################################################################################
### Useful functions
################################################################################

# This module contains useful functions for processing, analyzing and visualizing multidimensional (esp. spatiotemporal) data.

################################################################################

# Approximate the distance between two locations from the latitudes and longitudes:
dist.latlon = function(lat1, lon1, lat2, lon2) {
   # This function calculates the distance (in km) between two locations on Earth, given their latitudes and longitudes.
   # lat/lon is in angle
   lat1 = lat1/180*pi
   lat2 = lat2/180*pi
   lon1 = lon1/180*pi
   lon2 = lon2/180*pi
   arc.angle = acos(sin(lat1)*sin(lat2) + cos(lat1)*cos(lat2)*cos(abs(lon2-lon1)))
   dist = arc.angle*6371.22
   # 6371.22 km is the mean Earth radius used in the Community Earth System Model.
   return(dist)
} 

################################################################################

# Approximate the area of a spatial grid square from the latitudes and longitudes of the diagonal vertices:
area.latlon = function(lat1, lon1, lat2, lon2) {
   # This function calculates the area (in km^2) of a spatial grid square, given the latitudes and longitudes of the two diagonal vertices of the grid square.
   # lat/lon is in angle; lat: [-90:90]; lon:[-180:180].
   # lat1/lon1 and lat2/lon2 are thus the diagonal vertices of the square grid.
   lat1 = lat1/180*pi
   lat2 = lat2/180*pi
   lon1 = lon1/180*pi
   lon2 = lon2/180*pi
   A = abs(6371.22^2*(sin(lat2)-sin(lat1))*(lon2-lon1))
   # 6371.22 km is the mean Earth radius used in the Community Earth System Model.
   return(A)
}

################################################################################

# Finding the indices for a particular location given its longitude and latitude:
find.lon.lat = function(lonspec, latspec, lon, lat) {
   # This function finds the indices for a particular location (given "lonspec" and "latspec") in the regularly spaced vectors of longitude ("lon") and latitude ("lat") frames.
   lat.orig = lat
   if (lat[length(lat)] > lat[1]) {
      # 'lat' is in increasing order. Need to reverse it to proceed further.
      lat = rev(lat)
   }
   dlon = abs(lon[2]-lon[1])
   dlat = abs(lat[2]-lat[1])
   if (lonspec < (min(lon)-dlon/2) | lonspec > (max(lon)+dlon/2) | latspec > (max(lat)+dlat/2) | latspec < (min(lat)-dlat/2)) {
      return(NaN)
   } else {
      if (lonspec < min(lon)) lon1 = 1 else lon1 = max(which(lon <= lonspec))
      if (lonspec > max(lon)) lon2 = length(lon) else lon2 = min(which(lon >= lonspec))
      if (latspec < min(lat)) lat1 = length(lat) else lat1 = min(which(lat <= latspec))
      if (latspec > max(lat)) lat2 = 1 else lat2 = max(which(lat >= latspec))
      if (abs(lon[lon1]-lonspec) > abs(lon[lon2]-lonspec)) lon.ind = lon2 else lon.ind = lon1
      if (abs(lat[lat1]-latspec) > abs(lat[lat2]-latspec)) lat.ind = lat2 else lat.ind = lat1
      if (lat.orig[length(lat.orig)] > lat.orig[1]) lat.ind = length(lat.orig) - lat.ind + 1
      return(c(lon.ind, lat.ind))
   }
}

################################################################################

# Regridding spatial field data into another resolution, e.g., dissolving spatial field data into lower resolution (larger pixel size):

sp.regrid = function(spdata, lon.in, lat.in, lon.out, lat.out, method="mean", full.lon=FALSE, dlon.in=NULL, dlat.in=NULL, dlon.out=NULL, dlat.out=NULL) {
   
   # This function dissolves spatial field data (matrix or 3D array) into lower resolution with larger grid squares; it can also be used to regrid spatial field data into any other resolution.
   # "spdata" is a matrix with rows corresponding to longitudes ("lon.in") and columns corresponding to latitudes ("lat.in"), or a 3D array with the first two dimensions between longitude and latitude.
   # "lon.out" and "lat.out" are the lon/lat frames for the output resolution.
   # "lon.in" and "lon.out" should be in increasing order; "lat.in" and "lat.out" can be in either order, but they should have the same order.
   # There are two supported methods: 1. default "mean", i.e. calculating the area-weighted average; 2. "mode", take the value occupying the largest area in each output grid. For 3D array, however, only "method='mean'" is supported.
   # This function utilizes another functions "area.latlon" and "find.lon.lat".
   
   ndim = length(dim(spdata))
   if (ndim < 2 | ndim > 3) stop('Dimensions of spdata are not correct! Only 2D and 3D arrays are supported.')
   lat.in.orig = lat.in
   lat.out.orig = lat.out
   
   if (lat.in[length(lat.in)] > lat.in[1] & lat.out[length(lat.out)] > lat.out[1]) {
      # 'lat.in' and 'lat.out' are in increasing order. Need to reverse them to proceed further.
      lat.in = rev(lat.in)
      if (ndim == 2) spdata = spdata[,length(lat.in):1]
      if (ndim == 3) spdata = spdata[,length(lat.in):1,]
      lat.out = rev(lat.out)
   } else {
      # Do nothing.
   }
   if ( (lat.in[length(lat.in)] > lat.in[1] & lat.out[length(lat.out)] < lat.out[1]) | (lat.in[length(lat.in)] < lat.in[1] & lat.out[length(lat.out)] > lat.out[1]) ) stop('lat.in and lat.out should have the same increasing/decreasing order.')
   
   if (is.null(dlon.in)) dlon.in = abs(lon.in[2]-lon.in[1])
   if (is.null(dlat.in)) dlat.in = abs(lat.in[2]-lat.in[1])
   if (is.null(dlon.out)) dlon.out = abs(lon.out[2]-lon.out[1])
   if (is.null(dlat.out)) dlat.out = abs(lat.out[2]-lat.out[1])
   
   if (full.lon) {
      # If the data has full longitude range (e.g. 0:360, -180:180), we add matching lon to head and tail of lon.in.
      library(abind)
      lon.in.new = c(lon.in[1]-dlon.in, lon.in, lon.in[length(lon.in)]+dlon.in)
      if (ndim == 2) spdata = abind(spdata[length(lon.in),], spdata, spdata[1,], along=1)
      if (ndim == 3) spdata = abind(spdata[length(lon.in),,], spdata, spdata[1,,], along=1)
      lon.in = lon.in.new
   }
   
   lon.in.lim = c((lon.in-dlon.in/2), (lon.in[length(lon.in)]+dlon.in/2))
   lat.in.lim = c((lat.in+dlat.in/2), (lat.in[length(lat.in)]-dlat.in/2))
   
   if (ndim == 2) {
      
      sp.out = matrix(0, nrow=length(lon.out), ncol=length(lat.out))
      for (i in 1:length(lon.out)) {
         for (j in 1:length(lat.out)) {
            lon.out1 = lon.out[i]-dlon.out/2
            lon.out2 = lon.out[i]+dlon.out/2
            lat.out1 = lat.out[j]+dlat.out/2
            lat.out2 = lat.out[j]-dlat.out/2
            lon.frame = c(lon.out1, lon.in.lim[which(lon.in.lim >= lon.out1 & lon.in.lim <= lon.out2)], lon.out2)
            lat.frame = c(lat.out1, lat.in.lim[which(lat.in.lim <= lat.out1 & lat.in.lim >= lat.out2)], lat.out2)
            area = matrix(0, nrow=(length(lat.frame)-1), ncol=(length(lon.frame)-1))
            for (x in 1:nrow(area)) {
               for (y in 1:ncol(area)) {
                  area[x,y] = area.latlon(lat.frame[x], lon.frame[y], lat.frame[x+1], lon.frame[y+1])
               }
            }
            data = matrix(0, nrow=(length(lat.frame)-1), ncol=(length(lon.frame)-1))
            for (x in 1:nrow(data)) {
               for (y in 1:ncol(data)) {
                  lat.centroid = (lat.frame[x]+lat.frame[x+1])/2
                  lon.centroid = (lon.frame[y]+lon.frame[y+1])/2
                  ind = find.lon.lat(lon.centroid, lat.centroid, lon.in, lat.in)
                  lat.ind = ind[2]
                  lon.ind = ind[1]
                  data[x,y] = spdata[lon.ind,lat.ind]
               }
            }
            MAT = cbind(as.vector(data), as.vector(area))
            MAT = na.omit(MAT)
            data = MAT[,1]
            area = MAT[,2]
            if (method == "mean") {
               sp.out[i,j] <- sum(data*area)/sum(area)
            } else {
               data.unique = unique(data)
               area.unique = rep(0, times=length(data.unique))
               for (n in 1:length(area.unique)) {
                  area.unique[n] = sum(area[which(data == data.unique[n])])
               }
               max.ind = which(area.unique == max(area.unique))
               if (length(max.ind) == 1) sp.out[i,j] = data.unique[max.ind] else sp.out[i,j] = data.unique[max.ind[1]]
            }
         }
      }
      if (lat.out.orig[length(lat.out.orig)] > lat.out.orig[1]) sp.out = sp.out[,length(lat.out.orig):1]
      
   } else {
      
      # "spdata" are 3D.
      if (method != 'mean') stop('For 3D array, only method=mean is supported.')
      
      sp.out = array(NaN, dim=c(length(lon.out), length(lat.out), dim(spdata)[3]))
      for (i in 1:length(lon.out)) {
         for (j in 1:length(lat.out)) {
            lon.out1 = lon.out[i]-dlon.out/2
            lon.out2 = lon.out[i]+dlon.out/2
            lat.out1 = lat.out[j]+dlat.out/2
            lat.out2 = lat.out[j]-dlat.out/2
            lon.frame = c(lon.out1, lon.in.lim[which(lon.in.lim >= lon.out1 & lon.in.lim <= lon.out2)], lon.out2)
            lat.frame = c(lat.out1, lat.in.lim[which(lat.in.lim <= lat.out1 & lat.in.lim >= lat.out2)], lat.out2)
            area = array(NaN, dim=c(length(lat.frame)-1, length(lon.frame)-1, dim(spdata)[3]))
            for (x in 1:dim(area)[1]) {
               for (y in 1:dim(area)[2]) {
                  for (t in 1:dim(area)[3]) {
                     area[x,y,t] = area.latlon(lat.frame[x], lon.frame[y], lat.frame[x+1], lon.frame[y+1])
                  }
               }
            }
            data = array(NaN, dim=c(length(lat.frame)-1, length(lon.frame)-1, dim(spdata)[3]))
            for (x in 1:dim(data)[1]) {
               for (y in 1:dim(data)[2]) {
                  lat.centroid = (lat.frame[x]+lat.frame[x+1])/2
                  lon.centroid = (lon.frame[y]+lon.frame[y+1])/2
                  ind = find.lon.lat(lon.centroid, lat.centroid, lon.in, lat.in)
                  lat.ind = ind[2]
                  lon.ind = ind[1]
                  data[x,y,] = spdata[lon.ind,lat.ind,]
               }
            }
            sp.out[i,j,] = apply(data*area, 3, sum, na.rm=TRUE)/apply(area*!is.na(data), 3, sum, na.rm=TRUE)
         }
      }
      if (lat.out.orig[length(lat.out.orig)] > lat.out.orig[1]) sp.out = sp.out[,length(lat.out.orig):1,]
      
   }
   
   return(sp.out)
   
}

###############################################################################

# Flip spatial data with longitude from 0-360 to -180-180:

flip.lon = function(spdata, lon) {
   # This function reorganizes a global spatial data array with longitude coordinates from 0 to 360 to the more conventional -180 to 180.
   # "spdata" can be of any dimensions until eighth, but the first dimension has to be longitude.
   # "lon" has to be in increasing order.
   # This function requires the package "abind" to be installed and loaded.
   ind.flip = which(lon >= 180)
   lon.out = c((lon[ind.flip] - 360), lon[1:(ind.flip[1]-1)])
   ndim = length(dim(spdata))
   if (ndim == 2) spdata.out = rbind(spdata[ind.flip,], spdata[1:(ind.flip[1]-1),])
   if (ndim == 3) spdata.out = abind(spdata[ind.flip,,], spdata[1:(ind.flip[1]-1),,], along=1)
   if (ndim == 4) spdata.out = abind(spdata[ind.flip,,,], spdata[1:(ind.flip[1]-1),,,], along=1)
   if (ndim == 5) spdata.out = abind(spdata[ind.flip,,,,], spdata[1:(ind.flip[1]-1),,,,], along=1)
   if (ndim == 6) spdata.out = abind(spdata[ind.flip,,,,,], spdata[1:(ind.flip[1]-1),,,,,], along=1)
   if (ndim == 7) spdata.out = abind(spdata[ind.flip,,,,,,], spdata[1:(ind.flip[1]-1),,,,,,], along=1)
   if (ndim == 8) spdata.out = abind(spdata[ind.flip,,,,,,,], spdata[1:(ind.flip[1]-1),,,,,,,], along=1)
   out = list(spdata=spdata.out, lon=lon.out)
   return(out)
}

################################################################################

# Turning zero or negative values into NaN, and vice versa:

# Turning zero of negative values into NaN:
zero.to.NaN = function(X) {
   X[which(X <= 0)] = NaN
   return(X)
}

# Turning NaN values into zero:
NaN.to.zero = function(X) {
   X[which(is.na(X))] = 0
   return(X)
}

################################################################################

# Functions dealing with dates:

# Convert between NCEP/NCAP (Julian) time format (hour since year 0000) and more common YYYYMMDD format:
# These two functions convert between the two formats, and uses a reference year 1948 with 1948/01/01 = 17067072.
# Please note that the date is strictly for 00:00-00:00 UTC. Make the necessary time zone shift for converting to and from local date.

from.yyyymmdd = function(date) {
   # "date" and "time" cannot be a vector or matrix, and cannot be in the year of 1948.
   year = round(date/10000)
   month = round((date-year*10000)/100)
   day = date-year*10000-month*100
   year.diff = year-1948
   year.ind = seq(1, year.diff)
   add.day = 0
   for (i in 1:length(year.ind)) {
      if (((i-1)/4) == round((i-1)/4)){
        add.day = add.day + 366
      }else add.day = add.day + 365
   }
   month.ind = seq(1,month)
   for (j in 1:length(month.ind)) {
      if (j == 1) add.day = add.day + 0
      else if (j == 2 | j == 4 | j == 6 | j == 8 | j == 9 | j == 11) add.day = add.day + 31
      else if (j == 5 | j == 7 | j == 10 | j == 12) add.day = add.day + 30
      else if (j == 3 & (year.diff/4) == round(year.diff/4)) add.day = add.day + 29
      else add.day = add.day + 28
   }
   add.day = add.day + day
   time = 17067072 + (add.day-1)*24
   return(time) 
}

to.yyyymmdd = function(time) {
   # "date" and "time" cannot be a vector or matrix, and cannot be in the year of 1948.
   add.day = floor((time-17067072)/24) + 1
   year4 = ceiling(add.day/(366+365*3)) - 1
   year = 1948
   if (year4 == 0)	{
      add.day = add.day
      year = year
   } else {
      add.day = add.day - year4*(366+365*3)
      year = year + year4*4
   }
   y = c(366, 365, 365, 365)
   y = cumsum(y)
   if (add.day <= y[1]) {
      add.day = add.day
      year = year
   }
   else if (add.day > y[1] & add.day <= y[2]) {
      add.day = add.day - y[1]
      year = year + 1
   }
   else if (add.day > y[2] & add.day <= y[3]) {
      add.day = add.day - y[2]
      year = year + 2
   }
   else {
      add.day = add.day - y[3]
      year = year + 3
   }
   if (((year-1948)/4) == round((year-1948)/4)) m = c(31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
   else m = c(31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
   m = cumsum(m)
   if (add.day <= m[1]) {
      add.day = add.day
      month = 1
   }
   else if (add.day > m[1] & add.day <= m[2]) {
      add.day = add.day - m[1]
      month = 2
   }
   else if (add.day > m[2] & add.day <= m[3]) {
      add.day = add.day - m[2]
      month = 3
   }
   else if (add.day > m[3] & add.day <= m[4]) {
      add.day = add.day - m[3]
      month = 4
   }
   else if (add.day > m[4] & add.day <= m[5]) {
      add.day = add.day - m[4]
      month = 5
   }
   else if (add.day > m[5] & add.day <= m[6]) {
      add.day = add.day - m[5]
      month = 6
   }
   else if (add.day > m[6] & add.day <= m[7]) {
      add.day = add.day - m[6]
      month = 7
   }
   else if (add.day > m[7] & add.day <= m[8]) {
      add.day = add.day - m[7]
      month = 8
   }
   else if (add.day > m[8] & add.day <= m[9]) {
      add.day = add.day - m[8]
      month = 9
   }
   else if (add.day > m[9] & add.day <= m[10]) {
      add.day = add.day - m[9]
      month = 10
   }
   else if (add.day > m[10] & add.day <= m[11]) {
      add.day = add.day - m[10]
      month = 11
   }
   else {
      add.day = add.day - m[11]
      month = 12
   }
   day = add.day
   date = year*10000 + month*100 + day
   return(date)
}


# Make a vector of dates:
make.date.vec = function(start.date, end.date, return="date.vec") {
   # This function creates a vector of dates specified by "start.date" and "end.date", both of which are numeric, e.g. 20040101 for Jan 1, 2004.
   # This function utilizes other functions "from.yyyymmdd" and "to.yyyymmdd".
   time.vec = seq(from.yyyymmdd(start.date), from.yyyymmdd(end.date), by=24)
   date.vec = time.vec
   for (t in 1:length(date.vec)) date.vec[t] = to.yyyymmdd(time.vec[t])
   if (return == "date.vec") return(date.vec)
   if (return == "time.vec") return(time.vec)
}


# Day of year from date and vice versa:
# This function finds the day of year given a standard date yyyymmdd, or finds the date in yyyymmdd from the day of year given a year.
# Input argument "yyyymmdd", "day.of.year" or "yyyy" can be a vector.

date.to.day = function(yyyymmdd, leap=FALSE) {
   if (leap) days.in.month = c(0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31) else days.in.month = c(0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
   day = yyyymmdd - signif(yyyymmdd, 6)
   month = (yyyymmdd - signif(yyyymmdd, 4) - day)/100
   day.of.year = day + cumsum(days.in.month)[month]
   return(day.of.year)
}

day.to.date = function(day.of.year, yyyy, leap=FALSE) {
   if (leap) days.in.month = c(0, 31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31) else days.in.month = c(0, 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
   month = NULL
   for (n in 1:length(day.of.year)) month = c(month, max(which(cumsum(days.in.month) < day.of.year[n])))
   day = day.of.year - cumsum(days.in.month)[month]
   yyyymmdd = yyyy*1e4 + month*1e2 + day
   return(yyyymmdd)
}


# Leap year:
# This function determine whether a given year is a leap year:
is.leap = function(yyyy) {
   if (yyyy/4 == round(yyyy/4)) leap = TRUE else leap = FALSE
   return(leap)
}


################################################################################

# Plotting maps:

# Red-white-blue color scale:
rwb.colors = function(n) {
   # This function creates a vector of colors from dark blue to white to dark red, thus useful for displaying and contrasting positive versus negative values of a variable when we do a field plot.
   # "n" is the number of color levels. "n" has to be an even number.
   if (n/2 != round(n/2)) stop('n has to be an even number!')
   red = rep(0, times=n)
   green = rep(0, times=n)
   blue = rep(0, times=n)
   red[1:(n/2)] = seq(0, 1, length=(n/2))
   red[(n/2 + 1):floor(n*3/4)] = 1
   red[(floor(n*3/4) + 1):n] = seq(1, 0.5, length=ceiling(n/4))
   green[1:(n/2)] = seq(0, 1, length=(n/2))
   green[(n/2 + 1):n] = seq(1, 0, length=(n/2))
   blue[1:ceiling(n/4)] = seq(0.5, 1, length=ceiling(n/4))
   blue[(ceiling(n/4) + 1):(n/2)] = 1
   blue[(n/2 + 1):n] = seq(1, 0, length=(n/2))
   rwb = rgb(red, green, blue)
   return(rwb)
}

# Field map plot:
plot.field = function(spdata, lon.map, lat.map, type=NULL, same=FALSE, zlim=NULL, col=NULL, nlevel=32, mai=c(0.2, 0.2, 0.1, 0.1), mgp=c(1.4, 0.5, 0), tcl=-0.3, ps=12, legend.mar=3, legend.width=1.2, xaxt="n", yaxt="n", map.region='world', Pacific.centric=FALSE, custom.breaks=NULL, map.xlim=NULL, map.ylim=NULL) {
   
   # This function plots spatial field data using function "image.plot" from package "fields".
   # Packages "fields" and "maps" must be pre-loaded.
   # "spdata" is a matrix or an array of spatial data.
   # "type" is a vector of intended types of presentation, as explained below.
   # Five types of presentation are supported:
   # 1. "sign": variable that has both positive and negative values, good for comparing signs of correlations/effects/deviations.
   # 2. "frac": variable that is a fraction or percentage, good for comparing proportions.
   # 3. "abs": variable that has absolute values only, good for comparing magtitudes.
   # 4. "def": user-defined scale and z-limits. (If chosen, must define the same same scale/limits for all plots.) 
   # 5. Default: no specification for display.
   # If you want all plots to have the same type, simply enter a string scalar for "type".
   
   if (lat.map[length(lat.map)] < lat.map[1]) {
      # 'lat.map' is in decreasing order. Need to reverse it to proceed further.
      lat.map = rev(lat.map)
      spdata = spdata[,length(lat.map):1]
   }	
   par(mai=mai, mgp=mgp, tcl=tcl, ps=ps)
   if (is.null(custom.breaks)) {
      if (is.null(type)) {
         zlim = c(min(na.omit(as.vector(spdata))), max(na.omit(as.vector(spdata))))
      } else if (type == 'sign') {
         if (is.null(zlim)) zlim = c(-max(abs(na.omit(as.vector(spdata)))), max(abs(na.omit(as.vector(spdata))))) else zlim = zlim
         if (is.null(col)) col = rwb.colors(nlevel)
      } else if (type == 'frac') {
         zlim = c(0,1)
      } else if (type == 'abs') {
         zlim = c(0, max(na.omit(as.vector(spdata))))
      } else if (type == 'def') {
         zlim = zlim
      } else {
         zlim = c(min(na.omit(as.vector(spdata))), max(na.omit(as.vector(spdata))))
      }
      if (is.null(col)) col = tim.colors(nlevel)
      spdata[which(spdata > zlim[2])] = zlim[2]
      spdata[which(spdata < zlim[1])] = zlim[1]
      image.plot(lon.map, lat.map, spdata, zlim=zlim, xlab='', ylab='', axis.args=list(mgp=c(0,0.5,0), tcl=-0.2), legend.mar=legend.mar, legend.width=legend.width, col=col, nlevel=nlevel, xaxt=xaxt, yaxt=yaxt)
   } else {
      nbreaks = length(custom.breaks)
      nlevel = nbreaks - 1
      if (!is.null(type)) {
         if (type == 'sign') col = rwb.colors(nlevel) else col = tim.colors(nlevel)
      } else col = tim.colors(nlevel)
      zval.breaks = 0:nlevel
      zval.center = 0:(nlevel - 1) + 0.5
      spdata.new = spdata
      for (n in 1:nlevel) spdata.new[which(spdata >= custom.breaks[n] & spdata <= custom.breaks[n + 1])] = zval.center[n]
      spdata.new[which(spdata < custom.breaks[1])] = zval.center[1]
      spdata.new[which(spdata > tail(custom.breaks, 1))] = tail(zval.center, 1)
      spdata = spdata.new
      zlim = c(0, nlevel)
      image.plot(lon.map, lat.map, spdata, zlim=zlim, xlab='', ylab='', axis.args=list(mgp=c(0,0.5,0), tcl=-0.2), legend.mar=legend.mar, legend.width=legend.width, col=col, nlevel=nlevel, xaxt=xaxt, yaxt=yaxt, breaks=zval.breaks, lab.breaks=custom.breaks)
   }
   if (Pacific.centric) map('world2', add=TRUE, xlim=map.xlim, ylim=map.ylim) else map(map.region, add=TRUE, xlim=map.xlim, ylim=map.ylim)
}

################################################################################

################################################################################
### Functions for reading meteorological fields (MERRA2 & GEOS-FP)
################################################################################

# Function to convert units of MERRA2 or GEOS-FP variables to units of TEMIR variables

f_met_unit_convert = function(met.name, TEMIR.var, met.var) {
   
   # 'met.name' is the name of the chosen meteorological fields for the current TEMIR simulation
   # 'TEMIR.var' is the name of the meteorological variable used in TEMIR simulation
   # 'met.var' is the name of meteorological variable found in meteorological ncdf file
   
   # Get meteorological field:
   met_grid = get(x = TEMIR.var)
   
   # Convert met unit to TEMIR variable unit if required:
   if (met.name == 'GEOSFP') {
      # hPa to Pa
      if (TEMIR.var == 'SLP') out_grid = data_grid * 100
   } else if (met.name == 'MERRA2') {
      # Do nothing.
   } else stop('Unidentified meteorlogical fields are chosen!')
   
   if (!exists('out_grid')) stop(paste0(met.name, ' ', met.var, ' has different units than the ', TEMIR.var, ' used in TEMIR but no conversion method is specified!'))
   
   return(out_grid)
   
}

################################################################################
### End of file
################################################################################
