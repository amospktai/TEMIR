# This file contains useful functions for analyzing and visualizing geographical (spatial) data.

# Contact: Amos P. K. Tai (amostai@cuhk.edu.hk)

# Update history:
# Nov 2017 (Tai, Leung): Added the capacity to customize color schemes in plot.field() when using "custom.breaks", and added functionalities to plot legend/image only, specify the position of the legend, and new country features. Also created inv.terrain.colors() to display soil-vegetation related variables. Also modified plot.site() so that the map can be overlayed onto an existing one.

###############################################################################

# Approximate the distance between two locations from the latitudes and longitudes:
dist.latlon = function(lat1, lon1, lat2, lon2) {
	# This function calculates the distance (in km) between two locations on Earth, given their latitudes and longitudes.
	# lat/lon is in angle
	lat1 = lat1/180*pi
	lat2 = lat2/180*pi
	lon1 = lon1/180*pi
	lon2 = lon2/180*pi
	arc.angle = acos(sin(lat1)*sin(lat2) + cos(lat1)*cos(lat2)*cos(abs(lon2-lon1)))
	dist = arc.angle*6371.009
	# 6371.009 km is the mean Earth radius defined by the International Union of Geodesy and Geophysics. Alternatively, you may use 6371.22 km as in the Community Earth System Model.
	return(dist)
} 

###############################################################################

# Approximate the area of a spatial grid square from the latitudes and longitudes of the diagonal vertices:
area.latlon = function(lat1, lon1, lat2, lon2) {
	# This function calculates the area (in km^2) of a spatial grid square, given the latitudes and longitudes of the two diagonal vertices of the grid square.
	# lat/lon is in angle; lat: [-90:90]; lon:[-180:180].
	# lat1/lon1 and lat2/lon2 are thus the diagonal vertices of the square grid.
	lat1 = lat1/180*pi
	lat2 = lat2/180*pi
	lon1 = lon1/180*pi
	lon2 = lon2/180*pi
	A = abs(6371.009^2*(sin(lat2)-sin(lat1))*(lon2-lon1))
	# 6371.009 km is the mean Earth radius defined by the International Union of Geodesy and Geophysics. Alternatively, you may use 6371.22 km as in the Community Earth System Model.
	return(A)
}

###############################################################################

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

###############################################################################

# Interpolation method 1 - spatial averaging:
spavg = function(spdata, latlon, lat.frame, lon.frame, dlat=NULL, dlon=NULL) {
   
	# This function performs interpolation by spatial averaging from a vector of point data into a specified lat-lon frame.
	# "spdata" is a vector of spatial (point) data.
	# "latlon" has two columns: 1. latitude; 2. longitude; each pair of lat/lon corresponding to each datum in "spdata".
	# "lon.frame" has to be in increasing order and "lat.frame" can be in either order.
	# This function utilizes another function "find.lat.lon".
   
	lat.frame.orig = lat.frame
	if (lat.frame[length(lat.frame)] > lat.frame[1]) {
		# 'lat.frame' is in increasing order. Need to reverse it to proceed further.
		lat.frame = rev(lat.frame)
	}
	if (is.null(dlat)) dlat = abs(lat.frame[2]-lat.frame[1])
	if (is.null(dlon)) dlon = abs(lon.frame[2]-lon.frame[1])
	interpol = matrix(0, nrow=length(lon.frame), ncol=length(lat.frame))
	MAT = cbind(spdata, latlon)
	MAT = na.omit(MAT)
	spdata = MAT[,1]
	latlon = MAT[,2:3]
	if (nrow(MAT) == 0) interpol = NaN
	else if (nrow(MAT) == 1) {
		ind = find.lon.lat(latlon[2], latlon[1], lon.frame, lat.frame)
		if (length(ind) == 1) interpol = NaN
		else interpol[ind[1],ind[2]] = spdata
	}
	else {
		for (i in 1:length(lon.frame)) {
			for (j in 1:length(lat.frame)) {
				ind = which(latlon[,1] >= (lat.frame[j]-dlat/2) & latlon[,1] <= (lat.frame[j]+dlat/2) & latlon[,2] >= (lon.frame[i]-dlon/2) & latlon[,2] <= (lon.frame[i]+dlon/2))
				interpol[i,j] = mean(spdata[ind], na.rm=TRUE)
			}
		}
	}
	if (lat.frame.orig[length(lat.frame.orig)] > lat.frame.orig[1]) interpol = interpol[,length(lat.frame.orig):1]
	return(interpol)
}

###############################################################################

# Interpolation method 2 - inverse distance weighting:
invdist = function(spdata, latlon, lat.frame, lon.frame, n=2, dlim=500) {
   
	# This function performs interpolation by IDW from a vector of point data into a specified lat-lon frame.
	# "spdata" is a vector of spatial (point) data.
	# "latlon" has two columns: 1. latitude; 2. longitude; each pair of lat/lon corresponding to each datum in "spdata".
	# "lon.frame" has to be in increasing order and "lat.frame" can be in either order.
	# "n" is the power parameter.
	# "dlim" (in km) is the maximum search distance only within which the data are counted for interpolation.
	# This function utilizes another functions "dist.latlon" and "find.lon.lat".
   
	lat.frame.orig = lat.frame
	# If 'lat.frame' is in increasing order, we need to reverse it to proceed further:
	if (lat.frame[length(lat.frame)] > lat.frame[1]) lat.frame = rev(lat.frame) 
	interpol = matrix(0, nrow=length(lon.frame), ncol=length(lat.frame))
	MAT = cbind(spdata, latlon)
	MAT = na.omit(MAT)
	spdata = MAT[,1]
	latlon = MAT[,2:3]
	if (nrow(MAT) == 0) {
	   interpol = NaN
	} else if (nrow(MAT) == 1) {
		ind = find.lon.lat(latlon[2], latlon[1], lon.frame, lat.frame)
		if (length(ind) == 1) interpol = NaN else interpol[ind[1],ind[2]] = spdata
	} else {
		for (i in 1:length(lon.frame)) {
			for (j in 1:length(lat.frame)) {
				dist = dist.latlon(lat.frame[j], lon.frame[i], latlon[,1], latlon[,2])
				weight = 1/dist^n
				ind = which(dist <= dlim)
				interpol[i,j] = spdata[ind]%*%weight[ind]/sum(weight[ind], na.rm=TRUE)
			}
		}
	}
	if (lat.frame.orig[length(lat.frame.orig)] > lat.frame.orig[1]) interpol = interpol[,length(lat.frame.orig):1]
	return(interpol)
}

###############################################################################

# Standard deviation of sample data in a square grid:
spsd = function(spdata, latlon, lat.frame, lon.frame) {
   
   # This function finds the standard deviation of spatially averaged interpolated fields (using function "spavg") for each grid square.
   # "spdata" is a vector of spatial (point) data.
   # "latlon" has two columns: 1. latitude; 2. longitude; each pair of lat/lon corresponding to each datum in "spdata".
   # "lon.frame" has to be in increasing order and "lat.frame" has to be in decreasing order.
   # This function utilizes another function "find.lat.lon".
   
	dlon = abs(lon.frame[2]-lon.frame[1])
	dlat = abs(lat.frame[2]-lat.frame[1])
	stdev = matrix(0, nrow=length(lon.frame), ncol=length(lat.frame))
	MAT = cbind(spdata, latlon)
	MAT = na.omit(MAT)
	spdata = MAT[,1]
	latlon = MAT[,2:3]
	if (nrow(MAT) <= 2) stdev = NaN else {
		for (i in 1:length(lon.frame)) {
			for (j in 1:length(lat.frame)) {
				ind = which(latlon[,1] >= (lat.frame[j]-dlat/2) & latlon[,1] <= (lat.frame[j]+dlat/2) & latlon[,2] >= (lon.frame[i]-dlon/2) & latlon[,2] <= (lon.frame[i]+dlon/2))
				if (length(ind) < 2) stdev[i,j] = NaN else stdev[i,j] = sd(spdata[ind], na.rm=TRUE)
			}
		}
	}
	return(stdev)
}

###############################################################################

# Calculating average quantity of neighboring grid squares:
neighbor.avg <- function(spdata) {
   # This function calculates, for each grid square in an interpolated field, the mean value of a variable in the neighboring grid squares.
   # "spdata" is a matrix of interpolated spatial data.
	navg <- matrix(0, nrow=nrow(spdata), ncol=ncol(spdata))
	for (i in 1:nrow(navg)) {
		for (j in 1:ncol(navg)) {
			if (i == 1 & j == 1)
				ndata <- spdata[i:(i+1),j:(j+1)]
			else if (i == nrow(navg) & j == ncol(navg))
				ndata <- spdata[(i-1):i,(j-1):j]
			else if (i == 1 & j == ncol(navg))
				ndata <- spdata[i:(i+1),(j-1):j]
			else if (i == nrow(navg) & j == 1)
				ndata <- spdata[(i-1):i,j:(j+1)]
			else if (i == 1)
				ndata <- spdata[i:(i+1),(j-1):(j+1)]
			else if (j == 1)
				ndata <- spdata[(i-1):(i+1),j:(j+1)]
			else if (i == nrow(navg))
				ndata <- spdata[(i-1):i,(j-1):(j+1)]
			else if (j == ncol(navg))
				ndata <- spdata[(i-1):(i+1),(j-1):j]
			else
				ndata <- spdata[(i-1):(i+1),(j-1):(j+1)]
			ndata <- na.omit(as.vector(ndata))
			navg[i,j] <- (mean(ndata)*length(ndata)-spdata[i,j])/(length(ndata)-1)
			}
		}
	return(navg)
	}
	
###############################################################################

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

# "sp.regrid" and "sp.regrid.3D" used to be separate functions. Now they have been merged into one "sp.regrid".
# Also, these functions used to be called "sp.dissolve" and "sp.dissolve.3D" because the primary purpose was to dissolve higher resolution into lower resolution.
# For backward compatibility:
sp.regrid.3D = sp.regrid
sp.dissolve = sp.regrid
sp.dissolve.3D = sp.regrid.3D


###############################################################################

# Flip longitude from 0-360 to -180-180:

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

###############################################################################

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

###############################################################################

# Moran's I spatial autocorrelation for gridded data:
moran.I = function(spdata) {
	
	# This function calculates the Moran's I spatial autocorrelation index for a matrix of gridded spatial data with an adjacency matrix (w_ij = 1 if i and j are adjacent; i.e., if a grid cell is one of the eight neiboring cells).
	
	i.ind = matrix(rep(1:nrow(spdata), times=ncol(spdata)), nrow=nrow(spdata), byrow=FALSE)
	j.ind = matrix(rep(1:ncol(spdata), times=nrow(spdata)), nrow=nrow(spdata), byrow=TRUE)
	i.ind.vec = as.vector(i.ind)
	j.ind.vec = as.vector(j.ind)
	spdata.vec = as.vector(spdata)
   MAT = na.omit(cbind(i.ind.vec, j.ind.vec, spdata.vec))
   i.ind.vec = MAT[,1]
   j.ind.vec = MAT[,2]
   spdata.vec = MAT[,3]
	weight.mat = matrix(0, nrow=length(spdata.vec), ncol=length(spdata.vec))
	for (I in 1:length(spdata.vec)) {
		for (J in 1:length(spdata.vec)) {
			if (abs(diff(i.ind.vec[I], i.ind.vec[J])) <= 1 & abs(diff(j.ind.vec[I], j.ind.vec[J])) <= 1) weight.mat[I,J] = 1
		}
	}
	diag(weight.mat) = 0
   sum.cov = 0
   y = spdata.vec
	for (I in 1:length(y)) {
	   for (J in 1:length(y)) {
         sum.cov = sum.cov + weight.mat[I,J]*(y[I] - mean(y))*(y[J] - mean(y))
	   }
	}
	I.value = length(y)/sum((y - mean(y))^2)*sum.cov/sum(weight.mat)
	return(I.value)
}

###############################################################################

# Calculating the corresponding threshold R-value for a given p-value:
find.Rlim = function(pval, n) {
	# "pval" is the desired p-value; "n" is the sample size.
	t0 = qt(p=(1-pval/2), df=(n-1))	# The corresponding t-statistics
	R = sqrt(t0^2/(n-2+t0^2))			# The corresponding R-value (correlation coefficient)
	return(R)
}

###############################################################################

# Finding correlation coefficients between two variables for multiple locations:
find.cor = function(X, Y, plim=0.05, normalized=FALSE) {
	
	# This function finds the correlation coefficients that are statistically significant (defined by a desired p-value "plim") for multiple locations along the third dimension of 3-dimensional arrays (usually time).
	# X and Y are arrays of data with the first two dimensions defining the spatial frame (e.g. dim1 = longitude; dim2 = latitude; dim3 = time), and the correlation coefficients are found for the vectors of data contained in the third dimension.
	# If "normalized" is TRUE, the inputs X and Y will be normalized to their standard deviations.
	# This function utilizes another function "find.Rlim" to define the threshold R-value corresponding to the desired p-value.
	
	cor.XY = matrix(0, nrow=nrow(X), ncol=ncol(Y))
	for (i in 1:nrow(cor.XY)) {
		for (j in 1:ncol(cor.XY)) {
			if (length(unique(X[i,j,])) <= 2)
				cor.XY[i,j] = NaN
			else if (length(unique(Y[i,j,])) <= 2)
				cor.XY[i,j] = NaN
			else {
				MAT = cbind(X[i,j,], Y[i,j,])
				MAT = na.omit(MAT)
				if (normalized) {
					MAT[,1] = (MAT[,1] - mean(MAT[,1]))/sd(MAT[,1])
					MAT[,2] = (MAT[,2] - mean(MAT[,2]))/sd(MAT[,2])
				}
				Rlim = find.Rlim(plim, nrow(MAT))
				cor.MAT = cor(MAT[,1], MAT[,2])
				if (abs(cor.MAT) < Rlim) cor.XY[i,j] = NaN
				else cor.XY[i,j] = cor.MAT
			}
		}
	}
	return(cor.XY)
}

###############################################################################

# Quiver (vector arrow) plot:

quiver = function(u, v, xaxis, yaxis, xlab="x", ylab="y", scale=1, length=0.1, add=FALSE) {
   
   # This function creates quiver (vector arrow) plots for 2-dimensional data, "u" being the x-component and "v" being the y-component.
   # Source code from "http://fawn.unibw-hamburg.de/cgi-bin/Rwiki.pl?QuiverPlot".
   # "u" and "v" are matrices of 2D data, and "xaxis" and "yaxis" correspond to the x- and y-axis for the data.
   # "quiver" utilizes another function "par.uin" to determine scale of inches/userunits in x and y.
   
	# First stab at MATLAB's quiver in R.
	# From "http://tolstoy.newcastle.edu.au/R/help/01c/2711.html".
	# Robin Hankin Tue 20 Nov 2001 - 13:10:28 EST.
	xpos = matrix(rep(xaxis, times=length(yaxis)), ncol=length(yaxis))
	ypos = matrix(rep(yaxis, times=length(xaxis)), nrow=length(xaxis), byrow=TRUE)
   speed = sqrt(u*u+v*v) 
   maxspeed = max(speed) 
   u = u*scale/maxspeed 
   v = v*scale/maxspeed 
   matplot(xpos, ypos, type="p", cex=0, xlab=xlab, ylab=ylab, add=add) 
   arrows(xpos, ypos, xpos+u, ypos+v, length=length*min(par.uin()))
   
}

par.uin = function() {
  	# Determine scale of inches/userunits in x and y.
  	# From "http://tolstoy.newcastle.edu.au/R/help/01c/2714.html".
  	# Brian Ripley Tue 20 Nov 2001 - 20:13:52 EST.
   u = par("usr")
   p = par("pin")
   c(p[1]/(u[2] - u[1]), p[2]/(u[4] - u[3]))
}

###############################################################################

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

###############################################################################

# Inverse terrain color scale:
inv.terrain.colors = function(n, alpha=1) {
   # This function creates a vector of colors from light gray red to yellow to deep green, thus useful for displaying soil-vegetation related variables when we do a field plot.
   # "n" is the number of color levels. "n" has to be an even number.
   if ((n <- as.integer(n[1L])) > 0) {
      k <- n%/%2
      h <- c(0/12, 2/12, 4/12)
      s <- c(0.1, 1, 1)
      v <- c(0.85, 0.95, 0.25)
      c(hsv(h = seq.int(h[1L], h[2L], length.out = k), s = seq.int(s[1L], s[2L], length.out = k), v = seq.int(v[1L], v[2L], length.out = k), alpha = alpha), hsv(h = seq.int(h[2L], h[3L], length.out = n - k + 1)[-1L], s = seq.int(s[2L], s[3L], length.out = n - k + 1)[-1L], v = seq.int(v[2L], v[3L], length.out = n - k + 1)[-1L], alpha = alpha))
   }
   else character()
}

###############################################################################

# Quick map plot:
# This function creates 1 to 9 map plots for spatial data using function "image.plot" from package "fields".
# Packages "fields" and "maps" must be pre-loaded.
# "spdata" is a matrix or an array of spatial data. The third dimension of "spdata" defines the number of plots.
# "type" is a vector of intended types of presentation, as explained below.
# Five types of presentation are supported:
# 1. "sign": variable that has both positive and negative values, good for comparing signs of correlations/effects/deviations.
# 2. "frac": variable that is a fraction or percentage, good for comparing proportions.
# 3. "abs": variable that has absolute values only, good for comparing magtitudes.
# 4. "def": user-defined scale and z-limits. (If chosen, must define the same same scale/limits for all plots.) 
# 5. Default: no specification for display.
# If you want all plots to have the same type, simply enter a string scalar for "type".
# When "same" is set TRUE, all plots will have the same scale ("zlim").

quick.image = function(spdata, lon.map, lat.map, type=NULL, same=FALSE, zlim=NULL, col=tim.colors(32), nlevel=32, title='', width=6.5, height=4, mai=c(0.5, 0.4, 0.2, 0.2), mgp=c(1.4, 0.5, 0), tcl=-0.3, ps=12, legend.mar=3, legend.width=1.2) {
	
	if (length(dim(spdata)) == 2) {
		num.plot = 1
		mfrow = c(1,1)
		quartz(title=title, width=width, height=height)
		par(mfrow=mfrow, mai=mai, mgp=mgp, tcl=tcl, ps=ps)
		if (is.null(type)) {
			zlim = c(min(na.omit(as.vector(spdata))), max(na.omit(as.vector(spdata))))
			}
		else if (type == "sign") {
			zlim = c(-max(abs(na.omit(as.vector(spdata)))), max(abs(na.omit(as.vector(spdata)))))
			col = rwb.colors(nlevel)
			}
		else if (type == "frac") {
			zlim = c(0,1)
			}
		else if (type == "abs") {
			zlim = c(0, max(na.omit(as.vector(spdata))))
			}
		else if (type == "def") {
			zlim = zlim
			}
		else {
			zlim = c(min(na.omit(as.vector(spdata))), max(na.omit(as.vector(spdata))))
			}
		image.plot(lon.map, lat.map, spdata, zlim=zlim, xlab="", ylab="", axis.args=list(mgp=c(0,0.3,0), tcl=-0.1), legend.mar=legend.mar, legend.width=legend.width, col=col, nlevel=nlevel)
		map("world", add=T)
		}
		
	else {
		num.plot = dim(spdata)[3]
		if (num.plot == 2)
			mfrow = c(2,1)
		else if (num.plot == 3)
			mfrow = c(3,1)
		else if (num.plot == 4)
			mfrow = c(2,2)
		else if (num.plot == 5)
			mfrow = c(3,2)
		else if (num.plot == 6)
			mfrow = c(3,2)
		else if (num.plot == 7)
			mfrow = c(4,2)
		else if (num.plot == 8)
			mfrow = c(4,2)
		else	# i.e. num.plot == 9
			mfrow = c(3,3)
		
		quartz(title=title, width=width, height=height)
		par(mfrow=mfrow, mai=mai, mgp=mgp, tcl=tcl, ps=ps)

		if (num.plot > 1 & length(type) == 1) # i.e. you want all the plots to have the same type.
			type = rep(type, times=num.plot)
		else { }
	
		for (k in 1:num.plot) {
			if (is.null(type[k])) {
				if (same == TRUE)
					zlim = c(min(na.omit(as.vector(spdata))), max(na.omit(as.vector(spdata))))
				else
					zlim = c(min(na.omit(as.vector(spdata[,,k]))), max(na.omit(as.vector(spdata[,,k]))))
				}
			else if (type[k] == "sign") {
				col = rwb.colors(nlevel)
				if (same == TRUE)
					zlim = c(-max(abs(na.omit(as.vector(spdata)))), max(abs(na.omit(as.vector(spdata)))))
				else
					zlim = c(-max(abs(na.omit(as.vector(spdata[,,k])))), max(abs(na.omit(as.vector(spdata[,,k])))))
				}
			else if (type[k] == "frac") {
				zlim = c(0,1)
				}
			else if (type[k] == "abs") {
				if (same == TRUE)
					zlim = c(0, max(na.omit(as.vector(spdata))))
				else
					zlim = c(0, max(na.omit(as.vector(spdata[,,k]))))
				}
			else if (type[k] == "def") {
				zlim = zlim
				}
			else {
				if (same == TRUE)
					zlim = c(min(na.omit(as.vector(spdata))), max(na.omit(as.vector(spdata))))
				else
					zlim = c(min(na.omit(as.vector(spdata[,,k]))), max(na.omit(as.vector(spdata[,,k]))))
				}
			image.plot(lon.map, lat.map, spdata[,,k], zlim=zlim, xlab="", ylab="", axis.args=list(mgp=c(0,0.3,0), tcl=-0.1), legend.mar=legend.mar, legend.width=legend.width, col=col, nlevel=nlevel)
			map("world", add=T)
			}
		}
	}

###############################################################################

# Field map plot:

plot.field = function(spdata, lon.map, lat.map, type=NULL, same=FALSE, zlim=NULL, image.only=FALSE, legend.only=FALSE, horizontal=FALSE, col=NULL, nlevel=32, mai=c(0.2, 0.2, 0.1, 0.1), mgp=c(1.4, 0.5, 0), tcl=-0.3, ps=12, legend.mar=3, legend.width=1.2, xaxt="n", yaxt="n", Pacific.centric=FALSE, country=NULL, custom.breaks=NULL, map.xlim=NULL, map.ylim=NULL) {
  
  # This function plots spatial field data using function "image.plot" from package "fields".
  # Packages "fields" and "maps" must be pre-loaded. If Chinese/Japan provinces, Package "mapdata" is needed.
  # "spdata" is a matrix or an array of spatial data.
  # "type" is a vector of intended types of presentation, as explained below.
  # Five types of presentation are supported:
  # 1. "sign": variable that has both positive and negative values, good for comparing signs of correlations/effects/deviations.
  # 2. "frac": variable that is a fraction or percentage, good for comparing proportions.
  # 3. "abs": variable that has absolute values only, good for comparing magtitudes.
  # 4. "def": user-defined scale and z-limits. 
  # 5. Default: no specification for display.
   
  # You can specify the color scheme for the legend by specifying "col" (a function for color scheme, e.g., "tim.colors"); apart from inputting a function, you can also input a character string of hex color codes.
  # You can also manually customize the numeric break points for the legend by specifying "custom.breaks" (a vector). For those who want to use custom.breaks together with hex color codes, you need to use a function colorRampPalette() from a package called "RColorBrewer".
  
  # Last modified by Danny Leung (30 Nov 2017):
   
  # New application: 1. plot.field() can now plot maps without legend by setting argument "image.only"; 2. plot.field() can now plot just the legend by setting argument "legend.only"; 3. the legend can be either horizontal at the bottom or vertical on the right by setting argument "horizontal=TRUE/FALSE".
  # When the argument "image.only=TRUE", plot.field() plots map without legend. Further specification of arguments about legends will be ignored in the function. It will be ignored if "legend.only=TRUE".
  # When the argument "legend.only=TRUE", plot.field() plots the legend without the map. "Please remember to set "image.only=FALSE"."image.only" will be ignored.
  # Maps for some countries are now available. 'china', 'states' for the US, 'japan', 'france', and so on. Set "country" to customize your country map.
  
  if (lat.map[length(lat.map)] < lat.map[1]) {
    # 'lat.map' is in decreasing order. Need to reverse it to proceed further.
    lat.map = rev(lat.map)
    spdata = spdata[,length(lat.map):1]
  }	
  par(mai=mai, mgp=mgp, tcl=tcl, ps=ps)
  
  if (is.null(custom.breaks)) {
    
    if (is.null(col)) col = tim.colors(nlevel) else if (is.function(col)) col = col(nlevel) else col = col
    if (is.null(type)) {
      zlim = c(min(na.omit(as.vector(spdata))), max(na.omit(as.vector(spdata))))
    } else if (type == 'sign') {
      if (is.null(zlim)) zlim = c(-max(abs(na.omit(as.vector(spdata)))), max(abs(na.omit(as.vector(spdata))))) else zlim = zlim
      # if (is.null(col)) col = rwb.colors(nlevel)
      # If "sign" is used, the color scheme is always rwb.colors(nlevel).
      col = rwb.colors(nlevel)
    } else if (type == 'frac') {
      zlim = c(0,1)
    } else if (type == 'abs') {
      #zlim = c(0, 3)
      zlim = c(0, max(na.omit(as.vector(spdata))))
    } else if (type == 'def') {
      zlim = zlim
    } else {
      zlim = c(min(na.omit(as.vector(spdata))), max(na.omit(as.vector(spdata))))
    }
    spdata[which(spdata > zlim[2])] = zlim[2]
    spdata[which(spdata < zlim[1])] = zlim[1]
    if (legend.only) image.only = FALSE
    if (image.only) {
       image(lon.map, lat.map, spdata, zlim=zlim, xlab='', ylab='', col=col, xaxt=xaxt, yaxt=yaxt)
    } else {
       image.plot(lon.map, lat.map, spdata, zlim=zlim, xlab='', ylab='', axis.args=list(mgp=mgp, tcl=tcl), legend.only=legend.only, horizontal=horizontal, legend.mar=legend.mar, legend.width=legend.width, col=col, nlevel=nlevel, xaxt=xaxt, yaxt=yaxt)
    }
    
  } else {
    
    nbreaks = length(custom.breaks)
    nlevel = nbreaks - 1
    if (!is.null(type)) {
      if (type == 'sign') col = rwb.colors(nlevel) else {
        if (is.null(col)) col = tim.colors(nlevel) else if (is.function(col)) col = col(nlevel) else col = colorRampPalette(col)(nlevel)
      }
    } else {
      if (is.null(col)) col = tim.colors(nlevel) else if (is.function(col)) col = col(nlevel) else col = colorRampPalette(col)(nlevel)
    }
    zval.breaks = 0:nlevel
    zval.center = 0:(nlevel - 1) + 0.5
    spdata.new = spdata
    for (n in 1:nlevel) spdata.new[which(spdata >= custom.breaks[n] & spdata <= custom.breaks[n + 1])] = zval.center[n]
    spdata.new[which(spdata < custom.breaks[1])] = zval.center[1]
    spdata.new[which(spdata > tail(custom.breaks, 1))] = tail(zval.center, 1)
    spdata = spdata.new
    zlim = c(0, nlevel)
    if (legend.only) image.only = FALSE
    if (image.only) {
       image(lon.map, lat.map, spdata, zlim=zlim, xlab='', ylab='', col=col, xaxt=xaxt, yaxt=yaxt, breaks=zval.breaks, lab.breaks=custom.breaks)
    } else {
       image.plot(lon.map, lat.map, spdata, zlim=zlim, xlab='', ylab='', axis.args=list(mgp=mgp, tcl=tcl), legend.only=legend.only, horizontal=horizontal, legend.mar=legend.mar, legend.width=legend.width, col=col, nlevel=nlevel, xaxt=xaxt, yaxt=yaxt, breaks=zval.breaks, lab.breaks=custom.breaks)
    }
    
  }
  
  # Map:
  if (Pacific.centric) map('world2', add=TRUE, xlim=map.xlim, ylim=map.ylim) else map('world', add=TRUE, xlim=map.xlim, ylim=map.ylim)
  if (!is.null(country)) map(country, add=TRUE, xlim=map.xlim, ylim=map.ylim, lwd=1)
  
}

###############################################################################

# Site map plot:

plot.site = function(spdata, latlon, xlim=NULL, ylim=NULL, pch=19, type=NULL, zlim=NULL, col=NULL, nlevel=32, mai=c(0.5, 0.4, 0.2, 0.2), mgp=c(1.4, 0.5, 0), tcl=-0.3, ps=12, legend.mar=3, legend.width=1.2, xaxt="n", yaxt="n", Pacific.centric=FALSE, add=FALSE) {
   
   # This function plots a map where individual site measurements are displayed as dots with colors indicating their values.
   # Packages "fields" and "maps" must be pre-loaded.
   # If there are more than one values for a given site, the mean will be calculated automatically.
   # "spdata" is a vector of site measurements; "latlon" is a matrix with the 1st and 2nd columns being the corresponding latitude and longitude, respectively.
   # "xlim" and "ylim" are the limits for the plot area; "col" is the vector of colors used to define values of "spdata"; "zlim" is the range of "spdata" values used to break down "col".
   # "type" is a vector of intended types of presentation, as explained below.
   # Four types of presentation are supported:
   # 1. "frac": variable that is a fraction or percentage, good for comparing proportions.
   # 2. "abs": variable that has absolute values only, good for comparing magtitudes.
   # 3. "def": user-defined scale and z-limits.
   # 4. Default: no specification for display.
   
	MAT = na.omit(cbind(spdata, latlon))
	spdata = MAT[,1]; latlon = MAT[,2:3]
   if (is.null(xlim)) xlim = range(latlon[,2])
   if (is.null(ylim)) ylim = range(latlon[,1])
	latlon.uniq = unique(latlon)
	spdata.uniq = rep(0, times=nrow(latlon.uniq))
	col.uniq = rep(0, times=nrow(latlon.uniq))
	for (i in 1:length(spdata.uniq)) {
		spdata.uniq[i] = mean(spdata[which(latlon[,1] == latlon.uniq[i,1] & latlon[,2] == latlon.uniq[i,2])], na.rm=TRUE)
	}
	if (is.null(type)) {
      zlim = c(min(spdata.uniq, na.rm=TRUE), max(spdata.uniq, na.rm=TRUE))
   } else if (type == 'sign') {
      if (is.null(zlim)) zlim = c(-max(abs(spdata.uniq)), max(abs(spdata.uniq))) else zlim = zlim
      if (is.null(col)) col = rwb.colors(nlevel)
   } else if (type == "frac") {
      zlim = c(0, 1)
   } else if (type == "abs") {
      zlim = c(0, max(spdata.uniq, na.rm=TRUE))
   } else if (type == "def") {
      zlim = zlim
   } else {
      zlim = c(min(spdata.uniq, na.rm=TRUE), max(spdata.uniq, na.rm=TRUE))
   }
	if (is.null(col)) col = tim.colors(nlevel)
	zbreak = seq(zlim[1], zlim[2], length=(length(col) + 1))
	for (i in 1:length(col)) col.uniq[which(spdata.uniq >= zbreak[i] & spdata.uniq <= zbreak[i+1])] = col[i]
	if (!add) {
	   par(mai=mai, mgp=mgp, tcl=tcl, ps=ps)
	   image.plot(seq(xlim[1], xlim[2], length=10), seq(ylim[1], ylim[2], length=10), matrix(zlim[1]-1, nrow=10, ncol=10), zlim=zlim, xlab="", ylab="", axis.args=list(mgp=c(0,0.5,0), tcl=-0.2), legend.mar=legend.mar, legend.width=legend.width, col=col, nlevel=nlevel, xaxt=xaxt, yaxt=yaxt)
	}
	points(latlon.uniq[,2], latlon.uniq[,1], col=col.uniq, pch=pch)
	if (!add) {
	   if (Pacific.centric) map('world2', add=TRUE) else map('world', add=TRUE)
	}
}

###############################################################################

# Model-observation comparison for concentrations:

plot.mod.obs.US = function(mod.data, site.data, lon.mod, lat.mod, method="spavg", n=2, dlim=500, xlab='Observed', ylab='Simulated', mai=c(0.4, 0.4, 0.2, 0.2), mgp=c(1.5, 0.5, 0), tcl=-0.25, ps=14, cex=1.4, divide.EW=TRUE) {
	
	# This function plots model results vs. observations for the contiguous U.S.; the regression line and 1-to-1 line are also plotted.
	# The plot also gives the R^2 and mean bias (b).
	# "mod.data" is a matrix of gridded model results with 1st dimension = "lon.mod" and 2nd dimension = "lat.mod"; "lon.mod" and "lat.mod" are the longitudes and latitudes of model results, both in increasing order.
	# "site.data" is an array of observations by site. It has 3 columns: 1st column = observed values; 2nd column = latitude of site; 3rd column = longitude of site.
	# "n" and "dlim" are parameters for external function "invdist".
	# The rest of parameters are for plotting.
	# The package/library "lmodel2" has to be loaded.
	
	library(lmodel2)

	if (method == "spavg") {
		obs.data = spavg(site.data[,1], site.data[,2:3], rev(lat.mod), lon.mod)
		obs.data = obs.data[,length(lat.mod):1]
		x.data = as.vector(obs.data)
		y.data = as.vector(mod.data)
		x.data.E = as.vector(obs.data[which(lon.mod >= -97.5),])
		x.data.W = as.vector(obs.data[which(lon.mod < -97.5),])
		y.data.E = as.vector(mod.data[which(lon.mod >= -97.5),])
		y.data.W = as.vector(mod.data[which(lon.mod < -97.5),])
		}
		
	else if (method == "invdist") {
		obs.data = invdist(site.data[,1], site.data[,2:3], rev(lat.mod), lon.mod, n=n, dlim=dlim)
		obs.data = obs.data[,length(lat.mod):1]
		x.data = as.vector(obs.data)
		y.data = as.vector(mod.data)
		x.data.E = as.vector(obs.data[which(lon.mod >= -97.5),])
		x.data.W = as.vector(obs.data[which(lon.mod < -97.5),])
		y.data.E = as.vector(mod.data[which(lon.mod >= -97.5),])
		y.data.W = as.vector(mod.data[which(lon.mod < -97.5),])
		}
		
	else if (method == "match.site") {
		mod.match = site.data
		for (i in 1:nrow(site.data)) {
			lon.lat.mod.ind = find.lon.lat(site.data[i,3], site.data[i,2], lon.mod, rev(lat.mod))
			lon.lat.mod.ind[2] = length(lat.mod) - lon.lat.mod.ind[2] + 1
			mod.match[i,1] = mod.data[lon.lat.mod.ind[1],lon.lat.mod.ind[2]]
			}
		x.data = site.data[,1]
		y.data = mod.match[,1]
		x.data.E = site.data[which(site.data[,3] >= -97.5),1]
		x.data.W = site.data[which(site.data[,3] < -97.5),1]
		y.data.E = mod.match[which(mod.match[,3] >= -97.5),1]
		y.data.W = mod.match[which(mod.match[,3] < -97.5),1]
		}
	
	xylim = c(min(c(x.data, y.data), na.rm=TRUE), max(c(x.data, y.data), na.rm=TRUE))
	
	par(mai=mai, mgp=mgp, tcl=tcl, ps=ps)
	
	if (divide.EW) {
		
		plot(x=x.data.E, y=y.data.E, pch=21, col="blue", xlim=xylim, ylim=xylim, xlab=xlab, ylab=ylab)
		points(x=x.data.W, y=y.data.W, pch=25, col="red")
		reg.E = lmodel2(y.data.E ~ x.data.E)
		reg.W = lmodel2(y.data.W ~ x.data.W)
		r.sq.E = reg.E$rsquare
		r.sq.W = reg.W$rsquare
		a.E = reg.E$regression.results[2,2]
		a.W = reg.W$regression.results[2,2]
		b.E = reg.E$regression.results[2,3]
		b.W = reg.W$regression.results[2,3]
		abline(a=0, b=1, col="black")
		abline(a=a.E, b=b.E, col="blue", lty=2)
		abline(a=a.W, b=b.W, col="red", lty=4)
		#abline(h=0, lty="dashed", col="gray"); abline(v=0, lty="dashed", col="gray")
	
		text(x=(xylim[1] + 0.000*diff(xylim)), y=(xylim[1] + 1.000*diff(xylim)), labels=expression(R^2), cex=cex, col='blue', adj=c(0,1))
		text(x=(xylim[1] + 0.080*diff(xylim)), y=(xylim[1] + 0.970*diff(xylim)), labels=paste('  = ', as.character(round(r.sq.E, digits=2)), sep=''), cex=cex, col='blue', adj=c(0,1))
		text(x=(xylim[1] + 0.050*diff(xylim)), y=(xylim[1] + 0.870*diff(xylim)), labels=paste('b = ', as.character(round(b.E, digits=3)), sep=''), cex=cex, col='blue', adj=c(0,1))
	
		text(x=(xylim[1] + 0.560*diff(xylim)), y=(xylim[1] + 0.430*diff(xylim)), labels=expression(R^2), cex=cex, col='red', adj=c(0,1))
		text(x=(xylim[1] + 0.640*diff(xylim)), y=(xylim[1] + 0.400*diff(xylim)), labels=paste('  = ', as.character(round(r.sq.W, digits=2)), sep=''), cex=cex, col='red', adj=c(0,1))
		text(x=(xylim[1] + 0.610*diff(xylim)), y=(xylim[1] + 0.300*diff(xylim)), labels=paste('b = ', as.character(round(b.W, digits=3)), sep=''), cex=cex, col='red', adj=c(0,1))
	
		LIST = list(eastern.US=reg.E, western.US=reg.W)
		
	} else {
		
		plot(x=x.data, y=y.data, pch=21, col="blue", xlim=xylim, ylim=xylim, xlab=xlab, ylab=ylab)
		reg = lmodel2(y.data ~ x.data)
		r.sq = reg$rsquare
		a = reg$regression.results[2,2]
		b = reg$regression.results[2,3]
		NMB = sum((y.data - x.data), na.rm=TRUE)/sum(x.data, na.rm=TRUE)
		abline(a=0, b=1, col="black")
		abline(a=a, b=b, col="blue", lty=2)
		#abline(h=0, lty="dashed", col="gray"); abline(v=0, lty="dashed", col="gray")
	
		text(x=(xylim[1] + 0.000*diff(xylim)), y=(xylim[1] + 1.000*diff(xylim)), labels=expression(R^2), cex=cex, col='blue', adj=c(0,1))
		text(x=(xylim[1] + 0.080*diff(xylim)), y=(xylim[1] + 0.970*diff(xylim)), labels=paste('  = ', as.character(round(r.sq, digits=2)), sep=''), cex=cex, col='blue', adj=c(0,1))
		text(x=(xylim[1] + 0.050*diff(xylim)), y=(xylim[1] + 0.870*diff(xylim)), labels=paste('b = ', as.character(round(b, digits=3)), sep=''), cex=cex, col='blue', adj=c(0,1))
		text(x=(xylim[1] + 0.050*diff(xylim)), y=(xylim[1] + 0.770*diff(xylim)), labels=paste('NMB = ', as.character(round(NMB, digits=3)), sep=''), cex=cex, col='blue', adj=c(0,1))
	
		LIST = list(reg=reg, NMB=NMB)
		
	}
		
	return(LIST)
	
}

###############################################################################

# Plot vertical profile:

plot.vprof = function(prof, lat, lev, type=NULL, lev.type='pres', same=FALSE, zlim=NULL, xlab='Latitude', ylab=NULL, col=NULL, nlevel=32, mai=c(0.5, 0.5, 0.2, 0.2), mgp=c(1.4, 0.5, 0), tcl=-0.3, ps=12, legend.mar=3, legend.width=1.2) {
   
   # This function plots vertical profile using function "image.plot" from package "fields".
   # Packages "fields" and "maps" must be pre-loaded.
   # "prof" is a matrix or an array of vertical profile data.
   # "lev" can either be pressure levels "pres" or altitudes "alt", specified by value of "lev.type".
   # "type" is a vector of intended types of presentation, as explained below.
   # Five types of presentation are supported:
   # 1. "sign": variable that has both positive and negative values, good for comparing signs of correlations/effects/deviations.
   # 2. "frac": variable that is a fraction or percentage, good for comparing proportions.
   # 3. "abs": variable that has absolute values only, good for comparing magtitudes.
   # 4. "def": user-defined scale and z-limits. (If chosen, must define the same same scale/limits for all plots.) 
   # 5. Default: no specification for display.
   # If you want all plots to have the same type, simply enter a string scalar for "type".
   
	if (lat[length(lat)] < lat[1]) {
		# 'lat' is in decreasing order. Need to reverse it to proceed further.
		lat = rev(lat)
		prof = prof[length(lat):1,]
	}
	if (lev.type == 'pres') {
		if (lev[length(lev)] > lev[1]) {
			# 'lev' is pressure levels and is now in increasing order. Need to reverse it to proceed further.
			lev = rev(lev)
			prof = prof[,length(lev):1]
		}
	} else {
		if (lev[length(lev)] < lev[1]) {
			# 'lev' is altitudes and is now in decreasing order. Need to reverse it to proceed further.
			lev = rev(lev)
			prof = prof[,length(lev):1]
		}
	}
	par(mai=mai, mgp=mgp, tcl=tcl, ps=ps)
	if (is.null(type)) {
		zlim = c(min(na.omit(as.vector(prof))), max(na.omit(as.vector(prof))))
	} else if (type == "sign") {
		if (is.null(zlim)) zlim = c(-max(abs(na.omit(as.vector(prof)))), max(abs(na.omit(as.vector(prof))))) else zlim = zlim
		if (is.null(col)) col = rwb.colors(32)
	} else if (type == "frac") {
		zlim = c(0,1)
	} else if (type == "abs") {
		zlim = c(0, max(na.omit(as.vector(prof))))
	} else if (type == "def") {
		zlim = zlim
	} else {
		zlim = c(min(na.omit(as.vector(prof))), max(na.omit(as.vector(prof))))
	}
	if (is.null(col)) col = tim.colors(32)
	if (lev.type == 'pres') {
		if (is.null(ylab)) ylab = 'Pressure (hPa)'
		image.plot(lat, -lev, prof, zlim=zlim, xlab=xlab, ylab=ylab, axis.args=list(mgp=c(0,0.5,0), tcl=-0.2), legend.mar=legend.mar, legend.width=legend.width, col=col, nlevel=nlevel)
	} else {
		if (is.null(ylab)) ylab = 'Altitude (km)'
		image.plot(lat, lev, prof, zlim=zlim, xlab=xlab, ylab=ylab, axis.args=list(mgp=c(0,0.5,0), tcl=-0.2), legend.mar=legend.mar, legend.width=legend.width, col=col, nlevel=nlevel)
	}
}

###############################################################################

# Convert country-level data to longitude-latitude grid:

CountryData2LonLat = function(in.df, lon, lat, country.code='ISO3', dlon=NULL, dlat=NULL) {
   
   # This function takes in country-by-country data and maps them into a user-defined lon-lat grid for analysis and display.
   # "in.df" is the input data frame with two columns of data. The first column is the country codes, and the second column is the data.
   # "lon" and "lat" are the user-defined vectors of longitudes and latitudes.
   # The type of country codes is set by "country.code". Only two types are allowed: "ISO3" (default), which is the 3-letter ISO 3166 standard code; and "UN", which is the numerical UN country code. See http://www.unc.edu/~rowlett/units/codes/country.htm for a list of ISO3 and UN codes.
   # If country.code = "ISO3", the first column of "in.df" should be characters; if "UN", it should be numeric.
	
   # This function requires installation of two external packages:
   library(rworldmap); library(raster)
   
	if (is.null(dlon)) dlon = (tail(lon, 1) - lon[1])/(length(lon) - 1)
	if (is.null(dlat)) dlat = (tail(lat, 1) - lat[1])/(length(lat) - 1)
	lonlat.grid = raster(nrow=length(lat), ncol=length(lon))
	extent(lonlat.grid) = c((lon[1] - dlon/2), (tail(lon, 1) + dlon/2), (lat[1] - dlat/2), (tail(lat, 1) + dlat/2))
	
	if (country.code == 'UN') in.df[,1] = as.numeric(in.df[,1])
	if (country.code == 'ISO3') in.df[,1] = as.character(in.df[,1])
	colnames(in.df) = c('country', 'value')

	# Deal with Czechoslovakia:
	if (country.code == 'UN') {
		old.country = 200
		new.countries = c(203, 703)
	} else {
		old.country = 'CSK'
		new.countries = c('CZE', 'SVK')
	}
	ind.old = which(in.df$country == old.country)
	for (n in 1:length(new.countries)) {
		ind.new = which(in.df$country == new.countries[n])
		if (length(ind.old) > 0 & length(ind.new) == 0) {
			in.df = rbind(in.df, c(NaN, NaN))
			in.df[nrow(in.df),1] = new.countries[n]
			in.df[nrow(in.df),2] = in.df[ind.old,2]
		} else in.df = in.df	
	}
	if (length(ind.old) > 0) in.df = in.df[-ind.old,]
	
	# Deal with Pacific Islands Trust Territories:
	if (country.code == 'UN') {
		old.country = 582
		new.countries = c(584, 583, 580, 585)
	} else {
		old.country = 'PCI'
		new.countries = c('MHL', 'FSM', 'MNP', 'PLW')
	}
	ind.old = which(in.df$country == old.country)
	for (n in 1:length(new.countries)) {
		ind.new = which(in.df$country == new.countries[n])
		if (length(ind.old) > 0 & length(ind.new) == 0) {
			in.df = rbind(in.df, c(NaN, NaN))
			in.df[nrow(in.df),1] = new.countries[n]
			in.df[nrow(in.df),2] = in.df[ind.old,2]
		} else in.df = in.df	
	}
	if (length(ind.old) > 0) in.df = in.df[-ind.old,]
	
	# Deal with USSR:
	if (country.code == 'UN') {
		old.country = 810
		new.countries = c(51, 31, 233, 268, 398, 417, 428, 440, 498, 643, 762, 795, 860)
	} else {
		old.country = 'SUN'
		new.countries = c('ARM', 'AZE', 'EST', 'GEO', 'KAZ', 'KGZ', 'LVA', 'LTU', 'MDA', 'RUS', 'TJK', 'TKM', 'UZB')
	}
	ind.old = which(in.df$country == old.country)
	for (n in 1:length(new.countries)) {
		ind.new = which(in.df$country == new.countries[n])
		if (length(ind.old) > 0 & length(ind.new) == 0) {
			in.df = rbind(in.df, c(NaN, NaN))
			in.df[nrow(in.df),1] = new.countries[n]
			in.df[nrow(in.df),2] = in.df[ind.old,2]
		} else in.df = in.df	
	}
	if (length(ind.old) > 0) in.df = in.df[-ind.old,]
	
	# Deal with Yugoslavia:
	if (country.code == 'UN') {
		old.country = 890
		new.countries = c(191, 705, 807, 70, 499, 688)
	} else {
		old.country = 'YUG'
		new.countries = c('HRV', 'SVN', 'MKD', 'BIH', 'MNE', 'SRB')
	}
	ind.old = which(in.df$country == old.country)
	for (n in 1:length(new.countries)) {
		ind.new = which(in.df$country == new.countries[n])
		if (length(ind.old) > 0 & length(ind.new) == 0) {
			in.df = rbind(in.df, c(NaN, NaN))
			in.df[nrow(in.df),1] = new.countries[n]
			in.df[nrow(in.df),2] = in.df[ind.old,2]
		} else in.df = in.df	
	}
	if (length(ind.old) > 0) in.df = in.df[-ind.old,]
	
	# Deal with Serbia and Montenegro:
	if (country.code == 'UN') {
		old.country = 891
		new.countries = c(499, 688)
	} else {
		old.country = 'SCG'
		new.countries = c('MNE', 'SRB')
	}
	ind.old = which(in.df$country == old.country)
	for (n in 1:length(new.countries)) {
		ind.new = which(in.df$country == new.countries[n])
		if (length(ind.old) > 0 & length(ind.new) == 0) {
			in.df = rbind(in.df, c(NaN, NaN))
			in.df[nrow(in.df),1] = new.countries[n]
			in.df[nrow(in.df),2] = in.df[ind.old,2]
		} else in.df = in.df	
	}
	if (length(ind.old) > 0) in.df = in.df[-ind.old,]
	
	out.poly = joinCountryData2Map(in.df, country.code, 'country')
	for (n in 1:nrow(in.df)) {
		if (country.code == 'UN') {
			if (length(which(out.poly$ISO_N3 == in.df$country[n])) == 0) print(paste(as.character(in.df$country[n]), ' is not found!', sep=''))
		} else if (country.code == 'ISO3') {
			if (length(which(as.character(out.poly$ISO_A3) == in.df$country[n])) == 0) print(paste(in.df$country[n], ' is not found!', sep=''))
		} else stop('Error!')
	}
	out.raster = rasterize(out.poly, lonlat.grid, field=out.poly$value)
	out = matrix(out.raster@data@values, nrow=length(lon), ncol=length(lat))[,length(lat):1]
	return(out)
	
}

###############################################################################

########################################################################
### The functions in the section below is very model- and region-specific...
########################################################################

# Removing data in grid squares that do not cover the contiguous United States:
# This function replaces values in non-U.S. grid squares with "NaN". X is a spatial field (matrix or 3-D array with 3rd dimension being time/date) at 2.5x2.5 NCEP/NCAR lon/lat resolution, with longitude ranging in increasing order from -125 to -65 (25 intervals in total), and longitude ranging in decreasing order from 50 to 22.5 (12 intervals in total).

rm.nonUS.2.5x2.5 = function(X) {
	for (i in 1:25) {
		if (i == 1) j = c(1:12)
		if (i == 2) j = c(1,7:12)
		if (i == 3) j = c(1,8:12)
		if (i == 4) j = c(1,9:12)
		if (i == 5) j = c(1,9:12)
		if (i == 6) j = c(1,9:12)
		if (i == 7) j = c(1,9:12)
		if (i == 8) j = c(1,9:12)
		if (i == 9) j = c(1,10:12)
		if (i == 10) j = c(1,10:12)
		if (i == 11) j = c(1,11,12)
		if (i == 12) j = c(1,11,12)
		if (i == 13) j = c(1,10:12)
		if (i == 14) j = c(1,10:12)
		if (i == 15) j = c(1,10:12)
		if (i == 16) j = c(1,2,10:12)
		if (i == 17) j = c(1,2,10:12)
		if (i == 18) j = c(1:3,11,12)
		if (i == 19) j = c(1:3,9,11,12)
		if (i == 20) j = c(1:3,8:12)
		if (i == 21) j = c(1,2,6:12)
		if (i == 22) j = c(1,2,5:12)
		if (i == 23) j = c(1,2,5:12)
		if (i == 24) j = c(1,2,4:12)
		if (i == 25) j = c(1:12)
		if (length(dim(X)) == 2) X[i,j] = NaN
		if (length(dim(X)) == 3) X[i,j,] = NaN
	}
	return(X)
}

# Since many older analyses used "rm.nonUS" as "rm.nonUS.2.5x2.5", we equate the two functions here:
rm.nonUS = rm.nonUS.2.5x2.5

###############################################################################

# Removing data in grid squares that do not cover the contiguous United States:
# This function replaces values in non-U.S. grid squares with "NaN". X is a spatial field (matrix or 3-D array with 3rd dimension being time/date) at 2x2.5 GEOS-Chem lon/lat resolution, with longitude ranging in increasing order from -125 to -67.5 (24 intervals in total), and longitude ranging in increasing order from 24 to 50 (14 intervals in total).

rm.nonUS.2x2.5 = function(spdata) {
	for (i in 1:dim(spdata)[1]) {
		for (j in 1:dim(spdata)[2]) {
			if (length(dim(spdata)) == 2) {
				if (j == 1) spdata[i,j] = NaN
				else if (j == 2 & (i <= 11 | (i >= 13 & i <= 18) | i >= 20)) spdata[i,j] = NaN
				else if (j == 3 & (i <= 10 | (i >= 13 & i <= 17) | i >= 20)) spdata[i,j] = NaN
				else if (j == 4 & (i <= 8 | i >= 19)) spdata[i,j] = NaN
				else if (j == 5 & (i <= 4 | i >= 20)) spdata[i,j] = NaN
				else if (j == 6 & (i <= 2 | i >= 21)) spdata[i,j] = NaN
				else if (j == 7 & (i <= 2 | i >= 21)) spdata[i,j] = NaN
				else if (j == 8 & (i <= 1 | i >= 22)) spdata[i,j] = NaN
				else if (j == 9 & (i <= 1 | i >= 22)) spdata[i,j] = NaN
				else if (j == 10 & (i <= 1 | i >= 24)) spdata[i,j] = NaN
				else if (j == 11 & (i <= 1 | (i >= 19 & i <= 20) | (i >= 24))) spdata[i,j] = NaN
				else if (j == 12 & (i <= 1 | (i >= 18 & i <= 22))) spdata[i,j] = NaN
				else if (j == 13 & (i <= 1 | i >= 15)) spdata[i,j] = NaN
				else if (j == 14) spdata[i,j] = NaN
				else { }
			} else {
				if (j == 1) spdata[i,j,] = NaN
				else if (j == 2 & (i <= 11 | (i >= 13 & i <= 18) | i >= 20)) spdata[i,j,] = NaN
				else if (j == 3 & (i <= 10 | (i >= 13 & i <= 17) | i >= 20)) spdata[i,j,] = NaN
				else if (j == 4 & (i <= 8 | i >= 19)) spdata[i,j,] = NaN
				else if (j == 5 & (i <= 4 | i >= 20)) spdata[i,j,] = NaN
				else if (j == 6 & (i <= 2 | i >= 21)) spdata[i,j,] = NaN
				else if (j == 7 & (i <= 2 | i >= 21)) spdata[i,j,] = NaN
				else if (j == 8 & (i <= 1 | i >= 22)) spdata[i,j,] = NaN
				else if (j == 9 & (i <= 1 | i >= 22)) spdata[i,j,] = NaN
				else if (j == 10 & (i <= 1 | i >= 24)) spdata[i,j,] = NaN
				else if (j == 11 & (i <= 1 | (i >= 19 & i <= 20) | (i >= 24))) spdata[i,j,] = NaN
				else if (j == 12 & (i <= 1 | (i >= 18 & i <= 22))) spdata[i,j,] = NaN
				else if (j == 13 & (i <= 1 | i >= 15)) spdata[i,j,] = NaN
				else if (j == 14) spdata[i,j,] = NaN
				else { }
			}
		}
	}
	return(spdata)
}

###############################################################################

# Removing data in grid squares that do not cover the contiguous United States:
# This function replaces values in non-U.S. grid squares with "NaN". X is a spatial field (matrix or 3-D array with 3rd dimension being time/date) at 4x5 GISS lon/lat resolution, with longitude ranging in increasing order from -122.5 to -67.5 (12 intervals in total), and longitude ranging in increasing order from 26 to 50 (7 intervals in total).

rm.nonUS.GISS = function(spdata) {
	if (length(dim(spdata)) == 2) {
		spdata[1,1:2] = NaN
		spdata[2,1:2] = NaN
		spdata[3,1] = NaN
		spdata[4,1] = NaN
		spdata[5,1] = NaN
		spdata[7,1] = NaN
		spdata[8,c(1,7)] = NaN
		spdata[9,7] = NaN
		spdata[10,c(1,2,7)] = NaN
		spdata[11,c(1,2,3,7)] = NaN
		spdata[12,c(1,2,3,4,5,7)] = NaN
	} else {
		spdata[1,1:2,] = NaN
		spdata[2,1:2,] = NaN
		spdata[3,1,] = NaN
		spdata[4,1,] = NaN
		spdata[5,1,] = NaN
		spdata[7,1,] = NaN
		spdata[8,c(1,7),] = NaN
		spdata[9,7,] = NaN
		spdata[10,c(1,2,7),] = NaN
		spdata[11,c(1,2,3,7),] = NaN
		spdata[12,c(1,2,3,4,5,7),] = NaN
	}
	return(spdata)
}

###############################################################################

# Find regional subset for the contiguous United States:
# This function finds the [lon, lat] indices for the regional subset of spatial data over the contiguous U.S. The spatial data (matrix) has 2x2.5 GEOS-Chem lon/lat resolution, with longitude ranging in increasing order from -125 to -67.5 (24 intervals in total), and longitude ranging in increasing order from 24 to 50 (14 intervals in total).

find.region.2x2.5 = function(region.name) {

	pacificnw = rbind(c(2,10), c(2,11), c(2,12), c(2,13), c(3,10), c(3,11), c(3,12), c(3,13))
	pacificsw = rbind(c(2,8), c(2,9), c(3,6), c(3,7), c(3,8), c(3,9), c(4,6), c(4,7), c(5,5), c(5,6))
	innernw = rbind(c(4,10), c(4,11), c(4,12), c(4,13), c(5,10), c(5,11), c(5,12), c(5,13), c(6,10), c(6,11), c(6,12), c(6,13), c(7,10), c(7,11), c(7,12), c(7,13))
	innersw = rbind(c(4,8), c(4,9), c(5,7), c(5,8), c(5,9), c(6,5), c(6,6), c(6,7), c(6,8), c(6,9), c(7,5), c(7,6), c(7,7), c(7,8), c(7,9), c(8,5), c(8,6), c(8,7), c(8,8), c(8,9))
	greatplains = rbind(c(8,10), c(8,11), c(8,12), c(8,13), c(9,9), c(9,10), c(9,11), c(9,12), c(9,13), c(10,9), c(10,10), c(10,11), c(10,12), c(10,13), c(11,9), c(11,10), c(11,11), c(11,12), c(11,13), c(12,9), c(12,10), c(12,11), c(12,12), c(12,13))
	southcentral = rbind(c(9,4), c(9,5), c(9,6), c(9,7), c(9,8), c(10,4), c(10,5), c(10,6), c(10,7), c(10,8), c(11,3), c(11,4), c(11,5), c(11,6), c(11,7), c(11,8), c(12,2), c(12,3), c(12,4), c(12,5), c(12,6), c(12,7), c(12,8), c(13,4), c(13,5), c(13,6), c(13,7), c(13,8))
	midwest = rbind(c(13,9), c(13,10), c(13,11), c(13,12), c(13,13), c(14,9), c(14,10), c(14,11), c(14,12), c(14,13), c(15,9), c(15,10), c(15,11), c(15,12), c(16,9), c(16,10), c(16,11), c(16,12), c(17,9), c(17,10), c(17,11), c(17,12), c(18,9), c(18,10), c(18,11), c(19,9), c(19,10))
	southeast = rbind(c(14,4), c(14,5), c(14,6), c(14,7), c(14,8), c(15,4), c(15,5), c(15,6), c(15,7), c(15,8), c(16,4), c(16,5), c(16,6), c(16,7), c(16,8), c(17,4), c(17,5), c(17,6), c(17,7), c(17,8), c(18,3), c(18,4), c(18,5), c(18,6), c(18,7), c(18,8), c(19,2), c(19,3), c(19,5), c(19,6), c(19,7), c(19,8), c(20,6), c(20,7), c(20,8), c(21,8))
	northeast = rbind(c(20,9), c(20,10), c(21,9), c(21,10), c(21,11), c(22,10), c(22,11), c(23,10), c(23,11), c(23,12), c(24,12))
	
	if (region.name == "Northeast") region = northeast
	if (region.name == "Midwest") region = midwest
	if (region.name == "Southeast") region = southeast
	if (region.name == "South-central") region = southcentral
	if (region.name == "Great Plains") region = greatplains
	if (region.name == "Pacific NW") region = pacificnw
	if (region.name == "Pacific SW") region = pacificsw
	if (region.name == "Inner NW") region = innernw
	if (region.name == "Inner SW") region = innersw
	
	return(region)
	
	}

###############################################################################

# Find regional subset for the contiguous United States:

find.region.GISS = function(region.name) {

# This function finds the [lon, lat] indices for the regional subset of spatial data over the contiguous U.S. The spatial data (matrix) has 4x5 GISS lon/lat resolution, with longitude ranging in increasing order from -122.5 to -67.5 (13 intervals in total), and longitude ranging in increasing order from 26 to 50 (7 intervals in total).
	
	pacificnw = rbind(c(1, 5), c(1, 6), c(1, 7))
	pacificsw = rbind(c(1, 3), c(1, 4), c(2, 3), c(2, 4))
	innernw = rbind(c(2, 5), c(2, 6), c(2, 7), c(3, 5), c(3, 6), c(3, 7), c(4, 5), c(4, 6), c(4, 7))
	innersw = rbind(c(3, 2), c(3, 3), c(3, 4), c(4, 2), c(4, 3), c(4, 4))
	greatplains = rbind(c(5, 4), c(5, 5), c(5, 6), c(5, 7), c(6, 4), c(6, 5), c(6, 6), c(6, 7))
	southcentral = rbind(c(5, 2), c(5, 3), c(6, 1), c(6, 2), c(6, 3))
	midwest = rbind(c(7, 4), c(7, 5), c(7, 6), c(7, 7), c(8, 4), c(8, 5), c(8, 6), c(9, 4), c(9, 5), c(9, 6))
	southeast = rbind(c(7, 2), c(7, 3), c(8, 2), c(8, 3), c(9, 1), c(9, 2), c(9, 3), c(10, 3), c(10, 4))
	northeast = rbind(c(10, 5), c(10, 6), c(11, 4), c(11, 5), c(11, 6), c(12, 6))
		
	if (region.name == "Northeast") region = northeast
	if (region.name == "Midwest") region = midwest
	if (region.name == "Southeast") region = southeast
	if (region.name == "South-central") region = southcentral
	if (region.name == "Great Plains") region = greatplains
	if (region.name == "Pacific NW") region = pacificnw
	if (region.name == "Pacific SW") region = pacificsw
	if (region.name == "Inner NW") region = innernw
	if (region.name == "Inner SW") region = innersw
	
	return(region)
	
}

###############################################################################

# Find regional subset for the contiguous United States:

find.region.gcap = function(region.name) {

# This function finds the [lon, lat] indices for the regional subset of spatial data over the contiguous U.S. The spatial data (matrix) has 4x5 GCAP lon/lat resolution, with longitude ranging in increasing order from -125 to -65 (13 intervals in total), and longitude ranging in increasing order from 24 to 48 (7 intervals in total).
	
	pacificnw = rbind(c(1, 6), c(1, 7), c(2, 6), c(2, 7))
	pacificsw = rbind(c(1, 5), c(2, 4), c(2, 5), c(3, 3), c(3, 4))
	innernw = rbind(c(3, 6), c(3, 7), c(4, 6), c(4, 7))
	innersw = rbind(c(3, 5), c(4, 3), c(4, 4), c(4, 5))
	greatplains = rbind(c(5, 5), c(5, 6), c(5, 7), c(6, 5), c(6, 6), c(6, 7))
	southcentral = rbind(c(5, 3), c(5, 4), c(6, 2), c(6, 3), c(6, 4), c(7, 2), c(7, 3), c(7, 4))
	midwest = rbind(c(7, 5), c(7, 6), c(7, 7), c(8, 5), c(8, 6), c(8, 7), c(9, 5), c(9, 6), c(10, 5))
	southeast = rbind(c(8, 3), c(8, 4), c(9, 3), c(9, 4), c(10, 2), c(10, 3), c(10, 4), c(11, 4))
	northeast = rbind(c(11, 5), c(11, 6), c(12, 6))
		
	if (region.name == "Northeast") region = northeast
	if (region.name == "Midwest") region = midwest
	if (region.name == "Southeast") region = southeast
	if (region.name == "South-central") region = southcentral
	if (region.name == "Great Plains") region = greatplains
	if (region.name == "Pacific NW") region = pacificnw
	if (region.name == "Pacific SW") region = pacificsw
	if (region.name == "Inner NW") region = innernw
	if (region.name == "Inner SW") region = innersw
	
	return(region)
	
	}

###############################################################################

# Find regional subset for the contiguous United States:
# This function finds the [lon, lat] indices for the regional subset of spatial data over the contiguous U.S. The spatial data (matrix) has 2.5x2.5 NCEP/NCAR lon/lat resolution, with longitude ranging in increasing order from -125 to -65 (25 intervals in total), and longitude ranging in decreasing order from 50 to 22.5 (12 intervals in total).

find.region.2.5x2.5 = function(region.name) {
	
	pacificnw = rbind(c(2,2), c(2,3), c(2,4), c(3,2), c(3,3), c(3,4))
	pacificsw = rbind(c(2,5), c(2,6), c(3,5), c(3,6), c(3,7), c(4,7), c(4,8), c(5,7), c(5,8))
	innernw = rbind(c(4,2), c(4,3), c(4,4), c(5,2), c(5,3), c(5,4), c(6,2), c(6,3), c(6,4), c(7,2), c(7,3), c(7,4))
	innersw = rbind(c(4,5), c(4,6), c(5,5), c(5,6), c(6,5), c(6,6), c(6,7), c(6,8), c(7,5), c(7,6), c(7,7), c(7,8), c(8,5), c(8,6), c(8,7), c(8,8))
	greatplains = rbind(c(8,2), c(8,3), c(8,4), c(9,2), c(9,3), c(9,4), c(9,5), c(10,2), c(10,3), c(10,4), c(10,5), c(11,2), c(11,3), c(11,4), c(11,5), c(12,2), c(12,3), c(12,4), c(12,5))
	southcentral = rbind(c(9,6), c(9,7), c(9,8), c(9,9), c(10,6), c(10,7), c(10,8), c(10,9), c(11,6), c(11,7), c(11,8), c(11,9), c(11,10), c(12,6), c(12,7), c(12,8), c(12,9), c(12,10), c(13,6), c(13,7), c(13,8), c(13,9))
	midwest = rbind(c(13,2), c(13,3), c(13,4), c(13,5), c(14,2), c(14,3), c(14,4), c(14,5), c(15,2), c(15,3), c(15,4), c(15,5), c(16,3), c(16,4), c(16,5), c(17,3), c(17,4), c(17,5), c(18,4), c(18,5), c(19,4), c(19,5))
	southeast = rbind(c(14,6), c(14,7), c(14,8), c(14,9), c(15,6), c(15,7), c(15,8), c(15,9), c(16,6), c(16,7), c(16,8), c(16,9), c(17,6), c(17,7), c(17,8), c(17,9), c(18,6), c(18,7), c(18,8), c(18,9), c(18,10), c(19,6), c(19,7), c(19,8), c(19,10), c(20,6), c(20,7))
	northeast = rbind(c(20,4), c(20,5), c(21,3), c(21,4), c(21,5),  c(22,3), c(22,4), c(23,3), c(23,4), c(24,3))
	
	if (region.name == "Northeast") region = northeast
	if (region.name == "Midwest") region = midwest
	if (region.name == "Southeast") region = southeast
	if (region.name == "South-central") region = southcentral
	if (region.name == "Great Plains") region = greatplains
	if (region.name == "Pacific NW") region = pacificnw
	if (region.name == "Pacific SW") region = pacificsw
	if (region.name == "Inner NW") region = innernw
	if (region.name == "Inner SW") region = innersw
	
	return(region)
	
}


########################################################################
### The functions in the section above are very model- and region-specific...
########################################################################

###############################################################################

########################################################################
### The functions below are for quick handling of ncdf files
########################################################################

# Data extraction from .nc file:
load.nc = function(filename, varid) {
   # This function simply extracts an data array associated with "varid" from a given .nc file associated with "filename".
   # Both "varid" and "filename" are a string scalar. "varid" is the name given to the variable of interest stored in the .nc file associated with "filename".
   library(ncdf4)
   nc = nc_open(filename)
   data = ncvar_get(nc, varid=varid)
   nc_close(nc)
   return(data)
}

###############################################################################

# Quick ncdf plotting:
# These functions allow a quick look at ncdf files (.nc)
# Quick look at info:
get.info.nc = function(filename) {
   library(ncdf4)
	nc = nc_open(filename)
	print(nc)
	nc_close(nc)
}
# Quick plot (requires functions "plot.field" in "get_geo.R"), and also print out summary statistics:
plot.nc = function(filename, var.name, dim3=NULL, dim4=NULL, lon=NULL, lat=NULL, lon.name='lon', lat.name='lat', summary.only=FALSE, print.all=FALSE, mfrow=c(1,1), mai=c(0.5, 0.5, 0.2, 0.2), mgp=c(1.4, 0.4, 0), tcl=-0.25, ps=12, type=NULL, zlim=NULL) {
	# "file" = character string for the input ncdf file
	# "var.name" = character string for the ame of variable to plot
	# "dim3" = integer scalar or vector for the indices of the 3rd dimension, if any, of the vaiable to plot
	# "dim4" = integer scalar or vector for the indices of the 4th dimension, if any, of the vaiable to plot
	# "lon" = vector for the longitudes of the variable; if not provided, values will be taken directly from the nc file
	# "lat" = vector for the latitudes of the variable; if not provided, values will be taken directly from the nc file
	# "lon.name" = short name for longitudes given to the nc file
	# "lat.name" = short name for latitudes given to the nc file
	# "summary.only" = if TRUE, will only give summary info of the variable requested without making a plot
	# "print.all" = if TRUE, will print out the values for all entries in the variable requested
	# The rest of the parameters are for functions "plot.field" and "par".
   library(ncdf4); library(fields); library(maps)
	nc = nc_open(filename)
	if (is.null(lon)) lon = ncvar_get(nc, lon.name)
	if (is.null(lat)) lat = ncvar_get(nc, lat.name)
	if (is.null(dim3) & is.null(dim4)) vars = ncvar_get(nc, var.name) else if (!is.null(dim3) & is.null(dim4)) vars = ncvar_get(nc, var.name)[,,dim3] else vars = ncvar_get(nc, var.name)[,,dim3,dim4]
	print('Dim:', quote=FALSE); print(dim(vars))
	print('Quantile:', quote=FALSE); print(quantile(vars, c(0, 0.25, 0.5, 0.75, 1), na.rm=TRUE))
	print('Mean:', quote=FALSE); print(mean(vars, na.rm=TRUE))
	print('Std dev:', quote=FALSE); print(sd(vars, na.rm=TRUE))
	if (print.all) {
		print('Values:', quote=FALSE)
		print(vars)
	} else { }
	if (!summary.only) {
		par(mfrow=mfrow, mai=mai, mgp=mgp, tcl=tcl, ps=ps)
		if (length(dim(vars)) == 2) {
			plot.field(vars, lon, lat, type=type, zlim=zlim)
		} else if (length(dim(vars)) == 3) {
			for (i in 1:(dim(vars)[3])) plot.field(vars[,,i], lon, lat, type=type, zlim=zlim)
		} else {
			for (i in 1:(dim(vars)[3])) {
				for (j in 1:(dim(vars)[4])) {
					plot.field(vars[,,i,j], lon, lat, type=type, zlim=zlim)
				}
			}
		}
	} else { }
	nc_close(nc)
}
		
###############################################################################

