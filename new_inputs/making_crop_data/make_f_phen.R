# prepare f_phen input data
# algrithm follows Simpson et al. 2012, supplementary Table S17 and S16
# used parameters for IAM_CR, same for all crops

# read in plant and harvest date, as start and end of growing season
setwd('~/Dropbox/Research/DOSE/TEMIR_201905/TEMIR_inputs/crops/')
load('crop_calendar.RData')

days = 1:365

f_phen = array(NA, dim=c(dim(crop_planting_date),length(days)))

# same set of paramters as Temperate/Boreal crops
f_phen_a = 0.1
f_phen_b = 0.1
f_phen_c = 1
f_phen_d = 0.1
f_phen_e = 0
f_phen_f = 45

for(i in 1:dim(f_phen)[1])
  for(j in 1:dim(f_phen)[2])
    for(k in 18:dim(f_phen)[3]){
      
      if(is.na(crop_planting_date[i,j,k])) next
      if(is.na(crop_harvest_date[i,j,k])) next
      
      d_sgs = crop_planting_date[i,j,k]
      d_egs = crop_harvest_date[i,j,k]
      A_start = d_sgs
      A_end = d_egs
      
      if(d_sgs < d_egs){ # boreal/summer crop
      for(dn in days){
        if(dn <= d_sgs | dn > d_egs) f_phen[i,j,k,dn] = 0
        if(dn > d_sgs & dn <= A_start) f_phen[i,j,k,dn] = f_phen_a
        if(dn > A_start & dn <= A_start+f_phen_e) f_phen[i,j,k,dn] = f_phen_b + (f_phen_c - f_phen_b)*(dn - A_start)/f_phen_e
        if(dn > A_start+f_phen_e & dn <= A_end - f_phen_f) f_phen[i,j,k,dn] = f_phen_c
        if(dn > A_end - f_phen_f & dn <= A_end) f_phen[i,j,k,dn] = f_phen_b + (f_phen_c - f_phen_d)*(A_end - dn)/f_phen_f
        if(dn > A_end & dn <= d_egs) f_phen[i,j,k,dn] = f_phen_d
      }
      }
      if(d_sgs > d_egs){ # austral/winter crop, which grows over a calendar year
        d_egs = d_egs + 365
        A_end = A_end + 365
        f_phen_long = rep(0,365*2)
        for(dn in 1:(2*length(days))){
          if(dn <= d_sgs | dn > d_egs) f_phen_long[dn] = 0
          else if(dn > d_sgs & dn <= A_start) f_phen_long[dn] = f_phen_a
          else if(dn > A_start & dn <= A_start+f_phen_e) f_phen_long[dn] = f_phen_b + (f_phen_c - f_phen_b)*(dn - A_start)/f_phen_e
          else if(dn > A_start+f_phen_e & dn <= A_end - f_phen_f) f_phen_long[dn] = f_phen_c
          else if(dn > A_end - f_phen_f & dn <= A_end) f_phen_long[dn] = f_phen_b + (f_phen_c - f_phen_d)*(A_end - dn)/f_phen_f
          else if(dn > A_end & dn <= d_egs) f_phen_long[dn] = f_phen_d
        }
        f_phen_tmp = rep(0,365)
        for(iday in 1:(365*2)){
          if(iday <= 365) f_phen_tmp[iday] = f_phen_long[iday]
          else{
            if(f_phen_tmp[iday-365] == 0) f_phen_tmp[iday-365] = f_phen_long[iday]
          }
        }
        f_phen[i,j,k,] = f_phen_tmp
      }
    }
# follow processed LAI dataset, save daily f_phen, dim=c(144,91,25,365)
setwd('~/Dropbox/Research/DOSE/TEMIR_201905/TEMIR_inputs/crops/')
save(f_phen, file = "f_phen.RData")
