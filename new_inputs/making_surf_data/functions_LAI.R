# three functions needed in make_lai.R

# 1. winter/austral crop LAI calculation, use it when d_sgs > d_egs
# follow Simpson_2012 algrithm (Sect. 5, Table 3 and S7.1)
lai_winter_crop = function(d_sgs = 213, d_egs = 123, LAI_min = 0, LAI_max = 3.5, l_s = 70, l_e = 20){
  days=1:365 # days in a year
  d_egs = d_egs + 365 # extend it to the next year
  LAI_long = rep(0, 365*2)
  LAI_tmp = rep(0, 365) # output LAI
  
  for(dn in 1:(2*length(days))){
    if(dn <= d_sgs | dn >= d_egs) LAI_long[dn] = LAI_min
    if(dn >= d_sgs & dn < d_sgs + l_s) LAI_long[dn] = LAI_min + (LAI_max - LAI_min)*(dn - d_sgs)/l_s
    if(dn >= d_sgs + l_s & dn <= d_egs - l_e) LAI_long[dn] = LAI_max
    if(dn > d_egs - l_e & dn < d_egs) LAI_long[dn] = LAI_min + (LAI_max - LAI_min)*(d_egs - dn)/l_e
  }
  
  for(i in 1:(2*length(days))){
    if(i <= length(days)) LAI_tmp[i] = LAI_long[i]
    else{
      if(LAI_tmp[i-365] == 0) LAI_tmp[i-365] = LAI_long[i]
    }
  }
  return(LAI_tmp) # return 365x1 LAI 
}
#a=lai_austral_crop()

# 2. boreal summer crop LAI calculation function, works when d_sgs < d_egs
# follow Simpson_2012 algrithm (Sect. 5, Table 3 and S7.1)
lai_summer_crop = function(d_sgs = 123, d_egs = 213, LAI_min = 0, LAI_max = 3.5, l_s = 70, l_e = 20){
  days=1:365
  LAI_tmp = rep(0, 365)
  
  for(dn in 1:length(days)){
    if(dn <= d_sgs | dn >= d_egs) LAI_tmp[dn] = LAI_min
    if(dn >= d_sgs & dn < d_sgs + l_s) LAI_tmp[dn] = LAI_min + (LAI_max - LAI_min)*(dn - d_sgs)/l_s
    if(dn >= d_sgs + l_s & dn <= d_egs - l_e) LAI_tmp[dn] = LAI_max
    if(dn > d_egs - l_e & dn < d_egs) LAI_tmp[dn] = LAI_min + (LAI_max - LAI_min)*(d_egs - dn)/l_e
  }
  return(LAI_tmp)
  
}
#b=lai_summer_crop()

# 3. take in daily LAI, output monthly mean
lai_monthly_mean = function(input = a){
  LAI_month = rep(0,12)
  
  days_in_months=c(31,28,31,30,31,30,31,31,30,31,30,31)
  end_of_each_month = rep(0,12)
  for(i in 1:12) end_of_each_month[i] = sum(days_in_months[1:i])
  end_of_each_month = c(1, end_of_each_month)
  
  for(j in 1:12){
    tmp1 = end_of_each_month[j]
    tmp2 = end_of_each_month[j+1]
    LAI_month[j] = mean(input[tmp1:tmp2])
  }
  return(LAI_month)
}

#c=lai_monthly_mean(input = a)
#plot(1:12,c, type = 'l')

