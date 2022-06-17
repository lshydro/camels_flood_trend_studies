###########################################################-
# Objective: Prepare climate input data for all CAMELS 
# catchments for trend analysis
# Author: Lina Stein, University of Bristol
# note: 
###########################################################-

#packages
library(rgdal)
library(sf)
library(raster)
library(ncdf4)
library(data.table)


fpath = "C:/Users/ls16959/Data/01_Climate/MSWEP_V2.2_010deg_20180326/"
fname = list.files(fpath, full.names = T)[1]

ALL_shp = readOGR(dsn = "C:/Users/ls16959/Data/CAMELS", layer = "ALL_shp")
cat_name_shape = paste(ALL_shp@data[,"country"], ALL_shp@data[,"gauge_id"], sep = "_")
save(cat_name_shape, file = "C:/Users/ls16959/Data/CAMELS/cat_name_shape.Rdata")

ncin <- nc_open(fname)
missValue = ncatt_get(ncin,'precipitation',"_FillValue")$value
nc_close(ncin)

temp_brick = brick(fname, varname = 'precipitation')
#Na value for brick is automatically correctly identified

s_time = Sys.time()
catweights_MSWEP_camels = extract(temp_brick[[1]], ALL_shp,cellnumbers=TRUE, weights=TRUE, small = T)
e_time = Sys.time()
e_time-s_time
#9.7 min
save(catweights_MSWEP_camels, file = paste0(fpath, "/MSWEP2_2_P_catweights_camels.Rdata"))




catweights_nrow = unlist(lapply(catweights_MSWEP_camels, function(x){nrow(x)}))
catweights_MSWEP_df = do.call(rbind, catweights_MSWEP_camels)
tempweights = catweights_MSWEP_df[,3]
cat_id = rep(1:length(catweights_nrow), times = catweights_nrow)

MSWEP_flist = list.files(fpath, full.names = T, pattern = ".nc")
system.time({
  out_df_list = lapply(MSWEP_flist, function(fname){
    out = tryCatch(
      {
        temp_brick = brick(fname, varname = 'precipitation')
        select_cell_df = temp_brick[catweights_MSWEP_df[,1]] 
        weighted_cells = data.table(sweep(select_cell_df, MARGIN=1, tempweights, '*'))
        out_df = weighted_cells[, lapply(.SD, sum), by=cat_id]
        return(out_df)
      },
      error=function(cond) {
        message(cond)
        # Choose a return value in case of error
        return(NA)
      },
      warning=function(cond) {
        message(cond)
        # Choose a return value in case of warning
        return(NA)
      }
    )
    return(out)
  })
})
#nearly two hours
save(out_df_list, file = paste0(fpath, "/MSWEP2_2_camels_out_df_list.Rdata"))
test = lapply(out_df_list, length)
# no missing months

fpath = "C:/Users/ls16959/Data/01_Climate/MSWEP_V2.2_010deg_20180326/"
load(file = paste0(fpath, "/MSWEP2_2_camels_out_df_list.Rdata"))
MSWEPdf = do.call(rbind, lapply(out_df_list, function(df){t(as.data.frame(df)[,-1])}))
MSWEP_date = as.Date(sub("X","", rownames(MSWEPdf)), format = "%Y.%m.%d")


# ERA5 soil moisture ------------------------------------------------------

library(raster)
library(rgdal)
library(sf)
library(ncdf4)
library(data.table)

ALL_shp = readOGR(dsn = "C:/Users/ls16959/Data/CAMELS", layer = "ALL_shp")

fpath = "C:/Users/ls16959/Data/01_Climate/ERA5-land"
fpath_out = "C:/Users/ls16959/Data/01_Climate/ERA5-land/CAMELS"


flist = list.files(fpath, full.names = T, pattern = "ERA5-Land_SoilWaterLayer")
attr_nc_list = c("swvl1", "swvl2", "swvl3")


ncin <- nc_open(flist[1])
missValue = ncatt_get(ncin,attr_nc_list[1],"_FillValue")$value
timev <- ncvar_get(ncin,"time")
nc_close(ncin)

temp_brick = brick(flist[1], varname = attr_nc_list[1])

s_time = Sys.time()
e1 = extent(c(xmin=0,xmax=180,ymin=-90,ymax=90))
e2 = extent(c(xmin=180, xmax=360, ymin=-90, ymax=90))
d1 = crop(temp_brick[[1]],e1)
d2 = crop(temp_brick[[1]],e2)

dm = rotate(merge(d1,d2))
NAvalue(dm)<-missValue

catweights_ERA5_camels = extract(dm, ALL_shp,cellnumbers=TRUE, weights=TRUE, small = T)
e_time = Sys.time()
e_time-s_time
#54s
save(catweights_ERA5_camels, file = paste0(fpath_out, "/catweights_ERA5_camels.Rdata"))

load(file = paste0(fpath_out, "/catweights_ERA5_camels.Rdata"))

catweights_nrow = unlist(lapply(catweights_ERA5_camels, function(x){nrow(x)}))
catweights_ERA5_df = do.call(rbind, catweights_ERA5_camels)
tempweights = catweights_ERA5_df[,3]
cat_id = rep(1:length(catweights_nrow), times = catweights_nrow)



s_time = Sys.time()
#Need to run one after another
for(file_num in c(3)){
  temp_brick = brick(flist[file_num], varname = attr_nc_list[file_num])
  
  #raster longitude is from 0 to 360, need to shift values over 180 to other side, shift for whole brick at once
  e1 = extent(c(xmin=0,xmax=180,ymin=-90,ymax=90))
  e2 = extent(c(xmin=180, xmax=360, ymin=-90, ymax=90))
  d1 = crop(temp_brick,e1)
  d2 = crop(temp_brick,e2)
  
  dm = rotate(merge(d1,d2))
  NAvalue(dm)<-missValue
  
  save(dm, file = paste0(fpath_out, "/ERA5_swvl", file_num,"_camels_dm.Rdata"))
  
  select_cell_df = dm[catweights_ERA5_df[,1]] 
  weighted_cells = data.table(sweep(select_cell_df, MARGIN=1, tempweights, '*'))
  out_df = weighted_cells[, lapply(.SD, sum), by=cat_id]
  save(out_df, file = paste0(fpath_out, "/ERA5_swvl", file_num,"_camels_out_df_list.Rdata"))
}
e_time = Sys.time()
e_time-s_time
#2-3 hours to run all three


flist_out = list.files(fpath_out, full.names = T, pattern = "_camels_out")
load(flist_out[1])
lyr1df = t(as.data.frame(out_df))[-1,]
load(flist_out[2])
lyr2df = t(as.data.frame(out_df))[-1,]
load(flist_out[3])
lyr3df = t(as.data.frame(out_df))[-1,]
#Volume of water in soil layer 1 (0 - 7 cm) 
#Volume of water in soil layer 2 (7 - 28 cm) 
#Volume of water in soil layer 3 (28-100 cm)

mean_lyr_test=((lyr1df[,1]*7)+(lyr2df[,1]*21)+(lyr3df[,1]*72))/100

SM_cat_df = do.call(cbind, lapply(c(1:ncol(lyr1df)), function(cat){
  mean_lyr=((lyr1df[,cat]*7)+(lyr2df[,cat]*21)+(lyr3df[,cat]*72))/100
  return(mean_lyr)
}))
save(SM_cat_df, file = paste0(fpath_out, "/ERA5_SM_cat_df.Rdata"))
#filename was just ERA5_SM.Rdata (but saved file before adding AUS was always called ERA5_SM_cat_df.Rdata?)
ERA5_date = as.Date(as.POSIXct(timev*3600,origin='1900-01-01 00:00:00.0'))
save(ERA5_date, file = paste0(fpath_out, "/ERA5_date.Rdata"))


# Compare to AU SM data ---------------------------------------------------

fpath_AU = "C:/Users/ls16959/Data/streamflow_countries/Australia/Wasko_streamflow/data_baseflow/"
flist_TS_AU = list.files(fpath_AU, full.names = T)

gauge_no_temp = sub( "C:/Users/ls16959/Data/streamflow_countries/Australia/Wasko_streamflow/data_baseflow/", "", flist_TS_AU[[1]])
gauge_no = sub(".txt", "", gauge_no_temp)
#read in observed streamflow from CAMELS data
temp = read.csv(flist_TS_AU[[1]], sep="")
temp_date = as.Date(temp$Date, format = "%Y-%m-%d")

temp_SMdf = temp[,c("Date", "sm.avg")]
temp_SMdf[,"ym_date"] = substr(as.character(temp_SMdf$Date), 1, 7)

which(cat_name_shape == paste0("AU_", gauge_no))

tempERA_SM = SM_cat_df[,which(cat_name_shape == paste0("AU_", gauge_no))]
tempERA_SMdf = data.frame(ym_date = substr(ERA5_date, 1,7), tempERA_SM)

tempSM_both = merge(temp_SMdf, tempERA_SMdf, by = "ym_date")
plot(tempSM_both[,c(3,4)])
cor(tempSM_both$sm.avg, tempSM_both$tempERA_SM)

cor_AU_SM = lapply(flist_TS_AU, function(xfile){
  gauge_no_temp = sub( "C:/Users/ls16959/Data/streamflow_countries/Australia/Wasko_streamflow/data_baseflow/", "", xfile)
  gauge_no = sub(".txt", "", gauge_no_temp)
  #read in observed streamflow from CAMELS data
  temp = read.csv(xfile, sep="")
  temp_date = as.Date(temp$Date, format = "%Y-%m-%d")
  
  temp_SMdf = temp[,c("Date", "sm.avg")]
  temp_SMdf[,"ym_date"] = substr(as.character(temp_SMdf$Date), 1, 7)
  
  cat_ind = which(cat_name_shape == paste0("AU_", gauge_no))
  
  if(length(cat_ind)==0){
    out = NA
  }else{
    tempERA_SM = SM_cat_df[,which(cat_name_shape == paste0("AU_", gauge_no))]
    tempERA_SMdf = data.frame(ym_date = substr(ERA5_date, 1,7), tempERA_SM)
    
    tempSM_both = merge(temp_SMdf, tempERA_SMdf, by = "ym_date")
    #plot(tempSM_both[,c(3,4)])
    out = cor(tempSM_both$sm.avg, tempSM_both$tempERA_SM)
  }
  return(out)
})

plot(unlist(cor_AU_SM))




# Berkeley Earth Surface Temperature -------------------------------------------------

#packages
library(rgdal)
library(sf)
library(raster)
library(ncdf4)
library(data.table)

library(abind)
library(lubridate)


#use version 2.2?
fpath = "C:/Users/ls16959/Data/01_Climate/BerkleyEarth/original_T_avg_20200713"
fname = list.files(fpath, full.names = T, pattern = "nc")[1]
load(file = "C:/Users/ls16959/Data/CAMELS/cat_name_shape.Rdata")

ALL_shp = readOGR(dsn = "C:/Users/ls16959/Data/CAMELS", layer = "ALL_shp")

temp_brick_climatology = brick(fname, varname = 'climatology')

s_time = Sys.time()
tempr = temp_brick_climatology[[1]]
#climatology and anomaly are same raster, either can be used to calculate weights
crs(tempr) = crs(ALL_shp)
catweights_BEST_camels = extract(readAll(tempr), ALL_shp,cellnumbers=TRUE, weights=TRUE, small = T) #readAll reads entire object into memory, which it would otherwise not have done due to size (?)
e_time = Sys.time()
e_time-s_time
#2min
save(catweights_BEST_camels, file = paste0(fpath, "/BEST_catweights_camels.Rdata"))





catweights_nrow = unlist(lapply(catweights_BEST_camels, function(x){nrow(x)}))
catweights_BEST_df = do.call(rbind, catweights_BEST_camels)
tempweights = catweights_BEST_df[,3]
cat_id = rep(1:length(catweights_nrow), times = catweights_nrow)

BEST_flist = list.files(fpath, full.names = T, pattern = ".nc")
system.time({
  out_df_list = lapply(BEST_flist, function(fname){
    print(fname)
    
    ncin <- nc_open(fname)
    clim_array <- ncvar_get(ncin,"climatology")
    anom_array <- ncvar_get(ncin,"temperature")
    day_of_year <- ncvar_get(ncin,"day_of_year")
    year_val <- ncvar_get(ncin,"year")
    #month_val <- ncvar_get(ncin,"month")
    #day_val <- ncvar_get(ncin,"day")
    nc_close(ncin)
    
    #climatology is lacking leap day, while anomaly array has leap day
    #calculate interpolated leap day climatology
    leap_day_clim = (clim_array[,,59]+clim_array[,,60])/2
    leap_year_clim = abind(clim_array[,,1:59], leap_day_clim, clim_array[,,60:365])
    
    #Create climatology array the same size as anomaly array
    clim_ar_list = lapply(unique(year_val), function(year_v){
      if(leap_year(year_v)){
        out = leap_year_clim
      }else{
        out = clim_array
      }
      return(out)
    })
    clim_ar = do.call(abind, clim_ar_list)
    
    if(dim(anom_array)[3]<dim(clim_ar)[3]){
      clim_ar = clim_ar[,,c(1:dim(anom_array)[3])]
    }
    
    #calculate temperature array
    temperature_array =  clim_ar+anom_array
    #load anomaly as a brick for crs information
    temp_brick_temperature = brick(fname, varname = 'temperature')
    #fill with temperature values
    temp_brick_temperature_filled = setValues(temp_brick_temperature, temperature_array)
    #flip for correct alignment
    tempr_brick_temperature =  flip(temp_brick_temperature_filled, direction = "y")
    
    select_cell_df = tempr_brick_temperature[catweights_BEST_df[,1]] 
    weighted_cells = data.table(sweep(select_cell_df, MARGIN=1, tempweights, '*'))
    out_df = weighted_cells[, lapply(.SD, sum), by=cat_id]
    return(out_df)
  })
})
#2 min hours
save(out_df_list, file = paste0(fpath, "/BEST_camels_out_df_list.Rdata"))

BESTdf = do.call(rbind, lapply(out_df_list, function(df){t(as.data.frame(df)[,-1])}))
save(BESTdf, file = paste0(fpath, "/BEST_camels_df.Rdata"))



system.time({
  BEST_date_list = lapply(BEST_flist, function(fname){
    print(fname)
    
    ncin <- nc_open(fname)
    
    year_val <- ncvar_get(ncin,"year")
    month_val <- ncvar_get(ncin,"month")
    day_val <- ncvar_get(ncin,"day")
    nc_close(ncin)
    BEST_date = as.Date(paste0(year_val,"-", month_val, "-", day_val), format = "%Y-%m-%d")
    return(BEST_date)
  })
})
BEST_date = do.call(c, BEST_date_list)
save(BEST_date, file = paste0(fpath, "/BEST_date.Rdata"))



# GLEAM actual evapotranspiration -----------------------------------------
#packages
library(rgdal)
library(sf)
library(raster)
library(ncdf4)
library(data.table)

fpath = "C:/Users/ls16959/Data/01_Climate/GLEAM_v3.3a/GLEAM_v33a_AET"
fname = list.files(fpath, full.names = T)[1]

load(file = "C:/Users/ls16959/Data/CAMELS/cat_name_shape.Rdata")

ALL_shp = readOGR(dsn = "C:/Users/ls16959/Data/CAMELS", layer = "ALL_shp")

ncin <- nc_open(fname)
missValue = ncatt_get(ncin,'E',"_FillValue")$value
nc_close(ncin)

temp_brick = brick(fname, varname = 'E')

s_time = Sys.time()
tempr = t(flip(temp_brick[[1]], direction = 'y'))
#need to change direction first before calculation
crs(tempr) = crs(ALL_shp)
catweights_GLEAM_camels = extract(tempr, ALL_shp,cellnumbers=TRUE, weights=TRUE, small = T)
e_time = Sys.time()
e_time-s_time
#9.7 min 2.1 min
save(catweights_GLEAM_camels, file = paste0(fpath, "/GLEAM_3_3a_E_catweights_camels.Rdata"))


catweights_nrow = unlist(lapply(catweights_GLEAM_camels, function(x){nrow(x)}))
catweights_GLEAM_df = do.call(rbind, catweights_GLEAM_camels)
tempweights = catweights_GLEAM_df[,3]
cat_id = rep(1:length(catweights_nrow), times = catweights_nrow)

GLEAM_flist = list.files(fpath, full.names = T, pattern = ".nc")
system.time({
  out_df_list = lapply(GLEAM_flist, function(fname){
    print(fname)
    temp_brick = brick(fname, varname = 'E')
    tempr_brick = t(flip(temp_brick, direction = 'y'))
    crs(tempr_brick) = crs(ALL_shp)
    select_cell_df = tempr_brick[catweights_GLEAM_df[,1]] 
    weighted_cells = data.table(sweep(select_cell_df, MARGIN=1, tempweights, '*'))
    out_df = weighted_cells[, lapply(.SD, sum), by=cat_id]
    return(out_df)
  })
})
#nearly seven hours
save(out_df_list, file = paste0(fpath, "/GLEAM_3_3a_camels_out_df_list.Rdata"))
GLEAMdf = do.call(rbind, lapply(out_df_list, function(df){t(as.data.frame(df)[,-1])}))
save(GLEAMdf, file = paste0(fpath, "/GLEAM_camels_df.Rdata"))


system.time({
  GLEAM_date_list = lapply(GLEAM_flist, function(fname){
    print(fname)
    temp_brick = brick(fname, varname = 'E')
    GLEAM_date = as.Date(substr(names(temp_brick), 2, 11), format = "%Y.%m.%d")
    return(GLEAM_date)
  })
})
GLEAM_date = do.call(c, GLEAM_date_list)
save(GLEAM_date, file = paste0(fpath, "/GLEAM_date.Rdata"))


