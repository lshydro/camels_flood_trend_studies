###########################################################-
# Objective: Download and pre-process all available CAMELS 
# datasets
# Author: Lina Stein, University of Bristol
# note: 
###########################################################-

#for each station create txt or csv file with Date, Q, QF (if available), naming convention for the file country-ISO_gaugenum_Q_m3s, i.e. "CL_1001001_streamflow_m3s.csv"
#For CAMELS_BR two quality flags columns are needed
#CL does not have quality flags, but some catchments should be excluded based on attributes (e.g. too much abstraction, dams...)
#check if all ts have full dates
library(measurements)

ALL_Q_path = "C:/Users/ls16959/Data/CAMELS/CAMELS_ALL_Q/"

# CAMELS-BR ---------------------------------------------------------------
fpath_BR = "C:/Users/ls16959/Data/CAMELS/CAMELS_BR/"
BR_dirlist = list.dirs(fpath_BR)
BRattributes_flist = list.files(BR_dirlist[2], full.names = T)
#more than 3000 gauge records, but only 897 with high quality and attributes
#Quality checks data gives list of high quality stations
BR_QC = read.delim(BRattributes_flist[8], sep = " ")

BR_streamflow_flist = list.files(BR_dirlist[7], full.names = T)
BR_streamflow_flist_short = list.files(BR_dirlist[7])
#create gauge no from file names
BR_gaugeno_all = unlist(lapply(strsplit(BR_streamflow_flist_short, "_"), function(x){x[1]}))
#Compare to list of quality checked catchments
BR_Q_files = BR_streamflow_flist[BR_gaugeno_all %in% BR_QC[,1]]


lapply(BR_QC[,1], function(x){
  cat_ind = grep(x, BR_Q_files)
  temp = read.delim(BR_Q_files[cat_ind], sep = " ")
  temp_date = as.Date(paste(temp[,1], temp[,2], temp[,3], sep = "-"), format = "%Y-%m-%d")
  #qual_control_by_ana should be 1, qual_flag should be 1 (potentially accept 2 as well)
  temp_Q = temp$streamflow_m3s
  temp_Q[temp_Q<0] = NA
  temp_Q[is.nan(temp_Q)] = NA
  #quality flags. Accept quality flag 1 and 2, turn Q value NA for values 0, 3, 4
  temp_Q[temp$qual_flag>2] = NA
  temp_Q[temp$qual_flag==0] = NA
  temp_Qdf = data.frame(Q_date = temp_date, Q_m3s = temp_Q)
  write.csv(temp_Qdf, file = paste0(ALL_Q_path, "BR_", x, "_streamflow_m3s.csv"), row.names = F)
})


# CAMELS-CL ---------------------------------------------------------------
fpath_CL = "C:/Users/ls16959/Data/CAMELS/CAMELS_CL/"
CL_flist = list.files(fpath_CL, full.names = T, pattern = ".txt")
CL_Q_file = list.files(fpath_CL, full.names = T, pattern = "streamflow_m3s.txt")
#all stations within one txt file
CL_Q = read.table(CL_Q_file, header = T)

CL_date = as.Date(CL_Q[,1], format = "%Y-%m-%d")
gauge_no = sub("X", "", colnames(CL_Q[,-1]))
CL_Qdf = CL_Q[,-1]
lapply(c(1:ncol(CL_Qdf)), function(x){
  temp_Q =  CL_Qdf[,x]
  temp_Q[temp_Q<0] = NA
  
  temp_Qdf = data.frame(Q_date = CL_date, Q_m3s = temp_Q)
  
  write.csv(temp_Qdf, file = paste0(ALL_Q_path, "CL_", gauge_no[x], "_streamflow_m3s.csv"), row.names = F)
})
#no quality flags to be included



# CAMELS-GB ---------------------------------------------------------------

fpath_GB = "C:/Users/ls16959/Data/CAMELS/CAMELS_GB/"
flistGB = list.files(fpath_GB, full.names = T)
flist_TS_GB = list.files(flistGB[grep("timeseries", flistGB)], full.names = T)

#Download time series data with "downthemall" from https://catalogue.ceh.ac.uk/datastore/eidchub/8344e4f3-d2ea-44f5-8afa-86d2987543a9/timeseries/

lapply(flist_TS_GB, function(xfile){
  gauge_no_temp = sub( "C:/Users/ls16959/Data/CAMELS/CAMELS_GB/timeseries/CAMELS_GB_hydromet_timeseries_", "", xfile)
  gauge_no = sub("_19701001-20150930.csv", "", gauge_no_temp)
  temp = read.csv(xfile)
  temp_Q = temp$discharge_vol
  temp_Q[temp_Q<0] = NA
  temp_date = as.Date(temp$date, format = "%Y-%m-%d")
  
  temp_Qdf = data.frame(Q_date = temp_date, Q_m3s = temp_Q)
  
  write.csv(temp_Qdf, file = paste0(ALL_Q_path, "GB_", gauge_no, "_streamflow_m3s.csv"), row.names = F)
})

# CAMELS-US ---------------------------------------------------------------

fpath_US = 'C:/Users/ls16959/Data/streamflow_countries/US/Camels/'
camels_US_topo = read.csv(file =  paste0(fpath_US, "camels_attributes_v2.0/camels_topo.txt"), sep = ";", stringsAsFactors = F, encoding = 'latin1')
#Model ouput files have climate and observed runoff all in  mm/d
model_output_qpath = list.dirs(paste0(fpath_US, "basin_timeseries_v1p2_modelOutput_daymet/model_output_daymet/model_output/flow_timeseries/daymet"))
flist_TS_US =unlist(lapply(model_output_qpath, list.files, full.name = T, pattern = "05_model_output"))

lapply(flist_TS_US, function(xfile){
  print(xfile)
  gauge_no_temp = strsplit(xfile, "/")[[1]][14]
  gauge_no = sub("_05_model_output.txt", "", gauge_no_temp)
  #read in observed streamflow from CAMELS data
  temp = read.table(xfile, quote="\"", comment.char="", header = T)
  temp_date = as.Date(paste(temp[,1], temp[,2], temp[,3], sep = "-"), format = "%Y-%m-%d")

  temp_area = camels_US_topo[grep(gauge_search_no, camels_US_loctopo$gauge_id),"area_geospa_fabric"]
  #divide by geospatial fabric area, since CAMELS article says that one was used to transfer to mm
  temp_Q = temp$OBS_RUN
  temp_Q[temp_Q<0] = NA
  #convert from mm to m3/s
  temp_Q = ((temp_Q/1000)*(temp_area*1000000))/(24*3600)
  
  temp_Qdf = data.frame(Q_date = temp_date, Q_m3s = temp_Q)
  write.csv(temp_Qdf, file = paste0(ALL_Q_path, "US_", gauge_no, "_streamflow_m3s.csv"), row.names = F)
})





# CAMELS-AUS --------------------------------------------------------------
library(readr)

fpath_AUS = "C:/Users/ls16959/Data/CAMELS/CAMELS_AUS/03_streamflow/"
flistAUS = list.files(fpath_AUS, full.names = T)
#all streamflow stations within one file
AUS_Q = read.csv(flistAUS[3], header = T, na.strings = "-99.99")

#quality flags. Accept quality flag A and B, turn Q value NA for values C, E, G
AUS_QF = data.frame(read_csv(flistAUS[5]), stringsAsFactors = F)
AUS_QF = AUS_QF[,-c(1:3)]

AUS_date = as.Date(paste(AUS_Q[,1], AUS_Q[,2], AUS_Q[,3], sep = "-"), format = "%Y-%m-%d")

AUS_Q_MLd = AUS_Q[,-c(1:3)]

lapply(c(1:ncol(AUS_Q_MLd)), function(i){
  gauge_no = substr(colnames(AUS_Q_MLd)[i],start = 2, stop = 100)
  temp_Q = AUS_Q[,i]
  temp_QF = AUS_QF[,i]
  
  temp_Q[temp_Q<0] = NA
  temp_Q[!(temp_QF %in% c("A", "B"))] = NA
  #convert from ML/day to cubic meter per second
  temp_Q = (temp_Q*1000)/(24*60*60)
  temp_Qdf = data.frame(Q_date = AUS_date, Q_m3s = temp_Q)
  write.csv(temp_Qdf, file = paste0(ALL_Q_path, "AUS_", gauge_no, "_streamflow_m3s.csv"), row.names = F)
})


# Australian Streamflow data ----------------------------------------------
fpath_AU = "C:/Users/ls16959/Data/streamflow_countries/Australia/Wasko_streamflow/data_baseflow/"
flist_TS_AU = list.files(fpath_AU, full.names = T)


lapply(flist_TS_AU, function(xfile){
  gauge_no_temp = sub( "C:/Users/ls16959/Data/streamflow_countries/Australia/Wasko_streamflow/data_baseflow/", "", xfile)
  gauge_no = sub(".txt", "", gauge_no_temp)
  #read in observed streamflow from CAMELS data
  temp = read.csv(xfile, sep="")
  temp_date = as.Date(temp$Date, format = "%Y-%m-%d")
  
  temp_Q = temp$Flow_ML
  #convert -999 to NA before conversion, otherwise will be converted as well
  temp_Q[temp_Q<0] = NA
  #quality flags. Accept quality flag A and B, turn Q value NA for values C, E, G
  temp_Q[!(temp$Bureau.QCode %in% c("A", "B"))] = NA
  #convert from ML/day to cubic meter per second
  temp_Q = (temp_Q*1000)/(24*60*60)
  temp_Qdf = data.frame(Q_date = temp_date, Q_m3s = temp_Q)
  
  write.csv(temp_Qdf, file = paste0(ALL_Q_path, "AU_", gauge_no, "_streamflow_m3s.csv"), row.names = F)
})


#remove AU files, since CAMELS_AUS was published
AU_files = list.files(ALL_Q_path, pattern =  "AU_", full.names = T)
#file.remove(AU_files)

# Quality check streamflow ------------------------------------------------
#for each country read all streamflow records into a list
#merge to complete dates (sequence from 01.01. of first year to 31.12 of last year)
#check if sufficient records remain

####Not needed anymore#### included in POT calculation

#GSIM: yearly index is reliable if at least 350 values are available
library(chron)
library(data.table)


qfiles = list.files(ALL_Q_path, full.names = T, pattern = ".csv")
s_time = Sys.time()
QC_list = lapply(qfiles, function(xfile){
  print(xfile)
  gauge_no = sub("C:/Users/ls16959/Data/CAMELS/CAMELS_ALL_Q/", "", xfile)
  gauge_no = sub("_streamflow_m3s.csv", "", gauge_no)
  temp = read.csv(xfile, header = T, stringsAsFactors = F)
  temp[,1] = as.Date(temp[,1])
  start_date = as.Date(paste0(years(temp[1,1]), "-01-01"))
  end_date = as.Date(paste0(years(temp[nrow(temp),1]), "-12-31"))
  date_seq = seq(start_date, end_date, by = "days")
  date_seq_df = data.frame(Q_date = date_seq, Q_years = years(date_seq))
  tempdf = data.table(merge(date_seq_df, temp, all.x = T))
  
  countNAdf = data.frame(tempdf[, .(sumNA = sum(is.na(.SD))), by = Q_years])
  no_years = nrow(countNAdf)
  no_years_enough = sum(countNAdf[,2]<15)
  years_enough = countNAdf[,1][countNAdf[,2]<15]
  return(list(gauge_no, no_years, no_years_enough, years_enough))
})
e_time = Sys.time()
e_time-s_time
save(QC_list, file = "C:/Users/ls16959/Data/CAMELS/QC_list.Rdata")
QC_df = do.call(rbind, lapply(QC_list, function(x){x[1:3]}))
sum(QC_df[,3]>=20)

load(file = "C:/Users/ls16959/Data/CAMELS/QC_df.Rdata")
table(substr(QC_df[,1][QC_df[,3]>=20], 1,2))


# Combine catchment attributes --------------------------------------------
#combine all climate attributes+elevation + any exclusion criteria (dams)

snow_frac_thresh = 0.15
snow_frac_thresh = 999 #to not make snow fraction a deciding factor
#%%%%%%%%%%%%%%%
#Brazil
#%%%%%%%%%%%%%%%
fpath_BR = "C:/Users/ls16959/Data/CAMELS/CAMELS_BR/"
BR_dirlist = list.dirs(fpath_BR)
BRattributes_flist = list.files(BR_dirlist[2], full.names = T)
camels_BR_location = read.delim(BRattributes_flist[7], sep = " ")
camels_BR_climate = read.delim(BRattributes_flist[2], sep = " ")
camels_BR_human = read.delim(BRattributes_flist[4], sep = " ")
camels_BR_topo = read.delim(BRattributes_flist[10], sep = " ")

#Climate
camels_BR_gauge_no = paste0("BR_", camels_BR_climate[,"gauge_id"])
camels_BR_climate[,"gauge_id"] = camels_BR_gauge_no
#climate attributes same as camels_US

#Quality
BR_reservoir = camels_BR_human$vol_reservoirs
BR_reservoir_infl = camels_BR_human$degree_of_regulation
BR_QF = data.frame(gauge_id = camels_BR_gauge_no, QF = ifelse(BR_reservoir_infl<10 & camels_BR_climate$frac_snow<snow_frac_thresh, T, F),stringsAsFactors = F)
#degree of regulation generally small for the selected 800sth catchments already. Only one catchment with very high degree of regulation. 10% is inspired by camels BR article talking about "% of the catchments have a degree of regulation greater than 10 % "

#BR_46455000 has a jump in the timeline, need to exclude from calculations
BR_QF[grep("BR_46455000", BR_QF$gauge_id),2] = F


#Loco-topo
BR_elevation = camels_BR_topo$elev_mean

camels_BR_loctopo = merge(camels_BR_location[,c("gauge_id", "gauge_lat", "gauge_lon")], camels_BR_topo[,c("gauge_id", "elev_mean", "slope_mean", "area")])
camels_BR_loctopo[,"gauge_id"] = camels_BR_gauge_no

#%%%%%%%%%%%%%%%
#Chile
#%%%%%%%%%%%%%%%
fpath_CL = "C:/Users/ls16959/Data/CAMELS/CAMELS_CL/"

CL_attributes_file = list.files(fpath_CL, full.names = T, pattern = "attributes.txt")
CL_attributes = read.table(CL_attributes_file, header = F, row.names = 1,colClasses = "character")
#turn attributes to data frame
CL_attributes = data.frame(t(CL_attributes), stringsAsFactors = F)

#CL_gauge id
CL_gauge_no = sub(" ", "", paste0("CL_", CL_attributes[,"gauge_id"]))


#Climate
camels_CL_clim = data.frame(gauge_id = CL_gauge_no, p_mean = CL_attributes$p_mean_mswep, pet_mean = CL_attributes$pet_mean, p_seasonality = CL_attributes$p_seasonality_mswep, frac_snow = CL_attributes$frac_snow_mswep, aridity = CL_attributes$aridity_mswep, high_prec_freq = CL_attributes$high_prec_freq_mswep, high_prec_dur = CL_attributes$high_prec_dur_mswep, high_prec_timing = CL_attributes$low_prec_timing_mswep, low_prec_freq = CL_attributes$low_prec_freq_mswep, low_prec_dur = CL_attributes$low_prec_dur_mswep, low_prec_timing = CL_attributes$low_prec_timing_mswep)
#change all factors to characters
camels_CL_clim = data.frame(lapply(camels_CL_clim, as.character), stringsAsFactors=FALSE)
#change all characterst to numeric (where possible)
camels_CL_clim[,1:ncol(camels_CL_clim)]=lapply(1:ncol(camels_CL_clim),function(x) {
  tryCatch({
    as.numeric(camels_CL_clim[[x]])
  },warning = function(w) {
    camels_CL_clim[[x]]}
  )} )

#Quality
CL_big_dam = as.numeric(unlist(CL_attributes[,which(colnames(CL_attributes)=="big_dam")]))
CL_glacier_frac = as.numeric(unlist(CL_attributes[,which(colnames(CL_attributes)=="lc_glacier")]))
Cl_quality_df = data.frame(gauge_id = CL_gauge_no, dams = CL_big_dam, glacier_perc = CL_glacier_frac, snow_frac = camels_CL_clim$frac_snow)
CL_QF = data.frame(gauge_id = CL_gauge_no, QF = ifelse(CL_big_dam == 0 & CL_glacier_frac<5 & camels_CL_clim$frac_snow<snow_frac_thresh, T, F),stringsAsFactors = F)


#Loco-Topo
camels_CL_loctopo = data.frame(gauge_id = CL_gauge_no, gauge_lat = CL_attributes$gauge_lat, gauge_lon = CL_attributes$gauge_lon, elev_mean = CL_attributes$elev_mean, slope_mean = CL_attributes$slope_mean, area = CL_attributes$area)
#change all factors to characters
camels_CL_loctopo = data.frame(lapply(camels_CL_loctopo, as.character), stringsAsFactors=FALSE)
#change all characterst to numeric (where possible)
camels_CL_loctopo[,1:ncol(camels_CL_loctopo)]=lapply(1:ncol(camels_CL_loctopo),function(x) {
  tryCatch({
    as.numeric(camels_CL_loctopo[[x]])
  },warning = function(w) {
    camels_CL_loctopo[[x]]}
)} )




#%%%%%%%%%%%%%%%
#GB
#%%%%%%%%%%%%%%%
fpath_GB = "C:/Users/ls16959/Data/CAMELS/CAMELS_GB/"
flistGB = list.files(fpath_GB, full.names = T)
camels_GB_clim = read.csv(grep("climatic_", flistGB, value = T), stringsAsFactors = F)
camels_GB_topo = read.csv(grep("topographic_", flistGB, value = T), stringsAsFactors = F)
camels_GB_human = read.csv(grep("humaninfluence_", flistGB, value = T), stringsAsFactors = F)
camels_GB_hydrometry = read.csv(grep("hydrometry_", flistGB, value = T), stringsAsFactors = F)

#Climate
camels_GB_gauge_no = paste0("GB_", camels_GB_clim[,"gauge_id"])
camels_GB_clim[,"gauge_id"] = camels_GB_gauge_no
#climate attributes same as camels_US


#Quality
GB_QF = data.frame(gauge_id = camels_GB_gauge_no, QF = ifelse(camels_GB_human$num_reservoir == 0 & camels_GB_clim$frac_snow<snow_frac_thresh, T, F),stringsAsFactors = F)
#manually detected error in catchment
GB_QF[which(GB_QF$gauge_id == "GB_54052"),2] = F


#Loco-Topo
camels_GB_loctopo = camels_GB_topo[,c("gauge_id", "gauge_lat", "gauge_lon", "elev_mean","dpsbar", "area")]
colnames(camels_GB_loctopo) = colnames(camels_BR_loctopo) #match columns to BR, use same names
camels_GB_loctopo[,"gauge_id"] = camels_GB_gauge_no

#%%%%%%%%%%%%%%%
#US
#%%%%%%%%%%%%%%%
fpath_US = 'C:/Users/ls16959/Data/streamflow_countries/US/Camels/'
camels_US_name = read.csv(file = paste0(fpath_US, "camels_attributes_v2.0/camels_name.txt"), sep = ";", stringsAsFactors = F, encoding = 'latin1')
camels_US_topo = read.csv(file =  paste0(fpath_US, "camels_attributes_v2.0/camels_topo.txt"), sep = ";", stringsAsFactors = F, encoding = 'latin1')
camels_US_clim = read.csv(file =  paste0(fpath_US, "camels_attributes_v2.0/camels_clim.txt"), sep = ";", stringsAsFactors = F, encoding = 'latin1')

US_gauge_no = paste0("US_", camels_US_name[,"gauge_id"])
#Climate
camels_US_clim[,"gauge_id"] = US_gauge_no

#Quality

camels_US_areaerror = abs(1-camels_US_topo$area_gages2/camels_US_topo$area_geospa_fabric)
US_QF = data.frame(gauge_id = US_gauge_no, QF = ifelse(camels_US_clim$frac_snow<snow_frac_thresh & camels_US_areaerror<0.1, T, F),stringsAsFactors = F)


#Loco-Topo
camels_US_loctopo = camels_US_topo[,c(1:5, 7)] #use geospatial area since that was used to convert m3/s into mm
colnames(camels_US_loctopo) = colnames(camels_BR_loctopo)
camels_US_loctopo[,"gauge_id"] = US_gauge_no


#%%%%%%%%%%%%%%%
#Australia
#%%%%%%%%%%%%%%%
#Climate
camels_AUS_clim = data.frame(read_csv("C:/Users/ls16959/Data/CAMELS/CAMELS_AUS/05_hydrometeorology/ClimaticIndices.csv"))
head(camels_AUS_clim)

camels_AUS_gauge_no = paste0("AUS_", camels_AUS_clim[,"ID"])
camels_AUS_clim[,"gauge_id"] = camels_AUS_gauge_no
camels_AUS_clim = camels_AUS_clim[,-1] #remove ID column

#Quality
#In the catchment metadata file, some catchments have reported uncertainty in catchment boundary or other notes. Flag those. 
camels_AUS_metadata = read.csv("C:/Users/ls16959/Data/CAMELS/CAMELS_AUS/01_id_name_metadata/id_name_metadata.csv")

camels_AUS_gauge_no = paste0("AUS_", camels_AUS_metadata$station_id)
AUS_QF =  data.frame(gauge_id = camels_AUS_gauge_no, QF = camels_AUS_metadata$notes=="No notes",stringsAsFactors = F)



#Loco-Topo
camels_AUS_loco = data.frame(read_csv("C:/Users/ls16959/Data/CAMELS/CAMELS_AUS/02_location_boundary_area/location_boundary_area.csv"))
camels_AUS_gauge_no = paste0("AUS_", camels_AUS_loco$station_id)
camels_AUS_loco = data.frame(gauge_id = camels_AUS_gauge_no, gauge_lat = camels_AUS_loco$lat_outlet, gauge_lon = camels_AUS_loco$long_outlet, area = camels_AUS_loco$catchment_area)

camels_AUS_topo = data.frame(read_csv("C:/Users/ls16959/Data/CAMELS/CAMELS_AUS/04_attributes/CatchmentAttributes_02_Topography&Geometry.csv"))
camels_AUS_gauge_no = paste0("AUS_", camels_AUS_topo$station_id)
camels_AUS_topo = data.frame(gauge_id = camels_AUS_gauge_no, elev_mean= camels_AUS_topo$elev_mean, slope_mean = camels_AUS_topo$mean_slope_pct)

camels_AUS_locotopo = merge(camels_AUS_loco, camels_AUS_topo, by = "gauge_id")


#AU old----
fpath_AU = "C:/Users/ls16959/Data/streamflow_countries/Australia/Wasko_streamflow/data_baseflow/"
flist_TS_AU = list.files(fpath_AU, full.names = T)
camels_AU_gauge_no = sub("C:/Users/ls16959/Data/streamflow_countries/Australia/Wasko_streamflow/data_baseflow/", "", flist_TS_AU)
camels_AU_gauge_no = sub(".txt", "", camels_AU_gauge_no)
camels_AU_gauge_no = paste0("AU_", camels_AU_gauge_no)
AU_QF =  data.frame(gauge_id = camels_AU_gauge_no, QF = rep(T, length(camels_AU_gauge_no)),stringsAsFactors = F)


library(readxl)
hrs_station_details <- read_excel("C:/Users/ls16959/Data/streamflow_countries/Australia/Wasko_cat_boundaries/hrs_station_details.xls")


AU_gauge_no = sub("C:/Users/ls16959/Data/streamflow_countries/Australia/Wasko_streamflow/data_baseflow/", "", flist_TS_AU)
AU_gauge_no = sub(".txt", "", AU_gauge_no)

head(hrs_station_details)

AU_loctopo = data.frame(gauge_id = hrs_station_details$`AWRC Station Number`, gauge_lat = hrs_station_details$Latitude, gauge_lon = hrs_station_details$Longitude, elev_mean = rep(NA, nrow(hrs_station_details)), slope_mean = rep(NA, nrow(hrs_station_details)), area = hrs_station_details$`Catchment Area (km2)`)
AU_loctopo = merge(data.frame(gauge_id = AU_gauge_no), AU_loctopo, all.x = T)
AU_loctopo[,"gauge_id"] = camels_AU_gauge_no


#%%%%%%%%%%%%%%%
#All----
#%%%%%%%%%%%%%%%
#Climate
camels_all_clim = do.call(rbind, list(camels_BR_climate, camels_CL_clim, camels_GB_clim, camels_US_clim, camels_AUS_clim))
save(camels_all_clim, file = "C:/Users/ls16959/Data/CAMELS/camels_all_clim.Rdata")


load(file = "C:/Users/ls16959/Data/CAMELS/camels_all_clim.Rdata")
head(camels_all_clim)
library(ggplot2)
camels_clim_test = camels_all_clim
camels_clim_test[,"country"] = substr(camels_clim_test$gauge_id, 1,2)
ggplot(camels_clim_test, aes(aridity, p_seasonality, col = country))+geom_point()+xlim(0,5)


#Quality
All_QF = do.call(rbind, list(BR_QF, CL_QF, GB_QF, US_QF, AUS_QF))
All_QF[is.na(All_QF$QF),"QF"] = F
#save(All_QF, file = "C:/Users/ls16959/Data/CAMELS/All_QF.Rdata")

save(All_QF, file = "C:/Users/ls16959/Data/CAMELS/All_QF_nosnowrestr.Rdata")

#Loco-Topo
library(data.table)
camels_all_loctopo = rbindlist(list(camels_BR_loctopo, camels_CL_loctopo, camels_GB_loctopo, camels_US_loctopo, camels_AUS_locotopo), use.names = T)
save(camels_all_loctopo,  file = "C:/Users/ls16959/Data/CAMELS/camels_all_loctopo.Rdata")


# Combine catchment shapefiles --------------------------------------------
library(rgdal)
library(raster)
#bug in  original catchment shapefile, this is the updated version (supplied by authors)
BR_shp = readOGR(dsn = "C:/Users/ls16959/Data/CAMELS/CAMELS_BR", layer = "camels_br_catchments_fixed")
BR_shpdf = data.frame(BR_shp)
BR_shp@data$country <- rep("BR", nrow(BR_shpdf))
BR_shpdf = data.frame(BR_shp)

CL_shp = readOGR(dsn = "C:/Users/ls16959/Data/CAMELS/CAMELS_CL/CAMELScl_catchment_boundaries", layer = "catchments_camels_cl_v1.3")
CL_shpdf = data.frame(CL_shp)
CL_shp@data$country <- rep("CL", nrow(CL_shpdf))
CL_shpdf = data.frame(CL_shp)

GB_shp = readOGR(dsn = "C:/Users/ls16959/Data/CAMELS/CAMELS_GB", layer = "CAMELS_GB_catchment_boundaries")
GB_shpdf = data.frame(GB_shp)
names(GB_shp@data)[names(GB_shp@data)=="ID_STRING"] <- "gauge_id"
GB_shp@data$country <- rep("GB", nrow(GB_shpdf))
GB_shpdf = data.frame(GB_shp)

US_shp = readOGR(dsn = "C:/Users/ls16959/Data/streamflow_countries/US/Camels/basin_set_full_res", layer = "HCDN_nhru_final_671")
names(US_shp@data)[names(US_shp@data)=="hru_id"] <- "gauge_id"
US_shpdf = data.frame(US_shp)
US_shp@data$country <- rep("US", nrow(US_shpdf))
US_shpdf = data.frame(US_shp)


AUS_shp = readOGR(dsn = "C:/Users/ls16959/Data/CAMELS/CAMELS_AUS/02_location_boundary_area/shp", layer = "CAMELS_AUS_Boundaries_adopted")
names(AUS_shp@data)[names(AUS_shp@data)=="CatchID"] <- "gauge_id"
AUS_shpdf = data.frame(AUS_shp)
AUS_shp@data$country <- rep("AUS", nrow(AUS_shpdf))
AUS_shpdf = data.frame(AUS_shp)

#old
# AU_shp = readOGR(dsn = "C:/Users/ls16959/Data/streamflow_countries/Australia/Wasko_cat_boundaries/GIS/HRS_Boundaries_fromSRTM_v0.1_20140326", layer = "HRS_Boundaries_fromSRTM")
# names(AU_shp@data)[names(AU_shp@data)=="CatchID"] <- "gauge_id"
# AU_shpdf = data.frame(AU_shp)
# AU_shp@data$country <- rep("AU", nrow(AU_shpdf))
# AU_shpdf = data.frame(AU_shp)

shp_list = list(BR_shp, CL_shp, GB_shp, US_shp, AUS_shp)
#GB and US have different reference system
lapply(shp_list, crs)
standard_crs = crs(BR_shp)
GB_shp = spTransform(GB_shp, standard_crs)
US_shp = spTransform(US_shp, standard_crs)
#all shapefiles have same projection now
ALL_shp = union(BR_shp, CL_shp) 
ALL_shp = union(ALL_shp, GB_shp) 
ALL_shp = union(ALL_shp, US_shp) 
ALL_shp = union(ALL_shp, AUS_shp)
plot(ALL_shp)
writeOGR(ALL_shp, dsn = "C:/Users/ls16959/Data/CAMELS" , layer = "ALL_shp", driver="ESRI Shapefile")

rm(ALL_shp)

ALL_shp = readOGR(dsn = "C:/Users/ls16959/Data/CAMELS", layer = "ALL_shp")
plot(ALL_shp)
ALL_shp_df = data.frame(ALL_shp)



#are brazil catchments in there twice?
#do shape data.frame transfer?
head(BR_shpdf)
head(CL_shpdf)
head(GB_shpdf)
head(US_shpdf)
head(AU_shpdf)

dim(AU_shpdf)
length(flist_TS_AU)


# Test shapefile quality --------------------------------------------------


ALL_shp$area_sqkm <- conv_unit(area(ALL_shp), "m2", "km2")
ALL_shp_df = data.frame(ALL_shp, stringsAsFactors = F)
ALL_shp_df[,"gauge_id"] = as.character(ALL_shp_df$gauge_id)

# Brazil ####
camels_br_topography <- read.csv("C:/Users/ls16959/Data/CAMELS/CAMELS_BR/1_CAMELS_BR_attributes/camels_br_topography.txt", sep="")

BR_all_df = ALL_shp_df[which(ALL_shp_df$country == "BR"),]
BR_merge_df = merge(BR_all_df, camels_br_topography, by = "gauge_id", all = T)
hist(BR_merge_df$area_sqkm/BR_merge_df$area)
#less than 5% difference

#previous version shapefile identified erroneous catchments
BR_shp_error = BR_shpdf[!(BR_shpdf$gauge_id %in% BR_gaugeno_all),]
write.csv(BR_shp_error, file = "C:/Users/ls16959/Data/CAMELS/clarify_gauge_id.csv")


# CL ####
#IDs all match
sum(CL_shpdf$gauge_id %in% gauge_no) #gauge_no from streamflow reading part
sum(gauge_no %in% CL_shpdf$gauge_id)

#check catchment quality
#cannot merge data together but need to trust they are in the right order?
CL_all_df = ALL_shp_df[which(ALL_shp_df$country == "CL"),]
CL_attr_area = as.numeric(CL_attributes[,which(CL_attributes[1,]=="area")][-1])
range(CL_all_df$area_sqkm/CL_attr_area)
#order seems to be the same, since difference is very little



# GB ####
CAMELS_GB_topographic_attributes <- read.csv("C:/Users/ls16959/Data/CAMELS/CAMELS_GB/CAMELS_GB_topographic_attributes.csv", stringsAsFactors=FALSE)
#IDs all match
sum(GB_shpdf$gauge_id %in% CAMELS_GB_topographic_attributes$gauge_id)
sum(CAMELS_GB_topographic_attributes$gauge_id %in% GB_shpdf$gauge_id)

#check catchment quality

GB_all_df = ALL_shp_df[which(ALL_shp_df$country == "GB"),]
GB_merge_df = merge(GB_all_df, CAMELS_GB_topographic_attributes, by = "gauge_id", all = T)
#all catchments very little difference between calculated and actual area (might be due to coordindate system)
hist(GB_merge_df$area_sqkm/GB_merge_df$area)


# US ####
camels_US_topo <- read.csv("C:/Users/ls16959/Data/streamflow_countries/US/Camels/camels_attributes_v2.0/camels_topo.txt", sep=";")
#IDs all match
sum(US_shpdf$gauge_id %in% camels_US_topo$gauge_id)
sum(camels_US_topo$gauge_id %in% US_shpdf$gauge_id)

#check catchment quality

US_all_df = ALL_shp_df[which(ALL_shp_df$country == "US"),]
US_merge_df = merge(US_all_df, camels_US_topo, by = "gauge_id", all = T)
#all catchments very little difference between calculated and actual area (might be due to coordindate system)
hist(US_merge_df$area_sqkm/US_merge_df$area_gages2)
#geospatial fabric fits better
hist(US_merge_df$area_sqkm/US_merge_df$area_geospa_fabric)


# AU ####
#
library(readr)
camels_AUS_loco = data.frame(read_csv("C:/Users/ls16959/Data/CAMELS/CAMELS_AUS/02_location_boundary_area/location_boundary_area.csv"))
camels_AUS_gauge_no = paste0("AUS_", camels_AUS_loco$station_id)
camels_AUS_loco = data.frame(gauge_id = camels_AUS_gauge_no, gauge_lat = camels_AUS_loco$lat_outlet, gauge_lon = camels_AUS_loco$long_outlet, area = camels_AUS_loco$catchment_area)

AUS_all_df = ALL_shp_df[which(ALL_shp_df$country == "AUS"),]
AUS_all_df[,"gauge_id"] = paste("AUS", AUS_all_df$gauge_id, sep = "_")
AUS_merge_df = merge(AUS_all_df, camels_AUS_loco, by = "gauge_id", all = T)

hist(AUS_merge_df$area_sqkm/AUS_merge_df$area)
#difference is 5% max

