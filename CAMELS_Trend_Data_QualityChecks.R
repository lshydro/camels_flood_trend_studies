###########################################################-
# Objective: Prepare flood event data for all CAMELS 
# catchments for trend analysis
# Author: Lina Stein, University of Bristol
# note: 
###########################################################-



# Packages ----------------------------------------------------------------


library(quantmod)
library(parallel)
library(mblm)

library(lubridate)

# Functions ---------------------------------------------------------------



#Find peaks over threshold
POT_func = function(qobsdf, POT_thresh, eventsperyear){
  #1. set low threshold (median?)
  #2. Find Peaks
  #3. remove all peaks less than threshold
  #check independence
  #4. "the flow between them dropped to less than two thirds of the earlier of the two" Cunnane, 1979
  #5. time between peaks needs to be higher than 3*mean rising time (Cunnane, 1979)
  ###Input
  #qobsdf: discharge timeseries with date and flow
  #POT_thresh: threshold for POT analysis
  #eventsperyear: limit for in average how many events occur per year
  
  
  Qobsdf = qobsdf
  colnames(Qobsdf) = c("Qobs_date", "Qobs")
  peak_vec = findPeaks(Qobsdf[,2])-1 #is always off by one
  #find peaks with a large difference to next peak and use those for mean rising limb (but use single peak events)
  #identify difference between peaks
  peak_dif_vec = rep(NA, length(peak_vec))
  peak_dif_vec[1] = 0
  for(i in c(2:length(peak_vec))){
    peak_dif_vec[i] = peak_vec[i]-peak_vec[i-1]
  }
  temp = data.frame(peak_vec, peak_dif_vec, Qobsdf[peak_vec,])
  #order dataframe by Q magnitude
  temp = temp[order(temp$Qobs, decreasing = T),]
  #select 5 highest events that have a peak rising time over longer than the 75% percentile of peak rising time (min 5 days)
  temp5 = temp[temp$peak_dif_vec>max(quantile(peak_dif_vec)[4], 5),][1:5,]
  #find valley points
  valley_vec = findValleys(Qobsdf[,2])-1 #is always off by one
  #calculate mean difference between 10 peaks and the valley right before
  mean_rising_time = mean(unlist(lapply(c(1:5), FUN = function(x){
    peakp = temp5$peak_vec[x]
    valleyp = tail(valley_vec[which(valley_vec<peakp)],1)
    peakp - valleyp
  })))
  #difference between peaks needs to be higher than 3*mean rising time
  #highest peak is independent, two closest peaks need to be checked
  
  #select all POT above treshold, order by magnitude (highest first), then take that one and remove all peaks that are within the mean rising time
  #3. remove all peaks less than threshold
  POTdf = temp[temp$Qobs>POT_thresh,]
  
  #4. "the flow between them dropped to less than two thirds of the earlier of the two" Cunnane, 1979
  POTdf = POTdf[order(POTdf$Qobs_date),]
  i = 1
  while((i+1)<=dim(POTdf)[1]){
    qP1 = POTdf$Qobs[i]
    qP2 = POTdf$Qobs[i+1]
    P1_ind = POTdf$peak_vec[i] #get all days between peak days
    P2_ind = POTdf$peak_vec[i+1]
    inbetween = Qobsdf[c(P1_ind:P2_ind),"Qobs"]
    inbetween = inbetween[!is.na(inbetween)]
    #if inbetween Q is less than 2/3 of the first peak. Use while loop to keep checking the next peek, if a peek got eliminated
    #2/3 of the smaller peak: https://www.tandfonline.com/doi/full/10.1080/02626667.2018.1444766 mentions Mediero et al, 2015
    while(sum(inbetween<((2/3)*min(qP1, qP2)))==0){
      if(qP1<qP2){
        POTdf = POTdf[-(i),]
      }else{
        POTdf = POTdf[-(i+1),]
      }
      if((i+1)>dim(POTdf)[1]){
        break #need break conditions so that indexing does not go over the dimensions of the data frame
      }
      qP1 = POTdf$Qobs[i]
      qP2 = POTdf$Qobs[i+1]
      P1_ind = POTdf$peak_vec[i] #get all days between peak days
      P2_ind = POTdf$peak_vec[i+1]
      inbetween = Qobsdf[c(P1_ind:P2_ind),"Qobs"]
      inbetween = inbetween[!is.na(inbetween)]
    }
    i = i+1
  }
  
  #5. remove all peaks within peak rising time of a bigger peak
  #while loop (while tests if i is still bigger than the dimensions of the data frame)
  #sort by magnitude
  POTdf = POTdf[order(POTdf$Qobs, decreasing = T),]
  i = 1
  while(i<=dim(POTdf)[1]){
    #create time vector to catch all dates in +- rising time window of biggest peak
    date_vec = seq(POTdf[i,3]-(3*mean_rising_time), POTdf[i,3]+(3*mean_rising_time), by ="days")
    #which events are too close to the max peak (which is always the first)
    all_ind = as.character(POTdf[,3]) %in% as.character(date_vec) #change to character, as otherwise does not work
    #which one is the maximum peak
    keep_ind = which.max(as.character(POTdf[,3]) %in% as.character(date_vec)) #change to character, as otherwise does not work
    all_ind[keep_ind] = FALSE #reverse the maximum to keep this event
    #remove all other dependent events
    POTdf = POTdf[!all_ind,]
    i = i+1
  }
  POTdf = POTdf[,3:4]
  #events are ordered by magnitude, to get x events per year, calculate number of years and cut df accordingly 
  year_num = length(unique(chron::years(Qobsdf[!is.na(Qobsdf[,2]),1])))
  POTdfy = POTdf[1:min((year_num*eventsperyear), length(POTdf[,1])),] #if there are not enough events, limit is set to maximum number of events
  return(POTdfy)
}




theil_func = function(date_vals, trend_vals){
  if(sum(!is.na(trend_vals))<20){
    sensslope = NA
    senspval = NA
    sensslope_10 = NA
  } else{
    TS_slope = mblm(trend_vals~date_vals)
    TS_summary = summary(TS_slope)
    sensslope = TS_summary$coefficients[2,1]
    sensslope_10 = ((sensslope*10*365)/mean(trend_vals))*100
    senspval = TS_summary$coefficients[2,4]
  }
  return(c(sensslope, sensslope_10, senspval))
}




# Prepare files -----------------------------------------------------------

ALL_Q_path = "C:/Users/ls16959/Data/CAMELS/CAMELS_ALL_Q/"


fpath_out = "C:/Users/ls16959/Data/01_Climate/ERA5-land/CAMELS"
load(file = paste0(fpath_out, "/ERA5_date.Rdata"))
load(file = paste0(fpath_out, "/ERA5_SM_cat_df.Rdata"))

fpath = "C:/Users/ls16959/Data/01_Climate/MSWEP_V2.2_010deg_20180326/"
load(file = paste0(fpath, "/MSWEP2_2_camels_out_df_list.Rdata"))
MSWEPdf = do.call(rbind, lapply(out_df_list, function(df){t(as.data.frame(df)[,-1])}))
MSWEP_date = as.Date(sub("X","", rownames(MSWEPdf)), format = "%Y.%m.%d")

range(MSWEP_date)
range(ERA5_date)
start_date = max(c(min(MSWEP_date), min(ERA5_date)))
end_date = min(c(max(MSWEP_date), max(ERA5_date)))
date_seq = seq(start_date, end_date, by = "days")
date_seq_df = data.frame(Q_date = date_seq, Q_years = chron::years(date_seq))

#which catchments are included in the P and SM data:
load(file = "C:/Users/ls16959/Data/CAMELS/cat_name_shape.Rdata")



# Calculate POT -----------------------------------------------------------

qfiles = list.files(ALL_Q_path, full.names = T, pattern = ".csv")

qfiles_name = sub("C:/Users/ls16959/Data/CAMELS/CAMELS_ALL_Q/", "", qfiles)
qfiles_name = sub("_streamflow_m3s.csv", "", qfiles_name)
#US catchment shapefile numbers are missing zeros at the start of the number (but not for all). Therefore remove all leading zeros, to make match with catchment gauge id from shapefiles
qfiles_name_sub = sub("^US_0", "US_", qfiles_name)
length(qfiles_name_sub %in% cat_name_shape)


#What to do with years with insufficient data? Exclude completely by setting all to NA?
#Include for calculating threshold (mean Q)? If frequency is to be counted per year, all years with less than 350 values need to be excluded.  #GSIM: yearly index is reliable if at least 350 values are available
#What about missing values in the dry season?
#merge to complete dates (sequence from 01.01. of first year to 31.12 of last year)
#check if sufficient records remain

#For faster reading, resave all CVS files as Rdata files


s_time = Sys.time()
lapply(qfiles[qfiles_name_sub %in% cat_name_shape], function(xfile){
  temp = read.csv(xfile, header = T, stringsAsFactors = F)
  save(temp, file = sub(".csv", ".Rdata", xfile))
})
e_time = Sys.time()
e_time-s_time
#2min

qfilesR = list.files(ALL_Q_path, full.names = T, pattern = ".Rdata")

s_time = Sys.time()
# Calculate the number of cores
no_cores <- detectCores() - 3
# Initiate cluster
cl <- makeCluster(no_cores)
clusterExport(cl, varlist = c("qfilesR", "POT_func", "data.table", "date_seq_df",  "findPeaks", "findValleys", "years"))
QC_POT_list = parLapply(cl, qfilesR, function(xfile){
  gauge_no = sub("C:/Users/ls16959/Data/CAMELS/CAMELS_ALL_Q/", "", xfile)
  gauge_no = sub("_streamflow_m3s.Rdata", "", gauge_no)
  load(xfile)
  temp[,1] = as.Date(temp[,1])
  tempdf = merge(date_seq_df, temp, all.x = T)
  #check for continous (zero) values 
  #identify repeat values with more than a year no change (long period to not exclude zero values in arid catchments)
  test = rle(tempdf$Q_m3s)
  ident_test = rep(c(1:length(test$lengths)), times = test$lengths)
  tempdf[which(ident_test == which(test$lengths>365)),"Q_m3s"] = NA

  #Check if sufficient years are supplied
  tempdf_QC = data.table(tempdf)
  countNAdf = data.frame(tempdf_QC[, .(sumNA = sum(is.na(.SD))), by = Q_years])
  no_years = nrow(countNAdf)
  no_years_enough = sum(countNAdf[,2]<15)
  years_enough = countNAdf[,1][countNAdf[,2]<15]
  
  tempdf_POT = tempdf[,-2]
  tempdf_POT[!(tempdf$Q_years %in% years_enough),2] = NA #tunr all values to NA where not sufficient values per year are available
  
  if(no_years_enough>=20){
    #calculate Peaks-over-threshold
    POT_thresh = mean(tempdf_POT[,2], na.rm = T)
    POT_temp_list = lapply(c(1, 3, 10), FUN = function(x){POT_func(tempdf_POT, POT_thresh, eventsperyear = x)})
  } else{
    POT_temp_list = vector(mode = "list", length = 3)
  }
  return(list(gauge_no, no_years, no_years_enough, years_enough, POT_temp_list))
})
stopCluster(cl)
e_time = Sys.time()
e_time-s_time
#41 min

save(QC_POT_list, file = "C:/Users/ls16959/Data/CAMELS/QC_POT_list.Rdata")
load(file = "C:/Users/ls16959/Data/CAMELS/QC_POT_list.Rdata")
QC_df = do.call(rbind, lapply(QC_POT_list, function(x){x[1:3]}))
sum(QC_df[,3]>=20)


POT_names = unlist(lapply(QC_POT_list, function(x){x[[1]]}))
POT_names = sub("^US_0", "US_", POT_names)

