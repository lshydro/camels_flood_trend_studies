###########################################################-
# Objective: P-SM-Q Interactions and flood trends
# Author: Lina Stein, University of Bristol
# note: 
###########################################################-


# Packages ----------------------------------------------------------------

library(zoo)
library(ggplot2)
library(patchwork)
library(reshape)
library(tidyverse)
library(viridis)
library(raster)
library(metR)
library(zyp)
library(pals)
library(scales) #squish
library(data.table)
library(ggExtra)
library(hydroGOF)

plot_path =  'E:/phd/phd_20210618/Documents/8_plots/CAMELS_trends/'

# Prepare files -----------------------------------------------------------

#load for snowmelt function
source("E:/phd/phd_20210618/Documents/1_Data/Code_by_me/variable_method_flood_mech_functions.R")

#data files
#ERA5 Soil moisture
fpath_out = "C:/Users/ls16959/Data/01_Climate/ERA5-land/CAMELS"
load(file = paste0(fpath_out, "/ERA5_date.Rdata"))
load(file = paste0(fpath_out, "/ERA5_SM_cat_df.Rdata"))

#MSWEP precipitation
fpath = "C:/Users/ls16959/Data/01_Climate/MSWEP_V2.2_010deg_20180326/"
load(file = paste0(fpath, "/MSWEP2_2_camels_out_df_list.Rdata"))
MSWEPdf = do.call(rbind, lapply(out_df_list, function(df){t(as.data.frame(df)[,-1])}))
MSWEP_date = as.Date(sub("X","", rownames(MSWEPdf)), format = "%Y.%m.%d")

#actual evapotranspiration
fpath_GLEAM = "C:/Users/ls16959/Data/01_Climate/GLEAM_v3.3a/GLEAM_v33a_AET"
load(file = paste0(fpath_GLEAM, "/GLEAM_camels_df.Rdata"))
load(file = paste0(fpath_GLEAM, "/GLEAM_date.Rdata"))
#Temperature
fpath_BEST = "C:/Users/ls16959/Data/01_Climate/BerkleyEarth/original_T_avg_20200713"
load(file = paste0(fpath_BEST, "/BEST_camels_df.Rdata"))
load(file = paste0(fpath_BEST, "/BEST_date.Rdata"))

#Available Water Capacity
load(file = "E:/phd/phd_20210618/Data/CAMELS/AWC_catchment_list.rdata")


start_date = max(c(min(MSWEP_date), min(ERA5_date), min(GLEAM_date), min(BEST_date)))
end_date = min(c(max(MSWEP_date), max(ERA5_date),  max(GLEAM_date), max(BEST_date)))
date_seq = seq(start_date, end_date, by = "days")
date_seq_df = data.frame(Q_date = date_seq, Q_years = chron::years(date_seq))


#CAMELS files
ALL_Q_path = "E:/phd/phd_20210618/Data/CAMELS/CAMELS_ALL_Q/"

load(file = "E:/phd/phd_20210618/Data/CAMELS/camels_all_clim.Rdata")
#POT/Q files

#POT files calculated during previous analysis. Publication can be found here: https://doi.org/10.1029/2020WR028300
#Code repository can be found here: http://doi.org/10.5281/zenodo.4277642.


load(file = "E:/phd/phd_20210618/Data/CAMELS/cat_name_shape.Rdata")
load(file = "E:/phd/phd_20210618/Data/CAMELS/All_QF_nosnowrestr.Rdata")
load(file = "E:/phd/phd_20210618/Data/CAMELS/QC_POT_list.Rdata")

qfilesR = list.files(ALL_Q_path, full.names = T, pattern = ".Rdata")

POT_names = unlist(lapply(QC_POT_list, function(x){x[[1]]}))
POT_names_shape = sub("^US_0", "US_", POT_names)

isnullPOT = unlist(lapply(QC_POT_list, function(x){is.null(x[[5]][[1]])}))
#if quality is not good enough for POT, do not use for this analysis either
years_included_list = lapply(QC_POT_list, FUN = function(x){df = x[[4]]})
enough_years_list = lapply(QC_POT_list, FUN = function(x){df = x[[3]]})



# Calculate time series ---------------------------------------------------
#Combined for full timeseries of daily precipitation, soil moisture, ET, temperature, streamflow and the indices 7-day precipitation, snowmelt, 7-day snowmelt and snowstorage calculated for each CAMELS catchment with exception of the ones removed due to quality issues or lack of sufficient data. 


system.time({
  clim_df_list = lapply(c(1:length(POT_names)), function(catnum){
    print(catnum)
    if(enough_years_list[[catnum]]>=20){
      catid = POT_names[catnum]
      load(qfilesR[grep(paste0(catid,"_"), qfilesR)]) #need paste to avoid identifying multiple files
      
      catid = POT_names_shape[catnum]
      cat_shape_ind = which(cat_name_shape == catid)
      
      #Combine data
      MSWEP_tempdf = data.frame(Date = MSWEP_date, Prec = MSWEPdf[,cat_shape_ind], month_ind = substr(MSWEP_date, 1,7))
      ERA5_tempdf = data.frame(SM = SM_cat_df[,cat_shape_ind], month_ind = substr(ERA5_date, 1,7))
      Q_tempdf = data.frame(Date = as.Date(temp$Q_date, formate = "%Y-%m-%d"), Q = temp$Q_m3s)
      GLEAM_tempdf = data.frame(Date = GLEAM_date, ET = GLEAMdf[,cat_shape_ind])
      BEST_tempdf = data.frame(Date = BEST_date, Temp = BESTdf[,cat_shape_ind])
      
      #check for continous (zero) values 
      #identify repeat values with more than a year no change (long period to not exclude zero values in arid catchments)
      test = rle(Q_tempdf$Q)
      ident_test = rep(c(1:length(test$lengths)), times = test$lengths)
      Q_tempdf[which(ident_test == which(test$lengths>365)),"Q"] = NA
      
      
      #combine data
      climdf = merge(MSWEP_tempdf, ERA5_tempdf)
      climdf = list(climdf,GLEAM_tempdf, BEST_tempdf, Q_tempdf) %>% reduce(left_join, by = "Date")
      
      climdf[,"Prec7"] = rollsum(climdf$Prec, 7, fill = NA, align = "right")
      
      if(all(is.na(climdf$ET)) | all(is.na(climdf$Temp)) | all(is.na(climdf$Prec))){
        climdf = NULL
      }else{
        class_input_df_CAMELS = soil_snow_func(AWC_in = AWC_cat_vals[cat_shape_ind,1], ET_mat = climdf[,"ET"], rain_mat=climdf[,"Prec"], Temp_mat = climdf[,"Temp"], fdd_in = 2, Tcrit_in = 1, climdf$Date)
        Snow_vec = class_input_df_CAMELS[[1]][,"P_snow_vec"]
        
        climdf[,"snowmelt"] = class_input_df_CAMELS[[1]][,"P_snow_vec"]
        climdf[,"snow7"] = rollsum(climdf$snowmelt, 7, fill = NA, align = "right")
        
        climdf[,"snowstorage"] = class_input_df_CAMELS[[1]][,"Ssnow_vec"]
        
        climdf_years = chron::years(climdf$Date)
        years_incl = years_included_list[[catnum]]
        climdf[!(climdf_years %in% years_incl),] = rep(NA, 11)
        climdf[,"gauge_id"] = rep(POT_names_shape[catnum], nrow(climdf))
      }
    }else{
      climdf = NULL
    }
    return(climdf)
  })
})
#14min
save(clim_df_list, file = "E:/phd/phd_20210618/Data/CAMELS/clim_df_list.Rdata")
load(file = "E:/phd/phd_20210618/Data/CAMELS/clim_df_list.Rdata")


# Event values ----------------------------------------------------------
#Select Q events based on P, but keep all magnitudes
#Q events missed likely have other process than rainfall (for example snowmelt/glacier melt/ice jam...)

load(file = "E:/phd/phd_20210618/Data/CAMELS/clim_df_list.Rdata")
system.time({
  event_signatures_list = lapply(c(1:length(clim_df_list)), function(dfnum){
    print(dfnum)
    climdf = clim_df_list[[dfnum]]
    if(is.null(climdf)){
      event_signatures = NULL
    }else{
      NA_test = apply(climdf[,c("Prec", "SM", "Q")], 2, function(x){(sum(is.na(x))/length(x))>0.5})
      if(any(NA_test)){
        event_signatures = NULL
      }else{
        #calculate lag time
        lag_object = ccf(climdf$Prec, climdf$Q, lag.max = 100, na.action = na.pass, plot = F)
        lag_time = lag_object$lag[which.max(lag_object$acf)]
        
        ecdfP = ecdf(climdf$Prec)
        ecdfP7 = ecdf(climdf$Prec7)
        ecdfSM = ecdf(climdf$SM)
        ecdfQ = ecdf(climdf$Q)
        
        
        Ptest = climdf[,"Prec"]
        Ptest[Ptest<1] = 0
        test = rle(Ptest>=1)
        testdf = data.frame(lengthseq = test$lengths, valueseq = test$values)
        ident_test = rep(c(1:length(testdf[,1])), times = test$lengths)
        include_event = rep(testdf$valueseq, times = test$lengths)
        
        group = as.data.table(data.frame(EventID = ident_test, Prec = climdf$Prec, Ind = climdf$Date, SM = climdf$SM, snowmelt = climdf$snowmelt, include_event = include_event))
        group = group[group$include_event,] #only keep rainfall and not inbetween values
        SM_startval = group[group[, .I[which.min(Ind)], by=EventID]$V1][,c(1,4)] #only select ID, start SM 
        Pevent_max = group[group[, .I[which.max(Prec)], by=EventID]$V1][,c(1,2,3)]
        setnames(Pevent_max, c("EventID", "Pevent_max", "Pevent_maxdate"))
        Pvolume = group[, list(Pvolume = sum(Prec)), by=EventID]
        snowvolume = group[, list(snowvolume = sum(snowmelt)), by=EventID]
        Peventlength = group[, .N,by=list(EventID)]
        Pstartdate = group[group[, .I[which.min(Ind)], by=EventID]$V1][,c(1,3)]
        setnames(Pstartdate, c("EventID", "Pstartdate"))
        Penddate = group[group[, .I[which.max(Ind)], by=EventID]$V1][,c(1,3)]
        setnames(Penddate, c("EventID", "Penddate"))
        
        event_signatures = data.frame(list(Pevent_max, Pvolume, SM_startval, snowvolume, Pstartdate, Penddate, Peventlength) %>% reduce(left_join, by = "EventID"))
        
        
        Pmax_Qdf = do.call(rbind, lapply(c(1:nrow(event_signatures)),  function(x){
          event_range = climdf[climdf$Date %in% seq(from = event_signatures[x,"Pevent_maxdate"], to = event_signatures[x,"Penddate"]+5, by = "days"),]
          max_Q = max(event_range$Q)
          if(is.na(max_Q)){
            max_Q_date = NA
          }else{
            max_Q_date =event_range$Date[which.max(event_range$Q)]
          }
          return(data.frame(max_Q, max_Q_date))
        }))
        
        event_signatures = cbind(event_signatures, Pmax_Qdf)
        #remove events where Q is NA
        event_signatures = event_signatures[!is.na(event_signatures$max_Q),]
        #remove duplicated entries, where max Q is associated with multiple P events, only use the latest P event
        event_signatures = event_signatures[!duplicated(event_signatures$max_Q_date, fromLast = T),]
        #calculate quantiles
        
        
        #add P7 and snow7 based on Qmax date
        event_signatures[,"Prec7"] = climdf$Prec7[climdf$Date %in% event_signatures$max_Q_date]
        event_signatures[,"snow7"] = climdf$snow7[climdf$Date %in% event_signatures$max_Q_date]
        event_signatures[,"snowstorage"] = climdf$snowstorage[climdf$Date %in% event_signatures$max_Q_date]
        
        
        event_signatures[,"Pquant"] = ecdfP(event_signatures$Pevent_max)
        event_signatures[,"P7quant"] = ecdfP7(event_signatures$Prec7)
        event_signatures[,"Qquant"] = ecdfQ(event_signatures$max_Q)
        #event_signatures[,"POTquant"] = ecdfQ(event_signatures$Qobs)
        event_signatures[,"SMquant"] = ecdfSM(event_signatures$SM)
        event_signatures[,"lag_time"] = rep(lag_time, nrow(event_signatures))
        event_signatures[,"gauge_id"] = rep(climdf$gauge_id[1], nrow(event_signatures))
        
        #Test independence of events: if two events are closer or as close as the lag time, remove the smaller one
        i = 2
        while(i<nrow(event_signatures)){
          date_diff = event_signatures$max_Q_date[i]-event_signatures$max_Q_date[i-1]
          if(date_diff<=abs(lag_time)){
            if(event_signatures$max_Q[i]<event_signatures$max_Q[i-1]){
              event_signatures = event_signatures[-(i),]
              i = i
            }else{
              event_signatures = event_signatures[-(i-1),]
              i = i
            }
          }else{
            event_signatures = event_signatures
            i = i+1
          }
        }
        
        
      }
    }
    
    
    return(event_signatures)
  })
  
})
#52 min
save(event_signatures_list, file = "E:/phd/phd_20210618/Data/CAMELS/event_signatures_list.Rdata")
load(file = "E:/phd/phd_20210618/Data/CAMELS/event_signatures_list.Rdata")


Prec_dfall = do.call(rbind, lapply(event_signatures_list, function(df){
  if(!is.null(df)){
    df[,c("Pquant","P7quant", "Qquant", "SMquant","snow7", "snowstorage","Pevent_max","Pvolume", "Prec7","SM","max_Q", "lag_time", "gauge_id")]
  }
}))

#remove all events with snow influence
Prec_dfall = Prec_dfall[which(Prec_dfall$snow7<1),]
Prec_dfall = Prec_dfall[which(Prec_dfall$snowstorage<1),]
#remove all catchments with positive lag times (1 for uncertainty) and more than 5 days lag time
Prec_dfall = Prec_dfall[which(Prec_dfall$lag_time<=1&Prec_dfall$lag_time>=(-5)),]

Prec_dfall[,"country"] = substr(Prec_dfall$gauge_id,1,2)
#Combine with climate data
load(file = "E:/phd/phd_20210618/Data/CAMELS/camels_all_clim.Rdata")
Prec_dfallclim = merge(Prec_dfall, camels_all_clim[,c("gauge_id", "aridity")], by = "gauge_id", all.x = T)
Prec_dfallclim[,"arid_dir"] = cut(Prec_dfallclim$aridity, breaks=c(0,0.5, 1, 1.5,30), right = FALSE)

#total number of events
nrow(Prec_dfallclim)
#total number of catchments included
length(unique(Prec_dfallclim$gauge_id))

# Map Aridity -------------------------------------------------------------

isnotnullcatnames = data.frame(gauge_id =unique(Prec_dfall$gauge_id))

clim_aridity = camels_all_clim[,c("gauge_id", "aridity")]
clim_aridity = merge(isnotnullcatnames, clim_aridity)

clim_aridity[,"arid_dir"] = cut(clim_aridity$aridity, breaks=c(0,0.5, 1, 1.5,30), right = FALSE)
clim_aridity[,"country"] = substr(clim_aridity$gauge_id, 1,2)

load(file = "E:/phd/phd_20210618/Data/CAMELS/camels_all_loctopo.Rdata")
arid_mapdf = merge(clim_aridity, camels_all_loctopo, by = "gauge_id")
worldmap <- map_data("world")

p_map = ggplot() +  geom_map(data = worldmap, map = worldmap, mapping = aes(map_id = region), fill = "white", col = "gray90") +
  geom_point(data = arid_mapdf, aes(x = gauge_lon, y = gauge_lat, group= arid_dir, colour = arid_dir))+scale_colour_brewer(name = "Aridity", palette = "BrBG", direction = -1, labels = c("Very wet", "Wet", "Dry", "Very dry"))+theme_minimal(base_size = 13)+xlab("Longitude") +
  ylab("Latitude")+theme(legend.position = "none", text =  element_text(size = 16))+coord_fixed(ratio = 1.3)+scale_x_continuous(breaks = c(-100, 0, 100), labels = c(expression(paste(100^o,~'W')), expression(paste(0^o)), expression(paste(100^o,~'E'))))+scale_y_continuous(breaks = c(-30, 0, 30, 60), labels = c(expression(paste(30^o,~'S')), expression(paste(0^o)), expression(paste(30^o,~'N')), expression(paste(60^o,~'N'))))

p_bar = ggplot(arid_mapdf, aes(arid_dir, fill = arid_dir))+geom_bar(stat = "count")+scale_fill_brewer(palette = "BrBG", direction = -1,)+scale_x_discrete(labels=c("Very wet", "Wet", "Dry", "Very dry"))+theme_minimal(base_size = 13)+theme(legend.position = "none", panel.grid.major.x = element_blank(),axis.text=element_text(size = 14), axis.title = element_text(size = 16),plot.background = element_rect(fill = "white", colour = "gray90"))+xlab("")+ylab("No. Catchments")

p_map + inset_element(p_bar, left = 0.6, bottom = 0.6, right = 1, top = 1)

ggsave(paste0(plot_path, "/CAMELS_ariditymap.pdf"), width = 12, height = 7 , dpi = "retina")





# Raster focal method ----------------------------------------------------

#1-day precipitation 

testpts = Prec_dfallclim
# set up an 'empty' raster, here via an extent object derived from your data
r <- raster(xmn=0, xmx=1, ymx=1, ymn=0.5, ncol=100, nrow=50) #changing limit to 1.1 does not change the outcome
r <- rasterize(testpts[, c("SMquant", "Pquant")], r, testpts[,"Qquant"], fun=median)
r2 <- focal(r, w= matrix(1,5,5), mean)
r2dfPmax = data.frame(rasterToPoints(r2))
r2dfPmax[,"arid_dir"] = rep("All", nrow(r2dfPmax))


contour_dflist_Pmax = lapply(unique(Prec_dfallclim$arid_dir), function(arid_ind){
  testpts = Prec_dfallclim[Prec_dfallclim$arid_dir == arid_ind,]
  
  # set up an 'empty' raster, here via an extent object derived from your data
  r <- raster(xmn=0, xmx=1, ymx=1, ymn=0.5, ncol=100, nrow=50) #changing limit to 1.1 does not change the outcome
  r <- rasterize(testpts[, c("SMquant", "Pquant")], r, testpts[,"Qquant"], fun=median)
  r2 <- focal(r, w= matrix(1,5,5), mean)
  r2df = data.frame(rasterToPoints(r2))
  r2df[,"arid_dir"] = rep(arid_ind, nrow(r2df))
  
  return(r2df)
})
contour_df = do.call(rbind, contour_dflist_Pmax)
contour_df = rbind(r2dfPmax, contour_df)
contour_df[,"arid_dir"] = factor(contour_df$arid_dir, levels = c("All", "[0,0.5)", "[0.5,1)","[1,1.5)", "[1.5,30)"))
var.labs <- c("(a) All", "(b) Very wet", "(c) Wet", "(d) Dry", "(e) Very dry")
names(var.labs) <- c("All", "[0,0.5)", "[0.5,1)","[1,1.5)", "[1.5,30)")



#7-day precipitation ----

testpts = Prec_dfallclim
# set up an 'empty' raster, here via an extent object derived from your data
r <- raster(xmn=0, xmx=1, ymx=1, ymn=0.5, ncol=100, nrow=50) #changing limit to 1.1 does not change the outcome
r <- rasterize(testpts[, c("SMquant", "P7quant")], r, testpts[,"Qquant"], fun=median)
r2 <- focal(r, w= matrix(1,5,5), mean)
r2dfP7 = data.frame(rasterToPoints(r2))
r2dfP7[,"arid_dir"] = rep("All", nrow(r2dfP7))



contour_dflist_P7 = lapply(unique(Prec_dfallclim$arid_dir), function(arid_ind){
  testpts = Prec_dfallclim[Prec_dfallclim$arid_dir == arid_ind,]
  
  # set up an 'empty' raster, here via an extent object derived from your data
  r <- raster(xmn=0, xmx=1, ymx=1, ymn=0.5, ncol=100, nrow=50) #changing limit to 1.1 does not change the outcome
  r <- rasterize(testpts[, c("SMquant", "P7quant")], r, testpts[,"Qquant"], fun=median)
  r2 <- focal(r, w= matrix(1,5,5), mean)
  r2df = data.frame(rasterToPoints(r2))
  r2df[,"arid_dir"] = rep(arid_ind, nrow(r2df))
  
  return(r2df)
})
contour_dfP7 = do.call(rbind, contour_dflist_P7)
contour_dfP7 = rbind(r2dfP7, contour_dfP7)
contour_dfP7[,"arid_dir"] = factor(contour_dfP7$arid_dir, levels = c("All", "[0,0.5)", "[0.5,1)","[1,1.5)", "[1.5,30)"))




contour_df100 = contour_df
contour_df100[,1:3] = contour_df100[,1:3]*100
p1 = ggplot(contour_df100[contour_df100$arid_dir == "All",], aes(x, y, fill = layer)) +
  geom_raster() +
  geom_contour(aes(z = layer), col = "white", binwidth = 5)+
  geom_contour(aes(z = layer), col = "black", breaks = 90, lty = 2)+
  scale_fill_viridis_c(expression(paste("Med(Pr(Q<Qmax))"," [%]")), limits=c(0,100), direction = -1)+geom_text_contour(aes(z = layer), col = "white")+theme_minimal(base_size = 13)+xlab(expression(paste("Pr(SM<SMI)"," [%]")))+ylab(expression(paste("Pr(P<Pmax)"," [%]")))+facet_wrap(~arid_dir, labeller = labeller(arid_dir = var.labs), nrow = 2 )+theme(strip.text = element_text(size = 13), legend.position = "bottom", legend.margin=margin(-10, 0, 0, 0))

p2 = ggplot(contour_df100[contour_df100$arid_dir != "All",], aes(x, y, fill = layer)) +
  geom_raster() +
  geom_contour(aes(z = layer), col = "white", binwidth = 5)+
  geom_contour(aes(z = layer), col = "black", breaks = 90, lty = 2)+
  scale_fill_viridis_c(expression(paste("Med(Pr(Q<Qmax))"," [%]")), limits=c(0,100), direction = -1)+geom_text_contour(aes(z = layer), col = "white")+theme_minimal(base_size = 13)+xlab(expression(paste("Pr(SM<SMI)"," [%]")))+ylab(expression(paste("Pr(P<Pmax)"," [%]")))+facet_wrap(~arid_dir, labeller = labeller(arid_dir = var.labs), nrow = 2 )+theme(strip.text = element_text(size = 13), legend.position = "bottom", legend.margin=margin(-10, 0, 0, 0))

p1+p2+plot_layout(guides = "collect")& theme(legend.position = 'bottom')
ggsave(paste0(plot_path, "/CAMELS_SMPQquant_contour_aridity.pdf"), width = 10, height = 6 , dpi = "retina")

## mean/sd
fun_in = sd

testpts = Prec_dfallclim
# set up an 'empty' raster, here via an extent object derived from your data
r <- raster(xmn=0, xmx=1, ymx=1, ymn=0.5, ncol=100, nrow=50) #changing limit to 1.1 does not change the outcome
r <- rasterize(testpts[, c("SMquant", "Pquant")], r, testpts[,"Qquant"], fun=fun_in)
r2 <- focal(r, w= matrix(1,5,5), mean)
r2dfPmax = data.frame(rasterToPoints(r2))
r2dfPmax[,"arid_dir"] = rep("All", nrow(r2dfPmax))


contour_dflist_Pmax = lapply(unique(Prec_dfallclim$arid_dir), function(arid_ind){
  testpts = Prec_dfallclim[Prec_dfallclim$arid_dir == arid_ind,]
  
  # set up an 'empty' raster, here via an extent object derived from your data
  r <- raster(xmn=0, xmx=1, ymx=1, ymn=0.5, ncol=100, nrow=50) #changing limit to 1.1 does not change the outcome
  r <- rasterize(testpts[, c("SMquant", "Pquant")], r, testpts[,"Qquant"], fun=fun_in)
  r2 <- focal(r, w= matrix(1,5,5), mean)
  r2df = data.frame(rasterToPoints(r2))
  r2df[,"arid_dir"] = rep(arid_ind, nrow(r2df))
  
  return(r2df)
})
contour_df = do.call(rbind, contour_dflist_Pmax)
contour_df = rbind(r2dfPmax, contour_df)
contour_df[,"arid_dir"] = factor(contour_df$arid_dir, levels = c("All", "[0,0.5)", "[0.5,1)","[1,1.5)", "[1.5,30)"))


contour_df100 = contour_df
contour_df100[,1:3] = contour_df100[,1:3]*100
ggplot(contour_df100, aes(x, y, fill = layer)) +  
  geom_tile() +
  geom_contour(aes(z = layer), col = "white", binwidth = 5)+
  geom_contour(aes(z = layer), col = "black", breaks = 90, lty = 2)+
  scale_fill_viridis_c(expression(paste(bar("Pr(Q>Qmax)")," [%]")), limits=c(0,100), direction = -1)+geom_text_contour(aes(z = layer), col = "white")+theme_minimal(base_size = 13)+xlab(expression(paste("Pr(SM<SMI)"," [%]")))+ylab(expression(paste("Pr(P<Pmax)"," [%]")))+facet_wrap(~arid_dir, labeller = labeller(arid_dir = var.labs), nrow = 1 )+theme(strip.text = element_text(size = 13), legend.position = "bottom", legend.margin=margin(-10, 0, 0, 0))
ggsave(paste0(plot_path, "/CAMELS_SMPQquant_contour_aridity_mean.pdf"), width = 14, height = 4 , dpi = "retina")


ggplot(contour_df100, aes(x, y, fill = layer)) +  
  geom_tile() +
  scale_fill_viridis_c(option = "A", expression(paste(sigma, " Pr(Q>Qmax)"," [%]")), limits=c(0,30), direction = -1)+theme_minimal(base_size = 13)+xlab(expression(paste("Pr(SM<SMI)"," [%]")))+ylab(expression(paste("Pr(P<Pmax)"," [%]")))+facet_wrap(~arid_dir, labeller = labeller(arid_dir = var.labs), nrow = 1 )+theme(strip.text = element_text(size = 13), legend.position = "bottom", legend.margin=margin(-10, 0, 0, 0))
ggsave(paste0(plot_path, "/CAMELS_SMPQquant_contour_aridity_sd.pdf"), width = 14, height = 4 , dpi = "retina")



#mean vs median
test = Prec_dfallclim[Prec_dfallclim$arid_dir == "[1.5,30)",]
test = test[test$SMquant>0.99,]
dim(test)
hist(test$Qquant)
plot(test$Qquant, test$Pquant)
test[,"P_cat"] = cut(test$Pquant, breaks=seq(0.7, 1,0.05), right = FALSE)
ggplot(test, aes(P_cat, Qquant))+geom_boxplot()+stat_summary(fun=mean, geom="point", shape=20, size = 7, color="red", fill="red")

##




contour_dfP7_100 = contour_dfP7
contour_dfP7_100[,1:3] = contour_dfP7_100[,1:3]*100
ggplot(contour_dfP7_100, aes(x, y, fill = layer)) +  
  geom_tile() +
  geom_contour(aes(z = layer), col = "white", binwidth = 5)+
  geom_contour(aes(z = layer), col = "black", breaks = 90, lty = 2)+
  scale_fill_viridis_c(expression(paste("Med(Pr(Q<Qmax))"," [%]")), limits=c(0,100), direction = -1)+geom_text_contour(aes(z = layer), col = "white")+theme_minimal(base_size = 13)+xlab(expression(paste("Pr(SM<SMI)"," [%]")))+ylab(expression(paste("Pr(P7<P7i)"," [%]")))+facet_wrap(~arid_dir, labeller = labeller(arid_dir = var.labs), nrow = 1 )+theme(strip.text = element_text(size = 13), legend.position = "bottom", legend.margin=margin(-10, 0, 0, 0))
ggsave(paste0(plot_path, "/CAMELS_SMP7Qquant_contour_aridity.pdf"), width = 14, height = 4 , dpi = "retina")

#supplement: unsmoothed raster without contours:


testpts = Prec_dfallclim
# set up an 'empty' raster, here via an extent object derived from your data
r <- raster(xmn=0, xmx=1, ymx=1, ymn=0.5, ncol=100, nrow=50) #changing limit to 1.1 does not change the outcome
r2 <- rasterize(testpts[, c("SMquant", "Pquant")], r, testpts[,"Qquant"], fun=median)
r2dfPmax = data.frame(rasterToPoints(r2))
r2dfPmax[,"arid_dir"] = rep("All", nrow(r2dfPmax))


contour_dflist_Pmax = lapply(unique(Prec_dfallclim$arid_dir), function(arid_ind){
  testpts = Prec_dfallclim[Prec_dfallclim$arid_dir == arid_ind,]
  
  # set up an 'empty' raster, here via an extent object derived from your data
  r <- raster(xmn=0, xmx=1, ymx=1, ymn=0.5, ncol=100, nrow=50) #changing limit to 1.1 does not change the outcome
  r2 <- rasterize(testpts[, c("SMquant", "Pquant")], r, testpts[,"Qquant"], fun=median)
  r2df = data.frame(rasterToPoints(r2))
  r2df[,"arid_dir"] = rep(arid_ind, nrow(r2df))
  
  return(r2df)
})
contour_df = do.call(rbind, contour_dflist_Pmax)
contour_df = rbind(r2dfPmax, contour_df)
contour_df[,"arid_dir"] = factor(contour_df$arid_dir, levels = c("All", "[0,0.5)", "[0.5,1)","[1,1.5)", "[1.5,30)"))


contour_df100 = contour_df
contour_df100[,1:3] = contour_df100[,1:3]*100
ggplot(contour_df100, aes(x, y, fill = layer)) +  
  geom_tile() +
  scale_fill_viridis_c(expression(paste("Med(Pr(Q<Qmax))"," [%]")), limits=c(0,100), direction = -1)+theme_minimal(base_size = 13)+xlab(expression(paste("Pr(SM<SMI)"," [%]")))+ylab(expression(paste("Pr(P<Pmax)"," [%]")))+facet_wrap(~arid_dir, labeller = labeller(arid_dir = var.labs), nrow = 1 )+theme(strip.text = element_text(size = 13), legend.position = "bottom", legend.margin=margin(-10, 0, 0, 0))
ggsave(paste0(plot_path, "/CAMELS_SMPQquant_raster_aridity.pdf"), width = 14, height = 4 , dpi = "retina")





# High flow plots ---------------------------------------------------------


#Percent threshold
quant_thresh = 0.99

Precdf_highQ = Prec_dfallclim[Prec_dfallclim$Qquant>quant_thresh,]
Precdf_highQall = Precdf_highQ
Precdf_highQall[,"arid_dir"] = rep("All", nrow(Precdf_highQall))
Precdf_highQ = rbind(Precdf_highQall, Precdf_highQ)
Precdf_highQ[,"arid_dir"] = factor(Precdf_highQ$arid_dir, levels = c("All", "[0,0.5)", "[0.5,1)","[1,1.5)", "[1.5,30)"))

Precdf_highQ[,c("Pquant", "SMquant", "Qquant", "P7quant")] = (Precdf_highQ[,c("Pquant", "SMquant", "Qquant", "P7quant")]*100)

p1 = ggplot(Precdf_highQ[Precdf_highQ$arid_dir == "All",], aes(SMquant, Pquant, col = Qquant))+geom_point(pch = ".", size = 0.3)+scale_colour_viridis(expression(paste("Med(Pr(Q<Qmax))"," [%]")), limits = c(0, 100), direction = -1)+geom_hline(yintercept = 99, colour = "red", lty = 2)+geom_hline(yintercept = 95, colour = "red")+xlab(expression(paste("Pr(SM<SMI)"," [%]")))+ylab(expression(paste("Pr(P<Pmax)"," [%]")))+facet_wrap(~arid_dir, labeller = labeller(arid_dir = var.labs), nrow = 1)+theme_minimal(base_size = 13)+theme(legend.position = "none", strip.text = element_text(size = 13))+scale_y_continuous(limits = c(50, 100))

p2 = ggplot(Precdf_highQ[Precdf_highQ$arid_dir != "All",], aes(SMquant, Pquant, col = Qquant))+geom_point(pch = ".", size = 0.3)+scale_colour_viridis(expression(paste("Med(Pr(Q<Qmax))"," [%]")), limits = c(0, 100), direction = -1)+geom_hline(yintercept = 99, colour = "red", lty = 2)+geom_hline(yintercept = 95, colour = "red")+xlab(expression(paste("Pr(SM<SMI)"," [%]")))+ylab(expression(paste("Pr(P<Pmax)"," [%]")))+facet_wrap(~arid_dir, labeller = labeller(arid_dir = var.labs), nrow = 2)+theme_minimal(base_size = 13)+theme(legend.position = "none", strip.text = element_text(size = 13))+scale_y_continuous(limits = c(50, 100))

p1 + p2 + plot_layout()

#+scale_x_continuous(labels = function(x)(x*100))+scale_y_continuous(labels = function(x)(x*100), limits = c(0.5,1))
ggsave(paste0(plot_path, "/CAMELS_SMPQ99quant_contour_aridity.pdf"), width = 10, height = 6 , dpi = "retina")



ggplot(Precdf_highQ, aes(SMquant, P7quant, col = Qquant))+geom_point(pch = ".", size = 0.3)+scale_colour_viridis(expression(paste("Med(Pr(Q<Qmax))"," [%]")), limits = c(0, 100), direction = -1)+geom_hline(yintercept = 99, colour = "red", lty = 2)+geom_hline(yintercept = 95, colour = "red")+xlab(expression(paste("Pr(SM<SMI)"," [%]")))+ylab(expression(paste("Pr(P7<P7i)"," [%]")))+facet_wrap(~arid_dir, labeller = labeller(arid_dir = var.labs), nrow = 1)+theme_minimal(base_size = 13)+theme(legend.position = "none", strip.text = element_text(size = 11))+   scale_y_continuous(limits = c(50, 100))

ggsave(paste0(plot_path, "/CAMELS_SMP7Q99quant_contour_aridity.pdf"), width = 12, height = 3 , dpi = "retina")

#Compare event distributions for lower values "for instance, many more high flow events are reported for precipitation percentile $<$ XX and soil moisture percentile $<$ YY in panel B and C than in D and E"

lapply(c("All", "[0,0.5)", "[0.5,1)","[1,1.5)", "[1.5,30)"), function(clim){
  Precdf_highQtemp = Precdf_highQ[Precdf_highQ$arid_dir == clim,]
  nrow(Precdf_highQtemp[(Precdf_highQtemp$Pquant>1&Precdf_highQtemp$SMquant>25),])/nrow(Precdf_highQtemp)
})


#Quantile averages boxplots----

quant_thresh = 0.99
tempPselect = Prec_dfallclim[Prec_dfallclim$Pquant>quant_thresh,]
tempQselect = Prec_dfallclim[Prec_dfallclim$Qquant>quant_thresh,]

aggregatetodf = function(agg_out){cbind(agg_out[-ncol(agg_out)], agg_out[[ncol(agg_out)]])}

#aggregate median value per catchment
tempP = aggregatetodf(aggregate(tempPselect[,c("Pquant", "Qquant", "SMquant")], by = list(tempPselect$gauge_id), FUN=median))
tempP[,"selectby"] = rep("(a) P>P99", nrow(tempP))
tempQ = aggregatetodf(aggregate(tempQselect[,c("Pquant", "Qquant", "SMquant")], by = list(tempQselect$gauge_id), FUN=median))
tempQ[,"selectby"] = rep("(b) Q>Q99", nrow(tempQ))
tempdf = rbind(tempP, tempQ)
colnames(tempdf) = c("gauge_id", "P", "Q", "SM", "selectby")

tempdf = merge(tempdf, camels_all_clim[,c("gauge_id", "aridity")], by = "gauge_id", all.x = T)
tempdf[,"arid_dir"] = cut(tempdf$aridity, breaks=c(0,0.5, 1, 1.5,30), right = FALSE)

tempdf_melt = reshape2::melt(tempdf[,c("P", "Q", "SM","selectby", "arid_dir")], id = c("arid_dir", "selectby"))
tempdf_melt[,"value"] = (tempdf_melt[,"value"])*100


ggplot(tempdf_melt, aes(x = variable, value, col = arid_dir))+geom_boxplot(lwd = 1.2, fatten = 0.9)+facet_wrap(~selectby)+theme_minimal(base_size = 13)+scale_colour_brewer(name = "Aridity", palette = "BrBG", direction = -1, labels = c("Very wet", "Wet", "Dry", "Very dry"))+xlab("")+ylab("Median non-exceedance probability [%]")+ylim(0,100)+theme(strip.text.x = element_text(size = 14), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(), panel.spacing = unit(3, "lines"))
ggsave(paste0(plot_path, "/CAMELS_quantile_median_aridity.pdf"), width = 9, height = 4 , dpi = "retina")


# all events, not aggregated by catchments

tempPselect[,"selectby"] = rep("(a) P>P99", nrow(tempPselect))
tempQselect[,"selectby"] = rep("(b) Q>Q99", nrow(tempQselect))
tempdf = rbind(tempPselect[,c("Pquant", "Qquant", "SMquant", "selectby", "arid_dir")], tempQselect[,c("Pquant", "Qquant", "SMquant", "selectby", "arid_dir")])
colnames(tempdf) = c("P", "Q", "SM", "selectby", "arid_dir")
tempdf_melt = reshape2::melt(tempdf[,c("P", "Q", "SM","selectby", "arid_dir")], id = c("arid_dir", "selectby"))
tempdf_melt[,"value"] = (tempdf_melt[,"value"])*100

ggplot(tempdf_melt, aes(x = variable, value, col = arid_dir))+geom_boxplot(lwd = 1.2, fatten = 0.9)+facet_wrap(~selectby)+theme_minimal(base_size = 13)+scale_colour_brewer(name = "Aridity", palette = "BrBG", direction = -1, labels = c("Very wet", "Wet", "Dry", "Very dry"))+xlab("")+ylab("non-exceedance probability [%]")+theme(strip.text.x = element_text(size = 14), panel.grid = element_line())

ggsave(paste0(plot_path, "/CAMELS_quantile_boxplot_all_aridity.pdf"), width = 9, height = 4 , dpi = "retina")

# Quantile averages boxplots P7 ----------------------------------------------------

quant_thresh = 0.99
tempP7select = Prec_dfallclim[Prec_dfallclim$P7quant>quant_thresh,]
tempQselect = Prec_dfallclim[Prec_dfallclim$Qquant>quant_thresh,]

aggregatetodf = function(agg_out){cbind(agg_out[-ncol(agg_out)], agg_out[[ncol(agg_out)]])}


tempP7 = aggregatetodf(aggregate(tempP7select[,c("P7quant", "Qquant", "SMquant")], by = list(tempP7select$gauge_id), FUN=median))
tempP7[,"selectby"] = rep("(a) P7>P7 99", nrow(tempP7))
tempQ = aggregatetodf(aggregate(tempQselect[,c("P7quant", "Qquant", "SMquant")], by = list(tempQselect$gauge_id), FUN=median))
tempQ[,"selectby"] = rep("(b) Q>Q99", nrow(tempQ))
tempdf = rbind(tempP7, tempQ)
colnames(tempdf) = c("gauge_id", "P7", "Q", "SM", "selectby")

tempdf = merge(tempdf, camels_all_clim[,c("gauge_id", "aridity")], by = "gauge_id", all.x = T)
tempdf[,"arid_dir"] = cut(tempdf$aridity, breaks=c(0,0.5, 1, 1.5,30), right = FALSE)
tempdf_melt = reshape2::melt(tempdf[,c("P7", "Q", "SM","selectby", "arid_dir")], id = c("arid_dir", "selectby"))
tempdf_melt[,"value"] = (tempdf_melt[,"value"])*100

ggplot(tempdf_melt, aes(x = variable, value, col = arid_dir))+geom_boxplot(lwd = 1.2, fatten = 0.9)+facet_wrap(~selectby)+theme_minimal(base_size = 13)+scale_colour_brewer(name = "Aridity", palette = "BrBG", direction = -1, labels = c("Very wet", "Wet", "Dry", "Very dry"))+xlab("")+ylab("Median non-exceedance probability [%]")+ylim(0,100)+theme(strip.text.x = element_text(size = 14),  panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(), panel.spacing = unit(3, "lines"))
ggsave(paste0(plot_path, "/CAMELS_quantile_MedianP7_aridity.pdf"), width = 9, height = 4 , dpi = "retina")





# Compare sampling method for trends --------------------------------------



load(file = "E:/phd/phd_20210618/Data/CAMELS/camels_all_clim.Rdata")

#Percent threshold
quant_thresh = 0.99

sen_mk_func_slope = function(date_vals, trend_vals){
  slope = coef(zyp.sen(trend_vals~date_vals))[2]
  sensslope_10 = ((slope*10*365)/mean(trend_vals))*100
  mktest = MannKendall(trend_vals)
  pval = c(mktest$sl)
  return(c(sensslope_10, pval)) #sensslope_10
}


#Function to process and label output from lapply trend analysis
df_processing_func = function(list_in){
  df = do.call(rbind, lapply(list_in, function(df){
    if(is.null(df)){
      out = df
    }else{
      out = df
      colnames(out) = c("max_Q_slope", "max_Q_pval", "Pevent_max_slope","Pevent_max_pval","Pvolume_slope","Pvolume_pval",  "Prec7_slope","Prec7_pval", "SM_slope", "SM_pval", "Event_length_slope", "Event_length_pval", "gauge_id", "lag_time", "meanSM")
    }
    return(out)
  }))
  
  df = merge(df, camels_all_clim[,c("gauge_id", "aridity")])
  df[,"arid_dir"] = cut(df$aridity, breaks=c(0,0.5, 1, 1.5,30), right = FALSE)
  
  trend_dir = ifelse(df$max_Q_slope<0, "neg", "pos")
  trend_dir[df$max_Q_pval>0.05] = "non_sign"
  df[,"trend_dir"] = trend_dir
  df = df[abs(df$lag_time)<5,]
  
  return(df)
}


system.time({
  trend_quantiles_listQ = lapply(c(quant_thresh), function(Qquant_thresh){
    print(Qquant_thresh)
    trend_list = lapply(c(1:length(event_signatures_list)), function(catnum){
      tempdf = event_signatures_list[[catnum]]
      if(!is.null(tempdf)){
        if(tempdf$lag_time[1]<=1&tempdf$lag_time[1]>=(-5)){
          snow_ind = tempdf$snow7<1
          snow_ind[tempdf$snowstorage>1] = F
          if(!all(is.na(snow_ind))){
            #if there are snow values available include them, otherwise not
            tempdf = tempdf[snow_ind,]
          }
          #remove Na values for calculation
          tempdf = tempdf[!is.na(tempdf$max_Q_date),]
          #select quantile
          tempdf = tempdf[tempdf$Qquant>Qquant_thresh,]
          
          if(nrow(tempdf)>=20){
            #Theil Sen 
            theil_all = do.call(c,lapply(c("max_Q", "Pevent_max","Pvolume", "Prec7", "SM", "N"), function(vals){
              sen_mk_func_slope(as.numeric(tempdf$max_Q_date), tempdf[,vals])
            }))
            out = data.frame(t(theil_all))
          }else{
            out = data.frame(t(rep(NA, 12)))
          }
          
          
          
          out[,"gauge_id"] =  tempdf[!is.na(tempdf$gauge_id),"gauge_id"][1]
          out[,"lag_time"] =  tempdf[!is.na(tempdf$gauge_id),"lag_time"][1]
          out[,"meanSM_Q"] =  mean(tempdf[!is.na(tempdf$gauge_id),"SMquant"], na.rm = T)
        }else{
          out = NULL
        }
      }else{
        out = NULL
      }
      return(out)
    })
  })
})

trend_quantilesdfQ = df_processing_func(trend_quantiles_listQ[[1]])




#Pmax----

system.time({
  trend_quantiles_listPmax = lapply(c(quant_thresh), function(Pquant_thresh){
    print(Pquant_thresh)
    trend_list = lapply(c(1:length(event_signatures_list)), function(catnum){
      tempdf = event_signatures_list[[catnum]]
      if(!is.null(tempdf)){
        if(tempdf$lag_time[1]<=1&tempdf$lag_time[1]>=(-5)){
          snow_ind = tempdf$snow7<1
          snow_ind[tempdf$snowstorage>1] = F
          if(!all(is.na(snow_ind))){
            #if there are snow values available include them, otherwise not
            tempdf = tempdf[snow_ind,]
          }
          #remove Na values for calculation
          tempdf = tempdf[!is.na(tempdf$max_Q_date),]
          #select quantile
          tempdf = tempdf[tempdf$Pquant>Pquant_thresh,]
          
          if(nrow(tempdf)>=20){
            #Theil Sen 
            theil_all = do.call(c,lapply(c("max_Q", "Pevent_max","Pvolume", "Prec7", "SM", "N"), function(vals){
              sen_mk_func_slope(as.numeric(tempdf$max_Q_date), tempdf[,vals])
            }))
            out = data.frame(t(theil_all))
          }else{
            out = data.frame(t(rep(NA, 12)))
          }
          
          
          
          out[,"gauge_id"] =  tempdf[!is.na(tempdf$gauge_id),"gauge_id"][1]
          out[,"lag_time"] =  tempdf[!is.na(tempdf$gauge_id),"lag_time"][1]
          out[,"meanSM_Q"] =  mean(tempdf[!is.na(tempdf$gauge_id),"SMquant"], na.rm = T)
        }else{
          out = NULL
        }
        
      }else{
        out = NULL
      }
      return(out)
    })
  })
})


trend_quantilesdfP = df_processing_func(trend_quantiles_listPmax[[1]])



dfQP = data.frame(gauge_id = trend_quantilesdfQ$gauge_id, PmaxQslope = trend_quantilesdfQ$Pevent_max_slope, QPQslope = trend_quantilesdfQ$max_Q_slope, arid_dir = trend_quantilesdfQ$arid_dir)
dfPP = data.frame(gauge_id = trend_quantilesdfP$gauge_id, PmaxPslope = trend_quantilesdfP$Pevent_max_slope, PPQslope = trend_quantilesdfP$max_Q_slope)
compare_trendP = merge(dfQP, dfPP)
compare_trendPall = compare_trendP
compare_trendPall[,"arid_dir"] = rep("All", nrow(compare_trendPall))
compare_trendP = rbind(compare_trendP, compare_trendPall)


#Compare Q trends for all catchments
pall = ggplot(compare_trendP[compare_trendP$arid_dir == "All",], aes(QPQslope, PPQslope))+geom_point(alpha = 0.2)+xlim(-60,60)+geom_abline(slope = 1, intercept = 0, col = "red")+ylim(-60,60)+xlab(paste0("Streamflow trend select by Q>Q(", 100*quant_thresh,")"))+ylab(paste0("Streamflow trend select by P>P(",100*quant_thresh,")"))+theme_minimal(base_size = 13)+theme(legend.position = "none", text = element_text(size = 11), strip.text = element_text(size = 11))
ggMarginal(pall, type="histogram", binwidth = 2, colour = "gray50", fill = "gray90")


mean(compare_trendP[compare_trendP$arid_dir == "All","QPQslope"], na.rm = T)
mean(compare_trendP[compare_trendP$arid_dir == "All","PPQslope"], na.rm = T)


pdf(file = paste0(plot_path, "/CAMELS_trendcomparePmax_AllpointsQ99.pdf"), width = 5, height = 4 )
ggMarginal(pall, type="histogram", binwidth = 2, colour = "gray50", fill = "gray90")
dev.off()


#Evaluate comparison between trends
compare_trendP_allclim = compare_trendP[compare_trendP$arid_dir == "All",]
QPQslope_vals = compare_trendP_allclim$QPQslope
PPQslope_vals = compare_trendP_allclim$PPQslope
#RMSE for all points
hydroGOF::rmse(QPQslope_vals, PPQslope_vals, na.rm = T)
#RMSE with outliers removed
compare_trendP_allclim = tibble(compare_trendP_allclim)
outlier_rm_df = compare_trendP_allclim %>%
  filter(QPQslope>(-30),
         QPQslope<(30)) %>%
  filter(PPQslope>(-30),
         PPQslope<(30))

hydroGOF::rmse(outlier_rm_df$QPQslope, outlier_rm_df$PPQslope, na.rm = T)


#Plot trend outlier
plot(QPQslope_vals, PPQslope_vals)
compare_trendP_allclim[which(QPQslope_vals > 60),]
#GB_39088

event_sign_nonull = event_signatures_list[!unlist(lapply(event_signatures_list, is.null))]
event_sign_gaugeid = unlist(lapply(event_sign_nonull, function(x)(x$gauge_id[1])))
test = event_sign_nonull[[which(event_sign_gaugeid == "GB_39088")]]
plot(test$max_Q_date, test$max_Q)
points(test$max_Q_date[test$Qquant>quant_thresh], test$max_Q[test$Qquant>quant_thresh], col = "red")





#Compare P and Q trends for catchments split by aridity

var.labs <- c("(a) All", "(b) Very wet", "(c) Wet", "(d) Dry", "(e) Very dry")
names(var.labs) <- c("All", "[0,0.5)", "[0.5,1)","[1,1.5)", "[1.5,30)")

#streamflow
compare_trendP[,"arid_dir"] = factor(compare_trendP$arid_dir, levels = c("All", "[0,0.5)", "[0.5,1)","[1,1.5)", "[1.5,30)"))
ggplot(compare_trendP, aes(QPQslope, PPQslope))+geom_abline(slope = 1, intercept = 0, col = "gray90")+geom_point(alpha = 0.2)+xlim(-30,30)+ylim(-30,30)+xlab(paste0("Streamflow trend select by Q>Q(", quant_thresh,")"))+ylab(paste0("Streamflow trend select by P>P(",quant_thresh,")"))+facet_wrap(~arid_dir, labeller = labeller(arid_dir = var.labs), ncol = 5 )+theme_minimal(base_size = 13)+theme(legend.position = "none", text = element_text(size = 11), strip.text = element_text(size = 11))
ggsave(paste0(plot_path, "/CAMELS_trendcomparePmax_ariditypoints.pdf"), width = 12, height = 3 , dpi = "retina")




# Compare trend P7 --------------------------------------------------------

system.time({
  trend_quantiles_listP7 = lapply(c(quant_thresh), function(Pquant_thresh){
    print(Pquant_thresh)
    trend_list = lapply(c(1:length(event_signatures_list)), function(catnum){
      tempdf = event_signatures_list[[catnum]]
      if(!is.null(tempdf)){
        snow_ind = tempdf$snow7<1
        snow_ind[tempdf$snowstorage>1] = F
        if(!all(is.na(snow_ind))){
          #if there are snow values available include them, otherwise not
          tempdf = tempdf[snow_ind,]
        }
        #remove Na values for calculation
        tempdf = tempdf[!is.na(tempdf$max_Q_date),]
        #select quantile
        tempdf = tempdf[tempdf$P7quant>Pquant_thresh,]
        
        if(nrow(tempdf)>=20){
          #Theil Sen 
          theil_all = do.call(c,lapply(c("max_Q", "Pevent_max","Pvolume", "Prec7", "SM", "N"), function(vals){
            sen_mk_func_slope(as.numeric(tempdf$max_Q_date), tempdf[,vals])
          }))
          out = data.frame(t(theil_all))
        }else{
          out = data.frame(t(rep(NA, 12)))
        }
        
        
        
        out[,"gauge_id"] =  tempdf[!is.na(tempdf$gauge_id),"gauge_id"][1]
        out[,"lag_time"] =  tempdf[!is.na(tempdf$gauge_id),"lag_time"][1]
        out[,"meanSM_Q"] =  mean(tempdf[!is.na(tempdf$gauge_id),"SMquant"], na.rm = T)
        
      }else{
        out = NULL
      }
      return(out)
    })
  })
})


trend_quantilesdfP7 = df_processing_func(trend_quantiles_listP7[[1]])



dfQP = data.frame(gauge_id = trend_quantilesdfQ$gauge_id, P7Qslope = trend_quantilesdfQ$Pvolume_slope, QPQslope = trend_quantilesdfQ$max_Q_slope, arid_dir = trend_quantilesdfQ$arid_dir)
dfPP = data.frame(gauge_id = trend_quantilesdfP7$gauge_id, P7Pslope = trend_quantilesdfP7$Pvolume_slope, PPQslope = trend_quantilesdfP7$max_Q_slope)
compare_trendP7 = merge(dfQP, dfPP)
compare_trendP7all = compare_trendP7
compare_trendP7all[,"arid_dir"] = rep("All", nrow(compare_trendP7all))
compare_trendP7 = rbind(compare_trendP7, compare_trendP7all)


#Compare Q trends for all catchments
pall = ggplot(compare_trendP7[compare_trendP7$arid_dir == "All",], aes(QPQslope, PPQslope))+geom_point(alpha = 0.2)+geom_abline(slope = 1, intercept = 0, col = "red")+xlim(-65,65)+ylim(-65,65)+xlab(paste0("Streamflow trend select by Q>Q(", 100*quant_thresh,")"))+ylab(paste0("Streamflow trend select by P7>P7(",100*quant_thresh,")"))+theme_minimal(base_size = 13)+theme(legend.position = "none", text = element_text(size = 11), strip.text = element_text(size = 11))
ggMarginal(pall, type="histogram", binwidth = 2, colour = "gray50", fill = "gray90")

pdf(file = paste0(plot_path, "/CAMELS_trendcompareP7_AllpointsQ99.pdf"), width = 5, height = 4 )
ggMarginal(pall, type="histogram", binwidth = 2, colour = "gray50", fill = "gray90")
dev.off()




# Supplement: different splits besides aridity ----------------------------
load(file = "E:/phd/phd_20210618/Data/CAMELS/clim_df_list.Rdata")
load(file = "E:/phd/phd_20210618/Data/CAMELS/camels_all_loctopo.Rdata")


Prec_dfallclimS = merge(Prec_dfall, camels_all_clim[,c("gauge_id", "p_seasonality", "pet_mean")], by = "gauge_id", all.x = T)
Prec_dfallclimS = merge(Prec_dfallclimS, camels_all_loctopo[,c("gauge_id", "area")], by = "gauge_id", all.x = T)
Prec_dfallclimS[,"season_dir"] = cut(Prec_dfallclimS$p_seasonality, breaks=c(-1.5,-0.5, 0.5, 1.5), right = FALSE)
Prec_dfallclimS[,"area_dir"] = cut(Prec_dfallclimS$area, breaks=c(0, 10, 100, 1000, 10000, 100000), right = FALSE)
Prec_dfallclimS[,"pet_dir"] = cut(Prec_dfallclimS$pet_mean, breaks=c(1,2,3,4,8), right = FALSE)



var.labs_area <- c("(a) Cat size 0-10 km2", "(b) Cat size 10-100 km2" , "(c) Cat size 100-1e+03 km2" ,  "(d) Cat size 1e+03-1e+04 km2", "(e) Cat size 1e+04-1e+05 km2")
names(var.labs_area) <- levels(unique(Prec_dfallclimS$area_dir))

var.labs_season <- c("(a) Seasonality = (-1.5)-(-0.5)", "(b) Seasonality = (-0.5)-0.5" , "(c) Seasonality = 0.5-1.5")
names(var.labs_season) <- levels(unique(Prec_dfallclimS$season_dir))

#1-day precipitation 

contour_dflist_Pmax = lapply(unique(Prec_dfallclimS$area_dir), function(area_ind){
  testpts = Prec_dfallclimS[Prec_dfallclimS$area_dir == area_ind,]
  
  # set up an 'empty' raster, here via an extent object derived from your data
  r <- raster(xmn=0, xmx=1, ymx=1, ymn=0.5, ncol=100, nrow=50) #changing limit to 1.1 does not change the outcome
  r <- rasterize(testpts[, c("SMquant", "Pquant")], r, testpts[,"Qquant"], fun=median)
  r2 <- focal(r, w= matrix(1,5,5), mean)
  r2df = data.frame(rasterToPoints(r2))
  r2df[,"area_dir"] = rep(area_ind, nrow(r2df))
  
  return(r2df)
})
contour_df = do.call(rbind, contour_dflist_Pmax)


contour_df100 = contour_df
contour_df100[,1:3] = contour_df100[,1:3]*100
ggplot(contour_df100, aes(x, y, fill = layer)) +  
  geom_tile() +
  geom_contour(aes(z = layer), col = "white", binwidth = 5)+
  geom_contour(aes(z = layer), col = "black", breaks = 90, lty = 2)+
  scale_fill_viridis_c(expression(paste("Med(Pr(Q<Qmax))"," [%]")), limits=c(0,100), direction = -1)+geom_text_contour(aes(z = layer), col = "white")+theme_minimal(base_size = 13)+xlab(expression(paste("Pr(SM<SMI)"," [%]")))+ylab(expression(paste("Pr(P<Pmax)"," [%]")))+facet_wrap(~area_dir, labeller = labeller(area_dir = var.labs_area), nrow = 1 )+theme(strip.text = element_text(size = 13), legend.position = "bottom", legend.margin=margin(-10, 0, 0, 0))
ggsave(paste0(plot_path, "/CAMELS_SMPQquant_contour_area.pdf"), width = 14, height = 4 , dpi = "retina")

#Seasonality

contour_dflist_Pmax = lapply(unique(Prec_dfallclimS$season_dir), function(season_ind){
  testpts = Prec_dfallclimS[Prec_dfallclimS$season_dir == season_ind,]
  
  # set up an 'empty' raster, here via an extent object derived from your data
  r <- raster(xmn=0, xmx=1, ymx=1, ymn=0.5, ncol=100, nrow=50) #changing limit to 1.1 does not change the outcome
  r <- rasterize(testpts[, c("SMquant", "Pquant")], r, testpts[,"Qquant"], fun=median)
  r2 <- focal(r, w= matrix(1,5,5), mean)
  r2df = data.frame(rasterToPoints(r2))
  r2df[,"season_dir"] = rep(season_ind, nrow(r2df))
  
  return(r2df)
})
contour_df = do.call(rbind, contour_dflist_Pmax)


contour_df100 = contour_df
contour_df100[,1:3] = contour_df100[,1:3]*100
ggplot(contour_df100, aes(x, y, fill = layer)) +  
  geom_tile() +
  geom_contour(aes(z = layer), col = "white", binwidth = 5)+
  geom_contour(aes(z = layer), col = "black", breaks = 90, lty = 2)+
  scale_fill_viridis_c(expression(paste("Med(Pr(Q<Qmax))"," [%]")), limits=c(0,100), direction = -1)+geom_text_contour(aes(z = layer), col = "white")+theme_minimal(base_size = 13)+xlab(expression(paste("Pr(SM<SMI)"," [%]")))+ylab(expression(paste("Pr(P<Pmax)"," [%]")))+facet_wrap(~season_dir, labeller = labeller(season_dir = var.labs_season),  nrow = 1 )+theme(strip.text = element_text(size = 13), legend.position = "bottom", legend.margin=margin(-10, 0, 0, 0))
ggsave(paste0(plot_path, "/CAMELS_SMPQquant_contour_season.pdf"), width = 9, height = 4 , dpi = "retina")



#Q90/Q50 subsurface storage index 
#(Borga et al., 2007; Norbiato et al., 2009; Tarasova et al.,2018)

system.time({
  Q90Q50_list = lapply(c(1:length(clim_df_list)), function(dfnum){
    print(dfnum)
    climdf = clim_df_list[[dfnum]]
    if(is.null(climdf)){
      Q90Q50 = NULL
    }else{
      NA_test = apply(climdf[,c("Prec", "SM", "Q")], 2, function(x){(sum(is.na(x))/length(x))>0.5})
      if(any(NA_test)){
        Q90Q50 = NULL
      }else{
        ratiotemp = quantile(climdf$Q, 0.1, na.rm = T)/quantile(climdf$Q, 0.5, na.rm = T)
        Q90Q50 = data.frame(gauge_id = climdf[1,"gauge_id"], Q90Q50 = ratiotemp)
      }
    }
    
    return(Q90Q50)
  })
  
})

Q90Q50df = do.call(rbind, Q90Q50_list)

Prec_dfallclimS = merge(Prec_dfall, Q90Q50df, by = "gauge_id", all.x = T)
Prec_dfallclimS[,"storage_dir"] = cut(Prec_dfallclimS$Q90Q50, breaks=seq(0,0.8,0.2), right = FALSE)


var.labs_storage <- c("(a) Q90/Q50 = 0-0.2",  "(b) Q90/Q50 = 0.2-0.4", "(c) Q90/Q50 = 0.4-0.6", "(d) Q90/Q50 = 0.6-0.8")
names(var.labs_storage) <- levels(unique(Prec_dfallclimS$storage_dir))


contour_dflist_Pmax = lapply(unique(Prec_dfallclimS$storage_dir), function(area_ind){
  testpts = Prec_dfallclimS[Prec_dfallclimS$storage_dir == area_ind,]
  
  # set up an 'empty' raster, here via an extent object derived from your data
  r <- raster(xmn=0, xmx=1, ymx=1, ymn=0.5, ncol=100, nrow=50) #changing limit to 1.1 does not change the outcome
  r <- rasterize(testpts[, c("SMquant", "Pquant")], r, testpts[,"Qquant"], fun=median)
  r2 <- focal(r, w= matrix(1,5,5), mean)
  r2df = data.frame(rasterToPoints(r2))
  r2df[,"storage_dir"] = rep(area_ind, nrow(r2df))
  
  return(r2df)
})
contour_df = do.call(rbind, contour_dflist_Pmax)


contour_df100 = contour_df
contour_df100[,1:3] = contour_df100[,1:3]*100
ggplot(contour_df100, aes(x, y, fill = layer)) +  
  geom_tile() +
  geom_contour(aes(z = layer), col = "white", binwidth = 5)+
  geom_contour(aes(z = layer), col = "black", breaks = 90, lty = 2)+
  scale_fill_viridis_c(expression(paste("Med(Pr(Q<Qmax))"," [%]")), limits=c(0,100), direction = -1)+geom_text_contour(aes(z = layer), col = "white")+theme_minimal(base_size = 13)+xlab(expression(paste("Pr(SM<SMI)"," [%]")))+ylab(expression(paste("Pr(P<Pmax)"," [%]")))+facet_wrap(~storage_dir, labeller = labeller(storage_dir = var.labs_storage),  nrow = 1 )+theme(strip.text = element_text(size = 13), legend.position = "bottom", legend.margin=margin(-10, 0, 0, 0))
ggsave(paste0(plot_path, "/CAMELS_SMPQquant_contour_storage.pdf"), width = 14, height = 4 , dpi = "retina")



