Sys.setenv(TZ='UTC')
library(temperheic)
library("lubridate")
library("xts")
library("zoo")

dir = ("~/Documents/Meacham Creek/documents/Manuscripts/Chapter 1/plot code/")
umatilla_data = read.csv(paste0(dir,"umatilla_raw_data_2.csv"), head = TRUE)
umatilla_data = data.frame(Time = umatilla_data$Time, Temp = umatilla_data$Temp, SiteName = umatilla_data$SiteName,
                           Status = umatilla_data$Status, ActionIdx = umatilla_data$ActionIdx, ActionType = umatilla_data$ActionType)
umatilla_data$Time = mdy_hm(as.character(umatilla_data$Time))
umatilla_data$SiteName = as.character(umatilla_data$SiteName)
site_names = unique(umatilla_data$SiteName)

temp1 <- umatilla_data$Temp[umatilla_data$SiteName == site_names[1]]
temp2 <- umatilla_data$Temp[umatilla_data$SiteName == site_names[2]]
temp3 <- umatilla_data$Temp[umatilla_data$SiteName == site_names[3]]

time1 <- umatilla_data$Time[umatilla_data$SiteName == site_names[1]]
time2 <- umatilla_data$Time[umatilla_data$SiteName == site_names[2]]
time3 <- umatilla_data$Time[umatilla_data$SiteName == site_names[3]]

temp1.zoo = zoo(temp1, order.by = time1)
temp2.zoo = zoo(temp2, order.by = time2)
temp3.zoo = zoo(temp3, order.by = time3)

temp1.xts <- as.xts(temp1.zoo)
temp1.xts.daily <- period.apply(temp1.xts,INDEX=endpoints
                                (index(temp1.xts), "days"),FUN=mean)
foo.bar <- floor_date(index(temp1.xts.daily), "day")
index(temp1.xts.daily) <- foo.bar
temp1.zoo.daily <- as.zoo(temp1.xts.daily)

temp2.xts <- as.xts(temp2.zoo)
temp2.xts.daily <- period.apply(temp2.xts,INDEX=endpoints
                                (index(temp2.xts), "days"),FUN=mean)
foo.bar <- floor_date(index(temp2.xts.daily), "day")
index(temp2.xts.daily) <- foo.bar
temp2.zoo.daily <- as.zoo(temp2.xts.daily)

temp3.xts <- as.xts(temp3.zoo)
temp3.xts.daily <- period.apply(temp3.xts,INDEX=endpoints
                                (index(temp3.xts), "days"),FUN=mean)
foo.bar <- floor_date(index(temp3.xts.daily), "day")
index(temp3.xts.daily) <- foo.bar
temp3.zoo.daily <- as.zoo(temp3.xts.daily)
temp.3.daily.approx = na.approx(temp3.zoo.daily)

umatillaTemps = cbind(temp.3.daily.approx, temp2.zoo.daily, temp1.zoo.daily)
#colnames(umatillaTemps) <-c("channel", "well175m", "well955m")
colnames(umatillaTemps) <-c(0, 175, 955)

#plot(umatillaTemps, plot.type = "single")
empirical.stream = window(na.approx(umatillaTemps$"0"),
                          start = ymd_hms("2003-07-13 00:00:00"),
                          end = ymd_hms("2004-07-12 00:00:00"))
empirical.175 = window(na.approx(umatillaTemps$"175"),
                       start = ymd_hms("2003-07-13 00:00:00"),
                       end = ymd_hms("2004-07-12 00:00:00"))
empirical.955 = window(na.approx(umatillaTemps$"955"),
                       start = ymd_hms("2003-07-13 00:00:00"),
                       end = ymd_hms("2004-07-12 00:00:00"))

umatillaTemps2003to2004 = window(na.approx(umatillaTemps),
                                 start = ymd_hms("2003-07-14 00:00:00"),
                                 end = ymd_hms("2004-07-12 00:00:00"))

#plot(umatillaTemps2003to2004, plot.type = "single")

thObservedSeriesUmatilla = thObservedSeries(umatillaTemps2003to2004,
                                            xVals = c("0" = 0, "175" = 175, "955" = 955),
                                            aquifer = aquiferHydroList$annual[[1]]$aquifer,
                                            hydro = aquiferHydroList$annual[[1]],
                                            period = 365*86400,
                                            headGrad = 0.008,
                                            nmin = 180,
                                            specificUnits = attr(aquiferHydroList$annual[[1]]$aquifer,
                                                                 "specificUnits"),
                                            laggedLinearFit = F)
