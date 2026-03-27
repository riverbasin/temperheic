library(zoo)
library(lubridate)
library(temperheic)
library(RColorBrewer)

filePrefix = "fluxBoundaryTesto."
# datadir = "C:/Users/Byron/Documents/HGS/"
datadir = "/home/amerson/"

directories = c(
  paste0(datadir,"noBedrock100k_Copy/")) #,
  # paste0(datadir,"surfaceOnly100k/"),
  # paste0(datadir,"bedrockOnly100k/"),
  # paste0(datadir,"bedrockSurface100k/"),
  # paste0(datadir,"noBedrock400k/"),
  # paste0(datadir,"surfaceOnly400k/"),
  # paste0(datadir,"bedrockOnly400k/"),
  # paste0(datadir,"bedrockSurface400k/")
# )

#directories = c(paste0(datadir,"sandBox_surfaceAndSubHeatFlowWithBedrock_k400mday/"))

fileNamePrefix = "fluxBoundaryTesto.observation_well_conc.well_"
fileNameSuffix = ".temperature.dat"


distances = c(
  0,
  10,
  20,
  30,
  40,
  50,
  60,
  70,
  80,
  90,
  100,
  200,
  300,
  400,
  500,
  600,
  700,
  800,
  900,
  1000
)

fpLengths = c("temp0" = 0,
              "temp10" = 10,
              "temp20" = 20,
              "temp30" = 30,
              "temp40" = 40,
              "temp50" = 50,
              "temp60" = 60,
              "temp70" = 70,
              "temp80" = 80,
              "temp90" = 90,
              "temp100" = 100,
              "temp200" = 200,
              "temp300" = 300,
              "temp400" = 400,
              "temp500" = 500,
              "temp600" = 600,
              "temp700" = 700,
              "temp800" = 800,
              "temp900" = 900,
              "temp1000" = 1000)

#create a thUnits Object
units = thUnits()

#create a thAquifer object
aquifer = thAquifer(0.25, 0.0016, 0.000598, 0.84, 4.186, 3000, 998, units) ##thermCond_sed from  Burns et al. 2015

### plot HGS modeled time series using time series derived with plotTemps function

tempTimeSeries1 = mapply(
  modelTempsSubsurfaceDates,
  directory = directories,
  MoreArgs = list(
    distances = distances,
    fileNamePrefix = fileNamePrefix,
    fileNameSuffix = fileNameSuffix,
    mdyHMS = "6/1/2012 00:00:00"
  ),
  SIMPLIFY = F
)

tempTimeSeriesWindow = lapply(tempTimeSeries1, window,
                              start = mdy_hms("06/01/2015 00:00:00"),
                              end = mdy_hms("01/30/2017 01:00:00"))

# modelData_1D = tempTimeSeriesWindow[[5]]

tempTS_noTop = list(tempTimeSeriesWindow[[1]])

# plot(analyticData_1D, plot.type = "single", lwd = 2)
# lines(modelData_1D$temp0, col = "green")
# lines(modelData_1D$temp1000, col = "green")

# n <- 20
# qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
# col_vector = unlist(mapply(brewer.pal,
#                            qual_col_pals$maxcolors,
#                            rownames(qual_col_pals)))
# par(mfcol = c(4, 2))
# par(oma = c (3, 3, 2, 2))
# par(mar = c(2, 2.2, .15, .1))
# lapply(tempTimeSeriesWindow,
#        plot,
#        plot.type = "single",
#        lwd = seq(4,.7,-.17),
#        col=sample(col_vector, n))

empiricalDataPeriodList = list(
  c(rep(1, 11),1,1,1,1,2,2,2,2,2)) #,
  # c(rep(1, 11),1,1,1,2,2,2,2,2,2),
  # rep(1, 20),
  # rep(1, 20),
  # rep(1, 20),
  # rep(1, 20)
# )

# balls = thObservedSeries(modelData_1D,
#                  xVals = fpLengths,
#                  period = 365*86400,
#                  aquifer = aquifer,
#                  nmin = 180,
#                  specificUnits = units,
#                  laggedLinearFit = F,
#                  headGrad = 0.01,
#                  optimizeRange = c(0,1),
#                  empiricalDataPeriods = c(rep(1, 20))
#                  )


modelInverseSolutions = mapply(
  thObservedSeries,
  tempTS_noTop,
  empiricalDataPeriods = empiricalDataPeriodList,
  MoreArgs = list(
    xVals = fpLengths,
    boundaryMean = 10.5,
    aquifer = aquifer,
    period = 365 * 86400,
    headGrad = 0.01,
    nmin = 1000,
    specificUnits = attr(aquifer, "specificUnits"),
    optimizeRange = c(0, 1)
  ),
  SIMPLIFY = F
)




disps =  lapply(lapply(modelInverseSolutions, '[[', i = "dispersivity"),
                '[',
                1,
                T)
disps = lapply(disps, round, 3)


hydConds =  lapply(lapply(modelInverseSolutions, '[[', i = "hydraulicCond"),
                   '[',
                   1,
                   T)

hydConds = lapply(hydConds, '*', 86400)
hydConds = lapply(hydConds, round, 2)

