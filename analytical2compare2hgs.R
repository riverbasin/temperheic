library(tidyverse)
library(ggplot2)
library(lubridate)
library(cowplot)


# myUnits = thUnits("ft", "lb", "day", "degF", "ft-lb")
# myUnits = thUnits(t = "d")
# myUnits = "T"
library(zoo)
library(temperheic)
Sys.setenv(TZ="UTC")

#create a thUnits Object
myUnits = thUnits()

#create a thAquifer object
myAquifer = thAquifer(0.25, 0.0016, 0.000598, 0.84, 4.186, 3000, 998, myUnits) ##thermCond_sed from  Burns et al. 2015

#permutations for hydro objects
k = annual = c("100k" = 100/86400, "400k" = 400/86400)
disp = c(1)
grad = 0.01

#permutations for boundary objects. Derived from nls fit in umatilla_temps.r
mean = annual = 10.5
amp = c(9.6)
phase = c(18169179-16632000)
period = 86400*365

#locations and times for signal plots
xVals = c(0,10,20,30,40,50,60,70,80,90,100,200,300,400,500,600,700,800,900,1000)
tVals = seq(0, (5*365*86400)-86400, 86400)
#tVals = seq(0, (5*365*24*3600)-3600, 3600)

hydro = pmap(list(k, list(disp), list(grad), list(myAquifer), list(myUnits)),
           thHydro)
boundary = pmap(list(list(mean), list(amp), list(phase), list(period), list(myUnits)),
                thBoundary)
signal = pmap(list(hydro, boundary),
              thSignal)
series = pmap(list(signal, list(xVals), list(tVals), list(myUnits)),
              thSeries)
analyticData_1D = series$timeSeries
index(analyticData_1D) = index(analyticData_1D) + mdy_hms("06/01/2015 00:00:00")
analyticData_1D = window(analyticData_1D,
       start = mdy_hms("06/01/2015 00:00:00"),
       end = mdy_hms("01/30/2017 01:00:00"))

## create plot function using ggplot and zoo to make TS plots
plotF = function(tempTS){
  ggplot(aes(x = Index,
             y = Value,
             group = Series,
             colour = Series),
         data = fortify(tempTS, melt = TRUE)) + geom_line() +
    theme_cowplot(8) +
    scale_colour_viridis_d(name  ="Flow Path",
                           breaks=c("temp0", "temp10",
                                    "temp20", "temp30",
                                    "temp40", "temp50",
                                    "temp60", "temp70",
                                    "temp80", "temp90",
                                    "temp100", "temp200",
                                    "temp300", "temp400",
                                    "temp500", "temp600",
                                    "temp700", "temp800",
                                    "temp900", "temp1000"),
                           labels=c("0 m", "10 m",
                                    "20 m", "30 m",
                                    "40 m", "50 m",
                                    "60 m", "70 m",
                                    "80 m", "90 m",
                                    "100 m", "200 m",
                                    "300 m", "400 m",
                                    "500 m", "600 m",
                                    "700 m", "800 m",
                                    "900 m", "1000 m"))
}

## generate a list of plot objects
tsPlots = lapply(tempTS_noTop, plotF)

# use cowplot function plotgrids to create facet plot of timeseries
plotGrid = plot_grid(
  tsPlots[[1]] + xlab(NULL) + ylab("Temperature (C)") + theme(legend.position="none"),
  tsPlots[[5]] + xlab(NULL) + ylab(NULL) + theme(legend.position="none"),
  tsPlots[[2]] + xlab(NULL) + ylab("Temperature (C)") + theme(legend.position="none"),
  tsPlots[[6]] + xlab(NULL) + ylab(NULL) + theme(legend.position="none"),
  tsPlots[[3]] + xlab("Time")+ ylab("Temperature (C)") + theme(legend.position="none"),
  tsPlots[[7]] + xlab("Time") + ylab(NULL) + theme(legend.position="none"),
  tsPlots[[4]] + xlab("Time")+ ylab("Temperature (C)") + theme(legend.position="none"),
  tsPlots[[8]] + xlab("Time") + ylab(NULL) + theme(legend.position="none"),
  nrow=4,
  ncol = 2,
  align = "vh",
  labels = c("US Only K=100 m/d",
             "US Only K=400 m/d",
             "US & Top K=100 m/d",
             "US & Top K=400 m/d",
             "US & BR K=100 m/d",
             "US & Br K=400 m/d",
             "All Bounds K=100 m/d",
             "All Bounds K=400 m/d"),
  label_size = 8,
  label_fontface = "plain",
  label_x = 0,
  label_y = .13,
  hjust = -0.3)


# extract the legend from one of the plots
legend <- get_legend(
  # create some space to the left of the legend
  tsPlots[[1]] +
    theme(legend.box.margin = margin(0, 0, 0, 12)))

# again use cowplot to combine faceted timeseries and legend into a plot
# note very slow plot time and not very satisfying plot labels
# check
# http://www.sthda.com/english/articles/24-ggpubr-publication-ready-plots/81-ggplot2-easy-way-to-mix-multiple-graphs-on-the-same-page/
# for other faceting options that may have better annotation capability

png(file="C:/Users/Byron/Desktop/forward2D.png",
    res = 300,
    width=2100,
    height=1500)

plot_grid(plotGrid, legend, rel_widths = c(3, .4))

dev.off()



ggplot(aes(x = Index,
           y = Value,
           group = Series),
       data = fortify(analyticData_1D, melt = TRUE)) +
  geom_line(aes(colour = Series)) +
  theme_cowplot(8) +
  scale_colour_viridis_d(name  ="Flow Path",
                        breaks=c("x0", "x10",
                                 "x20", "x30",
                                 "x40", "x50",
                                 "x60", "x70",
                                 "x80", "x90",
                                 "x100", "x200",
                                 "x300", "x400",
                                 "x500", "x600",
                                 "x700", "x800",
                                 "x900", "x1000"),
                        labels=c("0 m", "10 m",
                                 "20 m", "30 m",
                                 "40 m", "50 m",
                                 "60 m", "70 m",
                                 "80 m", "90 m",
                                 "100 m", "200 m",
                                 "300 m", "400 m",
                                 "500 m", "600 m",
                                 "700 m", "800 m",
                                 "900 m", "1000 m"))

testObsSeries = thObservedSeries(analyticData_1D,
                                 xVals = c(x0 = 0,
                                           x10 = 10,
                                           x20=20,
                                           x30=30,
                                           x40=40,
                                           x50=50,
                                           x60=60,
                                           x70=70,
                                           x80=80,
                                           x90=90,
                                           x100=100,
                                           x200=200,
                                           x300=300,
                                           x400=400,
                                           x500=500,
                                           x600=600,
                                           x700=700,
                                           x800=800,
                                           x900=900,
                                           x1000=1000),
                                 period = 365*86400,
                                 aquifer = myAquifer,
                                 nmin = 180,
                                 specificUnits = myUnits,
                                 laggedLinearFit = F,
                                 headGrad = 0.01,
                                 optimizeRange = c(0,1),
                                 empiricalDataPeriods = c(rep(1, 11),1,1,1,2,2,2,2,2,2))
