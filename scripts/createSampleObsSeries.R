# myUnits = thUnits("ft", "lb", "day", "degF", "ft-lb")
# myUnits = thUnits(t = "d")
# myUnits = "T"

Sys.setenv(TZ="UTC")

#create a thUnits Object
myUnits = thUnits()

#create a thAquifer object
myAquifer = thAquifer(0.25, 0.00837, 0.000598, 0.84, 4.186, 2650, 998, myUnits)

#permutations for hydro objects
k = list(annual = c(10,100,1000)/86400, daily = c(.1, 1, 10)/86400)
disp = list(annual = c(10, 100, 1000, 10000), daily = c(0.01, 0.1, 1, 10))
grad = list(annual = 0.01, daily = 0.1)

#permutations for boundary objects
mean = list(annual = 10, daily = 22)
amp = list(annual = 8, daily = 6)
phase = list(annual = 734040, daily = 86400*3/5)
period = list(annual = 86400*365, daily = 86400)

#locations and times for signal plots
xVals = list(annual = c(0, 100, 200, 500, 1000, 1500, 2000), daily = c(0, 2, 3, 6, 8))
tVals = list(annual = seq(0, (730*86400), 86400), daily = seq(0, (3600*48), 3600))

# set these lists to the index of the singal you want to plot.  Permutations are
# in the "sensitivity" object calculated next
plotList = list(annual = T, daily = T)

sensitivity = mapply(function (x, y, z) structure(expand.grid(x,y,z), names = c("disp", "k", "grad")), disp, k, grad, SIMPLIFY = F)


#____________________________________________________

permuteHydros = function(sa) {
  thHydroList = (mapply(function (k, disp, grad) thHydro(k, disp, grad, myAquifer, myUnits), sa$k, sa$disp, sa$grad, SIMPLIFY = F))
  return(thHydroList)
}

aquiferHydroList = sapply(sensitivity, permuteHydros, simplify = F)

sensitivity = mapply(function (w, x, y, z) structure(expand.grid(w, x,y,z), names = c("mean", "amp", "phase", "period")), mean, amp, phase, period, SIMPLIFY = F)

permuteBoundarys = function(sa) {
  thBoundaryList = (mapply(function (mean, amp, phase, period) thBoundary(mean, amp, phase, period, myUnits), sa$mean, sa$amp, sa$phase, sa$period, SIMPLIFY = F))
  return(thBoundaryList)
}

aquiferBoundaryList = sapply(sensitivity, permuteBoundarys, simplify = F)

nBoundary = sapply(aquiferBoundaryList, length)
nHydro = sapply(aquiferHydroList, length)

sensitivity = mapply(function(nH, nB) structure(expand.grid(1:nH, 1:nB), names = c("hydro", "boundary")), nHydro, nBoundary, SIMPLIFY = F)

permuteSignals = function(sa, hydroList, boundaryList) {
  thSignalList = (mapply(function (hIdx, bIdx) thSignal(hydroList[[hIdx]], boundaryList[[bIdx]]), sa$hydro, sa$boundary, SIMPLIFY = F))
  return(thSignalList)
}

aquiferSignalList = mapply(permuteSignals, sensitivity, aquiferHydroList, aquiferBoundaryList, SIMPLIFY = F)

permuteSeries = function(sig, x, t) {
  seriesList = sapply(sig, thSeries, xVals = x, tVals = t, specificUnits = myUnits, simplify = F)
}

aquiferSeriesList = mapply(permuteSeries, aquiferSignalList, xVals, tVals, SIMPLIFY = F )

# lapply(aquiferSeriesList$annual[plotList$annual], htPlot)
# lapply(aquiferSeriesList$daily[plotList$daily], htPlot)

chosenSeries = 3

testObsSeries = thObservedSeries(empiricalData = aquiferSeriesList$annual[[chosenSeries]]$timeSeries[,1:4],
                                 xVals = c(x0 = 0, x100 = 100, x200 = 200, x500 = 500),
                                 boundaryMean = aquiferSeriesList$annual[[chosenSeries]]$signal$boundary$mean,
                                 period = 365*86400,
                                 aquifer = myAquifer,
                                 nmin = 180,
                                 specificUnits = myUnits,
                                 laggedLinearFit = F,
                                 headGrad = 0.01,
                                 optimizeRange = c(0,1),
                                 empiricalDataPeriods = c(1,1,1,1))
