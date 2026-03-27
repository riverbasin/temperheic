#' @export
generateExampleThAquifer <- function(
  porosity = 0.25,
  thermCond_sed = 0.0016, #Columbia River basalt from Burns et al. 2015 in GeoFluids v. 15
  thermCond_h2o = 0.000598,
  spHeat_sed = 0.84, #Columbia River basalt from Burns et al. 2015 in GeoFluids v. 15
  spHeat_h2o = 4.186,
  density_sed = 3000, #Columbia River basalt from Burns et al. 2015 in GeoFluids v. 15
  density_h2o = 998,
  specificUnits = thUnits()
){
  #create a thAquifer object
  thAquifer(porosity, thermCond_sed, thermCond_h2o, spHeat_sed, spHeat_h2o, density_sed, density_h2o, specificUnits)
}

#' @export
generateExampleThSeries <-
  function(
    myAquifer,
    k = list(annual = c(100)/86400, daily = c(10)/86400),
    disp = list(annual = c(100, 1000, 10000), daily = c(0.001, 0.01, 0.1, 1)),
    grad = list(annual = 0.01, daily = 0.1),
    mean = list(annual = 11.7798, daily = 22),
    amp = list(annual = c(8.869333), daily = 6),
    phase = list(annual = c(1702226), daily = 86400 + 3600),
    period = list(annual = 86400*365, daily = 86400),
    xVals = list(annual = c(0, 120, 240, 360), daily = c(0, .33,.66, 1)),
    tVals = list(annual = seq(0, (2*365*86400), 86400), daily = seq(0, (3600*48), 3600)),
    specificUnits = thUnits()
  ){

    storedSysEnv <- Sys.getenv("TZ")
    Sys.setenv(TZ="UTC")

    sensitivity = mapply(function (x, y, z) structure(expand.grid(x,y,z), names = c("disp", "k", "grad")), disp, k, grad, SIMPLIFY = F)

    #____________________________________________________

    permuteHydros = function(sa) {
      thHydroList = (mapply(function (k, disp, grad) thHydro(k, disp, grad, myAquifer, specificUnits), sa$k, sa$disp, sa$grad, SIMPLIFY = F))
      return(thHydroList)
    }

    aquiferHydroList = sapply(sensitivity, permuteHydros, simplify = F)

    sensitivity = mapply(function (w, x, y, z) structure(expand.grid(w, x,y,z), names = c("mean", "amp", "phase", "period")), mean, amp, phase, period, SIMPLIFY = F)

    permuteBoundarys = function(sa) {
      thBoundaryList = (mapply(function (mean, amp, phase, period) thBoundary(mean, amp, phase, period, specificUnits), sa$mean, sa$amp, sa$phase, sa$period, SIMPLIFY = F))
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
      seriesList = sapply(sig, thSeries, xVals = x, tVals = t, specificUnits = specificUnits, simplify = F)
    }

    aquiferSeriesList = mapply(permuteSeries, aquiferSignalList, xVals, tVals, SIMPLIFY = F )

    Sys.setenv(TZ = storedSysEnv)

    return(aquiferSeriesList)

  }

#' @export
generateExampleThObservedSeries <- function(
  aquiferSeriesList,
  nmin = c(annual = 180, daily = 12),
  periodOptimizeRange = c(-1/8, 7/8),
  laggedLinearFit = T
)
{
  makeThObservedSeries <- function(x, nmin, periodOptimizeRange)
  {
    empiricalData = x$timeSeries
    xVals = x$xVals
    boundaryMean = x$signal$boundary$mean
    period = x$signal$boundary$period
    hydro = x$signal$hydro
    freq = (2*pi)/period
    optimizeRange = periodOptimizeRange
    specificUnits = attr(x, "specificUnits")
    thObservedSeries(empiricalData, xVals, hydro$aquifer, boundaryMean, period, hydro$headGrad, nmin, freq, optimizeRange, specificUnits, laggedLinearFit)
  }

  mapply(
    function(annualOrDailyList, annualOrDailyNmin){
      mapply(
        makeThObservedSeries,
        annualOrDailyList,
        MoreArgs = list(periodOptimizeRange = periodOptimizeRange, nmin = annualOrDailyNmin),
        SIMPLIFY = F
      )
    },
    aquiferSeriesList,
    nmin
  )
}

#' @export
testForwardBackward <- function(laggedLinearFit = T)
{
  thSeriesList = generateExampleThSeries(generateExampleThAquifer())
  thObservedSeriesList = generateExampleThObservedSeries(thSeriesList, laggedLinearFit = laggedLinearFit)

  matrixListNames <- lapply(
    thSeriesList,
    function(x) {
      lapply(
        x,
        function(x2) {
          c(
            paste0("K = ", x2$signal$hydro$hydCond),
            paste0("dispersivity = ", x2$signal$hydro$dispersivity)
          )
        }
      )
    }
  )

  matrixList <- mapply(
    function(x, n) {
      mapply(
        function(x2, n2) {
          empVals <- list(
            KEmpirical = x2$hydraulicCond,
            dispersivityEmpirical = x2$dispersivity
          )
          names(empVals) <- n2
          return(empVals)
        },
        x,
        n,
        SIMPLIFY = F
      )
    },
    thObservedSeriesList,
    matrixListNames,
    SIMPLIFY = F
  )

  return(matrixList)

}

