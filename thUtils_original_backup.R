

#' @export
as.xts.thSeries = function(x, POSIXct.origin) {
  theOrder = as.POSIXct(x$time, origin = POSIXct.origin)
  x = x[-match("time", names(x))]
  x = xts(x, order.by = theOrder)
  return(x)
}

#' @export
as.zoo.thSeries = function(x, POSIXct.origin) {
  theOrder = as.POSIXct(x$time, origin = POSIXct.origin)
  x = x[-match("time", names(x))]
  x = zoo(x, order.by = theOrder)
  return(x)
}

#' @export
htPlot = function(myThSeries, POSIXct.origin = "2014-01-01 00:00:00") {
  myXTS = as.xts.thSeries(myThSeries, POSIXct.origin)
  plot(myXTS$x0, main = "", ylab = "Temperature (degC)", type = 'n', auto.grid = F, minor.ticks = F)
  lines(myXTS$x0, col = "grey70", lwd = 3, lty=2)
  for(col in names(myXTS)[2:length(names(myXTS))]) {
    lines(myXTS[,col], lwd = 3)
  }
  seriesAttr = attributes(myThSeries)
  sg = attr(myThSeries, "signal")
  hy = attr(sg, "hydro")
  k = hy$hydCond
  d = hy$dispersivity
  text(as.numeric(min(index(myXTS))), max(myXTS$x0), pos = 4, paste0("k = ", round(k,6), "  B = ", round(d, 6)))
}



# Replaces general units (e.g., "E t-1 L-3 T-1") with specific units specified
# by the user (e.g., "kJ s-1 m-3 degC-1")
.thSpecificUnits = function(generalUnits, specificUnits = thUnits()) {
  atomicUnits = unlist(lapply(generalUnits, strsplit, split = " "), recursive = F)

  #regular expression for an optional "-" and any number of digits at the end of
  #a string
  pattern = "[-]?[[:digit:]]+$"
  #strip off any match to pattern
  generalSymbols = lapply(atomicUnits, sub, pattern = pattern, replacement = "")
  #return any substring that matches to pattern
  generalExponents =
    lapply(
      atomicUnits,
      function(x) {
        locations = regexpr(pattern = pattern, x)
        locations[locations == -1] = 1000000L
        substring(x, locations)
      }
    )

  #return any generalSymbols that are not expected.  If any, throw error and
  #report to user
  unexpectedGenUnits = sapply(generalSymbols, function(x) {any(!(x %in% names(specificUnits)))})
  if(any(unexpectedGenUnits)) {
    cat("Unexpected general units found in: ", paste0("'", generalUnits[unexpectedGenUnits], "'", collapse = ", "), '\n  Expected units are: ', paste0(attributes(specificUnits)$longUnitName, "='", names(specificUnits), "'", collapse = ", "))
    stop("A general unit can be followed by positive or negative number (representing and exponent).\nThere can be no space between a unit and its exponent (e.g., 'M-3' not 'M -3').\nThere must be spaces between unit/exponent pairs (e.g. 'E L-3', not 'EL-3')")
  }

  #replace the generalSymbols with corresponding specific ones
  generalSymbols = lapply(generalSymbols, function(x) unlist(specificUnits)[x])
  #concatinate specific units with exponents and return vector specific units and exponents
  return(mapply(paste0, generalSymbols, generalExponents, collapse = " "))
}

#' Mean Squared Residuals and linear regression between two lagged time series
#'
#' Functions for regression two time series against one another (y ~ x), while
#' accounting for any lag time of the pattern in y realative to similar patterns
#' in x.
#'
#' When two time series of the same length are subjected to a lag, the number of
#' x,y pairs is reduced in proportion to the size of the lag because the series
#' are become more and more offset in time (analogous to the reduction of the
#' amount of overlap between two meter-sticks, which start out aligned, but are
#' then then slid in opposite directions).  The nmin ensures that the estimate
#' of mean squared residuals is based on at least nmin x,y pairs, once the time
#' lag in y is accounted for.
#'
#' \code{laggedMSR()} is desigend to be passed to optimize() in order to find
#' the lag with the minimum mean squared residuals between time series x and y.
#'
#' @return \code{laggedMSR} Returns the mean squared residuals of a linear model
#'   (y ~ x) given a time series x and y, assuming that time series y lags time
#'   series x by lag time units.  NOTE: When regressed agains one another, two
#'   cos waves with a lag of pi radians will yield a mean squared residual of
#'   zero and a slope of -1.0.  This is an undesirable solution.  The prefered
#'   solution is a lag of zero, which will yield a MSR of 0 and a slope of 1.0.
#'   This, in this function, residuals are calculated using the absolute value
#'   of the regression slope. This ensures that the prefered solution (where low
#'   MSR is associated with lags that are in phase rather than out of phase) is
#'   always returned.
#' @param lag Time that time series y lags time series x.
#' @param thSeriesPair Zoo object with two columns -- the starting sine wave and ending sine wave
#' @param t A vector of times of observations of values in x and y
#' @param nmin Minimum number of x.y pairs desired (see Details)
#' @export
laggedMSR = function(lag, thSeriesPair, nmin) {
  lData = laggedData(lag, thSeriesPair)
  if(nrow(lData) < nmin) {
    result = NULL
  } else {
    fit = laggedModel(lData)
    a = coefficients(fit)
    result = mean((lData[,2] - (a[1] + abs(a[2]) * lData[,1]))^2)
  }
  return(result)
}

#' @rdname laggedMSR
#' @return \code{laggedModel} runs \code{\link{lm}()} on lData and returns
#'   the results.  Usually, lData is generated by calling \code{laggedData()}
#' @param lData A zoo object, typcally returned by calling \code{laggedData()}
#' @export
laggedModel = function(lData) {
  return(lm(lData[,2] ~ lData[,1]))
}

#'@rdname laggedMSR
#'@return \code{laggedData()} creates a \code{\link{zoo}} object with two
#'  columns (x and y).  Each row in the zoo object contains a pair of
#'  observations, after accounting for the lag -- the amount of time y lags x.
#'  \code{\link{na.spline}()} is used to calculate the y column in the zoo
#'  object if lag is not an even multiple of the times between observations.
#'@export
laggedData = function(lag, thSeriesPair) {
  # Create zoo objects
  thSeriesX = thSeriesPair[,1]
  thSeriesY = thSeriesPair[,2]
  zoo::index(thSeriesY) = zoo::index(thSeriesY) - lag
  # Merge series into one object
  storeNames = names(thSeriesPair)
  thSeriesPair <- merge(thSeriesX, thSeriesY)
  names(thSeriesPair) = storeNames
  # Interpolate calibration data (na.spline could also be used)
  thSeriesPair[,2] <- zoo::na.approx(object = thSeriesPair[,2], na.rm = F)
  # Only keep index values from sample data
  thSeriesPair <- thSeriesPair[!(is.na(thSeriesPair[,1]) | is.na(thSeriesPair[,2])) ,]
  return(thSeriesPair)
}

# converts passed values into 3d array
derivedArray = function(ampRatio, deltaPhaseRadians, eta, seriesNames) {
  derivedVals = array(
    data = c(ampRatio, deltaPhaseRadians, eta),
    dim = c(length(seriesNames), length(seriesNames), 3),
    dimnames = list(from = seriesNames, to = seriesNames, value = c("ampRatio", "deltaPhaseRadians", "eta"))
  )
  return(derivedVals)
}

derived2DArray = function(x, seriesNames) {
  derivedVals = array(
    data = x,
    dim = c(length(seriesNames), length(seriesNames)),
    dimnames = list(from = seriesNames, to = seriesNames)
  )
  return(derivedVals)
}


fitCosine = function(empiricalData, boundaryMean, periodInSeconds, optimizeRange, nmin, empiricalDataPeriods) {

  # means = sapply(empiricalData, mean)
  # grandMean = mean(means)
  ## using nls to fit the data
  ampEst = sapply(empiricalData, function(x) (max(x) - min(x))/2)
  seconds = as.numeric(zoo::index(empiricalData))

  fits = mapply(
    function(eD, ampEst, obsTime, boundaryMean, period){
      if(length(na.omit(eD)) >= nmin) {
        # put in teeny tiny bit of noise so nls will fit a cosine to data that represent a perfect wave
        eD = eD + runif(length(eD), min = -(mean(eD) * 0.0001), max = mean(eD)* 0.0001)
        fitResult = nls(
          formula = eD ~ boundaryMean + AmpY._ * cos((2*pi/period) * obsTime - PhaY._),
          start = list(AmpY._ = ampEst, PhaY._ = 0),
          na.action = "na.exclude"
        )
        # if the fit amplitude is negative, negate it and add pi radians to the phase
        # modulus phase by 2pi to make sure the phase is between 0 and 2pi radians
        if(coef(fitResult)[1] < 0){
          results = c(fitAmp = -coef(fitResult)[1], fitPhase = (coef(fitResult)[2] + pi)%%(2*pi))
        }else{
          results = c(fitAmp = coef(fitResult)[1], fitPhase = coef(fitResult)[2]%%(2*pi))
        }

      } else {
        results = c(fitAmp = NA, fitPhase = NA)
      }
      return(results)
    },
    eD = as.data.frame(empiricalData),
    ampEst = ampEst,
    MoreArgs = list(
      obsTime = seconds - seconds[1],
      boundaryMean = boundaryMean,
      period = periodInSeconds
    ),
    SIMPLIFY = T
  )

  fitAmp = fits[1,]
  fitPhase = fits[2,]
  initialPhase = fitPhase[1]
  relativeRange = 2*pi*optimizeRange + initialPhase
  fitPhase[fitPhase < relativeRange[1]] = fitPhase[fitPhase < relativeRange[1]] + 2*pi
  fitPhase = fitPhase + (empiricalDataPeriods - 1)*2*pi
  #outOfRange = fitPhase > optimizeRange[[2]]*2*pi + fitPhase[1]
  #fitPhase[which(outOfRange)] = fitPhase[which(outOfRange)] - 2*pi
  ## find permutative combintation of differences for fitPhases
  permuteCombos = expand.grid(1:length(fitPhase), 1:length(fitAmp))
  deltaPhaseRadians = apply(permuteCombos, 1, function(x) fitPhase[x[2]] - fitPhase[x[1]])
  ampRatio = apply(permuteCombos, 1, function(x) fitAmp[x[2]] / fitAmp[x[1]])

  # 1.0000000  1.1128460         NA  0.8985648  1.0000000         NA  0.5578810  0.6208431  1.0000000
  #1.666759e-16 -2.262049e-01            NA  2.262049e-01  1.563496e-14            NA  1.234432e+00  1.008227e+00 -1.000026e-16

  return(structure(list(deltaPhaseRadians = deltaPhaseRadians, ampRatio = ampRatio), amplitudes = fitAmp, phases = fitPhase))
}


lagLinFit <- function(empiricalData, periodInSeconds, optimizeRange, nmin) {
  combos = expand.grid(from = names(empiricalData), to = names(empiricalData), stringsAsFactors = F)
  # calculate phases using phase shifting approach
  dphase = t(
    sapply(
      1:nrow(combos),
      function(row) {
        if(row %in% diagonalLocs) {
          result = list(minimum = 0, objective = 0)
        } else {
          result = optimize(f = laggedMSR, interval = periodInSeconds * optimizeRange, thSeriesPair = empiricalData[,as.character(combos[row,])], nmin = nmin)
        }
        return(result)
      }
    )
  )

  dphase[dphase[,"objective"] < 0,"minimum"] = NA
  dphase = unlist(dphase[,"minimum"])
  deltaPhaseRadians = dphase*2*pi/periodInSeconds

  ampRatio = sapply(
    1:nrow(combos),
    function(rowIndex) {
      if(is.na(dphase[rowIndex])) {
        result = NA
      } else {
        result = laggedData(dphase[rowIndex], empiricalData[,as.character(combos[rowIndex,])])
        result = laggedModel(result)
        result = coefficients(result)[2]
      }
      return(result)
    }
  )

  return(list(deltaPhaseRadians = deltaPhaseRadians, ampRatio = ampRatio))

}
