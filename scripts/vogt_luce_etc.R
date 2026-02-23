
#####################################################################
#some tests

## daily conductive thermal decay distance
aquiferSeriesList$daily[[1]]$signal$thermDecayDist_cond

## annual conductive thermal decay distance
aquiferSeriesList$annual[[4]]$signal$thermDecayDist_cond

## daily conductive thermal decay distance * 19
aquiferSeriesList$daily[[1]]$signal$thermDecayDist_cond * 19

## annual conductive-advective thermal decay distance
z = aquiferSeriesList$annual[[4]]$zdAdvective

## amplitude at thermal decay distance
8.869333/exp(1)

## thSeries amplitude at thermal decay distance
(max(aquiferSeriesList$annual[[4]]$timeSeries[,4]) - min(aquiferSeriesList$annual[[4]]$timeSeries[,4]))/2
max(aquiferSeriesList$annual[[4]]$timeSeries$x955) - mean(aquiferSeriesList$annual[[4]]$timeSeries$x955)

t = aquiferSeriesList$annual[[1]]$signal$thermDecayDist / aquiferSeriesList$annual[[1]]$signal$phaseVel

P = aquiferSeriesList$annual[[1]]$signal$boundary$period
f = aquiferSeriesList$annual[[1]]$signal$boundary$frequency
w = 2*pi*f

z = aquiferSeriesList$annual[[1]]$zdAdvective

vt <- aquiferSeriesList$annual[[1]]$signal$hydro$advectiveThermVel
ke <- aquiferSeriesList$annual[[1]]$signal$hydro$diffusivity_effective

# after Luce 2013;
vd = aquiferSeriesList$annual[[1]]$vdAdvective

# after Vogt 2014; i.e. c
num1 <- 64*(pi^2)*(f^2)*(ke^2)
rad1 <- (1 + (num1/(vt^4)))^0.5
rad2 <- (0.5 + (0.5*rad1))^0.5
vd2 <- vt * rad2

z1 <- (2*ke)/(vd2-vt)
