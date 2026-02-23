
library(numDeriv)
signal = aquiferSignalList$daily[[4]]
tvalsDaily = 54000 #seq(0, 172800, length.out = 731)-seq(0, 172800, length.out = 731)[1]
func0 = function(x, tvals) {
  signal$boundary$mean + signal$boundary$amplitude * exp(-x/signal$thermDecayDist) * cos((2 * pi / signal$boundary$period) * ((tvals) - signal$boundary$phase - (x/signal$phaseVel)))
}
func1 = function(tvals, x){grad(func0, x = x, tvals = tvals)}
distance = seq(0,5, length.out = 100)
plot(distance, unlist(lapply(distance, func1, tvals = tvalsDaily)), col = "red", type = "l")
# lines(unlist(lapply(tvalsDaily, func1, x = .33)), col = "green")
# lines(unlist(lapply(tvalsDaily, func1, x = 5)), col = "blue")

signal = aquiferSignalList$annual[[5]]
tvalsAnn = yday(ymd('2012-8-15'))*86400 #tVals$annual-tVals$annual[2]
func0 = function(x, tvals) {
  signal$boundary$mean + signal$boundary$amplitude * exp(-x/signal$thermDecayDist) * cos((2 * pi / signal$boundary$period) * ((tvals) - signal$boundary$phase - (x/signal$phaseVel)))
}
func1 = function(tvals, x){grad(func0, x = x, tvals = tvals)}
#lines(unlist(lapply(tvalsAnn, func1, x = 120)))
distance = seq(0, 5000, length.out = 100)
plot(distance, unlist(lapply(distance, func1, tvals = tvalsAnn)), col = "green", type = "l")
# lines(unlist(lapply(tvalsAnn, func1, x = 120)), col = "green")
# lines(unlist(lapply(tvalsAnn, func1, x = 2000)), col = "blue")

###### trying out some Peclet numbers
#### assume d50 of 65mm or 0.065 m

PeVelocity = (0.0065*aquiferSeriesList$annual[[1]]$signal$hydro$velocity_h2o)/aquiferSeriesList$annual[[1]]$signal$hydro$diffusivity_cond
PeThermVel = (0.0065*aquiferSeriesList$annual[[1]]$signal$hydro$advectiveThermVel)/aquiferSeriesList$annual[[1]]$signal$hydro$diffusivity_cond

# or using diffusive decay depth a la Luce or Goto
zd = aquiferSeriesList$annual[[1]]$signal$thermDecayDist
(zd*aquiferSeriesList$annual[[4]]$signal$hydro$velocity_h2o)/aquiferSeriesList$annual[[4]]$signal$hydro$diffusivity_cond

### trying out Bons stuff
### there are two unknowns: disperivity and the exponent
#
## from meacham data we know ke needs to be about 0.001158744 and we know tha kcond is 4.6e-07
## sum of kdisp and kcond = 0.0012, so kcond + disp * Pe^n must = 0.0012

## we can either set disp to solve for exponent, or alternatively set exponent to solve for disp
disp = 0.0001
ke = c(0.0011592042, 0.0025806947, 0.0015080708, 0.0009331538, 0.0045749676, 0.0014509814, 0.0020399125,
       0.0128507885, 0.0062586404, 0.0084686338, 0.0166172580)
bonsExponent = log(ke/disp)/log(Pe)

estDisp = (ke/Pe^1.2)-aquiferSeriesList$annual[[4]]$signal$hydro$diffusivity_cond

Pe = c(1e-1, 1e0, 1e1, 1e2, 1e3, 1e4)
n = c(.0001, .3, 1.1, 1.05, 1.01, 1.0)
aquiferSeriesList$annual[[4]]$signal$hydro$diffusivity_cond + 0.0001*Pe^n

v = (Pe*aquiferSeriesList$annual[[5]]$signal$hydro$diffusivity_cond)/.05
KA = (v*.25)/0.01
KD = (v*.25)/0.1

## classic Pe solutions
PeVelocity = (0.0065*aquiferSeriesList$annual[[1]]$signal$hydro$velocity_h2o)/aquiferSeriesList$annual[[1]]$signal$hydro$diffusivity_cond
PeThermVel = (0.0065*aquiferSeriesList$annual[[1]]$signal$hydro$advectiveThermVel)/aquiferSeriesList$annual[[1]]$signal$hydro$diffusivity_cond
## our updated Pe using effective Ke and thermal decay length

meachamPe = thObservedSeriesMeacham$pecletNumber[1,][which(thObservedSeriesMeacham$dispersivity[1,] > 0 & thObservedSeriesMeacham$dispersivity[1,] < 1000)]

thObservedSeriesUmatilla$pecletNumber

(thObservedSeriesUmatilla$advectiveThermVelEmpirical*0.0065)/aquiferSeriesList$annual[[1]]$signal$hydro$diffusivity_cond

meachamKe = thObservedSeriesMeacham$diffusivity_effective_empirical[1,][which(thObservedSeriesMeacham$dispersivity[1,] > 0 & thObservedSeriesMeacham$dispersivity[1,] < 1000)]
meachamBons = log(meachamKe/0.003)/log(meachamPe)
log(meachamKe/0.003)/log(PeVelocity)
log(meachamKe/aquiferSeriesList$annual[[1]]$signal$hydro$diffusivity_cond)/log(meachamPe)


estDisp = (meachamKe/meachamPe^exponentBons)/aquiferSeriesList$annual[[1]]$signal$hydro$diffusivity_cond


# using exponents based on Pe to solve for dispersion coeff to comapre to our solution from Ke
exponentBons = c(.1, .2, .3, .4, .5, .6, .7, .8, .9, 1, 1.1, 1.2, 1.1, 1)
pecletNums = c(.1,.2,.4,.6,.8,1,4,6,8,10,20,100,200, 500)
cond = aquiferSeriesList$annual[[1]]$signal$hydro$diffusivity_cond
Deff = exp(exponentBons*log(pecletNums)+log(cond))
Deff/cond

