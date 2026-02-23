##Characteristic time stuff

((100^2)/posKe)/86400
((10^2)/posKe)/86400

((meachamThermDecayDist^2)/posKe)/86400
((2^2)/posKe)/86400

##meacham
(0.01*aquiferSeriesList$annual[[1]]$signal$hydro$velocity_h2o)/aquiferSeriesList$annual[[1]]$signal$hydro$diffusivity_cond
(0.01*aquiferSeriesList$annual[[1]]$signal$hydro$advectiveThermVel)/aquiferSeriesList$annual[[1]]$signal$hydro$diffusivity_cond
(0.01*mean(thObservedSeriesMeacham$advectiveThermVelEmpirical[1,][which(thObservedSeriesMeacham$dispersivity[1,] > 0 & thObservedSeriesMeacham$dispersivity[1,] < 1000)]))/aquiferSeriesList$annual[[1]]$signal$hydro$diffusivity_cond
(0.01*mean(thObservedSeriesMeacham$velocity_h2o[1,][which(thObservedSeriesMeacham$dispersivity[1,] > 0 & thObservedSeriesMeacham$dispersivity[1,] < 1000)]))/aquiferSeriesList$annual[[1]]$signal$hydro$diffusivity_cond

##umatilla
(thObservedSeriesUmatilla$advectiveThermVelEmpirical*0.0065)/aquiferSeriesList$annual[[1]]$signal$hydro$diffusivity_cond
(thObservedSeriesUmatilla$velocity_h2o*0.0065)/aquiferSeriesList$annual[[1]]$signal$hydro$diffusivity_cond
thObservedSeriesUmatilla$pecletNumber

kay = 10
vel = ((kay/86400)*.01)/.25
d50 = 6.5/1000
(d50 * vel)/4.6e-7

dispRatio = posKe/aquiferSeriesList$annual[[1]]$signal$hydro$diffusivity_cond

Pe = (aquiferSeriesList$annual[[1]]$signal$thermDecayDist*aquiferSeriesList$annual[[1]]$signal$hydro$advectiveThermVel)/aquiferSeriesList$annual[[1]]$signal$hydro$diffusivity_cond
log(dispRatio)/log(Pe)
10^3

