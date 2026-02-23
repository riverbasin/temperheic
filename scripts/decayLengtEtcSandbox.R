freq = aquiferBoundaryList$annual[[1]]$frequency  ## 1/P in seconds
ke = aquiferHydroList$annual[[2]]$diffusivity_effective
vt = aquiferHydroList$annual[[2]]$advectiveThermVel
w = (2*pi)*freq
pw = aquiferHydroList$annual[[2]]$aquifer$density_h2o
cw = aquiferHydroList$annual[[2]]$aquifer$spHeat_h2o
q = aquiferHydroList$annual[[2]]$darcyFlux
k = aquiferHydroList$annual[[2]]$aquifer$thermCond_bulk
zdDiff = aquiferSignalList$annual[[1]]$thermDecayDist
zdVd = aquiferSignalList$annual[[2]]$phaseVel/w
zdAdv = aquiferSignalList$annual[[2]]$thermDecayDist


num1 <- 64*(pi^2)*(freq^2)*(ke^2)
rad1 <- (1 + (num1/(ke^4)))^0.5
rad2 <- (0.5 + (0.5*rad1))^0.5
phaseVel <- vt * rad2
thermDecayDist <- (2*ke)/(phaseVel-vt)

rad1A = (vt^4 + (4*w*ke)^2)^0.5
rad2A = ((rad1A + vt^2)/2)^0.5 - vt
aye = (1/(2*ke)) * rad2A

rad22A = sqrt(2)/(((rad1A + vt^2))^0.5 - sqrt(2)*vt)
aye2 = 2*ke * rad22A

rad1B = (vt^4 + (4*w*ke)^2)^0.5
rad2B = ((rad1B - vt^2)/2)^0.5
bee = (1/(2*ke)) * rad2B

rad22B = sqrt(2)/((rad1B - vt^2))^0.5
bee2 = 2*ke * rad22B

a = (1 / (2 * ke)) * (sqrt(((sqrt((vt^4) + (4 * 2 * pi * freq * ke)^2) + vt^2)/2)) - vt)
b = (1 / (2 * ke)) * sqrt(((sqrt((vt^4) + (4 * 2 * pi * freq * ke)^2) - vt^2)/2))


PeT = (pw*cw*q) /  ke
peTLuce = vt / ke


Dref = aquiferSeriesList$annual[[1]]$signal$hydro$diffusivity_cond
Ds = aquiferHydroList$annual[[1]]$aquifer$thermCond_sed/aquiferHydroList$annual[[1]]$aquifer$volHeatCap_sed
Df = aquiferHydroList$annual[[1]]$aquifer$thermCond_h2o/aquiferHydroList$annual[[1]]$aquifer$volHeatCap_h2o

d50 = 0.005 ## pretend D50 for Umatilla 5mm converted to meters
R = ((1-.25) * aquiferSeriesList$annual[[1]]$signal$hydro$aquifer$volHeatCap_sed) / (.25 * aquiferSeriesList$annual[[1]]$signal$hydro$aquifer$volHeatCap_h2o)
peRef = ((aquiferSeriesList$annual[[1]]$signal$hydro$darcyFlux/.25) * d50) / Dref
vc = aquiferSeriesList$annual[[1]]$signal$hydro$diffusivity_cond/d50
peCrit = (vc * d50) / aquiferSeriesList$annual[[1]]$signal$hydro$diffusivity_cond
t = 1  #tortuosity, usually larger than 1
Dh = 1/(1+R) * (Dhs*R Dh + (R+1)*Dref  )
