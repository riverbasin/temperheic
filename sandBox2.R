# myUnits = thUnits("ft", "lb", "day", "degF", "ft-lb")
# myUnits = thUnits(t = "d")
# myUnits = "T"

Sys.setenv(TZ="UTC")

#create a thUnits Object
myUnits = thUnits()

#create a thAquifer object
myAquifer = thAquifer(0.25, 0.00837, 0.000598, 0.84, 4.186, 3000, 998, myUnits)

#permutations for hydro objects
k = list(annual = c(100)/86400)
disp = list(annual = c(1))
grad = list(annual = 0.01)

#permutations for boundary objects
mean = list(annual = 10.5)
amp = list(annual = 8)
phase = list(annual = 734040)
period = list(annual = 86400*365)

#locations and times for signal plots
xVals = list(annual = c(0,10,20,30,40,50,60,70,80,90,100,200,300,400,500))
tVals = list(annual = seq(0, (730*86400), 86400))

sensitivity = mapply(function (x, y, z) structure(expand.grid(x,y,z), names = c("disp", "k", "grad")), disp, k, grad, SIMPLIFY = F)

permuteHydros = function(sa) {
  thHydroList = (mapply(function (k, disp, grad) thHydro(k, disp, grad, myAquifer, myUnits), sa$k, sa$disp, sa$grad, SIMPLIFY = F))
  return(thHydroList)
}

aquiferHydroList = sapply(sensitivity, permuteHydros, simplify = F)

R = (0.75*(myAquifer$spHeat_sed*myAquifer$density_sed))/(.25*(myAquifer$spHeat_h2o*myAquifer$density_h2o))
Ds = myAquifer$thermCond_sed/(myAquifer$spHeat_sed*myAquifer$density_sed)
Df = myAquifer$thermCond_h2o/(myAquifer$spHeat_h2o*myAquifer$density_h2o)
Dref = aquiferHydroList$annual[[1]]$diffusivity_cond
betaH = 1.33
L = c(1)
PeH = (aquiferHydroList$annual[[1]]$advectiveThermVel*L)/aquiferHydroList$annual[[1]]$diffusivity_cond
PeC = 12
Disp = (1/(R+1)) * (Ds*R + (Df + betaH*(R+1)*Dref*PeH*(PeH/PeC))/(1+(PeH/PeC)))
Deff = aquiferHydroList$annual[[1]]$diffusivity_effective
Disp
Deff
(Disp - aquiferHydroList$annual[[1]]$diffusivity_cond)/aquiferHydroList$annual[[1]]$advectiveThermVel



plot(PeH[2:9], Disp[2:9]/Disp[1], log = "xy")
lines(PeH[2:9], Disp[2:9]/Disp[1])


(Disp - aquiferHydroList$annual[[1]]$diffusivity_cond)/aquiferHydroList$annual[[1]]$advectiveThermVel


R = (0.75*(myAquifer$spHeat_sed*myAquifer$density_sed))/(.25*(myAquifer$spHeat_h2o*myAquifer$density_h2o))
Ds = myAquifer$thermCond_sed/(myAquifer$spHeat_sed*myAquifer$density_sed)
Df = myAquifer$thermCond_h2o/(myAquifer$spHeat_h2o*myAquifer$density_h2o)
Dref = aquiferHydroList$annual[[1]]$diffusivity_effective
betaH = 1.33
L = 1
PeH = (aquiferHydroList$annual[[1]]$advectiveThermVel*L)/aquiferHydroList$annual[[1]]$diffusivity_effective
PeC = 12
Disp = (1/(R+1)) * (Ds*R + (Df + betaH*(R+1)*Dref*PeH*(PeH/PeC))/(1+(PeH/PeC)))
