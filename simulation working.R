############## Change to installed folder containing the scripts ##############
setwd("C:/Users/Notandi/Dropbox/Projects/simulation/")
setwd("C:/Users/au384062/Documents/GitHub/simulation") #AU workstation
# source("first_time.r")
source("definitions.r")
source("generatebasegen.r")
source("utility_functions.r")
source("runsimulation.r")
source("simulation_function.r")



Simulation(
  n = 10,
  nruns = 200,
  selection.method = phenotypic ,
  mblup = 0,
  make.obs.file = 0,
  cull.ratio = 0.85,
  purebreeding = 0,
  quantile.setting.ls = 0.4,
  quantile.setting.bw,
  crossmating = 1
  )

