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
  n =10,
  nruns = 200,
  selection.method = 1 ,
  mblup = 0,
  make.obs.file = 0,
  cull.ratio = 0.85,
  purebreeding = 0,
  quantile.setting.ls = 0.4,
  quantile.setting.bw =0.4,
  prop.oldfemales = 0.6,
  max.age.females = 6,
  crossmating = 0
  # true.sire.chance = 0.85,
  # mating.will.yearling.1st   =0.95, # probability of yearling being mated
  # mating.will.yearling.2nd   =0.80, # probability of yearling being remated
  # mating.will.old.1st        =0.98, # probability of old female being mated
  # mating.will.old.2nd        =0.84, # probability of old female being remated
  # pr.barren.one.mating.yearling     = 0.73, # probability of single mated yearling being barren
  # pr.barren.double.mating.yearling  = 0.82, # probability of double mated yearling being barren
  # pr.barren.one.mating.old   =0.80, # probability of single mated old female being barren
  # pr.barren.double.mating.old=0.87
  )

