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
  n =4,
  nruns = 15,
  selection.method = 2 ,
  mblup = 1,
  make.obs.file = 1,
  cull.ratio = 0.85,
  purebreeding = 0,
  quantile.setting.ls = 0.4,
  quantile.setting.bw =0.4,
  prop.oldfemales = 0.4,
  max.age.females = 3,
  crossmating = 0,
  n.cages = 7000,
  n.females = 2000,
  use.true.sire = 0,
  pseudo.import = 0,
  use.blup.to.assort.mat = 1,
  use.comb.ind.for.males = 1
  
  # weight.fert.old.females = 0.5,
  # weight.bw.old.females   = 0.250,
  # weight.qual.old.females = 0.25,
  # weight.fert.kits        = 0.10,
  # weight.bw.kits          = 0.5,
  # weight.qual.kits        = 0.4
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

