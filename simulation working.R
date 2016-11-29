############## Change to installed folder containing the scripts ##############
setwd("C:/Users/Notandi/Dropbox/Projects/simulation/")
setwd("C:/Users/au384062/Documents/GitHub/simulation") #AU workstation
# source("first_time.r")
source("definitions.r")
source("generatebasegen.r")
source("utility_functions.r")
source("runsimulation.r")
source("simulation_function.r")


# library(profvis)
# profvis({ 
  Simulation(
  n =5,
  nruns = 1,
  selection.method = 2 ,
  mblup = 1,
  make.obs.file = 0,
  cull.ratio = 0.85,
  purebreeding = 0,
  quantile.setting.ls = 0.4,
  quantile.setting.bw =0.4,
  prop.oldfemales = 0.4,
  max.age.females = 3,
  crossmating = 0,
  n.cages = 3500,
  n.females = 1000,
  use.true.sire = 0,
  pseudo.import = 0,
  use.blup.to.assort.mat = 0,
  use.comb.ind.for.males = 0,
  genetic.means=c(
    0,                  # live.qual
    0,                  # h.length
    0,                  # skin.qual
    0,                  # skin.length.male
    0,                  # skin.length.female
    0.01,               # litter.size
    0,              # body weight females
    0,              # body weight males
    0,     # rfi1.m
    0,           # rfi2.m
    0,     # rfi1.f
    0      # rfi2.f
  )
  
  
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
# }
# )