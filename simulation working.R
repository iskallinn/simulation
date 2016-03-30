

############## Connection for output ##############
setwd("C:/Users/Notandi/Dropbox/Projects/simulation/")
source("definitions.r")
source("utility function poisson.r")
source("generatebasegen.r")
source("runsimulation.r")
source("simulation_function.r")

system.time(Simulation( n = 3, nruns = 1, use.true.sire = 1))

closeAllConnections()
