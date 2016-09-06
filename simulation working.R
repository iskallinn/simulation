############## Change to installed folder containing the scripts ##############
setwd("C:/Users/Notandi/Dropbox/Projects/simulation/")
setwd("C:/Users/au384062/Documents/GitHub/simulation") #AU workstation
# source("first_time.r")
source("definitions.r")
source("generatebasegen.r")
source("utility_functions.r")
source("runsimulation.r")
source("simulation_function.r")


Simulation( n = 2, nruns = 1, selection.method = phenotypic ,mblup=0, make.obs.file = 0)

