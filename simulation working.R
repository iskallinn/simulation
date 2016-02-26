############## Connection for output ##############
setwd("C:/Users/Notandi/Dropbox/Projects/simulation/")
source("definitions.r")
source("utility function poisson.r")
source("generatebasegen.r")
source("runsimulation.r")

setwd("C:/Users/Notandi/Dropbox/Projects/simulation of mink farm/Output/DMU analysis/")
con <- file(description="results",open="w")
cat("Gen","Gmean","Gvar","Fis","Obs.fert", "mean.phenotype.bs.females", "gen.value.bs", "mean.phenotype.bs.males","bw.var",sep="\t",file=con)
cat("\n",file=con)

for (p in 1:nruns) {
if (make.obs.file == 1) {
  output <- file(description = paste("Replicate_",p, sep=""), open="w")
  phenotypes <- file(description = paste("Phenotypes",p, sep=""), open="w")
}
year <- 1
l <- RunFirstYear()
for (y in 1:n) {
  year <- 1+y
  l <-RunSimulation(l)
  
} 
if (make.obs.file == 1) {
  close(con=output) 
  close(con = phenotypes)
}
}
closeAllConnections()
