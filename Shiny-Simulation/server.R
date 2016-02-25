#
# This is the server logic of a Shiny web application. You can run the 
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
   
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
    # runcounter <- p
    if (make.obs.file == 1) {
      output <- file(description = paste("Replicate_",p, sep=""), open="w")
      phenotypes <- file(description = paste("Phenotypes",p, sep=""), open="w")
    }
    year <- 1
    l <- generate.gen0(p)
    for (y in 1:n) {
      year <- 1+y
      l <-run.simulation(l)
      
    } 
    if (make.obs.file == 1) {
      close(con=output) 
      close(con = phenotypes)
    }
  }
  closeAllConnections()
  
  results <- read.table("results", header = TRUE)
  output$plot1 <- renderPlot(results$Gen)
})
