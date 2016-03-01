#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(ggplot2)
library(doBy)

# Define server logic required to draw a histogram
shinyServer(function(input, output) {
  setwd("C:/Users/Notandi/Dropbox/Projects/simulation/")
  source("definitions.r")
  source("utility function poisson.r")
  source("generatebasegen.r")
  source("runsimulation.r")
  
  setwd("C:/Users/Notandi/Dropbox/Projects/simulation of mink farm/Output/DMU analysis/")
  con <- file(description = "results", open = "w")
  cat(
    "Gen",
    "Gmean",
    "Gvar",
    "Fis",
    "Obs.fert",
    "mean.phenotype.bs.females",
    "gen.value.bs",
    "mean.phenotype.bs.males",
    "bw.var",
    sep = "\t",
    file = con
  )
  cat("\n", file = con)
  close(con = con)
  for (p in 1:nruns) {
    # runcounter <- p
    # if (make.obs.file == 1) {
    #   output <- file(description = paste("Replicate_",p, sep=""), open="w")
    #   phenotypes <- file(description = paste("Phenotypes",p, sep=""), open="w")
    # }
    year <- 1
    l <- RunFirstYear(p, year)
    for (y in 1:n) {
      year <- 1 + y
      l <- RunSimulation(l, year, p)
      
    }
    # if (make.obs.file == 1) {
    #   close(con=output)
    #   close(con = phenotypes)
    # }
    print("reached end")
    if (p == nruns) {
    results <- read.table("results", header = TRUE)
    summarized <-
      summaryBy(
        Gen + Gmean + Gvar + Fis + Obs.fert + mean.phenotype.bs.females + gen.value.bs +
          mean.phenotype.bs.males + bw.var
        ~ Gen,
        data = results,
        FUN = c(mean, var),
        na.rm = T
      )
    
  output$plot1 <-
    reactivePlot(function() {
      pl1 <- (ggplot(data = summarized, aes(x=Gen, y =Gmean.mean,fill = Gmean.mean))
            + geom_line(linetype =2)
            + geom_errorbar(
              aes(
                x = Gen,
                ymin = Gmean.mean - sqrt(Gmean.var),
                ymax = Gmean.mean + sqrt(Gmean.var)
              ),
              width = 0.25,
              data = summarized
            ))
      print(pl1)
    })
  output$plot2 <-
    reactivePlot(function() {
      pl2 <- (ggplot(data = summarized, aes(x=Gen, y =Obs.fert.mean))
              + geom_line(linetype =2)
              + geom_errorbar(
                aes(
                  x = Gen,
                  ymin = Obs.fert.mean - sqrt(Obs.fert.var),
                  ymax = Obs.fert.mean + sqrt(Obs.fert.var)
                ),
                width = 0.25,
                data = summarized
              ))
      print(pl2)
    })
  
  }
      }
  closeAllConnections()
  
  # results <- read.table("results", header = TRUE)
  
})
