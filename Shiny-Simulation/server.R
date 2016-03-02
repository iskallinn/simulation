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

setwd("C:/Users/Notandi/Dropbox/Projects/simulation/")
source("definitions.r")
source("utility function poisson.r")
source("generatebasegen.r")
source("runsimulation.r")
source("simulation_function.r")


# Define server logic required to draw a histogram
shinyServer(function(input, output) {
  source("definitions.r")
  
sim <- reactive({
  Simulation(nruns = input$n_runs )
})

# sim ()
  
 
   output$plot1 <-
     reactivePlot(function() {
       if (input$run == 0) return ()
       # browser()
       sim ()
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
       })
  
 
 
  closeAllConnections()
  
  # 
  
})  
