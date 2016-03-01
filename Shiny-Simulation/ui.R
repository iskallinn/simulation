#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
# 
#    http://shiny.rstudio.com/
#

library(shiny)

# Define UI for application that draws a histogram
shinyUI(fluidPage(

  # Application title
  titlePanel("Mink Sim 0.02"),

  # Sidebar with a slider input for number of bins
  sidebarLayout(
    sidebarPanel(
      numericInput("n",
                   "Number of generations in each replicate:",
                   min = 1,
                   max = 15,
                   value = 10),
      numericInput("n_runs",
                   "Number of Replicates",
                   min = 1,
                   max = 15,
                   value = 10),
      numericInput("n.females",
                   "Number females at farm",
                   min = 100,
                   max = 3000,
                   value = 1000),
      numericInput("prop.oldfemales",
                   "Proportion of old females on farm",
                   min = 0,
                   max = 1,
                   value = 0.4),
      # numericInput("weight.fert.old.females",
      #              "Weight on fertility index if selection index is used, NOTE the weights must sum to 1",
      #              min = 0,
      #              max = 1,
      #              value = 0.35),
      # numericInput("weight.bw.old.females",
      #              "Weight on body weight index if selection index is used,NOTE the weights must sum to 1",
      #              min = 0,
      #              max = 1,
      #              value = 0.65),
      checkboxGroupInput("selection.method", 
                         label = h3("Selection method"), 
                         choices = list("Selection" = 1, 
                                        "BLUP" = 2),
                         selected = 1),
      submitButton()
      
      
      
    ),

    # Show a plot of the generated distribution
    mainPanel(
      tabsetPanel(
      tabPanel("Genetic trend, litter size",plotOutput("plot1"))
      ,
      tabPanel("Litter size, 1st counting",plotOutput("plot2"))
    )
    )
  )
)) 

# renderInputs <- function(prefix) {
#   wellPanel(
#     fluidRow(
#       column(6,
#              sliderInput(paste0(prefix, "_", "n_obs"), "Number of observations (in Years):", min = 0, max = 40, value = 20),
#              sliderInput(paste0(prefix, "_", "start_capital"), "Initial capital invested :", min = 100000, max = 10000000, value = 2000000, step = 100000, pre = "$", sep = ","),
#              sliderInput(paste0(prefix, "_", "annual_mean_return"), "Annual investment return (in %):", min = 0.0, max = 30.0, value = 5.0, step = 0.5),
#              sliderInput(paste0(prefix, "_", "annual_ret_std_dev"), "Annual investment volatility (in %):", min = 0.0, max = 25.0, value = 7.0, step = 0.1)
#       ),
#       column(6,
#              sliderInput(paste0(prefix, "_", "annual_inflation"), "Annual inflation (in %):", min = 0, max = 20, value = 2.5, step = 0.1),
#              sliderInput(paste0(prefix, "_", "annual_inf_std_dev"), "Annual inflation volatility. (in %):", min = 0.0, max = 5.0, value = 1.5, step = 0.05),
#              sliderInput(paste0(prefix, "_", "monthly_withdrawals"), "Monthly capital withdrawals:", min = 1000, max = 100000, value = 10000, step = 1000, pre = "$", sep = ","),
#              sliderInput(paste0(prefix, "_", "n_sim"), "Number of simulations:", min = 0, max = 2000, value = 200)
#       )
#     ),
#     p(actionButton(paste0(prefix, "_", "recalc"),
#                    "Re-run simulation", icon("random")
#     ))
#   )
# }
