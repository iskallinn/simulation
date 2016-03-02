# simulation
development of a simulation of a mink farm

The file definition.r contains all of the variables created so far

Utility function poisson.r contains all of the minor functions that are run with
in the simulatino

Generatebasegen.r is the first year of the simulation
runsimulation.r is the nth year of the simulation
simulation_function.r is the "overall" simulation function

Shiny app is in the "Shiny app" folder and contains the server.r and ui.r 
Currently very much an alpha. 

It will not work since a few folders are needed in the right spots, that has to do 
with DMU. To make it work, a the "output" folder is needed, which I can provide
Then the "setwd" needs to be changed at two spots to point to where the source 
files are located and secondly where the "output" folder is located wherein are the 
appropriate DMU .dir files and .bat to run DMU.
