############## Trace pedigree ###############
# this is a function to trim the pedigree which gets too big and unwieldy after
# a certain number of generations
# it takes a list of animals to trace and returns a trimmed pedigree that is 
# then used by DMU to calculate breeding values

TracePed <- function(list.trace) {
 trace <- file(description = "trace", open="w")
  write.table(format(trace, scientific = FALSE),file= list.trace[,1], 
              quote = FALSE, col.names = FALSE, sep = "\t", row.names = FALSE  )
  close(trace)
system2("run_dmutrace", "trace")
  }