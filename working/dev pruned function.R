############## Trace pedigree ###############
# this is a function to trim the pedigree which gets too big and unwieldy after
# a certain number of generations
# it takes a list of animals to trace and returns a trimmed pedigree that is 
# then used by DMU to calculate breeding values
# currently unused, controlled by switch since it did not improve functionality

TracePed <- function(list.trace, next.gen) {
 trace <- file(description = "trace", open="w")
 write.table(rbind(list.trace[,.(id)],next.gen[,.(id)]), 
             file= trace,
             col.names=FALSE, row.names=FALSE, quote=FALSE)
 close(trace)
system2("run_DmuTrace.bat", "trace")
}