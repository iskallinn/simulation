run.simulation <- function (x) { # x = output from gen0
next.gen <- rbindlist(x[1])
next.gen.males <- rbindlist(x[2])
if (selection.method == blup) {
  big.pedfile <- rbindlist(x[4])
  pedfile <- rbindlist(x[3])
} else if (selection.method == selection) {
  pedfile <- rbindlist(x[3])
}
  # print( y )
temp <- copy(next.gen)
set( temp, j=which(colnames(temp) %in% c("obs_fert","perm.env"))  , value=NULL ) 
t <- rbind(next.gen.males,temp, fill=T) # needed down the road 
remove(temp) 
############### Make barren males ##########################
next.gen.males <- make.barren.males(next.gen.males)
next.gen <- mat.will(next.gen)
# ############### assign each female a male,  based on his mating willingness #####
mating.list <- mate(next.gen.males,next.gen)
# ############### Update pedigree ########################
pedfile <- update.pedigree(pedfile, next.gen, next.gen.males)
# ############### Inbreeding Coefficient and calculation of breeding value#############################
  # mating.list <- subset( mating.list,  sire.id>0 ) # remove the females who were not mated
mating.list <- dam.age (mating.list)  # checks the dam age and puts the effect for yearlings

mating.list = transform( mating.list,  obs_fert =  rpois(nrow(mating.list),
                                                         lambda = exp(1.95 + perm.env.ls + dam.fert+dam.age))
                         *barren*semen.quality )
set( mating.list, j=which(colnames(mating.list) %in% c("semen.quality","dam.age")) , value=NULL )
# 
# stat.crate[year+(runcounter -1)*(n+1),2] <- nrow(mating.list)
# 
mating.list <- subset(mating.list, obs_fert > 0)
if (make.obs.file == 1) {
  write.output(mating.list)
}
# 
if (selection.method == blup){
kit.list <- bv.n(mating.list,pedfile, big.pedfile)
} else if (selection.method == selection) {
  kit.list <- bv.n(mating.list,pedfile, pedfile)
}
if (selection.method== blup) {
  big.pedfile <- make.big.pedfile(kit.list, big.pedfile) # this makes the big pedigree with all animals in the pedigree
  create.obs.phenotypes (kit.list)
  solutions <- calculate.selection.index()
  solutions.bw <- blup.bw.nov()
}
# ############### Selection of next generation    #############
# # See utility functions for method
if (selection.method == selection) {
  old.females <- sel.old.females ( next.gen,mating.list)
  next.gen <- sel.yearlings.females (kit.list, old.females)
  next.gen.males <- sel.males (kit.list)
  ############### Index selection of next generation #############
} else if (selection.method == blup ) {
  big.pedfile <- update.big.pedigree(big.pedfile, next.gen, next.gen.males)
  old.females <- selind.old.females (next.gen, solutions, solutions.bw)
  next.gen <- indexsel.yearlings.females (kit.list,solutions,solutions.bw,old.females)
  next.gen.males <- indsel.males (kit.list,solutions,solutions.bw)
}
if("f0.dam" %in% colnames(old.females)) {
  set( old.females, j=which(colnames(old.females) %in%
                              "f0.dam")  , value=NULL )}
if("f0.dam" %in% colnames(next.gen.males)) {
  set( next.gen.males, j=which(colnames(next.gen.males) %in% "f0.dam")  , value=NULL )
}

next.gen <- rbind(next.gen,old.females)
# # gather up mean number of true sires
# stat.crate[year+(runcounter -1)*(n+1),1] <- c(mean(kit.list$true.sire == kit.list$sire.assumed))
# stat.crate[year+(runcounter -1)*(n+1),3] <- nrow(kit.list)
cat (y, mean(next.gen$fert),var(next.gen$fert),mean(mating.list$f0.dam)
     ,mean(mating.list$obs_fert), mean(next.gen$bs.phenotype), mean(next.gen$direct.genetic.body.size)
     ,mean(next.gen.males$bs.phenotype),var(next.gen$direct.genetic.body.size), sep="\t",file=con)
cat("\n",file=con)
if (selection.method == blup){
return (list(next.gen, next.gen.males,pedfile,big.pedfile))

} else if (selection.method == selection) {
  return (list(next.gen, next.gen.males,pedfile))  
}
}
run.simulation <-compiler::cmpfun(run.simulation,options= c(suppressAll=TRUE)) # performance boost
