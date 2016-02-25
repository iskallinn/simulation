generate.gen0 <- function ()  {
  #cat("ID","prodyear","damage","obs.fert", sep="\t", file=output) 
  #cat("\n", file= output) 
  
  
  if (make.obs.file == 1) {
    pedigree <- file(description = paste("pedigree_",p, sep=""), open="w")
  }
  year <- 1
  # runcounter <- sum(p)
  # if (selection.method == blup) {
   modify.dir.file ()
  # }
  ############### Create base population ############
  gen0.females <- generate.base.females() # create females
  effgen0.males <- generate.base.males() # create males
  ################## assign each female a male,  based on his mating willingness ####
   mating.list <- mate (effgen0.males,gen0.females)
  # ############### make gen 0 pedigree & and some bookkeeping ##############
  # # This is not really needed at this point since all gen0 animals are unrelated
  # 
  mating.list <- dam.age (mating.list)
  # 
  mating.list = transform( mating.list,  obs_fert =  rpois(nrow(mating.list),
                                                           lambda = exp(1.95 + perm.env.ls + dam.fert+dam.age))*barren*semen.quality)
  # # TODO how to implement inbreeding depression in this poisson model??
  # 
  set( mating.list, j=which(colnames(mating.list) %in% c("semen.quality","dam.age")) , value=NULL )
  # 
  # # if (runcounter == 1 ) {
  # #   # stat.crate[1,2] <- nrow(mating.list)
  # # } else if (runcounter > 1) { 
  #   # stat.crate[year+(runcounter -1)*(n+1),2] <- nrow(mating.list)}
  mating.list <-  subset( mating.list,  obs_fert >  0 ) # remove females who are barren or mated with barren male
  # # it is important to get rid of negative litter sizes! otherwise code crashes at many points
  # ############### Breeding value of offspring ###############################
  # # put in dam id's and make genetic value of each kit for fertility
  # # see utility functions for bv function
  gen1 <- bv(mating.list, effgen0.males)
  # 
  pedfile <- make.pedfile.gen0(gen0.females,effgen0.males)
  if (selection.method == blup) {
    big.pedfile <- make.big.pedfile(gen1, pedfile)
    create.obs.phenotypes (gen1)}
  # # At this point I think it is safe to delete some stuff from memory
  # # Note that I skip the construction of pedfile1 here. I don't think it is needed. Will check on that later
  # ############### Selection of first generation #########################################
  # 
  # # We choose females to fit n.females, based on prop of old females
  # # Note: currently all females are replaced in year 1
  # 
  # 
  old.females <- sel.old.females ( gen0.females,mating.list)
  next.gen <- sel.yearlings.females (gen1, old.females)
  next.gen.males <- sel.males (gen1)
  if("f0.dam" %in% colnames(old.females)) {
    set( old.females, j=which(colnames(old.females) %in%
                                "f0.dam")  , value=NULL )
  }
  next.gen <- rbind(next.gen, old.females)
  # # add in next gen and gen1 to big pedigree
  if (selection.method ==blup){
    big.pedfile <- update.big.pedigree (big.pedfile, next.gen, next.gen.males)
  }
  # ############## First year statistics #######################
  cat (0, mean(mating.list$dam.fert),var(mating.list$dam.fert),0,mean(mating.list$obs_fert), mean(next.gen$bs.phenotype)
       , mean(next.gen$direct.genetic.body.size),mean(next.gen.males$bs.phenotype), var(next.gen$direct.genetic.body.size), sep="\t",file=con)
  cat("\n",file=con)
  # 
  # # if (runcounter == 1) {
  # #   # stat.crate[[1,1]] <- (mean(gen1$true.sire == gen1$sire.assumed))
  # #   # stat.crate[1,3] <- nrow(gen1)
  # # } else if (runcounter > 1) {
  #   # stat.crate[[year+(runcounter -1)*(n+1),1]] <- (mean(gen1$true.sire == gen1$sire.assumed))
  #   # stat.crate[year+(runcounter -1)*(n+1),3] <-nrow(gen1)
  # # }
  if (make.obs.file == 1) {
    write.output(mating.list)
  }
  # 
  # remove(mating.list,gen0.females,gen1,effgen0.males,old.females) #remove all the stuff I don't need anymore
  # return( list(next.gen, next.gen.males,mating.list))
  if (selection.method == selection){ 
  return( list(next.gen,next.gen.males,pedfile))
  } else if (selection.method == blup) {
    return( list(next.gen,next.gen.males,pedfile,big.pedfile))
    }
}
generate.gen0 <-compiler::cmpfun(generate.gen0,options= c(suppressAll=TRUE)) # performance boost
