############### Documentation ###################
# This function connects all the underlying functions and goes through the first
# generation and makes and selects offspring, creates a new pedigree
RunFirstYear <- function (p,year)  { # p = is the loopcounter for the replicates
  #cat("ID","prodyear","damage","obs.fert", sep="\t", file=output) 
  #cat("\n", file= output) 
  p <- p
  stat.crate <- c(0,0,0,0,0,0)
  
  
  if (make.obs.file == 1) {
    pedigree <- file(description = paste("pedigree_",p, sep=""), open="w")
  }
  year <- 1
  # runcounter <- sum(p)
  # if (selection.method == blup) {
  ModifyDIRFile (p)
  # }
  ############### Create base population ############
  gen0.females <- GenerateBaseFemales() # create females
  effgen0.males <- GenerateBaseMales() # create males
  ################## assign each female a male,  based on his mating willingness ####
   mating.list <- mate (effgen0.males,gen0.females,year)
  # ############### make gen 0 pedigree & and some bookkeeping ##############
  # # This is not really needed at this point since all gen0 animals are unrelated
  # 
  mating.list <- YearlingEffectOnFertility (mating.list, year)
  # 
  mating.list = transform( mating.list,  obs_fert =  rpois(nrow(mating.list),
                                                           lambda = exp(1.95 + perm.env.ls + dam.fert+dam.age))*barren*semen.quality)
  # # TODO how to implement inbreeding depression in this poisson model??
  # 
  set( mating.list, j=which(colnames(mating.list) %in% c("semen.quality","dam.age")) , value=NULL )
  # 
  # if (runcounter == 1 ) {
  #   # stat.crate[1,2] <- nrow(mating.list)
  # } else if (runcounter > 1) {
  # stat.crate[year+(runcounter -1)*(n+1),2] <- nrow(mating.list)}
  stat.crate[1] <- nrow(mating.list)
  mating.list <-  subset( mating.list,  obs_fert >  0 ) # remove females who are barren or mated with barren male
  stat.crate[2] <- (stat.crate[1] - nrow(mating.list))
  stat.crate[5] <- mean(mating.list$sire.id.1st == mating.list$sire.id.2nd) # percentage of females mated with 1st male
  stat.crate[6] <- ( mean(mating.list$sire.id.2nd == 0)) # percentage of single mated females
  # # it is important to get rid of negative litter sizes! otherwise code crashes at many points
  # ############### Breeding value of offspring ###############################
  # # put in dam id's and make genetic value of each kit for fertility
  # # see utility functions for bv function
  kit.list <- MakeKitsGen0(mating.list, effgen0.males,year)
  # 
  stat.crate[3] <- mean(kit.list$true.sire == kit.list$sire.assumed)
  stat.crate[4] <-nrow(kit.list)
  pedfile <- MakePedfileGen0(gen0.females,effgen0.males)
  if (selection.method == blup) {
    big.pedfile <- WriteBigPedigree(kit.list, pedfile,year,p)
    WriteObservationFileBodyWeight (kit.list, year,p)}
  kit.list$birthyear.dam <- NULL
  # # At this point I think it is safe to delete some stuff from memory
  # # Note that I skip the construction of pedfile1 here. I don't think it is needed. Will check on that later
  # ############### Selection of first generation #########################################
  # 
  # # We choose females to fit n.females, based on prop of old females
  # # Note: currently all females are replaced in year 1
  # 
  # 
  old.females <- PhenoSelectionOldFemales ( gen0.females,mating.list,year)
  next.gen <- PhenoSelectionFemaleKits (kit.list, old.females)
  next.gen.males <- PhenoSelectionMaleKits (kit.list)
  if("f0.dam" %in% colnames(old.females)) {
    set( old.females, j=which(colnames(old.females) %in%
                                "f0.dam")  , value=NULL )
  }
  next.gen <- rbind(next.gen, old.females,fill=TRUE)
  # # add in next gen and kit.list to big pedigree
  if (selection.method ==blup){
    big.pedfile <- update.big.pedigree (big.pedfile, next.gen, next.gen.males)
  }
  # ############## First year statistics #######################
  con <- file(description="results",open="a")
  if (selection.method == blup) {
    stat <- summaryBy(phenotype.bw.oct ~ sex, data = kit.list, FUN= c(mean))
    
    cat (
      year, #simulation year
      mean(kit.list$litter.size), #avg genetic value of litter size
      var(kit.list$litter.size), # variance of genetic value of litter size
      0, #avg inbreeding
      mean(mating.list$obs_fert), #observed fertility
      stat[[2,2]], #avg phenotype, october females
      mean(kit.list$bw.oct.male), #avg gen val oct weight
      stat[[1,2]], #avg phenotype oct males
      var(kit.list$bw.oct.male), #variance oct weight
      0, #correlation bw blup and phenotype
      cor(kit.list$bw.oct.male, kit.list$phenotype.bw.oct), #correlation bw phenotype and genetic value
      mean(kit.list$skin.length.male), # avg genetic value for skin length
      var(kit.list$skin.length.male),  # var of skin length
      mean(kit.list$skin.qual), # avg genetic value of skin qual
      var(kit.list$skin.qual), # var of skin qual
      0, # correlation of gen value litter size to blup
      0, # correlation of own litter size to genetic value 
      stat.crate[1], # number of mated females
      stat.crate[2], # number of barren females
      stat.crate[3], # proportion of kits with true sires
      stat.crate[4], # number of kits
      stat.crate[5], # percentages of females mated with "own" male
      stat.crate[6], # number of females mated once
      mean(kit.list$live.qual), # avg live quality
      var(kit.list$live.qual),  #variance of live quality
      0, # correlation bw blup and gen value live qual
      0,
      0,
      sep = "\t",
      file = con
    )  } else if (selection.method == phenotypic) {
    stat <- summaryBy(phenotype.bw.oct ~ sex, data = kit.list, FUN= c(mean))
    stat1 <- subset(kit.list, sex==1)#males
    cat (
      year,
      mean(mating.list$dam.fert),
      var(mating.list$dam.fert),
      0,
      mean(mating.list$obs_fert),
      stat[[2,2]],
      mean(kit.list$bw.oct.male),
      stat[[1,2]],
      var(kit.list$bw.oct.male),
      cor(stat1$bw.oct.male, stat1$phenotype.bw.oct),
      mean(kit.list$skin.length.male),
      var(kit.list$skin.length.male),
      mean(kit.list$skin.qual),
      var(kit.list$skin.qual),
      0,
      stat.crate[1],
      stat.crate[2],
      stat.crate[3],
      stat.crate[4],
      stat.crate[5],
      stat.crate[6],
      mean(kit.list$live.qual), # avg live quality
      var(kit.list$live.qual),  #variance of live quality
      sep = "\t",
      file = con
    )
    
    }
  cat("\n", file = con)
  close(con=con)
  # 
  # # if (runcounter == 1) {
  # #   # stat.crate[[1,1]] <- (mean(kit.list$true.sire == kit.list$sire.assumed))
  # #   # stat.crate[1,3] <- nrow(kit.list)
  # # } else if (runcounter > 1) {
  #   # stat.crate[[year+(runcounter -1)*(n+1),1]] <- (mean(kit.list$true.sire == kit.list$sire.assumed))
  #   # stat.crate[year+(runcounter -1)*(n+1),3] <-nrow(kit.list)
  # # }
  if (make.obs.file == 1) {
    WriteFertObservations(mating.list,year,p)
  }
  # 
  # remove(mating.list,gen0.females,kit.list,effgen0.males,old.females) #remove all the stuff I don't need anymore
  # return( list(next.gen, next.gen.males,mating.list))
  if (selection.method == phenotypic){ 
  return( list(next.gen,next.gen.males,pedfile))
  } else if (selection.method == blup) {
    return( list(next.gen,next.gen.males,pedfile,big.pedfile))
    }
}
RunFirstYear <-compiler::cmpfun(RunFirstYear,options= c(suppressAll=TRUE)) # performance boost
