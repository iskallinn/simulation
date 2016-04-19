RunSimulation <- function (x, year, p) {
  # x = output from gen0
  stat.crate <- c(0,0,0,0,0,0)
  next.gen <- rbindlist(x[1])
  next.gen.males <- rbindlist(x[2])
  if (selection.method == blup) {
    big.pedfile <- rbindlist(x[4])
    pedfile <- rbindlist(x[3])
  } else if (selection.method == phenotypic) {
    pedfile <- rbindlist(x[3])
  }
  # print( y )
  temp <- copy(next.gen)
  
  set(temp, j = which(colnames(temp) %in% c("obs_fert", "perm.env"))  , value =
        NULL)
  t <- rbind(next.gen.males, temp, fill = T) # needed down the road
  remove(temp)
  ############### Make barren males ##########################
  next.gen.males <- PrepareMalesForMating(next.gen.males,year)
  next.gen <- PrepareFemalesForMating(next.gen, year)
  # ############### assign each female a male,  based on his mating willingness #####
  mating.list <- mate(next.gen.males, next.gen, year)
  # ############### Update pedigree ########################
  pedfile <- UpdatePedigree(pedfile, next.gen, next.gen.males, year)
  # ############### Inbreeding Coefficient and calculation of breeding value#############################
  # mating.list <- subset( mating.list,  sire.id>0 ) # remove the females who were not mated
  mating.list <-
    YearlingEffectOnFertility (mating.list, year)  # checks the dam age and puts the effect for yearlings
  
  mating.list = transform(mating.list,
                          obs_fert =  rpois(
                            nrow(mating.list),
                            lambda = exp(1.95 + perm.env.ls + dam.fert +
                                           dam.age)
                          )
                          * barren * semen.quality)
  set(mating.list,
      j = which(colnames(mating.list) %in% c("semen.quality", "dam.age")) ,
      value = NULL)
  #
  # stat.crate[year+(runcounter -1)*(n+1),2] <- nrow(mating.list)
  #
  stat.crate[1] <- nrow(mating.list)
  mating.list <- subset(mating.list, obs_fert > 0)
  stat.crate[2] <- (stat.crate[1] - nrow(mating.list))
  stat.crate[5] <- mean(mating.list$sire.id.1st == mating.list$sire.id.2nd) # percentage of females mated with 1st male
  stat.crate[6] <- ( mean(mating.list$sire.id.2nd == 0)) # percentage of single mated females
  if (make.obs.file == 1) {
    WriteFertObservations(mating.list, year, p)
  }
  #
  if (selection.method == blup) {
    kit.list <- MakeKitsGenN(mating.list, pedfile, big.pedfile, year, p)
    
  } else if (selection.method == phenotypic) {
    kit.list <- MakeKitsGenN(mating.list, pedfile, pedfile, year, p)
  }
  if (selection.method == blup) {
    big.pedfile <- WriteBigPedigree (kit.list, big.pedfile, year, p)
    # this makes the big pedigree with all animals in the pedigree
    # dirfile <- readLines("reml_bwnov.PAROUT")
    # if (weighing.method == oct ) {
    # dirfile[2] <- c(paste("  2  1  1    ",var(next.gen$bw.oct),sep=""))
    # } else if (weighing.method == sept) {
    #   dirfile[2] <- c(paste("  2  1  1    ",var(next.gen$bw.sept),sep=""))
    # }
    stat.crate[3] <- mean(kit.list$true.sire == kit.list$sire.assumed)
    stat.crate[4] <-nrow(kit.list)
    # writeLines(dirfile,"reml_bwnov.PAROUT")
    if(trace.ped == 1 ){TracePed(kit.list,next.gen)}
    solutions.littersize <- CalculateBLUPLitterSize ()
    WriteObservationFileBodyWeight (kit.list, year, p, solutions.littersize)
    solutions.bw.nov     <- CalculateBLUPBodyWeightNov ()
    solutions.qual <- CalculateBLUPQuality ()
  }
  stat.crate[3] <- mean(kit.list$true.sire == kit.list$sire.assumed)
  stat.crate[4] <-nrow(kit.list)
  kit.list$birthyear.dam <- NULL
  # ############### Selection of next generation    #############
  # # See utility functions for method
  if (selection.method == phenotypic) {
    old.females <-
      PhenoSelectionOldFemales (next.gen, mating.list, year)
    next.gen <- PhenoSelectionFemaleKits (kit.list, old.females)
    next.gen.males <- PhenoSelectionMaleKits (kit.list)
    ############### Index selection of next generation #############
  } else if (selection.method == blup) {
    big.pedfile    <-
      update.big.pedigree (big.pedfile, next.gen, next.gen.males)
    old.females    <-
      IndSelectionOldFemales (next.gen, solutions.littersize, solutions.bw.nov,solutions.qual, year)
    next.gen       <-
      IndSelFemaleKits (kit.list, solutions.littersize, solutions.bw.nov,solutions.qual ,old.females)
    next.gen.males <-
      IndSelMaleKits (kit.list, solutions.littersize, solutions.bw.nov,solutions.qual)
  }
  if ("f0.dam" %in% colnames(old.females)) {
    set(old.females,
        j = which(colnames(old.females) %in%
                    "f0.dam")  ,
        value = NULL)
  }
  if ("f0.dam" %in% colnames(next.gen.males)) {
    set(next.gen.males,
        j = which(colnames(next.gen.males) %in% "f0.dam")  ,
        value = NULL)
  }
  
  next.gen <- rbind(next.gen, old.females)
  # # gather up mean number of true sires
  # stat.crate[year+(runcounter -1)*(n+1),1] <- c(mean(kit.list$true.sire == kit.list$sire.assumed))
  # stat.crate[year+(runcounter -1)*(n+1),3] <- nrow(kit.list)
  con <- file(description = "results", open = "a")
  if (selection.method == blup) {
    kit.list <- merge(kit.list, solutions.littersize, by= "id", all.x=TRUE) # merge to solutions of blup of fertility
    kit.list <- merge(kit.list, solutions.bw.nov, by="id", all.x=TRUE)
    kit.list <- merge(kit.list, solutions.qual, by="id", all.x=TRUE)
    stat <- summaryBy(phenotype.bw.oct ~ sex, data = kit.list, FUN= c(mean))
    stat1 <- subset(kit.list, sex==1)#males
    
    cat (
    year,
    mean(kit.list$litter.size),
    var(kit.list$litter.size),
    mean(kit.list$f0),
    mean(mating.list$obs_fert),
    stat[[2,2]],
    mean(kit.list$bw.oct.male),
    stat[[1,2]],
    var(kit.list$bw.oct.male),
    cor(kit.list$bw.oct.male, kit.list$blup.bwnov),
    cor(stat1$bw.oct.male, stat1$phenotype.bw.oct), #only males
    mean(kit.list$skin.length.male),
    var(kit.list$skin.length.male),
    mean(kit.list$skin.qual),
    var(kit.list$skin.qual),
    cor(kit.list$blup.fert, kit.list$litter.size),
    cor(kit.list$own_littersize, kit.list$litter.size),
    stat.crate[1],
    stat.crate[2],
    stat.crate[3],
    stat.crate[4],
    stat.crate[5],
    stat.crate[6],
    mean(kit.list$live.qual),
    var(kit.list$live.qual),
    cor(kit.list$blup.qual, kit.list$live.qual),
    cor(kit.list$blup.qual, kit.list$live.score),
    cor(kit.list$blup.qual, kit.list$skin.qual),
    sep = "\t",
    file = con
  )
  } else if (selection.method == phenotypic) {
    stat <- summaryBy(phenotype.bw.oct ~ sex, data = kit.list, FUN= c(mean))
    stat1 <- subset(kit.list, sex==1)#males
    
        cat (
      year,
      mean(next.gen$litter.size),
      var(next.gen$litter.size),
      mean(mating.list$f0.dam),
      mean(mating.list$obs_fert),
      stat[[2,2]], # avg oct weight of females
      mean(kit.list$bw.oct.male),
      stat[[1,2]], # avg oct weight of males
      var(kit.list$bw.oct.male),
      cor(stat1$bw.oct.male, stat1$phenotype.bw.oct),
      mean(kit.list$skin.length.male),
      var(kit.list$skin.length.male),
      mean(kit.list$skin.qual),
      var(kit.list$skin.qual),
      cor(kit.list$own_littersize, kit.list$litter.size),
      stat.crate[1],
      stat.crate[2],
      stat.crate[3],
      stat.crate[4],
      stat.crate[5],
      stat.crate[6],
      mean(kit.list$live.qual),
      var(kit.list$live.qual),
      cor(kit.list$skin.qual, kit.list$blup.qual),
      cor(kit.list$live.score, kit.list$skin.qual),
      sep = "\t",
      file = con
    )
    }
  cat("\n", file = con)
  close(con = con)
  if (selection.method == blup) {
    return (list(next.gen, next.gen.males, pedfile, big.pedfile))
    
  } else if (selection.method == phenotypic) {
    return (list(next.gen, next.gen.males, pedfile))
  }
}
RunSimulation <-
  compiler::cmpfun(RunSimulation, options = c(suppressAll = TRUE)) # performance boost
