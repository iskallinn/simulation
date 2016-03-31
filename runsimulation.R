RunSimulation <- function (x, year, p) {
  # x = output from gen0
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
  
  set(temp, j = which(colnames(temp) %in% c("obs_fert", "perm.env"))  , value =
        NULL)
  t <- rbind(next.gen.males, temp, fill = T) # needed down the road
  remove(temp)
  ############### Make barren males ##########################
  next.gen.males <- PrepareMalesForMating(next.gen.males)
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
  mating.list <- subset(mating.list, obs_fert > 0)
  if (make.obs.file == 1) {
    WriteFertObservations(mating.list, year, p)
  }
  #
  if (selection.method == blup) {
    kit.list <- MakeKitsGenN(mating.list, pedfile, big.pedfile, year, p)
  } else if (selection.method == selection) {
    kit.list <- bv.n(mating.list, pedfile, pedfile, year, p)
  }
  if (selection.method == blup) {
    big.pedfile <- WriteBigPedigree (kit.list, big.pedfile, year, p)
    # this makes the big pedigree with all animals in the pedigree
    dirfile <- readLines("reml_bwnov.PAROUT")
    dirfile[2] <- c(paste("  2  1  1    ",var(next.gen$direct.genetic.body.size),sep=""))
    writeLines(dirfile,"reml_bwnov.PAROUT")
    if(trace.ped == 1 ){TracePed(kit.list,next.gen)}
    solutions.littersize <- CalculateBLUPLitterSize ()
    WriteObservationFileBodyWeight (kit.list, year, p, solutions.littersize)
    solutions.bw.nov     <- CalulateBLUPBodyWeightNov ()
  }
  # ############### Selection of next generation    #############
  # # See utility functions for method
  if (selection.method == selection) {
    old.females <-
      PhenoSelectionOldFemales (next.gen, mating.list, year)
    next.gen <- PhenoSelectionFemaleKits (kit.list, old.females)
    next.gen.males <- PhenoSelectionMaleKits (kit.list)
    ############### Index selection of next generation #############
  } else if (selection.method == blup) {
    big.pedfile    <-
      update.big.pedigree (big.pedfile, next.gen, next.gen.males)
    old.females    <-
      IndSelectionOldFemales (next.gen, solutions.littersize, solutions.bw.nov, year)
    next.gen       <-
      IndSelFemaleKits (kit.list, solutions.littersize, solutions.bw.nov, old.females)
    next.gen.males <-
      IndSelMaleKits (kit.list, solutions.littersize, solutions.bw.nov)
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
  cat (
    year,
    mean(next.gen$fert),
    var(next.gen$fert),
    mean(mating.list$f0.dam)
    ,
    mean(mating.list$obs_fert),
    mean(next.gen$bs.phenotype),
    mean(next.gen$direct.genetic.body.size)
    ,
    mean(next.gen.males$bs.phenotype),
    var(next.gen$direct.genetic.body.size),
    cor(next.gen$direct.genetic.body.size, next.gen$blup.bwnov),
    cor(next.gen$direct.genetic.body.size, next.gen$bs.phenotype),
    sep = "\t",
    file = con
  )
  cat("\n", file = con)
  close(con = con)
  if (selection.method == blup) {
    return (list(next.gen, next.gen.males, pedfile, big.pedfile))
    
  } else if (selection.method == selection) {
    return (list(next.gen, next.gen.males, pedfile))
  }
}
RunSimulation <-
  compiler::cmpfun(RunSimulation, options = c(suppressAll = TRUE)) # performance boost
