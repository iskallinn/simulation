############### Documentation ###################
# This function connects all the underlying functions and goes through the first
# generation and makes and selects offspring, creates a new pedigree
RunFirstYear <-
  function (p, # p = is the loopcounter for the replicates
            year,
            selection.method,
            mblup,
            trace.ped,
            crossmating,
            make.obs.file,
            leg2,
            t,
            n.females,
            mating.will.yearling.1st,
            mating.will.yearling.2nd,
            n.males,
            male.ratio,
            male.inf,
            intensity.remating,
            use.blup.to.assort.mat,
            purebreeding,
            pr.barren.double.mating.yearling,
            pr.barren.double.mating.old,
            pr.barren.one.mating.yearling,
            pr.barren.one.mating.old,
            yearling.effect,
            cull.ratio,
            leg1,
            max.age.females,
            prop.oldfemales,
            quantile.setting.ls,
            quantile.setting.bw,
            true.sire.chance,
            sorting.prop
  )
  {
    p <- p
    stat.crate <- c(0,0,0,0,0,0) # temporary holding space for values
  if (make.obs.file == 1) {
    pedigree <- file(description = paste("pedigree_",p, sep=""), open="w")
  }
  year <- 1
  ModifyDIRFile (p, mblup, trace.ped)
  ############### Create base population ############
  gen0.females <- GenerateBaseFemales( leg2,
                                       t,
                                       n.females,
                                       mating.will.yearling.1st,
                                       mating.will.yearling.2nd,
                                       qual.classes ) # create females
  effgen0.males <- GenerateBaseMales(leg2,
                                     t,
                                     n.males,
                                     male.ratio,
                                     male.inf,
                                     qual.classes,
                                     intensity.remating) # create males
  ################## assign each female a male,  based on his mating willingness ####
  mating.list <- mate (
    effgen0.males,
    gen0.females,
    year,
    use.blup.to.assort.mat,
    selection.method,
    crossmating,
    purebreeding,
    pr.barren.double.mating.yearling,
    pr.barren.double.mating.old,
    pr.barren.one.mating.yearling,
    pr.barren.one.mating.old
  )
  # ############### make gen 0 pedigree & and some bookkeeping ##############
  # # This is not really needed at this point since all gen0 animals are unrelated
  # 
  mating.list <- YearlingEffectOnFertility (mating.list, year,yearling.effect)
  if (crossmating == 1) {
  mating.list = transform( mating.list,  obs_fert =  rpois(nrow(mating.list),
                                                           lambda = exp(1.99 + perm.env.ls + dam.fert+dam.age))*barren*semen.quality)
  } else if (crossmating == 0 ) {
    mating.list = transform( mating.list,  obs_fert =  rpois(nrow(mating.list),
                                                             lambda = exp(1.95 + perm.env.ls + dam.fert+dam.age))*barren*semen.quality)
    
  }
  set( mating.list, j=which(colnames(mating.list) %in% c("semen.quality","dam.age")) , value=NULL )
  stat.crate[1] <- nrow(mating.list)
  mating.list <-  subset( mating.list,  obs_fert >  0 ) # remove females who are barren or mated with barren male
  stat.crate[2] <- (stat.crate[1] - nrow(mating.list))
  stat.crate[5] <- mean(mating.list$sire.id.1st == mating.list$sire.id.2nd) # percentage of females mated with 1st male
  stat.crate[6] <- ( mean(mating.list$sire.id.2nd == 0)) # percentage of single mated females
  # ############### Breeding value of offspring ###############################
  # # put in dam id's and make genetic value of each kit for fertility
  # # see utility functions for bv function
  kit.list <- MakeKitsGen0(mating.list, effgen0.males,year,leg2, qual.classes,true.sire.chance)
  # 
  stat.crate[3] <- mean(kit.list$true.sire == kit.list$sire.assumed)
  stat.crate[4] <-nrow(kit.list)
  kit.list <- RandCull(kit.list,cull.ratio)
  kit.list.for.stats <- kit.list
  stat.crate[7] <- nrow(kit.list)
  kit.list <- RFI(kit.list, leg2, leg1, t)
  kit.list.nomasked <- kit.list
  kit.list <- MaskKits(kit.list)
  
  pedfile <- MakePedfileGen0(gen0.females,effgen0.males)
  if (selection.method == blup) {
    big.pedfile <- WriteBigPedigree(kit.list, pedfile,year,p)
    WriteObservations(mating.list, gen0.females,effgen0.males,kit.list,year,p,sorting.prop)
    WriteMBLUPObservations(mating.list, gen0.females, effgen0.males, kit.list, year,p)
    }
  kit.list$birthyear.dam <- NULL # to  do figure out error
  # ############### Selection of first generation #########################################
  # # We choose females to fit n.females, based on prop of old females
  # # Note: currently all females are replaced in year 1
  old.females <- PhenoSelectionOldFemales ( gen0.females,mating.list,year,max.age.females,
                                            n.females,
                                            prop.oldfemales )
  next.gen <- PhenoSelectionFemaleKits (kit.list, old.females, quantile.setting.ls,
                                        quantile.setting.bw)
  next.gen.males <- PhenoSelectionMaleKits (kit.list,quantile.setting.ls,quantile.setting.bw)
  if("f0.dam" %in% colnames(old.females)) {
    set( old.females, j=which(colnames(old.females) %in%
                                "f0.dam")  , value=NULL )
  }
  next.gen <- rbind(next.gen, old.females,fill=TRUE)
  feed.intake <- sum(kit.list$FI)
  kit.list <- SkinPrices(kit.list.nomasked, next.gen, next.gen.males,year)
  # # add in next gen and kit.list to big pedigree
  if (selection.method ==blup){
    big.pedfile <- update.big.pedigree (big.pedfile, next.gen, next.gen.males)
  }
  # ############## First year statistics #######################
  con <- file(description="results",open="a")
  if (selection.method == blup) {
    stat <- summaryBy(phenotype.bw.oct + phenotype.skin.length ~ sex, data = kit.list, FUN= c(mean))
    stat1 <- subset(kit.list, sex==1)#males
    cat (
      year, #simulation year
      mean(kit.list.for.stats$litter.size), #avg genetic value of litter size
      var(kit.list.for.stats$litter.size), # variance of genetic value of litter size
      0, #avg inbreeding
      mean(mating.list$obs_fert), #observed fertility
      stat[[2,2]], #avg phenotype, october females
      mean(kit.list.for.stats$add.gen.bw.m), #avg gen val oct weight
      stat[[1,2]], #avg phenotype oct males
      var(kit.list.for.stats$add.gen.bw.m), #variance oct weight
      0, #correlation bw blup and phenotype
      cor(stat1$add.gen.bw.m, stat1$phenotype.bw.oct), #correlation bw phenotype and genetic value
      mean(kit.list.for.stats$skin.length.male), # avg genetic value for skin length
      var(kit.list.for.stats$skin.length.male),  # var of skin length
      mean(kit.list.for.stats$skin.qual), # avg genetic value of skin qual
      var(kit.list.for.stats$skin.qual), # var of skin qual
      0, # correlation of gen value litter size to blup
      0, # correlation of own litter size to genetic value 
      stat.crate[1], # number of mated females
      stat.crate[2], # number of barren females
      stat.crate[3], # proportion of kits with true sires
      stat.crate[4], # number of kits
      stat.crate[5], # percentages of females mated with "own" male
      stat.crate[6], # number of females mated once
      stat.crate[7], # survived kits
      mean(kit.list.for.stats$live.qual), # avg live quality
      var(kit.list.for.stats$live.qual),  #variance of live quality
      0, # correlation bw blup and gen value live qual
      0,
      0,
      stat[[1,3]], # skin length phenotype, male
      stat[[2,3]], # skin length phenotype, female
      #sum(kit.list$FI)/(nrow(kit.list)-(n.females*(1-prop.oldfemales)+n.males)),
      feed.intake/(stat.crate[7]-(1-prop.oldfemales)*n.females-n.males),
      sum(kit.list$skin.price, na.rm =T)/n.females,
      sum(kit.list$skin.price, na.rm =T)-(nrow(kit.list)*variable.costs)-feed.intake*feed.price, #pr farm margin
      feed.intake*feed.price/nrow(kit.list.nomasked),
      mean(kit.list$skin.price),
      sep = "\t",
      file = con
    )  } else if (selection.method == phenotypic) {
      stat <- summaryBy(phenotype.bw.oct + phenotype.skin.length ~ sex, data = kit.list, FUN= c(mean))
      stat1 <- subset(kit.list, sex==1)#males
    cat (
      year,
      mean(mating.list$dam.fert),
      var(mating.list$dam.fert),
      0,
      mean(mating.list$obs_fert),
      stat[[2,2]],
      mean(kit.list.for.stats$add.gen.bw.m),
      stat[[1,2]],
      var(kit.list.for.stats$add.gen.bw.m),
      cor(stat1$add.gen.bw.m, stat1$phenotype.bw.oct),
      mean(kit.list.for.stats$skin.length.male),
      var(kit.list.for.stats$skin.length.male),
      mean(kit.list.for.stats$skin.qual),
      var(kit.list.for.stats$skin.qual),
      0,
      stat.crate[1],
      stat.crate[2],
      stat.crate[3],
      stat.crate[4],
      stat.crate[5],
      stat.crate[6],
      stat.crate[7], # survived kits
      mean(kit.list.for.stats$live.qual), # avg live quality
      var(kit.list.for.stats$live.qual),  #variance of live quality
      stat[[1,3]], # skin length phenotype, male
      stat[[2,3]],
      #sum(kit.list$FI)/(nrow(kit.list)-(n.females*(1-prop.oldfemales)+n.males)),
      feed.intake/(stat.crate[7]-(1-prop.oldfemales)*n.females-n.males),
      sum(kit.list$skin.price, na.rm =T)/n.females,
      sum(kit.list$skin.price, na.rm =T)-(nrow(kit.list)*variable.costs)-feed.intake*feed.price, #pr farm margin
      feed.intake*feed.price/nrow(kit.list.nomasked),
      mean(kit.list$skin.price),
      sep = "\t",
      file = con
    )
    
    }
  cat("\n", file = con)
  close(con=con)
  if (make.obs.file == 1) {
    WriteFertObservations(mating.list,year,p)
  }
  if (selection.method == phenotypic){ 
  return( list(next.gen,next.gen.males,pedfile))
  } else if (selection.method == blup) {
    return( list(next.gen,next.gen.males,pedfile,big.pedfile))
    }
}
RunFirstYear <-compiler::cmpfun(RunFirstYear,options= c(suppressAll=TRUE)) # performance boost
