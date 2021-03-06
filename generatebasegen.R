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
            sorting.prop,
            n,
            n.cages,
            variable.costs,
            fixed.costs,
            pelting.costs,
            price.sold.kit,
            cheat,
            genetic.means,
            risktaking,
            root,
            fileoutputpath,
            feed.price,
            maintenance.energy
  )
  {
    number.of.females.start.of.year <- n.females
    resultfile <- paste(root, fileoutputpath, 'raw_data/results', sep = '/')
    p <- p
    stat.crate <- c(0,0,0,0,0,0,0,0,0) # temporary holding space for values
    fert.memory <- rep(0, times= n)
  if (make.obs.file == 1) {
    pedigree <- file(description = paste("pedigree_",p, sep=""), open="w")
  }
  year <- 1
  # browser ()
  if(selection.method == 2) {
  ModifyDIRFile (p, mblup, trace.ped)}
  ############### Create base population ############
  gen0.females <- GenerateBaseFemales( leg2,
                                       t,
                                       n.females,
                                       mating.will.yearling.1st,
                                       mating.will.yearling.2nd,
                                       qual.classes,
                                       genetic.means) # create females
  effgen0.males <- GenerateBaseMales(leg2,
                                     t,
                                     n.males,
                                     male.ratio,
                                     male.inf,
                                     qual.classes,
                                     intensity.remating,
                                     n.females,
                                     genetic.means) # create males
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
  
  feed.usage.breeders <- FeedUsageBreeders(mating.list, effgen0.males, gen0.females,maintenance.energy )
  stat.crate[2] <- (stat.crate[1] - nrow(mating.list))
  stat.crate[5] <- mean(mating.list$sire.id.1st == mating.list$sire.id.2nd) # percentage of females mated with 1st male
  stat.crate[6] <- ( mean(mating.list$sire.id.2nd == 0)) # percentage of single mated females
  # ############### Breeding value of offspring ###############################
  # # put in dam id's and make genetic value of each kit for fertility
  # # see utility functions for bv function
  kit.list <- MakeKitsGen0(mating.list, effgen0.males,year,leg2, qual.classes,true.sire.chance)
  # 
  # browser()
  stat.crate[3] <- mean(kit.list$true.sire == kit.list$sire.assumed)
  stat.crate[4] <-nrow(kit.list)
  kit.list <- RandCull(kit.list,cull.ratio)
  kit.list.for.stats <- kit.list
  stat.crate[7] <- nrow(kit.list)
  fert.memory[year] <- stat.crate[7]/n.females
  kit.list <- SellExtraKits(kit.list, stat.crate[1], n.cages)
  numb.sold.kits <- stat.crate[7]-nrow(kit.list)
  
  stat.crate[8] <- (nrow(kit.list)/2+stat.crate[1])/n.cages
  kit.list <- RFI(kit.list, leg2, leg1, t)
  kit.list.nomasked <- kit.list
  kit.list <- MaskKits(kit.list)
  # guess number of females and adjust male numbers
  n.females <- NumberofBreeders(fert.memory,n.cages,year)
  n.males <- ceiling( n.females/male.ratio )
  pedfile <- MakePedfileGen0(gen0.females,effgen0.males)
  if (selection.method == blup) {
    big.pedfile <- WriteBigPedigree(kit.list, pedfile,year,p)
    WriteObservations(mating.list, gen0.females,effgen0.males,kit.list,year,p,sorting.prop)
    WriteMBLUPObservations(mating.list, gen0.females, effgen0.males, kit.list, year,p,cheat)
    }
  kit.list$birthyear.dam <- NULL # to  do figure out error
  # ############### Selection of first generation #########################################
  # # We choose females to fit n.females, based on prop of old females
  # # Note: currently all females are replaced in year 1
# browser()
  if (selection.method != 3 ) 
  {  
  old.females <- PhenoSelectionOldFemales ( gen0.females,mating.list,year,max.age.females,
                                            n.females,
                                            prop.oldfemales )
  next.gen <- PhenoSelectionFemaleKits (kit.list, old.females, quantile.setting.ls,
                                        quantile.setting.bw,n.females)
  next.gen.males <- PhenoSelectionMaleKits (kit.list,quantile.setting.ls,quantile.setting.bw,n.males)
} else if (selection.method == random) {
  kit.list <- EqualizeLitters(kit.list.nomasked)
  old.females <- RandomSelectionOldFemales(gen0.females,n.females,prop.oldfemales,mating.list,year)
    next.gen <- RandomSelectionYearlings(kit.list, old.females,n.females)
    next.gen.males <- RandomSelectionMales(kit.list,n.males) 
}
if("f0.dam" %in% colnames(old.females)) {
    set( old.females, j=which(colnames(old.females) %in%
                                "f0.dam")  , value=NULL )
  }
  next.gen <- rbind(next.gen, old.females,fill=TRUE)
  feed.intake.kits <- sum(kit.list.nomasked$FI)
  feed.intake.kits.pr.kit <- feed.intake.kits/nrow(kit.list.nomasked)
  # truncs <- StartPosSkins(kit.list,next.gen,next.gen.males) #defunct
  # htruncs <- StartPosSkinsVelvet(kit.list,next.gen,next.gen.males)
  kit.list <- SkinPrices(kit.list.nomasked, next.gen, next.gen.males,year,root,fileoutputpath,truncs,htruncs)
  income <- sum(kit.list$skin.price, na.rm =T)
  # # add in next gen and kit.list to big pedigree
  if (selection.method ==blup){
    big.pedfile <- update.big.pedigree (big.pedfile, next.gen, next.gen.males)
  }
  # browser()
  labcosts <-  LaborCosts(number.of.females.start.of.year, ceiling(stat.crate[7]-n.females*(1-prop.oldfemales)-n.males) )
  # ############## First year statistics #######################
  con <- file(description=resultfile,open="a")
  if (selection.method == blup) {
    stat <- summaryBy(phenotype.bw.oct + phenotype.skin.length ~ sex, data = kit.list, FUN= c(mean))
    stat1 <- subset(kit.list, sex==1)#males
    cat (
      year,                                      # simulation year
      mean(kit.list.for.stats$litter.size),      # avg genetic value of litter size
      var(kit.list.for.stats$litter.size),       # variance of genetic value of litter size
      0,                                         # avg inbreeding
      mean(mating.list$obs_fert),                # observed fertility
      stat[[2,2]],                               # avg phenotype, october females
      mean(kit.list.for.stats$bw_m),             # avg gen val oct weight
      stat[[1,2]],                               # avg phenotype oct males
      var(kit.list.for.stats$bw_m),              # variance oct weight
      0,                                         # correlation bw blup and phenotype
      cor(stat1$bw_m, stat1$phenotype.bw.oct),   #correlation bw phenotype and genetic value
      mean(kit.list.for.stats$skin.length.male), # avg genetic value for skin length
      var(kit.list.for.stats$skin.length.male),  # var of skin length
      mean(kit.list.for.stats$skin.qual),        # avg genetic value of skin qual
      var(kit.list.for.stats$skin.qual),         # var of skin qual
      0,                                         # correlation of gen value litter size to blup
      0,                                         # correlation of own litter size to genetic value 
      stat.crate[1],                             # number of mated females
      stat.crate[2],                             # number of barren females
      stat.crate[3],                             # proportion of kits with true sires
      stat.crate[4],                             # number of kits
      stat.crate[5],                             # percentages of females mated with "own" male
      stat.crate[6],                             # number of females mated once
      stat.crate[7],                             # survived kits
      mean(kit.list.for.stats$live.qual),        # avg live quality
      var(kit.list.for.stats$live.qual),         # variance of live quality
      0,                                         # correlation bw blup and gen value live qual
      0,
      0,
      stat[[1,3]],                               # skin length phenotype, male
      stat[[2,3]],                               # skin length phenotype, female
      #sum(kit.list$FI)/(nrow(kit.list)-(n.females*(1-prop.oldfemales)+n.males)),
      feed.intake.kits/(stat.crate[7]-(1-prop.oldfemales)*n.females-n.males),
      sum(kit.list$skin.price, na.rm =T)/n.females,
      income - number.of.females.start.of.year *
        variable.costs -
        (feed.intake.kits+feed.usage.breeders) * feed.price - fixed.costs - nrow(kit.list) * pelting.costs +
        numb.sold.kits * price.sold.kit,         # pr farm margin
      income,                                    # income from skins
      number.of.females.start.of.year *variable.costs+labcosts, #variable costs
      fixed.costs,                               # fixed costs
      ceiling(stat.crate[7] - n.females *
                (1 - prop.oldfemales) - n.males) * 
        pelting.costs,                           # pelting costs
      (feed.intake.kits+feed.usage.breeders)*feed.price, # feeding costs
      numb.sold.kits*price.sold.kit,             # income from sold kits
      (feed.intake.kits+feed.usage.breeders)*feed.price/nrow(kit.list.nomasked),
      mean(kit.list$skin.price),
      0,
      0,
      0,
      stat.crate[8],
      number.of.females.start.of.year,
      ceiling(stat.crate[7]-n.females*(1-prop.oldfemales)-n.males),
      (number.of.females.start.of.year *
          variable.costs + labcosts +
          feed.intake.kits * feed.price + fixed.costs +  ceiling(stat.crate[7] -
          n.females * (1 - prop.oldfemales) - n.males) * pelting.costs
      ) /ceiling(stat.crate[7]-n.females*(1-prop.oldfemales)-n.males),
      (number.of.females.start.of.year *
         variable.costs +labcosts +
         feed.intake.kits * feed.price + fixed.costs +  ceiling(stat.crate[7]-n.females*(1-prop.oldfemales)-n.males) * pelting.costs)/number.of.females.start.of.year,
      feed.intake.kits.pr.kit,
      labcosts,
      feed.usage.breeders,
      feed.intake.kits,
      sep = "\t",
      file = con
    )  
    } else if (selection.method != blup) {
      stat <- summaryBy(phenotype.bw.oct + phenotype.skin.length ~ sex, data = kit.list, FUN= c(mean))
      stat1 <- subset(kit.list, sex==1)#males
    cat (
      year,
      mean(mating.list$dam.fert),
      var(mating.list$dam.fert),
      0,
      mean(mating.list$obs_fert),
      stat[[2,2]],
      mean(kit.list.for.stats$bw_m),
      stat[[1,2]],
      var(kit.list.for.stats$bw_m),
      cor(stat1$bw_m, stat1$phenotype.bw.oct),
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
      stat.crate[7],
      mean(kit.list.nomasked$live.qual),
      var(kit.list.nomasked$live.qual),
      stat[[1,3]],
      stat[[2,3]],
      #sum(kit.list$FI)/(nrow(kit.list)-(n.females*(1-prop.oldfemales)+n.males)),
      feed.intake.kits/(stat.crate[7]-(1-prop.oldfemales)*n.females-n.males),
      sum(kit.list$skin.price, na.rm =T)/n.females,
      sum(kit.list$skin.price, na.rm = T) - number.of.females.start.of.year *
        variable.costs -
        (feed.intake.kits+feed.usage.breeders) * feed.price - fixed.costs - nrow(kit.list) * pelting.costs +
        numb.sold.kits * price.sold.kit-labcosts, #pr farm margin
      sum(kit.list$skin.price, na.rm = T), #income from skins
      number.of.females.start.of.year *variable.costs+labcosts, #variable costs
      fixed.costs, #fixed costs
      ceiling(stat.crate[7]-n.females*(1-prop.oldfemales)-n.males) * pelting.costs, #pelting costs
      (feed.intake.kits+feed.usage.breeders)*feed.price, #feeding costs
      numb.sold.kits*price.sold.kit,
      (feed.intake.kits+feed.usage.breeders)*feed.price/nrow(kit.list.nomasked),  
      mean(kit.list$skin.price),  
      stat.crate[8],  
      number.of.females.start.of.year,
      ceiling(stat.crate[7]-n.females*(1-prop.oldfemales)-n.males),
      (number.of.females.start.of.year *
         variable.costs +labcosts +
         feed.intake.kits * feed.price + fixed.costs +  ceiling(stat.crate[7]-n.females*(1-prop.oldfemales)-n.males) * pelting.costs) /ceiling(stat.crate[7]-n.females*(1-prop.oldfemales)-n.males),
      (number.of.females.start.of.year *
         variable.costs +labcosts +
         feed.intake.kits * feed.price + fixed.costs +  ceiling(stat.crate[7]-n.females*(1-prop.oldfemales)-n.males) * pelting.costs)/number.of.females.start.of.year,
      feed.intake.kits.pr.kit,
      labcosts,
      feed.usage.breeders,
      feed.intake.kits,
      sep = "\t",
      file = con
    )
    
    }
  cat("\n", file = con)
  close(con=con)
  if (make.obs.file == 1) {
    WriteFertObservations(mating.list,year,p)
  }
  if (selection.method != blup){ 
  return( list(next.gen,next.gen.males,pedfile,fert.memory,n.females))
  } else if (selection.method == blup) {
    return( list(next.gen,next.gen.males,pedfile,big.pedfile,fert.memory,n.females))
    }
}
RunFirstYear <-compiler::cmpfun(RunFirstYear,options= c(suppressAll=TRUE)) # performance boost
