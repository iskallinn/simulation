RunSimulation <-
  function (x,   # x = output from gen0
            year,
            p,
            selection.method,
            mblup,
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
            use.comb.ind.for.males,
            weighing.method,
            true.sire.chance,
            use.true.sire,
            weight.bw.old.females,
            weight.fert.old.females,
            weight.qual.old.females,
            weight.bw.kits,
            weight.fert.kits,
            weight.qual.kits,
            sorting.prop,
            pseudo.import,
            pseudo.import.prop,
            n,
            n.cages,
            variable.costs,
            fixed.costs,
            pelting.costs,
            price.sold.kit,
            cheat) 
{
    # browser()
  stat.crate <- c(0,0,0,0,0,0,0)
  next.gen <- rbindlist(x[1])
  number.of.females.start.of.year <- nrow(next.gen)
  
  next.gen.males <- rbindlist(x[2])
  if (selection.method == blup) {
    big.pedfile <- rbindlist(x[4])
    pedfile <- rbindlist(x[3])
    fert.memory <- unlist(x[5])
    n.females <- unlist(x[6])
    } else if (selection.method != blup) {
    pedfile <- rbindlist(x[3])
    fert.memory <- unlist(x[4])
    n.females <- unlist(x[5])
    }
  temp <- copy(next.gen)
  
  set(temp, j = which(colnames(temp) %in% c("obs_fert"))  , value =
        NULL)
  t <- rbind(next.gen.males, temp, fill = T) # needed down the road
  remove(temp)
  ############### Make barren males ##########################
  next.gen.males <-
    PrepareMalesForMating(
      next.gen.males,
      year,
      use.comb.ind.for.males,
      selection.method,
      weighing.method,
      intensity.remating,
      pseudo.import,
      pseudo.import.prop,
      weight.bw.kits,
      weight.fert.kits,
      weight.qual.kits
    )
  # browser()
  next.gen <-
    PrepareFemalesForMating(
      next.gen,
      year,
      mating.will.old.1st,
      mating.will.yearling.1st,
      mating.will.old.2nd,
      mating.will.yearling.2nd
    )
  # ############### assign each female a male,  based on his mating willingness #####
  mating.list <- mate(next.gen.males, 
                      next.gen, 
                      year,
                      use.blup.to.assort.mat,
                      selection.method,
                      crossmating,
                      purebreeding,
                      pr.barren.double.mating.yearling,
                      pr.barren.double.mating.old,
                      pr.barren.one.mating.yearling,
                      pr.barren.one.mating.old ) 
  # ############### Update pedigree ########################
  pedfile <- UpdatePedigree(pedfile, next.gen, next.gen.males, year)
  # ############### Inbreeding Coefficient and calculation of breeding value#############################
  mating.list <-
    YearlingEffectOnFertility (mating.list, year,yearling.effect)  # checks the dam age and puts the effect for yearlings
  # browser()
  if (crossmating == 1) { #arbitrarily have more kits if there is systematic 1+1 crossmating
    mating.list = transform( mating.list,  obs_fert =  rpois(nrow(mating.list),
                                                             lambda = exp(1.99 + perm.env.ls + dam.fert+dam.age))*barren*semen.quality)
  } else if (crossmating == 0 ) {
    mating.list = transform( mating.list,  obs_fert =  rpois(nrow(mating.list),
                                                             lambda = exp(1.95 + perm.env.ls + dam.fert+dam.age))*barren*semen.quality)
    
  }
  feed.usage.breeders <- FeedUsageBreeders(mating.list, next.gen.males, next.gen )
  set(mating.list,
      j = which(colnames(mating.list) %in% c("semen.quality", "dam.age")) ,
      value = NULL)
  stat.crate[1] <- nrow(mating.list) #number of mated females
  mating.list <- subset(mating.list, obs_fert > 0)
  stat.crate[2] <- (stat.crate[1] - nrow(mating.list))
  stat.crate[5] <- mean(mating.list$sire.id.1st == mating.list$sire.id.2nd) # percentage of females mated with 1st male
  stat.crate[6] <- ( mean(mating.list$sire.id.2nd == 0)) # percentage of single mated females
  # browser()
  
  if (selection.method == blup) {
    kit.list <-
      MakeKitsGenN(
        mating.list,
        pedfile,
        big.pedfile,
        year,
        p,
        leg2,
        t,
        true.sire.chance,
        use.true.sire,
        make.obs.file,
        qual.classes
      )
    stat.crate[4] <-nrow(kit.list)
    kit.list <- RandCull(kit.list,cull.ratio)
    stat.crate[7] <- nrow(kit.list)
    fert.memory[year] <- stat.crate[7]/n.females
    kit.list <- SellExtraKits(kit.list, stat.crate[1], n.cages)
    numb.sold.kits <- stat.crate[7]-nrow(kit.list)
    stat.crate[8] <- (nrow(kit.list)/2+n.females)/n.cages
    
    kit.list <- RFI(kit.list, leg2, leg1, t)
    kit.list.nomasked <- kit.list
    kit.list <- MaskKits(kit.list)
    # guess number of females and adjust male numbers
    n.females <- NumberofBreeders(fert.memory,n.cages,year)
    n.males <- ceiling( n.females/male.ratio )
    } else if (selection.method != blup) {
      kit.list <-
        MakeKitsGenN(
          mating.list,
          pedfile,
          pedfile,
          year,
          p,
          leg2,
          t,
          true.sire.chance,
          use.true.sire,
          make.obs.file,
          qual.classes
        )
      stat.crate[4] <-nrow(kit.list)
    kit.list <- RandCull(kit.list,cull.ratio)
    stat.crate[7] <- nrow(kit.list) #number of kits after culling
    fert.memory[year] <- stat.crate[7]/n.females
    kit.list <- SellExtraKits(kit.list, stat.crate[1], n.cages)
    numb.sold.kits <- stat.crate[7]-nrow(kit.list)
    stat.crate[8] <- (nrow(kit.list)/2+stat.crate[1])/n.cages
    kit.list <- RFI(kit.list, leg2, leg1, t)
    kit.list.nomasked <- kit.list
    kit.list <- MaskKits(kit.list)
    n.females <- NumberofBreeders(fert.memory,n.cages,year)
    n.males <- ceiling( n.females/male.ratio )
    }
  if (selection.method == blup) {
    big.pedfile <- WriteBigPedigree (kit.list, big.pedfile, year, p)
    stat.crate[3] <- mean(kit.list$true.sire == kit.list$sire.assumed)
    if(trace.ped == 1 ){TracePed(kit.list,next.gen)}
    if (mblup == 0) {
     WriteObservations(mating.list, next.gen,next.gen.males,kit.list,year,p,sorting.prop)
      solutions <- CalculateBLUP ()
    } else if (mblup ==1 ) { WriteMBLUPObservations(mating.list, next.gen, next.gen.males, kit.list, year,p,cheat)
     solutions <- CalculateMBLUP ()
      }
  }
  
  stat.crate[3] <- mean(kit.list$true.sire == kit.list$sire.assumed)
  kit.list$birthyear.dam <- NULL
  # ############### Selection of next generation    #############
  # # See utility functions for method
  
  if (selection.method == phenotypic) {
    old.females <-
      PhenoSelectionOldFemales (next.gen,
                                mating.list,
                                year,
                                max.age.females,
                                n.females,
                                prop.oldfemales)
    next.gen <-
      PhenoSelectionFemaleKits (kit.list,
                                old.females,
                                quantile.setting.ls,
                                quantile.setting.bw,
                                n.females)
    next.gen.males <-
      PhenoSelectionMaleKits (kit.list, quantile.setting.ls, quantile.setting.bw,n.males)
    ############### Index selection of next generation #############
  } else if (selection.method == random ) {
    old.females <- RandomSelectionOldFemales(next.gen,n.females,prop.oldfemales,mating.list,year)
    next.gen <- RandomSelectionYearlings(kit.list, old.females,n.females)
    next.gen.males <- RandomSelectionMales(kit.list,n.males) 
  }  else if (selection.method == blup) {
    big.pedfile    <-
      update.big.pedigree (big.pedfile, next.gen, next.gen.males)
    old.females    <-
      IndSelectionOldFemales (next.gen, solutions, year,
                              weight.bw.old.females,
                              weight.fert.old.females,
                              weight.qual.old.females,n.females,prop.oldfemales)
    next.gen       <-
      IndSelFemaleKits (
        kit.list,
        solutions,
        old.females,
        mblup,
        weight.bw.kits,
        weight.fert.kits,
        weight.qual.kits,
        n.females
      )
    next.gen.males <-
      IndSelMaleKits (kit.list,
                      solutions,
                      weight.bw.kits,
                      weight.fert.kits,
                      weight.qual.kits,
                      n.males)
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
  if ("obs_fert" %in% colnames(next.gen)) {
    set(next.gen,
        j = which(colnames(next.gen) %in% c("obs_fert","dam.age"))  ,
        value = NULL)
  }
  
  next.gen$FI <- NULL
  next.gen <- rbind(next.gen, old.females, fill = T)
  feed.intake <- sum(kit.list$FI)
  feed.intake.pr.kit <- feed.intake/nrow(kit.list)
  kit.list.masked <- kit.list
  kit.list <- SkinPrices(kit.list.nomasked, next.gen, next.gen.males,year)
  skin.price <- sum(kit.list$skin.price, na.rm =T)/number.of.females.start.of.year
  income <- sum(kit.list$skin.price, na.rm =T)
  # if (year == n) {
  #   save(kit.list, file="last.year.Rdata")
  # }
  if(mean(is.na(next.gen$id)) > 0) {
    stop("NA in id's")
  }
  con <- file(description = "results", open = "a")
  # browser()
  if (selection.method == blup) {
    kit.list <- merge(kit.list.masked, solutions, by= "id", all.x=TRUE) # merge to solutions of blup of fertility
    stat <- summaryBy(phenotype.bw.oct + phenotype.skin.length ~ sex, data = kit.list, FUN= c(mean))
    stat1 <- subset(kit.list, sex==1)#males
    lm1  <- lm(data=kit.list, litter.size~blup.fert )
    lm2  <- lm(data=kit.list, live.qual~blup.qual )
    lm3  <- lm(data=kit.list, bw_m~blup.bw.male )
    
    cat (
    year,
    mean(kit.list$litter.size),
    var(kit.list$litter.size),
    mean(kit.list$f0),
    mean(mating.list$obs_fert),
    stat[[2,2]],
    mean(kit.list$bw_m),
    stat[[1,2]],
    var(kit.list$bw_m,na.rm=T),
    cor(kit.list$bw_m, kit.list$blup.bw.male,use="complete"),
    cor(stat1$bw_m, stat1$phenotype.bw.oct), #only males
    mean(kit.list$skin.length.male),
    var(kit.list$skin.length.male),
    mean(kit.list$skin.qual),
    var(kit.list$skin.qual),
    cor(kit.list$blup.fert, kit.list$litter.size, use="complete"),
    cor(kit.list$own_littersize, kit.list$litter.size),
    stat.crate[1],
    stat.crate[2],
    stat.crate[3],
    stat.crate[4],
    stat.crate[5],
    stat.crate[6],
    stat.crate[7],
    mean(kit.list$live.qual),
    var(kit.list$live.qual),
    cor(kit.list$blup.qual, kit.list$live.qual, use="complete"),
    cor(kit.list$blup.qual, kit.list$live.score, use="complete"),
    cor(kit.list$blup.qual, kit.list$skin.qual, use="complete"),
    stat[[1,3]],
    stat[[2,3]],
    #sum(kit.list$FI)/(nrow(kit.list)-(n.females*(1-prop.oldfemales)+n.males)),
    feed.intake/(stat.crate[7]-(1-prop.oldfemales)*n.females-n.males),
    skin.price,
   income - number.of.females.start.of.year *
      variable.costs -
      feed.intake * feed.price - fixed.costs - ceiling(stat.crate[7]-n.females*(1-prop.oldfemales)-n.males-numb.sold.kits) * pelting.costs +
      numb.sold.kits * price.sold.kit, #pr farm margin
    income, #income from skins
    number.of.females.start.of.year *variable.costs, #variable costs
    fixed.costs, #fixed costs
    ceiling(stat.crate[7]-n.females*(1-prop.oldfemales)-n.males) * pelting.costs, #pelting costs
    (feed.intake+feed.usage.breeders)*feed.price, #feeding costs
    numb.sold.kits*price.sold.kit, #sold kits
    (feed.intake+feed.usage.breeders)*feed.price/nrow(kit.list.nomasked),   
    skin.price*n.females/nrow(kit.list.nomasked), 
    coef(lm1)[2],
    coef(lm2)[2],
    coef(lm3)[2],
    stat.crate[8],
    number.of.females.start.of.year,
    ceiling(stat.crate[7]-n.females*(1-prop.oldfemales)-n.males-numb.sold.kits),
    (number.of.females.start.of.year *
       variable.costs +
       feed.intake * feed.price + fixed.costs + nrow(kit.list) * pelting.costs) /ceiling(stat.crate[7]-n.females*(1-prop.oldfemales)-n.males),
    (number.of.females.start.of.year *
       variable.costs +
       feed.intake * feed.price + fixed.costs + nrow(kit.list) * pelting.costs)/number.of.females.start.of.year,
   feed.intake.pr.kit,
    sep = "\t",
    file = con
  )
  } else if (selection.method != blup) {
    # browser()
    stat <- summaryBy(phenotype.bw.oct + phenotype.skin.length ~ sex, data = kit.list, FUN= c(mean))
    stat1 <- subset(kit.list, sex==1)#males
    
        cat ( 
      year,
      mean(next.gen$litter.size),
      var(next.gen$litter.size),
      mean(mating.list$f0.dam, na.rm=TRUE),
      mean(mating.list$obs_fert),
      stat[[2,2]], # avg oct weight of females
      mean(kit.list.nomasked$bw_m),
      stat[[1,2]], # avg oct weight of males
      var(kit.list.nomasked$bw_m),
      cor(stat1$bw_m, stat1$phenotype.bw.oct),
      mean(kit.list.nomasked$skin.length.male),
      var(kit.list.nomasked$skin.length.male),
      mean(kit.list.nomasked$skin.qual),
      var(kit.list.nomasked$skin.qual),
      cor(kit.list.nomasked$own_littersize, kit.list.nomasked$litter.size),
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
      feed.intake/(stat.crate[7]-(1-prop.oldfemales)*n.females-n.males),
      sum(kit.list$skin.price, na.rm =T)/n.females,
      sum(kit.list$skin.price, na.rm = T) - number.of.females.start.of.year *
        variable.costs -
        feed.intake * feed.price - fixed.costs - ceiling(stat.crate[7]-n.females*(1-prop.oldfemales)-n.males-numb.sold.kits) * pelting.costs +
        numb.sold.kits * price.sold.kit, #pr farm margin
      sum(kit.list$skin.price, na.rm = T), #income from skins
      number.of.females.start.of.year *variable.costs, #variable costs
        fixed.costs, #fixed costs
      ceiling(stat.crate[7]-n.females*(1-prop.oldfemales)-n.males) * pelting.costs, #pelting costs
      (feed.intake+feed.usage.breeders)*feed.price, #feeding costs
      numb.sold.kits*price.sold.kit,
      (feed.intake+feed.usage.breeders)*feed.price/nrow(kit.list.nomasked),  
      mean(kit.list$skin.price),  
      stat.crate[8],  
      number.of.females.start.of.year,   
      ceiling(stat.crate[7]-n.females*(1-prop.oldfemales)-n.males-numb.sold.kits),
      (number.of.females.start.of.year *
        variable.costs +
        feed.intake * feed.price + fixed.costs + nrow(kit.list) * pelting.costs) /ceiling(stat.crate[7]-n.females*(1-prop.oldfemales)-n.males),
      (number.of.females.start.of.year *
         variable.costs +
         feed.intake * feed.price + fixed.costs + nrow(kit.list) * pelting.costs)/number.of.females.start.of.year,
      feed.intake.pr.kit,
      sep = "\t",  
      file = con
    )  
    } 
  cat("\n", file = con)
  # close(con = con)
  closeAllConnections()
  if (selection.method == blup) {
    return (list(next.gen, next.gen.males, pedfile, big.pedfile,fert.memory,n.females))
    
  } else if (selection.method != blup) {
    return (list(next.gen, next.gen.males, pedfile,fert.memory,n.females))
  }
}
RunSimulation <-
  compiler::cmpfun(RunSimulation, options = c(suppressAll = TRUE)) # performance boost
