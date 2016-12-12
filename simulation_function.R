Simulation <- function (  
  phenotypic = 1,  
  blup = 2,
  random = 3,
  cheat = 0,# to make bw obs for both males and females, only for estimating var comps
  sept = 1, # weigh kits in september CURRENTLY DEFUNCT, ONLY SEPT = 1 SUPPORTED
  oct = 0,  # weigh kits in october CURRENTLY DEFUNCT
  use.blup.to.assort.mat = 1, # if 1 then the males and females will be sorted on 
  #their combined index before mating
  use.comb.ind.for.males = 1, # if 1 then the usage of the males will depend on 
  # their combined index, not quality & weight
  risktaking = 0.4,
  mblup = 0, # If == 1 then MBLUP will be used for selection index, note that 
  n.females =  1000,             # NUMBER OF FEMALES
  nruns = 1,                     # how many replicates of the simulation
  n = 5,                         # number of generation per replicate
  assortative = 1,
  pseudo.import = 0,             # Delete pseudo imported male pedigree
  pseudo.import.prop = 0.2,      # proportion of best quality males to mask
  mating.method = assortative,   # mating method, random or assortative
  selection.method = phenotypic,       # selection strategy, 
  weighing.method = oct,         # control for when to weigh kits for selection cands
  crossmating = 0,               # 1 = systematic crossmating
  purebreeding = 1,
  yearling.effect = -0.07,       # Change at own risk, currently gives ~0.6 fewer kits on yearlings
  qual.classes = 5,              # quality classes, 5 or 10 are supported
  intensity.remating = 0.2,      # this controls how many males are chosen not to be remated
  make.obs.file = 1,             # 1 = make observation file, 0 otherwise
  use.true.sire = 0,             # 1 if true sire of kits is wanted for BV prediction, 0 otherwise
  feed.price    = 2.85,          # feed price per kg feed
  true.sire.chance = 0.85,       # probability of kits being sired by 2nd male
  trace.ped = 0,       # if 1 then DmuTrace is used to prunce ped 5 gens back
  mask.phenotypes = 1, # 0 no masking, 1 mask phenotypes of kits from small litters
  weight.fert.old.females = 0.175,
  weight.bw.old.females   = 0.425,
  weight.qual.old.females = 0.4,
  weight.fert.kits        = 0.175,
  weight.bw.kits          = 0.425,
  weight.qual.kits        = 0.40,
  prop.oldfemales = 0.4,         # Proportion of older females
  max.age.females = 2,           # define the max age of females
  quantile.setting.ls = 0.4,     # amount of kits to throw away because of too low littersize
  quantile.setting.bw = 0.4,     # prop of kits to disqualify due to weight
  male.ratio      = 6,           # MALE TO FEMALE RATIO
  male.inf        = 0.98,        # % Proportion OF MALE BEING NOT BARREN
  mating.will.yearling.1st          = 0.95, # probability of yearling being mated
  mating.will.yearling.2nd          = 0.92, # probability of yearling being remated
  mating.will.old.1st               = 0.98, # probability of old female being mated
  mating.will.old.2nd               = 0.93, # probability of old female being remated
  pr.barren.one.mating.yearling     = 0.8, # probability of single mated yearling being barren
  pr.barren.double.mating.yearling  = 0.9, # probability of double mated yearling being barren
  pr.barren.one.mating.old          = 0.9, # probability of single mated old female being barren
  pr.barren.double.mating.old       = 0.95, # probability of double mated old female being barren
  n.males =  ceiling( n.females/male.ratio ), # calculates needed amount of males 
  cull.ratio                        = 0.85, # survival rate of kits, farmwise from 2nd cnt to pelting 
  sorting.prop                      = 1, # proportion of animals to live grade
  n.cages = 5500,
  variable.costs = 165,          # variable costs per female
  pelting.costs =  12,           # pelting costs pr skin
  fixed.costs = 286*n.females,   # fixed costs at start of simulation 
  price.sold.kit = 80,
  genetic.means = c(
    0,                  # live.qual
    0,                  # h.length
    0,                  # skin.qual
    0,                  # skin.length.male
    0,                  # skin.length.female
    0,                  # litter.size
    0,                  # body weight females
    0,                  # body weight males
    0,                  # rfi1.m
    0,                  # rfi2.m
    0,                  # rfi1.f
    0                   # rfi2.f
  ),
  fileoutputpath = "simulation of mink farm/Output/DMU analysis/",
  root = 'C:/Users/au384062/Dropbox/Projects'
  # price per kit sold
) # closing paranthesis for definitions 
  { # opening curly brace for function 
  # browser()
  if (selection.method != blup) {
    use.blup.to.assort.mat <- 0
  }
  n.males =  ceiling( n.females/male.ratio ) # calculates needed amount of males 
  #setwd("C:/Users/Notandi/Dropbox/Projects/simulation of mink farm/Output/DMU analysis/")
  setwd("C:/Users/au384062/Dropbox/Projects/simulation of mink farm/Output/DMU analysis/")
  WriteLogFile( n.females,
                n,
                nruns,
                phenotypic,
                blup,
                mblup,
                selection.method,
                quantile.setting.ls,
                quantile.setting.bw,
                weight.fert.kits,
                weight.bw.kits,
                weight.qual.kits, 
                weight.fert.old.females,
                weight.bw.old.females,
                weight.qual.old.females,
                male.ratio,
                n.males,
                male.inf,
                prop.oldfemales,
                max.age.females,
                crossmating,
                purebreeding,
                cull.ratio,
                feed.price,
                variable.costs,
                use.true.sire,
                use.blup.to.assort.mat,
                trace.ped,
                intensity.remating,
                fileoutputpath,
                root
  )
  
  skin.metrics.males <- file(description = paste(root, fileoutputpath,"raw_data/skin_metrics_males", sep = '/'), open ="w")
  skin.metrics.females <- file(description = paste(root, fileoutputpath,"raw_data/skin_metrics_females", sep = '/'), open ="w")
  cat(
    "Gen",
    "S50", #50
    "S40",#40
    "S30",#30
    "S00",#00
    "S0",#0
    "S1",#1
    "S2",#2
    "S3",#3
    "S4",#4
    "S5",#5
    "purple",#purple
    "platinum",#platinum
    "burgundy",#burgundy
    "ivory",#ivory
    "vel3",#velv3
    "vel2",#velv2
    "vel1",#vel1
    "kl",#kl
    "long.nap",#long nap 
    "avg.price.males",
    "numb.animals",
    "\n",
  sep="\t",  
    file = skin.metrics.males
  ) 
  cat(
    "Gen",
    "S50", #50
    "S40",#40
    "S30",#30
    "S00",#00
    "S0",#0
    "S1",#1
    "S2",#2
    "S3",#3
    "S4",#4
    "S5",#5
    "purple",#purple
    "platinum",#platinum
    "burgundy",#burgundy
    "ivory",#ivory
    "vel3",#velv3
    "vel2",#velv2
    "vel1",#vel1
    "kl",#kl
    "long.nap",#long nap 
    "avg.price.females",
    "numb.animals",
    "\n",
    sep="\t",  
    file = skin.metrics.females
  ) 
  if (selection.method == 2) { # 2 = blup
  con <- file(description = paste(root, fileoutputpath,"raw_data/results", sep = '/'), open = "w")
  cat(
    "Gen",
    "Gmean",
    "Gvar",
    "Fis",
    "Obs.fert",
    "mean.phenotype.bw.females",
    "gen.value.bs",
    "mean.phenotype.bw.males",
    "bw.var",
    "cor.bw.to.blup",
    "cor.bw.phenotype",
    "skin.length.mean",
    "skin.length.var",
    "skin.qual.mean",
    "skin.qual.var",
    "cor.ls.blup",
    "cor.ls.own.to.ls",
    "mated.females",
    "barren.females",
    "numb.false.sires",
    "numb.kits",
    "remating.perc",
    "perc.single.mat",
    "survived.kits",
    "mean.gen.val.qual",
    "var.gen.val.qual",
    "cor.blup.qual.gen.val.qual",
    "cor.live.score.skin.qa",
    "cor.blup.qual.to.skin.qual",
    "avg.skin.length.male",
    "avg.skin.length.female", 
    "feed.per.skin",
    "skin.price",
    "margin",
    "income.from.skins",
    "variable.costs",
    "fixed.costs",
    "pelting.costs",
    "gross.feeding.costs",
    "income.fr.sold.kits",
    "feeding.cost.pr.skin",
    "avg.skin.price",
    "reg_EBV_TBV_LS",
    "reg_EBV_TBV_qual",
    "reg_EBV_TBV_size",
    "cage_utilization",
    "n.females",
    "sold.skins",
    "costs.pr.sold.skin",
    "costs.pr.female",
    "FI.pr.kit",
    "labor.costs",
    sep = "\t",
    file = con
  )  } else if (selection.method != blup) { # phenotypic or random
    
    con <- file(description = paste(root, fileoutputpath,"raw_data/results", sep = '/'), open = "w")
    cat(
      "Gen",
      "Gmean",
      "Gvar",
      "Fis",
      "Obs.fert",
      "mean.phenotype.bw.females",
      "gen.value.bs",
      "mean.phenotype.bw.males",
      "bw.var",
      "cor.bw.phenotype",
      "skin.length.mean",
      "skin.length.var",
      "skin.qual.mean",
      "skin.qual.var",
      "cor.ls.own.to.ls",
      "mated.females",
      "barren.females",
      "numb.false.sires",
      "numb.kits",
      "remating.perc",
      "perc.single.mat",
      "survived.kits",
      "mean.gen.val.qual",
      "var.gen.val.qual",
      "avg.skin.length.male",
      "avg.skin.length.female",
      "feed.per.skin",
      "skin.price",
      "margin",
      "income.from.skins",
      "variable.costs",
      "fixed.costs",
      "pelting.costs",
      "gross.feeding.costs",
      "income.fr.sold.kits",
      "feeding.cost.pr.skin",
      "avg.skin.price",
      "cage_utilization",
      "n.females",
      "sold.skins",
      "costs.pr.sold.skin",
      "costs.pr.female",
      "FI.pr.kit",
          "labor.costs",
      sep = "\t",
      file = con
    )
  }
  cat("\n", file = con)
  close(con = con)
  # browser()
  for (p in 1:nruns) {
    year <- 1
    l <- RunFirstYear(p,
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
                      fileoutputpath
    )
    
    for (y in 1:n) {
      year <- 1 + y
      l <-
        RunSimulation(
          l,
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
          cheat,
          risktaking,
          root,
          fileoutputpath
          )
      
    }
  }
  templog <- readLines(paste(root, fileoutputpath,"log.log", sep = '/'))
  templog[3] <- c(paste(  "Simulation ended",
                          format(Sys.time(), " %b %d %X"), sep="")) 
  writeLines(templog, paste(root, fileoutputpath,"log.log", sep = '/'))
  ReadAndSummarize (fileoutputpath,root)
  closeAllConnections()
}