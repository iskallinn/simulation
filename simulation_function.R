Simulation <- function (
  weight.fert.old.females = 1,
  weight.bw.old.females   = 0,
  weight.qual.old.females = 0,
  weight.fert.kits        = 1,
  weight.bw.kits          = 0,
  weight.qual.kits        = 0,
  trace.ped = 0,       # if 1 then DmuTrace is used to prunce ped 5 gens back
  mask.phenotypes = 1, # 0 no masking, 1 mask phenotypes of kits from small litters
  sept = 1, # weigh kits in september CURRENTLY DEFUNCT, ONLY SEPT = 1 SUPPORTED
  oct = 0,  # weigh kits in october CURRENTLY DEFUNCT
  use.blup.to.assort.mat = 0, # if 1 then the males and females will be sorted on 
  use.comb.ind.for.males = 1, # if 1 then the usage of the males will depend on 
  mblup = 0, # If == 1 then MBLUP will be used for selection index, note that 
  n.females =  1000,             # NUMBER OF FEMALES
  nruns = 1,                     # how many replicates of the simulation
  n = 5     ,                    # number of generation per replicate
  selection.method = 2,       # selection strategy, 
  qual.classes = 5      ,        # quality classes, 5 or 10 are supported
  intensity.remating = 0.2,      # this controls how many males are chosen not to be remated
  make.obs.file = 1  ,           # 1 = make observation file, 0 otherwise
  use.true.sire = 0  ,           # 1 if true sire of kits is wanted for BV prediction, 0 otherwise
  feed.price    = 2.85,          # feed price per kg feed
  male.ratio      = 6  ,         # MALE TO FEMALE RATIO
  male.inf        = 0.98,        # % Proportion OF MALE BEING NOT BARREN
  prop.oldfemales = 0.4  ,       # Proportion of older females
  max.age.females = 2     ,      # define the max age of females
  quantile.setting.ls = 0.6,     # amount of kits to throw away because of too low littersize
  quantile.setting.bw = 0.1,     # prop of kits to disqualify due to weight
  mating.will.yearling.1st          = 0.95, # probability of yearling being mated
  mating.will.yearling.2nd          = 0.92, # probability of yearling being remated
  mating.will.old.1st               = 0.98, # probability of old female being mated
  mating.will.old.2nd               = 0.93, # probability of old female being remated
  pr.barren.one.mating.yearling     = 0.8, # probability of single mated yearling being barren
  pr.barren.double.mating.yearling  = 0.9, # probability of double mated yearling being barren
  pr.barren.one.mating.old          = 0.9, # probability of single mated old female being barren
  pr.barren.double.mating.old       = 0.95, # probability of double mated old female being barren
  n.males =  ceiling( n.females/male.ratio ), # calculates needed amount of males 
  cull.ratio                        = 0.8 # survival rate of kits, farmwise from 2nd cnt to pelting
 ) # closing paranthesis for definitions 
  { # opening curly brace for function 
  if (selection.method == phenotypic) {
    use.blup.to.assort.mat <- 0
  }
  #setwd("C:/Users/Notandi/Dropbox/Projects/simulation of mink farm/Output/DMU analysis/")
  setwd("C:/Users/au384062/Dropbox/Projects/simulation of mink farm/Output/DMU analysis/")
  WriteLogFile(n.females, n,nruns,cull.ratio)
  
  skin.metrics.males <- file(description = "skin_metrics_males", open ="w")
  skin.metrics.females <- file(description = "skin_metrics_females", open ="w")
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
  if (selection.method == blup) {
  con <- file(description = "results", open = "w")
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
    "pr.farm.margin",
    "feeding.cost",
    "avg.skin.price",
    sep = "\t",
    file = con
  )  } else if (selection.method == phenotypic) {
    con <- file(description = "results", open = "w")
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
      "pr.farm.margin",
      "feeding.cost",
      "avg.skin.price",
      sep = "\t",
      file = con
    )
  }
  cat("\n", file = con)
  close(con = con)
  
  for (p in 1:nruns) {
    year <- 1
    l <- RunFirstYear(p , year, selection.method,mblup,trace.ped)
    for (y in 1:n) {
      year <- 1 + y
      l <- RunSimulation(l, year, p,selection.method, mblup)
      
    }
  }
  templog <- readLines("log.log")
  templog[3] <- c(paste(  "Simulation ended",
                          format(Sys.time(), " %b %d %X"), sep="")) 
  writeLines(templog, "log.log")
  
  closeAllConnections()
}