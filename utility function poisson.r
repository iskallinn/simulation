############### Creation of Gen0 females ############
# This function creates the base population of females
GenerateBaseFemales <- function () {
  id        <-  seq(1:n.females)
  #   fert      <-  rnorm(n.females)*sqrt(variance.fertility)
  add.gen <- as.data.table(rmvnorm(n.females, sigma = sigma.new,method="svd" ))
  colnames(add.gen)     <-
    c(
      "live.qual",
      "h.length",
      "bw.oct.male",
      "bw.oct.female",
      "skin.qual",
      "skin.length.male",
      "skin.length.female",
      "litter.size",
      "bw.sept.male",
      "bw.sept.female",
      "bw.june.maternal",
      "bw.june"
    )
  perm.env.ls           <- rnorm(n.females) * sqrt(var.perm.env.ls)
  sex                   <-  rep(2, times = n.females)
  dam.id                <-  rep(0, times = n.females)
  sire.id               <-  rep(0, times = n.females)
  birthyear             <-  rep(0, times = n.females)
  mating.will.1st.round <-
    rbinom(n.females, 1, mating.will.yearling.1st)
  mating.will.2nd.round <- numeric(n.females)
  #create breeding value of fertility for females
  
  gen0.females <-  data.table(
    id,
    add.gen,
    perm.env.ls,
    sex,
    sire.id,
    dam.id,
    birthyear,
    mating.will.1st.round,
    mating.will.2nd.round
  )
  gen0.females$mating.will.2nd.round <-
    ifelse (gen0.females$mating.will.1st.round == 1,
            rbinom(sum(gen0.females$mating.will.1st.round),1,mating.will.yearling.2nd
      ), 0 )
  
  
  gen0.females[, c("litter.size")
               := lapply(.SD, function(x)
                 x * sqrt(variance.fertility)), .SDcols = c("litter.size")] # makes the variance proper
  
  gen0.females[, c("bw.oct.female")
               := lapply(.SD, function(x)
                 x * sqrt(var.bw.oct.female)), .SDcols = c("bw.oct.female")] # makes variance proper
  gen0.females[, c("bw.oct.male")
               := lapply(.SD, function(x)
                 x * sqrt(var.bw.oct.male)), .SDcols = c("bw.oct.male")] # makes variance proper
  
  gen0.females[, c("bw.sept.female")
               := lapply(.SD, function(x)
                 x * sqrt(var.bw.sept.female)), .SDcols = c("bw.sept.female")] # makes variance proper
  gen0.females[, c("bw.sept.male")
               := lapply(.SD, function(x)
                 x * sqrt(var.bw.sept.male)), .SDcols = c("bw.sept.male")] # makes variance proper
  
  gen0.females[, c("live.qual")
               := lapply(.SD, function(x)
                 x * sqrt(var.live.qual.gen)), .SDcols = c("live.qual")] # makes variance proper
  gen0.females[, c("skin.qual")
               := lapply(.SD, function(x)
                 x * sqrt(var.skin.qual.gen)), .SDcols = c("skin.qual")] # makes variance proper
  gen0.females[, c("skin.length.male")
               := lapply(.SD, function(x)
                 x * sqrt(var.skin.length.male)), .SDcols = c("skin.length.male")] # makes variance proper
  gen0.females[, c("skin.length.female")
               := lapply(.SD, function(x)
                 x * sqrt(var.skin.length.female)), .SDcols = c("skin.length.female")] # makes variance proper
  
  gen0.females[,`:=`(phenotype.bw.oct = 
                      mean.body.size.female.oct + bw.oct.female +
                      rnorm(nrow(gen0.females))*sqrt(var.c.bw.oct.female)+
                       rnorm(nrow(gen0.females))*sqrt(var.res.bw.oct.female),
                    phenotype.bw.sept = 
                      mean.body.size.female.sept + bw.sept.female +
                      rnorm(nrow(gen0.females))*sqrt(var.c.bw.sept.female)+
                      rnorm(nrow(gen0.females))*sqrt(var.res.bw.sept.female),
                    phenotype.live.qual = live.qual + 
                      rnorm(nrow(gen0.females))*sqrt(var.live.qual.res),
                    phenotype.skin.length = mean.skin.length.female + skin.length.female + 
                      rnorm(nrow(gen0.females))*sqrt(var.skin.length.c.female)+
                    rnorm(nrow(gen0.females))*sqrt(var.skin.length.res.female),
                    phenotype.skin.qual = skin.qual + rnorm(nrow(gen0.females))*
                      sqrt(var.skin.qual.res)
  )] 
  if (qual.classes == 5) {
    truncs <- qnorm(
      p = c(0.05, 0.3, 0.7, 0.95),
      mean = mean(gen0.females$phenotype.live.qual),
                  sd = sqrt(var(
                    gen0.females$phenotype.live.qual
                  )),
      lower.tail = TRUE
    )
    gen0.females[, `:=`(live.score = ifelse(
      phenotype.live.qual >= truncs[4],
      5,
      ifelse(
        truncs[3] < phenotype.live.qual & phenotype.live.qual <= truncs[4],
        4,
        ifelse(
          phenotype.live.qual > truncs[2] & phenotype.live.qual <= truncs[3],
          3,
          ifelse(
            phenotype.live.qual > truncs[1] & phenotype.live.qual <= truncs[2],
            2,
            ifelse(phenotype.live.qual <=
                     truncs[1], 1, 0)
          )
        )
      )
    ))]
  } else if (qual.classes == 10) {
    truncs <-
      qnorm(
        p = c(0.01, 0.05, 0.15, 0.3, 0.5, 0.7, 0.85, 0.95, 0.99),
        mean = mean(gen0.females$phenotype.live.qual),
                    sd = sqrt(var(
                      gen0.females$phenotype.live.qual
                    )),
        lower.tail = TRUE
      )
    gen0.females[, `:=`(live.score =
                         ifelse(
                           phenotype.live.qual >= truncs[9],
                           10,
                           ifelse(
                             truncs[8] < phenotype.live.qual & phenotype.live.qual <= truncs[9],
                             9,
                             ifelse(
                               phenotype.live.qual > truncs[7] & phenotype.live.qual <= truncs[8],
                               8,
                               ifelse(
                                 phenotype.live.qual > truncs[6] & phenotype.live.qual <= truncs[7],
                                 7,
                                 ifelse(
                                   phenotype.live.qual > truncs[5] & phenotype.live.qual <= truncs[6],
                                   6,
                                   ifelse(
                                     phenotype.live.qual > truncs[4] & phenotype.live.qual <= truncs[5],
                                     5,
                                     ifelse(
                                       phenotype.live.qual > truncs[3] & phenotype.live.qual <= truncs[4],
                                       4,
                                       ifelse(
                                         phenotype.live.qual > truncs[2] & phenotype.live.qual <= truncs[3],
                                         3,
                                         ifelse(
                                           phenotype.live.qual > truncs[1] & phenotype.live.qual <= truncs[2],
                                           2,
                                           ifelse(phenotype.live.qual <=
                                                    truncs[1], 1, 0)
                                         )
                                       )
                                     )
                                   )
                                 )
                               )
                             )
                           )
                         ))]
  }
  
  return (gen0.females)
}
############### Creation of Gen0 males#######################################
# This function creates  the base population of males
GenerateBaseMales <- function () {
  mating.willingness.1st <-  numeric( n.males )  
  mating.willingness.2nd  <-  numeric( n.males)
  semen.quality.1st      <-  numeric( n.males )  
  semen.quality.2nd      <-  numeric( n.males )  
  id                 <-  numeric( n.males )  
  add.gen <- as.data.table(rmvnorm(n.males, sigma = sigma.new,method="svd" ))
  colnames(add.gen)     <-
    c(
      "live.qual",
      "h.length",
      "bw.oct.male",
      "bw.oct.female",
      "skin.qual",
      "skin.length.male",
      "skin.length.female",
      "litter.size",
      "bw.sept.male",
      "bw.sept.female",
      "bw.june.maternal",
      "bw.june"
    )
  sex                <-  numeric( n.males )  
  dam.id             <-  numeric( n.males )  
  sire.id            <-  numeric( n.males ) 
  birthyear          <-  numeric( n.males ) 
  can.remate         <-  rep(0,times=n.males) 
  
  
  #TODO make this into something more meaningful once I have quality or something to rank the males on
  
  for ( i in 1:n.males )  {
    mating.willingness.1st[i]     <-  rZIP( 1,  mu = male.ratio,  sigma = 0.05 )               
    # this is guesswork
    mating.willingness.2nd [i]     <-  rZIP( 1,  mu = 8,  sigma = 0.05 )                       
    # too high, since this will make sure all females are mated 2
    id[i]       <-  n.females + i                                             
    # create ID for males
    semen.quality.1st[i] <-  rbinom( 1,  1,  male.inf )                           
    # create on off for barren males
    semen.quality.2nd [i] <- semen.quality.1st [i]
    sex [i]         <-  1                                                         
    # 1 = male, 2 = female
    dam.id [i]     <-  0                                                         
    # this is gen0 so unknown parents
    sire.id [i]   <-  0                                                         
    # this is gen0 so unknown parents
  }
  
  # make data table out of the males 
  gen0.males <-  data.table( id, mating.willingness.1st,mating.willingness.2nd, add.gen
                             , semen.quality.1st,semen.quality.2nd, sex, sire.id, dam.id, birthyear,can.remate ) 
  
  gen0.males[, c("litter.size")
               := lapply(.SD, function(x)
                 x * sqrt(variance.fertility)), .SDcols = c("litter.size")] # makes the variance proper
  
  gen0.males[, c("bw.oct.female")
               := lapply(.SD, function(x)
                 x * sqrt(var.bw.oct.female)), .SDcols = c("bw.oct.male")] # makes variance proper
  gen0.males[, c("bw.oct.male")
               := lapply(.SD, function(x)
                 x * sqrt(var.bw.oct.male)), .SDcols = c("bw.oct.male")] # makes variance proper
  
  gen0.males[, c("bw.sept.female")
               := lapply(.SD, function(x)
                 x * sqrt(var.bw.sept.female)), .SDcols = c("bw.sept.male")] # makes variance proper
  gen0.males[, c("bw.sept.male")
               := lapply(.SD, function(x)
                 x * sqrt(var.bw.sept.male)), .SDcols = c("bw.sept.male")] # makes variance proper
  
  gen0.males[, c("live.qual")
               := lapply(.SD, function(x)
                 x * sqrt(var.live.qual.gen)), .SDcols = c("live.qual")] # makes variance proper
  gen0.males[, c("skin.qual")
               := lapply(.SD, function(x)
                 x * sqrt(var.skin.qual.gen)), .SDcols = c("skin.qual")] # makes variance proper
  gen0.males[, c("skin.length.male")
               := lapply(.SD, function(x)
                 x * sqrt(var.skin.length.male)), .SDcols = c("skin.length.male")] # makes variance proper
  gen0.males[, c("skin.length.male")
               := lapply(.SD, function(x)
                 x * sqrt(var.skin.length.male)), .SDcols = c("skin.length.male")] # makes variance proper
  
  gen0.males[,`:=`(phenotype.bw.oct = 
                       mean.body.size.male.oct + bw.oct.male +
                       rnorm(nrow(gen0.males))*sqrt(var.c.bw.oct.male)+
                     rnorm(nrow(gen0.males))*sqrt(var.res.bw.oct.male),
                     phenotype.bw.sept = 
                       mean.body.size.male.sept + bw.sept.male +
                       rnorm(nrow(gen0.males))*sqrt(var.c.bw.oct.male)+
                     rnorm(nrow(gen0.males))*sqrt(var.res.bw.sept.male),
                     phenotype.live.qual = live.qual + 
                       rnorm(nrow(gen0.males))*sqrt(var.live.qual.res),
                   phenotype.skin.length = mean.skin.length.male + skin.length.male + 
                     rnorm(nrow(gen0.males))*sqrt(var.skin.length.c.male)+
                     rnorm(nrow(gen0.males))*sqrt(var.skin.length.res.male),
                   phenotype.skin.qual = skin.qual + rnorm(nrow(gen0.males))*
                     sqrt(var.skin.qual.res)
  )]
  
  if (qual.classes == 5) {
    truncs <- qnorm(
      p = c(0.05, 0.3, 0.7, 0.95),
      mean = mean(gen0.males$phenotype.live.qual),
                  sd = sqrt(var(
                    gen0.males$phenotype.live.qual
                  )),
      lower.tail = TRUE
    )
    gen0.males[, `:=`(live.score = ifelse(
      phenotype.live.qual >= truncs[4],
      5,
      ifelse(
        truncs[3] < phenotype.live.qual & phenotype.live.qual <= truncs[4],
        4,
        ifelse(
          phenotype.live.qual > truncs[2] & phenotype.live.qual <= truncs[3],
          3,
          ifelse(
            phenotype.live.qual > truncs[1] & phenotype.live.qual <= truncs[2],
            2,
            ifelse(phenotype.live.qual <=
                     truncs[1], 1, 0)
          )
        )
      )
    ))]
  } else if (qual.classes == 10) {
    truncs <-
      qnorm(
        p = c(0.01, 0.05, 0.15, 0.3, 0.5, 0.7, 0.85, 0.95, 0.99),
        mean = mean(gen0.males$phenotype.live.qual),
                    sd = sqrt(var(
                      gen0.males$phenotype.live.qual
                    )),
        lower.tail = TRUE
      )
    gen0.males[, `:=`(live.score =
                          ifelse(
                            phenotype.live.qual >= truncs[9],
                            10,
                            ifelse(
                              truncs[8] < phenotype.live.qual & phenotype.live.qual <= truncs[9],
                              9,
                              ifelse(
                                phenotype.live.qual > truncs[7] & phenotype.live.qual <= truncs[8],
                                8,
                                ifelse(
                                  phenotype.live.qual > truncs[6] & phenotype.live.qual <= truncs[7],
                                  7,
                                  ifelse(
                                    phenotype.live.qual > truncs[5] & phenotype.live.qual <= truncs[6],
                                    6,
                                    ifelse(
                                      phenotype.live.qual > truncs[4] & phenotype.live.qual <= truncs[5],
                                      5,
                                      ifelse(
                                        phenotype.live.qual > truncs[3] & phenotype.live.qual <= truncs[4],
                                        4,
                                        ifelse(
                                          phenotype.live.qual > truncs[2] & phenotype.live.qual <= truncs[3],
                                          3,
                                          ifelse(
                                            phenotype.live.qual > truncs[1] & phenotype.live.qual <= truncs[2],
                                            2,
                                            ifelse(phenotype.live.qual <=
                                                     truncs[1], 1, 0)
                                          )
                                        )
                                      )
                                    )
                                  )
                                )
                              )
                            )
                          ))]
  }
  setorder(gen0.males, -live.score, -phenotype.bw.oct)
  for (i in 1:ceiling((1-intensity.remating)*nrow(gen0.males))) {
    gen0.males$can.remate[i] <- 1
  }
  
  #make a subset of the males which will mate, this must be moved into the mating function 
  effgen0.males <- subset( gen0.males,  mating.willingness.1st > 0 ) 
  return (effgen0.males)
}
############### Mating list and mate function #################
# This functions performs the mating in both the first and subsequent generations
# It works on a simple principle, it first expands the male list by their mating 
# mating willingness and allocates each female a mate. Based on if the male 
# is allowed to remate, the male is then remated first with "their" females
# otherwise the males with "left over" matings are iterated through the list 
# of females using an inefficient looping function, it would speed up the program
# if a more lean method could be thought out to perform it. 

mate <- function (x, y, year) {
  # x = males, y = females
  if (year >2 & use.blup.to.assort.mat == 1 & selection.method ==blup) {
    setorder(x, -comb.ind)
    setorder(y, -comb.ind)
  }
  if (year <= 2) {
    x$comb.ind <- 0
  } else if( selection.method == phenotypic) {
    x$comb.ind <- 0
    
  }
  mating.list <-
    x[rep(seq(nrow(x)), mating.willingness.1st),  #expands the male list into a mating list, based on mat.will.1st
      c(
        "id",
        "litter.size",
        "bw.oct.male",
        "bw.oct.female",
        "bw.sept.male",
        "bw.sept.female",
        "live.qual",
        "skin.qual",
        "skin.length.male",
        "skin.length.female",
        "h.length",
        "bw.june.maternal",
        "bw.june",
        "semen.quality.1st",
        "semen.quality.2nd",
        "can.remate",
        "mating.willingness.1st",
        "mating.willingness.2nd",
        "live.score",
        "comb.ind"
      )
      , with = F] #specify which columns to incl.
  if (nrow(mating.list) > sum(y$mating.will.1st.round)) {
    mating.list <- mating.list[1:sum(y$mating.will.1st.round), ]
  }
  
  setnames(
    mating.list,
    c(
      "id",
      "bw.oct.male",
      "bw.oct.female",
      "litter.size",
      "bw.sept.male",
      "bw.sept.female",
      "live.qual",
      "skin.qual",
      "skin.length.male",
      "skin.length.female",
      "bw.june",
      "bw.june.maternal",
      "h.length",
      "live.score"
    ),
    c(
      "sire.id.1st",
      "sire.bw.oct.male.1st",
      "sire.bw.oct.female.1st",
      "sire.fert.1st",
      "sire.bw.sept.male.1st",
      "sire.bw.sept.female.1st",
      "sire.live.qual.1st",
      "sire.skin.qual.1st",
      "sire.skin.length.male.1st",
      "sire.skin.length.female.1st",
      "sire.bw.june.1st",
      "sire.bw.june.maternal.1st",
      "sire.h.length.1st",
      "sire.live.score.1st"
    )
  )
  setkey(mating.list, sire.id.1st)
  mating.list[mating.list, c(
    "dam.id",
    "dam.fert",
    "f0",
    "obs_fert",
    "dam.bw.oct.male",
    "dam.bw.oct.female",
    "dam.bw.sept.male",
    "dam.bw.sept.female",
    "dam.skin.qual",
    "dam.skin.length.male",
    "dam.skin.length.female",
    "dam.live.qual",
    "perm.env.ls",
    "birthyear.dam",
    "sire.id.2nd",
    "sire.bw.oct.male.2nd",
    "sire.bw.oct.female.2nd",
    "sire.bw.sept.male.2nd",
    "sire.bw.sept.female.2nd",
    "sire.fert.2nd",
    "sire.skin.qual.2nd",
    "sire.skin.length.male.2nd",
    "sire.skin.length.female.2nd",
    "sire.live.qual.2nd",
    "sire.h.length.2nd",
    "sire.bw.june.2nd",
    "sire.bw.june.maternal.2nd",
    "sire.live.score.2nd"
  ) := 0]
  # moved scaling of litter specific environment to the phenotype function later on
  mating.list[, `:=`(specific.env.bw = rnorm(nrow(mating.list)) 
  )]
  # here I subset the dam list to throw away those who will not mate on first round
  y <- subset (y, mating.will.1st.round == 1)
  mating.list$dam.id       <- y[1:nrow(mating.list), .(id)]
  mating.list$dam.fert     <- y[1:nrow(mating.list), .(litter.size)]
  mating.list$dam.bw.oct.male   <- y[1:nrow(mating.list), .(bw.oct.male)]
  mating.list$dam.bw.oct.female   <- y[1:nrow(mating.list), .(bw.oct.female)]
  mating.list$perm.env.ls  <- y[1:nrow(mating.list), .(perm.env.ls)]
  mating.list$dam.bw.sept.male   <- y[1:nrow(mating.list), .(bw.sept.male)]
  mating.list$dam.bw.sept.female   <- y[1:nrow(mating.list), .(bw.sept.female)]
  mating.list$dam.skin.length.male  <- y[1:nrow(mating.list), .(skin.length.male)]
  mating.list$dam.skin.length.female  <- y[1:nrow(mating.list), .(skin.length.female)]
  mating.list$dam.skin.qual  <- y[1:nrow(mating.list), .(skin.qual)]
  mating.list$dam.live.qual  <- y[1:nrow(mating.list), .(live.qual)]
  mating.list$dam.live.score <- y[1:nrow(mating.list), .(live.score)]
  mating.list$dam.h.length <- y[1:nrow(mating.list), .(h.length)]
  mating.list$dam.bw.june <- y[1:nrow(mating.list), .(bw.june)]
  mating.list$dam.bw.june.maternal <- y[1:nrow(mating.list), .(bw.june.maternal)]
  
  
  
  
  if ("f0" %in% colnames(y)) {
    mating.list$f0  <- y[1:nrow(mating.list), .(f0)]
  }
  if ("perm.env.bs" %in% colnames(y)) {
    mating.list$perm.env.bs  <- y[1:nrow(mating.list), .(perm.env.bs)]
  } # think this one is redundant
  mating.list$birthyear.dam  <- y[1:nrow(mating.list), .(birthyear)]
  setnames(mating.list, "f0", "f0.dam")
  if (year == 1) {
    setnames(mating.list, "f0.dam", "f0") 
    # this is because of in the breeding values of the first gen this is needed, silly really since its zero all over
  }
  # now to make the 2nd round of matings
  mating.list.allowed <- subset(mating.list, can.remate == 1)
  mating.list.allowed[, `:=`(IDX = 1:.N) , by = sire.id.1st]
  mating.list.allowed$sire.id.2nd       <-
    ifelse(
      mating.list.allowed$IDX <= mating.list.allowed$mating.willingness.2nd,
      mating.list.allowed$sire.id.1st   ,
      0
    )
  mating.list.allowed$test <- ifelse(mating.list.allowed$sire.id.1st == mating.list.allowed$sire.id.2nd,TRUE,FALSE)
  
  mating.list.allowed$semen.quality.2nd <-
    ifelse(
      mating.list.allowed$test == TRUE,
      mating.list.allowed$semen.quality.1st,
      0
    )
  mating.list.allowed$sire.fert.2nd     <-
    ifelse(
      mating.list.allowed$test == TRUE,
      mating.list.allowed$sire.fert.1st    ,
      0
    )
  mating.list.allowed$sire.bw.oct.male.2nd       <-
    ifelse(
      mating.list.allowed$test == TRUE,
      mating.list.allowed$sire.bw.oct.male.1st      ,
      0
    )
  mating.list.allowed$sire.bw.oct.female.2nd       <-
    ifelse(
      mating.list.allowed$test == TRUE,
      mating.list.allowed$sire.bw.oct.female.1st      ,
      0
    )
  
  mating.list.allowed$sire.bw.sept.male.2nd       <-
    ifelse(
      mating.list.allowed$test == TRUE,
      mating.list.allowed$sire.bw.sept.male.1st      ,
      0
    )
  mating.list.allowed$sire.bw.sept.female.2nd       <-
    ifelse(
      mating.list.allowed$test == TRUE,
      mating.list.allowed$sire.bw.sept.female.1st      ,
      0
    )
  
  mating.list.allowed$sire.skin.length.male.2nd       <-
    ifelse(
      mating.list.allowed$test == TRUE,
      mating.list.allowed$sire.skin.length.male.1st      ,
      0
    )
  mating.list.allowed$sire.skin.length.female.2nd       <-
    ifelse(
      mating.list.allowed$test == TRUE,
      mating.list.allowed$sire.skin.length.female.1st      ,
      0
    )
  
  mating.list.allowed$sire.skin.qual.2nd       <-
    ifelse(
      mating.list.allowed$test == TRUE,
      mating.list.allowed$sire.skin.qual.1st,
      0
    )
  mating.list.allowed$sire.live.qual.2nd       <-
    ifelse(
      mating.list.allowed$test == TRUE,
      mating.list.allowed$sire.live.qual.1st      ,
      0
    )
  mating.list.allowed$sire.h.length.2nd       <-
    ifelse(
      mating.list.allowed$test == TRUE,
      mating.list.allowed$sire.h.length.1st      ,
      0
    )
  mating.list.allowed$sire.bw.june.2nd       <-
    ifelse(
      mating.list.allowed$test == TRUE,
      mating.list.allowed$sire.bw.june.1st      ,
      0
    )
  
  mating.list.allowed$sire.bw.june.maternal.2nd       <-
    ifelse(
      mating.list.allowed$test == TRUE,
      mating.list.allowed$sire.bw.june.maternal.1st      ,
      0
    )
  
  
  
  
  #subset mating.list.allowed into two, those done and those left
  
  mating.list.remated <-
    subset (mating.list.allowed, sire.id.2nd != 0)
  mating.list.leftover <-
    subset (mating.list.allowed, sire.id.2nd == 0) 
  # those females who were allowed to
  # be remated with 1st male but male did not have enough to mate them all
  mating.list.notallowed <- subset(mating.list, can.remate == 0)
  x[, `:=`(matings.left = mating.willingness.2nd - mating.willingness.1st)] 
  # figures out how many matings the males performed that did have spares
  x$matings.left <-
    ifelse(x$matings.left < 0, 0, x$matings.left)  # if they have negs, they are done
  x <-
    subset(x, matings.left > 0 &
             can.remate == 1)  # only the males that have spare matings left
  
  set(mating.list.leftover,
      j = c("IDX","test"),
      value = NULL)
  mating.list.leftover <-
    rbind(mating.list.leftover, mating.list.notallowed)
  
  # here I should reorder the dams that are leftovers to prioritize younger dams and the ones mated with shitty males
  setkey(mating.list.leftover, dam.id)
  if (year >2 & use.blup.to.assort.mat == 1 & selection.method ==blup) {
    setorder(mating.list.leftover, -comb.ind)
  } else if (use.blup.to.assort.mat == 0) {
    setorder(mating.list.leftover, -birthyear.dam, -dam.live.score)}
  myvars <- c("dam.id",
              "sire.id.2nd",
              "semen.quality.2nd",
              "sire.fert.2nd",
              "sire.bw.oct.male.2nd",
              "sire.bw.oct.female.2nd",
              "sire.bw.sept.male.2nd",
              "sire.bw.sept.female.2nd",
              "sire.skin.length.male.2nd",
              "sire.skin.length.female.2nd",
              "sire.skin.qual.2nd",
              "sire.live.qual.2nd",
              "sire.bw.june.2nd",
              "sire.bw.june.maternal.2nd",
              "sire.h.length.2nd")
  mating.list.leftover.temp <-
    as.matrix(mating.list.leftover[,myvars,with=FALSE]) 
  myvars <- c("sire.id.2nd",
              "semen.quality.2nd",
              "sire.fert.2nd",
              "sire.bw.oct.male.2nd",
              "sire.bw.oct.female.2nd",
              "sire.bw.sept.male.2nd",
              "sire.bw.sept.female.2nd",
              "sire.skin.length.male.2nd",
              "sire.skin.length.female.2nd",
              "sire.skin.qual.2nd",
              "sire.live.qual.2nd",
              "sire.bw.june.2nd",
              "sire.bw.june.maternal.2nd",
              "sire.h.length.2nd")
  mating.list.leftover <- mating.list.leftover[,!myvars,with=FALSE]
  # way faster to loop through a matrix
  # this big loop has two modes, one if the selection method is blup and
  # another if the selection method is phenotypic
  # then each has three versions of the loop, one for year == 1, year ==2 and
  # year > 2. This is because of differing amounts of years in the matrix that
  # the loop accepts as input. Year 1 loops are always the same, regardless of 
  # selection method
  myvars <- c("id", # 1
              "semen.quality.2nd", #2
              "litter.size", # 3 
              "bw.oct.male", # 4
              "bw.oct.female", # 5
              "bw.sept.male", # 6
              "bw.sept.female", # 7
              "skin.length.male", # 8
              "skin.length.female", # 9
              "skin.qual", # 10
              "live.qual", # 11
              "live.score", # 12
              "h.length", # 13
              "bw.june", # 14
              "bw.june.maternal", # 15
              "matings.left") # 16
  setkey(x, id)
  if (year >2 & use.blup.to.assort.mat == 1 & selection.method ==blup) {
    setorder(x, -comb.ind)
  } else if (use.blup.to.assort.mat == 0 ) {
    setorder(x, -live.score)
  }
  x <- as.matrix(x[,myvars, with=FALSE])
  
  for (i in 1:nrow(x))  { #number of mating males 
    #print(i)
    for (j in 1:x[[i, 16]])  { # number of matings left
      s <-
        sum(x[1:i, 16]) 
      # keeps track of how many females the male has been assigned
      #         print(s)
      #     print(j)
      if (s < nrow(mating.list.leftover.temp)) {
        if ( i == 1) {
          # this is not a solution, need to make another controlling mechanism if the mating willingness
          #exceeds the number of females to be mated
          mating.list.leftover.temp[[s - (j - 1), 2]] <- x[[i, 1]]          # id of male
          mating.list.leftover.temp[[s - (j - 1), 3]] <- x[[i, 2]]          # semen.quality
          mating.list.leftover.temp[[s - (j - 1), 4]] <- x[[i, 3]]          # fertility of male
          mating.list.leftover.temp[[s - (j - 1), 5]] <- x[[i, 4]]          # body weight oct.male of male
          mating.list.leftover.temp[[s - (j - 1), 6]] <- x[[i, 5]]          # body weight oct.female of male
          mating.list.leftover.temp[[s - (j - 1), 7]] <- x[[i, 6]]          # body weight sept.male of male
          mating.list.leftover.temp[[s - (j - 1), 8]] <- x[[i, 7]]          # body weight sept.female of male
          mating.list.leftover.temp[[s - (j - 1), 9]] <- x[[i, 8]]          # body weight sept.male of male
          mating.list.leftover.temp[[s - (j - 1), 10]] <- x[[i, 9]]          # body weight sept.female of male
          mating.list.leftover.temp[[s - (j - 1), 11]] <- x[[i, 10]]          # skin.qual  of male
          mating.list.leftover.temp[[s - (j - 1), 12]] <- x[[i, 11]]          # live.qual  of male
          mating.list.leftover.temp[[s - (j - 1), 13]] <- x[[i, 14]]          # bw.june  of male
          mating.list.leftover.temp[[s - (j - 1), 14]] <- x[[i, 15]]          # bw.june.maternal  of male
          mating.list.leftover.temp[[s - (j - 1), 15]] <- x[[i, 13]]          # h.length  of male
          
        } else if ( i > 1) {
          t <- sum(x[1:(i-1), 16])+1
          mating.list.leftover.temp[[t+(j-1), 2]] <- x[[i, 1]]          # id of male
          mating.list.leftover.temp[[t+(j-1), 3]] <- x[[i, 2]]          # semen.quality
          mating.list.leftover.temp[[t+(j-1), 4]] <- x[[i, 3]]          # fertility of male
          mating.list.leftover.temp[[t+(j-1), 5]] <- x[[i, 4]]          # body weight oct.male of male
          mating.list.leftover.temp[[t+(j-1), 6]] <- x[[i, 5]]          # body weight oct.female of male
          mating.list.leftover.temp[[t+(j-1), 7]] <- x[[i, 6]]          # body weight sept.male of male
          mating.list.leftover.temp[[t+(j-1), 8]] <- x[[i, 7]]          # body weight sept.female of male
          mating.list.leftover.temp[[t+(j-1), 9]] <- x[[i, 8]]          # body weight sept.male of male
          mating.list.leftover.temp[[t+(j-1), 10]] <- x[[i, 9]]          # body weight sept.female of male
          mating.list.leftover.temp[[t+(j-1), 11]] <- x[[i, 10]]          # skin.qual  of male
          mating.list.leftover.temp[[t+(j-1), 12]] <- x[[i, 11]]          # live.qual  of male
          mating.list.leftover.temp[[t+(j-1), 13]] <- x[[i, 14]]          # bw.june  of male
          mating.list.leftover.temp[[t+(j-1), 14]] <- x[[i, 15]]          # bw.june.maternal  of male
          mating.list.leftover.temp[[t+(j-1), 15]] <- x[[i, 13]]          # h.length  of male
          
        }
      }
      if (s > nrow(mating.list.leftover.temp)) {
        while (sum(x[1:(i - 1) , 16]) + j <= nrow(mating.list.leftover.temp)) {
          # Here i solve the problem of if the males have more mating willingness than the number of females
          if (s < nrow(mating.list.leftover.temp)) {
            if ( i == 1) {
              # this is not a solution, need to make another controlling mechanism if the mating willingness
              #exceeds the number of females to be mated
              mating.list.leftover.temp[[s - (j - 1), 2]] <- x[[i, 1]]          # id of male
              mating.list.leftover.temp[[s - (j - 1), 3]] <- x[[i, 2]]          # semen.quality
              mating.list.leftover.temp[[s - (j - 1), 4]] <- x[[i, 3]]          # fertility of male
              mating.list.leftover.temp[[s - (j - 1), 5]] <- x[[i, 4]]          # body weight oct.male of male
              mating.list.leftover.temp[[s - (j - 1), 6]] <- x[[i, 5]]          # body weight oct.female of male
              mating.list.leftover.temp[[s - (j - 1), 7]] <- x[[i, 6]]          # body weight sept.male of male
              mating.list.leftover.temp[[s - (j - 1), 8]] <- x[[i, 7]]          # body weight sept.female of male
              mating.list.leftover.temp[[s - (j - 1), 9]] <- x[[i, 8]]          # body weight sept.male of male
              mating.list.leftover.temp[[s - (j - 1), 10]] <- x[[i, 9]]          # body weight sept.female of male
              mating.list.leftover.temp[[s - (j - 1), 11]] <- x[[i, 10]]          # skin.qual  of male
              mating.list.leftover.temp[[s - (j - 1), 12]] <- x[[i, 11]]          # live.qual  of male
              mating.list.leftover.temp[[s - (j - 1), 13]] <- x[[i, 14]]          # bw.june  of male
              mating.list.leftover.temp[[s - (j - 1), 14]] <- x[[i, 15]]          # bw.june.maternal  of male
              mating.list.leftover.temp[[s - (j - 1), 15]] <- x[[i, 13]]          # h.length  of male
            } else if ( i > 1) {
              t <- sum(x[1:(i-1), 16])+1
              mating.list.leftover.temp[[t+(j-1), 2]] <- x[[i, 1]]          # id of male
              mating.list.leftover.temp[[t+(j-1), 3]] <- x[[i, 2]]          # semen.quality
              mating.list.leftover.temp[[t+(j-1), 4]] <- x[[i, 3]]          # fertility of male
              mating.list.leftover.temp[[t+(j-1), 5]] <- x[[i, 4]]          # body weight oct.male of male
              mating.list.leftover.temp[[t+(j-1), 6]] <- x[[i, 5]]          # body weight oct.female of male
              mating.list.leftover.temp[[t+(j-1), 7]] <- x[[i, 6]]          # body weight sept.male of male
              mating.list.leftover.temp[[t+(j-1), 8]] <- x[[i, 7]]          # body weight sept.female of male
              mating.list.leftover.temp[[t+(j-1), 9]] <- x[[i, 8]]          # body weight sept.male of male
              mating.list.leftover.temp[[t+(j-1), 10]] <- x[[i, 9]]          # body weight sept.female of male
              mating.list.leftover.temp[[t+(j-1), 11]] <- x[[i, 10]]          # skin.qual  of male
              mating.list.leftover.temp[[t+(j-1), 12]] <- x[[i, 11]]          # live.qual  of male
              mating.list.leftover.temp[[t+(j-1), 13]] <- x[[i, 14]]          # bw.june  of male
              mating.list.leftover.temp[[t+(j-1), 14]] <- x[[i, 15]]          # bw.june.maternal  of male
              mating.list.leftover.temp[[t+(j-1), 15]] <- x[[i, 13]]          # h.length  of male
              
            }
          }
          
          break
        }
      }
    }
  }
  
  
  mating.list.leftover <- as.data.table(mating.list.leftover)
  mating.list.leftover.temp <- as.data.table(mating.list.leftover.temp)
  mating.list.leftover <- merge(mating.list.leftover, mating.list.leftover.temp, by="dam.id")
  x <- as.data.table(x)
  set(mating.list.remated, j = c("IDX","test"), value = NULL) # need to remove the counter to bind the mated ones to the rest
  
  mating.list <-
    rbind(mating.list.remated, mating.list.leftover) # merge the now completed mating list together
  mating.list[, `:=`(
    semen.quality = ifelse (
      mating.list$semen.quality.2nd == 0 & mating.list$sire.id.2nd ==0,
      mating.list$semen.quality.1st,
      ifelse(mating.list$semen.quality.2nd == 0 & mating.list$sire.id.2nd != 0, 0,mating.list$semen.quality.2nd)
    ),
    barren = # this rather cumbersome function is to have different probs of barren acc. to mating method
      ifelse (
        year - mating.list$birthyear.dam == 1 &
          mating.list$sire.id.2nd != 0 ,
        rbinom(nrow(mating.list), 1, pr.barren.double.mating.yearling),
        ifelse(
          year - mating.list$birthyear.dam != 1 &
            mating.list$sire.id.2nd != 0,
          rbinom(nrow(mating.list), 1, pr.barren.double.mating.old),
          ifelse(
            year - mating.list$birthyear.dam == 1 &
              mating.list$sire.id.2nd == 0,
            rbinom(nrow(mating.list), 1, pr.barren.one.mating.yearling),
            ifelse(
              year - mating.list$birthyear.dam != 1 &
                mating.list$sire.id.2nd == 0,
              rbinom(nrow(mating.list), 1, pr.barren.one.mating.old),
              0
            )
          )
        )
      )
  )]
  # remove columns who're not needed at this point
  set(
    mating.list,
    j = c(
      "mating.willingness.1st",
      "mating.willingness.2nd",
      "can.remate",
      "semen.quality.1st",
      "semen.quality.2nd",
      "sire.live.score.2nd",
      "sire.live.score.1st"
    ),
    value = NULL
  )
  
  return(mating.list)
}
mate <-compiler::cmpfun(mate,options= c(suppressAll=TRUE)) # performance boost
############### make gen 0 pedigree & and some bookkeeping ##############
# This function makes the pedfile for the first generation which is then expanded upon each generation
MakePedfileGen0 <- function(x,y) { #x = gen0 females, y = effgen0.males
setkey(x, id)
setkey(y, id)
  pedfile <- x[,.(id,sire.id,dam.id,sex,birthyear)] #use .(colnames) instead of c(colnames) in data.table
  l <- list (pedfile, y[,.(id,sire.id,dam.id,sex,birthyear)])
  pedfile <- rbindlist( l,use.names = TRUE ) 
  return (pedfile)
} 
############### Breeding value of offspring ###############################
# This function creates the kits for the first generation only, since the 
# subsequent generation is more complex this is the only instance of it
MakeKitsGen0 <- function (x,y, z) { #x = mating.list, y= effgen0.males, z = year
  kit.list <- x[rep(seq(nrow(x)), obs_fert), 
            c("dam.id",
              "sire.id.1st",
              "dam.fert",
              "sire.fert.1st",
              "sire.bw.oct.male.1st",
              "sire.bw.oct.female.1st",
              "sire.bw.sept.male.1st",
              "sire.bw.sept.female.1st",
              "sire.skin.length.male.1st",
              "sire.skin.length.female.1st",
              "sire.skin.qual.1st",
              "sire.live.qual.1st",
              "sire.bw.june.1st",
              "sire.bw.june.maternal.1st",
              "sire.h.length.1st",
              "f0",
              "obs_fert",
              "dam.bw.oct.male",
              "dam.bw.oct.female",
              "dam.bw.sept.male",
              "dam.bw.sept.female",
              "dam.skin.length.male",
              "dam.skin.length.female",
              "dam.skin.qual",
              "dam.live.qual",
              "dam.bw.june",
              "dam.bw.june.maternal",
              "dam.h.length",
              "specific.env.bw",
              "birthyear.dam",
              "sire.id.2nd",
              "sire.fert.2nd",
              "sire.bw.oct.male.2nd",
              "sire.bw.oct.female.2nd",
              "sire.bw.sept.male.2nd",
              "sire.bw.sept.female.2nd",
              "sire.skin.length.male.2nd",
              "sire.skin.length.female.2nd",
              "sire.skin.qual.2nd",
              "sire.live.qual.2nd",
              "sire.bw.june.2nd",
              "sire.bw.june.maternal.2nd",
              "sire.h.length.2nd"
            ) 
            , with=F] #specify which columns to incl.
  id <- seq(1:sum(x$obs_fert)) + max(y$id) # makes ID
  birthyear <- rep (1, sum(x$obs_fert)) # makes birthyear
  sex <- rbinom(sum(x$obs_fert),1,0.5)+1 # makes sex, TODO check if this is true
  true.sire <- numeric(nrow(kit.list))
  true.sire.check <- numeric(nrow(kit.list))
  mendelian <- as.data.table(rmvnorm(sum(x$obs_fert),sigma=sigma.new,method="svd"))
  colnames(mendelian)     <-
    c(
      "mend.live.qual",
      "mend.h.length",
      "mend.bw.oct.male",
      "mend.bw.oct.female",
      "mend.skin.qual",
      "mend.skin.length.male",
      "mend.skin.length.female",
      "mend.litter.size",
      "mend.bw.sept.male",
      "mend.bw.sept.female",
      "mend.bw.june.maternal",
      "mend.bw.june"
    )
  
  kit.list <- cbind (id,kit.list,  birthyear, sex, mendelian, true.sire,true.sire.check) # binds id, sex and birthyear to data.table
  
  kit.list$sire.id.2nd  <-
    ifelse(kit.list$sire.id.2nd == 0,
           kit.list$sire.id.1st,
           kit.list$sire.id.2nd)
  # puts the 1st sire into the 2nd sire column for single mated females
  # here i use the true sire check in order to speed up, replace it later
  kit.list$true.sire.check <- ifelse(kit.list$sire.id.2nd == kit.list$sire.id.1st,TRUE,FALSE)
  kit.list$sire.fert.2nd <-
    ifelse(
      kit.list$true.sire.check == TRUE,
      kit.list$sire.fert.1st,
      kit.list$sire.fert.2nd
    )
  kit.list$sire.bw.oct.male.2nd <-
    ifelse(
      kit.list$true.sire.check == TRUE,
      kit.list$sire.bw.oct.male.1st,
      kit.list$sire.bw.oct.male.2nd
    )
  kit.list$sire.bw.oct.female.2nd <-
    ifelse(
      kit.list$true.sire.check == TRUE,
      kit.list$sire.bw.oct.female.1st,
      kit.list$sire.bw.oct.female.2nd
    )
  
  kit.list$sire.bw.sept.male.2nd <-
    ifelse(
      kit.list$true.sire.check == TRUE,
      kit.list$sire.bw.sept.male.1st,
      kit.list$sire.bw.sept.male.2nd
    )
  kit.list$sire.bw.sept.female.2nd <-
    ifelse(
      kit.list$true.sire.check == TRUE,
      kit.list$sire.bw.sept.female.1st,
      kit.list$sire.bw.sept.female.2nd
    )
  
  kit.list$sire.skin.length.male.2nd <-
    ifelse(
      kit.list$true.sire.check == TRUE,
      kit.list$sire.skin.length.male.1st,
      kit.list$sire.skin.length.male.2nd
    )
  kit.list$sire.skin.length.female.2nd <-
    ifelse(
      kit.list$true.sire.check == TRUE,
      kit.list$sire.skin.length.female.1st,
      kit.list$sire.skin.length.female.2nd
    )
  kit.list$sire.skin.qual.2nd <-
    ifelse(
      kit.list$true.sire.check == TRUE,
      kit.list$sire.skin.qual.1st,
      kit.list$sire.skin.qual.2nd
    )
  kit.list$sire.live.qual.2nd <-
    ifelse(
      kit.list$true.sire.check == TRUE,
      kit.list$sire.live.qual.1st,
      kit.list$sire.live.qual.2nd
    )
  kit.list$sire.h.length.2nd <-
    ifelse(
      kit.list$true.sire.check == TRUE,
      kit.list$sire.h.length.1st,
      kit.list$sire.h.length.2nd
    )
  kit.list$sire.bw.june.2nd <-
    ifelse(
      kit.list$true.sire.check == TRUE,
      kit.list$sire.bw.june.1st,
      kit.list$sire.bw.june.2nd
    )
  kit.list$sire.bw.june.maternal.2nd <-
    ifelse(
      kit.list$true.sire.check == TRUE,
      kit.list$sire.bw.june.maternal.1st,
      kit.list$sire.bw.june.maternal.2nd
    )
  kit.list$true.sire.check <- ifelse( kit.list$sire.id.1st != kit.list$sire.id.2nd, rbinom(nrow(kit.list), 1, 0.85), 1) # 85% chance that the kits are sired by 2nd mating
  kit.list$true.sire <- ifelse( kit.list$true.sire.check == 0, kit.list$sire.id.1st, kit.list$sire.id.2nd)
  kit.list[, `:=`(true.sire.fert = 
                ifelse(kit.list$true.sire == kit.list$sire.id.2nd, 
                       kit.list$sire.fert.2nd, kit.list$sire.fert.1st), 
              true.sire.bw.oct.male = 
                ifelse(kit.list$true.sire == kit.list$sire.id.2nd, 
                       kit.list$sire.bw.oct.male.2nd, kit.list$sire.bw.oct.male.1st),
              true.sire.bw.oct.female = 
                ifelse(kit.list$true.sire == kit.list$sire.id.2nd, 
                       kit.list$sire.bw.oct.female.2nd, kit.list$sire.bw.oct.female.1st),
              
              true.sire.bw.sept.male = 
                ifelse(kit.list$true.sire == kit.list$sire.id.2nd, 
                       kit.list$sire.bw.sept.male.2nd, kit.list$sire.bw.sept.male.1st),
              true.sire.bw.sept.female = 
                ifelse(kit.list$true.sire == kit.list$sire.id.2nd, 
                       kit.list$sire.bw.sept.female.2nd, kit.list$sire.bw.sept.female.1st),
              
              true.sire.skin.length.male = 
                ifelse(kit.list$true.sire == kit.list$sire.id.2nd, 
                       kit.list$sire.skin.length.male.2nd, kit.list$sire.skin.length.male.1st),
              true.sire.skin.length.female = 
                ifelse(kit.list$true.sire == kit.list$sire.id.2nd, 
                       kit.list$sire.skin.length.female.2nd, kit.list$sire.skin.length.female.1st),
              
              true.sire.skin.qual = 
                ifelse(kit.list$true.sire == kit.list$sire.id.2nd, 
                       kit.list$sire.skin.qual.2nd, kit.list$sire.skin.qual.1st),
              true.sire.live.qual = 
                ifelse(kit.list$true.sire == kit.list$sire.id.2nd, 
                       kit.list$sire.live.qual.2nd, kit.list$sire.live.qual.1st),
              true.sire.h.length = 
                ifelse(kit.list$true.sire == kit.list$sire.id.2nd, 
                       kit.list$sire.h.length.2nd, kit.list$sire.h.length.1st),
              true.sire.bw.june = 
                ifelse(kit.list$true.sire == kit.list$sire.id.2nd, 
                       kit.list$sire.bw.june.2nd, kit.list$sire.bw.june.1st),
              true.sire.bw.june.maternal = 
                ifelse(kit.list$true.sire == kit.list$sire.id.2nd, 
                       kit.list$sire.bw.june.maternal.2nd, kit.list$sire.bw.june.maternal.1st)
              
  )]

  kit.list$mend.bw.june <-  ifelse(kit.list$sex==1, kit.list$mend.bw.june*sqrt(0.5*var.bw.june.male),
                                   kit.list$mend.bw.june*sqrt(0.5*var.bw.june.female))
  
  kit.list[, `:=`(litter.size = 0.5*(dam.fert + true.sire.fert) + 
                mend.litter.size*(sqrt(0.5*variance.fertility)),  # Breeding value of offspring, littersize
               perm.env.ls = rnorm(sum(x$obs_fert))*sqrt(var.perm.env.ls), # perm env for litter size
              bw.oct.male = 0.5*(dam.bw.oct.male + true.sire.bw.oct.male) + 
                mend.bw.oct.male*sqrt(0.5*var.bw.oct.male),
                bw.oct.female = 0.5*(dam.bw.oct.female + true.sire.bw.oct.female) + 
                mend.bw.oct.female*sqrt(0.5*var.bw.oct.female),
              bw.sept.male = 0.5*(dam.bw.sept.male + true.sire.bw.sept.male) + 
                mend.bw.sept.male*sqrt(0.5*var.bw.sept.male),
              bw.sept.female = 0.5*(dam.bw.sept.female + true.sire.bw.sept.female) + 
                mend.bw.sept.female*sqrt(0.5*var.bw.sept.female),
              skin.length.male = 0.5*(dam.skin.length.male + true.sire.skin.length.male) + 
                mend.skin.length.male*sqrt(0.5*var.skin.length.male),
              skin.length.female = 0.5*(dam.skin.length.female + true.sire.skin.length.female) + 
                mend.skin.length.female*sqrt(0.5*var.skin.length.female),
              live.qual = 0.5*(dam.live.qual + true.sire.live.qual)+ sqrt(0.5*var.live.qual.gen)*(mend.live.qual),
              skin.qual = 0.5*(dam.skin.qual + true.sire.skin.qual)+ sqrt(0.5)*(mend.skin.qual),
              h.length = 0.5*(dam.h.length + true.sire.h.length)+ sqrt(0.5)*(mend.h.length),
              bw.june= 0.5*(dam.bw.june + true.sire.bw.june + mend.bw.june),
              bw.june.maternal =0.5*(dam.bw.june.maternal + true.sire.bw.june.maternal + sqrt(0.5)*mend.bw.june.maternal)
              )]# Breeding value of offspring, body size
  #  kit.list <- count.sex.siblings(kit.list) # calls function to count offspring NOT NEEDED ATM
  
  setnames(kit.list, c("obs_fert","sire.id.2nd"), c("own_littersize","sire.assumed")) # renames obs_fert to own littersize of kits
  kit.list$dam.age <- ifelse( z - kit.list$birthyear.dam > 1, 1,0 )
  
  # kit.list$phenotype.bw.oct <- ifelse( kit.list$sex == 1,MakePhenotypesBWMalesOct(mean.body.size.male.oct , kit.list$bw.oct , kit.list$specific.env.bw , kit.list$own_littersize, kit.list$dam.age,x  )
  #                                  , MakePhenotypesBWFemalesOct(mean.body.size.female.oct , kit.list$bw.oct , kit.list$specific.env.bw,kit.list$own_littersize,x ))
  # kit.list$phenotype.bw.sept <- ifelse( kit.list$sex == 1,
  #                                   MakePhenotypesBWMalesSept(mean.body.size.male.sept , kit.list$bw.sept , kit.list$specific.env.bw , kit.list$own_littersize, kit.list$dam.age,x  )
  #                                   , MakePhenotypesBWFemalesOct(mean.body.size.female.sept , kit.list$bw.sept , kit.list$specific.env.bw,kit.list$own_littersize,x ))
  # kit.list$phenotype.live.qual <- kit.list$live.qual  + rnorm(nrow(kit.list))*sqrt(var.live.qual.res)
  kit.list[, `:=`( 
    phenotype.bw.oct = ifelse(
      kit.list$sex == 1,
      MakePhenotypesBWMalesOct(
        mean.body.size.male.oct ,
        kit.list$bw.oct.male ,
        kit.list$specific.env.bw ,
        kit.list$own_littersize,
        kit.list$dam.age,
        x
      )
      ,
      MakePhenotypesBWFemalesOct(
        mean.body.size.female.oct ,
        kit.list$bw.oct.fe ,
        kit.list$specific.env.bw,
        kit.list$own_littersize,
        x
      )
    ),
    phenotype.bw.sept = ifelse(
      kit.list$sex == 1,
      MakePhenotypesBWMalesSept(
        mean.body.size.male.sept ,
        kit.list$bw.sept.male ,
        kit.list$specific.env.bw ,
        kit.list$own_littersize,
        kit.list$dam.age,
        x
      )
      ,
      MakePhenotypesBWFemalesSept(
        mean.body.size.female.sept ,
        kit.list$bw.sept.female ,
        kit.list$specific.env.bw,
        kit.list$own_littersize,
        x
      )
    ),
    #
    phenotype.bw.june = ifelse(
      kit.list$sex == 1,
      MakePhenotypesBWMalesJune(
        mean.bw.june.male ,
        kit.list$bw.june ,
        kit.list$specific.env.bw ,
        kit.list$own_littersize,
        kit.list$dam.age,
        x,
        kit.list$dam.bw.june.maternal
      )
      ,
      MakePhenotypesBWFemalesJune(
        mean.bw.june.female ,
        kit.list$bw.june ,
        kit.list$specific.env.bw,
        kit.list$own_littersize,
        kit.list$dam.age,
        x,
        kit.list$dam.bw.june.maternal
      )
    ),
    phenotype.live.qual = live.qual  + rnorm(nrow(kit.list)) *
      sqrt(var.live.qual.res),
    phenotype.skin.length = ifelse(
      sex == 1,
      mean.skin.length.male + skin.length.male + rnorm(1)*sqrt(var.skin.length.res.male)+
        specific.env.bw*sqrt(var.skin.length.c.male),
      mean.skin.length.female + skin.length.female + rnorm(1)*sqrt(var.skin.length.res.female)+
        specific.env.bw*sqrt(var.skin.length.c.female) ),
    phenotype.skin.qual = 
      skin.qual + rnorm(nrow(kit.list))*sqrt(var.skin.qual.res)
    )
  ]  

    set( kit.list, j=which(colnames(kit.list) %in% c(
    "sire.fert.1st",
    "sire.bw.oct.male.1st",
    "sire.bw.oct.female.1st",
    "sire.bw.sept.male.1st",
    "sire.bw.sept.female.1st",
    "sire.skin.length.male.1st",
    "sire.skin.length.female.1st",
    "sire.skin.qual.1st",
    "sire.live.qual.1st",
    "sire.bw.june.1st",
    "sire.bw.june.maternal.1st",
    "sire.h.length.1st",
    "sire.id.1st",
    "dam.fert",
    "dam.bw.oct.male",
    "dam.bw.oct.female",
    "dam.bw.sept.male",
    "dam.bw.sept.female",
    "dam.skin.length.male",
    "dam.skin.length.female",
    "dam.skin.qual",
    "dam.live.qual",
    "dam.bw.june",
    "dam.bw.june.maternal",
    "dam.h.length",
    "sire.bw.oct.1st",
    "sire.bw.sept.1st",
    "specific.env.bw",
    "sire.fert.2nd",
    "sire.bw.oct.male.2nd",
    "sire.bw.oct.female.2nd",
    "sire.bw.sept.male.2nd",
    "sire.bw.sept.female.2nd",
    "sire.skin.length.male.2nd",
    "sire.skin.length.female.2nd",
    "sire.skin.qual.2nd",
    "sire.live.qual.2nd",
    "sire.bw.june.2nd",
    "sire.bw.june.maternal.2nd",
    "sire.h.length.2nd",
    "mend.bw.oct.male", 
    "mend.bw.oct.female", 
    "mend.bw.sept.male", 
    "mend.bw.sept.female",
    "mend.bw.june.maternal",
    "mend.bw.june",
    "mend.litter.size", 
    "mend.live.qual",
    "mend.skin.qual", 
    "mend.skin.length.male",
    "mend.skin.length.female",
    "true.sire.check",
    "true.sire.fert",
    "true.sire.bw.oct.male",
    "true.sire.bw.oct.female",
    "true.sire.bw.sept.male",
    "true.sire.bw.sept.female",
    "true.sire.skin.length.male",
    "true.sire.skin.length.female",
    "true.sire.skin.qual",
    "true.sire.live.qual",
    "true.sire.h.length",
    "true.sire.bw.june",
    "true.sire.bw.june.maternal",
    "dam.age"
  )) , value=NULL ) # removes bv of parents
  if (qual.classes == 5) {
    truncs <- qnorm(
      p = c(0.05, 0.3, 0.7, 0.95),
      mean = mean(kit.list$live.qual),
                  sd = sqrt(var(
                    kit.list$live.qual
                  )),
      lower.tail = TRUE
    )
    kit.list[, `:=`(live.score = ifelse(
      live.qual >= truncs[4],
      5,
      ifelse(
        truncs[3] < live.qual & live.qual <= truncs[4],
        4,
        ifelse(
          live.qual > truncs[2] & live.qual <= truncs[3],
          3,
          ifelse(
            live.qual > truncs[1] & live.qual <= truncs[2],
            2,
            ifelse(live.qual <=
                     truncs[1], 1, 0)
          )
        )
      )
    ))]
  } else if (qual.classes == 10) {
    truncs <-
      qnorm(
        p = c(0.01, 0.05, 0.15, 0.3, 0.5, 0.7, 0.85, 0.95, 0.99),
        mean = mean(kit.list$live.qual),
                    sd = sqrt(var(
                      kit.list$live.qual
                    )),
        lower.tail = TRUE
      )
    kit.list[, `:=`(live.score =
                      ifelse(
                        live.qual >= truncs[9],
                        10,
                        ifelse(
                          truncs[8] < live.qual & live.qual <= truncs[9],
                          9,
                          ifelse(
                            live.qual > truncs[7] & live.qual <= truncs[8],
                            8,
                            ifelse(
                              live.qual > truncs[6] & live.qual <= truncs[7],
                              7,
                              ifelse(
                                live.qual > truncs[5] & live.qual <= truncs[6],
                                6,
                                ifelse(
                                  live.qual > truncs[4] & live.qual <= truncs[5],
                                  5,
                                  ifelse(
                                    live.qual > truncs[3] & live.qual <= truncs[4],
                                    4,
                                    ifelse(
                                      live.qual > truncs[2] & live.qual <= truncs[3],
                                      3,
                                      ifelse(
                                        live.qual > truncs[1] & live.qual <= truncs[2],
                                        2,
                                        ifelse(live.qual <=
                                                 truncs[1], 1, 0)
                                      )
                                    )
                                  )
                                )
                              )
                            )
                          )
                        )
                      ))]
  }
  return(kit.list)
}
############### Phenotypic Selection of old females ###############################################
# currently they are only truncated on their litter size phenotype
PhenoSelectionOldFemales <- function (y,x, year) { # y = gen0.females, x = mating.list
  setkey(x, obs_fert)
  setorder(x, -obs_fert)
  if ("birthyear.dam" %in% colnames(x)) {
    x <- subset(x, year -birthyear.dam < max.age.females, 
                select=c("dam.id","obs_fert")) 
  }
  if("obs_fert" %in% colnames(y)) {
    set( y, j=which(colnames(y) %in% 
                                "obs_fert")  , value=NULL )
  }
  if("dam.age" %in% colnames(y)) {
    set( y, j=which(colnames(y) %in% 
                      "dam.age")  , value=NULL )
  }
  
  old.females <- x[1:(n.females*prop.oldfemales),]
  setnames(old.females, "dam.id", "id")
  old.females <- merge(old.females, y, by="id")
  
  # if("phenotype.bw.oct" %in% colnames(old.females)){

#   } else {
# #     if (qual.classes == 5){
# #       truncs <- qnorm(p=c(0.05,0.3,0.7,0.95), 
# #                       mean=mean(old.females$live.qual,
# #                                 sd= sqrt(var(old.females$live.qual))),
# #                       lower.tail=TRUE)
# #       old.females[,`:=`(live.score= ifelse(live.qual >= truncs[4],5,
# # ifelse(truncs[3] < live.qual & live.qual <= truncs[4],4,
# # ifelse(live.qual > truncs[2] & live.qual<=truncs[3],3,
# # ifelse(live.qual > truncs[1] & live.qual <=truncs[2],2,
# # ifelse(live.qual <=truncs[1],1,0
# # ))))))]
# #     } else if (qual.classes == 10) {
# #       truncs <- qnorm(p=c(0.01, 0.05, 0.15, 0.3, 0.5, 0.7, 0.85, 0.95, 0.99), 
# #                       mean=mean(old.females$live.qual,
# #                                 sd= sqrt(var(old.females$live.qual))),
# #                       lower.tail=TRUE)
# #       old.females[,`:=`(live.score= 
# # ifelse(live.qual >= truncs[9],10,
# # ifelse(truncs[8] < live.qual & live.qual <= truncs[9],9,
# # ifelse(live.qual > truncs[7] & live.qual<=truncs[8],8,
# # ifelse(live.qual > truncs[6] & live.qual <=truncs[7],7,
# # ifelse(live.qual > truncs[5] & live.qual <=truncs[6],6,
# # ifelse(live.qual > truncs[4] & live.qual <=truncs[5],5,
# # ifelse(live.qual > truncs[3] & live.qual <=truncs[4],4,
# # ifelse(live.qual > truncs[2] & live.qual <=truncs[3],3,
# # ifelse(live.qual > truncs[1] & live.qual <=truncs[2],2,
# # ifelse(live.qual <=truncs[1],1,0 
# # )))))))))))]
# #     }
#     
# old.females[,`:=`(phenotype.bw.oct = 
# mean.body.size.female.oct + bw.oct +
# rnorm(nrow(old.females))*sqrt(var.c.bw.oct.female),
# phenotype.bw.sept = 
# mean.body.size.female.sept + bw.sept +
# rnorm(nrow(old.females))*sqrt(var.c.bw.oct.female)
#     )]
#   }
  if("sire.id" %in% colnames(old.females)) {
    setnames(old.females, "sire.id", "sire.assumed")
    old.females[,`:=`(true.sire = old.females$sire.assumed)]
  }
  if ("own_littersize" %in% colnames(old.females)) {
    
  } else {old.females[,`:=`(own_littersize = 0)]}
  if ("mating.will.1st.round" %in% colnames(old.females)) {
    set( old.females, j=which(colnames(old.females) %in% 
                                c("mating.will.1st.round","mating.will.2nd.round"))  , value=NULL )
    
  }
  if (year == 1 ) {
    old.females$f0 <- 0
  }
  return(old.females)
}  
############### Phenotypic Selection of yearling females ############
PhenoSelectionFemaleKits <- function (x,y) { # x = kit.list, y = y  
  truncation.point <-  quantile( x$own_littersize,  probs =  quantile.setting ) 
  selection.candidates.females <- subset(x, own_littersize >= truncation.point) # throw away the smallest litters
  selection.candidates.females <-  subset( selection.candidates.females,  sex  ==   2) # take the female kits
  if (weighing.method == oct){
    truncation.point <-  quantile( selection.candidates.females$phenotype.bw.oct,  probs =  quantile.setting.bw ) 
    selection.candidates.females <- subset(selection.candidates.females, phenotype.bw.oct >= truncation.point) # throw away the smallest litters
  } else if (weighing.method == sept) {
    truncation.point <-  quantile( selection.candidates.females$phenotype.bw.sept,
                                   probs =  quantile.setting.bw ) 
    selection.candidates.females <- subset(selection.candidates.females, 
                                           phenotype.bw.sept >= truncation.point) 
  }
  if (qual.classes == 5){
    truncs <- qnorm(p=c(0.05,0.3,0.7,0.95), 
                    mean=mean(selection.candidates.females$live.qual,
                              sd= sqrt(var(selection.candidates.females$live.qual))),
                    lower.tail=TRUE)
    selection.candidates.females[,`:=`(live.score= ifelse(live.qual >= truncs[4],5,
     ifelse(truncs[3] < live.qual & live.qual <= truncs[4],4,
      ifelse(live.qual > truncs[2] & live.qual<=truncs[3],3,
       ifelse(live.qual > truncs[1] & live.qual <=truncs[2],2,
        ifelse(live.qual <=truncs[1],1,0
         ))))))]
  } else if (qual.classes == 10) {
    truncs <- qnorm(p=c(0.01, 0.05, 0.15, 0.3, 0.5, 0.7, 0.85, 0.95, 0.99), 
                    mean=mean(selection.candidates.females$live.qual,
                              sd= sqrt(var(selection.candidates.females$live.qual))),
                    lower.tail=TRUE)
    selection.candidates.females[,`:=`(live.score= 
     ifelse(live.qual >= truncs[9],10,
      ifelse(truncs[8] < live.qual & live.qual <= truncs[9],9,
       ifelse(live.qual > truncs[7] & live.qual<=truncs[8],8,
        ifelse(live.qual > truncs[6] & live.qual <=truncs[7],7,
         ifelse(live.qual > truncs[5] & live.qual <=truncs[6],6,
          ifelse(live.qual > truncs[4] & live.qual <=truncs[5],5,
           ifelse(live.qual > truncs[3] & live.qual <=truncs[4],4,
            ifelse(live.qual > truncs[2] & live.qual <=truncs[3],3,
             ifelse(live.qual > truncs[1] & live.qual <=truncs[2],2,
              ifelse(live.qual <=truncs[1],1,0 
               )))))))))))]
  }
  
  setkey(selection.candidates.females, live.score)           # in order to speed up ordering 
  setorder(selection.candidates.females,-live.score)         # order the kits according to body size 
  next.gen <- selection.candidates.females[1:(n.females-nrow(y)),]   # take kits with highest qual
  setkey(next.gen, id)
  next.gen[next.gen,obs_fert:=0]  
  if("f0.dam" %in% colnames(next.gen)) {
    set( next.gen, j=which(colnames(next.gen) %in% 
                             "f0.dam")  , value=NULL )
  }
  if("prodyear" %in% colnames(next.gen)) {
    set( next.gen, j=which(colnames(next.gen) %in% 
                             "prodyear")  , value=NULL )
  }
  
  
  return (next.gen)
}
############### Selection of yearling males in 1st gen ###############################
PhenoSelectionMaleKits <- function (x) { # x = kit.list 
  truncation.point <-  quantile( x$own_littersize,  probs =  quantile.setting ) 
  selection.candidates.males <- subset(x, own_littersize >= truncation.point) # throw away the smallest litters
  selection.candidates.males <-  subset( selection.candidates.males,  sex  ==   1) # take the male kits
  if (weighing.method == oct){
    truncation.point <-  quantile( selection.candidates.males$phenotype.bw.oct,  probs =  quantile.setting.bw ) 
    selection.candidates.males <- subset(selection.candidates.males, phenotype.bw.oct >= truncation.point) # throw away the smallest litters
  } else if (weighing.method == sept) {
    truncation.point <-  quantile( selection.candidates.males$phenotype.bw.sept,
                                   probs =  quantile.setting.bw ) 
    selection.candidates.males <- subset(selection.candidates.males, 
                                         phenotype.bw.sept >= truncation.point) 
  }
  if (qual.classes == 5){
    truncs <- qnorm(p=c(0.05,0.3,0.7,0.95), 
                    mean=mean(selection.candidates.males$live.qual,
                              sd= sqrt(var(selection.candidates.males$live.qual))),
                    lower.tail=TRUE)
    selection.candidates.males[,`:=`(live.score= ifelse(live.qual >= truncs[4],5,
     ifelse(truncs[3] < live.qual & live.qual <= truncs[4],4,
      ifelse(live.qual > truncs[2] & live.qual<=truncs[3],3,
       ifelse(live.qual > truncs[1] & live.qual <=truncs[2],2,
        ifelse(live.qual <=truncs[1],1,0
         ))))))]
  } else if (qual.classes == 10) {
    truncs <- qnorm(p=c(0.01, 0.05, 0.15, 0.3, 0.5, 0.7, 0.85, 0.95, 0.99), 
                    mean=mean(selection.candidates.males$live.qual,
                              sd= sqrt(var(selection.candidates.males$live.qual))),
                    lower.tail=TRUE)
    selection.candidates.males[,`:=`(live.score= 
      ifelse(live.qual >= truncs[9],10,
       ifelse(truncs[8] < live.qual & live.qual <= truncs[9],9,
        ifelse(live.qual > truncs[7] & live.qual<=truncs[8],8,
         ifelse(live.qual > truncs[6] & live.qual <=truncs[7],7,
          ifelse(live.qual > truncs[5] & live.qual <=truncs[6],6,
           ifelse(live.qual > truncs[4] & live.qual <=truncs[5],5,
            ifelse(live.qual > truncs[3] & live.qual <=truncs[4],4,
             ifelse(live.qual > truncs[2] & live.qual <=truncs[3],3,
              ifelse(live.qual > truncs[1] & live.qual <=truncs[2],2,
               ifelse(live.qual <=truncs[1],1,0 
               )))))))))))]
  }
  
  setkey(selection.candidates.males, live.score)           # in order to speed up ordering 
  if (weighing.method == sept){
  setorder(selection.candidates.males,-live.score,-phenotype.bw.sept)         # order the kits according to body size 
  } else if (weighing.method == oct) {
    setorder(selection.candidates.males,-live.score,-phenotype.bw.oct)         # order the kits according to body size 
    
} 
  next.gen <- selection.candidates.males[1:(n.males),]   # take kits with highest qual
  setkey(next.gen, id)
  if("f0.dam" %in% colnames(next.gen)) {
    set( next.gen, j=which(colnames(next.gen) %in% 
                             "f0.dam")  , value=NULL )
  }
  if("prodyear" %in% colnames(next.gen)) {
    set( next.gen, j=which(colnames(next.gen) %in% 
                             "prodyear")  , value=NULL )
  }
  
  
  return (next.gen)
}
############### Noselection of old females  ###############
nosel.old.females <- function (x) {
    old.females <- x[is.element(x$id, sample(x$id, n.females*prop.oldfemales)),]
    #safety valve to remove too old females
    # here i cbind the observation for litter size in this year to the old.females list
    old.females <- merge(old.females, mating.list, by.x="id", by.y="dam.id", all.x=TRUE)
    set( old.females, j=which(colnames(old.females) %in% 
                                c("sire.id.y","barren","dam.fert","sire.fert")) 
         , value=NULL )
     setnames(old.females, "sire.id.x", "sire.id")
     old.females[is.na(obs_fert), obs_fert:= 0]   
     old.females[is.na(f0), f0:= 0]   
        old.females[old.females,own_littersize:=0]  # adds a column named own littersize with value 0 
        return (old.females)
  }
############### No selection of yearlings ###################
nosel.yearlings.females <- function (x) {
    selection.candidates.females <-  subset( x,  sex  ==   2  ) 
    # kit.file_df in generation after first
    # We choose females to fit n.females, based on prop of old females
    
    next.gen <- sample(selection.candidates.females$id, (n.females-nrow(old.females)))
    next.gen <- selection.candidates.females[is.element(selection.candidates.females$id,next.gen),]
    obs_fert <- numeric(nrow(next.gen))
    next.gen <- cbind(next.gen, obs_fert)
     return (next.gen)
  }
############### No selection of males #############
nosel.males <- function (x) {  
    nfem <- ceiling((n.females/male.ratio))
    
    selection.candidates.males <-  subset( x,  sex  ==   1 ) 
    
    
    next.gen.males <- sample(selection.candidates.males$id, nfem)
    next.gen.males <- selection.candidates.males[is.element(selection.candidates.males$id,next.gen.males),]
    set( next.gen.males, j=which(colnames(next.gen.males) %in%
                                   "perm.env.ls"), value=NULL)
    
    return (next.gen.males)
  }
############### Index selection of old females ###############################################
# currently they are only truncated on their litter size phenotype
IndSelectionOldFemales <- function (x,y,z,solqual,year) { # x = next.gen, y = solutions, z = solutions.bw
   # delete last years index
   if("blup.fert" %in% colnames(x)) {
     set( x, j=which(colnames(x) %in% 
                              c("blup.fert", "sem.blup.fert"))  , value=NULL )
   }
  if("blup.bwnov" %in% colnames(x)) {
    set( x, j=which(colnames(x) %in% 
                             c("blup.bwnov", "sep.blup.bwnov"))  , value=NULL )
  }
  if("blup.qual" %in% colnames(x)) {
    set( x, j=which(colnames(x) %in% 
                      c("blup.qual", "sep.blup.qual"))  , value=NULL )
  }
  
  
   
   x <- merge(x, y, by ="id", all.x=TRUE) # merge solutions from blup to data.table containing the old females
   setkey(x, blup.fert)
   x <-merge(x, z, by="id", all.x=TRUE) # merge solutions from blup to data.table containing the old females
   x <- merge(x, solqual, by="id",all.x=TRUE)
   x[, `:=`(index.bw   = 100+ (blup.bwnov-mean(x$blup.bwnov))/(sqrt(var(x$blup.bwnov)))*10,
            index.fert = 100+ (blup.fert-mean(x$blup.fert))/(sqrt(var(x$blup.fert)))*10,
            index.qual = 100+ (blup.qual-mean(x$blup.qual))/(sqrt(var(x$blup.qual)))*10) ]
   x <- transform(x, comb.ind = index.bw*weight.bw.old.females+
                    index.fert*weight.fert.old.females+
                    weight.qual.old.females*index.qual)
   x <- subset(x, year -birthyear < max.age.females) 
   setkey(x, comb.ind)
      setorder(x, -comb.ind)
   #then I will need to make rangering
   old.females <- x[1:(n.females*prop.oldfemales),]
   set( old.females, j=which(colnames(old.females) %in% 
                               c("mating.will.1st.round",
                                 "mating.will.2nd.round"))  , value=NULL )
   if("dam.age" %in% colnames(old.females)) {
     set( old.females, j=which(colnames(old.females) %in% 
                                 c("dam.age"))  , value=NULL )
   }
   return(old.females)
 }  
 
############### Index selection for females ################
IndSelFemaleKits <- function (x,y,z,solqual,q) { # x = kit.list , y = solutions.fert, z = solutions.bw, q = old.females 
   x <- merge(x, y, by= "id", all.x=TRUE) # merge to solutions of blup of fertility
   x <- merge(x, z, by="id", all.x=TRUE)
   x <- merge(x, solqual, by="id",all.x=TRUE)
   x[, `:=`(index.bw   = 100+ (blup.bwnov-mean(x$blup.bwnov))/(sqrt(var(x$blup.bwnov)))*10,
            index.fert = 100+ (blup.fert-mean(x$blup.fert))/(sqrt(var(x$blup.fert)))*10,
            index.qual = 100+ (blup.qual-mean(x$blup.qual))/(sqrt(var(x$blup.qual)))*10) ]
   x <- transform(x, comb.ind = index.bw*weight.bw.kits+
                    index.fert*weight.fert.kits+
                    weight.qual.kits*index.qual)
   
   truncation.point <-  quantile( x$blup.fert,  probs =  quantile.setting ) 
   selection.candidates.females <- subset(x, blup.fert >= truncation.point) # throw away the smallest litters
   selection.candidates.females <-  subset( selection.candidates.females,  sex  ==   2) # take the female kits
   # selection.candidates.females <- subset( selection.candidates.females, phenotype.bw.oct < roof.body.size)
   # if (nrow(selection.candidates.females) < n.females*(1-prop.oldfemales)){ 
   #   truncation.point <-  quantile( x$own_littersize,  probs =  (quantile.setting  ) ) #can't change this since the kits won't have cards 
   #   selection.candidates.females <- subset(x, own_littersize >= truncation.point) # throw away the smallest litters
   #   selection.candidates.females <-  subset( selection.candidates.females,  sex  ==   2) # take the female kits
   #   selection.candidates.females <- subset( selection.candidates.females, phenotype.bw.oct < (roof.body.size+200)) #ease restrictions on size
   #   
   # }
   setkey(selection.candidates.females, comb.ind)           # in order to speed up ordering 
   setorder(selection.candidates.females,-comb.ind)         # order the kits according to body size 
   next.gen <- selection.candidates.females[1:(n.females-nrow(q)),]              # take the biggest kits
   setkey(next.gen, id)
   next.gen[next.gen,obs_fert:=0]  
   if("f0.dam" %in% colnames(next.gen)) {
     set( next.gen, j=which(colnames(next.gen) %in% 
                              "f0.dam")  , value=NULL )
   }
   if("prodyear" %in% colnames(next.gen)) {
     set( next.gen, j=which(colnames(next.gen) %in% 
                              "prodyear")  , value=NULL )
   }
   
   return (next.gen)
 }
############### Index selection for males ##############
 IndSelMaleKits <- function (x,y,z,solqual) { # x = kit.list, y = solutions, z = solutions.bw
   x <- merge(x, y, by= "id", all.x=TRUE) # merge to solutions of blup of fertility
   x <- merge(x, z, by="id", all.x=TRUE)
   x <- merge(x, solqual, by="id",all.x=TRUE)
   x[, `:=`(index.bw   = 100+ (blup.bwnov-mean(x$blup.bwnov))/(sqrt(var(x$blup.bwnov)))*10,
            index.fert = 100+ (blup.fert-mean(x$blup.fert))/(sqrt(var(x$blup.fert)))*10,
            index.qual = 100+ (blup.qual-mean(x$blup.qual))/(sqrt(var(x$blup.qual)))*10) ]
   x <- transform(x, comb.ind = index.bw*weight.bw.kits+
                    index.fert*weight.fert.kits+
                    weight.qual.kits*index.qual)
   truncation.point <-  quantile( x$blup.fert,  probs =  quantile.setting ) 
   selection.candidates.males <- subset(x, blup.fert >= truncation.point) # throw away the smallest litters
   selection.candidates.males <-  subset( selection.candidates.males,  sex  ==   1  ) # take the female kits
   setkey(selection.candidates.males, comb.ind)           # in order to speed up ordering 
   setorder(selection.candidates.males,-comb.ind)         # order the kits according to litter size 
   next.gen.males <- selection.candidates.males[1:n.males,]              # take the kits from biggest litters
   set( next.gen.males, j=which(colnames(next.gen.males) %in%
                                  "perm.env.ls"), value=NULL)
   if("f0.dam" %in% colnames(next.gen.males)) {
     set( next.gen.males, j=which(colnames(next.gen.males) %in% "f0.dam")  , value=NULL )
     
   }
   if("prodyear" %in% colnames(next.gen.males)) {
     set( next.gen.males, j=which(colnames(next.gen.males) %in% 
                              "prodyear")  , value=NULL )
   }
   # if("blup.fert" %in% colnames(next.gen.males)) {
   #   set( next.gen.males, j=which(colnames(next.gen.males) %in% 
   #                                  "blup.fert", "sep.blup.bwnov","index.bw","index.fert","blup.bwnov","comb.ind")  , value=NULL )
   # }
   
   
   return (next.gen.males)
   
 }
 
 
############### Make barren males ##################
PrepareMalesForMating <- function (x,year) { # x = next.gen.males 
  setkey(x, id)
  #   next.gen.males[next.gen.males,semen.quality:=rbinom( nrow(next.gen.males),  1,  male.inf )]
  #   next.gen.males[next.gen.males,mating.willingness:=rZIP( nrow(next.gen.males),  mu = male.ratio,  sigma = 0.05 )]
  x[,`:=`( semen.quality.1st=  rbinom( nrow(x),  1,  male.inf ),  
                        mating.willingness.1st = rZIP( nrow(x),  mu = male.ratio,  sigma = 0.05 ),
                        mating.willingness.2nd = rZIP( nrow(x), mu = 7.5, sigma = 0.05),
                        can.remate = rep(0, times=nrow(x)))]
  # truncation.point <- quantile(x$live.score, probs= 1/3 )
  x[,`:=`(semen.quality.2nd = semen.quality.1st)]
  if (use.comb.ind.for.males == 1 & year > 2 & selection.method ==blup) {
  setorder(x, -comb.ind)  
  } else if (use.comb.ind.for.males == 0 & weighing.method == sept){
    setorder(x,-live.score,-phenotype.bw.sept) 
  } else if (use.comb.ind.for.males == 0 & weighing.method == oct) {
    setorder(x,-live.score,-phenotype.bw.oct)
  } 
for (i in 1:ceiling((1-intensity.remating)*nrow(x))) {
  x$can.remate[i] <- 1
}
  x <- subset( x, mating.willingness.1st > 0 ) 
  return(x)
} 
############### mating will of females #############
PrepareFemalesForMating <- function (x,year) { # x = next.gen
  x[,`:=`(dam.age = ifelse(year-birthyear != 1, 0,1))] # 1 = older females, 0 = yearlings
  
  x[,`:=`( # this makes the mating will of the females in the first round of matings
    mating.will.1st.round = ifelse( dam.age == 1, rbinom(nrow(x), 1,mating.will.old.1st),
                                    ifelse(dam.age == 0, rbinom(nrow(x),1, mating.will.yearling.1st),0))
  )]
  
  x[,`:=`( # this makes the mating will of the females in the second round of matings
    mating.will.2nd.round = ifelse( dam.age == 1 & mating.will.1st.round == 1, rbinom(nrow(x), 1,mating.will.old.2nd),
                                    ifelse(dam.age == 0 & mating.will.1st.round ==1, rbinom(nrow(x),1, mating.will.yearling.2nd),0))
  )]
  return(x)
}
############### Update pedigree ########################
UpdatePedigree <- function (x,y,z,year) { # x = pedfile, y = next.gen, z = next.gen.males
if (year ==2) {
  x[,`:=`(true.sire = sire.id)]
  setnames(x, "sire.id", "sire.assumed")
}
  x <- rbind(x, # orginal pedfile
                   y[, .( id, sire.assumed, true.sire, dam.id, sex, birthyear ) ], #current gen females 
                   z[, .( id, sire.assumed, true.sire, dam.id, sex, birthyear )] ) # current gen males
  dups <- duplicated(x$id)
  x <- x[!dups,]
  return (x)
}
############### Breeding value of offspring in gen n###############################
# put in dam id's and make genetic value of each kit for fertility
# common cause for crash was a negative litter size
MakeKitsGenN <- function (x,y,z,year,p) { #x = mating.list, y = pedfile, z = big.pedfile
  kit.list <- x[rep(seq(nrow(x)), obs_fert), # makes the kit list by expanding the mating list
                c("dam.id",
                  "sire.id.1st",
                  "dam.fert",
                  "sire.fert.1st",
                  "sire.bw.oct.1st",
                  "sire.bw.sept.1st",
                  "sire.skin.length.1st",
                  "sire.skin.qual.1st",
                  "sire.live.qual.1st",
                  "f0.dam",
                  "obs_fert",
                  "dam.bw.oct",
                  "dam.bw.sept",
                  "dam.skin.length",
                  "dam.skin.qual",
                  "dam.live.qual",
                  "specific.env.bw",
                  "specific.env.bw",
                  "birthyear.dam",
                  "sire.id.2nd",
                  "sire.fert.2nd",
                  "sire.bw.oct.2nd",
                  "sire.bw.sept.2nd",
                  "sire.skin.length.2nd",
                  "sire.skin.qual.2nd",
                  "sire.live.qual.2nd"
                ) 
                , with=F] #specify which columns to incl.
  
  id <- seq(1:sum(x$obs_fert)) + max(z$id) # makes ID
  
  birthyear <- rep (year, sum(x$obs_fert)) # makes birthyear
  sex <- rbinom(sum(x$obs_fert),1,0.5)+1 # makes sex, TODO check if this is true
  true.sire <- numeric(nrow(kit.list))
  true.sire.check <- numeric(nrow(kit.list))
  mendelian <- as.data.table(rmvnorm(sum(x$obs_fert),sigma=sigma,method="svd"))
  colnames(mendelian) <- c("mend.bw.oct", "mend.bw.sept", "mend.litter.size", "mend.live.qual",
                           "mend.skin.qual", "mend.skin.length" )
  
  kit.list <-
    cbind (id,
           kit.list,
           birthyear,
           sex,
           mendelian,
           true.sire,
           true.sire.check)
  # now check the true sire before calculating inbreeding
  kit.list$sire.id.2nd  <-
    ifelse(kit.list$sire.id.2nd == 0,
           kit.list$sire.id.1st,
           kit.list$sire.id.2nd)
  # puts the 1st sire into the 2nd sire column for single mated females
  kit.list$sire.fert.2nd <-
    ifelse(
      kit.list$sire.id.2nd == kit.list$sire.id.1st,
      kit.list$sire.fert.1st,
      kit.list$sire.fert.2nd
    )
  kit.list$sire.bw.oct.2nd <-
    ifelse(
      kit.list$sire.id.2nd == kit.list$sire.id.1st,
      kit.list$sire.bw.oct.1st,
      kit.list$sire.bw.oct.2nd
    )
  kit.list$sire.bw.sept.2nd <-
    ifelse(
      kit.list$sire.id.2nd == kit.list$sire.id.1st,
      kit.list$sire.bw.sept.1st,
      kit.list$sire.bw.sept.2nd
    )
  
  kit.list$sire.skin.length.2nd <-
    ifelse(
      kit.list$sire.id.2nd == kit.list$sire.id.1st,
      kit.list$sire.skin.length.1st,
      kit.list$sire.skin.length.2nd
    )
  kit.list$sire.skin.qual.2nd <-
    ifelse(
      kit.list$sire.id.2nd == kit.list$sire.id.1st,
      kit.list$sire.skin.qual.1st,
      kit.list$sire.skin.qual.2nd
    )
  kit.list$sire.live.qual.2nd <-
    ifelse(
      kit.list$sire.id.2nd == kit.list$sire.id.1st,
      kit.list$sire.live.qual.1st,
      kit.list$sire.live.qual.2nd
    )
  
  kit.list$true.sire.check <-
    ifelse(kit.list$sire.id.1st != kit.list$sire.id.2nd,
           rbinom(nrow(kit.list), 1, 0.85), 1)
  # 85% chance that the kits are sired by 2nd mating
  kit.list$true.sire <-
    ifelse(kit.list$true.sire.check == 0,
           kit.list$sire.id.1st,
           kit.list$sire.id.2nd)
  # kit.list$true.sire <- ifelse( kit.list$sire.id.1st != 0 & kit.list$sire.id.2nd == 0, kit.list$sire.id.1st, kit.list$true.sire) # if the females are only mated once, the true sire is the first
  kit.list[, `:=`(
    true.sire.fert = ifelse(
      kit.list$true.sire == kit.list$sire.id.2nd,
      kit.list$sire.fert.2nd,
      kit.list$sire.fert.1st
    ),
    true.sire.bw.oct = ifelse(
      kit.list$true.sire == kit.list$sire.id.2nd,
      kit.list$sire.bw.oct.2nd,
      kit.list$sire.bw.oct.1st
    ),
    true.sire.bw.sept = ifelse(
      kit.list$true.sire == kit.list$sire.id.2nd,
      kit.list$sire.bw.sept.2nd,
      kit.list$sire.bw.sept.1st
    ),
    true.sire.skin.length = ifelse(
      kit.list$true.sire == kit.list$sire.id.2nd,
      kit.list$sire.skin.length.2nd,
      kit.list$sire.skin.length.1st
    ),
    true.sire.skin.qual = ifelse(
      kit.list$true.sire == kit.list$sire.id.2nd,
      kit.list$sire.skin.qual.2nd,
      kit.list$sire.skin.qual.1st
    ),
    true.sire.live.qual = ifelse(
      kit.list$true.sire == kit.list$sire.id.2nd,
      kit.list$sire.live.qual.2nd,
      kit.list$sire.live.qual.1st
    )
  )]
  
  setnames(kit.list, "sire.id.2nd", "sire.assumed")
  pedfile1 <- rbind(y, # orginal pedfile
                    kit.list[, .( id,sire.assumed, true.sire, dam.id, sex, birthyear ) ] ) #puts the kits into the pedfile
  
  ttt = trimPed(pedfile1,is.element( pedfile1$id, kit.list$id)) #trims the pedfile into only those related to kits
  t4 = cbind(pedfile1,ttt) #bind the logical vector to the temporary file
  if( make.obs.file == 1 & use.true.sire == 0) { 
    # This is done to get solutions for the litter size of the kits born this year
    pedigree <- file(description = c(paste("pedigree_",p, sep="")), open="w")
    
    write.table(pedfile1[,.(id,sire.assumed,dam.id,birthyear)], file= pedigree, col.names=FALSE, row.names=FALSE, quote=FALSE)
    close(pedigree)}
  if( make.obs.file == 1 & use.true.sire == 1) {
    pedigree <- file(description = c(paste("pedigree_",p, sep="")), open="w")
    
    write.table(pedfile1[,.(id,true.sire,dam.id,birthyear)], file= pedigree, col.names=FALSE, row.names=FALSE, quote=FALSE)
    close(pedigree)}
  
  
  
  t4 = subset(t4, ttt==TRUE) # delete everyone in the temporary pedigree who do not provide information on kits
  t4 <-as.data.frame(t4) # change temporary file into data frame since calcInbreeding expects a data frame
  f0 <- calcInbreeding(t4) # calculate inbreeding. NOTE This is what causes slowdown of program as the pedigree gets deeper
  t4 <- cbind(t4,f0) #cbind inbreeding to t4
  t4 <- as.data.table(t4)
  setkey(t4,birthyear)
  t4 <- subset(t4, birthyear ==year) # removes everything but the kits
  t4=t4[,.(id,f0)] # removes everything but the id and the inbreeding coefficiant
  kit.list=merge(kit.list, t4, by= "id") # merges the inbreeding coeff into the list
  
  ###### Now calculate the breeding values
  
  kit.list$mend.bw.oct <-  ifelse(kit.list$sex==1, kit.list$mend.bw.oct*sqrt(0.5*var.bw.oct.male),
                                  kit.list$mend.bw.oct*sqrt(0.5*var.bw.oct.female))
  kit.list$mend.bw.sept <-  ifelse(kit.list$sex==1, kit.list$mend.bw.sept*sqrt(0.5*var.bw.sept.male),
                                   kit.list$mend.bw.oct*sqrt(0.5*var.bw.sept.female))
  
  kit.list[, `:=`(
    litter.size = 0.5 * (dam.fert + true.sire.fert) +
      mend.litter.size * (sqrt(0.5 * variance.fertility) *
                            (1 - f0)) # Breeding value of offspring, littersize
    ,
    perm.env.ls = rnorm(sum(x$obs_fert)) * sqrt(var.perm.env.ls) # perm env for litter size
    ,
    bw.oct = 0.5 * (dam.bw.oct + true.sire.bw.oct) +
      mend.bw.oct * (1 - f0),
    bw.sept = 0.5 * (dam.bw.sept + true.sire.bw.sept) +
      mend.bw.oct * (1 - f0),
    live.qual = 0.5 * (dam.live.qual + true.sire.live.qual) + sqrt(0.5 *
         var.live.qual.gen) * (mend.live.qual) * (1 - f0),
    skin.qual = 0.5 * (dam.skin.qual + true.sire.skin.qual) + sqrt(0.5) *
      (mend.skin.qual) * (1 - f0),
    skin.length = 0.5 * (dam.skin.length + true.sire.skin.length) + sqrt(0.5) *
      (mend.skin.length) * (1 - f0)
  )]# Breeding value of offspring, body size
  
  setnames(kit.list, "obs_fert", "own_littersize") # changes obs fert into the littersize of the kit
  kit.list$dam.age <- ifelse( year - kit.list$birthyear.dam > 1, 1,0 )
  # generate phenotype for body size 
  # kit.list$phenotype.bw.oct <-
  #   ifelse(
  #     kit.list$sex == 1,
  #     MakePhenotypesBWMalesOct(
  #       mean.body.size.male.oct ,
  #       kit.list$bw.oct ,
  #       kit.list$specific.env.bw ,
  #       kit.list$own_littersize,
  #       kit.list$dam.age,
  #       x
  #     )
  #     ,
  #     MakePhenotypesBWFemalesOct(
  #       mean.body.size.female.oct ,
  #       kit.list$bw.oct ,
  #       kit.list$specific.env.bw,
  #       kit.list$own_littersize,
  #       x
  #     )
  #   )
  # # generate phenotype for body size
  # kit.list$phenotype.bw.sept <- ifelse(
  #   kit.list$sex == 1,
  #   MakePhenotypesBWMalesSept(
  #     mean.body.size.male.sept ,
  #     kit.list$bw.sept ,
  #     kit.list$specific.env.bw ,
  #     kit.list$own_littersize,
  #     kit.list$dam.age,
  #     x
  #   )
  #   ,
  #   MakePhenotypesBWFemalesOct(
  #     mean.body.size.female.sept ,
  #     kit.list$bw.sept ,
  #     kit.list$specific.env.bw,
  #     kit.list$own_littersize,
  #     x
  #   )
  # )
  
  kit.list[, `:=`(
    phenotype.bw.oct = ifelse(
      kit.list$sex == 1,
      MakePhenotypesBWMalesOct(
        mean.body.size.male.oct ,
        kit.list$bw.oct ,
        kit.list$specific.env.bw ,
        kit.list$own_littersize,
        kit.list$dam.age,
        x
      )
      ,
      MakePhenotypesBWFemalesOct(
        mean.body.size.female.oct ,
        kit.list$bw.oct ,
        kit.list$specific.env.bw,
        kit.list$own_littersize,
        x
      )
    ),
    phenotype.bw.sept = ifelse(
      kit.list$sex == 1,
      MakePhenotypesBWMalesSept(
        mean.body.size.male.sept ,
        kit.list$bw.sept ,
        kit.list$specific.env.bw ,
        kit.list$own_littersize,
        kit.list$dam.age,
        x
      )
      ,
      MakePhenotypesBWFemalesOct(
        mean.body.size.female.sept ,
        kit.list$bw.sept ,
        kit.list$specific.env.bw,
        kit.list$own_littersize,
        x
      )
    ),
    phenotype.live.qual = live.qual  + rnorm(nrow(kit.list)) *
      sqrt(var.live.qual.res),
    
    specific.env.bw = ifelse(
        sex == 1,
        rnorm(sum(x$obs_fert)) * sqrt(var.c.bw.sept.male),
        rnorm(sum(x$obs_fert)) *
          sqrt(var.c.bw.sept.female)
      ),
    specific.env.bw = ifelse(
      sex == 1,
      rnorm(sum(x$obs_fert)) * sqrt(var.c.bw.oct.male),
      rnorm(sum(x$obs_fert)) * sqrt(var.c.bw.oct.female)
    )
  )]
  
  # kit.list[,`:=`(specific.env.bw = ifelse(sex==1, rnorm(sum(x$obs_fert))*sqrt(var.c.bw.sept.male),
  #                                          rnorm(sum(x$obs_fert))*sqrt(var.c.bw.sept.female)),
  #                specific.env.bw = ifelse(sex==1, rnorm(sum(x$obs_fert))*sqrt(var.c.bw.oct.male),
  #                                         rnorm(sum(x$obs_fert))*sqrt(var.c.bw.oct.female)))] # generate specific env. for body size
  
  set( kit.list, j=which(colnames(kit.list) %in% c(
    "sire.fert.1st",
    "sire.bw.oct.1st",
    "sire.bw.sept.1st",
    "sire.skin.length.1st",
    "sire.skin.qual.1st",
    "sire.live.qual.1st",
    "sire.id.1st",
    "dam.fert",
    "dam.bw.oct",
    "dam.bw.sept",
    "dam.skin.length",
    "dam.skin.qual",
    "dam.live.qual",
    "specific.env.bw",
    "specific.env.bw",
    "sire.fert.2nd",
    "sire.bw.oct.2nd",
    "sire.bw.sept.2nd",
    "sire.skin.length.2nd",
    "sire.skin.qual.2nd",
    "sire.live.qual.2nd",
    "mend.bw.oct", 
    "mend.bw.sept", 
    "mend.litter.size", 
    "mend.live.qual",
    "mend.skin.qual", 
    "mend.skin.length",
    "true.sire.check",
    "true.sire.fert",
    "true.sire.bw.oct",
    "true.sire.bw.sept",
    "true.sire.skin.length",
    "true.sire.skin.qual",
    "true.sire.live.qual",
    "dam.age"
  )) , value=NULL )
  if (qual.classes == 5) {
    truncs <- qnorm(
      p = c(0.05, 0.3, 0.7, 0.95),
      mean = mean(kit.list$phenotype.live.qual),
                  sd = sqrt(var(
                    kit.list$phenotype.live.qual
                  )),
      lower.tail = TRUE
    )
    kit.list[, `:=`(live.score = ifelse(
      phenotype.live.qual >= truncs[4],
      5,
      ifelse(
        truncs[3] < phenotype.live.qual & phenotype.live.qual <= truncs[4],
        4,
        ifelse(
          phenotype.live.qual > truncs[2] & phenotype.live.qual <= truncs[3],
          3,
          ifelse(
            phenotype.live.qual > truncs[1] & phenotype.live.qual <= truncs[2],
            2,
            ifelse(phenotype.live.qual <=
                     truncs[1], 1, 0)
          )
        )
      )
    ))]
  } else if (qual.classes == 10) {
    truncs <-
      qnorm(
        p = c(0.01, 0.05, 0.15, 0.3, 0.5, 0.7, 0.85, 0.95, 0.99),
        mean = mean(kit.list$phenotype.live.qual),
                    sd = sqrt(var(
                      kit.list$phenotype.live.qual
                    )),
        lower.tail = TRUE
      )
    kit.list[, `:=`(live.score =
                      ifelse(
                        phenotype.live.qual >= truncs[9],
                        10,
                        ifelse(
                          truncs[8] < phenotype.live.qual & phenotype.live.qual <= truncs[9],
                          9,
                          ifelse(
                            phenotype.live.qual > truncs[7] & phenotype.live.qual <= truncs[8],
                            8,
                            ifelse(
                              phenotype.live.qual > truncs[6] & phenotype.live.qual <= truncs[7],
                              7,
                              ifelse(
                                phenotype.live.qual > truncs[5] & phenotype.live.qual <= truncs[6],
                                6,
                                ifelse(
                                  phenotype.live.qual > truncs[4] & phenotype.live.qual <= truncs[5],
                                  5,
                                  ifelse(
                                    phenotype.live.qual > truncs[3] & phenotype.live.qual <= truncs[4],
                                    4,
                                    ifelse(
                                      phenotype.live.qual > truncs[2] & phenotype.live.qual <= truncs[3],
                                      3,
                                      ifelse(
                                        phenotype.live.qual > truncs[1] & phenotype.live.qual <= truncs[2],
                                        2,
                                        ifelse(phenotype.live.qual <=
                                                 truncs[1], 1, 0)
                                      )
                                    )
                                  )
                                )
                              )
                            )
                          )
                        )
                      ))]
  }
  
  return(kit.list)
} 
MakeKitsGenN <-compiler::cmpfun(MakeKitsGenN,options= c(suppressAll=TRUE)) # performance boost
############### write observation file #########################
# This function writes to file in the manner that DMU requires the observed litter size
# in the year which it is performed and the replicate
WriteFertObservations <- function (x,year,p) { # x = mating.list
p <- p
# this is to format the mating list to write to the observation file for this replicate for the calculation of BV
  # fertility first and then write the observations for body weight of kits
  x[,`:=`(dam.age=0)]

    x[, c("dam.id","birthyear.dam")
            :=lapply(.SD, function(x) as.integer(x)), .SDcols=c("dam.id","birthyear.dam")]
  x[,`:=`(prodyear=as.integer(year), 
                    dam.age = ifelse( year-x$birthyear.dam == 1, as.integer(1), as.integer(2)))]
  setkey(x, dam.id)
x[, c("dam.id","dam.age")
            :=lapply(.SD, function(x) as.integer(x)), .SDcols=c("dam.id","dam.age")]
if (year ==1 ) {
  output <- file(description = paste("Replicate_",p, sep=""), open="w")
write.table(format(x[,.(dam.id,prodyear,dam.age,obs_fert)], nsmall=1), file= output, append= TRUE,col.names = FALSE, row.names = FALSE, quote = FALSE)
  close(con = output)
} else if (year > 1) {
  output <- file(description = paste("Replicate_",p, sep=""), open="a")
  write.table(format(x[,.(dam.id,prodyear,dam.age,obs_fert)], nsmall=1), file= output, append= TRUE,col.names = FALSE, row.names = FALSE, quote = FALSE)
  close(con = output)
}
} 
############### make dam.age checks #################
YearlingEffectOnFertility <- function (x,y){ # x = mating.list, y = year
  dam.age <- numeric(nrow(x))
  for (i in 1:nrow(x)) { # this is to change the parity into a two class thing
    if (y - x$birthyear.dam[i] == 1) {
      dam.age [i] <- yearling.effect
    }
    else {
      dam.age[i] <- 0
    }
  }
  x <- cbind(x, dam.age)
  return(x)
}

############  Write pedigree file to DMU
# Currently not in use, would be good idea to use this function instead of repeating the code within 
# MakeKitsGenN
 WritePedigreeForDMU <- function ( ) {
   if (make.obs.file == 1) {
     if (use.true.sire== 1) {

       write.table(pedfile[,.(id,true.sire,dam.id,birthyear)], file= pedigree, col.names=FALSE, row.names=FALSE, quote=FALSE)
       close(con=pedigree)
     }
   # set(pedfile, j = which(colnames(pedfile) %in%
   #                          "sex","true.sire"), value=NULL )
     write.table(pedfile[,.(id,sire.assumed,dam.id,birthyear)], file= pedigree, col.names=FALSE, row.names=FALSE, quote=FALSE)
   close(con=pedigree)
     }
 }
############### Make phenotype for body weight ######
 MakePhenotypesBWMalesOct <- function(x,y,z,t,u,mat) {
# x = mean, y = additive genetic, z = specific env, t = number of  sibs, mat =mating.list
   value <- x+ y + z*sqrt(var.c.bw.oct.male) + 
     rnorm( sum(mat$obs_fert))*(sqrt(var.res.bw.oct.male))+
     t*sib.effect.male+ u * bw.eff.damage
   return(value)
 } 
 MakePhenotypesBWFemalesOct <- function(x,y,z,t,mat) {
  # x = mean, y = additive genetic, z = specific env, t = number of sibs, mat=mating.list
   value <- x+ y + z*sqrt(var.c.bw.oct.female) + 
     rnorm( sum(mat$obs_fert))*(sqrt(var.res.bw.oct.female))+
     t*sib.effect.female 
   return(value)
 } 
 
 MakePhenotypesBWMalesSept <- function(x,y,z,t,u,mat) {
   # x = mean, y = additive genetic, z = specific env, t = number of  sibs, mat =mating.list
   value <- x+ y + z*sqrt(var.c.bw.sept.male) + 
     rnorm( sum(mat$obs_fert))*(sqrt(var.res.bw.oct.male))+
     t*sib.effect.male.sept+ u * bw.eff.damage
   return(value)
 } 
 MakePhenotypesBWFemalesSept <- function(x,y,z,t,mat) {
   # x = mean, y = additive genetic, z = specific env, t = number of sibs, mat=mating.list
   value <- x+ y + z*sqrt(var.c.bw.sept.female) + 
     rnorm( sum(mat$obs_fert))*(sqrt(var.res.bw.sept.female))+
     t*sib.effect.female.sept 
   return(value)
 } 

 MakePhenotypesBWMalesJune <- function(mean,a,c,sibs,damage,matinglist,maternal) {
   # mean = mean, a = additive genetic, c = specific env, sibs = number of  sibs, matinglist =mating.list
   value <- mean+ a + c*sqrt(var.bw.june.c.male) + 
     rnorm( sum(matinglist$obs_fert))*(sqrt(var.bw.june.res.male))+
     sibs*sib.effect.male.june+ damage * eff.dam.age.june +
     maternal*sqrt(var.maternal.male)
   return(value)
 } 
 MakePhenotypesBWFemalesJune <-  function(mean,a,c,sibs,damage,matinglist,maternal) {
   # x = mean, y = additive genetic, z = specific env, t = number of sibs, mat=mating.list
   value <- mean+ a + c*sqrt(var.bw.june.c.female) + 
     rnorm( sum(matinglist$obs_fert))*(sqrt(var.bw.june.res.female))+
     sibs*sib.effect.female.june + damage * eff.dam.age.june +
     maternal*sqrt(var.maternal.female)
   
   return(value)
 } 
 
 
 ############ Count sex of siblings ############
############### Function for counting sex of siblings ###########
# Currently unused
count.sex.siblings <- function (x) {
   temp <- subset(x, sex == 1)
   setkey(temp, dam.id)
   temp <- temp[, lapply(.SD, sum), by=dam.id, .SDcols=c("sex") ] 
 setnames(temp, "sex","sib.males")
 x <- merge(x, temp, by="dam.id", all.x=TRUE)
 x$sib.males[is.na(x$sib.males)] = 0 # change cases of NA into zeroes
 temp <- subset(x, sex ==2)
 setkey(temp, dam.id)
 temp <- temp[, lapply(.SD, function(p) sum(p)/2), by=dam.id, .SDcols=c("sex") ] 
 setnames(temp, "sex","sib.females")
 x <- merge(x, temp, by="dam.id", all.x=TRUE)
 x$sib.females[is.na(x$sib.females)] = 0 # change cases of NA into zeroes
 return(x)
 }
############### Calculate selection index ###########
CalculateBLUPLitterSize <- function () {
   system2("C:/Users/Notandi/Dropbox/Projects/simulation of mink farm/Output/DMU analysis/run_dmu4.bat", " bl_ass")
   # read the solutions and only keep the predictions of BV (they're not that right)
   if ( !file.exists('fort.70')) {
  now <-as.numeric(Sys.time()); howlong<-1; delt<-0; while(delt < howlong)
    { delt<-as.numeric(Sys.time())-now }
  }
  
   solutions <- as.matrix(read.table(file="bl_ass.SOL"))
   solutions <- as.data.table(solutions)
   solutions <- subset(solutions, V1 == 4 ) # throw away the estimation of permanent environment
   set (solutions, j=c("V1","V2","V3", "V4", "V6", "V7"), value= NULL)
   setnames(solutions, c("V5","V8", "V9"),c("id", "blup.fert", "sem.blup.fert"))
   
   
   # kit.list1 <- merge(kit.list, solutions, by= "id", all.x=TRUE)
   # next.gen <- merge(next.gen, solutions, by ="id", all.x=TRUE)
   
   # if (use.true.sire == 1) {
   #   system2("C:/Users/Notandi/Dropbox/Projects/simulation of mink farm/Output/DMU analysis/run_dmu4.bat", " bl_true")
   #   # read the solutions and only keep the predictions of BV (they're not that right)
   #   solutions <- as.matrix(read.table(file="bl_true.SOL"))
   #   solutions <- as.data.table(solutions)
   #   solutions <- subset(solutions, V1 == 4 ) # throw away the estimation of permanent environment
   #   set (solutions, j=c("V1","V2","V3", "V4", "V6", "V7"), value= NULL)
   #   setnames(solutions, c("V5","V8", "V9"),c("id", "blup.fert", "sem.blup.fert"))
   #   
   # }
   
   return(solutions)
 }
CalculateBLUPBodyWeightNov <- function () {
# run reml and then use the output as priors in BLUP
  # system2("C:/Users/Notandi/Dropbox/Projects/simulation of mink farm/Output/DMU analysis/run_dmuai.bat", " reml_bwnov")
  # run blup on BW  
system2("C:/Users/Notandi/Dropbox/Projects/simulation of mink farm/Output/DMU analysis/run_dmu4.bat", " bw_nov")
# read the solutions and only keep the predictions of BV (they're not that right)
   if ( !file.exists('fort.70')) {
  now <-as.numeric(Sys.time()); howlong<-1; delt<-0; while(delt < howlong)
    { delt<-as.numeric(Sys.time())-now }
  }
  solutions.bw <- as.matrix(read.table(file="C:/Users/Notandi/Dropbox/Projects/simulation of mink farm/Output/DMU analysis/bw_nov.SOL"))
  # }
solutions.bw <- as.data.table(solutions.bw)
solutions.bw <- subset(solutions.bw, V1 == 4 ) # throw away the estimation of permanent environment
set (solutions.bw, j=c("V1","V2","V3", "V4", "V6", "V7"), value= NULL)
setnames(solutions.bw, c("V5","V8", "V9"),c("id", "blup.bwnov", "sep.blup.bwnov"))
return (solutions.bw)
}
CalculateBLUPQuality <- function () {
  # run reml and then use the output as priors in BLUP
  # system2("C:/Users/Notandi/Dropbox/Projects/simulation of mink farm/Output/DMU analysis/run_dmuai.bat", " reml_bwnov")
  # run blup on BW  
  system2("C:/Users/Notandi/Dropbox/Projects/simulation of mink farm/Output/DMU analysis/run_dmu4.bat", " qual")
  # read the solutions and only keep the predictions of BV (they're not that right)
   if ( !file.exists('fort.70')) {
  now <-as.numeric(Sys.time()); howlong<-1; delt<-0; while(delt < howlong)
    { delt<-as.numeric(Sys.time())-now }
  }
  solutions.qual <- as.matrix(read.table(file="C:/Users/Notandi/Dropbox/Projects/simulation of mink farm/Output/DMU analysis/qual.SOL"))
  # }
  solutions.qual <- as.data.table(solutions.qual)
  solutions.qual <- subset(solutions.qual, V1 == 4, select= c("V5","V8", "V9")) # throw away the estimation of permanent environment
  # set (solutions.bw, j=c("V1","V2","V3", "V4", "V6", "V7"), value= NULL)
  setnames(solutions.qual, c("V5","V8", "V9"),c("id", "blup.qual", "sep.blup.qual"))
  return (solutions.qual)
}

############### Modify Dir file for fertility ################
 # this function changes the .DIR file for fertility and BW
# So it takes in the replicate from the current replicate of the simulation
 ModifyDIRFile <- function (p) {
   # fertility modification
   if (trace.ped == 0 ){
   dirfile <- readLines("bl_ass.DIR")
   dirfile[8] <- c(paste("$DATA  ASCII (3,1,-9999) Replicate_",p, sep=""))
   # change the input file for BLUP so it uses the next outputfile

      writeLines(dirfile, "bl_ass.DIR")
   dirfile[25] <- c(paste("$VAR_STR 1 PED 2 ASCII pedigree_",p, sep="")) 
   # change the input file for BLUP so it uses the next pedigree
   writeLines(dirfile,"bl_ass.DIR")
   
   # bw modification
   dirfile <- readLines("bw_nov.DIR")
   dirfile[8] <- c(paste("$DATA  ASCII (5,2,-9999) Phenotypes",p, sep="")) 
   # change the input file for BLUP so it uses the next outputfile
   writeLines(dirfile, "bw_nov.DIR")
   dirfile[25] <- c(paste("$VAR_STR 2 PED 2 ASCII Big_pedigree_",p, sep="")) 
   # change the input file for BLUP so it uses the next pedigree
   writeLines(dirfile,"bw_nov.DIR")
   # quality
   dirfile <- readLines("qual.DIR")
   dirfile[8] <- c(paste("$DATA  ASCII (5,2,-9999) Phenotypes",p, sep="")) 
   # change the input file for BLUP so it uses the next outputfile
   writeLines(dirfile, "qual.DIR")
   dirfile[25] <- c(paste("$VAR_STR 1 PED 2 ASCII Big_pedigree_",p, sep="")) 
   # change the input file for BLUP so it uses the next pedigree
   writeLines(dirfile,"qual.DIR")
   } else if (trace.ped == 1) {
   dirfile <- readLines("trace.DIR")
   dirfile[1] <- c(paste("Big_pedigree_",p, sep=""))
   writeLines(dirfile,"trace.DIR")
   # change the pedigree name in the dir file
   dirfile <- readLines("bl_ass.DIR")
   dirfile[25] <- c(paste("$VAR_STR 1 PED 2 ASCII trace.PRUNE")) # change the input file for BLUP so it uses the next pedigree
    writeLines(dirfile,"bl_ass.DIR")
  # change the pedigree name in case of trace being asked for
    dirfile <- readLines("bw_nov.DIR")
    dirfile[25] <- c(paste("$VAR_STR 1 PED 2 ASCII trace.PRUNE")) # change the input file for BLUP so it uses the next pedigree
    writeLines(dirfile,"bw_nov.DIR")
  # dir file for prune
    dirfile <- readLines("qual.DIR")
    dirfile[25] <- c(paste("$VAR_STR 1 PED 2 ASCII trace.PRUNE")) # change the input file for BLUP so it uses the next pedigree
    writeLines(dirfile,"qual.DIR")
    
   }
   
    }
 
############### Create observation file for phenotypes #################
# currently only size
WriteObservationFileBodyWeight <- function (x,year,p,solutions) {
  inter <- as.integer(rep(1, times=nrow(x))) # intercept for the quality regression
  x<- cbind(x, inter)
  x[, c("id","dam.id","own_littersize","sex")
              :=lapply(.SD, function(x) as.integer(x)), .SDcols=c("id","dam.id","own_littersize","sex")]
  # x[,`:=`(prodyear=as.integer(ifelse( sex == 1, year*(sex+6),year*(sex+8))))]
  # NOTE if simulation runs to 81 years, then year*sex will be the same for year 63 for dams and 81 for males
  if(mask.phenotypes == 1 & year > 1) {
    x <- merge(x, solutions, by= "id", all.x=TRUE) # merge to solutions of blup of fertility
    truncation.point <-  quantile( x$blup.fert,  probs =  quantile.setting ) 
    x <- subset(x, blup.fert >= truncation.point, 
                select= c("id",
                          "dam.id",
                          "own_littersize",
                          "sex",
                          "phenotype.bw.oct",
                          "phenotype.bw.sept",
                          "live.score",
                          "inter"))
  }
  if (year  == 1 ) {phenotypes <- file(description = paste("Phenotypes",p, sep=""), open="w")
  if (weighing.method == oct){
    x[, c("id","dam.id","own_littersize","sex","inter")
      :=lapply(.SD, function(x) as.integer(x)), .SDcols=c("id","dam.id","own_littersize","sex","inter")]
    
    write.table(format(x[,.(id,dam.id,sex,own_littersize,inter,phenotype.bw.oct,live.score)], nsmall=1, digits=2), 
                file= phenotypes, append= TRUE,col.names = FALSE, row.names = FALSE, quote = FALSE)
  } else if (weighing.method == sept ) {
    x[, c("id","dam.id","own_littersize","sex","inter")
      :=lapply(.SD, function(x) as.integer(x)), .SDcols=c("id","dam.id","own_littersize","sex","inter")]
    
    write.table(format(x[,.(id,dam.id,sex,own_littersize,inter,phenotype.bw.sept,live.score)], nsmall=1, digits=2), 
                file= phenotypes, append= TRUE,col.names = FALSE, row.names = FALSE, quote = FALSE)
  }
  close(con=phenotypes)
   } else if (year > 1) {
     x[, c("id","dam.id","own_littersize","sex","inter")
       :=lapply(.SD, function(x) as.integer(x)), .SDcols=c("id","dam.id","own_littersize","sex","inter")]
     
     phenotypes <- file(description = paste("Phenotypes",p, sep=""), open="a")
     if (weighing.method == oct){
      write.table(format(x[,.(id,dam.id,sex,own_littersize,inter,phenotype.bw.oct,live.score)], nsmall=1, digits=2), 
                  file= phenotypes, append= TRUE,col.names = FALSE, row.names = FALSE, quote = FALSE)
     } else if (weighing.method == sept ) {
       write.table(format(x[,.(id,dam.id,sex,own_littersize,inter,phenotype.bw.sept,live.score)], nsmall=1, digits=2), 
                   file= phenotypes, append= TRUE,col.names = FALSE, row.names = FALSE, quote = FALSE)
}
  close(con=phenotypes)
  
   }
  # if (year == 4) {
  #   quality <- file(description = "quality", open="w")
  #   write.table(format(x[,.(id,dam.id,sex,own_littersize,inter,phenotype.bw.oct,live.score)], nsmall=1, digits=2), 
  #               file= quality,append= TRUE ,col.names = FALSE, row.names = FALSE, quote = FALSE)
  #   close(quality)
  #   }
  }

# WriteObservationFileQuality <- function (kit.list, year, p, solutions) {
#   x[, c("id")
#     :=lapply(.SD, function(x) as.integer(x)), .SDcols=c("id")]
#   
# }
############### Make kit.list big.pedfile ##################################
WriteBigPedigree <- function (x,y,year,p) { # x = kit.list or kit.list , y = pedfile, year = year
  if (year == 1){
    big.pedfile <- copy(y)
    big.pedfile[,`:=`(true.sire = sire.id)]
    setnames(big.pedfile, "sire.id", "sire.assumed")
    big.pedfile <- rbindlist( list (big.pedfile, x[,.(id,true.sire,sire.assumed,dam.id,sex,birthyear)]), use.names=TRUE)
  } else if (year > 1 ) {
      big.pedfile <- rbindlist( list (y, x[,.(id,true.sire,sire.assumed,dam.id,sex,birthyear)]), use.names=TRUE)
  }
  big.pedigree.con <- file(description=paste("big_pedigree_",p, sep=""),open="w")
  if (use.true.sire == 1 ) {
  write.table(big.pedfile[,.(id,true.sire,dam.id,birthyear)], file= big.pedigree.con, col.names=FALSE, row.names=FALSE, quote=FALSE)
  } else if (use.true.sire==0) {
    write.table(big.pedfile[,.(id,sire.assumed,dam.id,birthyear)], file= big.pedigree.con, col.names=FALSE, row.names=FALSE, quote=FALSE)
  }
  close(con=big.pedigree.con)
  return(big.pedfile)
}
############### Update big pedigree ############
update.big.pedigree <-
  function (x, y, z) {
    # x = big.pedfile, y = next.gen, z = next.gen.males
    big.pedfile <- rbindlist(list (x,
                                   y[, .(id, true.sire, sire.assumed, dam.id, sex, birthyear)],
                                   z[, .(id, true.sire, sire.assumed, dam.id, sex, birthyear)]), use.names =
                               TRUE)
    dups <- duplicated(big.pedfile$id)
    big.pedfile <- big.pedfile[!dups, ]
    return(big.pedfile)
  }
############## Trace pedigree ###############
# this is a function to trim the pedigree which gets too big and unwieldy after
# a certain number of generations
# it takes a list of animals to trace and returns a trimmed pedigree that is 
# then used by DMU to calculate breeding values
# currently unused, controlled by switch since it did not improve functionality

TracePed <- function(list.trace, next.gen) {
  trace <- file(description = "trace", open="w")
  write.table(rbind(list.trace[,.(id)],next.gen[,.(id)]), 
              file= trace,
              col.names=FALSE, row.names=FALSE, quote=FALSE)
  close(trace)
  system2("run_DmuTrace.bat", "trace")
}
