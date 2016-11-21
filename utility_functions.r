# This file has all the subroutines used in the simulation, note that some
# functions are currently not used in this build of the simulation

############### Creation of Gen0 females ############
# This function creates the base population of females
GenerateBaseFemales <- function ( leg2,
                                  t,
                                  n.females,
                                 mating.will.yearling.1st,
                                 mating.will.yearling.2nd,
                                 qual.classes ) 
  {
  id        <-  seq(1:n.females)
  add.gen <- rmvnorm(n.females, sigma = G_sigma,method="svd" ) 
  # multiplies the random deviates by the variances
  add.gen <- as.data.table(t(t(add.gen) * sqrt(variances))) 
  
   colnames(add.gen)     <-
    c(
      "live.qual",
      "h.length",
      "skin.qual",
      "skin.length.male",
      "skin.length.female",
      "litter.size",
      "bw_f",
      "bw_m",
      "rfi1_m",
      "rfi2_m",
      "rfi1_f",
      "rfi2_f"
    )
  perm.env.ls           <- rnorm(n.females) * sqrt(var.perm.env.ls)
  sex                   <-  rep(2, times = n.females)
  dam.id                <-  rep(0, times = n.females)
  sire.id               <-  rep(0, times = n.females)
  birthyear             <-  rep(0, times = n.females)
  mating.will.1st.round <-
    rbinom(n.females, 1, mating.will.yearling.1st)
  mating.will.2nd.round <- numeric(n.females)

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

    gen0.females[,`:=`(phenotype.bw.oct =(BW.mean.females +
                   bw_f + rnorm(nrow(gen0.females))*sqrt(pe.var.bw.female)+
                     rnorm(nrow(gen0.females))*sqrt(bw.res.female))  ,
                    phenotype.live.qual = live.qual + 
                      rnorm(nrow(gen0.females))*sqrt(var.live.qual.res),
                    phenotype.skin.length = mean.skin.length.female + skin.length.female + 
                      rnorm(nrow(gen0.females))*sqrt(var.skin.length.c.female)+
                    rnorm(nrow(gen0.females))*sqrt(var.skin.length.res.female),
                    phenotype.skin.qual = skin.qual + rnorm(nrow(gen0.females))*
                      sqrt(var.skin.qual.res),
                    phenotype.h.length = h.length+ rnorm(nrow(gen0.females))*
                      sqrt(var.h.length.res)
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
GenerateBaseMales <- function (leg2,
                               t,
                               n.males,
                               male.ratio,
                               male.inf,
                               qual.classes,
                               intensity.remating,
                               n.females) 
  {
  mating.willingness.1st <-  numeric( n.males )  
  mating.willingness.2nd  <-  numeric( n.males)
  semen.quality.1st      <-  numeric( n.males )  
  semen.quality.2nd      <-  numeric( n.males )  
  id                 <-  numeric( n.males )  
  add.gen <- rmvnorm(n.males, sigma = G_sigma,method="svd" ) 
  # multiplies the random deviates by the variances
  add.gen <- as.data.table(t(t(add.gen) * sqrt(variances))) 
  
  colnames(add.gen)     <-
    c(
      "live.qual",
      "h.length",
      "skin.qual",
      "skin.length.male",
      "skin.length.female",
      "litter.size",
      "bw_f",
      "bw_m",
      "rfi1_m",
      "rfi2_m",
      "rfi1_f",
      "rfi2_f"
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

  gen0.males[,`:=`(phenotype.bw.oct = 
                     (BW.mean.males +
                        bw_f + rnorm(nrow(gen0.males))*sqrt(pe.var.bw.male)+
                        rnorm(nrow(gen0.males))*sqrt(bw.res.male)),
                     phenotype.live.qual = live.qual + 
                       rnorm(nrow(gen0.males))*sqrt(var.live.qual.res),
                   phenotype.skin.length = mean.skin.length.male + skin.length.male + 
                     rnorm(nrow(gen0.males))*sqrt(var.skin.length.c.male)+
                     rnorm(nrow(gen0.males))*sqrt(var.skin.length.res.male),
                   phenotype.skin.qual = skin.qual + rnorm(nrow(gen0.males))*
                     sqrt(var.skin.qual.res),
                   phenotype.h.length = h.length+ rnorm(nrow(gen0.males))*
                     sqrt(var.h.length.res)
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

mate <- function (x, 
                  y, 
                  year,
                  use.blup.to.assort.mat,
                  selection.method,
                  crossmating,
                  purebreeding,
                  pr.barren.double.mating.yearling,
                  pr.barren.double.mating.old,
                  pr.barren.one.mating.yearling,
                  pr.barren.one.mating.old ) 
{
  # x = males, y = females
  # browser()
  if (year >2 & use.blup.to.assort.mat == 1 & selection.method ==blup) {
    setorder(x, -comb.ind)
    setorder(y, -comb.ind)
  }
  if (year <= 2) {
    x$comb.ind <- 0
  } else if( selection.method == phenotypic) {
    x$comb.ind <- 0
    y[sample(nrow(y)),]
    
  }
  mating.list <-
    x[rep(seq(nrow(x)), mating.willingness.1st),  #expands the male list into a mating list, based on mat.will.1st
      c(
        "id",
        "litter.size",
        "bw_f",
        "bw_m",
        "rfi1_m",
        "rfi2_m",
        "rfi1_f",
        "rfi2_f",
        "live.qual",
        "skin.qual",
        "skin.length.male",
        "skin.length.female",
        "h.length",
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
      "litter.size",
      "bw_f",
      "bw_m",
      "rfi1_m",
      "rfi2_m",
      "rfi1_f",
      "rfi2_f",
      "live.qual",
      "skin.qual",
      "skin.length.male",
      "skin.length.female",
      "h.length",
      "live.score"
    ),
    c(
      "sire.id.1st",
      "sire.fert.1st",
      "sire.bw_f.1st",
      "sire.bw_m.1st",
      "sire.rfi1_m.1st",
      "sire.rfi2_m.1st",
      "sire.rfi1_f.1st",
      "sire.rfi2_f.1st",
      "sire.live.qual.1st",
      "sire.skin.qual.1st",
      "sire.skin.length.male.1st",
      "sire.skin.length.female.1st",
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
    "dam.bw_f",
    "dam.bw_m",
    "dam.rfi1_m",
    "dam.rfi2_m",
    "dam.rfi1_f",
    "dam.rfi2_f",
    "dam.skin.qual",
    "dam.skin.length.male",
    "dam.skin.length.female",
    "dam.live.qual",
    "perm.env.ls",
    "birthyear.dam",
    "sire.id.2nd",
    "sire.bw_f.2nd",
    "sire.bw_m.2nd",
    "sire.rfi1_m.2nd",
    "sire.rfi2_m.2nd",
    "sire.rfi1_f.2nd",
    "sire.rfi2_f.2nd",
    "sire.fert.2nd",
    "sire.skin.qual.2nd",
    "sire.skin.length.male.2nd",
    "sire.skin.length.female.2nd",
    "sire.live.qual.2nd",
    "sire.h.length.2nd",
    "sire.live.score.2nd"
  ) := 0]
  # moved scaling of litter specific environment to the phenotype function later on
  perm.env.bw.f <- rnorm(nrow(mating.list) )*sqrt(pe.var.bw.male)
  specific.env.skin <- rnorm(nrow(mating.list))
  perm.env.bw.m <- rnorm(nrow(mating.list) )*sqrt(pe.var.bw.male)
  perm.env.rfi <- rmvnorm(nrow(mating.list), sigma = P_RFI,method="svd" )
  colnames(perm.env.rfi) <- c("pe1.rfi.m", "pe2.rfi.m","pe1.rfi.f", "pe2.rfi.f")
  perm.env.rfi <- as.data.table(t(t(perm.env.rfi) * sqrt(pe.var.rfi))) 
  
  mating.list <- cbind(mating.list, perm.env.bw.f,perm.env.bw.m,specific.env.skin,perm.env.rfi)
  # perm.env.bw <- t(t(perm.env.bw)*pe.var.bw.female)
  # here I subset the dam list to throw away those who will not mate on first round
  y <- subset (y, mating.will.1st.round == 1)
  mating.list$dam.id       <- y[1:nrow(mating.list), .(id)]
  mating.list$dam.fert     <- y[1:nrow(mating.list), .(litter.size)]
  mating.list$dam.bw_f    <- y[1:nrow(mating.list), .(bw_f)]
  mating.list$dam.bw_m    <- y[1:nrow(mating.list), .(bw_m)]
  mating.list$dam.rfi1_m   <- y[1:nrow(mating.list), .(rfi1_m)]
  mating.list$dam.rfi2_m   <- y[1:nrow(mating.list), .(rfi2_m)]
  mating.list$dam.rfi1_f   <- y[1:nrow(mating.list), .(rfi1_f)]
  mating.list$dam.rfi2_f   <- y[1:nrow(mating.list), .(rfi2_f)]
  mating.list$perm.env.ls  <- y[1:nrow(mating.list), .(perm.env.ls)]
  mating.list$dam.skin.length.male  <- y[1:nrow(mating.list), .(skin.length.male)]
  mating.list$dam.skin.length.female  <- y[1:nrow(mating.list), .(skin.length.female)]
  mating.list$dam.skin.qual  <- y[1:nrow(mating.list), .(skin.qual)]
  mating.list$dam.live.qual  <- y[1:nrow(mating.list), .(live.qual)]
  mating.list$dam.live.score <- y[1:nrow(mating.list), .(live.score)]
  mating.list$dam.h.length   <- y[1:nrow(mating.list), .(h.length)]
  
  
  
  
  if ("f0" %in% colnames(y)) {
    mating.list$f0  <- y[1:nrow(mating.list), .(f0)]
  }
  mating.list$birthyear.dam  <- y[1:nrow(mating.list), .(birthyear)]
  setnames(mating.list, "f0", "f0.dam")
  if (year == 1) {
    setnames(mating.list, "f0.dam", "f0") 
    # this is because of in the breeding values of the first gen this is needed, silly really since its zero all over
  }
  # now to make the 2nd round of matings
  
  if (crossmating == 1 ) {
    mating.list <- mating.list[sample(nrow(mating.list)),] # randomize order of females
    myvars <- c("dam.id",
                "sire.id.2nd",
                "semen.quality.2nd",
                "sire.fert.2nd",
                "sire.bw_f.2nd",
                "sire.bw_m.2nd",
                "sire.rfi1_m.2nd",
                "sire.rfi2_m.2nd",
                "sire.rfi1_f.2nd",
                "sire.rfi2_f.2nd",
                "sire.skin.length.male.2nd",
                "sire.skin.length.female.2nd",
                "sire.skin.qual.2nd",
                "sire.live.qual.2nd",
                "sire.h.length.2nd")
    mating.list.temp <-
      as.matrix(mating.list[,myvars,with=FALSE]) 
    myvars <- c("sire.id.2nd",
                "semen.quality.2nd",
                "sire.fert.2nd",
                "sire.bw_f.2nd",
                "sire.bw_m.2nd",
                "sire.rfi1_m.2nd",
                "sire.rfi2_m.2nd",
                "sire.rfi1_f.2nd",
                "sire.rfi2_f.2nd",
                "sire.skin.length.male.2nd",
                "sire.skin.length.female.2nd",
                "sire.skin.qual.2nd",
                "sire.live.qual.2nd",
                "sire.h.length.2nd")
    mating.list <- mating.list[,!myvars,with=FALSE]
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
                "bw_f",# 4
                "bw_m",# 5
                "rfi1_m",# 6
                "rfi2_m",# 7
                "rfi1_f", #8
                "rfi2_f", # 9
                "skin.length.male", # 10
                "skin.length.female", # 11
                "skin.qual", # 12
                "live.qual", # 13
                "live.score", # 14
                "h.length", # 15
                "mating.willingness.2nd") #16 
    setkey(x, id)
    x <- subset(x, can.remate == 1 & mating.willingness.2nd > 0)
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
        if (s < nrow(mating.list.temp)) {
          if ( i == 1) {
            # this is not a solution, need to make another controlling mechanism if the mating willingness
            #exceeds the number of females to be mated
            mating.list.temp[[s - (j - 1), 2]] <- x[[i, 1]]          # id of male
            mating.list.temp[[s - (j - 1), 3]] <- x[[i, 2]]          # semen.quality
            mating.list.temp[[s - (j - 1), 4]] <- x[[i, 3]]          # fertility of male
            mating.list.temp[[s - (j - 1), 5]] <- x[[i, 4]]          # bw_f
            mating.list.temp[[s - (j - 1), 6]] <- x[[i, 5]]         # bw_m
            mating.list.temp[[s - (j - 1), 7]] <- x[[i, 6]]        # rfi1_m
            mating.list.temp[[s - (j - 1), 8]] <- x[[i, 7]]        # rfi2_m
            mating.list.temp[[s - (j - 1), 9]] <- x[[i, 8]]        # rfi1_f
            mating.list.temp[[s - (j - 1), 10]] <- x[[i, 9]]        # rfi2_f
            mating.list.temp[[s - (j - 1), 11]] <- x[[i, 10]]        # skin.length of male
            mating.list.temp[[s - (j - 1), 12]] <- x[[i, 11]]        # skin.length of female
            mating.list.temp[[s - (j - 1), 13]] <- x[[i, 12]]        # skin.qual  of male
            mating.list.temp[[s - (j - 1), 14]] <- x[[i, 13]]        # live.qual  of male
            mating.list.temp[[s - (j - 1), 15]] <- x[[i, 14]]        # h.length  of male
          } else if ( i > 1) {
            t <- sum(x[1:(i-1), 16])+1
            mating.list.temp[[t+(j-1), 2]] <- x[[i, 1]]          # id of male
            mating.list.temp[[t+(j-1), 3]] <- x[[i, 2]]          # semen.quality
            mating.list.temp[[t+(j-1), 4]] <- x[[i, 3]]          # fertility of male
            mating.list.temp[[t+(j-1), 5]] <- x[[i, 4]]          # bw_f
            mating.list.temp[[t+(j-1), 6]] <- x[[i, 5]]         # bw_m
            mating.list.temp[[t+(j-1), 7]] <- x[[i, 6]]        # rfi1_m
            mating.list.temp[[t+(j-1), 8]] <- x[[i, 7]]        # rfi2_m
            mating.list.temp[[t+(j-1), 9]] <- x[[i, 8]]        # rfi1_f
            mating.list.temp[[t+(j-1), 10]] <- x[[i, 9]]        # rfi2_f
            mating.list.temp[[t+(j-1), 11]] <- x[[i, 10]]        # skin.length of male
            mating.list.temp[[t+(j-1), 12]] <- x[[i, 11]]        # skin.length of female
            mating.list.temp[[t+(j-1), 13]] <- x[[i, 12]]        # skin.qual  of male
            mating.list.temp[[t+(j-1), 14]] <- x[[i, 13]]        # live.qual  of male
            mating.list.temp[[t+(j-1), 15]] <- x[[i, 14]]        # h.length  of male
          }
        }
        if (s > nrow(mating.list.temp)) {
          while (sum(x[1:(i - 1) , 16]) + j <= nrow(mating.list.temp)) {
            # Here i solve the problem of if the males have more mating willingness than the number of females
            if (s < nrow(mating.list.temp)) {
              if ( i == 1) {
                # this is not a solution, need to make another controlling mechanism if the mating willingness
                #exceeds the number of females to be mated
                mating.list.temp[[s - (j - 1), 2]] <- x[[i, 1]]          # id of male
                mating.list.temp[[s - (j - 1), 3]] <- x[[i, 2]]          # semen.quality
                mating.list.temp[[s - (j - 1), 4]] <- x[[i, 3]]          # fertility of male
                mating.list.temp[[s - (j - 1), 5]] <- x[[i, 4]]          # bw_f
                mating.list.temp[[s - (j - 1), 6]] <- x[[i, 5]]         # bw_m
                mating.list.temp[[s - (j - 1), 7]] <- x[[i, 6]]        # rfi1_m
                mating.list.temp[[s - (j - 1), 8]] <- x[[i, 7]]        # rfi2_m
                mating.list.temp[[s - (j - 1), 9]] <- x[[i, 8]]        # rfi1_f
                mating.list.temp[[s - (j - 1), 10]] <- x[[i, 9]]        # rfi2_f
                mating.list.temp[[s - (j - 1), 11]] <- x[[i, 10]]        # skin.length of male
                mating.list.temp[[s - (j - 1), 12]] <- x[[i, 11]]        # skin.length of female
                mating.list.temp[[s - (j - 1), 13]] <- x[[i, 12]]        # skin.qual  of male
                mating.list.temp[[s - (j - 1), 14]] <- x[[i, 13]]        # live.qual  of male
                mating.list.temp[[s - (j - 1), 15]] <- x[[i, 14]]        # h.length  of male
              } else if ( i > 1) {
                t <- sum(x[1:(i-1), 16])+1
                mating.list.temp[[s - (j - 1), 2]] <- x[[i, 1]]          # id of male
                mating.list.temp[[s - (j - 1), 3]] <- x[[i, 2]]          # semen.quality
                mating.list.temp[[s - (j - 1), 4]] <- x[[i, 3]]          # fertility of male
                mating.list.temp[[s - (j - 1), 5]] <- x[[i, 4]]          # bw_f
                mating.list.temp[[s - (j - 1), 6]] <- x[[i, 5]]         # bw_m
                mating.list.temp[[s - (j - 1), 7]] <- x[[i, 6]]        # rfi1_m
                mating.list.temp[[s - (j - 1), 8]] <- x[[i, 7]]        # rfi2_m
                mating.list.temp[[s - (j - 1), 9]] <- x[[i, 8]]        # rfi1_f
                mating.list.temp[[s - (j - 1), 10]] <- x[[i, 9]]        # rfi2_f
                mating.list.temp[[s - (j - 1), 11]] <- x[[i, 10]]        # skin.length of male
                mating.list.temp[[s - (j - 1), 12]] <- x[[i, 11]]        # skin.length of female
                mating.list.temp[[s - (j - 1), 13]] <- x[[i, 12]]        # skin.qual  of male
                mating.list.temp[[s - (j - 1), 14]] <- x[[i, 13]]        # live.qual  of male
                mating.list.temp[[s - (j - 1), 15]] <- x[[i, 14]]        # h.length  of male
              }
            }
            break
          }
        }
      }
    } 
    mating.list.temp <- as.data.table(mating.list.temp)
    mating.list <- merge(mating.list, mating.list.temp, by="dam.id")
    
  } else if (crossmating == 0 ) {
    if (purebreeding == 1 ) { # if purebreeding we will not use diff males in 2nd mating
      # we will not attempt to mate females with different male if 1st fails
      mating.list.allowed <- mating.list 
    } else if (purebreeding == 0){ 
      mating.list.allowed <- subset(mating.list, can.remate == 1)
    }
    
    mating.list.allowed[, `:=`(IDX = 1:.N) , by = sire.id.1st]
    mating.list.allowed$sire.id.2nd       <-
      ifelse(
        mating.list.allowed$IDX <= mating.list.allowed$mating.willingness.2nd,
        mating.list.allowed$sire.id.1st   ,
        0
      )
    mating.list.allowed$test <- ifelse(mating.list.allowed$sire.id.1st == mating.list.allowed$sire.id.2nd,TRUE,FALSE)
    ############### 2nd mating for animals allowed to remate with 1st ###########
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
    mating.list.allowed$sire.bw_f.2nd       <-
      ifelse(
        mating.list.allowed$test == TRUE,
        mating.list.allowed$sire.bw_f.1st      ,
        0
      )
    mating.list.allowed$sire.bw_m.2nd       <-
      ifelse(
        mating.list.allowed$test == TRUE,
        mating.list.allowed$sire.bw_m.1st      ,
        0
      )
    mating.list.allowed$sire.rfi1_m.2nd       <-
      ifelse(
        mating.list.allowed$test == TRUE,
        mating.list.allowed$sire.rfi1_m.1st      ,
        0
      )
    mating.list.allowed$sire.rfi2_m.2nd       <-
      ifelse(
        mating.list.allowed$test == TRUE,
        mating.list.allowed$sire.rfi2_m.1st      ,
        0
      )
    mating.list.allowed$sire.rfi1_f.2nd       <-
      ifelse(
        mating.list.allowed$test == TRUE,
        mating.list.allowed$sire.rfi1_f.1st      ,
        0
      )
    mating.list.allowed$sire.rfi2_f.2nd       <-
      ifelse(
        mating.list.allowed$test == TRUE,
        mating.list.allowed$sire.rfi2_f.1st      ,
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
    if (purebreeding == 0 ) {  

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
      # browser()
      # here I should reorder the dams that are leftovers to prioritize younger dams and the ones mated with shitty males
      setkey(mating.list.leftover, dam.id)
      if (year >2 & use.blup.to.assort.mat == 1 & selection.method ==blup) {
        # setorder(mating.list.leftover, -comb.ind)
      } else if (use.blup.to.assort.mat == 0) {
        # setorder(mating.list.leftover, -birthyear.dam, -dam.live.score)
      }
      myvars <- c("dam.id",
                  "sire.id.2nd",
                  "semen.quality.2nd",
                  "sire.fert.2nd",
                  "sire.bw_f.2nd",
                  "sire.bw_m.2nd",
                  "sire.rfi1_m.2nd",
                  "sire.rfi2_m.2nd",
                  "sire.rfi1_f.2nd",
                  "sire.rfi2_f.2nd",
                  "sire.skin.length.male.2nd",
                  "sire.skin.length.female.2nd",
                  "sire.skin.qual.2nd",
                  "sire.live.qual.2nd",
                  "sire.h.length.2nd")
      mating.list.leftover.temp <-
        as.matrix(mating.list.leftover[,myvars,with=FALSE]) 
      myvars <- c("sire.id.2nd",
                  "semen.quality.2nd",
                  "sire.fert.2nd",
                  "sire.bw_f.2nd",
                  "sire.bw_m.2nd",
                  "sire.rfi1_m.2nd",
                  "sire.rfi2_m.2nd",
                  "sire.rfi1_f.2nd",
                  "sire.rfi2_f.2nd",
                  "sire.skin.length.male.2nd",
                  "sire.skin.length.female.2nd",
                  "sire.skin.qual.2nd",
                  "sire.live.qual.2nd",
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
                  "bw_f",# 4
                  "bw_m",# 5
                  "rfi1_m",# 6
                  "rfi2_m",# 7
                  "rfi1_f", #8
                  "rfi2_f", # 9
                  "skin.length.male", # 10
                  "skin.length.female", # 11
                  "skin.qual", # 12
                  "live.qual", # 13
                  "live.score", # 14
                  "h.length", # 15
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
              mating.list.leftover.temp[[s - (j - 1), 5]] <- x[[i, 4]]          # bw_f
              mating.list.leftover.temp[[s - (j - 1), 6]] <- x[[i, 5]]         # bw_m
              mating.list.leftover.temp[[s - (j - 1), 7]] <- x[[i, 6]]        # rfi1_m
              mating.list.leftover.temp[[s - (j - 1), 8]] <- x[[i, 7]]        # rfi2_m
              mating.list.leftover.temp[[s - (j - 1), 9]] <- x[[i, 8]]        # rfi1_f
              mating.list.leftover.temp[[s - (j - 1), 10]] <- x[[i, 9]]        # rfi2_f
              mating.list.leftover.temp[[s - (j - 1), 11]] <- x[[i, 10]]        # skin.length of male
              mating.list.leftover.temp[[s - (j - 1), 12]] <- x[[i, 11]]        # skin.length of female
              mating.list.leftover.temp[[s - (j - 1), 13]] <- x[[i, 12]]        # skin.qual  of male
              mating.list.leftover.temp[[s - (j - 1), 14]] <- x[[i, 13]]        # live.qual  of male
              mating.list.leftover.temp[[s - (j - 1), 15]] <- x[[i, 14]]        # h.length  of male
            } else if ( i > 1) {
              t <- sum(x[1:(i-1), 16])+1
              mating.list.leftover.temp[[t+(j-1), 2]] <- x[[i, 1]]          # id of male
              mating.list.leftover.temp[[t+(j-1), 3]] <- x[[i, 2]]          # semen.quality
              mating.list.leftover.temp[[t+(j-1), 4]] <- x[[i, 3]]          # fertility of male
              mating.list.leftover.temp[[t+(j-1), 5]] <- x[[i, 4]]          # bw_f
              mating.list.leftover.temp[[t+(j-1), 6]] <- x[[i, 5]]         # bw_m
              mating.list.leftover.temp[[t+(j-1), 7]] <- x[[i, 6]]        # rfi1_m
              mating.list.leftover.temp[[t+(j-1), 8]] <- x[[i, 7]]        # rfi2_m
              mating.list.leftover.temp[[t+(j-1), 9]] <- x[[i, 8]]        # rfi1_f
              mating.list.leftover.temp[[t+(j-1), 10]] <- x[[i, 9]]        # rfi2_f
              mating.list.leftover.temp[[t+(j-1), 11]] <- x[[i, 10]]        # skin.length of male
              mating.list.leftover.temp[[t+(j-1), 12]] <- x[[i, 11]]        # skin.length of female
              mating.list.leftover.temp[[t+(j-1), 13]] <- x[[i, 12]]        # skin.qual  of male
              mating.list.leftover.temp[[t+(j-1), 14]] <- x[[i, 13]]        # live.qual  of male
              mating.list.leftover.temp[[t+(j-1), 15]] <- x[[i, 14]]        # h.length  of male
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
                  mating.list.leftover.temp[[s - (j - 1), 5]] <- x[[i, 4]]          # bw_f
                  mating.list.leftover.temp[[s - (j - 1), 6]] <- x[[i, 5]]           # bw_m
                  mating.list.leftover.temp[[s - (j - 1), 7]] <- x[[i, 6]]          # rfi1_m
                  mating.list.leftover.temp[[s - (j - 1), 8]] <- x[[i, 7]]          # rfi2_m
                  mating.list.leftover.temp[[s - (j - 1), 9]] <- x[[i, 8]]          # rfi1_f
                  mating.list.leftover.temp[[s - (j - 1), 10]] <- x[[i, 9]]         # rfi2_f
                  mating.list.leftover.temp[[s - (j - 1), 11]] <- x[[i, 10]]        # skin.length of male
                  mating.list.leftover.temp[[s - (j - 1), 12]] <- x[[i, 11]]        # skin.length of female
                  mating.list.leftover.temp[[s - (j - 1), 13]] <- x[[i, 12]]        # skin.qual  of male
                  mating.list.leftover.temp[[s - (j - 1), 14]] <- x[[i, 13]]        # live.qual  of male
                  mating.list.leftover.temp[[s - (j - 1), 15]] <- x[[i, 14]]        # h.length  of male
                } else if ( i > 1) {
                  t <- sum(x[1:(i-1), 16])+1
                  mating.list.leftover.temp[[t+(j-1), 2]] <- x[[i, 1]]          # id of male
                  mating.list.leftover.temp[[t+(j-1), 3]] <- x[[i, 2]]          # semen.quality
                  mating.list.leftover.temp[[t+(j-1), 4]] <- x[[i, 3]]          # fertility of male
                  mating.list.leftover.temp[[t+(j-1), 5]] <- x[[i, 4]]          # bw_f
                  mating.list.leftover.temp[[t+(j-1), 6]] <- x[[i, 5]]          # bw_m
                  mating.list.leftover.temp[[t+(j-1), 7]] <- x[[i, 6]]          # rfi1_m
                  mating.list.leftover.temp[[t+(j-1), 8]] <- x[[i, 7]]          # rfi2_m
                  mating.list.leftover.temp[[t+(j-1), 9]] <- x[[i, 8]]          # rfi1_f
                  mating.list.leftover.temp[[t+(j-1), 10]] <- x[[i, 9]]         # rfi2_f
                  mating.list.leftover.temp[[t+(j-1), 11]] <- x[[i, 10]]        # skin.length of male
                  mating.list.leftover.temp[[t+(j-1), 12]] <- x[[i, 11]]        # skin.length of female
                  mating.list.leftover.temp[[t+(j-1), 13]] <- x[[i, 12]]        # skin.qual  of male
                  mating.list.leftover.temp[[t+(j-1), 14]] <- x[[i, 13]]        # live.qual  of male
                  mating.list.leftover.temp[[t+(j-1), 15]] <- x[[i, 14]]        # h.length  of male
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
    } # end of no crossmating if
  } #end of purebreeding if
  if (purebreeding == 1 ) { 
    mating.list <- mating.list.allowed 
    set(
      mating.list,
      j = c(
        "IDX",
        "test"
      ),
      value = NULL
    )
    
  } 
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
  specific.env.skin <- rnorm(nrow(mating.list))
  cbind(mating.list, specific.env.skin)
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
MakeKitsGen0 <- function (x,y, z, leg2,qual.classes,true.sire.chance) { #x = mating.list, y= effgen0.males, z = year, leg2 = legendre poly
  kit.list <- x[rep(seq(nrow(x)), obs_fert), 
                c("dam.id",
                  "sire.id.1st",
                  "dam.fert",
                  "sire.fert.1st",
                  "sire.bw_f.1st",
                  "sire.bw_m.1st",
                  "sire.rfi1_m.1st",
                  "sire.rfi2_m.1st",
                  "sire.rfi1_f.1st",
                  "sire.rfi2_f.1st",
                  "sire.skin.length.male.1st",
                  "sire.skin.length.female.1st",
                  "sire.skin.qual.1st",
                  "sire.live.qual.1st",
                  "sire.h.length.1st",
                  "f0",
                  "obs_fert",
                  "dam.bw_f",
                  "dam.bw_m",
                  "dam.rfi1_m",
                  "dam.rfi2_m",
                  "dam.rfi1_f",
                  "dam.rfi2_f",
                  "dam.skin.length.male",
                  "dam.skin.length.female",
                  "dam.skin.qual",
                  "dam.live.qual",
                  "dam.h.length",
                  "perm.env.bw.f",
                  "perm.env.bw.m",
                  "birthyear.dam",
                  "sire.id.2nd",
                  "sire.fert.2nd",
                  "sire.bw_f.2nd",
                  "sire.bw_m.2nd",
                  "sire.rfi1_m.2nd",
                  "sire.rfi2_m.2nd",
                  "sire.rfi1_f.2nd",
                  "sire.rfi2_f.2nd",
                  "sire.skin.length.male.2nd",
                  "sire.skin.length.female.2nd",
                  "sire.skin.qual.2nd",
                  "sire.live.qual.2nd",
                  "sire.h.length.2nd",
                  "specific.env.skin",
                  "pe1.rfi.m", 
                  "pe2.rfi.m",
                  "pe1.rfi.f", 
                  "pe2.rfi.f"
                ) 
                , with=F] #specify which columns to incl.
  id <- seq(1:sum(x$obs_fert)) + max(y$id) # makes ID
  birthyear <- rep (1, sum(x$obs_fert)) # makes birthyear
  sex <- rbinom(sum(x$obs_fert),1,0.5)+1 # makes sex, TODO check if this is true
  true.sire <- numeric(nrow(kit.list))
  true.sire.check <- numeric(nrow(kit.list))
  mendelian <- rmvnorm(sum(x$obs_fert),sigma=G_sigma,method="svd")
  mendelian <- as.data.table(t(t(mendelian) * sqrt(0.5*variances))) 
  
  colnames(mendelian)     <-
    c(
      "mend.live.qual",
      "mend.h.length",
      "mend.skin.qual",
      "mend.skin.length.male",
      "mend.skin.length.female",
      "mend.litter.size",
      "mend.bw_f",
      "mend.bw_m",
      "mend.rfi1_m",
      "mend.rfi2_m",
      "mend.rfi1_f",
      "mend.rfi2_f"    )
  
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
  kit.list$sire.bw_f.2nd <-
    ifelse(
      kit.list$true.sire.check == TRUE,
      kit.list$sire.bw_f.1st,
      kit.list$sire.bw_f.2nd
    )
  kit.list$sire.bw_m.2nd <-
    ifelse(
      kit.list$true.sire.check == TRUE,
      kit.list$sire.bw_m.1st,
      kit.list$sire.bw_m.2nd
    )
  kit.list$sire.rfi1_m.2nd <-
    ifelse(
      kit.list$true.sire.check == TRUE,
      kit.list$sire.rfi1_m.1st,
      kit.list$sire.rfi1_m.2nd
    )
  kit.list$sire.rfi2_m.2nd <-
    ifelse(
      kit.list$true.sire.check == TRUE,
      kit.list$sire.rfi2_m.1st,
      kit.list$sire.rfi2_m.2nd
    )
  kit.list$sire.rfi1_f.2nd <-
    ifelse(
      kit.list$true.sire.check == TRUE,
      kit.list$sire.rfi1_f.1st,
      kit.list$sire.rfi1_f.2nd
    )
  kit.list$sire.rfi2_f.2nd <-
    ifelse(
      kit.list$true.sire.check == TRUE,
      kit.list$sire.rfi2_f.1st,
      kit.list$sire.rfi2_f.2nd
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
  kit.list$true.sire.check <- ifelse( kit.list$sire.id.1st != kit.list$sire.id.2nd, rbinom(nrow(kit.list), 1, true.sire.chance), 1) # 85% chance that the kits are sired by 2nd mating
  kit.list$true.sire <- ifelse( kit.list$true.sire.check == 0, kit.list$sire.id.1st, kit.list$sire.id.2nd)
  kit.list[, `:=`(true.sire.fert = 
                    ifelse(kit.list$true.sire == kit.list$sire.id.2nd, 
                           kit.list$sire.fert.2nd, kit.list$sire.fert.1st), 
                  true.sire.bw_m = 
                    ifelse(kit.list$true.sire == kit.list$sire.id.2nd, 
                           kit.list$sire.bw_m.2nd, kit.list$sire.bw_m.1st),
                  true.sire.bw_f = 
                    ifelse(kit.list$true.sire == kit.list$sire.id.2nd, 
                           kit.list$sire.bw_f.2nd, kit.list$sire.bw_f.1st),
                  true.sire.rfi1_m = 
                    ifelse(kit.list$true.sire == kit.list$sire.id.2nd, 
                           kit.list$sire.rfi1_m.2nd, kit.list$sire.rfi1_m.1st),
                  true.sire.rfi2_m = 
                    ifelse(kit.list$true.sire == kit.list$sire.id.2nd, 
                           kit.list$sire.rfi2_m.2nd, kit.list$sire.rfi2_m.1st),
                  true.sire.rfi1_f = 
                    ifelse(kit.list$true.sire == kit.list$sire.id.2nd, 
                           kit.list$sire.rfi1_f.2nd, kit.list$sire.rfi1_f.1st),
                  true.sire.rfi2_f = 
                    ifelse(kit.list$true.sire == kit.list$sire.id.2nd, 
                           kit.list$sire.rfi2_f.2nd, kit.list$sire.rfi2_f.1st),              
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
                           kit.list$sire.h.length.2nd, kit.list$sire.h.length.1st)
  )]
  
  
  kit.list[, `:=`(litter.size = 0.5*(dam.fert + true.sire.fert) + 
                    mend.litter.size,  # Breeding value of offspring, littersize
                  bw_f = 0.5*(true.sire.bw_f + dam.bw_f) + mend.bw_f,
                  bw_m = 0.5*(true.sire.bw_m + dam.bw_m) + mend.bw_m,
                  rfi1_m = 0.5*(true.sire.rfi1_m + dam.rfi1_m) + mend.rfi1_m,
                  rfi2_m = 0.5*(true.sire.rfi2_m + dam.rfi2_m) + mend.rfi2_m,
                  rfi1_f = 0.5*(true.sire.rfi1_f + dam.rfi1_f) + mend.rfi1_f,
                  rfi2_f = 0.5*(true.sire.rfi2_f + dam.rfi2_f) + mend.rfi2_f,
                  perm.env.ls = rnorm(sum(x$obs_fert))*sqrt(var.perm.env.ls), # perm env for litter size
                  skin.length.male = 0.5*(dam.skin.length.male + true.sire.skin.length.male) + 
                    mend.skin.length.male,
                  skin.length.female = 0.5*(dam.skin.length.female + true.sire.skin.length.female) + 
                    mend.skin.length.female,
                  live.qual = 0.5*(dam.live.qual + true.sire.live.qual)+ (mend.live.qual),
                  skin.qual = 0.5*(dam.skin.qual + true.sire.skin.qual)+ (mend.skin.qual),
                  h.length = 0.5*(dam.h.length + true.sire.h.length)+ (mend.h.length)
  )]# Breeding value of offspring, body size
  q <- as.matrix(as.data.frame(polynomial.values(polynomials = leg2, x =t[6])))
  
  
  setnames(kit.list, c("obs_fert","sire.id.2nd"), c("own_littersize","sire.assumed")) # renames obs_fert to own littersize of kits
  kit.list$dam.age <- ifelse( z - kit.list$birthyear.dam > 1, 1,0 )
  
  ## Make phenotypes for body weight
  kit.list$phenotype.bw.oct <- ifelse(
    kit.list$sex == 1,
    MakePhenotypesBWMalesOct(BW.mean.males , kit.list$bw_m , kit.list$perm.env.bw.m , kit.list$own_littersize, kit.list$dam.age,x  )
    , MakePhenotypesBWFemalesOct(BW.mean.females , kit.list$bw_f , kit.list$perm.env.bw.f, kit.list$own_littersize,x ))
  # kit.list$phenotype.bw.sept <- ifelse( kit.list$sex == 1,
  #                                   MakePhenotypesBWMalesSept(mean.body.size.male.sept , kit.list$bw.sept , kit.list$specific.env.bw , kit.list$own_littersize, kit.list$dam.age,x  )
  #                                   , MakePhenotypesBWFemalesOct(mean.body.size.female.sept , kit.list$bw.sept , kit.list$specific.env.bw,kit.list$own_littersize,x ))
  # kit.list$phenotype.live.qual <- kit.list$live.qual  + rnorm(nrow(kit.list))*sqrt(var.live.qual.res)
  kit.list[, `:=`( 
    phenotype.live.qual = live.qual  + rnorm(nrow(kit.list)) *
      sqrt(var.live.qual.res),
    phenotype.skin.length = ifelse(
      sex == 1,
      mean.skin.length.male + skin.length.male + rnorm(1)*sqrt(var.skin.length.res.male)+
        specific.env.skin*sqrt(var.skin.length.c.male),
      mean.skin.length.female + skin.length.female + rnorm(1)*sqrt(var.skin.length.res.female)+
        specific.env.skin*sqrt(var.skin.length.c.female) ),
    phenotype.skin.qual = 
      skin.qual + rnorm(nrow(kit.list))*sqrt(var.skin.qual.res),
    phenotype.h.length = h.length+ rnorm(nrow(kit.list))*
      sqrt(var.h.length.res)
  )
  ]  
  
  set( kit.list, j=which(colnames(kit.list) %in% c(
    "sire.fert.1st",
    "sire.bw_f.1st",
    "sire.bw_m.1st",
    "sire.rfi1_m.1st",
    "sire.rfi2_m.1st",
    "sire.rfi1_f.1st",
    "sire.rfi2_f.1st",
    "sire.skin.length.male.1st",
    "sire.skin.length.female.1st",
    "sire.skin.qual.1st",
    "sire.live.qual.1st",
    "sire.bw.june.1st",
    "sire.bw.june.maternal.1st",
    "sire.h.length.1st",
    "sire.id.1st",
    "dam.fert",
    "dam.bw_f",
    "dam.bw_m",
    "dam.rfi1_m",
    "dam.rfi2_m",
    "dam.rfi1_f",
    "dam.rfi2_f",
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
    "sire.bw_f.2nd",
    "sire.bw_m.2nd",
    "sire.rfi1_m.2nd",
    "sire.rfi2_m.2nd",
    "sire.rfi1_f.2nd",
    "sire.rfi2_f.2nd",
    "sire.skin.length.male.2nd",
    "sire.skin.length.female.2nd",
    "sire.skin.qual.2nd",
    "sire.live.qual.2nd",
    "sire.bw.june.2nd",
    "sire.bw.june.maternal.2nd",
    "sire.h.length.2nd",
    "mend.bw_f",
    "mend.bw_m",
    "mend.rfi1_m",
    "mend.rfi2_m",
    "mend.rfi1_f",
    "mend.rfi2_f",
    "mend.litter.size", 
    "mend.live.qual",
    "mend.skin.qual", 
    "mend.skin.length.male",
    "mend.skin.length.female",
    "true.sire.check",
    "true.sire.fert",
    "true.sire.bw_f",
    "true.sire.bw_m",
    "true.sire.rfi1_m",
    "true.sire.rfi2_m",
    "true.sire.rfi1_f",
    "true.sire.rfi2_f",
    "true.sire.skin.length.male",
    "true.sire.skin.length.female",
    "true.sire.skin.qual",
    "true.sire.live.qual",
    "true.sire.h.length",
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
PhenoSelectionOldFemales <- function (y,x, year,
                                      max.age.females,
                                      n.females,
                                      prop.oldfemales ) 
  { # y = gen0.females, x = mating.list
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
PhenoSelectionFemaleKits <- function (x, # x = kit.list
                                      y, # y = old females
                                      quantile.setting.ls,
                                      quantile.setting.bw,n.females) {   
  truncation.point <-  quantile( x$own_littersize,  probs =  quantile.setting.ls ) 
  selection.candidates.females <- subset(x, own_littersize >= truncation.point) # throw away the smallest litters
  selection.candidates.females <-  subset( selection.candidates.females,  sex  ==   2) # take the female kits

  truncation.point <-  quantile( selection.candidates.females$phenotype.bw.oct,  probs =  quantile.setting.bw ) 
selection.candidates.females <- subset(selection.candidates.females, phenotype.bw.oct >= truncation.point) # throw away the smallest litters
for (i in 1:4){
if (nrow(selection.candidates.females) <= (n.females - nrow(y)) ) {
  truncation.point <-  quantile( x$own_littersize,  probs =  quantile.setting.ls ) 
  selection.candidates.females <- subset(x, own_littersize >= truncation.point) # throw away the smallest litters
  selection.candidates.females <-  subset( selection.candidates.females,  sex  ==   2) # take the female kits
  
  truncation.point <-  quantile( selection.candidates.females$phenotype.bw.oct,  probs =  (quantile.setting.bw -i/10) ) 
  selection.candidates.females <- subset(selection.candidates.females, phenotype.bw.oct >= truncation.point) # throw away the smallest litters
} else if (nrow(selection.candidates.females) >= nrow(y) ) {
  break}
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
############### Selection of yearling males ###############################
PhenoSelectionMaleKits <- function (x,quantile.setting.ls,quantile.setting.bw,n.males ) { # x = kit.list 
  # browser()
  truncation.point <-  quantile( x$own_littersize,  probs =  quantile.setting.ls ) 
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
# nosel.old.females <- function (x) {
#     old.females <- x[is.element(x$id, sample(x$id, n.females*prop.oldfemales)),]
#     #safety valve to remove too old females
#     # here i cbind the observation for litter size in this year to the old.females list
#     old.females <- merge(old.females, mating.list, by.x="id", by.y="dam.id", all.x=TRUE)
#     set( old.females, j=which(colnames(old.females) %in% 
#                                 c("sire.id.y","barren","dam.fert","sire.fert")) 
#          , value=NULL )
#      setnames(old.females, "sire.id.x", "sire.id")
#      old.females[is.na(obs_fert), obs_fert:= 0]   
#      old.females[is.na(f0), f0:= 0]   
#         old.females[old.females,own_littersize:=0]  # adds a column named own littersize with value 0 
#         return (old.females)
#   }
# ############### No selection of yearlings ###################
# nosel.yearlings.females <- function (x) {
#     selection.candidates.females <-  subset( x,  sex  ==   2  ) 
#     # kit.file_df in generation after first
#     # We choose females to fit n.females, based on prop of old females
#     
#     next.gen <- sample(selection.candidates.females$id, (n.females-nrow(old.females)))
#     next.gen <- selection.candidates.females[is.element(selection.candidates.females$id,next.gen),]
#     obs_fert <- numeric(nrow(next.gen))
#     next.gen <- cbind(next.gen, obs_fert)
#      return (next.gen)
#   }
# ############### No selection of males #############
# nosel.males <- function (x) {  
#     nfem <- ceiling((n.females/male.ratio))
#     
#     selection.candidates.males <-  subset( x,  sex  ==   1 ) 
#     
#     
#     next.gen.males <- sample(selection.candidates.males$id, nfem)
#     next.gen.males <- selection.candidates.males[is.element(selection.candidates.males$id,next.gen.males),]
#     set( next.gen.males, j=which(colnames(next.gen.males) %in%
#                                    "perm.env.ls"), value=NULL)
#     
#     return (next.gen.males)
#   }
############### Index selection of old females ###############################################
# currently they are only truncated on their litter size phenotype

IndSelectionOldFemales <- function (x,
                                    y,
                                    year,
                                    weight.bw.old.females,
                                    weight.fert.old.females,
                                    weight.qual.old.females,n.females,prop.oldfemales) {
  # x = next.gen, y = solutions, z = solutions.bw
  # delete last years index
  if("blup.fert" %in% colnames(x)) {
    set( x, j=which(colnames(x) %in% 
                      c("blup.fert"))  , value=NULL )
  }
  if("blup.bwnov" %in% colnames(x)) {
    set( x, j=which(colnames(x) %in% 
                      c("blup.bwnov"))  , value=NULL )
  }
  if("blup.qual" %in% colnames(x)) {
    set( x, j=which(colnames(x) %in% 
                      c("blup.qual"))  , value=NULL )
  }
  x <- merge(x, y, by ="id", all.x=TRUE) # merge solutions from blup to data.table containing the old females
  setkey(x, blup.fert)
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
IndSelFemaleKits <-
  function (x,
            y,
            q,
            mblup,
            weight.bw.kits,
            weight.fert.kits,
            weight.qual.kits,n.females) {
    # x = kit.list , y = solutions q = old.females
    x <-
      merge(x, y, by = "id", all.x = TRUE) # merge to solutions of blup of fertility
    x[, `:=`(
      index.bw   = 100 + (blup.bwnov - mean(x$blup.bwnov)) / (sqrt(var(x$blup.bwnov))) *
        10,
      index.fert = 100 + (blup.fert - mean(x$blup.fert)) / (sqrt(var(x$blup.fert))) *
        10,
      index.qual = 100 + (blup.qual - mean(x$blup.qual)) / (sqrt(var(x$blup.qual))) *
        10
    )]
    x <- transform(
      x,
      comb.ind = index.bw * weight.bw.kits +
        index.fert * weight.fert.kits +
        weight.qual.kits * index.qual
    )
    if (mblup == 0) { # throw away the smallest litters
  truncation.point <-  quantile( x$blup.fert,  probs =  quantile.setting.ls ) 
  selection.candidates.females <- subset(x, blup.fert >= truncation.point)
  } else if (mblup == 1) { selection.candidates.females <- x 
  } 
  selection.candidates.females <-  subset( selection.candidates.females,  sex  ==   2) # take the female kits
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
IndSelMaleKits <- function (x, # x = kit.list
                            y, # y = solutions
                            weight.bw.kits,
                            weight.fert.kits,
                            weight.qual.kits,n.males) {

  x <- merge(x, y, by= "id", all.x=TRUE) # merge to solutions 
  x[, `:=`(index.bw   = 100+ (blup.bwnov-mean(x$blup.bwnov))/(sqrt(var(x$blup.bwnov)))*10,
           index.fert = 100+ (blup.fert-mean(x$blup.fert))/(sqrt(var(x$blup.fert)))*10,
           index.qual = 100+ (blup.qual-mean(x$blup.qual))/(sqrt(var(x$blup.qual)))*10) ]
  x <- transform(x, comb.ind = index.bw*weight.bw.kits+
                   index.fert*weight.fert.kits+
                   weight.qual.kits*index.qual)
  truncation.point <-  quantile( x$blup.fert,  probs =  quantile.setting.ls ) 
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
  
  return (next.gen.males)
  
}

 
############### Make barren males ##################
PrepareMalesForMating <-
  function (x,
            year,
            use.comb.ind.for.males,
            selection.method,
            weighing.method,
            intensity.remating,
            pseudo.import,
            pseudo.import.prop,
            weight.bw.kits,
            weight.fert.kits,
            weight.qual.kits) {
    # x = next.gen.males
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
  if (pseudo.import == 1) {
  x <- MaskPseudoImport (pseudo.import,
                         pseudo.import.prop,
                         x,
                         selection.method,
                         weight.bw.kits,
                         weight.fert.kits,
                         weight.qual.kits)
  }
  return(x)
} 
############### mating will of females #############
PrepareFemalesForMating <-
  function (x,
            year,
            mating.will.old.1st,
            mating.will.yearling.1st,
            mating.will.old.2nd,
            mating.will.yearling.2nd) {
    # x = next.gen
    x[, `:=`(dam.age = ifelse(year - birthyear != 1, 0, 1))] # 1 = older females, 0 = yearlings
    
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
MakeKitsGenN <- function (x,
                          y,
                          z,
                          year,
                          p,
                          leg2,
                          t,
                          true.sire.chance,
                          use.true.sire,
                          make.obs.file,
                          qual.classes
) { #x = mating.list, y = pedfile, z = big.pedfile
  kit.list <- x[rep(seq(nrow(x)), obs_fert), # makes the kit list by expanding the mating list
                c("dam.id",
                  "sire.id.1st",
                  "dam.fert",
                  "sire.fert.1st",
                  "sire.bw1_f.1st",
                  "sire.bw2_f.1st",
                  "sire.bw3_f.1st",
                  "sire.bw1_m.1st",
                  "sire.bw2_m.1st",
                  "sire.bw3_m.1st",
                  "sire.rfi1_m.1st",
                  "sire.rfi2_m.1st",
                  "sire.rfi1_f.1st",
                  "sire.rfi2_f.1st",
                  "sire.skin.length.male.1st",
                  "sire.skin.length.female.1st",
                  "sire.skin.qual.1st",
                  "sire.live.qual.1st",
                  "sire.h.length.1st",
                  "obs_fert",
                  "dam.bw1_f",
                  "dam.bw2_f",
                  "dam.bw3_f",
                  "dam.bw1_m",
                  "dam.bw2_m",
                  "dam.bw3_m",
                  "dam.rfi1_m",
                  "dam.rfi2_m",
                  "dam.rfi1_f",
                  "dam.rfi2_f",
                  "dam.skin.length.male",
                  "dam.skin.length.female",
                  "dam.skin.qual",
                  "dam.live.qual",
                  "dam.h.length",
                  "pe1.bw.f",
                  "pe2.bw.f",
                  "pe3.bw.f",
                  "pe1.bw.m", 
                  "pe2.bw.m", 
                  "pe3.bw.m",
                  "birthyear.dam",
                  "sire.id.2nd",
                  "sire.fert.2nd",
                  "sire.bw1_f.2nd",
                  "sire.bw2_f.2nd",
                  "sire.bw3_f.2nd",
                  "sire.bw1_m.2nd",
                  "sire.bw2_m.2nd",
                  "sire.bw3_m.2nd",
                  "sire.rfi1_m.2nd",
                  "sire.rfi2_m.2nd",
                  "sire.rfi1_f.2nd",
                  "sire.rfi2_f.2nd",
                  "sire.skin.length.male.2nd",
                  "sire.skin.length.female.2nd",
                  "sire.skin.qual.2nd",
                  "sire.live.qual.2nd",
                  "sire.h.length.2nd",
                  "specific.env.skin",
                  "pe1.rfi.m", 
                  "pe2.rfi.m",
                  "pe1.rfi.f", 
                  "pe2.rfi.f"
                )                  
                , with=F] #specify which columns to incl.
  
  id <- seq(1:sum(x$obs_fert)) + max(z$id) # makes ID
  
  birthyear <- rep (year, sum(x$obs_fert)) # makes birthyear
  sex <- rbinom(sum(x$obs_fert),1,0.5)+1 # makes sex, TODO check if this is true
  true.sire <- numeric(nrow(kit.list))
  true.sire.check <- numeric(nrow(kit.list))
  mendelian <- rmvnorm(sum(x$obs_fert),sigma=G_sigma,method="svd")
  mendelian <- as.data.table(t(t(mendelian) * sqrt(0.5*variances))) 
  
  colnames(mendelian)     <-
    c(
      "mend.live.qual",
      "mend.h.length",
      "mend.skin.qual",
      "mend.skin.length.male",
      "mend.skin.length.female",
      "mend.litter.size",
      "mend.bw1_f",
      "mend.bw2_f",
      "mend.bw3_f",
      "mend.bw1_m",
      "mend.bw2_m",
      "mend.bw3_m",
      "mend.rfi1_m",
      "mend.rfi2_m",
      "mend.rfi1_f",
      "mend.rfi2_f"    )
  
  kit.list <-
    cbind (id,
           kit.list,
           birthyear,
           sex,
           mendelian,
           true.sire,
           true.sire.check) # binds id, sex and birthyear to data.table
  
  # now check the true sire before calculating inbreeding
  # here i use the true sire check in order to speed up, replace it later
  kit.list$true.sire.check <- ifelse(kit.list$sire.id.2nd == kit.list$sire.id.1st,TRUE,FALSE)
  kit.list$sire.fert.2nd <-
    ifelse(
      kit.list$true.sire.check == TRUE,
      kit.list$sire.fert.1st,
      kit.list$sire.fert.2nd
    )
  kit.list$sire.bw1_f.2nd <-
    ifelse(
      kit.list$true.sire.check == TRUE,
      kit.list$sire.bw1_f.1st,
      kit.list$sire.bw1_f.2nd
    )
  kit.list$sire.bw2_f.2nd <-
    ifelse(
      kit.list$true.sire.check == TRUE,
      kit.list$sire.bw2_f.1st,
      kit.list$sire.bw2_f.2nd
    )
  kit.list$sire.bw3_f.2nd <-
    ifelse(
      kit.list$true.sire.check == TRUE,
      kit.list$sire.bw3_f.1st,
      kit.list$sire.bw3_f.2nd
    )
  kit.list$sire.bw1_m.2nd <-
    ifelse(
      kit.list$true.sire.check == TRUE,
      kit.list$sire.bw1_m.1st,
      kit.list$sire.bw1_m.2nd
    )
  kit.list$sire.bw2_m.2nd <-
    ifelse(
      kit.list$true.sire.check == TRUE,
      kit.list$sire.bw2_m.1st,
      kit.list$sire.bw2_m.2nd
    )
  kit.list$sire.bw3_m.2nd <-
    ifelse(
      kit.list$true.sire.check == TRUE,
      kit.list$sire.bw3_m.1st,
      kit.list$sire.bw3_m.2nd
    )
  kit.list$sire.rfi1_m.2nd <-
    ifelse(
      kit.list$true.sire.check == TRUE,
      kit.list$sire.rfi1_m.1st,
      kit.list$sire.rfi1_m.2nd
    )
  kit.list$sire.rfi2_m.2nd <-
    ifelse(
      kit.list$true.sire.check == TRUE,
      kit.list$sire.rfi2_m.1st,
      kit.list$sire.rfi2_m.2nd
    )
  kit.list$sire.rfi1_f.2nd <-
    ifelse(
      kit.list$true.sire.check == TRUE,
      kit.list$sire.rfi1_f.1st,
      kit.list$sire.rfi1_f.2nd
    )
  kit.list$sire.rfi2_f.2nd <-
    ifelse(
      kit.list$true.sire.check == TRUE,
      kit.list$sire.rfi2_f.1st,
      kit.list$sire.rfi2_f.2nd
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
  kit.list$true.sire.check <-
    ifelse(kit.list$sire.id.1st != kit.list$sire.id.2nd,
           rbinom(nrow(kit.list), 1, true.sire.chance),
           1) # 85% chance that the kits are sired by 2nd mating
  kit.list$true.sire <-
    ifelse(kit.list$true.sire.check == 0,
           kit.list$sire.id.1st,
           kit.list$sire.id.2nd)
  kit.list[, `:=`(true.sire.fert = 
                    ifelse(kit.list$true.sire == kit.list$sire.id.2nd, 
                           kit.list$sire.fert.2nd, kit.list$sire.fert.1st), 
                  true.sire.bw1_m = 
                    ifelse(kit.list$true.sire == kit.list$sire.id.2nd, 
                           kit.list$sire.bw1_m.2nd, kit.list$sire.bw1_m.1st),
                  true.sire.bw2_m = 
                    ifelse(kit.list$true.sire == kit.list$sire.id.2nd, 
                           kit.list$sire.bw2_m.2nd, kit.list$sire.bw2_m.1st),
                  true.sire.bw3_m = 
                    ifelse(kit.list$true.sire == kit.list$sire.id.2nd, 
                           kit.list$sire.bw3_m.2nd, kit.list$sire.bw3_m.1st),
                  true.sire.bw1_f = 
                    ifelse(kit.list$true.sire == kit.list$sire.id.2nd, 
                           kit.list$sire.bw1_f.2nd, kit.list$sire.bw1_f.1st),
                  true.sire.bw2_f = 
                    ifelse(kit.list$true.sire == kit.list$sire.id.2nd, 
                           kit.list$sire.bw2_f.2nd, kit.list$sire.bw2_f.1st),
                  true.sire.bw3_f = 
                    ifelse(kit.list$true.sire == kit.list$sire.id.2nd, 
                           kit.list$sire.bw3_f.2nd, kit.list$sire.bw3_f.1st),
                  true.sire.rfi1_m = 
                    ifelse(kit.list$true.sire == kit.list$sire.id.2nd, 
                           kit.list$sire.rfi1_m.2nd, kit.list$sire.rfi1_m.1st),
                  true.sire.rfi2_m = 
                    ifelse(kit.list$true.sire == kit.list$sire.id.2nd, 
                           kit.list$sire.rfi2_m.2nd, kit.list$sire.rfi2_m.1st),
                  true.sire.rfi1_f = 
                    ifelse(kit.list$true.sire == kit.list$sire.id.2nd, 
                           kit.list$sire.rfi1_f.2nd, kit.list$sire.rfi1_f.1st),
                  true.sire.rfi2_f = 
                    ifelse(kit.list$true.sire == kit.list$sire.id.2nd, 
                           kit.list$sire.rfi2_f.2nd, kit.list$sire.rfi2_f.1st),              
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
                           kit.list$sire.h.length.2nd, kit.list$sire.h.length.1st)
  )]
  kit.list$sire.id.2nd <-
    ifelse(
      kit.list$sire.id.2nd == 0,
      kit.list$sire.id.1st,
      kit.list$sire.id.2nd
    )
  
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
  

  kit.list[, `:=`(litter.size = 0.5*(dam.fert + true.sire.fert) + 
                    mend.litter.size,  # Breeding value of offspring, littersize
                  bw1_f = 0.5*(true.sire.bw1_f + dam.bw1_f) + mend.bw1_f*(1-f0),
                  bw2_f = 0.5*(true.sire.bw2_f + dam.bw2_f) + mend.bw2_f*(1-f0),
                  bw3_f = 0.5*(true.sire.bw3_f + dam.bw3_f) + mend.bw3_f*(1-f0),
                  bw1_m = 0.5*(true.sire.bw1_m + dam.bw1_m) + mend.bw1_m*(1-f0),
                  bw2_m = 0.5*(true.sire.bw2_m + dam.bw2_m) + mend.bw2_m*(1-f0),
                  bw3_m = 0.5*(true.sire.bw3_m + dam.bw3_m) + mend.bw3_m*(1-f0),
                  rfi1_m = 0.5*(true.sire.rfi1_m + dam.rfi1_m) + mend.rfi1_m*(1-f0),
                  rfi2_m = 0.5*(true.sire.rfi2_m + dam.rfi2_m) + mend.rfi2_m*(1-f0),
                  rfi1_f = 0.5*(true.sire.rfi1_f + dam.rfi1_f) + mend.rfi1_f*(1-f0),
                  rfi2_f = 0.5*(true.sire.rfi2_f + dam.rfi2_f) + mend.rfi2_f*(1-f0),
                  perm.env.ls = rnorm(sum(x$obs_fert))*sqrt(var.perm.env.ls), # perm env for litter size
                  skin.length.male = 0.5*(dam.skin.length.male + true.sire.skin.length.male) + 
                    mend.skin.length.male*(1-f0),
                  skin.length.female = 0.5*(dam.skin.length.female + true.sire.skin.length.female) + 
                    mend.skin.length.female*(1-f0),
                  live.qual = 0.5*(dam.live.qual + true.sire.live.qual)+ (mend.live.qual)*(1-f0),
                  skin.qual = 0.5*(dam.skin.qual + true.sire.skin.qual)+ (mend.skin.qual)*(1-f0),
                  h.length = 0.5*(dam.h.length + true.sire.h.length)+ (mend.h.length)*(1-f0)
  )]# Breeding value of offspring, body size
  t <- c(63,84,105,126,147,168,189,210)
  t <- StandardTime(t)
  t <- t[3:8]
  
  q <- as.matrix(as.data.frame(polynomial.values(polynomials = leg2, x =t[6])))
  
  kit.list[, `:=`(add.gen.bw.f =  q[1]*bw1_f+q[2]*bw2_f+q[3]*bw3_f,
                  add.gen.bw.m =  q[1]*bw1_m+q[2]*bw2_m+q[3]*bw3_m)]
  # kit.list$dam.age <- ifelse( z - kit.list$birthyear.dam > 1, 1,0 )
  
  
  setnames(kit.list, "obs_fert", "own_littersize") # changes obs fert into the littersize of the kit

  kit.list[, `:=`( 
    phenotype.live.qual = live.qual  + rnorm(nrow(kit.list)) *
      sqrt(var.live.qual.res),
    phenotype.skin.length = ifelse(
      sex == 1,
      mean.skin.length.male + skin.length.male + rnorm(1)*sqrt(var.skin.length.res.male)+
        specific.env.skin*sqrt(var.skin.length.c.male),
      mean.skin.length.female + skin.length.female + rnorm(1)*sqrt(var.skin.length.res.female)+
        specific.env.skin*sqrt(var.skin.length.c.female) ),
    phenotype.skin.qual = 
      skin.qual + rnorm(nrow(kit.list))*sqrt(var.skin.qual.res),
    phenotype.h.length = h.length+ rnorm(nrow(kit.list))*
      sqrt(var.h.length.res)
  )
  ]  
  set( kit.list, j=which(colnames(kit.list) %in% c(
    "sire.fert.1st",
    "sire.bw1_f.1st",
    "sire.bw2_f.1st",
    "sire.bw3_f.1st",
    "sire.bw1_m.1st",
    "sire.bw2_m.1st",
    "sire.bw3_m.1st",
    "sire.rfi1_m.1st",
    "sire.rfi2_m.1st",
    "sire.rfi1_f.1st",
    "sire.rfi2_f.1st",
    "sire.skin.length.male.1st",
    "sire.skin.length.female.1st",
    "sire.skin.qual.1st",
    "sire.live.qual.1st",
    "sire.bw.june.1st",
    "sire.bw.june.maternal.1st",
    "sire.h.length.1st",
    "sire.id.1st",
    "dam.fert",
    "dam.bw1_f",
    "dam.bw2_f",
    "dam.bw3_f",
    "dam.bw1_m",
    "dam.bw2_m",
    "dam.bw3_m",
    "dam.rfi1_m",
    "dam.rfi2_m",
    "dam.rfi1_f",
    "dam.rfi2_f",
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
    "sire.bw1_f.2nd",
    "sire.bw2_f.2nd",
    "sire.bw3_f.2nd",
    "sire.bw1_m.2nd",
    "sire.bw2_m.2nd",
    "sire.bw3_m.2nd",
    "sire.rfi1_m.2nd",
    "sire.rfi2_m.2nd",
    "sire.rfi1_f.2nd",
    "sire.rfi2_f.2nd",
    "sire.skin.length.male.2nd",
    "sire.skin.length.female.2nd",
    "sire.skin.qual.2nd",
    "sire.live.qual.2nd",
    "sire.bw.june.2nd",
    "sire.bw.june.maternal.2nd",
    "sire.h.length.2nd",
    "mend.bw1_f",
    "mend.bw2_f",
    "mend.bw3_f",
    "mend.bw1_m",
    "mend.bw2_m",
    "mend.bw3_m",
    "mend.rfi1_m",
    "mend.rfi2_m",
    "mend.rfi1_f",
    "mend.rfi2_f",
    "mend.litter.size", 
    "mend.live.qual",
    "mend.skin.qual", 
    "mend.skin.length.male",
    "mend.skin.length.female",
    "true.sire.check",
    "true.sire.fert",
    "true.sire.bw1_f",
    "true.sire.bw2_f",
    "true.sire.bw3_f",
    "true.sire.bw1_m",
    "true.sire.bw2_m",
    "true.sire.bw3_m",
    "true.sire.rfi1_m",
    "true.sire.rfi2_m",
    "true.sire.rfi1_f",
    "true.sire.rfi2_f",
    "true.sire.skin.length.male",
    "true.sire.skin.length.female",
    "true.sire.skin.qual",
    "true.sire.live.qual",
    "true.sire.h.length"
    # "dam.age"
  )) , value=NULL ) # removes bv of parents
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
YearlingEffectOnFertility <- function (x,y,yearling.effect){ # x = mating.list, y = year
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

CalculateBLUP <- function () {
  system2("run_dmu4.bat", " MBLUP")
  # read the solutions and only keep the predictions of BV (they're not that right)
  
  solutions <- as.matrix(read.table(file="MBLUP.SOL"))
  solutions <- as.data.table(solutions)
  solutions <- subset(solutions, V1 == 4 & V2==1 ) # BW
  set (solutions, j=c("V1","V2","V3", "V4", "V6", "V7","V9"), value= NULL)
  setnames(solutions, c("V5","V8"),c("id", "blup.bwnov"))
  temp <- as.matrix(read.table(file="MBLUP.SOL"))
  temp <- as.data.table(temp)
  temp <- subset(temp, V1 == 4 & V2==2 ) # qual
  set (temp, j=c("V1","V2","V3", "V4", "V6", "V7","V9"), value= NULL)
  setnames(temp, c("V5","V8"),c("id", "blup.qual"))
  solutions <- merge(solutions, temp, by="id")
  temp <- as.matrix(read.table(file="MBLUP.SOL"))
  temp <- as.data.table(temp)
  temp <- subset(temp, V1 == 4 & V2==3 ) # litter size
  set (temp, j=c("V1","V2","V3", "V4", "V6", "V7","V9"), value= NULL)
  setnames(temp, c("V5","V8"),c("id", "blup.fert"))
  solutions <- merge(solutions, temp, by="id")
  
  return(solutions)
}
############### Modify Dir file  ################
 # this function changes the .DIR file for fertility and BW
# So it takes in the replicate from the current replicate of the simulation
 ModifyDIRFile <- function (p, mblup, trace.ped) {
   # fertility modification
   if (mblup == 0) {
   if (trace.ped == 0 ){
     
     # bw modification
     dirfile <- readLines("MBLUP.DIR")
     dirfile[8] <- c(paste("$DATA  ASCII (8,3,-9999.0) Phenotypes",p, sep="")) 
     # change the input file for BLUP so it uses the next outputfile
     writeLines(dirfile, "MBLUP.DIR")
     dirfile[30] <- c(paste("$VAR_STR 1 PED 2 ASCII Big_pedigree_",p, sep="")) 
     # change the input file for BLUP so it uses the next pedigree
     writeLines(dirfile,"MBLUP.DIR")
     # change the input file for BLUP so it uses the next pedigree
   } else if (trace.ped == 1) {
     dirfile <- readLines("MBLUP.DIR")
     dirfile[8] <- c(paste("$DATA  ASCII (8,3,-9999.0) Phenotypes",p, sep="")) 
     # change the input file for BLUP so it uses the next outputfile
     writeLines(dirfile, "MBLUP.DIR")
     
     dirfile <- readLines("trace.DIR")
     dirfile[1] <- c(paste("Big_pedigree_",p, sep=""))
     writeLines(dirfile,"trace.DIR")
     # change the pedigree name in the dir file
     dirfile <- readLines("MBLUP.DIR")
     dirfile[30] <- c(paste("$VAR_STR 1 PED 2 ASCII trace.PRUNE")) # change the input file for BLUP so it uses the next pedigree
     writeLines(dirfile,"MBLUP.DIR")
     # change the pedigree name in case of trace being asked for
     
   }
   } else if (mblup ==1 ) {
     if (trace.ped == 0 ){
       
       # bw modification
       dirfile <- readLines("MBLUP_full.DIR")
       dirfile[8] <- c(paste("$DATA  ASCII (8,4,-9999.0) MBLUP_Y",p, sep="")) 
       # change the input file for BLUP so it uses the next outputfile
       writeLines(dirfile, "MBLUP_full.DIR")
       dirfile[36] <- c(paste("$VAR_STR 1 PED 2 ASCII Big_pedigree_",p, sep="")) 
       # change the input file for BLUP so it uses the next pedigree
       writeLines(dirfile,"MBLUP_full.DIR")
       # change the input file for BLUP so it uses the next pedigree
     } else if (trace.ped == 1) {
       dirfile <- readLines("MBLUP_full.DIR")
       dirfile[8] <- c(paste("$DATA  ASCII (8,4,-9999.0) MBLUP_Y",p, sep="")) 
       # change the input file for BLUP so it uses the next outputfile
       writeLines(dirfile, "MBLUP_full.DIR")
       
       dirfile <- readLines("trace.DIR")
       dirfile[1] <- c(paste("Big_pedigree_",p, sep=""))
       writeLines(dirfile,"trace.DIR")
       # change the pedigree name in the dir file
       dirfile <- readLines("MBLUP.DIR")
       dirfile[36] <- c(paste("$VAR_STR 1 PED 2 ASCII trace.PRUNE")) # change the input file for BLUP so it uses the next pedigree
       writeLines(dirfile,"MBLUP.DIR")
       # change the pedigree name in case of trace being asked for
       
     }
   
 } 
} # function
 
############### Create observation file for phenotypes #################
# currently one phenotype file is made
WriteObservations <- function (mating.list, next.gen, next.gen.males,kit.list,year,p,sorting.prop) {
  matinglist <- mating.list[,c("dam.id","obs_fert"),with=FALSE]
  setnames(matinglist, "dam.id", "id")
  if("obs_fert" %in% colnames(next.gen)) {
    set( next.gen, j=which(colnames(next.gen) %in% 
                             "obs_fert")  , value=NULL )
  }
  next.gen <- merge(next.gen, matinglist, by="id",all.x=TRUE)
  
  # if (year < 2) {
  # writefile[,`:=`(dam.id= as.integer(0))] 
  # }
  if (year == 1){
    next.gen.males$own_littersize= as.integer(0)
    next.gen$own_littersize = as.integer(0)
  } # make up littersizes for the base gen, probably ok to set it to a constant as it will then be regressed to zero
  
  writefile <- next.gen[, c("id","dam.id", "phenotype.bw.oct", "live.score", "sex", "birthyear", "obs_fert","own_littersize"), with=FALSE]
  writefile[,`:=`(dam.age = ifelse( year-writefile$birthyear == 1, as.integer(1), as.integer(2)))]
  
  
  next.gen.males[,`:=`(dam.age= as.integer(0),obs_fert = as.integer(-9999))] 
  temp <- next.gen.males[, c("id","dam.id", "phenotype.bw.oct", "live.score", "dam.age","sex", "birthyear", "obs_fert","own_littersize"),with=F]
  writefile <- rbind(writefile, temp) #rbind the males first
  if (year >1) { # this is to "remove" the body weight and quality measurement for females so they do not appear twice
    writefile$phenotype.bw.oct <- as.integer(-9999)
    writefile$live.score <- as.integer(-9999)
    
  }
  kit.list[,`:=`(dam.age= as.integer(0),obs_fert = as.integer(-9999))] 
  if (sorting.prop != 1) {
    # hingað
    numb.animals.to.sort <- ceiling(nrow(kit.list)*(1-sorting.prop))
   animals.to.mask <-  kit.list[,c("id"),with=FALSE][sample(.N,numb.animals.to.sort)]
  animals.to.mask <- merge(animals.to.mask,kit.list,by="id", all.x=TRUE)
  animals.to.mask[,`:=`(phenotype.bw.oct= as.integer(-9999),live.score = as.integer(-9999))]
  sd <- setdiff(kit.list$id, animals.to.mask$id)
  sd <- is.element(kit.list$id, sd)
  kit.list <- kit.list[sd, ]
  kit.list <- rbind(kit.list,animals.to.mask)
  
   }
  temp <- kit.list[, c("id", "dam.id","phenotype.bw.oct", "live.score", "dam.age","sex", "birthyear","own_littersize" ,"obs_fert"), with=FALSE]
  writefile <- rbind (writefile, temp)
  writefile[is.na(writefile)]	=	as.integer(-9999)
  writefile[,`:=`(inter= as.integer(1), year= as.integer(year))]
  writefile[, c("dam.id","sex","birthyear","id","own_littersize")
            :=lapply(.SD, function(x) as.integer(x)), .SDcols=c("dam.id","sex","birthyear","id","own_littersize")]
  
  if (year == 1) {
    phenotypes <- file(description = paste("Phenotypes",p, sep=""), open="w")
  } else if (year > 1) {
    phenotypes <- file(description = paste("Phenotypes",p, sep=""), open="a")
    
  }
  write.table(format(writefile[,.(id,dam.id,year,dam.age,birthyear,sex,own_littersize,inter,phenotype.bw.oct,live.score,obs_fert)], nsmall=1, digits=2), 
              file= phenotypes, append= TRUE,col.names = FALSE, row.names = FALSE, quote = FALSE)
  close(con=phenotypes)
  
}

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
############### Trace pedigree ###############
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
############### Mask kits from big litters ######
# the role of this function is to "hide" the kits from the big litters
# this is because farmers move kits from big litters to small litters
# this is essentially a simplification since the moved kits would have a portion
# of the "caretaker" permanent environment. That however is to conveluted to go into
# right now. This is a simple way to avoid the problem.
MaskKits <- function (kitlist) {
  setkey(kitlist, id)
  setorder(kitlist, id, sex)
  kitlist[, `:=`(IDX = 1:.N) , by = dam.id]
  kitlist[,`:=`(mask = ifelse( IDX> 9, 0,1))]
  kitlist <- subset(kitlist, mask == 1  )
  kitlist[,c("IDX", "mask"):=NULL]
  
  return(kitlist)
}
############### Random Culling ##############
# random culling is implemented since survival is not very heritable and a lot of effort
# TODO make more variability in 
RandCull <- function (kitlist,cull.ratio) {

    kitlist[, SO:=rbinom(nrow(kitlist), 1, cull.ratio)]
  kitlist <- subset(kitlist, SO == 1  )
  return(kitlist)
}

 
############### Residual feed intake ###########
 RFI <- function (kit.list, leg2,leg1, t) {
   t <- c(63,84,105,126,147,168,189,210)
   t <- StandardTime(t)
   t <- t[3:8]
   temp <-
     kit.list[, c("id",
                  "sex",
                  "bw1_f",
                  "bw2_f",
                  "bw3_f",
                  "bw1_m",
                  "bw2_m",
                  "bw3_m",
                  "pe1.bw.f",
                  "pe2.bw.f",
                  "pe3.bw.f",
                  "pe1.bw.m",
                  "pe2.bw.m",
                  "pe3.bw.m"), with = FALSE]
   # scale variances
   q <- as.matrix(as.data.frame(polynomial.values(polynomials = leg2, x =t)))
   FR.females <- c(1, 0.3, -0.12) 
   FR.males <- c(1, .94, 1.05)
   
   const.m <- q %*% FR.males
   const.f <- 0.9 + q %*% FR.females # this 0.9 is to get the body weight to a "realistic" level
   
   # this step makes the 6 needed BW for the RFI 
   # note the commented version has permanent environment 
   # temp[,`:=`( 
   #   phenotype.bw.oct = ifelse(sex==1, const.m[6] + bw1_m*q[[6,1]]+bw2_m*q[[6,2]]+bw3_m*q[[6,3]]+
   #                               pe1.bw.m*q[[6,1]]+pe2.bw.m*q[[6,2]]+pe3.bw.m*q[[6,3]]+rnorm(nrow(temp))*sqrt(bw.res.male[8]),
   #                             const.f[6] + bw1_f*q[[6,1]]+bw2_f*q[[6,2]]+bw3_f*q[[6,3]]+
   #                               pe1.bw.f*q[[6,1]]+pe2.bw.f*q[[6,2]]+pe3.bw.f*q[[6,3]]+rnorm(nrow(temp))*sqrt(bw.res.female[8])),
   #   perm.env = ifelse( sex == 1, pe1.bw.m*q[[6,1]]+pe2.bw.m*q[[6,2]]+pe3.bw.m*q[[6,3]], pe1.bw.f*q[[6,1]]+pe2.bw.f*q[[6,2]]+pe3.bw.f*q[[6,3]])
   #   
   # )]  
   temp[,`:=`( 
     phenotype.bw.oct = ifelse(sex==1, const.m[6] + bw1_m*q[[6,1]]+bw2_m*q[[6,2]]+bw3_m*q[[6,3]]+
                                 rnorm(nrow(temp))*sqrt(bw.res.male[8]),
                               const.f[6] + bw1_f*q[[6,1]]+bw2_f*q[[6,2]]+bw3_f*q[[6,3]]+
                                 rnorm(nrow(temp))*sqrt(bw.res.female[8])),
     perm.env = ifelse( sex == 1, pe1.bw.m*q[[6,1]]+pe2.bw.m*q[[6,2]]+pe3.bw.m*q[[6,3]], pe1.bw.f*q[[6,1]]+pe2.bw.f*q[[6,2]]+pe3.bw.f*q[[6,3]])
     
   )]  
   
   temp <- temp[,c("id","sex","phenotype.bw.oct"),with=F]
   temp <- merge(temp, 
                 kit.list[,
                          c("id",
                            "rfi1_m",
                            "rfi2_m",
                            "rfi1_f",
                            "rfi2_f",
                            "pe1.rfi.m",
                            "pe2.rfi.m",
                            "pe1.rfi.f",            
                            "pe2.rfi.f"  
                          ),with=F], by="id")
   leg1 <- legendre.polynomials(1, normalized = T)
   q <- as.matrix(as.data.frame(polynomial.values(polynomials = leg1, x =t)))
   const.rfi <- q %*% FR.RFI
   temp[,`:=`(
     FI = ifelse( sex == 1, 
                  phenotype.bw.oct*b.bw.male + const.rfi[6]+ 
                    q[[6,1]]*rfi1_m + q[[6,2]]*rfi2_m+ q[[6,1]]*pe1.rfi.m+q[[6,2]]*pe2.rfi.m+  
                    rnorm(nrow(temp))*sqrt(res.rfi[6]),
                  phenotype.bw.oct*b.bw.female+ const.rfi[6]+
                    q[[6,1]]*rfi1_f+ q[[6,2]]*rfi2_f+ q[[6,1]]*pe1.rfi.f+q[[6,2]]*pe2.rfi.f+
                    rnorm(nrow(temp))*sqrt(res.rfi[6]))
   )]
   kit.list <- merge(kit.list,temp[,c("id",
                                      "phenotype.bw.oct",
                                      "FI"), with=FALSE], by="id")
   # kit.list$phenotype.bw.oct <- kit.list$phenotype.bw.oct*1000
   set( kit.list, j=which(colnames(kit.list) %in% 
                     c("pe1.rfi.m",
                       "pe2.rfi.m",
                       "pe1.rfi.f",            
                       "pe2.rfi.f",  
                       "pe1.bw.f",
                       "pe2.bw.f",
                       "pe3.bw.f",
                       "pe1.bw.m",
                       "pe2.bw.m",
                       "pe3.bw.m"
                     ))  , value=NULL )
   
   return(kit.list) 
 }
 ############### Write observation file for body weight ####################
 # defunct
 WriteObservationFileBodyWeight <- function (x,year,p) {
   inter <- as.integer(rep(1, times=nrow(x))) # intercept for the quality regression
   x<- cbind(x, inter)
   x[, c("id","dam.id","sex")
     :=lapply(.SD, function(x) as.integer(x)), .SDcols=c("id","dam.id","sex")]
   # x[,`:=`(prodyear=as.integer(ifelse( sex == 1, year*(sex+6),year*(sex+8))))]
   # NOTE if simulation runs to 81 years, then year*sex will be the same for year 63 for dams and 81 for males
   if (year  == 1 ) {test <- file(description = paste("Test",p, sep=""), open="w")
   if (weighing.method == oct){
     x[, c("id","dam.id","own_littersize","sex","inter")
       :=lapply(.SD, function(x) as.integer(x)), .SDcols=c("id","dam.id","own_littersize","sex","inter")]
     
     write.table(format(x[,.(id,dam.id,sex,inter,phenotype.bw.oct)], nsmall=1, digits=2), 
                 file= test, append= TRUE,col.names = FALSE, row.names = FALSE, quote = FALSE)
   } else if (weighing.method == sept ) {
     x[, c("id","dam.id","own_littersize","sex","inter")
       :=lapply(.SD, function(x) as.integer(x)), .SDcols=c("id","dam.id","own_littersize","sex","inter")]
     
     write.table(format(x[,.(id,dam.id,sex,inter,phenotype.bw.sept)], nsmall=1, digits=2), 
                 file= test, append= TRUE,col.names = FALSE, row.names = FALSE, quote = FALSE)
   }
   close(con=test)
   } else if (year > 1) {
     x[, c("id","dam.id","own_littersize","sex","inter")
       :=lapply(.SD, function(x) as.integer(x)), .SDcols=c("id","dam.id","own_littersize","sex","inter")]
     
     test <- file(description = paste("Test",p, sep=""), open="a")
     if (weighing.method == oct){
       write.table(format(x[,.(id,dam.id,sex,inter,phenotype.bw.oct)], nsmall=1, digits=2), 
                   file= test, append= TRUE,col.names = FALSE, row.names = FALSE, quote = FALSE)
     } else if (weighing.method == sept ) {
       write.table(format(x[,.(id,dam.id,sex,inter,phenotype.bw.sept)], nsmall=1, digits=2), 
                   file= test, append= TRUE,col.names = FALSE, row.names = FALSE, quote = FALSE)
     }
     close(con=test)
     
   }
 }
 
 ################ MBLUP write observation file #####################
 WriteMBLUPObservations <- function (mating.list, next.gen, next.gen.males,kit.list,year,p) {
   matinglist <- mating.list[,c("dam.id","obs_fert"),with=FALSE]
   setnames(matinglist, "dam.id", "id")
   if("obs_fert" %in% colnames(next.gen)) {
     set( next.gen, j=which(colnames(next.gen) %in% 
                              "obs_fert")  , value=NULL )
   }
   next.gen <- merge(next.gen, matinglist, by="id",all.x=TRUE)
   
   # if (year < 2) {
   # writefile[,`:=`(dam.id= as.integer(0))] 
   # }
   if (year == 1){
     next.gen.males$own_littersize= as.integer(0)
     next.gen$own_littersize = as.integer(0)
   } # make up littersizes for the base gen, probably ok to set it to a constant as it will then be regressed to zero
   
   writefile <- next.gen[, c("id","dam.id", "phenotype.bw.oct", "live.score", "sex", "birthyear", "obs_fert","own_littersize"), with=FALSE]
   writefile[,`:=`(dam.age = ifelse( year-writefile$birthyear == 1, as.integer(1), as.integer(2)))]
   
   
   next.gen.males[,`:=`(dam.age= as.integer(0),obs_fert = as.integer(-9999))] 
   temp <- next.gen.males[, c("id","dam.id", "phenotype.bw.oct", "live.score", "dam.age","sex", "birthyear", "obs_fert","own_littersize"),with=F]
   writefile <- rbind(writefile, temp) #rbind the males first
   if (year >1) { # this is to "remove" the body weight and quality measurement for females so they do not appear twice
     writefile$phenotype.bw.oct <- as.integer(-9999)
     writefile$live.score <- as.integer(-9999)
     
   }
   kit.list[,`:=`(dam.age= as.integer(0), obs_fert = as.integer(-9999))] 
   temp <- kit.list[, c("id", "dam.id","phenotype.bw.oct", "live.score", "dam.age","sex", "birthyear","own_littersize" ,"obs_fert"), with=FALSE]
   writefile <- rbind (writefile, temp)
   writefile[is.na(writefile)]	=	as.integer(-9999)
   writefile[,`:=`(inter= as.integer(1), year= as.integer(year))]
   writefile[, c("dam.id","sex","birthyear","id","own_littersize")
             :=lapply(.SD, function(x) as.integer(x)), .SDcols=c("dam.id","sex","birthyear","id","own_littersize")]
   
   if (year == 1) {
     MBLUP_Y <- file(description = paste("MBLUP_Y",p, sep=""), open="w")
   } else if (year > 1) {
     MBLUP_Y <- file(description = paste("MBLUP_Y",p, sep=""), open="a")
     
   }
   writefile[,`:=`(bw.female = ifelse(sex == 2, phenotype.bw.oct, -9999.0),
                  phenotype.bw.oct = ifelse(sex ==1, phenotype.bw.oct, -9999.0))]
   write.table(format(writefile[,.(id,dam.id,year,dam.age,birthyear,sex,own_littersize,inter,phenotype.bw.oct,live.score,obs_fert,bw.female)], nsmall=1, digits=2), 
               file= MBLUP_Y, append= TRUE,col.names = FALSE, row.names = FALSE, quote = FALSE)
   close(con=MBLUP_Y)
   
 }
 
 ######################## Skin price function ##############################
 SkinPrices <- function (kitlist, next.gen, next.gen.males,y) {
   
   sd <-
     setdiff(kitlist$id, next.gen$id) # remove the next.gen females from kit.list
   sd <- is.element(kitlist$id, sd)
   kitlist <- kitlist[sd, ]
   sd <- setdiff(kitlist$id, next.gen.males$id)
   sd <- is.element(kitlist$id, sd)
   kitlist <- kitlist[sd, ]
   kitlist[, `:=`(
     
     P1 = ifelse (
       phenotype.skin.length >= 101,
       1,0),#50
     P2= ifelse(
       phenotype.skin.length < 101 & phenotype.skin.length >= 95,
       1,0),#40
     P3=ifelse(
       phenotype.skin.length < 95 & phenotype.skin.length >= 89,
       1,0),#30
     P4 = ifelse (
       phenotype.skin.length < 89 & phenotype.skin.length >= 83,
       1,0),#00
     P5=ifelse(
       phenotype.skin.length < 83 &
         phenotype.skin.length >= 77,
       1,0),#0
     P6=ifelse(
       phenotype.skin.length < 77 & phenotype.skin.length >= 71,
       1,0),#1
     P7 = ifelse (
       phenotype.skin.length < 71 & phenotype.skin.length >= 65,
       1,0),#2
     P8 = ifelse(
       phenotype.skin.length < 65 & phenotype.skin.length >= 59,
       1,0),#3
     P9=ifelse(
       phenotype.skin.length < 59 &
         phenotype.skin.length >= 53,
       1,0),#4
     P10 = ifelse(phenotype.skin.length < 53, 1, 0),#5
     P11 = ifelse(
       phenotype.skin.qual >= truncs[3],
       1,0),
     # purple
     P12 = ifelse(
       truncs[2] < phenotype.skin.qual & phenotype.skin.qual <= truncs[3],
       1,0),
     # platinum
     P13= ifelse(
       phenotype.skin.qual > truncs[1] & phenotype.skin.qual <= truncs[2],
       1,0),
     # burgundy
     P14 = ifelse(phenotype.skin.qual <=
                    truncs[1], 1, 0), # ivory
     P15 = ifelse(phenotype.h.length > htruncs[4] ,
                  1,0), #vel3
     P16 =  ifelse(phenotype.h.length > htruncs[3] & phenotype.h.length < htruncs[4],
                   1,0),
     # velv2
     P17 = ifelse(
       phenotype.h.length > htruncs[2] & phenotype.h.length < htruncs[3],
       1,0),
     # vel1
     P18=ifelse(
       phenotype.h.length > htruncs[1] & phenotype.h.length < htruncs[2],
       1,0),
     # klassik
     P19 = ifelse(phenotype.h.length < htruncs[1], 1, 0) # long nap
     
   )]
   # the loop function here below is to make sure that the skins are allocated
   # into legal classes, i.e. no purple skins in velvet 3 and no purple skins
   # in klassik
for (k in 1:nrow(kitlist)) {
  if (kitlist$P15[k] == 1) {
    set( kitlist,i=k, j=which(colnames(kitlist) %in% c("P11","P12","P14")), value=0)
    kitlist$P13[k] <- 1                       
  }
  if(kitlist$P18[k] == 1 & kitlist$P11[k] == 1) {
    
    set( kitlist,i=k, j=which(colnames(kitlist) %in% c("P11")), value=0)
    set( kitlist,i=k, j=which(colnames(kitlist) %in% c("P12")), value=1)
  }
  if (kitlist$P19[k] == 1 ) {
    set( kitlist,i=k, j=which(colnames(kitlist) %in% c("P11","P12","P13")), value=0)
    set( kitlist,i=k, j=which(colnames(kitlist) %in% c("P14")), value=1)
  }
}
kits.males <- subset(kitlist, sex == 1, select=
                          c("id",
                            "P1", #50
                            "P2",#40
                            "P3",#30
                            "P4",#00
                            "P5",#0
                            "P6",#1
                            "P7",#2
                            "P8",#3
                            "P9",#4
                            "P10",#5
                            "P11",#purple
                            "P12",#platinum
                            "P13",#burgundy
                            "P14",#ivory
                            "P15",#velv3
                            "P16",#velv2
                            "P17",#vel1
                            "P18",#kl
                            "P19")#long nap
   )
skin.metric.bin <- colSums(kits.males)
skin.metrics.males <- file(description = "skin_metrics_males", open ="a")
skin.metrics.females <- file(description = "skin_metrics_females", open ="a")

cat(
  y,
  skin.metric.bin[2],
  skin.metric.bin[3],
  skin.metric.bin[4],
  skin.metric.bin[5],
  skin.metric.bin[6],
  skin.metric.bin[7],
  skin.metric.bin[8],
  skin.metric.bin[9],
  skin.metric.bin[10],
  skin.metric.bin[11],
  skin.metric.bin[12],
  skin.metric.bin[13],
  skin.metric.bin[14],
  skin.metric.bin[15],
  skin.metric.bin[16],
  skin.metric.bin[17],
  skin.metric.bin[18],
  skin.metric.bin[19],
  skin.metric.bin[20],
  "\t", 
  sep="\t",
file = skin.metrics.males
  )
   # prices for 4 and 5 are made up (size)
   skin.prices.males <- intercept.bm+(as.matrix(kits.males))[,2:20] %*% (prices.bmales)
   
   cat(
     mean(skin.prices.males),
     nrow(skin.prices.males),
     "\n",
     sep="\t",
     file = skin.metrics.males
   )
   skin.prices.males<- cbind(kits.males$id, skin.prices.males)
    
    
   kits.females <- subset(kitlist, sex == 2, select=
                            c("id",
                              "P1", #50
                              "P2",#40
                              "P3",#30
                              "P4",#00
                              "P5",#0
                              "P6",#1
                              "P7",#2
                              "P8",#3
                              "P9",#4
                              "P10",#5
                              "P11",#purple
                              "P12",#platinum
                              "P13",#burgundy
                              "P14",#ivory
                              "P15",#velv3
                              "P16",#velv2
                              "P17",#vel1
                              "P18",#kl
                              "P19")#long nap
   )
   skin.metric.bin <- colSums(kits.females)
   
   cat(
     y,
     skin.metric.bin[2],
     skin.metric.bin[3],
     skin.metric.bin[4],
     skin.metric.bin[5],
     skin.metric.bin[6],
     skin.metric.bin[7],
     skin.metric.bin[8],
     skin.metric.bin[9],
     skin.metric.bin[10],
     skin.metric.bin[11],
     skin.metric.bin[12],
     skin.metric.bin[13],
     skin.metric.bin[14],
     skin.metric.bin[15],
     skin.metric.bin[16],
     skin.metric.bin[17],
     skin.metric.bin[18],
     skin.metric.bin[19],
     skin.metric.bin[20],
     "\t", 
     sep="\t",
     file = skin.metrics.females
   )
   # prices for 4 and 5 are made up (size)
   skin.prices.females <- intercept.bfemales+(as.matrix(kits.females))[,2:20] %*% (prices.bfemales)
   cat(
     mean(skin.prices.females),
     nrow(skin.prices.females),
     "\n",
     sep="\t",
     file = skin.metrics.females
   )
   skin.prices.females<- cbind(kits.females$id, skin.prices.females)
   skin.prices <- rbind(skin.prices.males, skin.prices.females)
   colnames(skin.prices) <- c("id","skin.price")
   kitlist <- merge(kitlist, skin.prices, by="id")
   set( kitlist, j=which(colnames(kitlist) %in% 
                           c("id",
                             "P1", #50
                             "P2",#40
                             "P3",#30
                             "P4",#00
                             "P5",#0
                             "P6",#1
                             "P7",#2
                             "P8",#3
                             "P9",#4
                             "P10",#5
                             "P11",#purple
                             "P12",#platinum
                             "P13",#burgundy
                             "P14",#ivory
                             "P15",#velv3
                             "P16",#velv2
                             "P17",#vel1
                             "P18",#kl
                             "P19"))  , value=NULL )
   return(kitlist)
 } 
############# calculate MBLUP ############
 CalculateMBLUP <- function () {
   system2("run_dmu4.bat", " MBLUP_full")
   # read the solutions and only keep the predictions of BV (they're not that right)
   
   solutions <- as.matrix(read.table(file="MBLUP_full.SOL"))
   solutions <- as.data.table(solutions)
   solutions <- subset(solutions, V1 == 4 & V2==1 ) # BW
   set (solutions, j=c("V1","V2","V3", "V4", "V6", "V7","V9"), value= NULL)
   setnames(solutions, c("V5","V8"),c("id", "blup.bwnov"))
   temp <- as.matrix(read.table(file="MBLUP_full.SOL"))
   temp <- as.data.table(temp)
   temp <- subset(temp, V1 == 4 & V2==2 ) # qual
   set (temp, j=c("V1","V2","V3", "V4", "V6", "V7","V9"), value= NULL)
   setnames(temp, c("V5","V8"),c("id", "blup.qual"))
   solutions <- merge(solutions, temp, by="id")
   temp <- as.matrix(read.table(file="MBLUP_full.SOL"))
   temp <- as.data.table(temp)
   temp <- subset(temp, V1 == 4 & V2==3 ) # litter size
   set (temp, j=c("V1","V2","V3", "V4", "V6", "V7","V9"), value= NULL)
   setnames(temp, c("V5","V8"),c("id", "blup.fert"))
   solutions <- merge(solutions, temp, by="id")
   
   return(solutions)
 }
 ################# Write the log file #####################
 WriteLogFile <-
   function ( n.females,
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
                 intensity.remating
   ) 
     {
     logfile <- file(description = "log.log", open = "w")
   cat("Logfile from MinkSim",file=logfile)
   cat("\n", file = logfile)
   cat("Simulation started",file=logfile,sep = "\t")
   cat(format(Sys.time(), " %b %d %X"),file=logfile,sep = "\t")
   cat("\n", file = logfile)
   cat("Simulation ended",file=logfile,sep = "\t")
   cat("\n", file = logfile)
   
   if (selection.method == phenotypic) {
     cat("Selection method is phenotypic", sep="\t",file=logfile)
     cat("\n", file = logfile)  
   } else if (selection.method == blup & mblup == 0) {
     cat("selection method is single-trait BLUP",sep="\t",file=logfile)
     cat("\n", file = logfile)
   } else if (selection.method== blup & mblup == 1) {
     cat("Selection method is multi-trait BLUP", sep="\t",file=logfile)
     cat("\n", file = logfile)
   }
   cat("Selection criterion", file=logfile, sep="\n")
   if (selection.method == phenotypic) {
     cat(
       " Proportion of animals deselected due to litter size ",
       quantile.setting.ls,
       file = logfile,
       sep = "\t"
     )
     cat("\n", file = logfile)
     cat(
       " Proportion of animals deselected due to body weight ",
       quantile.setting.bw,
       file = logfile,
       sep = "\t"
     )
     cat("\n", file = logfile)
     cat(
       " Proportion remaining for selection on quality",
       ((1 - quantile.setting.ls) * quantile.setting.bw),
       file = logfile,
       sep = "\t"
     )
     cat("\n", file = logfile)
   } else if (selection.method == blup) {
     cat("Index weights on kits", file=logfile)
     cat("\n", file = logfile)
     cat("Index weight on litter size ", weight.fert.kits, file =logfile, sep="\t")
     cat("\n", file = logfile)
     cat("Index weight on body weight ", weight.bw.kits, file =logfile, sep="\t")
     cat("\n", file = logfile)
     cat("Index weight on live graded quality ", weight.qual.kits, file =logfile, sep="\t")
     # old females    
     cat("\n", file = logfile)
     cat("Index weights on old females", file=logfile)
     cat("\n", file = logfile)
     cat("Index weight on litter size ", weight.fert.old.females, file =logfile, sep="\t")
     cat("\n", file = logfile)
     cat("Index weight on body weight ", weight.bw.old.females, file =logfile, sep="\t")
     cat("\n", file = logfile)
     cat("Index weight on live graded quality ", weight.qual.old.females, file =logfile, sep="\t")
     cat("\n", file = logfile)
   }
   
   cat("Controls for simulation", file = logfile, sep="\n")
   cat("Number of females on farm                    	/",n.females, file=logfile, sep="\t")
   cat("\n", file = logfile)
   cat("Number of replicates	                      /",nruns, file=logfile, sep="\t")
   cat("\n", file = logfile)
   cat("Number of years per replicates	               /",n, file=logfile, sep="\t")
   cat("\n", file = logfile)
   cat("Males per female                             	/",male.ratio, file=logfile, sep="\t")
   cat("\n", file = logfile)
   cat("Number of males	                               /",n.males, file=logfile, sep="\t")
   cat("\n", file = logfile)
   cat("Proportion of males barren                 	/",1-male.inf, file=logfile, sep="\t")
   cat("\n", file = logfile)
   cat("Proportion of females older than 1	           /",prop.oldfemales, file=logfile, sep="\t")
   cat("\n", file = logfile)
   cat("Max age of females	                           /",max.age.females, file=logfile, sep="\t")
   cat("\n", file = logfile)
   cat("Systematic crossmating, 0 = no, 1 = yes    	/",crossmating, file=logfile, sep="\t")
   cat("\n", file = logfile)
   cat("Systematic purebreeding, 0 = no, 1 = yes   	/",purebreeding, file=logfile, sep="\t")
   cat("\n", file = logfile)
   cat("Random culling ratio                      	/",cull.ratio, file=logfile, sep="\t")
   cat("\n", file = logfile)
   cat("Feed price per kg                            	/",feed.price, file=logfile, sep="\t")
   cat("\n", file = logfile)
   cat("Variable costs per skin                   	/",variable.costs, file=logfile, sep="\t")
   cat("\n", file = logfile)
   cat("Use true sire in pedigree, 0=no, 1 = yes    /",use.true.sire, file=logfile, sep="\t")
   cat("\n", file = logfile)
   cat("Use EBV to rank males in mating, 0=no, 1=yes/",use.blup.to.assort.mat, file=logfile, sep="\t")
   cat("\n", file = logfile)
   cat("Trace pedigree, 0=no, 1=yes                 /",trace.ped, file=logfile, sep="\t")
   cat("\n", file = logfile)
   cat("Proportion of males to deselect in 2nd mating	/",intensity.remating, file=logfile, sep="\t")
   cat("\n", file = logfile)
   
   
   close(con = logfile)
   
 }
 
 ################## Feed usage of adults #########################
 # this simple functions used a simple way to generate feeding costs for the breeders
 # males are considered to live for 110 days after pelting and females 365 days, use data from Tauson et al 2004
 # & mink nutrition book and assume K_lactation = 0.78 (like in sows)
 FeedUsageBreeders <- function (mating.list, next.gen.males, next.gen )   {
next.gen.males$feed.used.males <- (next.gen.males$phenotype.bw.oct^0.75*0.527*110)/5.0208
feed.used.males <- sum(next.gen.males$feed.used.males)
feed.used <- numeric(nrow(mating.list))
mating.list <- cbind(mating.list,feed.used)
transform(mating.list,feed.used, 13*obs_fert*0.78*4.53
                          +22.2*obs_fert*0.78*5.14+  
                            28.6*obs_fert*.78*5.86+
                            32.8*obs_fert*0.78*5.96)
feed.used.females <- sum(mating.list$feed.used)
next.gen$feed.used <- (next.gen$phenotype.bw.oct^0.75*0.527*365)/5.020
feed.used.females <- feed.used.females + sum(next.gen$feed.used)
feed.used.breeders <- (feed.used.females+feed.used.males)
return(feed.used.breeders)

} 
 
############### Pseudo imported males ##########################
 # Idea is to delete pedigree information about a certain proportion of 
 # males to use as a simple way to look at possible bias that comes about because
 # of missing pedigree information when importing from a better farm
 
 MaskPseudoImport <-
   function (pseudo.import,
             pseudo.import.prop,
             next.gen.males,
             selection.method,
             weight.bw.kits,
             weight.fert.kits,
             weight.qual.kits
) 
  {
     # browser()
     
         if (selection.method ==2  ) { # 2 == blup
           mean.SL <- mean(next.gen.males$skin.length.male)
           mean.SQ <- mean(next.gen.males$skin.qual)
           mean.LS <- mean(next.gen.males$litter.size)
           sd.SL <- sqrt(var(next.gen.males$skin.length.male))
           sd.SQ <- sqrt(var(next.gen.males$skin.qual))
           sd.LS <- sqrt(var(next.gen.males$litter.size))
           next.gen.males[, rank := 100 + (skin.length.male - mean.SL)/sd.SL*10 * weight.bw.kits +
                            (skin.qual-mean.SQ)/sd.SQ*10 * weight.qual.kits + (litter.size-mean.LS)/sd.LS*10 * weight.fert.kits]
           truncation <-
             quantile(next.gen.males$rank, probs = (1 - pseudo.import.prop))
           next.gen.males[, mask := ifelse (next.gen.males$rank < truncation, 0, 1)]
           next.gen.males$sire.assumed <-
             ifelse(next.gen.males$mask == 1, 0, next.gen.males$sire.assumed)
           next.gen.males$dam.id <-
             ifelse(next.gen.males$mask == 1, 0, next.gen.males$sire.assumed)
           set(next.gen.males, j = which(colnames(next.gen.males) %in%
                              c("mask", "rank"))  , value = NULL)
           return(next.gen.males)
    }
   }
 ####################### Sell of kits ##########################################
 # this function figures out if there are too many kits on the farm and then sells of
 # kits, picked at random until the number of kits fits with cage numbers 
 
 SellExtraKits <- function (kit.list, numb.of.old.females, n.cages) {
  # browser()
   if (nrow(kit.list)/2> (n.cages - numb.of.old.females) ) { 
     
     
     
     old.female.cages.needed <- numb.of.old.females # keeps track of how many cages are needed for old females
     
     
     sold.kit.numb <- (ceiling(nrow(kit.list)/2)-(n.cages-old.female.cages.needed) )*2
     sold.kit.id <- sample(kit.list$id, size = sold.kit.numb)
     
     sd <-
       setdiff(kit.list$id, sold.kit.id) # remove the sold kits from kit.list
     sd <- is.element(kit.list$id, sd)
     kit.list <- kit.list[sd, ] 
     
   }
   return(kit.list)
 }
 ################## Guess number of animals #############################
 NumberofBreeders <- function (fert.memory,n.cages,year) {
   # browser()
   if (year == 1 ) { 
     objective <- fert.memory[1]/2
   } else if (year > 1 ) {
       objective <- mean(fert.memory[year-1:year])*0.4
     }
   f.obj <- c(1,1)
   f.con <- matrix (c(objective,1,1,-1), nrow=2, byrow=TRUE)
   f.dir <- c("==","==")
   f.rhs <- c(n.cages,0)
   lp ("max", f.obj, f.con, f.dir, f.rhs)$solution
   n.females <- ceiling(lp ("max", f.obj, f.con, f.dir, f.rhs)$solution[1])
   return(n.females)
 }
 
 
 ############### Make phenotype for body weight ######
 MakePhenotypesBWMalesOct <- function(x,y,z,t,u,mat) {
   # x = mean, y = additive genetic, z = specific env, t = number of  sibs, mat =mating.list
   value <- x+ y + z + 
     rnorm( sum(mat$obs_fert))*(sqrt(bw.res.male))+
     t*sib.effect.male+ u * bw.eff.damage
   return(value)
 } 
 MakePhenotypesBWFemalesOct <- function(x,y,z,t,mat) {
   # x = mean, y = additive genetic, z = specific env, t = number of sibs, mat=mating.list
   value <- x+ y + z + 
     rnorm( sum(mat$obs_fert))*(sqrt(bw.res.female))+
     t*sib.effect.female 
   return(value)
 } 