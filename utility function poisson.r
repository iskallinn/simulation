############### Creation of Gen0 females ############
# This function creates the base population of females
GenerateBaseFemales <- function () {
  id        <-  seq(1:n.females)
  #   fert      <-  rnorm(n.females)*sqrt(variance.fertility)
  add.gen <- as.data.table(rmvnorm(n.females, sigma = sigma,method="svd" ))
  colnames(add.gen)     <- c("bw.oct", "bw.sept", "litter.size", "live.qual",
                             "skin.qual", "skin.length" )
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
  
  gen0.females[, c("bw.oct")
               := lapply(.SD, function(x)
                 x * sqrt(var.bw.oct.female)), .SDcols = c("bw.oct")] # makes variance proper
  gen0.females[, c("bw.sept")
               := lapply(.SD, function(x)
                 x * sqrt(var.bw.sept.female)), .SDcols = c("bw.sept")] # makes variance proper
  
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
  add.gen <- as.data.table(rmvnorm(n.males, sigma = sigma,method="svd" ))
  colnames(add.gen)     <- c("bw.oct", "bw.sept", "litter.size", "live.qual",
                             "skin.qual", "skin.length" )
  sex                <-  numeric( n.males )  
  dam.id             <-  numeric( n.males )  
  sire.id            <-  numeric( n.males ) 
  birthyear          <-  numeric( n.males ) 
  can.remate         <-  rbinom(n.males, 1,0.8) 
  
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
             :=lapply(.SD, function(x) x*sqrt(variance.fertility)), .SDcols=c("litter.size")] # makes the variance proper
  
  gen0.males[, c("bw.oct")
             :=lapply(.SD, function(x) x*sqrt(var.bw.oct.male)), .SDcols=c("bw.oct")] # makes variance proper
  gen0.males[, c("bw.sept")
             :=lapply(.SD, function(x) x*sqrt(var.bw.sept.male)), .SDcols=c("bw.sept")] # makes variance proper
  
  
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
  mating.list <-
    x[rep(seq(nrow(x)), mating.willingness.1st),  #expands the male list into a mating list, based on mat.will.1st
      c(
        "id",
        "litter.size",
        "bw.oct",
        "bw.sept",
        "live.qual",
        "skin.qual",
        "skin.length",
        "semen.quality.1st",
        "semen.quality.2nd",
        "can.remate",
        "mating.willingness.1st",
        "mating.willingness.2nd"
      )
      , with = F] #specify which columns to incl.
  if (nrow(mating.list) > sum(y$mating.will.1st.round)) {
    mating.list <- mating.list[1:sum(y$mating.will.1st.round), ]
  }
  
  setnames(
    mating.list,
    c("id", "bw.oct", "litter.size","bw.sept", "live.qual", "skin.qual", "skin.length"),
    c("sire.id.1st", "sire.bw.oct.1st", "sire.fert.1st","sire.bw.sept.1st", 
      "sire.live.qual.1st", "sire.skin.qual.1st", "sire.skin.length.1st")
  )
  setkey(mating.list, sire.id.1st)
  mating.list[mating.list, c(
    "dam.id",
    "dam.fert",
    "f0",
    "obs_fert",
    "dam.bw.oct",
    "dam.bw.sept",
    "dam.skin.qual",
    "dam.skin.length",
    "dam.live.qual",
    "perm.env.ls",
    "birthyear.dam",
    "sire.id.2nd",
    "sire.bw.oct.2nd",
    "sire.bw.sept.2nd",
    "sire.fert.2nd",
    "sire.skin.qual.2nd",
    "sire.skin.length.2nd",
    "sire.live.qual.2nd"
  ) := 0]
  # moved scaling of perm.env.bw to the phenotype function later on
  mating.list[, `:=`(perm.env.bw.oct = rnorm(nrow(mating.list)), 
                     perm.env.bw.sept= rnorm(nrow(mating.list)))]
  # here I subset the dam list to throw away those who will not mate on first round
  y <- subset (y, mating.will.1st.round == 1)
  mating.list$dam.id       <- y[1:nrow(mating.list), .(id)]
  mating.list$dam.fert     <- y[1:nrow(mating.list), .(litter.size)]
  mating.list$dam.bw.oct   <- y[1:nrow(mating.list), .(bw.oct)]
  mating.list$perm.env.ls  <- y[1:nrow(mating.list), .(perm.env.ls)]
  mating.list$dam.bw.sept  <- y[1:nrow(mating.list), .(bw.sept)]
  mating.list$dam.skin.length  <- y[1:nrow(mating.list), .(skin.length)]
  mating.list$dam.skin.qual  <- y[1:nrow(mating.list), .(skin.qual)]
  mating.list$dam.live.qual  <- y[1:nrow(mating.list), .(live.qual)]
  
    
  if ("f0" %in% colnames(y)) {
    mating.list$f0  <- y[1:nrow(mating.list), .(f0)]
  }
  if ("perm.env.bs" %in% colnames(y)) {
    mating.list$perm.env.bs  <- y[1:nrow(mating.list), .(perm.env.bs)]
  }
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
  mating.list.allowed$semen.quality.2nd <-
    ifelse(
      mating.list.allowed$sire.id.1st == mating.list.allowed$sire.id.2nd,
      mating.list.allowed$semen.quality.1st,
      0
    )
  mating.list.allowed$sire.fert.2nd     <-
    ifelse(
      mating.list.allowed$sire.id.1st == mating.list.allowed$sire.id.2nd,
      mating.list.allowed$sire.fert.1st    ,
      0
    )
  mating.list.allowed$sire.bw.oct.2nd       <-
    ifelse(
      mating.list.allowed$sire.id.1st == mating.list.allowed$sire.id.2nd,
      mating.list.allowed$sire.bw.oct.1st      ,
      0
    )
  mating.list.allowed$sire.bw.sept.2nd       <-
    ifelse(
      mating.list.allowed$sire.id.1st == mating.list.allowed$sire.id.2nd,
      mating.list.allowed$sire.bw.sept.1st      ,
      0
    )
  mating.list.allowed$sire.skin.length.2nd       <-
    ifelse(
      mating.list.allowed$sire.id.1st == mating.list.allowed$sire.id.2nd,
      mating.list.allowed$sire.skin.length.1st      ,
      0
    )
  mating.list.allowed$sire.skin.qual.2nd       <-
    ifelse(
      mating.list.allowed$sire.id.1st == mating.list.allowed$sire.id.2nd,
      mating.list.allowed$sire.skin.qual.1st,
      0
    )
  mating.list.allowed$sire.live.qual.2nd       <-
    ifelse(
      mating.list.allowed$sire.id.1st == mating.list.allowed$sire.id.2nd,
      mating.list.allowed$sire.live.qual.1st      ,
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
      j = c("IDX"),
      value = NULL)
  mating.list.leftover <-
    rbind(mating.list.leftover, mating.list.notallowed)
  
  # here I should reorder the dams that are leftovers to prioritize younger dams and the ones mated with shitty males
  mating.list.leftover <-
    as.matrix(mating.list.leftover) 
  # way faster to loop through a matrix
  # this big loop has two modes, one if the selection method is blup and
  # another if the selection method is phenotypic
  # then each has three versions of the loop, one for year == 1, year ==2 and
  # year > 2. This is because of differing amounts of years in the matrix that
  # the loop accepts as input. Year 1 loops are always the same, regardless of 
  # selection method
  x <- as.matrix(x)
  if (selection.method == blup) {
    if (year == 1) {
      for (i in 1:nrow(x))  { #number of mating males
        #print(i)
        for (j in 1:x[[i, 17]])  { # number of matings left
          s <-
            sum(x[1:i, 17]) 
          # keeps track of how many females the male has been assigned
          #         print(s)
          #     print(j)
          if (s < nrow(mating.list.leftover)) {
            # this is not a solution, need to make another controlling mechanism if the mating willingness
            #exceeds the number of females to be mated
            mating.list.leftover[[s - (j - 1), 24]] <- x[[i, 1]]          # id of male
            mating.list.leftover[[s - (j - 1), 9]]  <- x[[i, 11]]          # semen.quality
            mating.list.leftover[[s - (j - 1), 27]] <- x[[i, 6]]          # fertility of male
            mating.list.leftover[[s - (j - 1), 25]] <- x[[i, 4]]          # body weight oct of male
            mating.list.leftover[[s - (j - 1), 26]] <- x[[i, 5]]          # body weight sept of male
            mating.list.leftover[[s - (j - 1), 29]] <- x[[i, 9]]          # skin.length of male
            mating.list.leftover[[s - (j - 1), 28]] <- x[[i, 8]]          # skin.qual sept of male
            mating.list.leftover[[s - (j - 1), 30]] <- x[[i, 7]]          # live.qual sept of male
            
            
            
          }
          if (s > nrow(mating.list.leftover)) {
            while (sum(x[1:(i - 1) , 17]) + j <= nrow(mating.list.leftover)) {
              # Here i solve the problem of if the males have more mating willingness than the number of females
              mating.list.leftover[[s - (j - 1), 24]] <- x[[i, 1]]          # id of male
              mating.list.leftover[[s - (j - 1), 9]]  <- x[[i, 11]]          # semen.quality
              mating.list.leftover[[s - (j - 1), 27]] <- x[[i, 6]]          # fertility of male
              mating.list.leftover[[s - (j - 1), 25]] <- x[[i, 4]]          # body weight oct of male
              mating.list.leftover[[s - (j - 1), 26]] <- x[[i, 5]]          # body weight sept of male
              mating.list.leftover[[s - (j - 1), 29]] <- x[[i, 9]]          # skin.length of male
              mating.list.leftover[[s - (j - 1), 28]] <- x[[i, 8]]          # skin.qual sept of male
              mating.list.leftover[[s - (j - 1), 30]] <- x[[i, 7]]          # live.qual sept of male
              
              break
            }
          }
        }
      }
    }
    if (year == 2) {
      for (i in 1:nrow(x))  {
        #print(i)
        for (j in 1:x[[i, 18]])  {
          s <-
            sum(x[1:i, 18]) # keeps track of how many females the male has been assigned
          #         print(s)
          #     print(j)
          if (s < nrow(mating.list.leftover)) {
            # this is not a solution, need to make another controlling mechanism if the mating willingness
            #exceeds the number of females to be mated
            mating.list.leftover       [[s - (j - 1), 16]] <-
              x[[i, 1]]          # id of male
            mating.list.leftover       [[s - (j - 1), 5]]  <-
              x[[i, 17]]          # semen.quality
            mating.list.leftover       [[s - (j - 1), 18]] <-
              x[[i, 10]]          # fertility of male
            mating.list.leftover       [[s - (j - 1), 17]] <-
              x[[i, 11]]          # body size of male
          }
          if (s > nrow(mating.list.leftover)) {
            while (sum(x[[1:(i - 1) , 18]]) + j <= nrow(mating.list.leftover)) {
              # Here i solve the problem of if the males have more mating willingness than the number of females
              mating.list.leftover       [[sum (x$matings.left [1:(i - 1)]) + j, 16]] <-
                x[[i, 1]]          # id of male
              mating.list.leftover       [[sum (x$matings.left [1:(i - 1)]) + j,  5]] <-
                x[[i, 17]]         # semen.quality
              mating.list.leftover       [[sum (x$matings.left [1:(i - 1)]) + j, 18]] <-
                x[[i, 10]]         # fertility of male
              mating.list.leftover       [[sum (x$matings.left [1:(i - 1)]) + j, 17]] <-
                x[[i, 11]]         # body size of male
              
              break
            }
          }
        }
      }
    }
    
    if (year > 2) {
      for (i in 1:nrow(x))  {
        #print(i)
        for (j in 1:x[[i, 25]])  {
          s <-
            sum(x[1:i, 25]) # keeps track of how many females the male has been assigned
          #         print(s)
          #     print(j)
          if (s < nrow(mating.list.leftover)) {
            # this is not a solution, need to make another controlling mechanism if the mating willingness
            #exceeds the number of females to be mated
            mating.list.leftover       [[s - (j - 1), 16]] <-
              x[[i, 1]]          # id of male
            mating.list.leftover       [[s - (j - 1), 5]]  <-
              x[[i, 24]]          # semen.quality
            mating.list.leftover       [[s - (j - 1), 18]] <-
              x[[i, 10]]          # fertility of male
            mating.list.leftover       [[s - (j - 1), 17]] <-
              x[[i, 11]]          # body size of male
          }
          if (s > nrow(mating.list.leftover)) {
            while (sum(x[[1:(i - 1) , 25]]) + j <= nrow(mating.list.leftover)) {
              # Here i solve the problem of if the males have more mating willingness than the number of females
              mating.list.leftover       [[sum (x$matings.left [1:(i - 1)]) + j, 16]] <-
                x[[i, 1]]          # id of male
              mating.list.leftover       [[sum (x$matings.left [1:(i - 1)]) + j,  5]] <-
                x[[i, 24]]         # semen.quality
              mating.list.leftover       [[sum (x$matings.left [1:(i - 1)]) + j, 18]] <-
                x[[i, 10]]         # fertility of male
              mating.list.leftover       [[sum (x$matings.left [1:(i - 1)]) + j, 17]] <-
                x[[i, 11]]         # body size of male
              
              break
            }
          }
        }
      }
    }
  } else if (selection.method == phenotypic) {
    if (year == 1) {
      for (i in 1:nrow(x))  {
        #print(i)
        for (j in 1:x[[i, 17]])  {
          s <-
            sum(x[1:i, 17]) # keeps track of how many females the male has been assigned
          #         print(s)
          #     print(j)
          if (s < nrow(mating.list.leftover)) {
            # this is not a solution, need to make another controlling mechanism if the mating willingness
            #exceeds the number of females to be mated
            mating.list.leftover[[s - (j - 1), 24]] <- x[[i, 1]]          # id of male
            mating.list.leftover[[s - (j - 1), 9]]  <- x[[i, 11]]          # semen.quality
            mating.list.leftover[[s - (j - 1), 27]] <- x[[i, 6]]          # fertility of male
            mating.list.leftover[[s - (j - 1), 25]] <- x[[i, 4]]          # body weight oct of male
            mating.list.leftover[[s - (j - 1), 26]] <- x[[i, 5]]          # body weight sept of male
            mating.list.leftover[[s - (j - 1), 29]] <- x[[i, 9]]          # skin.length of male
            mating.list.leftover[[s - (j - 1), 28]] <- x[[i, 8]]          # skin.qual sept of male
            mating.list.leftover[[s - (j - 1), 30]] <- x[[i, 7]]          # live.qual sept of male
          }
          if (s > nrow(mating.list.leftover)) {
            while (sum(x[1:(i - 1) , 17]) + j <= nrow(mating.list.leftover)) {
              # Here i solve the problem of if the males have more mating willingness than the number of females
              mating.list.leftover[[s - (j - 1), 24]] <- x[[i, 1]]          # id of male
              mating.list.leftover[[s - (j - 1), 9]]  <- x[[i, 11]]          # semen.quality
              mating.list.leftover[[s - (j - 1), 27]] <- x[[i, 6]]          # fertility of male
              mating.list.leftover[[s - (j - 1), 25]] <- x[[i, 4]]          # body weight oct of male
              mating.list.leftover[[s - (j - 1), 26]] <- x[[i, 5]]          # body weight sept of male
              mating.list.leftover[[s - (j - 1), 29]] <- x[[i, 9]]          # skin.length of male
              mating.list.leftover[[s - (j - 1), 28]] <- x[[i, 8]]          # skin.qual sept of male
              mating.list.leftover[[s - (j - 1), 30]] <- x[[i, 7]]          # live.qual sept of male
              
              break
            }
          }
        }
      }
    }
    if (year > 1) {
      for (i in 1:nrow(x))  {
        #print(i)
        for (j in 1:x[[i, 18]])  {
          s <-
            sum(x[1:i, 18]) # keeps track of how many females the male has been assigned
          #         print(s)
          #     print(j)
          if (s < nrow(mating.list.leftover)) {
            # this is not a solution, need to make another controlling mechanism if the mating willingness
            #exceeds the number of females to be mated
            mating.list.leftover       [[s - (j - 1), 16]] <-
              x[[i, 1]]          # id of male
            mating.list.leftover       [[s - (j - 1), 5]]  <-
              x[[i, 17]]          # semen.quality
            mating.list.leftover       [[s - (j - 1), 18]] <-
              x[[i, 10]]          # fertility of male
            mating.list.leftover       [[s - (j - 1), 17]] <-
              x[[i, 11]]          # body size of male
          }
          if (s > nrow(mating.list.leftover)) {
            while (sum(x[[1:(i - 1) , 18]]) + j <= nrow(mating.list.leftover)) {
              # Here i solve the problem of if the males have more mating willingness than the number of females
              mating.list.leftover       [[sum (x$matings.left [1:(i - 1)]) + j, 16]] <-
                x[[i, 1]]          # id of male
              mating.list.leftover       [[sum (x$matings.left [1:(i - 1)]) + j,  5]] <-
                x[[i, 17]]         # semen.quality
              mating.list.leftover       [[sum (x$matings.left [1:(i - 1)]) + j, 18]] <-
                x[[i, 10]]         # fertility of male
              mating.list.leftover       [[sum (x$matings.left [1:(i - 1)]) + j, 17]] <-
                x[[i, 11]]         # body size of male
              
              break
            }
          }
        }
      }
    }
  }
  mating.list.leftover <- as.data.table(mating.list.leftover)
  x <- as.data.table(x)
  set(mating.list.remated, j = c("IDX"), value = NULL) # need to remove the counter to bind the mated ones to the rest
  
  mating.list <-
    rbind(mating.list.remated, mating.list.leftover) # merge the now completed mating list together
  mating.list[, `:=`(
    semen.quality = ifelse (
      mating.list$semen.quality.2nd == 0,
      0,
      mating.list$semen.quality.2nd
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
      "semen.quality.2nd"
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
  gen1 <- x[rep(seq(nrow(x)), obs_fert), 
            c("dam.id",
              "sire.id.1st",
              "dam.fert",
              "sire.fert.1st",
              "sire.bw.oct.1st",
              "sire.bw.sept.1st",
              "sire.skin.length.1st",
              "sire.skin.qual.1st",
              "sire.live.qual.1st",
              "f0",
              "obs_fert",
              "dam.bw.oct",
              "dam.bw.sept",
              "dam.skin.length",
              "dam.skin.qual",
              "dam.live.qual",
              "sire.bw.oct.1st",
              "sire.bw.sept.1st",
              "perm.env.bw.oct",
              "perm.env.bw.sept",
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
  id <- seq(1:sum(x$obs_fert)) + max(y$id) # makes ID
  birthyear <- rep (1, sum(x$obs_fert)) # makes birthyear
  sex <- rbinom(sum(x$obs_fert),1,0.5)+1 # makes sex, TODO check if this is true
  true.sire <- numeric(nrow(gen1))
  true.sire.check <- numeric(nrow(gen1))
  mendelian <- as.data.table(rmvnorm(sum(x$obs_fert),sigma=sigma,method="svd"))
  colnames(mendelian) <- c("mend.bw.oct", "mend.bw.sept", "mend.litter.size", "mend.live.qual",
                           "mend.skin.qual", "mend.skin.length" )
  
  gen1 <- cbind (id,gen1,  birthyear, sex, mendelian, true.sire,true.sire.check) # binds id, sex and birthyear to data.table
  
  gen1$true.sire.check <- ifelse( gen1$sire.id.1st != gen1$sire.id.2nd, rbinom(nrow(gen1), 1, 0.85), 1) # 85% chance that the kits are sired by 2nd mating
  gen1$true.sire <- ifelse( gen1$true.sire.check == 0, gen1$sire.id.1st, gen1$sire.id.2nd)
  gen1[, `:=`(true.sire.fert = 
                ifelse(gen1$true.sire == gen1$sire.id.2nd, 
                       gen1$sire.fert.2nd, gen1$sire.fert.1st), 
              true.sire.bw.oct = 
                ifelse(gen1$true.sire == gen1$sire.id.2nd, 
                       gen1$sire.bw.oct.2nd, gen1$sire.bw.oct.1st),
              true.sire.bw.sept = 
                ifelse(gen1$true.sire == gen1$sire.id.2nd, 
                       gen1$sire.bw.sept.2nd, gen1$sire.bw.sept.1st),
              true.sire.skin.length = 
                ifelse(gen1$true.sire == gen1$sire.id.2nd, 
                       gen1$sire.skin.length.2nd, gen1$sire.skin.length.1st),
              true.sire.skin.qual = 
                ifelse(gen1$true.sire == gen1$sire.id.2nd, 
                       gen1$sire.skin.qual.2nd, gen1$sire.skin.qual.1st),
              true.sire.live.qual = 
                ifelse(gen1$true.sire == gen1$sire.id.2nd, 
                       gen1$sire.live.qual.2nd, gen1$sire.live.qual.1st)
  )]
  gen1$mend.bw.oct <-  ifelse(gen1$sex==1, gen1$mend.bw.oct*sqrt(0.5*var.bw.oct.male),
                              gen1$mend.bw.oct*sqrt(0.5*var.bw.oct.female))
  gen1$mend.bw.sept <-  ifelse(gen1$sex==1, gen1$mend.bw.sept*sqrt(0.5*var.bw.sept.male),
                               gen1$mend.bw.oct*sqrt(0.5*var.bw.sept.female))
  
  gen1[, `:=`(litter.size = 0.5*(dam.fert + true.sire.fert) + 
                mend.litter.size*(sqrt(0.5*variance.fertility)) # Breeding value of offspring, littersize
              , perm.env.ls = rnorm(sum(x$obs_fert))*sqrt(var.perm.env.ls) # perm env for litter size
              , bw.oct = 0.5*(dam.bw.oct + true.sire.bw.oct) + 
                mend.bw.oct,
              bw.sept = 0.5*(dam.bw.sept + true.sire.bw.sept) + 
                mend.bw.oct)]# Breeding value of offspring, body size
  #  gen1 <- count.sex.siblings(gen1) # calls function to count offspring NOT NEEDED ATM
  
  setnames(gen1, c("obs_fert","sire.id.2nd"), c("own_littersize","sire.assumed")) # renames obs_fert to own littersize of kits
  gen1$dam.age <- ifelse( z - gen1$birthyear.dam > 1, 1,0 )
  
  gen1$phenotype.bw.oct <- ifelse( gen1$sex == 1,MakePhenotypesBWMalesOct(mean.body.size.male.oct , gen1$bw.oct , gen1$perm.env.bw.oct , gen1$own_littersize, gen1$dam.age,x  ) 
                                   , MakePhenotypesBWFemalesOct(mean.body.size.female.oct , gen1$bw.oct , gen1$perm.env.bw.oct,gen1$own_littersize,x )) 
  gen1$phenotype.bw.sept <- ifelse( gen1$sex == 1,
                                    MakePhenotypesBWMalesSept(mean.body.size.male.sept , gen1$bw.sept , gen1$perm.env.bw.sept , gen1$own_littersize, gen1$dam.age,x  ) 
                                    , MakePhenotypesBWFemalesOct(mean.body.size.female.sept , gen1$bw.sept , gen1$perm.env.bw.sept,gen1$own_littersize,x )) 
  
  # generate phenotype for body size
  
  gen1[,`:=`(perm.env.bw.sept = ifelse(sex==1, rnorm(sum(x$obs_fert))*sqrt(var.c.bw.sept.male),
                                       rnorm(sum(x$obs_fert))*sqrt(var.c.bw.sept.female)),
             perm.env.bw.oct = ifelse(sex==1, rnorm(sum(x$obs_fert))*sqrt(var.c.bw.oct.male),
                                      rnorm(sum(x$obs_fert))*sqrt(var.c.bw.oct.female)))] # generate specific env. for body size
  
  set( gen1, j=which(colnames(gen1) %in% c(
    "sire.fert.1st",
    "sire.bw.oct.1st",
    "sire.bw.sept.1st",
    "sire.skin.length.1st",
    "sire.skin.qual.1st",
    "sire.live.qual.1st",
    "dam.fert",
    "dam.bw.oct",
    "dam.bw.sept",
    "dam.skin.length",
    "dam.skin.qual",
    "dam.live.qual",
    "sire.bw.oct.1st",
    "sire.bw.sept.1st",
    "perm.env.bw.oct",
    "perm.env.bw.sept",
    "birthyear.dam",
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
  )) , value=NULL ) # removes bv of parents
  return(gen1)
}
############### Phenotypic Selection of old females ###############################################
# currently they are only truncated on their litter size phenotype
PhenoSelectionOldFemales <- function (y,x, year) { # y = gen0.females, x = x
    setkey(x, obs_fert)
    setorder(x, -obs_fert)
    if ("birthyear.dam" %in% colnames(x)) {
      x <- subset(x, year -birthyear.dam < max.age.females) 
    }
    old.females <- x[1:(n.females*prop.oldfemales),]
    set( old.females, j=which(colnames(old.females) %in% 
                                c("sire.id.1st","barren","dam.fert","sire.fert.1st","barren","sire.bs.1st", "dam.bs",
                                  "sire.id.2nd","sire.fert.2nd","sire.bs.2nd", "mating.will.1st.round",
                                  "mating.will.2nd.round"))  , value=NULL )
    setnames(old.females, "dam.id", "id")
    old.females <- merge(old.females, y, by="id")
    if("birthyear.dam" %in% colnames(old.females)) {
      set( old.females, j=which(colnames(old.females) %in% 
                                  c("birthyear.dam"))  , value=NULL )
    }
    if("f0.dam" %in% colnames(old.females)) {
      set( old.females, j=which(colnames(old.females) %in% 
                                  "f0.dam")  , value=NULL )
    }
    if("obs_fert.y" %in% colnames(old.females)) {
      set( old.females, j=which(colnames(old.females) %in% 
                                  "obs_fert.y")  , value=NULL )
    }
    if("obs_fert.x" %in% colnames(old.females)) {
      setnames(old.females, "obs_fert.x", "obs_fert")
    }
    if("obs_fert.x" %in% colnames(old.females)) {
    }  else {
      old.females[old.females,own_littersize:=0]   
    }
    if("perm.env.ls.x" %in% colnames(old.females)) {
      setnames(old.females, "perm.env.ls.x", "perm.env.ls")
      set( old.females, j=which(colnames(old.females) %in% 
                                  "perm.env.ls.y")  , value=NULL )
    }
    if("perm.env.bs.x" %in% colnames(old.females)) {
      setnames(old.females, "perm.env.bs.x", "perm.env.bs")
      set( old.females, j=which(colnames(old.females) %in% 
                                  "perm.env.bs.y")  , value=NULL )
    }
    if("phenotype.bw.oct" %in% colnames(old.females)){
      
    } else {
      old.females <-transform(old.females, phenotype.bw.oct = mean.body.size.female + bw.oct + rnorm(1)*sqrt(var.res.body.size))
    }
    if("mating.will.2nd.round" %in% colnames(old.females)) {
      set(old.females, j=c("mating.will.2nd.round", "mating.will.1st.round"),value=NULL)
    }
    if("sire.id" %in% colnames(old.females)) {
      setnames(old.females, "sire.id", "sire.assumed")
      old.females[,`:=`(true.sire = old.females$sire.assumed)]
    }
    if("dam.age" %in% colnames(old.females)) {
      set(old.females, j="dam.age", value = NULL)
    }
    if("dam.age.x" %in% colnames(old.females)) {
      set(old.females, j=c("dam.age.x","dam.age.y","prodyear"), value = NULL)
      }
    
    return(old.females)
    
    
  }  
############### Phenotypic Selection of yearling females in 1st gen ############
PhenoSelectionFemaleKits <- function (x,y) { # x = kit.list, y = y  
    truncation.point <-  quantile( x$own_littersize,  probs =  quantile.setting ) 
    selection.candidates.females <- subset(x, own_littersize >= truncation.point) # throw away the smallest litters
    selection.candidates.females <-  subset( selection.candidates.females,  sex  ==   2) # take the female kits
    # selection.candidates.females <- subset( selection.candidates.females, phenotype.bw.oct < roof.body.size)
    # if (nrow(selection.candidates.females) < n.females*(1-prop.oldfemales)){ 
    #   truncation.point <-  quantile( x$own_littersize,  probs =  (quantile.setting  ) ) #can't change this since the kits won't have cards 
    #   selection.candidates.females <- subset(x, own_littersize >= truncation.point) # throw away the smallest litters
    #   selection.candidates.females <-  subset( selection.candidates.females,  sex  ==   2) # take the female kits
    #   selection.candidates.females <- subset( selection.candidates.females, phenotype.bw.oct < (roof.body.size+200)) #ease restrictions on size
    #   
    # }
    setkey(selection.candidates.females, phenotype.bw.oct)           # in order to speed up ordering 
    setorder(selection.candidates.females,-phenotype.bw.oct)         # order the kits according to body size 
    next.gen <- selection.candidates.females[1:(n.females-nrow(y)),]              # take the biggest kits
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
PhenoSelectionMaleKits <- function (x) {
    truncation.point <-  quantile( x$own_littersize,  probs =  quantile.setting ) 
    selection.candidates.males <- subset(x, own_littersize >= truncation.point) # throw away the smallest litters
    selection.candidates.males <-  subset( selection.candidates.males,  sex  ==   1  ) # take the female kits
    setkey(selection.candidates.males, phenotype.bw.oct)           # in order to speed up ordering 
    setorder(selection.candidates.males,-phenotype.bw.oct)         # order the kits according to litter size 
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
IndSelectionOldFemales <- function (x,y,z,year) { # x = next.gen, y = solutions, z = solutions.bw
   
   if("blup.fert" %in% colnames(x)) {
     set( x, j=which(colnames(x) %in% 
                              c("blup.fert", "sem.blup.fert"))  , value=NULL )
   }
  if("blup.bwnov" %in% colnames(x)) {
    set( x, j=which(colnames(x) %in% 
                             c("blup.bwnov", "sep.blup.bwnov"))  , value=NULL )
  }
  
   
   x <- merge(x, y, by ="id", all.x=TRUE) # merge solutions from blup to data.table containing the old females
   setkey(x, blup.fert)
   x <-merge(x, z, by="id", all.x=TRUE) # merge solutions from blup to data.table containing the old females
   x[, `:=`(index.bw = 100+ (blup.bwnov-mean(x$blup.bwnov))/(sqrt(var(x$blup.bwnov)))*10,
                   index.fert =100+ (blup.fert-mean(x$blup.fert))/(sqrt(var(x$blup.fert)))*10)]
   x <- transform(x, comb.ind = index.bw*weight.bw.old.females+index.fert*weight.fert.old.females)
   x <- subset(x, year -birthyear < max.age.females) 
   setkey(x, comb.ind)
      setorder(x, -comb.ind)
   #then I will need to make rangering
   old.females <- x[1:(n.females*prop.oldfemales),]
   set( old.females, j=which(colnames(old.females) %in% 
                               c("sire.id.1st","barren","dam.fert","sire.fert.1st","barren","sire.bs.1st", "dam.bs",
                                 "sire.id.2nd","sire.fert.2nd","sire.bs.2nd", "mating.will.1st.round",
                                 "mating.will.2nd.round"))  , value=NULL )
   if("dam.age" %in% colnames(old.females)) {
     set( old.females, j=which(colnames(old.females) %in% 
                                 c("dam.age"))  , value=NULL )
   }
   return(old.females)
 }  
 
############### Index selection for females ################
IndSelFemaleKits <- function (x,y,z,q) { # x = kit.list , y = solutions.fert, z = solutions.bw, q = old.females 
   x <- merge(x, y, by= "id", all.x=TRUE) # merge to solutions of blup of fertility
   x <- merge(x, z, by="id", all.x=TRUE)
   x[, `:=`(index.bw = 100+ (blup.bwnov-mean(x$blup.bwnov))/(sqrt(var(x$blup.bwnov)))*10,
                   index.fert =100+ (blup.fert-mean(x$blup.fert))/(sqrt(var(x$blup.fert)))*10)]
   x <- transform(x, comb.ind = index.bw*weight.bw.kits+index.fert*weight.fert.kits)
   
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
 IndSelMaleKits <- function (x,y,z) { # x = kit.list, y = solutions, z = solutions.bw
   x <- merge(x, y, by= "id", all.x=TRUE) # merge to solutions of blup of fertility
   x <- merge(x, z, by="id", all.x=TRUE)
   x[, `:=`(index.bw = 100+ (blup.bwnov-mean(x$blup.bwnov))/(sqrt(var(x$blup.bwnov)))*10,
            index.fert =100+ (blup.fert-mean(x$blup.fert))/(sqrt(var(x$blup.fert)))*10)]
   x <- transform(x, comb.ind = index.bw*weight.bw.kits+index.fert*weight.fert.kits)
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
PrepareMalesForMating <- function (x) { # x = next.gen.males 
  setkey(x, id)
  #   next.gen.males[next.gen.males,semen.quality:=rbinom( nrow(next.gen.males),  1,  male.inf )]
  #   next.gen.males[next.gen.males,mating.willingness:=rZIP( nrow(next.gen.males),  mu = male.ratio,  sigma = 0.05 )]
  x[,`:=`( semen.quality.1st=  rbinom( nrow(x),  1,  male.inf ),  
                        mating.willingness.1st = rZIP( nrow(x),  mu = male.ratio,  sigma = 0.05 ),
                        mating.willingness.2nd = rZIP( nrow(x), mu = male.ratio, sigma = 0.05),
                        can.remate = rbinom(nrow(x), 1,0.8))]
  x[,`:=`(semen.quality.2nd = semen.quality.1st)]
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
                          c("dam.id","sire.id.1st","dam.fert","sire.fert.1st","f0.dam","obs_fert","dam.bs","sire.bs.1st"
                            ,"perm.env.bs","birthyear.dam","sire.id.2nd","sire.fert.2nd","sire.bs.2nd") , with=F]
 
  id <- seq(1:sum(x$obs_fert)) + max(z$id) # makes ID
  
  birthyear <- rep (year, sum(x$obs_fert)) # makes birthyear
  sex <- rbinom(sum(x$obs_fert),1,0.5)+1 # makes sex, TODO check if this is true
  true.sire <- numeric(nrow(kit.list))
  true.sire.check <- numeric(nrow(kit.list))
  mendelian <- as.data.table(rmvnorm(nrow(kit.list),sigma=sigma))
  colnames(mendelian) <- c("mend.fert","mend.bw")
  
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
  kit.list$sire.bs.2nd <-
    ifelse(
      kit.list$sire.id.2nd == kit.list$sire.id.1st,
      kit.list$sire.bs.1st,
      kit.list$sire.bs.2nd
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
    true.sire.bs = ifelse(
      kit.list$true.sire == kit.list$sire.id.2nd,
      kit.list$sire.bs.2nd,
      kit.list$sire.bs.1st
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
  kit.list[, `:=`(fert = 0.5*(dam.fert + true.sire.fert) + 
                    mend.fert*(sqrt(0.5*variance.fertility)*( 1- f0 )) # Breeding value of offspring, littersize
                  , perm.env.ls = rnorm(sum(x$obs_fert))*sqrt(var.perm.env.ls) # perm env for litter size
                  , bw.oct = 0.5*(dam.bs + true.sire.bs) + 
                    mend.bw*(sqrt(0.5*var.bw.oct)*( 1- f0 )))]# Breeding value of offspring, body size
  
  setnames(kit.list, "obs_fert", "own_littersize") # changes obs fert into the littersize of the kit
  kit.list$dam.age <- ifelse( year - kit.list$birthyear.dam > 1, 1,0 )
  kit.list$phenotype.bw.oct <- ifelse( kit.list$sex == 1,MakePhenotypesBWMales ( mean.body.size.male ,
  kit.list$bw.oct , kit.list$perm.env.bs , kit.list$own_littersize, kit.list$dam.age ,x ) ,
  MakePhenotypesBWFemales(mean.body.size.female ,kit.list$bw.oct , 
  kit.list$perm.env.bs, kit.list$own_littersize,x )) 
  # generate phenotype for body size
  
  kit.list[,`:=`(perm.env.bs = rnorm(sum(x$obs_fert))*sqrt(var.body.size.spec.env))] # generate specific env. for body size
  
  set( kit.list, j=which(colnames(kit.list) %in% c("dam.fert","sire.fert","sire.bs","dam.bs","mend.fert"
                                                   ,"mend.bw","birthyear.dam","dam.age","sire.id.1st","sire.bs.1st",
                                                   "sire.fert.1st","sire.fert.2nd","sire.bs.2nd"
                                                   ,"true.sire.fert","true.sire.bs", "true.sire.check")) , value=NULL ) # removes bv of parents
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
CalulateBLUPBodyWeightNov <- function () {
# run reml and then use the output as priors in BLUP
  # system2("C:/Users/Notandi/Dropbox/Projects/simulation of mink farm/Output/DMU analysis/run_dmuai.bat", " reml_bwnov")
  # run blup on BW  
system2("C:/Users/Notandi/Dropbox/Projects/simulation of mink farm/Output/DMU analysis/run_dmu4.bat", " bw_nov")
# read the solutions and only keep the predictions of BV (they're not that right)
  #  if ( !file.exists('fort.70')) {
  # now <-as.numeric(Sys.time()); howlong<-1; delt<-0; while(delt < howlong) 
  #   { delt<-as.numeric(Sys.time())-now }
  # }
  solutions.bw <- as.matrix(read.table(file="C:/Users/Notandi/Dropbox/Projects/simulation of mink farm/Output/DMU analysis/bw_nov.SOL"))
  # }
solutions.bw <- as.data.table(solutions.bw)
solutions.bw <- subset(solutions.bw, V1 == 4 ) # throw away the estimation of permanent environment
set (solutions.bw, j=c("V1","V2","V3", "V4", "V6", "V7"), value= NULL)
setnames(solutions.bw, c("V5","V8", "V9"),c("id", "blup.bwnov", "sep.blup.bwnov"))
return (solutions.bw)
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
   dirfile[8] <- c(paste("$DATA  ASCII (4,1,-9999) Phenotypes",p, sep="")) 
   # change the input file for BLUP so it uses the next outputfile
   writeLines(dirfile, "bw_nov.DIR")
   dirfile[25] <- c(paste("$VAR_STR 2 PED 2 ASCII Big_pedigree_",p, sep="")) 
   # change the input file for BLUP so it uses the next pedigree
   writeLines(dirfile,"bw_nov.DIR")
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
    
   }
   
    }
 
############### Create observation file for phenotypes #################
# currently only size
WriteObservationFileBodyWeight <- function (x,year,p,solutions) {
  x[, c("id","dam.id","own_littersize","sex")
              :=lapply(.SD, function(x) as.integer(x)), .SDcols=c("id","dam.id","own_littersize","sex")]
  # x[,`:=`(prodyear=as.integer(ifelse( sex == 1, year*(sex+6),year*(sex+8))))]
  # NOTE if simulation runs to 81 years, then year*sex will be the same for year 63 for dams and 81 for males
  if(mask.phenotypes == 1 & year > 1) {
    x <- merge(x, solutions, by= "id", all.x=TRUE) # merge to solutions of blup of fertility
    truncation.point <-  quantile( x$blup.fert,  probs =  quantile.setting ) 
    x <- subset(x, blup.fert >= truncation.point, select= c("id","dam.id","own_littersize","sex","phenotype.bw.oct"))
  }
  if (year  == 1 ) {phenotypes <- file(description = paste("Phenotypes",p, sep=""), open="w")
  write.table(format(x[,.(id,dam.id,sex,own_littersize,phenotype.bw.oct)], nsmall=1, digits=2), file= phenotypes, append= TRUE,col.names = FALSE, row.names = FALSE, quote = FALSE)
 close(con=phenotypes)
   } else if (year > 1) {
     phenotypes <- file(description = paste("Phenotypes",p, sep=""), open="a")
      write.table(format(x[,.(id,dam.id,sex,own_littersize,phenotype.bw.oct)], nsmall=1, digits=2), file= phenotypes, append= TRUE,col.names = FALSE, row.names = FALSE, quote = FALSE)
      close(con=phenotypes)
    }
  }
############### Make gen1 big.pedfile
WriteBigPedigree <- function (x,y,year,p) { # x = gen1 or kit.list , y = pedfile, year = year
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
