#backup of mating function monster
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
