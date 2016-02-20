############### Creation of Gen0 females ############

# create a list for the females
# create ID for females
generate.base.females <- function () {
  id        <-  seq(1:n.females)
  #   fert      <-  rnorm(n.females)*sqrt(variance.fertility)
  add.gen <- as.data.table(rmvnorm(n.females,sigma=sigma))
  colnames(add.gen) <- c("fert","direct.genetic.body.size")
  perm.env.ls <- rnorm(n.females)*sqrt(var.perm.env.ls)
  sex       <-  rep(2, times = n.females)
  dam.id    <-  rep(0, times = n.females) 
  sire.id   <-  rep(0, times = n.females)
  birthyear <-  rep(0, times = n.females)
  mating.will.1st.round <- rbinom( n.females, 1, mating.will.yearling.1st)
  mating.will.2nd.round <- numeric(n.females)
  #create breeding value of fertility for females
  
  gen0.females <-  data.table( id, add.gen,perm.env.ls, sex, sire.id
                               ,  dam.id, birthyear, mating.will.1st.round, mating.will.2nd.round )
  gen0.females$mating.will.2nd.round <- ifelse (gen0.females$mating.will.1st.round == 1
                                                , rbinom(sum(gen0.females$mating.will.1st.round),1,mating.will.yearling.2nd)
                                                ,0 ) 
  
   
  gen0.females[, c("fert")
               :=lapply(.SD, function(x) x*sqrt(variance.fertility)), .SDcols=c("fert")] # makes the variance proper
  
  gen0.females[, c("direct.genetic.body.size")
               :=lapply(.SD, function(x) x*sqrt(var.body.size.direct)), .SDcols=c("direct.genetic.body.size")] # makes variance proper
  
  
  return (gen0.females)
}
############### Creation of Gen0 males#######################################
# create a list for the males
generate.base.males <- function () {
  mating.willingness.1st <-  numeric( n.males )  # preallocate a numeric vector
  mating.willingness.2nd  <-  numeric( n.males)
  semen.quality.1st      <-  numeric( n.males )  # preallocate a numeric vector
  semen.quality.2nd      <-  numeric( n.males )  # preallocate a numeric vector
  id                 <-  numeric( n.males )  # preallocate a numeric vector
  #   fert               <-  numeric( n.males )  # preallocate a numeric vector
  add.gen <- as.data.table(rmvnorm(n.males,sigma=sigma))
  colnames(add.gen) <- c("fert","direct.genetic.body.size")
  sex                <-  numeric( n.males )  # preallocate a numeric vector
  dam.id             <-  numeric( n.males ) # preallocate a numeric vector 
  sire.id            <-  numeric( n.males ) # preallocate a numeric vector
  birthyear          <-  numeric( n.males ) # preallocate a numeric vector
  can.remate         <-  rbinom(n.males, 1,0.8) #TODO make this into something more meaningful once I have quality or something to rank the males on
  
  for ( i in 1:n.males )  {
    mating.willingness.1st[i]     <-  rZIP( 1,  mu = male.ratio,  sigma = 0.05 )               # this is guesswork
    mating.willingness.2nd [i]     <-  rZIP( 1,  mu = 8,  sigma = 0.05 )                       # too high, since this will make sure all females are mated 2
    id[i]       <-  n.females + i                                             # create ID for males
    #     fert[i]     <-  rnorm( 1 )*sqrt(variance.fertility)  # create breeding value of fertility for males
    semen.quality.1st[i] <-  rbinom( 1,  1,  male.inf )                           # create on off for barren males
    semen.quality.2nd [i] <- semen.quality.1st [i]
    sex [i]         <-  1                                                         # 1 = male, 2 = female
    dam.id [i]     <-  0                                                         # this is gen0 so unknown parents
    sire.id [i]   <-  0                                                         # this is gen0 so unknown parents
  }
  
  # make data table out of the males 
  gen0.males <-  data.table( id, mating.willingness.1st,mating.willingness.2nd, add.gen
                             , semen.quality.1st,semen.quality.2nd, sex, sire.id, dam.id, birthyear,can.remate ) 
  
  gen0.males[, c("fert")
             :=lapply(.SD, function(x) x*sqrt(variance.fertility)), .SDcols=c("fert")] # makes the variance proper
  
  gen0.males[, c("direct.genetic.body.size")
             :=lapply(.SD, function(x) x*sqrt(var.body.size.direct)), .SDcols=c("direct.genetic.body.size")] # makes variance proper
  
  
  #make a subset of the males which will mate, this must be moved into the mating function 
  effgen0.males <- subset( gen0.males,  mating.willingness.1st > 0 ) 
  return (effgen0.males)
}
############### Mating list and mate function #################
mate <- function (x,y) { # x = males, y = females
  mating.list <- x[rep(seq(nrow(x)), mating.willingness.1st),  #expands the male list into a mating list, based on mat.will.1st
                   c("id","fert", "direct.genetic.body.size", "semen.quality.1st","semen.quality.2nd","can.remate","mating.willingness.1st","mating.willingness.2nd") 
                   , with=F] #specify which columns to incl.
  if (nrow(mating.list) > sum(y$mating.will.1st.round)) {
    mating.list <- mating.list[1:sum(y$mating.will.1st.round),]
  }
  
  setnames(mating.list, c("id", "direct.genetic.body.size", "fert"), c("sire.id.1st", "sire.bs.1st", "sire.fert.1st"))
  setkey(mating.list, sire.id.1st)
  mating.list[mating.list, c("dam.id","dam.fert", "f0","obs_fert","dam.bs","perm.env.ls","birthyear.dam","sire.id.2nd","sire.bs.2nd","sire.fert.2nd"):=0]
  mating.list[, `:=`(perm.env.bs=rnorm(nrow(mating.list))*sqrt(var.body.size.spec.env))]
  # here I subset the dam list to throw away those who will not mate on first round
  y <- subset (y, mating.will.1st.round ==1)
  mating.list$dam.id       <- y[1:nrow(mating.list),.(id)] 
  mating.list$dam.fert     <- y[1:nrow(mating.list),.(fert)] 
  mating.list$dam.bs       <- y[1:nrow(mating.list),.(direct.genetic.body.size)] 
  mating.list$perm.env.ls  <- y[1:nrow(mating.list),.(perm.env.ls)] 
  if("f0" %in% colnames(y)) { 
    mating.list$f0  <- y[1:nrow(mating.list),.(f0)] 
  }
  if("perm.env.bs" %in% colnames(y)) { 
    mating.list$perm.env.bs  <- y[1:nrow(mating.list),.(perm.env.bs)] 
  }
  mating.list$birthyear.dam  <- y[1:nrow(mating.list),.(birthyear)] 
  setnames(mating.list, "f0", "f0.dam") 
  if (year == 1 ) {
    setnames(mating.list, "f0.dam", "f0") # this is because of in the breeding values of the first gen this is needed, silly really since its zero all over
  }
  # now to make the 2nd round of matings
  mating.list.allowed <- subset(mating.list, can.remate == 1)
  mating.list.allowed[ , `:=`( IDX = 1:.N ) , by = sire.id.1st ]
  mating.list.allowed$sire.id.2nd       <- ifelse(mating.list.allowed$IDX <= mating.list.allowed$mating.willingness.2nd, mating.list.allowed$sire.id.1st   , 0 )
  mating.list.allowed$semen.quality.2nd <- ifelse(mating.list.allowed$sire.id.1st == mating.list.allowed$sire.id.2nd, mating.list.allowed$semen.quality.1st, 0 )
  mating.list.allowed$sire.fert.2nd     <- ifelse(mating.list.allowed$sire.id.1st == mating.list.allowed$sire.id.2nd, mating.list.allowed$sire.fert.1st    , 0 )
  mating.list.allowed$sire.bs.2nd       <- ifelse(mating.list.allowed$sire.id.1st == mating.list.allowed$sire.id.2nd, mating.list.allowed$sire.bs.1st      , 0 )
  
  #subset mating.list.allowed into two, those done and those left
  
  mating.list.remated <- subset (mating.list.allowed, sire.id.2nd != 0)
  mating.list.leftover <- subset (mating.list.allowed, sire.id.2nd == 0) # those females who were allowed to
  # be remated with 1st male but male did not have enough to mate them all
  mating.list.notallowed <- subset(mating.list, can.remate == 0)
  x[,`:=`(matings.left = mating.willingness.2nd - mating.willingness.1st)] # figures out how many matings the males performed that did have spares
  x$matings.left <- ifelse( x$matings.left < 0, 0, x$matings.left  )  # if they have negs, they are done
  x <- subset(x, matings.left > 0 & can.remate ==1)  # only the males that have spare matings left
  
  set(mating.list.leftover, j=c("IDX"), value = NULL)
  mating.list.leftover <- rbind(mating.list.leftover, mating.list.notallowed)
  
  # here I should reorder the dams that are leftovers to prioritize younger dams and the ones mated with shitty males
  mating.list.leftover <- as.matrix(mating.list.leftover) # way faster to loop through a matrix
  x <- as.matrix(x)
  if (year == 1) {
 for ( i in 1:nrow(x) )  {
    #print(i)
    for ( j in 1:x[[i,13]] )  { 
      s <-  sum( x[1:i,13] ) # keeps track of how many females the male has been assigned
      #         print(s)
      #     print(j)
      if (s < nrow(mating.list.leftover)) { # this is not a solution, need to make another controlling mechanism if the mating willingness
        #exceeds the number of females to be mated
        mating.list.leftover       [[s- ( j- 1 ),16]] <-  x[[i,1]]          # id of male
        mating.list.leftover       [[s- ( j- 1 ),5]]  <-   x[[i,7]]          # semen.quality
        mating.list.leftover       [[s- ( j- 1 ),18]] <-  x[[i,4]]          # fertility of male
        mating.list.leftover       [[s- ( j- 1 ),17]] <-  x[[i,5]]          # body size of male
        }
      if (s> nrow(mating.list.leftover)) { 
        while (sum(x[ 1:(i-1) ,13] ) + j <= nrow(mating.list.leftover) ) {
          # Here i solve the problem of if the males have more mating willingness than the number of females
          mating.list.leftover       [[ sum ( x[ 1:( i - 1 ),13 ] ) + j,16]] <-  x[[i,1]]          # id of male
          mating.list.leftover       [[ sum ( x[ 1:( i - 1 ),13 ] ) + j, 5]] <-  x[[i,7]]         # semen.quality
          mating.list.leftover       [[ sum ( x[ 1:( i - 1 ),13 ] ) + j, 18]] <-  x[[i,4]]         # fertility of male
          mating.list.leftover       [[ sum ( x[ 1:( i - 1 ),13 ] ) + j, 17]] <-  x[[i,5]]         # body size of male
          
          break
        }
      }
    }
 }}
  if (year > 1) {
    for ( i in 1:nrow(x) )  {
      #print(i)
      for ( j in 1:x[[i,18]] )  { 
        s <-  sum( x[1:i,18] ) # keeps track of how many females the male has been assigned
        #         print(s)
        #     print(j)
        if (s < nrow(mating.list.leftover)) { # this is not a solution, need to make another controlling mechanism if the mating willingness
          #exceeds the number of females to be mated
          mating.list.leftover       [[s- ( j- 1 ),16]] <-  x[[i,1]]          # id of male
          mating.list.leftover       [[s- ( j- 1 ),5]]  <-   x[[i,17]]          # semen.quality
          mating.list.leftover       [[s- ( j- 1 ),18]] <-  x[[i,10]]          # fertility of male
          mating.list.leftover       [[s- ( j- 1 ),17]] <-  x[[i,11]]          # body size of male
        }
        if (s> nrow(mating.list.leftover)) { 
          while (sum(x[[ 1:(i-1) ,18]] ) + j <= nrow(mating.list.leftover) ) {
            # Here i solve the problem of if the males have more mating willingness than the number of females
            mating.list.leftover       [[ sum ( x$matings.left [ 1:( i - 1 ) ] ) + j, 16]] <-  x[[i,1]]          # id of male
            mating.list.leftover       [[ sum ( x$matings.left [ 1:( i - 1 ) ] ) + j,  5]] <-  x[[i,17]]         # semen.quality
            mating.list.leftover       [[ sum ( x$matings.left [ 1:( i - 1 ) ] ) + j, 18]] <-  x[[i,10]]         # fertility of male
            mating.list.leftover       [[ sum ( x$matings.left [ 1:( i - 1 ) ] ) + j, 17]] <-  x[[i,11]]         # body size of male
            
            break
          }
        }
      }
    }
  }
  mating.list.leftover <- as.data.table(mating.list.leftover)
  x <- as.data.table(x)
  set(mating.list.remated, j= c("IDX"), value= NULL) # need to remove the counter to bind the mated ones to the rest
  
  mating.list <- rbind(mating.list.remated, mating.list.leftover) # merge the now completed mating list together
  mating.list[,`:=`(semen.quality= ifelse ( mating.list$semen.quality.2nd == 0, 0, mating.list$semen.quality.2nd), 
                    barren = # this rather cumbersome function is to have different probs of barren acc. to mating method
                      ifelse (year-mating.list$birthyear.dam == 1 & mating.list$sire.id.2nd != 0 , rbinom(nrow(mating.list),1,pr.barren.double.mating.yearling),
                              ifelse(year-mating.list$birthyear.dam != 1 & mating.list$sire.id.2nd != 0, rbinom(nrow(mating.list),1,pr.barren.double.mating.old), 
                                     ifelse( year - mating.list$birthyear.dam ==1 & mating.list$sire.id.2nd == 0, rbinom(nrow(mating.list),1,pr.barren.one.mating.yearling),
                                             ifelse(year-mating.list$birthyear.dam != 1 & mating.list$sire.id.2nd == 0, rbinom(nrow(mating.list),1,pr.barren.one.mating.old),0)))  
                      ))] 
  # remove columns who're not needed at this point
  set(mating.list, j= c("mating.willingness.1st", "mating.willingness.2nd", "can.remate","semen.quality.1st", "semen.quality.2nd"), value=NULL)
  
  return(mating.list)
} 
mate <-compiler::cmpfun(mate,options= c(suppressAll=TRUE)) # performance boost

############### make gen 0 pedigree & and some bookkeeping ##############
# This is not really needed at this point since all gen0 animals are unrelated
make.pedfile.gen0 <- function() {
setkey(gen0.females, id)
setkey(effgen0.males, id)
  pedfile <- gen0.females[,.(id,sire.id,dam.id,sex,birthyear)] #use .(colnames) instead of c(colnames) in data.table
  l <- list (pedfile, effgen0.males[,.(id,sire.id,dam.id,sex,birthyear)])
  pedfile <- rbindlist( l,use.names = TRUE ) 
  return (pedfile)
} 
############### Breeding value of offspring ###############################
# put in dam id's and make genetic value of each kit for fertility
# common cause for crash was a negative litter size
bv <- function () {
  gen1 <- mating.list[rep(seq(nrow(mating.list)), obs_fert), 
                      c("dam.id","sire.id.1st","dam.fert","sire.fert.1st","f0","obs_fert","dam.bs","sire.bs.1st"
                        ,"perm.env.bs","birthyear.dam","sire.id.2nd","sire.fert.2nd","sire.bs.2nd") 
                      , with=F] #specify which columns to incl.
  id <- seq(1:sum(mating.list$obs_fert)) + max(effgen0.males$id) # makes ID
  birthyear <- rep (1, sum(mating.list$obs_fert)) # makes birthyear
  sex <- rbinom(sum(mating.list$obs_fert),1,0.5)+1 # makes sex, TODO check if this is true
  true.sire <- numeric(nrow(gen1))
  true.sire.check <- numeric(nrow(gen1))
  mendelian <- as.data.table(rmvnorm(sum(mating.list$obs_fert),sigma=sigma))
  colnames(mendelian) <- c("mend.fert","mend.bw")
  
  gen1 <- cbind (id,gen1,  birthyear, sex, mendelian, true.sire,true.sire.check) # binds id, sex and birthyear to data.table
  
  gen1$true.sire.check <- ifelse( gen1$sire.id.1st != gen1$sire.id.2nd, rbinom(nrow(gen1), 1, 0.85), 1) # 85% chance that the kits are sired by 2nd mating
  gen1$true.sire <- ifelse( gen1$true.sire.check == 0, gen1$sire.id.1st, gen1$sire.id.2nd)
  gen1[, `:=`(true.sire.fert = ifelse(gen1$true.sire == gen1$sire.id.2nd, gen1$sire.fert.2nd, gen1$sire.fert.1st), 
              true.sire.bs = ifelse(gen1$true.sire == gen1$sire.id.2nd, gen1$sire.bs.2nd, gen1$sire.bs.1st))]
  gen1[, `:=`(fert = 0.5*(dam.fert + true.sire.fert) + 
                mend.fert*(sqrt(0.5*variance.fertility)*( 1- f0 )) # Breeding value of offspring, littersize
              , perm.env.ls = rnorm(sum(mating.list$obs_fert))*sqrt(var.perm.env.ls) # perm env for litter size
              , direct.genetic.body.size = 0.5*(dam.bs + true.sire.bs) + 
                mend.bw*(sqrt(0.5*var.body.size.direct)*( 1- f0 )))]# Breeding value of offspring, body size
  #  gen1 <- count.sex.siblings(gen1) # calls function to count offspring NOT NEEDED ATM
  
  setnames(gen1, c("obs_fert","sire.id.2nd"), c("own_littersize","sire.assumed")) # renames obs_fert to own littersize of kits
  gen1$dam.age <- ifelse( year - gen1$birthyear.dam > 1, 1,0 )
  
  gen1$bs.phenotype <- ifelse( gen1$sex == 1,phenotype.bs.male(mean.body.size.male , gen1$direct.genetic.body.size , gen1$perm.env.bs , gen1$own_littersize, gen1$dam.age  ) 
                               , phenotype.bs.female(mean.body.size.female , gen1$direct.genetic.body.size , gen1$perm.env.bs,gen1$own_littersize )) 
  # generate phenotype for body size
  
  gen1[,`:=`(perm.env.bs = rnorm(sum(mating.list$obs_fert))*sqrt(var.body.size.spec.env))] # generate specific env. for body size
  
  set( gen1, j=which(colnames(gen1) %in% c("dam.fert","sire.fert","sire.bs","dam.bs","mend.fert"
                                           ,"mend.bw","birthyear.dam","dam.age","sire.id.1st","sire.bs.1st",
                                           "sire.fert.1st","sire.fert.2nd","sire.bs.2nd"
                                           ,"true.sire.fert","true.sire.bs", "true.sire.check")) , value=NULL ) # removes bv of parents
  return(gen1)
}
############### Selection of old females ###############################################
# currently they are only truncated on their litter size phenotype
sel.old.females <- function (y) {
  setkey(mating.list, obs_fert)
  setorder(mating.list, -obs_fert)
  if ("birthyear.dam" %in% colnames(mating.list)) {
    mating.list <- subset(mating.list, year -birthyear.dam < max.age.females) 
  }
  old.females <- mating.list[1:(n.females*prop.oldfemales),]
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
  if("bs.phenotype" %in% colnames(old.females)){
    
  } else {
    old.females <-transform(old.females, bs.phenotype = mean.body.size.female + direct.genetic.body.size + rnorm(1)*sqrt(var.res.body.size))
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
  return(old.females)
  
  
}  
############### Selection of yearling females in 1st gen ###############################

sel.yearlings.females <- function (x) { # x = kit.list  
  truncation.point <-  quantile( x$own_littersize,  probs =  quantile.setting ) 
  selection.candidates.females <- subset(x, own_littersize >= truncation.point) # throw away the smallest litters
  selection.candidates.females <-  subset( selection.candidates.females,  sex  ==   2) # take the female kits
  selection.candidates.females <- subset( selection.candidates.females, bs.phenotype < roof.body.size)
  if (nrow(selection.candidates.females) < n.females*(1-prop.oldfemales)){ 
    truncation.point <-  quantile( x$own_littersize,  probs =  (quantile.setting  ) ) #can't change this since the kits won't have cards 
    selection.candidates.females <- subset(x, own_littersize >= truncation.point) # throw away the smallest litters
    selection.candidates.females <-  subset( selection.candidates.females,  sex  ==   2) # take the female kits
    selection.candidates.females <- subset( selection.candidates.females, bs.phenotype < (roof.body.size+200)) #ease restrictions on size
    
  }
  setkey(selection.candidates.females, bs.phenotype)           # in order to speed up ordering 
  setorder(selection.candidates.females,-bs.phenotype)         # order the kits according to body size 
  next.gen <- selection.candidates.females[1:(n.females-nrow(old.females)),]              # take the biggest kits
  setkey(next.gen, id)
  next.gen[next.gen,obs_fert:=0]  
  if("f0.dam" %in% colnames(next.gen)) {
    set( next.gen, j=which(colnames(next.gen) %in% 
                             "f0.dam")  , value=NULL )
    
  }
  
  return (next.gen)
}
############### Selection of yearling males in 1st gen ###############################

sel.males <- function (x) {
  truncation.point <-  quantile( x$own_littersize,  probs =  quantile.setting ) 
  selection.candidates.males <- subset(x, own_littersize >= truncation.point) # throw away the smallest litters
  selection.candidates.males <-  subset( selection.candidates.males,  sex  ==   1  ) # take the female kits
  setkey(selection.candidates.males, bs.phenotype)           # in order to speed up ordering 
  setorder(selection.candidates.males,-bs.phenotype)         # order the kits according to litter size 
  next.gen.males <- selection.candidates.males[1:n.males,]              # take the kits from biggest litters
  set( next.gen.males, j=which(colnames(next.gen.males) %in%
                                 "perm.env.ls"), value=NULL)
  if("f0.dam" %in% colnames(next.gen.males)) {
    set( next.gen.males, j=which(colnames(next.gen.males) %in% "f0.dam")  , value=NULL )
    
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
############### Make barren males ##################
make.barren.males <- function () { 
  setkey(next.gen.males, id)
  #   next.gen.males[next.gen.males,semen.quality:=rbinom( nrow(next.gen.males),  1,  male.inf )]
  #   next.gen.males[next.gen.males,mating.willingness:=rZIP( nrow(next.gen.males),  mu = male.ratio,  sigma = 0.05 )]
  next.gen.males[,`:=`( semen.quality.1st=  rbinom( nrow(next.gen.males),  1,  male.inf ),  
                        mating.willingness.1st = rZIP( nrow(next.gen.males),  mu = male.ratio,  sigma = 0.05 ),
                        mating.willingness.2nd = rZIP( nrow(next.gen.males), mu = male.ratio, sigma = 0.05),
                        can.remate = rbinom(nrow(next.gen.males), 1,0.8))]
  next.gen.males[,`:=`(semen.quality.2nd = semen.quality.1st)]
  next.gen.males <- subset( next.gen.males, mating.willingness.1st > 0 ) 
  return(next.gen.males)
} 
############### mating will of females #############
mat.will <- function () {
  next.gen[,`:=`(dam.age = ifelse(year-birthyear != 1, 0,1))] # 1 = older females, 0 = yearlings
  
  next.gen[,`:=`( # this makes the mating will of the females in the first round of matings
    mating.will.1st.round = ifelse( dam.age == 1, rbinom(nrow(next.gen), 1,mating.will.old.1st),
                                    ifelse(dam.age == 0, rbinom(nrow(next.gen),1, mating.will.yearling.1st),0))
  )]
  
  next.gen[,`:=`( # this makes the mating will of the females in the second round of matings
    mating.will.2nd.round = ifelse( dam.age == 1 & mating.will.1st.round == 1, rbinom(nrow(next.gen), 1,mating.will.old.2nd),
                                    ifelse(dam.age == 0 & mating.will.1st.round ==1, rbinom(nrow(next.gen),1, mating.will.yearling.2nd),0))
  )]
  return(next.gen)
}
############### Update pedigree ########################
update.pedigree <- function () {
if (year ==2) {
  pedfile[,`:=`(true.sire = sire.id)]
  setnames(pedfile, "sire.id", "sire.assumed")
}
  pedfile <- rbind(pedfile, # orginal pedfile
                   next.gen[, .( id, sire.assumed, true.sire, dam.id, sex, birthyear ) ], #current gen females 
                   next.gen.males[, .( id, sire.assumed, true.sire, dam.id, sex, birthyear )] ) # current gen males
  dups <- duplicated(pedfile$id)
  pedfile <- pedfile[!dups,]
  return (pedfile)
}
############### Breeding value of offspring in gen n###############################
# put in dam id's and make genetic value of each kit for fertility
# common cause for crash was a negative litter size
bv.n <- function () {
  
  
  kit.list <- mating.list[rep(seq(nrow(mating.list)), obs_fert), # makes the kit list by expanding the mating list
                          c("dam.id","sire.id.1st","dam.fert","sire.fert.1st","f0.dam","obs_fert","dam.bs","sire.bs.1st"
                            ,"perm.env.bs","birthyear.dam","sire.id.2nd","sire.fert.2nd","sire.bs.2nd") , with=F]
  
  id <- seq(1:sum(mating.list$obs_fert)) + max(pedfile$id) # makes ID
  birthyear <- rep (year, sum(mating.list$obs_fert)) # makes birthyear
  sex <- rbinom(sum(mating.list$obs_fert),1,0.5)+1 # makes sex, TODO check if this is true
  true.sire <- numeric(nrow(kit.list))
  true.sire.check <- numeric(nrow(kit.list))
  mendelian <- as.data.table(rmvnorm(nrow(kit.list),sigma=sigma))
  colnames(mendelian) <- c("mend.fert","mend.bw")
  
  kit.list <- cbind (id,kit.list,  birthyear, sex,mendelian, true.sire, true.sire.check)
  # now check the true sire before calculating inbreeding
  kit.list$sire.id.2nd  <- ifelse( kit.list$sire.id.2nd == 0, kit.list$sire.id.1st, kit.list$sire.id.2nd) # puts the 1st sire into the 2nd sire column for single mated females
  kit.list$sire.fert.2nd <- ifelse( kit.list$sire.id.2nd == kit.list$sire.id.1st, kit.list$sire.fert.1st, kit.list$sire.fert.2nd)
  kit.list$sire.bs.2nd <- ifelse( kit.list$sire.id.2nd == kit.list$sire.id.1st, kit.list$sire.bs.1st, kit.list$sire.bs.2nd)
  
    kit.list$true.sire.check <- ifelse( kit.list$sire.id.1st != kit.list$sire.id.2nd, rbinom(nrow(kit.list), 1, 0.85), 1) # 85% chance that the kits are sired by 2nd mating
  kit.list$true.sire <- ifelse( kit.list$true.sire.check == 0, kit.list$sire.id.1st, kit.list$sire.id.2nd)
  # kit.list$true.sire <- ifelse( kit.list$sire.id.1st != 0 & kit.list$sire.id.2nd == 0, kit.list$sire.id.1st, kit.list$true.sire) # if the females are only mated once, the true sire is the first
  

  
  kit.list[, `:=`(true.sire.fert = ifelse(kit.list$true.sire == kit.list$sire.id.2nd, kit.list$sire.fert.2nd, kit.list$sire.fert.1st), 
              true.sire.bs = ifelse(kit.list$true.sire == kit.list$sire.id.2nd, kit.list$sire.bs.2nd, kit.list$sire.bs.1st))]
  
  setnames(kit.list, "sire.id.2nd", "sire.assumed")
  pedfile1 <- rbind(pedfile, # orginal pedfile
                    kit.list[, .( id,sire.assumed, true.sire, dam.id, sex, birthyear ) ] ) #puts the kits into the pedfile
  
  ttt = trimPed(pedfile1,is.element( pedfile1$id, kit.list$id)) #trims the pedfile into only those related to kits
  t4 = cbind(pedfile1,ttt) #bind the logical vector to the temporary file
  if( make.obs.file == 1) {
    pedigree <- file(description = paste("pedigree", sep=""), open="w")
    
    write.table(pedfile1[,.(id,sire.assumed,dam.id,birthyear)], file= pedigree, col.names=FALSE, row.names=FALSE, quote=FALSE)
    close(pedigree)}
  if( use.true.sire == 1){
    pedigree.true <- file(description = paste("pedigree.true", sep=""), open="w")
    
    write.table(pedfile1[,.(id,true.sire,dam.id,birthyear)], file= pedigree.true, col.names=FALSE, row.names=FALSE, quote=FALSE)
    close(pedigree.true)}
  
   
  
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
                  , perm.env.ls = rnorm(sum(mating.list$obs_fert))*sqrt(var.perm.env.ls) # perm env for litter size
                  , direct.genetic.body.size = 0.5*(dam.bs + true.sire.bs) + 
                    mend.bw*(sqrt(0.5*var.body.size.direct)*( 1- f0 )))]# Breeding value of offspring, body size
  
  setnames(kit.list, "obs_fert", "own_littersize") # changes obs fert into the littersize of the kit
  kit.list$dam.age <- ifelse( year - kit.list$birthyear.dam > 1, 1,0 )
  kit.list$bs.phenotype <- ifelse( kit.list$sex == 1,phenotype.bs.male(mean.body.size.male , kit.list$direct.genetic.body.size , kit.list$perm.env.bs , kit.list$own_littersize, kit.list$dam.age  ) 
                                   , phenotype.bs.female(mean.body.size.female , kit.list$direct.genetic.body.size , kit.list$perm.env.bs,kit.list$own_littersize )) 
  # generate phenotype for body size
  
  kit.list[,`:=`(perm.env.bs = rnorm(sum(mating.list$obs_fert))*sqrt(var.body.size.spec.env))] # generate specific env. for body size
  
  set( kit.list, j=which(colnames(kit.list) %in% c("dam.fert","sire.fert","sire.bs","dam.bs","mend.fert"
                                                   ,"mend.bw","birthyear.dam","dam.age","sire.id.1st","sire.bs.1st",
                                                   "sire.fert.1st","sire.fert.2nd","sire.bs.2nd"
                                                   ,"true.sire.fert","true.sire.bs", "true.sire.check")) , value=NULL ) # removes bv of parents
  return(kit.list)
} 
bv.n <-compiler::cmpfun(bv.n,options= c(suppressAll=TRUE)) # performance boost

############### write observation file #########################
write.output <- function () {

# this is to format the mating list to write to the observation file for this replicate for the calculation of BV
  # fertility first and then write the observations for body weight of kits
# set(mating.list, j = which(colnames(mating.list) %in% 
#                              c("sire.id", "f0.dam", "dam.fert", "perm.env.ls", "sire.fert", "barren")), value= NULL)
  mating.list[,`:=`(dam.age=0)]
  
  mating.list[, c("dam.id","birthyear.dam")
            :=lapply(.SD, function(x) as.integer(x)), .SDcols=c("dam.id","birthyear.dam")]
  mating.list[,`:=`(prodyear=as.integer(year), 
                    dam.age = ifelse( year-mating.list$birthyear.dam == 1, as.integer(1), as.integer(2)))]
  
# for (i in 1:nrow(mating.list)) { # this is to change the parity into a two class thing
#   if (year - mating.list$birthyear.dam[i] == 1) {
#     mating.list$dam.age [i] <- as.integer(1)
#   }
#   else {
#     mating.list$dam.age[i] <- as.integer(2)
#   }
# }
# set(mating.list, j = which(colnames(mating.list) %in%
#                              c("birthyear.dam")), value= NULL)
setkey(mating.list, dam.id)
mating.list[mating.list,prodyear:=as.integer(year)]
mating.list[, c("dam.id","dam.age")
            :=lapply(.SD, function(x) as.integer(x)), .SDcols=c("dam.id","dam.age")]

# if("f0" %in% colnames(mating.list)) {
#   set( mating.list, j=which(colnames(mating.list) %in%
#                               "f0")  , value=NULL )
# }
# setcolorder(mating.list, c("dam.id", "prodyear" , "dam.age" , "obs_fert"))
# x <- mating.list[,.(dam.id,prodyear,dam.age,obs_fert)]
write.table(format(mating.list[,.(dam.id,prodyear,dam.age,obs_fert)], nsmall=1), file= output, append= TRUE,col.names = FALSE, row.names = FALSE, quote = FALSE)

} 
############### make dam.age checks #################
dam.age <- function (){
  dam.age <- numeric(nrow(mating.list))
  for (i in 1:nrow(mating.list)) { # this is to change the parity into a two class thing
    if (year - mating.list$birthyear.dam[i] == 1) {
      dam.age [i] <- yearling.effect
    }
    else {
      dam.age[i] <- 0
    }
  }
  mating.list <- cbind(mating.list, dam.age)
  return(mating.list)
}

############  Write pedigree file to DMU
 write.pedigree <- function ( ) {
   if (make.obs.file == 1) {
     if (use.true.sire== 1) {
       pedigree.true <- file(description = paste("pedigree.true", sep=""), open="w")
       
       write.table(pedfile[,.(id,true.sire,dam.id,birthyear)], file= pedigree.true, col.names=FALSE, row.names=FALSE, quote=FALSE)
       close(pedigree.true)
     }
   # set(pedfile, j = which(colnames(pedfile) %in%
   #                          "sex","true.sire"), value=NULL )
     write.table(pedfile[,.(id,sire.assumed,dam.id,birthyear)], file= pedigree, col.names=FALSE, row.names=FALSE, quote=FALSE)
   close(con=pedigree)
     }
 }
############### Make phenotype for body size ######
 phenotype.bs.male <- function(x,y,z,t,u) { # x = mean, y = additive genetic, z = specific gen, t = number of  sibs
   value <- x+ y + z + rnorm( sum(mating.list$obs_fert))*(sqrt(var.res.body.size))+t*sib.effect.male+ u * bw.eff.damage
   return(value)
 } 
 phenotype.bs.female <- function(x,y,z,t) { # x = mean, y = additive genetic, z = specific gen, t = number of sibs
   value <- x+ y + z + rnorm( sum(mating.list$obs_fert))*(sqrt(var.res.body.size))+t*sib.effect.female 
   return(value)
 } 
 
 
 
 ############ Count sex of siblings ############
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
 