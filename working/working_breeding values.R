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
  
  kit.list[, `:=`(litter.size = 0.5*(dam.fert + true.sire.fert) + 
                mend.litter.size*(sqrt(0.5*variance.fertility)*(1-f0)) # Breeding value of offspring, littersize
              , perm.env.ls = rnorm(sum(x$obs_fert))*sqrt(var.perm.env.ls) # perm env for litter size
              , bw.oct = 0.5*(dam.bw.oct + true.sire.bw.oct) + 
                mend.bw.oct*(1-f0),
              bw.sept = 0.5*(dam.bw.sept + true.sire.bw.sept) + 
                mend.bw.oct*(1-f0),
              live.qual = 0.5*(dam.live.qual + true.sire.live.qual)+ sqrt(0.5)*(mend.live.qual)*(1-f0),
              skin.qual = 0.5*(dam.skin.qual + true.sire.skin.qual)+ sqrt(0.5)*(mend.skin.qual)*(1-f0),
              skin.length = 0.5*(dam.skin.length + true.sire.skin.length)+ sqrt(0.5)*(mend.skin.length)*(1-f0))]# Breeding value of offspring, body size
  
  setnames(kit.list, "obs_fert", "own_littersize") # changes obs fert into the littersize of the kit
  kit.list$dam.age <- ifelse( year - kit.list$birthyear.dam > 1, 1,0 )
  # generate phenotype for body size 
  kit.list$phenotype.bw.oct <-
    ifelse(
      kit.list$sex == 1,
      MakePhenotypesBWMalesOct(
        mean.body.size.male.oct ,
        kit.list$bw.oct ,
        kit.list$perm.env.bw.oct ,
        kit.list$own_littersize,
        kit.list$dam.age,
        x
      )
      ,
      MakePhenotypesBWFemalesOct(
        mean.body.size.female.oct ,
        kit.list$bw.oct ,
        kit.list$perm.env.bw.oct,
        kit.list$own_littersize,
        x
      )
    )
  # generate phenotype for body size
  kit.list$phenotype.bw.sept <- ifelse(
    kit.list$sex == 1,
    MakePhenotypesBWMalesSept(
      mean.body.size.male.sept ,
      kit.list$bw.sept ,
      kit.list$perm.env.bw.sept ,
      kit.list$own_littersize,
      kit.list$dam.age,
      x
    )
    ,
    MakePhenotypesBWFemalesOct(
      mean.body.size.female.sept ,
      kit.list$bw.sept ,
      kit.list$perm.env.bw.sept,
      kit.list$own_littersize,
      x
    )
  )
  
  
  kit.list[,`:=`(perm.env.bw.sept = ifelse(sex==1, rnorm(sum(x$obs_fert))*sqrt(var.c.bw.sept.male),
                                       rnorm(sum(x$obs_fert))*sqrt(var.c.bw.sept.female)),
             perm.env.bw.oct = ifelse(sex==1, rnorm(sum(x$obs_fert))*sqrt(var.c.bw.oct.male),
                                      rnorm(sum(x$obs_fert))*sqrt(var.c.bw.oct.female)))] # generate specific env. for body size
  
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
  )) , value=NULL )
  if (qual.classes == 5) {
    truncs <- qnorm(
      p = c(0.05, 0.3, 0.7, 0.95),
      mean = mean(kit.list$live.qual,
                  sd = sqrt(var(
                    kit.list$live.qual
                  ))),
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
        mean = mean(kit.list$live.qual,
                    sd = sqrt(var(
                      kit.list$live.qual
                    ))),
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
