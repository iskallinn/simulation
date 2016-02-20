############### Index selection of old females ###############################################
# currently they are only truncated on their litter size phenotype
selind.old.females <- function () {
  
  if("blup.fert" %in% colnames(next.gen)) {
    set( next.gen, j=which(colnames(next.gen) %in% 
                                c("blup.fert", "sem.blup.fert"))  , value=NULL )
  }
  
  next.gen <- merge(next.gen, solutions, by ="id", all.x=TRUE)
  setkey(next.gen, blup.fert)
  setorder(next.gen, -blup.fert)
    mating.list <- subset(mating.list, year -birthyear.dam < max.age.females) 
  #then I will need to make rangering
    old.females <- next.gen[1:(n.females*prop.oldfemales),]
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


indexsel.yearlings.females <- function (x) { # x = kit.list  
  x <- merge(x, solutions, by= "id", all.x=TRUE) # merge to solutions of blup of fertility
  truncation.point <-  quantile( x$blup.fert,  probs =  quantile.setting ) 
  selection.candidates.females <- subset(x, blup.fert >= truncation.point) # throw away the smallest litters
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

indsel.males <- function (x) {
  x <- merge(x, solutions, by= "id", all.x=TRUE) # merge to solutions of blup of fertility
  truncation.point <-  quantile( x$blup.fert,  probs =  quantile.setting ) 
  selection.candidates.males <- subset(x, blup.fert >= truncation.point) # throw away the smallest litters
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


# script for running BLUP
calculate.selection.index <- function () {
  system2("C:/Users/Notandi/Dropbox/Projects/simulation of mink farm/Output/DMU analysis/run_dmu4.bat", " bl_ass")
  # read the solutions and only keep the predictions of BV (they're not that right)
  solutions <- as.matrix(read.table(file="bl_ass.SOL"))
  solutions <- solutions[-1:-sum(year+3),] # throw away the first lines
  solutions <- as.data.table(solutions)
  solutions <- subset(solutions, V1 == 4 ) # throw away the estimation of permanent environment
  set (solutions, j=c("V1","V2","V3", "V4", "V6", "V7"), value= NULL)
  setnames(solutions, c("V5","V8", "V9"),c("id", "blup.fert", "sem.blup.fert"))
  
  
  # kit.list1 <- merge(kit.list, solutions, by= "id", all.x=TRUE)
  # next.gen <- merge(next.gen, solutions, by ="id", all.x=TRUE)
  
  if (use.true.sire == 1) {
    system2("C:/Users/Notandi/Dropbox/Projects/simulation of mink farm/Output/DMU analysis/run_dmu4.bat", " bl_true")
    # read the solutions and only keep the predictions of BV (they're not that right)
    solutions <- as.matrix(read.table(file="bl_true.SOL"))
    solutions <- solutions[-1:-sum(year+3),] # throw away the first lines
    solutions <- as.data.table(solutions)
    solutions <- subset(solutions, V1 == 4 ) # throw away the estimation of permanent environment
    set (solutions, j=c("V1","V2","V3", "V4", "V6", "V7"), value= NULL)
    setnames(solutions, c("V5","V8", "V9"),c("id", "blup.fert", "sem.blup.fert"))
    
    kit.list1 <- merge(kit.list1, solutions, by= "id", all.x=TRUE)
  }
  return(solutions)
}