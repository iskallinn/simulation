WriteObservations <- function (mating.list, next.gen, next.gen.males,kit.list,year,p) {
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

CalculateBLUP <- function () {
  system2("C:/Users/Notandi/Dropbox/Projects/simulation of mink farm/Output/DMU analysis/run_dmu4.bat", " MBLUP")
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

IndSelectionOldFemales <- function (x,y,year) { # x = next.gen, y = solutions, z = solutions.bw
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
IndSelFemaleKits <- function (x,y,q) { # x = kit.list , y = solutions.fert, z = solutions.bw, q = old.females 
  x <- merge(x, y, by= "id", all.x=TRUE) # merge to solutions of blup of fertility
  x[, `:=`(index.bw   = 100+ (blup.bwnov-mean(x$blup.bwnov))/(sqrt(var(x$blup.bwnov)))*10,
           index.fert = 100+ (blup.fert-mean(x$blup.fert))/(sqrt(var(x$blup.fert)))*10,
           index.qual = 100+ (blup.qual-mean(x$blup.qual))/(sqrt(var(x$blup.qual)))*10) ]
  x <- transform(x, comb.ind = index.bw*weight.bw.kits+
                   index.fert*weight.fert.kits+
                   weight.qual.kits*index.qual)
  
  truncation.point <-  quantile( x$blup.fert,  probs =  quantile.setting ) 
  selection.candidates.females <- subset(x, blup.fert >= truncation.point) # throw away the smallest litters
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
IndSelMaleKits <- function (x,y) { # x = kit.list, y = solutions, z = solutions.bw
  x <- merge(x, y, by= "id", all.x=TRUE) # merge to solutions 
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
  
  return (next.gen.males)
  
}

ModifyDIRFile <- function (p) {
  # fertility modification
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
  
}

ModifyPriors <- function (p) {
  # here I should update the priors at some specified point in time
    # bw modification
    dirfile <- readLines("MBLUP.DIR")
    dirfile[8] <- c(paste("$DATA  ASCII (5,2,-9999) Phenotypes",p, sep="")) 
    # change the input file for BLUP so it uses the next outputfile
    writeLines(dirfile, "MBLUP.DIR")
 
  
  
}
