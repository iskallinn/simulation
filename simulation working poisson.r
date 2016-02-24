source("C:/Users/Notandi/Dropbox/Projects/simulation/source_code.r")
# Rprof("profile2.out", line.profiling=TRUE)
# solutions <-numeric(n*nruns)
setwd("C:/Users/Notandi/Dropbox/Projects/simulation of mink farm/Output/DMU analysis/")
con <- file(description="results",open="w")
cat("Gen","Gmean","Gvar","Fis","Obs.fert", "mean.phenotype.bs.females", "gen.value.bs", "mean.phenotype.bs.males","bw.var",sep="\t",file=con)
cat("\n",file=con)

n.males <-  ceiling( n.females/male.ratio ) # calculates needed amount of males 


for (t in 1:nruns) {

  if (make.obs.file == 1) {
  output <- file(description = paste("Replicate", sep=""), open="w")
#cat("ID","prodyear","damage","obs.fert", sep="\t", file=output) 
#cat("\n", file= output) 
  }
  
  if (make.obs.file == 1) {
    pedigree <- file(description = paste("pedigree", sep=""), open="w")
  }
year <- 1
# runcounter <- sum(t)


############### Create base population ############
gen0.females <- generate.base.females() # create females
effgen0.males <- generate.base.males() # create males

############### Make mating list ##########################
################## assign each female a male,  based on his mating willingness ####
mating.list <- mate (effgen0.males,gen0.females)
############### make gen 0 pedigree & and some bookkeeping ##############
# This is not really needed at this point since all gen0 animals are unrelated

# mating.list <- data.table( mating.list ) 
# mating.list <- subset( mating.list,  sire.id>0 ) # remove the females who were not mated
mating.list <- dam.age ()  

mating.list = transform( mating.list,  obs_fert =  rpois(nrow(mating.list), 
                                                         lambda = exp(1.95 + perm.env.ls + dam.fert+dam.age))*barren*semen.quality)
# TODO how to implement inbreeding depression in this poisson model??

set( mating.list, j=which(colnames(mating.list) %in% c("semen.quality","dam.age")) , value=NULL )

mating.list <-  subset( mating.list,  obs_fert >  0 ) # remove females who are barren or mated with barren male
# it is important to get rid of negative litter sizes! otherwise code crashes at many points

############### Breeding value of offspring ###############################
# put in dam id's and make genetic value of each kit for fertility
# see utility functions for bv function
gen1 <- bv()

pedfile <- make.pedfile.gen0()

# At this point I think it is safe to delete some stuff from memory
# Note that I skip the construction of pedfile1 here. I don't think it is needed. Will check on that later
############### Selection of first generation #########################################

# We choose females to fit n.females, based on prop of old females
# Note: currently all females are replaced in year 1

if (selection.method == selection) {
  old.females <- sel.old.females ( gen0.females)
  next.gen <- sel.yearlings.females (gen1)
  next.gen.males <- sel.males (gen1)
  if("f0.dam" %in% colnames(old.females)) {
    set( old.females, j=which(colnames(old.females) %in% 
                                "f0.dam")  , value=NULL )
    }
  next.gen <- rbind(next.gen, old.females)
} else if (selection.method == no.selection) {
  old.females <- nosel.old.females(gen0.females) #old females
  next.gen <- nosel.yearlings.females(gen1) # females
  next.gen.males <- nosel.males(gen1)  # males
  next.gen <- rbind(next.gen, old.females)
}
############## First year statistics #######################
cat (0, mean(mating.list$dam.fert),var(mating.list$dam.fert),0,mean(mating.list$obs_fert), mean(next.gen$bs.phenotype)
     , mean(next.gen$direct.genetic.body.size),mean(next.gen.males$bs.phenotype), var(next.gen$direct.genetic.body.size), sep="\t",file=con)
cat("\n",file=con)

# solutions[[year+(runcounter-1)]] <- mean(next.gen$bs.phenotype)

if (make.obs.file == 1) {
  write.output()
}

remove(mating.list,gen0.females,gen1,effgen0.males,old.females) #remove all the stuff I don't need anymore

# Start of big loop 
for (y in 1:n) {
  print( y )
  temp <- copy(next.gen)
  set( temp, j=which(colnames(temp) %in% c("obs_fert","perm.env"))  , value=NULL ) 
  t <- rbind(next.gen.males,temp, fill=T) # needed down the road 
  remove(temp) 
  year <- 1+y  
############### Make barren males ##########################
  next.gen.males <- make.barren.males()
  next.gen <- mat.will()
############### assign each female a male,  based on his mating willingness #####
  mating.list <- mate(next.gen.males,next.gen)    
############### Update pedigree ########################
  pedfile <- update.pedigree()
############### Inbreeding Coefficient and calculation of breeding value#############################
#   mating.list <- subset( mating.list,  sire.id>0 ) # remove the females who were not mated
mating.list <- dam.age ()  # checks the dam age and puts the effect for yearlings
  
  mating.list = transform( mating.list,  obs_fert =  rpois(nrow(mating.list), 
                                                           lambda = exp(1.95 + perm.env.ls + dam.fert+dam.age)) 
                           *barren*semen.quality ) 
  set( mating.list, j=which(colnames(mating.list) %in% c("semen.quality","dam.age")) , value=NULL )
  
   kit.list <- bv.n()
############### Selection of next generation    #############
   # See utility functions for method
   if (selection.method == selection) {
     old.females <- sel.old.females ( next.gen)
     next.gen <- sel.yearlings.females (kit.list)
     next.gen.males <- sel.males (kit.list)
############### No selection of next generation #############
   } else if (selection.method == no.selection ) {
     old.females <- nosel.old.females ()
     next.gen <- nosel.yearlings.females (kit.file_df)
     next.gen.males <- nosel.males (kit.file_df)
   }
   if("f0.dam" %in% colnames(old.females)) {
     set( old.females, j=which(colnames(old.females) %in% 
                                 "f0.dam")  , value=NULL )}
   if("f0.dam" %in% colnames(next.gen.males)) {
     set( next.gen.males, j=which(colnames(next.gen.males) %in% "f0.dam")  , value=NULL )
   }
   
   next.gen <- rbind(next.gen,old.females)
      # solutions[year+(runcounter-1)] <- c(mean(next.gen$bs.phenotype))

    cat (y, mean(next.gen$fert),var(next.gen$fert),mean(mating.list$f0.dam)
         ,mean(mating.list$obs_fert), mean(next.gen$bs.phenotype), mean(next.gen$direct.genetic.body.size)
         ,mean(next.gen.males$bs.phenotype),var(next.gen$direct.genetic.body.size), sep="\t",file=con)
   cat("\n",file=con)
   
if (make.obs.file == 1) {
  write.output()
  }
   
   } 
write.pedigree()


}
close(con=con)
if (make.obs.file == 1) {
close(con=output) 
  }
closeAllConnections()
# Rprof(NULL)
