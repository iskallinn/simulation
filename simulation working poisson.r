############### Packages to use     ###############
library(gamlss.dist) # for making zero inflated poisson random deviates for mating willingness
library(pedigree)    # for calculating inbreeding and trimming pedigree
library(data.table)  # faster then data frames
library(mvtnorm)     # generating random deviates from a multivariate normal distribution
############### Controls for simulation ############
n.females <-  1000             # NUMBER OF FEMALES
nruns <- 2                    # how many replicates of the simulation
n <- 5                        # number of generation per replicate
mating.method <- assortative   # mating method, random or assortative
selection.method <- selection  # selection mating, selection = truncation,  no selection = next gen is chosen at random
make.obs.file <- 0 # 1 = make observation file, 0 otherwise
use.true.sire <- 0 # 1 if true sire of kits is wanted for BV prediction, 0 otherwise
############### Innter settings, change at own risk##########
male.ratio <-  6               # MALE TO FEMALE RATIO
male.inf <-  0.98              # % ODDS OF MALE BEING NOT BARREN
female.inf <- 0.9              # % ODDS OF FEMALE BEING NOT BARREN
prop.oldfemales <-  0.4        # Proportion of older females
ibd.fertility <- 0.06          # Inbreeding depression on fertility NOTE: not in use atm
max.age.females <- 3           # define how old the females can be
yearling.effect <- -0.07       # setting for yearling effect, current best guess is -0.07
sib.effect.male <- -13.4      # Effect on body size of (male) one extra kit in litter, Hansen(1997)
sib.effect.female <- -18.6    # Effect on body size of (female) one extra kit in litter, Hansen(1997)
quantile.setting <- 0.4        # amount of kits to throw away because of too low littersize
mating.will.yearling.1st          <- 0.95 # probability of yearling being mated
mating.will.yearling.2nd          <- 0.98 # probability of yearling being remated
mating.will.old.1st               <- 0.98 # probability of old female being mated
mating.will.old.2nd               <- 0.98 # probability of old female being remated
pr.barren.one.mating.yearling     <- 0.8 # probability of single mated yearling being barren
pr.barren.double.mating.yearling  <- 0.9 # probability of double mated yearling being barren
pr.barren.one.mating.old          <- 0.9 # probability of single mated old female being barren
pr.barren.double.mating.old       <- 0.95 # probaility of double mated old female being barren

############# Variance settings for traits ###################
variance.fertility     <-  0.0122738     # Genetic variance of fertility, live born
var.perm.env.ls        <-  0.0004464     # Variance of permanent environment of litter size of dam
var.body.size.direct   <-  40000         # Genetic variance of body size, direct effect 
var.body.size.spec.env <-  5300          # variance of specific environment (litter) on body size
var.res.body.size      <-  32500         # residual variance of body size
roof.body.size         <-  2400          # Maximum body weight of females to be allowed for selection
bw.eff.damage          <-  34.5          # effect of dam age on weight of males
############# (Co)Variance matrix ############################
sigma <-  matrix( 
  c(1, -0.5, -0.5, 1), # genetic correlations of traits
  nrow=2, 
  ncol=2) 
colnames(sigma) <- c("body.size.direct","litter.size")
rownames(sigma)<- c("body.size.direct","litter.size")

############# Mean settings for traits   #####################
mean.body.size.male <- 3000
mean.body.size.female <- 1650

##############Switches, best left alone ###########
assortative <- 1
random <- 0
no.selection <- 0
selection <- 1
blup <- 2
############## Connection for output ##############

# Rprof("profile2.out", line.profiling=TRUE)
stat.crate <-matrix(0, nrow = (n+1)*nruns, ncol = 3)
setwd("C:/Users/Notandi/Dropbox/Projects/simulation of mink farm/Output/DMU analysis/") 
con <- file(description="results",open="w")
cat("Gen","Gmean","Gvar","Fis","Obs.fert", "mean.phenotype.bs.females", "gen.value.bs", "mean.phenotype.bs.males","bw.var",sep="\t",file=con)
cat("\n",file=con)

n.males <-  ceiling( n.females/male.ratio ) # calculates needed amount of males 


for (p in 1:nruns) {

  if (make.obs.file == 1) {
  output <- file(description = paste("Replicate_",p, sep=""), open="w")
  } 
#cat("ID","prodyear","damage","obs.fert", sep="\t", file=output) 
#cat("\n", file= output) 
  
  
  # if (make.obs.file == 1) {
  #   pedigree <- file(description = paste("pedigree_",t, sep=""), open="w")
  # }
  
year <- 1
runcounter <- sum(p)
modify.dir.file ()

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

if (runcounter == 1 ) {
  stat.crate[1,2] <- nrow(mating.list)
} else if (runcounter > 1) { 
  stat.crate[year+(runcounter -1)*(n+1),2] <- nrow(mating.list)}
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


  old.females <- sel.old.females ( gen0.females)
  next.gen <- sel.yearlings.females (gen1)
  next.gen.males <- sel.males (gen1)
  if("f0.dam" %in% colnames(old.females)) {
    set( old.females, j=which(colnames(old.females) %in% 
                                "f0.dam")  , value=NULL )
  }
  next.gen <- rbind(next.gen, old.females)

############## First year statistics #######################
cat (0, mean(mating.list$dam.fert),var(mating.list$dam.fert),0,mean(mating.list$obs_fert), mean(next.gen$bs.phenotype)
     , mean(next.gen$direct.genetic.body.size),mean(next.gen.males$bs.phenotype), var(next.gen$direct.genetic.body.size), sep="\t",file=con)
cat("\n",file=con)

if (runcounter == 1) {
  stat.crate[[1,1]] <- (mean(gen1$true.sire == gen1$sire.assumed))
  stat.crate[1,3] <- nrow(gen1)
} else if (runcounter > 1) {
  stat.crate[[year+(runcounter -1)*(n+1),1]] <- (mean(gen1$true.sire == gen1$sire.assumed))
stat.crate[year+(runcounter -1)*(n+1),3] <-nrow(gen1)
  }
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
  
  stat.crate[year+(runcounter -1)*(n+1),2] <- nrow(mating.list)
  
    mating.list <- subset(mating.list, obs_fert > 0)
  if (make.obs.file == 1) {
    write.output()
  }
  
   kit.list <- bv.n()
   if (selection.method== blup) {
     solutions <- calculate.selection.index()
   }
############### Selection of next generation    #############
   # See utility functions for method
   if (selection.method == selection) {
     old.females <- sel.old.females ( next.gen)
     next.gen <- sel.yearlings.females (kit.list)
     next.gen.males <- sel.males (kit.list)
############### No selection of next generation #############
   } else if (selection.method == blup ) {
     old.females <- selind.old.females ()
     next.gen <- indexsel.yearlings.females (kit.list)
     next.gen.males <- indsel.males (kit.list)
   } 
   if("f0.dam" %in% colnames(old.females)) {
     set( old.females, j=which(colnames(old.females) %in% 
                                 "f0.dam")  , value=NULL )}
   if("f0.dam" %in% colnames(next.gen.males)) {
     set( next.gen.males, j=which(colnames(next.gen.males) %in% "f0.dam")  , value=NULL )
   }
   
   next.gen <- rbind(next.gen,old.females)
# gather up mean number of true sires
   stat.crate[year+(runcounter -1)*(n+1),1] <- c(mean(kit.list$true.sire == kit.list$sire.assumed))
   stat.crate[year+(runcounter -1)*(n+1),3] <- nrow(kit.list)
    cat (y, mean(next.gen$fert),var(next.gen$fert),mean(mating.list$f0.dam)
         ,mean(mating.list$obs_fert), mean(next.gen$bs.phenotype), mean(next.gen$direct.genetic.body.size)
         ,mean(next.gen.males$bs.phenotype),var(next.gen$direct.genetic.body.size), sep="\t",file=con)
   cat("\n",file=con)
   
   
   } 
# write.pedigree()
# rm(mating.list)

} 
# close(con=con)
if (make.obs.file == 1) {
  closeAllConnections()
# close(con=output) 
  }
closeAllConnections()
# Rprof(NULL)
