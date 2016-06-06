############### Packages to use     ###############
library(gamlss.dist) # for making zero inflated poisson random deviates for mating willingness
library(pedigree)    # for calculating inbreeding and trimming pedigree
library(data.table)  # faster then data frames
library(mvtnorm)     # generating random deviates from a multivariate normal distribution
library(doBy)        # for gathering up statistics along the way
library(orthopolynom) # for the random regression model for BW and RFI

############## Index weights ############
# weights must sum to 1
weight.fert.old.females <- 1
weight.bw.old.females   <- 0
weight.qual.old.females <- 0
weight.fert.kits        <- 1
weight.bw.kits          <- 0
weight.qual.kits        <- 0
##############Switches, best left alone ###########
assortative <- 1
random <- 0
no.selection <- 0
phenotypic <- 1
blup <- 2
trace.ped <- 0
mask.phenotypes <- 1 # 0 no masking, 1 mask phenotypes of kits from small litters
sept <- 1 # weigh kits in september
oct <- 0  # weigh kits in october
use.blup.to.assort.mat <- 0 # if 1 then the males and females will be sorted on 
#their combined index before mating
use.comb.ind.for.males <- 1 # if 1 then the usage of the males will depend on 
# their combined index, not quality
mblup <- 0
############### Controls for simulation ############
n.females <-  1000             # NUMBER OF FEMALES
nruns <- 1                    # how many replicates of the simulation
n <- 5                        # number of generation per replicate
mating.method <- assortative   # mating method, random or assortative
selection.method <- blup  # selection mating, 
weighing.method <- oct         # control for when to weigh kits for selection cands
qual.classes <- 5              # quality classes, 5 or 10 are supported
intensity.remating <- 0.2      # this controls how many males are chosen not to be remated
# selection = truncation,  no selection = next gen is chosen at random
make.obs.file <- 1 # 1 = make observation file, 0 otherwise
use.true.sire <- 0 # 1 if true sire of kits is wanted for BV prediction, 0 otherwise
feed.price    <- 2.85 # feed price per kg feed
use.MBLUP     <-  1
############# Mean settings for traits   #####################
mean.body.size.male.oct <- 3750
mean.body.size.female.oct <- 1960
mean.body.size.male.sept <- 2270
mean.body.size.female.sept <- 1280
############### Innter settings, change at own risk##########
male.ratio <-  6               # MALE TO FEMALE RATIO
male.inf <-  0.98              # % ODDS OF MALE BEING NOT BARREN
female.inf <- 0.9              # % ODDS OF FEMALE BEING NOT BARREN
prop.oldfemales <-  0.4        # Proportion of older females
ibd.fertility <- 0.06          # Inbreeding depression on fertility NOTE: not in use atm
max.age.females <- 2           # define how old the females can be
yearling.effect <- -0.07       # setting for yearling effect, current best guess is -0.07
sib.effect.male <- -13.4       # Effect on body size of (male) one extra kit in litter, Hansen(1997)
sib.effect.female <- -18.6      # Effect on body size of (female) one extra kit in litter, Hansen(1997)
sib.effect.male.sept   <- -21.5 # Effect on body size of (female) one extra kit in litter, Hansen(1997)
sib.effect.female.sept <- -12.6 # Effect on body size of (female) one extra kit in litter, Hansen(1997)
sib.effect.male.june   <- -19.9       
sib.effect.female.june <- -15.9
eff.dam.age.june       <- 29
quantile.setting.ls <- 0.6        # amount of kits to throw away because of too low littersize
quantile.setting.bw <- 0.1    # prop of kits to disqualify due to weight
mating.will.yearling.1st          <- 0.95 # probability of yearling being mated
mating.will.yearling.2nd          <- 0.92 # probability of yearling being remated
mating.will.old.1st               <- 0.98 # probability of old female being mated
mating.will.old.2nd               <- 0.93 # probability of old female being remated
pr.barren.one.mating.yearling     <- 0.8 # probability of single mated yearling being barren
pr.barren.double.mating.yearling  <- 0.9 # probability of double mated yearling being barren
pr.barren.one.mating.old          <- 0.9 # probability of single mated old female being barren
pr.barren.double.mating.old       <- 0.95 # probability of double mated old female being barren
n.males <-  ceiling( n.females/male.ratio ) # calculates needed amount of males 
cull.ratio                       <- 0.8 # survival rate of kits, farmwise from 2nd cnt to pelting

############# Variance settings for traits ###################
# fertility
variance.fertility     <-  0.0122738     # Genetic variance of fertility, live born
var.perm.env.ls        <-  0.0004464     # Variance of permanent environment of litter size of dam
# quality variances
var.live.qual.gen       <- 0.12
var.live.qual.res       <- 0.41
# skin quality
var.skin.qual.gen       <- 1.53
var.skin.qual.res       <- 3.42
# guard hair length
var.h.length            <- 0.07
var.h.length.res        <- 0.41
# skin sizes
mean.skin.length.male   <- 95.59
var.skin.length.male    <- 11.6
var.skin.length.c.male  <- 3.14
var.skin.length.res.male<- 26
mean.skin.length.female <- 75.62
var.skin.length.female  <- 8.94
var.skin.length.c.female<- 1.49
var.skin.length.res.female <- 9.29


# misc
roof.body.size          <-  2400          # Maximum body weight of females to be allowed for selection
bw.eff.damage           <-  34.5          # effect of dam age on weight of males

################# RR parameters ###################

G_BWF <- as.matrix(read.table("G_BWF.txt"))
P_BWF <- as.matrix(read.table("P_BWF.txt"))
P_BWM <- as.matrix(read.table("P_BWM.txt"))
P_RFI <- as.matrix(read.table("P_RFI.txt"))
# new sigma matrix
G_sigma <- as.matrix(read.table("G_sigma.txt"),header=T, row.names = 1)
# make a vector of the variances to quicken calculations later
variances <-
  c(
    0.12, # live.qual
    0.07, #h.length
    1.53, #skin.qual
    11.6, #skin.length.male
    8.94, #skin.length.female
    0.01222, #litter.size
    0.187570056965907E-01, #bw1.f
    0.105996930807828E-01, #bw2.f
    0.108406600737242E-02, #bw3.f
    0.499379849225065E-01, #bw1.m
    0.341965058424113E-01, #bw2.m
    0.375707107280707E-02, #bw3.m
    0.583035959863543, # rfi1.m
    1.015749852,       # rfi2.m
    0.502027738625064, #rfi1.f
    0.851297595421713  #rfi2.f
  )
pe.var.bw.female <- c(0.463272591471770E-02, 0.110671108702861E-02, 0.295290173247113E-03)
pe.var.bw.male <-
  c(0.140998088558170E-01, 0.310587884472959E-02, 0.721913463410233E-03)
bw.res.male <- as.matrix(read.table(file="RES_BWM.txt"))
bw.res.male[8] <- bw.res.male[8] *32 # this is the adjustment needed to get the proper heritability
bw.res.female <- as.matrix(read.table(file="RES_BWF.txt"))
bw.res.female[8] <- bw.res.female[8]*54 # this is the adjustment needed to get the proper heritability
pe.var.rfi <- c(
  0.686214714336812,
  0.744536664027022 ,
  0.660456914599594, 
  0.467616097053697 
)

FR.RFI <- c(15.4, 11.1) # this is the adjusted fixed regression from Shirali et. al (2016)

# function to make standardized time, kept here since this is sourced before the functions
StandardTime <- function (t) {
  t1 <- 2*(t- min(t))/(max(t)-min(t))-1 
  return(t1)
}
# solutions for fixed regression BW
FR.females <- c(1, 0.3, -0.12)
FR.males <- c(1, .94, 1.05)
# standardized time vector for use in bw and rfi
t <- c(63,84,105,126,147,168,189,210)
t <- StandardTime(t)
t <- t[3:8]
# 2nd order legendre polynomials
leg2 <- legendre.polynomials(2, normalized = T)
leg1 <- legendre.polynomials(1, normalized = T)

######### regression coefficients
b.bw.male   <- 1.86
b.bw.female <- 1.43
res.rfi <-  as.matrix(read.table("RES_RFI.txt"))/2

################ Skin price constants ################
truncs <- c(-4.11, -0.47, 4.25) # for qual categories
              # 50 40 30 00 0 1 2
htruncs <- c(-2.21672259, -0.74, 0.457, 1.723) # hair length categories

intercept.bm <- 419 #intercept plus average of auctions
prices.bmales <- c(0,-38,-80,-140,-204,-219,-216,-227,-240,-250,119,103,75,44,0,56,22,-33,-25)
prices.bfemales <- c(200,195,185,175,137,92,43,-7,16,0,113,104,85,42,0,145,129,76,92)
intercept.bfemales <- 18 #intercept + average of auctions
