############### Packages to use     ###############
library(gamlss.dist) # for making zero inflated poisson random deviates for mating willingness
library(pedigree)    # for calculating inbreeding and trimming pedigree
library(data.table)  # faster then data frames
library(mvtnorm)     # generating random deviates from a multivariate normal distribution
library(doBy)

############## Index weights ############
# weights must sum to 1
weight.fert.old.females <- 0.35
weight.bw.old.females   <- 0.35
weight.qual.old.females <- 0.3
weight.fert.kits        <- 0.35
weight.bw.kits          <- 0.35
weight.qual.kits        <- 0.3
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
use.comb.ind.for.males <- 0 # if 1 then the usage of the males will depend on 
# their combined index, not quality

############### Controls for simulation ############
n.females <-  1000             # NUMBER OF FEMALES
nruns <- 1                    # how many replicates of the simulation
n <- 4                        # number of generation per replicate
mating.method <- assortative   # mating method, random or assortative
selection.method <- blup  # selection mating, 
weighing.method <- oct         # control for when to weigh kits for selection cands
qual.classes <- 10              # quality classes, 5 or 10 are supported
intensity.remating <- 0.2      # this controls how many males are chosen not to be remated
# selection = truncation,  no selection = next gen is chosen at random
make.obs.file <- 1 # 1 = make observation file, 0 otherwise
use.true.sire <- 1 # 1 if true sire of kits is wanted for BV prediction, 0 otherwise
############# Mean settings for traits   #####################
mean.body.size.male.oct <- 3000
mean.body.size.female.oct <- 1650
mean.body.size.male.sept <- 2750
mean.body.size.female.sept <- 1450
############### Innter settings, change at own risk##########
male.ratio <-  6               # MALE TO FEMALE RATIO
male.inf <-  0.98              # % ODDS OF MALE BEING NOT BARREN
female.inf <- 0.9              # % ODDS OF FEMALE BEING NOT BARREN
prop.oldfemales <-  0.4        # Proportion of older females
ibd.fertility <- 0.06          # Inbreeding depression on fertility NOTE: not in use atm
max.age.females <- 3           # define how old the females can be
yearling.effect <- -0.07       # setting for yearling effect, current best guess is -0.07
sib.effect.male <- -13.4       # Effect on body size of (male) one extra kit in litter, Hansen(1997)
sib.effect.female <- -18.6      # Effect on body size of (female) one extra kit in litter, Hansen(1997)
sib.effect.male.sept   <- -21.5 # Effect on body size of (female) one extra kit in litter, Hansen(1997)
sib.effect.female.sept <- -12.6 # Effect on body size of (female) one extra kit in litter, Hansen(1997)
quantile.setting <- 0.4        # amount of kits to throw away because of too low littersize
quantile.setting.bw <- 0.4    # prop of kits to disqualify due to weight
mating.will.yearling.1st          <- 0.95 # probability of yearling being mated
mating.will.yearling.2nd          <- 0.98 # probability of yearling being remated
mating.will.old.1st               <- 0.98 # probability of old female being mated
mating.will.old.2nd               <- 0.98 # probability of old female being remated
pr.barren.one.mating.yearling     <- 0.8 # probability of single mated yearling being barren
pr.barren.double.mating.yearling  <- 0.9 # probability of double mated yearling being barren
pr.barren.one.mating.old          <- 0.9 # probability of single mated old female being barren
pr.barren.double.mating.old       <- 0.95 # probaility of double mated old female being barren
n.males <-  ceiling( n.females/male.ratio ) # calculates needed amount of males 


############# Variance settings for traits ###################
# fertility
variance.fertility     <-  0.0122738     # Genetic variance of fertility, live born
var.perm.env.ls        <-  0.0004464     # Variance of permanent environment of litter size of dam
# october body weight
var.bw.oct.female      <-  26116          # Genetic variance of body size, direct effect 
var.bw.oct.male        <-  83389         # from hansen et al 1992, ass h^2 = 0.51
var.c.bw.oct.male      <-  13036          # common env. from hansen et al 1992, ass c^2 = 0.068
var.c.bw.oct.female    <-  3306          # from hansen et al 1992, ass c^2 = 0.068
var.res.bw.oct.female  <-  9040          # from hansen et al 1992, ass c^2 andh^2 as above
var.res.bw.oct.male    <-  39021         # from hansen et al 1992, ass c^2 andh^2 as above
# september body weight
var.bw.sept.female      <-  7101         # from hansen et al 1992, ass h^2 = 0.51
var.bw.sept.male        <-  22705        # from hansen et al 1992, ass h^2 = 0.51
var.c.bw.sept.male      <-  3027         # from hansen et al 1992, ass c^2 = 0.068
var.c.bw.sept.female    <-  947          # from hansen et al 1992, ass c^2 = 0.068
var.res.bw.sept.female  <-  5876         # from hansen et al 1992, ass c^2 andh^2 as above
var.res.bw.sept.male    <-  18788        # from hansen et al 1992, ass c^2 andh^2 as above
# quality variances
var.live.qual.gen       <- 24.0
var.live.qual.res       <- 60.0
# misc
roof.body.size          <-  2400          # Maximum body weight of females to be allowed for selection
bw.eff.damage           <-  34.5          # effect of dam age on weight of males
############# (Co)Variance matrix ############################
# sigma <-  matrix( 
#   c(1, -0.5, -0.5, 1), # genetic correlations of traits
#   nrow=2, 
#   ncol=2) 
# colnames(sigma) <- c("bw.oct","litter.size")
# rownames(sigma)<- c("bw.oct","litter.size")

# bigger sigma matrix for more traits
# values are found in "Correlation between the development of mink kits in the
# lactation and growth perios, correlations to fur properties and heritability
# estimation, by Hansen, Lohi & Berg, 1992
sigma <- matrix(
  c(  1,    0.84, -0.5, -0.08,   -0.41,  0.73,
     0.84,  1,    -0.28,   0,     -0.4,   0.7,
    -0.5, -0.28,   1,     0,      0,     0,
    -0.08, 0,     0,     1,      0.6,   0,
    -0.41, -0.4, 0,     0.6,     1,     0,
    0.73,  0.7,  0,     0,       0,     1),
  nrow=6, ncol=6, dimnames= list(c("bw.oct", "bw.sept", "litter.size", 
  "live.qual", "skin.qual", "skin.length" ),
  c("bw.oct", "bw.sept", "litter.size", "live.qual", "skin.qual", "skin.length" )))

sigma.new <- matrix( 
  # qual velv   wMale wFemale sk.qual s.sm   s.sf    ls wsept.m  wsept.fe mat.w weaning
  c( 1,   0.86, -0.13,  0.18,   0.74,  -0.15, 0.12,   0,  -0.13,   0.18, 0,   0,           # qual
    0.86, 1,     0,     0,      0.55,   0,    0,      0,   0,      0,    0,   0,           # velv
   -0.13, 0,     1,     0.83,  -0.52,   0.79, 0.77,  -0.5, 0.84,   0.83, 0,   0.3,           # wMale
    0.18, 0,     0.83,  1,     -0.38,   0.64, 0.8,    0,   0.83,   0.94,-0.2, 0.3,           # wFemale
    0.74, 0.55, -0.52, -0.38,   1,     -0.55,-0.42,   0,   0,      0,    0,  -0.25,           # sk.qual
   -0.15, 0,     0.79,  0.64,  -0.55,  1,     0.89,   0,   0.7,    0.7,  0,   0.45,           # sk.size.male
    0.12, 0,     0.77,  0.8,   -0.42,  0.89,  1,      0,   0.7,    0.7,  0,   0.45,           # sk.size.fema
    0,    0,     0,    -0.5,    0,     0,     0,      1,  -0.28,  -0.28, 0,   0,           # litter.size
   -0.13, 0,     0.84,  0.83,  0,     0.7,   0.7,   -0.28,1,       0.83, 0,   0.57,           # bw.sept
    0.18, 0,     0.84,  0.83,  0,     0.7,   0.7,   -0.28,0.83,    1,    0,   0.57,           # bw.sept
    0,    0,    0,    -0.2,     0,     0,     0,      0,   0,      0,    1,   0,           # maternal 
    0,    0,    0,    -0.2,     0,     0,     0,      0,   0,      0,    0,   1              # weaning 
    
    ), nrow=12, ncol=12)

var.maternal <- matrix(c(25.8, 45, 35, 23,16,28,22,14),nrow=4,ncol=2,dimnames = list(
        c("maternal","additive","specific", "residual"),c("male","female")))
