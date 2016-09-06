# Skin analysis
setwd("C:/Users/Notandi/Dropbox/Projects/simulation of mink farm/Output/DMU analysis/")
skin.metrics.females <- read.table("skin_metrics_females", header=T, row.names = NULL)
skin.metrics.males <- read.table("skin_metrics_males", header=T, row.names = NULL)


test1 <- summaryBy(S50/numb.animals + 
                     S40/numb.animals +
                     S30/numb.animals+
                     S00/numb.animals+
                     S0/numb.animals+
                     S1/numb.animals+
                     S2/numb.animals+
                     S3/numb.animals+
                     S4/numb.animals+
                     S5/numb.animals+
                     purple/numb.animals+
                     platinum/numb.animals+
                     burgundy/numb.animals+
                     ivory/numb.animals+
                     vel3/numb.animals+
                     vel2/numb.animals+
                   vel1/numb.animals+
                   kl/numb.animals+
                  long.nap/numb.animals+
                    avg.price.females
                   ~ Gen, data=skin.metrics.females)
test1 <- summaryBy(S50/numb.animals + 
                     S40/numb.animals +
                     S30/numb.animals+
                     S00/numb.animals+
                     S0/numb.animals+
                     S1/numb.animals+
                     S2/numb.animals+
                     S3/numb.animals+
                     S4/numb.animals+
                     S5/numb.animals+
                     purple/numb.animals+
                     platinum/numb.animals+
                     burgundy/numb.animals+
                     ivory/numb.animals+
                     vel3/numb.animals+
                     vel2/numb.animals+
                     vel1/numb.animals+
                     kl/numb.animals+
                     long.nap/numb.animals+
                     avg.price.males
                   ~ Gen, data=skin.metrics.males)


write.table(test1, file="C:/Users/Notandi/Dropbox/Projects/Kopenhagen Fur/Analysis results/Scenarios/Crossmating/Summarized Results//prices_males", quote=F, row.names=F)
