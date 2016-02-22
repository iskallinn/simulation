### To compare different selection methods
# it is a bit of a bother to make the data fit the ggplot stuff

summarized.all <- rbind(summarized,summarized.blup)
keeps <- c("Gen", "Gmean.mean","method") # choose what you want to keep
summarized.all <- summarized.all[keeps]

summarized.blup$method <- "blup"
summarized$method <- "trad"

temp1 <- subset(summarized.all, method=="blup")
keeps <- c("Gmean.mean")
temp1 <- temp1 [keeps]
names( temp1 ) [names( temp1 )   ==   'Gmean.mean'] <-  'blup'

temp2 <- subset(summarized.all, method=="trad")
temp2 <- temp2 [keeps]
names( temp2 ) [names( temp2 )   ==   'Gmean.mean'] <-  'trad'
y <- seq(1:14)

df <- cbind(temp1,temp2,y)


mdf <- melt (df, id = "y", value.name="Gmean.mean",variable.name="Gen")

ggplot(data=mdf,
       aes(x=y, y=value, colour=variable)) +
  geom_line()

########### mean bs of females

summarized.all <- rbind(summarized,summarized.blup)
keeps <- c("Gen", "mean.phenotype.bs.females.mean","method") # choose what you want to keep
summarized.all <- summarized.all[keeps]

summarized.blup$method <- "blup"
summarized$method <- "trad"

temp1 <- subset(summarized.all, method=="blup")
keeps <- c("mean.phenotype.bs.females.mean")
temp1 <- temp1 [keeps]
names( temp1 ) [names( temp1 )   ==   'mean.phenotype.bs.females.mean'] <-  'blup'

temp2 <- subset(summarized.all, method=="trad")
temp2 <- temp2 [keeps]
names( temp2 ) [names( temp2 )   ==   'mean.phenotype.bs.females.mean'] <-  'trad'
y <- seq(1:14)

df <- cbind(temp1,temp2,y)


mdf <- melt (df, id = "y", value.name="mean.phenotype.bs.females.mean",variable.name="Gen")

ggplot(data=mdf,
       aes(x=y, y=value, colour=variable)) +
  geom_line()

################## obs fert

summarized.all <- rbind(summarized,summarized.blup)
keeps <- c("Gen", "Obs.fert.mean","method") # choose what you want to keep
summarized.all <- summarized.all[keeps]

summarized.blup$method <- "blup"
summarized$method <- "trad"

temp1 <- subset(summarized.all, method=="blup")
keeps <- c("Obs.fert.mean")
temp1 <- temp1 [keeps]
names( temp1 ) [names( temp1 )   ==   'Obs.fert.mean'] <-  'blup'

temp2 <- subset(summarized.all, method=="trad")
temp2 <- temp2 [keeps]
names( temp2 ) [names( temp2 )   ==   'Obs.fert.mean'] <-  'trad'
y <- seq(1:14)

df <- cbind(temp1,temp2,y)


mdf <- melt (df, id = "y", value.name="Obs.fert.mean",variable.name="Gen")

ggplot(data=mdf,
       aes(x=y, y=value, colour=variable)) +
  geom_line()


################## mean size of males

summarized.all <- rbind(summarized,summarized.blup)
keeps <- c("Gen", "mean.phenotype.bs.males.mean","method") # choose what you want to keep
summarized.all <- summarized.all[keeps]

summarized.blup$method <- "blup"
summarized$method <- "trad"

temp1 <- subset(summarized.all, method=="blup")
keeps <- c("mean.phenotype.bs.males.mean")
temp1 <- temp1 [keeps]
names( temp1 ) [names( temp1 )   ==   'mean.phenotype.bs.males.mean'] <-  'blup'

temp2 <- subset(summarized.all, method=="trad")
temp2 <- temp2 [keeps]
names( temp2 ) [names( temp2 )   ==   'mean.phenotype.bs.males.mean'] <-  'trad'
y <- seq(1:14)

df <- cbind(temp1,temp2,y)


mdf <- melt (df, id = "y", value.name="mean.phenotype.bs.males.mean",variable.name="Gen")

ggplot(data=mdf,
       aes(x=y, y=value, colour=variable)) +
  geom_line()

################## inbreeding

summarized.all <- rbind(summarized,summarized.blup)
keeps <- c("Gen", "Fis.mean","method") # choose what you want to keep
summarized.all <- summarized.all[keeps]

summarized.blup$method <- "blup"
summarized$method <- "trad"

temp1 <- subset(summarized.all, method=="blup")
keeps <- c("Fis.mean")
temp1 <- temp1 [keeps]
names( temp1 ) [names( temp1 )   ==   'Fis.mean'] <-  'blup'

temp2 <- subset(summarized.all, method=="trad")
temp2 <- temp2 [keeps]
names( temp2 ) [names( temp2 )   ==   'Fis.mean'] <-  'trad'
y <- seq(1:14)

df <- cbind(temp1,temp2,y)


mdf <- melt (df, id = "y", value.name="Fis.mean",variable.name="Gen")

ggplot(data=mdf,
       aes(x=y, y=value, colour=variable)) +
  geom_line()


################## variance of fertility

summarized.all <- rbind(summarized,summarized.blup)
keeps <- c("Gen", "Gvar.mean","method") # choose what you want to keep
summarized.all <- summarized.all[keeps]

summarized.blup$method <- "blup"
summarized$method <- "trad"

temp1 <- subset(summarized.all, method=="blup")
keeps <- c("Gvar.mean")
temp1 <- temp1 [keeps]
names( temp1 ) [names( temp1 )   ==   'Gvar.mean'] <-  'blup'

temp2 <- subset(summarized.all, method=="trad")
temp2 <- temp2 [keeps]
names( temp2 ) [names( temp2 )   ==   'Gvar.mean'] <-  'trad'
y <- seq(1:14)

df <- cbind(temp1,temp2,y)


mdf <- melt (df, id = "y", value.name="Gvar.mean",variable.name="Gen")

ggplot(data=mdf,
       aes(x=y, y=value, colour=variable)) +
  geom_line()
