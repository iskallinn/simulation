# truncation points of skins
trunc <- qnorm(
  p = c(0.04, 0.41, 0.96),
  mean = mean(gen0.females$phenotype.skin.qual),
  sd = sqrt(var(
    gen0.females$phenotype.skin.qual
  )),
  lower.tail = TRUE
)
truncs

htruncs  <- qnorm(
  p = c(0.00102, 0.143955, 0.73731,0.99244),
  mean = mean(gen0.females$phenotype.h.length),
  sd = sqrt(var(
    gen0.females$phenotype.h.length
  )),
  lower.tail = TRUE
)


colnames( results) <- c("Gen",
"Gmean",
"Gvar",
"Fis",
"Obs.fert",
"mean.phenotype.bw.females",
"gen.value.bs",
"mean.phenotype.bw.males",
"bw.var",
"cor.bw.phenotype",
"skin.length.mean",
"skin.length.var",
"skin.qual.mean",
"skin.qual.var",
"cor.ls.own.to.ls",
"mated.females",
"barren.females",
"numb.false.sires",
"numb.kits",
"remating.perc",
"perc.single.mat",
"survived.kits",
"mean.gen.val.qual",
"var.gen.val.qual",
"avg.skin.length.male",
"avg.skin.length.female",
"feed.per.skin",
"permenv.var.bw",
"var.bw.female",
"skin.price")


colnames(results) <- c(
  "Gen",
  "Gmean",
  "Gvar",
  "Fis",
  "Obs.fert",
  "mean.phenotype.bw.females",
  "gen.value.bs",
  "mean.phenotype.bw.males",
  "bw.var",
  "cor.bw.to.blup",
  "cor.bw.phenotype",
  "skin.length.mean",
  "skin.length.var",
  "skin.qual.mean",
  "skin.qual.var",
  "cor.ls.blup",
  "cor.ls.own.to.ls",
  "mated.females",
  "barren.females",
  "numb.false.sires",
  "numb.kits",
  "remating.perc",
  "perc.single.mat",
  "survived.kits",
  "mean.gen.val.qual",
  "var.gen.val.qual",
  "cor.blup.qual.gen.val.qual",
  "cor.live.score.skin.qa",
  "cor.blup.qual.to.skin.qual",
  "avg.skin.length.male",
  "avg.skin.length.female", 
  "feed.per.skin",
  "permenv.var.bw",
  "female.bw.var",
  "skin.price")
)