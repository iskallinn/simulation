## last year analysis and first year analysis
load(file="last.year.Rdata")
kit.list[, `:=`(hl.class = ifelse(
  phenotype.h.length > htruncs[4],
  0,
  # velv3
  ifelse(
    phenotype.h.length > htruncs[3] & phenotype.h.length < htruncs[4],
    1,
    # velv2
    ifelse(
      phenotype.h.length > htruncs[2] & phenotype.h.length < htruncs[3],
      2,
      # vel1
      ifelse(
        phenotype.h.length > htruncs[1] & phenotype.h.length < htruncs[2],
        3,
        # klassik
        ifelse(phenotype.h.length < htruncs[1], 4, NA) # long nap
      )
    )
  )
))]

kit.list[, `:=`(
  skin.qual.score = ifelse(
    phenotype.skin.qual >= trunc[3],
    4,
    # purple
    ifelse(
      truncs[2] < phenotype.skin.qual & phenotype.skin.qual <= truncs[3],
      3,
      # platinum
      ifelse(
        phenotype.skin.qual > truncs[1] & phenotype.skin.qual <= truncs[2],
        2,
        # burgundy
        ifelse(phenotype.skin.qual <=
                 truncs[1], 1, 0) # ivory
      )
    )
  ),
  skin.size = ifelse ( phenotype.skin.length >= 101,9,
                       ifelse(
                         phenotype.skin.length < 101 & phenotype.skin.length >= 95,
                         8,
                         ifelse(
                           phenotype.skin.length < 95 & phenotype.skin.length >= 89,
                           7,
                           ifelse (
                             phenotype.skin.length < 89 & phenotype.skin.length >= 83,
                             6,
                             ifelse(
                               phenotype.skin.length < 83 &
                                 phenotype.skin.length >= 77,
                               5,
                               ifelse(
                                 phenotype.skin.length < 77 & phenotype.skin.length >= 71,
                                 4,
                                 ifelse (
                                   phenotype.skin.length < 71 & phenotype.skin.length >= 65,
                                   3,
                                   ifelse(
                                     phenotype.skin.length < 65 & phenotype.skin.length >= 59,
                                     2,
                                     ifelse(
                                       phenotype.skin.length < 59 &
                                         phenotype.skin.length >= 53,
                                       1,
                                       ifelse(phenotype.skin.length < 53, 0, -999)
                                     )
                                   )
                                 )
                               )
                             )
                           )
                         )
                       )
  )
)]

hist(kit.list$skin.size)
hist(subset(kit.list, sex == 1)$skin.size)
hist(subset(kit.list, sex == 2)$skin.size)
hist(kit.list$skin.qual.score)
hist(kit.list$hl.class)
hist(kit.list$skin.price)
hist(kit.list$FI)
mean(kit.list$FI)
mean(kitlist$FI)


table(subset(kit.list, sex == 1)$skin.size)
table(subset(kit.list, sex == 2)$skin.size)
table(kit.list$skin.qual.score)
table(kit.list$hl.class)
table(kit.list$skin.price)
