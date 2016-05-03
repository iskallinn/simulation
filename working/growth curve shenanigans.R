############## Feed intake function ################
# this function fits a growth function to each kit to estimate 
# the amount of feed eaten during the growing season

growth_function_1 <- function( x, t ) {
  BW <- log((x)/(1+exp(1)^-(t-9.3)/1.388))
  return(exp(BW)/1000)^0.75)
}
test <- growth_function_1(test2[2],t1)
const %*% test
const <- 700/(144*4.13)*100*7/1000
growth_function_2 <- function( x, t ) {
  BW <- log((x)/(1+exp(1)^-(t-17.5)/2.682))
  return(exp(BW))
}

t1 <- c(seq(1:24))

FeedIntakeKits <- function (kitlist, weaning, per2, per3) {}
  t1 <- c(seq(1:28))
  sum(growth_function_1(per3,t1)+growth_function_2(per3,t1))+weaning

  test <- kit.list$phenotype.bw.sept-kit.list$phenotype.bw.june
  test2 <-kit.list$phenotype.bw.oct-kit.list$phenotype.bw.sept
  test3 <-kit.list$phenotype.bw.june
sum(((growth_function_1((test),t1)+growth_function_2((test2),t1)+(test3))/1000)^0.75*700/(144*4.13)*100*7/1000)

  +growth_function_2(per3,t1))+weaning
  summary(test3)
kit.list$FI <- ifelse( sex == 1, sum())

kit.list$test <- test

sum((growth_function_1(kit.list$test[1],t1)))
    