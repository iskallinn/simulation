RFI <- function (kit.list, leg2, t) {
  temp <-
    kit.list[, c("id",
                 "sex",
                 "bw1_f",
                 "bw2_f",
                 "bw3_f",
                 "bw1_m",
                 "bw2_m",
                 "bw3_m",
                 "pe1.bw.f",
                 "pe2.bw.f",
                 "pe3.bw.f",
                 "pe1.bw.m",
                 "pe2.bw.m",
                 "pe3.bw.m"), with = FALSE]
  # scale variances
  q <- as.matrix(as.data.frame(polynomial.values(polynomials = leg2, x =t)))
  FR.females <- c(1, 0.3, -0.12)
  FR.males <- c(1, .94, 1.05)
  
  const.m <- q %*% FR.males
  plot(const.f, type = "l", main= "Fixed regression, females", ylab= "Weight [kg]", xlab="time period")
  const.f <- 0.9 + q %*% FR.females
  
  # this step makes the 6 needed BW for the RFI
  temp[,`:=`( 
  bw1 = ifelse(sex==1, const.m[1] + bw1_m*q[[1,1]]+bw2_m*q[[1,2]]+bw3_m*q[[1,3]]+
                 pe1.bw.m*q[[1,1]]+pe2.bw.m*q[[1,2]]+pe3.bw.m*q[[1,3]]+rnorm(nrow(temp))*sqrt(bw.res.male[3]),
               const.f[1] + bw1_f*q[[1,1]]+bw2_f*q[[1,2]]+bw3_f*q[[1,3]]+
                 pe1.bw.f*q[[1,1]]+pe2.bw.f*q[[1,2]]+pe3.bw.f*q[[1,3]]+rnorm(nrow(temp))*sqrt(bw.res.female[3])),
  bw2 = ifelse(sex==1, const.m[2] + bw1_m*q[[2,1]]+bw2_m*q[[2,2]]+bw3_m*q[[2,3]]+
                 pe1.bw.m*q[[2,1]]+pe2.bw.m*q[[2,2]]+pe3.bw.m*q[[2,3]]+rnorm(nrow(temp))*sqrt(bw.res.male[4]),
               const.f[2] + bw1_f*q[[2,1]]+bw2_f*q[[2,2]]+bw3_f*q[[2,3]]+
                 pe1.bw.f*q[[2,1]]+pe2.bw.f*q[[2,2]]+pe3.bw.f*q[[2,3]]+rnorm(nrow(temp))*sqrt(bw.res.female[4])),
  bw3 = ifelse(sex==1, const.m[3] + bw1_m*q[[3,1]]+bw2_m*q[[3,2]]+bw3_m*q[[3,3]]+
                 pe1.bw.m*q[[3,1]]+pe2.bw.m*q[[3,2]]+pe3.bw.m*q[[3,3]]+rnorm(nrow(temp))*sqrt(bw.res.male[5]),
               const.f[3] + bw1_f*q[[3,1]]+bw2_f*q[[3,2]]+bw3_f*q[[3,3]]+
                 pe1.bw.f*q[[3,1]]+pe2.bw.f*q[[3,2]]+pe3.bw.f*q[[3,3]]+rnorm(nrow(temp))*sqrt(bw.res.female[5])),
  bw4 = ifelse(sex==1, const.m[4] + bw1_m*q[[4,1]]+bw2_m*q[[4,2]]+bw3_m*q[[4,3]]+
                 pe1.bw.m*q[[4,1]]+pe2.bw.m*q[[4,2]]+pe3.bw.m*q[[4,3]]+rnorm(nrow(temp))*sqrt(bw.res.male[5]),
               const.f[4] + bw1_f*q[[4,1]]+bw2_f*q[[4,2]]+bw3_f*q[[4,3]]+
                 pe1.bw.f*q[[4,1]]+pe2.bw.f*q[[4,2]]+pe3.bw.f*q[[4,3]]+rnorm(nrow(temp))*sqrt(bw.res.female[5])),
  bw5 = ifelse(sex==1, const.m[5] + bw1_m*q[[5,1]]+bw2_m*q[[5,2]]+bw3_m*q[[5,3]]+
                 pe1.bw.m*q[[5,1]]+pe2.bw.f*q[[5,2]]+pe3.bw.m*q[[5,3]]+rnorm(nrow(temp))*sqrt(bw.res.male[6]),
               const.f[5] + bw1_f*q[[5,1]]+bw2_f*q[[5,2]]+bw3_f*q[[5,3]]+
                 pe1.bw.f*q[[5,1]]+pe2.bw.f*q[[5,2]]+pe3.bw.f*q[[5,3]]+rnorm(nrow(temp))*sqrt(bw.res.female[6])),
  bw6 = ifelse(sex==1, const.m[6] + bw1_m*q[[6,1]]+bw2_m*q[[6,2]]+bw3_m*q[[6,3]]+
                 pe1.bw.m*q[[6,1]]+pe2.bw.m*q[[6,2]]+pe3.bw.m*q[[6,3]]+rnorm(nrow(temp))*sqrt(bw.res.male[7]),
               const.f[6] + bw1_f*q[[6,1]]+bw2_f*q[[6,2]]+bw3_f*q[[6,3]]+
                 pe1.bw.f*q[[6,1]]+pe2.bw.f*q[[6,2]]+pe3.bw.f*q[[6,3]]+rnorm(nrow(temp))*sqrt(bw.res.female[7]))
  
)]  
  temp1 <- temp[,c("id","sex","bw1", "bw2", "bw3", "bw4", "bw5", "bw6"),with=F]
  temp <- merge(temp1, 
                kit.list[,
                         c("id",
                           "rfi1_m",
                           "rfi2_m",
                           "rfi1_f",
                           "rfi2_f",
                           "pe1.rfi.m",
                           "pe2.rfi.m",
                           "pe1.rfi.f",            
                           "pe2.rfi.f"  
                         ),with=F], by="id")
  leg1 <- legendre.polynomials(1, normalized = T)
  q <- as.matrix(as.data.frame(polynomial.values(polynomials = leg1, x =t)))
  const.rfi <- q %*% FR.RFI
  temp[,`:=`(
    RFI1 = ifelse( sex == 1, 
                   bw1*b.bw.male + const.rfi[1]+ 
                     q[[1,1]]*rfi1_m+q[[1,2]]*rfi2_m+ q[[1,2]]*pe1.rfi.m+q[[1,1]]*pe2.rfi.m+  
                     rnorm(nrow(temp))*sqrt(res.rfi[1]),
                   bw1*b.bw.female+ const.rfi[1]+
                     q[[1,1]]*rfi1_f+ q[[1,2]]*rfi2_f + q[[1,1]]*pe1.rfi.f+q[[1,2]]*pe2.rfi.f+
                     rnorm(nrow(temp))*sqrt(res.rfi[1])),
    RFI6 = ifelse( sex == 1, 
                   bw6*b.bw.male + const.rfi[6]+ 
                     q[[6,1]]*rfi1_m + q[[6,2]]*rfi2_m+ q[[6,1]]*pe1.rfi.m+q[[6,2]]*pe2.rfi.m+  
                     rnorm(nrow(temp))*sqrt(res.rfi[6]),
                   bw6*b.bw.female+ const.rfi[6]+
                     q[[6,1]]*rfi1_f+ q[[6,2]]*rfi2_f+ q[[6,1]]*pe1.rfi.f+q[[6,2]]*pe2.rfi.f+
                     rnorm(nrow(temp))*sqrt(res.rfi[6]))
  )]
  # TO DO
  # 3) increase the residual variance to decrease heritability, need to reestimate heritability 
  # with dmu until I get something that looks ok
  # 4) test the RFI, do I get "enough" feed usage, does it matter if I change the fixed regression curve
  # so it looks like it should do (or i think it should)
  stat <- summaryBy(RFI6 ~ sex, data = temp, FUN= c(mean))
  p1 <-as.matrix( stat[1,c(2,3,4,5,6,7), with=F])
  p2 <-as.matrix( stat[2,c(2,3,4,5,6,7), with=F])
  
  plot(t(p1),type="l")
}
