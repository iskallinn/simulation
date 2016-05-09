RFI <- function (kit.list, leg2,leg1 t) {
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
  const.f <- 0.9 + q %*% FR.females # this 0.9 is to get the body weight to a "realistic" level
  
  # this step makes the 6 needed BW for the RFI
  temp[,`:=`( 
    phenotype.bw.oct = ifelse(sex==1, const.m[6] + bw1_m*q[[6,1]]+bw2_m*q[[6,2]]+bw3_m*q[[6,3]]+
                   pe1.bw.m*q[[6,1]]+pe2.bw.m*q[[6,2]]+pe3.bw.m*q[[6,3]]+rnorm(nrow(temp))*sqrt(bw.res.male[7]),
                 const.f[6] + bw1_f*q[[6,1]]+bw2_f*q[[6,2]]+bw3_f*q[[6,3]]+
                   pe1.bw.f*q[[6,1]]+pe2.bw.f*q[[6,2]]+pe3.bw.f*q[[6,3]]+rnorm(nrow(temp))*sqrt(bw.res.female[7]))
    
  )]  
  temp <- temp[,c("id","sex","phenotype.bw.oct"),with=F]
  temp <- merge(temp, 
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
    FI = ifelse( sex == 1, 
                 phenotype.bw.oct*b.bw.male + const.rfi[6]+ 
                     q[[6,1]]*rfi1_m + q[[6,2]]*rfi2_m+ q[[6,1]]*pe1.rfi.m+q[[6,2]]*pe2.rfi.m+  
                     rnorm(nrow(temp))*sqrt(res.rfi[6]),
                 phenotype.bw.oct*b.bw.female+ const.rfi[6]+
                     q[[6,1]]*rfi1_f+ q[[6,2]]*rfi2_f+ q[[6,1]]*pe1.rfi.f+q[[6,2]]*pe2.rfi.f+
                     rnorm(nrow(temp))*sqrt(res.rfi[6]))
  )]
  kit.list <- merge(kit.list,temp[,c("id",
                  "phenotype.bw.oct",
                  "FI"), with=FALSE], by="id")
 return(kit.list) 
}