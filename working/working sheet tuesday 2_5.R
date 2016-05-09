################## BW phenotype RR function ##############
library(orthopolynom)
StandardTime <- function (t) {

 t1 <- 2*(t- min(t))/(max(t)-min(t))-1 
  
  return(t1)
}
t <- c(63,84,105,126,147,168,189,210)
t <- StandardTime(t)
leg2 <- legendre.polynomials(2, normalized = T)
r1 <- c(1, 0.3, -0.12)
r2 <- c(1,0.17,-0.26)
q <- as.matrix(as.data.frame(polynomial.values(polynomials = leg2, x =t)))

colnames(q) <- c("leg0", "leg1", "leg2")
q <- q[, 1:ncol(q)]

bw <- (q) %*% r
# bw <- (bw-0.707)
# fm <- lm(bw ~ q-1)
# summary(fm)
# make the phenotype for each instance
a <- as.matrix(add.gen[,7:9, with=FALSE])
pe <- as.matrix(perm.env.bw[,])
bw <- matrix(nrow=1000, ncol = 8)
FR.males <- c(1, .34, .24)

bw <- (q) %*% (FR.males)
plot(bw, type = "l")
for (i in 1:nrow(a)) {
  bw[i,] <- 0.9 + q %*% r1 + q %*% (a[i,]) + q %*% pe[i,] 
}
plot(bw[50,], type = "l")
lines(bw[60,])
lines(bw[100,])
plot(colMeans(bw), type = "l")


# generate phenotypes
q <- as.matrix(as.data.frame(polynomial.values(polynomials = leg2, x =t)))
bw <- 0.9 + q %*% r1 + q %*% (a[i,]) + q %*% pe[i,]
