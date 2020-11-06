###############################################################
#############   PROGRAM FOR THE ESTIMATION OF  ################
######    IMPACTS IN NONLINEAR SAR MODELS    ##################
###############################################################


##### PKG
Packages <- c('spatialreg', 'effects', 'mgcv', 'sp', 'spdep', 'Matrix','MASS','truncnorm', 'stringr','MCMCpack')
lapply(Packages, library, character.only = TRUE)


##### DATA
data(boston)

DATA <- boston.c
NEIG <-  nb2listw(boston.soi, style="W")


##### FREQUENTIST ESTIMATION
DATA$SPLINES <- bs(DATA$DIS)

LR  <- lm(CMEDV  ~  SPLINES + CRIM + TAX + AGE + INDUS, data = DATA)
SAR  <- lagsarlm(LR, data = DATA, listw = NEIG)


##### BAYESIAN ESTIMATION

#### SETUP 
MAT <- listw2mat(NEIG) 
W_D <- MAT # interaction matrix

n <- nrow(DATA)
DATA$CONST <- c(rep(1, n)) # cst


#### ESTIMATION (MCMC)
Iterations <- 10000
m = ncol(as.matrix(DATA[, c('CONST', 'SPLINES','CRIM', 'TAX', 'AGE', 'INDUS')])) # number of exogeneous variables including the cst

### A PRIORI DISTRIBUTIONS FOR PARAMETERS
# FOR beta (normal distribution)
beta0 <- matrix(c(rep(0, m)))
var0 <- 10e+10*diag(m)

# FOR sigma (gamma inverse distribution)
a0 <- 0
b0 <- 0

# FOR rho: no known distribution form

### MATRIX FOR THE RESULT OF ITERATIONS
theta    <- matrix(nrow = Iterations, ncol = m + 2)  # + 2 is for rho and sigma2


### STARTING VALUES
cur.beta <- matrix(c(rep(0.5, m)))  
cur.sigma2 <- 0.5
cur.rho <- 0.1

acc.prob <- 0  # initialize acceptance probability counter (for rho)


### DATA
X <- DATA[, c("CONST", "SPLINES", "CRIM", "TAX", "AGE", "INDUS")]
X <- as.matrix(X)

Y <- DATA[, c("CMEDV")]


CG <- solve(t(X)%*%X + solve(var0), tol = 1e-40)
JG <- solve(var0)%*%beta0 


for (i in 1:Iterations){ 
  
  ## beta
  A <- (diag(n) - cur.rho*W_D)
  k1 <- CG                        # variance
  k2 <- k1%*%(t(X)%*%A%*%Y + JG)  # mean
  
  cur.beta <- mvrnorm(1,  k2,  cur.sigma2*k1)
  
  ## sigma
  a = a0 + n/2
  b = b0 + (t(A%*%Y - X%*%cur.beta)%*%(A%*%Y - X%*%cur.beta))/2
  
  cur.sigma2 <- rinvgamma(1, a, b)
  
  ## rho
  # Posterior distribution: with current value
  Q1 <- (diag(n) - cur.rho*W_D)
  Q2 <- log(det(Q1))
  
  Q3 <- t(A%*%Y - X%*%cur.beta)%*%(A%*%Y - X%*%cur.beta)
  Q4 <- -1/(2*cur.sigma2)
  Q5 <- Q4*Q3
  Q6 <- Q2 + Q5
  
  
  
  # Posterior distribution: with proposed value
  propo.rho  <- 0.02 # tuning parameter (help to control convergence rate. One should test many values)
  prop.rho <- rtruncnorm(1, -1, 1, cur.rho, propo.rho)
  
  A1 <- (diag(n) - prop.rho*W_D)
  
  P1 <- (diag(n) - prop.rho*W_D)
  P2 <- log(det(P1))
  
  P3 <- t(A1%*%Y - X%*%cur.beta)%*%(A1%*%Y - X%*%cur.beta)
  P4 <- -1/(2*cur.sigma2)
  P5 <- P4*P3
  P6 <- P2 + P5
  
  
  Test <- P6 - Q6
  
  
  ask <- runif(1)
  ask1 <- log(ask)
  if( ask1 < Test) {
    cur.rho <- prop.rho	 
    acc.prob <- acc.prob + 1
  }
    
  theta[i,] <- c(cur.beta, cur.sigma2, cur.rho) 
  
} 


RATE <- acc.prob/Iterations  # acceptation rate 



#### ANALYSIS OF THE RESULTS (PARAMETERS OF THE SAR MODEL)
MCMC <- theta

q1 <- 0.1*Iterations  # burn-in: discard the first 10%
MCMC.FINAL <- theta[(q1+1):Iterations, ]  

mcmc.median     <- apply(MCMC.FINAL,2, median) 
mcmc.mean       <- apply(MCMC.FINAL,2, mean) 
S.tD            <- apply(MCMC.FINAL,2, sd) 
Q2.5            <- apply(MCMC.FINAL, 2, quantile, probs=0.025) 
Q97.5           <- apply(MCMC.FINAL, 2, quantile, probs=0.975)
MCerror         <- S.tD/sqrt(0.1*Iterations)

TAB_SUM <- cbind(mcmc.mean, S.tD, MCerror, Q2.5, mcmc.median, Q97.5)
colnames(TAB_SUM) <- c("MEAN", "STD", "MCError", "Q2.5", "MEDIAN", "Q97.5")

rownames(TAB_SUM) <- c("CONST", "SPLINES1", "SPLINES2", "SPLINES3", "CRIM", "TAX", "AGE", "INDUS", "SIGMA_SQUARE", "RHO")



#### CONVERGENCE CHECK
png("CONVERGENCE.png",
    units="in",
    width = 10, 
    height =8,
    res = 400,
    pointsize=15)

par(mfrow = c(4,3), par(mai = c(0.7, 0.7, 0.2, 0.5))) 

plot(MCMC.FINAL[,1],  type="l", col = "darkblue", xlab = "", ylab = "CONST")
plot(MCMC.FINAL[,2],  type="l", col = "darkblue", xlab = "", ylab = "SPLINES1 (DIST)")
plot(MCMC.FINAL[,3],  type="l", col = "darkblue", xlab = "", ylab = "SPLINES2 (DIST)")

plot(MCMC.FINAL[,4],  type="l", col = "darkblue", xlab = "", ylab = "SPLINES3 (DIST)")
plot(MCMC.FINAL[,5],  type="l", col = "darkblue", xlab = "", ylab = "CRIM")
plot(MCMC.FINAL[,6],  type="l", col = "darkblue", xlab = "", ylab = "TAX")

plot(MCMC.FINAL[,7],  type="l", col = "darkblue", xlab = "", ylab = "AGE")
plot(MCMC.FINAL[,8],  type="l", col = "darkblue", xlab = "ITERATIONS", ylab = "INDUS")
plot(MCMC.FINAL[,9],  type="l", col = "darkblue",  xlab = "ITERATIONS", ylab = "SIGMA_SQUARE")

plot(MCMC.FINAL[,10],  type="l", col = "darkblue",  xlab = "ITERATIONS", ylab = "RHO")
dev.off()






#### COMPUTATION OF EFFECTS (DIRECT AND INDIRECT)

DATA$SPLINES <- bs(DATA$DIS,3, Boundary.knots = c(min(DATA$DIS), max(DATA$DIS)))
DATA$EPSILON <- DATA$DIS + 10e-10
DATA$SP <- predict(DATA$SPLINES, DATA$EPSILON)

XP1 <- DATA[, c("CONST", "SPLINES", "CRIM", "TAX", "AGE", "INDUS")]
XP2 <- DATA[, c("CONST", "SP", "CRIM", "TAX", "AGE", "INDUS")]

XP1 <- as.matrix(XP1)
XP2 <- as.matrix(XP2)



# SIMULATED EFFECTS
DIRECT_S   <- matrix(nrow = n, ncol = Iterations)
INDIRECT_S <- matrix(nrow = n, ncol = Iterations)
TOTAL_S    <- matrix(nrow = n, ncol = Iterations)


W    <- as(as_dgRMatrix_listw(NEIG), "CsparseMatrix")
yobs <- Y


for (i in 1:Iterations){ 
  
  
  cur.beta    <- theta[i, 1:8] 
  cur.sigma2  <- theta[i, 9]
  cur.rho     <- theta[i, 10]
  
  
  # BASE 1
  prdX_1   <- as.vector(XP1%*%cur.beta)
  
  
  # BASE 2 (with epsilon)
  prdX_2   <- as.vector(XP2%*%cur.beta)
  
  
 
  MAT_INV1 <- powerWeights(W_D, rho = cur.rho, X=as.matrix(diag(n)))
  MAT_INV2 <- MAT_INV1%*%diag(n)
  
  
  DERIVEE_1 <- (prdX_2 - prdX_1)/10e-10
  
  
  
  # DIRECT EFFECTS
  DIRECT_S[, i] <- (diag(MAT_INV2)*DERIVEE_1) 
  
  
  # TOTAL EFFECTS
  TOTAL_S[, i] <- (MAT_INV2%*%DERIVEE_1) 
  
  # INDIRECT EFFECTS
  INDIRECT_S[, i] <- TOTAL_S[, i] - DIRECT_S[, i]
  
} 




#### DESCRIPTIVE STATISTICS ON EFFECTS
# DIRECT
DIRECT_EFF_S   <- DIRECT_S[, (q1+1):Iterations]

Median_DS      <- apply(DIRECT_EFF_S,1, median) 
Mean_DS        <- apply(DIRECT_EFF_S,1, mean) 
S.tD_DS        <- apply(DIRECT_EFF_S,1, sd) 
Q2.5_DS        <- apply(DIRECT_EFF_S, 1, quantile, probs=0.025) 
Q97.5_DS       <- apply(DIRECT_EFF_S, 1, quantile, probs=0.975)
MCerror_DS     <- S.tD_DS/sqrt(0.1*Iterations)

DIRECT_EFFECT_S <- cbind(Mean_DS, S.tD_DS, MCerror_DS, Q2.5_DS, Median_DS, Q97.5_DS)


# INDIRECT
INDIRECT_EFF_S <- INDIRECT_S[, (q1+1):Iterations]

Median_IS      <- apply(INDIRECT_EFF_S,1, median) 
Mean_IS        <- apply(INDIRECT_EFF_S,1, mean) 
S.tD_IS        <- apply(INDIRECT_EFF_S,1, sd) 
Q2.5_IS        <- apply(INDIRECT_EFF_S, 1, quantile, probs=0.025) 
Q97.5_IS       <- apply(INDIRECT_EFF_S, 1, quantile, probs=0.975)
MCerror_IS     <- S.tD_IS/sqrt(0.1*Iterations)

INDIRECT_EFFECT_S <- cbind(Mean_IS, S.tD_IS, MCerror_IS, Q2.5_IS, Median_IS, Q97.5_IS)


# TOTAL
TOTAL_EFF_S <- TOTAL_S[, (q1+1):Iterations]

Median_TS      <- apply(TOTAL_EFF_S,1, median) 
Mean_TS        <- apply(TOTAL_EFF_S,1, mean) 
S.tD_TS        <- apply(TOTAL_EFF_S,1, sd) 
Q2.5_TS        <- apply(TOTAL_EFF_S, 1, quantile, probs=0.025) 
Q97.5_TS       <- apply(TOTAL_EFF_S, 1, quantile, probs=0.975)
MCerror_TS     <- S.tD_IS/sqrt(0.1*Iterations)

TOTAL_EFFECT_S <- cbind(Mean_TS, S.tD_TS, MCerror_TS, Q2.5_TS, Median_TS, Q97.5_TS)





#### GRAPHICAL REPRESENTATION
png("EFFECTS_OK.png",
    units="in",
    width = 10, 
    height = 10,
    res = 400,
    pointsize=15)
par(mfrow = c(2, 2), par(mai = c(0.9, 1.3, 0.3, 0.1))) 


matplot(DATA$DIS, DIRECT_EFFECT_S[, c("Q2.5_DS", "Mean_DS", "Q97.5_DS")], type="p", pch = 1, cex=0.5, ylab = "Individual direct impact", xlab = "")
legend("topleft", inset=.05, legend=c("Q2.5", "Mean", "Q97.5"), pch="-",  col = c("black", "red", "green"), horiz=FALSE)


matplot(DATA$DIS, INDIRECT_EFFECT_S[, c("Q2.5_IS", "Mean_IS", "Q97.5_IS")], type="p", pch = 1, cex=0.5, ylab = "Individual indirect impact", xlab = "")
legend("topleft", inset=.05, legend=c("Q2.5", "Mean", "Q97.5"), pch="-",  col = c("black", "red", "green"), horiz=FALSE)





mod_D1=loess(DIRECT_EFFECT_S[,"Q2.5_DS"] ~ DATA$DIS, span=0.3)
mod_D2=loess(DIRECT_EFFECT_S[, "Mean_DS"] ~ DATA$DIS, span=0.3)
mod_D3=loess(DIRECT_EFFECT_S[, "Q97.5_DS"] ~ DATA$DIS, span=0.3)


DATA$yfit_D1 = predict(mod_D1, newdata=DATA$DIS)
DATA$yfit_D2 = predict(mod_D2, newdata=DATA$DIS)
DATA$yfit_D3 = predict(mod_D3, newdata=DATA$DIS)

DATA1 <- DATA[order(DATA$DIS), ]
matplot(DATA1$DIS, DATA1[, c("yfit_D1", "yfit_D2", "yfit_D3")], type="l", pch = 1, cex=0.5, ylab = "Individual direct impact", xlab="DIS")
legend("topleft", inset=.05, legend=c("Q2.5", "Mean", "Q97.5"), pch="-",  col = c("black", "red", "green"), horiz=FALSE)




mod_I1=loess(INDIRECT_EFFECT_S[,"Q2.5_IS"] ~ DATA$DIS, span=0.3)
mod_I2=loess(INDIRECT_EFFECT_S[, "Mean_IS"] ~ DATA$DIS, span=0.3)
mod_I3=loess(INDIRECT_EFFECT_S[, "Q97.5_IS"] ~ DATA$DIS, span=0.3)

DATA$yfit_I1 = predict(mod_I1, newdata=DATA$DIS)
DATA$yfit_I2 = predict(mod_I2, newdata=DATA$DIS)
DATA$yfit_I3 = predict(mod_I3, newdata=DATA$DIS)

DATA2 <- DATA[order(DATA$DIS), ]
matplot(DATA2$DIS, DATA2[, c("yfit_I1", "yfit_I2", "yfit_I3")], type="l", pch = 1, cex=0.5, ylab = "Individual indirect impact", xlab="DIS")
legend("topleft", inset=.05, legend=c("Q2.5", "Mean", "Q97.5"), pch="-",  col = c("black", "red", "green"), horiz=FALSE)

dev.off()

