library("mombf")
source("normal_ig_functions.R")

set.seed(14)

# define parameters used -----------------
n <- 100
m <- 5
theta <- matrix(sample(c(0, 1), m, replace = TRUE))
Z <- matrix(rnorm(n * m), nrow = n, ncol = m)
y <- Z %*% theta + matrix(rnorm(n))
normid_V0 <- diag(m)
group_idx <- c(4,5)
Z_aux <- Z[, group_idx]
zellner_V0 <- solve(t(Z_aux) %*% Z_aux)
V0 <- normid_V0
V0[group_idx, group_idx] <- zellner_V0
g <- .3
a <- .1
b <- .1
print(paste("det(V0)", determinant(V0, log=TRUE)$modulus, sep=" "))
print(paste("det(S)", determinant(t(Z) %*% Z + solve(V0)/g, log=FALSE)$modulus, sep=" "))

# use analytical solution ----------------
log_p_y_given_gam(y, Z, V0, g, a, b)[1]

# use mombf solution ---------------------
nlpMarginal(seq(m), y, Z, priorCoef=normalidprior(tau=g),
            priorGroup=groupzellnerprior(tau=g), priorVar=igprior(a,b), groups=c(1,2,3,4,4))

fit <- modelSelection(y,Z, priorCoef=normalidprior(tau=g), priorGroup=groupzellnerprior(tau=g), priorVar=igprior(a,b), groups=c(1,2,3,4,4))
print(postProb(fit))

### Second case ##########################
cat("\nSecond case\n")
# define parameters used -----------------
n <- 100
m <- 7
theta <- matrix(sample(c(0, 1), m, replace = TRUE))
Z <- matrix(rnorm(n * m), nrow = n, ncol = m)
y <- Z %*% theta + matrix(rnorm(n))
tau <- 0.3
taugroup <- 2
normid_V0 <- diag(m)*tau
V0 <- normid_V0
group_idx1 <- c(3,4,5)
Z_aux1 <- Z[, group_idx1]
zellner_V01 <- solve(t(Z_aux1) %*% Z_aux1)
V0[group_idx1, group_idx1] <- zellner_V01*taugroup
group_idx2 <- c(6,7)
Z_aux2 <- Z[, group_idx2]
zellner_V02 <- solve(t(Z_aux2) %*% Z_aux2)
V0[group_idx2, group_idx2] <- zellner_V02*taugroup
g <- 1
a <- .1
b <- .1
print(paste("det(V0)", determinant(V0/c(tau,tau,taugroup,taugroup,taugroup,taugroup,taugroup), log=TRUE)$modulus, sep=" "))
print(paste("det(S)", determinant(t(Z) %*% Z + solve(V0)/g, log=FALSE)$modulus, sep=" "))

# use analytical solution ----------------
log_p_y_given_gam(y, Z, V0, g, a, b)[1]

# use mombf solution ---------------------
nlpMarginal(seq(m), y, Z, priorCoef=normalidprior(tau=tau),
            priorGroup=groupzellnerprior(tau=taugroup), priorVar=igprior(a,b), groups=c(1,2,3,3,3,4,4))

fit <- modelSelection(y, Z, priorCoef=normalidprior(tau=tau),
            priorGroup=groupzellnerprior(tau=taugroup), priorVar=igprior(a,b), groups=c(1,2,3,3,3,4,4))
print(postProb(fit))
