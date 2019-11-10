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
# zellner_V0 <- solve(t(Z) %*% Z)
V0 <- normid_V0
# V0[group_idx, group_idx] <- zellner_V0[group_idx, group_idx]
V0[group_idx, group_idx] <- zellner_V0
g <- .3
a <- .1
b <- .1
print(as.character(determinant(V0, log=TRUE)$modulus))
print(as.character(determinant(t(Z) %*% Z + solve(V0)/g, log=FALSE)$modulus))

# use analytical solution ----------------
log_p_y_given_gam(y, Z, V0, g, a, b)[1]

# use mombf solution ---------------------
nlpMarginal(seq(m), y, Z, priorCoef=normalidprior(tau=g),
            priorGroup=groupzellnerprior(tau=g), priorVar=igprior(a,b), groups=c(1,2,3,4,4))
