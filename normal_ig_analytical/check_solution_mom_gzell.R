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
V0 <- normid_V0
group_idxs <- list(c(1, 2, 3), c(4, 5))
for (group_idx in group_idxs) {
  Z_aux <- Z[, group_idx]
  zellner_V0 <- solve(t(Z_aux) %*% Z_aux)
  V0[group_idx, group_idx] <- zellner_V0 / length(group_idx)
}
g <- .3
a <- .1
b <- .1

cat("Only groups: exact first vs mom_gzell in mombf\n")
# use analytical solution ----------------
paste(log_p_y_given_gam(y, Z, V0, g, a, b)[1])

# use mombf solution ---------------------
paste(
  nlpMarginal(
    seq(m), y, Z, priorCoef=momprior(tau=g), priorGroup=groupzellnerprior(tau=g),
    priorVar=igprior(a,b), groups=c(1,1,1,4,4)
  )
)

cat("\nOnly single covariates: mom vs mom_gzell in mombf\n")
# use mombf solution ---------------------
paste(nlpMarginal(
  seq(m), y, Z, priorCoef=momprior(tau=g),
  priorGroup=momprior(tau=g), priorVar=igprior(a,b)
))
# use mombf solution ---------------------
paste(nlpMarginal(
  seq(m), y, Z, priorCoef=momprior(tau=g),
  priorGroup=groupzellnerprior(tau=g), priorVar=igprior(a,b)
))

### Second case ##########################
cat("\nSecond case\n")
# define parameters used -----------------
n <- 100
m <- 9
theta <- matrix(sample(c(0, 1), m, replace = TRUE))
Z <- matrix(rnorm(n * m), nrow = n, ncol = m)
y <- Z %*% theta + matrix(rnorm(n))
tau <- 0.3
taugroup <- 2
normid_V0 <- diag(m)*tau
V0 <- normid_V0
group_idxs <- list(c(1, 2), c(3, 4, 5, 6), c(7, 8, 9))
for (group_idx in group_idxs) {
  Z_aux <- Z[, group_idx]
  zellner_V0 <- solve(t(Z_aux) %*% Z_aux)
  V0[group_idx, group_idx] <- zellner_V0 * taugroup / length(group_idx)
}
g <- 1
a <- .1
b <- .1

# use analytical solution ----------------
paste(log_p_y_given_gam(y, Z, V0, g, a, b)[1])

# use mombf solution ---------------------
paste(
  nlpMarginal(
    seq(m), y, Z, priorCoef=momprior(tau=tau), priorGroup=groupzellnerprior(tau=taugroup),
    priorVar=igprior(a,b), groups=c(1,1,3,3,3,3,4,4,4)
  )
)


