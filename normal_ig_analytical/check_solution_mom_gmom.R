library("mombf")
source("normal_ig_functions.R")

set.seed(14)

# define parameters used -----------------
n <- 100
m <- 5
theta <- matrix(sample(c(0, 1), m, replace = TRUE))
print(t(theta))
Z <- matrix(rnorm(n * m), nrow = n, ncol = m)
y <- Z %*% theta + matrix(rnorm(n))
Zy <- t(Z) %*% y
taugroup <- .8
g <- 1
a <- 1
b <- 1
group_idxs <- list(c(1, 2, 3), c(4, 5))
p_j <- c(3,3,3,2,2)
V0 <- make_group_zellner_matrix(Z, diag(m), taugroup, group_idxs) * n / (p_j + 2)

cat("Only groups: orthogonal approx first vs mom_gmom in mombf\n")
# use analytical solution ----------------
mom_gmom <- log_p_y_given_gam_gmom(y, Z, V0, g, a, b, group_idxs, p_j)[1]
print("mom_gmom")
paste("R    :", mom_gmom)

# use mombf solution ---------------------
paste(
  "mombf:",
  nlpMarginal(
    seq(m), y, Z, priorCoef=momprior(tau=g), priorGroup=groupmomprior(tau=taugroup),
    priorVar=igprior(a,b), groups=c(1,1,1,4,4)
  )
)

cat("\nOnly single covariates: mom vs mom_gmom in mombf\n")
# use mombf solution ---------------------
paste(nlpMarginal(
  seq(m), y, Z, priorCoef=momprior(tau=g),
  priorGroup=momprior(tau=g), priorVar=igprior(a,b)
))
# use mombf solution ---------------------
paste(nlpMarginal(
  seq(m), y, Z, priorCoef=momprior(tau=g),
  priorGroup=groupmomprior(tau=g), priorVar=igprior(a,b)
))

### Second case ##########################
cat("\nSecond case\n")
# define parameters used -----------------
n <- 100
m <- 9
theta <- matrix(sample(c(0, 1), m, replace = TRUE))
print(t(theta))
Z <- matrix(rnorm(n * m), nrow = n, ncol = m)
y <- Z %*% theta + matrix(rnorm(n))
tau <- 0.3
taugroup <- 1
group_idxs <- list(c(1, 2), c(3, 4, 5, 6), c(7, 8, 9))
p_j <- c(2,2,4,4,4,4,3,3,3)
V0 <- make_group_zellner_matrix(Z, diag(m), taugroup, group_idxs) * n / (p_j + 2)
g <- 1
a <- 1
b <- 1

# use analytical solution ----------------
paste("R    :", log_p_y_given_gam_gmom(y, Z, V0, g, a, b, group_idxs, p_j)[1])

# use mombf solution ---------------------
paste(
  "mombf:",
  nlpMarginal(
    seq(m), y, Z, priorCoef=momprior(tau=tau), priorGroup=groupmomprior(tau=taugroup),
    priorVar=igprior(a,b), groups=c(1,1,3,3,3,3,4,4,4)
  )
)

paste("R gzell:", log_p_y_given_gam(y, Z, V0/p_j*(p_j+2), g, a, b)[1])

paste(
  "mombf gzell:",
  nlpMarginal(
    seq(m), y, Z, priorCoef=momprior(tau=tau), priorGroup=groupzellnerprior(tau=taugroup*n),
    priorVar=igprior(a,b), groups=c(1,1,3,3,3,3,4,4,4)
  )
)


