source("../integrate_2D.R")
source("normal_ig_functions.R")

set.seed(1234)

###########################
# define parameters used  #
###########################
n <- 10
m <- 1
y <- matrix(rnorm(n))
Z <- matrix(rnorm(n * m), nrow = n, ncol = m)
V0 <- diag(m)
g <- 2
a <- 6
b <- 4

###########################
# use analytical solution #
###########################
exp(log_p_y_given_gam(y, Z, V0, g, a, b)[1])

###########################
# use numerical approach  #
###########################
# use numerical integration (in 2D) to check analytical result
integrate2(joint_p_1d, lower = c(-Inf, 0), upper = c(Inf, Inf), y = y, Z = Z, V0 = V0, g = g, a = a, b = b)[1]
