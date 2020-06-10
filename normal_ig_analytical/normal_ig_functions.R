# Analytical marginal posterior probability
log_p_y_given_gam <- function(y, Z, V0, g, a, b) {
  # returns the posterior probability of y given gamma
  # gamma is the model identifier, which must be taken
  # into account when setting inputs Z and V0
  n <- length(y)
  m <- nrow(V0)
  V <- t(Z) %*% Z + solve(V0) / g
  const <- (-n * log(pi) - m * log(g) + a * log(b)) / 2
  const_gam <- lgamma((n + a) / 2) - lgamma(a / 2)
  dets <- -(
    determinant(V0, logarithm = TRUE)$modulus
      + determinant(V, logarithm = TRUE)$modulus
  ) / 2
  last_term <- -(n + a) / 2 * log(b +
    t(y) %*% y -
    t(y) %*% Z %*% solve(V) %*% t(Z) %*% y)
  const + const_gam + dets + last_term
}

# formulae for the joint posterior probability
log_joint_p_1d <- function(theta, phi, y, Z, V0, g, a, b) {
  n <- length(y)
  m <- nrow(V0)
  const1 <- -(n + m) / 2 * log(2 * pi) - m / 2 * log(g)
  const2 <- (a / 2) * log(b / 2) - lgamma(a / 2)
  v0_det <- -determinant(V0, logarithm = TRUE)$modulus / 2
  phi_term <- -((n + m + a) / 2 + 1) * log(phi)
  yzt <- y - Z %*% theta
  last_term <- (
    -(t(yzt) %*% yzt + t(theta) %*% solve(V0) %*% theta / g + b) / (2 * phi)
  )
  const1 + const2 + v0_det + phi_term + last_term
}

joint_p_1d <- function(theta, phi, y, Z, V0, g, a, b) {
  exp(log_joint_p_1d(theta, phi, y, Z, V0, g, a, b))
}

log_post_theta_phi <- function(theta, phi, y, Z, V0, g, a, b) {
  V_1 <- t(Z) %*% Z + solve(V0) / g
  V <- solve(V_1)
  m <- V %*% t(Z) %*% y
  b_hat <- (t(y) %*% y + b - t(m) %*% V_1 %*% m) / 2
  a_hat <- (length(y) + a) / 2
  ig_term <- invgamma::dinvgamma(phi, a_hat, b_hat, log = TRUE)
  norm_term <- mvtnorm::dmvnorm(theta, m, phi * V, log = TRUE)
  norm_term + ig_term
}

make_group_zellner_matrix <- function(Z, V0, taugroup, group_indexs) {
  for (group_idx in group_indexs) {
    Z_aux <- Z[, group_idx]
    zellner_V0 <- solve(t(Z_aux) %*% Z_aux)
    V0[group_idx, group_idx] <- zellner_V0 * taugroup
  }
  return(V0)
}

# Orthogonal approx marginal posterior probability (gmom only)
log_p_y_given_gam_gmom <- function(y, Z, V0, g, a, b, group_indexs, p_j) {
  # returns the posterior probability of y given gamma
  # gamma is the model identifier, which must be taken
  # into account when setting inputs Z and V0
  n <- length(y)
  m <- nrow(V0)
  a_hat <- (n + a) / 2
  Zty <- t(Z) %*% y
  V0_1 <- solve(V0)
  S <- t(Z) %*% Z + V0_1 / g
  Sinv <- solve(S)
  b_hat <- b + t(y) %*% y - t(Zty) %*% Sinv %*% Zty
  const <- (-n * log(pi) - m * log(g) + a * log(b)) / 2
  const_gam <- lgamma(a_hat) - lgamma(a / 2)
  dets <- -(
    determinant(V0, logarithm = TRUE)$modulus
      + determinant(S, logarithm = TRUE)$modulus
  ) / 2
  m_vec <- Sinv %*% Zty
  last_term <- -a_hat * log(b_hat)
  norm_marg_like <- (const + const_gam + dets + last_term)[1]
  print(paste("gzell: ", norm_marg_like))
  V0_1 <- V0_1 / p_j
  nu = as.integer(a_hat * 2)
  print(paste("term1:", b_hat/(nu-2)))
  print(paste("b_hat:", b_hat))
  for (group_idx in group_indexs) {
    W_j <- V0_1[group_idx, group_idx]
    S_j <- Sinv[group_idx, group_idx]
    m_j <- m_vec[group_idx]
    trace <- sum(diag(S_j %*% W_j))
    mm <- (nu - 2)/b_hat*t(m_j)%*%W_j%*%m_j
    aux <- g*length(group_idx) / (length(group_idx) +2)
    W_aux <- W_j / aux
    # print((nu - 2)/b_hat*t(m_j)%*%W_aux%*%m_j)
    # print(aux)
    print(paste("trace:", trace[1]))
    print(paste("aux:", t(m_j)%*%W_j%*%m_j))
    print(log(trace[1] + mm[1]))
    norm_marg_like <- norm_marg_like + log(trace[1] + mm[1])
  }
  return(norm_marg_like)
}
