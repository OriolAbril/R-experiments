logistic_post_optimizer <- function(par, fn.post, niter, X, y, mu.beta, sd.beta) {
  logistic_post_helper <- function(beta) {
    fn.post(beta, X, y, mu.beta, sd.beta, 0)
  }
  logistic_post_helper_g <- function(beta) {
    fn.post(beta, X, y, mu.beta, sd.beta, 1)$g
  }
  optim_res <- optim(
    par, logistic_post_helper, logistic_post_helper_g,
    method="BFGS", control=list(maxit=niter)
  )
  return(drop(optim_res$par))
}

laplace_logistic_post <- function(par, fn_post, niter, X, y, mu.beta, sd.beta) {
  beta_hatish <- logistic_post_optimizer(par, fn_post, niter, X, y, mu.beta, sd.beta)
  f_grad_hess <- fn_post(beta_hatish, X, y, mu.beta, sd.beta, 2)
  f <- f_grad_hess$f; hess <- f_grad_hess$h
  det_h <- - det(hess)
  # print(det_h)
  # print(f)
  lapprox <- f + log(sqrt(2 * pi / det_h))
  return(lapprox)
}

list_vars <- function(p) {
  gammas <- gtools::permutations(n = 2, r = p, v = c(FALSE, TRUE), repeats.allowed = TRUE)
  return(gammas)
}

gamma_to_str <- function(gamma) {
  return(paste(seq_along(gamma)[gamma], collapse=","))
}
