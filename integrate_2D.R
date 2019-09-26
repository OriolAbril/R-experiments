integrate2 <- function(f, lower, upper, ...) {
  # Bivariate integration via recursive use of R's integrate function
  # int f(z1,z2,...) dz1 dz2
  #
  # - f: function that we wish to integrate
  # - lower: vector of length(2) indicating the lower limits of integration for (z1,z2)
  # - upper: vector of length(2) indicating the upper limits of integration for (z1,z2)
  #
  # Example of usage
  # f  <-  function(z1,z2,sigma) dmvnorm(c(z1,z2),sigma=sigma)
  # integrate2(f, lower=c(-Inf,-Inf), upper=c(Inf,Inf), sigma=matrix(c(2,1,1,2),ncol=2))
  if (length(lower) == 1) lower <- rep(lower, 2)
  if (length(upper) == 1) upper <- rep(upper, 2)
  fint <- function(z1) {
    ans <- double(length(z1))
    for (i in 1:length(z1)) {
      f2 <- function(z2) sapply(z2, function(zz) f(z1[i], zz, ...))
      ans[i] <- integrate(f2, lower[2], upper[2])$value
    }
    return(ans)
  }
  integrate(fint, lower[1], upper[1])$value
}
