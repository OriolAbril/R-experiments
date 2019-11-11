m <- 7
vec <- matrix(rnorm(m))
mat <- matrix(rnorm(m^2), nrow = m, ncol = m)
print(drop(t(vec) %*% mat %*% vec))
print(determinant(mat, log=FALSE))

new_order <- matrix(sample(seq(m), m, replace = FALSE))

vec_reordered <- vec[new_order]
mat_reordered <- mat[new_order, new_order]
print(drop(t(vec_reordered) %*% mat_reordered %*% vec_reordered))
print(determinant(mat_reordered, log=FALSE))
