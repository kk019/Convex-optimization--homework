f <- function(x){
  -1 * sum(log(1 - A %*% x)) - 1 * sum(log(1 - x**2))
}
grad <- function(x){
  t(t(1 / (1 - A %*% x)) %*% A) + 2 * x / (1 - x**2)
}
hess <- function(x){
  t(A) %*% diag(1 / as.vector(1 - A %*% x)**2) %*% A - 2 * diag((x**2 + 1) / (1 - x**2)**2)
}


m <- 300
n <- 50
alpha <- 0.01
beta <- 0.5
gradtol <- 1e-3
maxiters <- 1e3



A <- matrix(rnorm(m*n), m, n)
# --------------------------------------------------------------------------------
# grad
x <- matrix(0, n, 1)
Val.gd <- c()
for (i in 1:maxiters){
  cat(i)
  val <- f(x)
  Val.gd <- c(Val.gd, val)
  g <- grad(x)
  if (sqrt(sum(g**2)) < gradtol){
    break
  }
  v <- -1 * g
  
  # initerior
  t <- 1
  while ((max(A %*% (x + t * v)) >= 1) | (max((x + t * v)**2) > 1)){
    t <- beta * t
  }
  # decrease
  while (f(x + t * v) > val - alpha * t * sum(g**2)){
    t <- beta * t
  }
  x <- x + t * v
}

plot(1:length(Val.gd), Val.gd, type = 'l')


# ---------------------------------------------------------------------------------------
# Newton
nttol <- 1e-8
x <- matrix(0, n, 1)
Val.nt <- c()
for (i in 1:maxiters){
  cat(i)
  val <- f(x)
  Val.nt <- c(Val.nt, val)
  g <- grad(x)
  # v <- -1 * solve(hess(x)) %*% g
  L <- chol(hess(x), lower=TRUE)
  v <- solve(L, solve(t(L), -1*g))
  
  fprime <- sum(g * v)
  if (abs(fprime) < nttol){
    break
  }
  
  # initerior
  t <- 1
  while ((max(A %*% (x + t * v)) >= 1) | (max((x + t * v)**2) > 1)){
    t <- beta * t
  }
  # decrease
  while (f(x + t * v) > val + alpha * t * fprime){
    t <- beta * t
  }
  x <- x + t * v
}

# -------------------------------------------------------------------------------
plot(1:length(Val.nt), Val.nt,  type = 'l')
