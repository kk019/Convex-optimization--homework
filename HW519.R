m <- 100
n <- 50
A <- vector(mode = 'list', length = m)
b <- vector(mode = 'list', length = m)
for (i in 1:m){
  tmp <- matrix(rnorm(n**2), n, n)
  A[[i]] <- tmp %*% t(tmp)
  tmp <- rnorm(n)
  t <- 1
  while(tmp %*% solve(A[[i]]) %*% tmp > 2){
    tmp <- tmp / 2
  }
  b[[i]] <- tmp
}

f1 <- function(x){
  r <- matrix(0, n, 1)
  for (i in 1:m){
    r <- r + c(t(x) %*% A[[i]] %*% x / 2 + sum(b[[i]] * x) + 1) * (A[[i]] %*% x + b[[i]])
  }
  r
}
f2 <- function(x){
  r <- matrix(0, n, n)
  for (i in 1:m){
    r <- r + (A[[i]] %*% x + b[[i]]) %*% t(A[[i]] %*% x + b[[i]])
  }
  r
}
f3 <- function(x){
  r <- 0
  for (i in 1:m){
    r <- r + t(x) %*% A[[i]] %*% x / 2 + sum(b[[i]] * x) + 1
  }
  r
}


beta <- 0.5
alpha <- 0.01
tol <- 0.1
maxiters <- 200
Val <- c()
x <- matrix(0, n, 1)
result <- optim(matrix(0, n, 1), f3, method="L-BFGS-B")
for (iter in 1:maxiters){
  cat(iter)
  val <- f3(x)
  Val <- c(Val, val)
  v <- -1 * solve(f2(x)) %*% f1(x)
  if (sqrt(sum(v**2)) < tol){
    break
  }
  while (f3(x + t * v) > val - alpha * t * sum(v**2)){
    t <- beta * t
  }
  x <- x + t * v
}
plot(1:length(Val), Val - result$value, type= 'l')


