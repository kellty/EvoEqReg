n_list <- c(200,700)
sig_list <- 0.2
s_list <- c(2,4)
Kmax <- 5
MC <- 200

# eigenpairs of -Laplacian in L2(0,1) with Neumann boundary conditions
r <- 2  # degree of operator / dimension of space
lambda <- function(k){ ((k-1)*pi)^r }
psi <- function(k,x) {
  if(k==1){ x^0 } else{ sqrt(2)*cos((k-1)*pi*x) }
}
Z <- function(K, obs_x, obs_t, lambda, psi) {
  sapply(1:K, function(k) exp(-lambda(k)*obs_t) * psi(k,obs_x))
}
# initial value
K0 <- 50
for (s in s_list) {# smoothness of function
ini.coef <- c(0.3, sapply(2:K0, function(k){4*(-1)^(k-1)/k^s}))
ini.func <- function(x){ sapply(1:K0, function(k) psi(k,x)) %*% ini.coef }
for (n in n_list) for (sig in sig_list) {
  print(paste('s',s,'n',n,'sig',sig))
result <- matrix(nrow=MC, ncol=2*Kmax)
for (mc in 1:MC) {
  set.seed(999+(mc+55)^2)
obs_x <- runif(n)
obs_t <- runif(n)
Z_ <- Z(K0, obs_x,obs_t, lambda,psi)
obs_u <- Z_ %*% ini.coef + rnorm(n, sd=sig)
ZTZ <- t(Z_) %*% Z_
ZTu <- t(Z_) %*% obs_u
BIC <- ISE <- NULL
for (K in 1:Kmax) {
  ini.coef_hat <- c(solve(ZTZ[1:K,1:K], ZTu[1:K]))
  RSS <- sum(obs_u^2) - sum(ZTu[1:K] * ini.coef_hat)
  BIC <- c(BIC, n * log(RSS / n) + log(n) * K)
  ISE <- c(ISE, sum((c(ini.coef_hat, rep(0,K0-K)) - ini.coef)^2))
}
result[mc,] <- c(BIC,ISE)
}
colnames(result) <- paste0(rep(c('BIC','ISE'),each=Kmax), 1:Kmax)
print(colMeans(result))
print(sqrt((colMeans(result^2)-colMeans(result)^2)/((MC-1)*MC)))
}}

