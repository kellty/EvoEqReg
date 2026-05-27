n_list <- c(100,200,350,700); l_n <- length(n_list)
sig_list <- (1:3)/10; l_sig <- length(sig_list)
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

# PDE regularization
library(splines2)
basis <- function(xt,flag=0,Kx=9,Kt=9) {
  knots_x <- (1:(Kx-4))/(Kx-3)
  knots_t <- (1:(Kt-4))/(Kt-3)
  if (flag==0) {# u
    b_x <- bSpline(xt[,1], knots=knots_x, degree=3, intercept=TRUE, Boundary.knots=c(0,1))
    b_t <- bSpline(xt[,2], knots=knots_t, degree=3, intercept=TRUE, Boundary.knots=c(0,1))
    t(sapply(1:nrow(xt), function(i) rep(b_x[i,],times=Kt)*rep(b_t[i,],each=Kx)))
  } else if (flag==1) {# u_t - u_xx
    b_x <- bSpline(xt[,1], knots=knots_x, degree=3, intercept=TRUE, Boundary.knots=c(0,1))
    b_t <- bSpline(xt[,2], knots=knots_t, degree=3, intercept=TRUE, Boundary.knots=c(0,1))
    b_x_2 <- deriv(b_x, 2)
    b_t_1 <- deriv(b_t)
    t(sapply(1:nrow(xt), function(i) 
      rep(b_x[i,],times=Kt)*rep(b_t_1[i,],each=Kx) - rep(b_x_2[i,],times=Kt)*rep(b_t[i,],each=Kx)))
  } else if (flag==-1) {# u_x
    b_x_1 <- deriv(bSpline(xt[,1], knots=knots_x, degree=3, intercept=TRUE, Boundary.knots=c(0,1)))
    b_t <- bSpline(xt[,2], knots=knots_t, degree=3, intercept=TRUE, Boundary.knots=c(0,1))
    t(sapply(1:nrow(xt), function(i) rep(b_x_1[i,],times=Kt)*rep(b_t[i,],each=Kx)))
  }
}
node <- (1:1000)/1001
b_bd0 <- basis(cbind(rep(0,1000),node), flag=-1)
b_bd1 <- basis(cbind(rep(1,1000),node), flag=-1)
pen_bd <- (t(b_bd0)%*%b_bd0 + t(b_bd1)%*%b_bd1) / 1000
pen_pde <- 0
for (x in node) {
  b_pde <- basis(cbind(rep(x,1000),node), flag=1)
  pen_pde <- pen_pde + t(b_pde)%*%b_pde
}
pen_pde <- pen_pde / 1000^2
b0 <- basis(cbind(node,rep(0,1000)))

# initial value
s <- 2  # smoothness of function
K0 <- 50
ini.coef <- c(0.3, sapply(2:K0, function(k){4*(-1)^(k-1)/k^s}))
ini.func <- function(x){ sapply(1:K0, function(k) psi(k,x)) %*% ini.coef }
ini_node <- c(ini.func(node))

for(idx_sig in 1:l_sig) {
  sig <- sig_list[idx_sig]
result <- matrix(nrow=l_n, ncol=4)
for(idx_n in 1:l_n) {
  n <- n_list[idx_n]
tmp <- matrix(nrow=MC, ncol=2)
for (mc in 1:MC) {
  set.seed(999+(mc+55)^2)
obs_x <- runif(n)
obs_t <- runif(n)
Z_ <- Z(K0, obs_x,obs_t, lambda,psi)
obs_u <- Z_ %*% ini.coef + rnorm(n, sd=sig)
ZTZ <- t(Z_) %*% Z_
ZTu <- t(Z_) %*% obs_u
# comparison
Kmax <- 5
tune_pde <- 10^(-8:3)
tune_bd <- 10^(-18:5)
tune <- cbind(rep(tune_pde,times=length(tune_bd)), rep(tune_bd,each=length(tune_pde)))
ISE_cut <- ISE_pen <- NULL
for (K in 1:Kmax) {
  ini.coef_hat <- c(solve(ZTZ[1:K,1:K], ZTu[1:K]))
  ISE_cut <- c(ISE_cut, sum((c(ini.coef_hat, rep(0,K0-K)) - ini.coef)^2))
}
b_xt <- basis(cbind(obs_x,obs_t))
bTb <- t(b_xt) %*% b_xt
bTu <- t(b_xt) %*% obs_u
for (idx_tune in 1:nrow(tune)) {
  pen <- tune[idx_tune,1]*pen_pde + tune[idx_tune,2]*pen_bd
  ini_hat_node <- b0 %*% solve(bTb + pen, bTu)
  ISE_pen <- c(ISE_pen, mean((ini_hat_node - ini_node)^2))
}
print(c(idx_sig, n, mc, which.min(ISE_cut), log10(tune[which.min(ISE_pen),])))
tmp[mc,] <- c(min(ISE_cut), min(ISE_pen))
}
result[idx_n,] <- c(mean(tmp[,1]), sqrt(var(tmp[,1])/MC),
                    mean(tmp[,2]), sqrt(var(tmp[,2])/MC))
}
pdf(paste0("simu_noise",idx_sig,".pdf"), height=6*1, width=6*1)
plot(n_list, result[,1], log='xy', type='b', lwd=2, main=bquote(sigma == .(sig)),
     xlab=expression(n), ylab='ISE') # expression(abs(abs(hat(g)-g))^2)
segments(n_list, result[,1]-result[,2], n_list, result[,1]+result[,2], lty=1)
arrows(n_list, result[,1]-result[,2], n_list, result[,1]+result[,2], code=3, angle=90, length=0.1)
lines(n_list, result[,3], type='b', lty=2, lwd=2)
segments(n_list, result[,3]-result[,4], n_list, result[,3]+result[,4], lty=1)
arrows(n_list, result[,3]-result[,4], n_list, result[,3]+result[,4], code=3, angle=90, length=0.1)
lines(1:max(n_list), 0.96*result[1,1]*(min(n_list)/1:max(n_list))^((2*s-1)/(2*s+r)), lty=6, lwd=2)
# expression(phantom() %prop% n^{-(2*s-1)/(2*s+r)})
legend('topright', c("proposed","penalized","theoretical"), lty=c(1,2,6), lwd=2)
dev.off()
}

