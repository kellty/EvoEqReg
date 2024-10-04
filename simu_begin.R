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
