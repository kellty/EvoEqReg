x_frap <- rep(dx * (1:nx), ny*nt)
y_frap <- rep(rep(dy * (1:ny), each=nx), times=nt)
t_frap <- rep(dt * (1:nt), each=nx*ny)
uTu <- sum(u_frap^2)

xnode <- xmax * (1:1000)/1001
ynode <- ymax * (1:1000)/1001
tnode <- tmax * (1:1000)/1001

# PDE regularization
library(splines2)
basis <- function(xyt,flag=0,Kx=9,Ky=9,Kt=9) {
  knots_x <- xmax * (1:(Kx-4))/(Kx-3)
  knots_y <- ymax * (1:(Ky-4))/(Ky-3)
  knots_t <- tmax * (1:(Kt-4))/(Kt-3)
  if (flag==0) {# u
    b_x <- bSpline(xyt[,1], knots=knots_x, degree=3, intercept=TRUE, Boundary.knots=c(0,xmax))
    b_y <- bSpline(xyt[,2], knots=knots_y, degree=3, intercept=TRUE, Boundary.knots=c(0,ymax))
    b_t <- bSpline(xyt[,3], knots=knots_t, degree=3, intercept=TRUE, Boundary.knots=c(0,tmax))
    t(sapply(1:nrow(xyt), function(i) 
      rep(b_x[i,],times=Ky*Kt)*rep(rep(b_y[i,],each=Kx),times=Kt)*rep(b_t[i,],each=Kx*Ky)))
  } else if (flag==1) {# u_t - D (u_xx + u_yy)
    b_x <- bSpline(xyt[,1], knots=knots_x, degree=3, intercept=TRUE, Boundary.knots=c(0,xmax))
    b_y <- bSpline(xyt[,2], knots=knots_y, degree=3, intercept=TRUE, Boundary.knots=c(0,ymax))
    b_t <- bSpline(xyt[,3], knots=knots_t, degree=3, intercept=TRUE, Boundary.knots=c(0,tmax))
    b_x_2 <- deriv(b_x, 2)
    b_y_2 <- deriv(b_y, 2)
    b_t_1 <- deriv(b_t)
    t(sapply(1:nrow(xyt), function(i) 
      rep(b_x[i,],times=Ky*Kt)*rep(rep(b_y[i,],each=Kx),times=Kt)*rep(b_t_1[i,],each=Kx*Ky)
      - Dpara * rep(b_x_2[i,],times=Ky*Kt)*rep(rep(b_y[i,],each=Kx),times=Kt)*rep(b_t[i,],each=Kx*Ky)
      - Dpara * rep(b_x[i,],times=Ky*Kt)*rep(rep(b_y_2[i,],each=Kx),times=Kt)*rep(b_t[i,],each=Kx*Ky)))
  } else if (flag==-1) {# u_x
    b_x_1 <- deriv(bSpline(xyt[,1], knots=knots_x, degree=3, intercept=TRUE, Boundary.knots=c(0,xmax)))
    b_y <- bSpline(xyt[,2], knots=knots_y, degree=3, intercept=TRUE, Boundary.knots=c(0,ymax))
    b_t <- bSpline(xyt[,3], knots=knots_t, degree=3, intercept=TRUE, Boundary.knots=c(0,tmax))
    t(sapply(1:nrow(xyt), function(i) 
      rep(b_x_1[i,],times=Ky*Kt)*rep(rep(b_y[i,],each=Kx),times=Kt)*rep(b_t[i,],each=Kx*Ky)))
  } else if (flag==-2) {# u_y
    b_x <- bSpline(xyt[,1], knots=knots_x, degree=3, intercept=TRUE, Boundary.knots=c(0,xmax))
    b_y_1 <- deriv(bSpline(xyt[,2], knots=knots_y, degree=3, intercept=TRUE, Boundary.knots=c(0,ymax)))
    b_t <- bSpline(xyt[,3], knots=knots_t, degree=3, intercept=TRUE, Boundary.knots=c(0,tmax))
    t(sapply(1:nrow(xyt), function(i) 
      rep(b_x[i,],times=Ky*Kt)*rep(rep(b_y_1[i,],each=Kx),times=Kt)*rep(b_t[i,],each=Kx*Ky)))
  }
}

library(foreach)
library(doParallel)
cl <- makeCluster(100)
registerDoParallel(cl)
pen_bd <- foreach(t_=tnode,.combine='+') %dopar% {
  library(splines2)
  b_bd_x0 <- basis(cbind(rep(xmax,1000),ynode,rep(t_,1000))) - basis(cbind(rep(0,1000),ynode,rep(t_,1000)))
  b_bd_x1 <- basis(cbind(rep(xmax,1000),ynode,rep(t_,1000)), flag=-1) - basis(cbind(rep(0,1000),ynode,rep(t_,1000)), flag=-1)
  b_bd_y0 <- basis(cbind(xnode,rep(ymax,1000),rep(t_,1000))) - basis(cbind(xnode,rep(0,1000),rep(t_,1000)))
  b_bd_y1 <- basis(cbind(xnode,rep(ymax,1000),rep(t_,1000)), flag=-2) - basis(cbind(xnode,rep(0,1000),rep(t_,1000)), flag=-2)
  pen_bdx <- (t(b_bd_x0)%*%b_bd_x0 + t(b_bd_x1)%*%b_bd_x1) * ymax
  pen_bdy <- (t(b_bd_y0)%*%b_bd_y0 + t(b_bd_y1)%*%b_bd_y1) * xmax
  pen_bdx + pen_bdy
}
pen_bd <- pen_bd * tmax / 1001^2
pen_pde <- 0
for (t_ in tnode) {
pen_pde <- pen_pde + foreach(y=ynode,.combine='+') %dopar% {
  library(splines2)
  b_pde <- basis(cbind(xnode,rep(y,1000),rep(t_,1000)), flag=1)
  t(b_pde)%*%b_pde
}}
pen_pde <- pen_pde * xmax*ymax*tmax/1001^3
stopCluster(cl)

b_frap <- basis(cbind(x_frap,y_frap,t_frap))
bTu <- bTb <- 0
for (it in 1:nt) {
  idx <- ((it-1)*nxy+1):(it*nxy)
  bTu <- bTu + t(b_frap[idx,]) %*% u_frap[idx]
  bTb <- bTb + t(b_frap[idx,]) %*% b_frap[idx,]
}

# estimation (penalized)
tune_pde <- 10^(-6:0)
tune_bd <- 10^(-9:-3)
tune <- cbind(rep(tune_pde,times=length(tune_bd)), rep(tune_bd,each=length(tune_pde)))
GCV_pen <- NULL
for (idx_tune in 1:nrow(tune)) {
  pen <- tune[idx_tune,1]*pen_pde + tune[idx_tune,2]*pen_bd
  # H <- b_xyt %*% solve(t(b_xyt)%*%b_xyt + pen, t(b_xyt))
  trH <- sum(diag(solve(bTb + pen, bTb)))
  b.coef <- solve(bTb + pen, bTu)
  RSS <- uTu + t(b.coef)%*%bTb%*%b.coef - 2*t(b.coef)%*%bTu
  # sum((u_frap - b_xyt %*% solve(bTb + pen, bTu))^2)
  GCV_pen <- c(GCV_pen, n_frap * RSS / (n_frap - trH)^2)
}
idx_tune <- which.min(GCV_pen)
b.coef <- solve(bTb + tune[idx_tune,1]*pen_pde + tune[idx_tune,2]*pen_bd, bTu)
upen_frap <- array(b_frap %*% b.coef, dim=nxyt)

