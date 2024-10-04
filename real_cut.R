x_frap <- rep(dx * (1:nx), ny*nt)
y_frap <- rep(rep(dy * (1:ny), each=nx), times=nt)
t_frap <- rep(dt * (1:nt), each=nx*ny)
uTu <- sum(u_frap^2)

# eigenpairs in L2([0,xmax]×[0,ymax]) with periodic (Born-von Karman) boundary conditions
lambda_x <- function(kx){ Dpara * (ceiling(kx/2)*pi/xmax)^2 }
lambda_y <- function(ky){ Dpara * (ceiling(ky/2)*pi/ymax)^2 }
psi_x <- function(kx,x) {
  if (kx==0) { x^0 }
  else if (kx%%2) { sqrt(2) * sin((kx+1)/2*pi*x/xmax) }
  else { sqrt(2) * cos(kx/2*pi*x/xmax) }
}
psi_y <- function(ky,y) {
  if (ky==0) { y^0 }
  else if (ky%%2) { sqrt(2) * sin((ky+1)/2*pi*y/ymax) }
  else { sqrt(2) * cos(ky/2*pi*y/ymax) }
}
Z <- function(k_list, obs_x, obs_y, obs_t, lambda_x, psi_x, lambda_y, psi_y) {
  apply(k_list, 1, function(k)
    exp(-(lambda_x(k[1])+lambda_y(k[2]))*obs_t) * psi_x(k[1],obs_x)*psi_y(k[2],obs_y))
}

# estimation (cutoff)
Kxmax <- Kymax <- 20
k_list <- cbind(rep(0:(Kxmax-1),times=Kymax), rep(0:(Kymax-1),each=Kxmax))
Z_ <- Z(k_list, x_frap,y_frap,t_frap, lambda_x,psi_x,lambda_y,psi_y)
ZTZ <- t(Z_) %*% Z_
ZTu <- t(Z_) %*% c(u_frap)
BIC_cut <- NULL
for (idx_K in 1:nrow(k_list)) {
  Kx <- k_list[idx_K,1]+1;  Ky <- k_list[idx_K,2]+1
  kx <- c(rep(TRUE,Kx), rep(FALSE,Kxmax-Kx))
  ky <- c(rep(TRUE,Ky), rep(FALSE,Kymax-Ky))
  ks <- rep(kx,times=Kymax) & rep(ky,each=Kxmax)
  RSS <- uTu - sum(ZTu[ks] * c(solve(ZTZ[ks,ks], ZTu[ks])))
  BIC_cut <- c(BIC_cut, n_frap * log(RSS / n_frap) + log(n_frap) * Kx*Ky)
}
idx_K <- which.min(BIC_cut)
Kx <- k_list[idx_K,1]+1;  Ky <- k_list[idx_K,2]+1
kx <- c(rep(TRUE,Kx), rep(FALSE,Kxmax-Kx))
ky <- c(rep(TRUE,Ky), rep(FALSE,Kymax-Ky))
ks <- rep(kx,times=Kymax) & rep(ky,each=Kxmax)
ini.coef <- c(solve(ZTZ[ks,ks], ZTu[ks]))
ucut_frap <- array(Z_[,ks] %*% ini.coef, dim=nxyt)

