library(R.matlab)
frap_post <- NULL
for (i in c(paste0(0,1:9),10:20)) {# post-bleach
  frap_post[[i]] <- readMat(paste0("frap/frap_32ww_0",i,".mat"))$experiment[[3]]
}
nxyt <- dim(frap_post[[01]][[1]])
nx <- nxyt[1];  ny <- nxyt[2];  nt <- nxyt[3]
n_frap <- nx * ny *nt;  nxy <- nx * ny
u_frap <- array(0, dim=c(nx,ny,nt))
for (i in c(paste0(0,1:9),10:20)) {
  u_frap <- u_frap + frap_post[[i]][[1]] / 20
}
dx <- as.numeric(frap_post[[01]][[2]])
dy <- as.numeric(frap_post[[01]][[3]])
dt <- as.numeric(frap_post[[01]][[7]])
xmax <- dx * (nx+1)
ymax <- dy * (ny+1)
tmax <- dt * (nt+1)
rm(list='frap_post')

Dpara <- 8.9e-11  # diffusion coefficient

save.image("real_prep.RData")
