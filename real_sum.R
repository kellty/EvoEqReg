x_frap <- rep(dx * (1:nx), ny*nt)
y_frap <- rep(rep(dy * (1:ny), each=nx), times=nt)
t_frap <- rep(dt * (1:nt), each=nx*ny)
uTu <- sum(u_frap^2)

library(plot3D)
pix <- 5

lim0 <- c(3e4,6e4)
lim1 <- c(-2e3,2e3)

pdf(paste0("real_o.pdf"), height=pix*2, width=pix*3)
par(mfrow=c(2,3))
for (it in c(1:3,10,50,100)) {
  idx <- ((it-1)*nxy+1):(it*nxy)
scatter3D(x_frap[idx], y_frap[idx], u_frap[idx], pch='.',
          xlab='x', ylab='y', zlab='u', zlim=lim0, clim=lim0, theta=-20, phi=20)
legend('bottomright', legend=bquote(t==.(it*dt)~s), box.col = "white")
}
dev.off()

pdf("real_fit.pdf", height=pix*5, width=pix*4)
par(mfrow=c(5,4))
for (it in c(1:3,10,50)) {
  idx <- ((it-1)*nxy+1):(it*nxy)
  idx_bd <- (it-1)*nxy + c(1,nx,nxy-nx+1,nxy)
scatter3D(x_frap[idx_bd], y_frap[idx_bd], u_frap[c(1,nx),c(1,ny),it], pch='.',
          xlab='x', ylab='y', zlab='u_cut', zlim=lim0, clim=lim0, theta=-20, phi=20,
          surf=list(x=matrix(x_frap[idx], nrow=nx),
                    y=matrix(y_frap[idx], nrow=nx),
                    z=ucut_frap[,,it]))
legend('bottomright', legend=bquote(t==.(it*dt)~s), box.col = "white")
scatter3D(x_frap[idx], y_frap[idx], u_frap[idx]-ucut_frap[idx], pch='.',
          xlab='x', ylab='y', zlab='u-u_cut', zlim=lim1, clim=lim1, theta=-20, phi=20)
legend('bottomright', legend=bquote(t==.(it*dt)~s), box.col = "white")
scatter3D(x_frap[idx_bd], y_frap[idx_bd], u_frap[c(1,nx),c(1,ny),it], pch='.',
          xlab='x', ylab='y', zlab='u_pen', zlim=lim0, clim=lim0, theta=-20, phi=20,
          surf=list(x=matrix(x_frap[idx], nrow=nx),
                    y=matrix(y_frap[idx], nrow=nx),
                    z=upen_frap[,,it]))
legend('bottomright', legend=bquote(t==.(it*dt)~s), box.col = "white")
scatter3D(x_frap[idx], y_frap[idx], u_frap[idx]-upen_frap[idx], pch='.',
          xlab='x', ylab='y', zlab='u-u_pen', zlim=lim1, clim=lim1, theta=-20, phi=20)
legend('bottomright', legend=bquote(t==.(it*dt)~s), box.col = "white")
}
dev.off()


RSS_cut <- RSS_pen <- RSS_0 <- rep(NA, nt)
for (it in 1:nt) {
  idx <- ((it-1)*nxy+1):(it*nxy)
  RSS_cut[it] <- sum((u_frap[idx]-ucut_frap[idx])^2)
  RSS_pen[it] <- sum((u_frap[idx]-upen_frap[idx])^2)
  RSS_0[it] <- sum((u_frap[idx]-mean(u_frap[idx]))^2)
}
