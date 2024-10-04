# model setting
sig <- .2  #noise level
n <- 200    #sample size
Kmax <- 5
MC <- 200
BIC.MC <- ISE.MC <- matrix(nrow=MC, ncol=Kmax)
set.seed(117592)
for (mc in 1:MC) {
obs_x <- runif(n)
obs_t <- runif(n)
Z_ <- Z(K0, obs_x,obs_t, lambda,psi)
obs_u <- Z_ %*% ini.coef + rnorm(n, sd=sig)
ZTZ <- t(Z_) %*% Z_
ZTu <- t(Z_) %*% obs_u
# selection
BIC <- ISE <- NULL
for (K in 1:Kmax) {
  ini.coef_hat <- c(solve(ZTZ[1:K,1:K], ZTu[1:K]))
  RSS <- sum(obs_u^2) - sum(ZTu[1:K] * ini.coef_hat)
  BIC <- c(BIC, n * log(RSS / n) + log(n) * K)
  ISE <- c(ISE, sum((c(ini.coef_hat, rep(0,K0-K)) - ini.coef)^2))
}
BIC.MC[mc,] <- BIC
ISE.MC[mc,] <- ISE
}
print(rbind(colMeans(BIC.MC), sqrt((colMeans(BIC.MC^2)-colMeans(BIC.MC)^2)/((MC-1)*MC)),
            colMeans(ISE.MC), sqrt((colMeans(ISE.MC^2)-colMeans(ISE.MC)^2)/((MC-1)*MC))))

# comparison
sig_list <- (1:4)/10; l_sig <- length(sig_list)
n_list <- c(100,200,350,700);  l_n <- length(n_list)
Kmax <- 5
tune_pde <- 10^(-8:3)
tune_bd <- 10^(-18:5)
tune <- cbind(rep(tune_pde,times=length(tune_bd)), rep(tune_bd,each=length(tune_pde)))
MC <- 200
ISE_cut.MC <- ISE_pen.MC <- array(dim=c(l_n,l_sig,MC), dimnames=list(n_list,sig_list,1:MC))
result <- array(dim=c(l_n,l_sig,4), dimnames=list(n_list,sig_list,c('cut1','cut2','pen1','pen2')))
for(idx_sig in 1:l_sig) for(idx_n in 1:l_n) {
  n <- n_list[idx_n]
  sig <- sig_list[idx_sig]
set.seed(996)
for (mc in 1:MC) {
obs_x <- runif(n)
obs_t <- runif(n)
Z_ <- Z(K0, obs_x,obs_t, lambda,psi)
obs_u <- Z_ %*% ini.coef + rnorm(n, sd=sig)
ZTZ <- t(Z_) %*% Z_
ZTu <- t(Z_) %*% obs_u
ISE_cut <- ISE_pen <- NULL
for (K in 1:Kmax) {
  ini.coef_hat <- c(solve(ZTZ[1:K,1:K], ZTu[1:K]))
  ISE_cut <- c(ISE_cut, sum((c(ini.coef_hat, rep(0,K0-K)) - ini.coef)^2))
}
ISE_cut.MC[idx_n,idx_sig,mc] <- min(ISE_cut)
b_xt <- basis(cbind(obs_x,obs_t))
bTb <- t(b_xt) %*% b_xt
bTu <- t(b_xt) %*% obs_u
for (idx_tune in 1:nrow(tune)) {
  pen <- tune[idx_tune,1]*pen_pde + tune[idx_tune,2]*pen_bd
  ini_hat_node <- b0 %*% solve(bTb + pen, bTu)
  ISE_pen <- c(ISE_pen, mean((ini_hat_node - ini_node)^2))
}
ISE_pen.MC[idx_n,idx_sig,mc] <- min(ISE_pen)
# print(c(idx_sig, n, mc, which.min(ISE_cut), log10(tune[which.min(ISE_pen),])))
}
result[idx_n,idx_sig,] <- c(mean(ISE_cut.MC[idx_n,idx_sig,]),
                            sqrt(var(ISE_cut.MC[idx_n,idx_sig,])/MC),
                            mean(ISE_pen.MC[idx_n,idx_sig,]),
                            sqrt(var(ISE_pen.MC[idx_n,idx_sig,])/MC))
}
for (idx in 1:l_sig) {
pdf(paste0("simu_rand_",idx,".pdf"), height=6*1, width=6*1)
plot(n_list, result[,idx,1], log='xy', type='b', lwd=2, main=bquote(sigma == .(sig_list[idx])),
     xlab=expression(n), ylab='ISE') # expression(abs(abs(hat(g)-g))^2)
segments(n_list, result[,idx,1]-result[,idx,2], n_list, result[,idx,1]+result[,idx,2], lty=1)
arrows(n_list, result[,idx,1]-result[,idx,2], n_list, result[,idx,1]+result[,idx,2], code=3, angle=90, length=0.1)
lines(n_list, result[,idx,3], type='b', lty=2, lwd=2)
segments(n_list, result[,idx,3]-result[,idx,4], n_list, result[,idx,3]+result[,idx,4], lty=1)
arrows(n_list, result[,idx,3]-result[,idx,4], n_list, result[,idx,3]+result[,idx,4], code=3, angle=90, length=0.1)
lines(1:max(n_list), 0.96*result[1,idx,1]*(min(n_list)/1:max(n_list))^((2*s-1)/(2*s+r)), lty=6, lwd=2)
# expression(phantom() %prop% n^{-(2*s-1)/(2*s+r)})
legend('topright', c("proposed","penalized","theoretical"), lty=c(1,2,6), lwd=2)
dev.off()
}

