rm(list=ls(all=TRUE))
library(TeachingSampling)
library(quadprog)
library(glmnet)
library(ncvreg)

mn1<-function(D=1000, N=2000, n=300, pho=0.6, r=0.6){
	rs1<-1:D
	for(i in 1:D){

		M0<-500
		M<-11
		Mt<-6
		c2<-r/(1-r)
		beta<-sqrt(c2)/(1:M0)
		beta[(Mt+1):M0]<-0

		x<-matrix(rnorm(N*(M0-1)), N, M0-1)
		x<-cbind(1,x)
		e<-rnorm(N)
		mu<-x%*%beta
		y<-mu+e
    
		z<-abs(as.vector(pho/sqrt((1-pho^2)*var(y)))*y+rnorm(N,10,1))
		z<-pmax(z, 1e-8)
    


		S1<-S.piPS(n, z)
		inxs1<-S1[ ,1]
		pis1<-S1[ ,2]
		
		ys<-y[inxs1]
		xs<-x[inxs1, 1:M]
    
    
		phi<-1/pis1 
		
		spn1<-round(0.5*n)
		sp1_inx <- sample.int(n, size = spn1, replace = FALSE)
		sp2_inx <- setdiff(seq_len(n), sp1_inx)
		spn2<-length(sp2_inx)

		xssp1<-xs[sp1_inx, ]
		yssp1<-ys[sp1_inx]
		phisp1<-phi[sp1_inx]


		p<-M	
		Bhatp <- matrix(0, p, M) 
		hm2<-1:M

		Bhatsp1 <- matrix(0, p, M)

		for(m in 1:M){
			Xm <- xs[ ,1:m]
			Aphi <- crossprod(Xm, phi * Xm)
			Gphi <- chol2inv(chol(Aphi))
			Bmp <- Gphi %*% crossprod(Xm, phi * ys)
			Bhatp[1:m, m] <- Bmp	
		}
    
		ebarp <- ys%*%t(rep(1,M))-xs%*%Bhatp	 
		

		#SMA1          
		sig <- -1
		if(sig == -1){
		  sig <- sum((ebarp[ , M])^2)/(n-p)  	
		}
		Dmat <- crossprod(ebarp, ebarp * phi)
		dvec <- as.vector(-(1:M)*sig*N/n)
		Amat <- t(rbind(matrix(1, nrow=1, ncol=M), diag(M)))
		bvec <- rbind(1,matrix(0, nrow=M, ncol=1))
		QP <- solve.QP(Dmat+diag(1e-8, M), dvec, Amat, bvec, meq = 1)
		w <- QP$solution
		w <- w*(w>0)
		w <- w/sum(w)
		what<- as.matrix(w)
		rs1[i] <- sum(what[Mt:M])
		print(i)
	}
	mean(rs1)
}


t1=Sys.time()
nx<-c(seq(30,500,50), seq(500,2000,100))
rs<-nx
for(pho in c(0.1, 0.6)){
	for(r in c(0.1, 0.3, 0.5, 0.7, 0.9)){
		rs0<-1:length(nx)
		for(i in 1:length(nx)){
			rs0[i]<-mn1(D=1000, N=2000, n=nx[i], pho, r)
		}
		rs<-cbind(rs, rs0)
	}
}
Sys.time()-t1

dat <- as.data.frame(rs)
colnames(dat)<-c('n', 'A01', 'A03', 'A05', 'A07', 'A09', 'B01', 'B03', 'B05', 'B07', 'B09')
cols <- c("#0072B2", "#D55E00", "#009E73", "#A020F0", "#000000")
ltys <- c(1, 2, 3, 4, 5)
pchs <- c(16, 17, 15, 18, 8)


draw_panel <- function(Y, rho_value) {
  matplot(
    dat$n, Y,
    type = "o",
    col = cols,
    lty = ltys,
    pch = pchs,
    lwd = 2,
    cex = 0.65,
    xlab = expression(n),
    ylab = "Sum of Weights on Correct Models",
    main = bquote(rho == .(rho_value)),
    xlim = c(0, 2000),
    ylim = c(0.2, 1),
    xaxt = "n",
    yaxt = "n",
    las = 1,
    cex.lab = 0.95
  )
abline(h = 1, col = "gray40", lty = 2, lwd = 1.5)
  axis(1, at = seq(0, 2000, 400))
  axis(2, at = seq(0, 1, 0.2), las = 1)

  legend(
    "bottomright",
    legend = expression(
      R^2 == 0.1,
      R^2 == 0.3,
      R^2 == 0.5,
      R^2 == 0.7,
      R^2 == 0.9
    ),
    col = cols,
    lty = ltys,
    pch = pchs,
    lwd = 2,
    cex = 0.85,
    bty = "n"
  )
}

plot_both <- function() {
  old_par <- par(
    mfrow = c(1, 2),
    mar = c(4.2, 4.6, 2.5, 1),
    mgp = c(2.8, 0.8, 0)
  )
  on.exit(par(old_par))

  draw_panel(as.matrix(dat[, 2:6]),  rho_value = 0.1)
  draw_panel(as.matrix(dat[, 7:11]), rho_value = 0.6)
}

plot_both()


