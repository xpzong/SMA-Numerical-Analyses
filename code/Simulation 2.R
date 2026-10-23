rm(list=ls(all=TRUE))
library(TeachingSampling)
library(quadprog)
library(glmnet)
library(ncvreg)


M0<-101;    N<-2000;    n<-300;
beta<-sqrt(1.5)/(1:M0)
beta[52:M0]<-0

mn1<-function(D=1000, N=2000, pho=0.6, M=5){
  
	med<-10
	rs1<-matrix(0, D, med+2)
	rs2<-matrix(0, D, med+2)
	
	for(i in 1:D){
		
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
			hm2[m] <- sum(phi^2 * rowSums((Xm %*% Gphi) * Xm))		
      	
      
			Xmsp1<-xssp1[, 1:m]
			Aphisp1 <- crossprod(Xmsp1, phisp1 * Xmsp1)
			Gphisp1 <- chol2inv(chol(Aphisp1))
			Bmsp1 <- Gphisp1 %*% crossprod(Xmsp1, phisp1 * yssp1)
			

			Bhatp[1:m, m] <- Bmp	
			Bhatsp1[1:m, m] <- Bmsp1 
		}
    
		ebarp <- ys%*%t(rep(1,M))-xs%*%Bhatp	 
    

		what<-matrix(0, M, med)
    
		#AIC, SAIC
		AIC_value <- n*log(colSums(ebarp^2))+2*(1:M)
		what[which.min(AIC_value), 3] <- 1
		AIC_value <- AIC_value - min(AIC_value)
		w <- exp(-0.5*AIC_value)/sum(exp(-0.5*AIC_value))
		what[ , 4] <- as.matrix(w)
    
		#BIC, SBIC
		BIC_value <- n*log(colSums(ebarp^2))+log(n)*(1:M)
		what[which.min(BIC_value), 5] <- 1
		BIC_value <- BIC_value - min(BIC_value)
		w <- exp(-0.5*BIC_value)/sum(exp(-0.5*BIC_value))
		what[ , 6] <- as.matrix(w)
    
		#Cp           
		sig <- -1
		if(sig == -1){
		  sig <- sum((ebarp[ , M])^2)/(n-p)	
		}
		Cp_value <- colSums(ebarp^2)/sig+2*(1:M)-n
		what[which.min(Cp_value), 7] <- 1
    
    
		#SMA1          
		Dmat <- crossprod(ebarp, ebarp * phi)   #(t(ebarp)%*%diag(phi)%*%ebarp)
		dvec <- as.vector(-(1:M)*sig*N/n)
		Amat <- t(rbind(matrix(1, nrow=1, ncol=M), diag(M)))
		bvec <- rbind(1,matrix(0, nrow=M, ncol=1))
		QP <- solve.QP(Dmat+diag(1e-8, M), dvec, Amat, bvec, meq = 1)
		w <- QP$solution
		w <- w*(w>0)
		w <- w/sum(w)
		what[ , 8] <- as.matrix(w)
    
		#SMA2
		Dmat <- crossprod(ebarp, ebarp * phi)
		dvec <- as.vector(-sig*hm2)
		Amat <- t(rbind(matrix(1, nrow=1, ncol=M), diag(M)))
		bvec <- rbind(1,matrix(0, nrow=M, ncol=1))
		QP <- solve.QP(Dmat+diag(1e-8, M), dvec, Amat, bvec, meq = 1)
		w <- QP$solution
		w <- w*(w>0)
		w <- w/sum(w)
		what[ , 9] <- as.matrix(w)
    
		#full
		what[ , 10] <- c(rep(0, M-1),1)

    
	
		#AGG1 
		yhatsp2<-xs[sp2_inx,]%*%Bhatsp1
		Dmat <- (t(yhatsp2)%*%yhatsp2)
		dvec <- as.vector(t(yhatsp2)%*%as.matrix(ys[sp2_inx]))
		Amat <- t(rbind(matrix(-1, nrow=1, ncol=M), diag(M)))
		bvec <- c(-1, rep(0, M))
		QP <- solve.QP(Dmat+diag(1e-8, M), dvec, Amat, bvec, meq = 0)
		w <- QP$solution
		w  <- w*(w>0)		
		what[ , 1] <- as.matrix(w)
    
		#AGG2
		sig <- -1
		if(sig == -1){
		  sig <-  sum((yssp1-xssp1%*%Bhatsp1[,M])^2)/(spn1-p)
		}
		f_j <- sqrt(colMeans(yhatsp2^2))
		rho <-  sqrt(2* sig)  * sqrt( (log(M) + log(spn2)) / spn2 )   #sigma=1
		lambda_glmnet <- rho 
		fit <- glmnet(x = yhatsp2, y = ys[sp2_inx], alpha = 1, intercept = FALSE, standardize = FALSE, penalty.factor = f_j, lambda = lambda_glmnet)
		w <- as.vector(coef(fit))[-1]
		what[ , 2] <- as.matrix(w)
    


		cvK <- min(10L, n)
		fold_id <- sample(rep(seq_len(cvK), length.out = n))

		# SCAD 
		cv_scad <- cv.ncvreg(X = xs[,2:M], y = ys, penalty = "SCAD", family = "gaussian", fold = fold_id)
		coef_scad_min <- drop(coef(cv_scad))
		tol <- 1e-8
		active_slope <- which(abs(coef_scad_min[-1]) > tol)
		ind_scad <- c(1L, active_slope + 1L)
		X_scad <- xs[, ind_scad, drop = FALSE]
		Aphi_scad <- crossprod(X_scad, phi * X_scad)
		Gphi_scad <- chol2inv(chol(Aphi_scad))
		B_scad <- Gphi_scad %*% crossprod(X_scad, phi * ys)
    
		# adaptive lasso
		fitols <- lm(ys~xs[,2:M])
		beta_init <- coef(fitols)[-1] 
		gamma <- 1
		wt <- 1 / (abs(beta_init) + 1e-6)^gamma
		cv_adlasso <- cv.glmnet(x = xs[,2:M], y = ys, alpha = 1, family = "gaussian", penalty.factor = wt,	foldid = fold_id,	intercept = TRUE, standardize = FALSE)
		coef_adlasso <- drop(as.matrix(coef(cv_adlasso, s = "lambda.min")))
		tol <- 1e-8
		active_slope <- which(abs(coef_adlasso[-1]) > tol)
		ind_ada <- c(1L, active_slope + 1L)

		X_ada <- xs[, ind_ada, drop = FALSE]
		Aphi_ada <- crossprod(X_ada, phi * X_ada)
		Gphi_ada <- chol2inv(chol(Aphi_ada))
		B_ada <- Gphi_ada %*% crossprod(X_ada, phi * ys)
    
    

		muhat <- cbind(x[,ind_scad, drop = FALSE]%*%matrix(B_scad), 
					x[,ind_ada, drop = FALSE]%*%matrix(B_ada),
					x[,1:M, drop = FALSE]%*%Bhatp%*%what)

		muccp <- mu%*%t(rep(1,med+2))-muhat 
		yccp <- ys%*%t(rep(1,med+2))-muhat[inxs1,]
		Y_hat <- as.vector(t(phi)%*%yccp) + colSums(muhat)

    	rs1[i, ] <- ((Y_hat-sum(y))/N)^2
    	rs2[i, ] <- apply(muccp, 2, function(t){mean(t^2)})
    	print(i)
	}

	return(cbind(rs1,rs2))
}


#rs<-mn1(D=10, N=2000, pho=0.6, M=11)
#colnames(rs)<-rep(c("scad","adlasso",'AGG1','AGG2','AIC', 'SAIC', 'BIC', 'SBIC',  'Cp', 'SMA1', 'SMA2', 'full'),2)
#apply(rs, 2, mean)




t1=Sys.time()
for(pho in c(0.1, 0.6)){
	rs01<-1:24
	for(M in seq(11,101,10)){
		rs0<-mn1(D=1000, N=2000, pho, M)
		colnames(rs0)<-rep(c("scad","adlasso",'AGG1','AGG2','AIC', 'SAIC', 'BIC', 'SBIC',  'Cp', 'SMA1', 'SMA2', 'full'),2)
		rs01<-rbind(rs01, apply(rs0, 2, mean))
	}
	rs01<-rs01[-1, ]
	write.csv( rs01, paste(pho, '300sq.csv') )
}
Sys.time()-t1
