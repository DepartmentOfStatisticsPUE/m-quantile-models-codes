glm.mq.binom <- function (x, y, case.weights = rep(1, nrow(x)), maxit = 100, acc = 1e-04,var.weights = rep(1, nrow(x)),weights.x=FALSE,q = 0.5,k=1.345)
{
	require(MASS)
	
#Stopping rule
	irls.delta <- function(old, new) abs(max(old-new))/abs(max(old))
	if (qr(x)$rank < ncol(x))
	stop("X matrix is singular")
	if (length(case.weights) != nrow(x))
	stop("Length of case.weights must equal number of observations")
	if (any(case.weights < 0))
	stop("Negative case.weights are not allowed")
	n<-length(case.weights)
	ifelse (weights.x,w.x <- sqrt(1-hat(x)),w.x <- rep(1,length=n))
	assign("w.x",w.x)
#We fit the glm.rob for computing the starting values
	
	temp.rob <-glm.rob.binom (X=x,y=y,weights.on.x=weights.x,chuber=k)
	resid.init <- y-temp.rob$fitted.values
	fit.init <- temp.rob$fitted.values
	
	done <- FALSE
	conv <- NULL
	qest <- matrix(0, nrow = ncol(x), ncol = length(q))
	qfit <- matrix(0, nrow = nrow(x), ncol = length(q))
	qres <- matrix(0, nrow = nrow(x), ncol = length(q))
	qvar <- matrix(0, nrow = ncol(x), ncol = length(q))
	for(i in 1:length(q)) {
		
#We define the starting values
		resid <- resid.init
		fit<-fit.init
		a.j<- case.weights
		w<-case.weights
		coef<-temp.rob$coef
		
		for (iiter in 1:maxit) {
			
			resid.old <- resid
			coef.old<-coef
			
			
# We define the probability mu=exp(xb)/(1+exp(xb))
			probab<-fit
			mu <- probab
			
#We define the variance
			V <- probab * (1 - probab)
			
#We define the scale 
			scale<- c(sqrt(V))
			
#We standardize the residuals
			r.stand <- (y-mu)/sqrt(V)
			
#we compute i1 and i2
			jinf <- floor(mu-k*sqrt(V))
			jsup <- floor(mu+k*sqrt(V))
			ni=rep(1,n)
			
#We compute the values of a_j(b)
			if(k==Inf)
			{
				a.j <- numeric(ni)
			}
			if(k!=Inf)
			{
				indic <- ifelse(jinf+1<=1 & jsup>=1,1,0)
				a.j <- -k*pbinom(jinf,1,probab)+
				k*(1-pbinom(pmin(jsup,ni),1,probab))+ 
				1/sqrt(V)*ifelse(ni==1,probab*indic,
								 mu*(pbinom(pmin(jsup-1,ni-1),pmax(0,1),probab)- 
									 pbinom(jinf-1,pmax(0,1),probab))) -
				mu/sqrt(V)*(pbinom(pmin(jsup,1),1,probab) - pbinom(jinf,1,probab))  }
			a.j<-2*a.j*(q[i]*(r.stand>0)+(1-q[i])*(r.stand<=0))
			
#we define a part of w_j
			w<-c((mu*(1-mu))/scale)*c(w.x)
			
#we compute psi_q(res)
			tmp <- psi.huber((resid)/scale,k=k) * case.weights*((resid)/scale)
			tmp1 <- 2 * (1 - q[i]) * tmp
			tmp1[resid > 0] <- 2 * q[i] * tmp[resid > 0]
			tmp <- tmp1
			
#we compute psi_q(r )-E(psi_q(r ))
			A<-(tmp-a.j)
			
#We compute the -E(psi_prime)
			res1<-(1-mu)
			res0<-mu
			tmp <- psi.huber((res1)/scale,k=k) * case.weights*((res1)/scale)
			tmp1 <- 2 * (1 - q[i]) * tmp
			tmp1[resid > 0] <- 2 * q[i] * tmp[resid > 0]
			tmp <- tmp1
			B.tmp<-tmp
			tmp <- psi.huber((res0)/scale,k=k) * case.weights*((res0)/scale)
			tmp1 <- 2 * (1 - q[i]) * tmp
			tmp1[resid > 0] <- 2 * q[i] * tmp[resid > 0]
			tmp <- tmp1
			B<-c(c(V)*(B.tmp+tmp))
			
#We estimate betas
			temp <- coef+solve(t(x*w*B)%*%x)%*%t(x*w)%*%A
			
			coef <- temp
			eta <- x%*%coef
			fit <- exp(eta)/(1+exp(eta))
			resid <- y-fit
			convi <- irls.delta(coef.old, coef)
			conv <- c(conv, convi)
			done <- (convi <= acc)
			if (done)
			break
		}
		if (!done)
		warning(paste("MQlogit failed to converge in", maxit, "steps at q = ", q[i]))
		
		
# Asymptotic estimated variance of the robust estimator
		
		probab<-fit
		mu <- probab
		deriv.mu<-exp(eta)/((1+exp(eta))^2)
		
#We define the variance
		V <- probab * (1 - probab)
		r.stand <- (y-mu)/sqrt(V)
		
		basepsiq<-function(res,q)
		{
			tmp <- psi.huber(res,k=k) *res
			tmp1 <- 2 * (1 - q) * tmp
			tmp1[res > 0] <- 2 * q * tmp[res > 0]
			tmp <- tmp1
		}
		
#p<-ncol(x)
# Asymptotic estimated variance of the robust estimator
		esp.cond <- basepsiq((1-mu)/sqrt(V),q=q[i])*mu+basepsiq(-mu/sqrt(V),q=q[i])*(1-mu)
		a.const <- apply(x*as.vector(1/n/sqrt(V)*esp.cond*deriv.mu), 2,sum)  
		esp.carre.cond <- (basepsiq((1-mu)/sqrt(V),q=q[i]))^2*mu+(basepsiq(-mu/sqrt(V),q=q[i]))^2*(1-mu)
		matQaux <- as.vector(esp.carre.cond/V*w.x^2*(mu*(1-mu))^2)
		matQ1 <- (1/n)*t(x)%*%(matQaux*x)
		matQ2 <- a.const%*%t(a.const)
		matQ <- matQ1-matQ2
		
		esp.psi.score <- basepsiq((1-mu)/sqrt(V),q=q[i])-basepsiq(-mu/sqrt(V),q=q[i])
		matMaux <- as.vector(esp.psi.score/sqrt(V)*w.x*V*(mu*(1-mu)))
		matM <- 1/n*t(x)%*%(matMaux*x)
		matMinv <- solve(matM)
		
		as.var <- 1/n*matMinv%*%matQ%*%matMinv
		
		
		qest[, i] <- coef
		qfit[, i] <- fit
		qres[,i] <- y-fit
		qvar[,i]<-as.numeric(round(diag(as.var),4))
	}
	list(fitted.values=qfit, var.beta=qvar,residuals=qres, q.values=q, coefficients=qest, matQ=matQ, matM=matM)
}