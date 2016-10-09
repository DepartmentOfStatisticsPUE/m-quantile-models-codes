glm.rob.binom <- function(X,y,weights.on.x=FALSE,chuber=1.345)
{
# Preliminary definitions
	mytol <- .Machine$double.eps^.25
	nb <- length(y)
	assign("nb",nb)
	
	
	basepsi <- function(x)
    {
		x*pmin(1,chuber/abs(x))
    }
	assign("basepsi",basepsi)
	
	basepsiprime <- function(x)
    {
		1*(abs(x)<chuber)
    }
	assign("basepsiprime",basepsiprime)
	
# Initializing....
	beta.old <- as.vector(glm(y~X-1,family=binomial)$coeff)
	Xwith <- X
	
	
	eta <- Xwith%*%beta.old
	probab <- exp(eta)/(1+exp(eta))
	mu <- probab
	V <- probab * (1 - probab)
	deriv.mu <- exp(eta)/((1+exp(eta))^2)
	r.stand <- (y-mu)/sqrt(V)
	ifelse (weights.on.x,w.x <- sqrt(1-hat(X)),w.x <- rep(1,length=nb))
	assign("w.x",w.x)
	assign("Xwith",Xwith)
	assign("y",y)
	
# pmin() on the argument of pbinom() is necessary, due to the 
# bad behavior of pbinom() when evaluated at values greater than size
#
# pmax() is used to avoid the warning messages due to the generation of NA
# in the ifelse procedure
	
	g.objective <- function(beta)
    {
		eta <- Xwith%*%beta
		probab <- exp(eta)/(1+exp(eta))
		mu <- probab
		V <- probab * (1 - probab)
		r.stand <- (y-mu)/sqrt(V)
		deriv.mu <- exp(eta)/((1+exp(eta))^2)
		jinf <- floor(mu-chuber*sqrt(V))
		jsup <- floor(mu+chuber*sqrt(V))
		ni=rep(1,nb)
		if(chuber==Inf)
		{
			esp.cond <- numeric(nb)
		}
		if(chuber!=Inf)
        {
			indic <- ifelse(jinf+1<=1 & jsup>=1,1,0)
			esp.cond <- -chuber*pbinom(jinf,1,probab)+
			chuber*(1-pbinom(pmin(jsup,ni),1,probab))+ 
			1/sqrt(V)*ifelse(ni==1,probab*indic,
							 mu*(pbinom(pmin(jsup-1,ni-1),pmax(0,1),probab)- 
								 pbinom(jinf-1,pmax(0,1),probab))) -
			mu/sqrt(V)*(pbinom(pmin(jsup,1),1,probab) - pbinom(jinf,1,probab))       }
		a.const <- apply(Xwith*as.vector(1/nb/sqrt(V)*w.x*esp.cond*deriv.mu), 2,sum)  
		apply(Xwith*as.vector(1/nb/sqrt(V)*w.x*basepsi(r.stand)*deriv.mu),2,sum)-a.const
    }
	assign("g.objective",g.objective)
	
	grad.g <- function(beta)
    {
		delta <- .Machine$double.eps^.5
		Ident <- diag(1,length(beta))
		1/delta*(apply(beta+delta*Ident,2,g.objective)-as.vector(g.objective(beta)))
    }
	
	tmp.maxit<-0
# Main  
	repeat
    {tmp.maxit<-tmp.maxit+1
		g.old <- g.objective(beta.old)
		grad.g.old <- grad.g(beta.old)
		csi <- solve(grad.g.old,-g.old)
		beta.new <- as.vector(beta.old+csi)
		if(abs(max(beta.old-beta.new))/abs(max(beta.old)) < mytol) break
		if(tmp.maxit>100) break
		beta.old <- beta.new
		NULL
    }
	
    eta <- Xwith%*%beta.old
    fit <- exp(eta)/(1+exp(eta))
	list(coef=beta.old,fitted.values=fit)
}