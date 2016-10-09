\name{BinaryMQ-package}
\alias{BinaryMQ-package}
\alias{BinaryMQ}
\docType{package}
\title{
BinaryMQ
}
\description{
M-quantile regression for Binary Data and Small Area Estimation
}
\details{
\tabular{ll}{
Package: \tab BinaryMQ\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2012-07-14\cr
License: \tab GPL (>= 2)\cr
}
}
\author{
Nicola Salvati, Ray Chambers, Nikos Tzavidis

Maintainer: Nicola Salvati <salvati@ec.unipi.it>
}
\references{
Chambers, Salvati, Tzavidis (2012). M-quantile Regression Models for Binary Data in Small Area Estimation.
}

\keyword{ package }
\seealso{
\code{\link[<robustbase>:<robustbase>-package]{<robustbase>}}
}
\examples{
#This is an example for using Binary MQ in small area estimation
library(MASS)
library(pps)
library(BinaryMQ)

m<-100
Est.MQNaive<-NULL
ni<-rep(10,m)
Ni<-rep(100,m)
N<-sum(Ni)
n<-sum(ni)
ar<-1:m

#Generate a population
h=1
set.seed(h)
X<-NULL
bb=8
for (i in 1:m)
{X<-c(X,runif(Ni[1],-1,ar[i]/bb))
}

gamma0=rnorm(m,0,sqrt(0.25))

beta<-1
eta<-X*beta+rep(gamma0,each=100) 
pd<-exp(eta)/(1+exp(eta))

Y<-rbinom(N,1,pd)
regioncode<-rep(1:m,each=100)
id<-seq(1:N)

pop<-cbind(id,Y,X,regioncode);dim(pop)

# Drawing a sample
s<-(sort(stratsrs(pop[,4],ni)))
data.s=cbind(pop[s,2],pop[s,3],pop[s,4])
data.r=cbind(pop[-s,2],pop[-s,3],pop[-s,4])

#M-quantile approach
Q<-c(0.25,0.30,0.4,seq(from=0.45,to=0.55,by=0.005),0.60,0.65,0.75)
temp<- try(QLogit.Scores(y=data.s[,1],x=cbind(rep(1,n),data.s[,2]),qgrid = Q,k.value=1.345))

scores=(cbind(data.s[,3],temp$qscores))

mq=rep(0,m)
for (i in ar){
mq[i]=try(mean(scores[,2][scores[,1]==ar[i]]))
}

for (i in ar){
x.s.tmp<-data.s[data.s[,3]==ar[i],2]
tmp=(glm.mq.binom(y=data.s[,1], x=cbind(rep(1,n),data.s[,2]),q = c(mq[i]),k=1.345))
tmp.eta<-try(exp(tmp$coefficients[1,1]+tmp$coefficients[2,1]*data.r[,2][data.r[,3]==ar[i]]))

Est.MQNaive[i]=try(1/Ni[i]*(sum(data.s[,1][data.s[,3]==ar[i]])+sum((tmp.eta/(1+tmp.eta)))))

}

}

#results
Est.MQNaive
