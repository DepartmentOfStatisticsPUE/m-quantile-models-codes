####Exemple for the use of the CD.quant function

source("CD.quant.R")

#GENERATE DATA
areas<-10
n.area<-30
n.area.pop<-300
X.pop<-rnorm(n.area.pop*areas,10,2)
X.pop<-cbind(rep(1,n.area.pop*areas),X.pop)
x<-rnorm(n.area*areas,10,2)
x.design<-cbind(rep(1,n.area*areas),x)
y<-5+2*x+rep(rnorm(areas,1,2),n.area)+rnorm(n.area*areas,0,4)

####VARIABLES DESCRIPTION
#areas <- the number of small area
#n.area <- the sample size of the areas
#y <- the target variable
#x <- auxiliary variable
#X.pop <- auxiliary variable for the population


####SET VARIABLES FOR THE ESTIMATION PROCEDURE (and for optional simulation)
qgrid<-c(0.25,0.50,0.75) #Set of quantiles to be estimated
regioncode<-rep(1:10,each=30)
regioncodepop<-rep(1:10,each=300)

# Quantile estimation Using Chambers Dunstan cdf estimator
cdf.cd.est<-CD.quant(qgrid,y,x.design,X.pop,regioncode,regioncodepop,
								myMSE=TRUE,B=1,R=40,method="eu")

