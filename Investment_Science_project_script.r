##SCENARIO 1

##install required packages
install.packages("matlib")
install.packages("autoimage")
install.packages("fPortfolio")
install.packages("timeSeries")
install.packages("ggplot2")
install.packages("reshape2")
install.packages("quadprog")
##Upload required packages
library(quantmod) 
library(TTR) 
library(matlib)
library(autoimage)
library(fPortfolio)
library(timeSeries) 
library(ggplot2)
library(reshape2)
library(quadprog)
##Clear all variables
rm(list = ls(all.names = TRUE))
##Create portfolio
#Stock 1
stockname='CPALL.BK'
startdate='2014-01-01'
enddate='2019-12-31'
stock.info<-getSymbols(c(stockname),src = 'yahoo',from=startdate,to=enddate,auto.assign=FALSE)
stock.info=na.omit(stock.info)
stock.adj.close.price<-stock.info[,6]
port.price<-stock.info[,6]
port.return<-dailyReturn(stock.info[,6])
#Stock 2
stockname='TMT.BK'
startdate='2014-01-01'
enddate='2019-12-31'
stock.info<-getSymbols(c(stockname),src = 'yahoo',from=startdate,to=enddate,auto.assign=FALSE)
stock.info=na.omit(stock.info)
port.price<-cbind(port.price,stock.info[,6])
port.return<-cbind(port.return,dailyReturn(stock.info[,6]))
#Stock 3
stockname='BTS.BK'
startdate='2014-01-01'
enddate='2019-12-31'
stock.info<-getSymbols(c(stockname),src = 'yahoo',from=startdate,to=enddate,auto.assign=FALSE)
stock.info=na.omit(stock.info)
port.price<-cbind(port.price,stock.info[,6])
port.return<-cbind(port.return,dailyReturn(stock.info[,6]))
port.annual.return<-cbind(annualReturn(port.price[,1]),annualReturn(port.price[,2]),annualReturn(port.price [,3]))
names(port.annual.return) <- c("CPALL", "TMT", "BTS")
covariance.mat<-252*cov(port.return)
rownames(covariance.mat) = c("CPALL", "TMT", "BTS")
colnames(covariance.mat) = c("CPALL", "TMT", "BTS")
##Define function
eff.frontier <- function (returns, short="no", max.allocation=NULL,
                          risk.premium.up=.5, risk.increment=.005){
  covariance <- cov(returns)
  n <- ncol(covariance)
  Amat <- matrix (1, nrow=n)
  bvec <- 1
  meq <- 1
  if(short=="no"){
    Amat <- cbind(1, diag(n))
    bvec <- c(bvec, rep(0, n))
  }
  if(!is.null(max.allocation)){
    if(max.allocation > 1 | max.allocation <0){
      stop("max.allocation must be greater than 0 and less than 1")
    }
    if(max.allocation * n < 1){
      stop("Need to set max.allocation higher; not enough assets to add to 1")
    }
    Amat <- cbind(Amat, -diag(n))
    bvec <- c(bvec, rep(-max.allocation, n))
  }
  loops <- risk.premium.up / risk.increment + 1
  loop <- 1
  eff <- matrix(nrow=loops, ncol=n+3)
  colnames(eff) <- c(colnames(returns), "Std.Dev", "Exp.Return", "sharpe")
  for (i in seq(from=0, to=risk.premium.up, by=risk.increment)){
    dvec <- colMeans(returns) * i
    sol <- solve.QP(covariance, dvec=dvec, Amat=Amat, bvec=bvec, meq=meq)
    eff[loop,"Std.Dev"] <- sqrt(sum(sol$solution*colSums((covariance*sol$solution))))
    eff[loop,"Exp.Return"] <- as.numeric(sol$solution %*% colMeans(returns))
    eff[loop,"sharpe"] <- eff[loop,"Exp.Return"] / eff[loop,"Std.Dev"]
    eff[loop,1:n] <- sol$solution
    loop <- loop+1
  }
  return(as.data.frame(eff))
}
##Find optimal point
#NO Short, Minimized risk, Return >= 0.02
data1<-as.timeSeries(port.annual.return)
eff <- eff.frontier(data1, short="no", max.allocation= NULL, risk.premium.up=0.5, risk.increment=.005)
eff<-subset(eff,eff$CPALL>=0) 
eff<-subset(eff,eff$TMT>=0) 
eff<-subset(eff,eff$BTS>=0) 
eff <- subset(eff,eff$Exp.Return>=0.02) 
eff.optimal.point <- subset(eff,eff$Std.Dev==min(eff$Std.Dev)) 
eff.optimal.point
##Plot Efficient Frontier
ealred <- "#7D110C"
ealtan <- "#CDC4B6"
eallighttan <- "#F7F6F0"
ealdark <- "#423C30"
ggplot(eff, aes(x=Std.Dev, y=Exp.Return)) + geom_point(alpha=.1, color=ealdark) +
  geom_point(data=eff.optimal.point, aes(x=Std.Dev, y=Exp.Return, label=sharpe),
             color=ealred, size=3) +
  annotate(geom="text", x=eff.optimal.point$Std.Dev,
           y=eff.optimal.point$Exp.Return,
           label=paste("Risk: ",
                       round(eff.optimal.point$Std.Dev*100, digits=3),"%\nReturn: ",
                       round(eff.optimal.point$Exp.Return*100, digits=4),"%\nSharpe: ",
                       round(eff.optimal.point$sharpe, digits=2), sep=""),
           hjust=-0.3, vjust=-0.3) +
  ggtitle("Efficient Frontier") +
  labs(x="Risk (standard deviation of portfolio)", y="Return") +
  theme(panel.background=element_rect(fill=eallighttan),
        text=element_text(color=ealdark),
        plot.title=element_text(size=24, color=ealred))
##One government bond
bond.info <- read.csv(file.choose(), head=TRUE, sep=",")
i<-which(bond.info[,1]=="T-BILL1Y")
bond.yield<-bond.info[i,2]
bond.yield<-bond.yield/100
bond.yield.vec<-matrix(bond.yield,nrow=nrow(eff),ncol=1)
tan.theta<-(eff[,6]-bond.yield)/eff[,5]
j<-which(tan.theta==max(tan.theta))
##Plot Efficient Frontier
ggplot(eff, aes(x=Std.Dev, y=Exp.Return)) + coord_cartesian(xlim = c(0.18,0.35),ylim = c(0.15,0.4)) + geom_point(alpha=.1, color=ealdark) +
  geom_point(data=eff[j,], aes(x=Std.Dev, y=Exp.Return),
             color=ealred, size=3) +
  annotate(geom="text", x=eff.optimal.point$Std.Dev,
           y=eff.optimal.point$Exp.Return,
           label=paste("Risk: ",
                       round(eff$Std.Dev[j]*100, digits=4),"%\nReturn: ",
                       round(eff$Exp.Return[j]*100, digits=4),"%\nSharpe: ",
                       round(eff$sharpe[j], digits=4), sep=""),
           hjust=1, vjust=1) +
  ggtitle("Efficient Frontier") +
  labs(x="Risk (standard deviation of portfolio)", y="Return") +
  theme(panel.background=element_rect(fill=eallighttan),
        text=element_text(color=ealdark),
        plot.title=element_text(size=24, color=ealred))
##Find Weight for Risky and Risk-free Assets
n.point=1000
port.stat.mat.w<-matrix(0,nrow=n.point,ncol=4)
for(i in 1:n.point) {
  w<-i/n.point
  port.stat.mat.w[i,1]<-w
  port.stat.mat.w[i,2]<-1-w  
  port.stat.mat.w[i,3]<-w*bond.yield+(1-w)*eff$Exp.Return[j]  
  port.stat.mat.w[i,4]<-(1-w)*eff$Std.Dev[j]
}
port.frame<-data.frame(Weight1<-port.stat.mat.w[,1],
                       Weight2<-port.stat.mat.w[,2],
                       Mean<-port.stat.mat.w[,3],
                       Sd<-port.stat.mat.w[,4])

##Plot Efficient Frontier
ggplot(port.frame, aes(x=Sd, y=Mean)) + geom_point(alpha=.1, color=ealdark) +
  geom_point(data=eff[j,], aes(x=Std.Dev, y=Exp.Return),
             color=ealred, size=3) +
  annotate(geom="text", x=eff$Std.Dev[j],
           y=eff$Exp.Return[j],
           label=paste("Point M    \n", "Risk: ",
                       round(eff$Std.Dev[j]*100, digits=3),"%\nReturn: ",
                       round(eff$Exp.Return[j]*100, digits=4), sep=""),
           hjust=1.5, vjust=1) +
  ggtitle("Efficient Frontier") +
  labs(x="Risk (standard deviation of portfolio)", y="Return") +
  theme(panel.background=element_rect(fill=eallighttan),
        text=element_text(color=ealdark),
        plot.title=element_text(size=24, color=ealred))
##Answer
wopt = (0.02-eff$Exp.Return[j])/(bond.yield-eff$Exp.Return[j])
Answer = data.frame("CPALL" = (1-wopt)*eff.optimal.point$CPALL,"AOT" = (1wopt)*eff.optimal.point$AOT,"BEM" = (1-wopt)*eff.optimal.point$BEM,"T-BILL1Y" = wopt,"Std.Dev" = (1-wopt)*eff$Std.Dev[j],"Exp.Return" = wopt*bond.yield+(1wopt)*eff$Exp.Return[j])
Answer


##SCENARIO 2

rm(list = ls(all.names = TRUE)) #clear objects

library(quantmod) 
library(TTR) 
library("ggplot2")
library(matlib)

library(autoimage)
library(fPortfolio)
library(timeSeries) 
library(ggplot2) # Used to graph efficient frontier
library(reshape2) # Used to melt the data
library(quadprog) # Needed for solve.QP

### graph efficient frontier
# Start with color scheme
ealred <- "#7D110C"
ealtan <- "#CDC4B6"
eallighttan <- "#F7F6F0"
ealdark <- "#423C30"

#par(no.readonly = TRUE)
####

#######################  stock info ############################

###extract stock info CPALL.BK
stock.info.CPALL.BK <- getSymbols(c('CPALL.BK'),src = 'yahoo',from='2016-02-01',to='2020-01-31',auto.assign=FALSE) 
###extract stock info ADVANC.BK
stock.info.ADVANC.BK <- getSymbols(c('ADVANC.BK'),src = 'yahoo',from='2016-02-01',to='2020-01-31',auto.assign=FALSE) 
###extract stock info BDMS.BK
stock.info.BDMS.BK <- getSymbols(c('BDMS.BK'),src = 'yahoo',from='2016-02-01',to='2020-01-31',auto.assign=FALSE) 
###extract stock info EGCO.BK
stock.info.EGCO.BK <- getSymbols(c('EGCO.BK'),src = 'yahoo',from='2016-02-01',to='2020-01-31',auto.assign=FALSE) 

###delete N/A (note: no N/A in this 4 stock)
stock.info.CPALL.BK = na.omit(stock.info.CPALL.BK)
stock.info.ADVANC.BK = na.omit(stock.info.ADVANC.BK)
stock.info.BDMS.BK = na.omit(stock.info.BDMS.BK)
stock.info.EGCO.BK = na.omit(stock.info.EGCO.BK)

### preview stock info
head(stock.info.CPALL.BK)
head(stock.info.ADVANC.BK)
head(stock.info.BDMS.BK)
head(stock.info.EGCO.BK)

###Extract adjusted close price and daily rate of return
#adjusted close price of each asset i.e. consider change in price, dividends and others as a return of a stock
port.price<-cbind(stock.info.CPALL.BK[,6],
                  stock.info.ADVANC.BK[,6], 
                  stock.info.BDMS.BK[,6], 
                  stock.info.EGCO.BK[,6])
head(port.price)

#dailyReturn of each asset
port.return<-cbind(dailyReturn(stock.info.CPALL.BK[,6]),
                   dailyReturn(stock.info.ADVANC.BK[,6]),
                   dailyReturn(stock.info.BDMS.BK[,6]),
                   dailyReturn(stock.info.EGCO.BK[,6]))
#rename column
names(port.return) <- c("CPALL.BK", "ADVANC.BK", "BDMS.BK", "EGCO.BK"); 
head(port.return)


################## plot price and dailyReturn of stocks #####################

#### for CPALL.BK of portfolio
plot(port.price[,1],main="Price",xlab='date',type='l')
plot(port.return[,1],main='Return',xlab='time',type='l')

#### for ADVANC.BK of portfolio
plot(port.price[,2],main="Price",xlab='date',type='l')
plot(port.return[,2],main='Return',xlab='time',type='l')

#### for BDMS.BK of portfolio
plot(port.price[,3],main="Price",xlab='date',type='l')
plot(port.return[,3],main='Return',xlab='time',type='l')

#### for EGCO.BK of portfolio
plot(port.price[,4],main="Price",xlab='date',type='l')
plot(port.return[,4],main='Return',xlab='time',type='l')

#### Volatility of asset price ####
#(2 years period from 2018-02-01 to 2020-01-30)
par(mar=c(1,1,1,1)) 
port.price.vol = port.price[(nrow(port.price)-488):(nrow(port.price)),] 
vol1=volatility(port.price.vol[,1],n=252,N=252,calc="close")
vol1<-na.omit(vol1)
vol2=volatility(port.price.vol[,2],n=252,N=252,calc="close")
vol2<-na.omit(vol2)
vol3=volatility(port.price.vol[,3],n=252,N=252,calc="close")
vol3<-na.omit(vol3)
vol4=volatility(port.price.vol[,4],n=252,N=252,calc="close")
vol4<-na.omit(vol4)
plot(vol1,ylim=c(0,500),main="Volatility",type="l")
lines(vol2, col="Red")
lines(vol3, col="Green")
lines(vol4, col="Blue")
addLegend(legend.loc = "topleft", legend.names = c("CPALL.BK", "ADVANC.BK", "BDMS.BK", "EGCO.BK"), fill = c("black","red", "green", "blue"),)


##################### prepare data for EF calculation ###########################

#annualized rate of return for each asset in portfolio
port.annual.return <- cbind(annualReturn(port.price[, 1]), 
                            annualReturn(port.price[, 2]), 
                            annualReturn(port.price[, 3]),
                            annualReturn(port.price[, 4]))
#rename column
names(port.annual.return) <- c("CPALL.BK", "ADVANC.BK", "BDMS.BK", "EGCO.BK"); 
head(port.annual.return)

###create a time series from port.annual.return
port.annual.return.ts <- as.timeSeries(port.annual.return);


############## statistics of stocks ##########################

#mean of annualized return of each asset in port
mean.asset<-c(mean(annualReturn(port.price[,1])),
              mean(annualReturn(port.price[,2])),
              mean(annualReturn(port.price[,3])),
              mean(annualReturn(port.price[,4])))
#rename column
names(mean.asset) <- c("CPALL.BK", "ADVANC.BK", "BDMS.BK", "EGCO.BK"); 
mean.asset

#(trading day = 252 days in 1 years) 
#SD = sqrt(var) = volatility
#Annualized Volatility = sqrt(252) * SD(daily rate of return, dailyReturn)
#variance and covariance of Annual return of each asset in port
covariance.mat = cov(port.annual.return)
covariance.mat

#Annualized Volatility of each asset in port
sd.asset = sqrt(diag(covariance.mat))
sd.asset


################ risk-free info ###############
bond.info <- read.csv(file.choose(), head=TRUE, sep=",")
bond.yield<-bond.info[which(bond.info[,1]=="T-BILL1Y"),2] # retrieve bond yield(%) of T-BILL1Y
bond.yield<-bond.yield/100 #convert value to percentage (i.e. Annualized risk-free rate of return)
bond.yield


#### Efficient Frontier function ####

eff.frontier <- function (returns, short="no", max.allocation=NULL,
                          risk.premium.up=.5, risk.increment=.005){
  # return argument should be a m x n matrix with one column per security
  # short argument is whether short-selling is allowed; default is no (short
  # selling prohibited)max.allocation is the maximum % allowed for any one
  # security (reduces concentration) risk.premium.up is the upper limit of the
  # risk premium modeled (see for loop below) and risk.increment is the
  # increment (by) value used in the for loop
  
  covariance <- cov(returns) #covariance matrix for return rate
  print(covariance) 
  n <- ncol(covariance) #number of assets in portfolio
  
  # Create initial Amat and bvec assuming only equality constraint
  # (short-selling is allowed, no allocation constraints)
  Amat <- matrix (1, nrow=n)
  bvec <- 1
  meq <- 1
  
  # Then modify the Amat and bvec if short-selling is prohibited
  if(short=="no"){
    Amat <- cbind(1, diag(n))
    bvec <- c(bvec, rep(0, n))
  }
  
  # And modify Amat and bvec if a max allocation (concentration) is specified
  if(!is.null(max.allocation)){
    if(max.allocation > 1 | max.allocation <0){
      stop("max.allocation must be greater than 0 and less than 1")
    }
    if(max.allocation * n < 1){
      stop("Need to set max.allocation higher; not enough assets to add to 1")
    }
    Amat <- cbind(Amat, -diag(n))
    bvec <- c(bvec, rep(-max.allocation, n))
  }
  
  # Calculate the number of loops
  loops <- risk.premium.up / risk.increment + 1
  loop <- 1
  
  # Initialize a matrix to contain allocation and statistics
  # This is not necessary, but speeds up processing and uses less memory
  eff <- matrix(nrow=loops, ncol=n+3)
  # Now I need to give the matrix column names
  colnames(eff) <- c(colnames(returns), "Std.Dev", "Exp.Return", "sharpe")
  
  # Loop through the quadratic program solver
  for (i in seq(from=0, to=risk.premium.up, by=risk.increment)){
    dvec <- colMeans(returns) * i # This moves the solution along the EF
    sol <- solve.QP(covariance, dvec=dvec, Amat=Amat, bvec=bvec, meq=meq)
    eff[loop,"Std.Dev"] <- sqrt(sum(sol$solution*colSums((covariance*sol$solution)))) #Annualized SD (Annualized volatility) of portfolio = sum of weight vector * sum(covariance matrix * weight vector)
    eff[loop,"Exp.Return"] <- as.numeric(sol$solution %*% colMeans(returns)) #Annualized Mean of portfolio
    eff[loop,"sharpe"] <- eff[loop,"Exp.Return"] / eff[loop,"Std.Dev"] #sharpe ratio where rate of return for risk-free asset = 0
    eff[loop,1:n] <- sol$solution #weight for each assets in port
    loop <- loop+1
  }
  
  return(as.data.frame(eff))
}
###### end of the function


################## use EF and shape ratio to find optimal portolio #######################

# Run the eff.frontier function based on no short and no allocation restrictions
eff.port.annual.return.ts <- eff.frontier(port.annual.return.ts, short="no", max.allocation= NULL, 
                                          risk.premium.up=0.5, risk.increment=.005)

#make sure no short selling is allowed (use only case which weight of asset >= 0)
eff.port.annual.return.ts<-subset(eff.port.annual.return.ts,eff.port.annual.return.ts$CPALL.BK>=0)
eff.port.annual.return.ts<-subset(eff.port.annual.return.ts,eff.port.annual.return.ts$ADVANC.BK>=0)
eff.port.annual.return.ts<-subset(eff.port.annual.return.ts,eff.port.annual.return.ts$BDMS.BK>=0)
eff.port.annual.return.ts<-subset(eff.port.annual.return.ts,eff.port.annual.return.ts$EGCO.BK>=0)
head(eff.port.annual.return.ts)

### find optimal porfolio ###

#find portfolio with min.target.mean and max.target.sd
# min.target.mean = 0.06
# max.target.sd = 0.08
target.eff.port.annual.return.ts = eff.port.annual.return.ts
target.eff.port.annual.return.ts = subset(target.eff.port.annual.return.ts, target.eff.port.annual.return.ts[,6] >= 0.06)
target.eff.port.annual.return.ts = subset(target.eff.port.annual.return.ts, target.eff.port.annual.return.ts[,5] <= 0.08)
target.eff.port.annual.return.ts

#recalculate sharpe ratio as risk-free rate = bond.yield for portfolio with min.target.mean and max.target.sd
tan.theta <- (target.eff.port.annual.return.ts[,6] - bond.yield) / target.eff.port.annual.return.ts[,5];#(expect return - bond.yeild)/Std.Dev
tan.theta

#find portfolio with min.target.mean and max.target.sd and have maximum sharpe ratio i.e. optimal portfolio
target.max.sharpe.index = which(tan.theta == max(tan.theta)) #index of port with max sharpe ratio
target.eff.port.annual.return.ts[target.max.sharpe.index,] #optimal port 
tan.theta[target.max.sharpe.index] #sharpe ratio of optimal port 


######### plot of efficient frontier and point of optimal portfolio ############

ggplot(eff.port.annual.return.ts, aes(x=Std.Dev, y=Exp.Return)) + geom_point(alpha=.1, color=ealdark) + #plot EF
  geom_point(data=eff.port.annual.return.ts[target.max.sharpe.index,], aes(x=Std.Dev, y=Exp.Return),
             color=ealred, size=3) + #plot point of portfolio M
  annotate(geom="text", x=eff.port.annual.return.ts$Std.Dev[target.max.sharpe.index],
           y=eff.port.annual.return.ts$Exp.Return[target.max.sharpe.index],
           label=paste("Optimal Portfolio    \n", "Risk: ",
                       round(eff.port.annual.return.ts$Std.Dev[target.max.sharpe.index]*100, digits=3),"% Return: ",
                       round(eff.port.annual.return.ts$Exp.Return[target.max.sharpe.index]*100, digits=4),"% Sharpe ratio: ",
                       round(tan.theta[target.max.sharpe.index], digits=4),"", sep=""),
           hjust=0, vjust=1.4, size=3.5) + #label for point M
  ggtitle("Efficient Frontier") +
  labs(x="Risk (standard deviation of portfolio)", y="Return") +
  theme(panel.background=element_rect(fill=eallighttan),
        text=element_text(color=ealdark),
        plot.title=element_text(size=24, color=ealred)) 
