library(gplots)
library(smoother)
library(matrixStats)
library(MASS)
#library(pheatmap)

args <- commandArgs(trailingOnly = TRUE)


data1=read.table(args[1],header=FALSE)
data2=read.table(args[2],header=FALSE)
data3=read.table(args[3],header=FALSE)

data11=as.matrix(as.matrix(data1)[-1,3:dim(data1)[2]])
data11=apply(data11, 2, as.numeric)
print(dim(data11))

data21=as.matrix(as.matrix(data2)[-1,3:dim(data2)[2]])
data21=apply(data21, 2, as.numeric)
print(dim(data21))

data31=as.matrix(as.matrix(data3)[-1,3:dim(data3)[2]])
data31=apply(data31, 2, as.numeric)
print(dim(data31))


use_id1=as.logical( (apply(data11,1,max)<=1000) )
use_id2=as.logical( (apply(data21,1,max)<=1000) )
use_id3=as.logical( (apply(data31,1,max)<=1000) )
#print(head(use_id))
print(dim(data11[use_id1,]))
data11_used=data11[use_id1,]
data21_used=data21[use_id2,]
data31_used=data31[use_id3,]

data11_composite=(colSums(data11_used))/dim(data11)[1]
data21_composite=(colSums(data21_used))/dim(data21)[1]
data31_composite=(colSums(data31_used))/dim(data31)[1]


ci_up <- function(t) {
  n <- length(t) # n is the sample size
  se <- sd(t)/sqrt(n); # Find the standard error of the sample
  m <- mean(t); # Find the sample mean
  cv <- qt(0.975,df=n-1) # cv is a critical value for the t distribution. P( t > cv ) = 0.025 = P( t < -cv )
  m+cv*se # Return the 95% confidence interval
}
ci_down <- function(t) {
  n <- length(t) # n is the sample size
  se <- sd(t)/sqrt(n); # Find the standard error of the sample
  m <- mean(t); # Find the sample mean
  #print(se)
  cv <- qt(0.975,df=n-1) # cv is a critical value for the t distribution. P( t > cv ) = 0.025 = P( t < -cv )
  #print(cv)
  m-cv*se # Return the 95% confidence interval
}

data11_composite_ci_up=apply(data11_used,2,ci_up)
#print(data11_composite_ci_up)
data11_composite_ci_down=apply(data11_used,2,ci_down)
data21_composite_ci_up=apply(data21_used,2,ci_up)
data21_composite_ci_down=apply(data21_used,2,ci_down)
data31_composite_ci_up=apply(data31_used,2,ci_up)
data31_composite_ci_down=apply(data31_used,2,ci_down)

#data11_composite=colQuantiles(data11,probs=0.95)#/dim(data11)[1]
#data21_composite=colQuantiles(data21,probs=0.95)#/dim(data11)[1]

#data11_composite=smth(data11_composite,method='gaussian',window=0.06,alpha=2.5,tails=TRUE)
#data21_composite=smth(data21_composite,method='gaussian',window=0.06,alpha=2.5,tails=TRUE)
#data31_composite=smth(data31_composite,method='gaussian',window=0.06,alpha=2.5,tails=TRUE)
#data41_composite=smth(data41_composite,method='gaussian',window=0.06,alpha=2.5,tails=TRUE)

x1=seq(-dim(data11)[2]/2,dim(data11)[2]/2-1)
x2=seq(-dim(data21)[2]/2,dim(data21)[2]/2-1)
x3=seq(-dim(data31)[2]/2,dim(data31)[2]/2-1)


print(length(data11_composite))
print(length(x1))
#class(data)
#colors = c(seq(0,quantile(data,0.95)[1],length=25),seq(quantile(data,0.95)[1]+0.1,(quantile(data,0.95)[1]+0.1)*2,length=25))
#colors = c(seq(0,8,length=25),seq(8+0.1,8*2,length=25))
y_uplim=max(c(max(data11_composite),max(data21_composite),max(data31_composite)))
y_downlim=min(c(min(data11_composite),min(data21_composite),min(data31_composite)))

ylim_use=c(y_downlim,y_uplim)

pdf(paste(args[1],'composite.pdf',sep='.'))#,res=300,width=500,height=800)
plot(x1,data11_composite,type="l",col="black",lwd=3,ylim=ylim_use,xlim=c(-24,23),xlab='distance from TATA box (bp)',ylab='')
lines(x1,data11_composite_ci_up,lwd=1.5,lty=2,col="black")
lines(x1,data11_composite_ci_down,lwd=1.5,lty=2,col="black")
dev.off()


pdf(paste(args[2],'composite.pdf',sep='.'))#,res=300,width=500,height=800)
plot(x2,data21_composite,type="l",col="black",lwd=3,ylim=ylim_use,xlim=c(-24,23),xlab='distance from TATA box (bp)',ylab='')
lines(x2,data21_composite_ci_up,lwd=1.5,lty=2,col="black")
lines(x2,data21_composite_ci_down,lwd=1.5,lty=2,col="black")
dev.off()


pdf(paste(args[3],'composite.pdf',sep='.'))#,res=300,width=500,height=800)
plot(x3,data31_composite,type="l",col="black",lwd=3,ylim=ylim_use,xlim=c(-24,23),xlab='distance from TATA box (bp)',ylab='')
lines(x3,data31_composite_ci_up,lwd=1.5,lty=2,col="black")
lines(x3,data31_composite_ci_down,lwd=1.5,lty=2,col="black")
dev.off()

