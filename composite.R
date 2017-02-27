library(gplots)
library(smoother)
library(matrixStats)
#library(pheatmap)

args <- commandArgs(trailingOnly = TRUE)

if (length(args)==3){


data1=read.table(args[1],header=T)
data2=read.table(args[2],header=T)

#dim(data)
data11=as.matrix(as.matrix(data1)[,3:dim(data1)[2]])
data11=apply(data11, 2, as.numeric)
print(dim(data11))

data21=as.matrix(as.matrix(data2)[,3:dim(data2)[2]])
data21=apply(data21, 2, as.numeric)
print(dim(data21))

use_id=as.logical( (apply(data11,1,max)<=1000) * (apply(data21,1,max)<=1000) )
#print(head(use_id))
print(dim(data11[use_id,]))
data11_composite=(colSums(data11[use_id,]))/dim(data21)[1]
data21_composite=(colSums(data21[use_id,]))/dim(data21)[1]

#data11_composite=colQuantiles(data11,probs=0.95)#/dim(data11)[1]
#data21_composite=colQuantiles(data21,probs=0.95)#/dim(data11)[1]

data11_composite=smth(data11_composite,method='gaussian',window=0.06,alpha=2.5,tails=TRUE)
data21_composite=smth(data21_composite,method='gaussian',window=0.06,alpha=2.5,tails=TRUE)

x=seq(-dim(data21)[2]/2,dim(data21)[2]/2-1)

print(length(data11_composite))
print(length(x))
#class(data)
#colors = c(seq(0,quantile(data,0.95)[1],length=25),seq(quantile(data,0.95)[1]+0.1,(quantile(data,0.95)[1]+0.1)*2,length=25))
#colors = c(seq(0,8,length=25),seq(8+0.1,8*2,length=25))
ylim_use=c(0,1.5)
png(paste(args[3],'png',sep='.'))#,res=300,width=500,height=800)
plot(x,data11_composite,type="l",col="blue",xlab='distance from TATA box',ylab='reads counts per gene',lwd=3,ylim=ylim_use,xlim=c(-500,499))
lines(x,data21_composite,col="red",lwd=3)
dev.off()

pdf(paste(args[3],'pdf',sep='.'))#,res=300,width=500,height=800)
plot(x,data11_composite,type="l",col="blue",lwd=3,ylim=ylim_use,xlim=c(-500,499),xlab='distance from TATA box',ylab='reads counts per gene')
lines(x,data21_composite,col="red",lwd=3)
dev.off()

} else {

data1=read.table(args[1],header=T)

#dim(data)
data11=as.matrix(as.matrix(data1)[,3:dim(data1)[2]])
data11=apply(data11, 2, as.numeric)
print(dim(data11))


data11_composite=log2(colSums(data11))#/dim(data11)[1]
#data11_composite=colQuantiles(data11,probs=0.9)#/dim(data11)[1]
#data11_composite=colQuantiles(data11,probs=0.95)#/dim(data11)[1]

print(head(data11_composite))
data11_composite=smth(data11_composite,method='gaussian',window=100,alpha=25,tails=TRUE)

x=seq(-dim(data11)[2]/2,dim(data11)[2]/2-1)

print(length(data11_composite))
print(length(x))
#class(data)
#colors = c(seq(0,quantile(data,0.95)[1],length=25),seq(quantile(data,0.95)[1]+0.1,(quantile(data,0.95)[1]+0.1)*2,length=25))
#colors = c(seq(0,8,length=25),seq(8+0.1,8*2,length=25))
ylim_use=c(5,15)
png(paste(args[2],'png',sep='.'))#,res=300,width=500,height=800)
plot(x,data11_composite,type="l",col="green",xlab='distance from TATA box',ylab='reads counts per gene',lwd=3,ylim=ylim_use,xlim=c(-500,499))
dev.off()

pdf(paste(args[2],'pdf',sep='.'))#,res=300,width=500,height=800)
plot(x,data11_composite,type="l",col="green",lwd=3,ylim=ylim_use,xlim=c(-500,499),xlab='distance from TATA box',ylab='reads counts per gene')
dev.off()


}

