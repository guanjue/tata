library(gplots)
library(smoother)
library(matrixStats)
library(MASS)
#library(pheatmap)

args <- commandArgs(trailingOnly = TRUE)


data1=read.table(args[1],header=FALSE)
data2=read.table(args[2],header=FALSE)
data3=read.table(args[3],header=FALSE)
data4=read.table(args[4],header=FALSE)

data11=as.matrix(as.matrix(data1)[-1,3:dim(data1)[2]])
data11_genename=as.vector(t(data1[1,]))
print(head(data11_genename))
data11=apply(data11, 2, as.numeric)
print(dim(data11))

data21=as.matrix(as.matrix(data2)[-1,3:dim(data2)[2]])
data21=apply(data21, 2, as.numeric)
print(dim(data21))

data31=as.matrix(as.matrix(data3)[-1,3:dim(data3)[2]])
data31=apply(data31, 2, as.numeric)
print(dim(data31))

data41=as.matrix(as.matrix(data4)[-1,3:dim(data4)[2]])
data41=apply(data41, 2, as.numeric)
print(dim(data41))

### remove genes that have more than 1000 reads at a single position
use_id=as.logical( (apply(data11,1,max)<=1000) * (apply(data21,1,max)<=1000) * (apply(data31,1,max)<=1000) * (apply(data41,1,max)<=1000))

print(dim(data11[use_id,]))
### get composite plots
data11_composite=(colSums(data11[use_id,]))/dim(data21)[1]
data21_composite=(colSums(data21[use_id,]))/dim(data21)[1]
data31_composite=(colSums(data31[use_id,]))/dim(data31)[1]
data41_composite=(colSums(data41[use_id,]))/dim(data41)[1]

### gaussian smoothing
data11_composite=smth(data11_composite,method='gaussian',window=0.06,alpha=2.5,tails=TRUE)
data21_composite=smth(data21_composite,method='gaussian',window=0.06,alpha=2.5,tails=TRUE)
data31_composite=smth(data31_composite,method='gaussian',window=0.06,alpha=2.5,tails=TRUE)
data41_composite=smth(data41_composite,method='gaussian',window=0.06,alpha=2.5,tails=TRUE)

x=seq(-dim(data21)[2]/2,dim(data21)[2]/2-1)

print(length(data11_composite))
print(length(x))

ylim_use=c(0,1.5)
png(paste(args[5],'procap.png',sep='.'))#,res=300,width=500,height=800)
plot(x,data11_composite,type="l",col="blue",xlab='distance from TATA box',ylab='reads counts per gene',lwd=3,ylim=ylim_use,xlim=c(-500,499))
lines(x,data21_composite,col="red",lwd=3)
dev.off()

pdf(paste(args[5],'procap.pdf',sep='.'))#,res=300,width=500,height=800)
plot(x,data11_composite,type="l",col="blue",lwd=3,ylim=ylim_use,xlim=c(-500,499),xlab='distance from TATA box',ylab='reads counts per gene')
lines(x,data21_composite,col="red",lwd=3)
dev.off()

ylim_use=c(0,11)
png(paste(args[5],'pipseq.png',sep='.'))#,res=300,width=500,height=800)
plot(x,data31_composite,type="l",col="green",xlab='distance from TATA box',ylab='reads counts per gene',lwd=3,ylim=ylim_use,xlim=c(-500,499))
dev.off()

pdf(paste(args[5],'pipseq.pdf',sep='.'))#,res=300,width=500,height=800)
plot(x,data31_composite,type="l",col="green",lwd=3,ylim=ylim_use,xlim=c(-500,499),xlab='distance from TATA box',ylab='reads counts per gene')
dev.off()

ylim_use=c(0,1.6)
png(paste(args[5],'chipexo.png',sep='.'))#,res=300,width=500,height=800)
plot(x,data41_composite,type="l",col="green",xlab='distance from TATA box',ylab='reads counts per gene',lwd=3,ylim=ylim_use,xlim=c(-500,499))
dev.off()

pdf(paste(args[5],'chipexo.pdf',sep='.'))#,res=300,width=500,height=800)
plot(x,data41_composite,type="l",col="green",lwd=3,ylim=ylim_use,xlim=c(-500,499),xlab='distance from TATA box',ylab='reads counts per gene')
dev.off()

### get new cdt files which the genes that have more than 1000 reads at a single position are removed
write.table(data1[-1,][use_id,], file=paste(args[1],'cdt',sep='.'), sep = "\t",row.names=FALSE,col.names=data11_genename,quote=FALSE)
write.table(data2[-1,][use_id,], file=paste(args[2],'cdt',sep='.'), sep = "\t",row.names=FALSE,col.names=data11_genename,quote=FALSE)
write.table(data3[-1,][use_id,], file=paste(args[3],'cdt',sep='.'), sep = "\t",row.names=FALSE,col.names=data11_genename,quote=FALSE)
write.table(data4[-1,][use_id,], file=paste(args[4],'cdt',sep='.'), sep = "\t",row.names=FALSE,col.names=data11_genename,quote=FALSE)

