library(gplots)
library(smoother)
library(matrixStats)
library(MASS)
library(psych)
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

use_id=as.logical( (apply(data11,1,max)<=1000) * (apply(data21,1,max)<=1000) * (apply(data31,1,max)<=1000) * (apply(data41,1,max)<=1000))
#print(head(use_id))
print(dim(data11[use_id,]))
data11_composite=(colSums(data11[use_id,]))/dim(data11)[1]
data21_composite=(colSums(data21[use_id,]))/dim(data21)[1]
data31_composite=(colSums(data31[use_id,]))/dim(data31)[1]
data41_composite=(colSums(data41[use_id,]))/dim(data41)[1]

### generate background model
t=qt(0.95, df=100-1)[1]

data11_composite_bg=matrix(0,100,dim(data11)[2])
data11_composite_bg = t(apply(data11_composite_bg,1, function(x) x+runif(dim(data11)[2])/sum(runif(dim(data11)[2]))*sum(data11[use_id,]) / dim(data11)[1] ))
data11_composite_bg_var=apply(data11_composite_bg, 2, var)
data11_composite_bg_mean=colMeans(data11_composite_bg)
data11_composite_bg_ci_up=data11_composite_bg_mean+data11_composite_bg_var*t/sqrt(1)
data11_composite_bg_ci_down=data11_composite_bg_mean-data11_composite_bg_var*t/sqrt(1)

use_id_2=data11_composite<=data11_composite_bg_ci_up

data11_composite_bg=matrix(0,100,dim(data11)[2])
data11_composite_bg = t(apply(data11_composite_bg,1, function(x) x+runif(dim(data11[use_id,use_id_2])[2])/sum(runif(dim(data11[use_id,use_id_2])[2]))*sum(data11[use_id,use_id_2]) / dim(data11[use_id,use_id_2])[1] ))
data11_composite_bg_var=apply(data11_composite_bg, 2, var)
data11_composite_bg_mean=colMeans(data11_composite_bg)
data11_composite_bg_ci_up=data11_composite_bg_mean+data11_composite_bg_var*t/sqrt(1)
data11_composite_bg_ci_down=data11_composite_bg_mean-data11_composite_bg_var*t/sqrt(1)

data11_composite=smth(data11_composite,method='gaussian',window=0.06,alpha=2.5,tails=TRUE)
data11_composite_bg_mean=smth(data11_composite_bg_mean,method='gaussian',window=0.06,alpha=2.5,tails=TRUE)
data11_composite_bg_ci_up=smth(data11_composite_bg_ci_up,method='gaussian',window=0.06,alpha=2.5,tails=TRUE)
data11_composite_bg_ci_down=smth(data11_composite_bg_ci_down,method='gaussian',window=0.06,alpha=2.5,tails=TRUE)


data21_composite_bg=matrix(0,100,dim(data21)[2])
data21_composite_bg = t(apply(data21_composite_bg,1, function(x) x+runif(dim(data21)[2])/sum(runif(dim(data21)[2]))*sum(data21[use_id,]) / dim(data21)[1] ))
data21_composite_bg_var=apply(data21_composite_bg, 2, var)
data21_composite_bg_mean=colMeans(data21_composite_bg)
data21_composite_bg_ci_up=data21_composite_bg_mean+data21_composite_bg_var*t/sqrt(1)
data21_composite_bg_ci_down=data21_composite_bg_mean-data21_composite_bg_var*t/sqrt(1)

use_id_2=data21_composite<=data21_composite_bg_ci_up

data21_composite_bg=matrix(0,100,dim(data21)[2])
data21_composite_bg = t(apply(data21_composite_bg,1, function(x) x+runif(dim(data21[use_id,use_id_2])[2])/sum(runif(dim(data21[use_id,use_id_2])[2]))*sum(data21[use_id,use_id_2]) / dim(data21[use_id,use_id_2])[1] ))
data21_composite_bg_var=apply(data21_composite_bg, 2, var)
data21_composite_bg_mean=colMeans(data21_composite_bg)
data21_composite_bg_ci_up=data21_composite_bg_mean+data21_composite_bg_var*t/sqrt(1)
data21_composite_bg_ci_down=data21_composite_bg_mean-data21_composite_bg_var*t/sqrt(1)

data21_composite=smth(data21_composite,method='gaussian',window=0.06,alpha=2.5,tails=TRUE)
data21_composite_bg_mean=smth(data21_composite_bg_mean,method='gaussian',window=0.06,alpha=2.5,tails=TRUE)
data21_composite_bg_ci_up=smth(data21_composite_bg_ci_up,method='gaussian',window=0.06,alpha=2.5,tails=TRUE)
data21_composite_bg_ci_down=smth(data21_composite_bg_ci_down,method='gaussian',window=0.06,alpha=2.5,tails=TRUE)

data31_composite_bg=matrix(0,100,dim(data31)[2])
data31_composite_bg = t(apply(data31_composite_bg,1, function(x) x+runif(dim(data31)[2])/sum(runif(dim(data31)[2]))*sum(data31[use_id,]) / dim(data31)[1] ))
data31_composite_bg_var=apply(data31_composite_bg, 2, var)
data31_composite_bg_mean=colMeans(data31_composite_bg)
data31_composite_bg_ci_up=data31_composite_bg_mean+data31_composite_bg_var*t/sqrt(1)
data31_composite_bg_ci_down=data31_composite_bg_mean-data31_composite_bg_var*t/sqrt(1)

use_id_2=data31_composite<=data31_composite_bg_ci_up

data31_composite_bg=matrix(0,100,dim(data31)[2])
data31_composite_bg = t(apply(data31_composite_bg,1, function(x) x+runif(dim(data31[use_id,use_id_2])[2])/sum(runif(dim(data31[use_id,use_id_2])[2]))*sum(data31[use_id,use_id_2]) / dim(data31[use_id,use_id_2])[1] ))
data31_composite_bg_var=apply(data31_composite_bg, 2, var)
data31_composite_bg_mean=colMeans(data31_composite_bg)
data31_composite_bg_ci_up=data31_composite_bg_mean+data31_composite_bg_var*t/sqrt(1)
data31_composite_bg_ci_down=data31_composite_bg_mean-data31_composite_bg_var*t/sqrt(1)

data31_composite=smth(data31_composite,method='gaussian',window=0.06,alpha=2.5,tails=TRUE)
data31_composite_bg_mean=smth(data31_composite_bg_mean,method='gaussian',window=0.06,alpha=2.5,tails=TRUE)
data31_composite_bg_ci_up=smth(data31_composite_bg_ci_up,method='gaussian',window=0.06,alpha=2.5,tails=TRUE)
data31_composite_bg_ci_down=smth(data31_composite_bg_ci_down,method='gaussian',window=0.06,alpha=2.5,tails=TRUE)

data41_composite_bg=matrix(0,100,dim(data41)[2])
data41_composite_bg = t(apply(data41_composite_bg,1, function(x) x+runif(dim(data41)[2])/sum(runif(dim(data41)[2]))*sum(data41[use_id,]) / dim(data41)[1] ))
data41_composite_bg_var=apply(data41_composite_bg, 2, var)
data41_composite_bg_mean=colMeans(data41_composite_bg)
data41_composite_bg_ci_up=data41_composite_bg_mean+data41_composite_bg_var*t/sqrt(1)
data41_composite_bg_ci_down=data41_composite_bg_mean-data41_composite_bg_var*t/sqrt(1)

use_id_2=data41_composite<=data41_composite_bg_ci_up

data41_composite_bg=matrix(0,100,dim(data41)[2])
data41_composite_bg = t(apply(data41_composite_bg,1, function(x) x+runif(dim(data41[use_id,use_id_2])[2])/sum(runif(dim(data41[use_id,use_id_2])[2]))*sum(data41[use_id,use_id_2]) / dim(data41[use_id,use_id_2])[1] ))
data41_composite_bg_var=apply(data41_composite_bg, 2, var)
data41_composite_bg_mean=colMeans(data41_composite_bg)
data41_composite_bg_ci_up=data41_composite_bg_mean+data41_composite_bg_var*t/sqrt(1)
data41_composite_bg_ci_down=data41_composite_bg_mean-data41_composite_bg_var*t/sqrt(1)

data41_composite=smth(data41_composite,method='gaussian',window=0.06,alpha=2.5,tails=TRUE)
data41_composite_bg_mean=smth(data41_composite_bg_mean,method='gaussian',window=0.06,alpha=2.5,tails=TRUE)
data41_composite_bg_ci_up=smth(data41_composite_bg_ci_up,method='gaussian',window=0.06,alpha=2.5,tails=TRUE)
data41_composite_bg_ci_down=smth(data41_composite_bg_ci_down,method='gaussian',window=0.06,alpha=2.5,tails=TRUE)




x=seq(-dim(data21)[2]/2,dim(data21)[2]/2-1)


#class(data)
#colors = c(seq(0,quantile(data,0.95)[1],length=25),seq(quantile(data,0.95)[1]+0.1,(quantile(data,0.95)[1]+0.1)*2,length=25))
#colors = c(seq(0,8,length=25),seq(8+0.1,8*2,length=25))
ylim_use=c(0,1.5)
png(paste(args[5],'procap.png',sep='.'))#,res=300,width=500,height=800)
plot(x,data11_composite,type="l",col="blue",xlab='distance from TATA box',ylab='reads counts per gene',lwd=3,ylim=ylim_use,xlim=c(-500,499))
lines(x,data21_composite,col="red",lwd=3)
#lines(x,data11_composite_bg,col="blue",lwd=1.5,lty=2)
#lines(x,data21_composite_bg,col="red",lwd=1.5,lty=2)
lines(x,data11_composite_bg_mean,col="blue",lwd=1.5,lty=2)
lines(x,data21_composite_bg_mean,col="red",lwd=1.5,lty=2)

#lines(x,data11_composite_bg_ci_up,col="black",lwd=1.5,lty=2)
#lines(x,data11_composite_bg_ci_down,col="black",lwd=1.5,lty=2)
dev.off()

pdf(paste(args[5],'procap.pdf',sep='.'))#,res=300,width=500,height=800)
plot(x,data11_composite,type="l",col="blue",lwd=3,ylim=ylim_use,xlim=c(-500,499),xlab='distance from TATA box',ylab='reads counts per gene')
lines(x,data21_composite,col="red",lwd=3)
lines(x,data11_composite_bg_mean,col="blue",lwd=1.5,lty=2)
lines(x,data21_composite_bg_mean,col="red",lwd=1.5,lty=2)

dev.off()

ylim_use=c(0,11)
png(paste(args[5],'pipseq.png',sep='.'))#,res=300,width=500,height=800)
plot(x,data31_composite,type="l",col="green",xlab='distance from TATA box',ylab='reads counts per gene',lwd=3,ylim=ylim_use,xlim=c(-500,499))
lines(x,data31_composite_bg_mean,col="green",lwd=1.5,lty=2)
dev.off()

pdf(paste(args[5],'pipseq.pdf',sep='.'))#,res=300,width=500,height=800)
plot(x,data31_composite,type="l",col="green",lwd=3,ylim=ylim_use,xlim=c(-500,499),xlab='distance from TATA box',ylab='reads counts per gene')
lines(x,data31_composite_bg_mean,col="green",lwd=1.5,lty=2)
dev.off()

ylim_use=c(0,1.6)
png(paste(args[5],'chipexo.png',sep='.'))#,res=300,width=500,height=800)
plot(x,data41_composite,type="l",col="green",xlab='distance from TATA box',ylab='reads counts per gene',lwd=3,ylim=ylim_use,xlim=c(-500,499))
lines(x,data41_composite_bg_mean,col="green",lwd=1.5,lty=2)
dev.off()

pdf(paste(args[5],'chipexo.pdf',sep='.'))#,res=300,width=500,height=800)
plot(x,data41_composite,type="l",col="green",lwd=3,ylim=ylim_use,xlim=c(-500,499),xlab='distance from TATA box',ylab='reads counts per gene')
lines(x,data41_composite_bg_mean,col="green",lwd=1.5,lty=2)
dev.off()

write.table(data1[-1,][use_id,], file=paste(args[1],'cdt',sep='.'), sep = "\t",row.names=FALSE,col.names=data11_genename,quote=FALSE)
write.table(data2[-1,][use_id,], file=paste(args[2],'cdt',sep='.'), sep = "\t",row.names=FALSE,col.names=data11_genename,quote=FALSE)
write.table(data3[-1,][use_id,], file=paste(args[3],'cdt',sep='.'), sep = "\t",row.names=FALSE,col.names=data11_genename,quote=FALSE)
write.table(data4[-1,][use_id,], file=paste(args[4],'cdt',sep='.'), sep = "\t",row.names=FALSE,col.names=data11_genename,quote=FALSE)

