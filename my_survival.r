#################################
## A list of functions for survival analysis and plot generation 
## This is an addition to survival package
#
# ...

library(survival)
library(gplots)


kmfunct <- function(dat,filename,colores,labs,coxpval=NULL,legend){
km.exp.plot<-survfit(Surv(dat[,1],dat[,2])~dat[,3])
km.exp<-survdiff(Surv(dat[,1],dat[,2])~dat[,3])
km.exp.pval<-sprintf("%.2e",1-pchisq(km.exp$chisq,length(unique(labs))-1))

png(file=paste(filename,".png",sep=""))
par(mar=c(5,5,1,1))
plot(km.exp.plot,col=colores,lwd=3, ylab="Frequency", xlab="Overall survival (days)",
                cex.axis=1.5,cex.lab=1.7,frame.plot = FALSE)
legend("topright",labs,lty=1,col=colores,cex=1.5,lwd=2,bg=rgb(1,1,1,0.8))
if(is.null(coxpval)){
	legend("bottomleft",paste("Logrank p-value: ",km.exp.pval,sep=""),
		cex=1.5,bg=rgb(1,1,1,0.8),lwd=2)
}else{
	legend("bottomleft",c(paste("Logrank p-value: ",km.exp.pval,sep=""),paste("Cox p-value: ",coxpval,sep="")),
		cex=1.5,bg=rgb(1,1,1,0.8))
	}
dev.off()
}


cox_km_plot <- function(dat,filename,gsize=50,ylabexp="Expression"){
## Designed to plot association between survival and contineus variables such as gene expression
## It will calculate the p-value using the cox proportional hazards model and the logrank test

# dat => data.frame with three columbs; 1- survival time; 2- censoring data; 3- conteneous variable
# filename => name of output png file
# group.size => define sample size per km curve
gs <- gsize
int <- gs/2
dat <- dat[order(dat[,3]),]
newdat <- cbind(dat[1:gs,1:3],rep(1,gs))
colnames(newdat) <- c("Surv","Event","Exp")
for(i in 2:(floor(nrow(dat)/int)-1) ){
	mult <- i-1
	slice <- cbind(dat[(int*mult):(gs-1+int*mult),1:3],rep(i,gs))
	colnames(slice) <- c("Surv","Event","Exp")
	newdat<- rbind(newdat,slice)
	}
newdat<-na.omit(newdat)

km.exp.plot <- survfit(Surv(newdat[,1],newdat[,2])~newdat[,4])
coxpval <- summary(coxph(Surv(dat[,1],dat[,2]) ~ dat[,3]))$waldtest["pvalue"]
cox.pval <- sprintf("%.2e",coxpval)

means <- sapply(1:max(newdat[,4]),function(i) mean(newdat[which(newdat[,4] == i),3]))

colorRampPalette(c("green","black","red"))
colores<- colorRampPalette(c("green","black","red"))(length(unique(newdat[,4])))
bars <- data.frame(seq(int,(max(newdat[,4])*int),int),means)

png(file=paste(filename,".png",sep=""),width=800,height=500)
layout(matrix(c(2,1),ncol=2,nrow=1),widths=c(300,500))
par(mar=c(5,0.2,1,5))
plot(km.exp.plot,col=colores,lwd=3, ylab="Frequency", xlab="Survival time (days)",
        cex.axis=1.5,cex.lab=1.7,yaxt='n')
mtext(side=4,"Frequency",cex=1.5,line=3.3)
axis(side=4,at= NULL,cex.lab=2,las=1,cex.axis=1.5)
legend("bottomright",paste("Wald test (P): ",cox.pval,sep=""),title="Cox model",cex=1.5)
par(mar=c(5,5,1,1))
plot(dat[,3],pch=20,cex=0.8,cex.axis=1.5,cex.lab=2,las=1,cex.axis=1.5,ylab=ylabexp)
for(i in 1:nrow(bars)){lines(c(bars[i,1]-int,bars[i,1]+int),c(bars[i,2],bars[i,2]),col=colores[i],lwd=3)}
dev.off()
}

##################################################
##################################################
## obsolete

continuous_km_plot <- function(dat,filename,breaks=2){
## Designed to plot association between survival and contineus variables such as gene expression
## It will calculate the p-value using the cox proportional hazards model and the logrank test

# 
# dat => data.frame with three columbs; 1- survival time; 2- censoring data; 3- conteneous variable
# filename => name of output png file
# brakes => number of groups to break the data (2-4)

dat<-cbind(dat,rep(1,nrow(dat)))
if(breaks==2){
	qtls <- quantile(imm_es[common])
	dat[which(dat[,3] > qtls[3]),4] <- 2
	colores <-c("blue","red")
	labs <- c("--","++")
}else if(breaks == 3){
	qtls <- quantile(imm_es[common],c(0.33,0.66))
	dat[which(dat[,3]> qtls[1]),4] <- 2
	dat[which(dat[,3]> qtls[2]),4] <- 3
	colores <-c("blue","orange","red")
	labs <- c("--","+-","++")
}else if(breaks == 4){
	qtls <- quantile(imm_es[common])
	dat[which(dat[,3]> qtls[2]),4] <- 2
	dat[which(dat[,3]> qtls[3]),4] <- 3
	dat[which(dat[,3]> qtls[4]),4] <- 4
	colores <-c("blue","green","orange","red")
	labs <- c("--","-","+","++")
	}
km.exp.plot<-survfit(Surv(dat[,1],dat[,2])~dat[,4])

coxpval <- summary(coxph(Surv(dat[,1],dat[,2]) ~ dat[,3]))$logtest["pvalue"]
cox.pval<-sprintf("%.2e",coxpval)
logrank<-survdiff(Surv(dat[,1],dat[,2])~dat[,4])
log.pval<-sprintf("%.2e",1-pchisq(logrank$chisq,3))

png(file=paste(filename,".png",sep=""))
par(mar=c(5,5,1,1))
plot(km.exp.plot,col=colores,lwd=3, ylab="Frequency", xlab="Overall survival (days)",
                cex.axis=1.5,cex.lab=1.7)
legend("topright",labs,lty=rep(1,length(labs)),lwd=3,col=colores,cex=1.5)
legend("bottomright",c(paste("Cox p-value: ",cox.pval,sep=""), paste("logrank p-value: ",log.pval,sep="")),cex=1.5)
dev.off()
}

