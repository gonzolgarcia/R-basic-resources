
######### functions
#########

top.pairs<-function(summ,featureMatrix.snp,numhits=40)
{
## this function returns a table with the most significan feature-TF pairs
# requires (entrez2symbol)
toppair<-head(sort(unlist(summ),decreasing=TRUE),numhits)

outtab <-matrix(ncol=6,nrow=length(toppair))
rownames(outtab) <- names(toppair)
outtab[,6] <- unname(toppair)
colnames(outtab)<-c("featureID","SNP geneid","SNP geneSymbol","TF geneID","TF geneSymbol","#pairs")
for (i in rownames(outtab)){
	snpid <- strsplit(i, "\\.")[[1]][1]
	tfgeneid <- strsplit(i, "\\.")[[1]][2]
	snpGeneID <-featureMatrix.snp[,"NAME"][which(featureMatrix.snp[,"ID"] == snpid)]

	outtab[i,1:5] <- c(snpid,snpGeneID,entrez2symbol(snpGeneID),tfgeneid,entrez2symbol(tfgeneid))
	}
return(outtab)

}


top.pairs.nbl<-function(summ,featureMatrix.snp,numhits=40)
{
## this function returns a table with the most significan feature-TF pairs
# requires (entrez2symbol)
toppair<-head(sort(unlist(summ),decreasing=TRUE),numhits)

outtab <-matrix(ncol=3,nrow=length(toppair))
rownames(outtab) <- names(toppair)
outtab[,3] <- unname(toppair)
colnames(outtab)<-c("featureID","TF geneSymbol","#pairs")
for (i in rownames(outtab)){
	featureID <- strsplit(i, "\\.")[[1]][1]
	tfid <- strsplit(i, "\\.")[[1]][2]
	outtab[i,1:2] <- c(featureID,tfid)
	}
return(outtab)

}

### 

poor.summary <- function(poor.mindy.ouput, pvalcut=1e-5, entrez2symbol=FALSE)
{
	res<-list()
	for(i in names(poor.mindy.ouput)){
		res[[i]]<-sapply(names(poor.mindy.ouput[[i]]["p.value",]), 
				function(j) length(poor.mindy.ouput[[i]]["p.value",j][[1]][which(poor.mindy.ouput[[i]]["p.value",j][[1]] < pvalcut)]))
	}
return(res)
}

# The input will be a featureMatrix in CiTRUS format, a list of TFs and an expression matrix
#

filterFeatures <- function(featureMatrix,expmat,minsize=10){
	# This takes care of the fact that sometimes there are -1
	feature.tmp <- as.matrix(featureMatrix[,intersect(colnames(featureMatrix),colnames(expmat))])
	# Number of non 0 events per genomic locus and type of mutation event
	events <- sapply(rownames(feature.tmp), function(i) sum(abs(na.omit(as.numeric(feature.tmp[i,])))))
	#rowSums(feature.tmp,na.rm=TRUE)
	valid.events1<-events[which(events > minsize)]
	len <- sapply(names(valid.events1),function(i) length(na.omit(feature.tmp[i,])))
	valid.events2<-(len - valid.events1)[which((len - valid.events1) > minsize)]
	valid.final<-rownames(unique(feature.tmp[names(valid.events2),]))
	return(valid.final)
}

poor.mindy <- function(featureMatrix,tflist,expmat,method="pearson",pval=0.0001,diff.exp.pval=0.05,diff.exp.tf.pval=0.05,minsize=20,variancep=1){
	## Given a list of binary features, we study associated to expression samples
# we study changes in the activity of TFs by computing changes in the correlation with 
# any possible target gene
	###############
# featureMatrix: a data.frame in the CITRUS featureMatrix format
# tflist: vector containing list of TF genes
# expmat: an expression matrix (sample names should match those in featureMatrix)
# method: "pearson", "spearman","kendall"
# pval: minimum p.value for the significance of the difference of two correlation to be considered
# diff.exp.pval: pvalue cut off for target genes differentially expressed
# diff.exp.tf.pval: pvalue cut off for TF genes differentially expressed
# minsize: group minimum size from featureMatrix
	
	validFeatures <- filterFeatures(featureMatrix,expmat,minsize=minsize)
	if(length(validFeatures) == 0){ b
		stop("No valid features available on featureMatrix")
	} else {
		message("Total number of valid features: ",length(validFeatures))
	}
	featureMatrix <- featureMatrix[validFeatures,]
	
	message("Iteration over features started at ",Sys.time())
	
	tfs<-intersect(tflist,rownames(expmat))
	if(length(tfs) == 0) stop("0 TFs defined in tflist")
	result <- list()
	for(i in rownames(featureMatrix)){
		featureName <- featureMatrix[i,"NAME"]
		print(paste("calculating",featureName,"feature",i,"out of",nrow(featureMatrix),sep=" "))
		
		# Split the expression matrix in two parts: with (1) and without (0) 
		pheno1 <- intersect(names(featureMatrix[i,][which(featureMatrix[i,] != 0)]),colnames(expmat))
		na <- length(pheno1)
		pheno2 <- intersect(names(featureMatrix[i,][which(featureMatrix[i,] == 0)]),colnames(expmat))
		nb <- length(pheno2)
		# generate lists of non-differentially expressed target genes and tfs
		diff.exp <- myttest(expmat[,pheno1],expmat[,pheno2])
		non.diff.genes <- names(diff.exp$p.value[which(diff.exp$p.value > diff.exp.pval)])
		non.diff.tfs <- intersect(names(diff.exp$p.value[which(diff.exp$p.value > diff.exp.tf.pval)]),tfs)
		if(variancep < 1){
			var.pvals<-myvartest(expmat[non.diff.tfs,pheno1],expmat[non.diff.tfs,pheno2])
			non.diff.tfs<-names(var.pvals[which(var.pvals > variancep)])
			}
		
		# calculate correlation for both non-diff-exp betwen tfs and targets 
		cor.a <- cor(t(expmat[non.diff.tfs,pheno1]),t(expmat[non.diff.genes,pheno1]), method = method) 
		cor.b <- cor(t(expmat[non.diff.tfs,pheno2]),t(expmat[non.diff.genes,pheno2]), method = method)
		zeta <- cor.dif(na,cor.a,nb,cor.b,pval=pval,minsize=minsize)
		tfres<-sapply(rownames(zeta$p.value), function(k) list(
							p.value=zeta$p.value[k,][which(zeta$p.value[k,] < pval)],
							score=zeta$score[k,][which(zeta$p.value[k,] < pval)] )
		)
		result[[i]]<-tfres
	}
	message("Iteration started at ",Sys.time())
	return(result)
}

## this function calculates pvalue for the change on variance for each row between two matrices
myvartest <- function(exp1,exp2){
sapply(intersect(rownames(exp1),rownames(exp2)), function(i) var.test(exp1[i,],exp2[i,])$p.value)
}


cor.dif <- function(na=NULL,ra=NULL,nb=NULL,rb=NULL, pval=0.05, minsize=6){
# Using the Fisher r-to-z transformation, this function calculates z values that can be applied to assess
# the significance of the difference between two correlation coefficients, ra and rb, found in two independent samples.
# If ra is greater than rb, the resulting value of z will have a positive sign; if ra is smaller than rb, the sign of z will be negative.
# na = length of sample set a
# ra = vector of correlations a
# nb = length of sample set b
# rb = vector of correlations b
        #if(na=NULL || ra=NULL || nb=NULL || rb=NULL){stop ("missing data! should provide 1:4 items")}
        if(max(round(ra, 7) ) > 1 || min(round(ra, 7)) < -1){stop("ra must fall between +1.0 and -1.0, inclusive.")}
        if(max(round(rb, 7)) > 1 || min(round(rb, 7)) < -1){stop("rb must fall between +1.0 and -1.0, inclusive.")}
        if(na < minsize || nb < minsize ){stop("n must be equal to or greater than 4.")}
        if(floor(na) < na) {stop("n_a must be an integer value.")}
        if(floor(nb) < nb) {stop("n_b must be an integer value.")}

# first we transform correlation matrix into a matrix of z values
        raplus = 1*ra+1 + 1e-5
        raminus = 1-ra + 1e-5
        rbplus = 1*rb+1 + 1e-5
        rbminus = 1-rb + 1e-5

        za = (log(raplus)-log(raminus))/2
        zb = (log(rbplus)-log(rbminus))/2

        se = sqrt((1/(na-3))+(1/(nb-3)))
        z = (za-zb)/se

#  pvalue estimation R version
        p <- 1-pnorm(abs(z))

        ## pvalue estimation javascript version
#z2 = abs(z)
#p2 =(((((.000005383*z2+.0000488906)*z2+.0000380036)*z2+.0032776263)*z2+.0211410061)*z2+.049867347)*z2+1
#p2 = p2^-16
#p1 = p2/2
#p<-p1[1,]
#z<-z[1,]

        return(list(
        score=z,
        p.value=p))

}


###  
val2gon<-function(z,nbreaks=256,c1="blue",c2="white",c3="red"){
	extreme=max(abs(z))+max(abs(z))/1000
	breaks <- seq(-extreme, extreme, length = nbreaks)
	ncol <- length(breaks)
	col <- colorpanel(ncol,c1,c2,c3)
	CUT <- cut(z, breaks=breaks)
	colorlevels <- col[match(CUT, levels(CUT))] # assign colors to heights for each point
	names(colorlevels)<-names(z)
	return(colorlevels)
}


## t.test function

myttest <-function (x, y, mu = 0, alternative = "two.sided", welch=T)
{
	lx <- ncol(x)
	ly <- ncol(y)
	x.var <- rowVars(x)
	y.var <- rowVars(y)
	if(welch) {
		t <- as.vector((rowMeans(x, na.rm=T) - rowMeans(y, na.rm=T))/sqrt((x.var/lx+y.var/ly)))
		df <- as.vector((x.var/lx+y.var/ly)^2/( (x.var/lx)^2/(lx-1) + (y.var/ly)^2/(ly-1)))
		df[which(df >  (lx+ly-2))] <- lx+ly-2
		df[which(df <  min(lx,ly))] <- min(lx, ly)
	}
	else{
		t <- as.vector((rowMeans(x, na.rm=T) - rowMeans(y, na.rm=T))/
						(sqrt(((lx - 1) * x.var + (ly - 1) * y.var)/(lx + ly - 2))*sqrt(1/lx + 1/ly)))
		df <- lx + ly -2
	}
	names(t) <- rownames(x)
	p <- as.vector(switch(pmatch(alternative,
							c("two.sided", "greater", "less")),
					pt(abs(t), df, lower.tail = F) * 2,
					pt(t, df, lower.tail = F),
					pt(t, df, lower.tail = T)))
	names(p) <- rownames(x)
	list(statistic = t,
			p.value = p)
}

# rowVars calcualtes variance per row
rowVars<-function(x){
	rowSums((x - rowMeans(x, na.rm=T))^2, na.rm=T)/ncol(x)
}
