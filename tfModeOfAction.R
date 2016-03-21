### collapse duplicate columns by mean

col_duplicates<-function(expmat){
dup <- colnames(expmat)[which(duplicated(colnames(expmat)))]
nodup <- setdiff(colnames(expmat),dup)
newnames <- unique(colnames(expmat))
mat <- matrix(ncol=length(newnames),nrow=nrow(expmat))
rownames(mat) <- rownames(expmat)
colnames(mat) <- newnames
mat[,nodup]<-expmat[,nodup]
for(j in dup){
	mat[,j]<-apply(expmat[,which(colnames(expmat)==j)],1,mean)
	}
return(mat)
}
  

###


peak.finder<-function(restab, span=0.2,cut=0.5){

len.tail<-ceiling(span*ncol(restab))
cerotail<-rep(0,len.tail)
getpeaks <- function(i,cut){
	p<-which(ppc.peaks(c(cerotail,restab[i,],cerotail),span)) - len.tail
	maxp<-restab[i,p]
	names(maxp)<-p
	p<-which(ppc.peaks(c(cerotail,-restab[i,],cerotail),span)) - len.tail
	minp<-restab[i,p]
	names(minp)<-p
	c(maxp[which(maxp > cut)],minp[which(minp < -cut)])
	}
loc<-sapply(rownames(restab), function(i) getpeaks(i,cut))
hits<-loc[sapply(names(loc), function(i) length(loc[[i]])) > 0]
max.peaks<-max(sapply(names(loc), function(i) length(loc[[i]])))
peak.tab<-c()
for(j in 1:max.peaks){
	peak.index <- as.numeric(sapply(names(hits), function(i) names(hits[[i]][j])))
	peak.names <- sapply(names(hits), function (i) paste(i,"_",j,sep=""))
	peak.corr <- as.numeric(sapply(names(hits), function(i) hits[[i]][j]))
	tab <- cbind(peak.index,peak.corr)
	rownames(tab) <- peak.names
	peak.tab <- rbind(peak.tab,na.omit(tab))
	}
res<-list(corr=restab,peak=peak.tab)
return(res)
}


namesplit<-function(l){
a<-sapply(l,function(i) strsplit(i,"_")[[1]][1])
return(unique(a))
}

nesguorquer<-function(tf,mexp,k){

mexp<-mexp[,names(sort(mexp[tf,]))]
l<-ncol(mexp)-k
genelist<-setdiff(rownames(mexp),tf)
restab<-matrix(ncol=l,nrow=length(genelist))
rownames(restab)<-genelist
for(j in 1:l){
	n<-j+k
	print (paste(j,n))
	restab[,j] <- cor(mexp[tf,j:n],t(mexp[genelist,j:n]))
	}
return(restab)
}


win.corr<-function(mexp,gene=NULL,lead=NULL,window=0.2,method="pearson"){

if(is.null(lead) && is.null(gene)){stop("you must provide a gene name or a vector!")}

if(is.null(lead)){
	if(gene %in% rownames(mexp)){
		genelist<-setdiff(rownames(mexp),gene)
		mexp<-mexp[,names(sort(mexp[gene,]))]
		lead<-sort(mexp[gene,])
		}
	else{	stop(paste("Gene name ",gene," does not exists!"))}
	}
else{	if(is.null(gene)){genelist<-rownames(mexp)}
	else{genelist<-setdiff(rownames(mexp),gene)}
	if(length(intersect(colnames(mexp),names(lead))) == 0){stop("Vector provided do not overlap with expression data")}
	mexp<-mexp[,names(sort(lead[intersect(colnames(mexp),names(lead))]))]
	lead<-sort(lead[intersect(colnames(mexp),names(lead))])
	}
pb <- txtProgressBar(style=3)
k<-floor(ncol(mexp)*window)
l<-ncol(mexp)-k
restab<-matrix(ncol=l,nrow=length(genelist))
rownames(restab)<-genelist
for(j in 1:l){
	restab[,j] <- cor(lead[j:(j+k)],t(mexp[genelist,j:(j+k)]), method = method)
	setTxtProgressBar(pb, j/l)
	}
close(pb)
return(restab)
}


nesploter.fit<-function(gene,restab,cut){
index<-c(1:ncol(restab))
plot(x=NULL,y=NULL,xlim=range(c(1,ncol(restab))),ylim=range(c(min(restab),max(restab))))
maxres <- which(apply(restab,1,max)  > cut)
minres <- which(apply(restab,1,min)  < -cut)
for(i in setdiff(names(maxres),gene)){
	if(max(abs(restab[i,])) == max(restab[i,]) && max(abs(restab[i,])) > cut){
		y.loess <- loess(restab[i,]~ index, span=0.25)
		y.predict <- predict(y.loess)
		lines(index, y.predict, lwd=0.3,col="red")
		print (i)
		}
	}
for(i in setdiff(names(minres),gene)){
	if(max(abs(restab[i,])) == -min(restab[i,]) && max(abs(restab[i,])) > cut){
		y.loess <- loess(restab[i,]~ index, span=0.25)
		y.predict <- predict(y.loess)
		lines(index, y.predict, lwd=0.3,col="blue")
		print (i)
		}
	}

}

nesploter<-function(filename,restab,cutoff){
index<-c(1:ncol(restab))
plot(x=NULL,y=NULL,xlim=range(c(1,ncol(restab))),ylim=range(c(min(restab),max(restab))))
maxres <- which(apply(restab,1,max)  > cut)
minres <- which(apply(restab,1,min)  < -cut)
for(i in names(maxres)){
	if(max(abs(restab[i,])) == max(restab[i,]) && max(abs(restab[i,])) > cut){
		lines(index, restab[i,], lwd=0.3,col="red")
		print (i)
		}
	}
for(i in names(minres)){
	if(max(abs(restab[i,])) == -min(restab[i,]) && max(abs(restab[i,])) > cut){
		lines(index, restab[i,], lwd=0.3,col="blue")
		print (i)
		}
	}

}


geneVSnet<-function(geneset,net,genelist){
res<-matrix(ncol=7,nrow=length(names(net)))
rownames(res)<-names(net)
s1 <- intersect(geneset,genelist)
for(i in names(net)){
	s2 <- intersect(net[[i]],genelist)
	common <- intersect(s1,s2)
	ft<-fisher.test(matrix(c(length(common),(length(s1)-length(common)),(length(s2)-length(common)),(length(genelist)-length(s1)-length(s1)+length(common))),2,2))
	res[i,5] <-ft$p.value
	res[i,4] <-ft$estimate
	res[i,1] <-length(common)
	res[i,2] <-length(s1)
	res[i,3] <-length(s2)
	res[i,7] <- paste(common,collapse=" ",sep="")
	}
res[,6] <- p.adjust(res[, 5],method="fdr")
res<-res[order(res[,5]),]
colnames(res)<-c("overlap","size.test","size.case","estimate","FET.P","fdr","Genes")
return(res)
}

geneVSnet.Hyp<-function(geneset,net,genelist){
res<-matrix(ncol=5,nrow=length(names(net)))
rownames(res)<-names(net)
for(i in names(net)){
	s1 <- intersect(geneset,genelist)
	n_A<-length(s1)
	s2 <- intersect(names(net[[i]]),genelist)
	n_B<-length(s2)
	n_A_B <- length(intersect(s1,s2))
	n_C <- length(genelist)
	p<-phype(n_A_B - 1, n_A, n_C-n_A, n_B, lower.tail = FALSE)
	res[i,4] <- p
	res[i,1] <- n_A_B
	res[i,2] <- n_A
	res[i,3] <- n_B
	}
res[,5] <- p.adjust(res[, 4],method="fdr")
res<-res[order(res[,4]),]
colnames(res)<-c("overlap","size.test","size.case","p.value","fdr")
return(res)
}


importAracne <- function(filename, type="adj") {
   if (type=="adj") {
     a <- readLines(filename)
     a <- a[-grep(">", a)]
     a <- strsplit(a, "\t")
     target <- lapply(a, function(x) list(x[1], matrix(x[-1], length(x[-1])/2, 2, byrow=T)))
     nom <- sapply(target, function(x) x[[1]])
     target <- lapply(target, function(x) x[[2]])
     target <- cbind(rep(nom, sapply(target, nrow)), unlist(lapply(target, function(x) x[, 1]), use.names=F), unlist(lapply(target, function(x) x[, 2]), use.names=F))
     return(split(target[, 2], target[, 1]))
   }
   if (type=="3col") {
     target <- t(sapply(strsplit(readLines(filename), "\t"), function(x) x[1:2]))
     return(split(target[, 2], target[, 1]))
   }
   stop("Only adj and 3col are acepted as types", call.=F)
 }


TFmodePos <- function(regulon, direction = 1){
newnet<-list()
for(i in names(regulon)){
	newnet[[i]] <- c(rep(direction, length(regulon[[i]])))
	names(newnet[[i]]) <- regulon[[i]]
	}
return(newnet)
}

 
 TFmode <- function(regulon, expset, method="spearman") {
   regulon <- regulon[names(regulon) %in% rownames(expset)]
   regulon <- lapply(regulon, function(x) {if (is.null(names(x))) return(x); names(x)})
   regulon <- lapply(regulon, function(x, genes) x[x %in% genes], genes=rownames(expset))
   tf <- rep(names(regulon), sapply(regulon, length))
   target <- unlist(regulon, use.names=F)
   score1 <- cor(t(expset[rownames(expset) %in% unique(tf), ]), t(expset[rownames(expset) %in% unique(target), ]), method=method)
   res <- lapply(1:length(regulon), function(i, regulon, score) score[rownames(score)==names(regulon)[i], colnames(score) %in% regulon[[i]]], 
 regulon=regulon, score=score1)
   names(res) <- names(regulon)
   return(res)
 }
#### color vector by its elements

vcolor <- function(x,cols){
if(length(unique(x)) != length(cols)) stop ("Wrong number of colors!")
y<-as.character(x)
items<-sort(unique(y))
for(i in 1:length(items)){
	y[which(y == items[i])]<-cols[i]
	}
names(y)<-names(x)
return(y)
}

## operating with networks


gene2net <- function(genelist=NULL,collection=NULL,pathways=NULL){
#### Returns the list of pathgways containing one or several genes (TRUE/FALSE matrix)
# genes => a vector of genes
# collection => a list of genesets
# pathway => a vector of pathway names; if empty then run the whole collection
if(pathways==NULL) pathways<-names(collection)
out<-matrix(ncol=length(genesets),nrow=length(genes))
colnames(out) <- genesets
rownames(out) <- genes
for(g in genes){
	out[g,] <- sapply(genesets, function(i) g %in% net[[i]])
	}
return(out)
}

translate.net<-function(net,ann){
netout<-list()
for(i in names(net)){
	tf <- ann[i]
	if(!is.na(tf)){
		netout[[tf]]<-unique(c(netout[[tf]],unname(ann[net[[i]]])))
		}
	}
return(netout)
}

net.unique<-function(net){
uninet<-list()
for(i in names(net)){
	uninet[[i]]<-unique(c(uninet[[i]],names(net[[i]])))
	}
return(uninet)
}


net.intersect<-function(net1,net2){
net1Anet2<-list()
for(i in intersect(names(net1),names(net2))){
	if(length(intersect(names(net1[[i]]),names(net2[[i]])))>0){
		net1Anet2[[i]] = intersect(names(net1[[i]]),names(net2[[i]]))
		}
	}
return(net1Anet2)
}

net.union<-function(net1,net2){
net1Unet2<-list()
for(i in unique(c(names(net1),names(net2)))){
	net1Unet2[[i]] = unique(c(net1Unet2[[i]],names(net1[[i]]),names(net2[[i]])))
	}
return(net1Unet2)
}

net.diff<-function(net1,net2){
net1Dnet2<-list()
for(i in names(net1)){
	if(length(setdiff(names(net1[[i]]),names(net2[[i]]))) > 0){
		net1Dnet2[[i]] = setdiff(names(net1[[i]]),names(net2[[i]]))
		}
	}
return(net1Dnet2)
}

net2col<-function(net=NULL){
netsize <- 0
for(i in names(net)){
	for(j in names(net[[i]])){
		netsize <- netsize +1
		}
	}
colnet<-matrix(nrow=netsize,ncol=3)
netsize<-0
for(i in names(net)){
	for(j in names(net[[i]])){
		netsize <-netsize +1
		colnet[netsize,1]<-i
		colnet[netsize,2]<-j
		}
	}
return(colnet)
}

### VIPER associated functions


viperTop <- function(vipermat,top=25,side="pos"){
# Obtain top MRs from each sample from a Viper matrix 
# side = c("pos","ned","abs")
if(side == "pos") vm <- - vipermat
if(side == "abs") vm <- abs(vipermat)
if(side == "neg") vm <- vipermat

rankmat<-apply(vm,2,rank, ties.method = "random")
topmat <- sapply(colnames(rankmat), function(i) 
	names(which(rankmat[,i] <= top))  	
	)
return (topmat)
}

viperTopP <- function(vipermat,nes=2,side="pos"){
# Obtain top MRs from each sample from a Viper matrix 
# side = c("pos","ned","abs")
if(side == "pos") vm <- vipermat
if(side == "abs") vm <- abs(vipermat)
if(side == "neg") vm <- -vipermat

toplist <- sapply(colnames(vm), function(i) 
	names(which(vm[,i] > nes))  	
	)
return (toplist)
}

vpFreq<-function(vtop,percent=5){
# agregate top MRs from viperTop
if(is.list(vtop)) {
cutoff <- round(length(vtop)*percent/100)
mrlist<-names(which(table(unlist(vtop)) > cutoff))
}else{
cutoff <- round(ncol(vtop)*percent/100)
mrlist<-names(which(table(c(vtop)) > cutoff))
}
return(mrlist)
}


