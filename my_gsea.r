###

gsea_es <- function(reflist, gs, w = 1, sizelim = 10)
{
  #GSEA Gene set enrichment analysis 
  #  reflist : named vector of reference scores
  #  gs   : gene set
  #  w     : weight
  #  es    : enrichment score
  #
  # Author: Wei Keat Lim (wl2131@columbia.edu)
  # Modified from Matlab to R by Celine Lefebvre (lefebvre@c2b2.columbia.edu)
  #
 
  # combine ranked list and score
  ix <- order(reflist, decreasing=T)
  reflist <- reflist[ix]
  
  es <- 0

  if(!is.null(gs) && length(gs) >= sizelim)
    {
      # check overlap of gene sets and ranked list
      isgs <- rep(0, length(reflist))
      isgs[which(names(reflist) %in% gs)] <- 1
      
      # compute ES
      score_hit <- cumsum((abs(reflist*isgs))^w)
      score_hit <- score_hit/tail(score_hit, 1)
      score_miss <- cumsum(1-isgs)
      score_miss <- score_miss/tail(score_miss, 1)
      es_all <- score_hit - score_miss
      #print(c(max(es_all), min(es_all)))
      es <- max(es_all) + min(es_all) 
      #if(max(es_all) > -min(es_all)) es <- max(es_all) 
      #else es <- min(es_all)
  }  
  return(es)
}


 

gsea <- function (reflist, set, method = "permutation",np = 1000, w = 1, gsea_null = NULL)
{
  #GSEA Gene set enrichment analysis 
  #  reflist : named vector of reference scores
  #  gs   : gene set
  #  w     : weight
  #  es    : enrichment score
  #  method = c("permutation", "pareto")
  # Author: Wei Keat Lim (wl2131@columbia.edu)
  # Modified from Matlab to R by Celine Lefebvre 
  # Modified from  Celine Lefebvre by Federico Giorgi 
    set <- intersect(names(reflist), set)
    ix <- order(reflist, decreasing = T)
    reflist <- reflist[ix]
    es <- 0
    nes <- 0
    p.value <- 1
    inSet <- rep(0, length(reflist))
    inSet[which(names(reflist) %in% set)] <- 1
    hits <- abs(reflist * inSet)
    hits <- hits^w
    score_hit <- cumsum(hits)
    score_hit <- score_hit/score_hit[length(score_hit)]
    score_miss <- cumsum(1 - inSet)
    score_miss <- score_miss/score_miss[length(score_miss)]
    running_score <- score_hit - score_miss
    if (all(is.na(running_score))) {
        running_score <- rep(0, length(running_score))
    }
#    if (abs(max(running_score)) > abs(min(running_score))) {
#        es <- max(running_score)
#    }else {
#        es <- min(running_score)
#    }
	es <- max(running_score)+ min(running_score)
    ledge_indeces <- rep(0, length(running_score))
    if (es < 0) {
        peak <- which(running_score == min(running_score))[1]
        ledge_indeces[peak:length(ledge_indeces)] <- 1
        ledge_indeces <- which(ledge_indeces == 1)
        ledge_names <- names(reflist[ledge_indeces])
    }else {
        peak <- which(running_score == max(running_score))
        ledge_indeces[1:peak] <- 1
        ledge_indeces <- which(ledge_indeces == 1)
        ledge_names <- names(reflist[ledge_indeces])
    }
    if (is.null(gsea_null)) {
        null_es <- null_gsea(set = set, reflist = reflist, np = np,w = w)
    }else {
        if (class(gsea_null) == "gsea_nullist") {
            null_es <- gsea_null[as.character(length(set))][[1]]
        }else {
            null_es <- gsea_null
        }
    }
    if (es < 0) {
        p.value <- sum(null_es <= es)/length(null_es)
    }else {
        p.value <- sum(null_es >= es)/length(null_es)
    }
    if (p.value < 0.05) {
        if (p.value == 0) {
            p.value <- 1/np
        }
        if (method == "pareto") {
            q95 <- as.numeric(quantile(abs(null_es), 0.95))
            fit <- pareto.fit(abs(null_es), threshold = q95)
            newp.value <- ppareto(abs(es), threshold = q95, exponent = fit$exponent,
                lower.tail = FALSE)/20
            if (is.na(newp.value)) {
                newp.value <- p.value
            }
            p.value <- newp.value
        }
    }
    nes <- p2z(p.value) * sign(es)
    gsea.obj <- list(es = es, nes = nes, p.value = p.value, ledge = ledge_names,
        running_score = running_score, set = set, reflist = reflist,
        inSet = inSet,nullDist=null_es)
    class(gsea.obj) <- "gsea"
    return(gsea.obj)
}



gsea2 <- function (reflist, set1,set2, method = c("permutation", "pareto"),
    np = 1000, w = 1, gsea_null = NULL)
{
  #GSEA Gene set enrichment analysis of two complementary gene sets using gsea (Federico Giorgi) function 
  #  reflist : named vector of reference scores
  #  set1   : gene set1
  #  set2   : gene set2
  #
  # Author: Gonzalo Lopez
g1 <- gsea(reflist,set1,method=method,np = np, w = w, gsea_null = gsea_null)
g2 <- gsea(reflist,set2,method=method,np = np, w = w, gsea_null = gsea_null)
ix <- order(reflist, decreasing = T)
reflist <- reflist[ix]

gsea.obj <- list(
es1 = g1$es,es2 = g2$es,
nes1 = g1$nes, nes2 =  g2$nes,
p.value1 = g1$p.value,p.value2 =  g2$p.value,
ledge1 = g1$ledge,ledge2 =  g2$ledge,
running_score1 = g1$running_score,running_score2 =  g2$running_score,
set1 = set1,set2 = set2,
reflist = reflist,
inSet1 = g1$inSet,inSet2 = g2$inSet)
return(gsea.obj)
}


pareto.fit <- function (data, threshold)
{
    return(pareto.fit.ml(data, threshold))
}

ppareto <- function (x, threshold = 1, exponent, lower.tail = TRUE)
{
    if (!lower.tail) {
        f <- function(x) {
            (x/threshold)^(1 - exponent)
        }
    }
    if (lower.tail) {
        f <- function(x) {
            1 - (x/threshold)^(1 - exponent)
        }
    }
    p <- ifelse(x < threshold, NA, f(x))
    return(p)
}

pareto.fit.ml <- function (data, threshold)
{
    data <- data[data >= threshold]
    n <- length(data)
    x <- data/threshold
    alpha <- 1 + n/sum(log(x))
    loglike = pareto.loglike(data, threshold, alpha)
    ks.dist <- ks.dist.fixed.pareto(data, threshold = threshold,
        exponent = alpha)
    fit <- list(type = "pareto", exponent = alpha, xmin = threshold,
        loglike = loglike, ks.dist = ks.dist, samples.over.threshold = n)
    return(fit)
}

pareto.loglike <- function (x, threshold, exponent)
{
    L <- sum(dpareto(x, threshold = threshold, exponent = exponent,
        log = TRUE))
    return(L)
}

dpareto <- function (x, threshold = 1, exponent, log = FALSE)
{
    if (!log) {
        prefactor <- (exponent - 1)/threshold
        f <- function(x) {
            prefactor * (x/threshold)^(-exponent)
        }
    }
    else {
        prefactor.log <- log(exponent - 1) - log(threshold)
        f <- function(x) {
            prefactor.log - exponent * (log(x) - log(threshold))
        }
    }
    d <- ifelse(x < threshold, NA, f(x))
    return(d)
}

ks.dist.fixed.pareto <- function (data, threshold, exponent)
{
    data <- data[data >= threshold]
    d <- suppressWarnings(ks.test(data, ppareto, threshold = threshold,
        exponent = exponent))
    return(as.vector(d$statistic))
}


p2z<- function (p)
{
    qnorm(p/2, lower.tail = F)
}

## generation of a null distriburion
null_gsea<-function (set, reflist, w = 1, np = 1000)
{
    gsea_null <- rep(0, np)
    gsea_null <- sapply(1:np, function(i) {
        inSet <- rep(0, length(reflist))
        inSet[which(names(reflist) %in% set)] <- 1
        null_inSet <- inSet[sample(1:length(inSet))]
        null_hit <- abs(reflist * null_inSet)
        null_hit <- null_hit^w
        null_hit <- cumsum(null_hit)
        null_hit <- null_hit/null_hit[length(null_hit)]
        null_miss <- cumsum(1 - null_inSet)
        null_miss <- null_miss/null_miss[length(null_miss)]
        null_running_score <- null_hit - null_miss
        if (abs(max(null_running_score)) > abs(min(null_running_score))) {
            null_es <- max(null_running_score) 
        }
        else {
            null_es <- min(null_running_score)
        }
        return(null_es)
    })
    class(gsea_null) <- "gsea_null"
    return(gsea_null)
}
## generation of a null distriburion for a matrix with ranking colomns

null_gsea_mat <- function (set, mrank, w = 1, np = 10000)
{
    gsea_null <- rep(0, np)
	refnums <- sample(1:ncol(mrank),np,replace=TRUE)    
	gsea_null <- sapply(1:np, function(i) {
        inSet <- rep(0, nrow(mrank))
        inSet[which(rownames(mrank) %in% set)] <- 1
        null_inSet <- inSet[sample(1:length(inSet))]
        null_hit <- abs(mrank[,refnums[i]] * null_inSet)
        null_hit <- null_hit^w
        null_hit <- cumsum(null_hit)
        null_hit <- null_hit/null_hit[length(null_hit)]
        null_miss <- cumsum(1 - null_inSet)
        null_miss <- null_miss/null_miss[length(null_miss)]
        null_running_score <- null_hit - null_miss
        if (abs(max(null_running_score)) > abs(min(null_running_score))) {
            null_es <- max(null_running_score) 
        }
        else {
            null_es <- min(null_running_score)
        }
        return(null_es)
    })
    class(gsea_null) <- "gsea_null"
    return(gsea_null)
}


# function to plot gsea object from gsea function

plot_gsea <- function (gsea.obj, twoColors = c("red", "blue"), plotNames = FALSE,
    title = "Running Enrichment Score", correctEntrez = FALSE,
    bottomYlabel = "Signature values")
{
    es <- gsea.obj$es
    nes <- gsea.obj$nes
    p.value <- gsea.obj$p.value
    ledge <- gsea.obj$ledge
    running_score <- gsea.obj$running_score
    set <- gsea.obj$set
    reflist <- gsea.obj$reflist
    inSet <- gsea.obj$inSet
    if (correctEntrez) {
        set <- entrez2symbol(set)
        names(reflist) <- entrez2symbol(names(reflist))
    }
    min.RES <- min(running_score)
    max.RES <- max(running_score)
    delta <- (max.RES - min.RES) * 0.5
    min.plot <- min.RES
    max.plot <- max.RES
    max.corr <- max(reflist)
    min.corr <- min(reflist)
    Obs.correl.vector.norm <- (reflist - min.corr)/(max.corr -
        min.corr) * 1.25 * delta + min.plot
    zero.corr.line <- (-min.corr/(max.corr - min.corr)) * 1.25 *
        delta + min.plot
    if (es < 0) {
        l.ledge.ref.plot <- length(reflist) - length(ledge)
    }
    else {
        l.ledge.ref.plot <- length(ledge)
    }
    if (nes > 0) {
        col.f <- twoColors[1]
    }
    else {
        col.f <- twoColors[2]
    }
    N <- length(reflist)
    ind <- 1:N
	
	setProp<-length(set)
	if(setProp < 100){	vlinelwd =2
	} else if(setProp < 300){	vlinelwd =1
	} else if(setProp < 600 ){vlinelwd =0.6
	} else if(setProp < 1000){vlinelwd =0.3
	} else{vlinelwd = 0.1}

    layoutMatrix <- rbind(1, 2, 3)
    layout(layoutMatrix, heights = c(1, 4, 1.5))
    par(mar = c(0, 4.1, 2, 2.1))
    plot(0, col = "white", xlim = c(1, N), ylim = c(0, 10), xaxt = "n",
        yaxt = "n", type = "n", frame.plot = FALSE, xlab = "",
        ylab = "", xaxs = "r", yaxs = "r", main = paste("Number of elements: ",
            N, " (in full list), ", length(set), " (in element set)",
            sep = "", collapse = ""))
    for (position in 1:N) {
        if (inSet[position] == 1) {
            if (N < 50 | length(set) <= 10) {
                rect(xleft = position - 0.2, ybottom = 0, xright = position +
                  0.2, ytop = 10, col = col.f, border = NA)
            }
            else {
                abline(v = position, lwd = vlinelwd, col = col.f)
            }
            if (plotNames) {
                text(labels = names(reflist[position]), x = position -
                  0.2, y = 0, srt = 90, offset = 0, pos = 4,
                  font = 2)
            }
        }
    }
    par(mar = c(0, 4.1, 2, 2.1))
    plot(ind, running_score, sub = "", xlab = "", ylab = "Running Enrichment Score",
        xlim = c(1, N), ylim = c(min.plot, max.plot), type = "l",
        pch = 20, lwd = 2, cex = 1, col = col.f, xaxt = "n",las=1,
        yaxs = "r", main = NULL)
    grid(col = "light grey", lty = 2)
    lines(c(1, N), c(0, 0), lwd = 1, lty = 1, cex = 1)
    lines(c(l.ledge.ref.plot, l.ledge.ref.plot), c(min.plot,
        max.plot), lwd = 1, lty = 1, cex = 1)
    if (es >= 0) {
        legend_position <- "topright"
    }
    else {
        legend_position <- "bottomleft"
    }
    legend(legend_position, legend = c(paste("ES = ", signif(es,
        3), sep = ""), paste("NES = ", signif(nes, 3), sep = ""),
        paste("p-val = ", sprintf("%.2e",p.value), sep = "")), bg = "white",cex=1.2)
    par(mar = c(2, 4.1, 2, 2.1))
    plot(x=NULL,y=NULL,xlim=range(c(0,max(ind))),ylim=range(c(min(reflist),max(reflist))), type = "b", lwd = 2,
	cex = 1, col = 1, xaxs = "r", yaxs = "r", main = NULL,las=1,
        ylab = bottomYlabel)
    #grid(col = "light grey", lty = 2)
    abline(h=0, lwd = 1,lty = 2, cex = 1, col = 1)
	lines(ind,reflist,lty=1,lwd=2)
	#legend("top", c("Reference list"),bty='n',cex=1.4)
	}

	
### function to plot gsea object from gsea2 function

plot_gsea2 <- function (gsea.obj, twoColors = c("red", "blue"), plotNames = FALSE,
    title = "Running Enrichment Score", bottomYlabel = "Signature values")
{
    es1 <- gsea.obj$es1; es2 <- gsea.obj$es2;
    nes1 <- gsea.obj$nes1;nes2 <- gsea.obj$nes2;
    p.value1 <- gsea.obj$p.value1; p.value2 <- gsea.obj$p.value2;
    ledge1 <- gsea.obj$ledge1; ledge2 <- gsea.obj$ledge2;
    running_score1 <- gsea.obj$running_score1;running_score2 <- gsea.obj$running_score2;
    set1 <- gsea.obj$set1;set2 <- gsea.obj$set2;
    inSet1 <- gsea.obj$inSet1; inSet2 <- gsea.obj$inSet2;
	reflist <- gsea.obj$reflist
	
    min.plot <- min(running_score1,running_score2)
    max.plot <- max(running_score1,running_score2)

    if (es1 < 0) {
        l.ledge.ref.plot1 <- length(reflist) - length(ledge1)
    }else {
        l.ledge.ref.plot1 <- length(ledge1)
    }

    if (es2 < 0) {
        l.ledge.ref.plot2 <- length(reflist) - length(ledge2)
    }else {
        l.ledge.ref.plot2 <- length(ledge2)
    }

    N <- length(reflist)
    ind <- 1:N

    layoutMatrix <- rbind(1, 2, 3, 4)
    layout(layoutMatrix, heights = c(1, 5, 1, 2))
	##-----------------------------------------
    par(mar = c(0, 4.1, 2, 2.1))
    plot(0, col = "white", xlim = c(1, N), ylim = c(0, 10), xaxt = "n",
        yaxt = "n", type = "n", frame.plot = FALSE, xlab = "",
        ylab = "", xaxs = "r", yaxs = "r", main = paste("Number of elements: ",
            N, " (in full list), ", length(set1), " (in element set1)",
            sep = "", collapse = ""))
    for (position in 1:N) {
        if (inSet1[position] == 1) {
            if (N < 50 | length(set1) <= 10) {
                rect(xleft = position - 0.2, ybottom = 0, xright = position +
                  0.2, ytop = 10, col = twoColors[1], border = NA)
            }
            else {
                abline(v = position, lwd = 1, col = twoColors[1])
            }
            if (plotNames) {
                text(labels = names(reflist[position]), x = position -
                  0.2, y = 0, srt = 90, offset = 0, pos = 4,
                  font = 2)
            }
        }
    }
	##--------------------------------
    par(mar = c(0, 4.1, 2, 2.1))
	plot(x=NULL,y=NULL,xlim=range(c(1,N)), ylim=range(c(min.plot, max.plot)), sub = "", xlab = "", ylab = "Running Enrichment Score",
        pch = 20, lwd = 2, cex = 1, xaxt = "n",
        yaxs = "r", main = NULL)
	    grid(col = "light grey", lty = 2)
    lines(ind, running_score1, sub = "", xlab = "", ylab = "Enrichment Score",
        lwd = 2, cex = 1, col = twoColors[1])
    lines(ind, running_score2, sub = "", xlab = "", ylab = "Enrichment Score",
        lwd = 2, cex = 1, col = twoColors[2])
    lines(c(1, N), c(0, 0), lwd = 1, lty = 1, cex = 1)
    lines(c(l.ledge.ref.plot1, l.ledge.ref.plot1), c(min(running_score1),
        max(running_score1)), lwd = 1, lty = 1, cex = 1)
    lines(c(l.ledge.ref.plot2, l.ledge.ref.plot2), c(min(running_score2),   
        max(running_score2)), lwd = 1, lty = 1, cex = 1)
	if(es1 >=0){legside1="topright"
	}else{legside1="topleft"}
    legend(legside1, legend = c(paste("NES1 = ", signif(nes1, 3), sep = ""),
        paste("P1 = ", sprintf("%.2e",p.value1), sep = "")), 
		bg = "white",box.col=twoColors[1],cex=1.5)
	if(es2 >=0){legside2="bottomright"
	}else{legside2="bottomleft"}
    legend(legside2, legend = c(paste("NES2 = ", signif(nes2, 3), sep = ""),
        paste("P2 = ", sprintf("%.2e",p.value2), sep = "")), 
		bg = "white",box.col=twoColors[2],cex=1.5)

	##----------------------------------
	par(mar = c(0, 4.1, 2, 2.1))
    plot(0, col = "white", xlim = c(1, N), ylim = c(0, 10), xaxt = "n",
        yaxt = "n", type = "n", frame.plot = FALSE, xlab = "",
        ylab = "", xaxs = "r", yaxs = "r", main = paste("Number of elements: ",
            N, " (in full list), ", length(set2), " (in element set2)",
            sep = "", collapse = ""))
    for (position in 1:N) {
        if (inSet2[position] == 1) {
            if (N < 50 | length(set2) <= 10) {
                rect(xleft = position - 0.2, ybottom = 0, xright = position +
                  0.2, ytop = 10, col = twoColors[2], border = NA)
            }
            else {
                abline(v = position, lwd = 1, col = twoColors[2])
            }
            if (plotNames) {
                text(labels = names(reflist[position]), x = position -
                  0.2, y = 0, srt = 90, offset = 0, pos = 4,
                  font = 2)
            }
        }
    }
	#-----------------------------------
    par(mar = c(2, 4.1, 2, 2.1))
    plot(x=NULL,y=NULL,xlim=range(c(0,N)),ylim=range(c(min(reflist),max(reflist))), type = "b", lwd = 2,
	cex = 1, col = 1, xaxs = "r", yaxs = "r", main = NULL, ylab = bottomYlabel)
    grid(col = "light grey", lty = 2)
    abline(h=0, lwd = 1,lty = 1, cex = 1, col = 1)
	lines(ind,reflist,lty=1,lwd=2)
	legend("top", c("Reference list"),bty='n',cex=1.4)
	}

######## PLOT multiple gseas of one set and multiple signatures

multi_sig_gsea <- function(signature,geneset){
# uses signature from viperSignature
# geneset
nsig <-null_gsea(geneset,signature$signature[,i])
out <- list()
for(i in colnames(signature$signature)){
gout <- gsea(signature$signature[,i],geneset,gsea_null=nsig,method="pareto")
out[[i]][["p_value"]] <- gout$p.value
out[[i]][["nes"]] <- gout$nes
}
return(out)
}



	
#############
## Gsea iterative functions

## iterative GSEA on a collection of pathways 

multi.gsea.gon <- function(collection,reflist,minsize=5, w = 1, np = 1000)
{
set.sizes<-sapply(names(collection), function(i) length(intersect(collection[[i]],names(reflist))))
non.empty.sets<-names(set.sizes[which(set.sizes > minsize)])

gsea.res<-data.frame(matrix(ncol=5,nrow=length(non.empty.sets)))
rownames(gsea.res)<-non.empty.sets
colnames(gsea.res)<-c("NES","P.value","fdr","LE.size","gs_size")

## we generate a list of null distr, starting by an NA fillecd list with top in the maximum possible length
toplen <- max(unique(unlist(lapply(collection,length))))
my_null <- lapply(1:toplen,function(i) NA)

pb <- txtProgressBar(style=3)
cc <-0
tot <- length(non.empty.sets)
for(i in non.empty.sets){
	cc <- cc+1
	n <- length(intersect(collection[[i]],names(reflist)))
#	if(is.na(my_null[[n]][1])) my_null[[n]] <- null_gsea(collection[[i]],reflist, w = w, np = np)
#	g <-gsea(reflist=reflist,set=collection[[i]], gsea_null= my_null[[n]], method="pareto")
	g <-gsea(reflist=reflist,set=collection[[i]],  method="pareto")

	le.size <- length(intersect(collection[[i]],g$ledge))
	gsea.res[i,] <- c(g$nes,g$p.value,NA,le.size,n)
	setTxtProgressBar(pb, cc/tot)
	}
close(pb)
gsea.res[,"fdr"] <- p.adjust(gsea.res[,"P.value"],method="fdr")
gsea.res<-gsea.res[order(gsea.res[,"P.value"]),]
return(gsea.res)
}


## compute GSEA for up and down regulons and gets the best/average NES 
# requires a collection with mode of action
multi.gsea.2t<-function(collection,reflist,minsize=15, w = 1, np = 1000,method="best"){
# method accepts "best" to select the biggest nes or "average" to use the get the 
# the average factorized by the proportion of negagive/positive targets 
newnet<-list()
for(i in names(collection)){
	pos_reg <-intersect(names(collection[[i]][which(collection[[i]] > 0)]),names(reflist))
	neg_reg <-intersect(names(collection[[i]][which(collection[[i]] < 0)]),names(reflist))
	pos_n<-length(pos_reg)
	neg_n<-length(neg_reg)
	if(pos_n > minsize & neg_n > minsize){
		newnet[[i]][["pos"]] <- c(pos_reg)
		newnet[[i]][["neg"]] <- c(neg_reg)
	}else if(pos_n > minsize & neg_n < minsize){
		newnet[[i]][["pos"]] <- c(pos_reg)
	}else if(pos_n < minsize & neg_n > minsize){
		newnet[[i]][["neg"]] <- c(neg_reg)
	}
}

toplen <- max(sapply(names(newnet), function(i) c(length(newnet[[i]][["neg"]]),length(newnet[[i]][["pos"]]) ) ))
my_null <- lapply(1:toplen,function(i) NA)

gsea.res<-data.frame(matrix(ncol=11,nrow=length(names(newnet))))
rownames(gsea.res)<-names(newnet)
colnames(gsea.res)<-c("NES+","P+","gs_size+","LE.size+","NES-","P-","gs_size-","LE.size-","NES","P.value","fdr")

pb <- txtProgressBar(style=3)
c<-0
for(i in names(newnet)){
	c <- c+1
	t <- length(names(newnet))
	n_pos<-length(newnet[[i]][["pos"]])
	n_neg<-length(newnet[[i]][["neg"]])
	if(n_pos > minsize & n_neg > minsize){
		if(is.na(my_null[[n_pos]][1])) my_null[[n_pos]] <- null_gsea(newnet[[i]][["pos"]],reflist, w = w, np = np)
		g1 <-gsea(reflist=reflist,set=newnet[[i]][["pos"]],method="pareto",w = w, np = np, gsea_null=my_null[[n_pos]])
		if(is.na(my_null[[n_neg]][1])) my_null[[n_neg]] <- null_gsea(newnet[[i]][["neg"]],reflist, w = w, np = np)		
		g2 <-gsea(reflist=reflist,set=newnet[[i]][["neg"]],method="pareto",w = w, np = np, gsea_null=my_null[[n_neg]])
		nes <- (n_pos/(n_pos+n_neg))*g1$nes + (n_neg/(n_pos+n_neg))*g2$nes
		p <- 1 - pnorm(abs(nes))
		gsea.res[i,] <- c(g1$nes,g1$p.value,n_pos,length(g1$ledge),-g2$nes,g2$p.value,n_neg,length(g2$ledge),nes,p,NA)
	}else if(n_pos > minsize & n_neg < minsize){
		if(is.na(my_null[[n_pos]][1])) my_null[[n_pos]] <- null_gsea(newnet[[i]][["pos"]],reflist, w = w, np = np)		
		g1 <-gsea(reflist=reflist,set=newnet[[i]][["pos"]],method="pareto",w = w, np = np, gsea_null=my_null[[n_pos]])
		gsea.res[i,]<-c(g1$nes,g1$p.value,n_pos,length(g1$ledge),NA,NA,NA,NA,g1$nes,g1$p.value,NA)
	}else if(n_pos < minsize & n_neg > minsize){
		if(is.na(my_null[[n_neg]][1])) my_null[[n_neg]] <- null_gsea(newnet[[i]][["neg"]],reflist, w = w, np = np)				
		g2 <-gsea(reflist=reflist,set=newnet[[i]][["neg"]],method="pareto",w = w, np = np, gsea_null=my_null[[n_neg]])
		gsea.res[i,]<-c(NA,NA,NA,NA,g2$nes,g2$p.value,n_neg,length(g2$ledge),-g2$nes,g2$p.value,NA)
	}
	setTxtProgressBar(pb, c/t)
}
close(pb)

gsea.res[,"fdr"] <- p.adjust(gsea.res[,"P.value"],method="fdr")
gsea.res<-gsea.res[order(gsea.res[,"P.value"]),]
return(gsea.res)

}


## compute GSEA for up and down regulons and gets the best/average NES 
# requires a collection with mode of action
multi.gsea.marina<-function(collection,reflist,minsize=15, w = 1, np = 1000,corcut=0){
# method accepts "best" to select the biggest nes or "average" to use the get the 
# the average factorized by the proportion of negagive/positive targets 
newnet<-list()
for(i in names(collection)){
	pos_reg <-intersect(names(collection[[i]][which(collection[[i]] > 0)]),names(reflist))
	neg_reg <-intersect(names(collection[[i]][which(collection[[i]] < 0)]),names(reflist))
	pos_n<-length(pos_reg)
	neg_n<-length(neg_reg)
	if(pos_n > minsize & neg_n > minsize){
		newnet[[i]][["pos"]] <- c(pos_reg)
		newnet[[i]][["neg"]] <- c(neg_reg)
	}else if(pos_n > minsize & neg_n < minsize){
		newnet[[i]][["pos"]] <- c(pos_reg)
	}else if(pos_n < minsize & neg_n > minsize){
		newnet[[i]][["neg"]] <- c(neg_reg)
	}
}

#toplen <- max(sapply(names(newnet), function(i) c(length(newnet[[i]][["neg"]]),length(newnet[[i]][["pos"]]) ) ))
#my_null <- lapply(1:toplen,function(i) NA)

gsea.res <- data.frame(matrix(ncol=8,nrow=length(names(newnet))))
rownames(gsea.res)<-names(newnet)
colnames(gsea.res)<-c("ES","NES","P.value","fdr","ES+","NES+","ES-","NES-")

pb <- txtProgressBar(style=3)
c<-0
for(i in names(newnet)){
	c <- c+1
	t <- length(names(newnet))
	n_pos<-length(newnet[[i]][["pos"]])
	n_neg<-length(newnet[[i]][["neg"]])
	if(n_pos > minsize & n_neg > minsize){
		#if(is.na(my_null[[n_pos]][1])) my_null[[n_pos]] <- null_gsea(newnet[[i]][["pos"]],reflist, w = w, np = np)
		g1 <-gsea(reflist=reflist,set=newnet[[i]][["pos"]],method="pareto",w = w, np = np)
		#if(is.na(my_null[[n_neg]][1])) my_null[[n_neg]] <- null_gsea(newnet[[i]][["neg"]],reflist, w = w, np = np)		
		g2 <-gsea(reflist=reflist,set=newnet[[i]][["neg"]],method="pareto",w = w, np = np)
		if(abs(g1$nes) >= abs(g2$nes)){
			gsea.res[i,] <- c(g1$es,g1$nes,g1$p.value,NA,g1$es,g1$nes,-g2$es,-g2$nes)
		}else{
			gsea.res[i,] <- c(-g2$es,-g2$nes,g1$p.value,NA,g1$es,g1$nes,-g2$es,-g2$nes)
			}
		
	}else if(n_pos > minsize & n_neg < minsize){
		#if(is.na(my_null[[n_pos]][1])) my_null[[n_pos]] <- null_gsea(newnet[[i]][["pos"]],reflist, w = w, np = np)		
		g1 <-gsea(reflist=reflist,set=newnet[[i]][["pos"]],method="pareto",w = w, np = np)
		gsea.res[i,]<-c(g1$es,g1$nes,g1$p.value,NA,g1$es,g1$nes,NA,NA)
	}else if(n_pos < minsize & n_neg > minsize){
		#if(is.na(my_null[[n_neg]][1])) my_null[[n_neg]] <- null_gsea(newnet[[i]][["neg"]],reflist, w = w, np = np)				
		g2 <-gsea(reflist=reflist,set=newnet[[i]][["neg"]],method="pareto",w = w, np = np)
		gsea.res[i,]<-c(g2$es,g2$nes,g2$p.value,NA,NA,NA,-g2$es,-g2$nes)
	}
	setTxtProgressBar(pb, c/t)
}
close(pb)
gsea.res[,"fdr"] <- p.adjust(gsea.res[,"P.value"],method="fdr")
gsea.res <- gsea.res[order(gsea.res[,"P.value"]),]
return(gsea.res)
}

## compute GSEA for up and down regulons separatedly 

multi.gsea.2tv2<-function(collection,reflist,minsize=10,corcut=0.1,w=1,np=1000){
newnet<-list()
for(i in names(collection)){
	pos_reg <-intersect(names(collection[[i]][which(collection[[i]] > corcut)]),names(reflist))
	neg_reg <- intersect(names(collection[[i]][which(collection[[i]] < -corcut)]),names(reflist))
	pos_n<-length(pos_reg)
	neg_n<-length(neg_reg)
	if(pos_n > minsize & neg_n > minsize){
		newnet[[i]][["pos"]] <- c(pos_reg)
		newnet[[i]][["neg"]] <- c(neg_reg)
	}else if(pos_n > minsize & neg_n < minsize){
		newnet[[i]][["pos"]] <- c(pos_reg)
	}else if(pos_n < minsize & neg_n > minsize){
		newnet[[i]][["neg"]] <- c(neg_reg)
	}
}

toplen <- max(sapply(names(newnet), function(i) c(length(newnet[[i]][["neg"]]),length(newnet[[i]][["pos"]]) ) ))
my_null <- lapply(1:toplen,function(i) NA)

regs.pos <- unlist(sapply(names(newnet), function(i) length(newnet[[i]][["pos"]])))
regs.neg <- unlist(sapply(names(newnet), function(i) length(newnet[[i]][["neg"]])))
reg.names<-c(paste(names(unlist(regs.pos)),".pos",sep=""),paste(names(unlist(regs.neg)),".neg",sep=""))
regs<-c(regs.pos,regs.neg)
names(regs)<-reg.names
regs <- regs[which(regs > minsize)]

toplen <- max(regs)
my_null <- lapply(1:toplen,function(i) NA)

geneId <- NES <- Mode <- P.value <- gs_size <- LE.size <- fdr <- rep(NA,length(names(regs)))
gsea.res <- data.frame(geneId , NES , Mode , P.value , gs_size , LE.size, fdr)
rownames(gsea.res)<-names(regs)

pb <- txtProgressBar(style=3)
c<-0
for(i in names(newnet)){
	c <- c+1
	t <- length(names(newnet))
	n_pos<-length(newnet[[i]][["pos"]])
	n_neg<-length(newnet[[i]][["neg"]])
	ipos<-paste(i,"pos",sep=".")
	ineg<-paste(i,"neg",sep=".")		
	if(ipos %in% names(regs)){
		if(is.na(my_null[[n_pos]][1])) my_null[[n_pos]] <- null_gsea(newnet[[i]][["pos"]],reflist, w = w, np = np)				
		g1 <-gsea(reflist=reflist,set=newnet[[i]][["pos"]],method="pareto",w=w,np=np)
		gsea.res[ipos,]<-c(i,g1$nes,"pos",g1$p.value,n_pos,length(g1$ledge),NA)
	}
	if(ineg %in% names(regs)){
		if(is.na(my_null[[n_neg]][1])) my_null[[n_neg]] <- null_gsea(newnet[[i]][["neg"]],reflist, w = w, np = np)				
		g2 <-gsea(reflist=reflist,set=newnet[[i]][["neg"]],method="pareto",w=w,np=np)
		gsea.res[ineg,]<-c(i,g2$nes,"neg",g2$p.value,n_neg,length(g2$ledge),NA)
	}
	setTxtProgressBar(pb, c/t)
}
close(pb)

gsea.res[,"fdr"] <- p.adjust(gsea.res[,"P.value"],method="fdr")
gsea.res<-gsea.res[order(as.numeric(gsea.res[,"P.value"])),]
return(gsea.res)
}




##############
# retrieve le genset genes

ledge_plot <- function(gobj,numcol=6,let.size=2.5,annot=NULL){

library(reshape2) ## melt
library(plyr)     ## round_any
library(ggplot2)

rounded <- ceiling((length(gobj$set)+2)/numcol)*numcol
numrow <- rounded/numcol
numna <- rounded-(length(gobj$set)+2)

leset <- intersect(gobj$ledge,gobj$set)
poset <- names(which(gobj$reflist[gobj$set] >= 0))
neset <- names(which(gobj$reflist[gobj$set] < 0))

if(gobj$es > 0 ){
s1 <- c(quantile(gobj$reflist,0.99),sort(gobj$reflist[leset],decreasing=TRUE))
s2 <- sort(gobj$reflist[setdiff(poset,leset)],decreasing=TRUE)
s3 <- c(sort(gobj$reflist[neset],decreasing=TRUE),quantile(gobj$reflist,0.01))
colv <- c(rep("black",length(s1)),rep("grey",length(s2)),rep("grey",numna),rep("grey",length(s3)))
}else if(gobj$es < 0){
s1 <- c(quantile(gobj$reflist,0.01),sort(gobj$reflist[leset]))
s2 <- sort(gobj$reflist[setdiff(neset,leset)])
s3 <- c(sort(gobj$reflist[poset]),quantile(gobj$reflist,0.99))
colv <- c(rep("black",length(s1)),rep("grey",length(s2)),rep("grey",numna),rep("grey",length(s3)))
}
dat <- expand.grid(var1=1:numcol, var2=numrow:1)
dat$value <- c(s1,s2,rep(0,numna),s3)
if(is.null(annot)){ dat$labels <- c(names(s1),names(s2),rep("",numna),names(s3)) 
}else{
	dat$labels <- paste(c(annot[names(s1)],annot[names(s2)],rep("",numna),annot[names(s3)]),c(names(s1),names(s2),rep("",numna),names(s3)),sep="\n")  }

dat$color <- colv
dat$colorScale  <- val2col(dat$value)

ggplot(dat, aes(x=var1,y=var2,label=labels),main="",colour='white')+ 
geom_tile(aes(fill = value),colour=colv) +
scale_fill_gradient2(low="blue",mid="white", high="red") +
geom_text(size=let.size)

}
  
############
es.gsea<-function(mexp.rank=NULL, set=NULL,sizelim=10){
res<-c(1:ncol(mexp.rank))
names(res)<-colnames(mexp.rank)
res<-sapply(colnames(mexp.rank), function(i) gsea_es(mexp.rank[,i],set,sizelim=sizelim))
return(res)
}

my.nes <- function(mrank=NULL, set=NULL,np=10000,w=1){
null <- null_gsea_mat(set,mrank,w=w,np=np)
res <- c(1:ncol(mrank))
names(res)<-colnames(mrank)
res <- sapply(colnames(mrank), function(i) gsea_es(mrank[,i],set))
out <- sapply(names(res),function(i) (res[i]-mean(null))/sd(null))
return(out)
}


my.nes2 <- function(mrank=NULL, set1=NULL,set2=NULL,np=10000,w=1){
null1 <- null_gsea_mat(set1,mrank,w=w,np=np)
null2 <- null_gsea_mat(set2,mrank,w=w,np=np)
res1 <- res2 <- c(1:ncol(mrank))
names(res1)<-names(res2) <- colnames(mrank)
res1 <- sapply(colnames(mrank), function(i) gsea_es(mrank[,i],set1))
res2 <- sapply(colnames(mrank), function(i) gsea_es(mrank[,i],set2))
out1<-sapply(names(res1),function(i) (res1[i]-mean(null))/sd(null))
out2<-sapply(names(res2),function(i) (res2[i]-mean(null))/sd(null))
pos_div <- length(set1)/(length(set1)+length(setneg))
neg_div <- length(set2)/(length(set2)+length(setneg))
out <- out1*pos_div - out2*neg_div
return(out)
}


nes<-function(exprank=NULL,set=NULL,np=100,null.type="collapsed",min.size=10){
## calculates GSEA Normalized enrichment Score for all signatures provided as columns in a matrix 
# exprank ; A matrix wich columns are ranked/reference list for GSEA
# set ; a list of elements (from the rownames cof exprank)
# null.type=c("collapsed","individual")
# np ; if small number ~ 100 better use collapsed method

set <- intersect(set,rownames(exprank))
if(length(set) < min.size) {warning(paste("Gene set length:",length(set),"is too short"))}

esSet <- es.gsea(exprank,set)
null_es<-list()
pb <- txtProgressBar(style=3)
for(i in 1:np){
	rSet <- sample(rownames(exprank))[1:length(set)]
	null_es[[i]] <- es.gsea(exprank,rSet)
	setTxtProgressBar(pb, i/np)
	}
close(pb)

if(null.type == "individual"){
null_es <- do.call(rbind,null_es)
nesSet <- sapply(colnames(null_es), function(i) (esSet[i]-mean(null_es[,i]))/sd(null_es[,i]) )
}else if(null.type == "collapsed"){
null_es <- unlist(null_es)
nesSet <- (esSet-mean(null_es))/sd(null_es)
}
return(nesSet)
}

es.gsea2<-function(mexp.rank=NULL, setpos=NULL,setneg=NULL){
setpos<-intersect(setpos,rownames(mexp.rank))
setneg<-intersect(setneg,rownames(mexp.rank))
pos_div <- length(setpos)/(length(setpos)+length(setneg))
neg_div <- length(setneg)/(length(setpos)+length(setneg))

respos<- sapply(colnames(mexp.rank), function(i) gsea_es(mexp.rank[,i],setpos))
resneg<- sapply(colnames(mexp.rank), function(i) gsea_es(mexp.rank[,i],setneg))
res <- respos*pos_div - resneg*neg_div
return(res)
}

nes2 <- function(mrank=NULL, setpos=NULL,setneg=NULL,np=100){
setpos <- intersect(setpos,rownames(mrank))
setneg <- intersect(setneg,rownames(mrank))
pos_div <- length(setpos)/(length(setpos)+length(setneg))
neg_div <- length(setneg)/(length(setpos)+length(setneg))
respos<- nes(mrank,setpos,np=np)
resneg<- nes(mrank,setneg,np=np)
res <- respos*pos_div - resneg*neg_div
return(res)
}


ES.matrix <- function (mexp.rank=NULL,net=NULL,tflist=tflist,corcut=0)
{
  # TF activity calculation
  # mexp.rank (matrix): a matrix of zrank transformed data
  # collection (list): a collection of regulons with mode of action
  # tflist (vector): a subset of names from collection to be ran
  #
es_mat<-matrix(nrow=length(tflist),ncol=ncol(mexp.rank))
colnames(es_mat)<-colnames(mexp.rank)
rownames(es_mat)<-tflist
genelist<-rownames(mexp.rank)
for(i in tflist){
	possig <-intersect(names(net[[i]][which(net[[i]] > corcut)]),genelist)
	negsig <-intersect(names(net[[i]][which(net[[i]] < corcut)]),genelist)
	print(paste(i,length(negsig),length(possig),sep="  "))
	if(length(possig) > 0 && length(negsig) > 0){
		es_pos <- es.gsea(mexp.rank,possig)*length(possig)/(length(possig)+length(negsig))
		es_neg <- -1* es.gsea(mexp.rank,negsig)*length(negsig)/(length(possig)+length(negsig))
		es_mat[i,] <- es_neg+es_pos
		}
	else if(length(possig) > 0 && length(negsig) == 0){
		es_mat[i,] <- es.gsea(mexp.rank,possig)
		}
	else if(length(possig) == 0 && length(negsig) > 0){
		es_mat[i,] <- -1* es.gsea(mexp.rank,negsig)
		}
	}
return(es_mat)
}


ES.geneset <- function (mexp.rank=NULL,collection=NULL,setnames=NULL)
{
  # Pathway activity calculation
  # mexp.rank (matrix): a matrix of zrank transformed data
  # collection (list): a collection of gene sets
  # setnames (vector): a subset of names from collection to be ran
  #
es_mat<-matrix(nrow=length(setnames),ncol=ncol(mexp.rank))
colnames(es_mat)<-colnames(mexp.rank)
rownames(es_mat)<-setnames
genelist<-rownames(mexp.rank)
for(i in setnames){
	possig <-intersect(collection[[i]],genelist)
	print(paste(i,length(possig),sep="  "))
	es_mat[i,] <- es.gsea(mexp.rank,possig)
	}
return(es_mat)
}


### Barplot for a pathway enrichment analisys (obsolet)

ea.barplot<-function(ea.table=ea.table,filename=filename, param="P.value",fdrcut=1,fdrline=0.05,pvalcut=0.05,maintitle="Enrichment analysis"){
# param = select the parameter to represent in the plot
#         if p.value the -logP with the sign of the NES will be showed

if(length(ea.table[,"NES"][which(ea.table[,"NES"] < 0 )]) >= 25 ){  resneg<-names(sort(ea.table[,"NES"])[1:25])
}else{ resneg<- names(sort(ea.table[,"NES"][which(ea.table[,"NES"] < 0 )]))}

if(length(ea.table[,"NES"][which(ea.table[,"NES"] > 0 )]) >= 25 ){  respos<-names(sort(ea.table[,"NES"],decreasing=TRUE)[1:25])
}else{ respos<-names(sort(ea.table[,"NES"][which(ea.table[,"NES"] > 0 )]))}
restrict<-c(respos,resneg)

sigres<-names(ea.table[restrict,"fdr"][which(ea.table[restrict,"fdr"] < fdrcut)])
sigres<-names(ea.table[restrict,"P.value"][which(ea.table[restrict,"P.value"] < pvalcut)])

if(length(sigres < 4)) sigres <- names(sort(ea.table[restrict,"P.value"])[1:4])
nes <- ea.table[sigres,"NES"]
color <- psign <- sign(nes)

color[which(color == 1)] <- "red"
color[which(color == -1)] <- "blue"

if(param == "P.value"){
	ppar <- -log10(ea.table[sigres,param])
	names(ppar)<-sigres
	param_lab<-"-log(p.value)"
}else if(param == "NES"){
	ppar <- abs(nes)
	param_lab<-"Normalized enrichment score"
}else{ stop("Wrong param !!!")}

label_length<-floor(max(sapply(sigres,function(i) nchar(i)))/3)
wp <- ceiling(50 + label_length)/10
hp <- ceiling(length(sigres)/6 )
if(hp < 2) hp <-2
pvcut<- -log10(max(ea.table[,"P.value"][which(ea.table[,"fdr"] < fdrline)]))

if(label_length > 25){ label_length = 25
}else if(label_length < 10) {label_length =4}

dev.new(file=filename,width=wp,height=hp)
par(mar=c(3,label_length,3,1))
barplot(ppar[names(sort(nes))], horiz=TRUE, col=color[names(sort(nes))], xlab=param_lab,cex.names=0.5,las=1,main=maintitle)
abline(a=NULL,b=NULL,v=pvcut)
savefont <- par(font=3)
legend(x=pvcut,y=length(ppar), paste("fdr = ",fdrline,sep=""),bty="n")
par(savefont)
axis(2,cex.axis=0.1)
dev.off()
}
