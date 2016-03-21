library("GO.db")


getgochildrens <- function(goterms , goroot="BP"){
# goroot accepted values: BP,CC,MF
#

# Convert the object to a list
if(goroot == "BP"){xx <- as.list(GOBPCHILDREN)
}else if(goroot == "CC"){xx <- as.list(GOCCCHILDREN)
}else if(goroot == "MF"){xx <- as.list(GOMFCHILDREN)
}else{message("error; wrong GO root; accepted values: BP,CC,MF");stop;}
# Remove GO IDs that do not have any children
xx2 <- xx[!is.na(xx)]

# retrieve GO childs for the provided GO term
golist <- unname(unlist(sapply(goterms,function(i)  unique(unname(xx[[i]][which(names(xx[[i]]) == "is_a")])))))
# retrieve GO childs with childs for the provided GO term
withchild <- unname(unlist(sapply(goterms,function(i)  unique(unname(xx2[[i]][which(names(xx2[[i]]) == "is_a")])))))

# repeat for "withchild" terms until no childs with childs remain and save the list in "resGoList"

resGoList <- c(goterms,golist)
if(length(withchild) > 0){
pb <- txtProgressBar(style=3)
cnt<-0
repeat{
	golist <- unname(unlist(sapply(withchild,function(i)  unique(unname(xx[[i]][which(names(xx[[i]]) == "is_a")])))))
	withchild <- unname(unlist(sapply(withchild,function(i)  unique(unname(xx2[[i]][which(names(xx2[[i]]) == "is_a")])))))
	#message(golist)
	resGoList<-c(resGoList,golist)
	cnt<-cnt+1
	setTxtProgressBar(pb, cnt/length(withchild))
  if(length(withchild) == 0){
    break
  }
}
close(pb)
}
return(unique(resGoList))
}

