library("GEOquery")
source("/ifs/data/c2b2/ac_lab/gl2377/soft/R/my_stat_functions.r")

## Downloading a GEO series RNA-seq data
# gets exp matrix and  
gse <- getGEO("GSE49710",GSEMatrix=FALSE)
show(gse)
head(Meta(gse))
gsmplatforms <- lapply(GSMList(gse),function(x){Meta(x)$platform})
Table(GSMList(gse)[[1]])[1:5,]
Columns(GSMList(gse)[[1]])[1:5,]
probesets <- Table(GPLList(gse)[[1]])$ID

data.matrix <- do.call('cbind',lapply(GSMList(gse),function(x){
	tab <-Table(x)
	mymatch <- match(probesets,tab$ID_REF)
	return(tab$VALUE[mymatch])}))
data.matrix <-apply(data.matrix,2,function(x) {as.numeric(as.character(x))})
rownames(data.matrix) <-probesets
colnames(data.matrix) <- names(GSMList(gse))

sample_info <- data.frame(sapply(pheno@data[,1],function(i)  
	sapply( gse@gsms[[i]]@header$characteristics_ch1, function(j)
		strsplit(j," ")[[1]][length( strsplit(j," ")[[1]])] 
		)))
colnames(sample_info) <- as.character(pheno@data[,1])
rownames(sample_info) <-  sapply(gse@gsms[[1]]@header$characteristics_ch1 , function(i) strsplit(i,":")[[1]][1])

platform <- names(GPLList(gse))
gpl <- getGEO(platform)
rna.sym <- data.matrix[which(gpl@dataTable@table$GeneSymbol != ""),]
rownames(rna.sym) <- gpl@dataTable@table$GeneSymbol[which(gpl@dataTable@table$GeneSymbol != "")]
rna.sym <- filterCV(rna.sym)


## Save data for NBL project
r_mexp.sym <- rna.sym
r_clinical<-list(
r_stage = as.numeric(sample_info["inss stage",]),
r_mycn = as.numeric(sample_info["mycn status",]) - 1,
r_risk = as.numeric(sample_info["high risk",]) -1,
r_death = as.numeric(sample_info["death from disease",]),
r_progr = as.numeric(sample_info["progression",]),
r_age = as.numeric(sample_info["age at diagnosis" ,])
)
expmat <- cbind(rownames(r_mexp.sym),r_mexp.sym)
colnames(expmat)<-c("gene",colnames(r_mexp.sym))
write.table(expmat,file="/ifs/scratch/c2b2/ac_lab/gl2377/Interactome/NBL/GSE49710_rna-seq/aracne/rnaseq_nbl_498.exp",quote=FALSE,row.names=FALSE,sep="\t")


names(r_clinical$r_stage) <- names(r_clinical$r_mycn) <- names(r_clinical$r_risk) <- names(r_clinical$r_death) <- names(r_clinical$r_progr) <- names(r_clinical$r_age) <- colnames(sample_info)
save(r_mexp.sym,r_clinical,sample_info,file="/ifs/scratch/c2b2/ac_lab/gl2377/Interactome/NBL/rdata/r_GSE49710.rda")

