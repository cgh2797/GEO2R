### Human+anmial brain disease mRNA data collected from geo

For ["Transient Potassium Channels: Therapeutic Targets for Brain Disorders"](https://www.frontiersin.org/articles/10.3389/fncel.2019.00265/full)

## GEO
GEO is a public functional genomics data repository supporting MIAME-compliant data submissions. [link](https://www.ncbi.nlm.nih.gov/geo/)

## GEO2R
gene data analysis tool

https://www.ncbi.nlm.nih.gov/geo/info/geo2r.html

[youtube guide](https://www.youtube.com/watch?v=EUPmGWS8ik0)

## example
```
#   Differential expression analysis with limma
library(Biobase)
library(GEOquery)
library(limma)

# load series and platform data from GEO

gset <- getGEO("GSE6834", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL4757", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

# make proper column names to match toptable 
fvarLabels(gset) <- make.names(fvarLabels(gset))

# group names for all samples
gsms <- "00000000001111111111XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX"
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }

# eliminate samples marked as "X"
sel <- which(sml != "X")
sml <- sml[sel]
gset <- gset[ ,sel]

# log2 transform
ex <- exprs(gset)
qx <- as.numeric(quantile(ex, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm=T))
LogC <- (qx[5] > 100) ||
          (qx[6]-qx[1] > 50 && qx[2] > 0) ||
          (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2)
if (LogC) { ex[which(ex <= 0)] <- NaN
  exprs(gset) <- log2(ex) }

# set up the data and proceed with analysis
sml <- paste("G", sml, sep="")    # set group names
fl <- as.factor(sml)
gset$description <- fl
design <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(fl)
fit <- lmFit(gset, design)
cont.matrix <- makeContrasts(G1-G0, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
# show top 10000 gene
tT <- topTable(fit2, adjust="fdr", sort.by="B", number=10000)

tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title"))
# write.table(tT, file=stdout(), row.names=F, sep="\t")
# stdout()

write.table(tT,file="C:/result1.txt" ,row.names=F, sep="\t")
# save as txt

```

## Limma
[Limma](https://bioconductor.org/packages/release/bioc/vignettes/limma/inst/doc/usersguide.pdf) : Linear Models for Microarray and RNA-Seq Data User’s Guide

Limma is a package for the analysis of gene expression data arising from microarray or RNA-Seq
technologies (Package for R)

## ETC
### [GeneCards](https://www.genecards.org/)

### [GeneMANIA](https://genemania.org/)
