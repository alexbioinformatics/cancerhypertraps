#BiocManager::install('GEOquery')
library(GEOquery)
library(limma)

# Load data from GEO
gset <- getGEO("GSE62932", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL570", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

#Extract Gene annotation
gene.annotation <- as.data.frame(gset@featureData@data)
gene.annotation <- subset(gene.annotation, select=c('ID','Gene symbol','Gene title', 'Gene ID'))

# Manual sample grouping (each number refers to which data group each column belongs to)
gsms <- "33224434223431434111241441323323112243243233224421234433131331220000"
sml <- c()
for (i in 1:nchar(gsms)) { sml[i] <- substr(gsms,i,i) }

# order samples by group
ex <- as.data.frame(exprs(gset)[, order(sml)])
sml <- sml[order(sml)]
sml <- as.factor(sml)
levels(sml) <- c("Control","Stage1","Stage2","Stage3","Stage4")
sample.annotation <- data.frame(Sample_ID=colnames(ex),Sample_Type=sml)

control <- ex[,sample.annotation$Sample_Type=='Control']
stage1 <- ex[,sample.annotation$Sample_Type=='Stage1']
stage2 <- ex[,sample.annotation$Sample_Type=='Stage2']
stage3 <- ex[,sample.annotation$Sample_Type=='Stage3']
stage4 <- ex[,sample.annotation$Sample_Type=='Stage4']


library(limma)
test1 <- cbind(control, stage1)
values <- c(rep('control', ncol(control)),rep('stage1',ncol(stage1)))
values <- as.factor(values)
design=model.matrix(~0+values)
colnames(design) <- levels(values)
fit <- lmFit(test1, design)
cont.matrix <- makeContrasts(stage1-control, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
stage1.tT <- topTable(fit2, adjust="fdr", sort.by="B", number=nrow(test1))

test2 <- cbind(control, stage2)
values <- c(rep('control', ncol(control)),rep('stage2',ncol(stage2)))
values <- as.factor(values)
design=model.matrix(~0+values)
colnames(design) <- levels(values)
fit <- lmFit(test2, design)
cont.matrix <- makeContrasts(stage2-control, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
stage2.tT <- topTable(fit2, adjust="fdr", sort.by="B", number=nrow(test2))

test3 <- cbind(control, stage3)
values <- c(rep('control', ncol(control)),rep('stage3',ncol(stage3)))
values <- as.factor(values)
design=model.matrix(~0+values)
colnames(design) <- levels(values)
fit <- lmFit(test3, design)
cont.matrix <- makeContrasts(stage3-control, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
stage3.tT <- topTable(fit2, adjust="fdr", sort.by="B", number=nrow(test3))

test4 <- cbind(control, stage4)
values <- c(rep('control', ncol(control)),rep('stage4',ncol(stage4)))
values <- as.factor(values)
design=model.matrix(~0+values)
colnames(design) <- levels(values)
fit <- lmFit(test4, design)
cont.matrix <- makeContrasts(stage4-control, levels=design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)
stage4.tT <- topTable(fit2, adjust="fdr", sort.by="B", number=nrow(test4))




