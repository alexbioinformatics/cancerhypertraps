#BiocManager::install('GEOquery')
library(GEOquery)

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

# Boxplots
palette(c("#c7ff9d","#f4dfdf","#f2cb98","#dfeaf4","#f3f388"))
#dev.new(width=4+dim(gset)[[2]]/5, height=6)
par(mar=c(2+round(max(nchar(sampleNames(gset)))/2),4,2,1))
title <- paste ("GSE62932", '/', annotation(gset), " All samples", sep ='')
boxplot(ex, boxwex=0.6, notch=T, main=title, outline=FALSE, las=2, col=sml)
legend("topleft", c("Control","Stage 1","Stage 2","Stage 3","Stage 4"), fill=palette(), bty="n")

# Define sub-groups
control <- ex[,sample.annotation$Sample_Type=='Control']
stage1 <- ex[,sample.annotation$Sample_Type=='Stage1']
stage2 <- ex[,sample.annotation$Sample_Type=='Stage2']
stage3 <- ex[,sample.annotation$Sample_Type=='Stage3']
stage4 <- ex[,sample.annotation$Sample_Type=='Stage4']

# Differential Analysis
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

# Extract the significant genes from each stage.
top.stage1 <- stage1.tT[which(stage1.tT$adj.P.Val<0.01),]
top.stage2 <- stage2.tT[which(stage2.tT$adj.P.Val<0.01),]
top.stage3 <- stage3.tT[which(stage3.tT$adj.P.Val<0.01),]
top.stage4 <- stage4.tT[which(stage4.tT$adj.P.Val<0.01),]

# Venn diagram
library("VennDiagram")
# Produce venn diagram to ascertain intersection between groups of significant genes. 
v1 <- venn.diagram(
  list('Stage 1'=row.names(top.stage1),
       'Stage 2'=row.names(top.stage2), 
       'Stage 3'=row.names(top.stage3), 
       'Stage 4'=row.names(top.stage4)),
  filename=NULL, fill=rainbow(4), fontfamily='Arial', cat.fontfamily = 'Arial')
grid.newpage()
grid.draw(v1)

# Store the genes attributed to each of the 15 intersections.
venn.partitions <- get.venn.partitions(list('Stage 1'=row.names(top.stage1),
                                            'Stage 2'=row.names(top.stage2), 
                                            'Stage 3'=row.names(top.stage3), 
                                            'Stage 4'=row.names(top.stage4)))

# Grouping
allsamples <- cbind(stage1,stage2,stage3,stage4) # Produce data frame of CRC samples.

# Produce dataframe of genes and groups they pertain to. 
groups=NULL
for (i in 1:nrow(venn.partitions)){
  tmp <- data.frame(gene=venn.partitions$..values..[[i]],group=i)
  groups=rbind(groups,tmp)
}

# Feature assignment function
HT_conversion <- function(data, features, stringency=0.5) {
  features$group <- as.factor(features$group) # Convert feature number to factor.
  no_samples <- ncol(data) # Extract number of samples.
  no_features <- length(levels(features$group)) # Number of features.
  
  # Produce matrix of zeros to store feature assignment. 
  binary <- matrix(rep(0,no_features*no_samples), ncol = no_features, nrow = no_samples)
  
  for (sample in 1:no_samples) {
    # Extract control samples and each CRC sample in turn.
    sampletest <- cbind(control, allsamples[,sample]) 
    
    # Run limma differntial expression pipeline on extracted data.
    values <- c(rep('control', ncol(control)),'sampletest')
    values <- as.factor(values)
    design=model.matrix(~0+values)
    colnames(design) <- levels(values)
    fit <- lmFit(sampletest, design)
    cont.matrix <- makeContrasts(sampletest-control, levels=design)
    fit2 <- contrasts.fit(fit, cont.matrix)
    fit2 <- eBayes(fit2, 0.01)
    sampletest.tT <- topTable(fit2, adjust="fdr", sort.by="B", number=nrow(sampletest))
    for (feature in 1:no_features) {
      tmp.group <- subset(features, group==feature) # Extract gene names for currently iterated feature. 
      tmp <- 0 # Define tmp as 0.
      for (genes in 1:nrow(tmp.group)) { # For each gene in currently iterated feature.
        # Extract the adjusted p-value for the current gene comparison.
        adj.p.val <- sampletest.tT$adj.P.Val[as.character(rownames(sampletest.tT)) %in% tmp.group$gene[genes]] # Extract expression data for each gene in group.
        if (adj.p.val<0.05) { # If this value is significant (<0.05) add 1 to tmp.
          tmp <- tmp + 1 
        } 
      }  
      if (tmp/nrow(tmp.group)>stringency) { # Check if prop. of genes in group significantly different from control is above stringency value.
        binary[sample,feature] <- 1 # If so store value of 1 for the currently iterated sample & feature.
        print(paste('Sample: ', sample,', Feature ', feature, ' Present', sep = ''))
      } else {
        binary[sample,feature] <- 0 # Otherwise store value of 0.
        print(paste('Sample: ', sample,', Feature ', feature, ' Absent', sep = ''))
      }
    }
  }
  out <- NULL
  for (i in 1:nrow(binary)) {
    # Append intial state of all zeros between each samples sample assignment. 
    ini <- rep(1, no_features)
    end <- binary[i,]
    tmp <- rbind(ini, out)
    out <- rbind(out, tmp)
  }
  # Save output file.
  write.table(out, file="out.txt", row.names=FALSE, col.names=FALSE)
  return(as.data.frame(out))
} 

# Run function. 
HT_conversion(allsamples,groups,0.75)

