#setwd("")

library(Biobase)
library(GEOquery)
library(limma)
library(pheatmap)
library(ggplot2)
library(plyr)
library(reshape2)
library(gplots)

series <-"GSE48558"
platform <- "GPL6244"

###Loading Data

gset <- getGEO(series, GSEMatrix =TRUE, AnnotGPL=TRUE, destdir = "Data/")
if (length(gset) > 1){idx <- grep(platform, attr(gset, "names"))
} else {idx <- 1 }
gset <- gset[[idx]]

gset<- gset[,which(gset$`phenotype:ch1` == "Normal" | gset$source_name_ch1 == "AML Patient")]

CreateLabel <- function(x) {
  if (gset$source_name_ch1[x] == "AML Patient") {
    return("Test")
  } else {
    return("Normal")
  }
}

gr <- sapply(1:length(gset$`phenotype:ch1`) , CreateLabel)

ex <-exprs(gset)

min_ex <- min(ex)
max_ex <- max(ex)

#print(c(min_ex, max_ex))


pdf("Results/boxplot.pdf", width = 20)
boxplot(ex)
dev.off()


#normalize
#ex<- normalizeQuantiles(ex)
#ex(gset) <- ex



###Correlation Heatmap
pdf("Results/CorHeatmap.pdf", width = 25, height = 15)
pheatmap(cor(ex), labels_row =  gr , labels_col =  gr, color = bluered(256))
dev.off()


###Principal Component Analysis (PCA1)
pc <- prcomp(ex)
pdf("Results/PC.pdf")
plot(pc)
plot(pc$x[,1:2])
dev.off()



###Principal Component Analysis (PC Scaled)
ex.scaled <- t(scale(t(ex), scale = F))
pcScaled <- prcomp(ex.scaled)
pdf("Results/PC_Scaled.pdf")
plot(pcScaled)
plot(pcScaled$x[,1:2])
dev.off()


pcr <- data.frame(pc$r[,1:3] , Group = gr)
pdf("Results/PCA_samples.pdf")
ggplot(pcr, aes(PC1, PC2, color=Group)) + geom_point(size=3)+ theme_bw()
dev.off()
#dim(pc$rotation)


### Differential Expression Analysis


gr <- c(rep("Test", 13), rep("other", 86), "CD34", rep("other", 3), "CD34", rep("other", 3), "CD34", rep("other", 36), rep("Test", 2), "other", rep("Test", 3), rep("other", 20))
gr <- factor(gr)

gset$description <- gr
design <- model.matrix(~ description + 0, gset)
colnames(design) <- levels(gr)
fit <- lmFit(gset, design)

cont.matrix <- makeContrasts(Test - CD34 , levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2, 0.01)

tT <- topTable(fit2, adjust="fdr", sort.by="B", number=Inf)
tT <- subset(tT, select=c("ID","adj.P.Val","P.Value","t","B","logFC","Gene.symbol","Gene.title"))
tT <- subset(tT , select = c("Gene.symbol" , "Gene.title", "adj.P.Val"  , "logFC"))

###aml up genes

aml.up <- subset(tT, logFC > 1 & adj.P.Val < 0.05)
aml.up.genes <-unique( as.character(strsplit2( (aml.up$Gene.symbol),"///")))
write.table(aml.up.genes, file = "Results/Test_CD34_Up.txt", row.names = F, col.names = F)

###aml down genes
aml.down <- subset(tT, logFC < -1 & adj.P.Val < 0.05)
aml.down.genes <-unique( as.character(strsplit2( (aml.down$Gene.symbol),"///")))
write.table(aml.down.genes, file = "Results/Test_CD34_down.txt", row.names = F, col.names = F )
