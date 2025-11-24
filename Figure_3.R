library("tidyverse")
library("ggpubr")
library("RColorBrewer")
library("Hmisc")
library("DESeq2")
library("METAFlux")

path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(path)

### DE analysis -----

### Input reading ----
sdata <- read.delim("0_Input_metadata_RNA-Seq.txt", row.names=1)
countdata <- read.delim("0_Input_Tissue_RNA-Seq_raw.tsv", row.names=1, check.names = F)
TPM <- read.delim("0_Input_Tissue_RNA-Seq_TPM.tsv", row.names=1, check.names = F)

countdata <- countdata[, row.names(sdata)]
TPM <- TPM[, row.names(sdata)]

### DESeq ----- 
refcov <- "Class"
covtomodel <- c("Subject")

#### Identification of the levels of the reference covariate
cov_class <- levels(as.factor(sdata[,refcov]))

# Definition of the models
fullmodel <- formula(paste("~", paste(refcov, paste(covtomodel, collapse = "+"), sep="+")))

## Sample selection
dds_complete <- DESeqDataSetFromMatrix(countData = countdata, 
                                       colData = sdata, design = fullmodel)

dds_complete <- DESeq(dds_complete, parallel=TRUE)

res_sub <- results(dds_complete, contrast=c(refcov, "CRC", "Adjacent"), cooksCutoff=FALSE)

write.table(res_sub, "2.06_Table_DEG.tsv", sep="\t", quote=F)

### Median levels ----
dat <- data.frame(t(TPM), Class = sdata$Class)
dat <- gather(dat, Gene, Level, -Class)
dat <- group_by(dat, Gene, Class) %>% summarize(med=median(Level))

out <- pivot_wider(dat, names_from = Class, values_from = med)

write.table(out, "2.2_Table_Median_TPM.tsv", sep="\t", quote=F, row.names=F)

## Figure 3A -- Volcano ----
colorh<-rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(256))

ggplot(a, aes(x=log2FC, y=-log10(q), fill=log2FC))+
  geom_point(pch=21)+
  theme_bw()+
  labs(x="log2FC(CRC / Adjacent)", y="-log10(Adjusted p-value)")+
  scale_fill_gradientn(colours = colorh)+
  geom_hline(yintercept = -log10(0.05), col="red", linetype="dashed")+
  ggtitle(label = "metaFlux pathway score")

### Metabolic Pathways analysis -----
data("human_blood")

exp <- TPM
scores<-calculate_reaction_score(exp)
flux<-compute_flux(mras=scores, mmedium=human_blood)

cbrt <- function(x) {
  sign(x) * abs(x)^(1/3)
}

flux=cbrt(flux)

pathway<-unique(unlist(human_gem$SUBSYSTEM))
pathway_score<-list()

for (i in pathway){
  path=i
  activity_score<-c()
  for (d in 1:ncol(flux)){
    activity_score[d]<-mean(abs(flux[which(unlist(human_gem$SUBSYSTEM)==i),d]))
  } 
  pathway_score[[i]]<-activity_score
}

all_pathway_score<-as.data.frame(do.call(rbind,pathway_score))

names(all_pathway_score) <- names(exp)

write.table(all_pathway_score, "0_Input_Tissue_metabolic_path.tsv", sep="\t", quote=F)

### Volcano plot
path <- read.delim("0_Input_Tissue_metabolic_path.tsv", row.names=1)

dat1 <- filter(dat, Sample %in% c("CRC", "Adjacent") & dat$ID %in% names(path)) %>%
  select("ID", agevec, Sample)

path <- path[, dat1$ID]

dat1 <- column_to_rownames(dat1, var = "ID")

CRCid <- row.names(filter(dat1, Sample == "CRC"))
ADid <- row.names(filter(dat1, Sample == "Adjacent"))

CRC <- path[, CRCid]
AD <- path[, ADid]

res <- matrix(nrow=nrow(path), ncol=2)
row.names(res) <- row.names(path)
colnames(res) <- c("log2FC", "p")

for(i in 1:nrow(res)){
  
  res[i,1] <- log2(rowMeans(CRC[i,])/rowMeans(AD[i,]))
  res[i,2] <- wilcox.test(as.numeric(CRC[i,]), as.numeric(AD[i,]))$p.value
  
}

res <- data.frame(res)

res$q <- p.adjust(res$p, method="BH")

colorh<-rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(256))

ggplot(res, aes(x=log2FC, y=-log10(q), fill=log2FC))+
  geom_point(pch=21)+
  theme_bw()+
  labs(x="log2FC(CRC / Adjacent)", y="-log10(Adjusted p-value)")+
  scale_fill_gradientn(colours = colorh)+
  geom_hline(yintercept = -log10(0.05), col="red", linetype="dashed")+
  ggtitle(label = "metaFlux pathway score")

## Figure 3B -- Volcano ----
CRCid <- row.names(filter(dat1, Sample == "CRC"))
ADid <- row.names(filter(dat1, Sample == "Adjacent"))

CRC <- path[, CRCid]
AD <- path[, ADid]

res <- matrix(nrow=nrow(path), ncol=2)
row.names(res) <- row.names(path)
colnames(res) <- c("log2FC", "p")

for(i in 1:nrow(res)){
  
  res[i,1] <- log2(rowMeans(CRC[i,])/rowMeans(AD[i,]))
  res[i,2] <- wilcox.test(as.numeric(CRC[i,]), as.numeric(AD[i,]))$p.value
  
}

res <- data.frame(res)

res$q <- p.adjust(res$p, method="BH")

colorh<-rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(256))

ggplot(res, aes(x=log2FC, y=-log10(q), fill=log2FC))+
  geom_point(pch=21)+
  theme_bw()+
  labs(x="log2FC(CRC / Adjacent)", y="-log10(Adjusted p-value)")+
  scale_fill_gradientn(colours = colorh)+
  geom_hline(yintercept = -log10(0.05), col="red", linetype="dashed")+
  ggtitle(label = "metaFlux pathway score")

write.table(res, "2.06_Table_Metaflux_diff_analysis.tsv", sep="\t", quote=F)

### Figure 3C -- Analysis gene of interest -----
library("Hmisc")

dat1 <- dplyr::filter(dat, Sample %in% c("Stool")) %>% 
  dplyr::select("Sample", "ID", agevec, "ID adjacent RNA-Seq")

goi <- c("GLO1", "TTPA", "MSR1", "SCARB1", "CD36", "SLC23A2", "GATD3", "AGER", "DDOST", "LGALS3", "TPI1", "HAGH", "AKR1B1", "MLXIPL", "MMP2", "SP1")

geneinfo <- read.delim("Gencode_v40_Gene_information.txt", header=F)

geneinfo <- filter(geneinfo, V3 %in% goi)

Exp <- read.delim("0_Input_Tissue_RNA-Seq_TPM.tsv", row.names=1, check.names = F)

id <- intersect(unique(dat1$`ID adjacent RNA-Seq`), names(Exp))

Exp <- Exp[, id]
dat1 <- filter(dat1, dat1$`ID adjacent RNA-Seq` %in% id)

Exp <- Exp[geneinfo$V1, ]

flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut])}

row.names(Exp) <- geneinfo$V3
row.names(dat1) <- dat1$`ID adjacent RNA-Seq`

dat1 <- dat1[, agevec]

res_p <- rcorr(as.matrix(dat1), as.matrix(t(Exp)), type="spearman")

out <- flattenCorrMatrix(res_p$r, res_p$P)
out <- filter(out, out$row %in% names(dat1) & out$column %in% row.names(Exp))

write.table(out, "2.06_Table_Correlation_Specific_genes_adjacent_stool.tsv", sep="\t", quote=F)

library("ComplexHeatmap")
library("circlize")

coefficient <- out[, c("row", "column", "cor")]
pval <- out[, c("row", "column", "p")]

coefficient$row <- factor(coefficient$row, levels = agevec)
pval$row <- factor(pval$row, levels = agevec)

coefficient_spread <- spread(coefficient, row, cor)

pval_spread <- spread(pval, row, p)

row.names(coefficient_spread) <- coefficient_spread$column
row.names(pval_spread) <- pval_spread$column
coefficient_spread$column <- NULL
pval_spread$column <- NULL

colorh<-rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(256))

agevec <- c("MGO", "GO", "3-DG","MGO/GO", "PB MG-H1", "PB CEL", "PB CML", "MG-H1+CEL/MGO", "CML/GO", "MG-H1+CEL/CML", "CEL/MG-H1")

pval_spread <- t(pval_spread[, agevec])

coefficient_spread <- t(coefficient_spread[, agevec])

Heatmap(coefficient_spread,
        row_names_gp = gpar(fontsize = 8),
        column_names_gp = gpar(fontsize = 8),
        column_title_gp = gpar(fontsize = 8),
        cluster_columns = F,
        cluster_rows = F,
        name = "Rho",
        border = "black",
        column_names_rot = 45,
        col=colorRamp2(breaks=c(-1,-0.5,0,0.5,1), 
                       colors=c(colorh[1], colorh[64], "white", colorh[192], colorh[256])),
        heatmap_legend_param = list(border = "black"),
        clustering_method_rows = "ward.D2",
        clustering_method_columns = "ward.D2",
        clustering_distance_columns =  "spearman", 
        clustering_distance_rows =  "spearman",
        cell_fun = function(j,i,x,y,w,h,fill) {
          if(pval_spread[i, j] < 0.001) {
            grid.text("***", x, y, vjust = 0.75)
          } else if(pval_spread[i, j] < 0.01) {
            grid.text("**", x, y, vjust = 0.75)
          } else if(pval_spread[i, j] < 0.05) {
            grid.text("*", x, y, vjust = 0.75)
          }})

### Figure 3D -- Pathway correlation -----
library("ComplexHeatmap")
library("circlize")
library("Hmisc")

path <- read.delim("0_Input_Tissue_metabolic_path.tsv", row.names=1)

dat1 <- filter(dat, Sample %in% c("Tumor", "Adjacent") & dat$ID %in% names(path)) %>%
  select("ID", agevec)

path <- path[, dat1$ID]

dat1 <- column_to_rownames(dat1, var = "ID")

flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut])}

res_p <- rcorr(as.matrix(dat1), as.matrix(t(path)), type="spearman")

out <- flattenCorrMatrix(res_p$r, res_p$P)
out <- filter(out, out$row %in% names(dat1) & out$column %in% row.names(path))

write.table(out, "2.06_Table_Correlation_Specific_pathway.tsv", sep="\t", quote=F, row.names=F)

out$q <- p.adjust(out$p, method="BH")
out$ID <- paste(out$row, out$column, sep="_")

sig <- filter(out, q<0.05)
out <- filter(out, column %in% sig$column)

coefficient <- out[, c("row", "column", "cor")]
pval <- out[, c("row", "column", "p")]

agevec <- c("MGO", "GO", "3-DG","MGO/GO", "PB MG-H1", "PB CEL", "PB CML", "MG-H1+CEL/MGO", "CML/GO", "MG-H1+CEL/CML", "CEL/MG-H1")

coefficient$row <- factor(coefficient$row, levels = agevec)

pval$row <- factor(pval$row, levels = agevec)

coefficient_spread <- spread(coefficient, row, cor)

pval_spread <- spread(pval, row, p)

row.names(coefficient_spread) <- coefficient_spread$column

row.names(pval_spread) <- pval_spread$column

coefficient_spread$column <- NULL
pval_spread$column <- NULL

colorh<-rev(colorRampPalette(RColorBrewer::brewer.pal(11,"RdBu"))(256))

coefficient_spread <- t(coefficient_spread)

pval_spread <- t(pval_spread)

coefficient_spread <- coefficient_spread[as.character(agevec),]
pval_spread <- pval_spread[as.character(agevec),]

Heatmap(coefficient_spread,
        row_names_gp = gpar(fontsize = 6),
        column_names_gp = gpar(fontsize = 4),
        column_title_gp = gpar(fontsize = 6),
        cluster_columns = T,
        cluster_rows = F,
        name = "Rho",
        border = "black",
        column_names_rot = 45,
        col=colorRamp2(breaks=c(-1,-0.5,0,0.5,1), 
                       colors=c(colorh[1], colorh[64], "white", colorh[192], colorh[256])),
        heatmap_legend_param = list(border = "black"),
        clustering_method_rows = "ward.D2",
        clustering_method_columns = "ward.D2",
        clustering_distance_columns =  "spearman", 
        clustering_distance_rows =  "spearman",
        cell_fun = function(j,i,x,y,w,h,fill) {
          if(pval_spread[i, j] < 0.001) {
            grid.text("***", x, y, vjust = 0.75)
          } else if(pval_spread[i, j] < 0.01) {
            grid.text("**", x, y, vjust = 0.75)
          } else if(pval_spread[i, j] < 0.05) {
            grid.text("*", x, y, vjust = 0.75)
          }}
)