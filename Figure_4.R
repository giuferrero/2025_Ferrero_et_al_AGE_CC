library("tidyverse")
library("ggpubr")
library("RColorBrewer")

path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(path)

### Plasma metabolite correlations ----
library("ComplexHeatmap")
library("circlize")
library("Hmisc")

dat <- read.delim("0_Input_data.txt", check.names = F)
met <- read.delim("0_Input_Plasma_metabolites_reduced.tsv", row.names=1)
sdata <- read.delim("0_Input_metadata_Metabolite.txt", row.names=1)

# Filtering metabolites
met <- met[rowSums(is.na(met) == T) < 1,]

dat1 <- filter(dat, Sample %in% c("Plasma") & dat$ID %in% names(met))

a.pca<-prcomp(t(log(met+1,2)), center=F, scale=F)

autoplot(a.pca, data = sdata, colour = 'Batch',  size=3.5, alpha=0.75) +
  theme_bw() + 
  theme(legend.position="top")

dat1 <- filter(dat, Sample %in% c("Plasma") & dat$ID %in% names(met))

dat2 <- filter(dat, Sample %in% c("Plasma") & dat$Patient %in% dat1$Patient)

met <- met[,dat2$`ID metabolite data`]

dat2 <- dat2 %>%
  dplyr::select("ID", agevec)

dat2 <- column_to_rownames(dat2, var = "ID")

flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut])}

res_p <- rcorr(as.matrix(dat2), as.matrix(t(met)), type="spearman")

out <- flattenCorrMatrix(res_p$r, res_p$P)
out <- filter(out, out$row %in% names(dat2) & out$column %in% row.names(met))

out$q <- p.adjust(out$p, method="BH")

write.table(out, "2.06_Table_Correlation_Specific_Metabolites_vs_AGE_Plasma.tsv", sep="\t", quote=F, row.names=F)

out$ID <- paste(out$row, out$column, sep="_")

sig <- filter(out, p<0.01)
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

### Fecal metagenome correlations ----
library("ComplexHeatmap")
library("circlize")
library("Hmisc")

met <- read.delim("0_Input_Stool_metagenome.tsv", row.names=1)

### Metagenome filer
met <- met[grep("s__", row.names(met)),]
met <- met[grep("t__", row.names(met), invert=T),]

row.names(met) <- str_extract(row.names(met), 's__.*$')
row.names(met) <- sub("s__", "", row.names(met))

prev <- matrixStats::rowMedians(as.matrix(met))

### Diversity
library("vegan")
output<- data.frame(Richness=rep(NA,19),
                  Diversity_Inverse_Simpson=rep(NA,19),
                  Diversity_Shannon=rep(NA,19),
                  Evenness_Simpson=rep(NA,19))

### Diversity
output$Richness <- colSums(met>0)
output$Diversity_Inverse_Simpson <- diversity(met, MARGIN = 2, index = "invsimpson")
output$Diversity_Shannon <- diversity(met, MARGIN = 2, index = "shannon")
output$Evenness_Simpson <- diversity(met, MARGIN = 2, index = "shannon") / log(output$Richness)

row.names(output) <- names(met)

write.table(output, "2.06_Table_Diversity_WMS.tsv", sep="\t", quote=F)

## Order relative abundance as sample data
oral <- read.delim("0_oral_species.txt")
SGB_oral <- met[oral$SGB_oral,]
SGB_oral <- na.omit(SGB_oral)
SGB_oral <- as.data.frame(t(SGB_oral))

## Oral_score
output$Oral_score_abundances <- rowSums(SGB_oral)
output$Oral_score_species <- rowSums(SGB_oral>0)

### Filtering
library("caret")
species_nzv <- nearZeroVar(t(met))
met <- met[-species_nzv, ]
met <- met[-species_nzv, ]
a <- rowSums(met>0)
ids <- a[a>5]
met <- met[names(ids), ]

####
dat1 <- filter(dat, Sample %in% c("Stool") & dat$ID %in% names(met))
dat2 <- filter(dat, Sample %in% c("Tumor") & dat$Patient %in% dat1$Patient) 

met <- met[,dat2$`ID stool sRNA-Seq data`]
output <- met[dat2$`ID stool sRNA-Seq data`,]

dat2 <- dat2 %>%
  dplyr::select("ID", agevec)

dat2 <- column_to_rownames(dat2, var = "ID")

flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut])}

res_p <- rcorr(as.matrix(dat2), as.matrix(t(met)), type="spearman")
res_p <- rcorr(as.matrix(dat2), as.matrix(output), type="spearman")

out <- flattenCorrMatrix(res_p$r, res_p$P)
out <- filter(out, out$row %in% names(dat2) & out$column %in% row.names(met))

out$q <- p.adjust(out$p, method="BH")

write.table(out, "2.06_Table_Correlation_Specific_Microbial_Species_vs_Adjacent.tsv", sep="\t", quote=F, row.names=F)

out$ID <- paste(out$row, out$column, sep="_")

sig <- filter(out, p<0.01)
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

### Figure 4A-B -- Barplot correlations ---
a$Class <- factor(a$Class, c("Tumor","Adjacent", "Stool", "Plasma"))

class.c<-colorRampPalette(brewer.pal(12,"Paired"))(12)

col <- c(class.c[6], class.c[2], class.c[8], class.c[4])

ggplot(a, aes(x=Class, y=n, fill=Class))+
  geom_bar(stat = "identity")+
  theme_bw()+
  ggtitle(label = "Correlations with stool species")+
  labs(y="Signficant correlations", x="Dicarbonyl/AGE source")+
  scale_fill_manual(values = col)

dat$Sample <- factor(dat$Sample, levels=c("Tumor", "Adjacent", "Plasma", "Stool"))

### Figure 4D -- Scatterplot oral score ----
class.c<-colorRampPalette(brewer.pal(12,"Paired"))(12)

col <- c(class.c[6], class.c[2], class.c[8], class.c[4])

ggplot(dat, aes(x=log10(Oral_score_species), y=log10(MO/GO), col=Sample))+
  geom_point()+
  geom_smooth(method="lm")+
  facet_wrap(~Sample)+
  scale_color_manual(values = col)+
  theme_bw()+
  stat_cor(method="spearman")+
  labs(x="log10(Oral species score)",
       y="log10(MGO levels)",
       col="Sample type")
