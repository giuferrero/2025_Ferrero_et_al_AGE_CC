library("tidyverse")
library("ggpubr")
library("RColorBrewer")
library("gtsummary")

path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(path)

agevec <- c("MGO", "GO", "3-DG", "PB CML", "PB CEL", "PB MG-H1", "MGO/GO", "MG-H1+CEL/CML", "MG-H1+CEL/MGO")

### Summary Table 1 ----
dat <- read.delim("0_Input_data.txt", check.names = F)
dat$OS.status <- ifelse(dat$OS.status == 0, "Alive", "Deceased")

table1 <- 
  dat %>%
  tbl_summary(include = c(Age, Sex, BMI, Smoking, Diabetes, Hypertriglyceridemia, Hypertension, Drugs, Familiarity, CMS, Stage, Grade, APC, KRAS, TP53, MSI), type = all_dichotomous() ~ "categorical",
              statistic = list(all_continuous() ~ "{mean} ({sd})"))

table1 <- tbl_split(table1, CMS)

table1 %>%
  as_gt() %>%
  gt::gtsave(filename = "Table1.docx")

### Summary Table 2 ----
dat <- read.delim("0_Input_data.txt", check.names = F)

table2 <- 
  dat %>%
  tbl_summary(include = c("MGO", "GO", "3-DG", "PB CML", "PB CEL", "PB MG-H1", "MGO/GO", "MG-H1+CEL/CML", "CEL/MG-H1", "MG-H1+CEL/MGO",	"CML/GO"), by=Sample, type = all_dichotomous() ~ "categorical",
              statistic = list(all_continuous() ~ "{mean} +/- {sd}")) %>% add_p(test = list(all_continuous() ~ "paired.wilcox.test"), group=Patient)

table2 %>%
  as_gt() %>%
  gt::gtsave(filename = "Table2.docx")

### Figure 1B - Boxplot tissue levels ----
dat$Sample <- factor(dat$Sample, levels=c("Adjacent", "Tumor", "Stool", "Plasma"))

dat1 <- filter(dat, Sample %in% c("Adjacent", "Tumor")) %>% 
  select("Sample", "ID",  "Subject", agevec) %>%
  gather(AGE, Level, -c(Sample, ID, Subject))

dat1$AGE <- factor(dat1$AGE, levels=agevec)

ggplot(dat1, aes(x=Sample, y=log10(Level), col=Sample))+
  geom_boxplot(width=0.4, outlier.alpha = 0, aes(group = Sample))+
  geom_jitter(width = 0, alpha=0.5)+
  geom_line(aes(group = Subject), col="grey", alpha=0.5)+
  facet_wrap(~AGE, nrow=1)+
  theme_bw()+
  labs(x="Sample type",  y="log10(AGE levels)") +
  scale_color_manual(values = c("darkblue", "darkred"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=1))+
  geom_pwc(method.args = list(paired=T))

### Figure 1C - Heatmap Ratios AGE levels ----
library("ComplexHeatmap")
library("circlize")

agevec <- c("MGO", "GO", "3-DG", "PB CML", "PB CEL", "PB MG-H1", "MGO/GO", "MG-H1+CEL/CML", "MG-H1+CEL/MGO")

dat <- filter(dat, Sample == "Ratio")
dat2 <- dplyr::select(dat, "ID", agevec)

pheno <- dplyr::select(dat, "ID", "Sample", "CMS", "Sex", "Age", "BMI", "Smoking", "Diabetes", "Hypertriglyceridemia", "CRIS", "OS.status", "Stage", "Grade", "Colon/Rectum", "APC", "KRAS", "TP53", "MSI", "Mucinous_class")

dat2 <- column_to_rownames(dat2, var = "ID")
pheno <- column_to_rownames(pheno, var = "ID")

colorh<-rev(colorRampPalette(brewer.pal(11,"RdBu"))(256))
class.c<-colorRampPalette(brewer.pal(12,"Paired"))(12)
stage.c<-colorRampPalette(brewer.pal(9,"YlOrRd"))(9)
grade.c<-colorRampPalette(brewer.pal(9,"Blues"))(9)
age.color <- colorRampPalette(brewer.pal(9,"Blues"))(256)

class.color <- c("CRC"=class.c[6], "Adjacent"=class.c[5])
stage.color <- c("I"=stage.c[2], "II"=stage.c[4], "III"=stage.c[6], "IV"=stage.c[8])
grade.color <- c("G1-G2"=grade.c[2], "G2"=grade.c[4], "G3"=grade.c[6])
msi.color <- c("High"="black", "Stable"="white")
mut.color <- c("Mutated"="black", "WT"="white")
sex.color <- c("Female"=class.c[5], "Male"=class.c[1])

cms.color <- c("#E89E33", "#0272AC", "#D079A4", "#009E76", "#8F8F8F", "lightsteelblue")
names(cms.color) <- c("CMS1", "CMS2", "CMS3", "CMS4", "Not assigned", "Healthy")
cris.color <- c("#FF6729", "#AD051E", "#03205B", "#13833F", "#25C9AF", "#8F8F8F", "lightsteelblue")
names(cris.color) <- c("CRIS-A", "CRIS-B", "CRIS-C", "CRIS-D", "CRIS-E", "Healthy")

pheno$CMS <- ifelse(pheno$CMS == "Not Assigned", NA, pheno$CMS)

column_ha = columnAnnotation(
  Sex=pheno$Sex,
  Age=pheno$Age,
  BMI=pheno$BMI,
  Stage=pheno$Stage,
  Grade=pheno$Grade,
  CMS=pheno$CMS,
  MSI=pheno$MSI,
  APC=pheno$APC,
  KRAS=pheno$KRAS,
  TP53=pheno$TP53,
  col = list(Type = class.color,
             Stage = stage.color,
             Grade = grade.color,
             MSI = msi.color,
             APC = mut.color,
             KRAS = mut.color,
             TP53 = mut.color,
             Sex = sex.color, 
             CMS = cms.color),
  border = TRUE,
  simple_anno_size = unit(0.2, "cm"),
  annotation_name_gp = gpar(fontsize =5),
  annotation_legend_param = list(border="black"))

Heatmap(t(dat2),
        top_annotation = column_ha,
        row_names_gp = gpar(fontsize = 6),
        column_names_gp = gpar(fontsize = 6),
        column_title_gp = gpar(fontsize = 4),
        cluster_columns = T,
        cluster_rows = T,
        name = "log2FC",
        border = "black",
        column_names_rot = 45,
        col=colorRamp2(breaks=c(-3,-1.5,0,1.5,3), 
                       colors=c(colorh[1], colorh[64], 
                                "white", colorh[192], colorh[256])),
        heatmap_legend_param = list(border = "black"),
        clustering_method_rows = "ward.D2",
        clustering_method_columns = "ward.D2")

### Figure 1C - Boxplot ratios ----
toplot <- gather(dat2, AGE, Level)

colorh<-rev(colorRampPalette(brewer.pal(11,"RdBu"))(256))

agevec <- c("MGO/GO", "MGO", "MG-H1+CEL/CML", "PB MG-H1", "GO", "3-DG", "PB CEL", "PB CML","MG-H1+CEL/MGO")

toplot$AGE <- factor(toplot$AGE, levels=rev(agevec))

ggplot(toplot, aes(y=AGE, x=Level, col=Level))+
  geom_boxplot(width=0.4, outlier.alpha = 0)+
  geom_jitter(height = 0.1, alpha=0.5)+
  theme_bw()+
  labs(y="AGE",  x="log2(Tumor / Adjacent levels)", col="log2(Ratio)") +
  scale_color_gradientn(colors = colorh, lim=c(-4.3,4.3))+
  theme(axis.text.x = element_text(angle = 0, hjust = 1, vjust=1))+
  geom_vline(xintercept = 0, linetype="dashed", col="red")