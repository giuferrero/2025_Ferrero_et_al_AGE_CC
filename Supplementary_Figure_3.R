## Code to generate the plots reported in the Supplementary Figure 3

library("tidyverse")
library("ggpubr")
library("RColorBrewer")

path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(path)

agevec <- c("MGO", "GO", "3-DG", "PB CML", "PB CEL", "PB MG-H1", "MGO/GO", "MG-H1+CEL/CML", "MG-H1+CEL/MGO")

### Supplementary Figure 3A -- Correlation Tumor-Stool ----
dat1 <- filter(dat, Sample %in% c("Stool", "Tumor") & Subject != "V216") %>% 
  select("Sample", "ID", agevec, "Subject", "CMS") %>%
  gather(AGE, Level, -c(Sample, ID, Subject, CMS))

dat1$AGE <- factor(dat1$AGE, levels=agevec)

tissue <- filter(dat1, Sample == "Tumor")
stool <- filter(dat1, Sample == "Stool")
id <- intersect(tissue$Subject, stool$Subject)

tissue <- filter(tissue, Subject %in% id)
stool <- filter(stool, Subject %in% id)
tissue$Stool <- stool$Level

ggplot(tissue, aes(x=log10(Level), y=log10(Stool)))+
  geom_point()+
  geom_smooth(method="lm", alpha=0.5)+
  facet_wrap(~AGE, nrow=3, scale="free")+
  theme_bw()+
  labs(x="log10(Tissue levels)",  y="log10(Stool levels)") +
  stat_cor(method = "spearman")+
  ggtitle("Correlation between Tumor and Stool AGE levels")

### Supplementary Figure 3B -- Correlation Adjacent-Stool ----
dat1 <- filter(dat, Sample %in% c("Stool", "Adjacent") & Subject != "V216") %>% 
  select("Sample", "ID", agevec, "Subject", "CMS") %>%
  gather(AGE, Level, -c(Sample, ID, Subject, CMS))

dat1$AGE <- factor(dat1$AGE, levels=agevec)

tissue <- filter(dat1, Sample == "Adjacent")
stool <- filter(dat1, Sample == "Stool")
id <- intersect(tissue$Subject, stool$Subject)

tissue <- filter(tissue, Subject %in% id)
stool <- filter(stool, Subject %in% id)
tissue$Stool <- stool$Level

ggplot(tissue, aes(x=log10(Level), y=log10(Stool)))+
  geom_point()+
  geom_smooth(method="lm", alpha=0.5)+
  facet_wrap(~AGE, nrow=3, scale="free")+
  theme_bw()+
  labs(x="log10(Tissue levels)",  y="log10(Stool levels)") +
  stat_cor(method = "spearman")+
  ggtitle("Correlation between Adjacent and Stool AGE levels")

### Supplementary Figure 3C-D -- Correlation Crossed ----
Plasma <- filter(dat, Sample %in% c("Plasma"))
Stool <- filter(dat, Sample %in% c("Stool"))
Tumor <- filter(dat, Sample %in% c("Tumor"))
Adj <- filter(dat, Sample %in% c("Adjacent"))

Tumor <- filter(Tumor, Subject %in% Stool$Subject)
Adj <- filter(Adj, Subject %in% Stool$Subject)

dat1 <- data.frame(Pre=Adj$MGO, Post=Stool$`PB MG-H1`)

ggplot(dat1, aes(x=log10(Pre), y=log10(Post)))+
  geom_point()+
  geom_smooth(method="lm", alpha=0.5)+
  theme_bw()+
  labs(x="log10(Adjacent tissue levels)",  y="log10(Stool levels)") +
  stat_cor(method = "spearman")+
  ggtitle("Correlation of PB MG-H1 levels")
