library("tidyverse")
library("ggpubr")
library("RColorBrewer")

path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(path)

agevec <- c("MGO", "GO", "3-DG", "PB CML", "PB CEL", "PB MG-H1", "MGO/GO", "MG-H1+CEL/CML", "MG-H1+CEL/MGO")

### Figure 2A -- Correlation Tumor-Plasma ----
dat1 <- filter(dat, Sample %in% c("Plasma", "Tumor") & Subject != "V216") %>% 
  select("Sample", "ID", agevec, "Subject", "CMS") %>%
  gather(AGE, Level, -c(Sample, ID, Subject, CMS))

dat1$AGE <- factor(dat1$AGE, levels=agevec)

tissue <- filter(dat1, Sample == "Tumor")
plasma <- filter(dat1, Sample == "Plasma")
id <- intersect(tissue$Subject, plasma$Subject)

tissue <- filter(tissue, Subject %in% id)
plasma <- filter(plasma, Subject %in% id)

tissue$Plasma <- plasma$Level

ggplot(tissue, aes(x=log10(Level), y=log10(Plasma)))+
  geom_point()+
  geom_smooth(method="lm", alpha=0.5)+
  facet_wrap(~AGE, nrow=3, scale="free")+
  theme_bw()+
  labs(x="log10(Tissue levels)",  y="log10(Plasma levels)") +
  stat_cor(method = "spearman")+
  ggtitle("Correlation between Tumor and Plasma AGE levels")

### Figure 2B -- Correlation Adjacent-Plasma ----
dat1 <- filter(dat, Sample %in% c("Plasma", "Adjacent")) %>% 
  select("Sample", "ID", agevec, "Subject", "CMS") %>%
  gather(AGE, Level, -c(Sample, ID, Subject, CMS))

dat1$AGE <- factor(dat1$AGE, levels=agevec)

tissue <- filter(dat1, Sample == "Adjacent")
plasma <- filter(dat1, Sample == "Plasma")
id <- intersect(tissue$Subject, plasma$Subject)

tissue <- filter(tissue, Subject %in% id)
plasma <- filter(plasma, Subject %in% id)
tissue$Plasma <- plasma$Level

ggplot(tissue, aes(x=log10(Level), y=log10(Plasma)))+
  geom_point()+
  geom_smooth(method="lm", alpha=0.5)+
  facet_wrap(~AGE, nrow=3, scale="free")+
  theme_bw()+
  labs(x="log10(Tissue levels)",  y="log10(Plasma levels)") +
  stat_cor(method = "spearman")+
  ggtitle("Correlation between Adjacent and Plasma AGE levels")