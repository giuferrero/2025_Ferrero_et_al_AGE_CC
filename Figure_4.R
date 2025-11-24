library("tidyverse")
library("ggpubr")
library("RColorBrewer")

path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(path)

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