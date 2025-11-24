### Supplementary Figure 1A -- Boxplot Diabetes ----
dat1 <- dat %>% 
  select("Sample", "Diabetes", "ID", agevec, "Subject") %>%
  gather(AGE, Level, -c(Sample, ID, Diabetes, Subject))

dat1$AGE <- factor(dat1$AGE, levels=c(agevec))

dat1$Sample <- factor(dat1$Sample, levels=c("Tumor", "Adjacent", "Stool", "Plasma"))

dat1 <- filter(dat1, Diabetes %in% c("No", "Yes"))

ggplot(dat1, aes(x=Diabetes, y=log10(Level), col=Diabetes))+
  geom_boxplot(width=0.4, outlier.alpha = 0)+
  geom_jitter(width = 0, alpha=0.5)+
  facet_grid(rows = vars(Sample), cols = vars(AGE))+
  theme_bw()+
  labs(x="Sample type",  y="log10(AGE levels)") +
  scale_color_manual(values = c("darkblue", "darkred"))+
  stat_compare_means(label.y = 3,method = "wilcox.test", paired=F, comparisons = list(c("No", "Yes")))

### Supplementary Figure 1B -- Hypertriglyceridemia ----
dat1 <- dat %>% 
  select("Sample", "Hypertriglyceridemia", "ID", agevec, "Subject") %>%
  gather(AGE, Level, -c(Sample, ID, Hypertriglyceridemia, Subject))

dat1$AGE <- factor(dat1$AGE, levels=agevec)

dat1$Sample <- factor(dat1$Sample, levels=c("Tumor", "Adjacent", "Stool", "Plasma"))

dat1 <- filter(dat1, Hypertriglyceridemia %in% c("No", "Yes"))

ggplot(dat1, aes(x=Hypertriglyceridemia, y=log10(Level), col=Hypertriglyceridemia))+
  geom_boxplot(width=0.4, outlier.alpha = 0)+
  geom_jitter(width = 0, alpha=0.5)+
  facet_grid(rows = vars(Sample), cols = vars(AGE))+
  theme_bw()+
  labs(x="Hypertriglyceridemia",  y="log10(AGE levels)") +
  scale_color_manual(values = c("darkblue", "darkred"))+
  stat_compare_means(label.y = 3,method = "wilcox.test", paired=F, comparisons = list(c("No", "Yes")))

### Supplementary Figure 2 -- Boxplot CMS ----
dat1 <- filter(dat, Sample %in% c("Adjacent", "Tumor")) %>%
  dplyr::select("Sample", "ID", agevec, "Subject", "CMS") %>%
  gather(AGE, Level, -c(Sample, ID, Subject, CMS))

dat1$AGE <- factor(dat1$AGE, levels=c("MGO", "GO", "3-DG", "PB CML", "PB CEL", "PB MG-H1", "MGO/GO", "MG-H1+CEL/CML", "CEL/MG-H1","MG-H1+CEL/MGO",	"CML/GO"))

dat1 <- filter(dat1, CMS %in% c("CMS1", "CMS2", "CMS3", "CMS4"))

ggplot(dat1, aes(x=Sample, y=log10(Level), col=Sample))+
  geom_boxplot(width=0.4, outlier.alpha = 0, aes(group = Sample))+
  geom_line(aes(group = Subject), col="grey", alpha=0.5)+
  geom_jitter(width = 0, alpha=0.5)+
  facet_grid(cols = vars(AGE), rows = vars(CMS))+
  theme_bw()+
  labs(x="Sample type",  y="log10(AGE levels)") +
  scale_color_manual(values = c("darkblue", "darkred"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=1))+
  ylim(c(-1, 4))+
  geom_pwc(label = "p={p.format}", label.size = 3,vjust = -0.2)