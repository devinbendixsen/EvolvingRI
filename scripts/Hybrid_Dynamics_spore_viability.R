#"`-:-.   ,-;"`-:-.   ,-;"`-:-.   ,-;"`-:-. -------------------------- ,-;"`-:-.
# `=`,'=/     `=`,'=/     `=`,'=/     `=`|  "If I have seen further |=/     `=
#   y==/        y==/        y==/        y| it is by standing on the |/
# ,=,-<=`.    ,=,-<=`.    ,=,-<=`.    ,=,|    shoulders of Giants"  |=`.    ,=
#,-'-'   `-=_,-'-'   `-=_,-'-'   `-=_,-'-'  -------------------------- `-=_,-'-'
# ==============================================================================
# PROJECT: Hybrid Dynamics
# PURPOSE: to plot the spore viability data
# DATE: 20 Dec 2023
# Devin P. Bendixsen, PhD
# Staff Bioinformatician | University of Edinburgh
# MRC Human Genetics Unit | Institute of Genetics and Cancer
# ==============================================================================

# ==============================================================================
# LOAD NEEDED PACKAGES
# ==============================================================================
library(readxl)

# ==============================================================================
# LOAD AND FORMAT DATA
# ==============================================================================
data <- read_excel('data/spore_viability/Tetrad dissections spore viability.xlsx',sheet='Table_F1_Parents')


# ==============================================================================
# plot
# ==============================================================================

H1 <- subset(data, (env == "NaCl" | env == "LiAc0.01" | env =='H1F1'))
H1$env <- factor(H1$env , levels=c('NaCl','LiAc0.01','H1F1'))
H1_plot <- ggplot(H1, aes(x=as.factor(env), y=Proportion,fill=env)) + 
  geom_boxplot(colour='black', alpha=1,outlier.shape = NA) +
  geom_jitter(aes(shape=color="black", size=1, alpha=0.9,width=0)) +
  facet_wrap(~factor(gen,c('G100','G200','G300','G400','G500','G700','G1000')), ncol=7)  +
  scale_fill_manual(values=c('#4472c4ff','#9e49e1ff','#808080'),name ='',labels = c("NaCl", "LiAc 0.01", "Hybrid F1"))+
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.grid.major.x = element_line(colour = "lightgrey", linetype = 1, linewidth = 0.5),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y =  element_line(colour = "lightgrey", linetype = 1, linewidth = 0.5),
        panel.grid.minor.y = element_blank(),
        panel.spacing = unit(0.5, "lines"),
        text=element_text(family="EB Garamond")) +
  ylim(0.2,1)

H1_total <- ggplot(H1, aes(x=as.factor(env), y=Proportion,fill=env)) + 
  geom_boxplot(colour='black', alpha=1,outlier.shape = NA) +
  geom_jitter(color="black", size=1, alpha=0.9,width=0) +
  scale_fill_manual(values=c('#4472c4ff','#9e49e1ff','#808080'),name ='',labels = c("NaCl", "LiAc 0.01", "Hybrid F1"))+
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.background = element_rect(fill = "lightgrey"),
        panel.grid.major.x = element_line(colour = "white", linetype = 1, linewidth = 0.5),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y =  element_line(colour = "white", linetype = 1, linewidth = 0.5),
        panel.grid.minor.y = element_blank(),
        panel.spacing = unit(0.5, "lines"),
        text=element_text(family="EB Garamond"),legend.position = "none",
        plot.tag.position=c(0.63,0.94)) + 
  ylab('spore viability') +
  ylim(0.2,1) +
  labs(tag='TOTAL')

H2 <- subset(data, (env == "NaCl" | env == "LiAc0.02" | env =='H2F1' ))
H2$env <- factor(H2$env , levels=c('NaCl','LiAc0.02','H2F1'))
H2_plot <- ggplot(H2, aes(x=as.factor(env), y=Proportion,fill=env)) + 
  geom_boxplot(colour='black', alpha=1,outlier.shape = NA) +
  geom_jitter(color="black", size=1, alpha=0.9,width=0) +
  facet_wrap(~factor(gen,c('G100','G200','G300','G400','G500','G700','G1000')), ncol=7)  +
  scale_fill_manual(values=c('#4472c4ff','#00b050ff','#808080'),name ='',labels = c("NaCl", "LiAc 0.02", "Hybrid F1"))+
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.grid.major.x = element_line(colour = "lightgrey", linetype = 1, linewidth = 0.5),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y =  element_line(colour = "lightgrey", linetype = 1, linewidth = 0.5),
        panel.grid.minor.y = element_blank(),
        panel.spacing = unit(0.5, "lines"),
        text=element_text(family="EB Garamond")) + 
  ylim(0.2,1)

H2_total <- ggplot(H2, aes(x=as.factor(env), y=Proportion,fill=env)) + 
  geom_boxplot(colour='black', alpha=1,outlier.shape = NA) +
  geom_jitter(color="black", size=1, alpha=0.9,width=0)  +
  scale_fill_manual(values=c('#4472c4ff','#00b050ff','#808080'),name ='',labels = c("NaCl", "LiAc 0.02", "Hybrid F1"))+
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.background = element_rect(fill = "lightgrey"),
        panel.grid.major.x = element_line(colour = "white", linetype = 1, linewidth = 0.5),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y =  element_line(colour = "white", linetype = 1, linewidth = 0.5),
        panel.grid.minor.y = element_blank(),
        panel.spacing = unit(0.5, "lines"),
        text=element_text(family="EB Garamond"),legend.position = "none",
        plot.tag.position=c(0.63,0.94)) + 
  ylab('spore viability') +
  ylim(0.2,1) +
  labs(tag='TOTAL')

H3 <- subset(data, (env == "NaCl" | env == "Ethanol" | env =='H3F1'))
H3$env <- factor(H3$env , levels=c('NaCl','Ethanol','H3F1'))
H3_plot <- ggplot(H3, aes(x=as.factor(env), y=Proportion,fill=env)) + 
  geom_boxplot(colour='black', alpha=1,outlier.shape = NA) +
  geom_jitter(color="black", size=1, alpha=0.9,width=0) +
  facet_wrap(~factor(gen,c('G100','G200','G300','G400','G500','G700','G1000')), ncol=7)  +
  scale_fill_manual(values=c('#4472c4ff','#ff0000ff','#808080'),name ='',labels = c("NaCl", "Ethanol", "Hybrid F1"))+
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        panel.background = element_rect(fill = "white"),
        panel.grid.major.x = element_line(colour = "lightgrey", linetype = 1, linewidth = 0.5),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y =  element_line(colour = "lightgrey", linetype = 1, linewidth = 0.5),
        panel.grid.minor.y = element_blank(),
        panel.spacing = unit(0.5, "lines"),
        text=element_text(family="EB Garamond")) + 
  ylim(0.2,1)

H3_total <- ggplot(H3, aes(x=as.factor(env), y=Proportion,fill=env)) + 
  geom_boxplot(colour='black', alpha=1,outlier.shape = NA) +
  geom_jitter(color="black", size=1, alpha=0.9,width=0) +
  scale_fill_manual(values=c('#4472c4ff','#ff0000ff','#808080'),name ='',labels = c("NaCl", "Ethanol", "Hybrid F1"))+
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.background = element_rect(fill = "lightgrey"),
        panel.grid.major.x = element_line(colour = "white", linetype = 1, linewidth = 0.5),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y =  element_line(colour = "white", linetype = 1, linewidth = 0.5),
        panel.grid.minor.y = element_blank(),
        panel.spacing = unit(0., "lines"),
        text=element_text(family="EB Garamond"),legend.position = "none",
        plot.tag.position=c(0.63,0.94)) + 
  ylab('spore viability') +
  ylim(0.2,1) +
  labs(tag='TOTAL')

H1_total + H1_plot + H2_total + H2_plot + H3_total + H3_plot +
  plot_layout(ncol=2,widths = c(1.3,5))

ggsave('figures/Hybrid_Dynamics_spore_viability.pdf',width=8,height=5,dpi = 900)

