#"`-:-.   ,-;"`-:-.   ,-;"`-:-.   ,-;"`-:-. -------------------------- ,-;"`-:-.
# `=`,'=/     `=`,'=/     `=`,'=/     `=`|  "If I have seen further |=/     `=
#   y==/        y==/        y==/        y| it is by standing on the |/
# ,=,-<=`.    ,=,-<=`.    ,=,-<=`.    ,=,|    shoulders of Giants"  |=`.    ,=
#,-'-'   `-=_,-'-'   `-=_,-'-'   `-=_,-'-'  -------------------------- `-=_,-'-'
# ==============================================================================
# PROJECT: Hybrid_Dynamics
# PURPOSE: 
# DATE:
# Devin P. Bendixsen, PhD
# Staff Bioinformatician | University of Edinburgh
# MRC Human Genetics Unit | Institute of Genetics and Cancer
# ==============================================================================

# ==============================================================================
# LOAD NEEDED PACKAGES
# ==============================================================================
library(tidyr)
library(dplyr)
library(ggplot2)
library(patchwork)
library(showtext)
font_add_google("EB Garamond")
showtext_auto()

# ==============================================================================
# IMPORT SOMATIC VCFS AND EXTRACT NEEDED DATA
# ==============================================================================
fileNames <- Sys.glob("results/FREEC/*_CNVs.p.value.txt") # identify all samples
cnv_data <- data.frame()
for (fileName in fileNames) {
  print(fileName)
  cnv_sample <- read.table(file=fileName,header=FALSE,skip = 1,col.names=c('chr','start','end','CN','type','Wilcoxon','KS'))
  cnv_sample$length <- (cnv_sample$end - cnv_sample$start)
  cnv_sample <- subset(cnv_sample,Wilcoxon<0.05 & KS <0.05)
  
  cnv <- split(cnv_sample,cnv_sample$type)
  amp = sum(cnv$gain$length)
  del = sum(cnv$loss$length)
  
  sample=sub('results/FREEC/','',fileName)
  sample=sub('.mapped.sort.picard.bam_CNVs.p.value.txt','',sample)
  env=strsplit(sample,split="_")[[1]][1]
  gen=strsplit(sample,split="_")[[1]][2]
  rep=strsplit(sample,split="_")[[1]][3]
  
  cnv_summary <- data.frame(row.names=c(sample))
  cnv_summary$sample = sample
  cnv_summary$env = env
  cnv_summary$gen = gen
  cnv_summary$rep = rep
  cnv_summary$amp = amp
  cnv_summary$del = del
  cnv_summary$total = amp+del
  cnv_data <- rbind(cnv_data,cnv_summary)
}
write.table(cnv_data,file='data/Hybrid_Dynamics_CNV_data.txt')
# ==============================================================================
# Environment Colors
# ==============================================================================
NaCl = '#4472c4ff'
LiAc0.01 = '#9e49e1ff'
LiAc0.02 = '#00b050ff'
Ethanol = '#ff0000ff'
F1 = '#808080'
F2 = '#d4a373'
edgewidth=0.4
markersize=2
width=0
alpha=0.7
# ==============================================================================
# 
# ==============================================================================
genome_size=12071326
cnv_data$amp <- (cnv_data$amp/genome_size)*100
cnv_data$del <- (cnv_data$del/genome_size)*100

founder <- subset(cnv_data, gen %in% "Founder")

H1 <- subset(cnv_data, (env == "NaCl" | env == "LiAc0.01" | env =='H1F1' | env =='H1F2'))
H1$env <- factor(H1$env , levels=c('NaCl','LiAc0.01','H1F1','H1F2'))
H1_plot <- ggplot(H1) + 
  geom_hline(yintercept=(mean(founder$amp)),color='black',linetype='dashed',linewidth=0.5)+
  geom_hline(yintercept=(-mean(founder$del)),color='black',linetype='dashed',linewidth=0.5)+
  geom_hline(yintercept=(0),color='black',linetype='solid',linewidth=0.7)+
  stat_summary(aes(x=as.factor(env),y=amp,color=env),fun=median, geom="point", shape=95, size=9,stroke=edgewidth,alpha=1)+
  geom_jitter(aes(x=as.factor(env),y=amp,fill=env),shape=21, colour='black',size=markersize, alpha=alpha,width=width,stroke=edgewidth) +
  stat_summary(aes(x=as.factor(env),y=-(del),color=env),fun=median, geom="point", shape=95, size=9,stroke=edgewidth,alpha=1)+
  geom_jitter(aes(x=as.factor(env),y=-(del),fill=env),shape=21,colour='black', size=markersize, alpha=alpha, width=width,stroke=edgewidth) +
  facet_wrap(~factor(gen,c('G100','G200','G300','G400','G500','G700','G1000')), ncol=7)+
  scale_fill_manual(values=c(NaCl,LiAc0.01,F1,F2),name ='',labels = c("NaCl", "LiAc0.01", "Hybrid F1",'Hybrid F2'))+
  scale_color_manual(values=c(NaCl,LiAc0.01,F1,F2),name ='',labels = c("NaCl", "LiAc0.01", "Hybrid F1",'Hybrid F2'))+
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
        text=element_text(family="EB Garamond"),
        legend.spacing.y = unit(-0.35, 'cm'))+
  guides(fill = guide_legend(byrow = TRUE))
H1_total <- ggplot(H1) + 
  geom_hline(yintercept=(mean(founder$amp)),color='black',linetype='dashed',linewidth=0.5)+
  geom_hline(yintercept=(-mean(founder$del)),color='black',linetype='dashed',linewidth=0.5)+
  geom_hline(yintercept=(0),color='black',linetype='solid',linewidth=0.7)+
  stat_summary(aes(x=as.factor(env),y=amp,color=env),fun=median, geom="point", shape=95, size=12,stroke=edgewidth,alpha=1)+
  geom_jitter(aes(x=as.factor(env),y=amp,fill=env),shape=21, colour='black',size=markersize, alpha=alpha,width=width,stroke=edgewidth) +
  stat_summary(aes(x=as.factor(env),y=-(del),color=env),fun=median, geom="point", shape=95, size=12,stroke=edgewidth,alpha=1)+
  geom_jitter(aes(x=as.factor(env),y=-(del),fill=env),shape=21,colour='black', size=markersize, alpha=alpha, width=width,stroke=edgewidth) +
  scale_fill_manual(values=c(NaCl,LiAc0.01,F1,F2),name ='',labels = c("NaCl", "LiAc0.01", "Hybrid F1",'Hybrid F2'))+
  scale_color_manual(values=c(NaCl,LiAc0.01,F1,F2),name ='',labels = c("NaCl", "LiAc0.01", "Hybrid F1",'Hybrid F2'))+
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.background = element_rect(fill = "lightgrey"),
        panel.grid.major.x = element_line(colour = "white", linetype = 1, linewidth = 0.5),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(colour = "white", linetype = 1, linewidth = 0.5),
        panel.grid.minor.y = element_blank(),
        panel.spacing = unit(0.5, "lines"), 
        text=element_text(family="EB Garamond"),
        legend.position = "none",
        plot.tag.position=c(0.62,0.94),
        plot.tag = element_text(hjust = 0.5)) + 
  ylab('genome altered (%)') +
  labs(tag='TOTAL')

H2 <- subset(cnv_data, (env == "NaCl" | env == "LiAc0.02" | env =='H2F1' | env =='H2F2'))
H2$env <- factor(H2$env , levels=c('NaCl','LiAc0.02','H2F1','H2F2'))
H2_plot <- ggplot(H2) + 
  geom_hline(yintercept=(mean(founder$amp)),color='black',linetype='dashed',linewidth=0.5)+
  geom_hline(yintercept=(-mean(founder$del)),color='black',linetype='dashed',linewidth=0.5)+
  geom_hline(yintercept=(0),color='black',linetype='solid',linewidth=0.7)+
  stat_summary(aes(x=as.factor(env),y=amp,color=env),fun=median, geom="point", shape=95, size=9,stroke=edgewidth,alpha=1)+
  geom_jitter(aes(x=as.factor(env),y=amp,fill=env),shape=21, colour='black',size=markersize, alpha=alpha, width=width,stroke=edgewidth) +
  stat_summary(aes(x=as.factor(env),y=-(del),color=env),fun=median, geom="point", shape=95, size=9,stroke=edgewidth,alpha=1)+
  geom_jitter(aes(x=as.factor(env),y=-(del),fill=env),shape=21,colour='black', size=markersize, alpha=alpha, width=width, stroke=edgewidth) +
  facet_wrap(~factor(gen,c('G100','G200','G300','G400','G500','G700','G1000')), ncol=7)+
  scale_fill_manual(values=c(NaCl,LiAc0.02,F1,F2),name ='',labels = c("NaCl", "LiAc 0.02", "Hybrid F1",'Hybrid F2'))+
  scale_color_manual(values=c(NaCl,LiAc0.02,F1,F2),name ='',labels = c("NaCl", "LiAc 0.02", "Hybrid F1",'Hybrid F2'))+
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
        text=element_text(family="EB Garamond"),
        legend.spacing.y = unit(-0.35, 'cm'))+
  guides(fill = guide_legend(byrow = TRUE))
H2_total <- ggplot(H2) + 
  geom_hline(yintercept=(mean(founder$amp)),color='black',linetype='dashed',linewidth=0.5)+
  geom_hline(yintercept=(-mean(founder$del)),color='black',linetype='dashed',linewidth=0.5)+
  geom_hline(yintercept=(0),color='black',linetype='solid',linewidth=0.7)+
  stat_summary(aes(x=as.factor(env),y=amp,color=env),fun=median, geom="point", shape=95, size=12,stroke=edgewidth,alpha=1)+
  geom_jitter(aes(x=as.factor(env),y=amp,fill=env),shape=21, colour='black',size=markersize, alpha=alpha, width=width,stroke=edgewidth) +
  stat_summary(aes(x=as.factor(env),y=-(del),color=env),fun=median, geom="point", shape=95, size=12,stroke=edgewidth,alpha=1)+
  geom_jitter(aes(x=as.factor(env),y=-(del),fill=env),shape=21,colour='black', size=markersize, alpha=alpha, width=width,stroke=edgewidth) +
  scale_fill_manual(values=c(NaCl,LiAc0.02,F1,F2),name ='',labels = c("NaCl", "LiAc 0.02", "Hybrid F1",'Hybrid F2'))+
  scale_color_manual(values=c(NaCl,LiAc0.02,F1,F2),name ='',labels = c("NaCl", "LiAc 0.02", "Hybrid F1",'Hybrid F2'))+
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.background = element_rect(fill = "lightgrey"),
        panel.grid.major.x = element_line(colour = "white", linetype = 1, linewidth = 0.5),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(colour = "white", linetype = 1, linewidth = 0.5),
        panel.grid.minor.y = element_blank(),
        panel.spacing = unit(0.5, "lines"), 
        text=element_text(family="EB Garamond"),
        legend.position = "none",
        plot.tag.position=c(0.62,0.94),
        plot.tag = element_text(hjust = 0.5)) + 
  ylab('genome altered (%)') +
  labs(tag='TOTAL')

H3 <- subset(cnv_data, (env == "NaCl" | env == "Ethanol" | env =='H3F1' | env =='H3F2'))
H3$env <- factor(H3$env , levels=c('NaCl','Ethanol','H3F1','H3F2'))
H3_plot <- ggplot(H3) + 
  geom_hline(yintercept=(mean(founder$amp)),color='black',linetype='dashed',linewidth=0.5)+
  geom_hline(yintercept=(-mean(founder$del)),color='black',linetype='dashed',linewidth=0.5)+
  geom_hline(yintercept=(0),color='black',linetype='solid',linewidth=0.7)+
  stat_summary(aes(x=as.factor(env),y=amp,color=env),fun=median, geom="point", shape=95, size=9,stroke=edgewidth,alpha=1)+
  geom_jitter(aes(x=as.factor(env),y=amp,fill=env),shape=21, colour='black',size=markersize, alpha=alpha, width=width,stroke=edgewidth) +
  stat_summary(aes(x=as.factor(env),y=-(del),color=env),fun=median, geom="point", shape=95, size=9,stroke=edgewidth,alpha=1)+
  geom_jitter(aes(x=as.factor(env),y=-(del),fill=env),shape=21,colour='black', size=markersize, alpha=alpha, width=width, stroke=edgewidth) +
  facet_wrap(~factor(gen,c('G100','G200','G300','G400','G500','G700','G1000')), ncol=7)+
  scale_fill_manual(values=c(NaCl,Ethanol,F1,F2),name ='',labels = c("NaCl", "Ethanol", "Hybrid F1",'Hybrid F2'))+
  scale_color_manual(values=c(NaCl,Ethanol,F1,F2),name ='',labels = c("NaCl", "Ethanol", "Hybrid F1",'Hybrid F2'))+
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
        text=element_text(family="EB Garamond"),
        legend.spacing.y = unit(-0.35, 'cm'))+
  guides(fill = guide_legend(byrow = TRUE))

H3_total <- ggplot(H3) + 
  geom_hline(yintercept=(mean(founder$amp)),color='black',linetype='dashed',linewidth=0.5)+
  geom_hline(yintercept=(-mean(founder$del)),color='black',linetype='dashed',linewidth=0.5)+
  geom_hline(yintercept=(0),color='black',linetype='solid',linewidth=0.7)+
  stat_summary(aes(x=as.factor(env),y=amp,color=env),fun=median, geom="point", shape=95, size=12,stroke=edgewidth,alpha=1)+
  geom_jitter(aes(x=as.factor(env),y=amp,fill=env),shape=21, colour='black',size=markersize, alpha=alpha,width=width,stroke=edgewidth) +
  stat_summary(aes(x=as.factor(env),y=-(del),color=env),fun=median, geom="point", shape=95, size=12,stroke=edgewidth,alpha=1)+
  geom_jitter(aes(x=as.factor(env),y=-(del),fill=env),shape=21,colour='black', size=markersize, alpha=alpha,width=width,stroke=edgewidth) +
  scale_fill_manual(values=c(NaCl,Ethanol,F1,F2),name ='',labels = c("NaCl", "Ethanol", "Hybrid F1",'Hybrid F2'))+
  scale_color_manual(values=c(NaCl,Ethanol,F1,F2),name ='',labels = c("NaCl", "Ethanol", "Hybrid F1",'Hybrid F2'))+
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.background = element_rect(fill = "lightgrey"),
        panel.grid.major.x = element_line(colour = "white", linetype = 1, linewidth = 0.5),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(colour = "white", linetype = 1, linewidth = 0.5),
        panel.grid.minor.y = element_blank(),
        panel.spacing = unit(0.5, "lines"), 
        text=element_text(family="EB Garamond"),
        legend.position = "none",
        plot.tag.position=c(0.62,0.94),
        plot.tag = element_text(hjust = 0.5)) + 
  ylab('genome altered (%)') +
  labs(tag='TOTAL')

H1_total + H1_plot + H2_total + H2_plot + H3_total + H3_plot + 
  plot_layout(ncol=2,widths = c(1.3,8))

ggsave('figures/Hybrid_Dynamics_CNV_burden.pdf',width=11,height=5,dpi = 900)

# ==============================================================================
# 
# ==============================================================================
cnv_data$total <- (cnv_data$total/genome_size)*100
founder <- subset(cnv_data, gen %in% "Founder")

H1 <- subset(cnv_data, (env == "NaCl" | env == "LiAc0.01" | env =='H1F1' | env =='H1F2'))
H1$env <- factor(H1$env , levels=c('NaCl','LiAc0.01','H1F1','H1F2'))
H1_plot <- ggplot(H1) + 
  geom_hline(yintercept=(mean(founder$total)),color='black',linetype='dashed',linewidth=0.5)+
  stat_summary(aes(x=as.factor(env),y=total,color=env),fun=median, geom="point", shape=95, size=9,stroke=edgewidth,alpha=1)+
  geom_jitter(aes(x=as.factor(env),y=total,fill=env),shape=21, colour='black',size=markersize, alpha=alpha,width=width,stroke=edgewidth) +
  facet_wrap(~factor(gen,c('G100','G200','G300','G400','G500','G700','G1000')), ncol=7)+
  scale_fill_manual(values=c(NaCl,LiAc0.01,F1,F2),name ='',labels = c("NaCl", "LiAc0.01", "Hybrid F1",'Hybrid F2'))+
  scale_color_manual(values=c(NaCl,LiAc0.01,F1,F2),name ='',labels = c("NaCl", "LiAc0.01", "Hybrid F1",'Hybrid F2'))+
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
        text=element_text(family="EB Garamond"),
        legend.spacing.y = unit(-0.35, 'cm'))+
  guides(fill = guide_legend(byrow = TRUE))

H1_total <- ggplot(H1) + 
  geom_hline(yintercept=(mean(founder$total)),color='black',linetype='dashed',linewidth=0.5)+
  stat_summary(aes(x=as.factor(env),y=total,color=env),fun=median, geom="point", shape=95, size=12,stroke=edgewidth,alpha=1)+
  geom_jitter(aes(x=as.factor(env),y=total,fill=env),shape=21, colour='black',size=markersize, alpha=alpha,width=width,stroke=edgewidth) +
  scale_fill_manual(values=c(NaCl,LiAc0.01,F1,F2),name ='',labels = c("NaCl", "LiAc0.01", "Hybrid F1",'Hybrid F2'))+
  scale_color_manual(values=c(NaCl,LiAc0.01,F1,F2),name ='',labels = c("NaCl", "LiAc0.01", "Hybrid F1",'Hybrid F2'))+
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.background = element_rect(fill = "lightgrey"),
        panel.grid.major.x = element_line(colour = "white", linetype = 1, linewidth = 0.5),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(colour = "white", linetype = 1, linewidth = 0.5),
        panel.grid.minor.y = element_blank(),
        panel.spacing = unit(0.5, "lines"), 
        text=element_text(family="EB Garamond"),
        legend.position = "none",
        plot.tag.position=c(0.62,0.94),
        plot.tag = element_text(hjust = 0.5)) + 
  ylab('genome altered (%)') +
  labs(tag='TOTAL')

H2 <- subset(cnv_data, (env == "NaCl" | env == "LiAc0.02" | env =='H2F1' | env =='H2F2'))
H2$env <- factor(H2$env , levels=c('NaCl','LiAc0.02','H2F1','H2F2'))
H2_plot <- ggplot(H2) + 
  geom_hline(yintercept=(mean(founder$total)),color='black',linetype='dashed',linewidth=0.5)+
  stat_summary(aes(x=as.factor(env),y=total,color=env),fun=median, geom="point", shape=95, size=9,stroke=edgewidth,alpha=1)+
  geom_jitter(aes(x=as.factor(env),y=total,fill=env),shape=21, colour='black',size=markersize, alpha=alpha,width=width,stroke=edgewidth) +
  facet_wrap(~factor(gen,c('G100','G200','G300','G400','G500','G700','G1000')), ncol=7)+
  scale_fill_manual(values=c(NaCl,LiAc0.02,F1,F2),name ='',labels = c("NaCl", "LiAc 0.02", "Hybrid F1",'Hybrid F2'))+
  scale_color_manual(values=c(NaCl,LiAc0.02,F1,F2),name ='',labels = c("NaCl", "LiAc 0.02", "Hybrid F1",'Hybrid F2'))+
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
        text=element_text(family="EB Garamond"),
        legend.spacing.y = unit(-0.35, 'cm'))+
  guides(fill = guide_legend(byrow = TRUE))

H2_total <- ggplot(H2) + 
  geom_hline(yintercept=(mean(founder$total)),color='black',linetype='dashed',linewidth=0.5)+
  stat_summary(aes(x=as.factor(env),y=total,color=env),fun=median, geom="point", shape=95, size=12,stroke=edgewidth,alpha=1)+
  geom_jitter(aes(x=as.factor(env),y=total,fill=env),shape=21, colour='black',size=markersize, alpha=alpha,width=width,stroke=edgewidth) +
  scale_fill_manual(values=c(NaCl,LiAc0.02,F1,F2),name ='',labels = c("NaCl", "LiAc 0.02", "Hybrid F1",'Hybrid F2'))+
  scale_color_manual(values=c(NaCl,LiAc0.02,F1,F2),name ='',labels = c("NaCl", "LiAc 0.02", "Hybrid F1",'Hybrid F2'))+
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.background = element_rect(fill = "lightgrey"),
        panel.grid.major.x = element_line(colour = "white", linetype = 1, linewidth = 0.5),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(colour = "white", linetype = 1, linewidth = 0.5),
        panel.grid.minor.y = element_blank(),
        panel.spacing = unit(0.5, "lines"), 
        text=element_text(family="EB Garamond"),
        legend.position = "none",
        plot.tag.position=c(0.62,0.94),
        plot.tag = element_text(hjust = 0.5)) + 
  ylab('genome altered (%)') +
  labs(tag='TOTAL')

H3 <- subset(cnv_data, (env == "NaCl" | env == "Ethanol" | env =='H3F1' | env =='H3F2'))
H3$env <- factor(H3$env , levels=c('NaCl','Ethanol','H3F1','H3F2'))
H3_plot <- ggplot(H3) + 
  geom_hline(yintercept=(mean(founder$total)),color='black',linetype='dashed',linewidth=0.5)+
  stat_summary(aes(x=as.factor(env),y=total,color=env),fun=median, geom="point", shape=95, size=9,stroke=edgewidth,alpha=1)+
  geom_jitter(aes(x=as.factor(env),y=total,fill=env),shape=21, colour='black',size=markersize, alpha=alpha,width=width,stroke=edgewidth) +
  facet_wrap(~factor(gen,c('G100','G200','G300','G400','G500','G700','G1000')), ncol=7)+
  scale_fill_manual(values=c(NaCl,Ethanol,F1,F2),name ='',labels = c("NaCl", "Ethanol", "Hybrid F1",'Hybrid F2'))+
  scale_color_manual(values=c(NaCl,Ethanol,F1,F2),name ='',labels = c("NaCl", "Ethanol", "Hybrid F1",'Hybrid F2'))+
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
        text=element_text(family="EB Garamond"),
        legend.spacing.y = unit(-0.35, 'cm'))+
  guides(fill = guide_legend(byrow = TRUE))

H3_total <- ggplot(H3) + 
  geom_hline(yintercept=(mean(founder$total)),color='black',linetype='dashed',linewidth=0.5)+
  stat_summary(aes(x=as.factor(env),y=total,color=env),fun=median, geom="point", shape=95, size=12,stroke=edgewidth,alpha=1)+
  geom_jitter(aes(x=as.factor(env),y=total,fill=env),shape=21, colour='black',size=markersize, alpha=alpha,width=width,stroke=edgewidth) +
  scale_fill_manual(values=c(NaCl,Ethanol,F1,F2),name ='',labels = c("NaCl", "Ethanol", "Hybrid F1",'Hybrid F2'))+
  scale_color_manual(values=c(NaCl,Ethanol,F1,F2),name ='',labels = c("NaCl", "Ethanol", "Hybrid F1",'Hybrid F2'))+
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.background = element_rect(fill = "lightgrey"),
        panel.grid.major.x = element_line(colour = "white", linetype = 1, linewidth = 0.5),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(colour = "white", linetype = 1, linewidth = 0.5),
        panel.grid.minor.y = element_blank(),
        panel.spacing = unit(0.5, "lines"), 
        text=element_text(family="EB Garamond"),
        legend.position = "none",
        plot.tag.position=c(0.62,0.94),
        plot.tag = element_text(hjust = 0.5)) + 
  ylab('genome altered (%)') +
  labs(tag='TOTAL')

H1_total + H1_plot + H2_total + H2_plot + H3_total + H3_plot + 
  plot_layout(ncol=2,widths = c(1.3,8))

ggsave('figures/Hybrid_Dynamics_GIN.pdf',width=11,height=5,dpi = 900)

# ==============================================================================
# ANOVA of CNV 
# ==============================================================================
res_aov <- aov(total ~ env,
               data = H1)
summary(res_aov)
post_test <- glht(res_aov,
                  linfct = mcp(env = "Tukey"))
summary(post_test)
plot(post_test)

res_aov <- aov(total ~ env,
               data = H2)
summary(res_aov)
post_test <- glht(res_aov,
                  linfct = mcp(env = "Tukey"))
summary(post_test)
plot(post_test)

res_aov <- aov(total ~ env,
               data = H3)
summary(res_aov)
post_test <- glht(res_aov,
                  linfct = mcp(env = "Tukey"))
summary(post_test)
plot(post_test)


# ==============================================================================
# 
# ==============================================================================

H1 <- subset(H1,env %in% c('H1F1','H1F2'))
H1$gen_rep <- paste0(H1$gen,'_',H1$rep)
p_H1 <- ggplot(H1, aes(x = env, y = total,fill=env))+
  geom_boxplot()+
  scale_fill_manual(values=c('#808080','#d4a373'))+
  geom_line(aes(group=gen_rep))+
  geom_point(shape=21,fill='white')+
  theme(panel.background = element_rect(fill = "white"),
        text=element_text(family="EB Garamond"),
        plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
        axis.title.x=element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(colour = "lightgrey", linetype = 1, linewidth = 0.5),
        panel.grid.minor.y = element_blank(),
        legend.position='none')+
  scale_x_discrete(labels=c('F1','F2'))+
  ylab('genome altered (%)')+
  ylim(3,20.5)
t.test(data=H1,total ~ env)

H2 <- subset(H2,env %in% c('H2F1','H2F2'))
H2$gen_rep <- paste0(H2$gen,'_',H2$rep)
p_H2 <- ggplot(H2, aes(x = env, y = total,fill=env))+
  geom_boxplot()+
  scale_fill_manual(values=c('#808080','#d4a373'))+
  geom_line(aes(group=gen_rep))+
  geom_point(shape=21,fill='white')+
  theme(panel.background = element_rect(fill = "white"),
        text=element_text(family="EB Garamond"),
        plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(colour = "lightgrey", linetype = 1, linewidth = 0.5),
        panel.grid.minor.y = element_blank(),
        legend.position='none')+
  scale_x_discrete(labels=c('F1','F2'))+
  ylab('genome altered (%)')+
  ylim(3,20.5)
t.test(data=H2,total ~ env)

H3 <- subset(H3,env %in% c('H3F1','H3F2'))
H3$gen_rep <- paste0(H3$gen,'_',H3$rep)
p_H3 <- ggplot(H3, aes(x = env, y = total,fill=env))+
  geom_boxplot()+
  scale_fill_manual(values=c('#808080','#d4a373'))+
  geom_line(aes(group=gen_rep))+
  geom_point(shape=21,fill='white')+
  theme(panel.background = element_rect(fill = "white"),
        text=element_text(family="EB Garamond"),
        plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
        axis.title.x=element_blank(),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(colour = "lightgrey", linetype = 1, linewidth = 0.5),
        panel.grid.minor.y = element_blank(),
        legend.position='none')+
  scale_x_discrete(labels=c('F1','F2'))+
  ylab('genome altered (%)')+
  ylim(3,20.5)
t.test(data=H3,total ~ env)

p_H1 + p_H2 + p_H3
ggsave('figures/Hybrid_Dynamics_F1-F2.pdf',width=6,height=2,dpi = 900)

# ==============================================================================
# 
# ==============================================================================
ggplot(c)


