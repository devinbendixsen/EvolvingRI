#"`-:-.   ,-;"`-:-.   ,-;"`-:-.   ,-;"`-:-. -------------------------- ,-;"`-:-.
# `=`,'=/     `=`,'=/     `=`,'=/     `=`|  "If I have seen further |=/     `=
#   y==/        y==/        y==/        y| it is by standing on the |/
# ,=,-<=`.    ,=,-<=`.    ,=,-<=`.    ,=,|    shoulders of Giants"  |=`.    ,=
#,-'-'   `-=_,-'-'   `-=_,-'-'   `-=_,-'-'  -------------------------- `-=_,-'-'
# ==============================================================================
# PROJECT:
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
quartz()
library(ggpubr)
library(multcomp)

# ==============================================================================
# READ IN AND SUMMAZE DATA FROM SNPEFF
# ==============================================================================
fileNames <- Sys.glob("results/SNV/vcf_snpeff/stats/*.stats.txt") # identify all samples
snv_data <- data.frame()
for (fileName in fileNames) {
  print(fileName)
  snv_sample_stat <- read.csv(file=fileName,header = FALSE, sep = ",", dec = ".",skip=4)
  
  sample=sub('results/SNV/vcf_snpeff/stats/','',fileName)
  sample=sub('.stats.txt','',sample)
  env=strsplit(sample,split="_")[[1]][1]
  gen=strsplit(sample,split="_")[[1]][2]
  rep=strsplit(sample,split="_")[[1]][3]
  
  snv_sample <- data.frame(row.names=c(sample))
  snv_sample$sample = sample
  snv_sample$env = env
  snv_sample$gen = gen
  snv_sample$rep = rep
  snv_sample$tot_burden = as.numeric(snv_sample_stat$V2[7])
  snv_sample$high = as.numeric(snv_sample_stat$V2[59])
  snv_sample$DEL = as.numeric(snv_sample_stat$V2[52])
  snv_sample$INS = as.numeric(snv_sample_stat$V2[53])
  snv_sample$MIXED = as.numeric(snv_sample_stat$V2[54])
  snv_sample$MNP = as.numeric(snv_sample_stat$V2[55])
  snv_sample$SNP = as.numeric(snv_sample_stat$V2[56])
  snv_data <- rbind(snv_data,snv_sample)
}
write.table(snv_data,file='data/Hybrid_Dynamics_SNV_data.txt')
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

founder <- subset(snv_data, gen %in% "Founder")

H1 <- subset(snv_data, (env == "NaCl" | env == "LiAc0.01" | env =='H1F1' | env =='H1F2'))
H1$env <- factor(H1$env , levels=c('NaCl','LiAc0.01','H1F1','H1F2'))
H1_plot <- ggplot(H1) + 
  geom_hline(yintercept=(mean(founder$tot_burden)),color='black',linetype='dashed',linewidth=0.7)+
  stat_summary(aes(x=as.factor(env),y=tot_burden,color=env),fun=median, geom="point", shape=95, size=9,stroke=edgewidth,alpha=1)+
  geom_jitter(aes(x=as.factor(env),y=tot_burden,fill=env),shape=21, colour='black',size=markersize, alpha=alpha,width=width,stroke=edgewidth) +
  facet_wrap(~factor(gen,c('G100','G200','G300','G400','G500','G700','G1000')), ncol=7)+
  scale_fill_manual(values=c(NaCl,LiAc0.01,F1,F2),name ='',labels = c("NaCl", "LiAc 0.01", "Hybrid F1",'Hybrid F2'))+
  scale_color_manual(values=c(NaCl,LiAc0.01,F1,F2),name ='',labels = c("NaCl", "LiAc 0.01", "Hybrid F1",'Hybrid F2'))+
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
  guides(fill = guide_legend(byrow = TRUE)) +
  ylim(70000, 97000)
H1_total <- ggplot(H1) + 
  geom_hline(yintercept=(mean(founder$tot_burden)),color='black',linetype='dashed',linewidth=0.7)+
  stat_summary(aes(x=as.factor(env),y=tot_burden,color=env),fun=median, geom="point", shape=95, size=12,stroke=edgewidth,alpha=1)+
  geom_jitter(aes(x=as.factor(env),y=tot_burden,fill=env),shape=21, colour='black',size=markersize, alpha=alpha,width=width,stroke=edgewidth) +
  scale_fill_manual(values=c(NaCl,LiAc0.01,F1,F2),name ='',labels = c("NaCl", "LiAc 0.01", "Hybrid F1",'Hybrid F2'))+
  scale_color_manual(values=c(NaCl,LiAc0.01,F1,F2),name ='',labels = c("NaCl", "LiAc 0.01", "Hybrid F1",'Hybrid F2'))+
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
  ylab('SNV mutational load') +
  labs(tag='TOTAL')+
  ylim(70000, 97000)

H2 <- subset(snv_data, (env == "NaCl" | env == "LiAc0.02" | env =='H2F1' | env =='H2F2'))
H2$env <- factor(H2$env , levels=c('NaCl','LiAc0.02','H2F1','H2F2'))
H2_plot <- ggplot(H2) + 
  geom_hline(yintercept=(mean(founder$tot_burden)),color='black',linetype='dashed',linewidth=0.7)+
  stat_summary(aes(x=as.factor(env),y=tot_burden,color=env),fun=median, geom="point", shape=95, size=9,stroke=edgewidth,alpha=1)+
  geom_jitter(aes(x=as.factor(env),y=tot_burden,fill=env),shape=21, colour='black',size=markersize, alpha=alpha,width=width,stroke=edgewidth) +
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
  guides(fill = guide_legend(byrow = TRUE))+
  ylim(70000, 97000)
H2_total <- ggplot(H2) + 
  geom_hline(yintercept=(mean(founder$tot_burden)),color='black',linetype='dashed',linewidth=0.7)+
  stat_summary(aes(x=as.factor(env),y=tot_burden,color=env),fun=median, geom="point", shape=95, size=12,stroke=edgewidth,alpha=1)+
  geom_jitter(aes(x=as.factor(env),y=tot_burden,fill=env),shape=21, colour='black',size=markersize, alpha=alpha,width=width,stroke=edgewidth) +
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
  ylab('SNV mutational load') +
  labs(tag='TOTAL')+
  ylim(70000, 97000)

H3 <- subset(snv_data, (env == "NaCl" | env == "Ethanol" | env =='H3F1' | env =='H3F2'))
H3$env <- factor(H3$env , levels=c('NaCl','Ethanol','H3F1','H3F2'))
H3_plot <- ggplot(H3) + 
  geom_hline(yintercept=(mean(founder$tot_burden)),color='black',linetype='dashed',linewidth=0.7)+
  stat_summary(aes(x=as.factor(env),y=tot_burden,color=env),fun=median, geom="point", shape=95, size=9,stroke=edgewidth,alpha=1)+
  geom_jitter(aes(x=as.factor(env),y=tot_burden,fill=env),shape=21, colour='black',size=markersize, alpha=alpha,width=width,stroke=edgewidth) +
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
  guides(fill = guide_legend(byrow = TRUE))+
  ylim(70000, 97000)
H3_total <- ggplot(H3) + 
  geom_hline(yintercept=(mean(founder$tot_burden)),color='black',linetype='dashed',linewidth=0.7)+
  stat_summary(aes(x=as.factor(env),y=tot_burden,color=env),fun=median, geom="point", shape=95, size=12,stroke=edgewidth,alpha=1)+
  geom_jitter(aes(x=as.factor(env),y=tot_burden,fill=env),shape=21, colour='black',size=markersize, alpha=alpha,width=width,stroke=edgewidth) +
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
  ylab('SNV mutational load') +
  labs(tag='TOTAL')+
  ylim(70000, 97000)

H1_total + H1_plot + H2_total + H2_plot + H3_total + H3_plot + 
  plot_layout(ncol=2,widths = c(1.3,8))

ggsave('figures/Hybrid_Dynamics_SNV_burden.pdf',width=11,height=5,dpi = 900)

# ==============================================================================
# HIGH IMPACT SNVs filtered by SnpEff
# ==============================================================================
genome_size=12071326

founder <- subset(snv_data, gen %in% "Founder")

H1 <- subset(snv_data, (env == "NaCl" | env == "LiAc0.01" | env =='H1F1' | env =='H1F2'))
H1$env <- factor(H1$env , levels=c('NaCl','LiAc0.01','H1F1','H1F2'))
H1_plot <- ggplot(H1) + 
  geom_hline(yintercept=(mean(founder$high)),color='black',linetype='dashed',linewidth=0.7)+
  stat_summary(aes(x=as.factor(env),y=high,color=env),fun=median, geom="point", shape=95, size=9,stroke=edgewidth,alpha=1)+
  geom_jitter(aes(x=as.factor(env),y=high,fill=env),shape=21, colour='black',size=markersize, alpha=alpha,width=width,stroke=edgewidth) +
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
  geom_hline(yintercept=(mean(founder$high)),color='black',linetype='dashed',linewidth=0.7)+
  stat_summary(aes(x=as.factor(env),y=high,color=env),fun=median, geom="point", shape=95, size=12,stroke=edgewidth,alpha=1)+
  geom_jitter(aes(x=as.factor(env),y=high,fill=env),shape=21, colour='black',size=markersize, alpha=alpha,width=width,stroke=edgewidth) +
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
  ylab('SNV load | high impact') +
  labs(tag='TOTAL')

H2 <- subset(snv_data, (env == "NaCl" | env == "LiAc0.02" | env =='H2F1' | env =='H2F2'))
H2$env <- factor(H2$env , levels=c('NaCl','LiAc0.02','H2F1','H2F2'))
H2_plot <- ggplot(H2) + 
  geom_hline(yintercept=(mean(founder$high)),color='black',linetype='dashed',linewidth=0.7)+
  stat_summary(aes(x=as.factor(env),y=high,color=env),fun=median, geom="point", shape=95, size=9,stroke=edgewidth,alpha=1)+
  geom_jitter(aes(x=as.factor(env),y=high,fill=env),shape=21, colour='black',size=markersize, alpha=alpha,width=width,stroke=edgewidth) +
  facet_wrap(~factor(gen,c('G100','G200','G300','G400','G500','G700','G1000')), ncol=7)+
  scale_fill_manual(values=c(NaCl,LiAc0.02,F1,F2),name ='',labels = c("NaCl", "LiAc0.02", "Hybrid F1",'Hybrid F2'))+
  scale_color_manual(values=c(NaCl,LiAc0.02,F1,F2),name ='',labels = c("NaCl", "LiAc0.02", "Hybrid F1",'Hybrid F2'))+
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
  geom_hline(yintercept=(mean(founder$high)),color='black',linetype='dashed',linewidth=0.7)+
  stat_summary(aes(x=as.factor(env),y=high,color=env),fun=median, geom="point", shape=95, size=12,stroke=edgewidth,alpha=1)+
  geom_jitter(aes(x=as.factor(env),y=high,fill=env),shape=21, colour='black',size=markersize, alpha=alpha,width=width,stroke=edgewidth) +
  scale_fill_manual(values=c(NaCl,LiAc0.02,F1,F2),name ='',labels = c("NaCl", "LiAc0.02", "Hybrid F1",'Hybrid F2'))+
  scale_color_manual(values=c(NaCl,LiAc0.02,F1,F2),name ='',labels = c("NaCl", "LiAc0.02", "Hybrid F1",'Hybrid F2'))+
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
  ylab('SNV load | high impact') +
  labs(tag='TOTAL')

H3 <- subset(snv_data, (env == "NaCl" | env == "Ethanol" | env =='H3F1' | env =='H3F2'))
H3$env <- factor(H3$env , levels=c('NaCl','Ethanol','H3F1','H3F2'))
H3_plot <- ggplot(H3) + 
  geom_hline(yintercept=(mean(founder$high)),color='black',linetype='dashed',linewidth=0.7)+
  stat_summary(aes(x=as.factor(env),y=high,color=env),fun=median, geom="point", shape=95, size=9,stroke=edgewidth,alpha=1)+
  geom_jitter(aes(x=as.factor(env),y=high,fill=env),shape=21, colour='black',size=markersize, alpha=alpha,width=width,stroke=edgewidth) +
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
  geom_hline(yintercept=(mean(founder$high)),color='black',linetype='dashed',linewidth=0.7)+
  stat_summary(aes(x=as.factor(env),y=high,color=env),fun=median, geom="point", shape=95, size=12,stroke=edgewidth,alpha=1)+
  geom_jitter(aes(x=as.factor(env),y=high,fill=env),shape=21, colour='black',size=markersize, alpha=alpha,width=width,stroke=edgewidth) +
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
  ylab('SNV load | high impact') +
  labs(tag='TOTAL')

H1_total + H1_plot + H2_total + H2_plot + H3_total + H3_plot + 
  plot_layout(ncol=2,widths = c(1.3,8))

ggsave('figures/Hybrid_Dynamics_SNV_burden-high.pdf',width=11,height=5,dpi = 900)

# ==============================================================================
# HIGH IMPACT SNVs filtered by SnpEff - panel for Figure 2
# ==============================================================================

founder <- subset(snv_data, gen %in% "Founder")

H1 <- subset(snv_data, (env == "NaCl" | env == "LiAc0.01" | env =='H1F1' | env =='H1F2'))
H1$env <- factor(H1$env , levels=c('NaCl','LiAc0.01','H1F1','H1F2'))
H1_total <- ggplot(H1) + 
  geom_hline(yintercept=(mean(founder$high)),color='black',linetype='dashed',linewidth=0.7)+
  stat_summary(aes(x=as.factor(env),y=high,color=env),fun=median, geom="point", shape=95, size=12,stroke=edgewidth,alpha=1)+
  geom_jitter(aes(x=as.factor(env),y=high,fill=env),shape=21, colour='black',size=markersize, alpha=alpha,width=width,stroke=edgewidth) +
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
        panel.spacing = unit(0., "lines"), 
        text=element_text(family="EB Garamond"),
        legend.position = "none",
        plot.tag.position=c(0.62,0.94),
        plot.tag = element_text(hjust = 0.5),
        plot.margin = margin(0,0,0,0, "cm")) +  
  ylab('') +
  ylim(550,765)

H2 <- subset(snv_data, (env == "NaCl" | env == "LiAc0.02" | env =='H2F1' | env =='H2F2'))
H2$env <- factor(H2$env , levels=c('NaCl','LiAc0.02','H2F1','H2F2'))
H2_total <- ggplot(H2) + 
  geom_hline(yintercept=(mean(founder$high)),color='black',linetype='dashed',linewidth=0.7)+
  stat_summary(aes(x=as.factor(env),y=high,color=env),fun=median, geom="point", shape=95, size=12,stroke=edgewidth,alpha=1)+
  geom_jitter(aes(x=as.factor(env),y=high,fill=env),shape=21, colour='black',size=markersize, alpha=alpha,width=width,stroke=edgewidth) +
  scale_fill_manual(values=c(NaCl,LiAc0.02,F1,F2),name ='',labels = c("NaCl", "LiAc0.02", "Hybrid F1",'Hybrid F2'))+
  scale_color_manual(values=c(NaCl,LiAc0.02,F1,F2),name ='',labels = c("NaCl", "LiAc0.02", "Hybrid F1",'Hybrid F2'))+
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        panel.background = element_rect(fill = "lightgrey"),
        panel.grid.major.x = element_line(colour = "white", linetype = 1, linewidth = 0.5),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(colour = "white", linetype = 1, linewidth = 0.5),
        panel.grid.minor.y = element_blank(),
        panel.spacing = unit(0., "lines"), 
        text=element_text(family="EB Garamond"),
        legend.position = "none",
        plot.tag.position=c(0.62,0.94),
        plot.tag = element_text(hjust = 0.5),
        plot.margin = margin(0,0,0,0, "cm")) +  
  ylab('') +
  ylim(550,765)

H3 <- subset(snv_data, (env == "NaCl" | env == "Ethanol" | env =='H3F1' | env =='H3F2'))
H3$env <- factor(H3$env , levels=c('NaCl','Ethanol','H3F1','H3F2'))
H3_total <- ggplot(H3) + 
  geom_hline(yintercept=(mean(founder$high)),color='black',linetype='dashed',linewidth=0.7)+
  stat_summary(aes(x=as.factor(env),y=high,color=env),fun=median, geom="point", shape=95, size=12,stroke=edgewidth,alpha=1)+
  geom_jitter(aes(x=as.factor(env),y=high,fill=env),shape=21, colour='black',size=markersize, alpha=alpha,width=width,stroke=edgewidth) +
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
        panel.spacing = unit(0., "lines"), 
        text=element_text(family="EB Garamond"),
        legend.position = "none",
        plot.tag.position=c(0.62,0.94),
        plot.tag = element_text(hjust = 0.5),
        plot.margin = margin(0,0,0,0, "cm")) + 
  ylab('') +
  ylim(550,765)

H1_total + H2_total + H3_total + 
  plot_layout(ncol=1)

ggsave('figures/Hybrid_Dynamics_SNV_burden-high_total.pdf',width=3,height=2.5,dpi = 900)
# ==============================================================================
# ANOVA of SNV mutational load
# ==============================================================================
res_aov <- aov(tot_burden ~ env,
               data = H1)
summary(res_aov)
post_test <- glht(res_aov,
                  linfct = mcp(env = "Tukey"))
summary(post_test)
plot(post_test)

res_aov <- aov(tot_burden ~ env,
               data = H2)
summary(res_aov)
post_test <- glht(res_aov,
                  linfct = mcp(env = "Tukey"))
summary(post_test)
plot(post_test)

res_aov <- aov(tot_burden ~ env,
               data = H3)
summary(res_aov)
post_test <- glht(res_aov,
                  linfct = mcp(env = "Tukey"))
summary(post_test)
plot(post_test)

# ==============================================================================
# ANOVA of SNV mutational load | HIGH IMPACT
# ==============================================================================
res_aov <- aov(high ~ env,
               data = H1)
summary(res_aov)
post_test <- glht(res_aov,
                  linfct = mcp(env = "Tukey"))
summary(post_test)
plot(post_test)

res_aov <- aov(high ~ env,
               data = H2)
summary(res_aov)
post_test <- glht(res_aov,
                  linfct = mcp(env = "Tukey"))
summary(post_test)
plot(post_test)

res_aov <- aov(high ~ env,
               data = H3)
summary(res_aov)
post_test <- glht(res_aov,
                  linfct = mcp(env = "Tukey"))
summary(post_test)
plot(post_test)

# ==============================================================================
# PLOT CORRELATIONS BETWEEN GENERATIONAL TIME AND MUTATIONAL BURDEN
# ==============================================================================
hybrid <- subset(snv_data, (env =='H1F1' | env =='H1F2'|env =='H2F1' | env =='H2F2'|env =='H3F1' | env =='H3F2'))
hybrid$gen_num <- as.numeric(sub('G','',hybrid$gen))
hybrid$env <- sub('H.','',hybrid$env)

hybrid_corr <- ggplot(hybrid, aes(x=gen_num, y=tot_burden)) + 
  geom_smooth(aes(color=env),method=lm, se=FALSE, fullrange=TRUE)+
  stat_cor(aes(color = env),label.x.npc = 0.4,label.y.npc = 0.2,p.accuracy=0.01,r.accuracy=0.01,size=2.5)+
  geom_point(alpha=alpha,aes(fill=env),shape=21, colour='black')+
  scale_fill_manual(values=c(F1,F2),name ='',labels = c("Hybrid F1",'Hybrid F2'))+ 
  scale_color_manual(values=c(F1,F2),name ='',labels = c("Hybrid F1",'Hybrid F2'))+
  theme(
        panel.background = element_rect(fill = "lightgrey"),
        panel.grid.major.x = element_blank(),#element_line(colour = "lightgrey", linetype = 1, linewidth = 0.5),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(colour = "white", linetype = 2, linewidth = 0.5),
        panel.grid.minor.y = element_blank(),
        panel.spacing = unit(0.5, "lines"), 
        text=element_text(family="EB Garamond"),
        legend.position = "none",
        plot.tag.position=c(0.62,0.94),
        plot.tag = element_text(hjust = 0.5)) + 
  ylab('SNV mutational load') +
  xlab('generations') +
  ylim(70000, 97000)
hybrid_corr
H1$gen_num <- as.numeric(sub('G','',H1$gen))
H1_corr <- ggplot(H1, aes(x=gen_num, y=tot_burden)) + 
  geom_smooth(aes(color=env),method=lm, se=FALSE, fullrange=TRUE)+
  stat_cor(aes(color = env),label.x.npc = 0.4,label.y.npc = 0.25,p.accuracy=0.01,r.accuracy=0.01,size=2.5)+
  geom_point(alpha=alpha,aes(fill=env),shape=21, colour='black')+
  scale_color_manual(values=c(NaCl,LiAc0.01,F1,F2),name ='',labels = c("NaCl", "LiAc0.01", "Hybrid F1",'Hybrid F2'))+ 
  scale_fill_manual(values=c(NaCl,LiAc0.01,F1,F2),name ='',labels = c("NaCl", "LiAc0.01", "Hybrid F1",'Hybrid F2'))+
  theme(axis.title.x=element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        panel.grid.major.x = element_blank(),#element_line(colour = "lightgrey", linetype = 1, linewidth = 0.5),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(colour = "lightgrey", linetype = 2, linewidth = 0.5),
        panel.grid.minor.y = element_blank(),
        panel.spacing = unit(0.5, "lines"), 
        text=element_text(family="EB Garamond"),
        legend.position = "none")+
  ylim(70000, 97000)
H1_corr
H2$gen_num <- as.numeric(sub('G','',H2$gen))
H2_corr <- ggplot(H2, aes(x=gen_num, y=tot_burden)) + 
  geom_smooth(aes(color=env),method=lm, se=FALSE, fullrange=TRUE)+
  stat_cor(aes(color = env),label.x.npc = 0.4,label.y.npc = 1,p.accuracy=0.01,r.accuracy=0.01,size=2.5)+
  geom_point(alpha=alpha,aes(fill=env),shape=21, colour='black')+
  scale_color_manual(values=c(NaCl,LiAc0.02,F1,F2),name ='',labels = c("NaCl", "LiAc0.02", "Hybrid F1",'Hybrid F2'))+ 
  scale_fill_manual(values=c(NaCl,LiAc0.02,F1,F2),name ='',labels = c("NaCl", "LiAc0.02", "Hybrid F1",'Hybrid F2'))+ 
  theme(axis.title.x=element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        panel.grid.major.x = element_blank(),#element_line(colour = "lightgrey", linetype = 1, linewidth = 0.5),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(colour = "lightgrey", linetype = 2, linewidth = 0.5),
        panel.grid.minor.y = element_blank(),
        panel.spacing = unit(0.5, "lines"), 
        text=element_text(family="EB Garamond"),
        legend.position = "none")+
  ylim(70000, 97000)

H3$gen_num <- as.numeric(sub('G','',H3$gen))
H3_corr <- ggplot(H3, aes(x=gen_num, y=tot_burden)) + 
  geom_smooth(aes(color=env),method=lm, se=FALSE, fullrange=TRUE)+
  stat_cor(aes(color = env),label.x.npc = 0.4,label.y.npc = 0.25,p.accuracy=0.01,r.accuracy=0.01,size=2.5)+
  geom_point(alpha=alpha,aes(fill=env),shape=21, colour='black')+
  scale_color_manual(values=c(NaCl,Ethanol,F1,F2),name ='',labels = c("NaCl", "Ethanol", "Hybrid F1",'Hybrid F2'))+ 
  scale_fill_manual(values=c(NaCl,Ethanol,F1,F2),name ='',labels = c("NaCl", "Ethanol", "Hybrid F1",'Hybrid F2'))+ 
  theme(axis.title.x=element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        panel.grid.major.x = element_blank(),#element_line(colour = "lightgrey", linetype = 1, linewidth = 0.5),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(colour = "lightgrey", linetype = 2, linewidth = 0.5),
        panel.grid.minor.y = element_blank(),
        panel.spacing = unit(0.5, "lines"), 
        text=element_text(family="EB Garamond"),
        legend.position = "none")+
  ylim(70000, 97000)

hybrid_corr | H1_corr | H2_corr | H3_corr
ggsave('figures/Hybrid_Dynamics_SNV_burden_corr.pdf',width=8,height=2.5,dpi = 900)


# ==============================================================================
# PLOT CORRELATIONS BETWEEN GENERATIONAL TIME AND MUTATIONAL BURDEN - HIGH
# ==============================================================================
hybrid <- subset(snv_data, (env =='H1F1' | env =='H1F2'|env =='H2F1' | env =='H2F2'|env =='H3F1' | env =='H3F2'))
hybrid$gen_num <- as.numeric(sub('G','',hybrid$gen))
hybrid$env <- sub('H.','',hybrid$env)

hybrid_corr <- ggplot(hybrid, aes(x=gen_num, y=high)) + 
  geom_point(alpha=1,aes(fill=env),shape=21, colour='black')+
  geom_smooth(aes(color=env),method=lm, se=FALSE, fullrange=TRUE)+
  stat_cor(aes(color = env),label.x.npc = 0.3,label.y.npc = 0.2,p.accuracy=0.01,r.accuracy=0.01,size=2.5)+
  scale_fill_manual(values=c(F1,F2),name ='',labels = c("Hybrid F1",'Hybrid F2'))+ 
  scale_color_manual(values=c(F1,F2),name ='',labels = c("Hybrid F1",'Hybrid F2'))+
  theme(
    panel.background = element_rect(fill = "lightgrey"),
    panel.grid.major.x = element_blank(),#element_line(colour = "lightgrey", linetype = 1, linewidth = 0.5),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(colour = "white", linetype = 2, linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    panel.spacing = unit(0.5, "lines"), 
    text=element_text(family="EB Garamond"),
    legend.position = "none",
    plot.tag.position=c(0.62,0.94),
    plot.tag = element_text(hjust = 0.5)) + 
  ylab('SNV load | high impact') +
  xlab('generations') +
  ylim(540,770)

H1$gen_num <- as.numeric(sub('G','',H1$gen))
H1_corr <- ggplot(H1, aes(x=gen_num, y=high)) + 
  geom_point(alpha=1,aes(fill=env),shape=21, colour='black')+
  geom_smooth(aes(color=env),method=lm, se=FALSE, fullrange=TRUE)+
  stat_cor(aes(color = env),label.x.npc = 0.3,label.y.npc = 0.2,p.accuracy=0.01,r.accuracy=0.01,size=2.5)+
  scale_color_manual(values=c(NaCl,LiAc0.01,F1,F2),name ='',labels = c("NaCl", "LiAc0.01", "Hybrid F1",'Hybrid F2'))+ 
  scale_fill_manual(values=c(NaCl,LiAc0.01,F1,F2),name ='',labels = c("NaCl", "LiAc0.01", "Hybrid F1",'Hybrid F2'))+
  theme(axis.title.x=element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        panel.grid.major.x = element_blank(),#element_line(colour = "lightgrey", linetype = 1, linewidth = 0.5),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(colour = "lightgrey", linetype = 2, linewidth = 0.5),
        panel.grid.minor.y = element_blank(),
        panel.spacing = unit(0.5, "lines"), 
        text=element_text(family="EB Garamond"),
        legend.position = "none")+
  ylim(540,770)
H2$gen_num <- as.numeric(sub('G','',H2$gen))
H2_corr <- ggplot(H2, aes(x=gen_num, y=high)) + 
  geom_point(alpha=1,aes(fill=env),shape=21, colour='black')+
  geom_smooth(aes(color=env),method=lm, se=FALSE, fullrange=TRUE)+
  stat_cor(aes(color = env),label.x.npc = 0.3,label.y.npc = 1,p.accuracy=0.01,r.accuracy=0.01,size=2.5)+
  scale_color_manual(values=c(NaCl,LiAc0.02,F1,F2),name ='',labels = c("NaCl", "LiAc0.02", "Hybrid F1",'Hybrid F2'))+ 
  scale_fill_manual(values=c(NaCl,LiAc0.02,F1,F2),name ='',labels = c("NaCl", "LiAc0.02", "Hybrid F1",'Hybrid F2'))+ 
  theme(axis.title.x=element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        panel.grid.major.x = element_blank(),#element_line(colour = "lightgrey", linetype = 1, linewidth = 0.5),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(colour = "lightgrey", linetype = 2, linewidth = 0.5),
        panel.grid.minor.y = element_blank(),
        panel.spacing = unit(0.5, "lines"), 
        text=element_text(family="EB Garamond"),
        legend.position = "none")+
  ylim(540,770)

H3$gen_num <- as.numeric(sub('G','',H3$gen))
H3_corr <- ggplot(H3, aes(x=gen_num, y=high)) + 
  geom_point(alpha=1,aes(fill=env),shape=21, colour='black')+
  geom_smooth(aes(color=env),method=lm, se=FALSE, fullrange=TRUE)+
  stat_cor(aes(color = env),label.x.npc = 0.3,label.y.npc = 0.2,p.accuracy=0.01,r.accuracy=0.01,size=2.5)+
  scale_color_manual(values=c(NaCl,Ethanol,F1,F2),name ='',labels = c("NaCl", "Ethanol", "Hybrid F1",'Hybrid F2'))+ 
  scale_fill_manual(values=c(NaCl,Ethanol,F1,F2),name ='',labels = c("NaCl", "Ethanol", "Hybrid F1",'Hybrid F2'))+ 
  theme(axis.title.x=element_blank(),
        panel.background = element_rect(fill = "white"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        panel.grid.major.x = element_blank(),#element_line(colour = "lightgrey", linetype = 1, linewidth = 0.5),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(colour = "lightgrey", linetype = 2, linewidth = 0.5),
        panel.grid.minor.y = element_blank(),
        panel.spacing = unit(0.5, "lines"), 
        text=element_text(family="EB Garamond"),
        legend.position = "none")+
  ylim(540,770)

hybrid_corr | H1_corr | H2_corr | H3_corr
ggsave('figures/Hybrid_Dynamics_SNV_burden_corr-high.pdf',width=8,height=2.5,dpi = 900)


# ==============================================================================
# FUNCTION TO BE ABLE TO CALCULATE THE SUMMARY (MEAN AND STANDARD DEVIATION)
# ==============================================================================
install.packages('plyr')
data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

# ==============================================================================
# PLOT THE PREVALENCE OF EACH MUTATIONAL TYPE IN THE TOTAL
# ==============================================================================

#H1 - NaCL x LiAc0.01
metric='SNP'
hybrid_type=H1
H1_stat <- data_summary(hybrid_type, varname=metric, groupnames=c('env'))
H1_stat$env <- factor(H1_stat$env , levels=c('H1F2','H1F1','LiAc0.01','NaCl'))
H1_SNP <- ggplot(H1_stat) + 
  geom_errorbar(aes(x=.data[[metric]],y=env,xmin=(.data[[metric]])-sd, xmax=(.data[[metric]])+sd), width=0)+
  geom_point( aes(x=.data[[metric]],y=env,fill=env),size=3, shape=21, color='black')+
  scale_fill_manual(values=c(F2,F1,LiAc0.01,NaCl),name ='',labels = c('Hybrid F2','Hybrid F1','LiAc 0.01','NaCl'))+
  scale_color_manual(values=c(F2,F1,LiAc0.01,NaCl),name ='',labels = c('Hybrid F2','Hybrid F1','LiAc 0.01','NaCl'))+
  theme(panel.background = element_rect(fill = "white"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major.y = element_blank(),#element_line(colour = "lightgrey", linetype = 1, linewidth = 0.5),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(colour = "lightgrey", linetype = 2, linewidth = 0.5),
        panel.grid.minor.x = element_blank(),
        panel.spacing = unit(0.5, "lines"), 
        text=element_text(family="EB Garamond"),
        legend.position = "none")

metric='DEL'
H1_stat <- data_summary(hybrid_type, varname=metric, groupnames=c('env'))
H1_stat$env <- factor(H1_stat$env , levels=c('H1F2','H1F1','LiAc0.01','NaCl'))
H1_DEL <- ggplot(H1_stat) + 
  geom_errorbar(aes(x=.data[[metric]],y=env,xmin=(.data[[metric]])-sd, xmax=(.data[[metric]])+sd), width=0)+
  geom_point( aes(x=.data[[metric]],y=env,fill=env),size=3, shape=21, color='black')+
  scale_fill_manual(values=c(F2,F1,LiAc0.01,NaCl),name ='',labels = c('Hybrid F2','Hybrid F1','LiAc 0.01','NaCl'))+
  scale_color_manual(values=c(F2,F1,LiAc0.01,NaCl),name ='',labels = c('Hybrid F2','Hybrid F1','LiAc 0.01','NaCl'))+
  theme(panel.background = element_rect(fill = "white"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major.y = element_blank(),#element_line(colour = "lightgrey", linetype = 1, linewidth = 0.5),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(colour = "lightgrey", linetype = 2, linewidth = 0.5),
        panel.grid.minor.x = element_blank(),
        panel.spacing = unit(0.5, "lines"), 
        text=element_text(family="EB Garamond"),
        legend.position = "none")

metric='INS'
H1_stat <- data_summary(hybrid_type, varname=metric, groupnames=c('env'))
H1_stat$env <- factor(H1_stat$env , levels=c('H1F2','H1F1','LiAc0.01','NaCl'))
H1_INS <- ggplot(H1_stat) + 
  geom_errorbar(aes(x=.data[[metric]],y=env,xmin=(.data[[metric]])-sd, xmax=(.data[[metric]])+sd), width=0)+
  geom_point( aes(x=.data[[metric]],y=env,fill=env),size=3, shape=21, color='black')+
  scale_fill_manual(values=c(F2,F1,LiAc0.01,NaCl),name ='',labels = c('Hybrid F2','Hybrid F1','LiAc 0.01','NaCl'))+
  scale_color_manual(values=c(F2,F1,LiAc0.01,NaCl),name ='',labels = c('Hybrid F2','Hybrid F1','LiAc 0.01','NaCl'))+
  theme(panel.background = element_rect(fill = "white"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major.y = element_blank(),#element_line(colour = "lightgrey", linetype = 1, linewidth = 0.5),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(colour = "lightgrey", linetype = 2, linewidth = 0.5),
        panel.grid.minor.x = element_blank(),
        panel.spacing = unit(0.5, "lines"), 
        text=element_text(family="EB Garamond"),
        legend.position = "none")

metric='MIXED'
H1_stat <- data_summary(hybrid_type, varname=metric, groupnames=c('env'))
H1_stat$env <- factor(H1_stat$env , levels=c('H1F2','H1F1','LiAc0.01','NaCl'))
H1_MIXED <- ggplot(H1_stat) + 
  geom_errorbar(aes(x=.data[[metric]],y=env,xmin=(.data[[metric]])-sd, xmax=(.data[[metric]])+sd), width=0)+
  geom_point( aes(x=.data[[metric]],y=env,fill=env),size=3, shape=21, color='black')+
  scale_fill_manual(values=c(F2,F1,LiAc0.01,NaCl),name ='',labels = c('Hybrid F2','Hybrid F1','LiAc 0.01','NaCl'))+
  scale_color_manual(values=c(F2,F1,LiAc0.01,NaCl),name ='',labels = c('Hybrid F2','Hybrid F1','LiAc 0.01','NaCl'))+
  theme(panel.background = element_rect(fill = "white"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major.y = element_blank(),#element_line(colour = "lightgrey", linetype = 1, linewidth = 0.5),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(colour = "lightgrey", linetype = 2, linewidth = 0.5),
        panel.grid.minor.x = element_blank(),
        panel.spacing = unit(0.5, "lines"), 
        text=element_text(family="EB Garamond"),
        legend.position = "none")

metric='MNP'
H1_stat <- data_summary(hybrid_type, varname=metric, groupnames=c('env'))
H1_stat$env <- factor(H1_stat$env , levels=c('H1F2','H1F1','LiAc0.01','NaCl'))
H1_MNP <- ggplot(H1_stat) + 
  geom_errorbar(aes(x=.data[[metric]],y=env,xmin=(.data[[metric]])-sd, xmax=(.data[[metric]])+sd), width=0)+
  geom_point( aes(x=.data[[metric]],y=env,fill=env),size=3, shape=21, color='black')+
  scale_fill_manual(values=c(F2,F1,LiAc0.01,NaCl),name ='',labels = c('Hybrid F2','Hybrid F1','LiAc 0.01','NaCl'))+
  scale_color_manual(values=c(F2,F1,LiAc0.01,NaCl),name ='',labels = c('Hybrid F2','Hybrid F1','LiAc 0.01','NaCl'))+
  theme(panel.background = element_rect(fill = "white"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major.y = element_blank(),#element_line(colour = "lightgrey", linetype = 1, linewidth = 0.5),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(colour = "lightgrey", linetype = 2, linewidth = 0.5),
        panel.grid.minor.x = element_blank(),
        panel.spacing = unit(0.5, "lines"), 
        text=element_text(family="EB Garamond"),
        legend.position = "right")

#H2 - NaCL x LiAc0.02
metric='SNP'
hybrid_type=H2
H2_stat <- data_summary(hybrid_type, varname=metric, groupnames=c('env'))
H2_stat$env <- factor(H2_stat$env , levels=c('H2F2','H2F1','LiAc0.02','NaCl'))
H2_SNP <- ggplot(H2_stat) + 
  geom_errorbar(aes(x=.data[[metric]],y=env,xmin=(.data[[metric]])-sd, xmax=(.data[[metric]])+sd), width=0)+
  geom_point( aes(x=.data[[metric]],y=env,fill=env),size=3, shape=21, color='black')+
  scale_fill_manual(values=c(F2,F1,LiAc0.02,NaCl),name ='',labels = c('Hybrid F2','Hybrid F1','LiAc 0.02','NaCl'))+
  scale_color_manual(values=c(F2,F1,LiAc0.02,NaCl),name ='',labels = c('Hybrid F2','Hybrid F1','LiAc 0.02','NaCl'))+
  theme(panel.background = element_rect(fill = "white"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major.y = element_blank(),#element_line(colour = "lightgrey", linetype = 1, linewidth = 0.5),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(colour = "lightgrey", linetype = 2, linewidth = 0.5),
        panel.grid.minor.x = element_blank(),
        panel.spacing = unit(0.5, "lines"), 
        text=element_text(family="EB Garamond"),
        legend.position = "none")

metric='DEL'
H2_stat <- data_summary(hybrid_type, varname=metric, groupnames=c('env'))
H2_stat$env <- factor(H2_stat$env , levels=c('H2F2','H2F1','LiAc0.02','NaCl'))
H2_DEL <- ggplot(H2_stat) + 
  geom_errorbar(aes(x=.data[[metric]],y=env,xmin=(.data[[metric]])-sd, xmax=(.data[[metric]])+sd), width=0)+
  geom_point( aes(x=.data[[metric]],y=env,fill=env),size=3, shape=21, color='black')+
  scale_fill_manual(values=c(F2,F1,LiAc0.02,NaCl),name ='',labels = c('Hybrid F2','Hybrid F1','LiAc 0.02','NaCl'))+
  scale_color_manual(values=c(F2,F1,LiAc0.02,NaCl),name ='',labels = c('Hybrid F2','Hybrid F1','LiAc 0.02','NaCl'))+
  theme(panel.background = element_rect(fill = "white"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major.y = element_blank(),#element_line(colour = "lightgrey", linetype = 1, linewidth = 0.5),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(colour = "lightgrey", linetype = 2, linewidth = 0.5),
        panel.grid.minor.x = element_blank(),
        panel.spacing = unit(0.5, "lines"), 
        text=element_text(family="EB Garamond"),
        legend.position = "none")

metric='INS'
H2_stat <- data_summary(hybrid_type, varname=metric, groupnames=c('env'))
H2_stat$env <- factor(H2_stat$env , levels=c('H2F2','H2F1','LiAc0.02','NaCl'))
H2_INS <- ggplot(H2_stat) + 
  geom_errorbar(aes(x=.data[[metric]],y=env,xmin=(.data[[metric]])-sd, xmax=(.data[[metric]])+sd), width=0)+
  geom_point( aes(x=.data[[metric]],y=env,fill=env),size=3, shape=21, color='black')+
  scale_fill_manual(values=c(F2,F1,LiAc0.02,NaCl),name ='',labels = c('Hybrid F2','Hybrid F1','LiAc 0.02','NaCl'))+
  scale_color_manual(values=c(F2,F1,LiAc0.02,NaCl),name ='',labels = c('Hybrid F2','Hybrid F1','LiAc 0.02','NaCl'))+
  theme(panel.background = element_rect(fill = "white"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major.y = element_blank(),#element_line(colour = "lightgrey", linetype = 1, linewidth = 0.5),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(colour = "lightgrey", linetype = 2, linewidth = 0.5),
        panel.grid.minor.x = element_blank(),
        panel.spacing = unit(0.5, "lines"), 
        text=element_text(family="EB Garamond"),
        legend.position = "none")

metric='MIXED'
H2_stat <- data_summary(hybrid_type, varname=metric, groupnames=c('env'))
H2_stat$env <- factor(H2_stat$env , levels=c('H2F2','H2F1','LiAc0.02','NaCl'))
H2_MIXED <- ggplot(H2_stat) + 
  geom_errorbar(aes(x=.data[[metric]],y=env,xmin=(.data[[metric]])-sd, xmax=(.data[[metric]])+sd), width=0)+
  geom_point( aes(x=.data[[metric]],y=env,fill=env),size=3, shape=21, color='black')+
  scale_fill_manual(values=c(F2,F1,LiAc0.02,NaCl),name ='',labels = c('Hybrid F2','Hybrid F1','LiAc 0.02','NaCl'))+
  scale_color_manual(values=c(F2,F1,LiAc0.02,NaCl),name ='',labels = c('Hybrid F2','Hybrid F1','LiAc 0.02','NaCl'))+
  theme(panel.background = element_rect(fill = "white"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major.y = element_blank(),#element_line(colour = "lightgrey", linetype = 1, linewidth = 0.5),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(colour = "lightgrey", linetype = 2, linewidth = 0.5),
        panel.grid.minor.x = element_blank(),
        panel.spacing = unit(0.5, "lines"), 
        text=element_text(family="EB Garamond"),
        legend.position = "none")

metric='MNP'
H2_stat <- data_summary(hybrid_type, varname=metric, groupnames=c('env'))
H2_stat$env <- factor(H2_stat$env , levels=c('H2F2','H2F1','LiAc0.02','NaCl'))
H2_MNP <- ggplot(H2_stat) + 
  geom_errorbar(aes(x=.data[[metric]],y=env,xmin=(.data[[metric]])-sd, xmax=(.data[[metric]])+sd), width=0)+
  geom_point( aes(x=.data[[metric]],y=env,fill=env),size=3, shape=21, color='black')+
  scale_fill_manual(values=c(F2,F1,LiAc0.02,NaCl),name ='',labels = c('Hybrid F2','Hybrid F1','LiAc 0.02','NaCl'))+
  scale_color_manual(values=c(F2,F1,LiAc0.02,NaCl),name ='',labels = c('Hybrid F2','Hybrid F1','LiAc 0.02','NaCl'))+
  theme(panel.background = element_rect(fill = "white"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major.y = element_blank(),#element_line(colour = "lightgrey", linetype = 1, linewidth = 0.5),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(colour = "lightgrey", linetype = 2, linewidth = 0.5),
        panel.grid.minor.x = element_blank(),
        panel.spacing = unit(0.5, "lines"), 
        text=element_text(family="EB Garamond"),
        legend.position = "right")

#H3 - NaCL x Ethanol
metric='SNP'
hybrid_type=H3
H3_stat <- data_summary(hybrid_type, varname=metric, groupnames=c('env'))
H3_stat$env <- factor(H3_stat$env , levels=c('H3F2','H3F1','Ethanol','NaCl'))
H3_SNP <- ggplot(H3_stat) + 
  geom_errorbar(aes(x=.data[[metric]],y=env,xmin=(.data[[metric]])-sd, xmax=(.data[[metric]])+sd), width=0)+
  geom_point( aes(x=.data[[metric]],y=env,fill=env),size=3, shape=21, color='black')+
  scale_fill_manual(values=c(F2,F1,Ethanol,NaCl),name ='',labels = c('Hybrid F2','Hybrid F1','Ethanol','NaCl'))+
  scale_color_manual(values=c(F2,F1,Ethanol,NaCl),name ='',labels = c('Hybrid F2','Hybrid F1','Ethanol','NaCl'))+
  theme(panel.background = element_rect(fill = "white"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major.y = element_blank(),#element_line(colour = "lightgrey", linetype = 1, linewidth = 0.5),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(colour = "lightgrey", linetype = 2, linewidth = 0.5),
        panel.grid.minor.x = element_blank(),
        panel.spacing = unit(0.5, "lines"), 
        text=element_text(family="EB Garamond"),
        legend.position = "none")

metric='DEL'
H3_stat <- data_summary(hybrid_type, varname=metric, groupnames=c('env'))
H3_stat$env <- factor(H3_stat$env , levels=c('H3F2','H3F1','Ethanol','NaCl'))
H3_DEL <- ggplot(H3_stat) + 
  geom_errorbar(aes(x=.data[[metric]],y=env,xmin=(.data[[metric]])-sd, xmax=(.data[[metric]])+sd), width=0)+
  geom_point( aes(x=.data[[metric]],y=env,fill=env),size=3, shape=21, color='black')+
  scale_fill_manual(values=c(F2,F1,Ethanol,NaCl),name ='',labels = c('Hybrid F2','Hybrid F1','Ethanol','NaCl'))+
  scale_color_manual(values=c(F2,F1,Ethanol,NaCl),name ='',labels = c('Hybrid F2','Hybrid F1','Ethanol','NaCl'))+
  theme(panel.background = element_rect(fill = "white"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major.y = element_blank(),#element_line(colour = "lightgrey", linetype = 1, linewidth = 0.5),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(colour = "lightgrey", linetype = 2, linewidth = 0.5),
        panel.grid.minor.x = element_blank(),
        panel.spacing = unit(0.5, "lines"), 
        text=element_text(family="EB Garamond"),
        legend.position = "none")

metric='INS'
H3_stat <- data_summary(hybrid_type, varname=metric, groupnames=c('env'))
H3_stat$env <- factor(H3_stat$env , levels=c('H3F2','H3F1','Ethanol','NaCl'))
H3_INS <- ggplot(H3_stat) + 
  geom_errorbar(aes(x=.data[[metric]],y=env,xmin=(.data[[metric]])-sd, xmax=(.data[[metric]])+sd), width=0)+
  geom_point( aes(x=.data[[metric]],y=env,fill=env),size=3, shape=21, color='black')+
  scale_fill_manual(values=c(F2,F1,Ethanol,NaCl),name ='',labels = c('Hybrid F2','Hybrid F1','Ethanol','NaCl'))+
  scale_color_manual(values=c(F2,F1,Ethanol,NaCl),name ='',labels = c('Hybrid F2','Hybrid F1','Ethanol','NaCl'))+
  theme(panel.background = element_rect(fill = "white"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major.y = element_blank(),#element_line(colour = "lightgrey", linetype = 1, linewidth = 0.5),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(colour = "lightgrey", linetype = 2, linewidth = 0.5),
        panel.grid.minor.x = element_blank(),
        panel.spacing = unit(0.5, "lines"), 
        text=element_text(family="EB Garamond"),
        legend.position = "none")

metric='MIXED'
H3_stat <- data_summary(hybrid_type, varname=metric, groupnames=c('env'))
H3_stat$env <- factor(H3_stat$env , levels=c('H3F2','H3F1','Ethanol','NaCl'))
H3_MIXED <- ggplot(H3_stat) + 
  geom_errorbar(aes(x=.data[[metric]],y=env,xmin=(.data[[metric]])-sd, xmax=(.data[[metric]])+sd), width=0)+
  geom_point( aes(x=.data[[metric]],y=env,fill=env),size=3, shape=21, color='black')+
  scale_fill_manual(values=c(F2,F1,Ethanol,NaCl),name ='',labels = c('Hybrid F2','Hybrid F1','Ethanol','NaCl'))+
  scale_color_manual(values=c(F2,F1,Ethanol,NaCl),name ='',labels = c('Hybrid F2','Hybrid F1','Ethanol','NaCl'))+
  theme(panel.background = element_rect(fill = "white"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major.y = element_blank(),#element_line(colour = "lightgrey", linetype = 1, linewidth = 0.5),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(colour = "lightgrey", linetype = 2, linewidth = 0.5),
        panel.grid.minor.x = element_blank(),
        panel.spacing = unit(0.5, "lines"), 
        text=element_text(family="EB Garamond"),
        legend.position = "none")

metric='MNP'
H3_stat <- data_summary(hybrid_type, varname=metric, groupnames=c('env'))
H3_stat$env <- factor(H3_stat$env , levels=c('H3F2','H3F1','Ethanol','NaCl'))
H3_MNP <- ggplot(H3_stat) + 
  geom_errorbar(aes(x=.data[[metric]],y=env,xmin=(.data[[metric]])-sd, xmax=(.data[[metric]])+sd), width=0)+
  geom_point( aes(x=.data[[metric]],y=env,fill=env),size=3, shape=21, color='black')+
  scale_fill_manual(values=c(F2,F1,Ethanol,NaCl),name ='',labels = c('Hybrid F2','Hybrid F1','Ethanol','NaCl'))+
  scale_color_manual(values=c(F2,F1,Ethanol,NaCl),name ='',labels = c('Hybrid F2','Hybrid F1','Ethanol','NaCl'))+
  theme(panel.background = element_rect(fill = "white"),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        panel.grid.major.y = element_blank(),#element_line(colour = "lightgrey", linetype = 1, linewidth = 0.5),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_line(colour = "lightgrey", linetype = 2, linewidth = 0.5),
        panel.grid.minor.x = element_blank(),
        panel.spacing = unit(0.5, "lines"), 
        text=element_text(family="EB Garamond"),
        legend.position = "right") 

(H1_SNP | H1_DEL | H1_INS | H1_MIXED | H1_MNP) / (H2_SNP | H2_DEL | H2_INS | H2_MIXED | H2_MNP) / (H3_SNP | H3_DEL | H3_INS | H3_MIXED | H3_MNP)

ggsave('figures/Hybrid_Dynamics_SNV_burden_type.pdf',width=11,height=5,dpi = 900)

