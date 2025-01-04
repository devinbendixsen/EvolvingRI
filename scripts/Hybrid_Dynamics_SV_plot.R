#"`-:-.   ,-;"`-:-.   ,-;"`-:-.   ,-;"`-:-. -------------------------- ,-;"`-:-.
# `=`,'=/     `=`,'=/     `=`,'=/     `=`|  "If I have seen further |=/     `=
#   y==/        y==/        y==/        y| it is by standing on the |/
# ,=,-<=`.    ,=,-<=`.    ,=,-<=`.    ,=,|    shoulders of Giants"  |=`.    ,=
#,-'-'   `-=_,-'-'   `-=_,-'-'   `-=_,-'-'  -------------------------- `-=_,-'-'
# ==============================================================================
# PROJECT: Hybrid_Dynamics
# PURPOSE: to plot the SV mutational landscape
# DATE: 6 November 2023
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
library(vcfR)
library(stringr)
library(showtext)
font_add_google("EB Garamond")
showtext_auto()
quartz()
# ==============================================================================
# 
# ==============================================================================
fileNames <- Sys.glob("results/SV/manta_INV/*.vcf") # identify all samples

sv_data <- data.frame()
for (fileName in fileNames) {
  print(fileName)
  vcf_sample <- read.vcfR(fileName,verbose = FALSE) # import VCF
  mut_data_sample <- vcfR2tidy(vcf_sample,single_frame=FALSE,info_only=TRUE,toss_INFO_column = TRUE) #Extract needed data including genotype AF
  mut_data_sample <- as.data.frame(mut_data_sample$fix)
  mut_data_sample <- subset(mut_data_sample,FILTER == "PASS")
  
  
  SV_burden <- mut_data_sample %>% filter(is.na(EVENT)) # filter to SVs that haven't been split
  SV_INV <- mut_data_sample %>% distinct(EVENT, .keep_all=TRUE) #filter to inversion events split
  SV_burden <- rbind(SV_burden,SV_INV) # recombine 
  mut_data_sample <- SV_burden %>% distinct(.keep_all=TRUE) #ensure there was no repeats
  
  
  sample=sub('results/SV/manta_INV/','',fileName)
  sample=sub('.manta.vcf','',sample)
  env=strsplit(sample,split="_")[[1]][1]
  gen=strsplit(sample,split="_")[[1]][2]
  rep=strsplit(sample,split="_")[[1]][3]
  
  sv_sample <- data.frame(row.names=c(sample))
  sv_sample$sample = sample
  sv_sample$env = env
  sv_sample$gen = gen
  sv_sample$rep = rep
  sv_sample$tot_burden = nrow(mut_data_sample)
  sv_sample$INV <- sum(str_count(mut_data_sample$SVTYPE, pattern="INV"))
  sv_sample$DEL <- sum(str_count(mut_data_sample$SVTYPE, pattern="DEL"))
  sv_sample$INS <- sum(str_count(mut_data_sample$SVTYPE, pattern="INS"))
  sv_sample$DUP <- sum(str_count(mut_data_sample$SVTYPE, pattern="DUP"))
  sv_sample$BND <- sum(str_count(mut_data_sample$SVTYPE, pattern="BND"))
  
  sv_data <- rbind(sv_data,sv_sample)
}
write.table(sv_data,file='data/Hybrid_Dynamics_SV_data.txt')
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
# PLOT TOTAL SV BURDEN
# ==============================================================================

founder <- subset(sv_data, gen %in% "Founder")

H1 <- subset(sv_data, (env == "NaCl" | env == "LiAc0.01" | env =='H1F1' | env =='H1F2'))
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
  guides(fill = guide_legend(byrow = TRUE))

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
  ylab('SV mutational load') +
  labs(tag='TOTAL')

H2 <- subset(sv_data, (env == "NaCl" | env == "LiAc0.02" | env =='H2F1' | env =='H2F2'))
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
  guides(fill = guide_legend(byrow = TRUE))
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
  ylab('SV mutational load') +
  labs(tag='TOTAL')

H3 <- subset(sv_data, (env == "NaCl" | env == "Ethanol" | env =='H3F1' | env =='H3F2'))
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
  guides(fill = guide_legend(byrow = TRUE))
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
  ylab('SV mutational load') +
  labs(tag='TOTAL')

H1_total + H1_plot + H2_total + H2_plot + H3_total + H3_plot + 
  plot_layout(ncol=2,widths = c(1.3,8))

ggsave('figures/Hybrid_Dynamics_SV_burden.pdf',width=11,height=5,dpi = 900)

# ==============================================================================
# ANOVA of SV mutational load
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
# PLOT CORRELATIONS BETWEEN GENERATIONAL TIME AND MUTATIONAL BURDEN
# ==============================================================================
hybrid <- subset(sv_data, (env =='H1F1' | env =='H1F2'|env =='H2F1' | env =='H2F2'|env =='H3F1' | env =='H3F2'))
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
  ylab('SV mutational load') +
  xlab('generations') +
  ylim(70,280)

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
  ylim(70,280)

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
  ylim(70,280)

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
  ylim(70,280)

hybrid_corr | H1_corr | H2_corr | H3_corr
ggsave('figures/Hybrid_Dynamics_SV_burden_corr.pdf',width=8,height=2.5,dpi = 900)

# ==============================================================================
# PLOT INSERTION SV BURDEN
# =============================================================================
H1 <- subset(sv_data, (env == "NaCl" | env == "LiAc0.01" | env =='H1F1' | env =='H1F2'))
H1$env <- factor(H1$env , levels=c('NaCl','LiAc0.01','H1F1','H1F2'))
H1_plot <- ggplot(H1) + 
  geom_hline(yintercept=(mean(founder$INS)),color='black',linetype='dashed',linewidth=0.7)+
  stat_summary(aes(x=as.factor(env),y=INS,color=env),fun=median, geom="point", shape=95, size=9,stroke=edgewidth,alpha=1)+
  geom_jitter(aes(x=as.factor(env),y=INS,fill=env),shape=21, colour='black',size=markersize, alpha=alpha,width=width,stroke=edgewidth) +
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
  guides(fill = guide_legend(byrow = TRUE))

H1_total <- ggplot(H1) + 
  geom_hline(yintercept=(mean(founder$INS)),color='black',linetype='dashed',linewidth=0.7)+
  stat_summary(aes(x=as.factor(env),y=INS,color=env),fun=median, geom="point", shape=95, size=12,stroke=edgewidth,alpha=1)+
  geom_jitter(aes(x=as.factor(env),y=INS,fill=env),shape=21, colour='black',size=markersize, alpha=alpha,width=width,stroke=edgewidth) +
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
  ylab('SV load | insertion') +
  labs(tag='TOTAL')

H2 <- subset(sv_data, (env == "NaCl" | env == "LiAc0.02" | env =='H2F1' | env =='H2F2'))
H2$env <- factor(H2$env , levels=c('NaCl','LiAc0.02','H2F1','H2F2'))
H2_plot <- ggplot(H2) + 
  geom_hline(yintercept=(mean(founder$INS)),color='black',linetype='dashed',linewidth=0.7)+
  stat_summary(aes(x=as.factor(env),y=INS,color=env),fun=median, geom="point", shape=95, size=9,stroke=edgewidth,alpha=1)+
  geom_jitter(aes(x=as.factor(env),y=INS,fill=env),shape=21, colour='black',size=markersize, alpha=alpha,width=width,stroke=edgewidth) +
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
  geom_hline(yintercept=(mean(founder$INS)),color='black',linetype='dashed',linewidth=0.7)+
  stat_summary(aes(x=as.factor(env),y=INS,color=env),fun=median, geom="point", shape=95, size=12,stroke=edgewidth,alpha=1)+
  geom_jitter(aes(x=as.factor(env),y=INS,fill=env),shape=21, colour='black',size=markersize, alpha=alpha,width=width,stroke=edgewidth) +
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
  ylab('SV load | insertion') +
  labs(tag='TOTAL')

H3 <- subset(sv_data, (env == "NaCl" | env == "Ethanol" | env =='H3F1' | env =='H3F2'))
H3$env <- factor(H3$env , levels=c('NaCl','Ethanol','H3F1','H3F2'))
H3_plot <- ggplot(H3) + 
  geom_hline(yintercept=(mean(founder$INS)),color='black',linetype='dashed',linewidth=0.7)+
  stat_summary(aes(x=as.factor(env),y=INS,color=env),fun=median, geom="point", shape=95, size=9,stroke=edgewidth,alpha=1)+
  geom_jitter(aes(x=as.factor(env),y=INS,fill=env),shape=21, colour='black',size=markersize, alpha=alpha,width=width,stroke=edgewidth) +
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
  geom_hline(yintercept=(mean(founder$INS)),color='black',linetype='dashed',linewidth=0.7)+
  stat_summary(aes(x=as.factor(env),y=INS,color=env),fun=median, geom="point", shape=95, size=12,stroke=edgewidth,alpha=1)+
  geom_jitter(aes(x=as.factor(env),y=INS,fill=env),shape=21, colour='black',size=markersize, alpha=alpha,width=width,stroke=edgewidth) +
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
  ylab('SV load | insertion') +
  labs(tag='TOTAL')

H1_total + H1_plot + H2_total + H2_plot + H3_total + H3_plot + 
  plot_layout(ncol=2,widths = c(1.3,8))

ggsave('figures/Hybrid_Dynamics_SV_burden-INS.pdf',width=11,height=5,dpi = 900)

# ==============================================================================
# PLOT DELETION SV BURDEN
# =============================================================================
H1 <- subset(sv_data, (env == "NaCl" | env == "LiAc0.01" | env =='H1F1' | env =='H1F2'))
H1$env <- factor(H1$env , levels=c('NaCl','LiAc0.01','H1F1','H1F2'))
H1_plot <- ggplot(H1) + 
  geom_hline(yintercept=(mean(founder$DEL)),color='black',linetype='dashed',linewidth=0.7)+
  stat_summary(aes(x=as.factor(env),y=DEL,color=env),fun=median, geom="point", shape=95, size=9,stroke=edgewidth,alpha=1)+
  geom_jitter(aes(x=as.factor(env),y=DEL,fill=env),shape=21, colour='black',size=markersize, alpha=alpha,width=width,stroke=edgewidth) +
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
  guides(fill = guide_legend(byrow = TRUE))

H1_total <- ggplot(H1) + 
  geom_hline(yintercept=(mean(founder$DEL)),color='black',linetype='dashed',linewidth=0.7)+
  stat_summary(aes(x=as.factor(env),y=DEL,color=env),fun=median, geom="point", shape=95, size=12,stroke=edgewidth,alpha=1)+
  geom_jitter(aes(x=as.factor(env),y=DEL,fill=env),shape=21, colour='black',size=markersize, alpha=alpha,width=width,stroke=edgewidth) +
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
  ylab('SV load | deletion') +
  labs(tag='TOTAL')

H2 <- subset(sv_data, (env == "NaCl" | env == "LiAc0.02" | env =='H2F1' | env =='H2F2'))
H2$env <- factor(H2$env , levels=c('NaCl','LiAc0.02','H2F1','H2F2'))
H2_plot <- ggplot(H2) + 
  geom_hline(yintercept=(mean(founder$DEL)),color='black',linetype='dashed',linewidth=0.7)+
  stat_summary(aes(x=as.factor(env),y=DEL,color=env),fun=median, geom="point", shape=95, size=9,stroke=edgewidth,alpha=1)+
  geom_jitter(aes(x=as.factor(env),y=DEL,fill=env),shape=21, colour='black',size=markersize, alpha=alpha,width=width,stroke=edgewidth) +
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
  geom_hline(yintercept=(mean(founder$DEL)),color='black',linetype='dashed',linewidth=0.7)+
  stat_summary(aes(x=as.factor(env),y=DEL,color=env),fun=median, geom="point", shape=95, size=12,stroke=edgewidth,alpha=1)+
  geom_jitter(aes(x=as.factor(env),y=DEL,fill=env),shape=21, colour='black',size=markersize, alpha=alpha,width=width,stroke=edgewidth) +
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
  ylab('SV load | deletion') +
  labs(tag='TOTAL')

H3 <- subset(sv_data, (env == "NaCl" | env == "Ethanol" | env =='H3F1' | env =='H3F2'))
H3$env <- factor(H3$env , levels=c('NaCl','Ethanol','H3F1','H3F2'))
H3_plot <- ggplot(H3) + 
  geom_hline(yintercept=(mean(founder$DEL)),color='black',linetype='dashed',linewidth=0.7)+
  stat_summary(aes(x=as.factor(env),y=DEL,color=env),fun=median, geom="point", shape=95, size=9,stroke=edgewidth,alpha=1)+
  geom_jitter(aes(x=as.factor(env),y=DEL,fill=env),shape=21, colour='black',size=markersize, alpha=alpha,width=width,stroke=edgewidth) +
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
  geom_hline(yintercept=(mean(founder$DEL)),color='black',linetype='dashed',linewidth=0.7)+
  stat_summary(aes(x=as.factor(env),y=DEL,color=env),fun=median, geom="point", shape=95, size=12,stroke=edgewidth,alpha=1)+
  geom_jitter(aes(x=as.factor(env),y=DEL,fill=env),shape=21, colour='black',size=markersize, alpha=alpha,width=width,stroke=edgewidth) +
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
  ylab('SV load | deletion') +
  labs(tag='TOTAL')

H1_total + H1_plot + H2_total + H2_plot + H3_total + H3_plot + 
  plot_layout(ncol=2,widths = c(1.3,8))

ggsave('figures/Hybrid_Dynamics_SV_burden-DEL.pdf',width=11,height=5,dpi = 900)


# ==============================================================================
# PLOT INVERSION SV BURDEN
# =============================================================================
H1 <- subset(sv_data, (env == "NaCl" | env == "LiAc0.01" | env =='H1F1' | env =='H1F2'))
H1$env <- factor(H1$env , levels=c('NaCl','LiAc0.01','H1F1','H1F2'))
H1_plot <- ggplot(H1) + 
  geom_hline(yintercept=(mean(founder$INV)),color='black',linetype='dashed',linewidth=0.7)+
  stat_summary(aes(x=as.factor(env),y=INV,color=env),fun=median, geom="point", shape=95, size=9,stroke=edgewidth,alpha=1)+
  geom_jitter(aes(x=as.factor(env),y=INV,fill=env),shape=21, colour='black',size=markersize, alpha=alpha,width=width,stroke=edgewidth) +
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
  guides(fill = guide_legend(byrow = TRUE))

H1_total <- ggplot(H1) + 
  geom_hline(yintercept=(mean(founder$INV)),color='black',linetype='dashed',linewidth=0.7)+
  stat_summary(aes(x=as.factor(env),y=INV,color=env),fun=median, geom="point", shape=95, size=12,stroke=edgewidth,alpha=1)+
  geom_jitter(aes(x=as.factor(env),y=INV,fill=env),shape=21, colour='black',size=markersize, alpha=alpha,width=width,stroke=edgewidth) +
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
  ylab('SV load | inversion') +
  labs(tag='TOTAL')

H2 <- subset(sv_data, (env == "NaCl" | env == "LiAc0.02" | env =='H2F1' | env =='H2F2'))
H2$env <- factor(H2$env , levels=c('NaCl','LiAc0.02','H2F1','H2F2'))
H2_plot <- ggplot(H2) + 
  geom_hline(yintercept=(mean(founder$INV)),color='black',linetype='dashed',linewidth=0.7)+
  stat_summary(aes(x=as.factor(env),y=INV,color=env),fun=median, geom="point", shape=95, size=9,stroke=edgewidth,alpha=1)+
  geom_jitter(aes(x=as.factor(env),y=INV,fill=env),shape=21, colour='black',size=markersize, alpha=alpha,width=width,stroke=edgewidth) +
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
  geom_hline(yintercept=(mean(founder$INV)),color='black',linetype='dashed',linewidth=0.7)+
  stat_summary(aes(x=as.factor(env),y=INV,color=env),fun=median, geom="point", shape=95, size=12,stroke=edgewidth,alpha=1)+
  geom_jitter(aes(x=as.factor(env),y=INV,fill=env),shape=21, colour='black',size=markersize, alpha=alpha,width=width,stroke=edgewidth) +
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
  ylab('SV load | inversion') +
  labs(tag='TOTAL')

H3 <- subset(sv_data, (env == "NaCl" | env == "Ethanol" | env =='H3F1' | env =='H3F2'))
H3$env <- factor(H3$env , levels=c('NaCl','Ethanol','H3F1','H3F2'))
H3_plot <- ggplot(H3) + 
  geom_hline(yintercept=(mean(founder$INV)),color='black',linetype='dashed',linewidth=0.7)+
  stat_summary(aes(x=as.factor(env),y=INV,color=env),fun=median, geom="point", shape=95, size=9,stroke=edgewidth,alpha=1)+
  geom_jitter(aes(x=as.factor(env),y=INV,fill=env),shape=21, colour='black',size=markersize, alpha=alpha,width=width,stroke=edgewidth) +
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
  geom_hline(yintercept=(mean(founder$INV)),color='black',linetype='dashed',linewidth=0.7)+
  stat_summary(aes(x=as.factor(env),y=INV,color=env),fun=median, geom="point", shape=95, size=12,stroke=edgewidth,alpha=1)+
  geom_jitter(aes(x=as.factor(env),y=INV,fill=env),shape=21, colour='black',size=markersize, alpha=alpha,width=width,stroke=edgewidth) +
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
  ylab('SV load | inversion') +
  labs(tag='TOTAL')

H1_total + H1_plot + H2_total + H2_plot + H3_total + H3_plot + 
  plot_layout(ncol=2,widths = c(1.3,8))

ggsave('figures/Hybrid_Dynamics_SV_burden-INV.pdf',width=11,height=5,dpi = 900)

# ==============================================================================
# PLOT DUPLICATION SV BURDEN
# =============================================================================
H1 <- subset(sv_data, (env == "NaCl" | env == "LiAc0.01" | env =='H1F1' | env =='H1F2'))
H1$env <- factor(H1$env , levels=c('NaCl','LiAc0.01','H1F1','H1F2'))
H1_plot <- ggplot(H1) + 
  geom_hline(yintercept=(mean(founder$DUP)),color='black',linetype='dashed',linewidth=0.7)+
  stat_summary(aes(x=as.factor(env),y=DUP,color=env),fun=median, geom="point", shape=95, size=9,stroke=edgewidth,alpha=1)+
  geom_jitter(aes(x=as.factor(env),y=DUP,fill=env),shape=21, colour='black',size=markersize, alpha=alpha,width=width,stroke=edgewidth) +
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
  guides(fill = guide_legend(byrow = TRUE))

H1_total <- ggplot(H1) + 
  geom_hline(yintercept=(mean(founder$DUP)),color='black',linetype='dashed',linewidth=0.7)+
  stat_summary(aes(x=as.factor(env),y=DUP,color=env),fun=median, geom="point", shape=95, size=12,stroke=edgewidth,alpha=1)+
  geom_jitter(aes(x=as.factor(env),y=DUP,fill=env),shape=21, colour='black',size=markersize, alpha=alpha,width=width,stroke=edgewidth) +
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
  ylab('SV load | duplication') +
  labs(tag='TOTAL')

H2 <- subset(sv_data, (env == "NaCl" | env == "LiAc0.02" | env =='H2F1' | env =='H2F2'))
H2$env <- factor(H2$env , levels=c('NaCl','LiAc0.02','H2F1','H2F2'))
H2_plot <- ggplot(H2) + 
  geom_hline(yintercept=(mean(founder$DUP)),color='black',linetype='dashed',linewidth=0.7)+
  stat_summary(aes(x=as.factor(env),y=DUP,color=env),fun=median, geom="point", shape=95, size=9,stroke=edgewidth,alpha=1)+
  geom_jitter(aes(x=as.factor(env),y=DUP,fill=env),shape=21, colour='black',size=markersize, alpha=alpha,width=width,stroke=edgewidth) +
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
  geom_hline(yintercept=(mean(founder$DUP)),color='black',linetype='dashed',linewidth=0.7)+
  stat_summary(aes(x=as.factor(env),y=DUP,color=env),fun=median, geom="point", shape=95, size=12,stroke=edgewidth,alpha=1)+
  geom_jitter(aes(x=as.factor(env),y=DUP,fill=env),shape=21, colour='black',size=markersize, alpha=alpha,width=width,stroke=edgewidth) +
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
  ylab('SV load | duplication') +
  labs(tag='TOTAL')

H3 <- subset(sv_data, (env == "NaCl" | env == "Ethanol" | env =='H3F1' | env =='H3F2'))
H3$env <- factor(H3$env , levels=c('NaCl','Ethanol','H3F1','H3F2'))
H3_plot <- ggplot(H3) + 
  geom_hline(yintercept=(mean(founder$DUP)),color='black',linetype='dashed',linewidth=0.7)+
  stat_summary(aes(x=as.factor(env),y=DUP,color=env),fun=median, geom="point", shape=95, size=9,stroke=edgewidth,alpha=1)+
  geom_jitter(aes(x=as.factor(env),y=DUP,fill=env),shape=21, colour='black',size=markersize, alpha=alpha,width=width,stroke=edgewidth) +
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
  geom_hline(yintercept=(mean(founder$DUP)),color='black',linetype='dashed',linewidth=0.7)+
  stat_summary(aes(x=as.factor(env),y=DUP,color=env),fun=median, geom="point", shape=95, size=12,stroke=edgewidth,alpha=1)+
  geom_jitter(aes(x=as.factor(env),y=DUP,fill=env),shape=21, colour='black',size=markersize, alpha=alpha,width=width,stroke=edgewidth) +
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
  ylab('SV load | duplication') +
  labs(tag='TOTAL')

H1_total + H1_plot + H2_total + H2_plot + H3_total + H3_plot + 
  plot_layout(ncol=2,widths = c(1.3,8))

ggsave('figures/Hybrid_Dynamics_SV_burden-DUP.pdf',width=11,height=5,dpi = 900)

# ==============================================================================
# PLOT BREAKEND SV BURDEN
# =============================================================================
H1 <- subset(sv_data, (env == "NaCl" | env == "LiAc0.01" | env =='H1F1' | env =='H1F2'))
H1$env <- factor(H1$env , levels=c('NaCl','LiAc0.01','H1F1','H1F2'))
H1_plot <- ggplot(H1) + 
  geom_hline(yintercept=(mean(founder$BND)),color='black',linetype='dashed',linewidth=0.7)+
  stat_summary(aes(x=as.factor(env),y=BND,color=env),fun=median, geom="point", shape=95, size=9,stroke=edgewidth,alpha=1)+
  geom_jitter(aes(x=as.factor(env),y=BND,fill=env),shape=21, colour='black',size=markersize, alpha=alpha,width=width,stroke=edgewidth) +
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
  guides(fill = guide_legend(byrow = TRUE))

H1_total <- ggplot(H1) + 
  geom_hline(yintercept=(mean(founder$BND)),color='black',linetype='dashed',linewidth=0.7)+
  stat_summary(aes(x=as.factor(env),y=BND,color=env),fun=median, geom="point", shape=95, size=12,stroke=edgewidth,alpha=1)+
  geom_jitter(aes(x=as.factor(env),y=BND,fill=env),shape=21, colour='black',size=markersize, alpha=alpha,width=width,stroke=edgewidth) +
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
  ylab('SV load | breakend') +
  labs(tag='TOTAL')

H2 <- subset(sv_data, (env == "NaCl" | env == "LiAc0.02" | env =='H2F1' | env =='H2F2'))
H2$env <- factor(H2$env , levels=c('NaCl','LiAc0.02','H2F1','H2F2'))
H2_plot <- ggplot(H2) + 
  geom_hline(yintercept=(mean(founder$BND)),color='black',linetype='dashed',linewidth=0.7)+
  stat_summary(aes(x=as.factor(env),y=BND,color=env),fun=median, geom="point", shape=95, size=9,stroke=edgewidth,alpha=1)+
  geom_jitter(aes(x=as.factor(env),y=BND,fill=env),shape=21, colour='black',size=markersize, alpha=alpha,width=width,stroke=edgewidth) +
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
  geom_hline(yintercept=(mean(founder$BND)),color='black',linetype='dashed',linewidth=0.7)+
  stat_summary(aes(x=as.factor(env),y=BND,color=env),fun=median, geom="point", shape=95, size=12,stroke=edgewidth,alpha=1)+
  geom_jitter(aes(x=as.factor(env),y=BND,fill=env),shape=21, colour='black',size=markersize, alpha=alpha,width=width,stroke=edgewidth) +
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
  ylab('SV load | breakend') +
  labs(tag='TOTAL')

H3 <- subset(sv_data, (env == "NaCl" | env == "Ethanol" | env =='H3F1' | env =='H3F2'))
H3$env <- factor(H3$env , levels=c('NaCl','Ethanol','H3F1','H3F2'))
H3_plot <- ggplot(H3) + 
  geom_hline(yintercept=(mean(founder$BND)),color='black',linetype='dashed',linewidth=0.7)+
  stat_summary(aes(x=as.factor(env),y=BND,color=env),fun=median, geom="point", shape=95, size=9,stroke=edgewidth,alpha=1)+
  geom_jitter(aes(x=as.factor(env),y=BND,fill=env),shape=21, colour='black',size=markersize, alpha=alpha,width=width,stroke=edgewidth) +
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
  geom_hline(yintercept=(mean(founder$BND)),color='black',linetype='dashed',linewidth=0.7)+
  stat_summary(aes(x=as.factor(env),y=BND,color=env),fun=median, geom="point", shape=95, size=12,stroke=edgewidth,alpha=1)+
  geom_jitter(aes(x=as.factor(env),y=BND,fill=env),shape=21, colour='black',size=markersize, alpha=alpha,width=width,stroke=edgewidth) +
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
  ylab('SV load | breakend') +
  labs(tag='TOTAL')

H1_total + H1_plot + H2_total + H2_plot + H3_total + H3_plot + 
  plot_layout(ncol=2,widths = c(1.3,8))

ggsave('figures/Hybrid_Dynamics_SV_burden-BND.pdf',width=11,height=5,dpi = 900)