#"`-:-.   ,-;"`-:-.   ,-;"`-:-.   ,-;"`-:-. -------------------------- ,-;"`-:-.
# `=`,'=/     `=`,'=/     `=`,'=/     `=`|  "If I have seen further |=/     `=
#   y==/        y==/        y==/        y| it is by standing on the |/
# ,=,-<=`.    ,=,-<=`.    ,=,-<=`.    ,=,|    shoulders of Giants"  |=`.    ,=
#,-'-'   `-=_,-'-'   `-=_,-'-'   `-=_,-'-'  -------------------------- `-=_,-'-'
# ==============================================================================
# PROJECT: Hybrid_Dynamics
# PURPOSE: to plot the FST calculations
# DATE: 1 Jan 2024
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
library(stringr)
library(reshape)
# ==============================================================================
# 
# ==============================================================================
H1_fst <- read.table(file='results/SNV/FST/H1_Fst.txt',row.names = NULL)
H2_fst <- read.table(file='results/SNV/FST/H2_Fst.txt',row.names = NULL)
H3_fst <- read.table(file='results/SNV/FST/H3_Fst.txt',row.names = NULL)
G200_400_fst <- read.table(file='results/SNV/FST/G200_G400_Fst.txt',row.names = NULL)
G1000_fst <- read.table(file='results/SNV/FST/G1000_Fst.txt',row.names = NULL)

fst <- rbind(H1_fst,H2_fst,H3_fst,G200_400_fst,G1000_fst)
fst <- separate_wider_delim(fst, cols = row.names, delim = ";", names = c("P1", "P2"))


fst$P1 <- str_replace_all(fst$P1, "N_Founder", "N_Founder_R1")
fst$P2 <- str_replace_all(fst$P2, "N_Founder", "N_Founder_R1")
fst$P2 <- str_replace_all(fst$P2, "LE_Founder", "LE_Founder_R1")
fst$P1 <- str_replace_all(fst$P1, "LE_Founder", "LE_Founder_R1")

fst <- separate_wider_delim(fst, cols = P1, delim = "_", names = c("p1env", "p1gen",'p1rep'),cols_remove = FALSE)
fst <- separate_wider_delim(fst, cols = P2, delim = "_", names = c("p2env", "p2gen",'p2rep'),cols_remove = FALSE)

fst$env <- ifelse(fst$p1env == fst$p2env, 'within', 'between')
fst$p1gen <- ifelse(fst$p1gen == 'Founder','G0',fst$p1gen)
fst$p2gen <- ifelse(fst$p2gen == 'Founder','G0',fst$p2gen)
fst$p1_gen_num <- as.numeric(sub('G','',fst$p1gen))
fst$p2_gen_num <- as.numeric(sub('G','',fst$p2gen))

fst <- subset(fst,p1gen==p2gen)
fst$p1_gen_num <- as.factor(fst$p1_gen_num)

H1_fst <- dplyr::filter(fst, grepl(paste(c('NaCl','LiAc0.01'),collapse="|"), p1env))
H1_fst <- dplyr::filter(H1_fst, grepl(paste(c('NaCl','LiAc0.01'),collapse="|"), p2env))
H1_fst$env_final <- ifelse(H1_fst$env=='within' & H1_fst$p1env == 'NaCl','NaCl',ifelse(H1_fst$env=='within' & H1_fst$p1env == 'LiAc0.01','LiAc0.01','Divergent'))

H2_fst <- dplyr::filter(fst, grepl(paste(c('NaCl','LiAc0.02'),collapse="|"), p1env))
H2_fst <- dplyr::filter(H2_fst, grepl(paste(c('NaCl','LiAc0.02'),collapse="|"), p2env))
H2_fst$env_final <- ifelse(H2_fst$env=='within' & H2_fst$p1env == 'NaCl','NaCl',ifelse(H2_fst$env=='within' & H2_fst$p1env == 'LiAc0.02','LiAc0.02','Divergent'))

H3_fst <- dplyr::filter(fst, grepl(paste(c('NaCl','Ethanol'),collapse="|"), p1env))
H3_fst <- dplyr::filter(H3_fst, grepl(paste(c('NaCl','Ethanol'),collapse="|"), p2env))
H3_fst$env_final <- ifelse(H3_fst$env=='within' & H3_fst$p1env == 'NaCl','NaCl',ifelse(H3_fst$env=='within' & H3_fst$p1env == 'Ethanol','Ethanol','Divergent'))


# ==============================================================================
# Environment Colors
# ==============================================================================
NaCl = '#4472c4ff'
LiAc0.01 = '#9e49e1ff'
LiAc0.02 = '#00b050ff'
Ethanol = '#ff0000ff'

# ==============================================================================
# 
# ==============================================================================


H1_overall<- ggplot(H1_fst, aes(x=env_final, Fst.Estimate)) +
  geom_hline(yintercept=0.01770595, linetype="dashed", 
             color = "black", size=0.5) +  
  geom_boxplot(aes(fill=env_final),outlier.alpha = 0.5,show.legend=FALSE)+
  scale_fill_manual(values=c('dimgrey',LiAc0.01,NaCl))+
  theme(
    panel.background = element_rect(fill = "lightgrey"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(colour = "white", linetype = 1, linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    panel.spacing = unit(0.5, "lines"), 
    axis.title.x=element_blank(),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  labs(y='Fst Estimate')+
  guides(fill=guide_legend(title="Environment"))+
  ylim(0,0.51)

H1_scatter<- ggplot(H1_fst, aes(x=p2_gen_num, metric)) +
  geom_hline(yintercept=0.01770595, linetype="dashed",  
             color = "black", size=0.5) + 
  stat_summary(aes(y=Fst.Estimate, group=env_final,color=env_final),fun=mean, geom="line",alpha=1,,show.legend=FALSE)+
  stat_summary(aes(x=p2_gen_num,y=Fst.Estimate,group=env_final),fun.data = mean_se,  geom = "errorbar",width=0,show.legend=FALSE) + 
  stat_summary(aes(x=p2_gen_num,y=Fst.Estimate,fill=env_final),fun=mean, geom="point", shape=21, size=2,stroke=0.7,alpha=1,show.legend=FALSE)+
  scale_fill_manual(values=c('dimgrey',LiAc0.01,NaCl))+
  scale_colour_manual(values=c('dimgrey',LiAc0.01,NaCl))+
  theme(
    panel.background = element_rect(fill = "white"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(colour = "lightgrey", linetype = 1, linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    axis.text.y=element_blank(),
    axis.title.y=element_blank(),
    panel.spacing = unit(0.5, "lines"), 
    legend.position = "bottom",
    axis.text.x = element_text())+
  labs(x='Generations')+
  guides(fill=guide_legend(title="Environment"))+
  ylim(0,0.51)


H2_overall<- ggplot(H2_fst, aes(x=env_final, Fst.Estimate)) +
  geom_hline(yintercept=0.01770595, linetype="dashed",  
             color = "black", size=0.5) +  
  geom_boxplot(aes(fill=env_final),outlier.alpha = 0.5,show.legend=FALSE)+
  scale_fill_manual(values=c('dimgrey',LiAc0.02,NaCl))+
  theme(
    panel.background = element_rect(fill = "lightgrey"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(colour = "white", linetype = 1, linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    panel.spacing = unit(0.5, "lines"), 
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  guides(fill=guide_legend(title="Environment"))+
  ylim(0,0.51)

H2_scatter<- ggplot(H2_fst, aes(x=p2_gen_num, metric)) +
  geom_hline(yintercept=0.01770595, linetype="dashed", 
             color = "black", size=0.5) + 
  stat_summary(aes(y=Fst.Estimate, group=env_final,color=env_final),fun=mean, geom="line",alpha=1,,show.legend=FALSE)+
  stat_summary(aes(x=p2_gen_num,y=Fst.Estimate,group=env_final),fun.data = mean_se,  geom = "errorbar",width=0,show.legend=FALSE) + 
  stat_summary(aes(x=p2_gen_num,y=Fst.Estimate,fill=env_final),fun=mean, geom="point", shape=21, size=2,stroke=0.7,alpha=1,show.legend=FALSE)+
  scale_fill_manual(values=c('dimgrey',LiAc0.02,NaCl))+
  scale_colour_manual(values=c('dimgrey',LiAc0.02,NaCl))+
  theme(
    panel.background = element_rect(fill = "white"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(colour = "lightgrey", linetype = 1, linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    axis.text.y=element_blank(),
    axis.title.y=element_blank(),
    panel.spacing = unit(0.5, "lines"), 
    legend.position = "bottom",
    axis.text.x = element_text())+
  labs(x='Generations')+
  guides(fill=guide_legend(title="Environment"))+
  ylim(0,0.51)


H3_overall<- ggplot(H3_fst, aes(x=env_final, Fst.Estimate)) +
  geom_hline(yintercept=0.01770595, linetype="dashed", 
             color = "black", size=0.5) +  
  geom_boxplot(aes(fill=env_final),outlier.alpha = 0.5,show.legend=FALSE)+
  scale_fill_manual(values=c('dimgrey',Ethanol,NaCl))+
  theme(
    panel.background = element_rect(fill = "lightgrey"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(colour = "white", linetype = 1, linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    panel.spacing = unit(0.5, "lines"), 
    axis.title.x=element_blank(),
    axis.title.y=element_blank(),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  guides(fill=guide_legend(title="Environment"))+
  ylim(0,0.51)

H3_scatter<- ggplot(H3_fst, aes(x=p2_gen_num, metric)) +
  geom_hline(yintercept=0.01770595, linetype="dashed", 
             color = "black", size=0.5) + 
  stat_summary(aes(y=Fst.Estimate, group=env_final,color=env_final),fun=mean, geom="line",alpha=1,,show.legend=FALSE)+
  stat_summary(aes(x=p2_gen_num,y=Fst.Estimate,group=env_final),fun.data = mean_se,  geom = "errorbar",width=0,show.legend=FALSE) + 
  stat_summary(aes(x=p2_gen_num,y=Fst.Estimate,fill=env_final),fun=mean, geom="point", shape=21, size=2,stroke=0.7,alpha=1,show.legend=FALSE)+
  scale_fill_manual(values=c('dimgrey',Ethanol,NaCl))+
  scale_colour_manual(values=c('dimgrey',Ethanol,NaCl))+
  theme(
    panel.background = element_rect(fill = "white"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.y = element_line(colour = "lightgrey", linetype = 1, linewidth = 0.5),
    panel.grid.minor.y = element_blank(),
    axis.text.y=element_blank(),
    axis.title.y=element_blank(),
    panel.spacing = unit(0.5, "lines"), 
    legend.position = "bottom",
    axis.text.x = element_text())+
  labs(x='Generations')+
  guides(fill=guide_legend(title="Environment"))+
  ylim(0,0.51)


H1_overall + H1_scatter +H2_overall + H2_scatter +H3_overall + H3_scatter +
  plot_layout(widths=c(2,5,2,5,2,5)) 

ggsave('figures/Hybrid_Dynamics_Fst-gen.pdf',width=7.5,height=2,dpi = 900)


# ==============================================================================
# Pairwise Kruskal-Wallis
# ==============================================================================
kruskal.test(Fst.Estimate ~ env_final, data = H1_fst)
pairwise.wilcox.test(H1_fst$Fst.Estimate, H1_fst$env_final,
                     p.adjust.method = "BH")
kruskal.test(Fst.Estimate ~ env_final, data = H2_fst)
pairwise.wilcox.test(H2_fst$Fst.Estimate, H2_fst$env_final,
                     p.adjust.method = "BH")
kruskal.test(Fst.Estimate ~ env_final, data = H3_fst)
pairwise.wilcox.test(H3_fst$Fst.Estimate, H3_fst$env_final,
                     p.adjust.method = "BH")
# ==============================================================================
# 
# ==============================================================================
library(multcomp)
H1_fst$env_final <- as.factor(H1_fst$env_final)
res_aov <- aov(Fst.Estimate ~ env_final,
               data = H1_fst)
summary(res_aov)
post_test <- glht(res_aov,
                  linfct = mcp(env_final = "Tukey"))
summary(post_test)




H2_fst$env_final <- as.factor(H2_fst$env_final)
res_aov <- aov(Fst.Estimate ~ env_final,
               data = H2_fst)
summary(res_aov)
post_test <- glht(res_aov,
                  linfct = mcp(env_final = "Tukey"))
summary(post_test)

H3_fst$env_final <- as.factor(H3_fst$env_final)
res_aov <- aov(Fst.Estimate ~ env_final,
               data = H3_fst)
summary(res_aov)
post_test <- glht(res_aov,
                  linfct = mcp(env_final = "Tukey"))
summary(post_test)
