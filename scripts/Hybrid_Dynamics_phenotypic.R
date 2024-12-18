#"`-:-.   ,-;"`-:-.   ,-;"`-:-.   ,-;"`-:-. -------------------------- ,-;"`-:-.
# `=`,'=/     `=`,'=/     `=`,'=/     `=`|  "If I have seen further |=/     `=
#   y==/        y==/        y==/        y| it is by standing on the |/
# ,=,-<=`.    ,=,-<=`.    ,=,-<=`.    ,=,|    shoulders of Giants"  |=`.    ,=
#,-'-'   `-=_,-'-'   `-=_,-'-'   `-=_,-'-'  -------------------------- `-=_,-'-'
# ==============================================================================
# PROJECT: Hybrid Dynamics
# PURPOSE: to plot the phenotypic distances between parental strains
# DATE: May 3rd, 2024
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
#quartz()
library(stringr)
library(reshape)

# ==============================================================================
# 
# ==============================================================================
pheno <- read.table(file='data/final_phenotypic.txt',header=TRUE)
pheno$env <- ifelse(pheno$P1 == pheno$P2, 'within', 'between')
founder <- dplyr::filter(pheno,Generation==0)
pheno <- dplyr::filter(pheno,Generation!=0)
pheno$Generation <- as.factor(pheno$Generation)

H1_pheno <- dplyr::filter(pheno, grepl(paste(c('NaCl','LiAc0.01'),collapse="|"), P1))
H1_pheno <- dplyr::filter(H1_pheno, grepl(paste(c('NaCl','LiAc0.01'),collapse="|"), P2))
H1_pheno$env_final <- ifelse(H1_pheno$env=='within' & H1_pheno$P1 == 'NaCl','NaCl',ifelse(H1_pheno$env=='within' & H1_pheno$P1 == 'LiAc0.01','LiAc0.01','divergent'))

H2_pheno <- dplyr::filter(pheno, grepl(paste(c('NaCl','LiAc0.02'),collapse="|"), P1))
H2_pheno <- dplyr::filter(H2_pheno, grepl(paste(c('NaCl','LiAc0.02'),collapse="|"), P2))
H2_pheno$env_final <- ifelse(H2_pheno$env=='within' & H2_pheno$P1 == 'NaCl','NaCl',ifelse(H2_pheno$env=='within' & H2_pheno$P1 == 'LiAc0.02','LiAc0.02','divergent'))

H3_pheno <- dplyr::filter(pheno, grepl(paste(c('NaCl','EtOH'),collapse="|"), P1))
H3_pheno <- dplyr::filter(H3_pheno, grepl(paste(c('NaCl','EtOH'),collapse="|"), P2))
H3_pheno$env_final <- ifelse(H3_pheno$env=='within' & H3_pheno$P1 == 'NaCl','NaCl',ifelse(H3_pheno$env=='within' & H3_pheno$P1 == 'EtOH','Ethanol','divergent'))


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

H1_overall<- ggplot(H1_pheno, aes(x=env_final, Prop)) +
  geom_hline(yintercept=0.0222, linetype="dashed", 
             color = "black", linewidth=0.5) +  
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
    text=element_text(family="EB Garamond"),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 60, vjust = 1, hjust=1))+
  labs(y='Parental\nPhenotypic Divergence')+
  guides(fill=guide_legend(title="Environment"))+
  ylim(0,1)

H1_scatter<- ggplot(H1_pheno, aes(x=Generation, Prop)) +
  geom_hline(yintercept=0.0222, linetype="dashed",  
             color = "black", linewidth=0.5) +
  stat_summary(aes(y=Prop, group=env_final,color=env_final),fun=mean, geom="line",alpha=1,,show.legend=FALSE)+
  stat_summary(aes(x=Generation,y=Prop,group=env_final),fun.data = mean_se,  geom = "errorbar",width=0,show.legend=FALSE) + 
  stat_summary(aes(x=Generation,y=Prop,fill=env_final),fun=mean, geom="point", shape=21, size=2,stroke=0.7,alpha=1,,show.legend=FALSE)+
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
    text=element_text(family="EB Garamond"),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 60, vjust = 1, hjust=1))+
  labs(x='generations')+
  guides(fill=guide_legend(title="Environment"))+
  ylim(0,1)

H2_overall<- ggplot(H2_pheno, aes(x=env_final, Prop)) +
  geom_hline(yintercept=0.0222, linetype="dashed", 
             color = "black", linewidth=0.5) +  
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
    text=element_text(family="EB Garamond"),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 60, vjust = 1, hjust=1))+
  guides(fill=guide_legend(title="Environment"))+
  ylim(0,1)

H2_scatter<- ggplot(H2_pheno, aes(x=Generation, metric),) +
  geom_hline(yintercept=0.0222, linetype="dashed",  
             color = "black", linewidth=0.5) +
  stat_summary(aes(y=Prop, group=env_final,color=env_final),fun=mean, geom="line",alpha=1,,show.legend=FALSE)+
  stat_summary(aes(x=Generation,y=Prop,group=env_final),fun.data = mean_se,  geom = "errorbar",width=0,show.legend=FALSE) + 
  stat_summary(aes(x=Generation,y=Prop,fill=env_final),fun=mean, geom="point", shape=21, size=2,stroke=0.7,alpha=1,show.legend=FALSE)+
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
    text=element_text(family="EB Garamond"),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 60, vjust = 1, hjust=1))+
  labs(x='generations')+
  guides(fill=guide_legend(title="Environment"))+
  ylim(0,1)

H3_overall<- ggplot(H3_pheno, aes(x=env_final, Prop)) +
  geom_hline(yintercept=0.0222, linetype="dashed", 
             color = "black", linewidth=0.5) +  
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
    text=element_text(family="EB Garamond"),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 60, vjust = 1, hjust=1))+
  guides(fill=guide_legend(title="Environment"))+
  ylim(0,1)

H3_scatter<- ggplot(H3_pheno, aes(x=Generation, Prop)) +
  geom_hline(yintercept=0.0222, linetype="dashed",  
             color = "black", linewidth=0.5) +
  stat_summary(aes(y=Prop, group=env_final,color=env_final),fun=mean, geom="line",alpha=1,,show.legend=FALSE)+
  stat_summary(aes(x=Generation,y=Prop,group=env_final),fun.data = mean_se,  geom = "errorbar",width=0,show.legend=FALSE) + 
  stat_summary(aes(x=Generation,y=Prop,fill=env_final),fun=mean, geom="point", shape=21, size=2,stroke=0.7,alpha=1,show.legend=FALSE)+
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
    text=element_text(family="EB Garamond"),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 60, vjust = 1, hjust=1))+
  labs(x='generations')+
  guides(fill=guide_legend(title="Environment"))+
  ylim(0,1)

H1_overall + H1_scatter +H2_overall + H2_scatter +H3_overall + H3_scatter +
  plot_layout(guides = 'collect',widths=c(1.5,5,1.5,5,1.5,5)) 

ggsave('figures/Hybrid_Dynamics_phenotypic-gen.pdf',width=10,height=2.5,dpi = 900)


# ==============================================================================
# 
# ==============================================================================
library(multcomp)
H1_pheno$env_final <- as.factor(H1_pheno$env_final)
res_aov <- aov(Prop ~ env_final,
               data = H1_pheno)
summary(res_aov)
post_test <- glht(res_aov,
                  linfct = mcp(env_final = "Tukey"))
summary(post_test)

H2_pheno$env_final <- as.factor(H2_pheno$env_final)
res_aov <- aov(Prop ~ env_final,
               data = H2_pheno)
summary(res_aov)
post_test <- glht(res_aov,
                  linfct = mcp(env_final = "Tukey"))
summary(post_test)

H3_pheno$env_final <- as.factor(H3_pheno$env_final)
res_aov <- aov(Prop ~ env_final,
               data = H3_pheno)
summary(res_aov)
post_test <- glht(res_aov,
                  linfct = mcp(env_final = "Tukey"))
summary(post_test)