#"`-:-.   ,-;"`-:-.   ,-;"`-:-.   ,-;"`-:-. -------------------------- ,-;"`-:-.
# `=`,'=/     `=`,'=/     `=`,'=/     `=`|  "If I have seen further |=/     `=
#   y==/        y==/        y==/        y| it is by standing on the |/
# ,=,-<=`.    ,=,-<=`.    ,=,-<=`.    ,=,|    shoulders of Giants"  |=`.    ,=
#,-'-'   `-=_,-'-'   `-=_,-'-'   `-=_,-'-'  -------------------------- `-=_,-'-'
# ==============================================================================
# PROJECT: Hybrid_Dynamics
# PURPOSE: to determine the similarity between copy number profiles
# DATE: 22 Apr 2024
# Devin P. Bendixsen, PhD
# Staff Bioinformatician | University of Edinburgh
# MRC Human Genetics Unit | Institute of Genetics and Cancer
# ==============================================================================

# ==============================================================================
# LOAD NEEDED MODULES
# ==============================================================================
library(GenomicRanges)
BiocManager::install("CNVMetrics")
library(CNVMetrics)
library(dplyr)
library(showtext)
library(patchwork)
font_add_google("EB Garamond",family="Garamond")
showtext_auto()
quartz()

fileNames <- Sys.glob("results/FREEC/*_CNVs.p.value.txt") # identify all samples
CNV_calls <- data.frame()
for (fileName in fileNames) {
  print(fileName)
  cnv_sample <- read.table(file=fileName,header=FALSE,skip = 1,col.names=c('chr','start','end','CN','type','Wilcoxon','KS'))
  cnv_sample$length <- (cnv_sample$end - cnv_sample$start)
  cnv_sample <- subset(cnv_sample,Wilcoxon<0.05 & KS <0.05)
  
  cnv_sample <- subset(cnv_sample,select=c(chr,start,end,type))
  
  sample_name <- sub("results/FREEC/","",fileName) # extract the sample name from the file pathway
  sample_name <- sub(".mapped.sort.picard.bam_CNVs.p.value.txt","",sample_name)
  
  cnv_sample$ID <- sample_name
  CNV_calls <- rbind(CNV_calls,cnv_sample)
}

CNV_calls[CNV_calls == "gain"] <- "AMPLIFICATION"
CNV_calls[CNV_calls == "loss"] <- "DELETION"

CNV_calls <- CNV_calls %>%
  dplyr::rename(state = type)


H1 <- dplyr::filter(CNV_calls, grepl(paste(c('NaCl','LiAc0.01'),collapse="|"), ID))
grl <- GenomicRanges::makeGRangesListFromDataFrame(H1, 
                                                   split.field="ID", keep.extra.columns=TRUE)
H1_metric <- calculateOverlapMetric(segmentData=grl, method="jaccard", nJobs=1)

plotMetric(H1_metric, type="AMPLIFICATION",
                 colorRange=c("white", "#9B1D20"),
                 show_colnames=TRUE, fontfamily='serif', fontsize=8,
                 filename='figures/Hybrid_Dynamics_H1_Jaccard-AMP.pdf')

plotMetric(H1_metric, type="DELETION",
                 colorRange=c("white", "#0B4F6C"),
                 show_colnames=TRUE, fontfamily='serif',fontsize=8,
                 filename='figures/Hybrid_Dynamics_H1_Jaccard-DEL.pdf')


H2 <- dplyr::filter(CNV_calls, grepl(paste(c('NaCl','LiAc0.02'),collapse="|"), ID))
grl <- GenomicRanges::makeGRangesListFromDataFrame(H2, 
                                                   split.field="ID", keep.extra.columns=TRUE)
H2_metric <- calculateOverlapMetric(segmentData=grl, method="jaccard", nJobs=1)

plotMetric(H2_metric, type="AMPLIFICATION",
           colorRange=c("white", "#9B1D20"),
           show_colnames=TRUE, fontfamily='serif', fontsize=8,
           filename='figures/Hybrid_Dynamics_H2_Jaccard-AMP.pdf')

plotMetric(H2_metric, type="DELETION",
           colorRange=c("white", "#0B4F6C"),
           show_colnames=TRUE, fontfamily='serif',fontsize=8,
           filename='figures/Hybrid_Dynamics_H2_Jaccard-DEL.pdf')


H3 <- dplyr::filter(CNV_calls, grepl(paste(c('NaCl','Ethanol'),collapse="|"), ID))
grl <- GenomicRanges::makeGRangesListFromDataFrame(H3, 
                                                   split.field="ID", keep.extra.columns=TRUE)
H3_metric <- calculateOverlapMetric(segmentData=grl, method="jaccard", nJobs=1)

plotMetric(H3_metric, type="AMPLIFICATION",
           colorRange=c("white", "#9B1D20"),
           show_colnames=TRUE, fontfamily='serif', fontsize=8,
           filename='figures/Hybrid_Dynamics_H3_Jaccard-AMP.pdf')

plotMetric(H3_metric, type="DELETION",
           colorRange=c("white", "#0B4F6C"),
           show_colnames=TRUE, fontfamily='serif',fontsize=8,
           filename='figures/Hybrid_Dynamics_H3_Jaccard-DEL.pdf')



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
# H1
# ==============================================================================

H1_amp <- as.data.frame(
  as.table(H1_metric$AMPLIFICATION),
  responseName = 'metric'
)

H1_amp <- na.omit(H1_amp)
H1_amp <- H1_amp %>%
  dplyr::rename(p1 = Var1,p2=Var2)

H1_amp$p1 <- str_replace_all(H1_amp$p1, "N_Founder", "N_Founder_R1")
H1_amp$p2 <- str_replace_all(H1_amp$p2, "N_Founder", "N_Founder_R1")
H1_amp$p2 <- str_replace_all(H1_amp$p2, "LE_Founder", "LE_Founder_R1")
H1_amp$p1 <- str_replace_all(H1_amp$p1, "LE_Founder", "LE_Founder_R1")
H1_amp <- separate_wider_delim(H1_amp, cols = p1, delim = "_", names = c("p1env", "p1gen",'p1rep'),cols_remove = FALSE)
H1_amp <- separate_wider_delim(H1_amp, cols = p2, delim = "_", names = c("p2env", "p2gen",'p2rep'),cols_remove = FALSE)
H1_amp$env <- ifelse(H1_amp$p1env == H1_amp$p2env, 'within', 'between')
H1_amp$p1_gen_num <- as.numeric(sub('G','',H1_amp$p1gen))
H1_amp$p2_gen_num <- as.numeric(sub('G','',H1_amp$p2gen))
H1_amp$dist <- ifelse(H1_amp$p1env == H1_amp$p2env & H1_amp$p1rep == H1_amp$p2rep, abs(H1_amp$p1_gen_num - H1_amp$p2_gen_num), H1_amp$p1_gen_num + H1_amp$p2_gen_num)
H1_amp$env[H1_amp$dist<=100]='within'

H1_amp <- subset(H1_amp,p1gen==p2gen)
H1_amp$dist <- H1_amp$p1_gen_num
#H1_amp$dist <- as.factor(H1_amp$dist)

H1_amp$env_final <- ifelse(H1_amp$env=='within' & H1_amp$p1env == 'NaCl','NaCl',ifelse(H1_amp$env=='within' & H1_amp$p1env == 'LiAc0.01','LiAc0.01','Divergent'))


H1_overall<- ggplot(H1_amp, aes(x=env_final, metric)) +
  geom_hline(yintercept=0.21147270, linetype="dashed", 
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
  labs(y='Jaccard similarity')+
  guides(fill=guide_legend(title="Environment"))+
  ylim(0,0.82)

H1_scatter<- ggplot(H1_amp, aes(x=dist, metric)) +
  geom_hline(yintercept=0.21147270, linetype="dashed", 
             color = "black", size=0.5) + 
  stat_summary(aes(y=metric, group=env_final,color=env_final),fun=mean, geom="line",alpha=1,,show.legend=FALSE)+
  stat_summary(aes(x=dist,y=metric,group=env_final),fun.data = mean_se,  geom = "errorbar",width=0,show.legend=FALSE) + 
  stat_summary(aes(x=dist,y=metric,fill=env_final),fun=mean, geom="point", shape=21, size=2,stroke=0.7,alpha=1,show.legend=FALSE)+
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
  ylim(0,0.82)


# ==============================================================================
# H2
# ==============================================================================

H2_amp <- as.data.frame(
  as.table(H2_metric$AMPLIFICATION),
  responseName = 'metric'
)

H2_amp <- na.omit(H2_amp)
H2_amp <- H2_amp %>%
  dplyr::rename(p1 = Var1,p2=Var2)

H2_amp$p1 <- str_replace_all(H2_amp$p1, "N_Founder", "N_Founder_R1")
H2_amp$p2 <- str_replace_all(H2_amp$p2, "N_Founder", "N_Founder_R1")
H2_amp$p2 <- str_replace_all(H2_amp$p2, "LE_Founder", "LE_Founder_R1")
H2_amp$p1 <- str_replace_all(H2_amp$p1, "LE_Founder", "LE_Founder_R1")
H2_amp <- separate_wider_delim(H2_amp, cols = p1, delim = "_", names = c("p1env", "p1gen",'p1rep'),cols_remove = FALSE)
H2_amp <- separate_wider_delim(H2_amp, cols = p2, delim = "_", names = c("p2env", "p2gen",'p2rep'),cols_remove = FALSE)
H2_amp$env <- ifelse(H2_amp$p1env == H2_amp$p2env, 'within', 'between')
H2_amp$p1_gen_num <- as.numeric(sub('G','',H2_amp$p1gen))
H2_amp$p2_gen_num <- as.numeric(sub('G','',H2_amp$p2gen))
H2_amp$dist <- ifelse(H2_amp$p1env == H2_amp$p2env & H2_amp$p1rep == H2_amp$p2rep, abs(H2_amp$p1_gen_num - H2_amp$p2_gen_num), H2_amp$p1_gen_num + H2_amp$p2_gen_num)
H2_amp$env[H2_amp$dist<=100]='within'

H2_amp <- subset(H2_amp,p1gen==p2gen)
H2_amp$dist <- H2_amp$p1_gen_num
#H2_amp$dist <- as.factor(H2_amp$dist)

H2_amp$env_final <- ifelse(H2_amp$env=='within' & H2_amp$p1env == 'NaCl','NaCl',ifelse(H2_amp$env=='within' & H2_amp$p1env == 'LiAc0.02','LiAc0.02','Divergent'))


H2_overall<- ggplot(H2_amp, aes(x=env_final, metric)) +
  geom_hline(yintercept=0.21147270, linetype="dashed", 
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
    axis.title.y=element_blank(),
    axis.title.x=element_blank(),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  guides(fill=guide_legend(title="Environment"))+
  ylim(0,0.82)

H2_scatter<- ggplot(H2_amp, aes(x=dist, metric)) +
  geom_hline(yintercept=0.21147270, linetype="dashed", 
             color = "black", size=0.5) +  
  stat_summary(aes(y=metric, group=env_final,color=env_final),fun=mean, geom="line",alpha=1,,show.legend=FALSE)+
  stat_summary(aes(x=dist,y=metric,group=env_final),fun.data = mean_se,  geom = "errorbar",width=0,show.legend=FALSE) + 
  stat_summary(aes(x=dist,y=metric,fill=env_final),fun=mean, geom="point", shape=21, size=2,stroke=0.7,alpha=1,show.legend=FALSE)+
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
  ylim(0,0.82)


# ==============================================================================
# H3
# ==============================================================================

H3_amp <- as.data.frame(
  as.table(H3_metric$AMPLIFICATION),
  responseName = 'metric'
)

H3_amp <- na.omit(H3_amp)
H3_amp <- H3_amp %>%
  dplyr::rename(p1 = Var1,p2=Var2)

H3_amp$p1 <- str_replace_all(H3_amp$p1, "N_Founder", "N_Founder_R1")
H3_amp$p2 <- str_replace_all(H3_amp$p2, "N_Founder", "N_Founder_R1")
H3_amp$p2 <- str_replace_all(H3_amp$p2, "LE_Founder", "LE_Founder_R1")
H3_amp$p1 <- str_replace_all(H3_amp$p1, "LE_Founder", "LE_Founder_R1")
H3_amp <- separate_wider_delim(H3_amp, cols = p1, delim = "_", names = c("p1env", "p1gen",'p1rep'),cols_remove = FALSE)
H3_amp <- separate_wider_delim(H3_amp, cols = p2, delim = "_", names = c("p2env", "p2gen",'p2rep'),cols_remove = FALSE)
H3_amp$env <- ifelse(H3_amp$p1env == H3_amp$p2env, 'within', 'between')
H3_amp$p1_gen_num <- as.numeric(sub('G','',H3_amp$p1gen))
H3_amp$p2_gen_num <- as.numeric(sub('G','',H3_amp$p2gen))
H3_amp$dist <- ifelse(H3_amp$p1env == H3_amp$p2env & H3_amp$p1rep == H3_amp$p2rep, abs(H3_amp$p1_gen_num - H3_amp$p2_gen_num), H3_amp$p1_gen_num + H3_amp$p2_gen_num)
H3_amp$env[H3_amp$dist<=100]='within'

H3_amp <- subset(H3_amp,p1gen==p2gen)
H3_amp$dist <- H3_amp$p1_gen_num
#H3_amp$dist <- as.factor(H3_amp$dist)


H3_amp$env_final <- ifelse(H3_amp$env=='within' & H3_amp$p1env == 'NaCl','NaCl',ifelse(H3_amp$env=='within' & H3_amp$p1env == 'Ethanol','Ethanol','Divergent'))


H3_overall<- ggplot(H3_amp, aes(x=env_final, metric)) +
  geom_hline(yintercept=0.21147270, linetype="dashed", 
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
    axis.title.y=element_blank(),
    axis.title.x=element_blank(),
    legend.position = "bottom",
    axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  guides(fill=guide_legend(title="Environment"))+
  ylim(0,0.82)

H3_scatter<- ggplot(H3_amp, aes(x=dist, metric)) +
  geom_hline(yintercept=0.21147270, linetype="dashed", 
             color = "black", size=0.5) + 
  stat_summary(aes(y=metric, group=env_final,color=env_final),fun=mean, geom="line",alpha=1,,show.legend=FALSE)+
  stat_summary(aes(x=dist,y=metric,group=env_final),fun.data = mean_se,  geom = "errorbar",width=0,show.legend=FALSE) + 
  stat_summary(aes(x=dist,y=metric,fill=env_final),fun=mean, geom="point", shape=21, size=2,stroke=0.7,alpha=1,show.legend=FALSE)+
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
  ylim(0,0.82)

H1_overall + H1_scatter +H2_overall + H2_scatter +H3_overall + H3_scatter +
  plot_layout(widths=c(2,5,2,5,2,5)) 

ggsave('figures/Hybrid_Dynamics_Jaccard-gen.pdf',width=7.5,height=2,dpi = 900)


# ==============================================================================
# Pairwise Kruskal-Wallis
# ==============================================================================
kruskal.test(metric ~ env_final, data = H1_amp)
pairwise.wilcox.test(H1_amp$metric, H1_amp$env_final,
                     p.adjust.method = "BH")
kruskal.test(metric ~ env_final, data = H2_amp)
pairwise.wilcox.test(H2_amp$metric, H2_amp$env_final,
                     p.adjust.method = "BH")
kruskal.test(metric ~ env_final, data = H3_amp)
pairwise.wilcox.test(H3_amp$metric, H3_amp$env_final,
                     p.adjust.method = "BH")
# ==============================================================================
# 
# ==============================================================================
library(multcomp)
H1_amp$env_final <- as.factor(H1_amp$env_final)
res_aov <- aov(metric ~ env_final,
               data = H1_amp)
summary(res_aov)
post_test <- glht(res_aov,
                  linfct = mcp(env_final = "Tukey"))
summary(post_test)

H2_amp$env_final <- as.factor(H2_amp$env_final)
res_aov <- aov(metric ~ env_final,
               data = H2_amp)
summary(res_aov)
post_test <- glht(res_aov,
                  linfct = mcp(env_final = "Tukey"))
summary(post_test)

H3_amp$env_final <- as.factor(H3_amp$env_final)
res_aov <- aov(metric ~ env_final,
               data = H3_amp)
summary(res_aov)
post_test <- glht(res_aov,
                  linfct = mcp(env_final = "Tukey"))
summary(post_test)

