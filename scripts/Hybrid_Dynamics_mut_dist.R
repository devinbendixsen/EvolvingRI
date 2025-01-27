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
library(tibble)
library(stringr)


genome_size=12071326
chr_data <- data.frame(chr=c('I','II','III','IV','V','VI','VII','VIII','IX','X','XI','XII','XIII','XIV','XV','XVI'),size=c(230218,	813184,	316620,	1531933,	576874,	270161,	1090940,	562643,	439888,	745751,	666816,	1078177,	924431,	784333,	1091291,	948066))
old_chr <- data.frame(chr=c("BK006935.2", "BK006936.2","BK006937.2",  "BK006938.2", "BK006939.2", "BK006940.2","BK006941.2",  "BK006934.2", "BK006942.2",  "BK006943.2", "BK006944.2",  "BK006945.2", "BK006946.2", "BK006947.3",  "BK006948.2", "BK006949.2"))

# ==============================================================================
# PLOT CHROMOSOME SIZES
# ==============================================================================
ggplot(chr_data) + 
  geom_bar(aes(x=chr))
ggplot(chr_data, aes(x = chr, y = size)) +
  geom_bar(stat = "identity")
# ==============================================================================
# LOAD AND FORMAT DATA
# ==============================================================================
fileNames <- Sys.glob("results/FREEC/*_CNVs.p.value.txt") # identify all samples
cnv_data <- data.frame()
for (fileName in fileNames) {
  print(fileName)
  cnv_sample <- read.table(file=fileName,header=FALSE,skip = 1,col.names=c('chr','start','end','CN','type','Wilcoxon','KS'))
  cnv_sample$length <- (cnv_sample$end - cnv_sample$start)
  cnv_sample <- subset(cnv_sample,Wilcoxon<0.05 & KS <0.05)
  cnv_sample <- cnv_sample %>%
    mutate(amp=if_else(CN>2,length,0)) %>%
    mutate(del=if_else(CN<2,length,0))
  
  for (chr in old_chr$chr) {
    cnv_sample <- rbind(cnv_sample,data.frame('chr'=chr,'start'=0,'end'=0,'CN'=2,'type'='x','Wilcoxon'=0,'KS'=0,'length'=0,'amp'=0,'del'=0))
  }
  
  cnv_summary <- data.frame(row.names=c(old_chr$chr))
  cnv <- split(cnv_sample,cnv_sample$chr)
  
  for (chro_num in unique(cnv_sample$chr)) {
    cnv_summary[chro_num,'cnv']=sum(cnv[[chro_num]]$length)
    cnv_summary[chro_num,'amp']=sum(cnv[[chro_num]]$amp)
    cnv_summary[chro_num,'del']=sum(cnv[[chro_num]]$del)
  }
  
  sample=sub('results/FREEC/','',fileName)
  sample=sub('.mapped.sort.picard.bam_CNVs.p.value.txt','',sample)
  env=strsplit(sample,split="_")[[1]][1]
  gen=strsplit(sample,split="_")[[1]][2]
  rep=strsplit(sample,split="_")[[1]][3]
  
  cnv_summary$sample = sample
  cnv_summary$env = env
  cnv_summary$gen = gen
  cnv_summary$rep = rep
  cnv_summary <- tibble::rownames_to_column(cnv_summary, "chr")

  cnv_data <- rbind(cnv_data,cnv_summary)
  
}

cnv_data <- cnv_data %>% 
  mutate(chr=recode(chr, "BK006935.2"="I", 
                    "BK006936.2"="II",
                    "BK006937.2"= "III", 
                    "BK006938.2"= "IV", 
                    "BK006939.2"= "V",
                    "BK006940.2"= "VI",
                    "BK006941.2"="VII",
                    "BK006934.2"="VIII",
                    "BK006942.2"="IX", 
                    "BK006943.2"= "X",
                    "BK006944.2"="XI",
                    "BK006945.2"="XII",
                    "BK006946.2"="XIII",
                    "BK006947.3"="XIV",
                    "BK006948.2"="XV",
                    "BK006949.2"="XVI"))


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
  
  sv_table <- data.frame(table(mut_data_sample$CHROM))
  
  sv_table$sample = sample
  sv_table$env = env
  sv_data <- rbind(sv_data,sv_table)
  
}

sv_data <- sv_data %>%
  dplyr::rename(chr=Var1,sv=Freq)
table(sv_data$chr)
sv_data <- sv_data %>% 
  mutate(chr=recode(chr, "BK006935.2"="I", 
                    "BK006936.2"="II",
                    "BK006937.2"= "III", 
                    "BK006938.2"= "IV", 
                    "BK006939.2"= "V",
                    "BK006940.2"= "VI",
                    "BK006941.2"="VII",
                    "BK006934.2"="VIII",
                    "BK006942.2"="IX", 
                    "BK006943.2"= "X",
                    "BK006944.2"="XI",
                    "BK006945.2"="XII",
                    "BK006946.2"="XIII",
                    "BK006947.3"="XIV",
                    "BK006948.2"="XV",
                    "BK006949.2"="XVI"))
table(sv_data$chr)

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
  
  snv_sample <- data.frame(row.names=c('snv'))

  
  snv_sample$I = as.numeric(snv_sample_stat$V3[18])
  snv_sample$II = as.numeric(snv_sample_stat$V3[20])
  snv_sample$III = as.numeric(snv_sample_stat$V3[22])
  snv_sample$IV = as.numeric(snv_sample_stat$V3[24])
  snv_sample$V = as.numeric(snv_sample_stat$V3[26])
  snv_sample$VI = as.numeric(snv_sample_stat$V3[28])
  snv_sample$VII = as.numeric(snv_sample_stat$V3[30])
  snv_sample$VIII = as.numeric(snv_sample_stat$V3[32])
  snv_sample$IX = as.numeric(snv_sample_stat$V3[34])
  snv_sample$X = as.numeric(snv_sample_stat$V3[36])
  snv_sample$XI = as.numeric(snv_sample_stat$V3[38])
  snv_sample$XII = as.numeric(snv_sample_stat$V3[40])
  snv_sample$XIII = as.numeric(snv_sample_stat$V3[42])
  snv_sample$XIV = as.numeric(snv_sample_stat$V3[44])
  snv_sample$XV = as.numeric(snv_sample_stat$V3[46])
  snv_sample$XVI = as.numeric(snv_sample_stat$V3[48])
  
  snv_sample <- data.frame(t(snv_sample))
  
  
  snv_sample <- tibble::rownames_to_column(snv_sample, "chr")
  snv_sample$sample = sample
  snv_sample$env = env
  snv_sample$gen = gen
  snv_sample$rep = rep
  snv_data <- rbind(snv_data,snv_sample)
}


# ==============================================================================
# Plot SNV distribution
# ==============================================================================
snv_mean <- snv_data %>%
  dplyr::group_by(env,chr) %>% 
  dplyr::summarise(mean_snv=mean(snv))

snv_mean <- left_join(snv_mean,chr_data)
snv_mean$mean_snv <- (snv_mean$mean_snv/snv_mean$size)*1000

snv_mean <- snv_mean %>%
  mutate(env=factor(env, levels=rev(c('N','LE','NaCl','LiAc0.01','LiAc0.02','Ethanol','H1F1','H1F2','H2F1','H2F2','H3F1','H3F2'))),
         chr=factor(chr,levels=c('I','II','III','IV','V','VI','VII','VIII','IX','X','XI','XII','XIII','XIV','XV','XVI')))

p_SNV<- ggplot(snv_mean, aes(x=chr, y=env)) +
  geom_tile( aes(fill=mean_snv),color = "white",
             lwd = 1,
             linetype = 1)+
  scale_fill_gradient(low='white',high='#136B6B',,
                      guide = guide_colorbar(frame.colour = "black", ticks.colour = "black",barwidth=5,barheight=1,title = 'SNV (muts/kb)',title.position='top'))+
  theme(panel.background = element_rect(fill = "white"),
        plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
        legend.position = "bottom",
        axis.title.y=element_blank()) +
  xlab('Chromosomes') + 
  scale_y_discrete(breaks=rev(c('N','LE','NaCl','LiAc0.01','LiAc0.02','Ethanol','H1F1','H1F2','H2F1','H2F2','H3F1','H3F2')),
                   labels=rev(c('N Founder','LE Founder','NaCl','LiAc 0.01','LiAc 0.02','Ethanol','NaCl x LiAc 0.01 F1','NaCl x LiAc 0.01 F2','NaCl x LiAc 0.02 F1','NaCl x LiAc 0.02 F2','NaCl x EtOH F1','NaCl x EtOH F2')))

p_SNV_right <- ggplot(snv_mean, aes(y = env, x = mean_snv,fill=env)) +
  geom_boxplot(color='black') +
  scale_fill_manual(values=rev(c("white", "white",'#4472c4ff','#9e49e1ff','#00b050ff','#ff0000ff','#808080','#d4a373','#808080','#d4a373','#808080','#d4a373')))+
  theme(panel.background = element_rect(fill = 'lightgrey'),
        plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
        legend.position = "none",
        axis.text.y=element_blank(),
        axis.title.y=element_blank()) +
  xlab('SNV load')+
  scale_x_continuous(minor_breaks = NULL) 

p_SNV_top <- ggplot(snv_mean, aes(y = mean_snv, x = chr)) +
  geom_boxplot(color='black',fill='white') +
  theme(panel.background = element_rect(fill = 'lightgrey'),
        plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
        legend.position = "none",
        axis.text.x=element_blank(),
        axis.title.x=element_blank()) +
  ylab('SNV load')+
  scale_y_continuous(minor_breaks = NULL) 


p_SNV_total <- p_SNV_top + guide_area() + p_SNV + p_SNV_right + plot_layout(heights=c(1,3),widths=c(6,2),guides='collect')
#ggsave('figures/Hybrid_Dynamics_SNV_dist.pdf',width=10,height=3,dpi = 900)

# ==============================================================================
# Plot SV distribution
# ==============================================================================
sv_mean <- sv_data %>%
  dplyr::group_by(env,chr) %>% 
  dplyr::summarise(mean_sv=mean(sv))

sv_mean <- left_join(sv_mean,chr_data)
sv_mean$mean_sv <- (sv_mean$mean_sv/sv_mean$size)*1000
sv_mean <- sv_mean %>%
  mutate(env=factor(env, levels=rev(c('N','LE','NaCl','LiAc0.01','LiAc0.02','Ethanol','H1F1','H1F2','H2F1','H2F2','H3F1','H3F2'))),
         chr=factor(chr,levels=c('I','II','III','IV','V','VI','VII','VIII','IX','X','XI','XII','XIII','XIV','XV','XVI')))

p_SV<- ggplot(sv_mean, aes(x=chr, y=env)) +
  geom_tile( aes(fill=mean_sv),color = "white",
             lwd = 1,
             linetype = 1)+
  scale_fill_gradient(low='white',high='#871F1F',,
                      guide = guide_colorbar(frame.colour = "black", ticks.colour = "black",barwidth=5,barheight=1,title = 'SV (muts/kb)',title.position='top'))+
  theme(panel.background = element_rect(fill = "white"),
        plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
        legend.position = "bottom",
        axis.title.y=element_blank()) +
  xlab('Chromosomes') + 
  scale_y_discrete(breaks=rev(c('N','LE','NaCl','LiAc0.01','LiAc0.02','Ethanol','H1F1','H1F2','H2F1','H2F2','H3F1','H3F2')),
                   labels=rev(c('N Founder','LE Founder','NaCl','LiAc 0.01','LiAc 0.02','Ethanol','NaCl x LiAc 0.01 F1','NaCl x LiAc 0.01 F2','NaCl x LiAc 0.02 F1','NaCl x LiAc 0.02 F2','NaCl x EtOH F1','NaCl x EtOH F2')))

p_SV_right <- ggplot(sv_mean, aes(y = env, x = mean_sv,fill=env)) +
  geom_boxplot(color='black') +
  scale_fill_manual(values=rev(c("white", "white",'#4472c4ff','#9e49e1ff','#00b050ff','#ff0000ff','#808080','#d4a373','#808080','#d4a373','#808080','#d4a373')))+
  theme(panel.background = element_rect(fill = 'lightgrey'),
        plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
        legend.position = "none",
        axis.text.y=element_blank(),
        axis.title.y=element_blank()) +
  xlab('SV load')+
  scale_x_continuous(minor_breaks = NULL) 

p_SV_top <- ggplot(sv_mean, aes(y = mean_sv, x = chr)) +
  geom_boxplot(color='black',fill='white') +
  theme(panel.background = element_rect(fill = 'lightgrey'),
        plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
        legend.position = "none",
        axis.text.x=element_blank(),
        axis.title.x=element_blank()) +
  ylab('SV load')+
  scale_y_continuous(minor_breaks = NULL) 

p_SV_total <- p_SV_top +guide_area() + p_SV + p_SV_right + plot_layout(heights=c(1,3),widths=c(6,2),guides='collect')
#ggsave('figures/Hybrid_Dynamics_SV_dist.pdf',width=10,height=3,dpi = 900)


# ==============================================================================
# Plot aneuploidy
# ==============================================================================
cnv_aneuploidy <- cnv_data
cnv_aneuploidy <- left_join(cnv_aneuploidy,chr_data)
cnv_aneuploidy$cnv <- (cnv_aneuploidy$cnv/cnv_aneuploidy$size)*100
cnv_aneuploidy$amp <- (cnv_aneuploidy$amp/cnv_aneuploidy$size)*100
cnv_aneuploidy$del <- (cnv_aneuploidy$del/cnv_aneuploidy$size)*100
cnv_aneuploidy <- cnv_aneuploidy %>%
  mutate(aneu=if_else(amp>=del,amp,del))
  
cnv_aneuploidy <- cnv_aneuploidy %>%
  mutate(aneu_count=if_else(aneu>60,1,0))
cnv_aneu_count <- cnv_aneuploidy %>%
  group_by(env,sample) %>%
  summarise(n=sum(aneu_count))

founder <- subset(cnv_aneu_count, sample %in% c("N_Founder","LE_Founder"))

H1 <- subset(cnv_aneu_count, (env == "NaCl" | env == "LiAc0.01" | env =='H1F1' | env =='H1F2'))
H1$env <- factor(H1$env , levels=c('NaCl','LiAc0.01','H1F1','H1F2'))

H1_aneu <- ggplot(H1) + 
  geom_hline(yintercept=(mean(founder$n)),color='black',linetype='dashed',linewidth=0.5)+
  geom_boxplot(aes(x=as.factor(env),y=n,fill=env)) +
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
        legend.position = "none") + 
  ylab('Aneuploidies') 
H1_aneu
H2 <- subset(cnv_aneu_count, (env == "NaCl" | env == "LiAc0.02" | env =='H2F1' | env =='H2F2'))
H2$env <- factor(H2$env , levels=c('NaCl','LiAc0.02','H2F1','H2F2'))

H2_aneu <- ggplot(H2) + 
  geom_hline(yintercept=(mean(founder$n)),color='black',linetype='dashed',linewidth=0.5)+
  geom_boxplot(aes(x=as.factor(env),y=n,fill=env)) +
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
        axis.title.y=element_blank(),
        panel.spacing = unit(0.5, "lines"), 
        legend.position = "none")
H2_aneu
H3 <- subset(cnv_aneu_count, (env == "NaCl" | env == "Ethanol" | env =='H3F1' | env =='H3F2'))
H3$env <- factor(H3$env , levels=c('NaCl','Ethanol','H3F1','H3F2'))

H3_aneu <- ggplot(H3) + 
  geom_hline(yintercept=(mean(founder$n)),color='black',linetype='dashed',linewidth=0.5)+
  geom_boxplot(aes(x=as.factor(env),y=n,fill=env)) +
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
        axis.title.y=element_blank(),
        panel.spacing = unit(0.5, "lines"), 
        legend.position = "none")
H1_aneu+H2_aneu+H3_aneu

kruskal.test(n ~ env, data = H1)
pairwise.wilcox.test(H1$n, H1$env,
                     p.adjust.method = "BH")
kruskal.test(n ~ env, data = H2)
pairwise.wilcox.test(H2$n, H2$env,
                     p.adjust.method = "BH")
kruskal.test(n ~ env, data = H3)
pairwise.wilcox.test(H3$n, H3$env,
                     p.adjust.method = "BH")

H1 %>%
  group_by(env) %>%
  summarise(mean=mean(n),n = n())
H2 %>%
  group_by(env) %>%
  summarise(mean=mean(n),n = n())
H3 %>%
  group_by(env) %>%
  summarise(mean=mean(n),n = n())
# ==============================================================================
# Plot CNV distribution
# ==============================================================================

cnv_mean <- cnv_data %>%
  dplyr::group_by(env,chr) %>% 
  dplyr::summarise(mean_cnv=mean(cnv))

cnv_mean <- left_join(cnv_mean,chr_data)
cnv_mean$mean_cnv <- (cnv_mean$mean_cnv/cnv_mean$size)*100
cnv_mean <- cnv_mean %>%
  mutate(env=factor(env, levels=rev(c('N','LE','NaCl','LiAc0.01','LiAc0.02','Ethanol','H1F1','H1F2','H2F1','H2F2','H3F1','H3F2'))),
         chr=factor(chr,levels=c('I','II','III','IV','V','VI','VII','VIII','IX','X','XI','XII','XIII','XIV','XV','XVI')))

p_CNV<- ggplot(cnv_mean, aes(x=chr, y=env)) +
  geom_tile( aes(fill=mean_cnv),color = "white",
             lwd = 1,
             linetype = 1)+
  scale_fill_gradient(low='white',high='#23395d',
                      guide = guide_colorbar(frame.colour = "black", ticks.colour = "black",barwidth=5,barheight=1,title = 'CNV (%)',title.position='top'))+
  theme(panel.background = element_rect(fill = "white"),
        plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
        legend.position = "bottom",
        axis.title.y=element_blank()) +
  xlab('Chromosomes') + 
  scale_y_discrete(breaks=rev(c('N','LE','NaCl','LiAc0.01','LiAc0.02','Ethanol','H1F1','H1F2','H2F1','H2F2','H3F1','H3F2')),
                   labels=rev(c('N Founder','LE Founder','NaCl','LiAc 0.01','LiAc 0.02','Ethanol','NaCl x LiAc 0.01 F1','NaCl x LiAc 0.01 F2','NaCl x LiAc 0.02 F1','NaCl x LiAc 0.02 F2','NaCl x EtOH F1','NaCl x EtOH F2')))


p_CNV_right <- ggplot(cnv_mean, aes(y = env, x = mean_cnv,fill = env)) +
  geom_boxplot(color='black') +
  scale_fill_manual(values=rev(c("white", "white",'#4472c4ff','#9e49e1ff','#00b050ff','#ff0000ff','#808080','#d4a373','#808080','#d4a373','#808080','#d4a373')))+
  theme(panel.background = element_rect(fill = 'lightgrey'),
        plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
        legend.position = "none",
        axis.text.y=element_blank(),
        axis.title.y=element_blank()) +
  xlab('CNV (%)') +
  scale_x_continuous(,minor_breaks = NULL) 

p_CNV_top <- ggplot(cnv_mean, aes(y = mean_cnv, x = chr)) +
  geom_boxplot(color='black',fill='white') +
  theme(panel.background = element_rect(fill = 'lightgrey'),
        plot.margin = margin(0.1, 0.1, 0.1, 0.1, "cm"),
        legend.position = "top",
        axis.text.x=element_blank(),
        axis.title.x=element_blank()) +
  ylab('CNV (%)') +
  scale_y_continuous(minor_breaks = NULL)

p_CNV_total <- p_CNV_top + guide_area() + p_CNV + p_CNV_right + plot_layout(heights=c(1,3),widths=c(6,2),guides='collect')
#ggsave('figures/Hybrid_Dynamics_CNV_dist.pdf',width=10,height=3,dpi = 900)

p_SNV_total/p_SV_total/p_CNV_total
ggsave('figures/Hybrid_Dynamics_mut_dist_SNV_SV_CNV.pdf',width=8,height=8,dpi = 900)


# ==============================================================================
# STATISTICS
# ==============================================================================
#SNV stats
snv_mean %>%
  dplyr::group_by(chr) %>% 
  dplyr::summarise(mean = mean(mean_snv))

kruskal.test(mean_snv ~ env, data = snv_mean)
pairwise.wilcox.test(snv_mean$mean_snv, snv_mean$env,
                     p.adjust.method = "BH")

kruskal.test(mean_snv ~ chr, data = snv_mean)
pairwise.wilcox.test(snv_mean$mean_snv, snv_mean$chr,
                     p.adjust.method = "BH")


#SV stats
sv_mean %>%
  dplyr::group_by(chr) %>% 
  dplyr::summarise(mean = mean(mean_sv))

kruskal.test(mean_sv ~ env, data = sv_mean)
pairwise.wilcox.test(sv_mean$mean_sv, sv_mean$env,
                     p.adjust.method = "BH")

kruskal.test(mean_sv ~ chr, data = sv_mean)
pairwise.wilcox.test(sv_mean$mean_sv, sv_mean$chr,
                     p.adjust.method = "BH")



#CNV stats
cnv_mean %>% 
  dplyr::group_by(env) %>% 
  dplyr::summarise(mean = mean(mean_cnv))
cnv_mean %>%
  dplyr::group_by(chr) %>% 
  dplyr::summarise(mean = mean(mean_cnv))

kruskal.test(mean_cnv ~ env, data = cnv_mean)
pairwise.wilcox.test(cnv_mean$mean_cnv, cnv_mean$env,
                     p.adjust.method = "BH")

kruskal.test(mean_cnv ~ chr, data = cnv_mean)
pairwise.wilcox.test(cnv_mean$mean_cnv, cnv_mean$chr,
                     p.adjust.method = "BH")



# ==============================================================================
# 
# ==============================================================================
chr_size <- ggplot(chr_data, aes(x = chr, y = size/1000)) +
  geom_bar(stat = "identity")+
  ylab('Chromosome Size (kb)')+
  xlab('Chromosome')+
  ylim(0,1600)
chr_size
cnv_mean <- cnv_mean %>%
  subset(!env %in% c('N','LE')) %>%
  mutate(env=factor(env, levels=rev(c('NaCl','LiAc0.01','LiAc0.02','Ethanol','H1F1','H1F2','H2F1','H2F2','H3F1','H3F2'))),
         chr=factor(chr,levels=c('I','II','III','IV','V','VI','VII','VIII','IX','X','XI','XII','XIII','XIV','XV','XVI')))
cnv_size <- ggplot(cnv_mean, aes(x=mean_cnv, y=size/1000)) + 
  geom_smooth(method=lm, se=TRUE, fullrange=TRUE,color='black')+
  stat_cor(label.x.npc = 0,label.y.npc = 0,p.accuracy=0.01,r.accuracy=0.01,size=2.5)+
  geom_point(alpha=alpha,aes(fill=env),shape=21, colour='black')+
  scale_fill_manual(values=rev(c('#4472c4ff','#9e49e1ff','#00b050ff','#ff0000ff','#808080','#d4a373','#808080','#d4a373','#808080','#d4a373')))+
  theme(legend.position='none',
        axis.title.y=element_blank(), 
        axis.text.y=element_blank())+
  xlab('CNV (%)') +
  ylab('Chromosome Size (kb)')+
  ylim(0,1600)

snv_mean <- snv_mean %>%
  subset(!env %in% c('N','LE')) %>%
  mutate(env=factor(env, levels=rev(c('NaCl','LiAc0.01','LiAc0.02','Ethanol','H1F1','H1F2','H2F1','H2F2','H3F1','H3F2'))),
         chr=factor(chr,levels=c('I','II','III','IV','V','VI','VII','VIII','IX','X','XI','XII','XIII','XIV','XV','XVI')))
snv_size <- ggplot(snv_mean, aes(x=mean_snv, y=size/1000)) + 
  geom_smooth(method=lm, se=TRUE, fullrange=TRUE,color='black')+
  stat_cor(label.x.npc = 0,label.y.npc = 0,p.accuracy=0.01,r.accuracy=0.01,size=2.5)+
  geom_point(alpha=alpha,aes(fill=env),shape=21, colour='black')+
  scale_fill_manual(values=rev(c('#4472c4ff','#9e49e1ff','#00b050ff','#ff0000ff','#808080','#d4a373','#808080','#d4a373','#808080','#d4a373'))) +
  theme(legend.position='none',
        axis.title.y=element_blank(), 
        axis.text.y=element_blank())+
  xlab('SNV (muts/kb)') +
  ylab('Chromosome Size (kb)')+
  ylim(0,1600)

sv_mean <- sv_mean %>%
  subset(!env %in% c('N','LE')) %>%
  mutate(env=factor(env, levels=rev(c('NaCl','LiAc0.01','LiAc0.02','Ethanol','H1F1','H1F2','H2F1','H2F2','H3F1','H3F2'))),
         chr=factor(chr,levels=c('I','II','III','IV','V','VI','VII','VIII','IX','X','XI','XII','XIII','XIV','XV','XVI')))
sv_size <- ggplot(sv_mean, aes(x=mean_sv, y=size/1000)) + 
  geom_smooth(method=lm, se=TRUE, fullrange=TRUE,color='black')+
  stat_cor(label.x.npc = 0,label.y.npc = 0,p.accuracy=0.01,r.accuracy=0.01,size=2.5)+
  geom_point(alpha=alpha,aes(fill=env),shape=21, colour='black')+
  scale_fill_manual(values=rev(c('#4472c4ff','#9e49e1ff','#00b050ff','#ff0000ff','#808080','#d4a373','#808080','#d4a373','#808080','#d4a373')))+
  xlab('SV (muts/kb)') +
  ylab('Chromosome Size (kb)')+
  ylim(0,1600) +
  theme(legend.position='none',
        axis.title.y=element_blank(), 
        axis.text.y=element_blank())

chr_size + snv_size + sv_size + cnv_size +plot_layout(widths=c(2.5,1,1,1))
ggsave('figures/Hybrid_Dynamics_SNV_SV_CNV_size.pdf',width=8,height=2,dpi = 900)


# ==============================================================================
# 
# ==============================================================================
chr_size <- ggplot(chr_data, aes(y = chr, x = size/1000)) +
  geom_bar(stat = "identity")+
  xlab('Chromosome Size (kb)')+
  ylab('Chromosome')+
  xlim(0,1600)+
  scale_x_continuous(position = "top") 
chr_size
cnv_mean <- cnv_mean %>%
  subset(!env %in% c('N','LE')) %>%
  mutate(env=factor(env, levels=rev(c('NaCl','LiAc0.01','LiAc0.02','Ethanol','H1F1','H1F2','H2F1','H2F2','H3F1','H3F2'))),
         chr=factor(chr,levels=c('I','II','III','IV','V','VI','VII','VIII','IX','X','XI','XII','XIII','XIV','XV','XVI')))
cnv_size <- ggplot(cnv_mean, aes(y=mean_cnv, x=size/1000)) + 
  geom_smooth(method=lm, se=TRUE, fullrange=TRUE,color='black')+
  stat_cor(label.x.npc = 0.5,label.y.npc = .75,p.accuracy=0.01,r.accuracy=0.01,size=2.5)+
  geom_point(alpha=alpha,aes(colour=env),shape=1)+
  scale_colour_manual(values=rev(c('#4472c4ff','#9e49e1ff','#00b050ff','#ff0000ff','#808080','#d4a373','#808080','#d4a373','#808080','#d4a373')))+
  theme(legend.position='none',
        axis.title.x=element_blank(), 
        axis.text.x=element_blank())+
  ylab('CNV (%)') +
  xlab('Chromosome Size (kb)')+
  xlim(0,1600)
cnv_size
snv_mean <- snv_mean %>%
  subset(!env %in% c('N','LE')) %>%
  mutate(env=factor(env, levels=rev(c('NaCl','LiAc0.01','LiAc0.02','Ethanol','H1F1','H1F2','H2F1','H2F2','H3F1','H3F2'))),
         chr=factor(chr,levels=c('I','II','III','IV','V','VI','VII','VIII','IX','X','XI','XII','XIII','XIV','XV','XVI')))
snv_size <- ggplot(snv_mean, aes(y=mean_snv, x=size/1000)) + 
  geom_smooth(method=lm, se=TRUE, fullrange=TRUE,color='black')+
  stat_cor(label.x.npc = 0.5,label.y.npc = 0.75,p.accuracy=0.01,r.accuracy=0.01,size=2.5)+
  geom_point(alpha=alpha,aes(colour=env),shape=1)+
  scale_colour_manual(values=rev(c('#4472c4ff','#9e49e1ff','#00b050ff','#ff0000ff','#808080','#d4a373','#808080','#d4a373','#808080','#d4a373'))) +
  theme(legend.position='none',
        axis.title.x=element_blank(), 
        axis.text.x=element_blank())+
  ylab('SNV (muts/kb)') +
  xlab('Chromosome Size (kb)')+
  xlim(0,1600)
snv_size
sv_mean <- sv_mean %>%
  subset(!env %in% c('N','LE')) %>%
  mutate(env=factor(env, levels=rev(c('NaCl','LiAc0.01','LiAc0.02','Ethanol','H1F1','H1F2','H2F1','H2F2','H3F1','H3F2'))),
         chr=factor(chr,levels=c('I','II','III','IV','V','VI','VII','VIII','IX','X','XI','XII','XIII','XIV','XV','XVI')))
sv_size <- ggplot(sv_mean, aes(y=mean_sv, x=size/1000)) + 
  geom_smooth(method=lm, se=TRUE, fullrange=TRUE,color='black')+
  stat_cor(label.x.npc = 0.5,label.y.npc = 0.75,p.accuracy=0.01,r.accuracy=0.01,size=2.5)+
  geom_point(alpha=alpha,aes(colour=env),shape=1)+
  scale_colour_manual(values=rev(c('#4472c4ff','#9e49e1ff','#00b050ff','#ff0000ff','#808080','#d4a373','#808080','#d4a373','#808080','#d4a373')))+
  ylab('SV (muts/kb)') +
  xlab('Chromosome Size (kb)')+
  xlim(0,1600) +
  theme(legend.position='none',
        axis.title.x=element_blank(), 
        axis.text.x=element_blank())

chr_size / snv_size / sv_size / cnv_size +plot_layout(heights=c(2.5,1,1,1))
ggsave('figures/Hybrid_Dynamics_SNV_SV_CNV_size.pdf',width=2,height=8,dpi = 900)



