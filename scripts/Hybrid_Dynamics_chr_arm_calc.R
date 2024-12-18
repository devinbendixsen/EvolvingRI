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
library(tidyverse)
library("readxl")

# ==============================================================================
# 
# ==============================================================================
fileNames <- Sys.glob("results/FREEC/*_CNVs") # identify all samples

CEN_data <- read_excel('data/reference_genome/genome_chr_data.xlsx',sheet="CEN")
rownames(CEN_data) <- CEN_data$chr

cnv_data <- data.frame()
for (fileName in fileNames) {
  print(fileName)
  chr_data <- read_excel('data/reference_genome/genome_chr_data.xlsx',sheet="Sheet3")
  cnv_sample <- read.table(file=fileName,header=FALSE,col.names=c('chr','start','end','CN','type'))
  cnv_sample$length <- (cnv_sample$end - cnv_sample$start)
  cnv_sample <- cnv_sample %>% 
    mutate(chr = str_replace(chr, "BK006935.2", "chr1"), 
           chr = str_replace(chr, "BK006936.2", "chr2"),
           chr = str_replace(chr, "BK006937.2", "chr3"),
           chr = str_replace(chr, "BK006938.2", "chr4"),
           chr = str_replace(chr, "BK006939.2", "chr5"),
           chr = str_replace(chr, "BK006940.2", "chr6"),
           chr = str_replace(chr, "BK006941.2", "chr7"),
           chr = str_replace(chr, "BK006934.2", "chr8"),
           chr = str_replace(chr, "BK006942.2", "chr9"),
           chr = str_replace(chr, "BK006943.2", "chr10"),
           chr = str_replace(chr, "BK006944.2", "chr11"),
           chr = str_replace(chr, "BK006945.2", "chr12"),
           chr = str_replace(chr, "BK006946.2", "chr13"),
           chr = str_replace(chr, "BK006947.3", "chr14"),
           chr = str_replace(chr, "BK006948.2", "chr15"),
           chr = str_replace(chr, "BK006949.2", "chr16"))
  chr_data$loss = 0
  chr_data$gain = 0
  for (variant in 1:nrow(cnv_sample)) {
    for (chr_arm in 1:nrow(chr_data)) {
      if (cnv_sample[variant,'chr'] == chr_data[chr_arm,'chr']) {
        if (cnv_sample[variant,'start'] >= chr_data[chr_arm,'start'] & cnv_sample[variant,'end'] <= chr_data[chr_arm,'end']) {
          cnv_sample[variant,'chr_arm'] <- chr_data[chr_arm,'chr_arm']
          type <- cnv_sample[variant,'type']
          chr_data[chr_arm,type] = chr_data[chr_arm,type]+cnv_sample[variant,'length']
        }
        else if (cnv_sample[variant,'start'] <= chr_data[chr_arm,'end'] & cnv_sample[variant,'end'] >= chr_data[chr_arm,'end']) {
          var1_start=cnv_sample[variant,'start']
          var1_end=(CEN_data[cnv_sample[variant,'chr'],'CEN'])
          type <- cnv_sample[variant,'type']
          chr_data[chr_arm,type] = chr_data[chr_arm,type]+abs(var1_start-var1_end)
        }
        if (cnv_sample[variant,'start'] < chr_data[chr_arm,'start'] & cnv_sample[variant,'end'] >= chr_data[chr_arm,'start']) {
          var2_start=CEN_data[cnv_sample[variant,'chr'],'CEN']
          var2_end=cnv_sample[variant,'end']
          type <- cnv_sample[variant,'type']
          chr_data[chr_arm,type] = chr_data[chr_arm,type]+abs(var2_start-var2_end)
        }
      }
    }
  }
  
  
  cnv <- split(cnv_sample,cnv_sample$type)
  amp = sum(cnv$gain$length)
  del = sum(cnv$loss$length)
  
  
  
  sample=sub('results/FREEC/','',fileName)
  sample=sub('.mapped.sort.bam_CNVs','',sample)
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



# ==============================================================================
# PRACTICE
# ==============================================================================
fileName = 'results/FREEC/Ethanol_G1000_R3.mapped.sort.picard.bam_CNVs'
cnv_sample <- read.table(file=fileName,header=FALSE,col.names=c('chr','start','end','CN','type'))
chr_data$loss = 0
chr_data$gain = 0
for (variant in 1:nrow(cnv_sample)) {
  for (chr_arm in 1:nrow(chr_data)) {
    if (cnv_sample[variant,'chr'] == chr_data[chr_arm,'chr']) {
      if (cnv_sample[variant,'start'] >= chr_data[chr_arm,'start'] & cnv_sample[variant,'end'] <= chr_data[chr_arm,'end']) {
        cnv_sample[variant,'chr_arm'] <- chr_data[chr_arm,'chr_arm']
        type <- cnv_sample[variant,'type']
        chr_data[chr_arm,type] = chr_data[chr_arm,type]+cnv_sample[variant,'length']
      }
      else if (cnv_sample[variant,'start'] <= chr_data[chr_arm,'end'] & cnv_sample[variant,'end'] >= chr_data[chr_arm,'end']) {
        var1_start=cnv_sample[variant,'start']
        var1_end=(CEN_data[cnv_sample[variant,'chr'],'CEN'])
        type <- cnv_sample[variant,'type']
        chr_data[chr_arm,type] = chr_data[chr_arm,type]+abs(var1_start-var1_end)
        print(abs(var1_start-var1_end))
      }
      if (cnv_sample[variant,'start'] < chr_data[chr_arm,'start'] & cnv_sample[variant,'end'] >= chr_data[chr_arm,'start']) {
        var2_start=CEN_data[cnv_sample[variant,'chr'],'CEN']
        var2_end=cnv_sample[variant,'end']
        type <- cnv_sample[variant,'type']
        chr_data[chr_arm,type] = chr_data[chr_arm,type]+abs(var2_start-var2_end)
      }
    }
  }
}




