#"`-:-.   ,-;"`-:-.   ,-;"`-:-.   ,-;"`-:-. -------------------------- ,-;"`-:-.
# `=`,'=/     `=`,'=/     `=`,'=/     `=`|  "If I have seen further |=/     `=
#   y==/        y==/        y==/        y| it is by standing on the |/
# ,=,-<=`.    ,=,-<=`.    ,=,-<=`.    ,=,|    shoulders of Giants"  |=`.    ,=
#,-'-'   `-=_,-'-'   `-=_,-'-'   `-=_,-'-'  -------------------------- `-=_,-'-'
# ==============================================================================
# PROJECT: Hybrid_Dynamics
# PURPOSE: to determine the linear fixed model that predicts hybridity
# DATE: 18 Dec 2024
# Devin P. Bendixsen, PhD
# Staff Bioinformatician | University of Edinburgh
# MRC Human Genetics Unit | Institute of Genetics and Cancer
# ==============================================================================

# ==============================================================================
# LOAD NEEDED PACKAGES
# ==============================================================================
library(nlme)
library(lme4)
library(MuMIn)
genome_size=12071326
require(plyr)

# ==============================================================================
# LOAD AND FORMAT DATA
# ==============================================================================
viability <- read_excel('data/spore_viability/Tetrad dissections spore viability.xlsx',sheet='Table_F1_Parents')
viability <- viability %>%
  subset(select=c('env','Replicate','gen','Proportion','Generation')) %>%
  dplyr::rename('rep'='Replicate','viability'='Proportion')
sv_data <- read.table(file='data/Hybrid_Dynamics_SV_data.txt')
sv_data <- sv_data %>%
  subset(select=c('env','gen','rep','tot_burden')) %>%
  dplyr::rename('SV'='tot_burden')

snv_data <- read.table(file='data/Hybrid_Dynamics_SNV_data.txt')
snv_data <- snv_data %>%
  subset(select=c('env','gen','rep','high')) %>%
  dplyr::rename('SNV'='high')

cnv_data <- read.table(file='data/Hybrid_Dynamics_CNV_data.txt')
cnv_data <- cnv_data %>%
  subset(select=c('env','gen','rep','total')) %>%
  dplyr::rename('CNV'='total')
cnv_data$CNV <- cnv_data$CNV/genome_size

data <- join_all(list(viability,snv_data,sv_data,cnv_data), by = c('env','gen','rep'), type='inner')
data$rep <- paste0(data$env,'_',data$rep)

# ==============================================================================
# 
# ==============================================================================
##  Upload the chromosome resolution data

fit <- lmer(viability ~ env*SNV*SV*CNV + (1|rep), data=data, REML=F,na.action=na.fail)
summary(fit)
models <- dredge(fit,rank="AIC")
print(summary(get.models(models, 1)[[1]]))
?dredge
final_model <- lmer(viability ~ CNV + env + SV + (1 | rep) + env:SV,data=data,REML=F,na.action=na.fail)
summary(final_model)
library(sjPlot)
install.packages('sjPlot')
sjPlot::plot_model(final_model)

model1 <- lmer(viability ~ env*SV*CNV*SNV + (1|rep), data=data, REML=F,na.action=na.fail)
summary(model1)
model2 <- lmer(viability ~ env*SV*CNV + (1|rep), data=data, REML=F,na.action=na.fail)
summary(model2)
model3 <- lmer(viability ~ env*SV*SNV + (1|rep), data=data, REML=F,na.action=na.fail)
summary(model3)
model4 <- lmer(viability ~ env*CNV*SNV + (1|rep), data=data, REML=F,na.action=na.fail)
summary(model4)
model5 <- lmer(viability ~ SV*SNV*CNV + (1|rep), data=data, REML=F,na.action=na.fail)
summary(model5)
model6 <- lmer(viability ~ env*SV + (1|rep), data=data, REML=F,na.action=na.fail)
summary(model6)
model7 <- lmer(viability ~ env*CNV + (1|rep), data=data, REML=F,na.action=na.fail)
summary(model7)
model8 <- lmer(viability ~ SV*CNV + (1|rep), data=data, REML=F,na.action=na.fail)
summary(model8)

# ==============================================================================
# 
# ==============================================================================
NaCl = '#4472c4ff'
LiAc0.01 = '#9e49e1ff'
LiAc0.02 = '#00b050ff'
Ethanol = '#ff0000ff'
H1F1 = '#808080'
H2F1 = '#808080'
H3F1 = '#808080'
F2 = '#d4a373'
ggplot(data, aes(x=SNV,y=viability)) +
  geom_point(shape=21,aes(size=Generation,colour=env))+
  scale_colour_manual(values=c(Ethanol,H1F1,H2F1,H3F1,LiAc0.01,LiAc0.02,NaCl))

ggplot(data, aes(x=SV,y=viability)) +
  geom_point(shape=21,aes(size=Generation,colour=env))+
  scale_colour_manual(values=c(Ethanol,H1F1,H2F1,H3F1,LiAc0.01,LiAc0.02,NaCl))

ggplot(data, aes(x=CNV*100,y=viability)) +
  geom_point(shape=21,aes(size=Generation,colour=env))+
  scale_colour_manual(values=c(Ethanol,H1F1,H2F1,H3F1,LiAc0.01,LiAc0.02,NaCl))
ggplot

ggplot(data, aes(x=CNV*100,y=SV)) +
  geom_point(shape=21,aes(size=Generation,colour=env))+
  scale_colour_manual(values=c(Ethanol,H1F1,H2F1,H3F1,LiAc0.01,LiAc0.02,NaCl))



data2 <- join_all(list(snv_data,sv_data,cnv_data), by = c('env','gen','rep'), type='inner')
ggplot(data2, aes(x=SNV,y=SV)) +
  geom_point(shape=21)
