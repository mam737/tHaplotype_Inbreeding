## Script Used to Generate Figures in the Manuscript

.libPaths('/path/to/Rlibs4')

library(dplyr)
library(stringr)
library(data.table)
library(tidyr)
library(rlist)
library(ggplot2)
library(ggpubr)
library(MetBrewer)
library(forcats)
library(viridis)

source('/path/to/generate_figures_functions.R')

# List out all possible values for each parameter
s_param <- seq(0,1,by=0.1)
t_param <- seq(0,1,by=0.1)
k_param <- sort(c(seq(0.5,1,by=0.1),0.99))
C_param <- seq(0,1,by=0.1)
f_param <- c(0.0000,0.0001,0.0010,0.0100,0.0250,0.0500,0.0750,0.1000,0.2000,0.4000,0.5000)

# Model 1
# 2-allele: wild-type + male-sterile
Model1_results_dir <- '/path/to/Model1/results'
Model1_results <- extract_results(Model1_results_dir)
Model1_faf <- obtain_faf(Model_results,model_type='Model1')

# Model 2
# 2-allele: wild-type + male-sterile w/ inbreeding
Model2_results_dir <- '/home/brandvai/mmunasin/LethalDistorter_Inbreeding/final_models/COI2_model/results'
Model2_results <- extract_results(Model2_results_dir)
Model2_faf <- obtain_faf(Model2_results,model_type='Model2')

# Model 3
# 3-allele: wild-type, male-sterile and lethal 
# with reproductive compensation without inbreeding
Model3_results_dir <- '/path/to/Model3/results'
Model3_results <- extract_results(Model3_results_dir)
Model3_faf <- obtain_faf(Model3_results,model_type='Model3')

# Model 4
# 3-allele: wild-type, male-sterile and lethal 
# with reproductive compensation with inbreeding
Model4_results_dir <- '/path/to/Model4/results'
Model4_results <- extract_results(Model4_results_dir)
Model4_faf <- obtain_faf(Model4_results,model_type='Model4')

# Special Case 1 - 
# Adapt Model 1 to introduce ONLY Lethal t into naive pop
# No reproductive compensation or inbreeding
# This is analogous to Dunn and Levene 1961
# 2-allele: wild-type + lethal
DL1961_results_dir <- '/path/to/ModelDL1961/results'
DL1961_results <- extract_results(DL1961_results_dir)
DL1961_faf <- obtain_faf(DL1961_results,model_type='DL1961')

# Special Case 2 -
# Adapt Model 2 to introduce ONLY Lethal t into naive pop w/ sib-mating
# Inbreeding via sib-mating but no reproductive compensation
# This is analogous to Petras 1967
# 2-allele: wild-type + lethal w/ inbreeding
P1967_results_dir <- '/path/to/P1967_model/results'
P1967_results <- extract_results(P1967_results_dir)
P1967_faf <- obtain_faf(P1967_results,model_type='P1967')

##### FIGURE 1a #####
Fig1_a <- Model2_faf %>% filter(t==1.0) %>% 
filter(f %in% c(0,0.0001,0.0100,0.050,0.1,0.2,0.4,0.5)) %>%
dplyr::mutate(f=factor(f,levels=sort(c(0,0.0001,0.0100,0.050,0.1,0.2,0.4,0.5)))) %>%
ggplot(aes(x=k,y=a_freq,group=f)) +
geom_line(aes(color=f)) +
geom_point(aes(color=f))+
scale_color_manual(values=met.brewer('Tam',n=8)) + theme_bw()+
labs(x='k: drive strength',y='Distorter Final Frequency',fill='f',main='s = 0, t = 1') +
theme(legend.position='bottom')
ggsave(Fig1_a,file='/path/to/Rplots//Fig1_a.pdf',height=4,width=3.5)

##### FIGURE 1b #####
Fig1_b <- Model2_faf %>% filter(f==0.1) %>%
filter(t %in% c(0,0.2,0.4,0.6,0.8,1.0)) %>%
dplyr::mutate(t=factor(t,levels=seq(0,1,by=0.2))) %>%
ggplot(aes(x=k,y=a_freq,group=t)) +
geom_line(aes(color=t)) +
geom_point(aes(color=t))+
scale_color_manual(values=met.brewer('Hokusai3',n=6)) + theme_bw()+
theme(axis.text.x=element_text(size=8,angle=90,vjust=0.5,hjust=1))+
labs(x='k: drive strength',y='Male Sterile Distorter Final Frequency',fill='t',main='s = 0, f = 0.1')+
theme(legend.position='bottom')
ggsave(Fig1_b,file='/path/to/Rplots//Fig1_b.pdf',height=4,width=3.5)

##### FIGURE 2 #####
Fig2 <- Model4_faf %>% filter(t==1.0) %>%
filter(k %in% c(0.9,0.99)) %>%
dplyr::mutate(k=factor(k,levels=c(0.99,0.9))) %>% 
filter(C %in% c(0.0,0.2,0.4,0.6,0.8,1.0)) %>%
dplyr::mutate(C=factor(C,levels=c(0.0,0.2,0.4,0.6,0.8,1.0))) %>%
ungroup() %>% 
select(k,f,C,gen,a_freq,updated_al_freq,t_freq) %>%
pivot_longer(cols=c('updated_al_freq','a_freq','t_freq'),values_to='allele_freq',names_to='allele_type') %>%
dplyr::mutate(allele_type=case_when(
	allele_type=='a_freq' ~ 'Male Sterile',
	allele_type=='updated_al_freq' ~ "Lethal",
	TRUE ~ 't-Haplotype')) %>%
dplyr::mutate(allele_type=factor(allele_type,levels=c('Male Sterile','Lethal','t-Haplotype'))) %>%
ggplot(aes(x=f,y=allele_freq,color=C))+
geom_line(aes(color=C))+
geom_point(aes(color=C))+
facet_grid(k~allele_type)+
scale_y_continuous(breaks=seq(0,1,by=0.2),limits=c(0,1.0))+
scale_color_manual(values=met.brewer('VanGogh3',n=6)) + theme_bw()+
labs(x='p_sib: Sib-Mating Proportion',y='Final Frequency',fill='f',main='s = 0, t = 1')+
theme(axis.text.x=element_text(size=8,angle=90,vjust=0.5,hjust=1))
ggsave(Fig2,file='/path/to/Rplots//Fig2.pdf',height=4,width=6)

##### FIGURE 3 #####
Fig3 <- Model4_faf %>% filter(s==0) %>% filter(t==1.0) %>%
dplyr::mutate(f=factor(f,levels=f_param)) %>%
dplyr::mutate(C=factor(C,levels=C_param)) %>%
filter(k %in% c(0.9,0.99)) %>%
dplyr::mutate(k=factor(k,levels=c(0.99,0.9))) %>% ungroup() %>%
select(k,f,C,a_freq,updated_al_freq,t_freq) %>%
dplyr::mutate(lethal_prop=updated_al_freq/t_freq) %>%
ggplot(aes(x=f,y=C,fill=lethal_prop))+
geom_tile()+
facet_wrap(~k,ncol=2,nrow=1)+
scale_fill_gradientn(colors=met.brewer('OKeeffe2'),na.value = 'grey90')+ theme_bw()+
labs(x='psib = Prop Sib Mating',y='C: compensation strength',fill='Lethal Prop',main='s = 0, t = 1')+
theme(axis.text.x=element_text(size=8,angle=90,vjust=0.5,hjust=1))
ggsave(Fig3,file='/path/to/Rplots//Fig3.pdf',height=3,width=6)

##### FIGURE S1 #####
FigS1 <- Model3_faf %>% filter(t==1.0) %>% 
filter(k %in% c(0.90,0.99)) %>%
dplyr::mutate(k=factor(k,levels=c(0.99,0.9))) %>% 
ungroup()%>%
select(t, k, C, updated_al_freq,t_freq) %>% 
pivot_longer(cols=c('updated_al_freq','t_freq'),values_to='allele_freq',names_to='allele_type') %>%
dplyr::mutate(allele_type=case_when(
	allele_type=='t_freq' ~ 't-Haplotype',
	TRUE ~ 'Lethal Haplotype')) %>%
dplyr::mutate(allele_type=factor(allele_type,levels=c('t-Haplotype','Lethal Haplotype')))%>%
ggplot(aes(x=C,y=allele_freq,group=allele_type))+
geom_line(aes(color=allele_type))+
geom_point(aes(color=allele_type,shape=allele_type))+
facet_wrap(~k,ncol=1)+
scale_color_manual(values=c('#BDC696','#13140f')) + 
scale_y_continuous(breaks=seq(0,1,by=0.2),limits=c(0,1.0))+
labs(x='C: Reproductive Compensation',y='Allele Frequency',main='Model 3: t = 1.0') +
theme(legend.position='right')
ggsave(FigS1,file='/path/to/Rplots/FigS1.pdf',height=4,width=4)

##### FIGURE S2 #####
FigS2 <- DL1961_faf %>% filter(t==1.0) %>%
filter(s==0) %>% select(k,al_freq) %>%
# Expected lethal t frequency given k from Dunn and Levene 1961
dplyr::mutate(DL_al_freq= (0.5-(sqrt(k * (1-k))/(2*k)))) %>%
pivot_longer(cols=c('al_freq','DL_al_freq'),values_to='allele_freq',names_to='allele_type') %>%
dplyr::mutate(allele_type=case_when(
	allele_type=='al_freq' ~ 'Observed t-Haplotype',
	TRUE ~ 'Predicted t-Haplotype')) %>%
ggplot(aes(x=k,y=allele_freq,group=allele_type))+
geom_line(aes(color=allele_type))+
geom_point(aes(shape=allele_type,color=allele_type))+
scale_y_continuous(breaks=seq(0,1,by=0.2),limits=c(0,1.0))+
scale_color_manual(values=c('#BDC696','#393b2d')) + 
labs(x='Segregation Distortion',y='Allele frequency',main='s=0, t= 1: Observed vs Predicted')
ggsave(FigS2,file='/path/to/Rplots/FigS2.pdf',height=2,width=4)

##### FIGURE S3 #####
FigS3 <- P1967_faf %>% 
filter(t==1.0) %>% 
filter(k==0.9) %>% 
select(s,t,k,f,C,al_freq) %>% 
# See supplement for how we derive transformation between
# our parameter f representing prop of sib mating
# and Petras 1967 inbreeding coefficient 
dplyr::mutate(transformed_f=f/(4- (3*f) )) %>% 
# Given transformed f we can calculate what expected
# lethal t frequency should be from Petras 1967
dplyr::mutate(petras_al_freq2=petras_prediction(m=k,F=transformed_f,type='minus')) %>%
pivot_longer(cols=c('al_freq','petras_al_freq2'),values_to='allele_freq',names_to='allele_type') %>%
dplyr::mutate(allele_type=case_when(
	allele_type=='al_freq' ~ 'Observed',
	TRUE ~ 'Predicted')) %>%
ggplot(aes(x=f,y=allele_freq,group=allele_type))+
geom_line(aes(color=allele_type))+
geom_point(aes(shape=allele_type,color=allele_type))+
scale_color_manual(values=c('#BDC696','#13140f')) + 
labs(x='Proportion of Sib-Mating',y='Allele frequency',main='s=0, t= 1,k=0.9: Observed vs Predicted',fill='Type')+
theme(legend.position="right")
ggsave(FigS3,file='/path/to/Rplots/FigS3.pdf',height=3,width=4)