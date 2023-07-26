# Model 3
# 3-allele Model w/ Wild-Type, Male-Sterile t-Haplotype, and Lethal t-Haplotype
# With Reproductive Compensation + No Inbreeding

.libPaths('/path/to/Rlibs4')
library(tidyr)
library(dplyr)
library(rlist)

source('/path/to/T_Haplotye_Functions.R')

# Suppress summarise info
options(dplyr.summarise.inform = FALSE)

args = commandArgs(trailingOnly=TRUE)
s <- 0
t <- as.numeric(args[1])
k <- as.numeric(args[2])
C <- as.numeric(args[3])

epsilon <- 1/1000

# Initialize Starting Genotype Frequencies
gen <- 1
AA_current <- 1 - epsilon
Aa_current <- epsilon
aa_current <- 0
Aal_current <- 0
aal_current <- 0
alal_current <- 0
a_allele_freq_current <- (0.5*Aa_current) + (0.5*aal_current) + aa_current
  
# Generate Empty Table for Results Each Generation
genotype_frequencies <- matrix(0,nrow=100000,ncol=11) 
genotype_frequencies[gen,] <- c(s,t,k,C,gen,AA_current,Aa_current,aa_current,Aal_current,aal_current,alal_current)
    
a_freq_dif <- 1

# Loop through until an equilibrium frequency for the male-sterile allele is reached
while (gen < 100 | a_freq_dif > (1/10000)) {
  # Generate Mating Table from Current Genotype Frequencies
	current_mating_table <- three_allele_generate_mating_table_with_compensation_random_mating(AA_current,Aa_current,aa_current,Aal_current, aal_current, alal_current,k,s,t,C)

  # From Mating Table Get Genotype Frequencies In The Next Generation
  genotype_freq_next_gen <- three_allele_get_genotype_freqs_nextgen(mating_table=current_mating_table)
    
	AA_nextgen <- as.numeric(genotype_freq_next_gen[1])
	Aa_nextgen <- as.numeric(genotype_freq_next_gen[2])
	aa_nextgen <- as.numeric(genotype_freq_next_gen[3])
	Aal_nextgen <- as.numeric(genotype_freq_next_gen[4])
	aal_nextgen <- as.numeric(genotype_freq_next_gen[5])
	alal_nextgen <- as.numeric(genotype_freq_next_gen[6])
    
  gen <- gen + 1
    
	genotype_frequencies[gen,] <- c(s,t,k,C,gen,AA_nextgen,Aa_nextgen,aa_nextgen,Aal_nextgen,aal_nextgen,alal_nextgen)
	a_allele_freq_nextgen <- (0.5 * Aa_nextgen) + (0.5 * aal_nextgen) + aa_nextgen
	a_freq_dif <- abs(a_allele_freq_nextgen - a_allele_freq_current)
  
  # Update Genotype Frequencies
	AA_current <- AA_nextgen
	Aa_current <- Aa_nextgen
	aa_current <- aa_nextgen
	Aal_current <- Aal_nextgen
	aal_current <- aal_nextgen
	alal_current <- alal_nextgen
	
	a_allele_freq_current <- a_allele_freq_nextgen
}

## Introduce Lethal Distorter
# Only proceed if male-sterile is stably polymorphic
# A > 0.01 and a > 0.01

if ( (a_allele_freq_current > 0.01) & ( (1-a_allele_freq_current) > 0.01) ) {
  introduced_al_check <- TRUE 

  # Introduce lethal (al) allele in heterozygotes
  al_introduction_gen <- gen
  al_gen_wiggle <- al_introduction_gen + 100
    
  if(Aa_current > 0.05){
    Aa_current <- Aa_current - 0.05
    Aal_current <- 0.05
  } else {
    Aal_current <- Aa_current/2
    Aa_current <- Aa_current/2      
  }
    
  al_allele_freq_current <- (0.5*Aal_current) + (0.5*aal_current) + alal_current
    
  al_freq_dif <- 1  

  # Loop through until an equilibrium frequency for the male-sterile allele is reached
  # of allele freq is sufficiently small (<0.001) to imply eventual loss of lethal

  while ( (gen < al_gen_wiggle) | ((al_freq_dif > (1.0e-08)) & (al_allele_freq_current > 0.005)) ) {
    # Generate Mating Table from Current Genotype Frequencies
    current_mating_table <- three_allele_generate_mating_table_with_compensation_random_mating(AA_current,Aa_current,aa_current,Aal_current, aal_current, alal_current,k,s,t,C)
    
    # From Mating Table Get Genotype Frequencies In The Next Generation
    genotype_freq_next_gen <- three_allele_get_genotype_freqs_nextgen(mating_table=current_mating_table)
      
    AA_nextgen <- as.numeric(genotype_freq_next_gen[1])
    Aa_nextgen <- as.numeric(genotype_freq_next_gen[2])
    aa_nextgen <- as.numeric(genotype_freq_next_gen[3])
    Aal_nextgen <- as.numeric(genotype_freq_next_gen[4])
    aal_nextgen <- as.numeric(genotype_freq_next_gen[5])
    alal_nextgen <- as.numeric(genotype_freq_next_gen[6])
      
    gen <- gen + 1
      
    genotype_frequencies[gen,] <- c(s,t,k,C,gen,AA_nextgen,Aa_nextgen,aa_nextgen,Aal_nextgen,aal_nextgen,alal_nextgen)
    al_allele_freq_nextgen <- (0.5 * Aal_nextgen) + (0.5 * aal_nextgen) + alal_nextgen
    al_freq_dif <- abs(al_allele_freq_nextgen - al_allele_freq_current)
    
    # Update Genotype Frequencies
    AA_current <- AA_nextgen
    Aa_current <- Aa_nextgen
    aa_current <- aa_nextgen
    Aal_current <- Aal_nextgen
    aal_current <- aal_nextgen
    alal_current <- alal_nextgen
      
    al_allele_freq_current <- al_allele_freq_nextgen            
  }
} else {
	introduced_al_check <- FALSE
}

# Remove Extra Generation Rows
# add col info letting us know if lethal was introduced
genotype_frequencies <- genotype_frequencies[genotype_frequencies[,5]!=0,]
genotype_frequencies <- data.frame(genotype_frequencies) %>% 
dplyr::mutate(introduce_al=introduced_al_check)

# Prepare Results Output File
colnames(genotype_frequencies) <- c("s","t","k","C","gen","AA_freq","Aa_freq","aa_freq",'Aal_freq','aal_freq','alal_freq','introduce_al')

filename <- paste('Model3_t=',t,'_k=',k,'_C=',C,'.txt',sep='')

output_path <- paste('/path/to/Model3/results/',filename,sep='')
write.table(genotype_frequencies,file=output_path,col.names=T,row.names=F,quote=F)