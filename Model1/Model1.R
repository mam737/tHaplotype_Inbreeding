# Model 1
# 2-allele Model w/ Wild-Type + Male-Sterile t-Haplotype
# No Reproductive Compensation + No Inbreeding

.libPaths('/path/to/Rlibs4')
library(tidyr)
library(dplyr)
library(rlist)

source('/path/to/T_Haplotye_Functions.R')

# Suppress summarise info
options(dplyr.summarise.inform = FALSE)

args = commandArgs(trailingOnly=TRUE)
s <- as.numeric(args[1])
t <- as.numeric(args[2])
k <- as.numeric(args[3])

epsilon <- 1/1000

# Initialize Starting Genotype Frequencies
gen <- 1
AA_current <- 1 - epsilon
Aa_current <- epsilon
aa_current <- 0
a_allele_freq_current <- 0.5*Aa_current + aa_current

# Generate Empty Table for Results Each Generation
genotype_frequencies <- matrix(0,nrow=100000,ncol=7) 
genotype_frequencies[gen,] <- c(s,t,k,gen,AA_current,Aa_current,aa_current)

a_freq_dif <- 1

# Loop through until an equilibrium frequency for the male-sterile allele is reached
while (gen < 100 | a_freq_dif > (1/10000)) {
	# Generate Mating Table from Current Genotype Frequencies
    current_mating_table <- two_allele_generate_mating_table_with_random_mating(AA_current,Aa_current,aa_current,k,s,t)

    # From Mating Table Get Genotype Frequencies In The Next Generation
	genotype_freq_next_gen <- two_allele_get_genotype_freqs_nextgen(mating_table=current_mating_table)

    AA_nextgen <- as.numeric(genotype_freq_next_gen[1])
    Aa_nextgen <- as.numeric(genotype_freq_next_gen[2])
    aa_nextgen <- as.numeric(genotype_freq_next_gen[3])
    gen <- gen + 1

    genotype_frequencies[gen,] <- c(s,t,k,gen,AA_nextgen,Aa_nextgen,aa_nextgen)

    a_allele_freq_nextgen <- (0.5 * Aa_nextgen) + aa_nextgen
    a_freq_dif <- abs(a_allele_freq_nextgen - a_allele_freq_current)

    # Update Genotype Frequencies
    AA_current <- AA_nextgen
    Aa_current <- Aa_nextgen
    aa_current <- aa_nextgen

    a_allele_freq_current <- a_allele_freq_nextgen
}

# Remove Extra Generation Rows
genotype_frequencies <- genotype_frequencies[genotype_frequencies[,4]!=0,]

# Prepare Results Output File
colnames(genotype_frequencies) <- c("s","t","k","gen","AA_freq","Aa_freq","aa_freq")

filename <- paste('Model1_s=',s,'_t=',t,'_k=',k,'.txt',sep='')

output_path <- paste('/path/to/Model1/results/',filename,sep='')
write.table(genotype_frequencies,file=output_path,col.names=T,row.names=F,quote=F)