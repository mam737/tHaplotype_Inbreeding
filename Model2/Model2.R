# Model 2
# 2-allele Model w/ Wild-Type + Male-Sterile t-Haplotype
# No Reproductive Compensation + With Sib-Mating

.libPaths('/path/to/Rlibs4')library(tidyr)
library(dplyr)
library(rlist)

source('/path/to/T_Haplotye_Functions.R')

# Suppress summarise info
options(dplyr.summarise.inform = FALSE)

args = commandArgs(trailingOnly=TRUE)
s <- 0.0
t <- as.numeric(args[1])
k <- as.numeric(args[2])
f <- as.numeric(args[3])

epsilon <- 1/1000

# Initialize Starting Genotype Frequencies
gen <- 1
AA_current <- 1 - epsilon
Aa_current <- epsilon
aa_current <- 0
a_allele_freq_current <- 0.5*Aa_current + aa_current

# Generate Empty Table for Results Each Generation
genotype_frequencies <- matrix(0,nrow=100000,ncol=8) 
genotype_frequencies[gen,] <- c(s,t,k,f,gen,AA_current,Aa_current,aa_current)

a_freq_dif <- 1

# Loop through until an equilibrium frequency for the male-sterile allele is reached
while (gen < 100 | a_freq_dif > (1/10000)) {
	if (gen==1) {
        # In first generation, we can only generate
        # the mating table under random mating
		current_mating_table <- two_allele_generate_mating_table_with_random_mating(AA_current,Aa_current,aa_current,k,s,t)
	} else {
        # Generate Mating Table from Current Genotype Frequencies
        # Proportion of Sib-Mating Relies on Knowing Prob Brother/Sis Genos
		current_mating_table <- two_allele_generate_mating_table_with_inbreeding(AA_current,Aa_current,aa_current,k,s,t,f,prob_brothergeno_given_female_offspring_geno=current_prob_brothergeno_given_female_offspring_geno)
	}
	
    # From Mating Table Get Genotype Frequencies In The Next Generation
	genotype_freq_next_gen <- two_allele_get_genotype_freqs_nextgen(mating_table=current_mating_table)

    # From Offspring Frequencies Determine
    # the probability of brother's geno given female geno
	current_prob_brothergeno_given_female_offspring_geno <- two_allele_calculate_inbreeding_genotype_probabilities(mating_table=current_mating_table,geno_freq_next_gen=genotype_freq_next_gen)

    # Update Genotype Frequencies
    AA_nextgen <- as.numeric(genotype_freq_next_gen[1])
    Aa_nextgen <- as.numeric(genotype_freq_next_gen[2])
    aa_nextgen <- as.numeric(genotype_freq_next_gen[3])
    gen <- gen + 1

    genotype_frequencies[gen,] <- c(s,t,k,f,gen,AA_nextgen,Aa_nextgen,aa_nextgen)

    a_allele_freq_nextgen <- (0.5 * Aa_nextgen) + aa_nextgen
    a_freq_dif <- abs(a_allele_freq_nextgen - a_allele_freq_current)

    AA_current <- AA_nextgen
    Aa_current <- Aa_nextgen
    aa_current <- aa_nextgen

    a_allele_freq_current <- a_allele_freq_nextgen
}
# Remove Extra Generation Rows
genotype_frequencies <- genotype_frequencies[genotype_frequencies[,5]!=0,]

# Prepare Results Output File
colnames(genotype_frequencies) <- c("s","t","k",'f',"gen","AA_freq","Aa_freq","aa_freq")

filename <- paste('Model2=',t,'_k=',k,'_f=',f,'.txt',sep='')

output_path <- paste('/path/to/Model2/results/',filename,sep='')
write.table(genotype_frequencies,file=output_path,col.names=T,row.names=F,quote=F)

