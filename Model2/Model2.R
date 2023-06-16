.libPaths('/home/brandvai/mmunasin/Rlibs4')
library(tidyr)
library(dplyr)
library(rlist)

source('/home/brandvai/mmunasin/LethalDistorter_Inbreeding/final_models/T_Haplotye_Functions.R')

# Suppress summarise info
options(dplyr.summarise.inform = FALSE)

args = commandArgs(trailingOnly=TRUE)
s <- 0.6
t <- as.numeric(args[1])
k <- as.numeric(args[2])
f <- as.numeric(args[3])

epsilon <- 1/1000

# Initialize starting genotype frequencies
gen <- 1
AA_current <- 1 - epsilon
Aa_current <- epsilon
aa_current <- 0
a_allele_freq_current <- 0.5*Aa_current + aa_current

genotype_frequencies <- matrix(0,nrow=100000,ncol=8) 
genotype_frequencies[gen,] <- c(s,t,k,f,gen,AA_current,Aa_current,aa_current)

a_freq_dif <- 1

# See whether you can obtain an equilibrium for the distorter allele 
# without the lethal component
while (gen < 100 | a_freq_dif > (1/10000)) {
	if (gen==1) {
		current_mating_table <- two_allele_generate_mating_table_with_random_mating(AA_current,Aa_current,aa_current,k,s,t)
	} else {
		current_mating_table <- two_allele_generate_mating_table_with_inbreeding(AA_current,Aa_current,aa_current,k,s,t,f,prob_brothergeno_given_female_offspring_geno=current_prob_brothergeno_given_female_offspring_geno)
	}
	
	genotype_freq_next_gen <- two_allele_get_genotype_freqs_nextgen(mating_table=current_mating_table)

	current_prob_brothergeno_given_female_offspring_geno <- two_allele_calculate_inbreeding_genotype_probabilities(mating_table=current_mating_table,geno_freq_next_gen=genotype_freq_next_gen)

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

genotype_frequencies <- genotype_frequencies[genotype_frequencies[,5]!=0,]
colnames(genotype_frequencies) <- c("s","t","k",'f',"gen","AA_freq","Aa_freq","aa_freq")

filename <- paste('COI2_t=',t,'_k=',k,'_f=',f,'.txt',sep='')

output_path <- paste('/home/brandvai/mmunasin/LethalDistorter_Inbreeding/final_models/COI2_model/results/',filename,sep='')
write.table(genotype_frequencies,file=output_path,col.names=T,row.names=F,quote=F)

