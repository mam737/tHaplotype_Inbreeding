################################################################################
################################################################################
##### TWO ALLELE MODELS ########################################################
##### THIS INCLUDES MODEL 1 AND MODEL 2 ########################################
################################################################################

# Consider a two allele model
# A = wild-type, a = distorter
# Generate the mating table under random mating
# Parameters of interest:
#    - k: probability a is transmitted in heterozygous males (drive strength)
#    - s: Fitness reduction in heterozygous Aa males
#    - t: Fitness reduction in homozygous aa males
two_allele_generate_mating_table_with_random_mating <- function(AA_current,Aa_current,aa_current,k,s,t) {
  mating_table <- crossing(female_geno=c("AA","Aa","aa"),male_geno=c("AA","Aa","aa")) %>%
    # Establish Female Genotype Frequencies
    dplyr::mutate(female_geno_freq=case_when(
      female_geno=='AA'~AA_current,
      female_geno=='Aa' ~ Aa_current,
      female_geno=='aa'~aa_current)) %>%
    # Establish Male Genotype Frequencies
    dplyr::mutate(male_geno_freq=case_when(
      male_geno=='AA'~ AA_current,
      male_geno=='Aa' ~ Aa_current,
      male_geno=='aa'~ aa_current)) %>%
    # Establish Mating Type Frequencies
    dplyr::mutate(mating_freq=female_geno_freq*male_geno_freq) %>%
    # Establish Frequencies of A allele in gametes (no drive in females)
    dplyr::mutate(female_A_gam_freq=case_when(
      female_geno=='AA' ~ 1,
      female_geno=='Aa' ~ 0.5,
      female_geno=='aa' ~ 0)) %>%
    # Establish Frequencies of a allele in gametes (no drive in females)
    dplyr::mutate(female_a_gam_freq=case_when(
      female_geno=='AA' ~ 0,
      female_geno=='Aa' ~ 0.5,
      female_geno=='aa' ~ 1)) %>%
    # Establish Frequencies of A allele in gametes (a allele drives in males)
    dplyr::mutate(male_A_gam_freq=case_when(
      male_geno=='AA' ~ 1,
      male_geno=='Aa' ~ (1-k),
      male_geno=='aa' ~ 0)) %>%
    # Establish Frequencies of a allele in gametes (a allele drives in males)
    dplyr::mutate(male_a_gam_freq=case_when(
      male_geno=='AA' ~ 0,
      male_geno=='Aa' ~ k,
      male_geno=='aa' ~ 1)) %>%
    # From each mating type, establish genotype freqs of offspring
    dplyr::mutate(
      freq_AA_geno=female_A_gam_freq * male_A_gam_freq,
      freq_Aa_geno=female_A_gam_freq*male_a_gam_freq + female_a_gam_freq*male_A_gam_freq,
      freq_aa_geno=female_a_gam_freq*male_a_gam_freq) %>%
    # Establish Family Fitness of Mating Type, here dependent on male fitness
    dplyr::mutate(fam_fit=case_when(
      male_geno=='AA' ~ 1,
      male_geno=='Aa' ~ (1-s),
      male_geno=='aa' ~ (1-t))) %>%
    # Determine Mean Population Fitness
    dplyr::mutate(mean_fit=sum(mating_freq*fam_fit)) %>%
    # Normalize Family Fitness by Mean Pop Fitness
    dplyr::mutate(rel_fam_fit=fam_fit/mean_fit) %>%
    # Extract Freq of Family After Selection
    dplyr::mutate(fam_freq_as=rel_fam_fit*mating_freq)
    return(mating_table)
}

# Consider a two allele model
# A = wild-type, a = distorter
# Pull genotype frequencies out from mating table
two_allele_get_genotype_freqs_nextgen <- function(mating_table) {
  geno_freq_next_gen <- mating_table %>% 
    dplyr::summarise(
      AA_ng=sum(freq_AA_geno*fam_freq_as),
      Aa_ng=sum(freq_Aa_geno*fam_freq_as),
      aa_ng=sum(freq_aa_geno*fam_freq_as))
  return(geno_freq_next_gen)  
}

# Consider a two allele model
# A = wild-type, a = distorter
# Determine probabilities given genotype frequencies
# for brother's genotype given sister's genotype
# this helps us determine different mating frequencies given inbreeding
# We consider inbreeding via sib-mating in Model 2
two_allele_calculate_inbreeding_genotype_probabilities <- function(mating_table,geno_freq_next_gen) {
  # Get probabilities offspring came from a given mating type given
  # offspring genotype
  mating_given_geno <- mating_table %>%
    dplyr::mutate(pmAA = freq_AA_geno * fam_freq_as / pull(geno_freq_next_gen,AA_ng),
                  pmAa = freq_Aa_geno * fam_freq_as / pull(geno_freq_next_gen,Aa_ng),
                  pmaa = freq_aa_geno * fam_freq_as / pull(geno_freq_next_gen,aa_ng)) %>%
    dplyr::mutate(pmAA=case_when(is.na(pmAA) ~ 0, TRUE ~ pmAA)) %>%
    dplyr::mutate(pmAa=case_when(is.na(pmAa) ~ 0, TRUE ~ pmAa)) %>%
    dplyr::mutate(pmaa=case_when(is.na(pmaa) ~ 0, TRUE ~ pmaa)) %>%
    dplyr::select(female_geno,male_geno,freq_AA_geno,freq_Aa_geno,freq_aa_geno,pmAA,pmAa,pmaa) %>% 
    dplyr::mutate(MT=paste(female_geno,male_geno,sep='_'))
  
  # Probability of a female offspring's brother having a given genotype
  # given the female offspring's genotype
  prob_brothergeno_given_female_offspring_geno <- mating_given_geno %>% 
    dplyr::select(-female_geno, - male_geno) %>%
    pivot_longer(cols = contains("freq"),values_to = "freq_male_offspring_geno",names_to = "male_offspring_geno")%>%
    separate(male_offspring_geno,sep='_',into=c('char_tag','male_offspring_geno','char_tag2'))  %>% select(-char_tag,-char_tag2) %>%
    pivot_longer(cols = contains("pm"),values_to = "prob_female_offspring_geno_fromMT",names_to = "female_offspring_geno") %>%
    separate(female_offspring_geno,sep='m',into=c('char_tag','female_offspring_geno')) %>% select(-char_tag) %>%
    group_by(female_offspring_geno,male_offspring_geno) %>%
    summarise(prob_brother_geno_given_female_offspring_geno = sum(prob_female_offspring_geno_fromMT *freq_male_offspring_geno))
  
  return(prob_brothergeno_given_female_offspring_geno)
}

# Consider a two allele model
# A = wild-type, a = distorter
# Generate the mating table under random mating
# Parameters of interest:
#    - k: probability a is transmitted in heterozygous males (drive strength)
#    - s: Fitness reduction in heterozygous Aa males
#    - t: Fitness reduction in homozygous aa males
#    - f: Probability of sib-mating
# We consider inbreeding via sib-mating in Model 2
two_allele_generate_mating_table_with_inbreeding <- function(AA_current,Aa_current,aa_current,k,s,t,f,prob_brothergeno_given_female_offspring_geno) {
  #Reset Colnames For Left-Join
  colnames(prob_brothergeno_given_female_offspring_geno) <- c('female_geno','male_geno','inbreeding_mating_prob')
  
  mating_table <- crossing(female_geno=c("AA","Aa","aa"),male_geno=c("AA","Aa","aa")) %>%
    # Establish Female Genotype Frequencies
    dplyr::mutate(female_geno_freq=case_when(
      female_geno=='AA'~AA_current,
      female_geno=='Aa' ~ Aa_current,
      female_geno=='aa'~aa_current)) %>%
    # Establish Male Genotype Frequencies
    dplyr::mutate(male_geno_freq=case_when(
      male_geno=='AA'~ AA_current,
      male_geno=='Aa' ~ Aa_current,
      male_geno=='aa'~ aa_current)) %>%
    # Random Mating Freq
    dplyr::mutate(mating_freq_random=female_geno_freq*male_geno_freq * (1-f)) %>%
    # Inbreed Mating Freq
    left_join(prob_brothergeno_given_female_offspring_geno,by=c('female_geno','male_geno')) %>%
    dplyr::mutate(mating_freq_inbred=inbreeding_mating_prob*f*female_geno_freq) %>%
    # Obtain Final Mating Freq
    dplyr::mutate(mating_freq=mating_freq_random+mating_freq_inbred) %>%
    # Establish Frequencies of A allele in gametes (no drive in females)
    dplyr::mutate(female_A_gam_freq=case_when(
      female_geno=='AA' ~ 1,
      female_geno=='Aa' ~ 0.5,
      female_geno=='aa' ~ 0)) %>%
    # Establish Frequencies of a allele in gametes (no drive in females)
    dplyr::mutate(female_a_gam_freq=case_when(
      female_geno=='AA' ~ 0,
      female_geno=='Aa' ~ 0.5,
      female_geno=='aa' ~ 1)) %>%
    # Establish Frequencies of A allele in gametes (a allele drives in males)
    dplyr::mutate(male_A_gam_freq=case_when(
      male_geno=='AA' ~ 1,
      male_geno=='Aa' ~ (1-k),
      male_geno=='aa' ~ 0)) %>%
    # Establish Frequencies of a allele in gametes (a allele drives in males)
    dplyr::mutate(male_a_gam_freq=case_when(
      male_geno=='AA' ~ 0,
      male_geno=='Aa' ~ k,
      male_geno=='aa' ~ 1)) %>%
    # From each mating type, establish genotype freqs of offspring
    dplyr::mutate(
      freq_AA_geno=female_A_gam_freq * male_A_gam_freq,
      freq_Aa_geno=female_A_gam_freq*male_a_gam_freq + female_a_gam_freq*male_A_gam_freq,
      freq_aa_geno=female_a_gam_freq*male_a_gam_freq) %>%
    # Establish Family Fitness of Mating Type, here dependent on male fitness
    dplyr::mutate(fam_fit=case_when(
      male_geno=='AA' ~ 1,
      male_geno=='Aa' ~ (1-s),
      male_geno=='aa' ~ (1-t))) %>%
    # Determine Mean Population Fitness
    dplyr::mutate(mean_fit=sum(mating_freq*fam_fit)) %>%
    # Normalize Family Fitness by Mean Pop Fitness
    dplyr::mutate(rel_fam_fit=fam_fit/mean_fit) %>%
    # Extract Freq of Family After Selection
    dplyr::mutate(fam_freq_as=rel_fam_fit*mating_freq)
  
  return(mating_table)
}
################################################################################
################################################################################
################################################################################
################################################################################

################################################################################
################################################################################
##### THREE ALLELE MODELS ######################################################
##### THIS INCLUDES MODEL 3 AND MODEL 4 ########################################
################################################################################
################################################################################

# Consider a three allele model
# A = wild-type, a = distorter, al = lethal distorter
# Generate the mating table under random mating
# WITHOUT COMPENSATION (i.e.,the benefit to female fitness by killing off alal homozygotes during embryogenesis)
# Parameters of interest:
#    - k: probability a is transmitted in heterozygous males (drive strength)
#    - s: Fitness reduction in heterozygous Aa males
#    - t: Fitness reduction in homozygous aa males
three_allele_generate_mating_table_with_random_mating <- function(AA_current,Aa_current,aa_current,Aal_current, aal_current, alal_current,k,s,t) {
  mating_table <- crossing(female_geno=c("AA","Aa","aa","Aal",'aal','alal'),male_geno=c("AA","Aa","aa",'Aal','aal','alal')) %>%
    # Establish Female Genotype Frequencies
    dplyr::mutate(female_geno_freq=case_when(
      female_geno=='AA'~AA_current,
      female_geno=='Aa' ~ Aa_current,
      female_geno=='aa'~aa_current,
      female_geno=='Aal' ~ Aal_current,
      female_geno=='aal' ~ aal_current,
      female_geno=='alal' ~ alal_current)) %>%
    # Establish Male Genotype Frequencies
    dplyr::mutate(male_geno_freq=case_when(
      male_geno=='AA'~AA_current,
      male_geno=='Aa' ~ Aa_current,
      male_geno=='aa'~aa_current,
      male_geno=='Aal' ~ Aal_current,
      male_geno=='aal' ~ aal_current,
      male_geno=='alal' ~ alal_current)) %>%
    # Establish Mating Type Frequencies
    dplyr::mutate(mating_freq=female_geno_freq*male_geno_freq) %>%
    # Establish Frequencies of A allele in gametes (no drive in females)
    dplyr::mutate(female_A_gam_freq=case_when(
      female_geno=='AA' ~ 1,
      female_geno=='Aa' ~ 0.5,
      female_geno=='aa' ~ 0,
      female_geno=='Aal' ~ 0.5,
      female_geno=='aal' ~ 0,
      female_geno=='alal' ~ 0)) %>%
    # Establish Frequencies of a allele in gametes (no drive in females)
    dplyr::mutate(female_a_gam_freq=case_when(
      female_geno=='AA' ~ 0,
      female_geno=='Aa' ~ 0.5,
      female_geno=='aa' ~ 1,
      female_geno=='Aal' ~ 0,
      female_geno=='aal' ~ 0.5,
      female_geno=='alal' ~ 0)) %>%
    # Establish Frequencies of al allele in gametes (no drive in females)
    dplyr::mutate(female_al_gam_freq=case_when(
      female_geno=='AA' ~ 0,
      female_geno=='Aa' ~ 0,
      female_geno=='aa' ~ 0,
      female_geno=='Aal' ~ 0.5,
      female_geno=='aal' ~ 0.5,
      female_geno=='alal' ~ 1)) %>%
    # Establish Frequencies of A allele in gametes (a allele drives in males)
    dplyr::mutate(male_A_gam_freq=case_when(
      male_geno=='AA' ~ 1,
      male_geno=='Aa' ~ (1-k),
      male_geno=='aa' ~ 0,
      male_geno=='Aal'~ (1-k),
      male_geno=='aal'~ 0,
      male_geno=='alal'~ 0)) %>%
    # Establish Frequencies of a allele in gametes (a allele drives in males)
    dplyr::mutate(male_a_gam_freq=case_when(
      male_geno=='AA' ~ 0,
      male_geno=='Aa' ~ k,
      male_geno=='aa' ~ 1,
      male_geno=='Aal'~ 0,
      male_geno=='aal'~ 0.5,
      male_geno=='alal'~ 0)) %>%
    # Establish Frequencies of al allele in gametes (a allele drives in males)
    dplyr::mutate(male_al_gam_freq=case_when(
      male_geno=='AA' ~ 0,
      male_geno=='Aa' ~ 0,
      male_geno=='aa' ~ 0,
      male_geno=='Aal'~ k,
      male_geno=='aal'~ 0.5,
      male_geno=='alal'~ 1)) %>%
    # From each mating type, establish genotype freqs of offspring
    dplyr::mutate(
      freq_AA_geno=female_A_gam_freq * male_A_gam_freq,
      freq_Aa_geno=female_A_gam_freq*male_a_gam_freq + female_a_gam_freq*male_A_gam_freq,
      freq_aa_geno=female_a_gam_freq*male_a_gam_freq,
      freq_Aal_geno=female_A_gam_freq*male_al_gam_freq + female_al_gam_freq*male_A_gam_freq,
      freq_aal_geno=female_a_gam_freq*male_al_gam_freq + female_al_gam_freq*male_a_gam_freq,
      freq_alal_geno=female_al_gam_freq * male_al_gam_freq) %>%
    # Establish Family Fitness of Mating Type, here dependent on male fitness
    dplyr::mutate(fam_fit=(1-freq_alal_geno)*case_when(
      male_geno=='AA' ~ 1,
      male_geno=='Aa' ~ (1-s),
      male_geno=='aa' ~ (1-t),
      male_geno=='Aal' ~ (1-s),
      male_geno=='aal' ~ (1-t),
      male_geno=='alal' ~ (1-t) )) %>%
    dplyr::mutate(updated_freq_alal_geno=0) %>% 
    dplyr::mutate(updated_freq_AA_geno=case_when(
      freq_alal_geno!=1 ~ freq_AA_geno/(1-freq_alal_geno),
      TRUE ~ 0)) %>%
    dplyr::mutate(updated_freq_Aa_geno=case_when(
      freq_alal_geno!=1 ~ freq_Aa_geno/(1-freq_alal_geno),
      TRUE ~ 0)) %>%
    dplyr::mutate(updated_freq_aa_geno=case_when(
      freq_alal_geno!=1 ~ freq_aa_geno/(1-freq_alal_geno),
      TRUE ~ 0)) %>%
    dplyr::mutate(updated_freq_Aal_geno=case_when(
      freq_alal_geno!=1 ~ freq_Aal_geno/(1-freq_alal_geno),
      TRUE ~ 0)) %>%
    dplyr::mutate(updated_freq_aal_geno=case_when(
      freq_alal_geno!=1 ~ freq_aal_geno/(1-freq_alal_geno),
      TRUE ~ 0)) %>%
    # Determine Mean Population Fitness
    dplyr::mutate(mean_fit=sum(mating_freq*fam_fit)) %>%
    # Normalize Family Fitness by Mean Pop Fitness
    dplyr::mutate(rel_fam_fit=fam_fit/mean_fit) %>%
    # Extract Freq of Family After Selection
    dplyr::mutate(fam_freq_as=rel_fam_fit*mating_freq)
  
  return(mating_table)
}

# Consider a three allele model
# A = wild-type, a = distorter, al = lethal distorter
# Generate the mating table under random mating
# WITH COMPENSATION (i.e.,the benefit to female fitness by killing off alal homozygotes during embryogenesis)
# Parameters of interest:
#    - k: probability a is transmitted in heterozygous males (drive strength)
#    - s: Fitness reduction in heterozygous Aa males
#    - t: Fitness reduction in homozygous aa males
#    - C: Strength of compensation
three_allele_generate_mating_table_with_compensation_random_mating <- function(AA_current,Aa_current,aa_current,Aal_current, aal_current, alal_current,k,s,t,C) {
  mating_table <- crossing(female_geno=c("AA","Aa","aa","Aal",'aal','alal'),male_geno=c("AA","Aa","aa",'Aal','aal','alal')) %>%
    # Establish Female Genotype Frequencies
    dplyr::mutate(female_geno_freq=case_when(
      female_geno=='AA'~AA_current,
      female_geno=='Aa' ~ Aa_current,
      female_geno=='aa'~aa_current,
      female_geno=='Aal' ~ Aal_current,
      female_geno=='aal' ~ aal_current,
      female_geno=='alal' ~ alal_current)) %>%
    # Establish Male Genotype Frequencies
    dplyr::mutate(male_geno_freq=case_when(
      male_geno=='AA'~AA_current,
      male_geno=='Aa' ~ Aa_current,
      male_geno=='aa'~aa_current,
      male_geno=='Aal' ~ Aal_current,
      male_geno=='aal' ~ aal_current,
      male_geno=='alal' ~ alal_current)) %>%
    # Establish Mating Type Frequencies
    dplyr::mutate(mating_freq=female_geno_freq*male_geno_freq) %>%
    # Establish Frequencies of A allele in gametes (no drive in females)
    dplyr::mutate(female_A_gam_freq=case_when(
      female_geno=='AA' ~ 1,
      female_geno=='Aa' ~ 0.5,
      female_geno=='aa' ~ 0,
      female_geno=='Aal' ~ 0.5,
      female_geno=='aal' ~ 0,
      female_geno=='alal' ~ 0)) %>%
    # Establish Frequencies of a allele in gametes (no drive in females)
    dplyr::mutate(female_a_gam_freq=case_when(
      female_geno=='AA' ~ 0,
      female_geno=='Aa' ~ 0.5,
      female_geno=='aa' ~ 1,
      female_geno=='Aal' ~ 0,
      female_geno=='aal' ~ 0.5,
      female_geno=='alal' ~ 0)) %>%
    # Establish Frequencies of al allele in gametes (no drive in females)
    dplyr::mutate(female_al_gam_freq=case_when(
      female_geno=='AA' ~ 0,
      female_geno=='Aa' ~ 0,
      female_geno=='aa' ~ 0,
      female_geno=='Aal' ~ 0.5,
      female_geno=='aal' ~ 0.5,
      female_geno=='alal' ~ 1)) %>%
    # Establish Frequencies of A allele in gametes (a allele drives in males)
    dplyr::mutate(male_A_gam_freq=case_when(
      male_geno=='AA' ~ 1,
      male_geno=='Aa' ~ (1-k),
      male_geno=='aa' ~ 0,
      male_geno=='Aal'~ (1-k),
      male_geno=='aal'~ 0,
      male_geno=='alal'~ 0)) %>%
    # Establish Frequencies of a allele in gametes (a allele drives in males)
    dplyr::mutate(male_a_gam_freq=case_when(
      male_geno=='AA' ~ 0,
      male_geno=='Aa' ~ k,
      male_geno=='aa' ~ 1,
      male_geno=='Aal'~ 0,
      male_geno=='aal'~ 0.5,
      male_geno=='alal'~ 0)) %>%
    # Establish Frequencies of al allele in gametes (a allele drives in males)
    dplyr::mutate(male_al_gam_freq=case_when(
      male_geno=='AA' ~ 0,
      male_geno=='Aa' ~ 0,
      male_geno=='aa' ~ 0,
      male_geno=='Aal'~ k,
      male_geno=='aal'~ 0.5,
      male_geno=='alal'~ 1)) %>%
    # From each mating type, establish genotype freqs of offspring
    dplyr::mutate(
      freq_AA_geno=female_A_gam_freq * male_A_gam_freq,
      freq_Aa_geno=female_A_gam_freq*male_a_gam_freq + female_a_gam_freq*male_A_gam_freq,
      freq_aa_geno=female_a_gam_freq*male_a_gam_freq,
      freq_Aal_geno=female_A_gam_freq*male_al_gam_freq + female_al_gam_freq*male_A_gam_freq,
      freq_aal_geno=female_a_gam_freq*male_al_gam_freq + female_al_gam_freq*male_a_gam_freq,
      freq_alal_geno=female_al_gam_freq * male_al_gam_freq) %>%
    # Establish Family Fitness of Mating Type, here dependent on male fitness
    dplyr::mutate(fam_fit=( ((1-freq_alal_geno) / (1-freq_alal_geno)^C))*case_when(
      male_geno=='AA' ~ 1,
      male_geno=='Aa' ~ (1-s),
      male_geno=='aa' ~ (1-t),
      male_geno=='Aal' ~ (1-s),
      male_geno=='aal' ~ (1-t),
      male_geno=='alal' ~ (1-t) )) %>%
    dplyr::mutate(fam_fit=case_when(is.na(fam_fit) ~ 0,TRUE ~ fam_fit)) %>%
    dplyr::mutate(updated_freq_alal_geno=0) %>% 
    dplyr::mutate(updated_freq_AA_geno=case_when(
      freq_alal_geno!=1 ~ freq_AA_geno/(1-freq_alal_geno),
      TRUE ~ 0)) %>%
    dplyr::mutate(updated_freq_Aa_geno=case_when(
      freq_alal_geno!=1 ~ freq_Aa_geno/(1-freq_alal_geno),
      TRUE ~ 0)) %>%
    dplyr::mutate(updated_freq_aa_geno=case_when(
      freq_alal_geno!=1 ~ freq_aa_geno/(1-freq_alal_geno),
      TRUE ~ 0)) %>%
    dplyr::mutate(updated_freq_Aal_geno=case_when(
      freq_alal_geno!=1 ~ freq_Aal_geno/(1-freq_alal_geno),
      TRUE ~ 0)) %>%
    dplyr::mutate(updated_freq_aal_geno=case_when(
      freq_alal_geno!=1 ~ freq_aal_geno/(1-freq_alal_geno),
      TRUE ~ 0)) %>%
    # Determine Mean Population Fitness
    dplyr::mutate(mean_fit=sum(mating_freq*fam_fit)) %>%
    # Normalize Family Fitness by Mean Pop Fitness
    dplyr::mutate(rel_fam_fit=fam_fit/mean_fit) %>%
    # Extract Freq of Family After Selection
    dplyr::mutate(fam_freq_as=rel_fam_fit*mating_freq)
  
  return(mating_table)
}

# Consider a three allele model
# A = wild-type, a = distorter, al = lethal distorter
# Pull genotype frequencies out from mating table
three_allele_get_genotype_freqs_nextgen <- function(mating_table) {
  geno_freq_next_gen <- mating_table %>% 
    dplyr::summarise(
      AA_ng=sum(updated_freq_AA_geno*fam_freq_as),
      Aa_ng=sum(updated_freq_Aa_geno*fam_freq_as),
      aa_ng=sum(updated_freq_aa_geno*fam_freq_as),
      Aal_ng=sum(updated_freq_Aal_geno*fam_freq_as),
      aal_ng=sum(updated_freq_aal_geno*fam_freq_as),
      alal_ng=sum(updated_freq_alal_geno*fam_freq_as))
  return(geno_freq_next_gen)  
}

# Consider a three allele model
# A = wild-type, a = distorter, al = lethal distorter
# Determine probabilities given genotype frequencies
# for brother's genotype given sister's genotype
# this helps us determine different mating frequencies given inbreeding
three_allele_calculate_inbreeding_genotype_probabilities <- function(mating_table,geno_freq_next_gen) {
  # Get probabilities offspring came from a given mating type given
  # offspring genotype
  mating_given_geno <- mating_table %>%
    dplyr::mutate(pmAA = updated_freq_AA_geno * fam_freq_as / pull(geno_freq_next_gen,AA_ng),
                  pmAa = updated_freq_Aa_geno * fam_freq_as / pull(geno_freq_next_gen,Aa_ng),
                  pmaa = updated_freq_aa_geno * fam_freq_as / pull(geno_freq_next_gen,aa_ng),
                  pmAal = updated_freq_Aal_geno * fam_freq_as / pull(geno_freq_next_gen,Aal_ng),
                  pmaal = updated_freq_aal_geno * fam_freq_as / pull(geno_freq_next_gen,aal_ng),
                  pmalal = updated_freq_alal_geno * fam_freq_as / pull(geno_freq_next_gen,alal_ng)) %>%
    dplyr::mutate(pmAA=case_when(is.na(pmAA) ~ 0, TRUE ~ pmAA)) %>%
    dplyr::mutate(pmAa=case_when(is.na(pmAa) ~ 0, TRUE ~ pmAa)) %>%
    dplyr::mutate(pmaa=case_when(is.na(pmaa) ~ 0, TRUE ~ pmaa)) %>%
    dplyr::mutate(pmAal=case_when(is.na(pmAal) ~ 0, TRUE ~ pmAal)) %>%
    dplyr::mutate(pmaal=case_when(is.na(pmaal) ~ 0, TRUE ~ pmaal)) %>%
    dplyr::mutate(pmalal=case_when(is.na(pmalal) ~ 0, TRUE ~ pmalal)) %>%    
    dplyr::select(female_geno,male_geno,updated_freq_AA_geno,updated_freq_Aa_geno,updated_freq_aa_geno,updated_freq_Aal_geno,updated_freq_aal_geno,updated_freq_alal_geno,pmAA,pmAa,pmaa,pmAal,pmaal,pmalal) %>% 
    dplyr::mutate(MT=paste(female_geno,male_geno,sep='_'))
  
  # Probability of a female offspring's brother having a given genotype
  # given the female offspring's genotype
  prob_brothergeno_given_female_offspring_geno <- mating_given_geno %>% 
    dplyr::select(-female_geno, - male_geno) %>%
    pivot_longer(cols = contains("freq"),values_to = "freq_male_offspring_geno",names_to = "male_offspring_geno")%>%
    separate(male_offspring_geno,sep='_',into=c('char_tag1','char_tag2','male_offspring_geno','char_tag3'))  %>% select(-char_tag1,-char_tag2,-char_tag3) %>%
    pivot_longer(cols = contains("pm"),values_to = "prob_female_offspring_geno_fromMT",names_to = "female_offspring_geno") %>%
    separate(female_offspring_geno,sep='m',into=c('char_tag','female_offspring_geno')) %>% select(-char_tag) %>%
    group_by(female_offspring_geno,male_offspring_geno) %>%
    summarise(prob_brother_geno_given_female_offspring_geno = sum(prob_female_offspring_geno_fromMT *freq_male_offspring_geno))
  
  return(prob_brothergeno_given_female_offspring_geno)
}

# Consider a three allele model
# A = wild-type, a = distorter, al = lethal distorter
# Generate the mating table under random mating
# WITH COMPENSATION (i.e.,the benefit to female fitness by killing off alal homozygotes during embryogenesis)
# Parameters of interest:
#    - k: probability a is transmitted in heterozygous males (drive strength)
#    - s: Fitness reduction in heterozygous Aa males
#    - t: Fitness reduction in homozygous aa males
#    - C: Strength of compensation
three_allele_generate_mating_table_with_compensation_inbreeding <- function(AA_current,Aa_current,aa_current,Aal_current, aal_current, alal_current,k,s,t,f,C,prob_brothergeno_given_female_offspring_geno) {
  
  #Reset Colnames For Left-Join
  colnames(prob_brothergeno_given_female_offspring_geno) <- c('female_geno','male_geno','inbreeding_mating_prob')

  mating_table <- crossing(female_geno=c("AA","Aa","aa","Aal",'aal','alal'),male_geno=c("AA","Aa","aa",'Aal','aal','alal')) %>%
    # Establish Female Genotype Frequencies
    dplyr::mutate(female_geno_freq=case_when(
      female_geno=='AA'~AA_current,
      female_geno=='Aa' ~ Aa_current,
      female_geno=='aa'~aa_current,
      female_geno=='Aal' ~ Aal_current,
      female_geno=='aal' ~ aal_current,
      female_geno=='alal' ~ alal_current)) %>%
    # Establish Male Genotype Frequencies
    dplyr::mutate(male_geno_freq=case_when(
      male_geno=='AA'~AA_current,
      male_geno=='Aa' ~ Aa_current,
      male_geno=='aa'~aa_current,
      male_geno=='Aal' ~ Aal_current,
      male_geno=='aal' ~ aal_current,
      male_geno=='alal' ~ alal_current)) %>%
    # Random Mating Freq
    dplyr::mutate(mating_freq_random=female_geno_freq*male_geno_freq * (1-f)) %>%
    # Inbreed Mating Freq
    left_join(prob_brothergeno_given_female_offspring_geno,by=c('female_geno','male_geno')) %>%
    dplyr::mutate(mating_freq_inbred=inbreeding_mating_prob*f*female_geno_freq) %>%
    # Obtain Final Mating Freq
    dplyr::mutate(mating_freq=mating_freq_random+mating_freq_inbred) %>%
    # Establish Frequencies of A allele in gametes (no drive in females)
    dplyr::mutate(female_A_gam_freq=case_when(
      female_geno=='AA' ~ 1,
      female_geno=='Aa' ~ 0.5,
      female_geno=='aa' ~ 0,
      female_geno=='Aal' ~ 0.5,
      female_geno=='aal' ~ 0,
      female_geno=='alal' ~ 0)) %>%
    # Establish Frequencies of a allele in gametes (no drive in females)
    dplyr::mutate(female_a_gam_freq=case_when(
      female_geno=='AA' ~ 0,
      female_geno=='Aa' ~ 0.5,
      female_geno=='aa' ~ 1,
      female_geno=='Aal' ~ 0,
      female_geno=='aal' ~ 0.5,
      female_geno=='alal' ~ 0)) %>%
    # Establish Frequencies of al allele in gametes (no drive in females)
    dplyr::mutate(female_al_gam_freq=case_when(
      female_geno=='AA' ~ 0,
      female_geno=='Aa' ~ 0,
      female_geno=='aa' ~ 0,
      female_geno=='Aal' ~ 0.5,
      female_geno=='aal' ~ 0.5,
      female_geno=='alal' ~ 1)) %>%
    # Establish Frequencies of A allele in gametes (a allele drives in males)
    dplyr::mutate(male_A_gam_freq=case_when(
      male_geno=='AA' ~ 1,
      male_geno=='Aa' ~ (1-k),
      male_geno=='aa' ~ 0,
      male_geno=='Aal'~ (1-k),
      male_geno=='aal'~ 0,
      male_geno=='alal'~ 0)) %>%
    # Establish Frequencies of a allele in gametes (a allele drives in males)
    dplyr::mutate(male_a_gam_freq=case_when(
      male_geno=='AA' ~ 0,
      male_geno=='Aa' ~ k,
      male_geno=='aa' ~ 1,
      male_geno=='Aal'~ 0,
      male_geno=='aal'~ 0.5,
      male_geno=='alal'~ 0)) %>%
    # Establish Frequencies of al allele in gametes (a allele drives in males)
    dplyr::mutate(male_al_gam_freq=case_when(
      male_geno=='AA' ~ 0,
      male_geno=='Aa' ~ 0,
      male_geno=='aa' ~ 0,
      male_geno=='Aal'~ k,
      male_geno=='aal'~ 0.5,
      male_geno=='alal'~ 1)) %>%
    # From each mating type, establish genotype freqs of offspring
    dplyr::mutate(
      freq_AA_geno=female_A_gam_freq * male_A_gam_freq,
      freq_Aa_geno=female_A_gam_freq*male_a_gam_freq + female_a_gam_freq*male_A_gam_freq,
      freq_aa_geno=female_a_gam_freq*male_a_gam_freq,
      freq_Aal_geno=female_A_gam_freq*male_al_gam_freq + female_al_gam_freq*male_A_gam_freq,
      freq_aal_geno=female_a_gam_freq*male_al_gam_freq + female_al_gam_freq*male_a_gam_freq,
      freq_alal_geno=female_al_gam_freq * male_al_gam_freq) %>%
    # Establish Family Fitness of Mating Type, here dependent on male fitness
    dplyr::mutate(fam_fit=( ((1-freq_alal_geno) / (1-freq_alal_geno)^C))*case_when(
      male_geno=='AA' ~ 1,
      male_geno=='Aa' ~ (1-s),
      male_geno=='aa' ~ (1-t),
      male_geno=='Aal' ~ (1-s),
      male_geno=='aal' ~ (1-t),
      male_geno=='alal' ~ (1-t) )) %>%
    dplyr::mutate(fam_fit=case_when(is.na(fam_fit) ~ 0,TRUE ~ fam_fit)) %>%
    dplyr::mutate(updated_freq_alal_geno=0) %>% 
    dplyr::mutate(updated_freq_AA_geno=case_when(
      freq_alal_geno!=1 ~ freq_AA_geno/(1-freq_alal_geno),
      TRUE ~ 0)) %>%
    dplyr::mutate(updated_freq_Aa_geno=case_when(
      freq_alal_geno!=1 ~ freq_Aa_geno/(1-freq_alal_geno),
      TRUE ~ 0)) %>%
    dplyr::mutate(updated_freq_aa_geno=case_when(
      freq_alal_geno!=1 ~ freq_aa_geno/(1-freq_alal_geno),
      TRUE ~ 0)) %>%
    dplyr::mutate(updated_freq_Aal_geno=case_when(
      freq_alal_geno!=1 ~ freq_Aal_geno/(1-freq_alal_geno),
      TRUE ~ 0)) %>%
    dplyr::mutate(updated_freq_aal_geno=case_when(
      freq_alal_geno!=1 ~ freq_aal_geno/(1-freq_alal_geno),
      TRUE ~ 0)) %>%
    # Determine Mean Population Fitness
    dplyr::mutate(mean_fit=sum(mating_freq*fam_fit)) %>%
    # Normalize Family Fitness by Mean Pop Fitness
    dplyr::mutate(rel_fam_fit=fam_fit/mean_fit) %>%
    # Extract Freq of Family After Selection
    dplyr::mutate(fam_freq_as=rel_fam_fit*mating_freq)
  
  return(mating_table)
}


################################################################################
################################################################################
################################################################################
################################################################################


obtain_nextgen_genotype_freqs_no_inbreeding <- function(AA_current,Aa_current, aa_current,k,s,t,epsilon=0.001) {
  mating_table <- crossing(female_geno=c("AA","Aa","aa"),male_geno=c("AA","Aa","aa")) %>%
    # Establish Female Genotype Frequencies
    dplyr::mutate(female_geno_freq=case_when(
      female_geno=='AA'~AA_current,
      female_geno=='Aa' ~ Aa_current,
      female_geno=='aa'~aa_current)) %>%
    # Establish Male Genotype Frequencies
    dplyr::mutate(male_geno_freq=case_when(
      male_geno=='AA'~ AA_current,
      male_geno=='Aa' ~ Aa_current,
      male_geno=='aa'~ aa_current)) %>%
    # Establish Mating Type Frequencies
    dplyr::mutate(mating_freq=female_geno_freq*male_geno_freq) %>%
    # Establish Frequencies of A allele in gametes (no drive in females)
    dplyr::mutate(female_A_gam_freq=case_when(
      female_geno=='AA' ~ 1,
      female_geno=='Aa' ~ 0.5,
      female_geno=='aa' ~ 0)) %>%
    # Establish Frequencies of a allele in gametes (no drive in females)
    dplyr::mutate(female_a_gam_freq=case_when(
      female_geno=='AA' ~ 0,
      female_geno=='Aa' ~ 0.5,
      female_geno=='aa' ~ 1)) %>%
    # Establish Frequencies of A allele in gametes (a allele drives in males)
    dplyr::mutate(male_A_gam_freq=case_when(
      male_geno=='AA' ~ 1,
      male_geno=='Aa' ~ (1-k),
      male_geno=='aa' ~ 0)) %>%
    # Establish Frequencies of a allele in gametes (a allele drives in males)
    dplyr::mutate(male_a_gam_freq=case_when(
      male_geno=='AA' ~ 0,
      male_geno=='Aa' ~ k,
      male_geno=='aa' ~ 1)) %>%
    # From each mating type, establish genotype freqs of offspring
    dplyr::mutate(
      freq_AA_geno=female_A_gam_freq * male_A_gam_freq,
      freq_Aa_geno=female_A_gam_freq*male_a_gam_freq + female_a_gam_freq*male_A_gam_freq,
      freq_aa_geno=female_a_gam_freq*male_a_gam_freq) %>%
    # Establish Family Fitness of Mating Type, here dependent on male fitness
    dplyr::mutate(fam_fit=case_when(
      male_geno=='AA' ~ 1,
      male_geno=='Aa' ~ (1-s),
      male_geno=='aa' ~ (1-t))) %>%
    # Determine Mean Population Fitness
    dplyr::mutate(mean_fit=sum(mating_freq*fam_fit)) %>%
    # Normalize Family Fitness by Mean Pop Fitness
    dplyr::mutate(rel_fam_fit=fam_fit/mean_fit) %>%
    # Extract Freq of Family After Selection
    dplyr::mutate(fam_freq_as=rel_fam_fit*mating_freq)
  
  geno_freq_next_gen <- mating_table %>% 
    dplyr::summarise(
      AA_ng=sum(freq_AA_geno*fam_freq_as),
      Aa_ng=sum(freq_Aa_geno*fam_freq_as),
      aa_ng=sum(freq_aa_geno*fam_freq_as))
  return(as.numeric(geno_freq_next_gen[1,]))
}

obtain_nextgen_genotype_freqs_with_inbreeding <- function(AA_current,Aa_current, aa_current,k,s,t,epsilon=0.001) {
  mating_table <- crossing(female_geno=c("AA","Aa","aa"),male_geno=c("AA","Aa","aa")) %>%
    # Establish Female Genotype Frequencies
    dplyr::mutate(female_geno_freq=case_when(
      female_geno=='AA'~AA_current,
      female_geno=='Aa' ~ Aa_current,
      female_geno=='aa'~aa_current)) %>%
    # Establish Male Genotype Frequencies
    dplyr::mutate(male_geno_freq=case_when(
      male_geno=='AA'~ AA_current,
      male_geno=='Aa' ~ Aa_current,
      male_geno=='aa'~ aa_current)) %>%
    # Establish Mating Type Frequencies
    dplyr::mutate(mating_freq=female_geno_freq*male_geno_freq) %>%
    # Establish Frequencies of A allele in gametes (no drive in females)
    dplyr::mutate(female_A_gam_freq=case_when(
      female_geno=='AA' ~ 1,
      female_geno=='Aa' ~ 0.5,
      female_geno=='aa' ~ 0)) %>%
    # Establish Frequencies of a allele in gametes (no drive in females)
    dplyr::mutate(female_a_gam_freq=case_when(
      female_geno=='AA' ~ 0,
      female_geno=='Aa' ~ 0.5,
      female_geno=='aa' ~ 1)) %>%
    # Establish Frequencies of A allele in gametes (a allele drives in males)
    dplyr::mutate(male_A_gam_freq=case_when(
      male_geno=='AA' ~ 1,
      male_geno=='Aa' ~ (1-k),
      male_geno=='aa' ~ 0)) %>%
    # Establish Frequencies of a allele in gametes (a allele drives in males)
    dplyr::mutate(male_a_gam_freq=case_when(
      male_geno=='AA' ~ 0,
      male_geno=='Aa' ~ k,
      male_geno=='aa' ~ 1)) %>%
    # From each mating type, establish genotype freqs of offspring
    dplyr::mutate(
      freq_AA_geno=female_A_gam_freq * male_A_gam_freq,
      freq_Aa_geno=female_A_gam_freq*male_a_gam_freq + female_a_gam_freq*male_A_gam_freq,
      freq_aa_geno=female_a_gam_freq*male_a_gam_freq) %>%
    # Establish Family Fitness of Mating Type, here dependent on male fitness
    dplyr::mutate(fam_fit=case_when(
      male_geno=='AA' ~ 1,
      male_geno=='Aa' ~ (1-s),
      male_geno=='aa' ~ (1-t))) %>%
    # Determine Mean Population Fitness
    dplyr::mutate(mean_fit=sum(mating_freq*fam_fit)) %>%
    # Normalize Family Fitness by Mean Pop Fitness
    dplyr::mutate(rel_fam_fit=fam_fit/mean_fit) %>%
    # Extract Freq of Family After Selection
    dplyr::mutate(fam_freq_as=rel_fam_fit*mating_freq)
  
  geno_freq_next_gen <- mating_table %>% 
    dplyr::summarise(
      AA_ng=sum(freq_AA_geno*fam_freq_as),
      Aa_ng=sum(freq_Aa_geno*fam_freq_as),
      aa_ng=sum(freq_aa_geno*fam_freq_as))
  
  # for  inbreeding stuff in the next gen
  mating_given_geno <- mating_table %>%
    dplyr::mutate(pmAA = freq_AA_geno * fam_freq_as / pull(geno_freq_next_gen,AA_ng),
                  pmAa = freq_Aa_geno * fam_freq_as / pull(geno_freq_next_gen,Aa_ng),
                  pmaa = freq_aa_geno * fam_freq_as / pull(geno_freq_next_gen,aa_ng)) %>%
    dplyr::select(female_geno,male_geno,freq_AA_geno,freq_Aa_geno,freq_aa_geno,pmAA,pmAa,pmaa) %>% 
    dplyr::mutate(MT=paste(female_geno,male_geno,sep='_'))
  
  # Probability of a female offspring's brother having a given genotype
  # given the female offspring's genotype
  prob_brothergeno_given_female_offspring_geno <- mating_given_geno %>% 
    dplyr::select(-female_geno, - male_geno) %>%
    pivot_longer(cols = contains("freq"),values_to = "freq_male_offspring_geno",names_to = "male_offspring_geno")%>%
    separate(male_offspring_geno,sep='_',into=c('char_tag','male_offspring_geno','char_tag2'))  %>% select(-char_tag,-char_tag2) %>%
    pivot_longer(cols = contains("pm"),values_to = "prob_female_offspring_geno_fromMT",names_to = "female_offspring_geno") %>%
    separate(female_offspring_geno,sep='m',into=c('char_tag','female_offspring_geno')) %>% select(-char_tag) %>%
    group_by(female_offspring_geno,male_offspring_geno) %>%
    summarise(prob_brother_geno_given_female_offspring_geno = sum(prob_female_offspring_geno_fromMT *freq_male_offspring_geno))
  
  return(as.numeric(geno_freq_next_gen[1,]))
}

