## Script Used to Store Functions Required to
# Generate Figures in the Manuscript

# Load and read in all of the result files
extract_results <- function(dir_name) {
	file_names <- list.files(dir_name,full.names=T)
	all_data <- lapply(file_names,fread)
	all_data <- do.call('rbind',all_data)
	return(all_data)
}

# From result files, extract 
# final alelle frequency
obtain_faf <- function(results,model_type) {
	if (model_type=='Model1') {
		faf <- results %>% group_by(s,t,k) %>%
		filter(gen==max(gen)) %>% 
		dplyr::mutate(a_freq=(0.5*Aa_freq)+aa_freq)
	}

	if (model_type=='DL1961') {
		faf <- results %>% group_by(s,t,k) %>%
		filter(gen==max(gen)) %>% 
		dplyr::mutate(a_freq=(0.5*Aa_freq)+aa_freq)%>%
		dplyr::mutate(al_freq=(0.5*Aal_freq)+(0.5*aal_freq)+alal_freq)
	}	

	if (model_type=='Model2') {
		faf <- results %>% group_by(s,t,k,f) %>%
		filter(gen==max(gen)) %>% 
		dplyr::mutate(a_freq=(0.5*Aa_freq)+aa_freq)
	}

	if (model_type=='P1967') {
		faf <- results %>% group_by(s,t,k,f,C) %>%
		filter(gen==max(gen)) %>% 
		dplyr::mutate(a_freq=(0.5*Aa_freq)+aa_freq)%>%
		dplyr::mutate(al_freq=(0.5*Aal_freq)+(0.5*aal_freq)+alal_freq)
	}
	
	if (model_type=='Model3') {
		faf <- results %>% group_by(s,t,k,C) %>%
		filter(gen==max(gen)) %>%
		dplyr::mutate(a_freq=(0.5*Aa_freq)+(0.5*aal_freq)+aa_freq) %>%
		dplyr::mutate(al_freq=(0.5*Aal_freq)+(0.5*aal_freq)+alal_freq) %>%
		dplyr::mutate(updated_al_freq=case_when(
			introduce_al==FALSE ~ NA_real_,
			TRUE ~ al_freq
			)) %>%
		dplyr::mutate(t_freq=a_freq+updated_al_freq)
	}

	if (model_type=='Model4') {
		faf <- results %>% group_by(s,t,k,f,C) %>%
		filter(gen==max(gen)) %>%
		dplyr::mutate(a_freq=(0.5*Aa_freq)+(0.5*aal_freq)+aa_freq) %>%
		dplyr::mutate(al_freq=(0.5*Aal_freq)+(0.5*aal_freq)+alal_freq) %>%
		dplyr::mutate(updated_al_freq=case_when(
			introduce_al==FALSE ~ NA_real_,
			TRUE ~ al_freq
			)) %>%
		dplyr::mutate(t_freq=a_freq+updated_al_freq)
	}

	return(faf)
}

# Get expected frequencies for Petras 1967
# from transformed prop sib-mating
petras_prediction <- function(m,F,type) {
	part1 <- ((4*m) - F - (6*F*m))
	part2 <- sqrt(((4*F*m) * (F - 8 + (F*m) + (10*m))) + ((16*m)*(1-m)) + (F*F))
	part3 <- (8*m)*(1-F)

	if (type=='add') {
		output_val <- (part1 + part2)/part3
	} else if (type=='minus') {
		output_val <- (part1 - part2)/part3
	}
	return(output_val)
}