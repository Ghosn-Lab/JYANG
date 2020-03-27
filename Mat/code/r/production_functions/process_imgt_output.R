#' process_imgt_output
#'
#' @param imgt_summary_file_directory Quoted file path of the IMGT summary file output folder
#' @param processed_cellranger_annotations_file Quoted file path of the CellRanger filtered and processed annotations .csv file from process_cellranger_output
#' @param imgt_individual_file_directory Quoted file path of the IMGT individual file output folder
#' @param return_df Returns the df to the console if TRUE else just writes the output .csv file
#'
#' @return
#' @export
#'
#' @examples

process_imgt_output <- function(imgt_summary_file_directory, processed_cellranger_annotations_file, imgt_individual_file_directory = FALSE, return_df = FALSE) {
  
  ###Check for dependencies
  require(plyr)
  require(dplyr)
  require(stringr)
  require(readr)
  require(Peptides)
  
###File Processing Functions
  
  ########################
  ##Process IMGT Summary File
  ########################
  
  process_imgt_summary_file <- function(summary_directory = imgt_summary_file_directory) {
    
    ##Set the working directory
    setwd(summary_directory)
    
    ##Read the summary file
    summary_file <- read_tsv("1_Summary.txt")
    
    ##Rename the summary file columns
    names(summary_file) <- tolower(names(summary_file))
    names(summary_file) <- gsub(" ", "_", names(summary_file))
    names(summary_file) <- gsub("-", "_", names(summary_file))
    names(summary_file) <- gsub("%", "pct", names(summary_file))
    
    ##Re-organize the V and J allele information
    summary_file$multiple_v_alleles_possible <- grepl("\\,", summary_file$v_gene_and_allele)
    summary_file$multiple_j_alleles_possible <- grepl("\\,", summary_file$j_gene_and_allele)
    summary_file$primary_v_allele <- tolower(gsub(".*(IG.*?) .*", "\\1", summary_file$v_gene_and_allele))
    summary_file$v_gene <- gsub("(.*)\\*.*", "\\1", summary_file$primary_v_allele)
    summary_file$primary_j_allele <- tolower(gsub(".*(IG.*?) .*", "\\1", summary_file$j_gene_and_allele))
    summary_file$j_gene <- gsub("(.*)\\*.*", "\\1", summary_file$primary_j_allele)
    summary_file$d_gene <- tolower(gsub(".*(IG.*)\\*.*", "\\1", summary_file$d_gene_and_allele))
    
    ##Process the summary file down to needed columns
    processed_summary_file <- summary_file %>%
      select(sequence_id,
             analysed_sequence_length,
             orientation,
             v_domain_functionality,
             v_gene,
             primary_v_allele,
             multiple_v_alleles_possible,
             v_region_score,
             v_region_identity_nt,
             v_region_identity_pct,
             j_gene,
             primary_j_allele,
             multiple_j_alleles_possible,
             j_region_score,
             j_region_identity_nt,
             j_region_identity_pct) %>%
      mutate(chain_type = gsub("(ig.).*", "\\1", v_gene))
    
    return(processed_summary_file)
  
    } ### End, process_summary_file()
  
  
  #########################
  ##Process IMGT Gapped nt Sequences file
  #########################
  
  process_imgt_gapped_nt_sequences_file <- function(summary_directory = imgt_summary_file_directory) {
    
    ##Set the working directory
    setwd(summary_directory)
    
    ##Read the IMGT Gapped nt Sequences file
    imgt_gapped_nt_sequences <- read_tsv("2_IMGT-gapped-nt-sequences.txt")
    
    ##Rename the IMGT Gapped nt Sequences file columns
    names(imgt_gapped_nt_sequences) <- tolower(names(imgt_gapped_nt_sequences))
    names(imgt_gapped_nt_sequences) <- gsub(" ", "_", names(imgt_gapped_nt_sequences))
    names(imgt_gapped_nt_sequences) <- gsub("-", "_", names(imgt_gapped_nt_sequences))
    
    ##Process the IMGT Gapped nt Sequences file columns
    processed_imgt_gapped_nt_sequences <- imgt_gapped_nt_sequences %>%
      select(sequence_id,
             nt_gapped_vdj = v_d_j_region,
             nt_gapped_vj = v_j_region,
             nt_gapped_v = v_region,
             nt_gapped_fr1 = fr1_imgt,
             nt_gapped_fr2 = fr2_imgt,
             nt_gapped_fr3 = fr3_imgt,
             nt_gapped_fr4 = fr4_imgt,
             nt_gapped_cdr1 = cdr1_imgt,
             nt_gapped_cdr2 = cdr2_imgt,
             nt_gapped_cdr3 = cdr3_imgt,
             nt_gapped_junction = junction,
             nt_gapped_j = j_region)
    
    return(processed_imgt_gapped_nt_sequences)  
  
  } ### End, process_imgt_gapped_nt_sequences_file()
  
  
  
  ########################
  ##Process IMGT NT Sequences file
  ########################
  
  process_imgt_nt_sequences_file <- function(summary_directory = imgt_summary_file_directory) {
    
    ##Set the working directory
    setwd(summary_directory)
    
    ##Read the nt Sequences file
    nt_sequences <- read_tsv("3_Nt-sequences.txt")
    
    ##Rename the nt Sequences file columns
    names(nt_sequences) <- tolower(names(nt_sequences))
    names(nt_sequences) <- gsub(" ", "_", names(nt_sequences))
    names(nt_sequences) <- gsub("-", "_", names(nt_sequences))
    names(nt_sequences) <- gsub("%", "pct", names(nt_sequences))
    
    ##Process the nt Sequences file
    processed_nt_sequences <- nt_sequences %>%
      select(sequence_id,
             nt_vdj = v_d_j_region,
             nt_vj = v_j_region,
             nt_v = v_region,
             nt_fr1 = fr1_imgt,
             nt_fr2 = fr2_imgt,
             nt_fr3 = fr3_imgt,
             nt_fr4 = fr4_imgt,
             nt_cdr1 = cdr1_imgt,
             nt_cdr2 = cdr2_imgt,
             nt_cdr3 = cdr3_imgt,
             nt_junction = junction,
             nt_j = j_region) %>%
      mutate(nt_cdr3_length = nchar(nt_cdr3))
    
    return(processed_nt_sequences)

  }###End, process_imgt_nt_sequences_file
  
  
  #########################
  ##Process the IMGT gapped AA sequences file
  #########################
  
  
  process_imgt_gapped_aa_sequences_file <- function(summary_directory = imgt_summary_file_directory) {
    
    ##Set the working directory
    setwd(summary_directory)
    
    ##Read the imgt_gapped_aa_sequences file
    imgt_gapped_aa_sequences <- read_tsv("4_IMGT-gapped-AA-sequences.txt")
    
    ##Rename the imgt_gapped_aa_sequences file columns
    names(imgt_gapped_aa_sequences) <- tolower(names(imgt_gapped_aa_sequences))
    names(imgt_gapped_aa_sequences) <- gsub(" ", "_", names(imgt_gapped_aa_sequences))
    names(imgt_gapped_aa_sequences) <- gsub("-", "_", names(imgt_gapped_aa_sequences))
    
    ##Process the imgt_gapped_aa_sequences file columns
    processed_imgt_gapped_aa_sequences <- imgt_gapped_aa_sequences %>%
      select(sequence_id,
             aa_gapped_vdj = v_d_j_region,
             aa_gapped_vj = v_j_region,
             aa_gapped_v = v_region,
             aa_gapped_fr1 = fr1_imgt,
             aa_gapped_fr2 = fr2_imgt,
             aa_gapped_fr3 = fr3_imgt,
             aa_gapped_fr4 = fr4_imgt,
             aa_gapped_cdr1 = cdr1_imgt,
             aa_gapped_cdr2 = cdr2_imgt,
             aa_gapped_cdr3 = cdr3_imgt,
             aa_gapped_junction = junction,
             aa_gapped_j = j_region)
    
    return(processed_imgt_gapped_aa_sequences)
    
  } ###End, process_imgt_gapped_aa_sequences_file()
  
  
  ########################
  ##Process the IMGT AA sequences file
  ########################
  
  
  process_imgt_aa_sequences_file <- function(summary_directory = imgt_summary_file_directory){
    
    ##Set the working directory
    setwd(summary_directory)
    
    ##Read the aa_sequences file
    aa_sequences <- read_tsv("5_AA-sequences.txt")
    
    ##Rename the aa_sequences file
    names(aa_sequences) <- tolower(names(aa_sequences))
    names(aa_sequences) <- gsub(" ", "_", names(aa_sequences))
    names(aa_sequences) <- gsub("-", "_", names(aa_sequences))
    names(aa_sequences) <- gsub("%", "pct", names(aa_sequences))
    
    ##Process the aa_sequences file
    processed_aa_sequences <- aa_sequences %>%
      select(sequence_id,
             aa_vdj = v_d_j_region,
             aa_vj = v_j_region,
             aa_v = v_region,
             aa_fr1 = fr1_imgt,
             aa_fr2 = fr2_imgt,
             aa_fr3 = fr3_imgt,
             aa_fr4 = fr4_imgt,
             aa_cdr1 = cdr1_imgt,
             aa_cdr2 = cdr2_imgt,
             aa_cdr3 = cdr3_imgt,
             aa_junction = junction,
             aa_j = j_region) 
    
    ##Establish a string of unique CDR3s
    unique_cdr3s <- unique(processed_aa_sequences$aa_cdr3)
    
    ##Make that string into a df
    cdr3_property_tibble <- tibble(aa_cdr3 = unique_cdr3s)
    
    ##Calculate the CDR3 peptide length
    cdr3_property_tibble$cdr3_aa_length <- lengthpep(unique_cdr3s)
    
    ##Calculate the CDR3 MW
    print("Calculating CDR3 Molecular Weights")
    cdr3_property_tibble$cdr3_molecular_weight <- mw(unique_cdr3s)
    
    ##Calculate the CDR3 charge
    print("Calculating CDR3 Charges")
    cdr3_property_tibble$cdr3_net_charge <- charge(unique_cdr3s)
    
    ##Calculate the CDR3 isoelectric point
    print("Calculating CDR3 Isoelectric Points")
    cdr3_property_tibble$cdr3_isoelectric_point <- pI(unique_cdr3s)
    
    ##Calculate the CDR3 hydrophobicity
    print("Calculating CDR3 Hydrophobicity")
    cdr3_property_tibble$cdr3_hydrophobicity <- hydrophobicity(unique_cdr3s)
    
    ##Calculate the CDR3 aliphatic index
    print("Calculating CDR3 Aliphatic Index")
    cdr3_property_tibble$cdr3_aliphatic_index <- aIndex(unique_cdr3s)
    
    ##Calculate the CDR3 MW
    print("Calculating CDR3 boman index")
    cdr3_property_tibble$cdr3_boman_index <- boman(unique_cdr3s)
    
    ##Calculate the CDR3 boman prediction
    cdr3_property_tibble$cdr3_boman_predicts_protein_interaction <- cdr3_property_tibble$cdr3_boman_index >= 2.48
    
    ##Calculate kidera factors
    print("Calculating Kidera Factors")
    kidera_factors <- as_tibble(t(bind_cols(kideraFactors(unique_cdr3s))))
    names(kidera_factors) <- c("cdr3_kf_helix_bend_preference",
                               "cdr3_kf_side_chain_size",
                               "cdr3_kf_extended_structure_preference",
                               "cdr3_kf_hydrophobicity",
                               "cdr3_kf_double_bend_preference",
                               "cdr3_kf_partial_specific_volume",
                               "cdr3_kf_flat_extended_preference",
                               "cdr3_kf_occurance_in_alpha_region",
                               "cdr3_kf_pK_c",
                               "cdr3_kf_surrounding_hydrophobicity")
    
    ##Join standard calculations and kidera factor calculations
    cdr3_property_tibble <- bind_cols(cdr3_property_tibble, kidera_factors)
    
    ##Join processed aa sequences and CDR3 properties
    processed_aa_sequences <- left_join(processed_aa_sequences, cdr3_property_tibble, by = "aa_cdr3")
    
    print("Done")
    
    return(processed_aa_sequences)
    
  }## End process_imgt_aa_sequences_file()
  
  
  #########################
  ##Process IMGT V-REGION-nt-mutation-statistics file
  #########################
  
  process_imgt_v_region_nt_mutation_statistics_file <- function(summary_directory = imgt_summary_file_directory) {
    
    ##Set the working directory
    setwd(summary_directory)
    
    ##Read the V-REGION-nt-mutation-statistics file
    v_region_nt_mut_stats <- read_tsv("8_V-REGION-nt-mutation-statistics.txt")
    
    ##Rename the V-REGION-nt-mutation-statistics file columns
    names(v_region_nt_mut_stats) <- tolower(names(v_region_nt_mut_stats))
    names(v_region_nt_mut_stats) <- gsub(" ", "_", names(v_region_nt_mut_stats))
    names(v_region_nt_mut_stats) <- gsub("-", "_", names(v_region_nt_mut_stats))
    names(v_region_nt_mut_stats) <- gsub("%", "pct", names(v_region_nt_mut_stats))
    names(v_region_nt_mut_stats) <- gsub(">", "_into_", names(v_region_nt_mut_stats))
    
    ## Select necessary columns
    processed_v_region_nt_mut_stats <- v_region_nt_mut_stats %>%
      select(2, 5:130)
    
    ## Clean the mutation data of IMGT gaps
    #Establish a progress bar
    print("Cleaning V region NT change statistics")
    pb <- txtProgressBar(min = 0, max = length(processed_v_region_nt_mut_stats), style = 3)
    
    #Clean the data of IMGT gaps
    for (i in 2:length(processed_v_region_nt_mut_stats)) {
      if (class(processed_v_region_nt_mut_stats[[i]]) == "character") {
        for (j in 1:nrow(processed_v_region_nt_mut_stats)) {
          if (grepl(" \\(.*\\)", processed_v_region_nt_mut_stats[j,i])) {
            processed_v_region_nt_mut_stats[j,i] <- gsub(" \\(.*\\)", "", processed_v_region_nt_mut_stats[j,i])
          }
        }
      } else {
        next()
      }
      setTxtProgressBar(pb, i)
    }
    
    #Convert cleaned columns to numeric values
    for (i in 2:length(processed_v_region_nt_mut_stats)) {
      if (class(processed_v_region_nt_mut_stats[[i]]) == "character") {
      processed_v_region_nt_mut_stats[[i]] <- as.numeric(processed_v_region_nt_mut_stats[[i]])
      }
    }
    
    ##Process v_region_nt_mut_stats file
    
    return(processed_v_region_nt_mut_stats)
    
  } ##End, process_imgt_v_region_nt_mutation_statistics_file()
  
  
  
  ######################
  ##Process V-REGION-AA-change-statistics file
  ######################
  
  process_imgt_v_region_aa_change_statistics_file <- function(summary_directory = imgt_summary_file_directory) {
    
    ##Set the working directory
    setwd(summary_directory)
    
    ##Read the V-REGION-AA-change-statistics file
    v_region_change_stats <- read_tsv("9_V-REGION-AA-change-statistics.txt")
    
    ##Read the V-REGION-AA-change-statistics columns
    names(v_region_change_stats) <- tolower(names(v_region_change_stats))
    names(v_region_change_stats) <- gsub(" ", "_", names(v_region_change_stats))
    names(v_region_change_stats) <- gsub("-", "_", names(v_region_change_stats))
    names(v_region_change_stats) <- gsub("___$", "---", names(v_region_change_stats))
    names(v_region_change_stats) <- gsub("__\\+$", "--+", names(v_region_change_stats))
    names(v_region_change_stats) <- gsub("_\\+_$", "-+-", names(v_region_change_stats))
    names(v_region_change_stats) <- gsub("\\+__$", "+--", names(v_region_change_stats))
    names(v_region_change_stats) <- gsub("_\\+\\+$", "-++", names(v_region_change_stats))
    names(v_region_change_stats) <- gsub("\\+\\+_$", "++-", names(v_region_change_stats))
    names(v_region_change_stats) <- gsub("\\+_\\+$", "+-+", names(v_region_change_stats))
    
    ##Select the V-REGION-AA-change-statistics columns
    processed_v_region_change_stats <- v_region_change_stats %>%
      select(2, 5:109)
    
    ##Clean the AA change statistics of inserted IMGT gap data
    #Establish a progress bar
    print("Cleaning V region AA change statistics")
    pb <- txtProgressBar(min = 0, max = length(processed_v_region_change_stats), style = 3)
    
    #Clean the data of IMGT gaps
    for (i in 2:length(processed_v_region_change_stats)) {
      if (class(processed_v_region_change_stats[[i]]) == "character") {
        for (j in 1:nrow(processed_v_region_change_stats)) {
          if (grepl(" \\(.*\\)", processed_v_region_change_stats[j,i])) {
            processed_v_region_change_stats[j,i] <- gsub(" \\(.*\\)", "", processed_v_region_change_stats[j,i])
          }
        }
      } else {
        next()
      }
      setTxtProgressBar(pb, i)
    }
    
    #Set cleaned data to numeric column classes
    for (i in 2:length(processed_v_region_change_stats)) {
      if (class(processed_v_region_change_stats[[i]]) == "character") {
      processed_v_region_change_stats[[i]] <- as.numeric(processed_v_region_change_stats[[i]])
      }
    }
    
    return(processed_v_region_change_stats)

  }##End, process_imgt_v_region_aa_change_statistics_file()
    
  
  ###########################
  ##Process IMGT V-REGION-mutation-hotspots files
  ###########################
  
  
  process_imgt_v_region_mutation_hotspots_files <- function(summary_directory = imgt_summary_file_directory, individual_directory = imgt_individual_file_directory){
    
    ##Set the working directory
    setwd(summary_directory)
    
    ##Read IMGT V-REGION-mutation-hotspots summary file
     v_region_hotspots <- read_tsv("10_V-REGION-mutation-hotspots.txt")
    
    ##Rename V-REGION-mutation-hotspots summary file columns
    names(v_region_hotspots) <- tolower(names(v_region_hotspots))
    names(v_region_hotspots) <- gsub(" ", "_", names(v_region_hotspots))
    names(v_region_hotspots) <- gsub("-", "_", names(v_region_hotspots))
    
    ##Select needed columns
    processed_v_region_hotspots <- v_region_hotspots %>%
      select(2, 5:8)
    
    ##Initialize hotspot mutation vactor columns for each region
    fr1_hotspot_mut_vector <- vector(length = nrow(v_region_hotspots))
    fr2_hotspot_mut_vector <- vector(length = nrow(v_region_hotspots))
    fr3_hotspot_mut_vector <- vector(length = nrow(v_region_hotspots))
    cdr1_hotspot_mut_vector <- vector(length = nrow(v_region_hotspots))
    cdr2_hotspot_mut_vector <- vector(length = nrow(v_region_hotspots))
    cdr3_hotspot_mut_vector <- vector(length = nrow(v_region_hotspots))
    
    ##Initialize a progress bar
    print("Identifying hotspot loci")
    pb <- txtProgressBar(min = 0, max = nrow(v_region_hotspots), style = 3)
    
    #Identify hotspot locations
    for(i in 1:nrow(v_region_hotspots)) {
      
      flattened_string <- str_flatten(v_region_hotspots[i,])
      
      fr1_hotspot_mut_vector[i] <- str_count(flattened_string, "FR1")
      fr2_hotspot_mut_vector[i] <- str_count(flattened_string, "FR2")
      fr3_hotspot_mut_vector[i] <- str_count(flattened_string, "FR3")
      cdr1_hotspot_mut_vector[i] <- str_count(flattened_string, "CDR1")
      cdr2_hotspot_mut_vector[i] <- str_count(flattened_string, "CDR2")
      cdr3_hotspot_mut_vector[i] <- str_count(flattened_string, "CDR3")
      
      setTxtProgressBar(pb, i)
    
    }
    
    #Condense hotspot locations into a data frame
    processed_v_region_hotspots <- tibble(sequence_id = v_region_hotspots$sequence_id,
                             fr1_hotspot_loci = fr1_hotspot_mut_vector,
                             fr2_hotspot_loci = fr2_hotspot_mut_vector,
                             fr3_hotspot_loci = fr3_hotspot_mut_vector,
                             cdr1_hotspot_loci = cdr1_hotspot_mut_vector,
                             cdr2_hotspot_loci = cdr2_hotspot_mut_vector,
                             cdr3_hotspot_loci = cdr3_hotspot_mut_vector,
                             total_non_junction_hotspot_loci = sum(fr1_hotspot_loci,
                                                  fr2_hotspot_loci,
                                                  fr3_hotspot_loci,
                                                  cdr1_hotspot_loci,
                                                  cdr2_hotspot_loci))
    
    ##Set working directory for individual files
    setwd(individual_directory)
    
    ##Get the individual file list
    file_list <- list.files(individual_directory)
    
    ##Initialize a list with an element for each sequence
    sequence_list <- vector("list", length = length(file_list))
    
    ##Initialize a status bar
    print("Identifying hotspot mutations")
    pb <- txtProgressBar(min = 0, max = length(file_list), style = 3)
    
    ##Find the hotspot mutations for each individual sequence
    for (i in 1:length(file_list)) {
      
      #Read in the individual file
      seq <- read_file(file_list[i])
      #Pull the hotspot text chunk
      hotspot_chunk <- gsub(".*and hotspots motifs\\\n(.*)\\\n\\\n\\\n12.*.", "\\1", seq)
      #Break the chunk into individual mutations
      mutation_string <- str_split(hotspot_chunk, "\n")[[1]]
      # #Get rid of the aggregate V mutations
      # mutation_string_vs_removed <- mutation_string[!grepl("V-REGION", mutation_string)]
      #Identify mutations in hotspots
      hotspot_mutations <- mutation_string[grepl("\\[", mutation_string)]
      #Flatten the hotspot mutations
      combined_hotspot_mutations <- str_flatten(hotspot_mutations)
      
      #If there are no hotspot mutations, fill with zeros
      if(length(combined_hotspot_mutations) == 0) {
        
      fr1_hotspot_muts <- 0
      fr2_hotspot_muts <- 0
      fr3_hotspot_muts <- 0
      cdr1_hotspot_muts <- 0
      cdr2_hotspot_muts <- 0
      cdr3_hotspot_muts <- 0
      total_non_junction_hotspot_muts <- 0
      
      #Else, count the hotspot mutation regions 
      } else {
      
      fr1_hotspot_muts <- str_count(combined_hotspot_mutations, "FR1")
      fr2_hotspot_muts <- str_count(combined_hotspot_mutations, "FR2")
      fr3_hotspot_muts <- str_count(combined_hotspot_mutations, "FR3")
      cdr1_hotspot_muts <- str_count(combined_hotspot_mutations, "CDR1")
      cdr2_hotspot_muts <- str_count(combined_hotspot_mutations, "CDR2")
      cdr3_hotspot_muts <- str_count(combined_hotspot_mutations, "CDR3")
      total_non_junction_hotspot_muts <- sum(fr1_hotspot_muts, fr2_hotspot_muts, fr3_hotspot_muts, cdr1_hotspot_muts, cdr2_hotspot_muts)
      
      }
      
      #Make a tibble with all of the hotspot mutations by region 
      hotspot_df <- tibble(sequence_id = gsub("(.*)_[0-9]+$", "\\1", file_list[i]),
                           fr1_hotspot_muts, fr2_hotspot_muts, fr3_hotspot_muts, cdr1_hotspot_muts, cdr2_hotspot_muts, cdr3_hotspot_muts, total_non_junction_hotspot_muts)
      
      #Add that tibble to the sequence list
      sequence_list[[i]] <- hotspot_df
      
      setTxtProgressBar(pb, i)
      
    }##End, FOR loop
    
    #Bind all of the sequences back together
    hotspot_df <- bind_rows(sequence_list)
    
    #Join hotspot locations and mutation data
    processed_v_region_hotspots <- left_join(processed_v_region_hotspots, hotspot_df, by = "sequence_id")
    
    return(processed_v_region_hotspots)
    
  }##End, process_imgt_v_region_mutation_hotspots_files()
  
  ############################
  ###Global processing
  ############################
  
  ##Create all of the processed files using the above functions
  processed_summary <- process_imgt_summary_file()
  processed_nt_sequences <- process_imgt_nt_sequences_file()
  processed_gapped_nt_sequences <- process_imgt_gapped_nt_sequences_file()
  processed_aa_sequences <- process_imgt_aa_sequences_file()
  processed_gapped_aa_sequences <- process_imgt_gapped_aa_sequences_file()
  processed_nt_mutation_summary <- process_imgt_v_region_nt_mutation_statistics_file()
  processed_aa_mutation_summary <- process_imgt_v_region_aa_change_statistics_file()
  
  ##Process hotspots if individual files are present
  processed_hotspot_mutation_summary <- FALSE
  if (imgt_individual_file_directory != FALSE) {
  processed_hotspot_mutation_summary <- process_imgt_v_region_mutation_hotspots_files()
  }
  
  ##Bind them together using the sequence_id 
  combined_tibble <- plyr::join_all(list(processed_summary,
                                         processed_nt_sequences,
                                         processed_gapped_nt_sequences,
                                         processed_aa_sequences,
                                         processed_gapped_aa_sequences,
                                         processed_nt_mutation_summary,
                                         processed_aa_mutation_summary),
                                    by = "sequence_id")
  
  if (processed_hotspot_mutation_summary != FALSE) {
    combined_tibble <- left_join(combined_tibble, processed_hotspot_mutation_summary, by = "sequence_id")
  }
  
  ##Add additional statistics of relavance to future study
  combined_tibble <- combined_tibble %>%
  mutate(fr1_nt_mut_freq = fr1_imgt_nb_of_mutations / fr1_imgt_nb_of_nucleotides * 100, ## Mutation frequencies
         fr2_nt_mut_freq = fr2_imgt_nb_of_mutations / fr2_imgt_nb_of_nucleotides * 100,
         fr3_nt_mut_freq = fr3_imgt_nb_of_mutations / fr3_imgt_nb_of_nucleotides * 100,
         cdr1_nt_mut_freq = cdr1_imgt_nb_of_mutations / cdr1_imgt_nb_of_nucleotides * 100,
         cdr2_nt_mut_freq = cdr2_imgt_nb_of_mutations / cdr2_imgt_nb_of_nucleotides * 100,
         cdr3_nt_mut_freq = cdr3_imgt_nb_of_mutations / cdr3_imgt_nb_of_nucleotides * 100,
         v_nt_mut_freq = v_region_nb_of_mutations / v_region_nb_of_nucleotides * 100,
         fr1_nonsilent_mut_ratio = fr1_imgt_nb_of_nonsilent_mutations / fr1_imgt_nb_of_mutations, ## Nonsilent mutation ratios
         fr2_nonsilent_mut_ratio = fr2_imgt_nb_of_nonsilent_mutations / fr2_imgt_nb_of_mutations,
         fr3_nonsilent_mut_ratio = fr3_imgt_nb_of_nonsilent_mutations / fr3_imgt_nb_of_mutations,
         cdr1_nonsilent_mut_ratio = cdr1_imgt_nb_of_nonsilent_mutations / cdr1_imgt_nb_of_mutations,
         cdr2_nonsilent_mut_ratio = cdr2_imgt_nb_of_nonsilent_mutations / cdr2_imgt_nb_of_mutations,
         cdr3_nonsilent_mut_ratio = cdr3_imgt_nb_of_nonsilent_mutations / cdr3_imgt_nb_of_mutations,
         v_nonsilent_mut_ratio = v_region_nb_of_nonsilent_mutations / v_region_nb_of_mutations,
         fr1_transitions = (fr1_imgt_a_into_g + fr1_imgt_g_into_a + fr1_imgt_c_into_t + fr1_imgt_t_into_c), ## Transition vs transversion counts
         fr1_transversions = fr1_imgt_nb_of_mutations - fr1_transitions,
         fr2_transitions = (fr2_imgt_a_into_g + fr2_imgt_g_into_a + fr2_imgt_c_into_t + fr2_imgt_t_into_c),
         fr2_transversions = fr2_imgt_nb_of_mutations - fr2_transitions,
         fr3_transitions = (fr3_imgt_a_into_g + fr3_imgt_g_into_a + fr3_imgt_c_into_t + fr3_imgt_t_into_c),
         fr3_transversions = fr3_imgt_nb_of_mutations - fr3_transitions,
         cdr1_transitions = (cdr1_imgt_a_into_g + cdr1_imgt_g_into_a + cdr1_imgt_c_into_t + cdr1_imgt_t_into_c),
         cdr1_transversions = cdr1_imgt_nb_of_mutations - cdr1_transitions,
         cdr2_transitions = (cdr2_imgt_a_into_g + cdr2_imgt_g_into_a + cdr2_imgt_c_into_t + cdr2_imgt_t_into_c),
         cdr2_transversions = cdr2_imgt_nb_of_mutations - cdr2_transitions,
         cdr3_transitions = (cdr3_imgt_a_into_g + cdr3_imgt_g_into_a + cdr3_imgt_c_into_t + cdr3_imgt_t_into_c),
         cdr3_transversions = cdr3_imgt_nb_of_mutations - cdr3_transitions,
         total_non_junction_v_transitions = fr1_transitions + fr2_transitions + fr3_transitions + cdr1_transitions + cdr2_transitions,
         total_non_junction_v_transversions = fr1_transversions + fr2_transversions + fr3_transversions + cdr1_transversions + cdr2_transversions,
         fr1_transition_ratio = fr1_transitions / fr1_imgt_nb_of_mutations, ## Transition ratios
         fr2_transition_ratio = fr2_transitions / fr2_imgt_nb_of_mutations,
         fr3_transition_ratio = fr3_transitions / fr3_imgt_nb_of_mutations,
         cdr1_transition_ratio = cdr1_transitions / cdr1_imgt_nb_of_mutations,
         cdr2_transition_ratio = cdr2_transitions / cdr2_imgt_nb_of_mutations,
         cdr3_transition_ratio = cdr3_transitions / cdr3_imgt_nb_of_mutations,
         non_junction_transition_ratio = total_non_junction_v_transitions / (fr1_imgt_nb_of_mutations + fr2_imgt_nb_of_mutations + fr3_imgt_nb_of_mutations + cdr1_imgt_nb_of_mutations + cdr2_imgt_nb_of_mutations),
         fr1_aa_mut_freq = fr1_imgt_nb_of_aa_changes / fr1_imgt_nb_of_aa, ## Amino acid conversion frequency
         fr2_aa__mut_freq = fr2_imgt_nb_of_aa_changes / fr2_imgt_nb_of_aa,
         fr3_aa_mut_freq = fr3_imgt_nb_of_aa_changes / fr3_imgt_nb_of_aa,
         cdr1_aa_mut_freq = cdr1_imgt_nb_of_aa_changes / cdr1_imgt_nb_of_aa,
         cdr2_aa_mut_freq = cdr2_imgt_nb_of_aa_changes / cdr2_imgt_nb_of_aa,
         cdr3_aa_mut_freq = cdr3_imgt_nb_of_aa_changes / cdr3_imgt_nb_of_aa,
         v_aa_mut_freq = (v_region_nb_of_aa / v_region_nb_of_aa_changes),
         contains_avy_in_fr1 = grepl("AVY", aa_fr1), ## Detects AVY in FR1
         contains_qw_in_fr1 = grepl("QW", aa_fr1)) ## Detects QW in FR1
  
  if (processed_hotspot_mutation_summary != FALSE) {
    
    combined_tibble <- combined_tibble %>% 
      mutate(fr1_pct_mutations_in_hotspots = fr1_hotspot_muts / fr1_imgt_nb_of_mutations * 100,
         fr2_pct_mutations_in_hotspots = fr2_hotspot_muts / fr2_imgt_nb_of_mutations * 100,
         fr3_pct_mutations_in_hotspots = fr3_hotspot_muts / fr3_imgt_nb_of_mutations * 100,
         cdr1_pct_mutations_in_hotspots = cdr1_hotspot_muts / cdr1_imgt_nb_of_mutations * 100,
         cdr2_pct_mutations_in_hotspots = cdr2_hotspot_muts / cdr2_imgt_nb_of_mutations * 100,
         cdr3_pct_mutations_in_hotspots = cdr3_hotspot_muts / cdr3_imgt_nb_of_mutations * 100,
         total_non_junction_pct_mutations_in_hotspots = total_non_junction_hotspot_muts / (fr1_imgt_nb_of_mutations + 
                                                                                             fr2_imgt_nb_of_mutations + 
                                                                                             fr3_imgt_nb_of_mutations + 
                                                                                             cdr1_imgt_nb_of_mutations + 
                                                                                             cdr2_imgt_nb_of_mutations),
         fr1_pct_hotspots_mutated = fr1_hotspot_muts / fr1_hotspot_loci * 100,
         fr2_pct_hotspots_mutated = fr2_hotspot_muts / fr2_hotspot_loci * 100,
         fr3_pct_hotspots_mutated = fr3_hotspot_muts / fr3_hotspot_loci * 100,
         cdr1_pct_hotspots_mutated = cdr1_hotspot_muts / cdr1_hotspot_loci * 100,
         cdr2_pct_hotspots_mutated = cdr2_hotspot_muts / cdr2_hotspot_loci * 100,
         cdr3_pct_hotspots_mutated = cdr3_hotspot_muts / cdr3_hotspot_loci * 100,
         total_non_junction_pct_hotspots_mutated = total_non_junction_hotspot_muts / total_non_junction_hotspot_loci * 100)
  }
  
  ##Bring in processed CellRanger annotations
  processed_cellranger_annotations <- read_csv(processed_cellranger_annotations_file)
  
  #Merge in the data
  merged_tibble <- left_join(processed_cellranger_annotations, combined_tibble, by = c("indexed_contig_id" = "sequence_id"))
  
  ##pull out the isotypes from the heavy chains
  isotypes <- merged_tibble %>% filter(chain_type == "igh") %>% select(indexed_barcode, isotype = c_gene)
  isotypes$isotype <- tolower(gsub("IGH(.*)", "\\1", isotypes$isotype))
  merged_tibble <- left_join(merged_tibble, isotypes, by = "indexed_barcode")
  
  ##Arrange the tibble for ease of future processing.
  if (processed_hotspot_mutation_summary != FALSE) {
    
    merged_tibble <- merged_tibble %>%
    select(indexed_contig_id:seq,
           chain_type, 
           isotype,
           reads:j_region_identity_pct,
           nt_vdj:aa_gapped_j,
           fr1_nt_mut_freq:contains_qw_in_fr1,
           "v_region_nb_of_positions_including_imgt_gaps_(nt)":cdr3_imgt_very_dissimilar,
           fr1_hotspot_loci:total_non_junction_hotspot_muts,
           fr1_pct_mutations_in_hotspots:total_non_junction_pct_hotspots_mutated)
    
  } else {
    
    merged_tibble <- merged_tibble %>%
    select(indexed_contig_id:seq,
           chain_type, 
           isotype,
           reads:j_region_identity_pct,
           nt_vdj:aa_gapped_j,
           fr1_nt_mut_freq:contains_qw_in_fr1,
           "v_region_nb_of_positions_including_imgt_gaps_(nt)":cdr3_imgt_very_dissimilar)
    
  }
  
  ##Return the table
  setwd(imgt_summary_file_directory)
  write_csv(merged_tibble, paste("superrseq_processed_imgt_output.index_", merged_tibble$sample_index[1], ".csv", sep = ""))
  
  if (return_df == TRUE) {
  return(merged_tibble)
  } else {
    print("Processing complete!")
  }
  
  }### End, process_imgt_output()