#' assign_lineages
#' 
#' The function takes a .db file by default, and adds logical columns to aid hybrid sequence filtering
#' @param data A flexible data input. Takes a single .csv file, a directory of .csv files, a single data frame, or a list of data frames
#' @param pct_homology The amount of CDR3 homology required to draw a lineage connection between two sequences
#' @param return_df Returns a data frame at the end of the function instead of writing a file
#' @export

assign_lineages <- function(data, pct_homology = 0.85, return_df = FALSE) {
   
  require(Biostrings)
  require(dplyr)
  require(readr)
  
  ###Identify type of data
  if ("data.frame" %in% class(data)) {
    data <- data
  } else if (class(data) == "list"){
    data <- do.call(rbind, data)
  } else if (grepl("\\.csv", data)) {
    
    ## Extract the directory for the single .csv
    subdirectory_path <- gsub("superrseq.*", "", data)
    setwd(subdirectory_path)
    
    data <- read_csv(data)
    
  } else if (dir.exists(data)) {
    setwd(data)
    file_list <- list.files(data)
    data_list <- vector("list", length = length(file_list))
    for (i in 1:length(file_list)) {
      data_list[[i]] <- read_csv(file_list[i])
    }
    data <- bind_rows(data_list)
  }

  ### Lineage assignment function
  
  lin_assign <- function(data) {
    
    ### Populate sequences
    cdr3_seqs <- data[["cdr3_seq"]]
    
    ### Create the initial distance matrix
    seq_dist_matrix <- stringDist(unique(cdr3_seqs), method = "hamming", upper = TRUE, diag = TRUE)
    percent_homology_matrix <- 1 - (as.matrix(seq_dist_matrix) / nchar(cdr3_seqs)[1])
    colnames(percent_homology_matrix) <- rownames(percent_homology_matrix) <- unique(cdr3_seqs)
    
    ### Create the sequence frequency table
    seq_freq_table <- sort(table(cdr3_seqs), decreasing = TRUE)
    
    ### iterative assignment of lineages
    
    while (length(seq_freq_table) >= 1) {
      
      ### identify most abundant sequence
      highest_freq_seq <- names(seq_freq_table)[1]
      
      ### Isolate the row of the highest-frequency seq
      seq_row <- as.vector(percent_homology_matrix[highest_freq_seq,,drop = FALSE])
      names(seq_row) <- colnames(percent_homology_matrix[highest_freq_seq,,drop = FALSE])
      
      ### Identify sequences within homology threshold
      homology_names <- names(which(seq_row >= pct_homology))
      
      ### Find homologous sequences in the data
      lineage_index <- which(data[["cdr3_seq"]] %in% homology_names)
      
      ### Assign lineage_id with the lineage id counter
      data$lineage_id[lineage_index] <- lineage_counter
      
      ### Prepare indexes for sequence removal from the table and matrix
      table_index <- which(names(seq_freq_table) %in% homology_names)
      matrix_index <- which(colnames(percent_homology_matrix) %in% homology_names)
      
      ### Remove those sequences from the table and matrix
      seq_freq_table <- seq_freq_table[-table_index]
      percent_homology_matrix <- percent_homology_matrix[-matrix_index,-matrix_index, drop = FALSE]
      
      ### Click over the lineage assignment counter
      
      lineage_counter <<- lineage_counter + 1
      
    } ### End "while" loop
    
    return(data)
    
  } ### End lin_assign function
  
  ### Make a trimmed data frame for processing (removing light chains)
  trimmed_data <- data %>% 
    filter(chain_type == "igh") %>%
    dplyr::select(indexed_barcode, v_gene, j_gene, cdr3_seq = nt_cdr3, cdr3_len = nt_cdr3_length) 

  ### Extract CDR3 length, and add a placeholder for lineage id  
  trimmed_data <- trimmed_data %>%
    mutate(lineage_id = rep(NA))
  
  ### Split the list on V, J, and CDR3 length
  split_list <- split(trimmed_data, list(factor(trimmed_data[["v_gene"]]), factor(trimmed_data[["j_gene"]]), factor(trimmed_data[["cdr3_len"]])), drop = TRUE)
  
  ### Initialize lineage counter
  lineage_counter <- 1

  ### Assignment visualizer
  for (i in seq_along(split_list)) {
    
  if (i %% 100 == 0) {
    
    print(paste(i, "out of", length(split_list), "v/j/cdr3 length families assigned"))
    
  }  
  
  ### assign the lineages  
  split_list[[i]] <- lin_assign(split_list[[i]])
    
  }
  
  ### collapse the list
  assigned_lins <- bind_rows(split_list)
  
  ### assign the new lineage ids
  assigned_lins <- assigned_lins %>% select(indexed_barcode, lineage_id)
  
  ### join the new lineage ids to the data
  data <- left_join(data, assigned_lins, by = "indexed_barcode")
  
  ### Add in pat_lins when available
  if ("patient_id" %in% names(data)) {
    
    data$pat_lin <- paste(data$patient_id, data$lineage_id, sep = "_")
    
  }

  ### Create an index .txt file
  split_index_list <- split(data, factor(data$sample_index))

  index_reference_key <- vector("list", length = length(split_index_list))
  
  split_lineage_list <- vector("list", length = length(split_index_list))
  names(split_lineage_list) <- unique(data$sample_index)
  
  start_metadata <- which(names(split_index_list[[1]]) == "sample_index")

  end_metadata <- which(names(split_index_list[[1]]) == "seq") - 1

  
  for (i in 1:length(split_index_list)) {
    
    index_reference_key[[i]] <- split_index_list[[i]][1,] %>% select(start_metadata:end_metadata) 
    
    split_lineage_list[[i]] <- unique(split_index_list[[i]]$lineage_id)
    data[[paste("lineage_in_index_", names(split_lineage_list)[i], sep = "")]] <- FALSE
    
  }
  
  indices <- unique(data$sample_index)
  
  for (i in 1:nrow(data)) {
    
    for (j in 1:length(indices)) {
      
      if (data[i,"lineage_id"] %in% split_lineage_list[[j]]) {
        
        data[i,paste("lineage_in_index_", names(split_lineage_list)[j], sep = "")] <- TRUE
        
      }
      
    }
    
  }
  
  index_reference_key <- bind_rows(index_reference_key)
  
  write_csv(index_reference_key, "index_reference_key.csv")
  
  write_csv(data, "lineage_processed_data.csv")
  
  ### return the data
   if (return_df == TRUE) {
  return(data)
  } else {
    print("Processing complete!")
  }

}