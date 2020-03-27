#' process_cellranger_outputs
#'
#' @param cellranger_all_contig_annotations_csv Quoted file path of the CellRanger unfiltered annotations .csv
#' @param cellranger_all_contig_fasta Quoted file path of the CellRanger unfiltered contig .fasta
#' @param sample_index Internal sample index to be appended as metadata to the dataset
#' @param filter_is_cell Use the is_cell CellRanger filter in generating the high-confidence .fasta
#' @param filter_high_confidence Use the high_confidence CellRanger filter in generating the high-confidence .fasta
#' @param filter_is_full_length Use the full_length CellRanger filter in generating the high-confidence .fasta
#' @param filter_is_productive Use the is_productive CellRanger filter in generating the high-confidence .fasta
#' @param filter_read_count Minimum reads requirement in generating the high-confidence .fasta
#' @param filter_umi_count Minimum UMI requirement in generating the high-confidence .fasta
#' @param require_light_chain Require both a single heavy chain and a single light chain - otherwise a lack of light chain is permitted
#'
#' @return
#' @export
#'
#' @examples


process_cellranger_outputs <- function(cellranger_all_contig_annotations_csv, 
                                       cellranger_all_contig_fasta,
                                       metadata_csv,
                                       sample_index = FALSE,
                                       filter_is_cell = TRUE, 
                                       filter_high_confidence = TRUE, 
                                       filter_is_full_length = TRUE, 
                                       filter_is_productive = TRUE, 
                                       filter_read_count = 50,
                                       filter_umi_count = 2,
                                       require_light_chain = FALSE) {
  
##Dependencies
require(dplyr)
require(readr)
require(stringr)
  
## Create a subdirectory for processed outputs
  subdirectory_path <- gsub("all_contig_annotations.csv", "superrseq_processed_cellranger_outputs/", cellranger_all_contig_annotations_csv)
  dir.create(subdirectory_path)
  setwd(subdirectory_path)

### Function reads a .fasta and returns a data framewith the sequence and related info
read_fasta <- function(fasta_file) {
  
  raw_file <- read_file(fasta_file) ###read raw file
  raw_matrix <- str_split_fixed(raw_file, ">", Inf) ###splits on the greater than
  raw_chr_strings <- raw_matrix[,2:ncol(raw_matrix)] ###trims the first empty entry
  
  info_strings <- gsub("^(.*?)\n.*", "\\1", raw_chr_strings) ### store the info line
  
  sequence_strings <- gsub("^.*?\n(.*)", "\\1", raw_chr_strings) ### store the sequence string
  returns_removed <- gsub("\n", "", sequence_strings) ### remove all returns
  carriage_returns_removed <- gsub("\r", "", returns_removed) ### remove carriage returns
  
  tibble("seq_info" = info_strings, "seq" = carriage_returns_removed) ### make a data frame
  
}

### Function takes a data frame containing ids and sequences, and returns a fasta formatted file
write_fasta <- function(df, seq_id, seq, file = "fasta.txt") {
  
  fasta_matrix <- matrix(nrow = nrow(df), ncol = 1)
  
  for (i in 1:nrow(df)) {
    fasta_matrix[i,1] <- paste(">", df[[seq_id]][i], "\n", df[[seq]][i], "\n", sep = "")
  }
  
  write(fasta_matrix, file = file)
  
}
  
##Load in the files
annotation_file <- read_csv(cellranger_all_contig_annotations_csv)
contigs <- read_fasta(cellranger_all_contig_fasta)
metadata_df <- read_csv(metadata_csv)

## Index the annotations file
annotation_file$sample_index <- sample_index

## Tack on the metadata to the annotations file
annotation_file <- left_join(annotation_file, metadata_df, by = "sample_index")

##Index the sample barcodes
raw_barcodes <- gsub("^([A-Z]*).*", "\\1", annotation_file$barcode)
annotation_file$indexed_barcode <- paste(raw_barcodes, annotation_file$sample_index, sep = "-")

##Extract the contig_id number
annotation_file$raw_contig_id <- gsub(".*(_contig_[0-9]*)$", "\\1", annotation_file$contig_id)
annotation_file$indexed_contig_id <- paste(annotation_file$indexed_barcode, annotation_file$raw_contig_id, sep = "")
  
##Update the contig fasta indexes
contigs$seq_info <- gsub("-[0-9]*_", paste("-", sample_index, "_", sep = ""), contigs$seq_info)

##Write the updated .fasta with indexed contig designations  
write_fasta(contigs, "seq_info", "seq", paste("processed_contigs.index_", sample_index, ".fasta", sep = "")) 

###Filter the annotations file

## is_cell filter
if (filter_is_cell == TRUE) {
  annotation_filter1 <- annotation_file %>% filter(is_cell == TRUE)
} else {
  annotation_filter1 <- annotation_file
}

## high_confidence filter
if (filter_high_confidence == TRUE) {
  annotation_filter2 <- annotation_filter1 %>% filter(high_confidence == TRUE)
} else {
  annotation_filter2 <- annotation_filter1
}

## full_length filter
if (filter_is_full_length == TRUE) {
  annotation_filter3 <- annotation_filter2 %>% filter(full_length == TRUE)
} else {
  annotation_filter3 <- annotation_filter2
}

## is_productive filter
if (filter_is_productive == TRUE) {
  annotation_filter4 <- annotation_filter3 %>% filter(productive == "True")
} else {
  annotation_filter4 <- annotation_filter3
}
 
## read count and UMI filter
annotation_filtered <- annotation_filter4 %>% filter(reads >= filter_read_count, umis >= filter_umi_count)

## Generate a summary table with the number of filtered heavy and light chains per barcode
annotation_chain_summary <- annotation_filtered %>%
  group_by(indexed_barcode) %>% 
  dplyr::summarize(n_seqs = n(),
            n_heavy_chains = sum(chain == "IGH"),
            n_light_chains = sum(chain != "IGH"))

## Return a df identifying the number of barcodes containing n heavy and n light chains
annotation_summary_df <- data.frame(table(annotation_chain_summary$n_heavy_chains, annotation_chain_summary$n_light_chains))

## Name the df appropriately
names(annotation_summary_df) <- c("n_heavy_chains", "n_light_chains", "n_barcodes")

## Identify barcodes with the correct number of heavy and light chains
if (require_light_chain == TRUE) {

  annotation_filtered_singlet_barcodes <- annotation_chain_summary %>% 
  filter(n_heavy_chains == 1, n_light_chains == 1) %>%
  .$indexed_barcode
  
} else {
  
  annotation_filtered_singlet_barcodes <- annotation_chain_summary %>% 
  filter(n_heavy_chains == 1, n_light_chains <= 1) %>%
  .$indexed_barcode
  
}
 
## Filter based on chain filter barcodes
annotation_filtered <- annotation_filtered %>% filter(indexed_barcode %in% annotation_filtered_singlet_barcodes)

## Add sequences to the annotations file
annotation_filtered <- left_join(annotation_filtered, contigs, by = c("indexed_contig_id" = "seq_info"))

## Select columns for later use following IMGT processing
final_annotations <- annotation_filtered %>%
  select(indexed_contig_id,
         indexed_barcode,
         sample_index:indexed_barcode,
         seq,
         reads,
         umis,
         c_gene)

### Write fasta file for filtered sequences
write_fasta(final_annotations, "indexed_contig_id", "seq", paste("filtered_singlet_sequences.index_", sample_index, ".fasta", sep = ""))

### Write filtered annotations file
write_csv(final_annotations, paste("filtered_singlet_annotations.index_", sample_index, ".csv", sep = ""))

### Write a .csv showing the heavy and light chain data for diagnostic purposes
write_csv(annotation_summary_df, paste("full_chain_summary.index_", sample_index, ".csv", sep = ""))

### Print the number of filtered unique barcodes to the console
print(paste(length(unique(final_annotations$indexed_barcode)), "high quality cells identified", sep = " "))

}