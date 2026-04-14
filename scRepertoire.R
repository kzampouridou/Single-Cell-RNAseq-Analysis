library(scRepertoire)
library(writexl)

# 1. SETUP
# Folder names containing your 'filtered_contig_annotations.csv' files
samples <- c("Sample_BCR")


# 2. DEFINE FUNCTION
# Function to process one sample
process_sample <- function(sample_name, threshold = 0.85) {
  
  message("--- Processing Sample: ", sample_name, " ---")
  
  # Construct file path to 10X filtered contigs
  input_path <- file.path(sample_name, "filtered_contig_annotations.csv")
  
  if (!file.exists(input_path)) {
    stop("Input file not found at: ", input_path)
  }
  
  # Load data
  contig_list <- loadContigs(input_path, format = "10X")
  
  # Run combineBCR "Matches Heavy and Light chains for each individual cell"
  combined <- combineBCR(contig_list,
                         samples = sample_name,
                         threshold = threshold,
                         )
  
  # Extract dataframe
  export_df <- combined[[1]]
  
    # Optional: preview key columns
      print(head(export_df[, c("barcode", "CTstrict", "IGH", "cdr3_aa1")]))
  
    # Optional: print cluster levels
      print(levels(as.factor(export_df$CTstrict)))
  
  # Save output
      output_fn <- paste0("BCR_Combined_", sample_name, ".xlsx")
      write_xlsx(export_df, output_fn)
      
      message("Successfully saved: ", output_fn)
  
  return(combined)
}

# 3. EXECUTION

# Apply to all samples
results_list <- lapply(samples, process_sample)

  # Optional: name the list
    names(results_list) <- samples


    