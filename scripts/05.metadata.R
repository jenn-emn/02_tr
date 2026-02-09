#!/usr/bin/env Rscript

# ===============================================
# Script to filter metadata of the HPRC samples
# ===============================================

library(tidyverse)



# ===============================================
# Verifying the input call
# ===============================================

option_list <- list(
  make_option(c("--out"), type = "character", default = NULL,
              help = "path to the output folder where the extracted MHC VCFs will be saved", 
              metavar = "path"),
  make_option(c("--name"), type = "character", default = NULL,
              help = "string to identify the job, e.g., 'trios' or 'trios_hla-mapper'", 
              metavar = "string")
)

usage_msg <- "%prog --out /path/to/output --name trios_analysis\n       %prog --out /home/jennifer/02_datas/01_intermediate --name test_hla-mapper"
opt_parser <- OptionParser(option_list = option_list, usage = usage_msg)
opt <- parse_args(opt_parser)

# Validate required arguments
if (is.null(opt$out) || is.null(opt$name)) {
  print_help(opt_parser)
  stop("Erro: Argumentos --out e --name sao obrigatorios.\n", call. = FALSE)
}

# Validate output directory and path
if (!dir.exists(opt$out)) {
  stop(paste("The directory path ('", opt$out, "') is not an existing directory."), call. = FALSE)
}

# Directory of actual Rscript
initial_options <- commandArgs(trailingOnly = FALSE)
script_path <- sub("--file=", "", initial_options[grep("--file=", initial_options)])
script_dir <- ifelse(length(script_path) > 0, dirname(normalizePath(script_path)), getwd())

# Assigning variables
path_out <- opt$out
name_job <- opt$name
cat("Job Name:", name_job, "\n")
cat("Output Path:", path_out, "\n")
cat("Script Dir:", script_dir, "\n")



# ===============================================
# Starting the analysis
# ===============================================

# intermediate paths
path_int <- file.path(path_out, name_job)
path_meta <- file.path(path_int, "metadata")
dir.create(path_meta, recursive = TRUE)


# download the metadata of the HPRC samples
cmd = "
    wget \
        -O ${path_meta}/hprc_release2_sample_metadata.csv \
        https://raw.githubusercontent.com/human-pangenomics/hprc_intermediate_assembly/main/data_tables/sample/hprc_release2_sample_metadata.csv"
system(cmd)

meta_df <- read_csv(file.path(path_meta, "hprc_release2_sample_metadata.csv")) %>%
  select(sample_id, population_descriptor, population_abbreviation, paternal_id, maternal_id, sex) %>%
  rename(Sample = sample_id, Population = population_descriptor, Pop_Abbr = population_abbreviation)

# join
samples_df <- data.frame(
    Sample = c("HG00733", "HG01109", "HG01243","HG02055", "HG02080", "HG02145",
    "HG02723", "HG02818", "HG03098","HG03486", "HG03492", "NA18906","NA19240", "NA20129") ) %>%
    left_join(meta_df, by = "Sample")

write.table(samples_df, file = file.path(path_meta, "hprc_samples_metadata.tsv"), quote = FALSE, row.names = FALSE, sep = "\t")


# end