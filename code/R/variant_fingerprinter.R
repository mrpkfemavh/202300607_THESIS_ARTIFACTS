library(tidyverse)
library(readr)
library(ggplot2)
library(stringr)
library(DescTools)
library(splitstackshape)
library(foreach)
library(doParallel)
library(lubridate)

dir_data_computed <-
  "/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_computed"

dir_data_references <-
  "/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_raw/reference_variants" # /mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_raw/reference_variants/MN908947.3.gff3

dir_data_generated <-
  "/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_generated"

file_all_snp <-
  paste0(dir_data_computed, "/", "MN908947.3_all_snp.data")
file_all_gff3_MN908947_3 <-
  paste0(dir_data_references, "/", "MN908947.3.gff3") # /mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_raw/reference_variants/MN908947.3.gff3

file_all_metadata <- paste0(dir_data_generated, "/", "metadata.csv") # "/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_generated/metadata.csv" # "/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_computed/column.header.csv"

file_all_collection <- paste0(dir_data_generated, "/", "collection_date.csv") # "/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_generated/collection_date.csv"

collection_date <- read_csv(file_all_collection,
                            col_types = cols(DATE_BEST_COLLECTION_SUBMISSION = col_character())) %>% tibble()

collection_date$DATE_BEST_COLLECTION_SUBMISSION <-
  as_date(collection_date$DATE_BEST_COLLECTION_SUBMISSION)

collection_date
write.csv(x = collection_date, file = paste0(dir_data_computed, "/", "collection_date.csv"))
saveRDS(collection_date, file = paste0(dir_data_computed, "/", "collection_date.rds"))

metadata <-
  read_csv(
    "/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_generated/metadata.csv",
    col_types = cols(
      METADATA_COMPUTED_TOTAL_COPIES = col_integer(),
      METADATA_COPIES_TOTAL = col_integer(),
      METADATA_VARIANT_FRACTION_C1_AUSTRALIA = col_double(),
      METADATA_VARIANT_FRACTION_C2_ANCESTRAL = col_double(),
      METADATA_VARIANT_FRACTION_C14_ALPHA = col_double(),
      METADATA_VARIANT_FRACTION_C29_DELTA = col_double(),
      METADATA_VARIANT_FRACTION_C48_OMICRON = col_double(),
      METADATA_VARIANT_COPIES_C1_MT007544.1 = col_integer(),
      METADATA_VARIANT_COPIES_C2_MN908947.3 = col_integer(),
      METADATA_VARIANT_COPIES_C14_EPI_ISL_710528 = col_integer(),
      METADATA_VARIANT_COPIES_C29_EPI_ISL_2693246 = col_integer(),
      METADATA_VARIANT_COPIES_C48_EPI_ISL_6841980 = col_integer(),
      METADATA_G_RUNS_HOW_MANY = col_integer(),
      METADATA_G_LANES_FORWARD = col_integer(),
      METADATA_G_LANES_REVERSE = col_integer(),
      METADATA_G_LANES_INDEX = col_integer()
    )
  ) %>% select(contains("METADATA_")) %>% tibble()

metadata$METADATA_COMPUTED_TOTAL_COPIES <-
  as.integer(metadata$METADATA_COMPUTED_TOTAL_COPIES)

metadata$METADATA_COPIES_TOTAL <-
  as.integer(metadata$METADATA_COPIES_TOTAL)

metadata$METADATA_VARIANT_COPIES_C1_MT007544.1 <-
  as.integer(metadata$METADATA_VARIANT_COPIES_C1_MT007544.1)

metadata$METADATA_VARIANT_COPIES_C2_MN908947.3 <-
  as.integer(metadata$METADATA_VARIANT_COPIES_C2_MN908947.3)

metadata$METADATA_VARIANT_COPIES_C14_EPI_ISL_710528 <-
  as.integer(metadata$METADATA_VARIANT_COPIES_C14_EPI_ISL_710528)

metadata$METADATA_VARIANT_COPIES_C29_EPI_ISL_2693246 <-
  as.integer(metadata$METADATA_VARIANT_COPIES_C29_EPI_ISL_2693246)

metadata$METADATA_VARIANT_COPIES_C48_EPI_ISL_6841980 <-
  as.integer(metadata$METADATA_VARIANT_COPIES_C48_EPI_ISL_6841980)

metadata$identifier <- metadata$METADATA_NAME_IDENTIFIER

metadata
write.csv(x = metadata, file = paste0(dir_data_computed, "/", "metadata_computed.csv"))
saveRDS(metadata, file = paste0(dir_data_computed, "/", "metadata_computed.rds"))

MN908947_3_all_snp <-
  read_delim(
    file_all_snp, # "/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/locate_snps/succesful_run/MN908947.3_all_snp.data", #/home/bgilot/Downloads/SarsCovControls/MN908947.3_all_snp_R.data",
    "\t",
    escape_double = FALSE,
    col_types = cols(
      BLANK_COLUMN = col_skip(),
      REFERENCE_POSITION = col_integer(),
      QUERY_POSITION = col_integer(),
      BUFF = col_integer(),
      DIST = col_integer(),
      REFERENCE_LENGTH = col_integer(),
      QUERY_LENGTH = col_integer(),
      REFERENCE_FRM = col_skip(),
      QUERY_FRM = col_skip(),
      SUBMISSION = col_character()
    ),
    trim_ws = TRUE
  ) %>% tibble()

MN908947_3_all_snp
write.csv(x = MN908947_3_all_snp, file = paste0(dir_data_computed, "/", "MN908947_3_all_snp.csv"))
saveRDS(MN908947_3_all_snp, file = paste0(dir_data_computed, "/", "MN908947_3_all_snp.rds"))

gff_MN908947_3 <-
  read_delim(
    file_all_gff3_MN908947_3, # "/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/locate_snps/variants_data/MN908947.3.gff3", #"Downloads/SarsCovControls/MN908947.3.gff3",
    "\t",
    escape_double = FALSE,
    col_names = FALSE,
    trim_ws = TRUE,
    skip = 2
  ) %>% tibble()

names(gff_MN908947_3) <-
  c("seqid",
    "source",
    "type",
    "start",
    "end",
    "score",
    "strand",
    "phase",
    "attributes")

MN908947_3 <- cSplit(gff_MN908947_3, "attributes", ";")


att_ref <-
  MN908947_3 %>% filter(type == "CDS" |
                          type == "five_prime_UTR") %>% select(start, end, attributes_02) %>% tibble()
att_ref$attributes_02 <-
  str_replace(att_ref$attributes_02, "gbkey=", "Parent=gene-")

att_ref <-
  cSplit(att_ref, "attributes_02", "=gene-") %>% select(start, end, attributes_02_2)

names(att_ref) <- c("start", "end", "CDS_gene")

att_ref
write.csv(x = att_ref, file = paste0(dir_data_computed, "/", "att_ref.csv"))
saveRDS(att_ref, file = paste0(dir_data_computed, "/", "att_ref.rds"))

SNPs_fingerprint <-
  MN908947_3_all_snp %>% filter(QUERY != ".") %>% filter(REFERENCE != ".") %>% arrange(REFERENCE_POSITION) %>% tibble()

SNPs_fingerprint <- left_join(SNPs_fingerprint, collection_date) %>% tibble()

SNPs_fingerprint$gff_start <- as.integer(NA)
SNPs_fingerprint$gff_end <- as.integer(NA)
SNPs_fingerprint$gff_gene <- as.character(NA)
SNPs_fingerprint$TRANSLATED <- as.character(NA)

for (i in 1:nrow(SNPs_fingerprint)) {
  check_position <- SNPs_fingerprint[i,]$REFERENCE_POSITION
  result <-
    att_ref %>% filter(start <= check_position, end >= check_position)
  how_many <- nrow(result)
  debut <-
    (att_ref %>% filter(start <= check_position, end >= check_position))$start
  fin <-
    (att_ref %>% filter(start <= check_position, end >= check_position))$end
  quoi <-
    (att_ref %>% filter(start <= check_position, end >= check_position))$CDS_gene
  
  translated <- as.character("TRANSLATED")
  if (how_many != 1) {
    print(check_position)
    print(result)
    
    debut <- as.integer(-1)
    fin <- as.integer(-2)
    quoi <- as.character("UNTRANSLATED")
    translated <- as.character("UNTRANSLATED")
  }
  SNPs_fingerprint[i,]$gff_start <- as.integer(debut)
  SNPs_fingerprint[i,]$gff_end <- as.integer(fin)
  SNPs_fingerprint[i,]$gff_gene <- as.character(quoi)
  SNPs_fingerprint[i,]$TRANSLATED <- as.character(translated)
}

SNPs_fingerprint

SNPs_MT007544_1_australia <-
  SNPs_fingerprint %>% filter(TRANSLATED == "TRANSLATED",
                              REFERENCE != ".",
                              QUERY != ".",
                              QUERY_TAG == "MT007544.1") %>% arrange(REFERENCE_POSITION) %>% tibble()

SNPs_EPI_ISL_710528_alpha <-
  SNPs_fingerprint %>% filter(TRANSLATED == "TRANSLATED",
                              REFERENCE != ".",
                              QUERY != ".",
                              QUERY_TAG == "EPI_ISL_710528") %>% arrange(REFERENCE_POSITION) %>% tibble()

SNPs_EPI_ISL_2693246_delta <-
  SNPs_fingerprint %>% filter(
    TRANSLATED == "TRANSLATED",
    REFERENCE != ".",
    QUERY != ".",
    QUERY_TAG == "EPI_ISL_2693246"
  ) %>% arrange(REFERENCE_POSITION) %>% tibble()

SNPs_EPI_ISL_omicron <-
  SNPs_fingerprint %>% filter(
    TRANSLATED == "TRANSLATED",
    REFERENCE != ".",
    QUERY != ".",
    QUERY_TAG == "EPI_ISL_6841980"
  ) %>% arrange(REFERENCE_POSITION) %>% tibble()

genes <- c("S")
variants <- SNPs_fingerprint$QUERY_TAG %>% sort() %>% unique()
coord_start <- 21563
coord_end <- 25384

for_plot <-
  SNPs_fingerprint %>% filter(TRANSLATED == "TRANSLATED", REFERENCE != ".", QUERY != ".") %>% filter(gff_gene %in% genes) %>% arrange(REFERENCE_POSITION, QUERY_TAG) %>% tibble()

for_plot

for_plot$REFERENCE_POSITION %>% sort() %>% unique()

write.csv(x = SNPs_fingerprint, file = paste0(dir_data_computed, "/", "SNPs_fingerprint_all.csv")) # "/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/locate_snps/variants_data/SNPs_fingerprint_all.csv")
write.csv(x = SNPs_MT007544_1_australia, file = paste0(dir_data_computed, "/", "SNPs_MT007544_1_australia.csv")) # "/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/locate_snps/variants_data/SNPs_MT007544_1_australia.csv")
write.csv(x = SNPs_EPI_ISL_710528_alpha, file = paste0(dir_data_computed, "/", "SNPs_EPI_ISL_710528_alpha.csv")) # "/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/locate_snps/variants_data/SNPs_EPI_ISL_710528_alpha.csv")
write.csv(x = SNPs_EPI_ISL_2693246_delta, file = paste0(dir_data_computed, "/", "SNPs_EPI_ISL_2693246_delta.csv")) # "/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/locate_snps/variants_data/SNPs_EPI_ISL_2693246_delta.csv")
write.csv(x = SNPs_EPI_ISL_omicron, file = paste0(dir_data_computed, "/", "SNPs_EPI_ISL_omicron.csv")) # "/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/locate_snps/variants_data/SNPs_EPI_ISL_omicron.csv")
write.csv(x = for_plot, file = paste0(dir_data_computed, "/", "SNPs_for_plot.csv")) # "/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/locate_snps/variants_data/SNPs_for_plot.csv")

saveRDS(SNPs_fingerprint, file = paste0(dir_data_computed, "/", "SNPs_fingerprint_all.rds")) # "/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/locate_snps/variants_data/SNPs_fingerprint_all.rds")
saveRDS(SNPs_MT007544_1_australia, file = paste0(dir_data_computed, "/", "SNPs_MT007544_1_australia.rds")) # "/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/locate_snps/variants_data/SNPs_MT007544_1_australia.rds")
saveRDS(SNPs_EPI_ISL_710528_alpha, file = paste0(dir_data_computed, "/", "SNPs_EPI_ISL_710528_alpha.rds")) # "/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/locate_snps/variants_data/SNPs_EPI_ISL_710528_alpha.rds")
saveRDS(SNPs_EPI_ISL_2693246_delta, file = paste0(dir_data_computed, "/", "SNPs_EPI_ISL_2693246_delta.rds")) # "/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/locate_snps/variants_data/SNPs_EPI_ISL_2693246_delta.rds")
saveRDS(SNPs_EPI_ISL_omicron, file = paste0(dir_data_computed, "/", "SNPs_EPI_ISL_omicron.rds")) # "/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/locate_snps/variants_data/SNPs_EPI_ISL_omicron.rds")
saveRDS(for_plot, file = paste0(dir_data_computed, "/", "SNPs_for_plot.rds")) # "/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/locate_snps/variants_data/SNPs_for_plot.rds")

SNPs_fingerprint$REFERENCE_POSITION
for_plot$REFERENCE_POSITION
# SNPs_delta <- SNPs_EPI_ISL_2693246
# SNPs_alpha <- SNPs_EPI_ISL_710528
# SNPs_omicron <- SNPs_EPI_ISL_6841980
# SNPs_Australia <- SNPs_MT007544_1

genes <- c("S")
variants <- SNPs_fingerprint$QUERY_TAG %>% sort() %>% unique()

for_plot <-
  SNPs_fingerprint %>% filter(TRANSLATED == "TRANSLATED", REFERENCE != ".", QUERY != ".") %>% filter(gff_gene %in% genes) %>% arrange(REFERENCE_POSITION, QUERY_TAG)

for_plot$REFERENCE_POSITION %>% sort() %>% unique()

# counts_experimental
input_dir <-
  "/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_raw/counts_experimental"
input_list <-
  list.files(input_dir, full.names = TRUE,  pattern = '.tsv$')
files_to_read <-
  list.files(path = input_dir,
             pattern = "\\.tsv$",
             full.names = TRUE)
all_files <- lapply(files_to_read, function(x) {
  read.table(file = x,
             sep = '\t',
             header = TRUE)
})

counts_experimental <- bind_rows(all_files)

# counts_controls
input_dir <-
  "/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_raw/counts_controls"
input_list <-
  list.files(input_dir, full.names = TRUE,  pattern = '.tsv$')
files_to_read <-
  list.files(path = input_dir,
             pattern = "\\.tsv$",
             full.names = TRUE)
all_files <- lapply(files_to_read, function(x) {
  read.table(file = x,
             sep = '\t',
             header = TRUE)
})

counts_controls <- bind_rows(all_files)

counts_all <- rbind(counts_experimental, counts_controls)
counts_metadata_all <- left_join(counts_all, metadata)

sample_group <-
  counts_metadata_all$sample_group %>% unique() %>% sort()
sample_method <-
  counts_metadata_all$sample_method %>% unique() %>% sort()
sample_name <-
  counts_metadata_all$sample_name %>% unique() %>% sort()

sample_depths <-
  data.frame(matrix(
    ncol = 9,
    nrow = 0,
    dimnames = list(
      NULL,
      c(
        "sample_group",
        "sample_method",
        "sample_name",
        "sample_min",
        "sample_mean",
        "sample_median",
        "sample_max",
        "sample_sd",
        "sample_length"
      )
    )
  ))
class(sample_depths)

for (group in sample_group) {
  for (method in sample_method) {
    for (name in sample_name) {
      #print(counts_metadata_all %>% filter(sample_group == group, sample_method == method, sample_name == name) %>% select(depth) %>% summary())
      current_depths <-
        (
          counts_metadata_all %>% filter(
            sample_group == group,
            sample_method == method,
            sample_name == name
          )
        )$depth
      if (identical(current_depths, integer(0))) {
        sample_min <- as.integer(NA)
        sample_mean <- as.double(NA)
        sample_median <- as.double(NA)
        sample_max <- as.integer(NA)
        sample_sd <- as.double(NA)
        sample_length <- as.integer(NA)
      } else {
        sample_min <- as.integer(min(current_depths))
        sample_mean <- (mean(current_depths))
        sample_median <- (median(current_depths))
        sample_max <- as.integer(max(current_depths))
        sample_sd <- (sd(current_depths))
        sample_length <- as.integer(length(current_depths))
        
        print(paste("The group is", group))
        print(paste("The method is", method))
        print(paste("The name is", name))
        print(sample_min)
        print(sample_mean)
        print(sample_median)
        print(sample_max)
        print(sample_sd)
        print(sample_length)
        print("")
        sample_depths[nrow(sample_depths) + 1, ] <-
          list(
            group,
            method,
            name,
            sample_min,
            sample_mean,
            sample_median,
            sample_max,
            sample_sd,
            sample_length
          )
      }
    }
  }
}

counts_metadata_all_depths <-
  left_join(counts_metadata_all, sample_depths)

sample_depth_stats <-
  counts_metadata_all_depths %>% select(
    sample_group,
    sample_method,
    sample_name,
    contains("_CT_"),
    contains("_TOTAL_"),
    contains("METADATA_VARIANT"),
    sample_min,
    sample_mean,
    sample_median,
    sample_max,
    sample_length
  ) %>% unique() %>% tibble()

sample_depth_stats
write.csv(x = sample_depth_stats, file = paste0(dir_data_computed, "/", "sample_depth_stats.csv"))
saveRDS(sample_depth_stats, file = paste0(dir_data_computed, "/", "sample_depth_stats.rds"))


















for (group in sample_group) {
  print(paste("The year is", group))
}

for (method in sample_method) {
  print(paste("The year is", method))
}

for (name in sample_name) {
  print(paste("The year is", name))
}

for (group in sample_group) {
  for (method in sample_method) {
    for (name in sample_name) {
      print(paste("The year is", group))
      print(paste("The year is", method))
      print(paste("The year is", name))
      
    }
  }
}

for (method in sample_method) {
  print(paste("The year is", method))
}

for (name in sample_name) {
  print(paste("The year is", name))
}




df_all_rows %>% colnames()

df_all_rows %>% filter(between(pos, 21563, 25384))

snp_positions <-
  df_all_rows %>% filter(pos %in% for_plot$REFERENCE_POSITION)



df_all_rows$gini <- as.double(NA)
for (i_r in 1:nrow(df_all_rows)) {
  vec_gini <-
    c(
      df_all_rows[i_r,]$acount,
      df_all_rows[i_r,]$ccount,
      df_all_rows[i_r,]$gcount,
      df_all_rows[i_r,]$tcount,
      df_all_rows[i_r,]$ncount,
      df_all_rows[i_r,]$indelcount
    )
  df_all_rows[i_r,]$gini <- Gini(vec_gini, unbiased = FALSE)
}



df_all_rows %>% filter(between(pos, 21563, 25384))

snp_positions <-
  df_all_rows %>% filter(pos %in% for_plot$REFERENCE_POSITION) # %>% filter(depth >= 100) # %>% filter(sample_method == "unblasted")
snp_positions$gini <- as.double(NA)
snp_positions$fraction_reference <-
  as.double(snp_positions$refcount / snp_positions$depth * 100)
snp_positions$fraction_alternate <-
  as.double(snp_positions$altcount / snp_positions$depth * 100)
snp_positions$fraction_uncalled <-
  as.double(snp_positions$ncount / snp_positions$depth * 100)
snp_positions$fraction_indel <-
  as.double(snp_positions$indelcount / snp_positions$depth * 100)
for (i_r in 1:nrow(snp_positions)) {
  vec_gini <-
    c(
      snp_positions[i_r,]$acount,
      snp_positions[i_r,]$ccount,
      snp_positions[i_r,]$gcount,
      snp_positions[i_r,]$tcount,
      snp_positions[i_r,]$ncount,
      snp_positions[i_r,]$indelcount
    )
  snp_positions[i_r,]$gini <- Gini(vec_gini, unbiased = FALSE)
}

write.csv(x = snp_positions, file = "/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/locate_snps/variants_data/snp_positions_all.csv")
saveRDS(snp_positions, file = "/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/locate_snps/variants_data/snp_positions_all.rds")

snp_positions <-
  df_all_rows %>% filter(pos %in% for_plot$REFERENCE_POSITION) %>% filter(depth >= 100) # %>% filter(sample_method == "unblasted")
snp_positions$gini <- as.double(NA)
snp_positions$fraction_reference <-
  as.double(snp_positions$refcount / snp_positions$depth * 100)
snp_positions$fraction_alternate <-
  as.double(snp_positions$altcount / snp_positions$depth * 100)
snp_positions$fraction_uncalled <-
  as.double(snp_positions$ncount / snp_positions$depth * 100)
snp_positions$fraction_indel <-
  as.double(snp_positions$indelcount / snp_positions$depth * 100)
for (i_r in 1:nrow(snp_positions)) {
  vec_gini <-
    c(
      snp_positions[i_r,]$acount,
      snp_positions[i_r,]$ccount,
      snp_positions[i_r,]$gcount,
      snp_positions[i_r,]$tcount,
      snp_positions[i_r,]$ncount,
      snp_positions[i_r,]$indelcount
    )
  snp_positions[i_r,]$gini <- Gini(vec_gini, unbiased = FALSE)
}

write.csv(x = snp_positions, file = "/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/locate_snps/variants_data/snp_positions_100.csv")
saveRDS(snp_positions, file = "/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/locate_snps/variants_data/snp_positions_100.rds")

snp_positions <-
  df_all_rows %>% filter(pos %in% for_plot$REFERENCE_POSITION) %>% filter(depth >= 250) # %>% filter(sample_method == "unblasted")
snp_positions$gini <- as.double(NA)
snp_positions$fraction_reference <-
  as.double(snp_positions$refcount / snp_positions$depth * 100)
snp_positions$fraction_alternate <-
  as.double(snp_positions$altcount / snp_positions$depth * 100)
snp_positions$fraction_uncalled <-
  as.double(snp_positions$ncount / snp_positions$depth * 100)
snp_positions$fraction_indel <-
  as.double(snp_positions$indelcount / snp_positions$depth * 100)
for (i_r in 1:nrow(snp_positions)) {
  vec_gini <-
    c(
      snp_positions[i_r,]$acount,
      snp_positions[i_r,]$ccount,
      snp_positions[i_r,]$gcount,
      snp_positions[i_r,]$tcount,
      snp_positions[i_r,]$ncount,
      snp_positions[i_r,]$indelcount
    )
  snp_positions[i_r,]$gini <- Gini(vec_gini, unbiased = FALSE)
}

write.csv(x = snp_positions, file = "/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/locate_snps/variants_data/snp_positions_250.csv")
saveRDS(snp_positions, file = "/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/locate_snps/variants_data/snp_positions_250.rds")

snp_positions <-
  df_all_rows %>% filter(pos %in% for_plot$REFERENCE_POSITION) %>% filter(depth >= 500) # %>% filter(sample_method == "unblasted")
snp_positions$gini <- as.double(NA)
snp_positions$fraction_reference <-
  as.double(snp_positions$refcount / snp_positions$depth * 100)
snp_positions$fraction_alternate <-
  as.double(snp_positions$altcount / snp_positions$depth * 100)
snp_positions$fraction_uncalled <-
  as.double(snp_positions$ncount / snp_positions$depth * 100)
snp_positions$fraction_indel <-
  as.double(snp_positions$indelcount / snp_positions$depth * 100)
for (i_r in 1:nrow(snp_positions)) {
  vec_gini <-
    c(
      snp_positions[i_r,]$acount,
      snp_positions[i_r,]$ccount,
      snp_positions[i_r,]$gcount,
      snp_positions[i_r,]$tcount,
      snp_positions[i_r,]$ncount,
      snp_positions[i_r,]$indelcount
    )
  snp_positions[i_r,]$gini <- Gini(vec_gini, unbiased = FALSE)
}

write.csv(x = snp_positions, file = "/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/locate_snps/variants_data/snp_positions_500.csv")
saveRDS(snp_positions, file = "/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/locate_snps/variants_data/snp_positions_500.rds")

snp_positions <-
  df_all_rows %>% filter(pos %in% for_plot$REFERENCE_POSITION) %>% filter(depth >= 1000) # %>% filter(sample_method == "unblasted")
snp_positions$gini <- as.double(NA)
snp_positions$fraction_reference <-
  as.double(snp_positions$refcount / snp_positions$depth * 100)
snp_positions$fraction_alternate <-
  as.double(snp_positions$altcount / snp_positions$depth * 100)
snp_positions$fraction_uncalled <-
  as.double(snp_positions$ncount / snp_positions$depth * 100)
snp_positions$fraction_indel <-
  as.double(snp_positions$indelcount / snp_positions$depth * 100)
for (i_r in 1:nrow(snp_positions)) {
  vec_gini <-
    c(
      snp_positions[i_r,]$acount,
      snp_positions[i_r,]$ccount,
      snp_positions[i_r,]$gcount,
      snp_positions[i_r,]$tcount,
      snp_positions[i_r,]$ncount,
      snp_positions[i_r,]$indelcount
    )
  snp_positions[i_r,]$gini <- Gini(vec_gini, unbiased = FALSE)
}

write.csv(x = snp_positions, file = "/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/locate_snps/variants_data/snp_positions_1000.csv")
saveRDS(snp_positions, file = "/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/locate_snps/variants_data/snp_positions_1000.rds")

snp_positions <-
  df_all_rows %>% filter(pos %in% for_plot$REFERENCE_POSITION) %>% filter(depth >= 5000) # %>% filter(sample_method == "unblasted")
snp_positions$gini <- as.double(NA)
snp_positions$fraction_reference <-
  as.double(snp_positions$refcount / snp_positions$depth * 100)
snp_positions$fraction_alternate <-
  as.double(snp_positions$altcount / snp_positions$depth * 100)
snp_positions$fraction_uncalled <-
  as.double(snp_positions$ncount / snp_positions$depth * 100)
snp_positions$fraction_indel <-
  as.double(snp_positions$indelcount / snp_positions$depth * 100)
for (i_r in 1:nrow(snp_positions)) {
  vec_gini <-
    c(
      snp_positions[i_r,]$acount,
      snp_positions[i_r,]$ccount,
      snp_positions[i_r,]$gcount,
      snp_positions[i_r,]$tcount,
      snp_positions[i_r,]$ncount,
      snp_positions[i_r,]$indelcount
    )
  snp_positions[i_r,]$gini <- Gini(vec_gini, unbiased = FALSE)
}

write.csv(x = snp_positions, file = "/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/locate_snps/variants_data/snp_positions_5000.csv")
saveRDS(snp_positions, file = "/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/locate_snps/variants_data/snp_positions_5000.rds")

snp_positions <-
  df_all_rows %>% filter(pos %in% for_plot$REFERENCE_POSITION) %>% filter(depth >= 10000) # %>% filter(sample_method == "unblasted")
snp_positions$gini <- as.double(NA)
snp_positions$fraction_reference <-
  as.double(snp_positions$refcount / snp_positions$depth * 100)
snp_positions$fraction_alternate <-
  as.double(snp_positions$altcount / snp_positions$depth * 100)
snp_positions$fraction_uncalled <-
  as.double(snp_positions$ncount / snp_positions$depth * 100)
snp_positions$fraction_indel <-
  as.double(snp_positions$indelcount / snp_positions$depth * 100)

snp_positions$fraction_a <-
  as.double(snp_positions$acount / snp_positions$depth * 100)
snp_positions$fraction_c <-
  as.double(snp_positions$ccount / snp_positions$depth * 100)
snp_positions$fraction_g <-
  as.double(snp_positions$gcount / snp_positions$depth * 100)
snp_positions$fraction_t <-
  as.double(snp_positions$tcount / snp_positions$depth * 100)

for (i_r in 1:nrow(snp_positions)) {
  vec_gini <-
    c(
      snp_positions[i_r,]$acount,
      snp_positions[i_r,]$ccount,
      snp_positions[i_r,]$gcount,
      snp_positions[i_r,]$tcount,
      snp_positions[i_r,]$ncount,
      snp_positions[i_r,]$indelcount
    )
  snp_positions[i_r,]$gini <- Gini(vec_gini, unbiased = FALSE)
}

write.csv(x = snp_positions, file = "/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/locate_snps/variants_data/snp_positions_10000.csv")
saveRDS(snp_positions, file = "/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/locate_snps/variants_data/snp_positions_10000.rds")

snp_positions <-
  snp_positions %>% arrange(sample_name,
                            sample_method,
                            desc(fraction_alternate),
                            pos,
                            gini)
write.csv(x = snp_positions, file = "/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/locate_snps/variants_data/snp_positions_10000_sorted_bases.csv")
saveRDS(snp_positions, file = "/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/locate_snps/variants_data/snp_positions_10000_sorted_bases.rds")


input_dir <-
  "/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data_raw/counts_controls/"
input_list <-
  list.files(input_dir, full.names = TRUE,  pattern = '.tsv$')
files_to_read <-
  list.files(path = input_dir,
             pattern = "\\.tsv$",
             full.names = TRUE)
all_files <- lapply(files_to_read, function(x) {
  read.table(file = x,
             sep = '\t',
             header = TRUE)
})

df_all_rows <- bind_rows(all_files)


foreach_snp_positions <-
  df_all_rows %>% filter(pos %in% for_plot$REFERENCE_POSITION) %>% filter(sample_method == "unblasted") %>% filter(depth >= 500)
d <- data.frame(x = 1:10, y = rnorm(10))
s <- foreach(d = iter(d, by = 'row'), .combine = rbind) %dopar% d

detectCores()
registerDoParallel(6)

start <- proc.time()
r <- foreach(icount(trials), .combine = rbind) %dopar% {
  ind <- sample(100, 100, replace = TRUE)
  result1 <- glm(x[ind, 2] ~ x[ind, 1], family = binomial(logit))
  coefficients(result1)
}
dopar_loop <- proc.time() - start

# Option to save shapefiles or not

# SNPs_EPI_ISL_710528 <- MN908947_3_all_snp %>% filter(QUERY != ".", QUERY_TAG == "EPI_ISL_710528") %>% arrange(REFERENCE_POSITION)
# SNPs_EPI_ISL_6841980 <- MN908947_3_all_snp %>% filter(QUERY != ".", QUERY_TAG == "EPI_ISL_6841980") %>% arrange(REFERENCE_POSITION)
# SNPs_MT007544_1 <- MN908947_3_all_snp %>% filter(QUERY != ".", QUERY_TAG == "MT007544.1") %>% arrange(REFERENCE_POSITION)
#
# SNPs_delta <- MN908947_3_all_snp %>% filter(QUERY != ".", QUERY_TAG == "EPI_ISL_2693246") %>% arrange(REFERENCE_POSITION)
# SNPs_alpha <- MN908947_3_all_snp %>% filter(QUERY != ".", QUERY_TAG == "EPI_ISL_710528") %>% arrange(REFERENCE_POSITION)
# SNPs_omicron <- MN908947_3_all_snp %>% filter(QUERY != ".", QUERY_TAG == "EPI_ISL_6841980") %>% arrange(REFERENCE_POSITION)
# SNPs_MT007544_1 <- MN908947_3_all_snp %>% filter(QUERY != ".", QUERY_TAG == "MT007544.1") %>% arrange(REFERENCE_POSITION)
