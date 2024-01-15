library(tidyverse)
library(readr)
library(ggplot2)
library(stringr)
library(DescTools)
library(splitstackshape)
library(foreach)
library(doParallel)
library(lubridate)
library(vcfR)
# library(VariantAnnotation)
library(genpwr)
library(Rsamtools)
library(pegas)
library(vroom)
library(ggvenn)
library(ggVennDiagram)
library(VennDiagram)

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

file_all_collection <- paste0(dir_data_generated, "/", "collection_date.csv") 
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

# "/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data_develo/20076a092_01_combined_from_vcf.rds"
# "/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data_develo/20076a092_01_ivar_from_vcf.rds"
# "/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data_develo/20076a092_01_lofreq_from_vcf.rds"
# "/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data_develo/20076a092_01_umivar_from_vcf.rds"
# "/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data_develo/20076a092_01_umivar_hq_from_vcf.rds"
# "/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data_develo/20076a092_01_varscan_from_vcf.rds"

combined_from_vcf.rds <- readRDS("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data_develo/20076a092_01_combined_from_vcf.rds")
readRDS("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data_develo/20076a092_01_ivar_from_vcf.rds")
readRDS("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data_develo/20076a092_01_lofreq_from_vcf.rds")
readRDS("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data_develo/20076a092_01_umivar_from_vcf.rds")
readRDS("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data_develo/20076a092_01_umivar_hq_from_vcf.rds")
readRDS("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data_develo/20076a092_01_varscan_from_vcf.rds")

"/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data_develo/20076a092_01_bowtie2_blasted_sorted_ivar_tsv.vcf"
"/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data_develo/20076a092_01_bowtie2_blasted_sorted_lofreq_corrected.vcf"
"/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data_develo/20076a092_01_bowtie2_blasted_sorted_varscan_corrected.vcf"
"/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data_develo/20076a092_01_bowtie2_unblasted_sorted_corrected_sorted_deduplicated_sorted_hq_umivar.vcf"
"/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data_develo/20076a092_01_bowtie2_unblasted_sorted_corrected_sorted_deduplicated_sorted_umivar.vcf"

ivar_tsv.vcf <- read.vcfR("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data_develo/20076a092_01_bowtie2_blasted_sorted_ivar_tsv.vcf")
lofreq_corrected.vcf <- read.vcfR("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data_develo/20076a092_01_bowtie2_blasted_sorted_lofreq_corrected.vcf")
varscan_corrected.vcf <- read.vcfR("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data_develo/20076a092_01_bowtie2_blasted_sorted_varscan_corrected.vcf")
hq_umivar.vcf <- read.vcfR("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data_develo/20076a092_01_bowtie2_unblasted_sorted_corrected_sorted_deduplicated_sorted_hq_umivar.vcf")
umivar.vcf <- read.vcfR("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data_develo/20076a092_01_bowtie2_unblasted_sorted_corrected_sorted_deduplicated_sorted_umivar.vcf")

ivar_tsv.tidy <- read.vcfR("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data_develo/20076a092_01_bowtie2_blasted_sorted_ivar_tsv.vcf") %>% vcfR2tidy()
lofreq_corrected.tidy <- read.vcfR("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data_develo/20076a092_01_bowtie2_blasted_sorted_lofreq_corrected.vcf") %>% vcfR2tidy()
varscan_corrected.tidy <- read.vcfR("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data_develo/20076a092_01_bowtie2_blasted_sorted_varscan_corrected.vcf") %>% vcfR2tidy()
hq_umivar.tidy <- read.vcfR("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data_develo/20076a092_01_bowtie2_unblasted_sorted_corrected_sorted_deduplicated_sorted_hq_umivar.vcf") %>% vcfR2tidy()
umivar.tidy <- read.vcfR("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data_develo/20076a092_01_bowtie2_unblasted_sorted_corrected_sorted_deduplicated_sorted_umivar.vcf") %>% vcfR2tidy()

ivar_tsv.tidy$gt
ivar_tsv.tidy$meta
ivar_tsv.tidy$fix
ivar_tsv.out <- cbind(ivar_tsv.tidy$fix, ivar_tsv.tidy$gt %>% select(-ChromKey, -POS)) %>% arrange(POS) %>% tibble()

lofreq_corrected.tidy$gt
lofreq_corrected.tidy$meta
lofreq_corrected.tidy$fix
lofreq_corrected.out <- cbind(lofreq_corrected.tidy$fix, lofreq_corrected.tidy$gt %>% select(-ChromKey, -POS)) %>% arrange(POS) %>% tibble()

varscan_corrected.tidy$gt
varscan_corrected.tidy$meta
varscan_corrected.tidy$fix
varscan_corrected.out <- cbind(varscan_corrected.tidy$fix, varscan_corrected.tidy$gt %>% select(-ChromKey, -POS)) %>% arrange(POS) %>% tibble()

hq_umivar.tidy$gt
hq_umivar.tidy$meta
hq_umivar.tidy$fix
hq_umivar.out <- cbind(hq_umivar.tidy$fix, hq_umivar.tidy$gt %>% select(-ChromKey, -POS)) %>% arrange(POS) %>% tibble()

bind_rows(hq_umivar.out, varscan_corrected.out)


combined_from_vcf.rds %>% filter(DP >= 100) %>% arrange(POS) %>% tibble()

combined_from_vcf.rds$THESIS_DP_TOTAL <- rep(0, nrow(combined_from_vcf.rds))
combined_from_vcf.rds$THESIS_DP_REF <- rep(0, nrow(combined_from_vcf.rds))
combined_from_vcf.rds$THESIS_DP_ALT <- rep(0, nrow(combined_from_vcf.rds))

combined_from_vcf.rds$THESIS_FREQ_TOTAL <- rep(0, nrow(combined_from_vcf.rds))
combined_from_vcf.rds$THESIS_FREQ_REF <- rep(0, nrow(combined_from_vcf.rds))
combined_from_vcf.rds$THESIS_FREQ_ALT <- rep(0, nrow(combined_from_vcf.rds))

combined_from_vcf.rds$THESIS_QUAL_TOTAL <- rep(0, nrow(combined_from_vcf.rds))
combined_from_vcf.rds$THESIS_QUAL_REF <- rep(0, nrow(combined_from_vcf.rds))
combined_from_vcf.rds$THESIS_QUAL_ALT <- rep(0, nrow(combined_from_vcf.rds))



str(combined_from_vcf.rds)


combined_from_vcf.rds

combined_from_vcf.rds$EXP_DP[combined_from_vcf.rds$CALLER == 'ivar']
combined_from_vcf.rds[CALLER == 'ivar']


combined_from_vcf.rds
combined_from_vcf.rds <- combined_from_vcf.rds %>% filter(FILTER == "PASS")
combined_from_vcf.rds <- combined_from_vcf.rds %>% filter(FILTER == "PASS")

combined_from_vcf.rds$THESIS_BASE_CHANGE <- paste(combined_from_vcf.rds$REF, combined_from_vcf.rds$ALT, sep = ">")
combined_from_vcf.rds$THESIS_VARIANT <- paste(combined_from_vcf.rds$POS, combined_from_vcf.rds$THESIS_BASE_CHANGE, sep = ":")

combined_from_vcf.rds$VARIANT <- paste(combined_from_vcf.rds$POS, combined_from_vcf.rds$REF, combined_from_vcf.rds$ALT, sep = ":")

#ALT_FREQ
combined_from_vcf.rds[combined_from_vcf.rds$CALLER == 'ivar',]$gt_ALT_FREQ %>% parse_number()
combined_from_vcf.rds[combined_from_vcf.rds$CALLER == 'varscan',]$gt_FREQ %>% parse_number()/100
combined_from_vcf.rds[combined_from_vcf.rds$CALLER == 'lofreq',]$AF
combined_from_vcf.rds[combined_from_vcf.rds$CALLER == 'umivar',]$gt_AF
combined_from_vcf.rds[combined_from_vcf.rds$CALLER == 'umivar_hq',]$gt_AF

#REF_FREQ
1 - (combined_from_vcf.rds[combined_from_vcf.rds$CALLER == 'ivar',]$gt_ALT_FREQ %>% parse_number())
1 - (combined_from_vcf.rds[combined_from_vcf.rds$CALLER == 'varscan',]$gt_FREQ %>% parse_number()/100)
1 - (combined_from_vcf.rds[combined_from_vcf.rds$CALLER == 'lofreq',]$AF)
1 - (combined_from_vcf.rds[combined_from_vcf.rds$CALLER == 'umivar',]$gt_AF)
1 - (combined_from_vcf.rds[combined_from_vcf.rds$CALLER == 'umivar_hq',]$gt_AF)

#VARIANTS_PASS
combined_from_vcf.rds[combined_from_vcf.rds$CALLER == 'ivar',]$POS %>% as.numeric()
combined_from_vcf.rds[combined_from_vcf.rds$CALLER == 'varscan',]$POS %>% as.numeric()
combined_from_vcf.rds[combined_from_vcf.rds$CALLER == 'lofreq',]$POS %>% as.numeric()
combined_from_vcf.rds[combined_from_vcf.rds$CALLER == 'umivar',]$POS %>% as.numeric()
combined_from_vcf.rds[combined_from_vcf.rds$CALLER == 'umivar_hq',]$POS %>% as.numeric()

#VARIANTS_PASS
combined_from_vcf.rds[combined_from_vcf.rds$CALLER == 'ivar',]$POS %>% length()
combined_from_vcf.rds[combined_from_vcf.rds$CALLER == 'varscan',]$POS %>% length()
combined_from_vcf.rds[combined_from_vcf.rds$CALLER == 'lofreq',]$POS %>% length()
combined_from_vcf.rds[combined_from_vcf.rds$CALLER == 'umivar',]$POS %>% length()
combined_from_vcf.rds[combined_from_vcf.rds$CALLER == 'umivar_hq',]$POS %>% length()

#VARIANTS_PASS_UNIQUE_POS
combined_from_vcf.rds[combined_from_vcf.rds$CALLER == 'ivar',]$POS %>% unique() %>% length()
combined_from_vcf.rds[combined_from_vcf.rds$CALLER == 'varscan',]$POS %>% unique() %>% length()
combined_from_vcf.rds[combined_from_vcf.rds$CALLER == 'lofreq',]$POS %>% unique() %>% length()
combined_from_vcf.rds[combined_from_vcf.rds$CALLER == 'umivar',]$POS %>% unique() %>% length()
combined_from_vcf.rds[combined_from_vcf.rds$CALLER == 'umivar_hq',]$POS %>% unique() %>% length()

#VARIANTS_PASS_MULTIPLE_VARIANTS_UNIQUE_POS
(combined_from_vcf.rds[combined_from_vcf.rds$CALLER == 'ivar',]$POS %>% length()) - (combined_from_vcf.rds[combined_from_vcf.rds$CALLER == 'ivar',]$POS %>% unique() %>% length())
(combined_from_vcf.rds[combined_from_vcf.rds$CALLER == 'varscan',]$POS %>% length()) - (combined_from_vcf.rds[combined_from_vcf.rds$CALLER == 'varscan',]$POS %>% unique() %>% length())
(combined_from_vcf.rds[combined_from_vcf.rds$CALLER == 'lofreq',]$POS %>% length()) - (combined_from_vcf.rds[combined_from_vcf.rds$CALLER == 'lofreq',]$POS %>% unique() %>% length())
(combined_from_vcf.rds[combined_from_vcf.rds$CALLER == 'umivar',]$POS %>% length()) - (combined_from_vcf.rds[combined_from_vcf.rds$CALLER == 'umivar',]$POS %>% unique() %>% length())
(combined_from_vcf.rds[combined_from_vcf.rds$CALLER == 'umivar_hq',]$POS%>% length()) - (combined_from_vcf.rds[combined_from_vcf.rds$CALLER == 'umivar_hq',]$POS %>% unique() %>% length())

combined_from_vcf.rds

pos_list_ivar <- combined_from_vcf.rds[combined_from_vcf.rds$CALLER == 'ivar',]$POS %>% as.numeric()
pos_list_varscan <- combined_from_vcf.rds[combined_from_vcf.rds$CALLER == 'varscan',]$POS %>% as.numeric()
pos_list_lofreq <- combined_from_vcf.rds[combined_from_vcf.rds$CALLER == 'lofreq',]$POS %>% as.numeric()
pos_list_umivar <- combined_from_vcf.rds[combined_from_vcf.rds$CALLER == 'umivar',]$POS %>% as.numeric()
pos_list_umivar_hq <- combined_from_vcf.rds[combined_from_vcf.rds$CALLER == 'umivar_hq',]$POS %>% as.numeric()

pos_list_ivar <- combined_from_vcf.rds[combined_from_vcf.rds$CALLER == 'ivar',]$THESIS_VARIANT %>% trimws() %>% toupper()
pos_list_varscan <- combined_from_vcf.rds[combined_from_vcf.rds$CALLER == 'varscan',]$THESIS_VARIANT %>% trimws() %>% toupper()
pos_list_lofreq <- combined_from_vcf.rds[combined_from_vcf.rds$CALLER == 'lofreq',]$THESIS_VARIANT %>% trimws() %>% toupper()
pos_list_umivar <- combined_from_vcf.rds[combined_from_vcf.rds$CALLER == 'umivar',]$THESIS_VARIANT %>% trimws() %>% toupper()
pos_list_umivar_hq <- combined_from_vcf.rds[combined_from_vcf.rds$CALLER == 'umivar_hq',]$THESIS_VARIANT %>% trimws() %>% toupper()

x <- list(
  ivar = pos_list_ivar, 
  varscan = pos_list_varscan, 
  lofreq = pos_list_lofreq,
  umivar = pos_list_umivar_hq
)

library(ggvenn)
ggvenn(
  x, 
  fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"),
  stroke_size = 0.5, set_name_size = 4
)

library("ggVennDiagram")
ggVennDiagram(x, label_alpha = 0)

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

# counts_experimental
input_dir <-
  "/mnt/bulk_data/THESIS/DATA_THESIS/GILOT_DATA_THESIS_CONTROL/rds"
input_list_ivar <-
  list.files(input_dir, full.names = TRUE,  pattern = '_ivar_from_vcf.rds$')
input_list_lofreq <-
  list.files(input_dir, full.names = TRUE,  pattern = '_lofreq_from_vcf.rds$')
input_list_varscan <-
  list.files(input_dir, full.names = TRUE,  pattern = '_varscan_from_vcf.rds$')
input_list_umivar <-
  list.files(input_dir, full.names = TRUE,  pattern = '_umivar_from_vcf.rds$')
input_list_umivar_hq <-
  list.files(input_dir, full.names = TRUE,  pattern = '_umivar_hq_from_vcf.rds$')

all_files_ivar <- do.call('rbind', lapply(input_list_ivar, readRDS))
all_files_lofreq <- do.call('rbind', lapply(input_list_lofreq, readRDS))
all_files_varscan <- do.call('rbind', lapply(input_list_varscan, readRDS))
all_files_umivar <- do.call('rbind', lapply(input_list_umivar, readRDS))
all_files_umivar_hq <- do.call('rbind', lapply(input_list_umivar_hq, readRDS))





all_files_ivar <- bind_rows(readRDS(input_list_ivar))

readRDS(input_list_ivar[2])

all_files_ivar <- do.call('rbind', lapply(input_list_ivar, readRDS))

all_files_ivar <- lapply(input_list_ivar, function(x) {
  bind_rows(x)
})


files_to_read <-
  list.files(path = input_dir,
             pattern = "\\.tsv$",
             full.names = TRUE)
all_files <- lapply(files_to_read, function(x) {
  read.table(file = x,
             sep = '\t',
             header = TRUE)
})
