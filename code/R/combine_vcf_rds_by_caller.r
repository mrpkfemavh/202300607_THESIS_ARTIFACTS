# library(DescTools)
# library(doParallel)
# library(foreach)
# library(genpwr)
# library(ggplot2)
# library(ggvenn)
# library(ggVennDiagram)
# library(lubridate)
# library(pegas)
# library(readr)
# library(Rsamtools)
# library(splitstackshape)
# library(stringr)
# library(VariantAnnotation)
# library(vcfR)
# library(VennDiagram)
# library(vroom)

rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() #free up memrory and report the memory usage.

library(tidyverse)
library(vroom)
library(lubridate)

github_data_unified <- "/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/"

##########################
##########################      SIGNATURES_OF_CONTROLS
##########################
########################################################################################

SNPs_nucmer_MT007544_1_australia <- readRDS("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_computed/SNPs_MT007544_1_australia.rds")
SNPs_nucmer_EPI_ISL_omicron <- readRDS("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_computed/SNPs_EPI_ISL_omicron.rds")
SNPs_nucmer_EPI_ISL_2693246_delta <- readRDS("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_computed/SNPs_EPI_ISL_2693246_delta.rds")
SNPs_nucmer_EPI_ISL_710528_alpha <- readRDS("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_computed/SNPs_EPI_ISL_710528_alpha.rds")

SNPs_nucmer_from_synthetic_controls <- bind_rows(SNPs_nucmer_MT007544_1_australia, SNPs_nucmer_EPI_ISL_omicron, SNPs_nucmer_EPI_ISL_2693246_delta, SNPs_nucmer_EPI_ISL_710528_alpha) %>% arrange(REFERENCE_POSITION, QUERY_TAG) %>% unique()

SNPs_nucmer_from_synthetic_controls$SIGNATURE_CONTROL_IDENTIFIER <- SNPs_nucmer_from_synthetic_controls$QUERY_TAG %>% str_replace("\\.", "_")
SNPs_nucmer_from_synthetic_controls$SIGNATURE_POS <- SNPs_nucmer_from_synthetic_controls$REFERENCE_POSITION %>% as.integer()
SNPs_nucmer_from_synthetic_controls$SIGNATURE_BASE_CHANGE <- paste(SNPs_nucmer_from_synthetic_controls$REFERENCE, SNPs_nucmer_from_synthetic_controls$QUERY, sep = ">")
SNPs_nucmer_from_synthetic_controls$SIGNATURE_VARIANT <- paste(SNPs_nucmer_from_synthetic_controls$SIGNATURE_BASE_CHANGE, SNPs_nucmer_from_synthetic_controls$SIGNATURE_POS, sep = ":")
SNPs_nucmer_from_synthetic_controls$SIGNATURE_VARIANT_BY_POS <- paste(SNPs_nucmer_from_synthetic_controls$SIGNATURE_POS, SNPs_nucmer_from_synthetic_controls$SIGNATURE_BASE_CHANGE, sep = ":")
SNPs_nucmer_from_synthetic_controls$SIGNATURE_VARIANT_BY_HGVS <- paste0(SNPs_nucmer_from_synthetic_controls$SIGNATURE_CONTROL_IDENTIFIER, ".", SNPs_nucmer_from_synthetic_controls$SIGNATURE_POS, SNPs_nucmer_from_synthetic_controls$SIGNATURE_BASE_CHANGE)
colnames(SNPs_nucmer_from_synthetic_controls)

data_fwrite <- SNPs_nucmer_from_synthetic_controls
name_data_frwite <- "SNPs_nucmer_from_synthetic_controls"
new_filename <- paste0(github_data_unified, name_data_frwite, ".csv")
print(new_filename)
data.table::fwrite(data_fwrite, file = new_filename, append = FALSE, quote = "auto", showProgress = getOption("datatable.showProgress", interactive()), buffMB = 1024, nThread = 6)

new_filename <- paste0(github_data_unified, name_data_frwite, ".rds")
print(new_filename)
saveRDS(data_fwrite, file = new_filename, compress = TRUE)


##########################
##########################      Location Extra Info
##########################
########################################################################################
file_all_snp_data_header <- "/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_generated/column.header.csv"
file_all_metadata <- "/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_generated/metadata.csv" # "/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_computed/column.header.csv"

file_MN908947_3_fasta <- "/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_raw/reference_variants/MN908947.3.fasta"
file_MN908947_3_gff <- "/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_raw/reference_variants/MN908947.3.gff3"
file_MN908947_3_gtf <- "/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_raw/reference_variants/Sars_cov_2.ASM985889v3.101.gtf"

# dna_MN908947_3 <- ape::read.dna(file_MN908947_3_fasta, format = "fasta")
# gff_MN908947_3 <- read.table(file_MN908947_3_gff, sep="\t", quote="") %>% tibble()

##########################
##########################      Features GTF
##########################
########################################################################################
gtf_MN908947_3 <-
  read_delim(
    file_MN908947_3_gtf, # "/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/locate_snps/variants_data/MN908947.3.gff3", #"Downloads/SarsCovControls/MN908947.3.gff3",
    "\t",
    escape_double = FALSE,
    col_names = FALSE,
    trim_ws = TRUE,
    comment = "#"
    # skip = 4
  ) %>% tibble()

# colnames(gtf_MN908947_3)
names(gtf_MN908947_3) <-
  c("seqid",
    "source",
    "type",
    "start",
    "end",
    "score",
    "strand",
    "phase",
    "attributes")
# colnames(gtf_MN908947_3)

gtf_MN908947_3 <- cSplit(gtf_MN908947_3, "attributes", ";", drop = TRUE, stripWhite = TRUE)
gtf_MN908947_3 <- gtf_MN908947_3 %>% arrange(start, end)

data_unified <- gtf_MN908947_3
file_data_unified <- "gtf_MN908947_3"
dir_data_unified <- "/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified"
write.csv(data_unified, paste0(dir_data_unified, "/", file_data_unified, ".csv"), row.names=FALSE, fileEncoding = "UTF-8")
saveRDS(data_unified, paste0(dir_data_unified, "/", file_data_unified, ".rds"), compress = TRUE)
infoRDS(paste0(dir_data_unified, "/", file_data_unified, ".rds"))
# readRDS(paste0(dir_data_unified, "/", file_data_unified, ".rds"))

##########################
##########################      Features GFF
##########################
########################################################################################
gff_MN908947_3 <-
  read_delim(
    file_MN908947_3_gff, # "/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/locate_snps/variants_data/MN908947.3.gff3", #"Downloads/SarsCovControls/MN908947.3.gff3",
    "\t",
    escape_double = FALSE,
    col_names = FALSE,
    trim_ws = TRUE,
    comment = "#"
    # skip = 2
  ) %>% tibble()

colnames(gff_MN908947_3)

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

colnames(gff_MN908947_3)

gff_MN908947_3 <- cSplit(gff_MN908947_3, "attributes", ";", drop = TRUE, stripWhite = TRUE)
gff_MN908947_3 <- gff_MN908947_3 %>% arrange(start, end)

# gff_MN908947_3 %>% filter(!(type %in% c("region", "gene"))) %>% arrange(start, end) %>% view()
# gtf_MN908947_3 %>% filter(!(type %in% c("region", "gene", "exon", "transcript"))) %>% arrange(start, end) %>% view()

gff_MN908947_3 %>% filter(!(type %in% c("region", "gene"))) %>% arrange(start, end) %>% tibble()
gff_join <- gff_MN908947_3 %>% filter(!(type %in% c("region", "gene"))) %>% arrange(start, end) %>% tibble()

colnames(gff_join) <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "cds", "gene", "dbxref", "name", "note", "translation", "genomic", "product", "attributes_09", "attributes_10", "attributes_11")

gff_join$start <- gff_join$start %>% as.integer() 
gff_join$end <- gff_join$end %>% as.integer() 

gff_join$cds <- gff_join$cds %>% as.character()
gff_join$cds <- (sapply(strsplit(gff_join$cds,"="), `[`, 2))

gff_join$gene <- gff_join$gene %>% as.character()
gff_join$gene <- (sapply(strsplit(gff_join$gene,"="), `[`, 2))

gff_join$dbxref <- gff_join$dbxref %>% as.character()
gff_join$dbxref <- (sapply(strsplit(gff_join$dbxref,"="), `[`, 2))

gff_join$name <- gff_join$name %>% as.character()
gff_join$name <- (sapply(strsplit(gff_join$name,"="), `[`, 2))

gff_join$note <- gff_join$note %>% as.character()
gff_join$note <- (sapply(strsplit(gff_join$note,"="), `[`, 2))

gff_join$translation <- gff_join$translation %>% as.character()
gff_join$translation <- (sapply(strsplit(gff_join$translation,"="), `[`, 2))

gff_join$genomic <- gff_join$genomic %>% as.character()
gff_join$genomic <- (sapply(strsplit(gff_join$genomic,"="), `[`, 2))

gff_join$product <- gff_join$product %>% as.character()
gff_join$product <- (sapply(strsplit(gff_join$product,"="), `[`, 2))

gff_join$THESIS_FEATURE_LABEL <- (sapply(strsplit(gff_join$gene,"-"), `[`, 2))
gff_join$THESIS_FEATURE_PRODUCT <- gff_join$product
gff_join$THESIS_FEATURE_GENE <- gff_join$gene
gff_join$THESIS_FEATURE_GENOMIC <- gff_join$genomic
gff_join$THESIS_FEATURE_NCBI <- gff_join$dbxref
gff_join$THESIS_FEATURE_ID <- gff_join$cds

gff_join[gff_join$type == "five_prime_UTR",]$THESIS_FEATURE_LABEL <- "five_prime_UTR"
gff_join[gff_join$type == "three_prime_UTR",]$THESIS_FEATURE_LABEL <- "three_prime_UTR"

gff_join[gff_join$type == "five_prime_UTR",]$THESIS_FEATURE_PRODUCT <- "five_prime_UTR"
gff_join[gff_join$type == "three_prime_UTR",]$THESIS_FEATURE_PRODUCT <- "three_prime_UTR"

gff_join$THESIS_FEATURE_LABEL
gff_join$THESIS_FEATURE_PRODUCT
gff_join$THESIS_FEATURE_GENE
gff_join$THESIS_FEATURE_NCBI

print("gff file")

data_unified <- gff_join
file_data_unified <- "gff_MN908947_3"
dir_data_unified <- "/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified"
saveRDS(data_unified, paste0(dir_data_unified, "/", file_data_unified, ".rds"), compress = TRUE)
infoRDS(paste0(dir_data_unified, "/", file_data_unified, ".rds"))

data_fwrite <- data_unified
name_data_frwite <- "gff_MN908947_3"
new_filename <- paste0(github_data_unified, name_data_frwite, ".csv")
print(new_filename)
data.table::fwrite(data_fwrite, file = new_filename, append = FALSE, quote = "auto", showProgress = getOption("datatable.showProgress", interactive()), buffMB = 1024, nThread = 6)
#write.csv(data_unified, paste0(dir_data_unified, "/", file_data_unified, ".csv"), row.names=FALSE, fileEncoding = "UTF-8")
# readRDS(paste0(dir_data_unified, "/", file_data_unified, ".rds"))

##########################
##########################      LOAD VCFs
##########################
########################################################################################

input_dir <-
  "/mnt/bulk_data/THESIS/DATA_THESIS/GILOT_DATA_THESIS_EXPERIMENTAL/rds"
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

all_vcf_experimental_ivar <- do.call('rbind', lapply(input_list_ivar, readRDS))
all_vcf_experimental_lofreq <- do.call('rbind', lapply(input_list_lofreq, readRDS))
all_vcf_experimental_varscan <- do.call('rbind', lapply(input_list_varscan, readRDS))
all_vcf_experimental_umivar <- do.call('rbind', lapply(input_list_umivar, readRDS))
all_vcf_experimental_umivar_hq <- do.call('rbind', lapply(input_list_umivar_hq, readRDS))

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

all_vcf_control_ivar <- do.call('rbind', lapply(input_list_ivar, readRDS))
all_vcf_control_lofreq <- do.call('rbind', lapply(input_list_lofreq, readRDS))
all_vcf_control_varscan <- do.call('rbind', lapply(input_list_varscan, readRDS))
all_vcf_control_umivar <- do.call('rbind', lapply(input_list_umivar, readRDS))
all_vcf_control_umivar_hq <- do.call('rbind', lapply(input_list_umivar_hq, readRDS))


all_vcf_control_ivar
all_vcf_experimental_ivar
(colnames(all_vcf_control_ivar) == colnames(all_vcf_experimental_ivar))
colnames(all_vcf_experimental_ivar)
str(all_vcf_experimental_ivar)

##########################
##########################      IVAR VCFs
##########################
########################################################################################
called <- 'ivar'
temp_vcf_combine <- bind_rows(all_vcf_control_ivar, all_vcf_experimental_ivar) %>% unique %>% arrange(POS)
called <- temp_vcf_combine$CALLER %>% unique()
called

# BASE INFORMATION
temp_vcf_combine$THESIS_POS <- temp_vcf_combine[temp_vcf_combine$CALLER == called,]$POS %>% as.integer()
temp_vcf_combine$THESIS_BASE_CHANGE <- paste(temp_vcf_combine[temp_vcf_combine$CALLER == called,]$REF, temp_vcf_combine[temp_vcf_combine$CALLER == called,]$ALT, sep = ">")
temp_vcf_combine$THESIS_VARIANT <- paste(temp_vcf_combine[temp_vcf_combine$CALLER == called,]$THESIS_BASE_CHANGE, temp_vcf_combine[temp_vcf_combine$CALLER == called,]$POS, sep = ":")
temp_vcf_combine$THESIS_VARIANT_BY_POS <- paste(temp_vcf_combine[temp_vcf_combine$CALLER == called,]$POS, temp_vcf_combine[temp_vcf_combine$CALLER == called,]$THESIS_BASE_CHANGE, sep = ":")
temp_vcf_combine$THESIS_VARIANT_BY_HGVS <- paste0((temp_vcf_combine[temp_vcf_combine$CALLER == called,]$CHROM %>% str_replace("\\.", "_")), ".", temp_vcf_combine[temp_vcf_combine$CALLER == called,]$POS, temp_vcf_combine[temp_vcf_combine$CALLER == called,]$THESIS_BASE_CHANGE)
temp_vcf_combine$THESIS_FLAG_CONSVAR <- as.character(NA)
temp_vcf_combine$THESIS_FLAG_INDEL <- as.character(NA)
temp_vcf_combine$THESIS_FLAG_WILDTYPE <- as.character(NA)

temp_vcf_combine$THESIS_PVALUE_UNCORRECTED <- as.numeric(NA)
temp_vcf_combine$THESIS_QVALUE <- as.numeric(NA)
temp_vcf_combine$THESIS_DEPTH_TOTAL <- temp_vcf_combine[temp_vcf_combine$CALLER == called,]$DP %>% as.integer()

temp_vcf_combine$THESIS_DEPTH_REFERENCE <- (temp_vcf_combine[temp_vcf_combine$CALLER == called,]$gt_REF_DP + temp_vcf_combine[temp_vcf_combine$CALLER == called,]$gt_REF_RV) %>% as.integer()
temp_vcf_combine$THESIS_FREQUENCY_NOT_ALTERNATE <- (1 - (temp_vcf_combine[temp_vcf_combine$CALLER == called,]$gt_ALT_FREQ) %>% parse_number())
temp_vcf_combine$THESIS_QUALITY_REFERENCE <- temp_vcf_combine[temp_vcf_combine$CALLER == called,]$gt_REF_QUAL %>% as.numeric()

temp_vcf_combine$THESIS_DEPTH_ALTERNATE <- (temp_vcf_combine[temp_vcf_combine$CALLER == called,]$gt_ALT_DP + temp_vcf_combine[temp_vcf_combine$CALLER == called,]$gt_ALT_RV) %>% as.integer()
temp_vcf_combine$THESIS_FREQUENCY_ALTERNATE <- temp_vcf_combine[temp_vcf_combine$CALLER == called,]$gt_ALT_FREQ %>% parse_number()
temp_vcf_combine$THESIS_QUALITY_ALTERNATE <- temp_vcf_combine[temp_vcf_combine$CALLER == called,]$gt_ALT_QUAL %>% as.numeric()

temp_vcf_combine$THESIS_QC_DEPTH_TOTAL_DISCREPANCY <- (temp_vcf_combine$THESIS_DEPTH_TOTAL - (temp_vcf_combine$THESIS_DEPTH_REFERENCE + temp_vcf_combine$THESIS_DEPTH_ALTERNATE)) %>% as.integer()
temp_vcf_combine$THESIS_QC_DEPTH_TOTAL_CALC_LESSTHAN <- ((temp_vcf_combine$THESIS_DEPTH_REFERENCE + temp_vcf_combine$THESIS_DEPTH_ALTERNATE) < temp_vcf_combine$THESIS_DEPTH_TOTAL)
temp_vcf_combine$THESIS_QC_DEPTH_TOTAL_CALC_EQUAL <- ((temp_vcf_combine$THESIS_DEPTH_REFERENCE + temp_vcf_combine$THESIS_DEPTH_ALTERNATE) == temp_vcf_combine$THESIS_DEPTH_TOTAL)
temp_vcf_combine$THESIS_QC_CALLER_SPECIFIC_5_PERCENT <- as.logical(NA) # temp_vcf_combine$QUAL

# temp_vcf_combine[temp_vcf_combine$CALLER == called,]$POS %>% length()
# temp_vcf_combine[temp_vcf_combine$CALLER == called,]$POS %>% unique() %>% length()
# (temp_vcf_combine[temp_vcf_combine$CALLER == called,]$POS %>% length()) - (temp_vcf_combine[temp_vcf_combine$CALLER == called,]$POS %>% unique() %>% length())

names(temp_vcf_combine)[names(temp_vcf_combine) == 'IDENTIFIER'] <- 'THESIS_IDENTIFIER'
names(temp_vcf_combine)[names(temp_vcf_combine) == 'COHORT'] <- 'THESIS_COHORT'
names(temp_vcf_combine)[names(temp_vcf_combine) == 'BLASTER'] <- 'THESIS_BLASTER'
names(temp_vcf_combine)[names(temp_vcf_combine) == 'CALLER'] <- 'THESIS_CALLER'
names(temp_vcf_combine)[names(temp_vcf_combine) == 'BASENAME'] <- 'THESIS_BASENAME'

all_ivar <- temp_vcf_combine
colnames(all_ivar)
str(all_ivar %>% filter(FILTER == "PASS"))

data_unified <- all_ivar
file_data_unified <- "all_ivar"
dir_data_unified <- "/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified"
#write.csv(data_unified, paste0(dir_data_unified, "/", file_data_unified, ".csv"), row.names=FALSE, fileEncoding = "UTF-8")
saveRDS(data_unified, paste0(dir_data_unified, "/", file_data_unified, ".rds"), compress = TRUE)
infoRDS(paste0(dir_data_unified, "/", file_data_unified, ".rds"))
# readRDS(paste0(dir_data_unified, "/", file_data_unified, ".rds"))

data_fwrite <- data_unified
name_data_frwite <- "all_ivar"
new_filename <- paste0(github_data_unified, name_data_frwite, ".csv")
print(new_filename)
data.table::fwrite(data_fwrite, file = new_filename, append = FALSE, quote = "auto", showProgress = getOption("datatable.showProgress", interactive()), buffMB = 1024, nThread = 6)


##########################
##########################      LOFREQ VCFs
##########################
########################################################################################

# lofreq
all_vcf_control_lofreq
all_vcf_experimental_lofreq
(colnames(all_vcf_control_lofreq) == colnames(all_vcf_experimental_lofreq))
colnames(all_vcf_experimental_lofreq)
str(all_vcf_experimental_lofreq)

called <- 'lofreq'
temp_vcf_combine <- bind_rows(all_vcf_control_lofreq, all_vcf_experimental_lofreq) %>% unique %>% arrange(POS)
called <- temp_vcf_combine$CALLER %>% unique()
called

# BASE INFORMATION
temp_vcf_combine$THESIS_POS <- temp_vcf_combine[temp_vcf_combine$CALLER == called,]$POS %>% as.integer()
temp_vcf_combine$THESIS_BASE_CHANGE <- paste(temp_vcf_combine[temp_vcf_combine$CALLER == called,]$REF, temp_vcf_combine[temp_vcf_combine$CALLER == called,]$ALT, sep = ">")
temp_vcf_combine$THESIS_VARIANT <- paste(temp_vcf_combine[temp_vcf_combine$CALLER == called,]$THESIS_BASE_CHANGE, temp_vcf_combine[temp_vcf_combine$CALLER == called,]$POS, sep = ":")
temp_vcf_combine$THESIS_VARIANT_BY_POS <- paste(temp_vcf_combine[temp_vcf_combine$CALLER == called,]$POS, temp_vcf_combine[temp_vcf_combine$CALLER == called,]$THESIS_BASE_CHANGE, sep = ":")
temp_vcf_combine$THESIS_VARIANT_BY_HGVS <- paste0((temp_vcf_combine[temp_vcf_combine$CALLER == called,]$CHROM %>% str_replace("\\.", "_")), ".", temp_vcf_combine[temp_vcf_combine$CALLER == called,]$POS, temp_vcf_combine[temp_vcf_combine$CALLER == called,]$THESIS_BASE_CHANGE)
temp_vcf_combine$THESIS_FLAG_CONSVAR <- as.character(NA)
temp_vcf_combine$THESIS_FLAG_INDEL <- temp_vcf_combine[temp_vcf_combine$CALLER == called,]$INDEL %>% as.character()
temp_vcf_combine$THESIS_FLAG_WILDTYPE <- as.character(NA)

temp_vcf_combine$THESIS_PVALUE_UNCORRECTED <- as.numeric(NA)
temp_vcf_combine$THESIS_QVALUE <- as.numeric(NA)
temp_vcf_combine$THESIS_DEPTH_TOTAL <- temp_vcf_combine[temp_vcf_combine$CALLER == called,]$DP %>% as.integer()

# DP4     4      Integer Counts for ref-forward bases, ref-reverse, alt-forward and alt-reverse bases   
(sapply(strsplit(temp_vcf_combine$DP4,","), `[`, 1) %>% as.integer()) # ref-forward bases
(sapply(strsplit(temp_vcf_combine$DP4,","), `[`, 2) %>% as.integer()) # ref-reverse bases
(sapply(strsplit(temp_vcf_combine$DP4,","), `[`, 3) %>% as.integer()) # alt-forward
(sapply(strsplit(temp_vcf_combine$DP4,","), `[`, 4) %>% as.integer()) # alt-reverse

temp_vcf_combine$THESIS_DEPTH_REFERENCE <- (sapply(strsplit(temp_vcf_combine$DP4,","), `[`, 1) %>% as.integer()) + (sapply(strsplit(temp_vcf_combine$DP4,","), `[`, 2) %>% as.integer())
temp_vcf_combine$THESIS_FREQUENCY_NOT_ALTERNATE <- (1 - (temp_vcf_combine[temp_vcf_combine$CALLER == called,]$AF)) %>% as.numeric()
temp_vcf_combine$THESIS_QUALITY_REFERENCE <- as.numeric(NA) # temp_vcf_combine[temp_vcf_combine$CALLER == called,]$gt_REF_QUAL %>% as.numeric()

temp_vcf_combine$THESIS_DEPTH_ALTERNATE <- (sapply(strsplit(temp_vcf_combine$DP4,","), `[`, 3) %>% as.integer()) + (sapply(strsplit(temp_vcf_combine$DP4,","), `[`, 4) %>% as.integer())
temp_vcf_combine$THESIS_FREQUENCY_ALTERNATE <- temp_vcf_combine[temp_vcf_combine$CALLER == called,]$AF %>% as.numeric()
temp_vcf_combine$THESIS_QUALITY_ALTERNATE <- as.numeric(NA) # temp_vcf_combine[temp_vcf_combine$CALLER == called,]$gt_ALT_QUAL %>% as.numeric()

temp_vcf_combine$THESIS_QC_DEPTH_TOTAL_DISCREPANCY <- (temp_vcf_combine$THESIS_DEPTH_TOTAL - (temp_vcf_combine$THESIS_DEPTH_REFERENCE + temp_vcf_combine$THESIS_DEPTH_ALTERNATE)) %>% as.integer()
temp_vcf_combine$THESIS_QC_DEPTH_TOTAL_CALC_LESSTHAN <- ((temp_vcf_combine$THESIS_DEPTH_REFERENCE + temp_vcf_combine$THESIS_DEPTH_ALTERNATE) < temp_vcf_combine$THESIS_DEPTH_TOTAL)
temp_vcf_combine$THESIS_QC_DEPTH_TOTAL_CALC_EQUAL <- ((temp_vcf_combine$THESIS_DEPTH_REFERENCE + temp_vcf_combine$THESIS_DEPTH_ALTERNATE) == temp_vcf_combine$THESIS_DEPTH_TOTAL)
temp_vcf_combine$THESIS_QC_CALLER_SPECIFIC_5_PERCENT <- (temp_vcf_combine[temp_vcf_combine$CALLER == called,]$QUAL >= 13 ) # https://sourceforge.net/p/lofreq/discussion/general/thread/7b713493/

# temp_vcf_combine[temp_vcf_combine$CALLER == called,]$POS %>% length()
# temp_vcf_combine[temp_vcf_combine$CALLER == called,]$POS %>% unique() %>% length()
# (temp_vcf_combine[temp_vcf_combine$CALLER == called,]$POS %>% length()) - (temp_vcf_combine[temp_vcf_combine$CALLER == called,]$POS %>% unique() %>% length())

names(temp_vcf_combine)[names(temp_vcf_combine) == 'IDENTIFIER'] <- 'THESIS_IDENTIFIER'
names(temp_vcf_combine)[names(temp_vcf_combine) == 'COHORT'] <- 'THESIS_COHORT'
names(temp_vcf_combine)[names(temp_vcf_combine) == 'BLASTER'] <- 'THESIS_BLASTER'
names(temp_vcf_combine)[names(temp_vcf_combine) == 'CALLER'] <- 'THESIS_CALLER'

all_lofreq <- temp_vcf_combine
colnames(all_lofreq)
str(all_lofreq %>% filter(FILTER == "PASS"))

data_unified <- all_lofreq
file_data_unified <- "all_lofreq"
dir_data_unified <- "/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified"
# write.csv(data_unified, paste0(dir_data_unified, "/", file_data_unified, ".csv"), row.names=FALSE, fileEncoding = "UTF-8")
saveRDS(data_unified, paste0(dir_data_unified, "/", file_data_unified, ".rds"), compress = TRUE)
infoRDS(paste0(dir_data_unified, "/", file_data_unified, ".rds"))
# readRDS(paste0(dir_data_unified, "/", file_data_unified, ".rds"))

data_fwrite <- data_unified
name_data_frwite <- "all_lofreq"
new_filename <- paste0(github_data_unified, name_data_frwite, ".csv")
print(new_filename)
data.table::fwrite(data_fwrite, file = new_filename, append = FALSE, quote = "auto", showProgress = getOption("datatable.showProgress", interactive()), buffMB = 1024, nThread = 6)

##########################
##########################      VARSCAN VCFs
##########################
########################################################################################

# varscan
all_vcf_control_varscan
all_vcf_experimental_varscan
(colnames(all_vcf_control_varscan) == colnames(all_vcf_experimental_varscan))
colnames(all_vcf_experimental_varscan)
str(all_vcf_experimental_varscan)

called <- 'varscan'
temp_vcf_combine <- bind_rows(all_vcf_control_varscan, all_vcf_experimental_varscan) %>% unique %>% arrange(POS)
called <- temp_vcf_combine$CALLER %>% unique()
called

# BASE INFORMATION
temp_vcf_combine$THESIS_POS <- temp_vcf_combine[temp_vcf_combine$CALLER == called,]$POS %>% as.integer()
temp_vcf_combine$THESIS_BASE_CHANGE <- paste(temp_vcf_combine[temp_vcf_combine$CALLER == called,]$REF, temp_vcf_combine[temp_vcf_combine$CALLER == called,]$ALT, sep = ">")
temp_vcf_combine$THESIS_VARIANT <- paste(temp_vcf_combine[temp_vcf_combine$CALLER == called,]$THESIS_BASE_CHANGE, temp_vcf_combine[temp_vcf_combine$CALLER == called,]$POS, sep = ":")
temp_vcf_combine$THESIS_VARIANT_BY_POS <- paste(temp_vcf_combine[temp_vcf_combine$CALLER == called,]$POS, temp_vcf_combine[temp_vcf_combine$CALLER == called,]$THESIS_BASE_CHANGE, sep = ":")
temp_vcf_combine$THESIS_VARIANT_BY_HGVS <- paste0((temp_vcf_combine[temp_vcf_combine$CALLER == called,]$CHROM %>% str_replace("\\.", "_")), ".", temp_vcf_combine[temp_vcf_combine$CALLER == called,]$POS, temp_vcf_combine[temp_vcf_combine$CALLER == called,]$THESIS_BASE_CHANGE)
temp_vcf_combine$THESIS_FLAG_CONSVAR <- as.character(NA)
temp_vcf_combine$THESIS_FLAG_INDEL <- as.character(NA)
temp_vcf_combine$THESIS_FLAG_WILDTYPE <- as.character(NA)

temp_vcf_combine$THESIS_PVALUE_UNCORRECTED <- temp_vcf_combine[temp_vcf_combine$CALLER == called,]$gt_PVAL %>% as.numeric()
temp_vcf_combine$THESIS_QVALUE <- as.numeric(NA)
temp_vcf_combine$THESIS_DEPTH_TOTAL <- temp_vcf_combine[temp_vcf_combine$CALLER == called,]$ADP %>% as.integer()

temp_vcf_combine$THESIS_DEPTH_REFERENCE <- (temp_vcf_combine[temp_vcf_combine$CALLER == called,]$gt_RDF + temp_vcf_combine[temp_vcf_combine$CALLER == called,]$gt_RDR) %>% as.integer()
temp_vcf_combine$THESIS_FREQUENCY_NOT_ALTERNATE <- (1 - (temp_vcf_combine[temp_vcf_combine$CALLER == called,]$gt_FREQ %>% parse_number()/100))
temp_vcf_combine$THESIS_QUALITY_REFERENCE <- temp_vcf_combine[temp_vcf_combine$CALLER == called,]$gt_RBQ %>% as.numeric()

temp_vcf_combine$THESIS_DEPTH_ALTERNATE <- (temp_vcf_combine[temp_vcf_combine$CALLER == called,]$gt_ADF + temp_vcf_combine[temp_vcf_combine$CALLER == called,]$gt_ADR) %>% as.integer()
temp_vcf_combine$THESIS_FREQUENCY_ALTERNATE <- (temp_vcf_combine[temp_vcf_combine$CALLER == called,]$gt_FREQ %>% parse_number()/100)
temp_vcf_combine$THESIS_QUALITY_ALTERNATE <- temp_vcf_combine[temp_vcf_combine$CALLER == called,]$gt_ABQ %>% as.numeric()

temp_vcf_combine$THESIS_QC_DEPTH_TOTAL_DISCREPANCY <- (temp_vcf_combine$THESIS_DEPTH_TOTAL - (temp_vcf_combine$THESIS_DEPTH_REFERENCE + temp_vcf_combine$THESIS_DEPTH_ALTERNATE)) %>% as.integer()
temp_vcf_combine$THESIS_QC_DEPTH_TOTAL_CALC_LESSTHAN <- ((temp_vcf_combine$THESIS_DEPTH_REFERENCE + temp_vcf_combine$THESIS_DEPTH_ALTERNATE) < temp_vcf_combine$THESIS_DEPTH_TOTAL)
temp_vcf_combine$THESIS_QC_DEPTH_TOTAL_CALC_EQUAL <- ((temp_vcf_combine$THESIS_DEPTH_REFERENCE + temp_vcf_combine$THESIS_DEPTH_ALTERNATE) == temp_vcf_combine$THESIS_DEPTH_TOTAL)
temp_vcf_combine$THESIS_QC_CALLER_SPECIFIC_5_PERCENT <- as.logical(NA)

# temp_vcf_combine[temp_vcf_combine$CALLER == called,]$POS %>% length()
# temp_vcf_combine[temp_vcf_combine$CALLER == called,]$POS %>% unique() %>% length()
# (temp_vcf_combine[temp_vcf_combine$CALLER == called,]$POS %>% length()) - (temp_vcf_combine[temp_vcf_combine$CALLER == called,]$POS %>% unique() %>% length())

names(temp_vcf_combine)[names(temp_vcf_combine) == 'IDENTIFIER'] <- 'THESIS_IDENTIFIER'
names(temp_vcf_combine)[names(temp_vcf_combine) == 'COHORT'] <- 'THESIS_COHORT'
names(temp_vcf_combine)[names(temp_vcf_combine) == 'BLASTER'] <- 'THESIS_BLASTER'
names(temp_vcf_combine)[names(temp_vcf_combine) == 'CALLER'] <- 'THESIS_CALLER'

all_varscan <- temp_vcf_combine
colnames(all_varscan)
str(all_varscan %>% filter(FILTER == "PASS"))

data_unified <- all_varscan
file_data_unified <- "all_varscan"
dir_data_unified <- "/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified"
#write.csv(data_unified, paste0(dir_data_unified, "/", file_data_unified, ".csv"), row.names=FALSE, fileEncoding = "UTF-8")
saveRDS(data_unified, paste0(dir_data_unified, "/", file_data_unified, ".rds"), compress = TRUE)
infoRDS(paste0(dir_data_unified, "/", file_data_unified, ".rds"))
# readRDS(paste0(dir_data_unified, "/", file_data_unified, ".rds"))

data_fwrite <- data_unified
name_data_frwite <- "all_varscan"
new_filename <- paste0(github_data_unified, name_data_frwite, ".csv")
print(new_filename)
data.table::fwrite(data_fwrite, file = new_filename, append = FALSE, quote = "auto", showProgress = getOption("datatable.showProgress", interactive()), buffMB = 1024, nThread = 6)

##########################
##########################      UMIVAR VCFs
##########################
########################################################################################

# umivar
all_vcf_control_umivar
all_vcf_experimental_umivar
(colnames(all_vcf_control_umivar) == colnames(all_vcf_experimental_umivar))
colnames(all_vcf_experimental_umivar)
str(all_vcf_experimental_umivar)

called <- 'umivar'
temp_vcf_combine <- bind_rows(all_vcf_control_umivar, all_vcf_experimental_umivar) %>% unique %>% arrange(POS)
called <- temp_vcf_combine$CALLER %>% unique()
called

# BASE INFORMATION
temp_vcf_combine$THESIS_POS <- temp_vcf_combine[temp_vcf_combine$CALLER == called,]$POS %>% as.integer()
temp_vcf_combine$THESIS_BASE_CHANGE <- paste(temp_vcf_combine[temp_vcf_combine$CALLER == called,]$REF, temp_vcf_combine[temp_vcf_combine$CALLER == called,]$ALT, sep = ">")
temp_vcf_combine$THESIS_VARIANT <- paste(temp_vcf_combine[temp_vcf_combine$CALLER == called,]$THESIS_BASE_CHANGE, temp_vcf_combine[temp_vcf_combine$CALLER == called,]$POS, sep = ":")
temp_vcf_combine$THESIS_VARIANT_BY_POS <- paste(temp_vcf_combine[temp_vcf_combine$CALLER == called,]$POS, temp_vcf_combine[temp_vcf_combine$CALLER == called,]$THESIS_BASE_CHANGE, sep = ":")
temp_vcf_combine$THESIS_VARIANT_BY_HGVS <- paste0((temp_vcf_combine[temp_vcf_combine$CALLER == called,]$CHROM %>% str_replace("\\.", "_")), ".", temp_vcf_combine[temp_vcf_combine$CALLER == called,]$POS, temp_vcf_combine[temp_vcf_combine$CALLER == called,]$THESIS_BASE_CHANGE)
temp_vcf_combine$THESIS_FLAG_CONSVAR <- as.character(NA)
temp_vcf_combine$THESIS_FLAG_INDEL <- as.character(NA)
temp_vcf_combine$THESIS_FLAG_WILDTYPE <- as.character(NA)

temp_vcf_combine$THESIS_PVALUE_UNCORRECTED <- temp_vcf_combine[temp_vcf_combine$CALLER == called,]$gt_Pval %>% as.numeric()
temp_vcf_combine$THESIS_QVALUE <- temp_vcf_combine[temp_vcf_combine$CALLER == called,]$gt_FDR %>% as.numeric()
temp_vcf_combine$THESIS_DEPTH_TOTAL <- temp_vcf_combine[temp_vcf_combine$CALLER == called,]$gt_DP %>% as.integer()


# DP4     4      Integer Counts for ref-forward bases, ref-reverse, alt-forward and alt-reverse bases   
# (sapply(strsplit(temp_vcf_combine$gt_Strand,"-"), `[`, 1) %>% as.integer()) # alt-forward bases
# (sapply(strsplit(temp_vcf_combine$gt_Strand,"-"), `[`, 2) %>% as.integer()) # alt-reverse bases
# (sapply(strsplit(temp_vcf_combine$gt_Strand,"-"), `[`, 3) %>% as.integer()) # ref-forward
# (sapply(strsplit(temp_vcf_combine$gt_Strand,"-"), `[`, 4) %>% as.integer()) # ref-reverse
# rf <- sapply(strsplit(temp_vcf_combine[temp_vcf_combine$CALLER == called,]$gt_Strand,"-"), `[`, 3) %>% as.integer()
# rr <- sapply(strsplit(temp_vcf_combine[temp_vcf_combine$CALLER == called,]$gt_Strand,"-"), `[`, 4) %>% as.integer()
# rf+rr

temp_vcf_combine$THESIS_DEPTH_REFERENCE <- ((sapply(strsplit(temp_vcf_combine[temp_vcf_combine$CALLER == called,]$gt_Strand,"-"), `[`, 3)  %>% as.integer()) + (sapply(strsplit(temp_vcf_combine[temp_vcf_combine$CALLER == called,]$gt_Strand,"-"), `[`, 4)) %>% as.integer())
temp_vcf_combine$THESIS_FREQUENCY_NOT_ALTERNATE <- (1 - (temp_vcf_combine[temp_vcf_combine$CALLER == called,]$gt_AF %>% as.numeric()))
temp_vcf_combine$THESIS_QUALITY_REFERENCE <- as.numeric(NA) # temp_vcf_combine[temp_vcf_combine$CALLER == called,]$gt_REF_QUAL %>% as.numeric()

temp_vcf_combine$THESIS_DEPTH_ALTERNATE <- ((sapply(strsplit(temp_vcf_combine[temp_vcf_combine$CALLER == called,]$gt_Strand,"-"), `[`, 1)  %>% as.integer()) + (sapply(strsplit(temp_vcf_combine[temp_vcf_combine$CALLER == called,]$gt_Strand,"-"), `[`, 2)) %>% as.integer())
temp_vcf_combine$THESIS_FREQUENCY_ALTERNATE <- (temp_vcf_combine[temp_vcf_combine$CALLER == called,]$gt_AF %>% as.numeric())
temp_vcf_combine$THESIS_QUALITY_ALTERNATE <- as.numeric(NA) # temp_vcf_combine[temp_vcf_combine$CALLER == called,]$gt_ALT_QUAL %>% as.numeric()

temp_vcf_combine$THESIS_QC_DEPTH_TOTAL_DISCREPANCY <- (temp_vcf_combine$THESIS_DEPTH_TOTAL - (temp_vcf_combine$THESIS_DEPTH_REFERENCE + temp_vcf_combine$THESIS_DEPTH_ALTERNATE)) %>% as.integer()
temp_vcf_combine$THESIS_QC_DEPTH_TOTAL_CALC_LESSTHAN <- ((temp_vcf_combine$THESIS_DEPTH_REFERENCE + temp_vcf_combine$THESIS_DEPTH_ALTERNATE) < temp_vcf_combine$THESIS_DEPTH_TOTAL)
temp_vcf_combine$THESIS_QC_DEPTH_TOTAL_CALC_EQUAL <- ((temp_vcf_combine$THESIS_DEPTH_REFERENCE + temp_vcf_combine$THESIS_DEPTH_ALTERNATE) == temp_vcf_combine$THESIS_DEPTH_TOTAL)
temp_vcf_combine$THESIS_QC_CALLER_SPECIFIC_5_PERCENT <- as.logical(NA)

# temp_vcf_combine[temp_vcf_combine$CALLER == called,]$POS %>% length()
# temp_vcf_combine[temp_vcf_combine$CALLER == called,]$POS %>% unique() %>% length()
# (temp_vcf_combine[temp_vcf_combine$CALLER == called,]$POS %>% length()) - (temp_vcf_combine[temp_vcf_combine$CALLER == called,]$POS %>% unique() %>% length())

names(temp_vcf_combine)[names(temp_vcf_combine) == 'IDENTIFIER'] <- 'THESIS_IDENTIFIER'
names(temp_vcf_combine)[names(temp_vcf_combine) == 'COHORT'] <- 'THESIS_COHORT'
names(temp_vcf_combine)[names(temp_vcf_combine) == 'BLASTER'] <- 'THESIS_BLASTER'
names(temp_vcf_combine)[names(temp_vcf_combine) == 'CALLER'] <- 'THESIS_CALLER'

all_umivar <- temp_vcf_combine
colnames(all_umivar)
str(all_umivar %>% filter(FILTER == "PASS"))

data_unified <- all_umivar
file_data_unified <- "all_umivar"
dir_data_unified <- "/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified"
dir_data_unified_bulk <- "/mnt/bulk_data/THESIS/bulk_unified_data"
#write.csv(data_unified, paste0(dir_data_unified_bulk, "/", file_data_unified, ".csv"), row.names=FALSE, fileEncoding = "UTF-8")
saveRDS(data_unified, paste0(dir_data_unified, "/", file_data_unified, ".rds"), compress = TRUE)
infoRDS(paste0(dir_data_unified, "/", file_data_unified, ".rds"))
# readRDS(paste0(dir_data_unified, "/", file_data_unified, ".rds"))

data_fwrite <- data_unified
name_data_frwite <- "all_umivar"
new_filename <- paste0("/mnt/bulk_data/THESIS/bulk_unified_data/", name_data_frwite, ".csv")
print(new_filename)
data.table::fwrite(data_fwrite, file = new_filename, append = FALSE, quote = "auto", showProgress = getOption("datatable.showProgress", interactive()), buffMB = 1024, nThread = 6)

##########################
##########################      FILTER VCFs
##########################
########################################################################################
# 
# str(all_ivar %>% filter(FILTER == "PASS"))
# str(all_lofreq %>% filter(FILTER == "PASS"))
# str(all_varscan %>% filter(FILTER == "PASS"))
# str(all_umivar %>% filter(FILTER == "PASS"))
# 
# colnames(all_ivar %>% filter(FILTER == "PASS"))
# colnames(all_lofreq %>% filter(FILTER == "PASS"))
# colnames(all_varscan %>% filter(FILTER == "PASS"))
# colnames(all_umivar %>% filter(FILTER == "PASS"))
# 
# ncol(all_ivar %>% filter(FILTER == "PASS"))
# ncol(all_lofreq %>% filter(FILTER == "PASS"))
# ncol(all_varscan %>% filter(FILTER == "PASS"))
# ncol(all_umivar %>% filter(FILTER == "PASS"))
# 
# nrow(all_ivar %>% filter(FILTER == "PASS"))
# nrow(all_lofreq %>% filter(FILTER == "PASS"))
# nrow(all_varscan %>% filter(FILTER == "PASS"))
# nrow(all_umivar %>% filter(FILTER == "PASS"))

##########################
##########################      COMBINE ALL VCFs
##########################
########################################################################################
all_samples <- bind_rows(all_lofreq, all_ivar, all_umivar, all_varscan) %>% arrange(POS) # min ncol first

data_unified <- all_samples
file_data_unified <- "all_callers"
dir_data_unified <- "/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified"
dir_data_unified_bulk <- "/mnt/bulk_data/THESIS/bulk_unified_data"
#write.csv(data_unified, paste0(dir_data_unified_bulk, "/", file_data_unified, ".csv"), row.names=FALSE, fileEncoding = "UTF-8")
saveRDS(data_unified, paste0(dir_data_unified, "/", file_data_unified, ".rds"), compress = TRUE)
infoRDS(paste0(dir_data_unified, "/", file_data_unified, ".rds"))
# readRDS(paste0(dir_data_unified, "/", file_data_unified, ".rds"))

data_fwrite <- data_unified
name_data_frwite <- "all_callers"
new_filename <- paste0("/mnt/bulk_data/THESIS/bulk_unified_data/", name_data_frwite, ".csv")
print(new_filename)
data.table::fwrite(data_fwrite, file = new_filename, append = FALSE, quote = "auto", showProgress = getOption("datatable.showProgress", interactive()), buffMB = 1024, nThread = 6)

##########################
##########################      FEATURES
##########################
########################################################################################
features_gff <- gff_join %>% select(start, end, starts_with("THESIS_"))

data_unified <- features_gff
file_data_unified <- "features_gff"
dir_data_unified <- "/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified"
dir_data_unified_bulk <- "/mnt/bulk_data/THESIS/bulk_unified_data"
#write.csv(data_unified, paste0(dir_data_unified, "/", file_data_unified, ".csv"), row.names=FALSE, fileEncoding = "UTF-8")
saveRDS(data_unified, paste0(dir_data_unified, "/", file_data_unified, ".rds"), compress = TRUE)
infoRDS(paste0(dir_data_unified, "/", file_data_unified, ".rds"))
# readRDS(paste0(dir_data_unified, "/", file_data_unified, ".rds"))


data_fwrite <- data_unified
name_data_frwite <- "features_gff"
new_filename <- paste0("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/", name_data_frwite, ".csv")
print(new_filename)
data.table::fwrite(data_fwrite, file = new_filename, append = FALSE, quote = "auto", showProgress = getOption("datatable.showProgress", interactive()), buffMB = 1024, nThread = 6)

##########################
##########################      FEATURES PLUS CALLERS
##########################
########################################################################################

all_samples$THESIS_FEATURE_LABEL <- as.character(NA)
all_samples$THESIS_FEATURE_PRODUCT <- as.character(NA)
all_samples$THESIS_FEATURE_GENE <- as.character(NA)
all_samples$THESIS_FEATURE_GENOMIC <- as.character(NA)
all_samples$THESIS_FEATURE_NCBI <- as.character(NA)
all_samples$THESIS_FEATURE_ID <- as.character(NA)

for (range_start in gff_join$start) {
  bound_lower <- range_start
  bound_upper <- gff_join[gff_join$start == range_start,]$end
  # all_samples_PASS[all_samples_PASS$POS >= 1 & all_samples_PASS$POS <= 265,]$POS
  all_samples[all_samples$POS >= bound_lower & all_samples$POS <= bound_upper,]$THESIS_FEATURE_LABEL <- gff_join[gff_join$start == bound_lower,]$THESIS_FEATURE_LABEL
  all_samples[all_samples$POS >= bound_lower & all_samples$POS <= bound_upper,]$THESIS_FEATURE_PRODUCT <- gff_join[gff_join$start == range_start,]$THESIS_FEATURE_PRODUCT
  all_samples[all_samples$POS >= bound_lower & all_samples$POS <= bound_upper,]$THESIS_FEATURE_GENE <- gff_join[gff_join$start == range_start,]$THESIS_FEATURE_GENE
  all_samples[all_samples$POS >= bound_lower & all_samples$POS <= bound_upper,]$THESIS_FEATURE_GENOMIC <- gff_join[gff_join$start == range_start,]$THESIS_FEATURE_GENOMIC
  all_samples[all_samples$POS >= bound_lower & all_samples$POS <= bound_upper,]$THESIS_FEATURE_NCBI <- gff_join[gff_join$start == range_start,]$THESIS_FEATURE_NCBI
  all_samples[all_samples$POS >= bound_lower & all_samples$POS <= bound_upper,]$THESIS_FEATURE_ID <- gff_join[gff_join$start == range_start,]$THESIS_FEATURE_ID
}

all_samples$THESIS_FEATURE_LABEL %>% unique()
all_samples$THESIS_FEATURE_PRODUCT %>% unique()
all_samples$THESIS_FEATURE_GENE %>% unique()
all_samples$THESIS_FEATURE_GENOMIC %>% unique()
all_samples$THESIS_FEATURE_NCBI %>% unique()
all_samples$THESIS_FEATURE_ID %>% unique()

colnames(all_samples)
ncol(all_samples)
nrow(all_samples)
str(all_samples)

data_unified <- all_samples
file_data_unified <- "all_callers_with_features"
dir_data_unified <- "/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified"
dir_data_unified_bulk <- "/mnt/bulk_data/THESIS/bulk_unified_data"
#write.csv(data_unified, paste0(dir_data_unified_bulk, "/", file_data_unified, ".csv"), row.names=FALSE, fileEncoding = "UTF-8")
saveRDS(data_unified, paste0(dir_data_unified, "/", file_data_unified, ".rds"), compress = TRUE)
infoRDS(paste0(dir_data_unified, "/", file_data_unified, ".rds"))
# readRDS(paste0(dir_data_unified, "/", file_data_unified, ".rds"))

data_fwrite <- data_unified
name_data_frwite <- "all_callers_with_features"
new_filename <- paste0("/mnt/bulk_data/THESIS/bulk_unified_data/", name_data_frwite, ".csv")
print(new_filename)
data.table::fwrite(data_fwrite, file = new_filename, append = FALSE, quote = "auto", showProgress = getOption("datatable.showProgress", interactive()), buffMB = 1024, nThread = 6)

##########################
##########################      NEW METADATA
##########################
########################################################################################

all_samples_PASS <- all_samples %>% filter(FILTER == "PASS") %>% unique %>% arrange(THESIS_IDENTIFIER, POS)
all_samples_FAIL <- all_samples %>% filter(FILTER != "PASS") %>% unique %>% arrange(THESIS_IDENTIFIER, POS)

list_experimental <- (all_umivar %>% filter(THESIS_COHORT =="experimental"))$THESIS_IDENTIFIER %>% sort() %>% unique()
list_control <- (all_umivar %>% filter(THESIS_COHORT =="control"))$THESIS_IDENTIFIER %>% sort() %>% unique()

all_samples_experimental <- all_samples_PASS %>% filter(THESIS_IDENTIFIER %in% list_experimental) %>% unique %>% arrange(THESIS_IDENTIFIER, POS)
all_samples_control <- all_samples_PASS %>% filter(THESIS_IDENTIFIER %in% list_control) %>% unique %>% arrange(THESIS_IDENTIFIER, POS)



all_samples_PASS %>% filter(THESIS_FEATURE_LABEL == "S") %>% unique %>% arrange(POS)
# all_samples_PASS %>% filter(THESIS_FREQUENCY_ALTERNATE <0.5) %>% unique %>% arrange(THESIS_IDENTIFIER, POS)
# all_samples_PASS %>% filter(THESIS_FREQUENCY_ALTERNATE <0.1) %>% unique %>% arrange(THESIS_IDENTIFIER, POS)
# all_samples_PASS %>% filter(THESIS_DEPTH_TOTAL > 1000) %>% filter(THESIS_FREQUENCY_ALTERNATE <0.01) %>% unique %>% arrange(THESIS_IDENTIFIER, POS) %>% view()

##########################
##########################
########################################################################################
file_all_snp_data_header <- "/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_generated/column.header.csv"
file_all_metadata <- "/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_generated/metadata.csv" # "/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_computed/column.header.csv"
metadata_samples <- vroom(file_all_metadata)

metadata_new <- metadata_samples %>% select(METADATA_NAME_IDENTIFIER, METADATA_ID, METADATA_SUBMITTER, METADATA_NAME_DESCRIPTION, G_FULL_NAME, M_ABNAHMEDATUM, METADATA_PL_LINEAGE, METADATA_PL_TAXON, METADATA_CT_S, METADATA_CT_E, METADATA_COPIES_TOTAL, METADATA_VARIANT_FRACTION_C1_AUSTRALIA, METADATA_VARIANT_FRACTION_C2_ANCESTRAL,	METADATA_VARIANT_FRACTION_C14_ALPHA,	METADATA_VARIANT_FRACTION_C29_DELTA,	METADATA_VARIANT_FRACTION_C48_OMICRON, METADATA_J_RUN_NAME,	METADATA_J_RUN_FLOWCELL_ID,	METADATA_J_RUN_FLOWCELL_TYPE,	METADATA_J_RUN_RECIPE, ZZZ_PROJECT_J, ZZZ_PROJECT, G_REAL_FOLDER_PATH,	G_REAL_PATH_PROJECT,	METADATA_G_RUNS_HOW_MANY,	METADATA_G_LANES_FORWARD,	METADATA_G_LANES_REVERSE,	METADATA_G_LANES_INDEX)

metadata_new$THESIS_SAMPLE_IDENTIFIER <- metadata_new$G_FULL_NAME
metadata_new$THESIS_SAMPLE_CT_S <- metadata_new$METADATA_CT_S %>% as.numeric()
metadata_new$THESIS_SAMPLE_CT_E <- metadata_new$METADATA_CT_E %>% as.numeric()
metadata_new$THESIS_SAMPLE_COLLECTION_DATE <- metadata_new$M_ABNAHMEDATUM %>% ymd()
metadata_new$THESIS_SAMPLE_ASSIGNED_LINEAGE <- metadata_new$METADATA_PL_LINEAGE
metadata_new$THESIS_SAMPLE_ASSIGNED_TAXON <- metadata_new$METADATA_PL_TAXON

metadata_new$THESIS_SAMPLE_PROJECT <- metadata_new$G_REAL_PATH_PROJECT
metadata_new$THESIS_SAMPLE_SEQ_FLOWCELL_ID <- metadata_new$METADATA_J_RUN_FLOWCELL_ID
metadata_new$THESIS_SAMPLE_SEQ_FLOWCELL_TYPE <- metadata_new$METADATA_J_RUN_FLOWCELL_TYPE
metadata_new$THESIS_SAMPLE_SEQ_RECIPE <- metadata_new$METADATA_J_RUN_RECIPE
metadata_new$THESIS_SAMPLE_SEQ_RUNS <- metadata_new$METADATA_G_RUNS_HOW_MANY %>% as.integer()
metadata_new$THESIS_SAMPLE_SEQ_LANES_FORWARD <- metadata_new$METADATA_G_LANES_FORWARD %>% as.integer()
metadata_new$THESIS_SAMPLE_SEQ_LANES_REVERSE <- metadata_new$METADATA_G_LANES_REVERSE %>% as.integer()
metadata_new$THESIS_SAMPLE_SEQ_LANES_INDEX <- metadata_new$METADATA_G_LANES_INDEX %>% as.integer()
metadata_new$THESIS_SYNTHETIC_COPIES_TOTAL <- metadata_new$METADATA_COPIES_TOTAL %>% as.integer()

metadata_new$THESIS_SYNTHETIC_PRESENT_UHR <- as.logical(FALSE)
metadata_new$THESIS_SYNTHETIC_DILUTION <- as.logical(FALSE)
metadata_new$THESIS_SYNTHETIC_PRESENT_C1_AUSTRALIA_MT007544_1 <- as.logical(FALSE)
metadata_new$THESIS_SYNTHETIC_PRESENT_C1_AUSTRALIA_MT007544_1_FRACTION <- metadata_new$METADATA_VARIANT_FRACTION_C1_AUSTRALIA/100 %>% as.numeric()
metadata_new$THESIS_SYNTHETIC_PRESENT_C1_AUSTRALIA_MT007544_1_COPIES <- (metadata_new$THESIS_SYNTHETIC_COPIES_TOTAL * metadata_new$THESIS_SYNTHETIC_PRESENT_C1_AUSTRALIA_MT007544_1_FRACTION) %>% as.integer()

metadata_new$THESIS_SYNTHETIC_PRESENT_C2_ANCESTRAL_MN908947_3 <- as.logical(FALSE)
metadata_new$THESIS_SYNTHETIC_PRESENT_C2_ANCESTRAL_MN908947_3_FRACTION <- metadata_new$METADATA_VARIANT_FRACTION_C2_ANCESTRAL/100 %>% as.numeric()
metadata_new$THESIS_SYNTHETIC_PRESENT_C2_ANCESTRAL_MN908947_3_COPIES <- (metadata_new$THESIS_SYNTHETIC_COPIES_TOTAL * metadata_new$THESIS_SYNTHETIC_PRESENT_C2_ANCESTRAL_MN908947_3_FRACTION) %>% as.integer()

metadata_new$THESIS_SYNTHETIC_PRESENT_C14_ALPHA_EPI_ISL_710528 <- as.logical(FALSE)
metadata_new$THESIS_SYNTHETIC_PRESENT_C14_ALPHA_EPI_ISL_710528_FRACTION <- metadata_new$METADATA_VARIANT_FRACTION_C14_ALPHA/100 %>% as.numeric()
metadata_new$THESIS_SYNTHETIC_PRESENT_C14_ALPHA_EPI_ISL_710528_COPIES <-(metadata_new$THESIS_SYNTHETIC_COPIES_TOTAL * metadata_new$THESIS_SYNTHETIC_PRESENT_C14_ALPHA_EPI_ISL_710528_FRACTION) %>% as.integer()

metadata_new$THESIS_SYNTHETIC_PRESENT_C29_DELTA_EPI_ISL_2693246 <- as.logical(FALSE)
metadata_new$THESIS_SYNTHETIC_PRESENT_C29_DELTA_EPI_ISL_2693246_FRACTION <- metadata_new$METADATA_VARIANT_FRACTION_C29_DELTA/100 %>% as.numeric()
metadata_new$THESIS_SYNTHETIC_PRESENT_C29_DELTA_EPI_ISL_2693246_COPIES <- (metadata_new$THESIS_SYNTHETIC_COPIES_TOTAL * metadata_new$THESIS_SYNTHETIC_PRESENT_C29_DELTA_EPI_ISL_2693246_FRACTION) %>% as.integer()

metadata_new$THESIS_SYNTHETIC_PRESENT_C48_OMICRON_EPI_ISL_6841980 <- as.logical(FALSE)
metadata_new$THESIS_SYNTHETIC_PRESENT_C48_OMICRON_EPI_ISL_6841980_FRACTION <- metadata_new$METADATA_VARIANT_FRACTION_C48_OMICRON/100 %>% as.numeric()
metadata_new$THESIS_SYNTHETIC_PRESENT_C48_OMICRON_EPI_ISL_6841980_COPIES <- (metadata_new$THESIS_SYNTHETIC_COPIES_TOTAL * metadata_new$THESIS_SYNTHETIC_PRESENT_C48_OMICRON_EPI_ISL_6841980_FRACTION) %>% as.integer()

metadata_new$THESIS_SYNTHETIC_PRESENT_UHR <- grepl("UHR", metadata_new$METADATA_NAME_DESCRIPTION, ignore.case = TRUE)
metadata_new$THESIS_SYNTHETIC_PRESENT_UHR <- ifelse(grepl("without UHR", metadata_new$METADATA_NAME_DESCRIPTION, ignore.case = TRUE), FALSE, metadata_new$THESIS_SYNTHETIC_PRESENT_UHR) 

metadata_new$THESIS_SYNTHETIC_DILUTION <- grepl("DIL-", metadata_new$METADATA_NAME_DESCRIPTION, ignore.case = TRUE)

# THESIS_SYNTHETIC_PRESENT_C1_AUSTRALIA_MT007544_1
metadata_new$THESIS_SYNTHETIC_PRESENT_C1_AUSTRALIA_MT007544_1 <- ifelse(grepl("COVID-1 ", metadata_new$METADATA_NAME_DESCRIPTION, ignore.case = TRUE), TRUE, metadata_new$THESIS_SYNTHETIC_PRESENT_C1_AUSTRALIA_MT007544_1) 
metadata_new$THESIS_SYNTHETIC_PRESENT_C1_AUSTRALIA_MT007544_1 <- ifelse(grepl("MIX COVID-1/2", metadata_new$METADATA_NAME_DESCRIPTION, ignore.case = TRUE), TRUE, metadata_new$THESIS_SYNTHETIC_PRESENT_C1_AUSTRALIA_MT007544_1) 
metadata_new$THESIS_SYNTHETIC_PRESENT_C1_AUSTRALIA_MT007544_1 <- ifelse(grepl("MIX 1/2", metadata_new$METADATA_NAME_DESCRIPTION, ignore.case = TRUE), TRUE, metadata_new$THESIS_SYNTHETIC_PRESENT_C1_AUSTRALIA_MT007544_1) 
metadata_new$THESIS_SYNTHETIC_PRESENT_C1_AUSTRALIA_MT007544_1 <- ifelse(grepl("C1$", metadata_new$METADATA_NAME_DESCRIPTION, ignore.case = TRUE), TRUE, metadata_new$THESIS_SYNTHETIC_PRESENT_C1_AUSTRALIA_MT007544_1)
metadata_new$THESIS_SYNTHETIC_PRESENT_C1_AUSTRALIA_MT007544_1 <- ifelse(grepl("C1,", metadata_new$METADATA_NAME_DESCRIPTION, ignore.case = TRUE), TRUE, metadata_new$THESIS_SYNTHETIC_PRESENT_C1_AUSTRALIA_MT007544_1)
metadata_new$THESIS_SYNTHETIC_PRESENT_C1_AUSTRALIA_MT007544_1 <- ifelse(grepl("CTRL1", metadata_new$METADATA_NAME_DESCRIPTION, ignore.case = TRUE), TRUE, metadata_new$THESIS_SYNTHETIC_PRESENT_C1_AUSTRALIA_MT007544_1)
metadata_new$THESIS_SYNTHETIC_PRESENT_C1_AUSTRALIA_MT007544_1 <- ifelse(grepl("CTRL1$", metadata_new$METADATA_NAME_DESCRIPTION, ignore.case = TRUE), TRUE, metadata_new$THESIS_SYNTHETIC_PRESENT_C1_AUSTRALIA_MT007544_1)

# THESIS_SYNTHETIC_PRESENT_C2_ANCESTRAL_MN908947_3
metadata_new$THESIS_SYNTHETIC_PRESENT_C2_ANCESTRAL_MN908947_3 <- ifelse(grepl("MIX COVID-1/2", metadata_new$METADATA_NAME_DESCRIPTION, ignore.case = TRUE), TRUE, metadata_new$THESIS_SYNTHETIC_PRESENT_C2_ANCESTRAL_MN908947_3) 
metadata_new$THESIS_SYNTHETIC_PRESENT_C2_ANCESTRAL_MN908947_3 <- ifelse(grepl("MIX 1/2", metadata_new$METADATA_NAME_DESCRIPTION, ignore.case = TRUE), TRUE, metadata_new$THESIS_SYNTHETIC_PRESENT_C2_ANCESTRAL_MN908947_3) 
metadata_new$THESIS_SYNTHETIC_PRESENT_C2_ANCESTRAL_MN908947_3 <- ifelse(grepl("C2$", metadata_new$METADATA_NAME_DESCRIPTION, ignore.case = TRUE), TRUE, metadata_new$THESIS_SYNTHETIC_PRESENT_C2_ANCESTRAL_MN908947_3)
metadata_new$THESIS_SYNTHETIC_PRESENT_C2_ANCESTRAL_MN908947_3 <- ifelse(grepl("C2,", metadata_new$METADATA_NAME_DESCRIPTION, ignore.case = TRUE), TRUE, metadata_new$THESIS_SYNTHETIC_PRESENT_C2_ANCESTRAL_MN908947_3)
metadata_new$THESIS_SYNTHETIC_PRESENT_C2_ANCESTRAL_MN908947_3 <- ifelse(grepl("CTRL2", metadata_new$METADATA_NAME_DESCRIPTION, ignore.case = TRUE), TRUE, metadata_new$THESIS_SYNTHETIC_PRESENT_C2_ANCESTRAL_MN908947_3)
metadata_new$THESIS_SYNTHETIC_PRESENT_C2_ANCESTRAL_MN908947_3 <- ifelse(grepl("CTRL2$", metadata_new$METADATA_NAME_DESCRIPTION, ignore.case = TRUE), TRUE, metadata_new$THESIS_SYNTHETIC_PRESENT_C2_ANCESTRAL_MN908947_3)

# THESIS_SYNTHETIC_PRESENT_C14_ALPHA_EPI_ISL_710528
metadata_new$THESIS_SYNTHETIC_PRESENT_C14_ALPHA_EPI_ISL_710528 <- ifelse(grepl("Ctrl 14,", metadata_new$METADATA_NAME_DESCRIPTION, ignore.case = TRUE), TRUE, metadata_new$THESIS_SYNTHETIC_PRESENT_C14_ALPHA_EPI_ISL_710528) 
metadata_new$THESIS_SYNTHETIC_PRESENT_C14_ALPHA_EPI_ISL_710528 <- ifelse(grepl("C14,", metadata_new$METADATA_NAME_DESCRIPTION, ignore.case = TRUE), TRUE, metadata_new$THESIS_SYNTHETIC_PRESENT_C14_ALPHA_EPI_ISL_710528) 
metadata_new$THESIS_SYNTHETIC_PRESENT_C14_ALPHA_EPI_ISL_710528 <- ifelse(grepl("C14$", metadata_new$METADATA_NAME_DESCRIPTION, ignore.case = TRUE), TRUE, metadata_new$THESIS_SYNTHETIC_PRESENT_C14_ALPHA_EPI_ISL_710528)

# THESIS_SYNTHETIC_PRESENT_C29_DELTA_EPI_ISL_2693246
metadata_new$THESIS_SYNTHETIC_PRESENT_C29_DELTA_EPI_ISL_2693246 <- ifelse(grepl("DELTA", metadata_new$METADATA_NAME_DESCRIPTION, ignore.case = TRUE), TRUE, metadata_new$THESIS_SYNTHETIC_PRESENT_C29_DELTA_EPI_ISL_2693246) 
metadata_new$THESIS_SYNTHETIC_PRESENT_C29_DELTA_EPI_ISL_2693246 <- ifelse(grepl("C29,", metadata_new$METADATA_NAME_DESCRIPTION, ignore.case = TRUE), TRUE, metadata_new$THESIS_SYNTHETIC_PRESENT_C29_DELTA_EPI_ISL_2693246) 
metadata_new$THESIS_SYNTHETIC_PRESENT_C29_DELTA_EPI_ISL_2693246 <- ifelse(grepl("C29$", metadata_new$METADATA_NAME_DESCRIPTION, ignore.case = TRUE), TRUE, metadata_new$THESIS_SYNTHETIC_PRESENT_C29_DELTA_EPI_ISL_2693246)

# THESIS_SYNTHETIC_PRESENT_C48_OMICRON_EPI_ISL_6841980
metadata_new$THESIS_SYNTHETIC_PRESENT_C48_OMICRON_EPI_ISL_6841980 <- ifelse(grepl("omicr", metadata_new$METADATA_NAME_DESCRIPTION, ignore.case = TRUE), TRUE, metadata_new$THESIS_SYNTHETIC_PRESENT_C48_OMICRON_EPI_ISL_6841980) 
metadata_new$THESIS_SYNTHETIC_PRESENT_C48_OMICRON_EPI_ISL_6841980 <- ifelse(grepl("C48,", metadata_new$METADATA_NAME_DESCRIPTION, ignore.case = TRUE), TRUE, metadata_new$THESIS_SYNTHETIC_PRESENT_C48_OMICRON_EPI_ISL_6841980) 
metadata_new$THESIS_SYNTHETIC_PRESENT_C48_OMICRON_EPI_ISL_6841980 <- ifelse(grepl("C48$", metadata_new$METADATA_NAME_DESCRIPTION, ignore.case = TRUE), TRUE, metadata_new$THESIS_SYNTHETIC_PRESENT_C48_OMICRON_EPI_ISL_6841980)

metadata_new$THESIS_SYNTHETIC_COPIES_SUM_QC <- (metadata_new$THESIS_SYNTHETIC_PRESENT_C1_AUSTRALIA_MT007544_1_COPIES +  metadata_new$THESIS_SYNTHETIC_PRESENT_C2_ANCESTRAL_MN908947_3_COPIES +  metadata_new$THESIS_SYNTHETIC_PRESENT_C14_ALPHA_EPI_ISL_710528_COPIES + metadata_new$THESIS_SYNTHETIC_PRESENT_C29_DELTA_EPI_ISL_2693246_COPIES + metadata_new$THESIS_SYNTHETIC_PRESENT_C48_OMICRON_EPI_ISL_6841980_COPIES)

metadata_new$THESIS_SYNTHETIC_COPIES_SUM_QC_TRANSCRIPT <- paste("C1", metadata_new$THESIS_SYNTHETIC_PRESENT_C1_AUSTRALIA_MT007544_1_COPIES, "C2", metadata_new$THESIS_SYNTHETIC_PRESENT_C2_ANCESTRAL_MN908947_3_COPIES, "C14", metadata_new$THESIS_SYNTHETIC_PRESENT_C14_ALPHA_EPI_ISL_710528_COPIES, "C29", metadata_new$THESIS_SYNTHETIC_PRESENT_C29_DELTA_EPI_ISL_2693246_COPIES, "C48", metadata_new$THESIS_SYNTHETIC_PRESENT_C48_OMICRON_EPI_ISL_6841980_COPIES, "SUM_SYNTHETIC_COPIES",  metadata_new$THESIS_SYNTHETIC_COPIES_SUM_QC, sep = ":")

metadata_new$THESIS_SYNTHETIC_FRACTION_SUM_QC <- (metadata_new$THESIS_SYNTHETIC_PRESENT_C1_AUSTRALIA_MT007544_1_FRACTION +  metadata_new$THESIS_SYNTHETIC_PRESENT_C2_ANCESTRAL_MN908947_3_FRACTION +  metadata_new$THESIS_SYNTHETIC_PRESENT_C14_ALPHA_EPI_ISL_710528_FRACTION + metadata_new$THESIS_SYNTHETIC_PRESENT_C29_DELTA_EPI_ISL_2693246_FRACTION + metadata_new$THESIS_SYNTHETIC_PRESENT_C48_OMICRON_EPI_ISL_6841980_FRACTION)

metadata_new$THESIS_SYNTHETIC_FRACTION_SUM_QC_TRANSCRIPT <- paste("C1", metadata_new$THESIS_SYNTHETIC_PRESENT_C1_AUSTRALIA_MT007544_1_FRACTION, "C2", metadata_new$THESIS_SYNTHETIC_PRESENT_C2_ANCESTRAL_MN908947_3_FRACTION, "C14", metadata_new$THESIS_SYNTHETIC_PRESENT_C14_ALPHA_EPI_ISL_710528_FRACTION, "C29", metadata_new$THESIS_SYNTHETIC_PRESENT_C29_DELTA_EPI_ISL_2693246_FRACTION, "C48", metadata_new$THESIS_SYNTHETIC_PRESENT_C48_OMICRON_EPI_ISL_6841980_FRACTION, "SUM_SYNTHETIC_FRACTIONS",  metadata_new$THESIS_SYNTHETIC_FRACTION_SUM_QC, sep = ":")

metadata_new$THESIS_SYNTHETIC_VARIANT_CONTROL_QC_TRANSCRIPT <- paste("C1", metadata_new$THESIS_SYNTHETIC_PRESENT_C1_AUSTRALIA_MT007544_1, "C2", metadata_new$THESIS_SYNTHETIC_PRESENT_C2_ANCESTRAL_MN908947_3, "C14", metadata_new$THESIS_SYNTHETIC_PRESENT_C14_ALPHA_EPI_ISL_710528, "C29", metadata_new$THESIS_SYNTHETIC_PRESENT_C29_DELTA_EPI_ISL_2693246, "C48", metadata_new$THESIS_SYNTHETIC_PRESENT_C48_OMICRON_EPI_ISL_6841980, sep = ":")

#metadata_new %>% select(METADATA_NAME_DESCRIPTION, starts_with("THESIS")) %>% view()

data_unified <- metadata_new
file_data_unified <- "metadata_new"
dir_data_unified <- "/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified"
dir_data_unified_bulk <- "/mnt/bulk_data/THESIS/bulk_unified_data"
#write.csv(data_unified, paste0(dir_data_unified, "/", file_data_unified, ".csv"), row.names=FALSE, fileEncoding = "UTF-8")
saveRDS(data_unified, paste0(dir_data_unified, "/", file_data_unified, ".rds"), compress = TRUE)
infoRDS(paste0(dir_data_unified, "/", file_data_unified, ".rds"))
# readRDS(paste0(dir_data_unified, "/", file_data_unified, ".rds"))

data_fwrite <- data_unified
name_data_frwite <- "metadata_new"
# new_filename <- paste0("/mnt/bulk_data/THESIS/bulk_unified_data/", name_data_frwite, ".csv")
new_filename <- paste0("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/", name_data_frwite, ".csv")
print(new_filename)
data.table::fwrite(data_fwrite, file = new_filename, append = FALSE, quote = "auto", showProgress = getOption("datatable.showProgress", interactive()), buffMB = 1024, nThread = 6)

##########################
##########################     metadata_goldmaster
##########################
########################################################################################

metadata_goldmaster <- metadata_new %>% select(METADATA_NAME_DESCRIPTION, starts_with("THESIS")) 

data_unified <- metadata_goldmaster
file_data_unified <- "metadata_goldmaster"
dir_data_unified <- "/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified"
dir_data_unified_bulk <- "/mnt/bulk_data/THESIS/bulk_unified_data"
write.csv(data_unified, paste0(dir_data_unified, "/", file_data_unified, ".csv"), row.names=FALSE, fileEncoding = "UTF-8")
saveRDS(data_unified, paste0(dir_data_unified, "/", file_data_unified, ".rds"), compress = TRUE)
infoRDS(paste0(dir_data_unified, "/", file_data_unified, ".rds"))
# readRDS(paste0(dir_data_unified, "/", file_data_unified, ".rds"))

data_fwrite <- data_unified
name_data_frwite <- "metadata_goldmaster"
# new_filename <- paste0("/mnt/bulk_data/THESIS/bulk_unified_data/", name_data_frwite, ".csv")
new_filename <- paste0("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/", name_data_frwite, ".csv")
print(new_filename)
data.table::fwrite(data_fwrite, file = new_filename, append = FALSE, quote = "auto", showProgress = getOption("datatable.showProgress", interactive()), buffMB = 1024, nThread = 6)

data_fwrite <- data_unified
name_data_frwite <- "metadata_new"
# new_filename <- paste0("/mnt/bulk_data/THESIS/bulk_unified_data/", name_data_frwite, ".csv")
new_filename <- paste0("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/", name_data_frwite, ".csv")
print(new_filename)
data.table::fwrite(data_fwrite, file = new_filename, append = FALSE, quote = "auto", showProgress = getOption("datatable.showProgress", interactive()), buffMB = 1024, nThread = 6)

##########################
##########################     all_callers_with_features_and_metadata
##########################
########################################################################################
rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() #free up memrory and report the memory usage.

library(tidyverse)
library(vroom)
library(lubridate)

file_data_unified <- "metadata_goldmaster"
dir_data_unified <- "/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified"
dir_data_unified_bulk <- "/mnt/bulk_data/THESIS/bulk_unified_data"
infoRDS(paste0(dir_data_unified, "/", file_data_unified, ".rds"))
metadata_goldmaster <- readRDS(paste0(dir_data_unified, "/", file_data_unified, ".rds"))

file_data_unified <- "all_callers_with_features" # all_samples
dir_data_unified <- "/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified"
dir_data_unified_bulk <- "/mnt/bulk_data/THESIS/bulk_unified_data"
infoRDS(paste0(dir_data_unified, "/", file_data_unified, ".rds"))
all_samples <- readRDS(paste0(dir_data_unified, "/", file_data_unified, ".rds"))

a <- metadata_goldmaster$THESIS_SAMPLE_IDENTIFIER %>% `[`(!is.na(.)) %>% unique() %>% sort()
b <- all_samples$THESIS_IDENTIFIER %>% `[`(!is.na(.)) %>% unique() %>% sort()
samples_sequenced <- intersect(a,b) %>% `[`(!is.na(.)) %>% unique() %>% sort()

metadata_goldmaster_sequenced <- metadata_goldmaster %>% filter(THESIS_SAMPLE_IDENTIFIER %in% samples_sequenced)

all_callers_with_features_and_metadata <- left_join(all_samples, metadata_goldmaster_sequenced, by = c("THESIS_IDENTIFIER"="THESIS_SAMPLE_IDENTIFIER"))

data_unified <- all_callers_with_features_and_metadata
file_data_unified <- "all_callers_with_features_and_metadata"
dir_data_unified <- "/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified"
dir_data_unified_bulk <- "/mnt/bulk_data/THESIS/bulk_unified_data"
#write.csv(data_unified, paste0(dir_data_unified_bulk, "/", file_data_unified, ".csv"), row.names=FALSE, fileEncoding = "UTF-8")
saveRDS(data_unified, paste0(dir_data_unified, "/", file_data_unified, ".rds"), compress = TRUE)
infoRDS(paste0(dir_data_unified, "/", file_data_unified, ".rds"))
# readRDS(paste0(dir_data_unified, "/", file_data_unified, ".rds"))

data_fwrite <- data_unified
name_data_frwite <- "all_callers_with_features_and_metadata"
# new_filename <- paste0("/mnt/bulk_data/THESIS/bulk_unified_data/", name_data_frwite, ".csv")
new_filename <- paste0("/mnt/bulk_data/THESIS/bulk_unified_data/", name_data_frwite, ".csv")
print(new_filename)
data.table::fwrite(data_fwrite, file = new_filename, append = FALSE, quote = "auto", showProgress = getOption("datatable.showProgress", interactive()), buffMB = 1024, nThread = 6)

##########################
##########################     GINI STATS
##########################
########################################################################################
rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() #free up memrory and report the memory usage.

library(tidyverse)
library(vroom)
library(lubridate)

dir_mpileups_experimental <- "/mnt/bulk_data/THESIS/DATA_THESIS/GILOT_DATA_THESIS/EXPERIMENTAL/counts/"
dir_mpileups <- dir_mpileups_experimental
unblasted_experimental <- "*_unblasted_gilot_experimental_counts.tsv"
blasted_experimental <- "*_blasted_gilot_experimental_counts.tsv"
unblasted_experimental_list <- list.files(dir_mpileups, full.names = TRUE,  pattern = unblasted_experimental)
blasted_experimental_list <- list.files(dir_mpileups, full.names = TRUE,  pattern = blasted_experimental)
unblasted_experimental_counts <- vroom(unblasted_experimental_list)
blasted_experimental_counts <- vroom(blasted_experimental_list)

dir_mpileups_control <- "/mnt/bulk_data/THESIS/DATA_THESIS/GILOT_DATA_THESIS/CONTROL/counts/"
dir_mpileups <- dir_mpileups_control
unblasted_control <- "*_unblasted_gilot_all_controls_counts.tsv"
blasted_control <- "*_blasted_gilot_all_controls_counts.tsv"
unblasted_experimental_list <- list.files(dir_mpileups, full.names = TRUE,  pattern = unblasted_control)
blasted_experimental_list <- list.files(dir_mpileups, full.names = TRUE,  pattern = blasted_control)
unblasted_control_counts <- vroom(unblasted_experimental_list)
blasted_control_counts <- vroom(blasted_experimental_list)



dir_mpileups_experimental <- "/mnt/bulk_data/THESIS/DATA_THESIS/GILOT_DATA_THESIS/EXPERIMENTAL/counts"
dir_mpileups <- dir_mpileups_experimental
unblasted_experimental <- "*_unblasted_gilot_experimental_counts.tsv"
blasted_experimental <- "*_blasted_gilot_experimental_counts.tsv"
unblasted_experimental_list <- list.files(dir_mpileups, full.names = TRUE,  pattern = unblasted_experimental)
blasted_experimental_list <- list.files(dir_mpileups, full.names = TRUE,  pattern = blasted_experimental)

dir_mpileups_control <- "/mnt/bulk_data/THESIS/DATA_THESIS/GILOT_DATA_THESIS/CONTROL/counts"
dir_mpileups <- dir_mpileups_control
unblasted_control <- "*_unblasted_gilot_all_controls_counts.tsv"
blasted_control <- "*_blasted_gilot_all_controls_counts.tsv"
unblasted_control_list <- list.files(dir_mpileups, full.names = TRUE,  pattern = unblasted_control)
blasted_control_list <- list.files(dir_mpileups, full.names = TRUE,  pattern = blasted_control)

sample_depth_statistics <- data.frame(STAT_IDENTIFIER=character(0), STAT_SAMPLE_METHOD=character(0),STAT_DEPTH_MIN=numeric(0), STAT_DEPTH_QUANT_1=numeric(0), STAT_DEPTH_MEDIAN=numeric(0), STAT_DEPTH_MEAN=numeric(0), STAT_DEPTH_QUANT_3=numeric(0), STAT_DEPTH_MAX=numeric(0), STAT_DEPTH_SD=numeric(0), STAT_GINI_MIN=numeric(0), STAT_GINI_QUANT_1=numeric(0), STAT_GINI_MEDIAN=numeric(0), STAT_GINI_MEAN=numeric(0), STAT_GINI_QUANT_3=numeric(0), STAT_GINI_MAX=numeric(0), STAT_GINI_SD=numeric(0)) %>% tibble()

counts_full_list <- c(unblasted_experimental_list, blasted_experimental_list, unblasted_control_list, blasted_control_list) %>% unique() %>% sort()
# counts_full_list <- c(unblasted_experimental_list, blasted_experimental_list) %>% unique() %>% sort()
# counts_full_list <- c(unblasted_experimental_list, blasted_experimental_list) %>% unique() %>% sort()
dir_save <- "/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/counts_from_mpileups"
for (filename_count in sample(counts_full_list)) { 
  print(filename_count)
  df_all_rows <- vroom(filename_count)
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
  colnames(df_all_rows) <- c("SAMPLE_GROUP", "SAMPLE_METHOD", "SAMPLE_NAME", "CHR", "POS", "DEPTH", "REFERENCE", "REFERENCE_COUNT", "ALTERNATE_COUNT", "A", "C", "G", "T", "N", "INDEL", "IDENTIFIER", "GINI_INDEX")
  df_all_rows$FREQUENCY_REFERENCE <- df_all_rows$REFERENCE_COUNT / df_all_rows$DEPTH
  df_all_rows$FREQUENCY_ALTERNATE <- df_all_rows$ALTERNATE_COUNT / df_all_rows$DEPTH
  
  df_all_rows$FREQUENCY_A <- df_all_rows$A / df_all_rows$DEPTH
  df_all_rows$FREQUENCY_C <- df_all_rows$C / df_all_rows$DEPTH
  df_all_rows$FREQUENCY_G <- df_all_rows$G / df_all_rows$DEPTH
  df_all_rows$FREQUENCY_T <- df_all_rows$T / df_all_rows$DEPTH
  df_all_rows$FREQUENCY_N <- df_all_rows$N / df_all_rows$DEPTH
  df_all_rows$FREQUENCY_INDEL <- df_all_rows$INDEL / df_all_rows$DEPTH
  df_all_rows$FREQUENCY_BASE_MAX <- colnames(df_all_rows[,10:15])[max.col(df_all_rows[,10:15],ties.method="first")]
  
  STAT_IDENTIFIER <- df_all_rows$IDENTIFIER %>% unique()
  STAT_SAMPLE_METHOD <- df_all_rows$SAMPLE_METHOD %>% unique()
  
  STAT_DEPTH_MIN <- summary(df_all_rows$DEPTH)[[1]]
  STAT_DEPTH_QUANT_1 <- summary(df_all_rows$DEPTH)[[2]]
  STAT_DEPTH_MEDIAN <- summary(df_all_rows$DEPTH)[[3]]
  STAT_DEPTH_MEAN <- summary(df_all_rows$DEPTH)[[4]]
  STAT_DEPTH_QUANT_3 <- summary(df_all_rows$DEPTH)[[5]]
  STAT_DEPTH_MAX <- summary(df_all_rows$DEPTH)[[6]]
  STAT_DEPTH_SD <- sd(df_all_rows$DEPTH)
  
  STAT_GINI_MIN <- summary(df_all_rows$GINI_INDEX)[[1]]
  STAT_GINI_QUANT_1 <- summary(df_all_rows$GINI_INDEX)[[2]]
  STAT_GINI_MEDIAN <- summary(df_all_rows$GINI_INDEX)[[3]]
  STAT_GINI_MEAN <- summary(df_all_rows$GINI_INDEX)[[4]]
  STAT_GINI_QUANT_3 <- summary(df_all_rows$GINI_INDEX)[[5]]
  STAT_GINI_MAX <- summary(df_all_rows$GINI_INDEX)[[6]]
  STAT_GINI_SD <- sd(df_all_rows$GINI_INDEX)
  
  sample_depth_statistics[nrow(sample_depth_statistics) + 1,] <- list(STAT_IDENTIFIER, STAT_SAMPLE_METHOD, STAT_DEPTH_MIN, STAT_DEPTH_QUANT_1, STAT_DEPTH_MEDIAN, STAT_DEPTH_MEAN, STAT_DEPTH_QUANT_3, STAT_DEPTH_MAX, STAT_DEPTH_SD, STAT_GINI_MIN, STAT_GINI_QUANT_1, STAT_GINI_MEDIAN, STAT_GINI_MEAN, STAT_GINI_QUANT_3, STAT_GINI_MAX, STAT_GINI_SD)
  
  # colnames(DF)[apply(DF,1,which.max)]
  
  new_filename <- paste0(dir_save, "/", (basename(filename_count) %>% tools::file_path_sans_ext()),"_gini.csv")
  print(new_filename)
  data.table::fwrite(df_all_rows, file = new_filename, append = FALSE, quote = "auto", showProgress = getOption("datatable.showProgress", interactive()), buffMB = 1024, nThread = 6)
}

new_filename <- paste0(dir_save, "/", "sample_depth_statistics_from_mpileups_gini.csv")
print(new_filename)
data.table::fwrite(sample_depth_statistics, file = new_filename, append = FALSE, quote = "auto", showProgress = getOption("datatable.showProgress", interactive()), buffMB = 1024, nThread = 6)

new_filename <- paste0("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/", "sample_depth_statistics_from_mpileups_gini.csv")
print(new_filename)
data.table::fwrite(sample_depth_statistics, file = new_filename, append = FALSE, quote = "auto", showProgress = getOption("datatable.showProgress", interactive()), buffMB = 1024, nThread = 6)

##########################
##########################     BACKGROUND AND CONSENSUS
##########################
########################################################################################
rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() #free up memrory and report the memory usage.

library(tidyverse)
library(vroom)
library(lubridate)

all_variant_data <- readRDS("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/all_callers_with_features_and_metadata.rds") %>% filter(FILTER == "PASS", THESIS_DEPTH_TOTAL >= 100)
min(all_variant_data$THESIS_DEPTH_TOTAL)
positions_to_exclude <- all_variant_data$POS %>% unique() %>% sort()

dir_mpileups_experimental <- "/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/counts_from_mpileups"
dir_mpileups <- dir_mpileups_experimental
unblasted_experimental <- "*_unblasted_gilot_experimental_counts_gini.csv"
blasted_experimental <- "*_blasted_gilot_experimental_counts_gini.csv"
unblasted_experimental_list <- list.files(dir_mpileups, full.names = TRUE,  pattern = unblasted_experimental)
blasted_experimental_list <- list.files(dir_mpileups, full.names = TRUE,  pattern = blasted_experimental)
# unblasted_experimental_counts <- vroom(unblasted_experimental_list)
# blasted_experimental_counts <- vroom(blasted_experimental_list)

dir_mpileups_control <- "/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/counts_from_mpileups"
dir_mpileups <- dir_mpileups_control
unblasted_control <- "*_unblasted_gilot_all_controls_counts_gini.csv"
blasted_control <- "*_blasted_gilot_all_controls_counts_gini.csv"
unblasted_control_list <- list.files(dir_mpileups, full.names = TRUE,  pattern = unblasted_control)
blasted_control_list <- list.files(dir_mpileups, full.names = TRUE,  pattern = blasted_control)
# unblasted_control_counts <- vroom(unblasted_experimental_list)
# blasted_control_counts <- vroom(blasted_experimental_list)


# library(ggplot2)
# all_variant_data <- readRDS("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/all_callers_with_features_and_metadata.rds") %>% filter(FILTER == "PASS")
# 
# # Histogram with density plot
# ggplot(all_variant_data%>% filter(THESIS_DEPTH_TOTAL >= 100,  THESIS_DEPTH_TOTAL <= 10000), aes(x=THESIS_DEPTH_TOTAL)) + 
#   geom_histogram(aes(y=..density..), colour="black", fill="white")+
#   geom_density(alpha=.2, fill="#FF6666") 


background_and_consensus <- data.frame(background_blaster=character(0), background_cohort=character(0), background_filename=character(0), background_identifier=character(0), base_ref_a_read_a_fraction=numeric(0), base_ref_a_read_c_fraction=numeric(0), base_ref_a_read_count_a=numeric(0), base_ref_a_read_count_c=numeric(0), base_ref_a_read_count_g=numeric(0), base_ref_a_read_count_indel=numeric(0), base_ref_a_read_count_n=numeric(0), base_ref_a_read_count_t=numeric(0), base_ref_a_read_count_total=numeric(0), base_ref_a_read_g_fraction=numeric(0), base_ref_a_read_indel_fraction=numeric(0), base_ref_a_read_n_fraction=numeric(0), base_ref_a_read_t_fraction=numeric(0), base_ref_c_read_a_fraction=numeric(0), base_ref_c_read_c_fraction=numeric(0), base_ref_c_read_count_a=numeric(0), base_ref_c_read_count_c=numeric(0), base_ref_c_read_count_g=numeric(0), base_ref_c_read_count_indel=numeric(0), base_ref_c_read_count_n=numeric(0), base_ref_c_read_count_t=numeric(0), base_ref_c_read_count_total=numeric(0), base_ref_c_read_g_fraction=numeric(0), base_ref_c_read_indel_fraction=numeric(0), base_ref_c_read_n_fraction=numeric(0), base_ref_c_read_t_fraction=numeric(0), base_ref_g_read_a_fraction=numeric(0), base_ref_g_read_c_fraction=numeric(0), base_ref_g_read_count_a=numeric(0), base_ref_g_read_count_c=numeric(0), base_ref_g_read_count_g=numeric(0), base_ref_g_read_count_indel=numeric(0), base_ref_g_read_count_n=numeric(0), base_ref_g_read_count_t=numeric(0), base_ref_g_read_count_total=numeric(0), base_ref_g_read_g_fraction=numeric(0), base_ref_g_read_indel_fraction=numeric(0), base_ref_g_read_n_fraction=numeric(0), base_ref_g_read_t_fraction=numeric(0), base_ref_t_read_a_fraction=numeric(0), base_ref_t_read_c_fraction=numeric(0), base_ref_t_read_count_a=numeric(0), base_ref_t_read_count_c=numeric(0), base_ref_t_read_count_g=numeric(0), base_ref_t_read_count_indel=numeric(0), base_ref_t_read_count_n=numeric(0), base_ref_t_read_count_t=numeric(0), base_ref_t_read_count_total=numeric(0), base_ref_t_read_g_fraction=numeric(0), base_ref_t_read_indel_fraction=numeric(0), base_ref_t_read_n_fraction=numeric(0), base_ref_t_read_t_fraction=numeric(0), consensus_pileup_genome=character(0), consensus_pileup_genome_length=numeric(0)) %>% tibble()

# x_temp <- background_and_consensus
# single_temp <- background_and_consensus

counts_full_list <- c(unblasted_experimental_list, blasted_experimental_list, unblasted_control_list, blasted_control_list) %>% unique() %>% sort()
dir_save <- "/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/counts_from_mpileups"
for (filename_count in sample(counts_full_list)) { 
  print(filename_count)
  background_filename <- filename_count
  # 99% Reference, and Alternate not 0
  tibble_background_raw <- vroom(background_filename)

  tibble_background <- tibble_background_raw %>% filter(DEPTH >= 100, REFERENCE == FREQUENCY_BASE_MAX, FREQUENCY_REFERENCE >= 0.995, ALTERNATE_COUNT != 0) %>% filter(!(POS %in% positions_to_exclude)) %>% arrange(POS)
  
  tibble_background_raw$FREQUENCY_BASE_MAX[tibble_background_raw$FREQUENCY_BASE_MAX == 'INDEL'] <- '-'
  consensus_pileup_genome <- paste(tibble_background_raw$FREQUENCY_BASE_MAX, collapse = '')
  consensus_pileup_genome_length <- consensus_pileup_genome %>% nchar()
  
  background_cohort <- tibble_background$SAMPLE_GROUP %>% unique()
  background_blaster <- tibble_background$SAMPLE_METHOD %>% unique()
  background_identifier <- tibble_background$IDENTIFIER %>% unique()
  
  mpiles_base <- tibble_background %>% filter(FREQUENCY_BASE_MAX == "A")
  base_ref_a_read_count_a <- sum(mpiles_base$A)
  base_ref_a_read_count_c <- sum(mpiles_base$C)
  base_ref_a_read_count_g <- sum(mpiles_base$G)
  base_ref_a_read_count_t <- sum(mpiles_base$T)
  base_ref_a_read_count_n <- sum(mpiles_base$N)
  base_ref_a_read_count_indel <- sum(mpiles_base$INDEL)
  base_ref_a_read_count_total <- base_ref_a_read_count_a + base_ref_a_read_count_c + base_ref_a_read_count_g + base_ref_a_read_count_t + base_ref_a_read_count_n + base_ref_a_read_count_indel
  base_ref_a_read_a_fraction <- base_ref_a_read_count_a /  base_ref_a_read_count_total
  base_ref_a_read_c_fraction <- base_ref_a_read_count_c /  base_ref_a_read_count_total
  base_ref_a_read_g_fraction <- base_ref_a_read_count_g /  base_ref_a_read_count_total
  base_ref_a_read_t_fraction <- base_ref_a_read_count_t /  base_ref_a_read_count_total
  base_ref_a_read_n_fraction <- base_ref_a_read_count_n /  base_ref_a_read_count_total
  base_ref_a_read_indel_fraction <- base_ref_a_read_count_indel /  base_ref_a_read_count_total
  
  mpiles_base <- tibble_background %>% filter(FREQUENCY_BASE_MAX == "C")
  base_ref_c_read_count_a <- sum(mpiles_base$A)
  base_ref_c_read_count_c <- sum(mpiles_base$C)
  base_ref_c_read_count_g <- sum(mpiles_base$G)
  base_ref_c_read_count_t <- sum(mpiles_base$T)
  base_ref_c_read_count_n <- sum(mpiles_base$N)
  base_ref_c_read_count_indel <- sum(mpiles_base$INDEL)
  base_ref_c_read_count_total <- base_ref_c_read_count_a + base_ref_c_read_count_c + base_ref_c_read_count_g + base_ref_c_read_count_t + base_ref_c_read_count_n + base_ref_c_read_count_indel
  base_ref_c_read_a_fraction <- base_ref_c_read_count_a /  base_ref_c_read_count_total
  base_ref_c_read_c_fraction <- base_ref_c_read_count_c /  base_ref_c_read_count_total
  base_ref_c_read_g_fraction <- base_ref_c_read_count_g /  base_ref_c_read_count_total
  base_ref_c_read_t_fraction <- base_ref_c_read_count_t /  base_ref_c_read_count_total
  base_ref_c_read_n_fraction <- base_ref_c_read_count_n /  base_ref_c_read_count_total
  base_ref_c_read_indel_fraction <- base_ref_c_read_count_indel /  base_ref_c_read_count_total
  
  mpiles_base <- tibble_background %>% filter(FREQUENCY_BASE_MAX == "G")
  base_ref_g_read_count_a <- sum(mpiles_base$A)
  base_ref_g_read_count_c <- sum(mpiles_base$C)
  base_ref_g_read_count_g <- sum(mpiles_base$G)
  base_ref_g_read_count_t <- sum(mpiles_base$T)
  base_ref_g_read_count_n <- sum(mpiles_base$N)
  base_ref_g_read_count_indel <- sum(mpiles_base$INDEL)
  base_ref_g_read_count_total <- base_ref_g_read_count_a + base_ref_g_read_count_c + base_ref_g_read_count_g + base_ref_g_read_count_t + base_ref_g_read_count_n + base_ref_g_read_count_indel
  base_ref_g_read_a_fraction <- base_ref_g_read_count_a /  base_ref_g_read_count_total
  base_ref_g_read_c_fraction <- base_ref_g_read_count_c /  base_ref_g_read_count_total
  base_ref_g_read_g_fraction <- base_ref_g_read_count_g /  base_ref_g_read_count_total
  base_ref_g_read_t_fraction <- base_ref_g_read_count_t /  base_ref_g_read_count_total
  base_ref_g_read_n_fraction <- base_ref_g_read_count_n /  base_ref_g_read_count_total
  base_ref_g_read_indel_fraction <- base_ref_g_read_count_indel /  base_ref_g_read_count_total
  
  mpiles_base <- tibble_background %>% filter(FREQUENCY_BASE_MAX == "T")
  base_ref_t_read_count_a <- sum(mpiles_base$A)
  base_ref_t_read_count_c <- sum(mpiles_base$C)
  base_ref_t_read_count_g <- sum(mpiles_base$G)
  base_ref_t_read_count_t <- sum(mpiles_base$T)
  base_ref_t_read_count_n <- sum(mpiles_base$N)
  base_ref_t_read_count_indel <- sum(mpiles_base$INDEL)
  base_ref_t_read_count_total <- base_ref_t_read_count_a + base_ref_t_read_count_c + base_ref_t_read_count_g + base_ref_t_read_count_t + base_ref_t_read_count_n + base_ref_t_read_count_indel
  base_ref_t_read_a_fraction <- base_ref_t_read_count_a /  base_ref_t_read_count_total
  base_ref_t_read_c_fraction <- base_ref_t_read_count_c /  base_ref_t_read_count_total
  base_ref_t_read_g_fraction <- base_ref_t_read_count_g /  base_ref_t_read_count_total
  base_ref_t_read_t_fraction <- base_ref_t_read_count_t /  base_ref_t_read_count_total
  base_ref_t_read_n_fraction <- base_ref_t_read_count_n /  base_ref_t_read_count_total
  base_ref_t_read_indel_fraction <- base_ref_t_read_count_indel /  base_ref_t_read_count_total
  
  print("DATA LOADED FOR NEW ROW")
  print(background_filename)
  print("BEGIN BIND ROWS")
  
  # rapply( z, f=function(x) ifelse(is.nan(x),0,x), how="replace" )
  background_and_consensus[nrow(background_and_consensus) + 1,] <- rapply( list(background_blaster, background_cohort, background_filename, background_identifier, base_ref_a_read_a_fraction, base_ref_a_read_c_fraction, base_ref_a_read_count_a, base_ref_a_read_count_c, base_ref_a_read_count_g, base_ref_a_read_count_indel, base_ref_a_read_count_n, base_ref_a_read_count_t, base_ref_a_read_count_total, base_ref_a_read_g_fraction , base_ref_a_read_indel_fraction, base_ref_a_read_n_fraction, base_ref_a_read_t_fraction, base_ref_c_read_a_fraction, base_ref_c_read_c_fraction, base_ref_c_read_count_a, base_ref_c_read_count_c, base_ref_c_read_count_g, base_ref_c_read_count_indel, base_ref_c_read_count_n, base_ref_c_read_count_t, base_ref_c_read_count_total, base_ref_c_read_g_fraction, base_ref_c_read_indel_fraction, base_ref_c_read_n_fraction, base_ref_c_read_t_fraction, base_ref_g_read_a_fraction, base_ref_g_read_c_fraction, base_ref_g_read_count_a, base_ref_g_read_count_c, base_ref_g_read_count_g, base_ref_g_read_count_indel, base_ref_g_read_count_n, base_ref_g_read_count_t, base_ref_g_read_count_total, base_ref_g_read_g_fraction, base_ref_g_read_indel_fraction, base_ref_g_read_n_fraction, base_ref_g_read_t_fraction, base_ref_t_read_a_fraction, base_ref_t_read_c_fraction, base_ref_t_read_count_a, base_ref_t_read_count_c, base_ref_t_read_count_g, base_ref_t_read_count_indel, base_ref_t_read_count_n, base_ref_t_read_count_t, base_ref_t_read_count_total, base_ref_t_read_g_fraction, base_ref_t_read_indel_fraction, base_ref_t_read_n_fraction, base_ref_t_read_t_fraction, consensus_pileup_genome, consensus_pileup_genome_length), f=function(x) ifelse(is.nan(x),0,x), how="replace" )
  
 # background_and_consensus[nrow(background_and_consensus) + 1,] <- rapply( nonan, f=function(x) ifelse(is.na(x),"exclude_sample",x), how="replace" )
  
 #  background_and_consensus[nrow(background_and_consensus) + 1,] <- list(background_blaster, background_cohort, background_filename, background_identifier, base_ref_a_read_a_fraction, base_ref_a_read_c_fraction, base_ref_a_read_count_a, base_ref_a_read_count_c, base_ref_a_read_count_g, base_ref_a_read_count_indel, base_ref_a_read_count_n, base_ref_a_read_count_t, base_ref_a_read_count_total, base_ref_a_read_g_fraction , base_ref_a_read_indel_fraction, base_ref_a_read_n_fraction, base_ref_a_read_t_fraction, base_ref_c_read_a_fraction, base_ref_c_read_c_fraction, base_ref_c_read_count_a, base_ref_c_read_count_c, base_ref_c_read_count_g, base_ref_c_read_count_indel, base_ref_c_read_count_n, base_ref_c_read_count_t, base_ref_c_read_count_total, base_ref_c_read_g_fraction, base_ref_c_read_indel_fraction, base_ref_c_read_n_fraction, base_ref_c_read_t_fraction, base_ref_g_read_a_fraction, base_ref_g_read_c_fraction, base_ref_g_read_count_a, base_ref_g_read_count_c, base_ref_g_read_count_g, base_ref_g_read_count_indel, base_ref_g_read_count_n, base_ref_g_read_count_t, base_ref_g_read_count_total, base_ref_g_read_g_fraction, base_ref_g_read_indel_fraction, base_ref_g_read_n_fraction, base_ref_g_read_t_fraction, base_ref_t_read_a_fraction, base_ref_t_read_c_fraction, base_ref_t_read_count_a, base_ref_t_read_count_c, base_ref_t_read_count_g, base_ref_t_read_count_indel, base_ref_t_read_count_n, base_ref_t_read_count_t, base_ref_t_read_count_total, base_ref_t_read_g_fraction, base_ref_t_read_indel_fraction, base_ref_t_read_n_fraction, base_ref_t_read_t_fraction, consensus_pileup_genome, consensus_pileup_genome_length)
  
#  x_temp[nrow(x_temp) + 1,] <- list(background_blaster, background_cohort, background_filename, background_identifier, base_ref_a_read_a_fraction, base_ref_a_read_c_fraction, base_ref_a_read_count_a, base_ref_a_read_count_c, base_ref_a_read_count_g, base_ref_a_read_count_indel, base_ref_a_read_count_n, base_ref_a_read_count_t, base_ref_a_read_count_total, base_ref_a_read_g_fraction , base_ref_a_read_indel_fraction, base_ref_a_read_n_fraction, base_ref_a_read_t_fraction, base_ref_c_read_a_fraction, base_ref_c_read_c_fraction, base_ref_c_read_count_a, base_ref_c_read_count_c, base_ref_c_read_count_g, base_ref_c_read_count_indel, base_ref_c_read_count_n, base_ref_c_read_count_t, base_ref_c_read_count_total, base_ref_c_read_g_fraction, base_ref_c_read_indel_fraction, base_ref_c_read_n_fraction, base_ref_c_read_t_fraction, base_ref_g_read_a_fraction, base_ref_g_read_c_fraction, base_ref_g_read_count_a, base_ref_g_read_count_c, base_ref_g_read_count_g, base_ref_g_read_count_indel, base_ref_g_read_count_n, base_ref_g_read_count_t, base_ref_g_read_count_total, base_ref_g_read_g_fraction, base_ref_g_read_indel_fraction, base_ref_g_read_n_fraction, base_ref_g_read_t_fraction, base_ref_t_read_a_fraction, base_ref_t_read_c_fraction, base_ref_t_read_count_a, base_ref_t_read_count_c, base_ref_t_read_count_g, base_ref_t_read_count_indel, base_ref_t_read_count_n, base_ref_t_read_count_t, base_ref_t_read_count_total, base_ref_t_read_g_fraction, base_ref_t_read_indel_fraction, base_ref_t_read_n_fraction, base_ref_t_read_t_fraction, consensus_pileup_genome, consensus_pileup_genome_length)
  
  # print(background_filename)
  # print(x_temp)
  # 
}

background_and_consensus <- background_and_consensus %>% arrange(background_filename)
view(background_and_consensus)
background_and_consensus$background_blaster
background_and_consensus$background_blaster <- background_and_consensus$background_blaster %>% replace_na("unknown")

view(background_and_consensus)
background_and_consensus$background_cohort
background_and_consensus$background_cohort <- background_and_consensus$background_cohort %>% replace_na("unknown")

view(background_and_consensus)
background_and_consensus$background_identifier
background_and_consensus$background_identifier <- background_and_consensus$background_identifier %>% replace_na("unknown")

view(background_and_consensus)

#########    NOT SAVING EVERYTHING #########  FIXED by replacing NaN with NA in the lists with rapply
new_filename <- paste0(dir_save, "/", "background_and_consensus_for_marascuillo_procedure.csv")
print(new_filename)
data.table::fwrite(background_and_consensus, file = new_filename, append = FALSE, quote = "auto", showProgress = getOption("datatable.showProgress", interactive()), buffMB = 1024, nThread = 6)
#########    NOT SAVING EVERYTHING #########
new_filename <- paste0("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/", "background_and_consensus_for_marascuillo_procedure.csv")
print(new_filename)
data.table::fwrite(background_and_consensus, file = new_filename, append = FALSE, quote = "auto", showProgress = getOption("datatable.showProgress", interactive()), buffMB = 1024, nThread = 6)

metadata_goldmaster <- readRDS("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/metadata_goldmaster.rds")
background_and_consensus_labelled <- left_join(background_and_consensus, metadata_goldmaster, by = c('background_identifier' = 'THESIS_SAMPLE_IDENTIFIER'))
# left_join(background_and_consensus, metadata_goldmaster, by = c('background_identifier' = 'THESIS_SAMPLE_IDENTIFIER'), multiple = "first" )

#########    NOT SAVING EVERYTHING #########  FIXED by replacing NaN with NA in the lists with rapply
new_filename <- paste0(dir_save, "/", "background_and_consensus_for_marascuillo_procedure_labelled.csv")
print(new_filename)
data.table::fwrite(background_and_consensus_labelled, file = new_filename, append = FALSE, quote = "auto", showProgress = getOption("datatable.showProgress", interactive()), buffMB = 1024, nThread = 6)
#########    NOT SAVING EVERYTHING #########
new_filename <- paste0("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/", "background_and_consensus_for_marascuillo_procedure_labelled.csv")
print(new_filename)
data.table::fwrite(background_and_consensus_labelled, file = new_filename, append = FALSE, quote = "auto", showProgress = getOption("datatable.showProgress", interactive()), buffMB = 1024, nThread = 6)

view(background_and_consensus_labelled)


background_and_consensus$base_ref_a_read_a_fraction %>% median()
background_and_consensus$base_ref_a_read_c_fraction %>% median()
background_and_consensus$base_ref_a_read_g_fraction %>% median()
background_and_consensus$base_ref_a_read_t_fraction %>% median()
background_and_consensus$base_ref_a_read_n_fraction %>% median()
background_and_consensus$base_ref_a_read_indel_fraction %>% median()

background_and_consensus$base_ref_a_read_a_fraction %>% mean()
background_and_consensus$base_ref_a_read_c_fraction %>% mean()
background_and_consensus$base_ref_a_read_g_fraction %>% mean()
background_and_consensus$base_ref_a_read_t_fraction %>% mean()
background_and_consensus$base_ref_a_read_n_fraction %>% mean()
background_and_consensus$base_ref_a_read_indel_fraction %>% mean()

background_and_consensus$base_ref_a_read_a_fraction %>% sd()
background_and_consensus$base_ref_a_read_c_fraction %>% sd()
background_and_consensus$base_ref_a_read_g_fraction %>% sd()
background_and_consensus$base_ref_a_read_t_fraction %>% sd()
background_and_consensus$base_ref_a_read_n_fraction %>% sd()
background_and_consensus$base_ref_a_read_indel_fraction %>% sd()

background_and_consensus$base_ref_a_read_a_fraction %>% summary()
background_and_consensus$base_ref_a_read_c_fraction %>% summary()
background_and_consensus$base_ref_a_read_g_fraction %>% summary()
background_and_consensus$base_ref_a_read_t_fraction %>% summary()
background_and_consensus$base_ref_a_read_n_fraction %>% summary()
background_and_consensus$base_ref_a_read_indel_fraction %>% summary()

sum(background_and_consensus$base_ref_a_read_count_a ) / sum(background_and_consensus$base_ref_a_read_count_total)
sum(background_and_consensus$base_ref_a_read_count_c ) / sum(background_and_consensus$base_ref_a_read_count_total)
sum(background_and_consensus$base_ref_a_read_count_g ) / sum(background_and_consensus$base_ref_a_read_count_total)
sum(background_and_consensus$base_ref_a_read_count_t ) / sum(background_and_consensus$base_ref_a_read_count_total)
sum(background_and_consensus$base_ref_a_read_count_n ) / sum(background_and_consensus$base_ref_a_read_count_total)
sum(background_and_consensus$base_ref_a_read_count_indel ) / sum(background_and_consensus$base_ref_a_read_count_total)


##########################
##########################     VARIANT SIGNATURES
##########################
########################################################################################
rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() #free up memrory and report the memory usage.

library(tidyverse)
library(vroom)
library(lubridate)


all_variant_data <- readRDS("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/all_callers_with_features_and_metadata.rds") %>% filter(FILTER == "PASS", THESIS_DEPTH_TOTAL >= 100)
min(all_variant_data$THESIS_DEPTH_TOTAL)
positions_to_exclude <- all_variant_data$POS %>% unique() %>% sort()
length(positions_to_exclude)
colnames(all_variant_data)
#view(all_variant_data)
#colnames(all_variant_data %>% select(contains("_COPIES")))
#colnames(all_variant_data %>% select(contains("THESIS_SYNTHETIC_PRESENT_C")))

THESIS_SYNTHETIC_PRESENT_C1_AUSTRALIA_MT007544_1 <- all_variant_data %>% filter(THESIS_SYNTHETIC_PRESENT_C1_AUSTRALIA_MT007544_1 == TRUE) %>% select(THESIS_IDENTIFIER, contains("_SYNTHETIC_PRESENT_C1_AUSTRALIA_MT007544_1")) %>% arrange (THESIS_SYNTHETIC_PRESENT_C1_AUSTRALIA_MT007544_1_COPIES, THESIS_SYNTHETIC_PRESENT_C1_AUSTRALIA_MT007544_1_FRACTION) %>% unique() 
data_fwrite <- THESIS_SYNTHETIC_PRESENT_C1_AUSTRALIA_MT007544_1
name_data_frwite <- "THESIS_SYNTHETIC_PRESENT_C1_AUSTRALIA_MT007544_1"
new_filename <- paste0("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/", name_data_frwite, ".csv")
print(new_filename)
data.table::fwrite(data_fwrite, file = new_filename, append = FALSE, quote = "auto", showProgress = getOption("datatable.showProgress", interactive()), buffMB = 1024, nThread = 6)

THESIS_SYNTHETIC_PRESENT_C2_ANCESTRAL_MN908947_3 <- all_variant_data %>% filter(THESIS_SYNTHETIC_PRESENT_C2_ANCESTRAL_MN908947_3 == TRUE) %>% select(THESIS_IDENTIFIER, contains("_SYNTHETIC_PRESENT_C2_ANCESTRAL_MN908947_3")) %>% arrange (THESIS_SYNTHETIC_PRESENT_C2_ANCESTRAL_MN908947_3_COPIES, THESIS_SYNTHETIC_PRESENT_C2_ANCESTRAL_MN908947_3_FRACTION) %>% unique() 
data_fwrite <- THESIS_SYNTHETIC_PRESENT_C2_ANCESTRAL_MN908947_3
name_data_frwite <- "THESIS_SYNTHETIC_PRESENT_C2_ANCESTRAL_MN908947_3"
new_filename <- paste0("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/", name_data_frwite, ".csv")
print(new_filename)
data.table::fwrite(data_fwrite, file = new_filename, append = FALSE, quote = "auto", showProgress = getOption("datatable.showProgress", interactive()), buffMB = 1024, nThread = 6)

THESIS_SYNTHETIC_PRESENT_C14_ALPHA_EPI_ISL_710528 <- all_variant_data %>% filter(THESIS_SYNTHETIC_PRESENT_C14_ALPHA_EPI_ISL_710528 == TRUE) %>% select(THESIS_IDENTIFIER, contains("_SYNTHETIC_PRESENT_C14_ALPHA_EPI_ISL_710528")) %>% arrange (THESIS_SYNTHETIC_PRESENT_C14_ALPHA_EPI_ISL_710528_COPIES, THESIS_SYNTHETIC_PRESENT_C14_ALPHA_EPI_ISL_710528_FRACTION) %>% unique()
data_fwrite <- THESIS_SYNTHETIC_PRESENT_C14_ALPHA_EPI_ISL_710528
name_data_frwite <- "THESIS_SYNTHETIC_PRESENT_C14_ALPHA_EPI_ISL_710528"
new_filename <- paste0("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/", name_data_frwite, ".csv")
print(new_filename)
data.table::fwrite(data_fwrite, file = new_filename, append = FALSE, quote = "auto", showProgress = getOption("datatable.showProgress", interactive()), buffMB = 1024, nThread = 6)

THESIS_SYNTHETIC_PRESENT_C29_DELTA_EPI_ISL_2693246 <- all_variant_data %>% filter(THESIS_SYNTHETIC_PRESENT_C29_DELTA_EPI_ISL_2693246 == TRUE) %>% select(THESIS_IDENTIFIER, contains("_SYNTHETIC_PRESENT_C29_DELTA_EPI_ISL_2693246")) %>% arrange (THESIS_SYNTHETIC_PRESENT_C29_DELTA_EPI_ISL_2693246_COPIES, THESIS_SYNTHETIC_PRESENT_C29_DELTA_EPI_ISL_2693246_FRACTION) %>% unique()
data_fwrite <- THESIS_SYNTHETIC_PRESENT_C29_DELTA_EPI_ISL_2693246
name_data_frwite <- "THESIS_SYNTHETIC_PRESENT_C29_DELTA_EPI_ISL_2693246"
new_filename <- paste0("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/", name_data_frwite, ".csv")
print(new_filename)
data.table::fwrite(data_fwrite, file = new_filename, append = FALSE, quote = "auto", showProgress = getOption("datatable.showProgress", interactive()), buffMB = 1024, nThread = 6)

THESIS_SYNTHETIC_PRESENT_C48_OMICRON_EPI_ISL_6841980 <- all_variant_data %>% filter(THESIS_SYNTHETIC_PRESENT_C48_OMICRON_EPI_ISL_6841980 == TRUE) %>% select(THESIS_IDENTIFIER, contains("_SYNTHETIC_PRESENT_C48_OMICRON_EPI_ISL_6841980")) %>% arrange (THESIS_SYNTHETIC_PRESENT_C48_OMICRON_EPI_ISL_6841980_COPIES, THESIS_SYNTHETIC_PRESENT_C48_OMICRON_EPI_ISL_6841980_FRACTION) %>% unique()
data_fwrite <- THESIS_SYNTHETIC_PRESENT_C48_OMICRON_EPI_ISL_6841980
name_data_frwite <- "THESIS_SYNTHETIC_PRESENT_C48_OMICRON_EPI_ISL_6841980"
new_filename <- paste0("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/", name_data_frwite, ".csv")
print(new_filename)
data.table::fwrite(data_fwrite, file = new_filename, append = FALSE, quote = "auto", showProgress = getOption("datatable.showProgress", interactive()), buffMB = 1024, nThread = 6)

synthetic_combinations_copies <- all_variant_data %>% filter(THESIS_SYNTHETIC_COPIES_TOTAL != 0) %>% arrange(THESIS_SYNTHETIC_COPIES_TOTAL) %>% select(THESIS_IDENTIFIER, THESIS_SYNTHETIC_COPIES_TOTAL, contains("_SYNTHETIC_PRESENT_C")) %>% select(!ends_with("_COPIES"))  %>% select(!ends_with("_FRACTION")) %>% unique() %>% select(!ends_with("_FRACTION")) # %>% view()# %>% #!contains("_FRACTION"))
data_fwrite <- synthetic_combinations_copies
name_data_frwite <- "synthetic_combinations_copies"
new_filename <- paste0("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/", name_data_frwite, ".csv")
print(new_filename)
data.table::fwrite(data_fwrite, file = new_filename, append = FALSE, quote = "auto", showProgress = getOption("datatable.showProgress", interactive()), buffMB = 1024, nThread = 6)

synthetic_combinations_fractions <- all_variant_data %>% filter(THESIS_SYNTHETIC_COPIES_TOTAL != 0) %>% arrange(THESIS_SYNTHETIC_COPIES_TOTAL) %>% select(THESIS_IDENTIFIER, THESIS_SYNTHETIC_COPIES_TOTAL, contains("_SYNTHETIC_PRESENT_C")) %>% select(!ends_with("_COPIES")) %>% unique()
data_fwrite <- synthetic_combinations_fractions
name_data_frwite <- "synthetic_combinations_fractions"
new_filename <- paste0("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/", name_data_frwite, ".csv")
print(new_filename)
data.table::fwrite(data_fwrite, file = new_filename, append = FALSE, quote = "auto", showProgress = getOption("datatable.showProgress", interactive()), buffMB = 1024, nThread = 6)

synthetic_combinations_all <- all_variant_data %>% filter(THESIS_SYNTHETIC_COPIES_TOTAL != 0) %>% arrange(THESIS_SYNTHETIC_COPIES_TOTAL) %>% select(THESIS_IDENTIFIER, THESIS_SYNTHETIC_COPIES_TOTAL, contains("_SYNTHETIC_PRESENT_C")) %>% unique()
data_fwrite <- synthetic_combinations_all
name_data_frwite <- "synthetic_combinations_all"
new_filename <- paste0("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/", name_data_frwite, ".csv")
print(new_filename)
data.table::fwrite(data_fwrite, file = new_filename, append = FALSE, quote = "auto", showProgress = getOption("datatable.showProgress", interactive()), buffMB = 1024, nThread = 6)

nonsynthetic_combinations_all <- all_variant_data %>% filter(THESIS_SYNTHETIC_COPIES_TOTAL == 0) %>% arrange(THESIS_SYNTHETIC_COPIES_TOTAL) %>% select(THESIS_IDENTIFIER, THESIS_SYNTHETIC_COPIES_TOTAL, contains("_SYNTHETIC_PRESENT_C")) %>% unique()
data_fwrite <- nonsynthetic_combinations_all
name_data_frwite <- "nonsynthetic_combinations_all"
new_filename <- paste0("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/", name_data_frwite, ".csv")
print(new_filename)
data.table::fwrite(data_fwrite, file = new_filename, append = FALSE, quote = "auto", showProgress = getOption("datatable.showProgress", interactive()), buffMB = 1024, nThread = 6)

nonsynthetic_combinations_all$THESIS_IDENTIFIER %>% unique()
synthetic_combinations_all$THESIS_IDENTIFIER %>% unique()

#readRDS(file = "/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/locate_snps/variants_data/snp_positions_all.rds")

rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
gc() #free up memrory and report the memory usage.

library(tidyverse)
library(vroom)
library(lubridate)

github_data_unified <- "/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/"

all_callers_with_features_and_metadata <- readRDS("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/all_callers_with_features_and_metadata.rds")

SNPs_nucmer_from_synthetic_controls <- readRDS("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/SNPs_nucmer_from_synthetic_controls.rds")

features_gff_MN908947_3 <- readRDS("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/gff_MN908947_3.rds")

colnames(all_callers_with_features_and_metadata)
colnames(SNPs_nucmer_from_synthetic_controls)
all_callers_with_features_and_metadata$THESIS_VARIANT %>% unique()
SNPs_nucmer_from_synthetic_controls$SIGNATURE_VARIANT

variants_observed_signatures <- (all_callers_with_features_and_metadata %>% filter(FILTER == "PASS", THESIS_DEPTH_TOTAL > 99))$THESIS_VARIANT %>% sort() %>% unique() 
variants_observed_positions <- (all_callers_with_features_and_metadata %>% filter(FILTER == "PASS", THESIS_DEPTH_TOTAL > 99))$POS %>% sort() %>% unique()

variants_controls_signatures <- SNPs_nucmer_from_synthetic_controls$SIGNATURE_VARIANT %>% sort() %>% unique()
variants_controls_positions <- SNPs_nucmer_from_synthetic_controls$SIGNATURE_POS %>% sort() %>% unique()

intersect(variants_observed_signatures, variants_controls_signatures) 

setdiff(variants_observed_signatures, variants_controls_signatures) 
setdiff(variants_controls_signatures, variants_observed_signatures) # in this not that

SNPs_nucmer_from_synthetic_controls$QUERY_TAG %>% sort() %>% unique()
SNPs_nucmer_from_synthetic_controls %>% select(QUERY_TAG, SIGNATURE_POS, SIGNATURE_VARIANT)

signature_EPI_ISL_2693246 <- SNPs_nucmer_from_synthetic_controls %>% filter(QUERY_TAG == "EPI_ISL_2693246") %>% arrange(SIGNATURE_POS) %>% select(QUERY_TAG, SIGNATURE_POS, REFERENCE, QUERY, SIGNATURE_VARIANT)
all_signature_EPI_ISL_2693246 <- signature_EPI_ISL_2693246$SIGNATURE_VARIANT
data_fwrite <- signature_EPI_ISL_2693246
name_data_frwite <- "signature_EPI_ISL_2693246"
new_filename <- paste0("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/", name_data_frwite, ".csv")
print(new_filename)
data.table::fwrite(data_fwrite, file = new_filename, append = FALSE, quote = "auto", showProgress = getOption("datatable.showProgress", interactive()), buffMB = 1024, nThread = 6)
saveRDS(signature_EPI_ISL_2693246, file = paste0("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/", "signature_EPI_ISL_2693246", ".rds"), compress = TRUE)

signature_EPI_ISL_6841980 <- SNPs_nucmer_from_synthetic_controls %>% filter(QUERY_TAG == "EPI_ISL_6841980") %>% arrange(SIGNATURE_POS) %>% select(QUERY_TAG, SIGNATURE_POS, REFERENCE, QUERY, SIGNATURE_VARIANT)
all_signature_EPI_ISL_6841980 <- signature_EPI_ISL_6841980$SIGNATURE_VARIANT
data_fwrite <- signature_EPI_ISL_6841980
name_data_frwite <- "signature_EPI_ISL_6841980"
new_filename <- paste0("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/", name_data_frwite, ".csv")
print(new_filename)
data.table::fwrite(data_fwrite, file = new_filename, append = FALSE, quote = "auto", showProgress = getOption("datatable.showProgress", interactive()), buffMB = 1024, nThread = 6)
saveRDS(signature_EPI_ISL_6841980, file = paste0("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/", "signature_EPI_ISL_6841980", ".rds"), compress = TRUE)

signature_EPI_ISL_710528 <- SNPs_nucmer_from_synthetic_controls %>% filter(QUERY_TAG == "EPI_ISL_710528") %>% arrange(SIGNATURE_POS) %>% select(QUERY_TAG, SIGNATURE_POS, REFERENCE, QUERY, SIGNATURE_VARIANT)
all_signature_EPI_ISL_710528 <- signature_EPI_ISL_710528$SIGNATURE_VARIANT
data_fwrite <- signature_EPI_ISL_710528
name_data_frwite <- "signature_EPI_ISL_710528"
new_filename <- paste0("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/", name_data_frwite, ".csv")
print(new_filename)
data.table::fwrite(data_fwrite, file = new_filename, append = FALSE, quote = "auto", showProgress = getOption("datatable.showProgress", interactive()), buffMB = 1024, nThread = 6)
saveRDS(signature_EPI_ISL_710528, file = paste0("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/", "signature_EPI_ISL_710528", ".rds"), compress = TRUE)

signature_MT007544_1 <- SNPs_nucmer_from_synthetic_controls %>% filter(QUERY_TAG == "MT007544.1") %>% arrange(SIGNATURE_POS) %>% select(QUERY_TAG, SIGNATURE_POS, REFERENCE, QUERY, SIGNATURE_VARIANT)
all_signature_MT007544_1 <- signature_MT007544_1$SIGNATURE_VARIANT
data_fwrite <- signature_MT007544_1
name_data_frwite <- "signature_MT007544_1"
new_filename <- paste0("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/", name_data_frwite, ".csv")
print(new_filename)
data.table::fwrite(data_fwrite, file = new_filename, append = FALSE, quote = "auto", showProgress = getOption("datatable.showProgress", interactive()), buffMB = 1024, nThread = 6)
saveRDS(signature_MT007544_1, file = paste0("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/", "signature_MT007544_1", ".rds"), compress = TRUE)

#union(all_signature_MT007544_1, all_signature_EPI_ISL_710528, all_signature_EPI_ISL_6841980, all_signature_EPI_ISL_2693246)
variant_stage_1 <- union(all_signature_MT007544_1, all_signature_EPI_ISL_710528) #, all_signature_EPI_ISL_6841980, all_signature_EPI_ISL_2693246)
variant_not_in <- union(variant_stage_1, all_signature_EPI_ISL_6841980) #, a
variant_is_in_signature_EPI_ISL_2693246 <- setdiff(all_signature_EPI_ISL_2693246, variant_not_in)  %>% unique() # in this not that
variant_is_in_signature_EPI_ISL_2693246 %>% length()
saveRDS(variant_is_in_signature_EPI_ISL_2693246, file = paste0("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/", "variant_is_in_signature_EPI_ISL_2693246", ".rds"), compress = TRUE)
variant_is_in_signature_EPI_ISL_2693246
signature_EPI_ISL_2693246$SIGNATURE_VARIANT

#union(all_signature_MT007544_1, all_signature_EPI_ISL_710528, all_signature_EPI_ISL_6841980, all_signature_EPI_ISL_2693246)
variant_stage_1 <- union(all_signature_EPI_ISL_2693246, all_signature_EPI_ISL_710528) #, all_signature_EPI_ISL_6841980, all_signature_EPI_ISL_2693246)
variant_not_in <- union(variant_stage_1, all_signature_EPI_ISL_6841980) #, a
variant_is_in_signature_MT007544_1 <- setdiff(all_signature_MT007544_1, variant_not_in)  %>% unique() # in this not that
variant_is_in_signature_MT007544_1 %>% length()
saveRDS(variant_is_in_signature_MT007544_1, file = paste0("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/", "variant_is_in_signature_MT007544_1", ".rds"), compress = TRUE)
variant_is_in_signature_MT007544_1
signature_MT007544_1$SIGNATURE_VARIANT

#union(all_signature_MT007544_1, all_signature_EPI_ISL_710528, all_signature_EPI_ISL_6841980, all_signature_EPI_ISL_2693246)
variant_stage_1 <- union(all_signature_EPI_ISL_2693246, all_signature_MT007544_1) #, all_signature_EPI_ISL_6841980, all_signature_EPI_ISL_2693246)
variant_not_in <- union(variant_stage_1, all_signature_EPI_ISL_6841980) #, a
variant_is_in_signature_EPI_ISL_710528 <- setdiff(all_signature_EPI_ISL_710528, variant_not_in)  %>% unique() # in this not that
variant_is_in_signature_EPI_ISL_710528 %>% length()
saveRDS(variant_is_in_signature_EPI_ISL_710528, file = paste0("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/", "variant_is_in_signature_EPI_ISL_710528", ".rds"), compress = TRUE)
variant_is_in_signature_EPI_ISL_710528
signature_EPI_ISL_710528$SIGNATURE_VARIANT

#union(all_signature_MT007544_1, all_signature_EPI_ISL_710528, all_signature_EPI_ISL_6841980, all_signature_EPI_ISL_2693246)
variant_stage_1 <- union(all_signature_EPI_ISL_2693246, all_signature_MT007544_1) #, all_signature_EPI_ISL_6841980, all_signature_EPI_ISL_2693246)
variant_not_in <- union(variant_stage_1, all_signature_EPI_ISL_710528) #, a
variant_is_in_signature_EPI_ISL_6841980 <- setdiff(all_signature_EPI_ISL_6841980, variant_not_in)  %>% unique() # in this not that
variant_is_in_signature_EPI_ISL_6841980 %>% length()
saveRDS(variant_is_in_signature_EPI_ISL_6841980, file = paste0("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/", "variant_is_in_signature_EPI_ISL_6841980", ".rds"), compress = TRUE)
variant_is_in_signature_EPI_ISL_6841980
signature_EPI_ISL_6841980$SIGNATURE_VARIANT

all_uniq_variants <- union(union(variant_is_in_signature_EPI_ISL_6841980, variant_is_in_signature_EPI_ISL_710528), union(variant_is_in_signature_MT007544_1, variant_is_in_signature_EPI_ISL_2693246)) %>% unique()
data_fwrite <- all_uniq_variants
name_data_frwite <- "all_uniq_variants"
new_filename <- paste0("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/", name_data_frwite, ".csv")
print(new_filename)
data.table::fwrite(data_fwrite, file = new_filename, append = FALSE, quote = "auto", showProgress = getOption("datatable.showProgress", interactive()), buffMB = 1024, nThread = 6)
saveRDS(data_fwrite, file = paste0("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/", "all_uniq_variants", ".rds"), compress = TRUE)


all_vcfs_labelled <- readRDS("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/all_callers_with_features_and_metadata.rds")

calls_variant_is_in_signature_MT007544_1 <- all_vcfs_labelled %>% filter(FILTER == "PASS") %>% filter(THESIS_VARIANT %in% variant_is_in_signature_MT007544_1) %>% arrange(desc(THESIS_COHORT), THESIS_IDENTIFIER, THESIS_BLASTER, THESIS_POS, THESIS_CALLER)
data_fwrite <- calls_variant_is_in_signature_MT007544_1
name_data_frwite <- "calls_variant_is_in_signature_MT007544_1"
new_filename <- paste0("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/", name_data_frwite, ".csv")
print(new_filename)
data.table::fwrite(data_fwrite, file = new_filename, append = FALSE, quote = "auto", showProgress = getOption("datatable.showProgress", interactive()), buffMB = 1024, nThread = 6)

calls_variant_is_in_signature_EPI_ISL_710528 <- all_vcfs_labelled %>% filter(FILTER == "PASS") %>% filter(THESIS_VARIANT %in% variant_is_in_signature_EPI_ISL_710528) %>% arrange(desc(THESIS_COHORT), THESIS_IDENTIFIER, THESIS_BLASTER, THESIS_POS, THESIS_CALLER)
data_fwrite <- calls_variant_is_in_signature_EPI_ISL_710528
name_data_frwite <- "calls_variant_is_in_signature_EPI_ISL_710528"
new_filename <- paste0("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/", name_data_frwite, ".csv")
print(new_filename)
data.table::fwrite(data_fwrite, file = new_filename, append = FALSE, quote = "auto", showProgress = getOption("datatable.showProgress", interactive()), buffMB = 1024, nThread = 6)

# variant_is_in_signature_EPI_ISL_6841980
calls_variant_is_in_signature_EPI_ISL_6841980 <- all_vcfs_labelled %>% filter(FILTER == "PASS") %>% filter(THESIS_VARIANT %in% variant_is_in_signature_EPI_ISL_6841980) %>% arrange(desc(THESIS_COHORT), THESIS_IDENTIFIER, THESIS_BLASTER, THESIS_POS, THESIS_CALLER)
data_fwrite <- calls_variant_is_in_signature_EPI_ISL_6841980
name_data_frwite <- "calls_variant_is_in_signature_EPI_ISL_6841980"
new_filename <- paste0("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/", name_data_frwite, ".csv")
print(new_filename)
data.table::fwrite(data_fwrite, file = new_filename, append = FALSE, quote = "auto", showProgress = getOption("datatable.showProgress", interactive()), buffMB = 1024, nThread = 6)

#variant_is_in_signature_EPI_ISL_2693246
calls_variant_is_in_signature_EPI_ISL_2693246 <- all_vcfs_labelled %>% filter(FILTER == "PASS") %>% filter(THESIS_VARIANT %in% variant_is_in_signature_EPI_ISL_2693246) %>% arrange(desc(THESIS_COHORT), THESIS_IDENTIFIER, THESIS_BLASTER, THESIS_POS, THESIS_CALLER)
data_fwrite <- calls_variant_is_in_signature_EPI_ISL_2693246
name_data_frwite <- "calls_variant_is_in_signature_EPI_ISL_2693246"
new_filename <- paste0("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/", name_data_frwite, ".csv")
print(new_filename)
data.table::fwrite(data_fwrite, file = new_filename, append = FALSE, quote = "auto", showProgress = getOption("datatable.showProgress", interactive()), buffMB = 1024, nThread = 6)

variant_is_in_signature_MT007544_1

save_this_EPI_ISL_6841980 <- calls_variant_is_in_signature_EPI_ISL_6841980 %>% filter(THESIS_SYNTHETIC_PRESENT_C48_OMICRON_EPI_ISL_6841980 == TRUE, THESIS_CALLER == "umivar") %>% select(THESIS_IDENTIFIER,	THESIS_COHORT,	THESIS_BLASTER,	THESIS_CALLER, THESIS_POS, THESIS_VARIANT, THESIS_DEPTH_TOTAL, THESIS_SYNTHETIC_PRESENT_C48_OMICRON_EPI_ISL_6841980_COPIES, THESIS_SYNTHETIC_PRESENT_C48_OMICRON_EPI_ISL_6841980_FRACTION, THESIS_FREQUENCY_NOT_ALTERNATE, THESIS_FREQUENCY_ALTERNATE) %>% arrange(desc(THESIS_SYNTHETIC_PRESENT_C48_OMICRON_EPI_ISL_6841980_FRACTION), THESIS_IDENTIFIER, THESIS_POS, THESIS_VARIANT)

view(save_this_EPI_ISL_6841980)

all_vcfs_labelled %>% filter(THESIS_COHORT == "experimental")

variant_is_in_signature_EPI_ISL_710528

colnames(all_vcfs_labelled)

all_vcfs_labelled %>% select(POS, REF,THESIS_FEATURE_LABEL, THESIS_FEATURE_PRODUCT, THESIS_FEATURE_GENE, THESIS_FEATURE_GENOMIC, THESIS_FEATURE_NCBI) %>% unique() %>% view()

features_ancestral <- all_vcfs_labelled %>% select(POS,THESIS_FEATURE_LABEL, THESIS_FEATURE_PRODUCT, THESIS_FEATURE_GENE, THESIS_FEATURE_GENOMIC, THESIS_FEATURE_NCBI) %>% unique()

features_key <- features_ancestral %>% select(THESIS_FEATURE_LABEL, THESIS_FEATURE_PRODUCT, THESIS_FEATURE_GENE, THESIS_FEATURE_GENOMIC, THESIS_FEATURE_NCBI) %>% unique()

del_signature_EPI_ISL_6841980 <- signature_EPI_ISL_6841980
del_signature_EPI_ISL_2693246 <- signature_EPI_ISL_2693246
del_signature_EPI_ISL_710528 <- signature_EPI_ISL_710528
del_signature_MT007544_1 <- signature_MT007544_1

all_the_variants <- bind_rows(signature_EPI_ISL_6841980, signature_EPI_ISL_2693246, signature_EPI_ISL_710528, signature_MT007544_1)
data_fwrite <- all_the_variants
name_data_frwite <- "all_the_variants"
new_filename <- paste0("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/", name_data_frwite, ".csv")
print(new_filename)
data.table::fwrite(data_fwrite, file = new_filename, append = FALSE, quote = "auto", showProgress = getOption("datatable.showProgress", interactive()), buffMB = 1024, nThread = 6)
saveRDS(data_fwrite, file = paste0("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/", "all_the_variants", ".rds"), compress = TRUE)

names(del_signature_EPI_ISL_6841980) <- c("QUERY_TAG", "SIGNATURE_POS", "REFERENCE", "QUERY", "EPI_ISL_6841980")
names(del_signature_EPI_ISL_2693246) <- c("QUERY_TAG", "SIGNATURE_POS", "REFERENCE", "QUERY", "EPI_ISL_2693246")
names(del_signature_EPI_ISL_710528) <- c("QUERY_TAG", "SIGNATURE_POS", "REFERENCE", "QUERY", "EPI_ISL_710528")
names(del_signature_MT007544_1) <- c("QUERY_TAG", "SIGNATURE_POS", "REFERENCE", "QUERY", "MT007544_1")

del_signature_EPI_ISL_6841980 <- del_signature_EPI_ISL_6841980 %>% select(SIGNATURE_POS, REFERENCE, EPI_ISL_6841980)
del_signature_EPI_ISL_2693246 <- del_signature_EPI_ISL_2693246 %>% select(SIGNATURE_POS, REFERENCE, EPI_ISL_2693246)
del_signature_EPI_ISL_710528 <- del_signature_EPI_ISL_710528 %>% select(SIGNATURE_POS, REFERENCE, EPI_ISL_710528)
del_signature_MT007544_1 <- del_signature_MT007544_1 %>% select(SIGNATURE_POS, REFERENCE, MT007544_1)

tab_1 <- merge(del_signature_EPI_ISL_6841980, del_signature_EPI_ISL_2693246, by = c("SIGNATURE_POS", "REFERENCE"), all = TRUE) %>% arrange(SIGNATURE_POS) %>% tibble()
tab_2 <- merge(tab_1, del_signature_EPI_ISL_710528, by = c("SIGNATURE_POS", "REFERENCE"), all = TRUE) %>% arrange(SIGNATURE_POS) %>% tibble()
tab_3 <- merge(tab_2, del_signature_MT007544_1, by = c("SIGNATURE_POS", "REFERENCE"), all = TRUE) %>% arrange(SIGNATURE_POS) %>% tibble()

signatures_joined <- left_join(tab_3, features_ancestral, by = c('SIGNATURE_POS' = 'POS')) %>% arrange(SIGNATURE_POS) %>% tibble()
signatures_joined %>% view()

signatures_joined$EPI_ISL_6841980 <- as.character(map(strsplit(signatures_joined$EPI_ISL_6841980, split = ":"), 1))
signatures_joined$EPI_ISL_2693246 <- as.character(map(strsplit(signatures_joined$EPI_ISL_2693246, split = ":"), 1))
signatures_joined$EPI_ISL_710528 <- as.character(map(strsplit(signatures_joined$EPI_ISL_710528, split = ":"), 1))
signatures_joined$MT007544_1 <- as.character(map(strsplit(signatures_joined$MT007544_1, split = ":"), 1))
colnames(signatures_joined)
names(signatures_joined)[names(signatures_joined) == 'REFERENCE'] <- 'MN908947_3'
names(signatures_joined)[names(signatures_joined) == 'SIGNATURE_POS'] <- 'POSITION'

signatures_joined %>% select(POSITION, MN908947_3, EPI_ISL_6841980, EPI_ISL_2693246, EPI_ISL_710528, MT007544_1, THESIS_FEATURE_PRODUCT, THESIS_FEATURE_NCBI) %>% view()
appendix_variants <- signatures_joined %>% select(POSITION, MN908947_3, EPI_ISL_6841980, EPI_ISL_2693246, EPI_ISL_710528, MT007544_1, THESIS_FEATURE_PRODUCT, THESIS_FEATURE_NCBI)

names(appendix_variants)[names(appendix_variants) == 'THESIS_FEATURE_PRODUCT'] <- 'Product'
names(appendix_variants)[names(appendix_variants) == 'THESIS_FEATURE_NCBI'] <- 'NCBI'

print(xtable(appendix_variants), include.rownames=FALSE)

data_fwrite <- appendix_variants
name_data_frwite <- "appendix_variants"
new_filename <- paste0("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/", name_data_frwite, ".csv")
print(new_filename)
data.table::fwrite(data_fwrite, file = new_filename, append = FALSE, quote = "auto", showProgress = getOption("datatable.showProgress", interactive()), buffMB = 1024, nThread = 6)
view(signatures_joined)

data_fwrite <- signatures_joined
name_data_frwite <- "signatures_joined"
new_filename <- paste0("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/", name_data_frwite, ".csv")
print(new_filename)
data.table::fwrite(data_fwrite, file = new_filename, append = FALSE, quote = "auto", showProgress = getOption("datatable.showProgress", interactive()), buffMB = 1024, nThread = 6)
view(signatures_joined)

features_key <- features_ancestral %>% select(THESIS_FEATURE_LABEL, THESIS_FEATURE_PRODUCT, THESIS_FEATURE_GENE, THESIS_FEATURE_GENOMIC, THESIS_FEATURE_NCBI) %>% unique()
data_fwrite <- features_key
name_data_frwite <- "features_key"
new_filename <- paste0("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/", name_data_frwite, ".csv")
print(new_filename)
data.table::fwrite(data_fwrite, file = new_filename, append = FALSE, quote = "auto", showProgress = getOption("datatable.showProgress", interactive()), buffMB = 1024, nThread = 6)

gini <- vroom("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/sample_depth_statistics_from_mpileups_gini.csv") %>% unique()
mara <- vroom("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/background_and_consensus_for_marascuillo_procedure_labelled.csv") %>% select(background_cohort, background_identifier) %>% unique()

mara$STAT_IDENTIFIER <- mara$background_identifier

sample_depth_statistics_from_mpileups_gini_cohorted <- left_join(gini, mara)
data_fwrite <- sample_depth_statistics_from_mpileups_gini_cohorted
name_data_frwite <- "sample_depth_statistics_from_mpileups_gini_cohorted"
new_filename <- paste0("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/", name_data_frwite, ".csv")
print(new_filename)
data.table::fwrite(data_fwrite, file = new_filename, append = FALSE, quote = "auto", showProgress = getOption("datatable.showProgress", interactive()), buffMB = 1024, nThread = 6)

