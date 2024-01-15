library(tidyverse)
library(readr)
library(ggplot2)
library(stringr)
library(DescTools)
library(splitstackshape)

MN908947_3_all_snp <-
  read_delim(
    "/home/bgilot/Downloads/SarsCovControls/MN908947.3_all_snp_R.data",
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
  )

gff_MN908947_3 <-
  read_delim(
    "Downloads/SarsCovControls/MN908947.3.gff3",
    "\t",
    escape_double = FALSE,
    col_names = FALSE,
    trim_ws = TRUE,
    skip = 2
  )

names(gff_MN908947_3) <- c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes")
MN908947_3 <- cSplit(gff_MN908947_3, "attributes", ";")


att_ref <- MN908947_3 %>% filter(type=="CDS" | type=="five_prime_UTR") %>% select(start, end, attributes_02)
att_ref$attributes_02 <- str_replace(att_ref$attributes_02, "gbkey=", "Parent=gene-")
att_ref <- cSplit(att_ref, "attributes_02", "=gene-") %>% select(start, end, attributes_02_2)
names(att_ref) <- c("start", "end", "CDS_gene")

SNPs_fingerprint <- MN908947_3_all_snp %>% filter(QUERY != ".") %>% filter(REFERENCE != ".") %>% arrange(REFERENCE_POSITION)
SNPs_fingerprint$gff_start <- as.integer(NA)
SNPs_fingerprint$gff_end <- as.integer(NA)
SNPs_fingerprint$gff_gene <- as.character(NA)
SNPs_fingerprint$TRANSLATED <- as.character(NA)

for(i in 1:nrow(SNPs_fingerprint)) {
  check_position <- SNPs_fingerprint[i,]$REFERENCE_POSITION
  result <- att_ref %>% filter(start <= check_position, end >= check_position) 
  how_many <- nrow(result)
  debut <- (att_ref %>% filter(start <= check_position, end >= check_position))$start
  fin <- (att_ref %>% filter(start <= check_position, end >= check_position))$end
  quoi <- (att_ref %>% filter(start <= check_position, end >= check_position))$CDS_gene

  translated <- as.character("TRANSLATED")
  if (how_many != 1 ){
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
  # do stuff with row
}

SNPs_EPI_ISL_2693246 <- SNPs_fingerprint %>% filter(TRANSLATED == "TRANSLATED", REFERENCE != ".", QUERY != ".", QUERY_TAG == "EPI_ISL_2693246") %>% arrange(REFERENCE_POSITION)
SNPs_EPI_ISL_710528 <- SNPs_fingerprint %>% filter(TRANSLATED == "TRANSLATED", REFERENCE != ".", QUERY != ".", QUERY_TAG == "EPI_ISL_710528") %>% arrange(REFERENCE_POSITION)
SNPs_EPI_ISL_6841980 <- SNPs_fingerprint %>% filter(TRANSLATED == "TRANSLATED", REFERENCE != ".", QUERY != ".", QUERY_TAG == "EPI_ISL_6841980") %>% arrange(REFERENCE_POSITION)
SNPs_MT007544_1 <- SNPs_fingerprint %>% filter(TRANSLATED == "TRANSLATED", REFERENCE != ".", QUERY != ".", QUERY_TAG == "MT007544.1") %>% arrange(REFERENCE_POSITION)

SNPs_delta <- SNPs_EPI_ISL_2693246
SNPs_alpha <- SNPs_EPI_ISL_710528
SNPs_omicron <- SNPs_EPI_ISL_6841980
SNPs_Australia <- SNPs_MT007544_1

genes <- c("S")
variants <- SNPs_fingerprint$QUERY_TAG %>% sort() %>% unique()

for_plot <- SNPs_fingerprint %>% filter(TRANSLATED == "TRANSLATED", REFERENCE != ".", QUERY != ".") %>% filter(gff_gene %in% genes) %>% arrange(REFERENCE_POSITION, QUERY_TAG)
for_plot$REFERENCE_POSITION %>% sort() %>% unique()


# SNPs_EPI_ISL_710528 <- MN908947_3_all_snp %>% filter(QUERY != ".", QUERY_TAG == "EPI_ISL_710528") %>% arrange(REFERENCE_POSITION)
# SNPs_EPI_ISL_6841980 <- MN908947_3_all_snp %>% filter(QUERY != ".", QUERY_TAG == "EPI_ISL_6841980") %>% arrange(REFERENCE_POSITION)
# SNPs_MT007544_1 <- MN908947_3_all_snp %>% filter(QUERY != ".", QUERY_TAG == "MT007544.1") %>% arrange(REFERENCE_POSITION)
# 
# SNPs_delta <- MN908947_3_all_snp %>% filter(QUERY != ".", QUERY_TAG == "EPI_ISL_2693246") %>% arrange(REFERENCE_POSITION)
# SNPs_alpha <- MN908947_3_all_snp %>% filter(QUERY != ".", QUERY_TAG == "EPI_ISL_710528") %>% arrange(REFERENCE_POSITION)
# SNPs_omicron <- MN908947_3_all_snp %>% filter(QUERY != ".", QUERY_TAG == "EPI_ISL_6841980") %>% arrange(REFERENCE_POSITION)
# SNPs_MT007544_1 <- MN908947_3_all_snp %>% filter(QUERY != ".", QUERY_TAG == "MT007544.1") %>% arrange(REFERENCE_POSITION)

