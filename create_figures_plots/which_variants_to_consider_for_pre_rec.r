# library(tidyverse)
# library(vroom)
# library(xtable)
# library(kableExtra)
# library(ggvenn)
# library(ggVennDiagram)
# library(VennDiagram)
# library(gplots)
# library(RColorBrewer)

# run this first   /mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/create_figures_plots/wrapper_thesis.r
# run this first   /mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/create_figures_plots/depth_plots.r

#minimum_depth <- 500

all_uniq_variants <- readRDS(file = paste0("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/", "all_uniq_variants", ".rds"))
all_the_variants <- readRDS(file = paste0("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/", "all_the_variants", ".rds"))


# SAMPLE METADATA 1 -------------------------------------------------------


sample_depth_statistics_from_mpileups_gini <- vroom("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/counts_from_mpileups/sample_depth_statistics_from_mpileups_gini.csv")
data_depth <- sample_depth_statistics_from_mpileups_gini %>% select(STAT_IDENTIFIER, STAT_SAMPLE_METHOD, STAT_DEPTH_MEDIAN, STAT_DEPTH_MEAN, STAT_DEPTH_MAX) %>% unique()
names(data_depth)  <- c("SAMPLE_NAME", "METHOD", "DEPTH_MED", "DEPTH_MEA", "DEPTH_MAX")

background_and_consensus_for_marascuillo <- vroom("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/background_and_consensus_for_marascuillo_procedure_labelled.csv")
data_cts <- background_and_consensus_for_marascuillo %>% select(background_identifier, background_blaster, background_cohort, THESIS_SAMPLE_CT_S) %>% unique()
names(data_cts)  <- c("SAMPLE_NAME", "METHOD", "COHORT", "CT_S")

label_data <- merge(data_depth, data_cts)
# data_label_select <-
#   label_data %>% filter(METHOD == methodi)  %>% filter(COHORT == cohorti) %>% filter(SAMPLE_NAME %in% (
#     counts_depth_plots %>% filter(SAMPLE_METHOD == methodi) %>% filter(SAMPLE_GROUP == cohorti)
#   )$SAMPLE_NAME)
label_data$Identifier  <- as.character(map(strsplit(label_data$SAMPLE_NAME, split = "_"), 1))

label_data %>% select(Identifier, CT_S) %>% unique()
# SAMPLE METADATA 2 -------------------------------------------------------


depths_blasted <- vroom("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/sample_depth_statistics_from_mpileups_gini_cohorted.csv") %>% 
  filter(STAT_SAMPLE_METHOD == "blasted") %>% 
  select(STAT_IDENTIFIER, STAT_DEPTH_MIN, STAT_DEPTH_MEDIAN, STAT_DEPTH_MEAN, STAT_DEPTH_MAX) %>% 
  unique() %>%
  arrange(STAT_IDENTIFIER)
colnames(depths_blasted) <- c("Identifier", "Min_blind", "Median_blind", "Mean_blind", "Max_blind")
glimpse(depths_blasted)

depths_unblasted <- vroom("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/sample_depth_statistics_from_mpileups_gini_cohorted.csv") %>% 
  filter(STAT_SAMPLE_METHOD == "unblasted") %>% 
  select(STAT_IDENTIFIER , STAT_DEPTH_MIN, STAT_DEPTH_MEDIAN, STAT_DEPTH_MEAN, STAT_DEPTH_MAX) %>% 
  unique() %>%
  arrange(STAT_IDENTIFIER)
colnames(depths_unblasted) <- c("Identifier", "Min_aware", "Median_aware", "Mean_aware", "Max_aware")
glimpse(depths_unblasted)

depths_by_method <- left_join(depths_blasted, depths_unblasted)

depths_by_method$sort_mean <- rowMeans(cbind(depths_by_method$Min_blind, depths_by_method$Median_blind, depths_by_method$Mean_blind, depths_by_method$Max_blind, depths_by_method$Min_aware, depths_by_method$Median_aware, depths_by_method$Mean_aware, depths_by_method$Max_aware)) %>% as.numeric() %>% round(1)
depths_by_method$Identifier  <- as.character(map(strsplit(depths_by_method$Identifier, split = "_"), 1))

metadata_by_sample <- left_join((label_data %>% select(Identifier, CT_S) %>% unique()), depths_by_method) %>% unique()

metadata_by_sample$Mean_blind <- metadata_by_sample$Mean_blind %>% round(1)
metadata_by_sample$Mean_aware <- metadata_by_sample$Mean_blind %>% round(1)

# other data --------------------------------------------------------------



SNPs_nucmer_MT007544_1_australia <- readRDS("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_computed/SNPs_MT007544_1_australia.rds")
SNPs_nucmer_EPI_ISL_omicron <- readRDS("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_computed/SNPs_EPI_ISL_omicron.rds")
SNPs_nucmer_EPI_ISL_2693246_delta <- readRDS("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_computed/SNPs_EPI_ISL_2693246_delta.rds")
SNPs_nucmer_EPI_ISL_710528_alpha <- readRDS("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_computed/SNPs_EPI_ISL_710528_alpha.rds")
glimpse(SNPs_nucmer_EPI_ISL_710528_alpha)

only_in_signature_EPI_ISL_710528 <- readRDS(file = paste0("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/", "variant_is_in_signature_EPI_ISL_710528", ".rds"))
only_in_signature_MT007544_1 <- readRDS(file = paste0("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/", "variant_is_in_signature_MT007544_1", ".rds"))
only_in_signature_EPI_ISL_6841980 <- readRDS(file = paste0("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/", "variant_is_in_signature_EPI_ISL_6841980", ".rds"))
only_in_signature_EPI_ISL_2693246 <- readRDS(file = paste0("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/", "variant_is_in_signature_EPI_ISL_2693246", ".rds"))

# all
all_vcfs_labelled <- readRDS("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/all_callers_with_features_and_metadata.rds") %>% filter(THESIS_COHORT == "control") %>% filter(FILTER == "PASS") # %>% filter(THESIS_VARIANT %in% all_uniq_variants) 
all_vcfs_labelled$IDENTIFIER <- as.character(map(strsplit(all_vcfs_labelled$THESIS_IDENTIFIER, split = "_"), 1))
all_vcfs_labelled$Identifier <- as.character(map(strsplit(all_vcfs_labelled$THESIS_IDENTIFIER, split = "_"), 1))



#  %>% filter(SIGNATURE_POS >= 21563) %>% filter(SIGNATURE_POS <= 25384)  
C1_Australia <- (all_the_variants %>% filter(QUERY_TAG == "MT007544.1"))$SIGNATURE_VARIANT  %>% unique() %>% sort()# Australia_MT007544_1
C2_Ancestral <- as.character(NA)
C14_Alpha <- (all_the_variants %>% filter(QUERY_TAG == "EPI_ISL_710528"))$SIGNATURE_VARIANT %>% unique() %>% sort()# Alpha_EPI_ISL_710528
C29_Delta <- (all_the_variants %>% filter(QUERY_TAG == "EPI_ISL_2693246"))$SIGNATURE_VARIANT %>% unique() %>% sort()# Delta_EPI_ISL_2693246
C48_Omicron  <- (all_the_variants %>% filter(QUERY_TAG == "EPI_ISL_6841980"))$SIGNATURE_VARIANT %>% unique() %>% sort() #EPI_ISL_6841980
C_ALL <- c(C1_Australia, C14_Alpha, C29_Delta, C48_Omicron) %>% unique() %>% sort()
#list(C1_Australia, C2_Ancestral, C14_Alpha, C29_Delta, C48_Omicron)

C1_Australia_negative <- setdiff(C_ALL, C1_Australia) %>% unique() %>% sort()
C2_Ancestral_negative <- setdiff(C_ALL, C2_Ancestral) %>% unique() %>% sort()
C14_Alpha_negative <- setdiff(C_ALL, C14_Alpha) %>% unique() %>% sort()
C29_Delta_negative <- setdiff(C_ALL, C29_Delta) %>% unique() %>% sort()
C48_Omicron_negative <- setdiff(C_ALL, C48_Omicron) %>% unique() %>% sort()
C_ALL_negative <- c(C1_Australia_negative, C2_Ancestral_negative, C14_Alpha_negative, C29_Delta_negative, C48_Omicron_negative) %>% unique() %>% sort()
#list(C1_Australia, C2_Ancestral, C14_Alpha, C29_Delta, C48_Omicron)
#list(C1_Australia_negative, C2_Ancestral_negative, C14_Alpha_negative, C29_Delta_negative, C48_Omicron_negative)

all_vcf_C1 <- all_vcfs_labelled %>% filter(Identifier %in% (readRDS("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/synthetic_samples_C1.rds")$Identifier %>% unique()))
all_vcf_C2 <- all_vcfs_labelled %>% filter(Identifier %in% (readRDS("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/synthetic_samples_C2.rds")$Identifier %>% unique()))
all_vcf_C14 <- all_vcfs_labelled %>% filter(Identifier %in% (readRDS("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/synthetic_samples_C14.rds")$Identifier %>% unique()))
all_vcf_C29 <- all_vcfs_labelled %>% filter(Identifier %in% (readRDS("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/synthetic_samples_C29.rds")$Identifier %>% unique()))
all_vcf_C48 <- all_vcfs_labelled %>% filter(Identifier %in% (readRDS("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/synthetic_samples_C48.rds")$Identifier %>% unique())) # %>% filter(THESIS_CALLER == "umivar")

control_variant_SNP <- list(C1_Australia, C2_Ancestral, C14_Alpha, C29_Delta, C48_Omicron)
control_variant_SNP_negation <- list(C1_Australia_negative, C2_Ancestral_negative, C14_Alpha_negative, C29_Delta_negative, C48_Omicron_negative)
control_variant_VCF <- list(all_vcf_C1, all_vcf_C2, all_vcf_C14, all_vcf_C29, all_vcf_C48)
control_variant_text_label <- list("C1_Australia", "C2_Ancestral", "C14_Alpha", "C29_Delta", "C48_Omicron")
control_variant_column_copies <- list("THESIS_SYNTHETIC_PRESENT_C1_AUSTRALIA_MT007544_1_COPIES", "THESIS_SYNTHETIC_PRESENT_C2_ANCESTRAL_MN908947_3_COPIES", "THESIS_SYNTHETIC_PRESENT_C14_ALPHA_EPI_ISL_710528_COPIES", "THESIS_SYNTHETIC_PRESENT_C29_DELTA_EPI_ISL_2693246_COPIES", "THESIS_SYNTHETIC_PRESENT_C48_OMICRON_EPI_ISL_6841980_COPIES")
control_variant_column_fraction <- list("THESIS_SYNTHETIC_PRESENT_C1_AUSTRALIA_MT007544_1_FRACTION", "THESIS_SYNTHETIC_PRESENT_C2_ANCESTRAL_MN908947_3_FRACTION", "THESIS_SYNTHETIC_PRESENT_C14_ALPHA_EPI_ISL_710528_FRACTION", "THESIS_SYNTHETIC_PRESENT_C29_DELTA_EPI_ISL_2693246_FRACTION", "THESIS_SYNTHETIC_PRESENT_C48_OMICRON_EPI_ISL_6841980_FRACTION")




# Begin _ Create _ Table --------------------------------------------------

per_sample_results <- tibble(
  minimum_depth = integer(),
  minimum_depth_alternate = integer(),
  reference_length = integer(),
  twist_control = integer(),
  ident_sample = character(),
  control_by_caller = character(),
  
  loop_control_variant_text_label = character(),
  loop_control_variant_column_copies = character(),
  loop_control_variant_column_fraction = character(),
  
  working_vcf_rows = integer(),
  #expected_snp = character(),
  expected_snp_chr = character(),
  expected_snp_len = integer(),
  #observed_snp = character(),
  observed_snp_chr = character(),
  observed_snp_len = integer(),
  #negation_snp = character(),
  negation_snp_chr = character(),
  negation_snp_len = integer(),
  #corrected_snp = character(),
  corrected_snp_chr = character(),
  corrected_snp_len = integer(),
  #other_control_snp = character(),
  other_control_snp_chr = character(),
  other_control_snp_len = integer(),
  #true_positive_snp = character(),
  true_positive_snp_chr = character(),
  true_positive_snp_len = integer(),
  #false_positive_snp = character(),
  false_positive_snp_chr = character(),
  false_positive_snp_len = integer(),
  #false_negative_snp = character(),
  false_negative_snp_chr = character(),
  false_negative_snp_len = integer(),
  corrected_ref_len = integer(),
  human_cts = numeric(),
  blind_minimum_depth = numeric(),
  blind_median_depth = numeric(),
  blind_mean_depth = numeric(),
  blind_max_depth = numeric(),
  aware_minimum_depth = numeric(),
  aware_median_depth = numeric(),
  aware_mean_depth = numeric(),
  aware_max_depth = numeric(),
  sort_mean = numeric(),
  # metadata_synthetic = tibble(),
  tp_mean = numeric(),
  tp_median = numeric(),
  tp_sd = numeric(),
  tp_summary = character(),
  fp_mean = numeric(),
  fp_median = numeric(),
  fp_sd = numeric(),
  fp_summary = character()
)


minimum_depth <- 100
minimum_depth_alternate <- 10
reference_length <- 29903

twist_control_by_C_number <-
  seq(1, (control_variant_VCF %>% length()), 1)
for (twist_control in twist_control_by_C_number) {
  #control_variant_VCF[[twist_control]] %>% glimpse() %>% print()
  #control_variant_VCF[[twist_control]] %>% tibble() %>% print()
  sample_identifiers_by_control <-
    control_variant_VCF[[twist_control]]$Identifier %>% unique() %>% sort()
  for (ident_sample in sample_identifiers_by_control){
    
    control_caller_list <- (control_variant_VCF[[twist_control]] %>% filter(Identifier == ident_sample))$THESIS_CALLER %>% unique() %>% sort()
    for (control_by_caller in control_caller_list) {

      loop_control_variant_text_label <- control_variant_text_label[[twist_control]]
      loop_control_variant_column_copies <- control_variant_column_copies[[twist_control]]
      loop_control_variant_column_fraction <- control_variant_column_fraction[[twist_control]]
      
       
      working_vcf <- control_variant_VCF[[twist_control]] %>% filter(THESIS_DEPTH_ALTERNATE >= minimum_depth_alternate) %>% filter(THESIS_DEPTH_TOTAL >= minimum_depth) %>% filter(FILTER == "PASS") %>% filter(Identifier == ident_sample) %>% filter(THESIS_CALLER == control_by_caller) %>% unique()
      working_vcf_rows <- working_vcf %>% nrow()
      # glimpse(working_vcf)
      
      #working_vcf_variant_list <- (control_variant_VCF[[twist_control]] %>% filter(Identifier == ident_sample) %>% filter(THESIS_CALLER == control_by_caller) %>% unique())$THESIS_VARIANT %>% unique() %>% sort()
      
      
      
      
      expected_snp <- control_variant_SNP[[twist_control]] %>% unique()
      expected_snp_chr <- paste(as.character(expected_snp), sep="' '", collapse=", ") 
      (expected_snp_chr %>% str_split(", "))[[1]] == expected_snp
      expected_snp_len <- expected_snp %>% length()
      
      observed_snp <- working_vcf$THESIS_VARIANT %>% unique()
      observed_snp_chr <- paste(as.character(observed_snp), sep="' '", collapse=", ") 
      (observed_snp_chr %>% str_split(", "))[[1]] == observed_snp
      observed_snp_len <- observed_snp %>% length()
      
      negation_snp <- control_variant_SNP_negation[[twist_control]]
      negation_snp_chr <- paste(as.character(negation_snp), sep="' '", collapse=", ") 
      (negation_snp_chr %>% str_split(", "))[[1]] == negation_snp
      negation_snp_len <- negation_snp %>% length()
      
      corrected_snp <- setdiff(observed_snp, control_variant_SNP_negation[[twist_control]])
      corrected_snp_chr <- paste(as.character(corrected_snp), sep="' '", collapse=", ") 
      (corrected_snp_chr %>% str_split(", "))[[1]] == corrected_snp
      corrected_snp_len <- corrected_snp %>% length()
      
      other_control_snp <- setdiff(observed_snp, corrected_snp)
      other_control_snp_chr <- paste(as.character(other_control_snp), sep="' '", collapse=", ") 
      (other_control_snp_chr %>% str_split(", "))[[1]] == other_control_snp
      other_control_snp_len <- other_control_snp %>% length()
      
      true_positive_snp <- intersect(expected_snp, corrected_snp)
      true_positive_snp_chr <- paste(as.character(true_positive_snp), sep="' '", collapse=", ") 
      (true_positive_snp_chr %>% str_split(", "))[[1]] == true_positive_snp
      true_positive_snp_len <- true_positive_snp %>% length()
      
      false_positive_snp <- setdiff(corrected_snp, expected_snp)
      false_positive_snp_chr <- paste(as.character(false_positive_snp), sep="' '", collapse=", ") 
      (false_positive_snp_chr %>% str_split(", "))[[1]] == false_positive_snp
      false_positive_snp_len <- false_positive_snp %>% length()
      
      false_negative_snp <- setdiff(expected_snp, true_positive_snp)
      false_negative_snp_chr <- paste(as.character(false_negative_snp), sep="' '", collapse=", ") 
      (false_negative_snp_chr %>% str_split(", "))[[1]] == false_negative_snp
      false_negative_snp_len <- false_negative_snp %>% length()
      
      corrected_ref_len <- reference_length - negation_snp_len
      
      human_cts <- (metadata_by_sample %>% filter(Identifier == ident_sample))$CT_S
      blind_minimum_depth <- (metadata_by_sample %>% filter(Identifier == ident_sample))$Min_blind
      blind_median_depth <- (metadata_by_sample %>% filter(Identifier == ident_sample))$Median_blind
      blind_mean_depth <- (metadata_by_sample %>% filter(Identifier == ident_sample))$Mean_blind
      blind_max_depth <- (metadata_by_sample %>% filter(Identifier == ident_sample))$Max_blind
      
      aware_minimum_depth <- (metadata_by_sample %>% filter(Identifier == ident_sample))$Min_aware
      aware_median_depth <- (metadata_by_sample %>% filter(Identifier == ident_sample))$Median_aware
      aware_mean_depth <- (metadata_by_sample %>% filter(Identifier == ident_sample))$Mean_aware
      aware_max_depth <- (metadata_by_sample %>% filter(Identifier == ident_sample))$Max_aware
      sort_mean <- (metadata_by_sample %>% filter(Identifier == ident_sample))$sort_mean
      
      metadata_synthetic <- working_vcf %>% select(Identifier, starts_with("THESIS_SYNTHETIC_PRESENT_")) %>% unique()
      #(working_vcf %>% filter(THESIS_VARIANT %in% true_positive_snp))$THESIS_VARIANT
      
      # alelle frequency TP
      tp_mean <- (working_vcf %>% filter(THESIS_VARIANT %in% true_positive_snp))$THESIS_FREQUENCY_ALTERNATE %>% mean()
      tp_median <- (working_vcf %>% filter(THESIS_VARIANT %in% true_positive_snp))$THESIS_FREQUENCY_ALTERNATE %>% median()
      tp_sd <- (working_vcf %>% filter(THESIS_VARIANT %in% true_positive_snp))$THESIS_FREQUENCY_ALTERNATE %>% sd()
      tp_summary <- (working_vcf %>% filter(THESIS_VARIANT %in% true_positive_snp))$THESIS_FREQUENCY_ALTERNATE %>% summary() %>% round(6) %>% as.character() %>% paste(collapse=", ")
      
      # alelle frequency FP
      fp_mean <- (working_vcf %>% filter(THESIS_VARIANT %in% false_positive_snp))$THESIS_FREQUENCY_ALTERNATE %>% mean()
      fp_median <- (working_vcf %>% filter(THESIS_VARIANT %in% false_positive_snp))$THESIS_FREQUENCY_ALTERNATE %>% median()
      fp_sd <- (working_vcf %>% filter(THESIS_VARIANT %in% false_positive_snp))$THESIS_FREQUENCY_ALTERNATE %>% sd()
      fp_summary <- (working_vcf %>% filter(THESIS_VARIANT %in% false_positive_snp))$THESIS_FREQUENCY_ALTERNATE %>% summary() %>% round(6) %>% as.character() %>% paste(collapse=", ")
      
      
      # print("Expected:")
      # print(expected_snp)
      # print("Observed:")
      # print(observed_snp)
      # print("Corrected:")
      # print(corrected_snp)
      print(paste(twist_control, control_variant_text_label[[twist_control]], ident_sample, control_by_caller, paste("Rows:", working_vcf_rows),  sep = " : "))
      print(metadata_by_sample %>% filter(Identifier == ident_sample))
      
      print(paste(human_cts, blind_minimum_depth, blind_median_depth, blind_mean_depth, blind_max_depth, aware_minimum_depth, aware_median_depth, aware_mean_depth, aware_max_depth, sort_mean,  sep = " : "))
      
      print(paste("expected_snp_len", expected_snp_len, "observed_snp_len", observed_snp_len,  "negation_snp_len",  negation_snp_len, "corrected_snp_len", corrected_snp_len, "other_control_snp_len", other_control_snp_len, "true_positive_snp_len", true_positive_snp_len, "false_positive_snp_len", false_positive_snp_len, "false_negative_snp_len", false_negative_snp_len, sep = " : "))
    print(metadata_synthetic)
    print(minimum_depth)
    print(minimum_depth_alternate)
    

# add row inside loop -----------------------------------------------------

    
    new_row_list <- list(
      minimum_depth,
      minimum_depth_alternate,
      reference_length,
      twist_control,
      ident_sample,
      control_by_caller,
      
      loop_control_variant_text_label,
      loop_control_variant_column_copies,
      loop_control_variant_column_fraction,
      
      working_vcf_rows,
      #expected_snp,
      expected_snp_chr,
      expected_snp_len,
      #observed_snp,
      observed_snp_chr,
      observed_snp_len,
      #negation_snp,
      negation_snp_chr,
      negation_snp_len,
      #corrected_snp,
      corrected_snp_chr,
      corrected_snp_len,
      #other_control_snp,
      other_control_snp_chr,
      other_control_snp_len,
      #true_positive_snp,
      true_positive_snp_chr,
      true_positive_snp_len,
      #false_positive_snp,
      false_positive_snp_chr,
      false_positive_snp_len,
      #false_negative_snp,
      false_negative_snp_chr,
      false_negative_snp_len,
      corrected_ref_len,
      human_cts,
      blind_minimum_depth,
      blind_median_depth,
      blind_mean_depth,
      blind_max_depth,
      aware_minimum_depth,
      aware_median_depth,
      aware_mean_depth,
      aware_max_depth,
      sort_mean,
      # metadata_synthetic,
      tp_mean,
      tp_median,
      tp_sd,
      tp_summary,
      fp_mean,
      fp_median,
      fp_sd,
      fp_summary
    )
    
    per_sample_results[nrow(per_sample_results) + 1,] <- new_row_list    
    
    
      }
  }
}
#metadata_synthetic <- working_vcf %>% select(Identifier, starts_with("THESIS_SYNTHETIC_PRESENT_")) %>% unique()

# to Join  
sample_info <- all_vcfs_labelled %>% select(Identifier, THESIS_SYNTHETIC_COPIES_TOTAL, starts_with("THESIS_SYNTHETIC_PRESENT_")) %>% unique()
per_sample_results$Identifier <- per_sample_results$ident_sample

master_performance <- left_join(per_sample_results, sample_info)
# (reference - correction) - (FP + FN + TP)
master_performance$true_negative_snp_len <- (master_performance$reference_length - master_performance$negation_snp_len) - (master_performance$false_positive_snp_len + master_performance$false_negative_snp_len + master_performance$true_positive_snp_len)
master_performance$precision <- master_performance$true_positive_snp_len / (master_performance$true_positive_snp_len + master_performance$false_positive_snp_len)
master_performance$recall <- master_performance$true_positive_snp_len / (master_performance$true_positive_snp_len + master_performance$false_negative_snp_len)

master_performance %>% glimpse()


# add relevant data section -----------------------------------------------


master_performance$RELEVANT_COPIES <- as.integer(0)
master_performance$RELEVANT_FRACTION <- as.numeric(0)
for (relevant_row in seq(1,nrow(master_performance),1)) {
  
  colname_copies <- master_performance[relevant_row,]$loop_control_variant_column_copies
  colname_fraction <- master_performance[relevant_row,]$loop_control_variant_column_fraction
  
  master_performance[relevant_row,]$RELEVANT_COPIES <- master_performance[relevant_row,][[colname_copies]]
  master_performance[relevant_row,]$RELEVANT_FRACTION <- master_performance[relevant_row,][[colname_fraction]] %>% round(5)
}
master_performance %>% glimpse()

# create table concordant results -----------------------------------------


working_00 <- master_performance
master_performance_consensus <- master_performance

master_performance_consensus$concordant_label <- as.character(NA)

master_performance_consensus$corrected_snp_chr_ivar <- as.character(NA)
master_performance_consensus$corrected_snp_chr_ivar_len <- as.numeric(NA)

master_performance_consensus$corrected_snp_chr_lofreq <- as.character(NA)
master_performance_consensus$corrected_snp_chr_lofreq_len <- as.numeric(NA)

master_performance_consensus$corrected_snp_chr_umivar <- as.character(NA)
master_performance_consensus$corrected_snp_chr_umivar_len <- as.numeric(NA)

master_performance_consensus$corrected_snp_chr_varscan <- as.character(NA)
master_performance_consensus$corrected_snp_chr_varscan_len <- as.numeric(NA)

master_performance_consensus$corrected_snp_chr_test <- as.character(NA)
master_performance_consensus$corrected_snp_chr_test_len <- as.numeric(NA)

working_00 <- master_performance
{
  list_target <-
    working_00$loop_control_variant_text_label %>% unique() %>% sort()
for (loop_target in list_target) {
  working_01 <-
    working_00 %>% filter(loop_control_variant_text_label == loop_target)
  list_copies <- working_01$RELEVANT_COPIES %>% unique() %>% sort()
  for (loop_copies in list_copies) {
    working_02 <- working_01 %>% filter(RELEVANT_COPIES == loop_copies)
    list_fractions <-
      working_02$RELEVANT_FRACTION %>% unique() %>% sort()
    for (loop_fractions in list_fractions) {
      working_03 <-
        working_02 %>% filter(RELEVANT_FRACTION == loop_fractions)
      list_ident <- working_03$Identifier %>% unique() %>% sort()
      for (loop_ident in list_ident) {
        working_04 <- working_03 %>% filter(Identifier == loop_ident)
        list_caller <- working_04$control_by_caller %>% unique() %>% sort()
        
        for (loop_caller in list_caller) {
          working_05 <- working_04 %>% filter(control_by_caller == loop_caller)
          
          
          # working_05$expected_snp_chr
          # working_05$expected_snp_len
          # working_05$corrected_snp_chr
          # working_05$corrected_snp_len
          
          print(paste(
            loop_target,
            loop_copies,
            loop_fractions,
            loop_ident,
            loop_caller,
            sep = " : "
          ))
          
          concordant_label <- paste(
            loop_target,
            loop_copies,
            loop_fractions,
            loop_ident,
            loop_caller,
            sep = " : "
          )
          
          master_performance_consensus[master_performance_consensus$loop_control_variant_text_label == loop_target & master_performance_consensus$RELEVANT_COPIES == loop_copies & master_performance_consensus$RELEVANT_FRACTION == loop_fractions & master_performance_consensus$Identifier == loop_ident & master_performance_consensus$control_by_caller == loop_caller, "concordant_label"] <- concordant_label
          #print(master_performance_consensus[master_performance_consensus$loop_control_variant_text_label == loop_target & master_performance_consensus$RELEVANT_COPIES == loop_copies & master_performance_consensus$RELEVANT_FRACTION == loop_fractions & master_performance_consensus$Identifier == loop_ident,]) 
          
          if (loop_caller == "ivar") {
            print("ivar: yes")
            master_performance_consensus[master_performance_consensus$loop_control_variant_text_label == loop_target & master_performance_consensus$RELEVANT_COPIES == loop_copies & master_performance_consensus$RELEVANT_FRACTION == loop_fractions & master_performance_consensus$Identifier == loop_ident, "corrected_snp_chr_ivar"] <- working_05$corrected_snp_chr
            print(master_performance_consensus[master_performance_consensus$loop_control_variant_text_label == loop_target & master_performance_consensus$RELEVANT_COPIES == loop_copies & master_performance_consensus$RELEVANT_FRACTION == loop_fractions & master_performance_consensus$Identifier == loop_ident, "corrected_snp_chr_ivar"])
             } else {
            print("ivar: no")
          }
          
          if (loop_caller == "lofreq") {
            print("lofreq: yes")
            master_performance_consensus[master_performance_consensus$loop_control_variant_text_label == loop_target & master_performance_consensus$RELEVANT_COPIES == loop_copies & master_performance_consensus$RELEVANT_FRACTION == loop_fractions & master_performance_consensus$Identifier == loop_ident, "corrected_snp_chr_lofreq"] <- working_05$corrected_snp_chr
            print(master_performance_consensus[master_performance_consensus$loop_control_variant_text_label == loop_target & master_performance_consensus$RELEVANT_COPIES == loop_copies & master_performance_consensus$RELEVANT_FRACTION == loop_fractions & master_performance_consensus$Identifier == loop_ident, "corrected_snp_chr_lofreq"])
                      } else {
            print("lofreq: no")
          }
          
          if (loop_caller == "umivar") {
            print("umivar: yes")
            master_performance_consensus[master_performance_consensus$loop_control_variant_text_label == loop_target & master_performance_consensus$RELEVANT_COPIES == loop_copies & master_performance_consensus$RELEVANT_FRACTION == loop_fractions & master_performance_consensus$Identifier == loop_ident, "corrected_snp_chr_umivar"] <- working_05$corrected_snp_chr
            print(master_performance_consensus[master_performance_consensus$loop_control_variant_text_label == loop_target & master_performance_consensus$RELEVANT_COPIES == loop_copies & master_performance_consensus$RELEVANT_FRACTION == loop_fractions & master_performance_consensus$Identifier == loop_ident, "corrected_snp_chr_umivar"])
            #master_performance_consensus[master_performance_consensus$loop_control_variant_text_label == loop_target & master_performance_consensus$RELEVANT_COPIES == loop_copies, "corrected_snp_chr_test"] <- working_05$corrected_snp_chr
          } else {
            print("umivar: no")
          }
          
          if (loop_caller == "varscan") {
            print("varscan: yes")
            master_performance_consensus[master_performance_consensus$loop_control_variant_text_label == loop_target & master_performance_consensus$RELEVANT_COPIES == loop_copies & master_performance_consensus$RELEVANT_FRACTION == loop_fractions & master_performance_consensus$Identifier == loop_ident, "corrected_snp_chr_varscan"] <- working_05$corrected_snp_chr
            print(master_performance_consensus[master_performance_consensus$loop_control_variant_text_label == loop_target & master_performance_consensus$RELEVANT_COPIES == loop_copies & master_performance_consensus$RELEVANT_FRACTION == loop_fractions & master_performance_consensus$Identifier == loop_ident, "corrected_snp_chr_varscan"])
            #master_performance_consensus[master_performance_consensus$loop_control_variant_text_label == loop_target & master_performance_consensus$RELEVANT_COPIES == loop_copies, "corrected_snp_chr_test"] <- working_05$corrected_snp_chr
          } else {
            print("varscan: no")
          }
          
        }
        
        
        
      }
    }
  }
}
}


master_performance_consensus %>% glimpse()
data_fwrite <- master_performance_consensus
name_data_frwite <- "master_performance_consensus"
new_filename <- paste0(folder_data_unified, name_data_frwite, ".csv")
print(new_filename)
fwrite(data_fwrite, file = new_filename, append = FALSE, quote = "auto", showProgress = getOption("datatable.showProgress", interactive()), buffMB = 1024, nThread = 6)
new_filename <- paste0(folder_data_unified, name_data_frwite, ".rds")
print(new_filename)
saveRDS(data_fwrite, file = new_filename, compress = TRUE)
readRDS(new_filename)

master_performance_consensus <- readRDS(new_filename)
master_performance_consensus %>% glimpse()


# Create / save  venn diagrams ---------------------------------------------

master_performance_consensus_safety <- master_performance_consensus

if(!exists(".Random.seed")) set.seed(NULL)
set.seed(sample(seq(1,100000,1),1))
set.seed(sample(seq(1,10000,1),1))
set.seed(sample(seq(1,1000,1),1))
set.seed(sample(seq(1,1000000,1),1))
run_rows <- sample(seq(1, nrow(master_performance_consensus),1),1)
run_rows

for (run_rows in seq(1, nrow(master_performance_consensus),1)) {
  #master_performance_consensus$corrected_snp_chr_ivar_len <- 
  
  plot_identifier <- master_performance_consensus[run_rows,]$Identifier
  plot_target <- master_performance_consensus[run_rows,]$loop_control_variant_text_label
  plot_copies <- master_performance_consensus[run_rows,]$RELEVANT_COPIES
  plot_fraction <- master_performance_consensus[run_rows,]$RELEVANT_FRACTION
  umi_median_depth <- master_performance_consensus[run_rows,]$aware_median_depth
  umi_blind_median_depth <- master_performance_consensus[run_rows,]$blind_median_depth
  af_mean <- master_performance_consensus[run_rows,]$tp_mean
  plot_caller <- master_performance_consensus[run_rows,]$control_by_caller
  
  plot_tp <- master_performance_consensus[run_rows,]$true_positive_snp_len
  plot_tn <- master_performance_consensus[run_rows,]$true_negative_snp_len
  plot_fp <- master_performance_consensus[run_rows,]$false_positive_snp_len
  plot_fn <- master_performance_consensus[run_rows,]$false_negative_snp_len
  plot_precision <- master_performance_consensus[run_rows,]$precision
  plot_recall <- master_performance_consensus[run_rows,]$recall
  
  ex_len <- master_performance_consensus[run_rows,]$expected_snp_len
  ex <- master_performance_consensus[run_rows,]$expected_snp_chr %>% str_split(", ") %>% unlist()
  
  if (is.na(ex)) {
    ex_len <- 0
  }
  
  iv <- master_performance_consensus[run_rows,]$corrected_snp_chr_ivar %>% str_split(", ") %>% unlist()
  lo <- master_performance_consensus[run_rows,]$corrected_snp_chr_lofreq %>% str_split(", ") %>% unlist()
  um <- master_performance_consensus[run_rows,]$corrected_snp_chr_umivar %>% str_split(", ") %>% unlist()
  va <- master_performance_consensus[run_rows,]$corrected_snp_chr_varscan %>% str_split(", ") %>% unlist()
    iv_len <- master_performance_consensus[run_rows,]$corrected_snp_chr_ivar %>% str_split(", ") %>% unlist() %>% length()
    lo_len <- master_performance_consensus[run_rows,]$corrected_snp_chr_lofreq %>% str_split(", ") %>% unlist() %>% length()
    um_len <- master_performance_consensus[run_rows,]$corrected_snp_chr_umivar %>% str_split(", ") %>% unlist() %>% length()
    va_len <- master_performance_consensus[run_rows,]$corrected_snp_chr_varscan %>% str_split(", ") %>% unlist() %>% length()
    
    master_performance_consensus[run_rows,]$corrected_snp_chr_ivar_len <- iv_len
    master_performance_consensus[run_rows,]$corrected_snp_chr_lofreq_len <- lo_len
    master_performance_consensus[run_rows,]$corrected_snp_chr_umivar_len <- um_len
    master_performance_consensus[run_rows,]$corrected_snp_chr_varscan_len <- va_len

    iv_lo <- intersect(iv,lo)
    iv_um <- intersect(iv,um)
    iv_va <- intersect(iv,va)
    lo_um <- intersect(lo,um)
    lo_va <- intersect(lo,va)
    um_va <- intersect(um,va)
    
    iv_ex <- intersect(iv,ex)
    lo_ex <- intersect(lo,ex)
    um_ex <- intersect(um,ex)
    va_ex <- intersect(va,ex)
    
    iv_ex_len <- iv_ex %>% unique() %>% length()
    lo_ex_len <- lo_ex %>% unique() %>% length()
    um_ex_len <- um_ex %>% unique() %>% length()
    va_ex_len <- va_ex %>% unique() %>% length()
   
    
    
    
    concordant_calls <- intersect(iv_lo,um_va)
    concordant_calls_len <- concordant_calls %>% length()
    
   union_calls <- union(( union(iv,lo) %>% unique()),    (union(um,va) %>% unique())) %>% unique()
   union_expected <- intersect(ex,union_calls)
   union_expected_len <- union_expected %>% unique() %>% length()
   
   union_expected_len - concordant_calls_len
   
    iv_lo_len <- iv_lo %>% unique() %>% length()
    iv_um_len <- iv_um %>% unique() %>% length()
    iv_va_len <- iv_va %>% unique() %>% length()
    lo_um_len <- lo_um %>% unique() %>% length()
    lo_va_len <- lo_va %>% unique() %>% length()
    um_va_len <- um_va %>% unique() %>% length()
    
    expected_concordant <- intersect(concordant_calls, ex)
    expected_concordant_len <- expected_concordant %>% length()
    
    union_expected_len - expected_concordant_len
    
    print(paste(iv_len,lo_len,um_len,va_len,sep = ":"))
    print(paste(iv_lo_len,iv_um_len,iv_va_len,lo_um_len, lo_va_len,um_va_len, sep = ":"))
    print(paste(iv_ex_len,lo_ex_len,um_ex_len,va_ex_len, sep = ":"))
    #print(concordant_calls)
    print((paste("Expect",ex_len, "Concordant Calls", concordant_calls_len,  "Intersect Concordant_expected",expected_concordant_len, "union_expected_len", union_expected_len, sep = ":")))
    
    venn_twist <- list(
      ivar = iv,
      lofreq = lo,
      varscan = va,
      umivar = um
    )
    
    
    if (plot_target == "C2_Ancestral") {
      ex_len <- 0
    plot_fn <- 0
    af_mean <- 0
    }
    
    plot_OBJECT <-
      ggVennDiagram(venn_twist, label_alpha = 0) + scale_x_continuous(expand = expansion(mult = .2)) + scale_fill_distiller(palette = "RdBu") + scale_color_brewer(palette = "Set3")  + #labs(title = NULL,subtitle = NULL,caption = "Mapped against MN908947_3")+labs(
      labs(
        title = NULL,
        subtitle = NULL,
        caption = paste0(
          plot_identifier,
          # " :: ",
          # plot_target,
          " :: ",
          plot_caller,
          "\n",
          # "Copies = ",
          # comma(plot_copies, 0),
          # " :: ",
          # "Fraction = ",
          # comma(plot_fraction, 3),
          #           "\n",
          paste0(comma(plot_copies, 0), " copies of ", plot_target, " @ fraction of ", comma(plot_fraction, 3)),
          "\n",
          "UMI Aware Median Depth = ",
          comma(umi_median_depth, digits = 0),
          "\n",
          "UMI Blind Median Depth = ",
          comma(umi_blind_median_depth, digits = 0),
          "\n",

          #"FN : TP : FP : TN = ",
          paste(paste0("precision = ", comma(plot_precision, digits = 4)), paste0("recall = ", comma(plot_recall, digits = 4)) ,sep = " : "),
          "\n",
          paste(paste0("FN = ",plot_fn), paste0("TP = ",plot_tp), paste0("FP = ",comma(plot_fp, digits = 0)), paste0("TN = ",comma(plot_tn, digits = 0)), sep = " : "),
          "\n",
          
          #"precision : recall = ",
          "Mean Alternate Allele Frequency = ",
          comma(af_mean, 4),
          "\n",
          "Concordant = ",
          concordant_calls_len,
          " SNPs :: Concordant True Positives = ",
          expected_concordant_len,
          " SNPs",
          " :: Expect = ",
          ex_len,
          " SNPs"
        )
      )
    #plot_caller
    #file <- paste0(folder_plots,"consensus_", plot_target, "_target_", plot_copies, "_copies_", plot_fraction, "_fraction_", run_rows, "_row_", plot_identifier, ".png")
    file <- paste0(folder_plots,"concordance_copies_", plot_copies,  "_fraction_", plot_fraction,  "_", plot_target, "_", plot_identifier, "_datarow_", run_rows, "_", plot_caller, ".png")
    ggsave(file, device = "png", dpi = 1800,  type = "cairo", width = 210, height = 148.5, units = "mm") 
    print(file)
}



# add plot data to master -------------------------------------------------
#######################  
master_performance_consensus %>% glimpse()
master_performance_consensus_super <- master_performance_consensus

master_performance_consensus_super %>% glimpse()

master_performance_consensus_super$plot_identifier <- as.character(NA)
master_performance_consensus_super$plot_target <- as.character(NA)
master_performance_consensus_super$plot_copies <- as.numeric(NA)
master_performance_consensus_super$plot_fraction <- as.numeric(NA)
master_performance_consensus_super$plot_umi_median_depth <- as.numeric(NA)
master_performance_consensus_super$plot_umi_blind_median_depth <- as.numeric(NA)
master_performance_consensus_super$plot_af_mean <- as.numeric(NA)
master_performance_consensus_super$plot_caller <- as.character(NA)
master_performance_consensus_super$plot_tp <- as.numeric(NA)
master_performance_consensus_super$plot_tn <- as.numeric(NA)
master_performance_consensus_super$plot_fp <- as.numeric(NA)
master_performance_consensus_super$plot_fn <- as.numeric(NA)
master_performance_consensus_super$plot_precision <- as.numeric(NA)
master_performance_consensus_super$plot_recall <- as.numeric(NA)
master_performance_consensus_super$plot_ex_len <- as.numeric(NA)
master_performance_consensus_super$plot_ex <- as.character(NA)
master_performance_consensus_super$plot_iv <- as.character(NA)
master_performance_consensus_super$plot_lo <- as.character(NA)
master_performance_consensus_super$plot_um <- as.character(NA)
master_performance_consensus_super$plot_va <- as.character(NA)
master_performance_consensus_super$plot_iv_len <- as.numeric(NA)
master_performance_consensus_super$plot_lo_len <- as.numeric(NA)
master_performance_consensus_super$plot_um_len <- as.numeric(NA)
master_performance_consensus_super$plot_va_len <- as.numeric(NA)
master_performance_consensus_super$plot_iv_lo <- as.character(NA)
master_performance_consensus_super$plot_iv_um <- as.character(NA)
master_performance_consensus_super$plot_iv_va <- as.character(NA)
master_performance_consensus_super$plot_lo_um <- as.character(NA)
master_performance_consensus_super$plot_lo_va <- as.character(NA)
master_performance_consensus_super$plot_um_va <- as.character(NA)
master_performance_consensus_super$plot_iv_ex <- as.character(NA)
master_performance_consensus_super$plot_lo_ex <- as.character(NA)
master_performance_consensus_super$plot_um_ex <- as.character(NA)
master_performance_consensus_super$plot_va_ex <- as.character(NA)
master_performance_consensus_super$plot_iv_ex_len <- as.numeric(NA)
master_performance_consensus_super$plot_lo_ex_len <- as.numeric(NA)
master_performance_consensus_super$plot_um_ex_len <- as.numeric(NA)
master_performance_consensus_super$plot_va_ex_len <- as.numeric(NA)
master_performance_consensus_super$plot_concordant_calls <- as.character(NA)
master_performance_consensus_super$plot_concordant_calls_len <- as.numeric(NA)
master_performance_consensus_super$plot_union_calls <- as.character(NA)
master_performance_consensus_super$plot_union_expected <- as.character(NA)
master_performance_consensus_super$plot_union_expected_len <- as.numeric(NA)
master_performance_consensus_super$plot_iv_lo_len <- as.numeric(NA)
master_performance_consensus_super$plot_iv_um_len <- as.numeric(NA)
master_performance_consensus_super$plot_iv_va_len <- as.numeric(NA)
master_performance_consensus_super$plot_lo_um_len <- as.numeric(NA)
master_performance_consensus_super$plot_lo_va_len <- as.numeric(NA)
master_performance_consensus_super$plot_um_va_len <- as.numeric(NA)
master_performance_consensus_super$plot_expected_concordant <- as.character(NA)
master_performance_consensus_super$plot_expected_concordant_len <- as.numeric(NA)
master_performance_consensus_super$plot_ex_len <- as.numeric(NA)
master_performance_consensus_super$plot_fn <- as.numeric(NA)
master_performance_consensus_super$plot_af_mean <- as.numeric(NA)
master_performance_consensus_super$plot_file <- as.character(NA)
master_performance_consensus_super$plot_legend  <- as.character(NA)

##############################$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# iv <- master_performance_consensus[run_rows,]$corrected_snp_chr_ivar %>% str_split(", ") %>% unlist()
# paste(as.character(expected_snp), sep="' '", collapse=", ") 
# paste(as.character(expected_snp), sep="' '", collapse=", ") %>% str_split(", ")  %>% unlist()
##############################$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

if(!exists(".Random.seed")) set.seed(NULL)
set.seed(sample(seq(1,100000,1),1))
set.seed(sample(seq(1,10000,1),1))
set.seed(sample(seq(1,1000,1),1))
set.seed(sample(seq(1,1000000,1),1))
run_rows <- sample(seq(1, nrow(master_performance_consensus_super),1),1)
run_rows

## add additional columns

master_performance_consensus_super %>% glimpse()

# master_performance_consensus_super[run_rows,]$corrected_snp_chr_ivar_len <- iv_len
# master_performance_consensus_super[run_rows,]$corrected_snp_chr_lofreq_len <- lo_len
# master_performance_consensus_super[run_rows,]$corrected_snp_chr_umivar_len <- um_len
# master_performance_consensus_super[run_rows,]$corrected_snp_chr_varscan_len <- va_len

for (run_rows in seq(1, nrow(master_performance_consensus_super),1)) {
  plot_identifier <- master_performance_consensus_super[run_rows,]$Identifier
  plot_target <- master_performance_consensus_super[run_rows,]$loop_control_variant_text_label
  plot_copies <- master_performance_consensus_super[run_rows,]$RELEVANT_COPIES
  plot_fraction <- master_performance_consensus_super[run_rows,]$RELEVANT_FRACTION
  umi_median_depth <- master_performance_consensus_super[run_rows,]$aware_median_depth
  umi_blind_median_depth <- master_performance_consensus_super[run_rows,]$blind_median_depth
  af_mean <- master_performance_consensus_super[run_rows,]$tp_mean
  plot_caller <- master_performance_consensus_super[run_rows,]$control_by_caller
  
  plot_tp <- master_performance_consensus_super[run_rows,]$true_positive_snp_len
  plot_tn <- master_performance_consensus_super[run_rows,]$true_negative_snp_len
  plot_fp <- master_performance_consensus_super[run_rows,]$false_positive_snp_len
  plot_fn <- master_performance_consensus_super[run_rows,]$false_negative_snp_len
  plot_precision <- master_performance_consensus_super[run_rows,]$precision
  plot_recall <- master_performance_consensus_super[run_rows,]$recall
  
  ex_len <- master_performance_consensus_super[run_rows,]$expected_snp_len
  ex <- master_performance_consensus_super[run_rows,]$expected_snp_chr %>% str_split(", ") %>% unlist()
  
  if (is.na(ex)) {
    ex_len <- 0
  }
  
  iv <- master_performance_consensus_super[run_rows,]$corrected_snp_chr_ivar %>% str_split(", ") %>% unlist()
  lo <- master_performance_consensus_super[run_rows,]$corrected_snp_chr_lofreq %>% str_split(", ") %>% unlist()
  um <- master_performance_consensus_super[run_rows,]$corrected_snp_chr_umivar %>% str_split(", ") %>% unlist()
  va <- master_performance_consensus_super[run_rows,]$corrected_snp_chr_varscan %>% str_split(", ") %>% unlist()
  
  iv_len <- master_performance_consensus_super[run_rows,]$corrected_snp_chr_ivar %>% str_split(", ") %>% unlist() %>% length()
  lo_len <- master_performance_consensus_super[run_rows,]$corrected_snp_chr_lofreq %>% str_split(", ") %>% unlist() %>% length()
  um_len <- master_performance_consensus_super[run_rows,]$corrected_snp_chr_umivar %>% str_split(", ") %>% unlist() %>% length()
  va_len <- master_performance_consensus_super[run_rows,]$corrected_snp_chr_varscan %>% str_split(", ") %>% unlist() %>% length()
  
  iv_lo <- intersect(iv,lo)
  iv_um <- intersect(iv,um)
  iv_va <- intersect(iv,va)
  lo_um <- intersect(lo,um)
  lo_va <- intersect(lo,va)
  um_va <- intersect(um,va)
  
  iv_ex <- intersect(iv,ex)
  lo_ex <- intersect(lo,ex)
  um_ex <- intersect(um,ex)
  va_ex <- intersect(va,ex)
  
  iv_ex_len <- iv_ex %>% unique() %>% length()
  lo_ex_len <- lo_ex %>% unique() %>% length()
  um_ex_len <- um_ex %>% unique() %>% length()
  va_ex_len <- va_ex %>% unique() %>% length()
  
  concordant_calls <- intersect(iv_lo,um_va)
  concordant_calls_len <- concordant_calls %>% length()
  
  union_calls <- union(( union(iv,lo) %>% unique()),    (union(um,va) %>% unique())) %>% unique()
  union_expected <- intersect(ex,union_calls)
  union_expected_len <- union_expected %>% unique() %>% length()
  
  union_expected_len - concordant_calls_len
  
  iv_lo_len <- iv_lo %>% unique() %>% length()
  iv_um_len <- iv_um %>% unique() %>% length()
  iv_va_len <- iv_va %>% unique() %>% length()
  lo_um_len <- lo_um %>% unique() %>% length()
  lo_va_len <- lo_va %>% unique() %>% length()
  um_va_len <- um_va %>% unique() %>% length()
  
  expected_concordant <- intersect(concordant_calls, ex)
  expected_concordant_len <- expected_concordant %>% length()
  
  union_expected_len - expected_concordant_len
  
  print(paste(iv_len,lo_len,um_len,va_len,sep = ":"))
  print(paste(iv_lo_len,iv_um_len,iv_va_len,lo_um_len, lo_va_len,um_va_len, sep = ":"))
  print(paste(iv_ex_len,lo_ex_len,um_ex_len,va_ex_len, sep = ":"))
  print((paste("Expect",ex_len, "Concordant Calls", concordant_calls_len,  "Intersect Concordant_expected",expected_concordant_len, "union_expected_len", union_expected_len, sep = ":")))
  
  venn_twist <- list(
    ivar = iv,
    lofreq = lo,
    varscan = va,
    umivar = um
  )
  
  if (plot_target == "C2_Ancestral") {
    ex_len <- 0
    plot_fn <- 0
    af_mean <- 0
  }
  
  
  plot_legend <- paste0(
    plot_identifier,
    " :: ",
    plot_caller,
    "\n",
    paste0(comma(plot_copies, 0), " copies of ", plot_target, " @ fraction of ", comma(plot_fraction, 3)),
    "\n",
    "UMI Aware Median Depth = ",
    comma(umi_median_depth, digits = 0),
    "\n",
    "UMI Blind Median Depth = ",
    comma(umi_blind_median_depth, digits = 0),
    "\n",
    paste(paste0("precision = ", comma(plot_precision, digits = 4)), paste0("recall = ", comma(plot_recall, digits = 4)) ,sep = " : "),
    "\n",
    paste(paste0("FN = ",plot_fn), paste0("TP = ",plot_tp), paste0("FP = ",comma(plot_fp, digits = 0)), paste0("TN = ",comma(plot_tn, digits = 0)), sep = " : "),
    "\n",
    "Mean Alternate Allele Frequency = ",
    comma(af_mean, 4),
    "\n",
    "Concordant = ",
    concordant_calls_len,
    " SNPs :: Concordant True Positives = ",
    expected_concordant_len,
    " SNPs",
    " :: Expect = ",
    ex_len,
    " SNPs"
  )
  
  file <- paste0(folder_plots,"concordance_copies_", plot_copies,  "_fraction_", plot_fraction,  "_", plot_target, "_", plot_identifier, "_datarow_", run_rows, "_", plot_caller, ".png")
  
  # iv <- master_performance_consensus[run_rows,]$corrected_snp_chr_ivar %>% str_split(", ") %>% unlist()
  # paste(as.character(expected_snp), sep="' '", collapse=", ") 
  # paste(as.character(expected_snp), sep="' '", collapse=", ") %>% str_split(", ")  %>% unlist()
  
  master_performance_consensus_super[run_rows,]$plot_identifier <- plot_identifier
  master_performance_consensus_super[run_rows,]$plot_target <- plot_target
  master_performance_consensus_super[run_rows,]$plot_copies <- plot_copies
  master_performance_consensus_super[run_rows,]$plot_fraction <- plot_fraction
  master_performance_consensus_super[run_rows,]$plot_umi_median_depth <- umi_median_depth
  master_performance_consensus_super[run_rows,]$plot_umi_blind_median_depth <- umi_blind_median_depth
  master_performance_consensus_super[run_rows,]$plot_af_mean <- af_mean
  master_performance_consensus_super[run_rows,]$plot_caller <- plot_caller
  master_performance_consensus_super[run_rows,]$plot_tp <- plot_tp
  master_performance_consensus_super[run_rows,]$plot_tn <- plot_tn
  master_performance_consensus_super[run_rows,]$plot_fp <- plot_fp
  master_performance_consensus_super[run_rows,]$plot_fn <- plot_fn
  master_performance_consensus_super[run_rows,]$plot_precision <- plot_precision
  master_performance_consensus_super[run_rows,]$plot_recall <- plot_recall
  master_performance_consensus_super[run_rows,]$plot_ex_len <- ex_len
  master_performance_consensus_super[run_rows,]$plot_ex <- ex %>% as.character() %>% paste(sep="' '", collapse=", ")
  master_performance_consensus_super[run_rows,]$plot_iv <- iv %>% as.character() %>% paste(sep="' '", collapse=", ")
  master_performance_consensus_super[run_rows,]$plot_lo <- lo %>% as.character() %>% paste(sep="' '", collapse=", ")
  master_performance_consensus_super[run_rows,]$plot_um <- um %>% as.character() %>% paste(sep="' '", collapse=", ")
  master_performance_consensus_super[run_rows,]$plot_va <- va %>% as.character() %>% paste(sep="' '", collapse=", ")
  master_performance_consensus_super[run_rows,]$plot_iv_len <- iv_len
  master_performance_consensus_super[run_rows,]$plot_lo_len <- lo_len
  master_performance_consensus_super[run_rows,]$plot_um_len <- um_len
  master_performance_consensus_super[run_rows,]$plot_va_len <- va_len
  master_performance_consensus_super[run_rows,]$plot_iv_lo <- iv_lo %>% as.character() %>% paste(sep="' '", collapse=", ")
  master_performance_consensus_super[run_rows,]$plot_iv_um <- iv_um %>% as.character() %>% paste(sep="' '", collapse=", ")
  master_performance_consensus_super[run_rows,]$plot_iv_va <- iv_va %>% as.character() %>% paste(sep="' '", collapse=", ")
  master_performance_consensus_super[run_rows,]$plot_lo_um <- lo_um %>% as.character() %>% paste(sep="' '", collapse=", ")
  master_performance_consensus_super[run_rows,]$plot_lo_va <- lo_va %>% as.character() %>% paste(sep="' '", collapse=", ")
  master_performance_consensus_super[run_rows,]$plot_um_va <- um_va %>% as.character() %>% paste(sep="' '", collapse=", ")
  master_performance_consensus_super[run_rows,]$plot_iv_ex <- iv_ex %>% as.character() %>% paste(sep="' '", collapse=", ")
  master_performance_consensus_super[run_rows,]$plot_lo_ex <- lo_ex %>% as.character() %>% paste(sep="' '", collapse=", ")
  master_performance_consensus_super[run_rows,]$plot_um_ex <- um_ex %>% as.character() %>% paste(sep="' '", collapse=", ")
  master_performance_consensus_super[run_rows,]$plot_va_ex <- va_ex %>% as.character() %>% paste(sep="' '", collapse=", ")
  master_performance_consensus_super[run_rows,]$plot_iv_ex_len <- iv_ex_len
  master_performance_consensus_super[run_rows,]$plot_lo_ex_len <- lo_ex_len
  master_performance_consensus_super[run_rows,]$plot_um_ex_len <- um_ex_len
  master_performance_consensus_super[run_rows,]$plot_va_ex_len <- va_ex_len
  master_performance_consensus_super[run_rows,]$plot_concordant_calls <- concordant_calls %>% as.character() %>% paste(sep="' '", collapse=", ")
  master_performance_consensus_super[run_rows,]$plot_concordant_calls_len <- concordant_calls_len
  master_performance_consensus_super[run_rows,]$plot_union_calls <- union_calls %>% as.character() %>% paste(sep="' '", collapse=", ")
  master_performance_consensus_super[run_rows,]$plot_union_expected <- union_expected %>% as.character() %>% paste(sep="' '", collapse=", ")
  master_performance_consensus_super[run_rows,]$plot_union_expected_len <- union_expected_len
  master_performance_consensus_super[run_rows,]$plot_iv_lo_len <- iv_lo_len
  master_performance_consensus_super[run_rows,]$plot_iv_um_len <- iv_um_len
  master_performance_consensus_super[run_rows,]$plot_iv_va_len <- iv_va_len
  master_performance_consensus_super[run_rows,]$plot_lo_um_len <- lo_um_len
  master_performance_consensus_super[run_rows,]$plot_lo_va_len <- lo_va_len
  master_performance_consensus_super[run_rows,]$plot_um_va_len <- um_va_len
  master_performance_consensus_super[run_rows,]$plot_expected_concordant <- expected_concordant %>% as.character() %>% paste(sep="' '", collapse=", ")
  master_performance_consensus_super[run_rows,]$plot_expected_concordant_len <- expected_concordant_len
  master_performance_consensus_super[run_rows,]$plot_ex_len <- ex_len
  master_performance_consensus_super[run_rows,]$plot_fn <- plot_fn
  master_performance_consensus_super[run_rows,]$plot_af_mean <- af_mean
  master_performance_consensus_super[run_rows,]$plot_file <- file
  master_performance_consensus_super[run_rows,]$plot_legend <- plot_legend
  
  print(paste("Done:",file))
  
}



master_performance_consensus_super %>% glimpse()
data_fwrite <- master_performance_consensus_super
name_data_frwite <- "master_performance_consensus_super"
new_filename <- paste0(folder_data_unified, name_data_frwite, ".csv")
print(new_filename)
fwrite(data_fwrite, file = new_filename, append = FALSE, quote = "auto", showProgress = getOption("datatable.showProgress", interactive()), buffMB = 1024, nThread = 6)
new_filename <- paste0(folder_data_unified, name_data_frwite, ".rds")
print(new_filename)
saveRDS(data_fwrite, file = new_filename, compress = TRUE)
readRDS(new_filename)
master_performance_consensus_super <- readRDS(new_filename)
master_performance_consensus_super %>% glimpse()


# compute precision/recall concordant -------------------------------------


master_performance_consensus_super_concordant_confusion_matrix <- master_performance_consensus_super
master_performance_consensus_super_concordant_confusion_matrix %>% glimpse()

master_performance_consensus_super_concordant_confusion_matrix$concordant_sample_negated_snps <-  master_performance_consensus_super_concordant_confusion_matrix$negation_snp_len 
# EXPECTED SNPs
master_performance_consensus_super_concordant_confusion_matrix$concordant_sample_expected_snps <-  master_performance_consensus_super_concordant_confusion_matrix$plot_ex_len  
# CORRECTED REFERENCE LENGTH
master_performance_consensus_super_concordant_confusion_matrix$concordant_corrected_ref_len <- master_performance_consensus_super_concordant_confusion_matrix$corrected_ref_len
# POSITIVES / CONCORDANT
master_performance_consensus_super_concordant_confusion_matrix$concordant_positive <- master_performance_consensus_super_concordant_confusion_matrix$plot_concordant_calls_len # POSITIVES / CONCORDANT
# NEGATIVES / DISCORDANT
master_performance_consensus_super_concordant_confusion_matrix$concordant_negative <- master_performance_consensus_super_concordant_confusion_matrix$concordant_corrected_ref_len - master_performance_consensus_super_concordant_confusion_matrix$concordant_positive # NEGATIVES / DISCORDANT

# TEST
master_performance_consensus_super_concordant_confusion_matrix$concordant_corrected_ref_len == (master_performance_consensus_super_concordant_confusion_matrix$concordant_positive + master_performance_consensus_super_concordant_confusion_matrix$concordant_negative)
# TRUE POSITIVE
master_performance_consensus_super_concordant_confusion_matrix$concordant_tp <- master_performance_consensus_super_concordant_confusion_matrix$plot_expected_concordant_len
# FALSE POSITIVE
master_performance_consensus_super_concordant_confusion_matrix$concordant_fp <- master_performance_consensus_super_concordant_confusion_matrix$concordant_positive - master_performance_consensus_super_concordant_confusion_matrix$concordant_tp
# FALSE NEGATIVE
master_performance_consensus_super_concordant_confusion_matrix$concordant_fn <- master_performance_consensus_super_concordant_confusion_matrix$plot_ex_len - master_performance_consensus_super_concordant_confusion_matrix$concordant_tp
# TRUE NEGATIVE
master_performance_consensus_super_concordant_confusion_matrix$concordant_tn <- master_performance_consensus_super_concordant_confusion_matrix$concordant_negative - master_performance_consensus_super_concordant_confusion_matrix$concordant_fn

master_performance_consensus_super_concordant_confusion_matrix$concordant_sample_negated_snps  + master_performance_consensus_super_concordant_confusion_matrix$concordant_tp + master_performance_consensus_super_concordant_confusion_matrix$concordant_fp + master_performance_consensus_super_concordant_confusion_matrix$concordant_fn + master_performance_consensus_super_concordant_confusion_matrix$concordant_tn

master_performance_consensus_super_concordant_confusion_matrix$concordant_precision <- master_performance_consensus_super_concordant_confusion_matrix$concordant_tp / (master_performance_consensus_super_concordant_confusion_matrix$concordant_tp + master_performance_consensus_super_concordant_confusion_matrix$concordant_fp)
master_performance_consensus_super_concordant_confusion_matrix$concordant_recall <- master_performance_consensus_super_concordant_confusion_matrix$concordant_tp / (master_performance_consensus_super_concordant_confusion_matrix$concordant_tp + master_performance_consensus_super_concordant_confusion_matrix$concordant_fn)


master_performance_consensus_super_concordant_confusion_matrix %>% glimpse()
data_fwrite <- master_performance_consensus_super_concordant_confusion_matrix
name_data_frwite <- "master_performance_consensus_super_concordant_confusion_matrix"
new_filename <- paste0(folder_data_unified, name_data_frwite, ".csv")
print(new_filename)
fwrite(data_fwrite, file = new_filename, append = FALSE, quote = "auto", showProgress = getOption("datatable.showProgress", interactive()), buffMB = 1024, nThread = 6)
new_filename <- paste0(folder_data_unified, name_data_frwite, ".rds")
print(new_filename)
saveRDS(data_fwrite, file = new_filename, compress = TRUE)
readRDS(new_filename)
master_performance_consensus_super_concordant_confusion_matrix <- readRDS(new_filename)
master_performance_consensus_super_concordant_confusion_matrix %>% glimpse()


# 
# # precision recall --------------------------------------------------------
# master_performance_consensus_super <- readRDS("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/master_performance_consensus_super.rds")
# 
# master_performance_consensus_super %>% 
#   filter(RELEVANT_COPIES  >= 1000, RELEVANT_COPIES  <= 1000000) %>% 
#   filter(RELEVANT_FRACTION  >= 0.0001, RELEVANT_FRACTION  <= 0.5) %>% 
#   # filter(control_by_caller == "umivar") %>% 
#   # filter(false_negative_snp_len == 0) %>%
#   #filter(loop_control_variant_text_label == "C1_Australia") %>% 
#   ggplot(aes(x = recall,
#              y = precision)) +
#   geom_point(alpha = 15 / 20, size = 0.5) + 
#   facet_grid(control_by_caller ~ RELEVANT_FRACTION)  +
#   geom_smooth(
#     method = "auto",
#     se = TRUE,
#     fullrange = FALSE,
#     level = 0.95,
#     linewidth = 0.5,
#     alpha = 0.5,
#     color = "yellow"
#   )   + 
#   scale_y_continuous(expand = c(0, 0)) +
#   scale_x_continuous(expand = c(0, 0)) 
# 
# 
# master_performance %>% ggplot(aes(x = recall,
#                                   y = precision)) +
#   geom_point(alpha = 10 / 20, size = 0.5) + 
#   facet_wrap(control_by_caller ~ THESIS_SYNTHETIC_PRESENT_C48_OMICRON_EPI_ISL_6841980_COPIES)  +
#   geom_smooth(
#     method = "auto",
#     se = TRUE,
#     fullrange = FALSE,
#     level = 0.95,
#     linewidth = 0.5,
#     alpha = 0.5,
#     color = "yellow"
#   ) 


# paste(as.character(expected_snp), sep="' '", collapse=", ") 
# paste(as.character(expected_snp), sep="' '", collapse=", ") %>% str_split(", ")

# 
# new_row_list <- list(
#   minimum_depth,
#   minimum_depth_alternate,
#   reference_length,
#   twist_control,
#   ident_sample,
#   control_by_caller,
#   working_vcf_rows,
#   #expected_snp,
#   expected_snp_chr,
#   expected_snp_len,
#   #observed_snp,
#   observed_snp_chr,
#   observed_snp_len,
#   #negation_snp,
#   negation_snp_chr,
#   negation_snp_len,
#   #corrected_snp,
#   corrected_snp_chr,
#   corrected_snp_len,
#   #other_control_snp,
#   other_control_snp_chr,
#   other_control_snp_len,
#   #true_positive_snp,
#   true_positive_snp_chr,
#   true_positive_snp_len,
#   #false_positive_snp,
#   false_positive_snp_chr,
#   false_positive_snp_len,
#   #false_negative_snp,
#   false_negative_snp_chr,
#   false_negative_snp_len,
#   corrected_ref_len,
#   human_cts,
#   blind_minimum_depth,
#   blind_median_depth,
#   blind_mean_depth,
#   blind_max_depth,
#   aware_minimum_depth,
#   aware_median_depth,
#   aware_mean_depth,
#   aware_max_depth,
#   sort_mean,
#   # metadata_synthetic,
#   tp_mean,
#   tp_median,
#   tp_sd,
#   tp_summary,
#   fp_mean,
#   fp_median,
#   fp_sd,
#   fp_summary
# )
# 
# per_sample_results[nrow(per_sample_results) + 1,] <- new_row_list

# 
# minimum_depth,
# minimum_depth_alternate,
# twist_control,
# ident_sample,
# control_by_caller,
# working_vcf_rows,
# expected_snp,
# expected_snp_len,
# observed_snp,
# observed_snp_len,
# negation_snp,
# negation_snp_len,
# corrected_snp,
# corrected_snp_len,
# other_control_snp,
# other_control_snp_len,
# true_positive_snp,
# true_positive_snp_len,
# false_positive_snp,
# false_positive_snp_len,
# false_negative_snp,
# false_negative_snp_len,
# human_cts,
# blind_minimum_depth,
# blind_median_depth,
# blind_mean_depth,
# blind_max_depth,
# aware_minimum_depth,
# aware_median_depth,
# aware_mean_depth,
# aware_max_depth,
# sort_mean,
# metadata_synthetic,
# tp_mean,
# tp_median,
# tp_sd,
# tp_summary,
# fp_mean,
# fp_median,
# fp_sd,
# fp_summary






# glimpse(all_vcf_C1)
# all_vcfs_labelled$Identifier %>% unique() %>% sort()
# all_vcf_C2$THESIS_VARIANT %>% unique() %>% sort()
# all_vcf_C2 %>% filter(!(THESIS_VARIANT %in% C2_Ancestral_negative))
# 
# list()
# 
# 
# spike_C14 <- intersect(only_in_signature_EPI_ISL_710528, C14_Alpha)
# spike_C1 <-intersect(only_in_signature_MT007544_1, C1_Australia)
# spike_C48 <-intersect(only_in_signature_EPI_ISL_6841980, C48_Omicron)
# spike_C29 <-intersect(only_in_signature_EPI_ISL_2693246, C29_Delta)
# 
# venn_twist <- list(
#   C1_Australia = (all_the_variants %>% filter(QUERY_TAG == "MT007544.1"))$SIGNATURE_VARIANT, # Australia_MT007544_1
#   C14_Alpha = (all_the_variants %>% filter(QUERY_TAG == "EPI_ISL_710528"))$SIGNATURE_VARIANT, # Alpha_EPI_ISL_710528
#   C29_Delta = (all_the_variants %>% filter(QUERY_TAG == "EPI_ISL_2693246"))$SIGNATURE_VARIANT, # Delta_EPI_ISL_2693246
#   C48_Omicron = (all_the_variants %>% filter(QUERY_TAG == "EPI_ISL_6841980"))$SIGNATURE_VARIANT # Omicron_EPI_ISL_6841980
# )
# ggVennDiagram(venn_twist, label_alpha = 0) + scale_x_continuous(expand = expansion(mult = .2)) + scale_fill_distiller(palette = "RdBu") + scale_color_brewer(palette = "Set3") + labs(title = NULL,
#                                                                                                                                                                                       subtitle = NULL,
#                                                                                                                                                                                       caption = "Pangenome Mapped against MN908947_3")
# 
# file <- paste0("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/GILOT_THESIS/figures/","twist_controls_pangenome.png")
# ggsave(file, device = "png", dpi = "retina", type = "cairo")
# 
# # all_the_variants %>% filter(SIGNATURE_POS >= 21563) %>% filter(SIGNATURE_POS <= 25384) 
# 
# venn_twist <- list(
#   C1_Australia = (all_the_variants %>% filter(SIGNATURE_POS >= 21563) %>% filter(SIGNATURE_POS <= 25384)  %>% filter(QUERY_TAG == "MT007544.1"))$SIGNATURE_VARIANT, # Australia_MT007544_1
#   C14_Alpha = (all_the_variants %>% filter(SIGNATURE_POS >= 21563) %>% filter(SIGNATURE_POS <= 25384)  %>% filter(QUERY_TAG == "EPI_ISL_710528"))$SIGNATURE_VARIANT, # Alpha_EPI_ISL_710528
#   C29_Delta = (all_the_variants %>% filter(SIGNATURE_POS >= 21563) %>% filter(SIGNATURE_POS <= 25384)  %>% filter(QUERY_TAG == "EPI_ISL_2693246"))$SIGNATURE_VARIANT, # Delta_EPI_ISL_2693246
#   C48_Omicron = (all_the_variants %>% filter(SIGNATURE_POS >= 21563) %>% filter(SIGNATURE_POS <= 25384)  %>% filter(QUERY_TAG == "EPI_ISL_6841980"))$SIGNATURE_VARIANT # Omicron_EPI_ISL_6841980
# )
# 
# ggVennDiagram(venn_twist, label_alpha = 0) + scale_x_continuous(expand = expansion(mult = .2)) + scale_fill_distiller(palette = "RdBu") + scale_color_brewer(palette = "Set3") + labs(title = NULL,
#                                                                                                                                                                                       subtitle = NULL,
#                                                                                                                                                                                       caption = "Spike Protein Positions Mapped against MN908947_3")
# 
# file <- paste0("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/GILOT_THESIS/figures/","twist_controls_spike_protein.png")
# ggsave(file, device = "png", dpi = "retina", type = "cairo")
# 
# 
# C1_Australia <- (all_the_variants %>% filter(SIGNATURE_POS >= 21563) %>% filter(SIGNATURE_POS <= 25384)  %>% filter(QUERY_TAG == "MT007544.1"))$SIGNATURE_VARIANT # Australia_MT007544_1
# C14_Alpha <- (all_the_variants %>% filter(SIGNATURE_POS >= 21563) %>% filter(SIGNATURE_POS <= 25384)  %>% filter(QUERY_TAG == "EPI_ISL_710528"))$SIGNATURE_VARIANT # Alpha_EPI_ISL_710528
# C29_Delta <- (all_the_variants %>% filter(SIGNATURE_POS >= 21563) %>% filter(SIGNATURE_POS <= 25384)  %>% filter(QUERY_TAG == "EPI_ISL_2693246"))$SIGNATURE_VARIANT # Delta_EPI_ISL_2693246
# C48_Omicron  <- (all_the_variants %>% filter(SIGNATURE_POS >= 21563) %>% filter(SIGNATURE_POS <= 25384)  %>% filter(QUERY_TAG == "EPI_ISL_6841980"))$SIGNATURE_VARIANT #EPI_ISL_6841980
# 
# spike_C14 <- intersect(only_in_signature_EPI_ISL_710528, C14_Alpha)
# spike_C1 <-intersect(only_in_signature_MT007544_1, C1_Australia)
# spike_C48 <-intersect(only_in_signature_EPI_ISL_6841980, C48_Omicron)
# spike_C29 <-intersect(only_in_signature_EPI_ISL_2693246, C29_Delta)
# 
# 


# 
# ggVennDiagram(venn_twist, fill_color = c("#0073C2FF", "#EFC000FF", "#868686FF", "#CD534CFF"))
# 
# 
# display_venn(venn_twist, category.names = c("Omicron" , "Delta" , "Alpha", "Australia"),
#              fill = c("#999999", "#E69F00", "#56B4E9", "#009E73"))
# 
# v.table <- venn(venn_twist)