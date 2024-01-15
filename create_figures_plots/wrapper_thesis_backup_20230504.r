names(sessionInfo()$loadedOnly)

# rm(list = ls(all = TRUE))
# lapply(names(sessionInfo()$loadedOnly), require, character.only = TRUE)
# invisible(lapply(paste0('package:', names(sessionInfo()$otherPkgs)), detach, character.only=TRUE, unload=TRUE, force=TRUE))

rm(list = ls(all.names = TRUE))
gc()
#.rs.restartR()

#library(DescTools)
#library(Rsamtools)
#library(doParallel)
#library(foreach)
#library(vcfR)
library(RColorBrewer)
library(VennDiagram)
library(genpwr)
library(ggVennDiagram)
library(ggplot2)
library(ggthemes)
library(ggvenn)
#library(ggsignif)
library(gplots)
library(ggpmisc)
library(kableExtra)
library(lubridate)
library(pegas)
library(readr)
library(scales)
library(splitstackshape)
library(stringr)
library(tidyverse)
library(vroom)
library(xtable)
library(data.table)
# library(memisc)
#library(conflicted)

#conflicted::conflicts_prefer(dplyr::filter)

folder_data_unified <- "/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/"

folder_plots <- "/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/GILOT_THESIS/folder_plots/"
folder_tables <- "/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/GILOT_THESIS/folder_tables/"

all_uniq_variants <- readRDS(file = paste0("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/", "all_uniq_variants", ".rds"))
all_the_variants <- readRDS(file = paste0("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/", "all_the_variants", ".rds"))

SNPs_nucmer_MT007544_1_australia <- readRDS("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_computed/SNPs_MT007544_1_australia.rds")
SNPs_nucmer_EPI_ISL_omicron <- readRDS("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_computed/SNPs_EPI_ISL_omicron.rds")
SNPs_nucmer_EPI_ISL_2693246_delta <- readRDS("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_computed/SNPs_EPI_ISL_2693246_delta.rds")
SNPs_nucmer_EPI_ISL_710528_alpha <- readRDS("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_computed/SNPs_EPI_ISL_710528_alpha.rds")

only_in_signature_EPI_ISL_710528 <- readRDS(file = paste0("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/", "variant_is_in_signature_EPI_ISL_710528", ".rds"))
only_in_signature_MT007544_1 <- readRDS(file = paste0("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/", "variant_is_in_signature_MT007544_1", ".rds"))
only_in_signature_EPI_ISL_6841980 <- readRDS(file = paste0("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/", "variant_is_in_signature_EPI_ISL_6841980", ".rds"))
only_in_signature_EPI_ISL_2693246 <- readRDS(file = paste0("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/", "variant_is_in_signature_EPI_ISL_2693246", ".rds"))

SNPs_nucmer_from_synthetic_controls <- readRDS("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/SNPs_nucmer_from_synthetic_controls.rds")
features_gff_MN908947_3 <- readRDS("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/gff_MN908947_3.rds")

appendix_features <- vroom("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/gff_MN908947_3.csv") %>% select(start, end, THESIS_FEATURE_PRODUCT)
appendix_features$start <- appendix_features$start %>% as.integer()
appendix_features$end <- appendix_features$end %>% as.integer()
appendix_features$Bases <- as.integer(appendix_features$end - appendix_features$start)
appendix_features$SNPs <- as.integer(NA)
appendix_features$SNP_per_Kb  <- as.numeric(NA)

#colnames(SNPs_nucmer_from_synthetic_controls)
uniq_snp <- SNPs_nucmer_from_synthetic_controls %>% select(REFERENCE_POSITION, SIGNATURE_VARIANT) %>% unique()  # 108 total uniq from 123

for (row in seq(nrow(appendix_features))) {
  print(row)
  print(appendix_features[row,]$start)
  print(appendix_features[row,]$end)
  uniq_snp %>% filter(REFERENCE_POSITION >= appendix_features[row,]$start) %>% filter(REFERENCE_POSITION <= appendix_features[row,]$end) %>% print()
  uniq_snp %>% filter(REFERENCE_POSITION >= appendix_features[row,]$start) %>% filter(REFERENCE_POSITION <= appendix_features[row,]$end) %>% length() %>% print()
  #appendix_features[row,]$SNPs <- as.integer()
  appendix_features[row,]$SNPs <- uniq_snp %>% filter(REFERENCE_POSITION >= appendix_features[row,]$start) %>% filter(REFERENCE_POSITION <= appendix_features[row,]$end) %>% nrow()
}

names(appendix_features)[names(appendix_features) == 'start'] <- 'Start'
names(appendix_features)[names(appendix_features) == 'end'] <- 'End'
names(appendix_features)[names(appendix_features) == 'THESIS_FEATURE_PRODUCT'] <- 'Product'
names(appendix_features)[names(appendix_features) == 'SNPs'] <- 'SNPs'

appendix_features$SNP_per_Kb <- appendix_features$SNPs / appendix_features$Bases * 1000

colnames(appendix_features) 
appendix_features <-
  appendix_features %>% select("Start", "End", "Bases", "Product", "SNPs", "SNP_per_Kb")

# table: genome_features --------------------------------------------------
colnames(appendix_features)
file_table <- paste0(folder_tables,"genome_features", ".tex")
print(file_table)
sink(file_table)
print(xtable(appendix_features), include.rownames=FALSE)
sink()

# SNPS CONTROLS with controls data ------------------------------------------------------
all_the_variants <- readRDS(file = paste0("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/", "all_the_variants", ".rds"))
list_control <- all_the_variants$QUERY_TAG %>% unique() %>% sort()

for (iter_control in list_control) {
  appendix_features[[iter_control]] <- as.integer(NA)
}

for (row_in_df in seq(1,nrow(appendix_features))) {
  
  for (iter_control in list_control) {
    #print(paste(iter_control, appendix_features[row_in_df,]$Start, appendix_features[row_in_df,]$End))
    #print(all_the_variants %>% filter(SIGNATURE_POS >= appendix_features[row_in_df,]$Start) %>% filter(SIGNATURE_POS <= appendix_features[row_in_df,]$End))
    paste(iter_control)
    how_many <- all_the_variants %>% filter(QUERY_TAG == iter_control) %>% filter(SIGNATURE_POS >= appendix_features[row_in_df,]$Start) %>% filter(SIGNATURE_POS <= appendix_features[row_in_df,]$End) %>% nrow()
    how_many_start <- appendix_features[row_in_df,]$Start
    how_many_end <- appendix_features[row_in_df,]$End
    print(paste(iter_control, how_many_start, how_many_end, how_many))
    appendix_features[row_in_df,][[iter_control]] <- how_many %>% as.integer()
  }
}

appendix_features_controls <- appendix_features %>% select(Start, End, Bases, Product, MT007544.1, EPI_ISL_710528, EPI_ISL_2693246, EPI_ISL_6841980, SNPs, SNP_per_Kb)
colnames(appendix_features_controls) <- c("Start", "End", "Bases", "Product", "Australia", "Alpha", "Delta", "Omicron", "Unique SNPs", "SNP/Kb")

file_table <- paste0(folder_tables,"appendix_features_controls", ".tex")
print(file_table)
sink(file_table)
print(xtable(appendix_features_controls), include.rownames=FALSE)
sink()

# Table: Synthetic Viral Controls -----------------------------------------
twist <- vroom("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/twist_part_numbers_clean.csv")
twist_clean <- twist %>% filter(PART_NUMBER %in% c("102019", "102024", "103907", "104539", "105204"))
twist_clean$PART_NUMBER <- as.integer(twist_clean$PART_NUMBER)
colnames(twist_clean) <- c("Part #","Control","GenBank GISAID","GISAID NAME")
file_table <- paste0(folder_tables,"synthetic_viral_controls", ".tex")
print(file_table)
sink(file_table)
print(xtable(twist_clean), include.rownames=FALSE)
sink()

# Table: proiminent variants ----------------------------------------------
prom_variants <- vroom("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/prominent_variants.csv") %>% arrange(`Name by WHO`)
colnames(prom_variants) <- c("Variant","WHO Designation","Country of Origin", "RBD", "S-Glycoprotein" )
file_table <- paste0(folder_tables,"prominent_variants", ".tex")
print(file_table)
sink(file_table)
print(xtable(prom_variants), include.rownames=FALSE)
sink()

# mutations_in_controls ---------------------------------------------------
venn_twist <- list(
  C1_Australia = (all_the_variants %>% filter(QUERY_TAG == "MT007544.1"))$SIGNATURE_VARIANT, # Australia_MT007544_1
  C14_Alpha = (all_the_variants %>% filter(QUERY_TAG == "EPI_ISL_710528"))$SIGNATURE_VARIANT, # Alpha_EPI_ISL_710528
  C29_Delta = (all_the_variants %>% filter(QUERY_TAG == "EPI_ISL_2693246"))$SIGNATURE_VARIANT, # Delta_EPI_ISL_2693246
  C48_Omicron = (all_the_variants %>% filter(QUERY_TAG == "EPI_ISL_6841980"))$SIGNATURE_VARIANT # Omicron_EPI_ISL_6841980
)

plot_OBJECT <- ggVennDiagram(venn_twist, label_alpha = 0) + scale_x_continuous(expand = expansion(mult = .2)) + scale_fill_distiller(palette = "RdBu") + scale_color_brewer(palette = "Set3") + labs(title = NULL,subtitle = NULL,caption = "Mapped against MN908947_3")
file <- paste0(folder_plots,"mutations_in_controls", ".png")
ggsave(file, device = "png", dpi = 1800,  type = "cairo") # , dpi = "retina",  type = "cairo") # width = 210, height = 148.5 , units = c("mm")

# table: depth_post_mapping -----------------------------------------------
depths_blasted <- vroom("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/sample_depth_statistics_from_mpileups_gini_cohorted.csv") %>% filter(STAT_SAMPLE_METHOD == "blasted") %>% select(STAT_IDENTIFIER, STAT_DEPTH_MIN ,STAT_DEPTH_MEDIAN, STAT_DEPTH_MEAN, STAT_DEPTH_MAX) %>% arrange(STAT_IDENTIFIER) %>% unique()
colnames(depths_blasted) <- c("Identifier", "Min_Samblaster", "Median_Samblaster", "Mean_Samblaster", "Max_Samblaster")

depths_unblasted <- vroom("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/sample_depth_statistics_from_mpileups_gini_cohorted.csv") %>% filter(STAT_SAMPLE_METHOD == "unblasted") %>% select(STAT_IDENTIFIER , STAT_DEPTH_MIN, STAT_DEPTH_MEDIAN, STAT_DEPTH_MEAN, STAT_DEPTH_MAX) %>% arrange(STAT_IDENTIFIER) %>% unique()
colnames(depths_unblasted) <- c("Identifier", "Min_DeDup", "Median_DeDup", "Mean_DeDup", "Max_DeDup")

depth_post_mapping <- left_join(depths_blasted, depths_unblasted) %>% arrange(Identifier)
depth_post_mapping$Identifier <- as.character(map(strsplit(depth_post_mapping$Identifier, split = "_"), 1))

colnames(depth_post_mapping) <- c("Identifier", "Min_blind", "Median_blind", "Mean_blind", "Max_blind", "Min_aware", "Median_aware", "Mean_aware", "Max_aware")
depth_post_mapping$sort_mean <- rowMeans(cbind(depth_post_mapping$Min_blind, depth_post_mapping$Median_blind, depth_post_mapping$Mean_blind, depth_post_mapping$Max_blind, depth_post_mapping$Min_aware, depth_post_mapping$Median_aware, depth_post_mapping$Mean_aware, depth_post_mapping$Max_aware)) %>% as.numeric() %>% round(1)

depth_post_mapping <- depth_post_mapping %>% arrange(desc(sort_mean))
glimpse(depth_post_mapping)
depth_post_mapping
file_table <- paste0(folder_tables,"depth_post_mapping", ".tex")
print(file_table)
sink(file_table)
print(xtable(depth_post_mapping), floating = FALSE, tabular.environment="longtable", include.rownames=FALSE)
sink()


# table: human_samples ----------------------------------------------------
pangolin_lineage_appendix <-
  vroom(
    "/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/metadata_new.csv"
  ) %>% select(METADATA_ID,
               M_ABNAHMEDATUM,
               METADATA_PL_LINEAGE,
               METADATA_CT_S,
               METADATA_CT_E) %>% filter(
                 !is.na(METADATA_PL_LINEAGE),
                 !is.na(M_ABNAHMEDATUM),
                 !is.na(METADATA_CT_S),
                 !is.na(METADATA_CT_E)
               ) %>% arrange(METADATA_CT_S, METADATA_CT_E, M_ABNAHMEDATUM)

names(pangolin_lineage_appendix) <- c("Identifier", "Collection", "Lineage", "CT_S", "CT_E")

pangolin_lineage_appendix$Collection <- pangolin_lineage_appendix$Collection %>% ymd() %>% as.character()
print(xtable(pangolin_lineage_appendix), include.rownames=FALSE)

file_table <- paste0(folder_tables,"pangolin_lineage_appendix", ".tex")
print(file_table)
sink(file_table)
print(xtable(pangolin_lineage_appendix), floating = FALSE, tabular.environment="longtable", include.rownames=FALSE)
sink()


# ct_vs_mean_depth --------------------------------------------------------
ct_vs_mean_depth <- left_join(pangolin_lineage_appendix, depth_post_mapping) %>% filter(!is.na(Median_blind)) %>% filter(!is.na(Median_aware)) %>% arrange(CT_S, CT_E, desc(sort_mean))
glimpse(ct_vs_mean_depth)
human_samples <- ct_vs_mean_depth %>% arrange(desc(sort_mean)) %>% select(Identifier, CT_S, Median_blind, Median_aware)
colnames(human_samples) <- c("Identifier", "CT(S) Value", "Median_UMI_blind", "Median_UMI_aware")
file_table <- paste0(folder_tables,"human_samples", ".tex")
print(file_table)
sink(file_table)
print(xtable(human_samples), floating = FALSE, tabular.environment="longtable", include.rownames=FALSE)
sink()

# plot: cts vs cte R2--------------------------------------------------------
label_x <- "CT(S) Value"
label_y <- "CT(E) Value"
label_caption <- NULL
plot_OBJECT <- ct_vs_mean_depth %>% ggplot(aes(x = CT_S, y = CT_E))  +
  stat_poly_line(    linewidth = 0.5,
                     alpha=0.5,
                     color = "cyan") +
  stat_poly_eq() +
  geom_point(alpha = 15/20, size = 1.5, color='blue') + 
  scale_y_continuous(breaks = seq(7.5, 25, by = 2.5), labels = label_comma(), expand = c(0, 0)) +
  scale_x_continuous(breaks = seq(7.5, 25, by = 2.5), labels = label_comma(), expand = c(0, 0)) +
  labs(x = label_x) +
  labs(y = label_y) +
  labs(caption = label_caption)

file <- paste0(folder_plots,"cts_vs_cte_r2", ".png")
ggsave(file, device = "png", dpi = 1800,  type = "cairo") 

# plot: cts vs cte --------------------------------------------------------
label_x <- "CT(S) Value"
label_y <- "CT(E) Value"
label_caption <- NULL
plot_OBJECT <- ct_vs_mean_depth %>% ggplot(aes(x = CT_S, y = CT_E))  +
  geom_point(alpha = 15/20, size = 1.5, color='blue') +
  stat_poly_eq() +
  geom_smooth(
    method = "auto",
    se = TRUE,
    fullrange = FALSE,
    level = 0.95,
    linewidth = 0.5,
    alpha=0.5,
    color = "cyan"
  ) + 
  scale_y_continuous(breaks = seq(7.5, 25, by = 2.5), labels = label_comma(), expand = c(0, 0)) +
  scale_x_continuous(breaks = seq(7.5, 25, by = 2.5), labels = label_comma(), expand = c(0, 0)) +
  labs(x = label_x) +
  labs(y = label_y) +
  labs(caption = label_caption)

file <- paste0(folder_plots,"cts_vs_cte", ".png")
ggsave(file, device = "png", dpi = 1800,  type = "cairo") 

# plot: cts UMI blind mapping -------------------------------------------------
label_x <- "CT(S) Value"
label_y <- "Median Depth"
label_caption <- "UMI blind mapping"
plot_OBJECT <- ct_vs_mean_depth %>% ggplot(aes(x = CT_S, y = Median_blind))  +
  geom_point(alpha = 15/20, size = 1.5, color='blue') +
  geom_smooth(
    method = "auto",
    se = TRUE,
    fullrange = FALSE,
    level = 0.95,
    linewidth = 0.5,
    alpha=0.5,
    color = "cyan"
  ) + 
  scale_y_continuous(labels = unit_format(unit = "k", scale = 1e-3), expand = c(0, 0)) +
  scale_x_continuous(breaks = seq(7.5, 25, by = 2.5), labels = label_comma(), expand = c(0, 0)) +
  labs(x = label_x) +
  labs(y = label_y) +
  labs(caption = label_caption)

file <- paste0(folder_plots,"cts_umi_blind_mapping_depth", ".png")
ggsave(file, device = "png", dpi = 1800,  type = "cairo") 


# plot: cts UMI aware mapping -------------------------------------------------
label_x <- "CT(S) Value"
label_y <- "Median Depth"
label_caption <- "UMI aware mapping"
plot_OBJECT <- ct_vs_mean_depth %>% ggplot(aes(x = CT_S, y = Median_aware))  +
  geom_point(alpha = 15/20, size = 1.5, color='blue') +
  geom_smooth(
    method = "auto",
    se = TRUE,
    fullrange = FALSE,
    level = 0.95,
    linewidth = 0.5,
    alpha=0.5,
    color = "cyan"
  ) + 
  scale_y_continuous(labels = unit_format(unit = "k", scale = 1e-3), expand = c(0, 0)) +
  scale_x_continuous(breaks = seq(7.5, 25, by = 2.5), labels = label_comma(), expand = c(0, 0)) +
  labs(x = label_x) +
  labs(y = label_y) +
  labs(caption = label_caption)

file <- paste0(folder_plots,"cts_umi_aware_mapping_depth", ".png")
ggsave(file, device = "png", dpi = 1800,  type = "cairo") 


# plot: cts UMI blind mapping -------------------------------------------------
label_x <- "CT(E) Value"
label_y <- "Median Depth"
label_caption <- "UMI blind mapping"
plot_OBJECT <- ct_vs_mean_depth %>% ggplot(aes(x = CT_E, y = Median_blind))  +
  geom_point(alpha = 15/20, size = 1.5, color='blue') +
  geom_smooth(
    method = "auto",
    se = TRUE,
    fullrange = FALSE,
    level = 0.95,
    linewidth = 0.5,
    alpha=0.5,
    color = "cyan"
  ) + 
  scale_y_continuous(labels = unit_format(unit = "k", scale = 1e-3), expand = c(0, 0)) +
  scale_x_continuous(breaks = seq(7.5, 25, by = 2.5), labels = label_comma(), expand = c(0, 0)) +
  labs(x = label_x) +
  labs(y = label_y) +
  labs(caption = label_caption)

file <- paste0(folder_plots,"cte_umi_blind_mapping_depth", ".png")
ggsave(file, device = "png", dpi = 1800,  type = "cairo") 


# plot: cts UMI aware mapping -------------------------------------------------
label_x <- "CT(E) Value"
label_y <- "Median Depth"
label_caption <- "UMI aware mapping"
plot_OBJECT <- ct_vs_mean_depth %>% ggplot(aes(x = CT_E, y = Median_aware))  +
  geom_point(alpha = 15/20, size = 1.5, color='blue') +
  geom_smooth(
    method = "auto",
    se = TRUE,
    fullrange = FALSE,
    level = 0.95,
    linewidth = 0.5,
    alpha=0.5,
    color = "cyan"
  ) + 
  scale_y_continuous(labels = unit_format(unit = "k", scale = 1e-3), expand = c(0, 0)) +
  scale_x_continuous(breaks = seq(7.5, 25, by = 2.5), labels = label_comma(), expand = c(0, 0)) +
  labs(x = label_x) +
  labs(y = label_y) +
  labs(caption = label_caption)

file <- paste0(folder_plots,"cte_umi_aware_mapping_depth", ".png")
ggsave(file, device = "png", dpi = 1800,  type = "cairo") 



# NEW --------------------------------------------------------------------

synthetic_combinations_fractions <- vroom("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/synthetic_combinations_fractions.csv")
synthetic_combinations_fractions <- synthetic_combinations_fractions %>%
  select(THESIS_IDENTIFIER, THESIS_SYNTHETIC_COPIES_TOTAL, ends_with("_FRACTION")) %>%
  arrange(desc(THESIS_SYNTHETIC_COPIES_TOTAL)) %>%
  rename('Identifier' = 'THESIS_IDENTIFIER') %>% 
  rename('Viral_Copies' = 'THESIS_SYNTHETIC_COPIES_TOTAL') 
synthetic_combinations_fractions$Identifier <- as.character(map(strsplit(synthetic_combinations_fractions$Identifier, split = "_"), 1))
synthetic_combinations_fractions$Viral_Copies <- synthetic_combinations_fractions$Viral_Copies %>% as.integer()
synthetic_combinations_fractions[,3:6] <- synthetic_combinations_fractions[,3:6] %>% round(6)

synthetic_combinations_fractions$THESIS_SYNTHETIC_PRESENT_C1_AUSTRALIA_MT007544_1_FRACTION_COPIES <- synthetic_combinations_fractions$Viral_Copies * synthetic_combinations_fractions$THESIS_SYNTHETIC_PRESENT_C1_AUSTRALIA_MT007544_1_FRACTION
synthetic_combinations_fractions$THESIS_SYNTHETIC_PRESENT_C2_ANCESTRAL_MN908947_3_FRACTION_COPIES <- synthetic_combinations_fractions$Viral_Copies * synthetic_combinations_fractions$THESIS_SYNTHETIC_PRESENT_C2_ANCESTRAL_MN908947_3_FRACTION

synthetic_combinations_fractions$THESIS_SYNTHETIC_PRESENT_C14_ALPHA_EPI_ISL_710528_FRACTION_COPIES <- synthetic_combinations_fractions$Viral_Copies * synthetic_combinations_fractions$THESIS_SYNTHETIC_PRESENT_C14_ALPHA_EPI_ISL_710528_FRACTION
synthetic_combinations_fractions$THESIS_SYNTHETIC_PRESENT_C29_DELTA_EPI_ISL_2693246_FRACTION_COPIES <- synthetic_combinations_fractions$Viral_Copies * synthetic_combinations_fractions$THESIS_SYNTHETIC_PRESENT_C29_DELTA_EPI_ISL_2693246_FRACTION
synthetic_combinations_fractions$THESIS_SYNTHETIC_PRESENT_C48_OMICRON_EPI_ISL_6841980_FRACTION_COPIES <- synthetic_combinations_fractions$Viral_Copies * synthetic_combinations_fractions$THESIS_SYNTHETIC_PRESENT_C48_OMICRON_EPI_ISL_6841980_FRACTION

glimpse(synthetic_combinations_fractions)

colnames(synthetic_combinations_fractions) <- c(
  "Identifier",
  "Total Copies",
  "C1 Fraction",
  "C2 Fraction",
  "C14 Fraction",
  "C29 Fraction",
  "C48 Fraction",
  "C1 Copies",
  "C2 Copies",
  "C14 Copies",
  "C29 Copies",
  "C48 Copies"
)

glimpse(synthetic_combinations_fractions)
synthetic_samples <- left_join(synthetic_combinations_fractions, depth_post_mapping)
glimpse(synthetic_samples)

synthetic_samples <- synthetic_samples %>%
  select(
    Identifier,
    Median_blind,
    Median_aware,
    `Total Copies`,
    `C1 Fraction`,
    `C2 Fraction`,
    `C14 Fraction`,
    `C29 Fraction`,
    `C48 Fraction`,
    `C1 Copies`,
    `C2 Copies`,
    `C14 Copies`,
    `C29 Copies`,
    `C48 Copies`,
    sort_mean
  ) %>% arrange(desc(sort_mean))

glimpse(synthetic_samples)

synthetic_samples_all <- synthetic_samples %>% arrange(desc(`Total Copies`), desc(`C1 Copies`), desc(`C2 Copies`), desc(`C14 Copies`), desc(`C29 Copies`), desc(`C48 Copies`)) 
data_fwrite <- synthetic_samples_all
name_data_frwite <- "synthetic_samples_all"
new_filename <- paste0(folder_data_unified, name_data_frwite, ".csv")
print(new_filename)
fwrite(data_fwrite, file = new_filename, append = FALSE, quote = "auto", showProgress = getOption("datatable.showProgress", interactive()), buffMB = 1024, nThread = 6)

synthetic_samples_C1 <- synthetic_samples %>% filter(`C1 Copies` > 0) %>% arrange(desc(`C1 Copies`), desc(`C2 Copies`), desc(`C14 Copies`), desc(`C29 Copies`), desc(`C48 Copies`)) 
data_fwrite <- synthetic_samples_C1
name_data_frwite <- "synthetic_samples_C1"
new_filename <- paste0(folder_data_unified, name_data_frwite, ".csv")
print(new_filename)
fwrite(data_fwrite, file = new_filename, append = FALSE, quote = "auto", showProgress = getOption("datatable.showProgress", interactive()), buffMB = 1024, nThread = 6)

synthetic_samples_C2 <- synthetic_samples %>% filter(`C2 Copies` > 0) %>% arrange(desc(`C2 Copies`), desc(`C1 Copies`), desc(`C14 Copies`), desc(`C29 Copies`), desc(`C48 Copies`)) 
data_fwrite <- synthetic_samples_C2
name_data_frwite <- "synthetic_samples_C2"
new_filename <- paste0(folder_data_unified, name_data_frwite, ".csv")
print(new_filename)
fwrite(data_fwrite, file = new_filename, append = FALSE, quote = "auto", showProgress = getOption("datatable.showProgress", interactive()), buffMB = 1024, nThread = 6)

synthetic_samples_C14 <- synthetic_samples %>% filter(`C14 Copies` > 0) %>% arrange(desc(`C14 Copies`), desc(`C1 Copies`), desc(`C2 Copies`), desc(`C29 Copies`), desc(`C48 Copies`)) 
data_fwrite <- synthetic_samples_C14
name_data_frwite <- "synthetic_samples_C14"
new_filename <- paste0(folder_data_unified, name_data_frwite, ".csv")
print(new_filename)
fwrite(data_fwrite, file = new_filename, append = FALSE, quote = "auto", showProgress = getOption("datatable.showProgress", interactive()), buffMB = 1024, nThread = 6)

synthetic_samples_C29 <- synthetic_samples %>% filter(`C29 Copies` > 0) %>% arrange(desc(`C29 Copies`), desc(`C1 Copies`), desc(`C2 Copies`), desc(`C14 Copies`), desc(`C48 Copies`)) 
data_fwrite <- synthetic_samples_C29
name_data_frwite <- "synthetic_samples_C29"
new_filename <- paste0(folder_data_unified, name_data_frwite, ".csv")
print(new_filename)
fwrite(data_fwrite, file = new_filename, append = FALSE, quote = "auto", showProgress = getOption("datatable.showProgress", interactive()), buffMB = 1024, nThread = 6)

synthetic_samples_C48 <- synthetic_samples %>% filter(`C48 Copies` > 0) %>% arrange(desc(`C48 Copies`), desc(`C1 Copies`), desc(`C2 Copies`), desc(`C14 Copies`), desc(`C29 Copies`)) 
data_fwrite <- synthetic_samples_C48
name_data_frwite <- "synthetic_samples_C48"
new_filename <- paste0(folder_data_unified, name_data_frwite, ".csv")
print(new_filename)
fwrite(data_fwrite, file = new_filename, append = FALSE, quote = "auto", showProgress = getOption("datatable.showProgress", interactive()), buffMB = 1024, nThread = 6)

synthetic_samples_to_xtable <- synthetic_samples
synthetic_samples_to_xtable <- synthetic_samples_to_xtable %>% arrange(desc(`Total Copies`), desc(`C1 Copies`), desc(`C2 Copies`), desc(`C14 Copies`), desc(`C29 Copies`), desc(`C48 Copies`)) 

synthetic_samples_to_xtable$`C1 Fraction` <- synthetic_samples_to_xtable$`C1 Fraction` %>% round(6)  %>% as.character()
synthetic_samples_to_xtable$`C2 Fraction` <- synthetic_samples_to_xtable$`C2 Fraction` %>% round(6)  %>% as.character()
synthetic_samples_to_xtable$`C14 Fraction` <- synthetic_samples_to_xtable$`C14 Fraction` %>% round(6)  %>% as.character() 
synthetic_samples_to_xtable$`C29 Fraction` <- synthetic_samples_to_xtable$`C29 Fraction` %>% round(6)  %>% as.character()
synthetic_samples_to_xtable$`C48 Fraction` <- synthetic_samples_to_xtable$`C48 Fraction` %>% round(6)  %>% as.character()

synthetic_samples_to_xtable$Median_blind <- synthetic_samples_to_xtable$Median_blind %>% as.integer() %>% as.character()
synthetic_samples_to_xtable$Median_aware <- synthetic_samples_to_xtable$Median_aware %>% as.integer() %>% as.character()
synthetic_samples_to_xtable$`Total Copies` <- synthetic_samples_to_xtable$`Total Copies` %>% as.integer() %>% as.character()

synthetic_samples_to_xtable$`C1 Copies` <- synthetic_samples_to_xtable$`C1 Copies` %>% as.integer() %>% as.character()
synthetic_samples_to_xtable$`C2 Copies` <- synthetic_samples_to_xtable$`C2 Copies` %>% as.integer() %>% as.character()
synthetic_samples_to_xtable$`C14 Copies` <- synthetic_samples_to_xtable$`C14 Copies` %>% as.integer() %>% as.character()
synthetic_samples_to_xtable$`C29 Copies` <- synthetic_samples_to_xtable$`C29 Copies` %>% as.integer() %>% as.character()
synthetic_samples_to_xtable$`C48 Copies` <- synthetic_samples_to_xtable$`C48 Copies` %>% as.integer() %>% as.character()

colnames(synthetic_samples_to_xtable)
glimpse(synthetic_samples)
file_table <- paste0(folder_tables,"synthetic_samples", ".tex")
print(file_table)
sink(file_table)
print(xtable(synthetic_samples_to_xtable), include.rownames=FALSE)
sink()









depths_blasted <- vroom("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/sample_depth_statistics_from_mpileups_gini_cohorted.csv") %>% 
  filter(STAT_SAMPLE_METHOD == "blasted") %>% 
  select(STAT_IDENTIFIER, STAT_DEPTH_MEDIAN, STAT_DEPTH_MEAN, STAT_DEPTH_MAX) %>% 
  unique() %>%
  arrange(STAT_IDENTIFIER)
colnames(depths_blasted) <- c("Identifier", "Median_blind", "Mean_blind", "Max_blind")
glimpse(depths_blasted)

depths_unblasted <- vroom("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/sample_depth_statistics_from_mpileups_gini_cohorted.csv") %>% 
  filter(STAT_SAMPLE_METHOD == "unblasted") %>% 
  select(STAT_IDENTIFIER , STAT_DEPTH_MEDIAN, STAT_DEPTH_MEAN, STAT_DEPTH_MAX) %>% 
  unique() %>%
  arrange(STAT_IDENTIFIER)
colnames(depths_unblasted) <- c("Identifier", "Median_aware", "Mean_aware", "Max_aware")
glimpse(depths_unblasted)

depths_by_method <- left_join(depths_blasted, depths_unblasted)

depth_post_mapping$sort_mean <- rowMeans(cbind(depth_post_mapping$Min_blind, depth_post_mapping$Median_blind, depth_post_mapping$Mean_blind, depth_post_mapping$Max_blind, depth_post_mapping$Min_aware, depth_post_mapping$Median_aware, depth_post_mapping$Mean_aware, depth_post_mapping$Max_aware)) %>% as.numeric() %>% round(1)

glimpse(depths_by_method)


colnames(depths_blasted) <- c("Identifier", "Median_Samblaster", "Mean_Samblaster", "Max_SamBlaster")

rename('STAT_IDENTIFIER' = 'Identifier', 'STAT_DEPTH_MEDIAN' = 'Median_DeDup', 'STAT_DEPTH_MEAN' = 'Mean_DeDup', 'STAT_DEPTH_MAX' = 'Max_DeDup') %>% 
  arrange(Identifier)

depths_by_method <- left_join(depths_blasted, depths_unblasted)

depths_by_method <- left_join(depths_blasted, depths_unblasted)

depths_by_method$Median_Samblaster <- depths_by_method$Median_Samblaster %>% as.integer()
depths_by_method$Max_SamBlaster <- depths_by_method$Max_SamBlaster %>% as.integer()
depths_by_method$Median_DeDup <- depths_by_method$Median_DeDup %>% as.integer()
depths_by_method$Max_DeDup <- depths_by_method$Max_DeDup %>% as.integer()

print(xtable(depths_by_method), include.rownames=FALSE)







# noname yet --------------------------------------------------------------

pangolin_lineage <- vroom("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/metadata_new.csv") %>% select(METADATA_PL_LINEAGE) %>% arrange(METADATA_PL_LINEAGE) %>% unique()
colnames(pangolin_lineage) <- c("Lineage")

all_vcfs_labelled <- readRDS("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/all_callers_with_features_and_metadata.rds")
# observed_experimental <- all_vcfs_labelled %>% filter(THESIS_COHORT == "experimental") %>% filter(FILTER == "PASS") %>% filter(POS >= 21563) %>% filter(POS <= 25384) %>% select(POS, THESIS_VARIANT) %>% unique() %>% arrange(POS, THESIS_VARIANT)
# observed_control <- all_vcfs_labelled %>% filter(THESIS_COHORT == "control") %>% filter(FILTER == "PASS") %>% filter(POS >= 21563) %>% filter(POS <= 25384) %>% select(POS, THESIS_VARIANT) %>% unique() %>% arrange(POS, THESIS_VARIANT)
observed_experimental <- all_vcfs_labelled %>% filter(THESIS_COHORT == "experimental") %>% filter(FILTER == "PASS") %>% select(POS, THESIS_VARIANT) %>% unique() %>% arrange(POS, THESIS_VARIANT)
observed_control <- all_vcfs_labelled %>% filter(THESIS_COHORT == "control") %>% filter(FILTER == "PASS") %>% select(POS, THESIS_VARIANT) %>% unique() %>% arrange(POS, THESIS_VARIANT)
# 
# four_files <- vroom("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/deep_assembly_visual.csv") %>% select(STAT_IDENTIFIER, Description, STAT_SAMPLE_METHOD, STAT_DEPTH_MIN, STAT_DEPTH_MEDIAN, STAT_DEPTH_MEAN, STAT_DEPTH_MAX)
# four_files$STAT_DEPTH_MIN <- four_files$STAT_DEPTH_MIN %>% as.integer()
# four_files$STAT_DEPTH_MEDIAN <- four_files$STAT_DEPTH_MEDIAN %>% as.integer()
# four_files$STAT_DEPTH_MAX <- four_files$STAT_DEPTH_MAX %>% as.integer()
# 
# colnames(four_files)
# four_files %>% rename('Identifier' = 'STAT_IDENTIFIER') %>% rename('UMI' = 'STAT_SAMPLE_METHOD') %>% rename('Minimum' = 'STAT_DEPTH_MIN') %>% rename('Median' = 'STAT_DEPTH_MEDIAN') %>% rename('Mean' = 'STAT_DEPTH_MEAN') %>% rename('Maximum' = 'STAT_DEPTH_MAX') %>% arrange(desc(Identifier), UMI, Median)
# 
# four_files <- four_files %>% rename('STAT_IDENTIFIER' = 'Identifier') %>% rename('STAT_SAMPLE_METHOD' = 'Method') %>% rename('STAT_DEPTH_MIN' = 'Minimum') %>% rename('STAT_DEPTH_MEDIAN' = 'Median') %>% rename('STAT_DEPTH_MEAN' = 'Mean') %>% rename('STAT_DEPTH_MAX' = 'Maximum') %>% arrange(desc(Identifier), Method, Median)
# 
# print(xtable(four_files), include.rownames=FALSE)

synthetic_combinations_fractions <- vroom("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/synthetic_combinations_fractions.csv") %>% select(THESIS_IDENTIFIER, THESIS_SYNTHETIC_COPIES_TOTAL, ends_with("_FRACTION")) %>% arrange (THESIS_SYNTHETIC_COPIES_TOTAL, THESIS_SYNTHETIC_PRESENT_C1_AUSTRALIA_MT007544_1_FRACTION, THESIS_SYNTHETIC_PRESENT_C2_ANCESTRAL_MN908947_3_FRACTION, THESIS_SYNTHETIC_PRESENT_C14_ALPHA_EPI_ISL_710528_FRACTION, THESIS_SYNTHETIC_PRESENT_C29_DELTA_EPI_ISL_2693246_FRACTION, THESIS_SYNTHETIC_PRESENT_C48_OMICRON_EPI_ISL_6841980_FRACTION) %>% unique()
synthetic_combinations_copies <- vroom("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/synthetic_combinations_copies.csv")
colnames(synthetic_combinations_fractions)

