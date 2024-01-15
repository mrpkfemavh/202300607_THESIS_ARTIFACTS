library(RColorBrewer)
library(VennDiagram)
library(genpwr)
library(ggVennDiagram)
library(ggplot2)
library(ggthemes)
library(ggvenn)
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

folder_data_unified <- "/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/"
folder_plots <- "/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/GILOT_THESIS/folder_plots/"
folder_tables <- "/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/GILOT_THESIS/folder_tables/"

# Table: Synthetic Viral Controls -----------------------------------------
twist <- vroom("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/twist_part_numbers_clean.csv")
twist_clean <- twist %>% filter(PART_NUMBER %in% c("102019", "102024", "103907", "104539", "105204"))
twist_clean$PART_NUMBER <- as.integer(twist_clean$PART_NUMBER)
colnames(twist_clean) <- c("Part #","Control","GenBank GISAID","GISAID NAME")
file_table <- paste0(folder_tables,"mm_synthetic_viral_controls", ".tex")
print(file_table)
sink(file_table)
print(xtable(twist_clean), include.rownames=FALSE)
sink()

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
# 
# synthetic_samples_all <- synthetic_samples %>% arrange(desc(`Total Copies`), desc(`C1 Copies`), desc(`C2 Copies`), desc(`C14 Copies`), desc(`C29 Copies`), desc(`C48 Copies`)) 
# data_fwrite <- synthetic_samples_all
# name_data_frwite <- "synthetic_samples_all"
# new_filename <- paste0(folder_data_unified, name_data_frwite, ".csv")
# print(new_filename)
# fwrite(data_fwrite, file = new_filename, append = FALSE, quote = "auto", showProgress = getOption("datatable.showProgress", interactive()), buffMB = 1024, nThread = 6)
# 
# synthetic_samples_C1 <- synthetic_samples %>% filter(`C1 Copies` > 0) %>% arrange(desc(`C1 Copies`), desc(`C2 Copies`), desc(`C14 Copies`), desc(`C29 Copies`), desc(`C48 Copies`)) 
# data_fwrite <- synthetic_samples_C1
# name_data_frwite <- "synthetic_samples_C1"
# new_filename <- paste0(folder_data_unified, name_data_frwite, ".csv")
# print(new_filename)
# fwrite(data_fwrite, file = new_filename, append = FALSE, quote = "auto", showProgress = getOption("datatable.showProgress", interactive()), buffMB = 1024, nThread = 6)
# 
# synthetic_samples_C2 <- synthetic_samples %>% filter(`C2 Copies` > 0) %>% arrange(desc(`C2 Copies`), desc(`C1 Copies`), desc(`C14 Copies`), desc(`C29 Copies`), desc(`C48 Copies`)) 
# data_fwrite <- synthetic_samples_C2
# name_data_frwite <- "synthetic_samples_C2"
# new_filename <- paste0(folder_data_unified, name_data_frwite, ".csv")
# print(new_filename)
# fwrite(data_fwrite, file = new_filename, append = FALSE, quote = "auto", showProgress = getOption("datatable.showProgress", interactive()), buffMB = 1024, nThread = 6)
# 
# synthetic_samples_C14 <- synthetic_samples %>% filter(`C14 Copies` > 0) %>% arrange(desc(`C14 Copies`), desc(`C1 Copies`), desc(`C2 Copies`), desc(`C29 Copies`), desc(`C48 Copies`)) 
# data_fwrite <- synthetic_samples_C14
# name_data_frwite <- "synthetic_samples_C14"
# new_filename <- paste0(folder_data_unified, name_data_frwite, ".csv")
# print(new_filename)
# fwrite(data_fwrite, file = new_filename, append = FALSE, quote = "auto", showProgress = getOption("datatable.showProgress", interactive()), buffMB = 1024, nThread = 6)
# 
# synthetic_samples_C29 <- synthetic_samples %>% filter(`C29 Copies` > 0) %>% arrange(desc(`C29 Copies`), desc(`C1 Copies`), desc(`C2 Copies`), desc(`C14 Copies`), desc(`C48 Copies`)) 
# data_fwrite <- synthetic_samples_C29
# name_data_frwite <- "synthetic_samples_C29"
# new_filename <- paste0(folder_data_unified, name_data_frwite, ".csv")
# print(new_filename)
# fwrite(data_fwrite, file = new_filename, append = FALSE, quote = "auto", showProgress = getOption("datatable.showProgress", interactive()), buffMB = 1024, nThread = 6)
# 
# synthetic_samples_C48 <- synthetic_samples %>% filter(`C48 Copies` > 0) %>% arrange(desc(`C48 Copies`), desc(`C1 Copies`), desc(`C2 Copies`), desc(`C14 Copies`), desc(`C29 Copies`)) 
# data_fwrite <- synthetic_samples_C48
# name_data_frwite <- "synthetic_samples_C48"
# new_filename <- paste0(folder_data_unified, name_data_frwite, ".csv")
# print(new_filename)
# fwrite(data_fwrite, file = new_filename, append = FALSE, quote = "auto", showProgress = getOption("datatable.showProgress", interactive()), buffMB = 1024, nThread = 6)

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

# synthetic_samples_to_xtable <- synthetic_samples_to_xtable %>% select("Identifier","Median_blind", "Median_aware", "C1 Fraction",  "C2 Fraction",  "C14 Fraction", "C29 Fraction", "C48 Fraction", "C1 Copies","C2 Copies","C14 Copies","C29 Copies","C48 Copies", "Total Copies")
# colnames(synthetic_samples_to_xtable) <- c("Identifier","Median Depth blind", "Median Depth aware", "C1 Fraction",  "C2 Fraction",  "C14 Fraction", "C29 Fraction", "C48 Fraction", "C1 Copies","C2 Copies","C14 Copies","C29 Copies","C48 Copies", "Total Copies")
# 
# colnames(synthetic_samples_to_xtable)
# glimpse(synthetic_samples)
# file_table <- paste0(folder_tables,"mm_synthetic_samples", ".tex")
# print(file_table)
# sink(file_table)
# print(xtable(synthetic_samples_to_xtable), floating = FALSE, tabular.environment="longtable", include.rownames=FALSE)
# sink()

synthetic_samples_to_xtable_fraction <- synthetic_samples_to_xtable %>% select("Identifier","Median_blind", "Median_aware", "C1 Fraction",  "C2 Fraction",  "C14 Fraction", "C29 Fraction", "C48 Fraction")
colnames(synthetic_samples_to_xtable_fraction) <- c("Identifier","Median_Depth_blind", "Median_depth_aware", "C1 Fraction",  "C2 Fraction",  "C14 Fraction", "C29 Fraction", "C48 Fraction")
file_table <- paste0(folder_tables,"mm_synthetic_samples_fraction", ".tex")
print(file_table)
sink(file_table)
print(xtable(synthetic_samples_to_xtable_fraction), floating = FALSE, tabular.environment="longtable", include.rownames=FALSE)
sink()

synthetic_samples_to_xtable_copies <- synthetic_samples_to_xtable %>% select("Identifier","Median_blind", "Median_aware",   "C1 Copies",
                                                                             "C2 Copies",
                                                                             "C14 Copies",
                                                                             "C29 Copies",
                                                                             "C48 Copies",
                                                                             "Total Copies")
colnames(synthetic_samples_to_xtable_fraction) <- c("Identifier","Median_Depth_blind", "Median_depth_aware",   "C1 Copies",
                                                    "C2 Copies",
                                                    "C14 Copies",
                                                    "C29 Copies",
                                                    "C48 Copies",
                                                    "Total Copies")
file_table <- paste0(folder_tables,"mm_synthetic_samples_copies", ".tex")
print(file_table)
sink(file_table)
print(xtable(synthetic_samples_to_xtable_copies), floating = FALSE, tabular.environment="longtable", include.rownames=FALSE)
sink()
