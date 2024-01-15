# library(tidyverse)
# library(vroom)
# library(xtable)
# library(kableExtra)
# library(ggvenn)
# library(ggVennDiagram)
# library(VennDiagram)
# library(gplots)
# library(RColorBrewer)
# library(scales)
# library(ggthemes)
# #library(Cairo)
# 
# 
# rm(list = ls(all.names = TRUE)) #will clear all objects includes hidden objects.
# gc() #free up memrory and report the memory usage.

# run this first   /mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/create_figures_plots/wrapper_thesis.r

folder_plots <- "/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/GILOT_THESIS/folder_plots/"
folder_tables <- "/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/GILOT_THESIS/folder_tables/"

all_uniq_variants <- readRDS(file = paste0("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/", "all_uniq_variants", ".rds"))
all_the_variants <- readRDS(file = paste0("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/", "all_the_variants", ".rds")) %>% filter(SIGNATURE_VARIANT %in% all_uniq_variants)
unique_spike_positions <- (all_the_variants %>% filter(SIGNATURE_POS >= 21563) %>% filter(SIGNATURE_POS <= 25384))$SIGNATURE_POS %>% unique() %>% sort()
unique_pangenome_positions <- (all_the_variants)$SIGNATURE_POS %>% unique() %>% sort()
unique_nonspike_positions <- setdiff(unique_pangenome_positions, unique_spike_positions)


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

sample_depth_statistics_from_mpileups_gini <- vroom("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/counts_from_mpileups/sample_depth_statistics_from_mpileups_gini.csv")
data_depth <- sample_depth_statistics_from_mpileups_gini %>% select(STAT_IDENTIFIER, STAT_SAMPLE_METHOD, STAT_DEPTH_MEDIAN, STAT_DEPTH_MEAN, STAT_DEPTH_MAX) %>% unique()
names(data_depth)  <- c("SAMPLE_NAME", "METHOD", "DEPTH_MED", "DEPTH_MEA", "DEPTH_MAX")

background_and_consensus_for_marascuillo <- vroom("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/background_and_consensus_for_marascuillo_procedure_labelled.csv")
data_cts <- background_and_consensus_for_marascuillo %>% select(background_identifier, background_blaster, background_cohort, THESIS_SAMPLE_CT_S) %>% unique()
names(data_cts)  <- c("SAMPLE_NAME", "METHOD", "COHORT", "CT_S")

label_data <- merge(data_depth, data_cts)
# sample_depth_statistics_from_mpileups_gini %>% filter(STAT_DEPTH_MEDIAN > 250)

# filtering bad runs ------------------------------------------------------


sample_depth_statistics_from_mpileups_gini %>% ggplot(aes(x=STAT_DEPTH_MEDIAN)) + geom_density()

med_1 <- (sample_depth_statistics_from_mpileups_gini$STAT_DEPTH_MEDIAN %>% summary())[2]
min_1 <- (sample_depth_statistics_from_mpileups_gini$STAT_DEPTH_MIN %>% summary())[2]
mea_1 <- (sample_depth_statistics_from_mpileups_gini$STAT_DEPTH_MEAN %>% summary())[2]
max_1 <- (sample_depth_statistics_from_mpileups_gini$STAT_DEPTH_MAX %>% summary())[2]

worst_eliminated <- sample_depth_statistics_from_mpileups_gini # %>% filter(STAT_DEPTH_MIN >= min_1) %>% filter(STAT_DEPTH_MEDIAN >= med_1) %>% filter(STAT_DEPTH_MEAN >= mea_1) %>% filter(STAT_DEPTH_MAX >= max_1)

worst_eliminated$IDENTIFIER <- as.character(map(strsplit(worst_eliminated$STAT_IDENTIFIER, split = "_"), 1))

worst_eliminated$STAT_IDENTIFIER

files_blasted_control <- (background_and_consensus_for_marascuillo %>% filter(background_identifier %in% worst_eliminated$STAT_IDENTIFIER) %>% filter(background_blaster == "blasted") %>% filter(background_cohort == "gilot_all_controls"))$background_filename
files_blasted_experimental <- (background_and_consensus_for_marascuillo %>% filter(background_identifier %in% worst_eliminated$STAT_IDENTIFIER) %>% filter(background_blaster == "blasted") %>% filter(background_cohort == "gilot_experimental"))$background_filename

files_unblasted_control <- (background_and_consensus_for_marascuillo %>% filter(background_identifier %in% worst_eliminated$STAT_IDENTIFIER) %>% filter(background_blaster == "unblasted") %>% filter(background_cohort == "gilot_all_controls"))$background_filename
files_unblasted_experimental <- (background_and_consensus_for_marascuillo %>% filter(background_identifier %in% worst_eliminated$STAT_IDENTIFIER) %>% filter(background_blaster == "unblasted") %>% filter(background_cohort == "gilot_experimental"))$background_filename

counts_blasted_control <- vroom(files_blasted_control)
counts_blasted_control$POS <- counts_blasted_control$POS %>% as.integer()
counts_blasted_control$DEPTH <- counts_blasted_control$DEPTH %>% as.integer()
counts_blasted_control %>% select(POS, DEPTH)

counts_depth_plots <- vroom(c(files_blasted_control, files_blasted_experimental, files_unblasted_control, files_unblasted_experimental))
counts_depth_plots$POS <- counts_depth_plots$POS %>% as.integer()
counts_depth_plots$DEPTH <- counts_depth_plots$DEPTH %>% as.integer()
glimpse(counts_depth_plots)

# SKIP _ 1 ----------------------------------------------------------------


# print("experimental_depths")
# plot_OBJECT <- counts_depth_plots  %>% filter(SAMPLE_GROUP == "gilot_experimental") %>%   ggplot(aes(x = POS,
#                                                                                       y = DEPTH)) +
#   geom_point(alpha = 1 / 20, size = 0.05) +
#   geom_smooth(
#     method = "auto",
#     se = TRUE,
#     fullrange = FALSE,
#     level = 0.95,
#     linewidth = 0.5,
#     alpha=0.5
#   ) + 
#   facet_wrap(. ~ SAMPLE_METHOD, scales = "free_y") +
#   scale_y_continuous(labels = unit_format(unit = "k", scale = 1e-3), expand = c(0, 0)) +
#   scale_x_continuous(breaks = seq(0, 30000, by = 5000), labels = label_comma(), expand = c(0, 0)) +
#   geom_vline(xintercept = c(21563, 25384), linetype = "solid", color = "black", linewidth = 0.1) +
#   geom_vline(xintercept = unique_spike_positions, linetype = "dotted", color = "red", linewidth = 0.1, alpha=0.4) +
#   geom_vline(xintercept = unique_nonspike_positions, linetype = "dotted", color = "green", linewidth = 0.25) +
#   geom_hline(yintercept = 5000, linetype = "dotted", color = "red", linewidth = 0.25) +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
# 
# file <-
#   paste0(
#     "/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/GILOT_THESIS/figures/",
#     "depth_gilot_experimental_pangenome.png"
#   )
# ggsave(file,
#        device = "png",
#        dpi = 1800,
#        type = "cairo")
# 
# plot_OBJECT


# SKIP _ 2 ----------------------------------------------------------------



# print("gilot_all_controls")
# plot_OBJECT <- counts_depth_plots  %>% filter(SAMPLE_GROUP == "gilot_all_controls") %>%   ggplot(aes(x = POS,
#                                                                                       y = DEPTH)) +
#   geom_point(alpha = 1 / 20, size = 0.05) +
#   geom_smooth(
#     method = "auto",
#     se = TRUE,
#     fullrange = FALSE,
#     level = 0.95,
#     linewidth = 0.5,
#     alpha=0.5
#   ) + 
#   facet_wrap(. ~ SAMPLE_METHOD, scales = "free_y") +
#   scale_y_continuous(labels = unit_format(unit = "k", scale = 1e-3), expand = c(0, 0)) +
#   scale_x_continuous(breaks = seq(0, 30000, by = 5000), labels = label_comma(), expand = c(0, 0)) +
#   geom_vline(xintercept = c(21563, 25384), linetype = "solid", color = "black", linewidth = 0.1) +
#   geom_vline(xintercept = unique_spike_positions, linetype = "dotted", color = "red", linewidth = 0.1, alpha=0.4) +
#   geom_vline(xintercept = unique_nonspike_positions, linetype = "dotted", color = "green", linewidth = 0.25) +
#   geom_hline(yintercept = 5000, linetype = "dotted", color = "red", linewidth = 0.25) +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
# 
# file <-
#   paste0(
#     "/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/GILOT_THESIS/figures/",
#     "depth_gilot_all_controls_pangenome.png"
#   )
# ggsave(file,
#        device = "png",
#        dpi = 1800,
#        type = "cairo")
# plot_OBJECT


# SKIP _ 3 ----------------------------------------------------------------



# print("gilot_spike")
# plot_OBJECT <- counts_depth_plots %>% filter(POS >= 21563) %>% filter(POS <= 25384) %>%   ggplot(aes(x = POS,
#                                                                                       y = DEPTH)) +
#   geom_point(alpha = 1 / 20, size = 0.05) +
#   geom_smooth(
#     method = "auto",
#     se = TRUE,
#     fullrange = FALSE,
#     level = 0.95,
#     linewidth = 0.5,
#     alpha=0.5
#   ) + 
#   facet_wrap(SAMPLE_GROUP ~ SAMPLE_METHOD, scales = "free_y") +
#   scale_y_continuous(labels = unit_format(unit = "k", scale = 1e-3), expand = c(0, 0)) +
#   scale_x_continuous(breaks = seq(21500, 25500, by = 500), labels = label_comma(), expand = c(0, 0)) +
#   geom_vline(xintercept = c(21563, 25384), linetype = "solid", color = "black", linewidth = 0.1) +
#   geom_vline(xintercept = unique_spike_positions, linetype = "dotted", color = "red", linewidth = 0.2, alpha=0.6) +
#   # geom_vline(xintercept = unique_nonspike_positions, linetype = "dotted", color = "green", linewidth = 0.25) +
#   geom_hline(yintercept = 5000, linetype = "dotted", color = "red", linewidth = 0.25) +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
# 
# file <-
#   paste0(
#     "/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/GILOT_THESIS/figures/",
#     "depth_spike.png"
#   )
# ggsave(file,
#        device = "png",
#        dpi = 1800,
#        type = "cairo")
# 
# plot_OBJECT



# DEPTH PLOT - UMI blind Mapping - Human Samples --------------------------
#data_cts_select <- data_cts %>% filter(METHOD == "unblasted") %>% filter(SAMPLE_NAME %in% (counts_depth_plots %>% filter(SAMPLE_METHOD == "unblasted") %>% filter(SAMPLE_GROUP == "gilot_experimental"))$SAMPLE_NAME)
methodi <- "blasted"
cohorti <- "gilot_experimental"
pos_begin <- 1 # 21563 #filter(POS >= 21563) %>% filter(POS <= 25384)
pos_end <- 29903 # 25384
label_x <- "Genomic Position"
label_y <- "Observed Depth"
label_caption <- "UMI Blind : Human Samples : CT Value in Blue"
data_label_select <-
  label_data %>% filter(METHOD == methodi)  %>% filter(COHORT == cohorti) %>% filter(SAMPLE_NAME %in% (
    counts_depth_plots %>% filter(SAMPLE_METHOD == methodi) %>% filter(SAMPLE_GROUP == cohorti)
  )$SAMPLE_NAME)

depth_plot <-
  counts_depth_plots  %>% filter(POS >= pos_begin) %>% filter(POS <= pos_end) %>% filter(SAMPLE_METHOD == methodi) %>% filter(SAMPLE_GROUP == cohorti)

how_many <- depth_plot$SAMPLE_NAME %>% unique() %>% length()

print(paste(methodi, ":", cohorti, ": n =", how_many))

data_label_select$SAMPLE_NAME  <- as.character(map(strsplit(data_label_select$SAMPLE_NAME, split = "_"), 1))
depth_plot$SAMPLE_NAME  <- as.character(map(strsplit(depth_plot$SAMPLE_NAME, split = "_"), 1))

plot_object <- depth_plot %>% ggplot(aes(x = POS,
                                         y = DEPTH)) +
  geom_point(alpha = 1 / 20, size = 0.05) +
  geom_smooth(
    method = "auto",
    se = TRUE,
    fullrange = FALSE,
    level = 0.95,
    linewidth = 0.5,
    alpha=0.5,
    color = "yellow"
  ) + 
  facet_wrap(. ~ SAMPLE_NAME, scales = "free_y")  + 
  scale_y_continuous(labels = unit_format(unit = "k", scale = 1e-3), expand = c(0, 0)) +
  scale_x_continuous(breaks = seq(0, 30000, by = 5000), labels = label_comma(), expand = c(0, 0)) +
  geom_vline(xintercept = c(21563, 25384), linetype = "dotted", color = "blue", linewidth = 0.1 ) + 
  geom_hline(yintercept = 5000,    linetype = "dotted",    color = "red",    linewidth = 0.25  ) +   
  theme(axis.text.x = element_text(    angle = 45,    vjust = 1,    hjust = 1  ))  + 
  geom_text(    data    = data_label_select,    mapping = aes(x = -Inf, y = -Inf, label = CT_S),     hjust   = -0.1,     vjust   = -1 , color = "blue"  ) +
  labs(x = label_x) +
  labs(y = label_y) +
  labs(caption = label_caption)

file <-
  paste0(
    folder_plots,
    "depth_facet_wrap_",
    cohorti,
    "_",
    methodi,
    "_pangenome_identifier.png"
  )

print(file)

ggsave(file,
       device = "png",
       dpi = 1800,
       type = "cairo") # , dpi = "retina",  type = "cairo") # width = 210, height = 148.5 , units = c("mm")

plot_object

latex_table <- data_label_select %>% select(SAMPLE_NAME, DEPTH_MED, CT_S) %>% arrange(CT_S)
latex_table$DEPTH_MED <- latex_table$DEPTH_MED %>% as.integer()
names(latex_table) <- c("SAMPLE", "Median Depth Blind", "CT(S)")
print(xtable(latex_table), include.rownames=FALSE)
file_table <- paste0(folder_tables,file_path_sans_ext(basename(file)), ".tex")
print(file_table)
sink(file_table)
print(xtable(latex_table), include.rownames=FALSE)
sink()


# unblasted human ---------------------------------------------------------
#data_cts_select <- data_cts %>% filter(METHOD == "unblasted") %>% filter(SAMPLE_NAME %in% (counts_depth_plots %>% filter(SAMPLE_METHOD == "unblasted") %>% filter(SAMPLE_GROUP == "gilot_experimental"))$SAMPLE_NAME)
methodi <- "unblasted"
cohorti <- "gilot_experimental"
pos_begin <- 1 # 21563 #filter(POS >= 21563) %>% filter(POS <= 25384)
pos_end <- 29903 # 25384
label_x <- "Genomic Position"
label_y <- "Observed Depth"
label_caption <- "UMI Aware : Human Samples : CT Value in Blue"
data_label_select <-
  label_data %>% filter(METHOD == methodi)  %>% filter(COHORT == cohorti) %>% filter(SAMPLE_NAME %in% (
    counts_depth_plots %>% filter(SAMPLE_METHOD == methodi) %>% filter(SAMPLE_GROUP == cohorti)
  )$SAMPLE_NAME)

depth_plot <-
  counts_depth_plots  %>% filter(POS >= pos_begin) %>% filter(POS <= pos_end) %>% filter(SAMPLE_METHOD == methodi) %>% filter(SAMPLE_GROUP == cohorti)

how_many <- depth_plot$SAMPLE_NAME %>% unique() %>% length()

print(paste(methodi, ":", cohorti, ": n =", how_many))

data_label_select$SAMPLE_NAME  <- as.character(map(strsplit(data_label_select$SAMPLE_NAME, split = "_"), 1))
depth_plot$SAMPLE_NAME  <- as.character(map(strsplit(depth_plot$SAMPLE_NAME, split = "_"), 1))

plot_object <- depth_plot %>% ggplot(aes(x = POS,
                                         y = DEPTH)) +
  geom_point(alpha = 1 / 20, size = 0.05) +
  geom_smooth(
    method = "auto",
    se = TRUE,
    fullrange = FALSE,
    level = 0.95,
    linewidth = 0.5,
    alpha=0.5,
    color = "yellow"
  ) + 
  facet_wrap(. ~ SAMPLE_NAME, scales = "free_y")  + 
  scale_y_continuous(labels = unit_format(unit = "k", scale = 1e-3), expand = c(0, 0)) +
  scale_x_continuous(breaks = seq(0, 30000, by = 5000), labels = label_comma(), expand = c(0, 0)) +
  geom_vline(xintercept = c(21563, 25384), linetype = "dotted", color = "blue", linewidth = 0.1 ) + 
  geom_hline(yintercept = 5000,    linetype = "dotted",    color = "red",    linewidth = 0.25  ) +   
  theme(axis.text.x = element_text(    angle = 45,    vjust = 1,    hjust = 1  ))  + 
  geom_text(    data    = data_label_select,    mapping = aes(x = -Inf, y = -Inf, label = CT_S),     hjust   = -0.1,     vjust   = -1 , color = "blue"  ) +
  labs(x = label_x) +
  labs(y = label_y) +
  labs(caption = label_caption)

file <-
  paste0(
    folder_plots,
    "depth_facet_wrap_",
    cohorti,
    "_",
    methodi,
    "_pangenome_identifier.png"
  )

print(file)

ggsave(file,
       device = "png",
       dpi = 1800,
       type = "cairo") # , dpi = "retina",  type = "cairo") # width = 210, height = 148.5 , units = c("mm")

plot_object

latex_table <- data_label_select %>% select(SAMPLE_NAME, DEPTH_MED, CT_S) %>% arrange(CT_S)
latex_table$DEPTH_MED <- latex_table$DEPTH_MED %>% as.integer()
names(latex_table) <- c("SAMPLE", "Median Depth Aware", "CT(S)")
print(xtable(latex_table), include.rownames=FALSE)
file_table <- paste0(folder_tables,file_path_sans_ext(basename(file)), ".tex")
print(file_table)
sink(file_table)
print(xtable(latex_table), include.rownames=FALSE)
sink()
print(file)
print(file_table)

# unblasted synthetic UMI Aware ---------------------------------------------------------
#data_cts_select <- data_cts %>% filter(METHOD == "unblasted") %>% filter(SAMPLE_NAME %in% (counts_depth_plots %>% filter(SAMPLE_METHOD == "unblasted") %>% filter(SAMPLE_GROUP == "gilot_experimental"))$SAMPLE_NAME)
methodi <- "unblasted"
cohorti <- "gilot_all_controls"
pos_begin <- 1 # 21563 #filter(POS >= 21563) %>% filter(POS <= 25384)
pos_end <- 29903 # 25384
subplots_number <- 16
label_x <- "Genomic Position"
label_y <- "Observed Depth"
label_caption <- "UMI Aware : Synthetic Samples : Median Depth in Blue"
data_label_select <-
  label_data %>% filter(METHOD == methodi)  %>% filter(COHORT == cohorti) %>% filter(SAMPLE_NAME %in% (
    counts_depth_plots %>% filter(SAMPLE_METHOD == methodi) %>% filter(SAMPLE_GROUP == cohorti)
  )$SAMPLE_NAME)
data_label_select_latex <- data_label_select

data_label_select$sort_mean <- rowMeans(cbind(data_label_select$DEPTH_MED, data_label_select$DEPTH_MEA, data_label_select$DEPTH_MAX)) %>% as.numeric() %>% round(2)

data_label_select <- data_label_select %>% arrange(desc(sort_mean))

carry_over <- data_label_select[1:subplots_number,]$SAMPLE_NAME

depth_plot <-
  counts_depth_plots  %>% filter(POS >= pos_begin) %>% filter(POS <= pos_end) %>% filter(SAMPLE_METHOD == methodi) %>% filter(SAMPLE_GROUP == cohorti) %>% filter(SAMPLE_NAME %in% carry_over)

data_label_select <- data_label_select %>% filter(SAMPLE_NAME %in% carry_over) 

how_many <- depth_plot$SAMPLE_NAME %>% unique() %>% length()

print(paste(methodi, ":", cohorti, ": n =", how_many))

data_label_select$SAMPLE_NAME  <- as.character(map(strsplit(data_label_select$SAMPLE_NAME, split = "_"), 1))
depth_plot$SAMPLE_NAME  <- as.character(map(strsplit(depth_plot$SAMPLE_NAME, split = "_"), 1))

plot_object <- depth_plot %>% ggplot(aes(x = POS,
                                         y = DEPTH)) +
  geom_point(alpha = 1 / 20, size = 0.05) +
  geom_smooth(
    method = "auto",
    se = TRUE,
    fullrange = FALSE,
    level = 0.95,
    linewidth = 0.5,
    alpha=0.5,
    color = "yellow"
  ) + 
  facet_wrap(. ~ SAMPLE_NAME, scales = "free_y")  + 
  scale_y_continuous(labels = unit_format(unit = "k", scale = 1e-3), expand = c(0, 0)) +
  scale_x_continuous(breaks = seq(0, 30000, by = 5000), labels = label_comma(), expand = c(0, 0)) +
  geom_vline(xintercept = c(21563, 25384), linetype = "dotted", color = "blue", linewidth = 0.1 ) + 
  geom_hline(yintercept = 5000,    linetype = "dotted",    color = "red",    linewidth = 0.25  ) +   
  theme(axis.text.x = element_text(    angle = 45,    vjust = 1,    hjust = 1  ))  + 
  geom_text(    data    = data_label_select,    mapping = aes(x = -Inf, y = -Inf, label = DEPTH_MED),     hjust   = -0.1,     vjust   = -1 , color = "blue"  ) +
  labs(x = label_x) +
  labs(y = label_y) +
  labs(caption = label_caption)

file <-
  paste0(
    folder_plots,
    "depth_facet_wrap_",
    cohorti,
    "_",
    methodi,
    "_pangenome_identifier.png"
  )

print(file)

ggsave(file,
       device = "png",
       dpi = 3600,
       type = "cairo") # , dpi = "retina",  type = "cairo") # width = 210, height = 148.5 , units = c("mm")

plot_object

latex_table <- data_label_select_latex %>% select(SAMPLE_NAME, DEPTH_MED) %>% arrange(desc(DEPTH_MED))
latex_table$DEPTH_MED <- latex_table$DEPTH_MED %>% as.integer()
names(latex_table) <- c("SAMPLE", "Median Depth Aware")
print(xtable(latex_table), include.rownames=FALSE)
file_table <- paste0(folder_tables,file_path_sans_ext(basename(file)), ".tex")
print(file_table)
sink(file_table)
print(xtable(latex_table), floating = FALSE, tabular.environment="longtable", include.rownames=FALSE)
sink()
print(file)
print(file_table)



# blasted synthetic UMI blind ---------------------------------------------------------
#data_cts_select <- data_cts %>% filter(METHOD == "unblasted") %>% filter(SAMPLE_NAME %in% (counts_depth_plots %>% filter(SAMPLE_METHOD == "unblasted") %>% filter(SAMPLE_GROUP == "gilot_experimental"))$SAMPLE_NAME)
methodi <- "blasted"
cohorti <- "gilot_all_controls"
pos_begin <- 1 # 21563 #filter(POS >= 21563) %>% filter(POS <= 25384)
pos_end <- 29903 # 25384
subplots_number <- 16
label_x <- "Genomic Position"
label_y <- "Observed Depth"
label_caption <- "UMI Blind : Synthetic Samples : Median Depth in Blue"
data_label_select <-
  label_data %>% filter(METHOD == methodi)  %>% filter(COHORT == cohorti) %>% filter(SAMPLE_NAME %in% (
    counts_depth_plots %>% filter(SAMPLE_METHOD == methodi) %>% filter(SAMPLE_GROUP == cohorti)
  )$SAMPLE_NAME)
data_label_select_latex <- data_label_select

data_label_select$sort_mean <- rowMeans(cbind(data_label_select$DEPTH_MED, data_label_select$DEPTH_MEA, data_label_select$DEPTH_MAX)) %>% as.numeric() %>% round(2)

data_label_select <- data_label_select %>% arrange(desc(sort_mean))

# depth_plot <-
#   counts_depth_plots  %>% filter(POS >= pos_begin) %>% filter(POS <= pos_end) %>% filter(SAMPLE_METHOD == methodi) %>% filter(SAMPLE_GROUP == cohorti) %>% filter(SAMPLE_NAME %in% data_label_select[1:subplots_number,]$SAMPLE_NAME)
depth_plot <-
  counts_depth_plots  %>% filter(POS >= pos_begin) %>% filter(POS <= pos_end) %>% filter(SAMPLE_METHOD == methodi) %>% filter(SAMPLE_GROUP == cohorti) %>% filter(SAMPLE_NAME %in% carry_over)


#data_label_select <- data_label_select %>% filter(SAMPLE_NAME %in% data_label_select[1:subplots_number,]$SAMPLE_NAME)   # carry_over
data_label_select <- data_label_select %>% filter(SAMPLE_NAME %in% carry_over)

how_many <- depth_plot$SAMPLE_NAME %>% unique() %>% length()

print(paste(methodi, ":", cohorti, ": n =", how_many))

data_label_select$SAMPLE_NAME  <- as.character(map(strsplit(data_label_select$SAMPLE_NAME, split = "_"), 1))
depth_plot$SAMPLE_NAME  <- as.character(map(strsplit(depth_plot$SAMPLE_NAME, split = "_"), 1))

plot_object <- depth_plot %>% ggplot(aes(x = POS,
                                         y = DEPTH)) +
  geom_point(alpha = 1 / 20, size = 0.05) +
  geom_smooth(
    method = "auto",
    se = TRUE,
    fullrange = FALSE,
    level = 0.95,
    linewidth = 0.5,
    alpha=0.5,
    color = "yellow"
  ) + 
  facet_wrap(. ~ SAMPLE_NAME, scales = "free_y")  + 
  scale_y_continuous(labels = unit_format(unit = "k", scale = 1e-3), expand = c(0, 0)) +
  scale_x_continuous(breaks = seq(0, 30000, by = 5000), labels = label_comma(), expand = c(0, 0)) +
  geom_vline(xintercept = c(21563, 25384), linetype = "dotted", color = "blue", linewidth = 0.1 ) + 
  geom_hline(yintercept = 5000,    linetype = "dotted",    color = "red",    linewidth = 0.25  ) +   
  theme(axis.text.x = element_text(    angle = 45,    vjust = 1,    hjust = 1  ))  + 
  geom_text(    data    = data_label_select,    mapping = aes(x = -Inf, y = -Inf, label = DEPTH_MED),     hjust   = -0.1,     vjust   = -1 , color = "blue"  ) +
  labs(x = label_x) +
  labs(y = label_y) +
  labs(caption = label_caption)

file <-
  paste0(
    folder_plots,
    "depth_facet_wrap_",
    cohorti,
    "_",
    methodi,
    "_pangenome_identifier.png"
  )

print(file)

ggsave(file,
       device = "png",
       dpi = 3600,
       type = "cairo") # , dpi = "retina",  type = "cairo") # width = 210, height = 148.5 , units = c("mm")

plot_object

latex_table <- data_label_select_latex %>% select(SAMPLE_NAME, DEPTH_MED) %>% arrange(desc(DEPTH_MED))
latex_table$DEPTH_MED <- latex_table$DEPTH_MED %>% as.integer()
names(latex_table) <- c("SAMPLE", "Median Depth Blind")
print(xtable(latex_table), include.rownames=FALSE)
file_table <- paste0(folder_tables,file_path_sans_ext(basename(file)), ".tex")
print(file_table)
sink(file_table)
print(xtable(latex_table), floating = FALSE, tabular.environment="longtable", include.rownames=FALSE)
sink()
print(file)
print(file_table)



#########################################################################################################################################



# 
# # DEPTH PLOT: UMI aware Mapping - Human Samples ---------------------------
# 
# 
# #data_cts_select <- data_cts %>% filter(METHOD == "unblasted") %>% filter(SAMPLE_NAME %in% (counts_depth_plots %>% filter(SAMPLE_METHOD == "unblasted") %>% filter(SAMPLE_GROUP == "gilot_experimental"))$SAMPLE_NAME)
# methodi <- "unblasted"
# cohorti <- "gilot_experimental"
# pos_begin <- 1 # 21563 #filter(POS >= 21563) %>% filter(POS <= 25384)
# pos_end <- 29903 # 25384
# label_x <- "Genomic Position"
# label_y <- "Observed Depth"
# label_caption <- "UMI Aware:Human Samples : CT Value in Blue"
# data_label_select <-
#   label_data %>% filter(METHOD == methodi)  %>% filter(COHORT == cohorti) %>% filter(SAMPLE_NAME %in% (
#     counts_depth_plots %>% filter(SAMPLE_METHOD == methodi) %>% filter(SAMPLE_GROUP == cohorti)
#   )$SAMPLE_NAME)
# 
# depth_plot <-
#   counts_depth_plots  %>% filter(POS >= pos_begin) %>% filter(POS <= pos_end) %>% filter(SAMPLE_METHOD == methodi) %>% filter(SAMPLE_GROUP == cohorti)
# 
# how_many <- depth_plot$SAMPLE_NAME %>% unique() %>% length()
# 
# print(paste(methodi, ":", cohorti, ": n =", how_many))
# 
# data_label_select$SAMPLE_NAME  <- as.character(map(strsplit(data_label_select$SAMPLE_NAME, split = "_"), 1))
# depth_plot$SAMPLE_NAME  <- as.character(map(strsplit(depth_plot$SAMPLE_NAME, split = "_"), 1))
# 
# plot_object <- depth_plot %>% ggplot(aes(x = POS,
#                                          y = DEPTH)) +
#   geom_point(alpha = 1 / 20, size = 0.05) +
#   geom_smooth(
#     method = "auto",
#     se = TRUE,
#     fullrange = FALSE,
#     level = 0.95,
#     linewidth = 0.5,
#     alpha=0.5,
#     color = "yellow"
#   ) + 
#   facet_wrap(. ~ SAMPLE_NAME, scales = "free_y")  + 
#   scale_y_continuous(labels = unit_format(unit = "k", scale = 1e-3), expand = c(0, 0)) +
#   scale_x_continuous(breaks = seq(0, 30000, by = 5000), labels = label_comma(), expand = c(0, 0)) +
#   geom_vline(xintercept = c(21563, 25384), linetype = "dotted", color = "blue", linewidth = 0.1 ) + 
#   geom_hline(yintercept = 5000,    linetype = "dotted",    color = "red",    linewidth = 0.25  ) +   
#   theme(axis.text.x = element_text(    angle = 45,    vjust = 1,    hjust = 1  ))  + 
#   geom_text(    data    = data_label_select,    mapping = aes(x = -Inf, y = -Inf, label = CT_S),     hjust   = -0.1,     vjust   = -1 , color = "blue"  ) +
#   labs(x = label_x) +
#   labs(y = label_y) +
#   labs(caption = label_caption)
# 
# file <-
#   paste0(
#     "/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/GILOT_THESIS/figures/",
#     "depth_facet_wrap_",
#     cohorti,
#     "_",
#     methodi,
#     "_pangenome_identifier.png"
#   )
# print(file)
# 
# ggsave(file,
#        device = "png",
#        dpi = 1800,
#        type = "cairo") # , dpi = "retina",  type = "cairo") # width = 210, height = 148.5 , units = c("mm")
# 
# plot_object
# 
# latex_table <- data_label_select %>% select(SAMPLE_NAME, DEPTH_MED, CT_S)
# names(latex_table) <- c("SAMPLE", "Median Depth", "CT(S)")
# print(xtable(latex_table), include.rownames=FALSE)
# 
# 
# 
# 
# 
# 
# 
# #data_cts_select <- data_cts %>% filter(METHOD == "unblasted") %>% filter(SAMPLE_NAME %in% (counts_depth_plots %>% filter(SAMPLE_METHOD == "unblasted") %>% filter(SAMPLE_GROUP == "gilot_experimental"))$SAMPLE_NAME)
# methodi <- "unblasted"
# cohorti <- "gilot_experimental"
# pos_begin <- 1 # 21563 #filter(POS >= 21563) %>% filter(POS <= 25384)
# pos_end <- 29903 # 25384
# data_label_select <-
#   label_data %>% filter(METHOD == methodi)  %>% filter(COHORT == cohorti) %>% filter(SAMPLE_NAME %in% (
#     counts_depth_plots %>% filter(SAMPLE_METHOD == methodi) %>% filter(SAMPLE_GROUP == cohorti)
#   )$SAMPLE_NAME)
# 
# depth_plot <-
#   counts_depth_plots  %>% filter(POS >= pos_begin) %>% filter(POS <= pos_end) %>% filter(SAMPLE_METHOD == methodi) %>% filter(SAMPLE_GROUP == cohorti)
# 
# how_many <- depth_plot$SAMPLE_NAME %>% unique() %>% length()
# 
# print(paste(methodi, ":", cohorti, ": n =", how_many))
# 
# plot_object <- depth_plot %>% ggplot(aes(x = POS,
#                                          y = DEPTH)) +
#   geom_point(alpha = 1 / 20, size = 0.05) +
#   geom_smooth(
#     method = "auto",
#     se = TRUE,
#     fullrange = FALSE,
#     level = 0.95,
#     linewidth = 0.5,
#     alpha=0.5,
#     color = "yellow"
#   ) + 
#   facet_wrap(. ~ SAMPLE_NAME, scales = "free_y")  + 
#   scale_y_continuous(labels = unit_format(unit = "k", scale = 1e-3), expand = c(0, 0)) +
#   scale_x_continuous(breaks = seq(0, 30000, by = 5000), labels = label_comma(), expand = c(0, 0)) +
#   geom_vline(xintercept = c(21563, 25384), linetype = "dotted", color = "blue", linewidth = 0.1 ) + 
#   geom_hline(yintercept = 5000,    linetype = "dotted",    color = "red",    linewidth = 0.25  ) +   
#   theme(axis.text.x = element_text(    angle = 45,    vjust = 1,    hjust = 1  ))  + 
#   geom_text(    data    = data_label_select,    mapping = aes(x = -Inf, y = -Inf, label = paste("CT(S) value:", CT_S)),     hjust   = -0.1,     vjust   = -1 , color = "blue"  )
# 
# #geom_text(data = data_label_select,mapping = aes(x = -Inf, y = -Inf, label = paste("Median Depth:", format(DEPTH_MED, big.mark = ",", scientific = FALSE))), hjust = -0.1, vjust   = -1)
# 
# # + scale_x_continuous(labels = label_comma()) + scale_y_continuous(breaks = seq(0, 250000, by = 50000))
# #mapping = aes(x = -Inf, y = -Inf, label = paste("CT(S) value:", CT_S, "- Median Depth:",  DEPTH_MED)),
# # +  theme(axis.text.y = element_text(angle = 45, vjust = 1, hjust = 1))
# #    mapping = aes(x = -Inf, y = -Inf, label = paste("Median Depth:",  format(DEPTH_MED,big.mark=",",scientific=FALSE))), #mapping = aes(x = -Inf, y = -Inf, label = paste("CT(S) value:", CT_S, "- Median Depth:",  DEPTH_MED)),
# # hjust   = 1.05,
# # vjust   = 1.5
# 
# file <-
#   paste0(
#     "/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/GILOT_THESIS/figures/",
#     "depth_facet_wrap_",
#     cohorti,
#     "_",
#     methodi,
#     "_pangenome_identifier.png"
#   )
# ggsave(file,
#        device = "png",
#        dpi = 1200,
#        type = "cairo") # , dpi = "retina",  type = "cairo") # width = 210, height = 148.5 , units = c("mm")
# 
# 
# latex_table <- data_label_select %>% select(SAMPLE_NAME, DEPTH_MED, CT_S)
# latex_table$DEPTH_MED <- latex_table$DEPTH_MED %>% as.integer()
# names(latex_table) <- c("SAMPLE", "Median Depth", "CT(S)")
# print(xtable(latex_table), include.rownames=FALSE)
# 
# 
# 
# 
# 
# 
# 
# #data_cts_select <- data_cts %>% filter(METHOD == "unblasted") %>% filter(SAMPLE_NAME %in% (counts_depth_plots %>% filter(SAMPLE_METHOD == "unblasted") %>% filter(SAMPLE_GROUP == "gilot_experimental"))$SAMPLE_NAME)
# methodi <- "unblasted"
# cohorti <- "gilot_all_controls"
# pos_begin <- 1 # 21563 #filter(POS >= 21563) %>% filter(POS <= 25384)
# pos_end <- 29903 # 25384
# data_label_select <-
#   label_data %>% filter(METHOD == methodi)  %>% filter(COHORT == cohorti) %>% filter(SAMPLE_NAME %in% (
#     counts_depth_plots %>% filter(SAMPLE_METHOD == methodi) %>% filter(SAMPLE_GROUP == cohorti)
#   )$SAMPLE_NAME)
# 
# depth_plot <-
#   counts_depth_plots  %>% filter(POS >= pos_begin) %>% filter(POS <= pos_end) %>% filter(SAMPLE_METHOD == methodi) %>% filter(SAMPLE_GROUP == cohorti)
# 
# how_many <- depth_plot$SAMPLE_NAME %>% unique() %>% length()
# 
# print(paste(methodi, ":", cohorti, ": n =", how_many))
# 
# plot_object <- depth_plot %>% ggplot(aes(x = POS,
#                                          y = DEPTH)) +
#   geom_point(alpha = 1 / 20, size = 0.05) +
#   geom_smooth(
#     method = "auto",
#     se = TRUE,
#     fullrange = FALSE,
#     level = 0.95,
#     linewidth = 0.5,
#     alpha=0.5,
#     color = "yellow"
#   ) + 
#   facet_wrap(. ~ SAMPLE_NAME, scales = "free_y")  + 
#   scale_y_continuous(labels = unit_format(unit = "k", scale = 1e-3), expand = c(0, 0)) +
#   scale_x_continuous(breaks = seq(0, 30000, by = 5000), labels = label_comma(), expand = c(0, 0)) +
#   geom_vline(xintercept = c(21563, 25384), linetype = "dotted", color = "blue", linewidth = 0.1 ) + 
#   geom_hline(yintercept = 5000,    linetype = "dotted",    color = "red",    linewidth = 0.25  ) +   
#   theme(axis.text.x = element_text(    angle = 45,    vjust = 1,    hjust = 1  ))  + 
#   geom_text(    data    = data_label_select,    mapping = aes(x = -Inf, y = -Inf, label = format(DEPTH_MED,big.mark=",",scientific=FALSE)),     hjust   = -0.1,     vjust   = -1 , color = "blue"  )
# 
# #geom_text(data = data_label_select,mapping = aes(x = -Inf, y = -Inf, label = paste("Median Depth:", format(DEPTH_MED, big.mark = ",", scientific = FALSE))), hjust = -0.1, vjust   = -1)
# 
# # + scale_x_continuous(labels = label_comma()) + scale_y_continuous(breaks = seq(0, 250000, by = 50000))
# #mapping = aes(x = -Inf, y = -Inf, label = paste("CT(S) value:", CT_S, "- Median Depth:",  DEPTH_MED)),
# # +  theme(axis.text.y = element_text(angle = 45, vjust = 1, hjust = 1))
# #    mapping = aes(x = -Inf, y = -Inf, label = paste("Median Depth:",  format(DEPTH_MED,big.mark=",",scientific=FALSE))), #mapping = aes(x = -Inf, y = -Inf, label = paste("CT(S) value:", CT_S, "- Median Depth:",  DEPTH_MED)),
# # hjust   = 1.05,
# # vjust   = 1.5
# 
# file <-
#   paste0(
#     "/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/GILOT_THESIS/figures/",
#     "depth_facet_wrap_",
#     cohorti,
#     "_",
#     methodi,
#     "_pangenome_identifier.png"
#   )
# ggsave(file,
#        device = "png",
#        dpi = 1600,
#        type = "cairo") # , dpi = "retina",  type = "cairo") # width = 210, height = 148.5 , units = c("mm")
# 
# 
# latex_table <- data_label_select %>% select(SAMPLE_NAME, DEPTH_MED, CT_S)
# latex_table$DEPTH_MED <- latex_table$DEPTH_MED %>% as.integer()
# names(latex_table) <- c("SAMPLE", "Median Depth", "CT(S)")
# print(xtable(latex_table %>% select(SAMPLE, "Median Depth")), include.rownames=FALSE)
# 
# 
# 
# #data_cts_select <- data_cts %>% filter(METHOD == "unblasted") %>% filter(SAMPLE_NAME %in% (counts_depth_plots %>% filter(SAMPLE_METHOD == "unblasted") %>% filter(SAMPLE_GROUP == "gilot_experimental"))$SAMPLE_NAME)
# methodi <- "blasted"
# cohorti <- "gilot_all_controls"
# pos_begin <- 1 # 21563 #filter(POS >= 21563) %>% filter(POS <= 25384)
# pos_end <- 29903 # 25384
# data_label_select <-
#   label_data %>% filter(METHOD == methodi)  %>% filter(COHORT == cohorti) %>% filter(SAMPLE_NAME %in% (
#     counts_depth_plots %>% filter(SAMPLE_METHOD == methodi) %>% filter(SAMPLE_GROUP == cohorti)
#   )$SAMPLE_NAME)
# 
# depth_plot <-
#   counts_depth_plots  %>% filter(POS >= pos_begin) %>% filter(POS <= pos_end) %>% filter(SAMPLE_METHOD == methodi) %>% filter(SAMPLE_GROUP == cohorti)
# 
# how_many <- depth_plot$SAMPLE_NAME %>% unique() %>% length()
# 
# print(paste(methodi, ":", cohorti, ": n =", how_many))
# 
# plot_object <- depth_plot %>% ggplot(aes(x = POS,
#                                          y = DEPTH)) +
#   geom_point(alpha = 1 / 20, size = 0.05) +
#   geom_smooth(
#     method = "auto",
#     se = TRUE,
#     fullrange = FALSE,
#     level = 0.95,
#     linewidth = 0.5,
#     alpha=0.5,
#     color = "yellow"
#   ) + 
#   facet_wrap(. ~ SAMPLE_NAME, scales = "free_y")  + 
#   scale_y_continuous(labels = unit_format(unit = "k", scale = 1e-3), expand = c(0, 0)) +
#   scale_x_continuous(breaks = seq(0, 30000, by = 5000), labels = label_comma(), expand = c(0, 0)) +
#   geom_vline(xintercept = c(21563, 25384), linetype = "dotted", color = "blue", linewidth = 0.1 ) + 
#   geom_hline(yintercept = 5000,    linetype = "dotted",    color = "red",    linewidth = 0.25  ) +   
#   theme(axis.text.x = element_text(    angle = 45,    vjust = 1,    hjust = 1  ))  + 
#   geom_text(    data    = data_label_select,    mapping = aes(x = -Inf, y = -Inf, label = format(DEPTH_MED,big.mark=",",scientific=FALSE)),     hjust   = -0.1,     vjust   = -1 , color = "blue"  )
# 
# #geom_text(data = data_label_select,mapping = aes(x = -Inf, y = -Inf, label = paste("Median Depth:", format(DEPTH_MED, big.mark = ",", scientific = FALSE))), hjust = -0.1, vjust   = -1)
# 
# # + scale_x_continuous(labels = label_comma()) + scale_y_continuous(breaks = seq(0, 250000, by = 50000))
# #mapping = aes(x = -Inf, y = -Inf, label = paste("CT(S) value:", CT_S, "- Median Depth:",  DEPTH_MED)),
# # +  theme(axis.text.y = element_text(angle = 45, vjust = 1, hjust = 1))
# #    mapping = aes(x = -Inf, y = -Inf, label = paste("Median Depth:",  format(DEPTH_MED,big.mark=",",scientific=FALSE))), #mapping = aes(x = -Inf, y = -Inf, label = paste("CT(S) value:", CT_S, "- Median Depth:",  DEPTH_MED)),
# # hjust   = 1.05,
# # vjust   = 1.5
# 
# file <-
#   paste0(
#     "/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/GILOT_THESIS/figures/",
#     "depth_facet_wrap_",
#     cohorti,
#     "_",
#     methodi,
#     "_pangenome_identifier.png"
#   )
# ggsave(file,
#        device = "png",
#        dpi = 1600,
#        type = "cairo") # , dpi = "retina",  type = "cairo") # width = 210, height = 148.5 , units = c("mm")
# 
# 
# latex_table <- data_label_select %>% select(SAMPLE_NAME, DEPTH_MED, CT_S)
# latex_table$DEPTH_MED <- latex_table$DEPTH_MED %>% as.integer()
# names(latex_table) <- c("SAMPLE", "Median Depth", "CT(S)")
# print(xtable(latex_table %>% select(SAMPLE, "Median Depth")), include.rownames=FALSE)
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# ###############################################################################
# counts_depth_plots %>% filter(POS >= 21563) %>% filter(POS <= 25384) %>% filter(SAMPLE_GROUP == "gilot_experimental") %>%   ggplot(
#   aes(
#     x = POS,
#     y = DEPTH
#   )
# ) +
#   geom_point(alpha = 1/20, size = 0.05) +
#   geom_smooth(
#     method = "auto",
#     se = TRUE,
#     fullrange = FALSE,
#     level = 0.95
#   ) + facet_wrap(. ~ SAMPLE_METHOD, scales="free_y") + scale_y_continuous(labels = label_comma())  + scale_x_continuous(breaks = seq(21500, 25500, by = 1000))  + geom_vline(xintercept = c(21563,25384), linetype="dotted", color = "blue", size=0.25)
# # + scale_y_continuous(labels = label_comma())
# file <- paste0("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/GILOT_THESIS/figures/","depth_gilot_experimental_spike_protein.png")
# ggsave(file, device = "png", dpi = "retina", type = "cairo")
# 
# 
# 
# 
# 
# 
# 
# 
# counts_depth_plots  %>%   ggplot(
#   aes(
#     x = POS,
#     y = DEPTH # ,
#     # color = SAMPLE_GROUP # ,
#     #shape = SAMPLE_METHOD
#   )
# ) +
#   geom_point(alpha = 1/20, size = 0.05) +
#   geom_smooth(
#     method = "auto",
#     se = TRUE,
#     fullrange = FALSE,
#     level = 0.95
#   ) + scale_y_continuous(labels = label_comma())  + scale_x_continuous(breaks = seq(0, 30000, by = 5000), labels = label_comma())  + geom_vline(xintercept = c(21563,25384), linetype="dotted", color = "blue", size=0.25)  + facet_wrap(SAMPLE_GROUP ~ SAMPLE_METHOD, scales="free_y") +
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +  theme(axis.text.y = element_text(angle = 45, vjust = 1, hjust = 1))
# 
# # + facet_wrap(. ~ SAMPLE_METHOD, scales="free_y") + scale_y_continuous(labels = label_comma())  + scale_x_continuous(breaks = seq(21500, 25500, by = 1000))  + geom_vline(xintercept = c(21563,25384), linetype="dotted", color = "blue", size=0.25)
# 
# file <- paste0("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/GILOT_THESIS/figures/","depth_all_pangenome.png")
# ggsave(file, device = "png", dpi = "retina", type = "cairo")
# 
# counts_depth_plots %>% str()
# 
# counts_depth_plots %>% filter(SAMPLE_METHOD == "unblasted") %>% filter(SAMPLE_GROUP == "gilot_experimental") %>%   ggplot(
#   aes(
#     x = POS,
#     y = DEPTH
#   )
# ) +
#   geom_point(alpha = 1/20, size = 0.05) +
#   geom_smooth(
#     method = "auto",
#     se = TRUE,
#     fullrange = FALSE,
#     level = 0.95
#   ) + facet_wrap(. ~ SAMPLE_NAME, scales="free_y")  + scale_y_continuous(labels = label_comma()) + scale_x_continuous(labels = label_comma()) + geom_vline(xintercept = c(21563,25384), linetype="dotted", color = "blue", size=0.25) +   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +  theme(axis.text.y = element_text(angle = 45, vjust = 1, hjust = 1))
# # + scale_x_continuous(labels = label_comma()) + scale_y_continuous(breaks = seq(0, 250000, by = 50000))
# 
# file <- paste0("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/GILOT_THESIS/figures/","depth_gilot_experimental_unblasted_pangenome_identifier.png")
# ggsave(file, device = "png", dpi = "retina", type = "cairo")
# # 
# # data_cts_select <- data_cts %>% filter(SAMPLE_NAME %in% (counts_depth_plots %>% filter(SAMPLE_METHOD == "blasted") %>% filter(SAMPLE_GROUP == "gilot_experimental"))$SAMPLE_NAME)
# # 
# # counts_depth_plots %>% filter(SAMPLE_METHOD == "blasted") %>% filter(SAMPLE_GROUP == "gilot_experimental") %>%   ggplot(
# #   aes(
# #     x = POS,
# #     y = DEPTH
# #   )
# # ) +
# #   geom_point(alpha = 1/20, size = 0.05) +
# #   geom_smooth(
# #     method = "auto",
# #     se = TRUE,
# #     fullrange = FALSE,
# #     level = 0.95
# #   ) + facet_wrap(. ~ SAMPLE_NAME, scales="free_y")  + scale_y_continuous(labels = label_comma()) + scale_x_continuous(labels = label_comma()) + geom_vline(xintercept = c(21563,25384), linetype="dotted", color = "blue", size=0.25) +   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +  theme(axis.text.y = element_text(angle = 45, vjust = 1, hjust = 1)) + geom_text(
# #     data    = data_cts_select,
# #     mapping = aes(x = -Inf, y = -Inf, label = paste("ct value:", CT_S)),
# #     hjust   = -0.1,
# #     vjust   = -1
# #   )
# # # + scale_x_continuous(labels = label_comma()) + scale_y_continuous(breaks = seq(0, 250000, by = 50000))
# # 
# # 
# # file <- paste0("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/GILOT_THESIS/figures/","depth_gilot_experimental_blasted_pangenome_identifier_ctvalue.png")
# # ggsave(file, device = "png", dpi = "retina", type = "cairo")
# # 
# 
# 
# 
# #data_cts_select <- data_cts %>% filter(METHOD == "unblasted") %>% filter(SAMPLE_NAME %in% (counts_depth_plots %>% filter(SAMPLE_METHOD == "unblasted") %>% filter(SAMPLE_GROUP == "gilot_experimental"))$SAMPLE_NAME)
# data_label_select <- label_data %>% filter(METHOD == "unblasted")  %>% filter(COHORT == "gilot_experimental") %>% filter(SAMPLE_NAME %in% (counts_depth_plots %>% filter(SAMPLE_METHOD == "unblasted") %>% filter(SAMPLE_GROUP == "gilot_experimental"))$SAMPLE_NAME)
# counts_depth_plots %>% filter(SAMPLE_METHOD == "unblasted") %>% filter(SAMPLE_GROUP == "gilot_experimental") %>%   ggplot(
#   aes(
#     x = POS,
#     y = DEPTH
#   )
# ) +
#   geom_point(alpha = 1/20, size = 0.05) +
#   geom_smooth(
#     method = "auto",
#     se = TRUE,
#     fullrange = FALSE,
#     level = 0.95
#   ) + facet_wrap(. ~ SAMPLE_NAME, scales="free_y")  + scale_y_continuous(labels = label_comma()) + scale_x_continuous(labels = label_comma()) + geom_vline(xintercept = c(21563,25384), linetype="dotted", color = "blue", size=0.1) + geom_hline(yintercept = 5000, linetype="dotted", color = "red", size=0.25) +   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +  theme(axis.text.y = element_text(angle = 45, vjust = 1, hjust = 1)) + geom_text(
#     data    = data_label_select,
#     mapping = aes(x = -Inf, y = -Inf, label = paste("CT(S) value:", CT_S, "- Median Depth:",  format(DEPTH_MED,big.mark=",",scientific=FALSE))), #mapping = aes(x = -Inf, y = -Inf, label = paste("CT(S) value:", CT_S, "- Median Depth:",  DEPTH_MED)),
#     # hjust   = 1.05,
#     # vjust   = 1.5
#     hjust   = -0.1,
#     vjust   = -1
#   )
# # + scale_x_continuous(labels = label_comma()) + scale_y_continuous(breaks = seq(0, 250000, by = 50000))
# 
# file <- paste0("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/GILOT_THESIS/figures/","depth_gilot_experimental_unblasted_pangenome_identifier_ctvalue_depth.png")
# ggsave(file, device = "png", dpi = "retina",   width = 210, height = 148.5 , units = c("mm"), type = "cairo")
# 
# 
# #data_cts_select <- data_cts %>% filter(METHOD == "unblasted") %>% filter(SAMPLE_NAME %in% (counts_depth_plots %>% filter(SAMPLE_METHOD == "unblasted") %>% filter(SAMPLE_GROUP == "gilot_experimental"))$SAMPLE_NAME)
# data_label_select <- label_data %>% filter(METHOD == "blasted")  %>% filter(COHORT == "gilot_experimental") %>% filter(SAMPLE_NAME %in% (counts_depth_plots %>% filter(SAMPLE_METHOD == "blasted") %>% filter(SAMPLE_GROUP == "gilot_experimental"))$SAMPLE_NAME)
# counts_depth_plots %>% filter(SAMPLE_METHOD == "blasted") %>% filter(SAMPLE_GROUP == "gilot_experimental") %>% ggplot(
#   aes(
#     x = POS,
#     y = DEPTH
#   )
# ) +
#   geom_point(alpha = 1/20, size = 0.05) +
#   geom_smooth(
#     method = "auto",
#     se = TRUE,
#     fullrange = FALSE,
#     level = 0.95
#   ) + facet_wrap(. ~ SAMPLE_NAME, scales="free_y")  + scale_y_continuous(labels = label_comma()) + scale_x_continuous(labels = label_comma()) + geom_vline(xintercept = c(21563,25384), linetype="dotted", color = "blue", size=0.1) + geom_hline(yintercept = 5000, linetype="dotted", color = "red", size=0.25) +   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +  theme(axis.text.y = element_text(angle = 45, vjust = 1, hjust = 1)) + 
#   geom_text(    data    = data_label_select,    mapping = aes(x = -Inf, y = -Inf, label = paste("CT(S) value:", CT_S, "- Median Depth:",  format(DEPTH_MED,big.mark=",",scientific=FALSE))),     hjust   = -0.1,     vjust   = -1   )
# # + scale_x_continuous(labels = label_comma()) + scale_y_continuous(breaks = seq(0, 250000, by = 50000))
# 
# file <- paste0("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/GILOT_THESIS/figures/","depth_gilot_experimental_blasted_pangenome_identifier_ctvalue_depth.png")
# ggsave(file, device = "png", dpi = "retina",   width = 210, height = 148.5 , units = c("mm"), type = "cairo")
# 
# 
# #data_cts_select <- data_cts %>% filter(METHOD == "unblasted") %>% filter(SAMPLE_NAME %in% (counts_depth_plots %>% filter(SAMPLE_METHOD == "unblasted") %>% filter(SAMPLE_GROUP == "gilot_experimental"))$SAMPLE_NAME)
# data_label_select <- label_data %>% filter(METHOD == "unblasted")  %>% filter(COHORT == "gilot_all_controls") %>% filter(SAMPLE_NAME %in% (counts_depth_plots %>% filter(SAMPLE_METHOD == "unblasted") %>% filter(SAMPLE_GROUP == "gilot_all_controls"))$SAMPLE_NAME)
# counts_depth_plots %>% filter(SAMPLE_METHOD == "unblasted") %>% filter(SAMPLE_GROUP == "gilot_all_controls") %>%   ggplot(
#   aes(
#     x = POS,
#     y = DEPTH
#   )
# ) +
#   geom_point(alpha = 1/20, size = 0.05) +
#   geom_smooth(
#     method = "auto",
#     se = TRUE,
#     fullrange = FALSE,
#     level = 0.95
#   ) + facet_wrap(. ~ SAMPLE_NAME, scales="free_y")  + scale_y_continuous(labels = label_comma()) + scale_x_continuous(labels = label_comma()) + geom_vline(xintercept = c(21563,25384), linetype="dotted", color = "blue", size=0.1) + geom_hline(yintercept = 5000, linetype="dotted", color = "red", size=0.25) +   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +  theme(axis.text.y = element_text(angle = 45, vjust = 1, hjust = 1)) + geom_text(
#     data    = data_label_select,
#     mapping = aes(x = -Inf, y = -Inf, label = paste("Median Depth:",  format(DEPTH_MED,big.mark=",",scientific=FALSE))), #mapping = aes(x = -Inf, y = -Inf, label = paste("CT(S) value:", CT_S, "- Median Depth:",  DEPTH_MED)),
#     # hjust   = 1.05,
#     # vjust   = 1.5
#     hjust   = -0.1,
#     vjust   = -1
#   )
# # + scale_x_continuous(labels = label_comma()) + scale_y_continuous(breaks = seq(0, 250000, by = 50000))
# 
# file <- paste0("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/GILOT_THESIS/figures/","depth_gilot_all_controls_unblasted_pangenome_identifier_depth.png")
# ggsave(file, device = "png", dpi = "retina",   width = 210, height = 148.5 , units = c("mm"), type = "cairo")
# 
# 
# 
# 
# #data_cts_select <- data_cts %>% filter(METHOD == "unblasted") %>% filter(SAMPLE_NAME %in% (counts_depth_plots %>% filter(SAMPLE_METHOD == "unblasted") %>% filter(SAMPLE_GROUP == "gilot_experimental"))$SAMPLE_NAME)
# methodi <- "blasted"
# cohorti <- "gilot_all_controls"
# data_label_select <- label_data %>% filter(METHOD == methodi)  %>% filter(COHORT == cohorti) %>% filter(SAMPLE_NAME %in% (counts_depth_plots %>% filter(SAMPLE_METHOD == methodi) %>% filter(SAMPLE_GROUP == cohorti))$SAMPLE_NAME)
# depth_plot <- counts_depth_plots %>% filter(SAMPLE_METHOD == methodi) %>% filter(SAMPLE_GROUP == cohorti)
# how_many <- depth_plot$SAMPLE_NAME %>% unique() %>% length()
# print(paste(methodi, ":", cohorti, ": n =", how_many))
# plot_object <- depth_plot %>% ggplot(
#   aes(
#     x = POS,
#     y = DEPTH
#   )
# ) +
#   geom_point(alpha = 1/20, size = 0.05) +
#   geom_smooth(
#     method = "auto",
#     se = TRUE,
#     fullrange = FALSE,
#     level = 0.95
#   ) + facet_wrap(. ~ SAMPLE_NAME, scales="free_y")  + scale_y_continuous(labels = unit_format(unit = "k", scale = 1e-3)) + scale_x_continuous(labels = label_comma()) + geom_vline(xintercept = c(21563,25384), linetype="dotted", color = "blue", linewidth=0.1) + geom_hline(yintercept = 5000, linetype="dotted", color = "red", linewidth=0.25) +   theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))  + 
#   geom_text(    data = data_label_select,  mapping = aes(x = -Inf, y = -Inf, label = paste("Median Depth:",  format(DEPTH_MED,big.mark=",",scientific=FALSE))), #mapping = aes(x = -Inf, y = -Inf, label = paste("CT(S) value:", CT_S, "- Median Depth:",  DEPTH_MED)),
#                 # +  theme(axis.text.y = element_text(angle = 45, vjust = 1, hjust = 1))
#                 #    mapping = aes(x = -Inf, y = -Inf, label = paste("Median Depth:",  format(DEPTH_MED,big.mark=",",scientific=FALSE))), #mapping = aes(x = -Inf, y = -Inf, label = paste("CT(S) value:", CT_S, "- Median Depth:",  DEPTH_MED)),
#                 # hjust   = 1.05,
#                 # vjust   = 1.5
#                 hjust   = -0.1,
#                 vjust   = -1
#   )
# # + scale_x_continuous(labels = label_comma()) + scale_y_continuous(breaks = seq(0, 250000, by = 50000))
# 
# file <- paste0("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/GILOT_THESIS/figures/","depth_", cohorti, "_", methodi, "_pangenome_identifier.png")
# ggsave(file, device = "png", dpi = 900,  type = "cairo") # , dpi = "retina",  type = "cairo") # width = 210, height = 148.5 , units = c("mm")
# 
# 







# Cairo::Cairo(
#   30, #length
#   30, #width
#   file = file, #paste("nameofplot", ".png", sep = ""),
#   type = "png", #tiff
#   bg = "transparent", #white or transparent depending on your requirement 
#   dpi = 300,
#   units = "cm" #you can change to pixels etc 
# )
# plot(p) #p is your graph object 
# dev.off()
