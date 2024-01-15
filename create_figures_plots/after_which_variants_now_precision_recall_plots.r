# /mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/create_figures_plots/wrapper_thesis.r
# /mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/create_figures_plots/depth_plots.r
# /mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/create_figures_plots/which_variants_to_consider_for_pre_rec.r

# precision recall --------------------------------------------------------
master_performance_consensus_super_concordant_confusion_matrix <- readRDS("/mnt/bulk_data/bgilot/GitHub/2022_EKUT_THESIS_BIOINFO/data/data_unified/master_performance_consensus_super_concordant_confusion_matrix.rds")

master_performance_consensus_super_concordant_confusion_matrix %>% glimpse()


label_x <- "Concordant Recall"
label_y <- "Concordant Precision"
label_caption <- NULL
master_performance_consensus_super_concordant_confusion_matrix %>%
  filter(plot_target != "C2_Ancestral") %>% 
  #filter(RELEVANT_COPIES  >= 10000, RELEVANT_COPIES  <= 100000000) %>% 
  #filter(RELEVANT_FRACTION  >= 0.01, RELEVANT_FRACTION  <= 1) %>%
  #filter(plot_target != "C2_Ancestral") %>% 
  #filter(control_by_caller == "umivar") %>% 
  #filter(false_negative_snp_len == 0) %>%
  #filter(loop_control_variant_text_label == "C1_Australia") %>% 
  ggplot(aes(x = concordant_recall,
             y = concordant_precision)) +
  geom_point(alpha = 20 / 20, size = 0.75) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  facet_wrap(. ~ plot_target) +
  labs(x = label_x) +
  labs(y = label_y) +
  labs(caption = label_caption) +
  #stat_poly_line() +
  geom_smooth(
    method = "auto",
    se = TRUE,
    fullrange = FALSE,
    level = 0.95,
    linewidth = 0.5,
    alpha = 0.5,
    color = "yellow"
  )
#210 x 297 mm
file <- paste0(folder_plots,"concordant_by_target", ".png")
#ggsave(file, device = "png", dpi = 1800,  type = "cairo", width = 297, height = 210, units = "mm") 
ggsave(file, device = "png", dpi = 1800,  type = "cairo", width = 400, height = 400, units = "mm") 
print(file)


label_x <- "Recall"
label_y <- "Precision"
label_caption <- NULL
master_performance_consensus_super_concordant_confusion_matrix %>%
  #filter(RELEVANT_COPIES  >= 10000, RELEVANT_COPIES  <= 100000000) %>% 
  #filter(RELEVANT_FRACTION  >= 0.01, RELEVANT_FRACTION  <= 1) %>%
  filter(plot_target != "C2_Ancestral") %>% 
  #filter(control_by_caller == "umivar") %>% 
  #filter(false_negative_snp_len == 0) %>%
  #filter(loop_control_variant_text_label == "C1_Australia") %>% 
  ggplot(aes(x = recall,
             y = precision)) +
  geom_point(alpha = 20 / 20, size = 0.75) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  facet_wrap(. ~ control_by_caller) + #, scales = "free") +
  labs(x = label_x) +
  labs(y = label_y) +
  labs(caption = label_caption) +
  #stat_poly_line() +
  geom_smooth(
    method = "auto",
    se = TRUE,
    fullrange = FALSE,
    level = 0.95,
    linewidth = 0.5,
    alpha = 0.5,
    color = "yellow"
  )

file <- paste0(folder_plots,"non_concordant_by_caller", ".png")
ggsave(file, device = "png", dpi = 1800,  type = "cairo", width = 400, height = 400, units = "mm") 
print(file)




label_x <- "Recall"
label_y <- "Precision"
label_caption <- NULL
master_performance_consensus_super_concordant_confusion_matrix %>%
  #filter(RELEVANT_COPIES  >= 10000, RELEVANT_COPIES  <= 100000000) %>% 
  #filter(RELEVANT_FRACTION  >= 0.01, RELEVANT_FRACTION  <= 1) %>%
  filter(plot_target != "C2_Ancestral") %>% 
  #filter(control_by_caller == "umivar") %>% 
  #filter(false_negative_snp_len == 0) %>%
  #filter(loop_control_variant_text_label == "C1_Australia") %>% 
  ggplot(aes(x = recall,
             y = precision)) +
  geom_point(alpha = 20 / 20, size = 0.75) +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  facet_grid(plot_target ~ control_by_caller) + #, scales = "free") +
  labs(x = label_x) +
  labs(y = label_y) +
  labs(caption = label_caption) +
  stat_poly_line() +
  stat_poly_eq()
  # geom_smooth(
  #   method = "auto",
  #   se = TRUE,
  #   fullrange = FALSE,
  #   level = 0.95,
  #   linewidth = 0.5,
  #   alpha = 0.5,
  #   color = "yellow"
  # )

file <- paste0(folder_plots,"non_concordant_by_caller_by_target", ".png")
#ggsave(file, device = "png", dpi = 1800,  type = "cairo", width = 400, height = 400, units = "mm")
ggsave(file, device = "png", dpi = 1800,  type = "cairo", width = 297, height = 210, units = "mm") 
print(file)



master_performance_consensus_super_concordant_confusion_matrix %>% 
  filter(plot_target != "C2_Ancestral") %>% 
    filter(RELEVANT_COPIES  >= 1000, RELEVANT_COPIES  <= 1000000) %>%
    filter(RELEVANT_FRACTION  >= 0.0005, RELEVANT_FRACTION  <= 0.5) %>%
    # filter(control_by_caller == "umivar") %>%
    # filter(false_negative_snp_len == 0) %>%
    #filter(loop_control_variant_text_label == "C1_Australia") %>%
    ggplot(aes(x = recall,
               y = precision)) +
    geom_point(alpha = 15 / 20, size = 0.5) +
  stat_poly_line() +
    facet_grid(control_by_caller ~ RELEVANT_FRACTION) 

master_performance_consensus_super_concordant_confusion_matrix %>% 
  #filter(plot_target != "C2_Ancestral") %>% 
  #filter(RELEVANT_COPIES  >= 1000, RELEVANT_COPIES  <= 1000000) %>%
  #filter(RELEVANT_FRACTION  >= 0.0005, RELEVANT_FRACTION  <= 0.5) %>%
  # filter(control_by_caller == "umivar") %>%
  # filter(false_negative_snp_len == 0) %>%
  #filter(loop_control_variant_text_label == "C1_Australia") %>%
  ggplot(aes(x = aware_median_depth,
             y = precision)) +
  geom_point(alpha = 15 / 20, size = 0.5) +
  stat_poly_line() +
  facet_wrap(. ~ control_by_caller)

master_performance_consensus_super_concordant_confusion_matrix %>% 
  #filter(control_by_caller != "umivar") %>%
  #filter(plot_target != "C2_Ancestral") %>% 
  #filter(RELEVANT_COPIES  >= 1000, RELEVANT_COPIES  <= 1000000) %>%
  #filter(RELEVANT_FRACTION  >= 0.0005, RELEVANT_FRACTION  <= 0.5) %>%
  # filter(control_by_caller == "umivar") %>%
  # filter(false_negative_snp_len == 0) %>%
  #filter(loop_control_variant_text_label == "C1_Australia") %>%
  ggplot(aes(x = RELEVANT_COPIES,
             y = precision)) +
  geom_point(alpha = 15 / 20, size = 0.5) +
  stat_poly_line() +
  facet_wrap(. ~ control_by_caller)

master_performance_consensus_super_concordant_confusion_matrix %>% 
  #filter(control_by_caller != "umivar") %>%
  #filter(plot_target != "C2_Ancestral") %>% 
  #filter(RELEVANT_COPIES  >= 1000, RELEVANT_COPIES  <= 1000000) %>%
  #filter(RELEVANT_FRACTION  >= 0.0005, RELEVANT_FRACTION  <= 0.5) %>%
  # filter(control_by_caller == "umivar") %>%
  # filter(false_negative_snp_len == 0) %>%
  #filter(loop_control_variant_text_label == "C1_Australia") %>%
  ggplot(aes(x = concordant_precision,
             y = concordant_recall)) +
  geom_point(alpha = 15 / 20, size = 0.5) +
  stat_poly_line() +
  facet_wrap(. ~ control_by_caller)

master_performance_consensus_super_concordant_confusion_matrix %>% 
  filter(plot_target != "C2_Ancestral") %>% 
  #filter(RELEVANT_COPIES %in% (master_performance_consensus_super_concordant_confusion_matrix$RELEVANT_COPIES %>% unique() %>% sort() %>% tail())) %>%
  filter(RELEVANT_COPIES  >= 1000, RELEVANT_COPIES  <= 1000000) %>%
  filter(RELEVANT_FRACTION  >= 0.005, RELEVANT_FRACTION  <= 0.5) %>%
  # filter(control_by_caller == "umivar") %>%
  # filter(false_negative_snp_len == 0) %>%
  #filter(loop_control_variant_text_label == "C1_Australia") %>%
  ggplot(aes(x = recall,
             y = precision)) +
  geom_point(alpha = 15 / 20, size = 0.5) +
  stat_poly_line() +
  #stat_poly_eq() +
  facet_grid(control_by_caller ~ RELEVANT_FRACTION) 


master_performance_consensus_super_concordant_confusion_matrix %>% 
  filter(plot_target != "C2_Ancestral") %>% 
  filter(RELEVANT_FRACTION  >= 0.5, RELEVANT_FRACTION  <= 1) %>% nrow()

master_performance_consensus_super_concordant_confusion_matrix %>% 
  filter(plot_target != "C2_Ancestral") %>% 
  #filter(RELEVANT_COPIES %in% (master_performance_consensus_super_concordant_confusion_matrix$RELEVANT_COPIES %>% unique() %>% sort() %>% tail())) %>%
  #filter(RELEVANT_COPIES  >= 1000, RELEVANT_COPIES  <= 1000000) %>%
  filter(RELEVANT_FRACTION  >= 0.5, RELEVANT_FRACTION  <= 1.0) %>%
  # filter(control_by_caller == "umivar") %>%
  # filter(false_negative_snp_len == 0) %>%
  #filter(loop_control_variant_text_label == "C1_Australia") %>%
  ggplot(aes(x = RELEVANT_FRACTION,
             y = tp_median)) +
  geom_point(alpha = 15 / 20, size = 0.5) +
  stat_poly_line() +
  stat_poly_eq() +
  facet_wrap(. ~ control_by_caller) 

master_performance_consensus_super_concordant_confusion_matrix %>% 
  #filter(plot_target == "C48_Omicron") %>% 
  #filter(RELEVANT_COPIES %in% (master_performance_consensus_super_concordant_confusion_matrix$RELEVANT_COPIES %>% unique() %>% sort() %>% tail())) %>%
  # filter(RELEVANT_COPIES  >= 1000, RELEVANT_COPIES  <= 500000000) %>%
  # filter(RELEVANT_FRACTION  >= 0, RELEVANT_FRACTION  <= 1) %>%
  #filter(control_by_caller == "umivar") %>%
  # filter(false_negative_snp_len == 0) %>%
  #filter(loop_control_variant_text_label == "C1_Australia") %>%
  ggplot(aes(x = RELEVANT_COPIES,
             y = concordant_recall)) +
  geom_point(alpha = 15 / 20, size = 0.5) +
  stat_poly_line() +
  # stat_poly_eq() +
  facet_wrap(. ~ control_by_caller) 

master_performance_consensus_super_concordant_confusion_matrix %>% filter(RELEVANT_FRACTION  >= 0.005, RELEVANT_FRACTION  <= 1) %>% ggplot(aes(x=af_mean, y=tp_mean)) + 
  geom_boxplot() +
  facet_wrap(. ~ RELEVANT_FRACTION)

master_performance_consensus_super_concordant_confusion_matrix %>% filter(RELEVANT_FRACTION  >= 0.00, RELEVANT_FRACTION  <= 0.5) %>% ggplot(aes(RELEVANT_FRACTION, tp_mean)) + geom_boxplot() +
  facet_wrap(. ~ RELEVANT_FRACTION)


master_performance_consensus_super_concordant_confusion_matrix %>% 
  filter(plot_target != "C2_Ancestral") %>% 
  #filter(RELEVANT_COPIES %in% (master_performance_consensus_super_concordant_confusion_matrix$RELEVANT_COPIES %>% unique() %>% sort() %>% tail())) %>%
  #filter(RELEVANT_COPIES  >= 1000, RELEVANT_COPIES  <= 1000000) %>%
  filter(RELEVANT_FRACTION  >= 0.5, RELEVANT_FRACTION  <= 1) %>%
  # filter(control_by_caller == "umivar") %>%
  # filter(false_negative_snp_len == 0) %>%
  #filter(loop_control_variant_text_label == "C1_Australia") %>%
  ggplot(aes(x = RELEVANT_FRACTION,
             y = tp_median)) +
  geom_point(alpha = 15 / 20, size = 0.5) +
  stat_poly_line() +
  facet_wrap(. ~ control_by_caller) 



master_performance_consensus_super_concordant_confusion_matrix %>% 
  filter(RELEVANT_COPIES  >= 00, RELEVANT_COPIES  <= 100000000) %>% 
  filter(RELEVANT_FRACTION  >= 0, RELEVANT_FRACTION  <= 1) %>% 
  # filter(control_by_caller == "umivar") %>% 
  # filter(false_negative_snp_len == 0) %>%
  #filter(loop_control_variant_text_label == "C1_Australia") %>% 
  ggplot(aes(x = recall,
             y = precision)) +
  geom_point(alpha = 20 / 20, size = 0.75) + 
  facet_wrap(. ~ plot_caller)  +
  geom_smooth(
    method = "auto",
    se = TRUE,
    fullrange = FALSE,
    level = 0.95,
    linewidth = 0.5,
    alpha = 0.5,
    color = "yellow"
  )   + 
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) 


for_plot <- master_performance_consensus_super_concordant_confusion_matrix %>% select(Identifier, blind_median_depth, aware_median_depth, RELEVANT_COPIES) %>% unique()

for_plot_n <- for_plot  %>% nrow()

for_plot %>% ggplot(aes(x = blind_median_depth,
                        y = aware_median_depth)) +
  geom_point(alpha = 20 / 20, size = 0.75)   +
  geom_smooth(
    method = "auto",
    se = TRUE,
    fullrange = FALSE,
    level = 0.95,
    linewidth = 0.5,
    alpha = 0.5,
    color = "yellow"
  )   +   stat_poly_line() +
  stat_poly_eq() +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) 


master_performance_consensus_super %>% 
  filter(RELEVANT_COPIES  >= 1000, RELEVANT_COPIES  <= 1000000) %>% 
  #filter(RELEVANT_FRACTION  >= 0, RELEVANT_FRACTION  <= 0.5) %>% 
  #filter(plot_af_mean  >= 0, RELEVANT_FRACTION  <= 0.75) %>% 
  # filter(control_by_caller == "umivar") %>% 
  # filter(false_negative_snp_len == 0) %>%
  #filter(loop_control_variant_text_label == "C1_Australia") %>% 
  ggplot(aes(x = plot_af_mean,
             y = plot_fraction)) +
  geom_point(alpha = 20 / 20, size = 0.75) + 
  facet_grid(RELEVANT_COPIES ~ plot_caller)  +
  stat_poly_line() +
  stat_poly_eq() +
scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) 
+
  geom_smooth(
    method = "auto",
    se = TRUE,
    fullrange = FALSE,
    level = 0.95,
    linewidth = 0.5,
    alpha = 0.5,
    color = "yellow"
  )   + 
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) 



master_performance_consensus_super %>% 
  filter(RELEVANT_COPIES  >= 0, RELEVANT_COPIES  <= 1000000) %>% 
  #filter(RELEVANT_FRACTION  >= 0.0001, RELEVANT_FRACTION  <= 0.5) %>% 
  # filter(control_by_caller == "umivar") %>% 
  # filter(false_negative_snp_len == 0) %>%
  #filter(loop_control_variant_text_label == "C1_Australia") %>% 
  ggplot(aes(x = RELEVANT_FRACTION,
             y = plot_precision)) +
  geom_point(alpha = 10 / 20, size = 0.75) + 
  stat_poly_line() +
  stat_poly_eq() +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0))  +
  facet_wrap(. ~ plot_caller)
+
  geom_smooth(
    method = "auto",
    se = TRUE,
    fullrange = FALSE,
    level = 0.95,
    linewidth = 0.5,
    alpha = 0.5,
    color = "yellow"
  ) +
facet_wrap(. ~ plot_caller)  +
  stat_poly_eq()  




master_performance_consensus_super %>% ggplot(aes(x = recall,
                                  y = precision)) +
  geom_point(alpha = 10 / 20, size = 0.5) + 
  facet_wrap(control_by_caller ~ THESIS_SYNTHETIC_PRESENT_C48_OMICRON_EPI_ISL_6841980_COPIES)  +
  geom_smooth(
    method = "auto",
    se = TRUE,
    fullrange = FALSE,
    level = 0.95,
    linewidth = 0.5,
    alpha = 0.5,
    color = "yellow"
  ) 

master_performance_consensus_super$factor_RELEVANT_FRACTION <- master_performance_consensus_super$RELEVANT_FRACTION %>% as.factor()

master_performance_consensus_super  %>% 
  filter(RELEVANT_COPIES  >= 1000, RELEVANT_COPIES  <= 1000000) %>% 
  filter(RELEVANT_FRACTION  >= 0.0000, RELEVANT_FRACTION  <= 1) %>% ggplot( aes(x=RELEVANT_FRACTION, y=plot_af_mean)) + 
  geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=2, notch=FALSE) + stat_summary(fun.y=mean, geom="point", shape=23, size=4) + 
  geom_dotplot(binaxis='y', stackdir='center', dotsize=0.1) + scale_x_discrete(limits=c("0.0", "1.0")) +
  facet_wrap(. ~ RELEVANT_FRACTION)

file <- paste0(folder_plots,"af_accuracy_by_caller_bad_label", ".png")
ggsave(file, device = "png", dpi = 1800,  type = "cairo", width = 400, height = 400, units = "mm") 
print(file)


master_performance_consensus_super  %>% 
  filter(plot_target != "C2_Ancestral") %>% 
  #filter(RELEVANT_COPIES  >= 1000, RELEVANT_COPIES  <= 1000000) %>% 
  filter(RELEVANT_FRACTION  >= 0.0000, RELEVANT_FRACTION  <= 0.5) %>% ggplot( aes(x=RELEVANT_FRACTION, y=plot_af_mean)) + 
  geom_point(alpha = 10 / 20, size = 0.5) +
  stat_poly_line() +
  stat_poly_eq() +
  facet_wrap(. ~ control_by_caller) 

file <- paste0(folder_plots,"af_accuracy_by_caller_bad_label", ".png")
ggsave(file, device = "png", dpi = 1800,  type = "cairo", width = 400, height = 400, units = "mm") 
print(file)


master_performance_consensus_super  %>% filter(abs(plot_af_mean - RELEVANT_FRACTION ) < 0.1) %>% 
  #filter(RELEVANT_COPIES  >= 1000, RELEVANT_COPIES  <= 1000000) %>% 
  filter(plot_target != "C2_Ancestral") %>% 
  filter(RELEVANT_FRACTION  >= 0.0000, RELEVANT_FRACTION  <= 0.5) %>% ggplot( aes(x=RELEVANT_FRACTION, y=plot_af_mean)) + 
  geom_point(alpha = 10 / 20, size = 0.5) +
  stat_poly_line() +
  stat_poly_eq() +
  facet_wrap(. ~ control_by_caller) #+
  # geom_smooth(
  #   method = "auto",
  #   se = TRUE,
  #   fullrange = FALSE,
  #   level = 0.95,
  #   linewidth = 0.5,
  #   alpha = 0.5,
  #   color = "yellow"
  # ) 
  # geom_boxplot(outlier.colour="black", outlier.shape=16, outlier.size=2, notch=FALSE) + stat_summary(fun.y=mean, geom="point", shape=23, size=4) + 
  # geom_dotplot(binaxis='y', stackdir='center', dotsize=0.1) + scale_x_discrete(limits=c("0.0", "1.0")) +
  # facet_wrap(. ~ RELEVANT_FRACTION)

file <- paste0(folder_plots,"af_accuracy_by_caller_good_label", ".png")
ggsave(file, device = "png", dpi = 1800,  type = "cairo", width = 400, height = 400, units = "mm") 
print(file)


master_performance_consensus_super %>% filter(abs(plot_af_mean - RELEVANT_FRACTION ) < 0.1) %>% select(Identifier) %>% unique()


master_performance_consensus_super_concordant_confusion_matrix %>% select(blind_median_depth, aware_median_depth, THESIS_SYNTHETIC_COPIES_TOTAL) %>% unique() %>%
  ggplot(aes(x=blind_median_depth, y=aware_median_depth)) + geom_point(alpha = 10 / 20, size = 0.5) +
labs(x = "blind_median_depth") +
  labs(y = "aware_median_depth") +
  labs(caption = "label_caption") +
  stat_poly_line() +
  stat_poly_eq() +
  geom_smooth(
    method = "auto",
    se = TRUE,
    fullrange = FALSE,
    level = 0.95,
    linewidth = 0.5,
    alpha = 0.5,
    color = "yellow"
  ) 

master_performance_consensus_super_concordant_confusion_matrix %>% select(blind_median_depth, aware_median_depth, THESIS_SYNTHETIC_COPIES_TOTAL) %>% unique() %>%
  ggplot(aes(x=THESIS_SYNTHETIC_COPIES_TOTAL, y=aware_median_depth)) + geom_point(alpha = 10 / 20, size = 0.5) +
  labs(x = "THESIS_SYNTHETIC_COPIES_TOTAL") +
  labs(y = "aware_median_depth") +
  labs(caption = "label_caption") +
  stat_poly_line() +
  stat_poly_eq() +
  geom_smooth(
    method = "auto",
    se = TRUE,
    fullrange = FALSE,
    level = 0.95,
    linewidth = 0.5,
    alpha = 0.5,
    color = "yellow"
  ) 


master_performance_consensus_super_concordant_confusion_matrix %>% select(blind_median_depth, aware_median_depth, THESIS_SYNTHETIC_COPIES_TOTAL) %>% unique() %>%
  ggplot(aes(x=THESIS_SYNTHETIC_COPIES_TOTAL, y=blind_median_depth)) + geom_point(alpha = 10 / 20, size = 0.5) +
  labs(x = "THESIS_SYNTHETIC_COPIES_TOTAL") +
  labs(y = "blind_median_depth") +
  labs(caption = "label_caption") +
  stat_poly_line() +
  stat_poly_eq() +
  geom_smooth(
    method = "auto",
    se = TRUE,
    fullrange = FALSE,
    level = 0.95,
    linewidth = 0.5,
    alpha = 0.5,
    color = "yellow"
  ) 


human_samples
ct_vs_mean_depth

ct_vs_mean_depth%>%
  ggplot(aes(x=Median_blind , y=Median_aware)) + geom_point(alpha = 10 / 20, size = 0.5) +
  labs(x = "Median_blind") +
  labs(y = "Median_aware") +
  labs(caption = "label_caption") +
  stat_poly_line() +
  stat_poly_eq() +
  geom_smooth(
    method = "auto",
    se = TRUE,
    fullrange = FALSE,
    level = 0.95,
    linewidth = 0.5,
    alpha = 0.5,
    color = "yellow"
  ) 

#################################################################

label_x <- "CT(S) Value"
label_y <- bquote("Median Coverage"^aware)
label_caption <- paste("n =", nrow(ct_vs_mean_depth))
ct_vs_mean_depth%>%
  ggplot(aes(x=CT_S , y=Median_aware)) + geom_point(alpha = 10 / 20, size = 0.95) +
  labs(x = "CT(S) Value") +
  labs(y = label_y) +
  labs(caption = label_caption) +
  #stat_poly_line() +
  stat_poly_eq() +
  scale_y_continuous(labels = unit_format(unit = "k", scale = 1e-3), expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_x_reverse()  +
  geom_smooth(
    method = "auto",
    se = TRUE,
    fullrange = FALSE,
    level = 0.95,
    linewidth = 0.5,
    alpha = 0.5,
    color = "yellow"
  ) 


file <- paste0(folder_plots,"clinical_samples_cts_vs_aware", ".png")
ggsave(file, device = "png", dpi = 1800,  type = "cairo", width = 210, height = 210, units = "mm") 
print(file)

################################################################################
# Compare depth
################################################################################

label_x <- "CT(S) Value"
label_y <- bquote("Median Coverage"^aware)
label_caption <- paste0("Clinical Specimens (n = ", nrow(ct_vs_mean_depth),")")
ct_vs_mean_depth%>%
  ggplot(aes(x=CT_S , y=Median_aware)) + geom_point(alpha = 10 / 20, size = 0.5) +
  labs(x = label_x) +
  labs(y = label_y) +
  labs(caption = label_caption) +
  #stat_poly_line() +
  stat_poly_eq() +
  scale_y_continuous(labels = unit_format(unit = "k", scale = 1e-3), expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_x_reverse()  +
  geom_smooth(
    method = "auto",
    se = TRUE,
    fullrange = FALSE,
    level = 0.95,
    linewidth = 0.5,
    alpha = 0.5,
    color = "yellow"
  ) 

file <- paste0(folder_plots,"clinical_samples_cts_vs_aware", ".png")
ggsave(file, device = "png", dpi = 1800,  type = "cairo", width = 400, height = 400, units = "mm") 
print(file)

plot_tibble <- master_performance_consensus_super_concordant_confusion_matrix %>% select(blind_median_depth, aware_median_depth, THESIS_SYNTHETIC_COPIES_TOTAL) %>% unique()
label_x <- bquote("Million Viral Copies / ul")
label_y <- bquote("Median Coverage"^aware)
label_caption <- paste0("Synthetic Specimens (n = ", nrow(plot_tibble),")")
plot_tibble %>%
  ggplot(aes(x=THESIS_SYNTHETIC_COPIES_TOTAL, y=aware_median_depth)) + geom_point(alpha = 10 / 20, size = 0.5) +
  labs(x = label_x) +
  labs(y = label_y) +
  labs(caption = label_caption) +
  stat_poly_line() +
  stat_poly_eq() +
  scale_y_continuous(labels = unit_format(unit = "k", scale = 1e-3), expand = c(0, 0)) +
  scale_x_continuous(labels = unit_format(unit = "", scale = 1e-6), expand = c(0, 0))

file <- paste0(folder_plots,"synthetic_samples_copies_vs_aware", ".png")
ggsave(file, device = "png", dpi = 1800,  type = "cairo", width = 400, height = 400, units = "mm") 
print(file)

################################################################################

label_x <- "CT(S) Value"
label_y <- bquote("Median Coverage"^blind)
label_caption <- paste0("Clinical Specimens (n = ", nrow(ct_vs_mean_depth),")")
ct_vs_mean_depth%>%
  ggplot(aes(x=CT_S , y=Median_blind)) + geom_point(alpha = 10 / 20, size = 0.5) +
  labs(x = "CT(S) Value") +
  labs(y = label_y) +
  labs(caption = label_caption) +
  #stat_poly_line() +
  stat_poly_eq() +
  scale_y_continuous(labels = unit_format(unit = "k", scale = 1e-3), expand = c(0, 0)) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_x_reverse()  +
  geom_smooth(
    method = "auto",
    se = TRUE,
    fullrange = FALSE,
    level = 0.95,
    linewidth = 0.5,
    alpha = 0.5,
    color = "yellow"
  ) 

file <- paste0(folder_plots,"clinical_samples_cts_vs_blind", ".png")
ggsave(file, device = "png", dpi = 1800,  type = "cairo", width = 400, height = 400, units = "mm") 
print(file)

plot_tibble <- master_performance_consensus_super_concordant_confusion_matrix %>% select(blind_median_depth, aware_median_depth, THESIS_SYNTHETIC_COPIES_TOTAL) %>% unique()
label_x <- bquote("Million Viral Copies / ul")
label_y <- bquote("Median Coverage"^blind)
label_caption <- paste0("Synthetic Specimens (n = ", nrow(plot_tibble),")")
plot_tibble %>%
  ggplot(aes(x=THESIS_SYNTHETIC_COPIES_TOTAL, y=blind_median_depth)) + geom_point(alpha = 10 / 20, size = 0.5) +
  labs(x = label_x) +
  labs(y = label_y) +
  labs(caption = label_caption) +
  stat_poly_line() +
  stat_poly_eq() +
  scale_y_continuous(labels = unit_format(unit = "k", scale = 1e-3), expand = c(0, 0)) +
  scale_x_continuous(labels = unit_format(unit = "", scale = 1e-6), expand = c(0, 0))

file <- paste0(folder_plots,"synthetic_samples_copies_vs_blind", ".png")
ggsave(file, device = "png", dpi = 1800,  type = "cairo", width = 400, height = 400, units = "mm") 
print(file)

#############################################################################

label_x <- bquote("Median Coverage"^blind)
label_y <- bquote("Median Coverage"^aware)
label_caption <- paste0("Clinical Specimens (n = ", nrow(ct_vs_mean_depth),")")
ct_vs_mean_depth%>%
  ggplot(aes(x=Median_blind , y=Median_aware)) + geom_point(alpha = 10 / 20, size = 0.5) +
  labs(x = label_x) +
  labs(y = label_y) +
  labs(caption = label_caption) +
  stat_poly_line() +
  stat_poly_eq() +
  scale_y_continuous(labels = unit_format(unit = "k", scale = 1e-3), expand = c(0, 0)) +
  scale_x_continuous(labels = unit_format(unit = "k", scale = 1e-3), expand = c(0, 0))

file <- paste0(folder_plots,"clinical_samples_blind_vs_aware", ".png")
ggsave(file, device = "png", dpi = 1800,  type = "cairo", width = 400, height = 400, units = "mm") 
print(file)


plot_tibble <- master_performance_consensus_super_concordant_confusion_matrix %>% select(blind_median_depth, aware_median_depth, THESIS_SYNTHETIC_COPIES_TOTAL) %>% unique()
label_x <- bquote("Median Coverage"^blind)
label_y <- bquote("Median Coverage"^aware)
label_caption <- paste0("Synthetic Specimens (n = ", nrow(plot_tibble),")")
plot_tibble %>%
  ggplot(aes(x=blind_median_depth, y=aware_median_depth)) + geom_point(alpha = 10 / 20, size = 0.5) +
  labs(x = label_x) +
  labs(y = label_y) +
  labs(caption = label_caption) +
  stat_poly_line() +
  stat_poly_eq() +
  scale_y_continuous(labels = unit_format(unit = "k", scale = 1e-3), expand = c(0, 0)) +
  scale_x_continuous(labels = unit_format(unit = "k", scale = 1e-3), expand = c(0, 0))

file <- paste0(folder_plots,"synthetic_samples_blind_vs_aware", ".png")
ggsave(file, device = "png", dpi = 1800,  type = "cairo", width = 400, height = 400, units = "mm") 
print(file)


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