#######################  

master_performance_consensus_super <- master_performance_consensus

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
