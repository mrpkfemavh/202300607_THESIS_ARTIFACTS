library(tidyverse, lib.loc = "/usr/lib/R/site-library")
library(readr)
library(ggplot2)
library(stringr)
library(DescTools)


controls_data <- read_csv("Downloads/controls_data.csv", 
                          col_types = cols(position_position = col_integer(), 
                                           position_depth = col_integer(), corrected_count = col_integer(), 
                                           sum_Reference = col_integer(), sum_A = col_integer(), 
                                           sum_C = col_integer(), sum_T = col_integer(), 
                                           sum_G = col_integer(), position_ins = col_integer(), 
                                           position_del = col_integer()))
summary(controls_data$corrected_count)

controls_data$C1 <- as.logical(NA)
controls_data$C2 <- as.logical(NA)
controls_data$C14 <- as.logical(NA)

gilot_present_final <- read.csv("~/Downloads/gilot_present_final.csv")

joined_controls <- left_join(controls_data, gilot_present_final, by = c("file_identifier" = "j_name"))

joined_controls$corrected_count %>% summary()
 
plot_data <- joined_controls %>% filter(corrected_count >= 2944)  %>% select("file_identifier", "j_name_external", "C1", "C2", "C14", "position_position", "corrected_count", "position_reference", "sum_Reference", "fraction_Reference", "sum_A", "fraction_A", "sum_C", "fraction_C", "sum_T", "fraction_T", "sum_G", "fraction_G", "position_ins", "fraction_insertions", "position_del", "fraction_deletions" )%>% arrange(position_position ,desc(corrected_count))


#write.csv(x=plot_data, file="/home/bgilot/Downloads/to_fix.csv")

fixed_data <- read_csv("Downloads/to_fix.csv") # %>% filter(corrected_count >= 10000) %>% filter(percent_c1 >= 0) %>% filter(percent_c1 <= 10)
fixed_data$colour <- as.factor(fixed_data$position_position)
fixed_data$variant <- 100 - fixed_data$fraction_Reference

# theme(axis.title.x = element_text(vjust = 0, size = 15),axis.title.y = element_text(vjust = 2, size = 15)) +

fixed_data_plot <- fixed_data %>% filter(percent_c1 >= 0) %>% filter(percent_c1 <= 100) %>% filter(position_position != 29754)
ggplot(fixed_data_plot, aes(percent_c1, variant, color = colour))  + geom_point(color = fixed_data_plot$colour, shape = fixed_data_plot$colour, size = 5) +   stat_summary(fun.data=mean_cl_normal) + geom_smooth(method='lm') +
  labs(x = "% Admixture - [MT007544.1]", y = "% Observed - [MT007544.1]",
       title = "Twist Synthetic SARS-CoV-2 RNA Control 1 (MT007544.1)",
       subtitle = "Detection of [19065 :: T > C], [22303 :: T > G], [26144 :: G > T]",
       caption = "unpublished data",
       color = "Genomic Position") + facet_wrap(~ colour) 

fixed_data_plot <- fixed_data %>% filter(percent_c1 >= 95) %>% filter(percent_c1 <= 100) %>% filter(position_position != 29754)
ggplot(fixed_data_plot, aes(percent_c1, variant, color = colour))  + geom_point(color = fixed_data_plot$colour, shape = fixed_data_plot$colour, size = 5) +   stat_summary(fun.data=mean_cl_normal) + geom_smooth(method='lm') + 
  labs(x = "% Admixture - [MT007544.1]", y = "% Observed - [MT007544.1]",
       title = "Twist Synthetic SARS-CoV-2 RNA Control 1 (MT007544.1)",
       subtitle = "Detection of [19065 :: T > C], [22303 :: T > G], [26144 :: G > T]",
       caption = "unpublished data",
       color = "Genomic Position") + facet_wrap(~ colour) 


fixed_data_plot <- fixed_data %>% filter(percent_c1 >= 0) %>% filter(percent_c1 <= 1) %>% filter(position_position != 29754)
ggplot(fixed_data_plot, aes(percent_c1, variant, color = colour))  + geom_point(color = fixed_data_plot$colour, shape = fixed_data_plot$colour, size = 5) +   stat_summary(fun.data=mean_cl_normal) + geom_smooth(method='lm') + 
  labs(x = "% Admixture - [MT007544.1]", y = "% Observed - [MT007544.1]",
       title = "Twist Synthetic SARS-CoV-2 RNA Control 1 (MT007544.1)",
       subtitle = "Detection of [19065 :: T > C], [22303 :: T > G], [26144 :: G > T]",
       caption = "unpublished data",
      color = "Genomic Position") + facet_wrap(~ colour) 

c("/home/bgilot/Downloads/20076a138_01_bowtie2_unblasted_sorted_corrected_sorted_deduplicated_sorted_samtools_mpileup_depth.tsv",
  "/home/bgilot/Downloads/20076a139_01_bowtie2_unblasted_sorted_corrected_sorted_deduplicated_sorted_samtools_mpileup_depth.tsv",
  "/home/bgilot/Downloads/20076a140_01_bowtie2_unblasted_sorted_corrected_sorted_deduplicated_sorted_samtools_mpileup_depth.tsv",
  "/home/bgilot/Downloads/20076a141_01_bowtie2_unblasted_sorted_corrected_sorted_deduplicated_sorted_samtools_mpileup_depth.tsv",
  "/home/bgilot/Downloads/20076a142_01_bowtie2_unblasted_sorted_corrected_sorted_deduplicated_sorted_samtools_mpileup_depth.tsv",
  "/home/bgilot/Downloads/20076a143_01_bowtie2_unblasted_sorted_corrected_sorted_deduplicated_sorted_samtools_mpileup_depth.tsv",
  "/home/bgilot/Downloads/20076a144_01_bowtie2_unblasted_sorted_corrected_sorted_deduplicated_sorted_samtools_mpileup_depth.tsv",
  "/home/bgilot/Downloads/20076a145_01_bowtie2_unblasted_sorted_corrected_sorted_deduplicated_sorted_samtools_mpileup_depth.tsv",
  "/home/bgilot/Downloads/20076a146_01_bowtie2_unblasted_sorted_corrected_sorted_deduplicated_sorted_samtools_mpileup_depth.tsv",
  "/home/bgilot/Downloads/20076Pos_01_bowtie2_unblasted_sorted_corrected_sorted_deduplicated_sorted_samtools_mpileup_depth.tsv",
  "/home/bgilot/Downloads/21014a011_01_bowtie2_unblasted_sorted_corrected_sorted_deduplicated_sorted_samtools_mpileup_depth.tsv",
  "/home/bgilot/Downloads/21014a012_01_bowtie2_unblasted_sorted_corrected_sorted_deduplicated_sorted_samtools_mpileup_depth.tsv",
  "/home/bgilot/Downloads/21014a013_01_bowtie2_unblasted_sorted_corrected_sorted_deduplicated_sorted_samtools_mpileup_depth.tsv",
  "/home/bgilot/Downloads/21014a014_01_bowtie2_unblasted_sorted_corrected_sorted_deduplicated_sorted_samtools_mpileup_depth.tsv",
  "/home/bgilot/Downloads/21014a015_01_bowtie2_unblasted_sorted_corrected_sorted_deduplicated_sorted_samtools_mpileup_depth.tsv",
  "/home/bgilot/Downloads/21014a016_01_bowtie2_unblasted_sorted_corrected_sorted_deduplicated_sorted_samtools_mpileup_depth.tsv",
  "/home/bgilot/Downloads/21014a017_01_bowtie2_unblasted_sorted_corrected_sorted_deduplicated_sorted_samtools_mpileup_depth.tsv",
  "/home/bgilot/Downloads/21014a018_01_bowtie2_unblasted_sorted_corrected_sorted_deduplicated_sorted_samtools_mpileup_depth.tsv",
  "/home/bgilot/Downloads/21014a019_01_bowtie2_unblasted_sorted_corrected_sorted_deduplicated_sorted_samtools_mpileup_depth.tsv",
  "/home/bgilot/Downloads/21014a020_01_bowtie2_unblasted_sorted_corrected_sorted_deduplicated_sorted_samtools_mpileup_depth.tsv",
  "/home/bgilot/Downloads/21014a021_01_bowtie2_unblasted_sorted_corrected_sorted_deduplicated_sorted_samtools_mpileup_depth.tsv",
  "/home/bgilot/Downloads/21014a022_01_bowtie2_unblasted_sorted_corrected_sorted_deduplicated_sorted_samtools_mpileup_depth.tsv",
  "/home/bgilot/Downloads/21014a023_01_bowtie2_unblasted_sorted_corrected_sorted_deduplicated_sorted_samtools_mpileup_depth.tsv",
  "/home/bgilot/Downloads/21014a024_01_bowtie2_unblasted_sorted_corrected_sorted_deduplicated_sorted_samtools_mpileup_depth.tsv",
  "/home/bgilot/Downloads/21014a025_01_bowtie2_unblasted_sorted_corrected_sorted_deduplicated_sorted_samtools_mpileup_depth.tsv")


depth_20076Pos_01  <- read.csv("/home/bgilot/Downloads/20076Pos_01_bowtie2_unblasted_sorted_corrected_sorted_deduplicated_sorted_samtools_mpileup_depth.tsv", sep = "\t")
depth_20076a138_01 <- read.csv("/home/bgilot/Downloads/20076a138_01_bowtie2_unblasted_sorted_corrected_sorted_deduplicated_sorted_samtools_mpileup_depth.tsv", sep = "\t")
depth_20076a139_01 <- read.csv("/home/bgilot/Downloads/20076a139_01_bowtie2_unblasted_sorted_corrected_sorted_deduplicated_sorted_samtools_mpileup_depth.tsv", sep = "\t")
depth_20076a140_01 <- read.csv("/home/bgilot/Downloads/20076a140_01_bowtie2_unblasted_sorted_corrected_sorted_deduplicated_sorted_samtools_mpileup_depth.tsv", sep = "\t")
depth_20076a141_01 <- read.csv("/home/bgilot/Downloads/20076a141_01_bowtie2_unblasted_sorted_corrected_sorted_deduplicated_sorted_samtools_mpileup_depth.tsv", sep = "\t")
depth_20076a142_01 <- read.csv("/home/bgilot/Downloads/20076a142_01_bowtie2_unblasted_sorted_corrected_sorted_deduplicated_sorted_samtools_mpileup_depth.tsv", sep = "\t")
depth_20076a143_01 <- read.csv("/home/bgilot/Downloads/20076a143_01_bowtie2_unblasted_sorted_corrected_sorted_deduplicated_sorted_samtools_mpileup_depth.tsv", sep = "\t")
depth_20076a144_01 <- read.csv("/home/bgilot/Downloads/20076a144_01_bowtie2_unblasted_sorted_corrected_sorted_deduplicated_sorted_samtools_mpileup_depth.tsv", sep = "\t")
depth_20076a145_01 <- read.csv("/home/bgilot/Downloads/20076a145_01_bowtie2_unblasted_sorted_corrected_sorted_deduplicated_sorted_samtools_mpileup_depth.tsv", sep = "\t")
depth_20076a146_01 <- read.csv("/home/bgilot/Downloads/20076a146_01_bowtie2_unblasted_sorted_corrected_sorted_deduplicated_sorted_samtools_mpileup_depth.tsv", sep = "\t")
depth_21014a011_01 <- read.csv("/home/bgilot/Downloads/21014a011_01_bowtie2_unblasted_sorted_corrected_sorted_deduplicated_sorted_samtools_mpileup_depth.tsv", sep = "\t")
depth_21014a012_01 <- read.csv("/home/bgilot/Downloads/21014a012_01_bowtie2_unblasted_sorted_corrected_sorted_deduplicated_sorted_samtools_mpileup_depth.tsv", sep = "\t")
depth_21014a013_01 <- read.csv("/home/bgilot/Downloads/21014a013_01_bowtie2_unblasted_sorted_corrected_sorted_deduplicated_sorted_samtools_mpileup_depth.tsv", sep = "\t")
depth_21014a014_01 <- read.csv("/home/bgilot/Downloads/21014a014_01_bowtie2_unblasted_sorted_corrected_sorted_deduplicated_sorted_samtools_mpileup_depth.tsv", sep = "\t")
depth_21014a015_01 <- read.csv("/home/bgilot/Downloads/21014a015_01_bowtie2_unblasted_sorted_corrected_sorted_deduplicated_sorted_samtools_mpileup_depth.tsv", sep = "\t")
depth_21014a016_01 <- read.csv("/home/bgilot/Downloads/21014a016_01_bowtie2_unblasted_sorted_corrected_sorted_deduplicated_sorted_samtools_mpileup_depth.tsv", sep = "\t")
depth_21014a017_01 <- read.csv("/home/bgilot/Downloads/21014a017_01_bowtie2_unblasted_sorted_corrected_sorted_deduplicated_sorted_samtools_mpileup_depth.tsv", sep = "\t")
depth_21014a018_01 <- read.csv("/home/bgilot/Downloads/21014a018_01_bowtie2_unblasted_sorted_corrected_sorted_deduplicated_sorted_samtools_mpileup_depth.tsv", sep = "\t")
depth_21014a019_01 <- read.csv("/home/bgilot/Downloads/21014a019_01_bowtie2_unblasted_sorted_corrected_sorted_deduplicated_sorted_samtools_mpileup_depth.tsv", sep = "\t")
depth_21014a020_01 <- read.csv("/home/bgilot/Downloads/21014a020_01_bowtie2_unblasted_sorted_corrected_sorted_deduplicated_sorted_samtools_mpileup_depth.tsv", sep = "\t")
depth_21014a021_01 <- read.csv("/home/bgilot/Downloads/21014a021_01_bowtie2_unblasted_sorted_corrected_sorted_deduplicated_sorted_samtools_mpileup_depth.tsv", sep = "\t")
depth_21014a022_01 <- read.csv("/home/bgilot/Downloads/21014a022_01_bowtie2_unblasted_sorted_corrected_sorted_deduplicated_sorted_samtools_mpileup_depth.tsv", sep = "\t")
depth_21014a023_01 <- read.csv("/home/bgilot/Downloads/21014a023_01_bowtie2_unblasted_sorted_corrected_sorted_deduplicated_sorted_samtools_mpileup_depth.tsv", sep = "\t")
depth_21014a024_01 <- read.csv("/home/bgilot/Downloads/21014a024_01_bowtie2_unblasted_sorted_corrected_sorted_deduplicated_sorted_samtools_mpileup_depth.tsv", sep = "\t")
depth_21014a025_01 <- read.csv("/home/bgilot/Downloads/21014a025_01_bowtie2_unblasted_sorted_corrected_sorted_deduplicated_sorted_samtools_mpileup_depth.tsv", sep = "\t")

depth_big <- inner_join(depth_20076Pos_01, depth_20076a138_01, by = c("position" = "position"))
depth_big <- inner_join(depth_big, depth_20076a139_01, by = c("position" = "position"))
depth_big <- inner_join(depth_big, depth_20076a140_01, by = c("position" = "position"))
depth_big <- inner_join(depth_big, depth_20076a141_01, by = c("position" = "position"))
depth_big <- inner_join(depth_big, depth_20076a142_01, by = c("position" = "position"))
depth_big <- inner_join(depth_big, depth_20076a143_01, by = c("position" = "position"))
depth_big <- inner_join(depth_big, depth_20076a144_01, by = c("position" = "position"))
depth_big <- inner_join(depth_big, depth_20076a145_01, by = c("position" = "position"))
depth_big <- inner_join(depth_big, depth_20076a146_01, by = c("position" = "position"))
depth_big <- inner_join(depth_big, depth_21014a011_01, by = c("position" = "position"))
depth_big <- inner_join(depth_big, depth_21014a012_01, by = c("position" = "position"))
depth_big <- inner_join(depth_big, depth_21014a013_01, by = c("position" = "position"))
depth_big <- inner_join(depth_big, depth_21014a014_01, by = c("position" = "position"))
depth_big <- inner_join(depth_big, depth_21014a015_01, by = c("position" = "position"))
depth_big <- inner_join(depth_big, depth_21014a016_01, by = c("position" = "position"))
depth_big <- inner_join(depth_big, depth_21014a017_01, by = c("position" = "position"))
depth_big <- inner_join(depth_big, depth_21014a018_01, by = c("position" = "position"))
depth_big <- inner_join(depth_big, depth_21014a019_01, by = c("position" = "position"))
depth_big <- inner_join(depth_big, depth_21014a020_01, by = c("position" = "position"))
depth_big <- inner_join(depth_big, depth_21014a021_01, by = c("position" = "position"))
depth_big <- inner_join(depth_big, depth_21014a022_01, by = c("position" = "position"))
depth_big <- inner_join(depth_big, depth_21014a023_01, by = c("position" = "position"))
depth_big <- inner_join(depth_big, depth_21014a024_01, by = c("position" = "position"))
depth_big <- inner_join(depth_big, depth_21014a025_01, by = c("position" = "position"))

depth_big %>% clean_names()
write.csv(x=depth_big, file="/home/bgilot/Downloads/depth_big.csv")
summary(depth_big[-1])
summary(depth_big)


sample(colnames(depth_big[-1]), 1)

ncol(depth_big)

for (sample in 2:ncol(depth_big)){
  
  minitable <- depth_big[,c(1,sample)]
  sample_name <- colnames(minitable)[2] %>% str_sub(2, )
  sample_min <- min(minitable[,2])
  sample_1quart <- as.integer(summary(minitable[,2])[2])
  sample_3quart <- as.integer(summary(minitable[,2])[5])
  sample_median <- as.integer(median(minitable[,2]))
  sample_max <- max(minitable[,2])
  
  print(paste(sample_name, "1st Quartile:", sample_1quart, "Median Depth:", sample_median, "3rd Quartile:",  sample_3quart))
  sample_title <- paste(sample_name, "::", "Median Read Depth =", sample_median)
  
  plot_last <- ggplot(minitable, aes(x = minitable[,1], y = minitable[,2])) + geom_point(size = 0.05) + 
    geom_hline(yintercept = sample_1quart, size=0.25, alpha=0.5, color = "red") +
    geom_hline(yintercept = sample_median, size=0.25, alpha=0.5, color = "black") +
    geom_hline(yintercept = sample_3quart, size=0.25, alpha=0.5, color = "green") +
    labs(x = "Genomic Position", y = "Observed Depth",
         title = sample_title,
         caption = "unpublished data") + scale_x_continuous(breaks = seq(0, 30000, by = 2500)) + theme_light()
  
  # plot_filename <- paste0("~/Downloads/", sample_name, ".png")
  ggsave(
    paste0(sample_name, ".png"),
    plot = plot_last,
    device = "png",
    path = paste0("~/Downloads/"),
    scale = 1,
    dpi = "retina",
    limitsize = TRUE
  )
  
  ggsave(
    paste0(sample_name, ".pdf"),
    plot = plot_last,
    device = "pdf",
    path = paste0("~/Downloads/"),
    scale = 1,
    dpi = "retina",
    limitsize = TRUE
  )
}




ggplot(depth_big, aes(x = position, y = X21014a017_01)) + geom_point(size = 0.05) +
  labs(x = "Genomic Position", y = "Observed Depth",
       title = summary(depth_21014a017_01[,-1]),
       caption = "unpublished data",
       color = "Genomic Position") + scale_x_continuous(breaks = seq(0, 30000, by = 2500)) + theme_light()


depth_20076Pos_01  
summary(depth_20076Pos_01[-1])





blasted_25 <- read_delim("Downloads/blasted_25.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)

unblasted_25 <- read_delim("Downloads/unblasted_25.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)

str(blasted_25)

blasted_25[10,7:12]

blasted_25$gini <- as.double(NA)
  
  Gini(c(blasted_25[10,7:12]$acount, blasted_25$ccount, blasted_25$gcount, blasted_25$tcount, blasted_25$ncount, blasted_25$indelcount), unbiased=FALSE)

nrow(blasted_25)

blasted_25 <- read_delim("Downloads/blasted_25.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)

blasted_25$fraction_reference <- blasted_25$refcount/blasted_25$depth * 100 %>% round(3)
blasted_25$fraction_alternate <- blasted_25$altcount/blasted_25$depth * 100 %>% round(3)
blasted_25$fraction_a <- blasted_25$acount/blasted_25$depth * 100 %>% round(3)
blasted_25$fraction_c <- blasted_25$ccount/blasted_25$depth * 100 %>% round(3)
blasted_25$fraction_g <- blasted_25$gcount/blasted_25$depth * 100 %>% round(3)
blasted_25$fraction_t <- blasted_25$tcount/blasted_25$depth * 100 %>% round(3)
blasted_25$fraction_indel <- blasted_25$indelcount/blasted_25$depth * 100 %>% round(3)
blasted_25$fraction_sum <- blasted_25$fraction_a + blasted_25$fraction_c + blasted_25$fraction_g + blasted_25$fraction_t # + blasted_25$fraction_indel

blasted_25$gini <- as.double(NA)
for (i_r in 1:nrow(blasted_25)) {
vec_gini <- c(blasted_25[i_r,]$acount, blasted_25[i_r,]$ccount, blasted_25[i_r,]$gcount, blasted_25[i_r,]$tcount, blasted_25[i_r,]$ncount, blasted_25[i_r,]$indelcount)
  blasted_25[i_r,]$gini <- Gini(vec_gini, unbiased = FALSE)
}

unblasted_25 <- read_delim("Downloads/unblasted_25.tsv", "\t", escape_double = FALSE, trim_ws = TRUE)

unblasted_25$fraction_reference <- unblasted_25$refcount/unblasted_25$depth * 100 %>% round(3)
unblasted_25$fraction_alternate <- unblasted_25$altcount/unblasted_25$depth * 100 %>% round(3)
unblasted_25$fraction_a <- unblasted_25$acount/unblasted_25$depth * 100 %>% round(3)
unblasted_25$fraction_c <- unblasted_25$ccount/unblasted_25$depth * 100 %>% round(3)
unblasted_25$fraction_g <- unblasted_25$gcount/unblasted_25$depth * 100 %>% round(3)
unblasted_25$fraction_t <- unblasted_25$tcount/unblasted_25$depth * 100 %>% round(3)
unblasted_25$fraction_indel <- unblasted_25$indelcount/unblasted_25$depth * 100 %>% round(3)
unblasted_25$fraction_sum <- unblasted_25$fraction_a + unblasted_25$fraction_c + unblasted_25$fraction_g + unblasted_25$fraction_t # + unblasted_25$fraction_indel

unblasted_25$gini <- as.double(NA)
for (i_r in 1:nrow(unblasted_25)) {
  vec_gini <- c(unblasted_25[i_r,]$acount, unblasted_25[i_r,]$ccount, unblasted_25[i_r,]$gcount, unblasted_25[i_r,]$tcount, unblasted_25[i_r,]$ncount, unblasted_25[i_r,]$indelcount)
  unblasted_25[i_r,]$gini <- Gini(vec_gini, unbiased = FALSE)
}


ggplot() + geom_line(aes(x=blasted_25$pos ,y=blasted_25$gini,color='red'))  + scale_x_continuous(breaks = seq(0, 30000, by = 2500)) + theme_light()
ggplot() + geom_line(aes(x=unblasted_25$pos ,y=unblasted_25$gini,color='blue'))  + scale_x_continuous(breaks = seq(0, 30000, by = 2500)) + theme_light()
ggplot() + geom_line(aes(x=blasted_25$pos ,y=blasted_25$gini,color='red')) + geom_line(aes(x=unblasted_25$pos ,y=unblasted_25$gini,color='blue'))  + scale_x_continuous(breaks = seq(0, 30000, by = 2500)) + theme_light()

ggplot(blasted_25, aes(x = gini)) + geom_density() + theme_light()
ggplot(unblasted_25, aes(x = gini)) + geom_density() + theme_light()


ggplot() + geom_line(aes(x=blasted_25$depth ,y=blasted_25$gini,color='red'))  + scale_x_continuous(breaks = seq(0, 30000, by = 2500)) + theme_light()


colnames(unblasted_25)[colnames(unblasted_25) == 'gini'] <- 'inequity'
write.csv(x=unblasted_25, file="/home/bgilot/Downloads/unblasted_25_fractions.csv")

colnames(blasted_25)[colnames(blasted_25) == 'gini'] <- 'inequity'
write.csv(x=blasted_25, file="/home/bgilot/Downloads/blasted_25_fractions.csv")

unblasted_25$depth %>% summary()
blasted_25$depth %>% summary()
# 
# +   stat_summary(fun.data=mean_cl_normal) + geom_smooth(method='lm') + 
#   labs(x = "% Admixture - [MT007544.1]", y = "% Observed - [MT007544.1]",
#        title = "Twist Synthetic SARS-CoV-2 RNA Control 1 (MT007544.1)",
#        subtitle = "Detection of [19065 :: T > C], [22303 :: T > G], [26144 :: G > T]",
#        caption = "unpublished data",
#        color = "Genomic Position") + facet_wrap(~ colour) 


# 
# # + geom_text() + annotate("text",label = c("19065:T->C", "22303:T->G", "26144:G->T"))
# 
#   theme(legend.position = c(.5, .97),
#         legend.background = element_rect(fill = "transparent")) +
#   guides(color = guide_legend(direction = "horizontal")) + facet_wrap(~ colour) + scale_fill_discrete(name = "Dose", labels = c("A", "B", "C"))
# 
# 
# 
# fixed_data_plot <- fixed_data %>% filter(position_position != 29754)
# ggplot(fixed_data_plot, aes(percent_c1, fraction_Reference, color = fixed_data_plot$colour))  + geom_point(color = fixed_data_plot$colour, shape = fixed_data_plot$colour, size = 5) +   stat_summary(fun.data=mean_cl_normal) + geom_smooth(method='lm') + theme(axis.title.x = element_text(vjust = 0, size = 15),axis.title.y = element_text(vjust = 2, size = 15)) +
#   labs(x = "% Synthetic Mixture", y = "% Variant Detected",
#        title = "Controls: Synthetic Mixtures",
#        subtitle = "Detectability of C1 from 0% to 99.9%",
#        caption = "unpublished data",
#        tag = "Fig. 1", color = "Genome Position") +
#   theme(legend.position = c(.5, .97),
#         legend.background = element_rect(fill = "transparent")) +
#   guides(color = guide_legend(direction = "horizontal"))
# 
# fixed_data_plot <- fixed_data %>% filter(percent_c1 >= 0) %>% filter(percent_c1 <= 1) %>% filter(position_position != 29754)
# ggplot(fixed_data_plot, aes(percent_c1, fraction_Reference, color = fixed_data_plot$colour))  + geom_point(color = fixed_data_plot$colour, shape = fixed_data_plot$colour, size = 5) +   stat_summary(fun.data=mean_cl_normal) + geom_smooth(method='lm') + theme(axis.title.x = element_text(vjust = 0, size = 15),axis.title.y = element_text(vjust = 2, size = 15)) +
#   labs(x = "% Synthetic Mixture", y = "% Variant Detected",
#        title = "Controls: Synthetic Mixtures",
#        subtitle = "Detectability of C1 from 0% to 1%",
#        caption = "unpublished data",
#        tag = "Fig. 2", color = "Genome Position") +
#   theme(legend.position = c(.5, .97),
#         legend.background = element_rect(fill = "transparent")) +
#   guides(color = guide_legend(direction = "horizontal"))
# 
# 
# fixed_data_plot <- fixed_data %>% filter(percent_c1 >= 95) %>% filter(percent_c1 <= 100) %>% filter(position_position != 29754)
# ggplot(fixed_data_plot, aes(percent_c1, fraction_Reference))  + geom_point(color = fixed_data_plot$colour, shape = fixed_data_plot$colour, size = 5) +   stat_summary(fun.data=mean_cl_normal) + geom_smooth(method='lm') + theme(axis.title.x = element_text(vjust = 0, size = 15),axis.title.y = element_text(vjust = 2, size = 15)) +
#   labs(x = "% Synthetic Mixture", y = "% Reference Observed",
#        title = "Controls: Synthetic Mixtures",
#        subtitle = "Detectability of C1 from 95% to 100%",
#        caption = "unpublished data",
#        tag = "Fig. 2", color = "Genome Position") +
#   theme(legend.position = c(.5, .97),
#         legend.background = element_rect(fill = "transparent")) +
#   guides(color = guide_legend(direction = "horizontal")) + facet_wrap(~ colour)
# 
# 
# 
# fixed_data_plot <- fixed_data %>% filter(percent_c1 >= 99) %>% filter(percent_c1 <= 100) %>% filter(position_position != 29754)
# ggplot(fixed_data_plot, aes(percent_c1, fraction_Reference, color = colour))  + geom_point(color = fixed_data_plot$colour, shape = fixed_data_plot$colour, size = 5) +   stat_summary(fun.data=mean_cl_normal) + geom_smooth(method='lm') + theme(axis.title.x = element_text(vjust = 0, size = 15),axis.title.y = element_text(vjust = 2, size = 15)) +
#   labs(x = "% Synthetic Mixture", y = "% Variant Detected",
#        title = "Controls: Synthetic Mixtures",
#        subtitle = "Detectability of C1 from 99% to 100%",
#        caption = "unpublished data",
#        tag = "Fig. 3") +
#   theme(legend.position = c(.85, .85),
#         legend.background = element_rect(fill = "transparent")) 
# # +
# # labs(x = "% Synthetic Mixture", y = "% Variant Detected") +
# # ggtitle("Controls: Synthetic Mixtures")
# 
# fixed_data_low <-  fixed_data %>% filter(percent_c1 >= 0) %>% filter(percent_c1 <= 5)
# ggplot(fixed_data_low, aes(percent_c1, fraction_Reference), colour = colour)  + geom_point(size = 1) +   stat_summary(fun.data=mean_cl_normal) +   geom_smooth(method='lm')
# 
# fixed_data <- read_csv("Downloads/to_fix.csv")  %>% filter(corrected_count >= 10000) %>% filter(percent_c1 >= 98) %>% filter(percent_c1 <= 100)
# 
# p <- ggplot(fixed_data, aes(percent_c1, fraction_Reference), colour = position_position)  + geom_point(size = 1) +   stat_summary(fun.data=mean_cl_normal) +   geom_smooth(method='lm', formula= y~x)
# 
# colnames(joined_controls)
# str(joined_controls)