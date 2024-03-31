#Script to take in a sgRNA summary from MAGeCK and generate columns for
#sgRNA sequences and only keep LFCs as value columns
#Added a section to determine false-positive rate for deltaLFCs based on nt-nt controls

library(tidyverse)
library(stringr)
library(writexl)
library(plotly)

#read in my sgRNA dictionary and convert to a dataframe with 3 columns for sequences and array name
sg_dict <- read_tsv("../screen_raw_data/sgRNA_Dictionary.txt", col_names = F)
sg_dict <- sg_dict %>%
  separate(col = X1, into = c("sequences", "sgrna"), sep = " >>> ") %>%
  separate(col = sequences, into = c("sgrna_pos1", "sgrna_pos2"), sep = ";")

guide_anno <- read_tsv("../screen_raw_data/guide_anno.tsv") %>%
  separate(col = ...1, into = c("sgrna_pos1", "sgrna_pos2"), sep = ";")

sg_dict <- sg_dict %>%
  left_join(guide_anno, by = c("sgrna_pos1", "sgrna_pos2"))

#read in MAGeCK sgRNA summary
infile <- "T98G_TP1vsTP3"
sgrna_summary <- read_tsv(str_c("../MAGeCK Analyses/", infile, ".sgrna_summary.txt"))

#wrangle MAGeCK file to include sequences and keep only LFCs
sgrna_summary_wrangled <- sgrna_summary %>%
  left_join(sg_dict, by = "sgrna") %>%
  select(sgrna, Gene, sgrna_pos1, sgrna_pos2, spot1, spot2, LFC)

sgrna_summary_wrangled <-  sgrna_summary_wrangled %>%
  mutate(T_gate = str_extract(sgrna_pos1, "TTTT")) %>%
  mutate(T_gate = if_else(is.na(T_gate), "no", "yes"))

#Split summary file into nt-nt data frame and all other arrays
ntnt_controls <- sgrna_summary_wrangled %>%
  filter(Gene == "controlcontrol")

sgrna_summary_wrangled <-  sgrna_summary_wrangled %>%
  filter(Gene != "controlcontrol")

#Calculate mean LFC (with sd and CV) in nt-nt controls
ntnt_mean_LFC <- mean(ntnt_controls$LFC)
ntnt_sd_LFC <- sd(ntnt_controls$LFC)
ntnt_CV_LFC <- ntnt_sd_LFC/ntnt_mean_LFC

#Start section to work with ntnt controls to determine false positive events and rates ----
#Perform LFC correction as in rest of dataset
ntnt_controls_corrected <- ntnt_controls %>%
  mutate(LFC = LFC - ntnt_mean_LFC)

#Get mean LFC for each nt-crRNA in spot 1
ntnt_controls_spot1_phenotypes <- ntnt_controls_corrected %>%
  group_by(spot1) %>%
  summarise(LFC_mean_pos1 = mean(LFC),
            LFC_SD_pos1 = sd(LFC),
            LFC_CV_pos1 = LFC_SD_pos1/LFC_mean_pos1,
            n = n())
#Get mean LFC for each nt-crRNA in spot 2
ntnt_controls_spot2_phenotypes <- ntnt_controls_corrected %>%
  group_by(spot2) %>%
  summarise(LFC_mean_pos2 = mean(LFC),
            LFC_SD_pos2 = sd(LFC),
            LFC_CV_pos2 = LFC_SD_pos2/LFC_mean_pos2,
            n = n())

#Merge dfs with averaged LFCs for each nt-crRNA with original df
#and calculate expected LFC as sum of averages in spot1 and spot 2
#calculate deltaLFC as difference of observed LFC and expected LFC
ntnt_controls_corrected <-  left_join(ntnt_controls_corrected, ntnt_controls_spot1_phenotypes, by = "spot1") %>%
  left_join(ntnt_controls_spot2_phenotypes, by = "spot2") %>%
  mutate(LFC_expected = LFC_mean_pos1 + LFC_mean_pos2,
         deltaLFC = LFC - LFC_expected)

##Generate the false_positive_rates for given cutoffs in deltaLFC ----
cutoffs <- seq(-4,0,0.2)
false_positive_rates <- numeric(length(cutoffs))
false_positive_rates <- lapply(cutoffs, function(num) {
  fpr <- sum(ntnt_controls_corrected$deltaLFC < num) / nrow(ntnt_controls_corrected)
  return(fpr)
})

plot(x = cutoffs, y = false_positive_rates, type = "b",
     ylab = "False Positive Rate", xlab = "deltaLFC cut-off",
     main = infile)

fpr_df <- tibble("dLFC_cutoff" = cutoffs,
                 "False_Positive_Rate" = unlist(false_positive_rates))

write_csv2(fpr_df, str_c("FalsePositive_", infile, ".csv"))

###End Section on ntnt controls false positive events ----

#Generate nt-nt LFC corrected LFCs ----
sgrna_summary_LFC_corrected <- sgrna_summary_wrangled %>%
  mutate(LFC = LFC - ntnt_mean_LFC)


#Filter out observations for single sgRNA paired with control guides (single Phenotypes)
single_phenotypes <- sgrna_summary_LFC_corrected %>%
  mutate(ctrl = case_when(
    endsWith(Gene, "control") ~ "ctrl2",
    startsWith(Gene, "control") ~ "ctrl1"
  )) %>%
  filter(ctrl == "ctrl2" | ctrl == "ctrl1") %>%
  mutate(Gene = str_remove(Gene, "control"))

single_phenotypes_grouped <- single_phenotypes %>%
  group_by(Gene) %>%
  summarise(LFC_mean = mean(LFC),
            LFC_std = sd(LFC),
            LFC_CV = LFC_std/LFC_mean,
            n = n()
  )


#Calculate expected pair phenotypes
#Be aware that the matrix calculated here will have all combinations of genes, even those that are NOT
#included in the library, so dont get confused here! These calculations are simply based on the addition
#of each gene's LFC with another gene as taken from its combination with control guides! A later step (i.e. a join)
#must be used to eliminate those combinations which are not in the library.
expected_pair_matrix <- outer(single_phenotypes_grouped$LFC_mean,
                              single_phenotypes_grouped$LFC_mean,
                              FUN = "+")

colnames(expected_pair_matrix) <- single_phenotypes_grouped$Gene

expected_pair_phenotype <- single_phenotypes_grouped %>%
  select(Gene) %>%
  cbind(expected_pair_matrix) %>%
  pivot_longer(cols = !Gene, names_to = "spot2", values_to = "LFC_expected") %>%
  mutate(Gene = str_c(Gene, spot2))


#Filter out observations for paired phenotypes only (i.e. exclude all control guides)
pair_phenotypes <- sgrna_summary_LFC_corrected %>%
  mutate(ctrl = case_when(
    endsWith(Gene, "control") ~ "ctrl2",
    startsWith(Gene, "control") ~ "ctrl1"
  )) %>%
  filter(is.na(ctrl))


#Merge pair_phenotypes with expected pair phenotype to enable final plot

pairs_observed_expected <- pair_phenotypes %>%
  inner_join(expected_pair_phenotype, by = c("Gene", "spot2")) %>%
  mutate(dLFC = LFC - LFC_expected)

#Generate a dataframe that contains the LFCs of each single sgRNA when paired with a control guide
spot1_singles <- single_phenotypes %>%
  filter(ctrl == "ctrl2") %>%
  select(sgrna_pos1, LFC) %>%
  group_by(sgrna_pos1) %>%
  summarise(spot1_ctrl_LFC = mean(LFC),
            spot1_ctrl_CV = sd(LFC))

spot2_singles <- single_phenotypes %>%
  filter(ctrl == "ctrl1") %>%
  select(sgrna_pos2, LFC) %>%
  group_by(sgrna_pos2) %>%
  summarise(spot2_ctrl_LFC = mean(LFC),
            spot2_ctrl_CV = sd(LFC))

pairs_observed_expected <- pairs_observed_expected %>%
  left_join(spot1_singles, by = "sgrna_pos1") %>%
  left_join(spot2_singles, by = "sgrna_pos2") %>%
  select(-ctrl)


#Start working on false negative calculation ----
threshold <- -1.542
essential_crRNAs_spot2 <- single_phenotypes %>%
  filter(ctrl == "ctrl1" & LFC < threshold) %>%
  select(sgrna_pos2) %>%
  unique()


expected_bufferings <- pairs_observed_expected %>%
  filter(T_gate == "yes") %>%
  filter(sgrna_pos2 %in% unlist(essential_crRNAs_spot2))


buff_cutoffs <- seq(-3,5,0.05)

false_negative_rates <- numeric(length(buff_cutoffs))
false_negative_rates <- lapply(buff_cutoffs, function(num) {
  fnr <- 1 - (sum(expected_bufferings$dLFC > num) / nrow(expected_bufferings))
  return(fnr)
})

plot(x = buff_cutoffs, y = false_negative_rates, type = "b",
     ylab = "False Negative Rate", xlab = "deltaLFC cut-off",
     main = infile)

fnr_df <- tibble("dLFC_cutoff" = buff_cutoffs,
                 "False_Negative_Rate" = unlist(false_negative_rates))

write_xlsx(fnr_df, str_c("FalseNegative_", infile, ".xlsx"))
write_xlsx(expected_bufferings, str_c("expected_bufferings_", infile, ".xlsx"))

#Writing Excel Files ----
write_xlsx(pairs_observed_expected, str_c(infile, "_sgRNA_level_dLFCs_broad.xlsx"))


#Plotting ----
p1 <- pairs_observed_expected %>%
  ggplot(aes(x = LFC_expected, y = LFC)) +
  geom_point(aes(color = T_gate, alpha = 0.5)) +
  geom_abline(slope = 1, intercept = 0)



#ggplotly(p2)

df <- pairs_observed_expected %>%
  select(sgrna, T_gate, Gene, spot1, LFC, LFC_expected)

plot_ly(data = df,
        type = "scatter",
        mode = "markers",
        x = ~LFC_expected,
        y = ~LFC,
        color = ~T_gate,
        text = ~sgrna,
        hoverinfo = "text",
        hovertext = paste("Array: ", df$sgrna,
                          "<br> LFC: ", df$LFC )
        )
