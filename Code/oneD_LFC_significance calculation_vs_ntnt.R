#Script to take in a sgRNA summary from MAGeCK and generate columns for
#sgRNA sequences and only keep LFCs as value columns
#Analysis focused on nt-nt controls and single targeting arrays only
#Assigning p-values comparing each gene to the nt-nt controls using Welch's t-test

library(tidyverse)
library(stringr)
library(writexl)
#library(plotly)

#read in my sgRNA dictionary and convert to a dataframe with 3 columns for sequences and array name
sg_dict <- read_tsv("sgRNA_Dictionary.txt")
sg_dict <- sg_dict %>%
  separate(col = A, into = c("sequences", "sgrna"), sep = " >>> ") %>%
  separate(col = sequences, into = c("sgrna_pos1", "sgrna_pos2"), sep = ";")

guide_anno <- read_tsv("guide_anno.tsv") %>%
  separate(col = A, into = c("sgrna_pos1", "sgrna_pos2"), sep = ";")

sg_dict <- sg_dict %>%
  left_join(guide_anno, by = c("sgrna_pos1", "sgrna_pos2"))

#read in MAGeCK sgRNA summary
infile <- "T98G_TP1vsTP3_"
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

#Generate nt-nt LFC corrected LFCs
sgrna_summary_LFC_corrected <- sgrna_summary_wrangled %>%
  mutate(LFC = LFC - ntnt_mean_LFC)

ntnt_controls_LFC_corrected <- ntnt_controls %>%
  mutate(LFC = LFC - ntnt_mean_LFC)


#Filter out observations for single sgRNA paired with control guides (single Phenotypes)
single_phenotypes <- sgrna_summary_LFC_corrected %>%
  mutate(ctrl = case_when(
    endsWith(Gene, "control") ~ "ctrl2",
    startsWith(Gene, "control") ~ "ctrl1"
  )) %>%
  filter(ctrl == "ctrl2" | ctrl == "ctrl1") %>%
  mutate(Gene = str_remove(Gene, "control"))

#Create Gene summaries and calculate p-value by Welch's t-test comparing all LFCs for one gene
#against all LFCs of the ntnt controls
#Choice for statistical test was made based on
#https://towardsdatascience.com/how-to-compare-two-or-more-distributions-9b06ee4d30bf

single_phenotypes_grouped <- single_phenotypes %>%
  group_by(Gene) %>%
  summarise(LFC_mean = mean(LFC),
            LFC_std = sd(LFC),
            LFC_CV = LFC_std/LFC_mean,
            n = n(),
            p_value = (t.test(LFC, ntnt_controls_LFC_corrected$LFC, var.equal = F))$p.value,
            neg_log10p = -log10(p_value)
  )

#Save Results
write_xlsx(single_phenotypes_grouped, str_c(infile, "oneD_Screen_Analysis_with_pvalue.xlsx"))

#Plot result as Volcano Plot (insoired by RNAseq data anaysis)
single_phenotypes_grouped %>%
  ggplot(aes(x = LFC_mean, y = neg_log10p)) +
  geom_point()

# p <- single_phenotypes %>%
#   filter(Gene == "ALYREF") %>%
#   select(LFC) %>%
#   t.test(ntnt_controls_LFC_corrected$LFC, var.equal = F)
# p$p.value
