#Script to take in a sgRNA summary from MAGeCK and generate columns for
#sgRNA sequences and only keep LFCs as value columns

library(tidyverse)
library(stringr)
library(writexl)
library(plotly)

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
infile <- "PATU_TP1vsTP3"
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
