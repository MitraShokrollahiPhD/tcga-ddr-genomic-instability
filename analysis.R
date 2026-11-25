library(cBioPortalData)
library(MultiAssayExperiment)
library(tidyverse)

dir.create("results", showWarnings = FALSE)

brca <- cBioDataPack("brca_tcga_pan_can_atlas_2018", ask = FALSE)
mut  <- brca[["mutations"]]
clin <- as.data.frame(colData(brca))

# via GRangesList to avoid expanding the full mutation matrix 
mut_gr <- unlist(as(mut, "GRangesList"))
mut_df <- data.frame(
  sample = sub("\\..*$", "", names(mut_gr)),
  gene   = sub("^[^.]*\\.", "", names(mut_gr)),
  vclass = as.character(mcols(mut_gr)$Variant_Classification)
)
nonsilent <- c("Missense_Mutation", "Nonsense_Mutation", "Frame_Shift_Del",
               "Frame_Shift_Ins", "In_Frame_Del", "In_Frame_Ins",
               "Splice_Site", "Translation_Start_Site", "Nonstop_Mutation")
lof <- c("Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site")

ddr_genes <- c("BRCA1", "BRCA2", "PALB2", "ATM", "CHEK2", "RAD51C", "RAD51D", "BRIP1")

ddr_samples <- mut_df %>%
  filter(gene %in% ddr_genes, vclass %in% nonsilent) %>%
  pull(sample) %>% unique()
ddr_lof_samples <- mut_df %>%
  filter(gene %in% ddr_genes, vclass %in% lof) %>%
  pull(sample) %>% unique()
tp53_samples <- mut_df %>%
  filter(gene == "TP53", vclass %in% nonsilent) %>%
  pull(sample) %>% unique()

# only tumours with mutation data
info <- data.frame(sample = colnames(mut)) %>%
  left_join(dplyr::select(clin, SAMPLE_ID, SUBTYPE, ANEUPLOIDY_SCORE, TMB_NONSYNONYMOUS),
            by = c("sample" = "SAMPLE_ID")) %>%
  mutate(ddr  = factor(ifelse(sample %in% ddr_samples, "DDR-mutated", "No DDR mutation"),
                       levels = c("No DDR mutation", "DDR-mutated")),
         ddr_lof = factor(ifelse(sample %in% ddr_lof_samples, "DDR LoF", "No DDR LoF"),
                          levels = c("No DDR LoF", "DDR LoF")),
         tp53 = factor(ifelse(sample %in% tp53_samples, "Mutant", "WT"),
                       levels = c("WT", "Mutant")),
         aneuploidy = as.numeric(ANEUPLOIDY_SCORE),
         tmb = as.numeric(TMB_NONSYNONYMOUS)) %>%
  filter(!is.na(SUBTYPE)) %>%
  mutate(SUBTYPE = relevel(factor(SUBTYPE), ref = "BRCA_LumA"),
         log_tmb = log10(tmb))

# Unadjusted comparisons
print(wilcox.test(aneuploidy ~ ddr, data = info, conf.int = TRUE))
print(wilcox.test(tmb ~ ddr, data = info, conf.int = TRUE))
info %>% group_by(ddr) %>%
  summarise(n = n(),
            aneuploidy_median = median(aneuploidy, na.rm = TRUE),
            tmb_median = median(tmb, na.rm = TRUE)) %>% print()

# Adjusted models
m_an  <- lm(aneuploidy ~ ddr + SUBTYPE + tp53, data = info)
m_tmb <- lm(log_tmb ~ ddr + SUBTYPE + tp53, data = info)
res <- bind_rows(
  broom::tidy(m_an,  conf.int = TRUE) %>% mutate(outcome = "aneuploidy"),
  broom::tidy(m_tmb, conf.int = TRUE) %>% mutate(outcome = "log10 TMB")
)
print(res)
write_csv(res, "results/model_results.csv")
print(10^broom::tidy(m_tmb, conf.int = TRUE)[2, c("estimate", "conf.low", "conf.high")])

# hypermutated tumours cutoff (TMB > 10/Mb)
info_nh <- info %>% filter(tmb <= 10)
cat("Non-hypermutated tumours:", nrow(info_nh), "\n")
print(table(info_nh$ddr))
print(broom::tidy(lm(aneuploidy ~ ddr + SUBTYPE + tp53, data = info_nh), conf.int = TRUE))
m_tmb_nh <- lm(log_tmb ~ ddr + SUBTYPE + tp53, data = info_nh)
print(broom::tidy(m_tmb_nh, conf.int = TRUE))
print(10^broom::tidy(m_tmb_nh, conf.int = TRUE)[2, c("estimate", "conf.low", "conf.high")])

print(table(info$ddr_lof))
print(table(info$ddr_lof, info$tmb > 10))
lof_res <- bind_rows(
  broom::tidy(lm(aneuploidy ~ ddr_lof + SUBTYPE + tp53, data = info), conf.int = TRUE) %>% mutate(outcome = "aneuploidy"),
  broom::tidy(lm(log_tmb ~ ddr_lof + SUBTYPE + tp53, data = info), conf.int = TRUE) %>% mutate(outcome = "log10 TMB"),
  broom::tidy(lm(log_tmb ~ ddr_lof + SUBTYPE + tp53, data = filter(info, tmb <= 10)), conf.int = TRUE) %>% mutate(outcome = "log10 TMB, non-hypermutated")
) %>% filter(term == "ddr_lofDDR LoF")
print(lof_res)
write_csv(lof_res, "results/model_results_lof.csv")

p1 <- ggplot(info, aes(ddr, aneuploidy, fill = ddr)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.15, size = 0.8, alpha = 0.5) +
  labs(x = NULL, y = "Aneuploidy score") + theme_bw() + theme(legend.position = "none")
p2 <- ggplot(info, aes(ddr, tmb, fill = ddr)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6) +
  geom_jitter(width = 0.15, size = 0.8, alpha = 0.5) +
  scale_y_log10() +
  labs(x = NULL, y = "TMB (mutations/Mb, log scale)") + theme_bw() + theme(legend.position = "none")
p3 <- p1 + facet_wrap(~ SUBTYPE, nrow = 1) + theme(axis.text.x = element_text(angle = 30, hjust = 1))

ggsave("results/aneuploidy_by_ddr.png", p1, width = 4, height = 4, dpi = 300)
ggsave("results/tmb_by_ddr.png", p2, width = 4, height = 4, dpi = 300)
ggsave("results/aneuploidy_by_ddr_subtype.png", p3, width = 10, height = 4, dpi = 300)