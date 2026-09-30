## ------------------------------------------------------------
## Differential APA Analysis (REDu / REDi)
## Original date: 2025-02-28 by Shan Yu
## Objective: quantify 3′UTR APA (distal vs proximal) and intronic APA
## Inputs: PAS expression matrix (RDS) with counts, RPM, and usage columns
## Notes: This is a re-factored script that improves readability and safety 
##        checks without changing the analysis logic or thresholds.
##        Re-factored on 2025-09-05 by Celeste Nobrega
## ------------------------------------------------------------
# REDu measures the relative expression levels between the top two most abundant APA isoforms in the 3'-most exon only.
# REDi measures the relative expression levels between the top expressed isoform in the 3'-most exon and the top expressed isoform in an intron.
# DAPA: fisher.test

library(reshape2)
library(tidyverse)
library(ggrepel)
library(openxlsx)
library(dplyr)


## =========================
## 0) USER CONFIG
## =========================
# Statistical thresholds
p_cutoff    <- 0.05           # p-value threshold (nominal & BH-adj)
fc_cutoff   <- 1.2            # fold-change threshold in linear space
pseudoCount <- 0              # added to RPMs to avoid log2(0)

# 3′UTR analysis cutoff (minimum average usage)
UPAS_cutoff <- 0

# Intronic analysis cutoffs (minimum average usage)
UPAS_cutoff_IPA <- 10
IPAS_cutoff     <- 5

# File system paths
mainpath="/Volumes/biocore01/scratch/cnobrega/3_REAP/APA_analysis/MAAPER"
inpath   <- "results"      # input directory (relative to mainpath)
outdir   <- "results"      # output directory (relative to mainpath)

# Sample naming
samples <- c(
  "Control_3day_rep1", "Control_3day_rep2", "Control_3day_rep3",
  "Control_5day_rep1", "Control_5day_rep2", "Control_5day_rep3",
  "Treatment_3day_rep1", "Treatment_3day_rep2", "Treatment_3day_rep3",
  "Treatment_5day_rep1", "Treatment_5day_rep2", "Treatment_5day_rep3"
)

## =========================
## 1) LOAD DATA & PREPARE
## =========================
setwd(mainpath)
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# Load PAS table (RDS) expected to contain at least these columns:
# chromosome, start, end, Position, Peak_ID, PAS_ID, gene_id, strand,
# region, and per-sample *_count, *_rpm, *_usage columns.
PAS_aaa1 <- readRDS(file.path(inpath, "pas_exp.rds"))
backup   <- PAS_aaa1  # keep an untouched copy; we reset to this before branches

# --- Defensive checks ---
required_base_cols <- c("chromosome","start","end","Position","Peak_ID","PAS_ID",
                        "gene_id","strand","region")
missing_base <- setdiff(required_base_cols, names(PAS_aaa1))
if (length(missing_base) > 0) {
  stop(sprintf("Missing required columns: %s", paste(missing_base, collapse = ", "))) }

# Verify sample-derived columns exist
affix_exists <- function(base, suffix) paste0(base, suffix) %in% names(PAS_aaa1)
missing_count <- setdiff(paste0(samples, "_count"), names(PAS_aaa1))
missing_rpm   <- setdiff(paste0(samples, "_rpm"),   names(PAS_aaa1))
missing_use   <- setdiff(paste0(samples, "_usage"), names(PAS_aaa1))
if (length(missing_count) > 0) warning("Missing *_count columns: ", paste(missing_count, collapse=", "))
if (length(missing_rpm)   > 0) warning("Missing *_rpm columns: ", paste(missing_rpm,   collapse=", "))
if (length(missing_use)   > 0) warning(sprintf("Missing *_usage columns:", paste(missing_use, collapse=", ")))

# Coerce *_count and *_rpm to numeric (keep NA if absent)
num_coerce_cols <- intersect(paste0(samples, "_count"), names(PAS_aaa1))
PAS_aaa1[num_coerce_cols] <- lapply(PAS_aaa1[num_coerce_cols], as.numeric)
num_coerce_cols <- intersect(paste0(samples, "_rpm"), names(PAS_aaa1))
PAS_aaa1[num_coerce_cols] <- lapply(PAS_aaa1[num_coerce_cols], as.numeric)

# Calculate average across all *_usage columns
use_cols <- names(PAS_aaa1)[grepl(paste0("_usage$"), names(PAS_aaa1))]
if (length(use_cols) == 0) stop(sprintf("No columns end with '_usage'"))
PAS_aaa1$avUsage <- rowMeans(PAS_aaa1[, use_cols], na.rm = TRUE)

## -------------------------
## 1.5) Aggregate biological replicates by condition
## -------------------------
# Four conditions, each with 3 biological replicates:
#   - 3-day Control:   Control_3day_rep1..3
#   - 3-day Treatment: Treatment_3day_rep1..3
#   - 5-day Control:   Control_5day_rep1..3
#   - 5-day Treatment: Treatment_5day_rep1..3
# For Fisher tests (counts), we sum counts across replicates.
# For RE metrics (RPM), we take the mean RPM across replicates (RPM is already library-size normalized).

# Helper to safely get intersecting columns
get_cols <- function(bases, suffix) intersect(paste0(bases, suffix), names(PAS_aaa1))

ctrl3_bases  <- c("Control_3day_rep1","Control_3day_rep2","Control_3day_rep3")
treat3_bases <- c("Treatment_3day_rep1","Treatment_3day_rep2","Treatment_3day_rep3")
ctrl5_bases  <- c("Control_5day_rep1","Control_5day_rep2","Control_5day_rep3")
treat5_bases <- c("Treatment_5day_rep1","Treatment_5day_rep2","Treatment_5day_rep3")

# Counts: sum across reps
PAS_aaa1$Ctrl3d_count  <- rowSums(PAS_aaa1[ get_cols(ctrl3_bases,  "_count") ], na.rm = TRUE)
PAS_aaa1$Treat3d_count <- rowSums(PAS_aaa1[ get_cols(treat3_bases, "_count") ], na.rm = TRUE)
PAS_aaa1$Ctrl5d_count  <- rowSums(PAS_aaa1[ get_cols(ctrl5_bases,  "_count") ], na.rm = TRUE)
PAS_aaa1$Treat5d_count <- rowSums(PAS_aaa1[ get_cols(treat5_bases, "_count") ], na.rm = TRUE)

# RPM: mean across reps
PAS_aaa1$Ctrl3d_rpm  <- rowMeans(PAS_aaa1[ get_cols(ctrl3_bases,  "_rpm") ], na.rm = TRUE)
PAS_aaa1$Treat3d_rpm <- rowMeans(PAS_aaa1[ get_cols(treat3_bases, "_rpm") ], na.rm = TRUE)
PAS_aaa1$Ctrl5d_rpm  <- rowMeans(PAS_aaa1[ get_cols(ctrl5_bases,  "_rpm") ], na.rm = TRUE)
PAS_aaa1$Treat5d_rpm <- rowMeans(PAS_aaa1[ get_cols(treat5_bases, "_rpm") ], na.rm = TRUE)

# Aggregated sample basenames used downstream
agg_samples <- c("Ctrl3d","Treat3d","Ctrl5d","Treat5d")

# Redefine matched pairs to run *two* comparisons with replicates aggregated:
#   - Treatment: 5d vs 3d  (3d as denominator/reference)
#   - Control:   5d vs 3d  (3d as denominator/reference)
ctrl_samples <- c("Treat3d", "Ctrl3d")   # denominators (reference = 3 day)
test_samples <- c("Treat5d", "Ctrl5d")   # numerators  (5 day)

## =============================================================
## 2) 3′UTR APA (REDu: distal vs proximal within 3′-most exon)
## =============================================================
# Folder encodes parameters for provenance
outpath1 <- file.path(outdir, sprintf("DE-UPA_UPAS-%s%d", "avUsage", UPAS_cutoff))
if (!dir.exists(outpath1)) dir.create(outpath1, recursive = TRUE)

# Keep terminal exon PAS with sufficient evaluate
# region uses terms: "3' most exon" or "Single exon"
PAS_3UTR <- subset(PAS_aaa1, region %in% c("3' most exon", "Single exon") & avUsage >= UPAS_cutoff)

# Within each gene, select top 2 expressed PAS by |avUsage|
PAS_3UTR_2top <- PAS_3UTR %>%
  arrange(desc(abs(avUsage))) %>%
  group_by(gene_id) %>%
  slice(1:2) %>%            # genes with only one PAS can still appear here
  ungroup()

# Keep only genes with >= 2 PAS entries (strictly 2 in this selection)
PAS_3UTR_2top2 <- PAS_3UTR_2top %>%
  group_by(gene_id) %>%
  filter(n() > 1) %>%
  ungroup()

message(sprintf("[3'UTR] genes with 2 top PAS: %d rows", nrow(PAS_3UTR_2top2)))

# --- BED outputs for genome browser visualization ---
temp <- as.data.frame(PAS_3UTR_2top2)
write.table(temp[, c("chromosome","start","end","Peak_ID","avUsage","strand")],
            file = file.path(outpath1, "Peak_3UTR_2top2.bed"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)

temp$name <- paste0(temp$gene_id,"_",temp$PAS_ID)
temp2 <- temp[,c("chromosome","Position","Position","name","avUsage","strand")]
temp2$Position.1 <- temp2$Position + 1
write.table(temp2, file.path(outpath1, "PAS_3UTR_2top2.bed"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)

gene_name <- unique(PAS_3UTR_2top2$gene_id) # 2646
message(sprintf("[3'UTR] gene_num: %d", length(gene_name)))

# --- Label distal vs proximal within each gene ---
df <- temp %>%
  group_by(gene_id) %>%
  mutate(
    is_distal = case_when(
      strand == "+" & Position == max(Position) ~ TRUE,
      strand == "-" & Position == min(Position) ~ TRUE,
      TRUE ~ FALSE
    )
  ) %>%
  ungroup()

# --- Compute dPUI and RE per sample (counts-based) ---
samples_count <- paste0(agg_samples, "_count")
resdf <- data.frame(gene_id = gene_name)

for (sample in samples_count) {
  distal_counts <- df %>%
    dplyr::filter(is_distal) %>%
    dplyr::select(gene_id, !!rlang::sym(sample)) %>%
    dplyr::rename(distal_count = !!rlang::sym(sample))
  
  proximal_counts <- df %>%
    dplyr::filter(!is_distal) %>%
    dplyr::select(gene_id, !!rlang::sym(sample)) %>%
    dplyr::rename(proximal_count = !!rlang::sym(sample))
   
  distal_col    <- paste0(sample, "_distal_count")
  proximal_col  <- paste0(sample, "_proximal_count")
  dPUI_col      <- paste0(sample, "_dPUI")
  RE_col        <- paste0(sample, "_RE")
  
  combined <- distal_counts %>%
    inner_join(proximal_counts, by = "gene_id") %>%
    mutate(
      !!distal_col := distal_count,
      !!proximal_col := proximal_count,
      !!dPUI_col := distal_count / (distal_count + proximal_count),
      !!RE_col := distal_count / proximal_count
    ) %>%
    select(gene_id, !!sym(distal_col), !!sym(proximal_col), !!sym(dPUI_col),!!sym(RE_col))
  
  resdf <- resdf %>% left_join(combined, by = "gene_id")
}
# resdf <- as.data.frame(resdf);head(resdf,2)
write.xlsx(resdf,file.path(outpath1, paste0("pas_UPAStop2_dPUI.xlsx")))
write.csv(resdf,file.path(outpath1, paste0("pas_UPAStop2_dPUI.csv")), quote = FALSE, row.names = FALSE)

# --- REDu: per-pair Fisher test & RED metrics on RPMs ---
List.REDu <- list()

for (n in seq_along(test_samples)) {
  ctrl <- ctrl_samples[n]
  test <- test_samples[n]
  pair <- paste0(test, ".vs.", ctrl)
  
  REDu <- data.frame(gene_id = character(),
                     dPASid = character(),
                     pPASid = character(),
                     aUTR_size = numeric(),
                     RE_dPAS = numeric(),
                     RE_pPAS = numeric(),
                     RE_test = numeric(),
                     RE_ctrl = numeric(),
                     REDu1 = numeric(),
                     REDu2 = numeric(),
                     pval.fisher = numeric(),
                     stringsAsFactors = FALSE
  )
  
  for (k in gene_name) {
    data <- subset(PAS_3UTR_2top2, gene_id == k)
    data <- data[order(data$Position, decreasing = FALSE), ]
    
    # 2x2 Fisher table on counts
    contingency_table <- matrix(c(
      data[[paste0(test, "_count")]][2] + pseudoCount, data[[paste0(ctrl, "_count")]][2] + pseudoCount,
      data[[paste0(test, "_count")]][1] + pseudoCount, data[[paste0(ctrl, "_count")]][1] + pseudoCount
    ), nrow = 2)
    pval.fisher <- fisher.test(contingency_table)$p.value
    
    aUTR_size <- abs(as.numeric(data$Position[2]) - as.numeric(data$Position[1]))
    
    if (data$strand[1] == "+") {
      dPASid  <- data$Peak_ID[2]
      pPASid  <- data$Peak_ID[1]
      RE_dPAS <- log2((data[[paste0(test, "_rpm")]][2] + pseudoCount) / (data[[paste0(ctrl, "_rpm")]][2] + pseudoCount))
      RE_pPAS <- log2((data[[paste0(test, "_rpm")]][1] + pseudoCount) / (data[[paste0(ctrl, "_rpm")]][1] + pseudoCount))
      RE_test <- log2((data[[paste0(test, "_rpm")]][2] + pseudoCount) / (data[[paste0(test, "_rpm")]][1] + pseudoCount))
      RE_ctrl <- log2((data[[paste0(ctrl, "_rpm")]][2] + pseudoCount) / (data[[paste0(ctrl, "_rpm")]][1] + pseudoCount))
      REDu1   <- RE_dPAS - RE_pPAS
      REDu2   <- RE_test - RE_ctrl
      REDu[nrow(REDu) + 1, ] <- c(k, dPASid, pPASid, aUTR_size, RE_dPAS, RE_pPAS, RE_test, RE_ctrl, REDu1, REDu2, pval.fisher)
    }
    if (data$strand[1] == "-") {
      dPASid  <- data$Peak_ID[1]
      pPASid  <- data$Peak_ID[2]
      RE_dPAS <- log2((data[[paste0(test, "_rpm")]][1] + pseudoCount) / (data[[paste0(ctrl, "_rpm")]][1] + pseudoCount))
      RE_pPAS <- log2((data[[paste0(test, "_rpm")]][2] + pseudoCount) / (data[[paste0(ctrl, "_rpm")]][2] + pseudoCount))
      RE_test <- log2((data[[paste0(test, "_rpm")]][1] + pseudoCount) / (data[[paste0(test, "_rpm")]][2] + pseudoCount))
      RE_ctrl <- log2((data[[paste0(ctrl, "_rpm")]][1] + pseudoCount) / (data[[paste0(ctrl, "_rpm")]][2] + pseudoCount))
      REDu1   <- RE_dPAS - RE_pPAS
      REDu2   <- RE_test - RE_ctrl
      REDu[nrow(REDu) + 1, ] <- c(k, dPASid, pPASid, aUTR_size, RE_dPAS, RE_pPAS, RE_test, RE_ctrl, REDu1, REDu2, pval.fisher)
    }
  }
  
  # Casting numeric columns; call regulation by thresholds
  REDu_aaa <- REDu
  REDu_aaa[c(4:11)] <- as.numeric(unlist(REDu_aaa[c(4:11)]))
  REDu_aaa$regu <- ifelse(REDu_aaa$pval.fisher < p_cutoff & REDu_aaa$REDu1 > log2(fc_cutoff), "Lengthened",
                          ifelse(REDu_aaa$pval.fisher < p_cutoff & REDu_aaa$REDu1 < (-log2(fc_cutoff)), "Shortened", "NO"))
  
  REDu_aaa$pval.fisher.adj <- p.adjust(REDu_aaa$pval.fisher, method = "BH")
  REDu_aaa$regu.adj <- ifelse(REDu_aaa$pval.fisher.adj < p_cutoff & REDu_aaa$REDu1 > log2(fc_cutoff), "Lengthened",
                              ifelse(REDu_aaa$pval.fisher.adj < p_cutoff & REDu_aaa$REDu1 < (-log2(fc_cutoff)), "Shortened", "NO"))
  List.REDu[[pair]] <- REDu_aaa
}

saveRDS(List.REDu, file = file.path(outpath1, sprintf("List.REDu-pseudoC%d.RData", pseudoCount)))

# Clean INF/-INF and drop NAs, then stack with a Comparison label
List.REDu <- lapply(List.REDu, function(x) { x <- x %>% mutate_if(is.numeric, list(~na_if(., Inf))) %>% mutate_if(is.numeric, list(~na_if(., -Inf))); x })
List.REDu <- lapply(List.REDu, function(x) { x <- na.omit(x); x })
List.REDu2 <- lapply(names(List.REDu), function(name) {
  List.REDu[[name]] %>% mutate(Comparison = name)
})
REDu_aaa <- do.call(rbind, List.REDu2)

# Join dPUI/RE table (resdf) by gene_id and write outputs
REDu_aaa2 <- left_join(REDu_aaa, resdf, by = "gene_id")

write.xlsx(REDu_aaa2, file.path(outpath1, sprintf("APA_3UTR-pseudoC%d.xlsx", pseudoCount)))
write.csv(REDu_aaa2, file.path(outpath1, sprintf("APA_3UTR-pseudoC%d.csv", pseudoCount)), quote = FALSE, row.names = FALSE)
saveRDS(REDu_aaa2,  file.path(outpath1, sprintf("APA_3UTR-pseudoC%d.rds",  pseudoCount)))

# Summary of calls per comparison
sink(file.path(outpath1, sprintf("List.REDu-pseudoC%d.summary.txt", pseudoCount)), split = TRUE)
for (i in seq_along(List.REDu)) {
  cat(paste0(names(List.REDu)[i], ":\n"))
  print(table(List.REDu[[i]]$regu))
  print(table(List.REDu[[i]]$regu.adj))
}
sink()

## =============================================================
## 3) Intronic APA (REDi: terminal 3′UTR PAS vs intronic PAS)
## =============================================================
# Folder encodes parameters for provenance
outpath2 <- file.path(outdir, sprintf("DE-IPA_UPAS-%s%d_IPAS-%s%d", "avUsage", UPAS_cutoff_IPA, "avUsage", IPAS_cutoff))
if (!dir.exists(outpath2)) dir.create(outpath2, recursive = TRUE)

# Select top terminal exon PAS (UPAS) and top intronic PAS (IPAS) per gene
iPAS_3UTR <- subset(PAS_aaa1, region == "3' most exon" & avUsage >= UPAS_cutoff_IPA)
iPAS_3UTR_top <- iPAS_3UTR %>% arrange(desc(abs(avUsage))) %>% group_by(gene_id) %>% slice(1) %>% ungroup()

iPAS_UR <- subset(PAS_aaa1, region == "Intron" & avUsage >= IPAS_cutoff)
iPAS_UR_top <- iPAS_UR %>% arrange(desc(abs(avUsage))) %>% group_by(gene_id) %>% slice(1) %>% ungroup()

PAS_IPA <- rbind(iPAS_3UTR_top, iPAS_UR_top)
PAS_IPA2 <- PAS_IPA %>% group_by(gene_id) %>% filter(n() > 1) %>% ungroup()

# --- BED outputs for genome browser visualization ---
temp <- as.data.frame(PAS_IPA2)
temp$name <- paste0(temp$gene_id, "_", temp$PAS_ID)
write.table(temp[, c("chromosome","start","end","Peak_ID","avUsage","strand")],
            file = file.path(outpath2, "Peak_UPAS_IPAS.bed"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)

temp2 <- temp[, c("chromosome","Position","Position","name","avUsage","strand")]
temp2$Position.1 <- temp2$Position + 1
write.table(temp2, file.path(outpath2, "PAS_UPAS_IPAS.bed"),
            row.names = FALSE, col.names = FALSE, quote = FALSE)

gene_name <- unique(PAS_IPA2$gene_id)
message(sprintf("[IPA] gene_num: %d", length(gene_name)))

# Identify the intronic PAS row per gene
# (we'll mark it to compute iPUI and RE per sample)
df <- temp %>%
  group_by(gene_id) %>%
  mutate(is_IPAS = case_when(region == "Intron" ~ TRUE, TRUE ~ FALSE)) %>%
  ungroup()

# --- Compute iPUI and RE per sample (counts-based) ---
samples_count <- paste0(agg_samples, "_count")
resdf <- data.frame(gene_id = gene_name)

for (sample in samples_count) {
  IPAS_counts <- df %>% filter(is_IPAS) %>% select(gene_id, !!sym(sample)) %>% rename(IPAS_count = !!sym(sample))
  UPAS_counts <- df %>% filter(!is_IPAS) %>% select(gene_id, !!sym(sample)) %>% rename(UPAS_count = !!sym(sample))
  
  IPAS_col <- paste0(sample, "_IPAS_count")
  UPAS_col <- paste0(sample, "_UPAS_count")
  iPUI_col <- paste0(sample, "_iPUI")
  RE_col <- paste0(sample, "_RE")
  
  combined <- UPAS_counts %>%
    inner_join(IPAS_counts, by = "gene_id") %>%
    mutate(
      !!IPAS_col := IPAS_count,
      !!UPAS_col := UPAS_count,
      !!iPUI_col := IPAS_count / (UPAS_count + IPAS_count),
      !!RE_col := IPAS_count / UPAS_count
    ) %>%
    select(gene_id, !!sym(IPAS_col), !!sym(UPAS_col), !!sym(iPUI_col),!!sym(RE_col))
  
  resdf <- resdf %>% left_join(combined, by = "gene_id")
}

write.xlsx(as.data.frame(resdf), file.path(outpath2, "pas_UPAS_IPAS_iPUI.xlsx"))
write.csv(as.data.frame(resdf), file.path(outpath2, "pas_UPAS_IPAS_iPUI.csv"), quote = FALSE, row.names = FALSE)


# --- REDi: per-pair Fisher test & RED metrics on RPMs ---
List.REDi <- list()

for (n in seq_along(test_samples)) {
  ctrl <- ctrl_samples[n]
  test <- test_samples[n]
  pair <- paste0(test, ".vs.", ctrl)
  
  REDi <- data.frame(gene_id = character(), 
                     TPAS_ID = character(), 
                     IPAS_ID = character(), 
                     aUTR_size = numeric(),
                     RE_TPA = numeric(), 
                     RE_IPA = numeric(), 
                     RE_test = numeric(), 
                     RE_ctrl = numeric(),
                     REDi1 = numeric(), 
                     REDi2 = numeric(), 
                     pval.fisher = numeric(),
                     stringsAsFactors = FALSE
  )
  
  
  for (k in gene_name) {
    data <- subset(PAS_IPA2, gene_id == k)
    data <- data[order(data$Position, decreasing = FALSE), ]
    
    # 2x2 Fisher table on counts (row2: intronic vs row1: terminal after order)
    contingency_table <- matrix(c(
      data[[paste0(test, "_count")]][2] + pseudoCount, data[[paste0(ctrl, "_count")]][2] + pseudoCount,
      data[[paste0(test, "_count")]][1] + pseudoCount, data[[paste0(ctrl, "_count")]][1] + pseudoCount
    ), nrow = 2)
    pval.fisher <- fisher.test(contingency_table)$p.value
    
    aUTR_size <- abs(as.numeric(data$Position[2]) - as.numeric(data$Position[1]))
    
    if (data$strand[1] == "+") {
      TPAS_ID <- data$PAS_ID[2]
      IPAS_ID <- data$PAS_ID[1]
      RE_TPA  <- log2((data[[paste0(test, "_rpm")]][2] + pseudoCount) / (data[[paste0(ctrl, "_rpm")]][2] + pseudoCount))
      RE_IPA  <- log2((data[[paste0(test, "_rpm")]][1] + pseudoCount) / (data[[paste0(ctrl, "_rpm")]][1] + pseudoCount))
      RE_test <- log2((data[[paste0(test, "_rpm")]][2] + pseudoCount) / (data[[paste0(test, "_rpm")]][1] + pseudoCount))
      RE_ctrl <- log2((data[[paste0(ctrl, "_rpm")]][2] + pseudoCount) / (data[[paste0(ctrl, "_rpm")]][1] + pseudoCount))
      REDi1   <- RE_TPA - RE_IPA
      REDi2   <- RE_test - RE_ctrl
      REDi[nrow(REDi) + 1, ] <- c(k, TPAS_ID, IPAS_ID, aUTR_size, RE_TPA, RE_IPA, RE_test, RE_ctrl, REDi1, REDi2, pval.fisher)
    }
    if (data$strand[1] == "-") {
      TPAS_ID <- data$PAS_ID[1]
      IPAS_ID <- data$PAS_ID[2]
      RE_TPA  <- log2((data[[paste0(test, "_rpm")]][1] + pseudoCount) / (data[[paste0(ctrl, "_rpm")]][1] + pseudoCount))
      RE_IPA  <- log2((data[[paste0(test, "_rpm")]][2] + pseudoCount) / (data[[paste0(ctrl, "_rpm")]][2] + pseudoCount))
      RE_test <- log2((data[[paste0(test, "_rpm")]][1] + pseudoCount) / (data[[paste0(test, "_rpm")]][2] + pseudoCount))
      RE_ctrl <- log2((data[[paste0(ctrl, "_rpm")]][1] + pseudoCount) / (data[[paste0(ctrl, "_rpm")]][2] + pseudoCount))
      REDi1   <- RE_TPA - RE_IPA
      REDi2   <- RE_test - RE_ctrl
      REDi[nrow(REDi) + 1, ] <- c(k, TPAS_ID, IPAS_ID, aUTR_size, RE_TPA, RE_IPA, RE_test, RE_ctrl, REDi1, REDi2, pval.fisher)
    }
  }
  
  # Casting numeric columns; call regulation by thresholds
  REDi_aaa <- REDi
  REDi_aaa[c(4:11)] <- as.numeric(unlist(REDi_aaa[c(4:11)]))
  REDi_aaa$regu <- ifelse(REDi_aaa$pval.fisher < p_cutoff & REDi_aaa$REDi1 > log2(fc_cutoff), "IPAS_decreased",
                          ifelse(REDi_aaa$pval.fisher < p_cutoff & REDi_aaa$REDi1 < (-log2(fc_cutoff)), "IPAS_increased", "NO"))
  
  REDi_aaa$pval.fisher.adj <- p.adjust(REDi_aaa$pval.fisher, method = "BH")
  REDi_aaa$regu.adj <- ifelse(REDi_aaa$pval.fisher.adj < p_cutoff & REDi_aaa$REDi1 > log2(fc_cutoff), "IPAS_decreased",
                              ifelse(REDi_aaa$pval.fisher.adj < p_cutoff & REDi_aaa$REDi1 < (-log2(fc_cutoff)), "IPAS_increased", "NO"))
  List.REDi[[pair]] <- REDi_aaa
}

saveRDS(List.REDi, file=file.path(outpath2, paste0("List.REDi-pseudoC",pseudoCount,".RData")))

# Clean INF/-INF and drop NAs, then stack with a Comparison label
List.REDi <- lapply(List.REDi, function(x) { x <- x %>% mutate_if(is.numeric, list(~na_if(., Inf))) %>% mutate_if(is.numeric, list(~na_if(., -Inf))); x })
List.REDi <- lapply(List.REDi, function(x) { x <- na.omit(x); x })
List.REDi2 <- lapply(names(List.REDi), function(name) {
  List.REDi[[name]] %>% mutate(Comparison = name)
})
REDi_aaa <- do.call(rbind, List.REDi2)

# Join iPUI/RE table (resdf) by gene_id and write outputs
REDi_aaa2 <- left_join(REDi_aaa, resdf, by = "gene_id")

write.xlsx(REDi_aaa2, file.path(outpath2, paste0("APA_IPA-pseudoC",pseudoCount,".xlsx")))
write.csv(REDi_aaa2, file.path(outpath2, paste0("APA_IPA-pseudoC",pseudoCount,".csv")), quote = FALSE, row.names = FALSE)
saveRDS(REDi_aaa2, file.path(outpath2, paste0("APA_IPA-pseudoC",pseudoCount,".rds")))

# Summary of calls per comparison
sink(file.path(outpath2, sprintf("List.REDi-pseudoC%d.summary.txt", pseudoCount)), split = TRUE)
for (i in seq_along(List.REDi)) {
  cat(paste0(names(List.REDi)[i], ":\n"))
  print(table(List.REDi[[i]]$regu))
  print(table(List.REDi[[i]]$regu.adj))
}
sink()

## =========================
## END OF SCRIPT
## =========================
