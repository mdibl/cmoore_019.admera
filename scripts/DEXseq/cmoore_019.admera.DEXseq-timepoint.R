library(DEXSeq)
library(SummarizedExperiment)
library(GenomicRanges)
library(DESeq2)
library(ggplot2)
library(ggrepel)
library(limma)
library(dplyr)
library(tibble)
library(matrixStats)
library(readr)
library(stringr)
setwd("/Volumes/biocore01/scratch/cnobrega/3_REAP")


# ---------- build count matrix  ---------- 
# Read counts
counts_df <- read.csv("./APA_analysis/DEXseq/cluster.all.reads.csv")
samples <- scan("./APA_analysis/DEXseq/sample_name.txt", what=character())

# make sure all the columns needed exist in the counts_df
stopifnot("hit_PAS_ID" %in% names(counts_df))
stopifnot(all(samples %in% names(counts_df)))

# make a counts matrix from the counts_df
# integer matrix: rows = PAS, cols = samples
count_mat <- as.matrix(counts_df[, samples])
rownames(count_mat) <- counts_df$hit_PAS_ID
storage.mode(count_mat) <- "integer"

# ---------- build a PAS feature annotation map ---------- 
# read in PAS annotation
pas_anno <- read.csv("./APA_analysis/DEXseq/DEX_human.PAS.hg38.txt", header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)

# filter to only PAS TYPE == 3'UTR(*)
pas_anno <- subset(pas_anno, grepl("^3'UTR", PAS_type))

# HANDLING ANNOTATION DUPLICATES
# create one single new column: group_ID == gene name 
# to determine gene name: try first Gene Symbol, if NA try Ensembl ID, if NA try RefSeq Gene ID, if NA try FAMTOM ID. If all empty, set NA
# investigate then drop rows with NA group_ID
# check for duplicate PAS IDs with the same group_ID and keep only one (these are illegal duplicates)
# ***if there are duplicate PAS IDs but each has a different group_ID, that is OKAY

blank_to_na <- function(x) {
  x <- str_squish(x)
  x[tolower(x) %in% c("", "na")] <- NA_character_
  x
}
pas_anno <- pas_anno %>%
  mutate(
    across(c(`Gene Symbol`, `Ensemble ID`, `RefSeq Gene ID`, `FAMTOM ID`), blank_to_na),
    groupID = coalesce(`Gene Symbol`, `Ensemble ID`, `RefSeq Gene ID`, `FAMTOM ID`)
  )

# Find PAS_ID values that are in counts_df but not in pas_anno if you were to remove all na entries
pas_anno_dropna <- pas_anno %>%
  filter(!is.na(groupID))
mismatched_pas_ids <- setdiff(counts_df$hit_PAS_ID, pas_anno_dropna$PAS_ID)
# Assign groupID = PAS_ID for those NA rows that you need to keep
pas_anno$groupID[pas_anno$PAS_ID %in% mismatched_pas_ids & is.na(pas_anno$groupID)] <- 
  pas_anno$PAS_ID[pas_anno$PAS_ID %in% mismatched_pas_ids & is.na(pas_anno$groupID)]

# drop NA group_IDs (already handled the needed NA values above, the rest are not needed)
pas_anno <- pas_anno %>%
  filter(!is.na(groupID))

# now check if there are any illegal duplicates (same gene name and same PASS_ID)
illegal_dup_check <- pas_anno %>%
  group_by(PAS_ID, groupID) %>%
  tally() %>%
  filter(n > 1)
# If any exist, keep only the first occurrence of each PAS_ID + groupID pair
if (nrow(illegal_dup_check) > 0) {
  pas_anno <- pas_anno[!duplicated(pas_anno[c("PAS_ID", "groupID")]), ]
}

# ---------- align matrix & annotation table; set groupIDs & featureIDs ---------- 
# align count matrix to PAS annotation such that 
# they only contain the same set of PAS IDs, in the same order, and that rownames(count_mat) == pas_anno$PAS_ID
# then set groupID & featureID vars (for use when creating DEXSeqDataSet object)
common_pas <- intersect(rownames(count_mat), pas_anno$PAS_ID) # length(common_pas) == length(rownames(count_mat)) *ONLY IF USING WHOLE COUNT MATRIX; if filtering based on annot. data (like 3'UTR only), then will be FALSE
count_mat  <- count_mat[common_pas, , drop=FALSE]
pas_anno   <- pas_anno[match(common_pas, pas_anno$PAS_ID), ]
groupID <- pas_anno$groupID
featureID <- pas_anno$PAS_ID

# ---------- optional: add GRanges ---------
# make 1-nt or small windows around the PAS
if (all(c("Chromosome","Position","Strand") %in% names(pas_anno))) {
  gr <- GRanges(seqnames = pas_anno$Chromosome,
                ranges   = IRanges(start = as.integer(pas_anno$Position),
                                   width = 1L),
                strand   = pas_anno$Strand)
} else {
  gr <- NULL
}

# ---------- set metadata & design formula ---------
colData <- read.delim("./APA_analysis/DEXseq/design.txt", row.names = 1, check.names = FALSE)
# Ensure these are factors
colData$condition <- factor(colData$condition, levels = c("Control","Treatment"))
colData$timepoint <- factor(colData$timepoint, levels = c("3day","5day"))

stopifnot(all(colnames(count_mat) == rownames(colData)))  # safety

# ---------- APA ANALYSIS ---------
run_dexseq_condition <- function(cond, count_mat, colData, featureID, groupID, gr=NULL,
                                 min_total=10, min_per_time = 0) {
  # ---- Subset to one condition  ---- 
  message(sprintf("Subsetting to one condition..."))
  sub_samples <- rownames(colData)[colData$condition == cond]
  sub_cd  <- droplevels(colData[sub_samples, , drop = FALSE])
  sub_cnt <- count_mat[, sub_samples, drop = FALSE]
  
  # Safety: sample names must match exactly
  stopifnot(all(colnames(sub_cnt) == rownames(sub_cd)))
  
  # Ensure factors used by the design
  sub_cd$sample    <- factor(rownames(sub_cd))
  sub_cd$timepoint <- droplevels(factor(sub_cd$timepoint, levels = c("3day","5day")))
  
  # Optional: convert any remaining character cols to factors to quiet warnings
  is_char <- vapply(sub_cd, is.character, logical(1))
  sub_cd[is_char] <- lapply(sub_cd[is_char], factor)
  
  # Sanity: both conditions present for this timepoint?
  if (length(levels(sub_cd$timepoint)) < 2 || any(table(sub_cd$timepoint) == 0)) {
    stop(sprintf("Condition '%s' does not contain both 3day and 5day samples.", cond))
  }
  
  # ---- Create DEXSeqDataSet object for subsetted data ----
  message("Creating DEXSeqDataSet object for subsetted data (timepoint effect within condition)...")
  dxd <- DEXSeqDataSet(
    countData     = sub_cnt,
    sampleData    = sub_cd,
    design        = ~ sample + exon + timepoint:exon,
    featureID     = featureID,
    groupID       = groupID,
    featureRanges = gr
  )
  
  # ---- Prefilter on the subset (operate on the counts inside dxd to keep alignment) ---- 
  message("Prefiltering...")
  mat <- counts(dxd)  # nrow(mat) == nrow(dxd)
  keep <- rowSums(mat) >= min_total
  dxd  <- dxd[keep, , drop = FALSE]
  if (min_per_time > 0) {
    mat <- counts(dxd)
    tp  <- colData(dxd)$timepoint
    idx3 <- which(tp == "3day")
    idx5 <- which(tp == "5day")
    keep2 <- rowSums(mat[, idx3, drop = FALSE] > 0) >= min_per_time &
      rowSums(mat[, idx5, drop = FALSE] > 0) >= min_per_time
    dxd <- dxd[keep2, , drop = FALSE]
  }
  
  # Drop genes with <2 PAS after filtering (no within-gene contrast possible)
  tab  <- table(rowData(dxd)$groupID)
  dxd  <- dxd[rowData(dxd)$groupID %in% names(tab[tab >= 2]), ]
  
  # ---- Normalize, estimate dispersion, test ---- 
  message("Normalizing and testing (DEU across timepoints within condition)...")
  dxd <- estimateSizeFactors(dxd)
  dxd <- estimateDispersions(dxd)
  dxd <- testForDEU(dxd)               # DEU w.r.t. timepoint in this subset
  dxr <- DEXSeqResults(dxd)
  # Effect sizes for 5day vs 3day (3day as denominator)
  dxd <- estimateExonFoldChanges(dxd, fitExpToVar = "timepoint", denominator = "3day")
  
  # ---- Build a clean results data.frame ---- 
  message("Building a clean results data frame...")
  res   <- S4Vectors::as.data.frame(dxr)
  rd    <- S4Vectors::mcols(dxd)
  fc_cols <- grep("^log2fold_", colnames(rd), value = TRUE)
  if (length(fc_cols) == 0) {
    stop("No fold-change columns found in rowData(dxd). Check your DEXSeq version / column names.")
  }
  
  # coerce possible Rle → numeric
  to_num <- function(z) if (methods::is(z, "Rle")) as.numeric(z) else as.numeric(z)
  fc_df <- base::as.data.frame(
    lapply(S4Vectors::as.list(rd[, fc_cols, drop = FALSE]), to_num),
    check.names = FALSE, stringsAsFactors = FALSE
  )
  
  # create a stable join key from dxd rownames (DEXSeq’s composite name)
  key <- rownames(dxd)
  res$key <- rownames(res)
  fc_df$key <- key
  
  # keep res order
  res_tp <- merge(res, fc_df, by = "key", all.x = TRUE, sort = FALSE)
  
  # sanity: missing FC rows?
  if (anyNA(res_tp$key)) stop("Join produced NA keys unexpectedly.")
  if (any(!res_tp$key %in% key)) stop("Keys in results not found in dxd rownames.")
  
  # restore rownames and drop helper
  rownames(res_tp) <- res_tp$key
  res_tp$key <- NULL
  
  # Annotate with gene/PAS IDs
  res_tp$gene    <- rowData(dxd)$groupID[match(rownames(res_tp), rownames(dxd))]
  res_tp$feature <- rowData(dxd)$featureID[match(rownames(res_tp), rownames(dxd))]
  
  # ---- Arrange output columns in a stable, readable order ----
  # Identify dynamic columns
  fc_cols <- grep("^log2fold_", colnames(res_tp), value = TRUE)
  count_cols <- grep("^countData\\.", colnames(res_tp), value = TRUE)
  
  # Put 3day before 5day and keep natural rep order
  day3_cols <- grep("^countData.*3day", count_cols, value = TRUE)
  day5_cols <- grep("^countData.*5day", count_cols, value = TRUE)
  
  # add extra annotation columns
  ann <- pas_anno[, c("PAS_ID", "Intron_exon_location", "PAS_type")]
  names(ann) <- c("featureID", "Intron_exon_location", "PAS_type")
  ann <- ann[!duplicated(ann$featureID), , drop = FALSE]
  ann$featureID   <- as.character(ann$featureID)
  res_tp$featureID <- as.character(res_tp$featureID)
  res_tp <- merge(res_tp, ann, by = "featureID", all.x = TRUE, sort = FALSE)
  
  # add APA direction if GRanges provided
  message("Calculating APA direction...")
  res_tp$APA_direction <- NA_character_  # create the column for schema stability
  
  if (!is.null(gr)) {
    if (is.null(names(gr)) || !all(rownames(dxd) %in% names(gr))) {
      names(gr) <- rownames(dxd)
    }
    
    rrn    <- rownames(dxd)
    pos    <- as.integer(GenomicRanges::start(gr[rrn]))
    strand <- as.character(GenomicRanges::strand(gr[rrn]))
    
    meta <- data.frame(
      row    = rrn,
      pos    = pos,
      strand = strand,
      stringsAsFactors = FALSE
    )
    
    res_tp$row <- rownames(res_tp)
    res_tp <- left_join(res_tp, meta, by = "row")
    
    # Compute groupwise extrema with guards against all-NA groups
    res_tp <- res_tp |>
      group_by(gene) |>
      mutate(
        g_max = suppressWarnings(max(pos, na.rm = TRUE)),
        g_min = suppressWarnings(min(pos, na.rm = TRUE)),
        g_max = ifelse(is.finite(g_max), g_max, NA_real_),
        g_min = ifelse(is.finite(g_min), g_min, NA_real_),
        is_distal = case_when(
          strand == "+" & !is.na(pos) & !is.na(g_max) & pos == g_max ~ TRUE,
          strand == "-" & !is.na(pos) & !is.na(g_min) & pos == g_min ~ TRUE,
          TRUE ~ FALSE
        )
      ) |>
      ungroup()
    
    # Choose the first log2fold_ column for direction
    fc_col <- grep("^log2fold_", names(res_tp), value = TRUE)[1]
    
    if (!is.na(fc_col)) {
      res_tp <- mutate(
        res_tp,
        APA_direction = case_when(
          is_distal  & .data[[fc_col]] > 0 ~ "Lengthened",
          is_distal  & .data[[fc_col]] < 0 ~ "Shortened",
          !is_distal & .data[[fc_col]] > 0 ~ "Shortened",
          !is_distal & .data[[fc_col]] < 0 ~ "Lengthened",
          TRUE                              ~ "Ambiguous"
        )
      )
      # mark rows with missing inputs as Ambiguous
      bad <- is.na(res_tp[[fc_col]]) | is.na(res_tp$is_distal)
      res_tp$APA_direction[bad] <- "Ambiguous"
    }
    
    # drop helpers
    res_tp$row <- NULL
    res_tp$g_max <- NULL
    res_tp$g_min <- NULL
  }
                  
  # ---- Define the desired order ---- 
  desired_order <- c(
    "groupID",
    "featureID",
    "gene",
    "feature",
    "Intron_exon_location",   # from pas_anno
    "PAS_type",               # from pas_anno
    "APA_direction",          # with directional annot.
    "exonBaseMean",
    "dispersion",
    "stat",
    "pvalue",
    "padj",
    fc_cols,
    day3_cols,
    day5_cols,
    "genomicData.seqnames",
    "genomicData.start",
    "genomicData.end",
    "genomicData.width",
    "genomicData.strand"
  )
  
  # Keep only columns that exist, in that order
  desired_order <- desired_order[desired_order %in% colnames(res_tp)]
  res_tp <- res_tp[, desired_order, drop = FALSE]
  
  list(dxd = dxd, dxr = dxr, results = res_tp)
}

## Run the above function for each timepoint
run_ctrl <- run_dexseq_condition("Control",   count_mat, colData, featureID, groupID, gr = gr, min_total = 10, min_per_time = 1)
run_trt  <- run_dexseq_condition("Treatment", count_mat, colData, featureID, groupID, gr = gr, min_total = 10, min_per_time = 1)

# ---- save results ----
res_ctrl <- run_ctrl$results    # has padj & log2FC for 5day vs 3day @ control
res_trt  <- run_trt$results   # has padj & log2FC for 5day vs 3day @ treated

# ---- Rename L2FC column for export only (keep internal name unchanged) ----
out_ctrl <- res_ctrl
out_trt <- res_trt
names(out_ctrl) <- sub("^log2fold_5day_3day$", "log2fold_5day_v_3day", names(res_ctrl))
names(out_trt)  <- sub("^log2fold_5day_3day$", "log2fold_5day_v_3day", names(res_trt))

## Export the results
write.csv(out_ctrl, "./APA_analysis/DEXseq/results/dexseq/3utr/with-dir/dexseq.3utr.dir.5-v-3.control.csv", row.names = FALSE, quote = FALSE)
write.csv(out_trt, "./APA_analysis/DEXseq/results/dexseq/3utr/with-dir/dexseq.3utr.dir.5-v-3.treatment.csv", row.names = FALSE, quote = FALSE)

# ---------- PC plots ---------- 
make_pca_expression <- function(dxd, color_cols = c("timepoint","condition"),
                                label_samples = TRUE, Wcols = NULL,
                                outfile = NULL, ntop = 2000) {
  stopifnot(methods::is(dxd, "DEXSeqDataSet"))
  cd_full <- as.data.frame(SummarizedExperiment::colData(dxd))
  mat_full <- SummarizedExperiment::assay(dxd)
  
  # ---- DEXSeq quirk: keep only "exon == this" ----
  if ("exon" %in% names(cd_full)) {
    keep_idx <- which(cd_full$exon == "this")
    if (!length(keep_idx)) stop("colData has 'exon' but no rows == 'this'.")
    cd <- cd_full[keep_idx, , drop = FALSE]
    mat <- mat_full[, keep_idx, drop = FALSE]
  } else {
    cd <- cd_full
    mat <- mat_full
  }
  
  # Build a DESeq2 object with the same counts & size factors (subset if needed)
  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = mat,
    colData   = cd,
    design    = ~ 1
  )
  
  # VST on PAS/bin counts
  vst_mat <- SummarizedExperiment::assay(DESeq2::vst(dds, blind = FALSE))
  
  # Optional: remove nuisance covariates
  if (!is.null(Wcols) && length(Wcols) > 0) {
    stopifnot(all(Wcols %in% colnames(cd)))
    vst_mat <- limma::removeBatchEffect(
      vst_mat,
      covariates = as.matrix(cd[, Wcols, drop = FALSE])
    )
  }
  
  # Top-variable features
  if (!is.null(ntop) && ntop > 0 && nrow(vst_mat) > ntop) {
    rv <- matrixStats::rowVars(vst_mat)
    sel <- order(rv, decreasing = TRUE)[seq_len(ntop)]
    vst_mat <- vst_mat[sel, , drop = FALSE]
  }
  
  # PCA
  p <- prcomp(t(vst_mat), center = TRUE, scale. = TRUE)
  pct <- round(100 * (p$sdev^2 / sum(p$sdev^2))[1:2], 1)
  
  pcs <- as.data.frame(p$x[, 1:2]) |>
    tibble::rownames_to_column(var = "sample_id") |>
    dplyr::left_join(cd |> tibble::rownames_to_column(var = "sample_id"), by = "sample_id")
  
  # Plot
  g <- ggplot2::ggplot(pcs, ggplot2::aes(x = PC1, y = PC2)) +
    {
      if (length(color_cols) >= 1 && color_cols[1] %in% names(pcs)) {
        ggplot2::aes(color = .data[[color_cols[1]]], shape = .data[[color_cols[2]]])
      } else NULL
    } +
    ggplot2::geom_point(size = 2) +
    ggplot2::labs(
      x = paste0("PC1 (", pct[1], "%)"),
      y = paste0("PC2 (", pct[2], "%)"),
      title = "PCA on VST PAS/bin counts (exon == 'this')"
    ) +
    ggplot2::theme_minimal(base_size = 12)
  
  if (label_samples && "replicate" %in% names(pcs)) {
    g <- g + ggrepel::geom_text_repel(ggplot2::aes(label = replicate), size = 3)
  }
  
  if (!is.null(outfile)) {
    ggplot2::ggsave(outfile, g, width = 7, height = 5, dpi = 300)
  }
  
  invisible(list(plot = g, pcs = pcs, prcomp = p))
}
p_expr_ctrl <- make_pca_expression(run_ctrl$dxd,
                                   color_cols = c("timepoint","condition"), 
                                   Wcols = NULL,
                                   outfile = "./APA_analysis/DEXseq/results/plots/PCA_expression_ctrl.png")
p_expr_trt  <- make_pca_expression(run_trt$dxd,
                                   color_cols = c("timepoint","condition"), 
                                   Wcols = NULL,
                                   outfile = "./APA_analysis/DEXseq/results/plots/PCA_expression_trt.png")
# ---------- Gene Level Summary ---------
# gene level q from per-bin tests
gene_q_ctrl <- perGeneQValue(run_ctrl$dxr)
gene_q_trt <- perGeneQValue(run_trt$dxr)

# collapse a result table 
collapse_gene <- function(res_df, gene_q = NULL, padj_cut = 0.05, lfc_cut = 0.5) {
  # helpers: return NA if no finite values
  safe_min  <- function(x) { x <- x[is.finite(x)]; if (length(x)) min(x) else NA_real_ }
  safe_max  <- function(x) { x <- x[is.finite(x)]; if (length(x)) max(x) else NA_real_ }
  safe_mean <- function(x) { x <- x[is.finite(x)]; if (length(x)) mean(x) else NA_real_ }
  
  has_dir <- "APA_direction" %in% names(res_df)
  
  res_df %>%
    mutate(
      padj = ifelse(is.infinite(padj), NA_real_, padj),
      log2fold_5day_3day = ifelse(is.infinite(log2fold_5day_3day), NA_real_, log2fold_5day_3day)
    ) %>%
    mutate(
      sig = !is.na(padj) & padj <= padj_cut,
      big = sig & abs(log2fold_5day_3day) >= lfc_cut
    ) %>%
    group_by(groupID) %>%
    summarise(
      n_PAS = n(),
      n_sig = sum(sig, na.rm = TRUE),
      n_big = sum(big, na.rm = TRUE),
      min_padj = safe_min(padj),
      max_absL2FC = safe_max(abs(log2fold_5day_3day[sig])),
      mean_absL2FC = safe_mean(abs(log2fold_5day_3day[sig])),
      dir_consensus = if (has_dir) {
        d <- APA_direction[sig & !is.na(APA_direction)]
        if (!length(d)) NA_character_
        else if (all(d == "Lengthened")) "Lengthened_only"
        else if (all(d == "Shortened"))  "Shortened_only"
        else "Mixed"
      } else NA_character_,
      .groups = "drop"
    ) %>%
    mutate(
      perGeneQ = if (!is.null(gene_q)) unname(gene_q[groupID]) else NA_real_
    ) %>%
    # sort with non‑NA perGeneQ first
    arrange(is.na(perGeneQ), perGeneQ, min_padj, desc(max_absL2FC))
}

gene_summary_ctrl <- collapse_gene(res_ctrl, gene_q_ctrl)
gene_summary_trt <- collapse_gene(res_trt, gene_q_trt)

## Export the tables
write.csv(gene_summary_ctrl, "./APA_analysis/DEXseq/results/summaries/gene/gene-summary.3utr.5-v-3.control.csv", row.names = FALSE, quote=FALSE)
write.csv(gene_summary_trt, "./APA_analysis/DEXseq/results/summaries/gene/gene-summary.3utr.5-v-3.treatment.csv", row.names = FALSE, quote=FALSE)

# ---------- Average Usage ---------
# per condition: the mean usage across the condition's replicates
# per sample: usage = counts_PAS / sum(counts_allPAS_in_gene)
usage_by_condition <- function(res_df, colData) {
  # Grab count columns
  cnt_cols <- grep("^countData\\.", names(res_df), value = TRUE)
  counts <- as.matrix(res_df[, cnt_cols, drop=FALSE])
  colnames(counts) <- sub("^countData\\.", "", colnames(counts))  # sample names
  
  # Map samples -> gene (groupID)
  genes <- res_df$groupID
  
  # Per-gene totals per sample
  gene_totals <- rowsum(counts, genes)         # genes x samples
  # Map back to PAS order
  totals_per_row <- gene_totals[genes, colnames(counts), drop=FALSE]
  
  # Usage per PAS per sample
  usage <- counts / pmax(totals_per_row, 1)  # avoid div/0 (PAS with all-zero totals should be filtered earlier)
  
  # Condition means (using your colData)
  cond <- colData[colnames(counts), "timepoint", drop=TRUE]
  stopifnot(!any(is.na(cond)))
  day3_idx <- which(cond=="3day")
  day5_idx <- which(cond=="5day")
  
  mean_usage_3day <- if (length(day3_idx)>0) rowMeans(usage[, day3_idx, drop=FALSE]) else NA_real_
  mean_usage_5day  <- if (length(day5_idx)>0) rowMeans(usage[, day5_idx,  drop=FALSE]) else NA_real_
  mean_usage_all  <- rowMeans(usage, na.rm = TRUE)
  
  cbind(res_df,
        meanUsage_All = mean_usage_all,
        meanUsage_3day = mean_usage_3day,
        meanUsage_5day = mean_usage_5day)
}

cd_ctrl <- subset(colData, condition=="Control")
cd_trt <- subset(colData, condition=="Treatment")

res_ctrl_u <- usage_by_condition(res_ctrl, cd_ctrl)
res_trt_u <- usage_by_condition(res_trt, cd_trt)

# --- Rename L2FC column for export only (usage tables) ---
out_ctrl_u <- res_ctrl_u
out_trt_u <- res_trt_u
names(out_ctrl_u) <- sub("log2fold_5day_3day","log2fold_5day_v_3day", names(out_ctrl_u), fixed=TRUE)
names(out_trt_u) <- sub("log2fold_5day_3day","log2fold_5day_v_3day", names(out_trt_u), fixed=TRUE)

## Export the tables
write.csv(out_ctrl_u, "./APA_analysis/DEXseq/results/dexseq/3utr/with-usage/dexseq.3utr.usage.5-v-3.control.csv", row.names = FALSE, quote=FALSE)
write.csv(out_trt_u, "./APA_analysis/DEXseq/results/dexseq/3utr/with-usage/dexseq.3utr.usage.5-v-3.treatment.csv", row.names = FALSE, quote=FALSE)

library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)

plot_gene_usage <- function(res_u, gene_symbol, show_fc = TRUE, title_suffix = NULL) {
  stopifnot(all(c("gene","meanUsage_3day","meanUsage_5day") %in% names(res_u)))
  
  df <- res_u %>%
    filter(gene == gene_symbol) %>%
    # use the most reliable PAS identifier column
    mutate(PAS_raw = coalesce(.data$feature, .data$featureID)) %>%
    # keep only rows with usage info in at least one condition
    filter(!(is.na(meanUsage_3day) & is.na(meanUsage_5day))) %>%
    mutate(
      pos    = coalesce(.data$genomicData.start, NA_integer_),
      strand = as.character(coalesce(.data$genomicData.strand, NA))
    )
  
  if (nrow(df) == 0) {
    stop(sprintf("No rows for gene '%s' in the supplied table.", gene_symbol))
  }
  
  # Distal flag (vectorized; no isTRUE)
  if (!all(is.na(df$pos)) && !all(is.na(df$strand))) {
    df <- df %>%
      group_by(gene) %>%
      mutate(
        is_distal = case_when(
          strand == "+" & pos == max(pos, na.rm = TRUE) ~ TRUE,
          strand == "-" & pos == min(pos, na.rm = TRUE) ~ TRUE,
          TRUE ~ FALSE
        )
      ) %>%
      ungroup()
  } else {
    df$is_distal <- NA
  }
  
  # Ordering: by position if available; else by PAS id
  if (!all(is.na(df$pos))) {
    if (all(df$strand %in% c("+","-"), na.rm = TRUE) &&
        identical(unique(na.omit(df$strand)), "-")) {
      df <- arrange(df, desc(pos))
    } else {
      df <- arrange(df, pos)
    }
  } else {
    df <- arrange(df, PAS_raw)
  }
  
  # Label builder (vectorized) + optional L2FC
  base_lab <- ifelse(!is.na(df$is_distal) & df$is_distal,
                     paste0(df$PAS_raw, " (distal)"),
                     df$PAS_raw)
  
  lab <- base_lab
  if (show_fc && "log2fold_5day_3day" %in% names(df)) {
    lab <- sprintf("%s\nL2FC=%.2f", base_lab, df$log2fold_5day_3day)
  }
  
  # Ensure uniqueness of x labels if PAS IDs are still duplicated
  if (any(duplicated(lab))) {
    lab <- make.unique(lab, sep = "_")
  }
  
  long <- df %>%
    mutate(PAS_label = lab) %>%
    transmute(
      PAS_label,
      day3 = meanUsage_3day,
      day5 = meanUsage_5day
    ) %>%
    tidyr::pivot_longer(c(day3, day5),
                        names_to = "timepoint", values_to = "mean_usage") %>%
    mutate(
      PAS_label = factor(PAS_label, levels = unique(PAS_label)),
      mean_usage = pmax(pmin(mean_usage, 1), 0)
    )
  
  ttl <- paste0("Relative PAS usage (mean across replicates): ", gene_symbol,
                if (!is.null(title_suffix)) paste0(" — ", title_suffix) else "")
  
  ggplot(long, aes(x = PAS_label, y = mean_usage, fill = timepoint)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.75) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                                expand = expansion(mult = c(0, .05))) +
    labs(x = "PAS", y = "Mean usage (fraction of gene PAS reads)",
                  title = ttl, fill = "Timepoint") +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
                   panel.grid.minor = element_blank())
}

# ---- Plotting Genes of Interest ----
#genes <- c("FUT4", "JHY", "TGFBI", "IL20RA")
genes <- c("DVL3", "ING3", "FZD2", "BOK")

for (g in genes) {
  png(filename = paste0("./APA_analysis/DEXseq/results/plots/relative_usage/usage_bar_plot.", g, ".5-v-3.ctrl.png"), 
      width = 1000, height = 800, res=100)
  p <- plot_gene_usage(res_ctrl_u, g, show_fc = TRUE, title_suffix = "5day vs 3day, in control")
  print(p)
  dev.off()
  
  png(filename = paste0("./APA_analysis/DEXseq/results/plots/relative_usage/usage_bar_plot.", g, ".5-v-3.trt.png"), 
      width = 1000, height = 800, res=100)
  p <- plot_gene_usage(res_trt_u, g, show_fc = TRUE, title_suffix = "5day vs 3day, in treatment")
  print(p)
  dev.off()
}
