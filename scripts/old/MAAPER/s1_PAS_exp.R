i## 2025/5/27
## objective: a. annotate gene; b. calculate PAS RPM, PAS usage(%), 3'UTR RPM, 
## inputfiles: PAS count, PAS reference (such as polyAdb files)
library(dplyr)
library(openxlsx)
library(janitor)
#=======================================================
# set up 
#=======================================================

mainpath <- "/Volumes/biocore01/scratch/cnobrega/3_REAP/APA_analysis/MAAPER";setwd(mainpath)
outdir <- "results/"
if (!dir.exists(outdir)) dir.create(outdir)

# sample base names
sample_names <- read.table("sample_name.txt")[,1] 

# PAS count file
exp_file <- "cluster.all.reads.csv" 

# PAS reference (polyAdb files) 
pas_ref.file <- "human.PAS.hg38.txt"

#=======================================================
# step1 data preprocessing 
#=======================================================
message("step1 data preprocessing")

# PAS count table
exp <- read.csv(exp_file) # first col is hit_PAS_ID, other cols are sample name as vector sample_names.

# define count columns with suffix "_count"
colnames(exp)[match(sample_names, colnames(exp))] <- paste0(sample_names, "_count")

count_cols <- paste0(sample_names, "_count")

# remove all 0 rows for subsequent quantile and only use genetic reads
# remove all 0 rows for subsequent quantile and only use genetic reads
pas <- exp[rowSums(exp[, count_cols, drop = FALSE]) > 0, ]
pas <- pas[!grepl("^chrM", pas$hit_PAS_ID), ]

# PAS reference (polyAdbv4 files). Make some selection if u need.
pas_ref0 <- read.csv(pas_ref.file, sep="\t")
pas_ref0 <- clean_names(pas_ref0, case="none")

# ========== HANDLING ANNOTATION DUPLICATES ==========
# create one single new column: gene_id == any kind of gene identifier in the annotation table 
# to determine gene_id: try first Gene Symbol, if NA try Ensembl ID, if NA try RefSeq Gene ID, if NA try FAMTOM ID. If all empty, set NA, then drop rows with NA gene_id

blank_to_na <- function(x) {
  x <- str_squish(x)
  x[tolower(x) %in% c("", "na")] <- NA_character_
  x
}
pas_ref0 <- pas_ref0 %>%
  mutate(
    across(c(`Gene_Symbol`, `Ensemble_ID`, `RefSeq_Gene_ID`, `FAMTOM_ID`), blank_to_na),
    gene_id = coalesce(`Gene_Symbol`, `Ensemble_ID`, `RefSeq_Gene_ID`, `FAMTOM_ID`)
  )

# Find PAS_ID values that are in pas (counts table) but not in pas_ref0 if you were to remove all na entries
pas_ref0_dropna <- pas_ref0 %>%
  filter(!is.na(gene_id))
mismatched_pas_ids <- setdiff(pas$hit_PAS_ID, pas_ref0_dropna$PAS_ID)
# Assign gene_id = PAS_ID for those NA rows that you need to keep
pas_ref0$gene_id[pas_ref0$PAS_ID %in% mismatched_pas_ids & is.na(pas_ref0$gene_id)] <- 
  pas_ref0$PAS_ID[pas_ref0$PAS_ID %in% mismatched_pas_ids & is.na(pas_ref0$gene_id)]

# drop NA group_IDs (already handled the needed NA values above, the rest are not needed)
pas_ref0 <- pas_ref0 %>%
  filter(!is.na(gene_id))

# now check if there are any illegal duplicates (same gene name and same PASS_ID)
illegal_dup_check <- pas_ref0 %>%
  group_by(PAS_ID, gene_id) %>%
  tally() %>%
  filter(n > 1)
# If any exist, keep only the first occurrence of each PAS_ID + gene_id pair
if (nrow(illegal_dup_check) > 0) {
  pas_ref0 <- pas_ref0[!duplicated(pas_ref0[c("PAS_ID", "gene_id")]), ]
}

# add additional metadata columns for use in downstream scripts: 
# Chromosome -> rename to "chromosome" 
# keep Position
# add "start" & "end"-> value should be == Position
# keep PAS_ID
# add "Peak_ID" -> value should be == PAS_ID
# keep gene_id
# Strand -> rename to "strand"
# Intron_exon_location -> rename to "region"
needed_cols <- c("Chromosome", "Position", "PAS_ID", "Strand", "Intron_exon_location", "gene_id")
missing_cols <- setdiff(needed_cols, names(pas_ref0))
if (length(missing_cols) > 0) {
  stop(sprintf("Missing columns in pas_ref0: %s", paste(missing_cols, collapse = ", ")))
}

pas_ref <- pas_ref0 %>%
  transmute(
    chromosome = Chromosome,
    Position   = as.integer(Position),
    start      = as.integer(Position),
    end        = as.integer(Position),
    PAS_ID     = PAS_ID,
    Peak_ID    = PAS_ID,
    gene_id    = gene_id,
    strand     = Strand,
    region     = Intron_exon_location, 
    PAS_type = PAS_type
  )

head(sort(table(pas_ref$gene_id), decreasing = TRUE), 10) # check if there's invaild gene name such as NA,"", ".", unknown,etc.
pas_id="PAS_ID"
#=======================================================
# step2 Calculate PAS RPM, using all genetic reads 
#=======================================================
print("step2 Calculate PAS RPM")

for(sample_name in sample_names){
  #sample_name = sample_names[1]
  all_counts <- pas[[paste0(sample_name, "_count")]]
  ## between 0.5 and 0.95 quantile
  range = quantile(all_counts, c(0.05, 0.95))
  trimmed_counts = all_counts[all_counts >= range[1] & all_counts <= range[2]]
  pas[[paste0(sample_name, "_rpm")]] = 10^6*all_counts/sum(trimmed_counts)
}

print(colnames(pas))
#pas0=pas
# pas=pas0
#=======================================================
# step3 annotate PAS and Calculate read number/RPM per gene
#=======================================================
pas <- merge(pas_ref,pas,by.y="hit_PAS_ID",by.x = pas_id,all.y=T,sort=F)
print(nrow(pas))
print(head(pas,2))

#------------------------
# Total read number per gene, without including reads from UA pA sites
# Calculate total reads number per gene, without counting pAs in UA 
counts_per_gene = aggregate(pas[, count_cols, drop=FALSE], 
                            list(gene_id = pas$gene_id), sum)
print(paste0("gene number: ", nrow(counts_per_gene)))
names(counts_per_gene)[-1] = paste0(sample_names, "_geneCount")

pas = merge(pas, counts_per_gene,by="gene_id", all.x = T, sort = F)

#rm(counts_per_gene)
print(head(pas,2))
print(colnames(pas))

#------------------------
# Calculate genewise RPM
rpm_per_gene = aggregate(pas[,paste0(sample_names, "_rpm")], 
                         list(gene_id = pas$gene_id), sum)
names(rpm_per_gene)[-1] = paste0(sample_names, "_geneRPM")

pas = merge(pas, rpm_per_gene, by="gene_id", all.x = T, sort = F)

print(head(pas,2))
print(colnames(pas))
#------------------------
# pA site usage (fraction of PAS reads from each pA site among all pA sites of the same gene)
# Calculate pA usage
for(sample_name in sample_names){
  pas[[paste0(sample_name, "_usage")]] <- pas[[paste0(sample_name, "_count")]]*100/pas[[paste0(sample_name, "_geneCount")]]
}
# Replace NaN with 0 across all "_usage" columns
usage_cols <- grep("_usage$", names(pas))
pas[usage_cols] <- lapply(pas[usage_cols], function(x) ifelse(is.nan(x), 0, x))

print(head(pas,2))
print(colnames(pas))

#------------------------
# Total 3'UTR reads and 3'UTR RPMs per gene
#### Calculate total 3'UTR reads per gene
counts_in_3utr_per_gene = aggregate(pas[grepl("^3'UTR", pas$PAS_type), count_cols, drop = FALSE],
                                    list(gene_id = pas$gene_id[grepl("^3'UTR", pas$PAS_type)]),
                                    sum)
names(counts_in_3utr_per_gene)[-1] = paste0(names(counts_in_3utr_per_gene)[-1], "_3UTRcount")
pas = merge(pas, counts_in_3utr_per_gene,by="gene_id", all.x=T, sort = F)
#rm(counts_in_3utr_per_gene)
print(head(pas,2))
print(colnames(pas))

#------------------------
# Calculate total 3'UTR RPMs per gene
rpm_in_3utr_per_gene = aggregate(pas[grepl("^3'UTR", pas$PAS_type), paste0(sample_names, "_rpm")], 
                                 list(gene_id = pas$gene_id[grepl("^3'UTR", pas$PAS_type)]),
                                 sum)

names(rpm_in_3utr_per_gene)[-1] = paste0(names(rpm_in_3utr_per_gene)[-1], "_3UTRrpm")

pas = merge(pas, rpm_in_3utr_per_gene,by="gene_id", all.x=T, sort = F)
pas1=pas

print(head(pas,2))
print(colnames(pas))

#------------------------
saveRDS(pas, paste0(outdir,"pas_exp.rds"))
openxlsx::write.xlsx(pas,paste0(outdir,"pas_exp.xlsx"))
write.csv(pas, paste0(outdir, "pas_exp.csv"), quote = FALSE, row.names = FALSE)
