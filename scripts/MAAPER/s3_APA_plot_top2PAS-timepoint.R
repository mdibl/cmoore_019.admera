## 2025/2/28
## objective: plot 
## inputfiles: APA results

#### SCRIPT ----
## Definitions:
# REDu measures the relative expression levels between the top two most abundant APA isoforms in the 3'-most exon only.
# REDi measures the relative expression levels between the top expressed isoform in the 3'-most exon and the top expressed isoform in an intron.

#### Loading libraries: ----
library(reshape2)
library(tidyverse)
library(ggrepel)
library(openxlsx)

#pseudoCount = 1
pseudoCount = 0

mainpath = paste0("/Volumes/biocore01/scratch/cnobrega/3_REAP/APA_analysis/MAAPER/results/5d-v-3d");setwd(mainpath)
inpath1 = paste0(mainpath,"/DE-UPA_UPAS-avUsage0/")
inpath2 = paste0(mainpath,"/DE-IPA_UPAS-avUsage10_IPAS-avUsage5/")
outpath1 = paste0(inpath1,"pseudoC",pseudoCount,"/");dir.create(outpath1,r=T)
outpath2 = paste0(inpath2,"pseudoC",pseudoCount,"/");dir.create(outpath2,r=T)


#### APA plot: ----
colors = c("#E41A1C", "grey", "#1F78B4")


#========================================
# 3'UTR-APA
#========================================
# Loading files: ----
REDu_aaa <- readRDS(paste0(inpath1,"APA_3UTR-pseudoC0.rds"))

REDu_aaa$regu_select <- REDu_aaa$regu.adj;filename = "ReguAdj" # optimal, dependent on differential APA-gene number
#REDu_aaa$regu_select <- REDu_aaa$regu;filename = "Regu"

REDu_aaa$regu_select <- factor(REDu_aaa$regu_select, levels = c("Lengthened","NO","Shortened"))
# 首先按分组计算每个分面的统计信息
facet_labels <- REDu_aaa %>%
  group_by(Comparison) %>%
  summarise(
    Lengthened = sum(regu_select == "Lengthened"),
    Shortened = sum(regu_select == "Shortened"),
    NO = sum(regu_select == "NO"),
    ratio = ifelse(Lengthened > Shortened, 
                   round(Lengthened / Shortened, 1), 
                   -round(Shortened / Lengthened, 1)),
    label = paste(
      "red=", Lengthened,
      #", gray=", NO,
      ", blue=", Shortened,
      ", ratio=", ratio, 
      sep = ""
    )
  )

# 
facet_labels$label <- paste0(facet_labels$Comparison,"\n",facet_labels$label)
facet_label_map <- setNames(facet_labels$label, facet_labels$Comparison)
write.xlsx(as.data.frame(facet_labels),file.path(outpath1, paste0("APA_3UTR-pseudoC",pseudoCount,"_",filename,"_summary.xlsx")))
# 
comparisons <- c("Treat5d.vs.Treat3d","Ctrl5d.vs.Ctrl3d")
REDu_aaa$Comparison <- factor(REDu_aaa$Comparison,levels=comparisons)

lim <- max(abs(c(REDu_aaa$RE_pPAS, REDu_aaa$RE_dPAS)), na.rm = TRUE)
n_facets <- length(comparisons)

ggplot(data = REDu_aaa, 
       aes(x = as.numeric(RE_pPAS), 
           y = as.numeric(RE_dPAS), 
           color = regu_select, fill = regu_select, size = regu_select)) +
  geom_point(alpha = 0.3) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  scale_size_manual(values = c(1.5, 1, 1.5)) +
  facet_wrap(Comparison ~ ., nrow = 1,
             labeller = labeller(Comparison = facet_label_map)) + # 修改分面标题
  coord_fixed(ratio = 1) +
  scale_x_continuous(limits = c(-lim, lim)) +                 # same limits on both axes
  scale_y_continuous(limits = c(-lim, lim)) +
  theme_classic() +
  theme(
    axis.text = element_text(colour = "black"),
    axis.title = element_text(colour = "black"),
    legend.position = "none",
    plot.title = element_text(size = 8),
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(size = 1),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, colour = "black", size = 1)
  ) +
  labs(
    x = "Log2Ratio of pPAS isoforms", 
    y = "Log2Ratio of dPAS isoforms", 
    color = ""
  ) 
ggsave(file.path(outpath1, paste0("APA_3UTR-pseudoC", pseudoCount, "_", filename, "_v1.pdf")), width = 4 * n_facets, height = 4)


#========================================
# IPA
#========================================
REDi_aaa <- readRDS(paste0(inpath2,"APA_IPA-pseudoC0.rds"))

REDi_aaa$regu_select <- REDi_aaa$regu.adj;filename = "ReguAdj"
#REDi_aaa$regu_select <- REDi_aaa$regu;filename = "Regu"
REDi_aaa$regu_select <- factor(REDi_aaa$regu_select, levels = c("IPAS_increased","NO","IPAS_decreased"))

# 首先按分组计算每个分面的统计信息
REDi_aaa$Comparison <- factor(REDi_aaa$Comparison,levels=comparisons)
facet_labels <- REDi_aaa %>%
  group_by(Comparison) %>%
  summarise(
    IPAS_decreased = sum(regu_select == "IPAS_decreased"),
    IPAS_increased = sum(regu_select == "IPAS_increased"),
    NO = sum(regu_select == "NO"),
    ratio = ifelse(IPAS_decreased > IPAS_increased, 
                   -round(IPAS_decreased / IPAS_increased, 1), 
                   round(IPAS_increased / IPAS_decreased, 1)),
    label = paste(
      "red=", IPAS_increased,
      #", gray=", NO,
      ", blue=", IPAS_decreased,
      ", ratio=", ratio, 
      sep = ""
    )
  )

# 
facet_labels$label <- paste0(facet_labels$Comparison,"\n",facet_labels$label)
facet_label_map <- setNames(facet_labels$label, facet_labels$Comparison)
write.xlsx(as.data.frame(facet_labels),file.path(outpath2, paste0("APA_IPA-pseudoC",pseudoCount,"_",filename,"_summary.xlsx")))
# 
#REDi_aaa$Comparison <- factor(REDi_aaa$Comparison,levels=comparisons)
lim <- max(abs(c(REDi_aaa$RE_IPA, REDi_aaa$RE_TPA)), na.rm = TRUE)
n_facets <- length(comparisons)

ggplot(data = REDi_aaa, 
       aes(x = as.numeric(RE_IPA), 
           y = as.numeric(RE_TPA), 
           color = regu_select, fill = regu_select, size = regu_select)) +
  geom_point(alpha = 0.3) +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  scale_size_manual(values = c(1.5, 1, 1.5)) +
  facet_wrap(Comparison ~ ., nrow = 1,
             labeller = labeller(Comparison = facet_label_map)) + # 修改分面标题
  coord_fixed(ratio = 1) +
  scale_x_continuous(limits = c(-lim, lim)) +                 # same limits on both axes
  scale_y_continuous(limits = c(-lim, lim)) +
  theme_classic() +
  theme(
    axis.text = element_text(colour = "black"),
    axis.title = element_text(colour = "black"),
    legend.position = "none",
    plot.title = element_text(size = 8),
    axis.line = element_line(colour = "black"),
    axis.ticks = element_line(size = 1),
    panel.grid.major = element_blank(), 
    panel.grid.minor = element_blank(),
    panel.border = element_rect(fill = NA, colour = "black", size = 1)
  ) +
  #coord_fixed(ratio = 1) + 
  #scale_x_continuous(limits = c(-3, 3), breaks = c(-4, -2, 0, 2, 4)) +
  #scale_y_continuous(limits = c(-3, 3), breaks = c(-4, -2, 0, 2, 4)) +
  labs(x = "Log2Ratio of IPA isoforms", y = "Log2Ratio of TPA isoforms", color = "")
ggsave(file.path(outpath2, paste0("APA_IPA-pseudoC", pseudoCount, "_", filename, "_v1.pdf")), width = 4 * n_facets, height = 4)








