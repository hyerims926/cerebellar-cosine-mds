seed <- 1234
set.seed(seed)

# obj <- readRDS("Rdata/mds_slingshot_E10.5.RDS")
# obj <- readRDS("Rdata/mds_slingshot_E11.5.RDS")
obj <- readRDS("Rdata/mds_slingshot_E12.5.RDS")
# obj <- readRDS("Rdata/mds_slingshot_E13.5.RDS")
# obj <- readRDS("Rdata/mds_slingshot_E14.5.RDS")
# obj <- readRDS("Rdata/mds_slingshot_E15.5.RDS")
# obj <- readRDS("Rdata/mds_slingshot_E17.5.RDS")
obj <- readRDS("Rdata/mds_slingshot_P0.RDS")
# obj <- readRDS("Rdata/mds_slingshot_P4.RDS")
# obj <- readRDS("Rdata/mds_slingshot_P7.RDS")
# obj <- readRDS("Rdata/mds_slingshot_P14.RDS")
# obj <- readRDS("Rdata/mds_slingshot_adult (9 weeks).RDS")

counts <- GetAssayData(obj, assay = "RNA", layer = "counts")


# pt_mds_purk / pt_umaph_purk / pt_mds_gra / pt_umaph_gra
cells_use <- colnames(obj)[!is.na(obj$pt_mds_gra)]
pt <- obj$pt_mds_gra[cells_use]

cells_use <- colnames(obj)[!is.na(obj$pt_umaph_gra)]
pt <- obj$pt_umaph_gra[cells_use]

counts_use <- counts[, cells_use]

pseudotime <- matrix(pt, ncol = 1)
rownames(pseudotime) <- cells_use

#"gra" or"purk"
colnames(pseudotime) <- "gra"

cellWeights <- matrix(1, nrow = length(cells_use), ncol = 1)
rownames(cellWeights) <- cells_use

#"gra" or"purk"
colnames(cellWeights) <- "gra"

# --- marker genes ---
# gra
genes_plot_symbol <- c( "Top2a","Atoh1", "Barhl1", "Neurod1","Dcx","Cntn2")

# purk
# genes_plot_symbol <- c("Top2a", "Pax6", "Foxp2", "Esrrg", "Dpp10", "Fgf14")

gene_map_plot <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys = genes_plot_symbol,
  keytype = "SYMBOL",
  columns = c("SYMBOL", "ENSEMBL")
)

genes_plot_ens <- gene_map_plot$ENSEMBL[
  match(genes_plot_symbol, gene_map_plot$SYMBOL)
]
genes_plot_ens

fit <- fitGAM(
  counts = counts_use,
  genes = genes_plot_ens,
  pseudotime = pseudotime,
  cellWeights = cellWeights,
  nknots = 6,
  verbose = TRUE,
  BPPARAM = BiocParallel::MulticoreParam(workers = 4))


genes_plot_ens %in% rownames(fit)

p1 <- plotSmoothers(fit, counts = counts_use, gene = genes_plot_ens[1]) + NoLegend() + ggtitle(genes_plot_symbol[1]) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) + ylab("log1p count") + scale_color_manual(values = "grey40")

p2 <- plotSmoothers(fit, counts = counts_use, gene = genes_plot_ens[2]) + NoLegend() + ggtitle(genes_plot_symbol[2]) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) + ylab("log1p count") + scale_color_manual(values = "grey25")

p3 <- plotSmoothers(fit, counts = counts_use, gene = genes_plot_ens[3]) + NoLegend() + ggtitle(genes_plot_symbol[3]) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) + ylab("log1p count") + scale_color_manual(values = "grey25")

p4 <- plotSmoothers(fit, counts = counts_use, gene = genes_plot_ens[4]) + NoLegend() + ggtitle(genes_plot_symbol[4]) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) + ylab("log1p count") + scale_color_manual(values = "grey25")

p5 <- plotSmoothers(fit, counts = counts_use, gene = genes_plot_ens[5]) + NoLegend() + ggtitle(genes_plot_symbol[5]) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) + ylab("log1p count") + scale_color_manual(values = "grey25")

p6 <- plotSmoothers(fit, counts = counts_use, gene = genes_plot_ens[6]) + NoLegend() + ggtitle(genes_plot_symbol[6]) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold")) + ylab("log1p count") + scale_color_manual(values = "grey25")

(p1 | p3 | p5)
ggsave("figures/tradeseq_umaph_gra_P0.png", width = 6, height = 2.25)


(p1 | p2 | p3) / (p4 | p5 | p6)
ggsave("figures/tradeseq_mds_gra_P0.png", width = 8, height = 6)
