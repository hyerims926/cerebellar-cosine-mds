# obj <- readRDS("Rdata/mds_slingshot_E10.5.RDS")
# obj <- readRDS("Rdata/mds_slingshot_E11.5.RDS")
# obj <- readRDS("Rdata/mds_slingshot_E12.5.RDS")
obj <- readRDS("Rdata/mds_slingshot_E13.5.RDS")
# obj <- readRDS("Rdata/mds_slingshot_E14.5.RDS")
# obj <- readRDS("Rdata/mds_slingshot_E15.5.RDS")
obj <- readRDS("Rdata/mds_slingshot_E17.5.RDS")
# obj <- readRDS("Rdata/mds_slingshot_P0.RDS")
# obj <- readRDS("Rdata/mds_slingshot_P4.RDS")
# obj <- readRDS("Rdata/mds_slingshot_P7.RDS")
# obj <- readRDS("Rdata/mds_slingshot_P14.RDS")
# obj <- readRDS("Rdata/mds_slingshot_adult (9 weeks).RDS")

markers_A <- c(
  "Ptf1a", "Lhx1", "Lhx5", "Corl2", "Skor2", "Foxp2",
  "Rora", "Calb1", "Pcp2", "Grid2", "Itpr1", "Car8"
)  # Purkinje lineage markers

markers_A <- c(
  "Atoh1", "Sox4", "Sox11", "Dcx", "Pou3f2", "Neurod1",
  "Neurod2", "Tubb3", "Slc17a7", "Rbfox3", "Snap25", "Camk2a"
)  # granule / general neuronal differentiation markers

markers_A <- AnnotationDbi::select(
  org.Mm.eg.db,
  keys = markers_A,
  keytype = "SYMBOL",
  columns = c("ENSEMBL", "SYMBOL"))

markers_A <- markers_A[["ENSEMBL"]]
markers_A <- intersect(markers_A, rownames(obj))

mat <- LayerData(obj, assay = "RNA", layer = "data")
expr_A <- mat[markers_A, , drop = FALSE]


########purk #########

df_A <- data.frame(
  pt_mds_purk = obj$pt_mds_purk,
  pt_umaph_purk = obj$pt_umaph_purk,
  t(expr_A))

df_A$bin_mds <- cut(
  df_A$pt_mds_purk,
  breaks = unique(quantile(df_A$pt_mds_purk, probs = seq(0, 1, length.out = 100), na.rm = TRUE)),
  include.lowest = TRUE)

df_bin_mds <- df_A %>%
  group_by(bin_mds) %>%
  summarise(across(all_of(markers_A), mean))

df_A$bin_mds <- cut(
  df_A$pt_umaph_purk,
  breaks = unique(quantile(df_A$pt_umaph_purk, probs = seq(0, 1, length.out = 100), na.rm = TRUE)),
  include.lowest = TRUE)

df_bin_umaph <- df_A %>%
  group_by(bin_mds) %>%
  summarise(across(all_of(markers_A), mean))

#####################

######gra ##########

df_A <- data.frame(
  pt_mds_gra = obj$pt_mds_gra,
  pt_umaph_gra = obj$pt_umaph_gra,
  t(expr_A))

df_A$bin_mds <- cut(
  df_A$pt_mds_gra,
  breaks = unique(quantile(df_A$pt_mds_gra, probs = seq(0, 1, length.out = 100), na.rm = TRUE)),
  include.lowest = TRUE)

df_bin_mds <- df_A %>%
  group_by(bin_mds) %>%
  summarise(across(all_of(markers_A), mean))

df_A$bin_umaph <- cut(
  df_A$pt_umaph_gra,
  breaks = unique(quantile(df_A$pt_umaph_gra, probs = seq(0, 1, length.out = 100), na.rm = TRUE)),
  include.lowest = TRUE)

df_bin_umaph <- df_A %>%
  group_by(bin_umaph) %>%
  summarise(across(all_of(markers_A), mean))


######################

mat_mds <- t(as.matrix(df_bin_mds[, markers_A]))
mat_umaph <- t(as.matrix(df_bin_umaph[, markers_A]))



# matrix Pearson correlation
cor(as.vector(mat_mds),
    as.vector(mat_umaph),
    method = "pearson")

# [1] 0.9899481
# [1] 0.9812594
# [1] 0.9717628
# [1] 0.9100348



# Spearman
cor(as.vector(mat_mds),
    as.vector(mat_umaph),
    method = "spearman")

# [1] 0.9565038
# [1] 0.9600262
# [1] 0.8976137
# [1] 0.901805



p_mds <- pheatmap(
  mat_mds,
  cluster_cols = FALSE,
  cluster_rows = FALSE)

p_umaph <- pheatmap(
  mat_umaph,
  cluster_cols = FALSE,
  cluster_rows = FALSE)

