seed <- 1234

obj <- readRDS("Rdata/obj_all.RDS")

obj <- NormalizeData(obj) 
obj <- FindVariableFeatures(obj, nfeatures = 3000)

## stage

stage_levels <- c(
  "E11.5", "E12.5", "E13.5", "E14.5",
  "E15.5", "E17.5", "P0", "P4", "P7"
)

stage_pairs <- data.frame(
  stage_a = head(stage_levels, -1),
  stage_b = tail(stage_levels, -1)
) %>%
  mutate(
    comparison = paste0(stage_a, "-", stage_b)
  )

DefaultAssay(obj) <- "RNA"

expr_mat <- GetAssayData(
  obj,
  assay = "RNA",
  layer = "data"
)


genes_use <- intersect(VariableFeatures(obj), rownames(expr_mat))

expr_mat <- expr_mat[
  genes_use,
  ,
  drop = FALSE]

meta_use <- obj@meta.data %>%
  mutate(
    cell = rownames(.),
    stage = as.character(stage),
    annotation = as.character(anno)
  ) %>%
  filter(
    stage %in% stage_levels,
    !is.na(annotation)
  )

min_cells <- 3

group_count <- meta_use %>%
  dplyr::count(
    stage,
    annotation,
    name = "n_cells"
  ) %>%
  filter(n_cells >= min_cells)

meta_valid <- meta_use %>%
  inner_join(
    group_count,
    by = c("stage", "annotation")
  )


group_key <- interaction(
  meta_valid$annotation,
  meta_valid$stage,
  sep = "___",
  drop = TRUE
)

group_cells <- split(
  meta_valid$cell,
  group_key
)

pseudobulk_mat <- sapply(
  group_cells,
  function(cells) {
    Matrix::rowMeans(
      expr_mat[, cells, drop = FALSE]
    )
  }
)

stage_a <- "P0"
stage_b <- "P4"

suffix_a <- paste0("___", stage_a, "$")
suffix_b <- paste0("___", stage_b, "$")

cols_a <- grep(
  suffix_a,
  colnames(pseudobulk_mat),
  value = TRUE
)

cols_b <- grep(
  suffix_b,
  colnames(pseudobulk_mat),
  value = TRUE
)

anno_a <- sub(suffix_a, "", cols_a)
anno_b <- sub(suffix_b, "", cols_b)

common_anno <- intersect(anno_a, anno_b)

cols_a_use <- paste0(common_anno, "___", stage_a)
cols_b_use <- paste0(common_anno, "___", stage_b)

mat_a <- pseudobulk_mat[, cols_a_use, drop = FALSE]
mat_b <- pseudobulk_mat[, cols_b_use, drop = FALSE]

colnames(mat_a) <- common_anno
colnames(mat_b) <- common_anno


cosine_distance_matrix <- function(mat) {
  
  mat <- as.matrix(mat)
  storage.mode(mat) <- "double"
  
  norms <- sqrt(colSums(mat^2))
  valid <- is.finite(norms) & norms > 0
  
  mat <- mat[, valid, drop = FALSE]
  norms <- norms[valid]
  
  mat_unit <- sweep(
    mat,
    MARGIN = 2,
    STATS = norms,
    FUN = "/"
  )
  
  sim <- crossprod(mat_unit)
  
  sim[sim > 1] <- 1
  sim[sim < -1] <- -1
  
  dist_mat <- 1 - sim
  diag(dist_mat) <- 0
  
  dist_mat
}


dist_a <- cosine_distance_matrix(mat_a)
dist_b <- cosine_distance_matrix(mat_b)

common_after_norm <- intersect(
  colnames(dist_a),
  colnames(dist_b)
)

dist_a <- dist_a[
  common_after_norm,
  common_after_norm,
  drop = FALSE
]

dist_b <- dist_b[
  common_after_norm,
  common_after_norm,
  drop = FALSE
]

upper_idx <- upper.tri(dist_a)

df_pair <- data.frame(
  dissim_P0 = dist_a[upper_idx],
  dissim_P4 = dist_b[upper_idx]
)

df_pair <- df_pair %>%
  filter(
    is.finite(dissim_P0),
    is.finite(dissim_P4)
  )


# linear model
fit <- lm(dissim_P4 ~ dissim_P0, data = df_pair)

summary(fit)

coef(fit)

# 핵심 값
slope <- coef(fit)[2]
intercept <- coef(fit)[1]
r2 <- summary(fit)$r.squared

slope
intercept
r2

cor(
  df_pair$dissim_P0,
  df_pair$dissim_P4,
  method = "pearson"
)

cor(
  df_pair$dissim_P0,
  df_pair$dissim_P4,
  method = "spearman"
)


# plot
ggplot(
  df_pair,
  aes(
    x = dissim_P0,
    y = dissim_P4
  )
) +
  geom_point(
    alpha = 0.7,
    size = 2
  ) +
  geom_smooth(
    method = "lm",
    se = FALSE
  ) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed"
  ) +
  labs(
    x = "P0 inter-annotation cosine dissimilarity",
    y = "P4 inter-annotation cosine dissimilarity",
    title = paste0(
      "slope = ", round(slope, 3),
      ", intercept = ", round(intercept, 3),
      ", R² = ", round(r2, 3)
    )
  ) +
  theme_classic()

ggsave("figures/P0P4-correlation.png", height =4, width =5)
