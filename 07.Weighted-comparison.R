seed <- 1234

obj <- readRDS("Rdata/obj_all.RDS")

obj <- NormalizeData(obj) 
obj <- FindVariableFeatures(obj, nfeatures = 3000)
obj <- ScaleData(obj, features = VariableFeatures(obj))

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

min_cells <-3

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

stage_a <- "E11.5"; stage_b <- "E12.5"
stage_a <- "E12.5"; stage_b <- "E13.5"
stage_a <- "E13.5"; stage_b <- "E14.5"
stage_a <- "E14.5"; stage_b <- "E15.5"
stage_a <- "E15.5"; stage_b <- "E17.5"
stage_a <- "E17.5"; stage_b <- "P0"
stage_a <- "P0"; stage_b <- "P4"
stage_a <- "P4"; stage_b <- "P7"


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

pair_idx <- which(
  upper_idx,
  arr.ind = TRUE
)

anno_i <- rownames(dist_a)[pair_idx[, 1]]
anno_j <- colnames(dist_a)[pair_idx[, 2]]


group_count_all <- meta_use %>%
  dplyr::count(
    stage,
    annotation,
    name = "n_cells"
  )

prop_a <- group_count_all %>%
  filter(stage == stage_a) %>%
  mutate(
    prop = n_cells / sum(n_cells)
  ) %>%
  filter(
    annotation %in% common_after_norm
  ) %>%
  select(
    annotation,
    prop
  )


prop_b <- group_count_all %>%
  filter(stage == stage_b) %>%
  mutate(
    prop = n_cells / sum(n_cells)
  ) %>%
  filter(
    annotation %in% common_after_norm
  ) %>%
  select(
    annotation,
    prop
  )


prop_mean <- full_join(
  prop_a,
  prop_b,
  by = "annotation",
  suffix = c("_a", "_b")
) %>%
  mutate(
    prop_a = replace_na(prop_a, 0),
    prop_b = replace_na(prop_b, 0),
    mean_prop = (prop_a + prop_b) / 2
  )


prop_lookup <- setNames(
  prop_mean$mean_prop,
  prop_mean$annotation
)

pair_weight <- (
  prop_lookup[anno_i] *
    prop_lookup[anno_j]
)


df_pair <- data.frame(
  anno_i = anno_i,
  anno_j = anno_j,
  dissim_a = dist_a[upper_idx],
  dissim_b = dist_b[upper_idx],
  pair_weight = pair_weight
) %>%
  filter(
    is.finite(dissim_a),
    is.finite(dissim_b),
    is.finite(pair_weight),
    pair_weight > 0
  )


# linear model
fit <- lm(
  dissim_b ~ dissim_a,
  data = df_pair,
  weights = df_pair$pair_weight
)

summary(fit)

coef(fit)

# 핵심 값
slope <- coef(fit)[2]
intercept <- coef(fit)[1]
r2 <- summary(fit)$r.squared

slope
intercept
r2


ggplot(
  df_pair,
  aes(
    x = dissim_a,
    y = dissim_b,
    size = pair_weight
  )
) +
  geom_point(
    alpha = 0.7
  ) +
  geom_abline(
    intercept = intercept,
    slope = slope,
    color = "blue",
    linewidth = 0.8
  ) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed"
  ) +
  scale_size_continuous(
    range = c(1, 8),
    name = "Pair weight"
  ) +
  coord_cartesian(
    xlim = c(0, 0.6),
    ylim = c(0, 0.6)
  ) +
  labs(
    x = paste0(
      stage_a,
      " inter-annotation cosine dissimilarity"
    ),
    y = paste0(
      stage_b,
      " inter-annotation cosine dissimilarity"
    ),
    title = paste0(
      "weighted slope = ", round(slope, 3),
      ", weighted R² = ", round(r2, 3))) +
  theme_classic() + NoLegend()


ggsave(paste0("figures/",stage_a, stage_b, "correlation.png"), height = 4, width = 5)

