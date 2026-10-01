
obj <- readRDS("Rdata/obj_all.RDS")

obj <- NormalizeData(obj) 
obj <- FindVariableFeatures(obj, nfeatures = 3000)
obj <- ScaleData(obj, features = VariableFeatures(obj))



compare_stage_relationships <- function(
    obj,
    stage_pairs,
    stage_col = "stage",
    anno_col = "anno",
    assay = "RNA",
    layer = "data",
    min_cells = 0,
    weighted = TRUE
) {
  
  DefaultAssay(obj) <- assay
  
  expr_mat <- GetAssayData(
    obj,
    assay = assay,
    layer = layer)
  
  genes_use <- intersect(VariableFeatures(obj),rownames(expr_mat))
  
  expr_mat <- expr_mat[
    genes_use,
    ,
    drop = FALSE]
  
  meta <- obj@meta.data %>%
    mutate(
      cell = rownames(.),
      stage = as.character(.data[[stage_col]]),
      annotation = as.character(.data[[anno_col]])
    ) %>%
    filter(
      !is.na(annotation))
  
  group_count_all <- meta %>%
    dplyr::count(
      stage,
      annotation,
      name = "n_cells"
    )
  
  group_count <- group_count_all %>%
    filter(n_cells >= min_cells)
  
  meta_valid <- meta %>%
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
        expr_mat[, cells, drop = FALSE])})

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
    
    dist_mat}
  
  compare_two_stages <- function(stage_a, stage_b) {
    
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
    
    anno_a <- sub(
      paste0("___", stage_a, "$"),
      "",
      cols_a
    )
    
    anno_b <- sub(
      paste0("___", stage_b, "$"),
      "",
      cols_b
    )
    
    common_anno <- intersect(
      anno_a,
      anno_b
    )
  
    
    cols_a_use <- paste0(
      common_anno,
      "___",
      stage_a
    )
    
    cols_b_use <- paste0(
      common_anno,
      "___",
      stage_b
    )
    
    mat_a <- pseudobulk_mat[
      ,
      cols_a_use,
      drop = FALSE
    ]
    
    mat_b <- pseudobulk_mat[
      ,
      cols_b_use,
      drop = FALSE
    ]
    
    colnames(mat_a) <- common_anno
    colnames(mat_b) <- common_anno
    
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
    
    vec_a <- dist_a[upper_idx]
    vec_b <- dist_b[upper_idx]
    
    if (weighted == TRUE) {
      
      weighted_cor <- function(x, y, w) {
        
        valid <- is.finite(x) &
          is.finite(y) &
          is.finite(w) &
          w > 0
        
        x <- x[valid]
        y <- y[valid]
        w <- w[valid]
        
        w <- w / sum(w)
        
        mx <- sum(w * x)
        my <- sum(w * y)
        
        cov_xy <- sum(w * (x - mx) * (y - my))
        var_x <- sum(w * (x - mx)^2)
        var_y <- sum(w * (y - my)^2)
        
        cov_xy / sqrt(var_x * var_y)
      }
      
      
      weighted_spearman <- function(x, y, w) {
        
        weighted_cor(
          rank(x, ties.method = "average"),
          rank(y, ties.method = "average"),
          w
        )
      }
      
      
      weighted_median <- function(x, w) {
        
        valid <- is.finite(x) &
          is.finite(w) &
          w > 0
        
        x <- x[valid]
        w <- w[valid]
        
        ord <- order(x)
        
        x <- x[ord]
        w <- w[ord]
        
        w <- w / sum(w)
        
        x[which(cumsum(w) >= 0.5)[1]]
      }

      pair_idx <- which(upper_idx, arr.ind = TRUE)
      
      anno_i <- rownames(dist_a)[pair_idx[, 1]]
      anno_j <- colnames(dist_a)[pair_idx[, 2]]
      
      prop_a <- group_count_all %>%
        filter(stage == stage_a) %>%
        mutate(prop = n_cells / sum(n_cells)) %>%
        filter(annotation %in% common_after_norm) %>%
        select(annotation, prop)
      
      prop_b <- group_count_all %>%
        filter(stage == stage_b) %>%
        mutate(prop = n_cells / sum(n_cells)) %>%
        filter(annotation %in% common_after_norm) %>%
        select(annotation, prop)
      
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
      
      pair_weight <- prop_lookup[anno_i] * prop_lookup[anno_j]
      
      
      valid_weight <- is.finite(vec_a) &
        is.finite(vec_b) &
        is.finite(pair_weight) &
        pair_weight > 0
      
      vec_a <- vec_a[valid_weight]
      vec_b <- vec_b[valid_weight]
      pair_weight <- pair_weight[valid_weight]
      
      
      abs_change <- abs(
        vec_b - vec_a)
      
      
      matrix_cor <- weighted_spearman(
        vec_a,
        vec_b,
        pair_weight)
      
      
      pearson_cor <- weighted_cor(
        vec_a,
        vec_b,
        pair_weight)
      
      
      median_change <- weighted_median(
        abs_change,
        pair_weight)
      
      
      rmse <- sqrt(
        weighted.mean(
          (vec_b - vec_a)^2,
          pair_weight ))
      
    } else {
      
      valid <- is.finite(vec_a) &
        is.finite(vec_b)
      
      vec_a <- vec_a[valid]
      vec_b <- vec_b[valid]
      
      matrix_cor <- cor(
        vec_a,
        vec_b,
        method = "spearman")
      
      pearson_cor <- cor(
        vec_a,
        vec_b,
        method = "pearson")
      
      median_change <- median(
        abs(vec_b - vec_a))
      
      rmse <- sqrt(
        mean((vec_b - vec_a)^2)
      )}
      
      data.frame(
        stage_a = stage_a,
        stage_b = stage_b,
        comparison = paste0(stage_a, "-", stage_b),
        matrix_correlation = matrix_cor,
        pearson_correlation = pearson_cor,
        median_absolute_change = median_change,
        rmse_change = rmse)
       }
  
  
  result_list <- lapply(
    seq_len(nrow(stage_pairs)),
    function(i) {
      compare_two_stages(
        stage_pairs$stage_a[i],
        stage_pairs$stage_b[i]
      )
    }
  )
  
  bind_rows(result_list)
  
}
  


stage_pairs <- data.frame(
  stage_a = c(
    "E11.5", "E12.5", "E13.5", "E14.5",
    "E15.5", "E17.5", "P0", "P4"
  ),
  stage_b = c(
    "E12.5", "E13.5", "E14.5", "E15.5",
    "E17.5", "P0", "P4", "P7"
  )
)

relationship_summary <- compare_stage_relationships(
  obj = obj,
  stage_pairs = stage_pairs,
  stage_col = "stage",
  anno_col = "anno",
  min_cells = 3,
  weighted = TRUE
)

relationship_summary <- compare_stage_relationships(
  obj = obj,
  stage_pairs = stage_pairs,
  stage_col = "stage",
  anno_col = "anno",
  min_cells = 20,
  weighted = FALSE
)


relationship_summary

write.csv(
  relationship_summary,
  "results/inter-anno_relationship_weighted.csv",
  row.names = FALSE)

write.csv(
  relationship_summary,
  "results/inter-anno_relationship_unweighted.csv",
  row.names = FALSE)



base_theme <- theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

p_matrix_cor <- ggplot(
  relationship_summary,
  aes(
    x = comparison,
    y = matrix_correlation,
    group = 1
  )
) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  labs(
    x = NULL,
    y = "Matrix correlation"
  ) +
  base_theme

p_pearson_cor <- ggplot(
  relationship_summary,
  aes(
    x = comparison,
    y = pearson_correlation,
    group = 1
  )
) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  labs(
    x = NULL,
    y = "Pearson correlation"
  ) +
  base_theme

p_mean_change <- ggplot(
  relationship_summary,
  aes(
    x = comparison,
    y = mean_absolute_change,
    group = 1
  )
) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  labs(
    x = NULL,
    y = "Mean absolute change"
  ) +
  base_theme

p_median_change <- ggplot(
  relationship_summary,
  aes(
    x = comparison,
    y = median_absolute_change,
    group = 1
  )
) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  labs(
    x = NULL,
    y = "Median absolute change"
  ) +
  base_theme

p_rmse_change <- ggplot(
  relationship_summary,
  aes(
    x = comparison,
    y = rmse_change,
    group = 1
  )
) +
  geom_line(linewidth = 1) +
  geom_point(size = 3) +
  labs(
    x = NULL,
    y = "RMSE change"
  ) +
  base_theme

(p_median_change | p_rmse_change) / (p_matrix_cor | p_pearson_cor )

ggsave("figures/between_cluster_relationship_weighted.png", width = 6, height = 5)
ggsave("figures/between_cluster_relationship_unweighted.png", width = 5, height = 5)
