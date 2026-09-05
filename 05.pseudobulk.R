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

cosine_dissimilarity <- function(x, y) {
  
  valid <- is.finite(x) & is.finite(y)
  
  x <- x[valid]
  y <- y[valid]
  
  norm_x <- sqrt(sum(x^2))
  norm_y <- sqrt(sum(y^2))
  
  
  cos_sim <- sum(x * y) / (norm_x * norm_y)
  cos_sim <- pmax(-1, pmin(1, cos_sim))
  
  1 - cos_sim
}

annotations <- unique(
  sub(
    "___.*$",
    "",
    colnames(pseudobulk_mat)
  )
)

shift_all <- lapply(
  seq_len(nrow(stage_pairs)),
  function(i) {
    
    stage_a <- stage_pairs$stage_a[i]
    stage_b <- stage_pairs$stage_b[i]
    comparison <- stage_pairs$comparison[i]
    
    bind_rows(
      lapply(
        annotations,
        function(anno) {
          
          col_a <- paste0(anno, "___", stage_a)
          col_b <- paste0(anno, "___", stage_b)
          
          if (!all(
            c(col_a, col_b) %in%
            colnames(pseudobulk_mat)
          )) {
            return(NULL)
          }
          
          data.frame(
            annotation = anno,
            stage_a = stage_a,
            stage_b = stage_b,
            comparison = comparison,
            cosine_dissimilarity =
              cosine_dissimilarity(
                pseudobulk_mat[, col_a],
                pseudobulk_mat[, col_b]
              )
          )
        }
      )
    )
  }
) %>%
  bind_rows()

shift_all <- shift_all %>%
  mutate(
    comparison = factor(
      comparison,
      levels = stage_pairs$comparison
    )
  )

anno_levels <- c(
  # Progenitor
  "progenitor",
  "VZ_neuroblast_1",
  "VZ_neuroblast_2",
  "VZ_neuroblast_3",
  "NTZ_neuroblast_1",
  "NTZ_neuroblast_2",

  # Granule lineage
  "GCP",
  "GC_diff_1",
  "GC_diff_2",
  "GC_defined",
  "GCP/UBCP",
  "GC/UBC_diff",
  "UBC_diff",
  "UBC_defined",
  
  # Purkinje lineage
  "Purkinje_diff",
  "Purkinje_defined",
  "Purkinje_maturing",
  
  # Other neurons
  "interneuron_diff",
  "interneuron_defined",
  "glut_DN_defined",
  "GABA_DN_defined",
  "isth_N_diff",
  "isth_N_defined",
  
  # Glia
  "glioblast",
  "oligo_progenitor"
)

shift_plot <- shift_all |>
  dplyr::filter(annotation %in% anno_levels)

shift_plot$annotation <- factor(
  shift_plot$annotation,
  levels = anno_levels
)

ggplot(
  shift_plot,
  aes(
    x = comparison,
    y = cosine_dissimilarity,
    group = annotation,
    color = annotation
  )
) +
  geom_line(
    linewidth = 0.9,
    alpha = 0.8
  ) +
  geom_point(
    size = 2.1
  ) +
  scale_color_manual(
    values = cols
  ) +
  labs(
    x = NULL,
    y = "Pseudobulk cosine dissimilarity",
    color = "Annotation"
    # title = "Annotation-specific transcriptomic shift across development"
  ) +
  theme_classic() +
  scale_y_continuous(limits = c(0, 0.3))+
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),
    legend.position = "right"
  ) +NoLegend()

ggsave("figures/pseudobulk_all.png", width = 5, height = 5)

ggplot(
  shift_plot,
  aes(
    x = comparison,
    y = cosine_dissimilarity,
    group = 1,
    color = annotation
  )
) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 1.8) +
  scale_color_manual(
    values = cols,
    guide = "none"
  ) +
  facet_wrap(
    ~ annotation,
    ncol = 5,
    scales = "free_y",
    drop = FALSE
  )+
  labs(
    x = NULL,
    y = "Pseudobulk cosine dissimilarity"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      size = 7
    ),
    strip.text = element_text(size = 8)
  )

ggsave("figures/pseudobulk_sep.png", width = 8, height = 8)


shift_summary <- shift_all %>%
  group_by(comparison) %>%
  summarise(
    median_pseudo_dissim = median(
      cosine_dissimilarity,
      na.rm = TRUE
    ),
    mean_pseudo_dissim = mean(
      cosine_dissimilarity,
      na.rm = TRUE
    ),
    n_annotation = sum(
      !is.na(cosine_dissimilarity)
    ),
    .groups = "drop"
  )

s_color_df <- read.csv(
  "results/s_color_by_stage_pair.csv",
  stringsAsFactors = FALSE
)

compare_df <- shift_summary %>%
  inner_join(
    s_color_df,
    by = "comparison"
  ) %>%
  mutate(
    landscape_change = 1 - S_color
  )


df_no_p0p4 <- compare_df %>%
  filter(comparison != "P0-P4")


fit_all <- lm(
  landscape_change ~ median_pseudo_dissim,
  data = compare_df
)

fit_no_p0p4 <- lm(
  landscape_change ~ median_pseudo_dissim,
  data = df_no_p0p4
)

r_all <- cor(
  compare_df$landscape_change,
  compare_df$median_pseudo_dissim,
  use = "complete.obs"
)

r_no_p0p4 <- cor(
  df_no_p0p4$landscape_change,
  df_no_p0p4$median_pseudo_dissim,
  use = "complete.obs"
)

r_all
r_no_p0p4

ggplot(
  compare_df,
  aes(
    x = median_pseudo_dissim,
    y = landscape_change
  )
) +
  
  geom_smooth(
    data = compare_df,
    aes(color = "Including P0-P4"),
    method = "lm",
    se = FALSE,
    linetype = "dashed"
  ) +
  
  geom_smooth(
    data = df_no_p0p4,
    aes(color = "Excluding P0-P4"),
    method = "lm",
    se = FALSE,
    linetype = "solid"
  ) +
  
  geom_point(size = 3) +
  
  ggrepel::geom_text_repel(
    aes(label = comparison),
    size = 4
  ) +
  annotate(
    "text",
    x = 0.12, y = 0.30,
    label = paste0("Excl. P0-P4: R = ", round(r_no_p0p4, 2)),
    size = 4,
    color = "red"
  ) +
  annotate(
    "text",
    x = 0.15, y = 0.24,
    label = paste0("All: R = ", round(r_all, 2)),
    size = 4,
    color = "blue"
  ) + 
  theme_classic() +
  
  labs(
    x = "Median pseudobulk cosine dissimilarity",
    y = "1 − SSIM"
  ) +
  
  NoLegend()

ggsave("figures/shift.png", width = 4, height = 4)
