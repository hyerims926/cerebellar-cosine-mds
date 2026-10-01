seed <- 1234
set.seed(seed)
pcs_use <- 1:50
harmony_cols <- paste0("harmony_", 1:50)

#########unintegrated harmony ##############
# path <- ""
#df <- read.csv(paste0(path, "mds_all.csv"))
df<- read.csv("mds-csv/mds_all_harmony.csv")   
########################################

######integrated pca###########################################
obj <- readRDS("Rdata/obj_all.RDS")

obj <- NormalizeData(obj) 
obj <- FindVariableFeatures(obj, nfeatures = 3000)
obj <- ScaleData(obj, features = VariableFeatures(obj))
obj <- RunPCA(obj, features = VariableFeatures(obj), npcs = 50, seed.use = seed)

obj <- RunHarmony(obj, group.by.vars = "batch")

add_harmony <- function(obj, prefix = "harmony", dims = 1:50) {
  emb <- Embeddings(obj, reduction = prefix)
  for (i in dims) {
    obj[[paste0(prefix, "_", i)]] <- emb[, i]
  }
  return(obj)}

obj  <- add_harmony(obj)

harmony_cols <- paste0("harmony_", 1:50)

df <- data.frame(
  cell = colnames(obj),
  stage = obj$stage,
  annotation = obj$anno,
  MDS_1 = obj$x,
  MDS_2 = obj$y,
  obj@meta.data[, harmony_cols])

write.csv(df, "mds-csv/mds_all_integreted_harmony.csv")   

#############################################################


cell_count <- df %>%
  dplyr::count(stage, annotation, name = "n_cells")

cell_count

sampled_df <- df %>%
  filter(
    complete.cases(across(all_of(harmony_cols)))
  ) %>%
  group_by(stage, annotation) %>%
  filter(n() >= 100) %>%
  slice_sample(n = 100, replace = FALSE) %>%
  ungroup()

sampled_df %>%
  dplyr::count(stage, annotation) %>%
  arrange(stage, annotation)

compute_cosine_dispersion <- function(
    df,
    harmony_cols,
    stage_col = "stage",
    annotation_col = "annotation"
) {
  
  group_vars <- c(stage_col, annotation_col)
  
  cell_level <- df %>%
    group_by(across(all_of(group_vars))) %>%
    group_modify(~ {
      
      coords <- as.matrix(
        .x[, harmony_cols, drop = FALSE]
      )
      
      storage.mode(coords) <- "double"
      
      # annotation centroid
      centroid <- colMeans(coords)
      
      cell_norm <- sqrt(rowSums(coords^2))
      centroid_norm <- sqrt(sum(centroid^2))
      
      valid <- (
        is.finite(cell_norm) &
          cell_norm > 0 &
          is.finite(centroid_norm) &
          centroid_norm > 0
      )
      
      coords_valid <- coords[valid, , drop = FALSE]
      cell_norm_valid <- cell_norm[valid]
      
      cosine_similarity <- as.vector(
        coords_valid %*% centroid
      ) / (
        cell_norm_valid * centroid_norm
      )
      
      cosine_similarity <- pmax(
        -1,
        pmin(1, cosine_similarity))
      
      cosine_dispersion <- 1 - cosine_similarity
      
      data.frame(
        cell = .x$cell[valid],
        cosine_dispersion = cosine_dispersion
      )
    }) %>%
    ungroup()
  
  summary_level <- cell_level %>%
    group_by(across(all_of(group_vars))) %>%
    summarise(
      n_cells = n(),
      median_dispersion = median(
        cosine_dispersion,
        na.rm = TRUE
      ),
      mean_dispersion = mean(
        cosine_dispersion,
        na.rm = TRUE
      ),
      q25 = quantile(
        cosine_dispersion,
        0.25,
        na.rm = TRUE
      ),
      q75 = quantile(
        cosine_dispersion,
        0.75,
        na.rm = TRUE
      ),
      IQR_dispersion = IQR(
        cosine_dispersion,
        na.rm = TRUE
      ),
      .groups = "drop"
    )
  
  list(
    cell_level = cell_level,
    summary_level = summary_level
  )
}

disp_result <- compute_cosine_dispersion(
  df = sampled_df,
  harmony_cols = harmony_cols)

disp_summary <- disp_result$summary_level


stage_levels <- c(
  "E11.5", "E12.5", "E13.5", "E14.5",
  "E15.5", "E17.5", "P0", "P4")

disp_plot_df <- disp_summary %>%
  mutate(stage = factor(stage, levels = stage_levels)) %>%
  filter(!is.na(stage))

ggplot(
  disp_plot_df,
  aes(
    x = stage,
    y = median_dispersion,
    group = annotation,
    color = annotation)) +
  geom_line(
    linewidth = 0.8,
    alpha = 0.8) +
  scale_color_manual(values = cols) +
  geom_point(
    size = 2) +
  labs(
    x = NULL,
    y = "Median cosine dispersion",
    color = "Annotation",
    title = "Within-annotation dispersion across development") +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1),
    legend.position = "right")

annotation_use <- c(
  "progenitor",
  "GCP",
  "GC_diff_1",
  "GC_diff_2",
  "GC_defined",
  "Purkinje_diff",
  "Purkinje_defined",
  "Purkinje_maturing",
  "glioblast"
)

disp_plot_df %>%
  filter(annotation %in% annotation_use) %>%
  ggplot(
    aes(
      x = stage,
      y = median_dispersion,
      group = annotation,
      color = annotation
    )
  ) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = cols) +
  geom_point(size = 2.5) +
  labs(
    x = NULL,
    y = "Median cosine dispersion",
    color = "Annotation"
    # title = "Within-annotation dispersion across development"
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  ) + NoLegend()

ggsave("figures/dispersion_within-cluster_integrated.png", width = 5, height = 4)
