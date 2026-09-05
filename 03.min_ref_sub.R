# obj from 11.merge.R

obj <- readRDS("Rdata/obj_all.RDS")

seed<- 1234
set.seed(seed)
pcs_use <- 1:50


get_downsample_flag <- function(obj_list, stage_pair) {
  base_stage <- stage_pair[1]
  target_stage <- stage_pair[2]
  
  obj <- subset(obj_list, subset = stage %in% c(base_stage, target_stage))
  
  meta <- obj@meta.data %>%
    mutate(cell = rownames(.))
  
  target_n <- meta %>%
    filter(stage == base_stage) %>%
    dplyr::count(anno, name = "n_target")
  
  flag_df <- meta %>%
    filter(stage == target_stage) %>%
    dplyr::count(anno, name = "n_avail") %>%
    inner_join(target_n, by = "anno") %>%
    mutate(downsampled = n_avail > n_target) %>%
    filter(n_target >= 20) %>%
    mutate(downsample_status = ifelse(downsampled, "TRUE", as.character(round(n_avail / n_target, 3))))
  
  return(flag_df)
}

merge_and_downsample_pairwise <- function(obj_list, stage_pair, project = "dev_cerebellum") {

  base_stage <- stage_pair[1]
  target_stage <- stage_pair[2]
  
  obj <- subset(obj_list, subset = stage %in% c(base_stage, target_stage))
  
  ## Meta
  meta <- obj@meta.data %>%
    mutate(cell = rownames(.))
  
  ## Target n
  target_n <- meta %>%
    filter(stage == base_stage) %>%
    dplyr::count(anno, name = "n_target")
  
  
  cells_sampled <- meta %>%
    filter(stage == target_stage) %>%
    inner_join(target_n, by = "anno") %>%
    group_by(anno) %>%
    group_modify(~ {
      n_sample <- min(nrow(.x), unique(.x$n_target))
      slice_sample(.x, n = n_sample)
    }) %>%
    ungroup() %>%
    pull(cell)
  

  obj_sub <- subset(obj, subset = stage == target_stage)
  obj_sub <- obj[, cells_sampled]
  
  return(obj_sub)
}

#stage_pair = c(reference, sampling)

stage_pair = c("E11.5", "E12.5")
stage_pair = c("E12.5", "E11.5")

stage_pair = c("E12.5", "E13.5")
stage_pair = c("E13.5", "E12.5")

stage_pair = c("E13.5", "E14.5")
stage_pair = c("E14.5", "E13.5")

stage_pair = c("E14.5", "E15.5")
stage_pair = c("E15.5", "E14.5")

stage_pair = c("E15.5", "E17.5")
stage_pair = c("E17.5", "E15.5")

stage_pair = c("E17.5", "P0")
stage_pair = c("P0", "E17.5")

stage_pair = c("P0", "P4")
stage_pair = c("P4", "P0")

stage_pair = c("P4", "P7")
stage_pair = c("P7", "P4")
# 
# stage_pair = c("P7", "P14")
# stage_pair = c("P14", "P7")
# 
# stage_pair = c("P14", "adult")
# stage_pair = c("adult","P14")

csv <- get_downsample_flag(obj,stage_pair)
csv

write.csv(csv, paste0("results/",stage_pair[2],"to",stage_pair[1],".csv"))

obj_reference <- subset(obj, subset = stage == stage_pair[1])

obj_sampling <- merge_and_downsample_pairwise(
  obj_list = obj,
  stage_pair = stage_pair)

df_plot <- data.frame(
  x = obj_sampling$x,
  y = obj_sampling$y,
  anno = obj_sampling$anno)

lim <- 1


ggplot(df_plot, aes(x, y)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey10", alpha = 0.4, linewidth = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey10", alpha = 0.4, linewidth = 0.7) +
  
  geom_point(aes(color = anno), size = 0.35, alpha = 0.5) +
  
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  coord_equal(xlim = c(-lim, lim), ylim = c(-lim, lim), expand = FALSE) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 1)) +
  theme_classic() +
  NoLegend() + NoAxes() +
  theme(panel.border = element_rect(color = "grey10", fill = NA, linewidth = 0.8))

ggsave(paste0("figures/min_", stage_pair[1],"_", stage_pair[2], "_before.png"), width = 4, height = 4)

obj_merge <- merge(
  x = obj_reference,
  y = list(obj_sampling),
  project = "dev_cerebellum")


obj_sampling <- NormalizeData(obj_sampling) 
obj_sampling <- FindVariableFeatures(obj_sampling, nfeatures = 3000)
obj_sampling <- ScaleData(obj_sampling, features = VariableFeatures(obj_sampling))
obj_sampling <- RunPCA(obj_sampling, features = VariableFeatures(obj_sampling), 
                       npcs = 50, seed.use = seed)

obj_sampling <- RunHarmony(obj_sampling, group.by.vars = "batch")  


emb <- Embeddings(obj_sampling, "harmony")[, pcs_use]

dist_normal <- sqrt(rowSums(emb^2))

vec <- sweep(emb, 2, 0, "-")

dir_unit <- vec / pmax(dist_normal)
cos_mat <- tcrossprod(dir_unit)
cos_mat <- pmin(pmax(cos_mat, -1), 1)
D2 <- (1 - cos_mat)^2

n <- nrow(D2)
rmean <- rowMeans(D2)
tmean <- mean(D2)

B <- -0.5 * (D2 - rmean - rep(rmean, each = n) + tmean)


eig <- eigs_sym(B, k = 2, which = "LA")
eig$values

mds <- eig$vectors %*% diag(sqrt(pmax(eig$values, 0)))
colnames(mds) <- c("MDS1", "MDS2")

obj_sampling$x <- mds[,1]
obj_sampling$x <- -mds[,1]
obj_sampling$y <- mds[,2]
obj_sampling$y <- -mds[,2]



df_plot <- data.frame(
  x = obj_sampling$x,
  y = obj_sampling$y,
  anno = obj_sampling$anno)

lim <- 1

ggplot(df_plot, aes(x, y)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey10", alpha = 0.4, linewidth = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey10", alpha = 0.4, linewidth = 0.7) +
  
  geom_point(aes(color = anno), size = 0.35, alpha = 0.5) +
  
  scale_color_manual(values = cols) +
  scale_fill_manual(values = cols) +
  coord_equal(xlim = c(-lim, lim), ylim = c(-lim, lim), expand = FALSE) +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1), ncol = 1)) +
  theme_classic() +
  NoLegend() + NoAxes() +
  theme(panel.border = element_rect(color = "grey10", fill = NA, linewidth = 0.8))

ggsave(paste0("figures/min_", stage_pair[1],"_", stage_pair[2], "_after.png"), width = 4, height = 4)


obj_sampling$stage <- paste0(stage_pair[1],stage_pair[2])

obj_merge <- merge(
  x = obj_merge,
  y = list(obj_sampling),
  project = "dev_cerebellum")


df <- data.frame(
  cell = colnames(obj_merge),
  stage = obj_merge$stage,
  annotation = obj_merge$anno,
  MDS_1 = obj_merge$x,
  MDS_2 = obj_merge$y,
  UMAP_1 = obj_merge$umap_harmony_1,
  UMAP_2 = obj_merge$umap_harmony_2)


write.csv(df, paste0("mds-csv/mds_sampling_", stage_pair[1], "_", stage_pair[2],".csv"))


