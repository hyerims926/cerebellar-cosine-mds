seed <- 1111
seed <- 2222
seed <- 3333
seed <- 4444
seed <- 1234

set.seed(seed)

pcs_use <- 1:10
pcs_use <- 1:30
pcs_use <- 1:50

obj <- readRDS("Rdata/obj_E10.5.RDS")
obj <- readRDS("Rdata/obj_E11.5.RDS")
obj <- readRDS("Rdata/obj_E12.5.RDS")
obj <- readRDS("Rdata/obj_E13.5.RDS")
obj <- readRDS("Rdata/obj_E14.5.RDS")
obj <- readRDS("Rdata/obj_E15.5.RDS")
obj <- readRDS("Rdata/obj_E17.5.RDS")
obj <- readRDS("Rdata/obj_P0.RDS")
obj <- readRDS("Rdata/obj_P4.RDS")
obj <- readRDS("Rdata/obj_P7.RDS")
obj <- readRDS("Rdata/obj_P14.RDS")
obj <- readRDS("Rdata/obj_adult.RDS")


obj <- NormalizeData(obj) 
obj <- FindVariableFeatures(obj, nfeatures = 3000)
obj <- ScaleData(obj, features = VariableFeatures(obj))
obj <- RunPCA(obj, features = VariableFeatures(obj), 
                    npcs = 50, seed.use = seed)


obj <- RunHarmony(obj, group.by.vars = "batch")    
obj <- RunUMAP(obj, reduction = "harmony", dims = pcs_use, 
                     reduction.name = "umap_harmony",reduction.key  = "UMAPH_", 
                     min.dist = 0.5, seed.use = seed)


DimPlot(obj, reduction = "umap_harmony", group.by = "anno", cols = cols) + 
  theme_classic() + ggtitle(NULL) + NoLegend() +
  theme_void() +
  annotate("segment", x = -14, xend = -10, y = -14, yend = -14,
           arrow = arrow(length = unit(0.15, "cm"))) +
  annotate("segment", x = -14, xend = -14, y = -14, yend = -10,
           arrow = arrow(length = unit(0.15, "cm"))) +
  annotate("text", x = -13, y = -15, label = "UMAP", size = 6) + NoLegend()

ggsave(paste0("figures/umaph_", unique(obj$author_stage), "_pc" , max(pcs_use),"_", seed, ".png"),
       width = 6, height = 6)


emb <- Embeddings(obj, "harmony")[, pcs_use]

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

# eig50 <- eigs_sym(B, k = 50, which = "LA")
# explained_approx <- eig$values / sum(eig50$values[eig50$values > 0])
# explained_approx
# sum(explained_approx)

mds <- eig$vectors %*% diag(sqrt(pmax(eig$values, 0)))
colnames(mds) <- c("MDS1", "MDS2")

obj$x <- mds[,1]
obj$x <- -mds[,1]
obj$y <- mds[,2]
obj$y <- -mds[,2]

df_plot <- data.frame(
  x = obj$x,
  y = obj$y,
  anno = obj$anno)

# xr <- max(abs(df_plot$polar_x), na.rm = TRUE)
# yr <- max(abs(df_plot$polar_y), na.rm = TRUE)
# lim <- max(xr, yr)
lim <- 1.2

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


ggsave(paste0("figures/mds_", unique(obj$author_stage), "_pc" , max(pcs_use),"_", seed, ".png"),
       width = 4, height = 4)

mds.emb <- as.matrix(
  obj@meta.data[, c("x", "y")])

colnames(mds.emb) <- c("MDS_1", "MDS_2")

obj[["mds"]] <- CreateDimReducObject(
  embeddings = mds.emb,
  key = "MDS_",
  assay = DefaultAssay(obj)
)


saveRDS(obj,paste0("Rdata/mds_", unique(obj$author_stage), "_pc" , max(pcs_use), ".RDS"))


obj <- readRDS("Rdata/mds_E10.5_pc50.RDS")
obj <- readRDS("Rdata/mds_E11.5_pc50.RDS")
obj <- readRDS("Rdata/mds_E12.5_pc50.RDS")
obj <- readRDS("Rdata/mds_E13.5_pc50.RDS")
obj <- readRDS("Rdata/mds_E14.5_pc50.RDS")
obj <- readRDS("Rdata/mds_E15.5_pc50.RDS")
obj <- readRDS("Rdata/mds_E17.5_pc50.RDS")
obj <- readRDS("Rdata/mds_P0_pc50.RDS")
obj <- readRDS("Rdata/mds_P4_pc50.RDS")
obj <- readRDS("Rdata/mds_P7_pc50.RDS")
obj <- readRDS("Rdata/mds_P14_pc50.RDS")
obj <- readRDS("Rdata/mds_adult (9 weeks)_pc50.RDS")


