seed<- 1234
set.seed(seed)
pcs_use <- 1:50


obj_17.5 <- readRDS("Rdata/mds_E17.5_pc50.RDS")
obj_P0 <- readRDS("Rdata/mds_P0_pc50.RDS")
obj_P4 <- readRDS("Rdata/mds_P4_pc50.RDS")

obj_17.5$stage <- "E17.5"
obj_P0$stage <- "P0"
obj_P4$stage <- "P4"

obj <- merge(
  x = obj_17.5,
  y = list(obj_P0, obj_P4),
  add.cell.ids = c("E17.5", "P0", "P4"),
  project = "dev_cerebellum")

obj_P0 <- subset(obj_P0, subset = !anno %in% c("Purkinje_defined", "Purkinje_diff", "Purkinje_maturing"))

obj_P0 <- subset(obj_P0, subset = !anno %in% c("GCP", "GC_diff_1", "GC_diff_2", "GC_defined"))

obj_P0 <- subset(obj_P0, subset = !anno %in% c("interneuron_diff", "interneuron_defined"))

obj_P0 <- subset(obj_P0, subset = !anno %in% c("glut_DN_defined", "isth_N_defined"))

obj_P0 <- subset(obj_P0, subset = !anno %in% "GABA_DN_defined")


obj_P0 <- NormalizeData(obj_P0) 
obj_P0 <- FindVariableFeatures(obj_P0, nfeatures = 3000)
obj_P0 <- ScaleData(obj_P0, features = VariableFeatures(obj_P0))
obj_P0 <- RunPCA(obj_P0, features = VariableFeatures(obj_P0), 
                        npcs = 50, seed.use = seed)

obj_P0 <- RunHarmony(obj_P0, group.by.vars = "batch")  


emb <- Embeddings(obj_P0, "harmony")[, pcs_use]

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

obj_P0$x <- mds[,1]
obj_P0$x <- -mds[,1]
obj_P0$y <- mds[,2]
obj_P0$y <- -mds[,2]

df_plot <- data.frame(
  x = obj_P0$x,
  y = obj_P0$y,
  anno = obj_P0$anno)

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

ggsave("figures/capped_purk_P0.png", width =4, height =4)
ggsave("figures/capped_gra_P0.png", width =4, height =4)
ggsave("figures/capped_inter_P0.png", width =4, height =4)
ggsave("figures/capped_glut_P0.png", width =4, height =4)
ggsave("figures/capped_gaba_P0.png", width =4, height =4)

obj_P0$stage = "capped_P0"

obj <- merge(
  x = obj,
  y = obj_P0,
  project = "dev_cerebellum")

saveRDS(obj, "Rdata/obj_purk_P0.RDS")
saveRDS(obj, "Rdata/obj_gra_P0.RDS")
saveRDS(obj, "Rdata/obj_inter_P0.RDS")
saveRDS(obj, "Rdata/obj_glut_P0.RDS")
saveRDS(obj, "Rdata/obj_gaba_P0.RDS")


df <- data.frame(
  cell = colnames(obj),
  stage = obj$stage,
  annotation = obj$anno,
  MDS_1 = obj$x,
  MDS_2 = obj$y)


write.csv(df, "mds-csv/mds_capped_purk_P0.csv")
write.csv(df, "mds-csv/mds_capped_gra_P0.csv")
write.csv(df, "mds-csv/mds_capped_inter_P0.csv")
write.csv(df, "mds-csv/mds_capped_glut_P0.csv")
write.csv(df, "mds-csv/mds_capped_gaba_P0.csv")


