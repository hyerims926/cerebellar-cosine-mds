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


### purk
sce <- as.SingleCellExperiment(obj)

mds <- Embeddings(obj, "mds")
reducedDim(sce, "MDS") <- mds

umaph <- Embeddings(obj, "umap_harmony")
reducedDim(sce, "UMAPH") <- umaph

keep_lineage <- c(
  "progenitor",
  "VZ_neuroblast_1",
  "VZ_neuroblast_2",
  "VZ_neuroblast_3",
  "Purkinje_diff",
  "Purkinje_defined"
)

sce_purk <- sce[, as.character(sce$anno) %in% keep_lineage]
sce_purk$anno <- factor(as.character(sce_purk$anno))

sce_purk <- slingshot(
  sce_purk,
  clusterLabels = "anno",
  reducedDim = "MDS",
  start.clus = "progenitor",
  dist.method = "simple")

slingLineages(sce_purk)

pt1 <- slingPseudotime(sce_purk)[, 1]
names(pt1) <- colnames(sce_purk)

obj$pt_mds_purk <- NA_real_
obj$pt_mds_purk[names(pt1)] <- pt1

sce_purk <- slingshot(
  sce_purk,
  clusterLabels = "anno",
  reducedDim = "UMAPH",
  start.clus = "progenitor",
  dist.method = "simple")

slingLineages(sce_purk)

pt2 <- slingPseudotime(sce_purk)[, 1]
names(pt2) <- colnames(sce_purk)

obj$pt_umaph_purk <- NA_real_
obj$pt_umaph_purk[names(pt2)] <- pt2

FeaturePlot(obj, features = "pt_mds_purk", reduction = "mds") + ggtitle(NULL)+
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", alpha = 0.7) +NoLegend()

ggsave("figures/pt_mds_purk.png", width = 4, height = 4)

FeaturePlot(obj, features = "pt_umaph_purk", reduction = "umap_harmony") + ggtitle(NULL)+
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", alpha = 0.7) +NoLegend()

ggsave("figures/pt_umaph_purk.png", width = 4, height = 4)

### gra
sce <- as.SingleCellExperiment(obj)

mds <- Embeddings(obj, "mds")
reducedDim(sce, "MDS") <- mds

umaph <- Embeddings(obj, "umap_harmony")
reducedDim(sce, "UMAPH") <- umaph

keep_lineage <- c(
  "progenitor",
  "GCP",
  "GC_diff_1",
  "GC_diff_2",
  "GC_defined")

sce_gra <- sce[, as.character(sce$anno) %in% keep_lineage]
sce_gra$anno <- factor(as.character(sce_gra$anno))

sce_gra <- slingshot(
  sce_gra,
  clusterLabels = "anno",
  reducedDim = "MDS",
  start.clus = "progenitor",
  dist.method = "simple")

slingLineages(sce_gra)
# [1] "progenitor" "GCP"        "GC_diff_1"  "GC_diff_2" 

pt1 <- slingPseudotime(sce_gra)[, 1]
names(pt1) <- colnames(sce_gra)

obj$pt_mds_gra <- NA_real_
obj$pt_mds_gra[names(pt1)] <- pt1

sce_gra <- slingshot(
  sce_gra,
  clusterLabels = "anno",
  reducedDim = "UMAPH",
  start.clus = "progenitor",
  dist.method = "simple")

slingLineages(sce_gra)

pt2 <- slingPseudotime(sce_gra)[, 1]
names(pt2) <- colnames(sce_gra)

obj$pt_umaph_gra <- NA_real_
obj$pt_umaph_gra[names(pt2)] <- pt2

FeaturePlot(obj, features = "pt_mds_gra", reduction = "mds") + ggtitle(NULL)+
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", alpha = 0.7) +NoLegend()

ggsave("figures/pt_mds_gra.png", width = 4, height = 4)

FeaturePlot(obj, features = "pt_umaph_gra", reduction = "umap_harmony") + ggtitle(NULL)+
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", alpha = 0.7) +NoLegend()

ggsave("figures/pt_umaph_gra.png", width = 4, height = 4)

### glut
sce <- as.SingleCellExperiment(obj)

mds <- Embeddings(obj, "mds")
reducedDim(sce, "MDS") <- mds

umaph <- Embeddings(obj, "umap_harmony")
reducedDim(sce, "UMAPH") <- umaph

keep_lineage <- c(
  "progenitor",
  "NTZ_neuroblast_1",
  "NTZ_neuroblast_2",
  "NTZ_neuroblast_3",
  "glut_DN_defined"
)

sce_ntz_glut <- sce[, as.character(sce$anno) %in% keep_lineage]
sce_ntz_glut$anno <- factor(as.character(sce_ntz_glut$anno))

sce_ntz_glut <- slingshot(
  sce_ntz_glut,
  clusterLabels = "anno",
  reducedDim = "MDS",
  start.clus = "progenitor",
  dist.method = "simple")

slingLineages(sce_ntz_glut)
# [1] "progenitor"       "NTZ_neuroblast_1" "NTZ_neuroblast_2" "NTZ_neuroblast_3" "glut_DN_defined" 
# [1] "progenitor"       "NTZ_neuroblast_2" "glut_DN_defined" 

pt1 <- slingPseudotime(sce_ntz_glut)[, 1]
names(pt1) <- colnames(sce_ntz_glut)

obj$pt_mds_ntz_glut <- NA_real_
obj$pt_mds_ntz_glut[names(pt1)] <- pt1

sce_ntz_glut <- slingshot(
  sce_ntz_glut,
  clusterLabels = "anno",
  reducedDim = "UMAPH",
  start.clus = "progenitor",
  dist.method = "simple")

slingLineages(sce_ntz_glut)

pt2 <- slingPseudotime(sce_ntz_glut)[, 1]
names(pt2) <- colnames(sce_ntz_glut)

obj$pt_umaph_ntz_glut <- NA_real_
obj$pt_umaph_ntz_glut[names(pt2)] <- pt2


FeaturePlot(obj, features = "pt_mds_ntz_glut", reduction = "mds") + ggtitle(NULL)+
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", alpha = 0.7) +NoLegend()

FeaturePlot(obj, features = "pt_umaph_ntz_glut", reduction = "umap_harmony") + ggtitle(NULL)+
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", alpha = 0.7) +NoLegend()


### interneuron
sce <- as.SingleCellExperiment(obj)

mds <- Embeddings(obj, "mds")
reducedDim(sce, "MDS") <- mds

umaph <- Embeddings(obj, "umap_harmony")
reducedDim(sce, "UMAPH") <- umaph

keep_lineage <- c(
  "progenitor",
  "interneuron_defined",
  "interneuron_diff"
)

sce_interneuron <- sce[, as.character(sce$anno) %in% keep_lineage]
sce_interneuron$anno <- factor(as.character(sce_interneuron$anno))

sce_interneuron <- slingshot(
  sce_interneuron,
  clusterLabels = "anno",
  reducedDim = "MDS",
  start.clus = "progenitor",
  dist.method = "simple")

slingLineages(sce_interneuron)

pt1 <- slingPseudotime(sce_interneuron)[, 1]
names(pt1) <- colnames(sce_interneuron)

obj$pt_mds_interneuron <- NA_real_
obj$pt_mds_interneuron[names(pt1)] <- pt1

sce_interneuron <- slingshot(
  sce_interneuron,
  clusterLabels = "anno",
  reducedDim = "UMAPH",
  start.clus = "progenitor",
  dist.method = "simple")

slingLineages(sce_interneuron)

pt2 <- slingPseudotime(sce_interneuron)[, 1]
names(pt2) <- colnames(sce_interneuron)

obj$pt_umaph_interneuron <- NA_real_
obj$pt_umaph_interneuron[names(pt2)] <- pt2

FeaturePlot(obj, features = "pt_mds_interneuron", reduction = "mds") + ggtitle(NULL) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", alpha = 0.7) +NoLegend()

FeaturePlot(obj, features = "pt_umaph_interneuron", reduction = "umap_harmony") + ggtitle(NULL) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey50", alpha = 0.7) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", alpha = 0.7) +NoLegend()

saveRDS(obj,paste0("Rdata/mds_slingshot_", unique(obj$author_stage), ".RDS"))
