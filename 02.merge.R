obj_10.5 <- readRDS("Rdata/mds_E10.5_pc50.RDS")
obj_11.5 <- readRDS("Rdata/mds_E11.5_pc50.RDS")
obj_12.5 <- readRDS("Rdata/mds_E12.5_pc50.RDS")
obj_13.5 <- readRDS("Rdata/mds_E13.5_pc50.RDS")
obj_14.5 <- readRDS("Rdata/mds_E14.5_pc50.RDS")
obj_15.5 <- readRDS("Rdata/mds_E15.5_pc50.RDS")
obj_17.5 <- readRDS("Rdata/mds_E17.5_pc50.RDS")
obj_P0 <- readRDS("Rdata/mds_P0_pc50.RDS")
obj_P4 <- readRDS("Rdata/mds_P4_pc50.RDS")
obj_P7 <- readRDS("Rdata/mds_P7_pc50.RDS")
obj_P14 <- readRDS("Rdata/mds_P14_pc50.RDS")
obj_adult <- readRDS("Rdata/mds_adult (9 weeks)_pc50.RDS")

obj_10.5$stage <- "E10.5"
obj_11.5$stage <- "E11.5"
obj_12.5$stage <- "E12.5"
obj_13.5$stage <- "E13.5"
obj_14.5$stage <- "E14.5"
obj_15.5$stage <- "E15.5"
obj_17.5$stage <- "E17.5"
obj_P0$stage <- "P0"
obj_P4$stage <- "P4"
obj_P4$stage <- "P4"
obj_P7$stage <- "P7"
obj_P14$stage <- "P14"
obj_adult$stage <- "adult"

add_umap_harmony_meta <- function(obj, prefix = "umap_harmony") {
  emb <- Embeddings(obj, reduction = "umap_harmony")
  
  obj[[paste0(prefix, "_1")]] <- emb[, 1]
  obj[[paste0(prefix, "_2")]] <- emb[, 2]
  
  return(obj)
}

add_pca <- function(obj, prefix = "pca", dims = 1:50) {
  emb <- Embeddings(obj, reduction = prefix)
  for (i in dims) {
    obj[[paste0(prefix, "_", i)]] <- emb[, i]
  }
  return(obj)
}

obj_10.5 <- add_umap_harmony_meta(obj_10.5)
obj_11.5 <- add_umap_harmony_meta(obj_11.5)
obj_12.5 <- add_umap_harmony_meta(obj_12.5)
obj_13.5 <- add_umap_harmony_meta(obj_13.5)
obj_14.5 <- add_umap_harmony_meta(obj_14.5)
obj_15.5 <- add_umap_harmony_meta(obj_15.5)
obj_17.5 <- add_umap_harmony_meta(obj_17.5)
obj_P0   <- add_umap_harmony_meta(obj_P0)
obj_P4   <- add_umap_harmony_meta(obj_P4)
obj_P7   <- add_umap_harmony_meta(obj_P7)
obj_P14   <- add_umap_harmony_meta(obj_P14)
obj_adult   <- add_umap_harmony_meta(obj_adult)


obj_10.5 <- add_pca(obj_10.5, prefix = "harmony")
obj_11.5 <- add_pca(obj_11.5, prefix = "harmony")
obj_12.5 <- add_pca(obj_12.5, prefix = "harmony")
obj_13.5 <- add_pca(obj_13.5, prefix = "harmony")
obj_14.5 <- add_pca(obj_14.5, prefix = "harmony")
obj_15.5 <- add_pca(obj_15.5, prefix = "harmony")
obj_17.5 <- add_pca(obj_17.5, prefix = "harmony")
obj_P0   <- add_pca(obj_P0, prefix = "harmony")
obj_P4   <- add_pca(obj_P4, prefix = "harmony")
obj_P7   <- add_pca(obj_P7, prefix = "harmony")
obj_P14   <- add_pca(obj_P14, prefix = "harmony")
obj_adult   <- add_pca(obj_adult, prefix = "harmony")

obj <- merge(
  x = obj_10.5,
  y = list(obj_11.5 ,obj_12.5, obj_13.5, obj_14.5, obj_15.5, obj_17.5, obj_P0, obj_P4, obj_P7, obj_P14, obj_adult),
  add.cell.ids = c("E10.5", "E11.5","E12.5", "E13.5", "E14.5","E15.5","E17.5","P0", "P4","P7","P14", "adult"),
  project = "dev_cerebellum"
)

saveRDS(obj, "Rdata/obj_all.RDS")
obj<- readRDS( "Rdata/obj_all.RDS")

harmony_cols <- paste0("harmony_", 1:50)

df <- data.frame(
  cell = colnames(obj),
  stage = obj$stage,
  annotation = obj$anno,
  MDS_1 = obj$x,
  MDS_2 = obj$y,
  UMAP_1 = obj$umap_harmony_1,
  UMAP_2 = obj$umap_harmony_2,
  obj@meta.data[, harmony_cols])

# df <- df %>%
#   filter(!is.na(annotation))

save(df, file = "Rdata/df_all.RData")

write.csv(df, "mds-csv/mds_all_harmony.csv")
