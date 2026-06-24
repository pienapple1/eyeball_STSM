###Author Jinpei Lin
###Times 2024.07

library(Seurat)
library(loomR)

drds <- "/00data/"
dcoor <- "/00df_coor/"
dex_df <- "/01df_ex/"
nvg = 5000
for(s in Sys.glob(file.path(drds,"*.RDS"))){
    key <- gsub(".RDS","",basename(s))
    rds <- paste0(drds,key,".RDS")
    sp <- readRDS(rds)
    object <- NormalizeData(sp,normalization.method = "LogNormalize")
    pro_name <- key
    position <- data.frame(NAME = rownames(object@meta.data),
                       X = object@meta.data$x,
                       Y = object@meta.data$y)
    write.csv(position, file.path(dcoor, paste0(pro_name, "_cellpos.csv")), row.names = F)
    csv_fname <- NULL
    loom_fname <- NULL
    expmtx <- NULL
    object@active.assay <- names(object@assays)[1]
    object <- Seurat::FindVariableFeatures(object, nfeatures = nvg)
    spatial_mtx <- Seurat::GetAssayData(object)
    var_genes <- object@assays[[object@active.assay]]@var.features
    lf <- loomR::create(file.path(dex_df, paste0(pro_name, "_expmtx_nvg", nvg, ".loom")), spatial_mtx[var_genes, ], overwrite = T)
    lf$close_all()
    loom_fname <- file.path(dex_df, paste0(pro_name, "_expmtx_nvg", nvg, ".loom"))
    system(paste0("chmod 755 ", loom_fname))
}