library(Seurat)
mer <- readRDS("data/spatial/osmfish.rds")

sc <- readRDS("data/scrna/zeisel.rds")
zei_g = rownames(sc)
rm(sc)
gc()
length(zei_g)

###
sc <- readRDS("data/scrna/allen_SSp.rds")
ssp_g = rownames(sc)
rm(sc)
gc()
length(ssp_g)

###
sc <- readRDS("data/scrna/allen_Visp.rds")
visp_g = rownames(sc)
rm(sc)
gc()
length(visp_g)

###
common_scg = Reduce(intersect,list(zei_g,
                                   ssp_g,
                                   visp_g))
length(common_scg)

##############################################
sc <- readRDS("data/scrna/zeisel.rds")
sc <- sc[common_scg,]
sc <- CreateSeuratObject(sc,project = "sc_use")
sc <- NormalizeData(object = sc,normalization.method = "LogNormalize")
sc <- FindVariableFeatures(sc,nfeatures = 2000)
hvg_zei <- sc@assays$RNA@var.features
rm(sc)
gc()

#######################3
sc <- readRDS("data/scrna/allen_SSp.rds")
sc <- sc[common_scg,]
dim(sc)
sc <- CreateSeuratObject(sc,project = "sc_use")
sc <- NormalizeData(object = sc,normalization.method = "LogNormalize")
sc <- FindVariableFeatures(sc,nfeatures = 2000)
hvg_ssp <- sc@assays$RNA@var.features
rm(sc)
gc()

###############################
sc <- readRDS("data/scrna/allen_Visp.rds")
sc <- sc[common_scg,]
dim(sc)
sc <- CreateSeuratObject(sc,project = "sc_use")
sc <- NormalizeData(object = sc,normalization.method = "LogNormalize")
sc <- FindVariableFeatures(sc,nfeatures = 2000)
hvg_visp <- sc@assays$RNA@var.features
rm(sc)
gc()
###########################################
common_hvg = Reduce(union,list(hvg_zei,
                               hvg_ssp,
                               hvg_visp))
length(common_hvg)

owngene <- setdiff(common_hvg,rownames(mer))
length(owngene)

write.csv(owngene,file = '.../osm_hvg.csv') 

