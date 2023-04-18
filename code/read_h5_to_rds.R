library(Seurat)
#After download datasets from Bdbag
sc10v3a  <- Read10X_h5(".../sc10v3a/umi_counts.h5")
qccell = read.csv('.../sc10v3a/QC.csv')
metadata = read.csv('.../sc10v3a/sample_metadata.csv')
sc10v3aa <- CreateSeuratObject(sc10v3a,project = "sc10v3a") 
colnames(sc10v3a) <- metadata[,1]
sc10v3a <- sc10v3a[,qccell[,2]]
saveRDS(object = sc10v3a, file = paste0(".../sc10v3a/sc10v3a_dgc.rds"))

