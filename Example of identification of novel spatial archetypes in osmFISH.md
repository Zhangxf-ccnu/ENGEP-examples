# Identification of novel spatial archetypes in osmFISH

This example use ENGEP to predict unmeasured genes for osmFISH and detect novel spatial archetypes for osmFISH. The datasets used in this examples can be downloaded from Zenodo (https://doi.org/10.5281/zenodo.7825993)

#### Get the predicted expression levels of unmeasured genes

Load ENGEP package.

```
library(ENGEP)
```

Read the spatial dataset. 

```
osmfish <- readRDS(".../osmfish.rds")
```

Read multiple reference datasets.

```
#read the first reference
zeisel <- readRDS(".../zeisel.rds")

#read the second reference
allen_SSp <- readRDS(".../allen_SSp.rds")

#read the third reference
allen_Visp <- readRDS(".../allen_Visp.rds")
```

These multiple reference datasets should be put into a list as an input to ENGEP.

```
ref_list = list(zeisel,allen_SSp,allen_Visp)
```

Here, we predict  3433 unmeasured genes ( the union of 2000 highly variable genes in each reference dataset) for osmFISH. The code to find these 3433 unmeasured genes is also in this repository.

```
pre_genes = c('Pvrl3','Wfs1','Cux2')
common_hvg <- read.csv('.../osm_hvg.csv',row.names=1)
pre_genes <- common_hvg$x
```

The predicted results can be obtained by running this function.

```
engep_predict <- engep_predict(osmfish,ref_list,pre_genes)
```

We can plot the expression patterns of these genes. Here, we plot expressions of three genes as an example.

```
locations=read.csv('.../osmloc.csv',row.names = 1)
gene_expressions = as.data.frame(engep_predict)
plot_genes = c('Pvrl3','Wfs1','Cux2')
gene_expressions = gene_expressions[,plot_genes]
locations = locations[rownames(gene_expressions),]

plot_exp_pattern(ncol(gene_expressions),locations,gene_expressions)
```
![image](https://github.com/Zhangxf-ccnu/ENGEP-examples/blob/main/docs/genes_of_osmFISH.png)

#### Identification of known spatial archetypes

First, we use MERINGUE to determine the spatially variable genes.

```
suppressMessages(library(MERINGUE))
cell_loc = read.csv('.../osmloc.csv',row.names = 1)
cell_loc = cell_loc[colnames(osmfish),]
mat <- normalizeCounts(counts = osmfish,
                       log=FALSE,
                       verbose=TRUE)
w <- getSpatialNeighbors(cell_loc,verbose=TRUE)
test <- rownames(mat)
I <- getSpatialPatterns(mat[test,], w)
results.filter <- filterSpatialPatterns(mat = mat[test,],
                                        I = I,
                                        w = w,
                                        adjustPv = TRUE,
                                        alpha = 0.05,
                                        minPercentCells = 0.05,
                                        verbose = TRUE)
```

Next, we perform hierarchical clustering on the spatially variable gene expression data using 1-Pearson correlation as the gene-wise distance metric. 

```
osm<- CreateSeuratObject(osmfish,project = "sc_use",min.cells = 0, 
                         min.features = 0)
osm <- NormalizeData(object = osm,normalization.method = "LogNormalize")
osm<- ScaleData(object = osm)

osm_fil <- as.matrix(osm@assays$RNA@scale.data[resultsfilter,])

corMatt <- stats::cor(as.matrix(t(osm_fil)),
                      as.matrix(t(osm_fil)), method = "pearson")

dd <- 1-corMatt
dd<-as.dist(dd)

hc <- hclust(dd, method = "ward.D")
```

 The number of archetypes is determined through visual inspection of the cluster dendrogram.

```
plot(hc,hang = -1)
out.id = cutree(hc,k=6)
```

In this paper, we identify six known archetypes in osmFISH. We can plot the expression patterns of every known archetype.

```
arche_exp <- get_expression(out.id,osm@assays$RNA@scale.data[resultsfilter,])

arche_exp_log <- get_expression(out.id,osm@assays$RNA@data[resultsfilter,])
gene_expressions = t(as.data.frame(arche_exp_log))
locations = as.data.frame(cell_loc)
plot_exp_pattern(ncol(gene_expressions),locations,gene_expressions,n_col=3)
```

![image](https://github.com/Zhangxf-ccnu/ENGEP-examples/blob/main/docs/osm_known.png)

#### Identification of novel spatial archetypes

First, we use MERINGUE to determine the spatially variable genes in the predicted genes.

```
pre_un = t(as.matrix(engep_predict))
counts <- cleanCounts(counts = pre_un,plot=FALSE,verbose=TRUE)
cell_loc = read.csv('.../osmloc.csv',row.names = 1)
cell_loc <- cell_loc[colnames(counts),]
mat <- normalizeCounts(counts = counts,
                       log=FALSE,
                       verbose=TRUE)
w <- getSpatialNeighbors(cell_loc,verbose=TRUE)
test <- rownames(mat)
I <- getSpatialPatterns(mat[test,], w)
results.filter.un <- filterSpatialPatterns(mat = mat[test,],
                                        I = I,
                                        w = w,
                                        adjustPv = TRUE,
                                        alpha = 0.05,
                                        minPercentCells = 0.05,
                                        verbose = TRUE)
```

We then calculate similarity between the predicted unmeasured genes and the known archetypes.

```
pre_un <- CreateSeuratObject(pre_un,project = "sc_use",min.cells = 0,min.features = 0)
pre_un <- NormalizeData(object = pre_un,normalization.method = "LogNormalize")
pre_un <- ScaleData(object = pre_un)
                    
corMat <- stats::cor(as.matrix(t(as.matrix(pre_un@assays$RNA@scale.data)[resultsfilter_un,]),as.matrix(t(arche_exp)), method = "pearson")
```

We use Pearson correlation as a metric to calculate the likelihood score of the unmeasured gene.  Then, we apply Hartigan’s dip test to analyze the distribution of likelihood scores of all unmeasured genes in order to detect any novel archetypes. 

```
maxcor <- apply(corMat, 1, max)
diptest::dip.test(maxcor, simulate.p.value = FALSE, B = 2000)
#Hartigans' dip test for unimodality / multimodality
#data:  maxcor
#D = 0.017995, p-value = 0.0003566
#alternative hypothesis: non-unimodal, i.e., at least bimodal
```

We use manifold mixup to partition the unmeasured genes into two groups.

```
c <- apply(as.matrix(pre_un@assays$RNA@scale.data[resultsfilter_un,],
           1,
           partition_threshold)
diag(c) <- 0
threshold <- sum(c)/((dim(pre_un)[1])*(dim(pre_un)[1]-1))
plot_score(maxcor,threshold=0.3933951)
```

![image](https://github.com/Zhangxf-ccnu/ENGEP-examples/blob/main/docs/osm_scores.png)

For genes belonging to known archetypes, we assign them to the closest matching known archetype.

```
similar_genes = names(maxcor)[which(maxcor>=threshold)]
corMat_similar = corMat[,similar_genes ]

cellidx = apply(corMat_similar ,2,which.max)

arche_exp_n <- get_expression(cellidx,pre_un@assays$RNA@scale.data[resultsfilter_un,])

arche_exp_n_log <- get_expression(cellidx,pre_un@assays$RNA@data[resultsfilter_un,])
gene_expressions <- t(as.data.frame(arche_exp_n_log))
plot_exp_pattern(ncol(gene_expressions),locations,gene_expressions,n_col=3)
```

![image](https://github.com/Zhangxf-ccnu/ENGEP-examples/blob/main/docs/osm_u_known.png)

For genes belonging to novel archetypes, we utilize hierarchical clustering to group them into distinct archetypes.

```
unsimilar_genes  = names(maxcor)[which(maxcor<threshold)]
corMat_un <- 
    stats::cor(as.matrix(t(pre_un@assays$RNA@scale.data[unsimilar_genes,])), 
               as.matrix(t(pre_un@assays$RNA@scale.data[unsimilar_genes,])),                
               method = "pearson")
dist_un <- 1-corMat_un
dist_un <- as.dist(dist_un)
hc <- hclust(dist_un, method = "ward.D")
rect.hclust(hcc,k=2)
id_un = cutree(hcc,k=2)
arche_exp_un <- 
    get_expression(id_un,pre_un@assays$RNA@scale.data[resultsfilter_un,])
arche_exp_un_log <- 
    get_expression(id_un,pre_un@assays$RNA@data[resultsfilter_un,])
gene_expressions <- t(as.data.frame(arche_exp_un_log))
plot_exp_pattern(ncol(as.data.frame(gene_expressions)),locations,gene_expres sions,n_col=2)
```
![image](https://github.com/Zhangxf-ccnu/ENGEP-examples/blob/main/docs/osm_novel.png)
