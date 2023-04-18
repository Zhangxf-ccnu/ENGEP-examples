# Identification of novel spatial archetypes in MERFISH

This example use ENGEP to predict unmeasured genes and detect novel spatial archetypes for MERFISH.

#### Get the predicted expression levels of unmeasured genes

Load ENGEP package.

```
library(ENGEP)
```

Read the spatial dataset. 

```
merfish <- readRDS(".../merfish.rds")
```

Read multiple reference datasets and put them into a list.

```
#read the first reference
cellsmart <- readRDS(".../cellsmart_qc.rds")

#read the second reference
nucsmart <- readRDS(".../nucsmart_qc.rds")

#read the third reference
nu10v2 <- readRDS(".../nu10v2_dgc.rds")

#read the forth reference
nu10v3a <- readRDS(".../nu10v3a_dgc.rds")

#read the fifth reference
sc10v2 <- readRDS(".../sc10v2_dgc.rds")

#read the sixth reference
nu10v3b <- readRDS(".../nu10v3b_dgc.rds")

#read the svevnth reference
sc10v3a <- readRDS(".../sc10v3a_dgc.rds")

# get the reference list
ref_list = list(cellsmart,nucsmart,nu10v2,nu10v3a,sc10v2,nu10v3b,sc10v3a)
```

Here, we predict 5234 unmeasured genes ( the union of 2000 highly variable genes in each reference dataset) for MERFISH.

```
common_hvg <- read.csv('.../union_hvg.csv',row.names=1)
pre_genes <- common_hvg$x
```

The predicted results can be obtained by running this function in ENGEP package.

```
engep_predict <- engep_predict(merfish,ref_list,pre_genes)
```

#### Identification of known spatial archetypes

First, we use MERINGUE to determine the spatially variable genes.

```
suppressMessages(library(MERINGUE))
cell_loc = read.csv('.../m2s300loc.csv',row.names = 1)
mat <- normalizeCounts(counts = merfish,
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
mer <- CreateSeuratObject(merfish,project = "sc_use",min.cells = 0,min.features = 0)
mer <- NormalizeData(object = mer,normalization.method = "LogNormalize")
mer <- ScaleData(object = mer)

mer_fil <- as.matrix(mer@assays$RNA@scale.data[resultsfilter,])

corMatt <- stats::cor(as.matrix(t(mer_fil)),
                      as.matrix(t(mer_fil)), method = "pearson")

dd <- 1-corMatt
dd<-as.dist(dd)

hc <- hclust(dd, method = "ward.D")
```

 The number of archetypes is determined through visual inspection of the cluster dendrogram.

```
plot(hc,hang = -1)
out.id = cutree(hc,k=5)
```

In this paper, we identify five known archetypes in MERFISH. The spatially variable genes are grouped into five clusters. The average expression of genes in each cluster is considered to represent the expression profile of a known archetype.

```
arche_exp <- get_expression(out.id,mer@assays$RNA@scale.data[resultsfilter,])
```

We can get expression patterns of these known archetypes.

```
arche_exp_log <- get_expression(out.id,mer@assays$RNA@data[resultsfilter,])
gene_expressions = t(as.data.frame(arche_exp_log))
locations = as.data.frame(cell_loc)
locations = rename(locations, X0=X, X1=Y)
plot_exp_pattern(ncol(gene_expressions),locations,gene_expressions,n_col=5)
```
![image](https://github.com/Zhangxf-ccnu/ENGEP-examples/blob/master/docs/known_archetype.png)



#### Identification of novel spatial archetypes

For each gene that has not been spatially measured, we aim to identify its closest known spatial archetype based on its predicted gene expression and the expression profifile of known archetypes.

We use MERINGUE to determine the spatially variable genes in the predicted genes.

```
pre_un = t(as.matrix(engep_predict))
counts <- cleanCounts(counts = pre_un,plot=FALSE,verbose=TRUE)
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

corMat <- stats::cor(as.matrix(t(arche_exp)),
                     as.matrix(t(as.matrix(pre_un@assays$RNA@scale.data) [resultsfilter_un,])), method = "pearson")
```

We use Pearson correlation as a metric to calculate the likelihood score of the unmeasured gene.  Then, we apply Hartigan’s dip test to analyze the distribution of likelihood scores of all unmeasured genes in order to detect any novel archetypes. 

```
maxcor <- apply(corMat, 2, max)
diptest::dip.test(maxcor, simulate.p.value = FALSE, B = 2000)
#Hartigans' dip test for unimodality / multimodality
#data:  maxcor
#D = 0.032315, p-value < 2.2e-16
#alternative hypothesis: non-unimodal, i.e., at least bimodal
```

We use manifold mixup to partition the unmeasured genes into two groups.

```
c <- apply(as.matrix(pre_un@assays$RNA@scale.data[resultsfilter_un,],
           1,
           partition_threshold)

diag(c) <- 0
threshold <- sum(c)/((dim(pre_un)[1])*(dim(pre_un)[1]-1))
plot_score(maxcor,threshold=0.232884)
```
![image](https://github.com/Zhangxf-ccnu/ENGEP-examples/blob/master/docs/scores.png)

For genes belonging to known archetypes, we assign them to the closest matching known archetype.

```
similar_genes = names(maxcor)[which(maxcor>=threshold)]
corMat_similar = corMat[,similar_genes ]

cellidx = apply(corMat_similar ,2,which.max)

arche_exp_n <- get_expression(cellidx,pre_un@assays$RNA@scale.data[resultsfilter_un,])
#plot
arche_exp_n_log <- get_expression(cellidx,pre_un@assays$RNA@data[resultsfilter_un,])
gene_expressions <- t(as.data.frame(arche_exp_n_log))
plot_exp_pattern(ncol(gene_expressions),locations,gene_expressions,n_col=5)
```
![image](https://github.com/Zhangxf-ccnu/ENGEP-examples/blob/master/docs/known_n.png)

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
rect.hclust(hcc,k=3)
id_un = cutree(hcc,k=3)
arche_exp_un <- 
    get_expression(id_un,pre_un@assays$RNA@scale.data[resultsfilter_un,])
#plot
arche_exp_un_log <- 
    get_expression(id_un,pre_un@assays$RNA@data[resultsfilter_un,])
gene_expressions <- t(as.data.frame(arche_exp_un_log))
plot_exp_pattern(ncol(as.data.frame(gene_expressions)),locations,gene_expressions,n_col=3)
```
![image](https://github.com/Zhangxf-ccnu/ENGEP-examples/blob/master/docs/novel.png)

