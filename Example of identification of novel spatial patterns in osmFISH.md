# Identification of novel spatial patterns in osmFISH

This example use ENGEP to predict unmeasured genes for osmFISH and detect novel spatial patterns for osmFISH. The datasets used in this examples can be downloaded from Zenodo (https://doi.org/10.5281/zenodo.7825993)

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

Here, we predict 3433 unmeasured genes ( the union of 2000 highly variable genes in each reference dataset) for osmFISH. The code to find these 3433 unmeasured genes is also in this repository.

```
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
![image](https://github.com/Zhangxf-ccnu/ENGEP-examples/blob/master/docs/genes_of_osmFISH.png)

#### Identification of known spatial patterns

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
resultsfilter <- filterSpatialPatterns(mat = mat[test,],
                                        I = I,
                                        w = w,
                                        adjustPv = TRUE,
                                        alpha = 0.05,
                                        minPercentCells = 0.05,
                                        verbose = TRUE)
```

Next, we calculate a Spatial Cross-Correlation Index (SCI) for all remaining gene pairs (our approach differs from the original MERINGUE approach for that we consider self-loops when constructing the cell adjacency relationship network).  We retain elements in the similarity matrix that surpass a certain threshold (at the 10th percentile) and replace the rest with zeros.  We employ Louvain clustering on the constructed gene network to get known spatial patterns.

```
osm<- CreateSeuratObject(osmfish,project = "sc_use",min.cells = 0, 
                         min.features = 0)
osm <- NormalizeData(object = osm,normalization.method = "LogNormalize")
osm<- ScaleData(object = osm)

gexpA = t(as.matrix(osm@assays$RNA@scale.data[resultsfilter,]))
gexpB = t(as.matrix(osm@assays$RNA@scale.data[resultsfilter,]))
weight = unitize_R2(w,0.5)

scc2 = CellTypeSpatialCrossCor_quick(gexpA,gexpB,weight)
data = scc2
diag(data) = 0
tau = quantile(data[upper.tri(data, diag = FALSE)], probs = 0.9)
connection_matrix = data
connection_matrix[data<=tau] =0
connection_matrix = connection_matrix - diag(diag(connection_matrix))
D = apply(connection_matrix, 1, sum)
idx  = (D!=0)
connection_matrix1 = connection_matrix[idx,idx]

library(igraph)
graph <- graph.adjacency(connection_matrix1, mode = "undirected", weighted = TRUE)
set.seed(1234)
clu = cluster_leiden(graph,objective_function = "modularity",resolution_parameter = 1)
table(clu$membership)
# 1  2  3 
# 2  3  3
```

 In this paper, we identify six known patterns in osmFISH. Notably, three of these patterns consist of only one gene (considering that the number of genes in osmFISH is relatively small, we retain the identified patterns composed of just a single gene). The average expression of genes in each cluster is considered to represent the expression profile of a known pattern.

```
groups = clu$membership
names(groups) = rownames(connection_matrix1)
genesD0 = setdiff(resultsfilter,names(groups))
```

In this paper, we identify six known archetypes in osmFISH. We can plot the expression patterns of every known archetype.

```
arche_exp <-get_expression_d(groups,osm@assays$RNA@scale.data[resultsfilter,],genesD0)
gene_expressions = t(as.data.frame(arche_exp))
locations = as.data.frame(cell_loc)
plot_exp_pattern(ncol(gene_expressions),locations,gene_expressions,n_col=6)
```

![image](https://github.com/Zhangxf-ccnu/ENGEP-examples/blob/master/docs/osmfish_measured_pattern.png)

#### Identification of novel spatial patterns

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
resultsfilter_un <- filterSpatialPatterns(mat = mat[test,],
                                        I = I,
                                        w = w,
                                        adjustPv = TRUE,
                                        alpha = 0.05,
                                        minPercentCells = 0.05,
                                        verbose = TRUE)
```

We then calculate SCI scores between the predicted unmeasured genes and the known patterns.

```
pre_un <- CreateSeuratObject(pre_un,project = "sc_use",min.cells = 0,min.features = 0)
pre_un <- NormalizeData(object = pre_un,normalization.method = "LogNormalize")
pre_un <- ScaleData(object = pre_un)
                    
corMat <- CellTypeSpatialCrossCor_quick(t(arche_exp),
                    t(as.matrix(pre_un@assays$RNA@scale.data[resultsfilter_un,])),
                    weight)
```

We use SCI score as a metric to calculate the likelihood score of the unmeasured gene.  Then, we apply Hartigan’s dip test to analyze the distribution of likelihood scores of all unmeasured genes in order to detect any novel patterns. 

```
maxcor <- apply(corMat, 2, max)
corMat = t(corMat)
diptest::dip.test(maxcor, simulate.p.value = FALSE, B = 2000)
#Hartigans' dip test for unimodality / multimodality
#data:  maxcor
#D = 0.014091, p-value = 0.0185
#alternative hypothesis: non-unimodal, i.e., at least bimodal
```

We employ a kernel density estimation method to estimate the density function of the probability distribution for likelihood scores, and  classify unmeasured genes into three distinct groups. We identify the two high peaks within this estimated density function and and determine the likelihood score at the low peak between the two high peaks, denoted as *s*. Genes with likelihood scores greater than $s+0.05$ are classified as members of known patterns. Conversely, genes with likelihood scores lower than $s-0.05$ are designated as belonging to novel patterns. The remaining genes are considered weakly associated with known patterns.

```
density_estimation <- density(maxcor)
max_indices <- which(diff(sign(diff(density_estimation$y))) == -2) + 1

min_density_index <- which.min(density_estimation$y[max_indices[2]:max_indices[3]])
x_value_at_min_density <- density_estimation$x[max_indices[2] + min_density_index]

plot_score(maxcor,density_estimation,x_value_at_min_density)
```

![image](https://github.com/Zhangxf-ccnu/ENGEP-examples/blob/master/docs/osmfish_density.png)

For genes belonging to known patterns, we assign them to the closest matching known pattern.

```
mmax = x_value_at_min_density+0.05
similar_genes = names(maxcor)[which(maxcor>=mmax)]
corMat_similar = corMat[similar_genes,]

cellidx = apply(corMat_similar,1,which.max)
table(cellidx)
#cellidx
#  1   2   3 
#161 157 336

arche_exp_n <- get_expression(cellidx,pre_un@assays$RNA@scale.data[resultsfilter_un,])
gene_expressions <- t(as.data.frame(arche_exp_n))
plot_exp_pattern(ncol(gene_expressions),locations,gene_expressions,n_col=3)
```

![image](https://github.com/Zhangxf-ccnu/ENGEP-examples/blob/master/docs/osm_unmeasured_known_pattern.png)

For genes belonging to novel patterns, we follow a meticulous three-step process to unveil these novel spatial patterns. First, we apply a filter to exclude genes with lower biological variability.  Specifically, we eliminate the lowest 10% of genes based on their mean expression levels. Subsequently, we further refine this set by removing the lowest 10% of genes with the least variability, employing the Seurat R package (FindVariableFeatures function).

```
mmin = x_value_at_min_density-0.05
unsimi_gene = names(maxcor[which(maxcor<mmin)])
gexpA = t(as.matrix(pre_un@assays$RNA@scale.data[unsimi_gene,]))
gexpB = t(as.matrix(pre_un@assays$RNA@scale.data[unsimi_gene,]))
scc_un <- CellTypeSpatialCrossCor_quick(gexpA,gexpB,weight)

gene_left1 = apply(ensem[unsimi_gene,],1,mean)
a = sort(gene_left1,decreasing = TRUE)[1:round(0.9*length(gene_left1))]
gene_left2 = names(a)

ensemn <- CreateSeuratObject(ensem[gene_left2,])
ensemn <- FindVariableFeatures(ensemn)
geneuse <- VariableFeatures(ensemn)[1:round(0.9*length(gene_left2))]
```
Then, we utilize a graph-based clustering method to categorize the remaining genes into different clusters.

```
datau = scc_un[geneuse,geneuse]
diag(datau) = 0
tau = quantile(datau[upper.tri(datau, diag = FALSE)], probs = 0.90)
connection_matrix = datau
connection_matrix[datau<=tau] =0
connection_matrix = connection_matrix - diag(diag(connection_matrix))
D = apply(connection_matrix, 1, sum)
idx  = (D!=0)
connection_matrix1 = connection_matrix[idx,idx]
graph <- graph.adjacency(connection_matrix1, mode = "undirected", weighted = TRUE)
set.seed(1234)
cluu = cluster_leiden(graph,objective_function = "modularity",resolution_parameter = 0.6)
table(cluu$membership)

#  1   2   3   4 
#249 149 319   4
```

Second, we filter out clusters based on a three-times standard deviation rule. 

```
groupu = cluu$membership
names(groupu) = rownames(connection_matrix1)
datauu = datau[rownames(connection_matrix1),rownames(connection_matrix1)]
mean_cluster =rep(0,length(unique(groupu)))
three_sigma =rep(0,length(unique(groupu)))
for (i in 1:length(unique(groupu))){
  scc1 = datauu[names(groupu)[which(groupu==i)],
                names(groupu)[which(groupu==i)]]
  mean_cluster[i]=mean(scc1)
  three_sigma[i]=mean(datauu)+ 3*sd(datauu)/sum(groupu==i)#
}
mean_cluster
#[1] 0.3478508 0.1303187 0.2995420 0.3140484
three_sigma
#[1] 0.1416579 0.1435540 0.1410379 0.3147070

iddk = which(mean_cluster>three_sigma)
iddk
#[1] 1 3
```

Third,  we refine the retained clusters by iteratively eliminating genes with fewer than ten connections within their respective sub-gene networks. These refined clusters, characterized by their average gene expression, are considered as novel patterns. 

```
genelist_use = list()
a=0
for (i in 1:length(iddk)){
  a=a+1
  genes_f1 <- names(groupu)[which(groupu==iddk[i])]
  connection_matrix1_u = connection_matrix1[genes_f1,genes_f1]
  connection_matrix1_u[connection_matrix1_u!=0]=1
  genes_f1_n = remove_degree_one_nodes(connection_matrix1_u)
  genelist_use[[i]] = genes_f1_n
}

arche_exp_un <-get_expression_novel(genelist_use,pre_un@assays$RNA@scale.data[resultsfilter_un,])
gene_expressions <- t(as.data.frame(arche_exp_un))
plot_exp_pattern(ncol(as.data.frame(gene_expressions)),locations,gene_expressions,n_col=2)

```

![image](https://github.com/Zhangxf-ccnu/ENGEP-examples/blob/master/docs/osm_novel_pattern.png)

Following the identification of novel spatial patterns, we conduct a multifaceted analysis to explore their biological significance. First, we select representative genes from each pattern.

```
for (i in 1:length(genelist_use)){
  corMatu <- stats::cor(t(arche_exp_un)[,i],
                        t(as.matrix(pre_un@assays$RNA@scale.data[genelist_use[[i]],])),
                        method = "pearson")

  corMattt = as.vector(corMatu)
  names(corMattt) = colnames(corMatu)
  a = order(corMattt, decreasing = TRUE)[1:5]
  print(names(corMattt)[a])
}

#[1] "Plod1"   "Sorbs3"  "Hbp1"    "Pnpla2"  "Plekhg1"
#[1] "Gprasp2" "Dleu7"   "Dlx5"    "Gm14204" "Rpp25" 
```

Next, we delve deeper into the predicted expression of these novel spatial patterns, investigating their potential associations with specific cell types within the tissue.

```
label = read.csv('.../osm_label.csv',row.names = 1)
label = label[colnames(pre_un@assays$RNA@scale.data),]
labelus = label
z =  model.matrix(~labelus-1)

#take the novel pattern 1 for an example

y = arche_exp_un[1,]
celltype = CellTypeSelection(y,z)
Coryz = cor(y,z)
Coryz[,celltype]
#labelusAstrocyte Gfap labelusOligodendrocyte Mature 
#0.3688995                     0.2990140
```

Finally, we employ Gene Ontology enrichment analysis using the clusterProfiler R package to unveil the biological processes associated with these novel spatial patterns. The background gene list comprises the union of genes included in the reference datasets.

```
background_gene_merfish <- readRDS(".../background_gene_osmfish.rds")
transID_back = bitr(background_gene_merfish,
                    fromType="SYMBOL",
                    toType="ENTREZID",
                    OrgDb="org.Mm.eg.db")
BP_f1_un = list()
for (i in length(genelist_use)){
  genes_f1_n = genelist_use[[i]]
  transID_f1 = bitr(genes_f1_n,
                    fromType="SYMBOL",
                    toType="ENTREZID",
                    OrgDb="org.Mm.eg.db")

  BP_f1 <- enrichGO(transID_f1$ENTREZID, "org.Mm.eg.db",
                    universe= transID_back$ENTREZID,
                    keyType="ENTREZID", ont="BP", pvalueCutoff=0.05)
  BP_f1 <- setReadable(BP_f1, OrgDb=org.Mm.eg.db)
  BP_f1_un[[i]] = BP_f1
}
```

