unitize_R2 <- function(w,alpha) {
  N <- ncol(w)  
  rs <- rowSums(w)
  rs[rs == 0] <- 1
  w1 <- apply(w, 2, function(row) row / rs)
  w2 <- diag(N)
  weight = alpha*w1 + (1-alpha)*w2
  return(weight)
}

CellTypeSpatialCrossCor_quick <- function(gexpA, gexpB, weight) {
  n1 =nrow(gexpA)
  n2 =nrow(gexpB)
  colmeanA = colMeans(gexpA)
  colmeanB = colMeans(gexpB)
  
  N = n1+n2
  W = 2*sum(weight)
  
  x = gexpA- matrix(rep(colmeanA,n1),n1, byrow =TRUE)
  y = gexpB- matrix(rep(colmeanB,n2),n2, byrow =TRUE)
  
  cv1 = t(x) %*% (weight+t(weight)) %*% y
  cv = cv1/2
  
  v = sqrt(colSums(x^2)) %o% sqrt(colSums(y^2))
  
  iSCI <- (N/W) * (cv/v) 
  
  return(iSCI)
}


get_expression <- function(out.id,mea){
  metaexp = matrix(0,nrow = length(table(out.id)),ncol = dim(mea)[2])
  for (i in 1:length(table(out.id))){
    if (table(out.id)[[i]] == 1){
      metaexp[i,] =  mea[names(out.id)[which(out.id==i)],]
    }else{
      metaexp[i,] <- apply(mea[names(out.id)[which(out.id==i)],],2,mean)
    }
  }
  return (metaexp)
}

get_expression_novel <- function(genelist_use,mea){
  metaexp = matrix(0,nrow = length(genelist_use),ncol = dim(mea)[2])
  for (i in 1:length(genelist_use)){
    
    metaexp[i,] <- apply(mea[genelist_use[[i]],],2,mean)
  }
  return (metaexp)
}

get_expression_d <- function(out.id,mea,genesD0){
  ll = length(table(out.id))+length(genesD0)
  metaexp = matrix(0,nrow = ll,ncol = dim(mea)[2])
  for (i in 1:length(table(out.id))){
    if (table(out.id)[[i]] == 1){
      metaexp[i,] =  mea[names(out.id)[which(out.id==i)],]
    }else{
      metaexp[i,] <- apply(mea[names(out.id)[which(out.id==i)],],2,mean)
    }
  }
  for (j in 1:length(genesD0)){
    metaexp[(j+length(table(out.id))),] = mea[genesD0[j],]
    }
  return (metaexp)
}

remove_degree_one_nodes <- function(adjacency_matrix) {
  graph <- graph_from_adjacency_matrix(adjacency_matrix, mode = "undirected")
  while (any(igraph::degree(graph) <=10)) {
    degree_one_nodes <- which(igraph::degree(graph) <=10)
    graph <- delete_vertices(graph, degree_one_nodes)
  }
  
  if (vcount(graph) == 0)
    final_nodes = NULL
  else{
    final_adjacency_matrix <- as.matrix(as_adjacency_matrix(graph))
    final_nodes = colnames(final_adjacency_matrix)}
  final_nodes
}

CellTypeSelection = function(y, z, tau=1.05){
  Coryz = cor(y,z)
  idx = sort(Coryz, decreasing = TRUE, index.return=TRUE)$ix
  CC = rep(0,length(idx))
  CC[1] = cor(y, z[,idx[1]])
  label = z[,idx[1]]
  change = rep(length(idx)-1)
  for (i in 2:length(idx)){
    label = label + z[,idx[i]]
    CC[i] = cor(y, label)
    change[i-1] = CC[i]/CC[i-1]
    if (change[i-1]<tau){
      celltype = idx[1:(i-1)]
      break
    }
    
  }
  celltype
}
