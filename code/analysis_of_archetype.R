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


partition_threshold <- function(pre_un2){
  m=runif(dim(pre_un2)[1] , min = 0, max = 1)
  corMat <- stats::cor(as.matrix(metaexp),#metaexp
                       t(as.matrix(m) %*% as.matrix(t(pre_un2)) +
                           matrix(rep((1-m),dim(pre_un2)[2]),nrow=length(m)) * as.matrix(pre_un2)),
                       method = "pearson")
  maxcor <- apply(corMat, 2, max)
  return(maxcor)
}

