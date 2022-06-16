#' @name df_to_dgCMatrix
#' @aliases df_to_dgCMatrix
#' 
#' @title converts to dgCMatrix
#'
#' @description converts dataframe to  sparse matrix
#' 
#' @param mat
#' 
#' @return mdgCMatrix

df_to_dgCMatrix <- function(df, nsplit=1000) {
  
  splitMxList = lapply(split(df, cut(1:nrow(df), nsplit)), function(mx) {
    Matrix::Matrix(as.matrix(mx), sparse=T)
  })
  sparse.M = Reduce(rbind, splitMxList)
  return(sparse.M)
  
}