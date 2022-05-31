#' @name as_matrix_transpose
#' @aliases as_matrix_transpose
#' 
#' @title converts to transposed matrix
#'
#' @description converts sparse matrix to normal transposed matrix for sparse data
#' 
#' @param mat
#' 
#' @return matrix

as_matrix_transpose <- function(mat){
  print('.')
  tmp <- matrix(data=0L, nrow = mat@Dim[2], ncol = mat@Dim[1])
  print('.')
  row_pos <- mat@i+1
  print('.')
  col_pos <- findInterval(seq(mat@x)-1,mat@p[-1])+1
  print('.')
  val <- mat@x
  print('.')
  for (i in seq_along(val)){
    tmp[col_pos[i],row_pos[i]] <- val[i]
  }
  print('.')
  row.names(tmp) <- mat@Dimnames[[2]]
  print('.')
  colnames(tmp) <- mat@Dimnames[[1]]
  return(tmp)
}