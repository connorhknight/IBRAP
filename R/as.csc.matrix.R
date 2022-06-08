#' @name as.csc.matrix
#' @aliases as.csc.matrix
#' 
#' @title converts to sparse matrix
#'
#' @return dgCMatrix

as.csc.matrix <- function (x, binary = FALSE, logical = FALSE, sort = FALSE) 
{
  if (binary && logical) 
    stop("Can pass only one of 'binary' or 'logical'.")
  if ((inherits(x, "dgCMatrix") && !binary && !logical) || 
      (inherits(x, "ngCMatrix") && binary) || (inherits(x, 
                                                        "lgCMatrix") && logical)) {
    return(x)
  }
  if (inherits(x, c("numeric", "integer", "logical", "data.frame"))) 
    x <- as.matrix(x)
  if (!inherits(x, "CsparseMatrix")) 
    x <- as(x, "CsparseMatrix")
  if (inherits(x, c("symmetricMatrix", "triangularMatrix"))) {
    if (!inherits(x, "dsparseMatrix")) 
      x <- as(x, "dsparseMatrix")
    x <- as(x, "dgCMatrix")
  }
  if (!binary && !logical && !inherits(x, "dgCMatrix")) {
    X_attr <- attributes(x)
    X_attr$class <- "dgCMatrix"
    if (.hasSlot(x, "x")) 
      X_attr$x <- as.numeric(X_attr$x)
    else X_attr$x <- rep(1, length(X_attr$i))
    if ("diag" %in% names(X_attr)) 
      X_attr$diag <- NULL
    if ("uplo" %in% names(X_attr)) 
      X_attr$uplo <- NULL
    attributes(x) <- X_attr
  }
  if (logical && !inherits(x, "lgCMatrix")) {
    X_attr <- attributes(x)
    X_attr$class <- "lgCMatrix"
    if (.hasSlot(x, "x")) 
      X_attr$x <- as.logical(X_attr$x)
    else X_attr$x <- rep(TRUE, length(X_attr$i))
    if ("diag" %in% names(X_attr)) 
      X_attr$diag <- NULL
    if ("uplo" %in% names(X_attr)) 
      X_attr$uplo <- NULL
    attributes(x) <- X_attr
  }
  if (binary && !inherits(x, "ngCMatrix")) {
    X_attr <- attributes(x)
    X_attr$class <- "ngCMatrix"
    if ("x" %in% names(X_attr)) 
      X_attr$x <- NULL
    if ("diag" %in% names(X_attr)) 
      X_attr$diag <- NULL
    if ("uplo" %in% names(X_attr)) 
      X_attr$uplo <- NULL
    attributes(x) <- X_attr
  }
  if (sort) 
    X <- sort_sparse_indices(X, copy = TRUE)
  return(x)
}