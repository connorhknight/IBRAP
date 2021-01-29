#' Here performs Seurat::SCTransfrom  normalisation
#'
#' SCTransfrom, Seurats novel normalisation technique
#'
#' @import Seurat
#' @param object Please specify a SCE object produced using IBRAP functions.
#' @param split.by Which metadata should the matrix be split by. 
#' @param which.assay Which assay should be used to perform SCTransform. 
#' @param new.assay What name should the new assay be saved under. 
#' @param ... Specify parameters to be passed to the SCTransform function.
#'
#' @examples normalisation_SCT(object = sce_object, split.by = NULL, which.assay = 'decontXcounts', new.assay = 'sctransform')
#' @export

perform.sct.normalisation <- function(object, 
                                      split.by=NULL, 
                                      which.assay, 
                                      new.assay = 'sctransform', 
                                      ...) {
  if(is.null(split.by)) {
    seuratobj <- CreateSeuratObject(counts = assay(object, which.assay), project = 'NA')
    seuratobj <- SCTransform(object = seuratobj, ...)
    genes <- intersect(rownames(object), rownames(as.matrix(seuratobj@assays$SCT@data)))
    object <- object[genes,]
    assay(object, new.assay) <- as.matrix(seuratobj@assays$SCT@data)[genes,]
    return(object)
  } else {
    t <-unique(colData(object)[[split.by]])
    list.matrix <- list()
    cat(crayon::cyan('loading datasets\n'))
    for(o in t) {
      cat(crayon::cyan(paste0('calculating: ', o, '\n')))
      sub <- object[,colData(object)[[split.by]] == o]
      seuratobj <- CreateSeuratObject(counts = assay(sub, which.assay), project = 'NA')
      seuratobj <- SCTransform(object = seuratobj, ...)
      list.matrix[[o]] <- seuratobj@assays$SCT@data
    }
    cat(crayon::cyan('Aligning features: stage 1\n'))
    genes.length <- list()
    for(o in t) {
      genes.length[[o]] <- length(rownames(list.matrix[[o]]))
    }
    
    cat(crayon::cyan('Aligning features: stage 2\n'))
    highest.rows <- names(which.max(rank(x = unlist(genes.length))))
    gene.list <- rownames(list.matrix[[highest.rows]])
    
    for(o in t[!t %in% highest.rows]){
      gene.list <- intersect(gene.list, rownames(list.matrix[[o]]))
    }
    
    cat(crayon::cyan('Aligning features: stage 3\n'))
    
    for(o in t) {
      list.matrix[[o]] <- list.matrix[[o]][gene.list,]
    }
    
    cat.mat <- do.call('cbind', list.matrix)
    cat(crayon::cyan('Aligning features: complete\n'))
    empty <- object
    common_rows <- intersect(rownames(cat.mat), rownames(empty))
    empty <- empty[common_rows,]
    cat.mat <- cat.mat[common_rows,]
    
    cat(crayon::cyan('Matrix added to S4 object\n'))
    
    assay(empty, new.assay) <- cat.mat
    cat(crayon::cyan('Complete\n'))
    return(empty)
  }
}

