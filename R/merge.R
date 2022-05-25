#' @name merge
#' 
#' @title Method override for merge function regarding IBRAP S4 objects
#' 
#' @exportMethod merge

setMethod(f = 'merge', signature = 'IBRAP',
          function(x, 
                   y){
            
            items <- c(x,y)
            
            feature.list <- list()
            
            for(i in 1:length(items)) {
              
              feature.list[[i]] <- rownames(items[[i]]@methods[[1]]@counts)
              
            }
            
            mutual_features <- Reduce(f = intersect, feature.list)
            
            column.names <- list()
            counts.list <- list()
            sample.list <- list()
            
            for(i in 1:length(items)) {
              
              column.names[[i]] <- colnames(items[[i]]@methods[[1]]@counts)
              counts.list[[i]] <- as_matrix(items[[i]]@methods[[1]]@counts)[mutual_features,]
              sample.list[[i]] <- as.data.frame(items[[i]]@sample_metadata)
              
            }
            
            column.names <- unlist(column.names)
            column.names[!isUnique(column.names)] <- make.unique(column.names[!isUnique(column.names)])
            
            for(q in 1:length(items)) {
              
              colnames(counts.list[[q]]) <- column.names[1:length(colnames(counts.list[[q]]))]
              rownames(sample.list[[q]]) <- column.names[1:length(rownames(sample.list[[q]]))]
              column.names <- column.names[sum(length(colnames(counts.list[[q]]))+1):length(column.names)]
              
            }
            
            .counts <- counts.list[[1]]
            .sample_metadata <- sample.list[[1]]
            .feature_metadata <- as.data.frame(items[[1]]@methods[[1]]@feature_metadata[mutual_features,])
            
            pb <- progress::progress_bar$new(total = sum(length(y)*4)+1)
            
            for(t in 2:length(items)) {
              
              if(!isUnique(c(colnames(.counts), colnames(counts.list[[t]])))) {
                
                names <- make.unique(c(colnames(.counts), colnames(counts.list[[t]])))
                
                colnames(.counts) <- names[1:ncol(.counts)]
                colnames(counts.list[[t]]) <- names(sum(ncol(.counts)+1):ncol(counts.list[[t]]))
                
              }
              
              pb$tick()

              .counts <- merge(x = .counts, y = counts.list[[t]], by = 'row.names', all = T)
              rownames(.counts) <- .counts$Row.names
              .counts$Row.names <- NULL
              
              pb$tick()
  
              .sample_metadata <- merge(.sample_metadata, sample.list[[t]], all=T)
              
              pb$tick()
              
              .feature_metadata <- merge(.feature_metadata, 
                                         as.data.frame(items[[t]]@methods[[1]]@feature_metadata[mutual_features,]), 
                                         by = 'row.names', all = T)
              
              .feature_metadata[is.na(.feature_metadata)] <- 0
              rownames(.feature_metadata) <- .feature_metadata$Row.names
              .feature_metadata$Row.names <- NULL
              .feature_metadata$total.cells.x <- .feature_metadata[,1] + .feature_metadata[,3]
              .feature_metadata$total.counts.x <- .feature_metadata[,2] + .feature_metadata[,4]
              .feature_metadata <- .feature_metadata[,1:2]
              colnames(.feature_metadata) <- c('total.cells', 'total.counts')
              
              pb$tick()

            }
            
            feat.met <- feature_metadata(assay = .counts, col.prefix = 'RAW')
            samp.met <- cell_metadata(assay = .counts, col.prefix = 'RAW')
            
            .counts[is.na(.counts)] <- 0
            
            .counts <- Matrix::Matrix(data = as_matrix(.counts), sparse = T)
            
            .sample_metadata[match(colnames(.counts), rownames(.sample_metadata)),]
            .sample_metadata[,which(grepl(pattern = 'total.counts', x = colnames(.sample_metadata)))] <- samp.met[,1]
            .sample_metadata[,which(grepl(pattern = 'total.features', x = colnames(.sample_metadata)))] <- samp.met[,2]
            
            .feature_metadata[match(rownames(.counts), rownames(.feature_metadata)),]
            .feature_metadata[,which(grepl(pattern = 'total.counts', x = colnames(.sample_metadata)))] <- feat.met[,1]
            .feature_metadata[,which(grepl(pattern = 'total.cells', x = colnames(.sample_metadata)))] <- feat.met[,2]
            
            pb$tick()
            
            new.method <- list()
            
            new.method[[names(x@methods)[1]]] <- new(Class = 'methods',
                                                     counts = .counts,
                                                     normalised = as(matrix(nrow = 0, ncol = 0), 'dgCMatrix'),
                                                     feature_metadata = .feature_metadata)
            
            ibrap <- new(Class = 'IBRAP',
                         methods = new.method, 
                         sample_metadata = .sample_metadata)
            
            pb$finished
            
            return(ibrap)
            
          })



