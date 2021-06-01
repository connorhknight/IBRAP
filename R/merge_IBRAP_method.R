#' @title Method override for merge function regarding IBRAP S4 objects

setMethod(f = 'merge', signature = 'IBRAP',
          function(x, 
                   y){
            
            items <- c(x,y)
            
            for(i in items) {
              
              if(length(i@methods[[1]]) > 1) {
                
                cat(crayon::cyan('No analysis can be performed prior to merging\n'))
                return(NULL)
                
              }
              
            }
            
            column.names <- list()
            counts.list <- list()
            sample.list <- list()
            
            for(i in 1:length(items)) {
              
              column.names[[i]] <- colnames(items[[i]]@methods[[1]]@counts)
              counts.list[[i]] <- as.matrix(items[[i]]@methods[[1]]@counts)
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
            .feature_metadata <- as.data.frame(items[[1]]@methods[[1]]@feature_metadata)
            
            for(t in 2:length(items)) {
              
              .counts <- merge(x = .counts, y = counts.list[[t]], by = 'row.names', all = T)
              rownames(.counts) <- .counts$Row.names
              .counts$Row.names <- NULL
              
              .sample_metadata <- rbind(.sample_metadata, sample.list[[t]])
              
              .feature_metadata <- merge(.feature_metadata, 
                                         as.data.frame(items[[t]]@methods[[1]]@feature_metadata), 
                                         by = 'row.names', all = T)
              .feature_metadata[is.na(.feature_metadata)] <- 0
              rownames(.feature_metadata) <- .feature_metadata$Row.names
              .feature_metadata$Row.names <- NULL
              .feature_metadata$total.cells.x <- .feature_metadata[,1] + .feature_metadata[,3]
              .feature_metadata$total.counts.x <- .feature_metadata[,2] + .feature_metadata[,4]
              .feature_metadata <- .feature_metadata[,1:2]
              colnames(.feature_metadata) <- c('total.cells', 'total.counts')
              print('.')
            }
            
            .counts[is.na(.counts)] <- 0
            
            .counts <- Matrix::Matrix(data = as.matrix(.counts), sparse = T)
            .sample_metadata[match(colnames(.counts), rownames(.sample_metadata)),]
            .feature_metadata[match(rownames(.counts), rownames(.feature_metadata)),]
            
            new.method <- list()
            
            new.method[[names(x@methods)[1]]] <- new(Class = 'methods',
                                                     counts = .counts,
                                                     feature_metadata = .feature_metadata)
            
            ibrap <- new(Class = 'IBRAP',
                         methods = new.method, 
                         sample_metadata = .sample_metadata)
            
            return(ibrap)
            
          })