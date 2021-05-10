#' @title An S4 class object that contains method-assays
#'
#' @description An S4 class object that contains assays and metadata.
#' 
#' @slot methods A list of methods class objects containing produced results
#' @slot sample_metadata A data.frame containing cell level metadata
#'
#' @export

setClass(Class = 'IBRAP', 
         representation = representation(
           methods = 'list', 
           sample_metadata = 'data.frame'
         ))