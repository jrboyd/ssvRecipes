#' get_hclust_gls
#'
#' All of this is required to get cluster gene lists in the "as plotted" order with cluster 1 at the top.  
#' Cluster numbers when cutting trees are normally arbitrary. 
#' Also gene order when cutting tree would not normally be "as plotted"
#'
#' @param hclust_res result from hclust().
#' @param nclust number of clusters to cutree with.
#'
#' @return
#' @export
#'
#' @examples
get_hclust_gls = function(hclust_res, nclust){
    
    tr = hclust_res
    
    gl_as_plotted = tr$labels[tr$order]
    
    clust_idents = cutree(tr, nclust)
    
    clust_idents = clust_idents[gl_as_plotted]
    
    clust_arbitrary2plotted = seq_len(nclust)
    names(clust_arbitrary2plotted) = unique(clust_idents)
    
    clust_idents.plotted = clust_arbitrary2plotted[as.character(clust_idents)]
    names(clust_idents.plotted) = names(clust_idents)
    
    clust_idents = clust_idents.plotted
    
    table(clust_idents)
    
    clust_gls = split(names(clust_idents), clust_idents)
    
    stopifnot(all(gl_as_plotted == unlist(clust_gls)))
    clust_gls
}

#' write_geneList_matrix
#'
#' @param gl list of gene ids, must have cluster ids as naames.
#' @param file csv file to write to.
#'
#' @return
#' @export
#'
#' @examples
write_geneList_matrix = function(gl, file){
    if(is.null(names(gl))) stop("gene lists must be named")
    gmat = matrix("", nrow = max(lengths(gl)), ncol = length(gl))
    for(i in seq_along(gl)){
        x = gl[[i]]
        gmat[seq_along(x), i] = as.character(x)
    }
    colnames(gmat) = names(gl)
    write.csv(gmat, file, row.names = FALSE, quote = FALSE)
}
