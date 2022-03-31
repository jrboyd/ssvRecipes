# BiocManager::install("IdeoViz")
# library(IdeoViz)

#' bfcif
#'
#' @param bfc
#' @param rname
#' @param FUN
#' @param force_overwrite
#' @param return_path_only boolean, if TRUE, FUN is never run and cache_path
#' is returned.  if file.exist(cache_path) == TRUE, data has already been
#' cached.
#' @param verbose boolean, if TRUE status messages are displayed
#'
#' @return
#' @export
#' @import BiocFileCache
#'
#' @examples
bfcif = function(bfc, rname, FUN, force_overwrite = FALSE, return_path_only = FALSE, verbose = TRUE){
    # is rname in cache?
    if(nrow(BiocFileCache::bfcquery(bfc, query = rname, field = "rname", exact = TRUE)) == 0){
        if(verbose) message("results not in cache. ", appendLF = FALSE)
        cache_path = BiocFileCache::bfcnew(bfc, rname = rname)

    }else{
        if(verbose) message("previous cache results found. ", appendLF = FALSE)
        cache_path = BiocFileCache::bfcrpath(bfc, rname)
    }
    if(return_path_only){
        if(verbose) message("returning cache path.")
        return(cache_path)
    }
    # does cached file exist?
    if(file.exists(cache_path) && !force_overwrite){
        if(verbose) message("loading previous cache results...")
        load(BiocFileCache::bfcrpath(bfc, rname))
    }else{
        if(verbose) message("running function...", appendLF = FALSE)
        res = FUN()
        if(verbose) message("caching results...")
        save(res, file = cache_path)
    }
    # return either new results or cached results
    res
}
