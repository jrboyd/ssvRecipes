#' is_url
#'
#' check if character f is a valid URL
#'
#' @param f a character possibly containing a URL
#'
#' @return TRUE if f looks like a URL
#' @export
#'
#' @examples
#' is_url("not_a_url")
#' is_url("ftp://is_a_url")
is_url = function(f){
    grepl("://", f)
}

#' Title
#'
#' @param f
#'
#' @return
#' @export
#'
#' @examples
is_gz = function(f){
    grepl(".gz$", f)
}

#' cache_gz
#'
#' Automatically handles cacheing of (potentially) gzipped resources.  Resources
#' can be local or remote (will be downloaded and cached).  If resource is
#' gzipped (only if ends in .gz) a gunzipped version will be added to cache.
#'
#' for a given path
#' 1. checks if it's a url and downloads, and updates to local path.
#' 2. if path ends in .gz, gunzips and updates path
#'
#' @param x a BiocFileCache cache
#' @param path url or path to gzipped resource
#' @param rname rname for gzipped resource
#'
#' @return path to gunzipped resource in cache
#' @export
#' @import BiocFileCache
#'
#' @examples
cache_gz = function(x, path){
    if(is_url(path)){
        if(is_gz(path)){
            path = BiocFileCache::bfcrpath(x,
                                           rnames = paste0(basename(path), ",gzip"),
                                           fpath = path)
            rname = paste0(basename(path), ",unzip")
            if(!any(BiocFileCache::bfcinfo(x)$rname == rname)){ #dl if necessary
                raw_lines = readLines(gzfile(path))
                path = BiocFileCache::bfcnew(x, rname = rname)
                writeLines(raw_lines, con = path)
            }else{ # use existing resource
                path = BiocFileCache::bfcrpath(x, rnames = rname)
            }
        }else{
            path = BiocFileCache::bfcrpath(x,
                                           rnames = paste0(basename(path), ",unzip"),
                                           fpath = path)
        }

    }
    if(grepl(".gz$", path)){ #gunzip if necessary


    }
    return(path)
}

#' Title
#'
#' @param x a bioc cache
#' @param input_name path to file FUN works on
#' @param FUN function that accepts one file path as it's only argument.
#' the output of FUN(input_name) will be returned, either from running or cached
#' if available.
#' @param output_tag character.  Used to create cache result identifier (rname).
#' New rname
#' will be either 1) the result of adding this tag to rname of input_name
#' if input_name is cached or 2) the result of adding this tag to basename of
#' input_name if uncached.  If this new rname exists in cache, cached result
#' will be used and FUN will not be rerun.
#'
#' @return object returned by FUN(input_obj)
#' @export
#' @import BiocFileCache
#'
#' @examples
cache_FUN = function(x, input_name, output_tag, FUN, input_obj = input_name, saveFUN = saveRDS, loadFUN = readRDS){
    if(any(BiocFileCache::bfcinfo(x)$rpath == basename(input_name))){
        input_rname = BiocFileCache::bfcinfo(x)$rname[BiocFileCache::bfcinfo(x)$rpath == basename(input_name)]
        output_rname = paste0(sub(',.+', "", input_rname), ",", output_tag)
    }else{
        output_rname = paste0(basename(input_name), ",", output_tag)
    }
    if(output_rname %in% BiocFileCache::bfcinfo(x)$rname){
        message("using cached result from: ", output_rname)
        out = loadFUN(BiocFileCache::bfcrpath(x, output_rname))
    }else{
        message("calculating result and storing in: ", output_rname)
        out = FUN(input_obj)
        saveFUN(out, BiocFileCache::bfcnew(x, rname = output_rname))
    }
    return(out)
}

#' Title
#'
#' @return
#' @export
#'
#' @examples
running_jids = function(){
    qstat_res = system("qstat", intern = TRUE)
    jids = sapply(strsplit(qstat_res[-(1:2)], " +"), function(x)x[2])
    return(jids)
}

#' Title
#'
#' @param hold_jids
#'
#' @return
#' @export
#'
#' @examples
wait_jids = function(hold_jids, interval = 5){
    hold_jids = hold_jids[!is.na(hold_jids)]
    jids = running_jids()
    while(any(hold_jids %in% jids)){
        message("Waiting on ", paste(intersect(hold_jids, jids), collapse = ","))
        jids = running_jids()
        Sys.sleep(interval)
    }
    message("Finished!")
}

#' Title
#'
#' @param hold_jids
#' @param interval
#' @param watch_prefix
#'
#' @return
#' @export
#'
#' @examples
watch_jids = function(hold_jids, interval = 5, watch_prefix = "watching"){
    hold_jids = hold_jids[!is.na(hold_jids)]
    jids = running_jids()
    if(any(hold_jids %in% jids)){
        message(watch_prefix, " : ", paste(intersect(hold_jids, jids), collapse = ","))
        later::later(function()watch_jids(hold_jids, interval, watch_prefix), interval)
    }else{
        message(watch_prefix, " : ", "finished!")
    }
}

#' Write gr to bigbed file
#'
#' @param gr gr to save
#' @param bedf name of intermediate bedfiles
#' @param bbf name of final output bigbed file
#' @param chrSizes chrSizes file to use, default is hg38
#' @param cleanup_beds remove beds when done? TRUE
#'
#' @return path to big bed file
#' @export
#'
#' @examples
writeBigBed = function(gr, bbf, bedf = sub("\\.bb$", "\\.bed", bbf), chrSizes = "~/hg38_chrsizes.txt", cleanup_beds = TRUE){
    rtracklayer::export.bed(gr, bedf)
    system(paste("bedSort", bedf, paste0(bedf, ".sorted")), intern = TRUE)
    system(paste("bedToBigBed", paste0(bedf, ".sorted"), chrSizes, bbf), intern = TRUE)
    if(cleanup_beds){
        file.remove(bedf)
        file.remove(paste0(bedf, ".sorted"))
    }
    bbf
}
