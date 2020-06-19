#' Title
#'
#' @param root
#' @param dir_regex
#' @param file_regex
#' @param name_key
#'
#' @return
#'
#' @examples
.get_files = function(root, dir_regex, file_regex, name_key){
    out_files = root %>%
        dir(., pattern = dir_regex, full.names = TRUE) %>%
        dir(., pattern = file_regex, full.names = TRUE)
    names(out_files) = sapply(strsplit(basename(out_files), "[-_\\.]"), function(x)paste(x[name_key], collapse = "_"))
    out_files
}

#' returns files from a root directory stored in sub directories per sample
#'
#' @param root directory containing directories to search
#' @param repKey parts of reps directories to use for naming. [-_\\.] delimitted.
#' @param poolKey parts of pooled directories to use for naming. [-_\\.] delimitted.
#'
#' @return
#' @export
#'
#' @examples
get_my_files = function(root, repKey = 2:4, poolKey = 2:3){
    my_files = list()
    my_files$root = root
    my_files$config = dir(root, pattern = ".csv$", full.names = TRUE)
    #bams
    my_files$bam_pooled = .get_files(my_files$root, "_pooled$", "bam$", poolKey)
    # my_files$bam_reps_input = .get_files(my_files$root, "input_R[0-9a-z]+$", "bam$", repKey)
    # my_files$bam_reps = setdiff(.get_files(my_files$root, "_R[0-9a-z]+$", "bam$", repKey),
    #                             my_files$bam_reps_input)
    my_files$bam_reps = .get_files(my_files$root, "_R[0-9a-z]+$", "bam$", repKey)
    #bigwigs
    my_files$bigwig_pooled = .get_files(my_files$root, "_pooled$", "FE.bw", poolKey)
    my_files$bigwig_reps = .get_files(my_files$root, "_R[0-9a-z]+$", "FE.bw", repKey)
    #peaks
    my_files$peaks_idr = .get_files(my_files$root, "_pooled$", "IDR.+Peak$", poolKey)
    my_files$peaks_pooled = .get_files(my_files$root, "_pooled$", "peaks.narrowPeak$", poolKey)

    my_files$peaksXls_reps = .get_files(my_files$root, "_R[0-9a-z]+$", "R[0-9]_peaks.xls$", repKey)
    my_files$peaksXls_pooled = .get_files(my_files$root, "_pooled$", "pooled_peaks.xls$", repKey)
    return(my_files)
}

#' parse lines in xls
#'
#' @param xls_file xls file
#' @param key key to match
#'
#' @return
#' @export
#'
#' @examples
parse_xls = function(xls_file, key){
    str = read.table(xls_file, nrows = 30, comment.char = "", sep = "\n", stringsAsFactors = FALSE)[,1]
    str[grepl(key, str)] %>% sub(".+: ", "", .)
}

#' parse lines in xls for treatment and control
#'
#' @param xls_file xls file
#'
#' @return
#' @export
#'
#' @examples
parse_xls_reads = function(xls_file){
    # str = read.table(xls_file, nrows = 30, comment.char = "", sep = "\n", stringsAsFactors = FALSE)[,1]
    # in_treat = str[grepl("total tags in treatment", str)] %>% sub(".+: ", "", .) %>% as.numeric()
    # in_ctrl = str[grepl("tags after filtering in control", str)] %>% sub(".+: ", "", .) %>% as.numeric()
    in_treat = parse_xls(xls_file, "total tags in treatment") %>% as.numeric # str[grepl("total tags in treatment", str)] %>% sub(".+: ", "", .) %>% as.numeric()
    in_ctrl = parse_xls(xls_file, "total tags in control") %>% as.numeric #str[grepl("tags after filtering in control", str)] %>% sub(".+: ", "", .) %>% as.numeric()
    c(treatmet = in_treat, control = in_ctrl)
}
