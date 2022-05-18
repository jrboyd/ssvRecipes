#' truncate_peak_file
#'
#' @param np_file path to .narrowPeak file or any .bed file using the same format
#' @param order_VAR variable to order by, typically signalValue or pValue
#' @param out_dir output location of truncated peaks, default is input file directory with ".diffbind_input" appended.
#' @param top_fraction fraction of top ranking peaks to retain. Default is .1
#' @param top_n number of top ranking peaks to retain. supersedes top_fraction if set. Default is NULL.
#'
#' @return path to truncated file
#' @export
#' 
#' @examples
#' peaks2 = dir("/slipstream/galaxy/uploads/working/qc_framework/output/old_peaks/", 
#'              pattern = "H3K4ME3_pooled_peaks.narrowPeak$", full.names = TRUE)
#' peaks2.truncated_20percent = sapply(peaks2, truncate_peak_file, order_VAR= "signalValue", top_fraction = .2)
#' peaks2.truncated_5k = sapply(peaks2, truncate_peak_file, order_VAR= "pValue", top_n = 5e3)
truncate_peak_file = function(np_file, order_VAR,  out_dir = paste0(dirname(np_file), ".diffbind_input"), top_fraction = .1, top_n = NULL){
    df = read.table(np_file, col.names = c("seqnames", "start", "end", "id", "score", "strand", "signalValue", "pValue", "qValue", "summit"))
    stopifnot(order_VAR %in% colnames(df))
    if(is.null(top_n)){
        top_n = ceiling(nrow(df)*top_fraction)
    }
    df = df[order(df[[order_VAR]], decreasing = TRUE),]
    if(top_n < nrow(df)){
        df = df[seq(top_n),]
    }
    dir.create(out_dir, showWarnings = FALSE)
    out_file = file.path(out_dir, basename(np_file))
    write.table(df, out_file, col.names = FALSE, row.names = FALSE, quote = FALSE, sep = "\t")
    out_file
}



