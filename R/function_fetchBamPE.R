
#' ssvFetchBamPE
#'
#' wrapper to handle pileup of standard RNA-seq FR paired end data
#'
#' @param file_paths character vector of file_paths to load from. Alternatively,
#' file_paths can be a data.frame or data.table whose first column is a
#' character vector of paths and additial columns will be used as metadata.
#' @param qgr Set of GRanges to query.  For valid results the width of each
#'   interval should be identical and evenly divisible by \code{win_size}.
#' @param win_size The window size that evenly divides widths in \code{qgr}.
#' @param target_strand character. if one of "+" or "-", reads are filtered to
#'   match. ignored if any other value.
#' @param splice_strategy  character, one of c("none", "ignore", "add", "only",
#'   "splice_count"). Default is "none" and spliced alignment are asssumed not
#'   present. fragLen must be NA for any other value to be valid.  "ignore" will
#'   not count spliced regions.  add" counts spliced regions along with others,
#'   "only" will only count spliced regions and ignore others.
#' @param return_data.table logical. If TRUE the internal data.table is returned
#'   instead of GRanges.  Default is FALSE.
#' @return a GRanges (or data.table if return_data.table == TRUE)
#' @export
#'
#' @examples
#' if(Sys.info()['sysname'] != "Windows"){
#' library(GenomicRanges)
#' bam_f = system.file("extdata/testPE.bam",
#'     package = "seqsetvis", mustWork = TRUE)
#' bam_files = c("a" = bam_f, "b" = bam_f)
#' qgr = CTCF_in_10a_overlaps_gr[1:5]
#' bw_gr = ssvFetchBamPE(bam_files, qgr, win_size = 10)
#' bw_gr
#' }
ssvFetchBamPE.RNA = function(file_paths, qgr, win_size = 50, target_strand = "both", splice_strategy = "ignore",
                         return_data.table = FALSE, win_method = "sample", 
                         flip_strand = FALSE, sum_reads = TRUE, force_skip_centerFix = TRUE){
    y = cn = NULL #reserve bindings
    strand(qgr) = "*"
    bam_r1 = seqsetvis::ssvFetchBam(
        file_paths = file_paths,
        qgr = qgr,
        target_strand = target_strand,
        splice_strategy = splice_strategy,
        return_data.table = TRUE,
        fragLens = NA,
        win_size = win_size,
        flag = scanBamFlag(isFirstMateRead = TRUE),
        win_method = win_method,
        flip_strand = flip_strand, 
        force_skip_centerFix = force_skip_centerFix
    )
    cn = colnames(bam_r1)
    bam_r1$read = "r1"
    bam_r2 =
        seqsetvis::ssvFetchBam(
        # ssvFetchBam.R2(
            file_paths = file_paths,
            qgr = qgr,
            target_strand = target_strand,
            splice_strategy = splice_strategy,
            return_data.table = TRUE,
            fragLens = NA,
            win_size = win_size,
            flag = scanBamFlag(isSecondMateRead = TRUE),
            win_method = win_method,
            flip_strand = !flip_strand,
            force_skip_centerFix = force_skip_centerFix
        )
    bam_r2$read = "r2"

    # ggplot(rbind(bam_r1, bam_r2), aes(x = x, y = y, color = strand)) + geom_path() + facet_wrap("read")
    # bam_r1[, strand := ifelse(strand == "+", "-", "+") ]
    # ggplot(rbind(bam_r1, bam_r2), aes(x = x, y = y, color = strand)) + geom_path() + facet_wrap("read")
    if(sum_reads){
        bam_dt = rbind(bam_r1, bam_r2)[, cn, with = FALSE]
        bam_dt = bam_dt[, list(y = sum(y)), by = c(cn[cn != "y"])][, cn, with = FALSE]
    }else{
        bam_dt = rbind(bam_r1, bam_r2)
    }

    # if(as.character(strand(qgr)) == "+"){
    #     bam_dt[, strand := ifelse(strand == "+", "-", "+") ]
    # }
    # ggplot(bam_dt, aes(x = x, y = y, color = strand)) + geom_path()
    if(!return_data.table){
        bam_dt = GRanges(bam_dt)
    }
    bam_dt
}
#
# ssvFetchBam.R2 = function(file_paths,
#                           qgr,
#                           unique_names = NULL,
#                           win_size = 50,
#                           win_method = c("sample", "summary")[1],
#                           summary_FUN = stats::weighted.mean,
#                           fragLens = "auto",
#                           target_strand = c("*", "+", "-", "both")[1],
#                           flip_strand = FALSE,
#                           anchor = c("left", "left_unstranded", "center",
#                                      "center_unstranded")[3],
#                           names_variable = "sample",
#                           return_data.table = FALSE,
#                           max_dupes = Inf,
#                           splice_strategy = c("none", "ignore", "add",
#                                               "only", "splice_count")[1],
#                           n_cores = getOption("mc.cores", 1),
#                           ...){
#     stopifnot(all(is.character(fragLens) |
#                       is.numeric(fragLens) |
#                       is.na(fragLens)))
#     exp_fragLen = ifelse(is.data.frame(file_paths) || is.data.table(file_paths), nrow(file_paths), length(file_paths))
#     stopifnot(length(fragLens) == 1 || length(fragLens) == exp_fragLen)
#     if(length(fragLens) == 1){
#         if (is.data.frame(file_paths) || is.data.table(file_paths)) {
#             fragLens = rep(fragLens[1], nrow(file_paths))
#         }else{
#             fragLens = rep(fragLens[1], length(file_paths))
#         }
#
#     }
#
#     if (is.data.frame(file_paths) || is.data.table(file_paths)) {
#         if(is.null(colnames(file_paths))){
#             names(fragLens) = file_paths[[1]]
#         }else{
#             if(any(grepl("file", colnames(file_paths)))){
#                 k = which(grepl("file", colnames(file_paths)))[1]
#                 names(fragLens) = file_paths[[k]]
#             }else{
#                 names(fragLens) = file_paths[[1]]
#             }
#         }
#
#
#     }else{
#         names(fragLens) = file_paths
#     }
#
#
#     load_bam = function(f, nam, qgr) {
#         message("loading ", f, " ...")
#         if(!file.exists(paste0(f, ".bai"))){
#             warning("creating index for ", f)
#             Rsamtools::indexBam(f)
#         }
#         fl = fragLens[f]
#         if(!is.na(fl))
#             if(fl == "auto"){
#                 fl = NULL
#             }
#         dt = ssvFetchBam.single.R2(bam_f = f,
#                                 qgr = qgr,
#                                 win_size = win_size,
#                                 win_method = win_method,
#                                 summary_FUN = summary_FUN,
#                                 fragLen = fl,
#                                 target_strand = target_strand,
#                                 anchor = anchor,
#                                 return_data.table = TRUE,
#                                 max_dupes = max_dupes,
#                                 splice_strategy = splice_strategy,
#                                 flip_strand = flip_strand,
#                                 ...)
#         # dt[[names_variable]] = rep(nam, nrow(dt))
#         message("finished loading ", nam, ".")
#         dt
#     }
#
#     bdt = ssvFetchSignal(file_paths = file_paths,
#                          qgr = qgr,
#                          load_signal = load_bam,
#                          unique_names = unique_names,
#                          names_variable = names_variable,
#                          win_size = win_size,
#                          win_method = win_method,
#                          return_data.table = TRUE,
#                          n_cores = n_cores)
#
#     if(flip_strand){
#         if(target_strand != "*"){
#             bdt[, strand := ifelse(strand == "+", "-", "+")]
#         }
#     }
#
#     if(!return_data.table){
#         bdt = GRanges(bdt)
#     }
#     bdt
# }
#
# ssvFetchBam.single.R2 = function(bam_f,
#                                  qgr,
#                                  win_size = 50,
#                                  win_method = c("sample", "summary")[1],
#                                  summary_FUN = stats::weighted.mean,
#                                  fragLen = NULL,
#                                  target_strand = c("*", "+", "-", "both")[1],
#                                  anchor = c("left", "left_unstranded", "center",
#                                             "center_unstranded")[3],
#                                  return_data.table = FALSE,
#                                  max_dupes = Inf,
#                                  splice_strategy = c("none", "ignore", "add",
#                                                      "only", "splice_count")[1],
#                                  flip_strand = FALSE,
#                                  ...) {
#     x = id = y = NULL
#     stopifnot(is.character(win_method))
#     stopifnot(length(win_method) == 1)
#     stopifnot(is(qgr, "GRanges"))
#     stopifnot(win_method %in% c("sample", "summary"))
#     stopifnot(is.function(summary_FUN))
#     stopifnot(target_strand %in% c("*", "+", "-", "both"))
#     stopifnot(anchor %in% c("left", "left_unstranded", "center",
#                             "center_unstranded"))
#     stopifnot(splice_strategy %in% c("none", "ignore", "add",
#                                      "only", "splice_count"))
#     if(splice_strategy == "splice_count"){
#         return(fetchBam.R2(bam_f, qgr, NA, "*",
#                         max_dupes, splice_strategy,
#                         flip_strand = flip_strand, ...))
#     }
#     switch (
#         win_method,
#         sample = {
#             qgr = prepare_fetch_GRanges(qgr, win_size)
#             if(target_strand == "both"){
#                 pos_gr = fetchBam.R2(bam_f, qgr, fragLen, "+",
#                                   max_dupes, splice_strategy,
#                                   flip_strand = flip_strand, ...)
#                 neg_gr = fetchBam.R2(bam_f, qgr, fragLen, "-",
#                                   max_dupes, splice_strategy,
#                                   flip_strand = flip_strand, ...)
#                 pos_dt = viewGRangesWinSample_dt(pos_gr, qgr,
#                                                  win_size, anchor = anchor)
#                 neg_dt = viewGRangesWinSample_dt(neg_gr, qgr,
#                                                  win_size, anchor = anchor)
#                 pos_dt[, strand := "+"]
#                 neg_dt[, strand := "-"]
#                 out = rbind(
#                     pos_dt,
#                     neg_dt
#                 )
#             }else{
#                 score_gr = fetchBam.R2(bam_f, qgr, fragLen, target_strand,
#                                     max_dupes, splice_strategy,
#                                     flip_strand = flip_strand, ...)
#                 out = viewGRangesWinSample_dt(score_gr, qgr,
#                                               win_size, anchor = anchor)
#                 out[, strand := target_strand]
#             }
#
#
#         },
#         summary = {
#             if(target_strand == "both"){
#                 pos_gr = fetchBam.R2(bam_f, qgr, fragLen, "+",
#                                   max_dupes, splice_strategy,
#                                   flip_strand = flip_strand, ...)
#                 neg_gr = fetchBam.R2(bam_f, qgr, fragLen, "-",
#                                   max_dupes, splice_strategy,
#                                   flip_strand = flip_strand, ...)
#                 pos_dt = viewGRangesWinSummary_dt(pos_gr, qgr, win_size,
#                                                   summary_FUN = summary_FUN,
#                                                   anchor = anchor)
#
#                 neg_dt = viewGRangesWinSummary_dt(neg_gr, qgr, win_size,
#                                                   summary_FUN = summary_FUN,
#                                                   anchor = anchor)
#                 pos_dt[, strand := "+"]
#                 neg_dt[, strand := "-"]
#                 out = rbind(
#                     pos_dt,
#                     neg_dt
#                 )
#             }else{
#                 score_gr = fetchBam.R2(bam_f, qgr, fragLen, target_strand,
#                                     max_dupes, splice_strategy,
#                                     flip_strand = flip_strand, ...)
#                 out = viewGRangesWinSummary_dt(score_gr, qgr, win_size,
#                                                summary_FUN = summary_FUN,
#                                                anchor = anchor)
#                 out[, strand := target_strand]
#             }
#         }
#     )
#     # if(any(strand(qgr) == "-")){
#     #     toflip = names(subset(qgr, strand == "-"))
#     #     out[id %in% toflip & strand != "*", strand := ifelse(strand == "+", "-", "+")]
#     # }
#     out = out[order(x)][order(id)][order(strand)]
#
#     if(!return_data.table){
#         out = GRanges(out)
#     }
#     return(out)
# }
#
# fetchBam.R2 = function(bam_f,
#                        qgr,
#                        fragLen = NULL,
#                        target_strand = c("*", "+", "-")[1],
#                        max_dupes = Inf,
#                        splice_strategy = c("none", "ignore", "add",
#                                            "only", "splice_count")[1],
#                        flip_strand = FALSE,
#                        ...){
#     which_label = NULL #reserve binding
#     stopifnot(is.numeric(max_dupes))
#     stopifnot(max_dupes >= 1)
#     if(!is.na(fragLen) && splice_strategy != "none"){
#         stop("fragLen must be NA if splice_strategy is not 'none'.")
#     }
#     if( ! splice_strategy %in% c("none", "ignore", "add",
#                                  "only", "splice_count")){
#         stop('splice_strategy must be one of: "none", "ignore", "add", "only"')
#     }
#     if(is.null(fragLen)){
#         fragLen = fragLen_calcStranded(bam_f, qgr, flip_strand = flip_strand)
#         message("fragLen was calculated as: ", fragLen)
#     }
#     if(!is.na(fragLen)){
#         stopifnot(is.numeric(fragLen))
#         stopifnot(fragLen %% 1 == 0)
#         stopifnot(fragLen >= 1)
#     }
#     sbgr = qgr
#     strand(sbgr) = "*"
#     # browser()
#     sbParam = Rsamtools::ScanBamParam(
#         which = sbgr,
#         what = c("rname", "strand", "pos", "qwidth", "cigar"), ...)
#     bam_raw = Rsamtools::scanBam(bam_f, param = sbParam)
#     bam_dt = lapply(bam_raw, function(x){
#         data.table(seqnames = x$rname, strand = x$strand,
#                    start = x$pos-x$qwidth, width = x$qwidth, cigar = x$cigar)
#     })
#     bam_dt = data.table::rbindlist(bam_dt,
#                                    use.names = TRUE,
#                                    idcol = "which_label")
#     bam_dt = bam_dt[!is.na(width)]
#     toflip = sub(":-", "", as.character(subset(qgr, strand == "-")))
#     if(target_strand %in% c("+", "-")){
#         bam_dt = bam_dt[strand == target_strand & !(which_label %in% toflip) |
#                             strand != target_strand & which_label %in% toflip]
#     }
#
#     bam_dt[, end := start + width - 1L]
#
#     if(max_dupes < Inf){
#         bam_dt = .rm_dupes(bam_dt, max_dupes)
#     }
#     if(is.na(fragLen)){
#         bam_dt = switch(
#             splice_strategy,
#             none = {bam_dt},
#             ignore = {.expand_cigar_dt(bam_dt)},
#             add = {.expand_cigar_dt(bam_dt,
#                                     op_2count = c("M", "D", "=", "X", "N"))},
#             only = {.expand_cigar_dt(bam_dt, op_2count = c("N"))},
#             splice_count = {.expand_cigar_dt(bam_dt, op_2count = c("N"))}
#         )
#     }else{
#         # bam_dt[, end := start + width - 1L]
#         if(flip_strand){
#             bam_dt[strand == "-", end := start + as.integer(fragLen) - 1L]
#             bam_dt[strand == "+", start := end - as.integer(fragLen) + 1L]
#         }else{
#             bam_dt[strand == "+", end := start + as.integer(fragLen) - 1L]
#             bam_dt[strand == "-", start := end - as.integer(fragLen) + 1L]
#         }
#
#     }
#     if(splice_strategy == "splice_count"){
#         return(bam_dt[, .N,
#                       by = list(which_label, seqnames, start, end, strand)])
#     }
#     ext_cov = coverage(split(GRanges(bam_dt), bam_dt$which_label))
#     score_gr = GRanges(ext_cov)
#     if(target_strand == "+"){
#         strand(score_gr) = "+"
#     }
#     if(target_strand == "-"){
#         strand(score_gr) = "-"
#     }
#     return(score_gr)
# }
