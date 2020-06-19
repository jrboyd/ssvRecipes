
#' Title
#'
#' @param file_paths character vector of file_paths to load from. Alternatively,
#' file_paths can be a data.frame or data.table whose first column is a
#' character vector of paths and additial columns will be used as metadata.
#' @param qgr Set of GRanges to query.  For valid results the width of each
#'   interval should be identical and evenly divisible by \code{win_size}.
#' @param cell_cluster_assignments either a list of character vectors containing cell ids or a data.table.  If data.table one column bust be 'id' and contain cell barcodes, the other must match cell_cluster_id_ and contain cluster names.
#' @param cell_cluster_id_ Attribute name of cell_cluster_assignments where cell barcodes are stored.
#' @param id_prefixes Character prefix to prepend to cell barcodes retrieved
#'   from bam files such that they match id items in cell_cluster_assignments.
#'   Either a single item to reuse or character vector of length equal to number
#'   of bam files.  Default of "" adds nothing.
#' @param unique_names not used.  see cell_cluster_assignments for analogous functionality.
#' @param win_size
#' @param win_method
#' @param summary_FUN
#' @param target_strand
#' @param flip_strand
#' @param anchor
#' @param names_variable
#' @param return_data.table
#' @param max_dupes
#' @param splice_strategy
#' @param n_cores
#' @param return_unprocessed
#' @param force_skip_centerFix
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
ssvFetchBam.scRNA = function(file_paths,
                             qgr,
                             cell_cluster_assignments,
                             cell_cluster_id_ = "seurat_clusters",
                             id_prefixes = "",
                             unique_names = NULL,
                             win_size = 50,
                             win_method = c("sample", "summary")[1],
                             summary_FUN = stats::weighted.mean,
                             target_strand = c("*", "+", "-", "both")[1],
                             flip_strand = FALSE,
                             anchor = c("left", "left_unstranded", "center",
                                        "center_unstranded")[3],
                             names_variable = "sample",
                             return_data.table = FALSE,
                             max_dupes = Inf,
                             splice_strategy = c("none", "ignore", "add",
                                                 "only", "splice_count")[1],
                             n_cores = getOption("mc.cores", 1),
                             return_unprocessed = FALSE,
                             force_skip_centerFix = FALSE,
                             ...){

    if (is.data.frame(file_paths) || is.data.table(file_paths)) {
        if (ncol(file_paths) == 1) {
            file_attribs = data.frame(matrix(
                0, nrow = nrow(file_paths), ncol = 0
            ))
        } else{
            file_attribs = file_paths[,-1, drop = FALSE]
        }
        
        file_paths = file_paths[[1]]
        
    } else{
        #file_paths is assumed to be a character vector
        file_attribs = data.frame(data.frame(matrix(
            0, nrow = length(file_paths), ncol = 0
        )))
    }
    if (is.list(file_paths)) {
        file_paths = unlist(file_paths)
    }
    if (is.factor(file_paths))
        file_paths = as.character(file_paths)
    if (!is.null(unique_names)) {
        warnings("unique_names is not used for ssvFetchBam.scRNA and should be NULL. See cell_cluster_assignments description for analogous functionality.")
    }
    stopifnot(is.character(file_paths))
    stopifnot(is(qgr, "GRanges"))
    stopifnot(is.character(names_variable))
    stopifnot(is.numeric(win_size))
    stopifnot(file.exists(file_paths))
    #scRNA checks
    if(is.list(cell_cluster_assignments)){
        cell_cluster_assignments = rbindlist(lapply(cell_cluster_assignments, function(x){
            data.table(id = x)
        }), idcol = cell_cluster_id_)
    }
    setnames(cell_cluster_assignments, cell_cluster_id_, "cell_cluster_id")
    stopifnot("cell_cluster_id" %in% colnames(cell_cluster_assignments))
    stopifnot("id" %in% colnames(cell_cluster_assignments))
    
    
    score_gr_list = parallel::mclapply(seq_along(file_paths), 
                                          function(i){
                                              sc_fetch_pileup(file_paths[i],
                                                              qgr = qgr, 
                                                              cell_cluster_dt = cell_cluster_assignments,  
                                                              id_prefix = id_prfixes[i],
                                                              splice_strategy = splice_strategy)
                                          })
    
    if(win_method == "sample"){
        dt = lapply(score_gr_list, function(score_gr){
            viewGRangesWinSample_dt(score_gr, qgr, window_size = win_size, anchor = anchor)
        }) %>% rbindlist(., idcol = "cell_cluster_id")    
    }else if(win_method == "summary"){
        dt = lapply(score_gr_list, function(score_gr){
            viewGRangesWinSummary_dt(score_gr, qgr, n_tiles = win_size, anchor = anchor, summary_FUN = summary_FUN)
        }) %>% rbindlist(., idcol = "cell_cluster_id")    
    }else{
        stop("invalid win_method, must be one of 'sample' or 'summary'.")
    }
    if(return_data.table){
        return(dt[])
    }else{
        return(GRanges(dt))
    }
    
    
}


sc_fetch_pileup = function(bam_file, qgr, 
                           cell_cluster_dt,
                           id_prefix = NULL,
                           splice_strategy = c("none", "ignore", "add", "only", "splice_count")[2]){
    what = c("rname", "strand", "pos", "qwidth", "cigar", "CR")
    what = Rsamtools::scanBamWhat()
    sbParam = Rsamtools::ScanBamParam(
        which = qgr,
        what = what, 
        tag = c("CR", "UB"))
    bam_raw = Rsamtools::scanBam(bam_file, param = sbParam)
    
    bam_dt = lapply(bam_raw, function(x){
        data.table(seqnames = x$rname, strand = x$strand,
                   # start = x$pos, width = x$qwidth, cigar = x$cigar, id = x$tag$CR, umi = x$tag$UB, seq = as.character(x$seq))
                   start = x$pos, width = x$qwidth, cigar = x$cigar, id = x$tag$CR, umi = x$tag$UB)
    })
    bam_dt = data.table::rbindlist(bam_dt,
                                   use.names = TRUE,
                                   idcol = "which_label")
    
    bam_dt = bam_dt[!is.na(width)]
    toflip = sub(":-", "", as.character(subset(qgr, strand == "-")))
    target_strand = "+"
    if(target_strand %in% c("+", "-")){
        bam_dt = bam_dt[strand == target_strand & !(which_label %in% toflip) |
                            strand != target_strand & which_label %in% toflip]
    }
    if(!is.null(id_prefix)){
        bam_dt$id = paste0(id_prefix, bam_dt$id)
    }
    
    bam_dt[, end := start + width - 1L]
    bam_dt = bam_dt %>% unique
    
    cig_dt = bam_dt[, .(which_label = paste(which_label, id, umi), cigar, seqnames, start, strand, width)]
    .expand_cigar_dt = seqsetvis:::.expand_cigar_dt
    .expand_cigar_dt_recursive = seqsetvis:::.expand_cigar_dt_recursive
    cig_dt = switch(
        splice_strategy,
        none = {cig_dt},
        ignore = {.expand_cigar_dt(cig_dt)},
        add = {.expand_cigar_dt(cig_dt,
                                op_2count = c("M", "D", "=", "X", "N"))},
        only = {.expand_cigar_dt(cig_dt, op_2count = c("N"))},
        splice_count = {.expand_cigar_dt(cig_dt, op_2count = c("N"))}
    )
    
    cig_flat_cell = split(GRanges(cig_dt), paste(cig_dt$which_label)) %>% reduce %>% unlist
    
    # browser()
    cig_flat_cell_dt = as.data.table(cig_flat_cell)
    cig_flat_cell_dt$which_label = names(cig_flat_cell)
    
    
    cig_flat_cell_dt[, c('which_label', "id", "umi") := tstrsplit(which_label, " ")]
    cig_flat_cell_dt = merge(cig_flat_cell_dt, cell_cluster_dt[, .(id, cell_cluster_id)], by = "id")
    cig_flat_cell_dt[, which_label := paste(which_label, cell_cluster_id)]
    #reduce umi to pile of 1
    
    mc_uniq = cig_flat_cell_dt$cell_cluster_id %>% unique
    ext_cov_list = lapply(mc_uniq, function(mc){
        coverage(split(GRanges(cig_flat_cell_dt[cell_cluster_id == mc]), cig_flat_cell_dt[cell_cluster_id == mc]$which_label))
    })
    names(ext_cov_list) = mc_uniq
    ext_cov = coverage(split(GRanges(cig_flat_cell_dt), cig_flat_cell_dt$which_label))
    score_gr_list = lapply(ext_cov_list, function(ext_cov){
        score_gr = GRanges(ext_cov)
        score_gr = harmonize_seqlengths(score_gr, bam_f)
        
        if(length(score_gr) == 0){
            score_gr = sbgr
            mcols(score_gr) = NULL
            score_gr$score = 0
        }
        if(target_strand == "+"){
            strand(score_gr) = "+"
        }
        if(target_strand == "-"){
            strand(score_gr) = "-"
        }
        return(score_gr)
    })
    score_gr_list
}

# 
# qgr = sc_qgr[delta_assign_dt[cluster_id == 1]$id]
# subset(tx_gr, gene_name == "Psmb4")
# 
# 
# sc_score_gr_list_df_pile = sc_fetch_pileup(sc_bam_files[1], qgr = qgr, id_prefix = "df_", splice_strategy = "ignore")
# sc_score_gr_list_wt_pile = sc_fetch_pileup(sc_bam_files[2], qgr = qgr, id_prefix = "wt_", splice_strategy = "ignore")
# 
# sc_score_gr_list_df = sc_fetch_pileup(sc_bam_files[1], qgr = qgr, id_prefix = "df_", splice_strategy = "add")
# sc_score_gr_list_df_spliced = sc_fetch_pileup(sc_bam_files[1], qgr = qgr, id_prefix = "df_", splice_strategy = "only")
# sc_score_gr_list_wt = sc_fetch_pileup(sc_bam_files[2], qgr = qgr, id_prefix = "wt_", splice_strategy = "add")
# sc_score_gr_list_wt_spliced = sc_fetch_pileup(sc_bam_files[2], qgr = qgr, id_prefix = "wt_", splice_strategy = "only")
# 
# sc_view = function(score_gr_list, type = c("pileup", "splicing")[1], status = c("df", "wt")){
#     dt = lapply(score_gr_list, function(score_gr){
#         viewGRangesWinSample_dt(score_gr, qgr, window_size = 50)
#     }) %>% rbindlist(., idcol = "cell_cluster_id")
#     dt$type = type
#     dt$status = status
#     dt[]
# }
# 
# sc_view_dt_pile = rbind(
#     sc_view(sc_score_gr_list_df_pile, "pileup", "df4"),
#     sc_view(sc_score_gr_list_wt_pile, "pileup", "wt")
# )