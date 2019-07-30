#' returns correlation scores for all regions in qgr as specified fragLens
#'
#' @param bam_file
#' @param qgr
#' @param fragLens
#'
#' @return data.table with id and correlation
#' @export
#'
#' @examples
corr_fragLens.single = function(bam_file, qgr, fragLens = 200){
    bam_dt = myFetchStrandedBam(bam_file, qgr, return_data.table = TRUE, fragLens = fragLens)
    s_dt = dcast(bam_dt, id + x + sample ~ strand, value.var = "y")
    cor_dt = s_dt[, cor(`+`, `-`), by = .(id)]
}

#' compares correlation of specified fragLens to un-extended read pileup
#'
#' @param bam_file
#' @param qgr
#' @param fragLens
#'
#' @return data.table merge with qgr with read and fragment correlation
#' @export
#'
#' @examples
corr_fragLens = function(bam_file, qgr, fragLens = 200){
    bam_noExt = corr_fragLens.single(bam_file, qgr, fragLens = NA)
    bam_wExt = corr_fragLens.single(bam_file, qgr, fragLens = fragLens)

    qdt = as.data.table(qgr)
    if(is.null(qdt$id)){
        qdt$id = bam_noExt$id
    }

    mdt = merge(bam_noExt, bam_wExt)
    colnames(mdt)[2:3] = c("read", "fragment")
    mdt[, diff := fragment - read]
    merge(qdt, mdt, by = "id")
}

#' Plots knowresults from homer results directory
#'
#' @param homer_dir root directory of homer results
#' @param dir_regex regex to match sub-directories
#' @param target_consensus consensus to match from homer. default is AAACCACARM for RUNX1.
#'
#' @return
#' @export
#' @import magrittr
#'
#' @examples
plotHomerKnowResults = function(homer_dir = "homer_strand_corr_results",
                       dir_regex = "enriched.+Sorted",
                       target_consensus = "AAACCACARM"){
    knownRes = dir(homer_dir, pattern = dir_regex, full.names = TRUE) %>%
        sapply(., function(x)dir(x, pattern = "knownResults.txt", full.names = TRUE))
    knownRes = knownRes[lengths(knownRes) > 0]
    names(knownRes) = sub(".+enriched_motifs_", "", names(knownRes))
    knownRes_dt = lapply(knownRes, fread)
    knownRes_dt = lapply(knownRes_dt, function(x){
        colnames(x)[c(6,8)] = c("# of Target Sequences with Motif", "# of Background Sequences with Motif")
        x
    })
    # lapply(knownRes_dt, dim)
    knownRes_dt = rbindlist(knownRes_dt, use.names = TRUE, idcol = "group")

    knownRes_dt[, c("treatment", "corr_group", "metric") := tstrsplit(group, "[_\\.]")]
    knownRes_dt$corr_group = factor(knownRes_dt$corr_group, levels = c("worst", "worse", "bad", "neutral", "good", "better", "best"))
    knownRes_dt$treatment = factor(knownRes_dt$treatment, levels = c("blocked", "released", "dmso"))
    knownRes_dt[, fraction := as.numeric(sub("%", "", `% of Target Sequences with Motif`))/100]
    p = ggplot(knownRes_dt[Consensus == target_consensus][order(corr_group)],
           aes(x = corr_group,
               y = fraction,
               color = treatment,
               group = treatment)) +
        geom_path() +
        geom_point() +
        labs(x = "strand correlation quality",
             y = "fraction of sites with Runx1 motif (homer)",
             title = "Runx1 motif and (fragment - read) strand correlation") +
        facet_wrap("metric", nrow = 1) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
    return(list(plot = p, homer_results = knownRes_dt))
}
