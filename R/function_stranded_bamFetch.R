#' Title
#'
#' @param bam_files
#' @param qgr
#' @param win_size
#' @param fragLens
#' @param return_data.table
#' @param flipStrand
#'
#' @return
#' @export
#' @import seqsetvis, GenomicRanges, data.table
#'
#' @examples
myFetchStrandedBam = function(bam_files, qgr, win_size = 5, fragLens = NA, return_data.table = TRUE, flipStrand = FALSE, max_dupes = Inf){
    pos_dt = ssvFetchBam(bam_files,
                         return_data.table = return_data.table,
                         qgr = qgr,
                         win_size = win_size,
                         target_strand = "+",
                         fragLens = fragLens,
                         max_dupes = max_dupes)
    neg_dt = ssvFetchBam(bam_files,
                         return_data.table = return_data.table,
                         qgr = qgr,
                         win_size = win_size,
                         target_strand = "-",
                         fragLens = fragLens,
                         max_dupes = max_dupes)
    if(flipStrand){
        pos_dt$strand = "-"
        neg_dt$strand = "+"
    }else{
        pos_dt$strand = "+"
        neg_dt$strand = "-"
    }

    rbind(pos_dt, neg_dt)
}

myStrandDiff = function(stranded_dt, cap = 30){
    stranded_dtw = dcast(stranded_dt, id + x + sample ~ strand, value.var = "y")
    stranded_dtw[, diff := `+` - `-`]

    hist(stranded_dtw$diff, breaks = 500, xlim = c(-cap, cap))
    stranded_dtw[, cap_diff := diff]

    stranded_dtw[cap_diff > cap, cap_diff := cap]
    stranded_dtw[cap_diff < -cap, cap_diff := -cap]
    hist(stranded_dtw$cap_diff, breaks = cap, xlim = c(-cap, cap))
    stranded_dtw
}



myPlotDiffClusters = function(diff_dt, clust_dt, fcap = 5){
    all_blocked_dc = merge(diff_dt, clust_dt, by = "id")

    all_blocked_dc[, cap_fill := cap_diff]
    all_blocked_dc[cap_fill > fcap, cap_fill := fcap]
    all_blocked_dc[cap_fill < -fcap, cap_fill := -fcap]

    # ggplot(all_blocked_dc,
    #        aes(x = x, y = id, fill = cap_fill)) +
    #     geom_raster() +
    #     theme(axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
    #     scale_fill_gradientn(colors = c("darkblue", "white", "darkred"), limits = c(-fcap,fcap)) +
    #     facet_wrap("sample", nrow = 1)

    pFullHeatmap = ssvSignalHeatmap(all_blocked_dc, fill_ = "cap_fill", cluster_ = "cluster_id")

    med_all_blocked_dc = all_blocked_dc[, .(mval = mean(cap_diff)), by = .(cluster_id, x, sample)]
    # hist(med_all_blocked_dc$mval, breaks = 50)

    med_all_blocked_dc[, str_clust_id := factor(cluster_id)]
    med_all_blocked_dc[, str_clust_id := factor(str_clust_id, levels = rev(levels(str_clust_id)))]

    #heatmap agg by cluster
    pAggHeatmap = ggplot(med_all_blocked_dc, aes(x = x, y = str_clust_id, fill = mval)) +
        geom_raster() + facet_wrap("sample", nrow = 1) +
        scale_fill_gradientn(colors = c("darkblue", "white", "darkred")) + scale_y_discrete(scale = )

    #lineplot agg by cluster
    pAggLineplot = ggplot(med_all_blocked_dc, aes(x = x, y = mval)) +
        geom_path() +
        facet_grid("cluster_id~sample") +
        scale_fill_gradientn(colors = c("darkblue", "white", "darkred")) +
        # scale_y_reverse() +
        annotate("line", x = c(0,0), y = c(-7,7), color = "red")
    return(list(detailHeatmap = pFullHeatmap, aggHeatmap = pAggHeatmap, aggLineplot = pAggLineplot ))
}

# load('results_follow_good_qgr_pipeline_10a_all_blocked_comparison.save')
#
# root = "/slipstream/galaxy/uploads/working/qc_framework"
# all_blocked = c(
#     # feb_blocked_R1 = "output_MinguJosh_feb2018/MCF10A-blocked_Runx1_R1/MCF10A-blocked_Runx1_R1_FE.bw",
#     # feb_blocked_R2 = "output_MinguJosh_feb2018/MCF10A-blocked_Runx1_R2/MCF10A-blocked_Runx1_R2_FE.bw",
#     # jun23980_blocked_R1 = "output_JR_bookmarking_blocked_RUNX1/MCF10A-blocked_Runx1-23980_R1/MCF10A-blocked_Runx1-23980_R1_FE.bw",
#     # jun23980_blocked_R2 = "output_JR_bookmarking_blocked_RUNX1/MCF10A-blocked_Runx1-23980_R2/MCF10A-blocked_Runx1-23980_R2_FE.bw",
#     # jun4336BF_blocked_R1 = "output_JR_bookmarking_blocked_RUNX1/MCF10A-blocked_Runx1-4336BF_R1/MCF10A-blocked_Runx1-4336BF_R1_FE.bw",
#     # jun4336BF_blocked_R2 = "output_JR_bookmarking_blocked_RUNX1/MCF10A-blocked_Runx1-4336BF_R2/MCF10A-blocked_Runx1-4336BF_R2_FE.bw",
#     # july_blocked_R1 = "output_JR_bookmarking_full_RUNX1/MCF10A-blocked_Runx1_R1/MCF10A-blocked_Runx1_R1_FE.bw",
#     # july_blocked_R2 = "output_JR_bookmarking_full_RUNX1/MCF10A-blocked_Runx1_R2/MCF10A-blocked_Runx1_R2_FE.bw",
#     blocked_R1 = "output_JR_bookmarking_full_RUNX1/MCF10A-blocked_Runx1_R1/MCF10A-blocked_Runx1_R1_FE.bw",
#     blocked_R2 = "output_JR_bookmarking_full_RUNX1/MCF10A-blocked_Runx1_R2/MCF10A-blocked_Runx1_R2_FE.bw",
#     released_R1 = "output_JR_bookmarking_full_RUNX1/MCF10A-released_Runx1_R1/MCF10A-released_Runx1_R1_FE.bw",
#     released_R2 = "output_JR_bookmarking_full_RUNX1/MCF10A-released_Runx1_R2/MCF10A-released_Runx1_R2_FE.bw",
#     dmso_R1 = "output_JR_bookmarking_full_RUNX1/MCF10A-dmso_Runx1_R1/MCF10A-dmso_Runx1_R1_FE.bw",
#     dmso_R2 = "output_JR_bookmarking_full_RUNX1/MCF10A-dmso_Runx1_R2/MCF10A-dmso_Runx1_R2_FE.bw"
# )
# all_blocked_bw = file.path(root, all_blocked)
# names(all_blocked_bw) = names(all_blocked)
# all_blocked_bam = sub("_FE.bw", ".bam", all_blocked_bw)
#
# bam_files = all_blocked_bam
# qgr = really_good_qgr
#
# bam_dt = myFetchStrandedBam(all_blocked_bam, really_good_qgr)
# diff_dt = myStrandDiff(bam_dt)
# diff_plots = myPlotDiffClusters(diff_dt, clust_dt, fcap = 5)
#
#
# diff_plots$detailHeatmap
# diff_plots$aggHeatmap
# diff_plots$aggLineplot
