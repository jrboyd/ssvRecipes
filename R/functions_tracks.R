
#' Title
#'
#' @param bw_files
#' @param qgr
#' @param flip_x
#' @param nwin
#' @param win_FUN
#' @param nspline
#' @param debug
#' @param show_lines
#' @param show_fill
#'
#' @return
#' @export
#'
#' @examples
track_chip = function(bw_files, qgr,
                      flip_x = FALSE, nwin = 100,
                      win_FUN = c("mean", "max")[1],
                      nspline = 10, debug = FALSE,
                      show_lines = TRUE, show_fill = FALSE){
    if(debug) browser()
    stopifnot(win_FUN %in% c("mean", "max"))
    strand(qgr) = "*"
    rng = c(start(qgr), end(qgr))

    sum_FUN = switch (win_FUN,
                      max = function(x, w)max(x),
                      mean = weighted.mean
    )

    bw_dt = ssvFetchBigwig(bw_files, qgr,
                           win_method = "summary",  win_size = nwin,
                           summary_FUN = sum_FUN,
                           return_data.table = TRUE, anchor = "left")
    bw_dt$sample = factor(bw_dt$sample, levels = names(bw_files))
    bw_dt[grepl("logFE", sample), y := 10^y]
    bw_dt[y < 1, y := 1]

    if(nspline > 1){
        bw_dt = applySpline(bw_dt, n = nspline, by_ = c("sample", "id"))
        bw_dt[, x := (seq_len(.N)-.5) / (.N) * (max(end) - min(start)) + min(start), by = c("sample", "id")]
        bw_dt[y < 1, y := 1]
    }
    # bw_dt = applySpline(bw_dt, n = nspline, by_ = c("sample", "id"))



    p_chip = ggplot(bw_dt) #coord_cartesian(xlim = c(2.5e4, 3.3e4)) +
    if(show_lines){
        p_chip = p_chip + geom_path(aes_string(x = "x", y = "y", color = "sample"))
    }
    if(show_fill){
        p_chip = p_chip + geom_ribbon(aes_string(x = "x", ymin = 0, ymax = "y", fill = "sample"))
    }
    p_chip = p_chip +        # scale_x_reverse(labels = function(x)x/10^3, limits = rev(rng)) +
        guides(color = "none", fill = "none") +
        labs(x = "kbp", y = "FE") +
        facet_wrap("sample", ncol = 1) +
        theme_classic() +
        theme(strip.background = element_blank())
    if(flip_x){
        p_chip = p_chip +
            scale_x_reverse(labels = function(x)x/10^3, limits = rev(rng))

    }else{
        p_chip = p_chip +
            scale_x_continuous(labels = function(x)x/10^3, limits = rng)
    }
    p_chip
}

#' Title
#'
#' @param bams
#' @param qgr
#' @param flip_x
#' @param flip_strand
#' @param max_dupes
#' @param win_size
#' @param strand_upsidedown
#'
#' @return
#' @export
#'
#' @examples
track_rna = function(bams, qgr, flip_x = FALSE, flip_strand = FALSE, max_dupes = Inf, win_size = 100, strand_upsidedown = TRUE){
    # bams = c("NS gapmer" = "~/R/KZ_P01_runx_binding_KD_DE/bams/MDA231_NSgapmer_merged.bam",
    #          "MANCR gapmer" = "~/R/KZ_P01_runx_binding_KD_DE/bams/MDA231_MANCRgapmer_merged.bam")
    # ,
    # "MDA231 NSgapmer" = "~/R/KZ_P01_runx_binding_KD_DE/bams/MDA231_NSgapmer_merged.bam",
    # "MDA231 MANCRgapmer" = "~/R/KZ_P01_runx_binding_KD_DE/bams/MDA231_MANCRgapmer_merged.bam")

    bam_counts = paste0(bams, ".count")
    names(bam_counts) = names(bams)
    bam_counts = sapply(bam_counts, function(f)read.table(f)[1,1])

    bam_dt = ssvRecipes::myFetchStrandedBam(bams, qgr,
                                            return_data.table = TRUE, anchor = "left",
                                            win_size = win_size, flipStrand = flip_strand,
                                            splice_strategy = "add", max_dupes = max_dupes)
    for(nam in names(bam_counts)){
        bam_dt[sample == nam, ynorm := y / bam_counts[nam] * mean(bam_counts)]
    }


    splice_dt = ssvRecipes::myFetchStrandedBam(bams, qgr,
                                               return_data.table = TRUE, anchor = "left",
                                               win_size = win_size, flipStrand = flip_strand,
                                               splice_strategy = "only", max_dupes = max_dupes)

    for(nam in names(bam_counts)){
        splice_dt[sample == nam, ynorm := y / bam_counts[nam] * mean(bam_counts)]
    }



    nbam_dt = copy(bam_dt)
    nsplice_dt = copy(splice_dt)

    if(strand_upsidedown){
        nbam_dt[strand == "+", ynorm := -ynorm]
        nsplice_dt[strand == "+", ynorm := -ynorm]
    }


    p_rna = ggplot() + #coord_cartesian(xlim = c(5e4, 7e4)) +
        geom_ribbon(data = nsplice_dt, aes(x = (start + end) / 2, ymin = 0, ymax = ynorm,
                                           fill = sample, group = sample), alpha = .4) +
        geom_path(data = nbam_dt, aes(x = (start + end) / 2, y = ynorm, color = sample, group = sample), size = 1) +
        scale_fill_manual(values = c("black", "red")) +
        scale_color_manual(values = c("black", "darkred")) +
        guides(fill = "none") +
        # scale_x_reverse(labels = function(x)x/10^3, limits = rev(rng)) +
        labs(x = "", y = "read pileup\n(stranded)") +
        # facet_wrap("sample", ncol = 1) +
        theme_classic() +
        theme(strip.background = element_blank(), legend.position = "right") #+
    # annotate("line", x = rng, y = c(0,0), color = "gray", size = 2)
    if(!strand_upsidedown){
        p_rna = p_rna +  facet_wrap("strand", ncol = 1)
    }
    # p_rna

    if(flip_x){
        p_rna = p_rna + scale_x_reverse(labels = function(x)x/10^3, limits = rev(rng))
    }else{
        p_rna = p_rna + scale_x_continuous(labels = function(x)x/10^3, limits = rng)
    }
    p_rna
}

#' Title
#'
#' @param ref
#' @param qgr
#' @param flip_x
#' @param exon_height
#' @param debug
#'
#' @return
#' @export
#'
#' @examples
track_ref = function(ref = "~/gencode.v28.annotation.gtf.gz", qgr, flip_x = FALSE,
                     exon_height = .5, intron_thickness = 2, debug = FALSE){
    if(debug) browser()
    if(!class(ref) == "GRanges"){
        if(file.exists(ref)){
            ref = rtracklayer::import.gff(ref, format =  "gtf", feature.type = "exon")
        }else{
            stop("ref must be gtf loaded as GRanges or path to gtf")
        }
    }
    rng = c(start(qgr), end(qgr))
    # ref_gr = rtracklayer::import.gff("~/gencode.v28.annotation.gtf.gz", format =  "gtf", feature.type = "exon")
    ref_dt = as.data.table(subsetByOverlaps(ref, qgr, ignore.strand = TRUE))
    # ref_dt = ref_dt[gene_name %in% c("LINC00704", "LINC00705")]
    yvar = "gene_name"

    ref_dt[gene_name == "RP11-117P22.1", gene_name := "MANCR"]
    ref_dt[gene_name == "RP11-117P22.2", gene_name := "LINC00705"]
    ref_dt[[yvar]] = factor(ref_dt[[yvar]])
    ref_dt$ymin = as.numeric(ref_dt[[yvar]]) - .5 * exon_height
    ref_dt$ymax = as.numeric(ref_dt[[yvar]]) + .5 * exon_height


    ref_dt_base = ref_dt[, .(start = min(start), end = max(end), y = mean(c(ymin, ymax)), strand = unique(strand)), by = yvar]
    ref_dt_base = melt(ref_dt_base, id.vars = c("strand", "gene_name", "y"), value.name = "x")

    p_ref = ggplot() +
        geom_line(data = ref_dt_base, aes_string(x = "x", y = "y", color = "strand", group = "gene_name"), size = intron_thickness) +
        geom_rect(data = ref_dt, aes(fill = strand, color = strand, xmin = start, xmax = end, ymin = ymin, ymax = ymax)) +
        scale_y_continuous(breaks = seq_along(levels(ref_dt[[yvar]])),
                           labels = function(x)levels(ref_dt[[yvar]])[round(x)]) +
        # scale_x_reverse(labels = function(x)x/10^3, limits = rev(rng)) +
        theme_classic() +
        labs(y = "gene\nannotation", x = "kbp") #+
    # scale_fill_manual(values = c("-" = "black", "+" = "darkgray")) +
    # scale_color_manual(values = c("-" = "black", "+" = "darkgray"))
    if(flip_x){
        p_ref = p_ref +
            scale_x_reverse(labels = function(x)x/10^3) +
            coord_cartesian(xlim = rev(rng )) +
            scale_fill_manual(values = c("-" = "black", "+" = "darkgray")) +
            scale_color_manual(values = c("-" = "black", "+" = "darkgray"))

    }else{
        p_ref = p_ref +
            scale_x_continuous(labels = function(x)x/10^3) +
            coord_cartesian(xlim = rng ) +
            scale_fill_manual(values = c("+" = "black", "-" = "darkgray")) +
            scale_color_manual(values = c("+" = "black", "-" = "darkgray"))
    }
    tss_dt = rbind(
        ref_dt[strand == "+", .(tss =min(start), y = mean(c(ymin, ymax)), strand = unique(strand)), by = "transcript_id"],
        ref_dt[strand == "-", .(tss =max(end), y = mean(c(ymin, ymax)), strand = unique(strand)), by = "transcript_id"]
    )
    tss_dt
    p_ref
}

# KT_fig = function(qgr, pdf_name){
#     rng = c(start(qgr), end(qgr))
#     # qgr = GRanges(gsub(",", "", "chr10:4,700,808-4,708,322"))
#
#     flip_x = FALSE
#
#     strand(qgr) = "-"
#     runx1_bw = c(
#         #`MDA231 Runx1` = "/slipstream/galaxy/uploads/working/qc_framework/output_MK_MDA231_Runx/MDA231_Runx1_pooled/MDA231_Runx1_pooled_FE.bw",
#         `MDA231 Runx2` = "/slipstream/galaxy/uploads/working/qc_framework/output_MK_MDA231_Runx/MDA231_Runx2_pooled/MDA231_Runx2_pooled_FE.bw",
#         `MDA231 MANCR` = "bams/MDA231_MANCER_GRIDseq.bw"
#     )
#     # ,
#     #              "/slipstream/galaxy/uploads/working/qc_framework/output_MDA231_RunxChIP/MDA231_RUNX1_pooled/MDA231_RUNX1_pooled_logFE.bw",
#     #              "/slipstream/galaxy/uploads/working/qc_framework/output_MDA231_RunxChIP/MDA231_RUNX2_pooled/MDA231_RUNX2_pooled_logFE.bw")
#
#     # bw_dt = ssvFetchBigwig(runx1_bw, qgr, win_method = "summary",  win_size = 100, summary_FUN = function(x, w)mean(x), return_data.table = TRUE)
#     # bw_dt[, ymov := seqsetvis:::movingAverage(y, 5), by = .(sample)]
#     # bw_dt[ y < 1, y := 1]
#     # bw_dt = applySpline(bw_dt, 3)
#     # bw_dt[ y < 1, y := 1]
#     # ssvSignalLineplot(bw_dt, y_ = "y")
#
#     p_chip = track_chip(runx1_bw[1], qgr = qgr, nspline = 10, nwin = 80, win_FUN = "max",
#                         show_fill = TRUE, show_lines = FALSE, flip_x = flip_x) +
#         labs(x = "", y = "Runx2 FE", title = "") + scale_fill_manual(values = "purple") +
#         theme(plot.margin = margin(0, 0, 0, 0, unit = "cm"))
#     p_grid = track_chip(runx1_bw[2], qgr = qgr, nspline = 10, nwin = 80, win_FUN = "max",
#                         show_fill = TRUE, show_lines = FALSE, flip_x = flip_x) +
#         labs(x = "", y = "MANCR GRID-seq", title = "") + scale_fill_manual(values = "forestgreen") +
#         theme(plot.margin = margin(0, 0, 0, 0, unit = "cm"))
#
#     p_chip = p_chip + theme(strip.text = element_blank())
#     p_grid = p_grid + theme(strip.text = element_blank())
#     library(GenomicRanges)
#     library(data.table)
#     library(seqsetvis)
#     bams = c("NS gapmer" = "~/R/KZ_P01_runx_binding_KD_DE/bams/MDA231_NSgapmer_merged.bam",
#              "MANCR gapmer" = "~/R/KZ_P01_runx_binding_KD_DE/bams/MDA231_MANCRgapmer_merged.bam")
#     # ,
#     # "MDA231 NSgapmer" = "~/R/KZ_P01_runx_binding_KD_DE/bams/MDA231_NSgapmer_merged.bam",
#     # "MDA231 MANCRgapmer" = "~/R/KZ_P01_runx_binding_KD_DE/bams/MDA231_MANCRgapmer_merged.bam")
#
#     bam_counts = paste0(bams, ".count")
#     names(bam_counts) = names(bams)
#     bam_counts = sapply(bam_counts, function(f)read.table(f)[1,1])
#
#     bam_dt = ssvRecipes::myFetchStrandedBam(bams, qgr,
#                                             return_data.table = TRUE, anchor = "left",
#                                             win_size = 100, flipStrand = TRUE,
#                                             splice_strategy = "add", max_dupes = 5)
#     for(nam in names(bam_counts)){
#         bam_dt[sample == nam, ynorm := y / bam_counts[nam] * mean(bam_counts)]
#     }
#
#
#     splice_dt = ssvRecipes::myFetchStrandedBam(bams, qgr,
#                                                return_data.table = TRUE, anchor = "left",
#                                                win_size = 100, flipStrand = TRUE,
#                                                splice_strategy = "only", max_dupes = 5)
#
#     for(nam in names(bam_counts)){
#         splice_dt[sample == nam, ynorm := y / bam_counts[nam] * mean(bam_counts)]
#     }
#
#
#     strand_upsidedown = TRUE
#
#     nbam_dt = copy(bam_dt)
#     nsplice_dt = copy(splice_dt)
#
#     if(strand_upsidedown){
#         nbam_dt[strand == "+", ynorm := -ynorm]
#         nsplice_dt[strand == "+", ynorm := -ynorm]
#     }
#
#
#     p_rna = ggplot() + #coord_cartesian(xlim = c(5e4, 7e4)) +
#         geom_ribbon(data = nsplice_dt, aes(x = (start + end) / 2, ymin = 0, ymax = ynorm,
#                                            fill = sample, group = sample), alpha = .4) +
#         geom_path(data = nbam_dt, aes(x = (start + end) / 2, y = ynorm, color = sample, group = sample), size = 1) +
#         scale_fill_manual(values = c("black", "red")) +
#         scale_color_manual(values = c("black", "darkred")) +
#         guides(fill = "none") +
#         # scale_x_reverse(labels = function(x)x/10^3, limits = rev(rng)) +
#         labs(x = "", y = "read pileup\n(stranded)") +
#         # facet_wrap("sample", ncol = 1) +
#         theme_classic() +
#         theme(strip.background = element_blank(), legend.position = "right") #+
#     # annotate("line", x = rng, y = c(0,0), color = "gray", size = 2)
#     if(!strand_upsidedown){
#         p_rna = p_rna +  facet_wrap("strand", ncol = 1)
#     }
#     # p_rna
#
#     if(flip_x){
#         p_rna = p_rna + scale_x_reverse(labels = function(x)x/10^3, limits = rev(rng))
#     }else{
#         p_rna = p_rna + scale_x_continuous(labels = function(x)x/10^3, limits = rng)
#     }
#
#     # p_list = list(p_rna, p_chip)
#     # p_list = dthic::sync_width(p_list)
#     # cowplot::plot_grid(plotlist = p_list, ncol = 1)
#
#     if(!exists("ref_gr"))
#         ref_gr = rtracklayer::import.gff("~/gencode.v28.annotation.gtf.gz", format =  "gtf", feature.type = "exon")
#     ref_dt = as.data.table(subsetByOverlaps(ref_gr, qgr, ignore.strand = TRUE))
#     # ref_dt = ref_dt[gene_name %in% c("LINC00704", "LINC00705")]
#     yvar = "gene_name"
#
#     min_w = round(width(qgr) / 500)
#
#     # ref_dt[gene_name == "RP11-117P22.1", gene_name := "MANCR"]
#     # ref_dt[gene_name == "RP11-117P22.2", gene_name := "LINC00705"]
#     # ref_dt = ref_dt[transcript_support_level <= 2 & gene_name %in% c("RP1-186E20.2", "CELF2")]
#     ref_dt = ref_dt[transcript_id %in% c("ENST00000542579.5", "ENST00000446372.2")]
#     ref_dt[[yvar]] = factor(ref_dt[[yvar]])
#
#     ref_dt[[yvar]] = droplevels(ref_dt[[yvar]])
#
#
#     ref_dt$ymin = as.numeric(ref_dt[[yvar]]) - .5
#     ref_dt$ymax = as.numeric(ref_dt[[yvar]]) + .5
#
#     ref_dt[end - start < min_w, c("start", "end") := .(start - as.integer(ceiling((min_w - (end - start))/2)),
#                                                        end + as.integer(ceiling((min_w - (end - start))/2)))]
#
#
#
#     ref_dt_base = ref_dt[, .(start = min(start), end = max(end), y = mean(c(ymin, ymax)), strand = unique(strand)), by = yvar]
#     ref_dt_base = melt(ref_dt_base, id.vars = c("strand", yvar, "y"), value.name = "x")
#
#     p_ref = ggplot() +
#         geom_line(data = ref_dt_base, aes_string(x = "x", y = "y", color = "strand", group = yvar), size = 2) +
#         geom_rect(data = ref_dt, aes(fill = strand, xmin = start, xmax = end, ymin = ymin, ymax = ymax)) +
#         scale_y_continuous(breaks = seq_along(levels(ref_dt[[yvar]])),
#                            labels = function(x)levels(ref_dt[[yvar]])[round(x)]) +
#         # scale_x_reverse(labels = function(x)x/10^3, limits = rev(rng)) +
#         theme_classic() +
#         labs(x = "", y = "gene\nannotation")
#
#     if(flip_x){
#         p_ref = p_ref + scale_x_reverse(labels = function(x)x/10^3, limits = rev(rng))+
#             scale_fill_manual(values = c("-" = "black", "+" = "darkgray")) +
#             scale_color_manual(values = c("-" = "black", "+" = "darkgray"))
#     }else{
#         p_ref = p_ref + scale_x_continuous(labels = function(x)x/10^3, limits = rng)+
#             scale_fill_manual(values = c("+" = "black", "-" = "darkgray")) +
#             scale_color_manual(values = c("+" = "black", "-" = "darkgray"))
#     }
#     # p_ref
#
#     p_list = list(p_ref, p_rna + ylim(-100, 300), p_chip, p_grid)
#     o = c(3, 4, 2, 1)
#     p_list = p_list[o]
#     p_list = lapply(p_list, function(p){
#         p +
#             theme(
#                 # legend.text = element_text(size = 8),
#                 # legend.title = element_text(size = 10),
#                 # axis.text = element_text(size = 8),
#                 text = element_text(size = 8)
#             )
#     })
#
#     p_list[[length(p_list)]] = p_list[[length(p_list)]] + labs(x = paste(as.character(seqnames(qgr)), "kbp"))
#     p_list = sync_width(p_list)
#     pg = cowplot::plot_grid(plotlist = p_list, ncol = 1,
#                             rel_heights = c(1.5, 2, 2, 2)[o], scale = 1)
#
#
#     ggsave(pdf_name, plot = pg, width = 6, height = 4)
#     return(pg)
# }

#' Title
#'
#' @param p_list
#' @param qgr
#' @param rel_heights
#'
#' @return
#' @export
#'
#' @examples
track_assembly = function(p_list, qgr, rel_heights = rep(1, length(p_list))){
    p_list[[length(p_list)]] = p_list[[length(p_list)]] + labs(x = paste(as.character(seqnames(qgr)), "kbp"))
    p_list = sync_width(p_list)
    pg = cowplot::plot_grid(plotlist = p_list, ncol = 1,
                            rel_heights = rel_heights, scale = 1)
}
#' 
#' #' Title
#' #'
#' #' @param my_plots
#' #'
#' #' @return
#' #' @export
#' #' @import grid
#' #'
#' #' @examples
#' sync_width = function(my_plots){
#'     stopifnot(class(my_plots) == "list")
#'     stopifnot(all(sapply(my_plots, function(x)"ggplot" %in% class(x))))
#'     my_grobs = lapply(my_plots, function(x){
#'         ggplotGrob(x)
#'     })
#' 
#'     my_widths = lapply(my_grobs, function(gt){
#'         gt$widths
#'     })
#'     maxWidth = my_widths[[1]]
#'     if(length(my_widths) > 1){
#'         for(i in 2:length(my_widths)){
#'             maxWidth = grid::unit.pmax(maxWidth, my_widths[[i]])
#'         }
#'     }
#'     for(j in 1:length(my_grobs)){
#'         my_grobs[[j]]$widths = maxWidth
#'     }
#'     my_grobs
#' }

# library(GenomicRanges)
# library(data.table)
# library(ggplot2)
# library(seqsetvis)
# if(!exists("ref")){
#     ref = "~/gencode.v28.annotation.gtf.gz"
#     ref = rtracklayer::import.gff(ref, format =  "gtf", feature.type = "exon")
# }
#
# qgr = GRanges("chr6", IRanges(25e6, 25.002e6))
#
# track_ref(ref, qgr = qgr, flip_x = TRUE, debug = FALSE)
# track_ref(subset(ref, gene_type == "protein_coding"), qgr = qgr, flip_x = TRUE)
#
# runx1_bw = c(`MDA231 Runx1` = "/slipstream/galaxy/uploads/working/qc_framework/output_MK_MDA231_Runx/MDA231_Runx1_pooled/MDA231_Runx1_pooled_FE.bw",
#              `MDA231 Runx2` = "/slipstream/galaxy/uploads/working/qc_framework/output_MK_MDA231_Runx/MDA231_Runx2_pooled/MDA231_Runx2_pooled_FE.bw")
#
# track_chip(runx1_bw, qgr, flip_x = TRUE, debug = FALSE, nwin = 50)
