
#' ssvR_plot_bam_mismatch
#'
#' @param bamfile
#' @param qgr
#'
#' @return
#' @export
#' @import Rsamtools
#' @import biovizBase
#' @import GenomeInfoDb
#' @rawNamespace import(data.table, except = c(shift, first, second, last))
#' @examples
#' library(GenomicRanges)
#' qgr = GRanges("chr21:34815000-34835000")
#' bam_file = file.path("/slipstream/galaxy/uploads/working/qc_framework/output_drugs_with_merged_inputs/MCF10A_ctrl_input_pooled/MCF10A_ctrl_input_pooled.bam")
#' ssvR_plot_bam_mismatch(bam_file, qgr)
ssvR_plot_bam_mismatch = function(bamfile, qgr,
                                  gen_ref = BSgenome.Hsapiens.UCSC.hg38::Hsapiens,
                                  min_count = 5, as_fraction = FALSE){
    qgr = GenomeInfoDb::dropSeqlevels(qgr, setdiff(GenomeInfoDb::seqlevels(qgr), as.character(GenomeInfoDb::seqnames(qgr))))
    GenomeInfoDb::seqlengths(qgr) = Rsamtools::scanBamHeader(bamfile)[[1]]$targets[GenomeInfoDb::seqlevels(qgr)]

    test <- biovizBase::pileupAsGRanges(bamfile, region = qgr)
    test.match <- biovizBase::pileupGRangesAsVariantTable(test, gen_ref)
    tdt = data.table::as.data.table(test.match)
    tdt[, .N, start][order(N)]
    term_start = tdt[!start %in% (start+1)]$start - 1
    term_end = tdt[!end %in% (end-1)]$end + 1
    term_dt = data.table::data.table(seqnames = tdt$seqnames[1],
                                     start = c(term_start, term_end),
                                     end = c(term_start, term_end),
                                     width = 1, strand = "+", ref = "C", read = "C",
                                     count = 0, depth = 0,
                                     bam = tdt$bam[1], match = TRUE)
    tdt2 = rbind(tdt,
                 term_dt
    )

    xmin = start(qgr)
    xmax = end(qgr)

    tdt2$read = factor(tdt2$read, levels = names(biovizBase::getBioColor("DNA_BASES_N")))


    rdt = tdt2[match == FALSE]
    rdt[, .N, by = .(start)]
    shift_dt = data.table::dcast(rdt, "seqnames+start~read", value.var = "count", fill = 0)
    shift_dt[, `A` := `C` + `G` + `T`]
    shift_dt[, `C` := `G` + `T`]
    shift_dt[, `G` := `T`]
    shift_dt[, `T` := 0]
    shift_dt = data.table::melt(shift_dt, id.vars = c("seqnames", "start"), variable.name = "read", value.name = "shift")
    rdt = merge(rdt, shift_dt)

    if(as_fraction){
        p = ggplot(rdt[count > min_count],
                   aes(xmin = start, xmax = start+1,
                       ymin = shift/depth, ymax = (shift+count)/depth,
                       fill = read, color = read)) +
            geom_rect()
    }else{
        p = ggplot() +
            geom_ribbon(data = tdt2[match == TRUE],
                        aes(x = start,
                            ymin = 0,
                            ymax = depth), fill = "gray40") +
            coord_cartesian(xlim = c(xmin, xmax)) +
            labs(x = "bp", y = "coverage", color = "mismatch")
        p = p + geom_rect(data = rdt[count > min_count],
                          aes(xmin = start, xmax = start+1,
                              ymin = shift, ymax = shift+count,
                              fill = read, color = read),
                          size = 1.2)

    }
    p + scale_fill_manual(values = biovizBase::getBioColor("DNA_BASES_N"), drop = FALSE)+
        scale_color_manual(values = biovizBase::getBioColor("DNA_BASES_N")) +
        guides(color = "none") +
        coord_cartesian(xlim = c(xmin, xmax)) +
        theme_classic() +
        theme(panel.background = element_rect(fill= "black"))
}

# library(ggplot2)
# bamfile = "~/R/RNAseq_BRCA/alignment/MCF10A_RNA-Seq_R1_GSM1944515.fastq_STARAligned.sortedByCoord.out.bam"
# qgr = GenomicRanges::GRanges("chr21", IRanges::IRanges(34789000, 34791500))
#
# ssvR_plot_bam_mismatch(bamfile, qgr)
# ssvR_plot_bam_mismatch(bamfile, resize(qgr, 20000, fix = "center"))


