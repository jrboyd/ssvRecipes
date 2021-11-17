
#' gr_distance
#'
#' Find items in b within a certain distance (max_dist) all items in a.
#'
#' Items in a will be repeated for every item in b they overlap with.
#'
#' a is the left side of output and the regions that extended when applying max_dist.
#'
#' @param a The left side of output, regions that get extended when applying max_dist.
#' @param b The right side of output
#' @param max_dist Maximum distance allowed from a to b.
#'
#' @return a GRanges of b appended to regions they overlap with from a.
#' @export
#' @import GenomicRanges
#' @import rtracklayer
#' @import S4Vectors
#'
#' @examples
#'
#' max_distance = 5e4
#' np_file = "/slipstream/galaxy/uploads/working/qc_framework/output/MCF7_H3K4ME3_R1/MCF7_H3K4ME3_R1_peaks.narrowPeak"
#' np_gr = seqsetvis::easyLoad_narrowPeak(np_file)[[1]]
#'
#' ref_gr = rtracklayer::import.gff("~/gencode.v36.annotation.gtf", feature.type = "gene")
#' pr_gr = GenomicRanges::promoters(ref_gr, 1)
#'
#' olap_gr = gr_distance_overlap(pr_gr, np_gr, max_distance)
#' #not the impact of reversing the inputs
#' olap_gr.flip = gr_distance_overlap(np_gr, pr_gr, max_distance)
#'
#' wd = "/slipstream/home/sophiakogut/IK_peaks"
#' np_file2 = file.path(wd, "IK_common.narrowPeak")
#' mir_file = file.path(wd, "DE_mirna_primary.txt")
#'
#' # This ain't a narrowPeak file!
#' np_gr2 = seqsetvis::easyLoad_bed(np_file2)[[1]]
#' mir_gr = seqsetvis::easyLoad_bed(mir_file)[[1]]
#' mir_olaps = gr_distance_overlap(mir_gr, np_gr2, 5e4)
#' mir_olaps
#'
#' # writing the output
#' df = as.data.frame(mir_olaps)
#' df = df[, colnames(df)[c(1:3, 6, 7, 5, seq(8, ncol(df)))]]
#' # this file is very similar to a bed file and will work with bedtools if you remove the first line with column names
#' # the first 6 columns are standard bed format
#' write.table(df, "mir_nearby_peaks.txt", sep = "\t", quote = FALSE, row.names = FALSE)
#'
#' # load it back into R like this
#' loaded_gr = GRanges(read.table("mir_nearby_peaks.txt", header = TRUE))
gr_distance_overlap = function(a, b, max_dist){
    stopifnot(is(a, "GRanges"))
    stopifnot(is(b, "GRanges"))
    stopifnot(is(max_dist, "numeric"))
    
    query_gr = GenomicRanges::resize(a, 2*max_dist + width(a), fix = "center")
    olaps = as.data.frame(GenomicRanges::findOverlaps(b, query_gr))
    
    dist_vals = GenomicRanges::distance(b[olaps$queryHits], a[olaps$subjectHits])
    append_df = as.data.frame(GenomicRanges::mcols(b))
    colnames(append_df)[colnames(append_df) == "name"] = "peak_name"
    append_df$peak_seqnames = as.character(seqnames(b))
    append_df$peak_start = start(b)
    append_df$peak_end = end(b)
    append_df$peak_strand = as.character(strand(b))
    
    out_gr = a[olaps$subjectHits]
    
    out_gr$distance_to_peak = dist_vals
    
    GenomicRanges::mcols(out_gr) = cbind(GenomicRanges::mcols(out_gr), S4Vectors::DataFrame(append_df[olaps$queryHits,]))
    out_gr
}


