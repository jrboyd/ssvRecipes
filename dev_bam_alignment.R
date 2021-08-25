#from ssvRecipes::ssvR_plot_bam_mismatch
library(GenomicRanges)
bamfile = "~/R/RNAseq_AT1_sorted/alignment/MCF10A_AT1_R1_Bulk_CGATGT_AHCYCJBCXY_S1_L001_R1_001.fastq_STARAligned.sortedByCoord.out.bam"
qgr = GRanges("chr21:34886836-34887135")
qgr = GRanges("chr21:34788191-34794827")

undebug(ssvRecipes::ssvR_plot_bam_mismatch)
ssvRecipes::ssvR_plot_bam_mismatch(bamfile, qgr, as_fraction = TRUE)
ssvRecipes::ssvR_plot_bam_mismatch(bamfile, qgr, as_fraction = FALSE)

qgr = GenomeInfoDb::dropSeqlevels(qgr, setdiff(GenomeInfoDb::seqlevels(qgr), as.character(GenomeInfoDb::seqnames(qgr))))
GenomeInfoDb::seqlengths(qgr) = Rsamtools::scanBamHeader(bamfile)[[1]]$targets[GenomeInfoDb::seqlevels(qgr)]

test <- biovizBase::pileupAsGRanges(bamfile, region = qgr)
test.match <- biovizBase::pileupGRangesAsVariantTable(test, gen_ref)
tdt = data.table::as.data.table(test.match)
