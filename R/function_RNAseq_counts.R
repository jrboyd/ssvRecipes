# fastq_paths = dir("~/R/JR_FLC_RNAseq_cross_promoter_capture/data_JR_RNAseq/fastqs/", full.names = TRUE)
# bams = dir("~/R/JR_FLC_RNAseq_cross_promoter_capture/data_JR_RNAseq/alignments/", full.names = TRUE, pattern = ".bam$")

#' Title
#'
#' @param fastq_paths paths to fastq files
#' @param index_path path to star index
#' @param gtf_path path (local or url) to gtf
#' @param star_path path to star executable
#' @param cache_path cache location
#' @param p number of threads to use
#' @param hold_jids job ids to hold for, default NA is none.
#' @param out_path directory to output to
#' @param out_dirs custom output directory per fastq_paths. default is basename
#' of fastq_paths
#' @param do_submit if FALSE, qsub is skipped but submit scripts remain.
#' Default is TRUE.
#'
#' @return list length 2 of output dirs and job ids
#' @export
#' @import BiocFileCache
#'
#' @examples
star_align_fastq_SE = function(fastq_paths,
                               index_path = HG38_STAR_INDEX,
                               gtf_path = HG38_v28_GTF_URL,
                               star_path = STAR_PATH,
                               cache_path = "~/.cache",
                               p = 8,
                               hold_jids = NA,
                               out_path = getwd(),
                               out_dirs = paste0(sub(".fastq", "", basename(fastq_paths)), ".STAR."),
                               do_submit = TRUE
){
    bfc = BiocFileCache::BiocFileCache(cache_path)
    if(is.list(index_path)) index_path = index_path[[1]]
    stopifnot(length(out_path) == 1)
    dir.create(out_path, showWarnings = FALSE)
    out_dirs = paste0(out_path, "/", out_dirs)
    stopifnot(file.exists(index_path))
    stopifnot(all(file.exists(fastq_paths)))
    stopifnot(length(fastq_paths) == length(out_dirs))
    stopifnot(!any(duplicated(out_dirs)))
    stopifnot(length(hold_jids) == 1 | length(hold_jids) == length(fastq_paths))
    if(length(hold_jids) == 1){
        hold_jids = rep(hold_jids, length(fastq_paths))
    }
    cmd_align = "STAR_VAR \
--runThreadN 8 \
--genomeDir INDEX_VAR \
--readFilesIn FASTQ_VAR \
--outFileNamePrefix OUT_VAR \
--outSAMtype BAM SortedByCoordinate \
--outMultimapperOrder Random \
--outFilterType BySJout \
--outFilterMultimapNmax 20 \
--alignSJoverhangMin THREADS_VAR \
--alignSJDBoverhangMin 1 \
--outFilterMismatchNmax 999 \
--outFilterMismatchNoverReadLmax 0.04 \
--alignIntronMin 20 \
--alignIntronMax 1000000 \
--alignMatesGapMax 1000000 \
--sjdbGTFfile GTF_VAR"
    cmd_align = gsub("\n", "", cmd_align)
    gtf_path = cache_gz(bfc, gtf_path)
    cmd_align = sub("GTF_VAR", gtf_path, cmd_align)
    cmd_align = sub("STAR_VAR", star_path, cmd_align)
    cmd_align = sub("INDEX_VAR", index_path, cmd_align)

    all_hjid = character()
    for(i in seq_along(fastq_paths)){
        hjid = hold_jids[i]
        cmd_this = cmd_align
        cmd_this = sub("FASTQ_VAR", fastq_paths[i], cmd_this)
        cmd_this = sub("THREADS_VAR", p, cmd_this)
        cmd_this = sub("OUT_VAR", out_dirs[i], cmd_this)

        bamout = paste0(out_dirs[i], "Aligned.sortedByCoord.out.bam")
        cmd_index = paste("samtools index", bamout)

        bash_lines = c(
            "#!/bin/bash",
            paste0("#$ -N STAR_align_", i),
            "#$ -cwd",
            paste("#$ -e", out_path),
            paste("#$ -o", out_path),
            paste("#$ -pe threads", p),
            ifelse(file.exists(bamout), "#bam exists, skip alignment", cmd_this),
            ifelse(file.exists(paste0(bamout, ".bai")), "#bam.bai exists, skip index", cmd_index),
            "echo done!"
        )
        submit_file = paste0("submit_STAR_align_", i, ".sh")
        writeLines(bash_lines, submit_file)

        # if(dir.exists(out_dirs[i])){
        #     warning(out_dirs[i], " output already exists, delete or submit manually")
        #     all_hjid = c(all_hjid, NULL)
        # }else{
            if(do_submit){
                if(!is.na(hjid)){
                    sub_out = system(paste("qsub -hold_jid", hjid, submit_file), intern = TRUE)
                }else{
                    sub_out = system(paste("qsub", submit_file), intern = TRUE)
                }

                hjid = strsplit(sub_out, " ")[[1]]
                hjid = hjid[which(grepl("job", hjid))[1]+1]
            }else{
                hjid = NULL
            }
            all_hjid = c(all_hjid, hjid)
        # }

    }
    return(list(result_paths = out_dirs, job_ids = all_hjid))
}

bam_plot_dt = function(bam_paths, qgr, flip_strand = TRUE){
    pos_bam = ssvFetchBam(file_paths = bam_paths, qgr = qgr, target_strand = "+", return_data.table = TRUE, fragLens = NA)
    neg_bam = ssvFetchBam(file_paths = bam_paths, qgr = qgr, target_strand = "-", return_data.table = TRUE, fragLens = NA)

    if(flip_strand){
        #strands seem flipped
        pos_bam$strand = "-"
        neg_bam$strand = "+"
    }

    bam_dt = rbind(pos_bam, neg_bam)
    bam_dt
}

bam_plot = function(bam_paths, qgr){
    bam_dt = bam_plot_dt(bam_paths, qgr)
    bam_dt[, sample_label := sub(".fastq_STARAligned.sortedByCoord.out.bam", "", sample)]
    ggplot(bam_dt, aes(x = x, y = y, color = strand)) + geom_path() + facet_wrap("sample_label", ncol = 3)
}

#' Title
#'
#' @return
#' @export
#' @import GenomicRanges
#' @import seqsetvis
#' @import Rsamtools
#' @import GenomicFeatures
#' @import GenomicAlignments
#' @import BiocParallel
#' @examples
counts_from_bams = function(bam_paths,
                            p = 8,
                            gtf_path = HG38_v28_GTF_URL,
                            counts_tag = "counts",
                            cache_path = "~/.cache"){
    # library(GenomicRanges)
    # library(seqsetvis)
    # library("Rsamtools")
    # library(GenomicFeatures)
    # library(GenomicAlignments)
    # library("BiocParallel")
    # library(ssvRecipes)
    # library(DESeq2)
    # library(data.table)
    bfc = BiocFileCache::BiocFileCache(cache_path)
    stopifnot(all(file.exists(bam_paths)))
    bamfiles <- Rsamtools::BamFileList(bam_paths, yieldSize=2000000)
    gtf_path = cache_gz(bfc, gtf_path)

    ### TODO extract to functiont that conditionally applies function and stores result in cache or returns cached
    txdb = cache_FUN(bfc, gtf_path, "txdb", function(gtf_path){
        GenomicFeatures::makeTxDbFromGFF(gtf_path, format="gtf")
    }, saveFUN = AnnotationDbi::saveDb, loadFUN = AnnotationDbi::loadDb)

    ebg = cache_FUN(bfc, gtf_path, "ebg", function(txdb){
        GenomicFeatures::exonsBy(txdb, by="gene")
    }, input_obj = txdb)



    se = cache_FUN(bfc, gtf_path, counts_tag, function(in_list){
        BiocParallel::register(BiocParallel::MulticoreParam(p))
        ebg = in_list$ebg
        bamfiles = in_list$bamfiles
        #single-end reads
        GenomicAlignments::summarizeOverlaps(features=ebg,
                                             reads=bamfiles,
                                             mode="Union",
                                             singleEnd=TRUE,
                                             ignore.strand=FALSE)
    }, input_obj = list(ebg = ebg, bamfiles = bamfiles))
    return(se)
}

