# fastq_paths = dir("~/R/JR_FLC_RNAseq_cross_promoter_capture/data_JR_RNAseq/fastqs/", full.names = TRUE)
# bams = dir("~/R/JR_FLC_RNAseq_cross_promoter_capture/data_JR_RNAseq/alignments/", full.names = TRUE, pattern = ".bam$")

#' Title
#' @param FASTQ_VAR readFilesIn arg to supply STAR
#' @param fastq_paths paths to fastq files
#' @param index_path path to star index
#' @param gtf_path path (local or url) to gtf
#' @param star_path path to star executable
#' @param cache_path cache location
#' @param n_cores number of threads to use
#' @param hold_jids job ids to hold for, default NA is none.
#' @param out_path directory to output to
#' @param output_prefix custom output directory per fastq_paths. default is basename
#' of fastq_paths
#' @param do_submit if FALSE, qsub is skipped but submit scripts remain.
#' Default is TRUE.
#'
#' @return list length 2 of output dirs and job ids
#' @export
#' @import BiocFileCache
#'
#' @examples
star_align_fastq_core = function(FASTQ_VAR, 
                                 fastq_paths,
                                 index_path = HG38_STAR_INDEX,
                                 gtf_path = HG38_v28_GTF_URL,
                                 star_path = STAR_PATH,
                                 cache_path = "~/.cache",
                                 n_cores = 8,
                                 hold_jids = NA,
                                 out_path = file.path(getwd(), "alignment"),
                                 output_prefix = NULL,
                                 do_submit = TRUE){
    bfc = BiocFileCache::BiocFileCache(cache_path)
    if(is.list(index_path)) index_path = index_path[[1]]
    stopifnot(length(out_path) == 1)
    dir.create(out_path, showWarnings = FALSE, recursive = TRUE)
    if(is.null(output_prefix)){
        output_prefix = sub("\\.fastq.+", ".", basename(fastq_paths))
    }
    output_prefix = file.path(out_path, output_prefix)
    stopifnot(file.exists(index_path))
    stopifnot(all(file.exists(fastq_paths)))
    stopifnot(length(FASTQ_VAR) == length(output_prefix))
    stopifnot(!any(duplicated(output_prefix)))
    stopifnot(length(hold_jids) == 1 | length(hold_jids) == length(FASTQ_VAR))
    if(length(hold_jids) == 1){
        hold_jids = rep(hold_jids, length(FASTQ_VAR))
    }
    
    is_fastq = all(grepl("\\.fastq$", fastq_paths))
    is_fastqgz = all(grepl("\\.fastq.gz$", fastq_paths))
    stopifnot(is_fastq | is_fastqgz)
    
    log_dir = file.path(normalizePath(out_path), "alignment_logs")
    dir.create(log_dir, showWarnings = FALSE)
    submit_dir = file.path(normalizePath(out_path), "alignment_scripts")
    dir.create(submit_dir, showWarnings = FALSE)
    
    
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
    if(is_fastqgz){
        cmd_align = paste(cmd_align, "\
                          --readFilesCommand zcat")
    }
    cmd_align = gsub("\n", "", cmd_align)
    gtf_path = cache_gz(bfc, gtf_path)
    cmd_align = sub("GTF_VAR", gtf_path, cmd_align)
    cmd_align = sub("STAR_VAR", star_path, cmd_align)
    cmd_align = sub("INDEX_VAR", index_path, cmd_align)
    
    all_hjid = character()
    for(i in seq_along(FASTQ_VAR)){
        hjid = hold_jids[i]
        cmd_this = cmd_align
        cmd_this = sub("FASTQ_VAR", FASTQ_VAR[i], cmd_this)
        cmd_this = sub("THREADS_VAR", n_cores, cmd_this)
        cmd_this = sub("OUT_VAR", output_prefix[i], cmd_this)
        
        bamout = paste0(output_prefix[i], "Aligned.sortedByCoord.out.bam")
        cmd_index = paste("samtools index", bamout)
        
        bash_lines = c(
            "#!/bin/bash",
            paste0("#$ -N STAR_align_", i),
            "#$ -cwd",
            paste("#$ -e", log_dir),
            paste("#$ -o", log_dir),
            paste("#$ -pe threads", n_cores),
            ifelse(file.exists(bamout), "#bam exists, skip alignment", cmd_this),
            ifelse(file.exists(paste0(bamout, ".bai")), "#bam.bai exists, skip index", cmd_index),
            "echo done!"
        )
        j = 1
        submit_file = file.path(submit_dir, paste0("submit_STAR_align.", basename(output_prefix[i]), ".", j, ".sh"))
        
        while(file.exists(submit_file)){
            j = j + 1
            submit_file = file.path(submit_dir, paste0("submit_STAR_align.", basename(output_prefix[i]), ".", j, ".sh"))
        }
        writeLines(bash_lines, submit_file)
        
        # if(dir.exists(output_prefix[i])){
        #     warning(output_prefix[i], " output already exists, delete or submit manually")
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
    return(list(result_paths = output_prefix, job_ids = all_hjid))
    
}

# fastq_paths = dir("~/R/JR_FLC_RNAseq_cross_promoter_capture/data_JR_RNAseq/fastqs/", full.names = TRUE)
# bams = dir("~/R/JR_FLC_RNAseq_cross_promoter_capture/data_JR_RNAseq/alignments/", full.names = TRUE, pattern = ".bam$")

#' Title
#'
#' @param fastq_paths paths to fastq files
#' @param index_path path to star index
#' @param gtf_path path (local or url) to gtf
#' @param star_path path to star executable
#' @param cache_path cache location
#' @param n_cores number of threads to use
#' @param hold_jids job ids to hold for, default NA is none.
#' @param out_path directory to output to
#' @param output_prefix custom output directory per fastq_paths. default is basename
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
                               n_cores = 8,
                               hold_jids = NA,
                               out_path = file.path(getwd(), "alignment"),
                               output_prefix = NULL,
                               do_submit = TRUE
){
    FASTQ_VAR = fastq_paths
    star_align_fastq_core(FASTQ_VAR = FASTQ_VAR, 
                          fastq_paths = fastq_paths,
                          index_path = index_path, 
                          gtf_path = gtf_path, 
                          star_path = star_path, 
                          cache_path = cache_path, 
                          n_cores = n_cores, 
                          hold_jids = hold_jids, 
                          out_path = out_path, 
                          output_prefix = output_prefix, 
                          do_submit = do_submit)
}

#' Title
#'
#' @param r1_fastq_paths paths to r1 files, mandatory
#' @param r2_fastq_paths paths to r2 files, will try to derive by calling r1_to_r2_FUN(r1_fastq_paths)
#' @param r1_to_r2_FUN optional function if r2 files no supplied
#' @param gtf_path path (local or url) to gtf
#' @param star_path path to star executable
#' @param cache_path cache location
#' @param n_cores number of threads to use
#' @param hold_jids job ids to hold for, default NA is none.
#' @param out_path directory to output to
#' @param output_prefix custom output directory per fastq_paths. default is basename
#' of fastq_paths
#' @param do_submit if FALSE, qsub is skipped but submit scripts remain.
#' Default is TRUE.
#'
#' @return list length 2 of output dirs and job ids
#' @export
#' @import BiocFileCache
#'
#' @examples
star_align_fastq_PE = function(r1_fastq_paths,
                               r2_fastq_paths = NULL,
                               pair_key = "P",
                               r1_to_r2_FUN = function(r1){sub(paste0(pair_key, "1.fastq"), paste0(pair_key, "2.fastq"), r1)},
                               index_path = HG38_STAR_INDEX,
                               gtf_path = HG38_v28_GTF_URL,
                               star_path = STAR_PATH,
                               cache_path = "~/.cache",
                               n_cores = 8,
                               hold_jids = NA,
                               out_path = file.path(getwd(), "alignment"),
                               output_prefix = paste0(sub(paste0("_", pair_key, "1.fastq"), "", basename(r1_fastq_paths)), ".STAR."),
                               do_submit = TRUE
){
    if(is.null(r2_fastq_paths)){
        r2_fastq_paths = r1_to_r2_FUN(r1_fastq_paths)
    }
    FASTQ_VAR = paste(r1_fastq_paths[i], r2_fastq_paths[i])
    star_align_fastq_core(FASTQ_VAR = FASTQ_VAR, 
                          fastq_paths = c(r1_fastq_paths, r2_fastq_paths),
                          index_path = index_path, 
                          gtf_path = gtf_path, 
                          star_path = star_path, 
                          cache_path = cache_path, 
                          n_cores = n_cores, 
                          hold_jids = hold_jids, 
                          out_path = out_path, 
                          output_prefix = output_prefix, 
                          do_submit = do_submit)
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
                            n_cores = 8,
                            gtf_path = HG38_v28_GTF_URL,
                            counts_tag = "counts",
                            cache_path = "~/.cache",
                            flip_strand = FALSE,
                            singleEnd=TRUE,
                            inter.feature = FALSE,
                            ignore.strand=FALSE,
                            mode="IntersectionStrict"){
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
    rname = digest::digest(list(bam_paths, gtf_path, flip_strand, singleEnd, inter.feature, ignore.strand, mode))
    
    if(flip_strand){
        preprocess.reads = GenomicAlignments::invertStrand
    }else{
        preprocess.reads = NULL    
    }
    
    # se = cache_FUN(bfc, gtf_path, counts_tag, function(in_list){
    se = bfcif(bfc, rname, function(){
        
        
        BiocParallel::register(BiocParallel::MulticoreParam(n_cores))
        # ebg = in_list$ebg
        # bamfiles = in_list$bamfiles
        #single-end reads
        GenomicAlignments::summarizeOverlaps(features=ebg,
                                             reads=bamfiles,
                                             mode=mode,
                                             singleEnd=singleEnd,
                                             inter.feature = inter.feature,
                                             ignore.strand=ignore.strand,
                                             preprocess.reads=preprocess.reads)
    })
    # }, input_obj = list(ebg = ebg, bamfiles = bamfiles))
    return(se)
}

#' rnaseq_asses_strandedness
#'
#' @param bam_files path to bam files
#' @param gtf_path path to gtf with exon info
#' @param max_bams only this many bams are used
#'
#' @return a grob of plots
#' @export
#'
#' @examples
rnaseq_asses_strandedness = function(bam_files, gtf_path, max_bams = 3){
    theme_set(cowplot::theme_cowplot())
    if(length(bam_files) > max_bams){
        bam_files = sample(bam_files, max_bams)
    }
    
    ex_gr = rtracklayer::import.gff(gtf_path, feature.type = "exon", format = 'gtf')
    qgr = sample(subset(ex_gr, gene_type == "protein_coding" & width(ex_gr) > 300), 500)
    names(qgr) = paste0("id_", seq_along(qgr))
    no_flip_dt = seqsetvis::ssvFetchBam(bam_files, qgr = qgr, 
                                        win_method = "summary", win_size = 1, 
                                        summary_FUN = function(x,w)sum(x), fragLens = 1, 
                                        n_cores = 20, return_data.table = TRUE)
    yes_flip_dt = seqsetvis::ssvFetchBam(bam_files, qgr = qgr, 
                                         win_method = "summary", win_size = 1, 
                                         summary_FUN = function(x,w)sum(x), fragLens = 1, 
                                         n_cores = 20, return_data.table = TRUE, 
                                         flip_strand = TRUE)
    
    dt = cbind(no_flip_dt[, .(id, no_flip = y)], yes_flip_dt[, .(yes_flip = y)])
    lim = log10(range(dt$no_flip, dt$yes_flip) +1)
    p1 = ggplot(dt, aes(x = log10(no_flip+1), y = log10(yes_flip+1))) + geom_point() +
        annotate("line", x = lim, y = lim, color = "red") +
        labs(x = "not flipped (log10)", y = "flipped (log10)", title = "flipped vs not flipped")
    
    qid_plus = dt[id %in% names(subset(qgr, strand == "+")), .(total = sum(no_flip+yes_flip)), .(id)][order(total, decreasing = TRUE)][1:5]$id
    qid_neg = dt[id %in% names(subset(qgr, strand == "-")), .(total = sum(no_flip+yes_flip)), .(id)][order(total, decreasing = TRUE)][1:5]$id
    qid = c(qid_plus, qid_neg)
    strand(qgr) = "*"
    no_flip_dt = seqsetvis::ssvFetchBam(data.table(file = bam_files, sample = substr(basename(bam_files), 1, 8)), 
                                        qgr = resize(qgr[qid], width(qgr)*2, fix = "center"), target_strand = "both",
                                        win_method = "summary", win_size = 100, 
                                        summary_FUN = function(x,w)sum(x), fragLens = NA, 
                                        n_cores = 20,
                                        return_data.table = TRUE)
    
    p2 = ggplot(no_flip_dt[id %in% qid_plus], aes(x = x, y = y, color = strand)) + 
        geom_path() + facet_grid(id~sample, scales = "free_y") +
        scale_color_manual(values = c("+" = "red", "-" = "blue")) +
        labs(title = "not flippped (+)strand", y = "pileup", x = "exon") +
        theme(legend.position = "bottom", 
              strip.text = element_text(size = 6),
              axis.text = element_text(size= 6))
    p3 = ggplot(no_flip_dt[id %in% qid_neg], aes(x = x, y = y, color = strand)) + 
        geom_path() + facet_grid(id~sample, scales = "free_y") +
        scale_color_manual(values = c("+" = "red", "-" = "blue")) +
        labs(title = "not flippped (-)strand", y = "pileup", x = "exon") +
        theme(legend.position = "bottom", 
              strip.text = element_text(size = 6), 
              axis.text = element_text(size= 6))
    plot(cowplot::plot_grid(p1, p2, p3, nrow = 1))
    invisible(list(scatter = p1, tracks_plus = p2, tracks_negative = p3))
}
