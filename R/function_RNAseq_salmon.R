# SALMON_PATH = "/slipstream/home/joeboyd/bin/salmon"
# GTF_TO_FASTA_PATH = "/slipstream/home/joeboyd/bin/gtf_to_fasta"
#
# HG38_v28_GTF_URL = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/gencode.v28.annotation.gtf.gz"
# HG38_SEQ_URL = "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz"
# HG38_SEQ_PATH = "/slipstream/galaxy/data/hg38/seq/hg38full.fa"
# fastq_paths = dir("~/R/JR_FLC_RNAseq_cross_promoter_capture/data_JR_RNAseq/fastqs/", full.names = TRUE)



#' salmon_index_transcriptome
#'
#' creates the salmon index for transcriptome created from specified assembly
#' fasta and gene reference gtf.
#'
#' @param seq_path path to local file or url to assembly fasta (only tested with hg38).
#' @param gtf_path path to local file or url to gtf reference (only tested with gencode hg38 v28 gtf).
#' @param gtf_to_fasta_path local path to tophat2 gtf_to_fasta utility.
#' @param salmon_path path to main salmon executable.
#' @param output_transcriptome file path for transcriptome output.
#' @param output_index file path for index for transcriptome.
#' @param cache_path path to cache directory, will be created if necessary.
#' @param use_qsub whether to qsub a script. if FALSE will use a direct system command
#' @param do_submit if FALSE, will not perform qsub.
#'
#' @return list containing path to index and job_id if relevant
#' @export
#' @import BiocFileCache
#' @examples
#' # as this may download several gigs of data, do not run
#' \dontrun{
#' salmon_index_transcriptome()
#' }
salmon_index_transcriptome = function(seq_path = HG38_SEQ_PATH,
                                      gtf_path = HG38_v28_GTF_URL,
                                      gtf_to_fasta_path = GTF_TO_FASTA_PATH,
                                      salmon_path = SALMON_PATH,
                                      output_transcriptome = NULL,
                                      output_index = NULL,
                                      cache_path = "~/.cache",
                                      use_qsub = TRUE,
                                      do_submit = TRUE){
    tmp_trans = paste0("tmp.", output_transcriptome)
    bfc = BiocFileCache::BiocFileCache(cache_path)
    if(is.null(output_transcriptome)){
        f = basename(gtf_path)
        f = sub(".gz$", "", f)
        f = sub(".gtf$", "", f)
        f = paste0(f, ".transcriptome.fa")
        output_transcriptome = f
    }
    if(is.null(output_index)){
        f = output_transcriptome
        f = sub(".fa$", "", f)
        f = sub(".fasta$", "", f)
        f = paste0(f, ".index")
        output_index = f
    }


    # if(is_url(seq_path)){
    #     seq_path = BiocFileCache::bfcrpath(bfc, seq_path)
    # }
    print(seq_path)
    seq_path = cache_gz(bfc, seq_path)
    print(gtf_path)
    gtf_path = cache_gz(bfc, gtf_path)

    # if(is_url(gtf_path)){
    #     gtf_path = BiocFileCache::bfcrpath(bfc, gtf_path)
    # }
    # if(grepl(".gz$", gtf_path)){ #gunzip if necessary
    #     rname = paste0(basename(gtf_path), ",unzip")
    #     if(!any(BiocFileCache::bfcinfo(bfc)$rname == rname)){ #dl if necessary
    #         gtf_raw = readLines(gzfile(gtf_path))
    #         gtf_file = BiocFileCache::bfcnew(bfc, rname = paste0(basename(gtf_path), ",unzip"))
    #         writeLines(gtf_raw, con = gtf_file)
    #         gtf_path = gtf_file
    #     }else{ # use existing resource
    #         gtf_path = BiocFileCache::bfcrpath(bfc, rnames = rname)
    #     }
    #
    # }
    tmp_trans = paste0("tmp.", basename(output_transcriptome))

    # log_file = paste0(output_transcriptome, ".log")
    cmd_gtf_fasta = "UTIL_VAR GTF_VAR SEQ_VAR TMP_VAR"
    cmd_gtf_fasta = sub("UTIL_VAR", gtf_to_fasta_path, cmd_gtf_fasta)
    cmd_gtf_fasta = sub("GTF_VAR", gtf_path, cmd_gtf_fasta)
    cmd_gtf_fasta = sub("SEQ_VAR", seq_path, cmd_gtf_fasta)
    cmd_gtf_fasta = sub("TMP_VAR", tmp_trans, cmd_gtf_fasta)

    cmd_clean_fasta = "cat TMP_VAR | awk '{if ($0 ~ \"^>\"){ sub(\">[0-9]+ \", \">\");} print $0 }' > OUT_VAR"
    cmd_clean_fasta = sub("TMP_VAR", tmp_trans, cmd_clean_fasta)
    cmd_clean_fasta = sub("OUT_VAR", output_transcriptome, cmd_clean_fasta)

    output_transcriptome_gz = paste0(output_transcriptome, ".gz")
    cmd_index = "UTIL_VAR index -t OUT_VAR -i INDEX_VAR"
    cmd_index = sub("UTIL_VAR", salmon_path, cmd_index)
    cmd_index = sub("OUT_VAR", output_transcriptome_gz, cmd_index)
    cmd_index = sub("INDEX_VAR", output_index, cmd_index)

    # for(f in c(tmp_trans, output_transcriptome, output_transcriptome_gz, output_index)){
    #     if(file.exists(f) | dir.exists(f))warning(f, " exists, not remaking.  delete to override.")
    # }
    if(any(file.exists(c(tmp_trans, output_transcriptome, output_transcriptome_gz, output_index)))){
        warning("Delete old results to override.")
    }

    all_cmds =  c(
        ifelse(file.exists(tmp_trans) |
                   file.exists(output_transcriptome) |
                   file.exists(output_transcriptome_gz), NA, cmd_gtf_fasta),
        ifelse(file.exists(output_transcriptome) |
                   file.exists(output_transcriptome_gz),
               NA,   cmd_clean_fasta),
        paste("rm -f", tmp_trans),
        ifelse(file.exists(output_transcriptome_gz), NA, paste("gzip", output_transcriptome)),
        ifelse(file.exists(output_index), NA, cmd_index)
    )
    all_cmds = all_cmds[!is.na(all_cmds)]
    if(use_qsub){
        bash_lines = c("#!/bin/bash",
                       "#$ -N salmon_index",
                       "#$ -cwd",
                       paste("#$ -e", cache_path),
                       paste("#$ -o", cache_path),
                       # paste0("#", cmd_gtf_fasta))
                       all_cmds)

        submit_file = "submit_salmon_index.sh"
        writeLines(bash_lines, submit_file)
        if(do_submit){
            sub_out = system(paste("qsub", submit_file), intern = TRUE)
            hjid = strsplit(sub_out, " ")[[1]]
            hjid = hjid[which(grepl("job", hjid))[1]+1]
        }else{
            hjid = NULL
        }
    }else{
        sapply(all_cmds, system)
        hjid = NULL
    }
    return(list(index_path = output_index, job_id = hjid))
}

# index_path = salmon_index_transcriptome(do_submit = FALSE)

#' salmon_quant_fastq_SE
#'
#' runs salmon quant using specified transcriptome index on every fastq supplied
#' for single-end reads only
#'
#' @param index_path path to index file
#' @param fastq_paths paths to fastq files (ungzipped?)
#' @param out_dirs directories to put results in per fastq
#'
#' @return list with paths to output directories and job ids
#' @export
#'
#' @examples
#' #nonsense without some toy data
#' \dontrun{
#' salmon_quant_fastq_SE(index_path, fastq_paths)
#' }
#'
salmon_quant_fastq_SE = function(index_path,
                                 fastq_paths,
                                 p = 8,
                                 out_path = getwd(),
                                 out_dirs = paste0("quant_", sub(".fastq", "", basename(fastq_paths))),
                                 do_submit = TRUE
){
    if(is.list(index_path)) index_path = index_path[[1]]
    stopifnot(length(out_path) == 1)
    dir.create(out_path, showWarnings = FALSE)
    out_dirs = paste0(out_path, "/", out_dirs)
    stopifnot(file.exists(index_path))
    stopifnot(all(file.exists(fastq_paths)))
    stopifnot(length(fastq_paths) == length(out_dirs))
    stopifnot(!any(duplicated(out_dirs)))
    cmd_quant = "salmon quant -i INDEX_VAR -l A -r FASTQ_VAR -p THREADS_VAR -o OUT_VAR --gcBias --seqBias"
    all_hjid = character()
    for(i in seq_along(fastq_paths)){
        cmd_this = cmd_quant
        cmd_this = sub("INDEX_VAR", index_path, cmd_this)
        cmd_this = sub("FASTQ_VAR", fastq_paths[i], cmd_this)
        cmd_this = sub("THREADS_VAR", p, cmd_this)
        cmd_this = sub("OUT_VAR", out_dirs[i], cmd_this)

        bash_lines = c("#!/bin/bash",
                       paste0("#$ -N salmon_quant_", i),
                       "#$ -cwd",
                       paste("#$ -e", out_dirs[i]),
                       paste("#$ -o", out_dirs[i]),
                       paste("#$ -pe threads", p),
                       cmd_this)
        submit_file = paste0("submit_salmon_quant_", i, ".sh")
        writeLines(bash_lines, submit_file)
        if(dir.exists(out_dirs[i])){
            warning(out_dirs[i], " output already exists, delete or submit manually")
            all_hjid = c(all_hjid, NULL)
        }else{
            if(do_submit){
                sub_out = system(paste("qsub", submit_file), intern = TRUE)
                hjid = strsplit(sub_out, " ")[[1]]
                hjid = hjid[which(grepl("job", hjid))[1]+1]
            }else{
                hjid = NULL
            }
            all_hjid = c(all_hjid, hjid)
        }

    }
    return(list(quant_results = out_dirs, job_ids = all_hjid))
}

#' Title
#'
#' @param sf_files
#' @param gtf_path
#' @param cache_path
#'
#' @return
#' @export
#'
#' @examples
salmon_tx_quant = function(sf_files, gtf_path = HG38_v28_GTF_URL, cache_path = "~/.cache"){
    bfc = BiocFileCache::BiocFileCache(cache_path)
    gtf_path = cache_gz(bfc, gtf_path)
    ref_gr = rtracklayer::import.gff(gtf_path, feature.type = "transcript", format = "gtf")
    tx2gene = data.frame(TXNAME = ref_gr$transcript_id, GENEID = ref_gr$gene_id)
    txi <- tximport::tximport(sf_files, type = "salmon", tx2gene = tx2gene)
    return(txi)
}

