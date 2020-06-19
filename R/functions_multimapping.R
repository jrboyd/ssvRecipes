#' write a fasta file of specified read_size from reference sequence for alignment
#' @param gr region to generate reads for
#' @param fasta_file path to output fasta file
#' @param mc number of cores to use. default is half of max.
#' @param genome BSgenome to use
#' @param read_size size of reads to make
writeFastaFromGR = function(gr, fasta_file, mc = parallel::detectCores() - 4,
                            genome = BSgenome.Hsapiens.UCSC.hg38::Hsapiens, read_size = 75, wrap = 100){
    #create fasta of 75 bp in silico reads of Hist2 locus
    mySeq = Biostrings::getSeq(genome, gr)
    mySeq = as.character(mySeq)
    mySeq = sub("^N+", "", mySeq)
    fasta_file = writeFastaFromSeq(mySeq, fasta_file = fasta_file, mc = mc, read_size = read_size, wrap = wrap)
    return(fasta_file)
}

#' write a fasta file of specified read_size from reference sequence for alignment
#' @param seq the seq
#' @param fasta_file path to output fasta file
#' @param mc number of cores to use. default is half of max.
#' @param genome BSgenome to use
#' @param read_size size of reads to make
writeFastaFromSeq = function(mySeq, fasta_file, mc = parallel::detectCores() - 4,
                             read_size = 75, wrap = 100){
    
    outf = fasta_file
    suppressWarnings({
        file.remove(outf)
    })
    mySplit = strsplit(mySeq, "")[[1]]
    
    MAX = nchar(mySeq) - read_size + 1
    
    hidden = parallel::mclapply(seq_len(MAX), mc.cores = mc, function(itr){
        rn = paste0(">read_", itr)
        seq = paste(x = mySplit[itr:(itr + read_size - 1)], collapse = "")
        rem = nchar(seq) %% wrap
        wrap_count = floor(nchar(seq) / wrap)
        
        tmp = sapply(seq_len(wrap_count), function(wi){
            substr(seq, (wi-1)*wrap + 1, (wi-1)*wrap + wrap)
            
        })
        if(rem > 0){
            tmp = c(tmp, substr(seq, nchar(seq)-rem+1, nchar(seq)))
        }
        write.table(rbind(rn, cbind(tmp)), file = outf,
                    append = T, row.names = F, col.names = F, quote = F)
    })
    return(fasta_file)
}

#' calls STAR on fasta
#' @param fasta_file fasta or fastq file to align
#' @param prefix prefix for STAR output
#' @param cleanup suffixes of files to remove.  by default only sam is kept.
#' @param star_path path to STAR executable
#' @param index_path path to genome indexes for STAR
alignFasta = function(fasta_file, 
                      prefix = sub("fasta", "", fasta_file),
                      cleanup = c("Log.out", "Log.progress.out", "SJ.out.tab"),
                      star_path = "/slipstream/galaxy/production/galaxy-dist/tools/star/STAR", 
                      index_path = "/slipstream/galaxy/data/hg38/star_index"){
    
    if(prefix == fasta_file) stop("can't determine prefix, please set manually.")
    cmd = paste(star_path, 
                "--genomeDir", index_path, 
                "--readFilesIn", fasta_file, 
                "--outFileNamePrefix", prefix, 
                "--alignIntronMax=0", 
                "--runThreadN=8")
    system(cmd)
    
    to_remove = paste0(prefix, cleanup)
    if(any(file.exists(to_remove))){
        file.remove(to_remove[file.exists(to_remove)])    
    }
    
}


writeBB = function(gr, prefix, chrSizes_file = "~/hg38_chrsizes.canon.txt"){
    bed_file = paste0(prefix, ".bed")
    sort_file = paste0(prefix, ".sorted.bed")
    bb_file = paste0(prefix, ".bb")
    
    rtracklayer::export.bed(gr, bed_file)
    
    bedSort = paste("bedSort", bed_file, sort_file)
    bedToBigBed = paste("bedToBigBed", sort_file, chrSizes_file, bb_file)
    
    system(bedSort)
    system(bedToBigBed)
    
    file.remove(bed_file)
    file.remove(sort_file)
    invisible(bb_file)
}
