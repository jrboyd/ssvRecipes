# assess_lib = function(dt_stats,
#                       miss_cutoff = 1.2,
#                       hit_cutoff = 1.2,
#                       split_cutoff = .7){
#     dt_stats = dt_stats[grepl("N_", gene_id)]
#     tdt_stats = melt(dt_stats, variable.name = "lib", value.name = "count", id.vars = "gene_id")
#     
#     un_count = tdt_stats[gene_id == "N_noFeature" & lib == "unstranded"]$count
#     fir_count = tdt_stats[gene_id == "N_noFeature" & lib == "first"]$count
#     sec_count = tdt_stats[gene_id == "N_noFeature" & lib == "second"]$count
#     
#     fir_ratio = fir_count / un_count
#     sec_ratio = sec_count / un_count
#     data.table(un_count, fir_count, sec_count, fir_ratio, sec_ratio)
# }


guess_lib = function(dt_stats,
                     miss_cutoff = 1.2,
                     hit_cutoff = 1.2,
                     split_cutoff = 1.5){
    dt_stats = dt_stats[grepl("N_", gene_id)]
    tdt_stats = melt(dt_stats, variable.name = "lib", value.name = "count", id.vars = "gene_id")
    
    un_count = tdt_stats[gene_id == "N_noFeature" & lib == "unstranded"]$count
    fir_count = tdt_stats[gene_id == "N_noFeature" & lib == "first"]$count
    sec_count = tdt_stats[gene_id == "N_noFeature" & lib == "second"]$count
    
    fir_ratio = fir_count / un_count
    sec_ratio = sec_count / un_count
    
    lib = "unrecognized"
    
    if(sec_ratio > miss_cutoff & fir_ratio < hit_cutoff){
        lib = "first"
    }else if(fir_ratio > miss_cutoff & sec_ratio < hit_cutoff){
        lib = "second"
    }else if(fir_ratio > split_cutoff & sec_ratio > split_cutoff){#this could be wrong
        lib = "unstranded"
    }
    lib
}

#' guess_lib_from_file
#'
#' @param f 
#' @param show_plots 
#' @param miss_cutoff 
#' @param hit_cutoff 
#' @param split_cutoff 
#'
#' @return
#' @export
#' @import cowplot ggplot2 
#' @rawNamespace import(data.table, except = c(shift, first, second, last))
#' @examples
#' files = dir("inst/extdata", pattern = "^test_[0-9].ReadsPerGene.out.tab", full = TRUE)
#' #unstranded
#' guess_lib_from_file(files[1], show_plots = TRUE)
#' guess_lib_from_file(files[2], show_plots = TRUE)
#' guess_lib_from_file(files[3], show_plots = TRUE)
#' 
#' #second strand
#' file.second_strand = dir("inst/extdata", pattern = "test_second.ReadsPerGene.out.tab", full = TRUE)
#' guess_lib_from_file(file.second_strand[1], show_plots = TRUE)
guess_lib_from_file = function(f,
                               show_plots = FALSE,
                               miss_cutoff = 1.2,
                               hit_cutoff = 1.2,
                               split_cutoff = 1.5){
    dt = data.table::fread(f, col.names = c("gene_id", "unstranded", "first", "second"))
    dt_stats = dt[grepl("N_", gene_id)]
    tdt_stats = melt(dt_stats, variable.name = "lib", value.name = "count", id.vars = "gene_id")
    
    p_stats = ggplot(tdt_stats, aes(x = lib, y = count)) +
        geom_bar(stat = "identity") +
        facet_wrap(~gene_id, scales = "free_y") +
        scale_y_continuous(labels = function(x)paste(x/1e6, "M")) +
        labs(title = "Problematic reads per strandedness lib_type")
    
    dt_genes = dt[!grepl("N_", gene_id)]
    tdt_genes = melt(dt_genes, variable.name = "lib", value.name = "count", id.vars = "gene_id")
    
    p_genes = ggplot(tdt_genes, aes(x = lib, y = log10(count+1))) +
        geom_boxplot() +
        labs(title = "Reads mapped to genes per strandedness lib_type")
    
    if(show_plots){
        plot(
            cowplot::plot_grid(ncol = 1, rel_heights = c(1, 12),
                ggplot() + theme_void() + cowplot::draw_text("If first looks like second -> unstranded\nElse, one of first or second, whichever has fewer noFeature and more mapped to genes."),
                cowplot::plot_grid(p_genes, p_stats))    
            )
            
    }
    
    guess_lib(
        dt_stats,
        miss_cutoff = miss_cutoff,
        hit_cutoff = hit_cutoff,
        split_cutoff = split_cutoff
    )
}

load_ReadsPerGene.out.tab = function(f, lib_type){
    stopifnot(lib_type %in% c("unstranded", "first", "second"))
    obs_lib = guess_lib_from_file(f)
    if(obs_lib != lib_type){
        warning("file ", f, " had lib_type of ", obs_lib, ', not ', lib_type, " as expected.")
    }
    
    dt = fread(f, col.names = c("gene_id", "unstranded", "first", "second"))
    
    dt_genes = dt[!grepl("N_", gene_id)]
    dt_genes = dt_genes[, c("gene_id", lib_type), with = FALSE]
    colnames(dt_genes)[2] = sub(".ReadsPerGene.out.tab", "", basename(f))
    dt_genes
}

#' load_matrix_from_ReadsPerGene.out.tab
#'
#' @param files 
#' @param lib_type 
#'
#' @return
#' @export
#'
#' @examples
#' files = dir("inst/extdata", pattern = "^test_[0-9].ReadsPerGene.out.tab", full = TRUE)
#' load_matrix_from_ReadsPerGene.out.tab(files)
load_matrix_from_ReadsPerGene.out.tab = function(files, lib_type = "guess"){
    stopifnot(lib_type %in% c("guess", "unstranded", "first", "second"))
    if(lib_type == "guess"){
        test_files = head(files)
        guesses = unique(sapply(files, guess_lib_from_file))
        if(length(guesses) > 1){
            stop("could not uniquely guess library type, please determine using guess_lib_from_file and manually specify lib_type.")
        }
        lib_type = guesses
    }
    stopifnot(lib_type %in% c("unstranded", "first", "second"))
    
    message("loading files...")
    dtl = pbmcapply::pbmclapply(files, mc.cores = 20, load_ReadsPerGene.out.tab, lib_type = lib_type)
    
    message("assembling matrix...")
    mat = matrix(0, nrow = nrow(dtl[[1]]), ncol = length(dtl))
    rownames(mat) = dtl[[1]]$gene_id
    colnames(mat) = sapply(dtl, function(x)colnames(x)[2])
    mat[,1] = dtl[[1]][[2]]
    if(length(dtl) > 1){
        for(i in seq(2, length(dtl))){
            mat[,i] = dtl[[i]][[2]]
        }
    }
    mat
}

