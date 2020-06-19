
#' get_meta_dt
#' 
#' retrieves metadata data.table from Seurat object with reduction coordinates
#' added.  optional arguments all meta clustres to be specified by combining
#' an existing cluster definition and factor level specification for arbitrary
#' variables.
#'
#' @param seurat_obj an object of class Seurat to extract meta data from
#' @param to_combine 
#' @param cluster_vars cluster variables, will be made factors with levels ordered by decreasing size.
#' @param manual_factors named list, names specify variables and entries specify factor levels to apply.
#' @param reduction string to specify reduction, must be in seurat_obj@reductions
#'
#' @return
#' @export
#'
#' @examples
get_meta_dt = function(seurat_obj, 
                       combine_source = "seurat_clusters",
                       to_combine = list("west_island" = c(15, 6), 
                                         "west_projection" = c(13, 4), 
                                         "lonely_island" = c(11), 
                                         "north_archipelago" = c(9, 8, 14, 12), 
                                         "main_continent" = c(7, 0, 2, 16, 2, 1, 3, 5, 10)),
                       cluster_vars = c("seurat_clusters", "meta_cluster"),
                       manual_factors = list(orig.ident = c("wt", "df4")),
                       reduction = "umap"){
    meta_dt = as.data.table(seurat_obj@meta.data)
    meta_dt$id = seurat_obj@meta.data %>% rownames
    
    stopifnot(reduction %in% names(seurat_obj@reductions))
    umap_dt = as.data.table(seurat_obj@reductions[[reduction]]@cell.embeddings, keep.rownames = TRUE) %>% setnames(., "rn", "id")
    meta_dt = merge(meta_dt, umap_dt, by = "id")
    
    
    
    to_combine_dt = lapply(to_combine, function(x){
        data.table(indi_clust = factor(x, levels = levels(meta_dt[[combine_source]])))
    }) %>% rbindlist(idcol = "meta_cluster")
    setnames(to_combine_dt, "indi_clust", combine_source)
    meta_dt = merge(meta_dt, to_combine_dt, by = combine_source, allow.cartesian = TRUE) %>% unique
    
    for(i in seq_along(manual_factors)){
        meta_dt[[names(manual_factors)[i]]] = factor(meta_dt[[names(manual_factors)[i]]], levels = manual_factors[[i]])    
    }
    
    
    stopifnot(cluster_vars %in% colnames(meta_dt))
    for(cv in cluster_vars){
        lev = meta_dt[, .N, c(cv)][order(N, decreasing = TRUE)][[cv]] %>% as.character
        meta_dt[[cv]] = factor(meta_dt[[cv]], levels = lev)
    }
    meta_dt
}

#' get_rna_dt
#' 
#' retrieves RNA data from a Seurat object for selected genes in a tidy format.
#'
#' @param seurat_obj 
#' @param sel_genes 
#' @param raw_counts 
#' @param assay_name 
#'
#' @return
#' @export
#'
#' @examples
get_rna_dt = function(seurat_obj, sel_genes = NULL, raw_counts = FALSE, assay_name = "RNA"){
    if(is.null(sel_genes)){
        rna_dt = as.data.frame(seurat_obj@assays[[assay_name]]@counts)
    }else{
        len_input = length(sel_genes)
        sel_genes = intersect(sel_genes, seurat_obj@assays[[assay_name]] %>% rownames)
        if(length(sel_genes) != len_input){
            d = len_input - length(sel_genes)
            perc = round(100 * d / len_input, 2)
            warning(perc, "% (", d, " of ", len_input, ") genes discarded to match scRNAseq")
        }
        if(!raw_counts){
            rna_dt = as.data.frame(seurat_obj@assays[[assay_name]][sel_genes, ])    
        }else{
            rna_dt = as.data.frame(seurat_obj@assays[[assay_name]]@counts[sel_genes, ])    
        }
        
        
    }
    rna_dt$gene_name = rownames(rna_dt)
    rna_dt = melt(as.data.table(rna_dt), variable.name = "id", value.name = "expression", id.vars = "gene_name")
    rna_dt
}

#' gene_name_gs2mm
#' 
#' converts a human gene name to mouse gene name ie. YFG -> Yfg
#' this is an informal conversion that mostly works but not ideal.
#'
#' @param x 
#'
#' @return
#' @export
#'
#' @examples
gene_name_hs2mm = function(x){
    x = tolower(x)
    substr(x, 1, 1) = toupper(substr(x, 1, 1))
    x = sub("(?<=[0-9])rik$", "Rik", x, perl = TRUE)
    is_rik = grepl("Rik$", x)
    substr(x[is_rik], 8, 8) = toupper(substr(x[is_rik], 8, 8))
    x
}