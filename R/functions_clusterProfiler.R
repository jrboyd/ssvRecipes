if(FALSE){
    library(seqsetvis)    
    library(ssvRecipes)
    
    set.seed(0)
    clust_dt = ssvSignalClustering(CTCF_in_10a_profiles_dt)
    clust_dt.assign = unique(clust_dt[, .(id, cluster_id)])
    qgr = CTCF_in_10a_overlaps_gr
    
    
    tx_gr = rtracklayer::import.gff("~/gencode.v36.annotation.gtf",
                                    feature.type = "transcript", format = "gtf")
    tss_gr = GenomicRanges::promoters(tx_gr, 1, 0)
    
    bg_genes = unique(tx_gr$gene_name)
    bg_genes.uni = symbol2uniprot(bg_genes)
    
    
    annotate_gr(qgr, tss_gr, max_dist = 1e4)
    gls = annotate_gr_clusters(qgr, clust_dt.assign, tss_gr, max_dist = 1e4)
    
    cprof.go = my_clusterProfiler_GO(gls, bg_genes = bg_genes, pvalueCutoff = 1, qvalueCutoff = 1)
    cprof.kegg = my_clusterProfiler_KEGG(gls, bg_genes.uni, bg_genes.as_is = TRUE, pvalueCutoff = 1, qvalueCutoff = 1)
    cprof.msig = my_clusterProfiler_MSigDB(gls, pvalueCutoff = 1, qvalueCutoff = 1)
    
    cowplot::plot_grid(
        cprof.go$dotplot,
        cprof.kegg$dotplot,
        cprof.msig$dotplot
    )
}


#' my_annotate
#'
#' @param gr
#' @param tx_gr
#' @param gr_size
#'
#' @return
#' @import IRanges S4Vectors data.table GenomicRanges
#' @export 
#'
#' @examples
#' tx_gr = rtracklayer::import.gff("~/gencode.v36.annotation.gtf",
#' feature.type = "transcript", format = "gtf")
#' tss_gr = GenomicRanges::promoters(tx_gr, 1, 0)
#' 
#' qgr = CTCF_in_10a_overlaps_gr
#'     
#' annotate_gr(qgr, tss_gr, max_dist = 1e4)
annotate_gr = function(gr,
                       tss_gr,
                       max_dist = 5e3){
    
    dists = IRanges::distanceToNearest(tss_gr, gr)
    dists = subset(dists, distance <= max_dist)
    anno_dt = cbind(data.table::as.data.table(tss_gr[S4Vectors::queryHits(dists)])[, .(gene_name, gene_id, transcript_id)],
                    data.table::as.data.table(gr[S4Vectors::subjectHits(dists)]), distance = GenomicRanges::mcols(dists)[[1]])
    anno_dt
}


#' annotate_gr_clusters
#'
#' @param qgr 
#' @param clust_assign A data.frame derived from ssvSignalClustering output (clust_assign = unique(clust_dt[, .(id, cluster_id)]))
#'
#' @return
#' @export
#'
#' @examples
#' clust_dt = ssvSignalClustering(CTCF_in_10a_profiles_dt)
#' clust_dt.assign = unique(clust_dt[, .(id, cluster_id)])
#' qgr = CTCF_in_10a_overlaps_gr
#' 
#' tx_gr = rtracklayer::import.gff("~/gencode.v36.annotation.gtf",
#'                                 feature.type = "transcript", format = "gtf")
#' tss_gr = GenomicRanges::promoters(tx_gr, 1, 0)
#' 
#' annotate_gr(qgr, tss_gr, max_dist = 1e4)
#' gls = annotate_gr_clusters(qgr, clust_dt.assign, tss_gr, max_dist = 1e4)
annotate_gr_clusters = function(qgr, clust_assign, tss_gr, max_dist = 5e3){
    peak_dt = data.table::as.data.table(qgr)
    peak_dt$name = names(qgr)
    mdt = merge(peak_dt, clust_assign, by.x = "name", by.y = "id")
    anno_dt = annotate_gr(GenomicRanges::GRanges(mdt), tss_gr, max_dist = max_dist)
    
    anno_dt = anno_dt[, .(list(unique(gene_name))), by = .(cluster_id)]
    anno_lists = anno_dt$V1
    names(anno_lists) = paste0("cluster_", anno_dt$cluster_id)
    anno_lists = anno_lists[order(anno_dt$cluster_id)]
    names(anno_lists)
    lengths(anno_lists)
    anno_lists
}


#' symbol2uniprot
#'
#' @param x
#'
#' @return
#' @export
#' @import clusterProfiler
#'
#' @examples
symbol2uniprot = function(x, OrgDb = org.Hs.eg.db::org.Hs.eg.db){
    bres = clusterProfiler::bitr(x, fromType = "SYMBOL",
                                 toType = c("UNIPROT"),
                                 OrgDb = OrgDb)
    bres[!duplicated(bres$SYMBOL),]$UNIPROT
}

#' my_clusterProfiler_KEGG
#'
#' @param gene_lists KEGG requires UNIPROT, gene_lists is expected to be SYMBOL.  If already UNIPROT set gene_lists.as_is = TRUE.
#' @param bg_genes KEGG requires UNIPROT, bg_genes is expected to be SYMBOL.  If already UNIPROT set bg_genes.as_is = TRUE.
#' @param force_overwrite
#'
#' @return
#' @export
#'
#' @examples
my_clusterProfiler_KEGG = function(gene_lists, 
                                   bg_genes = NULL, 
                                   gene_lists.as_is = FALSE,
                                   bg_genes.as_is = FALSE,
                                   OrgDb = org.Hs.eg.db::org.Hs.eg.db,
                                   force_overwrite = FALSE, 
                                   bfc = BiocFileCache::BiocFileCache(), 
                                   organism = c("hsa", "mmu")[1],  
                                   pvalueCutoff  = 0.05,
                                   qvalueCutoff  = 0.1){
    if(!gene_lists.as_is){
        gene_lists = lapply(gene_lists, symbol2uniprot, OrgDb = OrgDb)    
    }
    
    if(is.null(bg_genes)){
        bg_genes = unique(unlist(gene_lists))
    }else{
        if(!bg_genes.as_is){
            bg_genes = unique(symbol2uniprot(bg_genes, OrgDb = OrgDb))    
        }
    }
    
    rname_go_dat = digest::digest(list(gene_lists, bg_genes, "kegg", pvalueCutoff, qvalueCutoff))
    rname_go_plot = digest::digest(list(gene_lists, bg_genes, "kegg", "plot", pvalueCutoff, qvalueCutoff))
    
    res = bfcif(bfc, rname_go_plot, force_overwrite = force_overwrite,
                function(){
                    message("calc KEGG res")
                    tryCatch(
                        expr = {
                            ck = bfcif(bfc, rname_go_dat, force_overwrite = force_overwrite,
                                       function(){
                                           message("calc compareCluster")
                                           compareCluster(geneCluster = gene_lists,
                                                          universe      = bg_genes,
                                                          fun = "enrichKEGG",
                                                          organism = organism,
                                                          keyType       = 'uniprot',
                                                          pAdjustMethod = "BH",
                                                          pvalueCutoff  = pvalueCutoff,
                                                          qvalueCutoff  = qvalueCutoff)
                                       })
                            p = dotplot(ck)
                            list(result = ck, dotplot = p)
                        }, error = {
                            function(e){
                                ck = NULL
                                p = ggplot() + annotate("text", x = 0, y = 0, label ="no KEGG results")
                                list(result = ck, dotplot = p)
                            }
                        })
                })
    res
}

#' make_msigdb_TERM2GENE
#' 
#' Consult msigdbr for valid inputs
#' msigdbr::msigdbr_species()
#' msigdbr::msigdbr_collections()
#' 
#' @param species 
#' @param category 
#' @param subcategory 
#'
#' @return
#' @export
#' @import msigdbr
#'
#' @examples
#' # Consult msigdbr for valid inputs
#' # msigdbr::msigdbr_species()
#' # msigdbr::msigdbr_collections()
#' 
#' make_msigdb_TERM2GENE("Mus musculus", "C8", NULL)
#' make_msigdb_TERM2GENE("Mus musculus", "C2", "CP:KEGG")
make_msigdb_TERM2GENE = function(species = "Mus musculus", category = "C8", subcategory = NULL){
    msigdbr_df =  msigdbr::msigdbr(species = species, category = category)
    msigdbr_t2g = as.data.frame(dplyr::distinct(msigdbr_df, gs_name, gene_symbol))
    msigdbr_t2g
}

#' my_clusterProfiler_MSigDB
#'
#' @param gene_lists 
#' @param bg_genes 
#' @param force_overwrite 
#' @param bfc 
#' @param organism 
#' @param pvalueCutoff 
#' @param qvalueCutoff 
#' @param msigdb_species 
#' @param msigdb_category 
#'
#' @return
#' @export
#'
#' @examples
my_clusterProfiler_MSigDB = function(gene_lists, 
                                     bg_genes = NULL, 
                                     force_overwrite = FALSE, 
                                     bfc = BiocFileCache::BiocFileCache(), 
                                     organism = c("hsa", "mmu")[1],  
                                     pvalueCutoff  = 0.05,
                                     qvalueCutoff  = 0.1,
                                     msigdb_species = c("Homo sapiens", "Mus musculus")[1],
                                     msigdb_category ="C8",
                                     msigdb_subcategory = NULL,
                                     msigdb_TERM2GENE = NULL){
    library(msigdbr)
    library(clusterProfiler)
    if(is.null(msigdb_TERM2GENE)){
        message("getting msigdb_TERM2GENE.")
        msigdbr_t2g = make_msigdb_TERM2GENE(species = msigdb_species, category = msigdb_category, subcategory = msigdb_subcategory)    
    }else{
        message("msigdb_TERM2GENE provided, ignoring any other msigb parameters.")
        msigdbr_t2g = msigdb_TERM2GENE
    }
    
    # enricher(gene = gene_symbols_vector, TERM2GENE = msigdbr_t2g, ...)
    
    
    
    rname_go_dat = digest::digest(list(gene_lists, bg_genes, "msigdb", pvalueCutoff, qvalueCutoff, msigdbr_t2g))
    rname_go_plot = digest::digest(list(gene_lists, bg_genes, "msigdb", "plot", pvalueCutoff, qvalueCutoff, msigdbr_t2g))
    
    res = bfcif(bfc, rname_go_plot, force_overwrite = force_overwrite,
                function(){
                    message("calc KEGG res")
                    tryCatch(
                        expr = {
                            ck = bfcif(bfc, rname_go_dat, force_overwrite = force_overwrite,
                                       function(){
                                           message("calc compareCluster")
                                           clusterProfiler::compareCluster(
                                               geneClusters = gene_lists,
                                               fun = "enricher",
                                               TERM2GENE = msigdbr_t2g,
                                               pvalueCutoff  = pvalueCutoff,
                                               qvalueCutoff  = qvalueCutoff
                                           )
                                           
                                       })
                            p = dotplot(ck)
                            list(result = ck, dotplot = p)
                        }, error = {
                            function(e){
                                ck = NULL
                                p = ggplot() + annotate("text", x = 0, y = 0, label ="no KEGG results")
                                list(result = ck, dotplot = p)
                            }
                        })
                })
    res
    
}


#' my_clusterProfiler_fromGenes
#'
#' @param gene_lists
#' @param bg_genes
#' @param force_overwrite
#'
#' @return
#' @export
#' @import org.Hs.eg.db clusterProfiler
#'
#' @examples
my_clusterProfiler_GO = function(gene_lists, 
                                 bg_genes = NULL, 
                                 OrgDb = org.Hs.eg.db::org.Hs.eg.db,
                                 force_overwrite = FALSE, 
                                 bfc = BiocFileCache::BiocFileCache(),  
                                 pvalueCutoff  = 0.05,
                                 qvalueCutoff  = 0.1){
    if(is.null(bg_genes)){
        bg_genes = unique(unlist(gene_lists))
    }
    rname_go_dat = digest::digest(list(gene_lists, bg_genes, "BP", pvalueCutoff, qvalueCutoff))
    rname_go_plot = digest::digest(list(gene_lists, bg_genes, "BP", "plot", pvalueCutoff, qvalueCutoff))
    res = bfcif(bfc, rname_go_plot, force_overwrite = force_overwrite, function(){
        message("calc go res")
        tryCatch(
            expr = {
                ck = bfcif(bfc, rname_go_dat, force_overwrite = force_overwrite, function(){
                    message("calc compareCluster")
                    compareCluster(geneCluster = gene_lists,
                                   universe      = bg_genes,
                                   fun = "enrichGO",
                                   OrgDb = OrgDb,
                                   keyType       = 'SYMBOL',
                                   ont           = "BP",
                                   pAdjustMethod = "BH",
                                   pvalueCutoff  = pvalueCutoff,
                                   qvalueCutoff  = qvalueCutoff)
                })
                p = ck %>% simplify %>% dotplot
                list(result = ck, dotplot = p)
            }, error = {
                function(e){
                    ck = NULL
                    p = ggplot() + annotate("text", x = 0, y = 0, label ="no GO results")
                    list(result = ck, dotplot = p)
                }
            })
    })
    res
}

#' clusterProfiler_table
#'
#' @param cp_res
#' @param clust_id
#'
#' @return
#' @export
#'
#' @examples
my_clusterProfiler_table = function(cp_res, clust_id = NULL){
    cdat = cp_res@compareClusterResult
    cdat = as.data.table(cdat)
    
    if(is.null(clust_id)){
        clust_id = levels(cdat$Cluster)[1]
    }
    if(!clust_id %in% levels(cdat$Cluster)){
        stop("clust_id (", clust_id, ") must be one of ", paste(levels(cdat$Cluster), collapse = ", "))
    }
    if(!clust_id %in% unique(cdat$Cluster)){
        warning("clust_id (", clust_id, ") is valid but has no significant enrichment")
        return(datatable(matrix("empty")))
    }
    
    gs_key = unique(cdat[, .(ID, Description)])
    setkey(gs_key, ID)
    
    col_order = cdat[Cluster == clust_id][order(qvalue)]$Description
    
    if(cp_res@fun == "enrichGO"){
        cmat = cdat[, .(gene_name = tstrsplit(geneID, "/")), by = .(ID,Cluster)]
        cmat$gene_name = unlist(cmat$gene_name)
        ctab =  dcast(cmat[Cluster == clust_id],
                      "ID~gene_name",
                      value.var = "gene_name")#, value.var = length)
    }else if(cp_res@fun == "enrichKEGG"){
        cmat = cdat[, .(uniprot = tstrsplit(geneID, "/")), by = .(ID,Cluster)]
        cmat$uniprot = unlist(cmat$uniprot)
        ctab =  dcast(cmat[Cluster == clust_id],
                      "ID~uniprot",
                      value.var = "gene_name")#, value.var = length)
    }
    
    
    dmat = as.matrix(ctab[,-1])
    dmat = ifelse(is.na(dmat), 0, -1)
    rownames(dmat) = gs_key[.(ctab$ID)]$Description
    
    if(cp_res@fun == "enrichGO"){
        
    }else if(cp_res@fun == "enrichKEGG"){
        colnames(dmat) = bitr(colnames(dmat), 
                              fromType = "UNIPROT", 
                              toType = "SYMBOL", 
                              OrgDb = org.Hs.eg.db::org.Hs.eg.db)$SYMBOL
    }
    
    dmat = t(dmat)
    dmat = dmat[, col_order]
    dmat = dmat[order(rowSums(dmat)),]
    library(DT)
    datatable(dmat, options = list(pageLength = nrow(dmat))) %>% formatStyle(
        T,
        backgroundColor = styleEqual(c(0, -1), c('white', 'blue')),
        color = "00000000"
        
    )
}

