library(GenomicRanges)
library(data.table)
library(ssvRecipes)
library(seqsetvis)
library(magrittr)
library(Seurat)
# source("functions_sc_bam_files.R")
bam_dir = "~/R/scRNA_HSV1/sc_tracks/bams/"
bam_files = dir(bam_dir, pattern = "is_cell.bam$", full.names = TRUE)
parallel::mclapply(bam_files, function(x){
    if(!file.exists(paste0(x, ".bai"))){
        Rsamtools::indexBam(x)
        return(TRUE)
    }else{
        return(FALSE)
    }
})

Rsamtools::idxstatsBam(bam_files[1])

lab_dt = data.table(sample = basename(bam_files))
lab_dt[, name := tstrsplit(sample, "_", keep = 2)]
lab_dt[name == "counts", name := tstrsplit(sample, "_", keep = 1)]
setkey(lab_dt, sample)

qgr_m = GRanges("GRCh38_________MT", IRanges(1, 16569))
qgr_hsv1 = GRanges("JN555585_1_10x_JN555585.1", IRanges(1, 152222))
qgr_hsv1.sm = qgr_hsv1
# start(qgr_hsv1.sm) = 100e3
# end(qgr_hsv1.sm) = 125e3

qgr_hsv1 = qgr_hsv1.sm

qgr_ex = qgr_m

barcode_files = c(
    "~/R/scRNA_HSV1/data/filtered/reseq_acute_counts/filtered_feature_bc_matrix/barcodes.tsv.gz",
    "~/R/scRNA_HSV1/data/filtered/reseq_latent_counts/filtered_feature_bc_matrix/barcodes.tsv.gz",
    "~/R/scRNA_HSV1/data/filtered/reseq_reactivation_counts/filtered_feature_bc_matrix/barcodes.tsv.gz",
    "~/R/scRNA_HSV1/data/filtered/uninfected-12-wm_counts/filtered_feature_bc_matrix/barcodes.tsv.gz",
    "~/R/scRNA_HSV1/data/filtered/uninfected_counts/filtered_feature_bc_matrix/barcodes.tsv.gz"
)

cell_assignments = lapply(barcode_files, function(bc_file){
    cells = fread(bc_file, header = FALSE)
    setnames(cells, "id")
    cells[, id := sub("-.+", "", id)]
    cells$seurat_clusters = "is_cell"
    cells
})

qdt = data.table(file = bam_files)
qdt[, sample := basename(bam_files)]
qdt[, name := tstrsplit(sample, "[\\._]", keep = 2)]


bam_m_l = 
    lapply(seq_along(bam_files), function(sel_bam){
        pbmcapply::pbmclapply(seq_along(qgr_m), 
                              mc.cores = 36, 
                              function(i){
                                  # ssvFetchBam.scRNA(bam_files[sel_bam], qgr_m[i], n_cores = 1,
                                  #                   cell_assignments[[sel_bam]], win_size = 1, 
                                  #                   flip_strand = FALSE, target_strand = "both", splice_strategy = "ignore",
                                  #                   force_barcodes = F, return_data.table = TRUE, min_seq_qual = 30)
                                  ssvFetchBam()
                              }) %>% rbindlist 
    })


bam_hsv_l = 
    lapply(seq_along(bam_files), function(sel_bam){
        pbmcapply::pbmclapply(seq_along(qgr_hsv1), 
                              mc.cores = 36, 
                              function(i){
                                  ssvFetchBam.scRNA(bam_files[sel_bam], qgr_hsv1[i], n_cores = 1, 
                                                    win_method = "summary",
                                                    win_size = 1000,
                                                    summary_FUN = function(x, w)max(x),
                                                    cell_assignments[[sel_bam]], 
                                                    flip_strand = FALSE, target_strand = "both", splice_strategy = "ignore",
                                                    force_barcodes = F, return_data.table = TRUE, min_seq_qual = 30)
                              }) %>% rbindlist 
    })
names(bam_hsv_l) = lab_dt[.(basename(bam_files))]$name
bam_hsv_dt = rbindlist(bam_hsv_l, idcol = "name")

hsv_gr = rtracklayer::import.gff("/slipstream/home/mmariani/projects/hsv1_scrna/genes.gtf")
hsv_dt = as.data.table(subset(hsv_gr, seqnames == "JN555585_1_10x_JN555585.1"))
hsv_dt[, gene_name := sub("JN555585_1_10x_", "", gene_name)]

xlim = c(start(qgr_hsv1), end(qgr_hsv1))
xlim = scales::expand_range(xlim, .05)



p1 = ggplot(bam_hsv_dt, aes(x = (start + end)/2, y = y, color = name, group = name)) + 
    geom_path(size = 1.2) +
    # theme(axis.text.x = element_blank()) +
    coord_cartesian(xlim = xlim, expand = FALSE) +
    facet_wrap(~name, ncol = 1, scales = "free_y") +
    scale_x_continuous(labels = function(x)paste(x/1e3, "kb")) +
    scale_color_brewer(palette = "Dark2") +
    labs(x = "")

p2 = ggplot(hsv_dt, aes(x = (start + end)/2, 
                        y = -.05, label = gene_name, 
                        angle = 90, hjust = 1, vjust = .5)) + 
    geom_rect(aes(xmin = start, xmax = end, ymin = 0, ymax = .1, fill = strand)) +
    geom_text(size = 2) +
    
    theme(panel.grid = element_blank(), 
          panel.background = element_blank(), 
          axis.text.y = element_blank(), 
          axis.ticks.y = element_blank()) +
    coord_cartesian(ylim = c(-1.3, .1), expand = FALSE, xlim = xlim) +
    scale_x_continuous(labels = function(x)paste(x/1e3, "kb")) +
    labs(x = "bp", y = "") 

gl = ssvRecipes::sync_width(list(p1, p2))
pg = cowplot::plot_grid(plotlist = gl, rel_heights = c(4, 1), ncol = 1)
start(qgr_hsv1)
end(qgr_hsv1)
ggsave(paste0("hsv_tracks.", start(qgr_hsv1), "_", end(qgr_hsv1), ".png"), pg, width = 6, height = 6)
