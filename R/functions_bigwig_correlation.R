# library(magrittr)
# library(seqsetvis)
# library(data.table)
# library(ggplot2)
# library(ggdendro)
# library(GenomicRanges)

if(FALSE){
    options(mc.cores = 20)
    bw_files = c(
        dir("/slipstream/galaxy/uploads/working/qc_framework/output_waldron_bivalency/", pattern = "R[0-9]$", full.names = TRUE) %>%
            dir(., "_FE.bw", full.names = TRUE),
        dir("/slipstream/galaxy/uploads/working/qc_framework/output_hESC_court/", pattern = "R[0-9]$", full.names = TRUE) %>%
            dir(., "_FE.bw", full.names = TRUE)
    )
    k = grepl("(^CD34)|(^H)|(^Kasumi)|(^DOHH)|(^U937)|(^Nalm)", basename(bw_files))
    bw_files = bw_files[k]
    
    np_files = c(
        dir("/slipstream/galaxy/uploads/working/qc_framework/output_waldron_bivalency/", pattern = "R[0-9]$", full.names = TRUE) %>%
            dir(., "R[0-9]_peaks.narrowPeak", full.names = TRUE),
        dir("/slipstream/galaxy/uploads/working/qc_framework/output_hESC_court/", pattern = "R[0-9]$", full.names = TRUE) %>%
            dir(., "R[0-9]_peaks.narrowPeak", full.names = TRUE)
    )
    k = grepl("(^CD34)|(^H)|(^Kasumi)|(^DOHH)|(^U937)|(^Nalm)", basename(np_files))
    np_files = np_files[k]
    
    olaps = consensus_peaks(np_files, max_peaks_kept = 5000)
    qgr = windowize_peaks(olaps, win_size = 500)
    
    qdt = data.table(file = bw_files)
    qdt[, c("cell", "mark", "rep") := tstrsplit(basename(file), "_", keep = 1:3)]
    qdt$mark = toupper(qdt$mark)
    
    bw_dt = ssvFetchBigwig(qdt, qgr = qgr, win_method = "summary", win_size = 1, return_data.table = TRUE)
    wide_dt = dcast(bw_dt, id~cell+mark+rep, value.var = "y")
    mat = as.matrix(wide_dt[,-1])
    head(mat)
    cor_mat = cor(mat)
    
    ggheatmap.dendro(cor_mat, fill = "correlation", title = "Consensus peak fold-erichment correlation")
}


consensus_peaks = function(np_files, min_coverage = .1, max_peaks_kept = 5000){
    np_grs = seqsetvis::easyLoad_narrowPeak(np_files)
    olaps = seqsetvis::ssvOverlapIntervalSets(np_grs)
    k = rowSums(as.data.frame(mcols(olaps))) > length(np_files)*min_coverage
    olaps = olaps[k]
    if(max_peaks_kept < length(olaps) & max_peaks_kept > 0){
        olaps = sample(olaps, 5000)    
    }
    olaps
}

windowize_peaks = function(rolaps, win_size = 500){
    new_width = ceiling(width(rolaps) / win_size) * win_size
    rolaps= resize(rolaps, new_width, fix = "center")
    names(rolaps) = paste0("region_", seq_along(rolaps))
    swin = slidingWindows(rolaps, win_size, win_size) 
    names(swin) = NULL
    swin = swin %>% unlist
    swin_dt = as.data.table(swin)
    swin_dt$id = names(swin)
    swin_dt[, id := paste0(id, ".window_", seq_len(.N)), .(id)]
    qgr = GRanges(swin_dt)
    qgr
}



ggheatmap.dendro = function(mat, fill = "", title = "", subtitle = ""){
    x <- mat
    dd.col <- as.dendrogram(hclust(dist(x)))
    dd.row <- as.dendrogram(hclust(dist(t(x))))
    dx <- dendro_data(dd.row)
    dy <- dendro_data(dd.col)
    
    ggdend <- function(df) {
        ggplot() +
            geom_segment(data = df, aes(x=x, y=y, xend=xend, yend=yend)) +
            labs(x = "", y = "") + theme_minimal() +
            theme(axis.text = element_blank(), axis.ticks = element_blank(),
                  panel.grid = element_blank())
    }
    
    # x/y dendograms
    px <- ggdend(dx$segments)
    py <- ggdend(dy$segments) + coord_flip()
    
    # heatmap
    col.ord <- order.dendrogram(dd.col)
    row.ord <- rev(order.dendrogram(dd.row))
    xx <- mat[col.ord, row.ord]
    xx_names <- attr(xx, "dimnames")
    df <- as.data.frame(xx)
    colnames(df) <- xx_names[[2]]
    df$row_var <- xx_names[[1]]
    df$row_var <- with(df, factor(row_var, levels=row_var, ordered=TRUE))
    mdf <- reshape2::melt(df, id.vars="row_var")
    mdf$col_var = mdf$variable
    p <- ggplot(mdf, aes(x = col_var, y = row_var)) + geom_tile(aes(fill = value))
    p = p + 
        coord_fixed() + 
        theme(axis.text.x = element_text(angle = 270, hjust = 0, vjust = .5)) +
        scale_fill_gradientn(colours = c("blue", "white", "red"), limits = c(-1, 1)) +
        labs(x = "", y = "", fill = fill, 
             title = title,
             subtitle = subtitle)
    p
}

