#' Title
#'
#' @param qdt 
#' @param qgr 
#'
#' @return
#' @export
#'
#' @examples
fetch_bw_summary = function(qdt, qgr){
    bw_dt = ssvFetchBigwig(qdt, qgr = qgr, win_method = "summary", win_size = 1, return_data.table = TRUE)
    wide_dt = dcast(bw_dt, id~sample, value.var = "y")
    mat = as.matrix(wide_dt[,-1])
    head(mat)
    cor_mat = cor(mat)
    cor_mat
}

ggdend <- function(df) {
    ggplot() +
        geom_segment(data = df, aes(x=x, y=y, xend=xend, yend=yend)) +
        labs(x = "", y = "") + theme_minimal() +
        theme(axis.text = element_blank(), axis.ticks = element_blank(),
              panel.grid = element_blank())
}

#' Title
#'
#' @param plot_list 
#'
#' @return
#' @export
#'
#' @examples
#' library(ssvRecipes)
#' mat = matrix(runif(100), ncol = 10)
#' plot_list = plot_hclust_heatmap(mat)
#' plot_list$heatmap = plot_list$heatmap + scale_fill_viridis_c(limits = c(-1, 1))
#' plot_hclust_heatmap.assemble(plot_list)
plot_hclust_heatmap.assemble = function(plot_list){
    p = plot_list$heatmap
    if(!is.null(plot_list$x_dendrogram)){
        p = cowplot::insert_xaxis_grob(p, plot_list$x_dendrogram)        
    }
    if(!is.null(plot_list$y_dendrogram)){
        p = cowplot::insert_yaxis_grob(p, plot_list$y_dendrogram)        
    }
    plot(p)
    invisible(p)
}

#' Title
#'
#' @param x 
#'
#' @return
#' @export
#'
#' @examples
#' library(ssvRecipes)
#' mat = matrix(runif(100), ncol = 10)
#' plot_hclust_heatmap(mat)
#' plot_hclust_heatmap(mat, Rowv = FALSE)
#' plot_hclust_heatmap(mat, Colv = FALSE)
#' plot_hclust_heatmap(mat, dendrogram = "none")
#' plot_hclust_heatmap(mat, dendrogram = "row")
#' plot_hclust_heatmap(mat, dendrogram = "column")
#' 
#' @import ggdendro
plot_hclust_heatmap = function(x, 
                               Rowv = TRUE, Colv = TRUE, 
                               dendrogram = c("both", "row", "column", "none")[1]){
    stopifnot(is.matrix(x))
    stopifnot(dendrogram %in% c("both", "row", "column", "none"))
    stopifnot(is.logical(Rowv))
    stopifnot(is.logical(Colv))
    px = NULL
    py = NULL
    if(!Rowv){
        if(dendrogram == "both") dendrogram = "column"
        if(dendrogram == "row") dendrogram = "none"
    }
    if(!Colv & dendrogram %in% c("both", "row")){
        if(dendrogram == "both") dendrogram = "row"
        if(dendrogram == "column") dendrogram = "none"
    }
    if(is.null(colnames(x))){
        colnames(x) = paste0("col", seq_len(ncol(x)))
    }
    if(is.null(rownames(x))){
        rownames(x) = paste0("row", seq_len(nrow(x)))
    }
    
    if(Rowv){
        dd.row <- as.dendrogram(hclust(dist(x)))
        dy <- ggdendro::dendro_data(dd.row)
        py <- ggdend(dy$segments) + coord_flip()
        row.ord <- order.dendrogram(dd.row)    
    }else{
        row.ord = seq_len(nrow(x))
    }
    
    if(Colv){
        dd.col <- as.dendrogram(hclust(dist(t(x))))
        dx <- ggdendro::dendro_data(dd.col)
        px <- ggdend(dx$segments)
        col.ord <- rev(order.dendrogram(dd.col))    
    }else{
        col.ord = seq_len(ncol(x))
    }
    
    # x/y dendograms
    
    # heatmap
    xx <- x[row.ord, col.ord]
    xx_names <- attr(xx, "dimnames")
    df <- as.data.frame(xx)
    colnames(df) <- xx_names[[2]]
    df$row_var <- xx_names[[1]]
    df$row_var <- with(df, factor(row_var, levels=row_var, ordered=TRUE))
    mdf <- reshape2::melt(df, id.vars="row_var")
    mdf$col_var = mdf$variable
    p <- ggplot(mdf, aes(x = col_var, y = row_var)) + geom_tile(aes(fill = value))
    p = p + 
        # coord_fixed() + 
        #scale_fill_gradientn(colours = c("blue", "white", "red"), limits = c(-1, 1)) +
        labs(x = "", y = "", fill = "value") + 
        theme(axis.text.x = element_text(angle = 270, hjust = 0, vjust = .5, size = 6),
              axis.text.y = element_text(size = 6),
              plot.subtitle = element_text(size = 8),
              plot.title = element_text(size = 10),
              panel.background = element_blank(), 
              panel.grid = element_blank()) 
    
    pg_assembly = p
    if(dendrogram %in% c("both", "row")){
        pg_assembly = cowplot::insert_yaxis_grob(pg_assembly, py)    
        
    }
    if(dendrogram %in% c("both", "column")){
        pg_assembly = cowplot::insert_xaxis_grob(pg_assembly, px)    
    }
    plot(p)
    gg_out = list(heatmap = p, x_dendrogram = px, y_dendrogram = py, assembled = pg_assembly)
    invisible(gg_out)
}
