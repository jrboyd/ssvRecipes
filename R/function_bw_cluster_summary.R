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
#' debug(plot_hclust_heatmap)
#' plot_hclust_heatmap(mat)
#' 
#' @import ggdendro
plot_hclust_heatmap = function(x){
    if(is.null(colnames(x))){
        colnames(x) = paste0("col", seq_len(ncol(x)))
    }
    if(is.null(rownames(x))){
        rownames(x) = paste0("row", seq_len(nrow(x)))
    }
    dd.col <- as.dendrogram(hclust(dist(x)))
    dd.row <- as.dendrogram(hclust(dist(t(x))))
    dx <- ggdendro::dendro_data(dd.row)
    dy <- ggdendro::dendro_data(dd.col)
    
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
    xx <- x[col.ord, row.ord]
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
        #scale_fill_gradientn(colours = c("blue", "white", "red"), limits = c(-1, 1)) +
        labs(x = "", y = "", fill = "value") + 
        theme(axis.text.x = element_text(angle = 270, hjust = 0, vjust = .5, size = 6),
              axis.text.y = element_text(size = 6),
              plot.subtitle = element_text(size = 8),
              plot.title = element_text(size = 10)) 
    p
}
