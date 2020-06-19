
#' returns x limits of a ggplot
#'
#' @param p a ggplot
#'
#' @return x limits
#' @export
#'
#' @examples
get_gg_xrange = function(p){
    ggplot2::layer_scales(p)$x$range$range 
}

#' returns y limits of a ggplot
#'
#' @param p a ggplot
#'
#' @return y limits
#' @export
#'
#' @examples
get_gg_yrange = function(p){
    ggplot2::layer_scales(p)$y$range$range 
}

#' synchronizes margins of central scatterplot to align density plots along x and y axis
#' places legend from primary scatter at top-right
#'
#' @param components list with ggplot items: scatter, x_density, and y_density
#' @param main_title optional title to add at top
#' @param main_title.x x coord of title, 0-1, default is 0.02
#' @param main_title.y y coord of title, 0-1, default is 0.5
#' @param main_title.hjust hjust of title, default is 0
#' @param main_title.vjust vjust of title, default is 0.5
#'
#' @return a grob from cowplot::plot_grid
#' @export
#' @import cowplot ggplot2
#' @examples
plot_scatter_side_density.assemble = function(components, main_title = "", main_title.x = .02, main_title.y = .5, main_title.hjust = 0, main_title.vjust = .5){
    p_scatter = components$scatter
    p_x_density = components$x_density
    p_y_density = components$y_density
    
    
    p_legend = cowplot::get_legend(p_scatter)
    
    grobs_y = sync_height(list(p_scatter+guides(color = "none", size = "none"), p_y_density+guides(color = "none")))
    grobs_x = sync_width(list(grobs_y[[1]], p_x_density+guides(color = "none")))
    
    pg = cowplot::plot_grid(plotlist = c(grobs_x[2], list(p_legend), grobs_y), 
                            rel_widths = c(2, 1), rel_heights = c(1,2))
    if(main_title != ""){
        pg = cowplot::plot_grid(
            cowplot::ggdraw() + 
                cowplot::draw_text(main_title, 
                                   x = main_title.x, 
                                   y = main_title.y, 
                                   hjust = main_title.hjust, 
                                   vjust = main_title.vjust),
            pg,
            rel_heights = c(1, 15),
            ncol = 1
        )
    }
    pg
}

library(data.table)
library(ggplot2)

#' plot  a scatterplot with sets in color
#'
#' @param xy_data a data.table with columns determined by x_, y_, id_, set_.
#' alternatively a data.frame or matrix whose rownames will be set to id_.
#' @param x_ string specifiying column in xy_data whose values map to x-axis, default is "x"
#' @param y_ string specifiying column in xy_data whose values map to y-axis, default is "y"
#' @param id_ string specifiying column in xy_data whose values uniquely identify items, default is "id"
#' @param set_ string specifiying column in xy_data whose values specify set membership, default is "set"
#' @param labs_x string labelling x-axis, default is value of x_.
#' @param labs_y string labelling y-axis, default is value of y_.
#' @param labs_sets string labelling color and size, default is value of set_.
#' @param main_title string for main title. default is NULL.
#' @param main_title.x x coord of title, 0-1, default is 0.02
#' @param main_title.y y coord of title, 0-1, default is 0.5
#' @param main_title.hjust hjust of title, default is 0
#' @param main_title.vjust vjust of title, default is 0.5
#' @param sets.colors colors to map to each set.
#' @param bg.string string specifying the background set. default is "background".
#' background is always plotted at the bottom and omitted from density.
#' @param bg.color string specifying color for background.
#' @param sets.sizes sizes of each set. default is 1. repeated to number of sets.
#' @param bg.size size of background points. default is 0.5.
#' @param xlim_ optional manual override of xlimits. default is NULL.
#' @param ylim_ optional manual override of ylimits. default is NULL.
#' @param n_auto_label number of points to autolabel.  picked based on distance from slope 1.
#' @param manual_label manual specification by id_ of points to label.
#' @param label_size size for item labels.
#' @param label_color colors for item labels.
#' @param label_GEOM geom function to use for labels.  default is shadowtext::geom_shadowtext.
#' consider, geom_label, ggrepel::geom_text_repel etc.
#' @param ref_line.x any reference lines to include along x axis
#' @param ref_line.x.color color of x axis reference lines
#' @param ref_line.y any reference lines to include along y axis
#' @param ref_line.y.color color of y axis reference lines
#' @param ref_line.slope partially impelemented.  if any numeric supplied, a line of slope 1 is drawn.
#' @param ref_line.slope.color color for slop line.
#' @param suppress_plot if TRUE, plot is not send to device.  invisible output should be captured.
#'
#' @return a list of two intems, grob of assembled plot and a list of grobs for all components.
#' @export
#' @import ggplot2 cowplot RColorBrewer shadowtext
#' @rawNamespace import(data.table, except = c(shift, first, second, last))
#' @examples
#' n = 50
#' xy_data = rbind(
#'   data.table(x = rnorm(10*n, 0, 1), y = rnorm(10*n, 0, 1), set = "background"),
#'   data.table(x = rnorm(2*n, 2, 1), y = rnorm(2*n, 0, 1), set = "set1"),
#'   data.table(x = rnorm(2*n, 0, 1), y = rnorm(2*n, 2, 1), set = "set2"),
#'   data.table(x = rnorm(2*n, 2, 1), y = rnorm(2*n, 2, 1), set = "set3")
#' )
#' xy_data$id = seq_len(nrow(xy_data))
#' 
#' #by default, an assembled plot is output to graphic device
#' plot_scatter_side_density.xy(xy_data, x_ = "x", y_ = "y")
#' 
#' #a list with assembled plots and components are also returned for extra customization
#' #here's an example with lots of extra options used
#' plots = plot_scatter_side_density.xy(xy_data, x_ = "x", y_ = "y", suppress_plot = TRUE,
#' ref_line.x = c(0, 2), ref_line.y = c(0, 2), 
#' ref_line.x.color = c("gray70", "forestgreen"),
#' ref_line.y.color = c("gray70", "forestgreen"),
#' labs_x = "fc x", labs_y = "fc y", labs_sets = "group", main_title = "an important plot")
#' plots$assembled
#' new_lim = c(-5, 10)
#' comp = plots$components
#' comp$scatter = comp$scatter + coord_cartesian(xlim = new_lim, ylim = new_lim)
#' comp$x_density = comp$x_density + coord_cartesian(xlim = new_lim)
#' comp$y_density = comp$y_density + coord_flip(xlim = new_lim)
#' plot_scatter_side_density.assemble(comp, main_title = "new limits")
#' 
#' xy_df = as.data.frame(xy_data[,1:3])
#' rownames(xy_df) = xy_data$id
#' plot_scatter_side_density.xy(xy_df)
#' 
plot_scatter_side_density.xy = function( xy_data,
                                         x_ = "x",
                                         y_ = "y",
                                         id_ = "id",
                                         set_ = "set",
                                         #labels
                                         labs_x = x_,
                                         labs_y = y_,
                                         labs_sets = set_,
                                         main_title = NULL,
                                         main_title.x = .02,
                                         main_title.y = .5,
                                         main_title.hjust = 0,
                                         main_title.vjust = .5,
                                         #point color and sizing
                                         sets.colors = NULL,
                                         bg.string = "background",
                                         bg.color = "gray70",
                                         sets.sizes = 1,
                                         bg.size = .5,
                                         #limits
                                         xlim_ = NULL,
                                         ylim_ = NULL,
                                         #point labelling
                                         n_auto_label = 8,
                                         manual_label = NULL,
                                         label_size = 6,
                                         label_color = 'black',
                                         label_GEOM = shadowtext::geom_shadowtext,
                                         #reference lines
                                         ref_line.x = 0,
                                         ref_line.x.color = "gray50",
                                         ref_line.y = 0,
                                         ref_line.y.color = "gray50",
                                         ref_line.slope = 1,
                                         ref_line.slope.color = "black",
                                         suppress_plot = FALSE){
    if(is.matrix(xy_data)){
        rn = rownames(xy_data)
        xy_data = as.data.table(xy_data)
        xy_data[[id_]] = rn
    }
    if(is.data.frame(xy_data) & !is.data.table(xy_data)){
        if(is.null(xy_data[[id_]]) & !is.null(rownames(xy_data))){
            rn = rownames(xy_data)
            xy_data = as.data.table(xy_data)
            xy_data[[id_]] = rn
        }else{
            xy_data = as.data.table(xy_data)
        }
    }
    if(is.data.frame(xy_data))
        if(is.null(xy_data[[set_]])){
            stop("set_ : '", set_, "' must be valid column in xy_data, not found!")
        }
    if(is.null(xy_data[[x_]])){
        stop("x_ : '", x_, "' must be valid column in xy_data, not found!")
    }
    if(is.null(xy_data[[y_]])){
        stop("y_ : '", y_, "' must be valid column in xy_data, not found!")
    }
    if(is.null(xy_data[[id_]])){
        stop("id_ : '", id_, "' must be valid column in xy_data, not found!")
    }
    #labels
    if(is.na(main_title) || is.null(main_title)){
        main_title = ""
    }
    stopifnot(is.character(main_title))
    
    if(!is.factor(xy_data[[set_]])){
        xy_data[[set_]] = factor(xy_data[[set_]])
        
    }
    sets.names = levels(xy_data[[set_]])
    sets.len = length(levels(xy_data[[set_]]))
    
    if(is.null(sets.colors)){
        sets.colors = RColorBrewer::brewer.pal(sets.len, "Dark2")
    }
    stopifnot(length(sets.colors) == sets.len)
    if(is.null(names(sets.colors))){
        names(sets.colors) = sets.names
    }
    if(bg.string %in% names(sets.colors) & is.character(bg.color)){
        sets.colors[bg.string] = bg.color
    }
    if(length(sets.sizes == 0)){
        sets.sizes = rep(sets.sizes, sets.len)
    }
    if(is.null(names(sets.sizes))){
        names(sets.sizes) = names(sets.colors)
    }
    stopifnot(names(sets.sizes) == sets.names)
    if(bg.string %in% names(sets.sizes) & is.numeric(bg.size)){
        sets.sizes[bg.string] = bg.size
    }
    
    stopifnot(is.numeric(xy_data[[x_]]))
    stopifnot(is.numeric(xy_data[[y_]]))
    
    if(is.null(xlim_)){
        xlim = range(xy_data[[x_]])
    }else{
        xlim = xlim_
    }
    if(is.null(ylim_)){
        ylim = range(xy_data[[y_]])
    }else{
        ylim = ylim_
    }
    
    xy_data = xy_data[order(get(set_), decreasing = TRUE)]
    gene_o = xy_data[set != bg.string,][order(get(set_), decreasing = TRUE)][order(abs(get(x_) - get(y_)), decreasing = TRUE)][[id_]]
    if(n_auto_label > length(gene_o)) n_auto_label = length(gene_o)
    to_label = c(manual_label, gene_o[seq_len(n_auto_label)])
    
    p_scatter = ggplot(xy_data, aes_string(x = x_, y = y_, 
                                           color = set_, size = set_,
                                           label = id_)) + 
        geom_point() +
        scale_color_manual(values = sets.colors, drop = FALSE) +
        scale_size_manual(values = sets.sizes, drop = FALSE) +
        coord_cartesian(xlim = xlim, ylim = ylim) +
        labs(x = labs_x, y = labs_y, color = labs_sets, size = labs_sets) +
        guides() +
        theme_classic()
    p_x_density = ggplot(mapping = aes_string(x = x_, color = set_)) +
        geom_density(data = xy_data, color = bg.color) +
        geom_density(data = xy_data[get(set_) != bg.string]) +
        scale_color_manual(values = sets.colors, drop = FALSE) +
        coord_cartesian(xlim = xlim) +
        labs(x = "")+
        theme_classic()
    p_y_density = ggplot(mapping = aes_string(x = y_, color = set_)) +
        geom_density(data = xy_data, color = bg.color) +
        geom_density(data = xy_data[get(set_) != bg.string]) +
        scale_color_manual(values = sets.colors, drop = FALSE) +
        coord_flip(xlim = ylim) +
        labs(x = "")+
        theme_classic()
    #add reference lines
    if(is.numeric(ref_line.x)){
        if(length(ref_line.x.color) == 1){
            ref_line.x.color = rep(ref_line.x.color, length(ref_line.x))
        }
        stopifnot(length(ref_line.x.color) == length(ref_line.x))
        for(i in seq_along(ref_line.x)){
            p_scatter = p_scatter + annotate("line", x = rep(ref_line.x[i], 2), y = ylim, color = ref_line.x.color[i], size = .5, linetype = "dashed")
            p_x_density = p_x_density + annotate("line", x = rep(ref_line.x[i], 2), y = get_gg_yrange(p_x_density), color = ref_line.x.color[i], size = .5, linetype = "dashed")
        }
        
    }
    if(is.numeric(ref_line.y)){
        if(length(ref_line.y.color) == 1){
            ref_line.y.color = rep(ref_line.y.color, length(ref_line.y))
        }
        stopifnot(length(ref_line.y.color) == length(ref_line.y))
        for(i in seq_along(ref_line.y)){
            p_scatter = p_scatter + annotate("line", x = xlim, y = rep(ref_line.y[i], 2), color = ref_line.y.color[i], size = .5, linetype = "dashed")
            p_y_density = p_y_density + annotate("line", x = rep(ref_line.y[i], 2), y = get_gg_yrange(p_x_density), color = ref_line.x.color[i], size = .5, linetype = "dashed")    
        }
        
    }
    if(is.numeric(ref_line.slope)){
        if(max(xlim) > max(ylim)){
            right_pt = max(ylim)
        }else{
            right_pt = max(xlim)
        }
        if(min(xlim) < min(ylim)){
            left_pt = min(ylim)
        }else{
            left_pt = min(xlim)
        }
        p_scatter = p_scatter + annotate("line", x = c(left_pt, right_pt), y = c(left_pt, right_pt), color = ref_line.slope.color, size = .5, linetype = "dashed")
    }
    
    #add labels
    if(length(to_label) > 0){
        p_scatter = p_scatter + label_GEOM(data = xy_data[get(id_) %in% to_label], show.legend = FALSE, size = label_size)        
    }
    
    
    
    components = list(scatter = p_scatter, x_density = p_x_density, y_density = p_y_density)
    
    pg = plot_scatter_side_density.assemble(components, 
                                            main_title = main_title,
                                            main_title.x = main_title.x,
                                            main_title.y = main_title.y,
                                            main_title.hjust = main_title.hjust,
                                            main_title.vjust = main_title.vjust)
    
    if(!suppress_plot)
        plot(pg)
    invisible(list(assembled = pg, components = components))
}
