#out goal is to make a figure ready heatmap with aggregated plots on one side



#' ssvH2
#'
#' Stores results from ssvHeatmap2
#'
#' @name ssvH2-class
#' @rdname ssvH2-class
#' @exportClass ssvH2
setClass("ssvH2",
         slots = c(
             final_plot = "list",
             cluster_members = "list",
             cluster_data = "data.table",
             plot_parts_grouped = "list",
             plot_parts_individual = "list"
         )
)

#' Constructor method of ssvH2 Class.
#'
#' @name ssvH2
#' @rdname ssvH2-class
setMethod("initialize", "ssvH2", function(.Object,
                                          final_plot,
                                          cluster_members,
                                          cluster_data,
                                          plot_parts_grouped,
                                          plot_parts_individual){
    .Object@final_plot = final_plot
    .Object@cluster_members = cluster_members
    .Object@cluster_data = cluster_data
    .Object@plot_parts_grouped = plot_parts_grouped
    .Object@plot_parts_individual = plot_parts_individual
    .Object
})

#' plot method for ssvH2
#'
#' @rdname plot-methods
#' @aliases plot,ssvH2-method
#' @export
setMethod("plot", "ssvH2", function(x){
    h.final_plot(x)
})

#' print method for ssvH2
#'
#' @rdname show-methods
#' @aliases show,ssvH2-method
#' @export
setMethod("show", "ssvH2", definition =
              function(object){
                  writeLines(paste(sep = "\n",
                                   "your_hmap = ssvHeatmap2(your_data), your_hmap is an S4 object of class ssvH2.",
                                   "Here's how to use it:",
                                   "  plot(your_hmap) \n    # draws the fully assembled heatmap",
                                   "  h.final_plot(your_hmap) \n    # same as plot(your_hmap)",
                                   "  h.cluster_members(your_hmap) \n    # list of cluster members,\n    # comes from rownames(your_data)",
                                   "  h.cluster_data(your_hmap) \n    # data.table of your_data as \n    # supplied to geom_raster() with cluster_id added",
                                   "  h.plot_parts_grouped(your_hmap) \n    # list of all intermediate\n    # assembly parts of final heatmap. \n    # use cowplot::plot_grid() for customized final assembly",
                                   "  h.plot_parts_individual(your_hmap) \n    # list of all of the \n    # component ggplots before any assembly"))
              }
)

#' final_plot accessor
#'
#' @param ssvH2 output from ssvHeatmap2()
#'
#' @return final assembly for heatmap figure, view with print()
#' @export
#'
#' @examples
#' h.final_plot(ssvHeatmap2(heatmap_demo_matrix))
h.final_plot = function(ssvH2){
    ssvH2@final_plot[[1]]
}

#' cluster members accessor
#'
#' @param ssvH2 output from ssvHeatmap2()
#'
#' @return list of cluster members
#' @export
#'
#' @examples
#' h.cluster_members(ssvHeatmap2(heatmap_demo_matrix))
h.cluster_members = function(ssvH2){
    ssvH2@cluster_members
}

#' clust data.table accessor
#'
#' @param ssvH2 output from ssvHeatmap2()
#'
#' @return data supplied to geom_raster() with cluster assignments
#' @export
#'
#' @examples
#' h.cluster_data(ssvHeatmap2(heatmap_demo_matrix))
h.cluster_data = function(ssvH2){
    ssvH2@cluster_data
}

#' h.plot_parts_grouped
#'
#' @param ssvH2 output from ssvHeatmap2()
#'
#' @return intermediate assembly of heatmap components
#' @export
#'
#' @examples
#' h.plot_parts_grouped(ssvHeatmap2(heatmap_demo_matrix))
h.plot_parts_grouped = function(ssvH2){
    ssvH2@plot_parts_grouped
}

#' h.plot_parts_individual
#'
#' @param ssvH2 output from ssvHeatmap2()
#'
#' @return individual ggplot2 components of heatmap figure.
#' @export
#'
#' @examples
#' h.plot_parts_individual(ssvHeatmap2(heatmap_demo_matrix))
h.plot_parts_individual = function(ssvH2){
    ssvH2@plot_parts_individual
}

#' seqsetvis improved Heatmap
#'
#' An upgrade to setsetvis' ssvSignalHeatmap() function
#' accepts a standard wide matrix with some colname an rowname assumptions.
#'
#' @details
#' colnames must contain a single "_".
#'
#' Condition or treatment groups should
#' precede the "_".
#'
#' A unique replicate identifier should follow the "_".
#' Replicate identifiers should be unique within each treatment group.
#'
#' valid colnames:
#'
#' A_1, A_2, B_1, B_2
#'
#' A_1, A_2, A_3
#'
#' A_1, B_1, C_1, C_2
#'
#' rownames must be set and be unique.
#'
#' @param mat numeric matrix with row and column names set.  see details.
#' @param treatment_ordering order used for columns and side plots.  default NULL
#' causes order of supplied data to be used.
#' @param replicate_ordering order used for sub columns and side plots ticks.
#' default NULL causes order of supplied data to be used.
#' @param treatment_space_size numeric >=0. size of white space separating heatmap
#' groups. Default is 0.1
#' @param treatment_label_size numeric > 0. size of text to label groups with, default is 8.
#' @param main_heights numeric of length 4.  controls size of heatmap
#' body, axis ticks, group labels, and scale portion.
#' @param main_widths numeric of length 5.  controls size of heatmap
#' body, spacer, cluster bars, cluster connector, and aggregate plots
#' @param treatment_ name of conditions/treatment groups.  appears in plots.
#' @param replicate_ name of replicate identifier.  appears in plots.
#' @param column_ name of heatmap column variable - no effect if supplying mat
#' as matrix.
#' @param row_ name of heatmap row variable - no effect if supplying mat as
#' matrix.
#' @param fill_ name of heatmap fill variable - no effect if supplying mat as
#' matrix.
#' @param nclust numeric number of clusters
#' @param main_title title appearing above assembled figure.
#' @param side_plot_type plotting method for aggregate side plots.
#' one of c("lines1", "lines2", "bars1", "bars2").
#' @param side_plot_FUN overrides sideplot type. supplied function is executed
#' once per cluster and must return a ggplot. call is:
#' side_plot_FUN(clust, i, fill_, replicate_, treatment_, cluster_, side_plot_colors)
#' Where clust is a tidy data.table containing 4 variables_, i is cluster id.
#' @param heatmap_colors colors for heatmap, default uses spectral palette
#' from color brewer.
#' @param side_plot_colors colors for side_plot.
#' must be same length as number of treatment groups. default uses dark2
#' palette from  color brewer.
#'
#' @return a ssvH2 object.  print() this object for help.
#' @export
#' @import seqsetvis
#' @rawNamespace import(data.table, except = c(shift, first, second, last))
#' @import cowplot
#' @examples
#' res = ssvHeatmap2(heatmap_demo_matrix)
#' print(res)
#' plot(res)
ssvHeatmap2 = function(
    mat,
    treatment_ordering = NULL,
    replicate_ordering = NULL,
    treatment_space_size = .1,
    treatment_label_size = 8,
    main_heights = c(8, 1, 2, 2),
    main_widths = c(2, .1, .2, .2, 1.8),
    cluster_ = "cluster_id",
    treatment_ = "sample",
    replicate_ = "x",
    column_ = "column",
    row_ = "id",
    fill_ = "y",
    nclust = 5,
    main_title = "Heatmap",
    side_plot_type = c("lines1", "lines2", "bars1", "bars2")[1],
    side_plot_FUN = NULL,
    heatmap_colors = rev(seqsetvis::safeBrew(5, pal = "spectral")),
    side_plot_colors = NULL
){
    if(length(unique(c(cluster_, treatment_, replicate_, column_, row_, fill_))) != 6){
        stop("variable names must be unique:",
             "\n  cluster_ = ",cluster_,
             "\n  treatment_ = ",treatment_,
             "\n  replicate_ = ",replicate_,
             "\n  column_ = ",column_,
             "\n  row_ = ",row_,
             "\n  fill_ = ",fill_
        )
    }
    if(is.data.frame(mat)){
        if(all(c(row_, column_) %in% colnames(mat))){
            #mat is a tidy data.frame
            dt = data.table::as.data.table(mat)
        }else if(all(c(row_, treatment_, replicate_) %in% colnames(mat))){
            dt = data.table::as.data.table(mat)
            dt[[column_]] = paste(sep = "_", dt[[treatment_]], dt[[replicate_]])
        }
        else{
            #mat is a wide data.frame
            dt = data.table::as.data.table(mat)
            stopifnot(is.character(rownames(mat)))
            if(is.null(dt[[row_]])){
                dt[[row_]] = rownames(mat)
            }

            # if(is.null(dt[[column_]])){
            #     stopifnot(all(c(treatment_, replicate_) %in% colnames(mat)))
            #     dt[[column_]] = paste(sep = "_", dt[[treatment_]], dt[[replicate_]])
            # }
            dt = data.table::melt(dt, id.vars = row_, variable.name = column_, value.name = fill_)
        }
    }
    if(is.matrix(mat)){
        #STEP 1 - start here with real data
        #reformat wide matrix to tidy/tall data.table
        dt = data.table::as.data.table(mat)
        dt[[row_]] = rownames(mat)
        dt = data.table::melt(dt, id.vars = row_, variable.name = column_, value.name = fill_)
        #dt is now a tidy data.table
    }

    stopifnot(class(dt)[1] == "data.table")
    stopifnot(all(c(row_, column_) %in% colnames(dt)))
    if(sum(duplicated(paste(dt[[column_]], dt[[row_]]))) > 0){
        stop("There are duplicated values when combining row_ and column_. ",
             "This is most likely due to duplicated row names or ",
             "column names in supplied matrix.")
    }
    #extract grouping info
    if(is.null(dt[[treatment_]])){
        dt[, c(treatment_) := data.table::tstrsplit(get(column_), "_", keep = 1)]
    }
    if(is.null(dt[[replicate_]])){
        dt[, c(replicate_) := data.table::tstrsplit(get(column_), "_", keep = 2)]
    }

    stopifnot(all(c(treatment_, replicate_) %in% colnames(dt)))
    stopifnot(fill_ %in% colnames(dt))



    #STEP 2 - clustering
    #perform clustering
    if(is.null(dt[[cluster_]])){
        clust = seqsetvis::ssvSignalClustering(dt,
                                               nclust = nclust,
                                               column_ = column_,
                                               fill_ = fill_,
                                               row_ = row_,
                                               max_rows = Inf,
                                               max_cols = Inf,
                                               facet_ = "")
    }else{
        clust = dt
    }



    # clust$xmin = as.numeric(clust[[column_]])-1
    # clust[, xmax := xmin + 1]
    #
    # clust$ymin = as.numeric(clust[[row_]])-1
    # clust[, ymax := ymin + 1]

    #STEP 3 - heatmap figure parts


    if(is.null(treatment_ordering)){
        treatment_ordering = unique(clust[[treatment_]])
        if(is.factor(treatment_ordering)){
            treatment_ordering = droplevels(treatment_ordering)
        }
    }else{
        #treatment_ordering must cover all treatment_
        stopifnot(all(sort(unique(clust[[treatment_]])) == sort(treatment_ordering)))
    }
    clust[[treatment_]] = factor(clust[[treatment_]], levels = treatment_ordering)

    if(is.character(clust[[replicate_]]) || is.factor(clust[[replicate_]])){
        if(is.null(replicate_ordering)){
            replicate_ordering = unique(clust[[replicate_]])
            if(is.factor(replicate_ordering)){
                replicate_ordering = droplevels(replicate_ordering)
            }
        }else{
            #replicate_ordering must cover all replicate_
            stopifnot(all(sort(unique(clust[[replicate_]])) == sort(replicate_ordering)))
        }
        clust[[replicate_]] = factor(clust[[replicate_]], levels = replicate_ordering)
    }

    fill_lim = range(clust[[fill_]])
    p_groups = lapply(treatment_ordering, function(g){
        p = ggplot(clust[get(treatment_) == g]) +
            geom_raster(aes_string(x = replicate_, y = row_, fill = fill_)) +
            scale_fill_gradientn(colours = heatmap_colors, limits = fill_lim) +
            theme(legend.position = "bottom") +
            coord_cartesian(expand = FALSE)
        p
    })
    names(p_groups) = paste0("heatmap_", treatment_ordering)
    plegend = get_legend(p_groups[[1]] +  theme(legend.position = "bottom",
                                                legend.direction = "horizontal",
                                                legend.justification = "center"))

    part_names = function(plot){
        gt <- plot_to_gtable(plot)
        gt$layout$name
    }

    get_ggpart = function(plot, part_name){
        gt <- plot_to_gtable(plot)
        panelIndex <- which(grepl(part_name, gt$layout$name))
        # if (length(panelIndex) == 1) {
        # panel <- gt$grobs[[panelIndex]]
        # }
        # else {
        panel <- gt$grobs[panelIndex]
        # }
    }

    ppanels = lapply(p_groups, get_panel)

    pxaxis = lapply(p_groups, function(p){
        get_ggpart(p + theme(axis.line = element_blank()), "axis-b")[[1]]
    })

    pblank = ggplot() + theme_nothing()

    #weight widths by group size
    unique(clust[, .(get(treatment_), get(replicate_))])
    rel_widths = vapply(treatment_ordering, function(grp_)length(unique(clust[get(treatment_) == grp_][[replicate_]])), FUN.VALUE = 1)
    rel_widths = rel_widths / sum(rel_widths)

    plot_spaced_grid = function(plotlist, spacer = .1){
        plist = append(list(pblank), plotlist)
        o = c(2, unlist(lapply(seq_len(length(plotlist) -1), function(i)(c(1, i+2)))))
        plist = plist[o]
        w = c(spacer, rel_widths)[o]
        # w = c(rel_widths[1], unlist(lapply(seq_len(length(plotlist) -1), function(i)(c(spacer, rel_widths)))))
        cowplot::plot_grid(plotlist = plist, nrow = 1, rel_widths = w)
    }

    plot_spaced_labels = function(plabels, spacer = .1, label_size = 8){
        plotlist = lapply(plabels, function(lab){
            ggplot() +
                annotate("text", x = .5, y = 1, label = lab, vjust = 1.5, hjust = .5, size = label_size) +
                annotate("line", x = c(0,1), y = c(1,1)) +
                theme_nothing()
        })
        pblank = ggplot() + theme_nothing()
        plist = append(list(pblank), plotlist)
        o = c(2, unlist(lapply(seq_len(length(plotlist) -1), function(i)(c(1, i+2)))))
        plist = plist[o]
        w = c(spacer, rel_widths)[o]
        cowplot::plot_grid(plotlist = plist, nrow = 1, rel_widths = w)
    }

    pg_plots = plot_spaced_grid(plotlist = ppanels, treatment_space_size)
    pg_xaxis = plot_spaced_grid(plotlist = pxaxis, treatment_space_size)
    pg_labels = plot_spaced_labels(plabels = treatment_ordering, treatment_space_size, treatment_label_size)

    pg_main = cowplot::plot_grid(pg_plots,
                        pg_xaxis,
                        pg_labels,
                        plegend,
                        ncol = 1, rel_heights = main_heights)

    #STEP 4 - cluster bars and aggregated plots

    #bars indicating cluster identity
    #an internal seqsetvis function does most of this already
    pclust = seqsetvis:::add_cluster_annotation(clust[[cluster_]]) +
        coord_cartesian(expand = FALSE) + theme_nothing()

    pg_clust = cowplot::plot_grid(pclust,
                         pblank,
                         pblank,
                         pblank,
                         rel_heights = main_heights, ncol = 1)

    #lines connect cluster bars to line plots
    line_start = c(0, cumsum(table(clust[[cluster_]])))

    line_end = sapply(0:nclust, function(i){
        i * max(line_start) / nclust
    })
    line_df = data.frame(start = line_start, end = line_end)
    pconnect = ggplot(line_df) + theme_nothing() + coord_cartesian(expand = FALSE) + scale_y_reverse()
    for(i in seq_len(nrow(line_df))){
        pconnect = pconnect + annotate("line", x = c(0, 1), y = as.numeric(line_df[i,]), lty = 2)
    }
    pg_connect = plot_grid(pconnect,
                           pblank,
                           pblank,
                           pblank,
                           rel_heights = main_heights, ncol = 1)

    #aggregated line plots
    if(is.null(side_plot_colors)){
        side_plot_colors = seqsetvis::safeBrew(length(treatment_ordering))
    }
    stopifnot(length(side_plot_colors) == length(treatment_ordering))
    if(is.null(names(side_plot_colors))){
        names(side_plot_colors) = treatment_ordering
    }


    gg_agg = function(clust, plot_type = c("lines1", "lines2", "bars1", "bars2")[1]){
        agg = clust[, .(agg_fill = mean(get(fill_))), by = c(treatment_, replicate_, cluster_)]

        agg_lim = c(min(c(0, agg$agg_fill)), max(agg$agg_fill))
        if(is.null(side_plot_FUN)){
            agg_plots = lapply(seq_len(nclust), function(i){
                p = switch(plot_type,
                           lines1 = {
                               ggplot(agg[get(cluster_) == i]) +
                                   geom_line(aes_string(x = replicate_,
                                                        y = "agg_fill",
                                                        color = treatment_,
                                                        group = treatment_)) +
                                   facet_wrap(treatment_, nrow = 1)
                           },
                           lines2 = {
                               ggplot(agg[get(cluster_) == i]) +
                                   geom_line(aes_string(x = replicate_,
                                                        y = "agg_fill",
                                                        color = treatment_,
                                                        group = treatment_))
                           },
                           bars1 = {
                               ggplot(agg[get(cluster_) == i]) +
                                   geom_bar(aes_string(x = replicate_,
                                                       y = "agg_fill",
                                                       fill = treatment_,
                                                       group = treatment_),
                                            stat = "identity") +
                                   facet_wrap(treatment_, nrow = 1)
                           },
                           bars2 = {
                               ggplot(agg[get(cluster_) == i]) +
                                   geom_bar(aes_string(x = replicate_,
                                                       y = "agg_fill",
                                                       fill = treatment_),
                                            stat = "identity", position = "dodge")
                           })
                p = p +
                    scale_y_continuous(limits = agg_lim) +
                    theme(legend.position = "bottom",
                          legend.direction = "horizontal",
                          legend.justification = "center") +
                    scale_color_manual(values = side_plot_colors)
            })
        }else{
            agg_plots = lapply(seq_len(nclust), function(i){
                side_plot_FUN(clust, i, fill_, replicate_, treatment_, cluster_, side_plot_colors)
            })
        }
        return(agg_plots)
    }
    pagg_parts = gg_agg(clust, plot_type = side_plot_type)
    names(pagg_parts) = paste0("clust_agg_", seq_len(nclust))
    pp_agg_leg = get_legend(pagg_parts[[1]])

    # plot = pagg_parts[[1]]
    # gt <- plot_to_gtable(plot)
    # npanels = sum(grepl("panel", gt$layout$name))

    pp_agg_axis = get_ggpart(pagg_parts[[1]] + theme(axis.line = element_blank()), "axis-b")
    pagg = lapply(pagg_parts, function(p){
        p + theme_nothing() + theme(plot.background = element_rect(color = "black"))
    })
    pagg = plot_grid(plotlist = pagg, ncol = 1)
    pg_agg = plot_grid(pagg,
                       plot_grid(plotlist = pp_agg_axis, nrow = 1),
                       pblank,
                       pp_agg_leg,
                       rel_heights = main_heights, ncol = 1)

    #STEP 5 - figure assembly
    pg_all = plot_grid(pg_main,
                       pblank,
                       pg_clust,
                       pg_connect,
                       pg_agg,
                       rel_widths = main_widths, nrow = 1)

    pg_final = plot_grid(pg_all, scale = .9, labels = main_title)

    #STEP 6 - output assembly

    clust_assignemnt = unique(clust[, .(get(row_), get(cluster_))])
    cluster_members = split(clust_assignemnt$V1, clust_assignemnt$V2)
    cluster_members = lapply(cluster_members, as.character)
    new("ssvH2",
        final_plot = list(pg_final),
        cluster_members = cluster_members,
        cluster_data = clust,
        plot_parts_grouped = list(
            heatmap = pg_main,
            blank = pblank,
            cluster_bars = pg_clust,
            cluster_connector = pg_connect,
            cluster_aggregation = pg_agg
        ),
        plot_parts_individual = append(append(p_groups,
                                              list(blank = pblank,
                                                   clust_bars = pclust,
                                                   clust_connector = pconnect)),
                                       pagg_parts)

    )
}


