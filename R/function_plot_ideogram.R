GIE2COL = c("grey100", "#C1C1C1", "#808080",
            "#404040", "#000000",   "grey0",
            "brown4",  "brown3")
names(GIE2COL) = c("gneg", "gpos25", "gpos50",
                   "gpos75", "gpos100", "gvar",
                   "acen", "stalk" )




#' ssvR_plot_ideogogram
#'
#' @param chr_to_show
#' @param gr_highlights
#' @param gen
#' @param bfc
#' @param ideo_ymin
#' @param ideo_ymax
#' @param highlight_fill
#' @param highlight_color
#' @param highlight_alpha
#'
#' @return
#' @export
#' @import GenomicRanges
#' @import biovizBase
#' @rawNamespace import(data.table, except = c(shift, first, second, last))
#' @import BiocFileCache
#' @examples
#' ssvR_plot_ideogogram(facet_cols = 2)
ssvR_plot_ideogogram = function(gen = "hg38",
                                chr_to_show = paste0("chr", c(1:22, "X", "Y")),
                                gr_highlights = NULL,
                                bfc = BiocFileCache::BiocFileCache(),
                                ideo_ymin = -1,
                                ideo_ymax = 0,
                                highlight_fill = "green",
                                highlight_color = NA,
                                highlight_alpha = .5,
                                facet_cols = 1,
                                facet_by_row = FALSE,
                                print_plot = TRUE){
    grIdeo = bfcif(bfc, rname = paste("ideo", gen), function(){
        biovizBase::getIdeogram(gen)
    })
    grIdeo = subset(grIdeo, as.character(GenomicRanges::seqnames(grIdeo)) %in% chr_to_show)

    full_gr = GenomicRanges::GRanges(
        chr_to_show,
        IRanges::IRanges(1,
                         GenomeInfoDb::seqlengths(grIdeo)[chr_to_show]))
    dfIdeo = data.table::as.data.table(grIdeo)
    chrIdeo = subset(dfIdeo, seqnames %in% chr_to_show)

    gie2col = biovizBase::getBioColor("CYTOBAND")
    chrIdeo$fill = gie2col[as.character(chrIdeo$gieStain)]
    chrIdeo$color = NA
    if(is.character(chr_to_show)){
        chr_to_show = factor(chr_to_show, levels = chr_to_show)
    }
    if(!facet_by_row){
        level.byrow <- function(fac, nc){
            # browser()
            # fac <- factor(vec) #if it is not a factor
            suppressWarnings({mlev <- c(matrix(levels(fac), nrow=nc, byrow=T))})
            mlev = mlev[!duplicated(mlev)]
            # mlev = levels(fac)[order(seq_along(fac) %% ceiling(length(fac)/nc), decreasing = TRUE)]
            factor(fac, levels= c(mlev))
        }
        chr_to_show = sort(level.byrow(chr_to_show, facet_cols))
    }

    chrIdeo$seqnames = factor(chrIdeo$seqnames, levels = chr_to_show)

    #https://stackoverflow.com/questions/12888302/facet-wrap-fill-by-column


    p = ggplot() +
        labs(x = "", y = "") +
        geom_rect(data = chrIdeo[gieStain != "acen"],
                  mapping = aes(xmin = start, xmax = end,
                                ymin = ideo_ymin, ymax = ideo_ymax,
                                fill = fill), color = NA) +
        scale_fill_identity() +
        scale_color_identity() +
        facet_wrap("seqnames", ncol = facet_cols, strip.position = "left") +
        theme(panel.background = element_blank(),
              axis.ticks.y =   element_blank(),
              axis.text.x = element_text(size = 8),
              axis.text.y = element_blank(),
              strip.background = element_blank(),
              plot.margin = margin(t = 1, r = 0, b = 0, l = 0, unit = "pt")) +
        scale_x_continuous(labels = function(x)paste(x/10^6, "Mb")) + guides(fill = "none")
    # add outline
    acen_gr = GenomicRanges::reduce(subset(grIdeo, gieStain == "acen"))
    outline_df = as.data.frame(GenomicRanges::setdiff(full_gr, acen_gr))

    acen_dt = data.table::as.data.table(acen_gr)
    poly_dt = acen_dt[, .(
        xs = c(start, (start+end)/2, end, end, (start+end)/2, start, start),
        ys = c(ideo_ymax, (ideo_ymax + ideo_ymin)/2, ideo_ymax, ideo_ymin, (ideo_ymax + ideo_ymin)/2, ideo_ymin, ideo_ymax)),
        by = .(seqnames)]
    poly_dt$fill = gie2col["acen"]# "darkred"
    poly_dt$color = gie2col["acen"]# "darkred"
    p = p + geom_polygon(data = poly_dt, aes(x = xs, y = ys, fill = fill, color = color))
    p = p + geom_rect(data = outline_df,
                      aes(xmin = start, xmax = end,
                          ymin = ideo_ymin, ymax = ideo_ymax), fill = NA, color = "black")
    # p = Ideogram(grIdeo, subchr = "chr1", which = c(distal, left, right), color = "red", fill = "red", alpha = 1)
    highlight_dt = NULL
    if(!is.null(gr_highlights)){
        highlight_dt = data.table::as.data.table(gr_highlights)
        p = p + geom_rect(data = highlight_dt,
                          aes(xmin = start, xmax = end,
                              ymin = ideo_ymin, ymax = ideo_ymax),
                          fill = highlight_fill, color = highlight_color, alpha = highlight_alpha)
    }

    stains = data.frame(stain = unique(chrIdeo$gieStain))
    levels(stains$stain) = c("gneg", "gpos25", "gpos50", "gpos75", "gpos100", "gvar", "acen", "stalk")
    p_leg = ggplot(stains, aes(xmin = 1, xmax = 1, ymin = 1, ymax = 1, fill = stain)) +
        geom_rect(color = "black")+
        scale_fill_manual(values = gie2col[as.character(stains$stain)])
    if(print_plot) plot(p)
    invisible(list(plot = p, legend = cowplot::get_legend(p_leg), data = list(cytobands = chrIdeo[gieStain != "acen"],
                                                                    centromere = poly_dt,
                                                                    outline = outline_df,
                                                                    highlight = highlight_dt)))
}



#' ssvR_plot_ideogogram_data
#'
#' @param data_dt
#' @param gen
#' @param chr_to_show
#' @param gr_highlights
#' @param bfc
#' @param ideo_ymin
#' @param ideo_ymax
#' @param highlight_fill
#' @param highlight_color
#' @param highlight_alpha
#' @param facet_cols
#' @param facet_by_row
#' @param data_ymin
#' @param data_ymax
#'
#' @return
#' @export
#'
#' @examples
#' library(data.table)
#' dat = fread("~/homer_hic/MCF10A_pooled_PE_tag_dir/pcaOut.PC1.txt")
#' dat = dat[, .(x = (start+end)/2, y = PC1, seqnames = chr)]
#' ssvR_plot_ideogogram_data(dat, facet_cols = 2)
#' ssvR_plot_ideogogram_data(dat, facet_cols = 4, data_ymax = 10)
ssvR_plot_ideogogram_data = function(data_dt,
                                     gen = "hg38",
                                     chr_to_show = paste0("chr", c(1:22, "X", "Y")),
                                     gr_highlights = NULL,
                                     bfc = BiocFileCache::BiocFileCache(),
                                     ideo_ymin = -1,
                                     ideo_ymax = 0,
                                     highlight_fill = "green",
                                     highlight_color = NA,
                                     highlight_alpha = .5,
                                     facet_cols = 1,
                                     facet_by_row = FALSE,
                                     data_ymin = 0, data_ymax = 1,
                                     print_plot = TRUE){
    ideo_res = ssvR_plot_ideogogram(gen = gen,
                             chr_to_show = chr_to_show,
                             gr_highlights = gr_highlights,
                             bfc = bfc,
                             ideo_ymin = ideo_ymin,
                             ideo_ymax = ideo_ymax,
                             highlight_fill = highlight_fill,
                             highlight_color = highlight_color,
                             highlight_alpha = highlight_alpha,
                             facet_cols = facet_cols,
                             facet_by_row = facet_by_row, print_plot = FALSE)
    ideo_res$plot
    dt = data.table::copy(data_dt)
    dt = dt[seqnames %in% chr_to_show]
    dt[, y := y - min(y)]
    dt[, y := y / max(y)]
    dt[, y := y * (data_ymax - data_ymin) + data_ymin]
    bin_size = as.numeric(names(-sort(-table(dt$x[-1] - dt$x[-nrow(dt)])))[1])
    dt[, bin := x / bin_size]

    zero_dt = unique(rbind(
        dt[, .(bin = bin[!bin %in% (bin+1)]-1), by = .(seqnames)],
        dt[, .(bin = bin[!bin %in% (bin-1)]+1), by = .(seqnames)]
    ))

    pdt = rbind(
        dt[, .(seqnames, x, y)],
        zero_dt[, .(seqnames, x = bin * bin_size + .5 * bin_size, y = NA)]
    )[order(x)][order(seqnames)]
    pdt$seqnames = factor(pdt$seqnames, chr_to_show)

    ideo_res$plot = ideo_res$plot + geom_path(data = pdt, aes(x = x, y = y))
    if(print_plot) plot(ideo_res$plot)
    ideo_res$data$line_data = pdt
    invisible(ideo_res)
}



chr_to_show = paste0("chr", c(1:22, "X", "Y"))
ssvR_plot_ideogogram(chr_to_show = chr_to_show, facet_cols = 3, facet_by_row = T)

dt = data.table::fread("~/homer_hic/MCF10A_pooled_PE_tag_dir/pcaOut.PC1.txt")
dt = dt[, .(seqnames = chr, y = PC1, x = (start + end) / 2)]
ssvR_plot_ideogogram_data(dt, chr_to_show = paste0("chr", 1:5), facet_by_row = F, facet_cols = 2)
ssvR_plot_ideogogram_data(dt, chr_to_show = paste0("chr", 10:15), facet_by_row = F, facet_cols = 2)
