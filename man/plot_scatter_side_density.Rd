% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions_scatterplot_with_side_density.R
\name{plot_scatter_side_density.assemble}
\alias{plot_scatter_side_density.assemble}
\alias{plot_scatter_side_density.xy}
\title{Title}
\usage{
plot_scatter_side_density.assemble(
  components,
  main_title = "",
  main_title.x = 0.02,
  main_title.y = 0.5,
  main_title.hjust = 0,
  main_title.vjust = 0.5,
  rel_widths = c(2, 1),
  rel_heights = c(1, 2)
)

plot_scatter_side_density.xy(
  xy_data,
  x_,
  y_,
  id_ = "id",
  set_ = "set",
  labs_x = x_,
  labs_y = y_,
  labs_sets = set_,
  main_title = NULL,
  main_title.x = 0.02,
  main_title.y = 0.5,
  main_title.hjust = 0,
  main_title.vjust = 0.5,
  sets.colors = NULL,
  bg.string = "background",
  bg.color = "gray70",
  sets.sizes = 1,
  bg.size = 0.5,
  xlim_ = NULL,
  ylim_ = NULL,
  n_auto_label = 8,
  manual_label = NULL,
  label_size = 2,
  label_color = "black",
  label_use_shadow = TRUE,
  label_use_ggrepel = TRUE,
  label_min.segment.length = 0.5,
  ref_line.x = 0,
  ref_line.x.color = "gray50",
  ref_line.y = 0,
  ref_line.y.color = "gray50",
  ref_line.slope = 1,
  ref_line.slope.color = "black",
  suppress_plot = FALSE,
  rel_widths = c(2, 1),
  rel_heights = c(1, 2)
)
}
\arguments{
\item{components}{list of components returned by plot_scatter_side_density.xy}

\item{main_title}{}

\item{main_title.x}{}

\item{main_title.y}{}

\item{main_title.hjust}{}

\item{main_title.vjust}{}

\item{rel_widths}{}

\item{rel_heights}{}

\item{xy_data}{data.frame with 4 required variables defined by x_, y_, id_, and set_}

\item{x_}{character specifying x attribute, default is "x"}

\item{y_}{character specifying y attribute, default is "y"}

\item{id_}{character specifying id attribute, default is "id"}

\item{set_}{character specifying set attribute, default is "set"}

\item{labs_x}{label for x-axis, default is value of x_}

\item{labs_y}{label for y-axis, default is value of y_}

\item{labs_sets}{label for color legend, default is value of set_}

\item{sets.colors}{named character vector specify color for sets. See scale_color_manual for details.  Default is to use Dark2 from RColorBrewer.}

\item{bg.string}{name for background set (all items)}

\item{bg.color}{color for backgorund set (all items)}

\item{sets.sizes}{point size for items in sets. Default is 1.}

\item{bg.size}{}

\item{xlim_}{x limits. Default of NULL is auto.}

\item{ylim_}{y limits. Default of NULL is auto.}

\item{n_auto_label}{number of most extreme items to label. Default is 8.}

\item{manual_label}{item id values to label. Default is NULL.}

\item{label_size}{text size for items in that get labelled (n_auto_label or manual_label). Default is 2.}

\item{label_color}{text color for items in that get labelled (n_auto_label or manual_label). Default is "black".}

\item{label_use_shadow}{logical to use shadow for items in that get labelled (n_auto_label or manual_label). Default is TRUE.}

\item{label_use_ggrepel}{}

\item{ref_line.x}{}

\item{ref_line.x.color}{}

\item{ref_line.y}{}

\item{ref_line.y.color}{}

\item{ref_line.slope}{}

\item{ref_line.slope.color}{}

\item{suppress_plot}{}
}
\value{
gg object created by cowplot::plot_grid


}
\description{
Title

plot  a scatterplot with sets in color
}
\examples{
library(data.table)
library(ggplot2)
n = 50
xy_data = rbind(
  data.table(x = rnorm(10*n, 0, 1), y = rnorm(10*n, 0, 1), set = "background"),
  data.table(x = rnorm(2*n, 2, 1), y = rnorm(2*n, 0, 1), set = "set1"),
  data.table(x = rnorm(2*n, 0, 1), y = rnorm(2*n, 2, 1), set = "set2"),
  data.table(x = rnorm(2*n, 2, 1), y = rnorm(2*n, 2, 1), set = "set3")
)
xy_data$id = seq_len(nrow(xy_data))

#by default, an assembled plot is output to graphic device
plot_scatter_side_density.xy(xy_data, x_ = "x", y_ = "y")

#a list with assembled plots and components are also returned invisibly for extra customization
#here's an example with lots of extra options used
plot_res = plot_scatter_side_density.xy(
  xy_data, 
  x_ = "x", 
  y_ = "y", 
  suppress_plot = TRUE,
  ref_line.x = c(0, 2), 
  ref_line.y = c(0, 2), 
  ref_line.x.color = c("gray70", "forestgreen"),
  ref_line.y.color = c("gray70", "forestgreen"),
  labs_x = "fc x", 
  labs_y = "fc y", 
  labs_sets = "group", 
  main_title = "an important plot")

plot_res$assembled

#as an example we override the limits for all plot components
new_lim = c(-5, 10)
comp = plot_res$components
comp$scatter = comp$scatter + coord_cartesian(xlim = new_lim, ylim = new_lim)
comp$x_density = comp$x_density + coord_cartesian(xlim = new_lim)
comp$y_density = comp$y_density + coord_flip(xlim = new_lim)
plot_scatter_side_density.assemble(comp, main_title = "an important plot : with new limits")

}
