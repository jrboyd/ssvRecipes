library(ssvRecipes)
library(data.table)
dt = readRDS("test.Rds")
#issue 1) there are duplicated Id that need to be made unique
dt$Id = paste0(dt$Id, "_", seq_len(nrow(dt)))
head(dt)
#issue 2) replicate needs to be delimited from sample
for(i in 1:3){
    colnames(dt) = sub(i, paste0("_", i), colnames(dt))
}
#check unpooled for sanity, no melt
hres1 = ssvHeatmap2(dt, row_ = "Id", fill_ = "val")
plot(hres1)

#check unpooled for sanity, after melt
dt = melt(dt, id.vars = "Id", variable.name = "sample", value.name = "val")
hres2 = ssvHeatmap2(dt, row_ = "Id", column_ = "sample", treatment_ = "treatment", fill_ = "val")
plot(hres2)

#aggregate
#first split sample into treatment and replicate
dt[, c("treatment", "replicate") := tstrsplit(sample, "_")]
dtagg = dt[, .(mean_val = mean(val)), by = .(Id, treatment)]
dtagg$replicate = "1"
hres3 = ssvHeatmap2(dtagg, row_ = "Id",
                    treatment_ = "treatment", fill_ = "mean_val",
                    replicate_ = "replicate", side_plot_type = "bars1")
plot(hres3)

#log2
dtagg[, lg_mean_val := log2(mean_val + 1)]
hres4 = ssvHeatmap2(dtagg, row_ = "Id",
                    treatment_ = "treatment", fill_ = "lg_mean_val",
                    replicate_ = "replicate", side_plot_type = "bars1")
plot(hres4)

#You might prefer swapping replicate and treatment.
dtagg$replicate = "My experiment"
hres5 = ssvHeatmap2(dtagg, row_ = "Id",
                    treatment_ = "replicate", fill_ = "lg_mean_val",
                    replicate_ = "treatment", side_plot_type = "lines1")
plot(hres5)

# a custom sideplot function (new feature i just added)
# plots replicates in heatmap but then pools values for sideplot

# myfun = function(clust, i, fill_, replicate_, treatment_, cluster_, side_plot_colors){
#     # clust = hres5@cluster_data; i = 1; fill_ = "lg_mean_val"; replicate_ = "replicate";
#     # cluster_ = "cluster_id"; side_plot_colors = "black"; treatment_ = "treatment"
#     # yrng = c(0 - 2, max(clust[[fill_]]) + 2)
#     yrng = range(clust[[fill_]])
#     yrng = extendrange(yrng, f = .2)
#     pdt = clust[get(cluster_) == i]
#     pdt = pdt[, .(agg_fill = mean(get(fill_))), by = c(treatment_)]
#     p = ggplot(pdt, aes_string(x = treatment_, y = "agg_fill", group = 1, color = treatment_)) +
#         geom_line(size = 1, lty= 1, color = "darkgray") +
#         geom_point(size = 3) +
#         scale_y_continuous(limits = yrng)
#     p
# }
dt[, lg_val := log2(val + 8)]
# hres6 = ssvHeatmap2(dt, row_ = "Id",
#                     treatment_ = "treatment", fill_ = "lg_val",
#                     replicate_ = "replicate", side_plot_FUN = myfun)
# plot(hres6)

# a custom sideplot function (new feature i just added)
# plots replicates in heatmap but then pools values for sideplot
mybars = function(clust, i, fill_, replicate_, treatment_, cluster_, side_plot_colors){
    # clust = hres5@cluster_data; i = 1; fill_ = "lg_mean_val"; replicate_ = "replicate";
    # cluster_ = "cluster_id"; side_plot_colors = rep("black", 3); treatment_ = "treatment"
    # yrng = c(0 - 2, max(clust[[fill_]]) + 2)
    # yrng = range(clust[[fill_]])
    yrng = c(0, 2^max(clust[[fill_]])-1)/3
    yrng = extendrange(yrng, f = .2)
    pdt = clust[get(cluster_) == i]
    pdt = pdt[, .(agg_fill = mean(2^get(fill_)-1)), by = c(treatment_)]
    p = ggplot(pdt, aes_string(x = treatment_, y = "agg_fill", group = 1, fill = treatment_)) +
        geom_bar(stat = "identity") + theme_void() +
        scale_fill_manual(values = side_plot_colors) +
        # geom_point(size = 3) +
        scale_y_continuous(limits = yrng) +
        annotate("line", c(.5, 3.5), c(0,0)) +
        annotate("line", c(.5, 3.5), rep(round(max(yrng)), 2)) +
        annotate("text", label = round(max(yrng)), x = 2, y = round(max(yrng))*.95, vjust = 1)
    p
}
hres6 = ssvHeatmap2(dt, row_ = "Id",
                    treatment_ = "treatment",
                    fill_ = "lg_val",
                    replicate_ = "replicate",
                    side_plot_FUN = mybars,
                    side_plot_colors = c("cadetblue3", "dodgerblue3", "darkorchid3"),
                    main_widths = c(5, .3, .5, 1, 2))
plot(hres6)

#writing cluster members to file
cmemb = h.cluster_members(hres6)
nr = max(sapply(cmemb, length))
m = matrix("", nrow = nr, ncol = length(cmemb))
for(i in seq_len(length(cmemb))){
    genes = cmemb[[i]]
    m[seq_along(genes), i] = genes
}
write.table(m, file = "cluster_genes.csv", col.names = F, row.names = F, quote = F, sep = ",")

#replotting only the first 3 clusters
cdt = h.cluster_data(hres6)
cdt = cdt[cluster_id %in% 1:3]
# cdt$cluster_id = NULL
hres7 = ssvHeatmap2(cdt, row_ = "Id", main_title = "cluster 1, 2 and 3 from hres6",
                    nclust = 3,
                    treatment_ = "treatment",
                    fill_ = "lg_val",
                    side_plot_colors = c("cadetblue3", "dodgerblue3", "darkorchid3"),
                    replicate_ = "replicate",
                    side_plot_FUN = mybars,
                    main_widths = c(5, .3, .5, 1, 2))
plot(hres7)

