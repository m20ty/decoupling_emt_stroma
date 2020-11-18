# For this we need to already have an all_deconvs object.  So, either run the script in
# deconv.R, or read in a previous one from file:

# all_deconvs <- readRDS('../data_and_figures/all_deconvs.rds')

# We'll also need these:

# library(data.table)
# library(ggplot2)

# expression_data <- fread('../../TCGA_data/tcga_expression_data.csv', key = 'id')
# meta_data <- fread('../../TCGA_data/tcga_meta_data.csv', key = 'id')





boxplot_data <- rbindlist(
    sapply(
        1:length(all_deconvs),
        function(i) {
            data.table(
                name = names(all_deconvs)[i],
                ave_emt = colMeans(
                    expression_data[
                        all_deconvs[[i]]$sample_ids,
                        all_deconvs[[i]]$genes_filtered,
                        with = FALSE
                    ]
                )
            )
        },
        simplify = FALSE
    )
)

ave_emt_boxplot <- ggplot(
    boxplot_data,
    aes(name, ave_emt)
) +
    geom_boxplot() +
    theme(
        axis.text.x = element_text(
            angle = 90,
            hjust = 1,
            vjust = 0.5
        )
    ) +
    labs(
        x = 'Analysis name',
        y = 'Average EMT gene expression'
    )





# Write to PDF:

pdf(
    '../data_and_figures/average_emt_gene_expression_boxplot.pdf',
    width = 12,
    height = 7
)

ave_emt_boxplot

dev.off()





# To look at immune marker expression (I don't think there's a pattern):

cell_type_markers <- fread('../../cell_type_markers.csv')

boxplot_data <- rbindlist(
    sapply(
        1:length(all_deconvs),
        function(i) {
            data.table(
                name = names(all_deconvs)[i],
                ave_immune = colMeans(
                    expression_data[
                        all_deconvs[[i]]$sample_ids,
                        cell_type_markers[
                            category == 'immune' & gene %in% names(expression_data),
                            gene
                        ],
                        with = FALSE
                    ]
                )
            )
        },
        simplify = FALSE
    )
)

ave_immune_boxplot <- ggplot(
    boxplot_data,
    aes(name, ave_immune)
) +
    geom_boxplot() +
    theme(
        axis.text.x = element_text(
            angle = 90,
            hjust = 1,
            vjust = 0.5
        )
    ) +
    labs(
        x = 'Analysis name',
        y = 'Average immune gene expression'
    )
