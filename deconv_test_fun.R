# This doesn't work. >:-( >:-( >:-( >:-( >:-( >:-(
# I guess it's to do with some scoping issues in the deconvolve_emt_caf_data() function itself.  Otherwise, this function should work, whether using dots = list(...) or
# dots_args <- as.list(match.call( etc.

# EDIT: I got the same error when trying to add CAFs to the <cell_types> argument and using manual cell type weights - I forgot to add CAFs to the list of cell type
# weights.  Does this mean the error is due to a missing cell type, or similar?

deconv_test <- function(deconv_args_per_ct, ct, n, ...) {

    # dots = list(...)
    dots_args <- as.list(match.call())[-(1:4)] # First element will be function name

    deconv_ct <- do.call(
        deconvolve_emt_caf_data,
        args = c(
            list(
                expression_data = expression_data,
                meta_data = meta_data,
                genes = emt_markers,
                cell_type_markers = cell_type_markers,
                # heatmap_annotations = heatmap_annotations,
                ccle_data = ccle,
                subtypes_data = subtypes_data,
                subtypes_var = 'subtype',
                extra_data = extra_data,
                genes_from_tcga_fun = function(x) top_cols_by_fun_cor(x, initial = initial_genes, threshold_fun = function(x) quantile(x, 0.99))$id,
                initial_gene_weights = FALSE,
                genes_filter_fun = function(x) 1:n
            ),
            deconv_args_per_ct[[ct]][!(names(deconv_args_per_ct[[ct]]) == 'plot_title')],
            # dots
            dots_args
        )
    )

    deconv_ct <- deconv_reorder(deconv_ct)

    deconv_plot_ct <- do.call(
        deconvolve_emt_caf_plots,
        args = c(
            list(
                data = deconv_ct,
                # Include the following only if you want epithelial scores (takes much longer):
                # expression_data = expression_data,
                heatmap_axis_title = '', # Change heat_map() function so I can put NULL here
                heatmap_legend_title = 'Coexpression',
                heatmap_colours = rev(colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))(50)),
                heatmap_annotations = heatmap_annotations,
                heatmap_annotations_nudge = 0.3,
                purity_colours = rev(colorRampPalette(RColorBrewer::brewer.pal(11, "PuOr"))(50)),
                purity_colour_limits = c(-0.3, 0.3),
                purity_legend_breaks = c(-0.3, 0, 0.3),
                purity_legend_title = 'Correlation\nwith purity',
                purity_legend_direction = 'vertical',
                purity_axis_title = NULL,
                ccle_colours = rev(colorRampPalette(RColorBrewer::brewer.pal(11, 'PiYG')[c(1:3, 6, 9:11)])(50)),
                ccle_fun = function(x) {caTools::runmean(x - mean(x), 30)/max(abs(caTools::runmean(x - mean(x), 30)))},
                ccle_legend_title = 'Tumours vs.\ncell lines',
                ccle_legend_direction = 'vertical',
                ccle_axis_title = NULL,
                extra_colours = colorRampPalette(c('turquoise4', 'turquoise', 'azure', 'gold2', 'gold4'))(50),
                extra_axis_title = NULL,
                extra_legend_title = 'scRNA-seq',
                extra_legend_direction = 'vertical',
                bar_legend_width = NULL,
                bar_legend_height = NULL
            ),
            paste(deconv_args_per_ct[[str_split(ct, ',', 2)[[1]][1]]]['plot_title'], str_split(ct, ',', 2)[[1]][2], sep = ',')
        )
    )

    pdf('temp.pdf', width = 10, height = 12)
    print(deconv_plot(list(deconv_plot_ct), legends_space = 0.2))
    dev.off()

}





deconv_ct <- do.call(
    deconvolve_emt_caf_data,
    args = c(
        list(
            expression_data = expression_data,
            meta_data = meta_data,
            genes = emt_markers,
            cell_type_markers = cell_type_markers,
            # heatmap_annotations = heatmap_annotations,
            ccle_data = ccle,
            subtypes_data = subtypes_data,
            subtypes_var = 'subtype',
            extra_data = extra_data,
            genes_from_tcga_fun = function(x) top_cols_by_fun_cor(x, initial = initial_genes, threshold_fun = function(x) quantile(x, 0.99))$id,
            initial_gene_weights = FALSE,
            genes_filter_fun = function(x) 1:n
        ),
        deconv_args_per_ct[[ct]][!(names(deconv_args_per_ct[[ct]]) == 'plot_title')],
        mapvalues(
            weights_arg,
            c('90th percentile', 'median', 'none'),
            list(list(gene_weights_fun = function(x) quantile(x, 0.9)), list(gene_weights_fun = median), list(cell_type_weights = FALSE)),
            warn_missing = FALSE
        ) %>% unlist(recursive = FALSE)
    )
)

deconv_ct <- deconv_reorder(deconv_ct)

deconv_plot_ct <- do.call(
    deconvolve_emt_caf_plots,
    args = c(
        list(
            data = deconv_ct,
            # Include the following only if you want epithelial scores (takes much longer):
            # expression_data = expression_data,
            heatmap_axis_title = '', # Change heat_map() function so I can put NULL here
            heatmap_legend_title = 'Coexpression',
            heatmap_colours = rev(colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))(50)),
            heatmap_annotations = heatmap_annotations,
            heatmap_annotations_nudge = 0.3,
            purity_colours = rev(colorRampPalette(RColorBrewer::brewer.pal(11, "PuOr"))(50)),
            purity_colour_limits = c(-0.3, 0.3),
            purity_legend_breaks = c(-0.3, 0, 0.3),
            purity_legend_title = 'Correlation\nwith purity',
            purity_legend_direction = 'vertical',
            purity_axis_title = NULL,
            ccle_colours = rev(colorRampPalette(RColorBrewer::brewer.pal(11, 'PiYG')[c(1:3, 6, 9:11)])(50)),
            ccle_fun = function(x) {caTools::runmean(x - mean(x), 30)/max(abs(caTools::runmean(x - mean(x), 30)))},
            ccle_legend_title = 'Tumours vs.\ncell lines',
            ccle_legend_direction = 'vertical',
            ccle_axis_title = NULL,
            extra_colours = colorRampPalette(c('turquoise4', 'turquoise', 'azure', 'gold2', 'gold4'))(50),
            extra_axis_title = NULL,
            extra_legend_title = 'scRNA-seq',
            extra_legend_direction = 'vertical',
            bar_legend_width = NULL,
            bar_legend_height = NULL
        ),
        paste(deconv_args_per_ct[[str_split(ct, ',', 2)[[1]][1]]]['plot_title'], str_split(ct, ',', 2)[[1]][2], sep = ',')
    )
)

print(deconv_plot(list(deconv_plot_ct), legends_space = 0.2))
