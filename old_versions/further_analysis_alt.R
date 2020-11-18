library(data.table) # 1.12.8
library(ggplot2) # This is version 3.3.1 on my laptop's WSL setup but 3.3.0 on WEXAC.
library(cowplot) # 1.0.0
library(magrittr) # 1.5
library(plyr) # 1.8.6
library(grid) # 3.6.3
library(ggrepel) # 0.8.2
library(latex2exp) # 0.4.0
library(RColorBrewer) # 1.1.2
library(stringr) # 1.4.0
library(egg) # 0.4.5
library(caTools) # 1.18.0
library(seriation) # 1.2.8
library(gridExtra) # 2.3
library(gtable) # 0.3.0
library(scales) # 1.1.1
library(grDevices) # 3.6.3

source('general_functions.R')
source('tcga_functions.R')

expression_data <- fread('../../TCGA_data/tcga_expression_data.csv', key = 'id')
meta_data <- fread('../../TCGA_data/tcga_meta_data.csv', key = 'id')

expression_data <- expression_data[meta_data[sample_type != 'normal', id]]
meta_data <- meta_data[sample_type != 'normal']

deconv_data <- readRDS('../data_and_figures/deconv_alt_data.rds')
deconv_plots <- readRDS('../data_and_figures/deconv_alt_plots.rds')





# Filtering of deconvs based on diagnostics, namely cell type correlations:

diagnostic_cell_type_lms <- sapply(
    deconv_data,
    function(deconv_ct) {
        with(
            deconv_ct,
            cor_with_initial_and_cell_types[genes_filtered][
                ordering,
                sapply(.SD, function(ct) setNames(lm(ct ~ I(.I/.N))$coeff['I(.I/.N)'], NULL), USE.NAMES = TRUE),
                .SDcols = cell_types
            ]
        )
    },
    USE.NAMES = TRUE
)

# Plot of regression slopes, showing cut-off, and heatmap of regression slopes for all cell types (the cancer types/subtypes commented out either didn't work or
# had too few samples):

plot_mat <- diagnostic_cell_type_lms[
    ,
    c(
        'blca_luminal_infiltrated',
        'blca_luminal_papillary',
        # 'blca_luminal',
        'blca_basal_squamous',
        # 'blca_neuronal',
        'brca_luminal_a',
        'brca_luminal_b',
        'brca_basal_like',
        'brca_her2_enriched',
        # 'brca_normal_like',
        'cesc',
        'coad',
        'esca_ac',
        'esca_escc',
        'hnsc_mesenchymal_basal',
        'hnsc_classical',
        'hnsc_atypical',
        'kich',
        # 'kirc',
        'kirp',
        'lihc',
        'luad_proximal_inflammatory',
        'luad_proximal_proliferative',
        'luad_terminal_respiratory_unit',
        # 'lusc_basal',
        'lusc_classical',
        # 'lusc_primitive',
        'lusc_secretory',
        'ov_differentiated',
        'ov_immunoreactive',
        'ov_mesenchymal',
        'ov_proliferative',
        'paad_basal_moffitt',
        'paad_classical_moffitt',
        'prad',
        'read',
        'stad_cin',
        'stad_ebv',
        # 'stad_gs',
        'stad_msi',
        # 'thca',
        'ucec'
    )
]

colnames(plot_mat) <- nice_names(colnames(plot_mat), TRUE)
rownames(plot_mat) <- gsub('_', ' ', sapply(rownames(plot_mat), function(w) gsub('^[A-Z]|^[a-z]', toupper(str_extract(w, '^[A-Z]|^[a-z]')), w)))

# In the following, we could apply a condition on positive correlation with cell types.  E.g.
# we could use x < -0.1 | x > 0.4, meaning we want to exclude cancer types which have a strong
# positive association with some cell type, with an increase of 0.4 in correlation over all
# the genes.  But I'm not convinced a positive association is a problem.  It means we're
# separating out TME components, even if they're not fibroblasts.

pdf('../data_and_figures/diagnostic_summary_alt.pdf', width = 10, height = 6)

ggarrange(
    ggplot(data.table(index = 1:ncol(plot_mat), slope = sort(apply(plot_mat, 2, min)))[, pass := switch((slope < -0.1) + 1, 'pass', 'fail'), by = index]) +
        geom_col(aes(index, slope, fill = pass, colour = NULL)) +
        scale_fill_manual(values = c('#E78AC3', '#66C2A5')) +
        scale_x_discrete(expand = c(0, 0)) +
        geom_hline(yintercept = -0.1, colour = '#377EB8', linetype = 'dashed', size = 0.75) +
        labs(y = 'Minimum regression slope', fill = NULL) +
        theme_test() +
        theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()),
    heat_map(
        t(plot_mat)[order(apply(plot_mat, 2, min)), ],
        colour_limits = c(-0.4, 0.4),
        legend_breaks = c(-0.4, -0.1, 0.1, 0.4),
        colours = c(colorRampPalette(brewer.pal(11, 'PuOr')[1:4])(15), rep('white', 10), colorRampPalette(brewer.pal(11, 'PuOr')[8:11])(15)),
        axis_text_size = 11,
        legend_title = 'Regression\nslope'
    ) + theme(plot.margin = unit(c(1, 5.5, 5.5, 5.5), 'pt'), panel.border = element_rect(size = 0.5, fill = NA)),
    ncol = 1,
    nrow = 2,
    newpage = FALSE
)

dev.off()





# After some additional manual filtering (examining figures by eye), we narrow it down to
# the below list:

ct_to_keep <- c(
    # 'blca_luminal_infiltrated',
    'blca_luminal_papillary',
    # 'blca_luminal',
    'blca_basal_squamous',
    # 'blca_neuronal',
    'brca_luminal_a',
    'brca_luminal_b',
    'brca_basal_like',
    'brca_her2_enriched',
    # 'brca_normal_like',
    # 'cesc',
    'coad',
    'esca_ac',
    'esca_escc',
    'hnsc_mesenchymal_basal',
    'hnsc_classical',
    'hnsc_atypical',
    'kich',
    # 'kirc',
    # 'kirp',
    'lihc',
    'luad_proximal_inflammatory',
    'luad_proximal_proliferative',
    'luad_terminal_respiratory_unit',
    # 'lusc_basal',
    'lusc_classical',
    # 'lusc_primitive',
    'lusc_secretory',
    'ov_differentiated',
    'ov_immunoreactive',
    'ov_mesenchymal',
    'ov_proliferative',
    'paad_basal_moffitt',
    'paad_classical_moffitt',
    'prad',
    'read',
    'stad_cin',
    'stad_ebv',
    # 'stad_gs',
    'stad_msi',
    # 'thca',
    'ucec'
)

nice_names_for_figure <- c(
    # 'BLCA - Luminal-Infiltrated'
    'BLCA - Luminal-Papillary',
    # 'BLCA - Luminal',
    'BLCA - Basal-Squamous',
    # 'BLCA - Neuronal',
    'BRCA - Luminal A',
    'BRCA - Luminal B',
    'BRCA - Basal-like',
    'BRCA - HER2-enriched',
    # 'BRCA - Normal-like',
    # 'CESC',
    'COAD',
    'ESCA - Adenocarcinoma',
    'ESCA - Squamous',
    'HNSC - Malignant-Basal',
    'HNSC - Classical',
    'HNSC - Atypical',
    'KICH',
    # 'KIRC',
    # 'KIRP',
    'LIHC',
    'LUAD - Squamoid',
    'LUAD - Magnoid',
    'LUAD - Bronchioid',
    # 'LUSC - Basal',
    'LUSC - Classical',
    # 'LUSC - Primitive',
    'LUSC - Secretory',
    'OV - Differentiated',
    'OV - Immunoreactive',
    'OV - Mesenchymal',
    'OV - Proliferative',
    'PAAD - Basal',
    'PAAD - Classical',
    'PRAD',
    'READ',
    'STAD - CIN',
    'STAD - EBV',
    # 'STAD - GS',
    'STAD - MSI',
    # 'THCA',
    'UCEC'
)

deconv_data <- deconv_data[ct_to_keep]
deconv_plots <- deconv_plots[ct_to_keep]

# Change to nice names:
names(deconv_data) <- mapvalues(names(deconv_data), names(deconv_data), nice_names_for_figure)
names(deconv_plots) <- mapvalues(names(deconv_plots), names(deconv_plots), nice_names_for_figure)





# Make figures:

pdf('../data_and_figures/deconv_figures_selected.pdf', width = 16, height = 18)

plot_grid(
    blank_plot(),
    big_deconv_plot(
        c('BRCA', 'COAD', 'HNSC', 'LIHC', 'LUAD', 'LUSC', 'PAAD', 'READ'),
        nice_names_for_figure,
        plots_rel_heights = c(4, 1, 1, 1, 15),
        title_font_size = 16
    ),
    blank_plot(),
    plot_grid(
        plot_grid(
            big_deconv_plot(
                c('BLCA', 'ESCA', 'OV', 'STAD', 'UCEC'),
                nice_names_for_figure,
                plots_rel_heights = c(4, 1, 1, 15),
                title_font_size = 16
            ),
            blank_plot(),
            plot_grid(
                big_deconv_plot('CESC', nice_names_for_figure, plots_rel_heights = c(4, 1, 15), title_font_size = 16),
                blank_plot(),
                nrow = 2,
                ncol = 1,
                rel_heights = c(0.91, 4.81) # 7.09)
            ),
            nrow = 1,
            ncol = 3,
            rel_widths = c(3, 0.2, 1.05)
        ),
        plot_grid(
            blank_plot(),
            plot_grid(
                plotlist = list(
                    get_legend(
                        deconv_plots$LIHC$plots$purity_bar +
							labs(fill = 'Correlation with purity') +
							theme(legend.direction = 'horizontal', legend.justification = 'left') +
							guides(fill = guide_colourbar(title.position = 'right'))
                    ),
                    get_legend(
                        deconv_plots$LIHC$plots$ccle_bar +
							labs(fill = 'Tumours vs. cell lines') +
							theme(legend.direction = 'horizontal', legend.justification = 'left') +
							guides(fill = guide_colourbar(title.position = 'right'))
                    ),
                    get_legend(
                        deconv_plots$LIHC$plots$extra_bar +
							labs(fill = 'scRNA-seq: CAF vs. cancer') +
							theme(legend.direction = 'horizontal', legend.justification = 'left') +
							guides(fill = guide_colourbar(title.position = 'right'))
                    ),
                    get_legend(
                        deconv_plots$LIHC$plots$heatmap +
							labs(fill = 'Correlation') +
							theme(legend.justification = 'left') +
							guides(fill = guide_colourbar(title.position = 'right'))
                    )
                ),
                nrow = 4,
                ncol = 1,
                rel_heights = c(1, 1, 1, 2)
            ),
            nrow = 1,
            ncol = 2,
            rel_widths = c(5, 4)
        ),
        nrow = 2,
        ncol = 1,
        rel_heights = c(5.72, 2.28)
    ),
    # blank_plot(),
    # plot_grid(
    #     big_deconv_plot(
    #         'CESC',
    #         nice_names_for_figure,
    #         plots_rel_heights = c(4, 1, 15),
    #         title_font_size = 16
    #     ),
    #     blank_plot(),
    #     nrow = 2,
    #     ncol = 1,
    #     rel_heights = c(0.91, 4.81) # 7.09)
    # ),
    nrow = 1,
    ncol = 4,
    rel_widths = c(0.2, 3.95, 0.2, 4.25)
    # ncol = 6,
    # rel_widths = c(0.2, 3.95, 0.2, 3, 0.2, 1.05)
)

dev.off()

# Subset of plots for main text figure:

# First remake the plots for these cancer types, so we can change the axis labels:

ct_subset <- c(
    'BRCA - Luminal A',
    'HNSC - Malignant-Basal',
    'PAAD',
    'BLCA - Basal-Squamous',
    'ESCA - Adenocarcinoma' ,
    'UCEC'
    # 'COAD',
    # 'LUSC - Classical',
    # 'OV - Immunoreactive',
    # 'STAD - CIN',
    # 'LIHC'
)

heatmap_annotations_subset <- lapply(
    list(
        `BRCA - Luminal A` = c('CALU', 'COL1A1', 'FAP', 'ITGAV', 'LAMC2', 'NOTCH2', 'PCOLCE', 'PVR', 'THY1'),
        `HNSC - Malignant-Basal` = c('ACTA2', 'CD44', 'COL6A1', 'ITGB1', 'LAMC2', 'POSTN', 'SERPINE1', 'TAGLN', 'TGFBI'),
        `PAAD` = c('AREG', 'CD44', 'COL1A2', 'COL8A1', 'ITGB5', 'LAMC2', 'MMP2', 'PDPN', 'TGFBI'),
        `BLCA - Basal-Squamous` = c('CD44', 'TGFBI', 'LAMC2', 'ITGA2', 'PVR', 'TAGLN', 'COL1A2', 'PCOLCE', 'COL1A1'),
        `ESCA - Adenocarcinoma` = c('AREG', 'COL1A1', 'COL1A2', 'CXCL8', 'FAP', 'ITGA2', 'LAMC2', 'MMP3', 'PVR'),
        `UCEC` = c('ACTA2', 'CALU', 'COL3A1', 'DCN', 'FAP', 'FBN2', 'ITGB1', 'NOTCH2', 'TNC')
        # `COAD` = c('CDH2', 'COL1A2', 'CXCL8', 'FAP', 'ITGA2', 'LAMC2', 'MMP3', 'POSTN', 'PVR'),
        # `LUSC - Classical` = c('CD44', 'COL1A2', 'FAP', 'LAMC2', 'PLOD3', 'POSTN', 'PVR', 'PVR', 'THY1'),
        # `OV - Immunoreactive` = c('ACTA2', 'CALU', 'CDH2', 'COL1A2', 'FAP', 'ITGB5', 'LAMC1', 'MMP3', 'TAGLN'),
        # `STAD - CIN` = c('ACTA2', 'AREG', 'CALD1', 'CXCL8', 'LAMC2', 'MMP3', 'PVR', 'TAGLN', 'THY1'),
        # `LIHC` = c('ACTA2', 'CALU', 'COL3A1', 'ITGAV', 'LAMC1', 'NOTCH2', 'TAGLN', 'THY1', 'WWTR1')
    ),
    c,
    c('SNAI1', 'SNAI2', 'TWIST1', 'VIM', 'ZEB1', 'ZEB2')
)

deconv_plots_subset <- sapply(
    ct_subset,
    function(ct) {
        do.call(
            deconvolve_emt_caf_plots,
            args = list(
                data = deconv_data[[ct]],
                plot_title = ct,
                heatmap_legend_title = 'Correlation',
                heatmap_annotations = heatmap_annotations_subset[[ct]],
                heatmap_annotations_nudge = 0.3,
                purity_colours = rev(
                    colorRampPalette(
                        brewer.pal(11, "PuOr")
                    )(50)
                ),
                purity_legend_title = 'Correlation with purity\n',
                purity_legend_direction = 'horizontal',
                purity_axis_title = NULL,
                ccle_colours = rev(
                    colorRampPalette(
                        brewer.pal(11, 'PiYG')[c(1:3, 6, 9:11)]
                    )(50)
                ),
                ccle_legend_title = 'Tumours vs. cell lines\n',
                ccle_legend_direction = 'horizontal',
                ccle_axis_title = NULL,
                extra_colours = colorRampPalette(
                    c('turquoise4', 'turquoise', 'azure', 'gold2', 'gold4')
                )(50),
                extra_legend_title = 'scRNA-seq: CAF vs. cancer\n',
                extra_legend_direction = 'horizontal',
                extra_axis_title = NULL,
                bar_legend_justification = 'left',
                bar_legend_width = NULL
            )
        )
    },
    simplify = FALSE,
    USE.NAMES = TRUE
)

# Summary figures:

# Scatterplot of within- versus between-cluster correlations:

within_between_clust_corr <- lapply(
    names(deconv_data),
    function(ct) {

        within_clust_corr <- with(
            deconv_data[[ct]],
            c(
                cor_mat[ordering, ordering][1:30, 1:30],
                cor_mat[rev(ordering), rev(ordering)][1:30, 1:30]
            )
        )

        within_clust_corr <- mean(within_clust_corr[within_clust_corr != 1])

        between_clust_corr <- with(
            deconv_data[[ct]],
            mean(cor_mat[ordering, rev(ordering)][1:30, 1:30])
        )

        list(
            cancer_type = ct,
            wthn = within_clust_corr,
            btw = between_clust_corr
        )

    }
) %>% rbindlist

within_between_clust_corr[
    ,
    n_annot := sum(
        c('purity_bar', 'ccle_bar', 'extra_bar') %in% names(deconv_plots[[cancer_type]]$plots)
    ),
    by = cancer_type
]

# Make the actual plots:

deconv_summary_scatterplot <- ggplot(
    within_between_clust_corr,
    aes(
        x = btw,
        y = wthn,
        shape = as.character(mapvalues(n_annot, c(1, 2, 3), c('No', 'No', 'Yes')))
    )
) +
    geom_point(size = 2, alpha = 0.75, colour = 'lightblue4') +
    # geom_text_repel(aes(label = cancer_type)) +
    geom_text_repel(
        aes(label = cancer_type),
        data = within_between_clust_corr[cancer_type %in% c('BLCA - Basal-Squamous', 'ESCA - Squamous')],
        nudge_x = 0.04, nudge_y = 0.01, direction = 'y', size = 3.5
    ) +
    geom_text_repel(
        aes(label = cancer_type),
        data = within_between_clust_corr[cancer_type %in% c('HNSC - Malignant-Basal', 'HNSC - Classical', 'UCEC')],
        direction = 'x',
        nudge_x = 0.02, nudge_y = 0.005, direction = 'x', size = 3.5
    ) +
    geom_text_repel(
        aes(label = cancer_type),
        data = within_between_clust_corr[cancer_type == 'COAD'],
        nudge_x = 0.02, nudge_y = 0.02, size = 3.5
    ) +
    geom_text_repel(
        aes(label = cancer_type),
        data = within_between_clust_corr[cancer_type == 'LIHC'],
        nudge_x = -0.02, nudge_y = 0.01, size = 3.5
    ) +
    geom_text_repel(
        aes(label = cancer_type),
        data = within_between_clust_corr[cancer_type == 'STAD - MSI'],
        nudge_x = -0.04, nudge_y = -0.02, size = 3.5
    ) +
    geom_text_repel(
        aes(label = cancer_type),
        data = within_between_clust_corr[cancer_type == 'BRCA - Luminal A'],
        nudge_x = -0.05, nudge_y = 0.01, size = 3.5
    ) +
    geom_text_repel(
        aes(label = cancer_type),
        data = within_between_clust_corr[cancer_type == 'ESCA - Adenocarcinoma'],
        nudge_x = -0.06, nudge_y = -0.01, size = 3.5
    ) +
    geom_text_repel(
        aes(label = cancer_type),
        data = within_between_clust_corr[cancer_type == 'PAAD'],
        nudge_x = -0.02, nudge_y = -0.005, size = 3.5
    ) +
    theme_test() +
    labs(
        x = 'Between-cluster correlation',
        y = 'Within-cluster correlation',
        shape = 'scRNA-seq\ndata\navailable'
    )

pdf(
    '../data_and_figures/deconv_figures_main.pdf',
    width = 14,
    height = 9
)

plot_grid(
    deconv_plot(
        deconv_plots_subset,
        n_row = 2,
        n_col = 3,
        plots_rel_heights = c(
            title = 2,
            purity_bar = 1,
            ccle_bar = 1,
            extra_bar = 1,
            heatmap = 15,
            axis_labels = 6
        ),
        rows_rel_heights = c(1.04, 1),
        left_plot_width = 1.06,
        legends = FALSE
    ),
    plot_grid(
        plotlist = c(
            list(
                deconv_summary_scatterplot +
                    theme(
                        legend.position = 'bottom',
                        legend.text = element_text(margin = margin(r = 15, unit = 'pt')),
                        legend.title = element_text(margin = margin(r = 10, unit = 'pt')),
                        legend.spacing.x = unit(1, 'pt')
                    ) +
                    labs(shape = 'scRNA-seq data available:'),
                blank_plot()
            ),
            lapply(
                c('purity_bar', 'ccle_bar', 'extra_bar', 'heatmap'),
                function(b) {
                    get_legend(
                        deconv_plots_subset$`BRCA - Luminal A`$plots[[b]] + theme(
                            legend.justification = 'left'
                        ) + guides(fill = guide_colourbar(title.position = 'right'))
                    )
                }
            ),
            list(blank_plot())
        ),
        nrow = 7,
        ncol = 1,
        rel_heights = c(26, 2, 3, 3, 3, 8, 6)
    ),
    nrow = 1,
    ncol = 2,
    rel_widths = c(2, 1)
)

dev.off()

# main_deconv_plot <- plot_grid(
#     plot_grid(
#         deconv_plot(
#             deconv_plots_subset,
#             n_row = 2,
#             n_col = 3,
#             plots_rel_heights = c(
#                 title = 2,
#                 purity_bar = 1,
#                 ccle_bar = 1,
#                 extra_bar = 1,
#                 heatmap = 15,
#                 axis_labels = 6
#             ),
#             rows_rel_heights = c(1.04, 1),
#             left_plot_width = 1.06,
#             legends = FALSE
#         ),
#         blank_plot(),
#         nrow = 1,
#         ncol = 2,
#         rel_widths = c(14, 1)
#     ),
#     plot_grid(
#         blank_plot(),
#         plot_grid(
#             plotlist = lapply(
#                 c('purity_bar', 'ccle_bar', 'extra_bar'),
#                 function(b) {
#                     get_legend(
#                         deconv_plots_subset$`BRCA - Luminal A`$plots[[b]] + theme(
#                             legend.justification = 'left'
#                         ) + guides(fill = guide_colourbar(title.position = 'right'))
#                     )
#                 }
#             ),
#             nrow = 3,
#             ncol = 1
#         ),
#         get_legend(
#             deconv_plots_subset$`BRCA - Luminal A`$plots$heatmap +
#                 theme(legend.justification = 'right') +
#                 guides(fill = guide_colourbar(title.position = 'right'))
#         ),
#         blank_plot(),
#         get_legend(
#             deconv_summary_barplot +
#                 theme(legend.justification = 'right') +
#                 guides(fill = guide_legend(title.hjust = 1, label.position = 'left'))
#         ),
#         nrow = 1,
#         ncol = 5,
#         rel_widths = c(0.1, 3, 1, 1, 2)
#     ),
#     nrow = 2,
#     ncol = 1,
#     rel_heights = c(5, 1)
# )

# pdf(
#     '../data_and_figures/deconv_figures_main.pdf',
#     width = 14,
#     height = 10
# )
#
# plot_grid(
#     main_deconv_plot,
#     plot_grid(
#         deconv_summary_scatterplot,
#         deconv_summary_barplot + theme(legend.position = 'none'),
#         nrow = 2,
#         ncol = 1,
#         rel_heights = c(1, 2)
#         # axis = 'l',
#         # align = 'v'
#     ),
#     nrow = 1,
#     ncol = 2,
#     rel_widths = c(2, 1)
# )
#
# dev.off()





# Barplots of agreement with annotations, for supplement:

annotation_agreement <- lapply(
    names(deconv_data),
    function(ct) {

        annots <- c(ccle = 'ccle_comp_diff', extra = 'extra_data_score')

        annots <- annots[annots %in% names(deconv_data[[ct]])]

        # I don't actually use the lm stuff in the end, but it's here in
        # case it becomes useful...

        lms_data <- with(
            deconv_data[[ct]],
            data.table(
                index = 1:length(ordering),
                pur = cor_with_purity$scale[ordering]
            )
        )

        if(length(annots) > 0) {
            lms_data <- cbind(
                lms_data,
                sapply(
                    annots,
                    function(annot) {
                        with(deconv_data[[ct]], get(annot)[ordering])
                    },
                    simplify = FALSE,
                    USE.NAMES = TRUE
                ) %>% as.data.table
            )
        }

        lm_coeffs <- lms_data[
            ,
            c(
                list(cancer_type = ct),
                sapply(
                    .SD,
                    function(x) lm(x ~ index)$coeff['index'],
                    simplify = FALSE,
                    USE.NAMES = TRUE
                )
            ),
            .SDcols = -'index'
        ]

        names(lm_coeffs) <- str_split_fixed(names(lm_coeffs), '\\.', 2)[, 1]

        annot_fun <- function(x) {runmean(x, 30)/max(abs(runmean(x, 30)))}

        corr_diff <- with(
            deconv_data[[ct]],
            list(
                cancer_type = ct,
                pur = mean(head(cor_with_purity$scale[ordering], 30)) -
                    mean(tail(cor_with_purity$scale[ordering], 30))
            )
        )

        if(length(annots) > 0) {
            corr_diff <- c(
                corr_diff,
                sapply(
                    annots,
                    function(annot) {
                        with(
                            deconv_data[[ct]],
                            mean(head(annot_fun(get(annot)[ordering]), 30)) -
                                mean(tail(annot_fun(get(annot)[ordering]), 30))
                        )
                    },
                    simplify = FALSE,
                    USE.NAMES = TRUE
                )
            )
        }

        list(
            lm_coeffs = lm_coeffs,
            corr_diff = corr_diff
        )

    }
)

annotation_agreement <- sapply(
    c('lm_coeffs', 'corr_diff'),
    function(summary_type) {
        rbindlist(
            sapply(annotation_agreement, `[[`, summary_type),
            fill = TRUE
        )[
            ,
            n_annot := sum(!is.na(as.numeric(.SD))),
            by = cancer_type
        ]
    },
    simplify = FALSE,
    USE.NAMES = TRUE
)

barplot_data <- copy(annotation_agreement$corr_diff)[
    ,
    c('pur', 'ccle', 'extra') := lapply(
        .SD,
        function(x) {
            # x[!is.na(x)] <- scale(x[!is.na(x)], center = FALSE)
            x[!is.na(x)] <- x[!is.na(x)]/mean(x[!is.na(x)])
            x
        }
    ),
    .SDcols = c('pur', 'ccle', 'extra')
]

barplot_data <- melt(
    barplot_data,
    id.vars = c('cancer_type', 'n_annot'),
    variable.name = 'annot',
    value.name = 'score'
)[
    ,
    n_annot := factor(n_annot, levels = c(3, 2, 1))
]

barplot_data <- barplot_data[!is.na(score)]

pdf(
    '../data_and_figures/deconv_summary_barplot.pdf',
    width = 7,
    height = 8
)

ggplot(
    barplot_data,
    aes(
        x = factor( # Have to do this here or else order isn't preserved in all facets
            cancer_type,
            levels = barplot_data[
                ,
                .(mean_score = mean(score)),
                by = cancer_type
            ][
                order(mean_score),
                cancer_type
            ]
        ),
        y = score,
        group = annot,
        colour = annot,
        fill = annot
    )
) +
    # geom_col(position = 'dodge') + # position_dodge2() makes bars in all facets the same width
    geom_col(position = position_dodge2(width = 0.9, preserve = 'single')) +
    geom_hline(yintercept = 0, colour = 'lightgrey') +
    coord_flip() +
    facet_grid(
        rows = vars(n_annot),
        space = 'free',
        scales = 'free_y'
    ) +
    theme_half_open() +
    theme(
        strip.text = element_blank(),
        strip.background = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.title.x = element_text(size = 12),
        legend.title = element_text(size = 12),
        legend.text = element_text(size = 11)
    ) +
    scale_colour_manual(
        values = c(
            'pur' = brewer.pal(11, "PuOr")[2],
            'ccle' = brewer.pal(11, "PiYG")[3],
            'extra' = 'gold2'
        ),
        guide = FALSE
    ) +
    scale_fill_manual(
        labels = c(
            'pur' = 'Correlation with purity',
            'ccle' = 'Tumours vs. cell lines',
            'extra' = 'scRNA-seq: CAF vs. cancer'
        ),
        values = c(
            'pur' = brewer.pal(11, "PuOr")[2],
            'ccle' = brewer.pal(11, "PiYG")[3],
            'extra' = 'gold2'
        )
    ) +
    labs( # Remember the axes are flipped!
        x = NULL,
        y = 'Between-cluster difference',
        fill = 'Annotation type'
    )

dev.off()





# Commonality heatmap:

rank_mat <- deconv_rank(deconv_data)

scores_data_transformed <- deconv_scores(
    expression_data,
    deconv_data,
    scale_fun = function(x) x/(3*sd(x)),
    transform_data = TRUE
)

deconv_names <- names(deconv_data)

# htmp_emt_caf_initial <- deconv_heatmap(
#     scores_data_transformed[
#         unique(c(names(head(sort(apply(rank_mat, 1, quantile, 0.25)), 50)), names(tail(sort(apply(rank_mat, 1, quantile, 0.75)), 50)))),
#         c('gene', ..deconv_names)
#     ],
#     order_genes_fun = 'hclust',
#     order_analyses_fun = 'hclust',
#     plot_title = 'Common EMT and stroma genes'
# )
#
# deconv_heatmap_dendro_plot(htmp_emt_caf_initial)

# Based on the above initial hierarchically-clustered heatmap, we'll remove KICH and PRAD, which look like outliers:
deconv_data <- deconv_data[!(names(deconv_data) %in% c('KICH', 'PRAD'))]
deconv_plots <- deconv_plots[!(names(deconv_plots) %in% c('KICH', 'PRAD'))]
deconv_names <- names(deconv_data)

# htmp_emt_caf_revised <- deconv_heatmap(
#     scores_data_transformed[
#         unique(c(names(head(sort(apply(rank_mat, 1, quantile, 0.25)), 50)), names(tail(sort(apply(rank_mat, 1, quantile, 0.75)), 50)))),
#         c('gene', ..deconv_names)
#     ],
#     order_genes_fun = 'hclust',
#     order_analyses_fun = 'hclust',
#     plot_title = 'Common EMT and stroma genes'
# )
#
# deconv_heatmap_dendro_plot(htmp_emt_caf_revised)

# set.seed(78548)
set.seed(44398)

htmp_emt_caf <- deconv_heatmap(
    scores_data_transformed[
        unique(c(names(head(sort(apply(rank_mat, 1, quantile, 0.25)), 50)), names(tail(sort(apply(rank_mat, 1, quantile, 0.75)), 50)))),
        c('gene', ..deconv_names)
    ],
    order_genes_fun = 'seriate',
    order_genes_method = 'SPIN_NH',
    order_analyses_fun = function(x) rev(get_order(seriate(dist(t(x)), method = 'SPIN_NH'))),
    # order_analyses_fun = 'seriate',
    # order_analyses_method = 'SPIN_NH',
    plot_title = 'Common EMT and stroma genes'
)

# PCA to describe the variation in EMT and CAF genes between cancer types:
scores_pca <- prcomp(t(scores_data_transformed[names(head(sort(apply(rank_mat, 1, quantile, 0.25)), 100)), ..deconv_names]))
pca_data <- as.data.table(scores_pca$x[, 1:2], keep.rownames = 'cancer_type')
setkey(pca_data, cancer_type)

# qplot(PC1, PC2, data = pca_data)

# k-means clustering with k = 3, just because to the eye it looks like there are 3:
set.seed(16918)
kclust <- kmeans(pca_data[, .(PC1, PC2)], 3)
pca_data[, kmeans_cluster := kclust$cluster]

# Make the names of the cluster match "1", "2", "3" in a consistent order:
kclust_map <- mapvalues(1:3, 1:3, pca_data[c('BRCA - Luminal A', 'HNSC - Malignant-Basal', 'ESCA - Adenocarcinoma'), kmeans_cluster])
pca_data[, kmeans_cluster_manual := mapvalues(kmeans_cluster, 1:3, kclust_map)]

# Compute segments to draw on a scatterplot to illustrate the 3 directions in the points:
kclust_centre <- colMeans(kclust$centers)
segment_endpoints <- as.data.table(cbind(index = c(1, 2, 3), kclust$centers[kclust_map, ]))[
    ,
    c('PC1_end', 'PC2_end') := as.list(extend_endpoint(kclust_centre, as.numeric(.SD), pca_data[, .(PC1, PC2)])),
    by = index
]

# Transform data once for each of the 3 clusters, so that that cluster is aligned with
# the x axis:

for(i in 1:3) {
    pca_data <- cbind(
        pca_data,
        transform_segment(kclust_centre, kclust$centers[kclust_map[i], ], pca_data[, .(PC1, PC2)], suffix = paste0('_', i))
    )
}

# Arbitrarily choose cutoffs for x axis such that by eye they appear to distinguish the
# cancer types most strongly associated with this "direction" in the PC space:
cutoffs <- c(0.3, 0.8, 0.8)
pca_data[
    ,
    thresh_pass := get(paste0('PC1_', unique(kmeans_cluster_manual))) > cutoffs[unique(kmeans_cluster_manual)],
    by = kmeans_cluster_manual
]

kmeans_clust_distinct_genes <- lapply(
    1:3,
    function(i) {

        cts <- list(
            in_i = pca_data[thresh_pass == TRUE & kmeans_cluster_manual == i, cancer_type],
            not_i = pca_data[thresh_pass == TRUE & kmeans_cluster_manual != i, cancer_type]
        )

        scores_data_transformed[
            names(head(sort(apply(rank_mat, 1, quantile, 0.25)), 100)),
            .(gene = gene, score_diff = rowMeans(.SD[, cts$in_i, with = FALSE]) - rowMeans(.SD[, cts$not_i, with = FALSE]))
        ][order(-score_diff)][1:20, gene] # 1:15

    }
)

# Final heatmap with hierarchical clustering:
htmp_all_kmeans_distinct <- deconv_heatmap(
    scores_data_transformed[
        # unique( # We're removing the CAF genes because we have the volcano plot for that.
        #     c(
        #         names(tail(sort(apply(rank_mat, 1, quantile, 0.75)), 25)),
        #         unique(unlist(kmeans_clust_distinct_genes))
        #     )
        # ),
        unique(unlist(kmeans_clust_distinct_genes)),
        c('gene', pca_data[thresh_pass == TRUE, cancer_type]),
        with = FALSE
    ],
    order_genes_fun = 'seriate',
    order_genes_method = 'GW_average',
    order_analyses_fun = 'seriate',
    order_analyses_method = 'GW_average',
    plot_title = NULL
)

# Including VIM and the EMT TFs:
htmp_all_kmeans_distinct_emttfs <- deconv_heatmap(
    scores_data_transformed[
        unique(
            c(
                unlist(kmeans_clust_distinct_genes),
                'SNAI1',
                'SNAI2', # This should already be in there
                'TWIST1',
                'VIM',
                'ZEB1',
                'ZEB2'
            )
        ),
        c('gene', pca_data[thresh_pass == TRUE, cancer_type]),
        with = FALSE
    ],
    order_genes_fun = 'seriate',
    order_genes_method = 'GW_average',
    order_analyses_fun = 'seriate',
    order_analyses_method = 'GW_average',
    plot_title = NULL
)

# Tables of cancer types in each cluster and genes scoring highly in each cluster
# relative to the others:

# Details on how to make and customise tables can be found in this vignette:
# https://cran.r-project.org/web/packages/gridExtra/vignettes/tableGrob.html

table_cancer_types_data <- setNames(
    as.data.table(
        lapply(
            1:3,
            function(i) {
                pca_data[
                    thresh_pass == TRUE & kmeans_cluster_manual == i,
                    c(cancer_type, rep('', max(table(pca_data[thresh_pass == TRUE, kmeans_cluster_manual])) - .N))
                ]
            }
        )
    ),
    c('Cluster 1\nGynaecological', 'Cluster 2\nSquamous-like', 'Cluster 3\nGastro-intestinal')
)

table_cancer_types <- do.call(
    gtable_combine,
    args = lapply(
        1:3,
        function(i) {
            # tableGrob(
            #     table_cancer_types_data[, ..i],
            #     theme = ttheme_minimal(
            #         padding = unit(c(25, 10), 'pt'),
            #         base_colour = brewer.pal(3, 'Dark2')[i]
            #     ),
            #     rows = NULL
            # )
            tableGrob(
                table_cancer_types_data[, ..i],
                theme = ttheme_minimal(
                    padding = unit(c(25, 10), 'pt'),
                    core = list(fg_params = list(col = brewer.pal(3, 'Dark2')[i])),
                    colhead = list(fg_params = list(col = 'white'), bg_params = list(fill = brewer.pal(3, 'Dark2')[i], col = 'white'))
                ),
                rows = NULL
            )
        }
    )
)

# Including intermediates:

table_cancer_types_data_all <- setNames(
    cbind(
        as.data.table(
            lapply(
                1:3,
                function(i) {
                    pca_data[
                        thresh_pass == TRUE & kmeans_cluster_manual == i,
                        c(cancer_type, rep('', max(table(pca_data[thresh_pass == TRUE, kmeans_cluster_manual])) - .N))
                    ]
                }
            )
        ),
        pca_data[
            thresh_pass == FALSE,
            c(cancer_type, rep('', max(table(pca_data[thresh_pass == TRUE, kmeans_cluster_manual])) - .N))
        ]
    ),
    c('Cluster 1\nGynaecological', 'Cluster 2\nSquamous-like', 'Cluster 3\nGastro-intestinal', 'Intermediates')
)

table_cancer_types_all <- do.call(
    gtable_combine,
    args = lapply(
        1:4,
        function(i) {
            tableGrob(
                table_cancer_types_data_all[, ..i],
                theme = ttheme_minimal(
                    padding = unit(c(25, 10), 'pt'),
                    core = list(fg_params = list(col = c(brewer.pal(3, 'Dark2'), 'darkgrey')[i])),
                    colhead = list(fg_params = list(col = 'white'), bg_params = list(fill = c(brewer.pal(3, 'Dark2'), 'darkgrey')[i], col = 'white'))
                ),
                rows = NULL
            )
        }
    )
)

table_genes_data <- setNames(
    as.data.table(kmeans_clust_distinct_genes),
    c('Cluster 1\nGynaecological', 'Cluster 2\nSquamous-like', 'Cluster 3\nGastro-intestinal')
)

# Stuff I've commented out is for if you want row numbers.
table_genes <- do.call(
    gtable_combine,
    args = lapply(
        1:3,
        function(i) {

            g <- tableGrob(
                table_genes_data[, ..i],
                theme = ttheme_default(
                    core = list(
                        # fg_params = list(col = brewer.pal(3, 'Dark2')[i]),
                        bg_params = list(fill = 'white')
                    ),
                    colhead = list(
                        fg_params = list(col = 'white'),
                        bg_params = list(fill = brewer.pal(3, 'Dark2')[i], col = 'black')
                    )
                ),
                rows = NULL
                # rows = switch((i == 1) + 1, NULL, 1:nrow(table_genes_data))
            )

            gtable_add_grob(g, grobs = grid::rectGrob(gp = grid::gpar(fill = NA)), t = 2, b = nrow(g), l = 1) # l = switch((i == 1) + 1, 1, 2)

        }
    )
)





# Write to PDFs:

# Initial heatmap of 100 most common EMT/CAF genes:

pdf('../data_and_figures/scores_heatmap_supp.pdf', width = 7.5, height = 14)

# Need to mess around with the ordering so that it's roughly the same in the clustered
# and SPIN_NH'd heatmaps.

htmp_emt_caf$heatmap + theme(
    axis.text.x = element_text(
        colour = mapvalues(
            pca_data[
                with(htmp_emt_caf, analyses[ordering_analyses]),
                .(new_clust = switch(thresh_pass + 1, 0, kmeans_cluster_manual)),
                by = cancer_type
            ]$new_clust,
            0:3,
            c('darkgrey', brewer.pal(3, 'Dark2'))
        )
    )
)

dev.off()

# PCA plots using segments to indicate directions and transformations/thresholds to
# identify ambiguous cases:

pca_segment_plot <- ggplot(pca_data, aes(PC1, PC2)) +
    geom_point(aes(colour = as.character(kmeans_cluster_manual))) +
    scale_colour_manual(name = 'k-means cluster', values = brewer.pal(3, 'Set2')) +
    geom_segment(aes(x = kclust_centre['PC1'], y = kclust_centre['PC2'], xend = PC1_end, yend = PC2_end), data = segment_endpoints) +
    theme_test() +
    theme(legend.text = element_text(size = 12), legend.title = element_text(size = 13))

pdf('../data_and_figures/scores_pca_supp_alt.pdf', width = 9, height = 15)

# The colours of the Set2 palette seem like lighter versions of the Dark2 palette, hence
# I'm using them in the below plot to make them look distinct from but similar enough to
# the final groups.

plot_grid(
    blank_plot(),
    plot_grid(
        pca_segment_plot + guides(colour = FALSE),
        get_legend(pca_segment_plot),
        nrow = 1,
        ncol = 2,
        rel_widths = c(3, 1)
    ),
    plot_grid(
        plotlist = lapply(
            1:3,
            function(i) {
                ggplot(pca_data, aes(get(paste0('PC1_', i)), get(paste0('PC2_', i)))) +
                    geom_point(aes(colour = get(paste0('PC1_', i)) >= cutoffs[i])) +
                    scale_colour_manual(values = c('lightgrey', brewer.pal(3, 'Dark2')[i])) +
                    geom_text(
                        aes(x = x, y = y, label = label),
                        data = data.frame(x = cutoffs[i], y = -3.25, label = paste('Cluster', i)),
                        hjust = -0.1,
                        vjust = 0
                    ) +
                    theme_test() +
                    geom_hline(yintercept = 0) +
                    geom_vline(xintercept = cutoffs[i], linetype = 'dashed') +
                    theme(
                        axis.text.y = switch((i == 1) + 1, element_blank(), element_text()),
                        axis.title.x = switch((i == 2) + 1, element_blank(), element_text()),
                        axis.title.y = switch((i == 1) + 1, element_blank(), element_text()),
                        plot.margin = unit(c(50, 5.5, 50, 5.5), 'pt')
                    ) +
                    labs(
                        # x = paste('PC1 transformed, cluster', i),
                        # y = paste('PC2 transformed, cluster', i)
                        x = 'PC1 transformed',
                        y = 'PC2 transformed',
                        title = paste('k-means cluster', i)
                    ) +
                    lims(x = c(-2.5, 3.25), y = c(-3.25, 2.5)) +
                    guides(colour = FALSE)
            }
        ),
        nrow = 1,
        ncol = 3,
        align = 'h',
        rel_widths = c(1.125, 1, 1)
    ),
    nrow = 3,
    ncol = 1,
    rel_heights = c(0.5, 5.5, 4),
    align = 'v',
    axis = 'l'
) %>% plot_grid(table_cancer_types_all, nrow = 2, ncol = 1, rel_heights = c(9, 3)) +
    draw_label('A', x = 0, y = 0.975, hjust = -0.1, vjust = 1.1, size = 16, fontface = 2) +
    draw_label('B', x = 0, y = 0.525, hjust = -0.1, vjust = 1.1, size = 16, fontface = 2) +
    draw_label('C', x = 0, y = 0.26, hjust = -0.1, vjust = 1.1, size = 16, fontface = 2)

dev.off()

# PCA plot identifying clusters which we define after excluding the ambiguous cases:

pdf('../data_and_figures/scores_pca_clust_alt.pdf', width = 4.5, height = 5)

ggplot(pca_data, aes(PC1, PC2)) +
    # geom_point(aes(colour = as.character(kmeans_cluster))) +
    geom_point(
        aes(
            colour = mapvalues(
                kmeans_cluster_manual,
                c(1, 2, 3),
                c('Cluster 1 - Gynaecological', 'Cluster 2 - Squamous-like', 'Cluster 3 - Gastro-intestinal')
            )
        )
    ) +
    # scale_colour_manual(values = brewer.pal(3, 'Set1')) +
    # scale_colour_manual(values = brewer.pal(5, 'Dark2')[3:5]) +
    scale_colour_manual(values = brewer.pal(3, 'Dark2')) +
    geom_point(data = pca_data[thresh_pass == FALSE], colour = 'lightgrey') +
    guides(shape = FALSE) +
    theme_test() +
    theme(legend.position = 'bottom', legend.direction = 'vertical', legend.title = element_blank())
    # labs(colour = 'Cluster')

dev.off()

# Tables of cancer types and signature genes for each cluster:
pdf('../data_and_figures/scores_pca_tables_alt.pdf', width = 10, height = 8)
ggdraw(table_cancer_types)
ggdraw(table_cancer_types_all)
ggdraw(table_genes)
dev.off()

# Final clustered heatmap of filtered list of cancer types:

# I think it's also possible to colour the dendrogram, but might be more work than it's
# worth.  Here's some info in case we decide to do this:

# https://stackoverflow.com/questions/21474388/colorize-clusters-in-dendogram-with-ggplot2

pdf(
    '../data_and_figures/scores_heatmap_alt.pdf',
    width = 6.5,
    height = 9
)

deconv_heatmap_dendro_plot(
    # htmp_all_kmeans_distinct,
    c(
        list(
            heatmap = htmp_all_kmeans_distinct$heatmap + theme(
                axis.text.x = element_text(
                    angle = 55,
                    hjust = 1,
                    colour = pca_data[
                        with(htmp_all_kmeans_distinct, analyses[ordering_analyses]),
                        mapvalues(kmeans_cluster_manual, c(1, 2, 3), brewer.pal(3, 'Dark2'))
                    ]
                )
            )
        ),
        htmp_all_kmeans_distinct[-1]
    ),
    direction = 'horizontal',
    title.position = 'top',
    title.hjust = 0.5,
    barwidth = unit(70, 'pt'),
    barheight = unit(10, 'pt'),
    rel_widths = c(9, 2),
    rel_heights = c(1, 10)
)

deconv_heatmap_dendro_plot(
    # htmp_all_kmeans_distinct,
    c(
        list(
            heatmap = htmp_all_kmeans_distinct_emttfs$heatmap + theme(
                axis.text.x = element_text(
                    angle = 55,
                    hjust = 1,
                    colour = pca_data[
                        with(htmp_all_kmeans_distinct, analyses[ordering_analyses]),
                        mapvalues(kmeans_cluster_manual, c(1, 2, 3), brewer.pal(3, 'Dark2'))
                    ]
                )
            )
        ),
        htmp_all_kmeans_distinct[-1]
    ),
    direction = 'horizontal',
    title.position = 'top',
    title.hjust = 0.5,
    barwidth = unit(70, 'pt'),
    barheight = unit(10, 'pt'),
    rel_widths = c(9, 2),
    rel_heights = c(1, 10)
)

dev.off()





# Additional volcano plot to fill in white space in figure 4...  The idea is that genes
# which occur in more lists will have a larger sample size and therefore greater chance
# of becoming significant.  The fold change is replaced by average EMT-CAF score.

# Can replace wilcox.test() with e.g. t.test() in the following (we check beforehand
# that we have a sample size greater than 1).

# plot_data <- scores_data_transformed[
#     ,
#     .(
#         ave_score = rowMeans(.SD),
#         signif_val = switch(
#             (sum(rank_mat[gene, ] != 0.5) > 1) + 1,
#             1,
#             wilcox.test(
#                 as.numeric(
#                     .SD[, rank_mat[gene, ] != 0.5, with = FALSE]
#                 )
#             )$p.value
#         )
#     ),
#     by = gene
# ]

plot_data <- scores_data_transformed[, .(ave_score = rowMeans(.SD), signif_val = t.test(as.numeric(.SD))$p.value), by = gene]

pdf('../data_and_figures/scores_volcano.pdf', width = 5, height = 5)

# In the following, I'm manually choosing which labels to nudge left and which to nudge
# right.  This is because ggrepel doesn't do it adequately by itself, but it has the
# disadvantage that ggrepel doesn't know to avoid overlaps between the labels arising
# from the two separate calls to geom_label_repel().  So, if we're not careful, the
# ones that are nudged left can overlap those that are nudged right.

ggplot(plot_data[signif_val < 1], aes(x = ave_score, y = -log10(signif_val))) +
    geom_vline(xintercept = 0, linetype = 'dashed', colour = 'grey', size = 0.25) +
    geom_point(colour = 'grey40', alpha = 0.75) +
    # geom_text_repel(aes(label = gene), data = plot_data[abs(ave_score) > 0.2 & -log10(signif_val) > 6]) +
    geom_text_repel( # Just the CAF labels
        aes(label = gene),
        data = plot_data[gene %in% c('DPT', 'ACTA2', 'BGN', 'DCN', 'COL1A1', 'COMP', 'VCAN', 'MMP2', 'COL3A1')],
        nudge_x = 0.17,
        size = 3,
        point.padding = 0.2,
        segment.size = 0.3
    ) +
    geom_text_repel( # Just the pEMT labels
        aes(label = gene),
        data = plot_data[gene %in% c('LAMC1', 'PVR', 'ITGA2', 'ITGB1', 'LAMC2', 'LAMA3', 'CALU')],
        nudge_x = 0.12,
        size = 3,
        point.padding = 0.2,
        segment.size = 0.3
    ) +
    geom_text_repel(
        aes(label = gene),
        data = plot_data[gene =='CD44'],
        nudge_x = 0.1,
        nudge_y = -1,
        size = 3,
        point.padding = 0.1,
        segment.size = 0.3
    ) +
    geom_text_repel(
        aes(label = gene),
        data = plot_data[gene %in% c('VEGFA', 'COLGALT1', 'TPM4', 'DST')],
        nudge_x = -0.15,
        size = 3,
        point.padding = 0.2,
        segment.size = 0.3
    ) +
    geom_text_repel(
        aes(label = gene),
        data = plot_data[gene %in% c('ASPN', 'COL1A2', 'LUM', 'TAGLN', 'THY1', 'PCOLCE')],
        nudge_x = -0.14,
        size = 3,
        point.padding = 0.2,
        segment.size = 0.3
    ) +
    geom_text_repel(
        aes(label = gene),
        data = plot_data[gene %in% c('FAP', 'POSTN', 'COL5A1', 'COL6A2', 'SPARC', 'COL5A2')],
        nudge_x = -0.175,
        size = 3,
        point.padding = 0.2,
        segment.size = 0.3
    ) +
    geom_point(
        data = plot_data[gene %in% c('SNAI1', 'SNAI2', 'TWIST1', 'ZEB1', 'ZEB2')],
        size = 3,
        colour = brewer.pal(4, 'Dark2')[4]
    ) +
    geom_label_repel(
        aes(label = gene),
        data = plot_data[gene == 'SNAI2'],
        nudge_x = 0.12,
        size = 3,
        point.padding = 0.3,
        segment.size = 0.3
    ) +
    geom_label_repel(
        aes(label = gene),
        data = plot_data[gene %in% c('SNAI1', 'TWIST1', 'ZEB1', 'ZEB2')],
        nudge_x = -0.1,
        size = 3,
        point.padding = 0.3,
        segment.size = 0.3
    ) +
    theme_test() +
    lims(x = c(-0.75, 0.55)) +
    labs(x = 'Average pEMT-CAF score', y = TeX('Significance (-log_{10}(p))'))

dev.off()

# Previous version, using Wilcoxon test and restricted samples:

# pdf('../data_and_figures/scores_volcano.pdf', width = 5, height = 5)

# In the following, I'm manually choosing which labels to nudge left and which to nudge
# right.  This is because ggrepel doesn't do it adequately by itself, but it has the
# disadvantage that ggrepel doesn't know to avoid overlaps between the labels arising
# from the two separate calls to geom_label_repel().  So, if we're not careful, the
# ones that are nudged left can overlap those that are nudged right.

# ggplot(plot_data[signif_val < 1], aes(x = ave_score, y = -log10(signif_val))) +
#     geom_vline(xintercept = 0, linetype = 'dashed', colour = 'grey', size = 0.25) +
#     geom_point(colour = 'darkgrey') +
#     # geom_text_repel(
#     #     aes(label = gene),
#     #     data = plot_data[abs(ave_score) > 0.2 & -log10(signif_val) > 5]
#     # ) +
#     geom_text_repel(
#         aes(label = gene),
#         data = plot_data[
#             gene %in% c(
#                 'COL1A1',
#                 'COL1A2',
#                 'ASPN',
#                 'PCOLCE',
#                 'BGN',
#                 'THY1',
#                 'CD44',
#                 'VEGFA',
#                 'PLOD2',
#                 'TPM4',
#                 'TAGLN',
#                 'SLC6A8'
#             )
#         ],
#         nudge_x = -0.15,
#         size = 3,
#         point.padding = 0.2
#     ) +
#     geom_text_repel(
#         aes(label = gene),
#         data = plot_data[
#             gene %in% c(
#                 'COL11A1',
#                 'COMP',
#                 'MMP2',
#                 'FAP',
#                 'VCAN',
#                 'ITGB1',
#                 'LAMC1',
#                 'PVR',
#                 'LAMC2',
#                 'LAMA3',
#                 'POSTN'
#             )
#         ],
#         nudge_x = 0.15,
#         size = 3,
#         point.padding = 0.2
#     ) +
#     geom_text_repel(
#         aes(label = gene),
#         data = plot_data[gene == 'COPA'],
#         nudge_y = 1,
#         size = 3,
#         point.padding = 0.2
#     ) +
#     geom_text_repel(
#         aes(label = gene),
#         data = plot_data[gene == 'ITGA2'],
#         nudge_x = 0.1,
#         nudge_y = -0.8,
#         size = 3,
#         point.padding = 0.2
#     ) +
#     geom_point(
#         data = plot_data[gene %in% c('SNAI1', 'SNAI2', 'TWIST1', 'ZEB1', 'ZEB2')],
#         size = 3,
#         colour = brewer.pal(4, 'Dark2')[4]
#     ) +
#     geom_label_repel(
#         aes(label = gene),
#         data = plot_data[gene %in% c('SNAI1', 'TWIST1', 'ZEB1', 'ZEB2')],
#         nudge_x = -0.1,
#         size = 3,
#         point.padding = 0.3
#     ) +
#     geom_label_repel(
#         aes(label = gene),
#         data = plot_data[gene =='SNAI2'],
#         nudge_x = 0.1,
#         size = 3,
#         point.padding = 0.3
#     ) +
#     theme_test() +
#     lims(x = c(-0.75, 0.6)) +
#     labs(x = 'Average pEMT-CAF score', y = latex2exp::TeX('Significance (-log_{10}(p))'))

dev.off()

# Alternative for fig. 4: heatmap of 100 most common pEMT genes.  I ditched this
# because it doesn't look great and a volcano plot is a bit different (plus
# significance is a useful measure that is otherwise lacking from this part of
# the analysis).

# htmp_emt_100 <- deconv_heatmap(
#     scores_data_transformed[
#         names(head(sort(apply(rank_mat, 1, quantile, 0.25)), 100)),
#         c('gene', ..deconv_names)
#     ],
#     order_genes_fun = 'seriate',
#     order_genes_method = 'SPIN_NH',
#     order_analyses_fun = 'seriate',
#     order_analyses_method = 'SPIN_NH',
#     plot_title = 'Common pEMT genes'
# )

# Note that using the method 'PCA_angle' works well.  This orders using the first two
# principal components (somehow) - I hadn't noticed this method before, but it's
# basically what I want!

# htmp_emt_100 <- deconv_heatmap(
#     scores_data_transformed[
#         names(head(sort(apply(rank_mat, 1, quantile, 0.25)), 100)),
#         c('gene', ..deconv_names)
#     ],
#     order_genes_fun = 'seriate',
#     order_genes_method = 'PCA_angle',
#     order_analyses_fun = 'seriate',
#     order_analyses_method = 'PCA_angle',
#     plot_title = 'Common pEMT genes'
# )

# Or just order genes by average:

# htmp_emt_100 <- heat_map(
#     scores_data_transformed[
#         names(head(sort(apply(rank_mat, 1, quantile, 0.25)), 100)),
#         as.matrix(.SD[order(rowMeans(.SD)), order(colMeans(.SD)), with = FALSE]),
#         .SDcols = deconv_names
#         ] %>% t,
#     axis_text_x = NULL,
#     axis_title_x = 'Cancer types',
#     plot_margin = c(5.5, 5.5, 5.5, 0)
# )
#
# y_axis_labels <- heat_map_labels_repel(
#     scores_data_transformed[
#         names(head(sort(apply(rank_mat, 1, quantile, 0.25)), 100)),
#         c(rep('', 80), tail(gene[order(rowMeans(.SD))], 20)),
#         .SDcols = deconv_names
#         ],
#     edge = 'right'
# )
#
# pdf('../data_and_figures/scores_high_ave_emt.pdf', width = 5, height = 6)
#
# ggarrange(
#     y_axis_labels,
#     htmp_emt_100,
#     nrow = 1,
#     ncol = 2,
#     widths = c(2, 4),
#     newpage = FALSE
# )
#
# dev.off()





# Combined figure:

pdf(
    '../data_and_figures/scores_combined_figure.pdf',
    width = 11.5,
    height = 11
)

plot_grid(
    blank_plot(),
    blank_plot(),
    plot_grid(
        ggplot(plot_data[signif_val < 1], aes(x = ave_score, y = -log10(signif_val))) +
            geom_vline(xintercept = 0, linetype = 'dashed', colour = 'grey', size = 0.25) +
            geom_point(colour = 'grey40', alpha = 0.75) +
            geom_text_repel( # Just the CAF labels
                aes(label = gene),
                data = plot_data[gene %in% c('DPT', 'ACTA2', 'BGN', 'DCN', 'COL1A1', 'COMP', 'VCAN', 'MMP2', 'COL3A1')],
                nudge_x = 0.17,
                size = 3,
                point.padding = 0.2,
                segment.size = 0.3
            ) +
            geom_text_repel( # Just the pEMT labels
                aes(label = gene),
                data = plot_data[gene %in% c('LAMC1', 'PVR', 'ITGA2', 'ITGB1', 'LAMC2', 'LAMA3', 'CALU')],
                nudge_x = 0.12,
                size = 3,
                point.padding = 0.2,
                segment.size = 0.3
            ) +
            geom_text_repel(
                aes(label = gene),
                data = plot_data[gene =='CD44'],
                nudge_x = 0.1,
                nudge_y = -1,
                size = 3,
                point.padding = 0.1,
                segment.size = 0.3
            ) +
            geom_text_repel(
                aes(label = gene),
                data = plot_data[gene %in% c('VEGFA', 'COLGALT1', 'TPM4', 'DST')],
                nudge_x = -0.15,
                size = 3,
                point.padding = 0.2,
                segment.size = 0.3
            ) +
            geom_text_repel(
                aes(label = gene),
                data = plot_data[gene %in% c('ASPN', 'COL1A2', 'LUM', 'TAGLN', 'THY1', 'PCOLCE')],
                nudge_x = -0.14,
                size = 3,
                point.padding = 0.2,
                segment.size = 0.3
            ) +
            geom_text_repel(
                aes(label = gene),
                data = plot_data[gene %in% c('FAP', 'POSTN', 'COL5A1', 'COL6A2', 'SPARC', 'COL5A2')],
                nudge_x = -0.175,
                size = 3,
                point.padding = 0.2,
                segment.size = 0.3
            ) +
            geom_point(
                data = plot_data[gene %in% c('SNAI1', 'SNAI2', 'TWIST1', 'ZEB1', 'ZEB2')],
                size = 3,
                colour = brewer.pal(4, 'Dark2')[4]
            ) +
            geom_label_repel(
                aes(label = gene),
                data = plot_data[gene == 'SNAI2'],
                nudge_x = 0.13,
                size = 3,
                point.padding = 0.3,
                segment.size = 0.3
            ) +
            geom_label_repel(
                aes(label = gene),
                data = plot_data[gene %in% c('SNAI1', 'TWIST1', 'ZEB1', 'ZEB2')],
                nudge_x = -0.1,
                size = 3,
                point.padding = 0.3,
                segment.size = 0.3
            ) +
            theme_test() +
            theme(plot.margin = unit(c(5.5, 15, 20, 5.5), 'pt')) +
            lims(x = c(-0.75, 0.55)) +
            labs(x = 'Average pEMT-CAF score', y = latex2exp::TeX('Significance (-log_{10}(p))')),
        ggplot(pca_data, aes(PC1, PC2)) +
            geom_point(
                aes(
                    colour = mapvalues(
                        kmeans_cluster_manual,
                        c(1, 2, 3),
                        c('Cluster 1 - Gynaecological', 'Cluster 2 - Squamous-like', 'Cluster 3 - Gastro-intestinal')
                    )
                )
            ) +
            scale_colour_manual(values = brewer.pal(3, 'Dark2')) +
            geom_point(data = pca_data[thresh_pass == FALSE], colour = 'lightgrey') +
            guides(shape = FALSE) +
            theme_test() +
            theme(
                legend.position = 'bottom',
                legend.direction = 'vertical',
                legend.title = element_blank(),
                legend.text = element_text(size = 11),
                plot.margin = unit(c(20, 15, 5.5, 5.5), 'pt')
            ),
        nrow = 2,
        ncol = 1,
        align = 'v'
    ),
    deconv_heatmap_dendro_plot(
        c(
            list(
                heatmap = htmp_all_kmeans_distinct$heatmap +
                    scale_fill_gradientn(
                        colours = rev(colorRampPalette(brewer.pal(11, "RdBu"))(50)),
                        limits = c(-1, 1),
                        breaks = c(-1, 0, 1),
                        oob = scales::squish
                    ) +
                    theme(
                        axis.text.x = element_text(
                            angle = 55,
                            hjust = 1,
                            colour = pca_data[
                                with(htmp_all_kmeans_distinct, analyses[ordering_analyses]),
                                mapvalues(kmeans_cluster_manual, c(1, 2, 3), brewer.pal(3, 'Dark2'))
                            ]
                        ),
                        legend.text = element_text(size = 9),
                        legend.title = element_text(size = 9),
                        legend.key.width = unit(5, 'pt'),
                        legend.key.height = unit(5, 'pt')
                    )
            ),
            htmp_all_kmeans_distinct[-1]
        ),
        direction = 'horizontal',
        title.position = 'top',
        title.hjust = 0.5,
        barwidth = unit(50, 'pt'),
        barheight = unit(7.5, 'pt'),
        rel_widths = c(9, 1.25),
        rel_heights = c(0.75, 10)
    ),
    nrow = 2,
    ncol = 2,
    rel_widths = c(5, 6.5),
    rel_heights = c(1, 20)
) +
    draw_label('A', x = 0, y = 0.975, hjust = -0.1, vjust = 1.1, size = 16, fontface = 2) +
    draw_label('B', x = 0, y = 0.475, hjust = -0.1, vjust = 1.1, size = 16, fontface = 2) +
    draw_label('C', x = 0.5, y = 0.975, hjust = -0.1, vjust = 1.1, size = 16, fontface = 2)

dev.off()





# Control clustering to see whether our 3 clusters arise naturally from the expression
# data and don't really reflect EMT:

gene_variances_top_n <- sapply(
    deconv_names,
    function(ct) {
        expression_data[deconv_data[[ct]]$sample_ids, apply(.SD, 2, var), .SDcols = -'id'] %>% sort(decreasing = TRUE) %>% head(5000)
    },
    simplify = FALSE,
    USE.NAMES = TRUE
)

ct_cor_mat <- sapply(
    deconv_names,
    function(ct1) {
        sapply(
            deconv_names,
            function(ct2) {
                top_variable_genes <- intersect(names(gene_variances_top_n[[ct1]]), names(gene_variances_top_n[[ct2]]))
                expression_data[
                    ,
                    lapply(c(ct1, ct2), function(ct) {
						.SD[deconv_data[[ct]]$sample_ids, colMeans(.SD), .SDcols = top_variable_genes]) %>% setNames(c('ct1means', 'ct2means'))
					}
                ][, cor(ct1means, ct2means)]
            },
            USE.NAMES = TRUE
        )
    },
    USE.NAMES = TRUE
)

saveRDS(ct_cor_mat, '../data_and_figures/ct_cor_mat.rds')

# ct_cor_mat <- readRDS('../data_and_figures/ct_cor_mat.rds')

# gplots::heatmap.2(1 - ct_cor_mat, distfun = function(x) x, trace = 'none')

ct_cor_mat_clust <- hclust(as.dist(1 - ct_cor_mat), method = 'average')

# Or:

# ct_cor_mat_clust <- seriate(as.dist(1 - ct_cor_mat), method = 'GW_average')[[1]]

ct_cor_htmp <- heat_map(
    ct_cor_mat,
    ordering = ct_cor_mat_clust$order,
    colour_limits = c(min(ct_cor_mat), 1),
    colours = grDevices::heat.colors(50),
    axis_text_size = 11
)

# The following uses the deconv_heatmap_dendro_plot() function to plot the heatmap with
# dendrograms.  I should really rename this function, because it doesn't really need
# the output of the deconv_heatmap() function: it only needs a heatmap, clustering
# objects and dendrograms.

pdf(
    '../data_and_figures/control_for_pEMT_groups.pdf',
    width = 10,
    height = 10
)

deconv_heatmap_dendro_plot(
    list(
        heatmap = ct_cor_htmp,
        hclust_genes = ct_cor_mat_clust,
        hclust_analyses = ct_cor_mat_clust,
        dendro_genes = dendro(ct_cor_mat_clust, edge = 'left'),
        dendro_analyses = dendro(ct_cor_mat_clust, edge = 'bottom')
    ),
    direction = 'horizontal',
    title.position = 'top',
    title.hjust = 0.5,
    barwidth = unit(70, 'pt'),
    barheight = unit(10, 'pt')
)

dev.off()

# The following shows that we don't get really low numbers of variable genes shared
# between cancer types (the lowest is actually 2510):

# gplots::heatmap.2(
#     sapply(
#         deconv_names,
#         function(ct1) {
#             sapply(
#                 deconv_names,
#                 function(ct2) {
#                     length(
#                         intersect(
#                             names(gene_variances_top_n[[ct1]]),
#                             names(gene_variances_top_n[[ct2]])
#                         )
#                     )
#                 }
#             )
#         }
#     ),
#     trace = 'none'
# )

# I considered an alternative method to more closely mimic the deconvolution, but I'm
# not sure exactly how it would work.  I wanted to take, for each cancer type, two
# random sets of 20 genes and calculate the correlation of each gene (from other
# random set, or maybe the set of all genes/most variable genes) with one minus that
# with the other, mimicking my EMT-CAF scores.  I thought about first taking a random
# set of ~250 genes, ranking by correlation with purity and taking the head and tail,
# but I think this would be too close to my deconvolution and not a true control,
# since you might pick up some EMT-related signal.  But one problem is getting a set
# of genes that is common to all cancer types, as in the ~70 genes I have in my
# EMT-CAF scores heatmap.  I don't want to pick one random set, because this would be
# very specific.  But if I average over lots of gene sets, I would probably lose all
# signal.  The same goes for the random sets of ~20 genes for each cancer type.  I
# think the above clustering of averages of variable genes is probably the most
# effective and least biased control.





# Correlation with clinical features:

clinical_data <- fread('../../TCGA_data/tcga_clinical_data.csv', key = 'id')

emt_types <- list(
    deconv_emt = deconv_names,
    deconv_caf = deconv_names
)

clin_cor_genes <- setNames(
    lapply(
        list(head, tail),
        function(FUN) {
            sapply(
                deconv_data[deconv_names],
                function(deconv_ct) {
                    FUN(deconv_ct$genes_filtered[deconv_ct$ordering], 20)
                },
                simplify = FALSE,
                USE.NAMES = TRUE
            )
        }
    ),
    c('deconv_emt', 'deconv_caf')
)

# If I use the binning procedure from scrabble::score() in the following, it seems
# as though the significance values all become negative - their strengths and their
# distributions in the heatmaps are basically the same, just negative.

clin_cor <- sapply(
    names(emt_types),
    function(emt_type) {

        clinical_test(
            expression_data,
            lapply(deconv_data[emt_types[[emt_type]]], `[[`, 'sample_ids'),
            clin_cor_genes[[emt_type]],
            clinical_data,
            clin_var = c(
                'number_of_lymphnodes_positive_by_he',
                'followup_treatment_success',
                'days_to_death',
                'lymphovascular_invasion',
                'pathologic_t',
                'pathologic_n',
                'pathologic_m',
                'neoplasm_histologic_grade'
            ),
            wilcox_test_x_expr = list(
                quote(variable > 1),
                quote(variable != 'complete remission/response'),
                quote(variable < quantile(variable, 0.4)),
                quote(variable %in% c('yes', 'present')),
                quote(startsWith(variable, 't3') | startsWith(variable, 't4')),
                # Maybe n1 should also be included in the following:
				quote(startsWith(variable, 'n2') | startsWith(variable, 'n3')),
                quote(startsWith(variable, 'm1')),
                quote(variable %in% c('gb', 'g1', 'g2', 'low grade'))
            ),
            wilcox_test_y_expr = list(
                quote(variable == 0),
                quote(variable == 'complete remission/response'),
                quote(variable > quantile(variable, 0.6)),
                quote(variable %in% c('no', 'absent')),
                quote(startsWith(variable, 't0') | startsWith(variable, 't1') | startsWith(variable, 't2')),
                quote(startsWith(variable, 'n0') | startsWith(variable, 'n1')),
                quote(startsWith(variable, 'm0') | startsWith(variable, 'cm0')),
                quote(variable %in% c('g3', 'g4', 'high grade'))
            ),
            min_samples = 10
        )[
            ,
            c('nice_variable_name', 'sigval', 'sigval_adj') := .(
                mapvalues(
                    variable_name,
                    c(
                        'number_of_lymphnodes_positive_by_he',
                        'followup_treatment_success',
                        'days_to_death',
                        'lymphovascular_invasion',
                        'pathologic_t',
                        'pathologic_n',
                        'pathologic_m',
                        'neoplasm_histologic_grade'
                    ),
                    c(
                        'Lymph node metastasis',
                        'Therapy resistance',
                        'Reduced survival',
                        'Lymphovascular invasion',
                        'T stage',
                        'N stage',
                        'M stage',
                        'Grade'
                    )
                ),
                -log10(pval)*sign(fold_change),
                -log10(pval_adj)*sign(fold_change)
            )
        ]

    },
    simplify = FALSE,
    USE.NAMES = TRUE
)





# Heatmaps of clinical test significance values:

# Get expressions from these commands:
# latex2exp::TeX('\\textbf{Significance of association}')
# latex2exp::TeX('-log_{10}(p) $\\times$ sign(fold change)')

# Except I deleted a space from the following.  Not sure why I didn't need to do this
# for the scatterplots.

sig_lab <- expression(
    atop(
        `\textbf{Significance of association}` = paste(
            "",
            bold(paste("Significance of association"))
        ),
        `-log_{10}(p) $\times$ sign(fold change)` = paste(
            "-log",
            phantom()[{paste("10")}],
            "(",
            "",
            "p",
            ")",
            " ",
            "",
            phantom() %*% phantom(),
            "sign",
            # " sign",
            "(",
            "fold change",
            ")",
            ""
        )
    )
)

clin_cor_heatmaps <- sapply(
    names(clin_cor),
    function(emt_type) {

        # We'll ignore M stage, because there are no significant cases.

        clinical_test_heatmap(
            clin_cor[[emt_type]][variable_name != 'pathologic_m'],
            x_var = 'test_name',
            y_var = 'nice_variable_name',
            fill_var = 'sigval',
            hclust_method = NULL,
            x_factor_levels = htmp_emt_caf$analyses[htmp_emt_caf$ordering_analyses], # Same ordering as in the commonality heatmap
            y_factor_levels = c('Lymph node metastasis', 'N stage', 'Lymphovascular invasion', 'Grade', 'T stage',
								'Reduced survival', 'Therapy resistance'),
            colours = c('#276419', '#4D9221', '#7FBC41', '#B8E186', '#E6F5D0', '#F7F7F7', '#F7F7F7', '#F7F7F7', '#F7F7F7', '#F7F7F7',
						'#FDE0EF', '#F1B6DA', '#DE77AE', '#C51B7D', '#8E0152'), # PiYG palette with fatter waist...
            limits = c(-3, 3),
            breaks = c(-3, -2, -1, 0, 1, 2, 3),
            labels = c('-3' = '\u2264 -3', '-2' = '-2', '-1' = '-1', '0' = '0', '1' = '1', '2' = '2', '3' = '\u2265 3'),
            x_lab = NULL,
            y_lab = NULL,
            legend_title = sig_lab,
            # legend_title = NULL,
            plot_title = mapvalues(
                emt_type,
                names(clin_cor),
                c('pEMT signature', 'CAF signature'),
                warn_missing = FALSE
            ),
            axis.ticks.length = unit(0, 'pt'),
            legend.key.height = unit(15, 'pt'),
            legend.key.width = unit(7, 'pt'),
            legend.title = element_text(size = 9),
            legend.text = element_text(size = 8)
        )

    },
    simplify = FALSE,
    USE.NAMES = TRUE
)

# Save to PDF:

cairo_pdf(
    '../data_and_figures/clinical_heatmaps.pdf',
    width = 9.5,
    height = 5
)

plot_grid(
    plot_grid(
        clin_cor_heatmaps$deconv_emt + theme(axis.text.x = element_blank(), legend.position = 'none'),
        clin_cor_heatmaps$deconv_caf +
            theme(axis.text.x = element_blank(), legend.position = 'none', plot.margin = unit(c(5.5, 5.5, 1, 5.5), 'pt')),
        ggdraw(
            get_x_axis(
                clin_cor_heatmaps$deconv_caf + theme(
                    axis.text.x = element_text(
                        colour = mapvalues(
                            pca_data[
                                with(htmp_emt_caf, analyses[ordering_analyses]),
                                .(new_clust = switch(thresh_pass + 1, 0, kmeans_cluster_manual)),
                                by = cancer_type
                            ]$new_clust,
                            c(0, 1, 2, 3),
                            c('darkgrey', brewer.pal(3, 'Dark2'))
                        )
                    )
                )
            )
        ),
        # ggdraw(get_x_axis(clin_cor_heatmaps$deconv_caf)),
        nrow = 3,
        ncol = 1,
        align = 'v',
        rel_heights = c(2, 2, 1.5)
    ),
    plot_grid(
        get_legend(clin_cor_heatmaps$deconv_emt + guides(fill = guide_colourbar(title.position = 'right'))),
        nrow = 2,
        ncol = 1,
        rel_heights = c(4, 1)
    ),
    nrow = 1,
    ncol = 2,
    rel_widths = c(9.5, 3)
)

dev.off()





# Scatterplots of significance of association for EMT vs. that for CAFs:

# latex2exp::TeX('-log_{10}(p_{CAF}) $\\times$ sign(fold change)')

sig_lab_caf <- expression(
    atop(
        `\textbf{Significance of association for CAFs}` = paste(
            "",
            bold(paste("Significance of association for CAFs"))
        ),
        `-log_{10}(p_{CAF}) $\times$ sign(fold change)` = paste(
            "-log",
            phantom()[{paste("10")}],
            "(",
            "",
            "p",
            phantom()[{paste("CAF")}],
            ")",
            "",
            " ",
            "",
            phantom() %*% phantom(),
            " sign",
            "(",
            "fold change",
            ")",
            ""
        )
    )
)

# latex2exp::TeX('-log_{10}(p_{pEMT}) $\\times$ sign(fold change)')

sig_lab_emt <- expression(
    atop(
        `\textbf{Significance of association for pEMT}` = paste(
            "",
            bold(paste("Significance of association for pEMT"))
        ),
        `-log_{10}(p_{pEMT}) $\times$ sign(fold change)` = paste(
            "-log",
            phantom()[{paste("10")}],
            "(",
            "",
            "p",
            phantom()[{paste("pEMT")}],
            ")",
            "",
            " ",
            "",
            phantom() %*% phantom(),
            " sign",
            "(",
            "fold change",
            ")",
            ""
        )
    )
)

# Scatterplot including all clinical features and cancer types:

emt_caf_sig_data <- merge(
    clin_cor$deconv_emt,
    clin_cor$deconv_caf,
    by = c('test_name', 'nice_variable_name')
)[
    ,
    .(
        test_name = test_name,
        variable_name = nice_variable_name,
        sig_emt = sigval.x,
        sig_caf = sigval.y,
        sig_adj_emt = sigval_adj.x, # Include adjusted values in case I want them
        sig_adj_caf = sigval_adj.y # (I don't use them here)
    )
]

pval_adj_threshold <- adjust_threshold_bh(
    c(
        clin_cor$deconv_emt[variable_name != 'pathologic_m', pval],
        clin_cor$deconv_caf[variable_name != 'pathologic_m', pval]
    )
)

pdf(
    '../data_and_figures/clinical_scatterplot_all_features.pdf',
    width = 7,
    height = 4.5
)

ggplot(
    emt_caf_sig_data[variable_name != 'M stage'],
    aes(x = sig_caf, y = sig_emt)
) +
    geom_point(
        aes(colour = variable_name),
        shape = 17,
        size = 2.5
    ) +
    geom_text_repel(
        aes(label = str_extract(test_name, '^[A-Z]+')),
        data = emt_caf_sig_data[
            variable_name != 'M stage' & (
                sig_emt > -log10(pval_adj_threshold) |
                    sig_emt < log10(pval_adj_threshold) |
                    sig_caf > -log10(pval_adj_threshold) |
                    sig_caf < log10(pval_adj_threshold)
            )
        ],
        point.padding = 0.1,
        nudge_x = 0.1
    ) +
    scale_colour_manual(
        values = setNames(
            brewer.pal(8, "Set1")[-6],
            c(
                'Lymph node metastasis',
                'Therapy resistance',
                'Reduced survival',
                'Lymphovascular invasion',
                'T stage',
                'N stage',
                'Grade'
            )
        )
    ) +
    labs(
        x = sig_lab_caf,
        y = sig_lab_emt,
        colour = 'Clinical feature'
    ) +
    theme_test()

dev.off()

# Below is a figure combining four plots showing my "top 4" clinical features, which
# are chosen because they have sufficiently many points and sufficiently many
# significant cases.  Survival and T stage show no significant cases after adjusting
# p values, while LVI shows one but has relatively few points anyway.

scatterplots <- sapply(
    c('Grade', 'N stage', 'Reduced survival', 'Therapy resistance'),
    function(clin_feat) {

        pval_adj_threshold <- adjust_threshold_bh(
            c(
                clin_cor$deconv_emt[nice_variable_name == clin_feat, pval],
                clin_cor$deconv_caf[nice_variable_name == clin_feat, pval]
            )
        )

        ggplot(
            emt_caf_sig_data[variable_name == clin_feat][
                ,
                sig := switch(
                    (abs(sig_emt) > -log10(pval_adj_threshold) | abs(sig_caf) > -log10(pval_adj_threshold)) + 1,
                    'not_significant',
                    'significant'
                ),
                by = test_name
            ],
            aes(x = sig_caf, y = sig_emt)
        ) +
            geom_hline(yintercept = 0, colour = 'grey', size = 0.5, linetype = 'dashed') +
            geom_vline(xintercept = 0, colour = 'grey', size = 0.5, linetype = 'dashed') +
            geom_point(aes(colour = sig), shape = 17, size = 2.5) +
            scale_colour_manual(
                values = c('not_significant' = 'dodgerblue4', 'significant' = 'darkgoldenrod1'),
                labels = c('not_significant' = 'Not significant', 'significant' = 'Significant')
            ) +
            labs(title = clin_feat, colour = NULL) +
            theme_test() +
            theme(axis.title = element_blank())

    },
    simplify = FALSE,
    USE.NAMES = TRUE
)

scatterplots$Grade <- scatterplots$Grade + geom_text_repel(
    aes(label = str_replace(test_name, ' - ', '\n')),
    data = emt_caf_sig_data[variable_name == 'Grade' & test_name == 'STAD - CIN'],
    point.padding = 0.1, size = 3.5, lineheight = 0.75, segment.size = 0.3,
    nudge_x = 0.1, nudge_y = 0.5
) + geom_text_repel(
    aes(label = str_replace(test_name, ' - ', '\n')),
    data = emt_caf_sig_data[variable_name == 'Grade' & test_name == 'ESCA - Adenocarcinoma'],
    point.padding = 0.1, size = 3.5, lineheight = 0.75, segment.size = 0.3,
    nudge_x = 0.1, nudge_y = 1.2
) + geom_text_repel(
    aes(label = str_replace(test_name, ' - ', '\n')),
    data = emt_caf_sig_data[variable_name == 'Grade' & test_name %in% c('HNSC - Atypical', 'HNSC - Malignant-Basal')],
    point.padding = 0.1, size = 3.5, lineheight = 0.75, segment.size = 0.3,
    nudge_x = -0.5, nudge_y = -1
) + geom_text_repel(
    aes(label = str_replace(test_name, ' - ', '\n')),
    data = emt_caf_sig_data[variable_name == 'Grade' & test_name == 'PAAD'],
    point.padding = 0.1, size = 3.5, lineheight = 0.75, segment.size = 0.3,
    nudge_x = -0.3, nudge_y = 0.1
) + geom_text_repel(
    aes(label = str_replace(test_name, ' - ', '\n')),
    data = emt_caf_sig_data[variable_name == 'Grade' & test_name == 'CESC'],
    point.padding = 0.1, size = 3.5, lineheight = 0.75, segment.size = 0.3,
    nudge_x = 0.1
)

scatterplots$`N stage` <- scatterplots$`N stage` + geom_text_repel(
    aes(label = str_replace(test_name, ' - ', '\n')),
    data = emt_caf_sig_data[variable_name == 'N stage' & test_name == 'LUAD - Magnoid'],
    point.padding = 0.1, size = 3.5, lineheight = 0.75, segment.size = 0.3,
    nudge_x = -0.75, nudge_y = -0.1
) + geom_text_repel(
    aes(label = str_replace(test_name, ' - ', '\n')),
    data = emt_caf_sig_data[variable_name == 'N stage' & test_name == 'HNSC - Malignant-Basal'],
    point.padding = 0.1, size = 3.5, lineheight = 0.75, segment.size = 0.3,
    nudge_x = -0.5, nudge_y = -0.5
) + geom_text_repel(
    aes(label = str_replace(test_name, ' - ', '\n')),
    data = emt_caf_sig_data[variable_name == 'N stage' & test_name == 'READ'],
    point.padding = 0.1, size = 3.5, lineheight = 0.75, segment.size = 0.3,
    nudge_y = 0.5
) + geom_text_repel(
    aes(label = str_replace(test_name, ' - ', '\n')),
    data = emt_caf_sig_data[variable_name == 'N stage' & test_name == 'BRCA - Luminal A'],
    point.padding = 0.1, size = 3.5, lineheight = 0.75, segment.size = 0.3,
    nudge_x = 0.1, nudge_y = -0.5
)

scatterplots$`Reduced survival` <- scatterplots$`Reduced survival` + geom_text_repel(
    aes(label = str_replace(test_name, ' - ', '\n')),
    data = emt_caf_sig_data[variable_name == 'Reduced survival' & test_name == 'LIHC'],
    point.padding = 0.1, size = 3.5, lineheight = 0.75, segment.size = 0.3,
    nudge_x = -0.2, nudge_y = -0.1
) + geom_text_repel(
    aes(label = str_replace(test_name, ' - ', '\n')),
    data = emt_caf_sig_data[variable_name == 'Reduced survival' & test_name %in% c('LUSC - Primitive', 'BRCA - Luminal B')],
    point.padding = 0.1, size = 3.5, lineheight = 0.75, segment.size = 0.3,
    nudge_y = 0.75
)

scatterplots$`Therapy resistance` <- scatterplots$`Therapy resistance` + geom_text_repel(
    aes(label = str_replace(test_name, ' - ', '\n')),
    data = emt_caf_sig_data[variable_name == 'Therapy resistance' & test_name == 'PAAD'],
    point.padding = 0.1, size = 3.5, lineheight = 0.75, segment.size = 0.3,
    nudge_x = 0.3, nudge_y = -0.1
) + geom_text_repel(
    aes(label = str_replace(test_name, ' - ', '\n')),
    data = emt_caf_sig_data[variable_name == 'Therapy resistance' & test_name == 'HNSC - Malignant-Basal'],
    point.padding = 0.1, size = 3.5, lineheight = 0.75, segment.size = 0.3,
    nudge_x = 0.6, nudge_y = -0.1
)

pdf(
    '../data_and_figures/clinical_4_scatterplots.pdf',
    width = 8,
    height = 6.5
)

plot_grid(
    textGrob(sig_lab_emt, rot = 90, gp = gpar(fontsize = 9)),
    plot_grid(
        plotlist = lapply(
            scatterplots,
            function(g) {g + theme(legend.position = 'none')}
        ),
        nrow = 2,
        ncol = 2
    ),
    get_legend(scatterplots[[1]]),
    blank_plot(),
    textGrob(sig_lab_caf, gp = gpar(fontsize = 9)),
    blank_plot(),
    nrow = 2,
    ncol = 3,
    rel_widths = c(1, 12, 3),
    rel_heights = c(11, 1)
)

dev.off()

# plot_grid(
#     blank_plot(),
#     plot_grid(
#         plotlist = lapply(
#             c(
#                 'Grade',
#                 # 'Lymph node metastasis',
#                 # 'Lymphovascular invasion',
#                 'N stage',
#                 'Survival',
#                 # 'T stage',
#                 'Therapy resistance'
#             ),
#             function(clin_feat) {
#
#                 pval_adj_threshold <- adjust_threshold_bh(
#                     c(
#                         clin_cor$deconv_emt[nice_variable_name == clin_feat, pval],
#                         clin_cor$deconv_caf[nice_variable_name == clin_feat, pval]
#                     )
#                 )
#
#                 ggplot(
#                     emt_caf_sig_data[variable_name == clin_feat][
#                         ,
#                         sig := switch(
#                             (
#                                 abs(sig_emt) > -log10(pval_adj_threshold) |
#                                     abs(sig_caf) > -log10(pval_adj_threshold)
#                             ) + 1,
#                             'not_significant',
#                             'significant'
#                         ),
#                         by = test_name
#                     ],
#                     aes(x = sig_caf, y = sig_emt)
#                 ) +
#                     geom_hline(yintercept = 0, colour = 'grey', size = 0.5, linetype = 'dashed') +
#                     geom_vline(xintercept = 0, colour = 'grey', size = 0.5, linetype = 'dashed') +
#                     geom_point(aes(colour = sig), shape = 17, size = 2.5) +
#                     scale_colour_manual(
#                         values = c(
#                             'not_significant' = 'dodgerblue4',
#                             'significant' = 'darkgoldenrod1'
#                         )
#                     ) +
#                     geom_text_repel(
#                         aes(label = str_extract(test_name, '^[A-Z]+')),
#                         data = emt_caf_sig_data[
#                             variable_name == clin_feat & (
#                                 abs(sig_emt) > 2 | abs(sig_caf) > 2
#                                 # sig_emt > -log10(pval_adj_threshold) |
#                                 #     sig_emt < log10(pval_adj_threshold) |
#                                 #     sig_caf > -log10(pval_adj_threshold) |
#                                 #     sig_caf < log10(pval_adj_threshold)
#                             )
#                         ],
#                         point.padding = 0.1,
#                         nudge_x = 0.1,
#                         nudge_y = 0.1
#                     ) +
#                     labs(title = clin_feat) +
#                     theme_test() +
#                     theme(
#                         axis.title = element_blank(),
#                         legend.position = 'none'
#                     )
#
#             }
#         ),
#         nrow = 2,
#         ncol = 2
#     ),
#     nrow = 2,
#     ncol = 2,
#     rel_widths = c(1, 10),
#     rel_heights = c(8, 1)
# ) +
#     cowplot::draw_label(
#         sig_lab_caf,
#         x = 0.5,
#         y = 0,
#         vjust = -0.5,
#         size = 12
#     ) +
#     cowplot::draw_label(
#         sig_lab_emt,
#         x = 0,
#         y = 0.5,
#         vjust = 1.3,
#         angle = 90,
#         size = 12
#     )





# For rare vs shared EMT clinical analysis:

# rare_shared_emt <- readRDS('../data_and_figures/rare_shared_emt.rds')
#
# rare_shared_emt_analysis_names <- names(deconv_data)[
#     grepl('^HNSC|^LUSC|^LUAD|^PAAD', names(deconv_data))
# ]

# emt_types <- list(
#     rare_emt = rare_shared_emt_analysis_names,
#     shared_emt = rare_shared_emt_analysis_names
# )
#
# clin_cor_genes <- sapply(
#     names(emt_types),
#     function(emt_type) {
#         sapply(
#             emt_types[[emt_type]],
#             function(ct) {
#                 rare_shared_emt[[
#                     mapvalues(
#                         strtrim(ct, 4),
#                         c('Head', 'Lung', 'Panc'),
#                         c('hnsc', 'lung', 'paad'),
#                         warn_missing = FALSE
#                     )
#                 ]]$rare_shared_emt_data[[paste0(emt_type, '_genes')]]
#             },
#             simplify = FALSE,
#             USE.NAMES = TRUE
#         )
#     },
#     simplify = FALSE,
#     USE.NAMES = TRUE
# )

# Or:

# emt_types <- list(
#     rare_emt = rare_shared_emt_analysis_names,
#     shared_emt = rare_shared_emt_analysis_names,
#     deconv_emt = deconv_names,
#     deconv_caf = deconv_names
# )
#
# clin_cor_genes <- c(
#     sapply(
#         c('rare_emt', 'shared_emt'),
#         function(emt_type) {
#             sapply(
#                 emt_types[[emt_type]],
#                 function(ct) {
#                     rare_shared_emt[[
#                         mapvalues(
#                             strtrim(ct, 4),
#                             c('Head', 'Lung', 'Panc'),
#                             c('hnsc', 'lung', 'paad'),
#                             warn_missing = FALSE
#                         )
#                     ]]$rare_shared_emt_data[[paste0(emt_type, '_genes')]]
#                 },
#                 simplify = FALSE,
#                 USE.NAMES = TRUE
#             )
#         },
#         simplify = FALSE,
#         USE.NAMES = TRUE
#     ),
#     setNames(
#         lapply(
#             list(head, tail),
#             function(FUN) {
#                 sapply(
#                     deconv_data[deconv_names],
#                     function(deconv_ct) {
#                         FUN(deconv_ct$genes_filtered[deconv_ct$ordering], 20)
#                     },
#                     simplify = FALSE,
#                     USE.NAMES = TRUE
#                 )
#             }
#         ),
#         c('deconv_emt', 'deconv_caf')
#     )
# )

# After this, I think you can safely run the same as for deconv_emt and deconv_caf.





# Stuff for volcano plots:

# clin_cor_labels_logicals <- list(
#     rare_emt = clin_cor$rare_emt[
#         ,
#         test_name == 'Head & Neck Malignant-Basal' & variable_name == 'number_of_lymphnodes_positive_by_he' |
#             test_name == 'Head & Neck Malignant-Basal' & variable_name == 'pathologic_n' |
#             test_name == 'Head & Neck Malignant-Basal' & variable_name == 'followup_treatment_success'
#             # test_name == 'Lung Adenocarcinoma Squamoid' & variable_name == 'days_to_death' |
#             # test_name == 'Pancreatic' & variable_name == 'followup_treatment_success' |
#             # test_name == 'Head & Neck Malignant-Basal' & variable_name == 'neoplasm_histologic_grade' |
#             # test_name == 'Lung Squamous Basal' & variable_name == 'pathologic_n' |
#             # test_name == 'Head & Neck Atypical' & variable_name == 'days_to_death' |
#             # test_name == 'Lung Squamous Primitive' & variable_name == 'pathologic_m' |
#             # test_name == 'Head & Neck Malignant-Basal' & variable_name == 'pathologic_m'
#     ],
#     shared_emt = clin_cor$shared_emt[
#         ,
#         test_name == 'Head & Neck Malignant-Basal' & variable_name == 'number_of_lymphnodes_positive_by_he' |
#             # test_name == 'Head & Neck Malignant-Basal' & variable_name == 'pathologic_n' |
#             test_name == 'Pancreatic' & variable_name == 'followup_treatment_success' |
#             test_name == 'Lung Adenocarcinoma Magnoid' & variable_name == 'days_to_death' |
#             test_name == 'Pancreatic' & variable_name == 'pathologic_t'
#             # test_name == 'Lung Squamous Primitive' & variable_name == 'pathologic_m' |
#             # test_name == 'Head & Neck Classical' & variable_name == 'number_of_lymphnodes_positive_by_he' |
#             # test_name == 'Head & Neck Classical' & variable_name == 'pathologic_n' |
#             # test_name == 'Head & Neck Atypical' & variable_name == 'days_to_death' |
#             # test_name == 'Head & Neck Malignant-Basal' & variable_name == 'pathologic_m'
#     ],
#     deconv_emt = clin_cor$deconv_emt[
#         ,
#         test_name == 'Head & Neck Malignant-Basal' & variable_name == 'number_of_lymphnodes_positive_by_he' |
#             test_name == 'Head & Neck Malignant-Basal' & variable_name == 'pathologic_n' |
#             # test_name == 'Head & Neck Malignant-Basal' & variable_name == 'followup_treatment_success' |
#             test_name == 'Pancreatic' & variable_name == 'followup_treatment_success' |
#             test_name == 'Liver' & variable_name == 'days_to_death' |
#             test_name == 'Lung Adenocarcinoma Magnoid' & variable_name == 'pathologic_n'
#             # test_name == 'Rectum' & variable_name == 'days_to_death' |
#             # test_name == 'Cervical' & variable_name == 'neoplasm_histologic_grade' |
#             # test_name == 'Bladder Basal-Squamous' & variable_name == 'pathologic_m' |
#             # test_name == 'Head & Neck Malignant-Basal' & variable_name == 'pathologic_m' |
#             # test_name == 'Lung Squamous Primitive' & variable_name == 'pathologic_m'
#     ],
#     deconv_caf = clin_cor$deconv_caf[
#         ,
#         test_name == 'Head & Neck Malignant-Basal' & variable_name == 'number_of_lymphnodes_positive_by_he' |
#             # test_name == 'Head & Neck Malignant-Basal' & variable_name == 'pathologic_n' |
#             test_name == 'Colon' & variable_name == 'lymphovascular_invasion' |
#             test_name == 'Rectum' & variable_name == 'pathologic_n' |
#             test_name == 'Stomach EBV' & variable_name == 'followup_treatment_success' |
#             # test_name == 'Ovarian Immunoreactive' & variable_name == 'followup_treatment_success' |
#             # test_name == 'Colon' & variable_name == 'number_of_lymphnodes_positive_by_he' |
#             test_name == 'Breast Luminal B' & variable_name == 'days_to_death'
#             # test_name == 'Pancreatic' & variable_name == 'number_of_lymphnodes_positive_by_he' |
#             # test_name == 'Head & Neck Malignant-Basal' & variable_name == 'neoplasm_histologic_grade' |
#             # test_name == 'Stomach CIN' & variable_name == 'neoplasm_histologic_grade' |
#             # test_name == 'Oesophageal Adenocarcinoma' & variable_name == 'neoplasm_histologic_grade' |
#             # test_name == 'Lung Squamous Primitive' & variable_name == 'days_to_death' |
#             # test_name == 'Head & Neck Malignant-Basal' & variable_name == 'pathologic_m'
#     ]
# )
#
# clin_cor_plots <- sapply(
#     names(clin_cor),
#     function(emt_type) {
#         clinical_test_volcano(
#             test_data = clin_cor[[emt_type]],
#             rows_for_labels = clin_cor_labels_logicals[[emt_type]],
#             legend_title = 'Prognostic feature',
#             legend_labels = c(
#                 'number_of_lymphnodes_positive_by_he' = 'Lymph node metastasis',
#                 'followup_treatment_success' = 'Therapy resistance',
#                 'days_to_death' = 'Survival',
#                 'lymphovascular_invasion' = 'Lymphovascular invasion',
#                 'pathologic_t' = 'T stage',
#                 'pathologic_n' = 'N stage',
#                 'pathologic_m' = 'M stage',
#                 'neoplasm_histologic_grade' = 'Grade'
#             ),
#             legend_colours = setNames(
#                 brewer.pal(9, "Set1")[-6],
#                 c(
#                     'number_of_lymphnodes_positive_by_he',
#                     'followup_treatment_success',
#                     'days_to_death',
#                     'lymphovascular_invasion',
#                     'pathologic_t',
#                     'pathologic_n',
#                     'pathologic_m',
#                     'neoplasm_histologic_grade'
#                 )
#             ),
#             plot_title = mapvalues(
#                 emt_type,
#                 names(clin_cor),
#                 c(
#                     'Rare EMT',
#                     'Shared EMT',
#                     'EMT signature from deconvolution',
#                     'CAF signature from deconvolution'
#                 ),
#                 warn_missing = FALSE
#             ),
#             labels_fun = geom_text_repel,
#             box.padding = 1,
#             max.iter = 10000,
#             seed = 1826,
#             # force = 1.5,
#             nudge_x = -0.05,
#             nudge_y = 0.1
#         )
#     },
#     simplify = FALSE,
#     USE.NAMES = TRUE
# )
#
# # Save to PDF:
#
# pdf(
#     '../data_and_figures/clinical_analyses_volcano.pdf',
#     width = 10.5,
#     height = 8
# )
#
# plot_grid(
#     plot_grid(
#         clin_cor_plots$rare_emt + theme(
#             axis.title.x = element_blank(),
#             legend.position = 'none'
#         ),
#         clin_cor_plots$shared_emt + theme(
#             axis.title = element_blank(),
#             legend.position = 'none'
#         ),
#         clin_cor_plots$deconv_emt + theme(legend.position = 'none'),
#         clin_cor_plots$deconv_caf + theme(
#             axis.title.y = element_blank(),
#             legend.position = 'none'
#         ),
#         nrow = 2,
#         ncol = 2,
#         rel_widths = c(21, 20),
#         rel_heights = c(20, 21),
#         align = 'hv'
#     ),
#     get_legend(clin_cor_plots[[1]]),
#     nrow = 1,
#     ncol = 2,
#     rel_widths = c(4.75, 1)
# )
#
# dev.off()
#
# # Just the signatures from the deconvolution:
#
# pdf(
#     '../data_and_figures/clinical_analyses_volcano_deconv.pdf',
#     width = 10.5,
#     height = 4.5
# )
#
# plot_grid(
#     plot_grid(
#         clin_cor_plots$deconv_emt +
#             theme(legend.position = 'none') +
#             labs(title = 'EMT signature'),
#         clin_cor_plots$deconv_caf +
#             theme(
#                 axis.title.y = element_blank(),
#                 legend.position = 'none'
#             ) +
#             labs(title = 'CAF signature'),
#         nrow = 1,
#         ncol = 2,
#         rel_widths = c(21, 20),
#         align = 'h'
#     ),
#     get_legend(clin_cor_plots$deconv_emt),
#     nrow = 1,
#     ncol = 2,
#     rel_widths = c(4.75, 1)
# )
#
# dev.off()





# Correlation with copy number variation:

gistic_data <- fread('../data_and_figures/gistic_data.csv')
gene_meta <- fread('../data_and_figures/gene_meta.csv')
sample_meta <- fread('../data_and_figures/sample_meta.csv')

# The following is still going to take way too long...  Perhaps we should limit the
# tests to the high-level amplifications/deletions, i.e. just values 2 and -2.  Or,
# we could try to build in use of parLapply etc, to use multiple threads.

cnv_test <- rbindlist(
    lapply(
        deconv_names,
        function(ct) {

            # We subset the IDs and genes before running clinical_test_2(), because
            # this makes it run much faster.  To do this, we first have to match
            # the IDs from the expression data and the CNV data: the code for this
            # is the same as that used in the definition of clinical_test_2().

            all_ids <- data.table(
                expression_data_id = deconv_data[[ct]]$sample_ids,
                patient_id = apply(
                    str_split_fixed(
                        deconv_data[[ct]]$sample_ids,
                        '\\.',
                        4
                    )[, 1:3],
                    1,
                    paste,
                    collapse = '.'
                )
            )[
                ,
                gistic_data_id := sapply(
                    patient_id,
                    function(x) {

                        matches <- names(gistic_data)[
                            grep(paste0('^', x), names(gistic_data))
                        ]

                        if(length(matches) == 1) {
                            return(matches)
                        } else {
                            return(NA)
                        }

                    }
                )
            ][
                !(
                    gistic_data_id %in% names(table(gistic_data_id))[
                        table(gistic_data_id) > 1
                    ]
                ) & !is.na(gistic_data_id)
            ]

            test_data <- gistic_data[
                ,
                c('symbol', all_ids$gistic_data_id),
                with = FALSE
            ]

            setkey(test_data, symbol)

            rbindlist(
                lapply(
                    c(-2, -1, 1, 2),
                    function(cnv_score) {

                        genes_to_test <- test_data[
                            ,
                            .(enough_samples = sum(as.numeric(.SD) == cnv_score) >= 10),
                            by = symbol
                        ][
                            enough_samples == TRUE,
                            symbol
                        ]

                        # Run statistical tests:

                        # clinical_test_2(
                        #     expression_data,
                        #     setNames(list(all_ids$expression_data_id), ct),
                        #     clin_cor_genes$deconv_emt[[ct]],
                        #     clinical_data = tdt(test_data[genes_to_test]),
                        #     clin_var = genes_to_test,
                        #     wilcox_test_x_fun = function(z) {z == cnv_score},
                        #     wilcox_test_y_fun = function(z) {z != cnv_score},
                        #     min_samples = 10,
                        #     amatch_max_dist = NULL
                        # )

                        # I am calculating the fold changes in the below, even though I
                        # don't need it for the immediate analysis, because I might need
                        # it later!  It takes more than 50% longer when we calculate the
                        # fold change, however.

                        if(length(genes_to_test) > 0) {
                            return(
                                cbind(
                                    clinical_test_2_concise(
                                        expression_data,
                                        setNames(list(all_ids$expression_data_id), ct),
                                        clin_cor_genes$deconv_emt[[ct]],
                                        clinical_data = tdt(test_data[genes_to_test]),
                                        clin_var = genes_to_test,
                                        wilcox_test_x_fun = function(z) {z == cnv_score},
                                        wilcox_test_y_fun = function(z) {z != cnv_score}
                                        # calculate_fold_change = FALSE
                                    ),
                                    cnv_score = cnv_score
                                )
                            )
                        } else {
                            return(NULL)
                        }

                    }
                )
            )

        }
    )
)

fwrite(cnv_test, '../data_and_figures/cnv_test.csv')





cnv_test <- fread('../data_and_figures/cnv_test.csv', key = 'variable_name')

setnames(cnv_test, c('test_name', 'variable_name'), c('cancer_type', 'symbol'))

# setkey(cnv_test, symbol)

source('../../chr.order.R')

library(biomaRt)

ordered_genes <- chr.order(gene_meta$symbol)$cna_genes

plot_data <- cnv_test[
    ,
    .SD[
        ordered_genes,
        .(
            cnv_score, # For this we need to manually include cnv_score in .SDcols
            symbol,
            pval,
            pval_adj
        )
    ],
    by = .(cancer_type, abs_cnv_score = abs(cnv_score)),
    .SDcols = c('symbol', 'cnv_score', 'pval', 'pval_adj')
][
    ,
    c('sig_val', 'sig_val_adj') := .(
        switch(is.na(pval) + 1, -log10(pval)*sign(cnv_score), 0),
        switch(is.na(pval_adj) + 1, -log10(pval_adj)*sign(cnv_score), 0)
    ),
    by = .(symbol, cancer_type, abs_cnv_score, cnv_score)
]

# The following takes ages to plot if you don't remove the x axis text.  Removing the
# x axis ticks also speeds it up quite a bit.

plot_grid(
    plotlist = lapply(
        0:3,
        function(i) {
            ggplot(
                plot_data[cancer_type %in% unique(cancer_type)[(i*8 + 1):(i*8 + 8)]],
                aes(
                    factor(symbol, levels = unique(symbol)),
                    sig_val_adj,
                    group = abs_cnv_score,
                    colour = as.character(abs_cnv_score)
                )
            ) +
                facet_grid(cancer_type ~ .) +
                geom_line() +
                theme_test() +
                theme(
                    axis.text.x = element_blank(),
                    axis.ticks.x = element_blank(),
                    legend.position = 'none'
                )
        }
    ),
    nrow = 1,
    ncol = 4
)

# Also try it after including all 1.1 million p values in the adjustment:

cnv_test[, pval_adj_all := p.adjust(pval, method = 'BH')]

plot_data <- cnv_test[
    ,
    .SD[
        ordered_genes,
        .(
            cnv_score, # For this we need to manually include cnv_score in .SDcols
            symbol,
            pval,
            pval_adj_all
        )
        ],
    by = .(cancer_type, abs_cnv_score = abs(cnv_score)),
    .SDcols = c('symbol', 'cnv_score', 'pval', 'pval_adj_all')
][
    ,
    c('sig_val', 'sig_val_adj_all') := .(
        switch(is.na(pval) + 1, -log10(pval)*sign(cnv_score), 0),
        switch(is.na(pval_adj_all) + 1, -log10(pval_adj_all)*sign(cnv_score), 0)
    ),
    by = .(symbol, cancer_type, abs_cnv_score, cnv_score)
]

pdf('../data_and_figures/cnv_test.pdf', width = 20, height = 14)

plot_grid(
    plotlist = lapply(
        0:3,
        function(i) {
            ggplot(
                plot_data[cancer_type %in% unique(cancer_type)[(i*8 + 1):(i*8 + 8)]],
                aes(
                    factor(symbol, levels = unique(symbol)),
                    sig_val_adj_all,
                    # sig_val,
                    group = abs_cnv_score,
                    colour = as.character(abs_cnv_score)
                )
            ) +
                facet_grid(cancer_type ~ ., scales = 'free_y') +
                geom_line() +
                geom_hline(
                    yintercept = -log10(0.05)*c(1, -1),
                    # yintercept = -log10(adjust_threshold_bh(cnv_test$pval))*c(1, -1),
                    linetype = 'dashed'
                ) +
                theme_test() +
                theme(
                    axis.text.x = element_blank(),
                    axis.ticks.x = element_blank(),
                    legend.position = 'none'
                ) +
                labs(x = 'genes', y = NULL)
        }
    ),
    nrow = 1,
    ncol = 4
)

dev.off()
