library(data.table)
library(ggplot2)
library(cowplot)
library(magrittr)
library(limma)

source('general_functions.R')
source('tcga_functions.R')
source('sc_functions.R')





expression_data <- fread('../../TCGA_data/tcga_expression_data.csv', key = 'id')
meta_data <- fread('../../TCGA_data/tcga_meta_data.csv', key = 'id')

expression_data <- expression_data[meta_data[sample_type != 'normal', id]]
meta_data <- meta_data[sample_type != 'normal']

emt_markers <- fread('../../emt_markers.csv')[
    ,
    gene := alias2SymbolTable(gene)
][
    source != 'GO',
    sort(unique(gene))
]





# We could select which cell types to keep and which to put under 'other' using a fraction
# of the maximum of the cell type frequencies.  I decided to select them manually instead,
# because there are certain cell types I want to show in some cases, even though thery're
# a bit rare.

single_cell_metadata <- list(
    
    hnsc = list(
        
        tcga_cancer_types = 'HNSC',
        
        read_quote = quote(
            fread('../data_and_figures/puram_hnscc_2017.csv')[
                ,
                c('cell_type', 'lymph_node', 'processed_by_maxima_enzyme') := .(
                    plyr::mapvalues(
                        cell_type,
                        c('', 'b_cell', 'dendritic', 'myocyte'),
                        rep('other', 4)
                    ),
                    NULL,
                    NULL
                )
                ]
        )
        
    ),
    
    lung = list(
        
        tcga_cancer_types = c('LUAD', 'LUSC'),
        
        read_quote = quote(
            fread('../data_and_figures/lambrechts_nsclc_2018.csv')[
                sample_type != 'normal',
                -c('cluster', 'disease', 'sample_type', 'annotation')
                ][
                    ,
                    cell_type := plyr::mapvalues(
                        cell_type,
                        c('alveolar', 'epithelial'),
                        rep('other', 2)
                    )
                    ]
        )
        
    ),
    
    paad = list(
        
        tcga_cancer_types = 'PAAD',
        
        read_quote = quote(
            fread('../data_and_figures/elyada_pdac_2019.csv')[
                patient %in% c('HT137', 'HT149', 'HN149', 'HT99', 'HN150')
                ][
                    ,
                    c('cell_type', 'cluster') := .(
                        plyr::mapvalues(
                            cell_type,
                            c(
                                'ductal',
                                'acinar',
                                'perivascular',
                                'dendritic',
                                'endothelial',
                                'mast'
                                # 'red_blood_cell' # None of these left are removing patients
                            ),
                            c(
                                'cancer',
                                rep('other', 5)
                            )
                        ),
                        NULL
                    )
                    ]
        )
        
    ),
    
    lihc = list(
        
        tcga_cancer_types = 'LIHC',
        
        # Here we remove the tumours with fewer than 50 cancer cells.  Also, there are
        # surprisingly few B and T cells, and I don't want to include "HPC-like" cells
        # as a separate type in the lineplots, so I'll put all these under 'other'.
        
        read_quote = quote(
            fread('../data_and_figures/ma_liver_2019.csv')[
                patient %in% c('C26', 'C25', 'H38', 'H37', 'C56', 'C46', 'H65', 'C66')
                ][
                    ,
                    c('cell_type', 'sample', 'disease') := .(
                        plyr::mapvalues(
                            cell_type,
                            c('b_cell', 'hpc-like', 't_cell', 'unclassified'),
                            rep('other', 4)
                        ),
                        NULL,
                        NULL
                    )
                    ]
        )
        
    )
    
)

genes_list <- readRDS('../data_and_figures/sc_genes_list.rds')

sc_cancer_fibroblast_args <- list(
    
    hnsc = list(
        seed = 8511,
        genes_filter_fun = function(x) {log2(mean(2^x - 1) + 1) >= 2.5},
        scores_filter_fun = function(x) {log2(mean(2^x - 1) + 1) >= 4}
    ),
    
    lung = list(
        seed = 2566,
        genes_filter_fun = function(x) {log2(mean(2^x - 1) + 1) >= 2.5},
        scores_filter_fun = function(x) {log2(mean(2^x - 1) + 1) >= 4}
    ),
    
    brca = list(
        seed = 3718,
        genes_filter_fun = function(x) {log2(mean(2^x - 1) + 1) >= 2},
        scores_filter_fun = function(x) {log2(mean(2^x - 1) + 1) >= 3.5}
    ),
    
    coadread = list(
        seed = 3361,
        genes_filter_fun = function(x) {log2(mean(2^x - 1) + 1) >= 2},
        scores_filter_fun = function(x) {log2(mean(2^x - 1) + 1) >= 3.5}
    ),
    
    paad = list(
        seed = 5368,
        genes_filter_fun = function(x) {log2(mean(2^x - 1) + 1) >= 2.5},
        scores_filter_fun = function(x) {log2(mean(2^x - 1) + 1) >= 4}
    ),
    
    lihc = list(
        seed = 4376,
        genes_filter_fun = function(x) {log2(mean(2^x - 1) + 1) >= 2},
        scores_filter_fun = function(x) {log2(mean(2^x - 1) + 1) >= 3.5}
    ),
    
    tnbc = list(
        seed = 456,
        genes_filter_fun = function(x) {log2(mean(2^x - 1) + 1) >= 2.5},
        scores_filter_fun = function(x) {log2(mean(2^x - 1) + 1) >= 4}
    )
    
)

sc_cancer_fibroblast <- readRDS('../data_and_figures/sc_cancer_fibroblast.rds')

sc_cancer_fibroblast_heatmaps_args <- list(
    
    hnsc = list(
        annotations = c(
            'ACTA2',
            # 'AXL',
            'COL1A1',
            'COL1A2',
            'FN1',
            'JUN',
            'LAMC2',
            'MMP2',
            'SNAI1',
            'SNAI2',
            'TNC',
            'TWIST1',
            'VIM',
            'ZEB1',
            'ZEB2'
        ),
        annotations_title = 'Head and Neck'
    ),
    
    lung = list(
        annotations = c(
            'AREG',
            'COL1A1',
            'COL1A2',
            'FN1',
            'ITGB1',
            'LAMC2',
            'LGALS',
            'SNAI1',
            'SNAI2',
            'TPM4',
            'TWIST1',
            'VIM',
            'ZEB1',
            'ZEB2'
        ),
        annotations_side = 'right',
        annotations_title = 'Lung'
    ),
    
    paad = list(
        annotations = c(
            'CCL2',
            'COL1A1',
            'COL1A2',
            'CTGF',
            'FN1',
            'ITGB1',
            # 'JUN',
            'SNAI1',
            'SNAI2',
            'SPP1',
            'TIMP1',
            'TWIST1',
            'VIM',
            'ZEB1',
            'ZEB2'
        ),
        annotations_title = 'Pancreatic'
    ),
    
    lihc = list(
        annotations = c(
            'ACTA2',
            'COL1A2',
            'FN1',
            'IL32',
            'ITGB1',
            'LAMC2',
            'PLOD2',
            'SNAI1',
            'SNAI2',
            'TWIST1',
            'VIM',
            'VIM',
            'ZEB1',
            'ZEB2'
        ),
        annotations_side = 'right',
        annotations_title = 'Liver'
    )
    
)





# Single cell heatmaps:

set.seed(sc_cancer_fibroblast_args$hnsc$seed)

sc_figures_hnsc <- sc_groups_heatmap(
    sc_groups_list = sc_cancer_fibroblast$hnsc,
    groups = c('cancer', 'fibroblast'),
    title = 'Head and Neck cancer',
    default_figure_widths = list(
        annotations = 1.5,
        cancer = 6,
        fibroblast = 1.2
    ),
    figure_spacing = 10,
    annotations = sc_cancer_fibroblast_heatmaps_args$hnsc$annotations,
    annotations_text_size = 25,
    annotations_segment_size = 3,
    annotations_nudge = 0.25,
    annotations_force = 2,
    annotations_box_padding = 2,
    es_line_size = 4,
    h_legend_title = 'Expression\nlevel',
    h_legend_width = 60,
    h_legend_height = 100,
    h_legend_direction = 'vertical',
    h_legend_title_position = 'top',
    h_legend_just = 'left',
    h_legend_text_size = 70,
    h_legend_title_size = 80,
    h_legend_box_margin = margin(l = 40),
    gd_legend_title = 'Genes\ndetected',
    gd_legend_width = 60,
    gd_legend_height = 60,
    gd_legend_direction = 'vertical',
    gd_legend_title_position = 'top',
    gd_legend_just = 'left',
    gd_legend_text_size = 70,
    gd_legend_title_size = 80,
    gd_legend_box_margin = margin(l = 40),
    ticks.linewidth = 3
)

sc_figures_hnsc$plots$genes_detected$cancer <- sc_figures_hnsc$plots$genes_detected$cancer +
    theme(
        plot.title = element_text(size = 120)
    )

# sc_figures_hnsc$plots$expression_summary$cancer <- sc_figures_hnsc$plots$expression_summary$cancer +
#     theme(
#         panel.border = element_rect(size = 3),
#         axis.ticks.y = element_line(size = 3)
#     )
# 
# sc_figures_hnsc$plots$expression_summary <- sapply(
#     sc_figures_hnsc$plots$expression_summary,
#     function(g) {
#         g + theme(
#             panel.border = element_rect(size = 5),
#             axis.title.x = element_text(size = 16),
#             axis.text.y = element_text(size = 14)
#         )
#     },
#     simplify = FALSE,
#     USE.NAMES = TRUE
# )

sc_figures_mini <- sapply(
    c('lung', 'lihc', 'paad'),
    function(ct) {
        set.seed(sc_cancer_fibroblast_args[[ct]]$seed)
        sc_groups_heatmap(
            sc_groups_list = sc_cancer_fibroblast[[ct]],
            groups = c('cancer', 'fibroblast'),
            default_figure_widths = list(
                annotations = 0,
                cancer = 6,
                fibroblast = 1.2
            ),
            figure_spacing = 7.5,
            es_line_size = 3
        )
    },
    simplify = FALSE,
    USE.NAMES = TRUE
)

# The size I use below (in pixels) is just enough for an A2 poster.

png('../data_and_figures/sc_plots_mini.png', width = 4000, height = 3200)

# cairo_pdf('../data_and_figures/sc_plots_mini.pdf', width = 7, height = 5.5)

plot_grid(
    cowplot_sc(
        sc_figures_hnsc,
        heights = c(2.6, 20, 4),
        legend_space = 0.2,
        legend_rel_heights = c(1, 5),
        es_x_axis_title_vjust = 1.3,
        es_x_axis_title_size = 80,
        es_y_axis_text_size = 70,
        es_y_axis_ticks_size = 4,
        es_y_axis_ticks_length = 10,
        es_y_axis_title_angle = 0,
        es_y_axis_title_xpos = 0.85,
        es_y_axis_title_size = 80,
        es_panel_border_size = 4
    ),
    blank_plot(),
    plot_grid(
        plotlist = c(
            list(blank_plot()),
            lapply(
                c('lung', 'lihc', 'paad'),
                function(ct) {
                    plot_grid(
                        blank_plot(),
                        plot_grid(
                            blank_plot() +
                                theme(plot.title = element_text(size = 80)) +
                                labs(
                                    title = sc_cancer_fibroblast_heatmaps_args[[ct]]$annotations_title
                                ),
                            plot_grid(
                                plotlist = c(
                                    unlist(
                                        sc_figures_mini[[ct]]$plots[c('genes_detected', 'heatmaps')],
                                        recursive = FALSE
                                    ),
                                    lapply(
                                        sc_figures_mini[[ct]]$plots$expression_summary,
                                        function(g) {
                                            g + theme(
                                                panel.border = element_rect(size = 3),
                                                axis.text.y = element_blank(),
                                                axis.title.x = element_blank(),
                                                axis.ticks.length = unit(0, 'pt')
                                            )
                                        }
                                    )
                                ),
                                nrow = 3,
                                ncol = 2,
                                rel_widths = unlist(sc_figures_mini[[ct]]$figure_widths[c('cancer', 'fibroblast')]),
                                rel_heights = c(1, 10, 2)
                            ),
                            nrow = 2,
                            ncol = 1,
                            rel_heights = c(1.5, 10)
                        ),
                        blank_plot(),
                        nrow = 1,
                        ncol = 3,
                        rel_widths = c(1, 20, 1)
                    )
                    
                }
            ),
            list(blank_plot())
        ),
        nrow = 1,
        ncol = 5,
        rel_widths = c(2, 5, 5, 5, 2.5)
    ),
    blank_plot(),
    nrow = 4,
    ncol = 1,
    rel_heights = c(3, 0.3, 1, 0.1)
)

dev.off()





# Lineplots:

lineplots <- readRDS('../data_and_figures/simulated_bulk_lineplots.rds')

lineplot_fun <- function(
    
    plot_data,
    x_axis_title = 'Proportion of tumour',
    y_axis_title = 'Proportion of gene expression',
    legend_title = 'Cell type',
    legend_labels = c(
        'b_cell' = 'B cell',
        'cancer' = 'Cancer',
        'endothelial' = 'Endothelial',
        'fibroblast' = 'Fibroblast',
        'macrophage' = 'Macrophage',
        'mast' = 'Mast',
        'myeloid' = 'Myeloid',
        'plasma' = 'Plasma',
        't_cell' = 'T cell',
        'other' = 'Other',
        'lymphoid' = 'T & NK cells'
    ),
    legend_colours = c(
        'b_cell' = '#8DD3C7',
        'cancer' = '#FB8072',
        'endothelial' = '#BC80BD',
        'fibroblast' = '#FDB462',
        'macrophage' = '#80B1D3',
        'mast' = '#FFED6F',
        'myeloid' = '#80B1D3',
        'plasma' = '#FCCDE5',
        't_cell' = '#B3DE69',
        'other' = 'grey',
        'lymphoid' = '#B3DE69'
    ),
    plot_title = '',
    point_size = NULL,
    point_stroke = 1,
    line_size = NULL,
    grid_line_size = 0.25,
    panel_border_size = NULL
    
) {
    
    cell_types = unique(plot_data$initial_type)
    
    ggplot(
        data = plot_data,
        aes(
            x = proportion_content,
            y = mean_proportion_contrib,
            colour = initial_type
        )
    ) +
        geom_line(size = line_size) +
        geom_point(shape = 0, size = point_size, stroke = point_stroke) +
        theme(
            panel.background = element_blank(),
            panel.border = element_rect(
                fill = NA,
                colour = 'black',
                size = panel_border_size
            ),
            panel.grid.minor = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.major.y = element_line(
                colour = 'grey',
                size = grid_line_size,
                linetype = 'dotted'
            )
        ) +
        scale_y_continuous(
            limits = c(0, 1),
            breaks = c(0, 0.25, 0.5, 0.75, 1),
            labels = c('0' = '0', '0.25' = '0.25', '0.5' = '0.5', '0.75' = '0.75', '1' = '1'),
            expand = c(0.02, 0.02)
        ) +
        scale_x_continuous(
            expand = c(0.02, 0.02)
        ) +
        scale_colour_manual(
            labels = switch(
                is.null(legend_labels) + 1,
                legend_labels,
                sort(unique(initial_type))
            ),
            values = switch(
                is.null(legend_colours) + 1,
                legend_colours,
                setNames(
                    hcl(
                        h = seq(
                            360/(2*length(cell_types)),
                            360 - 360/(2*length(cell_types)),
                            length.out = length(cell_types)
                        ),
                        c = 100,
                        l = 65
                    ),
                    sort(unique(initial_type))
                )
            )
        ) +
        labs(
            x = x_axis_title,
            y = y_axis_title,
            colour = legend_title,
            title = plot_title
        )
    
}

new_lineplots <- sapply(
    lineplots,
    function(l) {
        lineplot_fun(
            l$data,
            line_size = 3,
            point_size = 8,
            point_stroke = 3,
            panel_border_size = 4,
            grid_line_size = 1.25,
            plot_title = l$lineplot$labels$title
        ) +
            theme(
                legend.position = 'none',
                plot.title = element_text(size = 100),
                axis.ticks = element_line(size = 4),
                axis.ticks.length = unit(10, 'pt')
            ) +
            labs(x = NULL, y = NULL)
    },
    simplify = FALSE,
    USE.NAMES = TRUE
)

dummy_legend_plot <- ggplot(
    data = data.table(
        x = 1:9,
        y = 1,
        f = factor(
            c(
                'b_cell',
                'cancer',
                'endothelial',
                'fibroblast',
                'myeloid',
                'mast',
                'plasma',
                'lymphoid',
                'rare'
            ),
            levels = c(
                'b_cell',
                'cancer',
                'endothelial',
                'fibroblast',
                'myeloid',
                'mast',
                'plasma',
                'lymphoid',
                'rare'
            )
        )
    )
) + 
    geom_tile(
        aes(x = x, y = y, fill = f)
    ) +
    scale_fill_manual(
        labels = c(
            'b_cell' = 'B cell',
            'cancer' = 'Cancer',
            'endothelial' = 'Endothelial',
            'fibroblast' = 'Fibroblast',
            'myeloid' = 'Macrophage',
            'mast' = 'Mast',
            'plasma' = 'Plasma',
            'lymphoid' = 'T & NK cells',
            'rare' = 'Rare cell types'
        ),
        values = c(
            'b_cell' = '#8DD3C7',
            'cancer' = '#FB8072',
            'endothelial' = '#BC80BD',
            'fibroblast' = '#FDB462',
            'myeloid' = '#80B1D3',
            'mast' = '#FFED6F',
            'plasma' = '#FCCDE5',
            'lymphoid' = '#B3DE69',
            'rare' = 'grey'
        )
    ) +
    labs(fill = 'Cell type') +
    theme(
        legend.title = element_text(size = 80),
        legend.text = element_text(size = 70, margin = margin(r = 100)),
        legend.key.size = unit(75, 'pt'),
        legend.justification = 'bottom'
    ) +
    guides(fill = guide_legend(nrow = 3, ncol = 3))

legend <- cowplot::get_legend(dummy_legend_plot)

png('../data_and_figures/simulated_bulk_lineplots_2x2.png', width = 2600, height = 3200)

plot_grid(
    plot_grid(
        new_lineplots$hnsc +
            theme(
                axis.text.x = element_text(size = 70, colour = 'white'), # Keeps vertical margins consistent
                axis.text.y = element_text(size = 70),
                plot.margin = unit(c(30, 30, 0, 120), 'pt')
            ),
        new_lineplots$lung +
            theme(
                axis.text = element_blank(),
                plot.margin = unit(c(30, 30, 0, 30), 'pt')
            ),
        nrow = 1,
        ncol = 2,
        rel_widths = c(1.175, 1), # Because left plots are squashed by larger left margin
        align = 'h'
    ),
    plot_grid(
        new_lineplots$paad +
            theme(
                axis.text.x = element_text(size = 70),
                axis.text.y = element_text(size = 70),
                plot.margin = unit(c(30, 30, 0, 120), 'pt')
            ),
        new_lineplots$lihc +
            theme(
                axis.text.x = element_text(size = 70),
                axis.text.y = element_blank(),
                plot.margin = unit(c(30, 30, 0, 30), 'pt')
            ),
        nrow = 1,
        ncol = 2,
        rel_widths = c(1.175, 1),
        align = 'h'
    ),
    legend,
    nrow = 3,
    ncol = 1,
    rel_heights = c(1, 1, 0.5)
) +
    cowplot::draw_label(
        'Proportion of tumour',
        x = 0.53,
        y = 0.16,
        vjust = -0.5,
        size = 80
    ) +
    cowplot::draw_label(
        'Proportion of gene expression',
        x = 0,
        y = 0.58,
        vjust = 1.3,
        angle = 90,
        size = 80
    )

dev.off()





# Stuff from TCGA deconvs:

deconv_data <- readRDS('../data_and_figures/deconv_data.rds')
deconv_plots <- readRDS('../data_and_figures/deconv_plots.rds')

ct_to_keep <- c(
    # 'acc',
    'blca_luminal_papillary',
    'blca_basal_squamous',
    'brca_luminal_a',
    'brca_luminal_b',
    'brca_basal_like_karaayvaz',
    'brca_her2_enriched',
    # 'brca_normal_like',
    'cesc',
    # 'chol',
    'coad',
    'esca_ac',
    'esca_escc',
    'hnsc_mesenchymal_basal',
    'hnsc_classical',
    'hnsc_atypical',
    'lihc',
    'luad_proximal_inflammatory',
    'luad_proximal_proliferative',
    'lusc_basal',
    'lusc_classical',
    'lusc_primitive',
    'lusc_secretory',
    # 'meso',
    'ov_differentiated',
    'ov_immunoreactive',
    'ov_proliferative',
    'paad',
    'read',
    'skcm_immune',
    'skcm_keratin',
    'skcm_mitf_low',
    'stad_cin',
    'stad_ebv',
    'stad_msi',
    # 'thca',
    'ucec'
    # 'uvm'
)

nice_names_for_figure <- c(
    # 'ACC',
    'BLCA - Luminal-Papillary',
    'BLCA - Basal-Squamous',
    'BRCA - Luminal A',
    'BRCA - Luminal B',
    'BRCA - Basal-like',
    'BRCA - HER2-enriched',
    # 'BRCA - Normal-like',
    'CESC',
    # 'CHOL',
    'COAD',
    'ESCA - Adenocarcinoma',
    'ESCA - Squamous',
    'HNSC - Malignant-Basal',
    'HNSC - Classical',
    'HNSC - Atypical',
    'LIHC',
    'LUAD - Squamoid',
    'LUAD - Magnoid',
    'LUSC - Basal',
    'LUSC - Classical',
    'LUSC - Primitive',
    'LUSC - Secretory',
    # 'MESO',
    'OV - Differentiated',
    'OV - Immunoreactive',
    'OV - Proliferative',
    'PAAD',
    'READ',
    'SKCM - Immune',
    'SKCM - Keratin',
    'SKCM - MITF-low',
    'STAD - CIN',
    'STAD - EBV',
    'STAD - MSI',
    # 'THCA',
    'UCEC'
    # 'UVM'
)

deconv_data <- deconv_data[ct_to_keep]
deconv_plots <- deconv_plots[ct_to_keep]

# Change to nice names:

names(deconv_data) <- plyr::mapvalues(
    names(deconv_data),
    names(deconv_data),
    nice_names_for_figure
)

names(deconv_plots) <- plyr::mapvalues(
    names(deconv_plots),
    names(deconv_plots),
    nice_names_for_figure
)





# Commonality heatmap:

rank_mat <- deconv_rank(deconv_data)

scores_data_transformed <- deconv_scores(
    expression_data,
    deconv_data,
    scale_fun = function(x) x/(3*sd(x)),
    transform_data = TRUE
)

deconv_names <- names(deconv_data)

scores_pca <- prcomp(
    t(
        scores_data_transformed[
            names(head(sort(apply(rank_mat, 1, quantile, 0.25)), 100)),
            ..deconv_names
        ]
    )
)

pca_data <- as.data.table(
    scores_pca$x[, 1:2],
    keep.rownames = 'cancer_type'
)

setkey(pca_data, cancer_type)

# k-means clustering with k = 3, just because to the eye it looks like there are 3:

set.seed(16918)

kclust <- kmeans(pca_data[, .(PC1, PC2)], 3)

pca_data[
    ,
    kmeans_cluster := kclust$cluster
]

# Make the names of the cluster match "1", "2", "3" in a consistent order:

kclust_map <- plyr::mapvalues(
    c(1, 2, 3),
    c(1, 2, 3),
    pca_data[
        c('BRCA - Luminal A', 'HNSC - Malignant-Basal', 'ESCA - Adenocarcinoma'),
        kmeans_cluster
    ]
)

pca_data[
    ,
    kmeans_cluster_manual := plyr::mapvalues(
        kmeans_cluster,
        c(1, 2, 3),
        kclust_map
    )
]

# Compute segments to draw on a scatterplot to illustrate the 3 directions in the points:

kclust_centre <- colMeans(kclust$centers)

segment_endpoints <- as.data.table(
    cbind(index = c(1, 2, 3), kclust$centers[kclust_map, ])
)[
    ,
    c('PC1_end', 'PC2_end') := as.list(
        extend_endpoint(kclust_centre, as.numeric(.SD), pca_data[, .(PC1, PC2)])
    ),
    by = index
]

# Transform data once for each of the 3 clusters, so that that cluster is aligned with
# the x axis:

for(i in 1:3) {
    
    pca_data <- cbind(
        pca_data,
        transform_segment(
            kclust_centre,
            kclust$centers[kclust_map[i], ],
            pca_data[, .(PC1, PC2)],
            suffix = paste0('_', i)
        )
    )
    
}

# Arbitrarily choose cutoffs for x axis such that by eye they appear to distinguish the
# cancer types most strongly associated with this "direction" in the PC space:

cutoffs <- c(0.3, 0.8, 0.8)

pca_data[
    ,
    thresh_pass := get(paste0('PC1_', unique(kmeans_cluster_manual))) >
        cutoffs[unique(kmeans_cluster_manual)],
    by = kmeans_cluster_manual
]

kmeans_clust_distinct_genes <- lapply(
    1:3,
    function(i) {
        
        cts <- list(
            in_i = pca_data[
                thresh_pass == TRUE & kmeans_cluster_manual == i,
                cancer_type
            ],
            not_i = pca_data[
                thresh_pass == TRUE & kmeans_cluster_manual != i,
                cancer_type
            ]
        )
        
        scores_data_transformed[
            names(head(sort(apply(rank_mat, 1, quantile, 0.25)), 100)),
            .(
                gene = gene,
                score_diff = rowMeans(.SD[, cts$in_i, with = FALSE]) -
                    rowMeans(.SD[, cts$not_i, with = FALSE])
            )
        ][
            order(-score_diff)
        ][
            1:14, # This is a weird number but gives a nice cluster layout
            gene
        ]
        
    }
)

htmp_all_kmeans_distinct <- deconv_heatmap(
    scores_data_transformed[
        unique(
            c(
                names(tail(sort(apply(rank_mat, 1, quantile, 0.75)), 25)),
                unique(unlist(kmeans_clust_distinct_genes)),
                c('SNAI1', 'SNAI2', 'TWIST1', 'VIM', 'ZEB1', 'ZEB2')
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

htmp_all_kmeans_distinct$dendro_genes <- ggplot(htmp_all_kmeans_distinct$dendro_genes$data) +
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend), size = 4) +
    coord_flip() +
    scale_x_continuous(expand = c(0, 0.5)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme(
        panel.background = element_rect(fill = NA),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.ticks.length = unit(0, 'pt'),
        axis.title = element_blank(),
        plot.margin = unit(c(0, 30, 0, 0), 'pt')
    )

htmp_all_kmeans_distinct$dendro_analyses <- ggplot(htmp_all_kmeans_distinct$dendro_analyses$data) +
    geom_segment(aes(x = x, y = y, xend = xend, yend = yend), size = 4) +
    scale_x_continuous(expand = c(0, 0.5)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme(
        panel.background = element_rect(fill = NA),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.ticks.length = unit(0, 'pt'),
        axis.title = element_blank(),
        plot.margin = unit(c(30, 0, 0, 0), 'pt')
    )

png('../data_and_figures/scores_heatmap_emt_caf.png', width = 5600, height = 8400)

deconv_heatmap_dendro_plot(
    c(
        list(
            heatmap = htmp_all_kmeans_distinct$heatmap +
                scale_fill_gradientn(
                    colours = rev(colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))(50)),
                    limits = c(-1, 1),
                    breaks = c(-1, 0, 1),
                    oob = scales::squish
                ) +
                theme(
                    axis.text.x = element_text(
                        size = 110,
                        angle = 55,
                        hjust = 1,
                        colour = pca_data[
                            with(htmp_all_kmeans_distinct, analyses[ordering_analyses]),
                            plyr::mapvalues(
                                kmeans_cluster_manual,
                                c(1, 2, 3),
                                RColorBrewer::brewer.pal(3, 'Dark2')
                            )
                        ],
                        margin = margin(t = 20)
                    ),
                    axis.text.y = element_text(size = 110, margin = margin(r = 20)),
                    plot.margin = unit(c(0, 0, 10, 30), 'pt'),
                    legend.text = element_text(size = 90),
                    legend.title = element_text(size = 110),
                    legend.key.width = unit(100, 'pt'),
                    legend.key.height = unit(80, 'pt')
                )
        ),
        htmp_all_kmeans_distinct[-1]
    ),
    direction = 'horizontal',
    title.position = 'top',
    title.hjust = 0.5,
    rel_widths = c(9, 1.25),
    rel_heights = c(0.75, 10),
    ticks.linewidth = 5
)

dev.off()





# Main TCGA deconv figures:

# deconv_plot(
#     list(deconv_plots$`BRCA - Luminal A`),
#     legends = FALSE
# )

deconv_plot_brca <- deconvolve_emt_caf_plots(
    data = deconv_data$`BRCA - Luminal A`,
    heatmap_axis_title = '', # Change heat_map() function so I can put NULL here
    heatmap_legend_title = 'Coexpression',
    heatmap_colours = rev(
        colorRampPalette(
            RColorBrewer::brewer.pal(11, "RdBu")
        )(50)
    ),
    heatmap_annotations = c(
        'CALU',
        'CD44',
        'COL1A1',
        'COL1A2',
        'FAP',
        'ITGAV',
        'ITGB3',
        'LAMC2',
        # 'NOTCH2',
        'PCOLCE',
        'PVR',
        'SNAI1',
        'SNAI2',
        'SPARC',
        'TGFBR3',
        'THY1',
        'TWIST1',
        'VIM',
        'ZEB1',
        'ZEB2'
    ),
    heatmap_annotations_nudge = 0.35,
    heatmap_annotations_text_size = 30,
    heatmap_annotations_segment_size = 3,
    heatmap_annotations_force = 2,
    heatmap_annotations_box_padding = 1,
    purity_colours = rev(
        colorRampPalette(
            RColorBrewer::brewer.pal(11, "PuOr")
        )(50)
    ),
    purity_colour_limits = c(-0.3, 0.3),
    purity_legend_breaks = c(-0.3, 0, 0.3),
    purity_legend_title = 'Correlation with purity\n',
    purity_legend_direction = 'horizontal',
    purity_axis_title = NULL,
    ccle_colours = rev(
        colorRampPalette(
            RColorBrewer::brewer.pal(11, 'PiYG')[c(1:3, 6, 9:11)]
        )(50)
    ),
    ccle_legend_title = 'Tumours vs. cell lines\n',
    ccle_legend_direction = 'horizontal',
    ccle_axis_title = NULL,
    extra_colours = colorRampPalette(
        c('turquoise4', 'turquoise', 'azure', 'gold2', 'gold4')
    )(50),
    extra_axis_title = NULL,
    extra_legend_title = 'scRNA-seq: CAF vs. cancer\n',
    extra_legend_direction = 'horizontal',
    bar_legend_width = NULL,
    bar_legend_height = NULL,
    diagnostics = FALSE
)

png('../data_and_figures/deconv_brca.png', width = 4800, height = 4000)

plot_grid(
    plot_grid(
        blank_plot() +
            labs(title = 'Breast Cancer') +
            theme(
                plot.title = element_text(size = 180),
                plot.margin = unit(c(50, 20, 20, 20), 'pt')
            ),
        deconv_plot_brca$plots$purity_bar +
            theme(legend.position = 'none', plot.margin = unit(c(5, 20, 5, 20), 'pt')),
        deconv_plot_brca$plots$ccle_bar +
            theme(legend.position = 'none', plot.margin = unit(c(5, 20, 5, 20), 'pt')),
        deconv_plot_brca$plots$extra_bar +
            theme(legend.position = 'none', plot.margin = unit(c(5, 20, 5, 20), 'pt')),
        deconv_plot_brca$plots$heatmap +
            theme(legend.position = 'none', plot.margin = unit(c(5, 20, 0, 20), 'pt')),
        deconv_plot_brca$plots$axis_labels +
            theme(plot.margin = unit(c(0, 20, 20, 20), 'pt')),
        nrow = 6,
        ncol = 1,
        rel_heights = c(1.5, 1, 1, 1, 20, 5)
    ),
    plot_grid(
        blank_plot(),
        get_legend(
            deconv_plot_brca$plots$purity_bar +
                theme(
                    legend.justification = 'left',
                    legend.text = element_text(size = 80),
                    legend.title = element_text(size = 100),
                    legend.key.height = unit(75, 'pt'),
                    legend.key.width = unit(100, 'pt'),
                    legend.box.margin = margin(l = 70)
                ) +
                guides(
                    fill = guide_colourbar(
                        ticks.linewidth = 4,
                        title.position = 'right'
                    )
                )
        ),
        get_legend(
            deconv_plot_brca$plots$ccle_bar +
                theme(
                    legend.justification = 'left',
                    legend.text = element_text(size = 80),
                    legend.title = element_text(size = 100),
                    legend.key.height = unit(75, 'pt'),
                    legend.key.width = unit(100, 'pt'),
                    legend.box.margin = margin(l = 70)
                ) +
                guides(
                    fill = guide_colourbar(
                        ticks.linewidth = 4,
                        title.position = 'right'
                    )
                )
        ),
        get_legend(
            deconv_plot_brca$plots$extra_bar +
                theme(
                    legend.justification = 'left',
                    legend.text = element_text(size = 80),
                    legend.title = element_text(size = 100),
                    legend.key.height = unit(75, 'pt'),
                    legend.key.width = unit(100, 'pt'),
                    legend.box.margin = margin(l = 70)
                ) +
                guides(
                    fill = guide_colourbar(
                        ticks.linewidth = 4,
                        title.position = 'right'
                    )
                )
        ),
        get_legend(
            deconv_plot_brca$plots$heatmap +
                theme(
                    legend.justification = c(0, 0),
                    legend.text = element_text(size = 80),
                    legend.title = element_text(size = 100),
                    legend.key.height = unit(150, 'pt'),
                    legend.key.width = unit(100, 'pt'),
                    legend.box.margin = margin(l = 70)
                ) +
                guides(fill = guide_colourbar(ticks.linewidth = 4)) +
                labs(fill = 'Correlation')
        ),
        blank_plot(),
        nrow = 6,
        ncol = 1,
        rel_heights = c(1.5, 2.5, 2.5, 2.5, 15.5, 5)
    ),
    nrow = 1,
    ncol = 2,
    rel_widths = c(3, 2)
)

dev.off()

# Mini heatmaps:

png('../data_and_figures/deconv_2x3_mini.png', width = 3000, height = 2400)

plot_grid(
    plotlist = c(
        lapply(
            c('Lung Squamous', 'Pancreatic', 'Liver'),
            function(ct) {
                blank_plot() +
                    labs(title = ct) +
                    theme(title = element_text(size = 80), plot.margin = unit(c(40, 20, 5, 20), 'pt'))
            }
        ),
        lapply(
            c('purity_bar', 'ccle_bar', 'extra_bar', 'heatmap'),
            function(plot_type) {
                lapply(
                    deconv_plots[c('LUSC - Classical', 'PAAD', 'LIHC')],
                    function(l) {
                        l$plots[[plot_type]] + 
                            theme(
                                legend.position = 'none',
                                title = element_blank(),
                                plot.margin = unit(
                                    c(
                                        5,
                                        20,
                                        switch((plot_type == 'heatmap') + 1, 5, 20),
                                        20
                                    ),
                                    'pt'
                                )
                            )
                    }
                )
            }
        ) %>% unlist(recursive = FALSE),
        lapply(
            c('Bladder', 'Oesophageal', 'Endometrial'),
            function(ct) {
                blank_plot() +
                    labs(title = ct) +
                    theme(title = element_text(size = 80), plot.margin = unit(c(40, 20, 5, 20), 'pt'))
            }
        ),
        lapply(
            c('purity_bar', 'ccle_bar', 'heatmap'),
            function(plot_type) {
                lapply(
                    deconv_plots[c('BLCA - Basal-Squamous', 'ESCA - Adenocarcinoma', 'UCEC')],
                    function(l) {
                        l$plots[[plot_type]] + 
                            theme(
                                legend.position = 'none',
                                title = element_blank(),
                                plot.margin = unit(
                                    c(
                                        5,
                                        20,
                                        switch((plot_type == 'heatmap') + 1, 5, 20),
                                        20
                                    ),
                                    'pt'
                                )
                            )
                    }
                )
            }
        ) %>% unlist(recursive = FALSE)
    ),
    nrow = 9,
    ncol = 3,
    rel_heights = c(2, 1, 1, 1, 15, 2, 1, 1, 15)
)

dev.off()





# Correlation with clinical features:

clinical_data <- fread('../../TCGA_data/tcga_clinical_data.csv', key = 'id')

emt_types <- c('deconv_emt', 'deconv_caf')

clin_cor_genes <- setNames(
    lapply(
        list(head, tail),
        function(FUN) {
            sapply(
                deconv_data,
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
    emt_types,
    function(emt_type) {
        
        clinical_test(
            expression_data,
            lapply(deconv_data, `[[`, 'sample_ids'),
            clin_cor_genes[[emt_type]],
            clinical_data,
            clin_var = c(
                'followup_treatment_success',
                'days_to_death',
                'pathologic_n'
            ),
            wilcox_test_x_expr = list(
                quote(variable != 'complete remission/response'),
                quote(variable < quantile(variable, 0.4)),
                quote( # Maybe n1 should be included in here as well...
                    startsWith(variable, 'n2') |
                        startsWith(variable, 'n3')
                )
            ),
            wilcox_test_y_expr = list(
                quote(variable == 'complete remission/response'),
                quote(variable > quantile(variable, 0.6)),
                quote(
                    startsWith(variable, 'n0') |
                        startsWith(variable, 'n1')
                )
            ),
            min_samples = 10
        )[
            ,
            c('nice_variable_name', 'sigval', 'sigval_adj') := .(
                plyr::mapvalues(
                    variable_name,
                    c(
                        'followup_treatment_success',
                        'days_to_death',
                        'pathologic_n'
                    ),
                    c(
                        'Therapy resistance',
                        'Survival',
                        'N stage'
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

scatterplots <- sapply(
    c(
        'N stage',
        'Survival',
        'Therapy resistance'
    ),
    function(clin_feat) {
        ggplot(
            emt_caf_sig_data[variable_name == clin_feat],
            aes(x = sig_caf, y = sig_emt)
        ) +
            geom_hline(yintercept = 0, colour = 'grey', size = 2, linetype = 'dashed') +
            geom_vline(xintercept = 0, colour = 'grey', size = 2, linetype = 'dashed') +
            geom_point(shape = 17, size = 12, colour = 'dodgerblue4') +
            labs(title = clin_feat) +
            lims(x = range(emt_caf_sig_data$sig_caf), y = range(emt_caf_sig_data$sig_emt)) +
            theme_test() +
            theme(
                plot.title = element_text(size = 80, margin = margin(b = 20)),
                axis.title = element_blank(),
                axis.text.x = element_text(size = 60),
                axis.text.y = switch(
                    (clin_feat == 'N stage') + 1,
                    element_blank(),
                    element_text(size = 60)
                ),
                axis.ticks = element_line(size = 4),
                axis.ticks.length = unit(10, 'pt'),
                panel.border = element_rect(size = 4),
                plot.margin = unit(c(20, 20, 20, 20), 'pt')
            )
    },
    simplify = FALSE,
    USE.NAMES = TRUE
)

scatterplots$`N stage` <- scatterplots$`N stage` +
    ggrepel::geom_text_repel(
        aes(label = stringr::str_replace(test_name, ' - ', '\n')),
        data = emt_caf_sig_data[variable_name == 'N stage' & test_name == 'HNSC - Malignant-Basal'],
        segment.size = 2,
        size = 20,
        box.padding = 1,
        nudge_x = -1.5,
        nudge_y = 0.3
    ) +
    ggrepel::geom_text_repel(
        aes(label = stringr::str_replace(test_name, ' - ', '\n')),
        data = emt_caf_sig_data[variable_name == 'N stage' & test_name == 'LUAD - Magnoid'],
        segment.size = 2,
        size = 20,
        box.padding = 1,
        nudge_x = -1,
        nudge_y = -0.1
    ) +
    ggrepel::geom_text_repel(
        aes(label = stringr::str_replace(test_name, ' - ', '\n')),
        data = emt_caf_sig_data[variable_name == 'N stage' & test_name == 'READ'],
        segment.size = 2,
        size = 20,
        box.padding = 1,
        nudge_x = -0.7,
        nudge_y = 0.2
    )

scatterplots$Survival <- scatterplots$Survival + ggrepel::geom_text_repel(
    aes(label = test_name),
    data = emt_caf_sig_data[variable_name == 'Survival' & test_name == 'LIHC'],
    segment.size = 2,
    size = 20,
    box.padding = 1,
    nudge_x = -0.5,
    nudge_y = -0.2
)

scatterplots$`Therapy resistance` <- scatterplots$`Therapy resistance` +
    ggrepel::geom_text_repel(
        aes(label = stringr::str_replace(test_name, ' - ', '\n')),
        data = emt_caf_sig_data[variable_name == 'Therapy resistance' & test_name == 'HNSC - Malignant-Basal'],
        segment.size = 2,
        size = 20,
        box.padding = 1,
        nudge_x = -1.5,
        nudge_y = -0.1
    ) +
    ggrepel::geom_text_repel(
        aes(label = test_name),
        data = emt_caf_sig_data[variable_name == 'Therapy resistance' & test_name == 'PAAD'],
        segment.size = 2,
        size = 20,
        box.padding = 1,
        nudge_x = 0.7,
        nudge_y = -0.1
    )

png('../data_and_figures/clinical_3_scatterplots.png', width = 4000, height = 1600)

plot_grid(
    blank_plot(),
    plot_grid(
        plotlist = scatterplots,
        nrow = 1,
        ncol = 3
    ),
    nrow = 2, # Second row will be blank, even though I haven't explicitly made blank plots
    ncol = 2,
    rel_widths = c(1, 21),
    rel_heights = c(9, 1)
) +
    cowplot::draw_label(
        sig_lab_caf,
        x = 0.53,
        y = 0,
        vjust = -0.2,
        size = 60
    ) +
    cowplot::draw_label(
        sig_lab_emt,
        x = 0,
        y = 0.53,
        vjust = 1.3,
        angle = 90,
        size = 60
    )

dev.off()
