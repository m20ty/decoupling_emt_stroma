library(data.table)
library(ggplot2)
library(cowplot)
library(magrittr)
library(limma)
library(grid)

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





deconv_data <- readRDS('../data_and_figures/deconv_data.rds')
deconv_plots <- readRDS('../data_and_figures/deconv_plots.rds')





# Automatic filtering of deconvs based on diagnostics, namely cell type correlations:

diagnostic_cell_type_lms <- sapply(
    deconv_data,
    function(deconv_ct) {
        with(
            deconv_ct,
            cor_with_initial_and_cell_types[
                genes_filtered
            ][
                ordering,
                sapply(
                    .SD,
                    function(ct) setNames(lm(ct ~ I(.I/.N))$coeff['I(.I/.N)'], NULL),
                    USE.NAMES = TRUE
                ),
                .SDcols = cell_types
            ]
        )
    },
    USE.NAMES = TRUE
)

lms_plot_mat <- diagnostic_cell_type_lms[
    ,
    c(
        # 'acc',
        'blca_luminal_infiltrated',
        'blca_luminal_papillary',
        'blca_luminal',
        'blca_basal_squamous',
        'blca_neuronal', # Should I be including this one?
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
        'kich',
        'kirc',
        'kirp',
        'lihc',
        'luad_proximal_inflammatory',
        'luad_proximal_proliferative',
        'luad_terminal_respiratory_unit',
        'lusc_basal',
        'lusc_classical',
        'lusc_primitive',
        'lusc_secretory',
        # 'meso',
        'ov_differentiated',
        'ov_immunoreactive',
        'ov_mesenchymal',
        'ov_proliferative',
        'paad',
        'prad',
        'read',
        'skcm_immune',
        'skcm_keratin',
        'skcm_mitf_low',
        'stad_cin',
        'stad_ebv',
        'stad_gs',
        'stad_msi',
        # 'thca',
        'ucec'
        # 'uvm'
    )
]

colnames(lms_plot_mat) <- nice_names(colnames(lms_plot_mat), TRUE)

rownames(lms_plot_mat) <- gsub(
    '_',
    ' ',
    sapply(
        rownames(lms_plot_mat),
        function(w) {
            gsub('^[A-Z]|^[a-z]', toupper(stringr::str_extract(w, '^[A-Z]|^[a-z]')), w)
        }
    )
)

# Diagnostic summary (filtering deconvs by correlation with other cell types):

# In the following, we could apply a condition on positive correlation with cell types.  E.g.
# we could use x < -0.1 | x > 0.4, meaning we want to exclude cancer types which have a strong
# positive association with some cell type, with an increase of 0.4 in correlation over all
# the genes.  But I'm not convinced a positive association is a problem.  It means we're
# separating out TME components, even if they're not fibroblasts.

pdf(
    '../data_and_figures/final_figures/S2.pdf',
    width = 10,
    height = 6
)

egg::ggarrange(
    heat_map(
        t(lms_plot_mat)[order(apply(lms_plot_mat, 2, min)), ],
        colour_limits = c(-0.4, 0.4),
        legend_breaks = c(-0.4, -0.1, 0.1, 0.4),
        colours = c(
            colorRampPalette(RColorBrewer::brewer.pal(11, 'PuOr')[1:4])(15),
            rep('white', 10),
            colorRampPalette(RColorBrewer::brewer.pal(11, 'PuOr')[8:11])(15)
        ),
        axis_text_x = '',
        axis_text_size = 11,
        legend_title = 'Regression\nslope'
    ) + theme(
        plot.margin = unit(c(5.5, 5.5, 1, 5.5), 'pt'),
        panel.border = element_rect(size = 0.5, fill = NA)
    ),
    ggplot(
        data.table(
            # index = 1:ncol(lms_plot_mat),
            cancer_type = colnames(lms_plot_mat),
            # slope = sort(apply(lms_plot_mat, 2, min))
            slope = apply(lms_plot_mat, 2, min)
        )[
            ,
            pass := switch((slope < -0.1) + 1, 'pass', 'fail'),
            by = cancer_type
            # by = index
        ]
    ) +
        geom_col(
            aes(
                # index,
                factor(cancer_type, levels = cancer_type[order(slope)]),
                slope,
                fill = pass,
                colour = NULL
            )
        ) +
        scale_fill_manual(values = c('#E78AC3', '#66C2A5')) +
        scale_x_discrete(expand = c(0, 0)) +
        geom_hline(
            yintercept = -0.1,
            colour = '#377EB8',
            linetype = 'dashed',
            size = 0.75
        ) +
        labs(
            y = 'Minimum regression slope',
            fill = NULL
        ) +
        theme_test() +
        theme(
            axis.title.x = element_blank(),
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 11),
            # axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            axis.ticks.length.x = unit(0, 'pt')
        ),
    ncol = 1,
    nrow = 2,
    newpage = FALSE
)

dev.off()

# Free up some space:

rm(diagnostic_cell_type_lms)
rm(lms_plot_mat)





# After some additional manual filtering (examining figures by eye), we narrow it down to
# the below list:

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

deconv_names <- names(deconv_data)





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
        
        # In the following I previously put acinar cells under the 'cancer' umbrella, but
        # decided to remove them for two reasons: first, they seem to express epithelial
        # markers at much lower levels than the ductal cells (maybe this is obvious to a
        # real biologist); and second, I believe it is actually the ductal cells that the
        # cancer originates from (the clue is in the name - PDAC - also see
        # https://www.pancreaticcancer.org.uk/information-and-support/facts-about-pancreatic-cancer/types-of-pancreatic-cancer/,
        # where they claim that less than 1% of pancreatic cancers are acinar cell
        # cancers).  So we can assume that the acinar cells won't be malignant.  This
        # also means I'm giving the acinar cells the same treatment as the alveolar cells
        # in lung, which makes sense because the alveolar duct is similar in form to the
        # acinar duct in pancreas.
        
        # EDIT 15/01/2020: we'll also filter out tumours with fewer than 50 ductal cells.
        
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





genes_list <- readRDS('../data_and_figures/sc_genes_list.rds')[
    c('hnsc', 'lung', 'paad', 'lihc')
]





# Parameters for heatmap data:

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

# Heatmap data:

sc_cancer_fibroblast <- readRDS('../data_and_figures/sc_cancer_fibroblast.rds')[
    c('hnsc', 'lung', 'paad', 'lihc')
]

# Parameters for making heatmaps:

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





# Make the heatmaps (using EMT score for the lineplots):

sc_cancer_fibroblast_heatmaps <- sapply(
    names(single_cell_metadata),
    function(ct) {
        set.seed(sc_cancer_fibroblast_args[[ct]]$seed)
        do.call(
            sc_groups_heatmap,
            args = c(
                list(
                    sc_groups_list = sc_cancer_fibroblast[[ct]],
                    groups = c('cancer', 'fibroblast'),
                    default_figure_widths = list(
                        annotations = 2.5,
                        cancer = 6,
                        fibroblast = 1.2
                    ),
                    figure_spacing = 2.5,
                    annotations_title_size = 16,
                    annotations_nudge = 0.25,
                    es_fun = NULL,
                    es_title = 'EMT score',
                    h_legend_title = 'Expression level',
                    h_legend_width = 20,
                    h_legend_height = 10,
                    h_legend_direction = 'horizontal',
                    # h_legend_title_position = 'right',
                    h_legend_title_position = 'top',
                    h_legend_just = 'top',
                    gd_legend_title = 'Genes detected',
                    gd_legend_width = 20,
                    gd_legend_height = 10,
                    gd_legend_direction = 'horizontal',
                    # gd_legend_title_position = 'left',
                    gd_legend_title_position = 'top',
                    # gd_legend_just = 'right',
                    gd_legend_just = 'top',
                    x_axis_titles = c('Cancer cells', 'CAFs')
                ),
                sc_cancer_fibroblast_heatmaps_args[[ct]]
            )
        )
    },
    simplify = FALSE,
    USE.NAMES = TRUE
)

# Make heatmaps showing average and 95th percentile of EMT TF (+ VIM) expression:

emt_tf_exp <- sapply(
    c('mean', '95th percentile'),
    function(fct) {
        sapply(
            c('cancer', 'fibroblast'),
            function(x) {
                sapply(
                    names(sc_cancer_fibroblast),
                    function(ct) {
                        sc_cancer_fibroblast[[ct]]$data[
                            cell_type == x,
                            apply(.SD, 2, switch((fct == 'mean') + 1, function(y) quantile(y, 0.95), mean)),
                            .SDcols = c('SNAI1', 'SNAI2', 'TWIST1', 'VIM', 'ZEB1', 'ZEB2')
                        ]
                    },
                    USE.NAMES = TRUE
                )
            },
            simplify = FALSE,
            USE.NAMES = TRUE
        )
    },
    simplify = FALSE,
    USE.NAMES = TRUE
)

emt_tf_exp_plots <- sapply(
    names(emt_tf_exp),
    function(fct) {
        sapply(
            names(emt_tf_exp[[fct]]),
            function(cell_type) {
                ggplot(
                    melt(
                        as.data.table(emt_tf_exp[[fct]][[cell_type]], keep.rownames = 'gene'),
                        id.vars = 'gene',
                        variable.name = 'cancer_type',
                        value.name = 'expression_level'
                    ),
                    aes(
                        x = plyr::mapvalues(
                            cancer_type,
                            c('hnsc', 'lung', 'paad', 'lihc'),
                            c('Head and Neck', 'Lung', 'Pancreatic', 'Liver')
                        ),
                        y = gene,
                        fill = expression_level
                    )
                ) +
                    geom_raster() +
                    scale_fill_gradientn(
                        colours = rev(colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))(50)),
                        limits = c(0, 12),
                        oob = scales::squish
                    ) +
                    scale_x_discrete(expand = c(0, 0)) +
                    scale_y_discrete(expand = c(0, 0)) +
                    theme(
                        legend.position = 'none',
                        axis.text.x = element_text(angle = 55, hjust = 1),
                        axis.text.y = switch((cell_type == 'cancer') + 1, element_blank(), element_text()),
                        axis.title = element_blank(),
                        plot.title = element_text(size = 11, vjust = 0),
                        panel.border = element_rect(colour = 'black', fill = NA)
                    ) +
                    labs(
                        title = plyr::mapvalues(
                            cell_type,
                            c('cancer', 'fibroblast'),
                            c('Cancer cells', 'CAFs'),
                            warn_missing = FALSE
                        )
                    )
            },
            simplify = FALSE,
            USE.NAMES = TRUE
        )
    },
    simplify = FALSE,
    USE.NAMES = TRUE
)

# The following was an attempt to convert an SVG made in Powerpoint to a grob, which could
# then be plotted along with the single cell heatmaps using plot_grid.  But the quality is
# not that good, and it doesn't look as editable (e.g. can't highlight text), so I think
# it's easier to make two separate PDFs and combine them.

# # Read in the gene filter scheme SVG:
# 
# rsvg::rsvg_ps(
#     '../data_and_figures/gene_filter_scheme.svg',
#     '../data_and_figures/gene_filter_scheme.ps',
#     # width = 25930, # Original dimensions from Powerpoint figure, multiplied by 1000
#     # height = 12200
#     width = 24000,
#     height = 10000
# )
# 
# # Had to do this:
# 
# # Sys.setenv(R_GSCMD = normalizePath('C:/Program Files/gs/gs9.51/bin/gswin64c.exe'))
# 
# grImport::PostScriptTrace(
#     '../data_and_figures/gene_filter_scheme.ps',
#     '../data_and_figures/gene_filter_scheme.ps.xml'
# )
# 
# gene_filter_scheme <- grImport::readPicture('../data_and_figures/gene_filter_scheme.ps.xml')
# 
# gene_filter_scheme <- grImport::pictureGrob(gene_filter_scheme)
# 
# # ggdraw(gene_filter_scheme)
# 
# cairo_pdf(
#     '../data_and_figures/gene_filter_scheme_cairo.pdf',
#     width = 12,
#     height = 5
# )
# 
# ggdraw(gene_filter_scheme)
# 
# dev.off()

# Single cell heatmaps:

cairo_pdf(
    '../data_and_figures/final_figures/1BC.pdf',
    width = 11,
    height = 11.3,
    onefile = TRUE
)

plot_grid(
    blank_plot(),
    plot_grid(
        cowplot_sc(
            sc_cancer_fibroblast_heatmaps$hnsc,
            legend_space = 0,
            heights = c(1.5, 20, 4),
            es_x_axis_title_vjust = 1.4,
            es_y_axis_title_angle = 0,
            es_y_axis_title_xpos = 0.82
        ),
        blank_plot(),
        cowplot_sc(
            sc_cancer_fibroblast_heatmaps$lung,
            legend_space = 0,
            heights = c(1.5, 20, 4),
            es_x_axis_title_vjust = 1.4,
            es_y_axis_title_angle = 0,
            es_y_axis_title_xpos = 0.18
        ),
        ncol = 3,
        nrow = 1,
        rel_widths = c(10, 0.8, 10)
    ),
    blank_plot(),
    plot_grid(
        cowplot_sc(
            sc_cancer_fibroblast_heatmaps$paad,
            legend_space = 0,
            heights = c(1.5, 20, 4),
            es_x_axis_title_vjust = 1.4,
            es_y_axis_title_angle = 0,
            es_y_axis_title_xpos = 0.82
        ),
        blank_plot(),
        cowplot_sc(
            sc_cancer_fibroblast_heatmaps$lihc,
            legend_space = 0,
            heights = c(1.5, 20, 4),
            es_x_axis_title_vjust = 1.4,
            es_y_axis_title_angle = 0,
            es_y_axis_title_xpos = 0.18
        ),
        ncol = 3,
        nrow = 1,
        rel_widths = c(10, 0.8, 10)
    ),
    blank_plot(),
    # plot_grid(
    #     sc_cancer_fibroblast_heatmaps$paad$plots$genes_detected_legend,
    #     blank_plot(),
    #     sc_cancer_fibroblast_heatmaps$paad$plots$heatmap_legend,
    #     nrow = 1,
    #     ncol = 3,
    #     rel_widths = c(5, 1, 5)
    # ),
    # # blank_plot(),
    plot_grid(
        plot_grid(
            plot_grid(
                plotlist = lapply(1:4, function(i) blank_plot()),
                nrow = 1,
                ncol = 4,
                rel_widths = c(2.5, 5, 2, 5)
            ) +
                draw_label('Average expression', x = 2.5/14.5 + 2.5/14.5, y = 0.5, size = 11) +
                draw_label('95th percentile', x = 2.5/14.5 + 9.5/14.5, y = 0.5, size = 11),
            plot_grid(
                blank_plot(),
                get_title(emt_tf_exp_plots$mean$cancer),
                blank_plot(),
                get_title(emt_tf_exp_plots$mean$fibroblast),
                blank_plot(),
                get_title(emt_tf_exp_plots$`95th percentile`$cancer),
                blank_plot(),
                get_title(emt_tf_exp_plots$`95th percentile`$fibroblast),
                get_y_axis(emt_tf_exp_plots$mean$cancer),
                emt_tf_exp_plots$mean$cancer + theme(
                    axis.text.x = element_blank(),
                    axis.text.y = element_blank(),
                    plot.title = element_blank(),
                    axis.ticks = element_blank(),
                    axis.ticks.length = unit(0, 'pt'),
                    plot.margin = unit(c(0, 0, 0, 0), 'pt')
                ),
                get_y_axis(emt_tf_exp_plots$mean$fibroblast),
                emt_tf_exp_plots$mean$fibroblast + theme(
                    axis.text.x = element_blank(),
                    axis.text.y = element_blank(),
                    plot.title = element_blank(),
                    axis.ticks = element_blank(),
                    axis.ticks.length = unit(0, 'pt'),
                    plot.margin = unit(c(0, 0, 0, 0), 'pt')
                ),
                get_y_axis(emt_tf_exp_plots$`95th percentile`$cancer),
                emt_tf_exp_plots$`95th percentile`$cancer + theme(
                    axis.text.x = element_blank(),
                    axis.text.y = element_blank(),
                    plot.title = element_blank(),
                    axis.ticks = element_blank(),
                    axis.ticks.length = unit(0, 'pt'),
                    plot.margin = unit(c(0, 0, 0, 0), 'pt')
                ),
                get_y_axis(emt_tf_exp_plots$`95th percentile`$fibroblast),
                emt_tf_exp_plots$`95th percentile`$fibroblast + theme(
                    axis.text.x = element_blank(),
                    axis.text.y = element_blank(),
                    plot.title = element_blank(),
                    axis.ticks = element_blank(),
                    axis.ticks.length = unit(0, 'pt'),
                    plot.margin = unit(c(0, 0, 0, 0), 'pt')
                ),
                blank_plot(),
                get_x_axis(emt_tf_exp_plots$mean$cancer),
                blank_plot(),
                get_x_axis(emt_tf_exp_plots$mean$fibroblast),
                blank_plot(),
                get_x_axis(emt_tf_exp_plots$`95th percentile`$cancer),
                blank_plot(),
                get_x_axis(emt_tf_exp_plots$`95th percentile`$fibroblast),
                nrow = 3,
                ncol = 8,
                rel_widths = c(2.5, 2.2, 0.4, 2.2, 2, 2.2, 0.4, 2.2),
                rel_heights = c(0.5, 5, 3)
            ),
            nrow = 2,
            ncol = 1,
            rel_heights = c(1, 9)
        ),
        plot_grid(
            blank_plot(),
            sc_cancer_fibroblast_heatmaps$paad$plots$genes_detected_legend,
            blank_plot(),
            sc_cancer_fibroblast_heatmaps$paad$plots$heatmap_legend,
            nrow = 3,
            ncol = 2,
            rel_widths = c(1, 5),
            rel_heights = c(1, 1, 1)
        ),
        nrow = 1,
        ncol = 2,
        rel_widths = c(14.1, 6.076)
    ),
    ncol = 1,
    nrow = 6,
    rel_heights = c(1, 10, 1, 10, 1.9, 10)
) +
    draw_label('B', x = 0, y = 0.99, hjust = -0.1, vjust = 1.1, size = 20, fontface = 2) +
    draw_label('C', x = 0, y = 0.31, hjust = -0.1, vjust = 1.1, size = 20, fontface = 2)

dev.off()

# Free up some space:

rm(sc_cancer_fibroblast_heatmaps)
rm(emt_tf_exp)
rm(emt_tf_exp_plots)





# Plots to show average expression of mesenchymal and epithelial genes:

# For this I need the analysis of EMT genes in just the cancer cells.  I was hoping
# not to have to do this - probably I can reengineer everything so I don't.

sc_cancer_args <- list(
    
    hnsc = list(
        seed = 8511,
        genes_filter_fun = function(x) {
            log2(mean(2^x - 1) + 1) >= 3 |
                sum(x >= 7) >= length(x)/100
        },
        scores_filter_fun = function(x) {
            log2(mean(2^x - 1) + 1) >= 4 |
                sum(x >= 7) >= length(x)/100
        }
    ),
    
    lung = list(
        seed = 2566,
        genes_filter_fun = function(x) {
            log2(mean(2^x - 1) + 1) >= 3 |
                sum(x >= 7) >= length(x)/100
        },
        scores_filter_fun = function(x) {
            log2(mean(2^x - 1) + 1) >= 4 |
                sum(x >= 7) >= length(x)/100
        }
    ),
    
    paad = list(
        seed = 5368,
        genes_filter_fun = function(x) {
            log2(mean(2^x - 1) + 1) >= 3 |
                sum(x >= 7) >= length(x)/100
        },
        scores_filter_fun = function(x) {
            log2(mean(2^x - 1) + 1) >= 4 |
                sum(x >= 7) >= length(x)/100
        }
    ),
    
    lihc = list(
        seed = 4376,
        genes_filter_fun = function(x) {
            log2(mean(2^x - 1) + 1) >= 2.5 |
                sum(x >= 7) >= length(x)/100
        },
        scores_filter_fun = function(x) {
            log2(mean(2^x - 1) + 1) >= 3.5 |
                sum(x >= 7) >= length(x)/100
        }
    )
    
)

sc_cancer <- readRDS('../data_and_figures/sc_cancer.rds')[
    c('hnsc', 'lung', 'paad', 'lihc')
]

emt_epi_comparison <- sapply(
    names(single_cell_metadata),
    function(ct) {
        
        cat(paste0(ct, '...'))
        
        sc_data <- eval(single_cell_metadata[[ct]]$read_quote)
        
        setkey(sc_data, id)
        
        all_genes_filtered <- sc_data[
            sc_cancer[[ct]]$cells_filtered,
            names(.SD)[apply(.SD, 2, sc_cancer_args[[ct]]$genes_filter_fun)],
            .SDcols = -c('id', 'patient', 'cell_type')
        ]
        
        epi_markers <- names(sc_data)[
            grep('^CDH1$|^EPCAM$|^SFN$|^KRT[0-9]', names(sc_data))
        ]
        
        epi_markers <- epi_markers[epi_markers %in% all_genes_filtered]
        
        plot_data <- sc_data[
            sc_cancer[[ct]]$cells_filtered,
            c(
                .(
                    id = id,
                    emt_score = sc_cancer[[ct]]$scores[id, score],
                    epi_score = signature_score(
                        set_colnames(t(.SD[, ..all_genes_filtered]), id),
                        epi_markers,
                        nbin = length(all_genes_filtered) %/% 110
                    )
                    # epi_score = scrabble::score(
                    #     set_colnames(t(.SD[, ..all_genes_filtered]), id),
                    #     list(epi_markers),
                    #     bin.control = TRUE,
                    #     nbin = length(all_genes_filtered) %/% 110
                    # )[, 1]
                ),
                .SD[, ..epi_markers]
            ),
            by = patient
        ][
            order(-emt_score),
            epi_score_runmean := caTools::runmean(epi_score, .N/10)
        ]
        
        lineplots <- setNames(
            lapply(
                c('emt_score', 'epi_score_runmean'),
                function(v) {
                    ggplot(
                        plot_data[order(-emt_score)],
                        aes(
                            factor(id, levels = id),
                            get(v),
                            group = 1
                        )
                    ) +
                        geom_line() +
                        theme_test() +
                        theme(
                            axis.text.x = element_blank(),
                            axis.ticks.x = element_blank()
                        ) +
                        labs(x = 'Cells', y = 'Relative average expression')
                }
            ),
            c('emt', 'epi')
        )
        
        combined_lineplot <- ggplot(
            melt(
                plot_data[, .(id, emt_score, epi_score_runmean)],
                id.vars = 'id'
            ),
            aes(
                factor(id, levels = plot_data[order(-emt_score), id]),
                value,
                group = variable,
                colour = variable
            )
        ) +
            geom_line() +
            scale_colour_discrete(
                labels = c(
                    'emt_score' = 'EMT score',
                    'epi_score_runmean' = 'Epithelial score'
                )
            ) +
            theme_test() +
            theme(
                axis.text.x = element_blank(),
                axis.ticks.x = element_blank()
            ) +
            labs(x = 'Cells', y = 'Score', colour = NULL)
        
        # epi_emt_corr <- sort(
        #     plot_data[, cor(emt_score, .SD)[1, ], .SDcols = epi_markers]
        # )
        
        # I changed the above step to calculate the correlations per tumour, then average them:
        
        epi_emt_corr <- sort(
            plot_data[
                ,
                .(gene = epi_markers, corr = cor(emt_score, .SD)[1, ]),
                .SDcols = epi_markers,
                by = patient
            ][
                ,
                .(ave_corr = mean(corr[!is.na(corr)])),
                by = gene
            ][
                ,
                setNames(ave_corr, gene)
            ]
        )
        
        epi_heatmap <- ggplot(
            melt(
                plot_data[
                    , # Z scores instead of expression levels to make correlation with EMT score clearer
                    c(
                        .(id = id),
                        sapply(
                            .SD,
                            function(x) {(x - mean(x))/sd(x)},
                            simplify = FALSE,
                            USE.NAMES = TRUE
                        )
                    ),
                    .SDcols = epi_markers
                    ],
                # plot_data[, c('id', ..epi_markers)],
                id.vars = 'id',
                variable.name = 'gene',
                value.name = 'expression_level'
            ),
            aes(
                x = factor(
                    id,
                    levels = unique(id)[sc_cancer[[ct]]$ordering_cells$cancer]
                ),
                y = factor(
                    gene,
                    # levels = plot_data[, names(sort(colMeans(.SD))), .SDcols = epi_markers]
                    levels = names(epi_emt_corr)
                ),
                fill = expression_level
            )
        ) +
            geom_raster() +
            scale_fill_gradientn(
                limits = c(-2, 2),
                # limits = c(0, 12),
                colours = rev(colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))(50)),
                # colours = rev(colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))(50)),
                oob = scales::squish,
                breaks = c(-2, -1, 0, 1, 2),
                labels = c('-2' = '\u2264 -2', '-1' = '-1', '0' = '0', '1' = '1', '2' = '\u2265 2')
                # breaks = c(0, 3, 6, 9, 12),
                # labels = c('0' = '0', '3' = '3', '6' = '6', '9' = '9', '12' = '\u2265 12')
            ) +
            scale_x_discrete(expand = c(0, 0)) +
            scale_y_discrete(expand = c(0, 0)) +
            theme(
                axis.text.x = element_blank(),
                axis.ticks.x = element_blank()
            ) +
            labs(x = 'Cells', y = 'Epithelial markers', fill = 'Expression\nlevel Z-score')
        
        epi_emt_corr_barplot <- ggplot(
            data.table(
                gene = names(epi_emt_corr),
                corr = epi_emt_corr
            ),
            aes(x = factor(gene, levels = gene), y = corr)
        ) +
            geom_col(width = 0.8) +
            scale_x_discrete(expand = c(0, 0.5)) +
            theme_test() +
            labs(x = 'Gene', y = 'Correlation with EMT score') +
            coord_flip()
        
        cat('Done!\n')
        
        list(
            lineplots = lineplots,
            combined_lineplot = combined_lineplot,
            epi_heatmap = epi_heatmap,
            # epi_emt_corr_lineplot = epi_emt_corr_lineplot,
            epi_emt_corr_barplot = epi_emt_corr_barplot,
            data = plot_data,
            epi_markers = epi_markers,
            epi_emt_corr = epi_emt_corr
        )
        
    },
    simplify = FALSE,
    USE.NAMES = TRUE
)

aligned_plots_emt_epi <- sapply(
    names(emt_epi_comparison),
    function(ct) {
        
        aligned_1 <- align_plots(
            emt_epi_comparison[[ct]]$epi_heatmap +
                theme(
                    axis.title = element_blank(),
                    plot.margin = unit(c(5.5, 1, 1, 5.5), 'pt'),
                    legend.position = 'none'
                ) +
                labs(title = sc_cancer_fibroblast_heatmaps_args[[ct]]$annotations_title),
            emt_epi_comparison[[ct]]$epi_emt_corr_barplot +
                theme(
                    axis.text.y = element_blank(),
                    axis.title.y = element_blank(),
                    axis.ticks.y = element_blank(),
                    # Setting margin is a bodge to allow us to align everything:
                    axis.title.x = element_text(margin = margin(t = 50)),
                    plot.margin = unit(c(5.5, 5.5, 5.5, 1), 'pt')
                ) +
                ylim(range(unlist(lapply(emt_epi_comparison, `[[`, 'epi_emt_corr')))) +
                labs(y = 'Correlation with\nEMT score'), # It's confusing with coord_flip()...
            align = 'h'
        )
        
        aligned_2 <- align_plots(
            emt_epi_comparison[[ct]]$epi_heatmap +
                theme(
                    axis.title = element_blank(),
                    plot.margin = unit(c(5.5, 1, 1, 5.5), 'pt'),
                    legend.position = 'none'
                ) +
                labs(title = sc_cancer_fibroblast_heatmaps_args[[ct]]$annotations_title),
            emt_epi_comparison[[ct]]$combined_lineplot +
                theme(
                    axis.title.y = element_text(margin = margin(r = -20)),
                    plot.margin = unit(c(1, 1, 5.5, 5.5), 'pt'),
                    legend.position = 'none'
                ) +
                labs(
                    x = switch(
                        (which(names(emt_epi_comparison) == ct) == length(emt_epi_comparison)) + 1,
                        NULL,
                        'Cells'
                    )
                ),
            align = 'v'
        )
        
        aligned_1[[1]]$heights[8] <- sum(unit(2.75, 'pt'), unit(0, 'cm'))
        aligned_1[[1]]$heights[9] <- unit(0, 'cm')
        aligned_1[[1]]$heights[12] <- unit(1, 'pt')
        aligned_1[[1]]$widths[3] <- aligned_2[[1]]$widths[3]
        aligned_1[[2]]$heights[8] <- sum(unit(2.75, 'pt'), unit(0, 'cm'))
        aligned_1[[2]]$heights[9] <- unit(0, 'cm')
        aligned_1[[2]]$heights[12] <- unit(1, 'pt')
        
        list(aligned_1[[1]], aligned_1[[2]], aligned_2[[2]])
        
    },
    simplify = FALSE,
    USE.NAMES = TRUE
)

# Make PDF:

cairo_pdf(
    '../data_and_figures/final_figures/S1.pdf',
    width = 10,
    height = 12
)

plot_grid(
    plot_grid(
        plotlist = lapply(
            c('hnsc', 'lung', 'paad', 'lihc'),
            # c('hnsc', 'lung', 'paad', 'lihc', 'tnbc'),
            # names(emt_epi_comparison),
            function(ct) {
                plot_grid(
                    plotlist = aligned_plots_emt_epi[[ct]],
                    nrow = 2,
                    ncol = 2,
                    rel_heights = c(
                        2 + (length(emt_epi_comparison[[ct]]$epi_markers) - 12)*0.15,
                        1
                    ),
                    rel_widths = c(5, 1)
                )
            }
        ),
        nrow = 4,
        ncol = 1,
        # align = 'v',
        rel_heights = sapply(
            emt_epi_comparison[c('hnsc', 'lung', 'paad', 'lihc')],
            # emt_epi_comparison[c('hnsc', 'lung', 'paad', 'lihc', 'tnbc')],
            function(li) {
                3 + (length(li$epi_markers) - 12)*0.075
            }
        )
    ),
    plot_grid(
        blank_plot(),
        get_legend(emt_epi_comparison[[1]]$epi_heatmap),
        get_legend(emt_epi_comparison[[1]]$combined_lineplot),
        blank_plot(),
        nrow = 4,
        ncol = 1,
        rel_heights = c(2, 1, 1, 2)
    ),
    nrow = 1,
    ncol = 2,
    rel_widths = c(10, 2)
)

dev.off()

# Free up some space:

rm(emt_epi_comparison)
rm(aligned_plots_emt_epi)





# Lineplots of contributions of different cell types to mesenchymal signal in simulated
# tumours:

max_mean_counts <- list(
    hnsc = 300,
    lung = 1000,
    paad = 500,
    lihc = 500
)

lineplots <- readRDS('../data_and_figures/simulated_bulk_lineplots.rds')[
    c('hnsc', 'lung', 'paad', 'lihc')
]

lineplots_dummy_legend_plot <- ggplot(
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
    theme(legend.key.size = unit(15, 'pt'))

# lineplots_legend <- cowplot::get_legend(lineplots_dummy_legend_plot)





# Deconvolution for simulated bulk profiles:

deconv_plot_args_per_ct <- list(
    
    hnsc = list(
        heatmap_annotations = c(
            'ACTA2',
            'COL1A2',
            'COL6A1',
            'ITGA2',
            'LAMC2',
            'SNAI1',
            'SNAI2',
            'TGFBI',
            'THY1',
            'TNC',
            'TWIST1',
            'VIM',
            'ZEB1',
            'ZEB2'
        ),
        plot_title = 'Head and Neck'
    ),
    
    lung = list(
        heatmap_annotations = c(
            'ACTA2',
            # 'AREG',
            'COL1A1',
            'COL1A2',
            'ITGB1',
            'LAMC2',
            'QSOX1',
            'RHOB',
            'SNAI1',
            'SNAI2',
            'THY1',
            'TWIST1',
            'VIM',
            'ZEB1',
            'ZEB2'
        ),
        plot_title = 'Lung'
    ),
    
    paad = list(
        heatmap_annotations = c(
            'COL1A1',
            'COL1A2',
            'COL6A2',
            'ITGB1',
            'LAMC2',
            'MMP3',
            'SNAI1',
            'SNAI2',
            'SPP1',
            'THY1',
            'TWIST1',
            'VIM',
            'ZEB1',
            'ZEB2'
        ),
        plot_title = 'Pancreatic'
    ),
    
    lihc = list(
        heatmap_annotations = c(
            'ACTA2',
            'CDH2',
            'COL1A2',
            'EFEMP2',
            'ITGA2',
            'LAMC2',
            'PLOD2',
            'SNAI1',
            'SNAI2',
            'THY1',
            'TWIST1',
            'VIM',
            'ZEB1',
            'ZEB2'
        ),
        plot_title = 'Liver'
    )
    
)

simulated_deconvs <- readRDS('../data_and_figures/simulated_deconvs.rds')[
    c('hnsc', 'lung', 'paad', 'lihc')
]

# In the following, I am not opting to include epithelial cells in the diagnostics, mainly to save
# time (it saves a lot of time), but also because the epithelial cell correlations are (thankfully)
# as I would expect, and because conceptually I don't think it adds much.

simulated_deconv_plots <- sapply(
    
    names(deconv_plot_args_per_ct),
    
    function(ct) {
        
        cat(ct, '\b...')
        
        deconv_ct <- do.call(
            deconvolve_emt_caf_plots,
            args = c(
                list(
                    data = simulated_deconvs[[ct]],
                    # Include the following only if you want epithelial scores (takes much longer):
                    # expression_data = simulated_bulk_data,
                    heatmap_legend_title = 'Correlation',
                    heatmap_colours = rev(
                        colorRampPalette(
                            RColorBrewer::brewer.pal(11, "RdBu")
                        )(50)
                    ),
                    heatmap_colour_limits = c(-1, 1),
                    heatmap_legend_breaks = c(-1, -0.5, 0, 0.5, 1),
                    heatmap_annotations_nudge = 0.3,
                    purity_colours = rev(
                        colorRampPalette(
                            RColorBrewer::brewer.pal(11, "PuOr")
                        )(50)
                    ),
                    purity_colour_limits = c(-1, 1),
                    purity_legend_breaks = c(-1, 0, 1),
                    purity_legend_title = 'Correlation with purity\n',
                    purity_legend_direction = 'horizontal',
                    purity_axis_title = NULL,
                    ccle_colours = rev(
                        colorRampPalette(
                            RColorBrewer::brewer.pal(11, 'PiYG')#[c(1:3, 6, 9:11)]
                        )(50)
                    ),
                    ccle_legend_breaks = c(-1, 0, 1),
                    ccle_legend_title = 'Tumours vs. cell lines\n',
                    ccle_legend_direction = 'horizontal',
                    ccle_axis_title = NULL,
                    extra_colours = colorRampPalette(
                        c('turquoise4', 'turquoise', 'azure', 'gold2', 'gold4')
                    )(50),
                    extra_legend_title = 'scRNA-seq: CAF vs. cancer\n',
                    extra_legend_direction = 'horizontal',
                    extra_axis_title = NULL,
                    bar_legend_justification = 'left'
                    # bar_legend_width = NULL,
                    # bar_legend_height = NULL
                ),
                deconv_plot_args_per_ct[[ct]]
            )
        )
        
        cat('Done!\n')
        
        deconv_ct
        
    },
    
    simplify = FALSE,
    USE.NAMES = TRUE
    
)

pemt_caf_bracket_params <- list(
    hnsc = list(
        left_bracket_xmax = 34,
        right_bracket_xmin = 52
    ),
    lung = list(
        left_bracket_xmax = 27,
        right_bracket_xmin = 60
    ),
    paad = list(
        left_bracket_xmax = 27,
        right_bracket_xmin = 54
    ),
    lihc = list(
        left_bracket_xmax = 27,
        right_bracket_xmin = 54
    )
)

# Main figure with aligned lineplots and simulated deconv plots:

pdf(
    '../data_and_figures/final_figures/2.pdf',
    width = 14,
    height = 8.8
)

plot_grid(
    blank_plot(),
    plot_grid(
        blank_plot(),
        plot_grid(
            plotlist = c(
                lapply(
                    c('hnsc', 'lung', 'paad', 'lihc'),
                    function(ct) {
                        plot_grid(
                            plotlist = c(
                                list(
                                    lineplots[[ct]]$lineplot +
                                        theme(
                                            legend.position = 'none',
                                            plot.margin = unit(c(5.5, 3.5, 55, 3.5), 'pt'),
                                            axis.text.y = switch((ct == 'hnsc') + 1, element_blank(), element_text()),# NULL),
                                            plot.title = element_text(size = 14.5)
                                        ) +
                                        labs(x = NULL, y = switch((ct == 'hnsc') + 1, NULL, waiver()))
                                ),
                                lapply(
                                    simulated_deconv_plots[[ct]]$plots[
                                        c('purity_bar', 'ccle_bar', 'extra_bar')
                                        ],
                                    function(g) {
                                        g + theme(
                                            legend.position = 'none',
                                            plot.margin = unit(c(0, 3.5, 2, 3.5), 'pt')
                                        )
                                    }
                                ),
                                list(
                                    simulated_deconv_plots[[ct]]$plots$heatmap +
                                        theme(
                                            legend.position = 'none',
                                            axis.title.y = element_text(vjust = -7.2),
                                            plot.margin = unit(c(0, 3.5, 0, 3.5), 'pt')
                                        ) +
                                        labs(title = NULL, y = switch((ct == 'hnsc') + 1, NULL, waiver())),
                                    simulated_deconv_plots[[ct]]$plots$axis_labels,
                                    do.call(
                                        pemt_caf_brackets,
                                        args = c(
                                            list(brackets_y = 0.7),
                                            pemt_caf_bracket_params[[ct]]
                                        )
                                    )
                                )
                            ),
                            nrow = 7,
                            ncol = 1,
                            rel_heights = c(23, 1, 1, 1, 15, 5, 2),
                            align = 'v'
                        )
                    }
                ),
                list(
                    plot_grid(
                        plotlist = c(
                            list(get_legend(lineplots_dummy_legend_plot + theme(legend.justification = 'left'))),
                            lapply(
                                c('purity_bar', 'ccle_bar', 'extra_bar', 'heatmap'),
                                function(plot_type) {
                                    get_legend(
                                        simulated_deconv_plots$hnsc$plots[[plot_type]] +
                                            theme(
                                                legend.title = element_text(size = 9),
                                                legend.justification = c(0, switch((plot_type == 'heatmap') + 1, 1, 0))
                                            ) +
                                            guides(fill = guide_colourbar(title.position = 'right')) +
                                            labs(
                                                fill = plyr::mapvalues(
                                                    plot_type,
                                                    c('purity_bar', 'ccle_bar', 'extra_bar', 'heatmap'),
                                                    c(
                                                        'Correlation\nwith purity',
                                                        'Tumours vs.\ncell lines',
                                                        'scRNA-seq:\nCAF vs. cancer',
                                                        'Correlation'
                                                    ),
                                                    warn_missing = FALSE
                                                )
                                            )
                                    )
                                }
                            ),
                            list(blank_plot())
                        ),
                        nrow = 6,
                        ncol = 1,
                        rel_heights = c(23, 3, 3, 3, 9, 7)
                    )
                )
            ),
            nrow = 1,
            ncol = 5,
            rel_widths = c(1.15, 1, 1, 1, 0.6)
        ),
        nrow = 2,
        ncol = 1,
        rel_heights = c(1, 21)
    ),
    nrow = 1,
    ncol = 2,
    rel_widths = c(1, 50)
) +
    draw_label('A', x = 0, y = 0.975, hjust = -0.1, vjust = 1.1, size = 20, fontface = 2) +
    draw_label('B', x = 0, y = 0.525, hjust = -0.1, vjust = 1.1, size = 20, fontface = 2) +
    draw_label(
        'Proportion of tumour',
        x = 0.465,
        y = 0.565,
        size = 11
    ) 

dev.off()

# Free up some space:

rm(simulated_deconvs)
rm(simulated_deconv_plots)
rm(lineplots)
rm(pemt_caf_bracket_params)





# Make plots to compare EMT and CAF genes in single cell analysis vs. TCGA bulk deconvolution:

sc_to_bulk_names <- list(
    hnsc = c(
        'HNSC - Malignant-Basal',
        'HNSC - Classical',
        'HNSC - Atypical'
    ),
    lung = c(
        'LUAD - Squamoid',
        'LUAD - Magnoid',
        'LUSC - Basal',
        'LUSC - Classical',
        'LUSC - Primitive',
        'LUSC - Secretory'
    ),
    paad = 'PAAD',
    lihc = 'LIHC'
)

# Get genes from single cell data and TCGA deconv data:

genes_sc <- unique(
    unlist(
        lapply(
            sc_cancer_fibroblast,
            `[[`,
            'genes_filtered'
        )
    )
)

genes_tcga <- unique(
    unlist(
        lapply(
            deconv_data[unlist(sc_to_bulk_names)],
            `[[`,
            'genes_filtered'
        )
    )
)

# Calculate scores from transformed TCGA bulk expression data:

scores_data_tcga <- deconv_scores(
    expression_data,
    deconv_data[unlist(sc_to_bulk_names)],
    scale_fun = function(x) x/(3*sd(x[!is.na(x)])),
    transform_data = TRUE,
    additional_genes = genes_sc[genes_sc %in% names(expression_data)]
)

scores <- sapply(
    names(single_cell_metadata),
    function(ct) {
        
        sc_data <- eval(single_cell_metadata[[ct]]$read_quote)
        
        all_genes <- unique(
            c(
                sc_cancer_fibroblast[[ct]]$genes_filtered,
                unlist(
                    lapply(
                        deconv_data[sc_to_bulk_names[[ct]]],
                        `[[`,
                        'genes_filtered'
                    )
                )
            )
        )
        
        scores_tcga <- scores_data_tcga[
            all_genes,
            setNames(rowMeans(.SD), gene),
            .SDcols = sc_to_bulk_names[[ct]]
        ]
        
        # In the following, would it be sensible to reverse the log before taking means?
        
        # I was going to take the top and bottom 20 from the sorted vector of the following
        # gene scores, then score genes by correlation with these head and tail genes.  But 
        # it's probably enough to use these scores as they are (or after log, to make them
        # normally distributed).
        
        scores_sc <- sc_data[
            cell_type == 'cancer',
            colMeans(.SD),
            .SDcols = all_genes[all_genes %in% names(sc_data)]
        ]/sc_data[
            cell_type == 'fibroblast',
            colMeans(.SD),
            .SDcols = all_genes[all_genes %in% names(sc_data)]
        ]
        
        out <- data.table(
            gene = all_genes,
            sc_score = scores_sc[all_genes],
            tcga_score = scores_tcga
        )
        
        out[
            ,
            sc_score := log10(sc_score)
        ][
            ,
            c('sc_score', 'tcga_score') := lapply(
                .SD,
                function(x) {
                    sapply(
                        x,
                        function(y) {
                            switch(
                                (!is.na(y) & !is.infinite(y)) + 1,
                                NA,
                                (y - mean(x[!is.na(x) & !is.infinite(x)]))/
                                    sd(x[!is.na(x) & !is.infinite(x)])
                            )
                        }
                    )
                }
            ),
            .SDcols = -'gene'
        ]
        
        out
        
    },
    simplify = FALSE,
    USE.NAMES = TRUE
)

# Function to make combined plot for given pair of variable names:

plot_var_pair <- function(
    
    var_pair,
    scores_list,
    n_row,
    n_col,
    sc_genes = NULL,
    simulated_genes = NULL,
    tcga_genes = NULL,
    collate_genes_fun = intersect, # Could also use union
    annotate_coords = c(3/4, 20/21),
    plot_titles = names(scores_list),
    cor_method = 'spearman',
    var_axis_titles = var_pair,
    rel_widths = rep(1, n_col),
    rel_heights = rep(1, n_row)
    
) {
    
    collate_genes_fun <- match.fun(collate_genes_fun)
    
    n <- length(scores_list)
    
    if(
        !is.null(get(paste0(strsplit(var_pair[1], '_')[[1]][1], '_genes'))) &
        !is.null(get(paste0(strsplit(var_pair[2], '_')[[1]][1], '_genes')))
    ) {
        
        plot_data <- sapply(
            names(scores_list),
            function(ct) {
                scores_list[[ct]][
                    gene %in% collate_genes_fun(
                        get(paste0(strsplit(var_pair[1], '_')[[1]][1], '_genes'))[[ct]],
                        get(paste0(strsplit(var_pair[2], '_')[[1]][1], '_genes'))[[ct]]
                    ),
                    ..var_pair
                    ][
                        !is.na(get(var_pair[1])) & !is.na(get(var_pair[2]))
                        ]
            },
            simplify = FALSE,
            USE.NAMES = TRUE
        )
        
    } else {
        
        plot_data <- sapply(
            scores_list,
            function(dt) {
                dt[
                    ,
                    ..var_pair
                    ][
                        !is.na(get(var_pair[1])) & !is.na(get(var_pair[2]))
                        ]
            },
            simplify = FALSE,
            USE.NAMES = TRUE
        )
        
    }
    
    lims_data <- sapply(
        sapply(
            var_pair,
            function(v) {
                sapply(
                    plot_data,
                    function(dt) {
                        c(
                            floor(10*min(dt[[v]]))/10,
                            ceiling(10*max(dt[[v]]))/10
                        )
                    }
                )
            },
            simplify = FALSE,
            USE.NAMES = TRUE
        ),
        function(m) c(min(m[1, ]), max(m[2, ]))
    )
    
    if(n == n_row*n_col) {
        
        plot_grid(
            plotlist = lapply(
                1:n,
                function(i) {
                    
                    ggplot(
                        plot_data[[i]],
                        aes(get(var_pair[1]), get(var_pair[2]))
                    ) +
                        geom_vline(
                            xintercept = 0,
                            linetype = 'dashed',
                            colour = 'grey',
                            size = 0.25
                        ) +
                        geom_hline(
                            yintercept = 0,
                            linetype = 'dashed',
                            colour = 'grey',
                            size = 0.25
                        ) +
                        geom_point() +
                        geom_smooth(method = 'lm', se = FALSE) +
                        theme_test() +
                        lims(x = lims_data[, var_pair[1]], y = lims_data[, var_pair[2]]) +
                        annotate(
                            'text',
                            x = annotate_coords[1]*lims_data[, var_pair[1]][2],
                            y = annotate_coords[2]*lims_data[, var_pair[2]][1],
                            label = latex2exp::TeX(
                                sprintf(
                                    '$\\rho = %.2f$',
                                    plot_data[[i]][
                                        ,
                                        cor(
                                            get(var_pair[1]),
                                            get(var_pair[2]),
                                            method = cor_method
                                        )
                                    ]
                                )
                            )
                        ) +
                        labs(title = plot_titles[i]) +
                        theme(
                            axis.text.x = switch(
                                (i %in% (n - n_col + 1):n) + 1,
                                element_blank(),
                                NULL
                            ),
                            axis.text.y = switch(
                                (i %in% ((0:((n / n_col) - 1))*n_col + 1)) + 1,
                                element_blank(),
                                NULL
                            ),
                            axis.title = element_blank(),
                            plot.margin = unit(
                                c(
                                    5.5,
                                    5.5,
                                    switch(
                                        (i %in% (n - n_col + 1):n) + 1,
                                        5.5,
                                        20
                                    ),
                                    switch(
                                        (i %in% ((0:((n / n_col) - 1))*n_col + 1)) + 1,
                                        5.5,
                                        20
                                    )
                                ),
                                'pt'
                            )
                        )
                    
                }
            ),
            nrow = n_row,
            ncol = n_col,
            # Don't need align = 'h' any more because no plots that aren't on the bottom
            # row will have x axis text, and aligning introduces space where the x axis
            # text would go if it was there...
            rel_widths = rel_widths,
            rel_heights = rel_heights
        ) +
            cowplot::draw_label(
                var_axis_titles[1],
                x = 0.5,
                y = 0,
                vjust = -0.5,
                size = 12
            ) +
            cowplot::draw_label(
                var_axis_titles[2],
                x = 0,
                y = 0.5,
                vjust = 1.3,
                angle = 90,
                size = 12
            )
        
    } else {
        
        plot_grid(
            plot_grid(
                plotlist = lapply(
                    1:(n_col*(n %/% n_col)),
                    function(i) {
                        
                        ggplot(
                            plot_data[[i]],
                            aes(get(var_pair[1]), get(var_pair[2]))
                        ) +
                            geom_vline(
                                xintercept = 0,
                                linetype = 'dashed',
                                colour = 'grey',
                                size = 0.25
                            ) +
                            geom_hline(
                                yintercept = 0,
                                linetype = 'dashed',
                                colour = 'grey',
                                size = 0.25
                            ) +
                            geom_point() +
                            geom_smooth(method = 'lm', se = FALSE) +
                            theme_test() +
                            lims(x = lims_data[, var_pair[1]], y = lims_data[, var_pair[2]]) +
                            annotate(
                                'text',
                                x = annotate_coords[1]*lims_data[, var_pair[1]][2],
                                y = annotate_coords[2]*lims_data[, var_pair[2]][1],
                                label = latex2exp::TeX(
                                    sprintf(
                                        '$\\rho = %.2f$',
                                        plot_data[[i]][
                                            ,
                                            cor(
                                                get(var_pair[1]),
                                                get(var_pair[2]),
                                                method = cor_method
                                            )
                                            ]
                                    )
                                )
                            ) +
                            labs(title = plot_titles[i]) +
                            theme(
                                axis.text.x = switch(
                                    (i %in% (n_col*(n %/% n_col) - (n_col - (n %% n_col)) + 1):(n_col*(n %/% n_col))) + 1,
                                    element_blank(),
                                    NULL
                                ),
                                axis.text.y = switch(
                                    (i %in% ((0:((n %/% n_col) - 1))*n_col + 1)) + 1,
                                    element_blank(),
                                    NULL
                                ),
                                axis.title = element_blank(),
                                plot.margin = unit(
                                    c(
                                        5.5,
                                        5.5,
                                        5.5,
                                        switch(
                                            (i %in% ((0:((n %/% n_col) - 1))*n_col + 1)) + 1,
                                            5.5,
                                            20
                                        )
                                    ),
                                    'pt'
                                )
                            )
                        
                    }
                ),
                nrow = n %/% n_col,
                ncol = n_col,
                align = 'h',
                rel_widths = rel_widths,
                rel_heights = rel_heights[-length(rel_heights)]
            ),
            plot_grid(
                plotlist = lapply(
                    (n - (n %% n_col) + 1):n,
                    function(i) {
                        
                        ggplot(
                            plot_data[[i]],
                            aes(get(var_pair[1]), get(var_pair[2]))
                        ) +
                            geom_vline(
                                xintercept = 0,
                                linetype = 'dashed',
                                colour = 'grey',
                                size = 0.25
                            ) +
                            geom_hline(
                                yintercept = 0,
                                linetype = 'dashed',
                                colour = 'grey',
                                size = 0.25
                            ) +
                            geom_point() +
                            geom_smooth(method = 'lm', se = FALSE) +
                            theme_test() +
                            lims(x = lims_data[, var_pair[1]], y = lims_data[, var_pair[2]]) +
                            annotate(
                                'text',
                                x = annotate_coords[1]*lims_data[, var_pair[1]][2],
                                y = annotate_coords[2]*lims_data[, var_pair[2]][1],
                                label = latex2exp::TeX(
                                    sprintf(
                                        '$\\rho = %.2f$',
                                        plot_data[[i]][
                                            ,
                                            cor(
                                                get(var_pair[1]),
                                                get(var_pair[2]),
                                                method = cor_method
                                            )
                                            ]
                                    )
                                )
                            ) +
                            labs(title = plot_titles[i]) +
                            theme(
                                axis.text.y = switch(
                                    (i == n - (n %% n_col) + 1) + 1,
                                    element_blank(),
                                    NULL
                                ),
                                axis.title = element_blank(),
                                plot.margin = unit(
                                    c(
                                        5.5,
                                        5.5,
                                        20,
                                        switch((i == n - (n %% n_col) + 1) + 1, 5.5, 20)
                                    ),
                                    'pt'
                                )
                            )
                        
                    }
                ),
                nrow = 1,
                ncol = n_col,
                align = 'h',
                rel_widths = rel_widths
            ),
            nrow = 2,
            ncol = 1,
            rel_heights = c(
                sum(rel_heights[-length(rel_heights)]),
                rel_heights[length(rel_heights)]
            )
        ) +
            cowplot::draw_label(
                var_axis_titles[1],
                x = 0.5,
                y = 0,
                vjust = -0.5,
                size = 12
            ) +
            cowplot::draw_label(
                var_axis_titles[2],
                x = 0,
                y = 0.5,
                vjust = 1.3,
                angle = 90,
                size = 12
            )
        
    }
    
}

# Not sure why we need these as well, but hey ho:

sc_genes <- sapply(
    sc_cancer_fibroblast,
    `[[`,
    'genes_filtered',
    simplify = FALSE,
    USE.NAMES = TRUE
)

tcga_genes <- sapply(
    names(sc_to_bulk_names),
    function(ct) {
        unique(
            unlist(
                lapply(deconv_data[sc_to_bulk_names[[ct]]], `[[`, 'genes_filtered')
            )
        )
    },
    simplify = FALSE,
    USE.NAMES = TRUE
)

# Signature concordance (TCGA deconv vs. scRNAseq) for supplementary:

pdf(
    '../data_and_figures/final_figures/S5.pdf',
    width = 6,
    height = 6
)

plot_var_pair(
    c('sc_score', 'tcga_score'),
    scores,
    n_row = 2,
    n_col = 2,
    sc_genes = sc_genes,
    tcga_genes = tcga_genes,
    collate_genes_fun = union,
    plot_titles = sapply(
        sc_cancer_fibroblast_heatmaps_args,
        `[[`,
        'annotations_title'
    ),
    var_axis_titles = c('scRNA-seq score', 'TCGA deconvolution score'),
    rel_widths = c(1.115, 1),
    rel_heights = c(1, 1, 1.075)
)

dev.off()

# Free up some space:

rm(scores_data_tcga)
rm(scores)





# Inter- vs. intra-tumour pEMT heterogeneity (or rare vs. shared EMT):

inter_intra_emt_all_cts <- readRDS('../data_and_figures/inter_intra_emt_all_cts.rds')

# Plots combining cancer types:

emt_scores_all_cts <- rbindlist(
    lapply(
        names(inter_intra_emt_all_cts),
        function(ct) {
            merge(
                cbind(
                    cancer_type = ct,
                    cancer_type_nice = sc_cancer_fibroblast_heatmaps_args[[ct]]$annotations_title,
                    inter_intra_emt_all_cts[[ct]]$data$cell_emt_scores
                ),
                inter_intra_emt_all_cts[[ct]]$data$inter_intra_emt_scores,
                by = 'patient'
            )
        }
    )
)[
    ,
    patient_unique := paste(cancer_type, patient, sep = '.')
]

inter_intra_emt_profiles <- ggplot(
    emt_scores_all_cts,
    aes(
        pos_frac,
        emt_score,
        group = patient,
        # colour = cancer_type_nice
        colour = scale(intra_emt)
    )
) +
    facet_grid(cols = vars(cancer_type_nice)) +
    geom_line() +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0), breaks = c(-5, 0, 5)) +
    scale_colour_gradientn(colours = colorRamps::blue2green(50)[11:45]) +
    theme(
        panel.background = element_rect(fill = NA, colour = 'black'),
        panel.grid.major.y = element_line(colour = 'grey', size = 0.3, linetype = 'dotted'),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.x = element_text(hjust = 0),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_rect(colour = 'black'),
        strip.text = element_text(size = 13),
        # legend.position = 'bottom',
        legend.position = c(1, -0.15),
        legend.direction = 'horizontal',
        legend.justification = 'right',
        legend.key.height = unit(10, 'pt'),
        legend.title = element_text(margin = margin(0, 2.5, 0, 0)),
        plot.margin = unit(c(5.5, 5.5, 40, 5.5), 'pt')
    ) +
    labs(x = 'Cells', y = 'pEMT score', colour = 'Intra-tumour pEMT heterogeneity score\n')

ct_colours = c(
    'Head and Neck' = RColorBrewer::brewer.pal(12, 'Set3')[4],
    'Lung' = RColorBrewer::brewer.pal(12, 'Set3')[1],
    # 'Pancreatic' = RColorBrewer::brewer.pal(12, 'Set3')[10],
    'Pancreatic' = colorRampPalette(
        c(
            RColorBrewer::brewer.pal(12, 'Set3')[10],
            RColorBrewer::brewer.pal(12, 'Set3')[8]
        )
    )(51)[21],
    'Liver' = RColorBrewer::brewer.pal(8, 'Set2')[6]
)

inter_intra_emt_profiles_gtable <- ggplot_gtable(ggplot_build(inter_intra_emt_profiles))

stript <- which(grepl('strip-t', inter_intra_emt_profiles_gtable$layout$name))

for (i in stript) {
    
    inter_intra_emt_profiles_gtable$grobs[[i]]$grobs[[1]]$children[[
        which(grepl('rect', inter_intra_emt_profiles_gtable$grobs[[i]]$grobs[[1]]$childrenOrder))
    ]]$gp$fill <-
        # This sapply() loop makes the ct_colours lighter:
        sapply(ct_colours, function(x) colorRampPalette(c(x, 'white'))(10)[4])[
            inter_intra_emt_profiles_gtable$grobs[[i]]$grobs[[1]]$children[[
                which(grepl('titleGrob', inter_intra_emt_profiles_gtable$grobs[[i]]$grobs[[1]]$childrenOrder))
            ]]$children[[1]]$label
        ]
    
}

inter_intra_emt_jitterplot_data <- melt(
    unique(
        emt_scores_all_cts[
            ,
            .(
                cancer_type_nice = cancer_type_nice,
                `Inter-tumour pEMT heterogeneity` = scale(inter_emt)[, 1],
                `Intra-tumour pEMT heterogeneity` = scale(intra_emt)[, 1]
            )
        ]
    ),
    id.vars = 'cancer_type_nice'
)

# To get a p value for the difference between HNSCC/PDAC and lung/liver:

# inter_intra_emt_jitterplot_data[
#     variable == 'Intra-tumour pEMT heterogeneity',
#     apply(
#         combn(unique(cancer_type_nice), 2),
#         2,
#         function(ct_pair) {
#             list(
#                 cancer_types = ct_pair,
#                 pval = t.test(
#                     .SD[cancer_type_nice == ct_pair[1], value],
#                     .SD[cancer_type_nice == ct_pair[2], value]
#                 )$p.value
#             )
#         }
#     )
# ]

# inter_intra_emt_boxplots <- ggplot(
#     inter_intra_emt_boxplot_data,
#     aes(x = cancer_type_nice, y = value)
# ) +
#     # geom_boxplot(varwidth = TRUE, outlier.shape = NA) +
#     # stat_boxplot(geom = 'errorbar', width = 0.2) +
#     geom_boxplot(outlier.shape = NA, width = 0.5) +
#     geom_jitter(width = 0.2) +
#     facet_grid(cols = vars(variable)) +
#     theme_bw()

inter_intra_emt_jitterplot <- ggplot(
    inter_intra_emt_jitterplot_data[
        ,
        .(
            y0 = mean(value) - sd(value),#/sqrt(length(value)),
            y1 = mean(value),
            y2 = mean(value) + sd(value)#/sqrt(length(value))
        ),
        by = .(cancer_type_nice, variable)
    ],#[variable == 'intra_emt_scaled'],
    aes(x = cancer_type_nice)
) +
    # geom_boxplot(
    #     aes(ymin = y0, lower = y1, middle = y1, upper = y1, ymax = y2),
    #     stat = 'identity'
    # ) +
    geom_segment(
        aes(x = cancer_type_nice, xend = cancer_type_nice, y = y0, yend = y2)
    ) +
    geom_segment(
        aes(
            x = as.numeric(as.factor(cancer_type_nice)) - 0.05,
            xend = as.numeric(as.factor(cancer_type_nice)) + 0.05,
            y = y0,
            yend = y0
        )
    ) +
    geom_segment(
        aes(
            x = as.numeric(as.factor(cancer_type_nice)) - 0.05,
            xend = as.numeric(as.factor(cancer_type_nice)) + 0.05,
            y = y2,
            yend = y2
        )
    ) +
    geom_point(aes(x = cancer_type_nice, y = y1), shape = 10, size = 3) +
    geom_jitter(
        data = inter_intra_emt_jitterplot_data,
        aes(x = cancer_type_nice, y = value, colour = cancer_type_nice),
        width = 0.2,
        size = 2
    ) +
    scale_colour_manual(values = ct_colours) +
    facet_grid(cols = vars(variable)) +
    theme_bw() +
    theme(
        axis.text.x = element_text(size = 11, angle = 50, hjust = 1),
        strip.text = element_text(size = 11),
        strip.background = element_rect(fill = NA, colour = NA),
        legend.position = 'none'
    ) +
    labs(x = NULL, y = 'Score')

# grid::grid.draw(inter_intra_emt_profiles_gtable)

# inter_intra_emt_scatterplot <- ggplot(
#     unique(
#         emt_scores_all_cts[
#             ,
#             .(
#                 cancer_type_nice = cancer_type_nice,
#                 inter_emt_scaled = scale(inter_emt)[, 1],
#                 intra_emt_scaled = scale(intra_emt)[, 1]
#             )
#         ]
#     ),
#     aes(inter_emt_scaled, intra_emt_scaled, colour = cancer_type_nice)
# ) +
#     geom_hline(yintercept = 0, size = 0.3, colour = 'grey', linetype = 'dotted') +
#     geom_vline(xintercept = 0, size = 0.3, colour = 'grey', linetype = 'dotted') +
#     geom_point(size = 2) +
#     # geom_point(
#     #     data = unique(
#     #         emt_scores_all_cts[
#     #             ,
#     #             .(
#     #                 cancer_type_nice = cancer_type_nice,
#     #                 inter_emt_scaled = scale(inter_emt)[, 1],
#     #                 intra_emt_scaled = scale(intra_emt)[, 1]
#     #             )
#     #         ][
#     #             ,
#     #             .(
#     #                 inter_emt_scaled = mean(inter_emt_scaled),
#     #                 intra_emt_scaled = mean(intra_emt_scaled)
#     #             ),
#     #             by = cancer_type_nice
#     #         ]
#     #     ),
#     #     size = 3,
#     #     shape = 3,
#     #     stroke = 2
#     # ) +
#     scale_colour_manual(values = ct_colours) +
#     theme_test() +
#     theme(
#         legend.title = element_text(size = 13, hjust = 1),
#         legend.text = element_text(size = 12),
#         legend.justification = c(0.85, 0.3),
#         axis.title.y = element_text(vjust = 3)
#     ) +
#     guides(colour = guide_legend(label.position = 'left')) +
#     labs(
#         x = 'Inter-tumour pEMT heterogeneity score', # These are medians, but then scaled...
#         y = 'Intra-tumour pEMT heterogeneity score',
#         colour = 'Cancer type'
#     )

# To make a "scheme" figure, showing what the inter- and intra-tumour heterogeneity
# scores mean in terms of the profiles:

# In the following, the colour '#5B8BAC' is half way between skyblue3 and skyblue4, and
# can be obtained by the command colorRampPalette(c('skyblue3', 'skyblue4'))(3)[2].

inter_intra_emt_scheme <- ggplot(
    data.table(
        y = c(
            seq(qnorm(0.01, mean = 0, sd = 1), qnorm(0.99, mean = 0, sd = 1), length.out = 500),
            seq(qnorm(0.01, mean = -1.5, sd = 0.6), qnorm(0.99, mean = -1.5, sd = 0.6), length.out = 500)
        ),
        grp = c(rep(2, 500), rep(1, 500))
    )[
        ,
        x := c(
            pnorm(y[1:500], mean = 0, sd = 1),
            pnorm(y[501:1000], mean = -1.5, sd = 0.6)
        )
    ],
    aes(x = x, y = y)
) +
    geom_line(aes(group = as.character(grp), colour = as.character(grp))) +
    ggrepel::geom_text_repel(
        aes(x, y, label = l, colour = as.character(grp)),
        data = data.table(x = 0.5, y = qnorm(0.5, mean = -1.5, sd = 0.6), l = 'median'),
        nudge_x = 0.15,
        nudge_y = -0.4,
        segment.colour = 'darkgrey',
        size = 4,
        # colour = 'sienna3'
        colour = colorRamps::blue2green(50)[11]
    ) +
    ggrepel::geom_text_repel(
        aes(x, y, label = l),
        data = data.table(x = 0.5, y = qnorm(0.5), l = 'median'),
        nudge_x = 0.05,
        nudge_y = 0.7,
        segment.colour = 'darkgrey',
        size = 4,
        # colour = '#5B8BAC'
        colour = colorRamps::blue2green(50)[36]
    ) +
    ggrepel::geom_text_repel(
        aes(x, y, label = l),
        data = data.table(x = 0.95, y = qnorm(0.95), l = '95th percentile'),
        nudge_x = -0.2,
        nudge_y = 0.4,
        segment.colour = 'darkgrey',
        size = 4,
        # colour = '#5B8BAC'
        colour = colorRamps::blue2green(50)[36]
    ) +
    # scale_colour_manual(values = c('1' = 'sienna3', '2' = '#5B8BAC')) +
    scale_colour_manual(
        values = c(
            '1' = colorRamps::blue2green(50)[11],
            '2' = colorRamps::blue2green(50)[36]
        )
    ) +
    scale_x_continuous(
        expand = c(0.01, 0)
    ) +
    scale_y_continuous(
        limits = c(qnorm(0.01, mean = -1.5, sd = 0.6) - 0.2, qnorm(0.99, mean = 0, sd = 1) + 0.2),
        expand = c(0, 0)
    ) +
    geom_segment( # Arrow between profiles to indicate inter-tumour heterogeneity
        aes(
            x = 0.5,
            y = qnorm(0.5, mean = -1.5, sd = 0.6) + 0.075,
            xend = 0.5,
            yend = qnorm(0.5, mean = 0, sd = 1) - 0.075
        ),
        arrow = arrow(ends = 'both', length = unit(5, 'pt'))
    ) +
    geom_segment( # Lower horizontal dashed line
        aes(
            x = 0.4025,
            y = qnorm(0.5, mean = 0, sd = 1),
            xend = 0.4975,
            yend = qnorm(0.5, mean = 0, sd = 1)
        ),
        linetype = 'dashed',
        colour = 'grey',
        size = 0.25
    ) +
    geom_segment( # Upper horizontal dashed line
        aes(
            x = 0.4025,
            y = qnorm(0.95, mean = 0, sd = 1),
            xend = 0.9475,
            yend = qnorm(0.95, mean = 0, sd = 1)
        ),
        linetype = 'dashed',
        colour = 'grey',
        size = 0.25
    ) +
    geom_segment( # Arrow between horizontal dashed lines to indicate intra-tumour heterogeneity
        aes(
            x = 0.4,
            y = qnorm(0.5, mean = 0, sd = 1) + 0.05,
            xend = 0.4,
            yend = qnorm(0.95, mean = 0, sd = 1) - 0.05
        ),
        arrow = arrow(ends = 'both', length = unit(5, 'pt'))
    ) +
    geom_point(
        aes(x, y, colour = as.character(grp)),
        data = data.table(
            x = c(0.5, 0.5, 0.95),
            y = c(qnorm(0.5, mean = -1.5, sd = 0.6), qnorm(0.5), qnorm(0.95)),
            grp = c(1, 2, 2)
        )
    ) +
    annotate(
        geom = 'text',
        x = 0.15,
        y = qnorm(0.04, mean = -1.5, sd = 0.6),
        label = 'tumour 1',
        size = 4,
        colour = colorRamps::blue2green(50)[11]
    ) +
    annotate(
        geom = 'text',
        x = 0.19,
        y = qnorm(0.08, mean = 0, sd = 1),
        label = 'tumour 2',
        size = 4,
        colour = colorRamps::blue2green(50)[36]
    ) +
    annotate(
        geom = 'text',
        x = 0.38,
        y = qnorm(0.5, mean = 0, sd = 1) + (qnorm(0.95, mean = 0, sd = 1) - qnorm(0.5, mean = 0, sd = 1))/2,
        label = 'Intra-tumour\nheterogeneity',
        hjust = 1,
        size = 4
    ) +
    annotate(
        geom = 'text',
        x = 0.52,
        y = qnorm(0.5, mean = -1.5, sd = 0.6) + 1,
        label = 'Inter-tumour\nheterogeneity',
        hjust = 0,
        size = 4
    ) +
    theme_test() +
    theme(
        # axis.title.x = element_text(vjust = 6),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = 'none'
    ) +
    labs(x = 'Cells ranked by pEMT score', y = 'pEMT score')

# Align plots so we can combine them into one figure:

# aligned_plots_inter_intra_1 <- align_plots(
#     inter_intra_emt_scheme,
#     inter_intra_emt_scatterplot,
#     align = 'h'
# )
# 
# aligned_plots_inter_intra_2 <- align_plots(
#     inter_intra_emt_scheme,
#     inter_intra_emt_profiles,
#     align = 'v',
#     axis = 'l'
# )
# 
# aligned_plots_inter_intra_1[[1]]$widths[1] <- aligned_plots_inter_intra_2[[1]]$widths[1]

# I think it is possible to manually change the facet label colours to match those used in
# the scatterplot.  I think you would have to make the facet plot into a grob, then delve
# into the grob elements.  See this stackoverflow post:

# https://stackoverflow.com/questions/22457981/ggplot2-facet-grid-strip-text-x-different-colours-based-on-factor

# I don't think it's worth the effort, though.

# Inter- vs. intra-tumour pEMT heterogeneity:

pdf('../data_and_figures/final_figures/5.pdf', width = 11, height = 7.5)

plot_grid(
    blank_plot(),
    plot_grid(
        plot_grid(
            inter_intra_emt_scheme + theme(plot.margin = unit(c(5.5, 40, 5.5, 20), 'pt')),
            # get_legend(inter_intra_emt_scatterplot),
            nrow = 2,
            ncol = 1,
            rel_heights = c(2, 0.7)
        ),
        plot_grid(
            # inter_intra_emt_profiles,
            inter_intra_emt_profiles_gtable,
            # inter_intra_emt_scatterplot + theme(legend.position = 'none'),
            inter_intra_emt_jitterplot + theme(plot.margin = unit(c(20, 5.5, 5.5, 5.5), 'pt')),
            nrow = 2,
            ncol = 1,
            axis = 'l',
            align = 'v',
            rel_heights = c(5, 6)
        ),
        nrow = 1,
        ncol = 2,
        rel_widths = c(1.9, 3)
    ),
    nrow = 2,
    ncol = 1,
    rel_heights = c(1, 20)
) +
    draw_label('A', x = 0, y = 0.975, hjust = -0.1, vjust = 1.1, size = 20, fontface = 2) +
    draw_label('B', x = 0.39, y = 0.975, hjust = -0.1, vjust = 1.1, size = 20, fontface = 2) +
    draw_label('C', x = 0.39, y = 0.51, hjust = -0.1, vjust = 1.1, size = 20, fontface = 2)

# plot_grid(
#     blank_plot(),
#     plot_grid(
#         aligned_plots_inter_intra_1[[1]],
#         blank_plot(),
#         aligned_plots_inter_intra_1[[2]],
#         nrow = 1,
#         ncol = 3,
#         rel_widths = c(2.3, 0.7, 4)
#     ),
#     blank_plot(),
#     aligned_plots_inter_intra_2[[2]],
#     nrow = 4,
#     ncol = 1,
#     rel_heights = c(0.1, 1.05, 0.05, 1)
# ) +
#     draw_label('A', x = 0, y = 0.975, hjust = -0.1, vjust = 1.1, size = 16, fontface = 2) +
#     draw_label('B', x = 0.425, y = 0.975, hjust = -0.1, vjust = 1.1, size = 16, fontface = 2) +
#     draw_label('C', x = 0, y = 0.475, hjust = -0.1, vjust = 1.1, size = 16, fontface = 2)

dev.off()

# Free up some space:

rm(inter_intra_emt_all_cts)
rm(emt_scores_all_cts)
rm(inter_intra_emt_profiles)
rm(inter_intra_emt_profiles_gtable)
rm(inter_intra_emt_jitterplot_data)
rm(inter_intra_emt_jitterplot)
# rm(inter_intra_emt_scatterplot)
rm(inter_intra_emt_scheme)





# All deconv figures, for supplementary:

plots_rel_heights <- c(
    title = 4,
    purity_bar = 1,
    ccle_bar = 1,
    extra_bar = 1,
    heatmap = 15
)

pdf(
    '../data_and_figures/final_figures/S3.pdf',
    width = 16,
    height = 13.5
)

plot_grid(
    blank_plot(),
    plot_grid(
        big_deconv_plot(
            c('BRCA', 'LUSC', 'HNSC', 'LUAD'),
            nice_names_for_figure,
            plots_rel_heights = plots_rel_heights,
            title_font_size = 16
        ),
        plot_grid(
            big_deconv_plot(
                'COAD',
                nice_names_for_figure,
                plots_rel_heights = plots_rel_heights,
                title_font_size = 16
            ),
            blank_plot(),
            big_deconv_plot(
                'READ',
                nice_names_for_figure,
                plots_rel_heights = plots_rel_heights,
                title_font_size = 16
            ),
            blank_plot(),
            nrow = 1,
            ncol = 4,
            rel_widths = c(1.075, 0.625, 1.075, 1.225)
        ),
        plot_grid(
            big_deconv_plot(
                'LIHC',
                nice_names_for_figure,
                plots_rel_heights = plots_rel_heights,
                title_font_size = 16
            ),
            blank_plot(),
            big_deconv_plot(
                'PAAD',
                nice_names_for_figure,
                plots_rel_heights = plots_rel_heights,
                title_font_size = 16
            ),
            blank_plot(),
            nrow = 1,
            ncol = 4,
            rel_widths = c(1.075, 0.625, 1.075, 1.225)
        ),
        nrow = 3,
        ncol = 1,
        rel_heights = c(4, 1, 1)
    ),
    blank_plot(),
    plot_grid(
        plot_grid(
            big_deconv_plot(
                c('OV', 'SKCM', 'STAD'),
                nice_names_for_figure,
                plots_rel_heights = plots_rel_heights[c('title', 'purity_bar', 'ccle_bar', 'heatmap')],
                title_font_size = 16
            ),
            blank_plot(),
            plot_grid(
                big_deconv_plot(
                    'CESC',
                    nice_names_for_figure,
                    plots_rel_heights = plots_rel_heights[c('title', 'purity_bar', 'heatmap')],
                    title_font_size = 16
                ),
                blank_plot(),
                nrow = 2,
                ncol = 1,
                rel_heights = c(1, 2.15)
            ),
            nrow = 1,
            ncol = 3,
            rel_widths = c(3, 0.2, 1.05)
        ),
        plot_grid(
            big_deconv_plot(
                c('BLCA', 'ESCA', 'UCEC'),
                nice_names_for_figure,
                plots_rel_heights = plots_rel_heights[c('title', 'purity_bar', 'ccle_bar', 'heatmap')],
                title_font_size = 16
            ),
            blank_plot(),
            plot_grid(
                plotlist = list(
                    blank_plot(),
                    get_legend(
                        deconv_plots$LIHC$plots$purity_bar + labs(
                            fill = 'Correlation with purity'
                        ) + theme(
                            legend.direction = 'horizontal',
                            legend.justification = 'left'
                        ) + guides(fill = guide_colourbar(title.position = 'right'))
                    ),
                    get_legend(
                        deconv_plots$LIHC$plots$ccle_bar + labs(
                            fill = 'Tumours vs. cell lines'
                        ) + theme(
                            legend.direction = 'horizontal',
                            legend.justification = 'left'
                        ) + guides(fill = guide_colourbar(title.position = 'right'))
                    ),
                    get_legend(
                        deconv_plots$LIHC$plots$extra_bar + labs(
                            fill = 'scRNA-seq: CAF vs. cancer'
                        ) + theme(
                            legend.direction = 'horizontal',
                            legend.justification = 'left'
                        ) + guides(fill = guide_colourbar(title.position = 'right'))
                    ),
                    get_legend(
                        deconv_plots$LIHC$plots$heatmap + labs(
                            fill = 'Correlation'
                        ) + theme(
                            legend.justification = 'left'
                        ) + guides(fill = guide_colourbar(title.position = 'right'))
                    )
                ),
                nrow = 5,
                ncol = 1,
                rel_heights = c(6, 1, 1, 1, 3)
            ),
            nrow = 1,
            ncol = 3,
            rel_widths = c(2, 0.2, 2)
        ),
        blank_plot(),
        nrow = 3,
        ncol = 1,
        rel_heights = c(1, 1, 0.1)
    ),
    nrow = 1,
    ncol = 4,
    rel_widths = c(0.2, 3.95, 0.2, 4.25)
)

dev.off()





# Plots for main deconv figure:

# First remake the plots for these cancer types, so we can change the axis labels:

ct_subset <- c(
    'BRCA - Luminal A',
    'HNSC - Malignant-Basal',
    'PAAD',
    'BLCA - Basal-Squamous',
    'ESCA - Adenocarcinoma' ,
    'UCEC'
)

heatmap_annotations_subset <- lapply(
    list(
        `BRCA - Luminal A` = c('CALU', 'COL1A1', 'FAP', 'ITGAV', 'LAMC2', 'NOTCH2', 'PCOLCE', 'PVR', 'THY1'),
        `HNSC - Malignant-Basal` = c('ACTA2', 'CD44', 'COL6A1', 'ITGB1', 'LAMC2', 'POSTN', 'SERPINE1', 'TAGLN', 'TGFBI'),
        `PAAD` = c('AREG', 'CD44', 'COL1A2', 'COL8A1', 'ITGB5', 'LAMC2', 'MMP2', 'PDPN', 'TGFBI'),
        `BLCA - Basal-Squamous` = c('CD44', 'TGFBI', 'LAMC2', 'ITGA2', 'PVR', 'TAGLN', 'COL1A2', 'PCOLCE', 'COL1A1'),
        `ESCA - Adenocarcinoma` = c('AREG', 'COL1A1', 'COL1A2', 'CXCL8', 'FAP', 'ITGA2', 'LAMC2', 'MMP3', 'PVR'),
        `UCEC` = c('ACTA2', 'CALU', 'COL3A1', 'DCN', 'FAP', 'FBN2', 'ITGB1', 'NOTCH2', 'TNC')
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
                        RColorBrewer::brewer.pal(11, "PuOr")
                    )(50)
                ),
                purity_legend_title = 'Correlation\nwith purity',
                purity_legend_direction = 'horizontal',
                purity_axis_title = NULL,
                ccle_colours = rev(
                    colorRampPalette(
                        RColorBrewer::brewer.pal(11, 'PiYG')[c(1:3, 6, 9:11)]
                    )(50)
                ),
                ccle_legend_title = 'Tumours vs.\ncell lines',
                ccle_legend_direction = 'horizontal',
                ccle_axis_title = NULL,
                extra_colours = colorRampPalette(
                    c('turquoise4', 'turquoise', 'azure', 'gold2', 'gold4')
                )(50),
                extra_legend_title = 'scRNA-seq:\nCAF vs. cancer',
                extra_legend_direction = 'horizontal',
                extra_axis_title = NULL,
                bar_legend_justification = 'left'
                # bar_legend_width = NULL
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
        shape = as.character(
            plyr::mapvalues(
                n_annot,
                c(1, 2, 3),
                c('No', 'No', 'Yes'),
                warn_missing = FALSE
            )
        )
    )
) +
    geom_point(size = 2, alpha = 0.75, colour = 'lightblue4') +
    scale_y_continuous(breaks = c(0.2, 0.3, 0.4)) +
    ggrepel::geom_text_repel(
        aes(label = cancer_type),
        data = within_between_clust_corr[cancer_type %in% c('BLCA - Basal-Squamous', 'ESCA - Squamous')],
        nudge_x = 0.04, nudge_y = 0.01, direction = 'y', size = 3.5
    ) +
    ggrepel::geom_text_repel(
        aes(label = cancer_type),
        data = within_between_clust_corr[cancer_type %in% c('HNSC - Malignant-Basal', 'HNSC - Classical', 'UCEC')],
        nudge_x = 0.02, nudge_y = 0.005, direction = 'x', size = 3.5
    ) +
    ggrepel::geom_text_repel(
        aes(label = cancer_type),
        data = within_between_clust_corr[cancer_type == 'COAD'],
        nudge_x = 0.02, nudge_y = 0.02, size = 3.5
    ) +
    ggrepel::geom_text_repel(
        aes(label = cancer_type),
        data = within_between_clust_corr[cancer_type == 'LIHC'],
        nudge_x = -0.02, nudge_y = 0.01, size = 3.5
    ) +
    ggrepel::geom_text_repel(
        aes(label = cancer_type),
        data = within_between_clust_corr[cancer_type == 'STAD - MSI'],
        nudge_x = -0.04, nudge_y = -0.02, size = 3.5
    ) +
    ggrepel::geom_text_repel(
        aes(label = cancer_type),
        data = within_between_clust_corr[cancer_type == 'BRCA - Luminal A'],
        nudge_x = -0.05, nudge_y = 0.01, size = 3.5
    ) +
    ggrepel::geom_text_repel(
        aes(label = cancer_type),
        data = within_between_clust_corr[cancer_type == 'ESCA - Adenocarcinoma'],
        nudge_x = -0.06, nudge_y = -0.01, size = 3.5
    ) +
    ggrepel::geom_text_repel(
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

# Table providing a dictionary between TCGA disease codes and common cancer type names:

tcga_codes_table <- gridExtra::tableGrob(
    data.table(
        `TCGA code` = c(
            'BLCA', 'BRCA', 'CESC', 'COAD', 'ESCA', 'HNSC', 'LIHC', 'LUAD', 'LUSC',
            'OV', 'PAAD', 'READ', 'SKCM', 'STAD', 'UCEC'
        ),
        `Cancer type` = c(
            'Bladder', 'Breast', 'Cervical', 'Colon', 'Oesophageal', 'Head and Neck',
            'Liver', 'Lung (Adeno.)', 'Lung (Squamous)', 'Ovarian', 'Pancreatic',
            'Rectum', 'Melanoma', 'Stomach', 'Endometrial'
        )
    ),
    theme = gridExtra::ttheme_minimal(
        padding = unit(c(10, 7.5), 'pt'),
        core = list(fg_params = list(x = 0.05, hjust = 0, fontsize = 11), padding = unit(c(15, 7.5), 'pt')),
        colhead = list(padding = unit(c(15, 15), 'pt'))
    ),
    rows = NULL
)

tcga_codes_table <- gtable::gtable_add_grob(
    tcga_codes_table,
    rectGrob(gp = gpar(fill = NA, lwd = 2)),
    t = 1, b = nrow(tcga_codes_table), l = 1, r = 2
)

tcga_codes_table <- gtable::gtable_add_grob(
    tcga_codes_table,
    segmentsGrob(
        x0 = unit(0, 'npc'), x1 = unit(0, 'npc'), y0 = unit(0, 'npc'), y1 = unit(1, 'npc'),
        gp = gpar(lwd = 2)
    ),
    t = 1, b = nrow(tcga_codes_table), l = 2, r = 2
)

tcga_codes_table <- gtable::gtable_add_grob(
    tcga_codes_table,
    segmentsGrob(
        x0 = unit(0, 'npc'), x1 = unit(1, 'npc'), y0 = unit(0, 'npc'), y1 = unit(0, 'npc'),
        gp = gpar(lwd = 2)
    ),
    t = 1, b = 1, l = 1, r = ncol(tcga_codes_table)
)

tcga_codes_table_width = sum(tcga_codes_table$widths)
tcga_codes_table_height = sum(tcga_codes_table$heights)

# To align to top left, with left side aligned with scatterplot's left axis:

tcga_codes_table$vp <- viewport(
    x = unit(0, 'npc') +
        sum(ggplotGrob(deconv_summary_scatterplot)$widths[c(3, 4)]) +
        # Could include index 1 (left margin, 5.5pt) but I change this later (to 40pt) anyway
        unit(40, 'pt') +
        0.5*tcga_codes_table_width,
    y = unit(1, 'npc') - 0.5*tcga_codes_table_height
)

# To align to bottom right:

# tcga_codes_table$vp <- viewport(
#     x = unit(1, 'npc') - 0.5*tcga_codes_table_width - unit(5.5, 'pt'),
#     y = unit(0, 'npc') + 0.5*tcga_codes_table_height
# )

# Parameters for 'pEMT' and 'CAF' brackets under gene annotations:

pemt_caf_bracket_params <- list(
    `BRCA - Luminal A` = list(
        left_bracket_xmax = 32,
        right_bracket_xmin = 69
    ),
    `HNSC - Malignant-Basal` = list(
        left_bracket_xmax = 39,
        right_bracket_xmin = 56
    ),
    `PAAD` = list(
        left_bracket_xmax = 32,
        right_bracket_xmin = 56
    ),
    `BLCA - Basal-Squamous` = list(
        left_bracket_xmax = 39,
        right_bracket_xmin = 63
    ),
    `ESCA - Adenocarcinoma` = list(
        left_bracket_xmax = 39,
        right_bracket_xmin = 56
    ),
    `UCEC` = list(
        left_bracket_xmax = 32,
        right_bracket_xmin = 62
    )
)

# Main deconv plot:

pdf(
    '../data_and_figures/final_figures/3.pdf',
    width = 14.5,
    height = 11
)

plot_grid(
    blank_plot(),
    plot_grid(
        blank_plot(),
        plot_grid(
            plot_grid(
                deconv_plot(
                    deconv_plots_subset[c('BRCA - Luminal A', 'HNSC - Malignant-Basal', 'PAAD')],
                    n_row = 1,
                    n_col = 3,
                    plots_rel_heights = c(
                        title = 2,
                        purity_bar = 1,
                        ccle_bar = 1,
                        extra_bar = 1,
                        heatmap = 15,
                        axis_labels = 5
                    ),
                    left_plot_width = 1.06,
                    legends = FALSE
                ),
                plot_grid(
                    plotlist = c(
                        list(blank_plot()),
                        lapply(
                            c('BRCA - Luminal A', 'HNSC - Malignant-Basal', 'PAAD'),
                            function(ct) {
                                do.call(
                                    pemt_caf_brackets,
                                    args = c(
                                        list(brackets_y = 0.7),
                                        pemt_caf_bracket_params[[ct]]
                                    )
                                ) + theme(plot.margin = unit(c(0, 5.5, 0, 5.5), 'pt'))
                            }
                        )
                    ),
                    nrow = 1,
                    ncol = 4,
                    rel_widths = c(0.06, 1, 1, 1)
                ),
                blank_plot(),
                deconv_plot(
                    deconv_plots_subset[c('BLCA - Basal-Squamous', 'ESCA - Adenocarcinoma' , 'UCEC')],
                    n_row = 1,
                    n_col = 3,
                    plots_rel_heights = c(
                        title = 2,
                        purity_bar = 1,
                        ccle_bar = 1,
                        heatmap = 15,
                        axis_labels = 5
                    ),
                    left_plot_width = 1.06,
                    legends = FALSE
                ),
                plot_grid(
                    plotlist = c(
                        list(blank_plot()),
                        lapply(
                            c('BLCA - Basal-Squamous', 'ESCA - Adenocarcinoma' , 'UCEC'),
                            function(ct) {
                                do.call(
                                    pemt_caf_brackets,
                                    args = c(
                                        list(brackets_y = 0.7),
                                        pemt_caf_bracket_params[[ct]]
                                    )
                                ) + theme(plot.margin = unit(c(0, 5.5, 0, 5.5), 'pt'))
                            }
                        )
                    ),
                    nrow = 1,
                    ncol = 4,
                    rel_widths = c(0.06, 1, 1, 1)
                ),
                nrow = 5,
                ncol = 1,
                rel_heights = c(25, 2, 3, 24, 2)
            ),
            # deconv_plot(
            #     deconv_plots_subset,
            #     n_row = 2,
            #     n_col = 3,
            #     plots_rel_heights = c(
            #         title = 2,
            #         purity_bar = 1,
            #         ccle_bar = 1,
            #         extra_bar = 1,
            #         heatmap = 15,
            #         axis_labels = 6
            #     ),
            #     rows_rel_heights = c(1.04, 1),
            #     left_plot_width = 1.06,
            #     legends = FALSE
            # ),
            plot_grid(
                blank_plot(),
                tcga_codes_table,
                blank_plot(),
                plot_grid(
                    blank_plot(),
                    plot_grid(
                        plotlist = lapply(
                            c('purity_bar', 'ccle_bar', 'extra_bar'),
                            function(b) {
                                get_legend(
                                    deconv_plots_subset$`BRCA - Luminal A`$plots[[b]] +
                                        theme(
                                            legend.justification = 'left',
                                            legend.title = element_text(margin = margin(0, 0, 0, 4))
                                        ) +
                                        guides(fill = guide_colourbar(title.position = 'right'))
                                )
                            }
                        ),
                        nrow = 3,
                        ncol = 1
                    ),
                    get_legend(
                        deconv_plots_subset$`BRCA - Luminal A`$plots$heatmap + guides(
                            fill = guide_colourbar(
                                title.position = 'right',
                                title.hjust = -4
                            )
                        )
                    ),
                    nrow = 1,
                    ncol = 3,
                    rel_widths = c(1.3, 3, 2)
                ),
                # plot_grid(
                #     get_legend(
                #         deconv_plots_subset$`BRCA - Luminal A`$plots$heatmap + theme(
                #             legend.justification = c(0, 0.8)
                #         )
                #     ),
                #     tcga_codes_table,
                #     nrow = 1,
                #     ncol = 2,
                #     rel_widths = c(1, 2)
                # ),
                # plotlist = c(
                #     lapply(
                #         c('purity_bar', 'ccle_bar', 'extra_bar', 'heatmap'),
                #         function(b) {
                #             get_legend(
                #                 deconv_plots_subset$`BRCA - Luminal A`$plots[[b]] + theme(
                #                     legend.justification = 'left'
                #                 ) + guides(fill = guide_colourbar(title.position = 'right'))
                #             )
                #         }
                #     ),
                blank_plot(),
                deconv_summary_scatterplot +
                    theme(
                        legend.position = 'bottom',
                        legend.text = element_text(margin = margin(r = 15, unit = 'pt')),
                        legend.title = element_text(margin = margin(r = 10, unit = 'pt')),
                        legend.spacing.x = unit(1, 'pt'),
                        plot.margin = unit(c(0, 5.5, 5.5, 40), 'pt')
                    ) +
                    labs(shape = 'scRNA-seq data available:'),
                blank_plot(),
                nrow = 7,
                ncol = 1,
                rel_heights = c(2, 18.5, 1.5, 7.5, 2.5, 23, 1)
                # rel_heights = c(26, 2, 3, 3, 3, 8, 6)
            ),
            nrow = 1,
            ncol = 2,
            rel_widths = c(2, 1)
        ),
        nrow = 2,
        ncol = 1,
        rel_heights = c(1, 40)
    ),
    nrow = 1,
    ncol = 2,
    rel_widths = c(1, 50)
) +
    draw_label('A', x = 0, y = 0.99, hjust = -0.1, vjust = 1.1, size = 20, fontface = 2) +
    draw_label('B', x = 0.71, y = 0.971, hjust = -0.1, vjust = 1.1, size = 20, fontface = 2) +
    # draw_label('B', x = 0.78, y = 0.795, hjust = -0.1, vjust = 1.1, size = 20, fontface = 2) +
    draw_label('C', x = 0, y = 0.473, hjust = -0.1, vjust = 1.1, size = 20, fontface = 2) +
    draw_label('D', x = 0.71, y = 0.452, hjust = -0.1, vjust = 1.1, size = 20, fontface = 2)

dev.off()





# Data for deconv summary barplot (for supplementary):

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
        
        names(lm_coeffs) <- stringr::str_split_fixed(names(lm_coeffs), '\\.', 2)[, 1]
        
        annot_fun <- function(x) {caTools::runmean(x, 30)/max(abs(caTools::runmean(x, 30)))}
        
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

# Barplots of agreement with annotations, for supplement:

pdf(
    '../data_and_figures/final_figures/S4.pdf',
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
            'pur' = RColorBrewer::brewer.pal(11, "PuOr")[2],
            'ccle' = RColorBrewer::brewer.pal(11, "PiYG")[3],
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
            'pur' = RColorBrewer::brewer.pal(11, "PuOr")[2],
            'ccle' = RColorBrewer::brewer.pal(11, "PiYG")[3],
            'extra' = 'gold2'
        )
    ) +
    labs( # Remember the axes are flipped!
        x = NULL,
        y = 'Between-cluster difference',
        fill = 'Annotation type'
    )

dev.off()





# Commonality heatmap (including melanoma Keratin and Immune):

rank_mat <- deconv_rank(deconv_data)

scores_data_transformed <- deconv_scores(
    expression_data,
    deconv_data,
    scale_fun = function(x) x/(3*sd(x)),
    transform_data = TRUE
)

set.seed(44398)

htmp_emt_caf <- deconv_heatmap(
    scores_data_transformed[
        unique(
            c(
                names(head(sort(apply(rank_mat, 1, quantile, 0.25)), 50)),
                names(tail(sort(apply(rank_mat, 1, quantile, 0.75)), 50))
            )
        ),
        c('gene', ..deconv_names)
    ],
    order_genes_fun = 'seriate',
    order_genes_method = 'SPIN_NH',
    order_analyses_fun = function(x) {
        rev(
            seriation::get_order(
                seriation::seriate(
                    dist(t(x)),
                    method = 'SPIN_NH'
                )
            )
        )
    }
    # order_analyses_fun = 'seriate',
    # order_analyses_method = 'SPIN_NH',
    # plot_title = 'Common EMT and stroma genes'
)

# PCA to describe the variation in EMT and CAF genes between cancer types:

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

# kclust_map <- plyr::mapvalues(
#     c(1, 2, 3),
#     c(1, 2, 3),
#     pca_data[
#         c('BRCA - Luminal A', 'HNSC - Malignant-Basal', 'ESCA - Adenocarcinoma'),
#         kmeans_cluster
#     ]
# )
# 
# pca_data[
#     ,
#     kmeans_cluster_manual := plyr::mapvalues(
#         kmeans_cluster,
#         c(1, 2, 3),
#         kclust_map
#     )
# ]

pca_data[
    ,
    kmeans_cluster_manual := plyr::mapvalues(
        kmeans_cluster,
        .SD[
            c('BRCA - Luminal A', 'HNSC - Malignant-Basal', 'ESCA - Adenocarcinoma'),
            kmeans_cluster
        ],
        c(1, 2, 3)
    )
]

kclust_map <- pca_data[
    order(kmeans_cluster_manual),
    .(kclust_map = unique(kmeans_cluster)),
    by = kmeans_cluster_manual
]$kclust_map

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

cutoffs <- c(0.4, 0.9, 0.9)

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
            1:20,
            gene
        ]
        
    }
)

# Final heatmap with hierarchical clustering:

htmp_all_kmeans_distinct <- deconv_heatmap(
    scores_data_transformed[
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
                    c(
                        cancer_type,
                        rep(
                            '',
                            max(
                                table(
                                    pca_data[
                                        thresh_pass == TRUE,
                                        kmeans_cluster_manual
                                    ]
                                )
                            ) - .N
                        )
                    )
                ]
            }
        )
    ),
    c(
        'Cluster 1\nGynaecological',
        'Cluster 2\nSquamous-like',
        'Cluster 3\nGastro-intestinal'
    )
)

table_cancer_types <- do.call(
    gridExtra::gtable_combine,
    args = lapply(
        1:3,
        function(i) {
            
            gridExtra::tableGrob(
                table_cancer_types_data[, ..i],
                theme = gridExtra::ttheme_minimal(
                    padding = unit(c(25, 10), 'pt'),
                    core = list(
                        fg_params = list(
                            col = RColorBrewer::brewer.pal(3, 'Dark2')[i]
                        )
                    ),
                    colhead = list(
                        fg_params = list(col = 'white'),
                        bg_params = list(
                            fill = RColorBrewer::brewer.pal(3, 'Dark2')[i],
                            col = 'white'
                        )
                    )
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
                        c(
                            cancer_type,
                            rep(
                                '',
                                max(
                                    table(
                                        pca_data[
                                            thresh_pass == TRUE,
                                            kmeans_cluster_manual
                                        ]
                                    )
                                ) - .N
                            )
                        )
                    ]
                }
            )
        ),
        pca_data[
            thresh_pass == FALSE,
            c(
                cancer_type,
                rep(
                    '',
                    max(
                        table(
                            pca_data[
                                thresh_pass == TRUE,
                                kmeans_cluster_manual
                            ]
                        )
                    ) - .N
                )
            )
        ]
    ),
    c(
        'Cluster 1\nGynaecological',
        'Cluster 2\nSquamous-like',
        'Cluster 3\nGastro-intestinal',
        'Intermediates'
    )
)

table_cancer_types_all <- do.call(
    gridExtra::gtable_combine,
    args = lapply(
        1:4,
        function(i) {
            
            gridExtra::tableGrob(
                table_cancer_types_data_all[, ..i],
                theme = gridExtra::ttheme_minimal(
                    padding = unit(c(25, 10), 'pt'),
                    core = list(
                        fg_params = list(
                            col = c(RColorBrewer::brewer.pal(3, 'Dark2'), 'darkgrey')[i]
                        )
                    ),
                    colhead = list(
                        fg_params = list(col = 'white'),
                        bg_params = list(
                            fill = c(RColorBrewer::brewer.pal(3, 'Dark2'), 'darkgrey')[i],
                            col = 'white'
                        )
                    )
                ),
                rows = NULL
            )
            
        }
    )
)

table_genes_data <- setNames(
    as.data.table(kmeans_clust_distinct_genes),
    c(
        'Cluster 1\nGynaecological',
        'Cluster 2\nSquamous-like',
        'Cluster 3\nGastro-intestinal'
    )
)

# Stuff I've commented out is for if you want row numbers.

table_genes <- do.call(
    gridExtra::gtable_combine,
    args = lapply(
        1:3,
        function(i) {
            
            g <- gridExtra::tableGrob(
                table_genes_data[, ..i],
                theme = gridExtra::ttheme_default(
                    core = list(
                        # fg_params = list(
                        #     col = RColorBrewer::brewer.pal(3, 'Dark2')[i]
                        # ),
                        bg_params = list(fill = 'white')
                    ),
                    colhead = list(
                        fg_params = list(col = 'white'),
                        bg_params = list(
                            fill = RColorBrewer::brewer.pal(3, 'Dark2')[i],
                            col = 'black'
                        )
                    )
                ),
                rows = NULL
                # rows = switch((i == 1) + 1, NULL, 1:nrow(table_genes_data))
            )
            
            gtable::gtable_add_grob(
                g,
                grobs = grid::rectGrob(gp = grid::gpar(fill = NA)),
                t = 2,
                b = nrow(g),
                l = 1
                # l = switch((i == 1) + 1, 1, 2)
            )
            
        }
    )
)

# PCA plot with line segments:

pca_segment_plot <- ggplot(pca_data, aes(PC1, PC2)) +
    geom_point(aes(colour = as.character(kmeans_cluster_manual)), size = 2) +
    scale_colour_manual(
        name = 'k-means cluster',
        values = RColorBrewer::brewer.pal(3, 'Set2')
    ) +
    geom_segment(
        aes(
            x = kclust_centre['PC1'],
            y = kclust_centre['PC2'],
            xend = PC1_end,
            yend = PC2_end
        ),
        data = segment_endpoints
    ) +
    theme_test() +
    theme(
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 13)
    )

# Additional volcano plot to fill in white space in figure 4...  The idea is that genes
# which occur in more lists will have a larger sample size and therefore greater chance
# of becoming significant.  The fold change is replaced by average EMT-CAF score.

scores_volcano_plot_data <- scores_data_transformed[
    ,
    .(
        ave_score = rowMeans(.SD),
        signif_val = t.test(as.numeric(.SD))$p.value
    ),
    by = gene
]

# Initial scores heatmap for 100 most common EMT/CAF genes:

pdf(
    '../data_and_figures/final_figures/S6.pdf',
    width = 7.5,
    height = 13.6
)

# Need to mess around with the ordering so that it's roughly the same in the clustered
# and SPIN_NH'd heatmaps.

htmp_emt_caf$heatmap + theme(
    axis.text.x = element_text(
        colour = plyr::mapvalues(
            pca_data[
                with(htmp_emt_caf, analyses[ordering_analyses]),
                .(new_clust = switch(thresh_pass + 1, 0, kmeans_cluster_manual)),
                by = cancer_type
                ]$new_clust,
            c(0, 1, 2, 3),
            c('darkgrey', RColorBrewer::brewer.pal(3, 'Dark2'))
        )
    )
)

dev.off()

# PCA plots using segments to indicate directions and transformations/thresholds to
# identify ambiguous cases:

# The colours of the Set2 palette seem like lighter versions of the Dark2 palette, hence
# I'm using them in the below plot to make them look distinct from but similar enough to
# the final groups.

pca_rot_plotlist <- lapply(
    1:3,
    function(i) {
        ggplot(
            pca_data,
            aes(get(paste0('PC1_', i)), get(paste0('PC2_', i)))
        ) +
            geom_point(aes(colour = get(paste0('PC1_', i)) >= cutoffs[i])) +
            scale_colour_manual(
                values = c('lightgrey', RColorBrewer::brewer.pal(3, 'Dark2')[i])
            ) +
            geom_text(
                aes(x = x, y = y, label = label),
                data = data.frame(x = cutoffs[i], y = -2.5, label = paste('Cluster', i)),
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
                x = 'PC1 transformed',
                y = 'PC2 transformed',
                title = paste('k-means cluster', i)
            ) +
            lims(x = c(-2.3, 3.5), y = c(-2.6, 3.1)) +
            guides(colour = FALSE)
    }
)

aligned_plots_pca <- align_plots(
    plot_grid(
        pca_segment_plot + guides(colour = FALSE),
        get_legend(pca_segment_plot),
        nrow = 1,
        ncol = 2,
        rel_widths = c(3, 1)
    ),
    plot_grid(
        plotlist = pca_rot_plotlist,
        nrow = 1,
        ncol = 3,
        align = 'h',
        rel_widths = c(1.125, 1, 1)
    ),
    align = 'v',
    axis = 'l'
)

# tempgrob <- ggplotGrob(temp[[1]])

tempgrob <- ggplotGrob(pca_rot_plotlist[[1]])

# In the following, we're subtracting from 1 npc (which I think is basically the width of
# the viewport, i.e. of the entire plot) the widths associated with left and right margins
# of the PCA rotation scatterplots (by manual inspection, these are 1, 3, 4 and 9) and an
# extra 1 null, which was found by inspection of aligned_plots_pca[[2]]$widths (which only
# has one nonzero element).  We then divide this whole lot by 4, i.e. into the 4 columns
# of the table.

for(i in 1:4) {
    table_cancer_types_all$widths[i] <- 0.25*(unit(1, 'npc') -
                                                  sum(tempgrob$widths[c(1, 3, 4, 9)]) -
                                                  unit(1, 'null'))
}

table_cancer_types_all_width = sum(table_cancer_types_all$widths)
table_cancer_types_all_height = sum(table_cancer_types_all$heights)
table_cancer_types_all$vp <- viewport(
    x = unit(1, 'npc') - 0.5*table_cancer_types_all_width - unit(5.5, 'pt')
)

pdf(
    '../data_and_figures/final_figures/S7.pdf',
    width = 9,
    height = 15
)

plot_grid(
    blank_plot(),
    aligned_plots_pca[[1]],
    aligned_plots_pca[[2]],
    table_cancer_types_all,
    nrow = 4,
    ncol = 1,
    rel_heights = c(0.5, 5.5, 4, 3.3)
) +
    draw_label('A', x = 0, y = 0.975, hjust = -0.1, vjust = 1.1, size = 20, fontface = 2) +
    draw_label('B', x = 0, y = 0.525, hjust = -0.1, vjust = 1.1, size = 20, fontface = 2) +
    draw_label('C', x = 0, y = 0.26, hjust = -0.1, vjust = 1.1, size = 20, fontface = 2)

dev.off()

# Combined scores figure:

pdf(
    '../data_and_figures/final_figures/4.pdf',
    width = 11.5,
    height = 11
)

plot_grid(
    blank_plot(),
    blank_plot(),
    plot_grid(
        ggplot(scores_volcano_plot_data[signif_val < 1], aes(x = ave_score, y = -log10(signif_val))) +
            scale_x_continuous(limits = c(-0.75, 0.55)) +
            scale_y_continuous(expand = c(0, 0), limits = c(-1, 18)) +
            geom_vline(xintercept = 0, linetype = 'dashed', colour = 'grey', size = 0.25) +
            geom_point(colour = 'grey40', alpha = 0.75) +
            geom_point(
                data = scores_volcano_plot_data[gene %in% c('SNAI1', 'SNAI2', 'TWIST1', 'VIM', 'ZEB1', 'ZEB2')],
                size = 2.5,
                colour = RColorBrewer::brewer.pal(4, 'Dark2')[4]
            ) +
            ggrepel::geom_text_repel( # Just the CAF labels
                aes(label = gene),
                data = scores_volcano_plot_data[
                    gene %in% c(
                        'DPT',
                        'ACTA2',
                        'BGN',
                        'DCN',
                        'COL1A1',
                        'COMP',
                        'VCAN',
                        'MMP2',
                        'COL3A1'
                    )
                ],
                nudge_x = 0.17,
                size = 3,
                point.padding = 0.2,
                segment.size = 0.3
            ) +
            ggrepel::geom_text_repel( # Just the pEMT labels
                aes(label = gene),
                data = scores_volcano_plot_data[
                    gene %in% c(
                        'LAMC1',
                        'PVR',
                        'ITGA2',
                        'ITGB1',
                        'LAMC2',
                        'LAMA3',
                        'CALU'
                    )
                ],
                nudge_x = 0.12,
                size = 3,
                point.padding = 0.2,
                segment.size = 0.3
            ) +
            ggrepel::geom_text_repel(
                aes(label = gene),
                data = scores_volcano_plot_data[gene =='CD44'],
                nudge_x = 0.1,
                nudge_y = -1,
                size = 3,
                point.padding = 0.1,
                segment.size = 0.3
            ) +
            ggrepel::geom_text_repel(
                aes(label = gene),
                data = scores_volcano_plot_data[
                    gene %in% c(
                        'VEGFA',
                        'COLGALT1',
                        'TPM4',
                        'DST'
                    )
                ],
                nudge_x = -0.15,
                size = 3,
                point.padding = 0.2,
                segment.size = 0.3
            ) +
            ggrepel::geom_text_repel(
                aes(label = gene),
                data = scores_volcano_plot_data[
                    gene %in% c(
                        'ASPN',
                        'COL1A2',
                        'LUM',
                        'TAGLN',
                        'THY1',
                        'PCOLCE'
                    )
                ],
                nudge_x = -0.14,
                size = 3,
                point.padding = 0.2,
                segment.size = 0.3
            ) +
            ggrepel::geom_text_repel(
                aes(label = gene),
                data = scores_volcano_plot_data[
                    gene %in% c(
                        'FAP',
                        'POSTN',
                        'COL5A1',
                        'COL6A2',
                        'SPARC',
                        'COL5A2'
                    )
                ],
                nudge_x = -0.175,
                size = 3,
                point.padding = 0.2,
                segment.size = 0.3
            ) +
            ggrepel::geom_label_repel(
                aes(label = gene),
                data = scores_volcano_plot_data[gene == 'SNAI2'],
                nudge_x = 0.13,
                size = 3,
                point.padding = 0.3,
                segment.size = 0.3
            ) +
            ggrepel::geom_label_repel(
                aes(label = gene),
                data = scores_volcano_plot_data[gene %in% c('VIM', 'TWIST1')],
                nudge_x = 0.13,
                nudge_y = 0.7,
                size = 3,
                point.padding = 0.3,
                segment.size = 0.3
            ) +
            ggrepel::geom_label_repel(
                aes(label = gene),
                data = scores_volcano_plot_data[gene %in% c('SNAI1', 'ZEB1', 'ZEB2')],
                nudge_x = -0.17,
                size = 3,
                point.padding = 0.3,
                segment.size = 0.3
            ) +
            geom_segment(
                aes(x = 0.1, xend = 0.5, y = -0.5, yend = -0.5),
                arrow = arrow(ends = 'last', length = unit(5, 'pt')),
                colour = RColorBrewer::brewer.pal(11, 'RdBu')[2]
            ) +
            geom_segment(
                aes(x = -0.6, xend = -0.1, y = -0.5, yend = -0.5),
                arrow = arrow(ends = 'first', length = unit(5, 'pt')),
                colour = RColorBrewer::brewer.pal(11, 'RdBu')[10]
            ) +
            annotate(
                geom = 'text',
                x = 0.5,
                y = 0.1,
                label = 'pEMT',
                colour = RColorBrewer::brewer.pal(11, 'RdBu')[2]
            ) +
            annotate(
                geom = 'text',
                x = -0.6,
                y = 0.1,
                label = 'CAF',
                colour = RColorBrewer::brewer.pal(11, 'RdBu')[10]
            ) +
            theme_test() +
            theme(plot.margin = unit(c(5.5, 20, 20, 5.5), 'pt')) +
            labs(x = 'Average pEMT-CAF score', y = latex2exp::TeX('Significance (-log_{10}(p))')),
        ggplot(pca_data, aes(PC1, PC2)) +
            geom_point(
                aes(
                    colour = plyr::mapvalues(
                        kmeans_cluster_manual,
                        c(1, 2, 3),
                        c(
                            'Cluster 1 - Gynaecological',
                            'Cluster 2 - Squamous-like',
                            'Cluster 3 - Gastro-intestinal'
                        )
                    )
                ),
                size = 2
            ) +
            scale_colour_manual(values = RColorBrewer::brewer.pal(3, 'Dark2')) +
            geom_point(data = pca_data[thresh_pass == FALSE], colour = 'lightgrey', size = 2) +
            guides(shape = FALSE) +
            theme_test() +
            theme(
                legend.position = 'bottom',
                legend.direction = 'vertical',
                legend.title = element_blank(),
                legend.text = element_text(size = 11),
                plot.margin = unit(c(20, 20, 5.5, 5.5), 'pt')
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
                        colours = rev(colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))(50)),
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
                                plyr::mapvalues(
                                    kmeans_cluster_manual,
                                    c(1, 2, 3),
                                    RColorBrewer::brewer.pal(3, 'Dark2')
                                )
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
    draw_label('A', x = 0, y = 0.975, hjust = -0.1, vjust = 1.1, size = 20, fontface = 2) +
    draw_label('B', x = 0, y = 0.475, hjust = -0.1, vjust = 1.1, size = 20, fontface = 2) +
    draw_label('C', x = 0.5, y = 0.975, hjust = -0.1, vjust = 1.1, size = 20, fontface = 2)

dev.off()





# Control clustering to see whether our 3 clusters arise naturally from the expression
# data and don't really reflect EMT:

ct_cor_mat <- readRDS('../data_and_figures/ct_cor_mat.rds')

ct_cor_mat_clust <- hclust(as.dist(1 - ct_cor_mat), method = 'average')

ct_cor_htmp <- heat_map(
    ct_cor_mat,
    ordering = ct_cor_mat_clust$order,
    colour_limits = c(min(ct_cor_mat), 1),
    colours = grDevices::heat.colors(50),
    axis_text_size = 11
)

# Global similarity between cancer types:

pdf(
    '../data_and_figures/final_figures/S8.pdf',
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
                quote(
                    startsWith(variable, 't3') |
                        startsWith(variable, 't4')
                ),
                quote( # Maybe n1 should be included in here as well...
                    startsWith(variable, 'n2') |
                        startsWith(variable, 'n3')
                ),
                quote(startsWith(variable, 'm1')),
                quote(variable %in% c('g3', 'g4', 'high grade'))
            ),
            wilcox_test_y_expr = list(
                quote(variable == 0),
                quote(variable == 'complete remission/response'),
                quote(variable > quantile(variable, 0.6)),
                quote(variable %in% c('no', 'absent')),
                quote(
                    startsWith(variable, 't0') |
                        startsWith(variable, 't1') |
                        startsWith(variable, 't2')
                ),
                quote(
                    startsWith(variable, 'n0') |
                        startsWith(variable, 'n1')
                ),
                quote(
                    startsWith(variable, 'm0') |
                        startsWith(variable, 'cm0')
                ),
                quote(variable %in% c('gb', 'g1', 'g2', 'low grade'))
            ),
            min_samples = 10
        )[
            ,
            c('nice_variable_name', 'sigval', 'sigval_adj') := .(
                plyr::mapvalues(
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
            # "sign",
            " sign",
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
            x_factor_levels = htmp_emt_caf$analyses[
                htmp_emt_caf$ordering_analyses
                ], # Same ordering of cancer types as in the commonality heatmap
            y_factor_levels = c(
                'Lymph node metastasis',
                'N stage',
                'Lymphovascular invasion',
                'Grade',
                'T stage',
                'Reduced survival',
                'Therapy resistance'
            ),
            colours = c( # PiYG palette with fatter waist...
                '#276419',
                '#4D9221',
                '#7FBC41',
                '#B8E186',
                '#E6F5D0',
                '#F7F7F7',
                '#F7F7F7',
                '#F7F7F7',
                '#F7F7F7',
                '#F7F7F7',
                '#FDE0EF',
                '#F1B6DA',
                '#DE77AE',
                '#C51B7D',
                '#8E0152'
            ),
            limits = c(-3, 3),
            breaks = c(-3, -2, -1, 0, 1, 2, 3),
            labels = c(
                '-3' = '\u2264 -3',
                '-2' = '-2',
                '-1' = '-1',
                '0' = '0',
                '1' = '1',
                '2' = '2',
                '3' = '\u2265 3'
            ),
            x_lab = NULL,
            y_lab = NULL,
            legend_title = sig_lab,
            # legend_title = NULL,
            plot_title = plyr::mapvalues(
                emt_type,
                names(clin_cor),
                c('pEMT signature', 'CAF signature'),
                warn_missing = FALSE
            ),
            grid_lines = 'grey90',
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

# Below is a figure combining four plots showing my "top 4" clinical features, which
# are chosen because they have sufficiently many points and sufficiently many
# significant cases.  Survival and T stage show no significant cases after adjusting
# p values, while LVI shows one but has relatively few points anyway.

scatterplots <- sapply(
    c(
        'Grade',
        'N stage',
        'Reduced survival',
        'Therapy resistance'
    ),
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
                    (
                        abs(sig_emt) > -log10(pval_adj_threshold) |
                            abs(sig_caf) > -log10(pval_adj_threshold)
                    ) + 1,
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
                values = c(
                    'not_significant' = 'dodgerblue4',
                    'significant' = 'darkgoldenrod1'
                ),
                labels = c(
                    'not_significant' = 'Not significant',
                    'significant' = 'Significant'
                )
            ) +
            labs(title = clin_feat, colour = NULL) +
            theme_test() +
            theme(axis.title = element_blank())
        
    },
    simplify = FALSE,
    USE.NAMES = TRUE
)

scatterplots$Grade <- scatterplots$Grade + ggrepel::geom_text_repel(
    aes(label = stringr::str_replace(test_name, ' - ', '\n')),
    data = emt_caf_sig_data[variable_name == 'Grade' & test_name == 'STAD - CIN'],
    point.padding = 0.1, size = 3.5, lineheight = 0.75, segment.size = 0.3,
    nudge_x = -0.1, nudge_y = -1
) + ggrepel::geom_text_repel(
    aes(label = stringr::str_replace(test_name, ' - ', '\n')),
    data = emt_caf_sig_data[variable_name == 'Grade' & test_name == 'ESCA - Adenocarcinoma'],
    point.padding = 0.1, size = 3.5, lineheight = 0.75, segment.size = 0.3,
    nudge_x = -0.2, nudge_y = -1.2#1.2
) + ggrepel::geom_text_repel(
    aes(label = stringr::str_replace(test_name, ' - ', '\n')),
    data = emt_caf_sig_data[variable_name == 'Grade' & test_name %in% c('HNSC - Atypical', 'HNSC - Malignant-Basal')],
    point.padding = 0.1, size = 3.5, lineheight = 0.75, segment.size = 0.3,
    nudge_x = 0.5, nudge_y = 1
) + ggrepel::geom_text_repel(
    aes(label = stringr::str_replace(test_name, ' - ', '\n')),
    data = emt_caf_sig_data[variable_name == 'Grade' & test_name == 'PAAD'],
    point.padding = 0.1, size = 3.5, lineheight = 0.75, segment.size = 0.3,
    nudge_x = 0.5, nudge_y = 0.1
) + ggrepel::geom_text_repel(
    aes(label = stringr::str_replace(test_name, ' - ', '\n')),
    data = emt_caf_sig_data[variable_name == 'Grade' & test_name == 'CESC'],
    point.padding = 0.1, size = 3.5, lineheight = 0.75, segment.size = 0.3,
    nudge_x = -0.5#, nudge_y = 0.1
)

scatterplots$`N stage` <- scatterplots$`N stage` + ggrepel::geom_text_repel(
    aes(label = stringr::str_replace(test_name, ' - ', '\n')),
    data = emt_caf_sig_data[variable_name == 'N stage' & test_name == 'LUAD - Magnoid'],
    point.padding = 0.1, size = 3.5, lineheight = 0.75, segment.size = 0.3,
    nudge_x = -0.75, nudge_y = -0.1
) + ggrepel::geom_text_repel(
    aes(label = stringr::str_replace(test_name, ' - ', '\n')),
    data = emt_caf_sig_data[variable_name == 'N stage' & test_name == 'HNSC - Malignant-Basal'],
    point.padding = 0.1, size = 3.5, lineheight = 0.75, segment.size = 0.3,
    nudge_x = -0.5, nudge_y = -0.5
) + ggrepel::geom_text_repel(
    aes(label = stringr::str_replace(test_name, ' - ', '\n')),
    data = emt_caf_sig_data[variable_name == 'N stage' & test_name == 'READ'],
    point.padding = 0.1, size = 3.5, lineheight = 0.75, segment.size = 0.3,
    nudge_y = 0.5
) + ggrepel::geom_text_repel(
    aes(label = stringr::str_replace(test_name, ' - ', '\n')),
    data = emt_caf_sig_data[variable_name == 'N stage' & test_name == 'BRCA - Luminal A'],
    point.padding = 0.1, size = 3.5, lineheight = 0.75, segment.size = 0.3,
    nudge_x = 0.1, nudge_y = -0.5
)

scatterplots$`Reduced survival` <- scatterplots$`Reduced survival` + ggrepel::geom_text_repel(
    aes(label = stringr::str_replace(test_name, ' - ', '\n')),
    data = emt_caf_sig_data[variable_name == 'Reduced survival' & test_name == 'LIHC'],
    point.padding = 0.1, size = 3.5, lineheight = 0.75, segment.size = 0.3,
    nudge_x = -0.2, nudge_y = -0.1
) + ggrepel::geom_text_repel(
    aes(label = stringr::str_replace(test_name, ' - ', '\n')),
    data = emt_caf_sig_data[variable_name == 'Reduced survival' & test_name %in% c('LUSC - Primitive', 'BRCA - Luminal B')],
    point.padding = 0.1, size = 3.5, lineheight = 0.75, segment.size = 0.3,
    nudge_y = 0.75
)

scatterplots$`Therapy resistance` <- scatterplots$`Therapy resistance` + ggrepel::geom_text_repel(
    aes(label = stringr::str_replace(test_name, ' - ', '\n')),
    data = emt_caf_sig_data[variable_name == 'Therapy resistance' & test_name == 'PAAD'],
    point.padding = 0.1, size = 3.5, lineheight = 0.75, segment.size = 0.3,
    nudge_x = 0.3, nudge_y = -0.1
) + ggrepel::geom_text_repel(
    aes(label = stringr::str_replace(test_name, ' - ', '\n')),
    data = emt_caf_sig_data[variable_name == 'Therapy resistance' & test_name == 'HNSC - Malignant-Basal'],
    point.padding = 0.1, size = 3.5, lineheight = 0.75, segment.size = 0.3,
    nudge_x = 0.6, nudge_y = -0.1
)

# Combined figure of clinical correlation heatmaps and 4 scatterplots:

clinical_heatmap_1 <- clin_cor_heatmaps$deconv_emt +
    theme(
        axis.text.x = element_blank(),
        legend.position = 'none'
    )

clinical_heatmap_2 <- clin_cor_heatmaps$deconv_caf +
    theme(
        axis.text.x = element_text(
            colour = plyr::mapvalues(
                pca_data[
                    with(htmp_emt_caf, analyses[ordering_analyses]),
                    .(new_clust = switch(thresh_pass + 1, 0, kmeans_cluster_manual)),
                    by = cancer_type
                    ]$new_clust,
                c(0, 1, 2, 3),
                c('darkgrey', RColorBrewer::brewer.pal(3, 'Dark2'))
            )
        ),
        legend.position = 'none',
        plot.margin = unit(c(5.5, 5.5, 1, 5.5), 'pt')
    )

clinical_scatterplots <- lapply(
    scatterplots,
    function(g) {g + theme(legend.position = 'none')}
)

aligned_clinical_heatmaps <- align_plots(
    clinical_heatmap_1,
    clinical_heatmap_2,
    align = 'v'
)

aligned_clinical_scatterplots_1_3 <- align_plots(
    aligned_clinical_heatmaps[[1]],
    aligned_clinical_heatmaps[[2]],
    clinical_scatterplots[[1]],
    clinical_scatterplots[[3]],
    axis = 'l',
    align = 'v'
)

aligned_clinical_scatterplots_2_4 <- align_plots(
    aligned_clinical_heatmaps[[1]],
    aligned_clinical_heatmaps[[2]],
    clinical_scatterplots[[2]],
    clinical_scatterplots[[4]],
    axis = 'r',
    align = 'v'
)

cairo_pdf(
    '../data_and_figures/final_figures/6.pdf',
    width = 9.5,
    height = 12
)

plot_grid(
    plot_grid(
        blank_plot(),
        aligned_clinical_scatterplots_1_3[[1]],
        aligned_clinical_scatterplots_1_3[[2]],
        blank_plot(),
        plot_grid(
            aligned_clinical_scatterplots_1_3[[3]],
            aligned_clinical_scatterplots_2_4[[3]],
            aligned_clinical_scatterplots_1_3[[4]],
            aligned_clinical_scatterplots_2_4[[4]],
            nrow = 2,
            ncol = 2,
            rel_widths = c(5.85, 4.15)
        ),
        blank_plot(),
        nrow = 6,
        ncol = 1,
        rel_heights = c(0.4, 2.03, 3.27, 0.4, 6, 0.7)
    ),
    plot_grid(
        get_legend(
            clin_cor_heatmaps$deconv_emt +
                theme(legend.justification = c(0, 0.595)) +
                guides(fill = guide_colourbar(title.position = 'right'))
        ),
        get_legend(
            scatterplots[[1]] +
                theme(legend.justification = c(0.1, 0.55))
        ),
        nrow = 2,
        ncol = 1,
        rel_heights = c(5.9, 6.6)
    ),
    nrow = 1,
    ncol = 2,
    rel_widths = c(9.5, 3)
) +
    draw_label('A', x = 0, y = 0.98, hjust = -0.1, vjust = 1.1, size = 20, fontface = 2) +
    draw_label('B', x = 0, y = 0.535, hjust = -0.1, vjust = 1.1, size = 20, fontface = 2) +
    draw_label(sig_lab_emt, x = 0.1, y = 0.29, angle = 90, size = 11) +
    draw_label(sig_lab_caf, x = 0.45, y = 0.033, size = 11)

dev.off()

# Clinical correlation scatterplot for all features:

pdf(
    '../data_and_figures/final_figures/S9.pdf',
    width = 7,
    height = 4.5
)

ggplot(
    emt_caf_sig_data[variable_name != 'M stage'],
    aes(x = sig_caf, y = sig_emt)
) +
    geom_hline(yintercept = 0, linetype = 'dashed', colour = 'lightgrey') +
    geom_vline(xintercept = 0, linetype = 'dashed', colour = 'lightgrey') +
    geom_point(
        aes(colour = variable_name),
        shape = 17,
        size = 2.5
    ) +
    ggrepel::geom_text_repel(
        aes(label = stringr::str_extract(test_name, '^[A-Z]+')),
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
            RColorBrewer::brewer.pal(8, "Set1")[-6],
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
