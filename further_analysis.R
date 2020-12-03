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

deconv_data <- readRDS('../data_and_figures/deconv_data.rds')
deconv_plots <- readRDS('../data_and_figures/deconv_plots.rds')

# Cancer types I will leave out: cesc, esca_escc, hnsc_atypical, kich, lusc_basal, stad_ebv (this last one has only 30 samples).
# I'm leaving out kich because the "CAF" end doesn't look a bit like CAFs - I'm not sure what it is.  I don't think this one is trustworthy.
# I'm also using paad on its own, instead of splitting it up into classical and basal.  This was mainly motivated by the loss of clinical significance
# after separating the subtypes, but anyway paad basal is identified as an intermediate in the clustering, so gets filtered out.





# Filtering deconvs based on annotation measures and diagnostics (cell type correlations):

# Barplots of agreement with annotations:

# Here I'm including ESCA ESCC, HNSC Atypical and LUSC Basal as examples of cases where the annotations disagree.
ct_to_keep <- c('blca_luminal_infiltrated', 'blca_luminal_papillary', 'blca_basal_squamous', 'brca_luminal_a', 'brca_luminal_b', 'brca_basal_like',
    'brca_her2_enriched', 'coad', 'esca_ac', 'esca_escc', 'hnsc_atypical', 'hnsc_classical', 'hnsc_mesenchymal_basal', 'kirp', 'lihc',
    'luad_proximal_inflammatory', 'luad_proximal_proliferative', 'luad_terminal_respiratory_unit', 'lusc_basal', 'lusc_classical', 'lusc_secretory',
    'ov_differentiated', 'ov_immunoreactive', 'ov_mesenchymal', 'ov_proliferative', 'paad', 'prad', 'read', 'stad_cin', 'stad_gs', 'stad_msi', 'ucec')
nice_names_for_figure <- c('BLCA - Luminal-Infiltrated', 'BLCA - Luminal-Papillary', 'BLCA - Basal-Squamous', 'BRCA - Luminal A', 'BRCA - Luminal B',
    'BRCA - Basal-like', 'BRCA - HER2-enriched', 'COAD', 'ESCA - Adenocarcinoma', 'ESCA - Squamous', 'HNSC - Atypical', 'HNSC - Classical',
    'HNSC - Malignant-Basal', 'KIRP', 'LIHC', 'LUAD - Squamoid', 'LUAD - Magnoid', 'LUAD - Bronchioid', 'LUSC - Basal', 'LUSC - Classical',
    'LUSC - Secretory', 'OV - Differentiated', 'OV - Immunoreactive', 'OV - Mesenchymal', 'OV - Proliferative', 'PAAD', 'PRAD', 'READ', 'STAD - CIN',
    'STAD - GS', 'STAD - MSI', 'UCEC')

annotation_agreement <- lapply(
    ct_to_keep,
    function(ct) {
        annots <- c(ccle = 'ccle_comp_diff', extra = 'extra_data_score')
        annots <- annots[annots %in% names(deconv_data[[ct]])]
        annot_fun <- function(x) {runmean(x, 30)/max(abs(runmean(x, 30)))}
        corr_diff <- with(
            deconv_data[[ct]],
            list(
                cancer_type = ct,
                pur = mean(head(cor_with_purity$scale[ordering], length(genes_filtered)/3)) -
                    mean(tail(cor_with_purity$scale[ordering], length(genes_filtered)/3))
            )
        )
        if(length(annots) > 0) {
            corr_diff <- c(
                corr_diff,
                sapply(
                    annots,
                    function(annot) with(
                        deconv_data[[ct]],
                        mean(head(annot_fun(get(annot)[ordering]), length(genes_filtered)/3)) -
                            mean(tail(annot_fun(get(annot)[ordering]), length(genes_filtered)/3))
                    ),
                    simplify = FALSE,
                    USE.NAMES = TRUE
                )
            )
        }
        return(corr_diff)
    }
)
annotation_agreement <- rbindlist(annotation_agreement, fill = TRUE)[, n_annot := sum(!is.na(as.numeric(.SD))), by = cancer_type]

barplot_annotations_data <- copy(annotation_agreement)[
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
barplot_annotations_data <- melt(barplot_annotations_data, id.vars = c('cancer_type', 'n_annot'), variable.name = 'annot', value.name = 'score')[
    ,
    n_annot := factor(n_annot, levels = c(3, 2, 1))
]
barplot_annotations_data <- barplot_annotations_data[!is.na(score)]
barplot_annotations_data[, cancer_type := mapvalues(cancer_type, ct_to_keep, nice_names_for_figure)]
barplot_annotations_data[
    ,
    cancer_type := factor(cancer_type, levels = .SD[, .(mean_score = mean(score)), by = cancer_type][order(-mean_score), cancer_type])
]

# The following should work with axis.text.x = element_markdown(angle = 55, hjust = 1) in the ggplot() call, but sadly I can't install ggtext (the
# package from where the element_markdown() function comes):
# barplot_annotations_data[
#     ,
#     cancer_type := paste0(
#         "<span style = 'color: ",
#         ifelse(cancer_type %in% c("ESCA - Squamous", "HNSC - Atypical", "LUSC - Basal"), "red", "black"),
#         ";'>",
#         cancer_type,
#         "</span>"
#     )
# ]

barplot_annotations <- ggplot(barplot_annotations_data, aes(x = cancer_type, y = score, group = annot, colour = annot, fill = annot)) +
    geom_col(position = position_dodge2(width = 0.9, preserve = 'single')) +
    geom_hline(yintercept = 0, colour = 'lightgrey') +
    facet_grid(cols = vars(n_annot), space = 'free', scales = 'free_x') +
    theme_half_open() +
    theme(
        strip.text = element_blank(),
        strip.background = element_blank(),
        axis.text.x = element_text(angle = 55, hjust = 1, size = 10),
        axis.title.y = element_text(size = 11),
        plot.margin = unit(c(5.5, 5.5, 70, 15), 'pt'),
        legend.text = element_text(size = 10, margin = margin(r = 15, unit = 'pt')),
        legend.title = element_text(size = 11, margin = margin(r = 10, unit = 'pt')),
        legend.spacing.x = unit(1, 'pt'),
        legend.justification = c(0.98, 0.1)
    ) +
    scale_colour_manual(values = c('pur' = brewer.pal(11, "PuOr")[2], 'ccle' = brewer.pal(11, "PiYG")[3], 'extra' = 'gold2'), guide = FALSE) +
    scale_fill_manual(
        labels = c('pur' = 'Correlation with purity', 'ccle' = 'Tumours vs. cell lines', 'extra' = 'scRNA-seq: CAF vs. cancer'),
        values = c('pur' = brewer.pal(11, "PuOr")[2], 'ccle' = brewer.pal(11, "PiYG")[3], 'extra' = 'gold2')
    ) +
    guides(fill = guide_legend(direction = 'horizontal', title.position = 'left')) +
    labs(x = NULL, y = 'Between-cluster difference', fill = 'Annotation type:')

# Plot of regression slopes, showing cut-off, and heatmap of regression slopes for all cell types:

# Remove cancer types that didn't pass the annotations check, namely ESCA ESCC, HNSC Atypical and LUSC Basal:
ct_to_keep <- c('blca_luminal_infiltrated', 'blca_luminal_papillary', 'blca_basal_squamous', 'brca_luminal_a', 'brca_luminal_b', 'brca_basal_like',
    'brca_her2_enriched', 'coad', 'esca_ac', 'hnsc_mesenchymal_basal', 'hnsc_classical', 'kirp', 'lihc', 'luad_proximal_inflammatory',
    'luad_proximal_proliferative', 'luad_terminal_respiratory_unit', 'lusc_classical', 'lusc_secretory', 'ov_differentiated', 'ov_immunoreactive',
    'ov_mesenchymal', 'ov_proliferative', 'paad', 'prad', 'read', 'stad_cin', 'stad_gs', 'stad_msi', 'ucec')
nice_names_for_figure <- c('BLCA - Luminal-Infiltrated', 'BLCA - Luminal-Papillary', 'BLCA - Basal-Squamous', 'BRCA - Luminal A', 'BRCA - Luminal B',
    'BRCA - Basal-like', 'BRCA - HER2-enriched', 'COAD', 'ESCA - Adenocarcinoma', 'HNSC - Malignant-Basal', 'HNSC - Classical', 'KIRP', 'LIHC',
    'LUAD - Squamoid', 'LUAD - Magnoid', 'LUAD - Bronchioid', 'LUSC - Classical', 'LUSC - Secretory', 'OV - Differentiated', 'OV - Immunoreactive',
    'OV - Mesenchymal', 'OV - Proliferative', 'PAAD', 'PRAD', 'READ', 'STAD - CIN', 'STAD - GS', 'STAD - MSI', 'UCEC')

cell_type_lms <- sapply(
    deconv_data[ct_to_keep],
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

colnames(cell_type_lms) <- mapvalues(ct_to_keep, ct_to_keep, nice_names_for_figure)
rownames(cell_type_lms) <- gsub(
    '_',
    ' ',
    sapply(rownames(cell_type_lms), function(w) gsub('^[A-Z]|^[a-z]', toupper(str_extract(w, '^[A-Z]|^[a-z]')), w))
)

barplot_min_slope <- ggplot(
    data.table(index = 1:ncol(cell_type_lms), slope = sort(apply(cell_type_lms, 2, min)))[
        ,
        pass := switch((slope < -0.1) + 1, 'pass', 'fail'),
        by = index
    ]
) +
    geom_col(aes(index, slope, fill = pass, colour = NULL)) +
    scale_fill_manual(values = c('#E78AC3', '#66C2A5')) +
    scale_x_discrete(expand = c(0, 0)) +
    geom_hline(yintercept = -0.1, colour = '#377EB8', linetype = 'dashed', size = 0.75) +
    labs(y = 'Minimum regression slope', fill = NULL) +
    theme_test() +
    theme(
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10)
    )

heatmap_slopes_data <- melt(
    as.data.table(t(cell_type_lms), keep.rownames = 'cancer_type'),
    id.vars = 'cancer_type',
    variable.name = 'cell_type',
    value.name = 'slope'
)[, cancer_type := factor(cancer_type, levels = names(sort(apply(cell_type_lms, 2, min))))]

heatmap_slopes <- ggplot(heatmap_slopes_data, aes(x = cancer_type, y = cell_type, fill = slope)) +
    geom_raster() +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    scale_fill_gradientn(
        limits = c(-0.4, 0.4),
        breaks = c(-0.4, -0.1, 0.1, 0.4),
        colours = c(colorRampPalette(brewer.pal(11, 'PuOr')[1:4])(15), rep('white', 10), colorRampPalette(brewer.pal(11, 'PuOr')[8:11])(15)),
        oob = squish
    ) +
    theme(
        panel.border = element_rect(size = 0.5, fill = NA),
        axis.text.x = element_text(angle = 55, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 11),
        legend.text = element_text(size = 10)
    ) +
    labs(x = NULL, y = NULL, fill = 'Regression\nslope')

pdf('../data_and_figures/final_figures_resubmission/S9.pdf', width = 9, height = 12)
print(
    plot_grid(
        get_legend(barplot_annotations),
        barplot_annotations + theme(legend.position = 'none'),
        plot_grid(
            barplot_min_slope + theme(axis.title.y = element_text(vjust = -12)),
            heatmap_slopes,
            nrow = 2,
            ncol = 1,
            align = 'v',
            rel_heights = c(2.3, 3.7)
        ),
        nrow = 3,
        ncol = 1,
        rel_heights = c(1, 5.2, 5.8)
    ) +
        draw_label('A', x = 0, y = 0.98, hjust = -0.1, vjust = 1.1, size = 20, fontface = 2) +
        draw_label('B', x = 0, y = 0.52, hjust = -0.1, vjust = 1.1, size = 20, fontface = 2)
)
dev.off()





# Remove cancer types that didn't pass the diagnostics check, namely BLCA Luminal-Infiltrated, KIRP, KIHC, PRAD, STAD GS and OV Mesenchymal:
ct_to_keep <- c('blca_luminal_papillary', 'blca_basal_squamous', 'brca_luminal_a', 'brca_luminal_b', 'brca_basal_like', 'brca_her2_enriched',
    'coad', 'esca_ac', 'hnsc_mesenchymal_basal', 'hnsc_classical', 'luad_proximal_inflammatory', 'luad_proximal_proliferative',
    'luad_terminal_respiratory_unit', 'lusc_classical', 'lusc_secretory', 'ov_differentiated', 'ov_immunoreactive', 'ov_proliferative', 'paad',
    'read', 'stad_cin', 'stad_msi', 'ucec')
nice_names_for_figure <- c('BLCA - Luminal-Papillary', 'BLCA - Basal-Squamous', 'BRCA - Luminal A', 'BRCA - Luminal B', 'BRCA - Basal-like',
    'BRCA - HER2-enriched', 'COAD', 'ESCA - Adenocarcinoma', 'HNSC - Malignant-Basal', 'HNSC - Classical', 'LUAD - Squamoid', 'LUAD - Magnoid',
    'LUAD - Bronchioid', 'LUSC - Classical', 'LUSC - Secretory', 'OV - Differentiated', 'OV - Immunoreactive', 'OV - Proliferative', 'PAAD', 'READ',
    'STAD - CIN', 'STAD - MSI', 'UCEC')





# Scores heatmap:

rank_mat <- deconv_rank(deconv_data[ct_to_keep])
scores_data_transformed <- deconv_scores(
    expression_data,
    deconv_data[ct_to_keep],
    scale_fun = function(x) x/(3*sd(x)),
    scale_fun_margin = 2,
    transform_data = TRUE
)

# Clustering cancer types by pEMT genes:

scores_mat <- set_rownames(as.matrix(scores_data_transformed[, ..ct_to_keep]), scores_data_transformed$gene)
scores_mat[scores_mat > 1] <- 1
scores_mat[scores_mat < -1] <- -1

ct_cor <- cor(scores_mat[names(head(sort(apply(rank_mat, 1, quantile, 0.25)), 100)), ])
ct_hclust <- hclust(as.dist(1 - ct_cor), method = 'average')
ct_gw <- seriate(as.dist(1 - ct_cor), method = 'GW_average')

ct_hclust_heatmap_data <- melt(as.data.table(ct_cor, keep.rownames = 'ct1'), id.vars = 'ct1', variable.name = 'ct2', value.name = 'corr')[
    ,
    c('ct1', 'ct2') := .(
        mapvalues(factor(ct1, levels = with(ct_hclust, labels[order])), ct_to_keep, nice_names_for_figure),
        mapvalues(factor(ct2, levels = with(ct_hclust, labels[order])), ct_to_keep, nice_names_for_figure)
    )
]

ct_gw_heatmap_data <- melt(as.data.table(ct_cor, keep.rownames = 'ct1'), id.vars = 'ct1', variable.name = 'ct2', value.name = 'corr')[
    ,
    c('ct1', 'ct2') := .(
        mapvalues(factor(ct1, levels = rownames(ct_cor)[get_order(ct_gw)]), ct_to_keep, nice_names_for_figure),
        mapvalues(factor(ct2, levels = rownames(ct_cor)[get_order(ct_gw)]), ct_to_keep, nice_names_for_figure)
    )
]

# CARTO diverging colour palettes:
# Fall: c('#3d5941', '#778868', '#b5b991', '#f6edbd', '#edbb8a', '#de8a5a', '#ca562c')
# ArmyRose: c('#798234', '#a3ad62', '#d0d3a2', '#fdfbe4', '#f0c6c3', '#df91a3', '#d46780')
# Tropic: c('#009B9E', '#42B7B9', '#A7D3D4', '#F1F1F1', '#E4C1D9', '#D691C1', '#C75DAB')

ct_hclust_heatmap <- ggplot(ct_hclust_heatmap_data, aes(x = ct1, y = ct2, fill = corr)) +
    geom_raster() +
    scale_fill_gradientn(
        colours = c('#798234', '#a3ad62', '#d0d3a2', '#fdfbe4', '#f0c6c3', '#df91a3', '#d46780'),
        limits = c(-0.6, 0.6),
        breaks = c(-0.6, -0.3, 0, 0.3, 0.6),
        oob = squish
    ) +
    scale_y_discrete(expand = c(0, 0)) +
    scale_x_discrete(expand = c(0, 0)) +
    theme(
        axis.ticks = element_blank(),
        axis.ticks.length = unit(0, 'pt'),
        axis.text.x = element_text(angle = 55, hjust = 1),
        axis.title = element_blank(),
        panel.background = element_blank(),
    ) +
    labs(fill = 'Correlation')

ct_gw_heatmap <- ggplot(ct_gw_heatmap_data, aes(x = ct1, y = ct2, fill = corr)) +
    geom_raster() +
    scale_fill_gradientn(
        colours = c('#798234', '#a3ad62', '#d0d3a2', '#fdfbe4', '#f0c6c3', '#df91a3', '#d46780'),
        limits = c(-0.6, 0.6),
        breaks = c(-0.6, -0.3, 0, 0.3, 0.6),
        oob = squish
    ) +
    scale_y_discrete(expand = c(0, 0)) +
    scale_x_discrete(expand = c(0, 0)) +
    theme(
        axis.ticks = element_blank(),
        axis.ticks.length = unit(0, 'pt'),
        axis.text.x = element_text(angle = 55, hjust = 1),
        axis.title = element_blank(),
        panel.background = element_blank(),
    ) +
    labs(fill = 'Correlation')

cairo_pdf('../data_and_figures/ct_clust.pdf', width = 8, height = 7, onefile = TRUE)

plot_grid(
    plot_grid(
        blank_plot(),
        dendro(ct_hclust, 'bottom') + theme(plot.margin = unit(c(5.5, 0, 0, 0), 'pt')),
        blank_plot(),
        get_y_axis(ct_hclust_heatmap),
        ct_hclust_heatmap + theme(
            legend.position = 'none',
            plot.margin = unit(c(0, 0, 0, 0), 'pt'),
            axis.text.x = element_blank(),
            axis.text.y = element_blank()
        ),
        dendro(ct_hclust, 'left') + theme(plot.margin = unit(c(5.5, 0, 0, 0), 'pt')),
        blank_plot(),
        get_x_axis(ct_hclust_heatmap),
        blank_plot(),
        nrow = 3,
        ncol = 3,
        rel_heights = c(0.6, 4.2, 1.2),
        rel_widths = c(1.5, 4, 0.5)
    ),
    get_legend(ct_hclust_heatmap),
    nrow = 1,
    ncol = 2,
    rel_widths = c(7, 1)
) %>% print

plot_grid(
    plot_grid(
        blank_plot(),
        dendro(ct_gw[[1]], 'bottom') + theme(plot.margin = unit(c(5.5, 0, 0, 0), 'pt')),
        blank_plot(),
        get_y_axis(ct_gw_heatmap),
        ct_gw_heatmap + theme(
            legend.position = 'none',
            plot.margin = unit(c(0, 0, 0, 0), 'pt'),
            axis.text.x = element_blank(),
            axis.text.y = element_blank()
        ),
        dendro(ct_gw[[1]], 'left') + theme(plot.margin = unit(c(5.5, 0, 0, 0), 'pt')),
        blank_plot(),
        get_x_axis(ct_gw_heatmap),
        blank_plot(),
        nrow = 3,
        ncol = 3,
        rel_heights = c(0.6, 4.2, 1.2),
        rel_widths = c(1.5, 4, 0.5)
    ),
    get_legend(ct_gw_heatmap),
    nrow = 1,
    ncol = 2,
    rel_widths = c(7, 1)
) %>% print

dev.off()

ct_hclust_cut <- cutree(ct_hclust, 3)
ct_clust <- data.table(cancer_type = names(ct_hclust_cut), cluster = ct_hclust_cut)[, memb_strength := silhouette_vals(1 - ct_cor, ct_hclust)]

# ct_clust <- data.table(cancer_type = names(ct_hclust_cut), cluster = ct_hclust_cut)[
#     ,
#     memb_strength := unlist(
#         lapply(
#             1:3,
#             function(i) {
#                 in_i <- names(ct_hclust_cut[ct_hclust_cut == i])
#                 not_in_i <- names(ct_hclust_cut[ct_hclust_cut != i])
#                 rowMeans(ct_cor[in_i, in_i]) - rowMeans(ct_cor[in_i, not_in_i])
#             }
#         )
#     )[cancer_type]
# ]

# Make the names of the cluster match "1", "2", "3" in a consistent order:
setkey(ct_clust, cancer_type)
ct_clust_map <- mapvalues(1:3, 1:3, ct_clust[c('brca_luminal_a', 'hnsc_mesenchymal_basal', 'stad_cin'), cluster])
ct_clust[, cluster_manual := mapvalues(cluster, 1:3, ct_clust_map)]

score_diff_table <- data.table(symbol = names(head(sort(apply(rank_mat, 1, quantile, 0.25)), 100)))[
    ,
    (paste0('score_diff_', 1:3)) := lapply(
        1:3,
        function(i) {
            cts <- list(
                in_i = ct_clust[cluster_manual == i & memb_strength > 0.05, cancer_type],
                not_i = ct_clust[cluster_manual != i & memb_strength > 0.05, cancer_type]
            )
            scores_data_transformed[symbol, rowMeans(.SD[, cts$in_i, with = FALSE]) - rowMeans(.SD[, cts$not_i, with = FALSE])]
        }
    )
] %>% setnames('symbol', 'gene')
score_diff_table[, c('which_max', 'which_max_score') := .(which.max(as.numeric(.SD)), max(as.numeric(.SD))), by = gene]
# ct_clust_distinct_genes <- score_diff_table[order(which_max, -which_max_score), .(gene = gene[which_max_score > 0.22]), by = which_max]
ct_clust_distinct_genes <- score_diff_table[order(which_max, -which_max_score), .(gene = gene[1:20]), by = which_max]

scores_heatmap_data <- scores_data_transformed[
    unique(unlist(ct_clust_distinct_genes$gene)),
    c('gene', ct_clust[memb_strength > 0.05, cancer_type]),
    with = FALSE
]
names(scores_heatmap_data) <- mapvalues(names(scores_heatmap_data), ct_to_keep, nice_names_for_figure, warn_missing = FALSE)

scores_heatmap <- deconv_heatmap(
    scores_heatmap_data,
    order_genes_fun = 'seriate',
    order_genes_method = 'GW_ward',
    order_analyses_fun = 'seriate',
    order_analyses_method = 'GW_ward',
    # order_genes_fun = 'hclust',
    # order_genes_method = 'ward.D2',
    # order_analyses_fun = 'hclust',
    # order_analyses_method = 'ward.D2',
    colour_limits = c(-1, 1),
    legend_breaks = c(-1, 0, 1),
    legend_labels = c('-1' = '-1', '0' = '0', '1' = '1'),
    plot_title = NULL,
    legend.text = element_text(size = 9),
    legend.title = element_text(size = 9)
)

scores_heatmap <- c(
    list(
        heatmap = scores_heatmap$heatmap +
            theme(
                axis.text.x = element_text(
                    angle = 55,
                    hjust = 1,
                    colour = ct_clust[
                        with(scores_heatmap, mapvalues(analyses[ordering_analyses], nice_names_for_figure, ct_to_keep, warn_missing = FALSE)),
                        mapvalues(cluster_manual, 1:3, brewer.pal(3, 'Dark2'))
                    ]
                )
            )
    ),
    scores_heatmap[-1]
)

# deconv_heatmap_dendro_plot(
#     scores_heatmap,
#     direction = 'horizontal',
#     title.position = 'top',
#     title.hjust = 0.5,
#     barwidth = unit(50, 'pt'),
#     barheight = unit(7.5, 'pt'),
#     rel_widths = c(9, 1.25),
#     rel_heights = c(0.75, 10)
# )

# I think it's also possible to colour the dendrogram, but might be more work than it's worth.  Here's some info in case we decide to do this:
# https://stackoverflow.com/questions/21474388/colorize-clusters-in-dendogram-with-ggplot2

# Alternative ordering for scores heatmap, by weighted mean within each cluster:

# ct_clust_distinct_scores <- scores_data_transformed[
#     ct_clust_distinct_genes$gene,
#     c('gene', ct_clust[memb_strength > 0.05, cancer_type]),
#     with = FALSE
# ]
# ct_clust_distinct_scores <- melt(ct_clust_distinct_scores, variable.name = 'cancer_type', value.name = 'score')
# setkey(ct_clust_distinct_scores, gene)
# setkey(score_diff_table, gene)
#
# ct_clust_cancer_types_scores <- rbindlist(
#     lapply(
#         1:3,
#         function(i) {
#             genes <- ct_clust_distinct_genes[which_max == i, gene]
#             ct_clust_distinct_scores[
#                 cancer_type %in% ct_clust[cluster_manual == i, cancer_type] & gene %in% genes,
#                 .(weighted_mean = weighted.mean(.SD[genes, score], score_diff_table[genes, which_max_score]), cluster_manual = i),
#                 by = cancer_type
#             ]
#         }
#     )
# )
#
# ct_clust_distinct_scores[
#     ,
#     c('gene', 'cancer_type') := .(
#         factor(gene, levels = ct_clust_distinct_genes$gene),
#         # factor(cancer_type, levels = with(tempcor_clust, labels[order]))
#         factor(cancer_type, levels = ct_clust_cancer_types_scores[order(cluster_manual, -weighted_mean), cancer_type])
#     )
# ]
#
# ggplot(ct_clust_distinct_scores, aes(x = cancer_type, y = gene, fill = score)) +
#     geom_raster() +
#     scale_x_discrete(expand = c(0, 0)) +
#     scale_y_discrete(expand = c(0, 0)) +
#     scale_fill_gradientn(colours = rev(colorRampPalette(brewer.pal(11, "RdBu"))(50)), limits = c(-1, 1), oob = squish) +
#     theme(axis.text.x = element_text(angle = 55, hjust = 1))

# Alternative ordering by seriate() within each cluster:

# orderings <- lapply(
#     1:3,
#     function(i) list(
#         ordering_genes = get_order(
#             seriate(
#                 # dist(scores_data_transformed[ct_clust_distinct_genes[which_max == i, gene], ct_clust[cluster == i, cancer_type], with = FALSE]),
#                 dist(scores_data_transformed[ct_clust_distinct_genes[which_max == i, gene], -'gene']),
#                 method = 'SPIN_STS'
#             )
#         ),
#         ordering_cancer_type = get_order(
#             seriate(
#                 # dist(
#                 #     t(scores_data_transformed[ct_clust_distinct_genes[which_max == i, gene], ct_clust[cluster == i, cancer_type], with = FALSE])
#                 # ),
#                 dist(t(scores_data_transformed[, ct_clust[cluster_manual == i, cancer_type], with = FALSE])),
#                 method = 'SPIN_STS'
#             )
#         )
#     )
# )
#
# ct_clust_distinct_scores[
#     ,
#     c('gene', 'cancer_type') := .(
#         factor(
#             gene,
#             levels = unlist(
#                 lapply(1:3, function(i) ct_clust_distinct_genes[which_max == i, gene][orderings[[i]]$ordering_genes]),
#                 recursive = FALSE
#             )
#         ),
#         factor(
#             cancer_type,
#             levels = unlist(
#                 lapply(1:3, function(i) ct_clust[cluster_manual == i, cancer_type][orderings[[i]]$ordering_cancer_type]),
#                 recursive = FALSE
#             )
#         )
#     )
# ]
#
# ggplot(ct_clust_distinct_scores, aes(x = cancer_type, y = gene, fill = score)) +
#     geom_raster() +
#     scale_x_discrete(expand = c(0, 0)) +
#     scale_y_discrete(expand = c(0, 0)) +
#     scale_fill_gradientn(colours = rev(colorRampPalette(brewer.pal(11, "RdBu"))(50)), limits = c(-1, 1), oob = squish) +
#     theme(axis.text.x = element_text(angle = 55, hjust = 1))





# Heatmap of top 50 pEMT and CAF genes, ordering genes by SPIN_NH and cancer types by ct_hclust:
set.seed(44398)
htmp_emt_caf <- deconv_heatmap(
    setNames(
        scores_data_transformed[
            unique(c(names(head(sort(apply(rank_mat, 1, quantile, 0.25)), 50)), names(tail(sort(apply(rank_mat, 1, quantile, 0.75)), 50)))),
            c('gene', ..ct_to_keep)
        ],
        c('gene', nice_names_for_figure)
    ),
    order_genes_fun = 'seriate',
    order_genes_method = 'SPIN_NH',
    order_analyses_fun = function(x) ct_hclust$order,
    # order_analyses_fun = function(x) rev(get_order(seriate(dist(t(x)), method = 'SPIN_NH'))),
    # order_analyses_fun = 'seriate',
    # order_analyses_method = 'SPIN_NH',
    plot_title = 'Common EMT and stroma genes'
)

pdf('../data_and_figures/final_figures_resubmission/S11.pdf', width = 7.5, height = 14)
htmp_emt_caf$heatmap + theme(
    axis.text.x = element_text(
        colour = mapvalues(
            ct_clust[
                with(htmp_emt_caf, mapvalues(analyses[ordering_analyses], nice_names_for_figure, ct_to_keep)),
                .(new_clust = switch((memb_strength > 0.05) + 1, 0L, cluster_manual)),
                by = cancer_type
            ]$new_clust,
            0:3,
            c('darkgrey', brewer.pal(3, 'Dark2'))
        )
    )
)
dev.off()

# I think it's also possible to colour the dendrogram, but might be more work than it's worth.  Here's some info in case we decide to do this:
# https://stackoverflow.com/questions/21474388/colorize-clusters-in-dendogram-with-ggplot2





# Tables of cancer types in each cluster and genes scoring highly in each cluster relative to the others:

# Details on making and customising tables can be found in this vignette: https://cran.r-project.org/web/packages/gridExtra/vignettes/tableGrob.html

table_cancer_types_data <- setNames(
    as.data.table(
        lapply(
            1:3,
            function(i) {
                ct_clust[
                    memb_strength > 0.05 & cluster_manual == i,
                    c(
                        mapvalues(cancer_type, ct_to_keep, nice_names_for_figure, warn_missing = FALSE),
                        rep('', max(table(ct_clust[memb_strength > 0.05, cluster_manual])) - .N)
                    )
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
                    ct_clust[
                        memb_strength > 0.05 & cluster_manual == i,
                        c(
                            mapvalues(cancer_type, ct_to_keep, nice_names_for_figure, warn_missing = FALSE),
                            rep('', max(table(ct_clust[memb_strength > 0.05, cluster_manual])) - .N)
                        )
                    ]
                }
            )
        ),
        ct_clust[
            memb_strength <= 0.05,
            c(
                mapvalues(cancer_type, ct_to_keep, nice_names_for_figure, warn_missing = FALSE),
                rep('', max(table(ct_clust[memb_strength > 0.05, cluster_manual])) - .N)
            )
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
    # lapply(
    #     1:3,
    #     function(i) ct_clust_distinct_genes[, c(.SD[which_max == i, gene], rep('', max(.SD[, .N, by = which_max]$N) - .SD[which_max == i, .N]))]
    # )
    as.data.table(lapply(1:3, function(i) score_diff_table[which_max == i][order(-which_max_score), gene[1:20]])),
    c('Cluster 1\nGynaecological', 'Cluster 2\nSquamous-like', 'Cluster 3\nGastro-intestinal')
)

table_genes <- do.call(
    gtable_combine,
    args = lapply(
        1:3,
        function(i) {
            g <- tableGrob(
                table_genes_data[, ..i],
                theme = ttheme_default(
                    core = list(bg_params = list(fill = 'white')),
                    colhead = list(fg_params = list(col = 'white'), bg_params = list(fill = brewer.pal(3, 'Dark2')[i], col = 'black'))
                ),
                rows = NULL
            )
            gtable_add_grob(g, grobs = rectGrob(gp = gpar(fill = NA)), t = 2, b = nrow(g), l = 1) # l = switch((i == 1) + 1, 1, 2)
        }
    )
)

pdf('../data_and_figures/scores_tables.pdf', width = 10, height = 8)
ggdraw(table_cancer_types)
ggdraw(table_cancer_types_all)
ggdraw(table_genes)
dev.off()





# Volcano plot for figure 4: the idea is that genes which occur in more lists will have a larger sample size and therefore greater chance of becoming
# significant.  The fold change is replaced by average EMT-CAF score.

scores_volcano_plot_data <- scores_data_transformed[, .(ave_score = rowMeans(.SD), signif_val = t.test(as.numeric(.SD))$p.value), by = gene]

pdf('../data_and_figures/scores_volcano.pdf', width = 5, height = 5)

# ggplot(scores_volcano_plot_data[signif_val < 1], aes(x = ave_score, y = -log10(signif_val))) +
#     geom_vline(xintercept = 0, linetype = 'dashed', colour = 'grey', size = 0.25) +
#     geom_point(colour = 'grey40', alpha = 0.75) +
#     geom_text_repel(
#         aes(label = gene),
#         data = scores_volcano_plot_data[
#             abs(ave_score) > 0.2 & -log10(signif_val) > 5 | gene %in% c('SNAI1', 'SNAI2', 'TWIST1', 'VIM', 'ZEB1', 'ZEB2')
#         ]
#     )

# In the following, I'm manually choosing which labels to nudge left and which to nudge right.  This is because ggrepel doesn't do it adequately by
# itself, but it has the disadvantage that ggrepel doesn't know to avoid overlaps between the labels arising from the two separate calls to
# geom_label_repel().  So, if we're not careful, the ones that are nudged left can overlap those that are nudged right.

scores_volcano_plot <- ggplot(scores_volcano_plot_data[signif_val < 1], aes(x = ave_score, y = -log10(signif_val))) +
    scale_x_continuous(limits = c(-0.9, 0.6)) +
    scale_y_continuous(expand = c(0, 0), limits = c(-1, 14)) +
    geom_vline(xintercept = 0, linetype = 'dashed', colour = 'grey', size = 0.25) +
    geom_point(colour = 'grey40', alpha = 0.75) +
    geom_point(
        data = scores_volcano_plot_data[gene %in% c('SNAI1', 'SNAI2', 'TWIST1', 'VIM', 'ZEB1', 'ZEB2')],
        size = 2.5,
        colour = brewer.pal(4, 'Dark2')[4]
    ) +
    geom_text_repel(
        aes(label = gene),
        data = scores_volcano_plot_data[gene %in% c('ITGA11', 'OLFML2B', 'ASPN', 'LUM', 'ACTA2', 'TAGLN', 'POSTN')],
        nudge_x = 0.18,
        size = 3,
        point.padding = 0.2,
        segment.size = 0.3
    ) +
    geom_text_repel(
        aes(label = gene),
        data = scores_volcano_plot_data[gene %in% c('SPARC', 'COL1A2', 'COL3A1', 'FAP', 'COL1A1', 'THY1', 'VCAN', 'MMP2', 'ECM2', 'PDGFRB')],
        nudge_x = -0.18,
        size = 3,
        point.padding = 0.2,
        segment.size = 0.3
    ) +
    geom_text_repel(
        aes(label = gene),
        data = scores_volcano_plot_data[gene %in% c('LUZP1', 'LAMC1', 'LAMC2', 'LAMA3', 'ITGA2', 'VCL')],
        nudge_x = 0.12,
        size = 3,
        point.padding = 0.2,
        segment.size = 0.3
    ) +
    geom_text_repel(
        aes(label = gene),
        data = scores_volcano_plot_data[gene %in% c('MCM7', 'CD59')],
        nudge_x = -0.1,
        size = 3,
        point.padding = 0.2,
        segment.size = 0.3
    ) +
    geom_text_repel(
        aes(label = gene),
        data = scores_volcano_plot_data[gene %in% c('MSN', 'PVR', 'ITGB1', 'CD44')],
        nudge_x = -0.15,
        size = 3,
        point.padding = 0.2,
        segment.size = 0.3
    ) +
    geom_label_repel(
        aes(label = gene),
        data = scores_volcano_plot_data[gene %in% c('SNAI1', 'VIM')],
        nudge_x = 0.15,
        nudge_y = 0.5,
        size = 3,
        point.padding = 0.3,
        segment.size = 0.3
    ) +
    geom_label_repel(
        aes(label = gene),
        data = scores_volcano_plot_data[gene %in% c('ZEB2', 'TWIST1')],
        nudge_x = -0.21,
        nudge_y = 0.5,
        size = 3,
        point.padding = 0.3,
        segment.size = 0.3
    ) +
    geom_label_repel(
        aes(label = gene),
        data = scores_volcano_plot_data[gene %in% c('SNAI2', 'ZEB1')],
        nudge_x = -0.2,
        nudge_y = -0.1,
        size = 3,
        point.padding = 0.3,
        segment.size = 0.3
    ) +
    geom_segment(
        aes(x = 0.1, xend = 0.5, y = -0.5, yend = -0.5),
        arrow = arrow(ends = 'last', length = unit(5, 'pt')),
        colour = brewer.pal(11, 'RdBu')[2]
    ) +
    geom_segment(
        aes(x = -0.7, xend = -0.1, y = -0.5, yend = -0.5),
        arrow = arrow(ends = 'first', length = unit(5, 'pt')),
        colour = brewer.pal(11, 'RdBu')[10]
    ) +
    annotate(geom = 'text', x = 0.5, y = 0.1, label = 'pEMT', colour = brewer.pal(11, 'RdBu')[2]) +
    annotate(geom = 'text', x = -0.7, y = 0.1, label = 'CAF', colour = brewer.pal(11, 'RdBu')[10]) +
    theme_test() +
    theme(plot.margin = unit(c(5.5, 20, 20, 5.5), 'pt')) +
    labs(x = 'Average pEMT-CAF score', y = TeX('Significance (-log_{10}(p))'))

dev.off()





top_dendro <- dendro(ct_hclust, 'bottom') + theme(plot.margin = unit(c(20, 0, 0, 0), 'pt'))
right_dendro <- dendro(ct_hclust, 'left') + theme(plot.margin = unit(c(0, 20, 0, 0), 'pt'))
corr_heatmap <- ct_hclust_heatmap +
    theme(
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_text(),
        legend.position = 'none',
        plot.margin = unit(c(0, 0, 5.5, 5.5), 'pt')
    ) +
    labs(x = 'Cancer types', y = 'Cancer types')

aligned_volcano_heatmap <- align_plots(scores_volcano_plot, corr_heatmap + theme(legend.position = 'none'), align = 'v', axis = 'l')
aligned_top_dendro_heatmap <- align_plots(top_dendro, aligned_volcano_heatmap[[2]], align = 'v')
aligned_right_dendro_heatmap <- align_plots(aligned_top_dendro_heatmap[[2]], right_dendro, align = 'h')

# aligned_top_dendro_heatmap[[1]]$heights
# aligned_right_dendro_heatmap[[1]]$heights
# aligned_right_dendro_heatmap[[1]]$widths
# aligned_right_dendro_heatmap[[2]]$widths

# Remove extra spacing added between heatmap and dendrograms by align_plots():
aligned_right_dendro_heatmap[[1]]$widths[9] <- unit(0, 'pt')

pdf('../data_and_figures/final_figures_resubmission/4.pdf', width = 11.5, height = 11)

plot_grid(
    blank_plot(),
    blank_plot(),
    blank_plot(),
    plot_grid(
        plot_grid(
            aligned_volcano_heatmap[[1]],
            plot_grid(
                aligned_top_dendro_heatmap[[1]],
                blank_plot(),
                aligned_right_dendro_heatmap[[1]],
                aligned_right_dendro_heatmap[[2]],
                nrow = 2,
                ncol = 2,
                rel_heights = c(1, 5),
                rel_widths = c(5.1, 0.9)
            ),
            nrow = 2,
            ncol = 1
        ),
        get_legend(
            ct_hclust_heatmap +
                theme(legend.justification = c(0.15, 0.5)) +
                guides(
                    fill = guide_colourbar(
                        direction = 'horizontal',
                        title = 'Correlation\n',
                        title.position = 'right',
                        barwidth = unit(100, 'pt'),
                        barheight = unit(10, 'pt')
                    )
                )
        ),
        nrow = 2,
        ncol = 1,
        rel_heights = c(19, 1)
    ),
    blank_plot(),
    deconv_heatmap_dendro_plot(
        scores_heatmap,
        direction = 'horizontal',
        title.position = 'top',
        title.hjust = 0.5,
        barwidth = unit(50, 'pt'),
        barheight = unit(7.5, 'pt'),
        rel_widths = c(8.75, 1.5),
        rel_heights = c(0.75, 10)
    ),
    nrow = 2,
    ncol = 3,
    rel_widths = c(5, 0.5, 6),
    rel_heights = c(1, 20)
) +
    draw_label('A', x = 0, y = 0.975, hjust = -0.1, vjust = 1.1, size = 20, fontface = 2) +
    draw_label('B', x = 0, y = 0.5, hjust = -0.1, vjust = 1.1, size = 20, fontface = 2) +
    draw_label('C', x = 0.5, y = 0.975, hjust = -0.1, vjust = 1.1, size = 20, fontface = 2)

dev.off()





# Deconv plots without single cell data:

cairo_pdf('../data_and_figures/final_figures_resubmission/S10B.pdf', width = 11.5, height = 10)

print(
	plot_grid(
		blank_plot(),
		plot_grid(
			plotlist = lapply(
				c('blca_luminal_papillary', 'blca_basal_squamous', 'esca_ac'),
				function(ct) {

					deconv_figures <- sapply(
						deconv_plots[[ct]]$plots,
						function(x) x + theme(legend.position = 'none', plot.title = element_blank()),
						simplify = FALSE,
						USE.NAMES = TRUE
					)
					deconv_figures$heatmap <- deconv_figures$heatmap +
						theme(plot.margin = unit(c(0, 20, 5.5, 20), 'pt'), axis.title.x = element_text(), axis.title.y = element_text()) +
                        labs(x = 'Genes', y = 'Genes')

					plot_grid(
						plotlist = c(
							list(
								blank_plot() +
								theme(plot.margin = unit(c(11, 5.5, 5.5, 0), 'pt')) +
								labs(title = nice_names_for_figure[ct_to_keep == ct])
							),
							deconv_figures[c('purity_bar', 'ccle_bar', 'heatmap')]
						),
						nrow = 4,
						ncol = 1,
						align = 'v',
						rel_heights = c(2, 1, 1, 15)
					)

				}
			),
			nrow = 1,
			ncol = 3
		),
		blank_plot(),
		plot_grid(
			plotlist = lapply(
				c('stad_cin', 'stad_msi', 'ucec'),
				function(ct) {

                    deconv_figures <- sapply(
						deconv_plots[[ct]]$plots,
						function(x) x + theme(legend.position = 'none', plot.title = element_blank()),
						simplify = FALSE,
						USE.NAMES = TRUE
					)
                    deconv_figures$heatmap <- deconv_figures$heatmap +
						theme(plot.margin = unit(c(0, 20, 5.5, 20), 'pt'), axis.title.x = element_text(), axis.title.y = element_text()) +
                        labs(x = 'Genes', y = 'Genes')

					plot_grid(
						plotlist = c(
							list(
								blank_plot() +
								theme(plot.margin = unit(c(11, 5.5, 5.5, 0), 'pt')) +
								labs(title = nice_names_for_figure[ct_to_keep == ct])
							),
							deconv_figures[c('purity_bar', 'ccle_bar', 'heatmap')]
						),
						nrow = 4,
						ncol = 1,
						align = 'v',
						rel_heights = c(2, 1, 1, 15)
					)

				}
			),
			nrow = 1,
			ncol = 3
		),
		blank_plot(),
		plot_grid(
			get_legend(
				deconv_plots$brca_luminal_a$plots$purity_bar +
					guides(fill = guide_colourbar(title.position = 'right')) +
					theme(
						legend.justification = c(0, 0.5),
						legend.direction = 'horizontal',
						legend.key.width = NULL,
						legend.key.height = unit(10, 'pt'),
						legend.box.margin = margin(l = 10)
					) +
					labs(fill = 'Correlation with purity\n')
			),
			get_legend(
				deconv_plots$brca_luminal_a$plots$ccle_bar +
					guides(fill = guide_colourbar(title.position = 'right')) +
					theme(
						legend.justification = c(0, 0.5),
						legend.direction = 'horizontal',
						legend.key.width = NULL,
						legend.key.height = unit(10, 'pt'),
						legend.box.margin = margin(l = 10)
					) +
					labs(fill = 'Tumours vs. cell lines\n')
			),
			get_legend(
				deconv_plots$brca_luminal_a$plots$heatmap +
					guides(fill = guide_colourbar(title.position = 'right')) +
					theme(
						legend.justification = c(0, 0.5),
						legend.direction = 'horizontal',
						legend.key.width = unit(25, 'pt'),
						legend.key.height = unit(10, 'pt'),
						legend.box.margin = margin(l = 10)
					) +
					labs(fill = 'Correlation\n')
			),
			nrow = 1,
			ncol = 3
		),
		nrow = 6,
		ncol = 1,
		rel_heights = c(0.15, 1.2, 0.125, 1.2, 0.125, 0.2)
		# rel_heights = c(0.05, 1, 0.05, 1)
	) + draw_label('B', x = 0, y = 0.99, hjust = -0.1, vjust = 1.1, size = 20, fontface = 2)
)

dev.off()





# Subset of plots for main text figure:

# First remake the plots for these cancer types, so we can change the axis labels:

ct_subset <- c('brca_luminal_a', 'hnsc_mesenchymal_basal', 'luad_proximal_proliferative')

heatmap_annotations_subset <- lapply(
    list(
        brca_luminal_a = c('CALU', 'COL1A1', 'FAP', 'ITGAV', 'LAMC2', 'NOTCH2', 'PCOLCE', 'PVR', 'THY1'),
        hnsc_mesenchymal_basal = c('ACTA2', 'CD44', 'COL6A1', 'ITGB1', 'LAMC2', 'POSTN', 'SERPINE1', 'TAGLN', 'TGFBI'),
        luad_proximal_proliferative = c('AREG', 'CD44', 'COL1A2', 'COL8A1', 'ITGB5', 'LAMC2', 'MMP2', 'PDPN', 'TGFBI')
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
                purity_colours = rev(colorRampPalette(brewer.pal(11, "PuOr"))(50)),
                purity_legend_title = 'Correlation with purity\n',
                purity_legend_direction = 'horizontal',
                purity_axis_title = NULL,
                ccle_colours = rev(colorRampPalette(brewer.pal(11, 'PiYG')[c(1:3, 6, 9:11)])(50)),
                ccle_legend_title = 'Tumours vs. cell lines\n',
                ccle_legend_direction = 'horizontal',
                ccle_axis_title = NULL,
                extra_colours = colorRampPalette(c('turquoise4', 'turquoise', 'azure', 'gold2', 'gold4'))(50),
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
        within_clust_corr <- with(deconv_data[[ct]], c(cor_mat[ordering, ordering][1:30, 1:30], cor_mat[rev(ordering), rev(ordering)][1:30, 1:30]))
        within_clust_corr <- mean(within_clust_corr[within_clust_corr != 1])
        between_clust_corr <- with(deconv_data[[ct]], mean(cor_mat[ordering, rev(ordering)][1:30, 1:30]))
        list(cancer_type = ct, wthn = within_clust_corr, btw = between_clust_corr)
    }
) %>% rbindlist

within_between_clust_corr[, n_annot := sum(c('purity_bar', 'ccle_bar', 'extra_bar') %in% names(deconv_plots[[cancer_type]]$plots)), by = cancer_type]

# Make the actual plots:

deconv_summary_scatterplot <- ggplot(
    within_between_clust_corr,
    aes(x = btw, y = wthn, shape = as.character(mapvalues(n_annot, c(1, 2, 3), c('No', 'No', 'Yes'))))
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
    labs(x = 'Between-cluster correlation', y = 'Within-cluster correlation', shape = 'scRNA-seq\ndata\navailable')

pdf('../data_and_figures/deconv_figures_main.pdf', width = 14, height = 9)

plot_grid(
    deconv_plot(
        deconv_plots_subset,
        n_row = 2,
        n_col = 3,
        plots_rel_heights = c(title = 2, purity_bar = 1, ccle_bar = 1, extra_bar = 1, heatmap = 15, axis_labels = 6),
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
                        deconv_plots_subset$`BRCA - Luminal A`$plots[[b]] +
                            theme(legend.justification = 'left') +
                            guides(fill = guide_colourbar(title.position = 'right'))
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





# Control clustering to see whether our 3 clusters arise naturally from the expression
# data and don't really reflect EMT:

gene_variances_top_n <- sapply(
    ct_to_keep,
    function(ct) head(sort(expression_data[deconv_data[[ct]]$sample_ids, apply(.SD, 2, var), .SDcols = -'id'], decreasing = TRUE), 5000),
    simplify = FALSE,
    USE.NAMES = TRUE
)

ct_cor_mat <- sapply(
    ct_to_keep,
    function(ct1) {
        sapply(
            ct_to_keep,
            function(ct2) {
                top_variable_genes <- intersect(names(gene_variances_top_n[[ct1]]), names(gene_variances_top_n[[ct2]]))
                expression_data[
                    ,
                    setNames(
                        lapply(c(ct1, ct2), function(ct) .SD[deconv_data[[ct]]$sample_ids, colMeans(.SD), .SDcols = top_variable_genes]),
                        c('ct1means', 'ct2means')
                    )
                ][, cor(ct1means, ct2means)]
            },
            USE.NAMES = TRUE
        )
    },
    USE.NAMES = TRUE
)

saveRDS(ct_cor_mat, '../data_and_figures/ct_cor_mat.rds')

# ct_cor_mat <- readRDS('../data_and_figures/ct_cor_mat.rds')

ct_cor_mat_clust <- hclust(as.dist(1 - ct_cor_mat), method = 'average')
# ct_cor_mat_clust <- seriate(as.dist(1 - ct_cor_mat), method = 'GW_average')[[1]]

ct_cor_htmp <- heat_map(
    set_colnames(
        set_rownames(ct_cor_mat, mapvalues(rownames(ct_cor_mat), ct_to_keep, nice_names_for_figure)),
        mapvalues(colnames(ct_cor_mat), ct_to_keep, nice_names_for_figure)
    ),
    ordering = ct_cor_mat_clust$order,
    colour_limits = c(min(ct_cor_mat), 1),
    colours = heat.colors(50),
    axis_text_size = 11
)

# The following uses the deconv_heatmap_dendro_plot() function to plot the heatmap with dendrograms.  I should really rename this function, because
# it doesn't really need the output of the deconv_heatmap() function: it only needs a heatmap, clustering objects and dendrograms.

pdf('../data_and_figures/final_figures_resubmission/S12.pdf', width = 10, height = 10)
# pdf('../data_and_figures/control_for_pEMT_groups.pdf', width = 10, height = 10)
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

# The following shows that we don't get really low numbers of variable genes shared between cancer types:
# gplots::heatmap.2(
#     sapply(
#         ct_to_keep,
#         function(ct1) sapply(ct_to_keep, function(ct2) length(intersect(names(gene_variances_top_n[[ct1]]), names(gene_variances_top_n[[ct2]]))))
#     ),
#     trace = 'none'
# )

# I considered an alternative method to more closely mimic the deconvolution, but I'm not sure exactly how it would work.  I wanted to take, for each
# cancer type, two random sets of 20 genes and calculate the correlation of each gene (from other random set, or maybe the set of all genes/most
# variable genes) with one minus that with the other, mimicking my EMT-CAF scores.  I thought about first taking a random set of ~250 genes, ranking
# by correlation with purity and taking the head and tail, but I think this would be too close to my deconvolution and not a true control, since you
# might pick up some EMT-related signal.  But one problem is getting a set of genes that is common to all cancer types, as in the ~70 genes I have in
# my EMT-CAF scores heatmap.  I don't want to pick one random set, because this would be very specific.  But if I average over lots of gene sets, I
# would probably lose all signal.  The same goes for the random sets of ~20 genes for each cancer type.  I think the above clustering of averages of
# variable genes is probably the most effective and least biased control.





# Correlation with clinical features:

clinical_data <- fread('../../TCGA_data/tcga_clinical_data.csv', key = 'id')

emt_types <- list(deconv_emt = ct_to_keep, deconv_caf = ct_to_keep)

clin_cor_genes <- setNames(
    lapply(
        list(head, tail),
        function(FUN) {
            sapply(
                deconv_data[ct_to_keep],
                function(deconv_ct) FUN(deconv_ct$genes_filtered[deconv_ct$ordering], 20),
                simplify = FALSE,
                USE.NAMES = TRUE
            )
        }
    ),
    c('deconv_emt', 'deconv_caf')
)

# If I use the binning procedure from scrabble::score() in the following, it seems as though the significance values all become negative - their
# strengths and their distributions in the heatmaps are basically the same, just negative.

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

clin_cor$deconv_emt[, nice_test_name := mapvalues(test_name, ct_to_keep, nice_names_for_figure)]
clin_cor$deconv_caf[, nice_test_name := mapvalues(test_name, ct_to_keep, nice_names_for_figure)]





# Heatmaps of clinical test significance values:

# Get expressions from these commands:
# latex2exp::TeX('\\textbf{Significance of association}')
# latex2exp::TeX('-log_{10}(p) $\\times$ sign(fold change)')
# Except I deleted a space from the following.  Not sure why I didn't need to do this for the scatterplots.

sig_lab <- expression(
    atop(
        `\textbf{Significance of association}` = paste("", bold(paste("Significance of association"))),
        `-log_{10}(p) $\times$ sign(fold change)` = paste("-log", phantom()[{paste("10")}], "(", "", "p", ")", " ", "", phantom() %*% phantom(),
            "sign",
            # " sign",
            "(", "fold change", ")", "")
    )
)

clin_cor_heatmaps <- sapply(
    names(clin_cor),
    function(emt_type) {
        # We'll ignore M stage, because there are no significant cases.
        clinical_test_heatmap(
            clin_cor[[emt_type]][variable_name != 'pathologic_m'],
            x_var = 'nice_test_name',
            y_var = 'nice_variable_name',
            fill_var = 'sigval',
            hclust_method = NULL,
            x_factor_levels = with(
                ct_hclust,
                mapvalues(labels[order][labels[order] %in% clin_cor[[emt_type]]$test_name], ct_to_keep, nice_names_for_figure)
            ),
            y_factor_levels = c('Lymph node metastasis', 'N stage', 'Lymphovascular invasion', 'Grade', 'T stage', 'Reduced survival',
                'Therapy resistance'),
            colours = c('#276419', '#4D9221', '#7FBC41', '#B8E186', '#E6F5D0', '#F7F7F7', '#F7F7F7', '#F7F7F7', '#F7F7F7', '#F7F7F7', '#FDE0EF',
                '#F1B6DA', '#DE77AE', '#C51B7D', '#8E0152'), # PiYG palette with fatter waist...
            limits = c(-3, 3),
            breaks = c(-3, -2, -1, 0, 1, 2, 3),
            labels = c('-3' = '\u2264 -3', '-2' = '-2', '-1' = '-1', '0' = '0', '1' = '1', '2' = '2', '3' = '\u2265 3'),
            x_lab = NULL,
            y_lab = NULL,
            legend_title = sig_lab,
            # legend_title = NULL,
            plot_title = mapvalues(emt_type, names(clin_cor), c('pEMT signature', 'CAF signature'), warn_missing = FALSE),
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

cairo_pdf('../data_and_figures/clinical_heatmaps.pdf', width = 8.5, height = 5)
plot_grid(
    plot_grid(
        clin_cor_heatmaps$deconv_emt + theme(axis.text.x = element_blank(), legend.position = 'none'),
        clin_cor_heatmaps$deconv_caf + theme(axis.text.x = element_blank(), legend.position = 'none', plot.margin = unit(c(5.5, 5.5, 1, 5.5), 'pt')),
        ggdraw(
            get_x_axis(
                clin_cor_heatmaps$deconv_caf + theme(
                    axis.text.x = element_text(
                        colour = mapvalues(
                            ct_clust[
                                with(ct_hclust, labels[order][labels[order] %in% clin_cor$deconv_caf$test_name]),
                                .(new_clust = switch((memb_strength > 0.05) + 1, 0L, cluster_manual)),
                                by = cancer_type
                            ]$new_clust,
                            c(0, 1, 2, 3),
                            c('darkgrey', brewer.pal(3, 'Dark2'))
                        )
                    )
                )
            )
        ),
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
    rel_widths = c(6, 2.5)
)
dev.off()





# Scatterplots of significance of association for EMT vs. that for CAFs:

# latex2exp::TeX('-log_{10}(p_{CAF}) $\\times$ sign(fold change)')
sig_lab_caf <- expression(
    atop(
        `\textbf{Significance of association for CAFs}` = paste("", bold(paste("Significance of association for CAFs"))),
        `-log_{10}(p_{CAF}) $\times$ sign(fold change)` = paste("-log", phantom()[{paste("10")}], "(", "", "p", phantom()[{paste("CAF")}], ")", "",
            " ", "", phantom() %*% phantom(), " sign", "(", "fold change", ")", "")
    )
)

# latex2exp::TeX('-log_{10}(p_{pEMT}) $\\times$ sign(fold change)')
sig_lab_emt <- expression(
    atop(
        `\textbf{Significance of association for pEMT}` = paste("", bold(paste("Significance of association for pEMT"))),
        `-log_{10}(p_{pEMT}) $\times$ sign(fold change)` = paste("-log", phantom()[{paste("10")}], "(", "", "p", phantom()[{paste("pEMT")}], ")", "",
            " ", "", phantom() %*% phantom(), " sign", "(", "fold change", ")", "")
    )
)

# Scatterplot including all clinical features and cancer types:

emt_caf_sig_data <- merge(clin_cor$deconv_emt, clin_cor$deconv_caf, by = c('nice_test_name', 'nice_variable_name'))[
    ,
    .(
        test_name = nice_test_name,
        variable_name = nice_variable_name,
        sig_emt = sigval.x,
        sig_caf = sigval.y,
        sig_adj_emt = sigval_adj.x, # Include adjusted values in case I want them
        sig_adj_caf = sigval_adj.y # (I don't use them here)
    )
]

pval_adj_threshold <- adjust_threshold_bh(
    c(clin_cor$deconv_emt[variable_name != 'pathologic_m', pval], clin_cor$deconv_caf[variable_name != 'pathologic_m', pval])
)

pdf('../data_and_figures/final_figures_resubmission/S13.pdf', width = 7, height = 4.5)
ggplot(emt_caf_sig_data[variable_name != 'M stage'], aes(x = sig_caf, y = sig_emt)) +
    geom_hline(yintercept = 0, linetype = 'dashed', colour = 'lightgrey') +
    geom_vline(xintercept = 0, linetype = 'dashed', colour = 'lightgrey') +
    geom_point(aes(colour = variable_name), shape = 17, size = 2.5) +
    geom_text_repel(
        aes(label = str_extract(test_name, '^[A-Z]+')),
        data = emt_caf_sig_data[
            variable_name != 'M stage' & (
                sig_emt > -log10(pval_adj_threshold) | sig_emt < log10(pval_adj_threshold) | sig_caf > -log10(pval_adj_threshold) |
                    sig_caf < log10(pval_adj_threshold)
            )
        ],
        point.padding = 0.1,
        nudge_x = 0.1
    ) +
    scale_colour_manual(
        values = setNames(
            brewer.pal(8, "Set1")[-6],
            c('Lymph node metastasis', 'Therapy resistance', 'Reduced survival', 'Lymphovascular invasion', 'T stage', 'N stage', 'Grade')
        )
    ) +
    labs(x = sig_lab_caf, y = sig_lab_emt, colour = 'Clinical feature') +
    theme_test()
dev.off()

# Below is a figure combining four plots showing my "top 4" clinical features, which
# are chosen because they have sufficiently many points and sufficiently many
# significant cases.  Survival and T stage show no significant cases after adjusting
# p values, while LVI shows one but has relatively few points anyway.

scatterplots_data <- rbindlist(
    lapply(
        c('Grade', 'Lymph node metastasis', 'Lymphovascular invasion', 'M stage', 'N stage', 'Reduced survival', 'T stage', 'Therapy resistance'),
        function(clin_feat) {
            pthresh <- adjust_threshold_bh(
                c(clin_cor$deconv_emt[nice_variable_name == clin_feat, pval], clin_cor$deconv_caf[nice_variable_name == clin_feat, pval])
            )
            dt <- emt_caf_sig_data[variable_name == clin_feat]
            if(is.na(pthresh)) {
                dt[, sig := 'not_significant']
            } else {
                dt[
                    ,
                    sig := switch((abs(sig_emt) > -log10(pthresh) | abs(sig_caf) > -log10(pthresh)) + 1, 'not_significant', 'significant'),
                    by = test_name
                ]
            }
            return(dt)
        }
    )
)

scatterplots <- sapply(
    c('Grade', 'Lymph node metastasis', 'N stage', 'Therapy resistance'),
    function(clin_feat) {
        ggplot(scatterplots_data[variable_name == clin_feat], aes(x = sig_caf, y = sig_emt)) +
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
    data = scatterplots_data[variable_name == 'Grade' & test_name == 'ESCA - Adenocarcinoma'],
    point.padding = 0.1, size = 3.5, lineheight = 0.75, segment.size = 0.3,
    nudge_y = -0.4
) + geom_text_repel(
    aes(label = str_replace(test_name, ' - ', '\n')),
    data = scatterplots_data[variable_name == 'Grade' & test_name == 'STAD - CIN'],
    point.padding = 0.1, size = 3.5, lineheight = 0.75, segment.size = 0.3,
    nudge_y = -0.1
)

scatterplots$`Lymph node metastasis` <- scatterplots$`Lymph node metastasis` + geom_text_repel(
    aes(label = str_replace(test_name, ' - ', '\n')),
    data = emt_caf_sig_data[variable_name == 'Lymph node metastasis' & test_name == 'COAD'],
    point.padding = 0.1, size = 3.5, lineheight = 0.75, segment.size = 0.3,
    nudge_x = 0.1, nudge_y = -0.1
) + geom_text_repel(
    aes(label = str_replace(test_name, ' - ', '\n')),
    data = emt_caf_sig_data[variable_name == 'Lymph node metastasis' & test_name == 'HNSC - Classical'],
    point.padding = 0.1, size = 3.5, lineheight = 0.75, segment.size = 0.3,
    nudge_x = 0.9
) + geom_text_repel(
    aes(label = str_replace(test_name, ' - ', '\n')),
    data = emt_caf_sig_data[variable_name == 'Lymph node metastasis' & test_name == 'HNSC - Malignant-Basal'],
    point.padding = 0.1, size = 3.5, lineheight = 0.75, segment.size = 0.3,
    nudge_y = -1.2
)

scatterplots$`N stage` <- scatterplots$`N stage` + geom_text_repel(
    aes(label = str_replace(test_name, ' - ', '\n')),
    data = emt_caf_sig_data[variable_name == 'N stage' & test_name == 'HNSC - Classical'],
    point.padding = 0.1, size = 3.5, lineheight = 0.75, segment.size = 0.3,
    nudge_x = 0.5
) + geom_text_repel(
    aes(label = str_replace(test_name, ' - ', '\n')),
    data = emt_caf_sig_data[variable_name == 'N stage' & test_name == 'HNSC - Malignant-Basal'],
    point.padding = 0.1, size = 3.5, lineheight = 0.75, segment.size = 0.3,
    nudge_x = 0.1, nudge_y = -0.3
) + geom_text_repel(
    aes(label = str_replace(test_name, ' - ', '\n')),
    data = emt_caf_sig_data[variable_name == 'N stage' & test_name == 'READ'],
    point.padding = 0.1, size = 3.5, lineheight = 0.75, segment.size = 0.3,
    nudge_y = 0.1
)

scatterplots$`Therapy resistance` <- scatterplots$`Therapy resistance` + geom_text_repel(
    aes(label = str_replace(test_name, ' - ', '\n')),
    data = emt_caf_sig_data[variable_name == 'Therapy resistance' & test_name %in% c('HNSC - Malignant-Basal', 'PAAD')],
    point.padding = 0.1, size = 3.5, lineheight = 0.75, segment.size = 0.3,
    nudge_x = 0.1, nudge_y = -0.1
)

pdf('../data_and_figures/clinical_4_scatterplots.pdf', width = 8, height = 6.5)
plot_grid(
    textGrob(sig_lab_emt, rot = 90, gp = gpar(fontsize = 9)),
    plot_grid(plotlist = lapply(scatterplots, function(g) {g + theme(legend.position = 'none')}), nrow = 2, ncol = 2),
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





# Combined figure of clinical correlation heatmaps and 4 scatterplots:

clinical_heatmap_1 <- clin_cor_heatmaps$deconv_emt + theme(axis.text.x = element_blank(), legend.position = 'none')

clinical_heatmap_2 <- clin_cor_heatmaps$deconv_caf +
    theme(
        axis.text.x = element_text(
            colour = mapvalues(
                ct_clust[
                    with(ct_hclust, labels[order][labels[order] %in% clin_cor$deconv_caf$test_name]),
                    .(new_clust = switch((memb_strength > 0.05) + 1, 0L, cluster_manual)),
                    by = cancer_type
                ]$new_clust,
                c(0, 1, 2, 3),
                c('darkgrey', brewer.pal(3, 'Dark2'))
            )
        ),
        legend.position = 'none',
        plot.margin = unit(c(5.5, 5.5, 1, 5.5), 'pt')
    )

clinical_scatterplots <- lapply(scatterplots, function(g) {g + theme(legend.position = 'none')})

aligned_clinical_heatmaps <- align_plots(clinical_heatmap_1, clinical_heatmap_2, align = 'v')

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

cairo_pdf('../data_and_figures/final_figures_resubmission/5.pdf', width = 8.5, height = 11.5)
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
            rel_widths = c(6.02, 3.98)
        ),
        blank_plot(),
        nrow = 6,
        ncol = 1,
        rel_heights = c(0.4, 2.03, 3.27, 0.4, 5.5, 0.7)
    ),
    plot_grid(
        get_legend(
            clin_cor_heatmaps$deconv_emt + theme(legend.justification = c(0, 0.595)) + guides(fill = guide_colourbar(title.position = 'right'))
        ),
        get_legend(scatterplots[[1]] + theme(legend.justification = c(0.1, 0.55))),
        nrow = 2,
        ncol = 1,
        rel_heights = c(5.9, 6.6)
    ),
    nrow = 1,
    ncol = 2,
    rel_widths = c(6, 2.5)
) +
    draw_label('A', x = 0, y = 0.98, hjust = -0.1, vjust = 1.1, size = 20, fontface = 2) +
    draw_label('B', x = 0, y = 0.52, hjust = -0.1, vjust = 1.1, size = 20, fontface = 2) +
    draw_label(sig_lab_emt, x = 0.1, y = 0.29, angle = 90, size = 11) +
    draw_label(sig_lab_caf, x = 0.45, y = 0.033, size = 11)
dev.off()





# To investigate/check the strong negative association of pEMT with LN mets in HNSC Classical:

# Check the subtype assignments match:
hnsc_subtypes_data <- as.data.table(readxl::read_xlsx('../../TCGA_data/HNSC/7.2.xlsx'))[
    ,
    .(
        individual_id = as.character(Barcode),
        subtype_ref = 'doi:10.1038/nature14129',
        subtype = as.character(RNA)
    )
]
sum(hnsc_subtypes_data[subtype == 'Classical', individual_id] %in% meta_data[deconv_data$hnsc_classical$sample_ids, patient_id])
sum(meta_data[deconv_data$hnsc_classical$sample_ids, patient_id] %in% hnsc_subtypes_data[subtype == 'Classical', individual_id])
# All looks fine!

# Check differences in pEMT gene expression manually using boxplots, including swapping the pEMT signatures:

hnsc_lnmets_data <- rbindlist(
    lapply(
        c('hnsc_classical', 'hnsc_mesenchymal_basal'),
        function(ct) rbind(
            data.table(
                cancer_type = ct,
                lnmets = 'Zero',
                individual_id = clinical_data[meta_data[deconv_data[[ct]]$sample_ids, patient_id]][number_of_lymphnodes_positive_by_he == 0, id]
            ),
            data.table(
                cancer_type = ct,
                lnmets = 'Multiple',
                individual_id = clinical_data[meta_data[deconv_data[[ct]]$sample_ids, patient_id]][number_of_lymphnodes_positive_by_he > 1, id]
            )
        )
    )
)

hnsc_lnmets_data[
    ,
    pemt_expression := expression_data[
        meta_data[patient_id %in% individual_id, id],
        rowMeans(.SD),
        .SDcols = with(deconv_data[[cancer_type]], head(genes_filtered[ordering], 20))
    ],
    by = .(cancer_type, lnmets)
][, pval := wilcox.test(pemt_expression[lnmets == 'Zero'], pemt_expression[lnmets == 'Multiple'])$p.value, by = cancer_type][
    ,
    c('lnmets', 'cancer_type') := .(
        factor(lnmets, levels = c('Zero', 'Multiple')),
        mapvalues(cancer_type, c('hnsc_classical', 'hnsc_mesenchymal_basal'), c('HNSC - Classical', 'HNSC - Malignant-Basal'))
    )
]

pemt_boxplot <- ggplot(hnsc_lnmets_data, aes(x = lnmets, y = pemt_expression)) +
    facet_grid(cols = vars(cancer_type)) +
    geom_boxplot() +
    geom_jitter(width = 0.2) +
    geom_text(
        aes(label = pval),
        data = data.frame(
            lnmets = c(1.5, 1.5),
            pemt_expression = c(5, 5),
            cancer_type = c('HNSC - Classical', 'HNSC - Malignant-Basal'),
            pval = hnsc_lnmets_data[order(cancer_type), paste('p =', sprintf("%.5f", unique(pval)))]
        )
    ) +
    theme_test() +
    labs(x = 'LN metastases', y = 'Average pEMT gene expression', title = 'Respective signatures')

hnsc_lnmets_data[
    ,
    pemt_expression := expression_data[
        meta_data[patient_id %in% individual_id, id],
        rowMeans(.SD),
        .SDcols = with(
            deconv_data[[c('hnsc_classical', 'hnsc_mesenchymal_basal')[c('HNSC - Classical', 'HNSC - Malignant-Basal') != unique(cancer_type)]]],
            head(genes_filtered[ordering], 20)
        )
    ],
    by = .(cancer_type, lnmets)
][, pval := wilcox.test(pemt_expression[lnmets == 'Zero'], pemt_expression[lnmets == 'Multiple'])$p.value, by = cancer_type]

pemt_boxplot_swapped <- ggplot(hnsc_lnmets_data, aes(x = lnmets, y = pemt_expression)) +
    facet_grid(cols = vars(cancer_type)) +
    geom_boxplot() +
    geom_jitter(width = 0.2) +
    geom_text(
        aes(label = pval),
        data = data.frame(
            lnmets = c(1.5, 1.5),
            pemt_expression = c(5, 5),
            cancer_type = c('HNSC - Classical', 'HNSC - Malignant-Basal'),
            pval = hnsc_lnmets_data[order(cancer_type), paste('p =', sprintf("%.5f", unique(pval)))]
        )
    ) +
    theme_test() +
    labs(x = 'LN metastases', y = 'Average pEMT gene expression', title = 'Swapped signatures')

plot_grid(pemt_boxplot, pemt_boxplot_swapped, nrow = 2, ncol = 1)

# Check for differences between the subtypes in days to death or number of LN mets:

hnsc_survival_data <- rbindlist(
    lapply(
        c('hnsc_classical', 'hnsc_mesenchymal_basal'),
        function(ct) data.table(
            cancer_type = ct,
            individual_id = meta_data[deconv_data[[ct]]$sample_ids, patient_id],
            days_to_death = clinical_data[meta_data[deconv_data[[ct]]$sample_ids, patient_id], days_to_death],
            lnmets = clinical_data[meta_data[deconv_data[[ct]]$sample_ids, patient_id], number_of_lymphnodes_positive_by_he]
        )
    )
)[, cancer_type := mapvalues(cancer_type, c('hnsc_classical', 'hnsc_mesenchymal_basal'), c('HNSC - Classical', 'HNSC - Malignant-Basal'))]

ggplot(hnsc_survival_data, aes(x = cancer_type, y = days_to_death)) +
    geom_boxplot() +
    geom_jitter(width = 0.2) +
    ylim(0, 1750) +
    theme_test() +
    labs(x = 'Cancer type', y = 'Days to death')

wilcox.test(
    hnsc_survival_data[!is.na(days_to_death) & cancer_type == 'HNSC - Classical', days_to_death],
    hnsc_survival_data[!is.na(days_to_death) & cancer_type == 'HNSC - Malignant-Basal', days_to_death]
)

ggplot(hnsc_survival_data, aes(x = cancer_type, y = lnmets)) +
    geom_boxplot() +
    geom_jitter(width = 0.2) +
    ylim(0, 10) +
    theme_test() +
    labs(x = 'Cancer type', y = 'Lymph node metastases')

wilcox.test(
    hnsc_survival_data[!is.na(days_to_death) & cancer_type == 'HNSC - Classical', lnmets],
    hnsc_survival_data[!is.na(days_to_death) & cancer_type == 'HNSC - Malignant-Basal', lnmets]
)

# Try with inferred subtypes, both with and without centring:

hnsc_subtypes_data_list <- list(
    centred = fread('../../TCGA_data/tcga_subtypes_data_centred.csv')[subtype_ref == 'doi:10.1038/nature14129'],
    uncentred = fread('../../TCGA_data/tcga_subtypes_data.csv')[subtype_ref == 'doi:10.1038/nature14129']
)

pemt_boxplots_inferred_subtypes <- sapply(
    hnsc_subtypes_data_list,
    function(dt) {

        lnmets_data <- rbindlist(
            lapply(
                list('Classical', c('Mesenchymal', 'Basal')),
                function(ct) rbind(
                    data.table(
                        cancer_type = tolower(paste(ct, collapse = '_')),
                        lnmets = 'Zero',
                        individual_id = clinical_data[dt[inferred_subtype %in% ct, unique(patient_id)]][number_of_lymphnodes_positive_by_he == 0, id]
                    ),
                    data.table(
                        cancer_type = tolower(paste(ct, collapse = '_')),
                        lnmets = 'Multiple',
                        individual_id = clinical_data[dt[inferred_subtype %in% ct, unique(patient_id)]][number_of_lymphnodes_positive_by_he > 1, id]
                    )
                )
            )
        )[
            ,
            pemt_expression := expression_data[
                meta_data[patient_id %in% individual_id & sample_type == 'primary', id],
                rowMeans(.SD),
                .SDcols = with(deconv_data[[paste0('hnsc_', cancer_type)]], head(genes_filtered[ordering], 20))
            ],
            by = .(cancer_type, lnmets)
        ][, pval := wilcox.test(pemt_expression[lnmets == 'Zero'], pemt_expression[lnmets == 'Multiple'])$p.value, by = cancer_type][
            ,
            c('lnmets', 'cancer_type') := .(
                factor(lnmets, levels = c('Zero', 'Multiple')),
                mapvalues(cancer_type, c('classical', 'mesenchymal_basal'), c('HNSC - Classical', 'HNSC - Malignant-Basal'))
            )
        ]

        lnmets_boxplot <- ggplot(lnmets_data, aes(x = lnmets, y = pemt_expression)) +
            facet_grid(cols = vars(cancer_type)) +
            geom_boxplot() +
            geom_jitter(width = 0.2) +
            geom_text(
                aes(label = pval),
                data = data.frame(
                    lnmets = c(1.5, 1.5),
                    pemt_expression = c(4, 4),
                    cancer_type = c('HNSC - Classical', 'HNSC - Malignant-Basal'),
                    pval = lnmets_data[order(cancer_type), paste('p =', sprintf("%.5f", unique(pval)))]
                )
            ) +
            theme_test() +
            labs(x = 'LN metastases', y = 'Average pEMT gene expression')

        return(list(data = lnmets_data, boxplot = lnmets_boxplot))

    },
    simplify = FALSE,
    USE.NAMES = TRUE
)

pdf('../data_and_figures/hnsc_pemt_lnmets_boxplot.pdf')
plot_grid(pemt_boxplot + labs(title = NULL), pemt_boxplots_inferred_subtypes$centred$boxplot, nrow = 2, ncol = 1)
# plot_grid(pemt_boxplot + labs(title = NULL), pemt_boxplots_inferred_subtypes$uncentred$boxplot, nrow = 2, ncol = 1)
dev.off()
