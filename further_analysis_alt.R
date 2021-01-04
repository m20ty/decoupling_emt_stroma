library(data.table) # 1.12.8
library(ggplot2) # 3.3.0
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
library(colorspace) # 1.4.1
library(cluster) # 2.1.0
library(msigdbr) # 7.1.1
library(clusterProfiler) # 3.14.3

source('general_functions.R')
source('tcga_functions.R')

expression_data <- fread('../../TCGA_data/tcga_expression_data.csv', key = 'id')
meta_data <- fread('../../TCGA_data/tcga_meta_data.csv', key = 'id')

expression_data <- expression_data[meta_data[sample_type != 'normal', id]]
meta_data <- meta_data[sample_type != 'normal']

deconv_data <- readRDS('../data_and_figures/deconv_data.rds')
deconv_plots <- readRDS('../data_and_figures/deconv_plots.rds')





# Filtering deconvs based on annotation measures and diagnostics (cell type correlations):

# Barplots of agreement with annotations:

# Here I'm including ESCA ESCC, HNSC Atypical and LUSC Basal as examples of cases where the annotations disagree.
ct_to_keep <- c('blca_luminal_infiltrated', 'blca_luminal_papillary', 'blca_basal_squamous', 'brca_luminal_a', 'brca_luminal_b', 'brca_basal_like',
    'brca_her2_enriched', 'coad', 'esca_ac', 'esca_escc', 'hnsc_atypical', 'hnsc_classical', 'hnsc_mesenchymal_basal', 'kirp', 'lihc',
    'luad_proximal_inflammatory', 'luad_proximal_proliferative', 'luad_terminal_respiratory_unit', 'lusc_basal', 'lusc_classical', 'lusc_secretory',
    'ov_differentiated', 'ov_immunoreactive', 'ov_mesenchymal', 'ov_proliferative', 'paad_basal_moffitt', 'paad_classical_moffitt', 'prad', 'read',
    'stad_cin', 'stad_gs', 'stad_msi', 'ucec')
nice_names_for_figure <- c('BLCA - Luminal-Infiltrated', 'BLCA - Luminal-Papillary', 'BLCA - Basal-Squamous', 'BRCA - Luminal A', 'BRCA - Luminal B',
    'BRCA - Basal-like', 'BRCA - HER2-enriched', 'COAD', 'ESCA - Adenocarcinoma', 'ESCA - Squamous', 'HNSC - Atypical', 'HNSC - Classical',
    'HNSC - Malignant-Basal', 'KIRP', 'LIHC', 'LUAD - Squamoid', 'LUAD - Magnoid', 'LUAD - Bronchioid', 'LUSC - Basal', 'LUSC - Classical',
    'LUSC - Secretory', 'OV - Differentiated', 'OV - Immunoreactive', 'OV - Mesenchymal', 'OV - Proliferative', 'PAAD - Basal', 'PAAD - Classical',
    'PRAD', 'READ', 'STAD - CIN', 'STAD - GS', 'STAD - MSI', 'UCEC')

sc_metadata <- list(
	brca = quote(
		fread('../data_and_figures/qian_breast_2020_reclassified.csv')[
			cell_type != 'ambiguous' & id != 'sc5rJUQ064_CCATGTCCATCCCATC',
			-c('cell_type_author', 'cell_type_lenient')
		]
	),
	coadread = quote(
		fread('../data_and_figures/lee_crc_2020_smc_reclassified.csv')[cell_type != 'ambiguous', -c('cell_type_author', 'cell_type_lenient')]
	),
    hnsc = quote(fread('../data_and_figures/puram_hnscc_2017_reclassified.csv')[cell_type != 'ambiguous', -'cell_type_author']),
	lihc = quote(
		fread('../data_and_figures/ma_liver_2019_reclassified.csv')[cell_type != 'ambiguous', -c('cell_type_author', 'cell_type_lenient')]
	),
	luad = quote(fread('../data_and_figures/kim_luad_2020_reclassified.csv')[cell_type != 'ambiguous', -'cell_type_author']),
    lusc = quote(
		fread('../data_and_figures/qian_lung_2020_reclassified.csv')[
			disease == 'LUSC' & cell_type != 'ambiguous',
			-c('disease', 'cell_type_author', 'cell_type_lenient')
		]
	),
	ov = quote(
		fread('../data_and_figures/qian_ovarian_2020_reclassified.csv')[
			cell_type != 'ambiguous' & !(id %in% c('scrSOL001_TCATTTGTCTGTCAAG', 'scrSOL004_TTGCCGTTCTCCTATA')),
			-c('cell_type_author', 'cell_type_lenient')
		]
	),
    paad = quote(
		fread('../data_and_figures/peng_pdac_2019_reclassified.csv')[
			cell_type != 'ambiguous' & !(id %in% c('T8_TGGTTCCTCGCATGGC', 'T17_CGTGTAACAGTACACT')),
			-'cell_type_author'
		]
	)
)

tcga_sc_map <- c(
	brca_luminal_a = 'brca',
	brca_luminal_b = 'brca',
	brca_basal_like = 'brca',
	brca_her2_enriched = 'brca',
	coad = 'coadread',
	hnsc_atypical = 'hnsc',
	hnsc_classical = 'hnsc',
	hnsc_mesenchymal_basal = 'hnsc',
	lihc = 'lihc',
	luad_proximal_inflammatory = 'luad',
	luad_proximal_proliferative = 'luad',
	luad_terminal_respiratory_unit = 'luad',
	lusc_basal = 'lusc',
	lusc_classical = 'lusc',
	lusc_secretory = 'lusc',
	ov_differentiated = 'ov',
	ov_immunoreactive = 'ov',
    ov_mesenchymal = 'ov',
	ov_proliferative = 'ov',
	paad_basal_moffitt = 'paad',
    paad_classical_moffitt = 'paad',
	read = 'coadread'
)

sc_comp_diff <- sapply(
	names(tcga_sc_map),
	function(ct) {

		cat(paste0(ct, '\n'))

		sc_data <- eval(sc_metadata[[tcga_sc_map[ct]]])
        sc_data <- sc_data[, c('id', 'cell_type', with(deconv_data[[ct]], genes_filtered[genes_filtered %in% names(sc_data)])), with = FALSE]
		sc_data <- melt(sc_data, id.vars = c('id', 'cell_type'), variable.name = 'gene', value.name = 'expression_level')

		# Centre genes:
		gene_averages <- sc_data[
			,
			.(ave_exp = mean(c(mean(expression_level[cell_type == 'cancer']), mean(expression_level[cell_type == 'caf'])))),
			by = .(symbol = gene)
		]
		sc_data[, expression_level := expression_level - gene_averages[symbol == unique(gene), ave_exp], by = gene]

		# Centre cells:
		sc_data[, expression_level := expression_level - mean(expression_level), by = id]

		sc_data[, .(d = mean(expression_level[cell_type == 'cancer']) - mean(expression_level[cell_type == 'caf'])), by = gene][, setNames(d, gene)]

	},
	simplify = FALSE,
	USE.NAMES = TRUE
)

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

annotation_agreement[
    !is.na(extra),
    sc_comp_diff := do.call(
        function(x) {
            rmx <- runmean(x, 30)/max(abs(runmean(x, 30)))
            mean(head(rmx, length(x)/3)) - mean(tail(rmx, length(x)/3))
        },
        list(x = sc_comp_diff[[cancer_type]][with(deconv_data[[cancer_type]], genes_filtered[ordering])])
    ),
    by = cancer_type
]

barplot_annotations_data <- copy(annotation_agreement)[
    ,
    c('pur', 'ccle', 'extra') := lapply(
        .SD,
        function(x) {
            x[!is.na(x)] <- x[!is.na(x)]/mean(x[!is.na(x)])
            x
        }
    ),
    .SDcols = c('pur', 'ccle', 'sc_comp_diff') # 'extra')
]
barplot_annotations_data[, sc_comp_diff := NULL] # We essentially replaced extra with sc_comp_diff, so can delete the latter
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
    'ov_mesenchymal', 'ov_proliferative', 'paad_basal_moffitt', 'paad_classical_moffitt', 'prad', 'read', 'stad_cin', 'stad_gs', 'stad_msi', 'ucec')
nice_names_for_figure <- c('BLCA - Luminal-Infiltrated', 'BLCA - Luminal-Papillary', 'BLCA - Basal-Squamous', 'BRCA - Luminal A', 'BRCA - Luminal B',
    'BRCA - Basal-like', 'BRCA - HER2-enriched', 'COAD', 'ESCA - Adenocarcinoma', 'HNSC - Malignant-Basal', 'HNSC - Classical', 'KIRP', 'LIHC',
    'LUAD - Squamoid', 'LUAD - Magnoid', 'LUAD - Bronchioid', 'LUSC - Classical', 'LUSC - Secretory', 'OV - Differentiated', 'OV - Immunoreactive',
    'OV - Mesenchymal', 'OV - Proliferative', 'PAAD - Basal', 'PAAD - Classical', 'PRAD', 'READ', 'STAD - CIN', 'STAD - GS', 'STAD - MSI', 'UCEC')

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
    'luad_terminal_respiratory_unit', 'lusc_classical', 'lusc_secretory', 'ov_differentiated', 'ov_immunoreactive', 'ov_proliferative',
    'paad_basal_moffitt', 'paad_classical_moffitt', 'read', 'stad_cin', 'stad_msi', 'ucec')
nice_names_for_figure <- c('BLCA - Luminal-Papillary', 'BLCA - Basal-Squamous', 'BRCA - Luminal A', 'BRCA - Luminal B', 'BRCA - Basal-like',
    'BRCA - HER2-enriched', 'COAD', 'ESCA - Adenocarcinoma', 'HNSC - Malignant-Basal', 'HNSC - Classical', 'LUAD - Squamoid', 'LUAD - Magnoid',
    'LUAD - Bronchioid', 'LUSC - Classical', 'LUSC - Secretory', 'OV - Differentiated', 'OV - Immunoreactive', 'OV - Proliferative', 'PAAD - Basal',
    'PAAD - Classical', 'READ', 'STAD - CIN', 'STAD - MSI', 'UCEC')





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

ct_hclust_heatmap <- ggplot(ct_hclust_heatmap_data, aes(x = ct1, y = ct2, fill = corr)) +
    geom_raster() +
    scale_fill_gradientn(
        colours = c('#798234', '#a3ad62', '#d0d3a2', '#fdfbe4', '#f0c6c3', '#df91a3', '#d46780'), # CARTO ArmyRose diverging colour palette
        limits = c(-0.6, 0.6),
        breaks = c(-0.6, 0, 0.6),
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
        colours = c('#798234', '#a3ad62', '#d0d3a2', '#fdfbe4', '#f0c6c3', '#df91a3', '#d46780'), # CARTO ArmyRose diverging colour palette
        limits = c(-0.6, 0.6),
        breaks = c(-0.6, 0, 0.6),
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

ct_gw_cut <- cutree(ct_gw[[1]], 3)
ct_clust <- data.table(cancer_type = names(ct_gw_cut), cluster = ct_gw_cut)[, memb_strength := silhouette(ct_gw_cut, 1 - ct_cor)[, 'sil_width']]

# Make the names of the cluster match "1", "2", "3" in a consistent order:
setkey(ct_clust, cancer_type)
ct_clust_map <- mapvalues(1:3, 1:3, ct_clust[c('brca_luminal_a', 'hnsc_mesenchymal_basal', 'stad_cin'), cluster])
ct_clust[, cluster_manual := mapvalues(cluster, 1:3, ct_clust_map)]

# Define cluster labels (including intermediates):
ct_clust[
    ,
    cluster_final := switch(
        (memb_strength > 0.2) + 1,
        'Intermediate',
        mapvalues(
            cluster_manual,
            1:3,
            c('Cluster 1 - Gynaecological', 'Cluster 2 - Squamous-like', 'Cluster 3 - Gastro-intestinal'),
            warn_missing = FALSE
        )
    ),
    by = cancer_type
]

clust_memb_data <- ct_clust[
    ,
    .(
        cancer_type = factor(cancer_type, levels = ct_gw[[1]]$labels[ct_gw[[1]]$order]),
        memb_strength = memb_strength,
        cluster = factor(
            cluster_final,
            levels = c('Cluster 1 - Gynaecological', 'Cluster 2 - Squamous-like', 'Cluster 3 - Gastro-intestinal', 'Intermediate')
        )
    )
]

clust_memb <- ggplot(clust_memb_data, aes(x = cancer_type, y = 0, fill = cluster)) +
    geom_raster() +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_fill_manual(values = c(brewer.pal(3, 'Dark2'), 'lightgrey')) +
    theme(
        axis.ticks = element_blank(),
        axis.ticks.length = unit(0, 'pt'),
        axis.text = element_blank(),
        axis.title = element_blank(),
        legend.key = element_rect(size = 1, colour = 'white'),
        legend.key.size = unit(15, 'pt')
    ) +
    guides(fill = guide_legend(label.position = 'left')) +
    labs(fill = NULL)

silvals <- ggplot(clust_memb_data, aes(x = cancer_type, y = 0, fill = memb_strength)) +
    geom_raster() +
    scale_y_continuous(expand = c(0, 0)) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_fill_gradientn(
        colours = c(
            sequential_hcl(25, h = 190, c = 70, l = c(90, 99), power = 0.7),
            sequential_hcl(25, h = 60, c = 100, l = c(90, 99), power = 0.7, rev = TRUE)
        ),
        limits = c(0, 0.4),
        breaks = c(0, 0.2, 0.4),
        labels = c('0' = '0', '0.2' = '0.2', '0.4' = '0.4'),
        oob = squish
    ) +
    theme(axis.ticks = element_blank(), axis.ticks.length = unit(0, 'pt'), axis.text = element_blank(), axis.title = element_blank()) +
    guides(
        fill = guide_colourbar(
            direction = 'horizontal',
            title = 'Silhouette',
            title.position = 'top',
            barwidth = unit(50, 'pt'),
            barheight = unit(8, 'pt')
        )
    )

score_diff_table <- data.table(symbol = names(head(sort(apply(rank_mat, 1, quantile, 0.25)), 100)))[
    ,
    (paste0('score_diff_', 1:3)) := lapply(
        1:3,
        function(i) {
            cts <- list(
                in_i = ct_clust[cluster_manual == i & memb_strength > 0.2, cancer_type],
                not_i = ct_clust[cluster_manual != i & memb_strength > 0.2, cancer_type]
            )
            scores_data_transformed[symbol, rowMeans(.SD[, cts$in_i, with = FALSE]) - rowMeans(.SD[, cts$not_i, with = FALSE])]
        }
    )
] %>% setnames('symbol', 'gene')
score_diff_table[, c('which_max', 'which_max_score') := .(which.max(as.numeric(.SD)), max(as.numeric(.SD))), by = gene]
ct_clust_distinct_genes <- score_diff_table[order(which_max, -which_max_score), .(gene = gene[1:20]), by = which_max]

scores_heatmap_data <- scores_data_transformed[
    unique(unlist(ct_clust_distinct_genes$gene)),
    c('gene', ct_clust[memb_strength > 0.2, cancer_type]),
    with = FALSE
]
names(scores_heatmap_data) <- mapvalues(names(scores_heatmap_data), ct_to_keep, nice_names_for_figure, warn_missing = FALSE)

scores_heatmap <- deconv_heatmap(
    scores_heatmap_data,
    order_genes_fun = 'seriate',
    order_genes_method = 'OLO_ward',
    order_analyses_fun = 'seriate',
    order_analyses_method = 'OLO_ward',
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
            geom_vline(xintercept = c(8.5, 14.5), size = 0.8) +
            theme(
                panel.border = element_rect(fill = NA, size = 0.8),
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
    plot_title = 'Common EMT and stroma genes'
)

pdf('../data_and_figures/final_figures_resubmission/S11.pdf', width = 7.5, height = 14)
htmp_emt_caf$heatmap + theme(
    axis.text.x = element_text(
        colour = mapvalues(
            ct_clust[
                with(htmp_emt_caf, mapvalues(analyses[ordering_analyses], nice_names_for_figure, ct_to_keep)),
                .(new_clust = switch((memb_strength > 0.2) + 1, 0L, cluster_manual)),
                by = cancer_type
            ]$new_clust,
            0:3,
            c('darkgrey', brewer.pal(3, 'Dark2'))
        )
    )
)
dev.off()





# Volcano plot for figure 4: the idea is that genes which occur in more lists will have a larger sample size and therefore greater chance of becoming
# significant.  The fold change is replaced by average EMT-CAF score.

scores_volcano_plot_data <- scores_data_transformed[, .(ave_score = rowMeans(.SD), signif_val = t.test(as.numeric(.SD))$p.value), by = gene]

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
# geom_label_repel().

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
        data = scores_volcano_plot_data[gene %in% c('ITGA11', 'OLFML2B', 'LUM', 'ACTA2', 'TAGLN', 'POSTN')],
        nudge_x = 0.18,
        size = 3,
        point.padding = 0.2,
        segment.size = 0.3
    ) +
    geom_text_repel(
        aes(label = gene),
        data = scores_volcano_plot_data[gene == 'MMP2'],
        nudge_x = 0.18,
        nudge_y = 0.3,
        size = 3,
        point.padding = 0.2,
        segment.size = 0.3
    ) +
    geom_text_repel(
        aes(label = gene),
        data = scores_volcano_plot_data[gene %in% c('ASPN', 'SPARC', 'COL1A2', 'COL3A1', 'FAP', 'COL1A1', 'THY1', 'VCAN', 'ECM2')],
        nudge_x = -0.18,
        size = 3,
        point.padding = 0.2,
        segment.size = 0.3
    ) +
    geom_text_repel(
        aes(label = gene),
        data = scores_volcano_plot_data[gene %in% c('LUZP1', 'LAMC1', 'LAMC2', 'LAMA3', 'ITGA2', 'CD44')],
        nudge_x = 0.15,
        size = 3,
        point.padding = 0.2,
        segment.size = 0.3
    ) +
    geom_text_repel(
        aes(label = gene),
        data = scores_volcano_plot_data[gene == 'MCM7'],
        nudge_x = -0.1,
        size = 3,
        point.padding = 0.2,
        segment.size = 0.3
    ) +
    geom_text_repel(
        aes(label = gene),
        data = scores_volcano_plot_data[gene %in% c('VCL', 'MSN', 'PVR', 'ITGB1')],
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





top_dendro <- dendro(ct_hclust, 'bottom') + theme(plot.margin = unit(c(20, 0, 0, 0), 'pt'))
right_dendro <- dendro(ct_hclust, 'left') + theme(plot.margin = unit(c(0, 20, 0, 0), 'pt'))
corr_heatmap <- ct_gw_heatmap +
    geom_rect(xmin = 1.5, xmax = 9.5, ymin = 1.5, ymax = 9.5, size = 0.8, colour = brewer.pal(3, 'Dark2')[1], fill = NA) + # Gynaecological
    geom_rect(xmin = 17.5, xmax = 23.5, ymin = 17.5, ymax = 23.5, size = 0.8, colour = brewer.pal(3, 'Dark2')[2], fill = NA) + # Squamous-like
    geom_rect(xmin = 11.5, xmax = 16.5, ymin = 11.5, ymax = 16.5, size = 0.8, colour = brewer.pal(3, 'Dark2')[3], fill = NA) + # Gastro-Intestinal
    theme(
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title = element_text(),
        legend.position = 'none',
        plot.margin = unit(c(0, 0, 0, 5.5), 'pt')
    ) +
    labs(x = 'Cancer types', y = 'Cancer types')

aligned_volcano_heatmap <- align_plots(
    scores_volcano_plot,
    corr_heatmap + theme(legend.position = 'none'),
    align = 'v',
    axis = 'l'
)
aligned_top_dendro_heatmap <- align_plots(
    top_dendro,
    silvals + theme(legend.position = 'none', plot.margin = unit(c(0, 0, 3, 0), 'pt')),
    clust_memb + theme(legend.position = 'none', plot.margin = unit(c(0, 0, 3, 0), 'pt')),
    aligned_volcano_heatmap[[2]],
    align = 'v'
)
aligned_right_dendro_heatmap <- align_plots(
    aligned_top_dendro_heatmap[[4]],
    clust_memb + coord_flip() + theme(legend.position = 'none', plot.margin = unit(c(0, 0, 0, 3), 'pt')),
    silvals + coord_flip() + theme(legend.position = 'none', plot.margin = unit(c(0, 0, 0, 3), 'pt')),
    right_dendro,
    align = 'h'
)

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
                blank_plot(),
                get_legend(silvals + theme(legend.position = c(0.4, 0.2))),
                aligned_top_dendro_heatmap[[2]],
                blank_plot(),
                blank_plot(),
                blank_plot(),
                aligned_top_dendro_heatmap[[3]],
                blank_plot(),
                blank_plot(),
                blank_plot(),
                aligned_right_dendro_heatmap[[1]],
                aligned_right_dendro_heatmap[[2]],
                aligned_right_dendro_heatmap[[3]],
                aligned_right_dendro_heatmap[[4]],
                nrow = 4,
                ncol = 4,
                rel_heights = c(1, 0.2, 0.2, 5),
                rel_widths = c(5.1, 0.2, 0.2, 0.9)
            ),
            nrow = 2,
            ncol = 1,
            rel_heights = c(9.5, 8.5)
        ),
        plot_grid(
            get_legend(
                ct_hclust_heatmap +
                    theme(legend.justification = c(0.3, 0.8)) +
                    guides(
                        fill = guide_colourbar(
                            direction = 'horizontal',
                            title = 'Correlation',
                            title.position = 'top',
                            barwidth = unit(70, 'pt'),
                            barheight = unit(8, 'pt')
                        )
                    )
            ),
            get_legend(clust_memb + theme(legend.justification = c(0.5, 1))),
            nrow = 1,
            ncol = 2
        ),
        nrow = 2,
        ncol = 1,
        rel_heights = c(18, 2)
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





# Control clustering to see whether our 3 clusters arise naturally from the expression data and don't really reflect EMT:

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

pdf('../data_and_figures/final_figures_resubmission/S12.pdf', width = 10, height = 10)
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

set.seed(2423)

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
            test_x_expr = list(
                quote(variable > 1),
                quote(variable != 'complete remission/response'),
                quote(variable < quantile(variable, 0.4)),
                quote(variable %in% c('yes', 'present')),
                quote(startsWith(variable, 't3') | startsWith(variable, 't4')),
				quote(startsWith(variable, 'n2') | startsWith(variable, 'n3')),
                quote(startsWith(variable, 'm1')),
                quote(variable %in% c('gb', 'g1', 'g2', 'low grade'))
            ),
            test_y_expr = list(
                quote(variable == 0),
                quote(variable == 'complete remission/response'),
                quote(variable > quantile(variable, 0.6)),
                quote(variable %in% c('no', 'absent')),
                quote(startsWith(variable, 't0') | startsWith(variable, 't1') | startsWith(variable, 't2')),
                quote(startsWith(variable, 'n0') | startsWith(variable, 'n1')),
                quote(startsWith(variable, 'm0') | startsWith(variable, 'cm0')),
                quote(variable %in% c('g3', 'g4', 'high grade'))
            ),
            min_samples = 10,
            score_method = 'signature_score',
            n = 100,
            nbin = 20
        )
    },
    simplify = FALSE,
    USE.NAMES = TRUE
)

clin_cor_sample_ids <- sapply(clin_cor, `[[`, 'sample_ids', simplify = FALSE, USE.NAMES = TRUE)
clin_cor_patient_ids <- sapply(clin_cor, `[[`, 'patient_ids', simplify = FALSE, USE.NAMES = TRUE)

clin_cor <- sapply(
    clin_cor,
    function(x) {
        clin_cor_data <- x$data[
            ,
            c('nice_test_name', 'nice_variable_name', 'sigval', 'sigval_adj') := .(
                mapvalues(test_name, ct_to_keep, nice_names_for_figure, warn_missing = FALSE),
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
                -log10(pval)*sign(diff),
                -log10(pval_adj)*sign(diff)
            )
        ]
        return(clin_cor_data)
    },
    simplify = FALSE,
    USE.NAMES = TRUE
)





# Heatmaps of clinical test significance values:

# latex2exp::TeX('-log_{10}(p) $\\times$ sign($\\Delta$(score))')
sig_lab <- expression(
    atop(
        `\textbf{Significance of association}` = paste("", bold(paste("Significance of association"))),
        `-log_{10}(p) $\times$ sign($\\Delta$(score))` = paste("-log", phantom()[{paste("10")}], "(", "", "p", ")", " ", "", phantom() %*% phantom(),
        " sign", "(", "", "", Delta, , , , "", "(", "score", ")", "", ")", "")
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
                '#F1B6DA', '#DE77AE', '#C51B7D', '#8E0152'), # PiYG palette with fatter waist
            limits = c(-3, 3),
            breaks = c(-3, -2, -1, 0, 1, 2, 3),
            labels = c('-3' = '\u2264 -3', '-2' = '-2', '-1' = '-1', '0' = '0', '1' = '1', '2' = '2', '3' = '\u2265 3'),
            x_lab = NULL,
            y_lab = NULL,
            legend_title = sig_lab,
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





# Scatterplots of significance of association for EMT vs. that for CAFs:

# latex2exp::TeX('-log_{10}(p_{CAF}) $\\times$ sign($\\Delta$(score))')
sig_lab_caf <- expression(
    atop(
        `\textbf{Significance of association for CAFs}` = paste("", bold(paste("Significance of association for CAFs"))),
        `-log_{10}(p_{CAF}) $\times$ sign($\\Delta$(score))` = paste("-log", phantom()[{paste("10")}], "(", "", "p", phantom()[{paste("CAF")}],
            ")", "", " ", "", phantom() %*% phantom(), " sign", "(", "", "", Delta, , , , "", "(", "score", ")", "", ")", "")
    )
)

# latex2exp::TeX('-log_{10}(p_{pEMT}) $\\times$ sign($\\Delta$(score))')
sig_lab_emt <- expression(
    atop(
        `\textbf{Significance of association for pEMT}` = paste("", bold(paste("Significance of association for pEMT"))),
        `-log_{10}(p_{pEMT}) $\times$ sign($\\Delta$(score))` = paste("-log", phantom()[{paste("10")}], "(", "", "p", phantom()[{paste("pEMT")}],
            ")", "", " ", "", phantom() %*% phantom(), " sign", "(", "", "", Delta, , , , "", "(", "score", ")", "", ")", "")
    )
)

emt_caf_sig_data <- merge(clin_cor$deconv_emt, clin_cor$deconv_caf, by = c('nice_test_name', 'nice_variable_name'))[
    ,
    .(
        test_name = nice_test_name,
        variable_name = nice_variable_name,
        sig_emt = sigval.x,
        sig_caf = sigval.y,
        sig_adj_emt = sigval_adj.x,
        sig_adj_caf = sigval_adj.y
    )
]

pval_adj_threshold <- adjust_threshold_bh(
    c(clin_cor$deconv_emt[variable_name != 'pathologic_m', pval], clin_cor$deconv_caf[variable_name != 'pathologic_m', pval])
)

# Volcano plots:

clin_vol_data <- merge(clin_cor$deconv_emt, clin_cor$deconv_caf, by = c('nice_test_name', 'nice_variable_name'))[
    ,
    .(
        test_name = nice_test_name,
        variable_name = nice_variable_name,
        sig_emt = abs(sigval.x),
        sig_caf = abs(sigval.y),
        diff_emt = diff.x,
        diff_caf = diff.y
    )
]

# Melt the data so that we have caf and emt in one column, which we can colour/group by:
clin_vol_data <- merge(
    melt(
        clin_vol_data[, .(test_name, variable_name, sig_emt, sig_caf)],
        id.vars = c('test_name', 'variable_name'),
        variable.name = 'emt_type',
        value.name = 'sigval'
    )[, emt_type := gsub('sig_', '', emt_type)],
    melt(
        clin_vol_data[, .(test_name, variable_name, diff_emt, diff_caf)],
        id.vars = c('test_name', 'variable_name'),
        variable.name = 'emt_type',
        value.name = 'diff'
    )[, emt_type := gsub('diff_', '', emt_type)],
    by = c('test_name', 'variable_name', 'emt_type')
)

clin_vol_data_adjvar <- rbindlist(
    lapply(
        c('Grade', 'Lymph node metastasis', 'Lymphovascular invasion', 'M stage', 'N stage', 'Reduced survival', 'T stage', 'Therapy resistance'),
        function(clin_feat) {
            pthresh <- adjust_threshold_bh(
                c(clin_cor$deconv_emt[nice_variable_name == clin_feat, pval], clin_cor$deconv_caf[nice_variable_name == clin_feat, pval])
            )
            dt <- clin_vol_data[variable_name == clin_feat]
            if(is.na(pthresh)) {
                dt[, sig := 'not_significant']
            } else {
                dt[, sig := switch((sigval > -log10(pthresh)) + 1, 'not_significant', 'significant'), by = .(test_name, emt_type)]
            }
            return(dt)
        }
    )
)

clin_vol <- ggplot(
    clin_vol_data_adjvar[variable_name %in% c('Grade', 'Lymph node metastasis', 'N stage', 'Therapy resistance')],
    aes(x = diff, y = sigval)
) +
    facet_wrap(facets = vars(variable_name), nrow = 2, ncol = 2) +
    geom_vline(xintercept = 0, colour = 'grey', size = 0.25, linetype = 'dashed') +
    geom_text_repel(
        aes(label = test_name),
        data = clin_vol_data_adjvar[variable_name == 'Grade' & emt_type == 'caf' & test_name == 'ESCA - Adenocarcinoma'],
        point.padding = 0.1, size = 3.5, lineheight = 0.75, segment.size = 0.3, nudge_x = 0.7, nudge_y = 0.2
    ) +
    geom_text_repel(
        aes(label = test_name),
        data = clin_vol_data_adjvar[variable_name == 'Grade' & emt_type == 'caf' & test_name == 'STAD - CIN'],
        point.padding = 0.1, size = 3.5, lineheight = 0.75, segment.size = 0.3, nudge_x = 0.4, nudge_y = -0.1
    ) +
    geom_text_repel(
        aes(label = test_name),
        data = clin_vol_data_adjvar[variable_name == 'Grade' & emt_type == 'emt' & test_name == 'UCEC'],
        point.padding = 0.1, size = 3.5, lineheight = 0.75, segment.size = 0.3, nudge_x = -0.05, nudge_y = 0.05
    ) +
    geom_text_repel(
        aes(label = test_name),
        data = clin_vol_data_adjvar[variable_name == 'Lymph node metastasis' & emt_type == 'caf' & test_name == 'HNSC - Malignant-Basal'],
        point.padding = 0.1, size = 3.5, lineheight = 0.75, segment.size = 0.3, nudge_x = -0.7
    ) +
    geom_text_repel(
        aes(label = test_name),
        data = clin_vol_data_adjvar[variable_name == 'Lymph node metastasis' & emt_type == 'caf' & test_name == 'COAD'],
        point.padding = 0.1, size = 3.5, lineheight = 0.75, segment.size = 0.3, nudge_x = 0.3, nudge_y = -0.05
    ) +
    geom_text_repel(
        aes(label = test_name),
        data = clin_vol_data_adjvar[variable_name == 'Lymph node metastasis' & emt_type == 'emt' & sig == 'significant'],
        point.padding = 0.1, size = 3.5, lineheight = 0.75, segment.size = 0.3
    ) +
    geom_text_repel(
        aes(label = test_name),
        data = clin_vol_data_adjvar[variable_name == 'N stage' & emt_type == 'emt' & test_name == 'HNSC - Classical'],
        point.padding = 0.1, size = 3.5, lineheight = 0.75, segment.size = 0.3, nudge_x = -0.1, nudge_y = 0.1
    ) +
    geom_text_repel(
        aes(label = test_name),
        data = clin_vol_data_adjvar[variable_name == 'N stage' & test_name %in% c('HNSC - Malignant-Basal', 'READ') & sig == 'significant'],
        point.padding = 0.1, size = 3.5, lineheight = 0.75, segment.size = 0.3
    ) +
    geom_text_repel(
        aes(label = test_name),
        data = clin_vol_data_adjvar[variable_name == 'Therapy resistance' & sig == 'significant'],
        point.padding = 0.1, size = 3.5, lineheight = 0.75, segment.size = 0.3
    ) +
    # geom_text_repel(
    #     # aes(label = str_replace(test_name, ' - ', '\n')),
    #     aes(label = test_name),
    #     data = clin_vol_data_adjvar[variable_name %in% c('Grade', 'Lymph node metastasis', 'N stage', 'Therapy resistance') & sig == 'significant'],
    #     point.padding = 0.1, size = 3.5, lineheight = 0.75, segment.size = 0.3
    # ) +
    geom_point(aes(colour = emt_type), shape = 17, size = 2.5) +
    scale_colour_manual(values = c('caf' = '#2F8AC4', 'emt' = '#ED645A'), labels = c('caf' = 'CAF', 'emt' = 'pEMT')) +
    labs(x = latex2exp::TeX('$\\Delta$(score)'), y = latex2exp::TeX('-log_{10}(p value)'), colour = 'Signature') +
    theme_test() +
    theme(strip.text = element_text(size = 13), legend.text = element_text(size = 11), legend.justification = c(0.1, 0.5))





# Combined figure of clinical correlation heatmaps and 4 volcano plots:

cairo_pdf('../data_and_figures/final_figures_resubmission/5.pdf', width = 8.5, height = 11.6)
plot_grid(
    plot_grid(
        blank_plot(),
        clin_cor_heatmaps$deconv_emt + theme(axis.text.x = element_blank(), legend.position = 'none'),
        clin_cor_heatmaps$deconv_caf +
            theme(
                axis.text.x = element_text(
                    colour = mapvalues(
                        ct_clust[
                            with(ct_hclust, labels[order][labels[order] %in% clin_cor$deconv_caf$test_name]),
                            .(new_clust = switch((memb_strength > 0.2) + 1, 0L, cluster_manual)),
                            by = cancer_type
                        ]$new_clust,
                        c(0, 1, 2, 3),
                        c('darkgrey', brewer.pal(3, 'Dark2'))
                    )
                ),
                legend.position = 'none',
                plot.margin = unit(c(5.5, 5.5, 1, 5.5), 'pt')
            ),
        blank_plot(),
        clin_vol + theme(legend.position = 'none'),
        nrow = 5,
        ncol = 1,
        rel_heights = c(0.4, 1.87, 3.03, 0.7, 5.6)
    ),
    plot_grid(
        get_legend(
            clin_cor_heatmaps$deconv_emt + theme(legend.justification = c(0, 0.595)) + guides(fill = guide_colourbar(title.position = 'right'))
        ),
        get_legend(clin_vol),
        nrow = 2,
        ncol = 1,
        rel_heights = c(5.9, 6.6)
    ),
    nrow = 1,
    ncol = 2,
    rel_widths = c(6, 2.5)
) +
    draw_label('A', x = 0, y = 0.98, hjust = -0.1, vjust = 1.1, size = 20, fontface = 2) +
    draw_label('B', x = 0, y = 0.51, hjust = -0.1, vjust = 1.1, size = 20, fontface = 2)
dev.off()





# To investigate/check the strong negative association of pEMT with LN mets in HNSC Classical:

lnmets_data <- rbindlist(
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
)[
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

# Tumour site:
lnmets_data[, site := clinical_data[individual_id, anatomic_neoplasm_subdivision]]

# Look at DE genes:

scaled_mat <- expression_data[
    meta_data[patient_id %in% lnmets_data[cancer_type == 'HNSC - Classical', individual_id], id],
    set_rownames(scale(as.matrix(.SD)), id),
    .SDcols = -'id'
]

degenes <- sapply(
    c('Zero', 'Multiple'),
    function(x) sort(
        apply(scaled_mat[meta_data[patient_id %in% lnmets_data[cancer_type == 'HNSC - Classical' & lnmets == x, individual_id], id], ], 2, mean),
        decreasing = TRUE
    ),
    simplify = FALSE,
    USE.NAMES = TRUE
)

# Get MSigDB C2, C5 and H gene sets (C2 = Curated gene sets; C5 = GO gene sets):
msigdb_table <- rbindlist(lapply(c('C2', 'C5', 'H'), function(categ) as.data.table(msigdbr(category = categ))))

# We don't get anything significant if we set pvalueCutoff to 0.01.
degenes_gsea <- sapply(
    c('Zero', 'Multiple'),
    function(lnmets_status) {
        sapply(
        unique(msigdb_table$gs_cat),
            function(categ) cbind(
                gs_cat = categ,
                as.data.table(GSEA(degenes[[lnmets_status]], TERM2GENE = msigdb_table[gs_cat == categ, .(gs_name, human_gene_symbol)]))
            ),
            simplify = FALSE,
            USE.NAMES = TRUE
        )
    },
    simplify = FALSE,
    USE.NAMES = TRUE
)

gsea_table <- lapply(
    c('Zero', 'Multiple'),
    function(lnmets_status) rbindlist(
        lapply(degenes_gsea[[lnmets_status]], function(dt) dt[!is.na(ID)][, lnmets := lnmets_status][order(-NES)]),
        fill = TRUE
    )
) %>% rbindlist

# OXPHOS signature expression:
lnmets_data[
    ,
    oxphos_expression := expression_data[
        meta_data[patient_id %in% individual_id, id],
        rowMeans(.SD),
        .SDcols = gsea_table[lnmets == 'Multiple' & ID == 'HALLMARK_OXIDATIVE_PHOSPHORYLATION', strsplit(core_enrichment, '/')[[1]]]
    ],
    by = .(cancer_type, lnmets)
][, pval_oxphos := wilcox.test(oxphos_expression[lnmets == 'Zero'], oxphos_expression[lnmets == 'Multiple'])$p.value, by = cancer_type]

# Cluster the classical tumours:

cell_type_markers <- fread('../../cell_type_markers.csv')[source %in% c('Tirosh', 'Tirosh immune'), .(cell_type, gene)] %>% unique

classical_data_all <- expression_data[deconv_data$hnsc_classical$sample_ids]
gene_vars <- classical_data_all[, sapply(.SD, var), .SDcols = -'id']
classical_data <- classical_data_all[
    ,
    c('id', names(head(sort(gene_vars[!is.na(gene_vars) & gene_vars != 0], decreasing = TRUE), 2000))),
    .SDcols = -'id',
    with = FALSE
]
classical_data[, names(classical_data[, -'id']) := lapply(.SD, scale), .SDcols = -'id']
classical_cormat <- classical_data[, cor(set_colnames(t(.SD), id)), .SDcols = -'id']

set.seed(6707)
classical_spin <- seriate(as.dist(1 - classical_cormat), method = 'SPIN_STS')[[1]]

classical_heatmap <- heat_map(
    classical_cormat,
    get_order(classical_spin),
    colour_limits = c(-0.4, 0.4),
    axis_title_x = 'Tumours',
    axis_title_y = 'Tumours',
    axis_text_x = '',
    axis_text_y = '',
    legend_title = 'Correlation',
    plot_margin = c(1, 5.5, 5.5, 5.5)
)

setkey(lnmets_data, individual_id)
set.seed(3056)
classical_bar_data <- data.table(sample_id = factor(rownames(classical_cormat), levels = rownames(classical_cormat)[get_order(classical_spin)]))[
    ,
    patient_id := meta_data[as.character(sample_id), patient_id]
][
    ,
    c('lnmets', 'site', 'pemt', 'oxphos') := c(
        lnmets_data[patient_id, .(lnmets, mapvalues(site, c('buccal mucosa', 'floor of mouth', 'oral tongue'), rep('oral cavity', 3)))],
        expression_data[
            as.character(sample_id),
            .(
                signature_score(
                    set_colnames(t(as.matrix(.SD)), id),
                    with(deconv_data$hnsc_classical, head(genes_filtered[ordering], 20)),
                    nbin = 20,
                    n = 100
                ),
                signature_score(
                    set_colnames(t(as.matrix(.SD)), id),
                    gsea_table[lnmets == 'Multiple' & ID == 'HALLMARK_OXIDATIVE_PHOSPHORYLATION', strsplit(core_enrichment, '/')[[1]]],
                    nbin = 20,
                    n = 100
                )
            ),
            .SDcols = -'id'
        ]
    )
][, c('pemt', 'oxphos') := lapply(.SD, function(x) {x/(2*sd(x))}), .SDcols = c('pemt', 'oxphos')]

# We could also try heatmaps of pEMT/myogenesis/oxphos genes in individual tumours, instead of averaging.

lnmets_bar <- ggplot(classical_bar_data, aes(x = sample_id, y = 0, fill = lnmets)) +
    geom_raster() +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_fill_manual(values = setNames(brewer.pal(8, 'Set2')[1:2], c('Zero', 'Multiple')), na.value = 'lightgrey') +
    theme(
        axis.ticks = element_blank(),
        axis.ticks.length = unit(0, 'pt'),
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.margin = unit(c(2, 5.5, 2, 5.5), 'pt'),
        legend.justification = 'left'
    ) +
    labs(fill = 'LN metastasis: ')

site_bar <- ggplot(classical_bar_data, aes(x = sample_id, y = 0, fill = mapvalues(site, c('larynx', 'oral cavity'), c('Larynx', 'Oral cavity')))) +
    geom_raster() +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    scale_fill_manual(values = setNames(brewer.pal(8, 'Set2')[4:5], c('Oral cavity', 'Larynx')), na.value = 'lightgrey') +
    theme(
        axis.ticks = element_blank(),
        axis.ticks.length = unit(0, 'pt'),
        axis.text = element_blank(),
        axis.title = element_blank(),
        plot.margin = unit(c(2, 5.5, 2, 5.5), 'pt'),
        legend.justification = 'left'
    ) +
    labs(fill = 'Tumour site: ')

pemt_oxphos_bars <- sapply(
    c('pemt', 'oxphos'),
    function(barvar) {
        ggplot(classical_bar_data, aes(x = sample_id, y = 0, fill = get(barvar))) +
            geom_raster() +
            scale_x_discrete(expand = c(0, 0)) +
            scale_y_continuous(expand = c(0, 0)) +
            scale_fill_gradientn(colours = colorRampPalette(brewer.pal(9, 'YlOrRd'))(50), limits = c(-1, 1), breaks = c(-1, 0, 1), oob = squish) +
            theme(
                axis.ticks = element_blank(),
                axis.ticks.length = unit(0, 'pt'),
                axis.text = element_blank(),
                axis.title.x = element_blank(),
                axis.title.y = element_text(angle = 0, vjust = 0.5), # I want to put hjust = 1 here, but it doesn't affect the aligned plot.
                plot.margin = unit(c(2, 5.5, 2, 5.5), 'pt'),
                legend.justification = c(0, 0.3),
                legend.key.height = unit(10, 'pt')
            ) +
            labs(y = switch((barvar == 'pemt') + 1, 'OXPHOS', 'pEMT'), fill = 'Score\n')
    },
    simplify = FALSE,
    USE.NAMES = TRUE
)

pdf('../data_and_figures/final_figures_resubmission/S13.pdf', width = 8, height = 12)
plot_grid(
    blank_plot(),
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
        theme_test(),
    blank_plot(),
    plot_grid(
        plot_grid(
            plotlist = lapply(
                list(lnmets_bar, site_bar, pemt_oxphos_bars$pemt, pemt_oxphos_bars$oxphos, classical_heatmap),
                `+`,
                theme(legend.position = 'none')
            ),
            nrow = 5,
            ncol = 1,
            align = 'v',
            rel_heights = c(rep(1, 4), 10)
        ),
        plot_grid(
            plotlist = c(
                lapply(list(lnmets_bar, site_bar, pemt_oxphos_bars$pemt), function(x) get_legend(x + theme(legend.direction = 'horizontal'))),
                list(get_legend(classical_heatmap + theme(legend.justification = c(0, 0.1))))
            ),
            nrow = 4,
            ncol = 1,
            rel_heights = c(1, 1, 2, 10)
        ),
        nrow = 1,
        ncol = 2,
        rel_widths = c(1, 0.7)
    ),
    nrow = 4,
    ncol = 1,
    rel_heights = c(0.5, 5.4, 0.7, 5.4)
) +
    draw_label('A', x = 0, y = 0.98, hjust = -0.1, vjust = 1.1, size = 20, fontface = 2) +
    draw_label('B', x = 0, y = 0.48, hjust = -0.1, vjust = 1.1, size = 20, fontface = 2)
dev.off()
