# bsub -q tirosh -R rusage[mem=128000] -o sc_deconv_comp_highscore_100.o -e sc_deconv_comp_highscore_100.e Rscript sc_deconv_comp_highscore_100.R

library(data.table) # 1.12.8
library(ggplot2) # This is version 3.3.1 on my laptop's WSL setup but 3.3.0 on WEXAC.
library(cowplot) # 1.0.0
library(magrittr) # 1.5
library(limma) # 3.42.2
library(plyr) # 1.8.6
library(stringr) # 1.4.0
library(caTools) # 1.18.0
library(RColorBrewer) # 1.1.2
library(scales) # 1.1.1
library(latex2exp) # 0.4.0
library(ggrepel) # 0.8.2
library(egg) # 0.4.5
library(colorspace) # 1.4.1

source('general_functions.R')
source('tcga_functions.R')
source('sc_functions.R')





expression_data <- fread('../../TCGA_data/tcga_expression_data.csv', key = 'id')
meta_data <- fread('../../TCGA_data/tcga_meta_data.csv', key = 'id')

expression_data <- expression_data[meta_data[sample_type != 'normal', id]]
meta_data <- meta_data[sample_type != 'normal']

emt_markers <- fread('../../emt_markers.csv')[, gene := alias2SymbolTable(gene)][source != 'GO', sort(unique(gene))]





sc_metadata <- list(
	brca = list(
		tcga_cancer_types = 'BRCA',
		read_quote = quote(
			fread('../data_and_figures/qian_breast_2020_reclassified.csv')[
				cell_type != 'ambiguous' & id != 'sc5rJUQ064_CCATGTCCATCCCATC',
				-c('cell_type_author', 'cell_type_lenient')
			]
		),
		initial_cell_types = c('b_cell', 'endothelial', 'macrophage', 'mast', 't_cell', 'caf', 'cancer'),
		rare_cell_types = 'dendritic'
	),
	coadread = list(
		tcga_cancer_types = c('COAD', 'READ'),
		read_quote = quote(
			fread('../data_and_figures/lee_crc_2020_smc_reclassified.csv')[cell_type != 'ambiguous', -c('cell_type_author', 'cell_type_lenient')]
		),
		initial_cell_types = c('b_cell', 'endothelial', 'macrophage', 't_cell', 'caf', 'cancer'),
		rare_cell_types = c('epithelial', 'mast')
	),
    hnsc = list(
        tcga_cancer_types = 'HNSC',
        read_quote = quote(fread('../data_and_figures/puram_hnscc_2017_reclassified.csv')[cell_type != 'ambiguous', -'cell_type_author']),
		initial_cell_types = c('endothelial', 'macrophage', 'mast', 't_cell', 'caf', 'cancer'),
		rare_cell_types = c('b_cell', 'dendritic', 'myocyte')
    ),
	lihc = list(
        tcga_cancer_types = 'LIHC',
        read_quote = quote(
			fread('../data_and_figures/ma_liver_2019_reclassified.csv')[cell_type != 'ambiguous', -c('cell_type_author', 'cell_type_lenient')]
		),
		initial_cell_types = c('endothelial', 'macrophage', 't_cell', 'caf', 'cancer'),
		rare_cell_types = c('b_cell', 'hpc-like')
    ),
	luad = list(
		tcga_cancer_types = 'LUAD',
		read_quote = quote(fread('../data_and_figures/kim_luad_2020_reclassified.csv')[cell_type != 'ambiguous', -'cell_type_author']),
		initial_cell_types = c('b_cell', 'endothelial', 'macrophage', 'mast', 't_cell', 'caf', 'cancer'),
		rare_cell_types = 'dendritic'
	),
    lusc = list(
        tcga_cancer_types = 'LUSC',
        read_quote = quote(
			fread('../data_and_figures/qian_lung_2020_reclassified.csv')[
				disease == 'LUSC' & cell_type != 'ambiguous',
				-c('disease', 'cell_type_author', 'cell_type_lenient')
			]
		),
		initial_cell_types = c('b_cell', 'endothelial', 'macrophage', 'mast', 't_cell', 'caf', 'cancer'),
		rare_cell_types = c('dendritic', 'erythroblast')
    ),
	ov = list(
		tcga_cancer_types = 'OV',
		read_quote = quote(
			fread('../data_and_figures/qian_ovarian_2020_reclassified.csv')[
				cell_type != 'ambiguous' & !(id %in% c('scrSOL001_TCATTTGTCTGTCAAG', 'scrSOL004_TTGCCGTTCTCCTATA')),
				-c('cell_type_author', 'cell_type_lenient')
			]
		),
		initial_cell_types = c('b_cell', 'endothelial', 'macrophage', 't_cell', 'caf', 'cancer'),
		rare_cell_types = NULL
	),
    paad = list(
        tcga_cancer_types = 'PAAD',
        read_quote = quote(
			fread('../data_and_figures/peng_pdac_2019_reclassified.csv')[
				cell_type != 'ambiguous' & !(id %in% c('T8_TGGTTCCTCGCATGGC', 'T17_CGTGTAACAGTACACT')),
				-'cell_type_author'
			]
		),
		initial_cell_types = c('b_cell', 'endothelial', 'macrophage', 't_cell', 'caf', 'cancer'),
		rare_cell_types = c('acinar', 'ductal_2', 'endocrine')
    )
)





sc_cancer_caf_args <- list(
    brca = list(
        seed = 3718,
        genes_filter_fun = function(x) {log2(mean(10*(2^x - 1)) + 1) >= 4.5 | sum(x >= 7) >= length(x)/100},
        scores_filter_fun = function(x) {log2(mean(10*(2^x - 1)) + 1) >= 6 | sum(x >= 7) >= length(x)/100}
    ),
    coadread = list(
        seed = 3361,
        genes_filter_fun = function(x) {log2(mean(10*(2^x - 1)) + 1) >= 4 | sum(x >= 7) >= length(x)/100},
        scores_filter_fun = function(x) {log2(mean(10*(2^x - 1)) + 1) >= 5.5 | sum(x >= 7) >= length(x)/100}
    ),
    hnsc = list(
        seed = 8511,
        genes_filter_fun = function(x) {log2(mean(10*(2^x - 1)) + 1) >= 4 | sum(x >= 7) >= length(x)/100},
        scores_filter_fun = function(x) {log2(mean(10*(2^x - 1)) + 1) >= 6 | sum(x >= 7) >= length(x)/100}
    ),
    lihc = list(
        seed = 4376,
        genes_filter_fun = function(x) {log2(mean(10*(2^x - 1)) + 1) >= 4 | sum(x >= 7) >= length(x)/100},
        scores_filter_fun = function(x) {log2(mean(10*(2^x - 1)) + 1) >= 6 | sum(x >= 7) >= length(x)/100}
    ),
	luad = list(
        seed = 1096,
        genes_filter_fun = function(x) {log2(mean(10*(2^x - 1)) + 1) >= 4 | sum(x >= 7) >= length(x)/100},
        scores_filter_fun = function(x) {log2(mean(10*(2^x - 1)) + 1) >= 5.5 | sum(x >= 7) >= length(x)/100}
    ),
    lusc = list(
        seed = 2566,
        genes_filter_fun = function(x) {log2(mean(10*(2^x - 1)) + 1) >= 4 | sum(x >= 7) >= length(x)/100},
        scores_filter_fun = function(x) {log2(mean(10*(2^x - 1)) + 1) >= 5.5 | sum(x >= 7) >= length(x)/100}
    ),
    ov = list(
        seed = 456,
        genes_filter_fun = function(x) {log2(mean(10*(2^x - 1)) + 1) >= 4 | sum(x >= 7) >= length(x)/100},
        scores_filter_fun = function(x) {log2(mean(10*(2^x - 1)) + 1) >= 5.5 | sum(x >= 7) >= length(x)/100}
    ),
    paad = list(
        seed = 5368,
        genes_filter_fun = function(x) {log2(mean(10*(2^x - 1)) + 1) >= 4.5 | sum(x >= 7) >= length(x)/100},
        scores_filter_fun = function(x) {log2(mean(10*(2^x - 1)) + 1) >= 6 | sum(x >= 7) >= length(x)/100}
    )
)





sample_sizes <- list(
    brca = 100,
    coadread = 100,
    hnsc = 100,
    lihc = 100,
    luad = 100,
	lusc = 100,
    ov = 100,
    paad = 100
)

sc_cancer_caf <- readRDS('../data_and_figures/sc_cancer_caf.rds')

deconv_data <- readRDS('../data_and_figures/deconv_data.rds')
# deconv_plots <- readRDS('../data_and_figures/deconv_plots.rds')

cell_type_markers <- fread('../../cell_type_markers.csv')
subtypes_data <- fread('../../TCGA_data/tcga_subtypes_data.csv', key = 'id')
ccle <- fread('../../CCLE_cancer_type_Av.csv', key = 'gene_id')

ct_to_keep <- c('brca_luminal_a', 'brca_luminal_b', 'brca_basal_like', 'brca_her2_enriched', 'coad', 'hnsc_mesenchymal_basal', 'hnsc_classical',
	'luad_proximal_inflammatory', 'luad_proximal_proliferative', 'lusc_classical', 'lusc_secretory', 'ov_differentiated', 'ov_immunoreactive',
	'ov_proliferative', 'paad', 'read')

tcga_sc_map <- c(
	brca_luminal_a = 'brca',
	brca_luminal_b = 'brca',
	brca_basal_like = 'brca',
	brca_her2_enriched = 'brca',
	coad = 'coadread',
	hnsc_mesenchymal_basal = 'hnsc',
	hnsc_classical = 'hnsc',
	# hnsc_atypical = 'hnsc',
	# lihc = 'lihc',
	luad_proximal_inflammatory = 'luad',
	luad_proximal_proliferative = 'luad',
	# lusc_basal = 'lusc',
	lusc_classical = 'lusc',
	# lusc_primitive = 'lusc',
	lusc_secretory = 'lusc',
	ov_differentiated = 'ov',
	ov_immunoreactive = 'ov',
	ov_proliferative = 'ov',
	paad = 'paad',
	read = 'coadread'
)

deconv_titles <- list(
    brca_luminal_a = 'BRCA - Luminal A',
    brca_luminal_b = 'BRCA - Luminal B',
    brca_basal_like = 'BRCA - Basal-like',
    brca_her2_enriched = 'BRCA - HER2-enriched',
    coad = 'COAD',
    hnsc_mesenchymal_basal = 'HNSC - Mesenchymal & Basal',
    hnsc_classical = 'HNSC - Classical',
    # hnsc_atypical = 'HNSC - Atypical',
    # lihc = 'LIHC',
    luad_proximal_inflammatory = 'LUAD - Proximal-inflammatory',
    luad_proximal_proliferative = 'LUAD - Proximal-proliferative',
    # lusc_basal = 'LUSC - Basal',
    lusc_classical = 'LUSC - Classical',
    # lusc_primitive = 'LUSC - Primitive',
    lusc_secretory = 'LUSC - Secretory',
    ov_differentiated = 'OV - Differentiated',
    ov_immunoreactive = 'OV - Immunoreactive',
    ov_proliferative = 'OV - Proliferative',
    paad = 'PAAD',
    read = 'READ'
)

# Re-make deconv plots to use centred CCLE and single cell data:
deconv_plots <- sapply(
	ct_to_keep,
    # names(deconv_args_per_ct),
    function(ct) {
        cat(ct, '\b...')
        deconv_ct <- deconv_data[[ct]]
        deconv_plot_ct <- deconvolve_emt_caf_plots(
            data = deconv_ct,
            # Include the following only if you want epithelial scores (takes much longer):
            # expression_data = expression_data,
            heatmap_axis_title = '', # Change heat_map() function so I can put NULL here
            heatmap_legend_title = 'Coexpression',
            heatmap_colours = rev(colorRampPalette(brewer.pal(11, "RdBu"))(50)),
			heatmap_colour_limits = c(-0.6, 0.6),
		    heatmap_legend_breaks = c(-0.6, -0.3, 0, 0.3, 0.6),
            heatmap_annotations = c('SNAI1', 'SNAI2', 'TWIST1', 'ZEB1', 'ZEB2'),
            heatmap_annotations_nudge = 0.3,
            purity_colours = rev(colorRampPalette(brewer.pal(11, "PuOr"))(50)),
            purity_colour_limits = c(-0.3, 0.3),
            purity_legend_breaks = c(-0.3, 0, 0.3),
            purity_legend_title = 'Correlation\nwith purity',
            purity_legend_direction = 'vertical',
            purity_axis_title = NULL,
            ccle_colours = rev(colorRampPalette(brewer.pal(11, 'PiYG')[c(1:3, 6, 9:11)])(50)),
			ccle_fun = function(x) {runmean(x - mean(x), 30)/max(abs(runmean(x - mean(x), 30)))},
            ccle_legend_title = 'Tumours vs.\ncell lines',
            ccle_legend_direction = 'vertical',
            ccle_axis_title = NULL,
            extra_colours = colorRampPalette(c('turquoise4', 'turquoise', 'azure', 'gold2', 'gold4'))(50),
            extra_axis_title = NULL,
            extra_legend_title = 'scRNA-seq',
            extra_legend_direction = 'vertical',
            bar_legend_width = NULL,
            bar_legend_height = NULL,
			plot_title = deconv_titles[ct]
        )
        cat('Done!\n')
        deconv_plot_ct
    },
    simplify = FALSE,
    USE.NAMES = TRUE
)

# Add data for single cell comparison colour bar:
# sc_comp_diff <- sapply(
# 	names(deconv_data),
# 	function(ct) {
#
# 		cat(paste0(ct, '\n'))
#
# 		sc_data <- eval(sc_metadata[[tcga_sc_map[ct]]]$read_quote)[, c('id', 'cell_type', deconv_data[[ct]]$genes_filtered), with = FALSE]
# 		sc_data <- melt(sc_data, id.vars = c('id', 'cell_type'), variable.name = 'gene', value.name = 'expression_level')
#
# 		# Centre genes:
# 		gene_averages <- sc_data[
# 			,
# 			.(ave_exp = mean(c(mean(expression_level[cell_type == 'cancer']), mean(expression_level[cell_type == 'caf'])))),
# 			by = .(symbol = gene)
# 		]
# 		sc_data[, expression_level := expression_level - gene_averages[symbol == unique(gene), ave_exp], by = gene]
#
# 		# Centre cells:
# 		sc_data[, expression_level := expression_level - mean(expression_level), by = id]
#
# 		sc_data[, .(d = mean(expression_level[cell_type == 'cancer']) - mean(expression_level[cell_type == 'caf'])), by = gene][, setNames(d, gene)]
#
# 	},
# 	simplify = FALSE,
# 	USE.NAMES = TRUE
# )
#
# for(ct in names(deconv_data)) {deconv_data[[ct]]$extra_data_score <- sc_comp_diff[[ct]]}





set.seed(8542)

sc_tcga_deconv_comp <- sapply(
	ct_to_keep,
	function(ct) {

		cat(paste0(ct, '\n'))

        sc_data <- eval(sc_metadata[[tcga_sc_map[ct]]]$read_quote)

		sc_deconv_comp <- sapply(
			c('cancer', 'caf'),
			function(x) {

				plot_data <- sc_data[cell_type == x][
					id %in% sc_cancer_caf[[tcga_sc_map[ct]]]$scores[id][order(-score), cell_id[1:sample_sizes[[tcga_sc_map[ct]]]]],
					c('id', with(deconv_data[[ct]], genes_filtered[genes_filtered %in% names(sc_data)])),
					with = FALSE
				]

				# plot_data <- copy(sc_data[cell_type == x])[
				# 	,
				# 	complexity := apply(.SD, 1, function(x) sum(x > 0)),
				# 	.SDcols = -c('id', 'patient', 'cell_type')
				# ][
				# 	,
				# 	.SD[sample(1:.N, sample_sizes[[tcga_sc_map[ct]]], prob = complexity)]
				# ][, complexity := NULL][, c('id', with(deconv_data[[ct]], genes_filtered[genes_filtered %in% names(sc_data)])), with = FALSE]

				plot_data <- melt(plot_data, id.vars = 'id', variable.name = 'gene', value.name = 'expression_level')

				ordered_cell_ids <- plot_data[
					gene %in% with(deconv_data[[ct]], do.call(c('head', 'tail')[which(c('cancer', 'caf') == x)], list(genes_filtered[ordering], 20))),
					.(top_20_mean = mean(expression_level)),
					by = id
				][order(top_20_mean), id]

				list(plot_data = plot_data, ordered_cell_ids = ordered_cell_ids)

			},
			simplify = FALSE,
			USE.NAMES = TRUE
		)

		sc_deconv_comp_data <- rbind(sc_deconv_comp$cancer$plot_data[, cell_type := 'cancer'], sc_deconv_comp$caf$plot_data[, cell_type := 'caf'])

		# Calculate averages in cancer cells and CAFs separately, so we can check them afterwards:
		gene_averages_cancer_caf <- sc_deconv_comp_data[, .(ave_exp = mean(expression_level)), by = .(gene, cell_type)]

		# To centre genes w.r.t. the average of the averages of cancer and CAF:
		gene_averages <- sc_deconv_comp_data[
			,
			.(ave_exp = mean(c(mean(expression_level[cell_type == 'cancer']), mean(expression_level[cell_type == 'caf'])))),
			by = .(symbol = gene)
		]
		sc_deconv_comp_data[, expression_level := expression_level - gene_averages[symbol == unique(gene), ave_exp], by = gene]

		# To centre the cells as well:
		sc_deconv_comp_data[, expression_level_cc := expression_level - mean(expression_level), by = id]

		# Apply running average to each cell:
		# sc_deconv_comp_data[
		# 	,
		# 	expression_level_cc_rm := setNames(
		# 		runmean(setNames(expression_level_cc, gene)[with(deconv_data[[ct]], genes_filtered[ordering])], 30),
		# 		with(deconv_data[[ct]], genes_filtered[ordering])
		# 	)[as.character(gene)], # Melting made gene a factor, and this reordering doesn't work unless we coerce it to character
		# 	by = id
		# ]

		sc_heatmaps <- sapply(
			c('expression_level', 'expression_level_cc'),
			function(expr_var) sapply(
				c('cancer', 'caf'),
				function(x) {
					ggplot(
						sc_deconv_comp_data[cell_type == x],
						aes(
							x = factor(gene, levels = with(deconv_data[[ct]], genes_filtered[ordering])),
							y = factor(id, levels = sc_deconv_comp[[x]]$ordered_cell_ids),
							fill = get(expr_var)
						)
					) +
						geom_raster() +
						scale_x_discrete(expand = c(0, 0)) +
						scale_y_discrete(expand = c(0, 0)) +
						scale_fill_gradientn(
							colours = c(
								sequential_hcl(25, h = 190, c = 70, l = c(70, 99), power = 0.8),
								sequential_hcl(25, h = 60, c = 100, l = c(80, 99), power = 0.8, rev = TRUE)
							),
							# colours = rev(colorRampPalette(brewer.pal(11, "RdBu"))(50)),
							limits = switch((expr_var == 'expression_level_cc_rm') + 1, c(-4, 4), c(-2, 2)),
							breaks = switch((expr_var == 'expression_level_cc_rm') + 1, c(-4, -2, 0, 2, 4), c(-2, -1, 0, 1, 2)),
							labels = switch(
								(expr_var == 'expression_level_cc_rm') + 1,
								c('-4' = '\u2264 -4', '-2' = '-2', '0' = '0', '2' = '2', '4' = '\u2265 4'),
								c('-2' = '\u2264 -2', '-1' = '-1', '0' = '0', '1' = '1', '2' = '\u2265 2')
							),
							# limits = c(-4, 4),
							# breaks = c(-4, -2, 0, 2, 4),
							# labels = c('-4' = '\u2264 -4', '-2' = '-2', '0' = '0', '2' = '2', '4' = '\u2265 4'),
							oob = squish
						) +
						theme(
							axis.text = element_blank(),
							axis.title.x = element_blank(),
							axis.ticks = element_blank(),
							axis.ticks.length = unit(0, 'pt'),
							plot.margin = unit(c(5.5, 5.5, 1.25, 5.5), 'pt')
							# plot.margin = switch((x == 'cancer') + 1, unit(c(1.25, 5.5, 5.5, 5.5), 'pt'), unit(c(5.5, 5.5, 1.25, 5.5), 'pt'))
						) +
						labs(y = mapvalues(x, c('cancer', 'caf'), c('Cancer', 'CAF'), warn_missing = FALSE), fill = 'Relative\nexpression\nlevel')
				},
				simplify = FALSE,
				USE.NAMES = TRUE
			),
			simplify = FALSE,
			USE.NAMES = TRUE
		)

		# Make versions where we filter out genes with low expression in the single cell data:

		deconv_filtered <- deconv_data[[ct]]

		filtered_genes <- sc_data[
			,
			sapply(.SD[cell_type == 'cancer'], sc_cancer_caf_args[[tcga_sc_map[ct]]]$genes_filter_fun) |
				sapply(.SD[cell_type == 'caf'], sc_cancer_caf_args[[tcga_sc_map[ct]]]$genes_filter_fun),
			.SDcols = with(deconv_data[[ct]], genes_filtered[genes_filtered %in% names(sc_data)])
		]
		filtered_genes <- names(filtered_genes)[filtered_genes]

		ordered_filtered_genes <- with(deconv_filtered, genes_filtered[ordering][genes_filtered[ordering] %in% filtered_genes])

		deconv_filtered$ordering <- order(filtered_genes)[order(order(ordered_filtered_genes))]
		deconv_filtered$genes_filtered <- filtered_genes
		deconv_filtered$cor_mat <- deconv_filtered$cor_mat[filtered_genes, filtered_genes]
		deconv_filtered$cor_with_purity <- sapply(deconv_filtered$cor_with_purity, function(x) x[filtered_genes], simplify = FALSE, USE.NAMES = TRUE)
		deconv_filtered$ccle_comp_diff <- deconv_filtered$ccle_comp_diff[filtered_genes]

		deconv_filtered_plots <- deconvolve_emt_caf_plots(
			data = deconv_filtered,
			heatmap_legend_title = 'Correlation',
			heatmap_colours = rev(colorRampPalette(brewer.pal(11, "RdBu"))(50)),
			heatmap_colour_limits = c(-0.6, 0.6),
		    heatmap_legend_breaks = c(-0.6, -0.3, 0, 0.3, 0.6),
			heatmap_annotations_side = 'left',
			heatmap_annotations_nudge = 0.3,
			purity_colours = rev(colorRampPalette(brewer.pal(11, "PuOr"))(50)),
			purity_colour_limits = c(-0.3, 0.3),
			purity_legend_breaks = c(-0.3, 0, 0.3),
			purity_legend_title = 'Correlation with purity\n',
			purity_legend_direction = 'horizontal',
			purity_axis_title = NULL,
			ccle_colours = rev(colorRampPalette(brewer.pal(11, 'PiYG')[c(1:3, 6, 9:11)])(50)),
			ccle_fun = function(x) {runmean(x - mean(x), 30)/max(abs(runmean(x - mean(x), 30)))},
			ccle_legend_breaks = c(-1, 0, 1),
			ccle_legend_title = 'Tumours vs. cell lines\n',
			ccle_legend_direction = 'horizontal',
			ccle_axis_title = NULL,
			bar_legend_justification = 'left',
			bar_legend_width = unit(10, 'pt'),
			bar_legend_height = unit(10, 'pt'),
			plot_title = deconv_titles[[ct]]
		)

		plot_data <- sc_deconv_comp_data[gene %in% deconv_filtered$genes_filtered]

		# Re-centre genes and cells w.r.t. the filtered gene list:
		gene_averages <- plot_data[
			,
			.(ave_exp = mean(c(mean(expression_level[cell_type == 'cancer']), mean(expression_level[cell_type == 'caf'])))),
			by = .(symbol = gene)
		]
		plot_data[, expression_level := expression_level - gene_averages[symbol == unique(gene), ave_exp], by = gene]
        plot_data[, expression_level_cc := expression_level - mean(expression_level), by = id]

		# Re-do running average per cell:
		# sc_deconv_comp_data[
		# 	,
		# 	expression_level_cc_rm := setNames(
		# 		runmean(setNames(expression_level_cc, gene)[with(deconv_data[[ct]], genes_filtered[ordering])], 30),
		# 		with(deconv_data[[ct]], genes_filtered[ordering])
		# 	)[as.character(gene)], # Melting made gene a factor, and this reordering doesn't work unless we coerce it to character
		# 	by = id
		# ]

		plot_data <- sapply(
			c('cancer', 'caf'),
			function(x) {
				plot_data[
					cell_type == x,
					.(
						id = factor(id, levels = sc_deconv_comp[[x]]$ordered_cell_ids),
						gene = factor(gene, levels = with(deconv_filtered, genes_filtered[ordering])),
						expression_level = expression_level,
						expression_level_cc = expression_level_cc
					)
				]
			},
			simplify = FALSE,
			USE.NAMES = TRUE
		)

		sc_heatmaps_filtered <- sapply(
			c('expression_level', 'expression_level_cc'),
			function(expr_var) sapply(
				c('cancer', 'caf'),
				function(x) {
					ggplot(
						plot_data[[x]],
						# plot_data[cell_type == x],
						aes(
							x = gene,
							y = id,
							# x = factor(gene, levels = with(deconv_filtered, genes_filtered[ordering])),
							# y = factor(id, levels = sc_deconv_comp[[x]]$ordered_cell_ids),
							fill = get(expr_var)
						)
					) +
						geom_raster() +
						scale_x_discrete(expand = c(0, 0)) +
						scale_y_discrete(expand = c(0, 0)) +
						scale_fill_gradientn(
							colours = c(
								sequential_hcl(25, h = 190, c = 70, l = c(70, 99), power = 0.8),
								sequential_hcl(25, h = 60, c = 100, l = c(80, 99), power = 0.8, rev = TRUE)
							),
							# colours = rev(colorRampPalette(brewer.pal(11, "RdBu"))(50)),
							limits = switch((expr_var == 'expression_level_cc_rm') + 1, c(-4, 4), c(-2, 2)),
							breaks = switch((expr_var == 'expression_level_cc_rm') + 1, c(-4, -2, 0, 2, 4), c(-2, -1, 0, 1, 2)),
							labels = switch(
								(expr_var == 'expression_level_cc_rm') + 1,
								c('-4' = '\u2264 -4', '-2' = '-2', '0' = '0', '2' = '2', '4' = '\u2265 4'),
								c('-2' = '\u2264 -2', '-1' = '-1', '0' = '0', '1' = '1', '2' = '\u2265 2')
							),
							# limits = c(-4, 4),
							# breaks = c(-4, -2, 0, 2, 4),
							# labels = c('-4' = '\u2264 -4', '-2' = '-2', '0' = '0', '2' = '2', '4' = '\u2265 4'),
							oob = squish
						) +
						theme(
							axis.text = element_blank(),
							axis.title.x = element_blank(),
							axis.ticks = element_blank(),
							axis.ticks.length = unit(0, 'pt'),
							plot.margin = unit(c(5.5, 5.5, 1.25, 0), 'pt')
						) +
						labs(y = mapvalues(x, c('cancer', 'caf'), c('Cancer', 'CAF'), warn_missing = FALSE), fill = 'Relative\nexpression\nlevel')
				},
				simplify = FALSE,
				USE.NAMES = TRUE
			),
			simplify = FALSE,
			USE.NAMES = TRUE
		)

		list(
			sc_deconv_comp_data = sc_deconv_comp_data,
			filtered_deconv_data = deconv_filtered,
			filtered_deconv_figures = deconv_filtered_plots$plots,
			sc_heatmaps_unfiltered = sc_heatmaps,
			sc_heatmaps_filtered = sc_heatmaps_filtered,
			sc_heatmaps_filtered_data = plot_data,
			gene_averages_cancer_caf = gene_averages_cancer_caf,
			gene_averages_for_centring = gene_averages,
			ordered_cell_ids <- sapply(sc_deconv_comp, `[[`, 'ordered_cell_ids', simplify = FALSE, USE.NAMES = TRUE)
		)

	},
	simplify = FALSE,
	USE.NAMES = TRUE
)

for(expr_var in c('expression_level', 'expression_level_cc')) {

	cairo_pdf(
		paste0('../data_and_figures/sc_deconv_comp_highscore_100', gsub('expression_level', '', expr_var), '.pdf'),
		width = 6,
		height = 9,
		onefile = TRUE
	)

	for(ct in names(sc_tcga_deconv_comp)) {

		plot_grid(
			plot_grid(
				plotlist = c(
					list(blank_plot() + theme(plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), 'pt')) + labs(title = deconv_titles[[ct]])),
					lapply(
						deconv_plots[[ct]]$plots[c('purity_bar', 'ccle_bar', 'heatmap')],
						function(x) x + theme(legend.position = 'none', plot.title = element_blank())
					),
					lapply(sc_tcga_deconv_comp[[ct]]$sc_heatmaps_unfiltered[[expr_var]], function(x) x + theme(legend.position = 'none'))
				),
				nrow = 6,
				ncol = 1,
				align = 'v',
				rel_heights = c(1, 1, 1, 15, 5, 5)
			),
			plot_grid(
				plotlist = c(
					list(
						blank_plot(),
						get_legend(
							deconv_plots[[1]]$plots$purity_bar +
								theme(
									legend.justification = c(0, 1),
									legend.direction = 'vertical',
									legend.key.width = NULL,
									legend.key.height = NULL
								) +
								labs(fill = 'Correlation\nwith purity')
						),
						get_legend(
							deconv_plots[[1]]$plots$ccle_bar +
								theme(
									legend.justification = c(0, 1),
									legend.direction = 'vertical',
									legend.key.width = NULL,
									legend.key.height = NULL
								) +
								labs(fill = 'Tumours vs.\ncell lines')
						),
						get_legend(deconv_plots[[1]]$plots$heatmap + theme(legend.justification = c(0, 1))),
						get_legend(sc_tcga_deconv_comp[[1]]$sc_heatmaps_unfiltered[[expr_var]][[1]] + theme(legend.justification = c(0, 1)))
					)
				),
				nrow = 5,
				ncol = 1,
				rel_heights = c(1, 5, 5, 7, 10)
			),
			nrow = 1,
			ncol = 2,
			rel_widths = c(5, 1)
		) %>% print

	}

	dev.off()

}

for(expr_var in c('expression_level', 'expression_level_cc')) {

	cairo_pdf(
		paste0('../data_and_figures/sc_deconv_comp_filtered_highscore_100', gsub('expression_level', '', expr_var), '.pdf'),
		width = 6,
		height = 9,
		onefile = TRUE
	)

	for(ct in names(sc_tcga_deconv_comp)) {

		plot_grid(
			plot_grid(
				plotlist = c(
					list(blank_plot() + theme(plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), 'pt')) + labs(title = deconv_titles[[ct]])),
					lapply(
						c(
							sc_tcga_deconv_comp[[ct]]$filtered_deconv_figures[c('purity_bar', 'ccle_bar', 'heatmap')],
							sc_tcga_deconv_comp[[ct]]$sc_heatmaps_filtered[[expr_var]]
						),
						function(x) x + theme(legend.position = 'none', plot.title = element_blank())
					)
				),
				nrow = 6,
				ncol = 1,
				align = 'v',
				rel_heights = c(1, 1, 1, 15, 5, 5)
			),
			plot_grid(
				plotlist = c(
					list(
						blank_plot(),
						get_legend(
							deconv_plots[[1]]$plots$purity_bar +
								theme(
									legend.justification = c(0, 1),
									legend.direction = 'vertical',
									legend.key.width = NULL,
									legend.key.height = NULL
								) +
								labs(fill = 'Correlation\nwith purity')
						),
						get_legend(
							deconv_plots[[1]]$plots$ccle_bar +
								theme(
									legend.justification = c(0, 1),
									legend.direction = 'vertical',
									legend.key.width = NULL,
									legend.key.height = NULL
								) +
								labs(fill = 'Tumours vs.\ncell lines')
						),
						get_legend(deconv_plots[[1]]$plots$heatmap + theme(legend.justification = c(0, 1))),
						get_legend(sc_tcga_deconv_comp[[1]]$sc_heatmaps_filtered[[expr_var]][[1]] + theme(legend.justification = c(0, 1)))
					)
				),
				nrow = 5,
				ncol = 1,
				rel_heights = c(1, 5, 5, 7, 10)
			),
			nrow = 1,
			ncol = 2,
			rel_widths = c(5, 1)
		) %>% print

	}

	dev.off()

}

# cairo_pdf('../data_and_figures/sc_deconv_comp_centred_temp.pdf', width = 6, height = 9, onefile = TRUE)
#
# for(ct in names(sc_tcga_deconv_comp)) {
#
# 	plot_grid(
# 		plot_grid(
# 			plotlist = c(
# 				list(blank_plot() + theme(plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), 'pt')) + labs(title = deconv_titles[[ct]])),
# 				lapply(
# 					deconv_plots[[ct]]$plots[c('purity_bar', 'ccle_bar', 'heatmap')],
# 					function(x) x + theme(legend.position = 'none', plot.title = element_blank())
# 				),
# 				lapply(sc_tcga_deconv_comp[[ct]]$sc_heatmaps_unfiltered, function(x) x + theme(legend.position = 'none'))
# 			),
# 			nrow = 6,
# 			ncol = 1,
# 			align = 'v',
# 			rel_heights = c(1, 1, 1, 15, 5, 5)
# 		),
# 		plot_grid(
# 			plotlist = c(
# 				list(
# 					blank_plot(),
# 					get_legend(
# 						deconv_plots[[1]]$plots$purity_bar +
# 							theme(legend.justification = c(0, 1), legend.direction = 'vertical', legend.key.width = NULL, legend.key.height = NULL) +
# 							labs(fill = 'Correlation\nwith purity')
# 					),
# 					get_legend(
# 						deconv_plots[[1]]$plots$ccle_bar +
# 							theme(legend.justification = c(0, 1), legend.direction = 'vertical', legend.key.width = NULL, legend.key.height = NULL) +
# 							labs(fill = 'Tumours vs.\ncell lines')
# 					),
# 					get_legend(deconv_plots[[1]]$plots$heatmap + theme(legend.justification = c(0, 1))),
# 					get_legend(sc_tcga_deconv_comp[[1]]$sc_heatmaps_unfiltered[[1]] + theme(legend.justification = c(0, 1)))
# 				)
# 			),
# 			nrow = 5,
# 			ncol = 1,
# 			rel_heights = c(1, 5, 5, 7, 10)
# 		),
# 		nrow = 1,
# 		ncol = 2,
# 		rel_widths = c(5, 1)
# 	) %>% print
#
# }
#
# dev.off()

# cairo_pdf('../data_and_figures/sc_deconv_comp_filtered_temp.pdf', width = 6, height = 9, onefile = TRUE)
#
# for(ct in names(sc_tcga_deconv_comp)) {
#
# 	plot_grid(
# 		plot_grid(
# 			plotlist = c(
# 				list(blank_plot() + theme(plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), 'pt')) + labs(title = deconv_titles[[ct]])),
# 				lapply(
# 					c(
# 						sc_tcga_deconv_comp[[ct]]$filtered_deconv_figures[c('purity_bar', 'ccle_bar', 'heatmap')],
# 						sc_tcga_deconv_comp[[ct]]$sc_heatmaps_filtered
# 					),
# 					function(x) x + theme(legend.position = 'none', plot.title = element_blank())
# 				)
# 			),
# 			nrow = 6,
# 			ncol = 1,
# 			align = 'v',
# 			rel_heights = c(1, 1, 1, 15, 5, 5)
# 		),
# 		plot_grid(
# 			plotlist = c(
# 				list(
# 					blank_plot(),
# 					get_legend(
# 						deconv_plots[[1]]$plots$purity_bar +
# 							theme(legend.justification = c(0, 1), legend.direction = 'vertical', legend.key.width = NULL, legend.key.height = NULL) +
# 							labs(fill = 'Correlation\nwith purity')
# 					),
# 					get_legend(
# 						deconv_plots[[1]]$plots$ccle_bar +
# 							theme(legend.justification = c(0, 1), legend.direction = 'vertical', legend.key.width = NULL, legend.key.height = NULL) +
# 							labs(fill = 'Tumours vs.\ncell lines')
# 					),
# 					get_legend(deconv_plots[[1]]$plots$heatmap + theme(legend.justification = c(0, 1))),
# 					get_legend(sc_tcga_deconv_comp[[1]]$sc_heatmaps_filtered[[1]] + theme(legend.justification = c(0, 1)))
# 				)
# 			),
# 			nrow = 5,
# 			ncol = 1,
# 			rel_heights = c(1, 5, 5, 7, 10)
# 		),
# 		nrow = 1,
# 		ncol = 2,
# 		rel_widths = c(5, 1)
# 	) %>% print
#
# }
#
# dev.off()
