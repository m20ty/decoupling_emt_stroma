# bsub -q tirosh -R rusage[mem=128000] -o sc_deconv_comp.o -e sc_deconv_comp.e Rscript sc_deconv_comp.R

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
library(grid) # 3.6.3
library(gtable) # 0.3.0
library(gridExtra) # 2.3

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





# sample_sizes <- list(
#     brca = 1000,
#     coadread = 800,
#     hnsc = 400,
#     lihc = 150,
#     luad = 1000,
# 	lusc = 100,
#     ov = 1500,
#     paad = 1500
# )

deconv_data <- readRDS('../data_and_figures/deconv_data.rds')
# deconv_plots <- readRDS('../data_and_figures/deconv_plots.rds')

cell_type_markers <- fread('../../cell_type_markers.csv')
subtypes_data <- fread('../../TCGA_data/tcga_subtypes_data.csv', key = 'id')
ccle <- fread('../../CCLE_cancer_type_Av.csv', key = 'gene_id')

ct_to_keep <- c('brca_luminal_a', 'brca_luminal_b', 'brca_basal_like', 'brca_her2_enriched', 'coad', 'hnsc_mesenchymal_basal', 'hnsc_classical',
	'luad_proximal_inflammatory', 'luad_proximal_proliferative', 'luad_terminal_respiratory_unit', 'lusc_classical', 'lusc_secretory',
	'ov_differentiated', 'ov_immunoreactive', 'ov_proliferative', 'paad', 'read')

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
	luad_terminal_respiratory_unit = 'luad',
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
    hnsc_mesenchymal_basal = 'HNSC - Malignant-Basal',
    hnsc_classical = 'HNSC - Classical',
    # hnsc_atypical = 'HNSC - Atypical',
    # lihc = 'LIHC',
    luad_proximal_inflammatory = 'LUAD - Squamoid',
    luad_proximal_proliferative = 'LUAD - Magnoid',
	luad_terminal_respiratory_unit = 'LUAD - Bronchioid',
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
			heatmap_annotations_side = 'left',
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
            # bar_legend_width = NULL,
            # bar_legend_height = NULL,
			# bar_legend_width = unit(10, 'pt'),
			bar_legend_height = unit(10, 'pt'),
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

				# plot_data <- copy(sc_data[cell_type == x])[
				# 	,
				# 	complexity := apply(.SD, 1, function(x) sum(x > 0)),
				# 	.SDcols = -c('id', 'patient', 'cell_type')
				# ][
				# 	,
				# 	.SD[sample(1:.N, sample_sizes[[tcga_sc_map[ct]]], prob = complexity)]
				# ][, complexity := NULL][, c('id', with(deconv_data[[ct]], genes_filtered[genes_filtered %in% names(sc_data)])), with = FALSE]

				plot_data <- sc_data[
					cell_type == x,
					c('id',  with(deconv_data[[ct]], genes_filtered[genes_filtered %in% names(sc_data)])),
					with = FALSE
				]

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
		sc_deconv_comp_data[
			,
			expression_level_cc_rm := setNames(
				runmean(setNames(expression_level_cc, gene)[with(deconv_data[[ct]], genes_filtered[ordering])], 30),
				with(deconv_data[[ct]], genes_filtered[ordering])
			)[as.character(gene)], # Melting made gene a factor, and this reordering doesn't work unless we coerce it to character
			by = id
		]

		sc_heatmaps <- sapply(
			c('expression_level', 'expression_level_cc', 'expression_level_cc_rm'),
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
							plot.margin = unit(c(5.5, 5.5, 1.25, 0), 'pt')
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
		# EDIT: is this working?  If we're using the genes from deconv_data[[ct]], surely it's not the filtered list?
		plot_data[
			,
			expression_level_cc_rm := setNames(
				runmean(setNames(expression_level_cc, gene)[with(deconv_filtered, genes_filtered[ordering])], 30),
				with(deconv_filtered, genes_filtered[ordering])
			)[as.character(gene)], # Melting made gene a factor, and this reordering doesn't work unless we coerce it to character
			by = id
		]

		plot_data <- sapply(
			c('cancer', 'caf'),
			function(x) {
				plot_data[
					cell_type == x,
					.(
						id = factor(id, levels = sc_deconv_comp[[x]]$ordered_cell_ids),
						gene = factor(gene, levels = with(deconv_filtered, genes_filtered[ordering])),
						expression_level = expression_level,
						expression_level_cc = expression_level_cc,
						expression_level_cc_rm = expression_level_cc_rm
					)
				]
			},
			simplify = FALSE,
			USE.NAMES = TRUE
		)

		sc_heatmaps_filtered <- sapply(
			c('expression_level', 'expression_level_cc', 'expression_level_cc_rm'),
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
			ordered_cell_ids = sapply(sc_deconv_comp, `[[`, 'ordered_cell_ids', simplify = FALSE, USE.NAMES = TRUE)
		)

	},
	simplify = FALSE,
	USE.NAMES = TRUE
)

for(expr_var in c('expression_level', 'expression_level_cc', 'expression_level_cc_rm')) {

	cairo_pdf(paste0('../data_and_figures/sc_deconv_comp', gsub('expression_level', '', expr_var), '.pdf'), width = 6, height = 9, onefile = TRUE)

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

for(expr_var in c('expression_level', 'expression_level_cc', 'expression_level_cc_rm')) {

	cairo_pdf(
		paste0('../data_and_figures/sc_deconv_comp_filtered', gsub('expression_level', '', expr_var), '.pdf'),
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





cairo_pdf('../data_and_figures/final_figures_resubmission/S10A.pdf', width = 11.5, height = 17.5, onefile = TRUE)

print(
	plot_grid(
		blank_plot(),
		plot_grid(
			plotlist = lapply(
				c('brca_luminal_b', 'brca_basal_like', 'brca_her2_enriched'),
				function(ct) {

					sc_deconv_comp_figures <- sapply( # Put deconv_plots[[ct]]$plots in the following...
						c(deconv_plots[[ct]]$plots, sc_tcga_deconv_comp[[ct]]$sc_heatmaps_unfiltered$expression_level_cc_rm),
						function(x) x + theme(legend.position = 'none', plot.title = element_blank()),
						simplify = FALSE,
						USE.NAMES = TRUE
					)
					sc_deconv_comp_figures$heatmap <- sc_deconv_comp_figures$heatmap +
						theme(plot.margin = unit(c(0, 20, 1.25, 20), 'pt')) +
						labs(y = '\nGenes')
					sc_deconv_comp_figures$cancer <- sc_deconv_comp_figures$cancer +
						labs(y = paste0('Cancer\ncells\n(n = ', length(sc_tcga_deconv_comp[[ct]]$ordered_cell_ids$cancer), ')'))
					sc_deconv_comp_figures$caf <- sc_deconv_comp_figures$caf +
						labs(y = paste0('\nCAFs\n(n = ', length(sc_tcga_deconv_comp[[ct]]$ordered_cell_ids$caf), ')'))

					plot_grid(
						plotlist = c(
							list(blank_plot() + theme(plot.margin = unit(c(11, 5.5, 5.5, 0), 'pt')) + labs(title = deconv_titles[[ct]])),
							sc_deconv_comp_figures[c('purity_bar', 'ccle_bar', 'heatmap', 'cancer', 'caf')],
							list(
								blank_plot() +
									scale_x_continuous(position = 'top') +
									theme(axis.title.x = element_text(), plot.margin = unit(c(4.5, 0, 0, 0), 'pt')) +
									labs(x = 'Genes')
							)
						),
						nrow = 7,
						ncol = 1,
						align = 'v',
						rel_heights = c(2, 1, 1, 15, 5, 5, 1)
					)

				}
			),
			nrow = 1,
			ncol = 3
		),
		blank_plot(),
		plot_grid(
			plotlist = lapply(
				c('coad', 'hnsc_classical', 'luad_proximal_inflammatory'),
				function(ct) {

					sc_deconv_comp_figures <- sapply(
						c(deconv_plots[[ct]]$plots, sc_tcga_deconv_comp[[ct]]$sc_heatmaps_unfiltered$expression_level_cc_rm),
						function(x) x + theme(legend.position = 'none', plot.title = element_blank()),
						simplify = FALSE,
						USE.NAMES = TRUE
					)
					sc_deconv_comp_figures$heatmap <- sc_deconv_comp_figures$heatmap +
						theme(plot.margin = unit(c(0, 20, 1.25, 20), 'pt'), axis.title.y = element_text()) +
						labs(y = '\nGenes')
					sc_deconv_comp_figures$cancer <- sc_deconv_comp_figures$cancer +
						labs(y = paste0('Cancer\ncells\n(n = ', length(sc_tcga_deconv_comp[[ct]]$ordered_cell_ids$cancer), ')'))
					sc_deconv_comp_figures$caf <- sc_deconv_comp_figures$caf +
						labs(y = paste0('\nCAFs\n(n = ', length(sc_tcga_deconv_comp[[ct]]$ordered_cell_ids$caf), ')'))

					plot_grid(
						plotlist = c(
							list(blank_plot() + theme(plot.margin = unit(c(11, 5.5, 5.5, 0), 'pt')) + labs(title = deconv_titles[[ct]])),
							sc_deconv_comp_figures[c('purity_bar', 'ccle_bar', 'heatmap', 'cancer', 'caf')],
							list(
								blank_plot() +
									scale_x_continuous(position = 'top') +
									theme(axis.title.x = element_text(), plot.margin = unit(c(4.5, 0, 0, 0), 'pt')) +
									labs(x = 'Genes')
							)
						),
						nrow = 7,
						ncol = 1,
						align = 'v',
						rel_heights = c(2, 1, 1, 15, 5, 5, 1)
					)

				}
			),
			nrow = 1,
			ncol = 3
		),
		blank_plot(),
		plot_grid(
			plotlist = c(
				lapply(
					c('luad_terminal_respiratory_unit', 'lusc_classical'),
					function(ct) {

						sc_deconv_comp_figures <- sapply(
							c(deconv_plots[[ct]]$plots, sc_tcga_deconv_comp[[ct]]$sc_heatmaps_unfiltered$expression_level_cc_rm),
							function(x) x + theme(legend.position = 'none', plot.title = element_blank()),
							simplify = FALSE,
							USE.NAMES = TRUE
						)
						sc_deconv_comp_figures$heatmap <- sc_deconv_comp_figures$heatmap +
							theme(plot.margin = unit(c(0, 20, 1.25, 20), 'pt'), axis.title.y = element_text()) +
							labs(y = '\nGenes')
						sc_deconv_comp_figures$cancer <- sc_deconv_comp_figures$cancer +
							labs(y = paste0('Cancer\ncells\n(n = ', length(sc_tcga_deconv_comp[[ct]]$ordered_cell_ids$cancer), ')'))
						sc_deconv_comp_figures$caf <- sc_deconv_comp_figures$caf +
							labs(y = paste0('\nCAFs\n(n = ', length(sc_tcga_deconv_comp[[ct]]$ordered_cell_ids$caf), ')'))

						plot_grid(
							plotlist = c(
								list(blank_plot() + theme(plot.margin = unit(c(11, 5.5, 5.5, 0), 'pt')) + labs(title = deconv_titles[[ct]])),
								sc_deconv_comp_figures[c('purity_bar', 'ccle_bar', 'heatmap', 'cancer', 'caf')],
								list(
									blank_plot() +
										scale_x_continuous(position = 'top') +
										theme(axis.title.x = element_text(), plot.margin = unit(c(4.5, 0, 0, 0), 'pt')) +
										labs(x = 'Genes')
								)
							),
							nrow = 7,
							ncol = 1,
							align = 'v',
							rel_heights = c(2, 1, 1, 15, 5, 5, 1)
						)

					}
				),
				list(
					plot_grid(
						blank_plot(),
						get_legend(
							deconv_plots$brca_luminal_a$plots$purity_bar +
								guides(fill = guide_colourbar(title.position = 'right')) +
								theme(
									legend.justification = c(0, 1),
									legend.direction = 'horizontal',
									legend.key.width = NULL,
									legend.box.margin = margin(l = 10)
								) +
								labs(fill = 'Correlation with purity\n')
						),
						get_legend(
							deconv_plots$brca_luminal_a$plots$ccle_bar +
								guides(fill = guide_colourbar(title.position = 'right')) +
								theme(
									legend.justification = c(0, 1),
									legend.direction = 'horizontal',
									legend.key.width = NULL,
									legend.box.margin = margin(l = 10)
								) +
								labs(fill = 'Tumours vs. cell lines\n')
						),
						get_legend(
							deconv_plots$brca_luminal_a$plots$heatmap +
								guides(fill = guide_colourbar(title.position = 'right')) +
								theme(
									legend.justification = c(0, 1),
									legend.direction = 'horizontal',
									legend.key.width = unit(25, 'pt'),
									legend.key.height = unit(10, 'pt'),
									legend.box.margin = margin(l = 10)
								) +
								labs(fill = 'Correlation\n')
						),
						get_legend(
							sc_tcga_deconv_comp$brca_luminal_a$sc_heatmaps_unfiltered$expression_level_cc_rm[[1]] +
								guides(fill = guide_colourbar(title.position = 'right')) +
								theme(
									legend.justification = c(0, 1),
									legend.direction = 'horizontal',
									legend.key.width = unit(25, 'pt'),
									legend.key.height = unit(10, 'pt'),
									legend.box.margin = margin(l = 10)
								) +
								labs(fill = 'Relative expression level\n')
						),
						blank_plot(),
						nrow = 6,
						ncol = 1,
						rel_heights = c(1, 1, 1, 1, 1, 3)
					)
				)
			),
			nrow = 1,
			ncol = 3
		),
		blank_plot(),
		nrow = 7,
		ncol = 1,
		rel_heights = c(0.2, 1.8, 0.15, 1.8, 0.15, 1.8, 0.125)
		# rel_heights = c(0.05, 1, 0.05, 1)
	) + draw_label('A', x = 0, y = 0.99, hjust = -0.1, vjust = 1.1, size = 20, fontface = 2)
)

print(
	plot_grid(
		blank_plot(),
		plot_grid(
			plotlist = lapply(
				c('lusc_secretory', 'ov_differentiated', 'ov_immunoreactive'),
				function(ct) {

					sc_deconv_comp_figures <- sapply(
						c(deconv_plots[[ct]]$plots, sc_tcga_deconv_comp[[ct]]$sc_heatmaps_unfiltered$expression_level_cc_rm),
						function(x) x + theme(legend.position = 'none', plot.title = element_blank()),
						simplify = FALSE,
						USE.NAMES = TRUE
					)
					sc_deconv_comp_figures$heatmap <- sc_deconv_comp_figures$heatmap +
						theme(plot.margin = unit(c(0, 20, 1.25, 20), 'pt'), axis.title.y = element_text()) +
						labs(y = '\nGenes')
					sc_deconv_comp_figures$cancer <- sc_deconv_comp_figures$cancer +
						labs(y = paste0('Cancer\ncells\n(n = ', length(sc_tcga_deconv_comp[[ct]]$ordered_cell_ids$cancer), ')'))
					sc_deconv_comp_figures$caf <- sc_deconv_comp_figures$caf +
						labs(y = paste0('\nCAFs\n(n = ', length(sc_tcga_deconv_comp[[ct]]$ordered_cell_ids$caf), ')'))

					plot_grid(
						plotlist = c(
							list(blank_plot() + theme(plot.margin = unit(c(11, 5.5, 5.5, 0), 'pt')) + labs(title = deconv_titles[[ct]])),
							sc_deconv_comp_figures[c('purity_bar', 'ccle_bar', 'heatmap', 'cancer', 'caf')],
							list(
								blank_plot() +
									scale_x_continuous(position = 'top') +
									theme(axis.title.x = element_text(), plot.margin = unit(c(4.5, 0, 0, 0), 'pt')) +
									labs(x = 'Genes')
							)
						),
						nrow = 7,
						ncol = 1,
						align = 'v',
						rel_heights = c(2, 1, 1, 15, 5, 5, 1)
					)

				}
			),
			nrow = 1,
			ncol = 3
		),
		blank_plot(),
		plot_grid(
			plotlist = lapply(
				c('ov_proliferative', 'paad', 'read'),
				function(ct) {

					sc_deconv_comp_figures <- sapply(
						c(deconv_plots[[ct]]$plots, sc_tcga_deconv_comp[[ct]]$sc_heatmaps_unfiltered$expression_level_cc_rm),
						function(x) x + theme(legend.position = 'none', plot.title = element_blank()),
						simplify = FALSE,
						USE.NAMES = TRUE
					)
					sc_deconv_comp_figures$heatmap <- sc_deconv_comp_figures$heatmap +
						theme(plot.margin = unit(c(0, 20, 1.25, 20), 'pt'), axis.title.y = element_text()) +
						labs(y = '\nGenes')
					sc_deconv_comp_figures$cancer <- sc_deconv_comp_figures$cancer +
						labs(y = paste0('Cancer\ncells\n(n = ', length(sc_tcga_deconv_comp[[ct]]$ordered_cell_ids$cancer), ')'))
					sc_deconv_comp_figures$caf <- sc_deconv_comp_figures$caf +
						labs(y = paste0('\nCAFs\n(n = ', length(sc_tcga_deconv_comp[[ct]]$ordered_cell_ids$caf), ')'))

					plot_grid(
						plotlist = c(
							list(blank_plot() + theme(plot.margin = unit(c(11, 5.5, 5.5, 0), 'pt')) + labs(title = deconv_titles[[ct]])),
							sc_deconv_comp_figures[c('purity_bar', 'ccle_bar', 'heatmap', 'cancer', 'caf')],
							list(
								blank_plot() +
									scale_x_continuous(position = 'top') +
									theme(axis.title.x = element_text(), plot.margin = unit(c(4.5, 0, 0, 0), 'pt')) +
									labs(x = 'Genes')
							)
						),
						nrow = 7,
						ncol = 1,
						align = 'v',
						rel_heights = c(2, 1, 1, 15, 5, 5, 1)
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
					theme(legend.justification = c(1, 1), legend.direction = 'horizontal', legend.key.width = NULL) +
					labs(fill = 'Correlation with purity\n')
			),
			blank_plot(),
			get_legend(
				deconv_plots$brca_luminal_a$plots$heatmap +
					guides(fill = guide_colourbar(title.position = 'right')) +
					theme(
						legend.justification = c(0, 1),
						legend.direction = 'horizontal',
						legend.key.width = unit(25, 'pt'),
						legend.key.height = unit(10, 'pt')
					) +
					labs(fill = 'Correlation\n')
			),
			get_legend(
				deconv_plots$brca_luminal_a$plots$ccle_bar +
					theme(legend.justification = c(1, 1), legend.direction = 'horizontal', legend.key.width = NULL) +
					labs(fill = 'Tumours vs. cell lines\n')
			),
			blank_plot(),
			get_legend(
				sc_tcga_deconv_comp$brca_luminal_a$sc_heatmaps_unfiltered$expression_level_cc_rm[[1]] +
					theme(
						legend.justification = c(0, 1),
						legend.direction = 'horizontal',
						legend.key.width = unit(25, 'pt'),
						legend.key.height = unit(10, 'pt')
					) +
					guides(fill = guide_colourbar(title.position = 'right')) +
					labs(fill = 'Relative expression level\n')
			),
			nrow = 2,
			ncol = 3,
			rel_widths = c(5, 1, 5)
		),
		blank_plot(),
		nrow = 7,
		ncol = 1,
		rel_heights = c(0.2, 1.8, 0.15, 1.8, 0.15, 0.35, 1.575)
		# rel_heights = c(0.05, 1, 0.05, 1)
	)
)

dev.off()





# Subset of plots for main figure and figure for rebuttal letter (so I can manually set the gene labels):

heatmap_annotations_subset <- lapply(
    list(
        brca_luminal_a = c('MPDZ', 'ITGB1', 'TGFBR3', 'CALU', 'COL1A1', 'COL1A2', 'THY1', 'FAP'),
        hnsc_mesenchymal_basal = c('LAMC2', 'SERPINE1', 'TGFBI', 'ITGB1', 'COL1A2', 'COL1A1', 'LUM', 'POSTN'),
        luad_proximal_proliferative = c('DST', 'AXL', 'AREG', 'LAMC2', 'FAP', 'LUM', 'POSTN', 'COL1A2')
    ),
    c,
    c('SNAI1', 'SNAI2', 'TWIST1', 'VIM', 'ZEB1', 'ZEB2')
)

deconv_plots_subset <- sapply(
	names(heatmap_annotations_subset),
	function(ct) deconvolve_emt_caf_plots(
		data = deconv_data[[ct]],
		# data = sc_tcga_deconv_comp[[ct]]$filtered_deconv_data,
		heatmap_legend_title = 'Correlation',
		heatmap_colours = rev(colorRampPalette(brewer.pal(11, "RdBu"))(50)),
		heatmap_colour_limits = c(-0.6, 0.6),
		heatmap_legend_breaks = c(-0.6, -0.3, 0, 0.3, 0.6),
		heatmap_annotations = heatmap_annotations_subset[[ct]],
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
	),
	simplify = FALSE,
	USE.NAMES = TRUE
)

pemt_caf_brackets_params <- list(
	brca_luminal_a = list(pemt_bracket_max = 42, caf_bracket_min = 71, pemt_label_hjust = 0.5, caf_label_hjust = 0.5),
	hnsc_mesenchymal_basal = list(pemt_bracket_max = 34, caf_bracket_min = 65, pemt_label_hjust = 0.5, caf_label_hjust = 0.5),
	luad_proximal_proliferative = list(pemt_bracket_max = 35, caf_bracket_min = 71, pemt_label_hjust = 0.5, caf_label_hjust = 0.5)
)





cairo_pdf('../data_and_figures/final_figures_resubmission/R1.pdf', width = 13, height = 6.7)

plot_grid(
	plot_grid(
		plotlist = lapply(
			names(heatmap_annotations_subset),
			function(ct) {

				sc_deconv_comp_figures <- sapply(
					c(deconv_plots_subset[[ct]]$plots, sc_tcga_deconv_comp[[ct]]$sc_heatmaps_unfiltered$expression_level_cc_rm),
					function(x) x + theme(legend.position = 'none', plot.title = element_blank(), axis.title.y = element_blank()),
					simplify = FALSE,
					USE.NAMES = TRUE
				)
				sc_deconv_comp_figures$heatmap <- sc_deconv_comp_figures$heatmap +
					theme(plot.margin = unit(c(0, 5.5, 1.25, 0), 'pt'))
				sc_deconv_comp_figures$axis_labels <- sc_deconv_comp_figures$axis_labels +
					theme(plot.margin = unit(c(0, 0, 1.25, 0), 'pt'))

				plot_grid(
					plot_grid(
						blank_plot(),
						pemt_caf_brackets(
							edge = 'right',
							pemt_bracket_max = pemt_caf_brackets_params[[ct]]$pemt_bracket_max,
							caf_bracket_min = pemt_caf_brackets_params[[ct]]$caf_bracket_min,
							brackets_just = 0.8,
							labels = c('pEMT genes', 'CAF genes'),
							pemt_label_hjust = pemt_caf_brackets_params[[ct]]$pemt_label_hjust,
							caf_label_hjust = pemt_caf_brackets_params[[ct]]$caf_label_hjust,
							labels_vjust = -0.7,
							labels_angle = 90,
							plot_margin = c(0, 0, 1.25, 0)
						),
						blank_plot(),
						nrow = 3,
						ncol = 1,
						rel_heights = c(4, 15, 11)
					),
					plot_grid(
						blank_plot(),
						sc_deconv_comp_figures$axis_labels,
						ggplot(data.frame(x = 0:1, y = 0:1), aes(x = x, y = y)) +
							geom_point(size = 0, colour = 'white') +
							geom_segment(
								aes(x = 0.65, xend = 0.65, y = 0.05, yend = 0.95),
								arrow = arrow(ends = 'both', length = unit(5, 'pt'))
							) +
							scale_x_continuous(expand = c(0, 0)) +
							scale_y_continuous(expand = c(0, 0)) +
							theme(
								axis.text = element_blank(),
								axis.title.x = element_blank(),
								axis.title.y = element_text(vjust = -3),
								axis.ticks = element_blank(),
								axis.ticks.length = unit(0, 'pt'),
								plot.background = element_blank(),
								panel.background = element_blank(),
								plot.margin = unit(c(5.5, 0, 1.25, 5.5), 'pt')
							) +
							labs(y = paste0('Cancer\ncells\n(n = ', length(sc_tcga_deconv_comp[[ct]]$ordered_cell_ids$cancer), ')')),
						ggplot(data.frame(x = 0:1, y = 0:1), aes(x = x, y = y)) +
							geom_point(size = 0, colour = 'white') +
							geom_segment(
								aes(x = 0.65, xend = 0.65, y = 0.05, yend = 0.95),
								arrow = arrow(ends = 'both', length = unit(5, 'pt'))
							) +
							scale_x_continuous(expand = c(0, 0)) +
							scale_y_continuous(expand = c(0, 0)) +
							theme(
								axis.text = element_blank(),
								axis.title.x = element_blank(),
								axis.title.y = element_text(vjust = -3),
								axis.ticks = element_blank(),
								axis.ticks.length = unit(0, 'pt'),
								plot.background = element_blank(),
								panel.background = element_blank(),
								plot.margin = unit(c(5.5, 0, 1.25, 5.5), 'pt')
							) +
							labs(y = paste0('\nCAFs\n(n = ', length(sc_tcga_deconv_comp[[ct]]$ordered_cell_ids$caf), ')')),
						# blank_plot() +
						# 	labs(y = 'Cancer\ncells') +
						# 	scale_y_continuous(position = 'right') +
						# 	theme(plot.margin = unit(c(5.5, 5.5, 1.25, 5.5), 'pt'), axis.title.y.right = element_text(angle = 90)),
						# blank_plot() +
						# 	labs(y = 'CAFs') +
						# 	scale_y_continuous(position = 'right') +
						# 	theme(plot.margin = unit(c(5.5, 5.5, 1.25, 5.5), 'pt'), axis.title.y.right = element_text(angle = 90)),
						blank_plot(),
						nrow = 5,
						ncol = 1,
						rel_heights = c(4, 15, 5, 5, 1)
					),
					plot_grid(
						plotlist = c(
							list(
								blank_plot() +
								theme(plot.margin = unit(c(11, 5.5, 5.5, 0), 'pt')) +
								labs(title = deconv_titles[[ct]])
							),
							sc_deconv_comp_figures[c('purity_bar', 'ccle_bar', 'heatmap', 'cancer', 'caf')],
							list(
								blank_plot() +
									scale_x_continuous(position = 'top') +
									theme(axis.title.x = element_text(), plot.margin = unit(c(4.5, 0, 0, 0), 'pt')) +
									labs(x = 'Genes')
							)
						),
						nrow = 7,
						ncol = 1,
						align = 'v',
						rel_heights = c(2, 1, 1, 15, 5, 5, 1)
					),
					nrow = 1,
					ncol = 3,
					rel_widths = c(0.8, 1.2, 4)
				)

			}
		),
		nrow = 1,
		ncol = 3
	),
	blank_plot(),
	plot_grid(
		get_legend(
			deconv_plots_subset$brca_luminal_a$plots$purity_bar +
				theme(legend.justification = c(1, 1), legend.direction = 'horizontal', legend.key.width = NULL) +
				labs(fill = 'Correlation with purity\n')
		),
		blank_plot(),
		get_legend(
			deconv_plots_subset$brca_luminal_a$plots$heatmap +
				guides(fill = guide_colourbar(title.position = 'right')) +
				# scale_fill_gradientn(
					# colours = rev(colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))(50)),
					# limits = c(-1, 1),
					# breaks = c(-1, 0, 1)
				# ) +
				theme(
					legend.justification = c(0, 1),
					legend.direction = 'horizontal',
					legend.key.width = unit(25, 'pt'),
					legend.key.height = unit(10, 'pt')
				) +
				labs(fill = 'Correlation\n')
		),
		get_legend(
			deconv_plots_subset$brca_luminal_a$plots$ccle_bar +
				theme(legend.justification = c(1, 1), legend.direction = 'horizontal', legend.key.width = NULL) +
				labs(fill = 'Tumours vs. cell lines\n')
		),
		blank_plot(),
		get_legend(
			sc_tcga_deconv_comp$brca_luminal_a$sc_heatmaps_unfiltered$expression_level_cc_rm[[1]] +
				theme(
					legend.justification = c(0, 1),
					legend.direction = 'horizontal',
					legend.key.width = unit(25, 'pt'),
					legend.key.height = unit(10, 'pt')
				) +
				guides(fill = guide_colourbar(title.position = 'right')) +
				labs(fill = 'Relative expression level\n')
		),
		nrow = 2,
		ncol = 3,
		rel_widths = c(5, 1, 5)
	),
	nrow = 3,
	ncol = 1,
	rel_heights = c(1, 0.05, 0.2)
) %>% print

dev.off()





# Scatterplot of within- versus between-cluster correlations:

ct_to_keep_summary <- c('blca_luminal_papillary', 'blca_basal_squamous', 'brca_luminal_a', 'brca_luminal_b', 'brca_basal_like', 'brca_her2_enriched',
    'coad', 'esca_ac', 'hnsc_mesenchymal_basal', 'hnsc_classical', 'luad_proximal_inflammatory', 'luad_proximal_proliferative',
    'luad_terminal_respiratory_unit', 'lusc_classical', 'lusc_secretory', 'ov_differentiated', 'ov_immunoreactive', 'ov_proliferative', 'paad',
    'read', 'stad_cin', 'stad_msi', 'ucec')
nice_names_summary <- c('BLCA - Luminal-Papillary', 'BLCA - Basal-Squamous', 'BRCA - Luminal A', 'BRCA - Luminal B', 'BRCA - Basal-like',
    'BRCA - HER2-enriched', 'COAD', 'ESCA - Adenocarcinoma', 'HNSC - Malignant-Basal', 'HNSC - Classical', 'LUAD - Squamoid', 'LUAD - Magnoid',
    'LUAD - Bronchioid', 'LUSC - Classical', 'LUSC - Secretory', 'OV - Differentiated', 'OV - Immunoreactive', 'OV - Proliferative', 'PAAD', 'READ',
    'STAD - CIN', 'STAD - MSI', 'UCEC')

within_between_clust_corr <- lapply(
	ct_to_keep_summary,
    function(ct) {
        within_clust_corr <- with(deconv_data[[ct]], c(cor_mat[ordering, ordering][1:30, 1:30], cor_mat[rev(ordering), rev(ordering)][1:30, 1:30]))
        within_clust_corr <- mean(within_clust_corr[within_clust_corr != 1])
        between_clust_corr <- with(deconv_data[[ct]], mean(cor_mat[ordering, rev(ordering)][1:30, 1:30]))
        list(cancer_type = ct, wthn = within_clust_corr, btw = between_clust_corr)
    }
) %>% rbindlist

within_between_clust_corr[
	,
	n_annot := sum(c('cor_with_purity', 'ccle_comp_diff', 'extra_data_score') %in% names(deconv_data[[cancer_type]])),
	by = cancer_type
][
	,
	c('scrnaseq_avail', 'cancer_type_nice', 'n_annot') := .(
		mapvalues(n_annot, c(2, 3), c('No', 'Yes')),
		mapvalues(cancer_type, ct_to_keep_summary, nice_names_summary),
		NULL
	)
]

deconv_summary_scatterplot <- ggplot(
    within_between_clust_corr,
    aes(x = btw, y = wthn, shape = scrnaseq_avail)
) +
    geom_point(size = 2, alpha = 0.75, colour = 'lightblue4') +
    # geom_text_repel(aes(label = cancer_type)) +
	geom_text_repel(
        aes(label = cancer_type_nice),
        data = within_between_clust_corr[cancer_type_nice == 'HNSC - Malignant-Basal'],
        nudge_x = 0.055, nudge_y = 0.01, size = 3.5
    ) +
	geom_text_repel(
        aes(label = cancer_type_nice),
        data = within_between_clust_corr[cancer_type_nice %in% c('HNSC - Classical', 'LUAD - Magnoid')],
        nudge_x = 0.03, nudge_y = 0.01, size = 3.5
    ) +
	geom_text_repel(
        aes(label = cancer_type_nice),
        data = within_between_clust_corr[cancer_type_nice == 'UCEC'],
        nudge_x = -0.025, size = 3.5
    ) +
    geom_text_repel(
        aes(label = cancer_type_nice),
        data = within_between_clust_corr[cancer_type_nice %in% c('READ', 'PAAD', 'COAD')],
        nudge_x = -0.02, nudge_y = -0.01, size = 3.5
    ) +
	geom_text_repel(
        aes(label = cancer_type_nice),
        data = within_between_clust_corr[cancer_type_nice == 'BRCA - Luminal A'],
        nudge_x = -0.05, nudge_y = -0.01, size = 3.5
    ) +
    theme_test() +
	theme(
		legend.position = 'bottom',
		legend.text = element_text(margin = margin(r = 15, unit = 'pt')),
		legend.title = element_text(margin = margin(r = 10, unit = 'pt')),
		legend.spacing.x = unit(1, 'pt')
	) +
    labs(x = 'Between-cluster correlation', y = 'Within-cluster correlation', shape = 'scRNA-seq data available:')

# Table providing a dictionary between TCGA disease codes and common cancer type names:

tcga_codes_table <- tableGrob(
    data.table(
        `TCGA code` = c('BLCA', 'BRCA', 'COAD', 'ESCA', 'HNSC', 'LUAD', 'LUSC', 'OV', 'PAAD', 'READ', 'STAD', 'UCEC'),
        `Cancer type` = c('Bladder', 'Breast', 'Colon', 'Oesophageal', 'Head and Neck', 'Lung (Adeno.)', 'Lung (Squamous)', 'Ovarian', 'Pancreatic',
            'Rectum', 'Stomach', 'Endometrial')
    ),
    theme = ttheme_minimal(
        padding = unit(c(10, 7.5), 'pt'),
        core = list(fg_params = list(x = 0.05, hjust = 0, fontsize = 11), padding = unit(c(15, 7.5), 'pt')),
        colhead = list(padding = unit(c(15, 15), 'pt'))
    ),
    rows = NULL
)

tcga_codes_table <- gtable_add_grob(
    tcga_codes_table,
    rectGrob(gp = gpar(fill = NA, lwd = 2)),
    t = 1, b = nrow(tcga_codes_table), l = 1, r = 2
)

tcga_codes_table <- gtable_add_grob(
    tcga_codes_table,
    segmentsGrob(x0 = unit(0, 'npc'), x1 = unit(0, 'npc'), y0 = unit(0, 'npc'), y1 = unit(1, 'npc'), gp = gpar(lwd = 2)),
    t = 1, b = nrow(tcga_codes_table), l = 2, r = 2
)

tcga_codes_table <- gtable_add_grob(
    tcga_codes_table,
    segmentsGrob(x0 = unit(0, 'npc'), x1 = unit(1, 'npc'), y0 = unit(0, 'npc'), y1 = unit(0, 'npc'), gp = gpar(lwd = 2)),
    t = 1, b = 1, l = 1, r = ncol(tcga_codes_table)
)

# To align to top:
tcga_codes_table$vp <- viewport(x = unit(0.5, 'npc'), y = unit(1, 'npc') - 0.5*sum(tcga_codes_table$heights))

# Summary of agreement with scRNA-seq data:

sc_deconv_comp_barplot_data <- rbindlist(
	lapply(
		names(sc_tcga_deconv_comp),
		function(ct) data.table(
			cancer_type = deconv_titles[[ct]],
			sc_diff = sc_tcga_deconv_comp[[ct]]$sc_deconv_comp_data[
				cell_type == 'cancer',
				.SD[gene %in% with(deconv_data[[ct]], head(genes_filtered[ordering], 20)), mean(expression_level_cc)] -
					.SD[gene %in% with(deconv_data[[ct]], tail(genes_filtered[ordering], 20)), mean(expression_level_cc)]
				# by = cell_type # No point doing it per cell type, because the centring means answer for one cell type is just minus the other
			]
		)
	)
)

sc_deconv_comp_barplot_data[, cancer_type := factor(cancer_type, levels = .SD[order(-sc_diff), cancer_type])]

sc_deconv_comp_barplot <- ggplot(sc_deconv_comp_barplot_data, aes(x = cancer_type, y = sc_diff)) +
	geom_col(colour = 'lightblue4', fill = 'lightblue3') +
	theme_test() +
	theme(axis.text.x = element_text(angle = 55, hjust = 1), axis.title.x = element_blank()) +
	labs(y = TeX('$\\Delta$(Expression level)'))

# sc_deconv_comp_scatterplot_data <- rbindlist(
# 	lapply(
# 		names(sc_tcga_deconv_comp),
# 		function(ct) cbind(
# 			cancer_type = deconv_titles[[ct]],
# 			merge(
# 				sc_tcga_deconv_comp[[ct]]$sc_deconv_comp_data[
# 					gene %in% with(deconv_data[[ct]], head(genes_filtered[ordering], 20)),
# 					.(pemt_sig = mean(expression_level_cc)),
# 					by = cell_type
# 				],
# 				sc_tcga_deconv_comp[[ct]]$sc_deconv_comp_data[
# 					gene %in% with(deconv_data[[ct]], tail(genes_filtered[ordering], 20)),
# 					.(caf_sig = mean(expression_level_cc)),
# 					by = cell_type
# 				]
# 			)
# 		)
# 	)
# )
#
# sc_deconv_comp_scatterplot_data <- dcast(
# 	melt(sc_deconv_comp_scatterplot_data, id.vars = c('cancer_type', 'cell_type'), variable.name = 'sig', value.name = 'expr'),
# 	cancer_type + sig ~ cell_type,
# 	value.var = 'expr'
# )
#
# sc_deconv_comp_scatterplot <- ggplot(sc_deconv_comp_scatterplot_data, aes(x = caf, y = cancer, shape = sig)) + geom_point() + theme_test()

cairo_pdf('../data_and_figures/final_figures_resubmission/3.pdf', width = 13, height = 12.7)

print(
	plot_grid(
		blank_plot(),
		plot_grid(
			plotlist = lapply(
				names(heatmap_annotations_subset),
				function(ct) {

					sc_deconv_comp_figures <- sapply(
						c(deconv_plots_subset[[ct]]$plots, sc_tcga_deconv_comp[[ct]]$sc_heatmaps_unfiltered$expression_level_cc_rm),
						function(x) x + theme(legend.position = 'none', plot.title = element_blank(), axis.title.y = element_blank()),
						simplify = FALSE,
						USE.NAMES = TRUE
					)
					sc_deconv_comp_figures$heatmap <- sc_deconv_comp_figures$heatmap +
						theme(plot.margin = unit(c(0, 5.5, 1.25, 0), 'pt'))
					sc_deconv_comp_figures$axis_labels <- sc_deconv_comp_figures$axis_labels +
						theme(plot.margin = unit(c(0, 0, 1.25, 0), 'pt'))

					plot_grid(
						plot_grid(
							blank_plot(),
							pemt_caf_brackets(
								edge = 'right',
								pemt_bracket_max = pemt_caf_brackets_params[[ct]]$pemt_bracket_max,
								caf_bracket_min = pemt_caf_brackets_params[[ct]]$caf_bracket_min,
								brackets_just = 0.8,
								labels = c('pEMT genes', 'CAF genes'),
								pemt_label_hjust = pemt_caf_brackets_params[[ct]]$pemt_label_hjust,
								caf_label_hjust = pemt_caf_brackets_params[[ct]]$caf_label_hjust,
								labels_vjust = -0.7,
								labels_angle = 90,
								plot_margin = c(0, 0, 1.25, 0)
							),
							blank_plot(),
							nrow = 3,
							ncol = 1,
							rel_heights = c(4, 15, 11)
						),
						plot_grid(
							blank_plot(),
							sc_deconv_comp_figures$axis_labels,
							ggplot(data.frame(x = 0:1, y = 0:1), aes(x = x, y = y)) +
								geom_point(size = 0, colour = 'white') +
								geom_segment(
									aes(x = 0.65, xend = 0.65, y = 0.05, yend = 0.95),
									arrow = arrow(ends = 'both', length = unit(5, 'pt'))
								) +
								scale_x_continuous(expand = c(0, 0)) +
								scale_y_continuous(expand = c(0, 0)) +
								theme(
									axis.text = element_blank(),
									axis.title.x = element_blank(),
									axis.title.y = element_text(vjust = -3),
									axis.ticks = element_blank(),
									axis.ticks.length = unit(0, 'pt'),
									plot.background = element_blank(),
									panel.background = element_blank(),
									plot.margin = unit(c(5.5, 0, 1.25, 5.5), 'pt')
								) +
								labs(y = paste0('Cancer\ncells\n(n = ', length(sc_tcga_deconv_comp[[ct]]$ordered_cell_ids$cancer), ')')),
							ggplot(data.frame(x = 0:1, y = 0:1), aes(x = x, y = y)) +
								geom_point(size = 0, colour = 'white') +
								geom_segment(
									aes(x = 0.65, xend = 0.65, y = 0.05, yend = 0.95),
									arrow = arrow(ends = 'both', length = unit(5, 'pt'))
								) +
								scale_x_continuous(expand = c(0, 0)) +
								scale_y_continuous(expand = c(0, 0)) +
								theme(
									axis.text = element_blank(),
									axis.title.x = element_blank(),
									axis.title.y = element_text(vjust = -3),
									axis.ticks = element_blank(),
									axis.ticks.length = unit(0, 'pt'),
									plot.background = element_blank(),
									panel.background = element_blank(),
									plot.margin = unit(c(5.5, 0, 1.25, 5.5), 'pt')
								) +
								labs(y = paste0('\nCAFs\n(n = ', length(sc_tcga_deconv_comp[[ct]]$ordered_cell_ids$caf), ')')),
							# blank_plot() +
							# 	labs(y = 'Cancer\ncells') +
							# 	scale_y_continuous(position = 'right') +
							# 	theme(plot.margin = unit(c(5.5, 5.5, 1.25, 5.5), 'pt'), axis.title.y.right = element_text(angle = 90)),
							# blank_plot() +
							# 	labs(y = 'CAFs') +
							# 	scale_y_continuous(position = 'right') +
							# 	theme(plot.margin = unit(c(5.5, 5.5, 1.25, 5.5), 'pt'), axis.title.y.right = element_text(angle = 90)),
							blank_plot(),
							nrow = 5,
							ncol = 1,
							rel_heights = c(4, 15, 5, 5, 1)
						),
						plot_grid(
							plotlist = c(
								list(
									blank_plot() +
									theme(plot.margin = unit(c(11, 5.5, 5.5, 0), 'pt')) +
									labs(title = deconv_titles[[ct]])
								),
								sc_deconv_comp_figures[c('purity_bar', 'ccle_bar', 'heatmap', 'cancer', 'caf')],
								list(
									blank_plot() +
										scale_x_continuous(position = 'top') +
										theme(axis.title.x = element_text(), plot.margin = unit(c(4.5, 0, 0, 0), 'pt')) +
										labs(x = 'Genes')
								)
							),
							nrow = 7,
							ncol = 1,
							align = 'v',
							rel_heights = c(2, 1, 1, 15, 5, 5, 1)
						),
						nrow = 1,
						ncol = 3,
						rel_widths = c(0.8, 1.2, 4)
					)

				}
			),
			nrow = 1,
			ncol = 3
		),
		blank_plot(),
		plot_grid(
			get_legend(
				deconv_plots_subset$brca_luminal_a$plots$purity_bar +
					theme(legend.justification = c(1, 1), legend.direction = 'horizontal', legend.key.width = NULL) +
					labs(fill = 'Correlation with purity\n')
			),
			blank_plot(),
			get_legend(
				deconv_plots_subset$brca_luminal_a$plots$heatmap +
					guides(fill = guide_colourbar(title.position = 'right')) +
					theme(
						legend.justification = c(0, 1),
						legend.direction = 'horizontal',
						legend.key.width = unit(25, 'pt'),
						legend.key.height = unit(10, 'pt')
					) +
					labs(fill = 'Correlation\n')
			),
			get_legend(
				deconv_plots_subset$brca_luminal_a$plots$ccle_bar +
					theme(legend.justification = c(1, 1), legend.direction = 'horizontal', legend.key.width = NULL) +
					labs(fill = 'Tumours vs. cell lines\n')
			),
			blank_plot(),
			get_legend(
				sc_tcga_deconv_comp$brca_luminal_a$sc_heatmaps_unfiltered$expression_level_cc_rm[[1]] +
					theme(
						legend.justification = c(0, 1),
						legend.direction = 'horizontal',
						legend.key.width = unit(25, 'pt'),
						legend.key.height = unit(10, 'pt')
					) +
					guides(fill = guide_colourbar(title.position = 'right')) +
					labs(fill = 'Relative expression level\n')
			),
			nrow = 2,
			ncol = 3,
			rel_widths = c(5, 1, 5)
		),
		blank_plot(),
		plot_grid(
			tcga_codes_table,
			sc_deconv_comp_barplot + theme(plot.margin = unit(c(0, 5.5, 5.5, 30), 'pt')),
			deconv_summary_scatterplot + theme(plot.margin = unit(c(0, 5.5, 5.5, 30), 'pt')),
			nrow = 1,
			ncol = 3,
			align = 'hv',
			axis = 'bl',
			rel_widths = c(2, 3, 3)
		),
		nrow = 6,
		ncol = 1,
		rel_heights = c(0.05, 1, 0.05, 0.2, 0.15, 0.95)
	) +
		draw_label('A', x = 0, y = 0.99, hjust = -0.1, vjust = 1.1, size = 20, fontface = 2) +
		draw_label('B', x = 0, y = 0.42, hjust = -0.1, size = 20, fontface = 2) +
		draw_label('C', x = 0.28, y = 0.42, size = 20, fontface = 2) +
		draw_label('D', x = 0.67, y = 0.42, size = 20, fontface = 2)
)

dev.off()
