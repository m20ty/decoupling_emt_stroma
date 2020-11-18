# bsub -q tirosh -R rusage[mem=64000] -o deconv_alt_v2.o -e deconv_alt_v2.e Rscript deconv_alt_v2.R

library(data.table) # 1.12.8
library(ggplot2) # This is version 3.3.1 on my laptop's WSL setup but 3.3.0 on WEXAC.
library(egg) # 0.4.5
library(limma) # 3.42.2
library(cowplot) # 1.0.0

source('general_functions.R')
source('tcga_functions.R')

# I need cell_type_markers for the gene filtering:
cell_type_markers <- fread('../../cell_type_markers.csv')

emt_markers <- fread('../../emt_markers.csv')[, gene := alias2SymbolTable(gene)][source != 'GO', sort(unique(gene))]

heatmap_annotations <- c('ACTA2', 'AXL', 'COL1A1', 'COL1A2', 'COL3A1', 'COL4A1', 'COL4A2', 'COL5A1', 'COL5A2', 'COL5A3', 'COL6A1', 'COL6A2',
						 'COL6A3', 'COL7A1', 'CD44', 'CDH2', 'ECM1', 'ECM2', 'FAP', 'FN1', 'IL6', 'ITGA2', 'ITGA5', 'ITGA6', 'ITGB1', 'ITGB3',
						 'ITGB5', 'ITGB6', 'LAMA1', 'LAMA2', 'LAMA3', 'LAMA5', 'LAMB3', 'LAMC1', 'LAMC2', 'MMP1', 'MMP2', 'MMP3', 'MMP10',
						 'MMP14', 'PDPN', 'SNAI1', 'SNAI2', 'SPARC', 'TGFB1', 'TGFBI', 'THY1', 'TNC', 'TWIST1', 'VCAN', 'VIM', 'ZEB1', 'ZEB2')

initial_genes <- c('SNAI1', 'SNAI2', 'TWIST1', 'VIM', 'ZEB1', 'ZEB2')





# Read in data:

expression_data <- fread('../../TCGA_data/tcga_expression_data.csv', key = 'id')
meta_data <- fread('../../TCGA_data/tcga_meta_data.csv', key = 'id')
subtypes_data <- fread('../../TCGA_data/tcga_subtypes_data_centred.csv', key = 'id')
ccle <- fread('../../CCLE_cancer_type_Av.csv', key = 'gene_id')
extra_data <- fread('../data_and_figures/collated_extra_data.csv', key = 'gene')

# I won't need the normal tissue samples, so let's take them out now:

expression_data <- expression_data[meta_data[sample_type != 'normal', id]]
meta_data <- meta_data[sample_type != 'normal']





# Note LUAD Proximal-Inflammatory = Squamoid; Proximal-Proliferative = Magnoid; and Terminal Respiratory Unit (TRU) = Bronchioid.

deconv_args_per_ct <- list(

    blca_luminal_infiltrated = list(
        tcga_cancer_types = 'BLCA',
        subtypes = 'Luminal_infiltrated',
        ref_for_subtypes = 'doi:10.1016/j.cell.2017.09.007',
        ccle_cancer_type = 'urinary tract',
        seed = 6007,
        plot_title = 'BLCA - Luminal-infiltrated',
		genes_filter_fun = function(x) 1:200,
		gene_weights_fun = function(x) quantile(x, 0.9)
    ),

    blca_luminal_papillary = list(
        tcga_cancer_types = 'BLCA',
        subtypes = 'Luminal_papillary',
        ref_for_subtypes = 'doi:10.1016/j.cell.2017.09.007',
        ccle_cancer_type = 'urinary tract',
        seed = 1582,
        plot_title = 'BLCA - Luminal-papillary',
		genes_filter_fun = function(x) 1:230,
		gene_weights_fun = median
    ),

    blca_basal_squamous = list(
        tcga_cancer_types = 'BLCA',
        subtypes = 'Basal_squamous',
        ref_for_subtypes = 'doi:10.1016/j.cell.2017.09.007',
        ccle_cancer_type = 'urinary tract',
        seed = 1027,
        plot_title = 'BLCA - Basal-squamous',
		genes_filter_fun = function(x) 1:250,
		gene_weights_fun = function(x) quantile(x, 0.9)
    ),

    brca_luminal_a = list(
        tcga_cancer_types = 'BRCA',
        subtypes = 'Luminal A',
        ref_for_subtypes = 'doi:10.1038/nature11412',
        ccle_cancer_type = 'breast',
        extra_data_source = 'breast_qian',
        seed = 8216,
        plot_title = 'BRCA - Luminal A',
		genes_filter_fun = function(x) 1:175,
		gene_weights_fun = function(x) quantile(x, 0.9)
    ),

    brca_luminal_b = list(
        tcga_cancer_types = 'BRCA',
        subtypes = 'Luminal B',
        ref_for_subtypes = 'doi:10.1038/nature11412',
        ccle_cancer_type = 'breast',
        extra_data_source = 'breast_qian',
        seed = 9056,
        plot_title = 'BRCA - Luminal B',
		genes_filter_fun = function(x) 1:230,
		gene_weights_fun = function(x) quantile(x, 0.9)
    ),

    brca_basal_like = list(
        tcga_cancer_types = 'BRCA',
        subtypes = 'Basal-like',
        ref_for_subtypes = 'doi:10.1038/nature11412',
        ccle_cancer_type = 'breast',
        extra_data_source = 'breast_qian',
        seed = 1982,
        plot_title = 'BRCA - Basal-like',
		genes_filter_fun = function(x) 1:200,
		gene_weights_fun = function(x) quantile(x, 0.9)
    ),

    brca_her2_enriched = list(
        tcga_cancer_types = 'BRCA',
        subtypes = 'HER2-enriched',
        ref_for_subtypes = 'doi:10.1038/nature11412',
        ccle_cancer_type = 'breast',
        extra_data_source = 'breast_qian',
        seed = 6994,
        plot_title = 'BRCA - HER2-enriched',
		genes_filter_fun = function(x) 1:250,
		gene_weights_fun = median
    ),

    cesc = list(
        tcga_cancer_types = 'CESC',
        seed = 4825,
        plot_title = 'CESC',
		genes_filter_fun = function(x) 1:150,
		gene_weights_fun = median
    ),

    coad = list(
        tcga_cancer_types = 'COAD',
        ccle_cancer_type = 'large intestine',
        extra_data_source = 'crc_lee_smc',
        seed = 3260,
        plot_title = 'COAD',
		genes_filter_fun = function(x) 1:160,
		gene_weights_fun = function(x) quantile(x, 0.9)
    ),

    esca_ac = list(
        tcga_cancer_types = 'ESCA',
        subtypes = 'AC',
        ref_for_subtypes = 'doi:10.1038/nature20805',
        ccle_cancer_type = 'oesophagus',
        seed = 2022,
        plot_title = 'ESCA - AC',
		genes_filter_fun = function(x) 1:200,
		gene_weights_fun = function(x) quantile(x, 0.9)
    ),

    esca_escc = list(
        tcga_cancer_types = 'ESCA',
        subtypes = 'ESCC',
        ref_for_subtypes = 'doi:10.1038/nature20805',
        ccle_cancer_type = 'oesophagus',
        seed = 6358,
        plot_title = 'ESCA - ESCC',
		genes_filter_fun = function(x) 1:250,
		gene_weights_fun = function(x) quantile(x, 0.9)
    ),

    hnsc_mesenchymal_basal = list(
        tcga_cancer_types = 'HNSC',
        subtypes = c('Mesenchymal', 'Basal'),
        ref_for_subtypes = 'doi:10.1038/nature14129',
        ccle_cancer_type = 'HNSCC',
        extra_data_source = 'hnscc_puram',
        seed = 7466,
        plot_title = 'HNSC - Mesenchymal & Basal',
		genes_filter_fun = function(x) 1:250,
		gene_weights_fun = function(x) quantile(x, 0.9)
    ),

    hnsc_classical = list(
        tcga_cancer_types = 'HNSC',
        subtypes = 'Classical',
        ref_for_subtypes = 'doi:10.1038/nature14129',
        ccle_cancer_type = 'HNSCC',
        extra_data_source = 'hnscc_puram',
        seed = 9998,
        plot_title = 'HNSC - Classical',
		genes_filter_fun = function(x) 1:200,
		gene_weights_fun = function(x) quantile(x, 0.9)
    ),

    hnsc_atypical = list(
        tcga_cancer_types = 'HNSC',
        subtypes = 'Atypical',
        ref_for_subtypes = 'doi:10.1038/nature14129',
        ccle_cancer_type = 'HNSCC',
        extra_data_source = 'hnscc_puram',
        seed = 8172,
        plot_title = 'HNSC - Atypical',
		genes_filter_fun = function(x) 1:200,
		gene_weights_fun = function(x) quantile(x, 0.9)
    ),

    kich = list(
        tcga_cancer_types = 'KICH',
        ccle_cancer_type = 'kidney',
        seed = 8796,
        plot_title = 'KICH',
		genes_filter_fun = function(x) 1:250,
		cell_type_weights = FALSE
    ),

    kirp = list(
        tcga_cancer_types = 'KIRP',
        ccle_cancer_type = 'kidney',
        seed = 3031,
        plot_title = 'KIRP',
		genes_filter_fun = function(x) 1:200,
		cell_type_weights = FALSE
    ),

    lihc = list(
        tcga_cancer_types = 'LIHC',
        ccle_cancer_type = 'liver',
        extra_data_source = 'liver_ma',
        seed = 9349,
        plot_title = 'LIHC',
		# genes_filter_fun = function(x) 1:260,
		genes_filter_fun = function(x) 1:230,
		cell_type_weights = FALSE
    ),

    luad_proximal_inflammatory = list(
        tcga_cancer_types = 'LUAD',
        subtypes = 'Proximal-inflammatory',
        ref_for_subtypes = 'doi:10.1038/nature13385',
        ccle_cancer_type = 'lung',
        extra_data_source = 'luad_kim',
        seed = 9955,
        plot_title = 'LUAD - Proximal-inflammatory',
		genes_filter_fun = function(x) 1:250,
		cell_type_weights = FALSE
    ),

    luad_proximal_proliferative = list(
        tcga_cancer_types = 'LUAD',
        subtypes = 'Proximal-proliferative',
        ref_for_subtypes = 'doi:10.1038/nature13385',
        ccle_cancer_type = 'lung',
        extra_data_source = 'luad_kim',
        seed = 743,
        plot_title = 'LUAD - Proximal-proliferative',
		genes_filter_fun = function(x) 1:250,
		gene_weights_fun = function(x) quantile(x, 0.9)
    ),

    luad_terminal_respiratory_unit = list(
        tcga_cancer_types = 'LUAD',
        subtypes = 'Terminal respiratory unit',
        ref_for_subtypes = 'doi:10.1038/nature13385',
        ccle_cancer_type = 'lung',
        extra_data_source = 'luad_kim',
        seed = 6026,
        plot_title = 'LUAD - Terminal respiratory unit',
		genes_filter_fun = function(x) 1:250,
		gene_weights_fun = median
    ),

	lusc_basal = list(
        tcga_cancer_types = 'LUSC',
        subtypes = 'basal',
        ref_for_subtypes = 'doi:10.1038/nature11404',
        ccle_cancer_type = 'lung',
        extra_data_source = 'lusc_qian',
        seed = 2287,
        plot_title = 'LUSC - Basal',
		genes_filter_fun = function(x) 1:190,
		gene_weights_fun = function(x) quantile(x, 0.9)
    ),

    lusc_classical = list(
        tcga_cancer_types = 'LUSC',
        subtypes = 'classical',
        ref_for_subtypes = 'doi:10.1038/nature11404',
        ccle_cancer_type = 'lung',
        extra_data_source = 'lusc_qian',
        seed = 7106,
        plot_title = 'LUSC - Classical',
		genes_filter_fun = function(x) 1:250,
		gene_weights_fun = median
    ),

    lusc_secretory = list(
        tcga_cancer_types = 'LUSC',
        subtypes = 'secretory',
        ref_for_subtypes = 'doi:10.1038/nature11404',
        ccle_cancer_type = 'lung',
        extra_data_source = 'lusc_qian',
        seed = 9240,
        plot_title = 'LUSC - Secretory',
		genes_filter_fun = function(x) 1:270,
		gene_weights_fun = function(x) quantile(x, 0.9)
    ),

    ov_differentiated = list(
        tcga_cancer_types = 'OV',
        subtypes = 'Differentiated',
        ref_for_subtypes = 'doi:10.1172/JCI65833',
        ccle_cancer_type = 'ovary',
		extra_data_source = 'ovarian_qian',
        seed = 1498,
        plot_title = 'OV - Differentiated',
		genes_filter_fun = function(x) 1:260,
		gene_weights_fun = function(x) quantile(x, 0.9)
    ),

    ov_immunoreactive = list(
        tcga_cancer_types = 'OV',
        subtypes = 'Immunoreactive',
        ref_for_subtypes = 'doi:10.1172/JCI65833',
        ccle_cancer_type = 'ovary',
		extra_data_source = 'ovarian_qian',
        seed = 991,
        plot_title = 'OV - Immunoreactive',
		genes_filter_fun = function(x) 1:200,
		gene_weights_fun = function(x) quantile(x, 0.9)
    ),

    ov_mesenchymal = list(
        tcga_cancer_types = 'OV',
        subtypes = 'Mesenchymal',
        ref_for_subtypes = 'doi:10.1172/JCI65833',
        ccle_cancer_type = 'ovary',
		extra_data_source = 'ovarian_qian',
        seed = 7549,
        plot_title = 'OV - Mesenchymal',
		genes_filter_fun = function(x) 1:200,
		gene_weights_fun = median
    ),

    ov_proliferative = list(
        tcga_cancer_types = 'OV',
        subtypes = 'Proliferative',
        ref_for_subtypes = 'doi:10.1172/JCI65833',
        ccle_cancer_type = 'ovary',
		extra_data_source = 'ovarian_qian',
        seed = 8583,
        plot_title = 'OV - Proliferative',
		genes_filter_fun = function(x) 1:200,
		cell_type_weights = FALSE
    ),

    paad_basal_moffitt = list(
        tcga_cancer_types = 'PAAD',
        subtypes = 'Basal',
        ref_for_subtypes = 'doi:10.1038/ng.3398',
        ccle_cancer_type = 'pancreas',
        extra_data_source = 'pdac_peng',
        seed = 7839,
        plot_title = 'PAAD - Basal',
		genes_filter_fun = function(x) 1:150,
		gene_weights_fun = median
    ),

    paad_classical_moffitt = list(
        tcga_cancer_types = 'PAAD',
        subtypes = 'Classical',
        ref_for_subtypes = 'doi:10.1038/ng.3398',
        ccle_cancer_type = 'pancreas',
        extra_data_source = 'pdac_peng',
        seed = 7794,
        plot_title = 'PAAD - Classical',
		genes_filter_fun = function(x) 1:160,
		gene_weights_fun = function(x) quantile(x, 0.9)
    ),

    prad = list(
        tcga_cancer_types = 'PRAD',
        ccle_cancer_type = 'prostate',
        seed = 6204,
        plot_title = 'PRAD',
		genes_filter_fun = function(x) 1:250,
		gene_weights_fun = function(x) quantile(x, 0.9)
    ),

    read = list(
        tcga_cancer_types = 'READ',
        ccle_cancer_type = 'large intestine',
        seed = 3309,
        plot_title = 'READ',
        extra_data_source = 'crc_lee_smc',
		genes_filter_fun = function(x) 1:200,
		gene_weights_fun = function(x) quantile(x, 0.9)
    ),

	stad_cin = list(
        tcga_cancer_types = 'STAD',
        subtypes = 'CIN',
        ref_for_subtypes = 'doi:10.1038/nature20805 - STAD',
        ccle_cancer_type = 'stomach',
        seed = 2773,
        plot_title = 'STAD - CIN',
        genes_filter_fun = function(x) 1:150,
        cell_type_weights = list(B_plasma = 0, myocyte = 0, macrophage = 0, endothelial = 0, DC = 0, mast = 0, T = 0, B = 1)
    ),

    stad_ebv = list(
        tcga_cancer_types = 'STAD',
        subtypes = 'EBV',
        ref_for_subtypes = 'doi:10.1038/nature20805 - STAD',
        ccle_cancer_type = 'stomach',
        seed = 3407,
        plot_title = 'STAD - EBV',
		genes_filter_fun = function(x) 1:150,
		gene_weights_fun = median
    ),

	stad_gs = list(
        tcga_cancer_types = 'STAD',
        subtypes = 'GS',
        ref_for_subtypes = 'doi:10.1038/nature20805 - STAD',
        ccle_cancer_type = 'stomach',
        seed = 5017,
        plot_title = 'STAD - GS',
		genes_filter_fun = function(x) 1:150,
		cell_type_weights = FALSE
    ),

    stad_msi = list(
        tcga_cancer_types = 'STAD',
        subtypes = 'MSI',
        ref_for_subtypes = 'doi:10.1038/nature20805 - STAD',
        ccle_cancer_type = 'stomach',
        seed = 6868,
        plot_title = 'STAD - MSI',
		genes_filter_fun = function(x) 1:180, # 220 also works well
        cell_type_weights = list(B_plasma = 0, myocyte = 0, macrophage = 0, endothelial = 0, DC = 0, mast = 0, T = 0, B = 1)
    ),

    ucec = list(
        tcga_cancer_types = 'UCEC',
        ccle_cancer_type = 'endometrium',
        seed = 475,
        plot_title = 'UCEC',
		genes_filter_fun = function(x) 1:200,
		gene_weights_fun = function(x) quantile(x, 0.9)
    )

)





deconv_data <- sapply(
	names(deconv_args_per_ct),
	function(ct) {
		cat(paste0(ct, '\n'))
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
					initial_gene_weights = FALSE
				),
				deconv_args_per_ct[[ct]][!(names(deconv_args_per_ct[[ct]]) == 'plot_title')]
			)
		)
		deconv_ct
	},
	simplify = FALSE,
	USE.NAMES = TRUE
)

# Reorder genes in deconv results:
deconv_data <- sapply(deconv_data, deconv_reorder, simplify = FALSE, USE.NAMES = TRUE)

deconv_plots <- sapply(
	names(deconv_args_per_ct),
	function(ct) {
		cat(paste0(ct, '\n'))
		deconv_ct <- deconv_data[[ct]]
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
				deconv_args_per_ct[[ct]]['plot_title']
			)
		)
		deconv_plot_ct
	},
	simplify = FALSE,
	USE.NAMES = TRUE
)

saveRDS(deconv_data, '../data_and_figures/deconv_alt_v2_data.rds')
saveRDS(deconv_plots, '../data_and_figures/deconv_alt_v2_plots.rds')





# Write all figures to a single PDF:
pdf('../data_and_figures/deconv_alt_v2_figures.pdf', width = 10, height = 12)
for(deconv_plot_ct in deconv_plots) {print(deconv_plot(list(deconv_plot_ct), legends_space = 0.2))}
dev.off()

# Diagnostic figures:
# pdf('../data_and_figures/deconv_alt_v2_diagnostics.pdf', width = 10, height = 12)
# i <- 1
# for(deconv_plot_ct in deconv_plots) {
# 	ggarrange(
# 		plots = c(list(deconv_plot_ct$plots$heatmap), deconv_plot_ct$diagnostics$alternative_purity_cor, deconv_plot_ct$diagnostics$cell_type_bars),
# 		ncol = 1,
# 		nrow = 11,
# 		heights = c(8, rep(0.6, 10)),
# 		newpage = switch((i == 1) + 1, TRUE, FALSE)
# 	)
# 	i <- i + 1
# }
# dev.off()

# Diagnostic figures:
pdf('../data_and_figures/deconv_alt_v2_diagnostics.pdf', width = 10, height = 12)
for(ct in names(deconv_data)) {

	deconv_ct <- deconv_data[[ct]]
	deconv_plot_ct <- deconv_plots[[ct]]

	regression_slopes <- with(
		deconv_ct,
		cor_with_initial_and_cell_types[genes_filtered][
			ordering,
			sapply(.SD, function(ct) setNames(lm(ct ~ I(.I/.N))$coeff['I(.I/.N)'], NULL), USE.NAMES = TRUE),
			.SDcols = cell_types
		]
	)

	print(
		plot_grid(
			plot_grid(
				plotlist = c(
					list(deconv_plot_ct$plots$heatmap + theme(legend.position = 'none')),
					lapply(deconv_plot_ct$diagnostics$alternative_purity_cor, function(g) {g + theme(legend.position = 'none')}),
					lapply(deconv_plot_ct$diagnostics$cell_type_bars, function(g) {g + theme(legend.position = 'none')}),
					list(blank_plot())
				),
				nrow = 12,
				ncol = 1,
				align = 'v',
				rel_heights = c(8, rep(0.6, 11))
			),
			plot_grid(
				plotlist = c(
					list(
						get_legend(deconv_plot_ct$plots$heatmap),
						get_legend(
							deconv_plot_ct$diagnostics$alternative_purity_cor[[1]] +
								theme(legend.title = element_text()) +
								guides(fill = guide_colorbar(title.position = 'top')) +
								labs(fill = 'Purity')
						),
						get_legend(
							deconv_plot_ct$diagnostics$cell_type_bars[[1]] +
								theme(legend.title = element_text()) +
								guides(fill = guide_colorbar(title.position = 'top')) +
								labs(fill = 'Cell types')
						),
						blank_plot()
					),
					lapply(
						deconv_ct$cell_types,
						function(x) {
							ggplot(data = data.table(x = 0:1, y = 0:1), aes(x = x, y = y)) +
								annotate(geom = 'text', x = 0.5, y = 0.5, label = round(regression_slopes[x], 2)) +
								theme(
						            axis.text = element_blank(),
						            axis.ticks = element_blank(),
						            axis.ticks.length = unit(0, 'pt'),
						            axis.title = element_blank(),
						            plot.background = element_blank(),
						            panel.background = element_blank()
						        )
						}
					),
					list(
						ggplot(data = data.table(x = 0:1, y = 0:1), aes(x = x, y = y)) +
							annotate(geom = 'text', x = 0.5, y = 0.5, label = paste('Min:', round(min(regression_slopes), 2))) +
							theme(
								axis.text = element_blank(),
								axis.ticks = element_blank(),
								axis.ticks.length = unit(0, 'pt'),
								axis.title = element_blank(),
								plot.background = element_blank(),
								panel.background = element_blank()
							)
					)
				),
				nrow = 13,
				ncol = 1,
				rel_heights = c(4, 2, 2, 1.2, rep(0.6, 9))
			),
			nrow = 1,
			ncol = 2,
			rel_widths = c(5, 1)
		)
	)

}
dev.off()
