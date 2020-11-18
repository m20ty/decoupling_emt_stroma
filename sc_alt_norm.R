# bsub -q tirosh -R rusage[mem=64000] -o sc_alt_norm.o -e sc_alt_norm.e Rscript sc_alt_norm.R

# Note there's a minor problem with the code at the end of this script that makes the final plot.  It won't affect the output if the whole
# script is run at once, but should probably be corrected.

library(data.table) # 1.12.8
library(ggplot2) # This is version 3.3.1 on my laptop's WSL setup but 3.3.0 on WEXAC.
library(Matrix) # 1.2.18
library(stringr) # 1.4.0
library(plyr) # 1.8.6
library(limma) # 3.42.2
library(magrittr) # 1.5
library(cowplot) # 1.0.0
library(RColorBrewer) # 1.1.2

source('general_functions.R')
source('sparse_matrix_functions.R')
source('sc_functions.R')
source('tcga_functions.R')

# Get TCGA data:
expression_data <- fread('../../TCGA_data/tcga_expression_data.csv', key = 'id')
meta_data <- fread('../../TCGA_data/tcga_meta_data.csv', key = 'id')
expression_data <- expression_data[meta_data[sample_type != 'normal', id]]
meta_data <- meta_data[sample_type != 'normal']

# Get EMT markers:
emt_markers <- fread('../../emt_markers.csv')[, gene := alias2SymbolTable(gene)][source != 'GO', sort(unique(gene))]





sc_metadata <- list(
	brca = list(
		ref = 'qian_breast_2020',
		tcga_cancer_types = 'BRCA',
		seed = 9181,
		initial_cell_types = c('b_cell', 'endothelial', 'macrophage', 'mast', 't_cell', 'caf', 'cancer'),
		rare_cell_types = 'dendritic',
		genes_filter_fun = function(x) {mean(x) >= 0.14 | sum(x >= 3.5) >= length(x)/100},
		scores_filter_fun = function(x) {mean(x) >= 0.37 | sum(x >= 3.5) >= length(x)/100},
		annotations = c('COL1A2', 'ECM1', 'FN1', 'NOTCH2', 'RHOB', 'SNAI1', 'SNAI2', 'THY1', 'TWIST1', 'VCAN', 'VEGFA', 'VIM', 'ZEB1', 'ZEB2'),
		annotations_title = 'Breast',
		max_mean_count = 800,
		deconv_args = list(
			tcga_cancer_types = 'brca',
			ccle_cancer_type = 'breast',
			seed = 1087,
			genes_filter_fun = function(x) 1:200,
			plot_title = 'Breast',
			heatmap_annotations = c('COL1A2', 'COL5A2', 'ECM1', 'FN1', 'NOTCH2', 'RHOB', 'SNAI1', 'SNAI2', 'TWIST1', 'VCAN', 'VEGFA', 'VIM', 'ZEB1', 'ZEB2')
		),
		sample_size = 1000,
		sim_deconv_filter = 0.14,
		# sim_deconv_filter = 0.1,
		pemt_bracket_max = 27,
		caf_bracket_min = 80
	),
	coadread = list(
		ref = 'lee_crc_2020_smc',
		tcga_cancer_types = c('COAD', 'READ'),
		seed = 9275,
		initial_cell_types = c('b_cell', 'endothelial', 'macrophage', 't_cell', 'caf', 'cancer'),
		rare_cell_types = c('epithelial', 'mast'),
		genes_filter_fun = function(x) {mean(x) >= 0.17 | sum(x >= 4.25) >= length(x)/100},
		scores_filter_fun = function(x) {mean(x) >= 0.46 | sum(x >= 4.25) >= length(x)/100},
		annotations = c('AREG', 'CALU', 'COL1A1', 'CXCL1', 'FN1', 'GADD45B', 'SNAI1', 'SNAI2', 'THY1', 'TWIST1', 'VIM', 'ZEB1', 'ZEB2'),
		annotations_title = 'Colorectal',
		annotations_side = 'left',
		max_mean_count = 1000,
		deconv_args = list(
			tcga_cancer_types = 'coadread',
			ccle_cancer_type = 'large intestine',
			seed = 6185,
			genes_filter_fun = function(x) 1:200,
			plot_title = 'Colorectal',
			heatmap_annotations = c('AREG', 'COL1A1', 'COL3A1', 'CXCL1', 'DCN', 'FN1', 'GADD45A', 'PLOD3', 'SNAI1', 'SNAI2', 'TWIST1', 'VIM', 'ZEB1', 'ZEB2')
		),
		sample_size = 800,
		sim_deconv_filter = 0.17,
		# sim_deconv_filter = 0.1,
		pemt_bracket_max = 27,
		caf_bracket_min = 80
	),
	luad = list(
		ref = 'kim_luad_2020',
		tcga_cancer_types = 'LUAD',
		seed = 8837,
		initial_cell_types = c('b_cell', 'endothelial', 'macrophage', 'mast', 't_cell', 'caf', 'cancer'),
		rare_cell_types = 'dendritic',
		genes_filter_fun = function(x) {mean(x) >= 0.12 | sum(x >= 3.5) >= length(x)/100},
		scores_filter_fun = function(x) {mean(x) >= 0.32 | sum(x >= 3.5) >= length(x)/100},
		annotations = c('AREG', 'CALU', 'COL1A1', 'FN1', 'RHOB', 'SDC1', 'SNAI1', 'SNAI2', 'SPARC', 'TWIST1', 'VEGFA', 'VIM', 'ZEB1', 'ZEB2'),
		annotations_title = 'Lung Adenocarcinoma',
		annotations_side = 'right',
		max_mean_count = 1000,
		deconv_args = list(
			tcga_cancer_types = 'luad',
			ccle_cancer_type = 'lung',
			seed = 5395,
			genes_filter_fun = function(x) 1:200,
			plot_title = 'Lung Adenocarcinoma',
			heatmap_annotations = c('CALU', 'COL1A2', 'COL3A1', 'MMP2', 'QSOX1', 'SDC4', 'SNAI1', 'SNAI2', 'THY1', 'TWIST1', 'VEGFA', 'VIM', 'ZEB1', 'ZEB2')
		),
		sample_size = 1000,
		sim_deconv_filter = 0.12,
		# sim_deconv_filter = 0.1,
		pemt_bracket_max = 27,
		caf_bracket_min = 80
	),
    paad = list(
		ref = 'peng_pdac_2019',
        tcga_cancer_types = 'PAAD',
		seed = 7466,
		initial_cell_types = c('b_cell', 'endothelial', 'macrophage', 't_cell', 'caf', 'cancer'),
		rare_cell_types = c('acinar', 'ductal_2', 'endocrine'),
		genes_filter_fun = function(x) {mean(x) >= 0.18 | sum(x >= 3.75) >= length(x)/100},
		scores_filter_fun = function(x) {mean(x) >= 0.48 | sum(x >= 3.75) >= length(x)/100},
		annotations = c('COL1A1', 'COL1A2', 'FN1', 'LAMC2', 'QSOX1', 'SDC4', 'SNAI1', 'SNAI2', 'TWIST1', 'VEGFA', 'VIM', 'ZEB1', 'ZEB2'),
        annotations_title = 'Pancreatic',
		annotations_side = 'left',
		max_mean_count = 1000,
		deconv_args = list(
			tcga_cancer_types = 'paad',
			ccle_cancer_type = 'pancreas',
			seed = 303,
			genes_filter_fun = function(x) 1:200,
			plot_title = 'Pancreatic',
			heatmap_annotations = c('COL1A2', 'COL3A1', 'DCN', 'FAP', 'LAMC2', 'QSOX1', 'SDC4', 'SNAI1', 'SNAI2', 'TWIST1', 'VEGFA', 'VIM', 'ZEB1', 'ZEB2')
		),
		sample_size = 1500,
		sim_deconv_filter = 0.18,
		# sim_deconv_filter = 0.08,
		pemt_bracket_max = 27,
		caf_bracket_min = 80
    )
)





all_figures <- list()

for(ct in names(sc_metadata)) {
	
	sc_data <- readRDS(paste0('../data_and_figures/sc_alt_norm/', sc_metadata[[ct]]$ref, '_scran.rds'))
	sc_meta <- fread(paste0('../data_and_figures/sc_alt_norm/', sc_metadata[[ct]]$ref, '_meta.csv'))
	
	sc_data <- cbind(sc_meta, t(as.matrix(sc_data)))
	rm(sc_meta)
	
	if(paste0('genes_filtered_scran_', sc_metadata[[ct]]$ref, '.rds') %in% dir('../data_and_figures/sc_alt_norm')) {
		genes_filtered_scran <- readRDS(paste0('../data_and_figures/sc_alt_norm/genes_filtered_scran_', sc_metadata[[ct]]$ref, '.rds'))
	} else {
		genes_unfiltered <- unique(
			c(
				emt_markers[emt_markers %in% names(sc_data)],
				top_cols_by_fun_cor(
					expression_data[meta_data[cancer_type %in% sc_metadata[[ct]]$tcga_cancer_types, id], -'id'],
					threshold_fun = function(x) quantile(x, 0.99)
				)[id %in% names(sc_data), id]
			)
		)
		genes_filtered_scran <- filter_for_groups(sc_data[, c('cell_type', ..genes_unfiltered)], groups = c('caf', 'cancer'))
		saveRDS(genes_filtered_scran, paste0('../data_and_figures/sc_alt_norm/genes_filtered_scran_', sc_metadata[[ct]]$ref, '.rds'))
	}
	
	# We'll also want to try with the gene list defined using the TPM data:
	genes_filtered_tpm <- readRDS('../data_and_figures/sc_genes_list.rds')[[ct]]$filtered_cancer_caf
	
	set.seed(sc_metadata[[ct]]$seed)
	
	if(paste0('sc_cancer_caf_scran_', sc_metadata[[ct]]$ref, '.rds') %in% dir('../data_and_figures/sc_alt_norm')) {
		sc_cancer_caf_scran <- readRDS(paste0('../data_and_figures/sc_alt_norm/sc_cancer_caf_scran_', sc_metadata[[ct]]$ref, '.rds'))
	} else {
		sc_cancer_caf_scran <- sc_groups(
			genes = genes_filtered_scran,
			sc_data = sc_data[cell_type %in% c('cancer', 'caf')],
			groups = c('cancer', 'caf'),
			score_cells_nbin = 30,
			score_cells_n = 40,
			min_sig_size = 0,
			scores_filter_groups = 'cancer',
			genes_filter_fun = sc_metadata[[ct]]$genes_filter_fun,
			scores_filter_fun = sc_metadata[[ct]]$scores_filter_fun
		)
		saveRDS(sc_cancer_caf_scran, paste0('../data_and_figures/sc_alt_norm/sc_cancer_caf_scran_', sc_metadata[[ct]]$ref, '.rds'))
	}
	
	sc_cancer_caf_heatmap_scran <- sc_groups_heatmap(
		sc_groups_list = sc_cancer_caf_scran,
		groups = c('cancer', 'caf'),
		x_axis_titles = c('Cancer cells', 'CAFs'),
		default_figure_widths = list(annotations = 2.2, cancer = 6, caf = 1.5),
		figure_spacing = 2.5,
		annotations_nudge = 0.25,
		es_fun = NULL,
		es_title = 'EMT\nscore',
		h_limits = c(0, 8),
		h_legend_breaks = c(0, 2, 4, 6, 8),
		h_legend_labels = c('0' = '0', '2' = '2', '4' = '4', '6' = '6', '8' = '\u2265 8'),
		h_legend_title = 'Expression level\n',
		h_legend_width = 20,
		h_legend_height = 10,
		h_legend_direction = 'horizontal',
		h_legend_title_position = 'right',
		h_legend_just = 'left',
		gd_legend_title = 'Genes detected\n',
		gd_legend_width = 20,
		gd_legend_height = 10,
		gd_legend_direction = 'horizontal',
		gd_legend_title_position = 'left',
		gd_legend_just = 'right',
		annotations = sc_metadata[[ct]]$annotations,
		title = sc_metadata[[ct]]$annotations_title
	)
	
	if('annotations_side' %in% names(sc_metadata[[ct]])) {
		sc_cancer_caf_heatmap_combining <- sc_groups_heatmap(
			sc_groups_list = sc_cancer_caf_scran,
			groups = c('cancer', 'caf'),
			x_axis_titles = c('Cancer cells', 'CAFs'),
			default_figure_widths = list(annotations = 2.5, cancer = 6, caf = 1.2),
			figure_spacing = 2.5,
			annotations_nudge = 0.25,
			es_fun = NULL,
			es_title = 'EMT score',
			h_limits = c(0, 8),
			h_legend_breaks = c(0, 2, 4, 6, 8),
			h_legend_labels = c('0' = '0', '2' = '2', '4' = '4', '6' = '6', '8' = '\u2265 8'),
			annotations = sc_metadata[[ct]]$annotations,
			annotations_title = sc_metadata[[ct]]$annotations_title,
			annotations_side = sc_metadata[[ct]]$annotations_side
		)
	}
	
	# Now with the genes from the TPM analysis:
	
	if(paste0('sc_cancer_caf_tpm_', sc_metadata[[ct]]$ref, '.rds') %in% dir('../data_and_figures/sc_alt_norm')) {
		sc_cancer_caf_tpm <- readRDS(paste0('../data_and_figures/sc_alt_norm/sc_cancer_caf_tpm_', sc_metadata[[ct]]$ref, '.rds'))
	} else {
		sc_cancer_caf_tpm <- sc_groups(
			genes = genes_filtered_tpm,
			sc_data = sc_data[cell_type %in% c('cancer', 'caf')],
			groups = c('cancer', 'caf'),
			score_cells_nbin = 30,
			score_cells_n = 40,
			min_sig_size = 0,
			scores_filter_groups = 'cancer',
			genes_filter_fun = sc_metadata[[ct]]$genes_filter_fun,
			scores_filter_fun = sc_metadata[[ct]]$scores_filter_fun
		)
		saveRDS(sc_cancer_caf_tpm, paste0('../data_and_figures/sc_alt_norm/sc_cancer_caf_tpm_', sc_metadata[[ct]]$ref, '.rds'))
	}
	
	sc_cancer_caf_heatmap_tpm <- sc_groups_heatmap(
		sc_groups_list = sc_cancer_caf_tpm,
		groups = c('cancer', 'caf'),
		x_axis_titles = c('Cancer cells', 'CAFs'),
		default_figure_widths = list(annotations = 2.5, cancer = 6, caf = 1.2),
		figure_spacing = 2.5,
		annotations_nudge = 0.25,
		es_fun = NULL,
		es_title = 'EMT\nscore',
		h_limits = c(0, 8),
		h_legend_breaks = c(0, 2, 4, 6, 8),
		h_legend_labels = c('0' = '0', '2' = '2', '4' = '4', '6' = '6', '8' = '\u2265 8'),
		annotations = sc_metadata[[ct]]$annotations,
		title = sc_metadata[[ct]]$annotations_title
	)
	
	sc_data[, cell_type := mapvalues(cell_type, sc_metadata[[ct]]$rare_cell_types, rep('rare', length(sc_metadata[[ct]]$rare_cell_types)))]
	
	if(paste0('lineplot_scran_', sc_metadata[[ct]]$ref, '.rds') %in% dir('../data_and_figures/sc_alt_norm')) {
		lineplot_scran <- readRDS(paste0('../data_and_figures/sc_alt_norm/lineplot_scran_', sc_metadata[[ct]]$ref, '.rds'))
	} else {
		lineplot_scran <- simulated_tumours_lineplot(
			sc_data,
			genes_filtered_scran,
			initial_types = sc_metadata[[ct]]$initial_cell_types,
			plot_title = sc_metadata[[ct]]$annotations_title,
			max_mean_count = sc_metadata[[ct]]$max_mean_count,
			legend_labels = c(
				'b_cell' = 'B cell',
				'cancer' = 'Cancer',
				'endothelial' = 'Endothelial',
				'caf' = 'CAF',
				'macrophage' = 'Macrophage',
				'mast' = 'Mast',
				't_cell' = 'T cell'
			),
			legend_colours = c(
				'b_cell' = '#8DD3C7',
				'cancer' = '#FB8072',
				'endothelial' = '#BC80BD',
				'caf' = '#FDB462',
				'macrophage' = '#80B1D3',
				'mast' = '#FCCDE5',
				't_cell' = '#B3DE69'
			),
		)
		saveRDS(lineplot_scran, paste0('../data_and_figures/sc_alt_norm/lineplot_scran_', sc_metadata[[ct]]$ref, '.rds'))
	}
	
	if(paste0('lineplot_tpm_', sc_metadata[[ct]]$ref, '.rds') %in% dir('../data_and_figures/sc_alt_norm')) {
		lineplot_tpm <- readRDS(paste0('../data_and_figures/sc_alt_norm/lineplot_tpm_', sc_metadata[[ct]]$ref, '.rds'))
	} else {
		lineplot_tpm <- simulated_tumours_lineplot(
			sc_data,
			genes_filtered_tpm,
			initial_types = sc_metadata[[ct]]$initial_cell_types,
			plot_title = sc_metadata[[ct]]$annotations_title,
			max_mean_count = sc_metadata[[ct]]$max_mean_count,
			legend_labels = c(
				'b_cell' = 'B cell',
				'cancer' = 'Cancer',
				'endothelial' = 'Endothelial',
				'caf' = 'CAF',
				'macrophage' = 'Macrophage',
				'mast' = 'Mast',
				't_cell' = 'T cell'
			),
			legend_colours = c(
				'b_cell' = '#8DD3C7',
				'cancer' = '#FB8072',
				'endothelial' = '#BC80BD',
				'caf' = '#FDB462',
				'macrophage' = '#80B1D3',
				'mast' = '#FCCDE5',
				't_cell' = '#B3DE69'
			),
		)
		saveRDS(lineplot_tpm, paste0('../data_and_figures/sc_alt_norm/lineplot_tpm_', sc_metadata[[ct]]$ref, '.rds'))
	}
	
	dummy_legend_plot <- ggplot(
		data = data.table(
			x = 1:7,
			y = 1,
			f = factor(
				c('b_cell', 'cancer', 'endothelial', 'caf', 'macrophage', 'mast', 't_cell'),
				levels = c('b_cell', 'cancer', 'endothelial', 'caf', 'macrophage', 'mast', 't_cell')
			)
		)
	) + 
		geom_tile(aes(x = x, y = y, fill = f)) +
		scale_fill_manual(
			labels = c(
				'b_cell' = 'B cell',
				'cancer' = 'Cancer',
				'endothelial' = 'Endothelial',
				'caf' = 'CAF',
				'macrophage' = 'Macrophage',
				'mast' = 'Mast',
				't_cell' = 'T cell'
			),
			values = c(
				'b_cell' = '#8DD3C7',
				'cancer' = '#FB8072',
				'endothelial' = '#BC80BD',
				'caf' = '#FDB462',
				'macrophage' = '#80B1D3',
				'mast' = '#FCCDE5',
				't_cell' = '#B3DE69'
			)
		) +
		labs(fill = 'Cell type') +
		theme(legend.key = element_rect(size = 1, colour = 'white'), legend.key.size = unit(15, 'pt'))
	
	dummy_legend <- get_legend(dummy_legend_plot)
	
	if(paste0('simulated_bulk_data_', sc_metadata[[ct]]$ref, '.csv') %in% dir('../data_and_figures/sc_alt_norm')) {
		simulated_bulk_data <- fread(paste0('../data_and_figures/sc_alt_norm/simulated_bulk_data_', sc_metadata[[ct]]$ref, '.csv'), key = 'id')
		simulated_bulk_metadata <- fread(paste0('../data_and_figures/sc_alt_norm/simulated_bulk_metadata_', sc_metadata[[ct]]$ref, '.csv'), key = 'id')
	} else {
		simulated_bulk_data <- simulated_tumours_data(
			as.matrix(sc_data[, -c('id', 'patient', 'cell_type')]),
			types = sc_data$cell_type,
			id_prefix = ct,
			max_mean_count = sc_data[cell_type == 'cancer', .N, by = patient][, round(quantile(N, 0.9))]
			# max_mean_count = max_mean_counts[[ct]]
			# genes = genes_list[[ct]]
		)
		simulated_bulk_metadata <- as.data.table(simulated_bulk_data$meta_data, keep.rownames = 'id')
		simulated_bulk_data <- as.data.table(simulated_bulk_data$expression_data, keep.rownames = 'id')
		simulated_bulk_metadata[, cancer_type := gsub('[0-9]+', '', id)]
		setcolorder(simulated_bulk_metadata, c('id', 'cancer_type'))
		setkey(simulated_bulk_data, 'id')
		setkey(simulated_bulk_metadata, 'id')
		fwrite(simulated_bulk_data, paste0('../data_and_figures/sc_alt_norm/simulated_bulk_data_', sc_metadata[[ct]]$ref, '.csv'))
		fwrite(simulated_bulk_metadata, paste0('../data_and_figures/sc_alt_norm/simulated_bulk_metadata_', sc_metadata[[ct]]$ref, '.csv'))
	}
	
	ccle <- fread('../../CCLE_cancer_type_Av.csv', key = 'gene_id')
	cell_type_markers <- fread('../../cell_type_markers.csv')
	
	if(paste0('simulated_deconv_', sc_metadata[[ct]]$ref, '.rds') %in% dir('../data_and_figures/sc_alt_norm')) {
		simulated_deconv <- readRDS(paste0('../data_and_figures/sc_alt_norm/simulated_deconv_', sc_metadata[[ct]]$ref, '.rds'))
	} else {
		
		simulated_deconv <- do.call(
			deconvolve_emt_caf_data,
			args = c(
				list(
					expression_data = simulated_bulk_data,
					meta_data = simulated_bulk_metadata,
					genes = emt_markers,
					cell_type_markers = cell_type_markers,
					ccle_data = ccle,
					initial_gene_weights = FALSE
				),
				sc_metadata[[ct]]$deconv_args[!(names(sc_metadata[[ct]]$deconv_args) %in% c('plot_title', 'heatmap_annotations'))]
			)
		)
		
		# Reorder genes in deconv results:
		simulated_deconv <- deconv_reorder(simulated_deconv)
		
		# Add data for single cell comparison colour bar:
		
		sc_data_melted <- melt(
			sc_data[, c('id', 'cell_type', simulated_deconv$genes_filtered), with = FALSE],
			id.vars = c('id', 'cell_type'),
			variable.name = 'gene',
			value.name = 'expression_level'
		)
		
		# Centre genes:
		gene_averages <- sc_data_melted[
			,
			.(ave_exp = mean(c(mean(expression_level[cell_type == 'cancer']), mean(expression_level[cell_type == 'caf'])))),
			by = .(symbol = gene)
		]
		sc_data_melted[, expression_level := expression_level - gene_averages[symbol == unique(gene), ave_exp], by = gene]
		
		# Centre cells:
		sc_data_melted[, expression_level := expression_level - mean(expression_level), by = id]
		
		simulated_deconv$extra_data_score <- sc_data_melted[
			,
			.(d = mean(expression_level[cell_type == 'cancer']) - mean(expression_level[cell_type == 'caf'])),
			by = gene
		][, setNames(d, gene)]
		
		saveRDS(simulated_deconv, paste0('../data_and_figures/sc_alt_norm/simulated_deconv_', sc_metadata[[ct]]$ref, '.rds'))
		
	}
	
	simulated_deconv_plots <- do.call(
		deconvolve_emt_caf_plots,
		args = c(
			list(
				data = simulated_deconv,
				heatmap_legend_title = 'Correlation',
				heatmap_colours = rev(colorRampPalette(brewer.pal(11, "RdBu"))(50)),
				heatmap_colour_limits = c(-1, 1),
				heatmap_legend_breaks = c(-1, -0.5, 0, 0.5, 1),
				heatmap_legend_justification = 'left',
				heatmap_annotations_nudge = 0.3,
				purity_colours = rev(colorRampPalette(brewer.pal(11, "PuOr"))(50)),
				purity_colour_limits = c(-1, 1),
				purity_legend_breaks = c(-1, 0, 1),
				purity_legend_title = 'Correlation with purity\n',
				purity_legend_direction = 'horizontal',
				purity_axis_title = NULL,
				ccle_colours = rev(colorRampPalette(brewer.pal(11, 'PiYG')[c(1:3, 6, 9:11)])(50)),
				ccle_fun = function(x) {caTools::runmean(x - mean(x), 30)/max(abs(caTools::runmean(x - mean(x), 30)))},
				ccle_legend_breaks = c(-1, 0, 1),
				ccle_legend_title = 'Tumours vs. cell lines\n',
				ccle_legend_direction = 'horizontal',
				ccle_axis_title = NULL,
				# extra_colours = c(
					# colorspace::sequential_hcl(25, h = 190, c = 70, l = c(70, 99), power = 0.8),
					# colorspace::sequential_hcl(25, h = 60, c = 100, l = c(80, 99), power = 0.8, rev = TRUE)
				# ),
				extra_colours = colorRampPalette(c('turquoise4', 'turquoise', 'azure', 'gold2', 'gold4'))(50),
				extra_fun = function(x) caTools::runmean(x, 30),
				extra_colour_limits = c(-4, 4),
				extra_legend_breaks = c(-4, 0, 4),
				extra_legend_labels = c('-4' = '\u2264 -4', '0' = '0', '4' = '\u2265 4'),
				extra_legend_title = 'scRNA-seq: CAF vs. cancer\n',
				extra_legend_direction = 'horizontal',
				extra_axis_title = NULL,
				bar_legend_justification = 'left'
				# bar_legend_width = NULL,
				# bar_legend_height = NULL
			),
			sc_metadata[[ct]]$deconv_args[names(sc_metadata[[ct]]$deconv_args) %in% c('plot_title', 'heatmap_annotations')]
		)
	)
	
	sc_deconv_comp <- sapply(
		c('cancer', 'caf'),
		function(x) {
			plot_data <- copy(sc_data[cell_type == x])[, complexity := apply(.SD, 1, function(x) sum(x > 0)), .SDcols = -c('id', 'patient', 'cell_type')][
				,
				.SD[sample(1:.N, sc_metadata[[ct]]$sample_size, prob = complexity)]
			][, complexity := NULL][, c('id', simulated_deconv$genes_filtered), with = FALSE]
			plot_data <- melt(plot_data, id.vars = 'id', variable.name = 'gene', value.name = 'expression_level')
			ordered_cell_ids <- plot_data[
				gene %in% with(
					simulated_deconv,
					do.call(c('head', 'tail')[which(c('cancer', 'caf') == x)], list(genes_filtered[ordering], 20))
				),
				.(top_20_mean = mean(expression_level)),
				by = id
			][order(top_20_mean), id]
			list(plot_data = plot_data, ordered_cell_ids = ordered_cell_ids)
		},
		simplify = FALSE,
		USE.NAMES = TRUE
	)
	
	sc_deconv_comp_data <- rbind(
		sc_deconv_comp$cancer$plot_data[, cell_type := 'cancer'],
		sc_deconv_comp$caf$plot_data[, cell_type := 'caf']
	)
	
	# Calculate averages in cancer cells and CAFs separately, so we can check them afterwards:
	gene_averages_cancer_caf <- sc_deconv_comp_data[
		,
		.(ave_exp = mean(expression_level)),
		by = .(gene, cell_type)
	]
	
	# To centre genes w.r.t. the average of the averages of cancer and CAF:
	gene_averages <- sc_deconv_comp_data[
		,
		.(ave_exp = mean(c(mean(expression_level[cell_type == 'cancer']), mean(expression_level[cell_type == 'caf'])))),
		by = .(symbol = gene)
	]
	sc_deconv_comp_data[, expression_level := expression_level - gene_averages[symbol == unique(gene), ave_exp], by = gene]
	
	# To centre the cells as well:
	sc_deconv_comp_data[, expression_level := expression_level - mean(expression_level), by = id]
	
	# simulated_deconv_filtered <- do.call(
		# deconvolve_emt_caf_data,
		# args = c(
			# list(
				# expression_data = simulated_bulk_data,
				# meta_data = simulated_bulk_metadata,
				# genes = gene_averages_cancer_caf[
					# ,
					# .(pass = ave_exp[cell_type == 'cancer'] > sc_metadata[[ct]]$sim_deconv_filter |
						# ave_exp[cell_type == 'caf'] > sc_metadata[[ct]]$sim_deconv_filter),
					# by = gene
				# ][pass == TRUE, as.character(gene)],
				# cell_type_markers = cell_type_markers,
				# ccle_data = ccle,
				# genes_from_tcga_fun = NULL,
				# genes_filter_fun = NULL
			# ),
			# sc_metadata[[ct]]$deconv_args[!(names(sc_metadata[[ct]]$deconv_args) %in% c('genes_filter_fun', 'plot_title', 'heatmap_annotations'))]
		# )
	# )
	
	# simulated_deconv_filtered <- deconv_reorder(simulated_deconv_filtered)
	
	simulated_deconv_filtered <- simulated_deconv
	
	filtered_genes <- gene_averages_cancer_caf[
		,
		.(pass = ave_exp[cell_type == 'cancer'] > sc_metadata[[ct]]$sim_deconv_filter | ave_exp[cell_type == 'caf'] > sc_metadata[[ct]]$sim_deconv_filter),
		by = gene
	][pass == TRUE, as.character(gene)]
	
	ordered_filtered_genes <- with(simulated_deconv_filtered, genes_filtered[ordering][genes_filtered[ordering] %in% filtered_genes])
	
	simulated_deconv_filtered$ordering <- order(filtered_genes)[order(order(ordered_filtered_genes))]
	simulated_deconv_filtered$genes_filtered <- filtered_genes
	simulated_deconv_filtered$cor_mat <- simulated_deconv_filtered$cor_mat[filtered_genes, filtered_genes]
	simulated_deconv_filtered$cor_with_purity <- sapply(
		simulated_deconv_filtered$cor_with_purity,
		function(x) x[filtered_genes],
		simplify = FALSE,
		USE.NAMES = TRUE
	)
	simulated_deconv_filtered$ccle_comp_diff <- simulated_deconv_filtered$ccle_comp_diff[filtered_genes]
	
	simulated_deconv_filtered_plots <- do.call(
		deconvolve_emt_caf_plots,
		args = c(
			list(
				data = simulated_deconv_filtered,
				heatmap_legend_title = 'Correlation',
				heatmap_colours = rev(colorRampPalette(brewer.pal(11, "RdBu"))(50)),
				heatmap_colour_limits = c(-1, 1),
				heatmap_legend_breaks = c(-1, -0.5, 0, 0.5, 1),
				heatmap_annotations_side = 'left',
				heatmap_annotations_nudge = 0.3,
				purity_colours = rev(colorRampPalette(brewer.pal(11, "PuOr"))(50)),
				purity_colour_limits = c(-1, 1),
				purity_legend_breaks = c(-1, 0, 1),
				purity_legend_title = 'Correlation with purity\n',
				purity_legend_direction = 'horizontal',
				purity_axis_title = NULL,
				ccle_colours = rev(colorRampPalette(brewer.pal(11, 'PiYG')[c(1:3, 6, 9:11)])(50)),
				ccle_fun = function(x) {caTools::runmean(x - mean(x), 30)/max(abs(caTools::runmean(x - mean(x), 30)))},
				ccle_legend_breaks = c(-1, 0, 1),
				ccle_legend_title = 'Tumours vs. cell lines\n',
				ccle_legend_direction = 'horizontal',
				ccle_axis_title = NULL,
				bar_legend_justification = 'left',
				bar_legend_width = unit(10, 'pt'),
				bar_legend_height = unit(10, 'pt')
			),
			sc_metadata[[ct]]$deconv_args[names(sc_metadata[[ct]]$deconv_args) %in% c('plot_title', 'heatmap_annotations')]
		)
	)
	
	plot_data <- sc_deconv_comp_data[gene %in% simulated_deconv_filtered$genes_filtered]
	
	plot_data <- sapply(
		c('cancer', 'caf'),
		function(x) {
			plot_data[
				cell_type == x,
				.(
					id = factor(id, levels = sc_deconv_comp[[x]]$ordered_cell_ids),
					gene = factor(gene, levels = with(simulated_deconv_filtered, genes_filtered[ordering])),
					expression_level = expression_level
				)
			]
		},
		simplify = FALSE,
		USE.NAMES = TRUE
	)
	
	sc_heatmaps <- sapply(
		c('cancer', 'caf'),
		function(x) {
			ggplot(
				plot_data[[x]],
				# plot_data[cell_type == x],
				aes(
					x = gene,
					y = id,
					# x = factor(gene, levels = with(simulated_deconv_filtered, genes_filtered[ordering])),
					# y = factor(id, levels = sc_deconv_comp[[x]]$ordered_cell_ids),
					fill = expression_level
				)
			) +
				geom_raster() +
				scale_x_discrete(expand = c(0, 0)) +
				scale_y_discrete(expand = c(0, 0)) +
				scale_fill_gradientn(
					colours = c(
						colorspace::sequential_hcl(25, h = 190, c = 70, l = c(70, 99), power = 0.8),
						colorspace::sequential_hcl(25, h = 60, c = 100, l = c(80, 99), power = 0.8, rev = TRUE)
					),
					# colours = colorspace::diverge_hcl(50, h = c(190, 70), c = 70, l = c(85, 95), power = 0.75),
					# colours = rev(colorRampPalette(brewer.pal(11, "RdBu"))(50)),
					limits = c(-4, 4),
					oob = scales::squish,
					breaks = c(-4, -2, 0, 2, 4),
					labels = c('-4' = '\u2264 -4', '-2' = '-2', '0' = '0', '2' = '2', '4' = '\u2265 4')
				) +
				theme(
					axis.text = element_blank(),
					axis.title.x = element_blank(),
					axis.ticks = element_blank(),
					axis.ticks.length = unit(0, 'pt'),
					plot.margin = unit(c(5.5, 5.5, 1.25, 0), 'pt')
				) +
				labs(y = mapvalues(x, c('cancer', 'caf'), c('Cancer cells', 'CAFs'), warn_missing = FALSE), fill = 'Relative\nexpression\nlevel')
		},
		simplify = FALSE,
		USE.NAMES = TRUE
	)
	
	
	
	
	
	cairo_pdf(paste0('../data_and_figures/sc_alt_norm/sc_heatmaps_', ct, '.pdf'), width = 8, height = 5, onefile = TRUE)
	cowplot_sc(
		sc_cancer_caf_heatmap_scran,
		legend_space = 0.2,
		heights = c(1.5, 20, 4),
		legend_rel_heights = c(2, 5),
		es_x_axis_title_vjust = 1.3,
		es_y_axis_title_angle = 0,
		es_y_axis_title_xpos = 0.8,
		x_axis_titles_space = 0.3
	) %>% print
	cowplot_sc(
		sc_cancer_caf_heatmap_tpm,
		legend_space = 0.2,
		heights = c(1.5, 20, 4),
		legend_rel_heights = c(2, 5),
		es_x_axis_title_vjust = 1.3,
		es_y_axis_title_angle = 0,
		es_y_axis_title_xpos = 0.8,
		x_axis_titles_space = 0.3
	) %>% print
	dev.off()
	
	pdf(paste0('../data_and_figures/sc_alt_norm/sim_lineplots_', ct, '.pdf'), width = 6, height = 4)
	plot_grid(lineplot_scran$lineplot + theme(legend.position = 'none'), dummy_legend, nrow = 1, ncol = 2, rel_widths = c(2, 1)) %>% print
	plot_grid(lineplot_tpm$lineplot + theme(legend.position = 'none'), dummy_legend, nrow = 1, ncol = 2, rel_widths = c(2, 1)) %>% print
	dev.off()
	
	cairo_pdf(paste0('../data_and_figures/sc_alt_norm/sim_deconv_', ct, '.pdf'), width = 7, height = 6)
	deconv_plot(
		list(simulated_deconv_plots),
		legends_rel_size = c(0.5, 0.75, 0, 0.75, 0, 0.75, 0.75, 3, 2),
		legends_space = 0.7,
		title.position = 'right'
	) %>% print
	dev.off()
	
	sc_sim_deconv_comp_figures <- sapply(
		c(simulated_deconv_filtered_plots$plots, sc_heatmaps),
		function(x) x + theme(legend.position = 'none', plot.title = element_blank(), axis.title.y = element_blank()),
		simplify = FALSE,
		USE.NAMES = TRUE
	)
	sc_sim_deconv_comp_figures$heatmap <- sc_sim_deconv_comp_figures$heatmap + theme(plot.margin = unit(c(0, 5.5, 1.25, 0), 'pt'))
	sc_sim_deconv_comp_figures$axis_labels <- sc_sim_deconv_comp_figures$axis_labels + theme(plot.margin = unit(c(0, 0, 1.25, 0), 'pt'))
	
	cairo_pdf(paste0('../data_and_figures/sc_alt_norm/sc_sim_deconv_comp_', ct, '.pdf'), width = 6, height = 7)
	plot_grid(
		plot_grid(
			blank_plot(),
			pemt_caf_brackets(
				edge = 'right',
				pemt_bracket_max = sc_metadata[[ct]]$pemt_bracket_max,
				caf_bracket_min = sc_metadata[[ct]]$caf_bracket_min,
				brackets_just = 0.65,
				labels_just = -1,
				labels = c('pEMT genes', 'CAF genes'),
				labels_angle = 90
			),
			blank_plot(),
			nrow = 3,
			ncol = 1,
			rel_heights = c(3, 15, 10)
		),
		plot_grid(
			blank_plot(),
			sc_sim_deconv_comp_figures$axis_labels,
			blank_plot() +
				labs(y = 'Cancer cells') +
				scale_y_continuous(position = 'right') +
				theme(plot.margin = unit(c(5.5, 5.5, 1.25, 5.5), 'pt'), axis.title.y.right = element_text(angle = 90)),
			blank_plot() +
				labs(y = 'CAFs') +
				scale_y_continuous(position = 'right') +
				theme(plot.margin = unit(c(5.5, 5.5, 1.25, 5.5), 'pt'), axis.title.y.right = element_text(angle = 90)),
			nrow = 4,
			ncol = 1,
			rel_heights = c(3, 15, 5, 5)
		),
		plot_grid(
			plotlist = c(
				list(blank_plot() + theme(plot.margin = unit(c(11, 5.5, 5.5, 0), 'pt')) + labs(title = sc_metadata[[ct]]$deconv_args$plot_title)),
				sc_sim_deconv_comp_figures[c('purity_bar', 'ccle_bar', 'heatmap', 'cancer', 'caf')]
				# lapply(
					# c(simulated_deconv_filtered_plots$plots[c('purity_bar', 'ccle_bar', 'heatmap')], sc_heatmaps),
					# function(x) x + theme(legend.position = 'none', plot.title = element_blank())
				# )
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
						simulated_deconv_filtered_plots$plots$purity_bar +
							theme(legend.justification = c(0, 1), legend.direction = 'vertical') +#, legend.key.width = NULL, legend.key.height = NULL) +
							labs(fill = 'Correlation\nwith purity')
					),
					get_legend(
						simulated_deconv_filtered_plots$plots$ccle_bar +
							theme(legend.justification = c(0, 1), legend.direction = 'vertical') +#, legend.key.width = NULL, legend.key.height = NULL) +
							labs(fill = 'Tumours vs.\ncell lines')
					),
					get_legend(simulated_deconv_filtered_plots$plots$heatmap + theme(legend.justification = c(0, 1))),
					get_legend(sc_heatmaps[[1]] + theme(legend.justification = c(0, 1)))
				)
			),
			nrow = 5,
			ncol = 1,
			rel_heights = c(1, 5, 5, 7, 10)
		),
		nrow = 1,
		ncol = 4,
		rel_widths = c(0.5, 1, 5, 1.5)
	) %>% print
	dev.off()
	
	all_figures[[ct]] <- list(
		sc_cancer_caf_heatmap_scran = sc_cancer_caf_heatmap_scran,
		sc_cancer_caf_heatmap_tpm = sc_cancer_caf_heatmap_tpm,
		lineplot_scran = lineplot_scran$lineplot,
		lineplot_tpm = lineplot_tpm$lineplot,
		simulated_deconv_plots = simulated_deconv_plots$plots,
		simulated_deconv_filtered_plots = simulated_deconv_filtered_plots$plots,
		sc_heatmaps = sc_heatmaps
	)
	
	if('annotations_side' %in% names(sc_metadata[[ct]])) {
		all_figures[[ct]] <- c(all_figures[[ct]], list(sc_cancer_caf_heatmap_combining = sc_cancer_caf_heatmap_combining))
	}
	
}





# Final figures:

dummy_legend_plot <- ggplot(
	data = data.table(
		x = 1:7,
		y = 1,
		f = factor(
			c('b_cell', 'cancer', 'endothelial', 'caf', 'macrophage', 'mast', 't_cell'),
			levels = c('b_cell', 'cancer', 'endothelial', 'caf', 'macrophage', 'mast', 't_cell')
		)
	)
) + 
	geom_tile(aes(x = x, y = y, fill = f)) +
	scale_fill_manual(
		labels = c(
			'b_cell' = 'B cell',
			'cancer' = 'Cancer',
			'endothelial' = 'Endothelial',
			'caf' = 'CAF',
			'macrophage' = 'Macrophage',
			'mast' = 'Mast',
			't_cell' = 'T cell'
		),
		values = c(
			'b_cell' = '#8DD3C7',
			'cancer' = '#FB8072',
			'endothelial' = '#BC80BD',
			'caf' = '#FDB462',
			'macrophage' = '#80B1D3',
			'mast' = '#FCCDE5',
			't_cell' = '#B3DE69'
		)
	) +
	labs(fill = 'Cell type') +
	theme(legend.key = element_rect(size = 1, colour = 'white'), legend.key.size = unit(15, 'pt'))

dummy_legend <- get_legend(dummy_legend_plot)

cairo_pdf(
    '../data_and_figures/final_figures_resubmission/S7.pdf',
    width = 12,
    height = 15.5
)

# It looks like there's a minor problem with the following: when getting legends for the last row of the plot, I'm using a variable "ct"
# outside the lapply() loop.  The fact that it works is presumably because when we source this whole script, the last value of "ct" from
# the above for loop will remain in the environment, and is used in this getting legends step.  We should replace this "ct" with any one
# of 'coadread', 'luad' or 'paad' (since the legends should be the same for all three).  I haven't done this yet because I didn't want
# to spend the time checking for sure that it worked, but it really should work.

print(
	plot_grid(
		blank_plot(),
		plot_grid(
			plot_grid(
				cowplot_sc(
					all_figures$coadread$sc_cancer_caf_heatmap_scran,
					legend_space = 0,
					heights = c(4, 20, 6),
					es_x_axis_title_vjust = 1.3,
					es_y_axis_title_angle = 0,
					es_y_axis_title_xpos = 0.8
				),
				# blank_plot(),
				cowplot_sc(
					all_figures$luad$sc_cancer_caf_heatmap_scran,
					legend_space = 0,
					heights = c(4, 20, 6),
					es_x_axis_title_vjust = 1.3,
					es_y_axis_title_angle = 0,
					es_y_axis_title_xpos = 0.8
				),
				cowplot_sc(
					all_figures$paad$sc_cancer_caf_heatmap_scran,
					legend_space = 0,
					heights = c(4, 20, 6),
					es_x_axis_title_vjust = 1.3,
					es_y_axis_title_angle = 0,
					es_y_axis_title_xpos = 0.8
				),
				ncol = 3,
				nrow = 1
				# rel_widths = c(20, 1, 20)
			),
			# blank_plot(),
			# plot_grid(
				# cowplot_sc(
					# all_figures$paad$sc_cancer_caf_heatmap_combining,
					# legend_space = 0,
					# heights = c(1.5, 20, 4),
					# es_x_axis_title_vjust = 1.3,
					# es_y_axis_title_angle = 0,
					# es_y_axis_title_xpos = 0.8
				# ),
				blank_plot(),
				plot_grid(
					all_figures$paad$sc_cancer_caf_heatmap_scran$plots$genes_detected_legend,
					blank_plot(),
					all_figures$paad$sc_cancer_caf_heatmap_scran$plots$heatmap_legend,
					nrow = 1,
					ncol = 3,
					rel_widths = c(5, 1, 5)
				),
				# ncol = 3,
				# nrow = 1,
				# rel_widths = c(20, 1, 20)
			# ),
			ncol = 1,
			nrow = 3,
			rel_heights = c(20, 2, 3)
		),
		blank_plot(),
		plot_grid(
			plot_grid(
				plotlist = list(
					all_figures$coadread$lineplot_scran +
						theme(legend.position = 'none', plot.margin = unit(c(5.5, 5.5, 15, 30), 'pt')),# +
						# labs(x = NULL, y = NULL),
					all_figures$luad$lineplot_scran +
						theme(legend.position = 'none', plot.margin = unit(c(5.5, 5.5, 15, 15), 'pt')) +
						labs(y = NULL),
					all_figures$paad$lineplot_scran +
						theme(legend.position = 'none', plot.margin = unit(c(5.5, 5.5, 15, 15), 'pt')) +
						labs(y = NULL)
				),
				nrow = 1,
				ncol = 3,
				rel_widths = c(1.125, 1, 1, 1), # Left plots are squashed by larger left margin
				align = 'h'
			),
			dummy_legend,
			nrow = 1,
			ncol = 2,
			rel_widths = c(7.5, 1.5)
		),# +
			# draw_label('Proportion of tumour', x = 0.47, y = 0, vjust = -0.5, size = 12) +
			# draw_label('Proportion of gene expression', x = 0, y = 0.5, vjust = 1.3, angle = 90, size = 12),
		blank_plot(),
		plot_grid(
			plot_grid(
				plotlist = lapply(
					c('coadread', 'luad', 'paad'),
					function(ct) {
						
						sc_sim_deconv_comp_figures <- sapply(
							c(all_figures[[ct]]$simulated_deconv_filtered_plots, all_figures[[ct]]$sc_heatmaps),
							function(x) x + theme(legend.position = 'none', plot.title = element_blank(), axis.title.y = element_blank()),
							simplify = FALSE,
							USE.NAMES = TRUE
						)
						sc_sim_deconv_comp_figures$heatmap <- sc_sim_deconv_comp_figures$heatmap +
							theme(plot.margin = unit(c(0, 5.5, 1.25, 0), 'pt'))
						sc_sim_deconv_comp_figures$axis_labels <- sc_sim_deconv_comp_figures$axis_labels +
							theme(plot.margin = unit(c(0, 0, 1.25, 0), 'pt'))
						
						plot_grid(
							plot_grid(
								blank_plot(),
								sc_sim_deconv_comp_figures$axis_labels,
								blank_plot() +
									labs(y = 'Cancer\ncells') +
									scale_y_continuous(position = 'right') +
									theme(plot.margin = unit(c(5.5, 5.5, 1.25, 5.5), 'pt'), axis.title.y.right = element_text(angle = 90)),
								blank_plot() +
									labs(y = 'CAFs') +
									scale_y_continuous(position = 'right') +
									theme(plot.margin = unit(c(5.5, 5.5, 1.25, 5.5), 'pt'), axis.title.y.right = element_text(angle = 90)),
								nrow = 4,
								ncol = 1,
								rel_heights = c(4, 15, 5, 5)
							),
							plot_grid(
								plotlist = c(
									list(
										blank_plot() +
										theme(plot.margin = unit(c(11, 5.5, 5.5, 0), 'pt')) +
										labs(title = sc_metadata[[ct]]$deconv_args$plot_title)
									),
									sc_sim_deconv_comp_figures[c('purity_bar', 'ccle_bar', 'heatmap', 'cancer', 'caf')]
								),
								nrow = 6,
								ncol = 1,
								align = 'v',
								rel_heights = c(2, 1, 1, 15, 5, 5)
							),
							nrow = 1,
							ncol = 2,
							rel_widths = c(1.6, 5)
						)
						
					}
				),
				nrow = 1,
				ncol = 3
			),
			blank_plot(),
			plot_grid(
				get_legend(
					all_figures[[ct]]$simulated_deconv_filtered_plots$purity_bar +
						theme(legend.justification = c(1, 1), legend.direction = 'horizontal', legend.key.width = NULL) +
						labs(fill = 'Correlation with purity\n')
				),
				blank_plot(),
				get_legend(
					all_figures[[ct]]$simulated_deconv_filtered_plots$heatmap +
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
					all_figures[[ct]]$simulated_deconv_filtered_plots$ccle_bar +
						theme(legend.justification = c(1, 1), legend.direction = 'horizontal', legend.key.width = NULL) +
						labs(fill = 'Tumours vs. cell lines\n')
				),
				blank_plot(),
				get_legend(
					all_figures[[ct]]$sc_heatmaps[[1]] +
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
		),
		nrow = 6,
		ncol = 1,
		rel_heights =c(0.15, 1.2, 0.125, 0.95, 0.125, 1.8)
	) +
		draw_label('A', x = 0, y = 0.99, hjust = -0.1, vjust = 1.1, size = 20, fontface = 2) +
		draw_label('B', x = 0, y = 0.68, hjust = -0.1, vjust = 1.1, size = 20, fontface = 2) +
		draw_label('C', x = 0, y = 0.43, hjust = -0.1, vjust = 1.1, size = 20, fontface = 2)
)

dev.off()
