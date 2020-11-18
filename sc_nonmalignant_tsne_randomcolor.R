library(data.table) # 1.12.8
library(ggplot2) # This is version 3.3.1 on my laptop's WSL setup but 3.3.0 on WEXAC.
library(magrittr) # 1.5

source('general_functions.R')

cell_type_markers <- fread('../../cell_type_markers.csv')





cohort_data <- list(
	breast_qian = list(
		patients = c(42, 43, 47, 49, 51, 53, 54),
		read_quote = quote(fread('../data_and_figures/qian_breast_2020.csv'))
	),
	crc_lee_smc = list(
		patients = c('SMC01', 'SMC02', 'SMC04', 'SMC07', 'SMC08', 'SMC09', 'SMC11', 'SMC14', 'SMC15', 'SMC16', 'SMC18', 'SMC20', 'SMC21', 'SMC23', 'SMC25'),
		read_quote = quote(fread('../data_and_figures/lee_crc_2020_smc.csv')[
			sample_type != 'normal',
			-c('sample_id', 'sample_type', 'cell_subtype')
		])
	),
	hnscc_puram = list(
		patients = c(5, 6, 18, 20, 22, 25, 26),
		read_quote = quote(fread('../data_and_figures/puram_hnscc_2017.csv')[, -c('lymph_node', 'processed_by_maxima_enzyme')])
	),
	liver_ma = list(
		patients = c('C25', 'C26', 'C46', 'C56', 'C66', 'H37', 'H38', 'H65'),
		read_quote = quote(fread('../data_and_figures/ma_liver_2019.csv')[, -c('sample', 'disease')])
	),
	luad_kim = list(
		patients = c('P0006', 'P0008', 'P0018', 'P0019', 'P0020', 'P0025', 'P0028', 'P0030', 'P0031', 'P0034', 'P1006', 'P1028', 'P1049', 'P1058'),
		read_quote = quote(fread('../data_and_figures/kim_luad_2020.csv')[, -c('cell_type_refined', 'cell_subtype', 'sample_site', 'sample_id', 'sample_origin')])
	),
	lung_qian = list(
		patients = 1:8,
		read_quote = quote(fread('../data_and_figures/qian_lung_2020.csv')[sample_type != 'normal', -c('sample_type', 'disease')])
	),
	ovarian_qian = list(
		patients = 11:14,
		read_quote = quote(fread('../data_and_figures/qian_ovarian_2020.csv')[sample_type != 'normal', -c('sample_type', 'site')])
	),
	pdac_peng = list(
		patients = c('T2', 'T3', 'T6', 'T7', 'T8', 'T9', 'T11', 'T13', 'T14', 'T15', 'T16', 'T17', 'T19', 'T20', 'T21', 'T22'),
		read_quote = quote(fread('../data_and_figures/peng_pdac_2019.csv')[sample_type != 'normal', -'sample_type'])
	)
)





for(cohort in names(cohort_data)) {
	
	cat(paste0(cohort, '\n'))
	
	classification_data <- fread(paste0('../data_and_figures/sc_reclassify/classification_data_', cohort, '.csv'), key = 'cell_id')
	
	sc_data <- eval(cohort_data[[cohort]]$read_quote)[patient %in% cohort_data[[cohort]]$patients]
	
	sc_data[, classification := classification_data[id, classification]]
	
	sc_data <- sc_data[classification == 'nonmalignant', -'classification']
	
	gene_averages <- sapply(
		sc_data[, -c('id', 'patient', 'cell_type')],
		function(x) {log2(mean(10*(2^x - 1)) + 1)},
		USE.NAMES = TRUE
	)
	
	sc_data <- sc_data[, c('id', 'patient', 'cell_type', names(gene_averages)[gene_averages >= 4]), with = FALSE]
	
	sc_tsne <- readRDS(paste0('../data_and_figures/sc_reclassify/tsne_', cohort, '.rds'))
	
	sc_dbscan <- readRDS(paste0('../data_and_figures/sc_reclassify/dbscan_', cohort, '.rds'))
	
	plot_data <- setNames(
		cbind(as.data.table(sc_tsne$Y), as.character(sc_dbscan$cluster), sc_data$cell_type),
		c('x', 'y', 'dbscan_cluster', 'cell_type')
	)
	
	plot_data[dbscan_cluster == 0, dbscan_cluster := 'noise']
	
	# Scatter plot of t-SNE coordinates, coloured by DBSCAN cluster:
	tsne_plot_dbscan <- ggplot(plot_data, aes(x = x, y = y, colour = dbscan_cluster)) +
		geom_point() +
		scale_colour_manual(values = randomcoloR::distinctColorPalette(length(unique(sc_dbscan$cluster)))) +
		theme_minimal()
	
	tsne_plot_author_cell_types <- ggplot(plot_data, aes(x = x, y = y, colour = cell_type)) +
		geom_point() +
		scale_colour_manual(values = randomcoloR::distinctColorPalette(length(unique(sc_data$cell_type)))) +
		theme_minimal()
	
	tsne_plot_complexity <- ggplot(
		cbind(plot_data, complexity = apply(sc_data[, -c('id', 'patient', 'cell_type')], 1, function(x) sum(x > 0))),
		aes(x = x, y = y, colour = complexity)
	) +
		geom_point() +
		theme_minimal() +
		labs(title = 'Complexity')
	
	ct_ave_exp_plots <- sapply(
		cell_type_markers[cell_type !='mesenchymal', unique(cell_type)],
		function(ct) {
			
			ave_exp <- sc_data[, rowMeans(.SD), .SDcols = cell_type_markers[cell_type == ct & gene %in% names(sc_data), unique(gene)]]
			
			tsne_plot <- ggplot(cbind(plot_data, ave_exp = ave_exp), aes(x = x, y = y, colour = ave_exp)) +
				geom_point() +
				theme_minimal() +
				labs(title = ct)
			
			list(ave_exp = ave_exp, plot = tsne_plot)
			
		},
		simplify = FALSE,
		USE.NAMES = TRUE
	)
	
	epithelial_markers <- c(
		'CDH1',
		'EPCAM',
		'SFN',
		names(sc_data)[grepl('^KRT[0-9]', names(sc_data))]
	)
	
	epithelial_markers <- epithelial_markers[epithelial_markers %in% names(sc_data)]
	
	ave_exp_epithelial <- sc_data[, rowMeans(.SD), .SDcols = epithelial_markers]
	
	tsne_plot_epithelial <- ggplot(cbind(plot_data, ave_exp = ave_exp_epithelial), aes(x = x, y = y, colour = ave_exp)) +
		geom_point() +
		theme_minimal() +
		labs(title = 'epithelial')
	
	pdf(paste0('../data_and_figures/sc_reclassify/tsne_plots_normal_', cohort, '.pdf'))
	
	print(tsne_plot_dbscan)
	print(tsne_plot_author_cell_types)
	print(tsne_plot_complexity)
	print(tsne_plot_epithelial)
	for(ct_list in ct_ave_exp_plots) print(ct_list$plot)
	
	dev.off()
	
}
