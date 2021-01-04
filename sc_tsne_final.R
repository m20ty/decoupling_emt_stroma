library(data.table) # 1.12.8
library(ggplot2) # 3.3.0
library(magrittr) # 1.5
library(cowplot) # 1.0.0
library(Rtsne) # 0.15
library(plyr) # 1.8.6
library(stringr) # 1.4.0
library(colorspace) # 1.4.1

source('general_functions.R')

cell_type_markers <- fread('../../cell_type_markers.csv')[!(cell_type %in% c('mesenchymal', 'myocyte')) & source != 'TCGA_CCLE_comparison']





cohort_data <- list(
	breast_qian = list(
		read_quote = quote(fread('../data_and_figures/qian_breast_2020_reclassified.csv')),
		seeds = c(9699, 6592, 9218),
		cells_to_exclude = 'sc5rJUQ064_CCATGTCCATCCCATC'
	),
	crc_lee_smc = list(
		read_quote = quote(fread('../data_and_figures/lee_crc_2020_smc_reclassified.csv')),
		seeds = c(9128, 1076, 5329)
	),
    hnscc_puram = list(
		read_quote = quote(fread('../data_and_figures/puram_hnscc_2017_reclassified.csv')),
		seeds = c(8889, 9887, 9134)
	),
	liver_ma = list(
		read_quote = quote(fread('../data_and_figures/ma_liver_2019_reclassified.csv')),
		seeds = c(3081, 46, 7526)
	),
	luad_kim = list(
		read_quote = quote(fread('../data_and_figures/kim_luad_2020_reclassified.csv')),
		seeds = c(4573, 2106, 7294)
	),
	lung_qian = list(
		read_quote = quote(fread('../data_and_figures/qian_lung_2020_reclassified.csv')[, -'disease']),
		seeds = c(9201, 4923, 7716)
	),
	ovarian_qian = list(
		read_quote = quote(fread('../data_and_figures/qian_ovarian_2020_reclassified.csv')),
		seeds = c(9419, 3874, 2197),
		cells_to_exclude = c('scrSOL001_TCATTTGTCTGTCAAG', 'scrSOL004_TTGCCGTTCTCCTATA')
	),
    pdac_peng = list(
		read_quote = quote(fread('../data_and_figures/peng_pdac_2019_reclassified.csv')[startsWith(cell_type, 'ductal'), cell_type := 'ductal']),
		seeds = c(8374, 3685, 273),
		cells_to_exclude = c('T8_TGGTTCCTCGCATGGC', 'T17_CGTGTAACAGTACACT')
	)
)





cell_type_labels <- c(
	'acinar' = 'Acinar',
	'alveolar' = 'Alveolar',
	'b_cell' = 'B cell',
	'caf' = 'CAF',
	'caf_potential' = 'Potential CAF',
	'cancer' = 'Cancer',
	'dendritic' = 'Dendritic',
	'ductal' = 'Ductal',
	'endocrine' = 'Endocrine',
	'endothelial' = 'Endothelial',
	'epithelial' = 'Epithelial',
	'erythroblast' = 'Erythroblast',
	'hpc-like' = 'HPC-like',
	'macrophage' = 'Macrophage',
	'mast' = 'Mast',
	'myocyte' = 'Myocyte',
	'nk_cell' = 'NK cell',
	'stellate' = 'Stellate',
	't_cell' = 'T cell'
)

cell_type_colours <- c(
	'acinar' = '#C0DF84',
	'alveolar' = '#E3E050',
	'b_cell' = '#887DDA',
	'caf' = '#D8A354',
	'caf_potential' = lighten('#D8A354', 0.5),
	'cancer' = '#DD4C61',
	'dendritic' = '#7AE1D8',
	'ductal' = '#DBE2BB',
	'endocrine' = '#7AE759',
	'endothelial' = '#8C5581',
	'epithelial' = '#D1D2DC',
	'erythroblast' = '#DCA9DB',
	'hpc-like' = '#A83EE5',
	'macrophage' = '#75B0DC',
	'mast' = '#799575',
	'myocyte' = '#DD6BCF',
	'nk_cell' = '#01AD1E',
	'stellate' = '#D5978A',
	't_cell' = '#79E5A3'
)





all_plot_data <- list(
	breast_qian = list(),
	crc_lee_smc = list(),
	hnscc_puram = list(),
	liver_ma = list(),
	luad_kim = list(),
	lung_qian = list(),
	ovarian_qian = list(),
	pdac_peng = list()
)

for(cohort in names(cohort_data)) {

	cat(paste0(cohort, '\n'))

	sc_data <- eval(cohort_data[[cohort]]$read_quote)[, -'cell_type_author']

	if('cell_type_lenient' %in% names(sc_data)) {

		sc_data <- sc_data[cell_type_lenient != 'ambiguous']

		sc_data[, cell_type_for_plot := cell_type]

		# If potential CAFs are marked as endothelial, change to 'caf_potential':
		if(sc_data[cell_type == 'endothelial' & cell_type_lenient == 'caf', .N] > 0) {
			sc_data[cell_type == 'endothelial' & cell_type_lenient == 'caf', cell_type_for_plot := 'caf_potential']
		}

		gene_averages <- sapply(
			sc_data[, -c('id', 'patient', 'cell_type', 'cell_type_lenient', 'cell_type_for_plot')],
			function(x) {log2(mean(10*(2^x - 1)) + 1)},
			USE.NAMES = TRUE
		)

		sc_data <- sc_data[
			,
			c('id', 'patient', 'cell_type', 'cell_type_lenient', 'cell_type_for_plot', names(gene_averages)[gene_averages >= 4]),
			with = FALSE
		]

		# Run t-SNE:
		# set.seed(cohort_data[[cohort]]$seeds[1])
		# sc_tsne <- Rtsne(
		# 	as.matrix(
		# 		sc_data[
		# 			,
		# 			lapply(.SD, function(x) {x - mean(x)}), .SDcols = -c('id', 'patient', 'cell_type', 'cell_type_lenient', 'cell_type_for_plot')
		# 		]
		# 	)#,
		# 	# num_threads = 16
		# )

		# saveRDS(sc_tsne, paste0('../data_and_figures/sc_tsne_final/tsne_', cohort, '.rds'))

		sc_tsne <- readRDS(paste0('../data_and_figures/sc_tsne_final/tsne_', cohort, '.rds'))

		if('cells_to_exclude' %in% names(cohort_data[[cohort]])) {
			sc_tsne$Y <- sc_tsne$Y[-which(sc_data$id %in% cohort_data[[cohort]]$cells_to_exclude), ]
			sc_data <- sc_data[!(id %in% cohort_data[[cohort]]$cells_to_exclude)]
		}

		plot_data <- setNames(
			cbind(as.data.table(sc_tsne$Y), sc_data[, .(patient, cell_type, cell_type_lenient, cell_type_for_plot)]),
			c('x', 'y', 'patient', 'cell_type', 'cell_type_lenient', 'cell_type_for_plot')
		)

		plot_data[cell_type == 'ambiguous', c('cell_type', 'cell_type_for_plot') := .('caf_potential', 'caf_potential')]

		tsne_plot_patient <- ggplot(plot_data, aes(x = x, y = y, colour = as.character(patient))) +
			geom_point() +
			theme_minimal() +
			labs(title = 'Patients', colour = NULL, x = 't-SNE 1', y = 't-SNE 2')

		tsne_plot_cell_types <- ggplot(plot_data, aes(x = x, y = y, colour = cell_type_for_plot)) +
			geom_point() +
			scale_colour_manual(
				labels = cell_type_labels[sort(unique(plot_data$cell_type_for_plot))],
				values = cell_type_colours[sort(unique(plot_data$cell_type_for_plot))]
			) +
			theme_minimal() +
			labs(title = 'Cell types', colour = NULL, x = 't-SNE 1', y = 't-SNE 2')

		tsne_plot_cell_types_lenient <- ggplot(plot_data, aes(x = x, y = y, colour = cell_type_lenient)) +
			geom_point() +
			scale_colour_manual(
				labels = cell_type_labels[sort(unique(plot_data$cell_type_lenient))],
				values = cell_type_colours[sort(unique(plot_data$cell_type_lenient))]
			) +
			theme_minimal() +
			labs(title = 'Cell types (lenient)', colour = NULL, x = 't-SNE 1', y = 't-SNE 2')

		# CNA scores:

		cna_score_files <- dir(paste0('../data_and_figures/sc_find_malignant/', cohort))
		cna_score_files <- cna_score_files[endsWith(cna_score_files, '_data.csv')]

		cna_score_data <- lapply(
			cna_score_files,
			function(x) cbind(fread(paste0('../data_and_figures/sc_find_malignant/', cohort, '/', x)))#, patient = gsub('_data.csv', '', x))
		) %>% rbindlist(fill = TRUE)

		setkey(cna_score_data, cell_id)

		cna_score_data <- cna_score_data[sc_data$id]

		cna_score_data[
			,
			classification_final := switch(
				(length(classification_final) == 1) + 1,
				switch( # In this case the cell must have been reference, but could have been classified as ambiguous or malignant in some patients.
					('ambiguous' %in% classification_final | 'malignant' %in% classification_final) + 1,
					'nonmalignant',
					'ambiguous'
				),
				classification_final # In this case the cell was not reference.
			),
			by = cell_id
		]

		cna_score_data <- cna_score_data[
			,
			.(
				classification_final = unique(classification_final),
				cna_score = switch(
					startsWith(unique(classification_final), 'malignant') + 1,
					mean(cna_signal),
					.SD[classification == classification_final, cna_signal]
				)
			),
			by = cell_id
		]

		cna_score_data[, patient := sc_data$patient]
		cna_score_data[, cna_score := cna_score/quantile(cna_score[startsWith(classification_final, 'malignant')], 0.75), by = patient]

		# cell_id should still be in the same order as sc_data$id, so it's safe to do this:
		plot_data[, cna_score := cna_score_data$cna_score]

		tsne_plot_cna_score <- ggplot(plot_data, aes(x = x, y = y, colour = cna_score)) +
			geom_point() +
			theme_minimal() +
			scale_colour_gradientn(colours = colorRamps::matlab.like(50), limits = c(0, 1), oob = scales::squish) +
			labs(title = 'CNA score', colour = NULL, x = 't-SNE 1', y = 't-SNE 2')

		immune_cell_types <- c('B', 'B_plasma', 'DC', 'macrophage', 'mast', 'T')[
			sapply(
				c('B', 'B_plasma', 'DC', 'macrophage', 'mast', 'T'),
				function(ct) cell_type_markers[cell_type == ct, sum(gene %in% names(sc_data)) >= 3] # This was previously just > 0
			)
		]

		immune_cell_types <- immune_cell_types[
			mapvalues(
				immune_cell_types,
				c('B', 'B_plasma', 'DC', 'macrophage', 'mast', 'T'),
				c('b_cell', 'b_cell', 'dendritic', 'macrophage', 'mast', 't_cell'),
				warn_missing = FALSE
			) %in% sc_data$cell_type
		]

		set.seed(cohort_data[[cohort]]$seeds[2])

		immune_scores <- as.data.table(
			sapply(
				immune_cell_types,
				function(ct) {
					signature_score(
						sc_data[, set_colnames(t(.SD), id), .SDcols = -c('id', 'patient', 'cell_type', 'cell_type_lenient', 'cell_type_for_plot')],
						sig_genes = cell_type_markers[cell_type == ct & gene %in% names(sc_data), unique(gene)],
						nbin = 30,
						n = 100
					)
				}
			),
			keep.rownames = 'cell_id'
		)

		immune_sigs <- sapply(
			immune_cell_types,
			function(ct) {

				sig_cor <- cell_type_markers[
					cell_type == ct & gene %in% names(sc_data),
					cor(sc_data[, unique(gene), with = FALSE], immune_scores[, ..ct])[, 1]
				]

				if(sum(sig_cor > 0.6) < 10) {
					if(length(sig_cor) <= 10) {
						return(names(sig_cor))
					} else {
						return(names(sig_cor)[names(sig_cor) %in% names(sort(sig_cor, decreasing = TRUE)[1:10])])
					}
				} else {
					return(names(sig_cor)[sig_cor > 0.6])
				}

			}
		)

		immune_scores <- as.data.table(
			sapply(
				immune_cell_types,
				function(ct) {
					signature_score(
						sc_data[, set_colnames(t(.SD), id), .SDcols = -c('id', 'patient', 'cell_type', 'cell_type_lenient', 'cell_type_for_plot')],
						sig_genes = immune_sigs[[ct]],
						nbin = 30,
						n = 100
					)
				}
			),
			keep.rownames = 'cell_id'
		)

		immune_scores[
			,
			c('which_max', 'cell_type') := .(apply(.SD, 1, function(x) names(.SD)[which.max(x)]), sc_data$cell_type),
			.SDcols = -'cell_id'
		]

		immune_cell_type_scaling_factors <- sapply(
			immune_cell_types,
			function(ict) {
				immune_scores[
					cell_type == mapvalues(
						ict,
						c('B', 'B_plasma', 'DC', 'macrophage', 'mast', 'T'),
						c('b_cell', 'b_cell', 'dendritic', 'macrophage', 'mast', 't_cell'),
						warn_missing = FALSE
					),
					setNames(quantile(get(ict), 0.9), NULL)
				]
			},
			USE.NAMES = TRUE
		)

		immune_scores[, score := get(unique(which_max))/immune_cell_type_scaling_factors[unique(which_max)], by = which_max]

		set.seed(cohort_data[[cohort]]$seeds[3])
		plot_data[
			,
			c('immune_score', 'caf_score', 'endothelial_score') := .( # Should add pericyte score, but need to make my own signature
				immune_scores$score,
				signature_score(
					sc_data[, set_colnames(t(.SD), id), .SDcols = -c('id', 'patient', 'cell_type', 'cell_type_lenient', 'cell_type_for_plot')],
					sig_genes = cell_type_markers[cell_type == 'CAF' & gene %in% names(sc_data), unique(gene)],
					nbin = 30,
					n = 100
				),
				signature_score(
					sc_data[, set_colnames(t(.SD), id), .SDcols = -c('id', 'patient', 'cell_type', 'cell_type_lenient', 'cell_type_for_plot')],
					sig_genes = cell_type_markers[cell_type == 'endothelial' & gene %in% names(sc_data), unique(gene)],
					nbin = 30,
					n = 100
				)
			)
		]

		caf_endothelial_sigs <- sapply(
			c('caf', 'endothelial'),
			function(ct) {

				sig_cor <- cell_type_markers[
					cell_type == mapvalues(ct, 'caf', 'CAF', warn_missing = FALSE) & gene %in% names(sc_data),
					cor(sc_data[, unique(gene), with = FALSE], plot_data[, get(paste0(ct, '_score'))])[, 1]
				]

				if(sum(sig_cor > 0.6) < 10) {
					if(length(sig_cor) <= 10) {
						return(names(sig_cor))
					} else {
						return(names(sig_cor)[names(sig_cor) %in% names(sort(sig_cor, decreasing = TRUE)[1:10])])
					}
				} else {
					return(names(sig_cor)[sig_cor > 0.6])
				}

			}
		)

		plot_data[
			,
			c('caf_score', 'endothelial_score') := .(
				signature_score(
					sc_data[, set_colnames(t(.SD), id), .SDcols = -c('id', 'patient', 'cell_type', 'cell_type_lenient', 'cell_type_for_plot')],
					sig_genes = caf_endothelial_sigs$caf,
					nbin = 30,
					n = 100
				),
				signature_score(
					sc_data[, set_colnames(t(.SD), id), .SDcols = -c('id', 'patient', 'cell_type', 'cell_type_lenient', 'cell_type_for_plot')],
					sig_genes = caf_endothelial_sigs$endothelial,
					nbin = 30,
					n = 100
				)
			)
		]

		plot_data[
			,
			c('caf_score', 'endothelial_score') := .(
				caf_score/.SD[cell_type == 'caf', quantile(caf_score, 0.9)],
				endothelial_score/.SD[cell_type == 'endothelial', quantile(endothelial_score, 0.9)]
			)
		]

		score_tsnes <- sapply(
			c('immune_score', 'caf_score', 'endothelial_score'),
			function(score_type) {
				ggplot(plot_data, aes(x = x, y = y, colour = get(score_type))) +
					geom_point() +
					theme_minimal() +
					scale_colour_gradientn(colours = colorRamps::matlab.like(50), limits = c(0, 1), oob = scales::squish) +
					labs(
						title = mapvalues(
							score_type,
							c('immune_score', 'caf_score', 'endothelial_score'),
							c('Immune score', 'CAF score', 'Endothelial score'),
							warn_missing = FALSE
						),
						colour = NULL,
						x = 't-SNE 1',
						y = 't-SNE 2'
					)
			},
			simplify = FALSE,
			USE.NAMES = TRUE
		)

		pdf(paste0('../data_and_figures/sc_tsne_final/tsne_plots_all_', cohort, '.pdf'))
		print(tsne_plot_cell_types)
		print(tsne_plot_cell_types_lenient)
		print(tsne_plot_patient)
		for(p in score_tsnes) print(p)
		print(tsne_plot_cna_score)
		dev.off()

		pdf(paste0('../data_and_figures/sc_tsne_final/tsne_plots_combined_', cohort, '.pdf'), width = 15, height = 5)
		plot_grid(
			tsne_plot_cell_types + theme(legend.position = 'none'),
			get_legend(tsne_plot_cell_types),
			score_tsnes$caf_score + theme(legend.position = 'none', axis.title.y = element_blank()),
			tsne_plot_cna_score + theme(legend.position = 'none', axis.title.y = element_blank()),
			get_legend(tsne_plot_cna_score),
			nrow = 1,
			ncol = 5,
			rel_widths = c(3, 1, 3, 3, 0.5)
		) %>% print
		dev.off()

		all_plot_data[[cohort]]$plot_data <- plot_data
		all_plot_data[[cohort]]$sig_genes <- c(immune_sigs, caf_endothelial_sigs)
		all_plot_data[[cohort]]$immune_cell_types <- immune_cell_types

	} else {

		sc_data <- sc_data[cell_type != 'ambiguous']

		gene_averages <- sapply(sc_data[, -c('id', 'patient', 'cell_type')], function(x) {log2(mean(10*(2^x - 1)) + 1)}, USE.NAMES = TRUE)

		sc_data <- sc_data[, c('id', 'patient', 'cell_type', names(gene_averages)[gene_averages >= 4]), with = FALSE]

		# Run t-SNE:
		# set.seed(cohort_data[[cohort]]$seeds[1])
		# sc_tsne <- Rtsne(
			# as.matrix(sc_data[, lapply(.SD, function(x) {x - mean(x)}), .SDcols = -c('id', 'patient', 'cell_type')])#,
			# # num_threads = 16
		# )

		# saveRDS(sc_tsne, paste0('../data_and_figures/sc_tsne_final/tsne_', cohort, '.rds'))

		sc_tsne <- readRDS(paste0('../data_and_figures/sc_tsne_final/tsne_', cohort, '.rds'))

		if('cells_to_exclude' %in% names(cohort_data[[cohort]])) {
			sc_tsne$Y <- sc_tsne$Y[-which(sc_data$id %in% cohort_data[[cohort]]$cells_to_exclude), ]
			sc_data <- sc_data[!(id %in% cohort_data[[cohort]]$cells_to_exclude)]
		}

		plot_data <- setNames(cbind(as.data.table(sc_tsne$Y), sc_data[, .(patient, cell_type)]), c('x', 'y', 'patient', 'cell_type'))

		tsne_plot_patient <- ggplot(plot_data, aes(x = x, y = y, colour = as.character(patient))) +
			geom_point() +
			theme_minimal() +
			labs(title = 'Patients', colour = NULL, x = 't-SNE 1', y = 't-SNE 2')

		tsne_plot_cell_types <- ggplot(plot_data, aes(x = x, y = y, colour = cell_type)) +
			geom_point() +
			scale_colour_manual(
				labels = cell_type_labels[sort(unique(plot_data$cell_type))],
				values = cell_type_colours[sort(unique(plot_data$cell_type))]
			) +
			theme_minimal() +
			labs(title = 'Cell types', colour = NULL, x = 't-SNE 1', y = 't-SNE 2')

		# CNA scores:

		cna_score_files <- dir(paste0('../data_and_figures/sc_find_malignant/', cohort))
		cna_score_files <- cna_score_files[endsWith(cna_score_files, '_data.csv')]

		cna_score_data <- lapply(
			cna_score_files,
			function(x) cbind(fread(paste0('../data_and_figures/sc_find_malignant/', cohort, '/', x)))
		) %>% rbindlist(fill = TRUE)

		setkey(cna_score_data, cell_id)

		cna_score_data <- cna_score_data[sc_data$id]

		cna_score_data[
			,
			classification_final := switch(
				(length(classification_final) == 1) + 1,
				switch( # In this case the cell must have been reference, but could have been classified as ambiguous or malignant in some patients.
					('ambiguous' %in% classification_final | 'malignant' %in% classification_final) + 1,
					'nonmalignant',
					'ambiguous'
				),
				classification_final # In this case the cell was not reference.
			),
			by = cell_id
		]

		cna_score_data <- cna_score_data[
			,
			.(
				classification_final = unique(classification_final),
				cna_score = switch(
					startsWith(unique(classification_final), 'malignant') + 1,
					mean(cna_signal),
					.SD[classification == classification_final, cna_signal]
				)
			),
			by = cell_id
		]

		cna_score_data[, patient := sc_data$patient]
		cna_score_data[, cna_score := cna_score/quantile(cna_score[startsWith(classification_final, 'malignant')], 0.75), by = patient]

		# cell_id should still be in the same order as sc_data$id, so it's safe to do this:
		plot_data[, cna_score := cna_score_data$cna_score]

		tsne_plot_cna_score <- ggplot(plot_data, aes(x = x, y = y, colour = cna_score)) +
			geom_point() +
			theme_minimal() +
			scale_colour_gradientn(colours = colorRamps::matlab.like(50), limits = c(0, 1), oob = scales::squish) +
			labs(title = 'CNA score', colour = NULL, x = 't-SNE 1', y = 't-SNE 2')

		immune_cell_types <- c('B', 'B_plasma', 'DC', 'macrophage', 'mast', 'T')[
			sapply(
				c('B', 'B_plasma', 'DC', 'macrophage', 'mast', 'T'),
				function(ct) cell_type_markers[cell_type == ct, sum(gene %in% names(sc_data)) >= 3]
			)
		]

		immune_cell_types <- immune_cell_types[
			mapvalues(
				immune_cell_types,
				c('B', 'B_plasma', 'DC', 'macrophage', 'mast', 'T'),
				c('b_cell', 'b_cell', 'dendritic', 'macrophage', 'mast', 't_cell'),
				warn_missing = FALSE
			) %in% sc_data$cell_type
		]

		set.seed(cohort_data[[cohort]]$seeds[2])

		immune_scores <- as.data.table(
			sapply(
				immune_cell_types,
				function(ct) {
					signature_score(
						sc_data[, set_colnames(t(.SD), id), .SDcols = -c('id', 'patient', 'cell_type')],
						sig_genes = cell_type_markers[cell_type == ct & gene %in% names(sc_data), unique(gene)],
						nbin = 30,
						n = 100
					)
				}
			),
			keep.rownames = 'cell_id'
		)

		immune_sigs <- sapply(
			immune_cell_types,
			function(ct) {

				sig_cor <- cell_type_markers[
					cell_type == ct & gene %in% names(sc_data),
					cor(sc_data[, unique(gene), with = FALSE], immune_scores[, ..ct])[, 1]
				]

				if(sum(sig_cor > 0.6) < 10) {
					if(length(sig_cor) <= 10) {
						return(names(sig_cor))
					} else {
						return(names(sig_cor)[names(sig_cor) %in% names(sort(sig_cor, decreasing = TRUE)[1:10])])
					}
				} else {
					return(names(sig_cor)[sig_cor > 0.6])
				}

			}
		)

		immune_scores <- as.data.table(
			sapply(
				immune_cell_types,
				function(ct) {
					signature_score(
						sc_data[, set_colnames(t(.SD), id), .SDcols = -c('id', 'patient', 'cell_type')],
						sig_genes = immune_sigs[[ct]],
						nbin = 30,
						n = 100
					)
				}
			),
			keep.rownames = 'cell_id'
		)

		immune_scores[
			,
			c('which_max', 'cell_type') := .(apply(.SD, 1, function(x) names(.SD)[which.max(x)]), sc_data$cell_type),
			.SDcols = -'cell_id'
		]

		immune_cell_type_scaling_factors <- sapply(
			immune_cell_types,
			function(x) {
				immune_scores[
					cell_type == mapvalues(
						x,
						c('B', 'B_plasma', 'DC', 'macrophage', 'mast', 'T'),
						c('b_cell', 'b_cell', 'dendritic', 'macrophage', 'mast', 't_cell'),
						warn_missing = FALSE
					),
					setNames(quantile(get(x), 0.9), NULL)
				]
			},
			USE.NAMES = TRUE
		)

		immune_scores[, score := get(unique(which_max))/immune_cell_type_scaling_factors[unique(which_max)], by = which_max]

		plot_data[
			,
			c('immune_score', 'caf_score', 'endothelial_score') := .(
				immune_scores$score,
				signature_score(
					sc_data[, set_colnames(t(.SD), id), .SDcols = -c('id', 'patient', 'cell_type')],
					sig_genes = cell_type_markers[cell_type == 'CAF' & gene %in% names(sc_data), unique(gene)],
					nbin = 30,
					n = 100
				),
				signature_score(
					sc_data[, set_colnames(t(.SD), id), .SDcols = -c('id', 'patient', 'cell_type')],
					sig_genes = cell_type_markers[cell_type == 'endothelial' & gene %in% names(sc_data), unique(gene)],
					nbin = 30,
					n = 100
				)
			)
		]

		caf_endothelial_sigs <- sapply(
			c('caf', 'endothelial'),
			function(ct) {

				sig_cor <- cell_type_markers[
					cell_type == mapvalues(ct, 'caf', 'CAF', warn_missing = FALSE) & gene %in% names(sc_data),
					cor(sc_data[, unique(gene), with = FALSE], plot_data[, get(paste0(ct, '_score'))])[, 1]
				]

				if(sum(sig_cor > 0.6) < 10) {
					if(length(sig_cor) <= 10) {
						return(names(sig_cor))
					} else {
						return(names(sig_cor)[names(sig_cor) %in% names(sort(sig_cor, decreasing = TRUE)[1:10])])
					}
				} else {
					return(names(sig_cor)[sig_cor > 0.6])
				}

			}
		)

		plot_data[
			,
			c('caf_score', 'endothelial_score') := .(
				signature_score(
					sc_data[, set_colnames(t(.SD), id), .SDcols = -c('id', 'patient', 'cell_type')],
					sig_genes = caf_endothelial_sigs$caf,
					nbin = 30,
					n = 100
				),
				signature_score(
					sc_data[, set_colnames(t(.SD), id), .SDcols = -c('id', 'patient', 'cell_type')],
					sig_genes = caf_endothelial_sigs$endothelial,
					nbin = 30,
					n = 100
				)
			)
		]

		plot_data[
			,
			c('caf_score', 'endothelial_score') := .(
				caf_score/.SD[cell_type == 'caf', quantile(caf_score, 0.9)],
				endothelial_score/.SD[cell_type == 'endothelial', quantile(endothelial_score, 0.9)]
			)
		]

		score_tsnes <- sapply(
			c('immune_score', 'caf_score', 'endothelial_score'),
			function(score_type) {
				ggplot(plot_data, aes(x = x, y = y, colour = get(score_type))) +
					geom_point() +
					theme_minimal() +
					scale_colour_gradientn(colours = colorRamps::matlab.like(50), limits = c(0, 1), oob = scales::squish) +
					labs(
						title = mapvalues(
							score_type,
							c('immune_score', 'caf_score', 'endothelial_score'),
							c('Immune score', 'CAF score', 'Endothelial score'),
							warn_missing = FALSE
						),
						colour = NULL,
						x = 't-SNE 1',
						y = 't-SNE 2'
					)
			},
			simplify = FALSE,
			USE.NAMES = TRUE
		)

		pdf(paste0('../data_and_figures/sc_tsne_final/tsne_plots_all_', cohort, '.pdf'))
		print(tsne_plot_cell_types)
		print(tsne_plot_patient)
		for(p in score_tsnes) print(p)
		print(tsne_plot_cna_score)
		dev.off()

		pdf(paste0('../data_and_figures/sc_tsne_final/tsne_plots_combined_', cohort, '.pdf'), width = 15, height = 5)
		plot_grid(
			tsne_plot_cell_types + theme(legend.position = 'none'),
			get_legend(tsne_plot_cell_types),
			score_tsnes$caf_score + theme(legend.position = 'none', axis.title.y = element_blank()),
			tsne_plot_cna_score + theme(legend.position = 'none', axis.title.y = element_blank()),
			get_legend(tsne_plot_cna_score),
			nrow = 1,
			ncol = 5,
			rel_widths = c(3, 1, 3, 3, 0.5)
		) %>% print
		dev.off()

		all_plot_data[[cohort]]$plot_data <- plot_data
		all_plot_data[[cohort]]$sig_genes <- c(immune_sigs, caf_endothelial_sigs)
		all_plot_data[[cohort]]$immune_cell_types <- immune_cell_types

	}

}





all_plots <- sapply(
	names(all_plot_data),
	function(cohort) {

		if('cell_type_for_plot' %in% names(all_plot_data[[cohort]]$plot_data)) {
			cell_type_plot <- ggplot(all_plot_data[[cohort]]$plot_data, aes(x = x, y = y, colour = cell_type_for_plot)) +
				geom_point() +
				scale_colour_manual(
					labels = cell_type_labels[sort(unique(all_plot_data[[cohort]]$plot_data$cell_type_for_plot))],
					values = cell_type_colours[sort(unique(all_plot_data[[cohort]]$plot_data$cell_type_for_plot))]
				) +
				theme_minimal() +
				labs(title = 'Cell types', colour = NULL, x = 't-SNE 1', y = 't-SNE 2')
		} else {
			cell_type_plot <- ggplot(all_plot_data[[cohort]]$plot_data, aes(x = x, y = y, colour = cell_type)) +
				geom_point() +
				scale_colour_manual(
					labels = cell_type_labels[sort(unique(all_plot_data[[cohort]]$plot_data$cell_type))],
					values = cell_type_colours[sort(unique(all_plot_data[[cohort]]$plot_data$cell_type))]
				) +
				theme_minimal() +
				labs(title = 'Cell types', colour = NULL, x = 't-SNE 1', y = 't-SNE 2')
		}

		caf_score_plot <- ggplot(all_plot_data[[cohort]]$plot_data, aes(x = x, y = y, colour = caf_score)) +
			geom_point() +
			theme_minimal() +
			scale_colour_gradientn(colours = colorRamps::matlab.like(50), limits = c(0, 1), oob = scales::squish) +
			labs(title = 'CAF score', colour = NULL, x = 't-SNE 1', y = 't-SNE 2')

		cna_score_plot <- ggplot(all_plot_data[[cohort]]$plot_data, aes(x = x, y = y, colour = cna_score)) +
			geom_point() +
			theme_minimal() +
			scale_colour_gradientn(colours = colorRamps::matlab.like(50), limits = c(0, 1), oob = scales::squish) +
			labs(title = 'CNA score', colour = NULL, x = 't-SNE 1', y = 't-SNE 2')

		list(cell_type_plot = cell_type_plot, caf_score_plot = caf_score_plot, cna_score_plot = cna_score_plot)

	},
	simplify = FALSE,
	USE.NAMES = TRUE
)

pdf('../data_and_figures/final_figures_resubmission/S2.pdf', width = 14, height = 20)

print(
	plot_grid(
		plotlist = unlist(
			lapply(
				names(all_plots[1:4]),
				function(cohort) {
					c(
						list(
							blank_plot() +
								labs(
									title = mapvalues(
										cohort,
										c('breast_qian', 'crc_lee_smc', 'hnscc_puram', 'liver_ma'),
										c('Breast - Qian et al.', 'Colorectal - Lee et al. - SMC cohort', 'Head and Neck - Puram et al.',
										  'Liver - Ma et al.'),
										warn_missing = FALSE
									)
								) +
								theme(plot.title = element_text(margin = margin(t = 30), size = 18), plot.margin = unit(c(5.5, 5.5, 10, 5.5), 'pt'))
						),
						list(
							plot_grid(
								all_plots[[cohort]]$cell_type_plot + theme(legend.position = 'none'),
								get_legend(all_plots[[cohort]]$cell_type_plot),
								all_plots[[cohort]]$caf_score_plot + theme(legend.position = 'none', axis.title.y = element_blank()),
								all_plots[[cohort]]$cna_score_plot + theme(legend.position = 'none', axis.title.y = element_blank()),
								get_legend(all_plots[[cohort]]$cna_score_plot),
								nrow = 1,
								ncol = 5,
								rel_widths = c(3.75, 1.75, 3.75, 3.75, 1)
							)
						)
					)
				}
			),
			recursive = FALSE
		),
		nrow = 8,
		ncol = 1,
		rel_heights = c(1, 4, 1, 4, 1, 4, 1, 4)
	)
)

print(
	plot_grid(
		plotlist = unlist(
			lapply(
				names(all_plots[5:8]),
				function(cohort) {
					c(
						list(
							blank_plot() +
								labs(
									title = mapvalues(
										cohort,
										c('luad_kim', 'lung_qian', 'ovarian_qian', 'pdac_peng'),
										c('Lung Adenocarcinoma - Kim et al.', 'Lung - Qian et al.', 'Ovarian - Qian et al.',
										  'Pancreatic - Peng et al.'),
										warn_missing = FALSE
									)
								) +
								theme(plot.title = element_text(margin = margin(t = 30), size = 18), plot.margin = unit(c(5.5, 5.5, 10, 5.5), 'pt'))
						),
						list(
							plot_grid(
								all_plots[[cohort]]$cell_type_plot + theme(legend.position = 'none'),
								get_legend(all_plots[[cohort]]$cell_type_plot),
								all_plots[[cohort]]$caf_score_plot + theme(legend.position = 'none', axis.title.y = element_blank()),
								all_plots[[cohort]]$cna_score_plot + theme(legend.position = 'none', axis.title.y = element_blank()),
								get_legend(all_plots[[cohort]]$cna_score_plot),
								nrow = 1,
								ncol = 5,
								rel_widths = c(3.75, 1.75, 3.75, 3.75, 1)
								# rel_widths = c(3, 1, 3, 3, 0.5)
							)
						)
					)
				}
			),
			recursive = FALSE
		),
		nrow = 8,
		ncol = 1,
		rel_heights = c(1, 4, 1, 4, 1, 4, 1, 4)
	)
)

dev.off()
