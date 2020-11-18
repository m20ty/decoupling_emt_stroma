# bsub -q tirosh -R "rusage[mem=64000]" \
	# -o ../data_and_figures/sc_tsne_final/sc_tsne_final_log.o \
	# -e ../data_and_figures/sc_tsne_final/sc_tsne_final_log.e \
	# Rscript sc_tsne_final.R

Sys.time()
cat(R.Version()$version.string, '\n')

library(data.table) # 1.12.8
library(ggplot2) # This is version 3.3.1 on my laptop's WSL setup but 3.3.0 on WEXAC.
library(magrittr) # 1.5
library(cowplot) # 1.0.0
library(Rtsne) # 0.15
library(plyr) # 1.8.6
library(stringr) # 1.4.0
library(colorspace) # 1.4.1

source('general_functions.R')

cell_type_markers <- fread('../../cell_type_markers.csv')[!(cell_type %in% c('mesenchymal', 'myocyte')) & source != 'TCGA_CCLE_comparison']

# The following merges the B and B_plasma marker gene lists.  I did this because we only have 5 B cell markers, which maybe
# isn't enough for the signature_score() function.  I undid it because it results in some B cell clusters having low immune score.
# cell_type_markers[cell_type == 'B_plasma', cell_type := 'B']

# Note that when I first did this, the macrophages often had high CAF scores.  I tried correcting this by subtracting from each score the average
# of the second highest scores, but this got me into a relativity mess...  It turns out the problem was simply that I was using too many "CAF
# markers", which was fixed by removing 'mesenchymal' and 'TCGA_CCLE_comparison' from cell_type_markers above.  I also removed 'myocyte', because
# I don't use these here.  While trying to fix the problem, I also built in a step to refine the cell type signatures by correlation with an
# initial cell type score.  I don't think this is necessary but it proably good practice, so I kept it.





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





# The following colours (used in the simulated bulk lineplots) are from Set3:
# 'b_cell' = '#8DD3C7',
# 'caf' = '#FDB462',
# 'cancer' = '#FB8072',
# 'dendritic' = '#FFED6F',
# 'endothelial' = '#BC80BD',
# 'macrophage' = '#80B1D3',
# 'mast' = '#FCCDE5',
# 't_cell' = '#B3DE69'

# The only ones we haven't used are 2, 3, 9 and 11.  Probably don't want to use 2.  Remaining cell types:
# epithelial
# myocyte
# hpc-like
# alveolar
# erythroblast
# stellate
# ductal (combine ductal_1 and ductal_2 in the PDAC dataset)
# acinar

# Too many for Set3!  Let's use randomcoloR to fix a set:
# set.seed(4323)
# random_colours <- randomcoloR::distinctColorPalette(17)
# To see the output:
# cat(paste0("'", paste(random_colours, collapse = "', '"), "'\n"))
# This is what we get:
# random_colours <- c('#799575', '#887DDA', '#D5978A', '#8C5581', '#79E5A3', '#75B0DC', '#D1D2DC', '#D8A354', '#DCA9DB',
					# '#E3E050', '#A83EE5', '#7AE1D8', '#DD4C61', '#DD6BCF', '#7AE759', '#DBE2BB', '#C0DF84')
# See what they look like:
# ggplot(data.table(x = letters[1:17])) +
	# geom_tile(aes(x = x, y = 0, fill = x)) +
	# scale_fill_manual(values = random_colours) +
	# scale_x_discrete(labels = random_colours, expand = c(0, 0)) +
	# scale_y_continuous(expand = c(0, 0)) +
	# theme(
		# axis.text.y = element_blank(),
		# axis.title.y = element_blank(),
		# axis.ticks.y = element_blank(),
		# axis.text.x = element_text(angle = 90, hjust = 1),
		# legend.position = 'none'
	# )

# Note I added NK cells after, so chose an additional colour using:
# set.seed(4211)
# randomColor(hue = 'green')

# To see this extra colour along with all the others:
# ggplot(data.table(x = letters[1:18])) +
	# geom_tile(aes(x = x, y = 0, fill = x)) +
	# scale_fill_manual(values = c(random_colours, '#01AD1E')) +
	# scale_x_discrete(labels = c(random_colours, '#01AD1E'), expand = c(0, 0)) +
	# scale_y_continuous(expand = c(0, 0)) +
	# theme(
		# axis.text.y = element_blank(),
		# axis.title.y = element_blank(),
		# axis.ticks.y = element_blank(),
		# axis.text.x = element_text(angle = 90, hjust = 1),
		# legend.position = 'none'
	# )

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
		
		sc_data <- sc_data[, c('id', 'patient', 'cell_type', 'cell_type_lenient', 'cell_type_for_plot', names(gene_averages)[gene_averages >= 4]), with = FALSE]
		
		# Run t-SNE:
		# set.seed(cohort_data[[cohort]]$seeds[1])
		# sc_tsne <- Rtsne(
			# as.matrix(sc_data[, lapply(.SD, function(x) {x - mean(x)}), .SDcols = -c('id', 'patient', 'cell_type', 'cell_type_lenient', 'cell_type_for_plot')])#,
			# # num_threads = 16
		# )
		
		# saveRDS(sc_tsne, paste0('../data_and_figures/sc_tsne_final/tsne_', cohort, '.rds'))
		
		sc_tsne <- readRDS(paste0('../data_and_figures/sc_tsne_final/tsne_', cohort, '.rds'))
		
		if('cells_to_exclude' %in% names(cohort_data[[cohort]])) {
			sc_tsne$Y <- sc_tsne$Y[-which(sc_data$id %in% cohort_data[[cohort]]$cells_to_exclude), ]
			sc_data <- sc_data[!(id %in% cohort_data[[cohort]]$cells_to_exclude)]
		}
		
		# If cohort == 'breast_qian' and we haven't removed sc5rJUQ064_CCATGTCCATCCCATC, we can see it is the one cancer cell in the
		# potential CAF cluster by the following:
		# sc_data[cell_type == 'cancer' & sc_tsne$Y[, 1] < -20 & sc_tsne$Y[, 2] < 0, id]
		# Similarly, if cohort == 'ovarian_qian' and we haven't removed scrSOL001_TCATTTGTCTGTCAAG and scrSOL004_TTGCCGTTCTCCTATA,
		# we can see these are the cancer cell and macrophage in the wrong clusters by the following:
		# sc_data[
			# (cell_type == 'cancer' & sc_tsne$Y[, 1] < -15 & sc_tsne$Y[, 2] > -5) |
				# (cell_type == 'macrophage' & sc_tsne$Y[, 1] < 10 & sc_tsne$Y[, 2] < -40),
			# id
		# ]
		
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
				switch( # In this case the cell must have been reference, but could have been classified as ambiguous or malignant in a subset of patients.
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
		
		# The following calculation of immune score is a bit complicated, but it's designed to put all the immune cells on roughly the same scale
		# so that we don't get some immune cell types having much higher scores than the others.  It works OK but not brilliantly, possibly because
		# taking the maximum makes it more noisy.
		
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
				
				sig_cor <- cell_type_markers[cell_type == ct & gene %in% names(sc_data), cor(sc_data[, unique(gene), with = FALSE], immune_scores[, ..ct])[, 1]]
				
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
			# which_max := apply(.SD, 1, function(x) names(.SD)[which.max(x)]),
			c('which_max', 'cell_type') := .(
				apply(.SD, 1, function(x) names(.SD)[which.max(x)]),
				sc_data$cell_type
			),
			.SDcols = -'cell_id'
		]
		
		# plot_data <- cbind(plot_data, immune_scores)
		# setcolorder(plot_data, 'cell_id')
		
		# set.seed(cohort_data[[cohort]]$seeds[3])
		# plot_data[
			# ,
			# c('caf_score', 'endothelial_score') := .(
				# signature_score(
					# sc_data[, set_colnames(t(.SD), id), .SDcols = -c('id', 'patient', 'cell_type', 'cell_type_lenient', 'cell_type_for_plot')],
					# sig_genes = cell_type_markers[cell_type == 'CAF' & gene %in% names(sc_data), unique(gene)],
					# nbin = 30,
					# n = 100
				# ),
				# signature_score(
					# sc_data[, set_colnames(t(.SD), id), .SDcols = -c('id', 'patient', 'cell_type', 'cell_type_lenient', 'cell_type_for_plot')],
					# sig_genes = cell_type_markers[cell_type == 'endothelial' & gene %in% names(sc_data), unique(gene)],
					# nbin = 30,
					# n = 100
				# )
			# )
		# ]
		
		# cell_type_score_ave <- plot_data[
			# ,
			# lapply(.SD, mean),
			# .SDcols = c(immune_cell_types, 'caf_score', 'endothelial_score', 'cna_score'),
			# by = .(ct = cell_type)
		# ]
		
		# plot_data[
			# ,
			# (c(immune_cell_types, 'caf_score', 'endothelial_score', 'cna_score')) := lapply(
				# .SD,
				# function(x) pmax(x - cell_type_score_ave[ct == unique(cell_type), sort(as.numeric(.SD), decreasing = TRUE)[2], .SDcols = -'ct'], 0)
			# ),
			# .SDcols = c(immune_cell_types, 'caf_score', 'endothelial_score', 'cna_score'),
			# by = cell_type
		# ]
		
		immune_cell_type_scaling_factors <- sapply(
			immune_cell_types,
			function(ict) {
				immune_scores[
				# plot_data[
					cell_type == mapvalues(
						ict,
						c('B', 'B_plasma', 'DC', 'macrophage', 'mast', 'T'),
						c('b_cell', 'b_cell', 'dendritic', 'macrophage', 'mast', 't_cell'),
						warn_missing = FALSE
					),
					setNames(quantile(get(ict), 0.9), NULL)
					# mean(get(ict))
				]
			},
			USE.NAMES = TRUE
		)
		
		# plot_data[
			# ,
			# immune_score := get(unique(which_max))/immune_cell_type_scaling_factors[unique(which_max)],
			# by = which_max
		# ]
		
		immune_scores[
			,
			score := get(unique(which_max))/immune_cell_type_scaling_factors[unique(which_max)],
			by = which_max
		]
		
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
			c('caf_score', 'endothelial_score') := .( # Should add pericyte score, but need to make my own signature
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
		
		# For each cell type, take the average of each score.  Then subtract the 2nd highest average from each of the scores, and take the
		# maximum of these centred scores and 0.  Then lower scores should become basically zero, and the 2nd highest score will become
		# centred at zero.  This will hopefully address the issue of macrophages having a high CAF score - this high CAF score will be
		# recentred to zero, but the immune score shouldn't be affected.
		
		# cell_type_score_ave <- plot_data[
			# ,
			# lapply(.SD, mean),
			# .SDcols = c('immune_score', 'caf_score', 'endothelial_score'),
			# by = .(ct = cell_type)
		# ]
		
		# plot_data[
			# ,
			# c('immune_score', 'caf_score', 'endothelial_score') := lapply(
				# .SD,
				# function(x) pmax(x - cell_type_score_ave[ct == unique(cell_type), sort(as.numeric(.SD))[2], .SDcols = -'ct'], 0)
			# ),
			# .SDcols = c('immune_score', 'caf_score', 'endothelial_score'),
			# by = cell_type
		# ]
		
		plot_data[
			,
			c('caf_score', 'endothelial_score') := .(
				caf_score/.SD[cell_type == 'caf', quantile(caf_score, 0.9)],
				endothelial_score/.SD[cell_type == 'endothelial', quantile(endothelial_score, 0.9)]
				# caf_score/.SD[cell_type == 'caf', mean(caf_score)],
				# endothelial_score/.SD[cell_type == 'endothelial', mean(endothelial_score)]
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
		
		gene_averages <- sapply(
			sc_data[, -c('id', 'patient', 'cell_type')],
			function(x) {log2(mean(10*(2^x - 1)) + 1)},
			USE.NAMES = TRUE
		)
		
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
		
		# If cohort == 'pdac_peng' and we haven't removed T8_TGGTTCCTCGCATGGC and T17_CGTGTAACAGTACACT, we can see these are the T cell and
		# CAF in the wrong clusters by the following:
		# sc_data[(cell_type == 't_cell' & sc_tsne$Y[, 1] > 35) | (cell_type == 'caf' & sc_tsne$Y[, 1] < -35), id]
		
		plot_data <- setNames(
			cbind(as.data.table(sc_tsne$Y), sc_data[, .(patient, cell_type)]),
			c('x', 'y', 'patient', 'cell_type')
		)
		
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
			function(x) cbind(fread(paste0('../data_and_figures/sc_find_malignant/', cohort, '/', x)))#, patient = gsub('_data.csv', '', x))
		) %>% rbindlist(fill = TRUE)
		
		setkey(cna_score_data, cell_id)
		
		cna_score_data <- cna_score_data[sc_data$id]
		
		cna_score_data[
			,
			classification_final := switch(
				(length(classification_final) == 1) + 1,
				switch( # In this case the cell must have been reference, but could have been classified as ambiguous or malignant in a subset of patients.
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
		
		# The following calculation of immune score is a bit complicated, but it's designed to put all the immune cells on roughly the same scale
		# so that we don't get some immune cell types having much higher scores than the others.  It works OK but not brilliantly, possibly because
		# taking the maximum makes it more noisy.
		
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
				
				sig_cor <- cell_type_markers[cell_type == ct & gene %in% names(sc_data), cor(sc_data[, unique(gene), with = FALSE], immune_scores[, ..ct])[, 1]]
				
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
			# which_max := apply(.SD, 1, function(x) names(.SD)[which.max(x)]),
			c('which_max', 'cell_type') := .(
				apply(.SD, 1, function(x) names(.SD)[which.max(x)]),
				sc_data$cell_type
			),
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
					# mean(get(x))
				]
			},
			USE.NAMES = TRUE
		)
		
		immune_scores[
			,
			score := get(unique(which_max))/immune_cell_type_scaling_factors[unique(which_max)],
			by = which_max
		]
		
		plot_data[
			,
			c('immune_score', 'caf_score', 'endothelial_score') := .( # Should add pericyte score, but need to make my own signature
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
			c('caf_score', 'endothelial_score') := .( # Should add pericyte score, but need to make my own signature
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
				# caf_score/.SD[cell_type == 'caf', mean(caf_score)],
				# endothelial_score/.SD[cell_type == 'endothelial', mean(endothelial_score)]
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
										c('Breast - Qian et al.', 'Colorectal - Lee et al. - SMC cohort', 'Head and Neck - Puram et al.', 'Liver - Ma et al.'),
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
										c('Lung Adenocarcinoma - Kim et al.', 'Lung - Qian et al.', 'Ovarian - Qian et al.', 'Pancreatic - Peng et al.'),
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
