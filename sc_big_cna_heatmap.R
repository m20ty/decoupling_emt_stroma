library(data.table) # 1.12.8
library(ggplot2) # This is version 3.3.1 on my laptop's WSL setup but 3.3.0 on WEXAC.
library(magrittr) # 1.5
library(cowplot) # 1.0.0
library(infercna) # 1.0.0

source('general_functions.R')

# Function from infercna:
.chromosomeBreaks = function(genes = NULL, halfway = F, hide = NULL) {
    n = length(genes)
    chrsum = cumsum(lengths(splitGenes(genes, by = 'chr')))
    Breaks = chrsum/max(chrsum) * n
    if (halfway) {
        b = Breaks
        Breaks = c(b[1]/2, b[2:length(b)] - diff(b)/2)
    }
    if (!is.null(hide)) {
        names(Breaks)[which(names(Breaks) %in% hide)] = rep('', length(hide))
    }
    Breaks
}





cohort_data <- list(
	breast_qian = list(
		read_quote = quote(fread('../data_and_figures/qian_breast_2020_reclassified.csv')),
		seed = 5347
	),
	crc_lee_smc = list(
		read_quote = quote(fread('../data_and_figures/lee_crc_2020_smc_reclassified.csv')),
		seed = 252
	),
    hnscc_puram = list(
		read_quote = quote(fread('../data_and_figures/puram_hnscc_2017_reclassified.csv')),
		seed = 5101
	),
	liver_ma = list(
		read_quote = quote(fread('../data_and_figures/ma_liver_2019_reclassified.csv')),
		seed = 756
	),
	luad_kim = list(
		read_quote = quote(fread('../data_and_figures/kim_luad_2020_reclassified.csv')),
		seed = 5481
	),
	lung_qian = list(
		read_quote = quote(fread('../data_and_figures/qian_lung_2020_reclassified.csv')[, -'disease']),
		seed = 1495
	),
	ovarian_qian = list(
		read_quote = quote(fread('../data_and_figures/qian_ovarian_2020_reclassified.csv')),
		seed = 1768
	),
    pdac_peng = list(
		read_quote = quote(fread('../data_and_figures/peng_pdac_2019_reclassified.csv')[startsWith(cell_type, 'ductal'), cell_type := 'ductal']),
		seed = 4574
	)
)





for(cohort in names(cohort_data)) {
	
	cat(paste0(cohort, '\n'))
	
	sc_data <- eval(cohort_data[[cohort]]$read_quote)[, -'cell_type_author']
	
	if('cell_type_lenient' %in% names(sc_data)) {
		sc_data <- sc_data[cell_type_lenient %in% c('caf', 'cancer')]
	} else {
		sc_data <- sc_data[cell_type %in% c('caf', 'cancer')]
	}
	
	inferred_cna <- readRDS(paste0('../data_and_figures/inferred_cna_', cohort, '.rds'))
	inferred_cna <- inferred_cna[names(inferred_cna) %in% sc_data$patient]
	
	# To adjust genes on chromosome 6p (increasing noise to account for bias from immune cell reference):
	# genes_split_by_arm <- splitGenes(rownames(inferred_cna[[1]]), by = 'arm')
	# for(p in names(inferred_cna)) {
		# inferred_cna[[p]][genes_split_by_arm[['6p']], ][
			# inferred_cna[[p]][genes_split_by_arm[['6p']], ] > -0.3 & inferred_cna[[p]][genes_split_by_arm[['6p']], ] < 0
		# ] <- 0
	# }
	# I think it's better to do this later, though, directly in plot_data.
	
	x_line_breaks <- .chromosomeBreaks(rownames(inferred_cna[[1]]))
	x_text_breaks <- round(.chromosomeBreaks(rownames(inferred_cna[[1]]), halfway = TRUE, hide = c('13', '18', '21', 'Y')))
	
	cna_signal_files <- dir(paste0('../data_and_figures/sc_find_malignant/', cohort))
	cna_signal_files <- cna_signal_files[endsWith(cna_signal_files, '_data.csv')]
	
	classification_data <- lapply(
		cna_signal_files,
		function(x) cbind(fread(paste0('../data_and_figures/sc_find_malignant/', cohort, '/', x)))#, patient = gsub('_data.csv', '', x))
	) %>% rbindlist(fill = TRUE)
	
	setkey(classification_data, cell_id)
	
	classification_data <- classification_data[sc_data[cell_type == 'cancer', id]]
	classification_data <- classification_data[classification == classification_final]
	
	# classification_data[
		# ,
		# classification_final := switch(
			# (length(classification_final) == 1) + 1,
			# switch( # In this case the cell must have been reference, but could have been classified as ambiguous or malignant in a subset of patients.
				# # EDIT: I don't think this is true!  A cell could have not been reference and been tested against two subclones.  Then it might be classified
				# # as malignant for one clone but nonmalignant or ambiguous for another.  So actually, in this case, since I'm taking only the cells I
				# # eventually classified as malignant, classification_final will have length > 1 iff unique(clone) has length > 1.  So it should be enough to
				# # check that classification doesn't consist of both 'malignant clone 1' and 'malignant clone 2' (for example), in which case it should be
				# # ambiguous.  Hence I'm changing the following.  I should probably check all my classifications.  Perhaps the reason for the cancer cell in
				# # the CAF cluster of the breast_qian final t-SNE is a misclassification...
				# # ('ambiguous' %in% classification_final | 'malignant' %in% classification_final) + 1,
				# # 'nonmalignant',
				# # 'ambiguous'
				# (sum(startsWith(classification, 'malignant')) > 1) + 1,
				# ,
				# 'ambiguous'
			# ),
			# classification_final # In this case the cell was not reference.
		# ),
		# by = cell_id
	# ]
	
	# For some reason, the key sometimes disappears after the above two lines, so I'm setting it again:
	setkey(classification_data, cell_id)
	
	classification_data <- sc_data[
		,
		.(
			patient = patient,
			classification = switch((cell_type == 'cancer') + 1, 'caf', classification_data[id, classification_final])
		),
		by = id
	]
	
	# The following deals with cases where all cells of one subclone get changed to 'ambiguous' at some point, presumably because their
	# I disagree with the authors as to their malignancy.  This happened in HNSCC, where all cells in malignant clone 2 became 'ambiguous'
	# leaving only 'Patient 6: clone 1' in the final heatmap, which obviously looks odd.
	classification_data[
		startsWith(classification, 'malignant') & patient %in% classification_data[
			startsWith(classification, 'malignant'),
			.(dodgy = any(grepl('clone', classification)) & length(unique(classification)) == 1),
			by = patient
		][dodgy == TRUE, patient],
		classification := 'malignant'
	]
	
	plot_data <- rbindlist(
		lapply(
			unique(classification_data$patient),
			function(p) cbind(classification_data[patient == p], t(inferred_cna[[as.character(p)]][, classification_data[patient == p, id]]))
		)
	)
	
	# Downsample where a patient is overrepresented:
	# Note the final proportion won't be 2/#patients, because you're taking a sample of size 2*sum(N)/#patients, but then of course sum(N)
	# becomes smaller.
	set.seed(cohort_data[[cohort]]$seeds)
	downsampled_ids <- plot_data[, .N, by = .(patient, classification)][
		,
		setNames(lapply(list(patient, classification, N), function(x) x[N > 2*sum(N)/length(unique(patient))]), c('patient_id', 'cl_full', 'N')),
		by = .(cl = gsub(' clone [0-9]$', '', classification))
	][
		,
		.(
			sampled_ids = plot_data[
				startsWith(classification, cl),
				sample(id[patient == patient_id & classification == cl_full], 2*.N/length(unique(patient)))
			]
		),
		# plot_data[, sample(id[patient == patient_id], 3*sum(startsWith(classification, cl))/length(unique(patient)))],
		by = .(cl, patient_id, cl_full)
	]
	
	downsampled_ids <- plot_data[
		,
		.(
			sampled_ids = switch(
				(patient %in% downsampled_ids$patient_id && classification %in% downsampled_ids[patient_id == patient, cl_full]) + 1,
				id,
				downsampled_ids[patient_id == patient & cl_full == classification, sampled_ids]
			)
		),
		by = .(patient, classification)
	]$sampled_ids
	
	plot_data <- plot_data[id %in% downsampled_ids]
	classification_data <- classification_data[id %in% downsampled_ids]
	
	plot_data <- melt(
		plot_data,
		id.vars = c('id', 'patient', 'classification'),
		variable.name = 'gene',
		value.name = 'cna_score'
	)
	
	# Heatmap for cancer cells:
	
	y_breaks_cancer_major <- classification_data[classification != 'caf'][order(patient), .(brks = .N), by = patient]$brks %>% cumsum
	
	y_breaks_cancer_minor <- numeric(0)
	temp <- classification_data[classification != 'caf'][order(patient, classification)]
	for(i in 2:nrow(temp)) {
		if(
			startsWith(temp[i - 1, classification], 'malignant clone') &&
				temp[i - 1, classification] != temp[i, classification] &&
				startsWith(temp[i, classification], 'malignant clone') &&
				temp[i - 1, patient] == temp[i, patient]
		) {y_breaks_cancer_minor <- c(y_breaks_cancer_minor, i)}
	}
	rm(temp)
	
	all_y_breaks <- sort(c(0, y_breaks_cancer_major, y_breaks_cancer_minor))
	y_text_breaks_cancer <- sapply(2:length(all_y_breaks), function(i) round(all_y_breaks[i - 1] + (all_y_breaks[i] - all_y_breaks[i - 1])/2))
	
	y_breaks_cancer_major <- y_breaks_cancer_major[1:(length(y_breaks_cancer_major) - 1)]
	
	cna_heatmap_cancer <- ggplot(plot_data[classification != 'caf']) +
		geom_raster(
			aes(
				x = factor(gene, levels = rownames(inferred_cna[[1]])),
				y = factor(id, levels = classification_data[classification != 'caf'][order(patient, classification, id), unique(id)]),
				fill = cna_score
			)
		) +
		scale_fill_gradientn(
			colours = rev(colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))(50)),
			limits = c(-1, 1),
			oob = scales::squish
		) +
		scale_y_discrete(
			expand = c(0, 0),
			breaks = classification_data[classification != 'caf'][order(patient, classification, id), id[y_text_breaks_cancer]],
			labels = classification_data[classification != 'caf'][
				order(patient, classification),
				sapply(
					y_text_breaks_cancer,
					function(x) switch(
						(classification[x] == 'malignant') + 1,
						paste0('Patient ', patient[x], ':', gsub('malignant', '', classification[x])),
						paste0('Patient ', patient[x])
					)
				)
			]
		) +
		scale_x_discrete(
			expand = c(0, 0),
			breaks = rownames(inferred_cna[[1]])[x_text_breaks],
			labels = names(x_text_breaks)
		) +
		geom_vline(xintercept = x_line_breaks, size = 0.25) +
		geom_hline(yintercept = y_breaks_cancer_major, size = 0.25) +
		geom_hline(yintercept = y_breaks_cancer_minor, size = 0.25, linetype = 'dashed') +
		theme(
			axis.ticks = element_blank(),
			axis.ticks.length = unit(0, 'pt'),
			panel.border = element_rect(colour = 'black', size = 0.5, fill = NA)
		) +
		labs(x = 'Chromosome', y = 'Cells', fill = 'Inferred CNA')
	
	patient_annotations_cancer <- character(classification_data[classification != 'caf', .N])
	patient_annotations_cancer[y_text_breaks_cancer] <- classification_data[classification != 'caf'][
		order(patient, classification),
		sapply(
			y_text_breaks_cancer,
			function(x) switch(
				(classification[x] == 'malignant') + 1,
				paste0('Patient ', patient[x], ':', gsub('malignant', '', classification[x])),
				paste0('Patient ', patient[x])
			)
		)
	]
	annotations_cancer <- heat_map_labels_repel(
		patient_annotations_cancer,
		edge = 'right',
		nudge = 0.1,
		axis_title = 'Cancer cells',
		title_edge_margin = 5.5
	)
	
	# Heatmap for CAFs:
	
	y_breaks_caf <- classification_data[classification == 'caf'][order(patient), .(brks = .N), by = patient]$brks %>% cumsum
	
	y_text_breaks_caf <- sapply(
		2:(length(y_breaks_caf) + 1),
		function(i) round(c(0, y_breaks_caf)[i - 1] + (c(0, y_breaks_caf)[i] - c(0, y_breaks_caf)[i - 1])/2)
	)
	
	y_breaks_caf <- y_breaks_caf[1:(length(y_breaks_caf) - 1)]
	
	cna_heatmap_caf <- ggplot(plot_data[classification == 'caf']) +
		geom_raster(
			aes(
				x = factor(gene, levels = rownames(inferred_cna[[1]])),
				y = factor(id, levels = classification_data[classification == 'caf'][order(patient, classification, id), unique(id)]),
				fill = cna_score
			)
		) +
		scale_fill_gradientn(
			colours = rev(colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))(50)),
			limits = c(-1, 1),
			oob = scales::squish
		) +
		scale_y_discrete(
			expand = c(0, 0),
			breaks = classification_data[classification == 'caf'][order(patient, classification, id), id[y_text_breaks_caf]],
			labels = classification_data[classification == 'caf'][
				order(patient, classification),
				sapply(y_text_breaks_caf, function(x) paste0('Patient ', patient[x]))
			]
		) +
		scale_x_discrete(
			expand = c(0, 0),
			breaks = rownames(inferred_cna[[1]])[x_text_breaks],
			labels = names(x_text_breaks)
		) +
		geom_vline(xintercept = x_line_breaks, size = 0.25) +
		geom_hline(yintercept = y_breaks_caf, size = 0.25) +
		theme(
			axis.ticks = element_blank(),
			axis.ticks.length = unit(0, 'pt'),
			panel.border = element_rect(colour = 'black', size = 0.5, fill = NA)
		) +
		labs(x = 'Chromosome', y = 'Cells', fill = 'Inferred CNA')
	
	patient_annotations_caf <- character(classification_data[classification == 'caf', .N])
	patient_annotations_caf[y_text_breaks_caf] <- classification_data[classification == 'caf'][
		order(patient, classification),
		sapply(y_text_breaks_caf, function(x) paste0('Patient ', patient[x]))
	]
	annotations_caf <- heat_map_labels_repel(
		patient_annotations_caf,
		edge = 'right',
		nudge = 0.1,
		axis_title = 'CAFs',
		title_edge_margin = 5.5
	)
	
	# Combined plot:
	
	pdf(paste0('../data_and_figures/sc_big_cna_heatmap/cna_heatmap_', cohort, '.pdf'), width = 10, height = 10)
	
	# plot_grid(
		# plot_grid(
			# cna_heatmap_caf + theme(axis.text.x = element_blank(), axis.title.x = element_blank(), legend.position = 'none'),
			# cna_heatmap_cancer + theme(legend.position = 'none'),
			# nrow = 2,
			# ncol = 1,
			# align = 'v',
			# rel_heights = c(1, 2)
		# ),
		# get_legend(cna_heatmap_caf),
		# nrow = 1,
		# ncol = 2,
		# rel_widths = c(7, 1)
	# ) %>% print
	
	print(
		plot_grid(
			plot_grid(
				annotations_caf,
				cna_heatmap_caf + theme(
					axis.text = element_blank(),
					axis.title = element_blank(),
					legend.position = 'none',
					plot.margin = unit(c(5.5, 5.5, 0, 0), 'pt')
				),
				annotations_cancer,
				cna_heatmap_cancer + theme(
					axis.text.y = element_blank(),
					axis.title.y = element_blank(),
					legend.position = 'none',
					plot.margin = unit(c(0, 5.5, 5.5, 0), 'pt')
				),
				align = 'h',
				nrow = 2,
				ncol = 2,
				rel_widths = c(1, 4),
				rel_heights = c(1, 2)
			),
			get_legend(cna_heatmap_caf),
			nrow = 1,
			ncol = 2,
			rel_widths = c(7, 1)
		)
	)
	
	dev.off()
	
	# Make heatmaps again but adjusting genes on chromosome 6p (increasing noise to account for bias from immune cell reference):
	genes_split_by_arm <- splitGenes(rownames(inferred_cna[[1]]), by = 'arm')
	plot_data[gene %in% genes_split_by_arm[['6p']] & cna_score < 0 & cna_score > -0.35, cna_score := 0]
	
	cna_heatmap_cancer <- ggplot(plot_data[classification != 'caf']) +
		geom_raster(
			aes(
				x = factor(gene, levels = rownames(inferred_cna[[1]])),
				y = factor(id, levels = classification_data[classification != 'caf'][order(patient, classification, id), unique(id)]),
				fill = cna_score
			)
		) +
		scale_fill_gradientn(
			colours = rev(colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))(50)),
			limits = c(-1, 1),
			oob = scales::squish
		) +
		scale_y_discrete(
			expand = c(0, 0),
			breaks = classification_data[classification != 'caf'][order(patient, classification, id), id[y_text_breaks_cancer]],
			labels = classification_data[classification != 'caf'][
				order(patient, classification),
				sapply(
					y_text_breaks_cancer,
					function(x) switch(
						(classification[x] == 'malignant') + 1,
						paste0('Patient ', patient[x], ':', gsub('malignant', '', classification[x])),
						paste0('Patient ', patient[x])
					)
				)
			]
		) +
		scale_x_discrete(
			expand = c(0, 0),
			breaks = rownames(inferred_cna[[1]])[x_text_breaks],
			labels = names(x_text_breaks)
		) +
		geom_vline(xintercept = x_line_breaks, size = 0.25) +
		geom_hline(yintercept = y_breaks_cancer_major, size = 0.25) +
		geom_hline(yintercept = y_breaks_cancer_minor, size = 0.25, linetype = 'dashed') +
		theme(
			axis.ticks = element_blank(),
			axis.ticks.length = unit(0, 'pt'),
			panel.border = element_rect(colour = 'black', size = 0.5, fill = NA)
		) +
		labs(x = 'Chromosome', y = 'Cells', fill = 'Inferred CNA')
	
	cna_heatmap_caf <- ggplot(plot_data[classification == 'caf']) +
		geom_raster(
			aes(
				x = factor(gene, levels = rownames(inferred_cna[[1]])),
				y = factor(id, levels = classification_data[classification == 'caf'][order(patient, classification, id), unique(id)]),
				fill = cna_score
			)
		) +
		scale_fill_gradientn(
			colours = rev(colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))(50)),
			limits = c(-1, 1),
			oob = scales::squish
		) +
		scale_y_discrete(
			expand = c(0, 0),
			breaks = classification_data[classification == 'caf'][order(patient, classification, id), id[y_text_breaks_caf]],
			labels = classification_data[classification == 'caf'][
				order(patient, classification),
				sapply(y_text_breaks_caf, function(x) paste0('Patient ', patient[x]))
			]
		) +
		scale_x_discrete(
			expand = c(0, 0),
			breaks = rownames(inferred_cna[[1]])[x_text_breaks],
			labels = names(x_text_breaks)
		) +
		geom_vline(xintercept = x_line_breaks, size = 0.25) +
		geom_hline(yintercept = y_breaks_caf, size = 0.25) +
		theme(
			axis.ticks = element_blank(),
			axis.ticks.length = unit(0, 'pt'),
			panel.border = element_rect(colour = 'black', size = 0.5, fill = NA)
		) +
		labs(x = 'Chromosome', y = 'Cells', fill = 'Inferred CNA')
	
	pdf(paste0('../data_and_figures/sc_big_cna_heatmap/cna_heatmap_6padj_', cohort, '.pdf'), width = 10, height = 10)
	
	print(
		plot_grid(
			plot_grid(
				annotations_caf,
				cna_heatmap_caf + theme(
					axis.text = element_blank(),
					axis.title = element_blank(),
					legend.position = 'none',
					plot.margin = unit(c(5.5, 5.5, 0, 0), 'pt')
				),
				annotations_cancer,
				cna_heatmap_cancer + theme(
					axis.text.y = element_blank(),
					axis.title.y = element_blank(),
					legend.position = 'none',
					plot.margin = unit(c(0, 5.5, 5.5, 0), 'pt')
				),
				align = 'h',
				nrow = 2,
				ncol = 2,
				rel_widths = c(1, 4),
				rel_heights = c(1, 2)
			),
			get_legend(cna_heatmap_caf),
			nrow = 1,
			ncol = 2,
			rel_widths = c(7, 1)
		)
	)
	
	dev.off()
	
}





# The following doesn't work - need to save the plot data from each iteration of the for loop, and remake the figures outside the loop.
# But I'm not sure this is necessary anyway - I think the figures can stay as individual PDFs.

# Make combined figures all in a single PDF:

# pdf('../data_and_figures/sc_big_cna_heatmap/cna_heatmap_all.pdf', width = 10, height = 14)

# for(cohort_figs in figure_list) {
	# print(
		# plot_grid(
			# plot_grid(
				# cohort_figs$caf$annotations,
				# cohort_figs$caf$heatmap + theme(
					# axis.text = element_blank(),
					# axis.title = element_blank(),
					# legend.position = 'none',
					# plot.margin = unit(c(5.5, 5.5, 0, 0), 'pt')
				# ),
				# cohort_figs$cancer$annotations,
				# cohort_figs$cancer$heatmap + theme(
					# axis.text.y = element_blank(),
					# axis.title.y = element_blank(),
					# legend.position = 'none',
					# plot.margin = unit(c(0, 5.5, 5.5, 0), 'pt')
				# ),
				# align = 'h',
				# nrow = 2,
				# ncol = 2,
				# rel_widths = c(1, 4),
				# rel_heights = c(1, 2)
			# ),
			# get_legend(cohort_figs$caf$heatmap),
			# nrow = 1,
			# ncol = 2,
			# rel_widths = c(7, 1)
		# )
	# )
# }

# dev.off()

# # Same for 6p-adjusted:

# pdf('../data_and_figures/sc_big_cna_heatmap/cna_heatmap_all_6padj.pdf', width = 10, height = 14)

# for(cohort_figs in figure_list) {
	# print(
		# plot_grid(
			# plot_grid(
				# cohort_figs$caf$annotations,
				# cohort_figs$caf$heatmap_6padj + theme(
					# axis.text = element_blank(),
					# axis.title = element_blank(),
					# legend.position = 'none',
					# plot.margin = unit(c(5.5, 5.5, 0, 0), 'pt')
				# ),
				# cohort_figs$cancer$annotations,
				# cohort_figs$cancer$heatmap_6padj + theme(
					# axis.text.y = element_blank(),
					# axis.title.y = element_blank(),
					# legend.position = 'none',
					# plot.margin = unit(c(0, 5.5, 5.5, 0), 'pt')
				# ),
				# align = 'h',
				# nrow = 2,
				# ncol = 2,
				# rel_widths = c(1, 4),
				# rel_heights = c(1, 2)
			# ),
			# get_legend(cohort_figs$caf$heatmap_6padj),
			# nrow = 1,
			# ncol = 2,
			# rel_widths = c(7, 1)
		# )
	# )
# }

# dev.off()
