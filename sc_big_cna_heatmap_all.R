# bsub -q tirosh -R rusage[mem=128000] -o sc_big_cna_heatmap_all.o -e sc_big_cna_heatmap_all.e Rscript sc_big_cna_heatmap_all.R

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
		seed = 5347,
		plot_title = 'Breast - Qian et al.'
	),
	crc_lee_smc = list(
		read_quote = quote(fread('../data_and_figures/lee_crc_2020_smc_reclassified.csv')),
		seed = 252,
		plot_title = 'Colorectal - Lee et al. - SMC cohort'
	),
    hnscc_puram = list(
		read_quote = quote(fread('../data_and_figures/puram_hnscc_2017_reclassified.csv')),
		seed = 5101,
		plot_title = 'Head and Neck - Puram et al.'
	),
	liver_ma = list(
		read_quote = quote(fread('../data_and_figures/ma_liver_2019_reclassified.csv')),
		seed = 756,
		plot_title = 'Liver - Ma et al.'
	),
	luad_kim = list(
		read_quote = quote(fread('../data_and_figures/kim_luad_2020_reclassified.csv')),
		seed = 5481,
		plot_title = 'Lung Adenocarcinoma - Kim et al.'
	),
	lung_qian = list(
		read_quote = quote(fread('../data_and_figures/qian_lung_2020_reclassified.csv')[, -'disease']),
		seed = 1495,
		plot_title = 'Lung - Qian et al.'
	),
	ovarian_qian = list(
		read_quote = quote(fread('../data_and_figures/qian_ovarian_2020_reclassified.csv')),
		seed = 1768,
		plot_title = 'Ovarian - Qian et al.'
	),
    pdac_peng = list(
		read_quote = quote(fread('../data_and_figures/peng_pdac_2019_reclassified.csv')[startsWith(cell_type, 'ductal'), cell_type := 'ductal']),
		seed = 4574,
		plot_title = 'Pancreatic - Peng et al.'
	)
)





all_plot_data <- sapply(
	names(cohort_data),
	function(cohort) {
		
		cat(paste0(cohort, '\n'))
		
		sc_data <- eval(cohort_data[[cohort]]$read_quote)[, -'cell_type_author']
		
		if('cell_type_lenient' %in% names(sc_data)) {
			sc_data <- sc_data[cell_type_lenient %in% c('caf', 'cancer')]
		} else {
			sc_data <- sc_data[cell_type %in% c('caf', 'cancer')]
		}
		
		inferred_cna <- readRDS(paste0('../data_and_figures/inferred_cna_', cohort, '.rds'))
		inferred_cna <- inferred_cna[names(inferred_cna) %in% sc_data$patient]
		
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
		plot_data[, gene := factor(gene, levels = rownames(inferred_cna[[1]]))]
		
		# For cancer cell heatmap:
		
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
		
		# For CAF heatmap:
		
		y_breaks_caf <- classification_data[classification == 'caf'][order(patient), .(brks = .N), by = patient]$brks %>% cumsum
		
		y_text_breaks_caf <- sapply(
			2:(length(y_breaks_caf) + 1),
			function(i) round(c(0, y_breaks_caf)[i - 1] + (c(0, y_breaks_caf)[i] - c(0, y_breaks_caf)[i - 1])/2)
		)
		
		y_breaks_caf <- y_breaks_caf[1:(length(y_breaks_caf) - 1)]
		
		patient_annotations_caf <- character(classification_data[classification == 'caf', .N])
		patient_annotations_caf[y_text_breaks_caf] <- classification_data[classification == 'caf'][
			order(patient, classification),
			sapply(y_text_breaks_caf, function(x) paste0('Patient ', patient[x]))
		]
		
		# For adjustment of chromosome 6p:
		genes_split_by_arm <- splitGenes(rownames(inferred_cna[[1]]), by = 'arm')
		
		list(
			plot_data = plot_data,
			classification_data = classification_data,
			patient_annotations_cancer = patient_annotations_cancer,
			y_text_breaks_cancer = y_text_breaks_cancer,
			y_breaks_cancer_major = y_breaks_cancer_major,
			y_breaks_cancer_minor = y_breaks_cancer_minor,
			patient_annotations_caf = patient_annotations_caf,
			y_breaks_caf = y_breaks_caf,
			y_text_breaks_caf = y_text_breaks_caf,
			x_line_breaks = x_line_breaks,
			x_text_breaks = x_text_breaks,
			genes_split_by_arm = genes_split_by_arm
		)
		
	},
	simplify = FALSE,
	USE.NAMES = TRUE
)





pdf('../data_and_figures/final_figures_resubmission/S1.pdf', width = 10, height = 12)
# pdf('../data_and_figures/sc_big_cna_heatmap/cna_heatmap_all.pdf', width = 10, height = 12)

for(cohort in names(all_plot_data)) {
	
	cna_heatmap_cancer <- ggplot(all_plot_data[[cohort]]$plot_data[classification != 'caf']) +
		geom_raster(
			aes(
				x = gene,
				y = factor(id, levels = all_plot_data[[cohort]]$classification_data[classification != 'caf'][order(patient, classification, id), unique(id)]),
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
			breaks = all_plot_data[[cohort]]$classification_data[classification != 'caf'][order(patient, classification, id), id[all_plot_data[[cohort]]$y_text_breaks_cancer]],
			labels = all_plot_data[[cohort]]$classification_data[classification != 'caf'][
				order(patient, classification),
				sapply(
					all_plot_data[[cohort]]$y_text_breaks_cancer,
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
			breaks = all_plot_data[[cohort]]$plot_data[classification != 'caf', unique(gene)][all_plot_data[[cohort]]$x_text_breaks],
			labels = names(all_plot_data[[cohort]]$x_text_breaks)
		) +
		geom_vline(xintercept = all_plot_data[[cohort]]$x_line_breaks, size = 0.25) +
		geom_hline(yintercept = all_plot_data[[cohort]]$y_breaks_cancer_major, size = 0.25) +
		geom_hline(yintercept = all_plot_data[[cohort]]$y_breaks_cancer_minor, size = 0.25, linetype = 'dashed') +
		theme(
			axis.ticks = element_blank(),
			axis.ticks.length = unit(0, 'pt'),
			panel.border = element_rect(colour = 'black', size = 0.5, fill = NA)
		) +
		labs(x = 'Chromosome', y = 'Cells', fill = 'Inferred CNA')
	
	annotations_cancer <- heat_map_labels_repel(
		all_plot_data[[cohort]]$patient_annotations_cancer,
		edge = 'right',
		nudge = 0.1,
		axis_title = 'Cancer cells',
		title_edge_margin = 5.5
	)
	
	cna_heatmap_caf <- ggplot(all_plot_data[[cohort]]$plot_data[classification == 'caf']) +
		geom_raster(
			aes(
				x = gene,
				y = factor(id, levels = all_plot_data[[cohort]]$classification_data[classification == 'caf'][order(patient, classification, id), unique(id)]),
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
			breaks = all_plot_data[[cohort]]$classification_data[classification == 'caf'][order(patient, classification, id), id[all_plot_data[[cohort]]$y_text_breaks_caf]],
			labels = all_plot_data[[cohort]]$classification_data[classification == 'caf'][
				order(patient, classification),
				sapply(all_plot_data[[cohort]]$y_text_breaks_caf, function(x) paste0('Patient ', patient[x]))
			]
		) +
		scale_x_discrete(
			expand = c(0, 0),
			breaks = all_plot_data[[cohort]]$plot_data[classification == 'caf', unique(gene)][all_plot_data[[cohort]]$x_text_breaks],
			labels = names(all_plot_data[[cohort]]$x_text_breaks)
		) +
		geom_vline(xintercept = all_plot_data[[cohort]]$x_line_breaks, size = 0.25) +
		geom_hline(yintercept = all_plot_data[[cohort]]$y_breaks_caf, size = 0.25) +
		theme(
			axis.ticks = element_blank(),
			axis.ticks.length = unit(0, 'pt'),
			panel.border = element_rect(colour = 'black', size = 0.5, fill = NA),
			plot.title = element_text(margin = margin(b = 30), size = 18)
		) +
		labs(x = 'Chromosome', y = 'Cells', fill = 'Inferred CNA', title = cohort_data[[cohort]]$plot_title)
	
	annotations_caf <- heat_map_labels_repel(
		all_plot_data[[cohort]]$patient_annotations_caf,
		edge = 'right',
		nudge = 0.1,
		axis_title = 'CAFs',
		title_edge_margin = 5.5
	)
	
	print(
		plot_grid(
			plot_grid(
				plot_grid(
					annotations_caf,
					cna_heatmap_caf + theme(
						axis.text = element_blank(),
						axis.title = element_blank(),
						legend.position = 'none',
						plot.margin = unit(c(5.5, 5.5, 0, 0), 'pt')
					),
					nrow = 1,
					ncol = 2,
					align = 'h',
					rel_widths = c(1, 4)
				),
				plot_grid(
					annotations_cancer,
					cna_heatmap_cancer + theme(
						axis.text.y = element_blank(),
						axis.title.y = element_blank(),
						legend.position = 'none',
						plot.margin = unit(c(40, 5.5, 5.5, 0), 'pt')
					),
					nrow = 1,
					ncol = 2,
					align = 'h',
					rel_widths = c(1, 4)
				),
				nrow = 2,
				ncol = 1,
				rel_heights = c(2, 3)
				# rel_heights = c(1, 2)
			),
			get_legend(cna_heatmap_caf),
			nrow = 1,
			ncol = 2,
			rel_widths = c(7, 1)
		)
	)
	
}

dev.off()





# Make heatmaps again but adjusting genes on chromosome 6p (increasing noise to account for bias from immune cell reference):

pdf('../data_and_figures/sc_big_cna_heatmap/cna_heatmap_all_6padj.pdf', width = 10, height = 12)

for(cohort in names(all_plot_data)) {
	
	plot_data <- copy(all_plot_data[[cohort]]$plot_data)
	plot_data[gene %in% all_plot_data[[cohort]]$genes_split_by_arm[['6p']] & cna_score < 0 & cna_score > -0.35, cna_score := 0]
	
	cna_heatmap_cancer_6padj <- ggplot(plot_data[classification != 'caf']) +
		geom_raster(
			aes(
				x = gene,
				y = factor(id, levels = all_plot_data[[cohort]]$classification_data[classification != 'caf'][order(patient, classification, id), unique(id)]),
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
			breaks = all_plot_data[[cohort]]$classification_data[classification != 'caf'][order(patient, classification, id), id[all_plot_data[[cohort]]$y_text_breaks_cancer]],
			labels = all_plot_data[[cohort]]$classification_data[classification != 'caf'][
				order(patient, classification),
				sapply(
					all_plot_data[[cohort]]$y_text_breaks_cancer,
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
			breaks = plot_data[classification != 'caf', unique(gene)][all_plot_data[[cohort]]$x_text_breaks],
			labels = names(all_plot_data[[cohort]]$x_text_breaks)
		) +
		geom_vline(xintercept = all_plot_data[[cohort]]$x_line_breaks, size = 0.25) +
		geom_hline(yintercept = all_plot_data[[cohort]]$y_breaks_cancer_major, size = 0.25) +
		geom_hline(yintercept = all_plot_data[[cohort]]$y_breaks_cancer_minor, size = 0.25, linetype = 'dashed') +
		theme(
			axis.ticks = element_blank(),
			axis.ticks.length = unit(0, 'pt'),
			panel.border = element_rect(colour = 'black', size = 0.5, fill = NA)
		) +
		labs(x = 'Chromosome', y = 'Cells', fill = 'Inferred CNA')
	
	annotations_cancer <- heat_map_labels_repel(
		all_plot_data[[cohort]]$patient_annotations_cancer,
		edge = 'right',
		nudge = 0.1,
		axis_title = 'Cancer cells',
		title_edge_margin = 5.5
	)
	
	cna_heatmap_caf_6padj <- ggplot(plot_data[classification == 'caf']) +
		geom_raster(
			aes(
				x = gene,
				y = factor(id, levels = all_plot_data[[cohort]]$classification_data[classification == 'caf'][order(patient, classification, id), unique(id)]),
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
			breaks = all_plot_data[[cohort]]$classification_data[classification == 'caf'][order(patient, classification, id), id[all_plot_data[[cohort]]$y_text_breaks_caf]],
			labels = all_plot_data[[cohort]]$classification_data[classification == 'caf'][
				order(patient, classification),
				sapply(all_plot_data[[cohort]]$y_text_breaks_caf, function(x) paste0('Patient ', patient[x]))
			]
		) +
		scale_x_discrete(
			expand = c(0, 0),
			breaks = plot_data[classification == 'caf', unique(gene)][all_plot_data[[cohort]]$x_text_breaks],
			labels = names(all_plot_data[[cohort]]$x_text_breaks)
		) +
		geom_vline(xintercept = all_plot_data[[cohort]]$x_line_breaks, size = 0.25) +
		geom_hline(yintercept = all_plot_data[[cohort]]$y_breaks_caf, size = 0.25) +
		theme(
			axis.ticks = element_blank(),
			axis.ticks.length = unit(0, 'pt'),
			panel.border = element_rect(colour = 'black', size = 0.5, fill = NA),
			plot.title = element_text(margin = margin(b = 30), size = 18)
		) +
		labs(x = 'Chromosome', y = 'Cells', fill = 'Inferred CNA', title = cohort_data[[cohort]]$plot_title)
	
	annotations_caf <- heat_map_labels_repel(
		all_plot_data[[cohort]]$patient_annotations_caf,
		edge = 'right',
		nudge = 0.1,
		axis_title = 'CAFs',
		title_edge_margin = 5.5
	)
	
	print(
		plot_grid(
			plot_grid(
				plot_grid(
					annotations_caf,
					cna_heatmap_caf_6padj + theme(
						axis.text = element_blank(),
						axis.title = element_blank(),
						legend.position = 'none',
						plot.margin = unit(c(5.5, 5.5, 0, 0), 'pt')
					),
					nrow = 1,
					ncol = 2,
					align = 'h',
					rel_widths = c(1, 4)
				),
				plot_grid(
					annotations_cancer,
					cna_heatmap_cancer_6padj + theme(
						axis.text.y = element_blank(),
						axis.title.y = element_blank(),
						legend.position = 'none',
						plot.margin = unit(c(40, 5.5, 5.5, 0), 'pt')
					),
					nrow = 1,
					ncol = 2,
					align = 'h',
					rel_widths = c(1, 4)
				),
				nrow = 2,
				ncol = 1,
				rel_heights = c(2, 3)
				# rel_heights = c(1, 2)
			),
			get_legend(cna_heatmap_caf),
			nrow = 1,
			ncol = 2,
			rel_widths = c(7, 1)
		)
	)
	
}

dev.off()
