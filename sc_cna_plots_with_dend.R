library(data.table) # 1.12.8
library(ggplot2) # 3.3.0
library(magrittr) # 1.5
library(infercna) # 1.0.0
library(cowplot) # 1.0.0

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

cohorts <- c(
	'breast_chung',
	'breast_karaayvaz',
	'breast_qian',
	'crc_lee_kul3',
	'crc_lee_smc',
	'crc_qian',
	'hnscc_puram',
	'liver_ma',
	'liver_ma_hcc',
	'liver_ma_icca',
	'lung_lambrechts',
	'lung_lambrechts_luad',
	'lung_lambrechts_lusc',
	'lung_qian',
	'lung_qian_luad',
	'lung_qian_lusc',
	'lung_song',
	'lung_song_luad',
	'lung_song_lusc',
	'ovarian_izar_10x',
	'ovarian_izar_ss2',
	'ovarian_qian',
	'pdac_elyada',
	'pdac_peng'
)





for(cohort in cohorts) {

	cat(paste0(cohort, '\n'))

	cat('\tReading in data\n')

	inferred_cna <- readRDS(paste0('../data_and_figures/inferred_cna_', cohort, '.rds'))
	cell_clust <- readRDS(paste0('../data_and_figures/cna_clust_', cohort, '.rds'))

	# The genes should be the same in all of the inferred_cna matrices:
	x_line_breaks <- .chromosomeBreaks(rownames(inferred_cna[[1]]))
	x_text_breaks <- round(.chromosomeBreaks(rownames(inferred_cna[[1]]), halfway = TRUE, hide = c('13', '18', '21', 'Y')))

	cat('\tMaking plots\n')

	pdf(paste0('../data_and_figures/cna_htmp_dend_', cohort, '.pdf'), width = 10, height = 8)

	for(p in names(inferred_cna)) {

		if(!is.null(cell_clust[[p]])) {

			cna_dendro <- dendro(cell_clust[[p]], edge = 'right')

			cna_heatmap <- ggplot(
				reshape2::melt(inferred_cna[[p]][, cell_clust[[p]]$labels], varnames = c('gene', 'cell'), value.name = 'cna_score')
			) +
				geom_raster(
					aes(
						x = factor(gene, levels = rownames(inferred_cna[[p]])),
						y = factor(cell, levels = cell_clust[[p]]$labels[cell_clust[[p]]$order]),
						fill = cna_score
					)
				) +
				scale_fill_gradientn(
					colours = rev(colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))(50)),
					limits = c(-1, 1),
					oob = scales::squish
				) +
				scale_y_discrete(expand = c(0, 0)) +
				scale_x_discrete(expand = c(0, 0), breaks = rownames(inferred_cna[[p]])[x_text_breaks], labels = names(x_text_breaks)) +
				geom_vline(xintercept = x_line_breaks, size = 0.25) +
				theme(
					axis.text.y = element_blank(),
					axis.title.y = element_blank(),
					axis.ticks = element_blank(),
					axis.ticks.length = unit(0, 'pt'),
					panel.border = element_rect(colour = 'black', size = 0.5, fill = NA),
					plot.margin = unit(c(5.5, 5.5, 5.5, 0), 'pt')
				) +
				labs(x = 'Chromosome', fill = 'Inferred CNA', title = p)

			combined_plot <- plot_grid(cna_dendro, cna_heatmap, nrow = 1, ncol = 2, align = 'h', rel_widths = c(1, 5))

			print(combined_plot)

		}

	}

	dev.off()

	cat('\tDone!\n')

	if(paste0('inferred_cna_normalref_', cohort, '.rds') %in% dir('../data_and_figures')) {

		cat(paste0(cohort, ' - normal cell reference\n'))

		cat('\tReading in data\n')

		inferred_cna <- readRDS(paste0('../data_and_figures/inferred_cna_normalref_', cohort, '.rds'))
		cell_clust <- readRDS(paste0('../data_and_figures/cna_clust_normalref_', cohort, '.rds'))

		# The genes should be the same in all of the inferred_cna matrices:
		x_line_breaks <- .chromosomeBreaks(rownames(inferred_cna[[1]]))
		x_text_breaks <- round(.chromosomeBreaks(rownames(inferred_cna[[1]]), halfway = TRUE, hide = c('13', '18', '21', 'Y')))

		cat('\tMaking plots\n')

		pdf(paste0('../data_and_figures/cna_htmp_dend_normalref_', cohort, '.pdf'), width = 10, height = 8)

		for(p in names(inferred_cna)) {

			if(!is.null(cell_clust[[p]])) {

				cna_dendro <- dendro(cell_clust[[p]], edge = 'right')

				cna_heatmap <- ggplot(
					reshape2::melt(inferred_cna[[p]][, cell_clust[[p]]$labels], varnames = c('gene', 'cell'), value.name = 'cna_score')
				) +
					geom_raster(
						aes(
							x = factor(gene, levels = rownames(inferred_cna[[p]])),
							y = factor(cell, levels = cell_clust[[p]]$labels[cell_clust[[p]]$order]),
							fill = cna_score
						)
					) +
					scale_fill_gradientn(
						colours = rev(colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))(50)),
						limits = c(-1, 1),
						oob = scales::squish
					) +
					scale_y_discrete(expand = c(0, 0)) +
					scale_x_discrete(expand = c(0, 0), breaks = rownames(inferred_cna[[p]])[x_text_breaks], labels = names(x_text_breaks)) +
					geom_vline(xintercept = x_line_breaks, size = 0.25) +
					theme(
						axis.text.y = element_blank(),
						axis.title.y = element_blank(),
						axis.ticks = element_blank(),
						axis.ticks.length = unit(0, 'pt'),
						panel.border = element_rect(colour = 'black', size = 0.5, fill = NA),
						plot.margin = unit(c(5.5, 5.5, 5.5, 0), 'pt')
					) +
					labs(x = 'Chromosome', fill = 'Inferred CNA', title = p)

				combined_plot <- plot_grid(cna_dendro, cna_heatmap, nrow = 1, ncol = 2, align = 'h', rel_widths = c(1, 5))

				print(combined_plot)

			}

		}

		dev.off()

		cat('\tDone!\n')

	}

}
