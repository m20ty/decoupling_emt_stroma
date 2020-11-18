# bsub -q new-short -R "rusage[mem=32000]" -o sc_cna_cor_signal_log.o -e sc_cna_cor_signal_log.e Rscript sc_cna_cor_signal.R

library(data.table)
library(ggplot2)
library(magrittr)
library(infercna)

source('general_functions.R')





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
	
	inferred_cna <- readRDS(paste0('../data_and_figures/inferred_cna_', cohort, '.rds'))
	
	pdf(paste0('../data_and_figures/cna_cor_signal_scatterplot_', cohort, '.pdf'))
	
	for(p in names(inferred_cna)) {
		
		scatterplot <- ggplot(
			data.table(
				cell_id = colnames(inferred_cna[[p]]),
				cna_signal = cnaSignal(inferred_cna[[p]], gene.quantile = 0.9),
				cna_cor = cnaCor(inferred_cna[[p]], cell.quantile = 0.9, gene.quantile = 0.9)
			),
			aes(x = cna_cor, y = cna_signal)
		) +
		geom_point(alpha = 0.25) +
		theme_minimal() +
		labs(title = p)
		
		print(scatterplot)
		
	}
	
	dev.off()
	
	if(paste0('inferred_cna_normalref_', cohort, '.rds') %in% dir('../data_and_figures')) {
		
		cat(paste0(cohort, ' - normal cell reference\n'))
		
		inferred_cna <- readRDS(paste0('../data_and_figures/inferred_cna_normalref_', cohort, '.rds'))
		
		pdf(paste0('../data_and_figures/cna_cor_signal_scatterplot_normalref_', cohort, '.pdf'))
		
		for(p in names(inferred_cna)) {
			
			scatterplot <- ggplot(
				data.table(
					cell_id = colnames(inferred_cna[[p]]),
					cna_signal = cnaSignal(inferred_cna[[p]], gene.quantile = 0.9),
					cna_cor = cnaCor(inferred_cna[[p]], cell.quantile = 0.9, gene.quantile = 0.9)
				),
				aes(x = cna_cor, y = cna_signal)
			) +
			geom_point(alpha = 0.25) +
			theme_minimal() +
			labs(title = p)
			
			print(scatterplot)
			
		}
		
		dev.off()
		
	}
	
}





# I decided on the gene.quantile and cell.quantile values by looking at the case of HNSCC, as below.  I also tried different values for gene.quantile.for.cells
# and cell.quantile.for.genes, but this didn't seem to make a lot of difference.

# inferred_cna <- readRDS('../data_and_figures/inferred_cna_hnscc_puram.rds')

# for(p in names(inferred_cna)) {
	
	# cna_cor_signal_data <- rbindlist(
		# lapply(
			# c(0.25, 0.5, 0.75, 0.9),
			# function(x) {
				# lapply(
					# c(0.25, 0.5, 0.75, 0.9),
					# function(y) {
						# data.table(
							# cell_id = colnames(inferred_cna[[p]]),
							# # gene_quantile_for_cells = x,
							# # cell_quantile_for_genes = y,
							# # cna_signal = colMeans(inferred_cna[[p]]^2)
							# # cna_signal = cnaSignal(inferred_cna[[p]], gene.quantile = 0.9, cell.quantile.for.genes = y),
							# # cna_cor = cnaCor(inferred_cna[[p]], cell.quantile = 0.9, gene.quantile.for.cells = x, gene.quantile = 0.9, cell.quantile.for.genes = y)
							# gene_quantile = x,
							# cell_quantile = y,
							# cna_signal = cnaSignal(inferred_cna[[p]], gene.quantile = x),
							# cna_cor = cnaCor(inferred_cna[[p]], cell.quantile = y, gene.quantile = x)
						# )
					# }
				# )
			# }
		# ) %>% unlist(recursive = FALSE)
	# )
	
	# facet_scatterplot <- ggplot(cna_cor_signal_data, aes(x = cna_cor, y = cna_signal)) +
		# geom_point(alpha = 0.25) +
		# # facet_grid(rows = vars(gene_quantile_for_cells), cols = vars(cell_quantile_for_genes))
		# facet_grid(
			# rows = vars(gene_quantile),
			# cols = vars(cell_quantile),
			# labeller = labeller(
				# .cols = c(`0.25` = 'cell quantile: 0.25', `0.5` = 'cell quantile: 0.5', `0.75` = 'cell quantile: 0.75', `0.9` = 'cell quantile: 0.9'),
				# .rows = c(`0.25` = 'gene quantile: 0.25', `0.5` = 'gene quantile: 0.5', `0.75` = 'gene quantile: 0.75', `0.9` = 'gene quantile: 0.9')
			# )
		# ) +
		# theme_minimal()
	
	# print(facet_scatterplot)
	
# }
