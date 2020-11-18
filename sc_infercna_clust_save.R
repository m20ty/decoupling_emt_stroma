# bsub -q new-all.q -n 4 -R "rusage[mem=32000]" -o sc_infercna_clust_save_log.o -e sc_infercna_clust_save_log.e Rscript sc_infercna_clust_save.R

cat(Sys.time(), '\n')

# The following is just to check the R version:
cat(R.Version()$version.string, '\n')

library(data.table) # 1.12.8
library(ggplot2) # This is version 3.3.1 on my laptop's WSL setup but 3.3.0 on WEXAC.
library(magrittr) # 1.5
library(infercna) # 1.0.0
library(stringr) # 1.4.0

source('general_functions.R')





cohort_data <- list(
	breast_chung = list(
		read_quote = quote(fread('../data_and_figures/chung_breast_cancer_2017.csv')[, -c('cell_type', 'lymph_node')]),
		ref_cell_clusters = list(
			b_t_cells = c(5, 6, 7, 15),
			macrophage = 11
		)
	),
	breast_karaayvaz = list(
		read_quote = quote(fread('../data_and_figures/karaayvaz_tnbc_2018.csv')[, -'cell_type']),
		ref_cell_clusters = list(
			b_t_mast = c(10, 12),
			endothelial = 7,
			macrophage = 9
		)
	),
	breast_qian = list(
		read_quote = quote(fread('../data_and_figures/qian_breast_2020.csv')[, -'cell_type']),
		ref_cell_clusters = list(
			b_cell = 7,
			b_plasma = 6,
			endothelial = 5,
			macrophage = 2,
			mast = 8,
			t_cell = 1
		)
	),
	# Note the KUL3 dataset also has "border" samples, which I assume are part of the tumour and not adjacent normal.
	crc_lee_kul3 = list(
		read_quote = quote(fread('../data_and_figures/lee_crc_2020_kul3.csv')[, -c('sample_id', 'cell_type', 'cell_subtype')]),
		ref_cell_clusters = list(
			b_cell = c(15, 17),
			b_plasma = c(16, 18),
			endothelial = 10,
			macrophage = 13,
			mast = 19,
			t_cell = 14
		),
		ref_cell_clusters_normal = list(
			b_cell = 12,
			b_plasma = 11,
			endothelial = 4,
			epithelial = 1:2,
			fibroblast = 3,
			macrophage = 9,
			t_mast = 10
		)
	),
	lung_lambrechts = list(
		read_quote = quote(fread('../data_and_figures/lambrechts_nsclc_2018.csv')[, -c('cell_type', 'cluster', 'disease', 'annotation')]),
		ref_cell_clusters = list(
			b_cell = 10,
			b_plasma = c(7, 13),
			endothelial = 5,
			macrophage = c(3, 4, 14, 19),
			mast = 2,
			t_cell = c(1, 11, 23, 27)
		),
		ref_cell_clusters_normal = list(
			b_cell = 16,
			b_plasma = 8,
			dc = 23,
			endothelial = 11,
			epithelial = c(1, 4, 15, 19, 21),
			fibroblast = 3,
			macrophage = c(2, 6, 12, 17:18, 22),
			mast = 10,
			t_cell = c(9, 13:14, 20)
		)
	),
	lung_lambrechts_luad = list(
		read_quote = quote(fread('../data_and_figures/lambrechts_nsclc_2018.csv')[disease == 'LUAD', -c('cell_type', 'cluster', 'disease', 'annotation')]),
		ref_cell_clusters = list(
			b_cell = 2,
			b_plasma = c(7, 14),
			endothelial = 5,
			macrophage = c(1, 9),
			t_cell = c(10, 12, 16)
		),
		ref_cell_clusters_normal = list(
			b_cell = 8,
			endothelial = c(7, 9),
			epithelial = c(4, 6),
			fibroblast = 3,
			macrophage_dc = c(1, 5, 10),
			t_cell = 2
		)
	),
	# There are only 122 normal LUSC cells in the Lambrechts dataset, so we only use tumour cells as reference for this dataset.
	lung_lambrechts_lusc = list(
		read_quote = quote(fread('../data_and_figures/lambrechts_nsclc_2018.csv')[
			sample_type != 'normal' & disease == 'LUSC',
			-c('cell_type', 'cluster', 'disease', 'sample_type', 'annotation')
		]),
		ref_cell_clusters = list(
			b_cell = 6,
			macrophage_dc = c(3, 4, 9),
			mast = 2,
			t_cell = 1
		)
	),
	lung_qian = list(
		read_quote = quote(fread('../data_and_figures/qian_lung_2020.csv')[, -c('disease', 'cell_type')]),
		ref_cell_clusters = list(
			b_cell = c(14, 25),
			b_plasma = c(6, 7),
			endothelial = 11,
			macrophage = 4,
			mast = c(20, 32),
			t_cell = c(2, 8, 18, 21, 24, 28, 29)
		),
		ref_cell_clusters_normal = list(
			endothelial = 6,
			epithelial = c(3:4, 10:11, 14, 21),
			fibroblast = 1,
			macrophage_dc = c(2, 7, 13, 16, 18, 20),
			mast = 5,
			t_cell = c(8:9, 15, 17, 19)
		)
	),
	lung_qian_luad = list(
		read_quote = quote(fread('../data_and_figures/qian_lung_2020.csv')[disease == 'LUAD', -c('disease', 'cell_type')]),
		ref_cell_clusters = list(
			b_cell = c(12, 15, 24),
			b_plasma = c(7, 14, 25),
			endothelial = 4,
			macrophage_dc = c(1, 16),
			mast = 18,
			t_cell = c(6, 10, 11, 17, 20, 21)
		),
		ref_cell_clusters_normal = list(
			b_cell = 9,
			endothelial = c(6, 11, 18),
			epithelial = c(5, 8, 14, 16, 22),
			fibroblast = 3,
			macrophage_dc = c(1, 7, 10, 15, 20),
			mast = 12,
			t_cell = c(2, 4, 13, 17, 19, 21)
		)
	),
	lung_qian_lusc = list(
		read_quote = quote(fread('../data_and_figures/qian_lung_2020.csv')[disease == 'LUSC', -c('disease', 'cell_type')]),
		ref_cell_clusters = list(
			b_cell = 15,
			b_plasma = c(9, 10),
			dc = 8,
			endothelial = 13,
			macrophage = c(4, 5, 12, 16),
			mast = c(3, 19),
			t_cell = c(2, 17)
		),
		ref_cell_clusters_normal = list(
			b_endothelial_fibroblast_mast_t = 3,
			epithelial = c(1, 4, 6, 8),
			macrophage_dc = c(2, 5, 7)
		)
	),
	lung_song = list(
		read_quote = quote(fread('../data_and_figures/song_nsclc_2019.csv')[, -c('cell_type', 'disease')]),
		ref_cell_clusters = list(
			b_cell = 15,
			b_t_cell = 3,
			macrophage_dc = c(2, 6, 11, 12, 13)
		),
		ref_cell_clusters_normal = list( # Epithelial cells seem quite patient-specific, so I'll leave them out
			b_plasma = 8,
			endothelial_fibroblast_myocyte = 20,
			macrophage_dc = c(3:6, 13, 15, 19),
			t_cell = c(2, 11)
		)
	),
	lung_song_luad = list(
		read_quote = quote(fread('../data_and_figures/song_nsclc_2019.csv')[disease == 'LUAD', -c('cell_type', 'disease')]),
		ref_cell_clusters = list(
			b_cell = 11,
			b_t_cell = 7,
			macrophage_dc = c(2, 8, 9)
		),
		ref_cell_clusters_normal = list(
			endothelial = 11,
			epithelial = c(8, 10, 12),
			macrophage_dc = c(2:4, 6:7, 9),
			t_cell = 1
		)
	),
	lung_song_lusc = list(
		read_quote = quote(fread('../data_and_figures/song_nsclc_2019.csv')[disease == 'LUSC', -c('cell_type', 'disease')]),
		ref_cell_clusters = list(
			endothelial = 7,
			macrophage_dc = c(1, 4)
		),
		ref_cell_clusters_normal = list( # Not sure this is informative enough
			b_plasma = 4,
			macrophage_t = 1
		)
	),
	ovarian_izar_10x = list(
		read_quote = quote(fread('../data_and_figures/izar_ovarian_2020_10x.csv')[, -c('barcode', 'sample_id', 'time', 'cell_type', 'cluster')]),
		ref_cell_clusters = list(
			b_cell = c(5, 8, 14),
			dc = 17,
			macrophage = c(1, 9, 12, 16, 19),
			t_cell = c(11, 13)
		)
	),
	ovarian_izar_ss2 = list(
		read_quote = quote(fread('../data_and_figures/izar_ovarian_2020_ss2.csv')[, -c('cell_type', 'cluster')]),
		ref_cell_clusters = list( # This was really not successful - I can only identify macrophages, and there aren't many of them.
			macrophage = 3
		)
	),
	ovarian_qian = list(
		read_quote = quote(fread('../data_and_figures/qian_ovarian_2020.csv')[, -c('site', 'cell_type')]),
		ref_cell_clusters = list(
			b_plasma = 7,
			dc_mast_t = c(5, 19),
			endothelial = 3,
			macrophage = c(2, 13, 18)
		),
		# In the following, I'm not sure about the epithelial and fibroblast clusters, because they're both from the same patient, and there's a somewhat
		# intermediate cluster between them which makes them together look like EMT.
		ref_cell_clusters_normal = list(
			endothelial = 7,
			epithelial = 2,
			fibroblast = c(1, 8),
			macrophage = c(3, 5),
			t_mast_b_plasma = 6
		)
	),
	pdac_elyada = list(
		read_quote = quote(fread('../data_and_figures/elyada_pdac_2019.csv')[, -c('cell_type', 'cluster')]),
		ref_cell_clusters = list(
			b_cell = c(9, 12),
			b_plasma = 11,
			endothelial = 14,
			macrophage = c(2, 3, 15),
			t_cell = 4
		)
	),
	pdac_peng = list(
		read_quote = quote(fread('../data_and_figures/peng_pdac_2019.csv')[, -'cell_type']),
		ref_cell_clusters = list(
			b_cell = c(8, 12),
			b_plasma = 6,
			endothelial = 3,
			macrophage_mast = c(2, 14),
			t_cell = 4
		),
		# The following wasn't very successful, because there's no clear B or T cell signal, and I don't really trust the clusters that look epithelial,
		# because they also express various other cell type markers (there could also be e.g. acinar cells), and they're quite patient-specific.
		ref_cell_clusters_normal = list(
			endothelial = c(2, 4, 9, 11),
			fibroblast = 5,
			macrophage = 7
		)
	)
)





for(cohort in names(cohort_data)) {
	
	cat(paste0(cohort, ':\n'))
	
	cat('\tReading in data\n')
	
	sc_data <- eval(cohort_data[[cohort]]$read_quote)
	
	if(
		'sample_type' %in% colnames(sc_data) &
			'ref_cell_clusters_normal' %in% names(cohort_data[[cohort]]) &
			paste0('dbscan_normal_', cohort, '.rds') %in% dir('../data_and_figures')
	) { # Use each of tumour and normal cells as reference
		
		# First, using cells from tumour sample as reference:
		
		cat('\tInferred CNAs with tumour cells as reference\n')
		
		sc_dbscan <- readRDS(paste0('../data_and_figures/dbscan_', cohort, '.rds'))
		
		ref_cells <- sapply(
			cohort_data[[cohort]]$ref_cell_clusters,
			function(x) sc_data[sample_type != 'normal'][sc_dbscan$cluster %in% x, id],
			simplify = FALSE,
			USE.NAMES = TRUE
		)
		
		cat('\t\tReading in inferred CNA matrices\n')
		
		inferred_cna <- readRDS(paste0('../data_and_figures/inferred_cna_', cohort, '.rds'))
		
		cat('\t\tRunning clustering and saving hclust objects\n')
		
		clust_objects <- sapply(
			names(inferred_cna),
			function(p) {
				
				cell_ids <- colnames(inferred_cna[[p]])[!(colnames(inferred_cna[[p]]) %in% unlist(ref_cells))]
				
				if(length(cell_ids) < 2) {
					warning('Not enough non-reference cells in this tumour!')
					return(NULL)
				}
				
				hclust(dist(t(inferred_cna[[p]][, cell_ids])))
				
			},
			simplify = FALSE,
			USE.NAMES = TRUE
		)
		
		saveRDS(clust_objects, paste0('../data_and_figures/cna_clust_', cohort, '.rds'))
		
		cat('\t\tDone!\n')
		
		# Next using the adjacent normal samples as reference:
		
		cat('\tInferred CNAs with normal cells as reference\n')
		
		sc_dbscan <- readRDS(paste0('../data_and_figures/dbscan_normal_', cohort, '.rds'))
		
		ref_cells <- sapply(
			cohort_data[[cohort]]$ref_cell_clusters_normal,
			function(x) sc_data[sample_type == 'normal'][sc_dbscan$cluster %in% x, id],
			simplify = FALSE,
			USE.NAMES = TRUE
		)
		
		cat('\t\tReading in inferred CNA matrices\n')
		
		inferred_cna <- readRDS(paste0('../data_and_figures/inferred_cna_normalref_', cohort, '.rds'))
		
		cat('\t\tRunning clustering and saving hclust objects\n')
		
		clust_objects <- sapply(
			names(inferred_cna),
			function(p) {
				
				cell_ids <- colnames(inferred_cna[[p]])[!(colnames(inferred_cna[[p]]) %in% unlist(ref_cells))]
				
				if(length(cell_ids) < 2) {
					warning('Not enough non-reference cells in this tumour!')
					return(NULL)
				}
				
				hclust(dist(t(inferred_cna[[p]][, cell_ids])))
				
			},
			simplify = FALSE,
			USE.NAMES = TRUE
		)
		
		saveRDS(clust_objects, paste0('../data_and_figures/cna_clust_normalref_', cohort, '.rds'))
		
		cat('\t\tDone!\n')
		
	} else { # Here we only use tumour cells as reference
		
		sc_dbscan <- readRDS(paste0('../data_and_figures/dbscan_', cohort, '.rds'))
		
		ref_cells <- sapply(
			cohort_data[[cohort]]$ref_cell_clusters,
			function(x) sc_data$id[sc_dbscan$cluster %in% x],
			simplify = FALSE,
			USE.NAMES = TRUE
		)
		
		cat('\t\tReading in inferred CNA matrices\n')
		
		inferred_cna <- readRDS(paste0('../data_and_figures/inferred_cna_', cohort, '.rds'))
		
		cat('\t\tRunning clustering and saving hclust objects\n')
		
		clust_objects <- sapply(
			names(inferred_cna),
			function(p) {
				
				cell_ids <- colnames(inferred_cna[[p]])[!(colnames(inferred_cna[[p]]) %in% unlist(ref_cells))]
				
				if(length(cell_ids) < 2) {
					warning('Not enough non-reference cells in this tumour!')
					return(NULL)
				}
				
				hclust(dist(t(inferred_cna[[p]][, cell_ids])))
				
			},
			simplify = FALSE,
			USE.NAMES = TRUE
		)
		
		saveRDS(clust_objects, paste0('../data_and_figures/cna_clust_', cohort, '.rds'))
		
		cat('\tDone!\n')
		
	}
	
}
