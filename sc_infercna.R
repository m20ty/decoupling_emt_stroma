# bsub -q new-all.q -n 4 -R "rusage[mem=32000]" -o sc_infercna_log.o -e sc_infercna_log.e Rscript sc_infercna.R

cat(Sys.time(), '\n')

# The following is just to check the R version:
cat(R.Version()$version.string, '\n')

library(data.table) # 1.12.8
library(ggplot2) # This is version 3.3.1 on my laptop's WSL setup but 3.3.0 on WEXAC.
library(magrittr) # 1.5
library(infercna) # 1.0.0
library(stringr) # 1.4.0

source('general_functions.R')

# Check all the single cell datasets are there:

# all(
	# c(
		# 'chung_breast_cancer_2017.csv',
		# 'karaayvaz_tnbc_2018.csv',
		# 'qian_breast_2020.csv',
		# 'lee_crc_2020_kul3.csv',
		# 'lee_crc_2020_smc.csv',
		# 'qian_crc_2020.csv',
		# 'puram_hnscc_2017.csv',
		# 'ma_liver_2019.csv',
		# 'lambrechts_nsclc_2018.csv',
		# 'qian_lung_2020.csv',
		# 'song_nsclc_2019.csv',
		# 'kim_luad_2020.csv',
		# 'izar_ovarian_2020_ss2.csv',
		# 'izar_ovarian_2020_10x.csv',
		# 'qian_ovarian_2020.csv',
		# 'elyada_pdac_2019.csv',
		# 'peng_pdac_2019.csv'
	# ) %in% dir('../data_and_figures')
# )

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
	crc_lee_smc = list(
		read_quote = quote(fread('../data_and_figures/lee_crc_2020_smc.csv')[, -c('sample_id', 'cell_type', 'cell_subtype')]),
		ref_cell_clusters = list(
			b_cell = c(20, 25),
			endothelial = 19,
			macrophage = 10,
			t_cell = c(22, 23, 24)
		),
		ref_cell_clusters_normal = list(
			b_cell = 14,
			b_plasma = 10,
			endothelial= 4,
			epithelial = 1:2,
			fibroblast_myocyte = c(3, 5:7),
			macrophage = 9,
			t_cell = 11:13
		)
	),
	crc_qian = list(
		read_quote = quote(fread('../data_and_figures/qian_crc_2020.csv')[, -'cell_type']),
		ref_cell_clusters = list(
			b_cell = c(9, 13),
			b_plasma = c(8, 16),
			endothelial = 6,
			macrophage = 2,
			mast = 7,
			t_cell = c(3, 18)
		),
		ref_cell_clusters_normal = list(
			b_cell = 12,
			b_plasma = 2,
			endothelial = 5,
			epithelial = c(6, 13),
			fibroblast_myocyte = c(3, 7:9),
			macrophage = 1,
			myocyte = 11,
			t_cell = 4
		)
	),
	hnscc_puram = list(
		read_quote = quote(fread('../data_and_figures/puram_hnscc_2017.csv')[, -c('cell_type', 'lymph_node', 'processed_by_maxima_enzyme')]),
		ref_cell_clusters = list(
			b_cell = 4,
			endothelial = c(11, 14),
			macrophage_dc = 9,
			mast = 16,
			t_cell = c(12, 17)
		)
	),
	liver_ma = list(
		read_quote = quote(fread('../data_and_figures/ma_liver_2019.csv')[, -c('cell_type', 'sample', 'disease')]),
		ref_cell_clusters = list(
			b_cell = c(11, 18),
			endothelial = c(2, 14),
			macrophage = 5,
			t_cell = c(3, 6)
		)
	),
	liver_ma_hcc = list(
		read_quote = quote(fread('../data_and_figures/ma_liver_2019.csv')[disease == 'HCC', -c('cell_type', 'sample', 'disease')]),
		ref_cell_clusters = list(
			b_cell = 9,
			endothelial = c(2, 4, 5, 8),
			macrophage = 7,
			t_cell = 3
		)
	),
	liver_ma_icca = list(
		read_quote = quote(fread('../data_and_figures/ma_liver_2019.csv')[disease == 'iCCA', -c('cell_type', 'sample', 'disease')]),
		ref_cell_clusters = list(
			b_cell = 10,
			endothelial = c(1, 14),
			macrophage = c(3, 5, 11),
			t_cell = 6
		)
	),
	# The following dataset does have normal samples (saved separately), but I don't remember if they're matched, and anyway I added this
	# dataset to the analysis after doing all the others and deciding that there's no benefit to using matched normal as reference.
	luad_kim = list(
		read_quote = quote(fread('../data_and_figures/kim_luad_2020.csv')[, -c('cell_type', 'cell_type_refined', 'cell_subtype', 'sample_site', 'sample_id', 'sample_origin')]),
		ref_cell_clusters = list(
			b_cell = c(7, 8, 11), # 17 could also be B cells, but they also have T cell signal, so could be doublets
			b_plasma = 13,
			dc = 23,
			endothelial = 21,
			macrophage = c(1, 19),
			mast = 5,
			t_cell = 6 # As above, 17 could be T cells, and possibly also 22, but signal in 22 is weak
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
		
		cat('\tInferring CNAs using tumour cells as reference\n')
		
		sc_dbscan <- readRDS(paste0('../data_and_figures/dbscan_', cohort, '.rds'))
		
		ref_cells <- sapply(
			cohort_data[[cohort]]$ref_cell_clusters,
			function(x) sc_data[sample_type != 'normal'][sc_dbscan$cluster %in% x, id],
			simplify = FALSE,
			USE.NAMES = TRUE
		)
		
		# Take genes with log average TPM at least 4:
		
		cat('\t\tFiltering genes\n')
		
		gene_averages <- sapply(
			sc_data[sample_type != 'normal', -c('id', 'patient', 'sample_type')],
			function(x) {log2(mean(10*(2^x - 1)) + 1)},
			USE.NAMES = TRUE
		)
		
		sc_data_subset <- sc_data[sample_type != 'normal', c('id', 'patient', names(gene_averages)[gene_averages >= 4]), with = FALSE]
		
		# Running infercna() function:
		
		cat('\t\tInferring CNAs\n')
		
		inferred_cna <- sapply(
			as.character(unique(sc_data_subset$patient)),
			function(p) {
				infercna(
					t(
						set_rownames(
							as.matrix(sc_data_subset[patient == p | id %in% unlist(ref_cells), -c('id', 'patient')]),
							sc_data_subset[patient == p | id %in% unlist(ref_cells), id]
						)
					),
					refCells = ref_cells,
					isLog = TRUE
				)
			},
			simplify = FALSE,
			USE.NAMES = TRUE
		)
		
		cat('\t\tSaving CNA matrices\n')
		
		saveRDS(inferred_cna, paste0('../data_and_figures/inferred_cna_', cohort, '.rds'))
		
		cat('\t\tMaking plots\n')
		
		# The genes should be the same in all of the inferred_cna matrices:
		x_line_breaks <- .chromosomeBreaks(rownames(inferred_cna[[1]]))
		x_text_breaks <- round(.chromosomeBreaks(rownames(inferred_cna[[1]]), halfway = TRUE, hide = c('13', '18', '21', 'Y')))
		
		clust_plots <- sapply(
			names(inferred_cna),
			function(p) {
				
				cell_ids <- colnames(inferred_cna[[p]])[!(colnames(inferred_cna[[p]]) %in% unlist(ref_cells))]
				
				if(length(cell_ids) < 2) {
					warning('Not enough non-reference cells in this tumour!')
					return(NULL)
				}
				
				cell_clust <- hclust(dist(t(inferred_cna[[p]][, cell_ids])))
				
				cna_heatmap <- ggplot(
					reshape2::melt(
						inferred_cna[[p]][, cell_ids],
						varnames = c('gene', 'cell'),
						value.name = 'cna_score'
					)
				) +
					geom_raster(
						aes(
							x = factor(gene, levels = rownames(inferred_cna[[p]])),
							y = factor(cell, levels = cell_ids[cell_clust$order]),
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
						axis.ticks = element_blank(),
						axis.ticks.length = unit(0, 'pt'),
						panel.border = element_rect(colour = 'black', size = 0.5, fill = NA)
					) +
					labs(x = 'Chromosome', y = 'Cells', fill = 'Inferred CNA', title = p)
				
				# We don't need to return cell_ids, because these are stored in cell_clust$labels anyway.
				
				list(cell_clust = cell_clust, heatmap = cna_heatmap)
				
			},
			simplify = FALSE,
			USE.NAMES = TRUE
		)
		
		# for(x in clust_plots) saveRDS(x$cell_clust, paste0('../data_and_figures/cna_clust_', cohort, '.rds'))
		
		saveRDS(
			sapply(clust_plots, `[[`, 'cell_clust', simplify = FALSE, USE.NAMES = TRUE),
			paste0('../data_and_figures/cna_clust_', cohort, '.rds')
		)
		
		pdf(paste0('../data_and_figures/cna_plots_', cohort, '.pdf'))
		for(x in clust_plots) print(x$heatmap)
		dev.off()
		
		cat('\t\tDone!\n')
		
		# Next using the adjacent normal samples as reference:
		
		cat('\tInferring CNAs using normal cells as reference\n')
		
		sc_dbscan <- readRDS(paste0('../data_and_figures/dbscan_normal_', cohort, '.rds'))
		
		ref_cells <- sapply(
			cohort_data[[cohort]]$ref_cell_clusters_normal,
			function(x) sc_data[sample_type == 'normal'][sc_dbscan$cluster %in% x, id],
			simplify = FALSE,
			USE.NAMES = TRUE
		)
		
		# Take genes with log average TPM at least 4:
		
		cat('\t\tFiltering genes\n')
		
		gene_averages <- sapply(
			sc_data[, -c('id', 'patient', 'sample_type')],
			function(x) {log2(mean(10*(2^x - 1)) + 1)},
			USE.NAMES = TRUE
		)
		
		sc_data_subset <- sc_data[, c('id', 'patient', names(gene_averages)[gene_averages >= 4]), with = FALSE]
		
		# Running infercna() function:
		
		cat('\t\tInferring CNAs\n')
		
		inferred_cna <- sapply(
			as.character(unique(sc_data_subset$patient)),
			function(p) {
				infercna(
					t(
						set_rownames(
							as.matrix(sc_data_subset[patient == p | id %in% unlist(ref_cells), -c('id', 'patient')]),
							sc_data_subset[patient == p | id %in% unlist(ref_cells), id]
						)
					),
					refCells = ref_cells,
					isLog = TRUE
				)
			},
			simplify = FALSE,
			USE.NAMES = TRUE
		)
		
		cat('\t\tSaving CNA matrices\n')
		
		saveRDS(inferred_cna, paste0('../data_and_figures/inferred_cna_normalref_', cohort, '.rds'))
		
		cat('\t\tMaking plots\n')
		
		# The genes should be the same in all of the inferred_cna matrices:
		x_line_breaks <- .chromosomeBreaks(rownames(inferred_cna[[1]]))
		x_text_breaks <- round(.chromosomeBreaks(rownames(inferred_cna[[1]]), halfway = TRUE, hide = c('13', '18', '21', 'Y')))
		
		clust_plots <- sapply(
			names(inferred_cna),
			function(p) {
				
				cell_ids <- colnames(inferred_cna[[p]])[!(colnames(inferred_cna[[p]]) %in% unlist(ref_cells))]
				
				if(length(cell_ids) < 2) {
					warning('Not enough non-reference cells in this tumour!')
					return(NULL)
				}
				
				cell_clust <- hclust(dist(t(inferred_cna[[p]][, cell_ids])))
				
				cna_heatmap <- ggplot(
					reshape2::melt(
						inferred_cna[[p]][, cell_ids],
						varnames = c('gene', 'cell'),
						value.name = 'cna_score'
					)
				) +
					geom_raster(
						aes(
							x = factor(gene, levels = rownames(inferred_cna[[p]])),
							y = factor(cell, levels = cell_ids[cell_clust$order]),
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
						axis.ticks = element_blank(),
						axis.ticks.length = unit(0, 'pt'),
						panel.border = element_rect(colour = 'black', size = 0.5, fill = NA)
					) +
					labs(x = 'Chromosome', y = 'Cells', fill = 'Inferred CNA', title = p)
				
				# We don't need to return cell_ids, because these are stored in cell_clust$labels anyway.
				
				list(cell_clust = cell_clust, heatmap = cna_heatmap)
				
			},
			simplify = FALSE,
			USE.NAMES = TRUE
		)
		
		# for(x in clust_plots) saveRDS(x$cell_clust, paste0('../data_and_figures/cna_clust_normalref_', cohort, '.rds'))
		
		saveRDS(
			sapply(clust_plots, `[[`, 'cell_clust', simplify = FALSE, USE.NAMES = TRUE),
			paste0('../data_and_figures/cna_clust_normalref_', cohort, '.rds')
		)
		
		pdf(paste0('../data_and_figures/cna_plots_normalref_', cohort, '.pdf'))
		for(x in clust_plots) print(x$heatmap)
		dev.off()
		
		cat('\t\tDone!\n')
		
	} else { # Here we only use tumour cells as reference
		
		which_avail <- c(
			'sample_type' %in% colnames(sc_data),
			'ref_cell_clusters_normal' %in% names(cohort_data[[cohort]]),
			paste0('dbscan_normal_', cohort, '.rds') %in% dir('../data_and_figures')
		)
		
		if(any(which_avail)) {
			
			warning_message <- paste0(
				paste(
					c("column 'sample_type'", "data 'ref_cell_clusters_normal'", paste0("file 'dbscan_normal_", cohort, ".rds'"))[which_avail],
					collapse = ' and '
				),
				' found but ',
				paste(
					c("column 'sample_type'", "data 'ref_cell_clusters_normal'", paste0("file 'dbscan_normal_", cohort, ".rds'"))[!which_avail],
					collapse = ' and '
				),
				' unavailable. CNAs will be computed without using normal samples.'
			)
			
			warning_message <- gsub('^[a-z]', toupper(str_extract(warning_message, '^[a-z]')), warning_message)
			
			warning(warning_message)
			
		}
		
		sc_dbscan <- readRDS(paste0('../data_and_figures/dbscan_', cohort, '.rds'))
		
		ref_cells <- sapply(
			cohort_data[[cohort]]$ref_cell_clusters,
			function(x) sc_data$id[sc_dbscan$cluster %in% x],
			simplify = FALSE,
			USE.NAMES = TRUE
		)
		
		# Take genes with log average TPM at least 4:
		
		cat('\tFiltering genes\n')
		
		gene_averages <- sapply(
			sc_data[, -c('id', 'patient')],
			function(x) {log2(mean(10*(2^x - 1)) + 1)},
			USE.NAMES = TRUE
		)
		
		sc_data <- sc_data[, c('id', 'patient', names(gene_averages)[gene_averages >= 4]), with = FALSE]
		
		# Running infercna() function:
		
		cat('\tInferring CNAs\n')
		
		inferred_cna <- sapply(
			as.character(unique(sc_data$patient)),
			function(p) {
				infercna(
					t(
						set_rownames(
							as.matrix(sc_data[patient == p | id %in% unlist(ref_cells), -c('id', 'patient')]),
							sc_data[patient == p | id %in% unlist(ref_cells), id]
						)
					),
					refCells = ref_cells,
					isLog = TRUE
				)
			},
			simplify = FALSE,
			USE.NAMES = TRUE
		)
		
		cat('\tSaving CNA matrices\n')
		
		saveRDS(inferred_cna, paste0('../data_and_figures/inferred_cna_', cohort, '.rds'))
		
		cat('\tMaking plots\n')
		
		# The genes should be the same in all of the inferred_cna matrices:
		x_line_breaks <- .chromosomeBreaks(rownames(inferred_cna[[1]]))
		x_text_breaks <- round(.chromosomeBreaks(rownames(inferred_cna[[1]]), halfway = TRUE, hide = c('13', '18', '21', 'Y')))
		
		clust_plots <- sapply(
			names(inferred_cna),
			function(p) {
				
				cell_ids <- colnames(inferred_cna[[p]])[!(colnames(inferred_cna[[p]]) %in% unlist(ref_cells))]
				
				# The following was necessary for the liver HCC data from Ma et al., in which patient H18 had no non-reference cells
				# according to my classification by t-SNE.
				if(length(cell_ids) < 2) {
					warning('Not enough non-reference cells in this tumour!')
					return(NULL)
				}
				
				cell_clust <- hclust(dist(t(inferred_cna[[p]][, cell_ids])))
				
				cna_heatmap <- ggplot(
					reshape2::melt(
						inferred_cna[[p]][, cell_ids],
						varnames = c('gene', 'cell'),
						value.name = 'cna_score'
					)
				) +
					geom_raster(
						aes(
							x = factor(gene, levels = rownames(inferred_cna[[p]])),
							y = factor(cell, levels = cell_ids[cell_clust$order]),
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
						axis.ticks = element_blank(),
						axis.ticks.length = unit(0, 'pt'),
						panel.border = element_rect(colour = 'black', size = 0.5, fill = NA)
					) +
					labs(x = 'Chromosome', y = 'Cells', fill = 'Inferred CNA', title = p)
				
				# We don't need to return cell_ids, because these are stored in cell_clust$labels anyway.
				
				list(cell_clust = cell_clust, heatmap = cna_heatmap)
				
			},
			simplify = FALSE,
			USE.NAMES = TRUE
		)
		
		# for(x in clust_plots) saveRDS(x$cell_clust, paste0('../data_and_figures/cna_clust_', cohort, '.rds'))
		
		saveRDS(
			sapply(clust_plots, `[[`, 'cell_clust', simplify = FALSE, USE.NAMES = TRUE),
			paste0('../data_and_figures/cna_clust_', cohort, '.rds')
		)
		
		pdf(paste0('../data_and_figures/cna_plots_', cohort, '.pdf'))
		for(x in clust_plots) print(x$heatmap)
		dev.off()
		
		cat('\tDone!\n')
		
	}
	
}

# gene quantile, cell quantile, gene quantile for cells, cell quantile for genes
# - First two for cells to include in average profile

# Scatterplots: use density, and maybe log10 scale.
