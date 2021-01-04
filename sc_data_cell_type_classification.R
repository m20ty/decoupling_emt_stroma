library(data.table) # 1.12.8
library(ggplot2) # 3.3.0
library(magrittr) # 1.5
library(plyr) # 1.8.6
library(msigdbr) # 7.1.1
library(clusterProfiler) # 3.14.3

source('general_functions.R')

# Get MSigDB C2, C5 and H gene sets (C2 = Curated gene sets; C5 = GO gene sets):
msigdb <- rbindlist(lapply(c('C2', 'C5', 'H'), function(categ) as.data.table(msigdbr(category = categ))))

# Function to test one cluster against all others, find DE genes and run GSEA:
test_cluster <- function(
	expression_mat,
	clust_data,
	clust,
	msigdb_table,
	clust_var_name = 'cluster',
	cell_id_var_name = 'cell_id',
	seed = NULL
) {

	# Centre the matrix gene-wise:
	expression_mat_centred <- t(apply(expression_mat, 2, function(x) {x - mean(x)}))

	averages <- sapply(
		clust_data[get(clust_var_name) != 'noise', unique(get(clust_var_name))],
		function(clst) rowMeans(expression_mat_centred[, clust_data[get(clust_var_name) == clst, get(cell_id_var_name)]]),
		simplify = FALSE,
		USE.NAMES = TRUE
	)

	gene_test <- lapply(
		lapply(clust_data[get(clust_var_name) != clust, unique(get(clust_var_name))], function(ct) c(clust, ct)),
		function(clustpair) {

			cat(paste0(clustpair[1], ' vs. ', clustpair[2], '...'))

			mat_clustpair <- sapply(
				clustpair,
				function(clst) expression_mat_centred[, clust_data[get(clust_var_name) == clst, get(cell_id_var_name)]],
				simplify = FALSE,
				USE.NAMES = TRUE
			)

			out <- list(
				clust_pair = clustpair,
				test_results = data.table(gene_id = rownames(expression_mat_centred))[
					,
					.(pval = try_default(t.test(mat_clustpair[[1]][gene_id, ], mat_clustpair[[2]][gene_id, ])$p.value, default = 1, quiet = TRUE)),
					by = gene_id
				]
			)

			cat('Done!\n')

			out

		}
	)

	gene_test_results <- lapply(
		gene_test,
		function(x) {
			data <- copy(x$test_results)
			data[, rel_exp := averages[[x$clust_pair[1]]][gene_id] - averages[[x$clust_pair[2]]][gene_id]]
			list(other_clst = x$clust_pair[2], data = data)
		}
	)

	gene_test_results <- setNames(lapply(gene_test_results, `[[`, 'data'), sapply(gene_test_results, `[[`, 'other_clst'))

	table_list <- lapply(gene_test_results, function(dt) dt[p.adjust(pval, method ='BH') < 0.05])
	for(dt in table_list) {setkey(dt, gene_id)}

	degenes <- data.table(gene = Reduce(intersect, lapply(table_list, `[[`, 'gene_id')))
	degenes[, rel_exp := apply(as.data.table(lapply(table_list, function(dt) dt[gene, rel_exp])), 1, mean)]

	# Run GSEA:

	if(!is.null(seed)) set.seed(seed)

	degenes_gsea <- sapply(
		unique(msigdb_table$gs_cat),
		function(categ) {
			cbind(
				gs_cat = categ,
				as.data.table(
					GSEA(degenes[order(-rel_exp), setNames(rel_exp, gene)], TERM2GENE = msigdb_table[gs_cat == categ, .(gs_name, human_gene_symbol)])
				)
			)
		},
		simplify = FALSE,
		USE.NAMES = TRUE
	)

	gsea_table <- rbindlist(lapply(degenes_gsea, function(dt) dt[, Reduce(intersect, lapply(degenes_gsea, names)), with = FALSE]))

	list(test_results = gene_test_results, degenes = degenes, gsea_table = gsea_table)

}





# breast_qian:

classification_data <- fread(paste0('../data_and_figures/sc_reclassify/classification_data_breast_qian.csv'), key = 'cell_id')
sc_data <- fread('../data_and_figures/qian_breast_2020.csv')[patient %in% c(42, 43, 47, 49, 51, 53, 54)]
sc_data[, classification := classification_data[id, classification]]
sc_data <- sc_data[classification == 'nonmalignant', -'classification']
gene_averages <- sapply(sc_data[, -c('id', 'patient', 'cell_type')], function(x) {log2(mean(10*(2^x - 1)) + 1)}, USE.NAMES = TRUE)
sc_data <- sc_data[, c('id', 'patient', 'cell_type', names(gene_averages)[gene_averages >= 4]), with = FALSE]
sc_tsne <- readRDS(paste0('../data_and_figures/sc_reclassify/tsne_breast_qian.rds'))
sc_dbscan <- readRDS(paste0('../data_and_figures/sc_reclassify/dbscan_breast_qian.rds'))

plot_data <- setNames(
	cbind(as.data.table(sc_tsne$Y), as.character(sc_dbscan$cluster), sc_data$cell_type),
	c('x', 'y', 'dbscan_cluster', 'cell_type')
)

plot_data[dbscan_cluster == 0, dbscan_cluster := 'noise']

cell_type_data <- data.table(
	cell_id = sc_data$id,
	cluster = mapvalues(
		sc_dbscan$cluster,
		c(0, 8, 4, 6, 3, 9, 1, 7, 10, 2, 5),
		c('noise', 'b_cell', 'b_plasma', 'endothelial', 'macrophage', 'mast', 't_cell', 't_cell', 't_cell', 'caf_1', 'caf_2')
	)
)

setkey(sc_data, id)

# Try to idenfity cluster 11:

de_11 <- test_cluster(
	set_rownames(
		as.matrix(sc_data[cell_type_data[cluster != 'noise', cell_id], -c('id', 'patient', 'cell_type')]),
		cell_type_data[cluster != 'noise', cell_id]
	),
	cell_type_data[cluster != 'noise'],
	'11',
	msigdb,
	seed = 9694
)

de_11_top_50 <- de_11$degenes[order(-rel_exp), head(gene, 50)]

# Check for differences between the two CAF clusters:

caf_seeds <- c(4843, 6826, 6365, 878)

# Difference between each CAF cluster and all the other clusters:
de_caf <- lapply(
	1:2,
	function(i) {
		test_cluster(
			set_rownames(
				as.matrix(
                    sc_data[cell_type_data[!(cluster %in% c('noise', paste0('caf_', (1:2)[1:2 != i]))), cell_id], -c('id', 'patient', 'cell_type')]
                ),
				cell_type_data[!(cluster %in% c('noise', paste0('caf_', (1:2)[1:2 != i]))), cell_id]
			),
			cell_type_data[!(cluster %in% c('noise', paste0('caf_', (1:2)[1:2 != i])))],
			paste0('caf_', i),
			msigdb,
			seed = caf_seeds[i]
		)
	}
)

de_caf_top_50 <- set_colnames(sapply(1:2, function(i) de_caf[[i]]$degenes[order(-rel_exp), head(gene, 50)]), c('Cluster 2', 'Cluster 5'))

# Difference between the CAF clusters:
de_between_caf <- lapply(
	1:2,
	function(i) {
		test_cluster(
			set_rownames(
				as.matrix(sc_data[cell_type_data[startsWith(cluster, 'caf'), cell_id], -c('id', 'patient', 'cell_type')]),
				cell_type_data[startsWith(cluster, 'caf'), cell_id]
			),
			cell_type_data[startsWith(cluster, 'caf')],
			paste0('caf_', i),
			msigdb,
			seed = caf_seeds[i + 2]
		)
	}
)

de_between_caf_top_50 <- set_colnames(
    sapply(1:2, function(i) de_between_caf[[i]]$degenes[order(-rel_exp), head(gene, 50)]),
    c('Cluster 2', 'Cluster 5')
)

# Save stuff:
fwrite(plot_data, '../data_and_figures/sc_reclassify/files_for_rmd/breast_qian_plot_data.csv')
saveRDS(de_11_top_50, '../data_and_figures/sc_reclassify/files_for_rmd/breast_qian_de_11_top_50.rds')
saveRDS(de_caf_top_50, '../data_and_figures/sc_reclassify/files_for_rmd/breast_qian_de_caf_top_50.rds')
saveRDS(de_between_caf_top_50, '../data_and_figures/sc_reclassify/files_for_rmd/breast_qian_de_between_caf_top_50.rds')





# crc_lee_smc:

classification_data <- fread(paste0('../data_and_figures/sc_reclassify/classification_data_crc_lee_smc.csv'), key = 'cell_id')
sc_data <- fread('../data_and_figures/lee_crc_2020_smc.csv')[
	sample_type != 'normal',
	-c('sample_id', 'sample_type', 'cell_subtype')
][
	patient %in% c('SMC01', 'SMC02', 'SMC04', 'SMC07', 'SMC08', 'SMC09', 'SMC11', 'SMC14', 'SMC15', 'SMC16', 'SMC18', 'SMC20', 'SMC21',
				   'SMC23', 'SMC25')
]
sc_data[, classification := classification_data[id, classification]]
sc_data <- sc_data[classification == 'nonmalignant', -'classification']
gene_averages <- sapply(sc_data[, -c('id', 'patient', 'cell_type')], function(x) {log2(mean(10*(2^x - 1)) + 1)}, USE.NAMES = TRUE)
sc_data <- sc_data[, c('id', 'patient', 'cell_type', names(gene_averages)[gene_averages >= 4]), with = FALSE]
sc_tsne <- readRDS(paste0('../data_and_figures/sc_reclassify/tsne_crc_lee_smc.rds'))
sc_dbscan <- readRDS(paste0('../data_and_figures/sc_reclassify/dbscan_crc_lee_smc.rds'))

plot_data <- setNames(
	cbind(as.data.table(sc_tsne$Y), as.character(sc_dbscan$cluster), sc_data$cell_type),
	c('x', 'y', 'dbscan_cluster', 'cell_type')
)

plot_data[dbscan_cluster == 0, dbscan_cluster := 'noise']

# Initial cell type definitions:
cell_type_data <- data.table(
	cell_id = sc_data$id,
	cluster = mapvalues(
		sc_dbscan$cluster,
		c(0, 9, 8, 2, 3, 1, 4, 6, 7),
		c('noise', 'b_cell', 'b_plasma', 'caf', 'endothelial', 'epithelial', 'macrophage', 't_cell', 't_cell')
	)
)

cell_type_data[cluster == 'noise' & sc_tsne$Y[, 2] > 15 & sc_tsne$Y[, 2] < 25, cluster := 'noise_1']
cell_type_data[cluster == 'noise' & sc_tsne$Y[, 2] > 40, cluster := 'noise_2']

setkey(sc_data, id)

# Try to idenfity cluster 5:

de_5 <- test_cluster(
	set_rownames(
		as.matrix(sc_data[cell_type_data[!startsWith(cluster, 'noise'), cell_id], -c('id', 'patient', 'cell_type')]),
		cell_type_data[!startsWith(cluster, 'noise'), cell_id]
	),
	cell_type_data[!startsWith(cluster, 'noise')],
	'5',
	msigdb,
	seed = 9627
)

de_5_top_50 <- de_5$degenes[order(-rel_exp), head(gene, 50)]

# Check if endothelial cluster really is distinct from the CAF cluster:

caf_endothelial_seeds <- c(8702, 1256, 2881, 9041)

de_caf_endothelial <- lapply(
	1:2,
	function(i) {
		test_cluster(
			set_rownames(
				as.matrix(sc_data[cell_type_data[!startsWith(cluster, 'noise'), cell_id], -c('id', 'patient', 'cell_type')]),
				cell_type_data[!startsWith(cluster, 'noise'), cell_id]
			),
			cell_type_data[!startsWith(cluster, 'noise')],
			c('caf', 'endothelial')[i],
			msigdb,
			seed = caf_endothelial_seeds[i]
		)
	}
)

de_caf_endothelial_top_50 <- set_colnames(
	sapply(1:2, function(i) de_caf_endothelial[[i]]$degenes[order(-rel_exp), head(gene, 50)]),
	c('CAF', 'endothelial')
)

de_between_caf_endothelial <- lapply(
	1:2,
	function(i) {
		test_cluster(
			set_rownames(
				as.matrix(sc_data[cell_type_data[cluster %in% c('caf', 'endothelial'), cell_id], -c('id', 'patient', 'cell_type')]),
				cell_type_data[cluster %in% c('caf', 'endothelial'), cell_id]
			),
			cell_type_data[cluster %in% c('caf', 'endothelial')],
			c('caf', 'endothelial')[i],
			msigdb,
			seed = caf_endothelial_seeds[i + 2]
		)
	}
)

de_between_caf_endothelial_top_50 <- set_colnames(
	sapply(1:2, function(i) de_between_caf_endothelial[[i]]$degenes[order(-rel_exp), head(gene, 50)]),
	c('CAF', 'endothelial')
)

de_noise_2 <- test_cluster(
	set_rownames(
		as.matrix(sc_data[cell_type_data[cluster != 'noise', cell_id], -c('id', 'patient', 'cell_type')]),
		cell_type_data[cluster != 'noise', cell_id]
	),
	cell_type_data[cluster != 'noise'],
	'noise_2',
	msigdb,
	seed = 7693
)

de_noise_2_top_50 <- de_noise_2$degenes[order(-rel_exp), head(gene, 50)]

de_noise_2_caf_endothelial <- test_cluster(
	set_rownames(
		as.matrix(sc_data[cell_type_data[cluster %in% c('caf', 'endothelial', 'noise_2'), cell_id], -c('id', 'patient', 'cell_type')]),
		cell_type_data[cluster %in% c('caf', 'endothelial', 'noise_2'), cell_id]
	),
	cell_type_data[cluster %in% c('caf', 'endothelial', 'noise_2')],
	'noise_2',
	msigdb,
	seed = 5880
)

de_noise_2_caf_endothelial_top_50 <- de_noise_2_caf_endothelial$degenes[order(-rel_exp), head(gene, 50)]

# Save stuff:
fwrite(plot_data, '../data_and_figures/sc_reclassify/files_for_rmd/crc_lee_smc_plot_data.csv')
saveRDS(de_5_top_50, '../data_and_figures/sc_reclassify/files_for_rmd/crc_lee_smc_de_5_top_50.rds')
saveRDS(de_caf_endothelial_top_50, '../data_and_figures/sc_reclassify/files_for_rmd/crc_lee_smc_de_caf_endothelial_top_50.rds')
saveRDS(de_between_caf_endothelial_top_50, '../data_and_figures/sc_reclassify/files_for_rmd/crc_lee_smc_de_between_caf_endothelial_top_50.rds')
saveRDS(de_noise_2_top_50, '../data_and_figures/sc_reclassify/files_for_rmd/crc_lee_smc_de_noise_2_top_50.rds')
saveRDS(de_noise_2_caf_endothelial_top_50, '../data_and_figures/sc_reclassify/files_for_rmd/crc_lee_smc_de_noise_2_caf_endothelial_top_50.rds')





# liver_ma:

classification_data <- fread(paste0('../data_and_figures/sc_reclassify/classification_data_liver_ma.csv'), key = 'cell_id')
sc_data <- fread('../data_and_figures/ma_liver_2019.csv')[, -c('sample', 'disease')][
	patient %in% c('C25', 'C26', 'C46', 'C56', 'C66', 'H37', 'H38', 'H65')
]
sc_data[, classification := classification_data[id, classification]]
sc_data <- sc_data[classification == 'nonmalignant', -'classification']
gene_averages <- sapply(sc_data[, -c('id', 'patient', 'cell_type')], function(x) {log2(mean(10*(2^x - 1)) + 1)}, USE.NAMES = TRUE)
sc_data <- sc_data[, c('id', 'patient', 'cell_type', names(gene_averages)[gene_averages >= 4]), with = FALSE]
sc_tsne <- readRDS(paste0('../data_and_figures/sc_reclassify/tsne_liver_ma.rds'))
sc_dbscan <- readRDS(paste0('../data_and_figures/sc_reclassify/dbscan_liver_ma.rds'))

plot_data <- setNames(
	cbind(as.data.table(sc_tsne$Y), as.character(sc_dbscan$cluster), sc_data$cell_type),
	c('x', 'y', 'dbscan_cluster', 'cell_type')
)

plot_data[dbscan_cluster == 0, dbscan_cluster := 'noise']

# Initial cell type definitions:

cell_type_data <- data.table(
	cell_id = sc_data$id,
	cluster = mapvalues(
		sc_dbscan$cluster,
		c(0, 5, 6),
		c('noise', 'endothelial_1', 'endothelial_2')
	)
)

cell_type_data[cluster == 2 & sc_data$cell_type != 'fibroblast', cluster := 'endothelial_3']
cell_type_data[cluster == 2 & sc_data$cell_type == 'fibroblast', cluster := 'caf']

setkey(sc_data, id)

de_caf <- test_cluster(
	set_rownames(
		as.matrix(sc_data[cell_type_data[cluster != 'noise', cell_id], -c('id', 'patient', 'cell_type')]),
		cell_type_data[cluster != 'noise', cell_id]
	),
	cell_type_data[cluster != 'noise'],
	'caf',
	msigdb,
	seed = 6826
)

de_caf_top_50 <- de_caf$degenes[order(-rel_exp), head(gene, 50)]

endothelial_seeds <- c(6327, 3648, 1966, 7305, 2122, 3678)

de_endothelial <- lapply(
	1:3,
	function(i) {
		test_cluster(
			set_rownames(
				as.matrix(
                    sc_data[
                        cell_type_data[!(cluster %in% c('noise', paste0('endothelial', (1:3)[1:3 != i]))), cell_id],
                        -c('id', 'patient', 'cell_type')
                    ]
                ),
				cell_type_data[!(cluster %in% c('noise', paste0('endothelial', (1:3)[1:3 != i]))), cell_id]
			),
			cell_type_data[!(cluster %in% c('noise', paste0('endothelial', (1:3)[1:3 != i])))],
			paste0('endothelial_', i),
			msigdb,
			seed = endothelial_seeds[i]
		)
	}
)

de_endothelial_top_50 <- set_colnames(
	sapply(1:3, function(i) de_endothelial[[i]]$degenes[order(-rel_exp), head(gene, 50)]),
	c('Cluster 5', 'Cluster 6', 'Cluster 2 - endothelial')
)

de_between_endothelial <- lapply(
	1:3,
	function(i) {
		test_cluster(
			set_rownames(
				as.matrix(sc_data[cell_type_data[startsWith(cluster, 'endothelial'), cell_id], -c('id', 'patient', 'cell_type')]),
				cell_type_data[startsWith(cluster, 'endothelial'), cell_id]
			),
			cell_type_data[startsWith(cluster, 'endothelial')],
			paste0('endothelial_', i),
			msigdb,
			seed = endothelial_seeds[i + 3]
		)
	}
)

de_between_endothelial_top_50 <- set_colnames(
	sapply(1:3, function(i) de_between_endothelial[[i]]$degenes[order(-rel_exp), head(gene, 50)]),
	c('Cluster 5', 'Cluster 6', 'Cluster 2 - endothelial')
)

# Save stuff:
fwrite(plot_data, '../data_and_figures/sc_reclassify/files_for_rmd/liver_ma_plot_data.csv')
saveRDS(de_caf_top_50, '../data_and_figures/sc_reclassify/files_for_rmd/liver_ma_de_caf_top_50.rds')
saveRDS(de_endothelial_top_50, '../data_and_figures/sc_reclassify/files_for_rmd/liver_ma_de_endothelial_top_50.rds')
saveRDS(de_between_endothelial_top_50, '../data_and_figures/sc_reclassify/files_for_rmd/liver_ma_de_between_endothelial_top_50.rds')





# lung_qian:

classification_data <- fread(paste0('../data_and_figures/sc_reclassify/classification_data_lung_qian.csv'), key = 'cell_id')
sc_data <- fread('../data_and_figures/qian_lung_2020.csv')[sample_type != 'normal', -c('sample_type', 'disease')]
sc_data[, classification := classification_data[id, classification]]
sc_data <- sc_data[classification == 'nonmalignant', -'classification']
gene_averages <- sapply(sc_data[, -c('id', 'patient', 'cell_type')], function(x) {log2(mean(10*(2^x - 1)) + 1)}, USE.NAMES = TRUE)
sc_data <- sc_data[, c('id', 'patient', 'cell_type', names(gene_averages)[gene_averages >= 4]), with = FALSE]
sc_tsne <- readRDS(paste0('../data_and_figures/sc_reclassify/tsne_lung_qian.rds'))
sc_dbscan <- readRDS(paste0('../data_and_figures/sc_reclassify/dbscan_lung_qian.rds'))

plot_data <- setNames(
	cbind(as.data.table(sc_tsne$Y), as.character(sc_dbscan$cluster), sc_data$cell_type),
	c('x', 'y', 'dbscan_cluster', 'cell_type')
)

plot_data[dbscan_cluster == 0, dbscan_cluster := 'noise']

# Initial cell type definitions:

cell_type_data <- data.table(
	cell_id = sc_data$id,
	cluster = mapvalues(
		sc_dbscan$cluster,
		c(0, 13, 18, 4, 19, 3, 14, 15, 5, 12, 6, 10, 16, 2, 9, 11, 17, 21, 1, 7, 20, 22),
		c('noise', 'b_cell', 'b_cell', 'b_plasma', 'b_plasma', 'caf_1', 'caf_2', 'caf_3', 'dc', 'endothelial', 'epithelial',
		  'epithelial', 'epithelial', 'macrophage', 'macrophage', 'macrophage', 'mast', 'mast', 't_cell', 't_cell', 't_cell', 't_cell')
	)
)

setkey(sc_data, id)

de_8 <- test_cluster(
	set_rownames(
		as.matrix(sc_data[cell_type_data[cluster != 'noise', cell_id], -c('id', 'patient', 'cell_type')]),
		cell_type_data[cluster != 'noise', cell_id]
	),
	cell_type_data[cluster != 'noise'],
	'8',
	msigdb,
	seed = 9611
)

de_8_top_50 <- de_8$degenes[order(-rel_exp), head(gene, 50)]

caf_seeds <- c(4573, 2106, 7294, 5951, 9657, 4862)

de_caf <- lapply(
	1:3,
	function(i) {
		test_cluster(
			set_rownames(
				as.matrix(
                    sc_data[
                        cell_type_data[cluster != 'noise' & !(cluster %in% paste0('caf_', (1:3)[1:3 != i])), cell_id],
                        -c('id', 'patient', 'cell_type')
                    ]
                ),
				cell_type_data[cluster != 'noise' & !(cluster %in% paste0('caf_', (1:3)[1:3 != i])), cell_id]
			),
			cell_type_data[cluster != 'noise' & !(cluster %in% paste0('caf_', (1:3)[1:3 != i]))],
			paste0('caf_', i),
			msigdb,
			seed = caf_seeds[i]
		)
	}
)

de_caf_top_50 <- set_colnames(
	sapply(1:3, function(i) de_caf[[i]]$degenes[order(-rel_exp), head(gene, 50)]),
	c('Cluster 3', 'Cluster 14', 'Cluster 15')
)

de_between_caf <- lapply(
	1:3,
	function(i) {
		test_cluster(
			set_rownames(
				as.matrix(sc_data[cell_type_data[startsWith(cluster, 'caf'), cell_id], -c('id', 'patient', 'cell_type')]),
				cell_type_data[startsWith(cluster, 'caf'), cell_id]
			),
			cell_type_data[startsWith(cluster, 'caf')],
			paste0('caf_', i),
			msigdb,
			seed = caf_seeds[i + 3]
		)
	}
)

de_between_caf_top_50 <- set_colnames(
	sapply(1:3, function(i) de_between_caf[[i]]$degenes[order(-rel_exp), head(gene, 50)]),
	c('Cluster 3', 'Cluster 14', 'Cluster 15')
)

fwrite(plot_data, '../data_and_figures/sc_reclassify/files_for_rmd/lung_qian_plot_data.csv')
saveRDS(de_8_top_50, '../data_and_figures/sc_reclassify/files_for_rmd/lung_qian_de_8_top_50.rds')
saveRDS(de_caf_top_50, '../data_and_figures/sc_reclassify/files_for_rmd/lung_qian_de_caf_top_50.rds')
saveRDS(de_between_caf_top_50, '../data_and_figures/sc_reclassify/files_for_rmd/lung_qian_de_between_caf_top_50.rds')





# ovarian_qian:

classification_data <- fread(paste0('../data_and_figures/sc_reclassify/classification_data_ovarian_qian.csv'), key = 'cell_id')
sc_data <- fread('../data_and_figures/qian_ovarian_2020.csv')[sample_type != 'normal', -c('sample_type', 'site')][patient %in% 11:14]
sc_data[, classification := classification_data[id, classification]]
sc_data <- sc_data[classification == 'nonmalignant', -'classification']
gene_averages <- sapply(sc_data[, -c('id', 'patient', 'cell_type')], function(x) {log2(mean(10*(2^x - 1)) + 1)}, USE.NAMES = TRUE)
sc_data <- sc_data[, c('id', 'patient', 'cell_type', names(gene_averages)[gene_averages >= 4]), with = FALSE]
sc_tsne <- readRDS(paste0('../data_and_figures/sc_reclassify/tsne_ovarian_qian.rds'))
sc_dbscan <- readRDS(paste0('../data_and_figures/sc_reclassify/dbscan_ovarian_qian.rds'))

plot_data <- setNames(
	cbind(as.data.table(sc_tsne$Y), as.character(sc_dbscan$cluster), sc_data$cell_type),
	c('x', 'y', 'dbscan_cluster', 'cell_type')
)

plot_data[dbscan_cluster == 0, dbscan_cluster := 'noise']

# Initial cell type definitions:

cell_type_data <- data.table(
	cell_id = sc_data$id,
	cluster = mapvalues(
		sc_dbscan$cluster,
		c(0, 13, 5, 14, 6, 3, 19, 21, 7, 12, 11, 16, 1, 2, 4, 10, 15, 17, 18, 20, 8),
		c('noise', 'b_cell', 'endothelial_1', 'endothelial_1', 'endothelial_2', 'macrophage', 'macrophage', 'macrophage', 't_cell', 't_cell',
		  'epithelial', 'epithelial', 'caf_1', 'caf_1', 'caf_1', 'caf_1', 'caf_2', 'caf_3', 'caf_4', 'caf_5', 'caf_6')
	)
)

setkey(sc_data, id)

de_9 <- test_cluster(
	set_rownames(
		as.matrix(sc_data[cell_type_data[cluster != 'noise', cell_id], -c('id', 'patient', 'cell_type')]),
		cell_type_data[cluster != 'noise', cell_id]
	),
	cell_type_data[cluster != 'noise'],
	'9',
	msigdb,
	seed = 1575
)

de_9_top_50 <- de_9$degenes[order(-rel_exp), head(gene, 100)]

caf_seeds <- c(4105, 530, 4113, 6442, 6390, 8160, 6406, 9617, 9597, 9388, 1284, 9511)

de_caf <- lapply(
	1:6,
	function(i) {
		test_cluster(
			set_rownames(
				as.matrix(
                    sc_data[
                        cell_type_data[cluster != 'noise' & !(cluster %in% paste0('caf_', (1:6)[1:6 != i])), cell_id],
                        -c('id', 'patient', 'cell_type')
                    ]
                ),
				cell_type_data[cluster != 'noise' & !(cluster %in% paste0('caf_', (1:6)[1:6 != i])), cell_id]
			),
			cell_type_data[cluster != 'noise' & !(cluster %in% paste0('caf_', (1:6)[1:6 != i]))],
			paste0('caf_', i),
			msigdb,
			seed = caf_seeds[i]
		)
	}
)

de_caf_top_50 <- set_colnames(
	sapply(1:6, function(i) de_caf[[i]]$degenes[order(-rel_exp), head(gene, 50)]),
	c('Clusters 1, 2, 4 and 10', 'Cluster 15', 'Cluster 17', 'Cluster 18', 'Cluster 20', 'Cluster 8')
)

de_between_caf <- lapply(
	1:6,
	function(i) {
		test_cluster(
			set_rownames(
				as.matrix(sc_data[cell_type_data[startsWith(cluster, 'caf'), cell_id], -c('id', 'patient', 'cell_type')]),
				cell_type_data[startsWith(cluster, 'caf'), cell_id]
			),
			cell_type_data[startsWith(cluster, 'caf')],
			paste0('caf_', i),
			msigdb,
			seed = caf_seeds[i + 6]
		)
	}
)

de_between_caf_top_50 <- set_colnames(
	sapply(1:6, function(i) de_between_caf[[i]]$degenes[order(-rel_exp), head(gene, 50)]),
	c('Clusters 1, 2, 4 and 10', 'Cluster 15', 'Cluster 17', 'Cluster 18', 'Cluster 20', 'Cluster 8')
)

endothelial_seeds <- c(5697, 6179, 2029, 9568)

de_endothelial <- lapply(
	1:2,
	function(i) {
		test_cluster(
			set_rownames(
				as.matrix(
                    sc_data[
                        cell_type_data[cluster != 'noise' & !(cluster %in% paste0('endothelial_', (1:2)[1:2 != i])), cell_id],
                        -c('id', 'patient', 'cell_type')
                    ]
                ),
				cell_type_data[cluster != 'noise' & !(cluster %in% paste0('endothelial_', (1:2)[1:2 != i])), cell_id]
			),
			cell_type_data[cluster != 'noise' & !(cluster %in% paste0('endothelial_', (1:2)[1:2 != i]))],
			paste0('endothelial_', i),
			msigdb,
			seed = endothelial_seeds[i]
		)
	}
)

de_endothelial_top_50 <- set_colnames(
	sapply(1:2, function(i) de_endothelial[[i]]$degenes[order(-rel_exp), head(gene, 50)]),
	c('Clusters 5 and 14', 'Cluster 6')
)

de_between_endothelial <- lapply(
	1:2,
	function(i) {
		test_cluster(
			set_rownames(
				as.matrix(sc_data[cell_type_data[startsWith(cluster, 'endothelial'), cell_id], -c('id', 'patient', 'cell_type')]),
				cell_type_data[startsWith(cluster, 'endothelial'), cell_id]
			),
			cell_type_data[startsWith(cluster, 'endothelial')],
			paste0('endothelial_', i),
			msigdb,
			seed = endothelial_seeds[i + 2]
		)
	}
)

de_between_endothelial_top_50 <- set_colnames(
	sapply(1:2, function(i) de_between_endothelial[[i]]$degenes[order(-rel_exp), head(gene, 50)]),
	c('Clusters 5 and 14', 'Cluster 6')
)

fwrite(plot_data, '../data_and_figures/sc_reclassify/files_for_rmd/ovarian_qian_plot_data.csv')
saveRDS(de_9_top_50, '../data_and_figures/sc_reclassify/files_for_rmd/ovarian_qian_de_9_top_50.rds')
saveRDS(de_caf_top_50, '../data_and_figures/sc_reclassify/files_for_rmd/ovarian_qian_de_caf_top_50.rds')
saveRDS(de_between_caf_top_50, '../data_and_figures/sc_reclassify/files_for_rmd/ovarian_qian_de_between_caf_top_50.rds')
saveRDS(de_endothelial_top_50, '../data_and_figures/sc_reclassify/files_for_rmd/ovarian_qian_de_endothelial_top_50.rds')
saveRDS(de_between_endothelial_top_50, '../data_and_figures/sc_reclassify/files_for_rmd/ovarian_qian_de_between_endothelial_top_50.rds')





# pdac_peng:

classification_data <- fread(paste0('../data_and_figures/sc_reclassify/classification_data_pdac_peng.csv'), key = 'cell_id')
sc_data <- fread('../data_and_figures/peng_pdac_2019.csv')[sample_type != 'normal', -'sample_type'][
	patient %in% c('T2', 'T3', 'T6', 'T7', 'T8', 'T9', 'T11', 'T13', 'T14', 'T15', 'T16', 'T17', 'T19', 'T20', 'T21', 'T22')
]
sc_data[, classification := classification_data[id, classification]]
sc_data <- sc_data[classification == 'nonmalignant', -'classification']
gene_averages <- sapply(sc_data[, -c('id', 'patient', 'cell_type')], function(x) {log2(mean(10*(2^x - 1)) + 1)}, USE.NAMES = TRUE)
sc_data <- sc_data[, c('id', 'patient', 'cell_type', names(gene_averages)[gene_averages >= 4]), with = FALSE]
sc_tsne <- readRDS(paste0('../data_and_figures/sc_reclassify/tsne_pdac_peng.rds'))
sc_dbscan <- readRDS(paste0('../data_and_figures/sc_reclassify/dbscan_pdac_peng.rds'))

plot_data <- setNames(
	cbind(as.data.table(sc_tsne$Y), as.character(sc_dbscan$cluster), sc_data$cell_type),
	c('x', 'y', 'dbscan_cluster', 'cell_type')
)

plot_data[dbscan_cluster == 0, dbscan_cluster := 'noise']

# Initial cell type definitions:

cell_type_data <- data.table(
	cell_id = sc_data$id,
	cluster = mapvalues(
		sc_dbscan$cluster,
		c(0, 10, 3, 5, 9, 1, 4, 7, 12, 13, 14, 6, 11),
		c('noise', 'b_plasma', 'endothelial', 'epithelial', 'epithelial', 'macrophage', 'caf_1', 'caf_2', 'caf_3', 'caf_4', 'caf_5',
		  'b_cell', 'b_cell')
	)
)

cell_type_data[cluster == 2 & sc_data$cell_type == 't_cell', cluster := 't_cell']
cell_type_data[cluster == 2 & sc_data$cell_type == 'b_cell', cluster := 'b_cell']
cell_type_data[cluster == 2 & !(sc_data$cell_type %in% c('b_cell', 't_cell')), cluster := 't_cell']

setkey(sc_data, id)

caf_seeds <- c(1160, 1583, 2746, 9909, 7037, 2347, 9093, 7311, 9648, 7849)

de_caf <- lapply(
	1:5,
	function(i) {
		test_cluster(
			set_rownames(
				as.matrix(
                    sc_data[
                        cell_type_data[cluster != 'noise' & !(cluster %in% paste0('caf_', (1:5)[1:5 != i])), cell_id],
                        -c('id', 'patient', 'cell_type')
                    ]
                ),
				cell_type_data[cluster != 'noise' & !(cluster %in% paste0('caf_', (1:5)[1:5 != i])), cell_id]
			),
			cell_type_data[cluster != 'noise' & !(cluster %in% paste0('caf_', (1:5)[1:5 != i]))],
			paste0('caf_', i),
			msigdb,
			seed = caf_seeds[i]
		)
	}
)

de_caf_top_50 <- set_colnames(
	sapply(1:5, function(i) de_caf[[i]]$degenes[order(-rel_exp), head(gene, 50)]),
	c('Cluster 4', 'Cluster 7', 'Cluster 12', 'Cluster 13', 'Cluster 14')
)

de_between_caf <- lapply(
	1:5,
	function(i) {
		test_cluster(
			set_rownames(
				as.matrix(sc_data[cell_type_data[startsWith(cluster, 'caf'), cell_id], -c('id', 'patient', 'cell_type')]),
				cell_type_data[startsWith(cluster, 'caf'), cell_id]
			),
			cell_type_data[startsWith(cluster, 'caf')],
			paste0('caf_', i),
			msigdb,
			seed = caf_seeds[i + 5]
		)
	}
)

de_between_caf_top_50 <- set_colnames(
	sapply(1:5, function(i) de_between_caf[[i]]$degenes[order(-rel_exp), head(gene, 50)]),
	c('Cluster 4', 'Cluster 7', 'Cluster 12', 'Cluster 13', 'Cluster 14')
)

fwrite(plot_data, '../data_and_figures/sc_reclassify/files_for_rmd/pdac_peng_plot_data.csv')
saveRDS(de_caf_top_50, '../data_and_figures/sc_reclassify/files_for_rmd/pdac_peng_de_caf_top_50.rds')
saveRDS(de_between_caf_top_50, '../data_and_figures/sc_reclassify/files_for_rmd/pdac_peng_de_between_caf_top_50.rds')
