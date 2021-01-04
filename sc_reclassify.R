library(data.table) # 1.12.8
library(ggplot2) # 3.3.0
library(magrittr) # 1.5
library(Rtsne) # 0.15
library(plyr) # 1.8.6
library(msigdbr) # 7.1.1
library(clusterProfiler) # 3.14.3
library(cowplot) # 1.0.0

source('general_functions.R')

cell_type_markers <- fread('../../cell_type_markers.csv')





cohort_data <- list(
	breast_qian = list(
		patients = c(42, 43, 47, 49, 51, 53, 54),
		read_quote = quote(fread('../data_and_figures/qian_breast_2020.csv')),
		seed = 6876,
		min_cells = 50,
		epsilon = 3.1
	),
	crc_lee_smc = list(
		patients = c('SMC01', 'SMC02', 'SMC04', 'SMC07', 'SMC08', 'SMC09', 'SMC11', 'SMC14', 'SMC15', 'SMC16', 'SMC18', 'SMC20', 'SMC21', 'SMC23',
					 'SMC25'),
		read_quote = quote(fread('../data_and_figures/lee_crc_2020_smc.csv')[
			sample_type != 'normal',
			-c('sample_id', 'sample_type', 'cell_subtype')
		]),
		seed = 4291,
		min_cells = 50,
		epsilon = 3.1
	),
	hnscc_puram = list(
		patients = c(5, 6, 18, 20, 22, 25, 26),
		read_quote = quote(fread('../data_and_figures/puram_hnscc_2017.csv')[, -c('lymph_node', 'processed_by_maxima_enzyme')]),
		seed = 1521,
		min_cells = 20,
		epsilon = 3.2
	),
	liver_ma = list(
		patients = c('C25', 'C26', 'C46', 'C56', 'C66', 'H37', 'H38', 'H65'),
		read_quote = quote(fread('../data_and_figures/ma_liver_2019.csv')[, -c('sample', 'disease')]),
		seed = 7568,
		min_cells = 20,
		epsilon = 3.1
	),
	luad_kim = list(
		patients = c('P0006', 'P0008', 'P0018', 'P0019', 'P0020', 'P0025', 'P0028', 'P0030', 'P0031', 'P0034', 'P1006', 'P1028', 'P1049', 'P1058'),
		read_quote = quote(
			fread('../data_and_figures/kim_luad_2020.csv')[, -c('cell_type_refined', 'cell_subtype', 'sample_site', 'sample_id', 'sample_origin')]
		),
		seed = 9684,
		min_cells = 50,
		epsilon = 2.4
	),
	lung_qian = list(
		patients = 1:8,
		read_quote = quote(fread('../data_and_figures/qian_lung_2020.csv')[sample_type != 'normal', -c('sample_type', 'disease')]),
		seed = 824,
		min_cells = 50,
		epsilon = 2.4
	),
	ovarian_qian = list(
		patients = 11:14,
		read_quote = quote(fread('../data_and_figures/qian_ovarian_2020.csv')[sample_type != 'normal', -c('sample_type', 'site')]),
		seed = 2802,
		min_cells = 20,
		epsilon = 1.9
	),
	pdac_peng = list(
		patients = c('T2', 'T3', 'T6', 'T7', 'T8', 'T9', 'T11', 'T13', 'T14', 'T15', 'T16', 'T17', 'T19', 'T20', 'T21', 'T22'),
		read_quote = quote(fread('../data_and_figures/peng_pdac_2019.csv')[sample_type != 'normal', -'sample_type']),
		seed = 3997,
		min_cells = 50,
		epsilon = 2.6
	)
)





for(cohort in names(cohort_data)) {

	classification_data <- fread(paste0('../data_and_figures/sc_reclassify/classification_data_', cohort, '.csv'), key = 'cell_id')

	sc_data <- eval(cohort_data[[cohort]]$read_quote)[patient %in% cohort_data[[cohort]]$patients]

	sc_data[, classification := classification_data[id, classification]]

	sc_data <- sc_data[classification == 'nonmalignant', -'classification']

	gene_averages <- sapply(sc_data[, -c('id', 'patient', 'cell_type')], function(x) {log2(mean(10*(2^x - 1)) + 1)}, USE.NAMES = TRUE)

	sc_data <- sc_data[, c('id', 'patient', 'cell_type', names(gene_averages)[gene_averages >= 4]), with = FALSE]

	sc_tsne <- readRDS(paste0('../data_and_figures/sc_reclassify/tsne_', cohort, '.rds'))

	# Check t-SNE plot for obvious problems:
	# qplot(x, y, data = setNames(as.data.table(sc_tsne$Y), c('x', 'y')))

	# Choose minimum number of cells per cluster and check the k-NN distance plot:
	# dbscan::kNNdistplot(sc_tsne$Y, k = cohort_data[[cohort]]$min_cells - 1)
	# abline(h = cohort_data[[cohort]]$epsilon)

	# Run DBSCAN:
	set.seed(cohort_data[[cohort]]$seed)
	sc_dbscan <- fpc::dbscan(sc_tsne$Y, eps = cohort_data[[cohort]]$epsilon, MinPts = cohort_data[[cohort]]$min_cells)

	saveRDS(sc_dbscan, paste0('../data_and_figures/sc_reclassify/dbscan_', cohort, '.rds'))

	plot_data <- setNames(
		cbind(as.data.table(sc_tsne$Y), as.character(sc_dbscan$cluster), sc_data$cell_type),
		c('x', 'y', 'dbscan_cluster', 'cell_type')
	)

	plot_data[dbscan_cluster == 0, dbscan_cluster := 'noise']

	# Scatter plot of t-SNE coordinates, coloured by DBSCAN cluster:
	tsne_plot_dbscan <- ggplot(plot_data, aes(x = x, y = y, colour = dbscan_cluster)) +
		geom_point() +
		# scale_colour_manual(values = randomcoloR::distinctColorPalette(length(unique(sc_dbscan$cluster)))) +
		theme_minimal()

	tsne_plot_author_cell_types <- ggplot(plot_data, aes(x = x, y = y, colour = cell_type)) +
		geom_point() +
		# scale_colour_manual(values = randomcoloR::distinctColorPalette(length(unique(sc_data$cell_type)))) +
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

	epithelial_markers <- c('CDH1', 'EPCAM', 'SFN', names(sc_data)[grepl('^KRT[0-9]', names(sc_data))])

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
	for(x in ct_ave_exp_plots) print(x$plot)

	dev.off()

}





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

	# Centre the matrix gene-wise (not sure this is necessary):
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





cohort <- 'breast_qian'

classification_data <- fread(paste0('../data_and_figures/sc_reclassify/classification_data_', cohort, '.csv'), key = 'cell_id')
sc_data <- eval(cohort_data[[cohort]]$read_quote)[patient %in% cohort_data[[cohort]]$patients]
sc_data[, classification := classification_data[id, classification]]
sc_data <- sc_data[classification == 'nonmalignant', -'classification']
gene_averages <- sapply(sc_data[, -c('id', 'patient', 'cell_type')], function(x) {log2(mean(10*(2^x - 1)) + 1)}, USE.NAMES = TRUE)
sc_data <- sc_data[, c('id', 'patient', 'cell_type', names(gene_averages)[gene_averages >= 4]), with = FALSE]
sc_tsne <- readRDS(paste0('../data_and_figures/sc_reclassify/tsne_', cohort, '.rds'))
sc_dbscan <- readRDS(paste0('../data_and_figures/sc_reclassify/dbscan_', cohort, '.rds'))

# Initial cell type definitions:
cell_type_data <- data.table(
	cell_id = sc_data$id,
	cluster = mapvalues(
		sc_dbscan$cluster,
		c(0, 8, 4, 6, 3, 9, 1, 7, 10, 2, 5),
		c('noise', 'b_cell', 'b_plasma', 'endothelial', 'macrophage', 'mast', 't_cell', 't_cell', 't_cell', 'caf_candidate_1', 'caf_candidate_2')
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

de_11$degenes[order(-rel_exp), head(gene,100)]

# This shows markers for various immune cell types, so possibly these cells are doublets.

caf_candidate_seeds <- c(4843, 6826, 6365, 878)

de_caf_candidate <- lapply(
	1:2,
	function(i) {
		test_cluster(
			set_rownames(
				as.matrix(sc_data[cell_type_data[cluster != 'noise', cell_id], -c('id', 'patient', 'cell_type')]),
				cell_type_data[cluster != 'noise', cell_id]
			),
			cell_type_data[cluster != 'noise'],
			paste0('caf_candidate_', i),
			msigdb,
			seed = caf_candidate_seeds[i]
		)
	}
)

sapply(1:2, function(i) de_caf_candidate[[i]]$degenes[order(-rel_exp), head(gene, 50)])

de_between_caf_candidate <- lapply(
	1:2,
	function(i) {
		test_cluster(
			set_rownames(
				as.matrix(sc_data[cell_type_data[startsWith(cluster, 'caf_candidate'), cell_id], -c('id', 'patient', 'cell_type')]),
				cell_type_data[startsWith(cluster, 'caf_candidate'), cell_id]
			),
			cell_type_data[startsWith(cluster, 'caf_candidate')],
			paste0('caf_candidate_', i),
			msigdb,
			seed = caf_candidate_seeds[i + 2]
		)
	}
)

sapply(1:2, function(i) de_between_caf_candidate[[i]]$degenes[order(-rel_exp), head(gene, 50)])

# It looks like caf_candidate_2 (DBSCAN cluster 5) could be pericytes, since RGS5 comes up prominently in this cluster.





cohort <- 'crc_lee_smc'

classification_data <- fread(paste0('../data_and_figures/sc_reclassify/classification_data_', cohort, '.csv'), key = 'cell_id')
sc_data <- eval(cohort_data[[cohort]]$read_quote)[patient %in% cohort_data[[cohort]]$patients]
sc_data[, classification := classification_data[id, classification]]
sc_data <- sc_data[classification == 'nonmalignant', -'classification']
gene_averages <- sapply(sc_data[, -c('id', 'patient', 'cell_type')], function(x) {log2(mean(10*(2^x - 1)) + 1)}, USE.NAMES = TRUE)
sc_data <- sc_data[, c('id', 'patient', 'cell_type', names(gene_averages)[gene_averages >= 4]), with = FALSE]
sc_tsne <- readRDS(paste0('../data_and_figures/sc_reclassify/tsne_', cohort, '.rds'))
sc_dbscan <- readRDS(paste0('../data_and_figures/sc_reclassify/dbscan_', cohort, '.rds'))

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

de_5$degenes[order(-rel_exp), head(gene, 100)]

# Looks like something immune, similar to cluster 11 in breast_qian.  Very few gene sets come up in GSEA.

# Check that clusters 2 and 3 definitely are CAFs and endothelial cells, respectively:

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

sapply(1:2, function(i) de_caf_endothelial[[i]]$degenes[order(-rel_exp), head(gene, 50)])

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

sapply(1:2, function(i) de_between_caf_endothelial[[i]]$degenes[order(-rel_exp), head(gene, 50)])

noise_seeds <- c(7172, 7693)

de_noise <- lapply(
	1:2,
	function(i) {
		test_cluster(
			set_rownames(
				as.matrix(sc_data[cell_type_data[cluster != 'noise', cell_id], -c('id', 'patient', 'cell_type')]),
				cell_type_data[cluster != 'noise', cell_id]
			),
			cell_type_data[cluster != 'noise'],
			paste0('noise_', i),
			msigdb,
			seed = noise_seeds[i]
		)
	}
)

sapply(1:2, function(i) de_noise[[i]]$degenes[order(-rel_exp), head(gene, 50)])

de_noise_1_b_cell <- test_cluster(
	set_rownames(
		as.matrix(sc_data[cell_type_data[cluster %in% c('b_cell', 'b_plasma', 'noise_1'), cell_id], -c('id', 'patient', 'cell_type')]),
		cell_type_data[cluster %in% c('b_cell', 'b_plasma', 'noise_1'), cell_id]
	),
	cell_type_data[cluster %in% c('b_cell', 'b_plasma', 'noise_1')],
	'noise_1',
	msigdb,
	seed = 9527
)

de_noise_1_b_cell$degenes[order(-rel_exp), head(gene, 50)]

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

de_noise_2_caf_endothelial$degenes[order(-rel_exp), head(gene, 50)]

# I think noise_2 could be pericytes - RGS5 is top of both lists, and I see PECAM1 as well.

# noise_1 indeed looks like B cells, since the top genes consist of several immunoglobulins, but also the B plasma marker expression is fairly high in
# these cells.





cohort <- 'liver_ma'

classification_data <- fread(paste0('../data_and_figures/sc_reclassify/classification_data_', cohort, '.csv'), key = 'cell_id')
sc_data <- eval(cohort_data[[cohort]]$read_quote)[patient %in% cohort_data[[cohort]]$patients]
sc_data[, classification := classification_data[id, classification]]
sc_data <- sc_data[classification == 'nonmalignant', -'classification']
gene_averages <- sapply(sc_data[, -c('id', 'patient', 'cell_type')], function(x) {log2(mean(10*(2^x - 1)) + 1)}, USE.NAMES = TRUE)
sc_data <- sc_data[, c('id', 'patient', 'cell_type', names(gene_averages)[gene_averages >= 4]), with = FALSE]
sc_tsne <- readRDS(paste0('../data_and_figures/sc_reclassify/tsne_', cohort, '.rds'))
sc_dbscan <- readRDS(paste0('../data_and_figures/sc_reclassify/dbscan_', cohort, '.rds'))

# Initial cell type definitions:
cell_type_data <- data.table(
	cell_id = sc_data$id,
	cluster = mapvalues(sc_dbscan$cluster, c(0, 5, 6), c('noise', 'endothelial_1', 'endothelial_2'))
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

de_caf$degenes[order(-rel_exp), head(gene, 100)]

# Looks like myofibroblasts.

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

sapply(1:3, function(i) de_endothelial[[i]]$degenes[order(-rel_exp), head(gene, 50)])

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

sapply(1:3, function(i) de_between_endothelial[[i]]$degenes[order(-rel_exp), head(gene, 50)])





cohort <- 'luad_kim'

classification_data <- fread(paste0('../data_and_figures/sc_reclassify/classification_data_', cohort, '.csv'), key = 'cell_id')
sc_data <- eval(cohort_data[[cohort]]$read_quote)[patient %in% cohort_data[[cohort]]$patients]
sc_data[, classification := classification_data[id, classification]]
sc_data <- sc_data[classification == 'nonmalignant', -'classification']
gene_averages <- sapply(sc_data[, -c('id', 'patient', 'cell_type')], function(x) {log2(mean(10*(2^x - 1)) + 1)}, USE.NAMES = TRUE)
sc_data <- sc_data[, c('id', 'patient', 'cell_type', names(gene_averages)[gene_averages >= 4]), with = FALSE]
sc_tsne <- readRDS(paste0('../data_and_figures/sc_reclassify/tsne_', cohort, '.rds'))
sc_dbscan <- readRDS(paste0('../data_and_figures/sc_reclassify/dbscan_', cohort, '.rds'))

# Initial cell type definitions:
cell_type_data <- data.table(
	cell_id = sc_data$id,
	cluster = mapvalues(
		sc_dbscan$cluster,
		c(0, 1, 8, 11, 19, 9, 16, 7, 15, 13, 2, 17, 3, 10, 12, 20, 4, 5, 6, 14),
		c('noise', 'b_cell', 'b_cell', 'b_cell', 'b_cell', 'b_plasma', 'b_plasma', 'caf', 'dc', 'endothelial', 'epithelial', 'epithelial',
		  'macrophage', 'macrophage', 'macrophage', 'macrophage', 'mast', 't_cell', 't_cell', 't_cell')
	)
)

setkey(sc_data, id)

de_18 <- test_cluster(
	set_rownames(
		as.matrix(sc_data[cell_type_data[cluster != 'noise', cell_id], -c('id', 'patient', 'cell_type')]),
		cell_type_data[cluster != 'noise', cell_id]
	),
	cell_type_data[cluster != 'noise'],
	'18',
	msigdb,
	seed = 9611
)

de_18$degenes[order(-rel_exp), head(gene, 100)]

# Immune-like cluster that I can't identify





cohort <- 'lung_qian'

classification_data <- fread(paste0('../data_and_figures/sc_reclassify/classification_data_', cohort, '.csv'), key = 'cell_id')
sc_data <- eval(cohort_data[[cohort]]$read_quote)[patient %in% cohort_data[[cohort]]$patients]
sc_data[, classification := classification_data[id, classification]]
sc_data <- sc_data[classification == 'nonmalignant', -'classification']
gene_averages <- sapply(sc_data[, -c('id', 'patient', 'cell_type')], function(x) {log2(mean(10*(2^x - 1)) + 1)}, USE.NAMES = TRUE)
sc_data <- sc_data[, c('id', 'patient', 'cell_type', names(gene_averages)[gene_averages >= 4]), with = FALSE]
sc_tsne <- readRDS(paste0('../data_and_figures/sc_reclassify/tsne_', cohort, '.rds'))
sc_dbscan <- readRDS(paste0('../data_and_figures/sc_reclassify/dbscan_', cohort, '.rds'))

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

de_8$degenes[order(-rel_exp), head(gene, 100)]

# Immune-like cluster that I can't identify

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

sapply(1:3, function(i) de_caf[[i]]$degenes[order(-rel_exp), head(gene, 50)])

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

sapply(1:3, function(i) de_between_caf[[i]]$degenes[order(-rel_exp), head(gene, 50)])

# caf_3 may be pericytes, because RGS5 comes up in both the general and between-caf-cluster tests (also PDGFRB appears in the between test).





cohort <- 'ovarian_qian'

classification_data <- fread(paste0('../data_and_figures/sc_reclassify/classification_data_', cohort, '.csv'), key = 'cell_id')
sc_data <- eval(cohort_data[[cohort]]$read_quote)[patient %in% cohort_data[[cohort]]$patients]
sc_data[, classification := classification_data[id, classification]]
sc_data <- sc_data[classification == 'nonmalignant', -'classification']
gene_averages <- sapply(sc_data[, -c('id', 'patient', 'cell_type')], function(x) {log2(mean(10*(2^x - 1)) + 1)}, USE.NAMES = TRUE)
sc_data <- sc_data[, c('id', 'patient', 'cell_type', names(gene_averages)[gene_averages >= 4]), with = FALSE]
sc_tsne <- readRDS(paste0('../data_and_figures/sc_reclassify/tsne_', cohort, '.rds'))
sc_dbscan <- readRDS(paste0('../data_and_figures/sc_reclassify/dbscan_', cohort, '.rds'))

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

de_9$degenes[order(-rel_exp), head(gene, 100)]

# Unidentified immune cluster

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

sapply(1:6, function(i) de_caf[[i]]$degenes[order(-rel_exp), head(gene, 50)])

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

sapply(1:6, function(i) de_between_caf[[i]]$degenes[order(-rel_exp), head(gene, 50)])

# caf_6 (DBSCAN cluster 8) looks like pericytes or similar; caf_1, caf_4 and caf_5 all look like CAFs.  caf_2 and caf_3 are harder to identify.

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

sapply(1:2, function(i) de_endothelial[[i]]$degenes[order(-rel_exp), head(gene, 50)])

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

sapply(1:2, function(i) de_between_endothelial[[i]]$degenes[order(-rel_exp), head(gene, 50)])

# No clear differences between the endothelial clusters.





cohort <- 'pdac_peng'

classification_data <- fread(paste0('../data_and_figures/sc_reclassify/classification_data_', cohort, '.csv'), key = 'cell_id')
sc_data <- eval(cohort_data[[cohort]]$read_quote)[patient %in% cohort_data[[cohort]]$patients]
sc_data[, classification := classification_data[id, classification]]
sc_data <- sc_data[classification == 'nonmalignant', -'classification']
gene_averages <- sapply(sc_data[, -c('id', 'patient', 'cell_type')], function(x) {log2(mean(10*(2^x - 1)) + 1)}, USE.NAMES = TRUE)
sc_data <- sc_data[, c('id', 'patient', 'cell_type', names(gene_averages)[gene_averages >= 4]), with = FALSE]
sc_tsne <- readRDS(paste0('../data_and_figures/sc_reclassify/tsne_', cohort, '.rds'))
sc_dbscan <- readRDS(paste0('../data_and_figures/sc_reclassify/dbscan_', cohort, '.rds'))

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

de_8 <- test_cluster(
	set_rownames(
		as.matrix(sc_data[cell_type_data[cluster != 'noise', cell_id], -c('id', 'patient', 'cell_type')]),
		cell_type_data[cluster != 'noise', cell_id]
	),
	cell_type_data[cluster != 'noise'],
	'8',
	msigdb,
	seed = 7370
)

de_8$degenes[order(-rel_exp), head(gene, 100)]

# The authors call these cells "endocrine", which fits with the GSEA analysis.

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

sapply(1:5, function(i) de_caf[[i]]$degenes[order(-rel_exp), head(gene, 50)])

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

sapply(1:5, function(i) de_between_caf[[i]]$degenes[order(-rel_exp), head(gene, 50)])

# caf_1 and caf_3 (DBSCAN clusters 4 and 12) both have RGS5 among their DE genes. The authors call them stellate cells, which makes sense.  The other
# three clusters all look like CAFs.





# Final cell type assignments:





# breast_qian:

sc_data <- fread('../data_and_figures/qian_breast_2020.csv')[patient %in% c(42, 43, 47, 49, 51, 53, 54)]
classification_data <- fread(paste0('../data_and_figures/sc_reclassify/classification_data_breast_qian.csv'), key = 'cell_id')
sc_data[, classification := classification_data[id, classification]]
sc_dbscan <- readRDS(paste0('../data_and_figures/sc_reclassify/dbscan_breast_qian.rds'))
sc_data[classification == 'nonmalignant', dbscan_cluster := sc_dbscan$cluster]

# We want to make strict and lenient classifications for CAFs (the former excluding cluster 5 and the latter including it), so duplicate the
# classification column:
sc_data[, classification_lenient := classification]

setcolorder(sc_data, c('id', 'patient', 'cell_type', 'classification', 'classification_lenient', 'dbscan_cluster'))

# Here I'm excluding the nonmalignant cells that the authors classified as malignant along with the cells in cluster 11 (possible doublets)
# and cluster 5 (possible pericytes).
sc_data[classification == 'nonmalignant' & (cell_type == 'Cancer' | dbscan_cluster %in% c(5, 11)), classification := 'ambiguous']

# Similar for the lenient classification column but retaining cluster 5:
sc_data[classification_lenient == 'nonmalignant' & cell_type == 'Cancer' | dbscan_cluster == 11, classification_lenient := 'ambiguous']

# Change remaining nonmalignant classifications to the authors' classifications, modifying names to my taste:
sc_data[
	classification == 'nonmalignant',
	classification := mapvalues(
		cell_type,
		c('B_cell', 'DC', 'EC', 'Fibroblast', 'Mast', 'Myeloid', 'T_cell'),
		c('b_cell', 'dendritic', 'endothelial', 'caf', 'mast', 'macrophage', 't_cell')
	)
]

# Similarly for classification_lenient (which will include cluster 5 in the CAF cells):
sc_data[
	classification_lenient == 'nonmalignant',
	classification_lenient := mapvalues(
		cell_type,
		c('B_cell', 'DC', 'EC', 'Fibroblast', 'Mast', 'Myeloid', 'T_cell'),
		c('b_cell', 'dendritic', 'endothelial', 'caf', 'mast', 'macrophage', 't_cell')
	)
]

# Exclude malignant cells that the authors classified as nonmalignant cell types:
sc_data[classification == 'malignant' & cell_type != 'Cancer', classification := 'ambiguous']
sc_data[classification_lenient == 'malignant' & cell_type != 'Cancer', classification_lenient := 'ambiguous']

# Change "malignant" to "cancer":
sc_data[classification == 'malignant', classification := 'cancer']
sc_data[classification_lenient == 'malignant', classification_lenient := 'cancer']

sc_data[, dbscan_cluster := NULL]

setnames(sc_data, c('cell_type', 'classification', 'classification_lenient'), c('cell_type_author', 'cell_type', 'cell_type_lenient'))

fwrite(sc_data, '../data_and_figures/qian_breast_2020_reclassified.csv')





# crc_lee_smc:

sc_data <- fread('../data_and_figures/lee_crc_2020_smc.csv')[
	sample_type != 'normal',
	-c('sample_id', 'sample_type', 'cell_subtype')
][
	patient %in% c('SMC01', 'SMC02', 'SMC04', 'SMC07', 'SMC08', 'SMC09', 'SMC11', 'SMC14',
				   'SMC15', 'SMC16', 'SMC18', 'SMC20', 'SMC21', 'SMC23', 'SMC25')
]
classification_data <- fread(paste0('../data_and_figures/sc_reclassify/classification_data_crc_lee_smc.csv'), key = 'cell_id')
sc_data[, classification := classification_data[id, classification]]
sc_dbscan <- readRDS(paste0('../data_and_figures/sc_reclassify/dbscan_crc_lee_smc.rds'))
sc_data[classification == 'nonmalignant', dbscan_cluster := sc_dbscan$cluster]

# Duplicate column for lenient classifications for CAFs (including cluster 3 and the little noise cluster):
sc_data[, classification_lenient := classification]

setcolorder(sc_data, c('id', 'patient', 'cell_type', 'classification', 'classification_lenient', 'dbscan_cluster'))

# Merge the small cluster of noise cells near the endothelial cluster with the endothelial cluster:
sc_tsne <- readRDS(paste0('../data_and_figures/sc_reclassify/tsne_crc_lee_smc.rds'))
sc_data[classification == 'nonmalignant'][dbscan_cluster == 0 & sc_tsne$Y[, 2] > 40, dbscan_cluster := 3]

# Exclude cluster 5 (not sure this is appropriate, but we can always reverse it):
sc_data[dbscan_cluster == 5, classification := 'ambiguous']
sc_data[dbscan_cluster == 5, classification_lenient := 'ambiguous']

# Change cluster 3 to endothelial for the strict classification and caf for the lenient classification:
sc_data[dbscan_cluster == 3, c('classification', 'classification_lenient') := .('endothelial', 'caf')]

# Change remaining nonmalignant classifications to the authors' classifications, modifying names to my taste:
sc_data[
	classification == 'nonmalignant',
	classification := mapvalues(cell_type, c('fibroblast', 'myeloid'), c('caf', 'macrophage'))
]
sc_data[
	classification_lenient == 'nonmalignant',
	classification_lenient := mapvalues(cell_type, c('fibroblast', 'myeloid'), c('caf', 'macrophage'))
]

# Exclude malignant cells that the authors classified as nonmalignant cell types:
sc_data[classification == 'malignant' & cell_type != 'epithelial', classification := 'ambiguous']
sc_data[classification_lenient == 'malignant' & cell_type != 'epithelial', classification_lenient := 'ambiguous']

# Change "malignant" to "cancer":
sc_data[classification == 'malignant', classification := 'cancer']
sc_data[classification_lenient == 'malignant', classification_lenient := 'cancer']

sc_data[, dbscan_cluster := NULL]

setnames(sc_data, c('cell_type', 'classification', 'classification_lenient'), c('cell_type_author', 'cell_type', 'cell_type_lenient'))

fwrite(sc_data, '../data_and_figures/lee_crc_2020_smc_reclassified.csv')





# hnscc_puram:

sc_data <- fread('../data_and_figures/puram_hnscc_2017.csv')[, -c('lymph_node', 'processed_by_maxima_enzyme')][
	patient %in% c(5, 6, 18, 20, 22, 25, 26)
]
classification_data <- fread(paste0('../data_and_figures/sc_reclassify/classification_data_hnscc_puram.csv'), key = 'cell_id')
sc_data[, classification := classification_data[id, classification]]
setcolorder(sc_data, c('id', 'patient', 'cell_type', 'classification'))

# Exclude unclassified cells and cells that I believe are nonmalignant but that were labelled malignant in the original publication:
sc_data[classification == 'nonmalignant' & cell_type %in% c('', 'cancer'), classification := 'ambiguous']

# Use classifications from the publication for the remaining nonmalignant cells:
sc_data[classification == 'nonmalignant', classification := mapvalues(cell_type, 'fibroblast', 'caf')]

# Exclude cells that I marked as malignant but which are classified as nonmalignant cell types in the publication:
sc_data[classification == 'malignant' & cell_type != 'cancer', classification := 'ambiguous']

# Change "malignant" to "cancer":
sc_data[classification == 'malignant', classification := 'cancer']

setnames(sc_data, c('cell_type', 'classification'), c('cell_type_author', 'cell_type'))

fwrite(sc_data, '../data_and_figures/puram_hnscc_2017_reclassified.csv')





# liver_ma:

sc_data <- fread('../data_and_figures/ma_liver_2019.csv')[, -c('sample', 'disease')][
	patient %in% c('C25', 'C26', 'C46', 'C56', 'C66', 'H37', 'H38', 'H65')
]
classification_data <- fread(paste0('../data_and_figures/sc_reclassify/classification_data_liver_ma.csv'), key = 'cell_id')
sc_data[, classification := classification_data[id, classification]]
sc_dbscan <- readRDS(paste0('../data_and_figures/sc_reclassify/dbscan_liver_ma.rds'))
sc_data[classification == 'nonmalignant', dbscan_cluster := sc_dbscan$cluster]

# Duplicate column for lenient classifications for CAFs (including clusters 2 and 6):
sc_data[, classification_lenient := classification]

setcolorder(sc_data, c('id', 'patient', 'cell_type', 'classification', 'classification_lenient', 'dbscan_cluster'))

# Rename 'unclassified' cells as 'ambiguous'
sc_data[classification == 'nonmalignant' & cell_type == 'unclassified', classification := 'ambiguous']

# For the lenient classifications, change clusters 2 and 6 to CAF and make the remaining 'unclassified' cells 'ambiguous':
sc_data[dbscan_cluster %in% c(2, 6), classification_lenient := 'caf']
sc_data[classification_lenient == 'nonmalignant' & cell_type == 'unclassified', classification_lenient := 'ambiguous']

# Use the authors' classifications for the remaining nonmalignant cells:
sc_data[classification == 'nonmalignant', classification := mapvalues(cell_type, 'fibroblast', 'caf')]
sc_data[classification_lenient == 'nonmalignant', classification_lenient := cell_type]

# Exclude cells that I marked as malignant but which the authors classified as nonmalignant cell types:
sc_data[classification == 'malignant' & cell_type != 'cancer', classification := 'ambiguous']
sc_data[classification_lenient == 'malignant' & cell_type != 'cancer', classification_lenient := 'ambiguous']

# Change "malignant" to "cancer":
sc_data[classification == 'malignant', classification := 'cancer']
sc_data[classification_lenient == 'malignant', classification_lenient := 'cancer']

sc_data[, dbscan_cluster := NULL]

setnames(sc_data, c('cell_type', 'classification', 'classification_lenient'), c('cell_type_author', 'cell_type', 'cell_type_lenient'))

fwrite(sc_data, '../data_and_figures/ma_liver_2019_reclassified.csv')





# luad_kim:

sc_data <- fread('../data_and_figures/kim_luad_2020.csv')[
	patient %in% c('P0006', 'P0008', 'P0018', 'P0019', 'P0020', 'P0025', 'P0028', 'P0030', 'P0031', 'P0034', 'P1006', 'P1028', 'P1049', 'P1058'),
	-c('cell_type_refined', 'cell_subtype', 'sample_site', 'sample_id', 'sample_origin')
]
classification_data <- fread(paste0('../data_and_figures/sc_reclassify/classification_data_luad_kim.csv'), key = 'cell_id')
sc_data[, classification := classification_data[id, classification]]
sc_dbscan <- readRDS(paste0('../data_and_figures/sc_reclassify/dbscan_luad_kim.rds'))
sc_data[classification == 'nonmalignant', dbscan_cluster := sc_dbscan$cluster]
setcolorder(sc_data, c('id', 'patient', 'cell_type', 'classification', 'dbscan_cluster'))

# Change cluster 15 to dendritic cells:
sc_data[dbscan_cluster == 15, classification := 'dendritic']

# Remove the unidentified immune-like cluster:
sc_data[dbscan_cluster == 18, classification := 'ambiguous']

# Use my classifications for the remaining nonmalignant cells that aren't in the T cell or noise clusters - these classifications are nearly the same
# as the authors', but there is enough disagreement to justify using my own definitions:
sc_data[
	classification == 'nonmalignant' & !(dbscan_cluster %in% c(0, 5, 6, 14)),
	classification := mapvalues(
		dbscan_cluster,
		c(1, 8, 11, 19, 9, 16, 7, 13, 2, 17, 3, 10, 12, 20, 4),
		c('b_cell', 'b_cell', 'b_cell', 'b_cell', 'b_cell', 'b_cell', 'caf', 'endothelial', 'epithelial', 'epithelial',
		  'macrophage', 'macrophage', 'macrophage', 'macrophage', 'mast')
	)
]

# For the T cell cluster, use the authors' NK cell classifications but mark all other cells as T cells (using all the authors' classifications
# would put some macrophages, fibroblasts etc. in the T cell cluster):
sc_data[dbscan_cluster %in% c(5, 6, 14), classification := switch((cell_type == 'NK cells') + 1, 't_cell', 'nk_cell'), by = id]

# Use the authors' classifications for the noise cells:
sc_data[
	dbscan_cluster == 0,
	classification := mapvalues(
		cell_type,
		c('B lymphocytes', 'Endothelial cells', 'Epithelial cells', 'Fibroblasts', 'MAST cells', 'Myeloid cells', 'NK cells', 'T lymphocytes'),
		c('b_cell', 'endothelial', 'epithelial', 'caf', 'mast', 'macrophage', 'nk_cell', 't_cell'),
		warn_missing = FALSE
	)
]

# Exclude cells that I marked as malignant but which the authors classified as non-epithelial cell types:
sc_data[classification == 'malignant' & cell_type != 'Epithelial cells', classification := 'ambiguous']

# Change "malignant" to "cancer":
sc_data[classification == 'malignant', classification := 'cancer']

sc_data[, dbscan_cluster := NULL]

setnames(sc_data, c('cell_type', 'classification'), c('cell_type_author', 'cell_type'))

fwrite(sc_data, '../data_and_figures/kim_luad_2020_reclassified.csv')





# lung_qian:

# Keep the 'disease' column so we can separate in to LUSC and LUAD:
sc_data <- fread('../data_and_figures/qian_lung_2020.csv')[sample_type != 'normal', -'sample_type']
classification_data <- fread(paste0('../data_and_figures/sc_reclassify/classification_data_lung_qian.csv'), key = 'cell_id')
sc_data[, classification := classification_data[id, classification]]
sc_dbscan <- readRDS(paste0('../data_and_figures/sc_reclassify/dbscan_lung_qian.rds'))
sc_data[classification == 'nonmalignant', dbscan_cluster := sc_dbscan$cluster]

# Duplicate column for lenient classifications for CAFs (including cluster 15):
sc_data[, classification_lenient := classification]

setcolorder(sc_data, c('id', 'patient', 'disease', 'cell_type', 'classification', 'classification_lenient', 'dbscan_cluster'))

# Exclude cells that the authors denote as 'Cancer' but I think are nonmalignant, along with cells in cluster 8 (unidenfitied immune
# cluster) and cluster 15 (part of group classified as fibroblast by authors but which I think are pericytes):
sc_data[classification == 'nonmalignant' & (cell_type == 'Cancer' | dbscan_cluster %in% c(8, 15)), classification := 'ambiguous']

# Similarly for the lenient classification, but change cluster 15 to CAFs:
sc_data[dbscan_cluster == 15, classification_lenient := 'caf']
sc_data[classification_lenient == 'nonmalignant' & (cell_type == 'Cancer' | dbscan_cluster == 8), classification_lenient := 'ambiguous']

# Change cluster 5 to dendritic cells:
sc_data[dbscan_cluster == 5, c('classification', 'classification_lenient') := .('dendritic', 'dendritic')]

# Change to 'b_cell' those cells which appear in my B cell clusters but which are identified as "Myeloid" by the authors:
sc_data[dbscan_cluster %in% c(13, 18, 4, 19) & cell_type == 'Myeloid', c('classification', 'classification_lenient') := .('b_cell', 'b_cell')]

# Use the authors' classifications for the remaining nonmalignant cells, modifying names to my taste:
sc_data[
	classification == 'nonmalignant',
	classification := mapvalues(
		cell_type,
		c('Alveolar', 'B_cell', 'EC', 'Epithelial', 'Erythroblast', 'Fibroblast', 'Mast_cell', 'Myeloid', 'T_cell'),
		c('alveolar', 'b_cell', 'endothelial', 'epithelial', 'erythroblast', 'caf', 'mast', 'macrophage', 't_cell')
	)
]
sc_data[
	classification_lenient == 'nonmalignant',
	classification_lenient := mapvalues(
		cell_type,
		c('Alveolar', 'B_cell', 'EC', 'Epithelial', 'Erythroblast', 'Fibroblast', 'Mast_cell', 'Myeloid', 'T_cell'),
		c('alveolar', 'b_cell', 'endothelial', 'epithelial', 'erythroblast', 'caf', 'mast', 'macrophage', 't_cell')
	)
]

# Exclude cells that I marked as malignant but which the authors classified as nonmalignant cell types:
sc_data[classification == 'malignant' & cell_type != 'Cancer', classification := 'ambiguous']
sc_data[classification_lenient == 'malignant' & cell_type != 'Cancer', classification_lenient := 'ambiguous']

# Change "malignant" to "cancer":
sc_data[classification == 'malignant', classification := 'cancer']
sc_data[classification_lenient == 'malignant', classification_lenient := 'cancer']

sc_data[, dbscan_cluster := NULL]

setnames(sc_data, c('cell_type', 'classification', 'classification_lenient'), c('cell_type_author', 'cell_type', 'cell_type_lenient'))

fwrite(sc_data, '../data_and_figures/qian_lung_2020_reclassified.csv')





# ovarian_qian:

sc_data <- fread('../data_and_figures/qian_ovarian_2020.csv')[sample_type != 'normal', -c('sample_type', 'site')][patient %in% 11:14]
classification_data <- fread(paste0('../data_and_figures/sc_reclassify/classification_data_ovarian_qian.csv'), key = 'cell_id')
sc_data[, classification := classification_data[id, classification]]
sc_dbscan <- readRDS(paste0('../data_and_figures/sc_reclassify/dbscan_ovarian_qian.rds'))
sc_data[classification == 'nonmalignant', dbscan_cluster := sc_dbscan$cluster]

# Duplicate column for lenient classifications for CAFs (including clusters 2 and 6):
sc_data[, classification_lenient := classification]

setcolorder(sc_data, c('id', 'patient', 'cell_type', 'classification', 'classification_lenient', 'dbscan_cluster'))

# Exclude nonmalignant cells classified by the authors as 'Cancer', along with cluster 9 (unidenfitied immune cells), cluster 8 (looks
# more like pericytes than CAFs) and clusters 15 and 17 (classified as fibroblasts by the authors but appear ambiguous in my analysis):
sc_data[classification == 'nonmalignant' & (cell_type == 'Cancer' | dbscan_cluster %in% c(9, 8, 15, 17)), classification := 'ambiguous']

# Similar for the lenient classification, but change clusters 8, 15 and 17 to CAF:
sc_data[dbscan_cluster %in% c(8, 15, 17), classification_lenient := 'caf']
sc_data[classification_lenient == 'nonmalignant' & (cell_type == 'Cancer' | dbscan_cluster == 9), classification_lenient := 'ambiguous']

# Use the authors' classifications for the remaining nonmalignant cells, modifying names to my taste:
sc_data[
	classification == 'nonmalignant',
	classification := mapvalues(
		cell_type,
		c('B_cell', 'EC', 'Fibroblast', 'Myeloid', 'T_cell'),
		c('b_cell', 'endothelial', 'caf', 'macrophage', 't_cell')
	)
]
sc_data[
	classification_lenient == 'nonmalignant',
	classification_lenient := mapvalues(
		cell_type,
		c('B_cell', 'EC', 'Fibroblast', 'Myeloid', 'T_cell'),
		c('b_cell', 'endothelial', 'caf', 'macrophage', 't_cell')
	)
]

# Exclude cells that I marked as malignant but which the authors classified as nonmalignant cell types:
sc_data[classification == 'malignant' & cell_type != 'Cancer', classification := 'ambiguous']
sc_data[classification_lenient == 'malignant' & cell_type != 'Cancer', classification_lenient := 'ambiguous']

# Change "malignant" to "cancer":
sc_data[classification == 'malignant', classification := 'cancer']
sc_data[classification_lenient == 'malignant', classification_lenient := 'cancer']

sc_data[, dbscan_cluster := NULL]

setnames(sc_data, c('cell_type', 'classification', 'classification_lenient'), c('cell_type_author', 'cell_type', 'cell_type_lenient'))

fwrite(sc_data, '../data_and_figures/qian_ovarian_2020_reclassified.csv')





# pdac_peng:

sc_data <- fread('../data_and_figures/peng_pdac_2019.csv')[sample_type != 'normal', -'sample_type'][
	patient %in% c('T2', 'T3', 'T6', 'T7', 'T8', 'T9', 'T11', 'T13', 'T14', 'T15', 'T16', 'T17', 'T19', 'T20', 'T21', 'T22')
]
classification_data <- fread(paste0('../data_and_figures/sc_reclassify/classification_data_pdac_peng.csv'), key = 'cell_id')
sc_data[, classification := classification_data[id, classification]]
sc_dbscan <- readRDS(paste0('../data_and_figures/sc_reclassify/dbscan_pdac_peng.rds'))
sc_data[classification == 'nonmalignant', dbscan_cluster := sc_dbscan$cluster]
setcolorder(sc_data, c('id', 'patient', 'cell_type', 'classification', 'dbscan_cluster'))

# Change cluster 10 to B cells:
sc_data[dbscan_cluster == 10, classification := 'b_cell']

# Use the authors' classifications for the remaining nonmalignant cells:
sc_data[classification == 'nonmalignant', classification := mapvalues(cell_type, 'fibroblast', 'caf')]

# Exclude cells that I marked as malignant but which the authors classified as nonmalignant cell types:
sc_data[classification == 'malignant' & !(cell_type %in% c('ductal_1', 'ductal_2')), classification := 'ambiguous']

# Change "malignant" to "cancer":
sc_data[classification == 'malignant', classification := 'cancer']

sc_data[, dbscan_cluster := NULL]

setnames(sc_data, c('cell_type', 'classification'), c('cell_type_author', 'cell_type'))

fwrite(sc_data, '../data_and_figures/peng_pdac_2019_reclassified.csv')
