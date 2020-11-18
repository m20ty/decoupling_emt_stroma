library(data.table) # 1.12.8
library(ggplot2) # This is version 3.3.1 on my laptop's WSL setup but 3.3.0 on WEXAC.
library(magrittr) # 1.5
library(plyr) # 1.8.6
library(cowplot) # 1.0.0
library(colorspace) # 1.4.1
library(msigdbr) # 7.1.1
library(clusterProfiler) # 3.14.3

source('general_functions.R')

# Get MSigDB C2, C5 and H gene sets (C2 = Curated gene sets; C5 = GO gene sets):
msigdb <- rbindlist(
    lapply(
        c('C2', 'C5', 'H'),
        function(categ) {
            as.data.table(msigdbr(category = categ))
        }
    )
)

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
					.(
						pval = try_default(
							t.test(
								mat_clustpair[[1]][gene_id, ],
								mat_clustpair[[2]][gene_id, ]
							)$p.value,
							default = 1,
							quiet = TRUE
						)
					),
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
	
	gene_test_results <- setNames(
		lapply(gene_test_results, `[[`, 'data'),
		sapply(gene_test_results, `[[`, 'other_clst')
	)
	
	table_list <- lapply(gene_test_results, function(dt) dt[p.adjust(pval, method ='BH') < 0.05])
	for(dt in table_list) {setkey(dt, gene_id)}
	
	degenes <- data.table(gene = Reduce(intersect, lapply(table_list, `[[`, 'gene_id')))
	# degenes <- data.table(gene = intersect(table_list[[1]]$gene_id, table_list[[2]]$gene_id))
	degenes[
		,
		rel_exp := apply(
			as.data.table(lapply(table_list, function(dt) dt[gene, rel_exp])),
			# cbind(table_list[[1]][gene, rel_exp], table_list[[2]][gene, rel_exp]),
			1,
			mean
		)
	]
	
	# Run GSEA (note GSEA involves random permutations, I think in correcting for multiple hypothesis tests, so we need to set a seed):
	
	if(!is.null(seed)) set.seed(seed)
	
	degenes_gsea <- sapply(
		unique(msigdb_table$gs_cat),
		function(categ) {
			cbind(
				gs_cat = categ,
				as.data.table(
					GSEA(
						degenes[
							order(-rel_exp),
							setNames(rel_exp, gene)
						],
						TERM2GENE = msigdb_table[gs_cat == categ, .(gs_name, human_gene_symbol)]
					)
				)
			)
		},
		simplify = FALSE,
		USE.NAMES = TRUE
	)
	
	gsea_table <- rbindlist(
		lapply(
			degenes_gsea,
			function(dt) dt[, Reduce(intersect, lapply(degenes_gsea, names)), with = FALSE]
		)
	)
	
	list(test_results = gene_test_results, degenes = degenes, gsea_table = gsea_table)
	
}

cell_type_labels <- c(
	'acinar' = 'Acinar',
	'alveolar' = 'Alveolar',
	'b_cell' = 'B cell',
	'caf' = 'CAF',
	'caf_1' = 'Potential CAF 1',
	'caf_2' = 'Potential CAF 2',
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
	'stellate' = 'Stellate',
	't_cell' = 'T cell',
	'unidentified' = 'Unidentified'
)

cell_type_colours <- c(
	'acinar' = '#C0DF84',
	'alveolar' = '#E3E050',
	'b_cell' = '#887DDA',
	'caf' = '#D8A354',
	'caf_1' = lighten('#F8A354', 0.5),
	'caf_2' = lighten('#D8BB54', 0.5),
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
	'stellate' = '#D5978A',
	't_cell' = '#79E5A3',
	'unidentified' = 'grey20'
)





tsne_data <- list()
heatmap_data <- list()
heatmap_annotations <- list()





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
	cbind(sc_data$id, as.data.table(sc_tsne$Y), as.character(sc_dbscan$cluster), sc_data$cell_type),
	c('cell_id', 'x', 'y', 'cluster', 'cell_type')
)

plot_data <- plot_data[cluster != 0]

plot_data[
	,
	c('cluster', 'cell_type') := .(
		mapvalues(
			cluster,
			c(8, 4, 6, 3, 9, 1, 7, 10, 2, 5, 11),
			c('b_cell', 'b_cell', 'endothelial', 'macrophage', 'mast', 't_cell', 't_cell', 't_cell', 'caf', 'caf_potential', 'unidentified')
		),
		mapvalues(
			cell_type,
			c('B_cell', 'Cancer', 'DC', 'EC', 'Fibroblast', 'Mast', 'Myeloid', 'T_cell'),
			c('b_cell', 'cancer', 'dendritic', 'endothelial', 'caf', 'mast', 'macrophage', 't_cell')
		)
	)
]

tsne_data$breast_qian <- plot_data

# The following is out of date because I updated gene names in the single cell data:
# de_between_caf_top_50 <- readRDS('../data_and_figures/sc_reclassify/files_for_rmd/breast_qian_de_between_caf_top_50.rds')

# Difference between the CAF clusters (taken from the sc_data_cell_type_classification.R file):

seeds <- c(6365, 878)
setkey(sc_data, id)

de_between_caf <- lapply(
	c('caf', 'caf_potential'),
	function(x) {
		test_cluster(
			set_rownames(
				as.matrix(sc_data[plot_data[cluster %in% c('caf', 'caf_potential'), cell_id], -c('id', 'patient', 'cell_type')]),
				plot_data[cluster %in% c('caf', 'caf_potential'), cell_id]
			),
			plot_data[cluster %in% c('caf', 'caf_potential')],
			x,
			msigdb,
			seed = seeds[c('caf', 'caf_potential') == x]
		)
	}
)

de_between_caf_top_50 <- set_colnames(sapply(1:2, function(i) de_between_caf[[i]]$degenes[order(-rel_exp), head(gene, 50)]), c('CAF', 'Potential CAF'))

de_data <- rbindlist(
	lapply(
		c('caf', 'caf_potential'),
		function(ct) {
			dt <- cbind(
				sc_data[
					id %in% plot_data[cluster == ct, cell_id],
					c('id', de_between_caf_top_50[, 1], rev(de_between_caf_top_50[, 2])),
					with = FALSE
				],
				cell_type = ct
			)
			setcolorder(dt, c('id', 'cell_type'))
			dt
		}
	)
)

de_data[
	,
	names(de_data[, -c('id', 'cell_type')]) := lapply(
		.SD,
		function(x) {x - mean(c(mean(x[cell_type == 'caf']), mean(x[cell_type == 'caf_potential'])))}
	),
	.SDcols = -c('id', 'cell_type')
]

de_data[
	,
	diff_exp := mean(as.numeric(.SD[, de_between_caf_top_50[, 1], with = FALSE])) -
		mean(as.numeric(.SD[, de_between_caf_top_50[, 2], with = FALSE])),
	by = id
]

setcolorder(de_data, c('id', 'cell_type', 'diff_exp'))

de_data[, r := order(-diff_exp), by = cell_type]
de_data[cell_type == 'caf_potential', r := r + de_data[cell_type == 'caf', .N]]
de_data <- de_data[de_data$r]
de_data[, r := NULL]

heatmap_data$breast_qian <- de_data

de_genes_for_annotations <- names(de_data[, -c('id', 'cell_type', 'diff_exp')])
de_genes_for_annotations[
	!(
		de_genes_for_annotations %in% c('MMP2', 'DCN', 'COL1A1', 'VCAN', 'FN1', 'COL1A2', 'COL3A1',
										'RGS5', 'COL4A1', 'MGP', 'ESAM', 'MYO1B', 'ADGRF5', 'MEF2C')
	)
] <- ''

heatmap_annotations$breast_qian <- de_genes_for_annotations

# Clean up workspace:
rm(classification_data)
rm(sc_data)
rm(gene_averages)
rm(sc_tsne)
rm(sc_dbscan)





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
	cbind(sc_data$id, as.data.table(sc_tsne$Y), as.character(sc_dbscan$cluster), sc_data$cell_type),
	c('cell_id', 'x', 'y', 'cluster', 'cell_type')
)

plot_data <- plot_data[cluster != 0]

plot_data[
	,
	c('cluster', 'cell_type') := .(
		mapvalues(
			cluster,
			c(9, 8, 2, 3, 1, 4, 6, 7, 5),
			c('b_cell', 'b_cell', 'caf', 'endothelial', 'epithelial', 'macrophage', 't_cell', 't_cell', 'unidentified')
		),
		mapvalues(cell_type, c('fibroblast', 'myeloid'), c('caf', 'macrophage'))
	)
]

tsne_data$crc_lee_smc <- plot_data

# de_between_caf_endothelial_top_50 <- readRDS('../data_and_figures/sc_reclassify/files_for_rmd/crc_lee_smc_de_between_caf_endothelial_top_50.rds')

seeds <- c(2881, 9041)
setkey(sc_data, id)

de_between_caf_endothelial <- lapply(
	c('caf', 'endothelial'),
	function(x) {
		test_cluster(
			set_rownames(
				as.matrix(sc_data[plot_data[cluster %in% c('caf', 'endothelial'), cell_id], -c('id', 'patient', 'cell_type')]),
				plot_data[cluster %in% c('caf', 'endothelial'), cell_id]
			),
			plot_data[cluster %in% c('caf', 'endothelial')],
			x,
			msigdb,
			seed = seeds[c('caf', 'endothelial') == x]
		)
	}
)

de_between_caf_endothelial_top_50 <- set_colnames(
	sapply(1:2, function(i) de_between_caf_endothelial[[i]]$degenes[order(-rel_exp), head(gene, 50)]),
	c('CAF', 'Endothelial')
)

de_data <- rbindlist(
	lapply(
		c('caf', 'endothelial'),
		function(ct) {
			dt <- cbind(
				sc_data[
					id %in% plot_data[cluster == ct, cell_id],
					c('id', de_between_caf_endothelial_top_50[, 1], rev(de_between_caf_endothelial_top_50[, 2])),
					with = FALSE
				],
				cell_type = ct
			)
			setcolorder(dt, c('id', 'cell_type'))
			dt
		}
	)
)

de_data[
	,
	names(de_data[, -c('id', 'cell_type')]) := lapply(
		.SD,
		function(x) {x - mean(c(mean(x[cell_type == 'caf']), mean(x[cell_type == 'endothelial'])))}
	),
	.SDcols = -c('id', 'cell_type')
]

de_data[
	,
	diff_exp := mean(as.numeric(.SD[, de_between_caf_endothelial_top_50[, 1], with = FALSE])) -
		mean(as.numeric(.SD[, de_between_caf_endothelial_top_50[, 2], with = FALSE])),
	by = id
]

setcolorder(de_data, c('id', 'cell_type', 'diff_exp'))

de_data[, r := order(-diff_exp), by = cell_type]
de_data[cell_type == 'endothelial', r := r + de_data[cell_type == 'caf', .N]]
de_data <- de_data[de_data$r]
de_data[, r := NULL]

heatmap_data$crc_lee_smc <- de_data

de_genes_for_annotations <- names(de_data[, -c('id', 'cell_type', 'diff_exp')])
de_genes_for_annotations[
	!(
		de_genes_for_annotations %in% c('COL1A1', 'COL1A2', 'COL3A1', 'TAGLN', 'ACTA2', 'MYL9', 'THY1', 'VCAN',
										'RAMP2', 'PECAM1', 'RAMP3', 'EGFL7', 'ADGRL4', 'ENG', 'FLT1', 'SELE', 'ESAM')
	)
] <- ''

heatmap_annotations$crc_lee_smc <- de_genes_for_annotations

# Clean up workspace:
rm(classification_data)
rm(sc_data)
rm(gene_averages)
rm(sc_tsne)
rm(sc_dbscan)





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
	cbind(sc_data$id, as.data.table(sc_tsne$Y), as.character(sc_dbscan$cluster), sc_data$cell_type),
	c('cell_id', 'x', 'y', 'cluster', 'cell_type')
)

plot_data[
	,
	c('cluster', 'cell_type') := .(
		mapvalues(
			cluster,
			c(1, 4, 3, 5, 6),
			c('macrophage', 'macrophage', 't_cell', 'endothelial', 'caf_potential')
		),
		mapvalues(cell_type, c('fibroblast', 'unclassified'), c('caf', 'unidentified'))
	)
]

plot_data[cluster == 2 & sc_data$cell_type != 'fibroblast', cluster := 'caf_potential']
plot_data[cluster == 2 & sc_data$cell_type == 'fibroblast', cluster := 'caf']

plot_data <- plot_data[cluster != 0]

tsne_data$liver_ma <- plot_data

seeds <- c(7152, 7920)
setkey(sc_data, id)

de_between <- lapply(
	c('caf', 'caf_potential'),
	function(x) {
		test_cluster(
			set_rownames(
				as.matrix(sc_data[plot_data[cluster %in% c('caf', 'caf_potential'), cell_id], -c('id', 'patient', 'cell_type')]),
				plot_data[cluster %in% c('caf', 'caf_potential'), cell_id]
			),
			plot_data[cluster %in% c('caf', 'caf_potential')],
			x,
			msigdb,
			seed = seeds[c('caf', 'caf_potential') == x]
		)
	}
)

de_between_top_50 <- set_colnames(
	sapply(1:2, function(i) de_between[[i]]$degenes[order(-rel_exp), head(gene, 50)]),
	c('CAF', 'Potential CAF')
)

# This actually makes it pretty obvious that the potential CAFs are endothelial, though I guess they could be doublets.

de_data <- rbindlist(
	lapply(
		c('caf', 'caf_potential'),
		function(ct) {
			dt <- cbind(
				sc_data[
					id %in% plot_data[cluster == ct, cell_id],
					c('id', de_between_top_50[, 1], rev(de_between_top_50[, 2])),
					with = FALSE
				],
				cell_type = ct
			)
			setcolorder(dt, c('id', 'cell_type'))
			dt
		}
	)
)

de_data[
	,
	names(de_data[, -c('id', 'cell_type')]) := lapply(
		.SD,
		function(x) {x - mean(c(mean(x[cell_type == 'caf']), mean(x[cell_type == 'caf_potential'])))}
	),
	.SDcols = -c('id', 'cell_type')
]

de_data[
	,
	diff_exp := mean(as.numeric(.SD[, de_between_top_50[, 1], with = FALSE])) -
		mean(as.numeric(.SD[, de_between_top_50[, 2], with = FALSE])),
	by = id
]

setcolorder(de_data, c('id', 'cell_type', 'diff_exp'))

de_data[, r := order(-diff_exp), by = cell_type]
de_data[cell_type == 'caf_potential', r := r + de_data[cell_type == 'caf', .N]]
de_data <- de_data[de_data$r]
de_data[, r := NULL]

heatmap_data$liver_ma <- de_data

de_genes_for_annotations <- names(de_data[, -c('id', 'cell_type', 'diff_exp')])
de_genes_for_annotations[
	!(
		de_genes_for_annotations %in% c('ACTA2', 'TAGLN', 'MYL9', 'DCN', 'MYLK', 'COL1A1', 'COL1A2', 'COL3A1',
										'RAMP2', 'EGFL7', 'PECAM1', 'FLT1', 'ENG', 'RAMP3', 'CDH5')
	)
] <- ''

heatmap_annotations$liver_ma <- de_genes_for_annotations

# Clean up workspace:
rm(classification_data)
rm(sc_data)
rm(gene_averages)
rm(sc_tsne)
rm(sc_dbscan)





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
	cbind(sc_data$id, as.data.table(sc_tsne$Y), as.character(sc_dbscan$cluster), sc_data$cell_type),
	c('cell_id', 'x', 'y', 'cluster', 'cell_type')
)

plot_data <- plot_data[cluster != 0]

plot_data[
	,
	c('cluster', 'cell_type') := .(
		mapvalues(
			cluster,
			c(8, 13, 18, 4, 19, 3, 14, 15, 5, 12, 6, 10, 16, 2, 9, 11, 17, 21, 1, 7, 20, 22),
			c('unidentified', 'b_cell', 'b_cell', 'b_cell', 'b_cell', 'caf', 'caf', 'caf_potential', 'dendritic', 'endothelial', 'epithelial',
			  'epithelial', 'epithelial', 'macrophage', 'macrophage', 'macrophage', 'mast', 'mast', 't_cell', 't_cell', 't_cell', 't_cell')
		),
		mapvalues(
			cell_type,
			c('Alveolar', 'B_cell', 'Cancer', 'EC', 'Epithelial', 'Erythroblast', 'Fibroblast', 'Mast_cell', 'Myeloid', 'T_cell'),
			c('alveolar', 'b_cell', 'cancer', 'endothelial', 'epithelial', 'erythroblast', 'caf', 'mast', 'macrophage', 't_cell')
		)
	)
]

tsne_data$lung_qian <- plot_data

seeds <- c(3358, 1326)
setkey(sc_data, id)

de_between <- lapply(
	c('caf', 'caf_potential'),
	function(x) {
		test_cluster(
			set_rownames(
				as.matrix(sc_data[plot_data[cluster %in% c('caf', 'caf_potential'), cell_id], -c('id', 'patient', 'cell_type')]),
				plot_data[cluster %in% c('caf', 'caf_potential'), cell_id]
			),
			plot_data[cluster %in% c('caf', 'caf_potential')],
			x,
			msigdb,
			seed = seeds[c('caf', 'caf_potential') == x]
		)
	}
)

de_between_top_50 <- set_colnames(
	sapply(1:2, function(i) de_between[[i]]$degenes[order(-rel_exp), head(gene, 50)]),
	c('CAF', 'Potential CAF')
)

de_data <- rbindlist(
	lapply(
		c('caf', 'caf_potential'),
		function(ct) {
			dt <- cbind(
				sc_data[
					id %in% plot_data[cluster == ct, cell_id],
					c('id', de_between_top_50[, 1], rev(de_between_top_50[, 2])),
					with = FALSE
				],
				cell_type = ct
			)
			setcolorder(dt, c('id', 'cell_type'))
			dt
		}
	)
)

de_data[
	,
	names(de_data[, -c('id', 'cell_type')]) := lapply(
		.SD,
		function(x) {x - mean(c(mean(x[cell_type == 'caf']), mean(x[cell_type == 'caf_potential'])))}
	),
	.SDcols = -c('id', 'cell_type')
]

de_data[
	,
	diff_exp := mean(as.numeric(.SD[, de_between_top_50[, 1], with = FALSE])) -
		mean(as.numeric(.SD[, de_between_top_50[, 2], with = FALSE])),
	by = id
]

setcolorder(de_data, c('id', 'cell_type', 'diff_exp'))

de_data[, r := order(-diff_exp), by = cell_type]
de_data[cell_type == 'caf_potential', r := r + de_data[cell_type == 'caf', .N]]
de_data <- de_data[de_data$r]
de_data[, r := NULL]

heatmap_data$lung_qian <- de_data

de_genes_for_annotations <- names(de_data[, -c('id', 'cell_type', 'diff_exp')])
de_genes_for_annotations[
	!(
		de_genes_for_annotations %in% c('MMP2', 'DCN', 'COL1A1', 'VCAN', 'COL3A1', 'EFEMP1', 'PDPN',
										'RGS5', 'COL4A1', 'PDGFRB', 'COL4A2', 'MEF2C', 'ESAM', 'ADGRF5')
	)
] <- ''

heatmap_annotations$lung_qian <- de_genes_for_annotations

# Clean up workspace:
rm(classification_data)
rm(sc_data)
rm(gene_averages)
rm(sc_tsne)
rm(sc_dbscan)





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
	cbind(sc_data[, .(id, patient)], as.data.table(sc_tsne$Y), as.character(sc_dbscan$cluster), sc_data$cell_type),
	c('cell_id', 'patient', 'x', 'y', 'cluster', 'cell_type')
)

plot_data <- plot_data[cluster != 0]

plot_data[
	,
	c('cluster', 'cell_type') := .(
		mapvalues(
			cluster,
			c(9, 13, 5, 14, 6, 3, 19, 21, 7, 12, 11, 16, 1, 2, 4, 10, 15, 17, 18, 20, 8),
			c('unidentified', 'b_cell', 'endothelial', 'endothelial', 'endothelial', 'macrophage', 'macrophage', 'macrophage', 't_cell',
			  't_cell', 'epithelial', 'epithelial', 'caf', 'caf', 'caf', 'caf', 'caf_1', 'caf_1', 'caf', 'caf', 'caf_2')
		),
		mapvalues(
			cell_type,
			c('B_cell', 'Cancer', 'EC', 'Fibroblast', 'Myeloid', 'T_cell'),
			c('b_cell', 'cancer', 'endothelial', 'caf', 'macrophage', 't_cell')
		)
	)
]

tsne_data$ovarian_qian <- plot_data

seeds <- c(9064, 7002, 7175)
setkey(sc_data, id)

de_between <- lapply(
	c('caf', 'caf_1', 'caf_2'),
	function(x) {
		test_cluster(
			set_rownames(
				as.matrix(sc_data[plot_data[cluster %in% c('caf', 'caf_1', 'caf_2'), cell_id], -c('id', 'patient', 'cell_type')]),
				plot_data[cluster %in% c('caf', 'caf_1', 'caf_2'), cell_id]
			),
			plot_data[cluster %in% c('caf', 'caf_1', 'caf_2')],
			x,
			msigdb,
			seed = seeds[c('caf', 'caf_1', 'caf_2') == x]
		)
	}
)

de_between_top_50 <- set_colnames(
	sapply(1:3, function(i) de_between[[i]]$degenes[order(-rel_exp), head(gene, 50)]),
	c('CAF', 'Potential CAF 1', 'Potential CAF 2')
)

de_genes_ovarian_qian <- de_between_top_50

nondupl_genes <- c(de_between_top_50[, 1], de_between_top_50[, 2], de_between_top_50[, 3])
nondupl_genes <- nondupl_genes[!(nondupl_genes %in% names(table(nondupl_genes)[table(nondupl_genes) > 1]))]

de_data <- rbindlist(
	lapply(
		c('caf', 'caf_1', 'caf_2'),
		function(ct) {
			dt <- cbind(
				sc_data[
					id %in% plot_data[cluster == ct, cell_id],
					unique(c('id', nondupl_genes)),
					with = FALSE
				],
				cell_type = ct
			)
			setcolorder(dt, c('id', 'cell_type'))
			dt
		}
	)
)

de_data[
	,
	names(de_data[, -c('id', 'cell_type')]) := lapply(
		.SD,
		function(x) {x - mean(c(mean(x[cell_type == 'caf']), mean(x[cell_type == 'caf_1']), mean(x[cell_type == 'caf_2'])))}
	),
	.SDcols = -c('id', 'cell_type')
]

de_data[
	,
	c('mean_caf', 'mean_caf_1', 'mean_caf_2') := .(
		mean(as.numeric(.SD[, de_between_top_50[, 1][de_between_top_50[, 1] %in% nondupl_genes], with = FALSE])),
		mean(as.numeric(.SD[, de_between_top_50[, 2][de_between_top_50[, 2] %in% nondupl_genes], with = FALSE])),
		mean(as.numeric(.SD[, de_between_top_50[, 3][de_between_top_50[, 3] %in% nondupl_genes], with = FALSE]))
	),
	by = id
]

setcolorder(de_data, c('id', 'cell_type', 'mean_caf', 'mean_caf_1', 'mean_caf_2'))

de_data[, r := order(-get(paste0('mean_', unique(cell_type)))), by = cell_type]
de_data[cell_type == 'caf_1', r := r + de_data[cell_type == 'caf', .N]]
de_data[cell_type == 'caf_2', r := r + de_data[cell_type %in% c('caf', 'caf_1'), .N]]
de_data <- de_data[de_data$r]
de_data[, r := NULL]

heatmap_data$ovarian_qian <- de_data

de_genes_for_annotations <- names(de_data[, -c('id', 'cell_type', 'mean_caf', 'mean_caf_1', 'mean_caf_2')])
de_genes_for_annotations[
	!(
		de_genes_for_annotations %in% c('COL3A1', 'COL1A1', 'VCAN', 'COLA2', 'FN1', 'MMP2',
										'STAR', 'NR4A2', 'NR4A1', 'JUNB', 'NFKBIA', 'BTG2',
										'RGS5', 'ACTA2', 'MEF2C', 'PDGFRB', 'MYL9', 'MYLK', 'ESAM')
	)
] <- ''

heatmap_annotations$ovarian_qian <- de_genes_for_annotations

# Clean up workspace:
rm(classification_data)
rm(sc_data)
rm(gene_averages)
rm(sc_tsne)
rm(sc_dbscan)





tsne_plots <- sapply(
	names(tsne_data),
	function(cohort) {
		list(
			dbscan = ggplot(tsne_data[[cohort]], aes(x = x, y = y, colour = cluster)) +
				geom_point(size = 0.7) +
				scale_colour_manual(labels = cell_type_labels, values = cell_type_colours) +
				theme_minimal() +
				labs(x = 't-SNE 1', y = 't-SNE 2', title = 'Manual classification', colour = 'Cell type'),
			author_cell_types = ggplot(tsne_data[[cohort]], aes(x = x, y = y, colour = cell_type)) +
				geom_point(size = 0.7) +
				scale_colour_manual(labels = cell_type_labels, values = cell_type_colours) +
				theme_minimal() +
				labs(x = 't-SNE 1', y = 't-SNE 2', title = 'Original classifications', colour = 'Cell type')
		)
	},
	simplify = FALSE,
	USE.NAMES = TRUE
)

de_annotations <- sapply(
	names(heatmap_annotations),
	function(cohort) heat_map_labels_repel(heatmap_annotations[[cohort]], edge = 'right', nudge = 0.15),
	simplify = FALSE,
	USE.NAMES = TRUE
)

de_heatmaps <- list()

de_heatmaps$breast_qian <- ggplot(
	melt(heatmap_data$breast_qian, id.vars = c('id', 'cell_type', 'diff_exp'), variable.name = 'gene', value.name = 'expression_level'),
	aes(x = factor(id, levels = unique(id)), y = factor(gene, levels = unique(gene)), fill = expression_level)
) +
	geom_raster() +
	geom_hline(yintercept = 50.5) +
	geom_vline(xintercept = heatmap_data$breast_qian[cell_type == 'caf', .N + 0.5]) +
	scale_x_discrete(
		labels = c('CAFs', 'Potential\nCAFs'),
		breaks = heatmap_data$breast_qian$id[heatmap_data$breast_qian[, .(N = .N), by = cell_type][, c(N[1]/2, N[1] + N[2]/2)]],
		expand = c(0, 0)
	) +
	scale_y_discrete(expand = c(0, 0)) +
	scale_fill_gradientn(
		colours = rev(colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))(50)),
		limits = c(-6, 6),
		oob = scales::squish
	) +
	theme(
		axis.text.y = element_blank(),
		axis.ticks = element_blank(),
		axis.ticks.length.y = unit(0, 'pt'),
		axis.title = element_blank(),
		panel.border = element_rect(fill = NA),
		plot.margin = unit(c(5.5, 5.5, 5.5, 0), 'pt')
	) +
	labs(fill = 'Relative\nexpression')

de_heatmaps$crc_lee_smc <- ggplot(
	melt(heatmap_data$crc_lee_smc, id.vars = c('id', 'cell_type', 'diff_exp'), variable.name = 'gene', value.name = 'expression_level'),
	aes(x = factor(id, levels = unique(id)), y = factor(gene, levels = unique(gene)), fill = expression_level)
) +
	geom_raster() +
	geom_hline(yintercept = 50.5) +
	geom_vline(xintercept = heatmap_data$crc_lee_smc[cell_type == 'caf', .N + 0.5]) +
	scale_x_discrete(
		labels = c('CAFs', 'Endothelial\ncells'),
		breaks = heatmap_data$crc_lee_smc$id[heatmap_data$crc_lee_smc[, .(N = .N), by = cell_type][, c(N[1]/2, N[1] + N[2]/2)]],
		expand = c(0, 0)
	) +
	scale_y_discrete(expand = c(0, 0)) +
	scale_fill_gradientn(
		colours = rev(colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))(50)),
		limits = c(-6, 6),
		oob = scales::squish
	) +
	theme(
		axis.text.y = element_blank(),
		axis.ticks = element_blank(),
		axis.ticks.length.y = unit(0, 'pt'),
		axis.title = element_blank(),
		panel.border = element_rect(fill = NA),
		plot.margin = unit(c(5.5, 5.5, 5.5, 0), 'pt')
	) +
	labs(fill = 'Relative\nexpression')

de_heatmaps$liver_ma <- ggplot(
	melt(heatmap_data$liver_ma, id.vars = c('id', 'cell_type', 'diff_exp'), variable.name = 'gene', value.name = 'expression_level'),
	aes(x = factor(id, levels = unique(id)), y = factor(gene, levels = unique(gene)), fill = expression_level)
) +
	geom_raster() +
	geom_hline(yintercept = 50.5) +
	geom_vline(xintercept = heatmap_data$liver_ma[cell_type == 'caf', .N + 0.5]) +
	scale_x_discrete(
		labels = c('CAFs', 'Potential\nCAFs'),
		breaks = heatmap_data$liver_ma$id[heatmap_data$liver_ma[, .(N = .N), by = cell_type][, c(N[1]/2, N[1] + N[2]/2)]],
		expand = c(0, 0)
	) +
	scale_y_discrete(expand = c(0, 0)) +
	scale_fill_gradientn(
		colours = rev(colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))(50)),
		limits = c(-6, 6),
		oob = scales::squish
	) +
	theme(
		axis.text.y = element_blank(),
		axis.ticks = element_blank(),
		axis.ticks.length.y = unit(0, 'pt'),
		axis.title = element_blank(),
		panel.border = element_rect(fill = NA),
		plot.margin = unit(c(5.5, 5.5, 5.5, 0), 'pt')
	) +
	labs(fill = 'Relative\nexpression')

de_heatmaps$lung_qian <- ggplot(
	melt(heatmap_data$lung_qian, id.vars = c('id', 'cell_type', 'diff_exp'), variable.name = 'gene', value.name = 'expression_level'),
	aes(x = factor(id, levels = unique(id)), y = factor(gene, levels = unique(gene)), fill = expression_level)
) +
	geom_raster() +
	geom_hline(yintercept = 50.5) +
	geom_vline(xintercept = heatmap_data$lung_qian[cell_type == 'caf', .N + 0.5]) +
	scale_x_discrete(
		labels = c('CAFs', 'Potential\nCAFs'),
		breaks = heatmap_data$lung_qian$id[heatmap_data$lung_qian[, .(N = .N), by = cell_type][, c(N[1]/2, N[1] + N[2]/2)]],
		expand = c(0, 0)
	) +
	scale_y_discrete(expand = c(0, 0)) +
	scale_fill_gradientn(
		colours = rev(colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))(50)),
		limits = c(-6, 6),
		oob = scales::squish
	) +
	theme(
		axis.text.y = element_blank(),
		axis.ticks = element_blank(),
		axis.ticks.length.y = unit(0, 'pt'),
		axis.title = element_blank(),
		panel.border = element_rect(fill = NA),
		plot.margin = unit(c(5.5, 5.5, 5.5, 0), 'pt')
	) +
	labs(fill = 'Relative\nexpression')

de_heatmaps$ovarian_qian <- ggplot(
	melt(heatmap_data$ovarian_qian, id.vars = c('id', 'cell_type', 'mean_caf', 'mean_caf_1', 'mean_caf_2'), variable.name = 'gene', value.name = 'expression_level'),
	aes(x = factor(id, levels = unique(id)), y = factor(gene, levels = unique(gene)), fill = expression_level)
) +
	geom_raster() +
	geom_hline(yintercept = length(de_genes_ovarian_qian[, 1][de_genes_ovarian_qian[, 1] %in% nondupl_genes]) + 0.5) +
	geom_hline(
		yintercept = length(de_genes_ovarian_qian[, 1][de_genes_ovarian_qian[, 1] %in% nondupl_genes]) +
			length(de_genes_ovarian_qian[, 2][de_genes_ovarian_qian[, 2] %in% nondupl_genes]) +
			0.5
	) +
	geom_vline(xintercept = heatmap_data$ovarian_qian[cell_type == 'caf', .N + 0.5]) +
	geom_vline(xintercept = heatmap_data$ovarian_qian[cell_type %in% c('caf', 'caf_1'), .N + 0.5]) +
	scale_x_discrete(
		labels = c('CAFs', 'Potential\nCAF 1', 'Potential\nCAF 2'),
		breaks = heatmap_data$ovarian_qian$id[heatmap_data$ovarian_qian[, .(N = .N), by = cell_type][, c(N[1]/2, N[1] + N[2]/2, N[1] + N[2] + N[3]/2)]],
		expand = c(0, 0)
	) +
	scale_y_discrete(expand = c(0, 0)) +
	scale_fill_gradientn(
		colours = rev(colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))(50)),
		limits = c(-6, 6),
		oob = scales::squish
	) +
	theme(
		axis.text.y = element_blank(),
		axis.ticks = element_blank(),
		axis.ticks.length.y = unit(0, 'pt'),
		axis.title = element_blank(),
		panel.border = element_rect(fill = NA),
		plot.margin = unit(c(5.5, 5.5, 5.5, 0), 'pt')
	) +
	labs(fill = 'Relative\nexpression')





pdf('../data_and_figures/final_figures_resubmission/S3.pdf', width = 14, height = 20)
plot_grid(
	plotlist = unlist(
		lapply(
			names(tsne_data),
			function(cohort) c(
				list(
					blank_plot() +
						labs(
							title = mapvalues(
								cohort,
								c('breast_qian', 'crc_lee_smc', 'liver_ma', 'lung_qian', 'ovarian_qian'),
								c('Breast - Qian et al.', 'Colorectal - Lee et al. - SMC cohort', 'Liver - Ma et al.', 'Lung - Qian et al.', 'Ovarian - Qian et al.'),
								warn_missing = FALSE
							)
						) +
						theme(plot.title = element_text(margin = margin(t = 30), size = 18), plot.margin = unit(c(5.5, 5.5, 10, 5.5), 'pt'))
				),
				list(
					plot_grid(
						plot_grid(
							plotlist = c(
								lapply(tsne_plots[[cohort]], function(g) g + theme(plot.margin = unit(c(5.5, 40, 5.5, 5.5), 'pt'), legend.position = 'none')),
								list(
									get_legend(
										ggplot(tsne_data[[cohort]][, .(x = 0, y = 0, clust = unique(c(cluster, cell_type)))], aes(x = x, y = y, colour = clust)) +
											geom_point() +
											scale_colour_manual(labels = cell_type_labels, values = cell_type_colours) +
											theme_minimal() +
											theme(legend.justification = 'left') +
											labs(colour = 'Cell type')
									)
								)
							),
							nrow = 1,
							ncol = 3,
							rel_widths = c(6.75, 6.75, 2.5)
						),
						plot_grid(
							de_annotations[[cohort]],
							de_heatmaps[[cohort]],
							nrow = 1,
							ncol = 2,
							align = 'h',
							rel_widths = c(2, 7)
						),
						nrow = 1,
						ncol = 2,
						rel_widths = c(16, 9)
					)
				)
			)
		),
		recursive = FALSE
	),
	nrow = 10,
	ncol = 1,
	# align = 'h',
	rel_heights = c(1, 4, 1, 4, 1, 4, 1, 4, 1, 4)
)
dev.off()
