# bsub -q tirosh -n 8 -R "rusage[mem=8000]" -o sc_infercna_crc_lee_kul3_log.o -e sc_infercna_crc_lee_kul3_log.e Rscript sc_infercna_crc_lee_kul3.R

cat(R.Version()$version.string, '\n')

library(data.table) # 1.12.8
library(ggplot2) # 3.3.1
library(magrittr) # 1.5
library(infercna) # 1.0.0

source('../../chr.order.R')
source('general_functions.R')

sc_data <- fread('../data_and_figures/lee_crc_2020_kul3.csv')[, -c('sample_id', 'cell_type', 'cell_subtype')]

# First, using cells from tumour sample as reference:

sc_dbscan <- readRDS('../data_and_figures/dbscan_crc_lee_kul3.rds')

ref_cells <- list(
	b_cell = sc_data[sample_type != 'normal'][sc_dbscan$cluster %in% c(15, 17), id],
	b_plasma = sc_data[sample_type != 'normal'][sc_dbscan$cluster %in% c(16, 18), id],
	endothelial = sc_data[sample_type != 'normal'][sc_dbscan$cluster == 10, id],
	macrophage = sc_data[sample_type != 'normal'][sc_dbscan$cluster == 13, id],
	mast = sc_data[sample_type != 'normal'][sc_dbscan$cluster == 19, id],
	t_cell = sc_data[sample_type != 'normal'][sc_dbscan$cluster == 14, id]
)

gene_averages <- sapply(
	sc_data[sample_type != 'normal', -c('id', 'patient', 'sample_type')],
	function(x) {log2(mean(10*(2^x - 1)) + 1)},
	USE.NAMES = TRUE
)

sc_data_subset <- sc_data[sample_type != 'normal', c('id', 'patient', names(gene_averages)[gene_averages >= 4]), with = FALSE]

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

pdf('../data_and_figures/cna_plots_crc_lee_kul3.pdf')

for(p in names(inferred_cna)) {
	
	cna_scatterplot <- cnaScatterPlot(inferred_cna[[p]], gene.quantile = 0.9, main = p)
	
	# I'm still undecided on whether removing the reference cell is a good idea, though it obviously makes plotting less computationally expensive.
	inferred_cna_subset <- inferred_cna[[p]][, colnames(inferred_cna[[p]])[!(colnames(inferred_cna[[p]]) %in% unlist(ref_cells))]]
	
	# cell_order <- hclust(dist(t(inferred_cna_subset[names(gene_averages)[names(gene_averages) %in% rownames(inferred_cna_subset)][1:2000], ])))$order
	cell_order <- hclust(dist(t(inferred_cna_subset)))$order
	
	# gene_order <- chr.order(rownames(inferred_cna[[p]]))
	
	cna_heatmap <- ggplot(
		reshape2::melt(
			inferred_cna_subset,
			varnames = c('gene', 'cell'),
			value.name = 'cna_score'
		)
	) +
		geom_raster(
			aes(
				x = factor(gene, levels = rownames(inferred_cna_subset)),
				y = factor(cell, levels = colnames(inferred_cna_subset)[cell_order]),
				fill = cna_score
			)
		) +
		scale_fill_gradientn(
			colours = rev(colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))(50)),
			limits = c(-1, 1),
			oob = scales::squish
		) +
		scale_y_discrete(expand = c(0, 0)) +
		scale_x_discrete(expand = c(0, 0)) +
		theme(
			axis.text = element_blank(),
			axis.ticks = element_blank(),
			axis.ticks.length = unit(0, 'pt'),
			panel.border = element_rect(colour = 'black', size = 0.5, fill = NA)
		) +
		labs(x = 'genes', y = 'cells', title = p)
	
	print(cna_heatmap)
	
}

dev.off()

# Next using the adjacent normal samples as reference:

sc_dbscan_normal <- readRDS('../data_and_figures/dbscan_normal_crc_lee_kul3.rds')

# In the following, I'm excluding some clusters that I can't identify, but they look non-malignant because they're not patient-specific.
# Should I include them as 'unknown'?
ref_cells <- list(
	b_cell = sc_data[sample_type == 'normal'][sc_dbscan_normal$cluster == 12, id],
	b_plasma = sc_data[sample_type == 'normal'][sc_dbscan_normal$cluster == 11, id],
	endothelial = sc_data[sample_type == 'normal'][sc_dbscan_normal$cluster == 4, id],
	epithelial = sc_data[sample_type == 'normal'][sc_dbscan_normal$cluster %in% 1:2, id],
	fibroblast = sc_data[sample_type == 'normal'][sc_dbscan_normal$cluster == 3, id],
	macrophage = sc_data[sample_type == 'normal'][sc_dbscan_normal$cluster == 9, id],
	t_mast = sc_data[sample_type == 'normal'][sc_dbscan_normal$cluster == 10, id]
)

gene_averages <- sapply(
	sc_data[, -c('id', 'patient', 'sample_type')],
	function(x) {log2(mean(10*(2^x - 1)) + 1)},
	USE.NAMES = TRUE
)

sc_data_subset <- sc_data[, c('id', 'patient', names(gene_averages)[gene_averages >= 4]), with = FALSE]

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

pdf('../data_and_figures/cna_plots_crc_lee_kul3_normalref.pdf')

for(p in names(inferred_cna)) {
	
	cna_scatterplot <- cnaScatterPlot(inferred_cna[[p]], gene.quantile = 0.9, main = p)
	
	# I'm still undecided on whether removing the reference cell is a good idea, though it obviously makes plotting less computationally expensive.
	inferred_cna_subset <- inferred_cna[[p]][, colnames(inferred_cna[[p]])[!(colnames(inferred_cna[[p]]) %in% unlist(ref_cells))]]
	
	# cell_order <- hclust(dist(t(inferred_cna_subset[names(gene_averages)[names(gene_averages) %in% rownames(inferred_cna_subset)][1:2000], ])))$order
	cell_order <- hclust(dist(t(inferred_cna_subset)))$order
	
	# gene_order <- chr.order(rownames(inferred_cna[[p]]))
	
	cna_heatmap <- ggplot(
		reshape2::melt(
			inferred_cna_subset,
			varnames = c('gene', 'cell'),
			value.name = 'cna_score'
		)
	) +
		geom_raster(
			aes(
				x = factor(gene, levels = rownames(inferred_cna_subset)),
				y = factor(cell, levels = colnames(inferred_cna_subset)[cell_order]),
				fill = cna_score
			)
		) +
		scale_fill_gradientn(
			colours = rev(colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))(50)),
			limits = c(-1, 1),
			oob = scales::squish
		) +
		scale_y_discrete(expand = c(0, 0)) +
		scale_x_discrete(expand = c(0, 0)) +
		theme(
			axis.text = element_blank(),
			axis.ticks = element_blank(),
			axis.ticks.length = unit(0, 'pt'),
			panel.border = element_rect(colour = 'black', size = 0.5, fill = NA)
		) +
		labs(x = 'genes', y = 'cells', title = p)
	
	print(cna_heatmap)
	
}

dev.off()
