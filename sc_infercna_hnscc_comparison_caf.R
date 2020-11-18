# This script is for inferring CNAs using the cells annotated as non-malignant in the original publication as reference.
# This is intended to be a sort of sanity check.
# We're also including CAFs here, as we suspect they may have a big impact on the CNV signal.

# bsub -q tirosh -n 8 -R "rusage[mem=16000]" -o sc_infercna_hnscc_log.o -e sc_infercna_hnscc_log.e Rscript sc_infercna_hnscc.R

cat(R.Version()$version.string, '\n')

library(data.table) # 1.12.8
library(ggplot2) # 3.3.1
library(magrittr) # 1.5
library(infercna) # 1.0.0
library(biomaRt) # 2.42.1

source('../../chr.order.R')
source('general_functions.R')

sc_data <- fread('../data_and_figures/puram_hnscc_2017.csv')[, -c('lymph_node', 'processed_by_maxima_enzyme')]

ref_cells <- list(
	b_cell = sc_data[cell_type == 'b_cell', id],
	fibroblast = sc_data[cell_type == 'fibroblast', id],
	dendritic = sc_data[cell_type == 'dendritic', id],
	endothelial = sc_data[cell_type == 'endothelial', id],
	macrophage = sc_data[cell_type == 'macrophage', id],
	mast = sc_data[cell_type == 'mast', id],
	myocyte = sc_data[cell_type == 'myocyte', id],
	t_cell = sc_data[cell_type == 't_cell', id]
)

gene_averages <- sapply(
	sc_data[, -c('id', 'patient', 'cell_type')],
	function(x) {log2(mean(10*(2^x - 1)) + 1)},
	USE.NAMES = TRUE
)

sc_data <- sc_data[, c('id', 'patient', names(gene_averages)[gene_averages >= 4]), with = FALSE]

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

# The genes should be the same in all of the inferred_cna matrices:
# gene_order <- chr.order(rownames(inferred_cna[[1]]))

pdf('../data_and_figures/cna_plots_hnscc_comparison_caf.pdf')

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
		# geom_vline(xintercept = gene_order$chr_breaks[-c(1, length(gene_order$chr_breaks))], size = 0.25) +
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
