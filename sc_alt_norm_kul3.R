library(data.table) # 1.12.8
library(Matrix) # 1.2.18
library(stringr) # 1.4.0
library(plyr) # 1.8.6
library(limma) # 3.42.2
library(magrittr) # 1.5
library(SCnorm) # Currently on WEXAC but not my laptop
library(infercna) # 1.0.0
library(scran)
library(scater)
library(umap)
library(Rtsne)

source('general_functions.R')
source('sparse_matrix_functions.R')

cell_type_markers <- fread('../../cell_type_markers.csv')

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





sc_data <- readRDS('../../single_cell_data/lee_crc_2020/GSE144735_processed_KUL3_CRC_10X_raw_UMI_count_matrix.rds')
sc_meta <- fread('../../single_cell_data/lee_crc_2020/GSE144735_processed_KUL3_CRC_10X_annotation.txt', key = 'Index')

# Update gene names:
gn <- alias2SymbolTable(rownames(sc_data))
sc_data <- sc_data[!is.na(gn) & gn %in% names(table(gn))[table(gn) == 1], ]

# Remove cells with fewer than 1000 genes detected:
sc_data <- sc_data[, col_nnz(sc_data) >= 1000]

sc_meta <- sc_meta[colnames(sc_data)]

sc_meta <- sc_meta[
	,
	.(
		id = Index,
		patient = Patient,
		sample_id = Sample,
		sample_type = mapvalues(Class, c('Normal', 'Border', 'Tumor'), c('normal', 'border', 'tumour')),
		cell_type = mapvalues(
			Cell_type,
			c(
				'B cells',
				'Epithelial cells',
				'Mast cells',
				'Myeloids',
				'Stromal cells',
				'T cells'
			),
			c(
				'b_cell',
				'epithelial',
				'mast',
				'myeloid',
				'fibroblast',
				't_cell'
			)
		),
		cell_subtype = Cell_subtype
	)
]





# Try normalising using SCnorm, which may not be appropriate anyway without additional filtering (see first paragraph of
# 'SCnorm: UMI data' section in the vignette).  The vignette says it won't work well with very sparse data (in particular
# data with more than ~80% zero counts, which is true in this case), so perhaps we should filter genes, perhaps weakly,
# beforehand.  On the other hand, I think SCnorm already excludes genes with very low counts, so maybe there's no need to
# filter the genes...  And we already filtered the cells by number of genes detected.

# Could consider using patient numbers as conditions - then I think SCnorm will effectively normalise possible batch effects.

# For now, let's use just one patient.  We'll choose 'KUL21' because it has a decent amount of tumour cells, and among these,
# a decent number of cells classified as epithelial.

sc_data <- sc_data[, sc_meta[sample_type != 'normal' & patient == 'KUL21', id]]
sc_meta <- sc_meta[sample_type != 'normal' & patient == 'KUL21']

sc_data_scnorm <- SCnorm(as.matrix(sc_data), Conditions = rep(1, ncol(sc_data)), NCores = 4)

saveRDS(SingleCellExperiment::normcounts(sc_data_scnorm), '../data_and_figures/lee_crc_2020_kul3_KUL21_scnorm.rds')





# Try normalising using scran:

# Do we have to convert from sparse matrix?
# I use calculateSumFactors instead of computeSumFactors, because this is what I'm told to do via warning message, since I'm not
# using an object of class SingleCellExperiment.

scran_clusters <- quickCluster(as.matrix(sc_data))
scran_size_factors <- calculateSumFactors(as.matrix(sc_data), clusters = sc_data_scran_clusters)
sc_data_scran <- SingleCellExperiment(list(counts = sc_data))
sc_data_scran <- logNormCounts(sc_data_scran, size_factors = scran_size_factors)

saveRDS(sc_data_scran, '../data_and_figures/lee_crc_2020_kul3_KUL21_scran.rds')





# In the below, it doesn't make sense to use rowSums for filtering genes, because we'll have to adjust the threshold for
# different sizes of datasets.  It's better to use rowMeans, in which case maybe 0.15 or 0.2 would give about the right
# number of genes...

# I can think of two ways to mimic the usual log(average TPM/10 + 1) >= 4.  First, we could take the mean/median cell size
# (cell size meaning total transcript count in a cell), which in the SCnorm case is around 10000.  Then if
# log(average TPM/10 + 1) >= 4, average TPM/10 >= 15, so if the average size factor is a 10th of 100,000, then we want
# genes with averages above 1.5.  But this only leaves ~800-900 genes in this case.  We could also divide each cell's
# counts by the cell size, and look for genes with average above 15/100,000 = 0.00015.  This again leaves ~800-900 genes.
# This is obviously well short of the several thousand genes I retain using the usual TPM/10 threshold, and I'm not sure
# why this is.

# We could also look at the distribution of gene averages.  E.g. this density plot suggests that a threshold of 0.15 is
# about right if we just want to shave off the peak at zero (sc_data comes from the readRDS call below):

# plot(density(log2(rowMeans(sc_data) + 1)), xlim = c(0, 3))
# abline(v = 0.15, col = 'blue', lty = 'dashed')





# Try UMAP and t-SNE on SCnorm data:

sc_data <- readRDS('../data_and_figures/lee_crc_2020_kul3_KUL21_scnorm.rds')

# This density plot suggests that a threshold of 0.15 is about right if we just want to shave off the peak at zero:
# plot(density(log2(rowMeans(sc_data) + 1)), xlim = c(0, 3))
# abline(v = 0.15, col = 'blue', lty = 'dashed')

sc_data <- sc_data[log2(rowMeans(sc_data) + 1) >= 0.15, ]

# UMAP:
# set.seed(6389)
# UMAP works OK without the log transformation, but with it we see very distinct clusters.
# Would it make sense to centre the data?
# sc_umap <- umap(t(log2(sc_data[log2(rowMeans(sc_data) + 1) >= 0.15, ] + 1)))
# sc_umap <- umap(t(sc_data))
# sc_umap_data <- setNames(as.data.table(sc_umap$layout, keep.rownames = TRUE), c('cell_id', 'x', 'y'))
# ggplot(sc_umap_data, aes(x = x, y = y)) + geom_point() + theme_minimal()
# dbscan::kNNdistplot(sc_umap_data[, .(x, y)], k = 19)
# abline(h = 0.34)
# set.seed(5094)
# sc_dbscan <- fpc::dbscan(sc_umap_data[, .(x, y)], eps = 0.34, MinPts = 20)
# sc_umap_data <- cbind(sc_umap_data, dbscan_cluster = as.character(sc_dbscan$cluster))
# umap_plot_dbscan <- ggplot(sc_umap_data, aes(x = x, y = y, colour = dbscan_cluster)) +
	# geom_point() +
	# scale_colour_manual(values = randomcoloR::distinctColorPalette(length(unique(sc_umap_data$dbscan_cluster)))) +
	# theme_minimal()
# umap_ct_ave_exp_plots <- sapply(
	# cell_type_markers[cell_type !='mesenchymal', unique(cell_type)],
	# function(ct) {
		# ave_exp <- log2(colMeans(sc_data[cell_type_markers[cell_type == ct & gene %in% rownames(sc_data), unique(gene)], ]) + 1)
		# umap_plot <- ggplot(
			# cbind(sc_umap_data, ave_exp = ave_exp)[dbscan_cluster != 0],
			# aes(x = x, y = y, colour = ave_exp)
		# ) +
			# geom_point() +
			# theme_minimal() +
			# labs(title = ct)
		# list(ave_exp = ave_exp, plot = umap_plot)
	# },
	# simplify = FALSE,
	# USE.NAMES = TRUE
# )
# pdf(paste0('../data_and_figures/umap_plots_lee_crc_2020_kul3_KUL21_scnorm.pdf'))
# c(list(umap_plot_dbscan), lapply(umap_ct_ave_exp_plots, `[[`, 'plot'))
# dev.off()

# t-SNE:

set.seed(7665)

sc_tsne <- Rtsne(t(log2(sc_data + 1)))
sc_tsne_data <- cbind(cell_id = colnames(sc_data), setNames(as.data.table(sc_tsne$Y), c('x', 'y')))

# ggplot(sc_tsne_data, aes(x = x, y = y)) + geom_point() + theme_minimal()
dbscan::kNNdistplot(sc_tsne_data[, .(x, y)], k = 19)
abline(h = 2.7)

set.seed(9200)
sc_dbscan <- fpc::dbscan(sc_tsne_data[, .(x, y)], eps = 2.7, MinPts = 20)

sc_tsne_data <- cbind(sc_tsne_data, dbscan_cluster = as.character(sc_dbscan$cluster))

tsne_plot_dbscan <- ggplot(sc_tsne_data, aes(x = x, y = y, colour = dbscan_cluster)) +
	geom_point() +
	# scale_colour_manual(values = randomcoloR::distinctColorPalette(length(unique(sc_tsne_data$dbscan_cluster)))) +
	theme_minimal()

tsne_ct_ave_exp_plots <- sapply(
	cell_type_markers[cell_type !='mesenchymal', unique(cell_type)],
	function(ct) {
		
		ave_exp <- colMeans(log2(sc_data[cell_type_markers[cell_type == ct & gene %in% rownames(sc_data), unique(gene)], ] + 1))
		# ave_exp <- log2(colMeans(sc_data[cell_type_markers[cell_type == ct & gene %in% rownames(sc_data), unique(gene)], ]) + 1)
		
		tsne_plot <- ggplot(
			cbind(sc_tsne_data, ave_exp = ave_exp)[dbscan_cluster != 0],
			aes(x = x, y = y, colour = ave_exp)
		) +
			geom_point() +
			theme_minimal() +
			labs(title = ct)
		
		list(ave_exp = ave_exp, plot = tsne_plot)
		
	},
	simplify = FALSE,
	USE.NAMES = TRUE
)

pdf(paste0('../data_and_figures/tsne_plots_lee_crc_2020_kul3_KUL21_scnorm.pdf'))
c(list(tsne_plot_dbscan), lapply(tsne_ct_ave_exp_plots, `[[`, 'plot'))
dev.off()

ref_cell_clusters <- list(
	b_cell = 8,
	dc_t_cell = 7,
	endothelial = 2,
	macrophage = 6,
	mast = 9
)

ref_cells <- sapply(
	ref_cell_clusters,
	function(x) colnames(sc_data)[sc_dbscan$cluster %in% x],
	simplify = FALSE,
	USE.NAMES = TRUE
)

# This doesn't work if I just supply sc_data with isLog = FALSE.  Maybe raise this with Julie.
inferred_cna <- infercna(log2(sc_data + 1), refCells = ref_cells, isLog = TRUE)

x_line_breaks <- .chromosomeBreaks(rownames(inferred_cna))
x_text_breaks <- round(.chromosomeBreaks(rownames(inferred_cna), halfway = TRUE, hide = c('13', '18', '21', 'Y')))

cell_clust <- hclust(dist(t(inferred_cna)))

cna_heatmap <- ggplot(reshape2::melt(inferred_cna, varnames = c('gene', 'cell'), value.name = 'cna_score')) +
	geom_raster(
		aes(
			x = factor(gene, levels = rownames(inferred_cna)),
			y = factor(cell, levels = cell_clust$labels[cell_clust$order]),
			fill = cna_score
		)
	) +
	scale_fill_gradientn(
		colours = rev(colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))(50)),
		limits = c(-0.5, 0.5),
		oob = scales::squish
	) +
	scale_y_discrete(expand = c(0, 0)) +
	scale_x_discrete(expand = c(0, 0), breaks = rownames(inferred_cna)[x_text_breaks], labels = names(x_text_breaks)) +
	geom_vline(xintercept = x_line_breaks, size = 0.25) +
	theme(
		axis.text.y = element_blank(),
		axis.ticks = element_blank(),
		axis.ticks.length = unit(0, 'pt'),
		panel.border = element_rect(colour = 'black', size = 0.5, fill = NA)
	) +
	labs(x = 'Chromosome', y = 'Cells', fill = 'Inferred CNA', title = 'KUL21')

# This doesn't look as good as the TPM/10 version.  I had to reduce the range from (-1, 1) to (-0.5, 0.5), because the signal is much fainter.
# Also, there are clearly some cells where the signal is weaker than in others, though still there, which I would guess are the cells with
# lower overall counts.  TPM/10 would presumably help to make the signal more even between the cells.
