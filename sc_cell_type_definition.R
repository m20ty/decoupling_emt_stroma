# Package version numbers listed below for R version 3.6.3 on WEXAC and on my laptop's WSL setup - same version on both unless stated otherwise.

library(data.table) # 1.12.8
library(ggplot2) # This is version 3.3.1 on my laptop's WSL setup but 3.3.0 on WEXAC.
library(magrittr) # 1.5
library(Rtsne) # 0.15
library(fpc) # 2.2.5
library(dbscan) # 1.1.5
library(infercna) # 1.0.0, but only on my laptop's WSL setup (it's the reason I started using WSL)

source('general_functions.R')

cell_type_markers <- fread('../../cell_type_markers.csv')





cohort_data <- list(
	breast_chung = list(
		read_quote = quote(fread('../data_and_figures/chung_breast_cancer_2017.csv')[, -c('cell_type', 'lymph_node')]),
		seed = 4292,
		min_cells = 10,
		epsilon = 1.8
	),
	breast_karaayvaz = list(
		read_quote = quote(fread('../data_and_figures/karaayvaz_tnbc_2018.csv')[, -'cell_type']),
		seed = 7733,
		min_cells = 20,
		epsilon = 2.8
	),
	breast_qian = list(
		read_quote = quote(fread('../data_and_figures/qian_breast_2020.csv')[, -'cell_type']),
		seed = 1776,
		min_cells = 100,
		epsilon = 2.7
	),
	# The CRC data from Li et al. doesn't cluster at all by t-SNE (probably due to low cell number), so I'll leave it out.
	# crc_li = list(
		# read_quote = quote(fread('../data_and_figures/li_colorectal_2017.csv')[
			# sample_type != 'normal',
			# -c('cell_type', 'sample_type', 'epithelial_cluster')
		# ]),
		# seed = 2456,
		# min_cells = 20,
		# epsilon = 2
	# ),
	# Note the KUL3 dataset also has "border" samples, which I'm retaining because I assume they're still part of the tumour.
	crc_lee_kul3 = list(
		read_quote = quote(fread('../data_and_figures/lee_crc_2020_kul3.csv')[
			sample_type != 'normal',
			-c('sample_id', 'sample_type', 'cell_type', 'cell_subtype')
		]),
		seed = 499,
		min_cells = 50,
		epsilon = 3
	),
	crc_lee_smc = list(
		read_quote = quote(fread('../data_and_figures/lee_crc_2020_smc.csv')[
			sample_type != 'normal',
			-c('sample_id', 'sample_type', 'cell_type', 'cell_subtype')
		]),
		seed = 4480,
		min_cells = 50,
		epsilon = 1.9
	),
	crc_qian = list(
		read_quote = quote(fread('../data_and_figures/qian_crc_2020.csv')[sample_type != 'normal', -c('sample_type', 'cell_type')]),
		seed = 7745,
		min_cells = 50,
		epsilon = 2.7
	),
	hnscc_puram = list(
		read_quote = quote(fread('../data_and_figures/puram_hnscc_2017.csv')[, -c('cell_type', 'lymph_node', 'processed_by_maxima_enzyme')]),
		seed = 7969,
		min_cells = 20,
		epsilon = 2.3
	),
	liver_ma = list(
		read_quote = quote(fread('../data_and_figures/ma_liver_2019.csv')[, -c('cell_type', 'sample', 'disease')]),
		seed = 6049,
		min_cells = 20,
		epsilon = 2.3
	),
	liver_ma_hcc = list(
		read_quote = quote(fread('../data_and_figures/ma_liver_2019.csv')[disease == 'HCC', -c('cell_type', 'sample', 'disease')]),
		seed = 5242,
		min_cells = 20,
		epsilon = 2.8
	),
	liver_ma_icca = list(
		read_quote = quote(fread('../data_and_figures/ma_liver_2019.csv')[disease == 'iCCA', -c('cell_type', 'sample', 'disease')]),
		seed = 420,
		min_cells = 10,
		epsilon = 1.9
	),
	luad_kim = list(
		read_quote = quote(fread('../data_and_figures/kim_luad_2020.csv')[, -c('cell_type', 'cell_type_refined', 'cell_subtype', 'sample_site', 'sample_id', 'sample_origin')]),
		seed = 6865,
		min_cells = 50,
		epsilon = 1.8
	),
	lung_lambrechts = list(
		read_quote = quote(fread('../data_and_figures/lambrechts_nsclc_2018.csv')[
			sample_type != 'normal',
			-c('cell_type', 'cluster', 'disease', 'sample_type', 'annotation')
		]),
		seed = 3480,
		min_cells = 20,
		epsilon = 1.9
	),
	lung_lambrechts_luad = list(
		read_quote = quote(fread('../data_and_figures/lambrechts_nsclc_2018.csv')[
			sample_type != 'normal' & disease == 'LUAD',
			-c('cell_type', 'cluster', 'disease', 'sample_type', 'annotation')
		]),
		seed = 6647,
		min_cells = 20,
		epsilon = 2
	),
	lung_lambrechts_lusc = list(
		read_quote = quote(fread('../data_and_figures/lambrechts_nsclc_2018.csv')[
			sample_type != 'normal' & disease == 'LUSC',
			-c('cell_type', 'cluster', 'disease', 'sample_type', 'annotation')
		]),
		seed = 5819,
		min_cells = 20,
		epsilon = 3.2
	),
	lung_qian = list(
		read_quote = quote(fread('../data_and_figures/qian_lung_2020.csv')[sample_type != 'normal', -c('sample_type', 'disease', 'cell_type')]),
		seed = 3308,
		min_cells = 50,
		epsilon = 2
	),
	lung_qian_luad = list(
		read_quote = quote(fread('../data_and_figures/qian_lung_2020.csv')[sample_type != 'normal' & disease == 'LUAD', -c('sample_type', 'disease', 'cell_type')]),
		seed = 9855,
		min_cells = 20,
		epsilon = 2
	),
	lung_qian_lusc = list(
		read_quote = quote(fread('../data_and_figures/qian_lung_2020.csv')[sample_type != 'normal' & disease == 'LUSC', -c('sample_type', 'disease', 'cell_type')]),
		seed = 1931,
		min_cells = 20,
		epsilon = 2.2
	),
	lung_song = list(
		read_quote = quote(fread('../data_and_figures/song_nsclc_2019.csv')[sample_type == 'tumour', -c('cell_type', 'sample_type', 'disease')]),
		seed = 7050,
		min_cells = 20,
		epsilon = 2.5
	),
	lung_song_luad = list(
		read_quote = quote(fread('../data_and_figures/song_nsclc_2019.csv')[sample_type == 'tumour' & disease == 'LUAD', -c('cell_type', 'sample_type', 'disease')]),
		seed = 9062,
		min_cells = 20,
		epsilon = 2.7
	),
	# In the following, I get no data for average expression of B cell markers, presumably because they don't pass the threshold of log average TPM >= 4.
	lung_song_lusc = list(
		read_quote = quote(fread('../data_and_figures/song_nsclc_2019.csv')[sample_type == 'tumour' & disease == 'LUSC', -c('cell_type', 'sample_type', 'disease')]),
		seed = 2735,
		min_cells = 20,
		epsilon = 2.9
	),
	ovarian_izar_10x = list(
		read_quote = quote(fread('../data_and_figures/izar_ovarian_2020_10x.csv')[, -c('barcode', 'sample_id', 'time', 'cell_type', 'cluster')]),
		seed = 3957,
		min_cells = 20,
		epsilon = 2.1
	),
	ovarian_izar_ss2 = list(
		read_quote = quote(fread('../data_and_figures/izar_ovarian_2020_ss2.csv')[, -c('cell_type', 'cluster')]),
		seed = 9096,
		min_cells = 20,
		epsilon = 2.6
	),
	ovarian_qian = list(
		read_quote = quote(fread('../data_and_figures/qian_ovarian_2020.csv')[sample_type != 'normal', -c('sample_type', 'site', 'cell_type')]),
		seed = 3758,
		min_cells = 50,
		epsilon = 2.6
	),
	pdac_elyada = list(
		read_quote = quote(fread('../data_and_figures/elyada_pdac_2019.csv')[, -c('cell_type', 'cluster')]),
		seed = 8305,
		min_cells = 50,
		epsilon = 2.7
	),
	pdac_peng = list(
		read_quote = quote(fread('../data_and_figures/peng_pdac_2019.csv')[sample_type != 'normal', -c('cell_type', 'sample_type')]),
		seed = 2765,
		min_cells = 50,
		epsilon = 1.6
	)
)





cohort <- 'breast_chung'

sc_data <- eval(cohort_data[[cohort]]$read_quote)

# Take genes with log average TPM at least 4:

gene_averages <- sapply(
	sc_data[, -c('id', 'patient')],
	function(x) {log2(mean(10*(2^x - 1)) + 1)},
	USE.NAMES = TRUE
)

sc_data <- sc_data[, c('id', 'patient', names(gene_averages)[gene_averages >= 4]), with = FALSE]

# Read in t-SNE:

sc_tsne <- readRDS(paste0('../data_and_figures/tsne_', cohort, '.rds'))

# Check t-SNE plot for obvious problems:

qplot(x, y, data = setNames(as.data.table(sc_tsne$Y), c('x', 'y')))

# Choose minimum number of cells per cluster and check the k-NN distance plot:

dbscan::kNNdistplot(sc_tsne$Y, k = cohort_data[[cohort]]$min_cells - 1)
abline(h = cohort_data[[cohort]]$epsilon)

# Run DBSCAN:

set.seed(cohort_data[[cohort]]$seed)
sc_dbscan <- fpc::dbscan(sc_tsne$Y, eps = cohort_data[[cohort]]$epsilon, MinPts = cohort_data[[cohort]]$min_cells)

saveRDS(sc_dbscan, paste0('../data_and_figures/dbscan_', cohort, '.rds'))

# sc_dbscan <- readRDS(paste0('../data_and_figures/dbscan_', cohort, '.rds'))

# Scatter plot of t-SNE coordinates, coloured by DBSCAN cluster:

tsne_plot_dbscan <- ggplot(
    setNames(
        cbind(as.data.table(sc_tsne$Y), as.character(sc_dbscan$cluster)),
        c('x', 'y', 'dbscan_cluster')
    )[dbscan_cluster != 0],
    aes(x = x, y = y, colour = dbscan_cluster)
) +
    geom_point() +
    scale_colour_manual(values = randomcoloR::distinctColorPalette(max(sc_dbscan$cluster))) +
    theme_minimal()

# Coloured by patient:

tsne_plot_patient <- ggplot(
    setNames(
        cbind(as.data.table(sc_tsne$Y), as.character(sc_dbscan$cluster), as.character(sc_data$patient)),
        c('x', 'y', 'dbscan_cluster', 'patient')
    )[dbscan_cluster != 0],
    aes(x = x, y = y, colour = patient)
) +
    geom_point() +
	scale_colour_manual(values = randomcoloR::distinctColorPalette(length(unique(sc_data$patient)))) +
    theme_minimal()

# Colour by average expression of cell type markers:

ct_ave_exp_plots <- sapply(
	cell_type_markers[cell_type !='mesenchymal', unique(cell_type)],
	function(ct) {
		
		ave_exp <- sc_data[, rowMeans(.SD), .SDcols = cell_type_markers[cell_type == ct & gene %in% names(sc_data), unique(gene)]]
		
		tsne_plot <- ggplot(
			setNames(
				cbind(as.data.table(sc_tsne$Y), as.character(sc_dbscan$cluster), ave_exp),
				c('x', 'y', 'dbscan_cluster', 'ave_exp')
			)[dbscan_cluster != 0],
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

pdf(paste0('../data_and_figures/tsne_plots_', cohort, '.pdf'))

c(
	list(tsne_plot_patient, tsne_plot_dbscan),
	lapply(ct_ave_exp_plots, `[[`, 'plot')
)

dev.off()





# HNSCC:

sc_data <- fread('../data_and_figures/puram_hnscc_2017.csv')[, -c('lymph_node', 'processed_by_maxima_enzyme')]

# Start by taking genes with log average TPM at least 4:

gene_averages <- sapply(
	sc_data[, -c('id', 'patient', 'cell_type')],
	function(x) {log2(mean(10*(2^x - 1)) + 1)},
	USE.NAMES = TRUE
)

sc_data <- sc_data[, c('id', 'patient', 'cell_type', names(gene_averages)[gene_averages >= 4]), with = FALSE]

# Clustering by t-SNE and DBSCAN:

set.seed(2450)
sc_tsne <- Rtsne(as.matrix(sc_data[, lapply(.SD, function(x) {x - mean(x)}), .SDcols = -c('id', 'patient', 'cell_type')]))

saveRDS(sc_tsne, '../data_and_figures/tsne_hnscc.rds')

# qplot(x, y, data = setNames(as.data.table(sc_tsne$Y), c('x', 'y')))

# dbscan::kNNdistplot(sc_tsne$Y, k = 19)
# abline(h = 2.3)

set.seed(484)
sc_dbscan <- fpc::dbscan(sc_tsne$Y, eps = 2.3, MinPts = 20)

# Scatter plot of t-SNE coordinates, coloured by DBSCAN cluster:

tsne_plot_dbscan <- ggplot(
    setNames(
        cbind(as.data.table(sc_tsne$Y), as.character(sc_dbscan$cluster)),
        c('x', 'y', 'dbscan_cluster')
    )[dbscan_cluster != 0],
    aes(x = x, y = y, colour = dbscan_cluster)
) +
    geom_point() +
    scale_colour_manual(values = randomcoloR::distinctColorPalette(max(sc_dbscan$cluster))) +
    theme_minimal()

# The following is to create a scatterplot with plotly, so we can easily identify which cluster is which, and
# save it to html using htmlwidgets.  However, without pandoc, we can't save it as a self-contained widget: it
# has to save all the assocaiated files in a separate folder, which is a bit messy.  If we just install pandoc
# with sudo apt-get install pandoc, will this be enough?

# tsne_plotly_dbscan <- plot_ly(
	# data = setNames(
        # cbind(as.data.table(sc_tsne$Y), as.character(sc_dbscan$cluster)),
        # c('x', 'y', 'dbscan_cluster')
    # )[dbscan_cluster != 0],
    # x = ~x,
    # y = ~y,
    # color = ~dbscan_cluster,
	# colors = randomcoloR::distinctColorPalette(max(sc_dbscan$cluster))
# )

# htmlwidgets::saveWidget(tsne_plotly_dbscan, '../data_and_figures/tsne_hnscc_dbscan.html')

# Coloured by patient:

tsne_plot_patient <- ggplot(
    setNames(
        cbind(as.data.table(sc_tsne$Y), as.character(sc_dbscan$cluster), as.character(sc_data$patient)),
        c('x', 'y', 'dbscan_cluster', 'patient')
    )[dbscan_cluster != 0],
    aes(x = x, y = y, colour = patient)
) +
    geom_point() +
	scale_colour_manual(values = randomcoloR::distinctColorPalette(length(unique(sc_data$patient)))) +
    theme_minimal()

# Colour by average expression of cell type markers:

ct_ave_exp_plots <- sapply(
	cell_type_markers[cell_type !='mesenchymal', unique(cell_type)],
	function(ct) {
		
		ave_exp <- sc_data[, rowMeans(.SD), .SDcols = cell_type_markers[cell_type == ct & gene %in% names(sc_data), unique(gene)]]
		
		tsne_plot <- ggplot(
			setNames(
				cbind(as.data.table(sc_tsne$Y), as.character(sc_dbscan$cluster), ave_exp),
				c('x', 'y', 'dbscan_cluster', 'ave_exp')
			)[dbscan_cluster != 0],
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

pdf('../data_and_figures/tsne_hnscc_plots.pdf')

c(
	list(tsne_plot_dbscan),
	lapply(ct_ave_exp_plots, `[[`, 'plot')
)

dev.off()

# We see that cluster 4 is B or B plasma cells; 8 is macrophages and DC cells; 10 and 14 both seem to be endothelial; 17 is mast cells; and
# 11, 18 and 19 are T cells.  Myocytes seem very similar to CAFs, so I won't include them.

# Try inferring CNAs using the above clusters as reference (use all reference cells for each tumour):

ref_cells <- list(
	b_cell = sc_data$id[sc_dbscan$cluster == 4],
	macrophage_dc = sc_data$id[sc_dbscan$cluster == 8],
	endothelial = sc_data$id[sc_dbscan$cluster %in% c(10, 14)],
	mast = sc_data$id[sc_dbscan$cluster == 17],
	t_cell = sc_data$id[sc_dbscan$cluster %in% c(11, 18, 19)]
)

# ref_cells <- sapply(
	# as.character(unique(sc_data$patient)),
	# function(p) {
		# list(
			# b_cell = sc_data[patient == p & sc_dbscan$cluster == 4, id],
			# macrophage_dc = sc_data[patient == p & sc_dbscan$cluster == 8, id],
			# endothelial = sc_data[patient == p & sc_dbscan$cluster %in% c(10, 14), id],
			# mast = sc_data[patient == p & sc_dbscan$cluster == 17, id],
			# t_cell = sc_data[patient == p & sc_dbscan$cluster %in% c(11, 18, 19), id]
		# )
	# },
	# simplify = FALSE,
	# USE.NAMES = TRUE
# )

# inferred_cna <- infercna(
	# t(set_rownames(as.matrix(sc_data[, -c('id', 'patient', 'cell_type')]), sc_data$id)),
	# refCells = ref_cells,
	# # n = 5000,
	# # noise = 0.1,
	# isLog = TRUE
# )

inferred_cna <- sapply(
	as.character(unique(sc_data$patient)),
	function(p) {
		infercna(
			t(
				set_rownames(
					as.matrix(sc_data[patient == p | id %in% unlist(ref_cells), -c('id', 'patient', 'cell_type')]),
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

saveRDS(inferred_cna, '../data_and_figures/inferred_cna_hnscc.rds')

# This doesn't work for some reason:
# cna_plot <- cnaPlot(inferred_cna[, !(colnames(inferred_cna) %in% unlist(ref_cells))], order.cells = TRUE)

cnaScatterPlot(inferred_cna, gene.quantile = 0.9)

# malignant_cells <- findMalignant(inferred_cna, refCells = ref_cells)

malignant_cells <- sapply(
	inferred_cna,
	findMalignant,
	refCells = ref_cells,
	simplify = FALSE,
	USE.NAMES = TRUE
)

# None of this inferring CNAs stuff worked well, and I think it's because I need to do everything per tumour.  It's intuitive, because we know
# tumours will have very different CNAs, so when we take the average CNA profile (correlation with which defines the malignant cells) we will
# average out a lot of signal.  Fixing this will hopefully fix the error with cnaPlot().
