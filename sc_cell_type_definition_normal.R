library(data.table) # 1.12.8
library(ggplot2) # This is version 3.3.1 on my laptop's WSL setup but 3.3.0 on WEXAC.
library(magrittr) # 1.5
library(Rtsne) # 0.15
library(fpc) # 2.2.5
library(dbscan) # 1.1.5
library(infercna) # 1.0.0

source('general_functions.R')

cell_type_markers <- fread('../../cell_type_markers.csv')





cohort_data <- list(
	crc_lee_kul3 = list(
		read_quote = quote(fread('../data_and_figures/lee_crc_2020_kul3.csv')[
			sample_type == 'normal',
			-c('sample_id', 'sample_type', 'cell_type', 'cell_subtype')
		]),
		seed = 7034,
		min_cells = 50,
		epsilon = 4
	),
	crc_lee_smc = list(
		read_quote = quote(fread('../data_and_figures/lee_crc_2020_smc.csv')[
			sample_type == 'normal',
			-c('sample_id', 'sample_type', 'cell_type', 'cell_subtype')
		]),
		seed = 2906,
		min_cells = 50,
		epsilon = 3.2
	),
	crc_qian = list(
		read_quote = quote(fread('../data_and_figures/qian_crc_2020.csv')[sample_type == 'normal', -c('sample_type', 'cell_type')]),
		seed = 1635,
		min_cells = 50,
		epsilon = 4
	),
	lung_lambrechts = list(
		read_quote = quote(fread('../data_and_figures/lambrechts_nsclc_2018.csv')[
			sample_type == 'normal',
			-c('cell_type', 'cluster', 'disease', 'sample_type', 'annotation')
		]),
		seed = 7045,
		min_cells = 20,
		epsilon = 2.1
	),
	lung_lambrechts_luad = list(
		read_quote = quote(fread('../data_and_figures/lambrechts_nsclc_2018.csv')[
			sample_type == 'normal' & disease == 'LUAD',
			-c('cell_type', 'cluster', 'disease', 'sample_type', 'annotation')
		]),
		seed = 3312,
		min_cells = 20,
		epsilon = 2.4
	),
	# There are only 122 normal LUSC cells in the Lambrechts dataset, so we should ignore this one.
	# lung_lambrechts_lusc = list(
		# read_quote = quote(fread('../data_and_figures/lambrechts_nsclc_2018.csv')[
			# sample_type == 'normal' & disease == 'LUSC',
			# -c('cell_type', 'cluster', 'disease', 'sample_type', 'annotation')
		# ]),
		# seed = 5907,
		# min_cells = 20,
		# epsilon = 3.2
	# ),
	lung_qian = list(
		read_quote = quote(fread('../data_and_figures/qian_lung_2020.csv')[sample_type == 'normal', -c('sample_type', 'disease', 'cell_type')]),
		seed = 7248,
		min_cells = 50,
		epsilon = 2.9
	),
	lung_qian_luad = list(
		read_quote = quote(fread('../data_and_figures/qian_lung_2020.csv')[sample_type == 'normal' & disease == 'LUAD', -c('sample_type', 'disease', 'cell_type')]),
		seed = 7817,
		min_cells = 20,
		epsilon = 2.3
	),
	lung_qian_lusc = list(
		read_quote = quote(fread('../data_and_figures/qian_lung_2020.csv')[sample_type == 'normal' & disease == 'LUSC', -c('sample_type', 'disease', 'cell_type')]),
		seed = 445,
		min_cells = 20,
		epsilon = 2.2
	),
	lung_song = list(
		read_quote = quote(fread('../data_and_figures/song_nsclc_2019.csv')[sample_type == 'normal', -c('cell_type', 'sample_type', 'disease')]),
		seed = 1660,
		min_cells = 20,
		epsilon = 2.5
	),
	lung_song_luad = list(
		read_quote = quote(fread('../data_and_figures/song_nsclc_2019.csv')[sample_type == 'normal' & disease == 'LUAD', -c('cell_type', 'sample_type', 'disease')]),
		seed = 2184,
		min_cells = 20,
		epsilon = 2.8
	),
	lung_song_lusc = list(
		read_quote = quote(fread('../data_and_figures/song_nsclc_2019.csv')[sample_type == 'normal' & disease == 'LUSC', -c('cell_type', 'sample_type', 'disease')]),
		seed = 7212,
		min_cells = 20,
		epsilon = 2.9
	),
	ovarian_qian = list(
		read_quote = quote(fread('../data_and_figures/qian_ovarian_2020.csv')[sample_type == 'normal', -c('sample_type', 'site', 'cell_type')]),
		seed = 7962,
		min_cells = 50,
		epsilon = 3.4
	),
	pdac_peng = list(
		read_quote = quote(fread('../data_and_figures/peng_pdac_2019.csv')[sample_type == 'normal', -c('cell_type', 'sample_type')]),
		seed = 4327,
		min_cells = 50,
		epsilon = 2.6
	)
)





cohort <- 'crc_lee_kul3'

sc_data <- eval(cohort_data[[cohort]]$read_quote)

# Take genes with log average TPM at least 4:

gene_averages <- sapply(
	sc_data[, -c('id', 'patient')],
	function(x) {log2(mean(10*(2^x - 1)) + 1)},
	USE.NAMES = TRUE
)

sc_data <- sc_data[, c('id', 'patient', names(gene_averages)[gene_averages >= 4]), with = FALSE]

# Read in t-SNE:

sc_tsne <- readRDS(paste0('../data_and_figures/tsne_normal_', cohort, '.rds'))

# Check t-SNE plot for obvious problems:

qplot(x, y, data = setNames(as.data.table(sc_tsne$Y), c('x', 'y')))

# Choose minimum number of cells per cluster and check the k-NN distance plot:

dbscan::kNNdistplot(sc_tsne$Y, k = cohort_data[[cohort]]$min_cells - 1)
abline(h = cohort_data[[cohort]]$epsilon)

# Run DBSCAN:

set.seed(cohort_data[[cohort]]$seed)
sc_dbscan <- fpc::dbscan(sc_tsne$Y, eps = cohort_data[[cohort]]$epsilon, MinPts = cohort_data[[cohort]]$min_cells)

saveRDS(sc_dbscan, paste0('../data_and_figures/dbscan_normal_', cohort, '.rds'))

# sc_dbscan <- readRDS(paste0('../data_and_figures/dbscan_normal_', cohort, '.rds'))





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

epithelial_markers <- c(
	'CDH1',
	'EPCAM',
	'SFN',
	names(sc_data)[grepl('^KRT[0-9]', names(sc_data))]
)

epithelial_markers <- epithelial_markers[epithelial_markers %in% names(sc_data)]

ave_exp_epithelial <- sc_data[, rowMeans(.SD), .SDcols = epithelial_markers]

tsne_plot_epithelial <- ggplot(
	setNames(
		cbind(as.data.table(sc_tsne$Y), as.character(sc_dbscan$cluster), ave_exp_epithelial),
		c('x', 'y', 'dbscan_cluster', 'ave_exp')
	)[dbscan_cluster != 0],
	aes(x = x, y = y, colour = ave_exp)
) +
	geom_point() +
	theme_minimal() +
	labs(title = 'epithelial')

pdf(paste0('../data_and_figures/tsne_plots_normal_', cohort, '.pdf'))

c(
	list(tsne_plot_patient, tsne_plot_dbscan, tsne_plot_epithelial),
	lapply(ct_ave_exp_plots, `[[`, 'plot')
)

dev.off()
