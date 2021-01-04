library(data.table) # 1.12.8
library(Rtsne) # 0.15





cohort_data <- list(
	breast_chung = list(
		read_quote = quote(fread('../data_and_figures/chung_breast_cancer_2017.csv')[, -c('cell_type', 'lymph_node')]),
		seed = 2470
	),
	breast_karaayvaz = list(
		read_quote = quote(fread('../data_and_figures/karaayvaz_tnbc_2018.csv')[, -'cell_type']),
		seed = 2604
	),
	breast_qian = list(
		read_quote = quote(fread('../data_and_figures/qian_breast_2020.csv')[, -'cell_type']),
		seed = 6500
	),
	crc_li = list(
		read_quote = quote(fread('../data_and_figures/li_colorectal_2017.csv')[
			sample_type != 'normal',
			-c('cell_type', 'sample_type', 'epithelial_cluster')
		]),
		seed = 5788
	),
	crc_lee_smc = list(
		read_quote = quote(fread('../data_and_figures/lee_crc_2020_smc.csv')[
			sample_type != 'normal',
			-c('sample_id', 'sample_type', 'cell_type', 'cell_subtype')
		]),
		seed = 2853
	),
	# Note the KUL3 dataset also has "border" samples, which I'm retaining because I assume they're still part of the tumour.
	crc_lee_kul3 = list(
		read_quote = quote(fread('../data_and_figures/lee_crc_2020_kul3.csv')[
			sample_type != 'normal',
			-c('sample_id', 'sample_type', 'cell_type', 'cell_subtype')
		]),
		seed = 4172
	),
	crc_qian = list(
		read_quote = quote(fread('../data_and_figures/qian_crc_2020.csv')[sample_type != 'normal', -c('sample_type', 'cell_type')]),
		seed = 9154
	),
	hnscc_puram = list(
		read_quote = quote(fread('../data_and_figures/puram_hnscc_2017.csv')[, -c('cell_type', 'lymph_node', 'processed_by_maxima_enzyme')]),
		seed = 5083
	),
	liver_ma = list(
		read_quote = quote(fread('../data_and_figures/ma_liver_2019.csv')[, -c('cell_type', 'sample', 'disease')]),
		seed = 2553
	),
	liver_ma_hcc = list(
		read_quote = quote(fread('../data_and_figures/ma_liver_2019.csv')[disease == 'HCC', -c('cell_type', 'sample', 'disease')]),
		seed = 8972
	),
	liver_ma_icca = list(
		read_quote = quote(fread('../data_and_figures/ma_liver_2019.csv')[disease == 'iCCA', -c('cell_type', 'sample', 'disease')]),
		seed = 8988
	),
	luad_kim = list(
		read_quote = quote(
			fread('../data_and_figures/kim_luad_2020.csv')[
				,
				-c('cell_type', 'cell_type_refined', 'cell_subtype', 'sample_site', 'sample_id', 'sample_origin')
			]
		),
		seed = 3118
	),
	lung_lambrechts = list(
		read_quote = quote(fread('../data_and_figures/lambrechts_nsclc_2018.csv')[
			sample_type != 'normal',
			-c('cell_type', 'cluster', 'disease', 'sample_type', 'annotation')
		]),
		seed = 4031
	),
	lung_lambrechts_luad = list(
		read_quote = quote(fread('../data_and_figures/lambrechts_nsclc_2018.csv')[
			sample_type != 'normal' & disease == 'LUAD',
			-c('cell_type', 'cluster', 'disease', 'sample_type', 'annotation')
		]),
		seed = 576
	),
	lung_lambrechts_lusc = list(
		read_quote = quote(fread('../data_and_figures/lambrechts_nsclc_2018.csv')[
			sample_type != 'normal' & disease == 'LUSC',
			-c('cell_type', 'cluster', 'disease', 'sample_type', 'annotation')
		]),
		seed = 4394
	),
	lung_qian = list(
		read_quote = quote(fread('../data_and_figures/qian_lung_2020.csv')[sample_type != 'normal', -c('sample_type', 'disease', 'cell_type')]),
		seed = 9656
	),
	lung_qian_luad = list(
		read_quote = quote(
			fread('../data_and_figures/qian_lung_2020.csv')[sample_type != 'normal' & disease == 'LUAD', -c('sample_type', 'disease', 'cell_type')]
		),
		seed = 3615
	),
	lung_qian_lusc = list(
		read_quote = quote(
			fread('../data_and_figures/qian_lung_2020.csv')[sample_type != 'normal' & disease == 'LUSC', -c('sample_type', 'disease', 'cell_type')]
		),
		seed = 1684
	),
	lung_song = list(
		read_quote = quote(fread('../data_and_figures/song_nsclc_2019.csv')[sample_type == 'tumour', -c('sample_type', 'cell_type', 'disease')]),
		seed = 6114
	),
	lung_song_luad = list(
		read_quote = quote(
			fread('../data_and_figures/song_nsclc_2019.csv')[sample_type == 'tumour' & disease == 'LUAD', -c('sample_type', 'cell_type', 'disease')]
		),
		seed = 8619
	),
	lung_song_lusc = list(
		read_quote = quote(
			fread('../data_and_figures/song_nsclc_2019.csv')[sample_type == 'tumour' & disease == 'LUSC', -c('sample_type', 'cell_type', 'disease')]
		),
		seed = 7720
	),
	ovarian_izar_ss2 = list(
		read_quote = quote(fread('../data_and_figures/izar_ovarian_2020_ss2.csv')[, -c('cell_type', 'cluster')]),
		seed = 4294
	),
	ovarian_izar_10x = list(
		read_quote = quote(fread('../data_and_figures/izar_ovarian_2020_10x.csv')[, -c('barcode', 'sample_id', 'time', 'cell_type', 'cluster')]),
		seed = 391
	),
	ovarian_qian = list(
		read_quote = quote(fread('../data_and_figures/qian_ovarian_2020.csv')[sample_type != 'normal', -c('sample_type', 'site', 'cell_type')]),
		seed = 2624
	),
	pdac_elyada = list(
		read_quote = quote(fread('../data_and_figures/elyada_pdac_2019.csv')[, -c('cell_type', 'cluster')]),
		seed = 7616
	),
	pdac_peng = list(
		read_quote = quote(fread('../data_and_figures/peng_pdac_2019.csv')[sample_type != 'normal', -c('cell_type', 'sample_type')]),
		seed = 6776
	)
)





for(cohort in names(cohort_data)) {

	cat(paste0(cohort, ':\n'))

	cat('\tReading in data\n')

	sc_data <- eval(cohort_data[[cohort]]$read_quote)

	# Take genes with log average TPM at least 4:

	cat('\tFiltering genes\n')

	gene_averages <- sapply(sc_data[, -c('id', 'patient')], function(x) {log2(mean(10*(2^x - 1)) + 1)}, USE.NAMES = TRUE)

	sc_data <- sc_data[, c('id', 'patient', names(gene_averages)[gene_averages >= 4]), with = FALSE]

	# Clustering by t-SNE and DBSCAN:

	cat('\tRunning t-SNE\n')

	set.seed(cohort_data[[cohort]]$seed)
	sc_tsne <- Rtsne(as.matrix(sc_data[, lapply(.SD, function(x) {x - mean(x)}), .SDcols = -c('id', 'patient')]))

	saveRDS(sc_tsne, paste0('../data_and_figures/tsne_', cohort, '.rds'))

	cat('\tDone!\n')

}
