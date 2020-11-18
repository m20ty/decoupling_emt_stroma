# bsub -q tirosh -n 8 -R "rusage[mem=2000]" -o sc_tsne_normal_log.o -e sc_tsne_normal_log.e Rscript sc_tsne_normal.R

library(data.table) # 1.12.8
library(Rtsne) # 0.15

# The following is just to check the R version:
cat(R.Version()$version.string, '\n')

cohort_data <- list(
	crc_lee_kul3 = list(
		read_quote = quote(fread('../data_and_figures/lee_crc_2020_kul3.csv')[
			sample_type == 'normal',
			-c('sample_id', 'sample_type', 'cell_type', 'cell_subtype')
		]),
		seed = 4876
	),
	crc_lee_smc = list(
		read_quote = quote(fread('../data_and_figures/lee_crc_2020_smc.csv')[
			sample_type == 'normal',
			-c('sample_id', 'sample_type', 'cell_type', 'cell_subtype')
		]),
		seed = 7428
	),
	crc_qian = list(
		read_quote = quote(fread('../data_and_figures/qian_crc_2020.csv')[sample_type == 'normal', -c('sample_type', 'cell_type')]),
		seed = 3384
	),
	lung_lambrechts = list(
		read_quote = quote(fread('../data_and_figures/lambrechts_nsclc_2018.csv')[
			sample_type == 'normal',
			-c('cell_type', 'cluster', 'disease', 'sample_type', 'annotation')
		]),
		seed = 1567
	),
	lung_lambrechts_luad = list(
		read_quote = quote(fread('../data_and_figures/lambrechts_nsclc_2018.csv')[
			sample_type == 'normal' & disease == 'LUAD',
			-c('cell_type', 'cluster', 'disease', 'sample_type', 'annotation')
		]),
		seed = 8942
	),
	lung_lambrechts_lusc = list(
		read_quote = quote(fread('../data_and_figures/lambrechts_nsclc_2018.csv')[
			sample_type == 'normal' & disease == 'LUSC',
			-c('cell_type', 'cluster', 'disease', 'sample_type', 'annotation')
		]),
		seed = 1313
	),
	lung_qian = list(
		read_quote = quote(fread('../data_and_figures/qian_lung_2020.csv')[sample_type == 'normal', -c('sample_type', 'disease', 'cell_type')]),
		seed = 2815
	),
	lung_qian_luad = list(
		read_quote = quote(fread('../data_and_figures/qian_lung_2020.csv')[sample_type == 'normal' & disease == 'LUAD', -c('sample_type', 'disease', 'cell_type')]),
		seed = 4148
	),
	lung_qian_lusc = list(
		read_quote = quote(fread('../data_and_figures/qian_lung_2020.csv')[sample_type == 'normal' & disease == 'LUSC', -c('sample_type', 'disease', 'cell_type')]),
		seed = 5864
	),
	lung_song = list(
		read_quote = quote(fread('../data_and_figures/song_nsclc_2019.csv')[sample_type == 'normal', -c('cell_type', 'sample_type', 'disease')]),
		seed = 1929
	),
	lung_song_luad = list(
		read_quote = quote(fread('../data_and_figures/song_nsclc_2019.csv')[sample_type == 'normal' & disease == 'LUAD', -c('cell_type', 'sample_type', 'disease')]),
		seed = 2479
	),
	lung_song_lusc = list(
		read_quote = quote(fread('../data_and_figures/song_nsclc_2019.csv')[sample_type == 'normal' & disease == 'LUSC', -c('cell_type', 'sample_type', 'disease')]),
		seed = 6029
	),
	ovarian_qian = list(
		read_quote = quote(fread('../data_and_figures/qian_ovarian_2020.csv')[sample_type == 'normal', -c('sample_type', 'site', 'cell_type')]),
		seed = 930
	),
	pdac_peng = list(
		read_quote = quote(fread('../data_and_figures/peng_pdac_2019.csv')[sample_type == 'normal', -c('cell_type', 'sample_type')]),
		seed = 5478
	)
)





for(cohort in names(cohort_data)) {
	
	cat(paste0(cohort, ':\n'))
	
	cat('\tReading in data\n')
	
	sc_data <- eval(cohort_data[[cohort]]$read_quote)
	
	# Start by taking genes with log average TPM at least 4:
	
	cat('\tFiltering genes\n')
	
	gene_averages <- sapply(
		sc_data[, -c('id', 'patient')],
		function(x) {log2(mean(10*(2^x - 1)) + 1)},
		USE.NAMES = TRUE
	)
	
	sc_data <- sc_data[, c('id', 'patient', names(gene_averages)[gene_averages >= 4]), with = FALSE]
	
	# Clustering by t-SNE and DBSCAN:
	
	cat('\tRunning t-SNE\n')
	
	set.seed(cohort_data[[cohort]]$seed)
	sc_tsne <- Rtsne(as.matrix(sc_data[, lapply(.SD, function(x) {x - mean(x)}), .SDcols = -c('id', 'patient')]))
	
	saveRDS(sc_tsne, paste0('../data_and_figures/tsne_normal_', cohort, '.rds'))
	
	cat('\tDone!\n')
	
}
