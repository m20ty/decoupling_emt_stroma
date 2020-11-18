# bsub -q new-short -R "rusage[mem=64000]" -o sc_nonmalignant_tsne_log.o -e sc_nonmalignant_tsne_log.e Rscript sc_nonmalignant_tsne.R

Sys.time()
cat(R.Version()$version.string, '\n')

library(data.table) # 1.12.8
library(ggplot2) # This is version 3.3.1 on my laptop's WSL setup but 3.3.0 on WEXAC.
library(magrittr) # 1.5
library(Rtsne) # 0.15

source('general_functions.R')





cohort_data <- list(
	breast_qian = list(
		patients = c(42, 43, 47, 49, 51, 53, 54),
		read_quote = quote(fread('../data_and_figures/qian_breast_2020.csv')[, -'cell_type']),
		seed = 2159
	),
	crc_lee_smc = list(
		patients = c('SMC01', 'SMC02', 'SMC04', 'SMC07', 'SMC08', 'SMC09', 'SMC11', 'SMC14', 'SMC15', 'SMC16', 'SMC18', 'SMC20', 'SMC21', 'SMC23', 'SMC25'),
		read_quote = quote(fread('../data_and_figures/lee_crc_2020_smc.csv')[
			sample_type != 'normal',
			-c('sample_id', 'sample_type', 'cell_type', 'cell_subtype')
		]),
		seed = 8587
	),
	hnscc_puram = list(
		patients = c(5, 6, 18, 20, 22, 25, 26),
		read_quote = quote(fread('../data_and_figures/puram_hnscc_2017.csv')[, -c('cell_type', 'lymph_node', 'processed_by_maxima_enzyme')]),
		seed = 216
	),
	liver_ma = list(
		patients = c('C25', 'C26', 'C46', 'C56', 'C66', 'H37', 'H38', 'H65'),
		read_quote = quote(fread('../data_and_figures/ma_liver_2019.csv')[, -c('cell_type', 'sample', 'disease')]),
		seed = 6565
	),
	luad_kim = list(
		patients = c('P0006', 'P0008', 'P0018', 'P0019', 'P0020', 'P0025', 'P0028', 'P0030', 'P0031', 'P0034', 'P1006', 'P1028', 'P1049', 'P1058'),
		read_quote = quote(fread('../data_and_figures/kim_luad_2020.csv')[, -c('cell_type', 'cell_type_refined', 'cell_subtype', 'sample_site', 'sample_id', 'sample_origin')]),
		seed = 6367
	),
	lung_qian = list(
		patients = 1:8,
		read_quote = quote(fread('../data_and_figures/qian_lung_2020.csv')[sample_type != 'normal', -c('sample_type', 'disease', 'cell_type')]),
		seed = 8176
	),
	ovarian_qian = list(
		patients = 11:14,
		read_quote = quote(fread('../data_and_figures/qian_ovarian_2020.csv')[sample_type != 'normal', -c('sample_type', 'site', 'cell_type')]),
		seed = 6403
	),
	pdac_peng = list(
		patients = c('T2', 'T3', 'T6', 'T7', 'T8', 'T9', 'T11', 'T13', 'T14', 'T15', 'T16', 'T17', 'T19', 'T20', 'T21', 'T22'),
		read_quote = quote(fread('../data_and_figures/peng_pdac_2019.csv')[sample_type != 'normal', -c('cell_type', 'sample_type')]),
		seed = 2371
	)
)





for(cohort in names(cohort_data)) {
	
	classification_data <- rbindlist(
		lapply(
			cohort_data[[cohort]]$patients,
			function(p) unique(
				fread(paste0('../data_and_figures/sc_find_malignant/', cohort, '/', p, '_data.csv'))[
					,
					.(cell_id, classification_final)
				]
			)
		)
	)
	
	classification_data[
		,
		classification := switch(
			(length(classification_final) == 1) + 1,
			switch( # In this case the cell must have been reference, but could have been classified as ambiguous or malignant in a subset of patients.
				('ambiguous' %in% classification_final | 'malignant' %in% classification_final) + 1,
				'nonmalignant',
				'ambiguous'
			),
			switch( # In this case the cell was not reference.
				grepl('^malignant', classification_final) + 1, # For some patients we have 'malignant clone 1' etc.
				classification_final,
				'malignant'
			)
		),
		by = cell_id
	]
	
	classification_data <- unique(classification_data[, .(cell_id, classification)])
	
	setkey(classification_data, cell_id)
	
	fwrite(classification_data, paste0('../data_and_figures/sc_reclassify/classification_data_', cohort, '.csv'))
	
	sc_data <- eval(cohort_data[[cohort]]$read_quote)[patient %in% cohort_data[[cohort]]$patients]
	
	sc_data[, classification := classification_data[id, classification]]
	
	sc_data <- sc_data[classification == 'nonmalignant', -'classification']
	
	gene_averages <- sapply(
		sc_data[, -c('id', 'patient')],
		function(x) {log2(mean(10*(2^x - 1)) + 1)},
		USE.NAMES = TRUE
	)
	
	sc_data <- sc_data[, c('id', 'patient', names(gene_averages)[gene_averages >= 4]), with = FALSE]
	
	# Run t-SNE:
	set.seed(cohort_data[[cohort]]$seed)
	sc_tsne <- Rtsne(as.matrix(sc_data[, lapply(.SD, function(x) {x - mean(x)}), .SDcols = -c('id', 'patient')]))
	
	saveRDS(sc_tsne, paste0('../data_and_figures/sc_reclassify/tsne_', cohort, '.rds'))
	
}
