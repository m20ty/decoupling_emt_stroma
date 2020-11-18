# bsub -q new-short -R "rusage[mem=64000]" -o sc_qc_checks_log.o -e sc_qc_checks_log.e Rscript sc_qc_checks.R

cat(R.Version()$version.string, '\n')

library(data.table)
library(ggplot2)

cohort_data <- list(
	breast_chung = quote(fread('../data_and_figures/chung_breast_cancer_2017.csv')[, -c('cell_type', 'lymph_node')]),
	breast_karaayvaz = quote(fread('../data_and_figures/karaayvaz_tnbc_2018.csv')[, -'cell_type']),
	breast_qian = quote(fread('../data_and_figures/qian_breast_2020.csv')[, -'cell_type']),
	crc_li = quote(fread('../data_and_figures/li_colorectal_2017.csv')[sample_type != 'normal', -c('cell_type', 'sample_type', 'epithelial_cluster')]),
	crc_lee_smc = quote(fread('../data_and_figures/lee_crc_2020_smc.csv')[sample_type != 'normal', -c('sample_id', 'sample_type', 'cell_type', 'cell_subtype')]),
	crc_lee_kul3 = quote(fread('../data_and_figures/lee_crc_2020_kul3.csv')[sample_type != 'normal', -c('sample_id', 'sample_type', 'cell_type', 'cell_subtype')]),
	crc_qian = quote(fread('../data_and_figures/qian_crc_2020.csv')[sample_type != 'normal', -c('sample_type', 'cell_type')]),
	hnscc_puram = quote(fread('../data_and_figures/puram_hnscc_2017.csv')[, -c('cell_type', 'lymph_node', 'processed_by_maxima_enzyme')]),
	liver_ma = quote(fread('../data_and_figures/ma_liver_2019.csv')[, -c('cell_type', 'sample', 'disease')]),
	liver_ma_hcc = quote(fread('../data_and_figures/ma_liver_2019.csv')[disease == 'HCC', -c('cell_type', 'sample', 'disease')]),
	liver_ma_icca = quote(fread('../data_and_figures/ma_liver_2019.csv')[disease == 'iCCA', -c('cell_type', 'sample', 'disease')]),
	lung_lambrechts = quote(fread('../data_and_figures/lambrechts_nsclc_2018.csv')[
		sample_type != 'normal',
		-c('cell_type', 'cluster', 'disease', 'sample_type', 'annotation')
	]),
	lung_lambrechts_luad = quote(fread('../data_and_figures/lambrechts_nsclc_2018.csv')[
		sample_type != 'normal' & disease == 'LUAD',
		-c('cell_type', 'cluster', 'disease', 'sample_type', 'annotation')
	]),
	lung_lambrechts_lusc = quote(fread('../data_and_figures/lambrechts_nsclc_2018.csv')[
		sample_type != 'normal' & disease == 'LUSC',
		-c('cell_type', 'cluster', 'disease', 'sample_type', 'annotation')
	]),
	lung_qian = quote(fread('../data_and_figures/qian_lung_2020.csv')[sample_type != 'normal', -c('sample_type', 'disease', 'cell_type')]),
	lung_qian_luad = quote(fread('../data_and_figures/qian_lung_2020.csv')[sample_type != 'normal' & disease == 'LUAD', -c('sample_type', 'disease', 'cell_type')]),
	lung_qian_lusc = quote(fread('../data_and_figures/qian_lung_2020.csv')[sample_type != 'normal' & disease == 'LUSC', -c('sample_type', 'disease', 'cell_type')]),
	lung_song = quote(fread('../data_and_figures/song_nsclc_2019.csv')[sample_type == 'tumour', -c('sample_type', 'cell_type', 'disease')]),
	lung_song_luad = quote(fread('../data_and_figures/song_nsclc_2019.csv')[sample_type == 'tumour' & disease == 'LUAD', -c('sample_type', 'cell_type', 'disease')]),
	lung_song_lusc = quote(fread('../data_and_figures/song_nsclc_2019.csv')[sample_type == 'tumour' & disease == 'LUSC', -c('sample_type', 'cell_type', 'disease')]),
	ovarian_izar_ss2 = quote(fread('../data_and_figures/izar_ovarian_2020_ss2.csv')[, -c('cell_type', 'cluster')]),
	ovarian_izar_10x = quote(fread('../data_and_figures/izar_ovarian_2020_10x.csv')[, -c('barcode', 'sample_id', 'time', 'cell_type', 'cluster')]),
	ovarian_qian = quote(fread('../data_and_figures/qian_ovarian_2020.csv')[sample_type != 'normal', -c('sample_type', 'site', 'cell_type')]),
	pdac_elyada = quote(fread('../data_and_figures/elyada_pdac_2019.csv')[, -c('cell_type', 'cluster')]),
	pdac_peng = quote(fread('../data_and_figures/peng_pdac_2019.csv')[sample_type != 'normal', -c('cell_type', 'sample_type')])
)

cat('Calculating number of genes detected:\n')

genes_detected <- sapply(
	names(cohort_data),
	function(cohort) {
		
		cat(paste0('\t', cohort, '\n'))
		
		sc_data <- eval(cohort_data[[cohort]])
		
		setNames(apply(sc_data[, -c('id', 'patient')], 1, function(x) sum(x > 0)), sc_data$id)
		
	},
	simplify = FALSE,
	USE.NAMES = TRUE
)

cat('Saving...')

saveRDS(genes_detected, '../data_and_figures/sc_genes_detected.rds')

cat('Done!\nMaking plots:\n')

pdf('../data_and_figures/sc_genes_detected_plots.pdf')

for(cohort in names(genes_detected)) {
	
	cat(paste0('\t', cohort, '\n'))
	
	# I'd like to do density plots as well, but then I probably have to choose appropriate bandwidths separately for each cohort.
	
	genes_detected_lineplot <- qplot(
		x = 1:length(genes_detected[[cohort]]),
		y = sort(genes_detected[[cohort]]),
		geom = 'line'
	) +
		theme_minimal() +
		labs(
			title = cohort,
			subtitle = 'Genes detected',
			x = 'Sorted cells',
			y = 'Number of genes detected'
		)
	
	print(genes_detected_lineplot)
	
}

dev.off()
