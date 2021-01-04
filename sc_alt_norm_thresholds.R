library(data.table) # 1.12.8
library(ggplot2) # 3.3.0
library(Matrix) # 1.2.18
library(stringr) # 1.4.0
library(plyr) # 1.8.6
library(limma) # 3.42.2
library(magrittr) # 1.5
library(cowplot) # 1.0.0

source('general_functions.R')
source('sparse_matrix_functions.R')
source('sc_functions.R')





density_plots <- list(
	gene_ave = list(breast_qian = list(), crc_lee_smc = list(), luad_kim = list(), pdac_peng = list()),
	exp_level = list(breast_qian = list(), crc_lee_smc = list(), luad_kim = list(), pdac_peng = list())
)





# breast_qian:

sc_data <- readRDS('../data_and_figures/sc_alt_norm/qian_breast_2020_scran.rds')
sc_meta <- fread('../data_and_figures/sc_alt_norm/qian_breast_2020_meta.csv')
tpm_data <- fread('../data_and_figures/qian_breast_2020_reclassified.csv', showProgress = FALSE)[
    cell_type != 'ambiguous' & id != 'sc5rJUQ064_CCATGTCCATCCCATC',
    -c('cell_type_author', 'cell_type_lenient')
]

gene_ave_scran <- rowMeans(sc_data)
gene_ave_tpm <- log2(colMeans(10*(2^tpm_data[, -c('id', 'patient', 'cell_type')] - 1)) + 1)
gene_ave_scran_density <- density(gene_ave_scran)
gene_ave_tpm_density <- density(gene_ave_tpm)

# 0.07 would be a good threshold for gene expression, as it just shaves off the peak at zero.  However, to give a number of genes more similar to
# the threshold I used for the TPM/10 data (4.5 in this case), I'm using a higher threshold:
density_plots$gene_ave$breast_qian$scran <- ggplot(data.frame(x = gene_ave_scran_density$x, y = gene_ave_scran_density$y), aes(x = x, y = y)) +
	geom_line() +
	lims(x = c(NA, 1), y = c(NA, 10)) +
	geom_vline(xintercept = 0.14, colour = 'blue', linetype = 'dashed') +
	labs(x = 'gene averages', y = 'density', title = 'Breast, Qian et al., scran-normalised') +
	annotate('text', x = 0.8, y = 8, label = paste('Genes above\nthreshold:', sum(gene_ave_scran > 0.14))) +
	theme_minimal()

density_plots$gene_ave$breast_qian$tpm <- ggplot(data.frame(x = gene_ave_tpm_density$x, y = gene_ave_tpm_density$y), aes(x = x, y = y)) +
	geom_line() +
	lims(x = c(NA, 12), y = c(NA, 0.45)) +
	geom_vline(xintercept = 4.5, colour = 'blue', linetype = 'dashed') +
	labs(x = 'gene averages', y = 'density', title = 'Breast, Qian et al., TPM/10') +
	annotate('text', x = 9, y = 0.3, label = paste('Genes above\nthreshold:', sum(gene_ave_tpm > 4.5))) +
	theme_minimal()

exp_level_scran_density <- density(as.numeric(sc_data))
exp_level_tpm_density <- density(unlist(tpm_data[, -c('id', 'patient', 'cell_type')]))
exp_level_scran_pass <- sum(
	tapply(
		sc_data[, sc_meta[cell_type == 'cancer', id]]@x,
		sc_data[, sc_meta[cell_type == 'cancer', id]]@i,
		function(g) {sum(g >= 3.5) >= sc_meta[cell_type == 'cancer', .N/100]}
	)
)
exp_level_tpm_pass <- sum(sapply(tpm_data[cell_type == 'cancer', -c('id', 'patient', 'cell_type')], function(x) {sum(x >= 7) >= length(x)/100}))

# The following plots suggest that 3.5 could be a good cut-off for "very high expression", and is comparable to using 7 in TPM/10 data:
density_plots$exp_level$breast_qian$scran <- ggplot(
	data.frame(
		x = exp_level_scran_density$x,
		y = replace(exp_level_scran_density$y, exp_level_scran_density$y >= 0.2, 0.2)
	),
	aes(x = x, y = y)
) +
	geom_line() +
	lims(x = c(NA, 13), y = c(NA, 0.2)) +
	geom_vline(xintercept = 3.5, colour = 'blue', linetype = 'dashed') +
	labs(x = 'expression level', y = 'density', title = 'Breast, Qian et al., scran-normalised') +
	annotate(
		'text',
		x = 10,
		y = 0.15,
		label = paste0('Genes above\nthreshold in at\nleast 1% of\ncancer cells:\n', exp_level_scran_pass)
	) +
	theme_minimal()

density_plots$exp_level$breast_qian$tpm <- ggplot(
	data.frame(
		x = exp_level_tpm_density$x,
		y = replace(exp_level_tpm_density$y, exp_level_tpm_density$y >= 0.2, 0.2)
	),
	aes(x = x, y = y)
) +
	geom_line() +
	lims(x = c(NA, 13), y = c(NA, 0.2)) +
	geom_vline(xintercept = 7, colour = 'blue', linetype = 'dashed') +
	labs(x = 'expression level', y = 'density', title = 'Breast, Qian et al., TPM/10') +
	annotate(
		'text',
		x = 10,
		y = 0.15,
		label = paste0('Genes above\nthreshold in at\nleast 1% of\ncancer cells:\n', exp_level_tpm_pass)
	) +
	theme_minimal()





# crc_lee_smc:

sc_data <- readRDS('../data_and_figures/sc_alt_norm/lee_crc_2020_smc_scran.rds')
sc_meta <- fread('../data_and_figures/sc_alt_norm/lee_crc_2020_smc_meta.csv')
tpm_data <- fread('../data_and_figures/lee_crc_2020_smc_reclassified.csv')[cell_type != 'ambiguous', -c('cell_type_author', 'cell_type_lenient')]

gene_ave_scran <- rowMeans(sc_data)
gene_ave_tpm <- log2(colMeans(10*(2^tpm_data[, -c('id', 'patient', 'cell_type')] - 1)) + 1)
gene_ave_scran_density <- density(gene_ave_scran)
gene_ave_tpm_density <- density(gene_ave_tpm)

density_plots$gene_ave$crc_lee_smc$scran <- ggplot(data.frame(x = gene_ave_scran_density$x, y = gene_ave_scran_density$y), aes(x = x, y = y)) +
	geom_line() +
	lims(x = c(NA, 1), y = c(NA, 10)) +
	geom_vline(xintercept = 0.17, colour = 'blue', linetype = 'dashed') +
	labs(x = 'gene averages', y = 'density', title = 'CRC - SMC, Lee et al., scran-normalised') +
	annotate('text', x = 0.8, y = 8, label = paste('Genes above\nthreshold:', sum(gene_ave_scran > 0.17))) +
	theme_minimal()

density_plots$gene_ave$crc_lee_smc$tpm <- ggplot(data.frame(x = gene_ave_tpm_density$x, y = gene_ave_tpm_density$y), aes(x = x, y = y)) +
	geom_line() +
	lims(x = c(NA, 12), y = c(NA, 0.45)) +
	geom_vline(xintercept = 4, colour = 'blue', linetype = 'dashed') +
	labs(x = 'gene averages', y = 'density', title = 'CRC - SMC, Lee et al., TPM/10') +
	annotate('text', x = 9, y = 0.3, label = paste('Genes above\nthreshold:', sum(gene_ave_tpm > 4))) +
	theme_minimal()

exp_level_scran_density <- density(as.numeric(sc_data))
exp_level_tpm_density <- density(unlist(tpm_data[, -c('id', 'patient', 'cell_type')]))
exp_level_scran_pass <- sum(
	tapply(
		sc_data[, sc_meta[cell_type == 'cancer', id]]@x,
		sc_data[, sc_meta[cell_type == 'cancer', id]]@i,
		function(g) {sum(g >= 4.25) >= sc_meta[cell_type == 'cancer', .N/100]}
	)
)
exp_level_tpm_pass <- sum(sapply(tpm_data[cell_type == 'cancer', -c('id', 'patient', 'cell_type')], function(x) {sum(x >= 7) >= length(x)/100}))

density_plots$exp_level$crc_lee_smc$scran <- ggplot(
	data.frame(
		x = exp_level_scran_density$x,
		y = replace(exp_level_scran_density$y, exp_level_scran_density$y >= 0.2, 0.2)
	),
	aes(x = x, y = y)
) +
	geom_line() +
	lims(x = c(NA, 13), y = c(NA, 0.2)) +
	geom_vline(xintercept = 4.25, colour = 'blue', linetype = 'dashed') +
	labs(x = 'expression level', y = 'density', title = 'CRC - SMC, Lee et al., scran-normalised') +
	annotate(
		'text',
		x = 10,
		y = 0.15,
		label = paste0('Genes above\nthreshold in at\nleast 1% of\ncancer cells:\n', exp_level_scran_pass)
	) +
	theme_minimal()

density_plots$exp_level$crc_lee_smc$tpm <- ggplot(
	data.frame(
		x = exp_level_tpm_density$x,
		y = replace(exp_level_tpm_density$y, exp_level_tpm_density$y >= 0.2, 0.2)
	),
	aes(x = x, y = y)
) +
	geom_line() +
	lims(x = c(NA, 13), y = c(NA, 0.2)) +
	geom_vline(xintercept = 7, colour = 'blue', linetype = 'dashed') +
	labs(x = 'expression level', y = 'density', title = 'CRC - SMC, Lee et al., TPM/10') +
	annotate(
		'text',
		x = 10,
		y = 0.15,
		label = paste0('Genes above\nthreshold in at\nleast 1% of\ncancer cells:\n', exp_level_tpm_pass)
	) +
	theme_minimal()





# luad_kim:

sc_data <- readRDS('../data_and_figures/sc_alt_norm/kim_luad_2020_scran.rds')
sc_meta <- fread('../data_and_figures/sc_alt_norm/kim_luad_2020_meta.csv')
tpm_data <- fread('../data_and_figures/kim_luad_2020_reclassified.csv')[cell_type != 'ambiguous', -'cell_type_author']

gene_ave_scran <- rowMeans(sc_data)
gene_ave_tpm <- log2(colMeans(10*(2^tpm_data[, -c('id', 'patient', 'cell_type')] - 1)) + 1)
gene_ave_scran_density <- density(gene_ave_scran)
gene_ave_tpm_density <- density(gene_ave_tpm)

density_plots$gene_ave$luad_kim$scran <- ggplot(data.frame(x = gene_ave_scran_density$x, y = gene_ave_scran_density$y), aes(x = x, y = y)) +
	geom_line() +
	lims(x = c(NA, 1), y = c(NA, 10)) +
	geom_vline(xintercept = 0.12, colour = 'blue', linetype = 'dashed') +
	labs(x = 'gene averages', y = 'density', title = 'LUAD, Kim et al., scran-normalised') +
	annotate('text', x = 0.8, y = 8, label = paste('Genes above\nthreshold:', sum(gene_ave_scran > 0.12))) +
	theme_minimal()

density_plots$gene_ave$luad_kim$tpm <- ggplot(data.frame(x = gene_ave_tpm_density$x, y = gene_ave_tpm_density$y), aes(x = x, y = y)) +
	geom_line() +
	lims(x = c(NA, 12), y = c(NA, 0.45)) +
	geom_vline(xintercept = 4, colour = 'blue', linetype = 'dashed') +
	labs(x = 'gene averages', y = 'density', title = 'LUAD, Kim et al., TPM/10') +
	annotate('text', x = 9, y = 0.3, label = paste('Genes above\nthreshold:', sum(gene_ave_tpm > 4))) +
	theme_minimal()

exp_level_scran_density <- density(as.numeric(sc_data))
exp_level_tpm_density <- density(unlist(tpm_data[, -c('id', 'patient', 'cell_type')]))
exp_level_scran_pass <- sum(
	tapply(
		sc_data[, sc_meta[cell_type == 'cancer', id]]@x,
		sc_data[, sc_meta[cell_type == 'cancer', id]]@i,
		function(g) {sum(g >= 3.5) >= sc_meta[cell_type == 'cancer', .N/100]}
	)
)
exp_level_tpm_pass <- sum(sapply(tpm_data[cell_type == 'cancer', -c('id', 'patient', 'cell_type')], function(x) {sum(x >= 7) >= length(x)/100}))

density_plots$exp_level$luad_kim$scran <- ggplot(
	data.frame(
		x = exp_level_scran_density$x,
		y = replace(exp_level_scran_density$y, exp_level_scran_density$y >= 0.2, 0.2)
	),
	aes(x = x, y = y)
) +
	geom_line() +
	lims(x = c(NA, 13), y = c(NA, 0.2)) +
	geom_vline(xintercept = 3.5, colour = 'blue', linetype = 'dashed') +
	labs(x = 'expression level', y = 'density', title = 'LUAD, Kim et al., scran-normalised') +
	annotate(
		'text',
		x = 10,
		y = 0.15,
		label = paste0('Genes above\nthreshold in at\nleast 1% of\ncancer cells:\n', exp_level_scran_pass)
	) +
	theme_minimal()

density_plots$exp_level$luad_kim$tpm <- ggplot(
	data.frame(
		x = exp_level_tpm_density$x,
		y = replace(exp_level_tpm_density$y, exp_level_tpm_density$y >= 0.2, 0.2)
	),
	aes(x = x, y = y)
) +
	geom_line() +
	lims(x = c(NA, 13), y = c(NA, 0.2)) +
	geom_vline(xintercept = 7, colour = 'blue', linetype = 'dashed') +
	labs(x = 'expression level', y = 'density', title = 'LUAD, Kim et al., TPM/10') +
	annotate(
		'text',
		x = 10,
		y = 0.15,
		label = paste0('Genes above\nthreshold in at\nleast 1% of\ncancer cells:\n', exp_level_tpm_pass)
	) +
	theme_minimal()





# pdac_peng:

sc_data <- readRDS('../data_and_figures/sc_alt_norm/peng_pdac_2019_scran.rds')
sc_meta <- fread('../data_and_figures/sc_alt_norm/peng_pdac_2019_meta.csv')
tpm_data <- fread('../data_and_figures/peng_pdac_2019_reclassified.csv')[
	cell_type != 'ambiguous' & !(id %in% c('T8_TGGTTCCTCGCATGGC', 'T17_CGTGTAACAGTACACT')),
	-'cell_type_author'
]

gene_ave_scran <- rowMeans(sc_data)
gene_ave_tpm <- log2(colMeans(10*(2^tpm_data[, -c('id', 'patient', 'cell_type')] - 1)) + 1)
gene_ave_scran_density <- density(gene_ave_scran)
gene_ave_tpm_density <- density(gene_ave_tpm)

density_plots$gene_ave$pdac_peng$scran <- ggplot(data.frame(x = gene_ave_scran_density$x, y = gene_ave_scran_density$y), aes(x = x, y = y)) +
	geom_line() +
	lims(x = c(NA, 1), y = c(NA, 10)) +
	geom_vline(xintercept = 0.18, colour = 'blue', linetype = 'dashed') +
	labs(x = 'gene averages', y = 'density', title = 'PDAC, Peng et al., scran-normalised') +
	annotate('text', x = 0.8, y = 8, label = paste('Genes above\nthreshold:', sum(gene_ave_scran > 0.18))) +
	theme_minimal()

density_plots$gene_ave$pdac_peng$tpm <- ggplot(data.frame(x = gene_ave_tpm_density$x, y = gene_ave_tpm_density$y), aes(x = x, y = y)) +
	geom_line() +
	lims(x = c(NA, 12), y = c(NA, 0.45)) +
	geom_vline(xintercept = 4.5, colour = 'blue', linetype = 'dashed') +
	labs(x = 'gene averages', y = 'density', title = 'PDAC, Peng et al., TPM/10') +
	annotate('text', x = 9, y = 0.3, label = paste('Genes above\nthreshold:', sum(gene_ave_tpm > 4.5))) +
	theme_minimal()

exp_level_scran_density <- density(as.numeric(sc_data))
exp_level_tpm_density <- density(unlist(tpm_data[, -c('id', 'patient', 'cell_type')]))
exp_level_scran_pass <- sum(
	tapply(
		sc_data[, sc_meta[cell_type == 'cancer', id]]@x,
		sc_data[, sc_meta[cell_type == 'cancer', id]]@i,
		function(g) {sum(g >= 3.75) >= sc_meta[cell_type == 'cancer', .N/100]}
	)
)
exp_level_tpm_pass <- sum(sapply(tpm_data[cell_type == 'cancer', -c('id', 'patient', 'cell_type')], function(x) {sum(x >= 7) >= length(x)/100}))

density_plots$exp_level$pdac_peng$scran <- ggplot(
	data.frame(
		x = exp_level_scran_density$x,
		y = replace(exp_level_scran_density$y, exp_level_scran_density$y >= 0.2, 0.2)
	),
	aes(x = x, y = y)
) +
	geom_line() +
	lims(x = c(NA, 13), y = c(NA, 0.2)) +
	geom_vline(xintercept = 3.75, colour = 'blue', linetype = 'dashed') +
	labs(x = 'expression level', y = 'density', title = 'PDAC, Peng et al., scran-normalised') +
	annotate(
		'text',
		x = 10,
		y = 0.15,
		label = paste0('Genes above\nthreshold in at\nleast 1% of\ncancer cells:\n', exp_level_scran_pass)
	) +
	theme_minimal()

density_plots$exp_level$pdac_peng$tpm <- ggplot(
	data.frame(
		x = exp_level_tpm_density$x,
		y = replace(exp_level_tpm_density$y, exp_level_tpm_density$y >= 0.2, 0.2)
	),
	aes(x = x, y = y)
) +
	geom_line() +
	lims(x = c(NA, 13), y = c(NA, 0.2)) +
	geom_vline(xintercept = 7, colour = 'blue', linetype = 'dashed') +
	labs(x = 'expression level', y = 'density', title = 'PDAC, Peng et al., TPM/10') +
	annotate(
		'text',
		x = 10,
		y = 0.15,
		label = paste0('Genes above\nthreshold in at\nleast 1% of\ncancer cells:\n', exp_level_tpm_pass)
	) +
	theme_minimal()





pdf('../data_and_figures/sc_alt_norm/gene_thresholds.pdf', width = 8, height = 12)
print(plot_grid(plotlist = unlist(density_plots$gene_ave, recursive = FALSE), nrow = 4, ncol = 2, align = 'hv'))
print(plot_grid(plotlist = unlist(density_plots$exp_level, recursive = FALSE), nrow = 4, ncol = 2, align = 'hv'))
dev.off()
