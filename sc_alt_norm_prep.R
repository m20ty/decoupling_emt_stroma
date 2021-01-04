library(data.table) # 1.12.8
library(Matrix) # 1.2.18
library(stringr) # 1.4.0
library(plyr) # 1.8.6
library(limma) # 3.42.2
library(magrittr) # 1.5
library(cowplot) # 1.0.0
library(scran) # 1.14.6
library(scater) # 1.14.6
library(SingleCellExperiment) # 1.8.0

source('general_functions.R')
source('sparse_matrix_functions.R')
source('sc_functions.R')





# Breast:

cat('breast_qian\n')

sc_data <- readMM('../../single_cell_data/qian_2020/BC_counts/matrix.mtx')
barcodes <- fread('../../single_cell_data/qian_2020/BC_counts/barcodes.tsv', header = FALSE)$V1
genes <- fread('../../single_cell_data/qian_2020/BC_counts/genes.tsv', header = FALSE)$V1
gn <- alias2SymbolTable(genes)
sc_data <- sc_data[!is.na(gn) & gn %in% names(table(gn))[table(gn) == 1], ]
rownames(sc_data) <- gn[!is.na(gn) & gn %in% names(table(gn))[table(gn) == 1]]
colnames(sc_data) <- barcodes
tpm_data <- fread('../data_and_figures/qian_breast_2020_reclassified.csv')[cell_type != 'ambiguous' & id != 'sc5rJUQ064_CCATGTCCATCCCATC', -c('cell_type_author', 'cell_type_lenient')]

sc_data <- sc_data[, tpm_data$id]
sc_meta <- tpm_data[, .(id, patient, cell_type)]

rm(tpm_data)
rm(barcodes)
rm(genes)
rm(gn)

fwrite(sc_meta, '../data_and_figures/sc_alt_norm/qian_breast_2020_meta.csv')

sc_data_scran <- SingleCellExperiment(list(counts = sc_data))
scran_clusters <- quickCluster(sc_data_scran)
sc_data_scran <- computeSumFactors(sc_data_scran, clusters = scran_clusters)
sc_data_scran <- logNormCounts(sc_data_scran)
saveRDS(SingleCellExperiment::logcounts(sc_data_scran), '../data_and_figures/sc_alt_norm/qian_breast_2020_scran.rds')





# CRC:

cat('crc_lee_smc\n')

sc_data <- readRDS('../../single_cell_data/lee_crc_2020/GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.rds')
gn <- alias2SymbolTable(rownames(sc_data))
sc_data <- sc_data[!is.na(gn) & gn %in% names(table(gn))[table(gn) == 1], ]
rownames(sc_data) <- gn[!is.na(gn) & gn %in% names(table(gn))[table(gn) == 1]]
tpm_data <- fread('../data_and_figures/lee_crc_2020_smc_reclassified.csv')[cell_type != 'ambiguous', -c('cell_type_author', 'cell_type_lenient')]

sc_data <- sc_data[, tpm_data$id]
sc_meta <- tpm_data[, .(id, patient, cell_type)]

rm(tpm_data)
rm(gn)

fwrite(sc_meta, '../data_and_figures/sc_alt_norm/lee_crc_2020_smc_meta.csv')

sc_data_scran <- SingleCellExperiment(list(counts = sc_data))
scran_clusters <- quickCluster(sc_data_scran)
sc_data_scran <- computeSumFactors(sc_data_scran, clusters = scran_clusters)
sc_data_scran <- logNormCounts(sc_data_scran)
saveRDS(SingleCellExperiment::logcounts(sc_data_scran), '../data_and_figures/sc_alt_norm/lee_crc_2020_smc_scran.rds')





# LUAD:

cat('luad_kim\n')

sc_data <- readRDS('../../single_cell_data/kim_luad_2020/kim_luad_2020_matrix_counts.rds')
tpm_data <- fread('../data_and_figures/kim_luad_2020_reclassified.csv')[cell_type != 'ambiguous', -'cell_type_author']

sc_data <- sc_data[, tpm_data$id]
sc_meta <- tpm_data[, .(id, patient, cell_type)]

rm(tpm_data)

fwrite(sc_meta, '../data_and_figures/sc_alt_norm/kim_luad_2020_meta.csv')

sc_data_scran <- SingleCellExperiment(list(counts = sc_data))
scran_clusters <- quickCluster(sc_data_scran)
sc_data_scran <- computeSumFactors(sc_data_scran, clusters = scran_clusters)
sc_data_scran <- logNormCounts(sc_data_scran)
saveRDS(SingleCellExperiment::logcounts(sc_data_scran), '../data_and_figures/sc_alt_norm/kim_luad_2020_scran.rds')





# PDAC:

cat('pdac_peng\n')

sc_data <- readRDS('../../single_cell_data/peng_pdac_2019/count-matrix.rds')
sc_data <- sc_data[!(rownames(sc_data) %in% c('AREGB', 'DDC8')), ]
gn <- alias2SymbolTable(rownames(sc_data))
sc_data <- sc_data[!is.na(gn) & gn %in% names(table(gn))[table(gn) == 1], ]
rownames(sc_data) <- gn[!is.na(gn) & gn %in% names(table(gn))[table(gn) == 1]]
tpm_data <- fread('../data_and_figures/peng_pdac_2019_reclassified.csv')[cell_type != 'ambiguous' & !(id %in% c('T8_TGGTTCCTCGCATGGC', 'T17_CGTGTAACAGTACACT')), -'cell_type_author']

sc_data <- sc_data[, tpm_data$id]
sc_meta <- tpm_data[, .(id, patient, cell_type)]

rm(gn)
rm(tpm_data)

fwrite(sc_meta, '../data_and_figures/sc_alt_norm/peng_pdac_2019_meta.csv')

sc_data_scran <- SingleCellExperiment(list(counts = sc_data))
scran_clusters <- quickCluster(sc_data_scran)
sc_data_scran <- computeSumFactors(sc_data_scran, clusters = scran_clusters)
sc_data_scran <- logNormCounts(sc_data_scran)
saveRDS(SingleCellExperiment::logcounts(sc_data_scran), '../data_and_figures/sc_alt_norm/peng_pdac_2019_scran.rds')
