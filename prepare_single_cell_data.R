library(data.table) # 1.12.8
library(Matrix) # 1.2.18
library(stringr) # 1.4.0
library(plyr) # 1.8.6
library(limma) # 3.42.2
library(AnnotationDbi) # 1.48.0
library(org.Hs.eg.db) # 3.10.0
library(magrittr) # 1.5

source('general_functions.R')
source('sparse_matrix_functions.R')





# Lung dataset from Song et al., 2019:

tumour_file_names <- c(
    '../../single_cell_data/song_nsclc_2019/GSM3304007_P1_Tumor_processed_data.txt',
	'../../single_cell_data/song_nsclc_2019/GSM3304008_P1_Normal_processed_data.txt',
    '../../single_cell_data/song_nsclc_2019/GSM3304009_P2_Tumor_processed_data.txt',
	'../../single_cell_data/song_nsclc_2019/GSM3304010_P2_Normal_processed_data.txt',
    '../../single_cell_data/song_nsclc_2019/GSM3304011_P3_Tumor_processed_data.txt',
	'../../single_cell_data/song_nsclc_2019/GSM3304012_P3_Normal_processed_data.txt',
    '../../single_cell_data/song_nsclc_2019/GSM3304013_P4_Tumor_processed_data.txt',
	'../../single_cell_data/song_nsclc_2019/GSM3304014_P4_Normal_processed_data.txt'
)

sc_data <- rbindlist(
    lapply(
        tumour_file_names,
        function(fn) {

            # Read in data file:
            dt <- fread(fn)

            # Transpose and add patient IDs and sample type:
			dt <- tdt(dt)
            dt[
				,
				c('patient', 'sample_type') := .(
					str_extract(str_extract(fn, 'P[0-9]'), '[0-9]'),
					mapvalues(str_extract(fn, 'Tumor|Normal'), c('Tumor', 'Normal'), c('tumour', 'normal'), warn_missing = FALSE)
				)
			]

			setcolorder(dt, c('id', 'patient', 'sample_type'))

            # Output:
            dt

        }
    ),
    fill = TRUE
)

sc_data[is.na(sc_data)] <- 0

# Update gene symbols using alias2SymbolTable() and remove resulting NAs and repeated genes:
sample_data <- sc_data[, .(id, patient, sample_type)]
sc_data <- tdt(sc_data[, -c('patient', 'sample_type')])
sc_data[, id := alias2SymbolTable(id)]
sc_data <- tdt(sc_data[!is.na(id) & id %in% names(table(id))[table(id) == 1]])

# In the following, we can't use merge() because there are repeated cell IDs.  We can use cbind() because the cell IDs are in the same order in
# sample_data and sc_data.
sc_data <- cbind(sample_data, sc_data[, -'id'])

rm(sample_data)

# Filter out cells of low complexity:
sc_data <- sc_data[apply(sc_data[, -c('id', 'patient', 'sample_type')], 1, function(x) sum(x > 0)) >= 1000]

# The following shows that all the IDs end in '-1', so we can remove this suffix.  This is important because the cell IDs in the annotations table do
# not have this suffix, so they won't match unless we remove it.
# sc_data[, unique(str_split_fixed(id, '-', 2)[, 2])]

sc_data[, id := str_split_fixed(id, '-', 2)[, 1]]

tumour_annotations_file_names <- paste0('../../single_cell_data/song_nsclc_2019/', dir('../../single_cell_data/song_nsclc_2019'))
tumour_annotations_file_names <- tumour_annotations_file_names[grepl('Annotation', tumour_annotations_file_names)]

sc_meta <- rbindlist(
    lapply(
        tumour_annotations_file_names,
        function(fn) {
			dt <- fread(fn, header = TRUE)
            dt[
				,
				c('V1', 'patient', 'sample_type') := .(
					NULL,
					str_extract(str_extract(fn, 'P[0-9]'), '[0-9]'),
					mapvalues(str_extract(fn, 'Tumor|Normal'), c('Tumor', 'Normal'), c('tumour', 'normal'), warn_missing = FALSE)
				)
			]
            dt
        }
    )
)

# There are a few duplicate cell IDs.  We'll deal with them by adjoining the patient number as a suffix.
sc_data[id %in% names(table(id))[table(id) > 1], id := paste(id, patient, sep = '_')]
sc_meta[cells %in% names(table(cells))[table(cells) > 1], cells := paste(cells, patient, sep = '_')]

# Include cell type column:
setkey(sc_data, id)
setkey(sc_meta, cells)
sc_data[, cell_type := sc_meta[id, type]]

# Read in series matrix to get disease type:
series_mat <- tdt(fread('../../single_cell_data/song_nsclc_2019/GSE117570_series_matrix.txt', skip = 31, header = FALSE))[, id := NULL]

series_mat <- series_mat[
    ,
    .(
        patient_number = gsub('P', '', str_split_fixed(`!Sample_title`, '_', 2)[, 1]),
        disease_code = mapvalues(
            gsub('diagnosis: ', '', `!Sample_characteristics_ch1`),
            c('Lung adenocarcinoma', 'Lung squamous cell carcinoma'),
            c('LUAD', 'LUSC')
        )
    )
] %>% unique

setkey(series_mat, patient_number)
sc_data[, disease := series_mat[patient, disease_code]]
setcolorder(sc_data, c('id', 'patient', 'sample_type', 'cell_type', 'disease'))

# Convert to TPM (actually counts per 100,000), take logs and round to 4 decimal places:
gene_cols <- names(sc_data[, -c('id', 'patient', 'sample_type', 'cell_type', 'disease')])
sc_data[, (gene_cols) := transpose(lapply(transpose(.SD), function(x) round(log2(1e+05*x/sum(x) + 1), 4))), .SDcols = gene_cols]

# It might actually be better to use 1e+04 in place of 1e+05, since the maximum total expression in any cell is still just less than 60000, so 10000
# may be a more reasonable scale.  Actually, it seems this isn't the first time I've thought this: I wrote a comment on it in
# lung_Lambrechts_sc-RNAseq_data_prep.R.

fwrite(sc_data, '../data_and_figures/song_nsclc_2019.csv')





# Lung dataset from Lambrechts et al., 2018:

sc_data <- readRDS('../../single_cell_data/lambrechts_nsclc_2018/RawDataLung.table.rds') # Note this is a sparse matrix
sc_meta <- as.data.table(readxl::read_xlsx('../../single_cell_data/lambrechts_nsclc_2018/MetaData.xlsx')[, -1])

# Update gene names using alias2SymbolTable() and remove resulting NAs and repeated genes:
gene_ids <- limma::alias2SymbolTable(rownames(sc_data))
rownames(sc_data) <- gene_ids
sc_data <- sc_data[!is.na(gene_ids) & gene_ids %in% names(table(gene_ids))[table(gene_ids) == 1], ]

# Convert to dgTMatrix, which uses x, i and j instead of x, i and p:
sc_data <- as(sc_data, 'dgTMatrix')

# Find number of genes detected per cell, which is just the length of x for each j:
genes_detected <- tapply(sc_data@x, sc_data@j, length)

sc_data <- sc_data[, which(genes_detected >= 1000)] # Doesn't work without which()!

sc_data <- as.data.table(as.matrix(sc_data), keep.rownames = 'id')

# Convert to "TPM/10", take logs and round to 4 decimal places:
sc_data[, names(sc_data[, -'id']) := lapply(.SD, function(x) {round(log2(1e+05*x/sum(x) + 1), 4)}), .SDcols = -'id']

sc_data <- tdt(sc_data)

setkey(sc_meta, cell)

# In the following, we add a 'disease' column, the information for which can be found in the paper in Fig. 1A and Supplementary Table 1.  Patients 1 &
# 2 are LUSC; patients 3 & 4 are LUAD; and patient 5 is ambiguous - in Fig. 1A they call it simply "NSCLC", while in the supplementary table they call
# it "Large cell" and describe it as "poorly differentiated morphology, negative for TTF1 and p63".

sc_data[
    ,
    c('patient', 'cell_type', 'cluster', 'disease', 'sample_type', 'annotation') := sc_meta[
        id,
        .(
            `PatientNumber MS`,
            mapvalues(tolower(CellType), c('ec', 'fibro', 'tumor', 'epi'), c('endothelial', 'fibroblast', 'cancer', 'epithelial')),
            cluster,
            mapvalues(`PatientNumber MS`, 1:5, c('LUSC', 'LUSC', 'LUAD', 'LUAD', 'NSCLC')),
            str_split_fixed(Patient_piece, '_', 2)[, 2],
            Annotation
        )
    ]
]

setcolorder(sc_data, c('id', 'patient', 'cell_type', 'cluster', 'disease', 'sample_type', 'annotation'))

fwrite(sc_data, '../data_and_figures/lambrechts_nsclc_2018.csv')





# BRCA data from Chung et al:

# While the Ensembl IDs are supplied here, we opt for using alias2SymbolTable() on the gene names, since the Ensembl IDs actually give more NAs when
# converted into symbols, while the non-NA elements don't give any more useful genes.  We then use mapIds() only on the genes that end up as repeats
# after alias2SymbolTable().  This is useful because the repeated genes include some that may be important for EMT, such as AREG, ITGB3 and TIMP2.
# Genes which still occur more than once (or become NA) after this procedure are removed.

sc_data <- fread(paste0('../../single_cell_data/chung_breast_cancer_2017/', 'GSE75688_GEO_processed_Breast_Cancer_raw_TPM_matrix.txt'))[
    gene_type %in% c('protein_coding', 'pseudogene')
][, gene_name := limma::alias2SymbolTable(gene_name)][
    gene_name %in% names(table(gene_name))[table(gene_name) > 1],
    gene_name := mapIds(
        org.Hs.eg.db,
        stringr::str_split_fixed(gene_id, '\\.', 2)[, 1], # Don't need version number
        keytype = 'ENSEMBL',
        column = 'SYMBOL'
    )
][!is.na(gene_name) & gene_name %in% names(table(gene_name))[table(gene_name) == 1]][, c('gene_id', 'gene_type') := NULL]

# Remove the bulk samples, and also the ones labelled 'Tumor', which I don't trust (they don't appear in the metadata, so they could also be bulk
# samples - the numbers for these samples look a bit like bulk):
sc_data[, names(sc_data)[endsWith(names(sc_data), 'Pooled') | endsWith(names(sc_data), 'Tumor')] := NULL]

# Filter out cells with less than 1000 genes detected (there aren't many):
sc_data <- sc_data[, c(TRUE, apply(sc_data[, -'gene_name'], 2, function(x) sum(x > 0)) >= 1000), with = FALSE]

sc_data[, names(sc_data[, -'gene_name']) := lapply(.SD, function(x) {round(log2(1e+05*x/sum(x) + 1), 4)}), .SDcols = -'gene_name']

sc_data <- tdt(sc_data)

# Add metadata:

sc_meta <- fread('../../single_cell_data/chung_breast_cancer_2017/GSE75688_final_sample_information.txt', key = 'sample')

sc_data[
    ,
    c('patient', 'cell_type') := .(
        as.numeric(gsub('BC', '', gsub('LN', '', str_split_fixed(id, '_', 2)[, 1]))),
        mapvalues(tolower(sc_meta[id, index3]), c('tumor', 'stromal', 'tcell', 'bcell'), c('cancer', 'fibroblast', 't_cell', 'b_cell'))
    )
][, lymph_node := switch(endsWith(str_split_fixed(id, '_', 2)[, 1], 'LN') + 1, 0, 1), by = id]

setcolorder(sc_data, c('id', 'patient', 'cell_type', 'lymph_node'))

fwrite(sc_data, '../data_and_figures/chung_breast_cancer_2017.csv')





# BRCA (TNBC) data from Karaayvaz et al.:

sc_data <- fread('../../single_cell_data/karaayvaz_tnbc_2018/GSE118389_tpm_rsem.txt', key = 'V1')

# Update gene names:
sc_data[, V1 := alias2SymbolTable(V1)]
sc_data <- sc_data[!is.na(V1) & !(V1 %in% names(table(V1))[table(V1) > 1])]

sc_data <- tdt(sc_data)
sc_data[, patient := stringr::str_split_fixed(id, '_', 2)[, 1]]

# The data on cell type assignments can be obtained from Github:
# download.file(
#     # 'https://github.com/Michorlab/tnbc_scrnaseq/blob/master/data/cell_types_tab_S9.txt',
#     'https://raw.githubusercontent.com/Michorlab/tnbc_scrnaseq/master/data/cell_types_tab_S9.txt',
#     destfile = '../../single_cell_data/karaayvaz_tnbc_2018/cell_types_tab_S9.txt'
# )

# There are only 1112 cells in this file, fewer than in the counts matrix, but this is because these are the ones that survived their two-step
# clustering process, and which they presumably therefore have confidence in assigning to cell types.  See pp. 35-37 of the supplementary information
# file from the paper (in particular, it says at the bottom of p. 37 that they're left with 1112 cells).

cell_types_data <- fread('../../single_cell_data/karaayvaz_tnbc_2018/cell_types_tab_S9.txt')[, 2:3] %>% setNames(c('cell_id', 'cell_type'))
setkey(cell_types_data, cell_id)
sc_data[, cell_type := cell_types_data[id, cell_type]]
setcolorder(sc_data, c('id', 'patient', 'cell_type'))

# Filter out cells with fewer than 1000 genes detected:
sc_data <- sc_data[apply(sc_data[, -c('id', 'patient', 'cell_type')], 1, function(x) sum(x > 0)) >= 1000]

# Remove the cells that seem to have considerably lower total gene expression than the others (despite the data supposedly being TPM):
sc_data <- sc_data[rowSums(sc_data[, -c('id', 'patient', 'cell_type')]) >= 9e+05]

# Scale to TPM/10, take logs and round to 4 decimal places:
sc_data[
    ,
    names(sc_data[, -c('id', 'patient', 'cell_type')]) := transpose(lapply(transpose(.SD), function(x) {round(log2(1e+05*x/sum(x) + 1), 4)})),
    .SDcols = -c('id', 'patient', 'cell_type')
]

fwrite(sc_data, '../data_and_figures/karaayvaz_tnbc_2018.csv')





# Colorectal data from Li et al.:

# I'm using the same strategy for updating gene names as I used in the Chung et al. Breast dataset.  The method of dealing with repeats is important
# because I have potentially EMT-relevant genes among the repeats including AREG, ITGB3 and TIMP2, as before.

sc_data <- merge(
    fread('../../single_cell_data/li_colorectal_2017/GSE81861_CRC_tumor_all_cells_FPKM.csv'),
    fread('../../single_cell_data/li_colorectal_2017/GSE81861_CRC_NM_all_cells_FPKM.csv'),
    by = 'V1'
)[
    ,
    c('ensembl_id', 'symbol') := .(
        str_split_fixed(gsub('.+_ENSG', 'ENSG', V1), '\\.', 2)[, 1],
        gsub('^chr[A-Z]*[0-9]*:[0-9]+-[0-9]+_', '', gsub('_ENSG[A-Z]*[0-9]+\\.[0-9]+$', '', V1))
    )
][, symbol := alias2SymbolTable(symbol)][
    symbol %in% names(table(symbol))[table(symbol) > 1],
    symbol := mapIds(
        org.Hs.eg.db,
        stringr::str_split_fixed(ensembl_id, '\\.', 2)[, 1], # Don't need version number
        keytype = 'ENSEMBL',
        column = 'SYMBOL'
    )
][!is.na(symbol) & symbol %in% names(table(symbol))[table(symbol) == 1]][, c('V1', 'ensembl_id') := NULL]

setcolorder(sc_data, 'symbol')

# Filter out cells with fewer than 1000 genes detected (there aren't many):
sc_data <- sc_data[, c(TRUE, apply(sc_data[, -'symbol'], 2, function(x) sum(x > 0)) >= 1000), with = FALSE]

# Convert to TPM (the values are already in FPKM, so conversion is the same process as for length-normalised counts, or rates, and they're not on a
# log scale so we don't have to take (2^matrix - 1) before doing this):
sc_data[, names(sc_data[, -'symbol']) := lapply(.SD, function(x) {round(log2(1e+05*x/sum(x) + 1), 4)}), .SDcols = -'symbol']

sc_data <- tdt(sc_data)

# Extract cell type names from the cell IDs:
sc_data[
    ,
    cell_type := mapvalues(tolower(str_split_fixed(id, '__', 3)[, 2]), c('tcell', 'bcell', 'na', 'mastcell'), c('t_cell', 'b_cell', NA, 'mast'))
][, id := str_split_fixed(id, '__', 3)[, 1]]

# Get series matrix (can't seem to download it from GEO now, but I assume the one I downloaded previously is the one on GEO) - note I'm pasting .I to
# V1 because some of them are the same, which is not good for variabe names:

series_mat <- tdt(
    fread('../../single_cell_data/li_colorectal_2017/series_matrix.txt', skip = 38, nrows = 44, header = FALSE)[
        ,
        V1 := paste0(gsub('!', '', V1), paste0('_', .I))
    ]
)[, -'id']

setkey(series_mat, Sample_title_1)

sc_data[
    ,
    c('patient', 'sample_type') := series_mat[
        id,
        .(
            as.numeric(gsub('patient id: CRC', '', Sample_characteristics_ch1_10)),
            mapvalues(Sample_characteristics_ch1_12, c('tissue type: Colorectal Tumor', 'tissue type: Normal Mucosa'), c('tumour', 'normal'))
        )
    ]
]

# Get epithelial cluster assignments:
epi_clusters <- unlist(
    lapply(
        c('NM', 'tumor'),
        function(x) {
            as.character(
                fread(paste0('../../single_cell_data/li_colorectal_2017/GSE81861_CRC_', x, '_epithelial_cells_FPKM.csv'), nrows = 1, header = FALSE)[
                    ,
                    -'V1'
                ]
            )
        }
    )
)
epi_clusters <- setNames(as.data.table(str_split_fixed(epi_clusters, '__', 3)[, 1:2]), c('cell_id', 'epi_cluster'))
setkey(epi_clusters, cell_id)
sc_data[cell_type == 'epithelial', epithelial_cluster := epi_clusters[id, epi_cluster]]

setcolorder(sc_data, c('id', 'patient', 'cell_type', 'sample_type', 'epithelial_cluster'))

fwrite(sc_data, '../data_and_figures/li_colorectal_2017.csv')





# CRC data from Lee et al., 2020:

# This data is divided into two parts: 23 Korean patients from the Samsung Medical Center (SMC) and 6 Belgian patients from Katholieke Universiteit
# Leuven (KUL3).  It's probably best to treat the SMC and KUL3 datasets separately, as they do in the paper, to avoid batch effects.

# To download the count matrices in R:

# download.file(
	# 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE132465&format=file&file=GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt.gz',
	# destfile = 'GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt.gz'
# )

# R.utils::gunzip('GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt.gz')

# download.file(
	# 'https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE144735&format=file&file=GSE144735_processed_KUL3_CRC_10X_raw_UMI_count_matrix.txt.gz',
	# destfile = 'GSE144735_processed_KUL3_CRC_10X_raw_UMI_count_matrix.txt.gz'
# )

# R.utils::gunzip('GSE144735_processed_KUL3_CRC_10X_raw_UMI_count_matrix.txt.gz')

# I first ran the following (in the lee_crc_2020 folder) to convert to sparse matrices:

# sc_data_smc <- fread('../../single_cell_data/lee_crc_2020/GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.txt')
# sc_data_smc <- set_rownames(as.matrix(sc_data_smc[, -'Index']), sc_data_smc$Index)
# sc_data_smc <- as(sc_data_smc, 'dgCMatrix')
# saveRDS(sc_data_smc, 'GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.rds')

# Same for the KUL3 dataset:

# sc_data_kul3 <- fread('../../single_cell_data/lee_crc_2020/GSE144735_processed_KUL3_CRC_10X_raw_UMI_count_matrix.txt')
# sc_data_kul3 <- set_rownames(as.matrix(sc_data_kul3[, -'Index']), sc_data_kul3$Index)
# sc_data_kul3 <- as(sc_data_kul3, 'dgCMatrix')
# saveRDS(sc_data_kul3, 'GSE144735_processed_KUL3_CRC_10X_raw_UMI_count_matrix.rds')

# From sparse matrix RDS files:



# First, SMC:

sc_data <- readRDS('../../single_cell_data/lee_crc_2020/GSE132465_GEO_processed_CRC_10X_raw_UMI_count_matrix.rds')
sc_meta <- fread('../../single_cell_data/lee_crc_2020/GSE132465_GEO_processed_CRC_10X_cell_annotation.txt', key = 'Index')

# Update gene names:
gn <- alias2SymbolTable(rownames(sc_data))
rownames(sc_data) <- gn
sc_data <- sc_data[!is.na(gn) & gn %in% names(table(gn))[table(gn) == 1], ]

# Remove cells with fewer than 1000 genes detected:
sc_data <- sc_data[, col_nnz(sc_data) >= 1000]

sc_data <- t(sc_data)

# Convert to log TPM/10 and round to 4 decimal places:
sc_data <- round(log_transform(to_frac(sc_data, MARGIN = 1)*1e+05), 4)

# Add cell types and patient IDs, converting to data.table:

sc_meta <- sc_meta[rownames(sc_data)]

sc_data <- cbind(
	sc_meta[
		,
		.(
			id = Index,
			patient = Patient,
			sample_id = Sample,
			sample_type = mapvalues(Class, c('Normal', 'Tumor'), c('normal', 'tumour')),
			cell_type = mapvalues(
				Cell_type,
				c('B cells', 'Epithelial cells', 'Mast cells', 'Myeloids', 'Stromal cells', 'T cells'),
				c('b_cell', 'epithelial', 'mast', 'myeloid', 'fibroblast', 't_cell')
			),
			cell_subtype = Cell_subtype
		)
	],
	as.data.table(as.matrix(sc_data))
)

fwrite(sc_data, '../data_and_figures/lee_crc_2020_smc.csv')



# Now KUL3:

sc_data <- readRDS('../../single_cell_data/lee_crc_2020/GSE144735_processed_KUL3_CRC_10X_raw_UMI_count_matrix.rds')
sc_meta <- fread('../../single_cell_data/lee_crc_2020/GSE144735_processed_KUL3_CRC_10X_annotation.txt', key = 'Index')

# Update gene names:
gn <- alias2SymbolTable(rownames(sc_data))
rownames(sc_data) <- gn
sc_data <- sc_data[!is.na(gn) & gn %in% names(table(gn))[table(gn) == 1], ]

# Remove cells with fewer than 1000 genes detected:
sc_data <- sc_data[, col_nnz(sc_data) >= 1000]

sc_data <- t(sc_data)

# Convert to log TPM/10 and round to 4 decimal places:
sc_data <- round(log_transform(to_frac(sc_data, MARGIN = 1)*1e+05), 4)

# Add cell types and patient IDs, converting to data.table:

sc_meta <- sc_meta[rownames(sc_data)]

sc_data <- cbind(
	sc_meta[
		,
		.(
			id = Index,
			patient = Patient,
			sample_id = Sample,
			sample_type = mapvalues(Class, c('Normal', 'Border', 'Tumor'), c('normal', 'border', 'tumour')),
			cell_type = mapvalues(
				Cell_type,
				c('B cells', 'Epithelial cells', 'Mast cells', 'Myeloids', 'Stromal cells', 'T cells'),
				c('b_cell', 'epithelial', 'mast', 'myeloid', 'fibroblast', 't_cell')
			),
			cell_subtype = Cell_subtype
		)
	],
	as.data.table(as.matrix(sc_data))
)

fwrite(sc_data, '../data_and_figures/lee_crc_2020_kul3.csv')





# HNSCC data from Puram et al.:

sc_data <- fread('../../single_cell_data/puram_hnscc_2017/HNSCC_all_data.txt', skip = 6, header = FALSE)[, V1 := gsub("'", "", V1)]
sc_meta <- tdt(fread('../../single_cell_data/puram_hnscc_2017/HNSCC_all_data.txt', nrows = 5), new_id = 'cell_id')
setnames(sc_data, names(sc_data), c('id', sc_meta$cell_id))

# Convert gene IDs using alias2SymbolTable() and delete repeated genes:
sc_data[, id := limma::alias2SymbolTable(id)]
sc_data <- sc_data[!is.na(id) & id %in% names(table(id))[table(id) == 1]]

# Re-scale to log TPM/10 (need to reverse the log):
sc_data[, names(sc_data[, -'id']) := lapply(.SD, function(x) {round(log2(1e+05*(2^x - 1)/sum(2^x - 1) + 1), 4)}), .SDcols = -'id']

sc_data <- tdt(sc_data)

# There are no cells with fewer than 1000 genes detected, so no need to filter for this.

# Use tumour assignment file to save me decomposing the cell IDs to get patient numbers:
tumour_assignments <- fread('../../single_cell_data/puram_hnscc_2017/puram_2017_tumour_assignment.csv', key = 'cell')
sc_meta[, cell_type := switch((`classified  as cancer cell` == 1) + 1, `non-cancer cell type`, 'cancer'), by = cell_id]
sc_meta[, c('patient', 'cell_type') := .(tumour_assignments[cell_id, tumor], mapvalues(tolower(gsub('-', '', gsub(' ', '_', cell_type))), '0', NA))]

setkey(sc_meta, cell_id)

sc_data[
    ,
    c('patient', 'cell_type', 'lymph_node', 'processed_by_maxima_enzyme') := sc_meta[
        id,
        .(patient, cell_type, `Lymph node`, `processed by Maxima enzyme`)
    ]
]

setcolorder(sc_data, c('id', 'patient', 'cell_type', 'lymph_node', 'processed_by_maxima_enzyme'))

# Subset only the patients whose data were deemed high-quality in the Puram et al. paper:
sc_data <- sc_data[patient %in% c(5, 6, 16, 17, 18, 20, 22, 25, 26, 28)]

fwrite(sc_data, '../data_and_figures/puram_hnscc_2017.csv')





# PDAC dataset from Elyada et al.:

sc_data <- fread('../../single_cell_data/elyada_pdac_2019/elyada_pdac_2019.txt', key = 'V1')

# The data has been transformed by natural log: see the Series Matrix File from the GEO page for the mouse scRNA-seq data from the same study
# (https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE129455), which says in one of the rows, "Each filtered counts matrix was then normalized such
# that the number of UMIs in each cell is equal to the median UMI count across the dataset and natural log transformed."

cluster_data <- fread('../../single_cell_data/elyada_pdac_2019/combined_tsne.csv', key = 'V1')
names(cluster_data)[1] <- 'cell_id'
genes_data <- fread('../../single_cell_data/elyada_pdac_2019/genes.tsv', header = FALSE, key = 'V1')[V1 %in% sc_data$V1]

# Extract patient IDs from the cell IDs, according to info from Ela:

# Each barcode has a hyphen and then a number (1-10). The numbers refer to the following samples:

# "1" = "HT97",
# "2" = "HT99",
# "3" = "HT103",
# "4" = "HT137",
# "5" = "HT137TN",
# "6" = "HT143",
# "7" = "HN149",
# "8" = "HT149",
# "9" = "HN150",
# "10" = "HT150"

# However, the matrix only contains 8 samples and not 10. 6 was omitted from the dataset because it's not a PDAC and 5 is the fibroblast-enriched
# sample, which was analyzed in Fig. 3. Samples 7 and 9 are the adj. normals relating to 8 and 10, respectively, and are present in the matrix.

# The following gives a warning, "The following `from` values were not present in `x`: 5, 6", which agrees nicely with what Ela said.

cluster_data[
    ,
    patient_id := mapvalues(
        str_split_fixed(cell_id, '-', 2)[, 2],
        1:10,
        c('HT97', 'HT99', 'HT103', 'HT137', 'HT137TN', 'HT143', 'HN149', 'HT149', 'HN150', 'HT150')
    )
]

sc_data <- sc_data[!(V1 %in% c('ENSGGENES', 'ENSGUMI'))]

# Convert the ensembl gene IDs to symbols:
genes_data[, ann_dbi_symbols := AnnotationDbi::mapIds(org.Hs.eg.db, keys = V1, keytype = 'ENSEMBL', column = 'SYMBOL')]
genes_data <- genes_data[!is.na(ann_dbi_symbols) & !(ann_dbi_symbols %in% names(table(ann_dbi_symbols))[table(ann_dbi_symbols) > 1])]

sc_data[, V1 := genes_data[sc_data$V1, ann_dbi_symbols]]
sc_data <- sc_data[!is.na(V1)]

# Filter out cells with fewer than 1000 genes detected:
sc_data <- sc_data[, c(TRUE, apply(sc_data[, -'V1'], 2, function(x) sum(x > 0)) >= 1000), with = FALSE]

# Re-scale the data by exponentiating, converting to TPM/10, then taking log to the base 2:
sc_data[, names(sc_data[, -'V1']) := lapply(.SD, function(x) {round(log2(1e+05*(exp(x) - 1)/sum(exp(x) - 1) + 1), 4)}), .SDcols = -'V1']

sc_data <- tdt(sc_data)

# Add column of cluster assignments, changing a few of the cluster names:
sc_data[
    ,
    c('patient', 'cluster', 'cell_type') := cluster_data[
        sc_data$id,
        .(
            patient_id,
            cluster_name,
            mapvalues(
                cluster_name,
                c('Acinar cells', 'Activated DC', 'B cells', 'Ductal cells 1', 'Ductal cells 2', 'Ductal cells 3', 'Fenestrated EC', 'Fibroblasts',
                  'Mast cells', 'Myeloid & Macrophages', 'pDCs', 'Perivascular cells', 'Plasma cells', 'RBCs', 'T and NK cells'),
                c('acinar', 'dendritic', 'b_cell', 'ductal', 'ductal', 'ductal', 'endothelial', 'fibroblast', 'mast', 'myeloid', 'dendritic',
                  'perivascular', 'plasma', 'red_blood_cell', 'lymphoid')
            )
        )
    ]
]

setcolorder(sc_data, c('id', 'patient', 'cell_type', 'cluster'))

fwrite(sc_data, '../data_and_figures/elyada_pdac_2019.csv')





# PDAC data from Peng et al., 2019:

# Use the following to download the files in R:

# To see what the files are at the FTP address:
# RCurl::getURL('ftp://download.big.ac.cn/gsa/CRA001160/other_files', ftp.use.epsv = FALSE, dirlistonly = TRUE)

# To actually download them (from within the peng_pdac_2019 directory in single_cell_data):
# download.file('ftp://download.big.ac.cn/gsa/CRA001160/other_files/all_celltype.txt', destfile = 'all_celltype.txt')
# download.file('ftp://download.big.ac.cn/gsa/CRA001160/other_files/count-matrix.txt', destfile = 'count-matrix.txt')

# I then ran the following to convert the data matrix to a sparse matrix:
# sc_data <- fread('../../single_cell_data/peng_pdac_2019/count-matrix.txt')
# sc_data <- set_rownames(as.matrix(sc_data[, -'V1']), sc_data$V1)
# sc_data <- as(sc_data, 'dgCMatrix')
# saveRDS(sc_data, 'count-matrix.rds')

# Preprocessing from sparse matrix:

sc_data <- readRDS('../../single_cell_data/peng_pdac_2019/count-matrix.rds')
cell_types_data <- fread('../../single_cell_data/peng_pdac_2019/all_celltype.txt', key = 'cell.name')

# When we run the following, we see that two of the repeated gene names are AREG and TIMP2:
# updated_gene_names <- alias2SymbolTable(rownames(sc_data))
# names(table(updated_gene_names))[table(updated_gene_names) > 1]
# Since these are potentially important pEMT genes, I'll remove AREGB and DDC8 from the data, which are the two extra genes that map to AREG and TIMP2
# (respectively), as can be seen from the following:
# rownames(sc_data)[!is.na(updated_gene_names) & updated_gene_names == 'AREG']
# rownames(sc_data)[!is.na(updated_gene_names) & updated_gene_names == 'TIMP2']

sc_data <- sc_data[!(rownames(sc_data) %in% c('AREGB', 'DDC8')), ]
gn <- alias2SymbolTable(rownames(sc_data))
rownames(sc_data) <- gn
sc_data <- sc_data[!is.na(gn) & gn %in% names(table(gn))[table(gn) == 1], ]

# Remove cells with fewer than 1000 genes detected (we still retain a lot of cells...):
sc_data <- sc_data[, col_nnz(sc_data) >= 1000]

sc_data <- t(sc_data)

# Convert to log TPM/10 and round to 4 decimal places:
sc_data <- round(log_transform(to_frac(sc_data, MARGIN = 1)*1e+05), 4)

# Add cell types and patient IDs, converting to data.table:

cell_types_data <- cell_types_data[rownames(sc_data)]

cell_types_data[
	,
	c('cluster', 'patient', 'sample_type') := .(
		mapvalues(
			cluster,
			c('Acinar cell', 'B cell', 'Ductal cell type 1', 'Ductal cell type 2', 'Endocrine cell', 'Endothelial cell', 'Fibroblast cell',
              'Macrophage cell', 'Stellate cell', 'T cell'),
			c('acinar', 'b_cell', 'ductal_1', 'ductal_2', 'endocrine', 'endothelial', 'fibroblast', 'macrophage', 'stellate', 't_cell')
		),
		str_split_fixed(cell.name, '_', 2)[, 1],
		mapvalues(str_extract(cell.name, '^[A-Z]'), c('T', 'N'), c('tumour', 'normal'))
	)
]

setnames(cell_types_data, c('cell.name', 'cluster'), c('id', 'cell_type'))

sc_data <- cbind(cell_types_data, as.data.table(as.matrix(sc_data)))

fwrite(sc_data, '../data_and_figures/peng_pdac_2019.csv')





# Liver dataset from Ma et al. (2019):

sc_data <- lapply(
    1:2,
    function(i) {

        file_root <- paste0('../../single_cell_data/ma_liver_2019/GSE125449_Set', i, '_')

        barcodes <- fread(paste0(file_root, 'barcodes.tsv'), header = FALSE)$V1
        genes <- setNames(fread(paste0(file_root, 'genes.tsv'), header = FALSE), c('ensembl_id', 'symbol'))

        # Get gene symbols using AnnotationDbi on the Ensembl IDs.  I could then remove all genes where the symbol and ensembl_id columns disagree,
        # but it seems like doing that will get rid of some important genes, like CAVIN1.
        genes[, ann_dbi_symbols := AnnotationDbi::mapIds(org.Hs.eg.db, keys = ensembl_id, keytype = 'ENSEMBL', column = 'SYMBOL')]

        samples <- fread(paste0(file_root, 'samples.txt'), key = 'Cell Barcode')

        # Get count matrix and subset genes which are not NA or repeats:
        expression_matrix <- Matrix::readMM(paste0(file_root, 'matrix.mtx'))
        expression_matrix <- expression_matrix[
            genes[, !is.na(ann_dbi_symbols) & !(ann_dbi_symbols %in% names(table(ann_dbi_symbols))[table(ann_dbi_symbols) > 1])],
        ]

        expression_matrix <- as.matrix(t(expression_matrix))
        rownames(expression_matrix) <- barcodes
        colnames(expression_matrix) <- genes[
            !is.na(ann_dbi_symbols) & !(ann_dbi_symbols %in% names(table(ann_dbi_symbols))[table(ann_dbi_symbols) > 1]),
            ann_dbi_symbols
        ]

        # Convert to data table and add metadata columns:
        expression_matrix <- as.data.table(expression_matrix, keep.rownames = 'id')
        expression_matrix[, c('sample', 'cell_type') := samples[id, .(Sample, Type)]]
        setcolorder(expression_matrix, c('id', 'cell_type', 'sample'))

        expression_matrix

    }
)

sc_data <- rbindlist(sc_data, use.names = TRUE, fill = TRUE)
sc_data[is.na(sc_data)] <- 0

# Change cell types and add patient ID and disease type columns (patient IDs as they appear in the paper):
sc_data[
    ,
    c('cell_type', 'patient') := .(
        mapvalues(
            cell_type,
            c('B cell', 'CAF', 'HPC-like', 'Malignant cell', 'T cell', 'TAM', 'TEC'),
            c('b_cell', 'fibroblast', 'hpc-like', 'cancer', 't_cell', 'macrophage', 'endothelial')
        ),
        as.numeric(stringr::str_extract(sample, '[0-9]+$'))
    )
][
    ,
    disease := switch((patient %in% c(18, 21, 23, 28, 30, 34, 37, 38, 65)) + 1, 'iCCA', 'HCC'),
    by = patient
][, patient := paste0(mapvalues(disease, c('iCCA', 'HCC'), c('C', 'H')), patient)]

setcolorder(sc_data, c('id', 'patient', 'cell_type', 'sample', 'disease'))

# Filter out cells with fewer than 1000 genes detected:
sc_data <- sc_data[apply(sc_data[, -c('id', 'patient', 'cell_type', 'sample', 'disease')], 1, function(x) sum(x > 0)) >= 1000]

# Convert to log TPM/10:
sc_data <- cbind(
    sc_data[, .(id, patient, cell_type, sample, disease)],
    apply(sc_data[, -c('id', 'patient', 'cell_type', 'sample', 'disease')], 1, function(x) round(log2(1e+05*x/sum(x) + 1), 4)) %>% t
)

fwrite(sc_data, '../data_and_figures/ma_liver_2019.csv')





# Ovarian SS2 data from Izar et al., 2020:

# Had to first manually delete an extra tab at the end of line 4 (patient numbers) in Notepad++...

sc_data <- fread('../../single_cell_data/izar_ovarian_2020/Izar_HGSOC_ascites_SS2_log.tsv')

sc_meta <- tdt(sc_data[1:5], new_id = 'Cell_ID')
setkey(sc_meta, 'Cell_ID')

sc_data <- sc_data[-(1:5)]

# Update gene names:
gn <- alias2SymbolTable(sc_data$Cell_ID)
sc_data$Cell_ID <- gn
sc_data <- sc_data[!is.na(gn) & gn %in% names(table(gn))[table(gn) == 1], ]

# There are no cells with fewer than 1000 genes detected.  In fact, the lowest number of genes detected among these cells is over 3000.

# Rescale to log TPM/10 and round to 4 decimal places:
sc_data[, names(sc_data[, -'Cell_ID']) := lapply(.SD, function(x) round(log2(1e+05*(2^x - 1)/sum(2^x - 1) + 1), 4)), .SDcols = -'Cell_ID']

sc_data <- tdt(sc_data)

# Add metadata:
sc_data <- cbind(
	sc_data[, .(id)],
	sc_meta[sc_data$id, .(patient = Patient, cell_type = mapvalues(clst, c(1:6, 7), c(rep('cancer', 6), 'fibroblast')), cluster = clst)],
	sc_data[, -'id']
)

fwrite(sc_data, '../data_and_figures/izar_ovarian_2020_ss2.csv')





# Ovarian 10x data from Izar et al., 2020:

sc_data <- fread('../../single_cell_data/izar_ovarian_2020/Izar_HGSOC_ascites_10x_log.tsv')

sc_meta <- tdt(sc_data[1:7], new_id = 'Cell_ID')
setkey(sc_meta, 'Cell_ID')

sc_data <- sc_data[-(1:7)]

# Columns are character vectors, because of the 10x ID row, so change this:

sc_data[, names(sc_data[, -'Cell_ID']) := lapply(.SD, as.numeric), .SDcols = -'Cell_ID']

# Update gene names:
gn <- alias2SymbolTable(sc_data$Cell_ID)
sc_data$Cell_ID <- gn
sc_data <- sc_data[!is.na(gn) & gn %in% names(table(gn))[table(gn) == 1]]

# Take cells with at least 1000 genes detected:

sc_data <- sc_data[, c(TRUE, apply(sc_data[, -'Cell_ID'], 2, function(x) {sum(x > 0) >= 1000})), with = FALSE]

# Rescale to log TPM/10 and round to 4 decimal places:
sc_data[, names(sc_data[, -'Cell_ID']) := lapply(.SD, function(x) round(log2(1e+05*(2^x - 1)/sum(2^x - 1) + 1), 4)), .SDcols = -'Cell_ID']

sc_data <- tdt(sc_data)

# Add metadata:
sc_data <- cbind(
	sc_data[, .(id)],
	sc_meta[
		sc_data$id,
		.(
			barcode = `10x_barcode`,
			patient = patient,
			sample_id = sample_ID,
			time = time,
			cell_type = mapvalues(clst, 1:9, c(rep('cancer', 5), rep('fibroblast', 4))),
			cluster = clst
		)
	],
	sc_data[, -'id']
)

fwrite(sc_data, '../data_and_figures/izar_ovarian_2020_10x.csv')





# 4 datasets from Qian et al., 2020 - note the lung one is the same as the one from Lambrechts et al., 2018, but includes all 8 samples.



# Breast:

sc_data <- readMM('../../single_cell_data/qian_2020/BC_counts/matrix.mtx')
barcodes <- fread('../../single_cell_data/qian_2020/BC_counts/barcodes.tsv', header = FALSE)$V1

# The genes table has 2 columns, but they are identical
genes <- fread('../../single_cell_data/qian_2020/BC_counts/genes.tsv', header = FALSE)$V1

# Update gene names:
gn <- alias2SymbolTable(genes)
sc_data <- sc_data[!is.na(gn) & gn %in% names(table(gn))[table(gn) == 1], ]
genes <- gn[!is.na(gn) & gn %in% names(table(gn))[table(gn) == 1]]

# Take cells with at least 1000 genes detected:
genes_detected <- col_nnz(sc_data)
barcodes <- barcodes[genes_detected >= 1000]
sc_data <- sc_data[, genes_detected >= 1000]

# Add row and column names:
rownames(sc_data) <- genes
colnames(sc_data) <- barcodes

# Rescale to log TPM/10, round to 4 decimal places and transpose:
sc_data <- round(log_transform(1e+05*to_frac(sc_data)), 4)
sc_data <- t(sc_data)

# Add metadata:
sc_meta <- fread('../../single_cell_data/qian_2020/BC_metadata.csv', key = 'Cell')[
	rownames(sc_data),
	.(id = Cell, patient = PatientNumber, cell_type = CellType)
]

sc_data <- cbind(sc_meta, as.data.table(as.matrix(sc_data), keep.rownames = FALSE))

fwrite(sc_data, '../data_and_figures/qian_breast_2020.csv')



# CRC:

sc_data <- readMM('../../single_cell_data/qian_2020/CRC_counts/matrix.mtx')
barcodes <- fread('../../single_cell_data/qian_2020/CRC_counts/barcodes.tsv', header = FALSE)$V1

# The genes table has 2 columns, but they are identical
genes <- fread('../../single_cell_data/qian_2020/CRC_counts/genes.tsv', header = FALSE)$V1

# Update gene names:
gn <- alias2SymbolTable(genes)
sc_data <- sc_data[!is.na(gn) & gn %in% names(table(gn))[table(gn) == 1], ]
genes <- gn[!is.na(gn) & gn %in% names(table(gn))[table(gn) == 1]]

# Take cells with at least 1000 genes detected:
genes_detected <- col_nnz(sc_data)
barcodes <- barcodes[genes_detected >= 1000]
sc_data <- sc_data[, genes_detected >= 1000]

# Add row and column names:
rownames(sc_data) <- genes
colnames(sc_data) <- barcodes

# Rescale to log TPM/10, round to 4 decimal places and transpose:
sc_data <- round(log_transform(1e+05*to_frac(sc_data)), 4)
sc_data <- t(sc_data)

# Add metadata:
sc_meta <- fread('../../single_cell_data/qian_2020/CRC_metadata.csv', key = 'Cell')[
	rownames(sc_data),
	.(
		id = Cell,
		patient = PatientNumber,
		sample_type = mapvalues(TumorSite, c('C', 'B', 'N'), c('core', 'border', 'normal')),
		cell_type = CellType
	)
]

sc_data <- cbind(sc_meta, as.data.table(as.matrix(sc_data), keep.rownames = FALSE))

fwrite(sc_data, '../data_and_figures/qian_crc_2020.csv')



# Lung:

sc_data <- readMM('../../single_cell_data/qian_2020/LC_counts/matrix.mtx')
barcodes <- fread('../../single_cell_data/qian_2020/LC_counts/barcodes.tsv', header = FALSE)$V1

# The genes table has 2 columns, but they are identical
genes <- fread('../../single_cell_data/qian_2020/LC_counts/genes.tsv', header = FALSE)$V1

# Update gene names:
gn <- alias2SymbolTable(genes)
sc_data <- sc_data[!is.na(gn) & gn %in% names(table(gn))[table(gn) == 1], ]
genes <- gn[!is.na(gn) & gn %in% names(table(gn))[table(gn) == 1]]

# Take cells with at least 1000 genes detected:
genes_detected <- col_nnz(sc_data)
barcodes <- barcodes[genes_detected >= 1000]
sc_data <- sc_data[, genes_detected >= 1000]

# Add row and column names:
rownames(sc_data) <- genes
colnames(sc_data) <- barcodes

# Rescale to log TPM/10, round to 4 decimal places and transpose:
sc_data <- round(log_transform(1e+05*to_frac(sc_data)), 4)
sc_data <- t(sc_data)

# Add metadata:

# The designation of 'core', 'middle' and 'edge' is missing from some of these tumours: for 6 and 7, TumorSite takes only the values "N" and "U". I'm
# assuming from the numbers (comparing them with supplementary table 2 from Qian et al.) that these are the normal and tumour cells, respectively.
# But I don't know why they haven't given the more precise tumour site information, when the supplementary table (and files in ArrayExpress) suggest
# they did take samples from core, middle and edge.

sc_meta <- fread('../../single_cell_data/qian_2020/LC_metadata.csv', key = 'Cell')[
	rownames(sc_data),
	.(
		id = Cell,
		patient = PatientNumber,
		sample_type = mapvalues(TumorSite, c('I', 'M', 'O', 'N', 'U'), c('core', 'middle', 'border', 'normal', 'tumour')),
		disease = mapvalues(PatientNumber, 1:8, c('LUSC', 'LUSC', 'LUAD', 'LUAD', 'Large cell', 'LUAD', 'LUSC', 'Pleimorphic')),
		cell_type = CellType
	)
]

sc_data <- cbind(sc_meta, as.data.table(as.matrix(sc_data), keep.rownames = FALSE))

fwrite(sc_data, '../data_and_figures/qian_lung_2020.csv')



# Ovarian:

sc_data <- readMM('../../single_cell_data/qian_2020/OvC_counts/matrix.mtx')
barcodes <- fread('../../single_cell_data/qian_2020/OvC_counts/barcodes.tsv', header = FALSE)$V1

# The genes table has 2 columns, but they are identical
genes <- fread('../../single_cell_data/qian_2020/OvC_counts/genes.tsv', header = FALSE)$V1

# Update gene names:
gn <- alias2SymbolTable(genes)
sc_data <- sc_data[!is.na(gn) & gn %in% names(table(gn))[table(gn) == 1], ]
genes <- gn[!is.na(gn) & gn %in% names(table(gn))[table(gn) == 1]]

# Take cells with at least 1000 genes detected:
genes_detected <- col_nnz(sc_data)
barcodes <- barcodes[genes_detected >= 1000]
sc_data <- sc_data[, genes_detected >= 1000]

# Add row and column names:
rownames(sc_data) <- genes
colnames(sc_data) <- barcodes

# Rescale to log TPM/10, round to 4 decimal places and transpose:
sc_data <- round(log_transform(1e+05*to_frac(sc_data)), 4)
sc_data <- t(sc_data)

# Add metadata:

# There is a mismatch between the metadata and supplementary table 2 from the paper - the metadata says patient 11 has one tumour and one normal
# sample from each of omentum and peritoneum; the supplementary table says both normals are from the omentum, while there are two tumour samples from
# the peritoneum.  I'll assume the metadata is correct, because it's less likely to have been affected by e.g. copy-paste issues, and anyway it would
# make much more sense to take tumour and normal from each of the two sites.

sc_meta <- fread('../../single_cell_data/qian_2020/OvC_metadata.csv', key = 'Cell')[
	rownames(sc_data),
	.(
		id = Cell,
		patient = PatientNumber,
		sample_type = mapvalues(CellFromTumor, 0:1, c('normal', 'tumour')),
		site = TumorSite,
		cell_type = CellType
	)
]

sc_data <- cbind(sc_meta, as.data.table(as.matrix(sc_data), keep.rownames = FALSE))

fwrite(sc_data, '../data_and_figures/qian_ovarian_2020.csv')





# LUAD data from Kim et al.:

# sc_data1 <- fread('../../single_cell_data/kim_luad_2020/GSE131907_Lung_Cancer_raw_UMI_matrix.txt', select = 1:50000)
# sc_data1 <- as(set_rownames(as.matrix(sc_data1[, -'Index']), sc_data1$Index), 'dgCMatrix')
# sc_data1 <- sc_data1[, col_nnz(sc_data1) >= 1000]
# saveRDS(sc_data1, '../../single_cell_data/kim_luad_2020/matrix1.rds')

# sc_data2 <- fread('../../single_cell_data/kim_luad_2020/GSE131907_Lung_Cancer_raw_UMI_matrix.txt', select = 50001:100000)
# sc_data2 <- as(as.matrix(sc_data2), 'dgCMatrix')
# sc_data2 <- sc_data2[, col_nnz(sc_data2) >= 1000]
# saveRDS(sc_data2, '../../single_cell_data/kim_luad_2020/matrix2.rds')

# sc_data3 <- fread('../../single_cell_data/kim_luad_2020/GSE131907_Lung_Cancer_raw_UMI_matrix.txt', select = 100001:150000)
# sc_data3 <- as(as.matrix(sc_data3), 'dgCMatrix')
# sc_data3 <- sc_data3[, col_nnz(sc_data3) >= 1000]
# saveRDS(sc_data3, '../../single_cell_data/kim_luad_2020/matrix3.rds')

# sc_data4 <- fread('../../single_cell_data/kim_luad_2020/GSE131907_Lung_Cancer_raw_UMI_matrix.txt', select = 150001:208507)
# sc_data4 <- as(as.matrix(sc_data4), 'dgCMatrix')
# sc_data4 <- sc_data4[, col_nnz(sc_data4) >= 1000]
# saveRDS(sc_data4, '../../single_cell_data/kim_luad_2020/matrix4.rds')

sc_data1 <- readRDS('../../single_cell_data/kim_luad_2020/matrix1.rds')
sc_data2 <- readRDS('../../single_cell_data/kim_luad_2020/matrix2.rds')
sc_data3 <- readRDS('../../single_cell_data/kim_luad_2020/matrix3.rds')
sc_data4 <- readRDS('../../single_cell_data/kim_luad_2020/matrix4.rds')

sc_data <- cbind(sc_data1, sc_data2, sc_data3, sc_data4)

rm(sc_data1)
rm(sc_data2)
rm(sc_data3)
rm(sc_data4)

gn <- alias2SymbolTable(rownames(sc_data))
sc_data <- sc_data[!is.na(gn) & gn %in% names(table(gn))[table(gn) == 1], ]
rownames(sc_data) <- gn[!is.na(gn) & gn %in% names(table(gn))[table(gn) == 1]]

saveRDS(sc_data, '../../single_cell_data/kim_luad_2020/kim_luad_2020_matrix_counts.rds')

sc_data <- round(log_transform(1e+05*to_frac(sc_data)), 4)

saveRDS(sc_data, '../../single_cell_data/kim_luad_2020/kim_luad_2020_matrix_scaled.rds')

sc_data <- t(sc_data)

# Add metadata:

sc_meta <- fread('../../single_cell_data/kim_luad_2020/GSE131907_Lung_Cancer_cell_annotation.txt', key = 'Index')

sample_origin_data <- as.data.table(readxl::read_xlsx('../../single_cell_data/kim_luad_2020/41467_2020_16164_MOESM4_ESM.xlsx', skip = 2, n_max = 58))
setkey(sample_origin_data, Samples)

sc_meta <- sc_meta[
	rownames(sc_data),
	.(
		id = Index,
		patient = sample_origin_data[Sample, `Patient id`],
		sample_type = mapvalues(
			str_extract(Sample_Origin, '^[a-z]|^[A-Z]'),
			c('m', 'n', 'P', 't'),
			c('metastasis', 'normal', 'pleural effusion', 'tumour')
		),
		sample_site = mapvalues(
			gsub('^[a-z]', '', Sample_Origin),
			c('Brain', 'L/B', 'LN', 'Lung', 'PE'),
			c('brain', 'lung', 'lymph node', 'lung', 'pleural fluids')
		),
		sample_id = Sample,
		sample_origin = Sample_Origin,
		cell_type = Cell_type,
		cell_type_refined = Cell_type.refined,
		cell_subtype = Cell_subtype
	)
]

fwrite(sc_meta, '../../single_cell_data/kim_luad_2020/kim_luad_2020_metadata.csv')

# Now I need to save the expression data in parts, to save memory.  I will split it up by sample type:

fwrite(
	cbind(sc_meta[sample_type == 'tumour', -'sample_type'], as.matrix(sc_data[sc_meta[sample_type == 'tumour', id], ])),
	'../data_and_figures/kim_luad_2020.csv'
)

fwrite(
	cbind(sc_meta[sample_type == 'normal', -'sample_type'], as.matrix(sc_data[sc_meta[sample_type == 'normal', id], ])),
	'../data_and_figures/kim_luad_2020_normal.csv'
)

fwrite(
	cbind(sc_meta[sample_type == 'metastasis', -'sample_type'], as.matrix(sc_data[sc_meta[sample_type == 'metastasis', id], ])),
	'../data_and_figures/kim_luad_2020_metastasis.csv'
)

fwrite(
	cbind(sc_meta[sample_type == 'pleural effusion', -'sample_type'], as.matrix(sc_data[sc_meta[sample_type == 'pleural effusion', id], ])),
	'../data_and_figures/kim_luad_2020_pe.csv'
)
