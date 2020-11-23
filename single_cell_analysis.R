library(data.table) # 1.12.8
library(ggplot2) # This is version 3.3.1 on my laptop's WSL setup but 3.3.0 on WEXAC.
library(cowplot) # 1.0.0
library(magrittr) # 1.5
library(limma) # 3.42.2
library(plyr) # 1.8.6
library(stringr) # 1.4.0
library(caTools) # 1.18.0
library(RColorBrewer) # 1.1.2
library(scales) # 1.1.1
library(latex2exp) # 0.4.0
library(ggrepel) # 0.8.2
library(egg) # 0.4.5
library(colorspace) # 1.4.1

source('general_functions.R')
source('tcga_functions.R')
source('sc_functions.R')





expression_data <- fread('../../TCGA_data/tcga_expression_data.csv', key = 'id')
meta_data <- fread('../../TCGA_data/tcga_meta_data.csv', key = 'id')

expression_data <- expression_data[meta_data[sample_type != 'normal', id]]
meta_data <- meta_data[sample_type != 'normal']

emt_markers <- fread('../../emt_markers.csv')[, gene := alias2SymbolTable(gene)][source != 'GO', sort(unique(gene))]





# We could select which cell types to keep and which to put under 'other' using a fraction
# of the maximum of the cell type frequencies.  I decided to select them manually instead,
# because there are certain cell types I want to show in some cases, even though they're
# a bit rare.

sc_metadata <- list(
	brca = list(
		tcga_cancer_types = 'BRCA',
		read_quote = quote(
			fread('../data_and_figures/qian_breast_2020_reclassified.csv')[
				cell_type != 'ambiguous' & id != 'sc5rJUQ064_CCATGTCCATCCCATC',
				-c('cell_type_author', 'cell_type_lenient')
			]
		),
		initial_cell_types = c('b_cell', 'endothelial', 'macrophage', 'mast', 't_cell', 'caf', 'cancer'),
		rare_cell_types = 'dendritic'
	),
	brca_lenient = list(
		tcga_cancer_types = 'BRCA',
		read_quote = quote(
			fread('../data_and_figures/qian_breast_2020_reclassified.csv')[
				cell_type_lenient != 'ambiguous' & id != 'sc5rJUQ064_CCATGTCCATCCCATC',
				-c('cell_type_author', 'cell_type')
			] %>% setnames('cell_type_lenient', 'cell_type')
		),
		initial_cell_types = c('b_cell', 'endothelial', 'macrophage', 'mast', 't_cell', 'caf', 'cancer'),
		rare_cell_types = 'dendritic'
	),
	coadread = list(
		tcga_cancer_types = c('COAD', 'READ'),
		read_quote = quote(
			fread('../data_and_figures/lee_crc_2020_smc_reclassified.csv')[cell_type != 'ambiguous', -c('cell_type_author', 'cell_type_lenient')]
		),
		initial_cell_types = c('b_cell', 'endothelial', 'macrophage', 't_cell', 'caf', 'cancer'),
		rare_cell_types = c('epithelial', 'mast')
	),
	coadread_lenient = list(
		tcga_cancer_types = c('COAD', 'READ'),
		read_quote = quote(
			fread('../data_and_figures/lee_crc_2020_smc_reclassified.csv')[
				cell_type_lenient != 'ambiguous',
				-c('cell_type_author', 'cell_type')
			] %>% setnames('cell_type_lenient', 'cell_type')
		),
		initial_cell_types = c('b_cell', 'macrophage', 't_cell', 'caf', 'cancer'), # Endothelial cells become CAFs under lenient definition
		rare_cell_types = c('epithelial', 'mast')
	),
    hnsc = list(
        tcga_cancer_types = 'HNSC',
        read_quote = quote(fread('../data_and_figures/puram_hnscc_2017_reclassified.csv')[cell_type != 'ambiguous', -'cell_type_author']),
		initial_cell_types = c('endothelial', 'macrophage', 'mast', 't_cell', 'caf', 'cancer'),
		rare_cell_types = c('b_cell', 'dendritic', 'myocyte')
    ),
	lihc = list(
        tcga_cancer_types = 'LIHC',
        read_quote = quote(
			fread('../data_and_figures/ma_liver_2019_reclassified.csv')[cell_type != 'ambiguous', -c('cell_type_author', 'cell_type_lenient')]
		),
		initial_cell_types = c('endothelial', 'macrophage', 't_cell', 'caf', 'cancer'),
		rare_cell_types = c('b_cell', 'hpc-like')
    ),
	lihc_lenient = list(
        tcga_cancer_types = 'LIHC',
        read_quote = quote(
			fread('../data_and_figures/ma_liver_2019_reclassified.csv')[
				cell_type_lenient != 'ambiguous',
				-c('cell_type_author', 'cell_type')
			] %>% setnames('cell_type_lenient', 'cell_type')
		),
		initial_cell_types = c('endothelial', 'macrophage', 't_cell', 'caf', 'cancer'),
		rare_cell_types = c('b_cell', 'hpc-like')
    ),
	luad = list(
		tcga_cancer_types = 'LUAD',
		read_quote = quote(fread('../data_and_figures/kim_luad_2020_reclassified.csv')[cell_type != 'ambiguous', -'cell_type_author']),
		initial_cell_types = c('b_cell', 'endothelial', 'macrophage', 'mast', 't_cell', 'caf', 'cancer'),
		rare_cell_types = 'dendritic'
	),
	# luad = list(
		# tcga_cancer_types = 'LUAD',
		# read_quote = quote(
		# 	fread('../data_and_figures/qian_lung_2020_reclassified.csv')[
		# 		disease == 'LUAD' & cell_type != 'ambiguous',
		# 		-c('disease', 'cell_type_author', 'cell_type_lenient')
		# 	]
		# ),
		# initial_cell_types = c('b_cell', 'endothelial', 'macrophage', 't_cell', 'caf', 'cancer'),
		# rare_cell_types = c('alveolar', 'dendritic', 'epithelial', 'mast')
	# ),
	# luad_lenient = list(
		# tcga_cancer_types = 'LUAD',
		# read_quote = quote(
		# 	fread('../data_and_figures/qian_lung_2020_reclassified.csv')[
		# 		disease == 'LUAD' & cell_type_lenient != 'ambiguous',
		# 		-c('disease', 'cell_type_author', 'cell_type')
		# 	] %>% setnames('cell_type_lenient', 'cell_type')
		# ),
		# initial_cell_types = c('b_cell', 'endothelial', 'macrophage', 't_cell', 'caf', 'cancer'),
		# rare_cell_types = c('alveolar', 'dendritic', 'epithelial', 'mast')
	# ),
    lusc = list(
        tcga_cancer_types = 'LUSC',
        read_quote = quote(
			fread('../data_and_figures/qian_lung_2020_reclassified.csv')[
				disease == 'LUSC' & cell_type != 'ambiguous',
				-c('disease', 'cell_type_author', 'cell_type_lenient')
			]
		),
		initial_cell_types = c('b_cell', 'endothelial', 'macrophage', 'mast', 't_cell', 'caf', 'cancer'),
		rare_cell_types = c('dendritic', 'erythroblast')
    ),
	lusc_lenient = list(
        tcga_cancer_types = 'LUSC',
        read_quote = quote(
			fread('../data_and_figures/qian_lung_2020_reclassified.csv')[
				disease == 'LUSC' & cell_type_lenient != 'ambiguous',
				-c('disease', 'cell_type_author', 'cell_type')
			] %>% setnames('cell_type_lenient', 'cell_type')
		),
		initial_cell_types = c('b_cell', 'endothelial', 'macrophage', 'mast', 't_cell', 'caf', 'cancer'),
		rare_cell_types = c('dendritic', 'erythroblast')
    ),
	ov = list(
		tcga_cancer_types = 'OV',
		read_quote = quote(
			fread('../data_and_figures/qian_ovarian_2020_reclassified.csv')[
				cell_type != 'ambiguous' & !(id %in% c('scrSOL001_TCATTTGTCTGTCAAG', 'scrSOL004_TTGCCGTTCTCCTATA')),
				-c('cell_type_author', 'cell_type_lenient')
			]
		),
		initial_cell_types = c('b_cell', 'endothelial', 'macrophage', 't_cell', 'caf', 'cancer'),
		rare_cell_types = NULL
	),
	ov_lenient = list(
		tcga_cancer_types = 'OV',
		read_quote = quote(
			fread('../data_and_figures/qian_ovarian_2020_reclassified.csv')[
				cell_type_lenient != 'ambiguous' & !(id %in% c('scrSOL001_TCATTTGTCTGTCAAG', 'scrSOL004_TTGCCGTTCTCCTATA')),
				-c('cell_type_author', 'cell_type')
			] %>% setnames('cell_type_lenient', 'cell_type')
		),
		initial_cell_types = c('b_cell', 'endothelial', 'macrophage', 't_cell', 'caf', 'cancer'),
		rare_cell_types = NULL
	),
    paad = list(
        tcga_cancer_types = 'PAAD',
        read_quote = quote(
			fread('../data_and_figures/peng_pdac_2019_reclassified.csv')[
				cell_type != 'ambiguous' & !(id %in% c('T8_TGGTTCCTCGCATGGC', 'T17_CGTGTAACAGTACACT')),
				-'cell_type_author'
			]
		),
		initial_cell_types = c('b_cell', 'endothelial', 'macrophage', 't_cell', 'caf', 'cancer'),
		rare_cell_types = c('acinar', 'ductal_2', 'endocrine')
    )
)

# Some background info on PDAC from previous versions of this analysis:
# Acinar cells seem to express epithelial markers at much lower levels than the ductal cells (maybe this is obvious to a real
# biologist), and it is actually the ductal cells that the cancer originates from (the clue is in the name - PDAC - also see
# https://www.pancreaticcancer.org.uk/information-and-support/facts-about-pancreatic-cancer/types-of-pancreatic-cancer/,
# where they claim that less than 1% of pancreatic cancers are acinar cell cancers).  So we can assume that the acinar cells
# won't be malignant.  This is similar to alveolar cells in the lung, which makes sense because the alveolar duct is similar
# in form to the acinar duct in pancreas.





# Define gene lists:

# If already done, skip by reading from RDS:
# genes_list <- readRDS('../data_and_figures/sc_genes_list.rds')

genes_list <- sapply(
    names(sc_metadata),
    function(ct) {
		cat(paste0(ct, '\n'))
        sc_data <- eval(sc_metadata[[ct]]$read_quote)
        unfiltered <- unique(
            c(
                emt_markers[emt_markers %in% names(sc_data)],
                top_cols_by_fun_cor(
                    expression_data[meta_data[cancer_type %in% sc_metadata[[ct]]$tcga_cancer_types, id], -'id'],
                    threshold_fun = function(x) quantile(x, 0.99)
                )[id %in% names(sc_data), id]
            )
        )
		out <- list(
			unfiltered = unfiltered,
			filtered_cancer_caf = filter_for_groups(sc_data[, c('cell_type', ..unfiltered)], groups = c('caf', 'cancer'))
		)
		# if('cell_type_lenient' %in% names(sc_data)) {
			# out <- c(
				# out,
				# list(filtered_cancer_caf_lenient = filter_for_groups(sc_data[, c('cell_type_lenient', ..unfiltered)], groups = c('caf', 'cancer')))
			# )
		# }
		out
    },
    simplify = FALSE,
    USE.NAMES = TRUE
)

# Save to RDS:

saveRDS(genes_list, '../data_and_figures/sc_genes_list.rds')





# Define parameters for constructing heatmap data:

sc_cancer_caf_args <- list(
    brca = list(
        seed = 3718,
        genes_filter_fun = function(x) {log2(mean(10*(2^x - 1)) + 1) >= 4.5 | sum(x >= 7) >= length(x)/100},
        scores_filter_fun = function(x) {log2(mean(10*(2^x - 1)) + 1) >= 6 | sum(x >= 7) >= length(x)/100}
    ),
	brca_lenient = list(
        seed = 4924,
        genes_filter_fun = function(x) {log2(mean(10*(2^x - 1)) + 1) >= 4.5 | sum(x >= 7) >= length(x)/100},
        scores_filter_fun = function(x) {log2(mean(10*(2^x - 1)) + 1) >= 6 | sum(x >= 7) >= length(x)/100}
    ),
    coadread = list(
        seed = 3361,
        genes_filter_fun = function(x) {log2(mean(10*(2^x - 1)) + 1) >= 4 | sum(x >= 7) >= length(x)/100},
        scores_filter_fun = function(x) {log2(mean(10*(2^x - 1)) + 1) >= 5.5 | sum(x >= 7) >= length(x)/100}
    ),
	coadread_lenient = list(
        seed = 7231,
        genes_filter_fun = function(x) {log2(mean(10*(2^x - 1)) + 1) >= 4 | sum(x >= 7) >= length(x)/100},
        scores_filter_fun = function(x) {log2(mean(10*(2^x - 1)) + 1) >= 5.5 | sum(x >= 7) >= length(x)/100}
    ),
    hnsc = list(
        seed = 8511,
        genes_filter_fun = function(x) {log2(mean(10*(2^x - 1)) + 1) >= 4 | sum(x >= 7) >= length(x)/100},
        scores_filter_fun = function(x) {log2(mean(10*(2^x - 1)) + 1) >= 6 | sum(x >= 7) >= length(x)/100}
    ),
    lihc = list(
        seed = 4376,
        genes_filter_fun = function(x) {log2(mean(10*(2^x - 1)) + 1) >= 4 | sum(x >= 7) >= length(x)/100},
        scores_filter_fun = function(x) {log2(mean(10*(2^x - 1)) + 1) >= 6 | sum(x >= 7) >= length(x)/100}
    ),
	lihc_lenient = list(
        seed = 240,
        genes_filter_fun = function(x) {log2(mean(10*(2^x - 1)) + 1) >= 4 | sum(x >= 7) >= length(x)/100},
        scores_filter_fun = function(x) {log2(mean(10*(2^x - 1)) + 1) >= 6 | sum(x >= 7) >= length(x)/100}
    ),
	luad = list(
        seed = 1096,
        genes_filter_fun = function(x) {log2(mean(10*(2^x - 1)) + 1) >= 4 | sum(x >= 7) >= length(x)/100},
        scores_filter_fun = function(x) {log2(mean(10*(2^x - 1)) + 1) >= 5.5 | sum(x >= 7) >= length(x)/100}
    ),
	# luad_lenient = list(
        # seed = 5226,
        # genes_filter_fun = function(x) {log2(mean(10*(2^x - 1)) + 1) >= 4.5 | sum(x >= 7) >= length(x)/100},
        # scores_filter_fun = function(x) {log2(mean(10*(2^x - 1)) + 1) >= 5.5 | sum(x >= 7) >= length(x)/100}
    # ),
    lusc = list(
        seed = 2566,
        genes_filter_fun = function(x) {log2(mean(10*(2^x - 1)) + 1) >= 4 | sum(x >= 7) >= length(x)/100},
        scores_filter_fun = function(x) {log2(mean(10*(2^x - 1)) + 1) >= 5.5 | sum(x >= 7) >= length(x)/100}
    ),
	lusc_lenient = list(
        seed = 9944,
        genes_filter_fun = function(x) {log2(mean(10*(2^x - 1)) + 1) >= 4 | sum(x >= 7) >= length(x)/100},
        scores_filter_fun = function(x) {log2(mean(10*(2^x - 1)) + 1) >= 5.5 | sum(x >= 7) >= length(x)/100}
    ),
    ov = list(
        seed = 456,
        genes_filter_fun = function(x) {log2(mean(10*(2^x - 1)) + 1) >= 4 | sum(x >= 7) >= length(x)/100},
        scores_filter_fun = function(x) {log2(mean(10*(2^x - 1)) + 1) >= 5.5 | sum(x >= 7) >= length(x)/100}
    ),
	ov_lenient = list(
        seed = 5412,
        genes_filter_fun = function(x) {log2(mean(10*(2^x - 1)) + 1) >= 4 | sum(x >= 7) >= length(x)/100},
        scores_filter_fun = function(x) {log2(mean(10*(2^x - 1)) + 1) >= 5.5 | sum(x >= 7) >= length(x)/100}
    ),
    paad = list(
        seed = 5368,
        genes_filter_fun = function(x) {log2(mean(10*(2^x - 1)) + 1) >= 4.5 | sum(x >= 7) >= length(x)/100},
        scores_filter_fun = function(x) {log2(mean(10*(2^x - 1)) + 1) >= 6 | sum(x >= 7) >= length(x)/100}
    )
)





# Construct heatmap data:

# If already done, skip by reading from RDS:

# sc_cancer_caf <- readRDS('../data_and_figures/sc_cancer_caf.rds')

sc_cancer_caf <- sapply(
    names(sc_metadata),
    function(ct) {
        cat(paste0(ct, '\n'))
        sc_data <- eval(sc_metadata[[ct]]$read_quote)
        set.seed(sc_cancer_caf_args[[ct]]$seed)
        do.call(
            sc_groups,
            args = c(
                list(
                    genes = genes_list[[ct]]$filtered_cancer_caf,
                    sc_data = sc_data[cell_type %in% c('cancer', 'caf')],
					groups = c('cancer', 'caf'),
					score_cells_nbin = 30,
					score_cells_n = 40,
                    min_sig_size = 0,
                    scores_filter_groups = 'cancer' # Better for removing effect of complexity
                ),
                sc_cancer_caf_args[[ct]][-1]
            )
        )
    },
    simplify = FALSE,
    USE.NAMES = TRUE
)

# Use the following to check if we need to reverse the order of the genes to make the genes more highly expressed by CAFs to appear
# at the bottom of the heatmap.  For each cancer type, it gives the average expression of the 20 genes in the head and tail of the
# list.  We want the "head" avearge to be highest in the CAFs and the "tail" average to be higher in the cancer cells.
check_ends <- rbindlist(
	lapply(
		names(sc_cancer_caf),
		function(ct) {
			cbind(
				ct,
				melt(
					sc_cancer_caf[[ct]]$data[, -c('genes_detected', 'patient')],
					id.vars = c('id', 'cell_type'),
					variable.name = 'gene',
					value.name = 'expval'
				)[
					,
					.(meanexp = mean(expval)),
					by = .(gene, cell_type)
				][
					,
					.SD[sc_cancer_caf[[ct]]$ordering_genes],
					keyby = cell_type
				][
					c('caf', 'cancer'),
					.(end = c('head', 'tail'), ave_exp = c(mean(head(meanexp, 20)), mean(tail(meanexp, 20)))),
					by = cell_type
				]
			)
		}
	)
)

sc_cancer_caf$ov_lenient$ordering_genes <- rev(sc_cancer_caf$ov_lenient$ordering_genes)

# Save to RDS:

saveRDS(sc_cancer_caf, '../data_and_figures/sc_cancer_caf.rds')





# Define parameters for making heatmaps:

sc_cancer_caf_heatmaps_args <- list(
	brca = list(
		annotations = c('COL1A2', 'ECM1', 'FN1', 'NOTCH2', 'RHOB', 'SNAI1', 'SNAI2', 'THY1', 'TWIST1', 'VCAN', 'VEGFA', 'VIM', 'ZEB1', 'ZEB2'),
		annotations_title = 'Breast'
	),
	brca_lenient = list(
		annotations = c('COL1A2', 'ECM1', 'FN1', 'NOTCH2', 'RHOB', 'SNAI1', 'SNAI2', 'THY1', 'TWIST1', 'VCAN', 'VEGFA', 'VIM', 'ZEB1', 'ZEB2'),
		annotations_title = 'Breast'
	),
	coadread = list(
		annotations = c('AREG', 'CALU', 'COL1A1', 'CXCL1', 'FN1', 'GADD45B', 'PLOD3', 'SNAI1', 'SNAI2', 'THY1', 'TWIST1', 'VIM', 'ZEB1', 'ZEB2'),
		annotations_title = 'Colorectal',
		annotations_side = 'right'
	),
	coadread_lenient = list(
		annotations = c('AREG', 'CALU', 'COL1A1', 'CXCL1', 'FN1', 'GADD45B', 'PLOD3', 'SNAI1', 'SNAI2', 'THY1', 'TWIST1', 'VIM', 'ZEB1', 'ZEB2'),
		annotations_title = 'Colorectal',
		annotations_side = 'right'
	),
    hnsc = list(
		annotations = c('ACTA2', 'CALU', 'COL1A1', 'FN1', 'LAMC2', 'SDC1', 'SNAI1', 'SNAI2', 'SPARC', 'TNC', 'TWIST1', 'VIM', 'ZEB1', 'ZEB2'),
        annotations_title = 'Head and Neck'
    ),
	lihc = list(
		annotations = c('ACTA2', 'COL6A2', 'FN1', 'ITGB1', 'LAMC2', 'QSOX1', 'SNAI1', 'SNAI2', 'TIMP1', 'TWIST1', 'VEGFA', 'VIM', 'ZEB1', 'ZEB2'),
		annotations_title = 'Liver',
        annotations_side = 'right'
	),
	lihc_lenient = list(
		annotations = c('ACTA2', 'COL6A2', 'FN1', 'ITGB1', 'LAMC2', 'QSOX1', 'SNAI1', 'SNAI2', 'TIMP1', 'TWIST1', 'VEGFA', 'VIM', 'ZEB1', 'ZEB2'),
		annotations_title = 'Liver'
	),
	luad = list(
		annotations = c('AREG', 'CALU', 'COL1A1', 'FN1', 'RHOB', 'SDC4', 'SNAI1', 'SNAI2', 'SPARC', 'TWIST1', 'VEGFA', 'VIM', 'ZEB1', 'ZEB2'),
		annotations_title = 'Lung Adenocarcinoma'
	),
	# luad_lenient = list(
		# annotations = c('CALU', 'VEGFA', 'AREG', 'PVR', 'QSOX1', 'SNAI1', 'SNAI2', 'TWIST1', 'ZEB1', 'ZEB2', 'VIM', 'TPM2', 'FN1', 'COL1A2'),
		# annotations_title = 'Lung Adeno.',
		# annotations_side = 'right'
	# ),
    lusc = list(
        annotations = c('CALU', 'COL1A1', 'DCN', 'FN1', 'IGFBP2', 'SDC1', 'SNAI1', 'SNAI2', 'TAGLN', 'TNC', 'TWIST1', 'VIM', 'ZEB1', 'ZEB2'),
        annotations_title = 'Lung Squamous',
		annotations_side = 'right'
    ),
	lusc_lenient = list(
        annotations = c('CALU', 'COL1A1', 'DCN', 'FN1', 'IGFBP2', 'SDC1', 'SNAI1', 'SNAI2', 'TAGLN', 'TNC', 'TWIST1', 'VIM', 'ZEB1', 'ZEB2'),
        annotations_title = 'Lung Squamous'
    ),
	ov = list(
		annotations = c('CALU', 'COL1A1', 'COL1A2', 'FN1', 'ITGB1', 'LAMC1', 'QSOX1', 'SNAI1', 'SNAI2', 'TWIST1', 'VCAN', 'VIM', 'ZEB1', 'ZEB2'),
		annotations_title = 'Ovarian'
	),
	ov_lenient = list(
		annotations = c('CALU', 'COL1A1', 'COL1A2', 'FN1', 'ITGB1', 'LAMC1', 'QSOX1', 'SNAI1', 'SNAI2', 'TWIST1', 'VCAN', 'VIM', 'ZEB1', 'ZEB2'),
		annotations_title = 'Ovarian',
		annotations_side = 'right'
	),
    paad = list(
        annotations = c('COL1A1', 'COL1A2', 'FN1', 'LAMC2', 'QSOX1', 'SDC4', 'SNAI1', 'SNAI2', 'TIMP1', 'TWIST1', 'VEGFA', 'VIM', 'ZEB1', 'ZEB2'),
        annotations_title = 'Pancreatic',
		annotations_side = 'right'
    )
)





# Heatmaps to be plotted one page at a time:

sc_cancer_caf_heatmaps <- sapply(
    names(sc_metadata),
    function(ct) {
		cat(paste0(ct, '\n'))
        set.seed(sc_cancer_caf_args[[ct]]$seed)
        do.call(
            sc_groups_heatmap,
            args = c(
                list(
                    sc_groups_list = sc_cancer_caf[[ct]],
                    groups = c('cancer', 'caf'),
                    default_figure_widths = list(
                        annotations = 1.5,
                        cancer = 6,
                        caf = 1.2
                    ),
                    figure_spacing = 2.5,
                    annotations_nudge = 0.25,
                    es_fun = NULL,
                    es_title = 'EMT score',
                    h_legend_title = 'Expression level',
                    gd_legend_title = 'Genes detected'
                ),
                sc_cancer_caf_heatmaps_args[[ct]][1:2]
            )
        )
    },
    simplify = FALSE,
    USE.NAMES = TRUE
)

cairo_pdf(
    '../data_and_figures/sc_heatmaps.pdf',
    width = 8,
    height = 5,
    onefile = TRUE
)

for(ct in names(sc_cancer_caf_heatmaps)[!grepl('_lenient', names(sc_cancer_caf_heatmaps))]) {
	cowplot_sc(
		sc_cancer_caf_heatmaps[[ct]],
		legend_space = 0.2,
		heights = c(1.5, 20, 4),
		legend_rel_heights = c(2, 5),
		es_x_axis_title_vjust = 1.3,
		es_y_axis_title_angle = 0,
		es_y_axis_title_xpos = 0.8,
		x_axis_titles_space = 0.3
	) %>% print
}

dev.off()

cairo_pdf(
    '../data_and_figures/sc_heatmaps_lenient.pdf',
    width = 8,
    height = 5,
    onefile = TRUE
)

for(ct in names(sc_cancer_caf_heatmaps)[grepl('_lenient', names(sc_cancer_caf_heatmaps))]) {
	cowplot_sc(
		sc_cancer_caf_heatmaps[[ct]],
		legend_space = 0.2,
		heights = c(1.5, 20, 4),
		legend_rel_heights = c(2, 5),
		es_x_axis_title_vjust = 1.3,
		es_y_axis_title_angle = 0,
		es_y_axis_title_xpos = 0.8,
		x_axis_titles_space = 0.3
	) %>% print
}

dev.off()





# Heatmaps to be made into one combined figure:

sc_cancer_caf_heatmaps <- sapply(
    names(sc_metadata),
    function(ct) {
		cat(paste0(ct, '\n'))
        set.seed(sc_cancer_caf_args[[ct]]$seed)
        do.call(
            sc_groups_heatmap,
            args = c(
                list(
                    sc_groups_list = sc_cancer_caf[[ct]],
                    groups = c('cancer', 'caf'),
                    default_figure_widths = list(
                        annotations = 2.5,
                        cancer = 6,
                        caf = 1.2
                    ),
                    figure_spacing = 2.5,
                    annotations_nudge = 0.25,
                    es_fun = NULL,
                    es_title = 'EMT score',
                    h_legend_title = 'Expression level\n',
                    h_legend_width = 20,
                    h_legend_height = 10,
                    h_legend_direction = 'horizontal',
                    h_legend_title_position = 'right',
                    h_legend_just = 'left',
                    gd_legend_title = 'Genes detected\n',
                    gd_legend_width = 20,
                    gd_legend_height = 10,
                    gd_legend_direction = 'horizontal',
                    gd_legend_title_position = 'left',
                    gd_legend_just = 'right'
                ),
                sc_cancer_caf_heatmaps_args[[ct]]
            )
        )
    },
    simplify = FALSE,
    USE.NAMES = TRUE
)

cairo_pdf(
    '../data_and_figures/sc_heatmaps_combined.pdf',
    width = 10,
    height = 14
)

plot_grid(
    plot_grid(
        cowplot_sc(
            sc_cancer_caf_heatmaps$brca,
            legend_space = 0,
            heights = c(1.5, 20, 4),
            es_x_axis_title_vjust = 1.3,
            es_y_axis_title_angle = 0,
            es_y_axis_title_xpos = 0.8
        ),
        blank_plot(),
        cowplot_sc(
            sc_cancer_caf_heatmaps$coadread,
            legend_space = 0,
            heights = c(1.5, 20, 4),
            es_x_axis_title_vjust = 1.3,
            es_y_axis_title_angle = 0,
            es_y_axis_title_xpos = 0.2
        ),
        ncol = 3,
        nrow = 1,
        rel_widths = c(20, 1, 20)
    ),
    blank_plot(),
    plot_grid(
        cowplot_sc(
            sc_cancer_caf_heatmaps$hnsc,
            legend_space = 0,
            heights = c(1.5, 20, 4),
            es_x_axis_title_vjust = 1.3,
            es_y_axis_title_angle = 0,
            es_y_axis_title_xpos = 0.8
        ),
        blank_plot(),
        cowplot_sc(
            sc_cancer_caf_heatmaps$lihc,
            legend_space = 0,
            heights = c(1.5, 20, 4),
            es_x_axis_title_vjust = 1.3,
            es_y_axis_title_angle = 0,
            es_y_axis_title_xpos = 0.2
        ),
        ncol = 3,
        nrow = 1,
        rel_widths = c(20, 1, 20)
    ),
    blank_plot(),
	plot_grid(
        cowplot_sc(
            sc_cancer_caf_heatmaps$luad,
            legend_space = 0,
            heights = c(1.5, 20, 4),
            es_x_axis_title_vjust = 1.3,
            es_y_axis_title_angle = 0,
            es_y_axis_title_xpos = 0.8
        ),
        blank_plot(),
        cowplot_sc(
            sc_cancer_caf_heatmaps$lusc,
            legend_space = 0,
            heights = c(1.5, 20, 4),
            es_x_axis_title_vjust = 1.3,
            es_y_axis_title_angle = 0,
            es_y_axis_title_xpos = 0.2
        ),
        ncol = 3,
        nrow = 1,
        rel_widths = c(20, 1, 20)
    ),
	blank_plot(),
	plot_grid(
        cowplot_sc(
            sc_cancer_caf_heatmaps$ov,
            legend_space = 0,
            heights = c(1.5, 20, 4),
            es_x_axis_title_vjust = 1.3,
            es_y_axis_title_angle = 0,
            es_y_axis_title_xpos = 0.8
        ),
        blank_plot(),
        cowplot_sc(
            sc_cancer_caf_heatmaps$paad,
            legend_space = 0,
            heights = c(1.5, 20, 4),
            es_x_axis_title_vjust = 1.3,
            es_y_axis_title_angle = 0,
            es_y_axis_title_xpos = 0.2
        ),
        ncol = 3,
        nrow = 1,
        rel_widths = c(20, 1, 20)
    ),
	blank_plot(),
    plot_grid(
        sc_cancer_caf_heatmaps$paad$plots$genes_detected_legend,
        blank_plot(),
        sc_cancer_caf_heatmaps$paad$plots$heatmap_legend,
        nrow = 1,
        ncol = 3,
        rel_widths = c(5, 1, 5)
    ),
    ncol = 1,
    nrow = 9,
    rel_heights = c(20, 1, 20, 1, 20, 1, 20, 1, 4)
) %>% print

dev.off()

# For the lenient classifications:

cairo_pdf(
    '../data_and_figures/sc_heatmaps_lenient_combined.pdf',
    width = 10,
    height = 11
)

plot_grid(
    plot_grid(
        cowplot_sc(
            sc_cancer_caf_heatmaps$brca_lenient,
            legend_space = 0,
            heights = c(1.5, 20, 4),
            es_x_axis_title_vjust = 1.3,
            es_y_axis_title_angle = 0,
            es_y_axis_title_xpos = 0.8
        ),
        blank_plot(),
        cowplot_sc(
            sc_cancer_caf_heatmaps$coadread_lenient,
            legend_space = 0,
            heights = c(1.5, 20, 4),
            es_x_axis_title_vjust = 1.3,
            es_y_axis_title_angle = 0,
            es_y_axis_title_xpos = 0.2
        ),
        ncol = 3,
        nrow = 1,
        rel_widths = c(20, 1, 20)
    ),
    blank_plot(),
    plot_grid(
        cowplot_sc(
            sc_cancer_caf_heatmaps$lihc_lenient,
            legend_space = 0,
            heights = c(1.5, 20, 4),
            es_x_axis_title_vjust = 1.3,
            es_y_axis_title_angle = 0,
            es_y_axis_title_xpos = 0.8
        ),
        blank_plot(),
        cowplot_sc(
            sc_cancer_caf_heatmaps$lusc_lenient,
            legend_space = 0,
            heights = c(1.5, 20, 4),
            es_x_axis_title_vjust = 1.3,
            es_y_axis_title_angle = 0,
            es_y_axis_title_xpos = 0.8
        ),
        ncol = 3,
        nrow = 1,
        rel_widths = c(20, 1, 20)
    ),
    blank_plot(),
	plot_grid(
        cowplot_sc(
            sc_cancer_caf_heatmaps$ov_lenient,
            legend_space = 0,
            heights = c(1.5, 20, 4),
            es_x_axis_title_vjust = 1.3,
            es_y_axis_title_angle = 0,
            es_y_axis_title_xpos = 0.2
        ),
        blank_plot(),
        blank_plot(),
        ncol = 3,
        nrow = 1,
        rel_widths = c(20, 1, 20)
    ),
	blank_plot(),
    plot_grid(
        sc_cancer_caf_heatmaps$paad$plots$genes_detected_legend,
        blank_plot(),
        sc_cancer_caf_heatmaps$paad$plots$heatmap_legend,
        nrow = 1,
        ncol = 3,
        rel_widths = c(5, 1, 5)
    ),
    ncol = 1,
    nrow = 7,
    rel_heights = c(20, 1, 20, 1, 20, 1, 4)
) %>% print

dev.off()





# Final version for main figures 1B/C:

sc_cancer_caf_heatmaps_args_final <- list(
	brca = list(
		annotations = c('COL1A2', 'ECM1', 'FN1', 'NOTCH2', 'RHOB', 'SNAI1', 'SNAI2', 'THY1', 'TWIST1', 'VCAN', 'VEGFA', 'VIM', 'ZEB1', 'ZEB2'),
		annotations_title = 'Breast'
	),
	coadread = list(
		annotations = c('AREG', 'COL1A1', 'CXCL1', 'FN1', 'GADD45B', 'PLOD3', 'SNAI1', 'SNAI2', 'THY1', 'TWIST1', 'VIM', 'ZEB1', 'ZEB2'),
		annotations_title = 'Colorectal'
	),
    hnsc = list(
		annotations = c('ACTA2', 'COL1A1', 'FN1', 'LAMC2', 'SDC1', 'SNAI1', 'SNAI2', 'SPARC', 'TNC', 'TWIST1', 'VIM', 'ZEB1', 'ZEB2'),
        annotations_title = 'Head and Neck',
		annotations_side = 'right'
    ),
	lihc = list(
		annotations = c('ACTA2', 'COL6A2', 'FN1', 'ITGB1', 'LAMC2', 'QSOX1', 'SNAI1', 'SNAI2', 'TIMP1', 'TWIST1', 'VEGFA', 'VIM', 'ZEB1', 'ZEB2'),
		annotations_title = 'Liver',
        annotations_side = 'right'
	),
	luad = list(
		annotations = c('AREG', 'COL1A1', 'FN1', 'RHOB', 'SDC4', 'SNAI1', 'SNAI2', 'SPARC', 'TWIST1', 'VEGFA', 'VIM', 'ZEB1', 'ZEB2'),
		annotations_title = 'Lung Adenocarcinoma'
	),
    lusc = list(
        annotations = c('CALU', 'COL1A1', 'DCN', 'FN1', 'IGFBP2', 'SDC1', 'SNAI1', 'SNAI2', 'TAGLN', 'TNC', 'TWIST1', 'VIM', 'ZEB1', 'ZEB2'),
        annotations_title = 'Lung Squamous'
    ),
	ov = list(
		annotations = c('CALU', 'COL1A1', 'COL1A2', 'FN1', 'ITGB1', 'LAMC1', 'QSOX1', 'SNAI1', 'SNAI2', 'TWIST1', 'VCAN', 'VIM', 'ZEB1', 'ZEB2'),
		annotations_title = 'Ovarian',
		annotations_side = 'right'
	),
    paad = list(
        annotations = c('COL1A1', 'COL1A2', 'FN1', 'LAMC2', 'SDC4', 'SNAI1', 'SNAI2', 'TIMP1', 'TWIST1', 'VEGFA', 'VIM', 'ZEB1', 'ZEB2'),
        annotations_title = 'Pancreatic',
		annotations_side = 'right'
    )
)

sc_cancer_caf_heatmaps_final <- sapply(
    names(sc_cancer_caf_heatmaps_args_final),
    function(ct) {
		cat(paste0(ct, '\n'))
        set.seed(sc_cancer_caf_args[[ct]]$seed)
        do.call(
            sc_groups_heatmap,
            args = c(
                list(
                    sc_groups_list = sc_cancer_caf[[ct]],
                    groups = c('cancer', 'caf'),
					x_axis_titles = c('Cancer cells', 'CAFs'),
                    default_figure_widths = list(annotations = 2.5, cancer = 6, caf = 1.2),
                    figure_spacing = 2.5,
					annotations_title_size = 18,
                    annotations_nudge = 0.25,
                    es_fun = NULL,
                    es_title = 'EMT score',
                    h_legend_title = 'Expression level',
                    h_legend_width = 20,
                    h_legend_height = 10,
                    h_legend_direction = 'horizontal',
                    h_legend_title_position = 'top',
                    h_legend_just = 'top',
                    gd_legend_title = 'Genes detected',
                    gd_legend_width = 20,
                    gd_legend_height = 10,
                    gd_legend_direction = 'horizontal',
                    gd_legend_title_position = 'top',
                    gd_legend_just = 'top'
                ),
                sc_cancer_caf_heatmaps_args_final[[ct]]
            )
        )
    },
    simplify = FALSE,
    USE.NAMES = TRUE
)

# Make heatmaps showing average and 95th percentile of EMT TF (+ VIM) expression (in all cancer types):

emt_tf_exp <- sapply(
    c('mean', '95th percentile'),
    function(fct) {
        rbindlist(
			sapply(
				c('cancer', 'caf'),
				function(x) {
					melt(
						as.data.table(
							sapply(
								names(sc_cancer_caf)[!grepl('_lenient', names(sc_cancer_caf))],
								function(ct) {
									sc_cancer_caf[[ct]]$data[
										cell_type == x,
										apply(.SD, 2, switch((fct == 'mean') + 1, function(y) quantile(y, 0.95), mean)),
										.SDcols = c('SNAI1', 'SNAI2', 'TWIST1', 'VIM', 'ZEB1', 'ZEB2')
									]
								},
								USE.NAMES = TRUE
							),
							keep.rownames = 'gene'
						),
						id.var = 'gene',
						variable.name = 'cancer_type',
						value.name = 'expression_level'
					)[, cell_type := x]
				},
				simplify = FALSE,
				USE.NAMES = TRUE
			)
		)
    },
    simplify = FALSE,
    USE.NAMES = TRUE
)

# Make VIM appear at the top of the heatmaps:
for(dt in emt_tf_exp) {dt[, gene := factor(gene, levels = c('SNAI1', 'SNAI2', 'TWIST1', 'ZEB1', 'ZEB2', 'VIM'))]}

emt_tf_exp_plots <- sapply(
    names(emt_tf_exp),
    function(fct) {
        ggplot(
			emt_tf_exp[[fct]],
			aes(
				x = plyr::mapvalues(
					cancer_type,
					unique(cancer_type),
					sapply(unique(cancer_type), function(ct) sc_cancer_caf_heatmaps_args_final[[ct]]$annotations_title)
				),
				y = gene,
				fill = expression_level
			)
		) +
			facet_grid(cols = vars(cell_type), labeller = labeller(cell_type = c(cancer = 'Cancer cells', caf = 'CAFs'))) +
			geom_raster() +
			scale_fill_gradientn(
				colours = rev(colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))(50)),
				limits = c(0, 12),
				breaks = c(0, 3, 6, 9, 12),
                labels = c('0' = '0', '3' = '3', '6' = '6', '9' = '9', '12' = '\u2265 12'),
				oob = scales::squish
			) +
			scale_x_discrete(expand = c(0, 0)) +
			scale_y_discrete(expand = c(0, 0)) +
			theme(
				legend.key.width = unit(10, 'pt'),
				axis.text.x = element_text(angle = 55, hjust = 1),
				axis.title = element_blank(),
				plot.title = element_text(hjust = 0.5),
				panel.border = element_rect(colour = 'black', fill = NA),
				strip.background = element_rect(colour = NA, fill = NA),
				strip.text = element_text(size = 11)
			) +
			labs(title = mapvalues(fct, 'mean', 'Average expression', warn_missing = FALSE), fill = 'Expression\nlevel')
    },
    simplify = FALSE,
    USE.NAMES = TRUE
)

cairo_pdf('../data_and_figures/final_figures_resubmission/1BC.pdf', width = 11, height = 12)

plot_grid(
    blank_plot(),
	plot_grid(
		blank_plot(),
        sc_cancer_caf_heatmaps_final$paad$plots$genes_detected_legend,
        sc_cancer_caf_heatmaps_final$paad$plots$heatmap_legend,
		blank_plot(),
        nrow = 1,
        ncol = 4
    ),
	# blank_plot(),
    plot_grid(
        cowplot_sc(
            sc_cancer_caf_heatmaps_final$coadread,
            legend_space = 0,
            heights = c(1.5, 20, 4),
            es_x_axis_title_vjust = 1.4,
            es_y_axis_title_angle = 0,
            es_y_axis_title_xpos = 0.82
        ),
        blank_plot(),
        cowplot_sc(
            sc_cancer_caf_heatmaps_final$hnsc,
            legend_space = 0,
            heights = c(1.5, 20, 4),
            es_x_axis_title_vjust = 1.4,
            es_y_axis_title_angle = 0,
            es_y_axis_title_xpos = 0.18
        ),
        ncol = 3,
        nrow = 1,
        rel_widths = c(10, 0.8, 10)
    ),
    blank_plot(),
    plot_grid(
        cowplot_sc(
            sc_cancer_caf_heatmaps_final$luad,
            legend_space = 0,
            heights = c(1.5, 20, 4),
            es_x_axis_title_vjust = 1.4,
            es_y_axis_title_angle = 0,
            es_y_axis_title_xpos = 0.82
        ),
        blank_plot(),
        cowplot_sc(
            sc_cancer_caf_heatmaps_final$paad,
            legend_space = 0,
            heights = c(1.5, 20, 4),
            es_x_axis_title_vjust = 1.4,
            es_y_axis_title_angle = 0,
            es_y_axis_title_xpos = 0.18
        ),
        ncol = 3,
        nrow = 1,
        rel_widths = c(10, 0.8, 10)
    ),
    blank_plot(),
    plot_grid(
		plotlist = c(
			# list(blank_plot()),
			list(
				ggplot(data.frame(y = 1:100), aes(x = 0:1, y = y)) +
					theme(
						axis.text = element_blank(),
						axis.ticks = element_blank(),
						axis.ticks.length = unit(0, 'pt'),
						axis.title = element_blank(),
						plot.background = element_blank(),
						panel.background = element_blank(),
						plot.margin = unit(c(0, 0, 0, 0), 'pt')
					) +
					scale_x_continuous(expand = c(0, 0), limits = c(0, 1)) +
					scale_y_continuous(expand = c(0, 0), limits = c(0, 100)) +
					geom_segment(aes(x = 0.85, xend = 0.85, y = 38, yend = 74)) +
					geom_segment(aes(x = 0.85, xend = 1, y = 38, yend = 38)) +
					geom_segment(aes(x = 0.85, xend = 1, y = 74, yend = 74)) +
					annotate(geom = 'text', x = 0.6, y = 56, label = 'EMT TFs', angle = 90)
			),
			lapply(emt_tf_exp_plots, function(x) x + theme(legend.position = 'none')),
			list(get_legend(emt_tf_exp_plots[[1]] + theme(legend.justification = c(0.1, 0.75))))
		),
		nrow = 1,
		ncol = 4,
		rel_widths = c(0.5, 4, 4, 1.5)
	),
    ncol = 1,
    nrow = 7,
    rel_heights = c(1, 2, 10, 1, 10, 1.9, 10)
) +
    draw_label('B', x = 0, y = 0.99, hjust = -0.1, vjust = 1.1, size = 20, fontface = 2) +
    draw_label('C', x = 0, y = 0.3, hjust = -0.1, vjust = 1.1, size = 20, fontface = 2)

dev.off()





# Check heatmaps against patient IDs and subclone information:

patient_subclone_bars <- sapply(
	names(sc_metadata),
	function(ct) {

		cat(paste0(ct, '\n'))

		sc_data <- sc_cancer_caf[[ct]]$data
		cells <- sc_cancer_caf[[ct]]$cells_filtered
		ordering_cells <- sc_cancer_caf[[ct]]$ordering_cells
		setkey(sc_data, id)
		annotations_side <- sc_cancer_caf_heatmaps_args[[ct]]$annotations_side
		if(is.null(annotations_side)) {annotations_side <- 'left'}
		# random_colours <- randomcoloR::distinctColorPalette(length(unique(sc_data$patient)))

		pbars <- sapply(
			names(ordering_cells),
			function(grp) {
				ggplot(
					sc_data[cells][cell_type == grp, .(id, patient)],
					aes(x = factor(id, levels = id[ordering_cells[[grp]]]), y = 0, fill = as.character(patient))
				) +
					geom_raster() +
					# scale_fill_manual(values = random_colours) +
					scale_x_discrete(expand = c(0, 0)) +
					scale_y_continuous(expand = c(0, 0)) +
					theme(
						axis.text = element_blank(),
						axis.ticks.length = unit(0, 'pt'),
						axis.ticks = element_blank(),
						axis.title = element_blank(),
						plot.margin = unit(
							switch(
								(annotations_side == 'left') + 1,
								c(2.5, 0, 2.5, switch((grp == 'cancer') + 1, 2.5, 0)),
								c(2.5, switch((grp == 'caf') + 1, 2.5, 0), 2.5, 0)
							),
							'pt'
						),
						legend.justification = 'top'
					) +
					labs(fill = 'Patient')
			},
			simplify = FALSE,
			USE.NAMES = TRUE
		)

		cdata <- rbindlist(
			lapply(
				unique(sc_data$patient),
				function(p) unique(
					fread(
						paste0(
							'../data_and_figures/sc_find_malignant/',
							mapvalues(
								gsub('_lenient', '', ct),
								c('brca', 'coadread', 'hnsc', 'lihc', 'luad', 'lusc', 'ov', 'paad'),
								c('breast_qian', 'crc_lee_smc', 'hnscc_puram', 'liver_ma', 'luad_kim', 'lung_qian', 'ovarian_qian', 'pdac_peng'),
								warn_missing = FALSE
							),
							'/',
							p,
							'_data.csv'
						)
					)[, .(cell_id, classification_final)]
				)
			)
		)

		cdata[
			,
			classification := switch(
				(length(classification_final) == 1) + 1,
				switch(('ambiguous' %in% classification_final | 'malignant' %in% classification_final) + 1, 'nonmalignant', 'ambiguous'),
				classification_final
			),
			by = cell_id
		]

		cdata <- unique(cdata[, .(cell_id, classification)])
		setkey(cdata, cell_id)

		subclones <- cdata[, sort(unique(classification[grep('[0-9]$', classification)]))]

		if(length(subclones) == 0) {
			sbar = blank_plot()
			sbar_legend = blank_plot()
		} else {

			subclone_colours <- c(malignant = 'white', setNames(c('coral', 'cornflowerblue', 'seagreen3')[1:length(subclones)], subclones))

			sbar <- ggplot(
				cdata[sc_data[cells][cell_type == 'cancer', id], .(cell_id, classification)],
				aes(x = factor(cell_id, levels = cell_id[ordering_cells$cancer]), y = 0, fill = classification)
			) +
				geom_raster() +
				scale_fill_manual(values = subclone_colours) +
				# scale_fill_manual(values = c('malignant' = 'white', 'malignant clone 1' = 'coral', 'malignant clone 2' = 'cornflowerblue')) +
				scale_x_discrete(expand = c(0, 0)) +
				scale_y_continuous(expand = c(0, 0)) +
				theme(
					axis.text = element_blank(),
					axis.ticks.length = unit(0, 'pt'),
					axis.ticks = element_blank(),
					axis.title = element_blank(),
					plot.margin = unit(switch((annotations_side == 'left') + 1, c(2.5, 0, 2.5, 0), c(2.5, 2.5, 2.5, 0)), 'pt')
				)

			sbar_legend <- get_legend(
				ggplot(
					data.table(cell_id = letters[1:length(subclones)], Subclone = as.character(1:length(subclones))),
					aes(x = cell_id, y = 0, fill = Subclone)
				) +
					geom_raster() +
					scale_fill_manual(values = setNames(subclone_colours[subclones], as.character(1:length(subclones)))) +
					# scale_fill_manual(values = c('1' = 'coral', '2' = 'cornflowerblue')) +
					theme(legend.justification = 'top')
			)

		}

		list(patient_bars = pbars, subclone_bar = sbar, subclone_bar_legend = sbar_legend, classification_data = cdata)

	},
	simplify = FALSE,
	USE.NAMES = TRUE
)

cairo_pdf('../data_and_figures/sc_heatmaps_patient_subclone.pdf', width = 8, height = 5, onefile = TRUE)

for(ct in names(sc_metadata)[!grepl('lenient', names(sc_metadata))]) {
	plot_grid(
		plot_grid(
			plotlist = c(
				list(
					blank_plot() +
						theme(plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), 'pt')) +
						labs(title = sc_cancer_caf_heatmaps_args[[ct]]$annotations_title),
					blank_plot()
				),
				sc_cancer_caf_heatmaps[[ct]]$plots$genes_detected,
				lapply(patient_subclone_bars[[ct]]$patient_bars, function(x) x + theme(legend.position = 'none')),
				list(patient_subclone_bars[[ct]]$subclone_bar + theme(legend.position = 'none'), blank_plot()),
				sc_cancer_caf_heatmaps[[ct]]$plots$heatmaps
			),
			nrow = 5,
			ncol = 2,
			rel_heights = c(2, 1.5, 1.5, 1.5, 20),
			rel_widths = unlist(sc_cancer_caf_heatmaps[[ct]]$figure_widths[c('cancer', 'caf')])
		),
		plot_grid(
			get_legend(patient_subclone_bars[[ct]]$patient_bars[[1]]),
			patient_subclone_bars[[ct]]$subclone_bar_legend,
			nrow = 1,
			ncol = 2
		),
		nrow = 1,
		ncol = 2,
		rel_widths = c(4, 1.5)
	) %>% print
}

dev.off()

cairo_pdf('../data_and_figures/sc_heatmaps_lenient_patient_subclone.pdf', width = 8, height = 5, onefile = TRUE)

for(ct in names(sc_metadata)[grep('lenient', names(sc_metadata))]) {
	plot_grid(
		plot_grid(
			plotlist = c(
				list(
					blank_plot() +
						theme(plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), 'pt')) +
						labs(title = sc_cancer_caf_heatmaps_args[[ct]]$annotations_title),
					blank_plot()
				),
				sc_cancer_caf_heatmaps[[ct]]$plots$genes_detected,
				lapply(patient_subclone_bars[[ct]]$patient_bars, function(x) x + theme(legend.position = 'none')),
				list(patient_subclone_bars[[ct]]$subclone_bar + theme(legend.position = 'none'), blank_plot()),
				sc_cancer_caf_heatmaps[[ct]]$plots$heatmaps
			),
			nrow = 5,
			ncol = 2,
			rel_heights = c(2, 1.5, 1.5, 1.5, 20),
			rel_widths = unlist(sc_cancer_caf_heatmaps[[ct]]$figure_widths[c('cancer', 'caf')])
		),
		plot_grid(
			get_legend(patient_subclone_bars[[ct]]$patient_bars[[1]]),
			patient_subclone_bars[[ct]]$subclone_bar_legend,
			nrow = 1,
			ncol = 2
		),
		nrow = 1,
		ncol = 2,
		rel_widths = c(4, 1.5)
	) %>% print
}

dev.off()





# Supplementary figure showing remaining cancer types and patient/subclone bars for breast and ovarian:

patient_subclone_bars_final <- sapply(
	c('brca', 'ov'),
	function(ct) {

		cat(paste0(ct, '\n'))

		sc_data <- sc_cancer_caf[[ct]]$data[cell_type == 'cancer', .(id, patient)]
		ordering_cells <- sc_cancer_caf[[ct]]$ordering_cells$cancer
		setkey(sc_data, id)

		cdata <- rbindlist(
			lapply(
				unique(sc_data$patient),
				function(p) unique(
					fread(
						paste0(
							'../data_and_figures/sc_find_malignant/',
							mapvalues(
								gsub('_lenient', '', ct),
								c('brca', 'coadread', 'hnsc', 'lihc', 'luad', 'lusc', 'ov', 'paad'),
								c('breast_qian', 'crc_lee_smc', 'hnscc_puram', 'liver_ma', 'luad_kim', 'lung_qian', 'ovarian_qian', 'pdac_peng'),
								warn_missing = FALSE
							),
							'/',
							p,
							'_data.csv'
						)
					)[, .(cell_id, classification_final)]
				)
			)
		)

		cdata[
			,
			classification := switch(
				(length(classification_final) == 1) + 1,
				switch(('ambiguous' %in% classification_final | 'malignant' %in% classification_final) + 1, 'nonmalignant', 'ambiguous'),
				classification_final
			),
			by = cell_id
		]

		cdata <- unique(cdata[, .(cell_id, classification)])
		setkey(cdata, cell_id)

		sc_data[
			,
			patient_clone := switch(
				(cdata[id, classification] == 'malignant') + 1,
				paste(as.character(patient), gsub('malignant ', '', cdata[id, classification]), sep = ' - '),
				as.character(patient)
			),
			by = id
		]

		# The colours in the following are from the CARTO Bold and Pastel palettes.
		pbar <- ggplot(
			sc_data,
			aes(x = factor(id, levels = id[ordering_cells]), y = 0, fill = patient_clone)
		) +
			geom_raster() +
			scale_fill_manual(
				values = switch(
					(ct == 'brca') + 1,
					c(
						'11 - clone 1' = '#F2B701',
						'11 - clone 2' = '#E68310',
						'14' = '#E73F74',
						setNames(
							lighten(c('#66C5CC', '#DCB0F2', '#8BE0A4', '#B3B3B3'), 0.5),
							sc_data[!(patient_clone %in% c('11 - clone 1', '11 - clone 2', '14')), unique(patient_clone)]
						)
					),
					c(
						'47' = '#E73F74',
						'49 - clone 1' = '#E68310',
						'49 - clone 2' = '#F2B701',
						setNames(
							lighten(c('#66C5CC', '#DCB0F2', '#87C55F', '#9EB9F3', '#8BE0A4', '#B3B3B3'), 0.5), # '#B497E7'
							sc_data[!(patient_clone %in% c('47', '49 - clone 1', '49 - clone 2')), unique(patient_clone)]
						)
					)
				)
			) +
			scale_x_discrete(expand = c(0, 0)) +
			scale_y_continuous(expand = c(0, 0)) +
			theme(
				axis.text = element_blank(),
				axis.ticks.length = unit(0, 'pt'),
				axis.ticks = element_blank(),
				axis.title = element_blank(),
				plot.margin = unit(c(5.5, 5.5, 2.5, 5.5), 'pt'),
				legend.justification = 'top'
			) +
			labs(fill = 'Patient', title = sc_cancer_caf_heatmaps_args_final[[ct]]$annotations_title)

		list(patient_bar = pbar, classification_data = sc_data)

	},
	simplify = FALSE,
	USE.NAMES = TRUE
)

cairo_pdf('../data_and_figures/final_figures_resubmission/S4.pdf', width = 11, height = 12)

plot_grid(
    blank_plot(),
	plot_grid(
		blank_plot(),
        sc_cancer_caf_heatmaps_final$brca$plots$genes_detected_legend,
        sc_cancer_caf_heatmaps_final$brca$plots$heatmap_legend,
		blank_plot(),
        nrow = 1,
        ncol = 4
    ),
	# blank_plot(),
    plot_grid(
        cowplot_sc(
            sc_cancer_caf_heatmaps_final$brca,
            legend_space = 0,
            heights = c(1.5, 20, 4),
            es_x_axis_title_vjust = 1.4,
            es_y_axis_title_angle = 0,
            es_y_axis_title_xpos = 0.82
        ),
        blank_plot(),
        cowplot_sc(
            sc_cancer_caf_heatmaps_final$lihc,
            legend_space = 0,
            heights = c(1.5, 20, 4),
            es_x_axis_title_vjust = 1.4,
            es_y_axis_title_angle = 0,
            es_y_axis_title_xpos = 0.18
        ),
        ncol = 3,
        nrow = 1,
        rel_widths = c(10, 0.8, 10)
    ),
    blank_plot(),
    plot_grid(
        cowplot_sc(
            sc_cancer_caf_heatmaps_final$lusc,
            legend_space = 0,
            heights = c(1.5, 20, 4),
            es_x_axis_title_vjust = 1.4,
            es_y_axis_title_angle = 0,
            es_y_axis_title_xpos = 0.82
        ),
        blank_plot(),
        cowplot_sc(
            sc_cancer_caf_heatmaps_final$ov,
            legend_space = 0,
            heights = c(1.5, 20, 4),
            es_x_axis_title_vjust = 1.4,
            es_y_axis_title_angle = 0,
            es_y_axis_title_xpos = 0.18
        ),
        ncol = 3,
        nrow = 1,
        rel_widths = c(10, 0.8, 10)
    ),
    blank_plot(),
    plot_grid(
		plot_grid(
			patient_subclone_bars_final$brca$patient_bar + theme(legend.position = 'none'),
			sc_cancer_caf_heatmaps_final$brca$plots$heatmaps$cancer +
				theme(axis.title.x = element_text(), axis.title.y = element_text(), plot.margin = unit(c(0, 5.5, 5.5, 40), 'pt')) +
				labs(x = 'Cancer cells', y = 'Genes'),
			nrow = 2,
			ncol = 1,
			align = 'v',
			rel_heights = c(1, 5)
		),
		get_legend(patient_subclone_bars_final$brca$patient_bar),
		plot_grid(
			patient_subclone_bars_final$ov$patient_bar + theme(legend.position = 'none'),
			sc_cancer_caf_heatmaps_final$ov$plots$heatmaps$cancer +
				theme(axis.title.x = element_text(), axis.title.y = element_text(), plot.margin = unit(c(0, 5.5, 5.5, 40), 'pt')) +
				labs(x = 'Cancer cells', y = 'Genes'),
			nrow = 2,
			ncol = 1,
			align = 'v',
			rel_heights = c(1, 5)
		),
		get_legend(patient_subclone_bars_final$ov$patient_bar),
		nrow = 1,
		ncol = 4,
		rel_widths = c(3, 1, 3, 1)
	),
    ncol = 1,
    nrow = 7,
    rel_heights = c(1, 2, 10, 1, 10, 1.9, 10)
) +
    draw_label('A', x = 0, y = 0.99, hjust = -0.1, vjust = 1.1, size = 20, fontface = 2) +
    draw_label('B', x = 0, y = 0.3, hjust = -0.1, vjust = 1.1, size = 20, fontface = 2)

dev.off()





# Plots to show average expression of mesenchymal and epithelial genes:

# For this I need the analysis of EMT genes in just the cancer cells.  I was hoping
# not to have to do this - probably I can reengineer everything so I don't.

# sc_cancer <- readRDS('../data_and_figures/sc_cancer.rds')

sc_cancer <- sapply(
    names(sc_cancer_caf_args)[!grepl('_lenient', names(sc_cancer_caf_args))],
    function(ct) {
        cat(paste0(ct, '\n'))
        sc_data <- eval(sc_metadata[[ct]]$read_quote)
        set.seed(sc_cancer_caf_args[[ct]]$seed)
        do.call(
            sc_groups,
            args = c(
                list(
                    genes = genes_list[[ct]]$unfiltered,
                    sc_data = sc_data[cell_type == 'cancer'],
                    groups = 'cancer',
                    to_keep = character(0),
					score_cells_nbin = 30,
					score_cells_n = 40,
                    min_sig_size = 0
                ),
                sc_cancer_caf_args[[ct]][-1]
            )
        )
    },
    simplify = FALSE,
    USE.NAMES = TRUE
)

saveRDS(sc_cancer, '../data_and_figures/sc_cancer.rds')

emt_epi_comparison <- sapply(
    names(sc_cancer),
    function(ct) {

        cat(paste0(ct, '...'))

        sc_data <- eval(sc_metadata[[ct]]$read_quote)[, -'cell_type_lenient']
        setkey(sc_data, id)

        all_genes_filtered <- sc_data[
            sc_cancer[[ct]]$cells_filtered,
            names(.SD)[apply(.SD, 2, sc_cancer_caf_args[[ct]]$genes_filter_fun)],
            .SDcols = -c('id', 'patient', 'cell_type')
        ]

        epi_markers <- names(sc_data)[grep('^CDH1$|^EPCAM$|^SFN$|^KRT[0-9]', names(sc_data))]
        epi_markers <- epi_markers[epi_markers %in% all_genes_filtered]

        plot_data <- sc_data[
            sc_cancer[[ct]]$cells_filtered,
            c(
                .(
                    id = id,
                    emt_score = sc_cancer[[ct]]$scores[id, score],
					epi_score = signature_score(set_colnames(t(.SD[, ..all_genes_filtered]), id), epi_markers, nbin = 30, n = 40)
                ),
                .SD[, ..epi_markers]
            ),
            by = patient
        ][order(-emt_score), epi_score_runmean := runmean(epi_score, .N/10)]

        lineplots <- setNames(
            lapply(
                c('emt_score', 'epi_score_runmean'),
                function(v) {
                    ggplot( plot_data[order(-emt_score)], aes(factor(id, levels = id), get(v), group = 1)) +
                        geom_line() +
                        theme_test() +
                        theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
                        labs(x = 'Cells', y = 'Relative average expression')
                }
            ),
            c('emt', 'epi')
        )

        combined_lineplot <- ggplot(
            melt(plot_data[, .(id, emt_score, epi_score_runmean)], id.vars = 'id'),
            aes(factor(id, levels = plot_data[order(-emt_score), id]), value, group = variable, colour = variable)
        ) +
            geom_line() +
            scale_colour_discrete(labels = c('emt_score' = 'EMT score', 'epi_score_runmean' = 'Epithelial score')) +
            theme_test() +
            theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
            labs(x = 'Cells', y = 'Score', colour = NULL)

        epi_emt_corr <- sort(plot_data[, cor(emt_score, .SD)[1, ], .SDcols = epi_markers])

        epi_heatmap <- ggplot(
            melt(
                plot_data[
                    , # Z scores instead of expression levels to make correlation with EMT score clearer
                    c(.(id = id), sapply(.SD, function(x) {(x - mean(x))/sd(x)}, simplify = FALSE, USE.NAMES = TRUE)),
                    .SDcols = epi_markers
                ],
                # plot_data[, c('id', ..epi_markers)],
                id.vars = 'id',
                variable.name = 'gene',
                value.name = 'expression_level'
            ),
            aes(
                x = factor(id, levels = unique(id)[sc_cancer[[ct]]$ordering_cells$cancer]),
                y = factor(
                    gene,
                    # levels = plot_data[, names(sort(colMeans(.SD))), .SDcols = epi_markers]
                    levels = names(epi_emt_corr)
                ),
                fill = expression_level
            )
        ) +
            geom_raster() +
            scale_fill_gradientn(
                limits = c(-2, 2),
                # limits = c(0, 12),
                colours = rev(colorRampPalette(brewer.pal(11, "RdBu"))(50)),
                # colours = rev(colorRampPalette(brewer.pal(11, "Spectral"))(50)),
                oob = scales::squish,
                breaks = c(-2, -1, 0, 1, 2),
                labels = c('-2' = '\u2264 -2', '-1' = '-1', '0' = '0', '1' = '1', '2' = '\u2265 2')
                # breaks = c(0, 3, 6, 9, 12),
                # labels = c('0' = '0', '3' = '3', '6' = '6', '9' = '9', '12' = '\u2265 12')
            ) +
            scale_x_discrete(expand = c(0, 0)) +
            scale_y_discrete(expand = c(0, 0)) +
            theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
            labs(x = 'Cells', y = 'Epithelial markers', fill = 'Expression\nlevel Z-score')

        epi_emt_corr_barplot <- ggplot(
            data.table(gene = names(epi_emt_corr), corr = epi_emt_corr),
            aes(x = factor(gene, levels = gene), y = corr)
        ) +
            geom_col(width = 0.8) +
            scale_x_discrete(expand = c(0, 0.5)) +
            theme_test() +
            labs(x = 'Gene', y = 'Correlation with EMT score') +
            coord_flip()

        cat('Done!\n')

        list(
            lineplots = lineplots,
            combined_lineplot = combined_lineplot,
            epi_heatmap = epi_heatmap,
            # epi_emt_corr_lineplot = epi_emt_corr_lineplot,
            epi_emt_corr_barplot = epi_emt_corr_barplot,
            data = plot_data,
            epi_markers = epi_markers,
            epi_emt_corr = epi_emt_corr
        )

    },
    simplify = FALSE,
    USE.NAMES = TRUE
)

# Combining heatmaps:

aligned_plots <- sapply(
	names(emt_epi_comparison),
	function(ct) {

		aligned_1 <- align_plots(
			emt_epi_comparison[[ct]]$epi_heatmap +
				theme(
					axis.title = element_blank(),
					plot.title = element_text(size = 16),
					plot.margin = unit(c(20, 1, 1, 20), 'pt'),
					legend.position = 'none'
				) +
				labs(title = sc_cancer_caf_heatmaps_args_final[[ct]]$annotations_title),
			emt_epi_comparison[[ct]]$epi_emt_corr_barplot +
				scale_y_continuous(
					limits = range(unlist(lapply(emt_epi_comparison, `[[`, 'epi_emt_corr'))),
					labels = c('-0.3' = '-0.3', '-0.2' = '', '-0.1' = '', '0.0' = '0', '0.1' = '', '0.2' = '0.2')
				) +
				theme(
					axis.text.y = element_blank(),
					axis.title.y = element_blank(),
					axis.ticks.y = element_blank(),
					# Setting margin is a bodge to allow us to align everything:
					axis.title.x = element_text(margin = margin(t = 50)),
					plot.margin = unit(c(20, 20, 5.5, 1), 'pt')
				) +
				# ylim(range(unlist(lapply(emt_epi_comparison, `[[`, 'epi_emt_corr')))) +
				labs(y = 'Correlation with\nEMT score'), # It's confusing with coord_flip()...
			align = 'h'
		)

		aligned_2 <- align_plots(
			emt_epi_comparison[[ct]]$epi_heatmap +
				theme(
					axis.title = element_blank(),
					plot.title = element_text(size = 16),
					plot.margin = unit(c(20, 1, 1, 20), 'pt'),
					legend.position = 'none'
				) +
				labs(title = sc_cancer_caf_heatmaps_args_final[[ct]]$annotations_title),
			emt_epi_comparison[[ct]]$combined_lineplot +
				theme(
					axis.title.y = element_text(margin = margin(r = -20)),
					plot.margin = unit(c(1, 1, 5.5, 20), 'pt'),
					legend.position = 'none'
				) +
				labs(
					x = 'Cells'
					# x = switch(
						# (which(names(emt_epi_comparison) == ct) == length(emt_epi_comparison)) + 1,
						# NULL,
						# 'Cells'
					# )
				),
			align = 'v'
		)

		aligned_1[[1]]$heights[8] <- sum(unit(2.75, 'pt'), unit(0, 'cm'))
		aligned_1[[1]]$heights[9] <- unit(0, 'cm')
		aligned_1[[1]]$heights[12] <- unit(1, 'pt')
		aligned_1[[1]]$widths[3] <- aligned_2[[1]]$widths[3]
		aligned_1[[2]]$heights[8] <- sum(unit(2.75, 'pt'), unit(0, 'cm'))
		aligned_1[[2]]$heights[9] <- unit(0, 'cm')
		aligned_1[[2]]$heights[12] <- unit(1, 'pt')

		list(aligned_1[[1]], aligned_1[[2]], aligned_2[[2]])

	},
	simplify = FALSE,
	USE.NAMES = TRUE
)

cairo_pdf('../data_and_figures/final_figures_resubmission/S5.pdf', width = 12, height = 16)
# cairo_pdf('../data_and_figures/sc_epi_vs_emt.pdf', width = 12, height = 16)

plot_grid(
    plot_grid(
        blank_plot(),
        get_legend(
			emt_epi_comparison[[1]]$epi_heatmap +
				theme(legend.direction = 'horizontal', legend.text = element_text(size = 11), legend.key.width = unit(30, 'pt'), legend.justification = c(1, 1)) +
				labs(fill = 'Expression level Z-score\n')
		),
		blank_plot(),
        get_legend(
			emt_epi_comparison[[1]]$combined_lineplot +
				theme(legend.text = element_text(size = 11), legend.key.width = unit(40, 'pt'), legend.justification = c(0, 1))
		),
        blank_plot(),
        nrow = 1,
        ncol = 5,
        rel_widths = c(1, 2, 1, 2, 1)
    ),
    plot_grid(
		plot_grid(
			plotlist = lapply(
				names(emt_epi_comparison)[1:4],
				function(ct) {
					plot_grid(
						plotlist = aligned_plots[[ct]],
						nrow = 2,
						ncol = 2,
						rel_heights = c(3 + length(emt_epi_comparison[[ct]]$epi_markers), 8),
						# rel_heights = c(2 + (length(emt_epi_comparison[[ct]]$epi_markers) - 12)*0.15, 1),
						rel_widths = c(3, 1)
					)
				}
			),
			nrow = 4,
			ncol = 1,
			# align = 'v',
			rel_heights = sapply(emt_epi_comparison[1:4], function(li) {11 + length(li$epi_markers)})
			# rel_heights = sapply(emt_epi_comparison[1:4], function(li) {3 + (length(li$epi_markers) - 12)*0.075})
		),
		plot_grid(
			plotlist = lapply(
				names(emt_epi_comparison)[5:8],
				function(ct) {
					plot_grid(
						plotlist = aligned_plots[[ct]],
						nrow = 2,
						ncol = 2,
						rel_heights = c(3 + length(emt_epi_comparison[[ct]]$epi_markers), 8),
						# rel_heights = c(2 + (length(emt_epi_comparison[[ct]]$epi_markers) - 12)*0.15, 1),
						rel_widths = c(3, 1)
					)
				}
			),
			nrow = 4,
			ncol = 1,
			# align = 'v',
			rel_heights = sapply(emt_epi_comparison[5:8], function(li) {11 + length(li$epi_markers)})
			# rel_heights = sapply(emt_epi_comparison[5:8], function(li) {3 + (length(li$epi_markers) - 12)*0.075})
		),
		nrow = 1,
		ncol = 2
	),
    nrow = 2,
    ncol = 1,
    rel_heights = c(1, 19)
)

dev.off()

# I don't think we need to bother doing this for the lenient classifications, because this shouldn't change the analysis of the cancer cells.





# Lineplots of contributions of different cell types to mesenchymal signal in simulated tumours:

# Use the following to help decide what the maximum mean count should be:
# sapply(sc_metadata, function(x) eval(x$read_quote)[, .N, by = cell_type], simplify = FALSE, USE.NAMES = TRUE)

max_mean_counts <- list(
    brca = 800,
    coadread = 1000,
    hnsc = 100,
    lihc = 100,
	luad = 1000,
    # luad = 400,
	lusc = 200,
    ov = 500,
    paad = 1000
)

# lineplots <- readRDS('../data_and_figures/simulated_bulk_lineplots.rds')

# In the below, I'm setting normalise_fun to NULL, so that I don't divide each gene by its 95th percentile (the default behaviour).  I think I could
# also use normalise_fun to reverse the log before computing the contributions: try setting normalise_fun = function(x) {2^x - 1}

lineplots <- sapply(
    names(sc_metadata),
    function(ct) {
		cat(paste0(ct, '\n'))
        sc_data <- eval(sc_metadata[[ct]]$read_quote)
		# if('cell_type_lenient' %in% names(sc_data)) {sc_data[, cell_type_lenient := NULL]}
		if(!is.null(sc_metadata[[ct]]$rare_cell_types)) {
			sc_data[, cell_type := mapvalues(cell_type, sc_metadata[[ct]]$rare_cell_types, rep('rare', length(sc_metadata[[ct]]$rare_cell_types)))]
		}
        set.seed(sc_cancer_caf_args[[ct]]$seed)
        simulated_tumours_lineplot(
            sc_data,
            genes_list[[ct]]$filtered_cancer_caf,
			initial_types = sc_metadata[[ct]]$initial_cell_types,
			normalise_fun = NULL,
            plot_title = sc_cancer_caf_heatmaps_args[[ct]]$annotations_title,
            max_mean_count = max_mean_counts[[gsub('_lenient', '', ct)]],
			legend_labels = c(
				'b_cell' = 'B cell',
				'cancer' = 'Cancer',
				'endothelial' = 'Endothelial',
				'caf' = 'CAF',
				'macrophage' = 'Macrophage',
				'mast' = 'Mast',
				't_cell' = 'T cell'
			),
			legend_colours = c(
				'b_cell' = '#8DD3C7',
				'cancer' = '#FB8072',
				'endothelial' = '#BC80BD',
				'caf' = '#FDB462',
				'macrophage' = '#80B1D3',
				'mast' = '#FCCDE5',
				't_cell' = '#B3DE69'
			),
        )
    },
    simplify = FALSE,
    USE.NAMES = TRUE
)

saveRDS(lineplots, '../data_and_figures/simulated_bulk_lineplots.rds')

# Save to PDF with each figure on a separate page:

# pdf(
#     '../data_and_figures/simulated_bulk_lineplots.pdf',
#     width = 7,
#     height = 5
# )
# lapply(lineplots, `[[`, 'lineplot')
# dev.off()

# Save as a single combined figure on one page:

dummy_legend_plot <- ggplot(
    data = data.table(
        x = 1:7,
        y = 1,
        f = factor(
            c('b_cell', 'caf', 'cancer', 'endothelial', 'macrophage', 'mast', 't_cell'),
            levels = c('b_cell', 'caf', 'cancer', 'endothelial', 'macrophage', 'mast', 't_cell')
        )
    )
) +
    geom_tile(aes(x = x, y = y, fill = f)) +
    scale_fill_manual(
        labels = c(
            'b_cell' = 'B cell',
            'caf' = 'CAF',
            'cancer' = 'Cancer',
            'endothelial' = 'Endothelial',
            'macrophage' = 'Macrophage',
            'mast' = 'Mast',
            't_cell' = 'T cell'
        ),
        values = c(
            'b_cell' = '#8DD3C7',
            'caf' = '#FDB462',
            'cancer' = '#FB8072',
            'endothelial' = '#BC80BD',
            'macrophage' = '#80B1D3',
			'mast' = '#FCCDE5', # This is better than the yellow (#FFED6F) previously used for mast
            # 'mast' = '#FFED6F',
            't_cell' = '#B3DE69'
        )
    ) +
    labs(fill = 'Cell type') +
    theme(legend.key = element_rect(size = 1, colour = 'white'), legend.key.size = unit(15, 'pt'))

dummy_legend <- get_legend(dummy_legend_plot)

# To see just the legend:

# grid::grid.newpage()
# grid::grid.draw(dummy_legend)

pdf('../data_and_figures/simulated_bulk_lineplots.pdf', width = 12, height = 6)

plot_grid(
    plot_grid(
        plotlist = list(
            brca = lineplots$brca$lineplot +
                theme(axis.text.x = element_blank(), legend.position = 'none', plot.margin = unit(c(5.5, 5.5, 0, 20), 'pt')) +
                labs(x = NULL, y = NULL),
            coadread = lineplots$coadread$lineplot +
                theme(axis.text = element_blank(), legend.position = 'none', plot.margin = unit(c(5.5, 5.5, 0, 5.5), 'pt')) +
                labs(x = NULL, y = NULL),
            hnsc = lineplots$hnsc$lineplot +
                theme(axis.text = element_blank(), legend.position = 'none', plot.margin = unit(c(5.5, 5.5, 0, 5.5), 'pt')) +
                labs(x = NULL, y = NULL),
            lihc = lineplots$lihc$lineplot +
                theme(axis.text = element_blank(), legend.position = 'none', plot.margin = unit(c(5.5, 5.5, 0, 5.5), 'pt')) +
                labs(x = NULL, y = NULL)
        ),
        nrow = 1,
        ncol = 4,
        rel_widths = c(1.17, 1, 1, 1), # Because left plots are squashed by larger left margin
        align = 'h'
    ),
    plot_grid(
        plotlist = list(
            luad = lineplots$luad$lineplot +
                theme(legend.position = 'none', plot.margin = unit(c(5.5, 5.5, 20, 20), 'pt')) +
                labs(x = NULL, y = NULL),
            lusc = lineplots$lusc$lineplot +
                theme(axis.text.y = element_blank(), legend.position = 'none', plot.margin = unit(c(5.5, 5.5, 20, 5.5), 'pt')) +
                labs(x = NULL, y = NULL),
            ov = lineplots$ov$lineplot +
                theme(axis.text.y = element_blank(), legend.position = 'none', plot.margin = unit(c(5.5, 5.5, 20, 5.5), 'pt')) +
                labs(x = NULL, y = NULL),
            paad = lineplots$paad$lineplot +
                theme(axis.text.y = element_blank(), legend.position = 'none', plot.margin = unit(c(5.5, 5.5, 20, 5.5), 'pt')) +
                labs(x = NULL, y = NULL)
        ),
        nrow = 1,
        ncol = 4,
        rel_widths = c(1.17, 1, 1, 1)
    ),
    nrow = 2,
    ncol = 1,
    rel_heights = c(1, 1.1) # Because lower plots are squashed by larger bottom margin
) %>% plot_grid(dummy_legend, nrow = 1, ncol = 2, rel_widths = c(10.5, 1)) +
    draw_label('Proportion of tumour', x = 0.47, y = 0, vjust = -0.5, size = 12) +
    draw_label('Proportion of gene expression', x = 0, y = 0.5, vjust = 1.3, angle = 90, size = 12)

dev.off()

# For the lenient classifications:

pdf('../data_and_figures/simulated_bulk_lineplots_lenient.pdf', width = 9.5, height = 6)

plot_grid(
    plot_grid(
        plotlist = list(
            brca = lineplots$brca_lenient$lineplot +
                theme(axis.text.x = element_blank(), legend.position = 'none', plot.margin = unit(c(5.5, 5.5, 0, 20), 'pt')) +
                labs(x = NULL, y = NULL),
            coadread = lineplots$coadread_lenient$lineplot +
                theme(axis.text = element_blank(), legend.position = 'none', plot.margin = unit(c(5.5, 5.5, 0, 5.5), 'pt')) +
                labs(x = NULL, y = NULL),
            lihc = lineplots$lihc_lenient$lineplot +
                theme(axis.text.y = element_blank(), legend.position = 'none', plot.margin = unit(c(5.5, 5.5, 0, 5.5), 'pt')) +
                labs(x = NULL, y = NULL)
        ),
        nrow = 1,
        ncol = 3,
        rel_widths = c(1.15, 1, 1, 1), # Because left plots are squashed by larger left margin
        align = 'h'
    ),
    plot_grid(
        plotlist = list(
            # luad = lineplots$luad_lenient$lineplot +
                # theme(legend.position = 'none', plot.margin = unit(c(5.5, 5.5, 20, 20), 'pt')) +
                # labs(x = NULL, y = NULL),
            lusc = lineplots$lusc_lenient$lineplot +
                theme(
                    # axis.text.y = element_blank(),
                    legend.position = 'none',
					plot.margin = unit(c(5.5, 5.5, 20, 20), 'pt')
                    # plot.margin = unit(c(5.5, 5.5, 20, 5.5), 'pt')
                ) +
                labs(x = NULL, y = NULL),
            ov = lineplots$ov_lenient$lineplot +
                theme(axis.text.y = element_blank(), legend.position = 'none', plot.margin = unit(c(5.5, 5.5, 20, 5.5), 'pt')) +
                labs(x = NULL, y = NULL)
        ),
        nrow = 1,
        ncol = 3,
        rel_widths = c(1.15, 1, 1)
    ),
    nrow = 2,
    ncol = 1,
    rel_heights = c(1, 1.1) # Because lower plots are squashed by larger bottom margin
) %>% plot_grid(dummy_legend, nrow = 1, ncol = 2, rel_widths = c(8, 1)) +
    draw_label('Proportion of tumour', x = 0.47, y = 0, vjust = -0.5, size = 12) +
    draw_label('Proportion of gene expression', x = 0, y = 0.5, vjust = 1.3, angle = 90, size = 12)

dev.off()





# Deconvolution for simulated bulk profiles:

simulated_bulk_data <- sapply(
    names(sc_metadata),
    function(ct) {

        cat(paste0(ct, '\n'))
        sc_data <- eval(sc_metadata[[ct]]$read_quote)
		# if('cell_type_lenient' %in% names(sc_data)) {sc_data[, cell_type_lenient := NULL]}
		if(!is.null(sc_metadata[[ct]]$rare_cell_types)) { # Check if this is necessary here
			sc_data[, cell_type := mapvalues(cell_type, sc_metadata[[ct]]$rare_cell_types, rep('rare', length(sc_metadata[[ct]]$rare_cell_types)))]
		}
        set.seed(sc_cancer_caf_args[[ct]]$seed)
        simulated_tumours_data(
            as.matrix(sc_data[, -c('id', 'patient', 'cell_type')]),
            types = sc_data$cell_type,
            id_prefix = ct,
            max_mean_count = sc_data[cell_type == 'cancer', .N, by = patient][, round(quantile(N, 0.9))]
            # max_mean_count = max_mean_counts[[ct]]
            # genes = genes_list[[ct]]
        )
    },
    simplify = FALSE,
    USE.NAMES = TRUE
)

# Change the following to convert the expression data from a matrix into a data table.

simulated_bulk_data <- sapply(
    c('expression_data', 'meta_data'),
    function(x) {
        rbindlist(
            lapply(
                simulated_bulk_data,
                function(y) as.data.table(y[[x]], keep.rownames = 'id')
            ),
            fill = TRUE
        )
    },
    simplify = FALSE,
    USE.NAMES = TRUE
)

# Replace NAs arising from fill = TRUE to zeros:

simulated_bulk_data$expression_data[is.na(simulated_bulk_data$expression_data)] <- 0

# simulated_bulk_data$meta_data[, cancer_type := str_extract(id, '^[a-z]*')]
simulated_bulk_data$meta_data[, cancer_type := gsub('[0-9]+', '', id)]

setcolorder(simulated_bulk_data$meta_data, c('id', 'cancer_type'))

fwrite(simulated_bulk_data$expression_data, '../data_and_figures/simulated_bulk_data.csv')
fwrite(simulated_bulk_data$meta_data, '../data_and_figures/simulated_bulk_metadata.csv')





# Start from saved simulated bulk data:

simulated_bulk_data <- fread('../data_and_figures/simulated_bulk_data.csv', key = 'id')
simulated_bulk_metadata <- fread('../data_and_figures/simulated_bulk_metadata.csv', key = 'id')

ccle <- fread('../../CCLE_cancer_type_Av.csv', key = 'gene_id')
cell_type_markers <- fread('../../cell_type_markers.csv')
# extra_data <- fread('../data_and_figures/collated_extra_data.csv', key = 'gene')

# The name 'tcga_cancer_type' in this list is confusing - it's really the elements of simulated_bulk_data$meta_data$cancer_type (hence it's lower case).
# I guess I need to call it this for the deconv function.

deconv_args_per_ct <- list(
    brca = list(tcga_cancer_types = 'brca', ccle_cancer_type = 'breast', seed = 1087,
        plot_title = 'Breast', genes_filter_fun = function(x) 1:200),
	brca_lenient = list(tcga_cancer_types = 'brca_lenient', ccle_cancer_type = 'breast', seed = 3510,
		plot_title = 'Breast', genes_filter_fun = function(x) 1:200),
    coadread = list(tcga_cancer_types = 'coadread', ccle_cancer_type = 'large intestine', seed = 6185,
        plot_title = 'Colorectal', genes_filter_fun = function(x) 1:200),
	coadread_lenient = list(tcga_cancer_types = 'coadread_lenient', ccle_cancer_type = 'large intestine', seed = 5106,
        plot_title = 'Colorectal', genes_filter_fun = function(x) 1:200),
    hnsc = list(tcga_cancer_types = 'hnsc', ccle_cancer_type = 'HNSCC', seed = 6367,
        plot_title = 'Head and Neck', genes_filter_fun = function(x) 1:200),# extra_data_source = 'puram_hnscc_2017
    lihc = list(tcga_cancer_types = 'lihc', ccle_cancer_type = 'liver', seed = 5199,
        plot_title = 'Liver', genes_filter_fun = function(x) 1:200),# extra_data_source = 'ma_liver_2019'
	lihc_lenient = list(tcga_cancer_types = 'lihc_lenient', ccle_cancer_type = 'liver', seed = 1682,
        plot_title = 'Liver', genes_filter_fun = function(x) 1:200),
    luad = list(tcga_cancer_types = 'luad', ccle_cancer_type = 'lung', seed = 5395,
        plot_title = 'Lung Adenocarcinoma', genes_filter_fun = function(x) 1:200),
	# luad_lenient = list(tcga_cancer_types = 'luad_lenient', ccle_cancer_type = 'lung', seed = 5561,
        # plot_title = 'Lung Adenocarcinoma', genes_filter_fun = function(x) 1:200),
    lusc = list(tcga_cancer_types = 'lusc', ccle_cancer_type = 'lung', seed = 6212,
        plot_title = 'Lung squamous', genes_filter_fun = function(x) 1:200),
    lusc_lenient = list(tcga_cancer_types = 'lusc_lenient', ccle_cancer_type = 'lung', seed = 5596,
        plot_title = 'Lung squamous', genes_filter_fun = function(x) 1:200),
    ov = list(tcga_cancer_types = 'ov', ccle_cancer_type = 'ovary', seed = 3854,
        plot_title = 'Ovarian', genes_filter_fun = function(x) 1:200),
    ov_lenient = list(tcga_cancer_types = 'ov_lenient', ccle_cancer_type = 'ovary', seed = 8619,
        plot_title = 'Ovarian', genes_filter_fun = function(x) 1:200),
    paad = list(tcga_cancer_types = 'paad', ccle_cancer_type = 'pancreas', seed = 303,
        plot_title = 'Pancreatic', genes_filter_fun = function(x) 1:200)
)

# Defining gene list from the simulated bulk data:

# Note I tried it with the genes we got from the single cell data and the result was
# worse, especially for colorectal, where the pattern disappears.  If you do want to
# use the genes from single cell data, delete initial_gene_weights = FALSE and put
# genes_filter_fun = NULL and genes_from_tcga_fun = NULL.

# simulated_deconvs <- readRDS('../data_and_figures/simulated_deconvs.rds')

simulated_deconvs <- sapply(
    names(deconv_args_per_ct),
    function(ct) {
        cat(paste0(ct, '\n'))
        do.call(
            deconvolve_emt_caf_data,
            args = c(
                list(
                    expression_data = simulated_bulk_data,
                    meta_data = simulated_bulk_metadata,
                    # extra_data = extra_data,
                    genes = emt_markers,
                    cell_type_markers = cell_type_markers,
                    ccle_data = ccle,
                    initial_gene_weights = FALSE
                ),
                deconv_args_per_ct[[ct]][names(deconv_args_per_ct[[ct]]) != 'plot_title']
            )
        )
    },
    simplify = FALSE,
    USE.NAMES = TRUE
)

# Reorder genes in deconv results:
simulated_deconvs <- sapply(simulated_deconvs, deconv_reorder, simplify = FALSE, USE.NAMES = TRUE)

# Add data for single cell comparison colour bar:
sc_comp_diff <- sapply(
	names(deconv_args_per_ct),
	function(ct) {

		cat(paste0(ct, '\n'))

		sc_data <- eval(sc_metadata[[ct]]$read_quote)[, c('id', 'cell_type', simulated_deconvs[[ct]]$genes_filtered), with = FALSE]
		sc_data <- melt(sc_data, id.vars = c('id', 'cell_type'), variable.name = 'gene', value.name = 'expression_level')

		# Centre genes:
		gene_averages <- sc_data[
			,
			.(ave_exp = mean(c(mean(expression_level[cell_type == 'cancer']), mean(expression_level[cell_type == 'caf'])))),
			by = .(symbol = gene)
		]
		sc_data[, expression_level := expression_level - gene_averages[symbol == unique(gene), ave_exp], by = gene]

		# Centre cells:
		sc_data[, expression_level := expression_level - mean(expression_level), by = id]

		sc_data[, .(d = mean(expression_level[cell_type == 'cancer']) - mean(expression_level[cell_type == 'caf'])), by = gene][, setNames(d, gene)]

	},
	simplify = FALSE,
	USE.NAMES = TRUE
)

for(ct in names(deconv_args_per_ct)) {simulated_deconvs[[ct]]$extra_data_score <- sc_comp_diff[[ct]]}

saveRDS(simulated_deconvs, '../data_and_figures/simulated_deconvs.rds')

deconv_plot_args_per_ct <- list(
    brca = list(
        heatmap_annotations = c('COL1A2', 'COL5A2', 'ECM1', 'FN1', 'NOTCH2', 'RHOB', 'SNAI1', 'SNAI2', 'TWIST1', 'VCAN', 'VEGFA', 'VIM', 'ZEB1', 'ZEB2'),
        plot_title = 'Breast'
    ),
	brca_lenient = list(
        heatmap_annotations = c('COL1A2', 'COL5A2', 'ECM1', 'FN1', 'NOTCH2', 'RHOB', 'SNAI1', 'SNAI2', 'TWIST1', 'VCAN', 'VEGFA', 'VIM', 'ZEB1', 'ZEB2'),
        plot_title = 'Breast'
    ),
    coadread = list(
        heatmap_annotations = c('AREG', 'COL1A1', 'COL3A1', 'CXCL1', 'DCN', 'FN1', 'GADD45A', 'PLOD3', 'SNAI1', 'SNAI2', 'TWIST1', 'VIM', 'ZEB1', 'ZEB2'),
        plot_title = 'Colorectal'
    ),
	coadread_lenient = list(
        heatmap_annotations = c('AREG', 'COL1A1', 'COL3A1', 'CXCL1', 'DCN', 'FN1', 'GADD45A', 'PLOD3', 'SNAI1', 'SNAI2', 'TWIST1', 'VIM', 'ZEB1', 'ZEB2'),
        plot_title = 'Colorectal'
    ),
    hnsc = list(
        heatmap_annotations = c('ACTA2', 'COL1A2', 'COL3A1', 'LAMC2', 'PCOLCE', 'SDC1', 'SNAI1', 'SNAI2', 'TGFBI', 'TNC', 'TWIST1', 'VIM', 'ZEB1', 'ZEB2'),
        plot_title = 'Head and Neck'
    ),
    lihc = list(
        heatmap_annotations = c('ACTA2', 'BGN', 'COL1A2', 'LAMC2', 'QSOX1', 'SNAI1', 'SNAI2', 'THY1', 'TNFRSF12A', 'TWIST1', 'VEGFA', 'VIM', 'ZEB1', 'ZEB2'),
        plot_title = 'Liver'
    ),
	lihc_lenient = list(
        heatmap_annotations = c('ACTA2', 'BGN', 'COL1A2', 'LAMC2', 'QSOX1', 'SNAI1', 'SNAI2', 'THY1', 'TNFRSF12A', 'TWIST1', 'VEGFA', 'VIM', 'ZEB1', 'ZEB2'),
        plot_title = 'Liver'
    ),
    luad = list(
        heatmap_annotations = c('CALU', 'COL1A2', 'COL3A1', 'MMP2', 'QSOX1', 'SDC4', 'SNAI1', 'SNAI2', 'THY1', 'TWIST1', 'VEGFA', 'VIM', 'ZEB1', 'ZEB2'),
        plot_title = 'Lung Adenocarcinoma'
    ),
	# luad_lenient = list(
        # heatmap_annotations = c('ACTA2', 'COL1A1', 'COL1A2', 'ITGB1', 'LAMC2', 'QSOX1', 'RHOB', 'SNAI1', 'SNAI2', 'THY1', 'TWIST1', 'VIM', 'ZEB1', 'ZEB2'),
        # plot_title = 'Lung Adenocarcinoma'
    # ),
    lusc = list( # Need to change annotations from here downwards...
        heatmap_annotations = c('COL1A2', 'COL3A1', 'FAP', 'IGFBP2', 'PFN2', 'SDC1', 'SNAI1', 'SNAI2', 'THY1', 'TNC', 'TWIST1', 'VIM', 'ZEB1', 'ZEB2'),
        plot_title = 'Lung squamous'
    ),
	lusc_lenient = list(
        heatmap_annotations = c('COL1A2', 'COL3A1', 'FAP', 'IGFBP2', 'PFN2', 'SDC1', 'SNAI1', 'SNAI2', 'THY1', 'TNC', 'TWIST1', 'VIM', 'ZEB1', 'ZEB2'),
        plot_title = 'Lung squamous'
    ),
    ov = list(
        heatmap_annotations = c('ACTA2', 'CDH2', 'COL3A1', 'FAP', 'LAMC1', 'MMP1', 'PFN2', 'SNAI1', 'SNAI2', 'THY1', 'TWIST1', 'VIM', 'ZEB1', 'ZEB2'),
        plot_title = 'Ovarian'
    ),
	ov_lenient = list(
        heatmap_annotations = c('ACTA2', 'CDH2', 'COL3A1', 'FAP', 'LAMC1', 'MMP1', 'PFN2', 'SNAI1', 'SNAI2', 'THY1', 'TWIST1', 'VIM', 'ZEB1', 'ZEB2'),
        plot_title = 'Ovarian'
    ),
    paad = list(
        heatmap_annotations = c('COL1A2', 'COL3A1', 'DCN', 'FAP', 'LAMC2', 'QSOX1', 'SDC4', 'SNAI1', 'SNAI2', 'TWIST1', 'VEGFA', 'VIM', 'ZEB1', 'ZEB2'),
        plot_title = 'Pancreatic'
    )
)

# In the following, I am not opting to include epithelial cells in the diagnostics, mainly to save
# time (it saves a lot of time), but also because the epithelial cell correlations are (thankfully)
# as I would expect, and because conceptually I don't think it adds much.

simulated_deconv_plots <- sapply(
    names(deconv_args_per_ct),
    function(ct) {
        cat(paste0(ct, '\n'))
        do.call(
            deconvolve_emt_caf_plots,
            args = c(
                list(
                    data = simulated_deconvs[[ct]],
                    # Include the following only if you want epithelial scores (takes much longer):
                    # expression_data = simulated_bulk_data,
                    heatmap_legend_title = 'Correlation',
                    heatmap_colours = rev(colorRampPalette(brewer.pal(11, "RdBu"))(50)),
                    heatmap_colour_limits = c(-1, 1),
                    heatmap_legend_breaks = c(-1, -0.5, 0, 0.5, 1),
					heatmap_legend_justification = 'left',
                    heatmap_annotations_nudge = 0.3,
                    purity_colours = rev(colorRampPalette(brewer.pal(11, "PuOr"))(50)),
                    purity_colour_limits = c(-1, 1),
                    purity_legend_breaks = c(-1, 0, 1),
                    purity_legend_title = 'Correlation with purity\n',
                    purity_legend_direction = 'horizontal',
                    purity_axis_title = NULL,
                    ccle_colours = rev(colorRampPalette(brewer.pal(11, 'PiYG')[c(1:3, 6, 9:11)])(50)),
					ccle_fun = function(x) {caTools::runmean(x - mean(x), 30)/max(abs(caTools::runmean(x - mean(x), 30)))},
                    ccle_legend_breaks = c(-1, 0, 1),
                    ccle_legend_title = 'Tumours vs. cell lines\n',
                    ccle_legend_direction = 'horizontal',
                    ccle_axis_title = NULL,
                    extra_colours = colorRampPalette(c('turquoise4', 'turquoise', 'azure', 'gold2', 'gold4'))(50),
					extra_fun = function(x) caTools::runmean(x, 30),
					extra_colour_limits = c(-4, 4),
					extra_legend_breaks = c(-4, 0, 4),
                    extra_legend_title = 'scRNA-seq: CAF vs. cancer\n',
                    extra_legend_direction = 'horizontal',
                    extra_axis_title = NULL,
                    bar_legend_justification = 'left'
                    # bar_legend_width = NULL,
                    # bar_legend_height = NULL
                ),
                deconv_plot_args_per_ct[[ct]]
            )
        )
    },
    simplify = FALSE,
    USE.NAMES = TRUE
)

cairo_pdf('../data_and_figures/simulated_deconv_figures.pdf', width = 7, height = 6, onefile = TRUE)

for(ct in names(simulated_deconv_plots)[!grepl('lenient', names(simulated_deconv_plots))]) {
	deconv_plot(
		simulated_deconv_plots[ct],
		legends_rel_size = c(0.5, 0.75, 0, 0.75, 0, 0.75, 0.75, 3, 2),
		legends_space = 0.7,
		title.position = 'right'
	) %>% print
}

# deconv_plot(
    # simulated_deconv_plots[!grepl('lenient', names(simulated_deconv_plots))],
    # n_row = 2,
    # n_col = 4,
    # left_plot_width = 1.06,
    # legends_arrange = 'vertical',
    # legends_rel_size = c(1, 1, 0, 1, 0, 1, 0, 3, 2),
    # title.position = 'right'
# )

dev.off()

cairo_pdf('../data_and_figures/simulated_deconv_figures_lenient.pdf', width = 7, height = 6, onefile = TRUE)

for(ct in names(simulated_deconv_plots)[grepl('lenient', names(simulated_deconv_plots))]) {
	deconv_plot(
		simulated_deconv_plots[ct],
		legends_rel_size = c(0.5, 0.75, 0, 0.75, 0, 0.75, 0.75, 3, 2),
		legends_space = 0.7,
		title.position = 'right'
	) %>% print
}

# deconv_plot(
    # simulated_deconv_plots[grepl('lenient', names(simulated_deconv_plots))],
    # n_row = 2,
    # n_col = 4,
    # left_plot_width = 1.06,
    # legends_arrange = 'vertical',
    # legends_rel_size = c(1, 1, 0, 1, 0, 1, 0, 3, 2),
    # title.position = 'right'
# )

dev.off()

# Diagnostic figures:

# pdf(
    # '../data_and_figures/simulated_deconv_diagnostics.pdf',
    # width = 10,
    # height = 16.8
# )

# i <- 1

# for(deconv_ct in simulated_deconv_plots) {

    # ggarrange(
        # plots = c(
            # list(deconv_ct$plots$heatmap),
            # deconv_ct$diagnostics$alternative_purity_cor,
            # deconv_ct$diagnostics$cell_type_bars
        # ),
        # ncol = 1,
        # nrow = 11,
        # heights = c(8, rep(0.8, 10)),
        # newpage = switch((i == 1) + 1, TRUE, FALSE)
    # )

    # i <- i + 1

# }

# dev.off()

# This is a bit worrying because 4 of the 7 analyses fail the criterion I used for
# filtering the TCGA deconvs (regression slope < -0.1):

# diagnostic_cell_type_lms <- sapply(
#     simulated_deconvs,
#     function(deconv_ct) {
#         with(
#             deconv_ct,
#             cor_with_initial_and_cell_types[
#                 genes_filtered
#             ][
#                 ordering,
#                 sapply(
#                     .SD,
#                     function(ct) setNames(lm(ct ~ I(.I/.N))$coeff['I(.I/.N)'], NULL),
#                     USE.NAMES = TRUE
#                 ),
#                 .SDcols = cell_types
#             ]
#         )
#     },
#     USE.NAMES = TRUE
# )





# Heatmaps showing expression of genes in single cells data, ordered as in the deconv results:

# Can't use max_mean_counts because that means too few cells in some cases.
sample_sizes <- list(
    brca = 1000,
	brca_lenient = 1000,
    coadread = 800,
	coadread_lenient = 800,
    hnsc = 400,
    lihc = 150,
	lihc_lenient = 200,
    luad = 1000,
	# luad = 200,
	# luad_lenient = 300,
	lusc = 100,
	lusc_lenient = 100,
    ov = 1500,
	ov_lenient = 2000,
    paad = 1500
)

set.seed(4508) # Is it enough to set a single seed before the whole loop?

sc_sim_deconv_comp <- sapply(
	names(simulated_deconvs),
	function(ct) {

		cat(paste0(ct, '\n'))

        sc_data <- eval(sc_metadata[[ct]]$read_quote)

		sc_deconv_comp <- sapply(
			c('cancer', 'caf'),
			function(x) {
				plot_data <- copy(sc_data[cell_type == x])[, complexity := apply(.SD, 1, function(x) sum(x > 0)), .SDcols = -c('id', 'patient', 'cell_type')][
					,
					.SD[sample(1:.N, sample_sizes[[ct]], prob = complexity)]
				][, complexity := NULL][, c('id', simulated_deconvs[[ct]]$genes_filtered), with = FALSE]
				plot_data <- melt(plot_data, id.vars = 'id', variable.name = 'gene', value.name = 'expression_level')
				ordered_cell_ids <- plot_data[
					gene %in% with(
						simulated_deconvs[[ct]],
						do.call(c('head', 'tail')[which(c('cancer', 'caf') == x)], list(genes_filtered[ordering], 20))
					),
					.(top_20_mean = mean(expression_level)),
					by = id
				][order(top_20_mean), id]
				list(plot_data = plot_data, ordered_cell_ids = ordered_cell_ids)
			},
			simplify = FALSE,
			USE.NAMES = TRUE
		)

		sc_deconv_comp_data <- rbind(
			sc_deconv_comp$cancer$plot_data[, cell_type := 'cancer'],
			sc_deconv_comp$caf$plot_data[, cell_type := 'caf']
		)

		# Calculate averages in cancer cells and CAFs separately, so we can check them afterwards:
		gene_averages_cancer_caf <- sc_deconv_comp_data[
			,
			.(ave_exp = mean(expression_level)),
			by = .(gene, cell_type)
		]

		# To centre genes w.r.t. the average of the averages of cancer and CAF:
		gene_averages <- sc_deconv_comp_data[
			,
			.(ave_exp = mean(c(mean(expression_level[cell_type == 'cancer']), mean(expression_level[cell_type == 'caf'])))),
			by = .(symbol = gene)
		]
		sc_deconv_comp_data[, expression_level := expression_level - gene_averages[symbol == unique(gene), ave_exp], by = gene]

		# To centre the cells as well:
		sc_deconv_comp_data[, expression_level := expression_level - mean(expression_level), by = id]

		sc_heatmaps <- sapply(
			c('cancer', 'caf'),
			function(x) {
				ggplot(
					sc_deconv_comp_data[cell_type == x],
					aes(
						x = factor(gene, levels = with(simulated_deconvs[[ct]], genes_filtered[ordering])),
						y = factor(id, levels = sc_deconv_comp[[x]]$ordered_cell_ids),
						fill = expression_level
					)
				) +
					geom_raster() +
					scale_x_discrete(expand = c(0, 0)) +
					scale_y_discrete(expand = c(0, 0)) +
					scale_fill_gradientn(
						colours = c(
							colorspace::sequential_hcl(25, h = 190, c = 70, l = c(70, 99), power = 0.8),
							colorspace::sequential_hcl(25, h = 60, c = 100, l = c(80, 99), power = 0.8, rev = TRUE)
						),
						# colours = rev(colorRampPalette(brewer.pal(11, "RdBu"))(50)),
						limits = c(-6, 6),
						oob = scales::squish,
						breaks = c(-6, -3, 0, 3, 6),
						labels = c('-6' = '\u2264 -6', '-3' = '-3', '0' = '0', '3' = '3', '6' = '\u2265 6')
					) +
					theme(
						axis.text = element_blank(),
						axis.title.x = element_blank(),
						axis.ticks = element_blank(),
						axis.ticks.length = unit(0, 'pt'),
						plot.margin = unit(c(5.5, 5.5, 1.25, 5.5), 'pt')
						# plot.margin = switch((x == 'cancer') + 1, unit(c(1.25, 5.5, 5.5, 5.5), 'pt'), unit(c(5.5, 5.5, 1.25, 5.5), 'pt'))
					) +
					labs(
						y = mapvalues(x, c('cancer', 'caf'), c('Cancer', 'CAF'), warn_missing = FALSE),
						fill = 'Relative\nexpression\nlevel'
					)
			},
			simplify = FALSE,
			USE.NAMES = TRUE
		)

		# Make versions where we filter out genes with low expression in the single cell data:

		simulated_deconv_filtered <- simulated_deconvs[[ct]]

		filtered_genes <- sc_data[
			,
			sapply(.SD[cell_type == 'cancer'], sc_cancer_caf_args[[ct]]$genes_filter_fun) |
				sapply(.SD[cell_type == 'caf'], sc_cancer_caf_args[[ct]]$genes_filter_fun),
			.SDcols = simulated_deconvs[[ct]]$genes_filtered
		]
		filtered_genes <- names(filtered_genes)[filtered_genes]

		# filtered_genes <- gene_averages_cancer_caf[
			# ,
			# .(pass = ave_exp[cell_type == 'cancer'] > 0.25 | ave_exp[cell_type == 'caf'] > 0.25),
			# by = gene
		# ][pass == TRUE, as.character(gene)]

		ordered_filtered_genes <- with(simulated_deconv_filtered, genes_filtered[ordering][genes_filtered[ordering] %in% filtered_genes])

		simulated_deconv_filtered$ordering <- order(filtered_genes)[order(order(ordered_filtered_genes))]
		simulated_deconv_filtered$genes_filtered <- filtered_genes
		simulated_deconv_filtered$cor_mat <- simulated_deconv_filtered$cor_mat[filtered_genes, filtered_genes]
		simulated_deconv_filtered$cor_with_purity <- sapply(
			simulated_deconv_filtered$cor_with_purity,
			function(x) x[filtered_genes],
			simplify = FALSE,
			USE.NAMES = TRUE
		)
		simulated_deconv_filtered$ccle_comp_diff <- simulated_deconv_filtered$ccle_comp_diff[filtered_genes]

		simulated_deconv_filtered_plots <- do.call(
			deconvolve_emt_caf_plots,
			args = c(
				list(
					data = simulated_deconv_filtered,
					heatmap_legend_title = 'Correlation',
					heatmap_colours = rev(colorRampPalette(brewer.pal(11, "RdBu"))(50)),
					heatmap_colour_limits = c(-1, 1),
					heatmap_legend_breaks = c(-1, -0.5, 0, 0.5, 1),
					heatmap_annotations_side = 'left',
					heatmap_annotations_nudge = 0.3,
					purity_colours = rev(colorRampPalette(brewer.pal(11, "PuOr"))(50)),
					purity_colour_limits = c(-1, 1),
					purity_legend_breaks = c(-1, 0, 1),
					purity_legend_title = 'Correlation with purity\n',
					purity_legend_direction = 'horizontal',
					purity_axis_title = NULL,
					ccle_colours = rev(colorRampPalette(brewer.pal(11, 'PiYG')[c(1:3, 6, 9:11)])(50)),
					ccle_fun = function(x) {caTools::runmean(x - mean(x), 30)/max(abs(caTools::runmean(x - mean(x), 30)))},
					ccle_legend_breaks = c(-1, 0, 1),
					ccle_legend_title = 'Tumours vs. cell lines\n',
					ccle_legend_direction = 'horizontal',
					ccle_axis_title = NULL,
					bar_legend_justification = 'left',
					bar_legend_width = unit(10, 'pt'),
					bar_legend_height = unit(10, 'pt')
				),
				deconv_plot_args_per_ct[[ct]]
			)
		)

		plot_data <- sc_deconv_comp_data[gene %in% simulated_deconv_filtered$genes_filtered]

		plot_data <- sapply(
			c('cancer', 'caf'),
			function(x) {
				plot_data[
					cell_type == x,
					.(
						id = factor(id, levels = sc_deconv_comp[[x]]$ordered_cell_ids),
						gene = factor(gene, levels = with(simulated_deconv_filtered, genes_filtered[ordering])),
						expression_level = expression_level
					)
				]
			},
			simplify = FALSE,
			USE.NAMES = TRUE
		)

		sc_heatmaps_filtered <- sapply(
			c('cancer', 'caf'),
			function(x) {
				ggplot(
					plot_data[[x]],
					# plot_data[cell_type == x],
					aes(
						x = gene,
						y = id,
						# x = factor(gene, levels = with(simulated_deconv_filtered, genes_filtered[ordering])),
						# y = factor(id, levels = sc_deconv_comp[[x]]$ordered_cell_ids),
						fill = expression_level
					)
				) +
					geom_raster() +
					scale_x_discrete(expand = c(0, 0)) +
					scale_y_discrete(expand = c(0, 0)) +
					scale_fill_gradientn(
						colours = c(
							sequential_hcl(25, h = 190, c = 70, l = c(70, 99), power = 0.8),
							sequential_hcl(25, h = 60, c = 100, l = c(80, 99), power = 0.8, rev = TRUE)
						),
						# colours = rev(colorRampPalette(brewer.pal(11, "RdBu"))(50)),
						limits = c(-6, 6),
						oob = squish,
						breaks = c(-6, -3, 0, 3, 6),
						labels = c('-6' = '\u2264 -6', '-3' = '-3', '0' = '0', '3' = '3', '6' = '\u2265 6')
					) +
					theme(
						axis.text = element_blank(),
						axis.title.x = element_blank(),
						axis.ticks = element_blank(),
						axis.ticks.length = unit(0, 'pt'),
						plot.margin = unit(c(5.5, 5.5, 1.25, 0), 'pt')
					) +
					labs(y = mapvalues(x, c('cancer', 'caf'), c('Cancer', 'CAF'), warn_missing = FALSE), fill = 'Relative\nexpression\nlevel')
			},
			simplify = FALSE,
			USE.NAMES = TRUE
		)

		list(
			sc_deconv_comp_data = sc_deconv_comp_data,
			filtered_deconv_data = simulated_deconv_filtered,
			filtered_deconv_figures = simulated_deconv_filtered_plots$plots,
			sc_heatmaps_unfiltered = sc_heatmaps,
			sc_heatmaps_filtered = sc_heatmaps_filtered,
			sc_heatmaps_filtered_data = plot_data,
			gene_averages_cancer_caf = gene_averages_cancer_caf,
			gene_averages_for_centring = gene_averages,
			ordered_cell_ids <- sapply(sc_deconv_comp, `[[`, 'ordered_cell_ids', simplify = FALSE, USE.NAMES = TRUE)
		)

	},
	simplify = FALSE,
	USE.NAMES = TRUE
)

cairo_pdf('../data_and_figures/sc_sim_deconv_comp_centred.pdf', width = 6, height = 9, onefile = TRUE)

for(ct in names(sc_sim_deconv_comp)[!grepl('lenient', names(sc_sim_deconv_comp))]) {

	sim_deconv_figures <- sapply(
		simulated_deconv_plots[[ct]]$plots[c('purity_bar', 'ccle_bar', 'heatmap')],
		function(x) x + theme(legend.position = 'none', plot.title = element_blank()),
		simplify = FALSE,
		USE.NAMES = TRUE
	)
	sim_deconv_figures$heatmap <- sim_deconv_figures$heatmap + theme(plot.margin = unit(c(0, 5.5, 1.25, 5.5), 'pt'))

	plot_grid(
		plot_grid(
			plotlist = c(
				list(blank_plot() + theme(plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), 'pt')) + labs(title = deconv_plot_args_per_ct[[ct]]$plot_title)),
				sim_deconv_figures,
				# lapply(
					# simulated_deconv_plots[[ct]]$plots[c('purity_bar', 'ccle_bar', 'heatmap')],
					# function(x) x + theme(legend.position = 'none', plot.title = element_blank())
				# ),
				lapply(sc_sim_deconv_comp[[ct]]$sc_heatmaps_unfiltered, function(x) x + theme(legend.position = 'none'))
			),
			nrow = 6,
			ncol = 1,
			align = 'v',
			rel_heights = c(1, 1, 1, 15, 5, 5)
		),
		plot_grid(
			plotlist = c(
				list(
					blank_plot(),
					get_legend(
						simulated_deconv_plots[[1]]$plots$purity_bar +
							theme(legend.justification = c(0, 1), legend.direction = 'vertical', legend.key.width = unit(10, 'pt'), legend.key.height = unit(10, 'pt')) +
							labs(fill = 'Correlation\nwith purity')
					),
					get_legend(
						simulated_deconv_plots[[1]]$plots$ccle_bar +
							theme(legend.justification = c(0, 1), legend.direction = 'vertical', legend.key.width = unit(10, 'pt'), legend.key.height = unit(10, 'pt')) +
							labs(fill = 'Tumours vs.\ncell lines')
					),
					get_legend(
						simulated_deconv_plots[[1]]$plots$heatmap + theme(legend.justification = c(0, 1))
					),
					get_legend(sc_sim_deconv_comp[[1]]$sc_heatmaps_unfiltered[[1]] + theme(legend.justification = c(0, 1)))
				)
			),
			nrow = 5,
			ncol = 1,
			rel_heights = c(1, 5, 5, 7, 10)
		),
		nrow = 1,
		ncol = 2,
		rel_widths = c(5, 1)
	) %>% print

}

dev.off()

cairo_pdf('../data_and_figures/sc_sim_deconv_comp_lenient_centred.pdf', width = 6, height = 10, onefile = TRUE)

for(ct in names(sc_sim_deconv_comp)[grepl('lenient', names(sc_sim_deconv_comp))]) {

	sim_deconv_lenient_figures <- sapply(
		simulated_deconv_plots[[ct]]$plots[c('purity_bar', 'ccle_bar', 'heatmap')],
		function(x) x + theme(legend.position = 'none', plot.title = element_blank()),
		simplify = FALSE,
		USE.NAMES = TRUE
	)
	sim_deconv_lenient_figures$heatmap <- sim_deconv_lenient_figures$heatmap + theme(plot.margin = unit(c(0, 5.5, 1.25, 5.5), 'pt'))

	plot_grid(
		plot_grid(
			plotlist = c(
				list(blank_plot() + theme(plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), 'pt')) + labs(title = deconv_plot_args_per_ct[[ct]]$plot_title)),
				sim_deconv_lenient_figures,
				# lapply(
					# simulated_deconv_plots[[ct]]$plots[c('purity_bar', 'ccle_bar', 'heatmap')],
					# function(x) x + theme(legend.position = 'none', plot.title = element_blank())
				# ),
				lapply(sc_sim_deconv_comp[[ct]]$sc_heatmaps_unfiltered, function(x) x + theme(legend.position = 'none'))
			),
			nrow = 6,
			ncol = 1,
			align = 'v',
			rel_heights = c(1, 1, 1, 15, 5, 5)
		),
		plot_grid(
			plotlist = c(
				list(
					blank_plot(),
					get_legend(
						simulated_deconv_plots[[1]]$plots$purity_bar +
							theme(legend.justification = c(0, 1), legend.direction = 'vertical', legend.key.width = unit(10, 'pt'), legend.key.height = unit(10, 'pt')) +
							labs(fill = 'Correlation\nwith purity')
					),
					get_legend(
						simulated_deconv_plots[[1]]$plots$ccle_bar +
							theme(legend.justification = c(0, 1), legend.direction = 'vertical', legend.key.width = unit(10, 'pt'), legend.key.height = unit(10, 'pt')) +
							labs(fill = 'Tumours vs.\ncell lines')
					),
					get_legend(
						simulated_deconv_plots[[1]]$plots$heatmap + theme(legend.justification = c(0, 1))
					),
					get_legend(sc_sim_deconv_comp[[1]]$sc_heatmaps_unfiltered[[1]] + theme(legend.justification = c(0, 1)))
				)
			),
			nrow = 5,
			ncol = 1,
			rel_heights = c(1, 5, 5, 7, 10)
		),
		nrow = 1,
		ncol = 2,
		rel_widths = c(5, 1)
	) %>% print

}

dev.off()

cairo_pdf('../data_and_figures/sc_sim_deconv_comp_filtered.pdf', width = 6, height = 9, onefile = TRUE)

for(ct in names(sc_sim_deconv_comp)[!grepl('lenient', names(sc_sim_deconv_comp))]) {

	sc_sim_deconv_comp_figures <- sapply(
		c(
			sc_sim_deconv_comp[[ct]]$filtered_deconv_figures[c('purity_bar', 'ccle_bar', 'heatmap')],
			sc_sim_deconv_comp[[ct]]$sc_heatmaps_filtered
		),
		function(x) x + theme(legend.position = 'none', plot.title = element_blank()),
		simplify = FALSE,
		USE.NAMES = TRUE
	)
	sc_sim_deconv_comp_figures$heatmap <- sc_sim_deconv_comp_figures$heatmap + theme(plot.margin = unit(c(0, 5.5, 1.25, 5.5), 'pt'))

	plot_grid(
		plot_grid(
			plotlist = c(
				list(blank_plot() + theme(plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), 'pt')) + labs(title = deconv_plot_args_per_ct[[ct]]$plot_title)),
				sc_sim_deconv_comp_figures
				# lapply(
					# c(
						# sc_sim_deconv_comp[[ct]]$filtered_deconv_figures[c('purity_bar', 'ccle_bar', 'heatmap')],
						# sc_sim_deconv_comp[[ct]]$sc_heatmaps_filtered
					# ),
					# function(x) x + theme(legend.position = 'none', plot.title = element_blank())
				# )
			),
			nrow = 6,
			ncol = 1,
			align = 'v',
			rel_heights = c(1, 1, 1, 15, 5, 5)
		),
		plot_grid(
			plotlist = c(
				list(
					blank_plot(),
					get_legend(
						simulated_deconv_plots[[1]]$plots$purity_bar +
							theme(legend.justification = c(0, 1), legend.direction = 'vertical') +#, legend.key.width = NULL, legend.key.height = NULL) +
							labs(fill = 'Correlation\nwith purity')
					),
					get_legend(
						simulated_deconv_plots[[1]]$plots$ccle_bar +
							theme(legend.justification = c(0, 1), legend.direction = 'vertical') +#, legend.key.width = NULL, legend.key.height = NULL) +
							labs(fill = 'Tumours vs.\ncell lines')
					),
					get_legend(
						simulated_deconv_plots[[1]]$plots$heatmap + theme(legend.justification = c(0, 1))
					),
					get_legend(sc_sim_deconv_comp[[1]]$sc_heatmaps_filtered[[1]] + theme(legend.justification = c(0, 1)))
				)
			),
			nrow = 5,
			ncol = 1,
			rel_heights = c(1, 5, 5, 7, 10)
		),
		nrow = 1,
		ncol = 2,
		rel_widths = c(5, 1)
	) %>% print

}

dev.off()

cairo_pdf('../data_and_figures/sc_sim_deconv_comp_lenient_filtered.pdf', width = 6, height = 10, onefile = TRUE)

for(ct in names(sc_sim_deconv_comp)[grepl('lenient', names(sc_sim_deconv_comp))]) {

	sc_sim_deconv_comp_lenient_figures <- sapply(
		c(
			sc_sim_deconv_comp[[ct]]$filtered_deconv_figures[c('purity_bar', 'ccle_bar', 'heatmap')],
			sc_sim_deconv_comp[[ct]]$sc_heatmaps_filtered
		),
		function(x) x + theme(legend.position = 'none', plot.title = element_blank()),
		simplify = FALSE,
		USE.NAMES = TRUE
	)
	sc_sim_deconv_comp_lenient_figures$heatmap <- sc_sim_deconv_comp_lenient_figures$heatmap + theme(plot.margin = unit(c(0, 5.5, 1.25, 5.5), 'pt'))

	plot_grid(
		plot_grid(
			plotlist = c(
				list(blank_plot() + theme(plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), 'pt')) + labs(title = deconv_plot_args_per_ct[[ct]]$plot_title)),
				sc_sim_deconv_comp_lenient_figures
				# lapply(
					# c(
						# sc_sim_deconv_comp[[ct]]$filtered_deconv_figures[c('purity_bar', 'ccle_bar', 'heatmap')],
						# sc_sim_deconv_comp[[ct]]$sc_heatmaps_filtered
					# ),
					# function(x) x + theme(legend.position = 'none', plot.title = element_blank())
				# )
			),
			nrow = 6,
			ncol = 1,
			align = 'v',
			rel_heights = c(1, 1, 1, 15, 5, 5)
		),
		plot_grid(
			plotlist = c(
				list(
					blank_plot(),
					get_legend(
						simulated_deconv_plots[[1]]$plots$purity_bar +
							theme(legend.justification = c(0, 1), legend.direction = 'vertical') +#, legend.key.width = NULL, legend.key.height = NULL) +
							labs(fill = 'Correlation\nwith purity')
					),
					get_legend(
						simulated_deconv_plots[[1]]$plots$ccle_bar +
							theme(legend.justification = c(0, 1), legend.direction = 'vertical') +#, legend.key.width = NULL, legend.key.height = NULL) +
							labs(fill = 'Tumours vs.\ncell lines')
					),
					get_legend(
						simulated_deconv_plots[[1]]$plots$heatmap + theme(legend.justification = c(0, 1))
					),
					get_legend(sc_sim_deconv_comp[[1]]$sc_heatmaps_filtered[[1]] + theme(legend.justification = c(0, 1)))
				)
			),
			nrow = 5,
			ncol = 1,
			rel_heights = c(1, 5, 5, 7, 10)
		),
		nrow = 1,
		ncol = 2,
		rel_widths = c(5, 1)
	) %>% print

}

dev.off()





# Combine lineplots, deconv and single cell comparison to make Figure 2:

pemt_caf_brackets_params <- list(
	coadread = list(pemt_bracket_max = 30, caf_bracket_min = 71, pemt_label_hjust = 0.45, caf_label_hjust = 0.5),
	hnsc = list(pemt_bracket_max = 36, caf_bracket_min = 73, pemt_label_hjust = 0.5, caf_label_hjust = 0.55),
	luad = list(pemt_bracket_max = 29, caf_bracket_min = 71, pemt_label_hjust = 0.45, caf_label_hjust = 0.5)
)

cairo_pdf('../data_and_figures/final_figures_resubmission/2.pdf', width = 13, height = 11)

print(
	plot_grid(
		blank_plot(),
		plot_grid(
			plot_grid(
				plotlist = list(
					lineplots$coadread$lineplot +
						theme(legend.position = 'none', plot.margin = unit(c(5.5, 5.5, 15, 30), 'pt')),# +
						# labs(x = NULL, y = NULL),
					lineplots$hnsc$lineplot +
						theme(legend.position = 'none', plot.margin = unit(c(5.5, 5.5, 15, 15), 'pt')) +
						labs(y = NULL),
					lineplots$luad$lineplot +
						theme(legend.position = 'none', plot.margin = unit(c(5.5, 5.5, 15, 15), 'pt')) +
						labs(y = NULL)
				),
				nrow = 1,
				ncol = 3,
				rel_widths = c(1.125, 1, 1), # Left plots are squashed by larger left margin
				align = 'h'
			),
			dummy_legend,
			blank_plot(),
			nrow = 1,
			ncol = 3,
			rel_widths = c(8, 1.5, 0.5)
		),# +
			# draw_label('Proportion of tumour', x = 0.47, y = 0, vjust = -0.5, size = 12) +
			# draw_label('Proportion of gene expression', x = 0, y = 0.5, vjust = 1.3, angle = 90, size = 12),
		blank_plot(),
		plot_grid(
			plot_grid(
				plotlist = lapply(
					c('coadread', 'hnsc', 'luad'),
					function(ct) {

						sc_sim_deconv_comp_figures <- sapply(
							c(sc_sim_deconv_comp[[ct]]$filtered_deconv_figures, sc_sim_deconv_comp[[ct]]$sc_heatmaps_filtered),
							function(x) x + theme(legend.position = 'none', plot.title = element_blank(), axis.title.y = element_blank()),
							simplify = FALSE,
							USE.NAMES = TRUE
						)
						sc_sim_deconv_comp_figures$heatmap <- sc_sim_deconv_comp_figures$heatmap +
							theme(plot.margin = unit(c(0, 5.5, 1.25, 0), 'pt'))
						sc_sim_deconv_comp_figures$axis_labels <- sc_sim_deconv_comp_figures$axis_labels +
							theme(plot.margin = unit(c(0, 0, 1.25, 0), 'pt'))
						# sc_sim_deconv_comp_figures$sc_heatmaps_filtered <- sapply(
							# sc_sim_deconv_comp_figures$sc_heatmaps_filtered,
							# function(x) {x + theme(plot.margin = unit(c(5.5, 5.5, 1.25, 0), 'pt'))},
							# simplify = FALSE,
							# USE.NAMES = TRUE
						# )

						plot_grid(
							plot_grid(
								blank_plot(),
								pemt_caf_brackets(
									edge = 'right',
									pemt_bracket_max = pemt_caf_brackets_params[[ct]]$pemt_bracket_max,
									caf_bracket_min = pemt_caf_brackets_params[[ct]]$caf_bracket_min,
									brackets_just = 0.8,
									labels = c('pEMT genes', 'CAF genes'),
									pemt_label_hjust = pemt_caf_brackets_params[[ct]]$pemt_label_hjust,
									caf_label_hjust = pemt_caf_brackets_params[[ct]]$caf_label_hjust,
									labels_vjust = -0.7,
									labels_angle = 90,
									plot_margin = c(0, 0, 1.25, 0)
								),
								blank_plot(),
								nrow = 3,
								ncol = 1,
								rel_heights = c(4, 15, 10)
							),
							plot_grid(
								blank_plot(),
								sc_sim_deconv_comp_figures$axis_labels,
								blank_plot() +
									labs(y = 'Cancer\ncells') +
									scale_y_continuous(position = 'right') +
									theme(plot.margin = unit(c(5.5, 5.5, 1.25, 5.5), 'pt'), axis.title.y.right = element_text(angle = 90)),
								blank_plot() +
									labs(y = 'CAFs') +
									scale_y_continuous(position = 'right') +
									theme(plot.margin = unit(c(5.5, 5.5, 1.25, 5.5), 'pt'), axis.title.y.right = element_text(angle = 90)),
								nrow = 4,
								ncol = 1,
								rel_heights = c(4, 15, 5, 5)
							),
							plot_grid(
								plotlist = c(
									list(
										blank_plot() +
										theme(plot.margin = unit(c(11, 5.5, 5.5, 0), 'pt')) +
										labs(title = deconv_plot_args_per_ct[[ct]]$plot_title)
									),
									sc_sim_deconv_comp_figures[c('purity_bar', 'ccle_bar', 'heatmap', 'cancer', 'caf')]
								),
								nrow = 6,
								ncol = 1,
								align = 'v',
								rel_heights = c(2, 1, 1, 15, 5, 5)
							),
							nrow = 1,
							ncol = 3,
							rel_widths = c(0.8, 1.2, 4)
						)

					}
				),
				nrow = 1,
				ncol = 3
			),
			blank_plot(),
			plot_grid(
				get_legend(
					sc_sim_deconv_comp[[ct]]$filtered_deconv_figures$purity_bar +
						theme(legend.justification = c(1, 1), legend.direction = 'horizontal', legend.key.width = NULL) +
						labs(fill = 'Correlation with purity\n')
				),
				blank_plot(),
				get_legend(
					sc_sim_deconv_comp[[ct]]$filtered_deconv_figures$heatmap +
						guides(fill = guide_colourbar(title.position = 'right')) +
						# scale_fill_gradientn(
							# colours = rev(colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))(50)),
							# limits = c(-1, 1),
							# breaks = c(-1, 0, 1)
						# ) +
						theme(
							legend.justification = c(0, 1),
							legend.direction = 'horizontal',
							legend.key.width = unit(25, 'pt'),
							legend.key.height = unit(10, 'pt')
						) +
						labs(fill = 'Correlation\n')
				),
				get_legend(
					sc_sim_deconv_comp[[ct]]$filtered_deconv_figures$ccle_bar +
						theme(legend.justification = c(1, 1), legend.direction = 'horizontal', legend.key.width = NULL) +
						labs(fill = 'Tumours vs. cell lines\n')
				),
				blank_plot(),
				get_legend(
					sc_sim_deconv_comp[[ct]]$sc_heatmaps_filtered[[1]] +
						theme(
							legend.justification = c(0, 1),
							legend.direction = 'horizontal',
							legend.key.width = unit(25, 'pt'),
							legend.key.height = unit(10, 'pt')
						) +
						guides(fill = guide_colourbar(title.position = 'right')) +
						labs(fill = 'Relative expression level\n')
				),
				nrow = 2,
				ncol = 3,
				rel_widths = c(5, 1, 5)
			),
			nrow = 3,
			ncol = 1,
			rel_heights = c(1, 0.05, 0.2)
		),
		nrow = 4,
		ncol = 1,
		rel_heights =c(0.125, 0.95, 0.125, 1.8)
	) +
		draw_label('A', x = 0, y = 0.99, hjust = -0.1, vjust = 1.1, size = 20, fontface = 2) +
		draw_label('B', x = 0, y = 0.62, hjust = -0.1, vjust = 1.1, size = 20, fontface = 2)
)

dev.off()





# Combine remaining lineplots, deconv and single cell comparison to make Figure S8 (lineplots in A, deconv in B, on separate pages):

cairo_pdf('../data_and_figures/final_figures_resubmission/S8A.pdf', width = 10.4, height = 7.5)

print(
	plot_grid(
		blank_plot(),
		plot_grid(
			plotlist = list(
				lineplots$brca$lineplot +
					theme(legend.position = 'none', plot.margin = unit(c(5.5, 5.5, 15, 30), 'pt')),# +
					# labs(x = NULL, y = NULL),
				lineplots$lihc$lineplot +
					theme(legend.position = 'none', plot.margin = unit(c(5.5, 5.5, 15, 15), 'pt')) +
					labs(y = NULL),
				lineplots$lusc$lineplot +
					theme(legend.position = 'none', plot.margin = unit(c(5.5, 5.5, 15, 15), 'pt')) +
					labs(y = NULL)
			),
			nrow = 1,
			ncol = 3,
			rel_widths = c(1.125, 1, 1) # Left plots are squashed by larger left margin
			# align = 'h'
		),
		plot_grid(
			plotlist = list(
				lineplots$ov$lineplot +
					theme(legend.position = 'none', plot.margin = unit(c(5.5, 5.5, 15, 30), 'pt')),# +
					# labs(x = NULL, y = NULL),
				lineplots$paad$lineplot +
					theme(legend.position = 'none', plot.margin = unit(c(5.5, 5.5, 15, 15), 'pt')) +
					labs(y = NULL),
				dummy_legend
			),
			nrow = 1,
			ncol = 3,
			rel_widths = c(1.125, 1, 1) # Left plots are squashed by larger left margin
			# align = 'h'
		),
		nrow = 3,
		ncol = 1,
		rel_heights = c(0.125, 0.95, 0.95)
	) + draw_label('A', x = 0, y = 0.99, hjust = -0.1, vjust = 1.1, size = 20, fontface = 2)
)

dev.off()

pemt_caf_brackets_params <- list(
	brca = list(pemt_bracket_max = 36, caf_bracket_min = 71, pemt_label_hjust = 0.5, caf_label_hjust = 0.5),
	lihc = list(pemt_bracket_max = 25, caf_bracket_min = 71, pemt_label_hjust = 0.39, caf_label_hjust = 0.5),
	lusc = list(pemt_bracket_max = 36, caf_bracket_min = 73, pemt_label_hjust = 0.5, caf_label_hjust = 0.52),
	ov = list(pemt_bracket_max = 36, caf_bracket_min = 71, pemt_label_hjust = 0.5, caf_label_hjust = 0.5),
	paad = list(pemt_bracket_max = 29, caf_bracket_min = 71, pemt_label_hjust = 0.45, caf_label_hjust = 0.5)
)

cairo_pdf('../data_and_figures/final_figures_resubmission/S8B.pdf', width = 13, height = 11.5)

print(
	plot_grid(
		blank_plot(),
		plot_grid(
			plotlist = lapply(
				c('brca', 'lihc', 'lusc'),
				function(ct) {

					sc_sim_deconv_comp_figures <- sapply(
						c(sc_sim_deconv_comp[[ct]]$filtered_deconv_figures, sc_sim_deconv_comp[[ct]]$sc_heatmaps_filtered),
						function(x) x + theme(legend.position = 'none', plot.title = element_blank(), axis.title.y = element_blank()),
						simplify = FALSE,
						USE.NAMES = TRUE
					)
					sc_sim_deconv_comp_figures$heatmap <- sc_sim_deconv_comp_figures$heatmap +
						theme(plot.margin = unit(c(0, 5.5, 1.25, 0), 'pt'))
					sc_sim_deconv_comp_figures$axis_labels <- sc_sim_deconv_comp_figures$axis_labels +
						theme(plot.margin = unit(c(0, 0, 1.25, 0), 'pt'))

					plot_grid(
						plot_grid(
							blank_plot(),
							pemt_caf_brackets(
								edge = 'right',
								pemt_bracket_max = pemt_caf_brackets_params[[ct]]$pemt_bracket_max,
								caf_bracket_min = pemt_caf_brackets_params[[ct]]$caf_bracket_min,
								brackets_just = 0.8,
								labels = c('pEMT genes', 'CAF genes'),
								pemt_label_hjust = pemt_caf_brackets_params[[ct]]$pemt_label_hjust,
								caf_label_hjust = pemt_caf_brackets_params[[ct]]$caf_label_hjust,
								labels_vjust = -0.7,
								labels_angle = 90,
								plot_margin = c(0, 0, 1.25, 0)
							),
							blank_plot(),
							nrow = 3,
							ncol = 1,
							rel_heights = c(4, 15, 10)
						),
						plot_grid(
							blank_plot(),
							sc_sim_deconv_comp_figures$axis_labels,
							blank_plot() +
								labs(y = 'Cancer\ncells') +
								scale_y_continuous(position = 'right') +
								theme(plot.margin = unit(c(5.5, 5.5, 1.25, 5.5), 'pt'), axis.title.y.right = element_text(angle = 90)),
							blank_plot() +
								labs(y = 'CAFs') +
								scale_y_continuous(position = 'right') +
								theme(plot.margin = unit(c(5.5, 5.5, 1.25, 5.5), 'pt'), axis.title.y.right = element_text(angle = 90)),
							nrow = 4,
							ncol = 1,
							rel_heights = c(4, 15, 5, 5)
						),
						plot_grid(
							plotlist = c(
								list(
									blank_plot() +
									theme(plot.margin = unit(c(11, 5.5, 5.5, 0), 'pt')) +
									labs(title = deconv_plot_args_per_ct[[ct]]$plot_title)
								),
								sc_sim_deconv_comp_figures[c('purity_bar', 'ccle_bar', 'heatmap', 'cancer', 'caf')]
							),
							nrow = 6,
							ncol = 1,
							align = 'v',
							rel_heights = c(2, 1, 1, 15, 5, 5)
						),
						nrow = 1,
						ncol = 3,
						rel_widths = c(0.8, 1.2, 4)
					)

				}
			),
			nrow = 1,
			ncol = 3
		),
		blank_plot(),
		plot_grid(
			plotlist = c(
				lapply(
					c('ov', 'paad'),
					function(ct) {

						sc_sim_deconv_comp_figures <- sapply(
							c(sc_sim_deconv_comp[[ct]]$filtered_deconv_figures, sc_sim_deconv_comp[[ct]]$sc_heatmaps_filtered),
							function(x) x + theme(legend.position = 'none', plot.title = element_blank(), axis.title.y = element_blank()),
							simplify = FALSE,
							USE.NAMES = TRUE
						)
						sc_sim_deconv_comp_figures$heatmap <- sc_sim_deconv_comp_figures$heatmap +
							theme(plot.margin = unit(c(0, 5.5, 1.25, 0), 'pt'))
						sc_sim_deconv_comp_figures$axis_labels <- sc_sim_deconv_comp_figures$axis_labels +
							theme(plot.margin = unit(c(0, 0, 1.25, 0), 'pt'))

						plot_grid(
							plot_grid(
								blank_plot(),
								pemt_caf_brackets(
									edge = 'right',
									pemt_bracket_max = pemt_caf_brackets_params[[ct]]$pemt_bracket_max,
									caf_bracket_min = pemt_caf_brackets_params[[ct]]$caf_bracket_min,
									brackets_just = 0.8,
									labels = c('pEMT genes', 'CAF genes'),
									pemt_label_hjust = pemt_caf_brackets_params[[ct]]$pemt_label_hjust,
									caf_label_hjust = pemt_caf_brackets_params[[ct]]$caf_label_hjust,
									labels_vjust = -0.7,
									labels_angle = 90,
									plot_margin = c(0, 0, 1.25, 0)
								),
								blank_plot(),
								nrow = 3,
								ncol = 1,
								rel_heights = c(4, 15, 10)
							),
							plot_grid(
								blank_plot(),
								sc_sim_deconv_comp_figures$axis_labels,
								blank_plot() +
									labs(y = 'Cancer\ncells') +
									scale_y_continuous(position = 'right') +
									theme(plot.margin = unit(c(5.5, 5.5, 1.25, 5.5), 'pt'), axis.title.y.right = element_text(angle = 90)),
								blank_plot() +
									labs(y = 'CAFs') +
									scale_y_continuous(position = 'right') +
									theme(plot.margin = unit(c(5.5, 5.5, 1.25, 5.5), 'pt'), axis.title.y.right = element_text(angle = 90)),
								nrow = 4,
								ncol = 1,
								rel_heights = c(4, 15, 5, 5)
							),
							plot_grid(
								plotlist = c(
									list(
										blank_plot() +
										theme(plot.margin = unit(c(11, 5.5, 5.5, 0), 'pt')) +
										labs(title = deconv_plot_args_per_ct[[ct]]$plot_title)
									),
									sc_sim_deconv_comp_figures[c('purity_bar', 'ccle_bar', 'heatmap', 'cancer', 'caf')]
								),
								nrow = 6,
								ncol = 1,
								align = 'v',
								rel_heights = c(2, 1, 1, 15, 5, 5)
							),
							nrow = 1,
							ncol = 3,
							rel_widths = c(0.8, 1.2, 4)
						)

					}
				),
				list(
					plot_grid(
						blank_plot(),
						get_legend(
							sc_sim_deconv_comp$brca$filtered_deconv_figures$purity_bar +
								guides(fill = guide_colourbar(title.position = 'right')) +
								theme(
									legend.justification = c(0, 1),
									legend.direction = 'horizontal',
									legend.key.width = NULL,
									legend.box.margin = margin(l = 10)
								) +
								labs(fill = 'Correlation with purity\n')
						),
						get_legend(
							sc_sim_deconv_comp$brca$filtered_deconv_figures$ccle_bar +
								guides(fill = guide_colourbar(title.position = 'right')) +
								theme(
									legend.justification = c(0, 1),
									legend.direction = 'horizontal',
									legend.key.width = NULL,
									legend.box.margin = margin(l = 10)
								) +
								labs(fill = 'Tumours vs. cell lines\n')
						),
						get_legend(
							sc_sim_deconv_comp$brca$filtered_deconv_figures$heatmap +
								guides(fill = guide_colourbar(title.position = 'right')) +
								theme(
									legend.justification = c(0, 1),
									legend.direction = 'horizontal',
									legend.key.width = unit(25, 'pt'),
									legend.key.height = unit(10, 'pt'),
									legend.box.margin = margin(l = 10)
								) +
								labs(fill = 'Correlation\n')
						),
						get_legend(
							sc_sim_deconv_comp$brca$sc_heatmaps_filtered[[1]] +
								guides(fill = guide_colourbar(title.position = 'right')) +
								theme(
									legend.justification = c(0, 1),
									legend.direction = 'horizontal',
									legend.key.width = unit(25, 'pt'),
									legend.key.height = unit(10, 'pt'),
									legend.box.margin = margin(l = 10)
								) +
								labs(fill = 'Relative expression level\n')
						),
						blank_plot(),
						nrow = 6,
						ncol = 1,
						rel_heights = c(1, 1, 1, 1, 1, 3)
					)
				)
			),
			nrow = 1,
			ncol = 3
		),
		blank_plot(),
		nrow = 5,
		ncol = 1,
		rel_heights = c(0.125, 1.8, 0.125, 1.8, 0.125)
		# rel_heights = c(0.05, 1, 0.05, 1)
	) + draw_label('B', x = 0, y = 0.99, hjust = -0.1, vjust = 1.1, size = 20, fontface = 2)
)

dev.off()





# Figure to show single cell analysis using lenient CAF definitions, for 3 cancer types, in the same format as for analysis of
# scran-normalised data:

sc_cancer_caf_heatmap_lenient_final <- sapply(
    names(sc_cancer_caf)[grepl('_lenient', names(sc_cancer_caf))],
    function(ct) {
		cat(paste0(ct, '\n'))
        set.seed(sc_cancer_caf_args[[ct]]$seed)
        do.call(
            sc_groups_heatmap,
            args = list(
				sc_groups_list = sc_cancer_caf[[ct]],
				groups = c('cancer', 'caf'),
				x_axis_titles = c('Cancer cells', 'CAFs'),
				default_figure_widths = list(annotations = 2.2, cancer = 6, caf = 1.5),
				figure_spacing = 2.5,
				annotations_title_size = 18,
				annotations_nudge = 0.25,
				es_fun = NULL,
				es_title = 'EMT\nscore',
				h_legend_title = 'Expression level\n',
				h_legend_width = 20,
				h_legend_height = 10,
				h_legend_direction = 'horizontal',
				h_legend_title_position = 'right',
				h_legend_just = 'left',
				gd_legend_title = 'Genes detected\n',
				gd_legend_width = 20,
				gd_legend_height = 10,
				gd_legend_direction = 'horizontal',
				gd_legend_title_position = 'left',
				gd_legend_just = 'right',
				annotations = sc_cancer_caf_heatmaps_args[[ct]]$annotations,
				title = sc_cancer_caf_heatmaps_args[[ct]]$annotations_title
			)
        )
    },
    simplify = FALSE,
    USE.NAMES = TRUE
)

cairo_pdf(
    '../data_and_figures/final_figures_resubmission/S6.pdf',
    width = 12,
    height = 15.5
)

print(
	plot_grid(
		blank_plot(),
		plot_grid(
			plot_grid(
				cowplot_sc(
					sc_cancer_caf_heatmap_lenient_final$brca_lenient,
					legend_space = 0,
					heights = c(4, 20, 6),
					es_x_axis_title_vjust = 1.3,
					es_y_axis_title_angle = 0,
					es_y_axis_title_xpos = 0.8
				),
				# blank_plot(),
				cowplot_sc(
					sc_cancer_caf_heatmap_lenient_final$coadread_lenient,
					legend_space = 0,
					heights = c(4, 20, 6),
					es_x_axis_title_vjust = 1.3,
					es_y_axis_title_angle = 0,
					es_y_axis_title_xpos = 0.8
				),
				cowplot_sc(
					sc_cancer_caf_heatmap_lenient_final$lihc_lenient,
					legend_space = 0,
					heights = c(4, 20, 6),
					es_x_axis_title_vjust = 1.3,
					es_y_axis_title_angle = 0,
					es_y_axis_title_xpos = 0.8
				),
				ncol = 3,
				nrow = 1
				# rel_widths = c(20, 1, 20)
			),
			# blank_plot(),
			# plot_grid(
				# cowplot_sc(
					# all_figures$paad$sc_cancer_caf_heatmap_combining,
					# legend_space = 0,
					# heights = c(1.5, 20, 4),
					# es_x_axis_title_vjust = 1.3,
					# es_y_axis_title_angle = 0,
					# es_y_axis_title_xpos = 0.8
				# ),
				blank_plot(),
				plot_grid(
					sc_cancer_caf_heatmap_lenient_final$brca_lenient$plots$genes_detected_legend,
					# all_figures$paad$sc_cancer_caf_heatmap_scran$plots$genes_detected_legend,
					blank_plot(),
					sc_cancer_caf_heatmap_lenient_final$brca_lenient$plots$heatmap_legend,
					nrow = 1,
					ncol = 3,
					rel_widths = c(5, 1, 5)
				),
				# ncol = 3,
				# nrow = 1,
				# rel_widths = c(20, 1, 20)
			# ),
			ncol = 1,
			nrow = 3,
			rel_heights = c(20, 2, 3)
		),
		blank_plot(),
		plot_grid(
			plot_grid(
				plotlist = list(
					lineplots$brca_lenient$lineplot +
						theme(legend.position = 'none', plot.margin = unit(c(5.5, 5.5, 15, 30), 'pt')),# +
						# labs(x = NULL, y = NULL),
					lineplots$coadread_lenient$lineplot +
						theme(legend.position = 'none', plot.margin = unit(c(5.5, 5.5, 15, 15), 'pt')) +
						labs(y = NULL),
					lineplots$lihc_lenient$lineplot +
						theme(legend.position = 'none', plot.margin = unit(c(5.5, 5.5, 15, 15), 'pt')) +
						labs(y = NULL)
				),
				nrow = 1,
				ncol = 3,
				rel_widths = c(1.125, 1, 1, 1), # Left plots are squashed by larger left margin
				align = 'h'
			),
			dummy_legend,
			nrow = 1,
			ncol = 2,
			rel_widths = c(7.5, 1.5)
		),# +
			# draw_label('Proportion of tumour', x = 0.47, y = 0, vjust = -0.5, size = 12) +
			# draw_label('Proportion of gene expression', x = 0, y = 0.5, vjust = 1.3, angle = 90, size = 12),
		blank_plot(),
		plot_grid(
			plot_grid(
				plotlist = lapply(
					c('brca_lenient', 'coadread_lenient', 'lihc_lenient'),
					function(ct) {

						sc_sim_deconv_comp_figures <- sapply(
							c(sc_sim_deconv_comp[[ct]]$filtered_deconv_figures, sc_sim_deconv_comp[[ct]]$sc_heatmaps_filtered),
							function(x) x + theme(legend.position = 'none', plot.title = element_blank(), axis.title.y = element_blank()),
							simplify = FALSE,
							USE.NAMES = TRUE
						)
						sc_sim_deconv_comp_figures$heatmap <- sc_sim_deconv_comp_figures$heatmap +
							theme(plot.margin = unit(c(0, 5.5, 1.25, 0), 'pt'))
						sc_sim_deconv_comp_figures$axis_labels <- sc_sim_deconv_comp_figures$axis_labels +
							theme(plot.margin = unit(c(0, 0, 1.25, 0), 'pt'))

						plot_grid(
							plot_grid(
								blank_plot(),
								sc_sim_deconv_comp_figures$axis_labels,
								blank_plot() +
									labs(y = 'Cancer\ncells') +
									scale_y_continuous(position = 'right') +
									theme(plot.margin = unit(c(5.5, 5.5, 1.25, 5.5), 'pt'), axis.title.y.right = element_text(angle = 90)),
								blank_plot() +
									labs(y = 'CAFs') +
									scale_y_continuous(position = 'right') +
									theme(plot.margin = unit(c(5.5, 5.5, 1.25, 5.5), 'pt'), axis.title.y.right = element_text(angle = 90)),
								nrow = 4,
								ncol = 1,
								rel_heights = c(4, 15, 5, 5)
							),
							plot_grid(
								plotlist = c(
									list(
										blank_plot() +
										theme(plot.margin = unit(c(11, 5.5, 5.5, 0), 'pt')) +
										labs(title = sc_cancer_caf_heatmaps_args[[ct]]$annotations_title)
									),
									sc_sim_deconv_comp_figures[c('purity_bar', 'ccle_bar', 'heatmap', 'cancer', 'caf')]
								),
								nrow = 6,
								ncol = 1,
								align = 'v',
								rel_heights = c(2, 1, 1, 15, 5, 5)
							),
							nrow = 1,
							ncol = 2,
							rel_widths = c(1.6, 5)
						)

					}
				),
				nrow = 1,
				ncol = 3
			),
			blank_plot(),
			plot_grid(
				get_legend(
					sc_sim_deconv_comp$brca_lenient$filtered_deconv_figures$purity_bar +
						theme(legend.justification = c(1, 1), legend.direction = 'horizontal', legend.key.width = NULL) +
						labs(fill = 'Correlation with purity\n')
				),
				blank_plot(),
				get_legend(
					sc_sim_deconv_comp$brca_lenient$filtered_deconv_figures$heatmap +
						guides(fill = guide_colourbar(title.position = 'right')) +
						# scale_fill_gradientn(
							# colours = rev(colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))(50)),
							# limits = c(-1, 1),
							# breaks = c(-1, 0, 1)
						# ) +
						theme(
							legend.justification = c(0, 1),
							legend.direction = 'horizontal',
							legend.key.width = unit(25, 'pt'),
							legend.key.height = unit(10, 'pt')
						) +
						labs(fill = 'Correlation\n')
				),
				get_legend(
					sc_sim_deconv_comp$brca_lenient$filtered_deconv_figures$ccle_bar +
						theme(legend.justification = c(1, 1), legend.direction = 'horizontal', legend.key.width = NULL) +
						labs(fill = 'Tumours vs. cell lines\n')
				),
				blank_plot(),
				get_legend(
					sc_sim_deconv_comp$brca_lenient$sc_heatmaps_filtered[[1]] +
						theme(
							legend.justification = c(0, 1),
							legend.direction = 'horizontal',
							legend.key.width = unit(25, 'pt'),
							legend.key.height = unit(10, 'pt')
						) +
						guides(fill = guide_colourbar(title.position = 'right')) +
						labs(fill = 'Relative expression level\n')
				),
				nrow = 2,
				ncol = 3,
				rel_widths = c(5, 1, 5)
			),
			nrow = 3,
			ncol = 1,
			rel_heights = c(1, 0.05, 0.2)
		),
		nrow = 6,
		ncol = 1,
		rel_heights =c(0.15, 1.2, 0.125, 0.95, 0.125, 1.8)
	) +
		draw_label('A', x = 0, y = 0.99, hjust = -0.1, vjust = 1.1, size = 20, fontface = 2) +
		draw_label('B', x = 0, y = 0.68, hjust = -0.1, vjust = 1.1, size = 20, fontface = 2) +
		draw_label('C', x = 0, y = 0.43, hjust = -0.1, vjust = 1.1, size = 20, fontface = 2)
)

dev.off()





# Figures to show correlation of TCGA deconv results with patterns in single cell data:

scores_data <- deconv_scores(
    expression_data,
    deconv_data,
    scale_fun = function(x) x/(3*sd(x)),
    transform_data = TRUE
)

# These are depressingly bad:
sapply(names(deconv_data), function(ct) cor(scores_data[deconv_data[[ct]]$genes_filtered, ..ct], deconv_data[[ct]]$extra_data_score))
sapply(names(deconv_data), function(ct) cor(deconv_data[[ct]]$new_scores, deconv_data[[ct]]$extra_data_score))
# It's also a bit confusing, because HNSC mesenchymal-basal gets really low correlation, but the single cell pattern looks good for this one.





# Main figure with aligned lineplots and deconv plots:

pdf(
    '../data_and_figures/simulated_bulk_all_figures.pdf',
    width = 14,
    height = 8.5
)

plot_grid(
    blank_plot(),
    plot_grid(
        plotlist = c(
            lapply(
                c('hnsc', 'lung', 'paad', 'lihc'),
                function(ct) {
                    plot_grid(
                        plotlist = c(
                            list(
                                lineplots[[ct]]$lineplot +
                                    theme(
                                        legend.position = 'none',
                                        plot.margin = unit(c(5.5, 3.5, 55, 3.5), 'pt'),
                                        axis.text.y = switch((ct == 'hnsc') + 1, element_blank(), NULL),
                                        plot.title = element_text(size = 14.5)
                                    ) +
                                    labs(x = NULL, y = switch((ct == 'hnsc') + 1, NULL, waiver()))
                            ),
                            lapply(
                                simulated_deconv_plots[[ct]]$plots[
                                    c('purity_bar', 'ccle_bar', 'extra_bar')
                                    ],
                                function(g) {
                                    g + theme(
                                        legend.position = 'none',
                                        plot.margin = unit(c(0, 3.5, 2, 3.5), 'pt')
                                    )
                                }
                            ),
                            list(
                                simulated_deconv_plots[[ct]]$plots$heatmap +
                                    theme(
                                        legend.position = 'none',
                                        axis.title.y = element_text(vjust = -7.2),
                                        plot.margin = unit(c(0, 3.5, 0, 3.5), 'pt')
                                    ) +
                                    labs(title = NULL, y = switch((ct == 'hnsc') + 1, NULL, waiver())),
                                simulated_deconv_plots[[ct]]$plots$axis_labels
                            )
                        ),
                        nrow = 6,
                        ncol = 1,
                        rel_heights = c(23, 1, 1, 1, 15, 5),
                        align = 'v'
                    )
                }
            ),
            list(
                plot_grid(
                    plotlist = c(
                        list(get_legend(dummy_legend_plot + theme(legend.justification = 'left'))),
                        lapply(
                            c('purity_bar', 'ccle_bar', 'extra_bar', 'heatmap'),
                            function(plot_type) {
                                get_legend(
                                    simulated_deconv_plots$hnsc$plots[[plot_type]] +
                                        theme(
                                            legend.title = element_text(size = 9),
                                            legend.justification = c(0, switch((plot_type == 'heatmap') + 1, 1, 0))
                                        ) +
                                        guides(fill = guide_colourbar(title.position = 'right')) +
                                        labs(
                                            fill = mapvalues(
                                                plot_type,
                                                c('purity_bar', 'ccle_bar', 'extra_bar', 'heatmap'),
                                                c(
                                                    'Correlation\nwith purity',
                                                    'Tumours vs.\ncell lines',
                                                    'scRNA-seq:\nCAF vs. cancer',
                                                    'Correlation'
                                                ),
                                                warn_missing = FALSE
                                            )
                                        )
                                )
                            }
                        ),
                        list(blank_plot())
                    ),
                    nrow = 6,
                    ncol = 1,
                    rel_heights = c(23, 3, 3, 3, 9, 5)
                )
            )
        ),
        nrow = 1,
        ncol = 5,
        rel_widths = c(1.14, 1, 1, 1, 0.6)
    ),
    nrow = 2,
    ncol = 1,
    rel_heights = c(1, 20)
) +
    draw_label('A', x = 0, y = 0.975, hjust = -0.1, vjust = 1.1, size = 16, fontface = 2) +
    draw_label('B', x = 0, y = 0.515, hjust = -0.1, vjust = 1.1, size = 16, fontface = 2) +
    draw_label(
        'Proportion of tumour',
        x = 0.452,
        y = 0.55,
        size = 11
    )

dev.off()





# Make plots to compare EMT and CAF genes in single cell analysis vs simulated bulk
# deconvolution and vs TCGA bulk deconvolution:

simulated_deconvs <- readRDS('../data_and_figures/simulated_deconvs.rds')

simulated_bulk_data <- fread(
    '../data_and_figures/simulated_bulk_data.csv',
    key = 'id'
)

deconv_data <- readRDS('../data_and_figures/deconv_data.rds')

ct_to_keep <- c(
    # 'acc',
    'blca_luminal_papillary',
    'blca_basal_squamous',
    'brca_luminal_a',
    'brca_luminal_b',
    'brca_basal_like_karaayvaz',
    'brca_her2_enriched',
    # 'brca_normal_like',
    'cesc',
    # 'chol',
    'coad',
    'esca_ac',
    'esca_escc',
    'hnsc_mesenchymal_basal',
    'hnsc_classical',
    'hnsc_atypical',
    'lihc',
    'luad_proximal_inflammatory',
    'luad_proximal_proliferative',
    'lusc_basal',
    'lusc_classical',
    'lusc_primitive',
    'lusc_secretory',
    # 'meso',
    'ov_differentiated',
    'ov_immunoreactive',
    'ov_proliferative',
    'paad',
    'read',
    'skcm_immune',
    'skcm_keratin',
    'skcm_mitf_low',
    'stad_cin',
    'stad_ebv',
    'stad_msi',
    # 'thca',
    'ucec'
    # 'uvm'
)

deconv_data <- deconv_data[ct_to_keep]

nice_names_for_figure <- c(
    # 'ACC',
    'BLCA - Luminal-Papillary',
    'BLCA - Basal-Squamous',
    'BRCA - Luminal A',
    'BRCA - Luminal B',
    'BRCA - Basal-like',
    'BRCA - HER2-enriched',
    # 'BRCA - Normal-like',
    'CESC',
    # 'CHOL',
    'COADREAD - Colon',
    'ESCA - Adenocarcinoma',
    'ESCA - Squamous',
    'HNSC - Malignant-Basal',
    'HNSC - Classical',
    'HNSC - Atypical',
    'LIHC',
    'LUAD - Squamoid',
    'LUAD - Magnoid',
    'LUSC - Basal',
    'LUSC - Classical',
    'LUSC - Primitive',
    'LUSC - Secretory',
    # 'MESO',
    'OV - Differentiated',
    'OV - Immunoreactive',
    'OV - Proliferative',
    'PAAD',
    'COADREAD - Rectum',
    'SKCM - Immune',
    'SKCM - Keratin',
    'SKCM - MITF-low',
    'STAD - CIN',
    'STAD - EBV',
    'STAD - MSI',
    # 'THCA',
    'UCEC'
    # 'UVM'
)

names(deconv_data) <- mapvalues(
    names(deconv_data),
    names(deconv_data),
    nice_names_for_figure
)

# deconv_names <- names(deconv_data)[
#     !grepl('^UVM|^MESO|^ACC|^SKCM|^CHOL|^THCA', names(deconv_data))
# ]

deconv_names <- names(deconv_data)

sc_to_bulk_names <- list(
    hnsc = c(
        'HNSC - Malignant-Basal',
        'HNSC - Classical',
        'HNSC - Atypical'
    ),
    lung = c(
        'LUAD - Squamoid',
        'LUAD - Magnoid',
        'LUSC - Basal',
        'LUSC - Classical',
        'LUSC - Primitive',
        'LUSC - Secretory'
    ),
    brca = c(
        'BRCA - Luminal A',
        'BRCA - Luminal B',
        'BRCA - Basal-like',
        'BRCA - HER2-enriched'
        # 'BRCA - Normal-like'
    ),
    coadread = c(
        'COADREAD - Colon',
        'COADREAD - Rectum'
    ),
    paad = 'PAAD',
    lihc = 'LIHC',
    tnbc = c(
        'BRCA - Luminal A',
        'BRCA - Luminal B',
        'BRCA - Basal-like',
        'BRCA - HER2-enriched'
        # 'BRCA - Normal-like'
    )
)

# Get genes from single cell data and TCGA deconv data:

genes_sc <- unique(
    unlist(
        lapply(
            sc_cancer_caf,
            `[[`,
            'genes_filtered'
        )
    )
)

genes_tcga <- unique(
    unlist(
        lapply(
            deconv_data[unlist(sc_to_bulk_names)],
            `[[`,
            'genes_filtered'
        )
    )
)

genes_simulated <- unique(
    unlist(
        lapply(
            simulated_deconvs,
            `[[`,
            'genes_filtered'
        )
    )
)

# Calculate scores from transformed bulk (TCGA and simulated) expression data:

scores_data_tcga <- deconv_scores(
    expression_data,
    deconv_data[unlist(sc_to_bulk_names)],
    scale_fun = function(x) x/(3*sd(x[!is.na(x)])),
    transform_data = TRUE,
    additional_genes = unique(
        c(
            genes_sc[genes_sc %in% names(expression_data)],
            genes_simulated[genes_simulated %in% names(expression_data)]
        )
    )
)

# Replace NAs in simulated_bulk_data with zeros:

# simulated_bulk_data[
#     ,
#     names(simulated_bulk_data[, -'id']) := lapply(
#         .SD,
#         function(x) {mapvalues(x, NA, 0, warn_missing = FALSE)}
#     ),
#     .SDcols = -'id'
# ]

# There are some NAs in the following, which we allowed thanks to the altered scale_fun
# which ignores the NAs during the scaling.  Perhaps we should remove these genes
# altogether.

# Note that since the simulated bulk data was defined from the single cell data, all the
# genes in genes_sc will be in names(simulated_bulk_data).

scores_data_simulated <- deconv_scores(
    simulated_bulk_data,
    simulated_deconvs,
    scale_fun = function(x) x/(3*sd(x[!is.na(x)])),
    transform_data = TRUE,
    additional_genes = unique(
        c(
            genes_sc,
            genes_tcga[genes_tcga %in% names(simulated_bulk_data)]
        )
    )
)

scores <- sapply(
    names(sc_metadata),
    function(ct) {

        sc_data <- eval(sc_metadata[[ct]]$read_quote)

        all_genes <- unique(
            c(
                sc_cancer_caf[[ct]]$genes_filtered,
                simulated_deconvs[[ct]]$genes_filtered,
                unlist(
                    lapply(
                        deconv_data[sc_to_bulk_names[[ct]]],
                        `[[`,
                        'genes_filtered'
                    )
                )
            )
        )

        scores_tcga <- scores_data_tcga[
            all_genes,
            setNames(rowMeans(.SD), gene),
            .SDcols = sc_to_bulk_names[[ct]]
        ]

        scores_simulated <- scores_data_simulated[
            all_genes,
            setNames(get(ct), gene)
        ]

        # In the following, would it be sensible to reverse the log before taking means?

        # I was going to take the top and bottom 20 from the sorted vector of the following
        # gene scores, then score genes by correlation with these head and tail genes.  But
        # it's probably enough to use these scores as they are (or after log, to make them
        # normally distributed).

        scores_sc <- sc_data[
            cell_type == 'cancer',
            colMeans(.SD),
            .SDcols = all_genes[all_genes %in% names(sc_data)]
        ]/sc_data[
            cell_type == 'fibroblast',
            colMeans(.SD),
            .SDcols = all_genes[all_genes %in% names(sc_data)]
        ]

        out <- data.table(
            gene = all_genes,
            sc_score = scores_sc[all_genes],
            tcga_score = scores_tcga,
            simulated_score = scores_simulated
        )

        out[
            ,
            sc_score := log10(sc_score)
        ][
            ,
            c('sc_score', 'tcga_score', 'simulated_score') := lapply(
                .SD,
                function(x) {
                    sapply(
                        x,
                        function(y) {
                            switch(
                                (!is.na(y) & !is.infinite(y)) + 1,
                                NA,
                                (y - mean(x[!is.na(x) & !is.infinite(x)]))/
                                    sd(x[!is.na(x) & !is.infinite(x)])
                            )
                        }
                    )
                }
            ),
            .SDcols = -'gene'
        ]

        out

    },
    simplify = FALSE,
    USE.NAMES = TRUE
)





# Function to make combined plot for given pair of variable names:

plot_var_pair <- function(

    var_pair,
    scores_list,
    n_row,
    n_col,
    sc_genes = NULL,
    simulated_genes = NULL,
    tcga_genes = NULL,
    collate_genes_fun = intersect, # Could also use union
    annotate_coords = c(3/4, 20/21),
    plot_titles = names(scores_list),
    cor_method = 'spearman',
    var_axis_titles = var_pair,
    rel_widths = rep(1, n_col),
    rel_heights = rep(1, n_row)

) {

    collate_genes_fun <- match.fun(collate_genes_fun)

    n <- length(scores_list)

    if(
        !is.null(get(paste0(strsplit(var_pair[1], '_')[[1]][1], '_genes'))) &
        !is.null(get(paste0(strsplit(var_pair[2], '_')[[1]][1], '_genes')))
    ) {
        plot_data <- sapply(
            names(scores_list),
            function(ct) {
                scores_list[[ct]][
                    gene %in% collate_genes_fun(
                        get(paste0(strsplit(var_pair[1], '_')[[1]][1], '_genes'))[[ct]],
                        get(paste0(strsplit(var_pair[2], '_')[[1]][1], '_genes'))[[ct]]
                    ),
                    ..var_pair
                ][!is.na(get(var_pair[1])) & !is.na(get(var_pair[2]))]
            },
            simplify = FALSE,
            USE.NAMES = TRUE
        )
    } else {
        plot_data <- sapply(
            scores_list,
            function(dt) dt[, ..var_pair][!is.na(get(var_pair[1])) & !is.na(get(var_pair[2]))],
            simplify = FALSE,
            USE.NAMES = TRUE
        )
    }

    lims_data <- sapply(
        sapply(
            var_pair,
            function(v) sapply(plot_data, function(dt) c(floor(10*min(dt[[v]]))/10, ceiling(10*max(dt[[v]]))/10)),
            simplify = FALSE,
            USE.NAMES = TRUE
        ),
        function(m) c(min(m[1, ]), max(m[2, ]))
    )

    if(n == n_row*n_col) {

        plot_grid(
            plotlist = lapply(
                1:n,
                function(i) {
                    ggplot(plot_data[[i]], aes(get(var_pair[1]), get(var_pair[2]))) +
                        geom_vline(xintercept = 0, linetype = 'dashed', colour = 'grey', size = 0.25) +
                        geom_hline(yintercept = 0, linetype = 'dashed', colour = 'grey', size = 0.25) +
                        geom_point() +
                        geom_smooth(method = 'lm', se = FALSE) +
                        theme_test() +
                        lims(x = lims_data[, var_pair[1]], y = lims_data[, var_pair[2]]) +
                        annotate(
                            'text',
                            x = annotate_coords[1]*lims_data[, var_pair[1]][2],
                            y = annotate_coords[2]*lims_data[, var_pair[2]][1],
                            label = TeX(sprintf('$\\rho = %.2f$', plot_data[[i]][, cor(get(var_pair[1]), get(var_pair[2]), method = cor_method)]))
                        ) +
                        labs(title = plot_titles[i]) +
                        theme(
                            axis.text.x = switch((i %in% (n - n_col + 1):n) + 1, element_blank(), NULL),
                            axis.text.y = switch((i %in% ((0:((n / n_col) - 1))*n_col + 1)) + 1, element_blank(), NULL),
                            axis.title = element_blank(),
                            plot.margin = unit(
                                c(5.5, 5.5, switch((i %in% (n - n_col + 1):n) + 1, 5.5, 20), switch((i %in% ((0:((n / n_col) - 1))*n_col + 1)) + 1, 5.5, 20)),
                                'pt'
                            )
                        )
                }
            ),
            nrow = n_row,
            ncol = n_col,
            # Don't need align = 'h' any more because no plots that aren't on the bottom
            # row will have x axis text, and aligning introduces space where the x axis
            # text would go if it was there...
            rel_widths = rel_widths,
            rel_heights = rel_heights
        ) +
            draw_label(
                var_axis_titles[1],
                x = 0.5,
                y = 0,
                vjust = -0.5,
                size = 12
            ) +
            draw_label(
                var_axis_titles[2],
                x = 0,
                y = 0.5,
                vjust = 1.3,
                angle = 90,
                size = 12
            )

    } else {

        plot_grid(
            plot_grid(
                plotlist = lapply(
                    1:(n_col*(n %/% n_col)),
                    function(i) {
                        ggplot(plot_data[[i]], aes(get(var_pair[1]), get(var_pair[2]))) +
                            geom_vline(xintercept = 0, linetype = 'dashed', colour = 'grey', size = 0.25) +
                            geom_hline(yintercept = 0, linetype = 'dashed', colour = 'grey', size = 0.25) +
                            geom_point() +
                            geom_smooth(method = 'lm', se = FALSE) +
                            theme_test() +
                            lims(x = lims_data[, var_pair[1]], y = lims_data[, var_pair[2]]) +
                            annotate(
                                'text',
                                x = annotate_coords[1]*lims_data[, var_pair[1]][2],
                                y = annotate_coords[2]*lims_data[, var_pair[2]][1],
                                label = TeX(sprintf('$\\rho = %.2f$', plot_data[[i]][, cor(get(var_pair[1]), get(var_pair[2]), method = cor_method)]))
                            ) +
                            labs(title = plot_titles[i]) +
                            theme(
                                axis.text.x = switch((i %in% (n_col*(n %/% n_col) - (n_col - (n %% n_col)) + 1):(n_col*(n %/% n_col))) + 1, element_blank(), NULL),
                                axis.text.y = switch((i %in% ((0:((n %/% n_col) - 1))*n_col + 1)) + 1, element_blank(), NULL),
                                axis.title = element_blank(),
                                plot.margin = unit(c(5.5, 5.5, 5.5, switch((i %in% ((0:((n %/% n_col) - 1))*n_col + 1)) + 1, 5.5, 20)), 'pt')
                            )
                    }
                ),
                nrow = n %/% n_col,
                ncol = n_col,
                align = 'h',
                rel_widths = rel_widths,
                rel_heights = rel_heights[-length(rel_heights)]
            ),
            plot_grid(
                plotlist = lapply(
                    (n - (n %% n_col) + 1):n,
                    function(i) {
                        ggplot(plot_data[[i]], aes(get(var_pair[1]), get(var_pair[2]))) +
                            geom_vline(xintercept = 0, linetype = 'dashed', colour = 'grey', size = 0.25) +
                            geom_hline(yintercept = 0, linetype = 'dashed', colour = 'grey', size = 0.25) +
                            geom_point() +
                            geom_smooth(method = 'lm', se = FALSE) +
                            theme_test() +
                            lims(x = lims_data[, var_pair[1]], y = lims_data[, var_pair[2]]) +
                            annotate(
                                'text',
                                x = annotate_coords[1]*lims_data[, var_pair[1]][2],
                                y = annotate_coords[2]*lims_data[, var_pair[2]][1],
                                label = TeX(sprintf('$\\rho = %.2f$', plot_data[[i]][, cor(get(var_pair[1]), get(var_pair[2]), method = cor_method)]))
                            ) +
                            labs(title = plot_titles[i]) +
                            theme(
                                axis.text.y = switch((i == n - (n %% n_col) + 1) + 1, element_blank(), NULL),
                                axis.title = element_blank(),
                                plot.margin = unit(c(5.5, 5.5, 20, switch((i == n - (n %% n_col) + 1) + 1, 5.5, 20)), 'pt')
                            )
                    }
                ),
                nrow = 1,
                ncol = n_col,
                align = 'h',
                rel_widths = rel_widths
            ),
            nrow = 2,
            ncol = 1,
            rel_heights = c(
                sum(rel_heights[-length(rel_heights)]),
                rel_heights[length(rel_heights)]
            )
        ) +
            draw_label(
                var_axis_titles[1],
                x = 0.5,
                y = 0,
                vjust = -0.5,
                size = 12
            ) +
            draw_label(
                var_axis_titles[2],
                x = 0,
                y = 0.5,
                vjust = 1.3,
                angle = 90,
                size = 12
            )

    }

}





# Write plots to PDF:

pdf(
    '../data_and_figures/signature_concordance.pdf',
    width = 12,
    height = 6
)

plot_var_pair(
    c('sc_score', 'simulated_score'),
    scores,
    n_row = 2,
    n_col = 4,
    plot_titles = sapply(sc_cancer_caf_heatmaps_args, `[[`, 'annotations_title'),
    var_axis_titles = c('scRNA-seq score', 'Simulated deconvolution score'),
    rel_widths = c(1.115, 1, 1, 1),
    rel_heights = c(1, 1.075)
)

plot_var_pair(
    c('sc_score', 'tcga_score'),
    scores,
    n_row = 2,
    n_col = 4,
    plot_titles = sapply(sc_cancer_caf_heatmaps_args, `[[`, 'annotations_title'),
    var_axis_titles = c('scRNA-seq score', 'TCGA deconvolution score'),
    rel_widths = c(1.115, 1, 1, 1),
    rel_heights = c(1, 1.075)
)

plot_var_pair(
    c('simulated_score', 'tcga_score'),
    scores,
    n_row = 2,
    n_col = 4,
    plot_titles = sapply(sc_cancer_caf_heatmaps_args, `[[`, 'annotations_title'),
    var_axis_titles = c('Simulated deconvolution score', 'TCGA deconvolution score'),
    rel_widths = c(1.115, 1, 1, 1),
    rel_heights = c(1, 1.075)
)

# Alternative, using pairwise intersections of gene sets:

sc_genes <- sapply(
    sc_cancer_caf,
    `[[`,
    'genes_filtered',
    simplify = FALSE,
    USE.NAMES = TRUE
)

simulated_genes <- sapply(
    simulated_deconvs,
    `[[`,
    'genes_filtered',
    simplify = FALSE,
    USE.NAMES = TRUE
)

tcga_genes <- sapply(
    names(sc_to_bulk_names),
    function(ct) {
        unique(
            unlist(
                lapply(deconv_data[sc_to_bulk_names[[ct]]], `[[`, 'genes_filtered')
            )
        )
    },
    simplify = FALSE,
    USE.NAMES = TRUE
)

plot_var_pair(
    c('sc_score', 'simulated_score'),
    scores,
    n_row = 2,
    n_col = 4,
    sc_genes = sc_genes,
    simulated_genes = simulated_genes,
    plot_titles = sapply(sc_cancer_caf_heatmaps_args, `[[`, 'annotations_title'),
    var_axis_titles = c('scRNA-seq score', 'Simulated deconvolution score'),
    rel_widths = c(1.115, 1, 1, 1),
    rel_heights = c(1, 1.075)
)

plot_var_pair(
    c('sc_score', 'tcga_score'),
    scores,
    n_row = 2,
    n_col = 4,
    sc_genes = sc_genes,
    tcga_genes = tcga_genes,
    plot_titles = sapply(sc_cancer_caf_heatmaps_args, `[[`, 'annotations_title'),
    var_axis_titles = c('scRNA-seq score', 'TCGA deconvolution score'),
    rel_widths = c(1.115, 1, 1, 1),
    rel_heights = c(1, 1.075)
)

plot_var_pair(
    c('simulated_score', 'tcga_score'),
    scores,
    n_row = 2,
    n_col = 4,
    simulated_genes = simulated_genes,
    tcga_genes = tcga_genes,
    plot_titles = sapply(sc_cancer_caf_heatmaps_args, `[[`, 'annotations_title'),
    var_axis_titles = c('Simulated deconvolution score', 'TCGA deconvolution score'),
    rel_widths = c(1.115, 1, 1, 1),
    rel_heights = c(1, 1.075)
)

# Pairwise unions:

plot_var_pair(
    c('sc_score', 'simulated_score'),
    scores,
    n_row = 2,
    n_col = 4,
    sc_genes = sc_genes,
    simulated_genes = simulated_genes,
    collate_genes_fun = union,
    plot_titles = sapply(sc_cancer_caf_heatmaps_args, `[[`, 'annotations_title'),
    var_axis_titles = c('scRNA-seq score', 'Simulated deconvolution score'),
    rel_widths = c(1.115, 1, 1, 1),
    rel_heights = c(1, 1.075)
)

plot_var_pair(
    c('sc_score', 'tcga_score'),
    scores,
    n_row = 2,
    n_col = 4,
    sc_genes = sc_genes,
    tcga_genes = tcga_genes,
    collate_genes_fun = union,
    plot_titles = sapply(sc_cancer_caf_heatmaps_args, `[[`, 'annotations_title'),
    var_axis_titles = c('scRNA-seq score', 'TCGA deconvolution score'),
    rel_widths = c(1.115, 1, 1, 1),
    rel_heights = c(1, 1.075)
)

plot_var_pair(
    c('simulated_score', 'tcga_score'),
    scores,
    n_row = 2,
    n_col = 4,
    simulated_genes = simulated_genes,
    tcga_genes = tcga_genes,
    collate_genes_fun = union,
    plot_titles = sapply(sc_cancer_caf_heatmaps_args, `[[`, 'annotations_title'),
    var_axis_titles = c('Simulated deconvolution score', 'TCGA deconvolution score'),
    rel_widths = c(1.115, 1, 1, 1),
    rel_heights = c(1, 1.075)
)

dev.off()

# Single one (just TCGA deconv vs. scRNAseq) for paper supplementary:

pdf(
    '../data_and_figures/signature_concordance_sc_tcga.pdf',
    width = 6,
    height = 9
)

plot_var_pair(
    c('sc_score', 'tcga_score'),
    scores[c('hnsc', 'lung', 'paad', 'lihc', 'tnbc')],
    n_row = 3,
    n_col = 2,
    sc_genes = sc_genes,
    tcga_genes = tcga_genes,
    collate_genes_fun = union,
    plot_titles = sapply(
        sc_cancer_caf_heatmaps_args[c('hnsc', 'lung', 'paad', 'lihc', 'tnbc')],
        `[[`,
        'annotations_title'
    ),
    var_axis_titles = c('scRNA-seq score', 'TCGA deconvolution score'),
    rel_widths = c(1.115, 1),
    rel_heights = c(1, 1, 1.075)
)

dev.off()

pdf(
    '../data_and_figures/signature_concordance_sc_tcga_supp.pdf',
    width = 6,
    height = 3.2
)

plot_var_pair(
    c('sc_score', 'tcga_score'),
    scores[c('brca', 'coadread')],
    n_row = 1,
    n_col = 2,
    sc_genes = sc_genes,
    tcga_genes = tcga_genes,
    collate_genes_fun = union,
    plot_titles = sapply(
        sc_cancer_caf_heatmaps_args[c('brca', 'coadread')],
        `[[`,
        'annotations_title'
    ),
    var_axis_titles = c('scRNA-seq score', 'TCGA deconvolution score'),
    rel_widths = c(1.115, 1)
)

dev.off()

# gene_position_profiles <- sapply(
#     names(sc_to_bulk_names),
#     function(ct1) {
#         sapply(
#             sc_to_bulk_names[[ct1]],
#             function(ct2) {
#
#                 common_genes <- intersect(
#                     simulated_deconvs[[ct1]]$genes_filtered,
#                     deconv_data[[ct2]]$genes_filtered
#                 )
#
#                 ordered_genes_1 <- with(
#                     simulated_deconvs[[ct1]],
#                     genes_filtered[ordering][genes_filtered[ordering] %in% common_genes]
#                 )
#
#                 ordered_genes_2 <- with(
#                     deconv_data[[ct2]],
#                     genes_filtered[ordering][genes_filtered[ordering] %in% common_genes]
#                 )
#
#                 qplot(
#                     x = 1:length(common_genes),
#                     y = order(
#                         order(ordered_genes_2)[
#                             order(order(ordered_genes_1))
#                         ]
#                     ),
#                     geom = 'line',
#                     xlab = 'Position in TCGA deconv',
#                     ylab = 'Position in simulated deconv',
#                     main = paste(ct1, 'vs.', ct2)
#                 ) + theme_test()
#
#             },
#             simplify = FALSE,
#             USE.NAMES = TRUE
#         )
#     },
#     simplify = FALSE,
#     USE.NAMES = TRUE
# )
#
# plot_grid(
#     plotlist = unlist(gene_position_profiles, recursive = FALSE),
#     nrow = 7,
#     ncol = 3,
#     align = 'hv'
# )





# Another attempt at rare vs. shared EMT:

deconv_data <- readRDS('../data_and_figures/deconv_data.rds')

ct_to_keep <- c(
    'blca_luminal_papillary',
    'blca_basal_squamous',
    'brca_luminal_a',
    'brca_luminal_b',
    'brca_basal_like',
    'brca_her2_enriched',
    'cesc',
    'coad',
    'esca_ac',
    'esca_escc',
    'hnsc_mesenchymal_basal',
    'hnsc_classical',
    'hnsc_atypical',
    'lihc',
    'luad_proximal_inflammatory',
    'luad_proximal_proliferative',
    'lusc_basal',
    'lusc_classical',
    'lusc_primitive',
    'lusc_secretory',
    'ov_differentiated',
    'ov_immunoreactive',
    'ov_proliferative',
    'paad',
    'read',
    # 'skcm_immune',
    # 'skcm_keratin',
    # 'skcm_mitf_low',
    'stad_cin',
    'stad_ebv',
    'stad_msi',
    'ucec'
)

deconv_data <- deconv_data[ct_to_keep]

nice_names_for_figure <- c(
    'BLCA - Luminal-Papillary',
    'BLCA - Basal-Squamous',
    'BRCA - Luminal A',
    'BRCA - Luminal B',
    'BRCA - Basal-like',
    'BRCA - HER2-enriched',
    'CESC',
    'COAD',
    'ESCA - Adenocarcinoma',
    'ESCA - Squamous',
    'HNSC - Malignant-Basal',
    'HNSC - Classical',
    'HNSC - Atypical',
    'LIHC',
    'LUAD - Squamoid',
    'LUAD - Magnoid',
    'LUSC - Basal',
    'LUSC - Classical',
    'LUSC - Primitive',
    'LUSC - Secretory',
    'OV - Differentiated',
    'OV - Immunoreactive',
    'OV - Proliferative',
    'PAAD',
    'READ',
    # 'SKCM - Immune',
    # 'SKCM - Keratin',
    # 'SKCM - MITF-low',
    'STAD - CIN',
    'STAD - EBV',
    'STAD - MSI',
    'UCEC'
)

names(deconv_data) <- mapvalues(
    names(deconv_data),
    names(deconv_data),
    nice_names_for_figure
)





# The following takes a long time, so read from RDS if already done:

# inter_intra_emt_all_cts <- readRDS('../data_and_figures/inter_intra_emt_all_cts.rds')

inter_intra_emt_all_cts <- sapply(
    # names(sc_metadata),
    c('brca', 'coadread', 'hnsc', 'lihc', 'luad', 'paad'),
    function(ct) {

        sc_data <- eval(sc_metadata[[ct]]$read_quote)
        setkey(sc_data, id)

        # Take subset of tumours with at least 50 cancer cells:
        sc_data <- sc_data[patient %in% sc_data[cell_type == 'cancer', .(N = .N), by = patient][N >= 50, patient]]

        # Get all the sufficiently highly-expressed genes:
        all_genes_filtered <- sc_data[
            cell_type == 'cancer',
            names(.SD)[apply(.SD, 2, sc_cancer_caf_args[[ct]]$genes_filter_fun)],
            .SDcols = -c('id', 'patient', 'cell_type')
        ]

        # Compute bins for the genes based on average expression:
        sc_mat <- sc_data[cell_type == 'cancer', set_colnames(t(.SD), id), .SDcols = all_genes_filtered]
        gene_averages <- sort(rowMeans(sc_mat))
        bins <- setNames(
            cut(seq_along(gene_averages), breaks = length(gene_averages) %/% 110, labels = FALSE, include.lowest = TRUE),
            names(gene_averages)
        )

        # Get EMT markers:
        ct_emt_markers <- unique(
            unlist(
                lapply(
                    deconv_data[grepl(paste(paste0('^', sc_metadata[[ct]]$tcga_cancer_types), collapse = '|'), names(deconv_data))],
                    function(deconv) deconv$genes_filtered[deconv$ordering]
                    # function(deconv) head(deconv$genes_filtered[deconv$ordering], length(deconv$genes_filtered)/3)
                )
            )
        )
        ct_emt_markers <- ct_emt_markers[ct_emt_markers %in% names(sc_data)]

        # Filter EMT markers by expression levels:
        ct_emt_markers <- sc_data[
            cell_type == 'cancer',
            ct_emt_markers[apply(.SD, 2, sc_cancer_caf_args[[ct]]$scores_filter_fun)],
            .SDcols = ct_emt_markers
        ]

        # Define control gene sets for distribution of scores:
        comparable_gene_sets <- lapply(ct_emt_markers, function(g) sample(names(bins)[bins == bins[g]], 100))
        comparable_gene_sets <- lapply(1:100, function(i) sapply(comparable_gene_sets, `[`, i))

        # Calculate EMT scores, then filter the EMT markers for correlation with these EMT scores, and recalculate the EMT scores using the
		# filtered list.  I think it makes more sense to filter for correlation with the initial EMT score than with the Z score (below), since
		# the Z score won't be on a comparable scale to the original distribution of expression levels (though this probably doesn't matter
		# much), and the aim of the Z score is to remove variability arising from poor data quality, which might significantly affect the
		# correlation with the Z scores.  EDIT: I don't think there's any difference between the correlations with EMT scores and with Z scores.

        emt_scores <- rowMeans(
            sapply(
                ct_emt_markers,
                function(g) {sc_mat[g, ] - colMeans(sc_mat[sample(names(bins)[bins == bins[g]], 100), ])},
                USE.NAMES = TRUE
            )
        )

        # sc_data[cell_type == 'cancer', sapply(ct_emt_markers, function(g) cor(get(g), emt_scores))] %>% sort %>% plot

        ct_emt_markers <- sc_data[cell_type == 'cancer', ct_emt_markers[sapply(ct_emt_markers, function(g) cor(get(g), emt_scores)) > 0.3]]

        emt_scores <- rowMeans(
            sapply(
                ct_emt_markers,
                function(g) {sc_mat[g, ] - colMeans(sc_mat[sample(names(bins)[bins == bins[g]], 100), ])},
                USE.NAMES = TRUE
            )
        )

        # Calculate distribution of scores for control gene sets:

        comparable_gene_sets_scores <- lapply(
            comparable_gene_sets,
            function(gli) {
                rowMeans(
                    sapply(
                        gli,
                        function(g) {sc_mat[g, ] - colMeans(sc_mat[sample(names(bins)[bins == bins[g]], 100), ])},
                        USE.NAMES = TRUE
                    )
                )
            }
        )

        # Define the EMT score as a Z score:

        z_scores <- sapply(
            names(emt_scores),
            function(cell_id) {
                distrib <- c(emt_scores[cell_id], sapply(comparable_gene_sets_scores, `[`, cell_id))
                (emt_scores[cell_id] - mean(distrib))/sd(distrib)
            },
            USE.NAMES = FALSE
        )

        emt_scores <- sc_data[cell_type == 'cancer', .(id = id, patient = patient, emt_score = z_scores[id])][
            order(patient, emt_score),
            pos_frac := (1:.N)/.N,
            by = patient
        ]

        emt_score_line <- ggplot(
            emt_scores[pos_frac > 0.01 & pos_frac < 0.99],
            aes(pos_frac, emt_score, group = patient, colour = as.character(patient))
        ) +
            scale_colour_manual(values = brewer.pal(12, 'Set3')[c(1, 3:10, 12, 11)][1:length(unique(emt_scores$patient))]) +
            geom_line() +
            theme_test()

        inter_intra_emt <- emt_scores[
            ,
            .(
                # inter_emt = .SD[pos_frac > 0.25 & pos_frac < 0.75, mean(emt_score)],
                # intra_emt = .SD[pos_frac > 0.9 & pos_frac < 0.975, mean(emt_score)] - median(emt_score)
                inter_emt = median(emt_score),
                intra_emt = quantile(emt_score, 0.95) - median(emt_score)
            ),
            by = patient
        ]

        inter_intra_plot <- ggplot(inter_intra_emt, aes(inter_emt, intra_emt, colour = as.character(patient))) +
            scale_colour_manual(values = brewer.pal(12, 'Set3')[c(1, 3:10, 12, 11)][1:length(unique(emt_scores$patient))]) +
            geom_point() +
            theme_test()

        list(
            plots = list(lineplot = emt_score_line, scatterplot = inter_intra_plot),# violin = emt_score_violin),
            data = list(
                cell_emt_scores = emt_scores,
                inter_intra_emt_scores = inter_intra_emt,
                all_genes_filtered = all_genes_filtered,
                emt_markers_filtered = ct_emt_markers
            )
        )

    },
    simplify = FALSE,
    USE.NAMES = TRUE
)

saveRDS(inter_intra_emt_all_cts, '../data_and_figures/inter_intra_emt_all_cts.rds')

# emt_scores <- scrabble::score(
#     sc_data[
#         cell_type == 'cancer',
#         set_colnames(t(.SD), id),
#         .SDcols = all_genes_filtered
#     ],
#     list(ct_emt_markers),
#     bin.control = TRUE,
#     nbin = length(all_genes_filtered) %/% 110
# )

# emt_scores <- data.table(
#     id = rownames(emt_scores),
#     patient = sc_data[rownames(emt_scores), patient],
#     emt_score = emt_scores[, 1]
# )[
#     order(patient, emt_score),
#     pos_frac := (1:.N)/.N,
#     by = patient
# ]

# emt_score_violin <- ggplot(
	# emt_scores,
	# aes(
		# factor(patient, levels = emt_scores[, .(med = median(emt_score)), by = patient][order(med), patient]),
		# emt_score,
		# fill = as.character(patient)
	# )
# ) +
	# scale_fill_manual(values = brewer.pal(12, 'Set3')[1:length(unique(emt_scores$patient))]) +
	# geom_violin(draw_quantiles = 0.5) +
	# theme_test()





# Plots combining cancer types:

emt_scores_all_cts <- rbindlist(
    lapply(
        names(inter_intra_emt_all_cts),
        function(ct) {
            merge(
                cbind(
                    cancer_type = ct,
                    cancer_type_nice = sc_cancer_caf_heatmaps_args[[ct]]$annotations_title,
                    inter_intra_emt_all_cts[[ct]]$data$cell_emt_scores
                ),
                inter_intra_emt_all_cts[[ct]]$data$inter_intra_emt_scores,
                by = 'patient'
            )
        }
    )
)[
    ,
    patient_unique := paste(cancer_type, patient, sep = '.')
]

# inter_intra_emt_scores_all_cts <- rbindlist(
#     lapply(
#         names(inter_intra_emt_all_cts),
#         function(ct) {
#             cbind(
#                 cancer_type = ct,
#                 inter_intra_emt_all_cts[[ct]]$data$inter_intra_emt_scores
#             )
#         }
#     )
# )

# inter_intra_emt_profiles <- ggplot(
    # emt_scores_all_cts,#[pos_frac > 0.01 & pos_frac < 0.99],
    # aes(
        # pos_frac,
        # emt_score,
        # group = patient,
        # # group = interaction(cancer_type, patient),
        # colour = cancer_type_nice
        # # alpha = intra_emt
        # # size = intra_emt
    # )
# ) +
    # # scale_size_continuous(range = c(0.4, 0.8)) +
    # facet_grid(cols = vars(cancer_type_nice)) +
    # geom_line() +
    # scale_x_continuous(expand = c(0, 0)) +
    # scale_y_continuous(expand = c(0, 0), breaks = c(-5, 0, 5)) +
    # theme(
        # panel.background = element_rect(fill = NA, colour = 'black'),
        # panel.grid.major.y = element_line(colour = 'grey', size = 0.3, linetype = 'dotted'),
        # panel.grid.minor.y = element_blank(),
        # panel.grid.major.x = element_blank(),
        # panel.grid.minor.x = element_blank(),
        # axis.text.x = element_blank(),
        # axis.ticks.x = element_blank(),
        # strip.background = element_rect(colour = 'black'),
        # strip.text = element_text(size = 14),
        # legend.position = 'none'
    # ) +
    # labs(x = 'Cells', y = 'pEMT score')

# inter_intra_emt_scatterplot <- ggplot(
    # # unique(emt_scores_all_cts[, .(cancer_type_nice, inter_emt, intra_emt)]),
    # unique(
        # emt_scores_all_cts[
            # ,
            # .(
                # cancer_type_nice = cancer_type_nice,
                # inter_emt_scaled = scale(inter_emt),
                # intra_emt_scaled = scale(intra_emt)
            # )
        # ]
    # ),
    # aes(inter_emt_scaled, intra_emt_scaled, colour = cancer_type_nice)
# ) +
    # geom_hline(yintercept = 0, size = 0.3, colour = 'grey', linetype = 'dotted') +
    # geom_vline(xintercept = 0, size = 0.3, colour = 'grey', linetype = 'dotted') +
    # geom_point() +
    # theme_test() +
    # labs(
        # x = 'Inter-tumour pEMT heterogeneity score', # These are medians, but then scaled...
        # y = 'Intra-tumour pEMT heterogeneity score',
        # colour = 'Cancer type'
    # )

inter_intra_emt_profiles <- ggplot(
    emt_scores_all_cts,
    aes(
        pos_frac,
        emt_score,
        group = patient,
        # colour = cancer_type_nice
        colour = scale(intra_emt)
    )
) +
    facet_grid(cols = vars(cancer_type_nice)) +
    geom_line() +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0), breaks = c(-5, 0, 5)) +
    scale_colour_gradientn(colours = colorRamps::blue2green(50)[11:45]) +
    theme(
        panel.background = element_rect(fill = NA, colour = 'black'),
        panel.grid.major.y = element_line(colour = 'grey', size = 0.3, linetype = 'dotted'),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.x = element_text(hjust = 0),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_rect(colour = 'black'),
        strip.text = element_text(size = 13),
        # legend.position = 'bottom',
        legend.position = c(1, -0.15),
        legend.direction = 'horizontal',
        legend.justification = 'right',
        legend.key.height = unit(10, 'pt'),
        legend.title = element_text(margin = margin(0, 2.5, 0, 0)),
        plot.margin = unit(c(5.5, 5.5, 40, 5.5), 'pt')
    ) +
    labs(x = 'Cells', y = 'pEMT score', colour = 'Intra-tumour pEMT heterogeneity score\n')

# Choosing colours:
# set.seed(2030)
# random_colours <- randomcoloR::distinctColorPalette(20)
# This is what we get:
# random_colours <- c('#d31fc1', '#7b32c9', '#3aeaaa', '#94e522', '#ed8ba4', '#e295cb', '#1f7fa5', '#f7da96', '#db5f3d', '#f293da',
					# '#d82295', '#ef97ee', '#e512c2', '#ffccd8', '#99f7b5', '#199e8a', '#a3e22d', '#5eedc0', '#e88db0', '#c2f794')
# See what they look like:
# ggplot(data.table(x = letters[1:20])) +
	# geom_tile(aes(x = x, y = 0, fill = x)) +
	# scale_fill_manual(values = random_colours) +
	# scale_x_discrete(labels = random_colours, expand = c(0, 0)) +
	# scale_y_continuous(expand = c(0, 0)) +
	# theme(
		# axis.text.y = element_blank(),
		# axis.title.y = element_blank(),
		# axis.ticks.y = element_blank(),
		# axis.text.x = element_text(angle = 90, hjust = 1),
		# legend.position = 'none'
	# )

ct_colours = c(
	'Breast' = '#CC8D81',
	'Colorectal' = '#DCE144',
    'Head and Neck' = '#7F9E9B',
    'Liver' = '#D4527A',
    'Lung Adeno.' = '#DF984C',
    'Pancreatic' = '#8975DC'
)

inter_intra_emt_profiles_gtable <- ggplot_gtable(ggplot_build(inter_intra_emt_profiles))

stript <- which(grepl('strip-t', inter_intra_emt_profiles_gtable$layout$name))

for (i in stript) {

    inter_intra_emt_profiles_gtable$grobs[[i]]$grobs[[1]]$children[[
        which(grepl('rect', inter_intra_emt_profiles_gtable$grobs[[i]]$grobs[[1]]$childrenOrder))
    ]]$gp$fill <-
        # This sapply() loop makes the ct_colours lighter:
        sapply(ct_colours, function(x) colorRampPalette(c(x, 'white'))(10)[4])[
            inter_intra_emt_profiles_gtable$grobs[[i]]$grobs[[1]]$children[[
                which(grepl('titleGrob', inter_intra_emt_profiles_gtable$grobs[[i]]$grobs[[1]]$childrenOrder))
            ]]$children[[1]]$label
        ]

}

inter_intra_emt_boxplot_data <- melt(
    unique(
        emt_scores_all_cts[
            ,
            .(
                cancer_type_nice = cancer_type_nice,
                `Inter-tumour pEMT heterogeneity` = scale(inter_emt)[, 1],
                `Intra-tumour pEMT heterogeneity` = scale(intra_emt)[, 1]
            )
        ]
    ),
    id.vars = 'cancer_type_nice'
)

inter_intra_emt_boxplot <- ggplot(
    inter_intra_emt_boxplot_data,
    aes(x = cancer_type_nice, y = value)
) +
    # geom_boxplot(varwidth = TRUE, outlier.shape = NA) +
    # stat_boxplot(geom = 'errorbar', width = 0.2) +
    geom_boxplot(outlier.shape = NA, width = 0.5) +
    geom_jitter(
		data = inter_intra_emt_jitterplot_data,
		aes(x = cancer_type_nice, y = value, colour = cancer_type_nice),
        width = 0.2,
        size = 2
	) +
	scale_colour_manual(values = ct_colours) +
    facet_grid(cols = vars(variable)) +
    theme_bw()

# Add the following to see that the standard deviation really is the same for both x and
# y (it's hard to believe from looking at the extreme y values!):

# xlim(c(-4.5, 4.5)) +
#     ylim(c(-4.5, 4.5)) +
#     geom_hline(yintercept = 0) +
#     geom_vline(xintercept = 0)

pdf('../data_and_figures/inter_intra_emt.pdf', width = 10, height = 4)
ggdraw(inter_intra_emt_profiles_gtable)
print(inter_intra_emt_boxplot)
dev.off()

# pdf('../data_and_figures/inter_intra_emt_profiles.pdf', width = 10, height = 4)
# inter_intra_emt_profiles
# dev.off()

# pdf('../data_and_figures/inter_intra_emt_scatterplot.pdf', width = 6, height = 4)
# inter_intra_emt_scatterplot
# dev.off()

# pdf('../data_and_figures/inter_intra_emt_violin.pdf', width = 16, height = 6)

# ggplot(
    # emt_scores_all_cts,
    # aes(
        # factor(
            # patient_unique,
            # levels = emt_scores_all_cts[
                # ,
                # .(med = median(emt_score)),
                # by = patient_unique
            # ][
                # order(med),
                # patient_unique
            # ]
        # ),
        # emt_score,
        # fill = cancer_type
    # )
# ) +
    # geom_violin(draw_quantiles = 0.5) +
    # theme_test() +
    # theme(axis.text.x = element_blank()) +
    # labs(x = 'Patient', y = 'EMT score')

# dev.off()





# To make a "scheme" figure, showing what the inter- and intra-tumour heterogeneity
# scores mean in terms of the profiles:

# In the following, the colour '#5B8BAC' is half way between skyblue3 and skyblue4, and
# can be obtained by the command colorRampPalette(c('skyblue3', 'skyblue4'))(3)[2].

inter_intra_emt_scheme <- ggplot(
    data.table(
        y = c(
            seq(qnorm(0.01, mean = 0, sd = 1), qnorm(0.99, mean = 0, sd = 1), length.out = 500),
            seq(qnorm(0.01, mean = -1.5, sd = 0.6), qnorm(0.99, mean = -1.5, sd = 0.6), length.out = 500)
        ),
        grp = c(rep(2, 500), rep(1, 500))
    )[
        ,
        x := c(
            pnorm(y[1:500], mean = 0, sd = 1),
            pnorm(y[501:1000], mean = -1.5, sd = 0.6)
        )
    ],
    aes(x = x, y = y)
) +
    # geom_line(aes(group = as.character(grp), alpha = grp), col = 'dodgerblue3') +
    geom_line(aes(group = as.character(grp), colour = as.character(grp))) +
    geom_text_repel(
        aes(x, y, label = l, colour = as.character(grp)),
        data = data.table(x = 0.5, y = qnorm(0.5, mean = -1.5, sd = 0.6), l = 'median'),
        # nudge_x = 0.2,
        # nudge_y = -0.1,
        nudge_x = 0.15,
        nudge_y = -0.4,
        segment.colour = 'darkgrey',
        size = 3,
        colour = 'sienna3'
    ) +
    geom_text_repel(
        aes(x, y, label = l),
        data = data.table(x = 0.5, y = qnorm(0.5), l = 'median'),
        nudge_x = 0.05,
        nudge_y = 0.8,
        segment.colour = 'darkgrey',
        size = 3,
        colour = '#5B8BAC'
    ) +
    geom_text_repel(
        aes(x, y, label = l),
        data = data.table(x = 0.95, y = qnorm(0.95), l = '95th percentile'),
        nudge_x = -0.2,
        nudge_y = 0.4,
        segment.colour = 'darkgrey',
        size = 3,
        colour = '#5B8BAC'
    ) +
    # scale_alpha_continuous(range = c(0.5, 1)) +
    scale_colour_manual(values = c('1' = 'sienna3', '2' = '#5B8BAC')) +
    scale_x_continuous(
        # breaks = c(0.5, 0.95),
        # labels = c('0.5' = '50%', '0.95' = '95%'),
        expand = c(0.01, 0)
    ) +
    scale_y_continuous(
        limits = c(qnorm(0.01, mean = -1.5, sd = 0.6) - 0.2, qnorm(0.99, mean = 0, sd = 1) + 0.2),
        expand = c(0, 0)
    ) +
    geom_segment( # Arrow between profiles to indicate inter-tumour heterogeneity
        aes(
            x = 0.5,
            y = qnorm(0.5, mean = -1.5, sd = 0.6) + 0.075,
            xend = 0.5,
            yend = qnorm(0.5, mean = 0, sd = 1) - 0.075
        ),
        arrow = arrow(ends = 'both', length = unit(5, 'pt'))
    ) +
    geom_segment( # Lower horizontal dashed line
        aes(
            x = 0.4025,
            y = qnorm(0.5, mean = 0, sd = 1),
            xend = 0.4975,
            yend = qnorm(0.5, mean = 0, sd = 1)
        ),
        linetype = 'dashed',
        colour = 'grey',
        size = 0.25
    ) +
    geom_segment( # Upper horizontal dashed line
        aes(
            x = 0.4025,
            y = qnorm(0.95, mean = 0, sd = 1),
            xend = 0.9475,
            yend = qnorm(0.95, mean = 0, sd = 1)
        ),
        linetype = 'dashed',
        colour = 'grey',
        size = 0.25
    ) +
    geom_segment( # Arrow between horizontal dashed lines to indicate intra-tumour heterogeneity
        aes(
            x = 0.4,
            y = qnorm(0.5, mean = 0, sd = 1) + 0.05,
            xend = 0.4,
            yend = qnorm(0.95, mean = 0, sd = 1) - 0.05
        ),
        arrow = arrow(ends = 'both', length = unit(5, 'pt'))
    ) +
    # geom_segment( # Vertical dashed line from lower profile down to 0.5 on x axis
    #     aes(
    #         x = 0.5,
    #         y = qnorm(0.01, mean = -1.5, sd = 0.6) - 0.2,
    #         xend = 0.5,
    #         yend = qnorm(0.5, mean = -1.5, sd = 0.6)
    #     ),
    #     linetype = 'dashed',
    #     colour = 'grey',
    #     size = 0.25
    # ) +
    # geom_segment( # Vertical dashed line from upper profile down to 0.95 on x axis
    #     aes(
    #         x = 0.95,
    #         y = qnorm(0.01, mean = -1.5, sd = 0.6) - 0.2,
    #         xend = 0.95,
    #         yend = qnorm(0.95, mean = 0, sd = 1)
    #     ),
    #     linetype = 'dashed',
    #     colour = 'grey',
    #     size = 0.25
    # ) +
    geom_point(
        aes(x, y, colour = as.character(grp)),
        data = data.table(
            x = c(0.5, 0.5, 0.95),
            y = c(qnorm(0.5, mean = -1.5, sd = 0.6), qnorm(0.5), qnorm(0.95)),
            grp = c(1, 2, 2)
        )
    ) +
    annotate(
        geom = 'text',
        x = 0.38,
        y = qnorm(0.5, mean = 0, sd = 1) + (qnorm(0.95, mean = 0, sd = 1) - qnorm(0.5, mean = 0, sd = 1))/2,
        label = 'Intra-tumour\nheterogeneity',
        hjust = 1,
        size = 4
    ) +
    annotate(
        geom = 'text',
        x = 0.52,
        y = qnorm(0.5, mean = -1.5, sd = 0.6) + 1,
        label = 'Inter-tumour\nheterogeneity',
        hjust = 0,
        size = 4
    ) +
    theme_test() +
    theme(
        # axis.title.x = element_text(hjust = 0, vjust = 6),
        # axis.text.y = element_blank(),
        # axis.ticks.y = element_blank(),
        axis.title.x = element_text(vjust = 6),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = 'none'
    ) +
    labs(x = 'Cells', y = 'pEMT score')
    # labs(x = 'Sorted cells', y = 'pEMT score')

pdf('../data_and_figures/inter_intra_emt_scheme.pdf', width = 3, height = 4)
inter_intra_emt_scheme
dev.off()





# Combination in one aligned plot:

aligned_plots_1 <- align_plots(
    inter_intra_emt_scheme,
    inter_intra_emt_scatterplot,
    align = 'h'
)

aligned_plots_2 <- align_plots(
    inter_intra_emt_scheme,
    inter_intra_emt_profiles,
    align = 'v',
    axis = 'l'
)

aligned_plots_1[[1]]$widths[1] <- aligned_plots_2[[1]]$widths[1]

pdf('../data_and_figures/inter_intra_emt.pdf', width = 10, height = 8.5)

plot_grid(
    blank_plot(),
    plot_grid(
        aligned_plots_1[[1]],
        blank_plot(),
        aligned_plots_1[[2]],
        nrow = 1,
        ncol = 3,
        rel_widths = c(2.3, 0.7, 4)
    ),
    blank_plot(),
    aligned_plots_2[[2]],
    nrow = 4,
    ncol = 1,
    rel_heights = c(0.1, 1.05, 0.05, 1)
) +
    draw_label('A', x = 0, y = 0.975, hjust = -0.1, vjust = 1.1, size = 16, fontface = 2) +
    draw_label('B', x = 0.425, y = 0.975, hjust = -0.1, vjust = 1.1, size = 16, fontface = 2) +
    draw_label('C', x = 0, y = 0.475, hjust = -0.1, vjust = 1.1, size = 16, fontface = 2)

dev.off()





# Trial with just HNSC:

sc_data <- eval(sc_metadata$hnsc$read_quote)

setkey(sc_data, id)

ct_emt_markers <- with(
    deconv_data$`HNSC - Malignant-Basal`,
    head(genes_filtered[ordering], 50)
)

ct_emt_markers <- sc_data[
    cell_type == 'cancer',
    ct_emt_markers[
        apply(.SD, 2, sc_cancer_caf_args$hnsc$scores_filter_fun)
    ],
    .SDcols = ct_emt_markers
]

all_genes_filtered <- sc_data[
    cell_type == 'cancer',
    names(.SD)[
        apply(.SD, 2, sc_cancer_caf_args$hnsc$genes_filter_fun)
    ],
    .SDcols = -c('id', 'patient', 'cell_type')
]

# Using scrabble::score():

emt_scores <- scrabble::score(
    sc_data[cell_type == 'cancer', set_colnames(t(.SD), id), .SDcols = all_genes_filtered],
    list(ct_emt_markers),
    bin.control = TRUE,
    nbin = length(all_genes_filtered) %/% 110
)

emt_scores <- data.table(
    id = rownames(emt_scores),
    patient = sc_data[rownames(emt_scores), patient],
    emt_score = emt_scores[, 1]
)[
    order(patient, emt_score),
    pos_frac := (1:.N)/.N,
    by = patient
]

# I think it's good to exclude at least the top and bottom 1% to eliminate really extreme
# values.  We could even remove more, e.g. top and bottom 2.5%, so we have a sort of "95%
# confidence interval".

emt_score_plot <- ggplot(
    emt_scores[pos_frac > 0.01 & pos_frac < 0.99],
    aes(pos_frac, emt_score, group = patient, colour = as.character(patient))
) +
    scale_colour_manual(values = brewer.pal(12, 'Set3')[-c(2, 11)]) +
    geom_line() +
    theme_test()

# My own implementation of scrabble::score():

# mat <- sc_data[
#     cell_type == 'cancer',
#     set_colnames(t(.SD), id),
#     .SDcols = all_genes_filtered
# ]
#
# gene_averages <- sort(rowMeans(mat))
#
# bins <- setNames(
#     cut(seq_along(gene_averages), breaks = 20, labels = FALSE, include.lowest = TRUE),
#     names(gene_averages)
# )
#
# emt_scores <- rowMeans(
#     sapply(
#         ct_emt_markers,
#         function(g) {
#             mat[g, ] - colMeans(mat[sample(names(bins)[bins == bins[g]], 100), ])
#         },
#         USE.NAMES = TRUE
#     )
# )
#
# emt_scores <- data.table(
#     id = names(emt_scores),
#     patient = sc_data[names(emt_scores), patient],
#     emt_score = emt_scores
# )[
#     order(patient, emt_score),
#     pos_frac := (1:.N)/.N, # .I doesn't work - it goes from 1 to nrow(emt_scores).
#     by = patient
# ]

# Alternative implementation of scrabble::score() using data.table, which is actually slower:

# gene_info <- sc_data[
#     cell_type == 'cancer',
#     .(
#         gene = all_genes_filtered,
#         average = colMeans(.SD)
#     ),
#     .SDcols = all_genes_filtered
# ][
#     order(average),
#     bin := cut(seq_along(average), breaks = 20, labels = FALSE, include.lowest = TRUE)
# ]
#
# emt_scores <- sc_data[
#     cell_type == 'cancer',
#     .(
#         id = id,
#         patient = patient,
#         emt_score = rowMeans(
#             sapply(
#                 ct_emt_markers,
#                 function(g) {
#                     get(g) - rowMeans(
#                         .SD[
#                             ,
#                             gene_info[bin == gene_info[gene == g, bin], sample(gene, 100)],
#                             with = FALSE
#                         ]
#                     )
#                 }
#             )
#         )
#     )
# ][
#     order(patient, emt_score),
#     pos_frac := (1:.N)/.N, # .I doesn't work - it goes from 1 to nrow(emt_scores).
#     by = patient
# ]

# The following calculates the scores per tumour, but I'm not sure I want to do this for
# the current analysis, because I want to preserve inter-tumour heterogeneity.

# emt_scores <- sc_data[
#     cell_type == 'cancer',
#     .(
#         id = id,
#         emt_score = scrabble::score(
#             set_colnames(t(.SD), id),
#             list(ct_emt_markers),
#             bin.control = TRUE,
#             nbin = length(all_genes_filtered) %/% 110
#         )[, 1]
#     ),
#     by = patient,
#     .SDcols = all_genes_filtered
# ][
#     order(patient, emt_score),
#     pos_frac := (1:.N)/.N,
#     by = patient
# ]

# The following also calculates the scores per tumour, but goes one step further - it also
# filters the gene lists on a per-tumour basis, so you have different gene lists for each
# tumour.  This further reduces inter-tumour heterogeneity in the results.

# emt_scores <- sc_data[
#     cell_type == 'cancer',
#     .(
#         id = id,
#         emt_score = scrabble::score(
#             set_colnames(
#                 t(
#                     .SD[
#                         ,
#                         apply(.SD, 2, sc_cancer_caf_args$hnsc$genes_filter_fun),
#                         with = FALSE
#                     ]
#                 ),
#                 id
#             ),
#             list(
#                 .SD[
#                     ,
#                     ct_emt_markers[
#                         apply(.SD, 2, sc_cancer_caf_args$hnsc$scores_filter_fun)
#                     ],
#                     .SDcols = ct_emt_markers
#                 ]
#             ),
#             bin.control = TRUE,
#             nbin = 20,
#             n = 80
#         )[, 1]
#     ),
#     by = patient,
#     .SDcols = -c('id', 'cell_type')
# ][
#     order(patient, emt_score),
#     pos_frac := (1:.N)/.N,
#     by = patient
# ]

# Try to assess inter- and intra-tumour heterogeneity:

inter_intra_emt <- emt_scores[
    ,
    .(
        inter_emt = .SD[pos_frac > 0.25 & pos_frac < 0.75, mean(emt_score)],
        intra_emt = .SD[pos_frac > 0.9 & pos_frac < 0.975, mean(emt_score)] -
            median(emt_score)
    ),
    by = patient
]

inter_intra_plot <- ggplot(
    inter_intra_emt,
    aes(inter_emt, intra_emt, colour = as.character(patient))
) +
    scale_colour_manual(values = brewer.pal(12, 'Set3')[-c(2, 11)]) +
    geom_point() +
    theme_test()

# In the same figure:

ggarrange(
    emt_score_plot,
    inter_intra_plot,
    nrow = 2,
    ncol = 1
)
