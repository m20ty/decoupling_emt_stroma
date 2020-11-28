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





# We could select which cell types to keep and which to put under 'other' using a fraction of the maximum of the cell type frequencies.  I decided to
# select them manually instead, because there are certain cell types I want to show in some cases, even though they're a bit rare.

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
# Acinar cells seem to express epithelial markers at much lower levels than the ductal cells (maybe this is obvious to a real biologist), and it is
# actually the ductal cells that the cancer originates from (the clue is in the name - PDAC - also see
# https://www.pancreaticcancer.org.uk/information-and-support/facts-about-pancreatic-cancer/types-of-pancreatic-cancer/, where they claim that less
# than 1% of pancreatic cancers are acinar cell cancers).  So we can assume that the acinar cells won't be malignant.  This is similar to alveolar
# cells in the lung, which makes sense because the alveolar duct is similar in form to the acinar duct in pancreas.





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
                    default_figure_widths = list(annotations = 1.5, cancer = 6, caf = 1.2),
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

cairo_pdf('../data_and_figures/sc_heatmaps.pdf', width = 8, height = 5, onefile = TRUE)
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

cairo_pdf('../data_and_figures/sc_heatmaps_lenient.pdf', width = 8, height = 5, onefile = TRUE)
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
                    default_figure_widths = list(annotations = 2.5, cancer = 6, caf = 1.2),
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

cairo_pdf('../data_and_figures/sc_heatmaps_combined.pdf', width = 10, height = 14)

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

cairo_pdf('../data_and_figures/sc_heatmaps_lenient_combined.pdf', width = 10, height = 11)

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

# For this I need the analysis of EMT genes in just the cancer cells.  I hoped not to have to do this - I could reengineer everything so I don't.
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
				theme(
					legend.direction = 'horizontal',
					legend.text = element_text(size = 11),
					legend.key.width = unit(30, 'pt'),
					legend.justification = c(1, 1)
				) +
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
		    x_axis_title = 'Proportion of cell mixture',
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
			error_bars = TRUE
        )
    },
    simplify = FALSE,
    USE.NAMES = TRUE
)

saveRDS(lineplots, '../data_and_figures/simulated_bulk_lineplots.rds')

# Save to PDF with each figure on a separate page:

# pdf('../data_and_figures/simulated_bulk_lineplots.pdf', width = 7, height = 5)
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
    draw_label('Proportion of cell mixture', x = 0.47, y = 0, vjust = -0.5, size = 12) +
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
    draw_label('Proportion of cell mixture', x = 0.47, y = 0, vjust = -0.5, size = 12) +
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
    function(x) rbindlist(lapply(simulated_bulk_data, function(y) as.data.table(y[[x]], keep.rownames = 'id')), fill = TRUE),
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

# Note I tried it with the genes we got from the single cell data and the result was worse.  If you do want to use the genes from single cell data,
# delete initial_gene_weights = FALSE and put genes_filter_fun = NULL and genes_from_tcga_fun = NULL.

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
        heatmap_annotations = c('COL1A2', 'COL5A2', 'ECM1', 'FN1', 'NOTCH2', 'RHOB', 'VCAN', 'VEGFA'),
        plot_title = 'Breast'
    ),
	brca_lenient = list(
        heatmap_annotations = c('COL1A2', 'COL5A2', 'ECM1', 'FN1', 'NOTCH2', 'RHOB', 'VCAN', 'VEGFA'),
        plot_title = 'Breast'
    ),
    coadread = list(
        heatmap_annotations = c('AREG', 'COL1A1', 'COL3A1', 'CXCL1', 'DCN', 'FN1', 'GADD45A', 'PLOD3'),
        plot_title = 'Colorectal'
    ),
	coadread_lenient = list(
        heatmap_annotations = c('AREG', 'COL1A1', 'COL3A1', 'CXCL1', 'DCN', 'FN1', 'GADD45A', 'PLOD3'),
        plot_title = 'Colorectal'
    ),
    hnsc = list(
        heatmap_annotations = c('ACTA2', 'COL1A2', 'COL3A1', 'LAMC2', 'PCOLCE', 'SDC1', 'TGFBI', 'TNC'),
        plot_title = 'Head and Neck'
    ),
    lihc = list(
        heatmap_annotations = c('ACTA2', 'BGN', 'COL1A2', 'LAMC2', 'QSOX1', 'TNFRSF12A', 'VEGFA'),
        plot_title = 'Liver'
    ),
	lihc_lenient = list(
        heatmap_annotations = c('ACTA2', 'BGN', 'COL1A2', 'LAMC2', 'QSOX1', 'THY1', 'TNFRSF12A', 'VEGFA'),
        plot_title = 'Liver'
    ),
    luad = list(
        heatmap_annotations = c('CALU', 'COL1A2', 'COL3A1', 'MMP2', 'QSOX1', 'SDC4', 'THY1', 'VEGFA'),
        plot_title = 'Lung Adenocarcinoma'
    ),
	# luad_lenient = list(
        # heatmap_annotations = c('ACTA2', 'COL1A1', 'COL1A2', 'ITGB1', 'LAMC2', 'QSOX1', 'RHOB', 'THY1'),
        # plot_title = 'Lung Adenocarcinoma'
    # ),
    lusc = list( # Need to change annotations from here downwards...
        heatmap_annotations = c('COL1A2', 'COL3A1', 'FAP', 'IGFBP2', 'PFN2', 'SDC1', 'THY1', 'TNC'),
        plot_title = 'Lung squamous'
    ),
	lusc_lenient = list(
        heatmap_annotations = c('COL1A2', 'COL3A1', 'FAP', 'IGFBP2', 'PFN2', 'SDC1', 'THY1', 'TNC'),
        plot_title = 'Lung squamous'
    ),
    ov = list(
        heatmap_annotations = c('ACTA2', 'CDH2', 'COL3A1', 'FAP', 'LAMC1', 'MMP1', 'PFN2', 'THY1'),
        plot_title = 'Ovarian'
    ),
	ov_lenient = list(
        heatmap_annotations = c('ACTA2', 'CDH2', 'COL3A1', 'FAP', 'LAMC1', 'MMP1', 'PFN2', 'THY1'),
        plot_title = 'Ovarian'
    ),
    paad = list(
        heatmap_annotations = c('COL1A2', 'COL3A1', 'DCN', 'FAP', 'LAMC2', 'QSOX1', 'SDC4', 'VEGFA'),
        plot_title = 'Pancreatic'
    )
)

for(ct in names(deconv_plot_args_per_ct)) {
	deconv_plot_args_per_ct[[ct]]$heatmap_annotations <- c(
		deconv_plot_args_per_ct[[ct]]$heatmap_annotations,
		'SNAI1', 'SNAI2', 'TWIST1', 'VIM', 'ZEB1', 'ZEB2'
	)
}

# In the following, I am not opting to include epithelial cells in the diagnostics, mainly to save time (it saves a lot of time), but also because
# the epithelial cell correlations are (thankfully) as I would expect, and because conceptually I don't think it adds much.

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
					heatmap_annotations_side = 'left',
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

# The following won't work unless we remove heatmap_annotations_side = 'left' from the above.

# cairo_pdf('../data_and_figures/simulated_deconv_figures.pdf', width = 7, height = 6, onefile = TRUE)
# for(ct in names(simulated_deconv_plots)[!grepl('lenient', names(simulated_deconv_plots))]) {
# 	deconv_plot(
# 		simulated_deconv_plots[ct],
# 		legends_rel_size = c(0.5, 0.75, 0, 0.75, 0, 0.75, 0.75, 3, 2),
# 		legends_space = 0.7,
# 		title.position = 'right'
# 	) %>% print
# }
# dev.off()
#
# cairo_pdf('../data_and_figures/simulated_deconv_figures_lenient.pdf', width = 7, height = 6, onefile = TRUE)
# for(ct in names(simulated_deconv_plots)[grepl('lenient', names(simulated_deconv_plots))]) {
# 	deconv_plot(
# 		simulated_deconv_plots[ct],
# 		legends_rel_size = c(0.5, 0.75, 0, 0.75, 0, 0.75, 0.75, 3, 2),
# 		legends_space = 0.7,
# 		title.position = 'right'
# 	) %>% print
# }
# dev.off()

# Diagnostic figures:
# pdf('../data_and_figures/simulated_deconv_diagnostics.pdf', width = 10, height = 16.8)
# i <- 1
# for(deconv_ct in simulated_deconv_plots) {
#     ggarrange(
#         plots = c(list(deconv_ct$plots$heatmap), deconv_ct$diagnostics$alternative_purity_cor, deconv_ct$diagnostics$cell_type_bars),
#         ncol = 1,
#         nrow = 11,
#         heights = c(8, rep(0.8, 10)),
#         newpage = switch((i == 1) + 1, TRUE, FALSE)
#     )
#     i <- i + 1
# }
# dev.off()

# This is a bit worrying because some analyses fail the criterion I used for filtering the TCGA deconvs (regression slope < -0.1):
# diagnostic_cell_type_lms <- sapply(
#     simulated_deconvs,
#     function(deconv_ct) {
#         with(
#             deconv_ct,
#             cor_with_initial_and_cell_types[genes_filtered][
#                 ordering,
#                 sapply(.SD, function(ct) setNames(lm(ct ~ I(.I/.N))$coeff['I(.I/.N)'], NULL), USE.NAMES = TRUE),
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

				# plot_data <- copy(sc_data[cell_type == x])[
				# 	,
				# 	complexity := apply(.SD, 1, function(x) sum(x > 0)),
				# 	.SDcols = -c('id', 'patient', 'cell_type')
				# ][
				# 	,
				# 	.SD[sample(1:.N, sample_sizes[[ct]], prob = complexity)]
				# ][, complexity := NULL][, c('id', simulated_deconvs[[ct]]$genes_filtered), with = FALSE]

				plot_data <- sc_data[cell_type == x, c('id', simulated_deconvs[[ct]]$genes_filtered), with = FALSE]

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

		sc_deconv_comp_data <- rbind(sc_deconv_comp$cancer$plot_data[, cell_type := 'cancer'], sc_deconv_comp$caf$plot_data[, cell_type := 'caf'])

		# Calculate averages in cancer cells and CAFs separately, so we can check them afterwards:
		gene_averages_cancer_caf <- sc_deconv_comp_data[, .(ave_exp = mean(expression_level)), by = .(gene, cell_type)]

		# To centre genes w.r.t. the average of the averages of cancer and CAF:
		gene_averages <- sc_deconv_comp_data[
			,
			.(ave_exp = mean(c(mean(expression_level[cell_type == 'cancer']), mean(expression_level[cell_type == 'caf'])))),
			by = .(symbol = gene)
		]
		sc_deconv_comp_data[, expression_level := expression_level - gene_averages[symbol == unique(gene), ave_exp], by = gene]

		# To centre the cells as well:
		sc_deconv_comp_data[, expression_level_cc := expression_level - mean(expression_level), by = id]

		# Apply running average to each cell:
		sc_deconv_comp_data[
			,
			expression_level_cc_rm := setNames(
				runmean(setNames(expression_level_cc, gene)[with(simulated_deconvs[[ct]], genes_filtered[ordering])], 30),
				with(simulated_deconvs[[ct]], genes_filtered[ordering])
			)[as.character(gene)], # Melting made gene a factor, and this reordering doesn't work unless we coerce it to character
			by = id
		]

		sc_heatmaps <- sapply(
			c('expression_level', 'expression_level_cc', 'expression_level_cc_rm'),
			function(expr_var) sapply(
				c('cancer', 'caf'),
				function(x) {
					ggplot(
						sc_deconv_comp_data[cell_type == x],
						aes(
							x = factor(gene, levels = with(simulated_deconvs[[ct]], genes_filtered[ordering])),
							y = factor(id, levels = sc_deconv_comp[[x]]$ordered_cell_ids),
							fill = get(expr_var)
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
							limits = switch((expr_var == 'expression_level_cc_rm') + 1, c(-4, 4), c(-2, 2)),
							breaks = switch((expr_var == 'expression_level_cc_rm') + 1, c(-4, -2, 0, 2, 4), c(-2, -1, 0, 1, 2)),
							labels = switch(
								(expr_var == 'expression_level_cc_rm') + 1,
								c('-4' = '\u2264 -4', '-2' = '-2', '0' = '0', '2' = '2', '4' = '\u2265 4'),
								c('-2' = '\u2264 -2', '-1' = '-1', '0' = '0', '1' = '1', '2' = '\u2265 2')
							),
							oob = squish
						) +
						theme(
							axis.text = element_blank(),
							axis.title.x = element_blank(),
							axis.ticks = element_blank(),
							axis.ticks.length = unit(0, 'pt'),
							plot.margin = unit(c(5.5, 5.5, 1.25, 0), 'pt')
							# plot.margin = switch((x == 'cancer') + 1, unit(c(1.25, 5.5, 5.5, 5.5), 'pt'), unit(c(5.5, 5.5, 1.25, 5.5), 'pt'))
						) +
						labs(y = mapvalues(x, c('cancer', 'caf'), c('Cancer', 'CAF'), warn_missing = FALSE), fill = 'Relative\nexpression\nlevel')
				},
				simplify = FALSE,
				USE.NAMES = TRUE
			),
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

		# Re-centre genes and cells w.r.t. the filtered gene list:
		gene_averages <- plot_data[
			,
			.(ave_exp = mean(c(mean(expression_level[cell_type == 'cancer']), mean(expression_level[cell_type == 'caf'])))),
			by = .(symbol = gene)
		]
		plot_data[, expression_level := expression_level - gene_averages[symbol == unique(gene), ave_exp], by = gene]
        plot_data[, expression_level_cc := expression_level - mean(expression_level), by = id]

		# Re-do running average per cell:
		plot_data[
			,
			expression_level_cc_rm := setNames(
				runmean(setNames(expression_level_cc, gene)[with(simulated_deconv_filtered, genes_filtered[ordering])], 30),
				with(simulated_deconv_filtered, genes_filtered[ordering])
			)[as.character(gene)], # Melting made gene a factor, and this reordering doesn't work unless we coerce it to character
			by = id
		]

		plot_data <- sapply(
			c('cancer', 'caf'),
			function(x) {
				plot_data[
					cell_type == x,
					.(
						id = factor(id, levels = sc_deconv_comp[[x]]$ordered_cell_ids),
						gene = factor(gene, levels = with(simulated_deconv_filtered, genes_filtered[ordering])),
						expression_level = expression_level,
						expression_level_cc = expression_level_cc,
						expression_level_cc_rm = expression_level_cc_rm
					)
				]
			},
			simplify = FALSE,
			USE.NAMES = TRUE
		)

		sc_heatmaps_filtered <- sapply(
			c('expression_level', 'expression_level_cc', 'expression_level_cc_rm'),
			function(expr_var) sapply(
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
							fill = get(expr_var)
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
							limits = switch((expr_var == 'expression_level_cc_rm') + 1, c(-4, 4), c(-2, 2)),
							breaks = switch((expr_var == 'expression_level_cc_rm') + 1, c(-4, -2, 0, 2, 4), c(-2, -1, 0, 1, 2)),
							labels = switch(
								(expr_var == 'expression_level_cc_rm') + 1,
								c('-4' = '\u2264 -4', '-2' = '-2', '0' = '0', '2' = '2', '4' = '\u2265 4'),
								c('-2' = '\u2264 -2', '-1' = '-1', '0' = '0', '1' = '1', '2' = '\u2265 2')
							),
							# limits = c(-4, 4),
							# breaks = c(-4, -2, 0, 2, 4),
							# labels = c('-4' = '\u2264 -4', '-2' = '-2', '0' = '0', '2' = '2', '4' = '\u2265 4'),
							oob = squish
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
			),
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
			ordered_cell_ids = sapply(sc_deconv_comp, `[[`, 'ordered_cell_ids', simplify = FALSE, USE.NAMES = TRUE)
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
				lapply(sc_sim_deconv_comp[[ct]]$sc_heatmaps_unfiltered$expression_level_cc_rm, function(x) x + theme(legend.position = 'none'))
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
							theme(
								legend.justification = c(0, 1),
								legend.direction = 'vertical',
								legend.key.width = unit(10, 'pt'),
								legend.key.height = unit(10, 'pt')
							) +
							labs(fill = 'Correlation\nwith purity')
					),
					get_legend(
						simulated_deconv_plots[[1]]$plots$ccle_bar +
							theme(
								legend.justification = c(0, 1),
								legend.direction = 'vertical',
								legend.key.width = unit(10, 'pt'),
								legend.key.height = unit(10, 'pt')
							) +
							labs(fill = 'Tumours vs.\ncell lines')
					),
					get_legend(
						simulated_deconv_plots[[1]]$plots$heatmap + theme(legend.justification = c(0, 1))
					),
					get_legend(sc_sim_deconv_comp[[1]]$sc_heatmaps_unfiltered$expression_level_cc_rm[[1]] + theme(legend.justification = c(0, 1)))
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
				lapply(sc_sim_deconv_comp[[ct]]$sc_heatmaps_unfiltered$expression_level_cc_rm, function(x) x + theme(legend.position = 'none'))
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
							theme(
								legend.justification = c(0, 1),
								legend.direction = 'vertical',
								legend.key.width = unit(10, 'pt'),
								legend.key.height = unit(10, 'pt')
							) +
							labs(fill = 'Correlation\nwith purity')
					),
					get_legend(
						simulated_deconv_plots[[1]]$plots$ccle_bar +
							theme(
								legend.justification = c(0, 1),
								legend.direction = 'vertical',
								legend.key.width = unit(10, 'pt'),
								legend.key.height = unit(10, 'pt')
							) +
							labs(fill = 'Tumours vs.\ncell lines')
					),
					get_legend(
						simulated_deconv_plots[[1]]$plots$heatmap + theme(legend.justification = c(0, 1))
					),
					get_legend(sc_sim_deconv_comp[[1]]$sc_heatmaps_unfiltered$expression_level_cc_rm[[1]] + theme(legend.justification = c(0, 1)))
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
			sc_sim_deconv_comp[[ct]]$sc_heatmaps_filtered$expression_level_cc_rm
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
							theme(legend.justification = c(0, 1), legend.direction = 'vertical') +
							labs(fill = 'Correlation\nwith purity')
					),
					get_legend(
						simulated_deconv_plots[[1]]$plots$ccle_bar +
							theme(legend.justification = c(0, 1), legend.direction = 'vertical') +
							labs(fill = 'Tumours vs.\ncell lines')
					),
					get_legend(simulated_deconv_plots[[1]]$plots$heatmap + theme(legend.justification = c(0, 1))),
					get_legend(sc_sim_deconv_comp[[1]]$sc_heatmaps_filtered$expression_level_cc_rm[[1]] + theme(legend.justification = c(0, 1)))
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
			sc_sim_deconv_comp[[ct]]$sc_heatmaps_filtered$expression_level_cc_rm
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
							theme(legend.justification = c(0, 1), legend.direction = 'vertical') +
							labs(fill = 'Correlation\nwith purity')
					),
					get_legend(
						simulated_deconv_plots[[1]]$plots$ccle_bar +
							theme(legend.justification = c(0, 1), legend.direction = 'vertical') +
							labs(fill = 'Tumours vs.\ncell lines')
					),
					get_legend(simulated_deconv_plots[[1]]$plots$heatmap + theme(legend.justification = c(0, 1))),
					get_legend(sc_sim_deconv_comp[[1]]$sc_heatmaps_filtered$expression_level_cc_rm[[1]] + theme(legend.justification = c(0, 1)))
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

cairo_pdf('../data_and_figures/final_figures_resubmission/2.pdf', width = 13, height = 11.4)

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
							c(simulated_deconv_plots[[ct]]$plots, sc_sim_deconv_comp[[ct]]$sc_heatmaps_unfiltered$expression_level_cc_rm),
							function(x) x + theme(legend.position = 'none', plot.title = element_blank(), axis.title.y = element_blank()),
							simplify = FALSE,
							USE.NAMES = TRUE
						)
						sc_sim_deconv_comp_figures$heatmap <- sc_sim_deconv_comp_figures$heatmap +
							theme(plot.margin = unit(c(0, 5.5, 1.25, 0), 'pt'))
						sc_sim_deconv_comp_figures$axis_labels <- sc_sim_deconv_comp_figures$axis_labels +
							theme(plot.margin = unit(c(0, 0, 1.25, 0), 'pt'))
						# sc_sim_deconv_comp_figures$sc_heatmaps_unfiltered <- sapply(
							# sc_sim_deconv_comp_figures$sc_heatmaps_unfiltered,
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
								rel_heights = c(4, 15, 11)
							),
							plot_grid(
								blank_plot(),
								sc_sim_deconv_comp_figures$axis_labels,
								ggplot(data.frame(x = 0:1, y = 0:1), aes(x = x, y = y)) +
									geom_point(size = 0, colour = 'white') +
									geom_segment(
										aes(x = 0.65, xend = 0.65, y = 0.05, yend = 0.95),
										arrow = arrow(ends = 'both', length = unit(5, 'pt'))
									) +
									scale_x_continuous(expand = c(0, 0)) +
									scale_y_continuous(expand = c(0, 0)) +
									theme(
							            axis.text = element_blank(),
										axis.title.x = element_blank(),
										axis.title.y = element_text(vjust = -3),
							            axis.ticks = element_blank(),
							            axis.ticks.length = unit(0, 'pt'),
							            plot.background = element_blank(),
							            panel.background = element_blank(),
										plot.margin = unit(c(5.5, 0, 1.25, 5.5), 'pt')
							        ) +
									labs(y = paste0('Cancer\ncells\n(n = ', length(sc_sim_deconv_comp[[ct]]$ordered_cell_ids$cancer), ')')),
								ggplot(data.frame(x = 0:1, y = 0:1), aes(x = x, y = y)) +
									geom_point(size = 0, colour = 'white') +
									geom_segment(
										aes(x = 0.65, xend = 0.65, y = 0.05, yend = 0.95),
										arrow = arrow(ends = 'both', length = unit(5, 'pt'))
									) +
									scale_x_continuous(expand = c(0, 0)) +
									scale_y_continuous(expand = c(0, 0)) +
									theme(
							            axis.text = element_blank(),
										axis.title.x = element_blank(),
										axis.title.y = element_text(vjust = -3),
							            axis.ticks = element_blank(),
							            axis.ticks.length = unit(0, 'pt'),
							            plot.background = element_blank(),
							            panel.background = element_blank(),
										plot.margin = unit(c(5.5, 0, 1.25, 5.5), 'pt')
							        ) +
									labs(y = paste0('\nCAFs\n(n = ', length(sc_sim_deconv_comp[[ct]]$ordered_cell_ids$caf), ')')),
								# blank_plot() +
								# 	labs(y = 'Cancer\ncells') +
								# 	scale_y_continuous(position = 'right') +
								# 	theme(plot.margin = unit(c(5.5, 5.5, 1.25, 5.5), 'pt'), axis.title.y.right = element_text(angle = 90)),
								# blank_plot() +
								# 	labs(y = 'CAFs') +
								# 	scale_y_continuous(position = 'right') +
								# 	theme(plot.margin = unit(c(5.5, 5.5, 1.25, 5.5), 'pt'), axis.title.y.right = element_text(angle = 90)),
								blank_plot(),
								nrow = 5,
								ncol = 1,
								rel_heights = c(4, 15, 5, 5, 1)
							),
							plot_grid(
								plotlist = c(
									list(
										blank_plot() +
										theme(plot.margin = unit(c(11, 5.5, 5.5, 0), 'pt')) +
										labs(title = deconv_plot_args_per_ct[[ct]]$plot_title)
									),
									sc_sim_deconv_comp_figures[c('purity_bar', 'ccle_bar', 'heatmap', 'cancer', 'caf')],
									list(
										blank_plot() +
											scale_x_continuous(position = 'top') +
											theme(axis.title.x = element_text(), plot.margin = unit(c(4.5, 0, 0, 0), 'pt')) +
											labs(x = 'Genes')
									)
								),
								nrow = 7,
								ncol = 1,
								align = 'v',
								rel_heights = c(2, 1, 1, 15, 5, 5, 1)
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
					simulated_deconv_plots[[ct]]$plots$purity_bar +
						theme(legend.justification = c(1, 1), legend.direction = 'horizontal', legend.key.width = NULL) +
						labs(fill = 'Correlation with purity\n')
				),
				blank_plot(),
				get_legend(
					simulated_deconv_plots[[ct]]$plots$heatmap +
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
					simulated_deconv_plots[[ct]]$plots$ccle_bar +
						theme(legend.justification = c(1, 1), legend.direction = 'horizontal', legend.key.width = NULL) +
						labs(fill = 'Tumours vs. cell lines\n')
				),
				blank_plot(),
				get_legend(
					sc_sim_deconv_comp[[ct]]$sc_heatmaps_unfiltered$expression_level_cc_rm[[1]] +
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
		rel_heights = c(0.125, 0.95, 0.125, 1.8)
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





# Make plots to compare EMT and CAF genes in single cell analysis vs simulated bulk deconvolution and vs TCGA bulk deconvolution:

simulated_deconvs <- readRDS('../data_and_figures/simulated_deconvs.rds')
simulated_bulk_data <- fread('../data_and_figures/simulated_bulk_data.csv', key = 'id')
deconv_data <- readRDS('../data_and_figures/deconv_data.rds')

ct_to_keep <- c('blca_luminal_papillary', 'blca_basal_squamous', 'brca_luminal_a', 'brca_luminal_b', 'brca_basal_like', 'brca_her2_enriched',
    'coad', 'esca_ac', 'hnsc_mesenchymal_basal', 'hnsc_classical', 'luad_proximal_inflammatory', 'luad_proximal_proliferative',
    'luad_terminal_respiratory_unit', 'lusc_classical', 'lusc_secretory', 'ov_differentiated', 'ov_immunoreactive', 'ov_proliferative', 'paad',
    'read', 'stad_cin', 'stad_msi', 'ucec')
nice_names_for_figure <- c('BLCA - Luminal-Papillary', 'BLCA - Basal-Squamous', 'BRCA - Luminal A', 'BRCA - Luminal B', 'BRCA - Basal-like',
    'BRCA - HER2-enriched', 'COAD', 'ESCA - Adenocarcinoma', 'HNSC - Malignant-Basal', 'HNSC - Classical', 'LUAD - Squamoid', 'LUAD - Magnoid',
    'LUAD - Bronchioid', 'LUSC - Classical', 'LUSC - Secretory', 'OV - Differentiated', 'OV - Immunoreactive', 'OV - Proliferative', 'PAAD', 'READ',
    'STAD - CIN', 'STAD - MSI', 'UCEC')

# The following doesn't make sense any more...
sc_to_bulk_names <- list(
	brca = c('BRCA - Luminal A', 'BRCA - Luminal B', 'BRCA - Basal-like', 'BRCA - HER2-enriched'),
    coadread = c('COADREAD - Colon', 'COADREAD - Rectum'),
    hnsc = c('HNSC - Malignant-Basal', 'HNSC - Classical', 'HNSC - Atypical'),
    lihc = 'LIHC',
    lung = c('LUAD - Squamoid', 'LUAD - Magnoid', 'LUSC - Basal', 'LUSC - Classical', 'LUSC - Primitive', 'LUSC - Secretory'),
    paad = 'PAAD'
)

# Get genes from single cell data and TCGA deconv data:
genes_sc <- unique(unlist(lapply(sc_cancer_caf, `[[`, 'genes_filtered')))
genes_tcga <- unique(unlist(lapply(deconv_data[unlist(sc_to_bulk_names)], `[[`, 'genes_filtered')))
genes_simulated <- unique(unlist(lapply(simulated_deconvs, `[[`, 'genes_filtered')))

# Calculate scores from transformed bulk (TCGA and simulated) expression data:
scores_data_tcga <- deconv_scores(
    expression_data,
    deconv_data[unlist(sc_to_bulk_names)],
    scale_fun = function(x) x/(3*sd(x[!is.na(x)])),
    transform_data = TRUE,
    additional_genes = unique(c(genes_sc[genes_sc %in% names(expression_data)], genes_simulated[genes_simulated %in% names(expression_data)]))
)

# There are some NAs in the following, which we allowed thanks to the altered scale_fun which ignores the NAs during the scaling.  Perhaps we should
# remove these genes altogether.

# Note that since the simulated bulk data was defined from the single cell data, all the genes in genes_sc will be in names(simulated_bulk_data).

scores_data_simulated <- deconv_scores(
    simulated_bulk_data,
    simulated_deconvs,
    scale_fun = function(x) x/(3*sd(x[!is.na(x)])),
    transform_data = TRUE,
    additional_genes = unique(c(genes_sc, genes_tcga[genes_tcga %in% names(simulated_bulk_data)]))
)

scores <- sapply(
    names(sc_metadata),
    function(ct) {

        sc_data <- eval(sc_metadata[[ct]]$read_quote)

        all_genes <- unique(
            c(
                sc_cancer_caf[[ct]]$genes_filtered,
                simulated_deconvs[[ct]]$genes_filtered,
                unlist(lapply(deconv_data[sc_to_bulk_names[[ct]]], `[[`, 'genes_filtered'))
            )
        )

        scores_tcga <- scores_data_tcga[all_genes, setNames(rowMeans(.SD), gene), .SDcols = sc_to_bulk_names[[ct]]]

        scores_simulated <- scores_data_simulated[all_genes, setNames(get(ct), gene)]

        # In the following, would it be sensible to reverse the log before taking means?

        # I was going to take the top and bottom 20 from the sorted vector of the following gene scores, then score genes by correlation with these
		# head and tail genes.  But it's probably enough to use these scores as they are (or after log, to make them normally distributed).

        scores_sc <- sc_data[cell_type == 'cancer', colMeans(.SD), .SDcols = all_genes[all_genes %in% names(sc_data)]]/
			sc_data[cell_type == 'fibroblast', colMeans(.SD), .SDcols = all_genes[all_genes %in% names(sc_data)]]

        out <- data.table(gene = all_genes, sc_score = scores_sc[all_genes], tcga_score = scores_tcga, simulated_score = scores_simulated)

        out[, sc_score := log10(sc_score)][
            ,
            c('sc_score', 'tcga_score', 'simulated_score') := lapply(
                .SD,
                function(x) {
                    sapply(
                        x,
                        function(y) switch(
                            (!is.na(y) & !is.infinite(y)) + 1,
                            NA,
                            (y - mean(x[!is.na(x) & !is.infinite(x)]))/sd(x[!is.na(x) & !is.infinite(x)])
                        )
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
                                c(
									5.5,
									5.5,
									switch((i %in% (n - n_col + 1):n) + 1, 5.5, 20),
									switch((i %in% ((0:((n / n_col) - 1))*n_col + 1)) + 1, 5.5, 20)
								),
                                'pt'
                            )
                        )
                }
            ),
            nrow = n_row,
            ncol = n_col,
            # Don't need align = 'h' any more because no plots that aren't on the bottom row will have x axis text, and aligning introduces space
			# where the x axis text would go if it was there...
            rel_widths = rel_widths,
            rel_heights = rel_heights
        ) +
            draw_label(var_axis_titles[1], x = 0.5, y = 0, vjust = -0.5, size = 12) +
            draw_label(var_axis_titles[2], x = 0, y = 0.5, vjust = 1.3, angle = 90, size = 12)

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
                                axis.text.x = switch(
									(i %in% (n_col*(n %/% n_col) - (n_col - (n %% n_col)) + 1):(n_col*(n %/% n_col))) + 1,
									element_blank(),
									NULL
								),
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
            rel_heights = c(sum(rel_heights[-length(rel_heights)]), rel_heights[length(rel_heights)])
        ) +
            draw_label(var_axis_titles[1], x = 0.5, y = 0, vjust = -0.5, size = 12) +
            draw_label(var_axis_titles[2], x = 0, y = 0.5, vjust = 1.3, angle = 90, size = 12)

    }

}





# Write plots to PDF:

pdf('../data_and_figures/signature_concordance.pdf', width = 12, height = 6)

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

sc_genes <- sapply(sc_cancer_caf, `[[`, 'genes_filtered', simplify = FALSE, USE.NAMES = TRUE)
simulated_genes <- sapply(simulated_deconvs, `[[`, 'genes_filtered', simplify = FALSE, USE.NAMES = TRUE)
tcga_genes <- sapply(
    names(sc_to_bulk_names),
    function(ct) unique(unlist(lapply(deconv_data[sc_to_bulk_names[[ct]]], `[[`, 'genes_filtered'))),
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

pdf('../data_and_figures/signature_concordance_sc_tcga.pdf', width = 6, height = 9)
plot_var_pair(
    c('sc_score', 'tcga_score'),
    scores[c('hnsc', 'lung', 'paad', 'lihc', 'tnbc')],
    n_row = 3,
    n_col = 2,
    sc_genes = sc_genes,
    tcga_genes = tcga_genes,
    collate_genes_fun = union,
    plot_titles = sapply(sc_cancer_caf_heatmaps_args[c('hnsc', 'lung', 'paad', 'lihc', 'tnbc')], `[[`, 'annotations_title'),
    var_axis_titles = c('scRNA-seq score', 'TCGA deconvolution score'),
    rel_widths = c(1.115, 1),
    rel_heights = c(1, 1, 1.075)
)
dev.off()

pdf('../data_and_figures/signature_concordance_sc_tcga_supp.pdf', width = 6, height = 3.2)
plot_var_pair(
    c('sc_score', 'tcga_score'),
    scores[c('brca', 'coadread')],
    n_row = 1,
    n_col = 2,
    sc_genes = sc_genes,
    tcga_genes = tcga_genes,
    collate_genes_fun = union,
    plot_titles = sapply(sc_cancer_caf_heatmaps_args[c('brca', 'coadread')], `[[`, 'annotations_title'),
    var_axis_titles = c('scRNA-seq score', 'TCGA deconvolution score'),
    rel_widths = c(1.115, 1)
)
dev.off()
