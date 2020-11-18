# bsub -q tirosh -R rusage[mem=64000] -o deconv_new_all.o -e deconv_new_all.e Rscript deconv_new_all.R

library(data.table) # 1.12.8
library(ggplot2) # This is version 3.3.1 on my laptop's WSL setup but 3.3.0 on WEXAC.
library(magrittr) # 1.5
library(plyr) # 1.8.6
library(stringr) # 1.4.0
library(egg) # 0.4.5
library(limma) # 3.42.2
library(cowplot) # 1.0.0

source('general_functions.R')
source('tcga_functions.R')

# I need cell_type_markers for the gene filtering:
cell_type_markers <- fread('../../cell_type_markers.csv')

emt_markers <- fread('../../emt_markers.csv')[, gene := alias2SymbolTable(gene)][source != 'GO', sort(unique(gene))]

heatmap_annotations <- c('ACTA2', 'AXL', 'COL1A1', 'COL1A2', 'COL3A1', 'COL4A1', 'COL4A2', 'COL5A1', 'COL5A2', 'COL5A3', 'COL6A1', 'COL6A2',
						 'COL6A3', 'COL7A1', 'CD44', 'CDH2', 'ECM1', 'ECM2', 'FAP', 'FN1', 'IL6', 'ITGA2', 'ITGA5', 'ITGA6', 'ITGB1', 'ITGB3',
						 'ITGB5', 'ITGB6', 'LAMA1', 'LAMA2', 'LAMA3', 'LAMA5', 'LAMB3', 'LAMC1', 'LAMC2', 'MMP1', 'MMP2', 'MMP3', 'MMP10',
						 'MMP14', 'PDPN', 'SNAI1', 'SNAI2', 'SPARC', 'TGFB1', 'TGFBI', 'THY1', 'TNC', 'TWIST1', 'VCAN', 'VIM', 'ZEB1', 'ZEB2')

initial_genes <- c('SNAI1', 'SNAI2', 'TWIST1', 'VIM', 'ZEB1', 'ZEB2')





# Read in data:

expression_data <- fread('../../TCGA_data/tcga_expression_data.csv', key = 'id')
meta_data <- fread('../../TCGA_data/tcga_meta_data.csv', key = 'id')
# Using different version of subtypes data, where I centre and scale the data before finding DE genes:
subtypes_data <- fread('../../TCGA_data/tcga_subtypes_data_centred.csv', key = 'id')
ccle <- fread('../../CCLE_cancer_type_Av.csv', key = 'gene_id')
extra_data <- fread('../data_and_figures/collated_extra_data.csv', key = 'gene')

# I won't need the normal tissue samples, so let's take them out now:

expression_data <- expression_data[meta_data[sample_type != 'normal', id]]
meta_data <- meta_data[sample_type != 'normal']





# Cancer subtypes that have fewer than 30 samples (before inferring subtypes):
# LUSC Primitive: 27 samples
# BLCA Luminal: 26 samples
# BLCA Neuronal: 20 samples
# BRCA Normal-like: 8 samples
# I think if I raised this threshold to 40, I would only lose STAD EBV, which has exactly 30 samples, and maybe splitting up PAAD into Basal and
# Classical would make up for this (in the GI cluster). I'm also leaving out whole cancer types that have fewer than 100 samples, namely ACC, CHOL and
# MESO (also UVM, but I'm not including melanoma any more anyway).

# ACC, CHOL, MESO and UVM have fewer than 100 samples.  But READ only has 105 samples in my dataset (shrinks to 95 after removing normal tissue samples),
# but on GDAC it says there should be 166.  We might want to redo the TCGA data collation.  If READ still has fewer than 100 samples, we could justify
# keeping it in the analysis because it's a subtype of COADREAD.

# Note LUAD Proximal-Inflammatory = Squamoid; Proximal-Proliferative = Magnoid; and Terminal Respiratory Unit (TRU) = Bronchioid.

deconv_args_per_ct <- list(

	blca_basal_squamous = list(
		tcga_cancer_types = 'BLCA',
		subtypes = 'Basal_squamous',
		ref_for_subtypes = 'doi:10.1016/j.cell.2017.09.007',
		ccle_cancer_type = 'urinary tract',
		seed = 1027,
		plot_title = 'BLCA - Basal-squamous'
	),

	blca_luminal = list(
        tcga_cancer_types = 'BLCA',
        subtypes = 'Luminal',
        ref_for_subtypes = 'doi:10.1016/j.cell.2017.09.007',
        ccle_cancer_type = 'urinary tract',
        seed = 9147,
        plot_title = 'BLCA - Luminal'
    ),

    blca_luminal_infiltrated = list(
        tcga_cancer_types = 'BLCA',
        subtypes = 'Luminal_infiltrated',
        ref_for_subtypes = 'doi:10.1016/j.cell.2017.09.007',
        ccle_cancer_type = 'urinary tract',
        seed = 6007,
        plot_title = 'BLCA - Luminal-infiltrated'
    ),

    blca_luminal_papillary = list(
        tcga_cancer_types = 'BLCA',
        subtypes = 'Luminal_papillary',
        ref_for_subtypes = 'doi:10.1016/j.cell.2017.09.007',
        ccle_cancer_type = 'urinary tract',
        seed = 1582,
        plot_title = 'BLCA - Luminal-papillary'
    ),

	blca_neuronal = list(
        tcga_cancer_types = 'BLCA',
        subtypes = 'Neuronal',
        ref_for_subtypes = 'doi:10.1016/j.cell.2017.09.007',
        ccle_cancer_type = 'urinary tract',
        seed = 3939,
        plot_title = 'BLCA - Neuronal'
    ),

    brca_luminal_a = list(
        tcga_cancer_types = 'BRCA',
        subtypes = 'Luminal A',
        ref_for_subtypes = 'doi:10.1038/nature11412',
        ccle_cancer_type = 'breast',
        extra_data_source = 'breast_qian',
        seed = 8216,
        plot_title = 'BRCA - Luminal A'
    ),

    brca_luminal_b = list(
        tcga_cancer_types = 'BRCA',
        subtypes = 'Luminal B',
        ref_for_subtypes = 'doi:10.1038/nature11412',
        ccle_cancer_type = 'breast',
        extra_data_source = 'breast_qian',
        seed = 9056,
        plot_title = 'BRCA - Luminal B'
    ),

    brca_basal_like = list(
        tcga_cancer_types = 'BRCA',
        subtypes = 'Basal-like',
        ref_for_subtypes = 'doi:10.1038/nature11412',
        ccle_cancer_type = 'breast',
        extra_data_source = 'breast_qian',
        seed = 1982,
        plot_title = 'BRCA - Basal-like'
    ),

    brca_her2_enriched = list(
        tcga_cancer_types = 'BRCA',
        subtypes = 'HER2-enriched',
        ref_for_subtypes = 'doi:10.1038/nature11412',
        ccle_cancer_type = 'breast',
        extra_data_source = 'breast_qian',
        seed = 6994,
        plot_title = 'BRCA - HER2-enriched'
    ),

    cesc = list(
        tcga_cancer_types = 'CESC',
        seed = 4825,
        plot_title = 'CESC'
    ),

    coad = list(
        tcga_cancer_types = 'COAD',
        ccle_cancer_type = 'large intestine',
        extra_data_source = 'crc_lee_smc',
        seed = 3260,
        plot_title = 'COAD'
    ),

    esca_ac = list(
        tcga_cancer_types = 'ESCA',
        subtypes = 'AC',
        ref_for_subtypes = 'doi:10.1038/nature20805',
        ccle_cancer_type = 'oesophagus',
        seed = 2022,
        plot_title = 'ESCA - AC'
    ),

    esca_escc = list(
        tcga_cancer_types = 'ESCA',
        subtypes = 'ESCC',
        ref_for_subtypes = 'doi:10.1038/nature20805',
        ccle_cancer_type = 'oesophagus',
        seed = 6358,
        plot_title = 'ESCA - ESCC'
    ),

    hnsc_mesenchymal_basal = list(
        tcga_cancer_types = 'HNSC',
        subtypes = c('Mesenchymal', 'Basal'),
        ref_for_subtypes = 'doi:10.1038/nature14129',
        ccle_cancer_type = 'HNSCC',
        extra_data_source = 'hnscc_puram',
        seed = 7466,
        plot_title = 'HNSC - Mesenchymal & Basal'
    ),

    hnsc_classical = list(
        tcga_cancer_types = 'HNSC',
        subtypes = 'Classical',
        ref_for_subtypes = 'doi:10.1038/nature14129',
        ccle_cancer_type = 'HNSCC',
        extra_data_source = 'hnscc_puram',
        seed = 9998,
        plot_title = 'HNSC - Classical'
    ),

    hnsc_atypical = list(
        tcga_cancer_types = 'HNSC',
        subtypes = 'Atypical',
        ref_for_subtypes = 'doi:10.1038/nature14129',
        ccle_cancer_type = 'HNSCC',
        extra_data_source = 'hnscc_puram',
        seed = 8172,
        plot_title = 'HNSC - Atypical'
    ),

    kich = list(
        tcga_cancer_types = 'KICH',
        ccle_cancer_type = 'kidney',
        seed = 8796,
        plot_title = 'KICH'
    ),

    kirc = list(
        tcga_cancer_types = 'KIRC',
        ccle_cancer_type = 'kidney',
        seed = 247,
        plot_title = 'KIRC'
    ),

    kirp = list(
        tcga_cancer_types = 'KIRP',
        ccle_cancer_type = 'kidney',
        seed = 3031,
        plot_title = 'KIRP'
    ),

    lihc = list(
        tcga_cancer_types = 'LIHC',
        ccle_cancer_type = 'liver',
        extra_data_source = 'liver_ma',
        seed = 9349,
        plot_title = 'LIHC'
    ),

    luad_proximal_inflammatory = list(
        tcga_cancer_types = 'LUAD',
        subtypes = 'Proximal-inflammatory',
        ref_for_subtypes = 'doi:10.1038/nature13385',
        ccle_cancer_type = 'lung',
        extra_data_source = 'luad_kim',
        seed = 9955,
        plot_title = 'LUAD - Proximal-inflammatory'
    ),

    luad_proximal_proliferative = list(
        tcga_cancer_types = 'LUAD',
        subtypes = 'Proximal-proliferative',
        ref_for_subtypes = 'doi:10.1038/nature13385',
        ccle_cancer_type = 'lung',
        extra_data_source = 'luad_kim',
        seed = 743,
        plot_title = 'LUAD - Proximal-proliferative'
    ),

    luad_terminal_respiratory_unit = list(
        tcga_cancer_types = 'LUAD',
        subtypes = 'Terminal respiratory unit',
        ref_for_subtypes = 'doi:10.1038/nature13385',
        ccle_cancer_type = 'lung',
        extra_data_source = 'luad_kim',
        seed = 6026,
        plot_title = 'LUAD - Terminal respiratory unit'
    ),

    lusc_basal = list(
        tcga_cancer_types = 'LUSC',
        subtypes = 'basal',
        ref_for_subtypes = 'doi:10.1038/nature11404',
        ccle_cancer_type = 'lung',
        extra_data_source = 'lusc_qian',
        seed = 2287,
        plot_title = 'LUSC - Basal'
    ),

    lusc_classical = list(
        tcga_cancer_types = 'LUSC',
        subtypes = 'classical',
        ref_for_subtypes = 'doi:10.1038/nature11404',
        ccle_cancer_type = 'lung',
        extra_data_source = 'lusc_qian',
        seed = 7106,
        plot_title = 'LUSC - Classical'
    ),

	lusc_primitive = list(
        tcga_cancer_types = 'LUSC',
        subtypes = 'primitive',
        ref_for_subtypes = 'doi:10.1038/nature11404',
        ccle_cancer_type = 'lung',
        extra_data_source = 'lusc_qian',
        seed = 303,
        plot_title = 'LUSC - Primitive'
    ),

    lusc_secretory = list(
        tcga_cancer_types = 'LUSC',
        subtypes = 'secretory',
        ref_for_subtypes = 'doi:10.1038/nature11404',
        ccle_cancer_type = 'lung',
        extra_data_source = 'lusc_qian',
        seed = 9240,
        plot_title = 'LUSC - Secretory'
    ),

    ov_differentiated = list(
        tcga_cancer_types = 'OV',
        subtypes = 'Differentiated',
        ref_for_subtypes = 'doi:10.1172/JCI65833',
        ccle_cancer_type = 'ovary',
		extra_data_source = 'ovarian_qian',
        seed = 1498,
        plot_title = 'OV - Differentiated'
    ),

    ov_immunoreactive = list(
        tcga_cancer_types = 'OV',
        subtypes = 'Immunoreactive',
        ref_for_subtypes = 'doi:10.1172/JCI65833',
        ccle_cancer_type = 'ovary',
		extra_data_source = 'ovarian_qian',
        seed = 991,
        plot_title = 'OV - Immunoreactive'
    ),

    ov_mesenchymal = list(
        tcga_cancer_types = 'OV',
        subtypes = 'Mesenchymal',
        ref_for_subtypes = 'doi:10.1172/JCI65833',
        ccle_cancer_type = 'ovary',
		extra_data_source = 'ovarian_qian',
        seed = 7549,
        plot_title = 'OV - Mesenchymal'
    ),

    ov_proliferative = list(
        tcga_cancer_types = 'OV',
        subtypes = 'Proliferative',
        ref_for_subtypes = 'doi:10.1172/JCI65833',
        ccle_cancer_type = 'ovary',
		extra_data_source = 'ovarian_qian',
        seed = 8583,
        plot_title = 'OV - Proliferative'
    ),

    paad = list(
        tcga_cancer_types = 'PAAD',
        ccle_cancer_type = 'pancreas',
        extra_data_source = 'pdac_peng',
        seed = 88,
        plot_title = 'PAAD'
    ),

    paad_basal_moffitt = list(
        tcga_cancer_types = 'PAAD',
        subtypes = 'Basal',
        ref_for_subtypes = 'doi:10.1038/ng.3398',
        ccle_cancer_type = 'pancreas',
        extra_data_source = 'pdac_peng',
        seed = 7839,
        plot_title = 'PAAD - Basal'
    ),

    paad_classical_moffitt = list(
        tcga_cancer_types = 'PAAD',
        subtypes = 'Classical',
        ref_for_subtypes = 'doi:10.1038/ng.3398',
        ccle_cancer_type = 'pancreas',
        extra_data_source = 'pdac_peng',
        seed = 7794,
        plot_title = 'PAAD - Classical'
    ),

    prad = list(
        tcga_cancer_types = 'PRAD',
        ccle_cancer_type = 'prostate',
        seed = 6204,
        plot_title = 'PRAD'
    ),

    read = list(
        tcga_cancer_types = 'READ',
        ccle_cancer_type = 'large intestine',
        seed = 3309,
        plot_title = 'READ',
        extra_data_source = 'crc_lee_smc'
    ),

    stad_cin = list(
        tcga_cancer_types = 'STAD',
        subtypes = 'CIN',
        ref_for_subtypes = 'doi:10.1038/nature20805 - STAD',
        ccle_cancer_type = 'stomach',
        seed = 2773,
        plot_title = 'STAD - CIN'
    ),

    stad_ebv = list(
        tcga_cancer_types = 'STAD',
        subtypes = 'EBV',
        ref_for_subtypes = 'doi:10.1038/nature20805 - STAD',
        ccle_cancer_type = 'stomach',
        seed = 3407,
        plot_title = 'STAD - EBV'
    ),

    stad_gs = list(
        tcga_cancer_types = 'STAD',
        subtypes = 'GS',
        ref_for_subtypes = 'doi:10.1038/nature20805 - STAD',
        ccle_cancer_type = 'stomach',
        seed = 5017,
        plot_title = 'STAD - GS'
    ),

    stad_msi = list(
        tcga_cancer_types = 'STAD',
        subtypes = 'MSI',
        ref_for_subtypes = 'doi:10.1038/nature20805 - STAD',
        ccle_cancer_type = 'stomach',
        seed = 6868,
        plot_title = 'STAD - MSI'
    ),

    thca = list(
        tcga_cancer_types = 'THCA',
        ccle_cancer_type = 'thyroid',
        seed = 6798,
        plot_title = 'THCA'
    ),

    ucec = list(
        tcga_cancer_types = 'UCEC',
        ccle_cancer_type = 'endometrium',
        seed = 475,
        plot_title = 'UCEC'
    )

)





# In the below, I'm using inferred subtypes and genes_from_tcga_fun which is slightly different from that used for scRNA-seq analysis.

deconv_data <- sapply(
	names(deconv_args_per_ct),
	function(ct) {
		cat(paste0(ct, '\n'))
		sapply(
			c('90th percentile', 'median', 'none'),
			function(weights_arg) {
				cat(paste0('\t', weights_arg, '\n'))
				setNames(
					lapply(
						c(150, 200, 250),
						function(n) {
							cat(paste0('\t\tn = ', n, '\n'))
							deconv_ct <- do.call(
								deconvolve_emt_caf_data,
								args = c(
									list(
										expression_data = expression_data,
										meta_data = meta_data,
										genes = emt_markers,
										cell_type_markers = cell_type_markers,
										# heatmap_annotations = heatmap_annotations,
										ccle_data = ccle,
										subtypes_data = subtypes_data,
										subtypes_var = 'inferred_subtype',
										extra_data = extra_data,
										# genes_from_tcga_fun = function(x) top_cols_by_fun_cor(x, initial = initial_genes, threshold_fun = function(x) quantile(x, 0.99))$id,
										initial_gene_weights = FALSE,
										genes_filter_fun = function(x) 1:n
									),
									deconv_args_per_ct[[ct]][!(names(deconv_args_per_ct[[ct]]) == 'plot_title')],
									mapvalues(
										weights_arg,
										c('90th percentile', 'median', 'none'),
										list(list(gene_weights_fun = function(x) quantile(x, 0.9)), list(gene_weights_fun = median), list(cell_type_weights = FALSE)),
										warn_missing = FALSE
									) %>% unlist(recursive = FALSE)
								)
							)
							deconv_ct
						}
					),
					c('n = 150', 'n = 200', 'n = 250')
				)
			},
			simplify = FALSE,
			USE.NAMES = TRUE
		)
	},
	simplify = FALSE,
	USE.NAMES = TRUE
)

deconv_data <- unlist(unlist(deconv_data, recursive = FALSE), recursive = FALSE)

names(deconv_data) <- gsub('\\.', ', ', names(deconv_data))

# Reorder genes in deconv results:
deconv_data <- sapply(deconv_data, deconv_reorder, simplify = FALSE, USE.NAMES = TRUE)

deconv_plots <- sapply(
	names(deconv_data),
	function(ct) {
		cat(paste0(ct, '\n'))
		deconv_ct <- deconv_data[[ct]]
		deconv_plot_ct <- do.call(
			deconvolve_emt_caf_plots,
			args = c(
				list(
					data = deconv_ct,
					# Include the following only if you want epithelial scores (takes much longer):
					# expression_data = expression_data,
					heatmap_axis_title = '', # Change heat_map() function so I can put NULL here
					heatmap_legend_title = 'Coexpression',
					heatmap_colours = rev(colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))(50)),
					heatmap_annotations = heatmap_annotations,
					heatmap_annotations_nudge = 0.3,
					purity_colours = rev(colorRampPalette(RColorBrewer::brewer.pal(11, "PuOr"))(50)),
					purity_colour_limits = c(-0.3, 0.3),
					purity_legend_breaks = c(-0.3, 0, 0.3),
					purity_legend_title = 'Correlation\nwith purity',
					purity_legend_direction = 'vertical',
					purity_axis_title = NULL,
					ccle_colours = rev(colorRampPalette(RColorBrewer::brewer.pal(11, 'PiYG')[c(1:3, 6, 9:11)])(50)),
					ccle_fun = function(x) {caTools::runmean(x - mean(x), 30)/max(abs(caTools::runmean(x - mean(x), 30)))},
					ccle_legend_title = 'Tumours vs.\ncell lines',
					ccle_legend_direction = 'vertical',
					ccle_axis_title = NULL,
					extra_colours = colorRampPalette(c('turquoise4', 'turquoise', 'azure', 'gold2', 'gold4'))(50),
					extra_axis_title = NULL,
					extra_legend_title = 'scRNA-seq',
					extra_legend_direction = 'vertical',
					bar_legend_width = NULL,
					bar_legend_height = NULL
				),
				paste(deconv_args_per_ct[[str_split(ct, ',', 2)[[1]][1]]]['plot_title'], str_split(ct, ',', 2)[[1]][2], sep = ',')
			)
		)
		deconv_plot_ct
	},
	simplify = FALSE,
	USE.NAMES = TRUE
)





saveRDS(deconv_data, '../data_and_figures/deconv_new_data_all.rds')
saveRDS(deconv_plots, '../data_and_figures/deconv_new_plots_all.rds')





# Write all figures to a single PDF:
pdf('../data_and_figures/deconv_new_figures_all.pdf', width = 10, height = 12)
for(deconv_plot_ct in deconv_plots) {print(deconv_plot(list(deconv_plot_ct), legends_space = 0.2))}
dev.off()

# Diagnostic figures:
pdf('../data_and_figures/deconv_new_diagnostics_all.pdf', width = 10, height = 12)
i <- 1
for(deconv_plot_ct in deconv_plots) {
	ggarrange(
		plots = c(list(deconv_plot_ct$plots$heatmap), deconv_plot_ct$diagnostics$alternative_purity_cor, deconv_plot_ct$diagnostics$cell_type_bars),
		ncol = 1,
		nrow = 11,
		heights = c(8, rep(0.6, 10)),
		newpage = switch((i == 1) + 1, TRUE, FALSE)
	)
	i <- i + 1
}
dev.off()
