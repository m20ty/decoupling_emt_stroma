# bsub -q tirosh -n 8 -R rusage[mem=8000] -o sc_deconv_comp.o -e sc_deconv_comp.e Rscript sc_deconv_comp.R

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
library(plyr) # 1.8.6
library(latex2exp) # 0.4.0
library(ggrepel) # 0.8.2
library(egg) # 0.4.5

source('general_functions.R')
source('tcga_functions.R')
source('sc_functions.R')





expression_data <- fread('../../TCGA_data/tcga_expression_data.csv', key = 'id')
meta_data <- fread('../../TCGA_data/tcga_meta_data.csv', key = 'id')

expression_data <- expression_data[meta_data[sample_type != 'normal', id]]
meta_data <- meta_data[sample_type != 'normal']

emt_markers <- fread('../../emt_markers.csv')[, gene := alias2SymbolTable(gene)][source != 'GO', sort(unique(gene))]





sc_metadata <- list(
	brca = list(
		tcga_cancer_types = 'BRCA',
		read_quote = quote(fread('../data_and_figures/qian_breast_2020_reclassified.csv')[cell_type != 'ambiguous' & id != 'sc5rJUQ064_CCATGTCCATCCCATC', -c('cell_type_author', 'cell_type_lenient')]),
		initial_cell_types = c('b_cell', 'caf', 'cancer', 'endothelial', 'macrophage', 'mast', 't_cell'),
		rare_cell_types = 'dendritic'
	),
	brca_lenient = list(
		tcga_cancer_types = 'BRCA',
		read_quote = quote(fread('../data_and_figures/qian_breast_2020_reclassified.csv')[cell_type_lenient != 'ambiguous' & id != 'sc5rJUQ064_CCATGTCCATCCCATC', -c('cell_type_author', 'cell_type')] %>% setnames('cell_type_lenient', 'cell_type')),
		initial_cell_types = c('b_cell', 'caf', 'cancer', 'endothelial', 'macrophage', 'mast', 't_cell'),
		rare_cell_types = 'dendritic'
	),
	coadread = list(
		tcga_cancer_types = c('COAD', 'READ'),
		read_quote = quote(fread('../data_and_figures/lee_crc_2020_smc_reclassified.csv')[cell_type != 'ambiguous', -c('cell_type_author', 'cell_type_lenient')]),
		initial_cell_types = c('b_cell', 'caf', 'cancer', 'endothelial', 'macrophage', 't_cell'),
		rare_cell_types = c('epithelial', 'mast')
	),
	coadread_lenient = list(
		tcga_cancer_types = c('COAD', 'READ'),
		read_quote = quote(fread('../data_and_figures/lee_crc_2020_smc_reclassified.csv')[cell_type_lenient != 'ambiguous', -c('cell_type_author', 'cell_type')] %>% setnames('cell_type_lenient', 'cell_type')),
		initial_cell_types = c('b_cell', 'caf', 'cancer', 'macrophage', 't_cell'), # Endothelial cells become CAFs under lenient definition
		rare_cell_types = c('epithelial', 'mast')
	),
    hnsc = list(
        tcga_cancer_types = 'HNSC',
        read_quote = quote(fread('../data_and_figures/puram_hnscc_2017_reclassified.csv')[cell_type != 'ambiguous', -'cell_type_author']),
		initial_cell_types = c('caf', 'cancer', 'endothelial', 'macrophage', 'mast', 't_cell'),
		rare_cell_types = c('b_cell', 'dendritic', 'myocyte')
    ),
	lihc = list(
        tcga_cancer_types = 'LIHC',
        read_quote = quote(fread('../data_and_figures/ma_liver_2019_reclassified.csv')[cell_type != 'ambiguous', -c('cell_type_author', 'cell_type_lenient')]),
		initial_cell_types = c('caf', 'cancer', 'endothelial', 'macrophage', 't_cell'),
		rare_cell_types = c('b_cell', 'hpc-like')
    ),
	lihc_lenient = list(
        tcga_cancer_types = 'LIHC',
        read_quote = quote(fread('../data_and_figures/ma_liver_2019_reclassified.csv')[cell_type_lenient != 'ambiguous', -c('cell_type_author', 'cell_type')] %>% setnames('cell_type_lenient', 'cell_type')),
		initial_cell_types = c('caf', 'cancer', 'endothelial', 'macrophage', 't_cell'),
		rare_cell_types = c('b_cell', 'hpc-like')
    ),
	luad = list(
		tcga_cancer_types = 'LUAD',
		read_quote = quote(fread('../data_and_figures/kim_luad_2020_reclassified.csv')[cell_type != 'ambiguous', -'cell_type_author']),
		initial_cell_types = c('b_cell', 'caf', 'cancer', 'endothelial', 'macrophage', 'mast', 't_cell'),
		rare_cell_types = 'dendritic'
	),
    lusc = list(
        tcga_cancer_types = 'LUSC',
        read_quote = quote(fread('../data_and_figures/qian_lung_2020_reclassified.csv')[disease == 'LUSC' & cell_type != 'ambiguous', -c('disease', 'cell_type_author', 'cell_type_lenient')]),
		initial_cell_types = c('b_cell', 'caf', 'cancer', 'endothelial', 'macrophage', 'mast', 't_cell'),
		rare_cell_types = c('dendritic', 'erythroblast')
    ),
	lusc_lenient = list(
        tcga_cancer_types = 'LUSC',
        read_quote = quote(fread('../data_and_figures/qian_lung_2020_reclassified.csv')[disease == 'LUSC' & cell_type_lenient != 'ambiguous', -c('disease', 'cell_type_author', 'cell_type')] %>% setnames('cell_type_lenient', 'cell_type')),
		initial_cell_types = c('b_cell', 'caf', 'cancer', 'endothelial', 'macrophage', 'mast', 't_cell'),
		rare_cell_types = c('dendritic', 'erythroblast')
    ),
	ov = list(
		tcga_cancer_types = 'OV',
		read_quote = quote(fread('../data_and_figures/qian_ovarian_2020_reclassified.csv')[cell_type != 'ambiguous' & !(id %in% c('scrSOL001_TCATTTGTCTGTCAAG', 'scrSOL004_TTGCCGTTCTCCTATA')), -c('cell_type_author', 'cell_type_lenient')]),
		initial_cell_types = c('b_cell', 'caf', 'cancer', 'endothelial', 'macrophage', 't_cell'),
		rare_cell_types = NULL
	),
	ov_lenient = list(
		tcga_cancer_types = 'OV',
		read_quote = quote(fread('../data_and_figures/qian_ovarian_2020_reclassified.csv')[cell_type_lenient != 'ambiguous' & !(id %in% c('scrSOL001_TCATTTGTCTGTCAAG', 'scrSOL004_TTGCCGTTCTCCTATA')), -c('cell_type_author', 'cell_type')] %>% setnames('cell_type_lenient', 'cell_type')),
		initial_cell_types = c('b_cell', 'caf', 'cancer', 'endothelial', 'macrophage', 't_cell'),
		rare_cell_types = NULL
	),
    paad = list(
        tcga_cancer_types = 'PAAD',
        read_quote = quote(fread('../data_and_figures/peng_pdac_2019_reclassified.csv')[cell_type != 'ambiguous' & !(id %in% c('T8_TGGTTCCTCGCATGGC', 'T17_CGTGTAACAGTACACT')), -'cell_type_author']),
		initial_cell_types = c('b_cell', 'caf', 'cancer', 'endothelial', 'macrophage', 't_cell'),
		rare_cell_types = c('acinar', 'ductal_2', 'endocrine')
    )
)





# Simulated bulk data:

simulated_bulk_data <- fread('../data_and_figures/simulated_bulk_data.csv', key = 'id')
simulated_bulk_metadata <- fread('../data_and_figures/simulated_bulk_metadata.csv', key = 'id')

ccle <- fread('../../CCLE_cancer_type_Av.csv', key = 'gene_id')
cell_type_markers <- fread('../../cell_type_markers.csv')

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
        plot_title = 'Head and Neck', genes_filter_fun = function(x) 1:200),
    lihc = list(tcga_cancer_types = 'lihc', ccle_cancer_type = 'liver', seed = 5199,
        plot_title = 'Liver', genes_filter_fun = function(x) 1:200),
	lihc_lenient = list(tcga_cancer_types = 'lihc_lenient', ccle_cancer_type = 'liver', seed = 1682,
        plot_title = 'Liver', genes_filter_fun = function(x) 1:200),
    luad = list(tcga_cancer_types = 'luad', ccle_cancer_type = 'lung', seed = 5395,
        plot_title = 'Lung adeno.', genes_filter_fun = function(x) 1:200),
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

simulated_deconvs <- readRDS('../data_and_figures/simulated_deconvs.rds')

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
        plot_title = 'Lung adeno.'
    ),
	# luad_lenient = list(
        # heatmap_annotations = c('ACTA2', 'COL1A1', 'COL1A2', 'ITGB1', 'LAMC2', 'QSOX1', 'RHOB', 'SNAI1', 'SNAI2', 'THY1', 'TWIST1', 'VIM', 'ZEB1', 'ZEB2'),
        # plot_title = 'Lung adeno.'
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

sample_sizes <- list(
    brca = 1000,
	brca_lenient = 1000,
    coadread = 800,
	coadread_lenient = 800,
    hnsc = 400,
    lihc = 150,
	lihc_lenient = 200,
    luad = 1000,
	lusc = 100,
	lusc_lenient = 100,
    ov = 1500,
	ov_lenient = 2000,
    paad = 1500
)

set.seed(4508) # Is it enough to set a single seed before the whole loop?

sc_deconv_comparison <- sapply(
	names(simulated_deconvs),
	function(ct) {
		cat(paste0(ct, '\n'))
        sc_data <- eval(sc_metadata[[ct]]$read_quote)
		sapply(
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
				
				htmp <- ggplot(
					plot_data,
					aes(
						x = factor(gene, levels = with(simulated_deconvs[[ct]], genes_filtered[ordering])),
						y = factor(id, levels = ordered_cell_ids),
						fill = expression_level
					)
				) +
					geom_raster() +
					scale_x_discrete(expand = c(0, 0)) +
					scale_y_discrete(expand = c(0, 0)) +
					scale_fill_gradientn(
						# colours = c('turquoise4', 'turquoise', 'azure', 'gold2', 'gold4')
						colours = rev(colorRampPalette(brewer.pal(11, "Spectral"))(50)),
						limits = c(0, 12),
						oob = scales::squish,
						breaks = c(0, 3, 6, 9, 12),
						labels = c('0' = '0', '3' = '3', '6' = '6', '9' = '9', '12' = '\u2265 12')
					) +
					theme(
						axis.text = element_blank(),
						axis.title.x = element_blank(),
						axis.ticks = element_blank(),
						axis.ticks.length = unit(0, 'pt'),
						plot.margin = unit(c(5.5, 5.5, 1.25, 1.1), 'pt')
						# plot.margin = switch((x == 'cancer') + 1, unit(c(1.25, 5.5, 5.5, 5.5), 'pt'), unit(c(5.5, 5.5, 1.25, 5.5), 'pt'))
					) +
					labs(
						y = mapvalues(x, c('cancer', 'caf'), c('Cancer', 'CAF'), warn_missing = FALSE),
						fill = 'Expression\nlevel'
					)
				
				ave_exp_bar <- ggplot(
					plot_data[, .(ave_exp = mean(expression_level)), keyby = gene][
						with(simulated_deconvs[[ct]], genes_filtered[ordering]),
						ave_exp := runmean(ave_exp, 30)
					],
					aes(x = factor(gene, levels = with(simulated_deconvs[[ct]], genes_filtered[ordering])), y = 'a', fill = ave_exp)
				) +
					geom_raster() +
					scale_x_discrete(expand = c(0, 0)) +
					scale_y_discrete(expand = c(0, 0)) +
					scale_fill_gradientn(
						# colours = c('turquoise4', 'turquoise', 'azure', 'gold2', 'gold4')
						colours = colorRampPalette(brewer.pal(9, "RdPu"))(50),
						limits = c(0, 4),
						oob = scales::squish,
						breaks = c(0, 2, 4),
						labels = c('0' = '0', '2' = '2', '4' = '\u2265 4')
					) +
					theme(
						axis.text = element_blank(),
						axis.title = element_blank(),
						axis.ticks = element_blank(),
						axis.ticks.length = unit(0, 'pt'),
						plot.margin = switch(
							(x == 'cancer') + 1,
							unit(c(1.25, 5.5, 5.5, 5.5), 'pt'),
							unit(c(1.25, 5.5, 1.25, 5.5), 'pt')
						)
					) +
					labs(
						y = mapvalues(x, c('cancer', 'caf'), c('Cancer', 'CAF'), warn_missing = FALSE),
						fill = 'Average\nexpression\nlevel'
					)
				
				list(heatmap = htmp, ave_exp_bar = ave_exp_bar, plot_data = plot_data, ordered_cell_ids = ordered_cell_ids)
				
			},
			simplify = FALSE,
			USE.NAMES = TRUE
		)
	},
	simplify = FALSE,
	USE.NAMES = TRUE
)

# cairo_pdf('../data_and_figures/sc_sim_deconv_comp.pdf', width = 6, height = 10, onefile = TRUE)

# for(ct in names(sc_deconv_comparison)[!grepl('lenient', names(sc_deconv_comparison))]) {
	
	# plot_grid(
		# plot_grid(
			# plotlist = c(
				# list(blank_plot() + theme(plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), 'pt')) + labs(title = deconv_plot_args_per_ct[[ct]]$plot_title)),
				# lapply(
					# simulated_deconv_plots[[ct]]$plots[c('purity_bar', 'ccle_bar', 'heatmap')],
					# function(x) x + theme(legend.position = 'none', plot.title = element_blank())
				# ),
				# unlist(
					# lapply(
						# sc_deconv_comparison[[ct]],
						# function(x) list(x$heatmap + theme(legend.position = 'none'), x$ave_exp_bar + theme(legend.position = 'none'))
					# ),
					# recursive = FALSE
				# )
			# ),
			# nrow = 8,
			# ncol = 1,
			# align = 'v',
			# rel_heights = c(1, 1, 1, 15, 5, 1, 5, 1)
		# ),
		# plot_grid(
			# plotlist = c(
				# list(
					# blank_plot(),
					# get_legend(
						# simulated_deconv_plots[[1]]$plots$purity_bar +
							# theme(legend.justification = c(0, 1), legend.direction = 'vertical', legend.key.width = NULL, legend.key.height = NULL) +
							# labs(fill = 'Correlation\nwith purity')
					# ),
					# get_legend(
						# simulated_deconv_plots[[1]]$plots$ccle_bar +
							# theme(legend.justification = c(0, 1), legend.direction = 'vertical', legend.key.width = NULL, legend.key.height = NULL) +
							# labs(fill = 'Tumours vs.\ncell lines')
					# ),
					# get_legend(
						# simulated_deconv_plots[[1]]$plots$heatmap + theme(legend.justification = c(0, 1))
					# )
				# ),
				# list(
					# get_legend(sc_deconv_comparison[[1]][[1]]$heatmap + theme(legend.justification = c(0, 1))),
					# get_legend(sc_deconv_comparison[[1]][[1]]$ave_exp_bar + theme(legend.justification = c(0, 1)))
				# )
			# ),
			# nrow = 6,
			# ncol = 1,
			# rel_heights = c(1, 5, 5, 7, 6, 6)
		# ),
		# nrow = 1,
		# ncol = 2,
		# rel_widths = c(5, 1)
	# ) %>% print
	
# }

# dev.off()

# cairo_pdf('../data_and_figures/sc_sim_deconv_comp_lenient.pdf', width = 6, height = 10, onefile = TRUE)

# for(ct in names(sc_deconv_comparison)[grepl('lenient', names(sc_deconv_comparison))]) {
	
	# plot_grid(
		# plot_grid(
			# plotlist = c(
				# list(blank_plot() + theme(plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), 'pt')) + labs(title = deconv_plot_args_per_ct[[ct]]$plot_title)),
				# lapply(
					# simulated_deconv_plots[[ct]]$plots[c('purity_bar', 'ccle_bar', 'heatmap')],
					# function(x) x + theme(legend.position = 'none', plot.title = element_blank())
				# ),
				# unlist(
					# lapply(
						# sc_deconv_comparison[[ct]],
						# function(x) list(x$heatmap + theme(legend.position = 'none'), x$ave_exp_bar + theme(legend.position = 'none'))
					# ),
					# recursive = FALSE
				# )
			# ),
			# nrow = 8,
			# ncol = 1,
			# align = 'v',
			# rel_heights = c(1, 1, 1, 15, 5, 1, 5, 1)
		# ),
		# plot_grid(
			# plotlist = c(
				# list(
					# blank_plot(),
					# get_legend(
						# simulated_deconv_plots[[1]]$plots$purity_bar +
							# theme(legend.justification = c(0, 1), legend.direction = 'vertical', legend.key.width = NULL, legend.key.height = NULL) +
							# labs(fill = 'Correlation\nwith purity')
					# ),
					# get_legend(
						# simulated_deconv_plots[[1]]$plots$ccle_bar +
							# theme(legend.justification = c(0, 1), legend.direction = 'vertical', legend.key.width = NULL, legend.key.height = NULL) +
							# labs(fill = 'Tumours vs.\ncell lines')
					# ),
					# get_legend(
						# simulated_deconv_plots[[1]]$plots$heatmap + theme(legend.justification = c(0, 1))
					# )
				# ),
				# list(
					# get_legend(sc_deconv_comparison[[1]][[1]]$heatmap + theme(legend.justification = c(0, 1))),
					# get_legend(sc_deconv_comparison[[1]][[1]]$ave_exp_bar + theme(legend.justification = c(0, 1)))
				# )
			# ),
			# nrow = 6,
			# ncol = 1,
			# rel_heights = c(1, 5, 5, 7, 6, 6)
		# ),
		# nrow = 1,
		# ncol = 2,
		# rel_widths = c(5, 1)
	# ) %>% print
	
# }

# dev.off()

# Make centred versions:

sc_deconv_comparison_centred <- sapply(
	names(sc_deconv_comparison),
	function(ct) {
		
		cat(paste0(ct, '\n'))
		
		plot_data <- rbind(
			sc_deconv_comparison[[ct]]$cancer$plot_data[, cell_type := 'cancer'],
			sc_deconv_comparison[[ct]]$caf$plot_data[, cell_type := 'caf']
		)
		
		# Calculate averages in cancer cells and CAFs separately, so we can check them afterwards:
		gene_averages_cancer_caf <- plot_data[
			,
			.(ave_exp = mean(expression_level)),
			by = .(gene, cell_type)
		]
		
		# To centre genes w.r.t. the average of the averages of cancer and CAF:
		gene_averages <- plot_data[
			,
			.(ave_exp = mean(c(mean(expression_level[cell_type == 'cancer']), mean(expression_level[cell_type == 'caf'])))),
			by = .(symbol = gene)
		]
		
		plot_data[, expression_level := expression_level - gene_averages[symbol == unique(gene), ave_exp], by = gene]
		
		# To centre the cells as well:
		plot_data[, expression_level := expression_level - mean(expression_level), by = id]
		
		htmps <- sapply(
			c('cancer', 'caf'),
			function(x) {
				ggplot(
					plot_data[cell_type == x],
					aes(
						x = factor(gene, levels = with(simulated_deconvs[[ct]], genes_filtered[ordering])),
						y = factor(id, levels = sc_deconv_comparison[[ct]][[x]]$ordered_cell_ids),
						fill = expression_level
					)
				) +
					geom_raster() +
					scale_x_discrete(expand = c(0, 0)) +
					scale_y_discrete(expand = c(0, 0)) +
					scale_fill_gradientn(
						# colours = c('turquoise4', 'turquoise', 'azure', 'gold2', 'gold4')
						colours = rev(colorRampPalette(brewer.pal(11, "RdBu"))(50)),
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
						plot.margin = unit(c(5.5, 5.5, 1.25, 1.1), 'pt')
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
		
		list(heatmaps = htmps, plot_data = plot_data, gene_averages = gene_averages, gene_averages_cancer_caf = gene_averages_cancer_caf)
		
	},
	simplify = FALSE,
	USE.NAMES = TRUE
)

# cairo_pdf('../data_and_figures/sc_sim_deconv_comp_centred.pdf', width = 6, height = 9, onefile = TRUE)

# for(ct in names(sc_deconv_comparison_centred)[!grepl('lenient', names(sc_deconv_comparison_centred))]) {
	
	# plot_grid(
		# plot_grid(
			# plotlist = c(
				# list(blank_plot() + theme(plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), 'pt')) + labs(title = deconv_plot_args_per_ct[[ct]]$plot_title)),
				# lapply(
					# simulated_deconv_plots[[ct]]$plots[c('purity_bar', 'ccle_bar', 'heatmap')],
					# function(x) x + theme(legend.position = 'none', plot.title = element_blank())
				# ),
				# lapply(sc_deconv_comparison_centred[[ct]]$heatmaps, function(x) x + theme(legend.position = 'none'))
			# ),
			# nrow = 6,
			# ncol = 1,
			# align = 'v',
			# rel_heights = c(1, 1, 1, 15, 5, 5)
		# ),
		# plot_grid(
			# plotlist = c(
				# list(
					# blank_plot(),
					# get_legend(
						# simulated_deconv_plots[[1]]$plots$purity_bar +
							# theme(legend.justification = c(0, 1), legend.direction = 'vertical', legend.key.width = NULL, legend.key.height = NULL) +
							# labs(fill = 'Correlation\nwith purity')
					# ),
					# get_legend(
						# simulated_deconv_plots[[1]]$plots$ccle_bar +
							# theme(legend.justification = c(0, 1), legend.direction = 'vertical', legend.key.width = NULL, legend.key.height = NULL) +
							# labs(fill = 'Tumours vs.\ncell lines')
					# ),
					# get_legend(
						# simulated_deconv_plots[[1]]$plots$heatmap + theme(legend.justification = c(0, 1))
					# ),
					# get_legend(sc_deconv_comparison_centred[[1]]$heatmaps[[1]] + theme(legend.justification = c(0, 1)))
				# )
			# ),
			# nrow = 5,
			# ncol = 1,
			# rel_heights = c(1, 5, 5, 7, 10)
		# ),
		# nrow = 1,
		# ncol = 2,
		# rel_widths = c(5, 1)
	# ) %>% print
	
# }

# dev.off()

# cairo_pdf('../data_and_figures/sc_sim_deconv_comp_lenient_centred.pdf', width = 6, height = 10, onefile = TRUE)

# for(ct in names(sc_deconv_comparison_centred)[grepl('lenient', names(sc_deconv_comparison_centred))]) {
	
	# plot_grid(
		# plot_grid(
			# plotlist = c(
				# list(blank_plot() + theme(plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), 'pt')) + labs(title = deconv_plot_args_per_ct[[ct]]$plot_title)),
				# lapply(
					# simulated_deconv_plots[[ct]]$plots[c('purity_bar', 'ccle_bar', 'heatmap')],
					# function(x) x + theme(legend.position = 'none', plot.title = element_blank())
				# ),
				# lapply(sc_deconv_comparison_centred[[ct]]$heatmaps, function(x) x + theme(legend.position = 'none'))
			# ),
			# nrow = 6,
			# ncol = 1,
			# align = 'v',
			# rel_heights = c(1, 1, 1, 15, 5, 5)
		# ),
		# plot_grid(
			# plotlist = c(
				# list(
					# blank_plot(),
					# get_legend(
						# simulated_deconv_plots[[1]]$plots$purity_bar +
							# theme(legend.justification = c(0, 1), legend.direction = 'vertical', legend.key.width = NULL, legend.key.height = NULL) +
							# labs(fill = 'Correlation\nwith purity')
					# ),
					# get_legend(
						# simulated_deconv_plots[[1]]$plots$ccle_bar +
							# theme(legend.justification = c(0, 1), legend.direction = 'vertical', legend.key.width = NULL, legend.key.height = NULL) +
							# labs(fill = 'Tumours vs.\ncell lines')
					# ),
					# get_legend(
						# simulated_deconv_plots[[1]]$plots$heatmap + theme(legend.justification = c(0, 1))
					# ),
					# get_legend(sc_deconv_comparison_centred[[1]]$heatmaps[[1]] + theme(legend.justification = c(0, 1)))
				# )
			# ),
			# nrow = 5,
			# ncol = 1,
			# rel_heights = c(1, 5, 5, 7, 10)
		# ),
		# nrow = 1,
		# ncol = 2,
		# rel_widths = c(5, 1)
	# ) %>% print
	
# }

# dev.off()

# Make versions where we filter out genes with low expression in the single cell data:

sc_deconv_comparison_filtered <- sapply(
	names(sc_deconv_comparison_centred),
	function(ct) {
		
		cat(paste0(ct, '\n'))
		
		simulated_deconv_ct <- do.call(
			deconvolve_emt_caf_data,
			args = c(
				list(
					expression_data = simulated_bulk_data,
					meta_data = simulated_bulk_metadata,
					genes = sc_deconv_comparison_centred[[ct]]$gene_averages_cancer_caf[
						,
						.(pass = ave_exp[cell_type == 'cancer'] > 0.5 | ave_exp[cell_type == 'caf'] > 0.5),
						by = gene
					][pass == TRUE, as.character(gene)],
					# genes = sc_deconv_comparison_centred[[ct]]$gene_averages[ave_exp >= 0.5, as.character(symbol)], # symbol is a factor for some reason
					# genes = gene_averages[
						# ,
						# .(pass = ave_exp[cancer_caf == 'cancer'] > 0.5 | ave_exp[cancer_caf == 'caf'] > 0.5),
						# by = symbol
					# ][pass == TRUE, as.character(symbol)],
					cell_type_markers = cell_type_markers,
					ccle_data = ccle,
					genes_from_tcga_fun = NULL,
					genes_filter_fun = NULL
				),
				deconv_args_per_ct[[ct]][!(names(deconv_args_per_ct[[ct]]) %in% c('plot_title', 'genes_filter_fun'))]
			)
		)
		
		simulated_deconv_ct <- deconv_reorder(simulated_deconv_ct)
		
		simulated_deconv_plots_ct <- do.call(
			deconvolve_emt_caf_plots,
			args = c(
				list(
					data = simulated_deconv_ct,
					heatmap_legend_title = 'Correlation',
					heatmap_colours = rev(colorRampPalette(brewer.pal(11, "RdBu"))(50)),
					heatmap_colour_limits = c(-1, 1),
					heatmap_legend_breaks = c(-1, -0.5, 0, 0.5, 1),
					heatmap_annotations_nudge = 0.3,
					purity_colours = rev(colorRampPalette(brewer.pal(11, "PuOr"))(50)),
					purity_colour_limits = c(-1, 1),
					purity_legend_breaks = c(-1, 0, 1),
					purity_legend_title = 'Correlation with purity\n',
					purity_legend_direction = 'horizontal',
					purity_axis_title = NULL,
					ccle_colours = rev(colorRampPalette(brewer.pal(11, 'PiYG')[c(1:3, 6, 9:11)])(50)),
					ccle_legend_breaks = c(-1, 0, 1),
					ccle_legend_title = 'Tumours vs. cell lines\n',
					ccle_legend_direction = 'horizontal',
					ccle_axis_title = NULL,
					bar_legend_justification = 'left'
				),
				deconv_plot_args_per_ct[[ct]]
			)
		)
		
		plot_data <- sc_deconv_comparison_centred[[ct]]$plot_data[
			gene %in% sc_deconv_comparison_centred[[ct]]$gene_averages_cancer_caf[
				,
				.(pass = ave_exp[cell_type == 'cancer'] > 0.5 | ave_exp[cell_type == 'caf'] > 0.5),
				by = gene
			][pass == TRUE, as.character(gene)]
			# gene %in% sc_deconv_comparison_centred[[ct]]$gene_averages[ave_exp >= 0.5, symbol]
		]
		
		sc_heatmaps <- sapply(
			c('cancer', 'caf'),
			function(x) {
				ggplot(
					plot_data[cell_type == x],
					aes(
						x = factor(gene, levels = with(simulated_deconv_ct, genes_filtered[ordering])),
						y = factor(id, levels = sc_deconv_comparison[[ct]][[x]]$ordered_cell_ids),
						fill = expression_level
					)
				) +
					geom_raster() +
					scale_x_discrete(expand = c(0, 0)) +
					scale_y_discrete(expand = c(0, 0)) +
					scale_fill_gradientn(
						colours = rev(colorRampPalette(brewer.pal(11, "RdBu"))(50)),
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
					) +
					labs(y = mapvalues(x, c('cancer', 'caf'), c('Cancer', 'CAF'), warn_missing = FALSE), fill = 'Relative\nexpression\nlevel')
			},
			simplify = FALSE,
			USE.NAMES = TRUE
		)
		
		list(deconv_data = simulated_deconv_ct, deconv_figures = simulated_deconv_plots_ct$plots, sc_heatmaps = sc_heatmaps, plot_data = plot_data)
		
	},
	simplify = FALSE,
	USE.NAMES = TRUE
)

cairo_pdf('../data_and_figures/sc_sim_deconv_comp_filtered.pdf', width = 6, height = 9, onefile = TRUE)

for(ct in names(sc_deconv_comparison_filtered)[!grepl('lenient', names(sc_deconv_comparison_filtered))]) {
	
	plot_grid(
		plot_grid(
			plotlist = c(
				list(blank_plot() + theme(plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), 'pt')) + labs(title = deconv_plot_args_per_ct[[ct]]$plot_title)),
				lapply(
					c(
						sc_deconv_comparison_filtered[[ct]]$deconv_figures[c('purity_bar', 'ccle_bar', 'heatmap')],
						sc_deconv_comparison_filtered[[ct]]$sc_heatmaps
					),
					function(x) x + theme(legend.position = 'none', plot.title = element_blank())
				)
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
							theme(legend.justification = c(0, 1), legend.direction = 'vertical', legend.key.width = NULL, legend.key.height = NULL) +
							labs(fill = 'Correlation\nwith purity')
					),
					get_legend(
						simulated_deconv_plots[[1]]$plots$ccle_bar +
							theme(legend.justification = c(0, 1), legend.direction = 'vertical', legend.key.width = NULL, legend.key.height = NULL) +
							labs(fill = 'Tumours vs.\ncell lines')
					),
					get_legend(
						simulated_deconv_plots[[1]]$plots$heatmap + theme(legend.justification = c(0, 1))
					),
					get_legend(sc_deconv_comparison_filtered[[1]]$sc_heatmaps[[1]] + theme(legend.justification = c(0, 1)))
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

for(ct in names(sc_deconv_comparison_filtered)[grepl('lenient', names(sc_deconv_comparison_filtered))]) {
	
	plot_grid(
		plot_grid(
			plotlist = c(
				list(blank_plot() + theme(plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), 'pt')) + labs(title = deconv_plot_args_per_ct[[ct]]$plot_title)),
				lapply(
					c(
						sc_deconv_comparison_filtered[[ct]]$deconv_figures[c('purity_bar', 'ccle_bar', 'heatmap')],
						sc_deconv_comparison_filtered[[ct]]$sc_heatmaps
					),
					function(x) x + theme(legend.position = 'none', plot.title = element_blank())
				)
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
							theme(legend.justification = c(0, 1), legend.direction = 'vertical', legend.key.width = NULL, legend.key.height = NULL) +
							labs(fill = 'Correlation\nwith purity')
					),
					get_legend(
						simulated_deconv_plots[[1]]$plots$ccle_bar +
							theme(legend.justification = c(0, 1), legend.direction = 'vertical', legend.key.width = NULL, legend.key.height = NULL) +
							labs(fill = 'Tumours vs.\ncell lines')
					),
					get_legend(
						simulated_deconv_plots[[1]]$plots$heatmap + theme(legend.justification = c(0, 1))
					),
					get_legend(sc_deconv_comparison_filtered[[1]]$sc_heatmaps[[1]] + theme(legend.justification = c(0, 1)))
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

# Make versions where we also filter out genes that don't fit the trend:

# We will want different thresholds for "anomalous" genes in each cancer type:
anomalous_threshold <- c(
	brca = 0.5,
	brca_lenient = 0.5,
	coadread = 0.6,
	coadread_lenient = 0.75,
	hnsc = 0.75,
	lihc = 0.75,
	lihc_lenient = 0.75,
	luad = 0.75,
	lusc = 0.5,
	lusc_lenient = 0.5,
	ov = 0.75,
	ov_lenient = 0.75,
	paad = 0.6
)

sc_deconv_comparison_filtered_twice <- sapply(
	names(sc_deconv_comparison_centred),
	function(ct) {
		
		cat(paste0(ct, '\n'))
		
		simulated_deconv_ct <- do.call(
			deconvolve_emt_caf_data,
			args = c(
				list(
					expression_data = simulated_bulk_data,
					meta_data = simulated_bulk_metadata,
					genes = sc_deconv_comparison_centred[[ct]]$gene_averages_cancer_caf[
						,
						.(pass = ave_exp[cell_type == 'cancer'] > 0.5 | ave_exp[cell_type == 'caf'] > 0.5),
						by = gene
					][pass == TRUE, as.character(gene)],
					# genes = sc_deconv_comparison_centred[[ct]]$gene_averages[ave_exp >= 0.5, as.character(symbol)], # symbol is a factor for some reason
					# genes = gene_averages[
						# ,
						# .(pass = ave_exp[cancer_caf == 'cancer'] > 0.5 | ave_exp[cancer_caf == 'caf'] > 0.5),
						# by = symbol
					# ][pass == TRUE, as.character(symbol)],
					cell_type_markers = cell_type_markers,
					ccle_data = ccle,
					genes_from_tcga_fun = NULL,
					genes_filter_fun = NULL
				),
				deconv_args_per_ct[[ct]][!(names(deconv_args_per_ct[[ct]]) %in% c('plot_title', 'genes_filter_fun'))]
			)
		)
		
		# Calculate new scores as in deconv_reorder() function.  Genes with positive score should be on the left side and those with
		# negative score should be on the right side.
		head_genes <- with(simulated_deconv_ct, head(genes_filtered[ordering], 20))
		tail_genes <- with(simulated_deconv_ct, tail(genes_filtered[ordering], 20))
		new_scores <- sapply(
			colnames(simulated_deconv_ct$cor_mat),
			function(g) {
				meanvec <- c(mean(simulated_deconv_ct$cor_mat[, g][head_genes]), mean(simulated_deconv_ct$cor_mat[, g][tail_genes]))
				ifelse(
					sign(meanvec[1]) == sign(meanvec[2]),
					ifelse(g %in% c('SNAI1', 'SNAI2', 'TWIST1', 'VIM','ZEB1', 'ZEB2'), return(0), return(NA)),
					return(max(meanvec)*c(1, -1)[which.max(meanvec)])
				)
			}
		)
		new_scores <- new_scores[!is.na(new_scores)]
		
		# For each gene, calculate the difference between its average across cancer cells and its average across CAFs.  For genes with
		# positive score, this difference should be positive, and for those with negative score it should be negative.  There will be
		# some exceptions - we need to decide if an exception is "extreme" enough.
		exp_diff_table <- sc_deconv_comparison_centred[[ct]]$plot_data[
			,
			.(exp_diff = mean(expression_level[cell_type == 'cancer']) - mean(expression_level[cell_type == 'caf'])),
			by = gene
		]
		
		# For each gene, consider the distribution of these differences for the genes having score with opposite sign.  If the gene is
		# beyond a certain quantile of this distribution ("beyond" w.r.t. its "own" distribution) then it is "anomalous".
		new_scores <- sapply(
			names(new_scores),
			function(g) {
				ifelse(
					new_scores[g] > 0,
					ifelse(
						exp_diff_table[
							,
							exp_diff[gene == g] < quantile(exp_diff[gene %in% names(new_scores)[new_scores < 0]], anomalous_threshold[ct])
						] & !(g %in% c('SNAI1', 'SNAI2', 'TWIST1', 'VIM', 'ZEB1', 'ZEB2')),
						return(NA),
						return(setNames(new_scores[g], NULL))
					),
					ifelse(
						new_scores[g] < 0,
						ifelse(
							exp_diff_table[
								,
								exp_diff[gene == g] > quantile(exp_diff[gene %in% names(new_scores)[new_scores > 0]], 1 - anomalous_threshold[ct])
							] & !(g %in% c('SNAI1', 'SNAI2', 'TWIST1', 'VIM', 'ZEB1', 'ZEB2')),
							return(NA),
							return(setNames(new_scores[g], NULL))
						),
						return(0)
					)
				)
				# exp_diff_table[, exp_diff[names(new_scores)[sign(new_scores) != sign(new_scores[g])]]
			}
		)
		
		anomalous_genes <- names(new_scores)[is.na(new_scores)]
		new_scores <- new_scores[!is.na(new_scores)]
		
		simulated_deconv_ct$anomalous_genes <- anomalous_genes
		simulated_deconv_ct$new_scores <- new_scores
		simulated_deconv_ct$cor_mat <- simulated_deconv_ct$cor_mat[names(new_scores), names(new_scores)]
		simulated_deconv_ct$genes_filtered <- names(new_scores)
		simulated_deconv_ct$ordering <- order(-new_scores)
		simulated_deconv_ct$cor_with_purity <- sapply(
			simulated_deconv_ct$cor_with_purity,
			function(x) x[names(new_scores)],
			simplify = FALSE,
			USE.NAMES = TRUE
		)
		simulated_deconv_ct$ccle_comp_diff <- simulated_deconv_ct$ccle_comp_diff[names(new_scores)]
		
		simulated_deconv_plots_ct <- do.call(
			deconvolve_emt_caf_plots,
			args = c(
				list(
					data = simulated_deconv_ct,
					heatmap_legend_title = 'Correlation',
					heatmap_colours = rev(colorRampPalette(brewer.pal(11, "RdBu"))(50)),
					heatmap_colour_limits = c(-1, 1),
					heatmap_legend_breaks = c(-1, -0.5, 0, 0.5, 1),
					heatmap_annotations_nudge = 0.3,
					purity_colours = rev(colorRampPalette(brewer.pal(11, "PuOr"))(50)),
					purity_colour_limits = c(-1, 1),
					purity_legend_breaks = c(-1, 0, 1),
					purity_legend_title = 'Correlation with purity\n',
					purity_legend_direction = 'horizontal',
					purity_axis_title = NULL,
					ccle_colours = rev(colorRampPalette(brewer.pal(11, 'PiYG')[c(1:3, 6, 9:11)])(50)),
					ccle_legend_breaks = c(-1, 0, 1),
					ccle_legend_title = 'Tumours vs. cell lines\n',
					ccle_legend_direction = 'horizontal',
					ccle_axis_title = NULL,
					bar_legend_justification = 'left'
				),
				deconv_plot_args_per_ct[[ct]]
			)
		)
		
		plot_data <- sc_deconv_comparison_centred[[ct]]$plot_data[
			gene %in% sc_deconv_comparison_centred[[ct]]$gene_averages_cancer_caf[
				,
				.(pass = ave_exp[cell_type == 'cancer'] > 0.5 | ave_exp[cell_type == 'caf'] > 0.5),
				by = gene
			][pass == TRUE, as.character(gene)]
			# gene %in% sc_deconv_comparison_centred[[ct]]$gene_averages[ave_exp >= 0.5, symbol]
		]
		
		sc_heatmaps <- sapply(
			c('cancer', 'caf'),
			function(x) {
				ggplot(
					plot_data[cell_type == x],
					aes(
						x = factor(gene, levels = with(simulated_deconv_ct, genes_filtered[ordering])),
						y = factor(id, levels = sc_deconv_comparison[[ct]][[x]]$ordered_cell_ids),
						fill = expression_level
					)
				) +
					geom_raster() +
					scale_x_discrete(expand = c(0, 0)) +
					scale_y_discrete(expand = c(0, 0)) +
					scale_fill_gradientn(
						colours = rev(colorRampPalette(brewer.pal(11, "RdBu"))(50)),
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
					) +
					labs(y = mapvalues(x, c('cancer', 'caf'), c('Cancer', 'CAF'), warn_missing = FALSE), fill = 'Relative\nexpression\nlevel')
			},
			simplify = FALSE,
			USE.NAMES = TRUE
		)
		
		list(deconv_data = simulated_deconv_ct, deconv_figures = simulated_deconv_plots_ct$plots, sc_heatmaps = sc_heatmaps, plot_data = plot_data)
		
	},
	simplify = FALSE,
	USE.NAMES = TRUE
)

cairo_pdf('../data_and_figures/sc_sim_deconv_comp_filtered_twice.pdf', width = 6, height = 9, onefile = TRUE)

for(ct in names(sc_deconv_comparison_filtered_twice)[!grepl('lenient', names(sc_deconv_comparison_filtered_twice))]) {
	
	plot_grid(
		plot_grid(
			plotlist = c(
				list(blank_plot() + theme(plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), 'pt')) + labs(title = deconv_plot_args_per_ct[[ct]]$plot_title)),
				lapply(
					c(
						sc_deconv_comparison_filtered_twice[[ct]]$deconv_figures[c('purity_bar', 'ccle_bar', 'heatmap')],
						sc_deconv_comparison_filtered_twice[[ct]]$sc_heatmaps
					),
					function(x) x + theme(legend.position = 'none', plot.title = element_blank())
				)
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
							theme(legend.justification = c(0, 1), legend.direction = 'vertical', legend.key.width = NULL, legend.key.height = NULL) +
							labs(fill = 'Correlation\nwith purity')
					),
					get_legend(
						simulated_deconv_plots[[1]]$plots$ccle_bar +
							theme(legend.justification = c(0, 1), legend.direction = 'vertical', legend.key.width = NULL, legend.key.height = NULL) +
							labs(fill = 'Tumours vs.\ncell lines')
					),
					get_legend(
						simulated_deconv_plots[[1]]$plots$heatmap + theme(legend.justification = c(0, 1))
					),
					get_legend(sc_deconv_comparison_filtered_twice[[1]]$sc_heatmaps[[1]] + theme(legend.justification = c(0, 1)))
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

cairo_pdf('../data_and_figures/sc_sim_deconv_comp_lenient_filtered_twice.pdf', width = 6, height = 10, onefile = TRUE)

for(ct in names(sc_deconv_comparison_filtered_twice)[grepl('lenient', names(sc_deconv_comparison_filtered_twice))]) {
	
	plot_grid(
		plot_grid(
			plotlist = c(
				list(blank_plot() + theme(plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), 'pt')) + labs(title = deconv_plot_args_per_ct[[ct]]$plot_title)),
				lapply(
					c(
						sc_deconv_comparison_filtered_twice[[ct]]$deconv_figures[c('purity_bar', 'ccle_bar', 'heatmap')],
						sc_deconv_comparison_filtered_twice[[ct]]$sc_heatmaps
					),
					function(x) x + theme(legend.position = 'none', plot.title = element_blank())
				)
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
							theme(legend.justification = c(0, 1), legend.direction = 'vertical', legend.key.width = NULL, legend.key.height = NULL) +
							labs(fill = 'Correlation\nwith purity')
					),
					get_legend(
						simulated_deconv_plots[[1]]$plots$ccle_bar +
							theme(legend.justification = c(0, 1), legend.direction = 'vertical', legend.key.width = NULL, legend.key.height = NULL) +
							labs(fill = 'Tumours vs.\ncell lines')
					),
					get_legend(
						simulated_deconv_plots[[1]]$plots$heatmap + theme(legend.justification = c(0, 1))
					),
					get_legend(sc_deconv_comparison_filtered_twice[[1]]$sc_heatmaps[[1]] + theme(legend.justification = c(0, 1)))
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





# TCGA bulk data:

deconv_data <- readRDS('../data_and_figures/deconv_data.rds')
deconv_plots <- readRDS('../data_and_figures/deconv_plots.rds')
cell_type_markers <- fread('../../cell_type_markers.csv')
subtypes_data <- fread('../../TCGA_data/tcga_subtypes_data.csv', key = 'id')
ccle <- fread('../../CCLE_cancer_type_Av.csv', key = 'gene_id')

ct_to_keep <- c('brca_luminal_a', 'brca_luminal_b', 'brca_basal_like', 'brca_her2_enriched', 'coad', 'hnsc_mesenchymal_basal',
				'hnsc_classical', 'hnsc_atypical', 'lihc', 'luad_proximal_inflammatory', 'luad_proximal_proliferative', 'lusc_basal',
				'lusc_classical', 'lusc_primitive', 'lusc_secretory', 'ov_differentiated', 'ov_immunoreactive', 'ov_proliferative',
				'paad', 'read')

deconv_data <- deconv_data[ct_to_keep]
deconv_plots <- deconv_plots[ct_to_keep]

tcga_sc_map <- c(
	brca_luminal_a = 'brca',
	brca_luminal_b = 'brca',
	brca_basal_like = 'brca',
	brca_her2_enriched = 'brca',
	coad = 'coadread',
	hnsc_mesenchymal_basal = 'hnsc',
	hnsc_classical = 'hnsc',
	hnsc_atypical = 'hnsc',
	lihc = 'lihc',
	luad_proximal_inflammatory = 'luad',
	luad_proximal_proliferative = 'luad',
	lusc_basal = 'lusc',
	lusc_classical = 'lusc',
	lusc_primitive = 'lusc',
	lusc_secretory = 'lusc',
	ov_differentiated = 'ov',
	ov_immunoreactive = 'ov',
	ov_proliferative = 'ov',
	paad = 'paad',
	read = 'coadread'
)

deconv_args_per_ct <- list(
    
    brca_luminal_a = list(
        tcga_cancer_types = 'BRCA',
        subtypes = 'Luminal A',
        ref_for_subtypes = 'doi:10.1038/nature11412',
        ccle_cancer_type = 'breast',
        extra_data_source = 'breast_qian',
        seed = 8216,
        plot_title = 'BRCA - Luminal A',
        genes_filter_fun = function(x) 1:175,
        gene_weights_fun = function(x) quantile(x, 0.9)
    ),
    
    brca_luminal_b = list(
        tcga_cancer_types = 'BRCA',
        subtypes = 'Luminal B',
        ref_for_subtypes = 'doi:10.1038/nature11412',
        ccle_cancer_type = 'breast',
        extra_data_source = 'breast_qian',
        seed = 9056,
        plot_title = 'BRCA - Luminal B',
        genes_filter_fun = function(x) 1:230,
        gene_weights_fun = function(x) quantile(x, 0.75)
    ),
    
    brca_basal_like = list(
        tcga_cancer_types = 'BRCA',
        subtypes = 'Basal-like',
        ref_for_subtypes = 'doi:10.1038/nature11412',
        ccle_cancer_type = 'breast',
        extra_data_source = 'breast_qian',
        seed = 1982,
        plot_title = 'BRCA - Basal-like',
        genes_filter_fun = function(x) 1:230,
        gene_weights_fun = function(x) quantile(x, 0.75)
    ),
    
    brca_her2_enriched = list(
        tcga_cancer_types = 'BRCA',
        subtypes = 'HER2-enriched',
        ref_for_subtypes = 'doi:10.1038/nature11412',
        ccle_cancer_type = 'breast',
        extra_data_source = 'breast_qian',
        seed = 6994,
        plot_title = 'BRCA - HER2-enriched',
        genes_filter_fun = function(x) 1:180,
        gene_weights_fun = function(x) quantile(x, 0.75)
    ),
    
    coad = list(
        tcga_cancer_types = 'COAD',
        ccle_cancer_type = 'large intestine',
        extra_data_source = 'crc_lee_smc',
        seed = 3260,
        plot_title = 'COAD',
        genes_filter_fun = function(x) 1:200,
        cell_type_weights = list(
            B_plasma = 0,
            myocyte = 0,
            macrophage = 0,
            endothelial = 1,
            DC = 0,
            mast = 0,
            T = 0,
            B = 0
        )
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
        plot_title = 'HNSC - Classical',
        gene_weights_fun = function(x) quantile(x, 0.75)
    ),
    
    hnsc_atypical = list(
        tcga_cancer_types = 'HNSC',
        subtypes = 'Atypical',
        ref_for_subtypes = 'doi:10.1038/nature14129',
        ccle_cancer_type = 'HNSCC',
        extra_data_source = 'hnscc_puram',
        seed = 8172,
        plot_title = 'HNSC - Atypical',
        gene_weights_fun = function(x) quantile(x, 0.9)
    ),
    
    lihc = list(
        tcga_cancer_types = 'LIHC',
        ccle_cancer_type = 'liver',
        extra_data_source = 'liver_ma',
        seed = 9349,
        plot_title = 'LIHC',
        genes_filter_fun = function(x) 1:260,
        cell_type_weights = FALSE
    ),
    
    luad_proximal_inflammatory = list(
        tcga_cancer_types = 'LUAD',
        subtypes = 'Proximal-inflammatory',
        ref_for_subtypes = 'doi:10.1038/nature13385',
        ccle_cancer_type = 'lung',
		extra_data_source = 'luad_kim',
        # extra_data_source = 'luad_qian',
        seed = 9955,
        plot_title = 'LUAD - Proximal-inflammatory'
    ),
    
    luad_proximal_proliferative = list(
        tcga_cancer_types = 'LUAD',
        subtypes = 'Proximal-proliferative',
        ref_for_subtypes = 'doi:10.1038/nature13385',
        ccle_cancer_type = 'lung',
		extra_data_source = 'luad_kim',
        # extra_data_source = 'luad_qian',
        seed = 743,
        plot_title = 'LUAD - Proximal-proliferative',
        genes_filter_fun = function(x) 1:225,
        gene_weights_fun = function(x) quantile(x, 0.75)
    ),
    
    lusc_basal = list(
        tcga_cancer_types = 'LUSC',
        subtypes = 'basal',
        ref_for_subtypes = 'doi:10.1038/nature11404',
        ccle_cancer_type = 'lung',
        extra_data_source = 'lusc_qian',
        seed = 2287,
        plot_title = 'LUSC - Basal',
        gene_weights_fun = function(x) quantile(x, 0.9)
    ),
    
    lusc_classical = list(
        tcga_cancer_types = 'LUSC',
        subtypes = 'classical',
        ref_for_subtypes = 'doi:10.1038/nature11404',
        ccle_cancer_type = 'lung',
        extra_data_source = 'lusc_qian',
        seed = 7106,
        plot_title = 'LUSC - Classical',
        gene_weights_fun = function(x) quantile(x, 0.9)
    ),
    
    lusc_primitive = list(
        tcga_cancer_types = 'LUSC',
        subtypes = 'primitive',
        ref_for_subtypes = 'doi:10.1038/nature11404',
        ccle_cancer_type = 'lung',
        extra_data_source = 'lusc_qian',
        seed = 303,
        plot_title = 'LUSC - Primitive',
        gene_weights_fun = function(x) quantile(x, 0.9)
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
        plot_title = 'OV - Differentiated',
        genes_filter_fun = function(x) 1:270,
        gene_weights_fun = function(x) quantile(x, 0.9)
    ),
    
    ov_immunoreactive = list(
        tcga_cancer_types = 'OV',
        subtypes = 'Immunoreactive',
        ref_for_subtypes = 'doi:10.1172/JCI65833',
        ccle_cancer_type = 'ovary',
		extra_data_source = 'ovarian_qian',
        seed = 991,
        plot_title = 'OV - Immunoreactive',
        genes_filter_fun = function(x) 1:225,
        gene_weights_fun = function(x) quantile(x, 0.9)
    ),
    
    ov_proliferative = list(
        tcga_cancer_types = 'OV',
        subtypes = 'Proliferative',
        ref_for_subtypes = 'doi:10.1172/JCI65833',
        ccle_cancer_type = 'ovary',
		extra_data_source = 'ovarian_qian',
        seed = 8583,
        plot_title = 'OV - Proliferative',
        gene_weights_fun = function(x) quantile(x, 0.9)
    ),
    
    paad = list(
        tcga_cancer_types = 'PAAD',
        ccle_cancer_type = 'pancreas',
        extra_data_source = 'pdac_peng',
        seed = 88,
        plot_title = 'PAAD',
        genes_filter_fun = function(x) 1:150, # Also works OK with 175, but 150 is better
        cell_type_weights = list(
            B_plasma = 0,
            myocyte = 0,
            macrophage = 0,
            endothelial = 1,
            DC = 0,
            mast = 0,
            T = 0,
            B = 0
        )
    ),
    
    read = list(
        tcga_cancer_types = 'READ',
        ccle_cancer_type = 'large intestine',
        seed = 3309,
        plot_title = 'READ',
        extra_data_source = 'crc_lee_smc'
    )
	
)

set.seed(8542)

sc_tcga_deconv_comparison <- sapply(
	ct_to_keep,
	function(ct) {
		cat(paste0(ct, '\n'))
        sc_data <- eval(sc_metadata[[tcga_sc_map[ct]]]$read_quote)
		sapply(
			c('cancer', 'caf'),
			function(x) {
				
				plot_data <- copy(sc_data[cell_type == x])[
					,
					complexity := apply(.SD, 1, function(x) sum(x > 0)),
					.SDcols = -c('id', 'patient', 'cell_type')
				][
					,
					.SD[sample(1:.N, sample_sizes[[tcga_sc_map[ct]]], prob = complexity)]
				][
					,
					complexity := NULL
				][
					,
					c('id', deconv_data[[ct]]$genes_filtered[deconv_data[[ct]]$genes_filtered %in% names(sc_data)]),
					with = FALSE
				]
				
				plot_data <- melt(plot_data, id.var = 'id', variable.name = 'gene', value.name = 'expression_level')
				
				ordered_cell_ids <- plot_data[
					gene %in% with(
						deconv_data[[ct]],
						do.call(c('head', 'tail')[which(c('cancer', 'caf') == x)], list(genes_filtered[ordering], 20))
					),
					.(top_20_mean = mean(expression_level)),
					by = id
				][order(top_20_mean), id]
				
				htmp <- ggplot(
					plot_data,
					aes(
						x = factor(gene, levels = with(deconv_data[[ct]], genes_filtered[ordering])),
						y = factor(id, levels = ordered_cell_ids),
						fill = expression_level
					)
				) +
					geom_raster() +
					scale_x_discrete(expand = c(0, 0)) +
					scale_y_discrete(expand = c(0, 0)) +
					scale_fill_gradientn(
						colours = rev(colorRampPalette(brewer.pal(11, "Spectral"))(50)),
						limits = c(0, 12),
						oob = scales::squish,
						breaks = c(0, 3, 6, 9, 12),
						labels = c('0' = '0', '3' = '3', '6' = '6', '9' = '9', '12' = '\u2265 12')
					) +
					theme(
						axis.text = element_blank(),
						axis.title.x = element_blank(),
						axis.ticks = element_blank(),
						axis.ticks.length = unit(0, 'pt'),
						plot.margin = unit(c(5.5, 5.5, 1.25, 1.1), 'pt')
					) +
					labs(
						y = mapvalues(x, c('cancer', 'caf'), c('Cancer', 'CAF'), warn_missing = FALSE),
						fill = 'Expression\nlevel'
					)
				
				ave_exp_bar <- ggplot(
					plot_data[, .(ave_exp = mean(expression_level)), keyby = gene][
						with(deconv_data[[ct]], genes_filtered[ordering]),
						ave_exp := runmean(ave_exp, 30)
					],
					aes(x = factor(gene, levels = with(deconv_data[[ct]], genes_filtered[ordering])), y = 'a', fill = ave_exp)
				) +
					geom_raster() +
					scale_x_discrete(expand = c(0, 0)) +
					scale_y_discrete(expand = c(0, 0)) +
					scale_fill_gradientn(
						colours = colorRampPalette(brewer.pal(9, "RdPu"))(50),
						limits = c(0, 4),
						oob = scales::squish,
						breaks = c(0, 2, 4),
						labels = c('0' = '0', '2' = '2', '4' = '\u2265 4')
					) +
					theme(
						axis.text = element_blank(),
						axis.title = element_blank(),
						axis.ticks = element_blank(),
						axis.ticks.length = unit(0, 'pt'),
						plot.margin = switch(
							(x == 'cancer') + 1,
							unit(c(1.25, 5.5, 5.5, 5.5), 'pt'),
							unit(c(1.25, 5.5, 1.25, 5.5), 'pt')
						)
					) +
					labs(
						y = mapvalues(x, c('cancer', 'caf'), c('Cancer', 'CAF'), warn_missing = FALSE),
						fill = 'Average\nexpression\nlevel'
					)
				
				list(heatmap = htmp, ave_exp_bar = ave_exp_bar, plot_data = plot_data, ordered_cell_ids = ordered_cell_ids)
				
			},
			simplify = FALSE,
			USE.NAMES = TRUE
		)
	},
	simplify = FALSE,
	USE.NAMES = TRUE
)

# cairo_pdf('../data_and_figures/sc_deconv_comp.pdf', width = 6, height = 10, onefile = TRUE)

# for(ct in names(sc_tcga_deconv_comparison)) {
	
	# plot_grid(
		# plot_grid(
			# plotlist = c(
				# list(blank_plot() + theme(plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), 'pt')) + labs(title = deconv_args_per_ct[[ct]]$plot_title)),
				# lapply(
					# deconv_plots[[ct]]$plots[c('purity_bar', 'ccle_bar', 'heatmap')],
					# function(x) x + theme(legend.position = 'none', plot.title = element_blank())
				# ),
				# unlist(
					# lapply(
						# sc_tcga_deconv_comparison[[ct]],
						# function(x) list(x$heatmap + theme(legend.position = 'none'), x$ave_exp_bar + theme(legend.position = 'none'))
					# ),
					# recursive = FALSE
				# )
			# ),
			# nrow = 8,
			# ncol = 1,
			# align = 'v',
			# rel_heights = c(1, 1, 1, 15, 5, 1, 5, 1)
		# ),
		# plot_grid(
			# plotlist = c(
				# list(
					# blank_plot(),
					# get_legend(
						# deconv_plots[[1]]$plots$purity_bar +
							# theme(legend.justification = c(0, 1), legend.direction = 'vertical', legend.key.width = NULL, legend.key.height = NULL) +
							# labs(fill = 'Correlation\nwith purity')
					# ),
					# get_legend(
						# deconv_plots[[1]]$plots$ccle_bar +
							# theme(legend.justification = c(0, 1), legend.direction = 'vertical', legend.key.width = NULL, legend.key.height = NULL) +
							# labs(fill = 'Tumours vs.\ncell lines')
					# ),
					# get_legend(
						# deconv_plots[[1]]$plots$heatmap + theme(legend.justification = c(0, 1))
					# )
				# ),
				# list(
					# get_legend(sc_tcga_deconv_comparison[[1]][[1]]$heatmap + theme(legend.justification = c(0, 1))),
					# get_legend(sc_tcga_deconv_comparison[[1]][[1]]$ave_exp_bar + theme(legend.justification = c(0, 1)))
				# )
			# ),
			# nrow = 6,
			# ncol = 1,
			# rel_heights = c(1, 5, 5, 7, 6, 6)
		# ),
		# nrow = 1,
		# ncol = 2,
		# rel_widths = c(5, 1)
	# ) %>% print
	
# }

# dev.off()

# Make centred versions:

sc_tcga_deconv_comparison_centred <- sapply(
	names(sc_tcga_deconv_comparison),
	function(ct) {
		
		cat(paste0(ct, '\n'))
		
		plot_data <- rbind(
			sc_tcga_deconv_comparison[[ct]]$cancer$plot_data[, cell_type := 'cancer'],
			sc_tcga_deconv_comparison[[ct]]$caf$plot_data[, cell_type := 'caf']
		)
		
		# Calculate averages in cancer cells and CAFs separately, so we can check them afterwards:
		gene_averages_cancer_caf <- plot_data[
			,
			.(ave_exp = mean(expression_level)),
			by = .(gene, cell_type)
		]
		
		# To centre genes w.r.t. the average of the averages of cancer and CAF:
		gene_averages <- plot_data[
			,
			.(ave_exp = mean(c(mean(expression_level[cell_type == 'cancer']), mean(expression_level[cell_type == 'caf'])))),
			by = .(symbol = gene)
		]
		
		plot_data[, expression_level := expression_level - gene_averages[symbol == unique(gene), ave_exp], by = gene]
		
		# To centre the cells:
		plot_data[, expression_level := expression_level - mean(expression_level), by = id]
		
		htmps <- sapply(
			c('cancer', 'caf'),
			function(x) {
				ggplot(
					plot_data[cell_type == x],
					aes(
						x = factor(gene, levels = with(deconv_data[[ct]], genes_filtered[ordering])),
						y = factor(id, levels = sc_tcga_deconv_comparison[[ct]][[x]]$ordered_cell_ids),
						fill = expression_level
					)
				) +
					geom_raster() +
					scale_x_discrete(expand = c(0, 0)) +
					scale_y_discrete(expand = c(0, 0)) +
					scale_fill_gradientn(
						# colours = c('turquoise4', 'turquoise', 'azure', 'gold2', 'gold4')
						colours = rev(colorRampPalette(brewer.pal(11, "RdBu"))(50)),
						limits = c(-6, 6),
						# limits = c(-2, 2),
						oob = scales::squish,
						breaks = c(-6, -3, 0, 3, 6),
						labels = c('-6' = '\u2264 -6', '-3' = '-3', '0' = '0', '3' = '3', '6' = '\u2265 6')
						# breaks = c(-2, -1, 0, 1, 2),
						# labels = c('-2' = '\u2264 -2', '-1' = '-1', '0' = '0', '1' = '1', '2' = '\u2265 2')
					) +
					theme(
						axis.text = element_blank(),
						axis.title.x = element_blank(),
						axis.ticks = element_blank(),
						axis.ticks.length = unit(0, 'pt'),
						plot.margin = unit(c(5.5, 5.5, 1.25, 1.1), 'pt')
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
		
		list(heatmaps = htmps, plot_data = plot_data, gene_averages = gene_averages, gene_averages_cancer_caf = gene_averages_cancer_caf)
		
	},
	simplify = FALSE,
	USE.NAMES = TRUE
)

# cairo_pdf('../data_and_figures/sc_deconv_comp_centred.pdf', width = 6, height = 9, onefile = TRUE)

# for(ct in names(sc_tcga_deconv_comparison_centred)) {
	
	# plot_grid(
		# plot_grid(
			# plotlist = c(
				# list(blank_plot() + theme(plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), 'pt')) + labs(title = deconv_args_per_ct[[ct]]$plot_title)),
				# lapply(
					# deconv_plots[[ct]]$plots[c('purity_bar', 'ccle_bar', 'heatmap')],
					# function(x) x + theme(legend.position = 'none', plot.title = element_blank())
				# ),
				# lapply(sc_tcga_deconv_comparison_centred[[ct]]$heatmaps, function(x) x + theme(legend.position = 'none'))
			# ),
			# nrow = 6,
			# ncol = 1,
			# align = 'v',
			# rel_heights = c(1, 1, 1, 15, 5, 5)
		# ),
		# plot_grid(
			# plotlist = c(
				# list(
					# blank_plot(),
					# get_legend(
						# deconv_plots[[1]]$plots$purity_bar +
							# theme(legend.justification = c(0, 1), legend.direction = 'vertical', legend.key.width = NULL, legend.key.height = NULL) +
							# labs(fill = 'Correlation\nwith purity')
					# ),
					# get_legend(
						# deconv_plots[[1]]$plots$ccle_bar +
							# theme(legend.justification = c(0, 1), legend.direction = 'vertical', legend.key.width = NULL, legend.key.height = NULL) +
							# labs(fill = 'Tumours vs.\ncell lines')
					# ),
					# get_legend(
						# deconv_plots[[1]]$plots$heatmap + theme(legend.justification = c(0, 1))
					# ),
					# get_legend(sc_tcga_deconv_comparison_centred[[1]]$heatmaps[[1]] + theme(legend.justification = c(0, 1)))
				# )
			# ),
			# nrow = 5,
			# ncol = 1,
			# rel_heights = c(1, 5, 5, 7, 10)
		# ),
		# nrow = 1,
		# ncol = 2,
		# rel_widths = c(5, 1)
	# ) %>% print
	
# }

# dev.off()

# Make versions where we filter out genes with low expression in the single cell data:

sc_tcga_deconv_comparison_filtered <- sapply(
	names(sc_tcga_deconv_comparison_centred),
	function(ct) {
		
		cat(paste0(ct, '\n'))
		
		deconv_ct <- do.call(
			deconvolve_emt_caf_data,
			args = c(
				list(
					expression_data = expression_data,
					meta_data = meta_data,
					genes = sc_tcga_deconv_comparison_centred[[ct]]$gene_averages_cancer_caf[
						,
						.(pass = ave_exp[cell_type == 'cancer'] > 0.5 | ave_exp[cell_type == 'caf'] > 0.5),
						by = gene
					][pass == TRUE, as.character(gene)],
					# genes = sc_tcga_deconv_comparison_centred[[ct]]$gene_averages[ave_exp >= 0.5, as.character(symbol)], # symbol is a factor for some reason
					cell_type_markers = cell_type_markers,
					ccle_data = ccle,
					subtypes_data = subtypes_data,
					genes_from_tcga_fun = NULL,
					genes_filter_fun = NULL
				),
				deconv_args_per_ct[[ct]][!(names(deconv_args_per_ct[[ct]]) %in% c('plot_title', 'genes_filter_fun'))]
			)
		)
		
		deconv_ct <- deconv_reorder(deconv_ct)
		
		deconv_plots_ct <- do.call(
			deconvolve_emt_caf_plots,
			args = c(
				list(
					data = deconv_ct,
					heatmap_legend_title = 'Correlation',
					heatmap_colours = rev(colorRampPalette(brewer.pal(11, "RdBu"))(50)),
					heatmap_colour_limits = c(-0.6, 0.6),
					heatmap_legend_breaks = c(-0.6, -0.3, 0, 0.3, 0.6),
					heatmap_annotations_nudge = 0.3,
					purity_colours = rev(colorRampPalette(brewer.pal(11, "PuOr"))(50)),
					purity_colour_limits = c(-0.3, 0.3),
					purity_legend_breaks = c(-0.3, 0, 0.3),
					purity_legend_title = 'Correlation with purity\n',
					purity_legend_direction = 'horizontal',
					purity_axis_title = NULL,
					ccle_colours = rev(colorRampPalette(brewer.pal(11, 'PiYG')[c(1:3, 6, 9:11)])(50)),
					ccle_legend_breaks = c(-1, 0, 1),
					ccle_legend_title = 'Tumours vs. cell lines\n',
					ccle_legend_direction = 'horizontal',
					ccle_axis_title = NULL,
					bar_legend_justification = 'left'
				),
				deconv_args_per_ct[[ct]]['plot_title']
			)
		)
		
		plot_data <- sc_tcga_deconv_comparison_centred[[ct]]$plot_data[
			gene %in% sc_tcga_deconv_comparison_centred[[ct]]$gene_averages_cancer_caf[
				,
				.(pass = ave_exp[cell_type == 'cancer'] > 0.5 | ave_exp[cell_type == 'caf'] > 0.5),
				by = gene
			][pass == TRUE, as.character(gene)]
			# gene %in% sc_tcga_deconv_comparison_centred[[ct]]$gene_averages[ave_exp >= 0.5, symbol]
		]
		
		sc_heatmaps <- sapply(
			c('cancer', 'caf'),
			function(x) {
				ggplot(
					plot_data[cell_type == x],
					aes(
						x = factor(gene, levels = with(deconv_ct, genes_filtered[ordering])),
						y = factor(id, levels = sc_tcga_deconv_comparison[[ct]][[x]]$ordered_cell_ids),
						fill = expression_level
					)
				) +
					geom_raster() +
					scale_x_discrete(expand = c(0, 0)) +
					scale_y_discrete(expand = c(0, 0)) +
					scale_fill_gradientn(
						colours = rev(colorRampPalette(brewer.pal(11, "RdBu"))(50)),
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
					) +
					labs(y = mapvalues(x, c('cancer', 'caf'), c('Cancer', 'CAF'), warn_missing = FALSE), fill = 'Relative\nexpression\nlevel')
			},
			simplify = FALSE,
			USE.NAMES = TRUE
		)
		
		list(deconv_data = deconv_ct, deconv_figures = deconv_plots_ct$plots, sc_heatmaps = sc_heatmaps, plot_data = plot_data)
		
	},
	simplify = FALSE,
	USE.NAMES = TRUE
)

cairo_pdf('../data_and_figures/sc_deconv_comp_filtered.pdf', width = 6, height = 9, onefile = TRUE)

for(ct in names(sc_tcga_deconv_comparison_filtered)) {
	
	plot_grid(
		plot_grid(
			plotlist = c(
				list(blank_plot() + theme(plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), 'pt')) + labs(title = deconv_args_per_ct[[ct]]$plot_title)),
				lapply(
					c(
						sc_tcga_deconv_comparison_filtered[[ct]]$deconv_figures[c('purity_bar', 'ccle_bar', 'heatmap')],
						sc_tcga_deconv_comparison_filtered[[ct]]$sc_heatmaps
					),
					function(x) x + theme(legend.position = 'none', plot.title = element_blank())
				)
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
						deconv_plots[[1]]$plots$purity_bar +
							theme(legend.justification = c(0, 1), legend.direction = 'vertical', legend.key.width = NULL, legend.key.height = NULL) +
							labs(fill = 'Correlation\nwith purity')
					),
					get_legend(
						deconv_plots[[1]]$plots$ccle_bar +
							theme(legend.justification = c(0, 1), legend.direction = 'vertical', legend.key.width = NULL, legend.key.height = NULL) +
							labs(fill = 'Tumours vs.\ncell lines')
					),
					get_legend(
						deconv_plots[[1]]$plots$heatmap + theme(legend.justification = c(0, 1))
					),
					get_legend(sc_tcga_deconv_comparison_filtered[[1]]$sc_heatmaps[[1]] + theme(legend.justification = c(0, 1)))
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

# Make versions where we also filter out genes that don't fit the trend:

anomalous_threshold_tcga <- c(
	brca_luminal_a = 0.4,
	brca_luminal_b = 0.3,
	brca_basal_like = 0.3,
	brca_her2_enriched = 0.25,
	coad = 0.5,
	hnsc_mesenchymal_basal = 0.3,
	hnsc_classical = 0.3,
	hnsc_atypical = 0.5,
	lihc = 0.4,
	luad_proximal_inflammatory = 0.4,
	luad_proximal_proliferative = 0.5,
	lusc_basal = 0.3,
	lusc_classical = 0.6,
	lusc_primitive = 0.4,
	lusc_secretory = 0.5,
	ov_differentiated = 0.3,
	ov_immunoreactive = 0.5,
	ov_proliferative = 0.3,
	paad = 0.6,
	read = 0.3
)

sc_tcga_deconv_comparison_filtered_twice <- sapply(
	names(sc_tcga_deconv_comparison_centred),
	function(ct) {
		
		cat(paste0(ct, '\n'))
		
		deconv_ct <- do.call(
			deconvolve_emt_caf_data,
			args = c(
				list(
					expression_data = expression_data,
					meta_data = meta_data,
					genes = sc_tcga_deconv_comparison_centred[[ct]]$gene_averages_cancer_caf[
						,
						.(pass = ave_exp[cell_type == 'cancer'] > 0.5 | ave_exp[cell_type == 'caf'] > 0.5),
						by = gene
					][pass == TRUE, as.character(gene)],
					# genes = sc_tcga_deconv_comparison_centred[[ct]]$gene_averages[ave_exp >= 0.5, as.character(symbol)], # symbol is a factor for some reason
					cell_type_markers = cell_type_markers,
					ccle_data = ccle,
					subtypes_data = subtypes_data,
					genes_from_tcga_fun = NULL,
					genes_filter_fun = NULL
				),
				deconv_args_per_ct[[ct]][!(names(deconv_args_per_ct[[ct]]) %in% c('plot_title', 'genes_filter_fun'))]
			)
		)
		
		# Calculate new scores as in deconv_reorder() function.  Genes with positive score should be on the left side and those with
		# negative score should be on the right side.
		head_genes <- with(deconv_ct, head(genes_filtered[ordering], 20))
		tail_genes <- with(deconv_ct, tail(genes_filtered[ordering], 20))
		new_scores <- sapply(
			colnames(deconv_ct$cor_mat),
			function(g) {
				meanvec <- c(mean(deconv_ct$cor_mat[, g][head_genes]), mean(deconv_ct$cor_mat[, g][tail_genes]))
				ifelse(
					sign(meanvec[1]) == sign(meanvec[2]),
					ifelse(g %in% c('SNAI1', 'SNAI2', 'TWIST1', 'VIM','ZEB1', 'ZEB2'), return(0), return(NA)),
					return(max(meanvec)*c(1, -1)[which.max(meanvec)])
				)
			}
		)
		new_scores <- new_scores[!is.na(new_scores)]
		
		# For each gene, calculate the difference between its average across cancer cells and its average across CAFs.  For genes with
		# positive score, this difference should be positive, and for those with negative score it should be negative.  There will be
		# some exceptions - we need to decide if an exception is "extreme" enough.
		exp_diff_table <- sc_tcga_deconv_comparison_centred[[ct]]$plot_data[
			,
			.(exp_diff = mean(expression_level[cell_type == 'cancer']) - mean(expression_level[cell_type == 'caf'])),
			by = gene
		]
		
		# For each gene, consider the distribution of these differences for the genes having score with opposite sign.  If the gene is
		# beyond a certain quantile of this distribution ("beyond" w.r.t. its "own" distribution) then it is "anomalous".
		new_scores <- sapply(
			names(new_scores),
			function(g) {
				ifelse(
					new_scores[g] > 0,
					ifelse(
						exp_diff_table[
							,
							exp_diff[gene == g] < quantile(exp_diff[gene %in% names(new_scores)[new_scores < 0]], anomalous_threshold_tcga[ct])
						] & !(g %in% c('SNAI1', 'SNAI2', 'TWIST1', 'VIM', 'ZEB1', 'ZEB2')),
						return(NA),
						return(setNames(new_scores[g], NULL))
					),
					ifelse(
						new_scores[g] < 0,
						ifelse(
							exp_diff_table[
								,
								exp_diff[gene == g] > quantile(exp_diff[gene %in% names(new_scores)[new_scores > 0]], 1 - anomalous_threshold_tcga[ct])
							] & !(g %in% c('SNAI1', 'SNAI2', 'TWIST1', 'VIM', 'ZEB1', 'ZEB2')),
							return(NA),
							return(setNames(new_scores[g], NULL))
						),
						return(0)
					)
				)
				# exp_diff_table[, exp_diff[names(new_scores)[sign(new_scores) != sign(new_scores[g])]]
			}
		)
		
		anomalous_genes <- names(new_scores)[is.na(new_scores)]
		new_scores <- new_scores[!is.na(new_scores)]
		
		deconv_ct$anomalous_genes <- anomalous_genes
		deconv_ct$new_scores <- new_scores
		deconv_ct$cor_mat <- deconv_ct$cor_mat[names(new_scores), names(new_scores)]
		deconv_ct$genes_filtered <- names(new_scores)
		deconv_ct$ordering <- order(-new_scores)
		deconv_ct$cor_with_purity <- sapply(
			deconv_ct$cor_with_purity,
			function(x) x[names(new_scores)],
			simplify = FALSE,
			USE.NAMES = TRUE
		)
		deconv_ct$ccle_comp_diff <- deconv_ct$ccle_comp_diff[names(new_scores)]
		
		deconv_plots_ct <- do.call(
			deconvolve_emt_caf_plots,
			args = c(
				list(
					data = deconv_ct,
					heatmap_legend_title = 'Correlation',
					heatmap_colours = rev(colorRampPalette(brewer.pal(11, "RdBu"))(50)),
					heatmap_colour_limits = c(-0.6, 0.6),
					heatmap_legend_breaks = c(-0.6, -0.3, 0, 0.3, 0.6),
					heatmap_annotations_nudge = 0.3,
					purity_colours = rev(colorRampPalette(brewer.pal(11, "PuOr"))(50)),
					purity_colour_limits = c(-0.3, 0.3),
					purity_legend_breaks = c(-0.3, 0, 0.3),
					purity_legend_title = 'Correlation with purity\n',
					purity_legend_direction = 'horizontal',
					purity_axis_title = NULL,
					ccle_colours = rev(colorRampPalette(brewer.pal(11, 'PiYG')[c(1:3, 6, 9:11)])(50)),
					ccle_legend_breaks = c(-1, 0, 1),
					ccle_legend_title = 'Tumours vs. cell lines\n',
					ccle_legend_direction = 'horizontal',
					ccle_axis_title = NULL,
					bar_legend_justification = 'left'
				),
				deconv_args_per_ct[[ct]]['plot_title']
			)
		)
		
		plot_data <- sc_tcga_deconv_comparison_centred[[ct]]$plot_data[
			gene %in% sc_tcga_deconv_comparison_centred[[ct]]$gene_averages_cancer_caf[
				,
				.(pass = ave_exp[cell_type == 'cancer'] > 0.5 | ave_exp[cell_type == 'caf'] > 0.5),
				by = gene
			][pass == TRUE, as.character(gene)]
			# gene %in% sc_tcga_deconv_comparison_centred[[ct]]$gene_averages[ave_exp >= 0.5, symbol]
		]
		
		sc_heatmaps <- sapply(
			c('cancer', 'caf'),
			function(x) {
				ggplot(
					plot_data[cell_type == x],
					aes(
						x = factor(gene, levels = with(deconv_ct, genes_filtered[ordering])),
						y = factor(id, levels = sc_tcga_deconv_comparison[[ct]][[x]]$ordered_cell_ids),
						fill = expression_level
					)
				) +
					geom_raster() +
					scale_x_discrete(expand = c(0, 0)) +
					scale_y_discrete(expand = c(0, 0)) +
					scale_fill_gradientn(
						colours = rev(colorRampPalette(brewer.pal(11, "RdBu"))(50)),
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
					) +
					labs(y = mapvalues(x, c('cancer', 'caf'), c('Cancer', 'CAF'), warn_missing = FALSE), fill = 'Relative\nexpression\nlevel')
			},
			simplify = FALSE,
			USE.NAMES = TRUE
		)
		
		list(deconv_data = deconv_ct, deconv_figures = deconv_plots_ct$plots, sc_heatmaps = sc_heatmaps, plot_data = plot_data)
		
	},
	simplify = FALSE,
	USE.NAMES = TRUE
)

cairo_pdf('../data_and_figures/sc_deconv_comp_filtered_twice.pdf', width = 6, height = 9, onefile = TRUE)

for(ct in names(sc_tcga_deconv_comparison_filtered_twice)) {
	
	plot_grid(
		plot_grid(
			plotlist = c(
				list(blank_plot() + theme(plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), 'pt')) + labs(title = deconv_args_per_ct[[ct]]$plot_title)),
				lapply(
					c(
						sc_tcga_deconv_comparison_filtered_twice[[ct]]$deconv_figures[c('purity_bar', 'ccle_bar', 'heatmap')],
						sc_tcga_deconv_comparison_filtered_twice[[ct]]$sc_heatmaps
					),
					function(x) x + theme(legend.position = 'none', plot.title = element_blank())
				)
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
						deconv_plots[[1]]$plots$purity_bar +
							theme(legend.justification = c(0, 1), legend.direction = 'vertical', legend.key.width = NULL, legend.key.height = NULL) +
							labs(fill = 'Correlation\nwith purity')
					),
					get_legend(
						deconv_plots[[1]]$plots$ccle_bar +
							theme(legend.justification = c(0, 1), legend.direction = 'vertical', legend.key.width = NULL, legend.key.height = NULL) +
							labs(fill = 'Tumours vs.\ncell lines')
					),
					get_legend(
						deconv_plots[[1]]$plots$heatmap + theme(legend.justification = c(0, 1))
					),
					get_legend(sc_tcga_deconv_comparison_filtered_twice[[1]]$sc_heatmaps[[1]] + theme(legend.justification = c(0, 1)))
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
