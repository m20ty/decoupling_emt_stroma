library(data.table)
library(ggplot2)
library(egg)
library(limma)
library(cowplot)
library(stringr)
library(readxl)
library(plyr)
library(stringi)
library(org.Hs.eg.db)
library(AnnotationDbi)

source('general_functions.R')
source('tcga_functions.R')

# I need cell_type_markers for the gene filtering:
cell_type_markers <- fread('../../cell_type_markers.csv')[cell_type != 'mesenchymal' & source != 'TCGA_CCLE_comparison']
# cell_type_markers <- fread('../../cell_type_markers.csv')
emt_markers <- fread('../../emt_markers.csv')[, gene := alias2SymbolTable(gene)][source != 'GO', sort(unique(gene))]

heatmap_annotations <- c('ACTA2', 'AXL', 'COL1A1', 'COL1A2', 'COL3A1', 'COL4A1', 'COL4A2', 'COL5A1', 'COL5A2', 'COL5A3', 'COL6A1', 'COL6A2',
	'COL6A3', 'COL7A1', 'CD44', 'CDH2', 'ECM1', 'ECM2', 'FAP', 'FN1', 'IL6', 'ITGA2', 'ITGA5', 'ITGA6', 'ITGB1', 'ITGB3', 'ITGB5', 'ITGB6',
    'LAMA1', 'LAMA2', 'LAMA3', 'LAMA5', 'LAMB3', 'LAMC1', 'LAMC2', 'MMP1', 'MMP2', 'MMP3', 'MMP10', 'MMP14', 'PDPN', 'SNAI1', 'SNAI2', 'SPARC',
    'TGFB1', 'TGFBI', 'THY1', 'TNC', 'TWIST1', 'VCAN', 'VIM', 'ZEB1', 'ZEB2')

initial_genes <- c('SNAI1', 'SNAI2', 'TWIST1', 'VIM', 'ZEB1', 'ZEB2')





# Get UCS data:

if(all(c('tcga_expression_data_ucs.csv', 'tcga_meta_data_ucs.csv') %in% dir('../../TCGA_data'))) {
    expression_data <- fread('../../TCGA_data/tcga_expression_data_ucs.csv')
    meta_data <- fread('../../TCGA_data/tcga_meta_data_ucs.csv')
} else {

    if('UCS_illuminahiseq_rnaseqv2_RSEM_genes.txt' %in% dir('../../TCGA_data/UCS')) {
        expression_data <- fread('../../TCGA_data/UCS/UCS_illuminahiseq_rnaseqv2_RSEM_genes.txt', showProgress = FALSE)
    } else {

        download.file(
            paste0(
                'http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/UCS/20160128/gdac.broadinstitute.org_UCS',
                '.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.Level_3.2016012800.0.0.tar.gz'
            ),
            destfile = 'tmp.tar.gz',
            quiet = TRUE
        )

        file_names <- untar('tmp.tar.gz', list = TRUE)
        untar('tmp.tar.gz', files = file_names[endsWith(file_names, 'data.txt')], exdir = 'tmp')
        file.remove('tmp.tar.gz')

        # I have to go through the following renaming procedure because the file name seems too long for R to handle.

        file.rename(paste0('tmp/', str_split_fixed(file_names[endsWith(file_names, 'data.txt')], '/', 2)[, 1]), 'tmp/tempdir')

        # Move the data file up one directory and delete tempdir:
        file.rename(
            paste0('tmp/tempdir/', str_split_fixed(file_names[endsWith(file_names, 'data.txt')], '/', 2)[, 2]),
            paste0('tmp/', str_split_fixed(file_names[endsWith(file_names, 'data.txt')], '/', 2)[, 2])
        )

        unlink('tmp/tempdir', recursive = TRUE)

        # Read in the data:
        expression_data <- fread(
            paste0('tmp/', str_split_fixed(file_names[endsWith(file_names, 'data.txt')], '/', 2)[, 2]),
            showProgress = FALSE
        )

        # Remove the created directory:
        unlink('tmp', recursive = TRUE)

    }

    names(expression_data)[1] <- 'id'

    # Extract scaled estimates:
    expression_data <- expression_data[, expression_data[1] %in% c('gene_id', 'scaled_estimate'), with = FALSE][-1]

    # Make columns numeric:
    expression_data[, names(expression_data[, -'id']) := lapply(.SD, as.numeric), .SDcols = -'id']

    # The following obtains gene symbols by converting the Entrez IDs.  This mostly
    # agrees with the output of alias2SymbolTable() applied to the symbols already in
    # the id column, but in the few places where it disagrees, it seems more accurate
    # than alias2SymbolTable().  We remove NAs after the conversion.
    expression_data[, id := mapIds(org.Hs.eg.db, str_split_fixed(id, '\\|', 2)[, 2], keytype = 'ENTREZID', column = 'SYMBOL')]
    expression_data <- expression_data[!is.na(id)]

    # Transpose and split IDs for extraction of patient IDs and sample type codes:
    expression_data <- tdt(expression_data)
    expression_data$id <- gsub('-', '\\.', expression_data$id)

    # Make sure it's in alphabetical order:
    expression_data <- setkey(expression_data, id)

    sample_codes <- str_split_fixed(expression_data$id, '\\.', 5)
    rows_to_keep <- grep('^01|^02|^05|^06|^07|^11', sample_codes[, 4])
    expression_data <- expression_data[rows_to_keep]
    sample_codes <- sample_codes[rows_to_keep, ]

    # Convert to (something like) log TPM by multiplying by 1e+06 and taking log:
    expression_data[, names(expression_data[, -'id']) := log2(.SD*1e+06 + 1), .SDcols = -'id']

    # Construct metadata table:
    meta_data <- data.table(
        id = expression_data$id,
        patient_id = apply(sample_codes[, 1:3], 1, paste, collapse = '.'),
        cancer_type = 'UCS',
        sample_type = mapvalues(
            str_split_fixed(sample_codes[, 4], '[A-Z]', 2)[, 1],
            c('01', '02', '05', '06', '07', '11'),
            c('primary', 'recurrent', 'primary_additional', 'metastatic', 'metastatic_additional', 'normal'),
            warn_missing = FALSE
        )
    )

    # Get purity data:
    purity_data <- as.data.table(read_xlsx('../../TCGA_data/UCS/1-s2.0-S1535610817300533-mmc4.xlsx', skip = 1))[, .(patient_id = gsub('-', '\\.', sample), purity = purity)]
    meta_data <- merge(meta_data, purity_data, by = 'patient_id')

    setkey(expression_data, id)
    setkey(meta_data, id)

    # Write to file:
    fwrite(expression_data, '../../TCGA_data/tcga_expression_data_ucs.csv')
    fwrite(meta_data, '../../TCGA_data/tcga_meta_data_ucs.csv')

}





# This is pretty robust to changes in number of genes and to the filtering.  It doesn't make much difference if I use 90th percentile or median for gene_weights_fun.
# Putting cell_type_weights = FALSE makes a bit of a difference (especially visible in the diagnostics), but not a great deal.  Maybe try several parameters and plot
# all the results together, so we can decide which we prefer.  It looks like CAFs may be hard to recognise, and it could indeed be different mesenchymal components
# among the malignant cells.

deconv_ucs <- sapply(
	c(150, 200, 250),
	function(n) setNames(
		sapply(
			list(
				list(cell_type_weights = list(B_plasma = 0, myocyte = 1, macrophage = 0, endothelial = 0, DC = 0, mast = 0, T = 0, B = 0)),
				list(gene_weights_fun = function(x) quantile(x, 0.9)),
				list(gene_weights_fun = median),
				list(cell_type_weights = FALSE)
			),
			function(weights_arg) do.call(
				deconvolve_emt_caf_data,
				args = c(
					list(
						expression_data = expression_data,
					    meta_data = meta_data,
					    genes = emt_markers,
					    cell_type_markers = cell_type_markers,
					    initial_gene_weights = FALSE,
					    tcga_cancer_types = 'UCS',
						cell_types = c('B_plasma', 'myocyte', 'macrophage', 'endothelial', 'DC', 'mast', 'T', 'B'),
					    seed = 5522,
						genes_filter_fun = function(x) 1:n
					),
					weights_arg
				)
			),
		    simplify = FALSE,
		    USE.NAMES = TRUE
		),
		c('all on myocyte', '90th percentile', 'median', 'none')
	),
    simplify = FALSE,
    USE.NAMES = TRUE
) %>% setNames(c('n = 150', 'n = 200', 'n = 250'))

deconv_ucs <- unlist(deconv_ucs, recursive = FALSE)

names(deconv_ucs) <- gsub('\\.', ', ', names(deconv_ucs))

for(deconv_ct in deconv_ucs) {
	deconv_ct$cor_with_initial_and_cell_types[
		deconv_ct$genes_filtered,
		CAF := expression_data[
			,
			rowMeans(
				cor(
					.SD[, deconv_ct$genes_filtered, with = FALSE],
					.SD[, cell_type_markers[cell_type == 'CAF' & gene %in% names(expression_data), gene], with = FALSE]
				)
			)
		]
	]
}

for(ct in names(deconv_ucs)) {deconv_ucs[[ct]]$cell_types <- c('CAF', deconv_ucs[[ct]]$cell_types)}

# deconv_ucs <- deconvolve_emt_caf_data(
#     expression_data = expression_data,
#     meta_data = meta_data,
#     genes = emt_markers,
#     cell_type_markers = cell_type_markers,
#     initial_gene_weights = FALSE,
#     tcga_cancer_types = 'UCS',
#     seed = 5522,
#     genes_filter_fun = function(x) 1:250,
#     cell_type_weights = FALSE
#     # gene_weights_fun = median
#     # gene_weights_fun = function(x) quantile(x, 0.9)
# )

# deconv_ucs <- deconv_reorder(deconv_ucs)

deconv_ucs <- sapply(deconv_ucs, deconv_reorder, simplify = FALSE, USE.NAMES = TRUE)

deconv_plot_ucs <- sapply(
	names(deconv_ucs),
	function(deconv_name) do.call(
		deconvolve_emt_caf_plots,
		args = list(
			data = deconv_ucs[[deconv_name]],
		    heatmap_axis_title = '',
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
		    bar_legend_width = NULL,
		    bar_legend_height = NULL,
		    plot_title = paste0('UCS - ', deconv_name)
		)
	),
	simplify = FALSE,
	USE.NAMES = TRUE
)

# deconv_plot_ucs <- deconvolve_emt_caf_plots(
#     data = deconv_ucs,
#     heatmap_axis_title = '',
#     heatmap_legend_title = 'Coexpression',
#     heatmap_colours = rev(colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))(50)),
#     heatmap_annotations = heatmap_annotations,
#     heatmap_annotations_nudge = 0.3,
#     purity_colours = rev(colorRampPalette(RColorBrewer::brewer.pal(11, "PuOr"))(50)),
#     purity_colour_limits = c(-0.3, 0.3),
#     purity_legend_breaks = c(-0.3, 0, 0.3),
#     purity_legend_title = 'Correlation\nwith purity',
#     purity_legend_direction = 'vertical',
#     purity_axis_title = NULL,
#     bar_legend_width = NULL,
#     bar_legend_height = NULL,
#     plot_title = 'UCS'
# )

saveRDS(deconv_ucs, '../data_and_figures/deconv_ucs.rds')

pdf('../data_and_figures/deconv_ucs_figures.pdf', width = 10, height = 12)
for(figs in deconv_plot_ucs) print(deconv_plot(list(figs), legends_space = 0.2))
dev.off()

# Diagnostic figures:
pdf('../data_and_figures/deconv_ucs_diagnostics.pdf', width = 10, height = 12)
i <- 1
for(figs in deconv_plot_ucs) {
	ggarrange(
	    plots = c(list(figs$plots$heatmap), figs$diagnostics$alternative_purity_cor, figs$diagnostics$cell_type_bars),
	    ncol = 1,
	    nrow = 12,
	    heights = c(8, rep(0.6, 11)),
		newpage = switch((i == 1) + 1, TRUE, FALSE)
	)
	i <- i + 1
}
dev.off()





# Save one separately - not sure if it's the "best" one, but it's good and illustrates a certain point.  The point is that I can't really get rid of the correlation
# of the "cancer" end with myocytes (according to the diagnostics), so this end might represent a sarcoma-like component.  There are a couple of TPMs at this end,
# which would make sense, as these are muscle-related genes.  The right end looks more like the pEMT from squamous or GI cancer types, and more like what I think of
# as pEMT.  But, interestingly, some of the genes at the left end are found in the gynae cluster of my pEMT-CAF scores heatmap, including SACS, CAP2, TGFBR3 and QKI.
# What could this imply about the gynae cluster?  That it's more sarcoma-like?  There are a few gynae genes around the middle, though, so this isn't a strong trend.

pdf('../data_and_figures/deconv_ucs_final.pdf', width = 10, height = 12)
print(deconv_plot(list(deconv_plot_ucs[['n = 200, 90th percentile']]), legends_space = 0.2))
ggarrange(
	plots = c(
		list(deconv_plot_ucs[['n = 200, 90th percentile']]$plots$heatmap),
		# deconv_plot_ucs[['n = 200, 90th percentile']]$diagnostics$alternative_purity_cor,
		deconv_plot_ucs[['n = 200, 90th percentile']]$diagnostics$cell_type_bars
	),
	ncol = 1,
	nrow = 10,
	heights = c(8, rep(0.6, 9))
)
dev.off()
