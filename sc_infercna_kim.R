# bsub -q tirosh -n 4 -R "rusage[mem=32000]" -o sc_infercna_kim_log.o -e sc_infercna_kim_log.e Rscript sc_infercna_kim.R

cat(Sys.time(), '\n')

# The following is just to check the R version:
cat(R.Version()$version.string, '\n')

library(data.table) # 1.12.8
library(ggplot2) # This is version 3.3.1 on my laptop's WSL setup but 3.3.0 on WEXAC.
library(magrittr) # 1.5
library(infercna) # 1.0.0
library(stringr) # 1.4.0

source('general_functions.R')

# Function from infercna:
.chromosomeBreaks = function(genes = NULL, halfway = F, hide = NULL) {
    n = length(genes)
    chrsum = cumsum(lengths(splitGenes(genes, by = 'chr')))
    Breaks = chrsum/max(chrsum) * n
    if (halfway) {
        b = Breaks
        Breaks = c(b[1]/2, b[2:length(b)] - diff(b)/2)
    }
    if (!is.null(hide)) {
        names(Breaks)[which(names(Breaks) %in% hide)] = rep('', length(hide))
    }
    Breaks
}

cohort_data <-list(
	luad_kim = list(
		read_quote = quote(fread('../data_and_figures/kim_luad_2020.csv')[, -c('cell_type', 'cell_type_refined', 'cell_subtype', 'sample_site', 'sample_id', 'sample_origin')]),
		ref_cell_clusters = list(
			b_cell = c(7, 8, 11), # 17 could also be B cells, but they also have T cell signal, so could be doublets
			b_plasma = 13,
			dc = 23,
			endothelial = 21,
			macrophage = c(1, 19),
			mast = 5,
			t_cell = 6 # As above, 17 could be T cells, and possibly also 22, but signal in 22 is weak
		)
	)
)

cohort <- 'luad_kim'

cat('\tReading in data\n')

sc_data <- eval(cohort_data[[cohort]]$read_quote)

sc_dbscan <- readRDS(paste0('../data_and_figures/dbscan_', cohort, '.rds'))

ref_cells <- sapply(
	cohort_data[[cohort]]$ref_cell_clusters,
	function(x) sc_data$id[sc_dbscan$cluster %in% x],
	simplify = FALSE,
	USE.NAMES = TRUE
)

# Take genes with log average TPM at least 4:

cat('\tFiltering genes\n')

gene_averages <- sapply(
	sc_data[, -c('id', 'patient')],
	function(x) {log2(mean(10*(2^x - 1)) + 1)},
	USE.NAMES = TRUE
)

sc_data <- sc_data[, c('id', 'patient', names(gene_averages)[gene_averages >= 4]), with = FALSE]

# Running infercna() function:

cat('\tInferring CNAs\n')

inferred_cna <- sapply(
	as.character(unique(sc_data$patient)),
	function(p) {
		infercna(
			t(
				set_rownames(
					as.matrix(sc_data[patient == p | id %in% unlist(ref_cells), -c('id', 'patient')]),
					sc_data[patient == p | id %in% unlist(ref_cells), id]
				)
			),
			refCells = ref_cells,
			isLog = TRUE
		)
	},
	simplify = FALSE,
	USE.NAMES = TRUE
)

cat('\tSaving CNA matrices\n')

saveRDS(inferred_cna, paste0('../data_and_figures/inferred_cna_', cohort, '.rds'))

cat('\tMaking plots\n')

# The genes should be the same in all of the inferred_cna matrices:
x_line_breaks <- .chromosomeBreaks(rownames(inferred_cna[[1]]))
x_text_breaks <- round(.chromosomeBreaks(rownames(inferred_cna[[1]]), halfway = TRUE, hide = c('13', '18', '21', 'Y')))

clust_plots <- sapply(
	names(inferred_cna),
	function(p) {
		
		cell_ids <- colnames(inferred_cna[[p]])[!(colnames(inferred_cna[[p]]) %in% unlist(ref_cells))]
		
		# The following was necessary for the liver HCC data from Ma et al., in which patient H18 had no non-reference cells
		# according to my classification by t-SNE.
		if(length(cell_ids) < 2) {
			warning('Not enough non-reference cells in this tumour!')
			return(NULL)
		}
		
		cell_clust <- hclust(dist(t(inferred_cna[[p]][, cell_ids])))
		
		cna_heatmap <- ggplot(
			reshape2::melt(
				inferred_cna[[p]][, cell_ids],
				varnames = c('gene', 'cell'),
				value.name = 'cna_score'
			)
		) +
			geom_raster(
				aes(
					x = factor(gene, levels = rownames(inferred_cna[[p]])),
					y = factor(cell, levels = cell_ids[cell_clust$order]),
					fill = cna_score
				)
			) +
			scale_fill_gradientn(
				colours = rev(colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))(50)),
				limits = c(-1, 1),
				oob = scales::squish
			) +
			scale_y_discrete(expand = c(0, 0)) +
			scale_x_discrete(expand = c(0, 0), breaks = rownames(inferred_cna[[p]])[x_text_breaks], labels = names(x_text_breaks)) +
			geom_vline(xintercept = x_line_breaks, size = 0.25) +
			theme(
				axis.text.y = element_blank(),
				axis.ticks = element_blank(),
				axis.ticks.length = unit(0, 'pt'),
				panel.border = element_rect(colour = 'black', size = 0.5, fill = NA)
			) +
			labs(x = 'Chromosome', y = 'Cells', fill = 'Inferred CNA', title = p)
		
		# We don't need to return cell_ids, because these are stored in cell_clust$labels anyway.
		
		list(cell_clust = cell_clust, heatmap = cna_heatmap)
		
	},
	simplify = FALSE,
	USE.NAMES = TRUE
)

# for(x in clust_plots) saveRDS(x$cell_clust, paste0('../data_and_figures/cna_clust_', cohort, '.rds'))

saveRDS(
	sapply(clust_plots, `[[`, 'cell_clust', simplify = FALSE, USE.NAMES = TRUE),
	paste0('../data_and_figures/cna_clust_', cohort, '.rds')
)

pdf(paste0('../data_and_figures/cna_plots_', cohort, '.pdf'))
for(x in clust_plots) print(x$heatmap)
dev.off()

cat('\tDone!\n')
