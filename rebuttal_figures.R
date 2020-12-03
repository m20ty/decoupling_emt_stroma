# This script doesn't produce all the rebuttal figures: the code for making R1 can be found in sc_deconv_comp.R.

library(data.table) # 1.12.8
library(ggplot2) # This is version 3.3.1 on my laptop's WSL setup but 3.3.0 on WEXAC.
library(cowplot) # 1.0.0
library(magrittr) # 1.5
library(cowplot) # 1.0.0
library(Rtsne) # 0.15
library(plyr) # 1.8.6
library(stringr) # 1.4.0
library(colorspace) # 1.4.1
library(infercna) # 1.0.0
library(limma) # 3.42.2
library(caTools) # 1.18.0
library(RColorBrewer) # 1.1.2
library(scales) # 1.1.1
library(latex2exp) # 0.4.0
library(ggrepel) # 0.8.2
library(egg) # 0.4.5

source('general_functions.R')
source('tcga_functions.R')
source('sc_functions.R')

cell_type_markers <- fread('../../cell_type_markers.csv')[!(cell_type %in% c('mesenchymal', 'myocyte')) & source != 'TCGA_CCLE_comparison']

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

cell_type_labels <- c(
	'acinar' = 'Acinar',
	'alveolar' = 'Alveolar',
	'b_cell' = 'B cell',
	'caf' = 'CAF',
	'caf_potential' = 'Potential CAF',
	'cancer' = 'Cancer',
	'dendritic' = 'Dendritic',
	'ductal' = 'Ductal',
	'endocrine' = 'Endocrine',
	'endothelial' = 'Endothelial',
	'epithelial' = 'Epithelial',
	'erythroblast' = 'Erythroblast',
	'hpc-like' = 'HPC-like',
	'macrophage' = 'Macrophage',
	'mast' = 'Mast',
	'myocyte' = 'Myocyte',
	'nk_cell' = 'NK cell',
	'stellate' = 'Stellate',
	't_cell' = 'T cell'
)

cell_type_colours <- c(
	'acinar' = '#C0DF84',
	'alveolar' = '#E3E050',
	'b_cell' = '#887DDA',
	'caf' = '#D8A354',
	'caf_potential' = lighten('#D8A354', 0.5),
	'cancer' = '#DD4C61',
	'dendritic' = '#7AE1D8',
	'ductal' = '#DBE2BB',
	'endocrine' = '#7AE759',
	'endothelial' = '#8C5581',
	'epithelial' = '#D1D2DC',
	'erythroblast' = '#DCA9DB',
	'hpc-like' = '#A83EE5',
	'macrophage' = '#75B0DC',
	'mast' = '#799575',
	'myocyte' = '#DD6BCF',
	'nk_cell' = '#01AD1E',
	'stellate' = '#D5978A',
	't_cell' = '#79E5A3'
)





cohort_data <- list(
	luad_kim = list(
		read_quote = quote(fread('../data_and_figures/kim_luad_2020_reclassified.csv')),
		seed = 5481
	)
)

cohort <- 'luad_kim'

sc_data <- eval(cohort_data[[cohort]]$read_quote)[, -'cell_type_author']

sc_data <- sc_data[cell_type %in% c('caf', 'cancer')]

inferred_cna <- readRDS(paste0('../data_and_figures/inferred_cna_', cohort, '.rds'))
inferred_cna <- inferred_cna[names(inferred_cna) %in% sc_data$patient]

x_line_breaks <- .chromosomeBreaks(rownames(inferred_cna[[1]]))
x_text_breaks <- round(.chromosomeBreaks(rownames(inferred_cna[[1]]), halfway = TRUE, hide = c('13', '18', '21', 'Y')))

cna_signal_files <- dir(paste0('../data_and_figures/sc_find_malignant/', cohort))
cna_signal_files <- cna_signal_files[endsWith(cna_signal_files, '_data.csv')]

classification_data <- lapply(
    cna_signal_files,
    function(x) cbind(fread(paste0('../data_and_figures/sc_find_malignant/', cohort, '/', x)))#, patient = gsub('_data.csv', '', x))
) %>% rbindlist(fill = TRUE)

setkey(classification_data, cell_id)

classification_data <- classification_data[sc_data[cell_type == 'cancer', id]]
classification_data <- classification_data[classification == classification_final]

# For some reason, the key sometimes disappears after the above two lines, so I'm setting it again:
setkey(classification_data, cell_id)

classification_data <- sc_data[
    ,
    .(
        patient = patient,
        classification = switch((cell_type == 'cancer') + 1, 'caf', classification_data[id, classification_final])
    ),
    by = id
]

# The following deals with cases where all cells of one subclone get changed to 'ambiguous' at some point, presumably because their I disagree with
# the authors as to their malignancy.  This happened in HNSCC, where all cells in malignant clone 2 became 'ambiguous' leaving only 'Patient 6:
# clone 1' in the final heatmap, which obviously looks odd.
classification_data[
    startsWith(classification, 'malignant') & patient %in% classification_data[
        startsWith(classification, 'malignant'),
        .(dodgy = any(grepl('clone', classification)) & length(unique(classification)) == 1),
        by = patient
    ][dodgy == TRUE, patient],
    classification := 'malignant'
]

plot_data <- rbindlist(
    lapply(
        unique(classification_data$patient),
        function(p) cbind(classification_data[patient == p], t(inferred_cna[[as.character(p)]][, classification_data[patient == p, id]]))
    )
)

# Downsample where a patient is overrepresented:
# Note the final proportion won't be 2/#patients, because you're taking a sample of size 2*sum(N)/#patients, but then sum(N) becomes smaller.
set.seed(cohort_data[[cohort]]$seeds)
downsampled_ids <- plot_data[, .N, by = .(patient, classification)][
    ,
    setNames(lapply(list(patient, classification, N), function(x) x[N > 2*sum(N)/length(unique(patient))]), c('patient_id', 'cl_full', 'N')),
    by = .(cl = gsub(' clone [0-9]$', '', classification))
][
    ,
    .(
        sampled_ids = plot_data[
            startsWith(classification, cl),
            sample(id[patient == patient_id & classification == cl_full], 2*.N/length(unique(patient)))
        ]
    ),
    # plot_data[, sample(id[patient == patient_id], 3*sum(startsWith(classification, cl))/length(unique(patient)))],
    by = .(cl, patient_id, cl_full)
]

downsampled_ids <- plot_data[
    ,
    .(
        sampled_ids = switch(
            (patient %in% downsampled_ids$patient_id && classification %in% downsampled_ids[patient_id == patient, cl_full]) + 1,
            id,
            downsampled_ids[patient_id == patient & cl_full == classification, sampled_ids]
        )
    ),
    by = .(patient, classification)
]$sampled_ids

plot_data <- plot_data[id %in% downsampled_ids]
classification_data <- classification_data[id %in% downsampled_ids]

plot_data <- melt(
    plot_data,
    id.vars = c('id', 'patient', 'classification'),
    variable.name = 'gene',
    value.name = 'cna_score'
)
plot_data[, gene := factor(gene, levels = rownames(inferred_cna[[1]]))]

# For cancer cell heatmap:

y_breaks_cancer_major <- classification_data[classification != 'caf'][order(patient), .(brks = .N), by = patient]$brks %>% cumsum

y_breaks_cancer_minor <- numeric(0)
temp <- classification_data[classification != 'caf'][order(patient, classification)]
for(i in 2:nrow(temp)) {
    if(
        startsWith(temp[i - 1, classification], 'malignant clone') &&
            temp[i - 1, classification] != temp[i, classification] &&
            startsWith(temp[i, classification], 'malignant clone') &&
            temp[i - 1, patient] == temp[i, patient]
    ) {y_breaks_cancer_minor <- c(y_breaks_cancer_minor, i)}
}
rm(temp)

all_y_breaks <- sort(c(0, y_breaks_cancer_major, y_breaks_cancer_minor))
y_text_breaks_cancer <- sapply(2:length(all_y_breaks), function(i) round(all_y_breaks[i - 1] + (all_y_breaks[i] - all_y_breaks[i - 1])/2))

y_breaks_cancer_major <- y_breaks_cancer_major[1:(length(y_breaks_cancer_major) - 1)]

patient_annotations_cancer <- character(classification_data[classification != 'caf', .N])
patient_annotations_cancer[y_text_breaks_cancer] <- classification_data[classification != 'caf'][
    order(patient, classification),
    sapply(
        y_text_breaks_cancer,
        function(x) switch(
            (classification[x] == 'malignant') + 1,
            paste0('Patient ', patient[x], ':', gsub('malignant', '', classification[x])),
            paste0('Patient ', patient[x])
        )
    )
]

# For CAF heatmap:

y_breaks_caf <- classification_data[classification == 'caf'][order(patient), .(brks = .N), by = patient]$brks %>% cumsum

y_text_breaks_caf <- sapply(
    2:(length(y_breaks_caf) + 1),
    function(i) round(c(0, y_breaks_caf)[i - 1] + (c(0, y_breaks_caf)[i] - c(0, y_breaks_caf)[i - 1])/2)
)

y_breaks_caf <- y_breaks_caf[1:(length(y_breaks_caf) - 1)]

patient_annotations_caf <- character(classification_data[classification == 'caf', .N])
patient_annotations_caf[y_text_breaks_caf] <- classification_data[classification == 'caf'][
    order(patient, classification),
    sapply(y_text_breaks_caf, function(x) paste0('Patient ', patient[x]))
]

# For adjustment of chromosome 6p:
genes_split_by_arm <- splitGenes(rownames(inferred_cna[[1]]), by = 'arm')

cna_heatmap_cancer <- ggplot(plot_data[classification != 'caf']) +
    geom_raster(
        aes(
            x = gene,
            y = factor(id, levels = classification_data[classification != 'caf'][order(patient, classification, id), unique(id)]),
            fill = cna_score
        )
    ) +
    scale_fill_gradientn(
        colours = rev(colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))(50)),
        limits = c(-1, 1),
        oob = scales::squish
    ) +
    scale_y_discrete(
        expand = c(0, 0),
        breaks = classification_data[classification != 'caf'][order(patient, classification, id), id[y_text_breaks_cancer]],
        labels = classification_data[classification != 'caf'][
            order(patient, classification),
            sapply(
                y_text_breaks_cancer,
                function(x) switch(
                    (classification[x] == 'malignant') + 1,
                    paste0('Patient ', patient[x], ':', gsub('malignant', '', classification[x])),
                    paste0('Patient ', patient[x])
                )
            )
        ]
    ) +
    scale_x_discrete(
        expand = c(0, 0),
        breaks = plot_data[classification != 'caf', unique(gene)][x_text_breaks],
        labels = names(x_text_breaks)
    ) +
    geom_vline(xintercept = x_line_breaks, size = 0.25) +
    geom_hline(yintercept = y_breaks_cancer_major, size = 0.25) +
    geom_hline(yintercept = y_breaks_cancer_minor, size = 0.25, linetype = 'dashed') +
    theme(
        axis.ticks = element_blank(),
        axis.ticks.length = unit(0, 'pt'),
        panel.border = element_rect(colour = 'black', size = 0.5, fill = NA)
    ) +
    labs(x = 'Chromosome', y = 'Cells', fill = 'Inferred CNA')

annotations_cancer <- heat_map_labels_repel(
    patient_annotations_cancer,
    edge = 'right',
    nudge = 0.1,
    axis_title = 'Cancer cells',
    title_edge_margin = 5.5
)

cna_heatmap_caf <- ggplot(plot_data[classification == 'caf']) +
    geom_raster(
        aes(
            x = gene,
            y = factor(id, levels = classification_data[classification == 'caf'][order(patient, classification, id), unique(id)]),
            fill = cna_score
        )
    ) +
    scale_fill_gradientn(
        colours = rev(colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))(50)),
        limits = c(-1, 1),
        oob = scales::squish
    ) +
    scale_y_discrete(
        expand = c(0, 0),
        breaks = classification_data[classification == 'caf'][order(patient, classification, id), id[y_text_breaks_caf]],
        labels = classification_data[classification == 'caf'][
            order(patient, classification),
            sapply(y_text_breaks_caf, function(x) paste0('Patient ', patient[x]))
        ]
    ) +
    scale_x_discrete(
        expand = c(0, 0),
        breaks = plot_data[classification == 'caf', unique(gene)][x_text_breaks],
        labels = names(x_text_breaks)
    ) +
    geom_vline(xintercept = x_line_breaks, size = 0.25) +
    geom_hline(yintercept = y_breaks_caf, size = 0.25) +
    theme(
        axis.ticks = element_blank(),
        axis.ticks.length = unit(0, 'pt'),
        panel.border = element_rect(colour = 'black', size = 0.5, fill = NA),
        plot.title = element_text(margin = margin(b = 30), size = 18)
    ) +
    labs(x = 'Chromosome', y = 'Cells', fill = 'Inferred CNA', title = cohort_data[[cohort]]$plot_title)

annotations_caf <- heat_map_labels_repel(
    patient_annotations_caf,
    edge = 'right',
    nudge = 0.1,
    axis_title = 'CAFs',
    title_edge_margin = 5.5
)





cohort_data <- list(
    luad_kim = list(
		read_quote = quote(fread('../data_and_figures/kim_luad_2020_reclassified.csv')),
		seeds = c(4573, 2106, 7294)
	)
)

sc_data <- eval(cohort_data[[cohort]]$read_quote)[, -'cell_type_author']

sc_data <- sc_data[cell_type != 'ambiguous']

gene_averages <- sapply(
    sc_data[, -c('id', 'patient', 'cell_type')],
    function(x) {log2(mean(10*(2^x - 1)) + 1)},
    USE.NAMES = TRUE
)

sc_data <- sc_data[, c('id', 'patient', 'cell_type', names(gene_averages)[gene_averages >= 4]), with = FALSE]

# Run t-SNE:
# set.seed(cohort_data[[cohort]]$seeds[1])
# sc_tsne <- Rtsne(
    # as.matrix(sc_data[, lapply(.SD, function(x) {x - mean(x)}), .SDcols = -c('id', 'patient', 'cell_type')])#,
    # # num_threads = 16
# )

# saveRDS(sc_tsne, paste0('../data_and_figures/sc_tsne_final/tsne_', cohort, '.rds'))

sc_tsne <- readRDS(paste0('../data_and_figures/sc_tsne_final/tsne_', cohort, '.rds'))

if('cells_to_exclude' %in% names(cohort_data[[cohort]])) {
    sc_tsne$Y <- sc_tsne$Y[-which(sc_data$id %in% cohort_data[[cohort]]$cells_to_exclude), ]
    sc_data <- sc_data[!(id %in% cohort_data[[cohort]]$cells_to_exclude)]
}

# If cohort == 'pdac_peng' and we haven't removed T8_TGGTTCCTCGCATGGC and T17_CGTGTAACAGTACACT, we can see these are the T cell and CAF in the wrong
# clusters by the following:
# sc_data[(cell_type == 't_cell' & sc_tsne$Y[, 1] > 35) | (cell_type == 'caf' & sc_tsne$Y[, 1] < -35), id]

plot_data <- setNames(
    cbind(as.data.table(sc_tsne$Y), sc_data[, .(patient, cell_type)]),
    c('x', 'y', 'patient', 'cell_type')
)

tsne_plot_cell_types <- ggplot(plot_data, aes(x = x, y = y, colour = cell_type)) +
    geom_point(size = 0.7) +
    scale_colour_manual(
        labels = cell_type_labels[sort(unique(plot_data$cell_type))],
        values = cell_type_colours[sort(unique(plot_data$cell_type))]
    ) +
    theme_minimal() +
    labs(title = 'Cell types', colour = NULL, x = 't-SNE 1', y = 't-SNE 2')

# CNA scores:

cna_score_files <- dir(paste0('../data_and_figures/sc_find_malignant/', cohort))
cna_score_files <- cna_score_files[endsWith(cna_score_files, '_data.csv')]

cna_score_data <- lapply(
    cna_score_files,
    function(x) cbind(fread(paste0('../data_and_figures/sc_find_malignant/', cohort, '/', x)))#, patient = gsub('_data.csv', '', x))
) %>% rbindlist(fill = TRUE)

setkey(cna_score_data, cell_id)

cna_score_data <- cna_score_data[sc_data$id]

cna_score_data[
    ,
    classification_final := switch(
        (length(classification_final) == 1) + 1,
        switch( # In this case the cell must have been reference, but could have been classified as ambiguous or malignant in a subset of patients.
            ('ambiguous' %in% classification_final | 'malignant' %in% classification_final) + 1,
            'nonmalignant',
            'ambiguous'
        ),
        classification_final # In this case the cell was not reference.
    ),
    by = cell_id
]

cna_score_data <- cna_score_data[
    ,
    .(
        classification_final = unique(classification_final),
        cna_score = switch(
            startsWith(unique(classification_final), 'malignant') + 1,
            mean(cna_signal),
            .SD[classification == classification_final, cna_signal]
        )
    ),
    by = cell_id
]

cna_score_data[, patient := sc_data$patient]
cna_score_data[, cna_score := cna_score/quantile(cna_score[startsWith(classification_final, 'malignant')], 0.75), by = patient]

# cell_id should still be in the same order as sc_data$id, so it's safe to do this:
plot_data[, cna_score := cna_score_data$cna_score]

tsne_plot_cna_score <- ggplot(plot_data, aes(x = x, y = y, colour = cna_score)) +
    geom_point(size = 0.7) +
    theme_minimal() +
    scale_colour_gradientn(colours = colorRamps::matlab.like(50), limits = c(0, 1), oob = scales::squish) +
    labs(title = 'CNA score', colour = NULL, x = 't-SNE 1', y = 't-SNE 2')

# The following calculation of immune score is a bit complicated, but it's designed to put all the immune cells on roughly the same scale so that we
# don't get some immune cell types having much higher scores than the others.  It works OK but not brilliantly, possibly because taking the maximum
# makes it more noisy.

immune_cell_types <- c('B', 'B_plasma', 'DC', 'macrophage', 'mast', 'T')[
    sapply(
        c('B', 'B_plasma', 'DC', 'macrophage', 'mast', 'T'),
        function(ct) cell_type_markers[cell_type == ct, sum(gene %in% names(sc_data)) >= 3] # This was previously just > 0
    )
]

immune_cell_types <- immune_cell_types[
    mapvalues(
        immune_cell_types,
        c('B', 'B_plasma', 'DC', 'macrophage', 'mast', 'T'),
        c('b_cell', 'b_cell', 'dendritic', 'macrophage', 'mast', 't_cell'),
        warn_missing = FALSE
    ) %in% sc_data$cell_type
]

set.seed(cohort_data[[cohort]]$seeds[2])

immune_scores <- as.data.table(
    sapply(
        immune_cell_types,
        function(ct) {
            signature_score(
                sc_data[, set_colnames(t(.SD), id), .SDcols = -c('id', 'patient', 'cell_type')],
                sig_genes = cell_type_markers[cell_type == ct & gene %in% names(sc_data), unique(gene)],
                nbin = 30,
                n = 100
            )
        }
    ),
    keep.rownames = 'cell_id'
)

immune_sigs <- sapply(
    immune_cell_types,
    function(ct) {

        sig_cor <- cell_type_markers[
            cell_type == ct & gene %in% names(sc_data),
            cor(sc_data[, unique(gene), with = FALSE], immune_scores[, ..ct])[, 1]
        ]

        if(sum(sig_cor > 0.6) < 10) {
            if(length(sig_cor) <= 10) {
                return(names(sig_cor))
            } else {
                return(names(sig_cor)[names(sig_cor) %in% names(sort(sig_cor, decreasing = TRUE)[1:10])])
            }
        } else {
            return(names(sig_cor)[sig_cor > 0.6])
        }

    }
)

immune_scores <- as.data.table(
    sapply(
        immune_cell_types,
        function(ct) {
            signature_score(
                sc_data[, set_colnames(t(.SD), id), .SDcols = -c('id', 'patient', 'cell_type')],
                sig_genes = immune_sigs[[ct]],
                nbin = 30,
                n = 100
            )
        }
    ),
    keep.rownames = 'cell_id'
)

immune_scores[
    ,
    # which_max := apply(.SD, 1, function(x) names(.SD)[which.max(x)]),
    c('which_max', 'cell_type') := .(
        apply(.SD, 1, function(x) names(.SD)[which.max(x)]),
        sc_data$cell_type
    ),
    .SDcols = -'cell_id'
]

immune_cell_type_scaling_factors <- sapply(
    immune_cell_types,
    function(x) {
        immune_scores[
            cell_type == mapvalues(
                x,
                c('B', 'B_plasma', 'DC', 'macrophage', 'mast', 'T'),
                c('b_cell', 'b_cell', 'dendritic', 'macrophage', 'mast', 't_cell'),
                warn_missing = FALSE
            ),
            setNames(quantile(get(x), 0.9), NULL)
            # mean(get(x))
        ]
    },
    USE.NAMES = TRUE
)

immune_scores[
    ,
    score := get(unique(which_max))/immune_cell_type_scaling_factors[unique(which_max)],
    by = which_max
]

plot_data[
    ,
    c('immune_score', 'caf_score', 'endothelial_score') := .( # Should add pericyte score, but need to make my own signature
        immune_scores$score,
        signature_score(
            sc_data[, set_colnames(t(.SD), id), .SDcols = -c('id', 'patient', 'cell_type')],
            sig_genes = cell_type_markers[cell_type == 'CAF' & gene %in% names(sc_data), unique(gene)],
            nbin = 30,
            n = 100
        ),
        signature_score(
            sc_data[, set_colnames(t(.SD), id), .SDcols = -c('id', 'patient', 'cell_type')],
            sig_genes = cell_type_markers[cell_type == 'endothelial' & gene %in% names(sc_data), unique(gene)],
            nbin = 30,
            n = 100
        )
    )
]

caf_endothelial_sigs <- sapply(
    c('caf', 'endothelial'),
    function(ct) {

        sig_cor <- cell_type_markers[
            cell_type == mapvalues(ct, 'caf', 'CAF', warn_missing = FALSE) & gene %in% names(sc_data),
            cor(sc_data[, unique(gene), with = FALSE], plot_data[, get(paste0(ct, '_score'))])[, 1]
        ]

        if(sum(sig_cor > 0.6) < 10) {
            if(length(sig_cor) <= 10) {
                return(names(sig_cor))
            } else {
                return(names(sig_cor)[names(sig_cor) %in% names(sort(sig_cor, decreasing = TRUE)[1:10])])
            }
        } else {
            return(names(sig_cor)[sig_cor > 0.6])
        }

    }
)

plot_data[
    ,
    c('caf_score', 'endothelial_score') := .( # Should add pericyte score, but need to make my own signature
        signature_score(
            sc_data[, set_colnames(t(.SD), id), .SDcols = -c('id', 'patient', 'cell_type')],
            sig_genes = caf_endothelial_sigs$caf,
            nbin = 30,
            n = 100
        ),
        signature_score(
            sc_data[, set_colnames(t(.SD), id), .SDcols = -c('id', 'patient', 'cell_type')],
            sig_genes = caf_endothelial_sigs$endothelial,
            nbin = 30,
            n = 100
        )
    )
]

plot_data[
    ,
    c('caf_score', 'endothelial_score') := .(
        caf_score/.SD[cell_type == 'caf', quantile(caf_score, 0.9)],
        endothelial_score/.SD[cell_type == 'endothelial', quantile(endothelial_score, 0.9)]
        # caf_score/.SD[cell_type == 'caf', mean(caf_score)],
        # endothelial_score/.SD[cell_type == 'endothelial', mean(endothelial_score)]
    )
]

score_tsnes <- sapply(
    c('immune_score', 'caf_score', 'endothelial_score'),
    function(score_type) {
        ggplot(plot_data, aes(x = x, y = y, colour = get(score_type))) +
            geom_point(size = 0.7) +
            theme_minimal() +
            scale_colour_gradientn(colours = colorRamps::matlab.like(50), limits = c(0, 1), oob = scales::squish) +
            labs(
                title = mapvalues(
                    score_type,
                    c('immune_score', 'caf_score', 'endothelial_score'),
                    c('Immune score', 'CAF score', 'Endothelial score'),
                    warn_missing = FALSE
                ),
                colour = NULL,
                x = 't-SNE 1',
                y = 't-SNE 2'
            )
    },
    simplify = FALSE,
    USE.NAMES = TRUE
)





# Combined figure:

pdf('../data_and_figures/final_figures_resubmission/R3.pdf', width = 13, height = 12)

plot_grid(
    blank_plot() +
        theme(plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), 'pt'), plot.title = element_text(size = 16)) +
        labs(title = 'Lung Adenocarcinoma - Kim et al.'),
    plot_grid(
        plot_grid(
            plot_grid(
                plot_grid(
                    annotations_caf,
                    cna_heatmap_caf + theme(
                        axis.text = element_blank(),
                        axis.title = element_blank(),
                        legend.position = 'none',
                        plot.margin = unit(c(5.5, 5.5, 0, 0), 'pt')
                    ),
                    nrow = 1,
                    ncol = 2,
                    align = 'h',
                    rel_widths = c(1.5, 4)
                ),
                plot_grid(
                    annotations_cancer,
                    cna_heatmap_cancer + theme(
                        axis.text.y = element_blank(),
                        axis.title.y = element_blank(),
                        legend.position = 'none',
                        plot.margin = unit(c(40, 5.5, 5.5, 0), 'pt')
                    ),
                    nrow = 1,
                    ncol = 2,
                    align = 'h',
                    rel_widths = c(1.5, 4)
                ),
                nrow = 2,
                ncol = 1,
                rel_heights = c(1, 2)
                # rel_heights = c(1, 2)
            ),
            get_legend(cna_heatmap_caf),
            nrow = 1,
            ncol = 2,
            rel_widths = c(5, 1)
        ),
        plot_grid(
            tsne_plot_cell_types + theme(plot.margin = unit(c(10, 5.5, 10, 5.5), 'pt')),
            score_tsnes$caf_score + theme(plot.margin = unit(c(10, 5.5, 10, 5.5), 'pt')),
            tsne_plot_cna_score + theme(plot.margin = unit(c(10, 5.5, 10, 5.5), 'pt')),
            nrow = 3,
            ncol = 1,
            align = 'v'
        ),
        nrow = 1,
        ncol = 2,
        rel_widths = c(6.5, 3.5)
    ),
    nrow = 2,
    ncol = 1,
    rel_heights = c(1, 15)
)

dev.off()





expression_data <- fread('../../TCGA_data/tcga_expression_data.csv', key = 'id')
meta_data <- fread('../../TCGA_data/tcga_meta_data.csv', key = 'id')

expression_data <- expression_data[meta_data[sample_type != 'normal', id]]
meta_data <- meta_data[sample_type != 'normal']

emt_markers <- fread('../../emt_markers.csv')[, gene := alias2SymbolTable(gene)][source != 'GO', sort(unique(gene))]

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
        heatmap_annotations = c('ACTA2', 'BGN', 'COL1A2', 'LAMC2', 'QSOX1', 'THY1', 'TNFRSF12A', 'VEGFA'),
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

simulated_deconvs <- readRDS('../data_and_figures/simulated_deconvs.rds')

simulated_deconv_plots <- sapply(
    names(deconv_plot_args_per_ct),
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

lineplots <- readRDS('../data_and_figures/simulated_bulk_lineplots.rds')

# sample_sizes <- list(
#     brca = 1000,
# 	brca_lenient = 1000,
#     coadread = 800,
# 	coadread_lenient = 800,
#     hnsc = 400,
#     lihc = 150,
# 	lihc_lenient = 200,
#     luad = 1000,
# 	# luad = 200,
# 	# luad_lenient = 300,
# 	lusc = 100,
# 	lusc_lenient = 100,
#     ov = 1500,
# 	ov_lenient = 2000,
#     paad = 1500
# )

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

cairo_pdf('../data_and_figures/final_figures_resubmission/R4.pdf', width = 11, height = 12.25)

plot_grid(
    blank_plot() +
        labs(title = 'Breast cancer - Qian et al.') +
        theme(plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), 'pt'), plot.title = element_text(size = 16)),
    plot_grid(
        plotlist = c(
            lapply(
                c('brca', 'brca_lenient'),
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

                    plot_grid(
                        plot_grid(
                            blank_plot(),
                            plot_grid(get_y_axis(lineplots[[ct]]$lineplot)) +
                                draw_label('Proportion of gene expression', x = 0.4, y = 0.5, vjust = 1.3, angle = 90, size = 12),
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
                                    axis.title.y = element_text(vjust = -5),
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
                                    axis.title.y = element_text(vjust = -5),
                                    axis.ticks = element_blank(),
                                    axis.ticks.length = unit(0, 'pt'),
                                    plot.background = element_blank(),
                                    panel.background = element_blank(),
                                    plot.margin = unit(c(5.5, 0, 1.25, 5.5), 'pt')
                                ) +
                                labs(y = paste0('\nCAFs\n(n = ', length(sc_sim_deconv_comp[[ct]]$ordered_cell_ids$caf), ')')),
                            # blank_plot() +
                            #     labs(y = 'Cancer\ncells') +
                            #     scale_y_continuous(position = 'right') +
                            #     theme(plot.margin = unit(c(5.5, 5.5, 1.25, 5.5), 'pt'), axis.title.y.right = element_text(angle = 90)),
                            # blank_plot() +
                            #     labs(y = 'CAFs') +
                            #     scale_y_continuous(position = 'right') +
                            #     theme(plot.margin = unit(c(5.5, 5.5, 1.25, 5.5), 'pt'), axis.title.y.right = element_text(angle = 90)),
                            blank_plot(),
                            nrow = 7,
                            ncol = 1,
                            rel_heights = c(3, 15, 7, 15, 5, 5, 1)
                        ),
                        plot_grid(
                            plotlist = c(
                                list(
                                    blank_plot() +
                                        labs(title = switch(grepl('_lenient', ct) + 1, '"Strict" CAF definition', '"Lenient" CAF definition')) +
                                        theme(plot.title = element_text(hjust = 0.5, size = 16)),
                                    plot_grid(
                                        lineplots[[ct]]$lineplot + theme(
                                            legend.position = 'none',
                                            plot.margin = unit(c(0, 0, 0, 0), 'pt'),
                                            plot.title = element_blank(),
                                            axis.ticks.length = unit(0, 'pt'),
                                            axis.text = element_blank(),
                                            axis.title = element_blank()
                                        ),
                                        plot_grid(get_x_axis(lineplots[[ct]]$lineplot)) +
                                            draw_label('Proportion of cell mixture', x = 0.5, y = 0.55, vjust = -0.5, size = 12),
                                        nrow = 2,
                                        ncol = 1,
                                        rel_heights = c(15, 5)
                                    )
                                ),
                                sc_sim_deconv_comp_figures[c('purity_bar', 'ccle_bar', 'heatmap', 'cancer', 'caf')],
                                list(
                                    blank_plot() +
                                        scale_x_continuous(position = 'top') +
                                        theme(axis.title.x = element_text(), plot.margin = unit(c(4.5, 0, 0, 0), 'pt')) +
                                        labs(x = 'Genes')
                                )
                            ),
                            nrow = 8,
                            ncol = 1,
                            align = 'v',
                            rel_heights = c(3, 20, 1, 1, 15, 5, 5, 1)
                        ),
                        nrow = 1,
                        ncol = 2,
                        rel_widths = c(1.2, 4)
                    )

                }
            ),
            list(
                plot_grid(
                    blank_plot(),
                    dummy_legend,
                    blank_plot(),
                    get_legend(
                        simulated_deconv_plots$brca$plots$purity_bar +
                            guides(fill = guide_colourbar(title.position = 'right')) +
                            theme(
                                legend.justification = c(0, 1),
                                legend.direction = 'vertical',
                                legend.key.width = unit(10, 'pt'),
                                legend.key.height = unit(10, 'pt'),
                                legend.box.margin = margin(l = 10)
                            ) +
                            labs(fill = 'Correlation\nwith purity')
                    ),
                    get_legend(
                        simulated_deconv_plots$brca$plots$ccle_bar +
                            guides(fill = guide_colourbar(title.position = 'right')) +
                            theme(
                                legend.justification = c(0, 1),
                                legend.direction = 'vertical',
                                legend.key.width = unit(10, 'pt'),
                                legend.key.height = unit(10, 'pt'),
                                legend.box.margin = margin(l = 10)
                            ) +
                            labs(fill = 'Tumours vs.\ncell lines')
                    ),
                    get_legend(
                        simulated_deconv_plots$brca$plots$heatmap +
                            guides(fill = guide_colourbar(title.position = 'right')) +
                            theme(
                                legend.justification = c(0, 0),
                                legend.direction = 'vertical',
                                legend.key.width = NULL,
                                legend.key.height = NULL,
                                legend.box.margin = margin(l = 10)
                            ) +
                            labs(fill = 'Correlation')
                    ),
                    get_legend(
                        sc_sim_deconv_comp$brca$sc_heatmaps_unfiltered$expression_level_cc_rm[[1]] +
                            guides(fill = guide_colourbar(title.position = 'right')) +
                            theme(
                                legend.justification = c(0, 0),
                                legend.direction = 'vertical',
                                legend.key.width = NULL,
                                legend.key.height = NULL,
                                legend.box.margin = margin(l = 10)
                            ) +
                            labs(fill = 'Relative\nexpression\nlevel')
                    ),
                    blank_plot(),
                    nrow = 8,
                    ncol = 1,
                    rel_heights = c(3, 15, 5, 5, 5, 7, 10, 1)
                )
            )
        ),
        nrow = 1,
        ncol = 3,
        rel_widths = c(2, 2, 1)
    ),
    blank_plot(),
    nrow = 3,
    ncol = 1,
    rel_heights = c(1, 15, 0.5)
) %>% print

dev.off()





sc_metadata <- list(
    coadread = list(
		ref = 'lee_crc_2020_smc',
		tcga_cancer_types = c('COAD', 'READ'),
		seed = 9275,
		initial_cell_types = c('b_cell', 'endothelial', 'macrophage', 't_cell', 'caf', 'cancer'),
		rare_cell_types = c('epithelial', 'mast'),
		genes_filter_fun = function(x) {mean(x) >= 0.17 | sum(x >= 4.25) >= length(x)/100},
		scores_filter_fun = function(x) {mean(x) >= 0.46 | sum(x >= 4.25) >= length(x)/100},
		annotations = c('AREG', 'CALU', 'COL1A1', 'CXCL1', 'FN1', 'GADD45B', 'SNAI1', 'SNAI2', 'THY1', 'TWIST1', 'VIM', 'ZEB1', 'ZEB2'),
		annotations_title = 'Colorectal',
		annotations_side = 'left',
		max_mean_count = 1000,
		deconv_args = list(
			tcga_cancer_types = 'coadread',
			ccle_cancer_type = 'large intestine',
			seed = 6185,
			genes_filter_fun = function(x) 1:200,
			plot_title = 'Colorectal',
			heatmap_annotations = c('AREG', 'COL1A1', 'COL3A1', 'CXCL1', 'DCN', 'FN1', 'GADD45A', 'PLOD3', 'SNAI1', 'SNAI2', 'TWIST1', 'VIM',
                                    'ZEB1', 'ZEB2')
		),
		sample_size = 800,
		sim_deconv_filter = 0.17,
		# sim_deconv_filter = 0.1,
		pemt_bracket_max = 27,
		caf_bracket_min = 80
	)
)

all_figures <- list()

ct <- 'coadread'

sc_data <- readRDS(paste0('../data_and_figures/sc_alt_norm/', sc_metadata[[ct]]$ref, '_scran.rds'))
sc_meta <- fread(paste0('../data_and_figures/sc_alt_norm/', sc_metadata[[ct]]$ref, '_meta.csv'))

sc_data <- cbind(sc_meta, t(as.matrix(sc_data)))
rm(sc_meta)

if(paste0('genes_filtered_scran_', sc_metadata[[ct]]$ref, '.rds') %in% dir('../data_and_figures/sc_alt_norm')) {
    genes_filtered_scran <- readRDS(paste0('../data_and_figures/sc_alt_norm/genes_filtered_scran_', sc_metadata[[ct]]$ref, '.rds'))
} else {
    genes_unfiltered <- unique(
        c(
            emt_markers[emt_markers %in% names(sc_data)],
            top_cols_by_fun_cor(
                expression_data[meta_data[cancer_type %in% sc_metadata[[ct]]$tcga_cancer_types, id], -'id'],
                threshold_fun = function(x) quantile(x, 0.99)
            )[id %in% names(sc_data), id]
        )
    )
    genes_filtered_scran <- filter_for_groups(sc_data[, c('cell_type', ..genes_unfiltered)], groups = c('caf', 'cancer'))
    saveRDS(genes_filtered_scran, paste0('../data_and_figures/sc_alt_norm/genes_filtered_scran_', sc_metadata[[ct]]$ref, '.rds'))
}

set.seed(sc_metadata[[ct]]$seed)

if(paste0('sc_cancer_caf_scran_', sc_metadata[[ct]]$ref, '.rds') %in% dir('../data_and_figures/sc_alt_norm')) {
    sc_cancer_caf_scran <- readRDS(paste0('../data_and_figures/sc_alt_norm/sc_cancer_caf_scran_', sc_metadata[[ct]]$ref, '.rds'))
} else {
    sc_cancer_caf_scran <- sc_groups(
        genes = genes_filtered_scran,
        sc_data = sc_data[cell_type %in% c('cancer', 'caf')],
        groups = c('cancer', 'caf'),
        score_cells_nbin = 30,
        score_cells_n = 40,
        min_sig_size = 0,
        scores_filter_groups = 'cancer',
        genes_filter_fun = sc_metadata[[ct]]$genes_filter_fun,
        scores_filter_fun = sc_metadata[[ct]]$scores_filter_fun
    )
    saveRDS(sc_cancer_caf_scran, paste0('../data_and_figures/sc_alt_norm/sc_cancer_caf_scran_', sc_metadata[[ct]]$ref, '.rds'))
}

sc_cancer_caf_heatmap_scran <- sc_groups_heatmap(
    sc_groups_list = sc_cancer_caf_scran,
    groups = c('cancer', 'caf'),
    x_axis_titles = c('Cancer cells', 'CAFs'),
    default_figure_widths = list(annotations = 2.2, cancer = 6, caf = 1.5),
    figure_spacing = 2.5,
    annotations_nudge = 0.25,
    es_fun = NULL,
    es_title = 'EMT\nscore',
    h_limits = c(0, 8),
    h_legend_breaks = c(0, 2, 4, 6, 8),
    h_legend_labels = c('0' = '0', '2' = '2', '4' = '4', '6' = '6', '8' = '\u2265 8'),
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
    annotations = sc_metadata[[ct]]$annotations,
    title = sc_metadata[[ct]]$annotations_title
)

if('annotations_side' %in% names(sc_metadata[[ct]])) {
    sc_cancer_caf_heatmap_combining <- sc_groups_heatmap(
        sc_groups_list = sc_cancer_caf_scran,
        groups = c('cancer', 'caf'),
        x_axis_titles = c('Cancer cells', 'CAFs'),
        default_figure_widths = list(annotations = 2.5, cancer = 6, caf = 1.2),
        figure_spacing = 2.5,
        annotations_nudge = 0.25,
        es_fun = NULL,
        es_title = 'EMT score',
        h_limits = c(0, 8),
        h_legend_breaks = c(0, 2, 4, 6, 8),
        h_legend_labels = c('0' = '0', '2' = '2', '4' = '4', '6' = '6', '8' = '\u2265 8'),
        annotations = sc_metadata[[ct]]$annotations,
        annotations_title = sc_metadata[[ct]]$annotations_title,
        annotations_side = sc_metadata[[ct]]$annotations_side
    )
}

sc_data[, cell_type := mapvalues(cell_type, sc_metadata[[ct]]$rare_cell_types, rep('rare', length(sc_metadata[[ct]]$rare_cell_types)))]

if(paste0('lineplot_scran_', sc_metadata[[ct]]$ref, '.rds') %in% dir('../data_and_figures/sc_alt_norm')) {
    lineplot_scran <- readRDS(paste0('../data_and_figures/sc_alt_norm/lineplot_scran_', sc_metadata[[ct]]$ref, '.rds'))
} else {
    lineplot_scran <- simulated_tumours_lineplot(
        sc_data,
        genes_filtered_scran,
        initial_types = sc_metadata[[ct]]$initial_cell_types,
        normalise_fun = NULL,
        x_axis_title = 'Proportion of cell mixture',
        plot_title = sc_metadata[[ct]]$annotations_title,
        max_mean_count = sc_metadata[[ct]]$max_mean_count,
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
    saveRDS(lineplot_scran, paste0('../data_and_figures/sc_alt_norm/lineplot_scran_', sc_metadata[[ct]]$ref, '.rds'))
}

if(paste0('simulated_bulk_data_', sc_metadata[[ct]]$ref, '.csv') %in% dir('../data_and_figures/sc_alt_norm')) {
    simulated_bulk_data <- fread(paste0('../data_and_figures/sc_alt_norm/simulated_bulk_data_', sc_metadata[[ct]]$ref, '.csv'), key = 'id')
    simulated_bulk_metadata <- fread(paste0('../data_and_figures/sc_alt_norm/simulated_bulk_metadata_', sc_metadata[[ct]]$ref, '.csv'), key = 'id')
} else {
    simulated_bulk_data <- simulated_tumours_data(
        as.matrix(sc_data[, -c('id', 'patient', 'cell_type')]),
        types = sc_data$cell_type,
        id_prefix = ct,
        max_mean_count = sc_data[cell_type == 'cancer', .N, by = patient][, round(quantile(N, 0.9))]
        # max_mean_count = max_mean_counts[[ct]]
        # genes = genes_list[[ct]]
    )
    simulated_bulk_metadata <- as.data.table(simulated_bulk_data$meta_data, keep.rownames = 'id')
    simulated_bulk_data <- as.data.table(simulated_bulk_data$expression_data, keep.rownames = 'id')
    simulated_bulk_metadata[, cancer_type := gsub('[0-9]+', '', id)]
    setcolorder(simulated_bulk_metadata, c('id', 'cancer_type'))
    setkey(simulated_bulk_data, 'id')
    setkey(simulated_bulk_metadata, 'id')
    fwrite(simulated_bulk_data, paste0('../data_and_figures/sc_alt_norm/simulated_bulk_data_', sc_metadata[[ct]]$ref, '.csv'))
    fwrite(simulated_bulk_metadata, paste0('../data_and_figures/sc_alt_norm/simulated_bulk_metadata_', sc_metadata[[ct]]$ref, '.csv'))
}

ccle <- fread('../../CCLE_cancer_type_Av.csv', key = 'gene_id')
cell_type_markers <- fread('../../cell_type_markers.csv')

if(paste0('simulated_deconv_', sc_metadata[[ct]]$ref, '.rds') %in% dir('../data_and_figures/sc_alt_norm')) {
    simulated_deconv <- readRDS(paste0('../data_and_figures/sc_alt_norm/simulated_deconv_', sc_metadata[[ct]]$ref, '.rds'))
} else {

    simulated_deconv <- do.call(
        deconvolve_emt_caf_data,
        args = c(
            list(
                expression_data = simulated_bulk_data,
                meta_data = simulated_bulk_metadata,
                genes = emt_markers,
                cell_type_markers = cell_type_markers,
                ccle_data = ccle,
                initial_gene_weights = FALSE
            ),
            sc_metadata[[ct]]$deconv_args[!(names(sc_metadata[[ct]]$deconv_args) %in% c('plot_title', 'heatmap_annotations'))]
        )
    )

    # Reorder genes in deconv results:
    simulated_deconv <- deconv_reorder(simulated_deconv)

    # Add data for single cell comparison colour bar:

    sc_data_melted <- melt(
        sc_data[, c('id', 'cell_type', simulated_deconv$genes_filtered), with = FALSE],
        id.vars = c('id', 'cell_type'),
        variable.name = 'gene',
        value.name = 'expression_level'
    )

    # Centre genes:
    gene_averages <- sc_data_melted[
        ,
        .(ave_exp = mean(c(mean(expression_level[cell_type == 'cancer']), mean(expression_level[cell_type == 'caf'])))),
        by = .(symbol = gene)
    ]
    sc_data_melted[, expression_level := expression_level - gene_averages[symbol == unique(gene), ave_exp], by = gene]

    # Centre cells:
    sc_data_melted[, expression_level := expression_level - mean(expression_level), by = id]

    simulated_deconv$extra_data_score <- sc_data_melted[
        ,
        .(d = mean(expression_level[cell_type == 'cancer']) - mean(expression_level[cell_type == 'caf'])),
        by = gene
    ][, setNames(d, gene)]

    saveRDS(simulated_deconv, paste0('../data_and_figures/sc_alt_norm/simulated_deconv_', sc_metadata[[ct]]$ref, '.rds'))

}

simulated_deconv_plots_scran <- do.call(
    deconvolve_emt_caf_plots,
    args = c(
        list(
            data = simulated_deconv,
            heatmap_legend_title = 'Correlation',
            heatmap_colours = rev(colorRampPalette(brewer.pal(11, "RdBu"))(50)),
            heatmap_colour_limits = c(-1, 1),
            heatmap_legend_breaks = c(-1, -0.5, 0, 0.5, 1),
            heatmap_legend_justification = 'left',
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
            # extra_colours = c(
                # colorspace::sequential_hcl(25, h = 190, c = 70, l = c(70, 99), power = 0.8),
                # colorspace::sequential_hcl(25, h = 60, c = 100, l = c(80, 99), power = 0.8, rev = TRUE)
            # ),
            extra_colours = colorRampPalette(c('turquoise4', 'turquoise', 'azure', 'gold2', 'gold4'))(50),
            extra_fun = function(x) caTools::runmean(x, 30),
            extra_colour_limits = c(-4, 4),
            extra_legend_breaks = c(-4, 0, 4),
            extra_legend_labels = c('-4' = '\u2264 -4', '0' = '0', '4' = '\u2265 4'),
            extra_legend_title = 'scRNA-seq: CAF vs. cancer\n',
            extra_legend_direction = 'horizontal',
            extra_axis_title = NULL,
            bar_legend_justification = 'left'
            # bar_legend_width = NULL,
            # bar_legend_height = NULL
        ),
        sc_metadata[[ct]]$deconv_args[names(sc_metadata[[ct]]$deconv_args) %in% c('plot_title', 'heatmap_annotations')]
    )
)

sc_deconv_comp <- sapply(
    c('cancer', 'caf'),
    function(x) {

        # plot_data <- copy(sc_data[cell_type == x])[
        #     ,
        #     complexity := apply(.SD, 1, function(x) sum(x > 0)),
        #     .SDcols = -c('id', 'patient', 'cell_type')
        # ][
        #     ,
        #     .SD[sample(1:.N, sc_metadata[[ct]]$sample_size, prob = complexity)]
        # ][, complexity := NULL][, c('id', simulated_deconv$genes_filtered), with = FALSE]

        plot_data <- sc_data[cell_type == x, c('id', simulated_deconv$genes_filtered), with = FALSE]

        plot_data <- melt(plot_data, id.vars = 'id', variable.name = 'gene', value.name = 'expression_level')

        ordered_cell_ids <- plot_data[
            gene %in% with(
                simulated_deconv,
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
        runmean(setNames(expression_level_cc, gene)[with(simulated_deconv, genes_filtered[ordering])], 30),
        with(simulated_deconv, genes_filtered[ordering])
    )[as.character(gene)], # Melting made gene a factor, and this reordering doesn't work unless we coerce it to character
    by = id
]

sc_heatmaps_unfiltered <- sapply(
    c('expression_level', 'expression_level_cc', 'expression_level_cc_rm'),
    function(expr_var) sapply(
        c('cancer', 'caf'),
        function(x) {
            ggplot(
                sc_deconv_comp_data[cell_type == x],
                aes(
                    x = factor(gene, levels = with(simulated_deconv, genes_filtered[ordering])),
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
                    limits = switch((expr_var == 'expression_level_cc_rm') + 1, c(-3, 3), c(-1.5, 1.5)),
                    breaks = switch((expr_var == 'expression_level_cc_rm') + 1, c(-3, -1.5, 0, 1.5, 3), c(-1.5, 0, 1.5)),
                    labels = switch(
                        (expr_var == 'expression_level_cc_rm') + 1,
                        c('-3' = '\u2264 -3', '-1.5' = '-1.5', '0' = '0', '1.5' = '1.5', '3' = '\u2265 3'),
                        c('-1.5' = '\u2264 -1.5', '0' = '0', '1.5' = '\u2265 1.5')
                    ),
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

simulated_deconv_filtered <- simulated_deconv

filtered_genes <- gene_averages_cancer_caf[
    ,
    .(pass = ave_exp[cell_type == 'cancer'] > sc_metadata[[ct]]$sim_deconv_filter |
        ave_exp[cell_type == 'caf'] > sc_metadata[[ct]]$sim_deconv_filter),
    by = gene
][pass == TRUE, as.character(gene)]

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
        sc_metadata[[ct]]$deconv_args[names(sc_metadata[[ct]]$deconv_args) %in% c('plot_title', 'heatmap_annotations')]
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
            ggplot(plot_data[[x]], aes(x = gene, y = id, fill = get(expr_var))) +
                geom_raster() +
                scale_x_discrete(expand = c(0, 0)) +
                scale_y_discrete(expand = c(0, 0)) +
                scale_fill_gradientn(
                    colours = c(
                        sequential_hcl(25, h = 190, c = 70, l = c(70, 99), power = 0.8),
                        sequential_hcl(25, h = 60, c = 100, l = c(80, 99), power = 0.8, rev = TRUE)
                    ),
                    limits = switch((expr_var == 'expression_level_cc_rm') + 1, c(-3, 3), c(-1.5, 1.5)),
                    breaks = switch((expr_var == 'expression_level_cc_rm') + 1, c(-3, -1.5, 0, 1.5, 3), c(-1.5, 0, 1.5)),
                    labels = switch(
                        (expr_var == 'expression_level_cc_rm') + 1,
                        c('-3' = '\u2264 -3', '-1.5' = '-1.5', '0' = '0', '1.5' = '1.5', '3' = '\u2265 3'),
                        c('-1.5' = '\u2264 -1.5', '0' = '0', '1.5' = '\u2265 1.5')
                    ),
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

all_figures[[ct]] <- list(
    sc_cancer_caf_heatmap_scran = sc_cancer_caf_heatmap_scran,
    lineplot_scran = lineplot_scran$lineplot,
    simulated_deconv_plots = simulated_deconv_plots_scran$plots,
    simulated_deconv_filtered_plots = simulated_deconv_filtered_plots$plots,
    sc_heatmaps_unfiltered = sc_heatmaps_unfiltered,
    sc_heatmaps_filtered = sc_heatmaps_filtered
)

if('annotations_side' %in% names(sc_metadata[[ct]])) {
    all_figures[[ct]] <- c(all_figures[[ct]], list(sc_cancer_caf_heatmap_combining = sc_cancer_caf_heatmap_combining))
}

sc_sim_deconv_comp_figures_tpm <- sapply(
    c(simulated_deconv_plots$coadread$plots, sc_sim_deconv_comp$coadread$sc_heatmaps_unfiltered$expression_level_cc_rm),
    # c(sc_sim_deconv_comp$coadread$filtered_deconv_figures, sc_sim_deconv_comp$coadread$sc_heatmaps_filtered),
    function(x) x + theme(legend.position = 'none', plot.title = element_blank(), axis.title.y = element_blank()),
    simplify = FALSE,
    USE.NAMES = TRUE
)
sc_sim_deconv_comp_figures_tpm$heatmap <- sc_sim_deconv_comp_figures_tpm$heatmap + theme(plot.margin = unit(c(0, 5.5, 1.25, 0), 'pt'))
sc_sim_deconv_comp_figures_tpm$axis_labels <- sc_sim_deconv_comp_figures_tpm$axis_labels + theme(plot.margin = unit(c(0, 0, 1.25, 0), 'pt'))

sc_sim_deconv_comp_figures_scran <- sapply(
    c(all_figures$coadread$simulated_deconv_plots, all_figures$coadread$sc_heatmaps_unfiltered$expression_level_cc_rm),
    function(x) x + theme(legend.position = 'none', plot.title = element_blank(), axis.title.y = element_blank()),
    simplify = FALSE,
    USE.NAMES = TRUE
)
sc_sim_deconv_comp_figures_scran$heatmap <- sc_sim_deconv_comp_figures_scran$heatmap + theme(plot.margin = unit(c(0, 5.5, 1.25, 0), 'pt'))
sc_sim_deconv_comp_figures_scran$axis_labels <- sc_sim_deconv_comp_figures_scran$axis_labels + theme(plot.margin = unit(c(0, 0, 1.25, 0), 'pt'))

dummy_legend_plot <- ggplot(
	data = data.table(
		x = 1:6,
		y = 1,
		f = factor(
			c('b_cell', 'cancer', 'endothelial', 'caf', 'macrophage', 't_cell'),
			levels = c('b_cell', 'cancer', 'endothelial', 'caf', 'macrophage', 't_cell')
		)
	)
) +
	geom_tile(aes(x = x, y = y, fill = f)) +
	scale_fill_manual(
		labels = c(
			'b_cell' = 'B cell',
			'cancer' = 'Cancer',
			'endothelial' = 'Endothelial',
			'caf' = 'CAF',
			'macrophage' = 'Macrophage',
			't_cell' = 'T cell'
		),
		values = c(
			'b_cell' = '#8DD3C7',
			'cancer' = '#FB8072',
			'endothelial' = '#BC80BD',
			'caf' = '#FDB462',
			'macrophage' = '#80B1D3',
			't_cell' = '#B3DE69'
		)
	) +
	labs(fill = 'Cell type') +
	theme(legend.key = element_rect(size = 1, colour = 'white'), legend.key.size = unit(15, 'pt'))

dummy_legend <- get_legend(dummy_legend_plot)

cairo_pdf('../data_and_figures/final_figures_resubmission/R2.pdf', width = 11, height = 12)

plot_grid(
    blank_plot() +
        labs(title = 'Colorectal cancer - Lee et al. (SMC cohort)') +
        theme(plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), 'pt'), plot.title = element_text(size = 16)),
    plot_grid(
        plot_grid(
            plot_grid(
                blank_plot(),
                plot_grid(get_y_axis(lineplots$coadread$lineplot)) +
                    draw_label('Proportion of gene expression', x = 0.4, y = 0.5, vjust = 1.3, angle = 90, size = 12),
                blank_plot(),
                sc_sim_deconv_comp_figures_tpm$axis_labels,
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
                        axis.title.y = element_text(vjust = -5),
                        axis.ticks = element_blank(),
                        axis.ticks.length = unit(0, 'pt'),
                        plot.background = element_blank(),
                        panel.background = element_blank(),
                        plot.margin = unit(c(5.5, 0, 1.25, 5.5), 'pt')
                    ) +
                    labs(y = paste0('Cancer\ncells\n(n = ', length(sc_sim_deconv_comp$coadread$ordered_cell_ids$cancer), ')')),
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
                        axis.title.y = element_text(vjust = -5),
                        axis.ticks = element_blank(),
                        axis.ticks.length = unit(0, 'pt'),
                        plot.background = element_blank(),
                        panel.background = element_blank(),
                        plot.margin = unit(c(5.5, 0, 1.25, 5.5), 'pt')
                    ) +
                    labs(y = paste0('\nCAFs\n(n = ', length(sc_sim_deconv_comp$coadread$ordered_cell_ids$caf), ')')),
                # blank_plot() +
                #     labs(y = 'Cancer\ncells') +
                #     scale_y_continuous(position = 'right') +
                #     theme(plot.margin = unit(c(5.5, 5.5, 1.25, 5.5), 'pt'), axis.title.y.right = element_text(angle = 90)),
                # blank_plot() +
                #     labs(y = 'CAFs') +
                #     scale_y_continuous(position = 'right') +
                #     theme(plot.margin = unit(c(5.5, 5.5, 1.25, 5.5), 'pt'), axis.title.y.right = element_text(angle = 90)),
                blank_plot(),
                nrow = 7,
                ncol = 1,
                rel_heights = c(3, 15, 7, 15, 5, 5, 3)
            ),
            plot_grid(
                plotlist = c(
                    list(
                        blank_plot() +
                            labs(title = 'Normalisation by TPM/10') +
                            theme(plot.title = element_text(hjust = 0.5, size = 16)),
                        plot_grid(
                            lineplots$coadread$lineplot + theme(
                                legend.position = 'none',
                                plot.margin = unit(c(0, 0, 0, 0), 'pt'),
                                plot.title = element_blank(),
                                axis.ticks.length = unit(0, 'pt'),
                                axis.text = element_blank(),
                                axis.title = element_blank()
                            ),
                            plot_grid(get_x_axis(lineplots$coadread$lineplot)) +
                                draw_label('Proportion of cell mixture', x = 0.5, y = 0.55, vjust = -0.5, size = 12),
                            nrow = 2,
                            ncol = 1,
                            rel_heights = c(15, 5)
                        )
                    ),
                    sc_sim_deconv_comp_figures_tpm[c('purity_bar', 'ccle_bar', 'heatmap', 'cancer', 'caf')],
                    list(
                        plot_grid(
                            get_legend(
                                sc_sim_deconv_comp$coadread$sc_heatmaps_unfiltered$expression_level_cc_rm[[1]] +
                                    guides(fill = guide_colourbar(title.position = 'right')) +
                                    theme(
                                        legend.justification = c(0, 1),
                                        legend.direction = 'horizontal',
                                        legend.key.width = unit(25, 'pt'),
            							legend.key.height = unit(10, 'pt'),
                                        legend.title = element_text(margin = margin(l = 7.5))
                                    ) +
                                    labs(fill = 'Relative expression\n level')
                            )
                        )
                    )
                ),
                nrow = 8,
                ncol = 1,
                align = 'v',
                rel_heights = c(3, 20, 1, 1, 15, 5, 5, 3)
            ),
            nrow = 1,
            ncol = 2,
            rel_widths = c(1.2, 4)
        ),
        plot_grid(
            plot_grid(
                blank_plot(),
                plot_grid(get_y_axis(all_figures$coadread$lineplot_scran)) +
                    draw_label('Proportion of gene expression', x = 0.4, y = 0.5, vjust = 1.3, angle = 90, size = 12),
                blank_plot(),
                sc_sim_deconv_comp_figures_scran$axis_labels,
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
                        axis.title.y = element_text(vjust = -5),
                        axis.ticks = element_blank(),
                        axis.ticks.length = unit(0, 'pt'),
                        plot.background = element_blank(),
                        panel.background = element_blank(),
                        plot.margin = unit(c(5.5, 0, 1.25, 5.5), 'pt')
                    ) +
                    labs(y = paste0('Cancer\ncells\n(n = ', length(sc_deconv_comp$cancer$ordered_cell_ids), ')')),
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
                        axis.title.y = element_text(vjust = -5),
                        axis.ticks = element_blank(),
                        axis.ticks.length = unit(0, 'pt'),
                        plot.background = element_blank(),
                        panel.background = element_blank(),
                        plot.margin = unit(c(5.5, 0, 1.25, 5.5), 'pt')
                    ) +
                    labs(y = paste0('\nCAFs\n(n = ', length(sc_deconv_comp$caf$ordered_cell_ids), ')')),
                # blank_plot() +
                #     labs(y = 'Cancer\ncells') +
                #     scale_y_continuous(position = 'right') +
                #     theme(plot.margin = unit(c(5.5, 5.5, 1.25, 5.5), 'pt'), axis.title.y.right = element_text(angle = 90)),
                # blank_plot() +
                #     labs(y = 'CAFs') +
                #     scale_y_continuous(position = 'right') +
                #     theme(plot.margin = unit(c(5.5, 5.5, 1.25, 5.5), 'pt'), axis.title.y.right = element_text(angle = 90)),
                blank_plot(),
                nrow = 7,
                ncol = 1,
                rel_heights = c(3, 15, 7, 15, 5, 5, 3)
            ),
            plot_grid(
                plotlist = c(
                    list(
                        blank_plot() +
                            labs(title = 'Normalisation by scran') +
                            theme(plot.title = element_text(hjust = 0.5, size = 16)),
                        plot_grid(
                            all_figures$coadread$lineplot_scran + theme(
                                legend.position = 'none',
                                plot.margin = unit(c(0, 0, 0, 0), 'pt'),
                                plot.title = element_blank(),
                                axis.ticks.length = unit(0, 'pt'),
                                axis.text = element_blank(),
                                axis.title = element_blank()
                            ),
                            plot_grid(get_x_axis(all_figures$coadread$lineplot_scran)) +
                                draw_label('Proportion of cell mixture', x = 0.5, y = 0.55, vjust = -0.5, size = 12),
                            nrow = 2,
                            ncol = 1,
                            rel_heights = c(15, 5)
                        )
                    ),
                    sc_sim_deconv_comp_figures_scran[c('purity_bar', 'ccle_bar', 'heatmap', 'cancer', 'caf')],
                    list(
                        plot_grid(
                            get_legend(
                                all_figures$coadread$sc_heatmaps_unfiltered$expression_level_cc_rm$cancer +
                                    guides(fill = guide_colourbar(title.position = 'right')) +
                                    theme(
                                        legend.justification = c(0, 1),
                                        legend.direction = 'horizontal',
                                        legend.key.width = unit(20, 'pt'),
            							legend.key.height = unit(10, 'pt'),
                                        legend.title = element_text(margin = margin(l = 7.5))
                                    ) +
                                    labs(fill = ' Relative expression\n  level')
                            )
                        )
                    )
                ),
                nrow = 8,
                ncol = 1,
                align = 'v',
                rel_heights = c(3, 20, 1, 1, 15, 5, 5, 3)
            ),
            nrow = 1,
            ncol = 2,
            rel_widths = c(1.2, 4)
        ),
        plot_grid(
            blank_plot(),
            dummy_legend,
            blank_plot(),
            get_legend(
                sc_sim_deconv_comp$coadread$filtered_deconv_figures$purity_bar +
                    guides(fill = guide_colourbar(title.position = 'right')) +
                    theme(
                        legend.justification = c(0, 1),
                        legend.direction = 'vertical',
                        legend.key.width = unit(10, 'pt'),
                        legend.key.height = unit(10, 'pt'),
                        legend.box.margin = margin(l = 10)
                    ) +
                    labs(fill = 'Correlation\nwith purity')
            ),
            get_legend(
                sc_sim_deconv_comp$coadread$filtered_deconv_figures$ccle_bar +
                    guides(fill = guide_colourbar(title.position = 'right')) +
                    theme(
                        legend.justification = c(0, 1),
                        legend.direction = 'vertical',
                        legend.key.width = unit(10, 'pt'),
                        legend.key.height = unit(10, 'pt'),
                        legend.box.margin = margin(l = 10)
                    ) +
                    labs(fill = 'Tumours vs.\ncell lines')
            ),
            get_legend(
                sc_sim_deconv_comp$coadread$filtered_deconv_figures$heatmap +
                    guides(fill = guide_colourbar(title.position = 'right')) +
                    theme(
                        legend.justification = c(0, 0),
                        legend.direction = 'vertical',
                        legend.key.width = NULL,
                        legend.key.height = NULL,
                        legend.box.margin = margin(l = 10)
                    ) +
                    labs(fill = 'Correlation')
            ),
            # get_legend(
            #     sc_sim_deconv_comp$coadread$sc_heatmaps_filtered[[1]] +
            #         guides(fill = guide_colourbar(title.position = 'right')) +
            #         theme(
            #             legend.justification = c(0, 0),
            #             legend.direction = 'vertical',
            #             legend.key.width = NULL,
            #             legend.key.height = NULL,
            #             legend.box.margin = margin(l = 10)
            #         ) +
            #         labs(fill = 'Relative\nexpression\nlevel')
            # ),
            blank_plot(),
            nrow = 7,
            ncol = 1,
            rel_heights = c(3, 15, 5, 5, 5, 7, 13)
        ),
        nrow = 1,
        ncol = 3,
        rel_widths = c(2, 2, 1)
    ),
    nrow = 2,
    ncol = 1,
    rel_heights = c(1, 15)
) %>% print

dev.off()





# Correlations between pEMT and CAF signatures:

expression_data <- fread('../../TCGA_data/tcga_expression_data.csv', key = 'id')
meta_data <- fread('../../TCGA_data/tcga_meta_data.csv', key = 'id')

expression_data <- expression_data[meta_data[sample_type != 'normal', id]]
meta_data <- meta_data[sample_type != 'normal']

deconv_data <- readRDS('../data_and_figures/deconv_data.rds')
deconv_plots <- readRDS('../data_and_figures/deconv_plots.rds')

ct_to_keep <- c('blca_luminal_papillary', 'blca_basal_squamous', 'brca_luminal_a', 'brca_luminal_b', 'brca_basal_like', 'brca_her2_enriched',
    'coad', 'esca_ac', 'hnsc_mesenchymal_basal', 'hnsc_classical', 'luad_proximal_inflammatory', 'luad_proximal_proliferative',
    'luad_terminal_respiratory_unit', 'lusc_classical', 'lusc_secretory', 'ov_differentiated', 'ov_immunoreactive', 'ov_proliferative', 'paad',
    'read', 'stad_cin', 'stad_msi', 'ucec')
nice_names_for_figure <- c('BLCA - Luminal-Papillary', 'BLCA - Basal-Squamous', 'BRCA - Luminal A', 'BRCA - Luminal B', 'BRCA - Basal-like',
    'BRCA - HER2-enriched', 'COAD', 'ESCA - Adenocarcinoma', 'HNSC - Malignant-Basal', 'HNSC - Classical', 'LUAD - Squamoid', 'LUAD - Magnoid',
    'LUAD - Bronchioid', 'LUSC - Classical', 'LUSC - Secretory', 'OV - Differentiated', 'OV - Immunoreactive', 'OV - Proliferative', 'PAAD', 'READ',
    'STAD - CIN', 'STAD - MSI', 'UCEC')

scores_data <- sapply(
    ct_to_keep,
    function(ct) {

        cat(paste0(ct, ':\n'))

        sample_ids_mat <- expression_data[deconv_data[[ct]]$sample_ids, set_colnames(t(.SD), id), .SDcols = -'id']
        max_n <- floor(length(deconv_data[[ct]]$genes_filtered)/3)

        cat('\tComputing mes score\n')
        mes_score <- signature_score(sample_ids_mat, deconv_data[[ct]]$genes_filtered, nbin = 20, n = 100)

        cat('\tComputing scores for random samples\n')
        gene_samples <- replicate(50, sample(deconv_data[[ct]]$genes_filtered, 20, replace = FALSE), simplify = FALSE)
        gene_sample_scores <- lapply(gene_samples, function(smpl) signature_score(sample_ids_mat, smpl, nbin = 20, n = 100))
        ave_mes_corr <- mean(sapply(gene_sample_scores, function(smpl_scores) cor(smpl_scores, mes_score)))

        cat('\tComputing pEMT and CAF scores: ')
        pemt_caf_scores <- lapply(
            20:max_n,
            function(i) {
                if(i == 20) {
                    cat('1/', max_n - 19, sep = '')
                } else {
                    cat(rep('\b', ceiling(log10(i - 19)) + 1 + ceiling(log10(max_n - 19))), i - 19, '/', max_n - 19, sep = '')
                }
                caf_scores <- with(deconv_data[[ct]], signature_score(sample_ids_mat, tail(genes_filtered[ordering], i), nbin = 20, n = 100))
                pemt_scores <- with(deconv_data[[ct]], signature_score(sample_ids_mat, head(genes_filtered[ordering], i), nbin = 20, n = 100))
                if(i == floor(length(deconv_data[[ct]]$genes_filtered)/3)) cat('\n')
                list(caf_scores = caf_scores, pemt_scores = pemt_scores)
            }
        )

        cat('\tCompiling data and fitting model\n')
        pemt_caf_scores_corr <- data.table(
            ngenes = 20:(length(pemt_caf_scores) + 19),
            score_corr = sapply(pemt_caf_scores, function(x) cor(x$pemt_scores, x$caf_scores))
        )
        pemt_caf_scores_mod <- loess(score_corr ~ ngenes, pemt_caf_scores_corr, span = 0.5, degree = 1)
        pemt_caf_scores_corr[, loess_fit := predict(pemt_caf_scores_mod, .SD)]

        list(
            mes_score = mes_score,
            gene_samples = gene_samples,
            gene_sample_scores = gene_sample_scores,
            ave_mes_corr = ave_mes_corr,
            pemt_caf_scores = pemt_caf_scores,
            pemt_caf_scores_corr = pemt_caf_scores_corr,
            pemt_caf_scores_mod = pemt_caf_scores_mod
        )

    },
    simplify = FALSE,
    USE.NAMES = TRUE
)

barplot_data <- rbindlist(
    lapply(
        ct_to_keep,
        function(ct) data.table(
            cancer_type = ct,
            Strict = scores_data[[ct]]$pemt_caf_scores_corr[ngenes == 20, loess_fit],
            Lenient = scores_data[[ct]]$pemt_caf_scores_corr[ngenes == max(ngenes), loess_fit]
        )
    )
) %>% melt(id.vars = 'cancer_type', measure.vars = c('Strict', 'Lenient'), variable.name = 'Signature threshold', value.name = 'corr')

barplot_data[
    ,
    c('cancer_type', 'Signature threshold') := .(
        mapvalues(
            factor(
                cancer_type,
                levels = .SD[`Signature threshold` == 'Lenient', cancer_type[order(-corr)]]
                # levels = .SD[, .(corr_diff = abs(diff(corr))), by = cancer_type][order(-corr_diff), cancer_type]
            ),
            ct_to_keep,
            nice_names_for_figure
        ),
        factor(`Signature threshold`, levels = c('Lenient', 'Strict'))
    )
]

summary_barplot <- ggplot(
    data = barplot_data,
    aes(x = cancer_type, y = corr, group = `Signature threshold`, fill = `Signature threshold`, colour = `Signature threshold`)
) +
    geom_hline(yintercept = 0, size = 0.5, colour = 'grey') +
    geom_col(position = 'dodge', width = 0.7) +
    theme_test() +
    theme(axis.text.x = element_text(angle = 55, hjust = 1)) +
    labs(x = NULL, y = 'Correlation', fill = NULL, colour = NULL)

scatterplot_data <- rbindlist(
    lapply(
        c('blca_basal_squamous', 'brca_luminal_b', 'luad_proximal_inflammatory'),
        function(ct) cbind(cancer_type = ct, scores_data[[ct]]$pemt_caf_scores_corr)
    )
)

scatterplots <- ggplot(data = scatterplot_data) +
    geom_point(aes(x = ngenes, y = score_corr), alpha = 0.5) +
    geom_path(aes(x = ngenes, y = loess_fit), colour = 'royalblue') +
    facet_wrap(~mapvalues(cancer_type, ct_to_keep, nice_names_for_figure, warn_missing = FALSE), nrow = 1, scales = 'free') +
    theme_test() +
    theme(plot.title = element_text(margin = margin(b = 30))) +
    labs(x = 'Number of genes in signatures', y = 'Correlation between signature scores', title = 'Correlation of pEMT and CAF signatures')

pdf('../data_and_figures/final_figures_resubmission/R5.pdf', width = 10, height = 7)
plot_grid(
    plot_grid(
        scatterplots,
        summary_barplot + theme(plot.margin = unit(c(20, 5.5, 5.5, 5.5), 'pt'), legend.position = 'none'),
        nrow = 2,
        ncol = 1,
        axis = 'l',
        align = 'v'
    ),
    plot_grid(blank_plot(), get_legend(summary_barplot + theme(legend.justification = c(0, 0.9))), nrow = 2, ncol = 1),
    nrow = 1,
    ncol = 2,
    rel_widths = c(11, 1)
)
dev.off()





# Subset of lineplots with "error bars":

sc_metadata_subset <- list(
	coadread = list(
		tcga_cancer_types = c('COAD', 'READ'),
		read_quote = quote(
            fread('../data_and_figures/lee_crc_2020_smc_reclassified.csv')[cell_type != 'ambiguous', -c('cell_type_author', 'cell_type_lenient')]
        ),
		initial_cell_types = c('b_cell', 'endothelial', 'macrophage', 't_cell', 'caf', 'cancer'),
		rare_cell_types = c('epithelial', 'mast'),
        seed = 3361,
        max_mean_count = 1000,
        plot_title = 'Colorectal'
	),
    hnsc = list(
        tcga_cancer_types = 'HNSC',
        read_quote = quote(fread('../data_and_figures/puram_hnscc_2017_reclassified.csv')[cell_type != 'ambiguous', -'cell_type_author']),
		initial_cell_types = c('endothelial', 'macrophage', 'mast', 't_cell', 'caf', 'cancer'),
		rare_cell_types = c('b_cell', 'dendritic', 'myocyte'),
        seed = 8511,
        max_mean_count = 100,
        plot_title = 'Head and Neck'
    ),
	luad = list(
		tcga_cancer_types = 'LUAD',
		read_quote = quote(fread('../data_and_figures/kim_luad_2020_reclassified.csv')[cell_type != 'ambiguous', -'cell_type_author']),
		initial_cell_types = c('b_cell', 'endothelial', 'macrophage', 'mast', 't_cell', 'caf', 'cancer'),
		rare_cell_types = 'dendritic',
        seed = 1096,
        max_mean_count = 1000,
        plot_title = 'Lung Adenocarcinoma'
	)
)

genes_list_subset <- readRDS('../data_and_figures/sc_genes_list.rds')[c('coadread', 'hnsc', 'luad')]

lineplots_sd <- sapply(
    names(sc_metadata_subset),
    function(cohort) {

		cat(paste0(cohort, '\n'))
        sc_data <- eval(sc_metadata_subset[[cohort]]$read_quote)
		sc_data[
            ,
            cell_type := mapvalues(
                cell_type,
                sc_metadata_subset[[cohort]]$rare_cell_types,
                rep('rare', length(sc_metadata_subset[[cohort]]$rare_cell_types))
            )
        ]
        set.seed(sc_metadata_subset[[cohort]]$seed)
        cell_types <- unique(sc_data$cell_type)
    	initial_types <- sc_metadata_subset[[cohort]]$initial_cell_types

        proportions_table <- rbindlist(
            lapply(
    			initial_types,
                # cell_types,
                function(ct) {
                    simulated_counts <- simulate_counts(cell_types, initial_types = ct, max_mean_count = sc_metadata_subset[[cohort]]$max_mean_count)
                    sampled_indices <- sample_indices(sc_data$cell_type, simulated_counts$counts_table[, ..cell_types])
                    proportions_table <- type_contrib(
                        sc_data[, c('cell_type', genes_list_subset[[cohort]]$filtered_cancer_caf), with = FALSE],
                        genes_list_subset[[cohort]]$filtered_cancer_caf,
                        simulated_counts$counts_table,
                        sampled_indices,
                        initial_types = ct,
    					normalise_fun = NULL
                    )[, initial_type := ct]
                    return(proportions_table)
                }
            )
        )

        plot_data <- proportions_table[
            ,
            .(mean_proportion_contrib = mean(proportion_contrib), sd_proportion_contrib = sd(proportion_contrib)),
            by = .(initial_type, proportion_content)
        ]

        plot_data <- rbind(
            plot_data,
            data.table(
                initial_type = unique(plot_data$initial_type),
                proportion_content = rep(0, length(unique(plot_data$initial_type))),
                mean_proportion_contrib = rep(0, length(unique(plot_data$initial_type))),
                sd_proportion_contrib = rep(0, length(unique(plot_data$initial_type)))
            )
        )

        plot_data[, initial_type := factor(initial_type, levels = initial_types)]
    	plot_data <- plot_data[order(initial_type, proportion_content)]

        lineplot <- ggplot(data = plot_data, aes(x = proportion_content, y = mean_proportion_contrib, colour = initial_type)) +
            geom_line() +
            geom_point(shape = 0) +
            geom_errorbar(
                aes(ymin = mean_proportion_contrib - sd_proportion_contrib, ymax = mean_proportion_contrib + sd_proportion_contrib),
                width = 0.01,
                size = 0.25
            ) +
            theme(
                panel.background = element_blank(),
                panel.border = element_rect(fill = NA, colour = 'black'),
                panel.grid.minor = element_blank(),
                panel.grid.major.x = element_blank(),
                panel.grid.major.y = element_line(colour = 'grey', size = 0.25, linetype = 'dotted')
            ) +
            scale_y_continuous(
                limits = c(0, 1),
                breaks = c(0, 0.25, 0.5, 0.75, 1),
                labels = c('0' = '0', '0.25' = '0.25', '0.5' = '0.5', '0.75' = '0.75', '1' = '1'),
                expand = c(0.02, 0.02)
            ) +
            scale_x_continuous(expand = c(0.02, 0.02)) +
            scale_colour_manual(
                labels = c(
    				'b_cell' = 'B cell',
    				'cancer' = 'Cancer',
    				'endothelial' = 'Endothelial',
    				'caf' = 'CAF',
    				'macrophage' = 'Macrophage',
    				'mast' = 'Mast',
    				't_cell' = 'T cell'
    			),
                values = c(
    				'b_cell' = '#8DD3C7',
    				'cancer' = '#FB8072',
    				'endothelial' = '#BC80BD',
    				'caf' = '#FDB462',
    				'macrophage' = '#80B1D3',
    				'mast' = '#FCCDE5',
    				't_cell' = '#B3DE69'
    			)
            ) +
            labs(
                x = 'Proportion of cell mixture',
                y = 'Proportion of gene expression',
                colour = 'Cell type',
                title = sc_metadata_subset[[cohort]]$plot_title
            )

        return(list(lineplot = lineplot, plot_data = plot_data, proportions_table = proportions_table))

    },
    simplify = FALSE,
    USE.NAMES = TRUE
)

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
			'mast' = '#FCCDE5',
            't_cell' = '#B3DE69'
        )
    ) +
    labs(fill = 'Cell type') +
    theme(legend.key = element_rect(size = 1, colour = 'white'), legend.key.size = unit(15, 'pt'))

dummy_legend <- get_legend(dummy_legend_plot)

pdf('../data_and_figures/sim_bulk_lineplots_errorbars.pdf', width = 12, height = 3.5)
plot_grid(
    lineplots_sd$coadread$lineplot + theme(legend.position = 'none'),
    lineplots_sd$hnsc$lineplot + theme(legend.position = 'none') + labs(y = NULL),
    lineplots_sd$luad$lineplot + theme(legend.position = 'none') + labs(y = NULL),
    dummy_legend,
    nrow = 1,
    ncol = 4,
    rel_widths = c(1, 1, 1, 0.5),
    align = 'v',
    axis = 'l'
)
dev.off()
