infer_subtypes <- function(expression_data, subtypes_dt) {

    # expression_data should be in the form returned by prepare_expression_data(), but any transformations of the data to influence the selection of
    # marker genes and the calculation of subtype scores, such as scaling and taking logarithms, should be done beforehand - this function will do no
    # further transformations for these calculations.

    # subtypes_dt must be a data table with first column being tumour IDs and second being the predefined subtype assignments.

    names(subtypes_dt)[1:2] <- c('tumour_id', 'subtype')

    subtypes_dt$tumour_id <- as.character(subtypes_dt$tumour_id)
    subtypes_dt$subtype <- as.character(subtypes_dt$subtype)

    # Find gene signatures for each of the subtypes:
    expression_data$subtype <- sapply(expression_data$id, function(tumour) {
        if(tumour %in% subtypes_dt$tumour_id) {subtypes_dt[tumour_id == tumour, subtype]} else {''}
    })

    subtype_markers <- lapply(
        unique(subtypes_dt$subtype),
        function(sbtp) names(sort(colMeans(expression_data[subtype == sbtp, -c('id', 'subtype')]), decreasing = TRUE)[1:100])
    )

    names(subtype_markers) <- unique(subtypes_dt$subtype)

    # Score each tumour for how closely it corresponds to each of the subtypes.  These scores are just the average expression across all genes in each
    # of the subtype signatures.  Each tumour thus has as many scores as there are subtypes - assign each tumour to the subtype for which it has the
    # highest score.
    subtype_scores <- as.data.table(lapply(subtype_markers, function(markers_vec) {rowMeans(expression_data[, markers_vec, with = FALSE])}))

    # Create data table of inferred TCGA subtype assignments, along with the assignments given in the TCGA paper supplementary data.  Also calculate
    # the minimum difference between each tumour's scores, for use in filtering.
    expression_data[
        ,
        c('inferred_subtype', 'min_diff') := .(
            names(sapply(1:.N, function(i) which.max(subtype_scores[i]))),
            sapply(
                1:.N,
                function(i) {
                    min(
                        as.numeric(subtype_scores[i, which.max(subtype_scores[i]), with = FALSE]) -
                            as.numeric(subtype_scores[i, -which.max(subtype_scores[i]), with = FALSE])
                    )
                }
            )
        )
    ]

    # Output:
    list(
        call = as.list(match.call()),
        inferred_subtypes_dt = expression_data[, .(id, subtype, inferred_subtype, min_diff)],
        subtype_markers = subtype_markers,
        subtype_scores = subtype_scores
    )

}





gene_filter_correlations <- function(
    genes_list,
    expression_data,
    cell_type_markers,
    initial_genes = c('SNAI1', 'SNAI2', 'TWIST1', 'VIM', 'ZEB1', 'ZEB2'),
    cell_types = c('B_plasma', 'myocyte', 'macrophage', 'endothelial', 'DC', 'mast', 'T', 'B')
) {
    as.data.table(
        cbind(
            cor(expression_data[, ..genes_list], expression_data[, ..initial_genes]),
            sapply(
                cell_types,
                function(ct) {
                    rowMeans(
                        cor(
                            expression_data[, ..genes_list],
                            expression_data[, cell_type_markers[cell_type == ct & gene %in% names(expression_data), gene], with = FALSE]
                        )
                    )
                },
                USE.NAMES = TRUE
            )
        ),
        keep.rownames = 'id'
    )[
        ,
        c('max_initial_gene', 'which_max_initial_gene', 'max_cell_type', 'which_max_cell_type') := .(
            max(abs(unlist(mget(initial_genes)))),
            initial_genes[which.max(abs(unlist(mget(initial_genes))))],
            max(unlist(mget(cell_types))),
            cell_types[which.max(unlist(mget(cell_types)))]
        ),
        by = id
    ]
}





deconvolve_emt_caf_data <- function(

    expression_data,
    meta_data,
    genes,
    cell_type_markers,
    sample_ids = NULL,
    tcga_cancer_types = NULL,
    subtypes = NULL,
    subtypes_data = NULL,
	subtypes_var = 'inferred_subtype',
    ref_for_subtypes = NULL,
    ccle_cancer_type = NULL,
    ccle_data = NULL,
    extra_data_source = NULL,
    extra_data = NULL,
    initial_genes = c('SNAI1', 'SNAI2', 'TWIST1', 'VIM', 'ZEB1', 'ZEB2'),
    cell_types = c('B_plasma', 'myocyte', 'macrophage', 'endothelial', 'DC', 'mast', 'T', 'B'),
    genes_from_tcga_fun = function(x) {top_cols_by_fun_cor(x, initial = initial_genes)[1:min(.N, 300), id]},
    genes_filter_fun = function(x) 1:250,
    gene_weights_fun = median,
    cell_type_weights = NULL,
    initial_gene_weights = NULL,
    caf_markers = c('COL1A1', 'COL1A2', 'COL3A1', 'COL6A3', 'THY1'),
    seed = NULL,
    ordering_fun = function(x) {as.vector(seriation::seriate(as.dist(1 - x), method = 'SPIN_STS')[[1]])}

) {

    # <caf_markers> is used for tie-breaking in deciding whether to reverse the ordering for the correlation matrix.

    if(!is.null(genes_filter_fun)) {genes_filter_fun <- match.fun(genes_filter_fun)}
    gene_weights_fun <- match.fun(gene_weights_fun)
    ordering_fun <- match.fun(ordering_fun)

    setkey(expression_data, id)
    setkey(meta_data, id)

    # Determine which of sample_ids, tcga_cancer_types and subtypes can be used:
    if(is.null(sample_ids) & is.null(tcga_cancer_types)) {
        stop('Please specify sample IDs or TCGA cancer types.')
    } else if((!is.null(sample_ids) & !is.null(tcga_cancer_types)) | (!is.null(sample_ids) & !is.null(subtypes))) {
        warning('If specifying sample_ids, no need to specify tcga_cancer_types or subtypes.  Using only sample_ids.')
        tcga_cancer_types <- NULL
        subtypes <- NULL
    } else if(!is.null(tcga_cancer_types) & !is.null(subtypes) & is.null(subtypes_data)) {
        stop('Please supply subtypes_data.')
    }

    # If we're not using sample_ids by default, get the appropriate sample IDs using tcga_cancer_types, and subtypes if it is supplied:
    if(is.null(sample_ids)) { # Then we must have !is.null(tcga_cancer_types)
        if(is.null(subtypes)) { # Then get all IDs corresponding to these cancer types
            sample_ids <- meta_data[cancer_type %in% tcga_cancer_types, id]
        } else { # Then we must have !is.null(subtypes_data)
            setkey(subtypes_data, id)
            if(is.null(ref_for_subtypes)) {
                sample_ids <- subtypes_data[meta_data[cancer_type %in% tcga_cancer_types, id]][get(subtypes_var) %in% subtypes, id]
            } else {
                # Optionally supply reference for subtypes, in case the same subtype name appears in multiple references.
                sample_ids <- subtypes_data[meta_data[cancer_type %in% tcga_cancer_types, id]][
                    subtype_ref == ref_for_subtypes & get(subtypes_var) %in% subtypes,
                    id
                ]
            }

        }

    }

    # Make sure the samples we're using have corresponding purity data:
    sample_ids <- meta_data[sample_ids][!is.na(purity), id]

    # Subset expression_data and remove genes with zero standard deviation:
    expression_data <- cbind(
        expression_data[sample_ids, .(id)],
        expression_data[sample_ids, sapply(.SD[, -'id'], function(x) {switch((sd(x) > 0) + 1, NULL, x)})]
    )

    # Subset genes which occur in expression_data and, if provided, ccle_data and in any extra data:
    genes <- genes[genes %in% names(expression_data)]

    # Get EMT markers by correlation with initial genes:
    if(!is.null(genes_from_tcga_fun)) {
        top_genes <- genes_from_tcga_fun(expression_data[, -'id'])
        genes <- unique(sort(c(genes, top_genes)))
    }

    # In case ccle and any extra data were provided, subset genes which occur in these datasets:

    if(!is.null(ccle_cancer_type) & !is.null(ccle_data)) {
        genes <- genes[genes %in% ccle_data$gene_id]
    } else if(!is.null(ccle_cancer_type) & is.null(ccle_data) | is.null(ccle_cancer_type) & !is.null(ccle_data)) {
        warning('Only one of ccle_cancer_type and ccle_data has been supplied.  Please supply both.')
    }

    if(!is.null(extra_data_source) & !is.null(extra_data)) {
        genes <- genes[genes %in% extra_data[source == extra_data_source, gene]]
    } else if(!is.null(extra_data_source) & is.null(extra_data) | is.null(extra_data_source) & !is.null(extra_data)) {
        warning('Only one of extra_data_source and extra_data has been supplied.  Please supply both.')
    }





    # Optionally filter this gene list:

    if(!is.null(genes_filter_fun)) {

        genes_cor_with_initial_and_cell_types <- gene_filter_correlations(
            genes,
            expression_data,
            cell_type_markers,
            initial_genes = initial_genes,
            cell_types = cell_types
        )

        # The following is the actual filtering step.  Note that the scores are designed so that a high score is good - we want to keep genes with
        # high scores.

        if(isFALSE(initial_gene_weights) & isFALSE(cell_type_weights)) {
            initial_gene_weights <- NULL
            cell_type_weights <- NULL
            warning('At least one of initial_gene_weights and cell_type_weights should not be FALSE.  Calculating all weights.')
        }

        # The user can supply a named list of weights for the initial genes and the cell types (the names must correspond to the initial_genes, resp.
        # cell_types arguments).  If they don't, we construct them below.

        if(is.null(initial_gene_weights)) {

            # Define weights for initial genes as quantiles of their correlations with other cell types (note we use '1 -' because we want to give
            # more weight to those initial genes that don't correlate highly with other cell types):

            initial_gene_weights <- as.list(
                1 - apply(
                    sapply(
                        cell_types,
                        function(ct) {
                            expression_data[
                                sample_ids,
                                rowMeans(
                                    cor(
                                        .SD[, ..initial_genes],
                                        .SD[, cell_type_markers[cell_type == ct & gene %in% names(expression_data), gene], with = FALSE]
                                    )
                                )
                            ]
                        },
                        USE.NAMES = TRUE
                    ),
                    1,
                    gene_weights_fun
                )
            )

        }

        if(is.null(cell_type_weights)) {

            cell_type_weights <- as.list(
                apply(
                    sapply(
                        cell_types,
                        function(ct) {
                            expression_data[
                                sample_ids,
                                rowMeans(
                                    cor(
                                        .SD[, ..genes],
                                        .SD[, cell_type_markers[cell_type == ct & gene %in% names(expression_data), gene], with = FALSE]
                                    )
                                )
                            ]
                        },
                        USE.NAMES = TRUE
                    ),
                    2,
                    gene_weights_fun
                )
            )

        }

        if(isFALSE(cell_type_weights)) {

            gene_scores_table <- genes_cor_with_initial_and_cell_types[
                ,
                .(score = sum(as.numeric(.SD[, ..initial_genes])*as.numeric(initial_gene_weights[initial_genes]))),
                by = id
            ][order(-score)]

            genes_filtered <- gene_scores_table[genes_filter_fun(score), id]

        } else if(isFALSE(initial_gene_weights)) {

            gene_scores_table <- genes_cor_with_initial_and_cell_types[
                ,
                .(score = -sum(as.numeric(.SD[, ..cell_types])*as.numeric(cell_type_weights[cell_types]))),
                by = id
            ][order(-score)]

            genes_filtered <- gene_scores_table[genes_filter_fun(score), id]

        } else {

            gene_scores_table <- genes_cor_with_initial_and_cell_types[
                ,
                .(
                    score = sum(as.numeric(.SD[, ..initial_genes])*as.numeric(initial_gene_weights[initial_genes])) -
                        sum(as.numeric(.SD[, ..cell_types])*as.numeric(cell_type_weights[cell_types]))
                ),
                by = id
            ][order(-score)]

            genes_filtered <- gene_scores_table[genes_filter_fun(score), id]

        }

        genes_filtered <- unique(c(genes_filtered, initial_genes))

    } else {

        genes_filtered <- genes

    }





    # Regress genes in filtered list against sum of genes in filtered list, and take residuals:
    resid_data <- as.data.table(
        sapply(
            genes_filtered,
            function(g) {
                # In the following, I have to put the backticks `` around g because there are genes like 'NKX2-1', which gets interpreted as NKX2 - 1.
                lm(
                    formula(paste0('`', g, '` ~ row_sums')),
                    data = cbind(expression_data[, ..genes_filtered], row_sums = rowSums(expression_data[, ..genes_filtered]))
                )$residuals
            },
            simplify = FALSE,
            USE.NAMES = TRUE
        )
    )

    # Correlation matrix for resid_data:
    cor_mat <- cor(resid_data)

    # Ordering of correlation matrix for resid_data using SPIN_STS:
    if(!is.null(seed)) {set.seed(seed)}
    ordering <- ordering_fun(cor_mat)





    # Detect which genes are cancer EMT genes and which are CAF genes by looking at correlations with purity and, if available, CCLE data and extra
    # data, taking a consensus between them.  Break ties using a few reliable CAF markers, checking which end they're closest to.

    cor_with_purity_scale <- cor(
        scaledt(expression_data[, c('id', ..genes_filtered)], margin = 1, scale = FALSE)[, -'id'],
        meta_data[sample_ids, purity]
    )[, 1]
    cor_with_purity_resid <- cor(resid_data, meta_data[sample_ids, purity])[, 1]
    cor_with_purity_raw <- cor(expression_data[, ..genes_filtered], meta_data[sample_ids, purity])[, 1]

    # CCLE comparison:
    if(!is.null(ccle_cancer_type) & !is.null(ccle_data)) {
        ccle_data <- expression_data[
            ,
            .(id = names(.SD), tumours = colMeans(.SD), ccle = ccle_data[names(expression_data[, -'id']), ..ccle_cancer_type][[ccle_cancer_type]]),
            .SDcols = -'id'
        ]
        mod_loess <- loess(tumours ~ ccle, ccle_data, span = 0.25, degree = 1, family = 'symmetric')
        setkey(ccle_data, id)
        ccle_comp_diff <- predict(mod_loess, ccle_data[genes_filtered, ccle]) + 1 - colMeans(expression_data[, ..genes_filtered])
    }

    # Extra data:
    if(!is.null(extra_data_source) & !is.null(extra_data)) {
        extra_data_score <- extra_data[source == extra_data_source][genes_filtered, setNames(diff, gene)] # Need a named vector for heat_map_bar()
    }

    # Store variables for checking the ordering in a list:
    ordering_check <- list(cor_with_purity_raw)
    if(exists('ccle_comp_diff')) {ordering_check <- c(ordering_check, list(ccle_comp_diff))}
    if(exists('extra_data_score')) {ordering_check <- c(ordering_check, list(extra_data_score))}

    # This is the actual check:
    ordering_check <- sum(
        sapply(
            ordering_check,
            function(x) {as.numeric(lm(x[ordering] ~ seq(1, length(genes_filtered), 1))$coeff[2])}
        ) > 0 # Greater than zero suggests the ordering is wrong
    )/length(ordering_check)

    # If a consensus suggests the ordering is wrong, reverse the ordering; if there is no consensus, try to break the tie using the average position
    # of the CAF markers (if this still doesn't break the tie, or if these CAF markers don't appear in the filtered gene list, then nothing happens):
    if(ordering_check > 0.5) {
        ordering <- rev(ordering)
    } else if(ordering_check == 0.5) {
        if(sum(caf_markers %in% genes_filtered) > 0) {
            if(mean(which(genes_filtered[ordering] %in% caf_markers)) < length(genes_filtered)/2) {ordering <- rev(ordering)}
        }
    }





    # Output:

    out <- list(
        sample_ids = sample_ids,
        genes_filtered = genes_filtered,
        cor_mat = cor_mat,
        ordering = ordering,
        cor_with_purity = list(scale = cor_with_purity_scale, resid = cor_with_purity_resid, raw = cor_with_purity_raw)
    )

    if(!is.null(genes_filter_fun)) {
        setkey(genes_cor_with_initial_and_cell_types, id)
        out <- c(
            out,
            list(
                genes_filter_fun = genes_filter_fun,
                cor_with_initial_and_cell_types = genes_cor_with_initial_and_cell_types,
                scores_table = gene_scores_table,
                cell_type_lms = genes_cor_with_initial_and_cell_types[genes_filtered][
                    ordering,
                    sapply(.SD, function(ct) {setNames(lm(ct ~ .I)$coeff, c('intercept', 'slope'))}, simplify = TRUE, USE.NAMES = TRUE),
                    .SDcols = cell_types
                ]
            )
        )
    }

    if(!is.null(ccle_cancer_type) & !is.null(ccle_data)) {
        out <- c(out, list(ccle_data = ccle_data, mod_loess = mod_loess, ccle_comp_diff = ccle_comp_diff))
    }

    if(!is.null(extra_data_source) & !is.null(extra_data)) {
        out <- c(out, list(extra_data_score = extra_data_score))
    }

    out <- c(
        out,
        list(
            tcga_cancer_types = tcga_cancer_types,
            subtypes = subtypes,
            ccle_cancer_type = ccle_cancer_type,
            cell_types = cell_types,
            initial_genes = initial_genes,
            caf_markers = caf_markers
        )
    )

    out

}





deconvolve_emt_caf_plots <- function(

    data,
    plot_title = '',
    heatmap_axis_title = 'genes',
    heatmap_legend_title = 'correlation',
    heatmap_colours = rev(colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))(50)),
    heatmap_colour_limits = c(-0.6, 0.6),
    heatmap_legend_breaks = c(-0.6, -0.3, 0, 0.3, 0.6),
    heatmap_legend_labels = waiver(),
	heatmap_legend_justification = 'left',
    heatmap_annotations = data$initial_genes,
	heatmap_annotations_side = c('bottom', 'left', 'top', 'right'),
    heatmap_annotations_nudge = 0.35,
    heatmap_annotations_text_size = 3,
    heatmap_annotations_segment_size = 0.5,
    heatmap_annotations_force = 0.25,
    heatmap_annotations_box_padding = 0.25,
    purity_colours = c('darkslategrey', 'turquoise4', 'yellow', 'chocolate1', 'coral4'),
    purity_fun = function(x) caTools::runmean(x, 30),
    purity_cor_method = c('scale', 'resid', 'raw'),
    purity_axis_title = 'Correlation\nwith purity',
    purity_legend_title = NULL,
    purity_colour_limits = c(-0.3, 0.3),
    purity_legend_breaks = c(-0.3, 0, 0.3),
	purity_legend_labels = waiver(),
    purity_legend_direction = c('horizontal', 'vertical'),
    ccle_colours = purity_colours,
    ccle_fun = function(x) {caTools::runmean(x, 30)/max(abs(caTools::runmean(x, 30)))},
    ccle_axis_title = 'Tumours vs.\ncell lines',
    ccle_legend_title = NULL,
    ccle_colour_limits = c(-1, 1),
    ccle_legend_breaks = c(-1, 0, 1),
	ccle_legend_labels = waiver(),
    ccle_legend_direction = c('horizontal', 'vertical'),
    extra_colours = ccle_colours,
    extra_fun = ccle_fun,
    extra_colour_limits = ccle_colour_limits,
    extra_legend_breaks = ccle_legend_breaks,
	extra_legend_labels = ccle_legend_labels,
    extra_legend_direction = c('horizontal', 'vertical'),
    extra_axis_title = 'scRNA-seq',
    extra_legend_title = NULL,
    bar_legend_height = unit(10, 'pt'),
    bar_legend_width = unit(10, 'pt'),
    bar_legend_justification = 'center',
    diagnostics = TRUE,
    diagnostic_colours = grDevices::rainbow(50),
    expression_data = NULL,
    epithelial_markers = c(
        'CDH1',
        'EPCAM',
        'SFN',
        switch(is.null(expression_data) + 1, names(expression_data)[grep('^KRT[0-9]|^KRTD', names(expression_data))], NULL)
    )

) {

    purity_cor_method = match.arg(purity_cor_method)

	heatmap_annotations_side <- match.arg(heatmap_annotations_side)
    purity_legend_direction <- match.arg(purity_legend_direction)
    ccle_legend_direction <- match.arg(ccle_legend_direction)
    extra_legend_direction <- match.arg(extra_legend_direction)

    purity_fun <- match.fun(purity_fun)
    ccle_fun <- match.fun(ccle_fun)
    extra_fun <- match.fun(extra_fun)

    cor_mat <- data$cor_mat
    ordering <- data$ordering
    genes_filtered <- data$genes_filtered
    genes_cor_with_initial_and_cell_types <- data$cor_with_initial_and_cell_types

    # Heatmap, annotations and purity bar:

    htmp <- heat_map(
        cor_mat,
        ordering,
        axis_title_y = heatmap_axis_title,
        axis_title_x = '',
        axis_text_x = '',
        axis_text_y = '',
        legend_title = heatmap_legend_title,
        colours = heatmap_colours,
        colour_limits = heatmap_colour_limits,
        legend_breaks = heatmap_legend_breaks,
        legend_labels = heatmap_legend_labels,
        plot_margin = c(5.5, switch((heatmap_annotations_side == 'right') + 1, 5.5, 0), 0, switch((heatmap_annotations_side == 'left') + 1, 5.5, 0)),
        plot_title = plot_title,
		legend.justification = heatmap_legend_justification
    )

    interesting_labels <- genes_filtered[ordering]

    interesting_labels[!(interesting_labels %in% heatmap_annotations)] <- ''

    axis_labels <- heat_map_labels_repel(
        interesting_labels,
		edge = plyr::mapvalues(
            heatmap_annotations_side,
            c('bottom', 'left', 'top', 'right'),
            c('top', 'right', 'bottom', 'left'),
            warn_missing = FALSE
        ),
        nudge = heatmap_annotations_nudge,
        text_size = heatmap_annotations_text_size,
        segment_size = heatmap_annotations_segment_size,
        force = heatmap_annotations_force,
        box.padding = heatmap_annotations_box_padding
    )

    purity_bar <- heat_map_bar(
        data$cor_with_purity[[purity_cor_method]],
        ordering,
        fun = purity_fun,
        colours = purity_colours,
        colour_limits = purity_colour_limits,
        axis_title_y = purity_axis_title,
        legend_title = purity_legend_title,
        legend_breaks = purity_legend_breaks,
		legend_labels = purity_legend_labels,
        legend.direction = purity_legend_direction,
        legend.justification = bar_legend_justification,
        legend.key.height = bar_legend_height,
        legend.key.width = bar_legend_width,
        axis.title.y = element_text(angle = 0, vjust = 0.5),
        axis.title.x = element_blank(),
        plot_margin = c(0, switch((heatmap_annotations_side == 'right') + 1, 5.5, 0), 2, switch((heatmap_annotations_side == 'left') + 1, 5.5, 0))
    )

    plots <- list(heatmap = htmp, axis_labels = axis_labels, purity_bar = purity_bar)

    # CCLE bar:

    if('ccle_comp_diff' %in% names(data)) {

        ccle_bar <- heat_map_bar(
            data$ccle_comp_diff,
            ordering,
            fun = ccle_fun,
            colours = ccle_colours,
            colour_limits = ccle_colour_limits,
            axis_title_y = ccle_axis_title,
            legend_title = ccle_legend_title,
            legend_breaks = ccle_legend_breaks,
			legend_labels = ccle_legend_labels,
            legend.direction = ccle_legend_direction,
            legend.justification = bar_legend_justification,
            legend.key.height = bar_legend_height,
            legend.key.width = bar_legend_width,
            axis.title.y = element_text(angle = 0, vjust = 0.5),
            axis.title.x = element_blank(),
            plot_margin = c(0, switch((heatmap_annotations_side == 'right') + 1, 5.5, 0), 2, switch((heatmap_annotations_side == 'left') + 1, 5.5, 0))
        )

        plots <- c(plots, ccle_bar = list(ccle_bar))

    }

    # Extra data bar:

    if('extra_data_score' %in% names(data)) {

        extra_bar <- heat_map_bar(
            data$extra_data_score,
            ordering,
            fun = extra_fun,
            colours = extra_colours,
            colour_limits = extra_colour_limits,
            axis_title_y = extra_axis_title,
            legend_title = extra_legend_title,
            legend_breaks = extra_legend_breaks,
			legend_labels = extra_legend_labels,
            legend.direction = extra_legend_direction,
            legend.justification = bar_legend_justification,
            legend.key.height = bar_legend_height,
            legend.key.width = bar_legend_width,
            axis.title.y = element_text(angle = 0, vjust = 0.5),
            axis.title.x = element_blank(),
            plot_margin = c(0, switch((heatmap_annotations_side == 'right') + 1, 5.5, 0), 2, switch((heatmap_annotations_side == 'left') + 1, 5.5, 0))
        )

        plots = c(plots, extra_bar = list(extra_bar))

    }





    # Optionally make diagnostic plots:

    if(diagnostics) {

        diagnostics <- list(

            alternative_purity_cor = sapply(
                c('scale', 'resid', 'raw')[c('scale', 'resid', 'raw') != purity_cor_method],
                function(cor_type) {
                    heat_map_bar(
                        data$cor_with_purity[[cor_type]],
                        ordering,
                        fun = purity_fun,
                        colours = diagnostic_colours,
                        colour_limits = purity_colour_limits,
                        axis_title_y = paste0('Cor. with\npurity (', cor_type, ')'),
                        legend_breaks = purity_legend_breaks,
                        legend.direction = 'horizontal',
                        legend.title = element_blank(),
                        legend.key.height = bar_legend_height,
                        legend.key.width = bar_legend_width,
                        axis.title.y = element_text(angle = 0, vjust = 0.5),
                        axis.title.x = element_blank(),
                        plot_margin = c(0, 5.5, 2, 5.5)
                    )
                },
                simplify = FALSE,
                USE.NAMES = TRUE
            )

        )

        if(!is.null(data$genes_filter_fun)) {

            setkey(genes_cor_with_initial_and_cell_types, id)

            if(xor(is.null(epithelial_markers), is.null(expression_data))) {
                warning('To include epithelial cells in diagnostics, please provide both epithelial markers and expression data.')
            }

            if(!is.null(epithelial_markers) & !is.null(expression_data)) {

                # Remove columns from expression_data which contain NAs or have zero standard deviation:

                expression_data <- cbind(
                    expression_data[data$sample_ids, .(id)],
                    expression_data[data$sample_ids, sapply(.SD[, -'id'], function(x) {switch((any(is.na(x)) || sd(x) == 0) + 1, x, NULL)})]
                )

                # Subset those epithelial markers that still appear in expression_data after removing columns:

                epithelial_markers <- epithelial_markers[epithelial_markers %in% names(expression_data)]

                genes_cor_with_initial_and_cell_types[
                    genes_filtered,
                    epithelial := expression_data[, rowMeans(cor(.SD[, ..genes_filtered], .SD[, ..epithelial_markers]))]
                ]

                cell_types <- c(data$cell_types, 'epithelial')

            } else {

                cell_types <- data$cell_types

            }

            cell_type_bars_colour_limits <- genes_cor_with_initial_and_cell_types[
                genes_filtered,
                quantile(apply(.SD, 2, function(x) purity_fun(x[ordering])), c(0, 1)),
                .SDcols = cell_types
            ]

            diagnostics <- c(

                diagnostics,

                list(

                    cell_type_bars = sapply(
                        cell_types,
                        function(ct) {
                            heat_map_bar(
                                genes_cor_with_initial_and_cell_types[genes_filtered, setNames(get(ct), id)],
                                ordering,
                                fun = purity_fun,
                                colours = grDevices::rainbow(50),
                                colour_limits = as.numeric(cell_type_bars_colour_limits),
                                legend_breaks = c(
                                    ceiling(as.numeric(cell_type_bars_colour_limits)[1]*10)/10,
                                    floor(as.numeric(cell_type_bars_colour_limits)[2]*10)/10
                                ),
                                axis_title_y = ct,
                                legend.direction = 'horizontal',
                                legend.title = element_blank(),
                                legend.key.height = bar_legend_height,
                                legend.key.width = bar_legend_width,
                                axis.title.y = element_text(angle = 0, vjust = 0.5),
                                axis.title.x = element_blank(),
                                plot_margin = c(0, 5.5, 2, 5.5)
                            )
                        },
                        simplify = FALSE,
                        USE.NAMES = TRUE
                    ),

                    initial_genes_cell_types_cor = ggplot(
                        melt(
                            genes_cor_with_initial_and_cell_types[data$initial_genes, c('id', cell_types), with = FALSE],
                            id.vars = 'id',
                            variable.name = 'cell_type',
                            value.name = 'correlation'
                        )
                    ) +
                        geom_raster(aes(x = cell_type, y = id, fill = correlation)) +
                        scale_y_discrete(expand = c(0, 0)) +
                        scale_x_discrete(expand = c(0, 0)),

                    gene_cell_type_scores_density = ggplot(data$scores_table) + geom_density(aes(score)) + theme_test()

                )

            )

        }

        if(!is.null(data$ccle_cancer_type) & !is.null(data$ccle_data)) {
            diagnostics <- c(
                diagnostics,
                tumours_vs_cell_lines = list(
                    ggplot(cbind(data$ccle_data, in_genes_filtered = data$ccle_data$id %in% genes_filtered)) +
                        geom_point(aes(x = ccle, y = tumours, col = in_genes_filtered)) +
                        scale_colour_manual(values = c('black', 'deepskyblue')) +
                        geom_point(aes(x = ccle, y = tumours), data = data$ccle_data[id %in% genes_filtered], colour = 'deepskyblue') +
                        geom_line(aes(x = ccle, y = predict(data$mod_loess, newdata = data$ccle_data)), colour = 'red', size = 1.5) +
                        labs(colour = 'In filtered\ngenes list')
                )
            )
        }

    }

    # Output:
    list(plots = plots, diagnostics = diagnostics)

}





deconv_rank <- function(deconv_data) {
    genes_to_rank <- sort(unique(unlist(lapply(deconv_data, `[[`, 'genes_filtered'))))
    sapply(
        deconv_data,
        function(li) {
            li_genes <- li$genes_filtered[li$ordering]
            sapply(genes_to_rank, function(g) {if(g %in% li_genes) {which(li_genes == g)/length(li_genes)} else {0.5}})
        }
    )
}





deconv_scores <- function(

    expression_data,
    deconv_data,
    head_tail = c('both', 'head', 'tail'),
    n = 20,
    additional_genes = NULL,
    genes = NULL,
    transform_data = FALSE,
    scale_fun = NULL,
    scale_fun_margin = NULL

) {

    # By default, scores will be calculated for all genes in <deconv_data>.  The user can supply a character vector to the <additional_genes>
    # argument, in which case scores will be calculated for genes in <additional_genes> as well as for those in <deconv_data>. If the user supplies a
    # character vector to the <genes> argument, then scores will only be calculated for the genes in <genes> and not for those in <deconv_data> (or
    # for any supplied to <additional_genes>).  In every case, the scores will be calculated using average correlation with genes in <deconv_data>.
    # If <head_or_tail> is set to 'head' or 'tail', we use respectively the genes from the head and tail of the genes in each item in <deconv_data>.
    # If <head_or_tail> is 'both' (the default), we calculate the difference between the average correlations with the head and tail genes.  The
    # argument <n> is the n to use in the head() and tail() functions.

    # If <transform_data> is TRUE, correlations will be calculated in the transformed space defined by the genes and samples in <deconv_data>, that
    # is, for each cancer type, the gene vectors are replaced with the residuals of a linear regression against sample sums.

    # <scale_fun> is a function to be applied to the genes by cancer types matrix of scores.  For example, scale_fun = scale would scale the scores
    # relative to the global mean and standard deviation; scale_fun = function(x) apply(x, 2, scale) would scale per cancer type;
    # scale_fun = function(x) apply(x, 2, function(y) y/max(abs(y))) would scale the scores for each cancer type to the interval [-1, 1].  You can set
    # <scale_fun_margin> to 1 or 2 if you want to apply <scale_fun> along the margins of the matrix.  1 means applying to genes, 2 means applying to
    # cancer types.

    head_tail <- match.arg(head_tail)
    if(head_tail == 'both') head_tail <- c('head', 'tail')

    if(!is.null(scale_fun)) {scale_fun <- match.fun(scale_fun)}

    setkey(expression_data, id)

    if(!is.null(genes)) {
        genes_to_score <- genes
    } else if(!is.null(additional_genes)) {
        genes_to_score <- unique(c(unlist(lapply(deconv_data, `[[`, 'genes_filtered')), additional_genes))
    } else {
        genes_to_score <- unique(unlist(lapply(deconv_data, `[[`, 'genes_filtered')))
    }

    scores <- sapply(

        deconv_data,

        function(li) {

            li_sample_ids <- li$sample_ids
            li_ordered_genes <- li$genes_filtered[li$ordering]

            for(fun_name in head_tail) {assign(paste0(fun_name, '_genes'), match.fun(fun_name)(li_ordered_genes, n))}

            genes_for_scores_data <- unique(c(genes_to_score, unlist(mget(paste0(head_tail, '_genes')))))

            if(transform_data) {
                scores_data <- cbind(
                    expression_data[li_sample_ids, ..genes_for_scores_data],
                    row_sums = rowSums(expression_data[li_sample_ids, ..li_ordered_genes])
                )
                scores_data <- as.data.table(
                    sapply(
                        genes_for_scores_data,
                        function(g) {lm(formula(paste0('`', g, '` ~ row_sums')), data = scores_data)$residuals},
                        simplify = FALSE,
                        USE.NAMES = TRUE
                    )
                )
            } else {
                scores_data <- expression_data[li_sample_ids, ..genes_for_scores_data]
            }

            if(length(head_tail) == 2) {
                scores_li <- scores_data[
                    ,
                    (rowMeans(cor(.SD[, ..genes_to_score], .SD[, ..head_genes])) - rowMeans(cor(.SD[, ..genes_to_score], .SD[, ..tail_genes])))/2
                ]
            } else {
                scores_li <- scores_data[, rowMeans(cor(.SD[, ..genes_to_score], .SD[, get(paste0(head_tail, '_genes')), with = FALSE]))]
            }

            scores_li

        },

        simplify = TRUE,
        USE.NAMES = TRUE

    )

    if(!is.null(scale_fun)) {
        if(!is.null(scale_fun_margin)) {
            if(scale_fun_margin == 1) {
                scores <- t(apply(scores, 1, scale_fun))
            } else if(scale_fun_margin == 2) {
                scores <- apply(scores, 2, scale_fun)
            } else {
                stop('scale_fun_margin must be 1, 2 or NULL.')
            }
        } else {
            scores <- scale_fun(scores)
        }
    }
    scores <- as.data.table(scores, keep.rownames = 'gene')
    setkey(scores, gene)

    scores

}





deconv_heatmap <- function(

    scores_data,
    order_genes_fun = c('hclust', 'seriate'),
    order_genes_method = 'average',
    order_analyses_fun = c('hclust', 'seriate'),
    order_analyses_method = 'average',
    colours = rev(colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))(50)),
    colour_limits = c(-1, 1),
    legend_breaks = c(-1, -0.5, 0, 0.5, 1),
    legend_labels = c('-1' = '-1', '-0.5' = '-0.5', '0' = '0', '0.5' = '0.5', '1' = '1'),
    x_axis_text_angle = 55,
    legend_title = 'pEMT-CAF\nscore',
    draw_dendro = TRUE,
    plot_margin = c(5.5, 5.5, 5.5, 5.5),
    plot_title = NULL,
    ...

) {

    # <order_genes_fun> and <order_analyses_fun> can each be 'hclust', 'seriate' or a user-defined function.  In the former two cases, hclust() and
    # seriate() functions will be used respectively, and the corresponding methods may be selected using the <order_genes_method> and
    # <order_analyses_method> arguments.  An hclust object will not be generated if the user defines their own function.

    # The '...' argument is for extra arguments to the theme() function in the construction of the heatmap.

    if(!(typeof(order_genes_fun) == 'character' && order_genes_fun %in% c('hclust', 'seriate'))) {
        order_genes_fun <- match.fun(order_genes_fun)
    }

    if(!(typeof(order_analyses_fun) == 'character' && order_analyses_fun %in% c('hclust', 'seriate'))) {
        order_analyses_fun <- match.fun(order_analyses_fun)
    }

    setkey(scores_data, gene)

    # Since the genes in scores_data are in alphabetical order, ordering_genes is the ordering with respect to the alphabetical order.

    if(typeof(order_genes_fun) == 'character' && order_genes_fun == 'seriate') {

        if(order_genes_method %in% seriation::list_seriation_methods('matrix')) {
            ordering_obj_genes <- seriation::seriate(as.matrix(scores_data[, -'gene']), method = order_genes_method, margin = 1)
        } else if(order_genes_method %in% seriation::list_seriation_methods('dist')) {
            dist_mat_genes <- dist(scores_data[, -'gene'])
            ordering_obj_genes <- seriation::seriate(dist_mat_genes, method = order_genes_method)
        } else {
            stop("If using 'seriate' for ordering, please select a seriation method for either 'dist' or 'matrix'.")
        }

        ordering_genes <- seriation::get_order(ordering_obj_genes)

        if('hclust' %in% class(ordering_obj_genes[[1]])) {hclust_genes <- ordering_obj_genes[[1]]}

    } else if(typeof(order_genes_fun) == 'character' && order_genes_fun == 'hclust') {
        dist_mat_genes <- dist(scores_data[, -'gene'])
        hclust_genes <- hclust(dist_mat_genes, method = order_genes_method)
        ordering_genes <- hclust_genes$order
    } else {
        ordering_genes <- order_genes_fun(scores_data[, -'gene'])
    }

    # ordering_analyses is with respect to the ordering in scores_data, which should be  the same as the ordering in the deconv_data used as input to
    # the deconv_scores() function.

    if(typeof(order_analyses_fun) == 'character' && order_analyses_fun == 'seriate') {

        if(order_analyses_method %in% seriation::list_seriation_methods('matrix')) {
            ordering_obj_analyses <- seriation::seriate(as.matrix(scores_data[, -'gene']), method = order_analyses_method, margin = 2)
        } else if(order_analyses_method %in% seriation::list_seriation_methods('dist')) {
            dist_mat_analyses <- dist(t(scores_data[, -'gene']))
            ordering_obj_analyses <- seriation::seriate(dist_mat_analyses, method = order_analyses_method)
        } else {
            stop("If using 'seriate' for ordering, please select a seriation method for either 'dist' or 'matrix'.")
        }

        ordering_analyses <- seriation::get_order(ordering_obj_analyses)

        if('hclust' %in% class(ordering_obj_analyses[[1]])) {hclust_analyses <- ordering_obj_analyses[[1]]}

    } else if(typeof(order_analyses_fun) == 'character' && order_analyses_fun == 'hclust') {
        dist_mat_analyses <- dist(t(scores_data[, -'gene']))
        hclust_analyses <- hclust(dist_mat_analyses, method = order_analyses_method)
        ordering_analyses <- hclust_analyses$order
    } else {
        ordering_analyses <- order_analyses_fun(scores_data[, -'gene'])
    }

    scores_melted <- melt(
        scores_data,
        id.vars = 'gene',
        measure.vars = names(scores_data[, -'gene'])[ordering_analyses],
        variable.name = 'analysis',
        value.name = 'score',
        variable.factor = FALSE
    )

    setkey(scores_melted, gene)

    # Create heatmap:

    htmp <- ggplot(
        scores_melted[scores_data$gene[ordering_genes]],
        aes(x = factor(analysis, levels = unique(analysis)), y = factor(gene, levels = unique(gene)))
    ) +
        geom_raster(aes(fill = score)) +
        scale_x_discrete(expand = c(0, 0)) +
        scale_y_discrete(expand = c(0, 0)) +
        scale_fill_gradientn(colours = colours, limits = colour_limits, oob = scales::squish, breaks = legend_breaks, labels = legend_labels) +
        theme(
            axis.ticks = element_blank(),
            axis.text.x = element_text(angle = x_axis_text_angle, hjust = 1),
            axis.title.y = element_blank(),
            axis.ticks.length = unit(0, 'pt'),
            plot.margin = unit(plot_margin, 'pt'),
            ...
        ) +
        labs(title = plot_title, x = '', fill = legend_title)

    out <- list(
        heatmap = htmp,
        genes = scores_data$gene,
        analyses = names(scores_data[, -'gene']),
        ordering_genes = ordering_genes,
        ordering_analyses = ordering_analyses
    )

    if(exists('hclust_genes')) {
        if(draw_dendro) {out <- c(out, list(dendro_genes = dendro(hclust_genes, edge = 'left')))}
        out <- c(out, list(hclust_genes = hclust_genes))
    }

    if(exists('hclust_analyses')) {
        if(draw_dendro) {out <- c(out, list(dendro_analyses = dendro(hclust_analyses, edge = 'bottom')))}
        out <- c(out, list(hclust_analyses = hclust_analyses))
    }

    if(exists('dist_mat_genes')) {out <- c(out, list(dist_mat_genes = dist_mat_genes))}

    if(exists('dist_mat_analyses')) {out <- c(out, list(dist_mat_analyses = dist_mat_analyses))}

    out

}





# The following function plots the aligned heatmaps and dendrogram(s) from the output of the deconv_heatmap() function.

deconv_heatmap_dendro_plot <- function(
    deconv_heatmap_list,
    show_dendros = na.omit(sapply(strsplit(names(deconv_heatmap_list), 'dendro_'), `[`, 2)),
    rel_heights = c(1, 8),
    rel_widths = switch((length(show_dendros) == 2 || show_dendros == 'analyses') + 1, c(6, 1, 1), c(6, 1)),
    padding = unit(0, 'pt'),
    ...
) {

    # The '...' argument is for additional arguments for guide_colourbar(), which can be used to change the formatting of the legend.

    if(sum(startsWith(names(deconv_heatmap_list), 'dendro')) == 0) {stop('There are no dendrograms in this plot list.')}

    show_dendros <- match.arg(show_dendros, c('genes', 'analyses'), several.ok = TRUE)
    grob_htmp <- ggplotGrob(deconv_heatmap_list$heatmap + theme(legend.position = 'none'))
    leg <- get_legend(deconv_heatmap_list$heatmap + guides(fill = guide_colourbar(...)))

    if('genes' %in% show_dendros) {
        aligned_genes <- align_plots(deconv_heatmap_list$heatmap + theme(legend.position = 'none'), deconv_heatmap_list$dendro_genes, align = 'h')
        grob_dend_genes <- aligned_genes[[2]]
        grob_htmp$heights[1] <- padding
        grob_dend_genes$heights[1] <- padding
    }

    if('analyses' %in% show_dendros) {
        aligned_analyses <- align_plots(
            deconv_heatmap_list$dendro_analyses,
            deconv_heatmap_list$heatmap + theme(legend.position = 'none'),
            axis = 'rl',
            align = 'v'
        )
        grob_dend_analyses <- aligned_analyses[[1]]
        grob_htmp$widths[9] <- padding
        grob_dend_analyses$widths[9] <- padding
    }


    if(length(show_dendros) == 2) {
        plot_grid(grob_dend_analyses, leg, grob_htmp, grob_dend_genes, nrow = 2, ncol = 2, rel_heights = rel_heights, rel_widths = rel_widths)
    } else if(show_dendros == 'genes') {
        plot_grid(grob_htmp, grob_dend_genes, leg, nrow = 1, ncol = 3, rel_widths = rel_widths)
    } else { # Then show_dendros == 'analyses'
        plot_grid(
            plot_grid(grob_dend_analyses, grob_htmp, ncol = 1, nrow = 2, rel_heights = rel_heights),
            leg,
            nrow = 1,
            ncol = 2,
            rel_widths = rel_widths
        )
    }

}





clinical_test <- function(

    expression_data,
    sample_ids_list,
    genes,
    clinical_data,
    clin_var,
    test_x_expr,
    test_y_expr,
    test_fun = wilcox.test,
    amatch_max_dist = 10,
    p_adjust_method = c('BH', 'bonferroni'),
    min_samples = 1,
    score_method = c('average', 'average_z', 'signature_score'),
    scores_filter_fun = NULL,
    genes_for_data_transform = NULL,
    ...

) {

    # <sample_ids_list> should have names!

    # <genes> can either be a character vector or a list of the same length as <sample_ids_list> (with the same names) with each element being a
    # character vector.  If the former, the same genes are used for all the elements of <sample_ids_list>; if the latter, the genes in the n-th
    # element of <genes> are used for the n-the element of <sample_ids_list>.

    # Similarly, <expression_data> and <clinical_data> can either be single data tables or lists of data tables (with the same length names as
    # <sample_ids_list>), depending on whether you want to use the same table for all cancer types or a different one for each.  If <sample_ids_list>
    # is not given and <expression_data> is a list, the sample IDs will be taken from the ids in <expression_data>.

    # Each of <test_x_expr> and <test_y_expr> can be either a list or a single expression.  If the latter, the same expression is used for every
    # clinical variable; if the former, the n-th expression will be used for the n-th variable.

    # By default, matches for the elements of <clin_var> are searched for using amatch() from the stringdist package, with maxDist equal to 10.  The
    # <amatch_max_dist> argument is passed to the maxDist parameter of amatch(), so can be used to specify how closely you need the variable names in
    # the data to match those in <clin_var>.  If exact matches are required, set <amatch_max_dist> to NULL.

    # <score_method> is the choice of method for assigning scores to each sample based on <genes>.  There are three options: 'average' simply
    # calculates the average expression of the genes in <genes> for each sample; 'average_z' calculates the average Z score of the genes in <genes>;
    # and 'signature_score' uses the signature_score() function, in which case the '...' argument may be used for additional arguments to this
    # function, such as n and nbin.  If <scores_filter_fun> is not NULL, the genes in <expression_data> will be filtered by this function before
    # calculating scores.

    if(missing(sample_ids_list) & !is.data.frame(expression_data)) {sample_ids_list <- lapply(expression_data, `[[`, 'id')}
    p_adjust_method <- match.arg(p_adjust_method)
    score_method <- match.arg(score_method)
    test_fun <- match.fun(test_fun)

    clinical_analyses <- sapply(

        names(sample_ids_list),

        function(ct) {

            cat(paste0(ct, '...\n'))

            if(typeof(genes) == 'list') {
                genes_vec <- genes[[ct]]
            } else {
                genes_vec <- genes
            }

            if(!is.data.frame(expression_data)) {
                expression_data_table <- expression_data[[ct]]
            } else {
                expression_data_table <- expression_data
            }

            if(!is.data.frame(clinical_data)) {
                clinical_data_table <- clinical_data[[ct]]
            } else {
                clinical_data_table <- clinical_data
            }

            setkey(expression_data_table, id)
            setkey(clinical_data_table, id)

            sample_ids <- sample_ids_list[[ct]]

            # The following is for the (optional) data transform by calculating residuals.  I do this before refining sample_ids because I want the
            # transformation to be the same as in deconv_data.
            if(!is.null(genes_for_data_transform)) {
                if(typeof(genes_for_data_transform) == 'list') {
                    genes_transform <- genes_for_data_transform[[ct]]
                } else {
                    genes_transform <- genes_for_data_transform
                }
                expression_data_table <- expression_data_table[sample_ids][, row_sums := rowSums(.SD), .SDcols = genes_transform]
                expression_data_table <- as.data.table(
                    sapply(
                        genes_vec,
                        function(g) lm(formula(paste0('`', g, '` ~ row_sums')), data = expression_data_table)$residuals,
                        simplify = FALSE,
                        USE.NAMES = TRUE
                    )
                )[, id := sample_ids]
                setcolorder(expression_data_table, 'id')
                setkey(expression_data_table, id)
            }

            # Match IDs from clinical data and expression data on the level of patients, i.e. first extract the first 3 components of each element in
            # sample_ids and then look for matches for these patient IDs in clinical_data_table$id.  We exclude non-unique matches.

            all_ids <- data.table(
                expression_data_id = sample_ids,
                patient_id = apply(stringr::str_split_fixed(sample_ids, '\\.', 4)[, 1:3], 1, paste, collapse = '.')
            )[
                ,
                clinical_data_id := sapply(
                    patient_id,
                    function(x) {
                        matches <- clinical_data_table[, id[grep(paste0('^', x), id)]]
                        if(length(matches) == 1) {return(matches)} else {return(NA)}
                    }
                )
            ][!(clinical_data_id %in% names(table(clinical_data_id))[table(clinical_data_id) > 1])]

            # The following two lines allow consistency with the rest of the code for this function, which I wrote before the addition of all_ids.
            patient_ids <- all_ids$clinical_data_id
            sample_ids <- all_ids$expression_data_id

            # Look for a match for clin_var, with specified max. distance, among the variables that are not all NA or the empty string for this cancer
            # type/subtype, and return NA if we don't find one:
            non_na_names <- colnames(clinical_data_table[patient_ids, sapply(.SD, function(x) switch((sum(!is.na(x) & x != '') == 0) + 1, x, NULL))])

            # The following returns a list with names being the non-NA variable name matches and each element being a list consisting of the
            # corresponding x and y Wilcoxon test expressions along with the input variable name.

            ct_clin_vars <- unlist(
                lapply(
                    clin_var,
                    function(v) {
                        if(typeof(test_x_expr) == 'list') {
                            x_expr <- test_x_expr[[which(clin_var == v)]]
                        } else {
                            x_expr <- test_x_expr
                        }
                        if(typeof(test_y_expr) == 'list') {
                            y_expr <- test_y_expr[[which(clin_var == v)]]
                        } else {
                            y_expr <- test_y_expr
                        }
                        if(is.null(amatch_max_dist)) {
                            if(v %in% non_na_names) {
                                return(
                                    setNames(
                                        list( # List whose one element is a list with 3 elements
                                            list(variable_name = v, x_expr = x_expr, y_expr = y_expr)
                                        ),
                                        v # Name for the single element of the 'outer' list
                                    )
                                )
                            } else {
                                cat("Cannot find", v, "in data.\n")
                            }
                        } else {
                            var_match <- non_na_names[stringdist::amatch(v, non_na_names, maxDist = amatch_max_dist)]
                            if(is.na(var_match)) {
                                cat("No variable match for ", v, ".\n", sep = '')
                            } else {
                                cat("Found match for ", v, ": ", var_match, ".\n", sep = '')
                                return(
                                    setNames(
                                        list( # List whose one element is a list with 3 elements
                                            list(variable_name = v, x_expr = x_expr, y_expr = y_expr)
                                        ),
                                        var_match # Name for the single element of the 'outer' list
                                    )
                                )
                            }
                        }
                    }
                ),
                recursive = FALSE
            )

            # Calculate scores by chosen method:

            if(is.null(scores_filter_fun)) {
                if(score_method == 'average') {
                    ct_scores <- setNames(rowMeans(expression_data_table[sample_ids, ..genes_vec]), sample_ids)
                } else if(score_method == 'average_z') {
                    ct_scores <- setNames(rowMeans(apply(expression_data_table[sample_ids, ..genes_vec], 2, scale)), sample_ids)
                } else {
                    ct_scores <- signature_score(set_colnames(t(expression_data_table[sample_ids, -'id']), sample_ids), genes_vec, ...)
                }
            } else {
                filtered_genes <- expression_data_table[sample_ids, names(.SD)[apply(.SD, 2, scores_filter_fun)], .SDcols = -'id']
                if(score_method == 'average') {
                    ct_scores <- setNames(
                        rowMeans(expression_data_table[sample_ids, genes_vec[genes_vec %in% filtered_genes], with = FALSE]),
                        sample_ids
                    )
                } else if(score_method == 'average_z') {
                    ct_scores <- setNames(
                        rowMeans(apply(expression_data_table[sample_ids, genes_vec[genes_vec %in% filtered_genes], with = FALSE], 2, scale)),
                        sample_ids
                    )
                } else {
                    ct_scores <- signature_score(
                        set_colnames(t(expression_data_table[sample_ids, ..filtered_genes]), sample_ids),
                        genes_vec[genes_vec %in% filtered_genes],
                        ...
                    )
                }
            }

            test_results <- sapply(
                names(ct_clin_vars),
                function(v) {

                    # Check that we have enough observations for the Wilcoxon test (by default, we only ignore cases with no observations):
                    criterion <- clinical_data_table[patient_ids, .(variable = get(v))][
                        !is.na(variable) & variable != '',
                        c(sum(eval(ct_clin_vars[[v]]$x_expr)) < min_samples, sum(eval(ct_clin_vars[[v]]$y_expr)) < min_samples)
                    ]

                    if(sum(criterion) > 0) {

                        cat("Not enough ", paste(c('x', 'y')[criterion], collapse = ' or '), " observations for ", v, ".\n", sep = '')

                    } else {

                        # Restrict sample and patient IDs to those where the clinical variable is not NA:
                        sids <- sample_ids[clinical_data_table[patient_ids, !is.na(get(v)) & get(v) != '']]
                        pids <- patient_ids[clinical_data_table[patient_ids, !is.na(get(v)) & get(v) != '']]

                        return(
                            list(
                                dt = clinical_data_table[pids, .(id = id, variable = get(v), score = ct_scores[sids])][
                                    ,
                                    c(
                                        test_name = ct,
                                        variable_name = ct_clin_vars[[v]]$variable_name,
                                        variable_match = v,
                                        setNames(
                                            test_fun(
                                                .SD[eval(ct_clin_vars[[v]]$x_expr), score],
                                                .SD[eval(ct_clin_vars[[v]]$y_expr), score]
                                            )[c('statistic', 'p.value')],
                                            c('stat', 'pval')
                                        ),
                                        ratio = .SD[eval(ct_clin_vars[[v]]$x_expr), mean(score)]/.SD[eval(ct_clin_vars[[v]]$y_expr), mean(score)],
                                        diff = .SD[eval(ct_clin_vars[[v]]$x_expr), mean(score)] - .SD[eval(ct_clin_vars[[v]]$y_expr), mean(score)]
                                    )
                                ],
                                sample_ids = sids,
                                patient_ids = pids
                            )
                        )

                    }

                },
                simplify = FALSE,
                USE.NAMES = TRUE
            )

            cat('Done!\n')

            test_results

        },

        simplify = FALSE,
        USE.NAMES = TRUE

    )

    # Bind together into a single data table:
    clinical_analyses_table <- rbindlist(unlist(lapply(clinical_analyses, function(x) lapply(x, `[[`, 'dt')), recursive = FALSE))

    if(nrow(clinical_analyses_table) > 0) {
        clinical_analyses_table[, c('test_id', 'pval_adj') := .(.I, p.adjust(pval, p_adjust_method))]
        setcolorder(clinical_analyses_table, c('test_id', 'test_name', 'variable_name', 'variable_match', 'stat', 'pval', 'pval_adj'))
        return(
            list(
                data = clinical_analyses_table,
                sample_ids = sapply(
                    clinical_analyses,
                    function(x) sapply(x, `[[`, 'sample_ids', simplify = FALSE, USE.NAMES = TRUE),
                    simplify = FALSE,
                    USE.NAMES = TRUE
                ),
                patient_ids = sapply(
                    clinical_analyses,
                    function(x) sapply(x, `[[`, 'patient_ids', simplify = FALSE, USE.NAMES = TRUE),
                    simplify = FALSE,
                    USE.NAMES = TRUE
                )
            )
        )
    } else {
        cat('\n')
        return(NULL)
    }

}





adj_signif_threshold <- function(signif_threshold, p_adjust_method = c('BH', 'bonferroni'), test_data = NULL) {
    p_adjust_method = match.arg(p_adjust_method)
    if(p_adjust_method == 'BH') {
        if(is.null(test_data)) {stop("Please provide <test_data> if using method 'BH'.")}
        test_data[
            order(pval),
            switch((sum(pval <= signif_threshold*.I/.N) == 0) + 1, signif_threshold*max(which(pval <= signif_threshold*.I/.N))/.N, NA)
        ]
    } else { # Then p_adjust_method == 'bonferroni'
        signif_threshold/nrow(test_data)
    }
}





clinical_test_heatmap <- function(

    test_data,
    x_var = 'variable',
    y_var = 'test',
    fill_var = 'fill',
    hclust_method = 'average',
    x_factor_levels = NULL,
    y_factor_levels = NULL,
    colours = rev(colorRampPalette(RColorBrewer::brewer.pal(11, "PiYG"))(50)),
    limits = c(-2, 2),
    breaks = c(-2, -1, 0, 1, 2),
    labels = c('-2' = '\u2264 -2', '-1' = '-1', '0' = '0', '1' = '1', '2' = '\u2265 2'),
    na_value = '#F7F7F7',
    x_lab = 'Prognostic feature',
    y_lab = 'Cancer type',
    legend_title = latex2exp::TeX('-log_{10}(p value) $\\times$ sign(fold change)'),
    plot_title = NULL,
    grid_lines = NA,
    ...

) {

    # If <hclust_method> is NULL, no clustering will be done.  Use <x_factor_levels> and <y_factor_levels> to manually specify the order in which you
    # want the variable values to appear on the axes.  These arguments will be ignored if <hclust_method> is NULL. The '...' argument is for extra
    # arguments for theme().

    # Use keys to fill in absent entries with NAs:
    setkeyv(test_data, c(x_var, y_var))

    # Output depends on whether we want to do clustering:

    if(!is.null(hclust_method)) {

        clust_data <- dcast(
            test_data,
            get(y_var) ~ get(x_var), # We end up with a column called 'y_var'
            value.var = fill_var,
            fill = 0
        )

        hclust_x <- hclust(dist(t(clust_data[, -'y_var'])), method = hclust_method)
        hclust_y <- hclust(dist(clust_data[, magrittr::set_rownames(as.matrix(.SD), y_var), .SD = -'y_var']), method = hclust_method)

        htmp <- ggplot(
            test_data[expand.grid(unique(get(x_var)), unique(get(y_var)))],
            aes(
                x = factor(get(x_var), levels = hclust_x$labels[hclust_x$order]),
                y = factor(get(y_var), levels = hclust_y$labels[hclust_y$order]),
                fill = get(fill_var)
            )
        ) +
            geom_tile(colour = grid_lines) +
            scale_x_discrete(expand = c(0, 0)) +
            scale_y_discrete(expand = c(0, 0)) +
            scale_fill_gradientn(
                limits = limits,
                breaks = breaks,
                labels = labels,
                colours = colours,
                na.value = na_value,
                oob = scales::squish
            ) +
            theme(panel.border = element_rect(size = 0.5, fill = NA), axis.text.x = element_text(angle = 55, hjust = 1), ...) +
            labs(x = x_lab, y = y_lab, fill = legend_title, title = plot_title)

        list(heatmap = htmp, hclust_x = hclust_x, hclust_y = hclust_y)

    } else {

        ggplot(
            test_data[expand.grid(unique(get(x_var)), unique(get(y_var)))],
            aes(
                x = switch(is.null(x_factor_levels) + 1, factor(get(x_var), levels = x_factor_levels), get(x_var)),
                y = switch(is.null(y_factor_levels) + 1, factor(get(y_var), levels = y_factor_levels), get(y_var)),
                fill = get(fill_var)
            )
        ) +
            geom_tile(colour = grid_lines) +
            scale_x_discrete(expand = c(0, 0)) +
            scale_y_discrete(expand = c(0, 0)) +
            scale_fill_gradientn(
                limits = limits,
                breaks = breaks,
                labels = labels,
                colours = colours,
                na.value = na_value,
                oob = scales::squish
            ) +
            theme(panel.border = element_rect(size = 0.5, fill = NA), axis.text.x = element_text(angle = 55, hjust = 1), ...) +
            labs(x = x_lab, y = y_lab, fill = legend_title, title = plot_title)

    }

}





deconv_plot <- function(

    deconvs_list,
    n_row = 1,
    n_col = length(deconvs_list),
    plots_rel_heights = c(title = 2, purity_bar = 1, ccle_bar = 1, extra_bar = 1, heatmap = 15, axis_labels = 5),
    rows_rel_heights = rep(1, n_row),
    left_plot_width = 1,
    align = c('centre', 'center', 'left'),
    legends = TRUE,
    legends_arrange = switch((length(deconvs_list) == n_row*n_col) + 1, 'horizontal', 'vertical'),
    legends_rel_size = switch((legends_arrange == 'vertical') + 1, 1, 0.2),
    legends_space = 1,
    ...

) {

    # The '...' argument is for guide_colourbar(), and can be used to change the format of the legends.

    # <deconvs_list> should be a list with each element of the form returned by deconvolve_emt_caf().

    # <left_plot_width> specifies the desired width of the leftmost plot in each row as a fraction of the normal width of the other plots.  It is 1 by
    # default, but it can be increased manually if, for example, y axis text on the leftmost plots makes these plots slightly wider than the others.

    # The <legends> argument can be used to specify whether you want the legends or not.  It is TRUE by default.

    # The legends are arranged in a list as follows: a blank plot first and last, and between these are the legends alternating with blank plots.  So
    # every other element of the list is a blank plot, and both the first and last elements of the list are blank plots.  This is so that you can
    # adjust the spacing between, before and after the legends by adjusting the relative sizes of the blank plots.  To do this manually, it helps to
    # know how long the legend list will be.  There will always be one more blank plot than legends, so the total length of the legend list will be
    # 2*(number of legends) + 1.

    # The <align> argument refers to the alignment of the plots (and legends) on the bottom row, if indeed there is space to spare.  If <align> is
    # 'centre' (or 'center'), blank plots are inserted at the left and right ends to provide the padding to make the alignment central.  The plots
    # will be the same width as those on the upper row(s).  The user can specify the width to be taken up by the legends via the <legends_space>
    # argument.  The number passed to this is interpreted as a fraction of the width of one plot.  It is 1 by default, meaning the legends take up the
    # same width as one plot.  The widths of the blank padding plots are chosen automatically to take up all remaining space.  If <align> is 'left',
    # no blank plots are inserted, <legend_space> is ignored, and the legend space takes up all the space not occupied by the plots.  In this case,
    # the <legends_rel_size> argument can be used to change the space occupied by the legends.

    if(length(rows_rel_heights) != n_row) {stop('<rows_rel_heights> does not have length equal to <n_row>.')}

    align <- match.arg(align)
    legends_arrange <- match.arg(legends_arrange, choices = c('horizontal', 'vertical'))

    plot_lists <- lapply(
        1:length(deconvs_list),
        function(i) {

            if(i %in% seq(from = 1, to = length(deconvs_list), by = n_col)) {

                bars <- lapply(
                    deconvs_list[[i]]$plots[
                        c('purity_bar', 'ccle_bar', 'extra_bar')[c('purity_bar', 'ccle_bar', 'extra_bar') %in% names(deconvs_list[[i]]$plots)]
                    ],
                    function(x) {x + theme(legend.position = 'none')}
                )

                pl <- c(
                    list(
                        title = ggplot() +
                            labs(caption = deconvs_list[[i]]$plots$heatmap$labels$title) +
                            theme(panel.background = element_rect(fill = 'white'), plot.caption = element_text(size = 14, hjust = 0))
                    ),
                    bars,
                    list(
                        heatmap = deconvs_list[[i]]$plots$heatmap + theme(
                            legend.position = 'none',
                            plot.title = element_blank(),
                            plot.margin = unit(c(0, 5.5, switch(('axis_labels' %in% names(deconvs_list[[i]]$plots)) + 1, 5.5, 0), 5.5), 'pt')
                        )
                    )
                )

                if('axis_labels' %in% names(deconvs_list[[i]]$plots)) {
                    pl <- c(pl, list(axis_labels = deconvs_list[[i]]$plots$axis_labels))
                }

            } else {

                bars <- lapply(
                    deconvs_list[[i]]$plots[
                        c('purity_bar', 'ccle_bar', 'extra_bar')[c('purity_bar', 'ccle_bar', 'extra_bar') %in% names(deconvs_list[[i]]$plots)]
                    ],
                    function(x) {x + theme(legend.position = 'none', axis.title.y = element_blank())}
                )

                pl <- c(
                    list(
                        title = ggplot() +
                            labs(caption = deconvs_list[[i]]$plots$heatmap$labels$title) +
                            theme(panel.background = element_rect(fill = 'white'), plot.caption = element_text(size = 14, hjust = 0))
                    ),
                    bars,
                    list(
                        heatmap = deconvs_list[[i]]$plots$heatmap + theme(
                            legend.position = 'none',
                            title = element_blank(),
                            plot.margin = unit(c(0, 5.5, switch(('axis_labels' %in% names(deconvs_list[[i]]$plots)) + 1, 5.5, 0), 5.5), 'pt')
                        )
                    )
                )

                if('axis_labels' %in% names(deconvs_list[[i]]$plots)) {
                    pl <- c(pl, list(axis_labels = deconvs_list[[i]]$plots$axis_labels))
                }

            }

            pl

        }
    )

    # Optionally make legends plot:

    if(legends) {

        # Make list of legend plots with alternating legends and blank plots:
        legend_list <- lapply(
            deconvs_list[[1]]$plots[
                c('purity_bar', 'ccle_bar', 'extra_bar')[c('purity_bar', 'ccle_bar', 'extra_bar') %in% names(deconvs_list[[1]]$plots)]
            ],
            function(g) {get_legend(g + guides(fill = guide_colourbar(...)))}
        )
        legend_list <- unlist(lapply(legend_list, list, blank_plot()), recursive = FALSE)
        legend_list <- c(
            list(blank_plot()),
            legend_list,
            list(get_legend(deconvs_list[[1]]$plots$heatmap + guides(fill = guide_colourbar(...))), blank_plot())
        )

        if(length(legends_rel_size) == 1) {legends_rel_size <- rep(legends_rel_size, length(legend_list))}

        # Make the legend plot:
        if(legends_arrange == 'vertical') {
            legend_plot <- plot_grid(plotlist = legend_list, nrow = length(legend_list), ncol = 1, rel_heights = legends_rel_size)
        } else {
            legend_plot <- plot_grid(plotlist = legend_list, nrow = 1, ncol = length(legend_list), rel_widths = legends_rel_size )
        }

    } else {

        legend_plot <- blank_plot()

    }

    # Make the whole plot:

    if(length(deconvs_list) == n_row*n_col) {

        if(legends) {
            plot_grid(
                plot_grid(
                    plotlist = lapply(
                        plot_lists,
                        function(l) plot_grid(plotlist = l, rel_heights = plots_rel_heights[names(l)], ncol = 1, nrow = length(l), align = 'v')
                    ),
                    nrow = n_row,
                    ncol = n_col,
                    rel_widths = c(left_plot_width, rep(1, n_col - 1)),
                    rel_heights = rows_rel_heights,
                    align = 'hv'
                ),
                legend_plot,
                nrow = 1,
                ncol = 2,
				rel_widths = c(n_col, legends_space)
            )
        } else {
            plot_grid(
                plotlist = lapply(
                    plot_lists,
                    function(l) plot_grid(plotlist = l, rel_heights = plots_rel_heights[names(l)], ncol = 1, nrow = length(l), align = 'v')
                ),
                nrow = n_row,
                ncol = n_col,
                rel_widths = c(left_plot_width, rep(1, n_col - 1)),
                rel_heights = rows_rel_heights,
                align = 'hv'
            )
        }

    } else if(length(deconvs_list) == n_row*n_col - 1) {

        if(length(plot_lists) %/% n_col == 0) {
            plot_grid(
                plotlist = c(
                    lapply(
                        plot_lists,
                        function(l) plot_grid(plotlist = l, rel_heights = plots_rel_heights[names(l)], ncol = 1, nrow = length(l), align = 'v')
                    ),
                    list(legend_plot)
                ),
                nrow = 1,
                ncol = n_col,
                rel_widths = c(left_plot_width, rep(1, n_col - 1))
            )
        } else {
            plot_grid(
                plot_grid(
                    plotlist = lapply(
                        plot_lists[1:((length(plot_lists) %/% n_col)*n_col)],
                        function(l) plot_grid(plotlist = l, rel_heights = plots_rel_heights[names(l)], ncol = 1, nrow = length(l), align = 'v')
                    ),
                    nrow = length(plot_lists) %/% n_col,
                    ncol = n_col,
                    rel_widths = c(left_plot_width, rep(1, n_col - 1)),
                    rel_heights = rows_rel_heights[-n_row],
                    align = 'hv'
                ),
                plot_grid( # Bottom row
                    plotlist = c(
                        lapply(
                            plot_lists[(length(plot_lists) - (length(plot_lists) %% n_col) + 1):length(plot_lists)],
                            function(l) plot_grid(plotlist = l, rel_heights = plots_rel_heights[names(l)], ncol = 1, nrow = length(l), align = 'v')
                        ),
                        list(legend_plot)
                    ),
                    nrow = 1,
                    ncol = n_col,
                    rel_widths = c(left_plot_width, rep(1, n_col - 1))
                ),
                nrow = 2,
                ncol = 1,
                rel_heights = c(sum(rows_rel_heights[-n_row]), rows_rel_heights[n_row])
            )
        }

    } else {

        n_bottom_row <- length(plot_lists) %% n_col

        if(align %in% c('centre', 'center')) {
            plot_grid(
                plot_grid(
                    plotlist = lapply(
                        plot_lists[1:((length(plot_lists) %/% n_col)*n_col)],
                        function(l) plot_grid(plotlist = l, rel_heights = plots_rel_heights[names(l)], ncol = 1, nrow = length(l), align = 'v')
                    ),
                    nrow = length(plot_lists) %/% n_col,
                    ncol = n_col,
                    rel_widths = c(left_plot_width, rep(1, n_col - 1)),
                    rel_heights = rows_rel_heights[-n_row],
                    align = 'hv'
                ),
                plot_grid( # Bottom row
                    plotlist = c(
                        list(blank_plot()),
                        lapply(
                            plot_lists[(length(plot_lists) - n_bottom_row + 1):length(plot_lists)],
                            function(l) plot_grid(plotlist = l, rel_heights = plots_rel_heights[names(l)], ncol = 1, nrow = length(l), align = 'v')
                        ),
                        list(legend_plot, blank_plot())
                    ),
                    nrow = 1,
                    rel_widths = c(
                        (n_col - n_bottom_row - legends_space)/2,
                        left_plot_width,
                        rep(1, n_bottom_row - 1),
                        legends_space,
                        (n_col - n_bottom_row - legends_space)/2
                    )
                ),
                nrow = 2,
                ncol = 1,
                rel_heights = c(sum(rows_rel_heights[-n_row]), rows_rel_heights[n_row])
            )
        } else {
            plot_grid(
                plot_grid(
                    plotlist = lapply(
                        plot_lists[1:((length(plot_lists) %/% n_col)*n_col)],
                        function(l) plot_grid(plotlist = l, rel_heights = plots_rel_heights[names(l)], ncol = 1, nrow = length(l), align = 'v')
                    ),
                    nrow = length(plot_lists) %/% n_col,
                    ncol = n_col,
                    rel_widths = c(left_plot_width, rep(1, n_col - 1)),
                    rel_heights = rows_rel_heights[-n_row],
                    align = 'hv'
                ),
                plot_grid( # Bottom row
                    plotlist = c(
                        lapply(
                            plot_lists[(length(plot_lists) - n_bottom_row + 1):length(plot_lists)],
                            function(l) plot_grid(plotlist = l, rel_heights = plots_rel_heights[names(l)], ncol = 1, nrow = length(l), align = 'v')
                        ),
                        list(legend_plot)
                    ),
                    nrow = 1,
                    rel_widths = c(left_plot_width, rep(1, n_bottom_row - 1), n_col - n_bottom_row)
                ),
                nrow = 2,
                ncol = 1,
                rel_heights = c(sum(rows_rel_heights[-n_row]), rows_rel_heights[n_row])
            )
        }

    }

}





# Function to make brackets with which to annotate the deconv results as 'pEMT' and 'CAF' (plotted separately, to be aligned with the deconv plot):

pemt_caf_brackets <- function(
	edge = c('top', 'right', 'bottom', 'left'),
    pemt_bracket_max = 20,
    caf_bracket_min = 80,
    brackets_just = 0,
    tips_end = ifelse(edge %in% c('top', 'right'), 1, -1),
    labels = c('pEMT', 'CAF'),
	pemt_label_hjust = 0.5,
	caf_label_hjust = 0.5,
    labels_vjust = ifelse(edge %in% c('top', 'right'), 1.5, -1.5),
	labels_angle = 0,
    pemt_bracket_colour = 'black',
    caf_bracket_colour = 'black',
    pemt_label_colour = 'black',
    caf_label_colour = 'black',
	plot_margin = c(0, 0, 0, 0)
) {

	edge <- match.arg(edge)

	if(edge == 'top') {

		ggplot(data.frame(x = 1:100, y = seq(-1, 1, length.out = 100))) +
			theme(
				panel.background = element_blank(),
				axis.text = element_blank(),
				axis.title = element_blank(),
				axis.ticks = element_blank(),
				axis.ticks.length = unit(0, 'pt'),
				plot.margin = unit(plot_margin, 'pt')
			) +
			scale_x_continuous(expand = c(0, 0)) +
			scale_y_continuous(expand = c(0, 0), limits = c(-1, 1)) +
			geom_segment(aes(x = 0, xend = pemt_bracket_max, y = brackets_just, yend = brackets_just), colour = pemt_bracket_colour) +
			geom_segment(aes(x = caf_bracket_min, xend = 100, y = brackets_just, yend = brackets_just), colour = caf_bracket_colour) +
			geom_segment(aes(x = 0, xend = 0, y = brackets_just, yend = tips_end), colour = pemt_bracket_colour) +
			geom_segment(aes(x = pemt_bracket_max, xend = pemt_bracket_max, y = brackets_just, yend = tips_end), colour = pemt_bracket_colour) +
			geom_segment(aes(x = caf_bracket_min, xend = caf_bracket_min, y = brackets_just, yend = tips_end), colour = caf_bracket_colour) +
			geom_segment(aes(x = 100, xend = 100, y = brackets_just, yend = tips_end), colour = caf_bracket_colour) +
			annotate(
				geom = 'text',
				x = pemt_bracket_max/2,
				y = brackets_just,
				label = labels[1],
				hjust = pemt_label_hjust,
				vjust = labels_vjust,
				colour = pemt_label_colour,
				angle = labels_angle
			) +
			annotate(
				geom = 'text',
				x = caf_bracket_min + (100 - caf_bracket_min)/2,
				y = brackets_just,
				label = labels[2],
				hjust = caf_label_hjust,
				vjust = labels_vjust,
				colour = caf_label_colour,
				angle = labels_angle
			)

	} else if(edge == 'right') {

		ggplot(data.frame(x = seq(-1, 1, length.out = 100), y = 1:100)) +
			theme(
				panel.background = element_blank(),
				axis.text = element_blank(),
				axis.title = element_blank(),
				axis.ticks = element_blank(),
				axis.ticks.length = unit(0, 'pt'),
				plot.margin = unit(plot_margin, 'pt')
			) +
			scale_x_continuous(expand = c(0, 0), limits = c(-1, 1)) +
			scale_y_continuous(expand = c(0, 0)) +
			geom_segment(aes(x = brackets_just, xend = brackets_just, y = 0, yend = pemt_bracket_max), colour = pemt_bracket_colour) +
			geom_segment(aes(x = brackets_just, xend = brackets_just, y = caf_bracket_min, yend = 100), colour = caf_bracket_colour) +
			geom_segment(aes(x = brackets_just, xend = tips_end, y = 0, yend = 0), colour = pemt_bracket_colour) +
			geom_segment(aes(x = brackets_just, xend = tips_end, y = pemt_bracket_max, yend = pemt_bracket_max), colour = pemt_bracket_colour) +
			geom_segment(aes(x = brackets_just, xend = tips_end, y = caf_bracket_min, yend = caf_bracket_min), colour = caf_bracket_colour) +
			geom_segment(aes(x = brackets_just, xend = tips_end, y = 100, yend = 100), colour = caf_bracket_colour) +
			annotate(
				geom = 'text',
				x = brackets_just,
				y = pemt_bracket_max/2,
				label = labels[1],
				hjust = pemt_label_hjust,
				vjust = labels_vjust,
				colour = pemt_label_colour,
				angle = labels_angle
			) +
			annotate(
				geom = 'text',
				x = brackets_just,
				y = caf_bracket_min + (100 - caf_bracket_min)/2,
				label = labels[2],
				hjust = caf_label_hjust,
				vjust = labels_vjust,
				colour = caf_label_colour,
				angle = labels_angle
			)

	} else if(edge == 'bottom') {

		ggplot(data.frame(x = 1:100, y = seq(-1, 1, length.out = 100))) +
			theme(
				panel.background = element_blank(),
				axis.text = element_blank(),
				axis.title = element_blank(),
				axis.ticks = element_blank(),
				axis.ticks.length = unit(0, 'pt'),
				plot.margin = unit(plot_margin, 'pt')
			) +
			scale_x_continuous(expand = c(0, 0)) +
			scale_y_continuous(expand = c(0, 0), limits = c(-1, 1)) +
			geom_segment(aes(x = 0, xend = pemt_bracket_max, y = brackets_just, yend = brackets_just), colour = pemt_bracket_colour) +
			geom_segment(aes(x = caf_bracket_min, xend = 100, y = brackets_just, yend = brackets_just), colour = caf_bracket_colour) +
			geom_segment(aes(x = 0, xend = 0, y = brackets_just, yend = tips_end), colour = pemt_bracket_colour) +
			geom_segment(aes(x = pemt_bracket_max, xend = pemt_bracket_max, y = brackets_just, yend = tips_end), colour = pemt_bracket_colour) +
			geom_segment(aes(x = caf_bracket_min, xend = caf_bracket_min, y = brackets_just, yend = tips_end), colour = caf_bracket_colour) +
			geom_segment(aes(x = 100, xend = 100, y = brackets_just, yend = tips_end), colour = caf_bracket_colour) +
			annotate(
				geom = 'text',
				x = pemt_bracket_max/2,
				y = brackets_just,
				label = labels[1],
				hjust = pemt_label_hjust,
				vjust = labels_vjust,
				colour = pemt_label_colour,
				angle = labels_angle
			) +
			annotate(
				geom = 'text',
				x = brackets_just,
				y = caf_bracket_min + (100 - caf_bracket_min)/2,
				label = labels[2],
				hjust = caf_label_hjust,
				vjust = labels_vjust,
				colour = caf_label_colour,
				angle = labels_angle
			)

	} else if(edge == 'top') {

		ggplot(data.frame(x = seq(-1, 1, length.out = 100), y = 1:100)) +
			theme(
				panel.background = element_blank(),
				axis.text = element_blank(),
				axis.title = element_blank(),
				axis.ticks = element_blank(),
				axis.ticks.length = unit(0, 'pt'),
				plot.margin = unit(plot_margin, 'pt')
			) +
			scale_x_continuous(expand = c(0, 0), limits = c(-1, 1)) +
			scale_y_continuous(expand = c(0, 0)) +
			geom_segment(aes(x = brackets_just, xend = brackets_just, y = 0, yend = pemt_bracket_max), colour = pemt_bracket_colour) +
			geom_segment(aes(x = brackets_just, xend = brackets_just, y = caf_bracket_min, yend = 100), colour = caf_bracket_colour) +
			geom_segment(aes(x = brackets_just, xend = tips_end, y = 0, yend = 0), colour = pemt_bracket_colour) +
			geom_segment(aes(x = brackets_just, xend = tips_end, y = pemt_bracket_max, yend = pemt_bracket_max), colour = pemt_bracket_colour) +
			geom_segment(aes(x = brackets_just, xend = tips_end, y = caf_bracket_min, yend = caf_bracket_min), colour = caf_bracket_colour) +
			geom_segment(aes(x = brackets_just, xend = tips_end, y = 100, yend = 100), colour = caf_bracket_colour) +
			annotate(
				geom = 'text',
				x = brackets_just,
				y = pemt_bracket_max/2,
				label = labels[1],
				hjust = pemt_label_hjust,
				vjust = labels_vjust,
				colour = pemt_label_colour,
				angle = labels_angle
			) +
			annotate(
				geom = 'text',
				x = brackets_just,
				y = caf_bracket_min + (100 - caf_bracket_min)/2,
				label = labels[2],
				hjust = caf_label_hjust,
				vjust = labels_vjust,
				colour = caf_label_colour,
				angle = labels_angle
			)

	}

}





# Function to filter and reorder the genes in the deconv results based on correlation with head and tail genes:

deconv_reorder <- function(deconv_obj, to_keep = c('SNAI1', 'SNAI2', 'TWIST1', 'VIM', 'ZEB1', 'ZEB2')) {

	head_genes <- with(deconv_obj, head(genes_filtered[ordering], 20))
	tail_genes <- with(deconv_obj, tail(genes_filtered[ordering], 20))

	new_scores <- sapply(
		colnames(deconv_obj$cor_mat),
		function(g) {
			meanvec <- c(mean(deconv_obj$cor_mat[, g][head_genes]), mean(deconv_obj$cor_mat[, g][tail_genes]))
			ifelse(
				sign(meanvec[1]) == sign(meanvec[2]),
				ifelse(g %in% to_keep, return(0), return(NA)),
				return(max(meanvec)*c(1, -1)[which.max(meanvec)])
			)
		}
	)

	new_scores <- new_scores[!is.na(new_scores)]

	new_deconv <- deconv_obj
	new_deconv$new_scores <- new_scores
	new_deconv$cor_mat <- new_deconv$cor_mat[names(new_scores), names(new_scores)]
	new_deconv$genes_filtered <- names(new_scores)
	new_deconv$ordering <- order(-new_scores)
	new_deconv$cor_with_purity <- sapply(new_deconv$cor_with_purity, function(x) x[names(new_scores)], simplify = FALSE, USE.NAMES = TRUE)

	if('ccle_comp_diff' %in% names(new_deconv)) {new_deconv$ccle_comp_diff <- new_deconv$ccle_comp_diff[names(new_scores)]}
	if('extra_data_score' %in% names(new_deconv)) {new_deconv$extra_data_score <- new_deconv$extra_data_score[names(new_scores)]}

	new_deconv

}
