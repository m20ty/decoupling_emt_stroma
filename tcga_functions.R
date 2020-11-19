infer_subtypes <- function(expression_data, subtypes_dt) {

    #expression_data should be in the form returned by prepare_expression_data(), but
    #any transformations of the data to influence the selection of marker genes and
    #the calculation of subtype scores, such as scaling and taking logarithms, should
    #be done beforehand - this function will do no further transformations for these
    #calculations.  I have found that centring, scaling and taking logs improves the
    #classification into subtypes, but there isn't much difference between using raw
    #counts vs TPM.

    #subtypes_dt must be a data table with first column being tumour IDs and
    #second being the predefined subtype assignments.

    #Make sure column names are what I want them to be:
    names(subtypes_dt)[1:2] <- c('tumour_id', 'subtype')

    #I also need to make sure the subtype column is a character column, not a factor (as
    #happens when you read it in with read.xlsx), otherwise the subtype column in
    #expression_data_subtypes (see below) will consist of numbers instead of names.  I'll
    #also make the tumour IDs characters, while I'm at it.
    subtypes_dt$tumour_id <- as.character(subtypes_dt$tumour_id)
    subtypes_dt$subtype <- as.character(subtypes_dt$subtype)



    #Find gene signatures for each of the subtypes:

    expression_data$subtype <- sapply(expression_data$id, function(tumour) {
        if(tumour %in% subtypes_dt$tumour_id) {
            subtypes_dt[tumour_id == tumour, subtype]
        } else {
            ''
        }
    })

    subtype_markers <- lapply(
        unique(subtypes_dt$subtype),
        function(sbtp) {
            names(sort(
                colMeans(expression_data[subtype == sbtp, -c('id', 'subtype')]),
                decreasing = TRUE
            )[1:100])
        }
    )

    names(subtype_markers) <- unique(subtypes_dt$subtype)



    #Score each tumour for how closely it corresponds to each of the subtypes.
    #These scores are just the average expression across all genes in each of the subtype
    #signatures.  Each tumour thus has as many scores as there are subtypes - assign each
    #tumour to the subtype for which it has the highest score.

    subtype_scores <- as.data.table(lapply(subtype_markers, function(markers_vec) {
        rowMeans(expression_data[, markers_vec, with = FALSE])
    }))



    #Create data table of inferred TCGA subtype assignments, along with the
    #assignments given in the TCGA paper supplementary data.  Also calculate the
    #minimum difference between each tumour's scores, for use in filtering.

    expression_data[
        ,
        c('inferred_subtype', 'min_diff') := .(
            names(sapply(1:.N, function(i) {
                which.max(subtype_scores[i])
            })),
            sapply(1:.N, function(i) {
                min(as.numeric(
                    subtype_scores[i, which.max(subtype_scores[i]), with = FALSE]
                ) - as.numeric(
                    subtype_scores[i, -which.max(subtype_scores[i]), with = FALSE]
                ))
            })
        )
    ]



    #Output:

    list(
        call = as.list(match.call()),
        inferred_subtypes_dt = expression_data[, .(
            id,
            subtype,
            inferred_subtype,
            min_diff
        )],
        subtype_markers = subtype_markers,
        subtype_scores = subtype_scores
    )

}





gene_filter_correlations <- function(

    genes_list,
    expression_data,
    cell_type_markers,
    initial_genes = c(
        'SNAI1',
        'SNAI2',
        'TWIST1',
        'VIM',
        'ZEB1',
        'ZEB2'
    ),
    cell_types = c(
        'B_plasma',
        'myocyte',
        'macrophage',
        'endothelial',
        'DC',
        'mast',
        'T',
        'B'
    )

) {

    as.data.table(
        cbind(
            cor(
                expression_data[, ..genes_list],
                expression_data[, ..initial_genes]
            ),
            sapply(
                cell_types,
                function(ct) {
                    rowMeans(
                        cor(
                            expression_data[, ..genes_list],
                            expression_data[
                                ,
                                cell_type_markers[
                                    cell_type == ct &
                                        gene %in% names(expression_data),
                                    gene
                                ],
                                with = FALSE
                            ]
                        )
                    )
                },
                USE.NAMES = TRUE
            )
        ),
        keep.rownames = 'id'
    )[
        ,
        c(
            'max_initial_gene',
            'which_max_initial_gene',
            'max_cell_type',
            'which_max_cell_type'
        ) := .(
            max(abs(unlist(mget(initial_genes)))),
            initial_genes[which.max(abs(unlist(mget(initial_genes))))],
            max(unlist(mget(cell_types))),
            cell_types[which.max(unlist(mget(cell_types)))]
        ),
        by = id
    ]

}





deconvolve_emt_caf <- function(

    expression_data,
    meta_data,
    genes,
    cell_type_markers,
    sample_ids = NULL,
    tcga_cancer_types = NULL,
    subtypes = NULL,
    subtypes_data = NULL,
    ref_for_subtypes = NULL,
    ccle_cancer_type = NULL,
    ccle_data = NULL,
    extra_data_source = NULL,
    extra_data = NULL,
    initial_genes = c(
        'SNAI1',
        'SNAI2',
        'TWIST1',
        'VIM',
        'ZEB1',
        'ZEB2'
    ),
    cell_types = c(
        'B_plasma',
        'myocyte',
        'macrophage',
        'endothelial',
        'DC',
        'mast',
        'T',
        'B'
    ),
    genes_from_tcga_fun = function(x) {
        top_cols_by_fun_cor(x, initial = initial_genes)[ # Using defaults for <FUN> and <threshold>
            1:min(.N, 300), # Take no more than 300 genes
            id
        ]
    },
    # genes_from_tcga_fun = function(x) {
    #     score_top_cols_by_cor(x, initial_genes, 0.7)[1:300, id]
    # },
    # top_genes_n = 300,
    genes_filter_method = c('weighted_sum', 'max_cor'),
    genes_filter_fun = function(x) 1:250,
    # genes_filter_fun = function(x) {1:min(sum(x > 1), 250)}, # Take no more than 250 genes
    gene_weights_fun = median,
    max_initial_gene = TRUE,
    max_cell_type = TRUE,
    cell_type_weights = NULL,
    initial_gene_weights = NULL,
    seed = NULL,
    ordering_fun = function(x) {
        as.vector(
            seriation::seriate(
                as.dist(1 - x), # This is meant for a correlation matrix
                method = 'SPIN_STS'
            )[[1]]
        )
    },
    plot_title = '',
    heatmap_colour_limits = c(-0.6, 0.6),
    heatmap_legend_breaks = c(-0.6, -0.3, 0, 0.3, 0.6),
    heatmap_annotations = initial_genes,
    heatmap_bar_colours = c(
        'darkslategrey',
        'turquoise4',
        'yellow',
        'chocolate1',
        'coral4'
    ),
    purity_fun = function(x) caTools::runmean(x, 30),
    purity_cor_method = c('scale', 'resid'),
    purity_colour_limits = c(-0.3, 0.3),
    purity_legend_breaks = c(-0.3, 0, 0.3),
    ccle_fun = function(x) {caTools::runmean(x, 30)/max(abs(caTools::runmean(x, 30)))},
    ccle_colour_limits = c(-1, 1),
    ccle_legend_breaks = c(-1, 0, 1),
    extra_fun = ccle_fun,
    extra_colour_limits = ccle_colour_limits,
    extra_legend_breaks = ccle_legend_breaks,
    extra_title = '',
    bar_legend_height = unit(10, 'pt'),
    bar_legend_width = unit(10, 'pt')

) {

    genes_filter_method <- match.arg(genes_filter_method)

    if(!is.null(genes_filter_fun)) {
        genes_filter_fun <- match.fun(genes_filter_fun)
    }
    # If we want to filter the genes using an absolute constant, we can use a constant function,
    # e.g. genes_filter_fun <- function(x) 1

    gene_weights_fun <- match.fun(gene_weights_fun)

    ordering_fun <- match.fun(ordering_fun)
    # We could do hierarchical clustering using hclust(1 - x)$order, or something.

    purity_cor_method = match.arg(purity_cor_method)

    purity_fun <- match.fun(purity_fun)
    ccle_fun <- match.fun(ccle_fun)
    extra_fun <- match.fun(extra_fun)





    # Set keys, in case this hasn't already been done:

    setkey(expression_data, id)
    setkey(meta_data, id)

    # Determine which of sample_ids, tcga_cancer_types and subtypes can be used:

    if(

        is.null(sample_ids) &
        is.null(tcga_cancer_types)

    ) {

        stop('Please specify sample IDs or TCGA cancer types.')

    } else if(

        (
            !is.null(sample_ids) &
            !is.null(tcga_cancer_types)
        ) | (
            !is.null(sample_ids) &
            !is.null(subtypes)
        )


    ) {

        warning(
            'If specifying sample_ids, no need to specify tcga_cancer_types or subtypes.  ',
            'Using only sample_ids.'
        )

        tcga_cancer_types <- NULL
        subtypes <- NULL

    } else if(

        !is.null(tcga_cancer_types) &
        !is.null(subtypes) &
        is.null(subtypes_data)
        # In other cases where subtypes_data is not supplied, we'll be using sample_ids anyway
        # (other combinations eliminated above).

    ) {

        stop('Please supply subtypes_data.')

    }

    # Table showing with which of the above 3 if statements each combination was eliminated
    # ('x' means the argument is NULL; '/' means something was supplied to it):

    # sample_ids        x / x x x / / / x x x / / / x /
    # tcga_cancer_types x x / x x / x x / / x / / x / /
    # subtypes          x x x / x x / x / x / / x / / /
    # subtypes_data     x x x x / x x / x / / x / / / /
    # which if stmnt:   1     1 1 2 2   3   1 2 2 2   2





    # If we're not using sample_ids by default, get the appropriate sample IDs using
    # tcga_cancer_types, and subtypes if it is supplied:

    if(

        is.null(sample_ids) # Then we must have !is.null(tcga_cancer_types)

    ) {

        if(

            is.null(subtypes) # Then get all IDs corresponding to these cancer types

        ) {

            sample_ids <- meta_data[
                cancer_type %in% tcga_cancer_types,
                id
            ]

        } else { # Then we must have !is.null(subtypes_data)

            # Set key, in case this hasn't already been done:

            setkey(subtypes_data, id)

            if(is.null(ref_for_subtypes)) {

                sample_ids <- subtypes_data[
                    meta_data[
                        cancer_type %in% tcga_cancer_types,
                        id
                    ]
                ][
                    inferred_subtype %in% subtypes,
                    id
                ]

            } else {

                # Optionally supply reference for subtypes, in case the same subtype name
                # appears in multiple references.

                sample_ids <- subtypes_data[
                    meta_data[
                        cancer_type %in% tcga_cancer_types,
                        id
                    ]
                ][
                    subtype_ref == ref_for_subtypes & inferred_subtype %in% subtypes,
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
        expression_data[
            sample_ids,
            sapply(
                .SD[, -'id'],
                function(x) {
                    switch((sd(x) > 0) + 1, NULL, x)
                }
            )
        ]
    )

    # Subset genes which occur in expression_data and, if provided, ccle_data and in any extra
    # data:

    genes <- genes[genes %in% names(expression_data)]





    # Get EMT markers by correlation with initial genes:

    if(!is.null(genes_from_tcga_fun)) {

        top_genes <- genes_from_tcga_fun(expression_data[, -'id'])

        genes <- unique(
            sort(
                c(
                    genes,
                    top_genes
                )
            )
        )

    }

    # cor_with_initial_genes <- cor(
    #     expression_data[, -c('id', ..initial_genes)],
    #     expression_data[, ..initial_genes]
    # )

    # In the following, high_cor_threshold is used to specify what counts as a "high
    # correlation" - any correlation above this threshold is deemed "high".  This is used in
    # scoring the genes for correlation with the initial genes - specifically, we use the
    # number of "high" correlations with the initial genes.  I tried, using HNSC data,
    # replacing the number of "high" correlations with simply the average correlation with
    # the initial genes.  The gene set I got was more or less the same, though it seemed to
    # be missing some genes that I think are important, like FAP and ACTA2.  So I'll stick
    # with the first method, though perhaps I could make the averaging method available as an
    # argument of the function.

    # top_genes <- rbindlist(
    #     lapply(
    #         rownames(cor_with_initial_genes),
    #         function(g) {
    #
    #             nhc <- sum(
    #                 sapply(
    #                     initial_genes,
    #                     function(s) {
    #                         switch(
    #                             (
    #                                 cor_with_initial_genes[g, s] >= high_cor_threshold
    #                             ) + 1,
    #                             0,
    #                             1
    #                         )
    #                     }
    #                 )
    #             )
    #
    #             if(
    #                 nhc >= 1
    #             ) {
    #                 list(
    #                     gene = g,
    #                     number_high_correlations = nhc,
    #                     highest_correlation = max(cor_with_initial_genes[g, ])
    #                 )
    #             }
    #
    #         }
    #     )
    # )
    #
    # top_genes[
    #     ,
    #     score := (number_high_correlations - min(number_high_correlations))/
    #         max(number_high_correlations) +
    #         (highest_correlation - min(highest_correlation))/
    #         max(highest_correlation)
    #     ]
    #
    # top_genes <- rbind(
    #     data.table(gene = initial_genes, score = Inf),
    #     top_genes[order(-score), .(gene, score)]
    # )

    # Combine these genes with the genes from earlier:

    # genes <- unique(
    #     sort(
    #         c(
    #             genes,
    #             top_genes[1:top_genes_n, gene]
    #         )
    #     )
    # )

    # In case ccle and any extra data were provided, subset genes which occur in these
    # datasets:

    if(
        !is.null(ccle_cancer_type) & !is.null(ccle_data)
    ) {
        genes <- genes[genes %in% ccle_data$gene_id]
    } else if(
        !is.null(ccle_cancer_type) & is.null(ccle_data) |
        is.null(ccle_cancer_type) & !is.null(ccle_data)
    ) {
        warning(
            'Only one of ccle_cancer_type and ccle_data has been supplied.  ',
            'Please supply both.'
        )
    }

    if(
        !is.null(extra_data_source) & !is.null(extra_data)
    ) {
        genes <- genes[
            genes %in% extra_data[
                source == extra_data_source,
                gene
            ]
        ]
    } else if(
        !is.null(extra_data_source) & is.null(extra_data) |
        is.null(extra_data_source) & !is.null(extra_data)
    ) {
        warning(
            'Only one of extra_data_source and extra_data has been supplied.  ',
            'Please supply both.'
        )
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

        # The following is the actual filtering step.  Note that the scores are designed so that
        # a high score is good - we want to keep genes with high scores.

        if(genes_filter_method == 'max_cor') {

            if(max_initial_gene == max_cell_type) { # Either both are TRUE or both are FALSE

                if(!max_initial_gene & !max_cell_type) { # Give warning if they're both FALSE
                    warning(
                        'If genes_filter_method is "max_cor", at least one of max_initial_gene ',
                        'and max_cell_type should not be FALSE.  Using both for gene scoring.'
                    )
                }

                # genes_filtered <- genes_cor_with_initial_and_cell_types[
                #     max_initial_gene - max_cell_type >
                #         genes_filter_fun(max_initial_gene - max_cell_type),
                #     id
                # ]

                # I replaced the above with a process which uses a more flexible genes_filter_fun,
                # which should now act on the scores to return a logical or row numbers specifying
                # which ids to take.  E.g. if you want the scores which are greater than the mean score, set
                # genes_filter_fun = function(x) {x > mean(x)}.  If you just want to take the top
                # 250 genes, set genes_filter_fun = function(x) {1:250} - this works because we
                # put the scores in descending order, specifically for this case.

                gene_scores_table <- genes_cor_with_initial_and_cell_types[
                    ,
                    .(
                        id = id,
                        score = max_initial_gene - max_cell_type
                    )
                ][
                    order(-score)
                ]

                genes_filtered <- gene_scores_table[
                    genes_filter_fun(score),
                    id
                ]

            } else if(max_initial_gene) {

                gene_scores_table <- genes_cor_with_initial_and_cell_types[
                    ,
                    .(
                        id = id,
                        score = max_initial_gene
                    )
                ][
                    order(-score)
                ]

                genes_filtered <- gene_scores_table[
                    genes_filter_fun(score),
                    id
                ]

            } else { # Then max_cell_type must be TRUE (and max_initial_gene is FALSE)

                gene_scores_table <- genes_cor_with_initial_and_cell_types[
                    ,
                    .(
                        id = id,
                        score = -max_cell_type
                    )
                ][
                    order(-score)
                ]

                genes_filtered <- gene_scores_table[
                    genes_filter_fun(score),
                    id
                ]

            }

        } else { # Then genes_filter_method must be 'weighted_sum'

            # Check that the user hasn't specified genes_filter_method = 'weighted sum' but also
            # initial_gene_weights = FALSE and cell_type_weights = FALSE:

            if(isFALSE(initial_gene_weights) & isFALSE(cell_type_weights)) {
                initial_gene_weights <- NULL
                cell_type_weights <- NULL
                warning(
                    'If genes_filter_method is "weighted_sum", at least one of initial_gene_weights ',
                    'and cell_type_weights should not be FALSE.  Calculating all weights.'
                )
            }

            # The user can supply a named list of weights for the initial genes and the cell
            # types (the names must correspond to the initial_genes, resp. cell_types arguments).
            # If they don't, we construct them below.

            if(is.null(initial_gene_weights)) {

                # Define weights for initial genes as quantiles of their correlations with
                # other cell types (note we use '1 -' because we want to give more weight to
                # those initial genes that don't correlate highly with other cell types):

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
                                            .SD[
                                                ,
                                                cell_type_markers[
                                                    cell_type == ct &
                                                        gene %in% names(expression_data),
                                                    gene
                                                ],
                                                with = FALSE
                                            ]
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

                # Define weights for cell types as functions of their average correlations
                # with EMT/CAF markers: (EDIT: that's not really what this does...  It's
                # more like a function applied to the vector of average correlations of
                # the cell type markers with each of the EMT/CAF genes)

                cell_type_weights <- as.list(
                    apply(
                        sapply(
                            cell_types,
                            function(ct) {
                                expression_data[
                                    sample_ids,
                                    rowMeans(
                                        cor( # Should we swap these two arguments of cor()?
                                            .SD[, ..genes],
                                            .SD[
                                                ,
                                                cell_type_markers[
                                                    cell_type == ct & gene %in% names(expression_data),
                                                    gene
                                                ],
                                                with = FALSE
                                            ]
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
                    .(
                        score = sum(
                            as.numeric(
                                .SD[, ..initial_genes]
                            )*as.numeric(
                                initial_gene_weights[
                                    initial_genes # Makes the genes in the right order
                                ]
                            )
                        )
                    ),
                    by = id
                ][
                    order(-score)
                ]

                genes_filtered <- gene_scores_table[
                    genes_filter_fun(score),
                    id
                ]

            } else if(isFALSE(initial_gene_weights)) {

                gene_scores_table <- genes_cor_with_initial_and_cell_types[
                    ,
                    .(
                        score = - sum(
                            as.numeric(
                                .SD[, ..cell_types]
                            )*as.numeric(
                                cell_type_weights[
                                    cell_types # Makes the cell types in the right order
                                ]
                            )
                        )
                    ),
                    by = id
                ][
                    order(-score)
                ]

                genes_filtered <- gene_scores_table[
                    genes_filter_fun(score),
                    id
                ]

            } else {

                gene_scores_table <- genes_cor_with_initial_and_cell_types[
                    ,
                    .(
                        score = sum(
                            as.numeric(
                                .SD[, ..initial_genes]
                            )*as.numeric(
                                initial_gene_weights[
                                    initial_genes # Makes the genes in the right order
                                ]
                            )
                        ) - sum(
                            as.numeric(
                                .SD[, ..cell_types]
                            )*as.numeric(
                                cell_type_weights[
                                    cell_types # Makes the cell types in the right order
                                ]
                            )
                        )
                    ),
                    by = id
                ][
                    order(-score)
                ]

                genes_filtered <- gene_scores_table[
                    genes_filter_fun(score),
                    id
                ]

            }

        }

        # Just to make sure the initial genes are all in there:

        genes_filtered <- unique(
            c(
                genes_filtered,
                initial_genes
            )
        )

    } else {

        genes_filtered <- genes

    }





    # Regress genes in filtered list against sum of genes in filtered list, and take
    # residuals:

    resid_data <- as.data.table(
        sapply(
            genes_filtered,
            function(g) {

                # In the following, I have to put the backticks `` around g because
                # there are genes like 'NKX2-1', which gets interpreted as NKX2 - 1.

                lm(
                    formula(
                        paste0(
                            '`',
                            g,
                            '` ~ row_sums'
                        )
                    ),
                    data = cbind(
                        expression_data[, ..genes_filtered], # Subsetting genes_filtered
                        # makes it much faster!
                        row_sums = rowSums(
                            expression_data[, ..genes_filtered]
                        )
                    )
                )$residuals

            },
            simplify = FALSE,
            USE.NAMES = TRUE
        )
    )

    # Correlation matrix for resid_data:

    cor_mat <- cor(resid_data)

    # Ordering of correlation matrix for resid_data using SPIN_STS:

    if(!is.null(seed)) {
        set.seed(seed)
    }

    ordering <- ordering_fun(cor_mat)





    # Detect which genes are cancer EMT genes and which are CAF genes by looking at correlations
    # with purity and, if available, CCLE data and extra data, taking a consensus between them.
    # I will break ties with a very crude method - just use a few seemingly very reliable CAF
    # markers and see which end they're closest to.

    # Begin by calculating correlations with purity, with 3 possible methods.  Only one will
    # appear in the main figure, and this is also the one we will use to judge which are the
    # EMT and CAF ends.  The other two will be included with the diagnostic plots.

    cor_with_purity_scale <- cor(
        scaledt(
            expression_data[, c('id', ..genes_filtered)],
            margin = 1,
            scale = FALSE
        )[, -'id'],
        meta_data[sample_ids, purity]
    )[, 1]

    cor_with_purity_resid <- cor(
        resid_data,
        meta_data[sample_ids, purity]
    )[, 1]

    cor_with_purity_raw <- cor(
        expression_data[, ..genes_filtered],
        meta_data[sample_ids, purity]
    )[, 1]

    # CCLE comparison:

    if(
        !is.null(ccle_cancer_type) & !is.null(ccle_data)
    ) {

        # Should I change this to accept multiple CCLE cancer types?  E.g.:
        # ccle_cancer_type = c('stomach', 'oesophagus')

        ccle_data <- expression_data[
            ,
            .(
                id = names(.SD),
                tumours = colMeans(.SD),
                ccle = ccle_data[
                    names(expression_data[, -'id']),
                    ..ccle_cancer_type
                    ]
            ),
            .SDcols = -'id'
            ]

        mod_loess <- loess(
            tumours ~ ccle,
            ccle_data,
            span = 0.25,
            degree = 1,
            family = 'symmetric'
        )

        setkey(ccle_data, id)

        ccle_comp_diff <- predict(
            mod_loess,
            ccle_data[genes_filtered, ccle]
        ) + 1 -
            colMeans(
                expression_data[, ..genes_filtered]
            )

    }

    # Extra data:

    if(
        !is.null(extra_data_source) & !is.null(extra_data)
    ) {
        extra_data_score <- extra_data[
            source == extra_data_source
            ][
                genes_filtered,
                setNames(diff, gene) # Need a named vector for heat_map_bar() to work
                ]
    }

    # Store variables for checking the ordering in a list:

    ordering_check <- list(
        switch(
            (purity_cor_method == 'scale') + 1,
            cor_with_purity_resid,
            cor_with_purity_scale
        )
    )

    if(exists('ccle_comp_diff')) {
        ordering_check <- c(ordering_check, list(ccle_comp_diff))
    }

    if(exists('extra_data_score')) {
        ordering_check <- c(ordering_check, list(extra_data_score))
    }

    # This is the actual check:

    ordering_check <- sum(
        sapply(
            ordering_check,
            function(x) {
                as.numeric(
                    lm(x[ordering] ~ seq(1, length(genes_filtered), 1))$coeff[2]
                )
            }
        ) > 0 # Greater than zero suggests the ordering is wrong
    )/length(ordering_check)

    # If a consensus suggests the ordering is wrong, reverse the ordering; if there is no
    # consensus, try to break the tie using the average position of the CAF markers (if this
    # still doesn't break the tie, or if these CAF markers don't appear in the filtered gene
    # list, then nothing happens):

    if(ordering_check > 0.5) {
        ordering <- rev(ordering)
    } else if(ordering_check == 0.5) {
        if(
            sum(
                c(
                    'COL1A1',
                    'COL1A2',
                    'COL3A1',
                    'COL6A3',
                    'THY1'
                ) %in% genes_filtered
            ) > 0
        ) {
            if(
                mean(
                    which(
                        genes_filtered[ordering] %in% c(
                            'COL1A1',
                            'COL1A2',
                            'COL3A1',
                            'COL6A3',
                            'THY1'
                        )
                    )
                ) < length(genes_filtered)/2
            ) {
                ordering <- rev(ordering)
            }
        }
    }





    # Heatmap, annotations and purity bar:

    htmp <- heat_map(
        cor_mat,
        ordering,
        axis_title_y = 'genes',
        axis_title_x = '',
        axis_text_x = '',
        axis_text_y = '',
        colour_limits = heatmap_colour_limits,
        legend_breaks = heatmap_legend_breaks,
        plot_margin = c(5.5, 5.5, 0, 5.5),
        plot_title = plot_title
    )

    interesting_labels <- genes_filtered[ordering]

    interesting_labels[!(interesting_labels %in% heatmap_annotations)] <- ''

    axis_labels <- heat_map_labels_repel(interesting_labels)

    purity_bar <- heat_map_bar(
        switch(
            (purity_cor_method == 'scale') + 1,
            cor_with_purity_resid,
            cor_with_purity_scale
        ),
        ordering,
        fun = purity_fun,
        colours = heatmap_bar_colours,
        colour_limits = purity_colour_limits,
        axis_title_y = 'Correlation\nwith purity',
        legend_breaks = purity_legend_breaks,
        legend.direction = 'horizontal',
        legend.title = element_blank(),
        legend.key.height = bar_legend_height,
        legend.key.width = bar_legend_width,
        axis.title.y = element_text(angle = 0, vjust = 0.5),
        axis.title.x = element_blank(),
        plot_margin = c(0, 5.5, 2, 5.5)
    )

    plots <- list(
        heatmap = htmp,
        axis_labels = axis_labels,
        purity_bar = purity_bar
    )





    # CCLE bar:

    if(
        !is.null(ccle_cancer_type) & !is.null(ccle_data)
    ) {

        ccle_bar <- heat_map_bar(
            ccle_comp_diff,
            ordering,
            fun = ccle_fun,
            colours = heatmap_bar_colours,
            colour_limits = ccle_colour_limits,
            axis_title_y = 'Tumours vs.\ncell lines',
            legend_breaks = ccle_legend_breaks,
            legend.direction = 'horizontal',
            legend.title = element_blank(),
            legend.key.height = bar_legend_height,
            legend.key.width = bar_legend_width,
            axis.title.y = element_text(angle = 0, vjust = 0.5),
            axis.title.x = element_blank(),
            plot_margin = c(0, 5.5, 2, 5.5)
        )

        plots <- c(plots, ccle_bar = list(ccle_bar))

    }





    # Extra data bar:

    if(
        !is.null(extra_data_source) & !is.null(extra_data)
    ) {

        setkey(extra_data, gene)

        extra_bar <- heat_map_bar(
            extra_data_score,
            ordering,
            fun = extra_fun,
            colours = heatmap_bar_colours,
            colour_limits = extra_colour_limits,
            axis_title_y = extra_title,
            legend_breaks = extra_legend_breaks,
            legend.direction = 'horizontal',
            legend.title = element_blank(),
            legend.key.height = bar_legend_height,
            legend.key.width = bar_legend_width,
            axis.title.y = element_text(angle = 0, vjust = 0.5),
            axis.title.x = element_blank(),
            plot_margin = c(0, 5.5, 2, 5.5)
        )

        plots = c(plots, extra_bar = list(extra_bar))

    }





    # Make diagnostic plots:

    if(purity_cor_method == 'scale') {
        purity_cor_vecs <- list(
            resid = cor_with_purity_resid,
            raw = cor_with_purity_raw
        )
    } else {
        purity_cor_vecs <- list(
            scale = cor_with_purity_scale,
            raw = cor_with_purity_raw
        )
    }

    diagnostics <- list(

        alternative_purity_cor = sapply(
            names(purity_cor_vecs),
            function(li_name) {
                heat_map_bar(
                    purity_cor_vecs[[li_name]],
                    ordering,
                    fun = purity_fun,
                    colours = grDevices::rainbow(50),
                    colour_limits = purity_colour_limits,
                    axis_title_y = paste0('Cor. with\npurity (', li_name, ')'),
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

    if(!is.null(genes_filter_fun)) {

        setkey(genes_cor_with_initial_and_cell_types, id)

        genes_cor_with_initial_and_cell_types[
            genes_filtered,
            epithelial := expression_data[
                ,
                rowMeans(
                    cor(
                        .SD[, ..genes_filtered],
                        .SD[
                            ,
                            c(
                                'CDH1',
                                'EPCAM',
                                'SFN',
                                names(.SD)[
                                    grep('^KRT[0-9]|^KRTD', names(.SD))
                                    ]
                            ),
                            with = FALSE
                            ]
                    )
                )
                ]
            ]

        cell_type_bars_colour_limits <- genes_cor_with_initial_and_cell_types[
            genes_filtered,
            quantile(
                apply(
                    .SD,
                    2,
                    function(x) purity_fun(x[ordering])
                ),
                c(0, 1)
            ),
            .SDcols = c(cell_types, 'epithelial')
            ]

        diagnostics <- c(

            diagnostics,

            list(

                cell_type_bars = sapply(
                    c(cell_types, 'epithelial'),
                    function(ct) {
                        heat_map_bar(
                            genes_cor_with_initial_and_cell_types[
                                genes_filtered,
                                setNames(get(ct), id)
                                ],
                            ordering,
                            fun = purity_fun,
                            colours = grDevices::rainbow(50), # I think I prefer this to the Spectral palette
                            # colours = rev(colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))(50)),
                            colour_limits = as.numeric(cell_type_bars_colour_limits),
                            legend_breaks = c(
                                ceiling(
                                    as.numeric(cell_type_bars_colour_limits)[1]*10
                                )/10,
                                floor(
                                    as.numeric(cell_type_bars_colour_limits)[2]*10
                                )/10
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
                        genes_cor_with_initial_and_cell_types[
                            initial_genes,
                            c('id', cell_types, 'epithelial'),
                            with = FALSE
                        ],
                        id.vars = 'id',
                        variable.name = 'cell_type',
                        value.name = 'correlation'
                    )
                ) +
                    geom_raster(
                        aes(
                            x = cell_type,
                            y = id,
                            fill = correlation
                        )
                    ) +
                    scale_y_discrete(
                        expand = c(0, 0)
                    ) +
                    scale_x_discrete(
                        expand = c(0, 0)
                    ),

                genes_cell_types_correlations = genes_cor_with_initial_and_cell_types[
                    genes_filtered
                ],

                cell_type_lms = genes_cor_with_initial_and_cell_types[
                    genes_filtered
                ][
                    ordering,
                    sapply(
                        .SD,
                        function(ct) {
                            setNames(
                                lm(ct ~ .I)$coeff,
                                c('intercept', 'slope')
                            )
                        },
                        simplify = TRUE,
                        USE.NAMES = TRUE
                    ),
                    .SDcols = c(cell_types, 'epithelial')
                ],

                gene_cell_type_scores = gene_scores_table,

                gene_cell_type_scores_density = ggplot(gene_scores_table) +
                    geom_density(aes(score)) +
                    theme_test()

            )

        )

    }

    if(!is.null(ccle_cancer_type) & !is.null(ccle_data)) {

        diagnostics <- c(
            diagnostics,
            tumours_vs_cell_lines = list(
                ggplot(
                    cbind(
                        ccle_data,
                        in_genes_filtered = ccle_data$id %in% genes_filtered
                    )
                ) +
                    geom_point(
                        aes(x = ccle, y = tumours, col = in_genes_filtered)
                    ) +
                    scale_colour_manual(values = c('black', 'deepskyblue')) +
                    # The following just makes sure the blue points appear on top, for clarity:
                    geom_point(
                        aes(x = ccle, y = tumours),
                        data = ccle_data[id %in% genes_filtered],
                        colour = 'deepskyblue'
                    ) +
                    geom_line(
                        aes(x = ccle, y = predict(mod_loess, newdata = ccle_data)),
                        colour = 'red',
                        size = 1.5
                    ) +
                    labs(
                        colour = 'In filtered\ngenes list'
                    )
            )
        )

    }

    # The following is a simple diagnostic measure that says how often we have a "strong"
    # positive association between EMT and TME, where a "strong" association is one where the
    # regression slope is less than -1/2500 (corresponding to a change in Pearson correlation
    # of 0.1 over 250 genes):

    # sum(diagnostics$cell_type_lms['slope', cell_types] < -1/2500)

    # The following would give the cell types that have such associations:

    # names(
    #     diagnostics$cell_type_lms[
    #         'slope',
    #         diagnostics$cell_type_lms['slope', cell_types] < -1/2500
    #     ]
    # )





    # Output:

    list(
        plots = plots,
        sample_ids = sample_ids,
        genes_filtered = genes_filtered,
        ordering = ordering,
        diagnostics = diagnostics
    )

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
    genes_from_tcga_fun = function(x) { # The name "genes_from_tcga_fun" is misleading - they come from <expression_data>
        top_cols_by_fun_cor(x, initial = initial_genes)[ # Using defaults for <FUN> and <threshold>
            1:min(.N, 300), # Take no more than 300 genes
            id
        ]
    },
    genes_filter_fun = function(x) 1:250,
    gene_weights_fun = median,
    cell_type_weights = NULL,
    initial_gene_weights = NULL,
    caf_markers = c('COL1A1', 'COL1A2', 'COL3A1', 'COL6A3', 'THY1'),
    seed = NULL,
    ordering_fun = function(x) {as.vector(seriation::seriate(as.dist(1 - x), method = 'SPIN_STS')[[1]])} # This is meant for a correlation matrix

) {

    # <caf_markers> is used for tie-breaking in deciding whether to reverse the ordering for the correlation matrix.
    if(!is.null(genes_filter_fun)) {genes_filter_fun <- match.fun(genes_filter_fun)}
    # If we want to filter the genes using an absolute constant, we can use a constant function, e.g. genes_filter_fun <- function(x) 1
    gene_weights_fun <- match.fun(gene_weights_fun)
    ordering_fun <- match.fun(ordering_fun)
    # We could do hierarchical clustering using hclust(1 - x)$order, or something.

    # Set keys, in case this hasn't already been done:
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
        # In other cases where subtypes_data is not supplied, we'll be using sample_ids anyway
        # (other combinations eliminated above).
        stop('Please supply subtypes_data.')
    }

    # Table showing with which of the above 3 if statements each combination was eliminated
    # ('x' means the argument is NULL; '/' means something was supplied to it):

    # sample_ids        x / x x x / / / x x x / / / x /
    # tcga_cancer_types x x / x x / x x / / x / / x / /
    # subtypes          x x x / x x / x / x / / x / / /
    # subtypes_data     x x x x / x x / x / / x / / / /
    # which if stmnt:   1     1 1 2 2   3   1 2 2 2   2

    # If we're not using sample_ids by default, get the appropriate sample IDs using
    # tcga_cancer_types, and subtypes if it is supplied:
    if(is.null(sample_ids)) { # Then we must have !is.null(tcga_cancer_types)
        if(is.null(subtypes)) { # Then get all IDs corresponding to these cancer types
            sample_ids <- meta_data[cancer_type %in% tcga_cancer_types, id]
        } else { # Then we must have !is.null(subtypes_data)
            # Set key, in case this hasn't already been done:
            setkey(subtypes_data, id)
            if(is.null(ref_for_subtypes)) {
                sample_ids <- subtypes_data[meta_data[cancer_type %in% tcga_cancer_types, id]][get(subtypes_var) %in% subtypes, id]
            } else {
                # Optionally supply reference for subtypes, in case the same subtype name
                # appears in multiple references.
                sample_ids <- subtypes_data[meta_data[cancer_type %in% tcga_cancer_types, id]][subtype_ref == ref_for_subtypes & get(subtypes_var) %in% subtypes, id]
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

    # In case ccle and any extra data were provided, subset genes which occur in these
    # datasets:

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

        # The following is the actual filtering step.  Note that the scores are designed so that
        # a high score is good - we want to keep genes with high scores.

        if(isFALSE(initial_gene_weights) & isFALSE(cell_type_weights)) {
            initial_gene_weights <- NULL
            cell_type_weights <- NULL
            warning('At least one of initial_gene_weights and cell_type_weights should not be FALSE.  Calculating all weights.')
        }

        # The user can supply a named list of weights for the initial genes and the cell
        # types (the names must correspond to the initial_genes, resp. cell_types arguments).
        # If they don't, we construct them below.

        if(is.null(initial_gene_weights)) {

            # Define weights for initial genes as quantiles of their correlations with
            # other cell types (note we use '1 -' because we want to give more weight to
            # those initial genes that don't correlate highly with other cell types):

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

            # Define weights for cell types as functions of their average correlations
            # with EMT/CAF markers: (EDIT: that's not really what this does...  It's
            # more like a function applied to the vector of average correlations of
            # the cell type markers with each of the EMT/CAF genes)

            cell_type_weights <- as.list(
                apply(
                    sapply(
                        cell_types,
                        function(ct) {
                            expression_data[
                                sample_ids,
                                rowMeans(
                                    cor( # Should we swap these two arguments of cor()?
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

        # Just to make sure the initial genes are all in there:
        genes_filtered <- unique(c(genes_filtered, initial_genes))

    } else {

        genes_filtered <- genes

    }





    # Regress genes in filtered list against sum of genes in filtered list, and take
    # residuals:
    resid_data <- as.data.table(
        sapply(
            genes_filtered,
            function(g) {
                # In the following, I have to put the backticks `` around g because there are genes like 'NKX2-1', which gets interpreted as NKX2 - 1.
                lm(
                    formula(paste0('`', g, '` ~ row_sums')),
                    data = cbind(
                        expression_data[, ..genes_filtered], # Subsetting genes_filtered
                        # makes it much faster!
                        row_sums = rowSums(expression_data[, ..genes_filtered])
                    )
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





    # Detect which genes are cancer EMT genes and which are CAF genes by looking at correlations
    # with purity and, if available, CCLE data and extra data, taking a consensus between them.
    # I will break ties with a very crude method - just use a few seemingly very reliable CAF
    # markers and see which end they're closest to.

    # Begin by calculating correlations with purity, with 3 possible methods.  Only one will
    # appear in the main figure, and this is also the one we will use to judge which are the
    # EMT and CAF ends.  The other two will be included with the diagnostic plots.
    cor_with_purity_scale <- cor(scaledt(expression_data[, c('id', ..genes_filtered)], margin = 1, scale = FALSE)[, -'id'], meta_data[sample_ids, purity])[, 1]
    cor_with_purity_resid <- cor(resid_data, meta_data[sample_ids, purity])[, 1]
    cor_with_purity_raw <- cor(expression_data[, ..genes_filtered], meta_data[sample_ids, purity])[, 1]

    # CCLE comparison:
    if(!is.null(ccle_cancer_type) & !is.null(ccle_data)) {

        # Should I change this to accept multiple CCLE cancer types?  E.g.: ccle_cancer_type = c('stomach', 'oesophagus')

        ccle_data <- expression_data[
            ,
            .(
                id = names(.SD),
                tumours = colMeans(.SD),
                ccle = ccle_data[
                    names(expression_data[, -'id']),
                    ..ccle_cancer_type
                ][[ccle_cancer_type]] # Didn't need these double brackets before - I assume it's to do with R version 3.6.3 and data.table version 1.12.8
            ),
            .SDcols = -'id'
        ]

        mod_loess <- loess(tumours ~ ccle, ccle_data, span = 0.25, degree = 1, family = 'symmetric')

        setkey(ccle_data, id)

        ccle_comp_diff <- predict(mod_loess, ccle_data[genes_filtered, ccle]) + 1 - colMeans(expression_data[, ..genes_filtered])

    }

    # Extra data:
    if(!is.null(extra_data_source) & !is.null(extra_data)) {
        extra_data_score <- extra_data[source == extra_data_source][genes_filtered, setNames(diff, gene)] # Need a named vector for heat_map_bar() to work
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

    # If a consensus suggests the ordering is wrong, reverse the ordering; if there is no
    # consensus, try to break the tie using the average position of the CAF markers (if this
    # still doesn't break the tie, or if these CAF markers don't appear in the filtered gene
    # list, then nothing happens):
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





# Note the following is useful for displaying colours.  In this example I'm using a subset
# of colours from an RColorBrewer palette and interpolating them with colorRampPalette:

# ggplot(
#     data = data.frame(x = 1:50, y = 1, clr = factor(1:50, levels = 1:50))
# ) +
#     geom_raster(aes(x = x, y = y, fill = clr)) +
#     scale_fill_manual(
#         values = colorRampPalette(RColorBrewer::brewer.pal(11, 'PiYG')[c(1:3, 6, 9:11)])(50)
#     ) +
#     scale_x_continuous(expand = c(0, 0)) +
#     scale_y_continuous(expand = c(0, 0)) +
#     theme(
#         axis.text = element_blank(),
#         axis.title = element_blank(),
#         axis.ticks = element_blank()
#     )

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
    epithelial_markers = c('CDH1', 'EPCAM', 'SFN', switch(is.null(expression_data) + 1, names(expression_data)[grep('^KRT[0-9]|^KRTD', names(expression_data))], NULL))

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
        axis_title_x = '', # Could set this to heatmap_axis_title if heatmap_annotations is NULL or empty
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
		edge = plyr::mapvalues(heatmap_annotations_side, c('bottom', 'left', 'top', 'right'), c('top', 'right', 'bottom', 'left'), warn_missing = FALSE),
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

                # Remove columns from expression_data which contain NAs or have zero standard
                # deviation:

                expression_data <- cbind(
                    expression_data[data$sample_ids, .(id)],
                    expression_data[data$sample_ids, sapply(.SD[, -'id'], function(x) {switch((any(is.na(x)) || sd(x) == 0) + 1, x, NULL)})]
                )

                # Subset those epithelial markers that still appear in expression_data after
                # removing columns:

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

                    gene_cell_type_scores_density = ggplot(data$scores_table) +
                        geom_density(aes(score)) +
                        theme_test()

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
                        # The following just makes sure the blue points appear on top, for clarity:
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
    genes_to_rank <- sort(unique(unlist(lapply(deconv_data, `[[`, 'genes_filtered')))) # function(li) li[[pmatch('gene', names(li))]]
    sapply(
        deconv_data,
        function(li) {
            # li_genes <- li[[pmatch('gene', names(li))]][li[[pmatch('order', names(li))]]]
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

    # By default, scores will be calculated for all genes in <deconv_data>.  The user can
    # supply a character vector to the <additional_genes> argument, in which case scores
    # will be calculated for genes in <additional_genes> as well as for those in <deconv_data>.
    # If the user supplies a character vector to the <genes> argument, then scores will only
    # be calculated for the genes in <genes> and not for those in <deconv_data> (or for any
    # supplied to <additional_genes>).  In every case, the scores will be calculated using
    # average correlation with genes in <deconv_data>.  If <head_or_tail> is set to 'head' or
    # 'tail', we use respectively the genes from the head and tail of the genes in each item
    # in <deconv_data>.  If <head_or_tail> is 'both' (the default), we calculate the
    # difference between the average correlations with the head and tail genes.  The argument
    # <n> is the n to use in the head() and tail() functions.

    # If <transform_data> is TRUE, correlations will be calculated in the transformed space
    # defined by the genes and samples in <deconv_data>, that is, for each cancer type, the
    # gene vectors are replaced with the residuals of a linear regression against sample sums.

    # <scale_fun> is a function to be applied to the genes by cancer types matrix of scores.
    # For example, scale_fun = scale would scale the scores relative to the global mean and
    # standard deviation; scale_fun = function(x) apply(x, 2, scale) would scale per cancer
    # type; scale_fun = function(x) apply(x, 2, function(y) y/max(abs(y))) would scale the
    # scores for each cancer type to the interval [-1, 1].  You can set <scale_fun_margin> to
    # 1 or 2 if you want to apply <scale_fun> along the margins of the matrix.  1 means applying
    # to genes, 2 means applying to cancer types.

    head_tail <- match.arg(head_tail)
    if(head_tail == 'both') head_tail <- c('head', 'tail')

    if(!is.null(scale_fun)) {scale_fun <- match.fun(scale_fun)}

    setkey(expression_data, id)

    if(!is.null(genes)) {
        genes_to_score <- genes
    } else if(!is.null(additional_genes)) {
        genes_to_score <- unique(c(unlist(lapply(deconv_data, `[[`, 'genes_filtered')), additional_genes)) # function(li) li[[pmatch('gene', names(li))]]
    } else {
        genes_to_score <- unique(unlist(lapply(deconv_data, `[[`, 'genes_filtered'))) # function(li) li[[pmatch('gene', names(li))]]
    }

    scores <- sapply(

        deconv_data,

        function(li) {

            # li_sample_ids <- li[[pmatch('sample', names(li))]]
            # li_ordered_genes <- li[[pmatch('gene', names(li))]][li[[pmatch('order', names(li))]]]

            li_sample_ids <- li$sample_ids
            li_ordered_genes <- li$genes_filtered[li$ordering]

            for(fun_name in head_tail) {assign(paste0(fun_name, '_genes'), match.fun(fun_name)(li_ordered_genes, n))}

            genes_for_scores_data <- unique(c(genes_to_score, unlist(mget(paste0(head_tail, '_genes')))))

            if(transform_data) {
                scores_data <- cbind(
                    expression_data[li_sample_ids, ..genes_for_scores_data],
                    row_sums = rowSums(
                        expression_data[
                            li_sample_ids,
                            # Use only the filtered genes from the deconv for row sums:
                            ..li_ordered_genes # They don't have to be ordered, but no harm in it
                        ]
                    )
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
                    , # Divide by 2 in the following to scale to interval [-1, 1]
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
                scores <- t(apply(scores, 1, scale_fun)) # Seems that whichever margin you apply along becomes the columns of the result
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
    # scores <- scores[order(rownames(scores)),]

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

    # <order_genes_fun> and <order_analyses_fun> can each be 'hclust', 'seriate' or a user-
    # defined function.  In the former two cases, hclust() and seriate() functions will be
    # used respectively, and the corresponding methods may be selected using the
    # <order_genes_method> and <order_analyses_method> arguments.  An hclust object will
    # not be generated if the user defines their own function (should add functionality
    # for that).

    # The '...' argument is for extra arguments to the theme() function in the
    # construction of the heatmap.

    if(!(typeof(order_genes_fun) == 'character' && order_genes_fun %in% c('hclust', 'seriate'))) {
        order_genes_fun <- match.fun(order_genes_fun)
    }

    if(!(typeof(order_analyses_fun) == 'character' && order_analyses_fun %in% c('hclust', 'seriate'))) {
        order_analyses_fun <- match.fun(order_analyses_fun)
    }

    # order_genes_fun <- match.arg(order_genes_fun)
    # order_analyses_fun <- match.arg(order_analyses_fun)

    setkey(scores_data, gene)

    # Since the genes in scores_data are in alphabetical order, ordering_genes is the
    # ordering with respect to the alphabetical order.

    if(typeof(order_genes_fun) == 'character' && order_genes_fun == 'seriate') {

        if(order_genes_method %in% seriation::list_seriation_methods('matrix')) {
            ordering_obj_genes <- seriation::seriate(as.matrix(scores_data[, -'gene']), method = order_genes_method, margin = 1)
        } else if(order_genes_method %in% seriation::list_seriation_methods('dist')) {
            # Recall the dist() function calculates distances between rows.
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

    # ordering_analyses is with respect to the ordering in scores_data, which should be
    # the same as the ordering in the deconv_data used as input to the deconv_scores()
    # function.

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

    # If we wanted to keep scores_data as a matrix, note that melt does work for matrics as
    # well, though you have less choice over column names.

    scores_melted <- melt(
        scores_data,
        id.vars = 'gene',
        measure.vars = names(scores_data[, -'gene'])[ordering_analyses],
        variable.name = 'analysis',
        value.name = 'score',
        variable.factor = FALSE # This is REALLY important!  Otherwise the 'analysis'
        # column will be a factor, which messes up stuff down the line.
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
            # plot.margin = switch(
            #     (exists('hclust_genes') & draw_dendro) + 1,
            #     switch(
            #         (exists('hclust_analyses') & draw_dendro) + 1,
            #         unit(c(5.5, 5.5, 5.5, 5.5), 'pt'),
            #         unit(c(1, 5.5, 5.5, 5.5), 'pt')
            #     ),
            #     switch(
            #         (exists('hclust_analyses')) + 1, # Already know draw_dendro is TRUE
            #         unit(c(5.5, 1, 5.5, 5.5), 'pt'),
            #         unit(c(1, 1, 5.5, 5.5), 'pt')
            #     )
            # ),
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





# The following function plots the aligned heatmaps and dendrogram(s) from the output
# of the deconv_heatmap() function.  The egg and cowplot packages both have trouble
# aligning two dendrograms to the heatmap, because for some reason they add additional
# padding between the heatmap and the dendrograms.  This seems to be because the
# aligning functions always choose the maximum of 2 values, one being the padding I
# specified in constructing the heatmap (1pt) and the other being 5.5pt.  When you
# look at the grob widths or heights, all the entries contain an expression using
# max(), and the problematic entry in this case is of the form max(5.5pt, 1pt). The
# only way I found to fix this is to manually change these grob entries to 1pt.

deconv_heatmap_dendro_plot <- function(

    deconv_heatmap_list,
    show_dendros = na.omit(
        sapply(strsplit(names(deconv_heatmap_list), 'dendro_'), `[`, 2)
    ),
    rel_heights = c(1, 8),
    rel_widths = switch(
        (length(show_dendros) == 2 || show_dendros == 'analyses') + 1,
        c(6, 1, 1),
        c(6, 1)
    ),
    padding = unit(0, 'pt'),
    ...

) {

    # The '...' argument is for additional arguments for guide_colourbar(), which can
    # be used to change the formatting of the legend.

    if(sum(startsWith(names(deconv_heatmap_list), 'dendro')) == 0) {
        stop('There are no dendrograms in this plot list.')
    }

    show_dendros <- match.arg(
        show_dendros,
        c('genes', 'analyses'),
        several.ok = TRUE
    )

    grob_htmp <- ggplotGrob(
        deconv_heatmap_list$heatmap + theme(legend.position = 'none')
    )

    leg <- get_legend(
        deconv_heatmap_list$heatmap +
            guides(fill = guide_colourbar(...)) # Use this to change legend format
    )

    if('genes' %in% show_dendros) {

        aligned_genes <- align_plots(
            deconv_heatmap_list$heatmap + theme(legend.position = 'none'),
            deconv_heatmap_list$dendro_genes,
            align = 'h'
        )

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

    # Use this to show which elements of the grobs we need to adjust:

    # grob_dend_analyses$widths
    # grob_htmp$widths

    # grob_dend_genes$heights
    # grob_htmp$heights

    # The latter two options are a bit of a bodge because the legend is not aligned with
    # the heatmap - it is centred vertically relative to the entire plot.  I don't know
    # how to align legends with cowplot.  I think we might have to make a dummy legend
    # for the analyses dendrogram or a dummy plot for the heatmap legend.

    if(length(show_dendros) == 2) {

        plot_grid(
            grob_dend_analyses,
            leg,
            grob_htmp,
            grob_dend_genes,
            nrow = 2,
            ncol = 2,
            rel_heights = rel_heights,
            rel_widths = rel_widths
        )

    } else if(show_dendros == 'genes') {

        plot_grid(
            grob_htmp,
            grob_dend_genes,
            leg,
            nrow = 1,
            ncol = 3,
            rel_widths = rel_widths
        )

    } else { # Then show_dendros == 'analyses'

        plot_grid(
            plot_grid(
                grob_dend_analyses,
                grob_htmp,
                ncol = 1,
                nrow = 2,
                rel_heights = rel_heights
            ),
            leg,
            nrow = 1,
            ncol = 2,
            rel_widths = rel_widths
        )

    }

}





# I think the following two functions should return the entire hclust object (if the hclust
# option is selected), so that we can use it afterwards to create/examine the dendrogram.

commonality_heatmap <- function(

    expression_data,
    ids_genes_ordering_list,
    extra_ordered_genes_list = NULL,
    two_sided = TRUE,
    transform_data = TRUE,
    score_type = c('correlation', 'expression_level'),
    quantile_prob = 0.25,
    n_head = 25,
    n_tail = n_head,
    genes_to_include = c(
        'SNAI1',
        'SNAI2',
        'ZEB1',
        'ZEB2',
        'TWIST1',
        'VIM'
    ),
    normalise_scores = TRUE,
    seriate_genes_fun = c('hclust', 'seriate'),
    seriate_genes_method = 'average',
    seriate_analyses_fun = c('hclust', 'seriate'),
    seriate_analyses_method = 'average',
    x_axis_text_angle = 55,
    legend_title = 'pEMT-CAF\nscore'

) {

    # ids_genes_ordering_list should be a list where each element has a vector of TCGA sample
    # IDs, a vector of gene names and, optionally, a numerical vector giving an ordering for
    # the genes.  This is designed so that a list of outputs of the deconvolve_emt_caf()
    # function can be put straight into this function.

    score_type <- match.arg(score_type)

    seriate_genes_fun <- match.arg(seriate_genes_fun)
    seriate_analyses_fun <- match.arg(seriate_analyses_fun)

    for(li in ids_genes_ordering_list) {
        if(is.na(pmatch('order', names(li)))) {
            li$ordering <- 1:length(li[[pmatch('genes', names(li))]])
        }
    }

    # If supplied, extra_ordered_genes_list should have the same names as
    # ids_genes_ordering_list.  If they don't, take the ones that are common between them.

    if(
        !is.null(extra_ordered_genes_list)
    ) {
        if(
            sum(
                names(ids_genes_ordering_list) != names(extra_ordered_genes_list)
            ) > 0
        ) {
            common_names <- intersect(
                names(ids_genes_ordering_list),
                names(extra_ordered_genes_list)
            )
            if(length(common_names) > 0) {
                ids_genes_ordering_list <- ids_genes_ordering_list[common_names]
                extra_ordered_genes_list <- extra_ordered_genes_list[common_names]
                warning(
                    'ids_genes_ordering_list and extra_ordered_genes_list do not have the same ',
                    'names.  Taking only the names they have in common.'
                )
            } else {
                stop('ids_genes_ordering_list and extra_ordered_genes_list have no names in common.')
            }
        }
    }

    # Find genes common to the gene vectors in ids_genes_ordering_list:

    # Rank genes using a quantile of their positions in the ordered gene list (relative to the
    # length of the gene list):

    if(!is.null(extra_ordered_genes_list)) {
        genes_list_to_rank <- extra_ordered_genes_list
    } else {
        genes_list_to_rank <- sapply(
            ids_genes_ordering_list,
            function(li) {
                li[[pmatch('gene', names(li))]][
                    li[[pmatch('order', names(li))]]
                    ]
            },
            simplify = FALSE,
            USE.NAMES = TRUE
        )
    }

    gene_ranking <- data.table(
        gene = unique(
            unlist(
                genes_list_to_rank
            )
        )
    )

    if(two_sided) {

        gene_ranking[
            ,
            c('quantile_rank_emt', 'quantile_rank_caf') := as.list(
                apply(
                    sapply(
                        genes_list_to_rank,
                        function(li) {
                            if(gene %in% li) {
                                c(
                                    which(li == gene),
                                    which(rev(li) == gene)
                                )/length(
                                    li
                                )
                            } else {
                                c(1, 1)
                            }
                        }
                    ),
                    1,
                    quantile,
                    quantile_prob
                )
            ),
            by = gene
            ]

        # Take top n EMT and CAF genes based on this ranking and combine with the genes we want
        # to include unconditionally:

        genes_for_heatmap <- sort(
            unique(
                c(
                    gene_ranking[order(quantile_rank_emt)][1:n_head, gene],
                    gene_ranking[order(quantile_rank_caf)][1:n_tail, gene],
                    genes_to_include
                )
            )
        )

    } else {

        gene_ranking[
            ,
            quantile_rank := quantile(
                sapply(
                    genes_list_to_rank,
                    function(li) {
                        if(gene %in% li) {
                            which(li == gene)/length(li)
                        } else {
                            1
                        }
                    }
                ),
                quantile_prob
            ),
            by = gene
            ]

        # Take top n genes based on this ranking and combine with the genes we want to
        # include unconditionally:

        genes_for_heatmap <- sort(
            unique(
                c(
                    gene_ranking[order(quantile_rank)][1:n_head, gene],
                    genes_to_include
                )
            )
        )

    }

    # Take total gene set for scores data:

    genes_for_scores_data <- unique(
        c(
            genes_to_include,
            gene_ranking$gene
        )
    )

    # For each cancer type/subtype, calculate pEMT-CAF scores using correlations with the genes
    # at the extremes:

    setkey(expression_data, id)

    scores <- sapply(

        ids_genes_ordering_list,

        function(li) {

            li_sample_ids <- li[[pmatch('sample', names(li))]]

            li_ordered_genes <- li[[pmatch('gene', names(li))]][
                li[[pmatch('order', names(li))]]
                ]

            emt_genes <- head(li_ordered_genes, 20)
            caf_genes <- tail(li_ordered_genes, 20)

            if(transform_data) {

                scores_data <- cbind(
                    expression_data[
                        li_sample_ids,
                        ..genes_for_scores_data
                        ],
                    row_sums = rowSums(
                        expression_data[
                            li_sample_ids,
                            # Use only the filtered genes from the deconv for row sums:
                            ..li_ordered_genes # They don't have to be ordered, but no harm in it
                            ]
                    )
                )

                scores_data <- as.data.table(
                    sapply(
                        genes_for_scores_data,
                        function(g) {

                            # In the following, I have to put the backticks `` around g because
                            # there are genes like 'NKX2-1', which gets interpreted as NKX2 - 1.

                            lm(
                                formula(
                                    paste0(
                                        '`',
                                        g,
                                        '` ~ row_sums'
                                    )
                                ),
                                data = scores_data
                            )$residuals

                        },
                        simplify = FALSE,
                        USE.NAMES = TRUE
                    )
                )

            } else {

                scores_data <- expression_data[
                    li_sample_ids,
                    ..genes_for_scores_data
                    ]

            }

            if(score_type == 'correlation') {

                # Take averages of head and tail genes in resid_data, so that we can calculate
                # correlations of the "scaled" genes_for_heatmap genes with these averages:

                if(two_sided) {

                    scores_li <- scores_data[
                        ,
                        (
                            rowMeans(
                                cor(
                                    .SD,
                                    expression_data[ # Use expression_data instead of .SD here
                                        # in case emt_genes are not all in scores_data, which
                                        # may happen when we're using extra_ordered_genes_list
                                        li_sample_ids,
                                        ..emt_genes
                                        ]
                                )
                            ) - rowMeans(
                                cor(
                                    .SD,
                                    expression_data[
                                        li_sample_ids,
                                        ..caf_genes
                                        ]
                                )
                            )
                        )/2 # Divide by 2 to scale to interval [-1, 1]
                        ]

                    # If you do want to do correlation with average instead of average
                    # correlation, use the following:

                    # scores_li <- scores_data[
                    #     ,
                    #     (
                    #         cor(
                    #             .SD,
                    #             rowMeans(.SD[, ..emt_genes])
                    #         )[, 1] - cor(
                    #             .SD,
                    #             rowMeans(.SD[, ..caf_genes])
                    #         )[, 1]
                    #     )/2
                    # ]

                } else {

                    scores_li <- scores_data[
                        ,
                        rowMeans(
                            cor(
                                .SD,
                                expression_data[
                                    li_sample_ids,
                                    ..emt_genes
                                    ]
                            )
                        )
                        ]

                }

            } else {

                scores_li <- colMeans(scores_data)

            }

            scores_li

        },

        simplify = TRUE,
        USE.NAMES = TRUE

    )

    if(normalise_scores) {
        scores <- scores/max(abs(scores[genes_for_heatmap, ]))
    }

    #Make into a data.frame for use in ggplot():

    # scores <- cbind(
    #     gene = names(scores[[1]]),
    #     as.data.table(scores)
    # )

    scores <- as.data.table(scores, keep.rownames = 'gene')

    setkey(scores, gene)

    # Since genes_for_heatmap is in alphabetical order, ordering_genes is the ordering with
    # respect to the alphabetical order.

    if(seriate_genes_fun == 'seriate') {

        ordering_genes <- as.vector(
            seriation::seriate(
                as.matrix(scores[genes_for_heatmap, -'gene']),
                method = seriate_genes_method,
                margin = 1
            )[[1]]
        )

    } else { # Then seriate_genes_fun must be 'hclust'.

        # Recall the dist() function calculates distances between rows.

        ordering_genes <- hclust(
            dist(scores[genes_for_heatmap, -'gene']),
            method = seriate_genes_method
        )$order

    }

    # ordering_analyses is with respect to the ordering in ids_genes_ordering_list (Though
    # note that we put the analyses in order in the scores_melted table anyway, that is,
    # unique(scores_melted$analysis) is the same as
    # names(ids_genes_ordering_list)[ordering_analyses]).

    if(seriate_analyses_fun == 'seriate') {

        ordering_analyses <- as.vector(
            seriation::seriate(
                as.matrix(scores[genes_for_heatmap, -'gene']),
                method = seriate_analyses_method,
                margin = 2
            )[[1]]
        )

    } else { # Then seriate_analyses_fun must be 'hclust'.

        # Recall the dist() function calculates distances between rows.

        ordering_analyses <- hclust(
            dist(t(scores[genes_for_heatmap, -'gene'])),
            method = seriate_analyses_method
        )$order

    }

    scores_melted <- melt(
        # scores[ordering_genes],
        scores,
        id.vars = 'gene',
        measure.vars = names(ids_genes_ordering_list)[ordering_analyses],
        variable.name = 'analysis',
        value.name = 'score',
        variable.factor = FALSE # This is REALLY important!  Otherwise the 'analysis' column
        # will be a factor, which messes up stuff down the line.
    )

    setkey(scores_melted, gene)

    # Create heatmap:

    heatmap_commonality <- ggplot(
        scores_melted[genes_for_heatmap[ordering_genes]],
        aes(
            x = factor(analysis, levels = unique(analysis)),
            y = factor(gene, levels = unique(gene))
        )
    ) +
        geom_raster(aes(fill = score)) +
        scale_x_discrete(
            expand = c(0, 0)#,
            # labels = analyses$analysis_name[ordering_analyses]
        ) +
        scale_y_discrete(
            expand = c(0, 0)
        ) +
        theme(
            axis.ticks = element_blank(),
            axis.text.x = element_text(
                angle = x_axis_text_angle,
                hjust = 1
            ),
            axis.title.y = element_blank(),
            axis.ticks.length = unit(0, 'pt'),
            plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), 'pt')
        ) +
        labs(
            x = '',
            # fill = latex2exp::TeX('$\\frac{\\rho_{mes} - \\rho_{CAF}}{2}$')
            fill = legend_title
        )

    if(score_type == 'correlation') {
        heatmap_commonality <- heatmap_commonality +
            scale_fill_gradientn(
                colours = rev(colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))(50)),
                limits = c(-1, 1) #This makes zero appear as white.
            )
    } else {
        heatmap_commonality <- heatmap_commonality +
            scale_fill_gradientn(
                colours = rev(colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))(50))
            )
    }

    # It might be useful to include in the output the orders vectors for the genes and cancer
    # types, so that we can plot other data aligned with the heatmap and in the same order.

    list(
        heatmap = heatmap_commonality,
        scores = scores_melted,
        genes = genes_for_heatmap,
        analyses = names(ids_genes_ordering_list),
        ordering_genes = ordering_genes,
        ordering_analyses = ordering_analyses
    )

}





commonality_heatmap_from_scores <- function(

    scores_data,
    x_var = 'analysis',
    y_var = 'gene',
    fill_var = 'score',
    reselect_genes = FALSE,
    ids_genes_ordering_list = NULL,
    two_sided = TRUE,
    quantile_prob = 0.25,
    n_head = 25,
    n_tail = n_head,
    genes_to_include = c(
        'SNAI1',
        'SNAI2',
        'ZEB1',
        'ZEB2',
        'TWIST1',
        'VIM'
    ),
    score_type = c('correlation', 'expression_level'),
    normalise_scores = TRUE,
    seriate_x_fun = c('hclust', 'seriate'),
    seriate_x_method = 'average',
    seriate_y_fun = c('hclust', 'seriate'),
    seriate_y_method = 'average',
    x_axis_text_angle = 55,
    legend_title = 'pEMT-CAF\nscore'

) {

    score_type <- match.arg(score_type)
    seriate_x_fun <- match.arg(seriate_x_fun)
    seriate_y_fun <- match.arg(seriate_y_fun)

    if(reselect_genes & is.null(ids_genes_ordering_list)) {

        warning('Must supply ids_genes_ordering_list in order to reselect genes.  Using all genes.')

        genes_for_heatmap <- unique(scores_data[[y_var]])

    } else if(reselect_genes) {

        gene_ranking <- data.table(
            gene_id = unique(
                scores_data$gene
            )
        )

        if(two_sided) {

            setkeyv(scores_data, y_var)

            gene_ranking[
                ,
                c('quantile_rank_emt', 'quantile_rank_caf') := as.list(
                    apply(
                        sapply(
                            names(ids_genes_ordering_list),
                            function(li_name) {

                                li <- ids_genes_ordering_list[[li_name]]

                                li_genes <- li[[pmatch('gene', names(li))]]

                                if(gene_id %in% li_genes) {
                                    scores_data[
                                        analysis == li_name
                                        ][
                                            li_genes
                                            ][
                                                order(-get(fill_var)),
                                                c(
                                                    which(get(y_var) == gene_id),
                                                    which(rev(get(y_var)) == gene_id)
                                                )
                                                ]/length(
                                                    li_genes
                                                )
                                } else {
                                    c(1, 1)
                                }

                            }
                        ),
                        1,
                        quantile,
                        quantile_prob
                    )
                ),
                by = gene_id
                ]

            # Take top n EMT and CAF genes based on this ranking and combine with the genes we want
            # to include unconditionally:

            genes_for_heatmap <- sort(
                unique(
                    c(
                        gene_ranking[order(quantile_rank_emt)][1:n_head, gene_id],
                        gene_ranking[order(quantile_rank_caf)][1:n_tail, gene_id],
                        genes_to_include
                    )
                )
            )

        } else {

            gene_ranking[
                ,
                quantile_rank := quantile(
                    sapply(
                        names(ids_genes_ordering_list),
                        function(li_name) {

                            li <- ids_genes_ordering_list[[li_name]]

                            li_genes <- li[[pmatch('gene', names(li))]]

                            if(gene_id %in% li_genes) {

                                scores_data[
                                    analysis == li_name
                                    ][
                                        li_genes
                                        ][
                                            order(-get(fill_var)),
                                            which(get(y_var) == gene_id)
                                            ]/length(
                                                li_genes
                                            )
                            } else {
                                1
                            }

                        }
                    ),
                    quantile_prob
                ),
                by = gene_id
                ]

            # Take top n genes based on this ranking and combine with the genes we want to
            # include unconditionally:

            genes_for_heatmap <- sort(
                unique(
                    c(
                        gene_ranking[order(quantile_rank)][1:n_head, gene_id],
                        genes_to_include
                    )
                )
            )

        }

    } else {

        genes_for_heatmap <- unique(scores_data[[y_var]])

    }

    if(normalise_scores) {
        scores_data[[fill_var]] <- scores_data[[fill_var]]/max(
            abs(
                scores_data[
                    get(y_var) %in% genes_for_heatmap,
                    ..fill_var
                    ]
            )
        )
    }

    # Cast the melted scores_data to get a matrix for use in e.g. hclust():

    cast_data <- dcast(
        scores_data[get(y_var) %in% genes_for_heatmap],
        get(y_var) ~ get(x_var),
        value.var = fill_var
    )

    names(cast_data)[1] <- y_var

    if(seriate_y_fun == 'seriate') {

        ordering_y <- as.vector(
            seriation::seriate(
                as.matrix(cast_data[, -..y_var]),
                method = seriate_y_method,
                margin = 1
            )[[1]]
        )

    } else { # Then seriate_genes_fun must be 'hclust'.

        # Recall the dist() function calculates distances between rows.

        ordering_y <- hclust(
            dist(cast_data[, -..y_var]),
            method = seriate_y_method
        )$order

    }

    if(seriate_x_fun == 'seriate') {

        ordering_x <- as.vector(
            seriation::seriate(
                as.matrix(cast_data[, -..y_var]),
                method = seriate_x_method,
                margin = 2
            )[[1]]
        )

    } else { # Then seriate_analyses_fun must be 'hclust'.

        # Recall the dist() function calculates distances between rows.

        ordering_x <- hclust(
            dist(t(cast_data[, -..y_var])),
            method = seriate_x_method
        )$order

    }

    setkeyv(scores_data, y_var)

    # Create heatmap:

    heatmap_commonality <- ggplot(
        scores_data[cast_data[, get(y_var)][ordering_y]],
        aes(
            x = factor(
                get(x_var),
                levels = names(cast_data[, -..y_var])[ordering_x]
            ),
            y = factor(
                get(y_var),
                levels = unique(get(y_var))
            )
        )
    ) +
        geom_raster(aes(fill = get(fill_var))) +
        scale_x_discrete(
            expand = c(0, 0)#,
            # labels = analyses$analysis_name[ordering_analyses]
        ) +
        scale_y_discrete(
            expand = c(0, 0)
        ) +
        theme(
            axis.ticks = element_blank(),
            axis.text.x = element_text(
                angle = x_axis_text_angle,
                hjust = 1
            ),
            axis.title.y = element_blank(),
            axis.ticks.length = unit(0, 'pt'),
            plot.margin = unit(c(5.5, 10, 5.5, 0), 'pt')
        ) +
        labs(
            x = '',
            # fill = latex2exp::TeX('$\\frac{\\rho_{mes} - \\rho_{CAF}}{2}$')
            fill = legend_title
        )

    if(score_type == 'correlation') {
        heatmap_commonality <- heatmap_commonality +
            scale_fill_gradientn(
                colours = rev(colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))(50)),
                limits = c(-1, 1) #This makes zero appear as white.
            )
    } else {
        heatmap_commonality <- heatmap_commonality +
            scale_fill_gradientn(
                colours = rev(colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))(50))
            )
    }

    # Output:

    list(
        heatmap = heatmap_commonality,
        scores_data = scores_data[, c(..y_var, ..x_var, ..fill_var)],
        x = names(cast_data[, -..y_var]),
        y = cast_data[[y_var]],
        ordering_x = ordering_x,
        ordering_y = ordering_y
    )

}





extra_genes_by_cor <- function(

    expression_data,
    ids_genes_ordering_list,
    additional_genes_list = NULL,
    n = 300,
    head_or_tail = c('head', 'tail'),
    head_tail_size = 20,
    transform_data = TRUE,
    genes_to_include = c(
        'SNAI1',
        'SNAI2',
        'ZEB1',
        'ZEB2',
        'TWIST1',
        'VIM'
    )

) {

    # Transforming the data with all genes takes a very long time...  Add an option to input
    # a restricted gene set to be used for the transformation, e.g. set a threshold for
    # minimum expression level, or use the top N most variable genes.

    head_or_tail <- match.arg(head_or_tail)

    for(li in ids_genes_ordering_list) {
        if(is.na(pmatch('order', names(li)))) {
            li$ordering <- 1:length(li[[pmatch('genes', names(li))]])
        }
    }

    sapply(

        names(ids_genes_ordering_list),

        function(li_name) {

            cat(li_name, '\b...')

            li_sample_ids <- ids_genes_ordering_list[[
                li_name
                ]][[
                    pmatch(
                        'sample',
                        names(ids_genes_ordering_list[[li_name]])
                    )
                    ]]

            li_ordered_genes <- ids_genes_ordering_list[[
                li_name
                ]][[
                    pmatch(
                        'gene',
                        names(ids_genes_ordering_list[[li_name]])
                    )
                    ]][
                        ids_genes_ordering_list[[
                            li_name
                            ]][[
                                pmatch(
                                    'order',
                                    names(ids_genes_ordering_list[[li_name]])
                                )
                                ]]
                        ]

            marker_genes <- match.fun(head_or_tail)(
                li_ordered_genes,
                head_tail_size
            )

            if(transform_data) {

                if(!is.null(additional_genes_list)) {

                    genes_to_transform <- c(
                        additional_genes_list[[li_name]],
                        marker_genes
                    )

                } else {

                    genes_to_transform <- names(expression_data[, -'id'])

                }

                scores_data <- cbind(
                    expression_data[
                        li_sample_ids,
                        ..genes_to_transform
                        ],
                    row_sums = rowSums(
                        expression_data[
                            li_sample_ids,
                            ..li_ordered_genes
                            ]
                    )
                )

                scores_data <- as.data.table(
                    sapply(
                        genes_to_transform,
                        function(g) {

                            # In the following, I have to put the backticks `` around g because
                            # there are genes like 'NKX2-1', which gets interpreted as NKX2 - 1.

                            lm(
                                formula(
                                    paste0(
                                        '`',
                                        g,
                                        '` ~ row_sums'
                                    )
                                ),
                                data = scores_data
                            )$residuals

                        },
                        simplify = FALSE,
                        USE.NAMES = TRUE
                    )
                )

            } else {

                scores_data <- expression_data[li_sample_ids]

            }

            cors <- cor(
                scores_data[
                    ,
                    -c('id', ..marker_genes)
                    ],
                rowMeans(
                    scores_data[
                        ,
                        ..marker_genes
                        ]
                )
            )

            cors <- cors[!is.na(cors[, 1]), 1]

            output_genes <- unique(
                c(
                    genes_to_include,
                    names(sort(cors, decreasing = TRUE))[1:n]
                )
            )

            cat('Done!\n')

            output_genes

        },

        simplify = FALSE,
        USE.NAMES = TRUE

    )

}





cor_with_clinical <- function(

    expression_data,
    ids_genes_ordering_list,
    clinical_data,
    clin_var,
    wilcox_test_x_expr,
    wilcox_test_y_expr,
    scores_data = NULL,
    id_var = 'id',
    emt_score_var = 'emt_score',
    caf_score_var = 'caf_score',
    genes_for_scores = NULL,
    amatch_max_dist = 10,
    p_adjust_method = 'BH',
    signif_threshold = 0.05,
    adjust_pvals = FALSE,
    only_signif_labels = TRUE,
    show_adjusted_pval_threshold = TRUE,
    plot_title = ''

) {

    # The following was intended for using the function in do.call with quote = TRUE, so I could
    # dynamically pass conditions for the Wilcoxon tests in lists.  But it didn't work.

    # fun_args <- as.list(match.call())
    #
    # for(arg_name in names(fun_args)) {
    #     assign(arg_name, eval(fun_args[[arg_name]]))
    # }

    # I added the following so you don't have to pass a quoted condition in the function arguments:

    wilcox_test_x_expr <- substitute(wilcox_test_x_expr)
    wilcox_test_y_expr <- substitute(wilcox_test_y_expr)

    setkey(expression_data, id)
    setkey(clinical_data, id)

    clinical_analyses <- sapply(

        names(ids_genes_ordering_list),

        function(li_name) {

            cat(li_name, '\b...')

            # Extract patient IDs from the sample IDs given in ids_genes_ordering_list (though
            # this should still work if the IDs given are already patient IDs):

            li <- ids_genes_ordering_list[[li_name]]

            sample_ids <- li[[pmatch('sample', names(li))]]

            patient_ids <- apply(
                stringr::str_split_fixed(
                    sample_ids,
                    '\\.',
                    4
                )[, 1:3],
                1,
                paste,
                collapse = '.'
            )

            # Look for a match for clin_var, with specified max. distance, among the variables
            # that are not all NA or the empty string for this cancer type/subtype, and return
            # NA if we don't find one:

            non_na_names <- names(
                clinical_data[
                    patient_ids,
                    sapply(
                        .SD,
                        function(x) {
                            switch((sum(!is.na(x) & x != '') == 0) + 1, x, NULL)
                        }
                    )
                ]
            )

            li_clin_var <- non_na_names[
                stringdist::amatch(
                    clin_var,
                    non_na_names,
                    maxDist = amatch_max_dist
                )
            ]

            if(is.na(li_clin_var)) {
                cat('No variable match.\n')
                return(NA)
            }

            # Restrict sample and patient IDs to those where the clinical variable is not NA:

            sample_ids <- sample_ids[
                clinical_data[
                    patient_ids,
                    !is.na(get(li_clin_var)) & get(li_clin_var) != ''
                ]
            ]

            patient_ids <- patient_ids[
                clinical_data[
                    patient_ids,
                    !is.na(get(li_clin_var)) & get(li_clin_var) != ''
                ]
            ]

            if(!is.null(scores_data)) {

                setkeyv(scores_data, id_var)

                data_for_test <- clinical_data[
                    patient_ids,
                    .(
                        id = id,
                        variable = get(li_clin_var),
                        emt_score = scores_data[
                            sample_ids,
                            ..emt_score_var
                        ],
                        caf_score = scores_data[
                            sample_ids,
                            ..caf_score_var
                        ]
                    )
                ]

            } else {

                if(!is.null(genes_for_scores)) {

                    emt_score_genes <- genes_for_scores[[li_name]]$emt
                    caf_score_genes <- genes_for_scores[[li_name]]$caf

                } else {

                    emt_score_genes <- head(
                        li[[pmatch('gene', names(li))]][
                            li[[pmatch('order', names(li))]]
                        ],
                        20
                    )

                    caf_score_genes <- tail(
                        li[[pmatch('gene', names(li))]][
                            li[[pmatch('order', names(li))]]
                        ],
                        20
                    )

                }

                data_for_test <- clinical_data[
                    patient_ids,
                    .(
                        id = id,
                        variable = get(li_clin_var),
                        emt_score = rowMeans(
                            expression_data[
                                sample_ids,
                                ..emt_score_genes
                            ]
                        ),
                        caf_score = rowMeans(
                            expression_data[
                                sample_ids,
                                ..caf_score_genes
                            ]
                        )
                    )
                ]

            }

            # Check that we have enough observations for the Wilcoxon test, and return NA
            # if not:

            criterion <- data_for_test[
                ,
                c(
                    sum(eval(wilcox_test_x_expr)) == 0,
                    sum(eval(wilcox_test_y_expr)) == 0
                )
            ]

            if(sum(criterion) > 0) {
                cat(
                    'Not enough',
                    paste(c('x', 'y')[criterion], collapse = ' or '),
                    'observations.\n'
                )
                return(NA)
            }

            # Running Wilcoxon rank-sum test:

            data_for_test <- data_for_test[
                ,
                setNames(
                    c(
                        li_clin_var,
                        wilcox.test(
                            .SD[eval(wilcox_test_x_expr), emt_score],
                            .SD[eval(wilcox_test_y_expr), emt_score]
                        )[c('statistic', 'p.value')],
                        wilcox.test(
                            .SD[eval(wilcox_test_x_expr), caf_score],
                            .SD[eval(wilcox_test_y_expr), caf_score]
                        )[c('statistic', 'p.value')]
                    ),
                    c('variable_name', 'stat_emt', 'pval_emt', 'stat_caf', 'pval_caf')
                )
            ]

            cat('Done!\n')

            data_for_test

        },

        simplify = FALSE,
        USE.NAMES = TRUE

    )

    # Bind together into a single data table:

    clinical_analyses <- rbindlist(
        lapply(
            names(clinical_analyses),
            function(an) {
                if(!is.na(clinical_analyses[an])) {
                    cbind(
                        test_name = an,
                        clinical_analyses[[an]]
                    )
                }
            }
        )
    )

    clinical_analyses[
        ,
        c('pval_emt_adj', 'pval_caf_adj') := .(
            p.adjust(pval_emt, p_adjust_method),
            p.adjust(pval_caf, p_adjust_method)
        )
    ]

    setcolorder(
        clinical_analyses,
        c(
            'test_name',
            'variable_name',
            'stat_emt',
            'pval_emt',
            'pval_emt_adj',
            'stat_caf',
            'pval_caf',
            'pval_caf_adj'
        )
    )





    # Make scatterplot:

    if(only_signif_labels) {
        test_labels <- clinical_analyses[
            ,
            .(
                test_labels = switch(
                    (pval_emt <= signif_threshold | pval_caf <= signif_threshold) + 1,
                    '',
                    test_name
                )
            ),
            by = test_name
        ]$test_labels
    } else {
        test_labels <- clinical_analyses$test_name
    }

    g <- ggplot(
        clinical_analyses
    ) +
        geom_point(
            aes(
                x = -log10(switch(adjust_pvals + 1, pval_caf, pval_caf_adj)),
                y = -log10(switch(adjust_pvals + 1, pval_emt, pval_emt_adj))
            ),
            shape = 17,
            size = 2.5
        ) +
        ggrepel::geom_label_repel(
            aes(
                x = -log10(switch(adjust_pvals + 1, pval_caf, pval_caf_adj)),
                y = -log10(switch(adjust_pvals + 1, pval_emt, pval_emt_adj)),
                label = test_labels
            )
        ) +
        geom_vline(
            xintercept = -log10(signif_threshold),
            colour = 'darkorange',
            linetype = 'dashed'
        ) +
        geom_hline(
            yintercept = -log10(signif_threshold),
            colour = 'darkorange',
            linetype = 'dashed'
        )

    if(!adjust_pvals & show_adjusted_pval_threshold & p_adjust_method == 'BH') {

        bh_adjusted_signif_threshold_caf <- clinical_analyses[
            order(pval_caf),
            switch(
                (sum(pval_caf <= signif_threshold*.I/.N) == 0) + 1,
                signif_threshold*max(
                    which(
                        pval_caf <= signif_threshold*.I/.N
                    )
                )/.N,
                NA
            )
        ]

        bh_adjusted_signif_threshold_emt <- clinical_analyses[
            order(pval_emt),
            switch(
                (sum(pval_emt <= signif_threshold*.I/.N) == 0) + 1,
                signif_threshold*max(
                    which(
                        pval_emt <= signif_threshold*.I/.N
                    )
                )/.N,
                NA
            )
        ]

        if(!is.na(bh_adjusted_signif_threshold_caf)) {
            g <- g + geom_vline(
                xintercept = -log10(bh_adjusted_signif_threshold_caf),
                colour = 'red3',
                linetype = 'dashed'
            )
        }

        if(!is.na(bh_adjusted_signif_threshold_emt)) {
            g <- g + geom_hline(
                yintercept = -log10(bh_adjusted_signif_threshold_emt),
                colour = 'red3',
                linetype = 'dashed'
            )
        }

    }

    if(adjust_pvals) {

        g <- g + labs(
            # x = latex2exp::TeX('$-\\log_{10}(p^{adj}_{CAF})$'),
            # y = latex2exp::TeX('$-\\log_{10}(p^{adj}_{EMT})$'),
            # What follows is the output from latex2exp::TeX, put into expression(atop()).
            # It is very messy, but despite RStudio's warnings it does work.  The reason for
            # using this is that neither \n nor \newline work if used directly in latex2exp::TeX.
            x = expression(
                atop(
                    `\textbf{Significance of association for CAFs}` = paste(
                        "",
                        bold(paste("Significance of association for CAFs"))
                    ),
                    `$-log_{10}(p^{adj}_{CAF})$` = paste(
                        "",
                        phantom() - phantom(),
                        log,
                        ,
                        ,
                        ,
                        phantom()[{paste("10")}],
                        "(",
                        "",
                        "p",
                        phantom()[paste("", "CAF")]^{paste("adj")},
                        ")",
                        "",
                        ""
                    )
                )
            ),
            y = expression(
                atop(
                    `\textbf{Significance of association for EMT}` = paste(
                        "",
                        bold(paste("Significance of association for EMT"))
                    ),
                    `$-log_{10}(p^{adj}_{EMT})$` = paste(
                        "",
                        phantom() - phantom(),
                        log,
                        ,
                        ,
                        ,
                        phantom()[{paste("10")}],
                        "(",
                        "",
                        "p",
                        phantom()[paste("", "EMT")]^{paste("adj")},
                        ")",
                        "",
                        ""
                    )
                )
            ),
            title = plot_title
        )

    } else {

        g <- g + labs(
            # x = latex2exp::TeX('$-\\log_{10}(p^{adj}_{CAF})$'),
            # y = latex2exp::TeX('$-\\log_{10}(p^{adj}_{EMT})$'),
            # What follows is the output from latex2exp::TeX, put into expression(atop()).
            # It is very messy, but despite RStudio's warnings it does work.  The reason for
            # using this is that neither \n nor \newline work if used directly in latex2exp::TeX.
            x = expression(
                atop(
                    `\textbf{Significance of association for CAFs}` = paste(
                        "",
                        bold(paste("Significance of association for CAFs"))
                    ),
                    `$-log_{10}(p_{CAF})$` = paste(
                        "",
                        "-log",
                        phantom()[{paste("10")}],
                        "(",
                        "",
                        "p",
                        phantom()[{paste("CAF")}],
                        ")",
                        "",
                        ""
                    )
                )
            ),
            y = expression(
                atop(
                    `\textbf{Significance of association for EMT}` = paste(
                        "",
                        bold(paste("Significance of association for EMT"))
                    ),
                    `$-log_{10}(p_{EMT})$` = paste(
                        "",
                        "-log",
                        phantom()[{paste("10")}],
                        "(",
                        "",
                        "p",
                        phantom()[{paste("EMT")}],
                        ")",
                        "",
                        ""
                    )
                )
            ),
            title = plot_title
        )

    }

    # Output:

    list(
        plot = g,
        test_results = clinical_analyses
    )

}





volcano_clinical <- function(

    expression_data,
    sample_ids_list,
    genes,
    clinical_data,
    clin_var,
    wilcox_test_x_expr,
    wilcox_test_y_expr,
    amatch_max_dist = 10,
    p_adjust_method = c('BH', 'bonferroni'),
    signif_threshold = 0.05,
    labels_signif_threshold = 'non-adjusted',
    show_pval_thresholds = c('non-adjusted', 'adjusted'),
    fold_change_fun = function(fc) log2(fc + 1),
    xlab_TeX = 'log_2(fold change + 1)',
    ylab_TeX = '-log_{10}(p value)',
    plot_title = ''

) {

    # <sample_ids_list> should have names!

    # <genes> can either be a character vector or a list of the same length as <sample_ids_list>
    # (with the same names) with each element being a character vector.  If the former, the same
    # genes are used for all the elements of <sample_ids_list>; if the latter, the genes in the
    # n-th element of <genes> are used for the n-the element of <sample_ids_list>.

    # Similarly, <expression_data> and <clinical_data> can either be single data tables or lists
    # of data tables (with the same length names as <sample_ids_list>), depending on whether you
    # want to use the same table for all cancer types or a different one for each.  If
    # <sample_ids_list> is not given and <expression_data> is a list, the sample IDs will be
    # taken from the ids in <expression_data>.

    # <show_pval_threshold> can take a character vector of length > 1, namely
    # c('non-adjusted', 'adjusted') in which case both non-adjusted and adjusted significance
    # thresholds will be shown on the plot.  <labels_signif_threshold> can only be a character
    # vector of length 1.

    if(missing(sample_ids_list) & !is.data.frame(expression_data)) {
        sample_ids_list <- lapply(expression_data, `[[`, 'id')
    }

    p_adjust_method <- match.arg(p_adjust_method)

    labels_signif_threshold <- match.arg(
        labels_signif_threshold,
        choices = c('non-adjusted', 'adjusted', 'none'),
        several.ok = FALSE
    )

    show_pval_thresholds <- match.arg(
        show_pval_thresholds,
        choices = c('non-adjusted', 'adjusted', 'none'),
        several.ok = TRUE
    )

    fold_change_fun <- match.fun(fold_change_fun)

    wilcox_test_x_expr <- substitute(wilcox_test_x_expr)
    wilcox_test_y_expr <- substitute(wilcox_test_y_expr)

    clinical_analyses <- sapply(

        names(sample_ids_list),

        function(ct) {

            cat(paste0(ct, '...'))

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

            patient_ids <- apply(
                stringr::str_split_fixed(
                    sample_ids,
                    '\\.',
                    4
                )[, 1:3],
                1,
                paste,
                collapse = '.'
            )

            # Look for a match for clin_var, with specified max. distance, among the variables
            # that are not all NA or the empty string for this cancer type/subtype, and return
            # NA if we don't find one:

            non_na_names <- names(
                clinical_data_table[
                    patient_ids,
                    sapply(
                        .SD,
                        function(x) {
                            switch((sum(!is.na(x) & x != '') == 0) + 1, x, NULL)
                        }
                    )
                ]
            )

            li_clin_var <- non_na_names[
                stringdist::amatch(
                    clin_var,
                    non_na_names,
                    maxDist = amatch_max_dist
                )
            ]

            if(is.na(li_clin_var)) {
                cat('No variable match.\n')
                return(NA)
            }

            # Restrict sample and patient IDs to those where the clinical variable is not NA:

            sample_ids <- sample_ids[
                clinical_data_table[
                    patient_ids,
                    !is.na(get(li_clin_var)) & get(li_clin_var) != ''
                ]
            ]

            patient_ids <- patient_ids[
                clinical_data_table[
                    patient_ids,
                    !is.na(get(li_clin_var)) & get(li_clin_var) != ''
                ]
            ]

            # Calculate scores by average expression:

            data_for_test <- clinical_data_table[
                patient_ids,
                .(
                    id = id,
                    variable = get(li_clin_var),
                    score = rowMeans(
                        expression_data_table[
                            sample_ids,
                            ..genes_vec
                        ]
                    )
                )
            ]

            # Check that we have enough observations for the Wilcoxon test, and return NA
            # if not:

            criterion <- data_for_test[
                ,
                c(
                    sum(eval(wilcox_test_x_expr)) == 0,
                    sum(eval(wilcox_test_y_expr)) == 0
                )
            ]

            if(sum(criterion) > 0) {
                cat(
                    'Not enough',
                    paste(c('x', 'y')[criterion], collapse = ' or '),
                    'observations.\n'
                )
                return(NA)
            }

            # Running Wilcoxon rank-sum test and calculating fold change:

            out <- data_for_test[
                ,
                c(
                    variable_name = li_clin_var,
                    setNames(
                        wilcox.test(
                            .SD[eval(wilcox_test_x_expr), score],
                            .SD[eval(wilcox_test_y_expr), score]
                        )[c('statistic', 'p.value')],
                        c('stat', 'pval')
                    ),
                    fold_change = .SD[eval(wilcox_test_x_expr), mean(score)]/
                        .SD[eval(wilcox_test_y_expr), mean(score)] - 1
                )
            ]

            cat('Done!\n')

            out

        },

        simplify = FALSE,
        USE.NAMES = TRUE

    )

    # Bind together into a single data table:

    clinical_analyses <- rbindlist(
        lapply(
            names(clinical_analyses),
            function(an) {
                if(!is.na(clinical_analyses[an])) {
                    cbind(
                        test_name = an,
                        clinical_analyses[[an]]
                    )
                }
            }
        )
    )

    clinical_analyses[, pval_adj := p.adjust(pval, p_adjust_method)]

    setcolorder(
        clinical_analyses,
        c(
            'test_name',
            'variable_name',
            'stat',
            'pval',
            'pval_adj',
            'fold_change'
        )
    )

    # Optionally calculate adjusted significance threshold:

    names(signif_threshold) <- 'non-adjusted'

    if('adjusted' %in% labels_signif_threshold | 'adjusted' %in% show_pval_thresholds) {

        if(p_adjust_method == 'BH') {

            adjusted_signif_threshold <- clinical_analyses[
                order(pval),
                switch(
                    (sum(pval <= signif_threshold*.I/.N) == 0) + 1,
                    signif_threshold*max(
                        which(
                            pval <= signif_threshold*.I/.N
                        )
                    )/.N,
                    NA
                )
            ]

        } else { # Then p_adjust_method == 'bonferroni'

            adjusted_signif_threshold <- signif_threshold/nrow(clinical_analyses)

        }

        signif_threshold <- c(
            signif_threshold,
            setNames(adjusted_signif_threshold, 'adjusted')
        )

    }

    # Vector of labels for points:

    if(labels_signif_threshold != 'none') {

        test_labels <- clinical_analyses[
            ,
            .(
                test_labels = switch(
                    (pval <= signif_threshold[labels_signif_threshold]) + 1,
                    '',
                    test_name
                )
            ),
            by = test_name
        ]$test_labels

    } else {

        test_labels <- clinical_analyses$test_name

    }

    # Make volcano plot:

    g <- ggplot(
        clinical_analyses
    ) +
        geom_point(
            aes(
                x = fold_change_fun(fold_change),
                y = -log10(pval)
            ),
            shape = 17,
            size = 2.5
        ) +
        ggrepel::geom_label_repel(
            aes(
                x = fold_change_fun(fold_change),
                y = -log10(pval),
                label = test_labels
            )
        ) +
        geom_hline(
            aes(
                yintercept = value,
                colour = threshold
            ),
            data = data.frame(
                threshold = names(signif_threshold),
                value = -log10(signif_threshold)
            ),
            linetype = 'dashed'
        ) +
        scale_colour_manual(
            values = c('non-adjusted' = 'darkorange', 'adjusted' = 'red3')
        ) +
        labs(
            x = latex2exp::TeX(xlab_TeX),
            y = latex2exp::TeX(ylab_TeX),
            title = plot_title
        ) +
        theme_test()

    # Output:

    cat('\n')

    list(
        plot = g,
        test_results = clinical_analyses
    )

}





volcano_clinical_list <- function(

    expression_data,
    sample_ids_list,
    genes,
    clinical_data,
    clin_var,
    wilcox_test_x_expr,
    wilcox_test_y_expr,
    amatch_max_dist = 10,
    p_adjust_method = c('BH', 'bonferroni'),
    signif_threshold = 0.05,
    labels_signif_threshold = 'non-adjusted',
    show_pval_thresholds = c('non-adjusted', 'adjusted'),
    fold_change_fun = function(fc) log2(fc + 1),
    xlab_TeX = 'log_2(fold change + 1)',
    ylab_TeX = '-log_{10}(p value)',
    legend_labels = NULL,
    legend_colours = NULL,
    legend_title = waiver(),
    plot_title = waiver()

) {

    # <sample_ids_list> should have names!

    # <genes> can either be a character vector or a list of the same length as <sample_ids_list>
    # (with the same names) with each element being a character vector.  If the former, the same
    # genes are used for all the elements of <sample_ids_list>; if the latter, the genes in the
    # n-th element of <genes> are used for the n-the element of <sample_ids_list>.

    # Similarly, <expression_data> and <clinical_data> can either be single data tables or lists
    # of data tables (with the same length names as <sample_ids_list>), depending on whether you
    # want to use the same table for all cancer types or a different one for each.  If
    # <sample_ids_list> is not given and <expression_data> is a list, the sample IDs will be
    # taken from the ids in <expression_data>.

    # <show_pval_threshold> can take a character vector of length > 1, namely
    # c('non-adjusted', 'adjusted') in which case both non-adjusted and adjusted significance
    # thresholds will be shown on the plot.  <labels_signif_threshold> can only be a character
    # vector of length 1.

    if(missing(sample_ids_list) & !is.data.frame(expression_data)) {
        sample_ids_list <- lapply(expression_data, `[[`, 'id')
    }

    p_adjust_method <- match.arg(p_adjust_method)

    labels_signif_threshold <- match.arg(
        labels_signif_threshold,
        choices = c('non-adjusted', 'adjusted', 'none'),
        several.ok = FALSE
    )

    show_pval_thresholds <- match.arg(
        show_pval_thresholds,
        choices = c('non-adjusted', 'adjusted', 'none'),
        several.ok = TRUE
    )

    fold_change_fun <- match.fun(fold_change_fun)

    # wilcox_test_x_expr <- lapply(wilcox_test_x_expr, substitute)
    # wilcox_test_y_expr <- lapply(wilcox_test_y_expr, substitute)

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

            patient_ids <- apply(
                stringr::str_split_fixed(
                    sample_ids,
                    '\\.',
                    4
                )[, 1:3],
                1,
                paste,
                collapse = '.'
            )

            # Look for a match for clin_var, with specified max. distance, among the variables
            # that are not all NA or the empty string for this cancer type/subtype, and return
            # NA if we don't find one:

            non_na_names <- names(
                clinical_data_table[
                    patient_ids,
                    sapply(
                        .SD,
                        function(x) {
                            switch((sum(!is.na(x) & x != '') == 0) + 1, x, NULL)
                        }
                    )
                ]
            )

            # The following returns a list with names being the non-NA variable name matches and
            # each element being a list consisting of the corresponding x and y Wilcoxon test
            # expressions along with the input variable name.

            ct_clin_vars <- unlist(
                lapply(
                    clin_var,
                    function(v) {

                        var_match <- non_na_names[
                            stringdist::amatch(
                                v,
                                non_na_names,
                                maxDist = amatch_max_dist
                            )
                        ]

                        if(is.na(var_match)) {
                            cat("No variable match for ", v, ".\n", sep = '')
                        } else {
                            cat("Found match for ", v, ": ", var_match, ".\n", sep = '')
                            return(
                                setNames(
                                    list( # List whose one element is a list with 3 elements
                                        list(
                                            variable_name = v,
                                            x_expr = wilcox_test_x_expr[clin_var == v][[1]],
                                            y_expr = wilcox_test_y_expr[clin_var == v][[1]]
                                        )
                                    ),
                                    var_match # Name for the single element of the 'outer' list
                                )
                            )
                        }

                    }
                ),
                recursive = FALSE
            )

            test_results <- rbindlist(
                lapply(
                    names(ct_clin_vars),
                    function(v) {

                        # Restrict sample and patient IDs to those where the clinical variable
                        # is not NA:

                        sample_ids <- sample_ids[
                            clinical_data_table[
                                patient_ids,
                                !is.na(get(v)) & get(v) != ''
                            ]
                        ]

                        patient_ids <- patient_ids[
                            clinical_data_table[
                                patient_ids,
                                !is.na(get(v)) & get(v) != ''
                            ]
                        ]

                        # Calculate scores by average expression:

                        data_for_test <- clinical_data_table[
                            patient_ids,
                            .(
                                id = id,
                                variable = get(v),
                                score = rowMeans(
                                    expression_data_table[
                                        sample_ids,
                                        ..genes_vec
                                    ]
                                )
                            )
                        ]

                        # Check that we have enough observations for the Wilcoxon test, and
                        # return NA if not:

                        criterion <- data_for_test[
                            ,
                            c(
                                sum(eval(ct_clin_vars[[v]]$x_expr)) == 0,
                                sum(eval(ct_clin_vars[[v]]$y_expr)) == 0
                            )
                        ]

                        if(sum(criterion) > 0) {

                            cat(
                                "Not enough ",
                                paste(c('x', 'y')[criterion], collapse = ' or '),
                                " observations for ",
                                v,
                                ".\n",
                                sep = ''
                            )

                            # return(NA)

                        } else {

                            # Running Wilcoxon rank-sum test and calculating fold change:

                            return(
                                data_for_test[
                                    ,
                                    c(
                                        test_name = ct,
                                        variable_name = ct_clin_vars[[v]]$variable_name,
                                        variable_match = v,
                                        setNames(
                                            wilcox.test(
                                                .SD[eval(ct_clin_vars[[v]]$x_expr), score],
                                                .SD[eval(ct_clin_vars[[v]]$y_expr), score]
                                            )[c('statistic', 'p.value')],
                                            c('stat', 'pval')
                                        ),
                                        fold_change = .SD[
                                            eval(ct_clin_vars[[v]]$x_expr),
                                            mean(score)
                                        ]/.SD[
                                            eval(ct_clin_vars[[v]]$y_expr),
                                            mean(score)
                                        ] - 1
                                    )
                                ]
                            )

                        }

                    }
                )
            )

            cat('Done!\n')

            test_results

        },

        simplify = FALSE,
        USE.NAMES = TRUE

    )

    # Bind together into a single data table:

    # clinical_analyses <- rbindlist(
    #     lapply(
    #         names(clinical_analyses),
    #         function(an) {
    #             if(!is.na(clinical_analyses[an])) {
    #                 cbind(
    #                     test_name = an,
    #                     clinical_analyses[[an]]
    #                 )
    #             }
    #         }
    #     )
    # )

    clinical_analyses <- rbindlist(clinical_analyses)

    clinical_analyses[
        ,
        c('test_id', 'pval_adj') := .(
            .I,
            p.adjust(pval, p_adjust_method)
        )
    ]

    setcolorder(
        clinical_analyses,
        c(
            'test_id',
            'test_name',
            'variable_name',
            'variable_match',
            'stat',
            'pval',
            'pval_adj',
            'fold_change'
        )
    )

    # Optionally calculate adjusted significance threshold:

    names(signif_threshold) <- 'non-adjusted'

    if('adjusted' %in% labels_signif_threshold | 'adjusted' %in% show_pval_thresholds) {

        if(p_adjust_method == 'BH') {

            adjusted_signif_threshold <- clinical_analyses[
                order(pval),
                switch(
                    (sum(pval <= signif_threshold*.I/.N) == 0) + 1,
                    signif_threshold*max(
                        which(
                            pval <= signif_threshold*.I/.N
                        )
                    )/.N,
                    NA
                )
            ]

        } else { # Then p_adjust_method == 'bonferroni'

            adjusted_signif_threshold <- signif_threshold/nrow(clinical_analyses)

        }

        signif_threshold <- c(
            signif_threshold,
            setNames(adjusted_signif_threshold, 'adjusted')
        )

    }

    # Vector of labels for points:

    if(labels_signif_threshold != 'none') {

        test_labels <- clinical_analyses[
            ,
            .(
                test_labels = switch(
                    (pval <= signif_threshold[labels_signif_threshold]) + 1,
                    '',
                    test_name
                )
            ),
            by = test_id
        ]$test_labels

    } else {

        test_labels <- clinical_analyses$test_name

    }

    # Make volcano plot:

    g <- ggplot(
        clinical_analyses
    ) +
        geom_point(
            aes(
                x = fold_change_fun(fold_change),
                y = -log10(pval),
                colour = variable_name
            ),
            shape = 17,
            size = 2.5
        ) +
        ggrepel::geom_label_repel(
            aes(
                x = fold_change_fun(fold_change),
                y = -log10(pval),
                label = test_labels
            )
        ) +
        # geom_hline(
        #     aes(
        #         yintercept = value,
        #         colour = threshold
        #     ),
        #     data = data.frame(
        #         threshold = names(signif_threshold),
        #         value = -log10(signif_threshold)
        #     ),
        #     linetype = 'dashed'
        # ) +
        # scale_colour_manual(
        #     values = c('non-adjusted' = 'darkorange', 'adjusted' = 'red3')
        # ) +
        scale_colour_manual(
            labels = switch(
                is.null(legend_labels) + 1,
                legend_labels,
                sort(clin_var)
            ),
            values = switch(
                is.null(legend_colours) + 1,
                legend_colours,
                setNames(
                    hcl(
                        h = seq(
                            360/(2*length(clin_var)),
                            360 - 360/(2*length(clin_var)),
                            length.out = length(clin_var)
                        ),
                        c = 100,
                        l = 65
                    ),
                    sort(clin_var)
                )
            )
        ) +
        labs(
            x = latex2exp::TeX(xlab_TeX),
            y = latex2exp::TeX(ylab_TeX),
            colour = legend_title,
            title = plot_title
        ) +
        theme_test()

    if('non-adjusted' %in% names(signif_threshold)) {
        g <- g + geom_hline(
            aes(yintercept = -log10(signif_threshold['non-adjusted'])),
            colour = 'darkorange',
            linetype = 'dashed'
        )
    }

    if('adjusted' %in% names(signif_threshold)) {
        g <- g + geom_hline(
            aes(yintercept = -log10(signif_threshold['adjusted'])),
            colour = 'red3',
            linetype = 'dashed'
        )
    }

    # Output:

    cat('\n')

    list(
        plot = g,
        test_results = clinical_analyses
    )

}





clinical_test <- function(

    expression_data,
    sample_ids_list,
    genes,
    clinical_data,
    clin_var,
    wilcox_test_x_expr,
    wilcox_test_y_expr,
    amatch_max_dist = 10,
    p_adjust_method = c('BH', 'bonferroni'),
    min_samples = 1,
    score_method = c('average', 'scrabble'),
    scores_filter_fun = NULL, # function(x) {mean(x) >= 2} is probably a good alternative
    genes_for_data_transform = NULL

) {

    # <sample_ids_list> should have names!

    # <genes> can either be a character vector or a list of the same length as <sample_ids_list>
    # (with the same names) with each element being a character vector.  If the former, the same
    # genes are used for all the elements of <sample_ids_list>; if the latter, the genes in the
    # n-th element of <genes> are used for the n-the element of <sample_ids_list>.

    # Similarly, <expression_data> and <clinical_data> can either be single data tables or lists
    # of data tables (with the same length names as <sample_ids_list>), depending on whether you
    # want to use the same table for all cancer types or a different one for each.  If
    # <sample_ids_list> is not given and <expression_data> is a list, the sample IDs will be
    # taken from the ids in <expression_data>.

    # Each of <wilcox_test_x_expr> and <wilcox_test_y_expr> can be either a list or a single
    # expression.  If the latter, the same expression is used for every clinical variable; if the
    # former, the n-th expression will be used for the n-th variable.

    # By default, matches for the elements of <clin_var> are searched for using amatch() from the
    # stringdist package, with maxDist equal to 10.  The <amatch_max_dist> argument is passed to
    # the maxDist parameter of amatch(), so can be used to specify how closely you need the
    # variable names in the data to match those in <clin_var>.  If exact matches are required,
    # set <amatch_max_dist> to NULL.

    # <score_method> is the choice of method for assigning scores to each sample based on <genes>.
    # There are two options: 'average' simply calculates the average expression of the genes in
    # <genes> for each sample; 'scrabble' uses the score() function from the scrabble package
    # with bins proportional to the number of genes in the subset of <expression_data> defined
    # by <sample_ids> and any filtering of <genes>, which can be specified under
    # <scores_filter_fun>.

    if(missing(sample_ids_list) & !is.data.frame(expression_data)) {
        sample_ids_list <- lapply(expression_data, `[[`, 'id')
    }
    p_adjust_method <- match.arg(p_adjust_method)
    score_method <- match.arg(score_method)

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

            # The following is for the (optional) data transform by calculating residuals.  I do this before refining sample_ids because I want
            # the transformation to be the same as in deconv_data.
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

            # Match IDs from clinical data and expression data on the level of patients,
            # i.e. first extract the first 3 components of each element in sample_ids
            # and then look for matches for these patient IDs in clinical_data_table$id.
            # We exclude non-unique matches.

            # patient_ids <- apply(stringr::str_split_fixed(sample_ids, '\\.', 4)[, 1:3], 1, paste, collapse = '.')

            all_ids <- data.table(
                expression_data_id = sample_ids,
                patient_id = apply(stringr::str_split_fixed(sample_ids, '\\.', 4)[, 1:3], 1, paste, collapse = '.')
            )[
                ,
                clinical_data_id := sapply(
                    patient_id,
                    function(x) {
                        matches <- clinical_data_table[, id[grep(paste0('^', x), id)]]
                        if(length(matches) == 1) {
                            return(matches)
                        } else {
                            return(NA)
                        }
                    }
                )
            ][!(clinical_data_id %in% names(table(clinical_data_id))[table(clinical_data_id) > 1])]

            # The following two lines just allow consistency with the rest of the code for this function, which I wrote before the addition of all_ids.
            patient_ids <- all_ids$clinical_data_id
            sample_ids <- all_ids$expression_data_id

            # Look for a match for clin_var, with specified max. distance, among the variables that are not all NA or the empty string for this cancer
            # type/subtype, and return NA if we don't find one:
            non_na_names <- colnames(clinical_data_table[patient_ids, sapply(.SD, function(x) switch((sum(!is.na(x) & x != '') == 0) + 1, x, NULL))])

            # The following returns a list with names being the non-NA variable name matches and
            # each element being a list consisting of the corresponding x and y Wilcoxon test
            # expressions along with the input variable name.

            ct_clin_vars <- unlist(
                lapply(
                    clin_var,
                    function(v) {
                        if(typeof(wilcox_test_x_expr) == 'list') {
                            x_expr <- wilcox_test_x_expr[[which(clin_var == v)]]
                        } else {
                            x_expr <- wilcox_test_x_expr
                        }
                        if(typeof(wilcox_test_y_expr) == 'list') {
                            y_expr <- wilcox_test_y_expr[[which(clin_var == v)]]
                        } else {
                            y_expr <- wilcox_test_y_expr
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
                } else {
                    ct_scores <- scrabble::score(
                        set_colnames(t(expression_data_table[sample_ids, -'id']), sample_ids),
                        list(genes_vec),
                        bin.control = TRUE,
                        nbin = (length(expression_data_table) - 1) %/% 110
                    )[, 1]
                }
            } else {
                filtered_genes <- expression_data_table[sample_ids, names(.SD)[apply(.SD, 2, scores_filter_fun)], .SDcols = -'id']
                if(score_method == 'average') {
                    ct_scores <- setNames(rowMeans(expression_data_table[sample_ids, genes_vec[genes_vec %in% filtered_genes], with = FALSE]), sample_ids)
                } else {
                    ct_scores <- scrabble::score(
                        set_colnames(t(expression_data_table[sample_ids, ..filtered_genes]), sample_ids),
                        list(genes_vec[genes_vec %in% filtered_genes]),
                        bin.control = TRUE,
                        nbin = length(filtered_genes) %/% 110
                    )[, 1]
                }
            }

            test_results <- rbindlist(
                lapply(
                    names(ct_clin_vars),
                    function(v) {

                        # Check that we have enough observations for the Wilcoxon test (by default, we only ignore cases with no observations):

                        # criterion <- data_for_test[, c(sum(eval(ct_clin_vars[[v]]$x_expr)) < min_samples, sum(eval(ct_clin_vars[[v]]$y_expr)) < min_samples)]

                        # The following actually takes a surprisingly long time...
                        criterion <- clinical_data_table[patient_ids, .(variable = get(v))][
                            !is.na(variable) & variable != '',
                            c(sum(eval(ct_clin_vars[[v]]$x_expr)) < min_samples, sum(eval(ct_clin_vars[[v]]$y_expr)) < min_samples)
                        ]

                        if(sum(criterion) > 0) {

                            cat("Not enough ", paste(c('x', 'y')[criterion], collapse = ' or '), " observations for ", v, ".\n", sep = '')

                        } else {

                            # Restrict sample and patient IDs to those where the clinical variable is not NA:
                            sample_ids <- sample_ids[clinical_data_table[patient_ids, !is.na(get(v)) & get(v) != '']]
                            patient_ids <- patient_ids[clinical_data_table[patient_ids, !is.na(get(v)) & get(v) != '']]

                            # Calculate scores by chosen method:
                            # if(is.null(scores_filter_fun)) {
                            #     if(score_method == 'average') {
                            #         data_for_test <- clinical_data_table[
                            #             patient_ids,
                            #             .(id = id, variable = get(v), score = rowMeans(expression_data_table[sample_ids, ..genes_vec]))
                            #         ]
                            #     } else {
                            #         data_for_test <- clinical_data_table[
                            #             patient_ids,
                            #             .(
                            #                 id = id,
                            #                 variable = get(v),
                            #                 score = scrabble::score(
                            #                     set_colnames(t(expression_data_table[sample_ids, -'id']), sample_ids),
                            #                     list(genes_vec),
                            #                     bin.control = TRUE,
                            #                     nbin = (length(expression_data_table) - 1) %/% 110
                            #                 )[, 1]
                            #             )
                            #         ]
                            #     }
                            # } else {
                            #     filtered_genes <- expression_data_table[sample_ids, names(.SD)[apply(.SD, 2, scores_filter_fun)], .SDcols = -'id']
                            #     if(score_method == 'average') {
                            #         data_for_test <- clinical_data_table[
                            #             patient_ids,
                            #             .(
                            #                 id = id,
                            #                 variable = get(v),
                            #                 score = rowMeans(expression_data_table[sample_ids, genes_vec[genes_vec %in% filtered_genes], with = FALSE])
                            #             )
                            #         ]
                            #     } else {
                            #         data_for_test <- clinical_data_table[
                            #             patient_ids,
                            #             .(
                            #                 id = id,
                            #                 variable = get(v),
                            #                 score = scrabble::score(
                            #                     set_colnames(t(expression_data_table[sample_ids, ..filtered_genes]), sample_ids),
                            #                     list(genes_vec[genes_vec %in% filtered_genes]),
                            #                     bin.control = TRUE,
                            #                     nbin = length(filtered_genes) %/% 110
                            #                 )[, 1]
                            #             )
                            #         ]
                            #     }
                            # }

                            # Running Wilcoxon rank-sum test and calculating fold change:

                            # return(
                            #     data_for_test[
                            #         ,
                            #         c(
                            #             test_name = ct,
                            #             variable_name = ct_clin_vars[[v]]$variable_name,
                            #             variable_match = v,
                            #             setNames(
                            #                 wilcox.test(.SD[eval(ct_clin_vars[[v]]$x_expr), score], .SD[eval(ct_clin_vars[[v]]$y_expr), score])[c('statistic', 'p.value')],
                            #                 c('stat', 'pval')
                            #             ),
                            #             fold_change = .SD[eval(ct_clin_vars[[v]]$x_expr), mean(score)]/.SD[eval(ct_clin_vars[[v]]$y_expr), mean(score)] - 1
                            #         )
                            #     ]
                            # )

                            return(
                                clinical_data_table[
                                    patient_ids,
                                    .(id = id, variable = get(v), score = ct_scores[sample_ids])
                                ][
                                    ,
                                    c(
                                        test_name = ct,
                                        variable_name = ct_clin_vars[[v]]$variable_name,
                                        variable_match = v,
                                        setNames(
                                            wilcox.test(.SD[eval(ct_clin_vars[[v]]$x_expr), score], .SD[eval(ct_clin_vars[[v]]$y_expr), score])[c('statistic', 'p.value')],
                                            c('stat', 'pval')
                                        ),
                                        fold_change = .SD[eval(ct_clin_vars[[v]]$x_expr), mean(score)]/.SD[eval(ct_clin_vars[[v]]$y_expr), mean(score)] - 1
                                    )
                                ]
                            )

                        }

                    }
                )
            )

            cat('Done!\n')

            test_results

        },

        simplify = FALSE,
        USE.NAMES = TRUE

    )

    # Bind together into a single data table:

    clinical_analyses <- rbindlist(clinical_analyses)

    if(nrow(clinical_analyses) > 0) {
        clinical_analyses[, c('test_id', 'pval_adj') := .(.I, p.adjust(pval, p_adjust_method))]
        setcolorder(clinical_analyses, c('test_id', 'test_name', 'variable_name', 'variable_match', 'stat', 'pval', 'pval_adj', 'fold_change'))
        return(clinical_analyses)
    } else {
        return(NULL)
    }

    # clinical_analyses[, c('test_id', 'pval_adj') := .(.I, p.adjust(pval, p_adjust_method))]
    # setcolorder(clinical_analyses, c('test_id', 'test_name', 'variable_name', 'variable_match', 'stat', 'pval', 'pval_adj', 'fold_change'))

    # Output:
    cat('\n')
    clinical_analyses

}





# The following function is the same as the above one except that you specify the
# criteria for the Wilcoxon tests using functions instead of expressions.  This
# allows for more flexibility with variable names.

clinical_test_2 <- function(

    expression_data,
    sample_ids_list,
    genes,
    clinical_data,
    clin_var,
    wilcox_test_x_fun,
    wilcox_test_y_fun,
    amatch_max_dist = 10,
    p_adjust_method = c('BH', 'bonferroni'),
    min_samples = 1,
    score_method = c('average', 'scrabble'),
    scores_filter_fun = NULL, # function(x) {mean(x) >= 2} is probably a good alternative
    genes_for_data_transform = NULL

) {

    # <sample_ids_list> should have names!

    # <genes> can either be a character vector or a list of the same length as <sample_ids_list>
    # (with the same names) with each element being a character vector.  If the former, the same
    # genes are used for all the elements of <sample_ids_list>; if the latter, the genes in the
    # n-th element of <genes> are used for the n-the element of <sample_ids_list>.

    # Similarly, <expression_data> and <clinical_data> can either be single data tables or lists
    # of data tables (with the same length names as <sample_ids_list>), depending on whether you
    # want to use the same table for all cancer types or a different one for each.  If
    # <sample_ids_list> is not given and <expression_data> is a list, the sample IDs will be
    # taken from the ids in <expression_data>.

    # Each of <wilcox_test_x_expr> and <wilcox_test_y_expr> can be either a list or a single
    # expression.  If the latter, the same expression is used for every clinical variable; if the
    # former, the n-th expression will be used for the n-th variable.

    # By default, matches for the elements of <clin_var> are searched for using amatch() from the
    # stringdist package, with maxDist equal to 10.  The <amatch_max_dist> argument is passed to
    # the maxDist parameter of amatch(), so can be used to specify how closely you need the
    # variable names in the data to match those in <clin_var>.  If exact matches are required,
    # set <amatch_max_dist> to NULL.

    # <score_method> is the choice of method for assigning scores to each sample based on <genes>.
    # There are two options: 'average' simply calculates the average expression of the genes in
    # <genes> for each sample; 'scrabble' uses the score() function from the scrabble package
    # with bins proportional to the number of genes in the subset of <expression_data> defined
    # by <sample_ids> and any filtering of <genes>, which can be specified under
    # <scores_filter_fun>.

    if(missing(sample_ids_list) & !is.data.frame(expression_data)) {
        sample_ids_list <- lapply(expression_data, `[[`, 'id')
    }

    p_adjust_method <- match.arg(p_adjust_method)

    score_method <- match.arg(score_method)

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

            # The following is for the (optional) data transform by calculating
            # residuals.  I do this before refining sample_ids because I want
            # the transformation to be the same as in deconv_data.

            if(!is.null(genes_for_data_transform)) {

                if(typeof(genes_for_data_transform) == 'list') {
                    genes_transform <- genes_for_data_transform[[ct]]
                } else {
                    genes_transform <- genes_for_data_transform
                }

                expression_data_table <- expression_data_table[
                    sample_ids
                ][
                    ,
                    row_sums := rowSums(.SD),
                    .SDcols = genes_transform
                ]

                expression_data_table <- as.data.table(
                    sapply(
                        genes_vec,
                        function(g) {

                            lm(
                                formula(
                                    paste0(
                                        '`',
                                        g,
                                        '` ~ row_sums'
                                    )
                                ),
                                data = expression_data_table
                            )$residuals

                        },
                        simplify = FALSE,
                        USE.NAMES = TRUE
                    )
                )[
                    ,
                    id := sample_ids
                ]

                setcolorder(expression_data_table, 'id')
                setkey(expression_data_table, id)

            }

            # Match IDs from clinical data and expression data on the level of patients,
            # i.e. first extract the first 3 components of each element in sample_ids
            # and then look for matches for these patient IDs in clinical_data_table$id.
            # We exclude non-unique matches.

            all_ids <- data.table(
                expression_data_id = sample_ids,
                patient_id = apply(
                    stringr::str_split_fixed(
                        sample_ids,
                        '\\.',
                        4
                    )[, 1:3],
                    1,
                    paste,
                    collapse = '.'
                )
            )[
                ,
                clinical_data_id := sapply(
                    patient_id,
                    function(x) {

                        matches <- clinical_data_table[
                            ,
                            id[grep(paste0('^', x), id)]
                        ]

                        if(length(matches) == 1) {
                            return(matches)
                        } else {
                            return(NA)
                        }

                    }
                )
            ][
                !(
                    clinical_data_id %in%
                        names(table(clinical_data_id))[
                            table(clinical_data_id) > 1
                        ]
                )
            ]

            # The following two lines just allow consistency with the rest of the code for this
            # function, which I wrote before the addition of all_ids.

            patient_ids <- all_ids$clinical_data_id
            sample_ids <- all_ids$expression_data_id

            # The following will be useful later for speeding things up:

            clinical_data_table_filtered <- clinical_data_table[patient_ids]

            # Look for a match for clin_var, with specified max. distance, among the variables
            # that are not all NA or the empty string for this cancer type/subtype, and return
            # NA if we don't find one:

            non_na_names <- colnames(
                clinical_data_table_filtered[
                    ,
                    sapply(
                        .SD,
                        function(x) {
                            switch((sum(!is.na(x) & x != '') == 0) + 1, x, NULL)
                        }
                    )
                ]
            )

            # The following returns a list with names being the non-NA variable name matches and
            # each element being a list consisting of the corresponding x and y Wilcoxon test
            # expressions along with the input variable name.

            ct_clin_vars <- unlist(
                lapply(
                    clin_var,
                    function(v) {

                        if(typeof(wilcox_test_x_fun) == 'list') {
                            x_fun <- wilcox_test_x_fun[[which(clin_var == v)]]
                        } else {
                            x_fun <- wilcox_test_x_fun
                        }

                        if(typeof(wilcox_test_y_fun) == 'list') {
                            y_fun <- wilcox_test_y_fun[[which(clin_var == v)]]
                        } else {
                            y_fun <- wilcox_test_y_fun
                        }

                        if(is.null(amatch_max_dist)) {

                            if(v %in% non_na_names) {
                                return(
                                    setNames(
                                        list( # List whose one element is a list with 3 elements
                                            list(
                                                variable_name = v,
                                                x_fun = x_fun,
                                                y_fun = y_fun
                                            )
                                        ),
                                        v # Name for the single element of the 'outer' list
                                    )
                                )
                            } else {
                                cat("Cannot find", v, "in data.\n")
                            }

                        } else {

                            var_match <- non_na_names[
                                stringdist::amatch(
                                    v,
                                    non_na_names,
                                    maxDist = amatch_max_dist
                                )
                            ]

                            if(is.na(var_match)) {
                                cat("No variable match for ", v, ".\n", sep = '')
                            } else {
                                cat("Found match for ", v, ": ", var_match, ".\n", sep = '')
                                return(
                                    setNames(
                                        list( # List whose one element is a list with 3 elements
                                            list(
                                                variable_name = v,
                                                x_fun = x_fun,
                                                y_fun = y_fun
                                            )
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

                    ct_scores <- setNames(
                        rowMeans(
                            expression_data_table[
                                sample_ids,
                                ..genes_vec
                            ]
                        ),
                        sample_ids
                    )

                } else {

                    ct_scores <- scrabble::score(
                        set_colnames(
                            t(expression_data_table[sample_ids, -'id']),
                            sample_ids
                        ),
                        list(genes_vec),
                        bin.control = TRUE,
                        nbin = (length(expression_data_table) - 1) %/% 110
                    )[, 1]

                }

            } else {

                filtered_genes <- expression_data_table[
                    sample_ids,
                    names(.SD)[apply(.SD, 2, scores_filter_fun)],
                    .SDcols = -'id'
                ]

                if(score_method == 'average') {

                    ct_scores <- setNames(
                        rowMeans(
                            expression_data_table[
                                sample_ids,
                                genes_vec[genes_vec %in% filtered_genes],
                                with = FALSE
                            ]
                        ),
                        sample_ids
                    )

                } else {

                    ct_scores <- scrabble::score(
                        set_colnames(
                            t(
                                expression_data_table[
                                    sample_ids,
                                    ..filtered_genes
                                ]
                            ),
                            sample_ids
                        ),
                        list(genes_vec[genes_vec %in% filtered_genes]),
                        bin.control = TRUE,
                        nbin = length(filtered_genes) %/% 110
                    )[, 1]

                }

            }

            cat(length(ct_clin_vars), 'variables to test.\n')

            test_results <- rbindlist(
                lapply(
                    names(ct_clin_vars),
                    function(v) {

                        # Check that we have enough observations for the Wilcoxon test
                        # (by default, we only ignore cases with no observations):

                        criterion <- clinical_data_table_filtered[
                            !is.na(get(v)) & get(v) != '',
                            c(
                                sum(ct_clin_vars[[v]]$x_fun(get(v))) < min_samples,
                                sum(ct_clin_vars[[v]]$y_fun(get(v))) < min_samples
                            )
                        ]

                        i <- which(names(ct_clin_vars) == v)

                        if(sum(criterion) > 0) {

                            cat(
                                rep(
                                    '\b',
                                    switch(
                                        (i == 1) + 1,
                                        floor(log10(i - 1)) + 1,
                                        0
                                    )
                                ),
                                "Not enough ",
                                paste(c('x', 'y')[criterion], collapse = ' or '),
                                " observations for ",
                                v,
                                ".\n",
                                i,
                                sep = ''
                            )

                        } else {

                            cat(
                                rep(
                                    '\b',
                                    switch(
                                        (i == 1) + 1,
                                        floor(log10(i - 1)) + 1,
                                        0
                                    )
                                ),
                                i,
                                sep = ''
                            )

                            # Restrict sample and patient IDs to those where the clinical variable
                            # is not NA:

                            sample_ids <- sample_ids[
                                clinical_data_table_filtered[
                                    ,
                                    !is.na(get(v)) & get(v) != ''
                                ]
                            ]

                            patient_ids <- patient_ids[
                                clinical_data_table_filtered[
                                    ,
                                    !is.na(get(v)) & get(v) != ''
                                ]
                            ]

                            # Running Wilcoxon rank-sum test and calculating fold change:

                            return(
                                clinical_data_table_filtered[
                                    patient_ids,
                                    c(
                                        test_name = ct,
                                        variable_name = ct_clin_vars[[v]]$variable_name,
                                        variable_match = v,
                                        setNames(
                                            wilcox.test(
                                                ct_scores[sample_ids[ct_clin_vars[[v]]$x_fun(get(v))]],
                                                ct_scores[sample_ids[ct_clin_vars[[v]]$y_fun(get(v))]]
                                            )[c('statistic', 'p.value')],
                                            c('stat', 'pval')
                                        ),
                                        fold_change = mean(
                                            ct_scores[sample_ids[ct_clin_vars[[v]]$x_fun(get(v))]]
                                        )/mean(
                                            ct_scores[sample_ids[ct_clin_vars[[v]]$y_fun(get(v))]]
                                        ) - 1
                                    )
                                ]
                            )

                        }

                    }
                )
            )

            cat('\nDone!\n')

            test_results

        },

        simplify = FALSE,
        USE.NAMES = TRUE

    )

    # Bind together into a single data table:

    clinical_analyses <- rbindlist(clinical_analyses)

    clinical_analyses[
        ,
        c('test_id', 'pval_adj') := .(
            .I,
            p.adjust(pval, p_adjust_method)
        )
    ]

    setcolorder(
        clinical_analyses,
        c(
            'test_id',
            'test_name',
            'variable_name',
            'variable_match',
            'stat',
            'pval',
            'pval_adj',
            'fold_change'
        )
    )

    # Output:

    cat('\n')

    clinical_analyses

}





# The following is a more concise version of the above, which does not carry out any
# checks that variables are present and sample sizes ar sufficient - you have to know
# these beforehand.  It should speed things up in cases where you know your variables
# exist and have no NAs.

# Note I also added an option to specify the test function, which previously was
# limited to wilcox.test().  Using t.test() makes it run very slightly faster, but I
# don't think it's very noticeable.  I should change the names of the arguments
# <wilcox_test_x_fun> and <wilcox_test_y_fun> to suit.

clinical_test_2_concise <- function(

    expression_data,
    sample_ids_list,
    genes,
    clinical_data,
    clin_var,
    wilcox_test_x_fun,
    wilcox_test_y_fun,
    p_adjust_method = c('BH', 'bonferroni'),
    score_method = c('average', 'scrabble'),
    scores_filter_fun = NULL, # function(x) {mean(x) >= 2} is probably a good alternative
    genes_for_data_transform = NULL,
    test_fun = wilcox.test,
    calculate_fold_change = TRUE

) {

    # <sample_ids_list> should have names!

    # <genes> can either be a character vector or a list of the same length as <sample_ids_list>
    # (with the same names) with each element being a character vector.  If the former, the same
    # genes are used for all the elements of <sample_ids_list>; if the latter, the genes in the
    # n-th element of <genes> are used for the n-the element of <sample_ids_list>.

    # Similarly, <expression_data> and <clinical_data> can either be single data tables or lists
    # of data tables (with the same length names as <sample_ids_list>), depending on whether you
    # want to use the same table for all cancer types or a different one for each.  If
    # <sample_ids_list> is not given and <expression_data> is a list, the sample IDs will be
    # taken from the ids in <expression_data>.

    # Each of <wilcox_test_x_expr> and <wilcox_test_y_expr> can be either a list or a single
    # expression.  If the latter, the same expression is used for every clinical variable; if the
    # former, the n-th expression will be used for the n-th variable.

    # By default, matches for the elements of <clin_var> are searched for using amatch() from the
    # stringdist package, with maxDist equal to 10.  The <amatch_max_dist> argument is passed to
    # the maxDist parameter of amatch(), so can be used to specify how closely you need the
    # variable names in the data to match those in <clin_var>.  If exact matches are required,
    # set <amatch_max_dist> to NULL.

    # <score_method> is the choice of method for assigning scores to each sample based on <genes>.
    # There are two options: 'average' simply calculates the average expression of the genes in
    # <genes> for each sample; 'scrabble' uses the score() function from the scrabble package
    # with bins proportional to the number of genes in the subset of <expression_data> defined
    # by <sample_ids> and any filtering of <genes>, which can be specified under
    # <scores_filter_fun>.

    if(missing(sample_ids_list) & !is.data.frame(expression_data)) {
        sample_ids_list <- lapply(expression_data, `[[`, 'id')
    }

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

            # The following is for the (optional) data transform by calculating
            # residuals.  I do this before refining sample_ids because I want
            # the transformation to be the same as in deconv_data.

            if(!is.null(genes_for_data_transform)) {

                if(typeof(genes_for_data_transform) == 'list') {
                    genes_transform <- genes_for_data_transform[[ct]]
                } else {
                    genes_transform <- genes_for_data_transform
                }

                expression_data_table <- expression_data_table[
                    sample_ids
                ][
                    ,
                    row_sums := rowSums(.SD),
                    .SDcols = genes_transform
                ]

                expression_data_table <- as.data.table(
                    sapply(
                        genes_vec,
                        function(g) {

                            lm(
                                formula(
                                    paste0(
                                        '`',
                                        g,
                                        '` ~ row_sums'
                                    )
                                ),
                                data = expression_data_table
                            )$residuals

                        },
                        simplify = FALSE,
                        USE.NAMES = TRUE
                    )
                )[
                    ,
                    id := sample_ids
                ]

                setcolorder(expression_data_table, 'id')
                setkey(expression_data_table, id)

            }

            # Match IDs from clinical data and expression data on the level of patients,
            # i.e. first extract the first 3 components of each element in sample_ids
            # and then look for matches for these patient IDs in clinical_data_table$id.
            # We exclude non-unique matches.

            all_ids <- data.table(
                expression_data_id = sample_ids,
                patient_id = apply(
                    stringr::str_split_fixed(
                        sample_ids,
                        '\\.',
                        4
                    )[, 1:3],
                    1,
                    paste,
                    collapse = '.'
                )
            )[
                ,
                clinical_data_id := sapply(
                    patient_id,
                    function(x) {

                        matches <- clinical_data_table[
                            ,
                            id[grep(paste0('^', x), id)]
                        ]

                        if(length(matches) == 1) {
                            return(matches)
                        } else {
                            return(NA)
                        }

                    }
                )
            ][
                !(
                    clinical_data_id %in%
                        names(table(clinical_data_id))[
                            table(clinical_data_id) > 1
                        ]
                )
            ]

            # The following two lines just allow consistency with the rest of the code for this
            # function, which I wrote before the addition of all_ids.

            patient_ids <- all_ids$clinical_data_id
            sample_ids <- all_ids$expression_data_id

            # The following will be useful later for speeding things up:

            # clinical_data_table_filtered <- clinical_data_table[patient_ids]

            # I'm transposing it so I can key by variable name, making it much faster to search
            # for the variables in the lapply() application below.

            clinical_data_table_filtered <- tdt(clinical_data_table[patient_ids])
            setkey(clinical_data_table_filtered, id)

            # non_na_names <- colnames(clinical_data_table_filtered)

            # The following returns a list with names being the non-NA variable name matches and
            # each element being a list consisting of the corresponding x and y Wilcoxon test
            # expressions along with the input variable name.

            if(typeof(wilcox_test_x_fun) == 'list' | typeof(wilcox_test_y_fun) == 'list') {

                ct_clin_vars <- unlist(
                    lapply(
                        clin_var[clin_var %in% clinical_data_table_filtered$id],
                        function(v) {

                            if(typeof(wilcox_test_x_fun) == 'list') {
                                x_fun <- wilcox_test_x_fun[[which(clin_var == v)]]
                            } else {
                                x_fun <- wilcox_test_x_fun
                            }

                            if(typeof(wilcox_test_y_fun) == 'list') {
                                y_fun <- wilcox_test_y_fun[[which(clin_var == v)]]
                            } else {
                                y_fun <- wilcox_test_y_fun
                            }

                            setNames(
                                list( # List whose one element is a list with 3 elements
                                    list(
                                        x_fun = x_fun,
                                        y_fun = y_fun
                                    )
                                ),
                                v # Name for the single element of the 'outer' list
                            )

                        }
                    ),
                    recursive = FALSE
                )

            } else {

                x_fun <- wilcox_test_x_fun
                y_fun <- wilcox_test_y_fun
                ct_clin_vars <- clin_var[clin_var %in% clinical_data_table_filtered$id]

            }

            # Calculate scores by chosen method:

            if(is.null(scores_filter_fun)) {

                if(score_method == 'average') {

                    ct_scores <- setNames(
                        rowMeans(
                            expression_data_table[
                                sample_ids,
                                ..genes_vec
                            ]
                        ),
                        sample_ids
                    )

                } else {

                    ct_scores <- scrabble::score(
                        set_colnames(
                            t(expression_data_table[sample_ids, -'id']),
                            sample_ids
                        ),
                        list(genes_vec),
                        bin.control = TRUE,
                        nbin = (length(expression_data_table) - 1) %/% 110
                    )[, 1]

                }

            } else {

                filtered_genes <- expression_data_table[
                    sample_ids,
                    names(.SD)[apply(.SD, 2, scores_filter_fun)],
                    .SDcols = -'id'
                ]

                if(score_method == 'average') {

                    ct_scores <- setNames(
                        rowMeans(
                            expression_data_table[
                                sample_ids,
                                genes_vec[genes_vec %in% filtered_genes],
                                with = FALSE
                            ]
                        ),
                        sample_ids
                    )

                } else {

                    ct_scores <- scrabble::score(
                        set_colnames(
                            t(
                                expression_data_table[
                                    sample_ids,
                                    ..filtered_genes
                                ]
                            ),
                            sample_ids
                        ),
                        list(genes_vec[genes_vec %in% filtered_genes]),
                        bin.control = TRUE,
                        nbin = length(filtered_genes) %/% 110
                    )[, 1]

                }

            }

            cat(length(ct_clin_vars), 'variables to test.\n')

            if(typeof(wilcox_test_x_fun) == 'list' | typeof(wilcox_test_y_fun) == 'list') {

                test_results <- rbindlist(
                    lapply(
                        names(ct_clin_vars),
                        function(v) {

                            i <- which(names(ct_clin_vars) == v)

                            cat(
                                rep(
                                    '\b',
                                    switch(
                                        (i == 1) + 1,
                                        floor(log10(i - 1)) + 1,
                                        0
                                    )
                                ),
                                i,
                                sep = ''
                            )

                            # Running Wilcoxon rank-sum test and calculating fold change:

                            # out <- clinical_data_table_filtered[
                            #     ,
                            #     c(
                            #         test_name = ct,
                            #         variable_name = ct_clin_vars[[v]]$variable_name,
                            #         variable_match = v,
                            #         setNames(
                            #             test_fun(
                            #                 ct_scores[sample_ids[ct_clin_vars[[v]]$x_fun(get(v))]],
                            #                 ct_scores[sample_ids[ct_clin_vars[[v]]$y_fun(get(v))]]
                            #             )[c('statistic', 'p.value')],
                            #             c('stat', 'pval')
                            #         )
                            #     )
                            # ]

                            # if(calculate_fold_change) {
                            #     out[
                            #         ,
                            #         fold_change := mean(
                            #             ct_scores[sample_ids[ct_clin_vars[[v]]$x_fun(get(v))]]
                            #         )/mean(
                            #             ct_scores[sample_ids[ct_clin_vars[[v]]$y_fun(get(v))]]
                            #         ) - 1
                            #     ]
                            # }

                            # I thought the following would be quite elegant, but the normal square
                            # bracket can't accept a list.  I could easily write a version of `[`
                            # that can accept a list, perhaps using lapply(), but I'm not sure it's
                            # worth the hassle: I'm not sure it will speed things up noticeably.

                            # ct_scores %>% `[`(
                            #     sample_ids %>% `[`(
                            #         as.numeric(clinical_data_table_filtered[v, -'id']) %>% (
                            #             function(vec) {
                            #                 lapply(
                            #                     ct_clin_vars[[v]][c('x_fun', 'y_fun')],
                            #                     function(my_fun) my_fun(vec)
                            #                 )
                            #             }
                            #         )
                            #     )
                            # )

                            # Could use the following:

                            # sub_by_list <- function(x, l) {
                            #     lapply(l, function(li) x[li])
                            # }
                            #
                            # ct_scores %>% sub_by_list(
                            #     sample_ids %>% sub_by_list(
                            #         as.numeric(clinical_data_table_filtered[v, -'id']) %>% (
                            #             function(vec) {
                            #                 lapply(
                            #                     ct_clin_vars[[v]][c('x_fun', 'y_fun')],
                            #                     function(my_fun) my_fun(vec)
                            #                 )
                            #             }
                            #         )
                            #     )
                            # )

                            # But microbenchmark shows that there's no difference between this and
                            # the stuff below.

                            out <- c(
                                test_name = ct,
                                variable_name = v,
                                setNames(
                                    do.call(
                                        test_fun,
                                        args = setNames(
                                            as.numeric(clinical_data_table_filtered[v, -'id']) %>% (
                                                function(vec) {
                                                    list(
                                                        ct_scores[sample_ids[ct_clin_vars[[v]]$x_fun(vec)]],
                                                        ct_scores[sample_ids[ct_clin_vars[[v]]$y_fun(vec)]]
                                                    )
                                                }
                                            ),
                                            c('x', 'y')
                                        )
                                    )[c('statistic', 'p.value')],
                                    c('stat', 'pval')
                                )
                            )

                            if(calculate_fold_change) {
                                out <- c(
                                    out,
                                    fold_change = as.numeric(clinical_data_table_filtered[v, -'id']) %>% (
                                        function(vec) {
                                            mean(ct_scores[sample_ids[ct_clin_vars[[v]]$x_fun(vec)]])/
                                                mean(ct_scores[sample_ids[ct_clin_vars[[v]]$y_fun(vec)]]) - 1
                                        }
                                    )
                                )
                            }

                            out

                        }
                    )
                )

            } else {

                test_results <- rbindlist(
                    lapply(
                        ct_clin_vars,
                        function(v) {

                            i <- which(ct_clin_vars == v)

                            cat(
                                rep(
                                    '\b',
                                    switch(
                                        (i == 1) + 1,
                                        floor(log10(i - 1)) + 1,
                                        0
                                    )
                                ),
                                i,
                                sep = ''
                            )

                            out <- c(
                                test_name = ct,
                                variable_name = v,
                                setNames(
                                    do.call(
                                        test_fun,
                                        args = setNames(
                                            as.numeric(clinical_data_table_filtered[v, -'id']) %>% (
                                                function(vec) {
                                                    list(
                                                        ct_scores[sample_ids[x_fun(vec)]],
                                                        ct_scores[sample_ids[y_fun(vec)]]
                                                    )
                                                }
                                            ),
                                            c('x', 'y')
                                        )
                                    )[c('statistic', 'p.value')],
                                    c('stat', 'pval')
                                )
                            )

                            if(calculate_fold_change) {
                                out <- c(
                                    out,
                                    fold_change = as.numeric(clinical_data_table_filtered[v, -'id']) %>% (
                                        function(vec) {
                                            mean(ct_scores[sample_ids[x_fun(vec)]])/
                                                mean(ct_scores[sample_ids[y_fun(vec)]]) - 1
                                        }
                                    )
                                )
                            }

                            out

                        }
                    )
                )

            }

            cat('\nDone!\n')

            test_results

        },

        simplify = FALSE,
        USE.NAMES = TRUE

    )

    # Bind together into a single data table:

    clinical_analyses <- rbindlist(clinical_analyses)

    clinical_analyses[
        ,
        c('test_id', 'pval_adj') := .(
            .I,
            p.adjust(pval, p_adjust_method)
        )
    ]

    new_col_order <- c(
        'test_id',
        'test_name',
        'variable_name',
        'stat',
        'pval',
        'pval_adj'
    )

    if('fold_change' %in% names(clinical_analyses)) {
        new_col_order <- c(new_col_order, 'fold_change')
    }

    setcolorder(
        clinical_analyses,
        new_col_order
    )

    # Output:

    cat('\n')

    clinical_analyses

}





adj_signif_threshold <- function(

    signif_threshold,
    p_adjust_method = c('BH', 'bonferroni'),
    test_data = NULL

) {

    p_adjust_method = match.arg(p_adjust_method)

    if(p_adjust_method == 'BH') {

        if(is.null(test_data)) {
            stop("Please provide <test_data> if using method 'BH'.")
        }

        test_data[
            order(pval),
            switch(
                (sum(pval <= signif_threshold*.I/.N) == 0) + 1,
                signif_threshold*max(
                    which(
                        pval <= signif_threshold*.I/.N
                    )
                )/.N,
                NA
            )
        ]

    } else { # Then p_adjust_method == 'bonferroni'

        signif_threshold/nrow(test_data)

    }

}





clinical_test_volcano <- function(

    test_data,
    rows_for_labels = 1:nrow(test_data),
    fold_change_fun = function(fc) log2(fc + 1),
    xlab_TeX = 'log_2(fold change + 1)',
    ylab_TeX = '-log_{10}(p value)',
    legend_labels = NULL,
    legend_colours = NULL,
    legend_title = waiver(),
    plot_title = waiver(),
    show_signif_threshold = 0.05,
    signif_threshold_colour = 'darkorange',
    labels_fun = ggrepel::geom_label_repel,
    ...

) {

    # <rows_for_labels> should be an integer or logical indicating which rows to label
    # on the plot.  The points corresponding to those rows will be labelled using
    # ggrepel.

    # Any value supplied to <show_signif_threshold> will be converted to -log10 space.

    # The '...' argument is for additional arguments to <labels_fun>.

    fold_change_fun <- match.fun(fold_change_fun)
    labels_fun <- match.fun(labels_fun)

    test_labels <- test_data$test_name

    if(is.integer(rows_for_labels)) {
        test_labels[-rows_for_labels] <- ''
    } else if(is.logical(rows_for_labels)) {
        test_labels[!rows_for_labels] <- ''
    } else {
        stop('<rows_for_labels> must be integer or logical.')
    }

    var_names <- sort(unique(test_data$variable_name))

    g <- ggplot(
        test_data
    ) +
        geom_point(
            aes(
                x = fold_change_fun(fold_change),
                y = -log10(pval),
                colour = variable_name
            ),
            shape = 17,
            size = 2.5
        ) +
        labels_fun(
            aes(
                x = fold_change_fun(fold_change),
                y = -log10(pval),
                label = test_labels
            ),
            ...
        ) +
        scale_colour_manual(
            labels = switch(
                is.null(legend_labels) + 1,
                legend_labels,
                var_names
            ),
            values = switch(
                is.null(legend_colours) + 1,
                legend_colours,
                setNames(
                    hcl(
                        h = seq(
                            360/(2*length(var_names)),
                            360 - 360/(2*length(var_names)),
                            length.out = length(var_names)
                        ),
                        c = 100,
                        l = 65
                    ),
                    var_names
                )
            )
        ) +
        labs(
            x = latex2exp::TeX(xlab_TeX),
            y = latex2exp::TeX(ylab_TeX),
            colour = legend_title,
            title = plot_title
        ) +
        theme_test()

    if(!is.null(show_signif_threshold)) {
        g <- g + geom_hline(
            aes(yintercept = -log10(show_signif_threshold)),
            colour = signif_threshold_colour,
            linetype = 'dashed'
        )
    }

    g

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

    # If <hclust_method> is NULL, no clustering will be done.

    # Use <x_factor_levels> and <y_factor_levels> to manually specify the order in which
    # you want the variable values to appear on the axes.  These arguments will be
    # ignored if <hclust_method> is NULL (maybe this doesn't make sense - we could make
    # it so that it overrides the hclust ordering, as an easy way to apply the hclust
    # ordering to only one axis).

    # The '...' argument is for extra arguments for theme().

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

        hclust_x <- hclust(
            dist(t(clust_data[, -'y_var'])),
            method = hclust_method
        )

        hclust_y <- hclust(
            dist(
                clust_data[
                    ,
                    magrittr::set_rownames(as.matrix(.SD), y_var),
                    .SD = -'y_var'
                ]
            ),
            method = hclust_method
        )

        htmp <- ggplot(
            test_data[
                expand.grid(
                    unique(get(x_var)),
                    unique(get(y_var))
                )
            ],
            aes(
                x = factor(
                    get(x_var),
                    levels = hclust_x$labels[hclust_x$order]
                ),
                y = factor(
                    get(y_var),
                    levels = hclust_y$labels[hclust_y$order]
                ),
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
            theme(
                panel.border = element_rect(size = 0.5, fill = NA),
                axis.text.x = element_text(angle = 55, hjust = 1),
                ...
                # legend.title = element_text(angle = 90)
            ) +
            labs(
                x = x_lab,
                y = y_lab,
                fill = legend_title,
                title = plot_title
            )

        list(
            heatmap = htmp,
            hclust_x = hclust_x,
            hclust_y = hclust_y
        )

    } else {

        ggplot(
            test_data[
                expand.grid(
                    unique(get(x_var)),
                    unique(get(y_var))
                )
            ],
            aes(
                x = switch(
                    is.null(x_factor_levels) + 1,
                    factor(get(x_var), levels = x_factor_levels),
                    get(x_var)
                ),
                y = switch(
                    is.null(y_factor_levels) + 1,
                    factor(get(y_var), levels = y_factor_levels),
                    get(y_var)
                ),
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
            theme(
                panel.border = element_rect(size = 0.5, fill = NA),
                axis.text.x = element_text(angle = 55, hjust = 1),
                ...
            ) +
            labs(
                x = x_lab,
                y = y_lab,
                fill = legend_title,
                title = plot_title
            )

    }

}





nice_names <- function(

    names_vec,
    keep_tcga_abbreviations = FALSE

) {

    plyr::mapvalues(
        names_vec,
        c(
            'acc',
            'blca',
            'blca_luminal_infiltrated',
            'blca_luminal_papillary',
            'blca_luminal',
            'blca_luminal_luminal_infiltrated',
            'blca_basal_squamous',
            'blca_neuronal',
            'blca_basal_squamous_neuronal',
            'brca',
            'brca_luminal_a',
            'brca_luminal_b',
            'brca_basal_like',
            'brca_her2_enriched',
            'brca_normal_like',
            'brca_karaayvaz',
            'brca_luminal_a_karaayvaz',
            'brca_luminal_b_karaayvaz',
            'brca_basal_like_karaayvaz',
            'brca_her2_enriched_karaayvaz',
            'brca_normal_like_karaayvaz',
            'cesc',
            'chol',
            'coad',
            'coadread',
            'coadread_cris_a',
            'coadread_cris_b',
            'coadread_cris_c',
            'coadread_cris_d',
            'coadread_cris_e',
            'esca',
            'esca_ac',
            'esca_escc',
            'hnsc',
            'hnsc_mesenchymal_basal',
            'hnsc_classical',
            'hnsc_atypical',
            'kich',
            'kirc',
            'kirp',
            'lihc',
            'lihc_icluster_1',
            'lihc_icluster_2',
            'lihc_icluster_3',
            'luad',
            'luad_proximal_inflammatory',
            'luad_proximal_proliferative',
            'luad_terminal_respiratory_unit',
            'lusc',
            'lusc_basal',
            'lusc_classical',
            'lusc_primitive',
            'lusc_secretory',
            'meso',
            'ov',
            'ov_differentiated',
            'ov_immunoreactive',
            'ov_mesenchymal',
            'ov_proliferative',
            'paad',
            'paad_basal_moffitt',
            'paad_classical_moffitt',
            'paad_classical_collisson',
            'paad_exocrine_collisson',
            'paad_quasi_mesenchymal_collisson',
            'paad_squamous_bailey',
            'paad_immunogenic_bailey',
            'paad_progenitor_bailey',
            'paad_adex_bailey',
            'prad',
            'prad_1',
            'prad_2',
            'prad_3',
            'read',
            'skcm',
            'skcm_immune',
            'skcm_keratin',
            'skcm_mitf_low',
            'stad',
            'stad_cin',
            'stad_ebv',
            'stad_gs',
            'stad_msi',
            'thca',
            'ucec',
            'uvm'
        ),
        switch(
            keep_tcga_abbreviations + 1,
            c(
                'Adrenocortical',
                'Bladder',
                'Bladder Luminal-Infiltrated',
                'Bladder Luminal-Papillary',
                'Bladder Luminal',
                'Bladder Luminal/Luminal-Infiltrated',
                'Bladder Basal-Squamous',
                'Bladder Neuronal',
                'Bladder Basal-Squamous/Neuronal',
                'Breast',
                'Breast Luminal A',
                'Breast Luminal B',
                'Breast Basal-like',
                'Breast HER2-enriched',
                'Breast Normal-like',
                'Breast',
                'Breast Luminal A',
                'Breast Luminal B',
                'Breast Basal-like',
                'Breast HER2-enriched',
                'Breast Normal-like',
                'Cervical',
                'Cholangiocarcinoma',
                'Colon',
                'Colorectal',
                'Colorectal CRIS A',
                'Colorectal CRIS B',
                'Colorectal CRIS C',
                'Colorectal CRIS D',
                'Colorectal CRIS E',
                'Oesophageal',
                'Oesophageal Adenocarcinoma',
                'Oesophageal Squamous',
                'Head & Neck',
                'Head & Neck Malignant-Basal',
                'Head & Neck Classical',
                'Head & Neck Atypical',
                'Kidney Chromophobe',
                'Kidney Clear Cell',
                'Kidney Papillary Cell',
                'Liver',
                'Liver iCluster 1',
                'Liver iCluster 2',
                'Liver iCluster 3',
                'Lung Adenocarcinoma',
                'Lung Adenocarcinoma Squamoid',
                'Lung Adenocarcinoma Magnoid',
                'Lung Adenocarcinoma Bronchioid',
                'Lung Squamous',
                'Lung Squamous Basal',
                'Lung Squamous Classical',
                'Lung Squamous Primitive',
                'Lung Squamous Secretory',
                'Mesothelioma',
                'Ovarian',
                'Ovarian Differentiated',
                'Ovarian Immunoreactive',
                'Ovarian Mesenchymal',
                'Ovarian Proliferative',
                'Pancreatic',
                'Pancreatic Basal (Moffitt)',
                'Pancreatic Classical (Moffitt)',
                'Pancreatic Classical (Collisson)',
                'Pancreatic Exocrine (Collisson)',
                'Pancreatic Quasi-mesenchymal (Collisson)',
                'Pancreatic Squamous (Bailey)',
                'Pancreatic Immunogenic (Bailey)',
                'Pancreatic Progenitor (Bailey)',
                'Pancreatic ADEX (Bailey)',
                'Prostate',
                'Prostate 1',
                'Prostate 2',
                'Prostate 3',
                'Rectum',
                'Melanoma',
                'Melanoma Immune',
                'Melanoma Keratin',
                'Melanoma MITF-low',
                'Stomach',
                'Stomach CIN',
                'Stomach EBV',
                'Stomach GS',
                'Stomach MSI',
                'Thyroid',
                'Endometrial',
                'Uveal Melanoma'
            ),
            c(
                'ACC',
                'BLCA',
                'BLCA Luminal-Infiltrated',
                'BLCA Luminal-Papillary',
                'BLCA Luminal',
                'BLCA Luminal/Luminal-Infiltrated',
                'BLCA Basal-Squamous',
                'BLCA Neuronal',
                'BLCA Basal-Squamous/Neuronal',
                'BRCA',
                'BRCA Luminal A',
                'BRCA Luminal B',
                'BRCA Basal-like',
                'BRCA HER2-enriched',
                'BRCA Normal-like',
                'BRCA',
                'BRCA Luminal A',
                'BRCA Luminal B',
                'BRCA Basal-like',
                'BRCA HER2-enriched',
                'BRCA Normal-like',
                'CESC',
                'CHOL',
                'COAD',
                'COADREAD',
                'COADREAD CRIS A',
                'COADREAD CRIS B',
                'COADREAD CRIS C',
                'COADREAD CRIS D',
                'COADREAD CRIS E',
                'ESCA',
                'ESCA Adenocarcinoma',
                'ESCA Squamous',
                'HNSC',
                'HNSC Malignant-Basal',
                'HNSC Classical',
                'HNSC Atypical',
                'KICH',
                'KIRC',
                'KIRP',
                'LIHC',
                'LIHC iCluster 1',
                'LIHC iCluster 2',
                'LIHC iCluster 3',
                'LUAD',
                'LUAD Squamoid',
                'LUAD Magnoid',
                'LUAD Bronchioid',
                'LUSC',
                'LUSC Basal',
                'LUSC Classical',
                'LUSC Primitive',
                'LUSC Secretory',
                'MESO',
                'OV',
                'OV Differentiated',
                'OV Immunoreactive',
                'OV Mesenchymal',
                'OV Proliferative',
                'PAAD',
                'PAAD Basal (Moffitt)',
                'PAAD Classical (Moffitt)',
                'PAAD Classical (Collisson)',
                'PAAD Exocrine (Collisson)',
                'PAAD Quasi-mesenchymal (Collisson)',
                'PAAD Squamous (Bailey)',
                'PAAD Immunogenic (Bailey)',
                'PAAD Progenitor (Bailey)',
                'PAAD ADEX (Bailey)',
                'PRAD',
                'PRAD 1',
                'PRAD 2',
                'PRAD 3',
                'READ',
                'SKCM',
                'SKCM Immune',
                'SKCM Keratin',
                'SKCM MITF-low',
                'STAD',
                'STAD CIN',
                'STAD EBV',
                'STAD GS',
                'STAD MSI',
                'THCA',
                'UCEC',
                'UVM'
            )
        ),
        warn_missing = FALSE
    )

}





# It would make more sense to have a separate function that calculates EMT and CAF scores,
# or even make this part of the deconvolve_emt_caf() function.

score_cor_with_cell_types <- function(

    ids_genes_ordering_list,
    which_score = c('gene', 'tumour'),
    scores_data = NULL,
    expression_data = NULL,
    cell_type_markers = NULL,
    cell_types = c(
        'B_plasma',
        'myocyte',
        'macrophage',
        'endothelial',
        'DC',
        'mast',
        'T',
        'B'
    ),
    gene_var = 'gene',
    analysis_var = 'analysis',
    score_var = 'score',
    transform_data = TRUE,
    subtract_caf_score = FALSE,
    epi_markers_list = NULL,
    x_lab = 'Cell type',
    y_lab = 'Cancer type',
    title = waiver(),
    colours = rev(colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))(50)),
    colour_limits = c(-1, 1),
    x_axis_text_angle = 55,
    x_axis_hjust = 1

) {

    which_score <- match.arg(which_score)

    if(which_score == 'gene' & is.null(scores_data)) {

        stop('Please supply scores_data if using gene scores.')

    } else if(which_score == 'gene') {

        setkeyv(scores_data, gene_var)

        data_for_cor <- rbindlist(
            lapply(
                names(ids_genes_ordering_list),
                function(li_name) {

                    genes_cell_types_correlations <- ids_genes_ordering_list[[
                        li_name
                        ]]$diagnostics$genes_cell_types_correlations

                    setkey(
                        genes_cell_types_correlations,
                        id
                    )

                    genes <- ids_genes_ordering_list[[li_name]]$genes_filtered

                    data_for_cor <- cbind(
                        scores_data[
                            get(analysis_var) == li_name
                            ][
                                genes,
                                .(
                                    analysis = li_name,
                                    gene = get(gene_var),
                                    score = get(score_var)
                                )
                                ],
                        genes_cell_types_correlations[
                            genes,
                            ..cell_types
                            ]
                    )

                }

            )
        )

    } else if(

        which_score == 'tumour' & (is.null(cell_type_markers) | is.null(expression_data))

    ) {

        stop(
            'Please supply both cell_type_markers and expression_data if using ',
            'tumour scores.'
        )

    } else {

        tumour_emt_scores <- sapply(
            names(ids_genes_ordering_list),
            function(li_name) {

                li <- ids_genes_ordering_list[[li_name]]

                # Genes to be included in the data for use in calculating tumour scores:

                emt_markers <- head(li$genes_filtered[li$ordering], 20)

                if(subtract_caf_score) {
                    caf_markers <- tail(li$genes_filtered[li$ordering], 20)
                    genes_for_tumour_scores <- c(emt_markers, caf_markers)
                } else {
                    genes_for_tumour_scores <- emt_markers
                }

                # Define the data for use in calculating tumour scores, depending on whether
                # we want the row sums/residuals transformation or not:

                if(transform_data) {

                    data_for_tumour_scores <- as.data.table(
                        sapply(
                            genes_for_tumour_scores,
                            function(g) {

                                lm(
                                    formula(
                                        paste0(
                                            '`',
                                            g,
                                            '` ~ row_sums'
                                        )
                                    ),
                                    data = cbind(
                                        expression_data[
                                            li$sample_ids,
                                            ..genes_for_tumour_scores
                                            ],
                                        row_sums = rowSums(
                                            expression_data[
                                                li$sample_ids,
                                                li$genes_filtered,
                                                with = FALSE
                                                ]
                                        )
                                    )
                                )$residuals

                            },
                            USE.NAMES = TRUE
                        )
                    )

                } else {

                    data_for_tumour_scores <- expression_data[
                        li$sample_ids,
                        ..genes_for_tumour_scores
                        ]

                }

                # Calculate the scores from data_for_tumour_scores:

                if(subtract_caf_score) {

                    scores <- setNames(
                        data_for_tumour_scores[
                            ,
                            rowMeans(.SD[, ..emt_markers]) - rowMeans(.SD[, ..caf_markers])
                            ],
                        li$sample_ids
                    )

                } else {

                    scores <- setNames(
                        rowMeans(data_for_tumour_scores),
                        li$sample_ids
                    )

                }

                # Optionally adjust EMT scores for epi scores (maybe try loess instead of
                # lm, except we'd have to manually look at some plots to choose the right
                # span etc.):

                if(!is.null(epi_markers_list)) {

                    scores <- setNames(
                        lm(
                            score ~ epi_score,
                            data = cbind(
                                data.table(
                                    id = names(scores),
                                    score = scores
                                ),
                                epi_score = expression_data[
                                    li$sample_ids,
                                    setNames(rowMeans(.SD), id),
                                    .SDcols = epi_markers_list[[li_name]]
                                    ]
                            )
                        )$residuals,
                        li$sample_ids
                    )

                }

                scores

            }
        )

        data_for_cor <- rbindlist(
            lapply(
                names(ids_genes_ordering_list),
                function(li_name) {

                    li <- ids_genes_ordering_list[[li_name]]

                    expression_data[
                        li$sample_ids,
                        c(
                            .(
                                analysis = li_name,
                                tumour = id,
                                score = tumour_emt_scores[[li_name]]
                            ),
                            sapply(
                                cell_types,
                                function(ct) {
                                    rowMeans(
                                        .SD[
                                            ,
                                            cell_type_markers[
                                                cell_type == ct & gene %in% names(expression_data),
                                                gene
                                                ],
                                            with = FALSE
                                            ]
                                    )
                                },
                                simplify = FALSE,
                                USE.NAMES = TRUE
                            )
                        )
                        ]#[
                    #     ,
                    #     cor(emt_score, .SD),
                    #     .SDcols = cell_types
                    # ][1, ]

                }
            )
        )

    }





    # Calculate the correlations for the heatmap:

    score_cor <- data_for_cor[
        ,
        .(
            cell_type = cell_types,
            correlation = cor(score, .SD)[1, ]
        ),
        .SDcols = cell_types,
        by = analysis
        ]

    # Put in alphabetical order (because casting data does this anyway, which would mess up
    # the orderings if it's not already in alphabetical order):

    setkey(score_cor, analysis, cell_type)

    # Order the data using hierarchical clustering (need to cast data for this):

    cast_score_cor <- dcast(
        score_cor,
        analysis ~ cell_type,
        value.var = 'correlation'
    )

    ordering_cell_types <- hclust(
        dist(t(cast_score_cor[, -'analysis'])),
        method = 'average'
    )$order

    ordering_analyses <- hclust(
        dist(cast_score_cor[, -'analysis']),
        method = 'average'
    )$order

    # Make the plot:

    g <- ggplot(
        score_cor,
        aes(
            x = factor(cell_type, levels = unique(cell_type)[ordering_cell_types]),
            y = factor(analysis, levels = unique(analysis)[ordering_analyses])
        )
    ) +
        geom_raster(
            aes(fill = correlation)
        ) +
        scale_x_discrete(
            expand = c(0, 0)
        ) +
        scale_y_discrete(
            expand = c(0, 0)
        ) +
        scale_fill_gradientn(
            colours = colours,
            limits = colour_limits
        ) +
        theme(
            axis.ticks = element_blank(),
            axis.text.x = element_text(
                angle = x_axis_text_angle,
                hjust = x_axis_hjust
            ),
            axis.ticks.length = unit(0, 'pt')
        ) +
        labs(
            x = x_lab,
            y = y_lab,
            title = title
        )

    list(
        heatmap = g,
        scores = data_for_cor,
        heatmap_data = score_cor,
        analyses = unique(score_cor$analysis),
        cell_types = unique(score_cor$cell_type),
        ordering_analyses = ordering_analyses,
        ordering_cell_types = ordering_cell_types
    )

}





calculate_tumour_scores <- function(

    ids_genes_ordering_list,
    expression_data,
    transform_data = TRUE,
    epi_markers_list = NULL

) {

    rbindlist(
        lapply(
            names(ids_genes_ordering_list),
            function(li_name) {

                li <- ids_genes_ordering_list[[li_name]]

                # Genes to be included in the data for use in calculating tumour scores:

                emt_markers <- head(li$genes_filtered[li$ordering], 20)
                caf_markers <- tail(li$genes_filtered[li$ordering], 20)

                # It doesn't make sense to calculate the scores by subtracting CAF correlation from
                # EMT correlation, and vice versa, because then the EMT and CAF scores are just the
                # negatives of each other.

                if(transform_data) {

                    scores <- as.data.table(
                        sapply(
                            c(emt_markers, caf_markers),
                            function(g) {

                                lm(
                                    formula(
                                        paste0(
                                            '`',
                                            g,
                                            '` ~ row_sums'
                                        )
                                    ),
                                    data = cbind(
                                        expression_data[
                                            li$sample_ids,
                                            c(..emt_markers, ..caf_markers)
                                            ],
                                        row_sums = rowSums(
                                            expression_data[
                                                li$sample_ids,
                                                li$genes_filtered,
                                                with = FALSE
                                                ]
                                        )
                                    )
                                )$residuals

                            },
                            USE.NAMES = TRUE
                        )
                    )[
                        ,
                        .(
                            analysis = li_name,
                            id = li$sample_ids,
                            emt_score = rowMeans(.SD[, ..emt_markers]),
                            caf_score = rowMeans(.SD[, ..caf_markers])
                        )
                    ]

                } else {

                    scores <- expression_data[
                        li$sample_ids,
                        c(..emt_markers, ..caf_markers)
                    ][
                        ,
                        .(
                            analysis = li_name,
                            id = li$sample_ids,
                            emt_score = rowMeans(.SD[, ..emt_markers]),
                            caf_score = rowMeans(.SD[, ..caf_markers])
                        )
                    ]

                }

                # Optionally adjust EMT scores for epi scores (maybe try loess instead of
                # lm, except we'd have to manually look at some plots to choose the right
                # span etc.):

                if(!is.null(epi_markers_list)) {

                    if(transform_data) {

                        scores[
                            ,
                            epi_score := expression_data[
                                li$sample_ids,
                                rowMeans(.SD),
                                .SDcols = epi_markers_list[[li_name]]
                            ]
                        ]

                        scores[
                            ,
                            emt_score_adjusted := lm(
                                emt_score ~ epi_score,
                                data = .SD
                            )$residuals
                        ]

                    } else {

                        scores[
                            ,
                            epi_score := expression_data[
                                li$sample_ids,
                                rowMeans(.SD),
                                .SDcols = epi_markers_list[[li_name]]
                            ]
                        ]

                        scores[
                            ,
                            emt_score_adjusted := emt_score - epi_score
                        ]

                    }

                }

                scores

            }

        )
    )

}





deconv_plot <- function(

    deconvs_list,
    n_row = 1,
    n_col = length(deconvs_list),
    plots_rel_heights = c(
        title = 2,
        purity_bar = 1,
        ccle_bar = 1,
        extra_bar = 1,
        heatmap = 15,
        axis_labels = 5
    ),
    # plots_rel_heights = switch(
    #     ('axis_labels' %in% names(deconvs_list[[1]]$plots)) + 1,
    #     c(1.5, rep(1, length(deconvs_list[[1]]$plots) - 1), 15),
    #     c(1.5, rep(1, length(deconvs_list[[1]]$plots) - 2), 15, 5)
    # ),
    rows_rel_heights = rep(1, n_row),
    left_plot_width = 1,
    align = c('centre', 'center', 'left'),
    legends = TRUE,
    legends_arrange = switch(
        (length(deconvs_list) == n_row*n_col) + 1,
        'horizontal',
        'vertical'
    ),
    legends_rel_size = switch((legends_arrange == 'vertical') + 1, 1, 0.2),
    legends_space = 1,
    ...

) {

    # The '...' argument is for guide_colourbar(), and can be used to change the format of the
    # legends.

    # <deconvs_list> should be a list with each element of the form returned by deconvolve_emt_caf().

    # <left_plot_width> specifies the desired width of the leftmost plot in each row as a fraction
    # of the normal width of the other plots.  It is 1 by default, but it can be increased manually
    # if, for example, y axis text on the leftmost plots makes these plots slightly wider than the
    # others.

    # The <legends> argument can be used to specify whether you want the legends or not.  It is
    # TRUE by default.

    # The legends are arranged in a list as follows: a blank plot first and last, and between these
    # are the legends alternating with blank plots.  So every other element of the list is a blank
    # plot, and both the first and last elements of the list are blank plots.  This is so that you
    # can adjust the spacing between, before and after the legends by adjusting the relative sizes
    # of the blank plots.  To do this manually, it helps to know how long the legend list will be.
    # There will always be one more blank plot than legends, so the total length of the legend list
    # will be 2*(number of legends) + 1.

    # The <align> argument refers to the alignment of the plots (and legends) on the bottom row,
    # if indeed there is space to spare.  If <align> is 'centre' (or 'center'), blank plots are
    # inserted at the left and right ends to provide the padding to make the alignment central.
    # The plots will be the same width as those on the upper row(s).  The user can specify the
    # width to be taken up by the legends via the <legends_space> argument.  The number passed to
    # this is interpreted as a fraction of the width of one plot.  It is 1 by default, meaning the
    # legends take up the same width as one plot.  The widths of the blank padding plots are chosen
    # automatically to take up all remaining space.  If <align> is 'left', no blank plots are
    # inserted, <legend_space> is ignored, and the legend space takes up all the space not occupied
    # by the plots.  In this case, the <legends_rel_size> argument can be used to change the space
    # occupied by the legends.

    if(length(rows_rel_heights) != n_row) {
        stop('<rows_rel_heights> does not have length equal to <n_row>.')
    }

    align <- match.arg(align)
    legends_arrange <- match.arg(legends_arrange, choices = c('horizontal', 'vertical'))

    plot_lists <- lapply(
        1:length(deconvs_list),
        function(i) {

            if(i %in% seq(from = 1, to = length(deconvs_list), by = n_col)) {

                bars <- lapply(
                    deconvs_list[[i]]$plots[ # This subsetting returns the bars in the order: purity-ccle-extra.
                        c('purity_bar', 'ccle_bar', 'extra_bar')[
                            c('purity_bar', 'ccle_bar', 'extra_bar') %in% names(deconvs_list[[i]]$plots)
                        ]
                    ],
                    function(x) {
                        x + theme(legend.position = 'none')
                    }
                )

                pl <- c(
                    list(
                        title = ggplot() +
                            labs(caption = deconvs_list[[i]]$plots$heatmap$labels$title) +
                            theme(
                                panel.background = element_rect(fill = 'white'),
                                plot.caption = element_text(size = 14, hjust = 0)
                            )
                        # blank_dummy_plot() +
                        #     labs(title = deconvs_list[[i]]$plots$heatmap$labels$title) +
                        #     theme(plot.margin = unit(c(5.5, 5.5, 0, 5.5), 'pt'))
                    ),
                    bars,
                    list(
                        heatmap = deconvs_list[[i]]$plots$heatmap + theme(
                            legend.position = 'none',
                            plot.title = element_blank(), # Remove just plot title
                            plot.margin = unit(
                                c(
                                    0,
                                    5.5,
                                    switch(
                                        ('axis_labels' %in% names(deconvs_list[[i]]$plots)) + 1,
                                        5.5,
                                        0
                                    ),
                                    5.5
                                ),
                                'pt'
                            )
                        )
                    )
                )

                if('axis_labels' %in% names(deconvs_list[[i]]$plots)) {
                    pl <- c(pl, list(axis_labels = deconvs_list[[i]]$plots$axis_labels))
                }

            } else {

                bars <- lapply(
                    deconvs_list[[i]]$plots[
                        c('purity_bar', 'ccle_bar', 'extra_bar')[
                            c('purity_bar', 'ccle_bar', 'extra_bar') %in% names(deconvs_list[[i]]$plots)
                        ]
                    ],
                    function(x) {
                        x + theme(legend.position = 'none', axis.title.y = element_blank())
                    }
                )

                pl <- c(
                    list(
                        title = ggplot() +
                            labs(caption = deconvs_list[[i]]$plots$heatmap$labels$title) +
                            theme(
                                panel.background = element_rect(fill = 'white'),
                                plot.caption = element_text(size = 14, hjust = 0)
                            )
                        # blank_dummy_plot() +
                        #     labs(title = deconvs_list[[i]]$plots$heatmap$labels$title) +
                        #     theme(plot.margin = unit(c(5.5, 5.5, 0, 5.5), 'pt'))
                    ),
                    bars,
                    list(
                        heatmap = deconvs_list[[i]]$plots$heatmap + theme(
                            legend.position = 'none',
                            title = element_blank(), # Remove plot and y axis titles
                            plot.margin = unit(
                                c(
                                    0,
                                    5.5,
                                    switch(
                                        ('axis_labels' %in% names(deconvs_list[[i]]$plots)) + 1,
                                        5.5,
                                        0
                                    ),
                                    5.5
                                ),
                                'pt'
                            )
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
                c('purity_bar', 'ccle_bar', 'extra_bar')[
                    c('purity_bar', 'ccle_bar', 'extra_bar') %in% names(deconvs_list[[1]]$plots)
                ]
            ],
            function(g) {
                get_legend(g + guides(fill = guide_colourbar(...)))
            }
            # get_legend
        )

        legend_list <- unlist(
            lapply(
                legend_list,
                list,
                blank_plot()
            ),
            recursive = FALSE
        )

        legend_list <- c(
            list(blank_plot()),
            legend_list,
            list(
                get_legend(deconvs_list[[1]]$plots$heatmap + guides(fill = guide_colourbar(...))),
                blank_plot()
            )
        )

        if(length(legends_rel_size) == 1) {
            legends_rel_size <- rep(legends_rel_size, length(legend_list))
        }

        # Make the legend plot:

        if(legends_arrange == 'vertical') {

            legend_plot <- plot_grid(
                plotlist = legend_list,
                nrow = length(legend_list),
                ncol = 1,
                rel_heights = legends_rel_size
            )

        } else {

            legend_plot <- plot_grid(
                plotlist = legend_list,
                nrow = 1,
                ncol = length(legend_list),
                rel_widths = legends_rel_size
            )

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
                        function(l) {
                            plot_grid(
                                plotlist = l,
                                # rel_heights = plots_rel_heights,
                                rel_heights = plots_rel_heights[names(l)],
                                ncol = 1,
                                nrow = length(l),
                                align = 'v'
                            )
                        }
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
                # rel_widths = c(n_col, legends_rel_size)
            )

        } else {

            plot_grid(
                plotlist = lapply(
                    plot_lists,
                    function(l) {
                        plot_grid(
                            plotlist = l,
                            # rel_heights = plots_rel_heights,
                            rel_heights = plots_rel_heights[names(l)],
                            ncol = 1,
                            nrow = length(l),
                            align = 'v'
                        )
                    }
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
                        function(l) {
                            plot_grid(
                                plotlist = l,
                                # rel_heights = plots_rel_heights,
                                rel_heights = plots_rel_heights[names(l)],
                                ncol = 1,
                                nrow = length(l),
                                align = 'v'
                            )
                        }
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
                        function(l) {
                            plot_grid(
                                plotlist = l,
                                # rel_heights = plots_rel_heights,
                                rel_heights = plots_rel_heights[names(l)],
                                ncol = 1,
                                nrow = length(l),
                                align = 'v'
                            )
                        }
                    ),
                    nrow = length(plot_lists) %/% n_col, # This should be n_row - 1, right?
                    ncol = n_col,
                    rel_widths = c(left_plot_width, rep(1, n_col - 1)),
                    rel_heights = rows_rel_heights[-n_row],
                    align = 'hv'
                ),
                plot_grid( # Bottom row
                    plotlist = c(
                        lapply(
                            plot_lists[(length(plot_lists) - (length(plot_lists) %% n_col) + 1):length(plot_lists)],
                            function(l) {
                                plot_grid(
                                    plotlist = l,
                                    # rel_heights = plots_rel_heights,
                                    rel_heights = plots_rel_heights[names(l)],
                                    ncol = 1,
                                    nrow = length(l),
                                    align = 'v'
                                )
                            }
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
                        function(l) {
                            plot_grid(
                                plotlist = l,
                                # rel_heights = plots_rel_heights,
                                rel_heights = plots_rel_heights[names(l)],
                                ncol = 1,
                                nrow = length(l),
                                align = 'v'
                            )
                        }
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
                            function(l) {
                                plot_grid(
                                    plotlist = l,
                                    # rel_heights = plots_rel_heights,
                                    rel_heights = plots_rel_heights[names(l)],
                                    ncol = 1,
                                    nrow = length(l),
                                    align = 'v'
                                )
                            }
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

            # Is this case really any different from the one where
            # length(deconvs_list) == n_row*n_col - 1?

            plot_grid(
                plot_grid(
                    plotlist = lapply(
                        plot_lists[1:((length(plot_lists) %/% n_col)*n_col)],
                        function(l) {
                            plot_grid(
                                plotlist = l,
                                # rel_heights = plots_rel_heights,
                                rel_heights = plots_rel_heights[names(l)],
                                ncol = 1,
                                nrow = length(l),
                                align = 'v'
                            )
                        }
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
                            function(l) {
                                plot_grid(
                                    plotlist = l,
                                    # rel_heights = plots_rel_heights,
                                    rel_heights = plots_rel_heights[names(l)],
                                    ncol = 1,
                                    nrow = length(l),
                                    align = 'v'
                                )
                            }
                        ),
                        list(legend_plot)
                    ),
                    nrow = 1,
                    rel_widths = c(
                        left_plot_width,
                        rep(1, n_bottom_row - 1),
                        n_col - n_bottom_row
                    )
                ),
                nrow = 2,
                ncol = 1,
                rel_heights = c(sum(rows_rel_heights[-n_row]), rows_rel_heights[n_row])
            )

        }

    }

}





big_deconv_plot <- function(

    cancer_types,
    nice_names_for_figure,
    title_font_size = 14,
    title_rel_width = 0.1,
    ...

) {

    cts <- stringr::str_split_fixed(nice_names_for_figure, ' - ', 2)[, 1]

    number_of_columns <- max(table(cts[cts %in% cancer_types]))

    plot_grid(
        plotlist = lapply(
            cancer_types,
            function(ct) {

                ct_inds <- which(startsWith(nice_names_for_figure, ct))

                titles <- stringr::str_split_fixed(
                    nice_names_for_figure[startsWith(nice_names_for_figure, ct)],
                    ' - ',
                    2
                )[, 2]

                if(length(titles) == 1 && titles == '') {titles <- 'All tumours'}

                ct_plotlist <- lapply(
                    1:length(ct_inds),
                    function(i) {
                        list(
                            plots = c(
                                deconv_plots[[ct_inds[i]]]$plots[
                                    !(names(deconv_plots[[ct_inds[i]]]$plots) %in% c('heatmap', 'axis_labels'))
                                ],
                                list(
                                    heatmap = deconv_plots[[ct_inds[i]]]$plots$heatmap + labs(title = titles[i])
                                )
                            )
                        )
                    }
                )

                if(length(ct_inds) == number_of_columns) {

                    plot_grid(
                        ggplot() +
                            labs(y = ct) +
                            theme(
                                panel.background = element_rect(fill = 'white'),
                                axis.title.y = element_text(size = title_font_size)
                            ),
                        deconv_plot(
                            ct_plotlist,
                            n_row = 1,
                            n_col = number_of_columns,
                            legends = FALSE,
                            ...
                        ),
                        nrow = 1,
                        ncol = 2,
                        rel_widths = c(title_rel_width, number_of_columns)
                    )

                } else {

                    plot_grid(
                        ggplot() +
                            labs(y = ct) +
                            theme(
                                panel.background = element_rect(fill = 'white'),
                                axis.title.y = element_text(size = title_font_size)
                            ),
                        deconv_plot(
                            c(
                                ct_plotlist,
                                lapply(
                                    1:(number_of_columns - length(ct_inds)),
                                    function(i) {
                                        list(
                                            plots = setNames(
                                                lapply(
                                                    1:(length(deconv_plots[ct_inds][[1]]$plots) - 1),
                                                    function(j) blank_plot()
                                                ),
                                                names(ct_plotlist[[1]]$plots)
                                            )
                                        )
                                    }
                                )
                            ),
                            n_row = 1,
                            n_col = number_of_columns,
                            legends = FALSE,
                            ...
                        ),
                        nrow = 1,
                        ncol = 2,
                        rel_widths = c(title_rel_width, number_of_columns)
                    )

                }

            }
        ),
        nrow = length(cancer_types),
        ncol = 1
    )

}





# Function to make brackets with which to annotate the deconv results as 'pEMT' and 'CAF'
# (plotted separately, to be aligned with the deconv plot):

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

	# new_scores <- apply(
		# deconv_obj$cor_mat,
		# 2,
		# function(x) {
			# meanvec <- c(mean(x[head_genes]), mean(x[tail_genes]))
			# ifelse(sign(meanvec[1]) == sign(meanvec[2]), return(NA), return(max(meanvec)*c(1, -1)[which.max(meanvec)]))
		# }
	# )

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
