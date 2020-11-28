filter_for_groups <- function(

    grouped_data,
    group_var = 'cell_type',
    groups = c('fibroblast', 'cancer'),
    FUN = mean,
    n_req = 1,
    to_keep = c(
        'SNAI1',
        'SNAI2',
        'TWIST1',
        'VIM',
        'ZEB1',
        'ZEB2'
    )

) {

    # grouped_data should consist of columns whose 'sizes' we want to measure and compare
    # between groups, along with a single annotation column whose name should be supplied to
    # group_var and which gives the group assignment to each row.  FUN is the function that
    # you want to apply to each column (by group) to measure their 'size', and n_req is the
    # number of groups (in the character vector supplied to the groups argument) that you
    # require to be bigger than the other groups.  The function filters for those columns for
    # which this condition is satisfied, along with the elements of the argument to_keep,
    # which are put back in if they get filtered out.

    FUN <- match.fun(FUN)

    tdt(
        grouped_data[
            ,
            lapply(.SD, FUN = FUN),
            by = group_var
        ]
    )[
        ,
        n_bigger := sum(
            as.numeric(.SD[, ..groups]) > max(as.numeric(.SD[, -..groups]))
        ),
        by = id
    ][
        n_bigger >= n_req,
        unique(c(id, to_keep))
    ]

}





# For the ordering, I settled on SPIN_NH for the genes and 75th percentile for the cells.  Using SPIN
# for the cells takes far too long, because we have thousands of cells.  But also, I think that to look
# for structure among the cancer cells, e.g. the rare EMT population, it is easier to plot a separate
# heatmap of just the cancer cells, using e.g. hclust with ward.D2 (we can also use this function to do
# that, just using groups = 'cancer').  The point of this plot is to compare cancer cells and CAFs, not
# to see the structure among them.

# Other ideas for ordering:

# genes_order_fun = function(x) {as.vector(seriation::seriate(dist(t(x)), method = 'SPIN_NH')[[1]])}
# cells_order_fun = function(x) {order(-apply(x, 1, quantile, 0.75))}

# hclust() using Euclidean distance (can change method as you like):

# genes_order_fun = function(x) {hclust(dist(t(x)), method = 'ward.D2')$order}
# cells_order_fun = function(x) {hclust(dist(x), method = 'ward.D2')$order}

# Hierarchical clustering with algorithm by Gruvaeus and Wainer, which orders the clusters at each
# level so that the objects at the edge of each cluster are adjacent to that object outside the cluster
# to which it is nearest (can also use method = 'OLO' for Optimal Leaf Ordering, which is similar):

# cells_order_fun = function(x) {seriation::seriate(dist(x), method = 'GW')[[1]]$order}

# Based on correlations:

# genes_order_fun = function(x) {hclust(as.dist(1 - cor(x)), method = 'average')$order}
# cells_order_fun = function(x) {hclust(as.dist(1 - cor(t(x))), method = 'average')$order}

# We could also try biclustering...

order_spin_nh_gw_olo = function(x, method = c('GW', 'OLO')) {

    method <- match.arg(method)

    spin_nh_ordering <- as.vector(seriation::seriate(dist(t(x)), method = 'SPIN_NH')[[1]])

    .NC <- ncol(x)

    spin_nh_ordering[
        c(
            seriation::seriate(
                dist(
                    t(
                        x[
                            ,
                            spin_nh_ordering[1:(.NC %/% 2)],
                            with = FALSE
                        ]
                    )
                ),
                method = method
            )[[1]]$order,
            seriation::seriate(
                dist(
                    t(
                        x[
                            ,
                            spin_nh_ordering[(.NC %/% 2 + 1):.NC],
                            with = FALSE
                        ]
                    )
                ),
                method = method
            )[[1]]$order + x[, .NC %/% 2]
        )
    ]

}

single_cell_analysis <- function(

    single_cell_data,
    genes = NULL,
    id_var = 'id',
    group_var = 'cell_type',
    groups = c('cancer', 'fibroblast'),
    genes_filter_fun = function(x) {mean(x) > 0.2},
    cells_filter_fun = function(x) {mean(x) > 0.2},
    genes_detected_threshold = 500,
    genes_filtered = NULL,
    cells_filtered = NULL,
    genes_order = NULL,
    cells_order = NULL,
    genes_order_fun = order_spin_nh_gw_olo,
    cells_order_fun = function(x) {seriation::seriate(dist(x), method = 'GW')[[1]]$order},
    annotations = c(
        'SNAI1',
        'SNAI2',
        'TWIST1',
        'VIM',
        'ZEB1',
        'ZEB2'
    ),
    title = '',
    x_axis_titles = groups,
    default_figure_widths = c(1, 8, 1.6),
    min_figure_width = 0.2,
    h_colours = rev(colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))(50)),
    h_upper_limit = 12,
    h_legend_title = 'Expression\nlevel',
    h_legend_width = 10,
    h_legend_breaks = c(0, 3, 6, 9, 12),
    h_legend_labels = c('0' = '0', '3' = '3', '6' = '6', '9' = '9', '12' = '\u2265 12'),
    expression_summary_fun = function(x) quantile(x, 0.75),
    gd_colours = c('black', 'deeppink', 'cyan'),
    gd_limits = c(0, 8000),
    gd_breaks = c(0, 4000, 8000),
    gd_labels = c('0' = '0', '4000' = '4000', '8000' = '\u2265 8000'),
    gd_legend_title = 'Genes\ndetected',
    gd_legend_width = 10,
    gd_legend_height = 10,
    gd_legend_direction = 'vertical',
    gd_legend_margin = c(20, 0, 0, 0)

) {

    genes_filter_fun <- match.fun(genes_filter_fun)
    cells_filter_fun <- match.fun(cells_filter_fun)
    genes_order_fun <- match.fun(genes_order_fun)
    cells_order_fun <- match.fun(cells_order_fun)
    expression_summary_fun <- match.fun(expression_summary_fun)

    if(is.null(genes) & is.null(genes_filtered)) {
        stop('Please supply either genes or genes_filtered.')
    }

    if(!is.null(genes_filtered)) {
        genes <- genes_filtered
    } else {
        genes <- sort(
            unique(
                c(
                    genes[
                        apply(
                            single_cell_data[, ..genes],
                            2,
                            genes_filter_fun
                        )
                    ],
                    'SNAI1',
                    'SNAI2',
                    'TWIST1',
                    'VIM',
                    'ZEB1',
                    'ZEB2'
                )
            )
        )
    }

    single_cell_data[
        get(group_var) %in% groups,
        genes_detected := sum(as.numeric(.SD) > 0),
        by = id_var,
        .SDcols = -group_var
    ]

    setkeyv(single_cell_data, id_var)

    if(!is.null(cells_filtered)) {
        cells <- cells_filtered
    } else {
        cells <- single_cell_data[
            genes_detected >= genes_detected_threshold,
            get(id_var)[
                apply(
                    .SD[, ..genes],
                    1,
                    cells_filter_fun
                )
            ]
        ]
    }

    if(!is.null(cells_order)) {
        ordering_cells <- cells_order
    } else {
        ordering_cells <- sapply(
            groups,
            function(grp) {
                cells_order_fun(
                    single_cell_data[
                        cells
                    ][
                        get(group_var) == grp,
                        ..genes
                    ]
                )
            },
            simplify = FALSE,
            USE.NAMES = TRUE
        )
    }

    if(!is.null(genes_order)) {
        ordering_genes <- genes_order
    } else {
        ordering_genes <- genes_order_fun(
            single_cell_data[
                cells,
                ..genes
            ]
        )
    }

    heatmaps <- sapply(
        groups,
        function(grp) {
            ggplot(
                melt(
                    single_cell_data[
                        cells
                    ][
                        get(group_var) == grp,
                        c(..id_var, ..genes)
                    ],
                    id.vars = 'id', # Shouldn't this be id_var?
                    variable.name = 'gene',
                    value.name = 'expression_level'
                ),
                aes(
                    x = factor(id, levels = unique(id)[ordering_cells[[grp]]]),
                    y = factor(gene, levels = unique(gene)[ordering_genes])
                )
            ) +
                geom_raster(
                    aes(fill = expression_level)
                ) +
                scale_fill_gradientn(
                    limits = c(0, h_upper_limit),
                    breaks = h_legend_breaks,
                    labels = h_legend_labels,
                    colours = h_colours,
                    oob = scales::squish
                ) +
                scale_x_discrete(
                    expand = c(0, 0)
                ) +
                scale_y_discrete(
                    expand = c(0, 0)
                ) +
                theme(
                    axis.text = element_blank(),
                    axis.ticks.length = unit(0, 'pt'),
                    axis.ticks = element_blank(),
                    axis.title.y = element_blank(),
                    legend.position = switch(
                        (which(groups == grp) == length(groups)) + 1,
                        'none',
                        'right'
                    ),
                    legend.key.width = unit(h_legend_width, 'pt'),
                    plot.margin = unit(c(0, 5.5, 5.5, 0), 'pt')
                ) +
                labs(
                    x = x_axis_titles[groups == grp],
                    fill = h_legend_title
                )
        },
        simplify = FALSE,
        USE.NAMES = TRUE
    )

    heatmap_annotations <- genes[
        ordering_genes
    ]

    heatmap_annotations[!(heatmap_annotations %in% annotations)] <- ''

    axis_labels <- heat_map_labels_repel(
        heatmap_annotations,
        edge = 'right'
    )

    expression_summary <- sapply(
        groups,
        function(grp) {
            ggplot(
                single_cell_data[
                    cells
                ][
                    get(group_var) == grp,
                    .(id = get(id_var), ave_exp = apply(.SD, 1, expression_summary_fun)),
                    .SDcols = genes
                ],
                aes(
                    x = factor(id, levels = id[ordering_cells[[grp]]]),
                    y = 0
                )
            ) +
                geom_raster(
                    aes(fill = ave_exp)
                ) +
                scale_fill_gradientn(
                    limits = c(0, h_upper_limit),
                    colours = h_colours,
                    oob = scales::squish
                ) +
                scale_x_discrete(
                    expand = c(0, 0)
                ) +
                scale_y_discrete(
                    expand = c(0, 0)
                ) +
                theme(
                    axis.text = element_blank(),
                    axis.ticks.length = unit(0, 'pt'),
                    axis.ticks = element_blank(),
                    axis.title = element_blank(),
                    legend.position = 'none',
                    plot.margin = unit(c(0, 5.5, 5.5, 0), 'pt')
                )
        },
        simplify = FALSE,
        USE.NAMES = TRUE
    )

    genes_detected <- sapply(
        groups,
        function(grp) {
            ggplot(
                single_cell_data[
                    cells
                ][
                    get(group_var) == grp,
                    .(id = get(id_var), genes_detected)
                ],
                aes(
                    x = factor(id, levels = id[ordering_cells[[grp]]]),
                    y = 0
                )
            ) +
                geom_raster(
                    aes(fill = genes_detected)
                ) +
                scale_fill_gradientn(
                    limits = gd_limits,
                    breaks = gd_breaks,
                    labels = gd_labels,
                    oob = scales::squish,
                    colours = gd_colours
                ) +
                scale_x_discrete(
                    expand = c(0, 0)
                ) +
                scale_y_continuous(
                    expand = c(0, 0)
                ) +
                theme(
                    axis.text = element_blank(),
                    axis.ticks.length = unit(0, 'pt'),
                    axis.ticks = element_blank(),
                    axis.title = element_blank(),
                    legend.position = switch(
                        (which(groups == grp) == length(groups)) + 1,
                        'none',
                        'right'
                    ),
                    legend.box.margin = margin(gd_legend_margin),
                    legend.direction = gd_legend_direction,
                    legend.key.height = unit(gd_legend_height, 'pt'),
                    legend.key.width = unit(gd_legend_width, 'pt'),
                    plot.margin = unit(c(5.5, 5.5, 5.5, 0), 'pt')
                ) +
                guides(
                    fill = guide_colourbar(title.position = 'top')
                ) +
                labs(
                    title = switch(
                        (which(groups == grp) == 1) + 1,
                        waiver(),
                        title
                    ),
                    fill = gd_legend_title
                )
        },
        simplify = FALSE,
        USE.NAMES = TRUE
    )

    figure_widths <- switch(
        (
            single_cell_data[
                get(group_var) %in% groups,
                .(n = .N),
                by = group_var
            ][
                ,
                sum(n/sum(n) < min_figure_width) == 0
            ]
        ) + 1,
        default_figure_widths,
        c(
            nrow(single_cell_data)*default_figure_widths[1]/
                sum(default_figure_widths[2:length(default_figure_widths)]),
            sapply(
                groups,
                function(grp) {
                    single_cell_data[
                        get(group_var) == grp & genes_detected >= genes_detected_threshold,
                        .N
                    ]
                }
            )
        )
    )

    # Remove genes_detected column from single_cell_data before output:

    single_cell_data[, genes_detected := NULL]

    # Output:

    list(
        plots = c(
            setNames(
                genes_detected,
                paste0('genes_detected_', groups)
            ),
            setNames(
                expression_summary,
                paste0('expression_summary_', groups)
            ),
            list(annotations = axis_labels),
            setNames(
                heatmaps,
                paste0('heatmap_', groups)
            )
        ),
        figure_widths = figure_widths,
        cells_filtered = cells,
        genes_filtered = genes,
        ordering_cells = ordering_cells,
        ordering_genes = ordering_genes
    )

}





sc_groups <- function(

    genes,
    sc_data,
    id_var = 'id',
    group_var = 'cell_type',
    patient_var = 'patient',
    groups = c('cancer', 'fibroblast'),
    to_keep = c(
        'SNAI1',
        'SNAI2',
        'TWIST1',
        'VIM',
        'ZEB1',
        'ZEB2'
    ),
    genes_filter_fun = function(x) {log2(mean(2^x - 1) + 1) >= 2.5},
    cells_filter_fun = function(x) {mean(x) >= 0}, # I.e. no extra cell filtering
    genes_detected_threshold = 1000,
    genes_order_fun = function(x) {
        as.vector(seriation::seriate(dist(t(x)), method = 'SPIN_NH')[[1]])
    },
    cells_order_fun = function(x) {seriation::seriate(dist(x), method = 'GW')[[1]]$order},
    score_cells = TRUE,
	score_cells_nbin = 30,
	score_cells_n = 100,
    scores_filter_groups = groups,
    scores_filter_fun = function(x) {log2(mean(10*(2^x - 1)) + 1) >= 4},
    min_sig_size = 40,
    seed = NULL

) {

    # It's important that the only non-numeric (i.e. non-gene) variables in <sc_data> are <id_var>
    # and <group_var>.

    # Note that if <score_cells> is set to TRUE, scores for the cells are calculated using the
    # signature_score() function, and the cells are then ordered by these scores - in particular,
    # the <cells_order_fun> argument is ignored.  In this case the function also returns an object
	# called all_genes_filtered, which is the genes obtained by applying <scores_filter_fun> (and
	# NOT <genes_filter_fun> to the set of all genes in the dataset.  The signature genes for
	# scoring are taken from all_genes_filtered, as are the control genes used in defining scores.

    genes_filter_fun <- match.fun(genes_filter_fun)
    cells_filter_fun <- match.fun(cells_filter_fun)
    genes_order_fun <- match.fun(genes_order_fun)
    cells_order_fun <- match.fun(cells_order_fun)

    sc_data[
        get(group_var) %in% groups,
        genes_detected := sum(as.numeric(.SD) > 0),
        by = id_var,
        .SDcols = -c(group_var, patient_var)
    ]

    cells <- sc_data[
        get(group_var) %in% groups & genes_detected >= genes_detected_threshold,
        get(id_var)[
            apply(
                .SD[, ..genes],
                1,
                cells_filter_fun
            )
        ]
    ]

    setkeyv(sc_data, id_var)

    genes <- sort(
        unique(
            c(
                genes[
                    apply(
                        sc_data[
                            cells,
                            ..genes
                        ],
                        2,
                        genes_filter_fun
                    )
                ],
                to_keep
            )
        )
    )

    if(score_cells) {

        all_genes_filtered <- sc_data[
            cells
        ][
            get(group_var) %in% scores_filter_groups,
            names(.SD)[apply(.SD, 2, scores_filter_fun)],
            .SDcols = -c(id_var, group_var, patient_var, 'genes_detected')
        ]

        surviving_genes <- genes[genes %in% all_genes_filtered]

        if(length(surviving_genes) >= min_sig_size) {

            sig_genes <- surviving_genes

        } else {

            sig_genes <- c(
                surviving_genes,
                names(
                    sort(
                        sc_data[
                            cells,
                            colMeans(
                                cor(
                                    .SD[, ..surviving_genes],
                                    .SD[, -..surviving_genes]
                                )
                            ),
                            .SDcols = all_genes_filtered
                        ],
                        decreasing = TRUE
                    )[
                        1:(min_sig_size - length(surviving_genes))
                    ]
                )
            )

        }

        if(!is.null(patient_var)) {

            # There is an element of randomness in calculating the scores, presumably due to
            # sampling genes from the bins.  So we include an option to set a seed.

            if(!is.null(seed)) {set.seed(seed)}

            scores <- sc_data[
                cells,
                setNames(
                    as.data.table(
                        # scrabble::score(
						signature_score(
                            set_colnames(t(.SD[, -'id']), .SD$id),
                            # list(sig_genes),
							sig_genes,
                            # bin.control = TRUE,
                            # nbin = length(all_genes_filtered) %/% 110,
							nbin = score_cells_nbin,
							n = score_cells_n
                        ),
                        keep.rownames = TRUE
                    ),
                    c('cell_id', 'score')
                ),
                .SDcols = c('id', all_genes_filtered),
                by = patient_var
            ][, -..patient_var]

            setkey(scores, cell_id)

        } else {

            if(!is.null(seed)) {set.seed(seed)}

            scores <- setNames(
                as.data.table(
                    # scrabble::score(
					signature_score(
                        set_colnames(
                            t(
                                sc_data[
                                    cells,
                                    ..all_genes_filtered
                                ]
                            ),
                            sc_data[cells, id]
                        ),
                        # list(sig_genes),
						sig_genes,
                        # bin.control = TRUE,
                        # nbin = length(all_genes_filtered) %/% 110
                        # n = 50
						nbin = score_cells_nbin,
						n = score_cells_n
                    ),
                    keep.rownames = TRUE
                ),
                c('cell_id', 'score')
            )

            setkey(scores, cell_id)

        }

        ordering_cells <- sapply(
            groups,
            function(grp) {
                sc_data[
                    cells
                ][
                    get(group_var) == grp,
                    order(-scores[id, score])
                ]
            },
            simplify = FALSE,
            USE.NAMES = TRUE
        )

    } else {

        ordering_cells <- sapply(
            groups,
            function(grp) {
                cells_order_fun(
                    sc_data[
                        cells
                    ][
                        get(group_var) == grp,
                        ..genes
                    ]
                )
            },
            simplify = FALSE,
            USE.NAMES = TRUE
        )

    }

    ordering_genes <- genes_order_fun(
        sc_data[
            cells,
            ..genes
        ]
    )

    # Output:

    out <- list(
        data = sc_data[
            ,
            c(..id_var, ..group_var, 'genes_detected', ..patient_var, ..genes)
        ],
        cells_filtered = cells,
        genes_filtered = genes,
        ordering_cells = ordering_cells,
        ordering_genes = ordering_genes
    )

    if(score_cells) {
        out <- c(
            out,
            list(
                scores = scores,
                signature_genes_for_scoring = sig_genes,
				all_genes_filtered = all_genes_filtered
            )
        )
    }

    out

}





sc_groups_heatmap <- function(

    sc_groups_list,
    id_var = 'id',
    group_var = 'cell_type',
    genes_detected_var = 'genes_detected',
    groups = c('cancer', 'fibroblast'),
    annotations = c(
        'SNAI1',
        'SNAI2',
        'TWIST1',
        'VIM',
        'ZEB1',
        'ZEB2'
    ),
    annotations_side = c('left', 'right'),
    annotations_nudge = 0.35,
    annotations_text_size = 3,
    annotations_segment_size = 0.5,
    annotations_title = NULL,
    annotations_title_size = 14,
    annotations_margin = 3,
    annotations_force = 0.25,
    annotations_box_padding = 0.25,
    title = waiver(),
    x_axis_titles = groups,
    default_figure_widths = list(
        annotations = 1,
        cancer = 8,
        fibroblast = 1.6
    ),
    min_figure_width = 0.2,
    figure_spacing = 3,
    h_colours = rev(colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))(50)),
    h_limits = c(0, 12),
    h_legend_title = 'Expression\nlevel',
    h_legend_width = 10,
    h_legend_breaks = c(0, 3, 6, 9, 12),
    h_legend_labels = c('0' = '0', '3' = '3', '6' = '6', '9' = '9', '12' = '\u2265 12'),
    h_legend_height = 15,
    h_legend_direction = c('vertical', 'horizontal'),
    h_legend_title_position = NULL,
    h_legend_just = 'center',
    h_legend_text_size = NULL,
    h_legend_title_size = NULL,
    h_legend_box_margin = NULL,
    es_fun = mean,
    es_plot_type = c('line', 'bar'),
    es_colours = h_colours,
    es_upper_limit = 4,
    es_title = switch(is.null(es_fun) + 1, 'Average\nexpression', 'Score'),
    es_legend_width = 10,
    es_legend_breaks = c(0, 2, 4),
    es_legend_labels = c('0' = '0', '2' = '2', '4' = '\u2265 4'),
    es_min_break_interval = 1,
    es_runmean_window = NULL,
    es_line_size = 1,
    gd_colours = colorRampPalette(RColorBrewer::brewer.pal(9, "RdPu"))(50),
    gd_limits = c(0, 8000),
    gd_breaks = c(0, 4000, 8000),
    gd_labels = c('0' = '0', '4000' = '4000', '8000' = '\u2265 8000'),
    gd_legend_title = 'Genes\ndetected',
    gd_legend_width = 10,
    gd_legend_height = 10,
    gd_legend_direction = c('vertical', 'horizontal'),
    gd_legend_title_position = NULL,
    gd_legend_just = 'center',
    gd_legend_text_size = NULL,
    gd_legend_title_size = NULL,
    gd_legend_box_margin = NULL,
    # gd_legend_margin = c(25, 0, 0, 0)
    ...

) {

    # <...> is for extra arguments to guide_colourbar(), to tweak the legend format.

    annotations_side <- match.arg(annotations_side)
    gd_legend_direction <- match.arg(gd_legend_direction)
    h_legend_direction <- match.arg(h_legend_direction)

    if(!is.null(es_fun)) {
        es_fun <- match.fun(es_fun)
    }

    es_plot_type <- match.arg(es_plot_type)

    single_cell_data <- sc_groups_list$data
    cells <- sc_groups_list[[names(sc_groups_list)[startsWith(names(sc_groups_list), 'cells')]]]
    genes <- sc_groups_list[[names(sc_groups_list)[startsWith(names(sc_groups_list), 'genes')]]]
    ordering_cells <- sc_groups_list$ordering_cells
    ordering_genes <- sc_groups_list$ordering_genes

    if(is.null(groups)) {groups <- unique(single_cell_data$cell_type)}

    setkeyv(single_cell_data, id_var)

    heatmaps <- sapply(
        groups,
        function(grp) {
            ggplot(
                melt(
                    single_cell_data[
                        cells
                    ][
                        get(group_var) == grp,
                        c(..id_var, ..genes)
                    ],
                    id.vars = id_var,
                    variable.name = 'gene',
                    value.name = 'expression_level'
                ),
                aes(
                    x = factor(id, levels = unique(id)[ordering_cells[[grp]]]),
                    y = factor(gene, levels = unique(gene)[ordering_genes])
                )
            ) +
                geom_raster(
                    aes(fill = expression_level)
                ) +
                scale_fill_gradientn(
                    limits = h_limits,
                    colours = h_colours,
                    oob = scales::squish
                ) +
                scale_x_discrete(
                    expand = c(0, 0)
                ) +
                scale_y_discrete(
                    expand = c(0, 0)
                ) +
                theme(
                    axis.text = element_blank(),
                    axis.ticks.length = unit(0, 'pt'),
                    axis.ticks = element_blank(),
                    axis.title.x = switch(
                        (es_plot_type == 'bar') + 1,
                        element_blank(),
                        element_text()
                    ),
                    axis.title.y = element_blank(),
                    legend.position = 'none',
                    plot.margin = unit(
                        switch(
                            (annotations_side == 'left') + 1,
                            c(
                                0,
                                0,
                                figure_spacing,
                                switch(
                                    (which(groups == grp) == 1) + 1,
                                    figure_spacing,
                                    0
                                )
                            ),
                            c(
                                0,
                                switch(
                                    (which(groups == grp) == length(groups)) + 1,
                                    figure_spacing,
                                    0
                                ),
                                figure_spacing,
                                0
                            )
                        ),
                        'pt'
                    )
                ) +
                labs(
                    x = switch(
                        (es_plot_type == 'bar') + 1,
                        '',
                        x_axis_titles[groups == grp]
                    )
                )
        },
        simplify = FALSE,
        USE.NAMES = TRUE
    )

    heatmap_legend <- get_legend(
        ggplot(
            data.table(id = 'x', gene = 'x', expression_level = 1),
            aes(x = id, y = gene)
        ) +
            geom_raster(aes(fill = expression_level)) +
            scale_fill_gradientn(
                limits = h_limits,
                breaks = h_legend_breaks,
                labels = h_legend_labels,
                colours = h_colours,
                oob = scales::squish
            ) +
            theme(
                legend.key.width = unit(h_legend_width, 'pt'),
                legend.key.height = unit(h_legend_height, 'pt'),
                legend.direction = h_legend_direction,
                legend.justification = h_legend_just,
                legend.text = element_text(size = h_legend_text_size),
                legend.title = element_text(size = h_legend_title_size),
                legend.box.margin = h_legend_box_margin
            ) +
            guides(fill = guide_colourbar(title.position = h_legend_title_position, ...)) +
            labs(fill = h_legend_title)
    )

    heatmap_annotations <- genes[
        ordering_genes
    ]

    heatmap_annotations[!(heatmap_annotations %in% annotations)] <- ''

    # Specifying annotations_side = 'left' actually means we want the <edge> argument in
    # heat_map_labels_repel to be 'right', and vice versa, hence the call to mapvalues()
    # in the following.

    axis_labels <- heat_map_labels_repel(
        heatmap_annotations,
        edge = plyr::mapvalues(
            annotations_side,
            c('left', 'right'),
            c('right', 'left'),
            warn_missing = FALSE
        ),
        nudge = annotations_nudge,
        text_size = annotations_text_size,
        axis_title = annotations_title,
        axis_title_size = annotations_title_size,
        title_edge_margin = annotations_margin,
        segment_size = annotations_segment_size,
        force = annotations_force,
        box.padding = annotations_box_padding
    )

    if(es_plot_type == 'bar') {

        if(is.null(es_upper_limit)) {
            es_upper_limit <- single_cell_data[
                cells,
                max(apply(.SD, 1, es_fun)),
                .SDcols = genes
            ]
        }

        if(is.null(es_fun)) {

            expression_summary_data <- sapply(
                groups,
                function(grp) {
                    sc_groups_list$scores[
                        single_cell_data[cells][get(group_var) == grp, get(id_var)],
                        .(
                            id = cell_id,
                            es = score
                        )
                    ][
                        ordering_cells[[grp]]
                    ]
                },
                simplify = FALSE,
                USE.NAMES = TRUE
            )

        } else {

            expression_summary_data <- sapply(
                groups,
                function(grp) {
                    single_cell_data[
                        cells
                    ][
                        get(group_var) == grp,
                        .(id = get(id_var), es = apply(.SD, 1, es_fun)),
                        .SDcols = genes
                    ][
                        ordering_cells[[grp]]
                    ]
                },
                simplify = FALSE,
                USE.NAMES = TRUE
            )

        }

        expression_summary <- sapply(
            groups,
            function(grp) {
                ggplot(
                    expression_summary_data,
                    aes(
                        x = factor(id, levels = id),
                        y = 0
                    )
                ) +
                    geom_raster(
                        aes(fill = es)
                    ) +
                    scale_fill_gradientn(
                        limits = c(0, es_upper_limit),
                        colours = es_colours,
                        oob = scales::squish
                    ) +
                    scale_x_discrete(
                        expand = c(0, 0)
                    ) +
                    scale_y_discrete(
                        expand = c(0, 0)
                    ) +
                    theme(
                        axis.text = element_blank(),
                        axis.ticks.length = unit(0, 'pt'),
                        axis.ticks = element_blank(),
                        axis.title = element_blank(),
                        legend.position = 'none',
                        plot.margin = unit(
                            switch(
                                (annotations_side == 'left') + 1,
                                c(
                                    0,
                                    0,
                                    figure_spacing,
                                    switch(
                                        (which(groups == grp) == 1) + 1,
                                        figure_spacing,
                                        0
                                    )
                                ),
                                c(
                                    0,
                                    switch(
                                        (which(groups == grp) == length(groups)) + 1,
                                        figure_spacing,
                                        0
                                    ),
                                    figure_spacing,
                                    0
                                )
                            ),
                            'pt'
                        )
                    )
            },
            simplify = FALSE,
            USE.NAMES = TRUE
        )

    } else { # Then es_plot_type == 'line'.

        if(is.null(es_runmean_window)) {
            es_runmean_window <- nrow(single_cell_data)/20
        }

        if(is.null(es_fun)) {

            setkey(sc_groups_list$scores, cell_id)

            expression_summary_data <- sapply(
                groups,
                function(grp) {
                    sc_groups_list$scores[
                        single_cell_data[cells][get(group_var) == grp, get(id_var)],
                        .(
                            id = cell_id,
                            es = score
                        )
                    ][
                        ordering_cells[[grp]]
                    ][
                        ,
                        es := caTools::runmean(es, es_runmean_window)
                    ]
                },
                simplify = FALSE,
                USE.NAMES = TRUE
            )

        } else {

            expression_summary_data <- sapply(
                groups,
                function(grp) {
                    single_cell_data[
                        cells
                    ][
                        get(group_var) == grp,
                        .(
                            id = get(id_var),
                            es = apply(.SD, 1, es_fun)
                        ),
                        .SDcols = genes
                    ][
                        ordering_cells[[grp]]
                    ][
                        ,
                        es := caTools::runmean(es, es_runmean_window)
                    ]
                },
                simplify = FALSE,
                USE.NAMES = TRUE
            )

        }

        es_upper_limit <- ceiling(
            2*max(
                sapply(
                    expression_summary_data,
                    function(x) x[, max(es)]
                )
            )
        )/2

        # es_lower_limit <- -ceiling(
        #     2*max(
        #         sapply(
        #             expression_summary_data,
        #             function(x) x[, abs(min(es))]
        #         )
        #     )
        # )/2

        es_lower_limit <- min(
            floor(
                2*min(
                    sapply(
                        expression_summary_data,
                        function(x) x[, min(es)]
                    )
                )
            )/2,
            0
        )

        # es_linegraph_breaks <- switch(
        #     (es_upper_limit < 2) + 1,
        #     seq(
        #         from = es_lower_limit,
        #         to = es_upper_limit,
        #         by = floor(log2(es_upper_limit))
        #     ),
        #     seq(
        #         from = es_lower_limit,
        #         to = es_upper_limit,
        #         by = floor(log2(2*es_upper_limit))/2
        #     )
        # )

        dist_btw_breaks <- (((es_upper_limit - es_lower_limit)/es_min_break_interval) %/%
                                c(2, 3))*es_min_break_interval

        n_breaks <- 1 +
            (es_upper_limit %/% dist_btw_breaks) +
            (abs(es_lower_limit) %/% dist_btw_breaks)

        dist_btw_breaks <- switch(
            (n_breaks[1] == n_breaks[2]) + 1,
            dist_btw_breaks[which.max(n_breaks[n_breaks <= 4])],
            max(dist_btw_breaks)
        )

        es_linegraph_breaks <- dist_btw_breaks*seq(
            from = -(abs(es_lower_limit) %/% dist_btw_breaks),
            to = es_upper_limit %/% dist_btw_breaks,
            by = 1
        )

        expression_summary <- sapply(
            groups,
            function(grp) {
                ggplot(
                    expression_summary_data[[grp]],
                    aes(
                        x = factor(id, levels = id),
                        y = es,
                        group = 1
                    )
                ) +
                    geom_line(
                        colour = 'mediumorchid1',
                        size = es_line_size
                    ) +
                    scale_y_continuous(
                        limits = c(es_lower_limit, es_upper_limit),
                        expand = c(0, 0),
                        breaks = es_linegraph_breaks,
                        position = annotations_side
                    ) +
                    scale_x_discrete(
                        expand = c(0, 0)
                    ) +
                    theme(
                        axis.text.x = element_blank(),
                        axis.text.y = switch(
                            (annotations_side == 'left') + 1,
                            switch( # Annotations on the right - display y axis text only on rightmost panel
                                (which(groups == grp) == length(groups)) + 1,
                                element_blank(),
                                element_text()
                            ),
                            switch( # Annotations on the left - display y axis text only on leftmost panel
                                (which(groups == grp) == 1) + 1,
                                element_blank(),
                                element_text()
                            )
                        ),
                        axis.title.y = element_blank(),
                        axis.ticks.x = element_blank(),
                        panel.background = element_blank(),
                        panel.border = element_rect(fill = NA, colour = 'black', size = 0.5),
                        panel.grid = element_blank(),
                        plot.margin = unit(
                            switch(
                                (annotations_side == 'left') + 1,
                                c(
                                    0,
                                    0,
                                    0,
                                    switch(
                                        (which(groups == grp) == 1) + 1,
                                        figure_spacing,
                                        0
                                    )
                                ),
                                c(
                                    0,
                                    switch(
                                        (which(groups == grp) == length(groups)) + 1,
                                        figure_spacing,
                                        0
                                    ),
                                    0,
                                    0
                                )
                            ),
                            'pt'
                        )
                    ) +
                    labs(
                        x = x_axis_titles[groups == grp]
                    )
            },
            simplify = FALSE,
            USE.NAMES = TRUE
        )

    }

    genes_detected <- sapply(
        groups,
        function(grp) {
            ggplot(
                single_cell_data[
                    cells
                ][
                    get(group_var) == grp,
                    .(id = get(id_var), genes_detected = get(genes_detected_var))
                ],
                aes(
                    x = factor(id, levels = id[ordering_cells[[grp]]]),
                    y = 0
                )
            ) +
                geom_raster(
                    aes(fill = genes_detected)
                ) +
                scale_fill_gradientn(
                    limits = gd_limits,
                    breaks = gd_breaks,
                    labels = gd_labels,
                    oob = scales::squish,
                    colours = gd_colours
                ) +
                scale_x_discrete(
                    expand = c(0, 0)
                ) +
                scale_y_continuous(
                    expand = c(0, 0)
                ) +
                theme(
                    axis.text = element_blank(),
                    axis.ticks.length = unit(0, 'pt'),
                    axis.ticks = element_blank(),
                    axis.title = element_blank(),
                    legend.position = 'none',
                    plot.margin = unit(
                        switch(
                            (annotations_side == 'left') + 1,
                            c(
                                figure_spacing,
                                0,
                                figure_spacing,
                                switch(
                                    (which(groups == grp) == 1) + 1,
                                    figure_spacing,
                                    0
                                )
                            ),
                            c(
                                figure_spacing,
                                switch(
                                    (which(groups == grp) == length(groups)) + 1,
                                    figure_spacing,
                                    0
                                ),
                                figure_spacing,
                                0
                            )
                        ),
                        'pt'
                    )
                ) +
                labs(
                    title = switch(
                        (which(groups == grp) == 1) + 1,
                        waiver(),
                        title
                    )
                )
        },
        simplify = FALSE,
        USE.NAMES = TRUE
    )

    genes_detected_legend <- get_legend(
        ggplot(
            data.table(id = 'x', genes_detected = 1),
            aes(x = id, y = 0)
        ) +
            geom_raster(
                aes(fill = genes_detected)
            ) +
            scale_fill_gradientn(
                limits = gd_limits,
                breaks = gd_breaks,
                labels = gd_labels,
                oob = scales::squish,
                colours = gd_colours
            ) +
            theme(
                legend.key.height = unit(gd_legend_height, 'pt'),
                legend.key.width = unit(gd_legend_width, 'pt'),
                legend.direction = gd_legend_direction,
                legend.justification = gd_legend_just,
                legend.text = element_text(size = gd_legend_text_size),
                legend.title = element_text(size = gd_legend_title_size),
                legend.box.margin = gd_legend_box_margin
            ) +
            guides(fill = guide_colourbar(title.position = gd_legend_title_position, ...)) +
            labs(fill = gd_legend_title)
    )

    figure_widths <- switch(
        (
            single_cell_data[
                # get(group_var) %in% groups,
				cells,
                .(n = .N),
                by = group_var
            ][
                ,
                sum(n/sum(n) < min_figure_width) == 0
            ]
        ) + 1,
        default_figure_widths,
        c(
            list(
                # annotations = nrow(single_cell_data)*default_figure_widths$annotations/
				annotations = length(cells)*default_figure_widths$annotations/
                    sum(unlist(default_figure_widths[groups]))
            ),
            single_cell_data[
                cells,
            # ][
                # ,
                .(N = .N),
                keyby = group_var
            ][
                groups,
                setNames(as.list(N), groups)
            ]
        )
    )

    out <- list(
        plots = c(
            list(
                genes_detected = setNames(
                    genes_detected,
                    groups
                )
            ),
            list(genes_detected_legend = genes_detected_legend),
            list(
                expression_summary = setNames(
                    expression_summary,
                    groups
                )
            ),
            list(
                heatmaps = setNames(
                    heatmaps,
                    groups
                )
            ),
            list(
                heatmap_legend = heatmap_legend,
                annotations = axis_labels
            )
        ),
        groups = groups,
        figure_widths = figure_widths,
        annotations_side = annotations_side,
        es_plot_type = es_plot_type,
        es_title = es_title
    )

    if(es_plot_type == 'line') {
        out <- c(out, list(es_upper_limit = es_upper_limit))
    }

    out

}





cowplot_sc <- function(

    sc_groups_figures,
    heights = switch(
        (sc_groups_figures$es_plot_type == 'bar') + 1,
        c(2, 20, 4),
        c(2, 1, 20)
    ),
    legend_space = 0.1,
    legend_rel_heights = c(1, 50),
    es_x_axis_title_vjust = 1,
    es_x_axis_title_size = 11,
    es_y_axis_text_size = NULL,
    es_y_axis_ticks_size = NULL,
    es_y_axis_ticks_length = 3,
    es_y_axis_title = sc_groups_figures$es_title,
    es_y_axis_title_angle = switch(
        (sc_groups_figures$annotations_side == 'left') + 1,
        -90,
        90
    ),
    es_y_axis_title_hjust = switch(
        (es_y_axis_title_angle == 0) + 1,
        0.5,
        switch(
            (sc_groups_figures$annotations_side == 'left') + 1,
            0,
            1
        )
    ),
    es_y_axis_title_vjust = switch(
        (es_y_axis_title_angle %in% c(-90, 90)) + 1,
        0.5,
        -0.75
    ),
    es_y_axis_title_xpos = switch(
        (es_y_axis_title_angle %in% c(-90, 90)) + 1,
        switch(
            (sc_groups_figures$annotations_side == 'left') + 1,
            0.1,
            0.9
        ),
        switch(
            (sc_groups_figures$annotations_side == 'left') + 1,
            0,
            1
        )
    ),
    es_y_axis_title_size = 11,
    es_panel_border_size = NULL,
    x_axis_titles_space = 0.2

) {

    if(sc_groups_figures$es_plot_type == 'bar') {

        # Since the code for this case is out of date it probably won't work, so let's just
        # use stop():

        stop("Code for es_plot_type == 'bar' case is out of date and won't work.")

        plot_grid(
            plot_grid(
                plot_grid(
                    plotlist = list(
                        blank_plot(),
                        sc_groups_figures$plots$genes_detected_cancer,
                        sc_groups_figures$plots$genes_detected_fibroblast +
                            theme(legend.position = 'none') +
                            labs(title = '')
                    ),
                    ncol = 3,
                    nrow = 1,
                    rel_widths = sc_groups_figures$figure_widths,
                    align = 'h'
                ),
                legend_genes_detected,
                ncol = 2,
                nrow = 1,
                rel_widths = c(1, legend_space)
            ),
            plot_grid(
                blank_plot(),
                sc_groups_figures$plots$expression_summary_cancer,
                sc_groups_figures$plots$expression_summary_fibroblast,
                blank_plot(),
                ncol = 4,
                nrow = 1,
                rel_widths = c(
                    sc_groups_figures$figure_widths,
                    sum(sc_groups_figures$figure_widths)*legend_space
                ),
                align ='h'
            ),
            plot_grid(
                plot_grid(
                    plotlist = list(
                        sc_groups_figures$plots$annotations,
                        sc_groups_figures$plots$heatmap_cancer,
                        sc_groups_figures$plots$heatmap_fibroblast +
                            theme(
                                legend.position = 'none'
                            )
                    ),
                    ncol = 3,
                    nrow = 1,
                    rel_widths = sc_groups_figures$figure_widths,
                    align = 'h'
                ),
                legend_heatmap,
                ncol = 2,
                nrow = 1,
                rel_widths = c(1, legend_space)
            ),
            ncol = 1,
            nrow = 3,
            rel_heights = heights
        )

    } else { # Then sc_groups_figures$es_plot_type == 'line'.

        legend_space_scaled <- switch(
            (legend_space == 0) + 1,
            sum(unlist(sc_groups_figures$figure_widths))*legend_space,
            NULL
        )

        cols_without_legends <- length(sc_groups_figures$groups) + 1

        total_cols <- switch(
            (legend_space == 0) + 1,
            cols_without_legends + 1,
            cols_without_legends
        )

        rel_widths_all_cols <- switch(
            (sc_groups_figures$annotations_side == 'left') + 1,
            c(
                legend_space_scaled,
                unlist(
                    sc_groups_figures$figure_widths[
                        c(sc_groups_figures$groups, 'annotations')
                    ]
                )
            ),
            c(
                unlist(
                    sc_groups_figures$figure_widths[
                        c('annotations', sc_groups_figures$groups)
                    ]
                ),
                legend_space_scaled
            )
        )

        if(legend_space > 0) {

            genes_detected_panel_plotlist <- c(
                list(blank_plot()),
                sc_groups_figures$plots$genes_detected,
                list(blank_plot())
            )

        } else if(sc_groups_figures$annotations_side == 'right') {

            genes_detected_panel_plotlist <- c(
                sc_groups_figures$plots$genes_detected,
                list(blank_plot())
            )

        } else if(sc_groups_figures$annotations_side == 'left') {

            genes_detected_panel_plotlist <- c(
                list(blank_plot()),
                sc_groups_figures$plots$genes_detected
            )

        }

        genes_detected_panel <- plot_grid(
            plotlist = genes_detected_panel_plotlist,
            ncol = total_cols,
            nrow = 1,
            rel_widths = rel_widths_all_cols,
            align = 'h'
        )

        if(sc_groups_figures$annotations_side == 'left') {

            if(legend_space > 0) {

                heatmap_and_legends_panel <- plot_grid(
                    plot_grid(
                        plotlist = c(
                            list(sc_groups_figures$plots$annotations),
                            sc_groups_figures$plots$heatmaps
                        ),
                        ncol = cols_without_legends,
                        nrow = 1,
                        rel_widths = unlist(
                            sc_groups_figures$figure_widths[
                                c('annotations', sc_groups_figures$groups)
                            ]
                        ),
                        align = 'h'
                    ),
                    plot_grid(
                        sc_groups_figures$plots$genes_detected_legend,
                        sc_groups_figures$plots$heatmap_legend,
                        ncol = 1,
                        nrow = 2,
                        rel_heights = legend_rel_heights
                    ),
                    ncol = 2,
                    nrow = 1,
                    rel_widths = c(1, legend_space)
                )

            } else {

                heatmap_and_legends_panel <- plot_grid(
                    plotlist = c(
                        list(sc_groups_figures$plots$annotations),
                        sc_groups_figures$plots$heatmaps
                    ),
                    ncol = cols_without_legends,
                    nrow = 1,
                    rel_widths = unlist(
                        sc_groups_figures$figure_widths[
                            c('annotations', sc_groups_figures$groups)
                        ]
                    ),
                    align = 'h'
                )

            }

        } else {

            if(legend_space > 0) {

                heatmap_and_legends_panel <- plot_grid(
                    plot_grid(
                        sc_groups_figures$plots$genes_detected_legend,
                        sc_groups_figures$plots$heatmap_legend,
                        ncol = 1,
                        nrow = 2,
                        rel_heights = legend_rel_heights
                    ),
                    plot_grid(
                        plotlist = c(
                            sc_groups_figures$plots$heatmaps,
                            list(sc_groups_figures$plots$annotations)
                        ),
                        ncol = cols_without_legends,
                        nrow = 1,
                        rel_widths = unlist(
                            sc_groups_figures$figure_widths[
                                c(sc_groups_figures$groups, 'annotations')
                            ]
                        ),
                        align = 'h'
                    ),
                    ncol = 2,
                    nrow = 1,
                    rel_widths = c(legend_space, 1)
                )

            } else {

                heatmap_and_legends_panel <- plot_grid(
                    plotlist = c(
                        sc_groups_figures$plots$heatmaps,
                        list(sc_groups_figures$plots$annotations)
                    ),
                    ncol = cols_without_legends,
                    nrow = 1,
                    rel_widths = unlist(
                        sc_groups_figures$figure_widths[
                            c(sc_groups_figures$groups, 'annotations')
                        ]
                    ),
                    align = 'h'
                )

            }

        }

        # Linegraph panel:

        # The lineplots and their x axes are plotted on two separate rows, which are then
        # combined.  The following makes lists for each of these two rows, depending on
        # <annotations_side> and <legend_space>.

        linegraphs_x_axes_list <- lapply(
            sc_groups_figures$groups,
            function(grp) {
                plot_grid(
                    get_x_axis(sc_groups_figures$plots$expression_summary[[grp]])
                ) +
                    draw_label(
                        # grp,
                        sc_groups_figures$plots$expression_summary[[grp]]$labels$x,
                        size = es_x_axis_title_size,
                        y = 1,
                        vjust = es_x_axis_title_vjust
                    )
            }
        )

        if(sc_groups_figures$annotations_side == 'left') {

            linegraphs_list <- c(
                list(
                    ggdraw(
                        get_y_axis(
                            sc_groups_figures$plots$expression_summary[[1]] +
                                theme(
                                    axis.text.y = element_text(size = es_y_axis_text_size),
                                    axis.ticks.y = element_line(size = es_y_axis_ticks_size),
                                    axis.ticks.length = unit(es_y_axis_ticks_length, 'pt')
                                ),
                            position = 'left'
                        )
                    ) +
                        draw_label(
                            es_y_axis_title,
                            x = es_y_axis_title_xpos,
                            y = 0.5,
                            angle = es_y_axis_title_angle,
                            hjust = es_y_axis_title_hjust,
                            vjust = es_y_axis_title_vjust,
                            size = es_y_axis_title_size
                        )
                ),
                lapply(
                    sc_groups_figures$plots$expression_summary,
                    function(g) {
                        g + theme(
                            axis.text.y = element_blank(),
                            axis.title = element_blank(),
                            axis.ticks.y = element_blank(),
                            axis.ticks.length = unit(0, 'pt'),
                            panel.border = element_rect(size = es_panel_border_size)
                        )
                    }
                )
            )

            linegraphs_x_axes_list <- c(
                list(blank_plot()),
                linegraphs_x_axes_list
            )

            if(legend_space > 0) {

                linegraphs_list <- c(
                    linegraphs_list,
                    list(blank_plot())
                )

                linegraphs_x_axes_list <- c(
                    linegraphs_x_axes_list,
                    list(blank_plot())
                )

            }

        } else {

            linegraphs_list <- c(
                lapply(
                    sc_groups_figures$plots$expression_summary,
                    function(g) {
                        g + theme(
                            axis.text.y = element_blank(),
                            axis.title = element_blank(),
                            axis.ticks.y = element_blank(),
                            axis.ticks.length = unit(0, 'pt'),
                            panel.border = element_rect(size = es_panel_border_size)
                        )
                    }
                ),
                list(
                    ggdraw(
                        get_y_axis(
                            sc_groups_figures$plots$expression_summary[[
                                length(sc_groups_figures$groups)
                            ]] +
                                theme(
                                    axis.text.y = element_text(size = es_y_axis_text_size),
                                    axis.ticks.y = element_line(size = es_y_axis_ticks_size),
                                    axis.ticks.length = unit(es_y_axis_ticks_length, 'pt')
                                ),
                            position = 'right'
                        )
                    ) +
                        draw_label(
                            es_y_axis_title,
                            x = es_y_axis_title_xpos,
                            y = 0.5,
                            angle = es_y_axis_title_angle,
                            hjust = es_y_axis_title_hjust,
                            vjust = es_y_axis_title_vjust,
                            size = es_y_axis_title_size
                        )
                )
            )

            linegraphs_x_axes_list <- c(
                linegraphs_x_axes_list,
                list(blank_plot())
            )

            if(legend_space > 0) {

                linegraphs_list <- c(
                    list(blank_plot()),
                    linegraphs_list
                )

                linegraphs_x_axes_list <- c(
                    list(blank_plot()),
                    linegraphs_x_axes_list
                )

            }

        }

        # Now make the linegraphs panel by making two separate rows - one for the graphs and
        # one for their x axes - and then combining these:

        linegraph_panel <- plot_grid(
            plot_grid( # This call to plot_grid is for the actual graphs, without x axis titles
                plotlist = linegraphs_list,
                ncol = total_cols,
                nrow = 1,
                rel_widths = rel_widths_all_cols
            ),
            plot_grid( # This call to plot_grid is just for the x axis titles
                plotlist = linegraphs_x_axes_list,
                ncol = total_cols,
                nrow = 1,
                rel_widths = rel_widths_all_cols
            ),
            ncol = 1,
            nrow = 2,
            rel_heights = c(1, x_axis_titles_space)
        )

        plot_grid(
            genes_detected_panel,
            heatmap_and_legends_panel,
            linegraph_panel,
            ncol = 1,
            nrow = 3,
            rel_heights = heights
        )

    }

}





cowplot_sc_rare_shared <- function(

    sc_groups_figures,
    heights = c(2, 20, 4),
    legend_space = 0.1,
    legend_rel_heights = c(1, 50),
    es_y_axis_title = sc_groups_figures$es_title,
    es_y_axis_title_angle = switch(
        (sc_groups_figures$annotations_side == 'left') + 1,
        -90,
        90
    ),
    es_y_axis_title_hjust = switch(
        (es_y_axis_title_angle == 0) + 1,
        0.5,
        switch(
            (sc_groups_figures$annotations_side == 'left') + 1,
            0,
            1
        )
    ),
    es_y_axis_title_vjust = switch(
        (es_y_axis_title_angle %in% c(-90, 90)) + 1,
        0.5,
        -0.75
    ),
    es_y_axis_title_xpos = switch(
        (es_y_axis_title_angle %in% c(-90, 90)) + 1,
        switch(
            (sc_groups_figures$annotations_side == 'left') + 1,
            0.1,
            0.9
        ),
        switch(
            (sc_groups_figures$annotations_side == 'left') + 1,
            0,
            1
        )
    ),
    heatmap_title = '',
    heatmap_title_vjust = 1.2

) {

    legend_space_scaled <- switch(
        (legend_space == 0) + 1,
        sum(unlist(sc_groups_figures$figure_widths))*legend_space,
        NULL
    )

    total_cols <- switch((legend_space == 0) + 1, 3, 2)

    rel_widths_all_cols <- switch(
        (sc_groups_figures$annotations_side == 'left') + 1,
        c(
            legend_space_scaled,
            unlist(
                sc_groups_figures$figure_widths[
                    c('heatmap', 'annotations')
                ]
            )
        ),
        c(
            unlist(
                sc_groups_figures$figure_widths[
                    c('annotations', 'heatmap')
                ]
            ),
            legend_space_scaled
        )
    )

    # Genes detected panel:

    if(legend_space > 0) {

        genes_detected_panel_plotlist <- c(
            list(blank_plot()),
            list(sc_groups_figures$plots$genes_detected),
            list(blank_plot())
        )

    } else if(sc_groups_figures$annotations_side == 'right') {

        genes_detected_panel_plotlist <- c(
            list(sc_groups_figures$plots$genes_detected),
            list(blank_plot())
        )

    } else if(sc_groups_figures$annotations_side == 'left') {

        genes_detected_panel_plotlist <- c(
            list(blank_plot()),
            list(sc_groups_figures$plots$genes_detected)
        )

    }

    genes_detected_panel <- plot_grid(
        plotlist = genes_detected_panel_plotlist,
        ncol = total_cols,
        nrow = 1,
        rel_widths = rel_widths_all_cols,
        align = 'h'
    )

    # Heatmaps and legends panel:

    if(sc_groups_figures$annotations_side == 'left') {

        if(legend_space > 0) {

            heatmap_and_legends_panel <- plot_grid(
                plot_grid(
                    plotlist = list(
                        sc_groups_figures$plots$annotations$rare +
                            theme(axis.title = element_blank()),
                        sc_groups_figures$plots$heatmaps$rare,
                        sc_groups_figures$plots$annotations$shared +
                            theme(axis.title = element_blank()),
                        sc_groups_figures$plots$heatmaps$shared
                    ),
                    ncol = 2,
                    nrow = 2,
                    rel_widths = unlist(
                        sc_groups_figures$figure_widths[
                            c('annotations', 'heatmap')
                        ]
                    ),
                    rel_heights = c(
                        length(sc_groups_figures$genes$rare),
                        length(sc_groups_figures$genes$shared)
                    ),
                    align = 'h'
                ),
                plot_grid(
                    sc_groups_figures$plots$genes_detected_legend,
                    sc_groups_figures$plots$heatmap_legend,
                    ncol = 1,
                    nrow = 2,
                    rel_heights = legend_rel_heights
                ),
                ncol = 2,
                nrow = 1,
                rel_widths = c(1, legend_space)
            ) +
                draw_label(
                    heatmap_title,
                    x = 0,
                    y = 0.5,
                    angle = 90,
                    vjust = heatmap_title_vjust,
                    size = 14
                )

        } else {

            heatmap_and_legends_panel <- plot_grid(
                plotlist = list(
                    sc_groups_figures$plots$annotations$rare +
                        theme(axis.title = element_blank()),
                    sc_groups_figures$plots$heatmaps$rare,
                    sc_groups_figures$plots$annotations$shared +
                        theme(axis.title = element_blank()),
                    sc_groups_figures$plots$heatmaps$shared
                ),
                ncol = 2,
                nrow = 2,
                rel_widths = unlist(
                    sc_groups_figures$figure_widths[
                        c('annotations', 'heatmap')
                    ]
                ),
                rel_heights = c(
                    length(sc_groups_figures$genes$rare),
                    length(sc_groups_figures$genes$shared)
                ),
                align = 'h'
            ) +
                draw_label(
                    heatmap_title,
                    x = 0,
                    y = 0.5,
                    angle = 90,
                    vjust = heatmap_title_vjust,
                    size = 14
                )

        }

    } else {

        if(legend_space > 0) {

            heatmap_and_legends_panel <- plot_grid(
                plot_grid(
                    sc_groups_figures$plots$genes_detected_legend,
                    sc_groups_figures$plots$heatmap_legend,
                    ncol = 1,
                    nrow = 2,
                    rel_heights = legend_rel_heights
                ),
                plot_grid(
                    plotlist = list(
                        sc_groups_figures$plots$heatmaps$rare,
                        sc_groups_figures$plots$annotations$rare +
                            theme(axis.title = element_blank()),
                        sc_groups_figures$plots$heatmaps$shared,
                        sc_groups_figures$plots$annotations$shared +
                            theme(axis.title = element_blank())
                    ),
                    ncol = 2,
                    nrow = 2,
                    rel_widths = unlist(
                        sc_groups_figures$figure_widths[
                            c('heatmap', 'annotations')
                        ]
                    ),
                    rel_heights = c(
                        length(sc_groups_figures$genes$rare),
                        length(sc_groups_figures$genes$shared)
                    ),
                    align = 'h'
                ),
                ncol = 2,
                nrow = 1,
                rel_widths = c(legend_space, 1)
            ) +
                draw_label(
                    heatmap_title,
                    x = 1,
                    y = 0.5,
                    angle = -90,
                    vjust = heatmap_title_vjust,
                    size = 14
                )

        } else {

            heatmap_and_legends_panel <- plot_grid(
                plotlist = list(
                    sc_groups_figures$plots$heatmaps$rare,
                    sc_groups_figures$plots$annotations$rare +
                        theme(axis.title = element_blank()),
                    sc_groups_figures$plots$heatmaps$shared,
                    sc_groups_figures$plots$annotations$shared +
                        theme(axis.title = element_blank())
                ),
                ncol = 2,
                nrow = 2,
                rel_widths = unlist(
                    sc_groups_figures$figure_widths[
                        c('heatmap', 'annotations')
                    ]
                ),
                rel_heights = c(
                    length(sc_groups_figures$genes$rare),
                    length(sc_groups_figures$genes$shared)
                ),
                align = 'h'
            ) +
                draw_label(
                    heatmap_title,
                    x = 1,
                    y = 0.5,
                    angle = -90,
                    vjust = heatmap_title_vjust,
                    size = 14
                )

        }

    }

    # Linegraph panel:

    if(sc_groups_figures$annotations_side == 'left') {

        linegraphs_list <- c(
            list(
                ggdraw(
                    get_y_axis(
                        sc_groups_figures$plots$expression_summary,
                        position = 'left'
                    )
                ) +
                    draw_label(
                        es_y_axis_title,
                        x = es_y_axis_title_xpos,
                        y = 0.5,
                        angle = es_y_axis_title_angle,
                        hjust = es_y_axis_title_hjust,
                        vjust = es_y_axis_title_vjust,
                        size = 11
                    )
            ),
            list(
                sc_groups_figures$plots$expression_summary + theme(
                    axis.text.y = element_blank(),
                    axis.title = element_blank(),
                    axis.ticks.y = element_blank(),
                    axis.ticks.length = unit(0, 'pt')
                )
            )
        )

        if(legend_space > 0) {

            linegraphs_list <- c(
                linegraphs_list,
                list(blank_plot())
            )

        }

    } else {

        linegraphs_list <- c(
            list(
                sc_groups_figures$plots$expression_summary + theme(
                    axis.text.y = element_blank(),
                    axis.title = element_blank(),
                    axis.ticks.y = element_blank(),
                    axis.ticks.length = unit(0, 'pt')
                )
            ),
            list(
                ggdraw(
                    get_y_axis(
                        sc_groups_figures$plots$expression_summary,
                        position = 'right'
                    )
                ) +
                    draw_label(
                        es_y_axis_title,
                        x = es_y_axis_title_xpos,
                        y = 0.5,
                        angle = es_y_axis_title_angle,
                        hjust = es_y_axis_title_hjust,
                        vjust = es_y_axis_title_vjust,
                        size = 11
                    )
            )
        )

        if(legend_space > 0) {

            linegraphs_list <- c(
                list(blank_plot()),
                linegraphs_list
            )

        }

    }

    # Now make the linegraphs panel by making two separate rows - one for the graphs and
    # one for their x axes - and then combining these:

    linegraph_panel <- plot_grid(
        plotlist = linegraphs_list,
        ncol = total_cols,
        nrow = 1,
        rel_widths = rel_widths_all_cols
    )

    plot_grid(
        genes_detected_panel,
        heatmap_and_legends_panel,
        linegraph_panel,
        ncol = 1,
        nrow = 3,
        rel_heights = heights
    )

}





score_dist <- function(

    sc_data,
    cells,
    genes_filter_fun,
    scores_table,
    seed = NULL,
    control_method = c('shuffle', 'random_genes'),
    n_random_genes = 30,
    n_replicate = 50,
    null_cor_threshold = 0.95,
    plot_title = ''

) {

    genes_filter_fun <- match.fun(genes_filter_fun)

    control_method <- match.arg(control_method)

    setkey(sc_data, id)

    # In the following, can we be sure that sc_data has id, patient and cell_type columns?

    all_genes_filtered <- sc_data[
        cells,
        names(.SD)[apply(.SD, 2, genes_filter_fun)],
        .SDcols = -c('id', 'patient', 'cell_type')
    ]

    if(!is.null(seed)) set.seed(seed)

    if(control_method == 'shuffle') {

        # The following computes the null distribution by shuffling the elements of
        # each gene vector.

        dt <- sc_data[
            cells,
            .(
                gene = all_genes_filtered,
                emt_cor = cor(.SD, scores_table[id, score])[, 1],
                null_cor = cor(
                    as.data.table(lapply(.SD, sample)),
                    scores_table[id, score]
                )[, 1]
            ),
            by = patient,
            .SDcols = all_genes_filtered
        ][
            ,
            sapply(
                .SD,
                # The following function is basically the mean, but we use it instead
                # of mean() to stop us losing genes which are all zero in a subset of
                # patients (and hence give NAs for those patients):
                function(cor_type) {
                    sum(cor_type[!is.na(cor_type)])/.N
                },
                simplify = FALSE,
                USE.NAMES = TRUE
            ),
            by = gene,
            .SDcols = -'patient'
        ]

        minx <- min(dt[, -'gene'])
        maxx <- max(dt[, -'gene'])

        xvals <- seq(
            from = minx,
            to = maxx,
            length.out = nrow(dt)
        )

        pdf_funs <- sapply(
            names(dt[, -'gene']),
            function(cor_type) {

                density_fun <- approxfun(density(dt[[cor_type]]))

                df_fun <- approxfun(
                    data.table(
                        x = xvals
                    )[
                        ,
                        y := switch(
                            is.na(density_fun(x)) + 1,
                            density_fun(x),
                            0
                        ),
                        by = x
                    ]
                )

            },
            simplify = FALSE,
            USE.NAMES = TRUE
        )

    } else {

        # The following computes the null distribution by scoring for a random sample
        # of genes.

        null_scores <- replicate(
            n_replicate,
            scrabble::score(
                set_colnames(
                    t(
                        sc_data[
                            cells,
                            ..all_genes_filtered
                        ]
                    ),
                    cells
                ),
                list(
                    sample(
                        all_genes_filtered,
                        n_random_genes
                    )
                ),
                bin.control = TRUE,
                nbin = length(all_genes_filtered) %/% 110
            ),
            simplify = FALSE
        )

        dt <- sc_data[
            cells,
            c(
                .(
                    gene = all_genes_filtered,
                    emt_cor = cor(.SD, scores_table[id, score])[, 1]
                ),
                setNames(
                    lapply(null_scores, function(z) cor(.SD, z[id, 1])[, 1]),
                    paste0('null_', 1:length(null_scores))
                )
            ),
            by = patient,
            .SDcols = all_genes_filtered
        ][
            ,
            sapply(
                .SD,
                function(cor_type) {
                    sum(cor_type[!is.na(cor_type)])/.N
                },
                simplify = FALSE,
                USE.NAMES = TRUE
            ),
            by = gene,
            .SDcols = -'patient'
        ]

        minx <- min(dt[, -'gene'])
        maxx <- max(dt[, -'gene'])

        xvals <- seq(
            from = minx,
            to = maxx,
            length.out = nrow(dt)
        )

        pdf_funs <- list(
            emt_cor = pdf_funs$emt_cor,
            null_cor = approxfun(
                data.table(
                    x = xvals
                )[
                    ,
                    y := mean(
                        sapply(
                            pdf_funs[-1],
                            function(pfun) {
                                pfun(x)
                            }
                        )
                    ),
                    by = x
                ]
            )
        )

    }

    prob_funs <- sapply(
        pdf_funs,
        function(pfun) {

            # In the following, we divide by the total area under the curve so
            # that we get a true CDF, i.e. the maximum value it takes is 1.  I
            # think this shouldn't be necessary, because the density() function
            # should give a curve the area under which is 1, but I think we get
            # some error in the various applications of approxfun().

            cdf_fun <- approxfun(
                x = xvals,
                y = sapply(
                    xvals,
                    function(z) {
                        integrate(
                            pfun,
                            lower = minx,
                            upper = z
                        )$value
                    }
                )/integrate(
                    pfun,
                    lower = minx,
                    upper = maxx
                )$value
            )

            quantile_fun <- function(x) {
                uniroot(
                    function(y) cdf_fun(y) - x,
                    lower = minx,
                    upper = maxx
                )$root
            }

            list(
                pdf_fun = pfun,
                cdf_fun = cdf_fun,
                quantile_fun = quantile_fun
            )

        },
        simplify = FALSE,
        USE.NAMES = TRUE
    )

    # Take a quantile of the null distribution as a threshold for correlation
    # with EMT score:

    threshold <- prob_funs$null_cor$quantile_fun(null_cor_threshold)

    density_data <- data.table(
        x = xvals
    )[
        ,
        c('y_emt', 'y_null') := .(
            sapply(x, prob_funs$emt_cor$pdf_fun),
            sapply(x, prob_funs$null_cor$pdf_fun)
        )
    ]

    density_plot <- ggplot(data = density_data) +
        geom_line(
            aes(x = x, y = y_emt, colour = 'Corr. with EMT score'),
            size = 1
        ) +
        geom_line(
            aes(x = x, y = y_null, colour = 'Null distribution'),
            size = 1
        ) +
        geom_vline(
            aes(xintercept = threshold, colour = 'Threshold'),
            size = 1
            # linetype = 'dashed'
        ) +
        scale_colour_manual(
            name = '',
            values = c(
                'Corr. with EMT score' = 'red2',
                'Null distribution' = 'steelblue',
                'Threshold' = 'limegreen'
            ),
            breaks = c(
                'Corr. with EMT score',
                'Null distribution',
                'Threshold'
            )
        ) +
        labs(
            x = 'correlation',
            y = 'density',
            title = plot_title
        ) +
        theme_test()

    # The null distributions may be useful for determining when we should continue
    # with the analysis - if there is no meaningful difference between the true and
    # null distributions, then we can just put this cancer type aside.

    list(
        density_data = density_data,
        threshold = threshold,
        density_plot = density_plot,
        prob_funs = prob_funs,
        score_cor_data = dt
    )

}





rare_vs_shared_emt <- function(

    sc_data,
    genes,
    score_dist_data,
    args_list,
    scores_filter_fun_shared,
    seed = NULL,
    shared_threshold_quantile = 0.4,
    rare_threshold_quantile = 0.6

) {

    setkey(score_dist_data$score_cor_data, gene)

    set.seed(args_list$seed)

    # In the following, we're ordering the genes in the heatmaps using the correlations
    # with the EMT scores.  This doesn't look so good when the heatmaps are showing
    # expression levels, because of the variation in expression levels between genes.
    # It looks better when we convert to Z scores (i.e. scale the genes), because this
    # cancels out the effect of the expression levels.

    # We define the rare EMT genes using skewness of the distribution defined by the
    # ordered gene vectors in the Z score heatmaps.  To get an appropriate threshold
    # for the skewness, we re-calculate Z scores from data that has been shuffled cell-
    # wise, i.e. the expression levels for each cell have been shuffled.  When the cells
    # are ordered again by EMT score, all the genes have some skew, and we use quantiles
    # of the skewnesses across genes as thresholds for rare and shared EMT.

    sc_groups_list <- do.call(
        sc_groups,
        args = c(
            list(
                genes = genes,
                sc_data = sc_data[cell_type == 'cancer'],
                groups = 'cancer',
                to_keep = character(0),
                genes_order_fun = function(x) {
                    order(
                        score_dist_data$score_cor_data[
                            names(x),
                            emt_cor
                        ]
                    )
                },
                min_sig_size = 0
            ),
            args_list[names(args_list) != 'seed']
        )
    )

    # Calculate Z scores from true and shuffled data:

    z_score_data_shuffled <- copy(sc_groups_list$data)[
        ,
        (sc_groups_list$genes_filtered) := lapply(
            transpose(transpose(.SD)[, lapply(.SD, sample)]),
            function(x) {(x - mean(x))/sd(x)}
        ),
        .SDcols = sc_groups_list$genes_filtered
    ]

    # First take copy of the data from sc_groups_list, so we can change the original
    # by assignment:

    data_copy <- copy(sc_groups_list$data)

    # We've used a copy of the data for shuffling, so we can replace the original
    # data with Z scores without copying:

    sc_groups_list$data[
        ,
        (sc_groups_list$genes_filtered) := lapply(
            .SD,
            function(x) {(x - mean(x))/sd(x)}
        ),
        .SDcols = sc_groups_list$genes_filtered
    ]

    set.seed(args_list$seed)

    skewnesses <- sc_groups_list$data[
        sc_groups_list$cells_filtered
    ][
        sc_groups_list$ordering_cells$cancer,
        sort(
            sapply(
                .SD,
                function(x) {
                    e1071::skewness(
                        sample(
                            1:.N,
                            100000,
                            replace = TRUE,
                            prob = x - min(x)
                        )
                    )
                }
            ),
            decreasing = TRUE
        ),
        .SDcols = -c('id', 'cell_type', 'genes_detected', 'patient')
    ]

    thresholds <- setNames(
        z_score_data_shuffled[
            sc_groups_list$cells_filtered
        ][
            sc_groups_list$ordering_cells$cancer,
            quantile(
                sapply(
                    .SD,
                    function(x) {
                        e1071::skewness(
                            sample(
                                1:.N,
                                100000,
                                replace = TRUE,
                                prob = x - min(x)
                            )
                        )
                    }
                ),
                c(shared_threshold_quantile, rare_threshold_quantile)
            ),
            .SDcols = -c('id', 'cell_type', 'genes_detected', 'patient')
        ],
        c('shared', 'rare')
    )

    rare_emt_genes <- score_dist_data$score_cor_data[
        sc_groups_list$genes_filtered
    ][
        emt_cor >= score_dist_data$prob_funs$emt_cor$quantile_fun(0.95) &
            skewnesses[gene] >= thresholds['rare']
    ][
        order(-emt_cor),
        gene
    ]

    shared_emt_genes <- data_copy[
        sc_groups_list$cells_filtered,
        sc_groups_list$genes_filtered[
            apply(.SD, 2, scores_filter_fun_shared) &
                skewnesses[sc_groups_list$genes_filtered] <= thresholds['shared']
        ],
        .SDcols = sc_groups_list$genes_filtered
    ]

    shared_emt_genes <- score_dist_data$score_cor_data[
        shared_emt_genes
    ][
        order(-emt_cor),
        gene
    ]

    c(
        list(sc_groups_list = sc_groups_list),
        list(score_dist_data = score_dist_data),
        list(z_score_data_shuffled = z_score_data_shuffled),
        list(
            rare_shared_emt_data = list(
                skewnesses = skewnesses,
                thresholds = thresholds,
                rare_emt_genes = rare_emt_genes,
                shared_emt_genes = shared_emt_genes
            )
        )
    )

}





simulate_counts <- function(

    types,
    initial_types = NULL,
    proportions_vec = seq(0.1, 0.8, by = 0.1),
    n = 100,
    max_mean_count = 1000,
    sd_frac = 0.2,
    density_fun = runif,
    density_fun_args = sapply(
        types,
        function(x) { # This x isn't actually used.
            list(
                min = 0,
                max = max_mean_count
            )
        },
        simplify = FALSE,
        USE.NAMES = TRUE
    )

) {

    # This function simulates counts for a collection of types using a specified density function
    # (uniform density by default) and arguments.  The types argument can either be a number,
    # indicating the number of types for which you want to simulate counts, or a vector of type
    # names (e.g. names of cell types), in which case the simulated counts will be annotated with
    # these names.

    # Additionally, the user can specify "initial" types (again, either a number or vector of
    # type names) and a vector of proportions.  In this case, the total for the initial types will
    # be sampled first using a normal distribution (with mean and standard deviation derived from
    # the arguments max_mean_count and sd_frac), and the remaining cell types will be sampled
    # using density_fun and then scaled so that the combined proportion of the initial types is
    # approximately as specified in proportions_vec (it is not exact, due to integer rounding).
    # The counts for the individual elements of initial_types are again sampled using density_fun
    # and scaled so that their total equals the count sampled via normal distribution.

    density_fun <- match.fun(density_fun)

    if(is.numeric(types) & length(types) == 1) {
        types <- 1:types
    }

    if(!is.null(initial_types)) {

        if(is.numeric(initial_types) & length(initial_types) == 1) {
            initial_types <- 1:initial_types
        }

        if(sum(initial_types %in% types) < length(initial_types)) {
            stop('Not all of initial_types are in types.')
        }

    }

    counts_table <- as.data.table( # Generate counts for each type using density_fun:
        sapply(
            types,
            function(type) {
                round(
                    do.call(
                        density_fun,
                        c(
                            density_fun_args[[type]],
                            list(n = n*length(proportions_vec))
                        )
                    )
                )
            },
            simplify = FALSE,
            USE.NAMES = TRUE
        )
    )

    if(is.null(initial_types)) {

        return(counts_table)

    } else {

        other_types <- types[!(types %in% initial_types)]

        # Generate combined count for initial_types by sampling from the set of integers
        # 1:(2*max_mean_count) using a normal distribution with mean and sd derived from
        # max_mean_count and sd_frac, scaled using the supplied proportions:

        counts_table <- cbind(
            rbindlist(
                lapply(
                    proportions_vec,
                    function(p) {
                        sim_counts <- data.table(
                            proportion = p,
                            initial_types_count = sample(
                                1:(2*max_mean_count),
                                size = n,
                                replace = TRUE,
                                prob = dnorm(
                                    1:(2*max_mean_count),
                                    mean = max_mean_count*p/max(proportions_vec),
                                    sd = sd_frac*max_mean_count*p/max(proportions_vec)
                                )
                            )
                        )
                    }
                )
            ),
            counts_table[ # Change zeros in initial_types columns to 1s, to avoid NaNs
                , # during scaling:
                (initial_types) := lapply(
                    .SD,
                    function(x) {
                        x[x == 0] <- 1
                        x
                    }
                ),
                .SDcols = initial_types
                ]
        )[ # Scale counts for types not in initial_types to make proportion of
            , # initial_types_count equal to the specified proportion:
            (other_types) := as.list(
                round(
                    as.numeric(.SD)*initial_types_count*(1 - proportion)/
                        (proportion*sum(as.numeric(.SD)))
                )
            ),
            by = row.names(counts_table), # .I and row.names(.SD) don't work here
            .SDcols = other_types
            ][ # Scale counts for initial types so that they sum to initial_types_count:
                ,
                (initial_types) := as.list(
                    round(
                        as.numeric(.SD)*initial_types_count/
                            sum(as.numeric(.SD))
                    )
                ),
                by = row.names(counts_table),
                .SDcols = initial_types
                ][ # Calculate the actual median proportion achieved for each specified proportion:
                    ,
                    median_proportion := median(
                        initial_types_count/(initial_types_count + rowSums(.SD))
                    ),
                    by = proportion,
                    .SDcols = other_types
                    ]

        return(list(counts_table = counts_table, initial_types = initial_types))

    }

}





sample_indices <- function(types_vec, counts_table) {

    apply(
        counts_table,
        1,
        function(s) {
            sapply(
                names(s),
                function(type) {

                    N <- sum(types_vec == type)
                    inds <- which(types_vec == type)

                    # If we're trying to sample more indices than actually exist for a given
                    # type, take all indices for that type the appropriate number of times,
                    # then sample for the remaining number.  This means some indices will be
                    # chosen multiple times.

                    # switch(
                        # (s[[type]] > N) + 1,
                        # sample(inds, s[[type]]),
                        # c(rep(inds, s[[type]] %/% N), sample(inds, s[[type]] %% N))
                    # )

					# I changed the above to just sample with replacement, which is conceptually simpler and arguably more appropraite:
					sample(inds, s[[type]], replace = TRUE)

                },
                simplify = FALSE,
                USE.NAMES = TRUE
            )
        }
    )

}





type_contrib <- function(

    data,
    cols,
    counts_table,
    sampled_indices,
    types_var = 'cell_type',
    initial_types = NULL,
    normalise_fun = function(x) {x/quantile(x[x > 0], 0.95)}

) {

    # This function is much faster if the data contains only the relevant columns.  If it
    # contains many more, we'll take just the relevant subset, to save computation time.  This
    # involves taking a copy, so it would save memory if we supplied to the function only the
    # columns that we need, namely those defined by genes and types_var.

    if(ncol(data) > length(cols) + 100) {
        data <- data[, c(..types_var, ..cols)]
        warning(
            "It looks like the data contains a lot of unnecessary columns. ",
            "To save computation time, we're copying the relevant subset of the data. ",
            "It would be more efficient to subset the data beforehand."
        )
    }

    # Optionally normalise the columns of interest using normalise_fun:

    if(!is.null(normalise_fun)) {

        normalise_fun <- match.fun(normalise_fun)

        # Start by taking a copy of (the relevant subset of) the data, so that we don't change the columns of the original data table.
        data <- data[, c(..types_var, ..cols)]
        data[, (cols) := lapply(.SD, normalise_fun), .SDcols = cols]

    }

    if(is.null(initial_types)) {

        # If initial_types is not specified, we calculate the proportions and contributions of
        # all the types.  For this option, counts_table should have one column per type and no
        # other columns.

        types <- unique(data[[types_var]])

        data.table(
            simulation = unlist(lapply(1:nrow(counts_table), rep, length(types))),
            type = rep(types, nrow(counts_table)),
            proportion_content = unlist( # Calculate proportion for each type:
                as.data.table(t(counts_table))[, lapply(.SD, function(s) {s/sum(s)})]
            ),
            proportion_contrib = unlist(
                lapply(
                    sampled_indices,
                    function(s) {
                        data[
                            unlist(s),
                            .(totals = sum(rowSums(.SD))), # Can't we just use sum() here?
                            .SDcols = cols,
                            keyby = types_var
                        ][
							types,
							plyr::mapvalues(totals, NA, 0, warn_missing = FALSE)/sum(totals[!is.na(totals)])
                        ]
                    }
                )
            )
        )

    } else {

        # If initial_types is specified, we'll calculate only the combined proportion of contribution of
        # the initial types, and we'll take for the content proportions the median_proportion column of
        # counts_table.  Note that taking sum over a data table sums all the elements, like summing a
        # matrix.  We don't need to do e.g. sum(rowSums()).

        data.table(
            simulation = 1:nrow(counts_table), # nrow(counts_table) should equal length(sampled_indices)
            proportion_content = counts_table$median_proportion,
            proportion_contrib = unlist(
                lapply(
                    sampled_indices,
                    function(s) data[unlist(s, use.names = FALSE), sum(.SD[get(types_var) %in% initial_types])/sum(.SD), .SDcols = cols]
                )
            )
        )

    }

}





# In the following, the purpose of the legend_labels and legend_colours arguments are
# to manually set the colours that correspond to particular cell types.  I have put as
# the default a list of the cell types I think I will need for my lineplots, with
# colours derived from RColorBrewer's Set3 palette.  Some are the same - e.g. myeloid
# and macrophage are both blue - to make it easier to compare plots with slightly
# differing but comparable categories.

simulated_tumours_lineplot <- function(

    single_cell_data,
    genes,
    types_var = 'cell_type',
	initial_types = NULL,
	normalise_fun = function(x) {x/quantile(x[x > 0], 0.95)},
    add_zero_one = 0,
    x_axis_title = 'Proportion of tumour',
    y_axis_title = 'Proportion of gene expression',
    legend_title = 'Cell type',
    legend_labels = c(
        'b_cell' = 'B cell',
		'caf' = 'CAF',
        'cancer' = 'Cancer',
        'endothelial' = 'Endothelial',
        'fibroblast' = 'Fibroblast',
        'macrophage' = 'Macrophage',
        'mast' = 'Mast',
        'myeloid' = 'Myeloid',
        'plasma' = 'Plasma',
        't_cell' = 'T cell',
        'other' = 'Other',
        'lymphoid' = 'T & NK cells'
    ),
    legend_colours = c(
        'b_cell' = '#8DD3C7',
		'caf' = '#FDB462',
        'cancer' = '#FB8072',
        'endothelial' = '#BC80BD',
        'fibroblast' = '#FDB462',
        'macrophage' = '#80B1D3',
        'mast' = '#FFED6F',
        'myeloid' = '#80B1D3',
        'plasma' = '#FCCDE5',
        't_cell' = '#B3DE69',
        'other' = 'grey',
        'lymphoid' = '#B3DE69'
    ),
    legend_key_size = NULL,
    plot_title = '',
    point_size = NULL,
    point_stroke = NULL,
    line_size = NULL,
    error_bars = FALSE,
    error_bars_width = 0.01,
    error_bars_size = 0.25, # I think 0.5 is the ggplot2 default; I guess the same is true for grid_line_size.
    grid_line_size = 0.25,
    panel_border_size = NULL,
    ...

) {

    # This is a wrapper function that runs the sampling and aggreating functions above for
    # several cell types, then makes a lineplot.

	# Note that the order of <initial_types> (if given) determines the order in which the corresponding lines and points are plotted,
	# the last one being plotted last and therefore appearing on top.  This is done by supplying ggplot() with a factor of initial
	# types with levels in the order that they appear in <initial_types>.

    # The <add_zero_one> argument is used to specify that you want to add points (0, 0) and/
    # or (1, 1) to the plot.  Supply 0 (the default) if you want the (0, 0) point; 1 if
    # you want (1, 1); and c(0, 1) if you want both.

    # The ... is for additional arguments to the simulate_counts() function.  I anticipate
    # it being used mostly to adjust max_mean_count to suit the supplied dataset.

    cell_types <- unique(single_cell_data[[types_var]])
	if(is.null(initial_types)) {initial_types <- cell_types}
	if(!is.null(normalise_fun)) {normalise_fun <- match.fun(normalise_fun)}

    proportions_table <- rbindlist(
        lapply(
			initial_types,
            # cell_types,
            function(ct) {

                simulated_counts <- simulate_counts(cell_types, initial_types = ct, ...)

                sampled_indices <- sample_indices(single_cell_data[[types_var]], simulated_counts$counts_table[, ..cell_types])

                proportions_table <- type_contrib(
                    single_cell_data[, c(..types_var, ..genes)],
                    genes,
                    simulated_counts$counts_table,
                    sampled_indices,
                    initial_types = ct,
					normalise_fun = normalise_fun
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

    if(!is.null(add_zero_one)) {
        plot_data <- rbind(
            plot_data,
            data.table(
                initial_type = rep(unique(plot_data$initial_type), length(add_zero_one)),
                proportion_content = unlist(lapply(add_zero_one, rep, length(unique(plot_data$initial_type)))),
                mean_proportion_contrib = unlist(lapply(add_zero_one, rep, length(unique(plot_data$initial_type)))),
                sd_proportion_contrib = unlist(lapply(add_zero_one, rep, length(unique(plot_data$initial_type))))
            )
        )#[order(initial_type, proportion_content)]
    }

	plot_data[, initial_type := factor(initial_type, levels = initial_types)]
	plot_data <- plot_data[order(initial_type, proportion_content)]

    lineplot <- ggplot(
        data = plot_data,
        aes(x = proportion_content, y = mean_proportion_contrib, colour = initial_type) # colour = factor(initial_type, levels = initial_types)
    )

    if(is.null(line_size)) {
        lineplot <- lineplot + geom_line()
    } else {
        lineplot <- lineplot + geom_line(size = line_size)
    }

    if(is.null(point_size)) {
        if(is.null(point_stroke)) {
            lineplot <- lineplot + geom_point(shape = 0)
        } else {
            lineplot <- lineplot + geom_point(shape = 0, stroke = point_stroke)
        }
    } else {
        if(is.null(point_stroke)) {
            lineplot <- lineplot + geom_point(shape = 0, size = point_size)
        } else {
            lineplot <- lineplot + geom_point(shape = 0, size = point_size, stroke = point_stroke)
        }
    }

    if(error_bars) {
        lineplot <- lineplot + geom_errorbar(
            aes(ymin = mean_proportion_contrib - sd_proportion_contrib, ymax = mean_proportion_contrib + sd_proportion_contrib),
            width = error_bars_width,
            size = error_bars_size
        )
    }

    lineplot <- lineplot +
        # geom_line(size = line_size) +
        # geom_point(shape = 0, size = point_size, stroke = point_stroke) +
        theme(
            panel.background = element_blank(),
            panel.border = element_rect(fill = NA, colour = 'black', size = panel_border_size),
            panel.grid.minor = element_blank(),
            panel.grid.major.x = element_blank(),
            panel.grid.major.y = element_line(colour = 'grey', size = grid_line_size, linetype = 'dotted')
        ) +
        scale_y_continuous(
            limits = c(0, 1),
            breaks = c(0, 0.25, 0.5, 0.75, 1),
            labels = c('0' = '0', '0.25' = '0.25', '0.5' = '0.5', '0.75' = '0.75', '1' = '1'),
            expand = c(0.02, 0.02)
        ) +
        scale_x_continuous(expand = c(0.02, 0.02)) +
        scale_colour_manual(
            labels = switch(is.null(legend_labels) + 1, legend_labels, sort(unique(initial_type))),
            values = switch(
                is.null(legend_colours) + 1,
                legend_colours,
                setNames(
                    hcl(h = seq(360/(2*length(cell_types)), 360 - 360/(2*length(cell_types)), length.out = length(cell_types)), c = 100, l = 65),
                    sort(unique(initial_type))
                )
            )
        ) +
        labs(x = x_axis_title, y = y_axis_title, colour = legend_title, title = plot_title)

    if(!is.null(legend_key_size)) {
        lineplot <- lineplot + theme(legend.key.size = legend_key_size)
    }

    list(data = plot_data, lineplot = lineplot)

}





simulated_tumours_data <- function(

    sc_mat,
    types,
    initial_types = 'cancer',
    # genes = names(sc_data[, -c('id', 'patient', 'cell_type')]),
    proportions_vec = seq(0.1, 0.9, length.out = 40),
    n = 25,
    id_prefix = '',
    show_progress = TRUE,
    restore_log = TRUE,
    ...

) {

    # The ... is for additional arguments to the simulate_counts() function.

    cell_types <- unique(types)

    simulated_counts <- simulate_counts(
        cell_types,
        initial_types = initial_types,
        proportions_vec = proportions_vec,
        n = n,
        ...
    )

    sampled_indices <- sample_indices(
        types,
        simulated_counts$counts_table[, ..cell_types]
    )

    # The lapply() application in which we construct expr is much faster if we prepare the data
    # beforehand, i.e. subset genes and reverse the log.  Originally I simply used the following
    # very simple code to achieve this:

    # data_for_aggregating <- 2^sc_data[, ..genes] - 1

    # This is quite fast, but it requires taking a copy and thus uses a lot of memory - too much
    # in the case of Lambrechts's lung scRNA-seq data, when R runs out of memory and throws an
    # error.  The same is true for the following, in which lapply() presumably has to remember
    # everything it produces until the last iteration, when it returns it all at once:

    # sc_data[
    #     ,
    #     (genes) := lapply(
    #         .SD,
    #         function(x) {
    #             round(2^as.numeric(x) - 1, 4)
    #         }
    #     ),
    #     .SDcols = genes
    # ]

    # So we have no choice but to use a for loop.  I settled on the following, which uses set()
    # and thus should be much faster than a for loop which uses `:=` in every iteration.  It's
    # still pretty slow, so it would be ideal if we could detect how much memory we have
    # available, calculate how much we would need for 2^sc_data[, ..genes] - 1, and
    # run the for loop only when there isn't enough available memory.

    # for(g in genes) {
    #     set(
    #         sc_data,
    #         j = g,
    #         value = round(2^sc_data[[g]] - 1, 4)
    #     )
    # }

    sc_mat <- round(2^sc_mat - 1, 4)

    if(show_progress) {

        cat(
            length(sampled_indices),
            'profiles to build.  Currently on iteration:\n'
        )

        expr <- t( # Transpose because genes end up as rows after sapply()
            sapply(
                1:length(sampled_indices),
                function(i) {

                    cat(
                        rep('\b', ceiling(log10(i))),
                        i,
                        sep = ''
                    )

                    colSums(sc_mat[unlist(sampled_indices[[i]]), ])

                }
            )
        )

        # expr <- log2(
        #     rbindlist(
        #         lapply(
        #             1:length(sampled_indices),
        #             function(i) {
        #
        #                 cat(
        #                     rep('\b', ceiling(log10(i))),
        #                     i,
        #                     sep = ''
        #                 )
        #
        #                 as.list(
        #                     colSums(
        #                         sc_mat[unlist(sampled_indices[[i]]), ..genes]
        #                     )
        #                 )
        #
        #             }
        #         )
        #     ) + 1
        # )

        cat('\nDone!')

    } else {

        expr <- t( # Transpose because genes end up as rows after sapply()
            sapply(
                1:length(sampled_indices),
                function(i) colSums(sc_mat[unlist(sampled_indices[[i]]), ])
            )
        )

        # expr <- log2(
        #     rbindlist(
        #         lapply(
        #             1:length(sampled_indices),
        #             function(i) {
        #                 as.list(
        #                     colSums(
        #                         sc_mat[unlist(sampled_indices[[i]])]
        #                     )
        #                 )
        #             }
        #         )
        #     ) + 1
        # )

    }

    rownames(expr) <- paste0(id_prefix, 1:length(sampled_indices))

    # expr[
    #     ,
    #     id := paste0(id_prefix, .I)
    # ]

    # setcolorder(expr, 'id')

    meta <- simulated_counts$counts_table[
        ,
        .(purity = cancer/sum(as.numeric(.SD))),
        by = .(id = paste0(id_prefix, row.names(simulated_counts$counts_table))),
        .SDcols = cell_types
    ]

    # if(restore_log) {
    #     for(g in genes) {
    #         set(
    #             sc_data,
    #             j = g,
    #             value = round(log2(sc_data[[g]] + 1), 4)
    #         )
    #     }
    # }

    if(restore_log) {
        expr <- round(log2(expr + 1), 4)
    }

    list(
        expression_data = expr,
        meta_data = meta
    )

}





gene_stat <- function(

    sc_groups_list,
    x_fun = mean,
    y_fun = var,
    mod_fun = quantreg::rq,
    mod_fun_args = list(
        formula = y ~ x - 1,
        tau = 0.75
    ),
    cond_fun = NULL

) {

    x_fun <- match.fun(x_fun)
    y_fun <- match.fun(y_fun)
    mod_fun <- match.fun(mod_fun)

    setkey(sc_groups_list$data, id)

    stat_data <- sc_groups_list$data[
        sc_groups_list$cells_filtered
        ][
            cell_type == 'cancer',
            .(
                gene = names(.SD),
                x = sapply(.SD, x_fun),
                y = sapply(.SD, y_fun)
            ),
            .SDcols = sc_groups_list$genes_filtered
            ]

    mod <- do.call(
        mod_fun,
        args = c(
            mod_fun_args,
            switch(
                is.null(cond_fun) + 1,
                list(data = stat_data[cond_fun(x, y)]),
                list(data = stat_data)
            )
        )
    )

    list(
        data = stat_data,
        mod = mod
    )

}





rare_genes <- function(

    gene_stat_data,
    y_mod_cond_fun = `>`,
    x_cond_fun = function(v) {v > quantile(v, 0.25)},
    y_cond_fun = x_cond_fun,
    title = NULL,
    x_lab = NULL,
    y_lab = NULL,
    rare_genes_colour = 'darkorange',
    mod_line_colour = 'royalblue3'

) {

    y_mod_cond_fun <- match.fun(y_mod_cond_fun)
    x_cond_fun <- match.fun(x_cond_fun)
    y_cond_fun <- match.fun(y_cond_fun)

    rare_genes_data <- gene_stat_data$data[
        y_mod_cond_fun(y, predict(gene_stat_data$mod, newdata = gene_stat_data$data))
        ][
            x_cond_fun(x) & y_cond_fun(y)
            ]

    rare_genes_plot <- qplot(
        x,
        y,
        data = gene_stat_data$data,
        main = title,
        xlab = x_lab,
        ylab = y_lab
    ) +
        geom_point(
            aes(x = x, y = y),
            data = rare_genes_data,
            colour = rare_genes_colour
        ) +
        geom_line(
            aes(
                x = x,
                y = predict(
                    gene_stat_data$mod,
                    newdata = gene_stat_data$data
                )
            ),
            data = gene_stat_data$data,
            size = 1,
            colour = mod_line_colour
        ) +
        theme_test()

    list(
        genes = rare_genes_data$gene,
        plot = rare_genes_plot
    )

}
