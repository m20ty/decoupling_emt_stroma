head2 <- function(dt, m = 6L, n = 6L) {
    dt[1:m, 1:n]
}

tail2 <- function(dt, m = 6L, n = 6L) {
    dt[(.N - m + 1):.N, 1:n]
}

tdt <- function(dt, new_id = 'id') {
    
    #Transpose a data table, keeping the entries of the first column for new column names
    #and choosing a name for the new ID variable.
    
    out <- cbind(names(dt)[-1], transpose(dt[, -1]))
    names(out) <- c(new_id, as.character(dt[[1]]))
    out
    
}

scaledt <- function(
    
    dt,
    centre = TRUE,
    scale = TRUE,
    margin = c(1, 2),
    skip_cols = 'id'
    
) {
    
    margin <- match.arg(as.character(margin), c('1', '2'))
    
    if(margin == 1) {
        
        cbind(
            dt[, skip_cols, with = FALSE],
            as.data.table(t(scale(
                t(dt[, -skip_cols, with = FALSE]),
                center = centre,
                scale = scale
            )))
        )
        
    } else {
        
        cbind(
            dt[, skip_cols, with = FALSE],
            as.data.table(scale(
                dt[, -skip_cols, with = FALSE],
                center = centre,
                scale = scale
            ))
        )
        
    }
    
}

blank_plot <- function() {
    ggplot() +
        theme(
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.ticks.length = unit(0, 'pt'),
            axis.title = element_blank(),
            plot.background = element_blank(),
            panel.background = element_blank(),
            plot.margin = unit(c(0, 0, 0, 0), 'pt')
        )
}

blank_dummy_plot <- function() {
    
    ggplot(
        data = data.frame(a = 0, b = 0),
        aes(x = a, y = b)
    ) +
        geom_point(colour = 'white') +
        theme(
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank(),
            panel.background = element_rect(fill = 'white'),
            plot.margin = unit(c(0, 0, 0, 0), 'pt')
        )
    
}





#Converting to HGNC approved symbols:

symbol_to_hgnc <- function(gene_symbol, gene_names_df) {
    if(gene_symbol %in% gene_names_df$approved_symbol) {
        hgnc_translation <- gene_symbol
    } else {
        crit <- gene_names_df$alias_symbol == gene_symbol | 
            gene_names_df$previous_symbol == gene_symbol
        hgnc_translation <- unique(gene_names_df$approved_symbol[crit])
        hgnc_translation <- unique(hgnc_translation[!is.na(hgnc_translation)])
    }
    if(length(hgnc_translation) == 1) {
        hgnc_translation
    } else if(length(hgnc_translation) > 1) {
        warning('This symbol has multiple HGNC forms.')
        gene_symbol
    } else {
        warning('This has no HGNC form according to gene_names_df. Are you sure the symbol
            is correct?')
        gene_symbol
    }
}





#The following function uses the above symbol_to_hgnc() function:

symbol_to_hgnc_list <- function(symbol_list, gene_names_df) {
    new_symbols <- symbol_list
    for(i in 1:length(new_symbols)) {
        hgnc_translation <- symbol_to_hgnc(new_symbols[i], gene_names_df)
        if(hgnc_translation != new_symbols[i] & !(hgnc_translation %in% new_symbols[-i])) {
            new_symbols[i] <- hgnc_translation
        }
    }
    new_symbols
}





# The following does the same as the above two functions together, except that it calculates
# the average correlation with the cell type markers, rather than the correlation with the
# average expression of the cell type markers.

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





heat_map <- function(
    
    mat,
    ordering,
    axis_title_y = '',
    axis_title_x = axis_title_y,
    axis_text_y = rownames(mat),
    axis_text_x = colnames(mat),
    legend_title = 'correlation',
    colours = rev(colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))(50)),
    colour_limits = c(-1, 1),
    oob = scales::squish,
    legend_breaks = waiver(),
    axis_text_size = 7,
    plot_margin = c(5.5, 5.5, 5.5, 5.5),
    border = FALSE,
    plot_title = '',
    ...
    
) {
    
    #plot_margin should be a numeric vector specifying the margin around the entire
    #plot, with sizes of top, right, bottom and left margins, *in that order*.  This
    #is the order used in the plot.margin argument of ggplot2's theme() function.
    
    #The '...' parameter is for extra arguments to theme().
    
    #Note this function outputs a ggplot2 object, which you can change as with any other
    #such object, e.g. g <- g + ...
    
    g <- ggplot(
        melt(mat[ordering, ordering]),
        aes(x = Var1, y = Var2)
        
    ) +
        geom_raster(aes(fill = value)) +
        scale_fill_gradientn(
            colours = colours,
            limits = colour_limits, #With the default colours, this makes zero appear as white.
            oob = oob,
            breaks = legend_breaks
        ) +
        theme(
            axis.ticks = element_blank(),
            axis.ticks.length = unit(0, 'pt'),
            axis.text = element_text(size = axis_text_size),
            axis.text.y = element_text(vjust = 0.5), #vjust makes labels in centre of rows
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), #likewise columns
            plot.margin = unit(plot_margin, 'pt'),
            ...
        ) +
        scale_y_discrete(
            breaks = axis_text_y,
            expand = c(0, 0)
        ) +
        scale_x_discrete(
            breaks = axis_text_x,
            expand = c(0, 0)
        ) +
        labs(
            fill = 'correlation',
            x = axis_title_x,
            y = axis_title_y,
            title = plot_title
        )
    
    if(paste(axis_text_y, collapse = '') == '') {
        g <- g + theme(axis.text.y = element_blank())
    }
    
    if(paste(axis_text_x, collapse = '') == '') {
        g <- g + theme(axis.text.x = element_blank())
    }
    
    if(axis_title_y == '') {
        g <- g + theme(axis.title.y = element_blank())
    }
    
    if(axis_title_x == '') {
        g <- g + theme(axis.title.x = element_blank())
    }
    
    if(plot_title == '') {
        g <- g + theme(plot.title = element_blank())
    }
    
    #if(!is.null(axis_labs)) {
    #  axis_text <- colnames(mat)
    #  axis_text[!(axis_text %in% axis_labs)] <- ''
    #  g <- g + scale_y_discrete(breaks = axis_text)
    #}
    
    g
    
}





heat_map_labels_repel <- function(
    
    labels_vec,
    edge = c('top', 'right', 'bottom', 'left'),
    text_size = 3,
    nudge = 0.35,
    axis_title = NULL,
    axis_title_size = 11, # Use 14 if you want it to look like a plot title
    title_edge_margin = 0, # In pts
    ...
    
) {
    
    #The '...' parameter is for additional arguments to ggrepel::geom_text_repel().
    
    # nudge = 0.25 pushes the labels a quarter of the way up/down/across the plot area.  This is roughly
    # appropriate because I'm using vjust to align the labels along one edge, so a certain space will be
    # occupied by the segments - I think a quarter of the space for the segments should be about right.
    
    edge <- match.arg(edge, c('top', 'right', 'bottom', 'left'))
    
    g <- ggplot(
        data.frame(
            pos = 1:length(labels_vec),
            label = labels_vec
        )
    ) + theme(
        panel.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks.length = unit(0, 'pt'),
        axis.ticks = element_blank(),
        axis.title = element_text(size = axis_title_size)
    )
    
    if(edge == 'top') {
        
        g <- g +
            theme(
                plot.margin = unit(c(0, 0, title_edge_margin, 0), 'pt')
            ) +
            ggrepel::geom_text_repel(
                aes(x = pos, y = 0, label = label),
                size = text_size,
                angle = 90,
                # nudge_y pushes the labels in the y direction up or down by a proportion of the plot area height.
                # We're using the negative of nudge here, because we want the labels to be pushed down from the top.
                nudge_y = -nudge,
                direction = 'x',
                force = 0.25,
                vjust = 1,
                ...
            ) +
            scale_x_continuous(
                limits = c(1, length(labels_vec)),
                #First limit must be 1 or else first label will be 1 from left edge
                expand = c(0, 0.5),
                #Expand 2nd component must be 0.5 to push lines into the middle of the columns
                breaks = NULL,
                labels = NULL,
                name = axis_title
            ) +
            scale_y_continuous(
                limits = c(-1, 0),
                expand = c(0, 0),
                breaks = NULL,
                labels = NULL,
                name = NULL
            )
        
    }
    
    if(edge == 'right') {
        
        g <- g +
            theme(
                plot.margin = unit(c(0, 0, 0, title_edge_margin), 'pt')
            ) +
            ggrepel::geom_text_repel(
                aes(x = 0, y = pos, label = label),
                size = text_size,
                nudge_x = -nudge,
                direction = 'y',
                force = 0.25,
                hjust = 1,
                ...
            ) +
            scale_x_continuous(
                limits = c(-1, 0),
                expand = c(0, 0),
                breaks = NULL,
                labels = NULL,
                name = NULL
            ) +
            scale_y_continuous(
                limits = c(1, length(labels_vec)),
                expand = c(0, 0.5),
                breaks = NULL,
                labels = NULL,
                name = axis_title
            )
        
    }
    
    if(edge == 'bottom') {
        
        # To get an axis title, we need to duplicate the x axis to get a second one at the top,
        # then set everything to blank on the bottom so that the label lines still reach the
        # edge of the plot area.
        
        g <- g +
            theme(
                axis.title.x.bottom = element_blank(),
                plot.margin = unit(c(title_edge_margin, 0, 0, 0), 'pt')
            ) +
            ggrepel::geom_text_repel(
                aes(x = pos, y = 0, label = label),
                size = text_size,
                angle = 90,
                nudge_y = nudge,
                direction = 'x',
                force = 0.25,
                vjust = 0,
                ...
            ) +
            scale_x_continuous(
                limits = c(1, length(labels_vec)),
                expand = c(0, 0.5),
                breaks = 1:length(labels_vec),
                labels = rep('', length(labels_vec)),
                name = axis_title,
                position = 'top',
                sec.axis = dup_axis()
                
            ) +
            scale_y_continuous(
                limits = c(0, 1),
                expand = c(0, 0),
                breaks = NULL,
                labels = NULL,
                name = NULL
            )
        
    }
    
    if(edge == 'left') {
        
        # As in the 'bottom' case, we duplicate the y axis so that an axis title can appear.
        
        g <- g +
            theme(
                axis.title.y.left = element_blank(),
                plot.margin = unit(c(0, title_edge_margin, 0, 0), 'pt')
            ) +
            ggrepel::geom_text_repel(
                aes(x = 0, y = pos, label = label),
                size = text_size,
                nudge_x = nudge,
                direction = 'y',
                force = 0.25,
                hjust = 0,
                ...
            ) +
            scale_x_continuous(
                limits = c(0, 1),
                expand = c(0, 0),
                breaks = NULL,
                labels = NULL,
                name = NULL
            ) +
            scale_y_continuous(
                limits = c(1, length(labels_vec)),
                expand = c(0, 0.5),
                breaks = 1:length(labels_vec),
                labels = rep('', length(labels_vec)),
                name = axis_title,
                position = 'right',
                sec.axis = dup_axis()
            )
        
    }
    
    g
    
}





heat_map_bar <- function(
    
    named_vec,
    ordering,
    colours,
    colour_limits = c(NA, NA),
    oob = scales::squish,
    fun = function(x) x,
    axis_title_x = '',
    axis_title_y = '',
    legend_title = '',
    legend_breaks = waiver(),
    plot_margin = c(5.5, 5.5, 5.5, 5.5),
    ...
    
) {
    
    #The '...' is for extra arguments to theme(), allowing you to customise axis titles etc.
    
    ggplot(
        data.frame(
            x = factor(
                names(named_vec)[ordering],
                levels = names(named_vec)[ordering]
            ),
            value = as.numeric(fun(named_vec[ordering]))
        ),
        aes(x = x, y = 0)
    ) + 
        geom_raster(
            aes(fill = value)
        ) + 
        scale_fill_gradientn(
            #colours = rev(colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))(50)),
            colours = colours,
            limits = colour_limits,
            oob = oob,
            breaks = legend_breaks
        ) +
        theme(
            axis.ticks = element_blank(),
            axis.ticks.length = unit(0, 'pt'),
            axis.text = element_blank(),
            plot.margin = unit(plot_margin, 'pt'),
            ...
        ) +
        scale_y_discrete(
            #breaks = axis_text_y,
            expand = c(0, 0)
        ) +
        scale_x_discrete(
            #breaks = axis_text_x,
            expand = c(0, 0)
        ) +
        labs(
            x = axis_title_x,
            y = axis_title_y,
            fill = legend_title
        )
    
}





col_side_colours <- function(
    
    #Add a way to remove the key, and maybe customise the key title.
    
    category_vec,
    ordering = NULL,
    plot_margin = NULL,
    key = TRUE
    
) {
    
    if(!is.null(ordering)) {
        
        g <- ggplot(
            data.frame(
                name = factor(
                    as.character(1:length(category_vec))[ordering],
                    levels = as.character(1:length(category_vec))[ordering]
                ),
                category = category_vec[ordering]
            ),
            aes(
                x = name,
                y = 'a',
                color = category,
                fill = category
            )
        )
        
    } else {
        
        g <- ggplot(
            data.frame(
                name = factor(
                    as.character(1:length(category_vec)),
                    levels = as.character(1:length(category_vec))
                ),
                category = category_vec
            ),
            aes(
                x = name,
                y = 'a',
                color = category,
                fill = category
            )
        )
        
    }
    
    g <- g +
        geom_raster(show.legend = key) +
        
        theme(
            axis.text = element_blank(),
            axis.ticks = element_blank(),
            axis.title = element_blank(),
            panel.background = element_blank()
        ) +
        
        scale_y_discrete(
            expand = c(0, 0)
        ) +
        
        scale_x_discrete(
            expand = c(0, 0)
        )
    
    #if(key) {
    #  g <- g + labs(fill = 'category')
    #} else {
    #  g <- g + labs(fill = NULL)
    #}
    
    if(!is.null(plot_margin)) {
        #Optionally tweak margins.  Particularly useful for combined plots.
        g <- g + theme(plot.margin = unit(plot_margin, 'pt'))
    }
    
    g
    
}





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





infer_subtypes <- function(
    
    expression_data,
    subtypes_dt
    
) {
    
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





top_cols_by_fun_cor <- function(
    
    x,
    initial = c(
        'SNAI1',
        'SNAI2',
        'TWIST1',
        'VIM',
        'ZEB1',
        'ZEB2'
    ),
    FUN = mean,
    threshold = 0.6
    
) {
    
    # x should be a matrix or data frame/table, among whose columns we will look for high
    # correlations (so the elements of initial should also be columns).  If it is a data
    # frame/table, it's important that there are no ID or non-numeric columns.
    
    FUN <- match.fun(FUN)
    
    # Convert x to data.table, in case it isn't already:
    
    x <- as.data.table(x)
    
    # Remove columns with zero standard deviation:
    
    x[
        ,
        names(x) := lapply(
            .SD,
            function(x) {
                switch((sd(x) > 0) + 1, NULL, x)
            }
        )
    ]
    
    # Calculate correlations with elements of initial, and apply FUN:
    
    fun_cor_with_initial <- x[
        ,
        apply(
            cor(
                .SD[, -..initial],
                .SD[, ..initial]
            ),
            1,
            FUN
        )
    ]
    
    # Output elements which surpass threshold, along with their correlation values:
    
    rbind(
        data.table(id = initial, correlation = 1),
        data.table(
            id = names(fun_cor_with_initial),
            correlation = fun_cor_with_initial
        )[
            correlation > threshold
        ][
            order(-correlation)
        ]
    )
    
}





score_top_cols_by_cor <- function(
    
    x,
    initial = c(
        'SNAI1',
        'SNAI2',
        'TWIST1',
        'VIM',
        'ZEB1',
        'ZEB2'
    ),
    threshold = 0.7
    
) {
    
    # This function finds the columns of <x> that correlate highly (above <threshold>) with
    # a set of <initial> or 'seed' columns, then scores them according to the number of high
    # correlations and the maximum correlation.  The data x should be a matrix or data frame
    # with column names, and the output will consist of a subset of these column names and their
    # corresponding scores.
    
    # Convert x to data.table, in case it isn't already:
    
    x <- as.data.table(x)
    
    cor_with_initial <- x[
        ,
        cor(
            .SD[, -..initial],
            .SD[, ..initial]
        )
    ]
    
    # In the following, <threshold> is used to specify what counts as a "high
    # correlation" - any correlation above this threshold is deemed "high".  This is used in
    # scoring the genes for correlation with the initial genes - specifically, we use the
    # number of "high" correlations with the initial genes.  I tried, using HNSC data,
    # replacing the number of "high" correlations with simply the average correlation with
    # the initial genes.  The gene set I got was more or less the same, though it seemed to
    # be missing some genes that I think are important, like FAP and ACTA2.  So I'll stick
    # with the first method, though perhaps I could make the averaging method available as an
    # argument of the function.
    
    above_threshold <- rbindlist(
        lapply(
            rownames(cor_with_initial),
            function(g) {
                
                nhc <- sum(
                    sapply(
                        initial,
                        function(s) {
                            switch(
                                (
                                    cor_with_initial[g, s] >= threshold
                                ) + 1,
                                0,
                                1
                            )
                        }
                    )
                )
                
                if(
                    nhc >= 1
                ) {
                    list(
                        id = g,
                        number_high_correlations = nhc,
                        highest_correlation = max(cor_with_initial[g, ])
                    )
                }
                
            }
        )
    )
    
    above_threshold[
        ,
        score := (number_high_correlations - min(number_high_correlations))/
            max(number_high_correlations) +
            (highest_correlation - min(highest_correlation))/
            max(highest_correlation)
    ]
    
    rbind(
        data.table(id = initial, score = Inf),
        above_threshold[order(-score), .(id, score)]
    )
    
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
    patient_var = NULL,
    groups = c('cancer', 'fibroblast'),
    to_keep = c(
        'SNAI1',
        'SNAI2',
        'TWIST1',
        'VIM',
        'ZEB1',
        'ZEB2'
    ),
    genes_filter_fun = function(x) {mean(x) > 0.2},
    cells_filter_fun = function(x) {mean(x) > 0.2},
    genes_detected_threshold = 1000,
    genes_order_fun = function(x) {
        as.vector(seriation::seriate(dist(t(x)), method = 'SPIN_NH')[[1]])
    },
    cells_order_fun = function(x) {seriation::seriate(dist(x), method = 'GW')[[1]]$order},
    score_cells = TRUE,
    scores_filter_fun = genes_filter_fun,
    min_sig_size = 40,
    seed = NULL
    
) {
    
    # It's important that the only non-numeric (i.e. non-gene) variables in <sc_data> are <id_var>
    # and <group_var>.
    
    # Note that if <score_cells> is set to TRUE, scores for the cells are calculated using the
    # scrabble::score() function, and the cells are then ordered by these scores - in particular,
    # the <cells_order_fun> argument is ignored.
    
    genes_filter_fun <- match.fun(genes_filter_fun)
    cells_filter_fun <- match.fun(cells_filter_fun)
    genes_order_fun <- match.fun(genes_order_fun)
    cells_order_fun <- match.fun(cells_order_fun)
    
    genes <- sort(
        unique(
            c(
                genes[
                    apply(
                        sc_data[
                            get(group_var) %in% groups,
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
    
    sc_data[
        get(group_var) %in% groups,
        genes_detected := sum(as.numeric(.SD) > 0),
        by = id_var,
        .SDcols = -c(group_var, patient_var)
    ]
    
    setkeyv(sc_data, id_var)
    
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
    
    if(score_cells) {
        
        all_genes_filtered <- sc_data[
            cells,
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
                        scrabble::score(
                            set_colnames(t(.SD[, -'id']), .SD$id),
                            list(sig_genes),
                            bin.control = TRUE,
                            nbin = length(all_genes_filtered) %/% 110
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
                    scrabble::score(
                        set_colnames(
                            t(
                                sc_data[
                                    cells,
                                    ..all_genes_filtered
                                ]
                            ),
                            sc_data[cells, id]
                        ),
                        list(sig_genes),
                        bin.control = TRUE,
                        nbin = length(all_genes_filtered) %/% 110
                        # n = 50
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
                signature_genes_for_scoring = sig_genes
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
    annotations_title = NULL,
    annotations_title_size = 14,
    annotations_margin = 3,
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
    h_upper_limit = 12,
    h_legend_title = 'Expression\nlevel',
    h_legend_width = 10,
    h_legend_breaks = c(0, 3, 6, 9, 12),
    h_legend_labels = c('0' = '0', '3' = '3', '6' = '6', '9' = '9', '12' = '\u2265 12'),
    h_legend_height = 15,
    h_legend_direction = c('vertical', 'horizontal'),
    es_fun = mean,
    es_plot_type = c('line', 'bar'),
    es_colours = h_colours,
    es_upper_limit = 4,
    es_legend_title = 'Average\nexpression',
    es_legend_width = 10,
    es_legend_breaks = c(0, 2, 4),
    es_legend_labels = c('0' = '0', '2' = '2', '4' = '\u2265 4'),
    es_runmean_window = NULL,
    gd_colours = colorRampPalette(RColorBrewer::brewer.pal(9, "RdPu"))(50),
    gd_limits = c(0, 8000),
    gd_breaks = c(0, 4000, 8000),
    gd_labels = c('0' = '0', '4000' = '4000', '8000' = '\u2265 8000'),
    gd_legend_title = 'Genes\ndetected',
    gd_legend_width = 10,
    gd_legend_height = 10,
    gd_legend_direction = c('vertical', 'horizontal')
    # gd_legend_margin = c(25, 0, 0, 0)
    
) {
    
    annotations_side <- match.arg(annotations_side)
    gd_legend_direction <- match.arg(gd_legend_direction)
    h_legend_direction <- match.arg(h_legend_direction)
    
    es_fun <- match.fun(es_fun)
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
                limits = c(0, h_upper_limit),
                breaks = h_legend_breaks,
                labels = h_legend_labels,
                colours = h_colours,
                oob = scales::squish
            ) +
            theme(
                legend.key.width = unit(h_legend_width, 'pt'),
                legend.key.height = unit(h_legend_height, 'pt'),
                legend.direction = h_legend_direction
            ) +
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
        axis_title = annotations_title,
        axis_title_size = annotations_title_size,
        title_edge_margin = annotations_margin
    )
    
    if(es_plot_type == 'bar') {
        
        if(is.null(es_upper_limit)) {
            es_upper_limit <- single_cell_data[
                cells,
                max(apply(.SD, 1, es_fun)),
                .SDcols = genes
            ]
        }
        
        expression_summary <- sapply(
            groups,
            function(grp) {
                ggplot(
                    single_cell_data[
                        cells
                    ][
                        get(group_var) == grp,
                        .(id = get(id_var), es = apply(.SD, 1, es_fun)),
                        .SDcols = genes
                    ],
                    aes(
                        x = factor(id, levels = id[ordering_cells[[grp]]]),
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
        
        es_upper_limit <- ceiling(
            2*max(
                sapply(
                    expression_summary_data,
                    function(x) x[, max(es)]
                )
            )
        )/2
        
        es_linegraph_breaks <- switch(
            (es_upper_limit < 2) + 1,
            seq(from = 0, to = es_upper_limit, by = floor(log2(es_upper_limit))),
            seq(from = 0, to = es_upper_limit, by = floor(log2(2*es_upper_limit))/2)
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
                        size = 1
                    ) +
                    scale_y_continuous(
                        limits = c(0, es_upper_limit),
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
                        # axis.text.y = switch(
                        #     (which(groups == grp) == 1) + 1,
                        #     element_blank(),
                        #     element_text()
                        # ),
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
                legend.direction = gd_legend_direction
            ) +
            labs(fill = gd_legend_title)
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
            list(
                annotations = nrow(single_cell_data)*default_figure_widths$annotations/
                    sum(unlist(default_figure_widths[groups]))
            ),
            single_cell_data[
                cells
            ][
                ,
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
        es_plot_type = es_plot_type
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
    es_y_axis_title = 'Average\nexpression',
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
                        grp,
                        size = 11,
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
                            sc_groups_figures$plots$expression_summary[[1]],
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
                lapply(
                    sc_groups_figures$plots$expression_summary,
                    function(g) {
                        g + theme(
                            axis.text.y = element_blank(),
                            axis.title = element_blank(),
                            axis.ticks.y = element_blank(),
                            axis.ticks.length = unit(0, 'pt')
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
                            axis.ticks.length = unit(0, 'pt')
                        )
                    }
                ),
                list(
                    ggdraw(
                        get_y_axis(
                            sc_groups_figures$plots$expression_summary[[
                                length(sc_groups_figures$groups)
                            ]],
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





sample_indices <- function(
    
    types_vec,
    counts_table
    
) {
    
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
                    
                    switch(
                        (s[[type]] > N) + 1,
                        sample(inds, s[[type]]),
                        c(
                            rep(inds, s[[type]] %/% N),
                            sample(inds, s[[type]] %% N)
                        )
                    )
                    
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
        
        # Start by taking a copy of (the relevant subset of) the data, so that we don't change
        # the columns of the original data table.
        
        data <- data[, c(..types_var, ..cols)]
        
        data[
            ,
            (cols) := lapply(.SD, normalise_fun),
            .SDcols = cols
        ]
        
    }
    
    if(is.null(initial_types)) {
        
        # If initial_types is not specified, we calculate the proportions and contributions of
        # all the types.  For this option, counts_table should have one column per type and no
        # other columns.
        
        types <- unique(data[[types_var]])
        
        data.table(
            simulation = unlist(
                lapply(
                    1:nrow(counts_table),
                    rep,
                    length(types)
                )
            ),
            type = rep(
                types,
                nrow(counts_table)
            ),
            proportion_content = unlist( # Calculate proportion for each type:
                as.data.table(
                    t(counts_table)
                )[
                    ,
                    lapply(
                        .SD,
                        function(s) {
                            s/sum(s)
                        }
                    )
                ]
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
                            plyr::mapvalues(totals, NA, 0, warn_missing = FALSE)/
                                sum(totals[!is.na(totals)])
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
                    function(s) {
                        data[
                            unlist(s, use.names = FALSE),
                            sum(.SD[get(types_var) %in% initial_types])/sum(.SD),
                            .SDcols = cols
                        ]
                    }
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
    add_zero_one = 0,
    x_axis_title = 'Proportion of tumour',
    y_axis_title = 'Proportion of gene expression',
    legend_title = 'Cell type',
    legend_labels = c(
        'b_cell' = 'B cell',
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
    plot_title = '',
    ...
    
) {
    
    # This is a wrapper function that runs the sampling and aggreating functions above for
    # several cell types, then makes a lineplot.
    
    # The add_zero_one argument is used to specify that you want to add points (0, 0) and/
    # or (1, 1) to the plot.  Supply 0 (the default) if you want the (0, 0) point; 1 if
    # you want (1, 1); and c(0, 1) if you want both.
    
    # The ... is for additional arguments to the simulate_counts() function.  I anticipate
    # it being used mostly to adjust max_mean_count to suit the supplied dataset.
    
    cell_types <- unique(single_cell_data[[types_var]])
    
    plot_data <- rbindlist(
        lapply(
            cell_types,
            function(ct) {
                
                simulated_counts <- simulate_counts(
                    cell_types,
                    initial_types = ct,
                    ...
                )
                
                sampled_indices <- sample_indices(
                    single_cell_data[[types_var]],
                    simulated_counts$counts_table[, ..cell_types]
                )
                
                proportions_table <- type_contrib(
                    single_cell_data[, c(..types_var, ..genes)],
                    genes,
                    simulated_counts$counts_table,
                    sampled_indices,
                    initial_types = ct
                )[
                    ,
                    initial_type := ct
                ]
                
                proportions_table
                
            }
        )
    )[
        ,
        .(mean_proportion_contrib = mean(proportion_contrib)),
        by = .(initial_type, proportion_content)
    ]
    
    if(!is.null(add_zero_one)) {
        
        plot_data <- rbind(
            plot_data,
            data.table(
                initial_type = rep(
                    unique(
                        plot_data$initial_type,
                        length(add_zero_one)
                    )
                ),
                proportion_content = unlist(
                    lapply(
                        add_zero_one,
                        rep,
                        length(unique(plot_data$initial_type))
                    )
                ),
                mean_proportion_contrib = unlist(
                    lapply(
                        add_zero_one,
                        rep,
                        length(unique(plot_data$initial_type))
                    )
                )
            )
        )[
            order(initial_type, proportion_content)
        ]
        
    }
    
    lineplot <- ggplot(
        data = plot_data,
        aes(
            x = proportion_content,
            y = mean_proportion_contrib,
            colour = initial_type
        )
    ) +
        geom_line() +
        geom_point(shape = 0) +
        theme(
            panel.background = element_blank(),
            panel.border = element_rect(fill = NA, colour = 'black', size = 1),
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
        scale_x_continuous(
            expand = c(0.02, 0.02)
        ) +
        scale_colour_manual(
            labels = switch(
                is.null(legend_labels) + 1,
                legend_labels,
                sort(unique(initial_type))
            ),
            values = switch(
                is.null(legend_colours) + 1,
                legend_colours,
                setNames(
                    hcl(
                        h = seq(
                            360/(2*length(cell_types)),
                            360 - 360/(2*length(cell_types)),
                            length.out = length(cell_types)
                        ),
                        c = 100,
                        l = 65
                    ),
                    sort(unique(initial_type))
                )
            )
        ) +
        labs(
            x = x_axis_title,
            y = y_axis_title,
            colour = legend_title,
            title = plot_title
        )
    
    list(
        data = plot_data,
        lineplot = lineplot
    )
    
}





simulated_tumours_data <- function(
    
    single_cell_data,
    types_var = 'cell_type',
    genes = names(single_cell_data[, -c('id', 'patient', 'cell_type')]),
    proportions_vec = seq(0.1, 0.9, length.out = 40),
    n = 25,
    id_prefix = '',
    show_progress = TRUE,
    restore_log = TRUE,
    ...
    
) {
    
    # The ... is for additional arguments to the simulate_counts() function.
    
    cell_types <- unique(single_cell_data[[types_var]])
    
    simulated_counts <- simulate_counts(
        cell_types,
        initial_types = 'cancer',
        proportions_vec = proportions_vec,
        n = n,
        ...
    )
    
    sampled_indices <- sample_indices(
        single_cell_data[[types_var]],
        simulated_counts$counts_table[, ..cell_types]
    )
    
    # The lapply() application in which we construct expr is much faster if we prepare the data
    # beforehand, i.e. subset genes and reverse the log.  Originally I simply used the following
    # very simple code to achieve this:
    
    # data_for_aggregating <- 2^single_cell_data[, ..genes] - 1
    
    # This is quite fast, but it requires taking a copy and thus uses a lot of memory - too much
    # in the case of Lambrechts's lung scRNA-seq data, when R runs out of memory and throws an
    # error.  The same is true for the following, in which lapply() presumably has to remember
    # everything it produces until the last iteration, when it returns it all at once:
    
    # single_cell_data[
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
    # available, calculate how much we would need for 2^single_cell_data[, ..genes] - 1, and
    # run the for loop only when there isn't enough available memory.
    
    for(g in genes) {
        set(
            single_cell_data,
            j = g,
            value = round(2^single_cell_data[[g]] - 1, 4)
        )
    }
    
    if(show_progress) {
        
        cat(
            length(sampled_indices),
            'profiles to build.  Currently on iteration:\n'
        )
        
        expr <- log2(
            rbindlist(
                lapply(
                    1:length(sampled_indices),
                    function(i) {
                        
                        cat(
                            rep('\b', ceiling(log10(i))),
                            i,
                            sep = ''
                        )
                        
                        as.list(
                            colSums(
                                single_cell_data[unlist(sampled_indices[[i]]), ..genes]
                            )
                        )
                        
                    }
                )
            ) + 1
        )
        
        cat('\nDone!')
        
    } else {
        
        expr <- log2(
            rbindlist(
                lapply(
                    1:length(sampled_indices),
                    function(i) {
                        as.list(
                            colSums(
                                single_cell_data[unlist(sampled_indices[[i]])]
                            )
                        )
                    }
                )
            ) + 1
        )
        
    }
    
    expr[
        ,
        id := paste0(id_prefix, .I)
    ]
    
    setcolorder(expr, 'id')
    
    meta <- simulated_counts$counts_table[
        ,
        .(purity = cancer/sum(as.numeric(.SD))),
        by = .(id = paste0(id_prefix, row.names(simulated_counts$counts_table))),
        .SDcols = cell_types
    ]
    
    if(restore_log) {
        for(g in genes) {
            set(
                single_cell_data,
                j = g,
                value = round(log2(single_cell_data[[g]] + 1), 4)
            )
        }
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
    # genes_filter_fun = function(x) 1:250,
    genes_filter_fun = function(x) {1:min(sum(x > 1), 250)}, # Take no more than 250 genes
    genes_filter_quantile = 0.5,
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
                
                genes_filtered <- genes_cor_with_initial_and_cell_types[
                    ,
                    .(
                        id = id,
                        score = max_initial_gene - max_cell_type
                    )
                ][
                    order(-score)
                ][
                    genes_filter_fun(score),
                    id
                ]
                
            } else if(max_initial_gene) {
                
                genes_filtered <- genes_cor_with_initial_and_cell_types[
                    ,
                    .(
                        id = id,
                        score = max_initial_gene
                    )
                ][
                    order(-score)
                ][
                    genes_filter_fun(score),
                    id
                ]
                
            } else { # Then max_cell_type must be TRUE (and max_initial_gene is FALSE)
                
                genes_filtered <- genes_cor_with_initial_and_cell_types[
                    ,
                    .(
                        id = id,
                        score = -max_cell_type
                    )
                ][
                    order(-score)
                ][
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
                        quantile,
                        genes_filter_quantile
                    )
                )
                
            }
            
            if(is.null(cell_type_weights)) {
                
                # Define weights for cell types as quantiles of their average correlations
                # with EMT/CAF markers:
                
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
                        quantile,
                        genes_filter_quantile
                    )
                )
                
            }
            
            if(isFALSE(cell_type_weights)) {
                
                genes_filtered <- genes_cor_with_initial_and_cell_types[
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
                ][
                    genes_filter_fun(score),
                    id
                ]
                
            } else if(isFALSE(initial_gene_weights)) {
                
                genes_filtered <- genes_cor_with_initial_and_cell_types[
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
                ][
                    genes_filter_fun(score),
                    id
                ]
                
            } else {
                
                genes_filtered <- genes_cor_with_initial_and_cell_types[
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
                ][
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
    
    # Ordering of corelation matrix for resid_data using SPIN_STS:
    
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
                ]
                
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
    legend_title = 'EMT-CAF\nscore'
    
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
    
    # For each cancer type/subtype, calculate EMT-CAF scores using correlations with the genes
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
    legend_title = 'EMT-CAF\nscore'
    
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
            
            # Running Wilcoxon rank-sum test for each gene:
            
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





nice_names <- function(names_vec) {
    
    plyr::mapvalues(
        names(ids_genes_ordering_list),
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
