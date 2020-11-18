headl <- function(x, m = 6L, n = 6L) {
    x[1:m, 1:n]
}

headr <- function(x, m = 6L, n = 6L) {
    x[1:m, (ncol(x) - n + 1):ncol(x)]
}

taill <- function(x, m = 6L, n = 6L) {
    x[(nrow(x) - m + 1):nrow(x), 1:n]
}

tailr <- function(x, m = 6L, n = 6L) {
    x[(nrow(x) - m + 1):nrow(x), (ncol(x) - n + 1):ncol(x)]
}

head2 <- function(x, m = 6L, n = 6L) {
    x[1:m, 1:n]
}

tail2 <- function(x, m = 6L, n = 6L) {
    x[(nrow(x) - m + 1):nrow(x), 1:n]
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





heat_map <- function(
    
    mat,
    ordering = NULL,
    axis_title_y = '',
    axis_title_x = axis_title_y,
    axis_text_y = colnames(mat),
    axis_text_x = rownames(mat),
    legend_title = 'correlation',
    colours = rev(colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))(50)),
    colour_limits = c(-1, 1),
    oob = scales::squish,
    legend_breaks = waiver(),
    legend_labels = waiver(),
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
    
    if(!is.null(ordering)) {
        g <- ggplot(
            melt(mat[ordering, ordering]),
            aes(x = Var1, y = Var2)
        )
    } else {
        g <- ggplot(
            melt(mat),
            aes(x = Var1, y = Var2)
        )
    }
    
    g <- g +
        geom_raster(aes(fill = value)) +
        scale_fill_gradientn(
            colours = colours,
            limits = colour_limits, #With the default colours, this makes zero appear as white.
            oob = oob,
            breaks = legend_breaks,
            labels = legend_labels
        ) +
        theme(
            axis.ticks = element_blank(),
            axis.ticks.length = unit(0, 'pt'),
            axis.text = element_text(size = axis_text_size),
            axis.text.y = element_text(vjust = 0.5), #vjust makes labels in centre of rows
            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5), #likewise columns
            plot.margin = unit(plot_margin, 'pt'),
            # The following eliminates the grey from the background.  This shouldn't be
            # necessary but sometimes the grey is still visible without this.
            panel.background = element_blank(),
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
            fill = legend_title,
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





dendro <- function(
    
    clust,
    edge = c('bottom', 'left', 'top', 'right'),
    plot_margin = sapply(
        c('top', 'right', 'bottom', 'left'),
        function(x) ifelse(x == edge, 0, 5.5)
    )
    
) {
    
    edge <- match.arg(edge)
    
    dend_data <- ggdendro::dendro_data(as.dendrogram(clust))
    
    if(edge == 'bottom') {
        
        ggplot(dend_data$segments) +
            geom_segment(aes(x = x, y = y, xend = xend, yend = yend)) +
            scale_x_continuous(expand = c(0, 0.5)) +
            scale_y_continuous(expand = c(0, 0)) +
            theme(
                axis.text = element_blank(),
                axis.ticks = element_blank(),
                axis.ticks.length = unit(0, 'pt'),
                axis.title = element_blank(),
                panel.background = element_rect(fill = NA),
                plot.margin = unit(plot_margin, 'pt')
            )
        
    } else if(edge == 'left') {
        
        ggplot(dend_data$segments) +
            geom_segment(aes(x = y, y = x, xend = yend, yend = xend)) +
            scale_x_continuous(expand = c(0, 0)) +
            scale_y_continuous(expand = c(0, 0.5)) +
            theme(
                axis.text = element_blank(),
                axis.ticks = element_blank(),
                axis.ticks.length = unit(0, 'pt'),
                axis.title = element_blank(),
                panel.background = element_rect(fill = NA),
                plot.margin = unit(plot_margin, 'pt')
            )
        
    } else if(edge == 'top') {
        
        ggplot(dend_data$segments) +
            geom_segment(aes(x = x, y = -y, xend = xend, yend = -yend)) +
            scale_x_continuous(expand = c(0, 0.5)) +
            scale_y_continuous(expand = c(0, 0)) +
            theme(
                axis.text = element_blank(),
                axis.ticks = element_blank(),
                axis.ticks.length = unit(0, 'pt'),
                axis.title = element_blank(),
                panel.background = element_rect(fill = NA),
                plot.margin = unit(plot_margin, 'pt')
            )
        
    } else if(edge == 'right') {
        
        ggplot(dend_data$segments) +
            geom_segment(aes(x = -y, y = x, xend = -yend, yend = xend)) +
            scale_x_continuous(expand = c(0, 0)) +
            scale_y_continuous(expand = c(0, 0.5)) +
            theme(
                axis.text = element_blank(),
                axis.ticks = element_blank(),
                axis.ticks.length = unit(0, 'pt'),
                axis.title = element_blank(),
                panel.background = element_rect(fill = NA),
                plot.margin = unit(plot_margin, 'pt')
            )
        
    }
    
}





heat_map_labels_repel <- function(
    
    labels_vec,
    edge = c('top', 'right', 'bottom', 'left'),
    text_size = 3,
    segment_size = 0.5,
    nudge = 0.35,
    axis_title = NULL,
    axis_title_size = 11, # Use 14 if you want it to look like a plot title
    title_edge_margin = 0, # In pts,
    force = 0.25,
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
                segment.size = segment_size,
                angle = 90,
                # nudge_y pushes the labels in the y direction up or down by a proportion of the plot area height.
                # We're using the negative of nudge here, because we want the labels to be pushed down from the top.
                nudge_y = -nudge,
                direction = 'x',
                force = force,
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
                segment.size = segment_size,
                nudge_x = -nudge,
                direction = 'y',
                force = force,
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
                segment.size = segment_size,
                angle = 90,
                nudge_y = nudge,
                direction = 'x',
                force = force,
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
                segment.size = segment_size,
                nudge_x = nudge,
                direction = 'y',
                force = force,
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
	legend_labels = waiver(),
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
			labels = legend_labels,
            oob = oob,
            breaks = legend_breaks
        ) +
        theme(
            axis.ticks = element_blank(),
            axis.ticks.length = unit(0, 'pt'),
            axis.text = element_blank(),
            plot.margin = unit(plot_margin, 'pt'),
            # The following eliminates the grey from the background.  This shouldn't be
            # necessary but sometimes the grey is still visible without this.
            panel.background = element_blank(),
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
    threshold_fun = function(x) 0.6
    
) {
    
    # x should be a matrix or data frame/table, among whose columns we will look for high
    # correlations (so the elements of initial should also be columns).  If it is a data
    # frame/table, it's important that there are no ID or non-numeric columns.
    
    FUN <- match.fun(FUN)
    threshold_fun <- match.fun(threshold_fun)
    
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
            correlation > threshold_fun(correlation)
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





# The following object takes an hclust object, cuts the tree at every height, then
# calculates the silhouette values for the resulting clusters for each height, producing
# a matrix with a silhouette value for each label and each height.  It returns either a
# list consisting of the cutree matrix and the matrix of silhouette values, or a numeric
# vector obtained by applying a function (<apply_fun>, mean by default) to the rows of
# the silhouette value matrix.  This vector is intended to characterise in general how
# well each label clusters with the others.  Note that we include the zero values in
# the application of <apply_fun>, since these also reflect how well the label clusters.
# This is because if a label gets a lot of zero silhouette values, it means it stays as
# it's own one-element cluster for a long time, never being close enough to another
# cluster to be combined with it during the agglomerative algorithm.  Such a label
# therefore clusters poorly in general, and should get a low score.  Including all the
# zeros in the application of <apply_fun> should reflect this.  By the way, we're
# assuming the hierarchical clustering is agglomerative, which it is in hclust().

silhouette_vals <- function(
    
    dist_mat,
    clust,
    labels = clust$labels,
    apply_fun = mean
    
) {
    
    if(!is.null(apply_fun)) {apply_fun <- match.fun(apply_fun)}
    
    # In the following, we remove the last element of hclust_genes$height because this is
    # a vector with one unique element, representing the largest cluster which includes
    # all genes.  This gives NA in silhouette_vals if we don't take it out.
    
    cutree_all <- cutree(
        clust,
        h = clust$height[-length(clust$height)]
    )
    
    if(is.null(clust$labels)) {
        rownames(cutree_all) <- labels
    }
    
    sil_vals <- apply(
        cutree_all,
        2,
        function(x) {
            cluster::silhouette(
                x,
                dist_mat
            )[, 'sil_width']
        }
    )
    
    # Could use magrittr::set_rownames() here, but prefer to avoid package dependency:
    
    rownames(sil_vals) <- rownames(cutree_all)
    
    if(is.null(apply_fun)) {
        list(
            sil_vals = sil_vals,
            cutree_clust_assign = cutree_all
        )
    } else {
        apply(sil_vals, 1, apply_fun)
    }
    
}





cov_mat_ellipse <- function(
    
    centre,
    data,
    level = 0.95,
    npoints = 100
    
) {
    
    # This function computes the (x, y) coordinates for a confidence ellipse, for use in
    # plotting.  Currently, it only works with 2-dimensional data.
    
    cov_mat <- cov(data)
    
    # Take the eigen decomposition of the covariance matrix.  The matrix of eigenvectors
    # then corresponds to the directions along which the data varies the most - the
    # eigenvectors wil be perpendicular, and, since R normalises them to length 1, the
    # matrix they define is also a rotation matrix.  The eigenvalues then represent the
    # "components" of the variance along the directions defined by the eigenvectors.
    
    eigen_decomp <- eigen(cov_mat)
    
    # Get coordinates for the ellipse centred at (0, 0) whose radii are aligned with the
    # x and y axes and have sizes defined by the eigenvalues (representing variances)
    # multiplied by a scaling factor (taken from a chi-squared distribution).  Do this
    # using the standard parametric representation of an ellipse applied to a vector of
    # angles.  Then transform these coordinates by applying the rotation and translation.
    
    angle_range <- seq(0, 2*pi, length.out = npoints)
    
    t(
        eigen_decomp$vectors %*% matrix(
            c(
                cos(angle_range)*sqrt(qchisq(level, 2)*eigen_decomp$values[1]),
                sin(angle_range)*sqrt(qchisq(level, 2)*eigen_decomp$values[2])
            ),
            nrow = 2,
            ncol = npoints,
            byrow = TRUE
        ) + centre
    )
    
}





check_point_ellipse <- function(
    
    point,
    centre,
    data,
    level = 0.95
    
) {
    
    # This function calculates computes a confidence ellipse from data and checks if a
    # given point is contained within the ellipse.  It returns TRUE or FALSE accordingly.
    
    # Get info for ellipse:
    
    cov_mat <- cov(data)
    
    eigen_decomp <- eigen(cov_mat)
    
    # Transform point back to base space:
    
    transformed_point <- as.numeric(
        solve(eigen_decomp$vectors) %*% (point - centre)
    )
    
    # Check if it's within transformed ellipse:
    
    (transformed_point[1]/sqrt(qchisq(level, 2)*eigen_decomp$values[1]))^2 +
        (transformed_point[2]/sqrt(qchisq(level, 2)*eigen_decomp$values[2]))^2
    
}





extend_endpoint <- function(centre, endpoint, data) {
    
    # This function computes a segment from <centre> to <endpoint>, extends it to the extreme
    # x or y values of <data>, and returns the end point of the extended segment.
    
    m_x <- (endpoint[1] - centre[1])/(endpoint[2] - centre[2])
    m_y <- 1/m_x
    
    x_fun <- function(yval) {centre[1] + (yval - centre[2])*m_x}
    y_fun <- function(xval) {centre[2] + (xval - centre[1])*m_y}
    
    x_sign <- sign(endpoint[1] - centre[1])
    y_sign <- sign(endpoint[2] - centre[2])
    
    x_end <- x_sign*max(x_sign*data[, 1])
    y_end <- y_sign*max(y_sign*data[, 2])
    
    cands <- matrix(
        c(x_end, y_fun(x_end), x_fun(y_end), y_end),
        2,
        2,
        byrow = TRUE
    )
    
    cands[
        apply(
            cands,
            1,
            function(v) {c(x_sign*v[1] <= x_sign*x_end && y_sign*v[2] <= y_sign*y_end)}
        ),
    ]
    
}





transform_segment <- function(centre, endpoint, data, suffix = '_t') {
    
    # This function transforms <data> via a rotation and translation, so that the segment
    # between <centre> and <endpoint> becomes the x axis in the transformed space.
    
    hypotenuse <- sqrt((endpoint[1] - centre[1])^2 + (endpoint[2] - centre[2])^2)
    
    rotation_matrix <- matrix(
        c(
            endpoint[1] - centre[1],
            endpoint[2] - centre[2],
            centre[2] - endpoint[2],
            endpoint[1] - centre[1]
        ),
        2,
        2,
        byrow = TRUE
    )/hypotenuse
    
    setNames(
        as.data.table(t(rotation_matrix %*% t(data) - centre)),
        paste0(names(data), suffix)
    )
    
}





fdr_bh <- function(pvals, threshold = 0.05) {
    
    # This function calculates the false discovery rate given a vector of p values and
    # a chosen threshold.  It answers the question: if I have a vector <pvals> and I
    # decide to declare all p values below <threshold> significant, then what is my
    # false discovery rate?  The calculation is based on the Benjamini-Hochberg
    # correction procedure.
    
    length(pvals)*max(pvals[pvals < threshold])/sum(pvals < threshold)
    
}

adjust_threshold_bh <- function(pvals, threshold = 0.05) {
    
    # This function takes a significance threshold and adjusts it via the Benjamini-
    # Hochberg method.  The output significance level is the value below which all
    # values of p.adjust(pvals, method = 'BH) would be considered significant.
    
    data.table(
        index = 1:length(pvals),
        pval = sort(pvals)
    )[
        ,
        switch(
            (sum(pval <= threshold*index/.N) == 0) + 1,
            threshold*max(which(pval <= threshold*index/.N))/.N,
            NA
        )
    ]
    
}





# The following is a simplified version of Julie's old scrabble::score() function.

signature_score <- function(
    mat,
    sig_genes,
    nbin = nrow(mat) %/% 110,
    n = 100,
    replace = FALSE,
    return_control_sets = FALSE
) {
    
    # <mat> should have genes for rows and samples for columns.  Set <return_control_sets>
    # to TRUE if you want to also return the control gene sets used to normalise the
    # expression levels.
    
    gene_averages <- sort(rowMeans(mat))
    
    bins <- setNames(
        cut(
            seq_along(gene_averages),
            breaks = nbin,
            labels = FALSE,
            include.lowest = TRUE
        ),
        names(gene_averages)
    )
    
    if(return_control_sets) {
        
        # Define control gene sets for distribution of scores:
        
        sig_gene_controls <- sapply(
            sig_genes,
            function(g) {
                sample(names(bins)[bins == bins[g]], n, replace = replace)
            },
            simplify = FALSE,
            USE.NAMES = TRUE
        )
        
        comparable_gene_sets <- lapply(
            1:n,
            function(i) sapply(sig_gene_controls, `[`, i)
        )
        
        sig_scores <- rowMeans(
            sapply(
                sig_genes,
                function(g) {
                    mat[g, ] - colMeans(
                        mat[sig_gene_controls[[g]], ]
                    )
                },
                USE.NAMES = TRUE
            )
        )
        
        return(
            list(
                scores = sig_scores,
                controls = sig_gene_controls,
                comparable_gene_sets = comparable_gene_sets # Redundant but possibly helpful
            )
        )
        
    } else {
        
        sig_scores <- rowMeans(
            sapply(
                sig_genes,
                function(g) {
                    mat[g, ] - colMeans(
                        mat[sample(names(bins)[bins == bins[g]], n, replace = replace), ]
                    )
                },
                USE.NAMES = TRUE
            )
        )
        
        return(sig_scores)
        
    }
    
}
