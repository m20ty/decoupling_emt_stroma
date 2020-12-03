# I don't think the following is as good as the ellipses method:

# Find genes associated with each cluster by transforming the PCA space (rotation and
# translation) so that one of the clusters is aligned with the x axis, then taking the
# genes that correlate highly with this new axis:

kclust_centre <- colMeans(kclust$centers)

find_endpoints <- function(centres, endpoints, data) {

    m_x <- (endpoints[1] - centres[1])/(endpoints[2] - centres[2])
    m_y <- 1/m_x

    x_fun <- function(yval) {centres[1] + (yval - centres[2])*m_x}
    y_fun <- function(xval) {centres[2] + (xval - centres[1])*m_y}

    x_sign <- sign(endpoints[1] - centres[1])
    y_sign <- sign(endpoints[2] - centres[2])

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

segment_endpoints <- as.data.table(
    kclust$centers,
    keep.rownames = 'index'
)[
    ,
    c('PC1_end', 'PC2_end') := as.list(
        find_endpoints(kclust_centre, as.numeric(.SD), plot_data[, .(PC1, PC2)])
    ),
    by = index
]

ggplot(plot_data, aes(PC1, PC2)) +
    geom_point(aes(colour = as.character(kmeans_cluster))) +
    geom_segment(
        aes(
            x = kclust_centre['PC1'],
            y = kclust_centre['PC2'],
            xend = PC1_end,
            yend = PC2_end
        ),
        data = segment_endpoints
    ) +
    theme_test() +
    labs(colour = 'k-means\ncluster')

transform_function <- function(centres, endpoints, data, suffix = '_t') {

    hypotenuse <- sqrt((endpoints[1] - centres[1])^2 + (endpoints[2] - centres[2])^2)

    rotation_matrix <- matrix(
        c(
            endpoints[1] - centres[1],
            endpoints[2] - centres[2],
            centres[2] - endpoints[2],
            endpoints[1] - centres[1]
        ),
        2,
        2,
        byrow = TRUE
    )/hypotenuse

    setNames(
        as.data.table(t(rotation_matrix %*% t(data) - centres)),
        paste0(names(data), suffix)
    )

}

plot_grid(
    plotlist = lapply(
        1:3,
        function(i) {
            ggplot(
                cbind(
                    plot_data,
                    transform_function(
                        kclust_centre,
                        kclust$centers[i, ],
                        plot_data[, .(PC1, PC2)]
                    )
                ),
                aes(PC1_t, PC2_t)
            ) +
                geom_point(aes(colour = as.character(kmeans_cluster))) +
                theme_test() +
                labs(colour = 'k-means\ncluster') +
                geom_hline(yintercept = 0)
        }
    ),
    nrow = 2,
    ncol = 2
)

kclust_gene_cor <- scores_data_transformed[
    names(head(sort(apply(rank_mat, 1, quantile, 0.25)), 100)),
    setNames(
        lapply(
            1:3,
            function(i) {
                cor(
                    as.numeric(.SD),
                    transform_function(
                        kclust_centre,
                        kclust$centers[i, ],
                        plot_data[, .(PC1, PC2)]
                    )$PC1_t
                )
            }
        ),
        paste0('clust', 1:3)
    ),
    by = gene,
    .SDcols = deconv_names
]

kclust_top_genes <- unique(
    unlist(
        lapply(
            names(kclust_gene_cor[, -'gene']),
            function(v) {
                kclust_gene_cor[
                    order(-get(v)),
                    head(gene, 20)
                ]
            }
        )
    )
)

set.seed(55413)

htmp_all_pca <- deconv_heatmap(
    scores_data_transformed[
        unique(
            c(
                names(tail(sort(apply(rank_mat, 1, quantile, 0.75)), 25)),
                kclust_top_genes
            )
        ),
        c('gene', ..deconv_names)
        ],
    order_genes_fun = 'seriate',
    order_genes_method = 'SPIN_NH',
    order_analyses_fun = 'seriate',
    order_analyses_method = 'SPIN_NH',
    plot_title = 'Common EMT and CAF genes'
)





# PCA to describe the variation in EMT and CAF genes between cancer types:

# Note if you calculate vr_av_score over all the genes and take the head and tail 100 genes,
# the k-means clustering perfectly separates squamous-like from non-squamous-like...

scores_pca <- prcomp(
    t(
        scores_data_transformed[
            names(head(sort(apply(rank_mat, 1, quantile, 0.25)), 100)),
            # vr_av_score[order(-score), head(gene, 50)],
            ..deconv_names
            ]
    )
)

pc_cor <- scores_data_transformed[
    names(head(sort(apply(rank_mat, 1, quantile, 0.25)), 100)),
    # vr_av_score[order(-score), head(gene, 50)],
    .(
        pc1_cor = cor(as.numeric(.SD), scores_pca$x[, 'PC1']),
        pc2_cor = cor(as.numeric(.SD), scores_pca$x[, 'PC2'])
    ),
    by = gene,
    .SDcols = deconv_names
    ]

plot_data <- as.data.table(
    scores_pca$x[, 1:2],
    keep.rownames = 'cancer_type'
)

squamoid_cancer_types <- deconv_names[grep('^[A-Z][A-Z]SC|Squam', deconv_names)]

plot_data[
    ,
    sa := switch(
        (cancer_type %in% squamoid_cancer_types) + 1,
        'Adeno-like',
        'Squamous-like'
    ),
    by = cancer_type
    ]

set.seed(16918)

kclust <- kmeans(plot_data[, .(PC1, PC2)], 3)

plot_data[
    ,
    kmeans_cluster := kclust$cluster
    ]

# ggplot(plot_data, aes(-PC1, PC2)) +
#     geom_point() +
#     geom_text(aes(label = cancer_type)) +
#     theme_test()

# In the following, it looks like PC1 + PC2 would do a very good job of separating Squamous
# and adeno:

pdf(
    '../data_and_figures/emt_caf_scores_pca_plot.pdf',
    width = 5,
    height = 3
)

ggplot(plot_data, aes(-PC1, PC2)) + # Minus sign to be consistent with heat map bars
    geom_point(aes(colour = sa)) +
    # geom_text(aes(label = cancer_type)) +
    theme_test() +
    labs(colour = NULL)

dev.off()

ggplot(plot_data, aes(-PC1, PC2)) +
    geom_point(aes(colour = as.character(kmeans_cluster))) +
    theme_test() +
    labs(colour = NULL)

kclust_centre <- colMeans(kclust$centers)

find_endpoints <- function(centres, endpoints, data) {
    
    m_x <- (endpoints[1] - centres[1])/(endpoints[2] - centres[2])
    m_y <- 1/m_x
    
    x_fun <- function(yval) {centres[1] + (yval - centres[2])*m_x}
    y_fun <- function(xval) {centres[2] + (xval - centres[1])*m_y}
    
    x_sign <- sign(endpoints[1] - centres[1])
    y_sign <- sign(endpoints[2] - centres[2])
    
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

segment_endpoints <- as.data.table(
    kclust$centers,
    keep.rownames = 'index'
)[
    ,
    c('PC1_end', 'PC2_end') := as.list(
        find_endpoints(kclust_centre, as.numeric(.SD), plot_data[, .(PC1, PC2)])
    ),
    by = index
    ]

ggplot(plot_data, aes(PC1, PC2)) +
    geom_point(aes(colour = as.character(kmeans_cluster))) +
    geom_segment(
        aes(
            x = kclust_centre['PC1'],
            y = kclust_centre['PC2'],
            xend = PC1_end,
            yend = PC2_end
        ),
        data = segment_endpoints
    ) +
    theme_test() +
    labs(colour = 'k-means\ncluster')

transform_function <- function(centres, endpoints, data, suffix = '_t') {
    
    hypotenuse <- sqrt((endpoints[1] - centres[1])^2 + (endpoints[2] - centres[2])^2)
    
    rotation_matrix <- matrix(
        c(
            endpoints[1] - centres[1],
            endpoints[2] - centres[2],
            centres[2] - endpoints[2],
            endpoints[1] - centres[1]
        ),
        2,
        2,
        byrow = TRUE
    )/hypotenuse
    
    setNames(
        as.data.table(t(rotation_matrix %*% t(data) - centres)),
        paste0(names(data), suffix)
    )
    
}

plot_grid(
    plotlist = lapply(
        1:3,
        function(i) {
            ggplot(
                cbind(
                    plot_data,
                    transform_function(
                        kclust_centre,
                        kclust$centers[i, ],
                        plot_data[, .(PC1, PC2)]
                    )
                ),
                aes(PC1_t, PC2_t)
            ) +
                geom_point(aes(colour = as.character(kmeans_cluster))) +
                theme_test() +
                labs(colour = 'k-means\ncluster') +
                geom_hline(yintercept = 0)
        }
    ),
    nrow = 2,
    ncol = 2
)

kclust_gene_cor <- scores_data_transformed[
    names(head(sort(apply(rank_mat, 1, quantile, 0.25)), 100)),
    setNames(
        lapply(
            1:3,
            function(i) {
                cor(
                    as.numeric(.SD),
                    transform_function(
                        kclust_centre,
                        kclust$centers[i, ],
                        plot_data[, .(PC1, PC2)]
                    )$PC1_t
                )
            }
        ),
        paste0('clust', 1:3)
    ),
    by = gene,
    .SDcols = deconv_names
    ]

kclust_top_genes <- unique(
    unlist(
        lapply(
            names(kclust_gene_cor[, -'gene']),
            function(v) {
                kclust_gene_cor[
                    order(-get(v)),
                    head(gene, 20)
                    ]
            }
        )
    )
)

set.seed(55413)

# Could just use the EMT genes here, and leave out the CAF ones...

htmp_all_pca <- deconv_heatmap(
    scores_data_transformed[
        unique(
            c(
                names(tail(sort(apply(rank_mat, 1, quantile, 0.75)), 25)),
                kclust_top_genes
            )
        ),
        c('gene', ..deconv_names)
        ],
    order_genes_fun = 'seriate',
    order_genes_method = 'SPIN_NH',
    order_analyses_fun = 'seriate',
    order_analyses_method = 'SPIN_NH',
    plot_title = 'Common EMT and CAF genes'
)

# Can use ggalt::geom_encircle(aes(group = kmeans_cluster)) to draw boundaries around the
# clusters.  Alternatively, there's an example of how to manually draw ellipses around the
# clusters on this page:

# https://stackoverflow.com/questions/23580095/how-to-plot-clusters-with-a-matrix

# Here's an implementation:

kclust_ellipses <- rbindlist(
    lapply(
        1:3,
        function(i) {
            data.table(
                kmeans_cluster = i,
                ellipse::ellipse(
                    plot_data[kmeans_cluster == i, cov(.SD), .SDcols = c('PC1', 'PC2')],
                    centre = kclust$centers[i, ]
                )
            )
        }
    )
)

ggplot(plot_data, aes(PC1, PC2)) +
    geom_point(aes(colour = as.character(kmeans_cluster))) +
    geom_segment(
        aes(
            x = kclust_centre['PC1'],
            y = kclust_centre['PC2'],
            xend = PC1_end,
            yend = PC2_end
        ),
        data = segment_endpoints
    ) +
    geom_path(
        aes(x = PC1, y = PC2, colour = as.character(kmeans_cluster)),
        data = kclust_ellipses
    ) +
    theme_test() +
    labs(colour = 'k-means\ncluster')

# I could abandon the rotation business and just use these ellipses.  I could take the cancer
# types unique to each ellipse (i.e. remove those in the overlaps) and just to a statistical
# test to find genes which score higher in one group than in the other two, regardless of
# correlation with any axis in the PC space.

# See the following for useful info on how to draw a confidence ellipse:

# https://www.visiondummy.com/2014/04/draw-error-ellipse-representing-covariance-matrix/

cov_mat_ellipse <- function(
    
    centres,
    data,
    level = 0.95,
    npoints = 100
    
) {
    
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
        ) + centres
    )
    
}

check_point_ellipse <- function(
    
    point,
    centres,
    data,
    level = 0.95
    
) {
    
    # Get info for ellipse:
    
    cov_mat <- cov(data)
    
    eigen_decomp <- eigen(cov_mat)
    
    # Transform point back to base space:
    
    transformed_point <- as.numeric(
        solve(eigen_decomp$vectors) %*% (point - centres)
    )
    
    # Check if it's within transformed ellipse:
    
    (transformed_point[1]/sqrt(qchisq(level, 2)*eigen_decomp$values[1]))^2 +
        (transformed_point[2]/sqrt(qchisq(level, 2)*eigen_decomp$values[2]))^2
    
}

ellipses_manual <- rbindlist(
    lapply(
        1:3,
        function(i) {
            cbind(
                kmeans_cluster = i,
                setNames(
                    as.data.table(
                        cov_mat_ellipse(
                            kclust$centers[i, ],
                            plot_data[kmeans_cluster == i, .(PC1, PC2)]
                        )
                    ),
                    c('X', 'Y')
                )
            )
        }
    )
)

ggplot(plot_data, aes(PC1, PC2)) +
    geom_point(aes(colour = as.character(kmeans_cluster))) +
    geom_segment(
        aes(
            x = kclust_centre['PC1'],
            y = kclust_centre['PC2'],
            xend = PC1_end,
            yend = PC2_end
        ),
        data = segment_endpoints
    ) +
    geom_path(
        aes(x = X, y = Y, colour = as.character(kmeans_cluster)),
        data = ellipses_manual
    ) +
    theme_test() +
    labs(colour = 'k-means\ncluster')

# I've never used the pipe operator in the following context before, but it's awesome!
# It allows me to keep the value of kmeans_cluster for each cancer type and then put
# it into plot_data again, via a custom function (which data.table calls an
# "anonymous" function and says it should be parenthesised).

plot_data[
    ,
    in_multiple_clusters := kmeans_cluster %>% (
        function(i) {
            
            other_clusters <- (1:3)[i != 1:3]
            
            check_point_ellipse(
                c(PC1, PC2),
                kclust$centers[other_clusters[1], ],
                plot_data[kmeans_cluster == other_clusters[1], .(PC1, PC2)]
            ) <= 1 | check_point_ellipse(
                c(PC1, PC2),
                kclust$centers[other_clusters[2], ],
                plot_data[kmeans_cluster == other_clusters[2], .(PC1, PC2)]
            ) <= 1
            
        }
    ),
    by = cancer_type
    ]

ggplot(plot_data, aes(PC1, PC2)) +
    geom_point() +
    geom_point(data = plot_data[in_multiple_clusters == TRUE], colour = 'red') +
    geom_segment(
        aes(
            x = kclust_centre['PC1'],
            y = kclust_centre['PC2'],
            xend = PC1_end,
            yend = PC2_end
        ),
        data = segment_endpoints
    ) +
    geom_path(
        aes(x = X, y = Y, colour = as.character(kmeans_cluster)),
        data = ellipses_manual
    ) +
    theme_test() +
    labs(colour = 'k-means\ncluster')

cancer_types_distinct <- plot_data[
    in_multiple_clusters == FALSE,
    cancer_type
    ]

kmeans_clust_distinct_genes <- lapply(
    1:3,
    function(i) {
        
        cancer_types_i <- plot_data[
            in_multiple_clusters == FALSE & kmeans_cluster == i,
            cancer_type
            ]
        
        cancer_types_not_i <- cancer_types_distinct[
            !(cancer_types_distinct %in% cancer_types_i)
            ]
        
        scores_data_transformed[
            names(head(sort(apply(rank_mat, 1, quantile, 0.25)), 100)),
            .(
                gene = gene,
                score_diff = rowMeans(.SD[, ..cancer_types_i]) -
                    rowMeans(.SD[, ..cancer_types_not_i])
            )
            ][
                order(-score_diff)
                ][
                    1:20,
                    gene
                    ]
        
    }
)

# One more heatmap...  I'm doing this one with hierarchical clustering, because we've
# filtered the genes and cancer types so that this will work well.  But using the same
# genes with all of deconv_names and SPIN_NH also works well.

htmp_all_kmeans_distinct <- deconv_heatmap(
    scores_data_transformed[
        unique(
            c(
                names(tail(sort(apply(rank_mat, 1, quantile, 0.75)), 25)),
                unique(unlist(kmeans_clust_distinct_genes))
            )
        ),
        c('gene', ..cancer_types_distinct)
        ],
    order_genes_fun = 'seriate',
    order_genes_method = 'GW_average',
    order_analyses_fun = 'seriate',
    order_analyses_method = 'GW_average',
    plot_title = NULL
)

deconv_heatmap_dendro_plot(htmp_all_kmeans_distinct)





# A rotation of the PCA space by -45 degrees (i.e. x goes to (PC1 + PC2)/sqrt(2) and y to
# (PC2 - PC1)/sqrt(2)) will make the x axis very effectively separate Squamous from adeno,
# while the y axis separates gynaecological from gastro-intestinal among the adenocarcinomas.
# I think the colour bars look better when using just the principal components, without
# rotation.  If we do the rotation and then use x and y for the colour bars, they look pretty
# much like reflections of each other, while the PCs alone do a good job of highlighting each
# of the three groups.  Note colour bars definitely work better than line charts, which are
# very jagged.

htmp_bar_pc1 <- heat_map_bar(
    -scores_pca$x[htmp_all$analyses, 'PC1'], # Minus sign makes green low and pink high
    htmp_all$ordering_analyses,
    colours = colorspace::terrain_hcl(12, c = c(80, 30), l = c(30, 90), power = c(0.2, 1.5)),
    colour_limits = c(-2, 2)
)

htmp_bar_pc2 <- heat_map_bar(
    scores_pca$x[htmp_all$analyses, 'PC2'],
    htmp_all$ordering_analyses,
    colours = colorspace::terrain_hcl(12, c = c(80, 30), l = c(30, 90), power = c(0.2, 1.5)),
    colour_limits = c(-2, 2)
)

pdf(
    '../data_and_figures/commonality_heatmap_pca_annotated.pdf',
    width = 7.5,
    height = 11
)

egg::ggarrange(
    htmp_bar_pc1 + theme(
        legend.position = 'none',
        plot.margin = unit(c(5.5, 5.5, 1, 5.5), 'pt')
    ) + labs(x = NULL),
    htmp_bar_pc2 + theme(plot.margin = unit(c(1, 5.5, 1, 5.5), 'pt')) + labs(x = NULL),
    htmp_all$heatmap + theme(plot.title = element_blank(), plot.margin = unit(c(1, 5.5, 5.5, 5.5), 'pt')),
    nrow = 3,
    ncol = 1,
    heights = c(1, 1, 40),
    newpage = FALSE
)

dev.off()

knitr::kable(
    data.frame(
        squamous = pc_cor[order(-pc1_cor), head(gene, 15)],
        gynae = pc_cor[order(pc1_cor), head(gene, 15)],
        gastro = pc_cor[order(-pc2_cor), head(gene, 15)]
    )
)

pc_genes <- unique(
    c(
        pc_cor[order(-pc1_cor), head(gene, 20)],
        pc_cor[order(pc1_cor), head(gene, 20)],
        pc_cor[order(-pc2_cor), head(gene, 20)],
        pc_cor[order(pc2_cor), head(gene, 20)]
    )
)





# I'm leaving out the following for now, because I need to figure out how to plot the
# colour bars and the dendrograms together on the same heatmap.  I'm not sure it's
# necessary, though, because we have the separate PC plot.

# Make colour bars of PC values with which to annotate the commonality heatmap:

htmp_bar_pc1 <- heat_map_bar(
    scores_pca$x[htmp_all_kmeans_distinct$analyses, 'PC1'],
    htmp_all_kmeans_distinct$ordering_analyses,
    colours = colorspace::terrain_hcl(
        12,
        c = c(80, 30),
        l = c(30, 90),
        power = c(0.2, 1.5)
    ),
    colour_limits = c(-2, 2)
)

htmp_bar_pc2 <- heat_map_bar(
    scores_pca$x[htmp_all_kmeans_distinct$analyses, 'PC2'],
    htmp_all_kmeans_distinct$ordering_analyses,
    colours = colorspace::terrain_hcl(
        12,
        c = c(80, 30),
        l = c(30, 90),
        power = c(0.2, 1.5)
    ),
    colour_limits = c(-2, 2)
)

egg::ggarrange(
    htmp_bar_pc1 + theme(
        legend.position = 'none',
        plot.margin = unit(c(5.5, 5.5, 1, 5.5), 'pt')
    ) + labs(x = NULL),
    htmp_bar_pc2 +
        theme(plot.margin = unit(c(1, 5.5, 1, 5.5), 'pt')) +
        labs(x = NULL),
    htmp_all_kmeans_distinct$heatmap +
        theme(
            plot.title = element_blank(),
            plot.margin = unit(c(1, 5.5, 5.5, 5.5), 'pt')
        ),
    nrow = 3,
    ncol = 1,
    heights = c(1, 1, 40),
    newpage = FALSE
)





set.seed(61937)

htmp_emt_caf <- deconv_heatmap(
    scores_data_transformed[
        unique(
            c(
                names(head(sort(apply(rank_mat, 1, quantile, 0.25)), 30)),
                names(tail(sort(apply(rank_mat, 1, quantile, 0.75)), 30))
                # 'SNAI1',
                # 'SNAI2',
                # 'ZEB1',
                # 'ZEB2',
                # 'TWIST1',
                # 'VIM'
            )
        ),
        c('gene', ..deconv_names)
        ],
    order_genes_fun = 'seriate',
    order_genes_method = 'SPIN_NH',
    order_analyses_fun = 'seriate',
    order_analyses_method = 'SPIN_NH',
    plot_title = 'Common EMT and stroma genes'
)

vr_av_score <- scores_data_transformed[
    names(head(sort(apply(rank_mat, 1, quantile, 0.25)), 100)),
    c('gene', ..deconv_names)
    ][
        ,
        .(
            vr = var(as.numeric(.SD)),
            av = sum(as.numeric(.SD)[as.numeric(.SD) > 0])/ncol(.SD)
        ),
        by = gene
        ][
            ,
            score := scale(vr) + scale(av)
            ]

set.seed(948136)

# In the following, I'm using my variance-average scores to choose the top 60 EMT genes,
# but actually I get similar genes with the command:

# names(head(sort(apply(rank_mat, 1, quantile, 0.25)), 60))

htmp_emt <- deconv_heatmap(
    scores_data_transformed[
        vr_av_score[order(-score), head(gene, 60)],
        c('gene', ..deconv_names)
        ],
    order_genes_fun = 'seriate',
    order_genes_method = 'SPIN_NH',
    order_analyses_fun = 'seriate',
    order_analyses_method = 'SPIN_NH',
    plot_title = 'Common EMT genes'
)

# Write heatmaps to one PDF:

pdf(
    '../data_and_figures/commonality_heatmaps.pdf',
    width = 7.5,
    height = 9
)

htmp_emt_caf$heatmap
htmp_emt$heatmap

dev.off()

# Together on the same page:

pdf(
    '../data_and_figures/commonality_heatmaps_onepage.pdf',
    width = 14,
    height = 9
)

plot_grid(
    htmp_emt_caf$heatmap + theme(legend.position = 'none'),
    htmp_emt$heatmap,
    nrow = 1,
    ncol = 2,
    rel_widths = c(0.865, 1),
    align = 'h'
)

dev.off()

# Heatmap with 50 EMT genes and 25 CAF genes:

set.seed(39384)

htmp_all <- deconv_heatmap(
    scores_data_transformed[
        unique(
            c(
                names(tail(sort(apply(rank_mat, 1, quantile, 0.75)), 25)),
                vr_av_score[order(-score), head(gene, 50)]
            )
        ),
        c('gene', ..deconv_names)
        ],
    order_genes_fun = 'seriate',
    order_genes_method = 'SPIN_NH',
    order_analyses_fun = 'seriate',
    order_analyses_method = 'SPIN_NH',
    plot_title = 'Common EMT and CAF genes'
)

pdf(
    '../data_and_figures/commonality_heatmap_all.pdf',
    width = 7.5,
    height = 11
)

htmp_all$heatmap

dev.off()
