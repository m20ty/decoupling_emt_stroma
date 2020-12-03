# This document is a bin for all my attempts at making the perfect commonality heatmaps.
# Some worked well, most did not.





# The following are the most recent attempts, which were quite decent, but I chose the
# version using SPIN_NH instead of hierarchical clustering.

htmp_emt_prelim <- deconv_heatmap(
    scores_data_transformed[
        names(head(sort(apply(rank_mat, 1, quantile, 0.25)), 100)),
        c('gene', ..deconv_names)
        ]
)

htmp_emt_prelim_gw <- deconv_heatmap(
    scores_data_transformed[
        names(head(sort(apply(rank_mat, 1, quantile, 0.25)), 100)),
        c('gene', ..deconv_names)
        ],
    order_genes_fun = 'seriate',
    order_genes_method = 'GW_complete',
    order_analyses_fun = 'seriate',
    order_analyses_method = 'GW_complete'
)

# Just to examine these preliminary heatmaps:

deconv_heatmap_dendro_plot(
    htmp_emt_prelim,
    direction = 'horizontal',
    title = 'EMT-CAF\nscore',
    title.position = 'top',
    title.hjust = 0.5,
    barwidth = unit(70, 'pt'),
    barheight = unit(10, 'pt')
)

deconv_heatmap_dendro_plot(
    htmp_emt_prelim_gw,
    direction = 'horizontal',
    title = 'EMT-CAF\nscore',
    title.position = 'top',
    title.hjust = 0.5,
    barwidth = unit(70, 'pt'),
    barheight = unit(10, 'pt')
)

# Silhouette values to determine optimal number of clusters:

sil_vals_emt <- silhouette_vals(
    htmp_emt_prelim$dist_mat_genes,
    htmp_emt_prelim$hclust_genes,
    labels = htmp_emt_prelim$genes,
    apply_fun = NULL
)

sil_data <- as.data.table(
    merge(
        melt(
            sil_vals_emt$sil_vals,
            varnames = c('gene', 'height'),
            value.name = 'sil_val'
        ),
        melt(
            sil_vals_emt$cutree_clust_assign,
            varnames = c('gene', 'height'),
            value.name = 'clust_assign'
        )
    )
)

sil_vals_emt_gw <- silhouette_vals(
    htmp_emt_prelim_gw$dist_mat_genes,
    htmp_emt_prelim_gw$hclust_genes,
    labels = htmp_emt_prelim_gw$genes,
    apply_fun = NULL
)

sil_data_gw <- as.data.table(
    merge(
        melt(
            sil_vals_emt_gw$sil_vals,
            varnames = c('gene', 'height'),
            value.name = 'sil_val'
        ),
        melt(
            sil_vals_emt_gw$cutree_clust_assign,
            varnames = c('gene', 'height'),
            value.name = 'clust_assign'
        )
    )
)

# To view the silhouette value plots:

# plot_grid(
#     plotlist = lapply(
#         1:10,
#         function(i) {
#             qplot(
#                 factor(gene, levels = gene),
#                 sil_val,
#                 data = sil_data[ # Or sil_data_gw
#                     height == sort(unique(height), decreasing = TRUE)[i]
#                 ][
#                     order(clust_assign, sil_val)
#                 ],
#                 geom = 'col'
#             ) +
#                 facet_grid(
#                     cols = vars(clust_assign),
#                     scales = 'free',
#                     space = 'free'
#                 ) +
#                 geom_hline(
#                     yintercept = sil_data[
#                         height == sort(unique(height), decreasing = TRUE)[i],
#                         mean(sil_val)
#                     ]
#                 )
#         }
#     ),
#     nrow = 5,
#     ncol = 2
# )

# I used these plots to decide the 'optimal' number of clusters for each case.

genes_for_heatmap <- sil_data[
    height == sort(unique(height), decreasing = TRUE)[2]
    ][
        order(-sil_val),
        head(as.character(gene), 50)
        ]

genes_for_heatmap_gw <- sil_data_gw[
    height == sort(unique(height), decreasing = TRUE)[3]
    ][
        order(-sil_val),
        head(as.character(gene), 50)
        ]

# Remove genes that stand on their own even at a very high level, then find the optimum
# number of clusters, then filter again:

sil_data_filtered <- sil_data[
    gene %in% rownames(sil_vals_emt$sil_vals)[
        apply(sil_vals_emt$sil_vals, 1, function(x) sum(x == 0)) < 50
        ]
    ]

sil_data_filtered_gw <- sil_data_gw[
    gene %in% rownames(sil_vals_emt_gw$sil_vals)[
        apply(sil_vals_emt_gw$sil_vals, 1, function(x) sum(x == 0)) < 50
        ]
    ]

# I chose the value of 50 based on the following plots, in which the trend seems to get
# steeper at around y = 50:

# plot(sort(apply(sil_vals_emt$sil_vals, 1, function(x) sum(x == 0))))
# plot(sort(apply(sil_vals_emt_gw$sil_vals, 1, function(x) sum(x == 0))))

# After this, different values for i in the following don't always give different numbers
# of clusters, since we've got rid of some genes.

genes_for_heatmap_filtered <- sil_data_filtered[
    height == sort(unique(height), decreasing = TRUE)[4]
    ][
        order(-sil_val),
        head(as.character(gene), 50)
        ]

genes_for_heatmap_filtered_gw <- sil_data_filtered_gw[
    height == sort(unique(height), decreasing = TRUE)[3]
    ][
        order(-sil_val),
        head(as.character(gene), 50)
        ]

# Filtering by taking highest-scoring genes in each of the x axis clusters:

dt <- scores_data_transformed[
    names(head(sort(apply(rank_mat, 1, quantile, 0.25)), 100)),
    c('gene', ..deconv_names)
    ]

cutree_analyses <- cutree(
    htmp_emt_prelim$hclust_analyses,
    k = 3
)

cutree_analyses_gw <- cutree(
    htmp_emt_prelim_gw$hclust_analyses,
    k = 3
)

genes_for_heatmap_3 <- unique(
    unlist(
        lapply(
            lapply(
                1:3,
                function(i) {
                    dt[
                        ,
                        rowMeans(magrittr::set_rownames(as.matrix(.SD), gene)),
                        .SDcols = names(cutree_analyses[cutree_analyses == i])
                        ]
                }
            ),
            function(x) names(x[x >= quantile(x, 0.75)])
        )
    )
)

genes_for_heatmap_3_gw <- unique(
    unlist(
        lapply(
            lapply(
                1:3,
                function(i) {
                    dt[
                        ,
                        rowMeans(magrittr::set_rownames(as.matrix(.SD), gene)),
                        .SDcols = names(cutree_analyses_gw[cutree_analyses_gw == i])
                        ]
                }
            ),
            function(x) names(x[x >= quantile(x, 0.75)])
        )
    )
)

# Also try just taking the most variable genes:

genes_for_heatmap_var <- scores_data_transformed[
    names(head(sort(apply(rank_mat, 1, quantile, 0.25)), 100)),
    c('gene', ..deconv_names)
    ][
        ,
        .(vr = var(as.numeric(.SD))),
        by = gene
        ][
            order(-vr),
            head(gene, 60)
            ]

# This doesn't depend on the clustering, so it doesn't matter if we use GW_complete or not, but
# just to make it fit with the other cases we'll define genes_for_heatmap_var_gw as well:

genes_for_heatmap_var_gw <- genes_for_heatmap_var

# Make heatmaps:

heatmaps_emt <- lapply(
    list(
        'genes_for_heatmap',
        'genes_for_heatmap_gw',
        'genes_for_heatmap_filtered',
        'genes_for_heatmap_filtered_gw',
        'genes_for_heatmap_3',
        'genes_for_heatmap_3_gw',
        'genes_for_heatmap_var',
        'genes_for_heatmap_var_gw'
    ),
    function(gene_set) {
        if(endsWith(gene_set, 'gw')) {
            deconv_heatmap(
                scores_data_transformed[
                    get(gene_set),
                    c('gene', ..deconv_names)
                    ],
                order_genes_fun = 'seriate',
                order_genes_method = 'GW_complete',
                order_analyses_fun = 'seriate',
                order_analyses_method = 'GW_complete'
            )
        } else {
            deconv_heatmap(
                scores_data_transformed[
                    get(gene_set),
                    c('gene', ..deconv_names)
                    ]
            )
        }
    }
)

# I'm beginning to think it might be best not to use hierarchical clustering, since I
# can't seem to get a robust one.  I could instead use PCA or SPIN_NH and present it as
# more of a spectrum.  At least, if we don't use a clustering algorithm then we're not
# obliged to present a dendrogram or anything - that is, we're not committing to any
# clustering.  SPIN_NH works quite well with genes chosen by a score based on variance
# and average (positive) EMT-CAF score, which we calculate in the following:

vr_av_score <- scores_data_transformed[
    names(head(sort(apply(rank_mat, 1, quantile, 0.25)), 100)),
    c('gene', ..deconv_names)
    ][
        ,
        .(
            vr = var(as.numeric(.SD)),
            av = mean(as.numeric(.SD)[as.numeric(.SD) > 0])
            # Maybe av should be:
            # av = mean(as.numeric(.SD)[sum(as.numeric(.SD))/ncol(.SD) > 0])
            # This effectively counts the negative scores as zeros.
        ),
        by = gene
        ][
            ,
            score := scale(vr) + scale(av)
            ]

# An example of one without hierarchical clustering (I think using score works better
# than either vr or av):

set.seed(948136)

deconv_heatmap(
    scores_data_transformed[
        vr_av_score[order(-score), head(gene, 50)],
        c('gene', ..deconv_names)
        ],
    order_genes_fun = 'seriate',
    order_genes_method = 'SPIN_NH',
    order_analyses_fun = 'seriate',
    order_analyses_method = 'SPIN_NH'
)$heatmap





# Below are the earlier, rougher attempts.

deconv_names <- names(deconv_data)[
    !grepl('^UVM|^MESO|^ACC|^SKCM|^CHOL|^THCA', names(deconv_data))
    ]

# In the following, I think the 25th percentile is sensible to use, because it prevents
# exclusion of genes which appear a lot at either end.  For example, SNAI2 commonly
# appears at the cancer end, but it also commonly appears at the CAF end.  We still want
# to include it because of its association with the cancer end, but if we used e.g. mean
# then it would be excluded because the average position would be roughly in the middle.
# This is not a problem for e.g. LAMC2, which appears at the cancer end much more often,
# so the mean would still pick up LAMC2.  But the 25th percentile should encourage
# inclusion of genes which appear at the cancer end often enough to be interesting but
# not necessarily the majority of the time.  

htmp_emt <- deconv_heatmap(
    scores_data_transformed[
        names(head(sort(apply(rank_mat, 1, quantile, 0.25)), 100)),
        c('gene', ..deconv_names)
        ]
    # order_genes_fun = 'seriate',
    # order_genes_method = 'GW_complete',
    # order_analyses_fun = 'seriate',
    # order_analyses_method = 'GW_complete'
)

deconv_heatmap_dendro_plot(
    htmp_emt,
    direction = 'horizontal',
    title = 'EMT-CAF\nscore',
    title.position = 'top',
    title.hjust = 0.5,
    barwidth = unit(70, 'pt'),
    barheight = unit(10, 'pt')
)

sil_vals_emt <- silhouette_vals(
    htmp_emt$dist_mat_genes,
    htmp_emt$hclust_genes,
    labels = htmp_emt$genes,
    apply_fun = NULL
)

sil_data <- as.data.table(
    merge(
        melt(
            sil_vals_emt$sil_vals,
            varnames = c('gene', 'height'),
            value.name = 'sil_val'
        ),
        melt(
            sil_vals_emt$cutree_clust_assign,
            varnames = c('gene', 'height'),
            value.name = 'clust_assign'
        )
    )
)

qplot(
    factor(gene, levels = gene),
    sil_val,
    data = sil_data[
        height == sort(unique(height), decreasing = TRUE)[2] # Use 3 if using GW_complete
        ][
            order(clust_assign, sil_val)
            ],
    geom = 'col'
) +
    facet_grid(
        cols = vars(clust_assign),
        scales = 'free',
        space = 'free'
    ) +
    geom_hline(
        yintercept = sil_data[
            height == sort(unique(height), decreasing = TRUE)[2],
            mean(sil_val)
            ]
    )

genes_filtered_sil <- sil_data[
    height == sort(unique(height), decreasing = TRUE)[2]
    ][
        order(-sil_val),
        head(as.character(gene), 50)
        ]

htmp_all <- deconv_heatmap(
    scores_data_transformed[
        genes_filtered_sil,
        c('gene', ..deconv_names)
        ]
    # order_genes_fun = 'seriate',
    # order_genes_method = 'GW_complete',
    # order_analyses_fun = 'seriate',
    # order_analyses_method = 'GW_complete'
)

deconv_heatmap_dendro_plot(
    htmp_all,
    direction = 'horizontal',
    title = 'EMT-CAF\nscore',
    title.position = 'top',
    title.hjust = 0.5,
    barwidth = unit(70, 'pt'),
    barheight = unit(10, 'pt')
)

# Filtering based on silhouette values for the x axis clusters:

dt <- scores_data_transformed[
    names(head(sort(apply(rank_mat, 1, quantile, 0.25)), 100)),
    c('gene', ..deconv_names)
    ]

cutree_analyses <- cutree(
    htmp_emt$hclust_analyses,
    k = 3
)

# The following gives, for each gene g, the average silhouette value of each of the 3 clusters
# when g is removed from the set of 100 genes.  So, a high value should mean that g adds
# noise to the cluster, because the cluster is stronger without g.

ave_sil_per_gene <- sapply(
    dt$gene,
    function(g) {
        
        dist_mat <- dist(t(dt[gene != g, -'gene']))
        
        sil <- cluster::silhouette(
            cutree_analyses,
            dist = dist_mat
        )
        
        sapply( # For USE.NAMES = TRUE to work, need to convert cluster number to character
            unique(as.character(sil[, 'cluster'])),
            function(i) {
                mean(sil[sil[, 'cluster'] == as.numeric(i), 'sil_width'])
            },
            USE.NAMES = TRUE
        )
        
    },
    USE.NAMES = TRUE
)

# High value for gene g means clusters are, on average, better off without g.  So when we
# sort, we want to take the head():

# genes_for_heatmap <- names(
#     head(
#         sort(
#             rowMeans(
#                 magrittr::set_rownames(
#                     apply(
#                         ave_sil_per_gene,
#                         1,
#                         order
#                     ),
#                     colnames(ave_sil_per_gene)
#                 )
#             )
#         ),
#         50
#     )
# )

# A "punishment" approach:

gene_ranks <- magrittr::set_rownames(
    apply(
        ave_sil_per_gene,
        1,
        order
    ),
    colnames(ave_sil_per_gene)
)

while(nrow(gene_ranks) > 50) {
    
    row_to_remove <- which.max(gene_ranks) %% nrow(gene_ranks)
    
    if(row_to_remove == 0) {row_to_remove <- nrow(gene_ranks)}
    
    gene_ranks <- gene_ranks[-row_to_remove, ]
    
}

genes_for_heatmap <- rownames(gene_ranks)

# Or, a "reward" approach:

gene_ranks <- magrittr::set_rownames(
    apply(
        ave_sil_per_gene,
        1,
        order
    ),
    colnames(ave_sil_per_gene)
)

genes_for_heatmap <- character()

while(length(genes_for_heatmap) < 50) {
    
    row_to_add <- which.min(gene_ranks) %% nrow(gene_ranks)
    
    if(row_to_add == 0) {row_to_add <- nrow(gene_ranks)}
    
    genes_for_heatmap <- c(
        genes_for_heatmap,
        rownames(gene_ranks)[row_to_add]
    )
    
    gene_ranks <- gene_ranks[-row_to_add, ]
    
}

# # Or, just take the highest scoring genes in each cluster:
# 
# temp <- lapply(
#     1:3,
#     function(i) {
#         dt[
#             ,
#             rowMeans(set_rownames(as.matrix(.SD), gene)),
#             .SDcols = names(cutree_analyses[cutree_analyses == i])
#         ]
#     }
# )
# 
# genes_for_heatmap <- unique(unlist(lapply(temp, function(x) names(x[x > 0.3]))))
# # Or
# genes_for_heatmap <- unique(unlist(lapply(temp, function(x) tail(names(sort(x)), 20))))
# 
# # It's a bit crude, but I think it actually works the best...

htmp_all <- deconv_heatmap(
    scores_data_transformed[
        genes_for_heatmap,
        c('gene', ..deconv_names)
        ],
    order_genes_fun = 'seriate',
    order_genes_method = 'GW_complete',
    order_analyses_fun = 'seriate',
    order_analyses_method = 'GW_complete'
)

deconv_heatmap_dendro_plot(
    htmp_all,
    direction = 'horizontal',
    title = 'EMT-CAF\nscore',
    title.position = 'top',
    title.hjust = 0.5,
    barwidth = unit(70, 'pt'),
    barheight = unit(10, 'pt')
)

# This doesn't work either...

# I tried filtering for score as well as just how common the genes are to the EMT
# component, but I'm not sure it adds much beyond just filtering for commonality.

sil_vals_emt <- silhouette_vals(
    htmp_emt$dist_mat_genes,
    htmp_emt$hclust_genes,
    labels = htmp_emt$genes,
    apply_fun = function(x) quantile(x, 0.25)
)

# I'm using the 75th percentile below for similar reasons to above.

htmp_caf <- deconv_heatmap(
    scores_data_transformed[
        names(tail(sort(apply(rank_mat, 1, quantile, 0.75)), 40)),
        c('gene', ..deconv_names)
        ]
)

sil_vals_caf <- silhouette_vals(
    htmp_caf$dist_mat_genes,
    htmp_caf$hclust_genes,
    labels = htmp_caf$genes,
    apply_fun = function(x) quantile(x, 0.25)
)

htmp_all <- deconv_heatmap(
    scores_data_transformed[
        unique(
            c(
                names(tail(sort(sil_vals_emt), 40)),
                names(tail(sort(sil_vals_caf), 25))
            )
        ),
        c('gene', ..deconv_names)
        ]
)

deconv_heatmap_dendro_plot(
    htmp_all,
    direction = 'horizontal',
    title = 'EMT-CAF\nscore',
    title.position = 'top',
    title.hjust = 0.5,
    barwidth = unit(70, 'pt'),
    barheight = unit(10, 'pt')
)

# Or...

scores <- rowMeans(
    -rank_mat + 1 + (as.matrix(scores_data_transformed[, -'gene']) + 1)/2
) # High score means EMT-associated; low means CAF-associated.

htmp_emt <- deconv_heatmap(
    scores_data_transformed[
        names(tail(sort(scores), 100)),
        c('gene', ..deconv_names)
        ]
)

sil_vals_emt <- silhouette_vals(
    htmp_emt$dist_mat_genes,
    htmp_emt$hclust_genes,
    labels = htmp_emt$genes
)

htmp_all <- deconv_heatmap(
    scores_data_transformed[
        names(tail(sort(sil_vals_emt), 50)),
        c('gene', ..deconv_names)
        ]
)

deconv_heatmap_dendro_plot(
    htmp_all,
    direction = 'horizontal',
    title = 'EMT-CAF\nscore',
    title.position = 'top',
    title.hjust = 0.5,
    barwidth = unit(70, 'pt'),
    barheight = unit(10, 'pt')
)

# Silhouette values to determine optimal number of clusters:

sil_vals_emt <- silhouette_vals(
    htmp_emt$dist_mat_genes,
    htmp_emt$hclust_genes,
    labels = htmp_emt$genes,
    apply_fun = NULL
)

sil_data <- as.data.table(
    merge(
        melt(
            sil_vals_emt$sil_vals,
            varnames = c('gene', 'height'),
            value.name = 'sil_val'
        ),
        melt(
            sil_vals_emt$cutree_clust_assign,
            varnames = c('gene', 'height'),
            value.name = 'clust_assign'
        )
    )
)

# I'm not sure this works very well either, but it does seem like 3 is the optimal
# number of clusters (try the following for other numbers besides 2 - if you increase
# the number you just end up with tiny clusters with zero silhouette value).

qplot(
    factor(gene, levels = gene),
    sil_val,
    data = sil_data[
        height == sort(unique(height), decreasing = TRUE)[2] # Gives 3 clusters
        ][
            order(clust_assign, sil_val)
            ],
    geom = 'col'
    # facets = sil_val ~ clust_assign
) +
    facet_grid(
        cols = vars(clust_assign),
        scales = 'free',
        space = 'free'
    ) +
    geom_hline(
        yintercept = sil_data[
            height == sort(unique(height), decreasing = TRUE)[2],
            mean(sil_val)
            ]
    )

genes_filtered_sil <- sil_data[
    height == sort(unique(height), decreasing = TRUE)[2] # Gives 3 clusters
    ][
        order(-sil_val),
        head(as.character(gene), 50)
        ]

htmp_all <- deconv_heatmap(
    scores_data_transformed[
        genes_filtered_sil,
        c('gene', ..deconv_names)
        ],
    order_genes_fun = 'seriate',
    order_genes_method = 'GW_complete'
    # order_analyses_fun = 'seriate',
    # order_analyses_method = 'GW_complete'
)

deconv_heatmap_dendro_plot(
    htmp_all,
    direction = 'horizontal',
    title = 'EMT-CAF\nscore',
    title.position = 'top',
    title.hjust = 0.5,
    barwidth = unit(70, 'pt'),
    barheight = unit(10, 'pt')
)

# I think this does work quite well, though you have to take the top 60 to keep SNAI2...

# Maybe I could remove genes that stand on their own even at a very high level, then
# find the optimum number of clusters, then filter again.

sil_data <- sil_data[
    gene %in% rownames(sil_vals_emt$sil_vals)[
        apply(sil_vals_emt$sil_vals, 1, function(x) sum(x == 0)) < 70
        ]
    ]

# Changing the number in the following doesn't always increase the number of clusters,
# since we've got rid of some genes.  But it looks like 4 clusters might now be the
# optimal number:

qplot(
    factor(gene, levels = gene),
    sil_val,
    data = sil_data[
        height == sort(unique(height), decreasing = TRUE)[4] # Gives 3 clusters
        ][
            order(clust_assign, sil_val)
            ],
    geom = 'col'
    # facets = sil_val ~ clust_assign
) +
    facet_grid(
        cols = vars(clust_assign),
        scales = 'free',
        space = 'free'
    ) +
    geom_hline(
        yintercept = sil_data[
            height == sort(unique(height), decreasing = TRUE)[4],
            mean(sil_val)
            ]
    )

genes_filtered_sil <- sil_data[
    height == sort(unique(height), decreasing = TRUE)[4] # Gives 3 clusters
    ][
        order(-sil_val),
        head(as.character(gene), 50)
        ]

htmp_all <- deconv_heatmap(
    scores_data_transformed[
        genes_filtered_sil,
        c('gene', ..deconv_names)
        ],
    order_genes_fun = 'seriate',
    order_genes_method = 'GW_complete'
    # order_analyses_fun = 'seriate',
    # order_analyses_method = 'GW_complete'
)

deconv_heatmap_dendro_plot(
    htmp_all,
    direction = 'horizontal',
    title = 'EMT-CAF\nscore',
    title.position = 'top',
    title.hjust = 0.5,
    barwidth = unit(70, 'pt'),
    barheight = unit(10, 'pt')
)

# Maybe it would be best to filter by silhouette values on a larger set of genes.





# Although I can't seem to improve on the following method, I'm not sure I like
# it that much, because I'm not sure taking the 25th percentile of the silhouette
# values is a very natural thing to do.

heatmaps <- setNames(
    lapply(
        list(scores_data, scores_data_transformed),
        function(dt) {
            
            setNames(
                lapply(
                    list(
                        names(deconv_data),
                        names(deconv_data)[
                            !(
                                names(deconv_data) %in% c(
                                    'Uveal Melanoma',
                                    'Mesothelioma',
                                    'Adrenocortical',
                                    'Melanoma Keratin'
                                )
                            )
                            ],
                        names(deconv_data)[
                            !(
                                names(deconv_data) %in% c(
                                    'Uveal Melanoma',
                                    'Mesothelioma',
                                    'Adrenocortical',
                                    'Melanoma Keratin',
                                    'Cholangiocarcinoma',
                                    'Melanoma MITF-low',
                                    'Melanoma Immune',
                                    'Thyroid'
                                )
                            )
                            ]
                    ),
                    function(deconv_names) {
                        
                        emt_caf <- deconv_heatmap(
                            dt[
                                unique(
                                    c(
                                        names(head(sort(apply(rank_mat, 1, quantile, 0.25)), 25)),
                                        names(tail(sort(apply(rank_mat, 1, quantile, 0.75)), 25)),
                                        'SNAI1',
                                        'SNAI2',
                                        'ZEB1',
                                        'ZEB2',
                                        'TWIST1',
                                        'VIM'
                                    )
                                ),
                                c('gene', ..deconv_names)
                                ]
                        )
                        
                        # Take top 75 most "common" EMT genes, make a heatmap, then calculate
                        # silhouette values and remake the heatmap with the top 50 most
                        # "clusterable" genes:
                        
                        emt <- deconv_heatmap(
                            dt[
                                names(head(sort(apply(rank_mat, 1, quantile, 0.25)), 75)),
                                c('gene', ..deconv_names)
                                ]
                        )
                        
                        sil_vals <- silhouette_vals(
                            emt$dist_mat_genes,
                            emt$hclust_genes,
                            labels = emt$genes
                        )
                        
                        emt <- deconv_heatmap(
                            dt[
                                unique(
                                    c(
                                        names(tail(sort(sil_vals), 50)),
                                        'SNAI1',
                                        'SNAI2',
                                        'ZEB1',
                                        'ZEB2',
                                        'TWIST1',
                                        'VIM'
                                    )
                                ),
                                c('gene', ..deconv_names)
                                ]
                        )
                        
                        # Now try silhouette values with EMT and CAF genes together:
                        
                        htmp_emt <- deconv_heatmap(
                            dt[
                                names(head(sort(apply(rank_mat, 1, quantile, 0.25)), 60)),
                                c('gene', ..deconv_names)
                                ]
                        )
                        
                        sil_vals_emt <- silhouette_vals(
                            htmp_emt$dist_mat_genes,
                            htmp_emt$hclust_genes,
                            labels = htmp_emt$genes,
                            apply_fun = function(x) quantile(x, 0.25)
                        )
                        
                        htmp_caf <- deconv_heatmap(
                            dt[
                                names(tail(sort(apply(rank_mat, 1, quantile, 0.75)), 40)),
                                c('gene', ..deconv_names)
                                ]
                        )
                        
                        sil_vals_caf <- silhouette_vals(
                            htmp_caf$dist_mat_genes,
                            htmp_caf$hclust_genes,
                            labels = htmp_caf$genes,
                            apply_fun = function(x) quantile(x, 0.25)
                        )
                        
                        htmp_all <- deconv_heatmap(
                            dt[
                                unique(
                                    c(
                                        names(tail(sort(sil_vals_emt), 40)),
                                        names(tail(sort(sil_vals_caf), 20)),
                                        'SNAI1',
                                        'SNAI2',
                                        'ZEB1',
                                        'ZEB2',
                                        'TWIST1',
                                        'VIM'
                                    )
                                ),
                                c('gene', ..deconv_names)
                                ]
                        )
                        
                        list(
                            
                            emt_caf = emt_caf,
                            emt = emt,
                            all = htmp_all
                            
                        )
                        
                    }
                ),
                c('unfiltered', 'light_filter', 'strong_filter')
            )
            
        }
    ),
    c('untransformed', 'transformed')
)





# Save plots in two separate PDFs, because we want to plot heatmap_commonality_emt with
# the dendrograms, which need more space:

pdf(
    '../data_and_figures/commonality_heatmap_emt_vs_caf.pdf',
    width = 9,
    height = 8.5
)

lapply(
    unlist(heatmaps, recursive = FALSE),
    function(li) li$emt_caf$heatmap
)

dev.off()

pdf(
    '../data_and_figures/commonality_heatmap_emt.pdf',
    width = 9.5,
    height = 10
)

# Plot the heatmap with dendrograms (the additional arguments go to guide_colourbar()):

lapply(
    unlist(heatmaps, recursive = FALSE),
    function(li) {
        lapply(
            c('emt', 'all'),
            function(type) {
                deconv_heatmap_dendro_plot(
                    li[[type]],
                    direction = 'horizontal',
                    title = 'EMT-CAF score',
                    title.position = 'top',
                    title.hjust = 0.5
                )
            }
        )
    }
)

dev.off()

# Optimally-sized version of my favourite heatmap:

pdf(
    '../data_and_figures/commonality_heatmap_all.pdf',
    width = 7.5,
    height = 10
)

deconv_heatmap_dendro_plot(
    heatmaps$transformed$strong_filter$all,
    direction = 'horizontal',
    title = 'EMT-CAF\nscore',
    title.position = 'top',
    title.hjust = 0.5,
    barwidth = unit(70, 'pt'),
    barheight = unit(10, 'pt')
)

dev.off()

# The same version without the dendrograms:

pdf(
    '../data_and_figures/commonality_heatmap_all_simple.pdf',
    width = 7,
    height = 8.5
)

heatmaps$transformed$strong_filter$all$heatmap

dev.off()
