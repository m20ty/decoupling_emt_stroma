# The following tries to distinguish rare from generic EMT.  I think it also works quite
# well (possibly better, though it's hard to define what is better) when using the 99th
# percentile instead of the 95th, or indeed the variance, though for this I think our
# loess parameters might not be so appropriate - a linear model is better, though I
# think quantile regression using quantreg::rq(y ~ x, tau = 0.5) works best, perhaps
# using different quantiles (i.e. values for tau) for each cancer type.

# We could also try with weighed least squares regression, which you can do with the lm()
# function, using the weights argument (there's also MASS::rlm(), 'robust linear regression').
# For the weights, we could use the correlations of the genes with the number of genes detected,
# so we give less weight to genes whose variation seems to come primarily from the variation
# in the complexity.

# I could also try replacing the loess curve with a linear model (perhaps using quantile
# regression) which is fitted using only the values where x is greater than the chosen
# threshold - that is, fit a linear model only to the right of the dashed vertical line.

for(ct in names(sc_groups_data)) {
    setkey(sc_groups_data[[ct]]$data, id)
}

# Rare and generic EMT genes using variance against mean (note using robustbase::Qn instead
# of variance produces nearly identical results):

rare_emt_genes_params <- list(
    hnsc = list(
        tau = 0.75,
        title = 'Head and Neck'
    ),
    lung = list(
        tau = 0.85,
        x_cond_fun = function(v) {v > quantile(v, 0.2)},
        title = 'Lung'
    ),
    brca = list(
        tau = 0.8,
        title = 'Breast'
    ),
    coadread = list(
        tau = 0.75,
        title = 'Colorectal'
    ),
    paad = list(
        tau = 0.65,
        x_cond_fun = function(v) {v > quantile(v, 0.15)},
        title = 'Pancreas'
    )
)

generic_emt_genes_params <- list(
    hnsc = list(
        tau = 0.35,
        y_mod_cond_fun = `<`,
        title = 'Head and Neck'
    ),
    lung = list(
        tau = 0.45,
        y_mod_cond_fun = `<`,
        x_cond_fun = function(v) {v > quantile(v, 0.05)},
        title = 'Lung'
    ),
    brca = list(
        tau = 0.5,
        y_mod_cond_fun = `<`,
        title = 'Breast'
    ),
    coadread = list(
        tau = 0.4,
        y_mod_cond_fun = `<`,
        title = 'Colorectal'
    ),
    paad = list(
        tau = 0.35,
        y_mod_cond_fun = `<`,
        title = 'Pancreas'
    )
)

rare_vs_generic_emt_genes <- setNames(
    lapply(
        list(
            rare_emt_genes_params,
            generic_emt_genes_params
        ),
        function(params) {
            
            sapply(
                names(sc_groups_data),
                function(ct) {
                    
                    gene_stat_data <- do.call(
                        gene_stat,
                        args = c(
                            list(sc_groups_list = sc_groups_data[[ct]]),
                            list(
                                mod_fun_args = list(
                                    formula = y ~ x - 1,
                                    tau = params[[ct]]$tau
                                )
                            )
                        )
                    )
                    
                    list(
                        gene_stat_data = gene_stat_data,
                        genes = do.call(
                            rare_genes,
                            args = c(
                                list(gene_stat_data = gene_stat_data),
                                params[[ct]][-1]
                            )
                        )
                    )
                    
                },
                simplify = FALSE,
                USE.NAMES = TRUE
            )
            
        }
    ),
    c('rare', 'generic')
)





# Using 95th percentile against mean:

rare_emt_genes_params <- list(
    hnsc = list(
        x_cond_fun = function(v) {v > 0.5},
        y_cond_fun = function(v) {v > 0},
        title = 'Head and Neck'
    ),
    lung = list(
        x_cond_fun = function(v) {v > 0.5},
        y_cond_fun = function(v) {v > 0},
        title = 'Lung'
    ),
    brca = list(
        x_cond_fun = function(v) {v > 0.75},
        y_cond_fun = function(v) {v > 0},
        title = 'Breast'
    ),
    coadread = list(
        x_cond_fun = function(v) {v > 0.5},
        y_cond_fun = function(v) {v > 0},
        title = 'Colorectal'
    ),
    paad = list(
        x_cond_fun = function(v) {v > 0.125},
        y_cond_fun = function(v) {v > 0},
        title = 'Pancreas'
    )
)

generic_emt_genes_params <- sapply(
    rare_emt_genes_params,
    c,
    list(y_mod_cond_fun = `<`),
    simplify = FALSE,
    USE.NAMES = TRUE
)

rare_vs_generic_emt_genes <- setNames(
    lapply(
        list(
            rare_emt_genes_params,
            generic_emt_genes_params
        ),
        function(params) {
            
            sapply(
                names(sc_groups_data),
                function(ct) {
                    
                    gene_stat_data <- do.call(
                        gene_stat,
                        args = c(
                            list(
                                sc_groups_list = sc_groups_data[[ct]],
                                y_fun = function(v) quantile(v, 0.95),
                                mod_fun = loess,
                                cond_fun = function(v, w) {w > 0}
                            ),
                            list(
                                mod_fun_args = list(
                                    formula = y ~ x,
                                    span = 0.85,
                                    degree = 2,
                                    family = 'symmetric'
                                )
                            )
                        )
                    )
                    
                    list(
                        gene_stat_data = gene_stat_data,
                        genes = do.call(
                            rare_genes,
                            args = c(
                                list(gene_stat_data = gene_stat_data),
                                params[[ct]]
                            )
                        )
                    )
                    
                },
                simplify = FALSE,
                USE.NAMES = TRUE
            )
            
        }
    ),
    c('rare', 'generic')
)





# Heatmap showing the genes for all cancer types:

membership_heatmaps <- sapply(
    names(rare_vs_generic_emt_genes),
    function(w) {
        
        genes_table <- dcast(
            rbindlist(
                lapply(
                    names(rare_vs_generic_emt_genes[[w]]),
                    function(ct) {
                        data.table(
                            analysis = ct,
                            gene = rare_vs_generic_emt_genes[[w]][[ct]]$genes$genes,
                            membership = 1
                        )
                    }
                )
            ),
            analysis ~ gene,
            value.var = 'membership',
            fill = 0
        )
        
        # The following restricts to just the genes that occur in more than one cancer type.
        
        genes_table[
            ,
            names(genes_table[, -'analysis']) := lapply(
                .SD,
                function(x) {switch((sum(x) == 1) + 1, x, NULL)}
            ),
            .SDcols = -'analysis'
            ]
        
        ordering_analyses <- hclust(dist(genes_table[, -'analysis']))$order
        ordering_genes <- hclust(dist(t(genes_table[, -'analysis'])))$order
        
        genes_table_melted <- melt(
            genes_table,
            id.vars = 'analysis',
            variable.name = 'gene',
            value.name = 'membership'
        )[
            ,
            analysis_name := plyr::mapvalues(
                analysis,
                c('hnsc', 'lung', 'brca', 'coadread', 'paad'),
                c('Head and Neck', 'Lung', 'Breast', 'Colorectal', 'Pancreas')
            )
            ]
        
        htmp <- ggplot(
            genes_table_melted,
            aes(
                x = factor(analysis_name, levels = unique(analysis_name)[ordering_analyses]),
                y = factor(gene, levels = unique(gene)[ordering_genes])
            )
        ) +
            geom_raster(
                aes(
                    fill = as.factor(membership),
                    alpha = plyr::mapvalues(membership, c(0, 1), c(0, 1)))
            ) +
            scale_x_discrete(expand = c(0, 0)) +
            scale_y_discrete(expand = c(0, 0)) +
            scale_fill_manual(values = c('white', 'navy')) +
            labs(x = 'Cancer type', y = 'Gene') +
            theme_bw() +
            theme(
                legend.position = 'none',
                axis.text.x = element_text(angle = 55, hjust = 1),
                axis.text.y = element_text(size = 8)
                # plot.margin = unit(c(2, 30, 2, 30), 'mm')
            ) +
            labs(title = paste(stringr::str_to_title(w), 'EMT markers'))
        
        list(
            heatmap = htmp,
            data = genes_table,
            melted_data = genes_table_melted,
            ordering_analyses = ordering_analyses,
            ordering_genes = ordering_genes
        )
        
    },
    simplify = FALSE,
    USE.NAMES = TRUE
)





# Save figures to PDF (this is for variance against mean; I similarly saved the analysis
# with 95th percentile against mean under file name 'rare_vs_generic_emt_95per.pdf'):

pdf('../data_and_figures/rare_vs_generic_emt_var.pdf', width = 12, height = 8)

egg::ggarrange(
    plots = lapply(
        rare_vs_generic_emt_genes$rare,
        function(li) {li$genes$plot}
    ),
    ncol = 3,
    nrow = 2,
    bottom = 'Mean',
    left = 'Variance',
    top = 'Rare EMT genes'
)

egg::ggarrange(
    plots = lapply(
        rare_vs_generic_emt_genes$generic,
        function(li) {li$genes$plot}
    ),
    ncol = 3,
    nrow = 2,
    bottom = 'Mean',
    left = 'Variance',
    top = 'Generic EMT genes'
)

egg::ggarrange(
    plots = lapply(membership_heatmaps, `[[`, 'heatmap'),
    ncol = 2,
    nrow = 1
)

dev.off()





# Try to widen rare and generic EMT programs:

expanded <- sapply(
    names(single_cell_metadata),
    function(ct) {
        
        sc_data <- eval(single_cell_metadata[[ct]]$read_quote)
        
        sc_data[
            cell_type == 'cancer',
            c(
                list(gene = names(.SD)),
                sapply(
                    names(rare_vs_generic_emt_genes),
                    function(w) {
                        rowMeans(
                            cor(
                                .SD,
                                .SD[
                                    ,
                                    rare_vs_generic_emt_genes[[w]][[ct]]$genes$genes,
                                    with = FALSE
                                    ]
                            )
                        )
                    },
                    simplify = FALSE,
                    USE.NAMES = TRUE
                )
            ),
            .SDcols = -switch(
                (ct == 'paad') + 1,
                c('id', 'patient', 'cell_type'),
                c('id', 'cell_type')
            )
            ][
                ,
                rare_generic_score := rare - generic
                ]
        
    },
    simplify = FALSE,
    USE.NAMES = TRUE
)

expanded_list <- sapply(
    expanded,
    function(li) {
        list(
            rare = li[rare_generic_score > 0.08, gene],
            generic = li[rare_generic_score < -0.08, gene]
        )
    },
    simplify = FALSE,
    USE.NAMES = TRUE
)

# There are no genes common to all the cancer types:

# Reduce(intersect, rare_emt_genes_expanded_list) # character(0)

# But let's take the genes which are common to at least two cancer types:

expanded_list <- sapply





# We can also try just scoring by variance over mean (not sure how accurate it is, but
# probably at least as accurate as the above):

rare_vs_generic_emt_scores <- sapply(
    sc_groups_data,
    function(li) {
        li$data[
            li$cells_filtered
            ][
                cell_type == 'cancer',
                .(
                    gene = names(.SD),
                    score = sapply(.SD, function(x) {var(x)/mean(x)})
                ),
                .SDcols = li$genes_filtered
                ]
    },
    simplify = FALSE,
    USE.NAMES = TRUE
)

all_rare_emt_genes <- sort(
    unique(
        unlist(
            lapply(
                rare_vs_generic_emt_scores,
                function(li) {
                    li[order(-score)][1:50, gene]
                }
            )
        )
    )
)

heatmap_data <- rbindlist(
    lapply(
        names(single_cell_metadata),
        function(ct) {
            
            sc_data <- eval(single_cell_metadata[[ct]]$read_quote)
            
            sc_data[
                cell_type == 'cancer',
                .(
                    gene = all_rare_emt_genes[all_rare_emt_genes %in% names(sc_data)],
                    analysis = ct,
                    score = sapply(.SD, function(x) {var(x)/mean(x)})
                ),
                .SDcols = all_rare_emt_genes[all_rare_emt_genes %in% names(sc_data)]
                ]
            
        }
    )
)

heatmap_data <- heatmap_data[complete.cases(heatmap_data)]

heatmap_data[
    ,
    normalised_score := score/quantile(score, 0.9),
    by = analysis
    ]

setkey(heatmap_data, gene, analysis)

top_genes <- sort(heatmap_data[, .(mean_score = mean(normalised_score)), by = gene][order(-mean_score)][1:50, gene])

cast_data <- dcast(heatmap_data, gene ~ analysis, value.var = 'normalised_score')

ordering_top_genes <- hclust(dist(cast_data[gene %in% top_genes, -'gene']), method = 'average')$order
ordering_analyses_top_genes <- hclust(dist(t(cast_data[gene %in% top_genes, -'gene'])), method = 'average')$order

rare_emt_genes_commonality <- ggplot(
    heatmap_data[top_genes],
    aes(
        x = factor(analysis, levels = unique(analysis)[ordering_analyses_top_genes]),
        y = factor(gene, levels = unique(gene)[ordering_top_genes])
    )
) +
    geom_raster(aes(fill = normalised_score)) +
    scale_fill_gradientn(
        colours = rev(colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))(50))
    ) +
    scale_x_discrete(
        expand = c(0, 0)
    ) +
    scale_y_discrete(
        expand = c(0, 0)
    ) +
    theme(
        axis.ticks = element_blank(),
        axis.text.x = element_text(
            angle = 55,
            hjust = 1
        ),
        axis.title.y = element_blank(),
        axis.ticks.length = unit(0, 'pt'),
        plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), 'pt')
    )

# I'm not sure this is useful...  I think I should at least have something bi-directional,
# perhaps with rare and generic EMT in the same heatmap.





# Attempt at biclustering for finding rare EMT genes:

sc_data <- eval(single_cell_metadata$hnsc$read_quote)

setkey(sc_data, id)

dat_mat <- as.matrix(
    sc_data[
        sc_groups_data$hnsc$cells_filtered
        ][
            cell_type == 'cancer',
            sc_groups_data$hnsc$genes_filtered,
            with = FALSE
            ]
)

biclust_result <- biclust::biclust(
    dat_mat,
    method = biclust::BCCC()
)

# Try the following also with number = 2 and number = 3:

biclust::drawHeatmap(
    dat_mat,
    bicResult = biclust_result,
    number = 1,
    local = FALSE
)

# Look at the genes in each cluster:

colnames(dat_mat)[biclust_result@NumberxCol[1, ]] # Replace 1 with 2 and 3

# Looks like none of the clusters picked up VIM!  But there are parameters in the BCCC method that
# we can tweak.  We could also try some other method: BCXmotifs, BCPlaid and BCSpectral seem like
# they should work, but I don't have a good reason to pick any one over another.  But before I
# delve into this too deeply, let's just try separate clustering of the cells and genes.

sc_data_cancer <- sc_data[
    sc_groups_data$hnsc$cells_filtered
    ][
        cell_type == 'cancer',
        c('id', sc_groups_data$hnsc$genes_filtered),
        with = FALSE
        ]

ordering_cells <- seriation::seriate(dist(sc_data_cancer[, -'id']), method = 'GW')[[1]]$order
ordering_genes <- seriation::seriate(dist(t(sc_data_cancer[, -'id'])), method = 'GW')[[1]]$order

ggplot(
    melt(
        sc_data_cancer,
        id.vars = 'id',
        variable.name = 'gene',
        value.name = 'expression_level'
    )
) +
    geom_raster(
        aes(
            x = factor(id, levels = unique(id)[ordering_cells]),
            y = factor(gene, levels = unique(gene)[ordering_genes]),
            fill = expression_level
        )
    ) +
    scale_fill_gradientn(
        colours = rev(colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))(50))
    )

# I'm not sure this works anyway.  It seems the problem is that there is no distinct cluster of
# "rare EMT cells" - I believe that rare EMT is a thing, but it is not seen uniquely in one set
# of cancer cells.  I think some of the genes are expressed in some cancer cells, others in
# others etc.  The other problem is that if we do clustering on the cells then most of the
# clustering will be into patients, so the only way to recover EMT via clustering is to cluster
# per patient and then cluster the resulting programs, as in the HNSCC paper.  But this is a
# technique that won't be available to us in bulk expression data, and also may not detect EMT
# in other cancer types where EMT is even rarer.  Thinking about it this way, the correlation
# with EMT score may still be the best method.  I suppose it might be worth applying the
# metaprograms clustering method to the other scRNA-seq datasets, to see if you find some EMT:
# then you can look at what separates rare from generic EMT in these cancer types, to hopefully
# find something that we can use in the bulk data.





# The following was some playing around I did with the Z score heatmaps/data:

# In the following, we could use an alternative function in the sapply loop, to return
# the average nonzero pairwise difference:

# function(x) {
#     temp <- as.numeric(dist(x))
#     mean(temp[temp > 0])
# }

# I'm not sure it's any more useful than the mean of the nonzero values, though.

sc_cancer_z_score[[ct]]$data[
    ,
    sort(sapply(.SD, function(x) mean(x[x > 0])), decreasing = TRUE),
    .SDcols = -c('id', 'cell_type', 'genes_detected', 'patient')
    ]

rare_emt_genes <- sapply(
    names(sc_cancer),
    function(ct) {
        
        rare_shared_emt[[ct]]$score_cor_data[
            sc_cancer_z_score[[ct]]$genes_filtered
            ][
                emt_cor >= rare_shared_emt[[ct]]$prob_funs$emt_cor$quantile_fun(0.975) &
                    sc_cancer_z_score[[ct]]$data[
                        ,
                        sapply(.SD, function(x) mean(x[x > 0])),
                        .SDcols = gene
                        ] >= 0.8
                ][
                    order(-emt_cor),
                    gene
                    ]
        
    },
    simplify = FALSE,
    USE.NAMES = TRUE
)

z_score_data_shuffled <- sapply(
    sc_cancer,
    function(li) {
        copy(li$data)[
            ,
            (li$genes_filtered) := lapply(
                transpose(transpose(.SD)[, lapply(.SD, sample)]),
                function(x) {(x - mean(x))/sd(x)}
            ),
            .SDcols = li$genes_filtered
            ]
    },
    simplify = FALSE,
    USE.NAMES = TRUE
)

z_score_data_shuffled[[ct]][
    ,
    mean(sapply(.SD, function(x) mean(x[x > 0])))
    ]

sc_cancer_z_score$hnsc$data[
    sc_cancer_z_score$hnsc$cells_filtered
    ][
        sc_cancer_z_score$hnsc$ordering_cells$cancer,
        sort(sapply(.SD, function(x) {var(caTools::runmean(x, 10))}), decreasing = TRUE),
        .SDcols = sc_cancer_z_score$hnsc$genes_filtered
        ]

sc_cancer_z_score[[ct]]$data[
    sc_cancer_z_score$hnsc$cells_filtered
    ][
        sc_cancer_z_score$hnsc$ordering_cells$cancer,
        sort(sapply(.SD, function(x) {sum((1:.N)*x/.N)})), # Could replace *x with *(x^2)
        .SDcols = -c('id', 'cell_type', 'genes_detected', 'patient')
        ]

# I guess I could also try using a t test, or similar, to see if the mean of the true
# distribution differs significantly from the shuffled distribution.  When I say "distribution"
# here, I mean the distribution defined by the PDF determined by the ordered gene vector (i.e.
# the vector of Z scores).

sort(
    sapply(
        sc_cancer[[ct]]$genes_filtered,
        function(g) {
            wilcox.test(
                sc_cancer_z_score[[ct]]$data[
                    sc_cancer_z_score[[ct]]$cells_filtered
                    ][
                        sc_cancer_z_score[[ct]]$ordering_cells$cancer,
                        sample(
                            1:.N,
                            1000,
                            replace = TRUE,
                            prob = get(g) - min(get(g))
                        )
                        ],
                z_score_data_shuffled[[ct]][
                    sc_cancer_z_score[[ct]]$cells_filtered
                    ][
                        sc_cancer_z_score[[ct]]$ordering_cells$cancer,
                        sample(
                            1:.N,
                            1000,
                            replace = TRUE,
                            prob = get(g) - min(get(g))
                        )
                        ]
            )$p.value
        }
    )
)

# This works OK, but still doesn't filter out some genes that I want to filter out, like
# CD44 (though it would depend on your significance threshold).
