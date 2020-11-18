# EMT minus CAF scores (but I'm not sure I want to do this):

# data_for_test <- sapply(
#     deconv_data,
#     function(li) {
#         expression_data[
#             li$sample_ids,
#             .SD[, head(li$genes_filtered[li$ordering], 20), with = FALSE] - rowMeans(
#                 .SD[, tail(li$genes_filtered[li$ordering], 20), with = FALSE]
#             )
#         ]
#     },
#     simplify = FALSE,
#     USE.NAMES = TRUE
# )

# Could do something similar with epithelial markers:

epi_genes <- c(
    'CDH1',
    'EPCAM',
    'SFN',
    names(expression_data)[
        grep('^KRT[0-9]+$', names(expression_data))
        ]
)

keratins <- names(expression_data)[
    grep('^KRT[0-9]+$', names(expression_data))
    ]

epi_genes_list <- sapply(
    deconv_names,
    function(ct) {
        c(
            'CDH1',
            'EPCAM',
            'SFN',
            expression_data[
                deconv_data[[ct]]$sample_ids,
                .(
                    krt = keratins,
                    av = colMeans(.SD)
                ),
                .SDcols = keratins
                ][
                    av >= 3,
                    krt
                    ]
        )
    },
    simplify = FALSE,
    USE.NAMES = TRUE
)

# data_for_test <- sapply(
#     deconv_data[deconv_names],
#     function(li) {
#         expression_data[
#             li$sample_ids,
#             c(
#                 .(id = id),
#                 .SD[, head(li$genes_filtered[li$ordering], 20), with = FALSE] -
#                     rowMeans(.SD[, ..epi_genes])
#             )
#         ]
#     },
#     simplify = FALSE,
#     USE.NAMES = TRUE
# )

data_for_clin_tests <- setNames(
    lapply(
        list(head, tail),
        function(FUN) {
            sapply(
                deconv_names,
                function(ct) {
                    expression_data[
                        deconv_data[[ct]]$sample_ids,
                        c(
                            .(id = id),
                            .SD[
                                ,
                                FUN(deconv_data[[ct]]$genes_filtered[deconv_data[[ct]]$ordering], 20),
                                with = FALSE
                                ] - rowMeans(
                                    .SD[, epi_genes_list[[ct]], with = FALSE]
                                )
                            # ] - rowMeans(.SD[, ..epi_genes])
                        )
                        ]
                },
                simplify = FALSE,
                USE.NAMES = TRUE
            )
        }
    ),
    c('deconv_emt', 'deconv_caf')
)

clin_cor <- sapply(
    names(emt_types),
    function(emt_type) {
        
        clinical_test(
            expression_data,
            # data_for_clin_tests[[emt_type]],
            lapply(deconv_data[emt_types[[emt_type]]], `[[`, 'sample_ids'),
            clin_cor_genes[[emt_type]],
            clinical_data,
            clin_var = c(
                'number_of_lymphnodes_positive_by_he',
                'followup_treatment_success',
                'days_to_death',
                'lymphovascular_invasion',
                'pathologic_t',
                'pathologic_n',
                'pathologic_m',
                'neoplasm_histologic_grade'
            ),
            wilcox_test_x_expr = list(
                quote(variable > 1),
                quote(variable != 'complete remission/response'),
                quote(variable < quantile(variable, 0.4)),
                quote(variable %in% c('yes', 'present')),
                quote(
                    startsWith(variable, 't3') |
                        startsWith(variable, 't4')
                ),
                quote( # Maybe n1 should be included in here as well...
                    startsWith(variable, 'n2') |
                        startsWith(variable, 'n3')
                ),
                quote(startsWith(variable, 'm1')),
                quote(variable %in% c('gb', 'g1', 'g2', 'low grade'))
            ),
            wilcox_test_y_expr = list(
                quote(variable == 0),
                quote(variable == 'complete remission/response'),
                quote(variable > quantile(variable, 0.6)),
                quote(variable %in% c('no', 'absent')),
                quote(
                    startsWith(variable, 't0') |
                        startsWith(variable, 't1') |
                        startsWith(variable, 't2')
                ),
                quote(
                    startsWith(variable, 'n0') |
                        startsWith(variable, 'n1')
                ),
                quote(
                    startsWith(variable, 'm0') |
                        startsWith(variable, 'cm0')
                ),
                quote(variable %in% c('g3', 'g4', 'high grade'))
            ),
            min_samples = 10
            # genes_for_data_transform = clin_cor_genes[[
            #     names(emt_types)[names(emt_types) != emt_type]
            # ]]
        )
        
    },
    simplify = FALSE,
    USE.NAMES = TRUE
)





# I tried to do something clever here to tidy up my code - I tried storing arguments for the
# cor_with_clinical() function in a list, so I could apply this function over the list using
# do.call.  The problem is that the wilcox_test_x_expr and wilcox_test_y_expr arguments have
# to be quoted expressions, which leads to an error in eval() inside the cor_with_clinical()
# function, I think because the expression is not quoted in the same environment that the
# data_for_test data table is created (which is where we look for 'variable').  I can't store
# an unquoted condition (like variable > 0) in a list, but if I quote it the function won't
# look for 'variable' in the right place.  So I had to give up on this idea.

# Read more about this problem here, especially in the 'Calling from another function'
# section:

# http://adv-r.had.co.nz/Computing-on-the-language.html

clin_cor <- list(
    
    lymph_node_metastasis_positive = cor_with_clinical(
        expression_data,
        ids_genes_ordering_list,
        clinical_data,
        clin_var = 'number_of_lymphnodes_positive_by_he',
        wilcox_test_x_expr = variable > 0,
        wilcox_test_y_expr = variable == 0,
        plot_title = 'Correlation with lymph node metastasis'
    ),
    
    lymph_node_metastasis_multiple = cor_with_clinical(
        expression_data,
        ids_genes_ordering_list,
        clinical_data,
        clin_var = 'number_of_lymphnodes_positive_by_he',
        wilcox_test_x_expr = variable > 1,
        wilcox_test_y_expr = variable == 0,
        plot_title = 'Correlation with lymph node metastasis (0 vs. >1)'
    ),
    
    therapy_resistance = cor_with_clinical(
        expression_data,
        ids_genes_ordering_list,
        clinical_data,
        clin_var = 'followup_treatment_success',
        wilcox_test_x_expr = variable != 'complete remission/response',
        wilcox_test_y_expr = variable == 'complete remission/response',
        plot_title = 'Correlation with resistance to therapy'
    ),
    
    mortality = cor_with_clinical(
        expression_data,
        ids_genes_ordering_list,
        clinical_data,
        clin_var = 'days_to_death',
        wilcox_test_x_expr = variable < median(variable),
        wilcox_test_y_expr = variable > median(variable),
        plot_title = 'Correlation with mortality (median)'
    ),
    
    mortality_strict = cor_with_clinical(
        expression_data,
        ids_genes_ordering_list,
        clinical_data,
        clin_var = 'days_to_death',
        wilcox_test_x_expr = variable < quantile(variable, 0.4),
        wilcox_test_y_expr = variable > quantile(variable, 0.6),
        plot_title = 'Correlation with mortality (40th vs. 60th percentile)'
    ),
    
    lymphovascular_invasion = cor_with_clinical(
        expression_data,
        ids_genes_ordering_list,
        clinical_data,
        clin_var = 'lymphovascular_invasion',
        wilcox_test_x_expr = variable %in% c('yes', 'present'),
        wilcox_test_y_expr = variable %in% c('no', 'absent'),
        plot_title = 'Correlation with lymphovascular invasion'
    ),
    
    t_stage = cor_with_clinical(
        expression_data,
        ids_genes_ordering_list,
        clinical_data,
        clin_var = 'pathologic_t',
        wilcox_test_x_expr = startsWith(variable, 't3') |
            startsWith(variable, 't4'),
        wilcox_test_y_expr = startsWith(variable, 't0') |
            startsWith(variable, 't1') |
            startsWith(variable, 't2'),
        plot_title = 'Correlation with T stage',
        amatch_max_dist = 2
    ),
    
    n_stage = cor_with_clinical(
        expression_data,
        ids_genes_ordering_list,
        clinical_data,
        clin_var = 'pathologic_n',
        wilcox_test_x_expr = startsWith(variable, 'n2') |
            startsWith(variable, 'n3'),
        wilcox_test_y_expr = startsWith(variable, 'n0') |
            startsWith(variable, 'n1'),
        plot_title = 'Correlation with N stage',
        amatch_max_dist = 2
    ),
    
    m_stage = cor_with_clinical(
        expression_data,
        ids_genes_ordering_list,
        clinical_data,
        clin_var = 'pathologic_m',
        wilcox_test_x_expr = startsWith(variable, 'm1'),
        wilcox_test_y_expr = startsWith(variable, 'm0') |
            startsWith(variable, 'cm0'),
        plot_title = 'Correlation with M stage',
        amatch_max_dist = 2
    ),
    
    grade = cor_with_clinical(
        expression_data,
        ids_genes_ordering_list,
        clinical_data,
        clin_var = 'neoplasm_histologic_grade',
        wilcox_test_x_expr = variable %in% c('gb', 'g1', 'g2', 'low grade'),
        wilcox_test_y_expr = variable %in% c('g3', 'g4', 'high grade'),
        plot_title = 'Correlation with tumour grade'
    )
    
)

pdf('../data_and_figures/clinical_analyses.pdf', width = 35, height = 14)

egg::ggarrange(
    plots = lapply(
        clin_cor,
        function(li) li$plot
    ),
    nrow = 2,
    ncol = 5
)

dev.off()

# To find other partial variable names we can use:

# clinical_data[
#     ,
#     .(
#         non_na_names = names(
#             .SD[
#                 ,
#                 sapply(
#                     .SD,
#                     function(x) {
#                         switch(
#                             (sum(is.na(x) | x == '') == .N) + 1,
#                             x,
#                             NULL
#                         )
#                     }
#                 )
#             ]
#         )
#     ),
#     by = cancer_type
# ][
#     ,
#     .(
#         closest_match = non_na_names[
#             stringdist::amatch(
#                 'lymphovascular_invasion',
#                 non_na_names,
#                 maxDist = 10
#             )
#         ]
#     ),
#     by = cancer_type
# ]





# Look for more genes that correlate with our EMT genes, using absolute expression levels:

genes_list <- extra_genes_by_cor(
    expression_data,
    ids_genes_ordering_list,
    n = 1000,
    transform_data = FALSE
)

# Use these genes for transforming the gene expression space, then look again for genes
# that correlate with EMT genes in this transformed space (we do this because it takes
# too long to transform the whole space, i.e. with all the genes):

genes_list_transformed <- extra_genes_by_cor(
    expression_data,
    ids_genes_ordering_list,
    additional_genes_list = genes_list,
    transform_data = TRUE
)

# Make commonality heatmaps using correlation:

heatmap_commonality_genes_list <- commonality_heatmap(
    expression_data,
    ids_genes_ordering_list,
    extra_ordered_genes_list = genes_list,
    two_sided = FALSE,
    n_head = 60,
    legend_title = 'Correlation'
)

heatmap_commonality_genes_list_transformed <- commonality_heatmap(
    expression_data,
    ids_genes_ordering_list,
    extra_ordered_genes_list = genes_list_transformed,
    two_sided = FALSE,
    n_head = 60,
    legend_title = 'Correlation'
)

# We can try with expression levels, though I'm not sure it helps:

heatmap_commonality_genes_list_exp <- commonality_heatmap(
    expression_data,
    ids_genes_ordering_list,
    extra_ordered_genes_list = genes_list,
    two_sided = FALSE,
    score_type = 'expression_level',
    transform_data = FALSE,
    normalise_scores = FALSE,
    n_head = 60,
    legend_title = 'Expression\nlevel'
)

heatmap_commonality_genes_list_transformed_exp <- commonality_heatmap(
    expression_data,
    ids_genes_ordering_list,
    extra_ordered_genes_list = genes_list_transformed,
    two_sided = FALSE,
    score_type = 'expression_level',
    transform_data = FALSE,
    normalise_scores = FALSE,
    n_head = 60,
    legend_title = 'Expression\nlevel'
)

# Save to PDF:

pdf('../data_and_figures/extra_genes_heatmaps.pdf', width = 10, height = 10)

heatmap_commonality_genes_list$heatmap
heatmap_commonality_genes_list_transformed$heatmap
heatmap_commonality_genes_list_exp$heatmap
heatmap_commonality_genes_list_transformed_exp$heatmap

dev.off()





# Look for genes correlating highly with epithelial marker genes:

epi_markers <- c(
    'CDH1',
    'EPCAM',
    'SFN',
    names(expression_data)[grep('^KRT[0-9]|^KRTD', names(expression_data))]
)

epi_markers <- sapply(
    
    ids_genes_ordering_list,
    
    function(li) {
        
        epi_levels <- expression_data[
            li$sample_ids,
            colMeans(.SD),
            .SDcols = epi_markers
            ]
        
        # I want to make sure I have at least E-cadherin and EPCAM, for cases where the
        # expression levels are extremely low:
        
        li_epi_markers <- unique(
            c(
                'CDH1',
                'EPCAM',
                names(epi_levels[epi_levels > 6])
            )
        )
        
        li_epi_cors <- expression_data[
            li$sample_ids,
            cor(.SD, rowMeans(.SD[, ..li_epi_markers])),
            .SDcols = -'id'
            ]
        
        li_epi_markers <- sort(
            unique(
                c(
                    li_epi_markers,
                    # names(li_epi_cors[order(-li_epi_cors[, 1]), ])[1:20]
                    names(li_epi_cors[order(-li_epi_cors[, 1]), ])[
                        li_epi_cors[order(-li_epi_cors[, 1]), ] > 0.75
                        ]
                )
            )
        )
        
        # The following finds, instead of correlation with the average, the average of the
        # individual correlations.  But I found it gave way too many genes in some cases,
        # especially Pancreas.
        
        # li_epi_cors <- expression_data[
        #     li$sample_ids,
        #     rowMeans(cor(.SD, .SD[, ..li_epi_markers])),
        #     .SDcols = -'id'
        # ]
        
        # li_epi_markers <- sort(
        #     unique(
        #         c(
        #             li_epi_markers,
        #             names(sort(li_epi_cors, decreasing = TRUE))[
        #                 sort(li_epi_cors, decreasing = TRUE) > 0.5
        #             ]
        #         )
        #     )
        # )
        
    },
    
    simplify = FALSE,
    USE.NAMES = TRUE
    
)

# The following doesn't work if the analysis column is a factor - have to coerce to character.
# I prevented this by adding 'variable.factor = FALSE' into the call to melt() in the
# commonality_heatmap() function.

heatmap_commonality$scores[
    ,
    epi_score := expression_data[
        ids_genes_ordering_list[[analysis]]$sample_ids,
        as.numeric(
            cor(
                .SD[, gene, with = FALSE],
                rowMeans(.SD[, epi_markers[[analysis]], with = FALSE])
            )
        )
        ],
    by = analysis
    ]

pdf('../data_and_figures/emt_epi_cor.pdf')

for(an in unique(heatmap_commonality$scores$analysis)) {
    
    print(
        qplot(
            epi_score,
            score,
            data = heatmap_commonality$scores[analysis == an],
            main = an
        ) +
            geom_smooth(method = 'loess', span = 0.95, method.args = list(degree = 1))
    )
    
}

dev.off()





# Calculate EMT score predicted by loess model of EMT score ~ epi score, then take "adjusted"
# EMT score, which is the residual of this model:

heatmap_commonality$scores[
    ,
    emt_score_predicted := predict(
        loess(
            score ~ epi_score,
            span = 0.95,
            degree = 1
        )
    ),
    by = analysis
    ][
        ,
        emt_score_adjusted := score - emt_score_predicted
        ]

# heatmap_commonality$scores[
#     ,
#     emt_score_adjusted := score - predict(
#         loess(
#             score ~ epi_score,
#             span = 0.95,
#             degree = 1
#         )
#     ),
#     by = analysis
# ]

# Make heatmaps using the same genes as in the original commonality heatmap:

heatmap_commonality_predicted <- commonality_heatmap_from_scores(
    heatmap_commonality$scores[
        gene %in% heatmap_commonality$genes
        ],
    fill_var = 'emt_score_predicted',
    legend_title = 'EMT-CAF\nscore\npredicted\nby loess'
)

heatmap_commonality_adjusted <- commonality_heatmap_from_scores(
    heatmap_commonality$scores[
        gene %in% heatmap_commonality$genes
        ],
    fill_var = 'emt_score_adjusted',
    legend_title = 'EMT-CAF\nscore\nadjusted\nby loess'
)

# Alternatively, reselect the common genes:

heatmap_commonality_predicted_reselected <- commonality_heatmap_from_scores(
    heatmap_commonality$scores,
    fill_var = 'emt_score_predicted',
    reselect_genes = TRUE,
    ids_genes_ordering_list = ids_genes_ordering_list,
    legend_title = 'EMT-CAF\nscore\npredicted\nby loess'
)

heatmap_commonality_adjusted_reselected <- commonality_heatmap_from_scores(
    heatmap_commonality$scores,
    fill_var = 'emt_score_adjusted',
    reselect_genes = TRUE,
    ids_genes_ordering_list = ids_genes_ordering_list,
    legend_title = 'EMT-CAF\nscore\nadjusted\nby loess'
)

# Write to PDF:

pdf(
    '../data_and_figures/commonality_heatmap_adjusted.pdf',
    width = 10,
    height = 10
)

heatmap_commonality_adjusted$heatmap +
    labs(title = 'Adjusted for epithelial score using loess')

heatmap_commonality_predicted$heatmap +
    labs(title = 'EMT score predited by loess')

heatmap_commonality_adjusted_reselected$heatmap +
    labs(title = 'Adjusted for epithelial score using loess - reselected genes')

heatmap_commonality_predicted_reselected$heatmap +
    labs(title = 'EMT score predicted by loess - reselected genes')

dev.off()





# Check for correlations between gene and tumour EMT scores and scores for other cell types
# (for adjusted and non-adjusted scores):

# Gene scores:

score_cor_genes <- score_cor_with_cell_types(
    ids_genes_ordering_list,
    scores_data = heatmap_commonality$scores,
    score_var = 'score',
    title = 'Correlation between gene scores and cell types'
)

score_cor_genes_adjusted <- score_cor_with_cell_types(
    ids_genes_ordering_list,
    scores_data = heatmap_commonality$scores,
    score_var = 'emt_score_adjusted',
    title = 'Correlation between adjusted gene scores and\ncell types'
)

# Tumour scores:

# If you choose to transform the data, it doesn't make much difference subtracting the
# CAF score, and I think there isn't much conceptual motivation for doing this because
# the transformation should already account for the CAF content.

score_cor_tumours <- score_cor_with_cell_types(
    ids_genes_ordering_list,
    which_score = 'tumour',
    expression_data = expression_data,
    cell_type_markers = cell_type_markers,
    transform_data = TRUE,
    title = 'Correlation between tumour scores and cell\ntypes - transformed data'
)

score_cor_tumours_adjusted <- score_cor_with_cell_types(
    ids_genes_ordering_list,
    which_score = 'tumour',
    expression_data = expression_data,
    cell_type_markers = cell_type_markers,
    transform_data = TRUE,
    epi_markers_list = epi_markers,
    title = 'Correlation between adjusted tumour scores\nand cell types - transformed data'
)

# If you choose not to transform the data, you should subtract the CAF score.
# Otherwise, you get a lot of positive correlations, because you're also picking up
# the CAF content.  Subtracting the CAF scores turns most of these into negative
# correlations, which I guess is what you want - it means the CAFs correlate better
# with other cell types than your supposed cancer cell signature.

score_cor_tumours_untransformed <- score_cor_with_cell_types(
    ids_genes_ordering_list,
    which_score = 'tumour',
    expression_data = expression_data,
    cell_type_markers = cell_type_markers,
    transform_data = FALSE,
    subtract_caf_score = TRUE,
    title = 'Correlation between tumour scores and cell\ntypes - no transformation'
)

score_cor_tumours_untransformed_adjusted <- score_cor_with_cell_types(
    ids_genes_ordering_list,
    which_score = 'tumour',
    expression_data = expression_data,
    cell_type_markers = cell_type_markers,
    transform_data = FALSE,
    subtract_caf_score = TRUE,
    epi_markers_list = epi_markers,
    title = 'Correlation between adjusted tumour scores\nand cell types - no transformation'
)

# Save to PDF:

pdf(
    '../data_and_figures/cor_scores_with_cell_types.pdf',
    height = 7,
    width = 7
)

score_cor_genes$heatmap
score_cor_genes_adjusted$heatmap
score_cor_tumours$heatmap
score_cor_tumours_adjusted$heatmap
score_cor_tumours_untransformed$heatmap
score_cor_tumours_untransformed_adjusted$heatmap

dev.off()





# Make colour bars to plot with commonality heatmap to indicate correlations with other cell
# types:

# First get the ordering that makes the analyses from heatmap_commonality into the same order
# as those in score_cor_tumours:

order_heatmap_commonality <- order(order(heatmap_commonality$analyses))

# Now make colour bars:

colour_bars <- sapply(
    cell_types,
    function(ct) {
        col_side_colours(
            score_cor_tumours$heatmap_data[
                cell_type == ct,
                .(category = switch((correlation > 0) + 1, 0, 1)),
                by = analysis
                ]$category[
                    order_heatmap_commonality
                    ],
            ordering = heatmap_commonality$ordering_analyses,
            plot_margin = c(0, 5.5, 0, 5.5)
        ) +
            theme(axis.title.y = element_text(angle = 0, hjust = 1, vjust = 0.5)) +
            labs(y = ct)
    },
    simplify = FALSE,
    USE.NAMES = TRUE
)

# Make plot:

# I'm not sure how useful this is.  In any case, it might be better to have non-binary colour
# bars, showing the actual correlations instead of just 1 or 0, corresponding to +ve or -ve.

pdf(
    '../data_and_figures/commonality_heatmap_cell_type_cor.pdf',
    width = 10,
    height = 10
)

egg::ggarrange(
    plots = c(
        colour_bars,
        list(heatmap_commonality$heatmap)
    ),
    ncol = 1,
    nrow = length(colour_bars) + 1,
    heights = c(rep(1, length(colour_bars)), 50)
)

dev.off()





# Re-run clinical analyses with alternative genes/scores:

# First, use the gene EMT scores adjusted for epithelial markers to define new EMT markers
# for use in the clinical analyses (I think this is a bit primitive - it would be better
# to use actual adjusted EMT and CAF gene scores, but I don't already have these in
# heatmap_commonality$scores, where the scores are EMT minus CAF, so useless for this
# clinical analysis):

genes_for_scores <- sapply(
    names(ids_genes_ordering_list),
    function(li_name) {
        list(
            emt = head(
                heatmap_commonality$scores[
                    analysis == li_name &
                        gene %in%  ids_genes_ordering_list[[li_name]]$genes_filtered
                    ][
                        order(-emt_score_adjusted),
                        gene
                        ],
                20
            ),
            caf = tail(
                heatmap_commonality$scores[
                    analysis == li_name &
                        gene %in%  ids_genes_ordering_list[[li_name]]$genes_filtered
                    ][
                        order(-emt_score_adjusted),
                        gene
                        ],
                20
            )
        )
    },
    simplify = FALSE,
    USE.NAMES = TRUE
)

clin_cor_genes_adjusted = list(
    
    lymph_node_metastasis_positive = cor_with_clinical(
        expression_data,
        ids_genes_ordering_list,
        clinical_data,
        genes_for_scores = genes_for_scores,
        clin_var = 'number_of_lymphnodes_positive_by_he',
        wilcox_test_x_expr = variable > 0,
        wilcox_test_y_expr = variable == 0,
        plot_title = 'Correlation with lymph node metastasis'
    ),
    
    lymph_node_metastasis_multiple = cor_with_clinical(
        expression_data,
        ids_genes_ordering_list,
        clinical_data,
        genes_for_scores = genes_for_scores,
        clin_var = 'number_of_lymphnodes_positive_by_he',
        wilcox_test_x_expr = variable > 1,
        wilcox_test_y_expr = variable == 0,
        plot_title = 'Correlation with lymph node metastasis (0 vs. >1)'
    ),
    
    therapy_resistance = cor_with_clinical(
        expression_data,
        ids_genes_ordering_list,
        clinical_data,
        genes_for_scores = genes_for_scores,
        clin_var = 'followup_treatment_success',
        wilcox_test_x_expr = variable != 'complete remission/response',
        wilcox_test_y_expr = variable == 'complete remission/response',
        plot_title = 'Correlation with resistance to therapy'
    ),
    
    mortality = cor_with_clinical(
        expression_data,
        ids_genes_ordering_list,
        clinical_data,
        genes_for_scores = genes_for_scores,
        clin_var = 'days_to_death',
        wilcox_test_x_expr = variable < median(variable),
        wilcox_test_y_expr = variable > median(variable),
        plot_title = 'Correlation with mortality (median)'
    ),
    
    mortality_strict = cor_with_clinical(
        expression_data,
        ids_genes_ordering_list,
        clinical_data,
        genes_for_scores = genes_for_scores,
        clin_var = 'days_to_death',
        wilcox_test_x_expr = variable < quantile(variable, 0.4),
        wilcox_test_y_expr = variable > quantile(variable, 0.6),
        plot_title = 'Correlation with mortality (40th vs. 60th percentile)'
    ),
    
    lymphovascular_invasion = cor_with_clinical(
        expression_data,
        ids_genes_ordering_list,
        clinical_data,
        genes_for_scores = genes_for_scores,
        clin_var = 'lymphovascular_invasion',
        wilcox_test_x_expr = variable %in% c('yes', 'present'),
        wilcox_test_y_expr = variable %in% c('no', 'absent'),
        plot_title = 'Correlation with lymphovascular invasion'
    ),
    
    t_stage = cor_with_clinical(
        expression_data,
        ids_genes_ordering_list,
        clinical_data,
        genes_for_scores = genes_for_scores,
        clin_var = 'pathologic_t',
        wilcox_test_x_expr = startsWith(variable, 't3') |
            startsWith(variable, 't4'),
        wilcox_test_y_expr = startsWith(variable, 't0') |
            startsWith(variable, 't1') |
            startsWith(variable, 't2'),
        plot_title = 'Correlation with T stage',
        amatch_max_dist = 2
    ),
    
    n_stage = cor_with_clinical(
        expression_data,
        ids_genes_ordering_list,
        clinical_data,
        genes_for_scores = genes_for_scores,
        clin_var = 'pathologic_n',
        wilcox_test_x_expr = startsWith(variable, 'n2') |
            startsWith(variable, 'n3'),
        wilcox_test_y_expr = startsWith(variable, 'n0') |
            startsWith(variable, 'n1'),
        plot_title = 'Correlation with N stage',
        amatch_max_dist = 2
    ),
    
    m_stage = cor_with_clinical(
        expression_data,
        ids_genes_ordering_list,
        clinical_data,
        genes_for_scores = genes_for_scores,
        clin_var = 'pathologic_m',
        wilcox_test_x_expr = startsWith(variable, 'm1'),
        wilcox_test_y_expr = startsWith(variable, 'm0') |
            startsWith(variable, 'cm0'),
        plot_title = 'Correlation with M stage',
        amatch_max_dist = 2
    ),
    
    grade = cor_with_clinical(
        expression_data,
        ids_genes_ordering_list,
        clinical_data,
        genes_for_scores = genes_for_scores,
        clin_var = 'neoplasm_histologic_grade',
        wilcox_test_x_expr = variable %in% c('gb', 'g1', 'g2', 'low grade'),
        wilcox_test_y_expr = variable %in% c('g3', 'g4', 'high grade'),
        plot_title = 'Correlation with tumour grade'
    )
    
)

# Now using tumour scores:

tumour_scores <- calculate_tumour_scores(
    ids_genes_ordering_list,
    expression_data,
    epi_markers_list = epi_markers
)

# Non-adjusted tumour scores:

clin_cor_tumour_scores = list(
    
    lymph_node_metastasis_positive = cor_with_clinical(
        expression_data,
        ids_genes_ordering_list,
        clinical_data,
        scores_data = tumour_scores,
        clin_var = 'number_of_lymphnodes_positive_by_he',
        wilcox_test_x_expr = variable > 0,
        wilcox_test_y_expr = variable == 0,
        plot_title = 'Correlation with lymph node metastasis'
    ),
    
    lymph_node_metastasis_multiple = cor_with_clinical(
        expression_data,
        ids_genes_ordering_list,
        clinical_data,
        scores_data = tumour_scores,
        clin_var = 'number_of_lymphnodes_positive_by_he',
        wilcox_test_x_expr = variable > 1,
        wilcox_test_y_expr = variable == 0,
        plot_title = 'Correlation with lymph node metastasis (0 vs. >1)'
    ),
    
    therapy_resistance = cor_with_clinical(
        expression_data,
        ids_genes_ordering_list,
        clinical_data,
        scores_data = tumour_scores,
        clin_var = 'followup_treatment_success',
        wilcox_test_x_expr = variable != 'complete remission/response',
        wilcox_test_y_expr = variable == 'complete remission/response',
        plot_title = 'Correlation with resistance to therapy'
    ),
    
    mortality = cor_with_clinical(
        expression_data,
        ids_genes_ordering_list,
        clinical_data,
        scores_data = tumour_scores,
        clin_var = 'days_to_death',
        wilcox_test_x_expr = variable < median(variable),
        wilcox_test_y_expr = variable > median(variable),
        plot_title = 'Correlation with mortality (median)'
    ),
    
    mortality_strict = cor_with_clinical(
        expression_data,
        ids_genes_ordering_list,
        clinical_data,
        scores_data = tumour_scores,
        clin_var = 'days_to_death',
        wilcox_test_x_expr = variable < quantile(variable, 0.4),
        wilcox_test_y_expr = variable > quantile(variable, 0.6),
        plot_title = 'Correlation with mortality (40th vs. 60th percentile)'
    ),
    
    lymphovascular_invasion = cor_with_clinical(
        expression_data,
        ids_genes_ordering_list,
        clinical_data,
        scores_data = tumour_scores,
        clin_var = 'lymphovascular_invasion',
        wilcox_test_x_expr = variable %in% c('yes', 'present'),
        wilcox_test_y_expr = variable %in% c('no', 'absent'),
        plot_title = 'Correlation with lymphovascular invasion'
    ),
    
    t_stage = cor_with_clinical(
        expression_data,
        ids_genes_ordering_list,
        clinical_data,
        scores_data = tumour_scores,
        clin_var = 'pathologic_t',
        wilcox_test_x_expr = startsWith(variable, 't3') |
            startsWith(variable, 't4'),
        wilcox_test_y_expr = startsWith(variable, 't0') |
            startsWith(variable, 't1') |
            startsWith(variable, 't2'),
        plot_title = 'Correlation with T stage',
        amatch_max_dist = 2
    ),
    
    n_stage = cor_with_clinical(
        expression_data,
        ids_genes_ordering_list,
        clinical_data,
        scores_data = tumour_scores,
        clin_var = 'pathologic_n',
        wilcox_test_x_expr = startsWith(variable, 'n2') |
            startsWith(variable, 'n3'),
        wilcox_test_y_expr = startsWith(variable, 'n0') |
            startsWith(variable, 'n1'),
        plot_title = 'Correlation with N stage',
        amatch_max_dist = 2
    ),
    
    m_stage = cor_with_clinical(
        expression_data,
        ids_genes_ordering_list,
        clinical_data,
        scores_data = tumour_scores,
        clin_var = 'pathologic_m',
        wilcox_test_x_expr = startsWith(variable, 'm1'),
        wilcox_test_y_expr = startsWith(variable, 'm0') |
            startsWith(variable, 'cm0'),
        plot_title = 'Correlation with M stage',
        amatch_max_dist = 2
    ),
    
    grade = cor_with_clinical(
        expression_data,
        ids_genes_ordering_list,
        clinical_data,
        scores_data = tumour_scores,
        clin_var = 'neoplasm_histologic_grade',
        wilcox_test_x_expr = variable %in% c('gb', 'g1', 'g2', 'low grade'),
        wilcox_test_y_expr = variable %in% c('g3', 'g4', 'high grade'),
        plot_title = 'Correlation with tumour grade'
    )
    
)

# Adjusted tumour scores:

clin_cor_tumour_scores_adjusted = list(
    
    lymph_node_metastasis_positive = cor_with_clinical(
        expression_data,
        ids_genes_ordering_list,
        clinical_data,
        scores_data = tumour_scores,
        emt_score_var = 'emt_score_adjusted',
        clin_var = 'number_of_lymphnodes_positive_by_he',
        wilcox_test_x_expr = variable > 0,
        wilcox_test_y_expr = variable == 0,
        plot_title = 'Correlation with lymph node metastasis'
    ),
    
    lymph_node_metastasis_multiple = cor_with_clinical(
        expression_data,
        ids_genes_ordering_list,
        clinical_data,
        scores_data = tumour_scores,
        emt_score_var = 'emt_score_adjusted',
        clin_var = 'number_of_lymphnodes_positive_by_he',
        wilcox_test_x_expr = variable > 1,
        wilcox_test_y_expr = variable == 0,
        plot_title = 'Correlation with lymph node metastasis (0 vs. >1)'
    ),
    
    therapy_resistance = cor_with_clinical(
        expression_data,
        ids_genes_ordering_list,
        clinical_data,
        scores_data = tumour_scores,
        emt_score_var = 'emt_score_adjusted',
        clin_var = 'followup_treatment_success',
        wilcox_test_x_expr = variable != 'complete remission/response',
        wilcox_test_y_expr = variable == 'complete remission/response',
        plot_title = 'Correlation with resistance to therapy'
    ),
    
    mortality = cor_with_clinical(
        expression_data,
        ids_genes_ordering_list,
        clinical_data,
        scores_data = tumour_scores,
        emt_score_var = 'emt_score_adjusted',
        clin_var = 'days_to_death',
        wilcox_test_x_expr = variable < median(variable),
        wilcox_test_y_expr = variable > median(variable),
        plot_title = 'Correlation with mortality (median)'
    ),
    
    mortality_strict = cor_with_clinical(
        expression_data,
        ids_genes_ordering_list,
        clinical_data,
        scores_data = tumour_scores,
        emt_score_var = 'emt_score_adjusted',
        clin_var = 'days_to_death',
        wilcox_test_x_expr = variable < quantile(variable, 0.4),
        wilcox_test_y_expr = variable > quantile(variable, 0.6),
        plot_title = 'Correlation with mortality (40th vs. 60th percentile)'
    ),
    
    lymphovascular_invasion = cor_with_clinical(
        expression_data,
        ids_genes_ordering_list,
        clinical_data,
        scores_data = tumour_scores,
        emt_score_var = 'emt_score_adjusted',
        clin_var = 'lymphovascular_invasion',
        wilcox_test_x_expr = variable %in% c('yes', 'present'),
        wilcox_test_y_expr = variable %in% c('no', 'absent'),
        plot_title = 'Correlation with lymphovascular invasion'
    ),
    
    t_stage = cor_with_clinical(
        expression_data,
        ids_genes_ordering_list,
        clinical_data,
        scores_data = tumour_scores,
        emt_score_var = 'emt_score_adjusted',
        clin_var = 'pathologic_t',
        wilcox_test_x_expr = startsWith(variable, 't3') |
            startsWith(variable, 't4'),
        wilcox_test_y_expr = startsWith(variable, 't0') |
            startsWith(variable, 't1') |
            startsWith(variable, 't2'),
        plot_title = 'Correlation with T stage',
        amatch_max_dist = 2
    ),
    
    n_stage = cor_with_clinical(
        expression_data,
        ids_genes_ordering_list,
        clinical_data,
        scores_data = tumour_scores,
        emt_score_var = 'emt_score_adjusted',
        clin_var = 'pathologic_n',
        wilcox_test_x_expr = startsWith(variable, 'n2') |
            startsWith(variable, 'n3'),
        wilcox_test_y_expr = startsWith(variable, 'n0') |
            startsWith(variable, 'n1'),
        plot_title = 'Correlation with N stage',
        amatch_max_dist = 2
    ),
    
    m_stage = cor_with_clinical(
        expression_data,
        ids_genes_ordering_list,
        clinical_data,
        scores_data = tumour_scores,
        emt_score_var = 'emt_score_adjusted',
        clin_var = 'pathologic_m',
        wilcox_test_x_expr = startsWith(variable, 'm1'),
        wilcox_test_y_expr = startsWith(variable, 'm0') |
            startsWith(variable, 'cm0'),
        plot_title = 'Correlation with M stage',
        amatch_max_dist = 2
    ),
    
    grade = cor_with_clinical(
        expression_data,
        ids_genes_ordering_list,
        clinical_data,
        scores_data = tumour_scores,
        emt_score_var = 'emt_score_adjusted',
        clin_var = 'neoplasm_histologic_grade',
        wilcox_test_x_expr = variable %in% c('gb', 'g1', 'g2', 'low grade'),
        wilcox_test_y_expr = variable %in% c('g3', 'g4', 'high grade'),
        plot_title = 'Correlation with tumour grade'
    )
    
)

# Simpler EMT score minus CAF score:

tumour_scores_untransformed <- calculate_tumour_scores(
    ids_genes_ordering_list,
    expression_data,
    transform_data = FALSE,
    epi_markers_list = epi_markers
)

clin_cor_tumour_scores_untransformed_adjusted = list(
    
    lymph_node_metastasis_positive = cor_with_clinical(
        expression_data,
        ids_genes_ordering_list,
        clinical_data,
        scores_data = tumour_scores_untransformed,
        emt_score_var = 'emt_score_adjusted',
        clin_var = 'number_of_lymphnodes_positive_by_he',
        wilcox_test_x_expr = variable > 0,
        wilcox_test_y_expr = variable == 0,
        plot_title = 'Correlation with lymph node metastasis'
    ),
    
    lymph_node_metastasis_multiple = cor_with_clinical(
        expression_data,
        ids_genes_ordering_list,
        clinical_data,
        scores_data = tumour_scores_untransformed,
        emt_score_var = 'emt_score_adjusted',
        clin_var = 'number_of_lymphnodes_positive_by_he',
        wilcox_test_x_expr = variable > 1,
        wilcox_test_y_expr = variable == 0,
        plot_title = 'Correlation with lymph node metastasis (0 vs. >1)'
    ),
    
    therapy_resistance = cor_with_clinical(
        expression_data,
        ids_genes_ordering_list,
        clinical_data,
        scores_data = tumour_scores_untransformed,
        emt_score_var = 'emt_score_adjusted',
        clin_var = 'followup_treatment_success',
        wilcox_test_x_expr = variable != 'complete remission/response',
        wilcox_test_y_expr = variable == 'complete remission/response',
        plot_title = 'Correlation with resistance to therapy'
    ),
    
    mortality = cor_with_clinical(
        expression_data,
        ids_genes_ordering_list,
        clinical_data,
        scores_data = tumour_scores_untransformed,
        emt_score_var = 'emt_score_adjusted',
        clin_var = 'days_to_death',
        wilcox_test_x_expr = variable < median(variable),
        wilcox_test_y_expr = variable > median(variable),
        plot_title = 'Correlation with mortality (median)'
    ),
    
    mortality_strict = cor_with_clinical(
        expression_data,
        ids_genes_ordering_list,
        clinical_data,
        scores_data = tumour_scores_untransformed,
        emt_score_var = 'emt_score_adjusted',
        clin_var = 'days_to_death',
        wilcox_test_x_expr = variable < quantile(variable, 0.4),
        wilcox_test_y_expr = variable > quantile(variable, 0.6),
        plot_title = 'Correlation with mortality (40th vs. 60th percentile)'
    ),
    
    lymphovascular_invasion = cor_with_clinical(
        expression_data,
        ids_genes_ordering_list,
        clinical_data,
        scores_data = tumour_scores_untransformed,
        emt_score_var = 'emt_score_adjusted',
        clin_var = 'lymphovascular_invasion',
        wilcox_test_x_expr = variable %in% c('yes', 'present'),
        wilcox_test_y_expr = variable %in% c('no', 'absent'),
        plot_title = 'Correlation with lymphovascular invasion'
    ),
    
    t_stage = cor_with_clinical(
        expression_data,
        ids_genes_ordering_list,
        clinical_data,
        scores_data = tumour_scores_untransformed,
        emt_score_var = 'emt_score_adjusted',
        clin_var = 'pathologic_t',
        wilcox_test_x_expr = startsWith(variable, 't3') |
            startsWith(variable, 't4'),
        wilcox_test_y_expr = startsWith(variable, 't0') |
            startsWith(variable, 't1') |
            startsWith(variable, 't2'),
        plot_title = 'Correlation with T stage',
        amatch_max_dist = 2
    ),
    
    n_stage = cor_with_clinical(
        expression_data,
        ids_genes_ordering_list,
        clinical_data,
        scores_data = tumour_scores_untransformed,
        emt_score_var = 'emt_score_adjusted',
        clin_var = 'pathologic_n',
        wilcox_test_x_expr = startsWith(variable, 'n2') |
            startsWith(variable, 'n3'),
        wilcox_test_y_expr = startsWith(variable, 'n0') |
            startsWith(variable, 'n1'),
        plot_title = 'Correlation with N stage',
        amatch_max_dist = 2
    ),
    
    m_stage = cor_with_clinical(
        expression_data,
        ids_genes_ordering_list,
        clinical_data,
        scores_data = tumour_scores_untransformed,
        emt_score_var = 'emt_score_adjusted',
        clin_var = 'pathologic_m',
        wilcox_test_x_expr = startsWith(variable, 'm1'),
        wilcox_test_y_expr = startsWith(variable, 'm0') |
            startsWith(variable, 'cm0'),
        plot_title = 'Correlation with M stage',
        amatch_max_dist = 2
    ),
    
    grade = cor_with_clinical(
        expression_data,
        ids_genes_ordering_list,
        clinical_data,
        scores_data = tumour_scores_untransformed,
        emt_score_var = 'emt_score_adjusted',
        clin_var = 'neoplasm_histologic_grade',
        wilcox_test_x_expr = variable %in% c('gb', 'g1', 'g2', 'low grade'),
        wilcox_test_y_expr = variable %in% c('g3', 'g4', 'high grade'),
        plot_title = 'Correlation with tumour grade'
    )
    
)

# Save to PDF:

pdf('../data_and_figures/clinical_analyses_alternative.pdf', width = 35, height = 14)

egg::ggarrange(
    plots = lapply(
        clin_cor_genes_adjusted,
        function(li) li$plot
    ),
    nrow = 2,
    ncol = 5,
    top = 'Correlations of gene scores with clinical features (genes chosen from adjusted EMT-CAF scores)'
)

egg::ggarrange(
    plots = lapply(
        clin_cor_tumour_scores,
        function(li) li$plot
    ),
    nrow = 2,
    ncol = 5,
    top = 'Correlations of tumour scores with clinical features'
)

egg::ggarrange(
    plots = lapply(
        clin_cor_tumour_scores_adjusted,
        function(li) li$plot
    ),
    nrow = 2,
    ncol = 5,
    top = 'Correlations of adjusted tumour scores with clinical features'
)

egg::ggarrange(
    plots = lapply(
        clin_cor_tumour_scores_untransformed_adjusted,
        function(li) li$plot
    ),
    nrow = 2,
    ncol = 5,
    top = 'Correlations of tumour EMT minus Epi scores with clinical features'
)

dev.off()





# Use the following to make boxplots of epi, EMT and CAF gene expression in single cell data (change
# cancer type names and single cell dataset names as appropriate):

# (Probably heatmaps would be better, of the sort that I previously made for single cell data.)

li <- fread('../li_colon_adenocarcinoma_2017_cells_as_rows.csv')

ggplot(
    data = li[
        ,
        .(ave_epi = rowMeans(.SD)),
        .SDcols = epi_markers$Colon[
            epi_markers$Colon %in% names(li)
            ],
        by = cell_type
        ]
) +
    geom_boxplot(aes(x = cell_type, y = ave_epi))

ggplot(
    data = li[
        ,
        .(ave_epi = rowMeans(.SD)),
        .SDcols = head(all_deconvs$coadread$genes_filtered[all_deconvs$coadread$ordering], 20)[
            head(all_deconvs$coadread$genes_filtered[all_deconvs$coadread$ordering], 20) %in% names(li)
            ],
        by = cell_type
        ]
) +
    geom_boxplot(aes(x = cell_type, y = ave_epi))

ggplot(
    data = li[
        ,
        .(ave_epi = rowMeans(.SD)),
        .SDcols = tail(all_deconvs$coadread$genes_filtered[all_deconvs$coadread$ordering], 20)[
            tail(all_deconvs$coadread$genes_filtered[all_deconvs$coadread$ordering], 20) %in% names(li)
            ],
        by = cell_type
        ]
) +
    geom_boxplot(aes(x = cell_type, y = ave_epi))





# Alternative EMT scores:

scores_data <- cbind(
    expression_data[
        li_sample_ids,
        ..genes_for_scores_data
        ],
    emt_sum = rowSums(
        expression_data[
            li_sample_ids,
            ..emt_genes
            ]
    ),
    caf_sum = rowSums(
        expression_data[
            li_sample_ids,
            ..caf_genes
            ]
    )
)

scores <- rbindlist(
    sapply(
        
        names(ids_genes_ordering_list),
        
        function(li_name) {
            
            li <- ids_genes_ordering_list[[li_name]]
            
            li_sample_ids <- li[[pmatch('sample', names(li))]]
            
            li_ordered_genes <- li[[pmatch('gene', names(li))]][
                li[[pmatch('order', names(li))]]
                ]
            
            emt_genes <- head(li_ordered_genes, 20)
            caf_genes <- tail(li_ordered_genes, 20)
            
            scores_data <- cbind(
                expression_data[
                    li_sample_ids,
                    ..genes_for_scores_data
                    ],
                emt_sum = rowSums(
                    expression_data[
                        li_sample_ids,
                        ..emt_genes
                        ]
                ),
                caf_sum = rowSums(
                    expression_data[
                        li_sample_ids,
                        ..caf_genes
                        ]
                )
            )
            
            scores_data <- data.table(
                gene = genes_for_scores_data,
                analysis = li_name,
                emt_score = as.data.table(
                    sapply(
                        genes_for_scores_data,
                        function(g) {
                            
                            lm(
                                formula(
                                    paste0(
                                        '`',
                                        g,
                                        '` ~ caf_sum'
                                    )
                                ),
                                data = scores_data
                            )$residuals
                            
                        },
                        USE.NAMES = TRUE
                    )
                )[
                    ,
                    rowMeans(cor(.SD, .SD[, ..emt_genes]))
                    ],
                caf_score = as.data.table(
                    sapply(
                        genes_for_scores_data,
                        function(g) {
                            
                            lm(
                                formula(
                                    paste0(
                                        '`',
                                        g,
                                        '` ~ emt_sum'
                                    )
                                ),
                                data = scores_data
                            )$residuals
                            
                        },
                        USE.NAMES = TRUE
                    )
                )[
                    ,
                    rowMeans(cor(.SD, .SD[, ..caf_genes]))
                    ]
            )
            
            scores_data
            
        },
        
        simplify = FALSE,
        USE.NAMES = TRUE
        
    )
)

scores[, score := emt_score - caf_score]

cast_data <- dcast(
    scores[gene %in% heatmap_commonality$genes],
    gene ~ analysis,
    value.var = 'score'
)

ordering_y <- hclust(
    dist(cast_data[, -'gene']),
    method = 'average'
)$order

ordering_x <- hclust(
    dist(t(cast_data[, -'gene'])),
    method = 'average'
)$order

setkey(scores, gene)

g <- ggplot(
    scores[cast_data[, gene][ordering_y]],
    aes(
        x = factor(
            analysis,
            levels = names(cast_data[, -'gene'])[ordering_x]
        ),
        y = factor(
            gene,
            levels = unique(gene)
        )
    )
) +
    geom_raster(aes(fill = score)) +
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
        plot.margin = unit(c(5.5, 10, 5.5, 0), 'pt')
    ) +
    labs(
        x = '',
        fill = 'EMT-CAF\nscore'
    )

g <- g +
    scale_fill_gradientn(
        colours = rev(colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))(50)),
        limits = c(-1, 1)
    )
