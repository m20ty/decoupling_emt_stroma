inter_intra_emt_all_cts <- sapply(
    # names(single_cell_metadata),
    c('hnsc', 'lung', 'paad', 'lihc'), # Exclude brca and coadread (too few cancer cells)
    function(ct) {
        
        sc_data <- eval(single_cell_metadata[[ct]]$read_quote)
        
        setkey(sc_data, id)
        
        # Take subset of tumours with at least 50 cancer cells:
        
        sc_data <- sc_data[
            patient %in% sc_data[
                cell_type == 'cancer',
                .(N = .N),
                by = patient
                ][
                    N >= 50,
                    patient
                    ]
            ]
        
        # Get all the sufficiently highly-expressed genes:
        
        all_genes_filtered <- sc_data[
            cell_type == 'cancer',
            names(.SD)[
                apply(.SD, 2, sc_cancer_args[[ct]]$genes_filter_fun)
                ],
            .SDcols = -c('id', 'patient', 'cell_type')
            ]
        
        # Compute bins for the genes based on average expression:
        
        sc_mat <- sc_data[
            cell_type == 'cancer',
            set_colnames(t(.SD), id),
            .SDcols = all_genes_filtered
            ]
        
        gene_averages <- sort(rowMeans(sc_mat))
        
        bins <- setNames(
            cut(
                seq_along(gene_averages),
                breaks = length(gene_averages) %/% 110,
                labels = FALSE,
                include.lowest = TRUE
            ),
            names(gene_averages)
        )
        
        # Get EMT markers:
        
        ct_emt_markers <- unique(
            unlist(
                lapply(
                    deconv_data[
                        grepl(
                            paste(
                                paste0('^', single_cell_metadata[[ct]]$tcga_cancer_types),
                                collapse = '|'
                            ),
                            names(deconv_data)
                        )
                        ],
                    function(deconv) deconv$genes_filtered[deconv$ordering]
                    # function(deconv) {
                    #     head(
                    #         deconv$genes_filtered[deconv$ordering],
                    #         length(deconv$genes_filtered)/3
                    #     )
                    # }
                )
            )
        )
        
        ct_emt_markers <- ct_emt_markers[ct_emt_markers %in% names(sc_data)]
        
        # Filter EMT markers by expression levels:
        
        ct_emt_markers <- sc_data[
            cell_type == 'cancer',
            ct_emt_markers[
                apply(.SD, 2, sc_cancer_args[[ct]]$scores_filter_fun)
                ],
            .SDcols = ct_emt_markers
            ]
        
        # Define control gene sets for distribution of scores:
        
        comparable_gene_sets <- lapply(
            ct_emt_markers,
            function(g) {
                sample(names(bins)[bins == bins[g]], 100)
            }
        )
        
        comparable_gene_sets <- lapply(
            1:100,
            function(i) sapply(comparable_gene_sets, `[`, i)
        )
        
        # Calculate EMT scores, then filter the EMT markers for correlation with these
        # EMT scores, and recalculate the EMT scores using the filtered list.  I think
        # it makes more sense to filter for correlation with the initial EMT score than
        # with the Z score (below), since the Z score won't be on a comparable scale
        # to the original distribution of expression levels (though this probably
        # doesn't matter much), and the aim of the Z score is to remove variability
        # arising from poor data quality, which might significantly affect the
        # correlation with the Z scores.  EDIT: I don't think there's any difference
        # between the correlations with EMT scores and with Z scores.
        
        emt_scores <- rowMeans(
            sapply(
                ct_emt_markers,
                function(g) {
                    sc_mat[g, ] - colMeans(
                        sc_mat[sample(names(bins)[bins == bins[g]], 100), ]
                    )
                },
                USE.NAMES = TRUE
            )
        )
        
        # sc_data[
        #     cell_type == 'cancer',
        #     sapply(ct_emt_markers, function(g) cor(get(g), emt_scores))
        # ] %>%
        #     sort %>%
        #     plot
        
        ct_emt_markers <- sc_data[
            cell_type == 'cancer',
            ct_emt_markers[
                sapply(ct_emt_markers, function(g) cor(get(g), emt_scores)) > 0.3
                ]
            ]
        
        # Calculate the EMT scores per tumour:
        
        emt_scores <- lapply(
            unique(sc_data$patient),
            function(p) {
                
                gene_averages <- sort(rowMeans(sc_mat[, sc_data[cell_type == 'cancer' & patient == p, id]]))
                
                bins <- setNames(
                    cut(
                        seq_along(gene_averages),
                        breaks = length(gene_averages) %/% 110,
                        labels = FALSE,
                        include.lowest = TRUE
                    ),
                    names(gene_averages)
                )
                
                comparable_gene_sets <- lapply(
                    ct_emt_markers,
                    function(g) {
                        sample(names(bins)[bins == bins[g]], 100)
                    }
                )
                
                comparable_gene_sets <- lapply(
                    1:100,
                    function(i) sapply(comparable_gene_sets, `[`, i)
                )
                
                emt_scores <- rowMeans(
                    sapply(
                        ct_emt_markers,
                        function(g) {
                            sc_mat[g, sc_data[cell_type == 'cancer' & patient == p, id]] - colMeans(
                                sc_mat[sample(names(bins)[bins == bins[g]], 100), sc_data[cell_type == 'cancer' & patient == p, id]]
                            )
                        },
                        USE.NAMES = TRUE
                    )
                )
                
                comparable_gene_sets_scores <- lapply(
                    comparable_gene_sets,
                    function(gli) {
                        rowMeans(
                            sapply(
                                gli,
                                function(g) {
                                    sc_mat[g, sc_data[cell_type == 'cancer' & patient == p, id]] - colMeans(
                                        sc_mat[sample(names(bins)[bins == bins[g]], 100), sc_data[cell_type == 'cancer' & patient == p, id]]
                                    )
                                },
                                USE.NAMES = TRUE
                            )
                        )
                    }
                )
                
                z_scores <- sapply(
                    names(emt_scores),
                    function(cell_id) {
                        
                        distrib <- c(
                            emt_scores[cell_id],
                            sapply(comparable_gene_sets_scores, `[`, cell_id)
                        )
                        
                        (emt_scores[cell_id] - mean(distrib))/sd(distrib)
                        
                    },
                    USE.NAMES = FALSE
                )
                
                emt_scores <- sc_data[
                    cell_type == 'cancer' & patient == p,
                    .(
                        id = id,
                        patient = patient,
                        emt_score = z_scores[id]
                    )
                ][
                    order(patient, emt_score),
                    pos_frac := (1:.N)/.N,
                    by = patient
                ]
                
                emt_scores
                
            }
        ) %>% rbindlist
        
        # emt_scores <- rowMeans(
        #     sapply(
        #         ct_emt_markers,
        #         function(g) {
        #             sc_mat[g, ] - colMeans(
        #                 sc_mat[sample(names(bins)[bins == bins[g]], 100), ]
        #             )
        #         },
        #         USE.NAMES = TRUE
        #     )
        # )
        
        # Calculate distribution of scores for control gene sets:
        
        # comparable_gene_sets_scores <- lapply(
        #     comparable_gene_sets,
        #     function(gli) {
        #         rowMeans(
        #             sapply(
        #                 gli,
        #                 function(g) {
        #                     sc_mat[g, ] - colMeans(
        #                         sc_mat[sample(names(bins)[bins == bins[g]], 100), ]
        #                     )
        #                 },
        #                 USE.NAMES = TRUE
        #             )
        #         )
        #     }
        # )
        
        # Define the EMT score as a Z score:
        
        # z_scores <- sapply(
        #     names(emt_scores),
        #     function(cell_id) {
        #         
        #         distrib <- c(
        #             emt_scores[cell_id],
        #             sapply(comparable_gene_sets_scores, `[`, cell_id)
        #         )
        #         
        #         (emt_scores[cell_id] - mean(distrib))/sd(distrib)
        #         
        #     },
        #     USE.NAMES = FALSE
        # )
        # 
        # emt_scores <- sc_data[
        #     cell_type == 'cancer',
        #     .(
        #         id = id,
        #         patient = patient,
        #         emt_score = z_scores[id]
        #     )
        #     ][
        #         order(patient, emt_score),
        #         pos_frac := (1:.N)/.N,
        #         by = patient
        #         ]
        
        # emt_scores <- scrabble::score(
        #     sc_data[
        #         cell_type == 'cancer',
        #         set_colnames(t(.SD), id),
        #         .SDcols = all_genes_filtered
        #     ],
        #     list(ct_emt_markers),
        #     bin.control = TRUE,
        #     nbin = length(all_genes_filtered) %/% 110
        # )
        
        # emt_scores <- data.table(
        #     id = rownames(emt_scores),
        #     patient = sc_data[rownames(emt_scores), patient],
        #     emt_score = emt_scores[, 1]
        # )[
        #     order(patient, emt_score),
        #     pos_frac := (1:.N)/.N,
        #     by = patient
        # ]
        
        emt_score_line <- ggplot(
            emt_scores[pos_frac > 0.01 & pos_frac < 0.99],
            aes(pos_frac, emt_score, group = patient, colour = as.character(patient))
        ) +
            scale_colour_manual(
                values = RColorBrewer::brewer.pal(12, 'Set3')[
                    c(1, 3:10, 12, 11)
                    ][
                        1:length(unique(emt_scores$patient))
                        ]
            ) +
            geom_line() +
            theme_test()
        
        # I think a violin plot might be more illustrative than a lineplot:
        
        emt_score_violin <- ggplot(
            emt_scores,
            aes(
                factor(
                    patient,
                    levels = emt_scores[
                        ,
                        .(med = median(emt_score)),
                        by = patient
                        ][
                            order(med),
                            patient
                            ]
                ),
                emt_score,
                fill = as.character(patient)
            )
        ) +
            scale_fill_manual(
                values = RColorBrewer::brewer.pal(12, 'Set3')[
                    1:length(unique(emt_scores$patient))
                    ]
            ) +
            geom_violin(draw_quantiles = 0.5) +
            theme_test()
        
        inter_intra_emt <- emt_scores[
            ,
            .(
                # inter_emt = .SD[pos_frac > 0.25 & pos_frac < 0.75, mean(emt_score)],
                # intra_emt = .SD[pos_frac > 0.9 & pos_frac < 0.975, mean(emt_score)] -
                #     median(emt_score)
                inter_emt = median(emt_score),
                intra_emt = quantile(emt_score, 0.95) - median(emt_score)
            ),
            by = patient
            ]
        
        inter_intra_plot <- ggplot(
            inter_intra_emt,
            aes(inter_emt, intra_emt, colour = as.character(patient))
        ) +
            scale_colour_manual(
                values = RColorBrewer::brewer.pal(12, 'Set3')[
                    c(1, 3:10, 12, 11)
                    ][
                        1:length(unique(emt_scores$patient))
                        ]
            ) +
            geom_point() +
            theme_test()
        
        list(
            plots = list(
                lineplot = emt_score_line,
                violin = emt_score_violin,
                scatterplot = inter_intra_plot
            ),
            data = list(
                cell_emt_scores = emt_scores,
                inter_intra_emt_scores = inter_intra_emt,
                all_genes_filtered = all_genes_filtered,
                emt_markers_filtered = ct_emt_markers
            )
        )
        
    },
    simplify = FALSE,
    USE.NAMES = TRUE
)
