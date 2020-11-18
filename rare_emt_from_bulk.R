# Run the first part of further_analysis.R to get all_deconvs.

# Try first with HNSC:

deconv_ct <- all_deconvs$hnsc

mes_data <- expression_data[
    deconv_ct$sample_ids,
    deconv_ct$genes_filtered[deconv_ct$ordering],
    with = FALSE
]

# Remove CAF contribution:

caf_score <- mes_data[
    ,
    rowMeans(.SD),
    .SDcols = tail(names(mes_data), 20)
]

resid_data <- as.data.table(
    sapply(
        names(mes_data),
        function(g) {
            lm(
                formula(
                    paste0(
                        '`',
                        g,
                        '` ~ caf_score'
                    )
                ),
                data = mes_data
            )$residuals
        },
        simplify = FALSE,
        USE.NAMES = TRUE
    )
)

# Find correlations of the genes with the row sums of resid_data.  This could be considered
# an analogy of what we do with the single cell data, where we find the correlation among
# the cancer cells of the mes genes with the mes signature scores.

mes_cor <- resid_data[, cor(.SD, rowSums(.SD))]

sort(mes_cor[, 1], decreasing = TRUE)

# This kind of works, but we still get SNAI2 near the top, and VIM is not as high up as I
# would expect.  We also don't get a very good correspondance with Itay's EMT signature:

sum(
    names(sort(temp_cor[, 1], decreasing = TRUE))[1:50] %in%
        cell_type_markers[source == 'Tirosh' & cell_type == 'mesenchymal', gene]
)
