# Test of decomposing a density into two normal components:

# This was applied to the part of the deconvove_emt_caf() function just before we define
# genes_filtered, i.e.:

# scores_table <- genes_cor_with_initial_and_cell_types[
#     ,
#     .(
#         score = - sum(
#             as.numeric(
#                 .SD[, ..cell_types]
#             )*as.numeric(
#                 cell_type_weights[
#                     cell_types # Makes the cell types in the right order
#                     ]
#             )
#         )
#     ),
#     by = id
# ][
#     order(-score)
# ]

# It seems to work well in the HNSC case, but if you follow through with the filtered
# gene set that you end up with, it doesn't give a good deconvolution.  This and the fact
# that we'd probably have to check each cancer type for appropriate parameters suggests
# that it may not be the best method to use.  But it's interesting and I feel I should be
# able to apply it somewhere, perhaps in the rare vs. shared EMT analysis.

mix_mod <- mixtools::normalmixEM(scores_table$score)

densities <- lapply(
    1:2,
    function(i) {
        mix_mod$lambda[i]*dnorm(
            mix_mod$x,
            mean = mix_mod$mu[i],
            sd = mix_mod$sigma[i]
        )
    }
)

# In the following, we have to use min() instead of max() because the values of min_mod$x
# are in descending order, and densities[[1]] and densities[[2]] have the corresponding
# order.  We want to take the scores which are strictly greater than this cutoff.

cutoff <- mix_mod$x[min(which(densities[[2]] - densities[[1]] <= 0))]

# To visualise:

qplot(
    x = x,
    y = y,
    col = col,
    data = data.table(
        x = rep(mix_mod$x, 3),
        y = c(
            unlist(densities),
            pmax(0, densities[[2]] - densities[[1]])
        ),
        col = unlist(
            lapply(
                c('a', 'b', 'c'),
                rep,
                length(mix_mod$x)
            )
        )
    ),
    geom = 'line'
) +
    geom_vline(aes(xintercept = cutoff), colour = 'gold', size = 1) +
    theme_test()

# Define filtered gene set by :

# genes_filtered <- scores_table[score > cutoff, id]





# Discussion:

# I think perhaps the reason this didn't work so well in the deconv analysis was that the
# genes all look more like cancer markers than fibroblast markers.  Perhaps we could use
# the mixture model to actually define the cancer cell component and the CAF component,
# or more generally the TME component - the second peak could be interpreted as the cancer
# cell markers, and the first as the TME components.  Perhaps, instead of taking the
# cutoff at the intersection of the two distributions, we could take it at, say, the mean
# of the first distribution.  Then we should keep two components - two distributions -
# while filtering out genes that seem to be strongly associated with other components of
# the TME (we could interpret genes below the mean of the first peak to be strongly
# associated with other TME components). An illustration of this:

qplot(
    x = x,
    y = y,
    col = col,
    data = data.table(
        x = rep(mix_mod$x, 3),
        y = c(
            unlist(densities),
            pmax(0, densities[[2]] - densities[[1]])
        ),
        col = unlist(
            lapply(
                c('a', 'b', 'c'),
                rep,
                length(mix_mod$x)
            )
        )
    ),
    geom = 'line'
) +
    geom_vline(
        aes(xintercept = cutoff),
        colour = 'gold',
        size = 1
    ) +
    geom_vline(
        aes(xintercept = weighted.mean(mix_mod$x, densities[[1]])),
        colour = 'gold',
        size = 1
    ) +
    theme_test()

genes_filtered <- scores_table[score > weighted.mean(mix_mod$x, densities[[1]]), id]

# The weighted.mean() function doesn't perfectly find the peak of the first distribution,
# but it's probably good enough.  The final deconv result looks pretty good!
