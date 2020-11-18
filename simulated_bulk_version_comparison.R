# This script is about comparing two versions of my code for simulating bulk profiles.
# The motivation was to understand why the deconvolution of simulated bulk for
# colorectal seemed much weaker in the later version.  I think in the end the two
# versions of the code are equivalent: they both sample a count for an "initial" cell
# type ('cancer' in both versions - no other cell types are used as initial ones)
# using a normal distribution, then sample the rest of the cell types using a
# uniform distribution and scale them appropriately so that the initial cell type
# count is the right proportion.  The difference between the two versions is not in
# the sampling method - it seems to be simply a change of parameter, with
# max_mean_count being 200 in the earlier version and 20 in the later version.
# Perhaps I should experiment a bit with the parameters.  The only other differences
# between the versions are to do with speed (I'm pretty sure the newer version is
# faster) and things that are not relevant to my present analysis (e.g. in the
# newer version you can set multiple initial cell types together by sampling multiple
# random uniforms and scaling them to match a single random normal, before again
# sampling the other cell types with random uniforms and also scaling them to fit
# this one random normal - this is irrelevant for me since I only use one initial
# cell type).





# New version:

# sim_data <- simulated_tumours_data(
#     as.matrix(sc_data[, -c('id', 'patient', 'cell_type')]),
#     types = sc_data$cell_type,
#     id_prefix = ct,
#     max_mean_count = max_mean_counts[[ct]]
# )

proportions_vec <- seq(0.1, 0.9, length.out = 40)
n <- 25
sd_frac <- 0.2

density_fun <- runif
initial_types <- 'cancer'

max_mean_count <- max_mean_counts[[ct]] # This is 20 in the case of COADREAD

types <- sc_data$cell_type
cell_types <- unique(sc_data$cell_type)

density_fun_args <- sapply(
    cell_types,
    function(x) { # This x isn't actually used.
        list(
            min = 0,
            max = max_mean_count
        )
    },
    simplify = FALSE,
    USE.NAMES = TRUE
)

sc_mat <- as.matrix(sc_data[, -c('id', 'patient', 'cell_type')])

# simulated_counts <- simulate_counts(
#     cell_types,
#     initial_types = initial_types,
#     proportions_vec = proportions_vec,
#     n = n,
#     ...
# )

# For each cell type, generate n*length(proportions_vec) (= 25*40) random uniforms,
# and store the results as a data table (so one column for each cell type):

counts_table <- as.data.table( # Generate counts for each type using density_fun:
    sapply(
        cell_types,
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

other_types <- cell_types[!(cell_types %in% initial_types)]

# Generate combined count for initial_types by sampling from the set of integers
# 1:(2*max_mean_count) using a normal distribution with mean and sd derived from
# max_mean_count and sd_frac, scaled using the supplied proportions:

# Sample n random uniforms for each proportion (so n*length(proportions_vec) = 25*40
# random uniforms in total), bind it to the table of random uniforms, and scale
# counts appropriately.  The random normal will become the combined count of the
# initial types (which in my case is just one count, so this is irrelevant), the
# individual counts for which come from the random uniforms (which are scaled, but
# I guess if we just have one initial type the random uniform and random normal
# values become the same).

counts_table <- cbind(
    rbindlist(
        lapply(
            proportions_vec,
            function(p) { # Is anything actually returned from this function?
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
][
    ,
    ..cell_types
]

# Don't need this line any more, because we only use counts_table:

# simulated_counts <- list(counts_table = counts_table, initial_types = initial_types)

# sampled_indices <- sample_indices(
#     types,
#     counts_table
# )

sampled_indices <- apply(
    counts_table,
    1,
    function(s) {
        sapply(
            names(s),
            function(type) {
                
                N <- sum(types == type)
                inds <- which(types == type)
                
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

sc_mat <- round(2^sc_mat - 1, 4)

expr <- t( # Transpose because genes end up as rows after sapply()
    sapply(
        1:length(sampled_indices),
        function(i) colSums(sc_mat[unlist(sampled_indices[[i]]), ])
    )
)

meta <- counts_table[
    ,
    .(purity = cancer/sum(as.numeric(.SD))),
    by = .(id = row.names(counts_table)),
    .SDcols = cell_types
]

expr <- round(log2(expr + 1), 4)

sim_data <- list(
    expression_data = expr,
    meta_data = meta
)





# Old version:

# simulated_counts <- simulate_counts_from_type_proportions(
#     unique(single_cell_data$cell_type),
#     initial_cell_types = 'cancer',
#     proportions_vec = seq(0.1, 0.9, length.out = 40),
#     n = 25,
#     max_mean_count = 200
# )$cancer[
#     ,
#     -c('proportion', 'median_proportion')
# ]

proportions_vec <- seq(0.1, 0.9, length.out = 40)
n <- 25
sd_frac <- 0.2

densities <- runif
initial_cell_types <- 'cancer'

max_mean_count <- 200

types <- unique(sc_data$cell_type)

density_fun_args <- sapply(
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

simulated_counts <- sapply(
    initial_cell_types,
    function(ct) {
        
        other_types <- types[types != ct]
        
        initial <- rbindlist(
            lapply(
                proportions_vec,
                function(p) {
                    
                    # sim_counts <- simulate_counts(
                    #     types = ct,
                    #     densities = function(n, max_mean_count, sd_frac) {
                    #         sample(
                    #             # I'm choosing 1 as the minimum because it makes no sense
                    #             # to sample zero of the initial cell type.
                    #             1:(2*max_mean_count),
                    #             size = n,
                    #             replace = TRUE,
                    #             prob = dnorm(
                    #                 1:(2*max_mean_count),
                    #                 mean = max_mean_count*p/max(proportions_vec),
                    #                 sd = sd_frac*max_mean_count*p/max(proportions_vec)
                    #             )
                    #         )
                    #     },
                    #     n = n,
                    #     density_fun_args = setNames(
                    #         list(
                    #             list(
                    #                 max_mean_count = max_mean_count,
                    #                 sd_frac = sd_frac
                    #             )
                    #         ),
                    #         ct
                    #     )
                    # )[
                    #     ,
                    #     proportion := p
                    # ]
                    
                    # Here we sample n random normals for each proportion, meaning
                    # n*40 = 25*40 random normals in total.  These will be the counts
                    # of the initial cell type(s).
                    
                    densities <- function(n, max_mean_count, sd_frac) {
                        sample(
                            # I'm choosing 1 as the minimum because it makes no sense
                            # to sample zero of the initial cell type.
                            1:(2*max_mean_count),
                            size = n,
                            replace = TRUE,
                            prob = dnorm(
                                1:(2*max_mean_count),
                                mean = max_mean_count*p/max(proportions_vec),
                                sd = sd_frac*max_mean_count*p/max(proportions_vec)
                            )
                        )
                    }
                    
                    density_fun_args <- setNames(
                        list(
                            list(
                                max_mean_count = max_mean_count,
                                sd_frac = sd_frac
                            )
                        ),
                        ct
                    )
                    
                    sim_counts <- setNames(
                        as.data.table(
                            round(
                                do.call(
                                    densities,
                                    c(
                                        density_fun_args[[tp]],
                                        list(n = n)
                                    )
                                )
                            )
                        ),
                        ct
                    )
                    
                    sim_counts <- as.data.table(
                        sapply(
                            ct,
                            function(tp) {
                                if('n' %in% names(density_fun_args[[tp]])) {
                                    round(
                                        do.call(
                                            densities,
                                            density_fun_args[[tp]]
                                        )
                                    )
                                } else {
                                    round(
                                        do.call(
                                            densities,
                                            c(
                                                density_fun_args[[tp]],
                                                list(n = n)
                                            )
                                        )
                                    )
                                }
                            },
                            USE.NAMES = TRUE
                        )
                    )
                    
                    sim_counts[, proportion := p]
                    
                    sim_counts
                    
                }
            )
        )
        
        # In the following, we sample n*length(proportions_vec) = 25*40 random uniforms
        # for each cel type not in initial_cell_types, then scales them so that the
        # random normal for the initial cell type is the right proportion:
        
        initial[
            ,
            # eval(other_types) := simulate_counts(
            #     other_types,
            #     densities = densities,
            #     n = n*length(proportions_vec),
            #     density_fun_args = density_fun_args[
            #         other_types
            #     ]
            # )
            eval(other_types) := as.data.table(
                sapply(
                    other_types,
                    function(tp) {
                        round(
                            do.call(
                                densities,
                                c(
                                    density_fun_args[[tp]],
                                    list(n = n*length(proportions_vec))
                                )
                            )
                        )
                    },
                    USE.NAMES = TRUE
                )
            )
            ][
                ,
                eval(other_types) := as.list(
                    round(
                        as.numeric(.SD)*get(ct)*(1 - get('proportion'))/
                            (get('proportion')*sum(as.numeric(.SD)))
                    )
                ),
                by = row.names(initial), # .I and row.names(.SD) don't work here
                .SDcols = other_types
                ]
        
        initial
        
    },
    simplify = FALSE,
    USE.NAMES = TRUE
)$cancer[
    ,
    -'median_proportion'
    ]

# simulated_tumours_row_numbers <- simulate_tumours_row_numbers(
#     sc_data,
#     simulated_counts
# )

# Now sample rows from the single cell data of the predefined counts (just return the
# row numbers, not the actual rows):

types_var <- 'cell_type'

sc_data[
    ,
    row_number := 1:.N
    ]

setcolorder(sc_data, 'row_number')

simulated_tumours_row_numbers <- apply(
    simulated_counts,
    1,
    function(s) {
        sapply(
            names(s),
            function(tp) {
                # If we're trying to sample more rows than actually exist for a given type,
                # take all rows for that type the appropriate number of times, then sample
                # for the remaining number.  This means some rows will be chosen multiple
                # times.
                single_cell_data[
                    get(types_var) == tp,
                    switch(
                        (s[[tp]] > .N) + 1,
                        sample(row_number, s[[tp]]),
                        c(
                            rep(row_number, s[[tp]] %/% .N),
                            sample(row_number, s[[tp]] %% .N)
                        )
                    )
                    ]
            },
            simplify = FALSE,
            USE.NAMES = TRUE
        )
    }
)

# Construct bulk expression matrix for simulated tumours:

single_cell_data_for_aggregating <- 2^sc_data[
    ,
    ccle[
        gene_id %in% names(sc_data),
        gene_id
        ],
    with = FALSE
    ] - 1

# Sum over the row numbers selected for each bulk sample:

bulk_tumours <- log2(
    rbindlist(
        lapply(
            simulated_tumours_row_numbers,
            function(st) {
                as.list(
                    colSums(
                        single_cell_data_for_aggregating[unlist(st)]
                    )
                )
            }
        )
    ) + 1
)
