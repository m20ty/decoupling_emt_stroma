library(data.table)
library(ggplot2)
library(cowplot)
library(magrittr)
library(limma)

source('general_functions.R')
source('tcga_functions.R')
source('sc_functions.R')





expression_data <- fread('../../TCGA_data/tcga_expression_data.csv', key = 'id')
meta_data <- fread('../../TCGA_data/tcga_meta_data.csv', key = 'id')

expression_data <- expression_data[meta_data[sample_type != 'normal', id]]
meta_data <- meta_data[sample_type != 'normal']

emt_markers <- fread('../../emt_markers.csv')[
    ,
    gene := alias2SymbolTable(gene)
][
    source != 'GO',
    sort(unique(gene))
]





# We could select which cell types to keep and which to put under 'other' using a fraction
# of the maximum of the cell type frequencies.  I decided to select them manually instead,
# because there are certain cell types I want to show in some cases, even though thery're
# a bit rare.

single_cell_metadata <- list(

    hnsc = list(

        tcga_cancer_types = 'HNSC',

        read_quote = quote(
            fread('../data_and_figures/puram_hnscc_2017.csv')[
                ,
                c('cell_type', 'lymph_node', 'processed_by_maxima_enzyme') := .(
                    plyr::mapvalues(
                        cell_type,
                        c('', 'b_cell', 'dendritic', 'myocyte'),
                        rep('other', 4)
                    ),
                    NULL,
                    NULL
                )
            ]
        )

    ),

    lung = list(

        tcga_cancer_types = c('LUAD', 'LUSC'),

        read_quote = quote(
            fread('../data_and_figures/lambrechts_nsclc_2018.csv')[
                sample_type != 'normal',
                -c('cluster', 'disease', 'sample_type', 'annotation')
            ][
                ,
                cell_type := plyr::mapvalues(
                    cell_type,
                    c('alveolar', 'epithelial'),
                    rep('other', 2)
                )
            ]
        )

    ),

    brca = list(

        tcga_cancer_types = 'BRCA',

        # Here we only have 3 tumours with more than 50 cancer cells.

        read_quote = quote(
            fread('../data_and_figures/chung_breast_cancer_2017.csv')[
                ,
                c('cell_type', 'lymph_node') := .(
                    plyr::mapvalues(cell_type, '', 'other'),
                    NULL
                )
            ]
        )

    ),

    coadread = list(

        tcga_cancer_types = c('COAD', 'READ'),

        # Note there's only one tumour with more than 50 cancer cells in this
        # dataset, so I'm not sure we should really be using it.  Perhaps it
        # should go into a supplementary figure, without controlling for number
        # of genes detected, to see if we can see any pattern at all.

        read_quote = quote(
            fread('../data_and_figures/li_colorectal_2017.csv')[
                sample_type == 'tumour',
                -c('sample_type', 'epithelial_cluster')
            ][
                ,
                cell_type := plyr::mapvalues(
                    cell_type,
                    c('', 'endothelial', 'mast', 'epithelial'),
                    c(rep('other', 3), 'cancer')
                )
            ]
        )

    ),

    paad = list(

        tcga_cancer_types = 'PAAD',

        # In the following I previously put acinar cells under the 'cancer' umbrella, but
        # decided to remove them for two reasons: first, they seem to express epithelial
        # markers at much lower levels than the ductal cells (maybe this is obvious to a
        # real biologist); and second, I believe it is actually the ductal cells that the
        # cancer originates from (the clue is in the name - PDAC - also see
        # https://www.pancreaticcancer.org.uk/information-and-support/facts-about-pancreatic-cancer/types-of-pancreatic-cancer/,
        # where they claim that less than 1% of pancreatic cancers are acinar cell
        # cancers).  So we can assume that the acinar cells won't be malignant.  This
        # also means I'm giving the acinar cells the same treatment as the alveolar cells
        # in lung, which makes sense because the alveolar duct is similar in form to the
        # acinar duct in pancreas.

        # EDIT 15/01/2020: we'll also filter out tumours with fewer than 50 ductal cells.

        read_quote = quote(
            fread('../data_and_figures/elyada_pdac_2019.csv')[
                patient %in% c('HT137', 'HT149', 'HN149', 'HT99', 'HN150')
            ][
                ,
                c('cell_type', 'cluster') := .(
                    plyr::mapvalues(
                        cell_type,
                        c(
                            'ductal',
                            'acinar',
                            'perivascular',
                            'dendritic',
                            'endothelial',
                            'mast'
                            # 'red_blood_cell' # None of these left are removing patients
                        ),
                        c(
                            'cancer',
                            rep('other', 5)
                        )
                    ),
                    NULL
                )
            ]
        )

    ),

    lihc = list(

        tcga_cancer_types = 'LIHC',

        # Here we remove the tumours with fewer than 50 cancer cells.  Also, there are
        # surprisingly few B and T cells, and I don't want to include "HPC-like" cells
        # as a separate type in the lineplots, so I'll put all these under 'other'.

        read_quote = quote(
            fread('../data_and_figures/ma_liver_2019.csv')[
                patient %in% c('C26', 'C25', 'H38', 'H37', 'C56', 'C46', 'H65', 'C66')
            ][
                ,
                c('cell_type', 'sample', 'disease') := .(
                    plyr::mapvalues(
                        cell_type,
                        c('b_cell', 'hpc-like', 't_cell', 'unclassified'),
                        rep('other', 4)
                    ),
                    NULL,
                    NULL
                )
            ]
        )

    ),

    tnbc = list(

        tcga_cancer_types = 'BRCA',

        # Remove the tumour with fewer than 50 epithelial cells:

        read_quote = quote(
            fread('../data_and_figures/karaayvaz_tnbc_2018.csv')[
                patient != 'PT058'
            ][
                ,
                cell_type := plyr::mapvalues(
                    cell_type,
                    c('epithelial', '', 'Bcell', 'endothelial', 'stroma', 'Tcell'),
                    c('cancer', 'other', 'other', 'other', 'fibroblast', 't_cell')
                )
            ]
        )

    )

)





# See single_cell_analysis.R for the definition of genes_list.

genes_list <- readRDS('../data_and_figures/sc_genes_list.rds')





# Check that the EMT-like cells in TNBC don't look like doublets:

sc_data <- eval(single_cell_metadata$tnbc$read_quote)
setkey(sc_data, id)
sc_data <- sc_data[sc_cancer_fibroblast$tnbc$cells_filtered]
ordering_cancer_tnbc <- sc_cancer_fibroblast$tnbc$ordering_cells$cancer

all_genes_filtered_tnbc <- sc_data[
    ,
    names(.SD)[apply(.SD, 2, sc_cancer_fibroblast_args$tnbc$genes_filter_fun)],
    .SDcols = -c('id', 'patient', 'cell_type')
]

set.seed(74)

# Take difference between average expression of genes in EMT-like cells and that in a
# subset of "normal" cancer cells which don't show such EMT characteristics:

norm_emt <- sc_data[cell_type == 'cancer'][
    ordering_cancer_tnbc,
    colMeans(.SD[1:70]) - colMeans(.SD[sample(1:200, 70) + 200]),#(.N/2 - 34):(.N/2 + 35)
    .SDcols = all_genes_filtered_tnbc
]

# Similarly, take difference between fibroblasts and "normal" cancer cells:

norm_fibroblast <- sc_data[cell_type == 'fibroblast'][
    ,
    colMeans(.SD),
    .SDcols = all_genes_filtered_tnbc
] - sc_data[cell_type == 'cancer'][
    ordering_cancer_tnbc,
    colMeans(.SD[sample(1:200, 70) + 200]),#(.N/2 - 34):(.N/2 + 35)])
    .SDcols = all_genes_filtered_tnbc
]

# Simulate doublets by adding profiles of fibroblasts to those of "normal" cancer cells":

fib_inds <- sample(1:length(sc_cancer_fibroblast$tnbc$ordering_cells$fibroblast))

cancer_inds <- sample(1:200, round(length(fib_inds)/2)) + 200

doublet_indices <- mapply(
    c,
    fib_inds[1:round(length(fib_inds)/2)],
    cancer_inds
)

sim_doublets <- apply(
    doublet_indices,
    2,
    function(inds) {
        sc_data[
            ,
            log2(
                2^(.SD[cell_type == 'fibroblast'][inds[1]])
                + 2^(.SD[cell_type == 'cancer'][ordering_cancer_tnbc][inds[2]])
                - 1
            ),
            .SDcols = all_genes_filtered_tnbc
        ]
    }
) %>% rbindlist

# Normalise the simulated doublets by the "normal" cancer cells:

norm_doublets <- colMeans(sim_doublets) - sc_data[cell_type == 'cancer'][
    ordering_cancer_tnbc,
    colMeans(.SD[sample((1:200)[-(cancer_inds - 200)], round(length(fib_inds)/2)) + 200]),
    .SDcols = all_genes_filtered_tnbc
]

# Normalise the fibroblasts not used in the doublet simulations:

norm_fibroblasts_not_for_doublets <- sc_data[cell_type == 'fibroblast'][
    fib_inds[-(1:round(length(fib_inds)/2))],
    colMeans(.SD),
    .SDcols = all_genes_filtered_tnbc
] - sc_data[cell_type == 'cancer'][
    ordering_cancer_tnbc,
    colMeans(.SD[sample((1:200)[-(cancer_inds - 200)], round(length(fib_inds)/2)) + 200]),
    .SDcols = all_genes_filtered_tnbc
]

# Make plots:

pdf('../data_and_figures/emt_program_tnbc_doublets.pdf', width = 12, height = 4)

plot_grid(
    plotlist = lapply(
        list(
            list(
                norm_fibroblast,
                norm_emt,
                names = c('Fibroblasts', 'EMT cells')
            ),
            list(
                norm_fibroblasts_not_for_doublets,
                norm_doublets,
                names = c('Fibroblasts', 'Simualted doublets')
            ),
            list(
                norm_doublets,
                norm_emt,
                names = c('Simulated doublets', 'EMT cells')
            )
        ),
        function(vp) {
            qplot(
                vp[[1]],
                vp[[2]],
                xlab = vp$names[1],
                ylab = vp$names[2]
            ) +
                geom_smooth(method = 'lm', se = FALSE) +
                geom_hline(yintercept = 0, colour = 'grey', linetype = 'dashed') +
                geom_vline(xintercept = 0, colour = 'grey', linetype = 'dashed') +
                annotate(
                    'text',
                    x = 7,
                    y = -1,
                    label = latex2exp::TeX(
                        sprintf(
                            '$\\rho = %.2f$',
                            cor(
                                vp[[1]],
                                vp[[2]],
                                method = 'spearman'
                            )
                        )
                    )
                ) +
                theme_test()
        }
    ),
    nrow = 1,
    ncol = 3
)

dev.off()

# To me this suggests it is unlikely that the EMT cells are doublets.  It shows that the
# EMT cells correlate significantly but not very strongly with fibroblasts (0.58), and
# this correlation is the same as with the simulated doublets (0.59).  Meanwhile, the
# simulated doublets correlate highly with the fibroblasts (0.89), much more highly than
# the EMT cells do.  I think that if the EMT cells were doublets then their correlation
# with the simulated doublets would be somewhere between 0.59 and 0.89.

# Try with all cancer types, for a comparison of the correlation values:

sim_doublets_data <- sapply(
    c('hnsc', 'lung', 'paad', 'liver', 'tnbc'),
    function(ct) {

        cat(paste0(ct, '...'))

        sc_data <- eval(single_cell_metadata[[ct]]$read_quote)
        setkey(sc_data, id)
        sc_data <- sc_data[sc_cancer_fibroblast[[ct]]$cells_filtered]
        ordering_cancer <- sc_cancer_fibroblast[[ct]]$ordering_cells$cancer
        ordering_fibroblast <- sc_cancer_fibroblast[[ct]]$ordering_cells$fibroblast

        all_genes_filtered <- sc_data[
            ,
            names(.SD)[apply(.SD, 2, sc_cancer_fibroblast_args[[ct]]$genes_filter_fun)],
            .SDcols = -c('id', 'patient', 'cell_type')
        ]

        set.seed(sc_cancer_fibroblast_args[[ct]]$seed)

        # Take difference between average expression of genes in EMT-like cells and that in a
        # subset of "normal" cancer cells which don't show such EMT characteristics:

        norm_emt <- sc_data[cell_type == 'cancer'][
            ordering_cancer,
            colMeans(.SD[1:70]) - colMeans(.SD[sample(1:200, 70) + .N/2 - 99]),
            .SDcols = all_genes_filtered
        ]

        # Similarly, take difference between fibroblasts and "normal" cancer cells:

        norm_fibroblast <- sc_data[cell_type == 'fibroblast'][
            ,
            colMeans(.SD),
            .SDcols = all_genes_filtered
        ] - sc_data[cell_type == 'cancer'][
            ordering_cancer,
            colMeans(.SD[sample(1:200, 70) + .N/2 - 99]),
            .SDcols = all_genes_filtered
        ]

        # Simulate doublets by adding profiles of fibroblasts to those of "normal" cancer cells":

        fib_inds <- sample(1:length(ordering_fibroblast), min(100, length(ordering_fibroblast)))

        cancer_inds <- sample(1:200, round(length(fib_inds)/2)) + length(ordering_cancer)/2 - 99

        doublet_inds <- mapply(
            c,
            fib_inds[1:round(length(fib_inds)/2)],
            cancer_inds
        )

        sim_doublets <- apply(
            doublet_inds,
            2,
            function(inds) {
                sc_data[
                    ,
                    log2(
                        2^(.SD[cell_type == 'fibroblast'][inds[1]])
                        + 2^(.SD[cell_type == 'cancer'][ordering_cancer][inds[2]])
                        - 1
                    ),
                    .SDcols = all_genes_filtered
                ]
            }
        ) %>% rbindlist

        # Normalise the simulated doublets by the "normal" cancer cells:

        norm_doublets <- colMeans(sim_doublets) - sc_data[
            cell_type == 'cancer'
        ][
            ordering_cancer,
            colMeans(
                .SD[
                    sample(
                        (1:200)[-(cancer_inds - length(ordering_cancer)/2 + 99)],
                        round(length(fib_inds)/2)
                    ) + length(ordering_cancer)/2 - 99
                ]
            ),
            .SDcols = all_genes_filtered
        ]

        # Normalise the fibroblasts not used in the doublet simulations:

        norm_fibroblasts_not_for_doublets <- sc_data[cell_type == 'fibroblast'][
            fib_inds[-(1:round(length(fib_inds)/2))],
            colMeans(.SD),
            .SDcols = all_genes_filtered
        ] - sc_data[cell_type == 'cancer'][
            ordering_cancer,
            colMeans(
                .SD[
                    sample(
                        (1:200)[-(cancer_inds - length(ordering_cancer)/2 + 99)],
                        round(length(fib_inds)/2)
                    ) + length(ordering_cancer)/2 - 99
                ]
            ),
            .SDcols = all_genes_filtered
        ]

        cat('Done!\n')

        list(
            norm_emt = norm_emt,
            norm_fibroblast = norm_fibroblast,
            norm_doublets = norm_doublets,
            norm_fibroblasts_not_for_doublets = norm_fibroblasts_not_for_doublets,
            fibroblast_indices = fib_inds,
            doublet_indices = doublet_inds
        )

    },
    simplify = FALSE,
    USE.NAMES = TRUE
)

plot_grid(
    plotlist = lapply(
        sim_doublets_data,
        function(li) {

            plot_grid(
                plotlist = lapply(
                    list(
                        list(
                            li$norm_fibroblast,
                            li$norm_emt,
                            names = c('Fibroblasts', 'EMT cells')
                        ),
                        list(
                            li$norm_fibroblasts_not_for_doublets,
                            li$norm_doublets,
                            names = c('Fibroblasts', 'Simualted doublets')
                        ),
                        list(
                            li$norm_doublets,
                            li$norm_emt,
                            names = c('Simulated doublets', 'EMT cells')
                        )
                    ),
                    function(vp) {
                        qplot(
                            vp[[1]],
                            vp[[2]],
                            xlab = vp$names[1],
                            ylab = vp$names[2]
                        ) +
                            geom_smooth(method = 'lm', se = FALSE) +
                            geom_hline(yintercept = 0, colour = 'grey', linetype = 'dashed') +
                            geom_vline(xintercept = 0, colour = 'grey', linetype = 'dashed') +
                            annotate(
                                'text',
                                x = 7,
                                y = -1,
                                label = latex2exp::TeX(
                                    sprintf(
                                        '$\\rho = %.2f$',
                                        cor(
                                            vp[[1]],
                                            vp[[2]],
                                            method = 'spearman'
                                        )
                                    )
                                )
                            ) +
                            theme_test()
                    }
                ),
                nrow = 1,
                ncol = 3
            )

        }
    ),
    nrow = length(sim_doublets_data)
)

# Just the EMT vs. fibroblast plots, without the simulated doublets stuff:

pdf('../data_and_figures/emt_program_doublets.pdf', width = 17, height = 8)

plot_grid(
    plotlist = lapply(
        names(sim_doublets_data),
        function(ct) {

            qplot(
                sim_doublets_data[[ct]]$norm_fibroblast,
                sim_doublets_data[[ct]]$norm_emt,
                xlab = 'Fibroblasts',
                ylab = 'EMT cells',
                main = sc_cancer_fibroblast_heatmaps_args[[ct]]$annotations_title
            ) +
                geom_smooth(method = 'lm', se = FALSE) +
                geom_hline(yintercept = 0, colour = 'grey', linetype = 'dashed') +
                geom_vline(xintercept = 0, colour = 'grey', linetype = 'dashed') +
                annotate(
                    'text',
                    x = floor(max(sim_doublets_data[[ct]]$norm_fibroblast)) - 1,
                    y = ceiling(min(sim_doublets_data[[ct]]$norm_emt)) + 0.5,
                    label = latex2exp::TeX(
                        sprintf(
                            '$\\rho = %.2f$',
                            cor(
                                sim_doublets_data[[ct]]$norm_fibroblast,
                                sim_doublets_data[[ct]]$norm_emt,
                                method = 'spearman'
                            )
                        )
                    )
                ) +
                theme_test()

        }
    ),
    nrow = 2,
    ncol = 4
)

dev.off()

# Itay: This doesn't make sense to me, so I suspect that the analysis you have done may
# differ from what I imagined: With X being the cancer cell average, if (EMT-X) is
# correlated with (CAF-X) then it must correlate more highly with
# ([EMT+CAF]/2-X)=(EMT-X)/2 + (CAF-X)/2. I do understand that my notation may not be
# exactly correct because the doublets are simulated one at a time rather than by
# combining their averages, but still, I would expect higher correlation with the
# simulated doublets rather than lower correlations...

# Me: My method differed slightly from this in that the doublets were not simulated as
# "EMT + CAF" but rather "normal cancer cell + CAF".  That is, I took for the doublets
# a sample of "ordinary" cancer cells which didn't have high EMT scores.  So, in your
# notation, I think I'm comparing cor(EMT-X, CAF-X) with cor(EMT - X, [X + CAF]/2 - X)...

# Shit!  This equals cor(EMT-X, [CAF-X]/2)!  So it's even weirder that there is worse
# correlation, except that, as Itay said, in reality you're not combining averages.
# Still, maybe something is wrong with the method.





# Heatmaps of EMT genes for just cancer cells:

sc_cancer_args <- list(

    hnsc = list(
        seed = 8511,
        genes_filter_fun = function(x) {
            log2(mean(2^x - 1) + 1) >= 3 |
                sum(x >= 7) >= length(x)/100
        },
        scores_filter_fun = function(x) {
            log2(mean(2^x - 1) + 1) >= 4 |
                sum(x >= 7) >= length(x)/100
        }
    ),

    lung = list(
        seed = 2566,
        genes_filter_fun = function(x) {
            log2(mean(2^x - 1) + 1) >= 3 |
                sum(x >= 7) >= length(x)/100
        },
        scores_filter_fun = function(x) {
            log2(mean(2^x - 1) + 1) >= 4 |
                sum(x >= 7) >= length(x)/100
        }
    ),

    brca = list(
        seed = 3718,
        genes_filter_fun = function(x) {
            log2(mean(2^x - 1) + 1) >= 2.5 |
                sum(x >= 7) >= length(x)/100
        },
        scores_filter_fun = function(x) {
            log2(mean(2^x - 1) + 1) >= 3.5 |
                sum(x >= 7) >= length(x)/100
        }
    ),

    coadread = list(
        seed = 3361,
        genes_filter_fun = function(x) {
            log2(mean(2^x - 1) + 1) >= 2.5 |
                sum(x >= 7) >= length(x)/100
        },
        scores_filter_fun = function(x) {
            log2(mean(2^x - 1) + 1) >= 3.5 |
                sum(x >= 7) >= length(x)/100
        }
    ),

    paad = list(
        seed = 5368,
        genes_filter_fun = function(x) {
            log2(mean(2^x - 1) + 1) >= 3 |
                sum(x >= 7) >= length(x)/100
        },
        scores_filter_fun = function(x) {
            log2(mean(2^x - 1) + 1) >= 4 |
                sum(x >= 7) >= length(x)/100
        }
    ),

    lihc = list(
        seed = 4376,
        genes_filter_fun = function(x) {
            log2(mean(2^x - 1) + 1) >= 2.5 |
                sum(x >= 7) >= length(x)/100
        },
        scores_filter_fun = function(x) {
            log2(mean(2^x - 1) + 1) >= 3.5 |
                sum(x >= 7) >= length(x)/100
        }
    ),

    tnbc = list(
        seed = 456,
        genes_filter_fun = function(x) {
            log2(mean(2^x - 1) + 1) >= 3 |
                sum(x >= 7) >= length(x)/100
        },
        scores_filter_fun = function(x) {
            log2(mean(2^x - 1) + 1) >= 4 |
                sum(x >= 7) >= length(x)/100
        }
    )

)

# sc_cancer <- readRDS('../data_and_figures/sc_cancer.rds')

sc_cancer <- sapply(
    names(single_cell_metadata),
    function(ct) {

        cat(paste0(ct, '...'))

        sc_data <- eval(single_cell_metadata[[ct]]$read_quote)

        set.seed(sc_cancer_args[[ct]]$seed)

        out <- do.call(
            sc_groups,
            args = c(
                list(
                    genes = genes_list[[ct]]$unfiltered,
                    sc_data = sc_data[cell_type == 'cancer'],
                    groups = 'cancer',
                    to_keep = character(0),
                    min_sig_size = 0
                ),
                sc_cancer_args[[ct]][-1]
            )
        )

        cat('Done!\n')

        out

    },
    simplify = FALSE,
    USE.NAMES = TRUE
)

saveRDS(sc_cancer, '../data_and_figures/sc_cancer.rds')

sc_cancer_heatmaps_args <- list(

    hnsc = list(
        annotations = c(
            'CD44',
            'FN1',
            'IL32',
            'ITGB1',
            'LAMC2',
            'MMP1',
            'PLAUR',
            'PPIB',
            'SAT1',
            'SDC1',
            'SNAI2',
            'TGFBI',
            'TNC',
            'VIM'
        ),
        annotations_title = 'Head and Neck'
    ),

    lung = list(
        annotations = c(
            'AREG',
            'CAPG',
            'CD44',
            'IGFBP3',
            'IL32',
            'JUN',
            'LAMC2',
            'PLAUR',
            'QSOX1',
            'RHOB',
            'SAT1',
            'SDC4',
            'SPP1',
            'TPM4'
        ),
        annotations_side = 'right',
        annotations_title = 'Lung'
    ),

    brca = list(
        annotations = c(
            'AREG',
            'CALU',
            'CD44',
            'CRYAB',
            'CXCL1',
            'ITGB1',
            'PPIB',
            'QKI',
            'RHOB',
            'SAT1',
            'SDC1',
            'SDC4',
            'TIMP3',
            'TPM1'
        ),
        annotations_title = 'Breast'
    ),

    coadread = list(
        annotations = c(
            'AREG',
            'CD44',
            'CXCL1',
            'DST',
            'IL32',
            'ITGB1',
            'LAMA3',
            'PLAUR',
            'PLOD3',
            'PPIB',
            'SAT1',
            'TGFBI',
            'TIMP1',
            'TPM4'
        ),
        annotations_side = 'right',
        annotations_title = 'Colorectal'
    ),

    paad = list(
        annotations = c(
            'CD44',
            'CTGF',
            'CXCL1',
            'IGFBP7',
            'IL32',
            'ITGB1',
            'JUN',
            'PPIB',
            'SAT1',
            'SDC4',
            'SPP1',
            'TFPI2',
            'TPM4',
            'VIM'
        ),
        annotations_title = 'Pancreatic'
    ),

    lihc = list(
        annotations = c(
            'CD44',
            'CXCL1',
            'CXCL8',
            'FN1',
            'IGFBP3',
            'IL32',
            'ITGB1',
            'LAMC2',
            'PLOD2',
            'SERPINE1',
            'SPP1',
            'TPM1',
            'VEGFA',
            'VIM'
        ),
        annotations_side = 'right',
        annotations_title = 'Liver'
    ),

    tnbc = list(
        annotations = c(
            'ACTA2',
            'AXL',
            'CD44',
            'COL3A1',
            'CYR61',
            'FN1',
            'MGP',
            'POSTN',
            'SERPING1',
            'SPARC',
            'SPARCL1',
            'TAGLN',
            'TIMP1',
            'VIM'
        ),
        annotations_title = 'Breast (TNBC)'
    )

)

sc_cancer_heatmaps <- sapply(
    names(single_cell_metadata),
    function(ct) {
        set.seed(sc_cancer_args[[ct]]$seed)
        do.call(
            sc_groups_heatmap,
            args = c(
                list(
                    sc_groups_list = sc_cancer[[ct]],
                    groups = 'cancer',
                    default_figure_widths = list(
                        annotations = 2,
                        cancer = 7.2
                    ),
                    figure_spacing = 2.5,
                    annotations_nudge = 0.25,
                    es_fun = NULL,
                    es_title = 'EMT score'
                ),
                sc_cancer_heatmaps_args[[ct]]
            )
        )
    },
    simplify = FALSE,
    USE.NAMES = TRUE
)





# Make PDFs:

cairo_pdf(
    '../data_and_figures/single_cell_heatmaps_cancer.pdf',
    width = 10,
    height = 10
)

plot_grid(
    plot_grid(
        cowplot_sc(
            sc_cancer_heatmaps$hnsc,
            legend_space = 0,
            heights = c(1.5, 20, 4),
            es_x_axis_title_vjust = 1.3,
            es_y_axis_title_angle = 0,
            es_y_axis_title_xpos = 0.8
        ),
        blank_plot(),
        cowplot_sc(
            sc_cancer_heatmaps$lung,
            legend_space = 0,
            heights = c(1.5, 20, 4),
            es_x_axis_title_vjust = 1.3,
            es_y_axis_title_angle = 0,
            es_y_axis_title_xpos = 0.2
        ),
        ncol = 3,
        nrow = 1,
        rel_widths = c(20, 1, 20)
    ),
    blank_plot(),
    plot_grid(
        cowplot_sc(
            sc_cancer_heatmaps$paad,
            legend_space = 0,
            heights = c(1.5, 20, 4),
            es_x_axis_title_vjust = 1.3,
            es_y_axis_title_angle = 0,
            es_y_axis_title_xpos = 0.8
        ),
        blank_plot(),
        cowplot_sc(
            sc_cancer_heatmaps$lihc,
            legend_space = 0,
            heights = c(1.5, 20, 4),
            es_x_axis_title_vjust = 1.3,
            es_y_axis_title_angle = 0,
            es_y_axis_title_xpos = 0.2
        ),
        ncol = 3,
        nrow = 1,
        rel_widths = c(20, 1, 20)
    ),
    blank_plot(),
    plot_grid(
        cowplot_sc(
            sc_cancer_heatmaps$tnbc,
            legend_space = 0,
            heights = c(1.5, 20, 4),
            es_x_axis_title_vjust = 1.3,
            es_y_axis_title_angle = 0,
            es_y_axis_title_xpos = 0.8
        ),
        blank_plot(),
        sc_cancer_heatmaps$paad$plots$heatmap_legend,
        sc_cancer_heatmaps$paad$plots$genes_detected_legend,
        blank_plot(),
        ncol = 5,
        nrow = 1,
        rel_widths = c(20, 1, 6, 4, 10)
    ),
    blank_plot(),
    ncol = 1,
    nrow = 6,
    rel_heights = c(20, 1, 20, 1, 20, 1)
)

dev.off()

# Just the Chung breast and Li colorectal data:

brca_cancer <- sc_cancer_heatmaps$brca
coadread_cancer <- sc_cancer_heatmaps$coadread

brca_cancer$plots$genes_detected$cancer <- brca_cancer$plots$genes_detected$cancer +
    labs(title = 'Breast')
brca_cancer$plots$annotations <- brca_cancer$plots$annotations +
    theme(axis.title = element_blank())

coadread_cancer$plots$genes_detected$cancer <- coadread_cancer$plots$genes_detected$cancer +
    labs(title = 'Colorectal')
coadread_cancer$plots$annotations <- coadread_cancer$plots$annotations +
    theme(axis.title = element_blank())

cairo_pdf(
    '../data_and_figures/single_cell_heatmaps_cancer_supp.pdf',
    width = 11,
    height = 3.5,
    onefile = TRUE
)

plot_grid(
    plot_grid(
        cowplot_sc(
            brca_cancer,
            legend_space = 0,
            heights = c(3.5, 20, 4),
            es_x_axis_title_vjust = 1.3,
            es_y_axis_title_angle = 0,
            es_y_axis_title_xpos = 0.8
        ),
        blank_plot(),
        cowplot_sc(
            coadread_cancer,
            legend_space = 0,
            heights = c(3.5, 20, 4),
            es_x_axis_title_vjust = 1.3,
            es_y_axis_title_angle = 0,
            es_y_axis_title_xpos = 0.2
        ),
        plot_grid(
            brca_cancer$plots$heatmap_legend,
            coadread_cancer$plots$genes_detected_legend,
            ncol = 1,
            nrow = 2
        ),
        ncol = 4,
        nrow = 1,
        rel_widths = c(20, 1, 20, 6)
    ),
    blank_plot(),
    ncol = 1,
    nrow = 2,
    rel_heights = c(20, 1)
)

dev.off()





# Mini un-annotated plots (not for head and neck):

titles <- list(
    lung = 'Lung',
    brca = 'Breast',
    coadread = 'Colorectal',
    paad = 'Pancreatic'
)

cairo_pdf(
    '../data_and_figures/single_cell_heatmaps_mini_examples.pdf',
    width = 16,
    height = 3.5,
    onefile = TRUE
)

plot_grid(
    plotlist = lapply(
        names(sc_cancer_fibroblast_heatmaps[c('lung', 'brca', 'coadread', 'paad')]),
        function(ct) {
            with(
                sc_cancer_fibroblast_heatmaps[[ct]],
                plot_grid(
                    plots$genes_detected$cancer +
                        labs(title = titles[[ct]]) +
                        theme(
                            plot.margin = unit(c(2.5, 1, 1.5, 4), 'pt'),
                            plot.title = element_text(
                                margin = unit(c(1, 5.5, 2, 5.5), 'pt'),
                                size = 24
                            )
                        ),
                    plots$genes_detected$fibroblast +
                        labs(title = '') +
                        theme(
                            plot.margin = unit(c(2.5, 4, 1, 1), 'pt'),
                            plot.title = element_text(
                                margin = unit(c(1, 5.5, 2, 5.5), 'pt'),
                                size = 24
                            )
                        ),
                    plots$heatmaps$cancer +
                        theme(plot.margin = unit(c(1, 1, 1, 4), 'pt')),
                    plots$heatmaps$fibroblast +
                        theme(plot.margin = unit(c(1, 4, 1, 1.5), 'pt')),
                    plots$expression_summary$cancer +
                        theme(
                            axis.text.y.left = element_blank(),
                            axis.ticks.length = unit(0, 'pt'),
                            axis.title = element_blank(),
                            plot.margin = unit(c(1, 1, 5.5, 4), 'pt')
                        ),
                    plots$expression_summary$fibroblast +
                        theme(
                            axis.text.y.right = element_blank(),
                            axis.ticks.length = unit(0, 'pt'),
                            axis.title = element_blank(),
                            plot.margin = unit(c(1, 4, 5.5, 1), 'pt')
                        ),
                    ncol = 2,
                    nrow = 3,
                    rel_widths = unlist(figure_widths[c('cancer', 'fibroblast')]),
                    rel_heights = c(1, 4, 1)
                )
            )
        }
    ),
    ncol = 4,
    nrow = 1
)

dev.off()





# Distinguishing rare from shared EMT:

# If already done, skip by reading from RDS:

# rare_shared_emt <- readRDS('../data_and_figures/rare_shared_emt.rds')

rare_shared_emt <- sapply(
    names(single_cell_metadata),
    function(ct) {

        cat(paste0(ct, '...'))

        sc_data <- eval(single_cell_metadata[[ct]]$read_quote)

        score_dist_data <- score_dist(
            sc_data,
            sc_cancer[[ct]]$cells_filtered,
            sc_cancer_args[[ct]]$genes_filter_fun,
            sc_cancer[[ct]]$scores,
            seed = sc_cancer_fibroblast_args[[ct]]$seed,
            plot_title = ct,
            null_cor_threshold = 0.99
        )

        out <- rare_vs_shared_emt(
            sc_data,
            genes_list[[ct]]$unfiltered,
            score_dist_data,
            args_list = sc_cancer_args[[ct]],
            scores_filter_fun_shared = sc_cancer_fibroblast_args[[ct]]$scores_filter_fun,
            seed = sc_cancer_args[[ct]]$seed,
            shared_threshold_quantile = 0.3,
            rare_threshold_quantile = 0.7
        )

        cat('Done!\n')

        out

    },
    simplify = FALSE,
    USE.NAMES = TRUE
)

saveRDS(rare_shared_emt, '../data_and_figures/rare_shared_emt.rds')





# Density plots:

pdf(
    '../data_and_figures/rare_emt_densities.pdf',
    width = 11,
    height = 7
)

plot_grid(
    plotlist = c(
        lapply(
            rare_shared_emt,
            function(li) {
                li$score_dist_data$density_plot +
                    theme(legend.position = 'none')
            }
        ),
        list(get_legend(rare_shared_emt[[1]]$score_dist_data$density_plot))
    ),
    ncol = 3,
    nrow = 2
)

dev.off()





sc_cancer_z_score_heatmaps_args_list <- list(

    hnsc = list(
        annotations = c(
            'CAPG',
            'CD44',
            'FN1',
            'ITGB1',
            'LAMC2',
            'MMP1',
            'PPIB',
            'SAT1',
            'SERPINE1',
            'SNAI2',
            'TGFBI',
            'TNC',
            'TPM4',
            'VIM'
        ),
        annotations_side = 'left',
        annotations_title = 'Head and Neck'
    ),

    lung = list(
        annotations = c(
            'AREG',
            'CAPG',
            'CD44',
            'IGFBP3',
            'IL32',
            'JUN',
            'LAMC2',
            'PPIB',
            'QSOX1',
            'RHOB',
            'SDC4',
            'SPP1',
            'TPM1',
            'VEGFA'
        ),
        annotations_side = 'right',
        annotations_title = 'Lung'
    ),

    brca = list(
        annotations = c(
            'AREG',
            'CALU',
            'CD44',
            'CRYAB',
            'CXCL1',
            'ITGB1',
            'PPIB',
            'QKI',
            'RHOB',
            'SAT1',
            'SDC1',
            'SDC4',
            'TIMP3',
            'TPM1'
        ),
        annotations_side = 'left',
        annotations_title = 'Breast'
    ),

    coadread = list(
        annotations = c(
            'AREG',
            'CD44',
            'CTGF',
            'CXCL1',
            'DST',
            'IL32',
            'ITGB1',
            'LAMA3',
            'PLAUR',
            'PLOD3',
            'PPIB',
            'SAT1',
            'TGFBI',
            'TIMP1',
            'TPM4'
        ),
        annotations_side = 'right',
        annotations_title = 'Colorectal'
    ),

    paad = list(
        annotations = c(
            'CCL2',
            'CTGF',
            'CXCL1',
            'CXCL6',
            'IGFBP7',
            'IL32',
            'ITGB1',
            'JUN',
            'PPIB',
            'SAT1',
            'SPP1',
            'TFPI2',
            'TPM4',
            'VIM'
        ),
        annotations_side = 'left',
        annotations_title = 'Pancreas'
    ),

    lihc = list(
        annotations = c(
            'SNAI1',
            'SNAI2',
            'TWIST1',
            'VIM',
            'ZEB1',
            'ZEB2'
        ),
        annotations_side = 'right',
        annotations_title = 'Liver'
    )

)

sc_cancer_heatmaps_z_score <- sapply(
    names(single_cell_metadata),
    function(ct) {

        set.seed(sc_cancer_args[[ct]]$seed)

        do.call(
            sc_groups_heatmap,
            args = c(
                list(
                    sc_groups_list = rare_shared_emt[[ct]]$sc_groups_list,
                    groups = 'cancer',
                    default_figure_widths = list(
                        annotations = 2,
                        cancer = 8
                    ),
                    figure_spacing = 2.5,
                    es_fun = NULL,
                    es_title = 'EMT score',
                    h_colours = rev(colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))(50)),
                    h_limits = c(-2, 2),
                    h_legend_breaks = c(-2, -1, 0, 1, 2),
                    h_legend_labels = c('-2' = -2, '-1' = -1, '0' = 0, '1' = 1, '2' = 2),
                    h_legend_title = 'Z score'
                ),
                sc_cancer_z_score_heatmaps_args_list[[ct]]
            )
        )

    },
    simplify = FALSE,
    USE.NAMES = TRUE
)

# Make the figure:

cairo_pdf(
    '../data_and_figures/single_cell_heatmaps_z_score.pdf',
    width = 10,
    height = 10
)

plot_grid(
    plot_grid(
        cowplot_sc(
            sc_cancer_heatmaps_z_score$hnsc,
            legend_space = 0,
            heights = c(1.5, 20, 4),
            es_x_axis_title_vjust = 1.3,
            es_y_axis_title_angle = 0,
            es_y_axis_title_xpos = 0.8
        ),
        blank_plot(),
        cowplot_sc(
            sc_cancer_heatmaps_z_score$lung,
            legend_space = 0,
            heights = c(1.5, 20, 4),
            es_x_axis_title_vjust = 1.3,
            es_y_axis_title_angle = 0,
            es_y_axis_title_xpos = 0.2
        ),
        ncol = 3,
        nrow = 1,
        rel_widths = c(20, 1, 20)
    ),
    blank_plot(),
    plot_grid(
        cowplot_sc(
            sc_cancer_heatmaps_z_score$brca,
            legend_space = 0,
            heights = c(1.5, 20, 4),
            es_x_axis_title_vjust = 1.3,
            es_y_axis_title_angle = 0,
            es_y_axis_title_xpos = 0.8
        ),
        blank_plot(),
        cowplot_sc(
            sc_cancer_heatmaps_z_score$coadread,
            legend_space = 0,
            heights = c(1.5, 20, 4),
            es_x_axis_title_vjust = 1.3,
            es_y_axis_title_angle = 0,
            es_y_axis_title_xpos = 0.2
        ),
        ncol = 3,
        nrow = 1,
        rel_widths = c(20, 1, 20)
    ),
    blank_plot(),
    plot_grid(
        cowplot_sc(
            sc_cancer_heatmaps_z_score$paad,
            legend_space = 0,
            heights = c(1.5, 20, 4),
            es_x_axis_title_vjust = 1.3,
            es_y_axis_title_angle = 0,
            es_y_axis_title_xpos = 0.8
        ),
        blank_plot(),
        sc_cancer_heatmaps_z_score$paad$plots$heatmap_legend,
        sc_cancer_heatmaps_z_score$paad$plots$genes_detected_legend,
        blank_plot(),
        ncol = 5,
        nrow = 1,
        rel_widths = c(20, 1, 6, 4, 10)
    ),
    blank_plot(),
    ncol = 1,
    nrow = 6,
    rel_heights = c(20, 1, 20, 1, 20, 1)
)

dev.off()

# I think this better illustrates the problem I was having with some of my rare EMT gene
# analysis.  E.g. in pancreas, ITGB1 has a relatively high correlation with the EMT scores,
# but it is expressed in a lot of cells, whereas VIM ranks lower in correlation with EMT
# scores but is relatively more highly expressed in fewer cells, so may better reflect
# truly "rare" EMT.





# Alternative, only using rare and shared EMT genes and restricting to head and neck, lung
# and pancreas (since breast and colorectal don't have a strong signal):

sc_cancer_heatmaps_z_score_restricted <- sapply(
    c('hnsc', 'lung', 'paad'),
    function(ct) {

        htmps <- sapply(
            c('rare', 'shared'),
            function(emt_type) {

                sc_groups_list <- rare_shared_emt[[ct]]$sc_groups_list

                sc_groups_list$genes_filtered <- rare_shared_emt[[ct]]$rare_shared_emt_data[[
                    paste0(emt_type, '_emt_genes')
                    ]]

                sc_groups_list$ordering_genes <- length(sc_groups_list$genes_filtered):1

                sc_groups_heatmap(
                    sc_groups_list,
                    groups = 'cancer',
                    default_figure_widths = list(
                        annotations = 2,
                        cancer = 8
                    ),
                    annotations_nudge = 0.25,
                    figure_spacing = 2.5,
                    es_fun = NULL,
                    es_title = 'EMT score',
                    h_colours = rev(colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))(50)),
                    h_limits = c(-2, 2),
                    h_legend_breaks = c(-2, -1, 0, 1, 2),
                    h_legend_labels = c('-2' = -2, '-1' = -1, '0' = 0, '1' = 1, '2' = 2),
                    h_legend_title = 'Z score',
                    annotations = sc_cancer_z_score_heatmaps_args_list[[ct]]$annotations,
                    annotations_title = sc_cancer_z_score_heatmaps_args_list[[ct]]$annotations_title,
                    annotations_side = sc_cancer_z_score_heatmaps_args_list[[ct]]$annotations_side
                )

            },
            simplify = FALSE,
            USE.NAMES = TRUE
        )

        c(
            list(
                plots = list(
                    genes_detected = htmps[[1]]$plots$genes_detected$cancer,
                    genes_detected_legend = htmps[[1]]$plots$genes_detected_legend,
                    expression_summary = htmps[[1]]$plots$expression_summary$cancer,
                    heatmaps = list(
                        rare = htmps$rare$plots$heatmaps$cancer,
                        shared = htmps$shared$plots$heatmaps$cancer
                    ),
                    heatmap_legend = htmps[[1]]$plots$heatmap_legend,
                    annotations = list(
                        rare = htmps$rare$plots$annotations,
                        shared = htmps$shared$plots$annotations
                    )
                )
            ),
            list(
                genes = list(
                    rare = rare_shared_emt[[ct]]$rare_shared_emt_data$rare_emt_genes,
                    shared = rare_shared_emt[[ct]]$rare_shared_emt_data$shared_emt_genes
                )
            ),
            list(figure_widths = setNames(htmps[[1]]$figure_widths, c('annotations', 'heatmap'))),
            htmps[[1]][
                c(
                    'annotations_side',
                    'es_plot_type',
                    'es_title',
                    'es_upper_limit'
                )
                ]
        )

    },
    simplify = FALSE,
    USE.NAMES = TRUE
)

# Make the figure:

cairo_pdf(
    '../data_and_figures/single_cell_heatmaps_rare_shared.pdf',
    width = 10,
    height = 7
)

plot_grid(
    plot_grid(
        cowplot_sc_rare_shared(
            sc_cancer_heatmaps_z_score_restricted$hnsc,
            legend_space = 0,
            heights = c(1.5, 20, 4),
            es_y_axis_title_angle = 0,
            es_y_axis_title_xpos = 0.8,
            heatmap_title = 'Head and Neck'
        ),
        blank_plot(),
        cowplot_sc_rare_shared(
            sc_cancer_heatmaps_z_score_restricted$lung,
            legend_space = 0,
            heights = c(1.5, 20, 4),
            es_y_axis_title_angle = 0,
            es_y_axis_title_xpos = 0.2,
            heatmap_title = 'Lung'
        ),
        ncol = 3,
        nrow = 1,
        rel_widths = c(20, 1, 20)
    ),
    blank_plot(),
    plot_grid(
        cowplot_sc_rare_shared(
            sc_cancer_heatmaps_z_score_restricted$paad,
            legend_space = 0,
            heights = c(1.5, 20, 4),
            es_y_axis_title_angle = 0,
            es_y_axis_title_xpos = 0.8,
            heatmap_title = 'Pancreas'
        ),
        blank_plot(),
        sc_cancer_heatmaps_z_score_restricted$paad$plots$heatmap_legend,
        sc_cancer_heatmaps_z_score_restricted$paad$plots$genes_detected_legend,
        blank_plot(),
        ncol = 5,
        nrow = 1,
        rel_widths = c(20, 1, 6, 4, 10)
    ),
    blank_plot(),
    ncol = 1,
    nrow = 4,
    rel_heights = c(20, 1, 20, 1, 20, 1)
)

dev.off()





# Make cancer-caf heatmaps using average expression as expression summary:

sc_cancer_fibroblast_heatmaps <- sapply(
    names(single_cell_metadata),
    function(ct) {
        set.seed(sc_cancer_fibroblast_args[[ct]]$seed)
        do.call(
            sc_groups_heatmap,
            args = c(
                list(
                    sc_groups_list = sc_cancer_fibroblast[[ct]],
                    groups = c('cancer', 'fibroblast'),
                    default_figure_widths = list(
                        annotations = 2,
                        cancer = 6,
                        fibroblast = 1.2
                    ),
                    figure_spacing = 2.5,
                    annotations_nudge = 0.25,
                    h_legend_title = 'Expression level\n',
                    h_legend_width = 20,
                    h_legend_height = 10,
                    h_legend_direction = 'horizontal',
                    h_legend_title_position = 'right',
                    h_legend_just = 'left',
                    gd_legend_title = 'Genes detected\n',
                    gd_legend_width = 20,
                    gd_legend_height = 10,
                    gd_legend_direction = 'horizontal',
                    gd_legend_title_position = 'left',
                    gd_legend_just = 'right'
                )#,
                # sc_cancer_fibroblast_heatmaps_args[[ct]]
            )
        )
    },
    simplify = FALSE,
    USE.NAMES = TRUE
)





# To see the genes (in table form):

all_rare_shared <- unlist(
    lapply(
        rare_shared_emt,
        function(li) li$rare_shared_emt_data[c('rare_emt_genes', 'shared_emt_genes')]
    ),
    recursive = FALSE
)

knitr::kable(
    as.data.frame(
        lapply(
            all_rare_shared,
            function(x) {
                c(x, rep('', max(sapply(all_rare_shared, length)) - length(x)))
            }
        ),
        col.names = gsub('_emt_genes', '', names(all_rare_shared))
    )
)

# Now it would be good to remake the heatmaps using (some of) the rare and shared EMT genes
# as labels, with a different colour for each...  Or, I could label the rare EMT genes on
# one side and the shared ones on the other side.





# Try a scatterplot of scores:

rare_shared_emt_scores <- sapply(
    names(sc_cancer),
    function(ct) {

        dt <- data.table(
            gene_id = sc_cancer[[ct]]$genes_filtered
        )[
            ,
            c('ave_exp', 'emt_cor', 'skewness') := .(
                sc_cancer[[ct]]$data[
                    sc_cancer[[ct]]$cells_filtered,
                    # apply(.SD, 2, mean),
                    scale(apply(.SD, 2, mean)),
                    .SDcols = gene_id
                    ],
                rare_shared_emt[[ct]]$score_dist_data$score_cor_data[
                    gene_id,
                    # emt_cor
                    scale(emt_cor)
                    ],
                # rare_shared_emt[[ct]]$rare_shared_emt_data$skewnesses[gene_id]
                scale(rare_shared_emt[[ct]]$rare_shared_emt_data$skewnesses[gene_id])
            )
            ][
                ,
                emt_type := switch(
                    (gene_id %in% rare_shared_emt[[ct]]$rare_shared_emt_data$rare_emt_genes) + 1,
                    switch(
                        (gene_id %in% rare_shared_emt[[ct]]$rare_shared_emt_data$shared_emt_genes) + 1,
                        'unclassified',
                        'shared'
                    ),
                    'rare'
                ),
                by = gene_id
                ]

        # fig <- qplot(
        #     x = emt_cor + skewness,
        #     y = ave_exp - skewness,
        #     colour = emt_type,
        #     data = dt
        # ) +
        #     labs(
        #         x = 'Rare EMT score',
        #         y = 'Shared EMT score',
        #         title = ct
        #     ) +
        #     theme_test()

        # Take logs to reign in the more extreme values:

        fig <- qplot(
            x = log2(abs(emt_cor + skewness) + 1)*sign(emt_cor + skewness),
            y = log2(abs(ave_exp - skewness) + 1)*sign(ave_exp - skewness),
            colour = emt_type,
            data = dt
        ) +
            labs(
                x = 'Rare EMT score',
                y = 'Shared EMT score',
                title = ct
            ) +
            theme_test()

        list(
            data = dt,
            plot = fig
        )

    },
    simplify = FALSE,
    USE.NAMES = TRUE
)

# Save to PDF (just hnsc, lung and paad):

pdf('../data_and_figures/rare_emt_scores.pdf', width = 5, height = 5)

plot_grid(
    plotlist = c(
        lapply(
            rare_shared_emt_scores[c('hnsc', 'lung', 'paad')],
            function(li) {
                li$plot + theme(legend.position = 'none')
            }
        ),
        list(get_legend(rare_shared_emt_scores[[1]]$plot))
    ),
    nrow = 2,
    ncol = 2
)

dev.off()





# Check correlations of rare EMT genes with clinical features:

deconv_data <- readRDS('../data_and_figures/all_deconvs.rds')

clinical_data <- fread('../../TCGA_data/tcga_clinical_data.csv', key = 'id')

# Restrict to those cancer types which we took in further_analysis.R and for which we
# have single cell data:

sample_ids_list <- lapply(
    deconv_data[
        c(
            # 'brca_luminal_a',
            # 'brca_luminal_b',
            # 'brca_basal_like',
            # 'brca_her2_enriched',
            # 'brca_normal_like',
            # 'coad',
            'hnsc_mesenchymal_basal',
            'hnsc_classical',
            'hnsc_atypical',
            'luad_proximal_inflammatory',
            'luad_proximal_proliferative',
            'lusc_basal',
            'lusc_classical',
            'lusc_primitive',
            'lusc_secretory',
            'paad'
            # 'read'
        )
        ],
    `[[`,
    'sample_ids'
)

genes <- sapply(
    names(sample_ids_list),
    function(ct) {
        rare_shared_emt[[
            plyr::mapvalues(
                stringr::str_extract(ct, '^[a-z][a-z][a-z][a-z]'),
                c('coad', 'luad', 'lusc', 'read'),
                c('coadread', 'lung', 'lung', 'coadread'),
                warn_missing = FALSE
            )
            ]]$rare_shared_emt_data$rare_emt_genes
    },
    simplify = FALSE,
    USE.NAMES = TRUE
)

clin_cor <- volcano_clinical_list(
    expression_data,
    sample_ids_list,
    genes,
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
        quote(
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
    labels_signif_threshold = 'non-adjusted',
    legend_title = 'Prognostic feature',
    legend_labels = c(
        'number_of_lymphnodes_positive_by_he' = 'Lymph node metastasis',
        'followup_treatment_success' = 'Therapy resistance',
        'days_to_death' = 'Survival',
        'lymphovascular_invasion' = 'Lymphovascular invasion',
        'pathologic_t' = 'T stage',
        'pathologic_n' = 'N stage',
        'pathologic_m' = 'M stage',
        'neoplasm_histologic_grade' = 'Grade'
    ),
    legend_colours = setNames(
        RColorBrewer::brewer.pal(8, "Set1"),
        c(
            'number_of_lymphnodes_positive_by_he',
            'followup_treatment_success',
            'days_to_death',
            'lymphovascular_invasion',
            'pathologic_t',
            'pathologic_n',
            'pathologic_m',
            'neoplasm_histologic_grade'
        )
    )
)

clin_cor <- list(

    lymph_node_metastasis_positive = volcano_clinical(
        expression_data,
        sample_ids_list,
        genes,
        clinical_data,
        clin_var = 'number_of_lymphnodes_positive_by_he',
        wilcox_test_x_expr = variable > 0,
        wilcox_test_y_expr = variable == 0,
        labels_signif_threshold = 'none',
        plot_title = 'Correlation with lymph node metastasis'
    ),

    lymph_node_metastasis_multiple = volcano_clinical(
        expression_data,
        sample_ids_list,
        genes,
        clinical_data,
        clin_var = 'number_of_lymphnodes_positive_by_he',
        wilcox_test_x_expr = variable > 1,
        wilcox_test_y_expr = variable == 0,
        labels_signif_threshold = 'none',
        plot_title = 'Correlation with lymph node metastasis (0 vs. >1)'
    ),

    therapy_resistance = volcano_clinical(
        expression_data,
        sample_ids_list,
        genes,
        clinical_data,
        clin_var = 'followup_treatment_success',
        wilcox_test_x_expr = variable != 'complete remission/response',
        wilcox_test_y_expr = variable == 'complete remission/response',
        labels_signif_threshold = 'none',
        plot_title = 'Correlation with resistance to therapy'
    ),

    mortality = volcano_clinical(
        expression_data,
        sample_ids_list,
        genes,
        clinical_data,
        clin_var = 'days_to_death',
        wilcox_test_x_expr = variable < median(variable),
        wilcox_test_y_expr = variable > median(variable),
        labels_signif_threshold = 'none',
        plot_title = 'Correlation with mortality (median)'
    ),

    mortality_strict = volcano_clinical(
        expression_data,
        sample_ids_list,
        genes,
        clinical_data,
        clin_var = 'days_to_death',
        wilcox_test_x_expr = variable < quantile(variable, 0.4),
        wilcox_test_y_expr = variable > quantile(variable, 0.6),
        labels_signif_threshold = 'none',
        plot_title = 'Correlation with mortality (40th vs. 60th percentile)'
    ),

    lymphovascular_invasion = volcano_clinical(
        expression_data,
        sample_ids_list,
        genes,
        clinical_data,
        clin_var = 'lymphovascular_invasion',
        wilcox_test_x_expr = variable %in% c('yes', 'present'),
        wilcox_test_y_expr = variable %in% c('no', 'absent'),
        labels_signif_threshold = 'none',
        plot_title = 'Correlation with lymphovascular invasion'
    ),

    t_stage = volcano_clinical(
        expression_data,
        sample_ids_list,
        genes,
        clinical_data,
        clin_var = 'pathologic_t',
        wilcox_test_x_expr = startsWith(variable, 't3') |
            startsWith(variable, 't4'),
        wilcox_test_y_expr = startsWith(variable, 't0') |
            startsWith(variable, 't1') |
            startsWith(variable, 't2'),
        labels_signif_threshold = 'none',
        plot_title = 'Correlation with T stage',
        amatch_max_dist = 2
    ),

    n_stage = volcano_clinical(
        expression_data,
        sample_ids_list,
        genes,
        clinical_data,
        clin_var = 'pathologic_n',
        wilcox_test_x_expr = startsWith(variable, 'n2') |
            startsWith(variable, 'n3'),
        wilcox_test_y_expr = startsWith(variable, 'n0') |
            startsWith(variable, 'n1'),
        labels_signif_threshold = 'none',
        plot_title = 'Correlation with N stage',
        amatch_max_dist = 2
    ),

    m_stage = volcano_clinical(
        expression_data,
        sample_ids_list,
        genes,
        clinical_data,
        clin_var = 'pathologic_m',
        wilcox_test_x_expr = startsWith(variable, 'm1'),
        wilcox_test_y_expr = startsWith(variable, 'm0') |
            startsWith(variable, 'cm0'),
        labels_signif_threshold = 'none',
        plot_title = 'Correlation with M stage',
        amatch_max_dist = 2
    ),

    grade = volcano_clinical(
        expression_data,
        sample_ids_list,
        genes,
        clinical_data,
        clin_var = 'neoplasm_histologic_grade',
        wilcox_test_x_expr = variable %in% c('gb', 'g1', 'g2', 'low grade'),
        wilcox_test_y_expr = variable %in% c('g3', 'g4', 'high grade'),
        labels_signif_threshold = 'none',
        plot_title = 'Correlation with tumour grade'
    )

)

pdf('../data_and_figures/clinical_analyses_rare_emt.pdf', width = 35, height = 14)

plot_grid(
    plot_grid(
        plotlist = lapply(
            clin_cor,
            function(li) {li$plot + theme(legend.position = 'none')}
        ),
        nrow = 2,
        ncol = 5
    ),
    plot_grid(
        get_legend(clin_cor[[1]]$plot),
        get_legend(clin_cor[[1]]$plot),
        nrow = 2,
        ncol = 1
    ),
    nrow = 1,
    ncol = 2,
    rel_widths = c(20, 1)
)

dev.off()

# The fold changes might have a better scale if we centred the data.  But then, if we're
# going to do this we might as well centre relative to the expected contribution of these
# genes from the CAFs, i.e. take a linear regression against average CAF marker expression
# and take residuals.

expression_data_list <- sapply(
    names(sample_ids_list),
    function(ct) {

        dt <- cbind(
            expression_data[sample_ids_list[[ct]], genes[[ct]], with = FALSE],
            caf_score = rowSums(
                expression_data[
                    sample_ids_list[[ct]],
                    tail(
                        all_deconvs[[ct]]$genes_filtered[all_deconvs[[ct]]$ordering],
                        20
                    ),
                    with = FALSE
                    ]
            )
        )

        as.data.table(
            c(
                list(id = sample_ids_list[[ct]]),
                sapply(
                    genes[[ct]],
                    function(g) {
                        lm(
                            formula(paste0('`', g, '` ~ caf_score')),
                            data = dt
                        )$residuals
                    },
                    simplify = FALSE,
                    USE.NAMES = TRUE
                )
            )
        )

    },
    simplify = FALSE,
    USE.NAMES = TRUE
)

clin_cor <- list(

    lymph_node_metastasis_positive = volcano_clinical(
        expression_data_list,
        sample_ids_list,
        genes,
        clinical_data,
        clin_var = 'number_of_lymphnodes_positive_by_he',
        wilcox_test_x_expr = variable > 0,
        wilcox_test_y_expr = variable == 0,
        labels_signif_threshold = 'none',
        fold_change_fun = function(fc) {log2(abs(fc) + 1)*sign(fc)},
        xlab_TeX = '$\\log_2$(|fold change| + 1)$\\times$sign(fold change)',
        plot_title = 'Correlation with lymph node metastasis'
    ),

    lymph_node_metastasis_multiple = volcano_clinical(
        expression_data_list,
        sample_ids_list,
        genes,
        clinical_data,
        clin_var = 'number_of_lymphnodes_positive_by_he',
        wilcox_test_x_expr = variable > 1,
        wilcox_test_y_expr = variable == 0,
        labels_signif_threshold = 'none',
        fold_change_fun = function(fc) {log2(abs(fc) + 1)*sign(fc)},
        xlab_TeX = '$\\log_2$(|fold change| + 1)$\\times$sign(fold change)',
        plot_title = 'Correlation with lymph node metastasis (0 vs. >1)'
    ),

    therapy_resistance = volcano_clinical(
        expression_data_list,
        sample_ids_list,
        genes,
        clinical_data,
        clin_var = 'followup_treatment_success',
        wilcox_test_x_expr = variable != 'complete remission/response',
        wilcox_test_y_expr = variable == 'complete remission/response',
        labels_signif_threshold = 'none',
        fold_change_fun = function(fc) {log2(abs(fc) + 1)*sign(fc)},
        xlab_TeX = '$\\log_2$(|fold change| + 1)$\\times$sign(fold change)',
        plot_title = 'Correlation with resistance to therapy'
    ),

    mortality = volcano_clinical(
        expression_data_list,
        sample_ids_list,
        genes,
        clinical_data,
        clin_var = 'days_to_death',
        wilcox_test_x_expr = variable < median(variable),
        wilcox_test_y_expr = variable > median(variable),
        labels_signif_threshold = 'none',
        fold_change_fun = function(fc) {log2(abs(fc) + 1)*sign(fc)},
        xlab_TeX = '$\\log_2$(|fold change| + 1)$\\times$sign(fold change)',
        plot_title = 'Correlation with mortality (median)'
    ),

    mortality_strict = volcano_clinical(
        expression_data_list,
        sample_ids_list,
        genes,
        clinical_data,
        clin_var = 'days_to_death',
        wilcox_test_x_expr = variable < quantile(variable, 0.4),
        wilcox_test_y_expr = variable > quantile(variable, 0.6),
        labels_signif_threshold = 'none',
        fold_change_fun = function(fc) {log2(abs(fc) + 1)*sign(fc)},
        xlab_TeX = '$\\log_2$(|fold change| + 1)$\\times$sign(fold change)',
        plot_title = 'Correlation with mortality (40th vs. 60th percentile)'
    ),

    lymphovascular_invasion = volcano_clinical(
        expression_data_list,
        sample_ids_list,
        genes,
        clinical_data,
        clin_var = 'lymphovascular_invasion',
        wilcox_test_x_expr = variable %in% c('yes', 'present'),
        wilcox_test_y_expr = variable %in% c('no', 'absent'),
        labels_signif_threshold = 'none',
        fold_change_fun = function(fc) {log2(abs(fc) + 1)*sign(fc)},
        xlab_TeX = '$\\log_2$(|fold change| + 1)$\\times$sign(fold change)',
        plot_title = 'Correlation with lymphovascular invasion'
    ),

    t_stage = volcano_clinical(
        expression_data_list,
        sample_ids_list,
        genes,
        clinical_data,
        clin_var = 'pathologic_t',
        wilcox_test_x_expr = startsWith(variable, 't3') |
            startsWith(variable, 't4'),
        wilcox_test_y_expr = startsWith(variable, 't0') |
            startsWith(variable, 't1') |
            startsWith(variable, 't2'),
        labels_signif_threshold = 'none',
        fold_change_fun = function(fc) {log2(abs(fc) + 1)*sign(fc)},
        xlab_TeX = '$\\log_2$(|fold change| + 1)$\\times$sign(fold change)',
        plot_title = 'Correlation with T stage',
        amatch_max_dist = 2
    ),

    n_stage = volcano_clinical(
        expression_data_list,
        sample_ids_list,
        genes,
        clinical_data,
        clin_var = 'pathologic_n',
        wilcox_test_x_expr = startsWith(variable, 'n2') |
            startsWith(variable, 'n3'),
        wilcox_test_y_expr = startsWith(variable, 'n0') |
            startsWith(variable, 'n1'),
        labels_signif_threshold = 'none',
        fold_change_fun = function(fc) {log2(abs(fc) + 1)*sign(fc)},
        xlab_TeX = '$\\log_2$(|fold change| + 1)$\\times$sign(fold change)',
        plot_title = 'Correlation with N stage',
        amatch_max_dist = 2
    ),

    m_stage = volcano_clinical(
        expression_data_list,
        sample_ids_list,
        genes,
        clinical_data,
        clin_var = 'pathologic_m',
        wilcox_test_x_expr = startsWith(variable, 'm1'),
        wilcox_test_y_expr = startsWith(variable, 'm0') |
            startsWith(variable, 'cm0'),
        labels_signif_threshold = 'none',
        fold_change_fun = function(fc) {log2(abs(fc) + 1)*sign(fc)},
        xlab_TeX = '$\\log_2$(|fold change| + 1)$\\times$sign(fold change)',
        plot_title = 'Correlation with M stage',
        amatch_max_dist = 2
    ),

    grade = volcano_clinical(
        expression_data_list,
        sample_ids_list,
        genes,
        clinical_data,
        clin_var = 'neoplasm_histologic_grade',
        wilcox_test_x_expr = variable %in% c('gb', 'g1', 'g2', 'low grade'),
        wilcox_test_y_expr = variable %in% c('g3', 'g4', 'high grade'),
        labels_signif_threshold = 'none',
        fold_change_fun = function(fc) {log2(abs(fc) + 1)*sign(fc)},
        xlab_TeX = '$\\log_2$(|fold change| + 1)$\\times$sign(fold change)',
        plot_title = 'Correlation with tumour grade'
    )

)

plot_grid(
    plot_grid(
        plotlist = lapply(
            clin_cor,
            function(li) {li$plot + theme(legend.position = 'none')}
        ),
        nrow = 2,
        ncol = 5
    ),
    plot_grid(
        get_legend(clin_cor[[1]]$plot),
        get_legend(clin_cor[[1]]$plot),
        nrow = 2,
        ncol = 1
    ),
    nrow = 1,
    ncol = 2,
    rel_widths = c(20, 1)
)

# This doesn't seem to work well...  A lot of the significant cases or either lost or
# have negative fold change, meaning that lower expression of the rare EMT genes relative
# to the CAF frequency is associated with poor clinical outcome...  This suggests it's
# really the fibroblasts themselves that lead to poor prognosis, and any effect they have
# of driving EMT is not significant.  Or, it means that the linear regression method is
# just rubbish, which it probably is, or that my marker genes are a poor proxy for CAF
# frequency (or TME more generally).





# For single-cell-simulated-deconv comparison: make versions where we also filter out genes that don't fit the trend:

# We will want different thresholds for "anomalous" genes in each cancer type:
anomalous_threshold <- c(
	brca = 0.5,
	brca_lenient = 0.5,
	coadread = 0.6,
	coadread_lenient = 0.75,
	hnsc = 0.75,
	lihc = 0.75,
	lihc_lenient = 0.75,
	luad = 0.75,
	lusc = 0.5,
	lusc_lenient = 0.5,
	ov = 0.75,
	ov_lenient = 0.75,
	paad = 0.6
)

sc_deconv_comparison_filtered_twice <- sapply(
	names(sc_deconv_comparison_centred),
	function(ct) {

		cat(paste0(ct, '\n'))

		simulated_deconv_ct <- do.call(
			deconvolve_emt_caf_data,
			args = c(
				list(
					expression_data = simulated_bulk_data,
					meta_data = simulated_bulk_metadata,
					genes = sc_deconv_comparison_centred[[ct]]$gene_averages_cancer_caf[
						,
						.(pass = ave_exp[cell_type == 'cancer'] > 0.5 | ave_exp[cell_type == 'caf'] > 0.5),
						by = gene
					][pass == TRUE, as.character(gene)],
					# genes = sc_deconv_comparison_centred[[ct]]$gene_averages[ave_exp >= 0.5, as.character(symbol)], # symbol is a factor for some reason
					# genes = gene_averages[
						# ,
						# .(pass = ave_exp[cancer_caf == 'cancer'] > 0.5 | ave_exp[cancer_caf == 'caf'] > 0.5),
						# by = symbol
					# ][pass == TRUE, as.character(symbol)],
					cell_type_markers = cell_type_markers,
					ccle_data = ccle,
					genes_from_tcga_fun = NULL,
					genes_filter_fun = NULL
				),
				deconv_args_per_ct[[ct]][!(names(deconv_args_per_ct[[ct]]) %in% c('plot_title', 'genes_filter_fun'))]
			)
		)

		# Calculate new scores as in deconv_reorder() function.  Genes with positive score should be on the left side and those with
		# negative score should be on the right side.
		head_genes <- with(simulated_deconv_ct, head(genes_filtered[ordering], 20))
		tail_genes <- with(simulated_deconv_ct, tail(genes_filtered[ordering], 20))
		new_scores <- sapply(
			colnames(simulated_deconv_ct$cor_mat),
			function(g) {
				meanvec <- c(mean(simulated_deconv_ct$cor_mat[, g][head_genes]), mean(simulated_deconv_ct$cor_mat[, g][tail_genes]))
				ifelse(
					sign(meanvec[1]) == sign(meanvec[2]),
					ifelse(g %in% c('SNAI1', 'SNAI2', 'TWIST1', 'VIM','ZEB1', 'ZEB2'), return(0), return(NA)),
					return(max(meanvec)*c(1, -1)[which.max(meanvec)])
				)
			}
		)
		new_scores <- new_scores[!is.na(new_scores)]

		# For each gene, calculate the difference between its average across cancer cells and its average across CAFs.  For genes with
		# positive score, this difference should be positive, and for those with negative score it should be negative.  There will be
		# some exceptions - we need to decide if an exception is "extreme" enough.
		exp_diff_table <- sc_deconv_comparison_centred[[ct]]$plot_data[
			,
			.(exp_diff = mean(expression_level[cell_type == 'cancer']) - mean(expression_level[cell_type == 'caf'])),
			by = gene
		]

		# For each gene, consider the distribution of these differences for the genes having score with opposite sign.  If the gene is
		# beyond a certain quantile of this distribution ("beyond" w.r.t. its "own" distribution) then it is "anomalous".
		new_scores <- sapply(
			names(new_scores),
			function(g) {
				ifelse(
					new_scores[g] > 0,
					ifelse(
						exp_diff_table[
							,
							exp_diff[gene == g] < quantile(exp_diff[gene %in% names(new_scores)[new_scores < 0]], anomalous_threshold[ct])
						] & !(g %in% c('SNAI1', 'SNAI2', 'TWIST1', 'VIM', 'ZEB1', 'ZEB2')),
						return(NA),
						return(setNames(new_scores[g], NULL))
					),
					ifelse(
						new_scores[g] < 0,
						ifelse(
							exp_diff_table[
								,
								exp_diff[gene == g] > quantile(exp_diff[gene %in% names(new_scores)[new_scores > 0]], 1 - anomalous_threshold[ct])
							] & !(g %in% c('SNAI1', 'SNAI2', 'TWIST1', 'VIM', 'ZEB1', 'ZEB2')),
							return(NA),
							return(setNames(new_scores[g], NULL))
						),
						return(0)
					)
				)
				# exp_diff_table[, exp_diff[names(new_scores)[sign(new_scores) != sign(new_scores[g])]]
			}
		)

		anomalous_genes <- names(new_scores)[is.na(new_scores)]
		new_scores <- new_scores[!is.na(new_scores)]

		simulated_deconv_ct$anomalous_genes <- anomalous_genes
		simulated_deconv_ct$new_scores <- new_scores
		simulated_deconv_ct$cor_mat <- simulated_deconv_ct$cor_mat[names(new_scores), names(new_scores)]
		simulated_deconv_ct$genes_filtered <- names(new_scores)
		simulated_deconv_ct$ordering <- order(-new_scores)
		simulated_deconv_ct$cor_with_purity <- sapply(
			simulated_deconv_ct$cor_with_purity,
			function(x) x[names(new_scores)],
			simplify = FALSE,
			USE.NAMES = TRUE
		)
		simulated_deconv_ct$ccle_comp_diff <- simulated_deconv_ct$ccle_comp_diff[names(new_scores)]

		simulated_deconv_plots_ct <- do.call(
			deconvolve_emt_caf_plots,
			args = c(
				list(
					data = simulated_deconv_ct,
					heatmap_legend_title = 'Correlation',
					heatmap_colours = rev(colorRampPalette(brewer.pal(11, "RdBu"))(50)),
					heatmap_colour_limits = c(-1, 1),
					heatmap_legend_breaks = c(-1, -0.5, 0, 0.5, 1),
					heatmap_annotations_nudge = 0.3,
					purity_colours = rev(colorRampPalette(brewer.pal(11, "PuOr"))(50)),
					purity_colour_limits = c(-1, 1),
					purity_legend_breaks = c(-1, 0, 1),
					purity_legend_title = 'Correlation with purity\n',
					purity_legend_direction = 'horizontal',
					purity_axis_title = NULL,
					ccle_colours = rev(colorRampPalette(brewer.pal(11, 'PiYG')[c(1:3, 6, 9:11)])(50)),
					ccle_legend_breaks = c(-1, 0, 1),
					ccle_legend_title = 'Tumours vs. cell lines\n',
					ccle_legend_direction = 'horizontal',
					ccle_axis_title = NULL,
					bar_legend_justification = 'left'
				),
				deconv_plot_args_per_ct[[ct]]
			)
		)

		plot_data <- sc_deconv_comparison_centred[[ct]]$plot_data[
			gene %in% sc_deconv_comparison_centred[[ct]]$gene_averages_cancer_caf[
				,
				.(pass = ave_exp[cell_type == 'cancer'] > 0.5 | ave_exp[cell_type == 'caf'] > 0.5),
				by = gene
			][pass == TRUE, as.character(gene)]
			# gene %in% sc_deconv_comparison_centred[[ct]]$gene_averages[ave_exp >= 0.5, symbol]
		]

		sc_heatmaps <- sapply(
			c('cancer', 'caf'),
			function(x) {
				ggplot(
					plot_data[cell_type == x],
					aes(
						x = factor(gene, levels = with(simulated_deconv_ct, genes_filtered[ordering])),
						y = factor(id, levels = sc_deconv_comparison[[ct]][[x]]$ordered_cell_ids),
						fill = expression_level
					)
				) +
					geom_raster() +
					scale_x_discrete(expand = c(0, 0)) +
					scale_y_discrete(expand = c(0, 0)) +
					scale_fill_gradientn(
						colours = rev(colorRampPalette(brewer.pal(11, "RdBu"))(50)),
						limits = c(-6, 6),
						oob = scales::squish,
						breaks = c(-6, -3, 0, 3, 6),
						labels = c('-6' = '\u2264 -6', '-3' = '-3', '0' = '0', '3' = '3', '6' = '\u2265 6')
					) +
					theme(
						axis.text = element_blank(),
						axis.title.x = element_blank(),
						axis.ticks = element_blank(),
						axis.ticks.length = unit(0, 'pt'),
						plot.margin = unit(c(5.5, 5.5, 1.25, 5.5), 'pt')
					) +
					labs(y = mapvalues(x, c('cancer', 'caf'), c('Cancer', 'CAF'), warn_missing = FALSE), fill = 'Relative\nexpression\nlevel')
			},
			simplify = FALSE,
			USE.NAMES = TRUE
		)

		list(deconv_data = simulated_deconv_ct, deconv_figures = simulated_deconv_plots_ct$plots, sc_heatmaps = sc_heatmaps, plot_data = plot_data)

	},
	simplify = FALSE,
	USE.NAMES = TRUE
)

cairo_pdf('../data_and_figures/sc_sim_deconv_comp_filtered_twice.pdf', width = 6, height = 9, onefile = TRUE)

for(ct in names(sc_deconv_comparison_filtered_twice)[!grepl('lenient', names(sc_deconv_comparison_filtered_twice))]) {

	plot_grid(
		plot_grid(
			plotlist = c(
				list(blank_plot() + theme(plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), 'pt')) + labs(title = deconv_plot_args_per_ct[[ct]]$plot_title)),
				lapply(
					c(
						sc_deconv_comparison_filtered_twice[[ct]]$deconv_figures[c('purity_bar', 'ccle_bar', 'heatmap')],
						sc_deconv_comparison_filtered_twice[[ct]]$sc_heatmaps
					),
					function(x) x + theme(legend.position = 'none', plot.title = element_blank())
				)
			),
			nrow = 6,
			ncol = 1,
			align = 'v',
			rel_heights = c(1, 1, 1, 15, 5, 5)
		),
		plot_grid(
			plotlist = c(
				list(
					blank_plot(),
					get_legend(
						simulated_deconv_plots[[1]]$plots$purity_bar +
							theme(legend.justification = c(0, 1), legend.direction = 'vertical', legend.key.width = NULL, legend.key.height = NULL) +
							labs(fill = 'Correlation\nwith purity')
					),
					get_legend(
						simulated_deconv_plots[[1]]$plots$ccle_bar +
							theme(legend.justification = c(0, 1), legend.direction = 'vertical', legend.key.width = NULL, legend.key.height = NULL) +
							labs(fill = 'Tumours vs.\ncell lines')
					),
					get_legend(
						simulated_deconv_plots[[1]]$plots$heatmap + theme(legend.justification = c(0, 1))
					),
					get_legend(sc_deconv_comparison_filtered_twice[[1]]$sc_heatmaps[[1]] + theme(legend.justification = c(0, 1)))
				)
			),
			nrow = 5,
			ncol = 1,
			rel_heights = c(1, 5, 5, 7, 10)
		),
		nrow = 1,
		ncol = 2,
		rel_widths = c(5, 1)
	) %>% print

}

dev.off()

cairo_pdf('../data_and_figures/sc_sim_deconv_comp_lenient_filtered_twice.pdf', width = 6, height = 10, onefile = TRUE)

for(ct in names(sc_deconv_comparison_filtered_twice)[grepl('lenient', names(sc_deconv_comparison_filtered_twice))]) {

	plot_grid(
		plot_grid(
			plotlist = c(
				list(blank_plot() + theme(plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), 'pt')) + labs(title = deconv_plot_args_per_ct[[ct]]$plot_title)),
				lapply(
					c(
						sc_deconv_comparison_filtered_twice[[ct]]$deconv_figures[c('purity_bar', 'ccle_bar', 'heatmap')],
						sc_deconv_comparison_filtered_twice[[ct]]$sc_heatmaps
					),
					function(x) x + theme(legend.position = 'none', plot.title = element_blank())
				)
			),
			nrow = 6,
			ncol = 1,
			align = 'v',
			rel_heights = c(1, 1, 1, 15, 5, 5)
		),
		plot_grid(
			plotlist = c(
				list(
					blank_plot(),
					get_legend(
						simulated_deconv_plots[[1]]$plots$purity_bar +
							theme(legend.justification = c(0, 1), legend.direction = 'vertical', legend.key.width = NULL, legend.key.height = NULL) +
							labs(fill = 'Correlation\nwith purity')
					),
					get_legend(
						simulated_deconv_plots[[1]]$plots$ccle_bar +
							theme(legend.justification = c(0, 1), legend.direction = 'vertical', legend.key.width = NULL, legend.key.height = NULL) +
							labs(fill = 'Tumours vs.\ncell lines')
					),
					get_legend(
						simulated_deconv_plots[[1]]$plots$heatmap + theme(legend.justification = c(0, 1))
					),
					get_legend(sc_deconv_comparison_filtered_twice[[1]]$sc_heatmaps[[1]] + theme(legend.justification = c(0, 1)))
				)
			),
			nrow = 5,
			ncol = 1,
			rel_heights = c(1, 5, 5, 7, 10)
		),
		nrow = 1,
		ncol = 2,
		rel_widths = c(5, 1)
	) %>% print

}

dev.off()





# Similarly for single-cell-TCGA-deconv comparison, make versions where we also filter out genes that don't fit the trend:

anomalous_threshold_tcga <- c(
	brca_luminal_a = 0.4,
	brca_luminal_b = 0.3,
	brca_basal_like = 0.3,
	brca_her2_enriched = 0.25,
	coad = 0.5,
	hnsc_mesenchymal_basal = 0.3,
	hnsc_classical = 0.3,
	hnsc_atypical = 0.5,
	lihc = 0.4,
	luad_proximal_inflammatory = 0.4,
	luad_proximal_proliferative = 0.5,
	lusc_basal = 0.3,
	lusc_classical = 0.6,
	lusc_primitive = 0.4,
	lusc_secretory = 0.5,
	ov_differentiated = 0.3,
	ov_immunoreactive = 0.5,
	ov_proliferative = 0.3,
	paad = 0.6,
	read = 0.3
)

sc_tcga_deconv_comparison_filtered_twice <- sapply(
	names(sc_tcga_deconv_comparison_centred),
	function(ct) {

		cat(paste0(ct, '\n'))

		deconv_ct <- do.call(
			deconvolve_emt_caf_data,
			args = c(
				list(
					expression_data = expression_data,
					meta_data = meta_data,
					genes = sc_tcga_deconv_comparison_centred[[ct]]$gene_averages_cancer_caf[
						,
						.(pass = ave_exp[cell_type == 'cancer'] > 0.5 | ave_exp[cell_type == 'caf'] > 0.5),
						by = gene
					][pass == TRUE, as.character(gene)],
					# genes = sc_tcga_deconv_comparison_centred[[ct]]$gene_averages[ave_exp >= 0.5, as.character(symbol)], # symbol is a factor for some reason
					cell_type_markers = cell_type_markers,
					ccle_data = ccle,
					subtypes_data = subtypes_data,
					genes_from_tcga_fun = NULL,
					genes_filter_fun = NULL
				),
				deconv_args_per_ct[[ct]][!(names(deconv_args_per_ct[[ct]]) %in% c('plot_title', 'genes_filter_fun'))]
			)
		)

		# Calculate new scores as in deconv_reorder() function.  Genes with positive score should be on the left side and those with
		# negative score should be on the right side.
		head_genes <- with(deconv_ct, head(genes_filtered[ordering], 20))
		tail_genes <- with(deconv_ct, tail(genes_filtered[ordering], 20))
		new_scores <- sapply(
			colnames(deconv_ct$cor_mat),
			function(g) {
				meanvec <- c(mean(deconv_ct$cor_mat[, g][head_genes]), mean(deconv_ct$cor_mat[, g][tail_genes]))
				ifelse(
					sign(meanvec[1]) == sign(meanvec[2]),
					ifelse(g %in% c('SNAI1', 'SNAI2', 'TWIST1', 'VIM','ZEB1', 'ZEB2'), return(0), return(NA)),
					return(max(meanvec)*c(1, -1)[which.max(meanvec)])
				)
			}
		)
		new_scores <- new_scores[!is.na(new_scores)]

		# For each gene, calculate the difference between its average across cancer cells and its average across CAFs.  For genes with
		# positive score, this difference should be positive, and for those with negative score it should be negative.  There will be
		# some exceptions - we need to decide if an exception is "extreme" enough.
		exp_diff_table <- sc_tcga_deconv_comparison_centred[[ct]]$plot_data[
			,
			.(exp_diff = mean(expression_level[cell_type == 'cancer']) - mean(expression_level[cell_type == 'caf'])),
			by = gene
		]

		# For each gene, consider the distribution of these differences for the genes having score with opposite sign.  If the gene is
		# beyond a certain quantile of this distribution ("beyond" w.r.t. its "own" distribution) then it is "anomalous".
		new_scores <- sapply(
			names(new_scores),
			function(g) {
				ifelse(
					new_scores[g] > 0,
					ifelse(
						exp_diff_table[
							,
							exp_diff[gene == g] < quantile(exp_diff[gene %in% names(new_scores)[new_scores < 0]], anomalous_threshold_tcga[ct])
						] & !(g %in% c('SNAI1', 'SNAI2', 'TWIST1', 'VIM', 'ZEB1', 'ZEB2')),
						return(NA),
						return(setNames(new_scores[g], NULL))
					),
					ifelse(
						new_scores[g] < 0,
						ifelse(
							exp_diff_table[
								,
								exp_diff[gene == g] > quantile(exp_diff[gene %in% names(new_scores)[new_scores > 0]], 1 - anomalous_threshold_tcga[ct])
							] & !(g %in% c('SNAI1', 'SNAI2', 'TWIST1', 'VIM', 'ZEB1', 'ZEB2')),
							return(NA),
							return(setNames(new_scores[g], NULL))
						),
						return(0)
					)
				)
				# exp_diff_table[, exp_diff[names(new_scores)[sign(new_scores) != sign(new_scores[g])]]
			}
		)

		anomalous_genes <- names(new_scores)[is.na(new_scores)]
		new_scores <- new_scores[!is.na(new_scores)]

		deconv_ct$anomalous_genes <- anomalous_genes
		deconv_ct$new_scores <- new_scores
		deconv_ct$cor_mat <- deconv_ct$cor_mat[names(new_scores), names(new_scores)]
		deconv_ct$genes_filtered <- names(new_scores)
		deconv_ct$ordering <- order(-new_scores)
		deconv_ct$cor_with_purity <- sapply(
			deconv_ct$cor_with_purity,
			function(x) x[names(new_scores)],
			simplify = FALSE,
			USE.NAMES = TRUE
		)
		deconv_ct$ccle_comp_diff <- deconv_ct$ccle_comp_diff[names(new_scores)]

		deconv_plots_ct <- do.call(
			deconvolve_emt_caf_plots,
			args = c(
				list(
					data = deconv_ct,
					heatmap_legend_title = 'Correlation',
					heatmap_colours = rev(colorRampPalette(brewer.pal(11, "RdBu"))(50)),
					heatmap_colour_limits = c(-0.6, 0.6),
					heatmap_legend_breaks = c(-0.6, -0.3, 0, 0.3, 0.6),
					heatmap_annotations_nudge = 0.3,
					purity_colours = rev(colorRampPalette(brewer.pal(11, "PuOr"))(50)),
					purity_colour_limits = c(-0.3, 0.3),
					purity_legend_breaks = c(-0.3, 0, 0.3),
					purity_legend_title = 'Correlation with purity\n',
					purity_legend_direction = 'horizontal',
					purity_axis_title = NULL,
					ccle_colours = rev(colorRampPalette(brewer.pal(11, 'PiYG')[c(1:3, 6, 9:11)])(50)),
					ccle_legend_breaks = c(-1, 0, 1),
					ccle_legend_title = 'Tumours vs. cell lines\n',
					ccle_legend_direction = 'horizontal',
					ccle_axis_title = NULL,
					bar_legend_justification = 'left'
				),
				deconv_args_per_ct[[ct]]['plot_title']
			)
		)

		plot_data <- sc_tcga_deconv_comparison_centred[[ct]]$plot_data[
			gene %in% sc_tcga_deconv_comparison_centred[[ct]]$gene_averages_cancer_caf[
				,
				.(pass = ave_exp[cell_type == 'cancer'] > 0.5 | ave_exp[cell_type == 'caf'] > 0.5),
				by = gene
			][pass == TRUE, as.character(gene)]
			# gene %in% sc_tcga_deconv_comparison_centred[[ct]]$gene_averages[ave_exp >= 0.5, symbol]
		]

		sc_heatmaps <- sapply(
			c('cancer', 'caf'),
			function(x) {
				ggplot(
					plot_data[cell_type == x],
					aes(
						x = factor(gene, levels = with(deconv_ct, genes_filtered[ordering])),
						y = factor(id, levels = sc_tcga_deconv_comparison[[ct]][[x]]$ordered_cell_ids),
						fill = expression_level
					)
				) +
					geom_raster() +
					scale_x_discrete(expand = c(0, 0)) +
					scale_y_discrete(expand = c(0, 0)) +
					scale_fill_gradientn(
						colours = rev(colorRampPalette(brewer.pal(11, "RdBu"))(50)),
						limits = c(-6, 6),
						oob = scales::squish,
						breaks = c(-6, -3, 0, 3, 6),
						labels = c('-6' = '\u2264 -6', '-3' = '-3', '0' = '0', '3' = '3', '6' = '\u2265 6')
					) +
					theme(
						axis.text = element_blank(),
						axis.title.x = element_blank(),
						axis.ticks = element_blank(),
						axis.ticks.length = unit(0, 'pt'),
						plot.margin = unit(c(5.5, 5.5, 1.25, 5.5), 'pt')
					) +
					labs(y = mapvalues(x, c('cancer', 'caf'), c('Cancer', 'CAF'), warn_missing = FALSE), fill = 'Relative\nexpression\nlevel')
			},
			simplify = FALSE,
			USE.NAMES = TRUE
		)

		list(deconv_data = deconv_ct, deconv_figures = deconv_plots_ct$plots, sc_heatmaps = sc_heatmaps, plot_data = plot_data)

	},
	simplify = FALSE,
	USE.NAMES = TRUE
)

cairo_pdf('../data_and_figures/sc_deconv_comp_filtered_twice.pdf', width = 6, height = 9, onefile = TRUE)

for(ct in names(sc_tcga_deconv_comparison_filtered_twice)) {

	plot_grid(
		plot_grid(
			plotlist = c(
				list(blank_plot() + theme(plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), 'pt')) + labs(title = deconv_args_per_ct[[ct]]$plot_title)),
				lapply(
					c(
						sc_tcga_deconv_comparison_filtered_twice[[ct]]$deconv_figures[c('purity_bar', 'ccle_bar', 'heatmap')],
						sc_tcga_deconv_comparison_filtered_twice[[ct]]$sc_heatmaps
					),
					function(x) x + theme(legend.position = 'none', plot.title = element_blank())
				)
			),
			nrow = 6,
			ncol = 1,
			align = 'v',
			rel_heights = c(1, 1, 1, 15, 5, 5)
		),
		plot_grid(
			plotlist = c(
				list(
					blank_plot(),
					get_legend(
						deconv_plots[[1]]$plots$purity_bar +
							theme(legend.justification = c(0, 1), legend.direction = 'vertical', legend.key.width = NULL, legend.key.height = NULL) +
							labs(fill = 'Correlation\nwith purity')
					),
					get_legend(
						deconv_plots[[1]]$plots$ccle_bar +
							theme(legend.justification = c(0, 1), legend.direction = 'vertical', legend.key.width = NULL, legend.key.height = NULL) +
							labs(fill = 'Tumours vs.\ncell lines')
					),
					get_legend(
						deconv_plots[[1]]$plots$heatmap + theme(legend.justification = c(0, 1))
					),
					get_legend(sc_tcga_deconv_comparison_filtered_twice[[1]]$sc_heatmaps[[1]] + theme(legend.justification = c(0, 1)))
				)
			),
			nrow = 5,
			ncol = 1,
			rel_heights = c(1, 5, 5, 7, 10)
		),
		nrow = 1,
		ncol = 2,
		rel_widths = c(5, 1)
	) %>% print

}

dev.off()





sc_deconv_comparison <- sapply(
	names(simulated_deconvs),
	function(ct) {
		cat(paste0(ct, '\n'))
        sc_data <- eval(sc_metadata[[ct]]$read_quote)
		sapply(
			c('cancer', 'caf'),
			function(x) {

				plot_data <- copy(sc_data[cell_type == x])[, complexity := apply(.SD, 1, function(x) sum(x > 0)), .SDcols = -c('id', 'patient', 'cell_type')][
					,
					.SD[sample(1:.N, sample_sizes[[ct]], prob = complexity)]
				][, complexity := NULL][, c('id', simulated_deconvs[[ct]]$genes_filtered), with = FALSE]

				plot_data <- melt(plot_data, id.vars = 'id', variable.name = 'gene', value.name = 'expression_level')

				ordered_cell_ids <- plot_data[
					gene %in% with(
						simulated_deconvs[[ct]],
						do.call(c('head', 'tail')[which(c('cancer', 'caf') == x)], list(genes_filtered[ordering], 20))
					),
					.(top_20_mean = mean(expression_level)),
					by = id
				][order(top_20_mean), id]

				htmp <- ggplot(
					plot_data,
					aes(
						x = factor(gene, levels = with(simulated_deconvs[[ct]], genes_filtered[ordering])),
						y = factor(id, levels = ordered_cell_ids),
						fill = expression_level
					)
				) +
					geom_raster() +
					scale_x_discrete(expand = c(0, 0)) +
					scale_y_discrete(expand = c(0, 0)) +
					scale_fill_gradientn(
						# colours = c('turquoise4', 'turquoise', 'azure', 'gold2', 'gold4')
						colours = rev(colorRampPalette(brewer.pal(11, "Spectral"))(50)),
						limits = c(0, 12),
						oob = scales::squish,
						breaks = c(0, 3, 6, 9, 12),
						labels = c('0' = '0', '3' = '3', '6' = '6', '9' = '9', '12' = '\u2265 12')
					) +
					theme(
						axis.text = element_blank(),
						axis.title.x = element_blank(),
						axis.ticks = element_blank(),
						axis.ticks.length = unit(0, 'pt'),
						plot.margin = unit(c(5.5, 5.5, 1.25, 1.1), 'pt')
						# plot.margin = switch((x == 'cancer') + 1, unit(c(1.25, 5.5, 5.5, 5.5), 'pt'), unit(c(5.5, 5.5, 1.25, 5.5), 'pt'))
					) +
					labs(
						y = mapvalues(x, c('cancer', 'caf'), c('Cancer', 'CAF'), warn_missing = FALSE),
						fill = 'Expression\nlevel'
					)

				ave_exp_bar <- ggplot(
					plot_data[, .(ave_exp = mean(expression_level)), keyby = gene][
						with(simulated_deconvs[[ct]], genes_filtered[ordering]),
						ave_exp := runmean(ave_exp, 30)
					],
					aes(x = factor(gene, levels = with(simulated_deconvs[[ct]], genes_filtered[ordering])), y = 'a', fill = ave_exp)
				) +
					geom_raster() +
					scale_x_discrete(expand = c(0, 0)) +
					scale_y_discrete(expand = c(0, 0)) +
					scale_fill_gradientn(
						# colours = c('turquoise4', 'turquoise', 'azure', 'gold2', 'gold4')
						colours = colorRampPalette(brewer.pal(9, "RdPu"))(50),
						limits = c(0, 4),
						oob = scales::squish,
						breaks = c(0, 2, 4),
						labels = c('0' = '0', '2' = '2', '4' = '\u2265 4')
					) +
					theme(
						axis.text = element_blank(),
						axis.title = element_blank(),
						axis.ticks = element_blank(),
						axis.ticks.length = unit(0, 'pt'),
						plot.margin = switch(
							(x == 'cancer') + 1,
							unit(c(1.25, 5.5, 5.5, 5.5), 'pt'),
							unit(c(1.25, 5.5, 1.25, 5.5), 'pt')
						)
					) +
					labs(
						y = mapvalues(x, c('cancer', 'caf'), c('Cancer', 'CAF'), warn_missing = FALSE),
						fill = 'Average\nexpression\nlevel'
					)

				list(heatmap = htmp, ave_exp_bar = ave_exp_bar, plot_data = plot_data, ordered_cell_ids = ordered_cell_ids)

			},
			simplify = FALSE,
			USE.NAMES = TRUE
		)
	},
	simplify = FALSE,
	USE.NAMES = TRUE
)

cairo_pdf('../data_and_figures/sc_sim_deconv_comp.pdf', width = 6, height = 10, onefile = TRUE)

for(ct in names(sc_deconv_comparison)[!grepl('lenient', names(sc_deconv_comparison))]) {

	sim_deconv_figures <- sapply(
		simulated_deconv_plots[[ct]]$plots[c('purity_bar', 'ccle_bar', 'heatmap')],
		function(x) x + theme(legend.position = 'none', plot.title = element_blank()),
		simplify = FALSE,
		USE.NAMES = TRUE
	)
	sim_deconv_figures$heatmap <- sim_deconv_figures$heatmap + theme(plot.margin = unit(c(0, 5.5, 1.25, 5.5), 'pt'))

	plot_grid(
		plot_grid(
			plotlist = c(
				list(blank_plot() + theme(plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), 'pt')) + labs(title = deconv_plot_args_per_ct[[ct]]$plot_title)),
				sim_deconv_figures,
				# lapply(
					# simulated_deconv_plots[[ct]]$plots[c('purity_bar', 'ccle_bar', 'heatmap')],
					# function(x) x + theme(legend.position = 'none', plot.title = element_blank())
				# ),
				unlist(
					lapply(
						sc_deconv_comparison[[ct]],
						function(x) list(x$heatmap + theme(legend.position = 'none'), x$ave_exp_bar + theme(legend.position = 'none'))
					),
					recursive = FALSE
				)
			),
			nrow = 8,
			ncol = 1,
			align = 'v',
			rel_heights = c(1, 1, 1, 15, 5, 1, 5, 1)
		),
		plot_grid(
			plotlist = c(
				list(
					blank_plot(),
					get_legend(
						simulated_deconv_plots[[1]]$plots$purity_bar +
							theme(legend.justification = c(0, 1), legend.direction = 'vertical', legend.key.width = unit(10, 'pt'), legend.key.height = unit(10, 'pt')) +
							labs(fill = 'Correlation\nwith purity')
					),
					get_legend(
						simulated_deconv_plots[[1]]$plots$ccle_bar +
							theme(legend.justification = c(0, 1), legend.direction = 'vertical', legend.key.width = unit(10, 'pt'), legend.key.height = unit(10, 'pt')) +
							labs(fill = 'Tumours vs.\ncell lines')
					),
					get_legend(
						simulated_deconv_plots[[1]]$plots$heatmap + theme(legend.justification = c(0, 1))
					)
				),
				list(
					get_legend(sc_deconv_comparison[[1]][[1]]$heatmap + theme(legend.justification = c(0, 1))),
					get_legend(sc_deconv_comparison[[1]][[1]]$ave_exp_bar + theme(legend.justification = c(0, 1)))
				)
			),
			nrow = 6,
			ncol = 1,
			rel_heights = c(1, 5, 5, 7, 6, 6)
		),
		nrow = 1,
		ncol = 2,
		rel_widths = c(5, 1)
	) %>% print

}

dev.off()

cairo_pdf('../data_and_figures/sc_sim_deconv_comp_lenient.pdf', width = 6, height = 10, onefile = TRUE)

for(ct in names(sc_deconv_comparison)[grepl('lenient', names(sc_deconv_comparison))]) {

	sim_deconv_figures <- sapply(
		simulated_deconv_plots[[ct]]$plots[c('purity_bar', 'ccle_bar', 'heatmap')],
		function(x) x + theme(legend.position = 'none', plot.title = element_blank()),
		simplify = FALSE,
		USE.NAMES = TRUE
	)
	sim_deconv_figures$heatmap <- sim_deconv_figures$heatmap + theme(plot.margin = unit(c(0, 5.5, 1.25, 5.5), 'pt'))

	plot_grid(
		plot_grid(
			plotlist = c(
				list(blank_plot() + theme(plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), 'pt')) + labs(title = deconv_plot_args_per_ct[[ct]]$plot_title)),
				sim_deconv_figures,
				# lapply(
					# simulated_deconv_plots[[ct]]$plots[c('purity_bar', 'ccle_bar', 'heatmap')],
					# function(x) x + theme(legend.position = 'none', plot.title = element_blank())
				# ),
				unlist(
					lapply(
						sc_deconv_comparison[[ct]],
						function(x) list(x$heatmap + theme(legend.position = 'none'), x$ave_exp_bar + theme(legend.position = 'none'))
					),
					recursive = FALSE
				)
			),
			nrow = 8,
			ncol = 1,
			align = 'v',
			rel_heights = c(1, 1, 1, 15, 5, 1, 5, 1)
		),
		plot_grid(
			plotlist = c(
				list(
					blank_plot(),
					get_legend(
						simulated_deconv_plots[[1]]$plots$purity_bar +
							theme(legend.justification = c(0, 1), legend.direction = 'vertical', legend.key.width = unit(10, 'pt'), legend.key.height = unit(10, 'pt')) +
							labs(fill = 'Correlation\nwith purity')
					),
					get_legend(
						simulated_deconv_plots[[1]]$plots$ccle_bar +
							theme(legend.justification = c(0, 1), legend.direction = 'vertical', legend.key.width = unit(10, 'pt'), legend.key.height = unit(10, 'pt')) +
							labs(fill = 'Tumours vs.\ncell lines')
					),
					get_legend(
						simulated_deconv_plots[[1]]$plots$heatmap + theme(legend.justification = c(0, 1))
					)
				),
				list(
					get_legend(sc_deconv_comparison[[1]][[1]]$heatmap + theme(legend.justification = c(0, 1))),
					get_legend(sc_deconv_comparison[[1]][[1]]$ave_exp_bar + theme(legend.justification = c(0, 1)))
				)
			),
			nrow = 6,
			ncol = 1,
			rel_heights = c(1, 5, 5, 7, 6, 6)
		),
		nrow = 1,
		ncol = 2,
		rel_widths = c(5, 1)
	) %>% print

}

dev.off()

# Make centred versions:

sc_deconv_comparison_centred <- sapply(
	names(sc_deconv_comparison),
	function(ct) {

		cat(paste0(ct, '\n'))

		plot_data <- rbind(
			sc_deconv_comparison[[ct]]$cancer$plot_data[, cell_type := 'cancer'],
			sc_deconv_comparison[[ct]]$caf$plot_data[, cell_type := 'caf']
		)

		# Calculate averages in cancer cells and CAFs separately, so we can check them afterwards:
		gene_averages_cancer_caf <- plot_data[
			,
			.(ave_exp = mean(expression_level)),
			by = .(gene, cell_type)
		]

		# To centre genes w.r.t. the average of the averages of cancer and CAF:
		gene_averages <- plot_data[
			,
			.(ave_exp = mean(c(mean(expression_level[cell_type == 'cancer']), mean(expression_level[cell_type == 'caf'])))),
			by = .(symbol = gene)
		]

		plot_data[, expression_level := expression_level - gene_averages[symbol == unique(gene), ave_exp], by = gene]

		# To centre the cells as well:
		plot_data[, expression_level := expression_level - mean(expression_level), by = id]

		htmps <- sapply(
			c('cancer', 'caf'),
			function(x) {
				ggplot(
					plot_data[cell_type == x],
					aes(
						x = factor(gene, levels = with(simulated_deconvs[[ct]], genes_filtered[ordering])),
						y = factor(id, levels = sc_deconv_comparison[[ct]][[x]]$ordered_cell_ids),
						fill = expression_level
					)
				) +
					geom_raster() +
					scale_x_discrete(expand = c(0, 0)) +
					scale_y_discrete(expand = c(0, 0)) +
					scale_fill_gradientn(
						colours = c(
							colorspace::sequential_hcl(25, h = 190, c = 70, l = c(70, 99), power = 0.8),
							colorspace::sequential_hcl(25, h = 60, c = 100, l = c(80, 99), power = 0.8, rev = TRUE)
						),
						# colours = rev(colorRampPalette(brewer.pal(11, "RdBu"))(50)),
						limits = c(-6, 6),
						oob = scales::squish,
						breaks = c(-6, -3, 0, 3, 6),
						labels = c('-6' = '\u2264 -6', '-3' = '-3', '0' = '0', '3' = '3', '6' = '\u2265 6')
					) +
					theme(
						axis.text = element_blank(),
						axis.title.x = element_blank(),
						axis.ticks = element_blank(),
						axis.ticks.length = unit(0, 'pt'),
						plot.margin = unit(c(5.5, 5.5, 1.25, 5.5), 'pt')
						# plot.margin = switch((x == 'cancer') + 1, unit(c(1.25, 5.5, 5.5, 5.5), 'pt'), unit(c(5.5, 5.5, 1.25, 5.5), 'pt'))
					) +
					labs(
						y = mapvalues(x, c('cancer', 'caf'), c('Cancer', 'CAF'), warn_missing = FALSE),
						fill = 'Relative\nexpression\nlevel'
					)
			},
			simplify = FALSE,
			USE.NAMES = TRUE
		)

		list(heatmaps = htmps, plot_data = plot_data, gene_averages = gene_averages, gene_averages_cancer_caf = gene_averages_cancer_caf)

	},
	simplify = FALSE,
	USE.NAMES = TRUE
)

cairo_pdf('../data_and_figures/sc_sim_deconv_comp_centred.pdf', width = 6, height = 9, onefile = TRUE)

for(ct in names(sc_deconv_comparison_centred)[!grepl('lenient', names(sc_deconv_comparison_centred))]) {

	sim_deconv_figures <- sapply(
		simulated_deconv_plots[[ct]]$plots[c('purity_bar', 'ccle_bar', 'heatmap')],
		function(x) x + theme(legend.position = 'none', plot.title = element_blank()),
		simplify = FALSE,
		USE.NAMES = TRUE
	)
	sim_deconv_figures$heatmap <- sim_deconv_figures$heatmap + theme(plot.margin = unit(c(0, 5.5, 1.25, 5.5), 'pt'))

	plot_grid(
		plot_grid(
			plotlist = c(
				list(blank_plot() + theme(plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), 'pt')) + labs(title = deconv_plot_args_per_ct[[ct]]$plot_title)),
				sim_deconv_figures,
				# lapply(
					# simulated_deconv_plots[[ct]]$plots[c('purity_bar', 'ccle_bar', 'heatmap')],
					# function(x) x + theme(legend.position = 'none', plot.title = element_blank())
				# ),
				lapply(sc_deconv_comparison_centred[[ct]]$heatmaps, function(x) x + theme(legend.position = 'none'))
			),
			nrow = 6,
			ncol = 1,
			align = 'v',
			rel_heights = c(1, 1, 1, 15, 5, 5)
		),
		plot_grid(
			plotlist = c(
				list(
					blank_plot(),
					get_legend(
						simulated_deconv_plots[[1]]$plots$purity_bar +
							theme(legend.justification = c(0, 1), legend.direction = 'vertical', legend.key.width = unit(10, 'pt'), legend.key.height = unit(10, 'pt')) +
							labs(fill = 'Correlation\nwith purity')
					),
					get_legend(
						simulated_deconv_plots[[1]]$plots$ccle_bar +
							theme(legend.justification = c(0, 1), legend.direction = 'vertical', legend.key.width = unit(10, 'pt'), legend.key.height = unit(10, 'pt')) +
							labs(fill = 'Tumours vs.\ncell lines')
					),
					get_legend(
						simulated_deconv_plots[[1]]$plots$heatmap + theme(legend.justification = c(0, 1))
					),
					get_legend(sc_deconv_comparison_centred[[1]]$heatmaps[[1]] + theme(legend.justification = c(0, 1)))
				)
			),
			nrow = 5,
			ncol = 1,
			rel_heights = c(1, 5, 5, 7, 10)
		),
		nrow = 1,
		ncol = 2,
		rel_widths = c(5, 1)
	) %>% print

}

dev.off()

cairo_pdf('../data_and_figures/sc_sim_deconv_comp_lenient_centred.pdf', width = 6, height = 10, onefile = TRUE)

for(ct in names(sc_deconv_comparison_centred)[grepl('lenient', names(sc_deconv_comparison_centred))]) {

	sim_deconv_lenient_figures <- sapply(
		simulated_deconv_plots[[ct]]$plots[c('purity_bar', 'ccle_bar', 'heatmap')],
		function(x) x + theme(legend.position = 'none', plot.title = element_blank()),
		simplify = FALSE,
		USE.NAMES = TRUE
	)
	sim_deconv_lenient_figures$heatmap <- sim_deconv_lenient_figures$heatmap + theme(plot.margin = unit(c(0, 5.5, 1.25, 5.5), 'pt'))

	plot_grid(
		plot_grid(
			plotlist = c(
				list(blank_plot() + theme(plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), 'pt')) + labs(title = deconv_plot_args_per_ct[[ct]]$plot_title)),
				sim_deconv_lenient_figures,
				# lapply(
					# simulated_deconv_plots[[ct]]$plots[c('purity_bar', 'ccle_bar', 'heatmap')],
					# function(x) x + theme(legend.position = 'none', plot.title = element_blank())
				# ),
				lapply(sc_deconv_comparison_centred[[ct]]$heatmaps, function(x) x + theme(legend.position = 'none'))
			),
			nrow = 6,
			ncol = 1,
			align = 'v',
			rel_heights = c(1, 1, 1, 15, 5, 5)
		),
		plot_grid(
			plotlist = c(
				list(
					blank_plot(),
					get_legend(
						simulated_deconv_plots[[1]]$plots$purity_bar +
							theme(legend.justification = c(0, 1), legend.direction = 'vertical', legend.key.width = unit(10, 'pt'), legend.key.height = unit(10, 'pt')) +
							labs(fill = 'Correlation\nwith purity')
					),
					get_legend(
						simulated_deconv_plots[[1]]$plots$ccle_bar +
							theme(legend.justification = c(0, 1), legend.direction = 'vertical', legend.key.width = unit(10, 'pt'), legend.key.height = unit(10, 'pt')) +
							labs(fill = 'Tumours vs.\ncell lines')
					),
					get_legend(
						simulated_deconv_plots[[1]]$plots$heatmap + theme(legend.justification = c(0, 1))
					),
					get_legend(sc_deconv_comparison_centred[[1]]$heatmaps[[1]] + theme(legend.justification = c(0, 1)))
				)
			),
			nrow = 5,
			ncol = 1,
			rel_heights = c(1, 5, 5, 7, 10)
		),
		nrow = 1,
		ncol = 2,
		rel_widths = c(5, 1)
	) %>% print

}

dev.off()

# Make versions where we filter out genes with low expression in the single cell data:

# Note I previously re-ran the deconvolution with the filtered genes.  The code for doing that (with some comments) is:

# simulated_deconv_ct <- do.call(
	# deconvolve_emt_caf_data,
	# args = c(
		# list(
			# expression_data = simulated_bulk_data,
			# meta_data = simulated_bulk_metadata,
			# genes = sc_deconv_comparison_centred[[ct]]$gene_averages_cancer_caf[
				# ,
				# .(pass = ave_exp[cell_type == 'cancer'] > 0.25 | ave_exp[cell_type == 'caf'] > 0.25),
				# by = gene
			# ][pass == TRUE, as.character(gene)],
			# # genes = sc_deconv_comparison_centred[[ct]]$gene_averages[ave_exp >= 0.5, as.character(symbol)], # symbol is a factor for some reason
			# # genes = gene_averages[
				# # ,
				# # .(pass = ave_exp[cancer_caf == 'cancer'] > 0.5 | ave_exp[cancer_caf == 'caf'] > 0.5),
				# # by = symbol
			# # ][pass == TRUE, as.character(symbol)],
			# cell_type_markers = cell_type_markers,
			# ccle_data = ccle,
			# genes_from_tcga_fun = NULL,
			# genes_filter_fun = NULL
		# ),
		# deconv_args_per_ct[[ct]][!(names(deconv_args_per_ct[[ct]]) %in% c('plot_title', 'genes_filter_fun'))]
	# )
# )

# simulated_deconv_ct <- deconv_reorder(simulated_deconv_ct)

sc_deconv_comparison_filtered <- sapply(
	names(sc_deconv_comparison_centred),
	function(ct) {

		cat(paste0(ct, '\n'))

		simulated_deconv_ct <- simulated_deconvs[[ct]]

		filtered_genes <- sc_deconv_comparison_centred[[ct]]$gene_averages_cancer_caf[
			,
			.(pass = ave_exp[cell_type == 'cancer'] > 0.25 | ave_exp[cell_type == 'caf'] > 0.25),
			by = gene
		][pass == TRUE, as.character(gene)]

		ordered_filtered_genes <- with(simulated_deconv_ct, genes_filtered[ordering][genes_filtered[ordering] %in% filtered_genes])

		simulated_deconv_ct$ordering <- order(filtered_genes)[order(order(ordered_filtered_genes))]
		simulated_deconv_ct$genes_filtered <- filtered_genes
		simulated_deconv_ct$cor_mat <- simulated_deconv_ct$cor_mat[filtered_genes, filtered_genes]
		simulated_deconv_ct$cor_with_purity <- sapply(
			simulated_deconv_ct$cor_with_purity,
			function(x) x[filtered_genes],
			simplify = FALSE,
			USE.NAMES = TRUE
		)
		simulated_deconv_ct$ccle_comp_diff <- simulated_deconv_ct$ccle_comp_diff[filtered_genes]

		simulated_deconv_plots_ct <- do.call(
			deconvolve_emt_caf_plots,
			args = c(
				list(
					data = simulated_deconv_ct,
					heatmap_legend_title = 'Correlation',
					heatmap_colours = rev(colorRampPalette(brewer.pal(11, "RdBu"))(50)),
					heatmap_colour_limits = c(-1, 1),
					heatmap_legend_breaks = c(-1, -0.5, 0, 0.5, 1),
					heatmap_annotations_side = 'left',
					heatmap_annotations_nudge = 0.3,
					purity_colours = rev(colorRampPalette(brewer.pal(11, "PuOr"))(50)),
					purity_colour_limits = c(-1, 1),
					purity_legend_breaks = c(-1, 0, 1),
					purity_legend_title = 'Correlation with purity\n',
					purity_legend_direction = 'horizontal',
					purity_axis_title = NULL,
					ccle_colours = rev(colorRampPalette(brewer.pal(11, 'PiYG')[c(1:3, 6, 9:11)])(50)),
					ccle_fun = function(x) {caTools::runmean(x - mean(x), 30)/max(abs(caTools::runmean(x - mean(x), 30)))},
					ccle_legend_breaks = c(-1, 0, 1),
					ccle_legend_title = 'Tumours vs. cell lines\n',
					ccle_legend_direction = 'horizontal',
					ccle_axis_title = NULL,
					bar_legend_justification = 'left',
					bar_legend_width = unit(10, 'pt'),
					bar_legend_height = unit(10, 'pt')
				),
				deconv_plot_args_per_ct[[ct]]
			)
		)

		plot_data <- sc_deconv_comparison_centred[[ct]]$plot_data[gene %in% simulated_deconv_ct$genes_filtered]

		sc_heatmaps <- sapply(
			c('cancer', 'caf'),
			function(x) {
				ggplot(
					plot_data[cell_type == x],
					aes(
						x = factor(gene, levels = with(simulated_deconv_ct, genes_filtered[ordering])),
						y = factor(id, levels = sc_deconv_comparison[[ct]][[x]]$ordered_cell_ids),
						fill = expression_level
					)
				) +
					geom_raster() +
					scale_x_discrete(expand = c(0, 0)) +
					scale_y_discrete(expand = c(0, 0)) +
					scale_fill_gradientn(
						colours = c(
							colorspace::sequential_hcl(25, h = 190, c = 70, l = c(70, 99), power = 0.8),
							colorspace::sequential_hcl(25, h = 60, c = 100, l = c(80, 99), power = 0.8, rev = TRUE)
						),
						# colours = rev(colorRampPalette(brewer.pal(11, "RdBu"))(50)),
						limits = c(-6, 6),
						oob = scales::squish,
						breaks = c(-6, -3, 0, 3, 6),
						labels = c('-6' = '\u2264 -6', '-3' = '-3', '0' = '0', '3' = '3', '6' = '\u2265 6')
					) +
					theme(
						axis.text = element_blank(),
						axis.title.x = element_blank(),
						axis.ticks = element_blank(),
						axis.ticks.length = unit(0, 'pt'),
						plot.margin = unit(c(5.5, 5.5, 1.25, 5.5), 'pt')
					) +
					labs(y = mapvalues(x, c('cancer', 'caf'), c('Cancer', 'CAF'), warn_missing = FALSE), fill = 'Relative\nexpression\nlevel')
			},
			simplify = FALSE,
			USE.NAMES = TRUE
		)

		list(deconv_data = simulated_deconv_ct, deconv_figures = simulated_deconv_plots_ct$plots, sc_heatmaps = sc_heatmaps, plot_data = plot_data)

	},
	simplify = FALSE,
	USE.NAMES = TRUE
)

# plot_data <- sc_deconv_comparison_centred[[ct]]$plot_data[
	# gene %in% sc_deconv_comparison_centred[[ct]]$gene_averages_cancer_caf[
		# ,
		# .(pass = ave_exp[cell_type == 'cancer'] > 0.5 | ave_exp[cell_type == 'caf'] > 0.5),
		# by = gene
	# ][pass == TRUE, as.character(gene)]
	# # gene %in% sc_deconv_comparison_centred[[ct]]$gene_averages[ave_exp >= 0.5, symbol]
# ]

# In the above, use the following if you want to recentre the genes and cells:

# plot_data <- rbind(
	# sc_deconv_comparison[[ct]]$cancer$plot_data[, cell_type := 'cancer'],
	# sc_deconv_comparison[[ct]]$caf$plot_data[, cell_type := 'caf']
# )[gene %in% sc_deconv_comparison_centred[[ct]]$gene_averages[ave_exp >= 0.5, symbol]]

# # To centre genes w.r.t. the average of the averages of cancer and CAF:
# gene_averages <- plot_data[
	# ,
	# .(ave_exp = mean(c(mean(expression_level[cell_type == 'cancer']), mean(expression_level[cell_type == 'caf'])))),
	# by = .(symbol = gene)
# ]

# plot_data[, expression_level := expression_level - gene_averages[symbol == unique(gene), ave_exp], by = gene]

# # To centre the cells as well:
# plot_data[, expression_level := expression_level - mean(expression_level), by = id]

cairo_pdf('../data_and_figures/sc_sim_deconv_comp_filtered.pdf', width = 6, height = 9, onefile = TRUE)

for(ct in names(sc_deconv_comparison_filtered)[!grepl('lenient', names(sc_deconv_comparison_filtered))]) {

	sc_sim_deconv_comp_figures <- sapply(
		c(
			sc_deconv_comparison_filtered[[ct]]$deconv_figures[c('purity_bar', 'ccle_bar', 'heatmap')],
			sc_deconv_comparison_filtered[[ct]]$sc_heatmaps
		),
		function(x) x + theme(legend.position = 'none', plot.title = element_blank()),
		simplify = FALSE,
		USE.NAMES = TRUE
	)
	sc_sim_deconv_comp_figures$heatmap <- sc_sim_deconv_comp_figures$heatmap + theme(plot.margin = unit(c(0, 5.5, 1.25, 5.5), 'pt'))

	plot_grid(
		plot_grid(
			plotlist = c(
				list(blank_plot() + theme(plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), 'pt')) + labs(title = deconv_plot_args_per_ct[[ct]]$plot_title)),
				sc_sim_deconv_comp_figures
				# lapply(
					# c(
						# sc_deconv_comparison_filtered[[ct]]$deconv_figures[c('purity_bar', 'ccle_bar', 'heatmap')],
						# sc_deconv_comparison_filtered[[ct]]$sc_heatmaps
					# ),
					# function(x) x + theme(legend.position = 'none', plot.title = element_blank())
				# )
			),
			nrow = 6,
			ncol = 1,
			align = 'v',
			rel_heights = c(1, 1, 1, 15, 5, 5)
		),
		plot_grid(
			plotlist = c(
				list(
					blank_plot(),
					get_legend(
						simulated_deconv_plots[[1]]$plots$purity_bar +
							theme(legend.justification = c(0, 1), legend.direction = 'vertical') +#, legend.key.width = NULL, legend.key.height = NULL) +
							labs(fill = 'Correlation\nwith purity')
					),
					get_legend(
						simulated_deconv_plots[[1]]$plots$ccle_bar +
							theme(legend.justification = c(0, 1), legend.direction = 'vertical') +#, legend.key.width = NULL, legend.key.height = NULL) +
							labs(fill = 'Tumours vs.\ncell lines')
					),
					get_legend(
						simulated_deconv_plots[[1]]$plots$heatmap + theme(legend.justification = c(0, 1))
					),
					get_legend(sc_deconv_comparison_filtered[[1]]$sc_heatmaps[[1]] + theme(legend.justification = c(0, 1)))
				)
			),
			nrow = 5,
			ncol = 1,
			rel_heights = c(1, 5, 5, 7, 10)
		),
		nrow = 1,
		ncol = 2,
		rel_widths = c(5, 1)
	) %>% print

}

dev.off()

cairo_pdf('../data_and_figures/sc_sim_deconv_comp_lenient_filtered.pdf', width = 6, height = 10, onefile = TRUE)

for(ct in names(sc_deconv_comparison_filtered)[grepl('lenient', names(sc_deconv_comparison_filtered))]) {

	sc_sim_deconv_comp_lenient_figures <- sapply(
		c(
			sc_deconv_comparison_filtered[[ct]]$deconv_figures[c('purity_bar', 'ccle_bar', 'heatmap')],
			sc_deconv_comparison_filtered[[ct]]$sc_heatmaps
		),
		function(x) x + theme(legend.position = 'none', plot.title = element_blank()),
		simplify = FALSE,
		USE.NAMES = TRUE
	)
	sc_sim_deconv_comp_lenient_figures$heatmap <- sc_sim_deconv_comp_lenient_figures$heatmap + theme(plot.margin = unit(c(0, 5.5, 1.25, 5.5), 'pt'))

	plot_grid(
		plot_grid(
			plotlist = c(
				list(blank_plot() + theme(plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), 'pt')) + labs(title = deconv_plot_args_per_ct[[ct]]$plot_title)),
				sc_sim_deconv_comp_lenient_figures
				# lapply(
					# c(
						# sc_deconv_comparison_filtered[[ct]]$deconv_figures[c('purity_bar', 'ccle_bar', 'heatmap')],
						# sc_deconv_comparison_filtered[[ct]]$sc_heatmaps
					# ),
					# function(x) x + theme(legend.position = 'none', plot.title = element_blank())
				# )
			),
			nrow = 6,
			ncol = 1,
			align = 'v',
			rel_heights = c(1, 1, 1, 15, 5, 5)
		),
		plot_grid(
			plotlist = c(
				list(
					blank_plot(),
					get_legend(
						simulated_deconv_plots[[1]]$plots$purity_bar +
							theme(legend.justification = c(0, 1), legend.direction = 'vertical') +#, legend.key.width = NULL, legend.key.height = NULL) +
							labs(fill = 'Correlation\nwith purity')
					),
					get_legend(
						simulated_deconv_plots[[1]]$plots$ccle_bar +
							theme(legend.justification = c(0, 1), legend.direction = 'vertical') +#, legend.key.width = NULL, legend.key.height = NULL) +
							labs(fill = 'Tumours vs.\ncell lines')
					),
					get_legend(
						simulated_deconv_plots[[1]]$plots$heatmap + theme(legend.justification = c(0, 1))
					),
					get_legend(sc_deconv_comparison_filtered[[1]]$sc_heatmaps[[1]] + theme(legend.justification = c(0, 1)))
				)
			),
			nrow = 5,
			ncol = 1,
			rel_heights = c(1, 5, 5, 7, 10)
		),
		nrow = 1,
		ncol = 2,
		rel_widths = c(5, 1)
	) %>% print

}

dev.off()





deconv_args_per_ct <- list(

    brca_luminal_a = list(
        tcga_cancer_types = 'BRCA',
        subtypes = 'Luminal A',
        ref_for_subtypes = 'doi:10.1038/nature11412',
        ccle_cancer_type = 'breast',
        extra_data_source = 'breast_qian',
        seed = 8216,
        plot_title = 'BRCA - Luminal A',
        genes_filter_fun = function(x) 1:175,
        gene_weights_fun = function(x) quantile(x, 0.9)
    ),

    brca_luminal_b = list(
        tcga_cancer_types = 'BRCA',
        subtypes = 'Luminal B',
        ref_for_subtypes = 'doi:10.1038/nature11412',
        ccle_cancer_type = 'breast',
        extra_data_source = 'breast_qian',
        seed = 9056,
        plot_title = 'BRCA - Luminal B',
        genes_filter_fun = function(x) 1:230,
        gene_weights_fun = function(x) quantile(x, 0.75)
    ),

    brca_basal_like = list(
        tcga_cancer_types = 'BRCA',
        subtypes = 'Basal-like',
        ref_for_subtypes = 'doi:10.1038/nature11412',
        ccle_cancer_type = 'breast',
        extra_data_source = 'breast_qian',
        seed = 1982,
        plot_title = 'BRCA - Basal-like',
        genes_filter_fun = function(x) 1:230,
        gene_weights_fun = function(x) quantile(x, 0.75)
    ),

    brca_her2_enriched = list(
        tcga_cancer_types = 'BRCA',
        subtypes = 'HER2-enriched',
        ref_for_subtypes = 'doi:10.1038/nature11412',
        ccle_cancer_type = 'breast',
        extra_data_source = 'breast_qian',
        seed = 6994,
        plot_title = 'BRCA - HER2-enriched',
        genes_filter_fun = function(x) 1:180,
        gene_weights_fun = function(x) quantile(x, 0.75)
    ),

    coad = list(
        tcga_cancer_types = 'COAD',
        ccle_cancer_type = 'large intestine',
        extra_data_source = 'crc_lee_smc',
        seed = 3260,
        plot_title = 'COAD',
        genes_filter_fun = function(x) 1:200,
        cell_type_weights = list(
            B_plasma = 0,
            myocyte = 0,
            macrophage = 0,
            endothelial = 1,
            DC = 0,
            mast = 0,
            T = 0,
            B = 0
        )
    ),

    hnsc_mesenchymal_basal = list(
        tcga_cancer_types = 'HNSC',
        subtypes = c('Mesenchymal', 'Basal'),
        ref_for_subtypes = 'doi:10.1038/nature14129',
        ccle_cancer_type = 'HNSCC',
        extra_data_source = 'hnscc_puram',
        seed = 7466,
        plot_title = 'HNSC - Mesenchymal & Basal'
    ),

    hnsc_classical = list(
        tcga_cancer_types = 'HNSC',
        subtypes = 'Classical',
        ref_for_subtypes = 'doi:10.1038/nature14129',
        ccle_cancer_type = 'HNSCC',
        extra_data_source = 'hnscc_puram',
        seed = 9998,
        plot_title = 'HNSC - Classical',
        gene_weights_fun = function(x) quantile(x, 0.75)
    ),

    hnsc_atypical = list(
        tcga_cancer_types = 'HNSC',
        subtypes = 'Atypical',
        ref_for_subtypes = 'doi:10.1038/nature14129',
        ccle_cancer_type = 'HNSCC',
        extra_data_source = 'hnscc_puram',
        seed = 8172,
        plot_title = 'HNSC - Atypical',
        gene_weights_fun = function(x) quantile(x, 0.9)
    ),

    lihc = list(
        tcga_cancer_types = 'LIHC',
        ccle_cancer_type = 'liver',
        extra_data_source = 'liver_ma',
        seed = 9349,
        plot_title = 'LIHC',
        genes_filter_fun = function(x) 1:260,
        cell_type_weights = FALSE
    ),

    luad_proximal_inflammatory = list(
        tcga_cancer_types = 'LUAD',
        subtypes = 'Proximal-inflammatory',
        ref_for_subtypes = 'doi:10.1038/nature13385',
        ccle_cancer_type = 'lung',
		extra_data_source = 'luad_kim',
        # extra_data_source = 'luad_qian',
        seed = 9955,
        plot_title = 'LUAD - Proximal-inflammatory'
    ),

    luad_proximal_proliferative = list(
        tcga_cancer_types = 'LUAD',
        subtypes = 'Proximal-proliferative',
        ref_for_subtypes = 'doi:10.1038/nature13385',
        ccle_cancer_type = 'lung',
		extra_data_source = 'luad_kim',
        # extra_data_source = 'luad_qian',
        seed = 743,
        plot_title = 'LUAD - Proximal-proliferative',
        genes_filter_fun = function(x) 1:225,
        gene_weights_fun = function(x) quantile(x, 0.75)
    ),

    lusc_basal = list(
        tcga_cancer_types = 'LUSC',
        subtypes = 'basal',
        ref_for_subtypes = 'doi:10.1038/nature11404',
        ccle_cancer_type = 'lung',
        extra_data_source = 'lusc_qian',
        seed = 2287,
        plot_title = 'LUSC - Basal',
        gene_weights_fun = function(x) quantile(x, 0.9)
    ),

    lusc_classical = list(
        tcga_cancer_types = 'LUSC',
        subtypes = 'classical',
        ref_for_subtypes = 'doi:10.1038/nature11404',
        ccle_cancer_type = 'lung',
        extra_data_source = 'lusc_qian',
        seed = 7106,
        plot_title = 'LUSC - Classical',
        gene_weights_fun = function(x) quantile(x, 0.9)
    ),

    lusc_primitive = list(
        tcga_cancer_types = 'LUSC',
        subtypes = 'primitive',
        ref_for_subtypes = 'doi:10.1038/nature11404',
        ccle_cancer_type = 'lung',
        extra_data_source = 'lusc_qian',
        seed = 303,
        plot_title = 'LUSC - Primitive',
        gene_weights_fun = function(x) quantile(x, 0.9)
    ),

    lusc_secretory = list(
        tcga_cancer_types = 'LUSC',
        subtypes = 'secretory',
        ref_for_subtypes = 'doi:10.1038/nature11404',
        ccle_cancer_type = 'lung',
        extra_data_source = 'lusc_qian',
        seed = 9240,
        plot_title = 'LUSC - Secretory'
    ),

    ov_differentiated = list(
        tcga_cancer_types = 'OV',
        subtypes = 'Differentiated',
        ref_for_subtypes = 'doi:10.1172/JCI65833',
        ccle_cancer_type = 'ovary',
		extra_data_source = 'ovarian_qian',
        seed = 1498,
        plot_title = 'OV - Differentiated',
        genes_filter_fun = function(x) 1:270,
        gene_weights_fun = function(x) quantile(x, 0.9)
    ),

    ov_immunoreactive = list(
        tcga_cancer_types = 'OV',
        subtypes = 'Immunoreactive',
        ref_for_subtypes = 'doi:10.1172/JCI65833',
        ccle_cancer_type = 'ovary',
		extra_data_source = 'ovarian_qian',
        seed = 991,
        plot_title = 'OV - Immunoreactive',
        genes_filter_fun = function(x) 1:225,
        gene_weights_fun = function(x) quantile(x, 0.9)
    ),

    ov_proliferative = list(
        tcga_cancer_types = 'OV',
        subtypes = 'Proliferative',
        ref_for_subtypes = 'doi:10.1172/JCI65833',
        ccle_cancer_type = 'ovary',
		extra_data_source = 'ovarian_qian',
        seed = 8583,
        plot_title = 'OV - Proliferative',
        gene_weights_fun = function(x) quantile(x, 0.9)
    ),

    paad = list(
        tcga_cancer_types = 'PAAD',
        ccle_cancer_type = 'pancreas',
        extra_data_source = 'pdac_peng',
        seed = 88,
        plot_title = 'PAAD',
        genes_filter_fun = function(x) 1:150, # Also works OK with 175, but 150 is better
        cell_type_weights = list(
            B_plasma = 0,
            myocyte = 0,
            macrophage = 0,
            endothelial = 1,
            DC = 0,
            mast = 0,
            T = 0,
            B = 0
        )
    ),

    read = list(
        tcga_cancer_types = 'READ',
        ccle_cancer_type = 'large intestine',
        seed = 3309,
        plot_title = 'READ',
        extra_data_source = 'crc_lee_smc'
    )

)

set.seed(8542)

sc_tcga_deconv_comparison <- sapply(
	ct_to_keep,
	function(ct) {
		cat(paste0(ct, '\n'))
        sc_data <- eval(sc_metadata[[tcga_sc_map[ct]]]$read_quote)
		sapply(
			c('cancer', 'caf'),
			function(x) {

				plot_data <- copy(sc_data[cell_type == x])[
					,
					complexity := apply(.SD, 1, function(x) sum(x > 0)),
					.SDcols = -c('id', 'patient', 'cell_type')
				][
					,
					.SD[sample(1:.N, sample_sizes[[tcga_sc_map[ct]]], prob = complexity)]
				][
					,
					complexity := NULL
				][
					,
					c('id', deconv_data[[ct]]$genes_filtered[deconv_data[[ct]]$genes_filtered %in% names(sc_data)]),
					with = FALSE
				]

				plot_data <- melt(plot_data, id.var = 'id', variable.name = 'gene', value.name = 'expression_level')

				ordered_cell_ids <- plot_data[
					gene %in% with(
						deconv_data[[ct]],
						do.call(c('head', 'tail')[which(c('cancer', 'caf') == x)], list(genes_filtered[ordering], 20))
					),
					.(top_20_mean = mean(expression_level)),
					by = id
				][order(top_20_mean), id]

				htmp <- ggplot(
					plot_data,
					aes(
						x = factor(gene, levels = with(deconv_data[[ct]], genes_filtered[ordering])),
						y = factor(id, levels = ordered_cell_ids),
						fill = expression_level
					)
				) +
					geom_raster() +
					scale_x_discrete(expand = c(0, 0)) +
					scale_y_discrete(expand = c(0, 0)) +
					scale_fill_gradientn(
						colours = rev(colorRampPalette(brewer.pal(11, "Spectral"))(50)),
						limits = c(0, 12),
						oob = scales::squish,
						breaks = c(0, 3, 6, 9, 12),
						labels = c('0' = '0', '3' = '3', '6' = '6', '9' = '9', '12' = '\u2265 12')
					) +
					theme(
						axis.text = element_blank(),
						axis.title.x = element_blank(),
						axis.ticks = element_blank(),
						axis.ticks.length = unit(0, 'pt'),
						plot.margin = unit(c(5.5, 5.5, 1.25, 1.1), 'pt')
					) +
					labs(
						y = mapvalues(x, c('cancer', 'caf'), c('Cancer', 'CAF'), warn_missing = FALSE),
						fill = 'Expression\nlevel'
					)

				ave_exp_bar <- ggplot(
					plot_data[, .(ave_exp = mean(expression_level)), keyby = gene][
						with(deconv_data[[ct]], genes_filtered[ordering]),
						ave_exp := runmean(ave_exp, 30)
					],
					aes(x = factor(gene, levels = with(deconv_data[[ct]], genes_filtered[ordering])), y = 'a', fill = ave_exp)
				) +
					geom_raster() +
					scale_x_discrete(expand = c(0, 0)) +
					scale_y_discrete(expand = c(0, 0)) +
					scale_fill_gradientn(
						colours = colorRampPalette(brewer.pal(9, "RdPu"))(50),
						limits = c(0, 4),
						oob = scales::squish,
						breaks = c(0, 2, 4),
						labels = c('0' = '0', '2' = '2', '4' = '\u2265 4')
					) +
					theme(
						axis.text = element_blank(),
						axis.title = element_blank(),
						axis.ticks = element_blank(),
						axis.ticks.length = unit(0, 'pt'),
						plot.margin = switch(
							(x == 'cancer') + 1,
							unit(c(1.25, 5.5, 5.5, 5.5), 'pt'),
							unit(c(1.25, 5.5, 1.25, 5.5), 'pt')
						)
					) +
					labs(
						y = mapvalues(x, c('cancer', 'caf'), c('Cancer', 'CAF'), warn_missing = FALSE),
						fill = 'Average\nexpression\nlevel'
					)

				list(heatmap = htmp, ave_exp_bar = ave_exp_bar, plot_data = plot_data, ordered_cell_ids = ordered_cell_ids)

			},
			simplify = FALSE,
			USE.NAMES = TRUE
		)
	},
	simplify = FALSE,
	USE.NAMES = TRUE
)

cairo_pdf('../data_and_figures/sc_deconv_comp.pdf', width = 6, height = 10, onefile = TRUE)

for(ct in names(sc_tcga_deconv_comparison)) {

	plot_grid(
		plot_grid(
			plotlist = c(
				list(blank_plot() + theme(plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), 'pt')) + labs(title = deconv_args_per_ct[[ct]]$plot_title)),
				lapply(
					deconv_plots[[ct]]$plots[c('purity_bar', 'ccle_bar', 'heatmap')],
					function(x) x + theme(legend.position = 'none', plot.title = element_blank())
				),
				unlist(
					lapply(
						sc_tcga_deconv_comparison[[ct]],
						function(x) list(x$heatmap + theme(legend.position = 'none'), x$ave_exp_bar + theme(legend.position = 'none'))
					),
					recursive = FALSE
				)
			),
			nrow = 8,
			ncol = 1,
			align = 'v',
			rel_heights = c(1, 1, 1, 15, 5, 1, 5, 1)
		),
		plot_grid(
			plotlist = c(
				list(
					blank_plot(),
					get_legend(
						deconv_plots[[1]]$plots$purity_bar +
							theme(legend.justification = c(0, 1), legend.direction = 'vertical', legend.key.width = NULL, legend.key.height = NULL) +
							labs(fill = 'Correlation\nwith purity')
					),
					get_legend(
						deconv_plots[[1]]$plots$ccle_bar +
							theme(legend.justification = c(0, 1), legend.direction = 'vertical', legend.key.width = NULL, legend.key.height = NULL) +
							labs(fill = 'Tumours vs.\ncell lines')
					),
					get_legend(
						deconv_plots[[1]]$plots$heatmap + theme(legend.justification = c(0, 1))
					)
				),
				list(
					get_legend(sc_tcga_deconv_comparison[[1]][[1]]$heatmap + theme(legend.justification = c(0, 1))),
					get_legend(sc_tcga_deconv_comparison[[1]][[1]]$ave_exp_bar + theme(legend.justification = c(0, 1)))
				)
			),
			nrow = 6,
			ncol = 1,
			rel_heights = c(1, 5, 5, 7, 6, 6)
		),
		nrow = 1,
		ncol = 2,
		rel_widths = c(5, 1)
	) %>% print

}

dev.off()

# Make centred versions:

sc_tcga_deconv_comparison_centred <- sapply(
	names(sc_tcga_deconv_comparison),
	function(ct) {

		cat(paste0(ct, '\n'))

		plot_data <- rbind(
			sc_tcga_deconv_comparison[[ct]]$cancer$plot_data[, cell_type := 'cancer'],
			sc_tcga_deconv_comparison[[ct]]$caf$plot_data[, cell_type := 'caf']
		)

		# Calculate averages in cancer cells and CAFs separately, so we can check them afterwards:
		gene_averages_cancer_caf <- plot_data[
			,
			.(ave_exp = mean(expression_level)),
			by = .(gene, cell_type)
		]

		# To centre genes w.r.t. the average of the averages of cancer and CAF:
		gene_averages <- plot_data[
			,
			.(ave_exp = mean(c(mean(expression_level[cell_type == 'cancer']), mean(expression_level[cell_type == 'caf'])))),
			by = .(symbol = gene)
		]

		plot_data[, expression_level := expression_level - gene_averages[symbol == unique(gene), ave_exp], by = gene]

		# To centre the cells:
		plot_data[, expression_level := expression_level - mean(expression_level), by = id]

		htmps <- sapply(
			c('cancer', 'caf'),
			function(x) {
				ggplot(
					plot_data[cell_type == x],
					aes(
						x = factor(gene, levels = with(deconv_data[[ct]], genes_filtered[ordering])),
						y = factor(id, levels = sc_tcga_deconv_comparison[[ct]][[x]]$ordered_cell_ids),
						fill = expression_level
					)
				) +
					geom_raster() +
					scale_x_discrete(expand = c(0, 0)) +
					scale_y_discrete(expand = c(0, 0)) +
					scale_fill_gradientn(
						colours = c(
							colorspace::sequential_hcl(25, h = 190, c = 70, l = c(70, 99), power = 0.8),
							colorspace::sequential_hcl(25, h = 60, c = 100, l = c(80, 99), power = 0.8, rev = TRUE)
						),
						# colours = rev(colorRampPalette(brewer.pal(11, "RdBu"))(50)),
						limits = c(-6, 6),
						# limits = c(-2, 2),
						oob = scales::squish,
						breaks = c(-6, -3, 0, 3, 6),
						labels = c('-6' = '\u2264 -6', '-3' = '-3', '0' = '0', '3' = '3', '6' = '\u2265 6')
						# breaks = c(-2, -1, 0, 1, 2),
						# labels = c('-2' = '\u2264 -2', '-1' = '-1', '0' = '0', '1' = '1', '2' = '\u2265 2')
					) +
					theme(
						axis.text = element_blank(),
						axis.title.x = element_blank(),
						axis.ticks = element_blank(),
						axis.ticks.length = unit(0, 'pt'),
						plot.margin = unit(c(5.5, 5.5, 1.25, 1.1), 'pt')
						# plot.margin = switch((x == 'cancer') + 1, unit(c(1.25, 5.5, 5.5, 5.5), 'pt'), unit(c(5.5, 5.5, 1.25, 5.5), 'pt'))
					) +
					labs(
						y = mapvalues(x, c('cancer', 'caf'), c('Cancer', 'CAF'), warn_missing = FALSE),
						fill = 'Relative\nexpression\nlevel'
					)
			},
			simplify = FALSE,
			USE.NAMES = TRUE
		)

		list(heatmaps = htmps, plot_data = plot_data, gene_averages = gene_averages, gene_averages_cancer_caf = gene_averages_cancer_caf)

	},
	simplify = FALSE,
	USE.NAMES = TRUE
)

cairo_pdf('../data_and_figures/sc_deconv_comp_centred.pdf', width = 6, height = 9, onefile = TRUE)

for(ct in names(sc_tcga_deconv_comparison_centred)) {

	plot_grid(
		plot_grid(
			plotlist = c(
				list(blank_plot() + theme(plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), 'pt')) + labs(title = deconv_args_per_ct[[ct]]$plot_title)),
				lapply(
					deconv_plots[[ct]]$plots[c('purity_bar', 'ccle_bar', 'heatmap')],
					function(x) x + theme(legend.position = 'none', plot.title = element_blank())
				),
				lapply(sc_tcga_deconv_comparison_centred[[ct]]$heatmaps, function(x) x + theme(legend.position = 'none'))
			),
			nrow = 6,
			ncol = 1,
			align = 'v',
			rel_heights = c(1, 1, 1, 15, 5, 5)
		),
		plot_grid(
			plotlist = c(
				list(
					blank_plot(),
					get_legend(
						deconv_plots[[1]]$plots$purity_bar +
							theme(legend.justification = c(0, 1), legend.direction = 'vertical', legend.key.width = NULL, legend.key.height = NULL) +
							labs(fill = 'Correlation\nwith purity')
					),
					get_legend(
						deconv_plots[[1]]$plots$ccle_bar +
							theme(legend.justification = c(0, 1), legend.direction = 'vertical', legend.key.width = NULL, legend.key.height = NULL) +
							labs(fill = 'Tumours vs.\ncell lines')
					),
					get_legend(
						deconv_plots[[1]]$plots$heatmap + theme(legend.justification = c(0, 1))
					),
					get_legend(sc_tcga_deconv_comparison_centred[[1]]$heatmaps[[1]] + theme(legend.justification = c(0, 1)))
				)
			),
			nrow = 5,
			ncol = 1,
			rel_heights = c(1, 5, 5, 7, 10)
		),
		nrow = 1,
		ncol = 2,
		rel_widths = c(5, 1)
	) %>% print

}

dev.off()

# Make versions where we filter out genes with low expression in the single cell data:

sc_tcga_deconv_comparison_filtered <- sapply(
	names(sc_tcga_deconv_comparison_centred),
	function(ct) {

		cat(paste0(ct, '\n'))

		deconv_ct <- do.call(
			deconvolve_emt_caf_data,
			args = c(
				list(
					expression_data = expression_data,
					meta_data = meta_data,
					genes = sc_tcga_deconv_comparison_centred[[ct]]$gene_averages_cancer_caf[
						,
						.(pass = ave_exp[cell_type == 'cancer'] > 0.5 | ave_exp[cell_type == 'caf'] > 0.5),
						by = gene
					][pass == TRUE, as.character(gene)],
					# genes = sc_tcga_deconv_comparison_centred[[ct]]$gene_averages[ave_exp >= 0.5, as.character(symbol)], # symbol is a factor for some reason
					cell_type_markers = cell_type_markers,
					ccle_data = ccle,
					subtypes_data = subtypes_data,
					genes_from_tcga_fun = NULL,
					genes_filter_fun = NULL
				),
				deconv_args_per_ct[[ct]][!(names(deconv_args_per_ct[[ct]]) %in% c('plot_title', 'genes_filter_fun'))]
			)
		)

		deconv_ct <- deconv_reorder(deconv_ct)

		deconv_plots_ct <- do.call(
			deconvolve_emt_caf_plots,
			args = c(
				list(
					data = deconv_ct,
					heatmap_legend_title = 'Correlation',
					heatmap_colours = rev(colorRampPalette(brewer.pal(11, "RdBu"))(50)),
					heatmap_colour_limits = c(-0.6, 0.6),
					heatmap_legend_breaks = c(-0.6, -0.3, 0, 0.3, 0.6),
					heatmap_annotations_nudge = 0.3,
					purity_colours = rev(colorRampPalette(brewer.pal(11, "PuOr"))(50)),
					purity_colour_limits = c(-0.3, 0.3),
					purity_legend_breaks = c(-0.3, 0, 0.3),
					purity_legend_title = 'Correlation with purity\n',
					purity_legend_direction = 'horizontal',
					purity_axis_title = NULL,
					ccle_colours = rev(colorRampPalette(brewer.pal(11, 'PiYG')[c(1:3, 6, 9:11)])(50)),
					ccle_legend_breaks = c(-1, 0, 1),
					ccle_legend_title = 'Tumours vs. cell lines\n',
					ccle_legend_direction = 'horizontal',
					ccle_axis_title = NULL,
					bar_legend_justification = 'left'
				),
				deconv_args_per_ct[[ct]]['plot_title']
			)
		)

		plot_data <- sc_tcga_deconv_comparison_centred[[ct]]$plot_data[
			gene %in% sc_tcga_deconv_comparison_centred[[ct]]$gene_averages_cancer_caf[
				,
				.(pass = ave_exp[cell_type == 'cancer'] > 0.5 | ave_exp[cell_type == 'caf'] > 0.5),
				by = gene
			][pass == TRUE, as.character(gene)]
			# gene %in% sc_tcga_deconv_comparison_centred[[ct]]$gene_averages[ave_exp >= 0.5, symbol]
		]

		sc_heatmaps <- sapply(
			c('cancer', 'caf'),
			function(x) {
				ggplot(
					plot_data[cell_type == x],
					aes(
						x = factor(gene, levels = with(deconv_ct, genes_filtered[ordering])),
						y = factor(id, levels = sc_tcga_deconv_comparison[[ct]][[x]]$ordered_cell_ids),
						fill = expression_level
					)
				) +
					geom_raster() +
					scale_x_discrete(expand = c(0, 0)) +
					scale_y_discrete(expand = c(0, 0)) +
					scale_fill_gradientn(
						colours = c(
							colorspace::sequential_hcl(25, h = 190, c = 70, l = c(70, 99), power = 0.8),
							colorspace::sequential_hcl(25, h = 60, c = 100, l = c(80, 99), power = 0.8, rev = TRUE)
						),
						# colours = rev(colorRampPalette(brewer.pal(11, "RdBu"))(50)),
						limits = c(-6, 6),
						oob = scales::squish,
						breaks = c(-6, -3, 0, 3, 6),
						labels = c('-6' = '\u2264 -6', '-3' = '-3', '0' = '0', '3' = '3', '6' = '\u2265 6')
					) +
					theme(
						axis.text = element_blank(),
						axis.title.x = element_blank(),
						axis.ticks = element_blank(),
						axis.ticks.length = unit(0, 'pt'),
						plot.margin = unit(c(5.5, 5.5, 1.25, 5.5), 'pt')
					) +
					labs(y = mapvalues(x, c('cancer', 'caf'), c('Cancer', 'CAF'), warn_missing = FALSE), fill = 'Relative\nexpression\nlevel')
			},
			simplify = FALSE,
			USE.NAMES = TRUE
		)

		list(deconv_data = deconv_ct, deconv_figures = deconv_plots_ct$plots, sc_heatmaps = sc_heatmaps, plot_data = plot_data)

	},
	simplify = FALSE,
	USE.NAMES = TRUE
)

cairo_pdf('../data_and_figures/sc_deconv_comp_filtered.pdf', width = 6, height = 9, onefile = TRUE)

for(ct in names(sc_tcga_deconv_comparison_filtered)) {

	plot_grid(
		plot_grid(
			plotlist = c(
				list(blank_plot() + theme(plot.margin = unit(c(5.5, 5.5, 5.5, 5.5), 'pt')) + labs(title = deconv_args_per_ct[[ct]]$plot_title)),
				lapply(
					c(
						sc_tcga_deconv_comparison_filtered[[ct]]$deconv_figures[c('purity_bar', 'ccle_bar', 'heatmap')],
						sc_tcga_deconv_comparison_filtered[[ct]]$sc_heatmaps
					),
					function(x) x + theme(legend.position = 'none', plot.title = element_blank())
				)
			),
			nrow = 6,
			ncol = 1,
			align = 'v',
			rel_heights = c(1, 1, 1, 15, 5, 5)
		),
		plot_grid(
			plotlist = c(
				list(
					blank_plot(),
					get_legend(
						deconv_plots[[1]]$plots$purity_bar +
							theme(legend.justification = c(0, 1), legend.direction = 'vertical', legend.key.width = NULL, legend.key.height = NULL) +
							labs(fill = 'Correlation\nwith purity')
					),
					get_legend(
						deconv_plots[[1]]$plots$ccle_bar +
							theme(legend.justification = c(0, 1), legend.direction = 'vertical', legend.key.width = NULL, legend.key.height = NULL) +
							labs(fill = 'Tumours vs.\ncell lines')
					),
					get_legend(
						deconv_plots[[1]]$plots$heatmap + theme(legend.justification = c(0, 1))
					),
					get_legend(sc_tcga_deconv_comparison_filtered[[1]]$sc_heatmaps[[1]] + theme(legend.justification = c(0, 1)))
				)
			),
			nrow = 5,
			ncol = 1,
			rel_heights = c(1, 5, 5, 7, 10)
		),
		nrow = 1,
		ncol = 2,
		rel_widths = c(5, 1)
	) %>% print

}

dev.off()





# Another attempt at rare vs. shared EMT:

deconv_data <- readRDS('../data_and_figures/deconv_data.rds')

ct_to_keep <- c(
    'blca_luminal_papillary',
    'blca_basal_squamous',
    'brca_luminal_a',
    'brca_luminal_b',
    'brca_basal_like',
    'brca_her2_enriched',
    'cesc',
    'coad',
    'esca_ac',
    'esca_escc',
    'hnsc_mesenchymal_basal',
    'hnsc_classical',
    'hnsc_atypical',
    'lihc',
    'luad_proximal_inflammatory',
    'luad_proximal_proliferative',
    'lusc_basal',
    'lusc_classical',
    'lusc_primitive',
    'lusc_secretory',
    'ov_differentiated',
    'ov_immunoreactive',
    'ov_proliferative',
    'paad',
    'read',
    # 'skcm_immune',
    # 'skcm_keratin',
    # 'skcm_mitf_low',
    'stad_cin',
    'stad_ebv',
    'stad_msi',
    'ucec'
)

deconv_data <- deconv_data[ct_to_keep]

nice_names_for_figure <- c(
    'BLCA - Luminal-Papillary',
    'BLCA - Basal-Squamous',
    'BRCA - Luminal A',
    'BRCA - Luminal B',
    'BRCA - Basal-like',
    'BRCA - HER2-enriched',
    'CESC',
    'COAD',
    'ESCA - Adenocarcinoma',
    'ESCA - Squamous',
    'HNSC - Malignant-Basal',
    'HNSC - Classical',
    'HNSC - Atypical',
    'LIHC',
    'LUAD - Squamoid',
    'LUAD - Magnoid',
    'LUSC - Basal',
    'LUSC - Classical',
    'LUSC - Primitive',
    'LUSC - Secretory',
    'OV - Differentiated',
    'OV - Immunoreactive',
    'OV - Proliferative',
    'PAAD',
    'READ',
    # 'SKCM - Immune',
    # 'SKCM - Keratin',
    # 'SKCM - MITF-low',
    'STAD - CIN',
    'STAD - EBV',
    'STAD - MSI',
    'UCEC'
)

names(deconv_data) <- mapvalues(
    names(deconv_data),
    names(deconv_data),
    nice_names_for_figure
)





# The following takes a long time, so read from RDS if already done:

# inter_intra_emt_all_cts <- readRDS('../data_and_figures/inter_intra_emt_all_cts.rds')

inter_intra_emt_all_cts <- sapply(
    # names(sc_metadata),
    c('brca', 'coadread', 'hnsc', 'lihc', 'luad', 'paad'),
    function(ct) {

        sc_data <- eval(sc_metadata[[ct]]$read_quote)
        setkey(sc_data, id)

        # Take subset of tumours with at least 50 cancer cells:
        sc_data <- sc_data[patient %in% sc_data[cell_type == 'cancer', .(N = .N), by = patient][N >= 50, patient]]

        # Get all the sufficiently highly-expressed genes:
        all_genes_filtered <- sc_data[
            cell_type == 'cancer',
            names(.SD)[apply(.SD, 2, sc_cancer_caf_args[[ct]]$genes_filter_fun)],
            .SDcols = -c('id', 'patient', 'cell_type')
        ]

        # Compute bins for the genes based on average expression:
        sc_mat <- sc_data[cell_type == 'cancer', set_colnames(t(.SD), id), .SDcols = all_genes_filtered]
        gene_averages <- sort(rowMeans(sc_mat))
        bins <- setNames(
            cut(seq_along(gene_averages), breaks = length(gene_averages) %/% 110, labels = FALSE, include.lowest = TRUE),
            names(gene_averages)
        )

        # Get EMT markers:
        ct_emt_markers <- unique(
            unlist(
                lapply(
                    deconv_data[grepl(paste(paste0('^', sc_metadata[[ct]]$tcga_cancer_types), collapse = '|'), names(deconv_data))],
                    function(deconv) deconv$genes_filtered[deconv$ordering]
                    # function(deconv) head(deconv$genes_filtered[deconv$ordering], length(deconv$genes_filtered)/3)
                )
            )
        )
        ct_emt_markers <- ct_emt_markers[ct_emt_markers %in% names(sc_data)]

        # Filter EMT markers by expression levels:
        ct_emt_markers <- sc_data[
            cell_type == 'cancer',
            ct_emt_markers[apply(.SD, 2, sc_cancer_caf_args[[ct]]$scores_filter_fun)],
            .SDcols = ct_emt_markers
        ]

        # Define control gene sets for distribution of scores:
        comparable_gene_sets <- lapply(ct_emt_markers, function(g) sample(names(bins)[bins == bins[g]], 100))
        comparable_gene_sets <- lapply(1:100, function(i) sapply(comparable_gene_sets, `[`, i))

        # Calculate EMT scores, then filter the EMT markers for correlation with these EMT scores, and recalculate the EMT scores using the
		# filtered list.  I think it makes more sense to filter for correlation with the initial EMT score than with the Z score (below), since
		# the Z score won't be on a comparable scale to the original distribution of expression levels (though this probably doesn't matter
		# much), and the aim of the Z score is to remove variability arising from poor data quality, which might significantly affect the
		# correlation with the Z scores.  EDIT: I don't think there's any difference between the correlations with EMT scores and with Z scores.

        emt_scores <- rowMeans(
            sapply(
                ct_emt_markers,
                function(g) {sc_mat[g, ] - colMeans(sc_mat[sample(names(bins)[bins == bins[g]], 100), ])},
                USE.NAMES = TRUE
            )
        )

        # sc_data[cell_type == 'cancer', sapply(ct_emt_markers, function(g) cor(get(g), emt_scores))] %>% sort %>% plot

        ct_emt_markers <- sc_data[cell_type == 'cancer', ct_emt_markers[sapply(ct_emt_markers, function(g) cor(get(g), emt_scores)) > 0.3]]

        emt_scores <- rowMeans(
            sapply(
                ct_emt_markers,
                function(g) {sc_mat[g, ] - colMeans(sc_mat[sample(names(bins)[bins == bins[g]], 100), ])},
                USE.NAMES = TRUE
            )
        )

        # Calculate distribution of scores for control gene sets:

        comparable_gene_sets_scores <- lapply(
            comparable_gene_sets,
            function(gli) {
                rowMeans(
                    sapply(
                        gli,
                        function(g) {sc_mat[g, ] - colMeans(sc_mat[sample(names(bins)[bins == bins[g]], 100), ])},
                        USE.NAMES = TRUE
                    )
                )
            }
        )

        # Define the EMT score as a Z score:

        z_scores <- sapply(
            names(emt_scores),
            function(cell_id) {
                distrib <- c(emt_scores[cell_id], sapply(comparable_gene_sets_scores, `[`, cell_id))
                (emt_scores[cell_id] - mean(distrib))/sd(distrib)
            },
            USE.NAMES = FALSE
        )

        emt_scores <- sc_data[cell_type == 'cancer', .(id = id, patient = patient, emt_score = z_scores[id])][
            order(patient, emt_score),
            pos_frac := (1:.N)/.N,
            by = patient
        ]

        emt_score_line <- ggplot(
            emt_scores[pos_frac > 0.01 & pos_frac < 0.99],
            aes(pos_frac, emt_score, group = patient, colour = as.character(patient))
        ) +
            scale_colour_manual(values = brewer.pal(12, 'Set3')[c(1, 3:10, 12, 11)][1:length(unique(emt_scores$patient))]) +
            geom_line() +
            theme_test()

        inter_intra_emt <- emt_scores[
            ,
            .(
                # inter_emt = .SD[pos_frac > 0.25 & pos_frac < 0.75, mean(emt_score)],
                # intra_emt = .SD[pos_frac > 0.9 & pos_frac < 0.975, mean(emt_score)] - median(emt_score)
                inter_emt = median(emt_score),
                intra_emt = quantile(emt_score, 0.95) - median(emt_score)
            ),
            by = patient
        ]

        inter_intra_plot <- ggplot(inter_intra_emt, aes(inter_emt, intra_emt, colour = as.character(patient))) +
            scale_colour_manual(values = brewer.pal(12, 'Set3')[c(1, 3:10, 12, 11)][1:length(unique(emt_scores$patient))]) +
            geom_point() +
            theme_test()

        list(
            plots = list(lineplot = emt_score_line, scatterplot = inter_intra_plot),# violin = emt_score_violin),
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

saveRDS(inter_intra_emt_all_cts, '../data_and_figures/inter_intra_emt_all_cts.rds')

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

# emt_score_violin <- ggplot(
	# emt_scores,
	# aes(
		# factor(patient, levels = emt_scores[, .(med = median(emt_score)), by = patient][order(med), patient]),
		# emt_score,
		# fill = as.character(patient)
	# )
# ) +
	# scale_fill_manual(values = brewer.pal(12, 'Set3')[1:length(unique(emt_scores$patient))]) +
	# geom_violin(draw_quantiles = 0.5) +
	# theme_test()





# Plots combining cancer types:

emt_scores_all_cts <- rbindlist(
    lapply(
        names(inter_intra_emt_all_cts),
        function(ct) {
            merge(
                cbind(
                    cancer_type = ct,
                    cancer_type_nice = sc_cancer_caf_heatmaps_args[[ct]]$annotations_title,
                    inter_intra_emt_all_cts[[ct]]$data$cell_emt_scores
                ),
                inter_intra_emt_all_cts[[ct]]$data$inter_intra_emt_scores,
                by = 'patient'
            )
        }
    )
)[
    ,
    patient_unique := paste(cancer_type, patient, sep = '.')
]

# inter_intra_emt_scores_all_cts <- rbindlist(
#     lapply(
#         names(inter_intra_emt_all_cts),
#         function(ct) {
#             cbind(
#                 cancer_type = ct,
#                 inter_intra_emt_all_cts[[ct]]$data$inter_intra_emt_scores
#             )
#         }
#     )
# )

# inter_intra_emt_profiles <- ggplot(
    # emt_scores_all_cts,#[pos_frac > 0.01 & pos_frac < 0.99],
    # aes(
        # pos_frac,
        # emt_score,
        # group = patient,
        # # group = interaction(cancer_type, patient),
        # colour = cancer_type_nice
        # # alpha = intra_emt
        # # size = intra_emt
    # )
# ) +
    # # scale_size_continuous(range = c(0.4, 0.8)) +
    # facet_grid(cols = vars(cancer_type_nice)) +
    # geom_line() +
    # scale_x_continuous(expand = c(0, 0)) +
    # scale_y_continuous(expand = c(0, 0), breaks = c(-5, 0, 5)) +
    # theme(
        # panel.background = element_rect(fill = NA, colour = 'black'),
        # panel.grid.major.y = element_line(colour = 'grey', size = 0.3, linetype = 'dotted'),
        # panel.grid.minor.y = element_blank(),
        # panel.grid.major.x = element_blank(),
        # panel.grid.minor.x = element_blank(),
        # axis.text.x = element_blank(),
        # axis.ticks.x = element_blank(),
        # strip.background = element_rect(colour = 'black'),
        # strip.text = element_text(size = 14),
        # legend.position = 'none'
    # ) +
    # labs(x = 'Cells', y = 'pEMT score')

# inter_intra_emt_scatterplot <- ggplot(
    # # unique(emt_scores_all_cts[, .(cancer_type_nice, inter_emt, intra_emt)]),
    # unique(
        # emt_scores_all_cts[
            # ,
            # .(
                # cancer_type_nice = cancer_type_nice,
                # inter_emt_scaled = scale(inter_emt),
                # intra_emt_scaled = scale(intra_emt)
            # )
        # ]
    # ),
    # aes(inter_emt_scaled, intra_emt_scaled, colour = cancer_type_nice)
# ) +
    # geom_hline(yintercept = 0, size = 0.3, colour = 'grey', linetype = 'dotted') +
    # geom_vline(xintercept = 0, size = 0.3, colour = 'grey', linetype = 'dotted') +
    # geom_point() +
    # theme_test() +
    # labs(
        # x = 'Inter-tumour pEMT heterogeneity score', # These are medians, but then scaled...
        # y = 'Intra-tumour pEMT heterogeneity score',
        # colour = 'Cancer type'
    # )

inter_intra_emt_profiles <- ggplot(
    emt_scores_all_cts,
    aes(
        pos_frac,
        emt_score,
        group = patient,
        # colour = cancer_type_nice
        colour = scale(intra_emt)
    )
) +
    facet_grid(cols = vars(cancer_type_nice)) +
    geom_line() +
    scale_x_continuous(expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0), breaks = c(-5, 0, 5)) +
    scale_colour_gradientn(colours = colorRamps::blue2green(50)[11:45]) +
    theme(
        panel.background = element_rect(fill = NA, colour = 'black'),
        panel.grid.major.y = element_line(colour = 'grey', size = 0.3, linetype = 'dotted'),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        axis.title.x = element_text(hjust = 0),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.background = element_rect(colour = 'black'),
        strip.text = element_text(size = 13),
        # legend.position = 'bottom',
        legend.position = c(1, -0.15),
        legend.direction = 'horizontal',
        legend.justification = 'right',
        legend.key.height = unit(10, 'pt'),
        legend.title = element_text(margin = margin(0, 2.5, 0, 0)),
        plot.margin = unit(c(5.5, 5.5, 40, 5.5), 'pt')
    ) +
    labs(x = 'Cells', y = 'pEMT score', colour = 'Intra-tumour pEMT heterogeneity score\n')

# Choosing colours:
# set.seed(2030)
# random_colours <- randomcoloR::distinctColorPalette(20)
# This is what we get:
# random_colours <- c('#d31fc1', '#7b32c9', '#3aeaaa', '#94e522', '#ed8ba4', '#e295cb', '#1f7fa5', '#f7da96', '#db5f3d', '#f293da',
					# '#d82295', '#ef97ee', '#e512c2', '#ffccd8', '#99f7b5', '#199e8a', '#a3e22d', '#5eedc0', '#e88db0', '#c2f794')
# See what they look like:
# ggplot(data.table(x = letters[1:20])) +
	# geom_tile(aes(x = x, y = 0, fill = x)) +
	# scale_fill_manual(values = random_colours) +
	# scale_x_discrete(labels = random_colours, expand = c(0, 0)) +
	# scale_y_continuous(expand = c(0, 0)) +
	# theme(
		# axis.text.y = element_blank(),
		# axis.title.y = element_blank(),
		# axis.ticks.y = element_blank(),
		# axis.text.x = element_text(angle = 90, hjust = 1),
		# legend.position = 'none'
	# )

ct_colours = c(
	'Breast' = '#CC8D81',
	'Colorectal' = '#DCE144',
    'Head and Neck' = '#7F9E9B',
    'Liver' = '#D4527A',
    'Lung Adeno.' = '#DF984C',
    'Pancreatic' = '#8975DC'
)

inter_intra_emt_profiles_gtable <- ggplot_gtable(ggplot_build(inter_intra_emt_profiles))

stript <- which(grepl('strip-t', inter_intra_emt_profiles_gtable$layout$name))

for (i in stript) {

    inter_intra_emt_profiles_gtable$grobs[[i]]$grobs[[1]]$children[[
        which(grepl('rect', inter_intra_emt_profiles_gtable$grobs[[i]]$grobs[[1]]$childrenOrder))
    ]]$gp$fill <-
        # This sapply() loop makes the ct_colours lighter:
        sapply(ct_colours, function(x) colorRampPalette(c(x, 'white'))(10)[4])[
            inter_intra_emt_profiles_gtable$grobs[[i]]$grobs[[1]]$children[[
                which(grepl('titleGrob', inter_intra_emt_profiles_gtable$grobs[[i]]$grobs[[1]]$childrenOrder))
            ]]$children[[1]]$label
        ]

}

inter_intra_emt_boxplot_data <- melt(
    unique(
        emt_scores_all_cts[
            ,
            .(
                cancer_type_nice = cancer_type_nice,
                `Inter-tumour pEMT heterogeneity` = scale(inter_emt)[, 1],
                `Intra-tumour pEMT heterogeneity` = scale(intra_emt)[, 1]
            )
        ]
    ),
    id.vars = 'cancer_type_nice'
)

inter_intra_emt_boxplot <- ggplot(
    inter_intra_emt_boxplot_data,
    aes(x = cancer_type_nice, y = value)
) +
    # geom_boxplot(varwidth = TRUE, outlier.shape = NA) +
    # stat_boxplot(geom = 'errorbar', width = 0.2) +
    geom_boxplot(outlier.shape = NA, width = 0.5) +
    geom_jitter(
		data = inter_intra_emt_jitterplot_data,
		aes(x = cancer_type_nice, y = value, colour = cancer_type_nice),
        width = 0.2,
        size = 2
	) +
	scale_colour_manual(values = ct_colours) +
    facet_grid(cols = vars(variable)) +
    theme_bw()

# Add the following to see that the standard deviation really is the same for both x and
# y (it's hard to believe from looking at the extreme y values!):

# xlim(c(-4.5, 4.5)) +
#     ylim(c(-4.5, 4.5)) +
#     geom_hline(yintercept = 0) +
#     geom_vline(xintercept = 0)

pdf('../data_and_figures/inter_intra_emt.pdf', width = 10, height = 4)
ggdraw(inter_intra_emt_profiles_gtable)
print(inter_intra_emt_boxplot)
dev.off()

# pdf('../data_and_figures/inter_intra_emt_profiles.pdf', width = 10, height = 4)
# inter_intra_emt_profiles
# dev.off()

# pdf('../data_and_figures/inter_intra_emt_scatterplot.pdf', width = 6, height = 4)
# inter_intra_emt_scatterplot
# dev.off()

# pdf('../data_and_figures/inter_intra_emt_violin.pdf', width = 16, height = 6)

# ggplot(
    # emt_scores_all_cts,
    # aes(
        # factor(
            # patient_unique,
            # levels = emt_scores_all_cts[
                # ,
                # .(med = median(emt_score)),
                # by = patient_unique
            # ][
                # order(med),
                # patient_unique
            # ]
        # ),
        # emt_score,
        # fill = cancer_type
    # )
# ) +
    # geom_violin(draw_quantiles = 0.5) +
    # theme_test() +
    # theme(axis.text.x = element_blank()) +
    # labs(x = 'Patient', y = 'EMT score')

# dev.off()





# To make a "scheme" figure, showing what the inter- and intra-tumour heterogeneity
# scores mean in terms of the profiles:

# In the following, the colour '#5B8BAC' is half way between skyblue3 and skyblue4, and
# can be obtained by the command colorRampPalette(c('skyblue3', 'skyblue4'))(3)[2].

inter_intra_emt_scheme <- ggplot(
    data.table(
        y = c(
            seq(qnorm(0.01, mean = 0, sd = 1), qnorm(0.99, mean = 0, sd = 1), length.out = 500),
            seq(qnorm(0.01, mean = -1.5, sd = 0.6), qnorm(0.99, mean = -1.5, sd = 0.6), length.out = 500)
        ),
        grp = c(rep(2, 500), rep(1, 500))
    )[
        ,
        x := c(
            pnorm(y[1:500], mean = 0, sd = 1),
            pnorm(y[501:1000], mean = -1.5, sd = 0.6)
        )
    ],
    aes(x = x, y = y)
) +
    # geom_line(aes(group = as.character(grp), alpha = grp), col = 'dodgerblue3') +
    geom_line(aes(group = as.character(grp), colour = as.character(grp))) +
    geom_text_repel(
        aes(x, y, label = l, colour = as.character(grp)),
        data = data.table(x = 0.5, y = qnorm(0.5, mean = -1.5, sd = 0.6), l = 'median'),
        # nudge_x = 0.2,
        # nudge_y = -0.1,
        nudge_x = 0.15,
        nudge_y = -0.4,
        segment.colour = 'darkgrey',
        size = 3,
        colour = 'sienna3'
    ) +
    geom_text_repel(
        aes(x, y, label = l),
        data = data.table(x = 0.5, y = qnorm(0.5), l = 'median'),
        nudge_x = 0.05,
        nudge_y = 0.8,
        segment.colour = 'darkgrey',
        size = 3,
        colour = '#5B8BAC'
    ) +
    geom_text_repel(
        aes(x, y, label = l),
        data = data.table(x = 0.95, y = qnorm(0.95), l = '95th percentile'),
        nudge_x = -0.2,
        nudge_y = 0.4,
        segment.colour = 'darkgrey',
        size = 3,
        colour = '#5B8BAC'
    ) +
    # scale_alpha_continuous(range = c(0.5, 1)) +
    scale_colour_manual(values = c('1' = 'sienna3', '2' = '#5B8BAC')) +
    scale_x_continuous(
        # breaks = c(0.5, 0.95),
        # labels = c('0.5' = '50%', '0.95' = '95%'),
        expand = c(0.01, 0)
    ) +
    scale_y_continuous(
        limits = c(qnorm(0.01, mean = -1.5, sd = 0.6) - 0.2, qnorm(0.99, mean = 0, sd = 1) + 0.2),
        expand = c(0, 0)
    ) +
    geom_segment( # Arrow between profiles to indicate inter-tumour heterogeneity
        aes(
            x = 0.5,
            y = qnorm(0.5, mean = -1.5, sd = 0.6) + 0.075,
            xend = 0.5,
            yend = qnorm(0.5, mean = 0, sd = 1) - 0.075
        ),
        arrow = arrow(ends = 'both', length = unit(5, 'pt'))
    ) +
    geom_segment( # Lower horizontal dashed line
        aes(
            x = 0.4025,
            y = qnorm(0.5, mean = 0, sd = 1),
            xend = 0.4975,
            yend = qnorm(0.5, mean = 0, sd = 1)
        ),
        linetype = 'dashed',
        colour = 'grey',
        size = 0.25
    ) +
    geom_segment( # Upper horizontal dashed line
        aes(
            x = 0.4025,
            y = qnorm(0.95, mean = 0, sd = 1),
            xend = 0.9475,
            yend = qnorm(0.95, mean = 0, sd = 1)
        ),
        linetype = 'dashed',
        colour = 'grey',
        size = 0.25
    ) +
    geom_segment( # Arrow between horizontal dashed lines to indicate intra-tumour heterogeneity
        aes(
            x = 0.4,
            y = qnorm(0.5, mean = 0, sd = 1) + 0.05,
            xend = 0.4,
            yend = qnorm(0.95, mean = 0, sd = 1) - 0.05
        ),
        arrow = arrow(ends = 'both', length = unit(5, 'pt'))
    ) +
    # geom_segment( # Vertical dashed line from lower profile down to 0.5 on x axis
    #     aes(
    #         x = 0.5,
    #         y = qnorm(0.01, mean = -1.5, sd = 0.6) - 0.2,
    #         xend = 0.5,
    #         yend = qnorm(0.5, mean = -1.5, sd = 0.6)
    #     ),
    #     linetype = 'dashed',
    #     colour = 'grey',
    #     size = 0.25
    # ) +
    # geom_segment( # Vertical dashed line from upper profile down to 0.95 on x axis
    #     aes(
    #         x = 0.95,
    #         y = qnorm(0.01, mean = -1.5, sd = 0.6) - 0.2,
    #         xend = 0.95,
    #         yend = qnorm(0.95, mean = 0, sd = 1)
    #     ),
    #     linetype = 'dashed',
    #     colour = 'grey',
    #     size = 0.25
    # ) +
    geom_point(
        aes(x, y, colour = as.character(grp)),
        data = data.table(
            x = c(0.5, 0.5, 0.95),
            y = c(qnorm(0.5, mean = -1.5, sd = 0.6), qnorm(0.5), qnorm(0.95)),
            grp = c(1, 2, 2)
        )
    ) +
    annotate(
        geom = 'text',
        x = 0.38,
        y = qnorm(0.5, mean = 0, sd = 1) + (qnorm(0.95, mean = 0, sd = 1) - qnorm(0.5, mean = 0, sd = 1))/2,
        label = 'Intra-tumour\nheterogeneity',
        hjust = 1,
        size = 4
    ) +
    annotate(
        geom = 'text',
        x = 0.52,
        y = qnorm(0.5, mean = -1.5, sd = 0.6) + 1,
        label = 'Inter-tumour\nheterogeneity',
        hjust = 0,
        size = 4
    ) +
    theme_test() +
    theme(
        # axis.title.x = element_text(hjust = 0, vjust = 6),
        # axis.text.y = element_blank(),
        # axis.ticks.y = element_blank(),
        axis.title.x = element_text(vjust = 6),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        legend.position = 'none'
    ) +
    labs(x = 'Cells', y = 'pEMT score')
    # labs(x = 'Sorted cells', y = 'pEMT score')

pdf('../data_and_figures/inter_intra_emt_scheme.pdf', width = 3, height = 4)
inter_intra_emt_scheme
dev.off()





# Combination in one aligned plot:

aligned_plots_1 <- align_plots(
    inter_intra_emt_scheme,
    inter_intra_emt_scatterplot,
    align = 'h'
)

aligned_plots_2 <- align_plots(
    inter_intra_emt_scheme,
    inter_intra_emt_profiles,
    align = 'v',
    axis = 'l'
)

aligned_plots_1[[1]]$widths[1] <- aligned_plots_2[[1]]$widths[1]

pdf('../data_and_figures/inter_intra_emt.pdf', width = 10, height = 8.5)

plot_grid(
    blank_plot(),
    plot_grid(
        aligned_plots_1[[1]],
        blank_plot(),
        aligned_plots_1[[2]],
        nrow = 1,
        ncol = 3,
        rel_widths = c(2.3, 0.7, 4)
    ),
    blank_plot(),
    aligned_plots_2[[2]],
    nrow = 4,
    ncol = 1,
    rel_heights = c(0.1, 1.05, 0.05, 1)
) +
    draw_label('A', x = 0, y = 0.975, hjust = -0.1, vjust = 1.1, size = 16, fontface = 2) +
    draw_label('B', x = 0.425, y = 0.975, hjust = -0.1, vjust = 1.1, size = 16, fontface = 2) +
    draw_label('C', x = 0, y = 0.475, hjust = -0.1, vjust = 1.1, size = 16, fontface = 2)

dev.off()





# Trial with just HNSC:

sc_data <- eval(sc_metadata$hnsc$read_quote)

setkey(sc_data, id)

ct_emt_markers <- with(
    deconv_data$`HNSC - Malignant-Basal`,
    head(genes_filtered[ordering], 50)
)

ct_emt_markers <- sc_data[
    cell_type == 'cancer',
    ct_emt_markers[
        apply(.SD, 2, sc_cancer_caf_args$hnsc$scores_filter_fun)
    ],
    .SDcols = ct_emt_markers
]

all_genes_filtered <- sc_data[
    cell_type == 'cancer',
    names(.SD)[
        apply(.SD, 2, sc_cancer_caf_args$hnsc$genes_filter_fun)
    ],
    .SDcols = -c('id', 'patient', 'cell_type')
]

# Using scrabble::score():

emt_scores <- scrabble::score(
    sc_data[cell_type == 'cancer', set_colnames(t(.SD), id), .SDcols = all_genes_filtered],
    list(ct_emt_markers),
    bin.control = TRUE,
    nbin = length(all_genes_filtered) %/% 110
)

emt_scores <- data.table(
    id = rownames(emt_scores),
    patient = sc_data[rownames(emt_scores), patient],
    emt_score = emt_scores[, 1]
)[
    order(patient, emt_score),
    pos_frac := (1:.N)/.N,
    by = patient
]

# I think it's good to exclude at least the top and bottom 1% to eliminate really extreme
# values.  We could even remove more, e.g. top and bottom 2.5%, so we have a sort of "95%
# confidence interval".

emt_score_plot <- ggplot(
    emt_scores[pos_frac > 0.01 & pos_frac < 0.99],
    aes(pos_frac, emt_score, group = patient, colour = as.character(patient))
) +
    scale_colour_manual(values = brewer.pal(12, 'Set3')[-c(2, 11)]) +
    geom_line() +
    theme_test()

# My own implementation of scrabble::score():

# mat <- sc_data[
#     cell_type == 'cancer',
#     set_colnames(t(.SD), id),
#     .SDcols = all_genes_filtered
# ]
#
# gene_averages <- sort(rowMeans(mat))
#
# bins <- setNames(
#     cut(seq_along(gene_averages), breaks = 20, labels = FALSE, include.lowest = TRUE),
#     names(gene_averages)
# )
#
# emt_scores <- rowMeans(
#     sapply(
#         ct_emt_markers,
#         function(g) {
#             mat[g, ] - colMeans(mat[sample(names(bins)[bins == bins[g]], 100), ])
#         },
#         USE.NAMES = TRUE
#     )
# )
#
# emt_scores <- data.table(
#     id = names(emt_scores),
#     patient = sc_data[names(emt_scores), patient],
#     emt_score = emt_scores
# )[
#     order(patient, emt_score),
#     pos_frac := (1:.N)/.N, # .I doesn't work - it goes from 1 to nrow(emt_scores).
#     by = patient
# ]

# Alternative implementation of scrabble::score() using data.table, which is actually slower:

# gene_info <- sc_data[
#     cell_type == 'cancer',
#     .(
#         gene = all_genes_filtered,
#         average = colMeans(.SD)
#     ),
#     .SDcols = all_genes_filtered
# ][
#     order(average),
#     bin := cut(seq_along(average), breaks = 20, labels = FALSE, include.lowest = TRUE)
# ]
#
# emt_scores <- sc_data[
#     cell_type == 'cancer',
#     .(
#         id = id,
#         patient = patient,
#         emt_score = rowMeans(
#             sapply(
#                 ct_emt_markers,
#                 function(g) {
#                     get(g) - rowMeans(
#                         .SD[
#                             ,
#                             gene_info[bin == gene_info[gene == g, bin], sample(gene, 100)],
#                             with = FALSE
#                         ]
#                     )
#                 }
#             )
#         )
#     )
# ][
#     order(patient, emt_score),
#     pos_frac := (1:.N)/.N, # .I doesn't work - it goes from 1 to nrow(emt_scores).
#     by = patient
# ]

# The following calculates the scores per tumour, but I'm not sure I want to do this for
# the current analysis, because I want to preserve inter-tumour heterogeneity.

# emt_scores <- sc_data[
#     cell_type == 'cancer',
#     .(
#         id = id,
#         emt_score = scrabble::score(
#             set_colnames(t(.SD), id),
#             list(ct_emt_markers),
#             bin.control = TRUE,
#             nbin = length(all_genes_filtered) %/% 110
#         )[, 1]
#     ),
#     by = patient,
#     .SDcols = all_genes_filtered
# ][
#     order(patient, emt_score),
#     pos_frac := (1:.N)/.N,
#     by = patient
# ]

# The following also calculates the scores per tumour, but goes one step further - it also
# filters the gene lists on a per-tumour basis, so you have different gene lists for each
# tumour.  This further reduces inter-tumour heterogeneity in the results.

# emt_scores <- sc_data[
#     cell_type == 'cancer',
#     .(
#         id = id,
#         emt_score = scrabble::score(
#             set_colnames(
#                 t(
#                     .SD[
#                         ,
#                         apply(.SD, 2, sc_cancer_caf_args$hnsc$genes_filter_fun),
#                         with = FALSE
#                     ]
#                 ),
#                 id
#             ),
#             list(
#                 .SD[
#                     ,
#                     ct_emt_markers[
#                         apply(.SD, 2, sc_cancer_caf_args$hnsc$scores_filter_fun)
#                     ],
#                     .SDcols = ct_emt_markers
#                 ]
#             ),
#             bin.control = TRUE,
#             nbin = 20,
#             n = 80
#         )[, 1]
#     ),
#     by = patient,
#     .SDcols = -c('id', 'cell_type')
# ][
#     order(patient, emt_score),
#     pos_frac := (1:.N)/.N,
#     by = patient
# ]

# Try to assess inter- and intra-tumour heterogeneity:

inter_intra_emt <- emt_scores[
    ,
    .(
        inter_emt = .SD[pos_frac > 0.25 & pos_frac < 0.75, mean(emt_score)],
        intra_emt = .SD[pos_frac > 0.9 & pos_frac < 0.975, mean(emt_score)] -
            median(emt_score)
    ),
    by = patient
]

inter_intra_plot <- ggplot(
    inter_intra_emt,
    aes(inter_emt, intra_emt, colour = as.character(patient))
) +
    scale_colour_manual(values = brewer.pal(12, 'Set3')[-c(2, 11)]) +
    geom_point() +
    theme_test()

# In the same figure:

ggarrange(
    emt_score_plot,
    inter_intra_plot,
    nrow = 2,
    ncol = 1
)
