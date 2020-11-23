library(data.table) # 1.12.8
library(ggplot2) # This is version 3.3.1 on my laptop's WSL setup but 3.3.0 on WEXAC.
library(cowplot) # 1.0.0
library(magrittr) # 1.5
library(plyr) # 1.8.6

source('general_functions.R')
source('tcga_functions.R')

expression_data <- fread('../../TCGA_data/tcga_expression_data.csv', key = 'id')
meta_data <- fread('../../TCGA_data/tcga_meta_data.csv', key = 'id')

expression_data <- expression_data[meta_data[sample_type != 'normal', id]]
meta_data <- meta_data[sample_type != 'normal']

deconv_data <- readRDS('../data_and_figures/deconv_data.rds')
deconv_plots <- readRDS('../data_and_figures/deconv_plots.rds')





ct_to_keep <- c('blca_luminal_papillary', 'blca_basal_squamous', 'brca_luminal_a', 'brca_luminal_b', 'brca_basal_like', 'brca_her2_enriched',
    'coad', 'esca_ac', 'hnsc_mesenchymal_basal', 'hnsc_classical', 'luad_proximal_inflammatory', 'luad_proximal_proliferative',
    'luad_terminal_respiratory_unit', 'lusc_classical', 'lusc_secretory', 'ov_differentiated', 'ov_immunoreactive', 'ov_proliferative', 'paad',
    'read', 'stad_cin', 'stad_msi', 'ucec')
nice_names_for_figure <- c('BLCA - Luminal-Papillary', 'BLCA - Basal-Squamous', 'BRCA - Luminal A', 'BRCA - Luminal B', 'BRCA - Basal-like',
    'BRCA - HER2-enriched', 'COAD', 'ESCA - Adenocarcinoma', 'HNSC - Malignant-Basal', 'HNSC - Classical', 'LUAD - Squamoid', 'LUAD - Magnoid',
    'LUAD - Bronchioid', 'LUSC - Classical', 'LUSC - Secretory', 'OV - Differentiated', 'OV - Immunoreactive', 'OV - Proliferative', 'PAAD', 'READ',
    'STAD - CIN', 'STAD - MSI', 'UCEC')

scores_data <- sapply(
    ct_to_keep,
    function(ct) {

        cat(paste0(ct, ':\n'))

        sample_ids_mat <- expression_data[deconv_data[[ct]]$sample_ids, set_colnames(t(.SD), id), .SDcols = -'id']
        max_n <- floor(length(deconv_data[[ct]]$genes_filtered)/3)

        cat('\tComputing mes score\n')
        mes_score <- signature_score(sample_ids_mat, deconv_data[[ct]]$genes_filtered, nbin = 20, n = 100)

        cat('\tComputing scores for random samples\n')
        gene_samples <- replicate(50, sample(deconv_data[[ct]]$genes_filtered, 20, replace = FALSE), simplify = FALSE)
        gene_sample_scores <- lapply(gene_samples, function(smpl) signature_score(sample_ids_mat, smpl, nbin = 20, n = 100))
        ave_mes_corr <- mean(sapply(gene_sample_scores, function(smpl_scores) cor(smpl_scores, mes_score)))

        cat('\tComputing pEMT and CAF scores: ')
        pemt_caf_scores <- lapply(
            20:max_n,
            function(i) {
                if(i == 20) {
                    cat('1/', max_n - 19, sep = '')
                } else {
                    cat(rep('\b', ceiling(log10(i - 19)) + 1 + ceiling(log10(max_n - 19))), i - 19, '/', max_n - 19, sep = '')
                }
                caf_scores <- with(deconv_data[[ct]], signature_score(sample_ids_mat, tail(genes_filtered[ordering], i), nbin = 20, n = 100))
                pemt_scores <- with(deconv_data[[ct]], signature_score(sample_ids_mat, head(genes_filtered[ordering], i), nbin = 20, n = 100))
                if(i == floor(length(deconv_data[[ct]]$genes_filtered)/3)) cat('\n')
                list(caf_scores = caf_scores, pemt_scores = pemt_scores)
            }
        )

        cat('\tCompiling data and fitting model\n')
        pemt_caf_scores_corr <- data.table(
            ngenes = 20:(length(pemt_caf_scores) + 19),
            score_corr = sapply(pemt_caf_scores, function(x) cor(x$pemt_scores, x$caf_scores))
        )
        pemt_caf_scores_mod <- loess(score_corr ~ ngenes, pemt_caf_scores_corr, span = 0.5, degree = 1)
        pemt_caf_scores_corr[, loess_fit := predict(pemt_caf_scores_mod, .SD)]

        list(
            mes_score = mes_score,
            gene_samples = gene_samples,
            gene_sample_scores = gene_sample_scores,
            ave_mes_corr = ave_mes_corr,
            pemt_caf_scores = pemt_caf_scores,
            pemt_caf_scores_corr = pemt_caf_scores_corr,
            pemt_caf_scores_mod = pemt_caf_scores_mod
        )

    },
    simplify = FALSE,
    USE.NAMES = TRUE
)

# scatterplot_data <- rbindlist(
#     lapply(
#         names(scores_data),
#         function(ct) {
#             data.table(
#                 cancer_type = ct,
#                 ngenes = 20:(length(scores_data[[ct]]$pemt_caf_scores) + 19),
#                 score_corr = sapply(scores_data[[ct]]$pemt_caf_scores, function(x) cor(x$pemt_scores, x$caf_scores))
#             )
#         }
#     )
# )

pdf('../data_and_figures/pemt_caf_corr.pdf')
for(ct in ct_to_keep) {
    print(
        ggplot(data = scores_data[[ct]]$pemt_caf_scores_corr) +
            geom_point(aes(x = ngenes, y = score_corr)) +
            geom_path(aes(x = ngenes, y = loess_fit), colour = 'royalblue') +
            theme_test() +
            labs(x = 'Number of genes in signatures', y = 'Correlation between signature scores', title = nice_names_for_figure[ct_to_keep == ct])
        # ggplot(data = plot_data[cancer_type == ct], aes(x = ngenes, y = score_corr)) +
        #     geom_point() +
        #     theme_test() +
        #     labs(x = 'Number of genes in signatures', y = 'Correlation between signature scores', title = nice_names_for_figure[ct_to_keep == ct])
    )
}
dev.off()

# barplot_data <- rbindlist(
#     lapply(
#         ct_to_keep,
#         function(ct) data.table(
#             cancer_type = ct,
#             ave_mes_corr = scores_data[[ct]]$ave_mes_corr,
#             pemt_caf_corr = scores_data[[ct]]$pemt_caf_scores_corr[ngenes == 20, loess_fit]
#         )
#     )
# ) %>% melt(id.vars = 'cancer_type', measure.vars = c('ave_mes_corr', 'pemt_caf_corr'), variable.name = 'corr_type', value.name = 'corr')
#
# barplot_data[
#     ,
#     cancer_type := mapvalues(
#         factor(
#             cancer_type,
#             levels = .SD[, .(corr_diff = abs(diff(corr))), by = cancer_type][order(-corr_diff), cancer_type]
#         ),
#         ct_to_keep,
#         nice_names_for_figure
#     )
# ]
#
# ggplot(data = barplot_data, aes(x = cancer_type, y = corr, group = corr_type, fill = corr_type, colour = corr_type)) +
#     geom_col(position = 'dodge') +
#     theme_test() +
#     theme(axis.text.x = element_text(angle = 55, hjust = 1)) +
#     labs(x = 'Cancer type', y = 'Correlation', fill = NULL, colour = NULL)

barplot_data <- rbindlist(
    lapply(
        ct_to_keep,
        function(ct) data.table(
            cancer_type = ct,
            Strict = scores_data[[ct]]$pemt_caf_scores_corr[ngenes == 20, loess_fit],
            Lenient = scores_data[[ct]]$pemt_caf_scores_corr[ngenes == max(ngenes), loess_fit]
        )
    )
) %>% melt(id.vars = 'cancer_type', measure.vars = c('Strict', 'Lenient'), variable.name = 'Signature threshold', value.name = 'corr')

barplot_data[
    ,
    c('cancer_type', 'Signature threshold') := .(
        mapvalues(
            factor(
                cancer_type,
                levels = .SD[`Signature threshold` == 'Lenient', cancer_type[order(-corr)]]
                # levels = .SD[, .(corr_diff = abs(diff(corr))), by = cancer_type][order(-corr_diff), cancer_type]
            ),
            ct_to_keep,
            nice_names_for_figure
        ),
        factor(`Signature threshold`, levels = c('Lenient', 'Strict'))
    )
]

summary_barplot <- ggplot(
    data = barplot_data,
    aes(x = cancer_type, y = corr, group = `Signature threshold`, fill = `Signature threshold`, colour = `Signature threshold`)
) +
    geom_hline(yintercept = 0, size = 0.5, colour = 'grey') +
    geom_col(position = 'dodge', width = 0.7) +
    theme_test() +
    theme(axis.text.x = element_text(angle = 55, hjust = 1)) +
    labs(x = NULL, y = 'Correlation', fill = NULL, colour = NULL)

scatterplot_data <- rbindlist(
    lapply(c('blca_basal_squamous', 'lusc_classical', 'read'), function(ct) cbind(cancer_type = ct, scores_data[[ct]]$pemt_caf_scores_corr))
)

scatterplots <- ggplot(data = scatterplot_data) +
    geom_point(aes(x = ngenes, y = score_corr), alpha = 0.5) +
    geom_path(aes(x = ngenes, y = loess_fit), colour = 'royalblue') +
    facet_wrap(~mapvalues(cancer_type, ct_to_keep, nice_names_for_figure, warn_missing = FALSE), nrow = 1, scales = 'free') +
    theme_test() +
    theme(plot.title = element_text(margin = margin(b = 30))) +
    labs(x = 'Number of genes in signatures', y = 'Correlation between signature scores', title = 'Correlation of pEMT and CAF signatures')

pdf('../data_and_figures/pemt_caf_corr_combined.pdf', width = 10, height = 7)
plot_grid(
    plot_grid(
        scatterplots,
        summary_barplot + theme(plot.margin = unit(c(20, 5.5, 5.5, 5.5), 'pt'), legend.position = 'none'),
        nrow = 2,
        ncol = 1,
        axis = 'l',
        align = 'v'
    ),
    plot_grid(blank_plot(), get_legend(summary_barplot + theme(legend.justification = c(0, 0.9))), nrow = 2, ncol = 1),
    nrow = 1,
    ncol = 2,
    rel_widths = c(11, 1)
)
dev.off()





ct <- 'brca_luminal_a'

sample_ids_mat <- expression_data[deconv_data[[ct]]$sample_ids, set_colnames(t(.SD), id), .SDcols = -'id']

temp <- lapply(
    10:floor(length(deconv_data[[ct]]$genes_filtered)/3),
    function(i) {

        if(i == 10) cat('\n')
        cat(rep('\b', ceiling(log10(i))), i, sep = '')

        caf_scores <- with(
            deconv_data[[ct]],
            signature_score(sample_ids_mat, tail(genes_filtered[ordering], i), nbin = 20, n = 100)
        )

        pemt_scores <- with(
            deconv_data[[ct]],
            signature_score(sample_ids_mat, head(genes_filtered[ordering], i), nbin = 20, n = 100)
        )

        # pemt_caf_corr <- cor(
        #     with(
        #         deconv_data[[ct]],
        #         expression_data[sample_ids, head(genes_filtered[ordering], i), with = FALSE]
        #     ),
        #     caf_scores
        # )[, 1]

        # list(caf_scores = caf_scores, pemt_caf_corr = pemt_caf_corr)

        list(caf_scores = caf_scores, pemt_scores = pemt_scores)

    }
)

# plot(1:length(temp), sapply(temp, function(x) mean(x$pemt_caf_corr)))
plot(1:length(temp), sapply(temp, function(x) cor(x$pemt_scores, x$caf_scores)))
