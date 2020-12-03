# Example of call to scrabble::score()

# scrabble::score(
#     set_colnames(t(.SD[, ..all_genes_filtered]), id),
#     list(epi_markers),
#     bin.control = TRUE,
#     nbin = length(all_genes_filtered) %/% 110
# )[, 1]





# Arguments of Julie's scrabble::score() function:

# mat,
# groups,
# binmat = NULL,
# bins = NULL,
# controls = NULL,
# bin.control = F,
# center = F,
# nbin = 30,
# n = 100,
# replace = F





signature_score <- function(
    mat,
    sig_genes,
    nbin = nrow(mat) %/% 110,
    n = 100,
    replace = FALSE,
    return_control_sets = FALSE
) {
    
    # <mat> should have genes for rows and samples for columns.  Set <return_control_sets>
    # to TRUE if you want to also return the control gene sets used to normalise the
    # expression levels.
    
    gene_averages <- sort(rowMeans(mat))
    
    bins <- setNames(
        cut(
            seq_along(gene_averages),
            breaks = nbin,
            labels = FALSE,
            include.lowest = TRUE
        ),
        names(gene_averages)
    )
    
    if(return_control_sets) {
        
        # Define control gene sets for distribution of scores:
        
        sig_gene_controls <- sapply(
            sig_genes,
            function(g) {
                sample(names(bins)[bins == bins[g]], n, replace = replace)
            },
            simplify = FALSE,
            USE.NAMES = TRUE
        )
        
        comparable_gene_sets <- lapply(
            1:n,
            function(i) sapply(sig_gene_controls, `[`, i)
        )
        
        sig_scores <- rowMeans(
            sapply(
                sig_genes,
                function(g) {
                    mat[g, ] - colMeans(
                        mat[sig_gene_controls[[g]], ]
                    )
                },
                USE.NAMES = TRUE
            )
        )
        
        return(
            list(
                scores = sig_scores,
                controls = sig_gene_controls,
                comparable_gene_sets = comparable_gene_sets # Redundant but possibly helpful
            )
        )
        
    } else {
        
        sig_scores <- rowMeans(
            sapply(
                sig_genes,
                function(g) {
                    mat[g, ] - colMeans(
                        mat[sample(names(bins)[bins == bins[g]], n, replace = replace), ]
                    )
                },
                USE.NAMES = TRUE
            )
        )
        
        return(sig_scores)
        
    }
    
}





# Code I used to manually calculate scores, without scrabble::score():

# sc_mat <- sc_data[
#     cell_type == 'cancer',
#     set_colnames(t(.SD), id),
#     .SDcols = all_genes_filtered
# ]
# 
# gene_averages <- sort(rowMeans(sc_mat))
# 
# bins <- setNames(
#     cut(
#         seq_along(gene_averages),
#         breaks = length(gene_averages) %/% 110,
#         labels = FALSE,
#         include.lowest = TRUE
#     ),
#     names(gene_averages)
# )
# 
# # Define control gene sets for distribution of scores:
# 
# comparable_gene_sets <- lapply(
#     ct_emt_markers,
#     function(g) {
#         sample(names(bins)[bins == bins[g]], 100)
#     }
# )
# 
# comparable_gene_sets <- lapply(
#     1:100,
#     function(i) sapply(comparable_gene_sets, `[`, i)
# )
# 
# # Calculate EMT scores:
# 
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
