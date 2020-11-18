library(data.table)
library(Seurat)
library(magrittr)

lambrechts <- fread('../data_and_figures/lambrechts_nsclc_2018.csv')[
    sample_type != 'normal'
][
    ,
    id := paste(id, 'l', sep = '_')
]

song <- fread('../data_and_figures/song_nsclc_2019.csv')[
    ,
    id := paste(id, 's', sep = '_')
]

meta_data <- rbind(
    lambrechts[, .(id, patient, cell_type, disease)][, study := 'l'],
    song[, .(id, patient, cell_type, disease)][, study := 's']
)

# The metadata has to be a data frame with rownames for the Seurat object:

meta_data <- set_rownames(
    as.data.frame(meta_data[, -'id']),
    meta_data$id
)

# The expression data should be a matrix with cells being columns:

expression_matrix <- rbindlist(
    list(
        lambrechts[, -c('patient', 'cell_type', 'cluster', 'disease', 'sample_type', 'annotation')],
        song[, -c('patient', 'cell_type', 'disease')]
    ),
    fill = TRUE
)

expression_matrix <- set_colnames(
    t(
        as.matrix(
            expression_matrix[, -'id']
        )
    ),
    expression_matrix$id
)

expression_matrix <- plyr::mapvalues(expression_matrix, NA, 0)

rm(lambrechts)
rm(song)





# LUSC:

# lusc <- CreateSeuratObject(
#     expression_matrix[, rownames(meta_data[meta_data$disease == 'LUSC', ])],
#     meta.data = meta_data[meta_data$disease == 'LUSC', ]
# )
# 
# lusc_list <- SplitObject(lusc, split.by = 'study')
# 
# for(i in 1:2) {
#     lusc_list[[i]] <- FindVariableFeatures(
#         lusc_list[[i]],
#         selection.method = "vst",
#         nfeatures = 2000,
#         verbose = FALSE
#     )
# }
# 
# lusc_anchors <- FindIntegrationAnchors(
#     object.list = lusc_list,
#     dims = 1:30
# )
# 
# lusc_integrated <- IntegrateData(
#     anchorset = lusc_anchors,
#     dims = 1:30
# )

# I think the following is the integrated data.  It's a sparse matrix, but it's relatively
# small - It has as many cells as the data supplied to CreateSeuratObject, but only 2000
# genes, I guess because of the nfeatures = 2000 argument.

# lusc_integrated@assays$integrated@data

# It's a sparse matrix, but it's relatively small - It has as many cells as the data
# supplied to CreateSeuratObject, but only 2000 genes.  I assume this is because of the
# nfeatures = 2000 argument, though it still only has 2000 genes even if you don't use
# the FindVariableFeatures() function, I think because the argument anchor.features in
# FindIntegrationAnchors is 2000 by default.  So there's probably not much point using
# FindVariableFeatures - the expression matrix you end up with seems to be the same with
# and without it.  Actually, looking at the messages the FindIntegrationAnchors()
# function displays as it's running, it runs FindVariableFeatures() anyway, if you don't
# call it first yourself.  If you set anchor.features = nrow(lusc_list$l@assays$RNA@data)
# in FindIntegrationAnchors() without having fun FindVariableFeatures(), it still only
# gives you 3466 genes/features in the resulting integrated matrix.

# The following gives a matrix of nearly the same dimensions as the input one: it has
# 17135 genes instead of 19071, but I think this is because there are 17135 non-zero
# genes in the input data - the remaining ones are zero for every cell.

lusc <- CreateSeuratObject(
    expression_matrix[, rownames(meta_data[meta_data$disease == 'LUSC', ])],
    meta.data = meta_data[meta_data$disease == 'LUSC', ]
)

lusc_list <- SplitObject(lusc, split.by = 'study')

for(i in 1:2) {
    lusc_list[[i]] <- FindVariableFeatures(
        lusc_list[[i]],
        selection.method = "vst",
        nfeatures = nrow(lusc_list$l@assays$RNA@data),
        verbose = FALSE
    )
}

lusc_anchors <- FindIntegrationAnchors(
    object.list = lusc_list,
    anchor.features = nrow(lusc_list$l@assays$RNA@data),
    dims = 1:30
)

lusc_integrated <- IntegrateData(
    anchorset = lusc_anchors,
    dims = 1:30
)

# Save as data.table:

lusc_integrated <- as.data.table(
    t(as.matrix(lusc_integrated@assays$integrated@data)),
    keep.rownames = 'id'
)[
    ,
    c('patient', 'cell_type') := meta_data[
        id,
        c('patient', 'cell_type')
    ]
]

setcolorder(lusc_integrated, c('id', 'patient', 'cell_type'))

setkey(lusc_integrated, id)

# Note the cell sums aren't anything like the same any more, but oh well.

lusc_integrated[
    ,
    names(lusc_integrated[, -c('id', 'patient', 'cell_type')]) := round(.SD, 4),
    .SDcols = -c('id', 'patient', 'cell_type')
]

lusc_integrated[
    ,
    c('patient', 'cell_type') := .(
        paste0(stringr::str_extract(id, '[a-z]$'), patient),
        plyr::mapvalues(
            tolower(cell_type),
            c('cd14+ monocyte', 't-cell', 'epithelial', 'm2', 'cd141+ dc'),
            c('myeloid', 't_cell', 'cancer', 'myeloid', 'myeloid')
        )
    )
]

fwrite(lusc_integrated, '../data_and_figures/lambrechts_song_lusc_integrated.csv')





# LUAD:

luad <- CreateSeuratObject(
    expression_matrix[, rownames(meta_data[meta_data$disease == 'LUAD', ])],
    meta.data = meta_data[meta_data$disease == 'LUAD', ]
)

luad_list <- SplitObject(luad, split.by = 'study')

for(i in 1:2) {
    luad_list[[i]] <- FindVariableFeatures(
        luad_list[[i]],
        selection.method = "vst",
        nfeatures = nrow(luad_list$l@assays$RNA@data),
        verbose = FALSE
    )
}

luad_anchors <- FindIntegrationAnchors(
    object.list = luad_list,
    anchor.features = nrow(luad_list$l@assays$RNA@data),
    dims = 1:30
)

luad_integrated <- IntegrateData(
    anchorset = luad_anchors,
    dims = 1:30
)

# Save as data.table:

luad_integrated <- as.data.table(
    t(as.matrix(luad_integrated@assays$integrated@data)),
    keep.rownames = 'id'
)[
    ,
    c('patient', 'cell_type') := meta_data[
        id,
        c('patient', 'cell_type')
    ]
]

setcolorder(luad_integrated, c('id', 'patient', 'cell_type'))

setkey(luad_integrated, id)

luad_integrated[
    ,
    names(luad_integrated[, -c('id', 'patient', 'cell_type')]) := round(.SD, 4),
    .SDcols = -c('id', 'patient', 'cell_type')
]

# In the following, we can't use tolower() on the cell types, because we have 'epithelial' from
# the Lambrechts dataset and 'Epithelial' from the Song dataset, the latter of which I want to
# change to 'cancer'.

luad_integrated[
    ,
    c('patient', 'cell_type') := .(
        paste0(stringr::str_extract(id, '[a-z]$'), patient),
        plyr::mapvalues(
            cell_type,
            c(
                'Epithelial',
                'CD14+ Monocyte',
                'NK',
                'Fibroblast',
                'Th',
                'CD1c+ DC',
                'CD8+ T-cell',
                'M2',
                'CD14- Monocyte',
                'B-cell',
                'T-cell'
            ),
            c(
                'cancer',
                'myeloid',
                'nk',
                'fibroblast',
                't_cell',
                'myeloid',
                't_cell',
                'myeloid',
                'myeloid',
                'b_cell',
                't_cell'
            )
        )
    )
]

fwrite(luad_integrated, '../data_and_figures/lambrechts_song_luad_integrated.csv')





# Clean up:

rm(lusc_integrated)
rm(luad_integrated)
rm(lusc)
rm(luad)
rm(lusc_list)
rm(luad_list)
rm(lusc_anchors)
rm(luad_anchors)
rm(expression_matrix)
rm(meta_data)





# Try the single cell analysis on these new datasets:

library(ggplot2)
library(cowplot)

source('functions.R')

expression_data <- fread('../../TCGA_data/tcga_expression_data.csv', key = 'id')
meta_data <- fread('../../TCGA_data/tcga_meta_data.csv', key = 'id')

expression_data <- expression_data[meta_data[sample_type != 'normal', id]]
meta_data <- meta_data[sample_type != 'normal']

# gene_names <- fread('../../gene_names.csv', showProgress = FALSE)

emt_markers <- fread('../../emt_markers.csv')[
    ,
    gene := symbol_to_hgnc(
        gene,
        gene_names
    ),
    by = gene
]

emt_markers <- emt_markers[
    source != 'GO',
    sort(unique(gene))
]





single_cell_metadata <- list(
    
    lusc = list(
        
        tcga_cancer_types = 'LUSC',
        
        read_quote = quote(
            fread('../data_and_figures/lambrechts_song_lusc_integrated.csv')[
                ,
                cell_type := plyr::mapvalues(
                    cell_type,
                    c('', 'alveolar', 'endothelial'),
                    rep('other', 3)
                )
            ]
        )
        
    ),
    
    luad = list(
        
        tcga_cancer_types = 'LUAD',
        
        read_quote = quote(
            fread('../data_and_figures/lambrechts_song_luad_integrated.csv')[
                ,
                cell_type := plyr::mapvalues(
                    cell_type,
                    c('', 'alveolar', 'epithelial', 'nk'),
                    rep('other', 4)
                )
            ]
        )
        
    )
    
)





genes_list <- sapply(
    names(single_cell_metadata),
    function(ct) {
        
        sc_data <- eval(single_cell_metadata[[ct]]$read_quote)
        
        filter_for_groups(
            sc_data[
                ,
                unique(
                    c(
                        'cell_type',
                        emt_markers[emt_markers %in% names(sc_data)],
                        top_cols_by_fun_cor(
                            expression_data[
                                meta_data[
                                    cancer_type %in% single_cell_metadata[[ct]]$tcga_cancer_types,
                                    id
                                ],
                                -'id'
                            ]
                        )[id %in% names(sc_data), id]
                    )
                ),
                with = FALSE
            ]
        )
        
    },
    simplify = FALSE,
    USE.NAMES = TRUE
)





sc_groups_args_list <- list(
    
    lusc = list(
        seed = 775,
        scores_filter_fun = function(x) {mean(x) > 1.5}
    ),
    
    luad = list(
        seed = 636,
        scores_filter_fun = function(x) {mean(x) > 1.5}
    )
    
)





sc_groups_data <- sapply(
    names(single_cell_metadata),
    function(ct) {
        
        sc_data <- eval(single_cell_metadata[[ct]]$read_quote)
        
        set.seed(sc_groups_args_list[[ct]]$seed)
        
        do.call(
            sc_groups,
            args = c(
                list(
                    genes = genes_list[[ct]],
                    sc_data = sc_data[cell_type %in% c('cancer', 'fibroblast')],
                    patient_var = 'patient'
                ),
                sc_groups_args_list[[ct]][-1]
            )
        )
        
    },
    simplify = FALSE,
    USE.NAMES = TRUE
)





heatmap_args_list <- list(
    
    lusc = list(
        annotations = c(
            'COL1A1',
            'COL1A2',
            'CXCL8',
            'FN1',
            'IGFBP3',
            'LAMC2',
            'SDC1',
            'SNAI1',
            'SNAI2',
            'TNC',
            'TPM4',
            'TWIST1',
            'VIM',
            'ZEB1',
            'ZEB2'
        ),
        annotations_title = 'Lung Squamous'
    ),
    
    luad = list(
        annotations = c(
            'COL1A1',
            'COL1A2',
            'VIM',
            'FN1',
            'LGALS1',
            'ITGB1',
            'AREG',
            'TUBA1A',
            'SNAI1',
            'SNAI2',
            'TWIST1',
            'ZEB1',
            'ZEB2',
            'ITGA2'
        ),
        annotations_side = 'right',
        annotations_title = 'Lung Adenocarcinoma'
    )
    
)

sc_heatmaps <- sapply(
    names(single_cell_metadata),
    function(ct) {
        set.seed(sc_groups_args_list[[ct]]$seed)
        do.call(
            sc_groups_heatmap,
            args = c(
                list(
                    sc_groups_list = sc_groups_data[[ct]],
                    groups = c('cancer', 'fibroblast'),
                    default_figure_widths = list(
                        annotations = 2.4,
                        cancer = 8,
                        fibroblast = 1.6
                    ),
                    figure_spacing = 2.5
                ),
                heatmap_args_list[[ct]]
            )
        )
    },
    simplify = FALSE,
    USE.NAMES = TRUE
)





# This is lacking legends, but oh well:

pdf(
    '../data_and_figures/single_cell_heatmaps_lambrechts_song_integrated.pdf',
    width = 11,
    height = 3.3
)

plot_grid(
    plot_grid(
        cowplot_sc(
            sc_heatmaps$lusc,
            legend_space = 0,
            heights = c(1.5, 20, 4),
            es_x_axis_title_vjust = 1.3,
            es_y_axis_title_angle = 0,
            es_y_axis_title_xpos = 0.8
        ),
        blank_plot(),
        cowplot_sc(
            sc_heatmaps$luad,
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
    ncol = 1,
    nrow = 2,
    rel_heights = c(20, 1)
)

dev.off()

# Not sure about this... It looks like there is a quite strong batch effect (at least in
# LUSC, maybe not in LUAD), so further investigation is needed to see if we can correct
# this during the integration.
