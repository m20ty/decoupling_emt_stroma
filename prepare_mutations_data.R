library(data.table)
library(magrittr)
library(Matrix)

source('general_functions.R')

tcga_cancer_types <- c(
    'ACC',
    'BLCA',
    'BRCA',
    'CESC',
    'CHOL',
    'COAD',
    'ESCA',
    'HNSC',
    'KICH',
    'KIRC',
    'KIRP',
    'LIHC',
    'LUAD',
    'LUSC',
    'MESO',
    'OV',
    'PAAD',
    'PRAD',
    'READ',
    'SKCM',
    'STAD',
    'THCA',
    'UCEC',
    'UVM'
)





# The following uses tryCatch() because I was getting the following error from getGistic():

# "Error while downloading CNV data"

# Based on the function code in here:

# https://rdrr.io/github/BioinformaticsFMRP/TCGAbiolinks/src/R/internal.R

# I think the error comes from checking MD5 checksums.  I don't know how to fix this!

all_gistic_data <- sapply(
    tcga_cancer_types,
    function(ct) {
        tryCatch(
            TCGAbiolinks::getGistic(ct),
            error = function(c) NULL
        )
    },
    simplify = FALSE,
    USE.NAMES = TRUE
)

# After first attempt, a couple of cancer types failed, as shown by the presence of NULLs:

# sapply(all_gistic_data, is.null)

# I downloaded these manually from the GDAC Firehose website.  They're a little hard to find,
# in that there are several CNV-related files and the ones I want are only in one place: you
# have to click on the grey "CopyNumber Analyses" block on the left, then "CopyNumber Gistic 2"
# in the resulting drop menu.  The file you want from the resulting tarball is called
# "all_thresholded.by_genes.txt".

saveRDS(all_gistic_data, '../data_and_figures/gistic_data_list.rds')





# To make a nice data file from the above RDS file:

all_gistic_data <- readRDS('../data_and_figures/gistic_data_list.rds')

sample_meta <- rbindlist(
    lapply(
        tcga_cancer_types,
        function(ct) {
            data.table(
                id = gsub('-', '\\.', names(all_gistic_data[[ct]][, -(1:3)])),
                cancer_type = ct
            )
        }
    )
)

gistic_data <- Reduce(
    merge,
    all_gistic_data
)

gene_meta <- as.data.table(
    setNames(
        gistic_data[, c('Gene Symbol', 'Locus ID', 'Cytoband')],
        c('symbol', 'locus_id', 'cytoband')
    )
)

# The columns of the data frame are characters, for some reason, hence the following:

gistic_data <- as.data.table(gistic_data)[, -c('Locus ID', 'Cytoband')]

setnames(gistic_data, 'Gene Symbol', 'symbol')

names(gistic_data) <- gsub('-', '\\.', names(gistic_data))

gistic_data[
    ,
    names(gistic_data[, -'symbol']) := lapply(.SD, as.numeric),
    .SDcols = -'symbol'
]





# Get gene annotations:

library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(org.Hs.eg.db)

gene_ids_aliases <- as.data.table(
    select(
        org.Hs.eg.db,
        genes(TxDb.Hsapiens.UCSC.hg19.knownGene)$gene_id,
        c('SYMBOL', 'ALIAS'),
        keytype = 'ENTREZID'
    )
)

gistic_data[
    !grepl('\\|', symbol),
    symbol := switch(
        (symbol %in% gene_ids_aliases$SYMBOL) + 1,
        switch(
            (symbol %in% gene_ids_aliases$ALIAS) + 1,
            'no_unique_match',
            gene_ids_aliases[
                ALIAS == symbol,
                switch(
                    (length(unique(SYMBOL)) > 1 || unique(SYMBOL) %in% gistic_data$symbol) + 1,
                    unique(SYMBOL),
                    'no_unique_match'
                )
            ]
        ),
        symbol
    ),
    by = symbol
]

gistic_data <- gistic_data[
    !grepl('\\|', symbol) & symbol != 'no_unique_match'
]

# Note in the following that the resulting entrez_id column is a character, not an integer
# vector.  This seems weird, but is actually useful because it can then be used to subset
# genes(TxDb.Hsapiens.UCSC.hg19.knownGene) in the right order.  The start and end positions
# do form integer columns, however.

setkey(gene_meta, symbol)

gene_meta <- gene_meta[gistic_data$symbol][
    ,
    entrez_id := unique(gene_ids_aliases[!is.na(SYMBOL), .(ENTREZID), keyby = SYMBOL])[
        symbol,
        ENTREZID
    ]
][
    ,
    c('chrom_start', 'chrom_end') := lapply(
        c(start, end),
        function(FUN) FUN(genes(TxDb.Hsapiens.UCSC.hg19.knownGene)[entrez_id])
    )
]

# Save:

fwrite(gistic_data, '../data_and_figures/gistic_data.csv')
fwrite(gene_meta, '../data_and_figures/gene_meta.csv')
fwrite(sample_meta, '../data_and_figures/sample_meta.csv')

# For saving in matrix format:

# gistic_data <- set_rownames(
#     apply(gistic_data[, -(1:3)], 2, as.numeric),
#     gistic_data$`Gene Symbol`
# )

# writeMM(
#     Matrix(gistic_data, sparse = TRUE),
#     '../data_and_figures/gistic_data.mtx'
# )

# gistic_data <- readMM('../data_and_figures/gistic_data.mtx')

# I gave up on the matrix format because the .mtx file is huge, and reading it in and
# doing operations on it is surprisingly slow...  It even loses the row and column
# names for some reason... Maybe I was doing something wrong.  It seems quicker to read
# in as a data.table and then convert to a normal matrix.

# When reading in gistic_data as a data.table, we might want to consider converting it
# into a matrix straight away:

# gistic_data <- set_rownames(
#     as.matrix(gistic_data[, -'symbol']),
#     gistic_data$symbol
# )

# This is because subsetting columns is a bit faster with matrices than with data.table...
# But not by that much, so it's probably not worth it.
