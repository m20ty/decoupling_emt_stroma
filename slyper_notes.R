library(rhdf5)
library(magrittr)
library(DropletUtils)

h5ls('../../single_cell_data/slyper_2020/GSM4186956_NSCLC14_fresh-C4_channel1_raw_gene_bc_matrices_h5.h5')

data1 <- h5read(
    '../../single_cell_data/slyper_2020/GSM4186956_NSCLC14_fresh-C4_channel1_raw_gene_bc_matrices_h5.h5',
    'GRCh38/data'
)

all_data <- H5Fopen(
    '../../single_cell_data/slyper_2020/GSM4186956_NSCLC14_fresh-C4_channel1_raw_gene_bc_matrices_h5.h5'
)

all_data$GRCh38$data %>% class

# This wasn't very successful - don't really know what to do with it.  The following at least gives me an
# object of class SingleCellExperiment:

all_data <- read10xCounts(
    '../../single_cell_data/slyper_2020/GSM4186956_NSCLC14_fresh-C4_channel1_raw_gene_bc_matrices_h5.h5'
)

features <- as.data.table(rowData(all_data))
setnames(features, c('gene_id', 'symbol'))

barcodes <- as.data.table(colData(all_data))[, .(barcode = Barcode)]

# The following produces an error "subscript contains NAs":
expression_matrix <- as(assays(all_data)$counts, 'dgCMatrix')

# We also get lots of warnings asking us to use this:
# h5closeAll()
