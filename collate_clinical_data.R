library(data.table) # 1.12.8

source('functions.R')

tcga_cancer_types <- c('ACC', 'BLCA', 'BRCA', 'CESC', 'CHOL', 'COAD', 'ESCA', 'HNSC', 'KICH', 'KIRC', 'KIRP', 'LIHC', 'LUAD', 'LUSC', 'MESO', 'OV',
    'PAAD', 'PRAD', 'READ', 'SKCM', 'STAD', 'THCA', 'THYM', 'UCEC', 'UVM')

clinical_datasets <- sapply(

    tcga_cancer_types,

    function(ct) {

        cat(ct, '\b...')

        if(ct %in% dir('../../TCGA_data') & 'All_CDEs.txt' %in% dir(paste0('../../TCGA_data/', ct))) {
            clinical_data <- fread(paste0('../../TCGA_data/', ct, '/All_CDEs.txt'), showProgress = FALSE)
        } else {

            download.file(
                paste0(
                    'http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/',
                    ct,
                    '/20160128/gdac.broadinstitute.org_',
                    ct,
                    '.Clinical_Pick_Tier1.Level_4.2016012800.0.0.tar.gz'
                ),
                destfile = 'tmp.tar.gz',
                quiet = TRUE
            )

            file_names <- untar('tmp.tar.gz', list = TRUE)
            untar('tmp.tar.gz', files = file_names[endsWith(file_names, 'All_CDEs.txt')], exdir = 'tmp')
            file.remove('tmp.tar.gz')

            # Read in the data:
            clinical_data <- fread(paste0('tmp/', file_names[endsWith(file_names, 'All_CDEs.txt')]), showProgress = FALSE)

            # Remove the created directory:
            unlink('tmp', recursive = TRUE)

        }

        clinical_data <- tdt(clinical_data)[order(id)][, c('id', 'cancer_type') := .(toupper(gsub('-', '\\.', id)), ct)]
        setcolorder(clinical_data, c('id', 'cancer_type'))

        cat('Done!\n')

        clinical_data

    },

    simplify = FALSE,
    USE.NAMES = TRUE

)

# Bind them all together using rbindlist() with fill argument:
clinical_data <- rbindlist(clinical_datasets, fill = TRUE)

rm(clinical_datasets)

# Write to file:
fwrite(clinical_data, '../../TCGA_data/tcga_clinical_data.csv')
