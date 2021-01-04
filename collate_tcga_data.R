library(data.table) # 1.12.8
library(stringr) # 1.4.0
library(readxl) # 1.3.1
library(plyr) # 1.8.6
library(stringi) # 1.5.3
library(org.Hs.eg.db) # 3.10.0
library(AnnotationDbi) # 1.48.0

source('general_functions.R')
source('tcga_functions.R')

gene_names <- fread('../../gene_names.csv', showProgress = FALSE)

tcga_cancer_types <- c('ACC', 'BLCA', 'BRCA', 'CESC', 'CHOL', 'COAD', 'ESCA', 'HNSC', 'KICH', 'KIRC', 'KIRP', 'LIHC', 'LUAD', 'LUSC', 'MESO', 'OV',
    'PAAD', 'PRAD', 'READ', 'SKCM', 'STAD', 'THCA', 'THYM', 'UCEC', 'UVM')





datasets_tpm <- sapply(

    tcga_cancer_types,

    function(ct) {

        cat(ct, '\b...')

        if(ct %in% dir('../../TCGA_data')) {
            expression_data <- fread(paste0('../../TCGA_data/', ct, '/', ct, '_illuminahiseq_rnaseqv2_RSEM_genes.txt'), showProgress = FALSE)
        } else {

            download.file(
                paste0(
                    'http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/',
                    ct,
                    '/20160128/gdac.broadinstitute.org_',
                    ct,
                    '.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes__data.Level_3.2016012800.0.0.tar.gz'
                ),
                destfile = 'tmp.tar.gz',
                quiet = TRUE
            )

            file_names <- untar('tmp.tar.gz', list = TRUE)
            untar('tmp.tar.gz', files = file_names[endsWith(file_names, 'data.txt')], exdir = 'tmp')
            file.remove('tmp.tar.gz')

            # I have to go through the following renaming procedure because the file name seems too long for R to handle.
            file.rename(paste0('tmp/', str_split_fixed(file_names[endsWith(file_names, 'data.txt')], '/', 2)[, 1]), 'tmp/tempdir')

            # Move the data file up one directory and delete tempdir:
            file.rename(
                paste0('tmp/tempdir/', str_split_fixed(file_names[endsWith(file_names, 'data.txt')], '/', 2)[, 2]),
                paste0('tmp/', str_split_fixed(file_names[endsWith(file_names, 'data.txt')], '/', 2)[, 2])
            )
            unlink('tmp/tempdir', recursive = TRUE)

            # Read in the data:
            expression_data <- fread(paste0('tmp/', str_split_fixed(file_names[endsWith(file_names, 'data.txt')], '/', 2)[, 2]), showProgress = FALSE)

            # Remove the created directory:
            unlink('tmp', recursive = TRUE)

        }

        names(expression_data)[1] <- 'id'

        # Extract scaled estimates:
        expression_data <- expression_data[, expression_data[1] %in% c('gene_id', 'scaled_estimate'), with = FALSE][-1]

        # Make columns numeric:
        expression_data[, names(expression_data[, -'id']) := lapply(.SD, as.numeric), .SDcols = -'id']

        # The following obtains gene symbols by converting the Entrez IDs.  This mostly agrees with the output of alias2SymbolTable() applied to the
        # symbols already in the id column, but in the few places where it disagrees, it seems more accurate than alias2SymbolTable().  We remove NAs
        # after the conversion.
        expression_data[, id := mapIds(org.Hs.eg.db, str_split_fixed(id, '\\|', 2)[, 2], keytype = 'ENTREZID', column = 'SYMBOL')]
        expression_data <- expression_data[!is.na(id)]

        # Transpose and split IDs for extraction of patient IDs and sample type codes:
        expression_data <- tdt(expression_data)
        expression_data$id <- gsub('-', '\\.', expression_data$id)

        # Make sure it's in alphabetical order:
        expression_data <- setkey(expression_data, id)

        sample_codes <- str_split_fixed(expression_data$id, '\\.', 5)
        rows_to_keep <- grep('^01|^02|^05|^06|^07|^11', sample_codes[, 4])
        expression_data <- expression_data[rows_to_keep]
        sample_codes <- sample_codes[rows_to_keep, ]

        # Convert to (something like) log TPM by multiplying by 1e+06 and taking log:
        expression_data[, names(expression_data[, -'id']) := log2(.SD*1e+06 + 1), .SDcols = -'id']

        # Construct metadata table:
        meta_data <- data.table(
            id = expression_data$id,
            patient_id = apply(sample_codes[, 1:3], 1, paste, collapse = '.'),
            cancer_type = ct,
            sample_type = mapvalues(
                str_split_fixed(sample_codes[, 4], '[A-Z]', 2)[, 1],
                c('01', '02', '05', '06', '07', '11'),
                c('primary', 'recurrent', 'primary_additional', 'metastatic', 'metastatic_additional', 'normal'),
                warn_missing = FALSE
            )
        )

        cat('Done!\n')

        # Output as list:
        list(expression_data = expression_data, meta_data = meta_data)

    },

    simplify = FALSE,

    USE.NAMES = TRUE

)





# Bind expression data and metadata tables together:
expression_data <- rbindlist(lapply(datasets_tpm, `[[`, 'expression_data'))
meta_data <- rbindlist(lapply(datasets_tpm, `[[`, 'meta_data'))

setkey(expression_data, id)
setkey(meta_data, id)

# Remove datasets_tpm to save memory:
rm(datasets_tpm)

# Write to file:
fwrite(expression_data, '../../TCGA_data/tcga_expression_data.csv')
fwrite(meta_data, '../../TCGA_data/tcga_meta_data.csv')





# Add purity data:

# Purity data from tcga_annotations file:

tcga_purity_annotations <- fread('../../tcga_annotations.csv', showProgress = FALSE)[
    ,
    .(individual_id = individual_id, cancer_type = tumor_type, purity = purity)
][!is.na(purity)]

# Fix empty cancer_type entries in tcga_purity_annotations (corresponding IDs are all in meta_data):
setkey(meta_data, patient_id)
tcga_purity_annotations[cancer_type == '', cancer_type := meta_data[sample_type == 'primary'][individual_id, cancer_type]]

# In the following, we're defining sample_type depending on the IDs in meta_data.  Since the IDs given are patient IDs and it doesn't say if they're
# metastasis, recurrent or whatever, we would assume the purity estimates are derived from primary samples.  So, if the patient ID matches that of a
# primary sample in meta_data, we'll assign it the sample type 'primary'.  If not, we look at what matches there are.  If there's a unique one, we'll
# take that.  If there is not a unique one - e.g. 'metastatic' and 'recurrent', then we'll bail out and just return 'cannot match sample type'.  This
# is because I don't want to make any more assumptions about the sample type (I only think it's safe to assume it's primary if a primary sample is
# available).  This is designed to account for e.g. melanoma cases where they only have samples from the metastases.

# Note there's no need to consider 'primary_additional', because if an ID has a 'primary_additional' entry then it also has a 'primary' entry.

# Note some NAs will appear here, if the patient IDs don't appear in meta_data. These will disappear when we merge the purity data with meta_data.

tcga_purity_annotations[
    ,
    sample_type := meta_data[
        individual_id,
        switch(
            ('primary' %in% sample_type) + 1,
            switch((length(unique(sample_type)) > 1) + 1, unique(sample_type), 'cannot match sample type'),
            'primary'
        )
    ],
    by = individual_id
]

tcga_purity_annotations <- tcga_purity_annotations[sample_type != 'cannot match sample type'] # Also gets rid of NAs
setnames(tcga_purity_annotations, 'individual_id', 'patient_id')
setcolorder(tcga_purity_annotations, c('patient_id', 'cancer_type', 'sample_type', 'purity'))

# Purity data from the paper by Taylor et al:

tcga_purity_taylor <- as.data.table(read_xlsx('../../TCGA_data/taylor_1-s2.0-S1535610818301119-mmc2.xlsx', skip = 1))[
    Type %in% tcga_cancer_types,
    .(
        patient_id = apply(str_split_fixed(Sample, '-', 4)[, 1:3], 1, paste, collapse = '.'),
        cancer_type = Type,
        sample_type = mapvalues(
            str_split_fixed(Sample, '-', 4)[, 4],
            c('01', '02', '05', '06'),
            c('primary', 'recurrent', 'primary_additional', 'metastatic'),
            warn_missing = FALSE
        ),
        purity = Purity
    )
][!is.na(purity)]

# Extra dataset for oesophageal and stomach:

esca_stad_purity <- as.data.table(read_xlsx('../../TCGA_data/ESCA/nature20805-s1-1.xlsx', skip = 1))[
    `Absolute extract purity` != 'NA', # Stops NAs being introduced by coercing 'NA' string to numeric
    .(
        individual_id = gsub('-', '\\.', as.character(barcode)),
        cancer_type = `Disease code`,
        purity = as.numeric(`Absolute extract purity`)
    )
][!is.na(purity)]

esca_stad_purity[
    ,
    sample_type := meta_data[
        individual_id,
        switch(
            ('primary' %in% sample_type) + 1,
            switch((length(unique(sample_type)) > 1) + 1, unique(sample_type), 'cannot match sample type'),
            'primary'
        )
    ],
    by = individual_id
]

# All the sample types end up being 'primary'.

setnames(esca_stad_purity, 'individual_id', 'patient_id')
setcolorder(esca_stad_purity, c('patient_id', 'cancer_type', 'sample_type', 'purity'))

# Merge these data tables (I think this automatically removes duplicate entries, since applying unique() doesn't change the dimensions):
tcga_purity_data <- merge(merge(tcga_purity_annotations, tcga_purity_taylor, all = TRUE), esca_stad_purity, all = TRUE)

rm(tcga_purity_annotations)
rm(tcga_purity_taylor)
rm(esca_stad_purity)

# Average purity values for repeated patient_id/sample_type combinations:
tcga_purity_data <- tcga_purity_data[, .(purity = round(mean(purity), 2)), by = .(patient_id, cancer_type, sample_type)]

# Add purity data to meta_data table:
setkey(tcga_purity_data, patient_id, sample_type)
meta_data[, purity := tcga_purity_data[.(meta_data$patient_id, meta_data$sample_type), purity]]
setkey(meta_data, id)

# Save changes:
fwrite(meta_data, '../../TCGA_data/tcga_meta_data.csv')





# Subtypes data:

# I decided to keep the subtypes data separate from the rest of the metadata.  This is because there are some cancer types with multiple subtype
# assignments (COADREAD and PAAD), which necessitates having multiple rows for the same TCGA sample ID, which would differ only by the subtypes.  To
# keep meta_data tidy and having the same number of rows as expression_data (one per sample ID in expression_data), I think it's easier to have the
# subtypes information in a separate table.

# When we want to infer subtypes, we will have to do it for every subtype assignment for a given cancer type.  So I can supply to the infer_subtypes()
# function a subtypes_dt table consisting of IDs and one set of subtype assignments, and expression_data with the all IDs for that cancer type.  Note
# I'll also want to make inferred_subtype equal subtype for those entries where subtype is not NA or ''.

tcga_subtypes_data <- rbindlist(

    list(

        BLCA = as.data.table(read_xlsx('../../TCGA_data/BLCA/1-s2.0-S0092867417310565-mmc1.xlsx', sheet = 'Master table'))[
            ,
            .(
                individual_id = gsub('-', '\\.', `Case ID`),
                # cancer_type = 'BLCA',
                subtype_ref = 'doi:10.1016/j.cell.2017.09.007',
                subtype = `mRNA cluster`
            )
        ],

        BRCA = as.data.table(read_xls('../../TCGA_data/BRCA/Supplementary Tables 1-4.xls', skip = 1))[
            ,
            .(
                individual_id = gsub('-', '\\.', `Complete TCGA ID`),
                # cancer_type = 'BRCA',
                subtype_ref = 'doi:10.1038/nature11412',
                subtype = `PAM50 mRNA`
            )
        ][subtype != 'NA'],

        COADREAD_isella = as.data.table(read_xlsx('../../TCGA_data/COADREAD/isella_2017_ncomms15107-s12.xlsx', skip = 2))[
            `Dataset ID` == 'TCGA',
            .(
                individual_id = gsub('-', '\\.', `Sample ID`),
                # cancer_type = 'COADREAD', # This won't do - need to separate COAD and READ.  Look in data from TCGA paper, or in expression_data.
                subtype_ref = 'doi:10.1038/ncomms15107',
                subtype = paste('CRIS', stri_extract_last_regex(`CRIS Assignment`, '[A-Z]'))
            )
        ],

        ESCA = as.data.table(read_xlsx('../../TCGA_data/ESCA/nature20805-s1-1.xlsx', skip = 1))[
            `Disease code` == 'ESCA' & `Histological Type` != 'NA',
            .(
                individual_id = gsub('-', '\\.', barcode),
                # cancer_type = `Disease code`,
                subtype_ref = 'doi:10.1038/nature20805',
                subtype = `Histological Type`
            )
        ],

        HNSC = as.data.table(read_xlsx('../../TCGA_data/HNSC/7.2.xlsx'))[
            ,
            .(
                individual_id = as.character(Barcode),
                # cancer_type = 'HNSC',
                subtype_ref = 'doi:10.1038/nature14129',
                subtype = as.character(RNA)
            )
        ],

        LIHC = as.data.table(read_xlsx('../../TCGA_data/LIHC/1-s2.0-S0092867417306396-mmc3.xlsx', skip = 1))[
            `iCluster clusters` != 'NA',
            .(
                individual_id = gsub('-', '\\.', short.id),
                # cancer_type = 'LIHC',
                subtype_ref = 'doi:10.1016/j.cell.2017.05.046',
                subtype = gsub(':', ' ', `iCluster clusters`)
            )
        ],

        LUAD = as.data.table(read_xlsx('../../TCGA_data/LUAD/nature13385-s2.xlsx', sheet = 7, skip = 4))[
            ,
            .(
                individual_id = gsub('-', '\\.', `Tumor ID`),
                # cancer_type = 'LUAD',
                subtype_ref = 'doi:10.1038/nature13385',
                subtype = mapvalues(
                    expression_subtype,
                    c('TRU', 'prox.-inflam', 'prox.-prolif.'),
                    c('Terminal respiratory unit', 'Proximal-inflammatory', 'Proximal-proliferative')
                )
            )
        ],

        LUSC = as.data.table(read_xls('../../TCGA_data/LUSC/data.file.S7.5.clinical.and.genomic.data.table.xls', skip = 2))[
            ,
            .(
                individual_id = gsub('LUSC', 'TCGA', gsub('-', '\\.', `Tumor ID`)),
                # cancer_type = 'LUSC',
                subtype_ref = 'doi:10.1038/nature11404',
                subtype = get(names(.SD)[length(.SD)]) #We want the last column.
            )
        ],

        OV = as.data.table(read_xls('../../TCGA_data/OV/JCI65833sd1.xls', skip = 1, sheet = 1))[
            DATASET == 'TCGA-discovery',
            .(
                individual_id = gsub('-', '\\.', ID),
                # cancer_type = 'OV',
                subtype_ref = 'doi:10.1172/JCI65833',
                subtype = SUBTYPE
            )
        ],

        # The actual data table I got the PAAD subtypes from is taken from the TCGA paper: doi:10.1016/j.ccell.2017.07.007
        PAAD = as.data.table(read_xlsx('../../TCGA_data/PAAD/mmc2.xlsx', skip = 1))[
            ,
            .(
                individual_id = rep(
                    sapply(
                        `Tumor Sample ID`,
                        function(tsid) {
                            paste(str_split_fixed(tsid, '-', 4)[, 1:3], collapse = '.')
                        } # The 4th part of the IDs is '01A' for all samples, so can remove this part
                    ),
                    3
                ),
                # cancer_type = 'PAAD',
                subtype_ref = c(
                    rep('doi:10.1038/ng.3398', .N), # Moffitt
                    rep('doi:10.1038/nm.2344', .N), # Collisson
                    rep('doi:10.1038/nature16965', .N) # Bailey
                ),
                subtype = c(
                    mapvalues(
                        `mRNA Moffitt clusters (All 150 Samples) 1basal  2classical`,
                        c(1, 2),
                        c('Basal', 'Classical')
                    ),
                    mapvalues(
                        `mRNA Collisson clusters (All 150 Samples) 1classical 2exocrine 3QM`,
                        c(1, 2, 3),
                        c('Classical', 'Exocrine', 'Quasi-mesenchymal')
                    ),
                    mapvalues(
                        `mRNA Bailey Clusters (All 150 Samples) 1squamous 2immunogenic 3progenitor 4ADEX`,
                        c(1, 2, 3, 4),
                        c('Squamous', 'Immunogenic', 'Progenitor', 'ADEX')
                    )
                )
            )
        ],

        PRAD = as.data.table(read_xlsx('../../TCGA_data/PRAD/1-s2.0-S0092867415013392-mmc2.xlsx'))[
            ,
            .(
                individual_id = gsub('-', '\\.', PATIENT_ID), # 4th parts of IDs are all 01, and all patient IDs are unique
                # cancer_type = 'PRAD',
                subtype_ref = 'doi:10.1016/j.cell.2015.10.025',
                subtype = mRNA_cluster
            )
        ],

        # Note the following Excel file for SKCM also contains purity data, which I haven't used.
        SKCM = as.data.table(read_xlsx('../../TCGA_data/SKCM/1-s2.0-S0092867415006340-mmc2.xlsx', sheet = 'Supplemental Table S1D', skip = 1))[
            `RNASEQ-CLUSTER_CONSENHIER` != '-', # All *patient* IDs in this subset are unique.
            .(
                individual_id = apply(str_split_fixed(Name, '-', 4)[, 1:3], 1, paste, collapse = '.'),
                # cancer_type = 'SKCM',
                subtype_ref = 'doi:10.1016/j.cell.2015.05.044',
                subtype = `RNASEQ-CLUSTER_CONSENHIER`
            )
        ],

        # We're using the same data table for the STAD subtypes that we used for ESCA, but we'll be using the 'Gastric classification' column.  I
        # changed the subtype_ref to make it distinct from the ESCA one, because we don't want to confuse the cancer types when inferring subtypes.
        # It's a bit ugly, though...  I hoped to use paper refs as unique identifiers of subtype classifications, but alas...
        STAD = as.data.table(read_xlsx('../../TCGA_data/ESCA/nature20805-s1-1.xlsx', skip = 1))[
            `Disease code` == 'STAD',
            .(
                individual_id = gsub('-', '\\.', barcode),
                # cancer_type = `Disease code`,
                subtype_ref = 'doi:10.1038/nature20805 - STAD',
                subtype = `Gastric classification`
            )
        ]

    )

)





# Infer subtypes:

# The following shows that I only ever have one sample for each patient which is designated 'primary' - presumably any other samples from that
# patient's primary tumour are called 'primary_additional'.  Similar code shows that each patient has only one sample designated 'metastatic'.
# meta_data[, .(repeats = sum(sample_type == 'primary')), by = patient_id][, sum(repeats > 1)]

# Since I only have patient IDs in the subtypes data, I can thus use the primary samples, or, if these don't exist, the metastatic ones, to infer
# subtypes, without any ambiguity.  If neither primary nor metastatic subtypes exist, I can use the recurrent one (actually there's only one such
# patient in our data).  The only really important thing is to exclude patients for which we only have normal samples.

# Maybe, just to be safe, check how many samples have just metastatic and recurrent samples.  Which do I want to use in this case?  Also check that
# the numbers add up, i.e. tumours that have primary + tumours that have only metastatic + only recurrent + only normal == all.
# EDIT: I checked the former, and there are only 2 patients which don't have primary samples and which have more than one unique sample type.  These
# have, respectively, metastatic and metastatic_additional samples and metastatic and normal samples.  So it's safe to take, in order of priority, the
# primary, metastatic and recurrent.

setkey(meta_data, patient_id)

tcga_subtypes_data_centred <- rbindlist(
    lapply(
        unique(tcga_subtypes_data$subtype_ref),
        function(ref) {

            cat(ref, '\b...')

            subtypes_dt <- data.table(
                id = meta_data[
                    tcga_subtypes_data[ # Condition on individual_id avoids NAs
                        subtype_ref == ref & individual_id %in% meta_data$patient_id,
                        individual_id
                    ],
                    .(
                        id = id[ # Gets IDs of to primary/metastatic/recurrent samples
                            sample_type == switch(
                                ('primary' %in% sample_type) + 1,
                                switch(('metastatic' %in% sample_type) + 1, 'recurrent', 'metastatic'),
                                'primary'
                            )
                        ]
                    ),
                    by = patient_id
                ][, id]
            )

            subtypes_dt$subtype <- tcga_subtypes_data[individual_id %in% meta_data[id %in% subtypes_dt$id, patient_id] & subtype_ref == ref, subtype]

			expression_data_centred <- copy(
				expression_data[
					meta_data[ # Gets all rows with cancer types corresponding to ref
						cancer_type %in% meta_data[
							tcga_subtypes_data[ # The condition on individual_id avoids NAs
								subtype_ref == ref & individual_id %in% meta_data$patient_id,
								individual_id
							],
							unique(cancer_type)
						],
						.(
							id = id[ # Gets IDs of primary/metastatic/recurrent samples
								sample_type == switch(
									('primary' %in% sample_type) + 1,
									switch(('metastatic' %in% sample_type) + 1, 'recurrent', 'metastatic'),
									'primary'
								)
							]
						),
						by = patient_id
					][, id]
				]
			)

			expression_data_centred[, (names(expression_data_centred[, -'id'])) := lapply(.SD, scale), .SDcols = -'id']

            inferred_subtypes <- infer_subtypes(expression_data_centred, subtypes_dt = subtypes_dt)$inferred_subtypes_dt

            # Make sure all the (non-normal) IDs from meta_data (for this cancer type) are in inferred_subtypes:

            setkey(inferred_subtypes, id)

            inferred_subtypes <- inferred_subtypes[
                meta_data[
                    cancer_type %in% meta_data[
                        tcga_subtypes_data[ # The condition on individual_id avoids NAs
                            subtype_ref == ref & individual_id %in% meta_data$patient_id,
                            individual_id
                        ],
                        unique(cancer_type)
                    ] & sample_type != 'normal',
                    id
                ]
            ]

            # Replace NAs in subtype column with '':
            inferred_subtypes[is.na(subtype), subtype := '']

            # Include patient IDs and delete min_diff:
            inferred_subtypes[
                ,
                c('patient_id', 'subtype_ref', 'min_diff') := .(
                    meta_data[
                        cancer_type %in% meta_data[
                            tcga_subtypes_data[ # The condition on individual_id avoids NAs
                                subtype_ref == ref & individual_id %in% meta_data$patient_id,
                                individual_id
                            ],
                            unique(cancer_type)
                        ] & sample_type != 'normal',
                        patient_id
                    ],
                    ref,
                    NULL
                )
            ]

            # Fill in empty entries in the inferred_subtype column with the inferred subtype already given to another sample with the same patient ID:
            inferred_subtypes[
                ,
                inferred_subtype := unique(inferred_subtype[!is.na(inferred_subtype)]),
                by = patient_id
            ]

            # Make incorrectly inferred subtypes equal "true" subtype:
            inferred_subtypes[subtype != '', inferred_subtype := subtype]

            setcolorder(inferred_subtypes, c('id', 'patient_id', 'subtype_ref', 'subtype', 'inferred_subtype'))

            any_warnings <- warnings()
            assign('last.warning', NULL, envir = baseenv())

            cat('Done!\n')

            inferred_subtypes

        }
    )
)

setkey(tcga_subtypes_data_centred, id)

tcga_subtypes_data <- rbindlist(
    lapply(
        unique(tcga_subtypes_data$subtype_ref),
        function(ref) {

            cat(ref, '\b...')

            subtypes_dt <- data.table(
                id = meta_data[
                    tcga_subtypes_data[ # Condition on individual_id avoids NAs
                        subtype_ref == ref & individual_id %in% meta_data$patient_id,
                        individual_id
                    ],
                    .(
                        id = id[ # Gets IDs of to primary/metastatic/recurrent samples
                            sample_type == switch(
                                ('primary' %in% sample_type) + 1,
                                switch(('metastatic' %in% sample_type) + 1, 'recurrent', 'metastatic'),
                                'primary'
                            )
                        ]
                    ),
                    by = patient_id
                ][, id]
            )

            subtypes_dt$subtype <- tcga_subtypes_data[individual_id %in% meta_data[id %in% subtypes_dt$id, patient_id] & subtype_ref == ref, subtype]

            inferred_subtypes <- infer_subtypes(
                expression_data[
                    meta_data[ # Gets all rows with cancer types corresponding to ref
                        cancer_type %in% meta_data[
                            tcga_subtypes_data[ # The condition on individual_id avoids NAs
                                subtype_ref == ref & individual_id %in% meta_data$patient_id,
                                individual_id
                            ],
                            unique(cancer_type)
                        ],
                        .(
                            id = id[ # Gets IDs of primary/metastatic/recurrent samples
                                sample_type == switch(
                                    ('primary' %in% sample_type) + 1,
                                    switch(('metastatic' %in% sample_type) + 1, 'recurrent', 'metastatic'),
                                    'primary'
                                )
                            ]
                        ),
                        by = patient_id
                    ][, id]
                ],
                subtypes_dt = subtypes_dt
            )$inferred_subtypes_dt

            # Make sure all the (non-normal) IDs from meta_data (for this cancer type) are in inferred_subtypes:

            setkey(inferred_subtypes, id)

            inferred_subtypes <- inferred_subtypes[
                meta_data[
                    cancer_type %in% meta_data[
                        tcga_subtypes_data[ # The condition on individual_id avoids NAs
                            subtype_ref == ref & individual_id %in% meta_data$patient_id,
                            individual_id
                        ],
                        unique(cancer_type)
                    ] & sample_type != 'normal',
                    id
                ]
            ]

            # Replace NAs in subtype column with '':
            inferred_subtypes[is.na(subtype), subtype := '']

            # Include patient IDs and delete min_diff:
            inferred_subtypes[
                ,
                c('patient_id', 'subtype_ref', 'min_diff') := .(
                    meta_data[
                        cancer_type %in% meta_data[
                            tcga_subtypes_data[ # The condition on individual_id avoids NAs
                                subtype_ref == ref & individual_id %in% meta_data$patient_id,
                                individual_id
                            ],
                            unique(cancer_type)
                        ] & sample_type != 'normal',
                        patient_id
                    ],
                    ref,
                    NULL
                )
            ]

            # Fill in empty entries in the inferred_subtype column with the inferred subtype already given to another sample with the same patient ID:
            inferred_subtypes[
                ,
                inferred_subtype := unique(inferred_subtype[!is.na(inferred_subtype)]),
                by = patient_id
            ]

            # Make incorrectly inferred subtypes equal "true" subtype:
            inferred_subtypes[subtype != '', inferred_subtype := subtype]

            setcolorder(inferred_subtypes, c('id', 'patient_id', 'subtype_ref', 'subtype', 'inferred_subtype'))

            any_warnings <- warnings()

            assign('last.warning', NULL, envir = baseenv())

            cat('Done!\n')

            inferred_subtypes

        }
    )
)

setkey(tcga_subtypes_data, id)





# Restore original key for meta_data:
setkey(meta_data, id)

# Save tcga_subtypes_data:
fwrite(tcga_subtypes_data, '../../TCGA_data/tcga_subtypes_data.csv')
fwrite(tcga_subtypes_data_centred, '../../TCGA_data/tcga_subtypes_data_centred.csv')
