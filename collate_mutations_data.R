library(data.table)

source('functions.R')





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
    'THYM',
    'UCEC',
    'UVM'
)





datasets_mutations <- sapply(
    
    tcga_cancer_types,
    
    function(ct) {
        
        cat(ct, '\b...')
        
        if(
            ct %in% dir('../../TCGA_data') &
            'mutation_packager_raw_calls' %in% dir(
                paste0(
                    '../../TCGA_data/',
                    ct
                )
            )
        ) {
            
            mutations_data_ct <- rbindlist(
                lapply(
                    dir(
                        paste0(
                            '../../TCGA_data/',
                            ct,
                            '/mutation_packager_raw_calls'
                        )
                    ),
                    function(x) {
                        fread(
                            paste0(
                                'TCGA_data/',
                                ct,
                                '/mutation_packager_raw_calls/',
                                x
                            )
                        )
                    }
                )
            )
            
        } else {
            
            download.file(
                paste0(
                    'http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/',
                    ct,
                    '/20160128/gdac.broadinstitute.org_',
                    ct,
                    '.Mutation_Packager_Raw_Calls.Level_3.2016012800.0.0.tar.gz'
                ),
                destfile = 'tmp.tar.gz',
                quiet = TRUE
            )
            
            file_names <- untar('tmp.tar.gz', list = TRUE)
            
            untar(
                'tmp.tar.gz',
                files = file_names[endsWith(file_names, 'maf.txt')],
                exdir = 'tmp'
            )
            
            file.remove('tmp.tar.gz')
            
            mutations_data_ct <- rbindlist(
                lapply(
                    dir(
                        paste0(
                            'tmp/gdac.broadinstitute.org_',
                            ct,
                            '.Mutation_Packager_Raw_Calls.Level_3.2016012800.0.0'
                        )
                    ),
                    function(x) {
                        fread(
                            paste0(
                                'tmp/gdac.broadinstitute.org_',
                                ct,
                                '.Mutation_Packager_Raw_Calls.Level_3.2016012800.0.0/',
                                x
                            )
                        )
                    }
                )
            )
            
            unlink('tmp', recursive = TRUE)
            
        }
        
        mutations_data_ct <- mutations_data_ct[
            ,
            .(
                id = gsub('-', '\\.', Tumor_Sample_Barcode),
                gene = Hugo_Symbol,
                chromosome = Chromosome,
                start_position = Start_position,
                end_position = End_position,
                strand = Strand,
                variant_classification = Variant_Classification,
                variant_type = Variant_Type,
                reference_allele = Reference_Allele,
                tumour_seq_allele_1 = Tumor_Seq_Allele1,
                tumour_seq_allele_2 = Tumor_Seq_Allele2,
                mutation_status = Mutation_Status,
                sequence_source = Sequence_Source,
                genome_change = Genome_Change,
                annotation_transcript = Annotation_Transcript,
                transcript_strand = Transcript_Strand,
                transcript_exon = Transcript_Exon,
                transcript_position = Transcript_Position,
                cdna_change = cDNA_Change,
                codon_change = Codon_Change,
                protein_change = Protein_Change,
                cosmic_total_alterations_in_gene = COSMIC_total_alterations_in_gene,
                ref_context,
                gc_content,
                t_alt_count,
                t_ref_count
            )
        ]
        
        # The problem with this is that the sample IDs from the mutations table don't match
        # the ones in expression_data.  I guess they used a different sample for the whole
        # exome sequencing, so there will be no matches.  We'll have to match based on only
        # the first 3 or 4 parts of the ID.
        
        # Maybe also add columns for EMT and CAF scores for each tumour.  We can use these
        # in statistical tests.
        
        # Note that previously I added another column consisting of 'y' or 'n', denoting
        # whether this tumour has a mutation in this gene.  This, of course, requires an
        # entry for every gene for each tumour, not just the genes which are mutated in
        # that tumour.  I did this to help with the wilcox.test() function, but I'm not
        # sure we need it - if we subset a gene-tumour combination that doesn't exist, we
        # just get NA, so I can test for the difference between tumours which have the
        # mutation vs those that return NA (for each gene).
        
        cat('Done!\n')
        
    }
    
)
