# bsub -q tirosh -R rusage[mem=64000] -o deconv_alt.o -e deconv_alt.e Rscript deconv_alt.R

library(data.table)
library(ggplot2)
library(egg)
library(limma)
library(cowplot)

source('general_functions.R')
source('tcga_functions.R')

# I need cell_type_markers for the gene filtering:
cell_type_markers <- fread('../../cell_type_markers.csv')
emt_markers <- fread('../../emt_markers.csv')[, gene := alias2SymbolTable(gene)][source != 'GO', sort(unique(gene))]

heatmap_annotations <- c('ACTA2', 'AXL', 'COL1A1', 'COL1A2', 'COL3A1', 'COL4A1', 'COL4A2', 'COL5A1', 'COL5A2', 'COL5A3', 'COL6A1', 'COL6A2',
	'COL6A3', 'COL7A1', 'CD44', 'CDH2', 'ECM1', 'ECM2', 'FAP', 'FN1', 'IL6', 'ITGA2', 'ITGA5', 'ITGA6', 'ITGB1', 'ITGB3', 'ITGB5', 'ITGB6',
    'LAMA1', 'LAMA2', 'LAMA3', 'LAMA5', 'LAMB3', 'LAMC1', 'LAMC2', 'MMP1', 'MMP2', 'MMP3', 'MMP10', 'MMP14', 'PDPN', 'SNAI1', 'SNAI2', 'SPARC',
    'TGFB1', 'TGFBI', 'THY1', 'TNC', 'TWIST1', 'VCAN', 'VIM', 'ZEB1', 'ZEB2')

initial_genes <- c('SNAI1', 'SNAI2', 'TWIST1', 'VIM', 'ZEB1', 'ZEB2')





# Read in data:
expression_data <- fread('../../TCGA_data/tcga_expression_data.csv', key = 'id')
meta_data <- fread('../../TCGA_data/tcga_meta_data.csv', key = 'id')
subtypes_data <- fread('../../TCGA_data/tcga_subtypes_data.csv', key = 'id')
ccle <- fread('../../CCLE_cancer_type_Av.csv', key = 'gene_id')
extra_data <- fread('../data_and_figures/collated_extra_data.csv', key = 'gene')

# I won't need the normal tissue samples, so let's take them out now:
expression_data <- expression_data[meta_data[sample_type != 'normal', id]]
meta_data <- meta_data[sample_type != 'normal']





deconv_args_per_ct <- list(
    
    # ACC doesn't seem to change much when you tweak the parameters, unless you
    # set cell_type_weights = FALSE, in which case a lot of genes shift sides...
    # This suggests this one isn't very trustworthy...
    acc = list(
        tcga_cancer_types = 'ACC',
        seed = 7135,
        plot_title = 'ACC'
    ),
    
    blca = list(
        tcga_cancer_types = 'BLCA',
        ccle_cancer_type = 'urinary tract',
        seed = 5499,
        plot_title = 'BLCA',
        genes_filter_fun = function(x) 1:225,
        gene_weights_fun = function(x) quantile(x, 0.9)
    ),
    
    blca_luminal_infiltrated = list(
        tcga_cancer_types = 'BLCA',
        subtypes = 'Luminal_infiltrated',
        ref_for_subtypes = 'doi:10.1016/j.cell.2017.09.007',
        ccle_cancer_type = 'urinary tract',
        seed = 6007,
        plot_title = 'BLCA - Luminal-infiltrated',
        genes_filter_fun = function(x) 1:190,
        gene_weights_fun = function(x) quantile(x, 0.9)
    ),
    
    blca_luminal_papillary = list(
        tcga_cancer_types = 'BLCA',
        subtypes = 'Luminal_papillary',
        ref_for_subtypes = 'doi:10.1016/j.cell.2017.09.007',
        ccle_cancer_type = 'urinary tract',
        seed = 1582,
        plot_title = 'BLCA - Luminal-papillary',
        genes_filter_fun = function(x) 1:220,
        gene_weights_fun = function(x) quantile(x, 0.9)
    ),
    
    # This one looks best without cell type weights, but it turns out there is a strong
    # correlation with endothelial cells.  I can't seem to get rid of this without
    # introducing other problems, so I think it's a lost cause.
    blca_luminal = list(
        tcga_cancer_types = 'BLCA',
        subtypes = 'Luminal',
        ref_for_subtypes = 'doi:10.1016/j.cell.2017.09.007',
        ccle_cancer_type = 'urinary tract',
        seed = 9147,
        plot_title = 'BLCA - Luminal',
        genes_filter_fun = function(x) 1:225,
        cell_type_weights = FALSE
    ),
    
    # This one suffers the same problems as the one above, but I have managed to get
    # something that looks OK.  There is still some correlation with endothelial cells
    # which may not pass QC.
    blca_luminal_luminal_infiltrated = list(
        tcga_cancer_types = 'BLCA',
        subtypes = c('Luminal', 'Luminal_infiltrated'),
        ref_for_subtypes = 'doi:10.1016/j.cell.2017.09.007',
        ccle_cancer_type = 'urinary tract',
        seed = 6867,
        plot_title = 'BLCA - Luminal & Luminal-infiltrated',
        genes_filter_fun = function(x) 1:200,
        cell_type_weights = list(
            B_plasma = 0,
            myocyte = 0,
            macrophage = 0,
            endothelial = 0.5,
            DC = 0,
            mast = 0,
            T = 0,
            B = 0.5
        )
    ),
    
    blca_basal_squamous = list(
        tcga_cancer_types = 'BLCA',
        subtypes = 'Basal_squamous',
        ref_for_subtypes = 'doi:10.1016/j.cell.2017.09.007',
        ccle_cancer_type = 'urinary tract',
        seed = 1027,
        plot_title = 'BLCA - Basal-squamous'
    ),
    
    # This one doesn't work that well.
    blca_neuronal = list(
        tcga_cancer_types = 'BLCA',
        subtypes = 'Neuronal',
        ref_for_subtypes = 'doi:10.1016/j.cell.2017.09.007',
        ccle_cancer_type = 'urinary tract',
        seed = 3939,
        plot_title = 'BLCA - Neuronal',
        gene_weights_fun = function(x) quantile(x, 0.9)
    ),
    
    # This one doesn't work that well - it's better to just take basal-squamous and
    # ignore neuronal.
    blca_basal_squamous_neuronal = list(
        tcga_cancer_types = 'BLCA',
        subtypes = c('Basal_squamous', 'Neuronal'),
        ref_for_subtypes = 'doi:10.1016/j.cell.2017.09.007',
        ccle_cancer_type = 'urinary tract',
        seed = 5278,
        plot_title = 'BLCA - Basal-squamous & Neuronal'
    ),
    
    brca = list(
        tcga_cancer_types = 'BRCA',
        ccle_cancer_type = 'breast',
        extra_data_source = 'chung_breast_cancer_2017',
        seed = 9138,
        plot_title = 'BRCA'
    ),
    
    # This looks quite good but might have problems with B plasma cells.
    brca_luminal_a = list(
        tcga_cancer_types = 'BRCA',
        subtypes = 'Luminal A',
        ref_for_subtypes = 'doi:10.1038/nature11412',
        ccle_cancer_type = 'breast',
        extra_data_source = 'chung_breast_cancer_2017',
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
        extra_data_source = 'chung_breast_cancer_2017',
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
        extra_data_source = 'chung_breast_cancer_2017',
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
        extra_data_source = 'chung_breast_cancer_2017',
        seed = 6994,
        plot_title = 'BRCA - HER2-enriched',
        genes_filter_fun = function(x) 1:180,
        gene_weights_fun = function(x) quantile(x, 0.75)
    ),
    
    brca_normal_like = list(
        tcga_cancer_types = 'BRCA',
        subtypes = 'Normal-like',
        ref_for_subtypes = 'doi:10.1038/nature11412',
        ccle_cancer_type = 'breast',
        extra_data_source = 'chung_breast_cancer_2017',
        seed = 3249,
        plot_title = 'BRCA - Normal-like',
        genes_filter_fun = function(x) 1:200,
        gene_weights_fun = function(x) quantile(x, 0.9)
    ),
    
    brca_karaayvaz = list(
        tcga_cancer_types = 'BRCA',
        ccle_cancer_type = 'breast',
        extra_data_source = 'karaayvaz_tnbc_2018',
        seed = 7198,
        plot_title = 'BRCA'
    ),
    
    brca_luminal_a_karaayvaz = list(
        tcga_cancer_types = 'BRCA',
        subtypes = 'Luminal A',
        ref_for_subtypes = 'doi:10.1038/nature11412',
        ccle_cancer_type = 'breast',
        extra_data_source = 'karaayvaz_tnbc_2018',
        seed = 3,
        plot_title = 'BRCA - Luminal A',
        genes_filter_fun = function(x) 1:190,
        gene_weights_fun = function(x) quantile(x, 0.9)
    ),
    
    brca_luminal_b_karaayvaz = list(
        tcga_cancer_types = 'BRCA',
        subtypes = 'Luminal B',
        ref_for_subtypes = 'doi:10.1038/nature11412',
        ccle_cancer_type = 'breast',
        extra_data_source = 'karaayvaz_tnbc_2018',
        seed = 3361,
        plot_title = 'BRCA - Luminal B',
        genes_filter_fun = function(x) 1:190,
        gene_weights_fun = function(x) quantile(x, 0.9)
    ),
    
    brca_basal_like_karaayvaz = list(
        tcga_cancer_types = 'BRCA',
        subtypes = 'Basal-like',
        ref_for_subtypes = 'doi:10.1038/nature11412',
        ccle_cancer_type = 'breast',
        extra_data_source = 'karaayvaz_tnbc_2018',
        seed = 7873,
        plot_title = 'BRCA - Basal-like',
        gene_weights_fun = function(x) quantile(x, 0.75)
    ),
    
    brca_her2_enriched_karaayvaz = list(
        tcga_cancer_types = 'BRCA',
        subtypes = 'HER2-enriched',
        ref_for_subtypes = 'doi:10.1038/nature11412',
        ccle_cancer_type = 'breast',
        extra_data_source = 'karaayvaz_tnbc_2018',
        seed = 4713,
        plot_title = 'BRCA - HER2-enriched',
        genes_filter_fun = function(x) 1:170,
        gene_weights_fun = function(x) quantile(x, 0.9)
    ),
    
    brca_normal_like_karaayvaz = list(
        tcga_cancer_types = 'BRCA',
        subtypes = 'Normal-like',
        ref_for_subtypes = 'doi:10.1038/nature11412',
        ccle_cancer_type = 'breast',
        extra_data_source = 'karaayvaz_tnbc_2018',
        seed = 9312,
        plot_title = 'BRCA - Normal-like'
    ),
    
    cesc = list(
        tcga_cancer_types = 'CESC',
        seed = 4825,
        plot_title = 'CESC',
        cell_type_weights = FALSE
    ),
    
    # This one doesn't work well - purity and cell lines disagree.  Putting
    # gene_weights_fun = function(x) quantile(x, 0.9) makes the disagreement more obvious.
    chol = list(
        tcga_cancer_types = 'CHOL',
        ccle_cancer_type = 'biliary tract',
        seed = 6688,
        plot_title = 'CHOL'
    ),
    
    # This one's quite sensitive to changes.
    coad = list(
        tcga_cancer_types = 'COAD',
        ccle_cancer_type = 'large intestine',
        extra_data_source = 'li_colorectal_2017',
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
    
    coadread = list(
        tcga_cancer_types = c('COAD', 'READ'),
        ccle_cancer_type = 'large intestine',
        extra_data_source = 'li_colorectal_2017',
        seed = 6578,
        plot_title = 'COADREAD',
        genes_filter_fun = function(x) 1:200,
        gene_weights_fun = function(x) quantile(x, 0.9)
    ),
    
    # coadread_cris_a = list(
    #     tcga_cancer_types = c('COAD', 'READ'),
    #     subtypes = 'CRIS A',
    #     ref_for_subtypes = 'doi:10.1038/ncomms15107',
    #     ccle_cancer_type = 'large intestine',
    #     extra_data_source = 'li_colorectal_2017',
    #     seed = 1441,
    #     plot_title = 'COADREAD - CRIS A',
    # ),
    
    # coadread_cris_b = list(
    #     tcga_cancer_types = c('COAD', 'READ'),
    #     subtypes = 'CRIS B',
    #     ref_for_subtypes = 'doi:10.1038/ncomms15107',
    #     ccle_cancer_type = 'large intestine',
    #     extra_data_source = 'li_colorectal_2017',
    #     seed = 3274,
    #     plot_title = 'COADREAD - CRIS B',
    # ),
    
    # coadread_cris_c = list(
    #     tcga_cancer_types = c('COAD', 'READ'),
    #     subtypes = 'CRIS C',
    #     ref_for_subtypes = 'doi:10.1038/ncomms15107',
    #     ccle_cancer_type = 'large intestine',
    #     extra_data_source = 'li_colorectal_2017',
    #     seed = 481,
    #     plot_title = 'COADREAD - CRIS C',
    # ),
    
    # coadread_cris_d = list(
    #     tcga_cancer_types = c('COAD', 'READ'),
    #     subtypes = 'CRIS D',
    #     ref_for_subtypes = 'doi:10.1038/ncomms15107',
    #     ccle_cancer_type = 'large intestine',
    #     extra_data_source = 'li_colorectal_2017',
    #     seed = 3365,
    #     plot_title = 'COADREAD - CRIS D',
    # ),
    
    # coadread_cris_e = list(
    #     tcga_cancer_types = c('COAD', 'READ'),
    #     subtypes = 'CRIS E',
    #     ref_for_subtypes = 'doi:10.1038/ncomms15107',
    #     ccle_cancer_type = 'large intestine',
    #     extra_data_source = 'li_colorectal_2017',
    #     seed = 3377,
    #     plot_title = 'COADREAD - CRIS E',
    # ),
    
    esca = list(
        tcga_cancer_types = 'ESCA',
        ccle_cancer_type = 'oesophagus',
        seed = 9507,
        plot_title = 'ESCA',
        genes_filter_fun = function(x) 1:260
    ),
    
    esca_ac = list(
        tcga_cancer_types = 'ESCA',
        subtypes = 'AC',
        ref_for_subtypes = 'doi:10.1038/nature20805',
        ccle_cancer_type = 'oesophagus',
        seed = 2022,
        plot_title = 'ESCA - AC',
        genes_filter_fun = function(x) 1:175,
        gene_weights_fun = function(x) quantile(x, 0.9)
    ),
    
    esca_escc = list(
        tcga_cancer_types = 'ESCA',
        subtypes = 'ESCC',
        ref_for_subtypes = 'doi:10.1038/nature20805',
        ccle_cancer_type = 'oesophagus',
        seed = 6358,
        plot_title = 'ESCA - ESCC',
        gene_weights_fun = function(x) quantile(x, 0.9)
    ),
    
    hnsc = list(
        tcga_cancer_types = 'HNSC',
        ccle_cancer_type = 'HNSCC',
        extra_data_source = 'puram_hnscc_2017',
        seed = 9138,
        plot_title = 'HNSC'
    ),
    
    hnsc_mesenchymal_basal = list(
        tcga_cancer_types = 'HNSC',
        subtypes = c('Mesenchymal', 'Basal'),
        ref_for_subtypes = 'doi:10.1038/nature14129',
        ccle_cancer_type = 'HNSCC',
        extra_data_source = 'puram_hnscc_2017',
        seed = 7466,
        plot_title = 'HNSC - Mesenchymal & Basal'
    ),
    
    hnsc_classical = list(
        tcga_cancer_types = 'HNSC',
        subtypes = 'Classical',
        ref_for_subtypes = 'doi:10.1038/nature14129',
        ccle_cancer_type = 'HNSCC',
        extra_data_source = 'puram_hnscc_2017',
        seed = 9998,
        plot_title = 'HNSC - Classical',
        gene_weights_fun = function(x) quantile(x, 0.75)
    ),
    
    hnsc_atypical = list(
        tcga_cancer_types = 'HNSC',
        subtypes = 'Atypical',
        ref_for_subtypes = 'doi:10.1038/nature14129',
        ccle_cancer_type = 'HNSCC',
        extra_data_source = 'puram_hnscc_2017',
        seed = 8172,
        plot_title = 'HNSC - Atypical',
        gene_weights_fun = function(x) quantile(x, 0.9)
    ),
    
    # The following doesn't work well - purity and cell lines disagree.  Tweaking
    # gene_weights_fun or genes_filter_fun doesn't do much, while setting
    # cell_type_weights = FALSE has a big impact but makes the whole cell lines bar
    # uniformly high.
    kich = list(
        tcga_cancer_types = 'KICH',
        ccle_cancer_type = 'kidney',
        seed = 8796,
        plot_title = 'KICH'
    ),
    
    # I can't get this one to work.
    kirc = list(
        tcga_cancer_types = 'KIRC',
        ccle_cancer_type = 'kidney',
        seed = 247,
        plot_title = 'KIRC',
        gene_weights_fun = function(x) quantile(x, 0.9)
        # Alternative:
        # genes_filter_fun = function(x) 1:150,
        # cell_type_weights = list(
        #     B_plasma = 0,
        #     myocyte = 0,
        #     macrophage = 0,
        #     endothelial = 1,
        #     DC = 0,
        #     mast = 0,
        #     T = 0,
        #     B = 0
        # )
    ),
    
    # This one actually works OK.
    kirp = list(
        tcga_cancer_types = 'KIRP',
        ccle_cancer_type = 'kidney',
        seed = 3031,
        plot_title = 'KIRP',
        cell_type_weights = FALSE
    ),
    
    # # This works surprisingly OK with cell_type_weights = FALSE.  May have problems with
    # # T cells and macrophages.
    # lihc = list(
    #     tcga_cancer_types = 'LIHC',
    #     ccle_cancer_type = 'liver',
    #     extra_data_source = 'ma_liver_2019',
    #     seed = 9349,
    #     plot_title = 'LIHC',
    #     genes_filter_fun = function(x) 1:260,
    #     cell_type_weights = FALSE
    #     # Alternative:
    #     # genes_filter_fun = function(x) 1:175,
    #     # gene_weights_fun = function(x) quantile(x, 0.9)
    #     # But the purity and cell lines disagree with this, and I think the cell lines
    #     # are more believable.
    # ),
    
    # # EDIT: the following works OK, but could be improved:
    # 
    # lihc = list(
    #     tcga_cancer_types = 'LIHC',
    #     ccle_cancer_type = 'liver',
    #     extra_data_source = 'ma_liver_2019',
    #     seed = 9349,
    #     plot_title = 'LIHC',
    #     genes_filter_fun = function(x) 1:260,
    #     gene_weights_fun = function(x) quantile(x, 0.9)
    # ),
    
    lihc = list(
        tcga_cancer_types = 'LIHC',
        ccle_cancer_type = 'liver',
        extra_data_source = 'ma_liver_2019',
        seed = 9349,
        plot_title = 'LIHC',
        genes_filter_fun = function(x) 1:260,
        cell_type_weights = FALSE
    ),
    
    # lihc_icluster_1 = list(
    #     tcga_cancer_types = 'LIHC',
    #     subtypes = 'iCluster 1',
    #     ref_for_subtypes = 'doi:10.1016/j.cell.2017.05.046',
    #     ccle_cancer_type = 'liver',
    #     seed = 6399,
    #     plot_title = 'LIHC - iCluster 1'
    # ),
    
    # lihc_icluster_2 = list(
    #     tcga_cancer_types = 'LIHC',
    #     subtypes = 'iCluster 2',
    #     ref_for_subtypes = 'doi:10.1016/j.cell.2017.05.046',
    #     ccle_cancer_type = 'liver',
    #     seed = 9030,
    #     plot_title = 'LIHC - iCluster 2'
    # ),
    
    # lihc_icluster_3 = list(
    #     tcga_cancer_types = 'LIHC',
    #     subtypes = 'iCluster 3',
    #     ref_for_subtypes = 'doi:10.1016/j.cell.2017.05.046',
    #     ccle_cancer_type = 'liver',
    #     seed = 1811,
    #     plot_title = 'LIHC - iCluster 3'
    # ),
    
    # The following doesn't work that well - it looks like the CAF genes might really be
    # in the middle, while the ones on the right might be related to something else,
    # maybe also cancer.  Using genes_filter_fun = function(x) 1:175 gave the "best" (I
    # think) result, but it's sensitive to changes.  Arguably, purity_cor_method = 'resid'
    # is slightly better.
    luad = list(
        tcga_cancer_types = 'LUAD',
        ccle_cancer_type = 'lung',
        extra_data_source = 'lambrechts_nsclc_2018',
        seed = 7999,
        plot_title = 'LUAD',
        genes_filter_fun = function(x) 1:175,
        gene_weights_fun = function(x) quantile(x, 0.9)
    ),
    
    luad_proximal_inflammatory = list(
        tcga_cancer_types = 'LUAD',
        subtypes = 'Proximal-inflammatory',
        ref_for_subtypes = 'doi:10.1038/nature13385',
        ccle_cancer_type = 'lung',
        extra_data_source = 'lambrechts_nsclc_2018',
        seed = 9955,
        plot_title = 'LUAD - Proximal-inflammatory'
    ),
    
    luad_proximal_proliferative = list(
        tcga_cancer_types = 'LUAD',
        subtypes = 'Proximal-proliferative',
        ref_for_subtypes = 'doi:10.1038/nature13385',
        ccle_cancer_type = 'lung',
        extra_data_source = 'lambrechts_nsclc_2018',
        seed = 743,
        plot_title = 'LUAD - Proximal-proliferative',
        genes_filter_fun = function(x) 1:225,
        gene_weights_fun = function(x) quantile(x, 0.75)
    ),
    
    # This one doesn't work.  There are problems with endothelial cells and macrophages
    # that I can't get rid of.
    luad_terminal_respiratory_unit = list(
        tcga_cancer_types = 'LUAD',
        subtypes = 'Terminal respiratory unit',
        ref_for_subtypes = 'doi:10.1038/nature13385',
        ccle_cancer_type = 'lung',
        extra_data_source = 'lambrechts_nsclc_2018',
        seed = 6026,
        plot_title = 'LUAD - Terminal respiratory unit',
        genes_filter_fun = function(x) 1:200,
        gene_weights_fun = mean
    ),
    
    lusc = list(
        tcga_cancer_types = 'LUSC',
        ccle_cancer_type = 'lung',
        extra_data_source = 'lambrechts_nsclc_2018',
        seed = 8446,
        plot_title = 'LUSC',
        gene_weights_fun = function(x) quantile(x, 0.9)
    ),
    
    lusc_basal = list(
        tcga_cancer_types = 'LUSC',
        subtypes = 'basal',
        ref_for_subtypes = 'doi:10.1038/nature11404',
        ccle_cancer_type = 'lung',
        extra_data_source = 'lambrechts_nsclc_2018',
        seed = 2287,
        plot_title = 'LUSC - Basal',
        gene_weights_fun = function(x) quantile(x, 0.9)
    ),
    
    lusc_classical = list(
        tcga_cancer_types = 'LUSC',
        subtypes = 'classical',
        ref_for_subtypes = 'doi:10.1038/nature11404',
        ccle_cancer_type = 'lung',
        extra_data_source = 'lambrechts_nsclc_2018',
        seed = 7106,
        plot_title = 'LUSC - Classical',
        gene_weights_fun = function(x) quantile(x, 0.9)
    ),
    
    lusc_primitive = list(
        tcga_cancer_types = 'LUSC',
        subtypes = 'primitive',
        ref_for_subtypes = 'doi:10.1038/nature11404',
        ccle_cancer_type = 'lung',
        extra_data_source = 'lambrechts_nsclc_2018',
        seed = 303,
        plot_title = 'LUSC - Primitive',
        gene_weights_fun = function(x) quantile(x, 0.9)
    ),
    
    lusc_secretory = list(
        tcga_cancer_types = 'LUSC',
        subtypes = 'secretory',
        ref_for_subtypes = 'doi:10.1038/nature11404',
        ccle_cancer_type = 'lung',
        extra_data_source = 'lambrechts_nsclc_2018',
        seed = 9240,
        plot_title = 'LUSC - Secretory'
    ),
    
    # This has an association between endothelial cells and the "CAF end", but otherwise no
    # problems.  There isn't a really obvious cancer cell signature, though, which I can't
    # improve no matter what I try.  There is no agreement between purity and cell lines.
    meso = list(
        tcga_cancer_types = 'MESO',
        ccle_cancer_type = 'pleura',
        seed = 4752,
        plot_title = 'MESO'
    ),
    
    ov = list(
        tcga_cancer_types = 'OV',
        ccle_cancer_type = 'ovary',
        seed = 8393,
        plot_title = 'OV'
    ),
    
    ov_differentiated = list(
        tcga_cancer_types = 'OV',
        subtypes = 'Differentiated',
        ref_for_subtypes = 'doi:10.1172/JCI65833',
        ccle_cancer_type = 'ovary',
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
        seed = 991,
        plot_title = 'OV - Immunoreactive',
        genes_filter_fun = function(x) 1:225,
        gene_weights_fun = function(x) quantile(x, 0.9)
    ),
    
    ov_mesenchymal = list(
        tcga_cancer_types = 'OV',
        subtypes = 'Mesenchymal',
        ref_for_subtypes = 'doi:10.1172/JCI65833',
        ccle_cancer_type = 'ovary',
        seed = 7549,
        plot_title = 'OV - Mesenchymal',
        genes_filter_fun = function(x) 1:175,
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
        # Alternative (may have problems with myocytes and endothelial cells):
        # genes_filter_fun = function(x) 1:175, # Maybe better with 150
        # gene_weights_fun = function(x) quantile(x, 0.9) # Also try 0.95
    ),
    
    ov_proliferative = list(
        tcga_cancer_types = 'OV',
        subtypes = 'Proliferative',
        ref_for_subtypes = 'doi:10.1172/JCI65833',
        ccle_cancer_type = 'ovary',
        seed = 8583,
        plot_title = 'OV - Proliferative',
        gene_weights_fun = function(x) quantile(x, 0.9)
    ),
    
    # For PAAD, I had such a problem with endothelial cells that I resorted to supplying
    # my own cell type weights that only punish endothelial cells.  It works quite well,
    # except the CAF end now looks more like endothelial cells than CAFs - the CAF genes
    # seem to be right of centre but not at the rightmost end.  So we should be careful
    # when we take the rightmost 20 genes and call them the 'CAF markers'.
    paad = list(
        tcga_cancer_types = 'PAAD',
        ccle_cancer_type = 'pancreas',
        extra_data_source = 'elyada_pdac_2019',
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
    
    paad_basal_moffitt = list(
        tcga_cancer_types = 'PAAD',
        subtypes = 'Basal',
        ref_for_subtypes = 'doi:10.1038/ng.3398',
        ccle_cancer_type = 'pancreas',
        extra_data_source = 'elyada_pdac_2019',
        seed = 7839,
        plot_title = 'PAAD - Basal (Moffitt et al.)',
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
    
    paad_classical_moffitt = list(
        tcga_cancer_types = 'PAAD',
        subtypes = 'Classical',
        ref_for_subtypes = 'doi:10.1038/ng.3398',
        ccle_cancer_type = 'pancreas',
        extra_data_source = 'elyada_pdac_2019',
        seed = 7794,
        plot_title = 'PAAD - Classical (Moffitt et al.)',
        genes_filter_fun = function(x) 1:175,
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
    
    # In contrast to previous paad subtypes, endothelial is not a big problem here, but it
    # looks like the CAF end might reflect macrophages.
    paad_classical_collisson = list(
        tcga_cancer_types = 'PAAD',
        subtypes = 'Classical',
        ref_for_subtypes = 'doi:10.1038/nm.2344',
        ccle_cancer_type = 'pancreas',
        extra_data_source = 'elyada_pdac_2019',
        seed = 9973,
        plot_title = 'PAAD - Classical (Collisson et al.)',
        genes_filter_fun = function(x) 1:200,
        gene_weights_fun = function(x) quantile(x, 0.9)
    ),
    
    paad_exocrine_collisson = list(
        tcga_cancer_types = 'PAAD',
        subtypes = 'Exocrine',
        ref_for_subtypes = 'doi:10.1038/nm.2344',
        ccle_cancer_type = 'pancreas',
        extra_data_source = 'elyada_pdac_2019',
        seed = 7809,
        plot_title = 'PAAD - Exocrine (Collisson et al.)',
        genes_filter_fun = function(x) 1:175,
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
    
    # I don't really trust this one, but the diagnostic plot doesn't show any problems,
    # and the patterns of purity, cell lines and scRNA-seq are OK.
    paad_quasi_mesenchymal_collisson = list(
        tcga_cancer_types = 'PAAD',
        subtypes = 'Quasi-mesenchymal',
        ref_for_subtypes = 'doi:10.1038/nm.2344',
        ccle_cancer_type = 'pancreas',
        extra_data_source = 'elyada_pdac_2019',
        seed = 1341,
        plot_title = 'PAAD - Quasi-mesenchymal (Collisson et al.)',
        genes_filter_fun = function(x) 1:175,
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
    
    # This one gives a nice pattern without any effort, and with no diagnostic problems.
    # It looks slightly nicer with gene_weights_fun = function(x) quantile(x, 0.9).
    paad_squamous_bailey = list(
        tcga_cancer_types = 'PAAD',
        subtypes = 'Squamous',
        ref_for_subtypes = 'doi:10.1038/nature16965',
        ccle_cancer_type = 'pancreas',
        extra_data_source = 'elyada_pdac_2019',
        seed = 1416,
        plot_title = 'PAAD - Squamous (Bailey et al.)',
        gene_weights_fun = function(x) quantile(x, 0.9)
    ),
    
    # I also don't trust this one.  It's like the Collisson quasi-mesenchymal one - there
    # don't seem to be any real problems, but the genes don't have much to do with EMT or
    # CAFs.
    paad_immunogenic_bailey = list(
        tcga_cancer_types = 'PAAD',
        subtypes = 'Immunogenic',
        ref_for_subtypes = 'doi:10.1038/nature16965',
        ccle_cancer_type = 'pancreas',
        extra_data_source = 'elyada_pdac_2019',
        seed = 6199,
        plot_title = 'PAAD - Immunogenic (Bailey et al.)',
        genes_filter_fun = function(x) 1:200
    ),
    
    # This one behaves similarly to Bailey's squamous subtype.  We may face problems with
    # B plasma cells, though, which seem to correlate with the cancer end.
    paad_progenitor_bailey = list(
        tcga_cancer_types = 'PAAD',
        subtypes = 'Progenitor',
        ref_for_subtypes = 'doi:10.1038/nature16965',
        ccle_cancer_type = 'pancreas',
        extra_data_source = 'elyada_pdac_2019',
        seed = 590,
        plot_title = 'PAAD - Progenitor (Bailey et al.)',
        genes_filter_fun = function(x) 1:220,
        gene_weights_fun = function(x) quantile(x, 0.9)
    ),
    
    # In ADEX, it seems like T cells are the biggest problem, and putting all the weight
    # on T cells seems to mostly solve it.  The 'CAF end' isn't very CAF-like, though.
    paad_adex_bailey = list(
        tcga_cancer_types = 'PAAD',
        subtypes = 'ADEX',
        ref_for_subtypes = 'doi:10.1038/nature16965',
        ccle_cancer_type = 'pancreas',
        extra_data_source = 'elyada_pdac_2019',
        seed = 3283,
        plot_title = 'PAAD - ADEX (Bailey et al.)',
        genes_filter_fun = function(x) 1:175,
        cell_type_weights = list(
            B_plasma = 0,
            myocyte = 0,
            macrophage = 0,
            endothelial = 0,
            DC = 0,
            mast = 0,
            T = 1,
            B = 0
        )
    ),
    
    # I don't trust this one, and I think we have problems with B plasma and DC cells.
    # putting all the weight on them doesn't seem to help.
    prad = list(
        tcga_cancer_types = 'PRAD',
        ccle_cancer_type = 'prostate',
        seed = 6204,
        plot_title = 'PRAD',
        genes_filter_fun = function(x) 1:200,
        gene_weights_fun = function(x) quantile(x, 0.9)
    ),
    
    # prad_1 = list(
    #     tcga_cancer_types = 'PRAD',
    #     subtypes = '1',
    #     ref_for_subtypes = 'doi:10.1016/j.cell.2015.10.025',
    #     ccle_cancer_type = 'prostate',
    #     seed = 687,
    #     plot_title = 'PRAD - 1'
    # ),
    
    # prad_2 = list(
    #     tcga_cancer_types = 'PRAD',
    #     subtypes = '2',
    #     ref_for_subtypes = 'doi:10.1016/j.cell.2015.10.025',
    #     ccle_cancer_type = 'prostate',
    #     seed = 2498,
    #     plot_title = 'PRAD - 2'
    # ),
    
    # prad_3 = list(
    #     tcga_cancer_types = 'PRAD',
    #     subtypes = '3',
    #     ref_for_subtypes = 'doi:10.1016/j.cell.2015.10.025',
    #     ccle_cancer_type = 'prostate',
    #     seed = 4716,
    #     plot_title = 'PRAD - 3'
    # ),
    
    # This works surprisingly well without interference.
    read = list(
        tcga_cancer_types = 'READ',
        ccle_cancer_type = 'large intestine',
        seed = 3309,
        plot_title = 'READ',
        extra_data_source = 'li_colorectal_2017'
    ),
    
    skcm = list(
        tcga_cancer_types = 'SKCM',
        ccle_cancer_type = 'skin',
        seed = 2536,
        plot_title = 'SKCM'
    ),
    
    skcm_immune = list(
        tcga_cancer_types = 'SKCM',
        subtypes = 'immune',
        ref_for_subtypes = 'doi:10.1016/j.cell.2015.05.044',
        ccle_cancer_type = 'skin',
        seed = 2137,
        plot_title = 'SKCM - Immune'
    ),
    
    skcm_keratin = list(
        tcga_cancer_types = 'SKCM',
        subtypes = 'keratin',
        ref_for_subtypes = 'doi:10.1016/j.cell.2015.05.044',
        ccle_cancer_type = 'skin',
        seed = 921,
        plot_title = 'SKCM - Keratin'
    ),
    
    skcm_mitf_low = list(
        tcga_cancer_types = 'SKCM',
        subtypes = 'MITF-low',
        ref_for_subtypes = 'doi:10.1016/j.cell.2015.05.044',
        ccle_cancer_type = 'skin',
        seed = 2091,
        plot_title = 'SKCM - MITF-low',
        genes_filter_fun = function(x) 1:175,
        cell_type_weights = FALSE
    ),
    
    # Quite extreme parameters for STAD - picks out cancer genes but CAF end is more
    # endothelial.
    stad = list(
        tcga_cancer_types = 'STAD',
        ccle_cancer_type = 'stomach',
        seed = 7658,
        plot_title = 'STAD',
        genes_filter_fun = function(x) 1:150,
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
    
    stad_cin = list(
        tcga_cancer_types = 'STAD',
        subtypes = 'CIN',
        ref_for_subtypes = 'doi:10.1038/nature20805 - STAD',
        ccle_cancer_type = 'stomach',
        seed = 2773,
        plot_title = 'STAD - CIN',
        genes_filter_fun = function(x) 1:200,
        gene_weights_fun = function(x) quantile(x, 0.9)
    ),
    
    stad_ebv = list(
        tcga_cancer_types = 'STAD',
        subtypes = 'EBV',
        ref_for_subtypes = 'doi:10.1038/nature20805 - STAD',
        ccle_cancer_type = 'stomach',
        seed = 3407,
        plot_title = 'STAD - EBV',
        genes_filter_fun = function(x) 1:170,
        cell_type_weights = list(
            B_plasma = 0,
            myocyte = 0,
            macrophage = 0,
            endothelial = 0,
            DC = 0,
            mast = 0,
            T = 1,
            B = 0
        )
    ),
    
    # This one still has a very strong endothelial signature despite the extreme parameters,
    # but at least the cancer end looks OK.
    stad_gs = list(
        tcga_cancer_types = 'STAD',
        subtypes = 'GS',
        ref_for_subtypes = 'doi:10.1038/nature20805 - STAD',
        ccle_cancer_type = 'stomach',
        seed = 5017,
        plot_title = 'STAD - GS',
        genes_filter_fun = function(x) 1:150,
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
    
    # This one's a bit different from the other stomach ones.  The 'CAF end' correlates
    # with macrophages.
    stad_msi = list(
        tcga_cancer_types = 'STAD',
        subtypes = 'MSI',
        ref_for_subtypes = 'doi:10.1038/nature20805 - STAD',
        ccle_cancer_type = 'stomach',
        seed = 6868,
        plot_title = 'STAD - MSI',
        genes_filter_fun = function(x) 1:150,
        gene_weights_fun = function(x) quantile(x, 0.9)
    ),
    
    # This one doesn't work that well.
    thca = list(
        tcga_cancer_types = 'THCA',
        ccle_cancer_type = 'thyroid',
        seed = 6798,
        plot_title = 'THCA',
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
    
    # Thymoma doesn't work because we don't have the purity data - purity column in meta_data
    # is all NAs.  The tcga_annotations.csv file doesn't have data for Thymoma, while the
    # file from the paper of Taylor et al. has entries for Thymoma but all its purity values
    # are NA.
    
    # thym = list(
    #     tcga_cancer_types = 'THYM',
    #     seed = 1207,
    #     plot_title = 'THYM'
    # ),
    
    # This one has a lot of cancer-relevant genes in the middle, but I can't make it better.
    ucec = list(
        tcga_cancer_types = 'UCEC',
        ccle_cancer_type = 'endometrium',
        seed = 475,
        plot_title = 'UCEC',
        cell_type_weights = FALSE
    ),
    
    # This one has a lot of cancer-relevant genes at the CAF end - not sure what the cancer
    # end really is.
    uvm = list(
        tcga_cancer_types = 'UVM',
        seed = 4158,
        plot_title = 'UVM',
        genes_filter_fun = function(x) 1:200
    )
    
)





deconv_data <- sapply(
    names(deconv_args_per_ct),
    function(ct) {
        cat(paste0(ct, '\n'))
        deconv_ct <- do.call(
            deconvolve_emt_caf_data,
            args = c(
                list(
                    expression_data = expression_data,
                    meta_data = meta_data,
                    genes = emt_markers,
                    cell_type_markers = cell_type_markers,
                    # heatmap_annotations = heatmap_annotations,
                    ccle_data = ccle,
                    subtypes_data = subtypes_data,
                    extra_data = extra_data,
                    initial_gene_weights = FALSE
                ),
                deconv_args_per_ct[[ct]][!(names(deconv_args_per_ct[[ct]]) == 'plot_title')]
            )
        )
        deconv_ct
    },
    simplify = FALSE,
    USE.NAMES = TRUE
)

# Reorder genes in deconv results:
deconv_data <- sapply(deconv_data, deconv_reorder, simplify = FALSE, USE.NAMES = TRUE)

deconv_plots <- sapply(
    names(deconv_args_per_ct),
    function(ct) {
        cat(paste0(ct, '\n'))
        deconv_ct <- deconv_data[[ct]]
        deconv_plot_ct <- do.call(
            deconvolve_emt_caf_plots,
            args = c(
                list(
                    data = deconv_ct,
                    # Include the following only if you want epithelial scores (takes much longer):
                    # expression_data = expression_data,
                    heatmap_axis_title = '', # Change heat_map() function so I can put NULL here
                    heatmap_legend_title = 'Coexpression',
                    heatmap_colours = rev(colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))(50)),
                    heatmap_annotations = heatmap_annotations,
                    heatmap_annotations_nudge = 0.3,
                    purity_colours = rev(colorRampPalette(RColorBrewer::brewer.pal(11, "PuOr"))(50)),
                    purity_colour_limits = c(-0.3, 0.3),
                    purity_legend_breaks = c(-0.3, 0, 0.3),
                    purity_legend_title = 'Correlation\nwith purity',
                    purity_legend_direction = 'vertical',
                    purity_axis_title = NULL,
                    ccle_colours = rev(colorRampPalette(RColorBrewer::brewer.pal(11, 'PiYG')[c(1:3, 6, 9:11)])(50)),
					ccle_fun = function(x) {caTools::runmean(x - mean(x), 30)/max(abs(caTools::runmean(x - mean(x), 30)))},
                    ccle_legend_title = 'Tumours vs.\ncell lines',
                    ccle_legend_direction = 'vertical',
                    ccle_axis_title = NULL,
                    extra_colours = colorRampPalette(c('turquoise4', 'turquoise', 'azure', 'gold2', 'gold4'))(50),
                    extra_axis_title = NULL,
                    extra_legend_title = 'scRNA-seq',
                    extra_legend_direction = 'vertical',
                    bar_legend_width = NULL,
                    bar_legend_height = NULL
                ),
                deconv_args_per_ct[[ct]]['plot_title']
            )
        )
        deconv_plot_ct
    },
    simplify = FALSE,
    USE.NAMES = TRUE
)

saveRDS(deconv_data, '../data_and_figures/deconv_data.rds')
saveRDS(deconv_plots, '../data_and_figures/deconv_plots.rds')

# deconv_data <- readRDS('../data_and_figures/deconv_data.rds')
# deconv_plots <- readRDS('../data_and_figures/deconv_plots.rds')





# Write all figures to a single PDF:
pdf('../data_and_figures/deconv_figures.pdf', width = 10, height = 12)
for(deconv_plot_ct in deconv_plots) {print(deconv_plot(list(deconv_plot_ct), legends_space = 0.2))}
dev.off()

# Diagnostic figures:
pdf('../data_and_figures/deconv_diagnostics.pdf', width = 10, height = 12)
i <- 1
for(deconv_plot_ct in deconv_plots) {
    ggarrange(
        plots = c(list(deconv_plot_ct$plots$heatmap), deconv_plot_ct$diagnostics$alternative_purity_cor, deconv_plot_ct$diagnostics$cell_type_bars),
        ncol = 1,
        nrow = 11,
        heights = c(8, rep(0.6, 10)),
        newpage = switch((i == 1) + 1, TRUE, FALSE)
    )
    i <- i + 1
}

dev.off()
