# bsub -q tirosh -R "rusage[mem=32000]" \
	# -o ../data_and_figures/deconv_bsub.o \
	# -e ../data_and_figures/deconv_bsub.e \
	# Rscript deconv_bsub.R

library(data.table) # 1.12.8
library(ggplot2) # This is version 3.3.1 on my laptop's WSL setup but 3.3.0 on WEXAC.
library(egg) # 0.4.5
library(limma) # 3.42.2
library(cowplot) # 1.0.0
library(RColorBrewer) # 1.1.2

source('general_functions.R')
source('tcga_functions.R')

cell_type_markers <- fread('../../cell_type_markers.csv')

emt_markers <- fread('../../emt_markers.csv')[, gene := alias2SymbolTable(gene)][source != 'GO', sort(unique(gene))]

heatmap_annotations <- c('ACTA2', 'AXL', 'COL1A1', 'COL1A2', 'COL3A1', 'COL4A1', 'COL4A2', 'COL5A1', 'COL5A2', 'COL5A3', 'COL6A1', 'COL6A2',
						 'COL6A3', 'COL7A1', 'CD44', 'CDH2', 'ECM1', 'ECM2', 'FAP', 'FN1', 'IL6', 'ITGA2', 'ITGA5', 'ITGA6', 'ITGB1', 'ITGB3',
						 'ITGB5', 'ITGB6', 'LAMA1', 'LAMA2', 'LAMA3', 'LAMA5', 'LAMB3', 'LAMC1', 'LAMC2', 'MMP1', 'MMP2', 'MMP3', 'MMP10',
						 'MMP14', 'PDPN', 'SNAI1', 'SNAI2', 'SPARC', 'TGFB1', 'TGFBI', 'THY1', 'TNC', 'TWIST1', 'VCAN', 'VIM', 'ZEB1', 'ZEB2')

initial_genes <- c('SNAI1', 'SNAI2', 'TWIST1', 'VIM', 'ZEB1', 'ZEB2')

expression_data <- fread('../../TCGA_data/tcga_expression_data.csv', key = 'id')
meta_data <- fread('../../TCGA_data/tcga_meta_data.csv', key = 'id')
subtypes_data <- fread('../../TCGA_data/tcga_subtypes_data.csv', key = 'id')
ccle <- fread('../../CCLE_cancer_type_Av.csv', key = 'gene_id')
extra_data <- fread('../data_and_figures/collated_extra_data.csv', key = 'gene')

expression_data <- expression_data[meta_data[sample_type != 'normal', id]]
meta_data <- meta_data[sample_type != 'normal']





deconv_args_per_ct <- list(
    
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
    
    blca_basal_squamous = list(
        tcga_cancer_types = 'BLCA',
        subtypes = 'Basal_squamous',
        ref_for_subtypes = 'doi:10.1016/j.cell.2017.09.007',
        ccle_cancer_type = 'urinary tract',
        seed = 1027,
        plot_title = 'BLCA - Basal-squamous'
    ),
    
    blca_neuronal = list(
        tcga_cancer_types = 'BLCA',
        subtypes = 'Neuronal',
        ref_for_subtypes = 'doi:10.1016/j.cell.2017.09.007',
        ccle_cancer_type = 'urinary tract',
        seed = 3939,
        plot_title = 'BLCA - Neuronal',
        gene_weights_fun = function(x) quantile(x, 0.9)
    ),
    
    brca = list(
        tcga_cancer_types = 'BRCA',
        ccle_cancer_type = 'breast',
        extra_data_source = 'breast_qian',
        seed = 9138,
        plot_title = 'BRCA'
    ),
    
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
    
    cesc = list(
        tcga_cancer_types = 'CESC',
        seed = 4825,
        plot_title = 'CESC',
        cell_type_weights = FALSE
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
    
    coadread = list(
        tcga_cancer_types = c('COAD', 'READ'),
        ccle_cancer_type = 'large intestine',
        extra_data_source = 'crc_lee_smc',
        seed = 6578,
        plot_title = 'COADREAD',
        genes_filter_fun = function(x) 1:200,
        gene_weights_fun = function(x) quantile(x, 0.9)
    ),
    
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
        extra_data_source = 'hnscc_puram',
        seed = 9138,
        plot_title = 'HNSC'
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
    
    kich = list(
        tcga_cancer_types = 'KICH',
        ccle_cancer_type = 'kidney',
        seed = 8796,
        plot_title = 'KICH'
    ),
    
    kirc = list(
        tcga_cancer_types = 'KIRC',
        ccle_cancer_type = 'kidney',
        seed = 247,
        plot_title = 'KIRC',
        gene_weights_fun = function(x) quantile(x, 0.9)
    ),
    
    kirp = list(
        tcga_cancer_types = 'KIRP',
        ccle_cancer_type = 'kidney',
        seed = 3031,
        plot_title = 'KIRP',
        cell_type_weights = FALSE
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
    
    luad = list(
        tcga_cancer_types = 'LUAD',
        ccle_cancer_type = 'lung',
        extra_data_source = 'luad_qian',
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
        extra_data_source = 'luad_qian',
        seed = 9955,
        plot_title = 'LUAD - Proximal-inflammatory'
    ),
    
    luad_proximal_proliferative = list(
        tcga_cancer_types = 'LUAD',
        subtypes = 'Proximal-proliferative',
        ref_for_subtypes = 'doi:10.1038/nature13385',
        ccle_cancer_type = 'lung',
        extra_data_source = 'luad_qian',
        seed = 743,
        plot_title = 'LUAD - Proximal-proliferative',
        genes_filter_fun = function(x) 1:225,
        gene_weights_fun = function(x) quantile(x, 0.75)
    ),
    
    luad_terminal_respiratory_unit = list(
        tcga_cancer_types = 'LUAD',
        subtypes = 'Terminal respiratory unit',
        ref_for_subtypes = 'doi:10.1038/nature13385',
        ccle_cancer_type = 'lung',
        extra_data_source = 'luad_qian',
        seed = 6026,
        plot_title = 'LUAD - Terminal respiratory unit',
        genes_filter_fun = function(x) 1:200,
        gene_weights_fun = mean
    ),
    
    lusc = list(
        tcga_cancer_types = 'LUSC',
        ccle_cancer_type = 'lung',
        extra_data_source = 'lusc_qian',
        seed = 8446,
        plot_title = 'LUSC',
        gene_weights_fun = function(x) quantile(x, 0.9)
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
    
    ov = list(
        tcga_cancer_types = 'OV',
        ccle_cancer_type = 'ovary',
		extra_data_source = 'ovarian_qian',
        seed = 8393,
        plot_title = 'OV'
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
    
    ov_mesenchymal = list(
        tcga_cancer_types = 'OV',
        subtypes = 'Mesenchymal',
        ref_for_subtypes = 'doi:10.1172/JCI65833',
        ccle_cancer_type = 'ovary',
		extra_data_source = 'ovarian_qian',
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
    
    paad_basal_moffitt = list(
        tcga_cancer_types = 'PAAD',
        subtypes = 'Basal',
        ref_for_subtypes = 'doi:10.1038/ng.3398',
        ccle_cancer_type = 'pancreas',
        extra_data_source = 'pdac_peng',
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
        extra_data_source = 'pdac_peng',
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
    
    prad = list(
        tcga_cancer_types = 'PRAD',
        ccle_cancer_type = 'prostate',
        seed = 6204,
        plot_title = 'PRAD',
        genes_filter_fun = function(x) 1:200,
        gene_weights_fun = function(x) quantile(x, 0.9)
    ),
    
    read = list(
        tcga_cancer_types = 'READ',
        ccle_cancer_type = 'large intestine',
        seed = 3309,
        plot_title = 'READ',
        extra_data_source = 'crc_lee_smc'
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
    
    ucec = list(
        tcga_cancer_types = 'UCEC',
        ccle_cancer_type = 'endometrium',
        seed = 475,
        plot_title = 'UCEC',
        cell_type_weights = FALSE
    )
    
)





deconv_data <- sapply(
    
    names(deconv_args_per_ct),
    
    function(ct) {
        
        cat(ct, '\b...')
        
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
        
        cat('Done!\n')
        
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
        
        cat(ct, '\b...')
        
        deconv_ct <- deconv_data[[ct]]
        
        deconv_plot_ct <- do.call(
            deconvolve_emt_caf_plots,
            args = c(
                list(
                    data = deconv_ct,
                    heatmap_axis_title = '',
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
        
        cat('Done!\n')
        
        deconv_plot_ct
        
    },
    
    simplify = FALSE,
    USE.NAMES = TRUE
    
)

saveRDS(deconv_data, '../data_and_figures/deconv_data.rds')
saveRDS(deconv_plots, '../data_and_figures/deconv_plots.rds')
