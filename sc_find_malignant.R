library(data.table) # 1.12.8
library(ggplot2) # This is version 3.3.1 on my laptop's WSL setup but 3.3.0 on WEXAC.
library(purrr) # 0.3.4
library(magrittr) # 1.5
library(infercna) # 1.0.0
library(cowplot) # 1.0.0
library(stringr) # 1.4.0

source('general_functions.R')

# Function from infercna:
.chromosomeBreaks = function(genes = NULL, halfway = F, hide = NULL) {
    n = length(genes)
    chrsum = cumsum(lengths(splitGenes(genes, by = 'chr')))
    Breaks = chrsum/max(chrsum) * n
    if (halfway) {
        b = Breaks
        Breaks = c(b[1]/2, b[2:length(b)] - diff(b)/2)
    }
    if (!is.null(hide)) {
        names(Breaks)[which(names(Breaks) %in% hide)] = rep('', length(hide))
    }
    Breaks
}





# Datasets that we won't include:
#	- Breast datasets from Chung and Karaayvaz - CNVs often indistinct and low cell number.  We have enough good samples from the Qian breast dataset.
#	- CRC datasets KUL3 and Qian are OK, but SMC gives more usable samples, and I think is better quality based on the genes detected curve.
#	- Lung data from Song: infercna works quite well, but we still have fewer samples and fewer cells than the Qian data.
#	- Izar ovarian 10x data (Qian 10x data is better and has at least 1 more usable sample) and SS2 data (I can't identify non-malignant cells well).
#	- PDAC data from Elyada - CNV inference didn't work well, and maybe have only 2 usable samples.

# Note that for Lung data from Qian et al., if we want to separate into LUAD and LUSC we might still want to use the CNVs from the full dataset to
# help identify the malignant cells for patient 1, since for this patient it looks a bit less clean when considering only the LUSC samples.  The
# normal reference might also help here, but it looks complicated because there are ambiguous cases.  Perhaps it's easier to leave patient 1 out.

# Samples to keep:

# Breast Qian: 42, 43, 47, 49, 51, 53, 54
# CRC Lee SMC: SMC01, SMC02, SMC04, SMC07, SMC08, SMC09, SMC11, SMC14, SMC15, SMC16, SMC17, SMC18, SMC20, SMC21, SMC23, SMC25
# HNSCC Puram: 26, 25, 22, 18, 20, 6, 5
# Liver Ma: C26, C25, H38, H37, C56, C46, H65, C66
# Lung Qian: All samples, except possibly 1
# Ovarian Qian: 11, 12, 13, 14
# PDAC Peng: T2, T3, T6, T7, T8, T9, T11, T13, T14, T15, T16, T17, T19, T20, T21, T22

cohort_data <- list(
	# breast_qian = list(
		# `42` = list(
			# candidate_malignant_branches = list(list(1, c(2, 1, 2))),
			# cna_cor_thresholds = list(c(0.4, 0.5)),
			# cna_signal_thresholds = list(c(50, 70)),
			# seed = 4939
		# ),
		# `43` = list(
			# candidate_malignant_branches = list(
				# list(c(2, 2, 2)),
				# list(c(2, 2, 1, 2))
			# ),
			# cna_cor_thresholds = list(
				# c(0.35, 0.45),
				# c(0.4, 0.5)
			# ),
			# cna_signal_thresholds = list(
				# c(45, 65),
				# c(50, 65)
			# ),
			# seed = 7160
		# ),
		# `47` = list(
			# candidate_malignant_branches = list(list(c(2, 2, 1))),
			# cna_cor_thresholds = list(c(0.4, 0.6)),
			# cna_signal_thresholds = list(c(70, 200)),
			# seed = 1383
		# ),
		# `49` = list(
			# candidate_malignant_branches = list(
				# list(1),
				# list(c(2, 1))
			# ),
			# cna_cor_thresholds = list(
				# c(0.4, 0.5),
				# c(0.4, 0.5)
			# ),
			# cna_signal_thresholds = list(
				# c(45, 60),
				# c(50, 70)
			# ),
			# seed = 8499
		# ),
		# `51` = list(
			# candidate_malignant_branches = list(list(c(2, 2))),
			# cna_cor_thresholds = list(c(0.4, 0.45)),
			# cna_signal_thresholds = list(c(45, 55)),
			# seed = 3821
		# ),
		# `53` = list(
			# candidate_malignant_branches = list(list(1)),
			# cna_cor_thresholds = list(c(0.45, 0.55)),
			# cna_signal_thresholds = list(c(55, 75)),
			# seed = 832
		# ),
		# `54` = list(
			# candidate_malignant_branches = list(list(c(2, 2, 2, 2))),
			# cna_cor_thresholds = list(c(0.3, 0.4)),
			# cna_signal_thresholds = list(c(40, 40)),
			# seed = 1909
		# )
	# ),
	# crc_lee_smc = list(
		# SMC01 = list(
			# candidate_malignant_branches = list(list(c(2, 2, 2, 2, 2, 1))),
			# cna_cor_thresholds = list(c(0.45, 0.55)),
			# cna_signal_thresholds = list(c(45, 50)),
			# seed = 5162
		# ),
		# SMC02 = list(
			# candidate_malignant_branches = list(list(c(2, 2, 2))),
			# cna_cor_thresholds = list(c(0.4, 0.5)),
			# cna_signal_thresholds = list(c(50, 65)),
			# seed = 8962
		# ),
		# SMC04 = list(
			# candidate_malignant_branches = list(list(c(2, 2, 2, 2))),
			# cna_cor_thresholds = list(c(0.4, 0.5)),
			# cna_signal_thresholds = list(c(40, 40)),
			# seed = 2334
		# ),
		# SMC07 = list(
			# candidate_malignant_branches = list(list(c(2, 2, 2, 2, 1, 1, 2), c(2, 2, 2, 2, 2), c(2, 2, 2, 2, 1, 2))),
			# cna_cor_thresholds = list(c(0.4, 0.5)),
			# cna_signal_thresholds = list(c(45, 45)),
			# seed = 8705
		# ),
		# SMC08 = list(
			# candidate_malignant_branches = list(list(c(2, 2, 2), c(2, 2, 1, 2))),
			# cna_cor_thresholds = list(c(0.4, 0.5)),
			# cna_signal_thresholds = list(c(50, 65)),
			# seed = 5152
		# ),
		# SMC09 = list(
			# candidate_malignant_branches = list(list(c(2, 2))),
			# cna_cor_thresholds = list(c(0.4, 0.5)),
			# cna_signal_thresholds = list(c(40, 40)),
			# seed = 6112
		# ),
		# SMC11 = list(
			# candidate_malignant_branches = list(list(c(2, 2, 2, 1, 2), c(2, 2, 2, 2))),
			# cna_cor_thresholds = list(c(0.4, 0.5)),
			# cna_signal_thresholds = list(c(50, 65)),
			# seed = 6490
		# ),
		# SMC14 = list(
			# candidate_malignant_branches = list(
				# list(c(2, 2, 2, 2, 2, 2)),
				# list(c(2, 2, 2, 2, 2, 1))
			# ),
			# cna_cor_thresholds = list(
				# c(0.45, 0.575),
				# c(0.45, 0.5)
			# ),
			# cna_signal_thresholds = list(
				# c(50, 60),
				# c(50, 60)
			# ),
			# seed = 1636
		# ),
		# SMC15 = list(
			# candidate_malignant_branches = list(list(c(2, 2, 2, 2, 2, 2, 2, 2))),
			# cna_cor_thresholds = list(c(0.4, 0.55)),
			# cna_signal_thresholds = list(c(50, 50)),
			# seed = 1253
		# ),
		# SMC16 = list(
			# candidate_malignant_branches = list(list(c(2, 2, 2, 2))),
			# cna_cor_thresholds = list(c(0.45, 0.525)),
			# cna_signal_thresholds = list(c(35, 35)),
			# seed = 1892
		# ),
		# # Maybe we should leave out SMC17 - there seems to be a problem with likely nonmalignant cells appearing to have a deletion on chromosome
        # # 19, which leads to negative correlations between these cells' CNA profiles and those of the candidate malignant cells.
		# SMC17 = list(
			# candidate_malignant_branches = list(list(c(2, 2, 2, 2), c(2, 2, 1))),
			# cna_cor_thresholds = list(c(0.4, 0.5)),
			# cna_signal_thresholds = list(c(30, 30)),
			# seed = 2597
		# ),
		# SMC18 = list(
			# candidate_malignant_branches = list(
				# list(c(1, 2, 2, 1, 2, 2, 1, 2, 1, 2, 2)),
				# list(c(1, 2, 2, 1, 2, 2, 1, 2, 1, 1), c(1, 2, 2, 1, 2, 2, 1, 2, 2, 2, 2))
			# ),
			# cna_cor_thresholds = list(
				# c(0.45, 0.55),
				# c(0.45, 0.55)
			# ),
			# cna_signal_thresholds = list(
				# c(55, 55),
				# c(55, 65)
			# ),
			# seed = 3854
		# ),
		# SMC20 = list(
			# candidate_malignant_branches = list(list(c(2, 2, 2, 2, 2))),
			# cna_cor_thresholds = list(c(0.45, 0.55)),
			# cna_signal_thresholds = list(c(45, 55)),
			# seed = 7843
		# ),
		# SMC21 = list(
			# candidate_malignant_branches = list(
				# list(c(2, 2, 2)),
				# list(c(2, 2, 1)),
				# list(c(2, 1, 2))
			# ),
			# cna_cor_thresholds = list(
				# c(0.45, 0.55),
				# c(0.475, 0.625),
				# c(0.5, 0.6)
			# ),
			# cna_signal_thresholds = list(
				# c(50, 60),
				# c(70, 90),
				# c(70, 90)
			# ),
			# seed = 4694
		# ),
		# SMC23 = list(
			# candidate_malignant_branches = list(list(c(2, 1, 2), c(2, 2))),
			# cna_cor_thresholds = list(c(0.45, 0.55)),
			# cna_signal_thresholds = list(c(40, 40)),
			# seed = 9108
		# ),
		# SMC25 = list(
			# candidate_malignant_branches = list(list(c(2, 2))),
			# cna_cor_thresholds = list(c(0.45, 0.55)),
			# cna_signal_thresholds = list(c(40, 40)),
			# seed = 3655
		# )
	# ),
	# hnscc_puram = list(
		# `5` = list(
			# candidate_malignant_branches = list(list(c(2, 2, 2, 2))),
			# cna_cor_thresholds = list(c(0.4, 0.5)),
			# cna_signal_thresholds = list(c(35, 45)),
			# seed = 3039
		# ),
		# `6` = list(
			# candidate_malignant_branches = list(
				# list(c(2, 2, 2)),
				# list(c(2, 2, 1))
			# ),
			# cna_cor_thresholds = list(
				# c(0.4, 0.55),
				# c(0.4, 0.5)
			# ),
			# cna_signal_thresholds = list(
				# c(40, 50),
				# c(40, 50)
			# ),
			# seed = 2448
		# ),
		# `18` = list(
			# candidate_malignant_branches = list(list(c(2, 2, 1))),
			# cna_cor_thresholds = list(c(0.4, 0.525)),
			# cna_signal_thresholds = list(c(40, 50)),
			# seed = 6782
		# ),
		# `20` = list(
			# candidate_malignant_branches = list(list(c(2, 2, 2, 1, 1), c(2, 2, 2, 1, 2, 2))),
			# cna_cor_thresholds = list(c(0.35, 0.5)),
			# cna_signal_thresholds = list(c(40, 40)),
			# seed = 137
		# ),
		# `22` = list(
			# candidate_malignant_branches = list(list(1)),
			# cna_cor_thresholds = list(c(0.375, 0.525)),
			# cna_signal_thresholds = list(c(35, 35)),
			# seed = 3581
		# ),
		# `25` = list(
			# candidate_malignant_branches = list(list(c(2, 2, 2, 2, 2))),
			# cna_cor_thresholds = list(c(0.375, 0.475)),
			# cna_signal_thresholds = list(c(30, 30)),
			# seed = 3542
		# ),
		# `26` = list(
			# candidate_malignant_branches = list(list(c(2, 2, 2, 2, 2, 2, 2, 2))),
			# cna_cor_thresholds = list(c(0.375, 0.475)),
			# cna_signal_thresholds = list(c(35, 45)),
			# seed = 4096
		# )
	# ),
	# liver_ma = list(
		# C25 = list(
			# candidate_malignant_branches = list(list(c(2, 1, 1))),
			# cna_cor_thresholds = list(c(0.4, 0.5)),
			# cna_signal_thresholds = list(c(35, 50)),
			# seed = 5864
		# ),
		# C26 = list(
			# candidate_malignant_branches = list(list(2)),
			# cna_cor_thresholds = list(c(0.4, 0.5)),
			# cna_signal_thresholds = list(c(50, 70)),
			# seed = 7863
		# ),
		# C46 = list(
			# candidate_malignant_branches = list(
				# list(c(2, 2, 1, 2, 2, 2, 2, 1, 2, 2, 2), c(2, 2, 1, 2, 2, 2, 2, 1, 1), c(2, 2, 1, 2, 2, 2, 2, 2)),
				# list(c(2, 2, 1, 2, 2, 1))
			# ),
			# cna_cor_thresholds = list(
				# c(0.375, 0.525),
				# c(0.4, 0.65)
			# ),
			# cna_signal_thresholds = list(
				# c(60, 70),
				# c(70, 100)
			# ),
			# seed = 3916
		# ),
		# C56 = list(
			# candidate_malignant_branches = list(list(c(2, 2))),
			# cna_cor_thresholds = list(c(0.45, 0.55)),
			# cna_signal_thresholds = list(c(50, 70)),
			# seed = 2522
		# ),
		# C66 = list(
			# candidate_malignant_branches = list(
				# list(c(1, 2)),
				# list(c(2, 2), c(1, 1))
			# ),
			# cna_cor_thresholds = list(
				# c(0.4, 0.5),
				# c(0.375, 0.5)
			# ),
			# cna_signal_thresholds = list(
				# c(50, 65),
				# c(40, 50)
			# ),
			# seed = 9478
		# ),
		# H37 = list(
			# candidate_malignant_branches = list(list(c(2, 2, 2), c(2, 2, 1, 2))),
			# cna_cor_thresholds = list(c(0.4, 0.5)),
			# cna_signal_thresholds = list(c(35, 50)),
			# seed = 5354
		# ),
		# H38 = list(
			# candidate_malignant_branches = list(list(2)),
			# cna_cor_thresholds = list(c(0.35, 0.55)),
			# cna_signal_thresholds = list(c(35, 45)),
			# seed = 4555
		# ),
		# H65 = list(
			# candidate_malignant_branches = list(list(2)),
			# cna_cor_thresholds = list(c(0.375, 0.5)),
			# cna_signal_thresholds = list(c(40, 50)),
			# seed = 6812
		# )
	# ),
	luad_kim = list(
		P0006 = list(
			candidate_malignant_branches = list(list(2)),
			cna_cor_thresholds = list(c(0.4, 0.5)),
			cna_signal_thresholds = list(c(45, 50)),
			seed = 2418
		),
		P0008 = list(
			candidate_malignant_branches = list(list(1)),
			cna_cor_thresholds = list(c(0.5, 0.6)),
			cna_signal_thresholds = list(c(45, 50)),
			seed = 1966
		),
		P0018 = list(
			candidate_malignant_branches = list(list(1, c(2, 2))),
			cna_cor_thresholds = list(c(0.45, 0.5)),
			cna_signal_thresholds = list(c(40, 50)),
			seed = 2133
		),
		P0019 = list(
			candidate_malignant_branches = list(list(2)),
			cna_cor_thresholds = list(c(0.45, 0.55)),
			cna_signal_thresholds = list(c(50, 60)),
			seed = 7146
		),
		P0020 = list(
			candidate_malignant_branches = list(list(c(2, 2))),
			cna_cor_thresholds = list(c(0.45, 0.55)),
			cna_signal_thresholds = list(c(45, 50)),
			seed = 300
		),
		P0025 = list(
			candidate_malignant_branches = list(list(2, c(1, 1))),
			cna_cor_thresholds = list(c(0.4, 0.5)),
			cna_signal_thresholds = list(c(40, 40)),
			seed = 3765
		),
		P0028 = list(
			candidate_malignant_branches = list(
				list(c(2, 2, 1)),
				list(c(2, 2, 2, 2, 1))
			),
			cna_cor_thresholds = list(
				c(0.4, 0.5),
				c(0.4, 0.5)
			),
			cna_signal_thresholds = list(
				c(40, 50),
				c(40, 50)
			),
			seed = 6318
		),
		P0030 = list(
			candidate_malignant_branches = list(list(c(2, 2))),
			cna_cor_thresholds = list(c(0.45, 0.525)),
			cna_signal_thresholds = list(c(40, 40)),
			seed = 3819
		),
		P0031 = list(
			candidate_malignant_branches = list(list(c(2, 2, 1, 2), c(2, 2, 2))),
			cna_cor_thresholds = list(c(0.45, 0.55)),
			cna_signal_thresholds = list(c(45, 55)),
			seed = 6547
		),
		P0034 = list(
			candidate_malignant_branches = list(
				list(c(2, 2, 2)),
				list(c(2, 2, 1, 2)),
				list(c(2, 2, 1, 1)),
				list(1)
			),
			cna_cor_thresholds = list(
				c(0.4, 0.5),
				c(0.4, 0.5),
				c(0.4, 0.5),
				c(0.45, 0.55)
			),
			cna_signal_thresholds = list(
				c(37, 37),
				c(35, 35),
				c(35, 35),
				c(50, 55)
			),
			seed = 9538
		),
		P1006 = list(
			candidate_malignant_branches = list(list(2)),
			cna_cor_thresholds = list(c(0.4, 0.5)),
			cna_signal_thresholds = list(c(40, 40)),
			seed = 7419
		),
		P1028 = list(
			candidate_malignant_branches = list(list(c(2, 2, 2, 2, 2, 2, 2))),
			cna_cor_thresholds = list(c(0.4, 0.5)),
			cna_signal_thresholds = list(c(40, 50)),
			seed = 7570
		),
		P1049 = list(
			candidate_malignant_branches = list(list(c(2, 2))),
			cna_cor_thresholds = list(c(0.45, 0.55)),
			cna_signal_thresholds = list(c(40, 40)),
			seed = 645
		),
		P1058 = list(
			candidate_malignant_branches = list(
				list(c(2, 2, 1, 2)),
				list(c(2, 1, 1)),
				list(c(2, 1, 2, 2)),
				list(c(2, 2, 2), c(2, 2, 1, 1))
			),
			cna_cor_thresholds = list(
				c(0.5, 0.6),
				c(0.5, 0.7),
				c(0.55, 0.65),
				c(0.45, 0.55)
			),
			cna_signal_thresholds = list(
				c(55, 65),
				c(55, 80),
				c(60, 120),
				c(45, 55)
			),
			seed = 2250
		)
	)#,
	# lung_qian = list(
		# `1` = list(
			# candidate_malignant_branches = list(list(c(2, 2, 2))),
			# cna_cor_thresholds = list(c(0.45, 0.575)),
			# cna_signal_thresholds = list(c(50, 70)),
			# seed = 1828
		# ),
		# `2` = list(
			# candidate_malignant_branches = list(list(2)),
			# cna_cor_thresholds = list(c(0.475, 0.6)),
			# cna_signal_thresholds = list(c(50, 70)),
			# seed = 656
		# ),
		# `3` = list(
			# candidate_malignant_branches = list(list(c(2, 2, 2, 2))),
			# cna_cor_thresholds = list(c(0.4, 0.55)),
			# cna_signal_thresholds = list(c(50, 70)),
			# seed = 531
		# ),
		# `4` = list(
			# candidate_malignant_branches = list(
				# list(c(2, 2)),
				# list(c(2, 1, 2))
			# ),
			# cna_cor_thresholds = list(
				# c(0.4, 0.525),
				# c(0.4, 0.525)
			# ),
			# cna_signal_thresholds = list(
				# c(40, 50),
				# c(40, 55)
			# ),
			# seed = 354
		# ),
		# `5` = list(
			# candidate_malignant_branches = list(list(c(2, 2, 2, 2, 1, 2, 2, 1), c(2, 2, 2, 2, 1, 1))),
			# cna_cor_thresholds = list(c(0.375, 0.5)),
			# cna_signal_thresholds = list(c(40, 40)),
			# seed = 1920
		# ),
		# `6` = list(
			# candidate_malignant_branches = list(list(c(2, 1, 2, 2, 2), c(2, 1, 2, 1))),
			# cna_cor_thresholds = list(c(0.4, 0.525)),
			# cna_signal_thresholds = list(c(40, 40)),
			# seed = 4757
		# ),
		# `7` = list(
			# candidate_malignant_branches = list(list(c(2, 2))),
			# cna_cor_thresholds = list(c(0.375, 0.475)),
			# cna_signal_thresholds = list(c(40, 40)),
			# seed = 8569
		# ),
		# `8` = list(
			# candidate_malignant_branches = list(
				# list(c(2, 2, 2, 2, 1)),
				# list(c(2, 2, 2, 2, 2, 1, 1), c(2, 2, 2, 2, 2, 1, 2, 1))
			# ),
			# cna_cor_thresholds = list(
				# c(0.45, 0.575),
				# c(0.45, 0.525)
			# ),
			# cna_signal_thresholds = list(
				# c(50, 50),
				# c(45, 45)
			# ),
			# seed = 2945
		# )
	# ),
	# ovarian_qian = list(
		# `11` = list(
			# candidate_malignant_branches = list(
				# list(c(1, 2, 2, 2, 2, 2, 2, 1, 2, 2, 2)),
				# list(c(1, 2, 2, 2, 2, 2, 2, 1, 2, 2, 1))
			# ),
			# cna_cor_thresholds = list(
				# c(0.4, 0.5),
				# c(0.4, 0.5)
			# ),
			# cna_signal_thresholds = list(
				# c(50, 50),
				# c(50, 70)
			# ),
			# seed = 7375
		# ),
		# `12` = list(
			# candidate_malignant_branches = list(
				# list(c(2, 1, 1), c(2, 2, 2)),
				# list(c(1, 1))
			# ),
			# cna_cor_thresholds = list(
				# c(0.3, 0.4),
				# c(0.4, 0.5)
			# ),
			# cna_signal_thresholds = list(
				# c(25, 25),
				# c(45, 65)
			# ),
			# seed = 7662
		# ),
		# `13` = list(
			# candidate_malignant_branches = list(
				# list(c(2, 2, 2, 2)),
				# list(c(2, 2, 2, 1))
			# ),
			# cna_cor_thresholds = list(
				# c(0.35, 0.45),
				# c(0.4, 0.5)
			# ),
			# cna_signal_thresholds = list(
				# c(40, 50),
				# c(40, 50)
			# ),
			# seed = 8283
		# ),
		# `14` = list(
			# candidate_malignant_branches = list(list(2)),
			# cna_cor_thresholds = list(c(0.4, 0.5)),
			# cna_signal_thresholds = list(c(40, 50)),
			# seed = 3969
		# )
	# ),
	# pdac_peng = list(
		# T2 = list(
			# candidate_malignant_branches = list(list(2)),
			# cna_cor_thresholds = list(c(0.375, 0.5)),
			# cna_signal_thresholds = list(c(50, 50)),
			# seed = 3912
		# ),
		# T3 = list(
			# candidate_malignant_branches = list(list(c(2, 2, 2, 2, 2, 1, 2, 1))),
			# cna_cor_thresholds = list(c(0.45, 0.525)),
			# cna_signal_thresholds = list(c(50, 65)),
			# seed = 2064
		# ),
		# T6 = list(
			# candidate_malignant_branches = list(list(c(2, 2, 2, 2, 2))),
			# cna_cor_thresholds = list(c(0.35, 0.475)),
			# cna_signal_thresholds = list(c(45, 45)),
			# seed = 6041
		# ),
		# T7 = list(
			# candidate_malignant_branches = list(list(c(2, 2))),
			# cna_cor_thresholds = list(c(0.35, 0.5)),
			# cna_signal_thresholds = list(c(40, 50)),
			# seed = 7922
		# ),
		# T8 = list(
			# candidate_malignant_branches = list(
				# list(c(2, 1, 2, 2)),
				# list(c(2, 1, 2, 1, 2))
			# ),
			# cna_cor_thresholds = list(
				# c(0.425, 0.525),
				# c(0.45, 0.55)
			# ),
			# cna_signal_thresholds = list(
				# c(50, 70),
				# c(50, 60)
			# ),
			# seed = 1207
		# ),
		# T9 = list(
			# candidate_malignant_branches = list(list(c(2, 2, 2))),
			# cna_cor_thresholds = list(c(0.4, 0.5)),
			# cna_signal_thresholds = list(c(45, 60)),
			# seed = 8948
		# ),
		# T11 = list(
			# candidate_malignant_branches = list(list(c(2, 2, 1))),
			# cna_cor_thresholds = list(c(0.4, 0.5)),
			# cna_signal_thresholds = list(c(45, 45)),
			# seed = 9533
		# ),
		# T13 = list(
			# candidate_malignant_branches = list(list(2)),
			# cna_cor_thresholds = list(c(0.45, 0.55)),
			# cna_signal_thresholds = list(c(50, 60)),
			# seed = 1866
		# ),
		# T14 = list(
			# candidate_malignant_branches = list(list(c(2, 2))),
			# cna_cor_thresholds = list(c(0.45, 0.525)),
			# cna_signal_thresholds = list(c(50, 70)),
			# seed = 8607
		# ),
		# T15 = list(
			# candidate_malignant_branches = list(list(1)),
			# cna_cor_thresholds = list(c(0.425, 0.55)),
			# cna_signal_thresholds = list(c(65, 85)),
			# seed = 4024
		# ),
		# T16 = list(
			# candidate_malignant_branches = list(list(c(2, 2, 2))),
			# cna_cor_thresholds = list(c(0.45, 0.55)),
			# cna_signal_thresholds = list(c(55, 70)),
			# seed = 845
		# ),
		# T17 = list(
			# candidate_malignant_branches = list(
				# list(c(2, 2, 1, 1)),
				# list(c(2, 2, 1, 2))
			# ),
			# cna_cor_thresholds = list(
				# c(0.45, 0.55),
				# c(0.45, 0.55)
			# ),
			# cna_signal_thresholds = list(
				# c(55, 70),
				# c(55, 70)
			# ),
			# seed = 8479
		# ),
		# T19 = list(
			# candidate_malignant_branches = list(list(c(2, 1))),
			# cna_cor_thresholds = list(c(0.45, 0.525)),
			# cna_signal_thresholds = list(c(45, 55)),
			# seed = 1801
		# ),
		# T20 = list(
			# candidate_malignant_branches = list(list(c(2, 2))),
			# cna_cor_thresholds = list(c(0.45, 0.525)),
			# cna_signal_thresholds = list(c(50, 50)),
			# seed = 8253
		# ),
		# T21 = list(
			# candidate_malignant_branches = list(list(c(2, 1, 2))),
			# cna_cor_thresholds = list(c(0.4, 0.5)),
			# cna_signal_thresholds = list(c(45, 45)),
			# seed = 1306
		# ),
		# T22 = list(
			# candidate_malignant_branches = list(
				# list(c(2, 2, 2, 2)),
				# list(c(2, 2, 1))
			# ),
			# cna_cor_thresholds = list(
				# c(0.45, 0.55),
				# c(0.425, 0.525)
			# ),
			# cna_signal_thresholds = list(
				# c(50, 65),
				# c(50, 65)
			# ),
			# seed = 4752
		# )
	# )
)





for(cohort in names(cohort_data)) {

	cat(paste0(cohort, ':\n\tReading data\n'))

	inferred_cna <- readRDS(paste0('../data_and_figures/inferred_cna_', cohort, '.rds'))
	cell_clust <- readRDS(paste0('../data_and_figures/cna_clust_', cohort, '.rds'))

	x_line_breaks <- .chromosomeBreaks(rownames(inferred_cna[[1]]))
	x_text_breaks <- round(.chromosomeBreaks(rownames(inferred_cna[[1]]), halfway = TRUE, hide = c('13', '18', '21', 'Y')))

	for(p in names(cohort_data[[cohort]])) {

		cat(paste0('\tPatient ', p, ':\n\t\tComputing CNA scores and correlation\n'))

		candidate_malignant <- lapply(
			cohort_data[[cohort]][[p]]$candidate_malignant_branches,
			function(x) lapply(x, function(y) labels(reduce(y, .f = `[[`, .init = as.dendrogram(cell_clust[[p]]))))
		) # This list will have length equal to the number of suspected subclones

		candidate_malignant <- lapply(candidate_malignant, unlist) # In case multiple branches correspond to the same subclone

		assumed_nonmalignant <- colnames(inferred_cna[[p]])[!(colnames(inferred_cna[[p]]) %in% unlist(candidate_malignant))]

		classification_analysis <- lapply(
			1:length(candidate_malignant), # Iterate over the different subclones (returns list of length 1 if there's only one clone)
			function(i) {

				candidate_malignant_mean_cna <- rowMeans(inferred_cna[[p]][, candidate_malignant[[i]]])

				cna_cor <- apply(
					inferred_cna[[p]][, c(assumed_nonmalignant, candidate_malignant[[i]])],
					2,
					cor,
					candidate_malignant_mean_cna
				)

				cna_signal <- apply(
					inferred_cna[[p]][, c(assumed_nonmalignant, candidate_malignant[[i]])],
					2,
					function(x) sum(x[candidate_malignant_mean_cna^2 >= quantile(candidate_malignant_mean_cna^2, 0.9)]^2)
				)

				annotation_data_ref <- data.table(
					cell_id = c(assumed_nonmalignant, candidate_malignant[[i]]),
					cna_cor = cna_cor,
					cna_signal = cna_signal
				)[
					,
					classification := switch( # If it's below both thresholds, then it's nonmalignant.  If not...
						(cna_cor < cohort_data[[cohort]][[p]]$cna_cor_thresholds[[i]][1] &
                            cna_signal < cohort_data[[cohort]][[p]]$cna_signal_thresholds[[i]][1]) + 1,
						switch( # ...Then if it's above both thresholds, it's malignant; otherwise, it's ambiguous.
							(cna_cor > cohort_data[[cohort]][[p]]$cna_cor_thresholds[[i]][2] &
                                cna_signal > cohort_data[[cohort]][[p]]$cna_signal_thresholds[[i]][2]) + 1,
							'ambiguous',
							switch( # Call it "malignant" if there's only one clone; "malignant clone i" if there are multiple subclones.
                                (length(candidate_malignant) > 1) + 1,
                                'malignant', # If candidate_malignant has length 1, there's only 1 clone
                                paste('malignant clone', i) # If it has length > 1, there are multiple subclones.
                            )
						),
						'nonmalignant'
					),
					by = cell_id
				]

				cna_cor_signal_scatterplot <- ggplot(annotation_data_ref, aes(x = cna_cor, y = cna_signal, colour = classification)) +
					geom_point(alpha = 0.5) +
					scale_colour_manual(
						values = setNames(
							c('orangered3', 'grey', 'palegreen3'),
							c(switch((length(candidate_malignant) > 1) + 1, 'malignant', paste('malignant clone', i)), 'ambiguous', 'nonmalignant')
						)
					) +
					geom_vline(xintercept = cohort_data[[cohort]][[p]]$cna_cor_thresholds[[i]], colour = c('palegreen4', 'orangered4'), size = 0.3,
                        linetype = 'dashed') +
					geom_hline(yintercept = cohort_data[[cohort]][[p]]$cna_signal_thresholds[[i]], colour = c('palegreen4', 'orangered4'), size = 0.3,
                        linetype = 'dashed') +
					theme_test() +
					labs(x = 'CNA correlation', y = 'CNA signal', title = p)

				list(data = annotation_data_ref, plot = cna_cor_signal_scatterplot)

			}
		)

		# for(x in classification_analysis) print(x$plot)
		# dev.off()

        # If there were multiple subclones, bind the multiple tables together:
		if(length(classification_analysis) == 1) {
			classification_data <- classification_analysis[[1]]$data
		} else {
			classification_data <- rbindlist(
				lapply(
					1:length(classification_analysis),
					function(i) cbind(classification_analysis[[i]]$data, clone = i)
				)
			)
		}

        # The classification data still only contains data for one patient.  So a cell should appear multiple times in this table if and only if
        # it was tested against multiple subclones, in which case, all cells will appear the same number of times, equal to the number of subclones.
        # It therefore seems wasteful to check length(classification) per cell_id in the below...

		cat('\t\tDefining final classifications\n')

		classification_data[
			,
			classification_final := switch(
				(length(classification) == 1) + 1,
				switch( # This should happen only where we're testing against multiple subclones.
					(sum(grep('^malignant', classification)) > 0) + 1,
					switch( # In this case, the cell isn't classified as malignant for any subclone but could have been classified as ambiguous.
						('ambiguous' %in% classification) + 1,
						switch(
							(!(unique(cell_id) %in% cell_clust[[p]]$labels)) + 1,
							'nonmalignant',
							'reference'
						),
						'ambiguous'
					),
					switch( # In this case, the cell was classified as malignant for at least one subclone.
						(sum(grep('^malignant', classification)) == 1) + 1,
						'ambiguous', # If it's classified as malignant w.r.t. more than one subclone, that makes it ambiguous.
						classification[grep('^malignant', classification)]
					)
				),
				switch( # In this case there's only one clone, and we just need to change classification to "refernece" where appropriate.
					(!(cell_id %in% cell_clust[[p]]$labels) && classification == 'nonmalignant') + 1,
					classification,
					'reference'
				)
			),
			keyby = cell_id
		]

        # Where multiple subclones were tested, classification_final will have multiple repeated entries, so take the unique ones:
		classification_data_unique <- unique(classification_data[, .(cell_id, classification_final)])

		# Cell IDs are in alphabetical order within each classification, including when manually ordering the classifications:
		# classification_data_unique[, sum(cell_id == sort(cell_id)) == .N, by = classification_final]
		# classification_data_unique[
        #     c('malignant', 'ambiguous', 'nonmalignant', 'reference'),
        #     sum(cell_id == sort(cell_id)) == .N,
        #     by = classification_final
        # ]

		setkey(classification_data_unique, classification_final)

		cat('\t\tMaking plots\n')

		set.seed(cohort_data[[cohort]][[p]]$seed)

		cell_ids_final <- classification_data_unique[
			,
			c(
				cell_id[classification_final != 'reference'],
				sample(
					cell_id[classification_final == 'reference'],
					min(sum(classification_final == 'reference'), round(sum(classification_final != 'reference')/2))
				)
			)
		]

		clone_names <- classification_data_unique[grep('^malignant', classification_final), sort(unique(classification_final))]

		cell_ids_final <- classification_data_unique[
			c(clone_names, 'ambiguous', 'nonmalignant', 'reference'),
			cell_id[cell_id %in% cell_ids_final]
		]

		y_line_breaks <- classification_data_unique[
			c(clone_names, 'ambiguous', 'nonmalignant', 'reference')
		][
			cell_id %in% cell_ids_final,
			.(brks = .N),
			by = classification_final
		]$brks %>% cumsum

		y_text_breaks <- classification_data_unique[
			c(clone_names, 'ambiguous', 'nonmalignant', 'reference')
		][
			cell_id %in% cell_ids_final,
			.(line_break = .N, text_break = .N/2),
			by = classification_final
		][
			,
			round(sapply(1:.N, function(i) {switch((i == 1) + 1, cumsum(line_break)[i - 1], 0) + text_break[i]}))
		]

		cna_heatmap_final <- ggplot(
			reshape2::melt(
				inferred_cna[[p]][, cell_ids_final],
				varnames = c('gene', 'cell'),
				value.name = 'cna_score'
			)
		) +
			geom_raster(
				aes(
					x = factor(gene, levels = rownames(inferred_cna[[p]])),
					y = factor(cell, levels = cell_ids_final),
					fill = cna_score
				)
			) +
			scale_fill_gradientn(
				colours = rev(colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))(50)),
				limits = c(-1, 1),
				oob = scales::squish
			) +
			scale_y_discrete(
				expand = c(0, 0),
				breaks = cell_ids_final[y_text_breaks],
				labels = c(str_to_title(clone_names), 'Ambiguous', 'Nonmalignant', 'Reference')
			) +
			scale_x_discrete(
				expand = c(0, 0),
				breaks = rownames(inferred_cna[[p]])[x_text_breaks],
				labels = names(x_text_breaks)
			) +
			geom_vline(xintercept = x_line_breaks, size = 0.25) +
			geom_hline(yintercept = y_line_breaks, size = 0.25) +
			theme(
				axis.ticks = element_blank(),
				axis.ticks.length = unit(0, 'pt'),
				panel.border = element_rect(colour = 'black', size = 0.5, fill = NA)
			) +
			labs(x = 'Chromosome', y = 'Cells', fill = 'Inferred CNA', title = p)

		cat('\t\tSaving data and plots\n')

        # Saving classification_data: I guess I opted to save classification_data instead of classification_data_unique, because the former contains
        # subclone and CNA score/correlation information, and we can always reconstruct classification_data_unique.
		fwrite(classification_data, paste0('../data_and_figures/sc_find_malignant/', cohort, '/', p, '_data.csv'))

		pdf(paste0('../data_and_figures/sc_find_malignant/', cohort, '/', p, '_figures.pdf'), width = 10, height = 8)
		for(x in classification_analysis) print(x$plot)
		print(cna_heatmap_final)
		dev.off()

		cat('\t\tDone!\n')

	}

}





# The following is for making a CNA heatmap with dendrogram and annotation bar labelling candidate malignant and nonmalignant cells.

# cna_dendro <- dendro(cell_clust[[p]], edge = 'right')

# cna_heatmap <- ggplot(
	# reshape2::melt(
		# inferred_cna[[p]][, cell_clust[[p]]$labels],
		# varnames = c('gene', 'cell'),
		# value.name = 'cna_score'
	# )
# ) +
	# geom_raster(
		# aes(
			# x = factor(gene, levels = rownames(inferred_cna[[p]])),
			# y = factor(cell, levels = cell_clust[[p]]$labels[cell_clust[[p]]$order]),
			# fill = cna_score
		# )
	# ) +
	# scale_fill_gradientn(
		# colours = rev(colorRampPalette(RColorBrewer::brewer.pal(11, "RdBu"))(50)),
		# limits = c(-1, 1),
		# oob = scales::squish
	# ) +
	# scale_y_discrete(expand = c(0, 0)) +
	# scale_x_discrete(expand = c(0, 0), breaks = rownames(inferred_cna[[p]])[x_text_breaks], labels = names(x_text_breaks)) +
	# geom_vline(xintercept = x_line_breaks, size = 0.25) +
	# theme(
		# axis.text.y = element_blank(),
		# axis.title.y = element_blank(),
		# axis.ticks = element_blank(),
		# axis.ticks.length = unit(0, 'pt'),
		# panel.border = element_rect(colour = 'black', size = 0.5, fill = NA),
		# plot.margin = unit(c(5.5, 5.5, 5.5, 0), 'pt')
	# ) +
	# labs(x = 'Chromosome', fill = 'Inferred CNA', title = p)

# annotation_data <- data.table(cell_id = cell_clust[[p]]$labels)[
	# ,
	# classification := switch(
		# (cell_id %in% candidate_malignant) + 1,
		# 'nonmalignant',
		# 'malignant'
	# ),
	# by = cell_id
# ]

# cna_annotation <- ggplot(annotation_data) +
	# geom_raster(aes(x = '0', y = factor(cell_id, levels = with(cell_clust[[p]], labels[order])), fill = classification)) +
	# scale_x_discrete(expand = c(0, 0)) +
	# scale_y_discrete(expand = c(0, 0)) +
	# theme(
		# axis.text = element_blank(),
		# axis.title = element_blank(),
		# axis.ticks = element_blank(),
		# axis.ticks.length = unit(0, 'pt'),
		# plot.margin = unit(c(5.5, 0, 5.5, 0), 'pt')
	# )

# combined_plot <- plot_grid(
	# cna_dendro,
	# cna_annotation + theme(legend.position = 'none'),
	# cna_heatmap + theme(legend.position = 'none'),
	# plot_grid(get_legend(cna_annotation), get_legend(cna_heatmap), nrow = 2, ncol = 1),
	# nrow = 1,
	# ncol = 4,
	# align = 'h',
	# rel_widths = c(1, 0.1, 5, 1)
# )
