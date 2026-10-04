#!/usr/bin/env Rscript
## ===========================================================================
## 02_dapa_postmortem.R
##
## Differential 3'UTR APA between SCZ and controls in adult postmortem DLPFC
## (CommonMind Consortium, CMC; MSSM, Penn and Pitt), using the iN 3'UTR
## reference of module 02 (lifted over to hg19).
##
## Manuscript: Fig. 2e-g, Fig. 3a (DLPFC), Supplementary Table 3
##             (DiffAPA_iPSC_PM, PMAPA columns).
##
## CMC data are controlled access (Synapse syn2759792) and are not distributed
## here. Obtain them from the CMC, place the files below in one directory and
## pass it with --cmc.
##
## USAGE
##   Rscript 03_differential_apa/02_dapa_postmortem.R --cmc <dir> [--cores <n>] [--out <dir>]
##          [--ethnicity <value>] [--cases <Dx,Dx,...>] [--tag <name>]
##          [--residuals-out <file>]
##
##   --ethnicity      restrict to donors with this CMC Ethnicity (e.g. Caucasian)
##   --cases          CMC diagnoses analysed as cases (default SCZ; SCZ,BP,AFF = all cases)
##   --prs            restrict to donors listed in this PRS file (column "sample" =
##                    RNA-seq sample ID), i.e. donors with imputed genotypes passing QC
##   --tag            output file stem (default SCZ_vs_Control)
##   --residuals-out  also write the per-library residuals of the covariate-only
##                    model (individual-level CMC-derived data; keep private)
##
##   Runs used in the manuscript:
##     default                                   Fig. 2e-g, Suppl. Table 3
##     --ethnicity Caucasian --prs <file> --tag SCZ_vs_Control_EUR_PRS
##                                               DLPFC dAPA gene set for the PRS analysis
##                                               (European donors with PRS: 150 / 173)
##     --cases SCZ,BP,AFF --tag allCases_vs_Control --residuals-out <file>
##                                               per-donor APA for the PRS analysis
##
## INPUTS (in --cmc <dir>)
##   fc_CMC_utr_counts_hg19.rds                    featureCounts of all CMC DLPFC RNA-seq
##                                                 libraries on dominant_pas_utrs (hg19
##                                                 liftover); counts, UTR regions x libraries
##   readsPerSample_mapped_in_proper_pair.txt      reads per library (header; col 1 BAM
##                                                 name, col 2 reads)
##   CMC_MSSM-Penn-Pitt_Clinical.txt               clinical data
##   CMC_MSSM-Penn-Pitt_DLPFC_DNA_IlluminaOmniExpressExome_GemToolsAncestry.tsv
##   CMC_MSSM-Penn-Pitt_DLPFC_mRNA-metaData.csv    RNA metadata (RIN, library batch)
##   CMC_MSSM-Penn-Pitt_DLPFC_clusteredLIB.txt     library batch -> library cluster
##   TechnicalCovariate_CommonMind.txt             Picard CollectRnaSeqMetrics per library
##                                                 (column sampleId = RNA-seq sample ID)
##
## OUTPUT (results/dapa/)
##   dAPA_DLPFC_SCZ_vs_Control.txt   ID, method, n_ref, n_comp, beta, se, p_wald,
##                                   p_lrt, delta, mean_expr, qvalue, padj_BH
##                                   dAPA at q <= 0.1; delta > 0 = higher distal usage
##                                   in controls
##
## MODEL
##   As for iNs (R/dapa_functions.R), with covariates institution, sex, age at
##   death, PMI, RIN, RIN^2, library cluster, ancestry EV1-5 and technical PC1-3
##   (PCA of Picard RNA-seq metrics); >= 50 libraries per group with 0 < R < 1;
##   genes below this threshold in either group are not tested; no expression
##   filter.
## ===========================================================================

.here <- local({ a <- commandArgs(FALSE); f <- sub("^--file=", "", a[grepl("^--file=", a)])
                 if (length(f)) dirname(normalizePath(f)) else getwd() })
REPO <- normalizePath(file.path(.here, ".."))
source(file.path(REPO, "R", "apa_pas_functions.R"))   # get_arg(), need_file()
source(file.path(REPO, "R", "dapa_functions.R"))

CMC     <- get_arg("--cmc"); if (is.null(CMC)) stop("set the CMC data directory with --cmc <dir>")
CORES   <- as.integer(get_arg("--cores", "1"))
OUT_DIR <- get_arg("--out", file.path(REPO, "results", "dapa"))
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
cmc <- function(f) need_file(file.path(CMC, f), "--cmc")
ETHN    <- get_arg("--ethnicity")
CASES   <- strsplit(get_arg("--cases", "SCZ"), ",")[[1]]
TAG     <- get_arg("--tag", "SCZ_vs_Control")
RES_OUT <- get_arg("--residuals-out")
PRS_FILE <- get_arg("--prs")

GROUPS        <- c("Control", "SCZ")
MIN_PER_GROUP <- 50
EXCLUDE       <- "MSSM_RNA_PFC_353"          # RNA metadata missing
TECH_METRICS  <- c("PF_ALIGNED_BASES", "PCT_CODING_BASES", "PCT_UTR_BASES", "PCT_INTRONIC_BASES",
                   "PCT_INTERGENIC_BASES", "MEDIAN_5PRIME_BIAS", "MEDIAN_3PRIME_BIAS")

## ------------------------------------------------------------ counts
cnt <- as.matrix(readRDS(cmc("fc_CMC_utr_counts_hg19.rds")))
colnames(cnt) <- gsub("\\.", "_", gsub("[._]accepted[._]hits[._]sort[._]coord[._]bam$", "", colnames(cnt)))
dup <- unique(rownames(cnt)[duplicated(rownames(cnt))])      # regions split by liftover
cnt <- cnt[!rownames(cnt) %in% dup, ]

reads <- read.delim(cmc("readsPerSample_mapped_in_proper_pair.txt"), stringsAsFactors = FALSE)
reads$sample <- gsub("\\.", "_", gsub("[._]accepted[._]hits[._]sort[._]coord[._]bam$", "", reads[[1]]))

## ------------------------------------------------------------ covariates
clin <- read.delim(cmc("CMC_MSSM-Penn-Pitt_Clinical.txt"), stringsAsFactors = FALSE)
anc  <- read.delim(cmc("CMC_MSSM-Penn-Pitt_DLPFC_DNA_IlluminaOmniExpressExome_GemToolsAncestry.tsv"),
                   stringsAsFactors = FALSE)
rna  <- read.csv(cmc("CMC_MSSM-Penn-Pitt_DLPFC_mRNA-metaData.csv"), stringsAsFactors = FALSE)
lib  <- read.delim(cmc("CMC_MSSM-Penn-Pitt_DLPFC_clusteredLIB.txt"), stringsAsFactors = FALSE)
tech <- read.delim(cmc("TechnicalCovariate_CommonMind.txt"), stringsAsFactors = FALSE)

ids  <- setdiff(intersect(clin$DLPFC_RNA_Sequencing_Sample_ID, colnames(cnt)), EXCLUDE)
clin <- clin[match(ids, clin$DLPFC_RNA_Sequencing_Sample_ID), ]
anc  <- anc[match(clin$Genotyping_Sample_ID, anc$Genotyping_Sample_ID), ]
rna  <- rna[match(ids, rna$DLPFC_RNA_Sequencing_Sample_ID), ]
tech <- tech[match(ids, tech$sampleId), ]
stopifnot(!anyNA(anc$Genotyping_Sample_ID), !anyNA(rna$DLPFC_RNA_Sequencing_Sample_ID),
          !anyNA(tech$sampleId))

## technical PCs: PCA of Picard metrics (log10 aligned bases), all libraries
tm <- as.matrix(tech[, TECH_METRICS]); tm[, 1] <- log10(tm[, 1])
tpc <- stats::prcomp(tm, center = TRUE, scale. = TRUE)$x

cov <- data.frame(condition = clin$Dx,
                  SITE = relevel(factor(clin$Institution), ref = "Pitt"),
                  SEX  = relevel(factor(clin$Gender), ref = "Female"),
                  AOD  = as.numeric(sub("90+", "90", clin$Age_of_Death, fixed = TRUE)),
                  PMI  = clin$PMI_hrs,
                  RIN  = rna$DLPFC_RNA_isolation_RIN, RIN2 = rna$DLPFC_RNA_isolation_RIN^2,
                  setNames(anc[, paste0("EV.", 1:5)], paste0("EV", 1:5)),
                  clustLIB = relevel(factor(lib$Cluster[match(rna$DLPFC_RNA_Sequencing_Library_Batch,
                                                              lib$Library_Batch)]), ref = "base"),
                  PC1 = tpc[, 1], PC2 = tpc[, 2], PC3 = tpc[, 3],
                  row.names = ids)
COVARIATES <- setdiff(names(cov), "condition")

## ------------------------------------------------------------ ratios + test
ann <- read.delim(need_file(file.path(REPO, "results", "apa_pas_library", "dominant_pas_utrs_hg38.txt"),
                            "module 02_apa_pas_library"), stringsAsFactors = FALSE)
cnt <- cnt[intersect(rownames(cnt), ann$ID), ids]
sf  <- setNames(reads[match(ids, reads$sample), 2] / 1e6, ids)
ap  <- apa_ratios(cnt, ann, sf)

cov$condition[cov$condition %in% CASES] <- "SCZ"            # cases analysed as one group
keep <- cov$condition %in% GROUPS
if (!is.null(ETHN)) keep <- keep & clin$Ethnicity == ETHN
if (!is.null(PRS_FILE)) keep <- keep & ids %in% read.delim(need_file(PRS_FILE, "--prs"))$sample
s  <- ids[keep]
si <- cov[s, ]
d  <- diff_apa(ap$R[, s], si, COVARIATES, GROUPS, min_per_group = MIN_PER_GROUP,
               ones_usable = FALSE, fallback = "none", cores = CORES)
r  <- finalize_dapa(d, ap$expr, si, min_expr = 0)
r  <- r[order(r$p_wald), ]
f  <- file.path(OUT_DIR, paste0("dAPA_DLPFC_", TAG, ".txt"))
write.table(r, f, sep = "\t", quote = FALSE, row.names = FALSE)
if (!is.null(RES_OUT))
  write.table(data.frame(ID = rownames(attr(d, "residuals")), attr(d, "residuals"), check.names = FALSE),
              RES_OUT, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("%s (DLPFC): %d / %d libraries, %d genes tested, %d at q <= 0.1 -> %s\n",
            TAG, sum(si$condition == GROUPS[1]), sum(si$condition == GROUPS[2]), nrow(r),
            sum(r$qvalue <= 0.1), f))
