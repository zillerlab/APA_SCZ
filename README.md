# Code for: Polygenic risk for schizophrenia is associated with alternative polyadenylation as molecular mechanism contributing to synaptic impairment

Raabe, Atella, Hausruckinger, Gagliardi, Almeida et al.

Analysis code for the study. Standard preprocessing (alignment, peak calling,
quantification) follows the Supplementary Methods and is not included.

## Layout

```
<NN>_<module>/      analysis scripts and a README per module
R/                  functions shared between modules
resources/<module>/ processed, de-identified input data
results/<module>/   outputs (figures in results/figures/)
```

Scripts resolve paths relative to the repository root and can be run from any
directory. Samples are identified by the iPSC line IDs of Supplementary Table 1.

## Modules

| Module | Content | Manuscript |
| :--- | :--- | :--- |
| [`01_eb_scorecard`](01_eb_scorecard) | Tri-lineage differentiation of iPSC lines (embryoid-body RNA-seq) | Ext. Data Fig. 1g,h |
| [`02_apa_pas_library`](02_apa_pas_library) | PAS detection and two-PAS 3'UTR reference | Fig. 2a,b; Suppl. Table 3 |
| [`03_differential_apa`](03_differential_apa) | Differential APA: iN cohort, postmortem DLPFC, RBP knockdowns | Fig. 2c,e-g, 5e-h; Ext. Data Fig. 6b; Suppl. Table 3 |
| [`04_prs_association`](04_prs_association) | Polygenic risk and cumulative differential APA | Fig. 4b,d,e; Ext. Data Fig. 6c-i |
| [`05_differential_expression`](05_differential_expression) | Differential expression (mRNA, microRNA), power analysis | Fig. 1c-e, 5d; Ext. Data Fig. 4c,d,f,g, 7d |
| [`06_proteomics_olink`](06_proteomics_olink) | Protein abundance (Olink) after RBP knockdown | Fig. 5i,j |
| [`07_western_blot`](07_western_blot) | Western blot quantification of RBP knockdowns and GRIA1 3'UTR KO | Fig. 5c, 7e |
| [`08_chromatin_accessibility`](08_chromatin_accessibility) | ATAC-seq PCA and differential accessibility | Fig. 1c; Ext. Data Fig. 4e |
| [`09_ptbp2_clip`](09_ptbp2_clip) | PTBP2 CLIP-seq and overlap with differential APA | Fig. 6a |
| [`10_mea`](10_mea) | Multielectrode array recordings | Fig. 6f,g, 7d |
| [`11_gria1_reporter`](11_gria1_reporter) | GRIA1 qPCR: GRIA1 3'UTR KO iNs and GFP reporter in HEK293T cells | Fig. 7c,h |

## Data availability

Raw sequencing data: European Genome-phenome
Archive [TBA]. Polygenic risk scores and CommonMind data are available under
controlled access (see modules 03 and 04).

## License

MIT
