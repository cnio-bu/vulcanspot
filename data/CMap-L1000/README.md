# CMap-L1000 data for KD-CP perturbations
The file sizes are very large to be included in the repository. The data can be downloaded from https://figshare.com/s/b34fc8e7f4d77ca6d580

The files contains the TAU scores from the signature maching between knock-down and compound cellular perturbations in cancer cell lines. In brief, transcriptional signatures of gene knock-down (KD) and compound (CP) perturbations in cancer cell lines that belong to the same context was obtained using limma (https://github.com/jperales/CLUE). Then, a Total Enrichment Score (TES) was calculated by each pair of KD-CP (Iorio et al. PNAS 2010). Finally, the TES were scaled using the TAU approach proposed in Subramanian et al 2017.

Files:
* `BREAST.TAU_matrix.rds`
* `CMap_contexts.txt`
* `KIDNEY.TAU_matrix.rds`
* `LARGE_INTESTINE.TAU_matrix.rds`
* `LIVER.TAU_matrix.rds`
* `LUNG.TAU_matrix.rds`
* `PANCANCER.TAU_matrix.rds`
* `PROSTATE.TAU_matrix.rds`
* `SKIN.TAU_matrix.rds`
