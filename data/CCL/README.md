# Molecular Profiles of cancer cell lines
This folder contains the molecular profiles of cancer cell lines in the raw and processed forms. File sizes are very large, so an ad-hoc script (`02_CCL_profiles.sh`) downloads and processed the data from the original data portals.

| Filename | Description |
| ----- | ----- |
| `CCLE_DepMap_18Q2_maf_20180502.txt` | Merged mutation calls (coding region, germline filtered) |
| `CCLE_copynumber_byGene_2013-12-03.txt` | Copy-number values per gene. Affymetrix SNP6.0 arrays |
| `CCLE_DepMap_18Q2_RNAseq_reads_20180502.gct` | CCLE RNAseq gene expression data (read count) |
| `portal-Avana-2018-05-10.csv` | CERES inferred gene effect matrix (Meyers et al. 2017) |
| `portal-RNAi_merged-2018-05-10.csv` | DEMETER2 inferred gene dependency scores (McFarland et al. 2018) |
