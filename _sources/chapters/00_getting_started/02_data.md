# Atacama Soil Microbiome

In our example, we showcase the application of our QIIME2 plugins for high-dimensional statistics using the Atacama soil microbiome dataset {cite}`neilson2017significant`. With q2-gglasso we solve various graphical lasso problems to identify microbial associations which we later assess by fitting sparse log-contrast models implemented in q2-classo. Microbiome bioinformatics analyses were conducted using QIIME 2 version 2022.4 {cite}`bolyen2019reproducible`. The processing of raw sequence data involved demultiplexing and quality filtering, facilitated by the q2-demux plugin. Subsequent denoising was performed using DADA2 {cite}`callahan2016dada2` through the q2-dada2 plugin. Taxonomic assignments for Amplicon Sequence Variants (ASVs) were accomplished using the q2-feature-classifier {cite}`bokulich2018q2` with the naive Bayes taxonomy classifier, trained on the Silva Database {cite}`quast2012silva`.

## Data Preprocessing Pipeline

Here is the original code for Atacama example data [preprocessing](https://docs.qiime2.org/2024.10/tutorials/atacama-soils/). One can simply follow the original Atacama tutorial, but we wanted to give these ASVs meaningful names for the sake of this example. That's why we rename the ASV keys as follows:

| ASV label | Real ASV ID                      |
| --------- | -------------------------------- |
| ASV-1     | 89cb1ddf89dcf11d86f725dbcaa9a5ce |
| ASV-2     | cdb9c0ee3bba4c3d8b9292eb575bd9e3 |
| ASV-3     | 4c8ff0ea98d2c0ebb486e33bd96c61f9 |
| ASV-4     | a56b903521b8ab33a7878b35582803a8 |
| ASV-5     | 5c78314ff92e6fec9aa07acc1fa0dc24 |
| ASV-6     | dc8a2f47b3d1dc2e1f5f805891976b29 |
| ASV-7     | 6b780e361cfc5f06def718518324bdcb |
| ASV-8     | f6c10a04d57159c0d64d6bc30c677471 |
| ASV-9     | ffd60d684f32e6fd5b47fe90095f9d34 |
| ASV-10    | ef3fdbe1dcde754d91130cde6a4b4d61 |
| ASV-11    | a36b38f754f6abd278aeb9dbc7696343 |
| ASV-12    | a7b877ae6d2f079a15b6b192a4425620 |
| ASV-13    | 409faa5f5353e543bf6d99125c7c0e83 |

## Data for Downstream analysis

We'll use the Atacama soil microbiome dataset {cite}`neilson2017significant`, which contains:
- 49 samples from Atacama Desert soil;
- 13 microbial taxa (ASVs)
- Environmental covariates: pH, elevation, temperature, humidity, and vegetayion

The original data is available through the European Nucleotide Archive under accession [ERP019482](https://www.ebi.ac.uk/ena/browser/view/PRJEB17617).

| Sample ID | ASV-1 | ASV-2 | ASV-3 | ... | ASV-11 | ASV-12 | ASV-13 | Elevation | pH | Avg Soil RH | Avg Soil Temp | Vegetation |
|-----------|-------|-------|-------|-----|--------|--------|--------|-----------|----|-----------|--------------|-----------| 
| BAQ2420.1.1 | 0.0 | 11.0 | 0.0 | ... | 115.0 | 0.0 | 0.0 | 2420 | 9.33 | 82.54 | 22.45 | no |
| BAQ2420.1.2 | 0.0 | 0.0 | 0.0 | ... | 0.0 | 0.0 | 0.0 | 2420 | 9.36 | 82.54 | 22.45 | no |
| BAQ2420.1.3 | 0.0 | 0.0 | 0.0 | ... | 0.0 | 0.0 | 0.0 | 2420 | 8.90 | 82.54 | 22.45 | no |
| ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... | ... |
| YUN3856.2 | 6.0 | 13.0 | 26.0 | ... | 0.0 | 0.0 | 104.0 | 3856 | 7.43 | 99.44 | 9.51 | yes |
| YUN3856.3 | 21.0 | 36.0 | 23.0 | ... | 33.0 | 0.0 | 0.0 | 3856 | 7.43 | 99.44 | 9.51 | yes |

QIIME2 .qza file can be downloaded from this [link](https://github.com/bio-datascience/q2-gglasso/tree/master/data/atacama-counts.qza), here is the snapshot of the count data and corresponding [metadata](https://data.qiime2.org/2024.10/tutorials/atacama-soils/sample_metadata.tsv) we're using in our example.
