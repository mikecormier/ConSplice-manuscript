# Scripts used for the ConSplice manuscript

This directory contains scripts used for data creation, collection, and anlaysis. 

## Script directory description

based on order of directories listed

| Directory | Description |
| --------- | ----------- |
| HGMD_Pathogenic_Splicing_DM | The HGMD SQL query used to download the disease causing mutations (DM) Pathogenic Splice Altering variants from the professional version of HGMD using in **Figure 5** | 
| benign_variants | Information on were to download the variants used as the benign truth set as seen in **Figure 5** of the ConSplice manuscript | 
| constitutive_cassette_exons | A set of scripts used to generate the data used to evaluate the regional ConSplice model at exon features seen in **Figure 3** |
| dosage_sensitivity | Information on how to download the V<sup>G</sup> values from [Mohammadi et al.](https://www.science.org/doi/10.1126/science.aay0256?) used in **Figure 2**. |
| gene_sets | A script to download the gene sets used for **Figure 2** and **Supplemental Figure S2, S3** in the manuscript. This script will download the gene set data files into the *data* directory of this repo. | 
| sOutliers | A script to download and extract the sOutliers identified by [Ferraro et al.](https://www.science.org/doi/10.1126/science.aaz5900) using RNA-seq from GTEx v8. **Figure 2**. sOutliers will be downloaded to the *data* directory of this repo. | 





