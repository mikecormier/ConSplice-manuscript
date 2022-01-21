# Splicing Outliers (sOutliers) 

This repo contains the script used to download the splicing outliers identified by [Ferraro et al.](https://www.science.org/doi/10.1126/science.aaz5900) and used to evaluate the ability of the genic ConSplice model to identify splicing specific constraint. **Figure 2** and **Supplemental Figure S3**



## Script definition 

| Directory | Description |
| --------- | ----------- |
| get_gtex_sOutliers.sh | A script to download the GTEx v8 sOutliers used for **Figure 2** and **Supplemental Figure S3** in the manuscript. This script will download the GTEx v8 Outliers, remove all but the sOutlier file, and move the sOutlier file to the *data* directory of this repo. The sOutlier file will be named `gtexV8.sOutlier.stats.globalOutliers.removed.txt.gz`. | 

## Running the scripts


```
bash get_gtex_sOutliers.sh
```

The files downloaded by this script will be downloaded to the [data](https://github.com/mikecormier/ConSplice-manuscript/tree/main/data) directory of this repo
