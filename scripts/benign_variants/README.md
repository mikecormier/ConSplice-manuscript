# Benign Variants

Benign variants used to evaluate the performance of CADD, SpliceAI, SQUIRLS, and ConSpliceML in **Figure 5** of the manuscript. 

Variants are a combination of de novo mutations (DNM) from [CEPH](https://elifesciences.org/articles/46922) and [deCODE](https://www.nature.com/articles/nature24018), and [benign alternative splicing variants](https://www.sciencedirect.com/science/article/pii/S0092867418316295?via%3Dihub) validated in the GTEx RNA-seq data, with each dataset downloaded from the respective paper.

## Benign variants download

- GRCh37 DNM from CEPH were downloaded from the GitHub repo associated with the [Sasani et al.](https://elifesciences.org/articles/46922) paper [here](https://github.com/quinlan-lab/ceph-dnm-manuscript).

    - 2-generation DNM
      ```
      wget https://raw.githubusercontent.com/quinlan-lab/ceph-dnm-manuscript/master/data/second_gen.dnms.txt
      ```
    - 3-generation DNM
      ```
      wget https://raw.githubusercontent.com/quinlan-lab/ceph-dnm-manuscript/master/data/third_gen.dnms.txt
      ```
  DNM were then lifted over to GRCh38 using [CrossMap](http://crossmap.sourceforge.net/) 
    
- GRCh38 DNM from deCODE were downloaded from Supplemental Table S4 of the [Jónssonet al.](https://www.nature.com/articles/nature2401) paper

- GRCh37 Benign validated splicing altering variants in GTEx were downloaded from the [illumina basespace](https://basespace.illumina.com/s/otSPW8hnhaZR) from the [Jaganathan et al.](https://www.sciencedirect.com/science/article/pii/S0092867418316295?via%3Dihub) paper.  

    - Supplemental files representing the data from the basespace download:

      * Variants seen in 1-4 individuals in GTEx: `all/SpliceAI_supplement_v2_ds.eebebfd43c2343518df86c448517ce79/validated_preds.csv`
      * Private variants in GTEx: `all/SupplementaryData_ds.b56bfe22e7394bf39780f4312d544683/validated_preds.csv`

  The unique set of variants were then lifted over to GRCh38 using [CrossMap](http://crossmap.sourceforge.net/). 



## Filtering benign variants

Before filtering, the combined set contained 141,835 unique variants. 

After filtering, the combined set contained 48,980 benign variants. Those variants are available in the [data](https://github.com/mikecormier/ConSplice-manuscript/data/) directory of this repo

Variants were filtered as follows:

  - Variants in protein-coding genes
  - Variants in the autosome 
  - Variants that had scores from CADD, SpliceAI, SQUIRLS, and ConSplice. 

The majority of variants filtered out are found in intragenic regions where most of these tools do not score variants. 



## Preprocessed Data

To improve reproducability, we have preprocessed this benign set of variants and make it available in the [data](https://github.com/mikecormier/ConSplice-manuscript/tree/main/data) directory of this repo under the name [benign.combined.txt](https://github.com/mikecormier/ConSplice-manuscript/blob/main/data/benign.combined.txt).

