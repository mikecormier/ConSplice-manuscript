
# Data donwloaded into the data directory of this repo


## HI, AD, AR, OR Gene lists

### Haploinsufficient Genes hosted on the MaCarthur lab gene list GitHub repo
wget -O ../../data/HI.genes.tsv https://raw.githubusercontent.com/macarthur-lab/gene_lists/master/lists/clingen_level3_genes_2018_09_13.tsv

### Autosomal Dominant Genes hosted on the MaCarthur lab gene list GitHub repo
wget -O ../../data/AD.genes.tsv https://raw.githubusercontent.com/macarthur-lab/gene_lists/master/lists/all_ad.tsv

### Autosomal Recessive Genes hosted on the MaCarthur lab gene list GitHub repo
wget -O ../../data/AR.genes.tsv https://raw.githubusercontent.com/macarthur-lab/gene_lists/master/lists/all_ar.tsv

### Olfactory Receptor Genes hosted on the MaCarthur lab gene list GitHub repo
wget -O ../../data/OR.genes.tsv https://raw.githubusercontent.com/macarthur-lab/gene_lists/master/lists/olfactory_receptors.tsv


## CRISPR Essential and Non-Essential Gene lists

### CRISPR Essential Genes hosted on the MaCarthur lab gene list GitHub repo
wget -O ../../data/CRISPR.Essential.genes.tsv https://raw.githubusercontent.com/macarthur-lab/gene_lists/master/lists/CEGv2_subset_universe.tsv

### CRISPR Non-Essential Genes hosted on the MaCarthur lab gene list GitHub repo
wget -O ../../data/CRISPR.Nonessential.genes.tsv https://raw.githubusercontent.com/macarthur-lab/gene_lists/master/lists/NEGv1_subset_universe.tsv


## Developmental Delay and Intelectual Disability Gene Lists

### DD/ID genes from the Deciphering Developmental Delay studies hosted on the Gene2Phenotype website
wget -O ../../data/DD.ID.genes.csv.gz https://www.ebi.ac.uk/gene2phenotype/downloads/DDG2P.csv.gz


## Homozygous LoF Tolerant genes from gnomAD
### Data from the Karczewki et al., Nature, 2020 on gnomAD. Supplemental Data Table 7: "supplementary_dataset_7_hom_ko_genes.txt" 
#### Download supplmentary data zip file
wget https://static-content.springer.com/esm/art%3A10.1038%2Fs41586-020-2308-7/MediaObjects/41586_2020_2308_MOESM4_ESM.zip

## Unzip it
unzip 41586_2020_2308_MOESM4_ESM.zip

## move the homozygous LoF genes to the data directory 
mv supplement/supplementary_dataset_7_hom_ko_genes.txt ../../data/Hom.LoF_tolerant.genes.tsv

## Remove unnecessary files. 
rm 41586_2020_2308_MOESM4_ESM.zip 
rm -r supplement/
rm -r __MACOSX/
