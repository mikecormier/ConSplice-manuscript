
# Data donwloaded into the data directory of this repo


## HI, AD, AR, OR Gene lists

### Haploinsufficient Genes hosted on the MaCarthur lab gene list GitHub repo
wget -O ../../data/HI.genes.tsv https://raw.githubusercontent.com/macarthur-lab/gene_lists/master/lists/clingen_level3_genes_2018_09_13.tsv

### Autosomal Dominant Genes hosted on the MaCarthur lab gene list GitHub repo
wget -O ../../data/AD.genes.tsv https://raw.githubusercontent.com/macarthur-lab/gene_lists/master/lists/berg_ad.tsv

### Autosomal Recessive Genes hosted on the MaCarthur lab gene list GitHub repo
wget -O ../../data/AR.genes.tsv https://raw.githubusercontent.com/macarthur-lab/gene_lists/master/lists/all_ar.tsv

### Olfactory Receptor Genes hosted on the MaCarthur lab gene list GitHub repo
wget -O ../../data/OR.genes.tsv https://raw.githubusercontent.com/macarthur-lab/gene_lists/master/lists/olfactory_receptors.tsv


## Homozygous Tolerant genes hosted on the MaCarthur lab gene list GitHub repo
wget -O ../../data/Hom.LoF_tolerant.genes.tsv https://raw.githubusercontent.com/macarthur-lab/gene_lists/master/lists/homozygous_lof_tolerant_twohit.tsv


## CRISPR Essential and Non-Essential Gene lists

### CRISPR Essential Genes hosted on the MaCarthur lab gene list GitHub repo
wget -O ../../data/CRISPR.Essential.genes.tsv https://raw.githubusercontent.com/macarthur-lab/gene_lists/master/lists/CEGv2_subset_universe.tsv

### CRISPR Non-Essential Genes hosted on the MaCarthur lab gene list GitHub repo
wget -O ../../data/CRISPR.Nonessential.genes.tsv https://raw.githubusercontent.com/macarthur-lab/gene_lists/master/lists/NEGv1_subset_universe.tsv


## Developmental Delay and Intelectual Disability Gene Lists

### DD/ID genes from the Deciphering Developmental Delay studies hosted on the Gene2Phenotype website
wget -O ../../data/DD.ID.genes.csv.gz https://www.ebi.ac.uk/gene2phenotype/downloads/DDG2P.csv.gz
