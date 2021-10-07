# ConSplice Data Files

Data files used for the ConSplice manuscript


## Gene sets:

The script used to download these files can be found in the [scripts](https://github.com/mikecormier/ConSplice-manuscript/scripts/) directory of this repo under the name [gene_sets.sh](https://github.com/mikecormier/ConSplice-manuscript/scripts/gene_sets.sh)

| Data File Name | Description | Data Host Link |
| -------------- | ----------- | -------------- |
| **AD.genes.tsv** | Autosomal Dominant genes | [AD MaCarthur Gene List](https://github.com/macarthur-lab/gene_lists/blob/master/lists/berg_ad.tsv) |
| **AR.genes.tsv** | Autosomal Recessive genes | [AR MaCarthur Gene List](https://github.com/macarthur-lab/gene_lists/blob/master/lists/all_ar.tsv) |
| **CRISPR.Essential.genes.tsv** | Genes determined to be essential in cell culture CRISPR screens | [CRISPER Ess. MaCarthur Gene List](https://github.com/macarthur-lab/gene_lists/blob/master/lists/CEGv2_subset_universe.tsv) |
| **CRISPR.Nonessential.genes.tsv** | Genes determined to be non-essential in cell culture CRISPR screens | [CRISPER NonEss. MaCarthur Gene List](https://github.com/macarthur-lab/gene_lists/blob/master/lists/NEGv1_subset_universe.tsv) |
| **DD.ID.genes.csv.gz** | Genes associated with Developmental Delay and Intelectual Disability from the Deciphering Developmental Disease studies | [Gene2Phenotype](https://www.ebi.ac.uk/gene2phenotype/) and [Gene2Phenotype donwloads](https://www.ebi.ac.uk/gene2phenotype/downloads/) |
| **HI.genes.tsv** | Haplosinsufficient genes | [HI MaCarthur Gene List](https://github.com/macarthur-lab/gene_lists/blob/master/lists/clingen_level3_genes_2018_09_13.tsv) |
| **Hom.LoF_tolerant.genes.tsv** | Homozogyous loss-of-function (LoF) tolerant genes | [Hom. LoF Tolerant MaCarthur Gene List](https://github.com/macarthur-lab/gene_lists/blob/master/lists/homozygous_lof_tolerant_twohit.tsv) |
| **OR.genes.tsv** | Olfactory Receptor genes | [OR MaCarthur Gene List](https://github.com/macarthur-lab/gene_lists/blob/master/lists/olfactory_receptors.tsv) |



## Constitutive and Cassette exons:

**GRCh38** *Constitutive* and *Cassete* exon definitions were downloaded from the [HEXEvent](http://hexevent.mmg.uci.edu/cgi-bin/HEXEvent/HEXEventWEB.cgi) web portal and stored in this repo. 

| Data File Name | Description | Data Host Link |
| -------------- | ----------- | -------------- |
| HEXEvent.GRCh38.cassette_exons.txt | Cassette exons defined with a constitutive level (constitLevel) and inclusion level (incLevel) < 1.0. (Score ranges between 0.0 - 1.0, with a 1.0 == a constitutively expressed exon and a score < 1.0 a cassette expressed exon) | [HEXEvent](http://hexevent.mmg.uci.edu/cgi-bin/HEXEvent/HEXEventWEB.cgi) | 
| HEXEvent.GRCh38.consitutive_exons.txt | Constitutive exons defined with a constitutive level (constitLevel) and inclusion level (incLevel) == 1.0. (Score ranges between 0.0 - 1.0, with a 1.0 score == a constitutively expressed exon and a score < 1.0 a cassette expressed exon) | [HEXEvent](http://hexevent.mmg.uci.edu/cgi-bin/HEXEvent/HEXEventWEB.cgi) | 
| constitutive_and_cassette_exons.gtf | A gtf file filtered for constitutive and cassette exons using the [constitutive_cassette_exon_gtf_filter.py](https://github.com/mikecormier/ConSplice-manuscript/scripts/constitutive_cassette_exon_gtf_filter.py) script in the script directory of this repo | [ConSplice Scripts](https://github.com/mikecormier/ConSplice-manuscript/scripts) |





