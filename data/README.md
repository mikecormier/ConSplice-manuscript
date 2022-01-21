# ConSplice Data Files

Data files used for the ConSplice manuscript


## Gene sets:

The script used to download these files can be found in the [scripts](https://github.com/mikecormier/ConSplice-manuscript/scripts/) directory of this repo under the name [gene_sets.sh](https://github.com/mikecormier/ConSplice-manuscript/scripts/gene_sets.sh)

| Data File Name | Description | Data Host Link |
| -------------- | ----------- | -------------- |
| [AD.genes.tsv](https://github.com/mikecormier/ConSplice-manuscript/blob/main/data/AD.genes.tsv) | Autosomal Dominant genes | [AD MaCarthur Gene List](https://github.com/macarthur-lab/gene_lists/blob/master/lists/berg_ad.tsv) |
| [AR.genes.tsv](https://github.com/mikecormier/ConSplice-manuscript/blob/main/data/AR.genes.tsv) | Autosomal Recessive genes | [AR MaCarthur Gene List](https://github.com/macarthur-lab/gene_lists/blob/master/lists/all_ar.tsv) |
| [CRISPR.Essential.genes.tsv](https://github.com/mikecormier/ConSplice-manuscript/blob/main/data/CRISPR.Essential.genes.tsv) | Genes determined to be essential in cell culture CRISPR screens | [CRISPER Ess. MaCarthur Gene List](https://github.com/macarthur-lab/gene_lists/blob/master/lists/CEGv2_subset_universe.tsv) |
| [CRISPR.Nonessential.genes.tsv](https://github.com/mikecormier/ConSplice-manuscript/blob/main/data/CRISPR.Nonessential.genes.tsv) | Genes determined to be non-essential in cell culture CRISPR screens | [CRISPER NonEss. MaCarthur Gene List](https://github.com/macarthur-lab/gene_lists/blob/master/lists/NEGv1_subset_universe.tsv) |
| [DD.ID.genes.csv.gz](https://github.com/mikecormier/ConSplice-manuscript/blob/main/data/DD.ID.genes.csv.gz) | Genes associated with Developmental Delay and Intelectual Disability from the Deciphering Developmental Disease studies | [Gene2Phenotype](https://www.ebi.ac.uk/gene2phenotype/) and [Gene2Phenotype donwloads](https://www.ebi.ac.uk/gene2phenotype/downloads/) |
| [HI.genes.tsv](https://github.com/mikecormier/ConSplice-manuscript/blob/main/data/HI.genes.tsv) | Haplosinsufficient genes | [HI MaCarthur Gene List](https://github.com/macarthur-lab/gene_lists/blob/master/lists/clingen_level3_genes_2018_09_13.tsv) |
| [Hom.LoF_tolerant.genes.tsv](https://github.com/mikecormier/ConSplice-manuscript/blob/main/data/Hom.LoF_tolerant.genes.tsv) | Homozogyous loss-of-function (LoF) tolerant genes currated by the gnomAD group in their 2020 Nature paper entitled "*The mutational constraint spectrum quantified from variation in 141,456 humans*" | Supplementary Data Table 7 from [Karczewski et al.](https://www.nature.com/articles/s41586-020-2308-7) |
| [OR.genes.tsv](https://github.com/mikecormier/ConSplice-manuscript/blob/main/data/OR.genes.tsv) | Olfactory Receptor genes | [OR MaCarthur Gene List](https://github.com/macarthur-lab/gene_lists/blob/master/lists/olfactory_receptors.tsv) |



## Constitutive and Cassette exons:

**GRCh38** *Constitutive* and *Cassete* exon definitions were downloaded from the [HEXEvent](http://hexevent.mmg.uci.edu/cgi-bin/HEXEvent/HEXEventWEB.cgi) web portal and stored in this repo. 

| Data File Name | Description | Data Host Link |
| -------------- | ----------- | -------------- |
| **HEXEvent.GRCh38.cassette_exons.txt** | Cassette exons defined with a constitutive level (constitLevel) and inclusion level (incLevel) < 1.0. (Score ranges between 0.0 - 1.0, with a 1.0 == a constitutively expressed exon and a score < 1.0 a cassette expressed exon) | [HEXEvent](http://hexevent.mmg.uci.edu/cgi-bin/HEXEvent/HEXEventWEB.cgi) | 
| **HEXEvent.GRCh38.consitutive_exons.txt** | Constitutive exons defined with a constitutive level (constitLevel) and inclusion level (incLevel) == 1.0. (Score ranges between 0.0 - 1.0, with a 1.0 score == a constitutively expressed exon and a score < 1.0 a cassette expressed exon) | [HEXEvent](http://hexevent.mmg.uci.edu/cgi-bin/HEXEvent/HEXEventWEB.cgi) | 
| **constitutive_and_cassette_exons.gtf** | A gtf file filtered for constitutive and cassette exons using the [constitutive_cassette_exon_gtf_filter.py](https://github.com/mikecormier/ConSplice-manuscript/blob/main/scripts/constitutive_cassette_exons/constitutive_cassette_exon_gtf_filter.py) script in the script directory of this repo | [ConSplice Scripts](https://github.com/mikecormier/ConSplice-manuscript/tree/main/scripts/constitutive_cassette_exons) |



## Dosage Sensitivity (V<sup>G</sup>)

V<sup>G</sup> values from the [Mohammadi et al.](https://www.science.org/doi/10.1126/science.aay0256?) paper generated using GTEx data. 

| Data File Name | Description | Data Host Link |
| -------------- | ----------- | -------------- |
| [Vg.gene_level.ANEVA_tableS1.weighted_harmonic_mean.txt](https://github.com/mikecormier/ConSplice-manuscript/blob/main/data/Vg.gene_level.ANEVA_tableS1.weighted_harmonic_mean.txt) | The gene level averaged harmonic mean across tissues V<sup>G</sup> scores from Mohammadi et al. representing expression variation by gene. | Supplemental Table S1 from [Mohammadi et al.](https://www.science.org/doi/10.1126/science.aay0256?) | 


## Manually Currated Pathogenic Splice Altering Variants

A set of 376 pathogenic splice altering variants manually currated from literature with functional evidence for each variants effect on splicing.

This set excludes any splice altering variant at the canonical acceptor or donor sites, (A-2, A-1, D+1, D+2), only including variants likely to be overlooked by an analyst during rare disease clinical diagnosis 

| Data File Name | Description | Data Host Link |
| -------------- | ----------- | -------------- |
| [scored_patho_vars.from.lit.txt](https://github.com/mikecormier/ConSplice-manuscript/blob/main/data/scored_patho_vars.from.lit.txt) | The set of manually currated pathogenic splice altering variants used to compare performance of ConSplice to other splicing prediction methods. Genomic coordinates for each variant are for hg38/GRCh38 | [ConSplice](https://github.com/mikecormier/ConSplice-manuscript/tree/main/data) (this repo) |


## Benign Variants

A set of 48,978 variants that either do not impact splicing or are benign splice altering variants. The set is a combination of de novo mutations from 3-generation pedigrees in the Utah CEPH cohort, de novo mutation from 2-generation pedigrees in the Iceland deCODE cohort, and benign alternative splicing variants from GTEx. 

The variants represent a filtered set of variants that are in protein-coding genes, in the autosome, with a CADD score, SpliceAI score, SQUIRLS score, and a ConSplice score. 

| Data File Name | Description | Data Host Link |
| -------------- | ----------- | -------------- |
| [benign.combined.txt](https://github.com/mikecormier/ConSplice-manuscript/blob/main/data/benign.combined.txt) | The set of benign variants to train and test ConSpliceML | [CEPH from Sasani et al.](https://elifesciences.org/articles/46922), [deCODE from Jónsson et al.](https://www.nature.com/articles/nature24018), and [GTEx from Jaganathan et al.](https://www.sciencedirect.com/science/article/pii/S0092867418316295?via%3Dihub) |

## pLI and LOEUF

The pLI and LOEUF scores used to compared constraint metrics against the genic splicing constraint metric. 

The original pLI and LOEUF files were filtered to the constraint score column and the gene column using the `pLI` or `oe_lof_upper_bin` column names from the pLI or LOEUF file, respectively.  

| Data File Name | Description | Data Host Link |
| -------------- | ----------- | -------------- |
| [gene_to_pLI.txt](https://github.com/mikecormier/ConSplice-manuscript/blob/main/data/gene_to_pLI.txt) | By gene ExAC pLI score | [website](https://gnomad.broadinstitute.org/downloads#exac-constraint), [file](https://gnomad-public-us-east-1.s3.amazonaws.com/legacy/exac_browser/forweb_cleaned_exac_r03_march16_z_data_pLI_CNV-final.txt.gz) |
| [gene_to_LOEUF.txt](https://github.com/mikecormier/ConSplice-manuscript/blob/main/data/gene_to_LOEUF.txt) | By gene gnomAD LOEUF score | [website](https://gnomad.broadinstitute.org/downloads#v2-constraint), [file](https://gnomad-public-us-east-1.s3.amazonaws.com/release/2.1.1/constraint/gnomad.v2.1.1.lof_metrics.by_gene.txt.bgz) |





