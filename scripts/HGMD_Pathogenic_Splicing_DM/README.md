## HGMD Pathogenic Splice Altering variants 

We used the professional version of the [Human Gene Mutation Database (HGMD)](http://www.hgmd.cf.ac.uk/) to download the alternative splicing disease causing mutations (DM). 

Because of HGMD's licensing, we are not allowed to share this dataset publicly. The user downloaded this dataset will need the professional license to HGMD.

Here we share the SQL query used to download the Splice Altering DM variants used in the ConSplice manuscript. 

These variants are used train and test the ConSpliceML model as seen in **Figure 5**, as well as to train the final ConSpliceML model for use to score variants.

### Downloading disease causing mutations from HGMD

If you have an HGMD license you can use the following command to download the pathogenic splice-altering variants used in the ConSplice manuscript.


Variants were downloaded using the following MySQL command from the MySQL import of HGMD 2021q1 SQL dump: 

```
SELECT * FROM splice INNER JOIN hgmd_hg38_vcf ON acc_num = hgmd_hg38_vcf.id
```
