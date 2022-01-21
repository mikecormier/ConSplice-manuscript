# Constitutive and Cassette exon scripts

Constitutive and Cassette exon definitions were downloaded from [HEXEvent](http://hexevent.mmg.uci.edu/cgi-bin/HEXEvent/HEXEventWEB.cgi)
  - A constitutive exon is defined as having a constitutive level (constLevel) and an inclusion level (incLevel) of 100% (== 1.0)
  - A cassette exon is defined as having a constitutive level (constLevel) and an inclusion level (incLevel) of  less than 100% (< 1.0 [i.e. 0.0 - 0.99])

The following scripts are used to create a gtf file of constittutive and cassette exon features, and score those exon features with the ConSplice regional scores. 

The scoring script will score a region around each exon to give the upstream and downstream ConSplice scores along with the ConSplice scores at the exon. 

The scripts and data associated with these scripts are used in **Figure 3** of the ConSplice manuscript.



## Script definition 

| Script | Description |
| ------ | ----------- |
| constitutive_cassette_exon_gtf_filter.py | A script to filter a GRCh38 gtf file to those exon features designated as *constitutive* and *cassette* exons as defined by [HEXEvent](http://hexevent.mmg.uci.edu/cgi-bin/HEXEvent/HEXEventWEB.cgi). The gtf file created from this script will be used to check ConSplice scores around exons as see in **Figure 3** | 
| score_constitutive_cassette_exons.py  | A script to score the gtf file created using the *constitutive_cassette_exon_gtf_filter.py* script with ConSplice scores. The output file from this script is used to see the pattern of constraint around constittutive and cassette exons as see in **Figure 3** | 



## Running the scripts

Prior to running these scripts the following need to be set up:

  - Python 3 is installed and active. 
  - conda is installed. 
  - The [ConSplice](https://github.com/mikecormier/ConSplice) CLI is installed via conda.
  - The HEXEvent constitutive and cassette exons files have been downloaded from HEXEvent. (These files can be found in the [data](https://github.com/mikecormier/ConSplice-manuscript/tree/main/data) directory of this repo under the names `HEXEvent.GRCh38.cassette_exons.txt` and `HEXEvent.GRCh38.consitutive_exons.txt`)
  - The data recipes in describe in the [ConSplice CLI Data Recipes](https://github.com/mikecormier/ConSplice/tree/main/data_recipes) directory are installed in the current environment 


 1) Run constitutive_cassette_exon_gtf_filter.py

 ```
 ## Get the ggd 'grch38-canonical-transcript-features-gencode-v1' data file 
 ### This assumes that the `grch38-canonical-transcript-features-gencode-v1` data packages has been installed in the current environment 
 canonical_gtf=$(ggd get-files grch38-canonical-transcript-features-gencode-v1 -p '*.gtf.gz')

 ## For now, we assume the path to the "data" directory of this repo with the constitutive and cassette exons data files is stored as $data

 ## We will also assume the data recipes defined in the ConSplice CLI repo are stored at $data

 constitutive_cassette_exon_gtf_filter.py \
    --gtf $canonical_gtf  \
    --constitutive $data/HEXEvent.GRCh38.consitutive_exons.txt \
    --cassette $data/HEXEvent.GRCh38.cassette_exons.txt \
    --alt-gene-symbol $data/alt_gene_names.txt \
    --out-file constitutive_and_cassette_exons.gtf 

 ```

 This script will create a gtf file named **constitutive_and_cassette_exons.gtf** 


 2) Run score_constitutive_cassette_exons.py
 ```
 ## Get the ggd 'grch38-canonical-transcript-features-gencode-v1' data file 
 ### This assumes that the `grch38-canonical-transcript-features-gencode-v1` data packages has been installed in the current environment 
 canonical_gtf=$(ggd get-files grch38-canonical-transcript-features-gencode-v1 -p '*.gtf.gz')

 ## For now, we assume the path to the "data" directory of this repo with the constitutive and cassette exons data files is stored as $data

 ## We will also assume the data recipes defined in the ConSplice CLI repo are stored at $data

 ## The '--exon-gtf' argument uses the gtf file created from the 'constitutive_cassette_exon_gtf_filter.py' in step 1

 ## We will assume the the ConSplice score file is in the $data dir

 ## We will assume a  region size of 50

 ## We will assume the percentile score name is "ConSplice_Percentile"

 python score_constitutive_cassette_exons.py \
    --gtf $canonical_gtf \
    --exon-gtf constitutive_and_cassette_exons.gtf \
    --out-file scored_constitutive_and_cassette_exons.txt \
    --constraint-scores $data/<ConSplice_Score_file> \
    --region-size 50 \
    --score-col "ConSplice_Percentile"

 ```

 This script will create a txt file name **scored_constitutive_and_cassette_exons.txt** 



