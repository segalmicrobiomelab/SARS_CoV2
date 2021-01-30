*Metagenome Analysis*
## Running the Analysis
* [Bacterial](https://github.com/segalmicrobiomelab/SARS_CoV2/tree/main/Metatranscriptome/Bacterial)
* [DNA Virome](https://github.com/segalmicrobiomelab/SARS_CoV2/tree/main/Metatranscriptome/DNA_Virome)
* [Phages](https://github.com/segalmicrobiomelab/SARS_CoV2/tree/main/Metatranscriptome/Phages)
* [Fungal](https://github.com/segalmicrobiomelab/SARS_CoV2/tree/main/Metatranscriptome/Fungal)
* [Resistome](https://github.com/segalmicrobiomelab/SARS_CoV2/tree/main/Metatranscriptome/Resistome)
* [Functional](https://github.com/segalmicrobiomelab/SARS_CoV2/tree/main/Metatranscriptome/Functional)

## For the Bacterial Analysis
The Required Input Data
* Microbe Raw Abundance (taxonomy_bacteria_metatranscriptome_G.txt):An abundance table with microbial species in rows and samples in columns.
* Microbe Relative Abundance (relab_taxonomy_bacteria_metatranscriptome_S.txt):A relative abundance table with microbial species in rows and samples in columns.
* MetaData (201104.WGS.UPDATE.txt): A spreadsheet with samples in rows and metadata in columns

1. Open R
1. Follow Code COVID19_Functional_Trial.r 
  1. At the start of the code you will need to load packages required for analysis, for example:
```
library(vegan)
```
* If you don't have this package you will need to install it, for example:
```
install.pakcages("vegan")
```

## For the DNA Virome Analysis
The Required Input Data
* Microbe Raw Abundance #1 (taxonomy_DNAvertebratevirus_metagenome_S.txt):An abundance table with microbial species in rows and samples in columns.
* Microbe Raw Abundance #2 (counts_taxonomy_DNAvertebratevirus_metagenome_extra_S.txt):An abundance table with microbial species in rows and samples in columns.
* Microbe Relative Abundance (relab_taxonomy_DNAvertebratevirus_metagenome_S.txt):A relative abundance table with microbial species in rows and samples in columns.
* MetaData (201104.DNA.UPDATE.txt): A spreadsheet with samples in rows and metadata in columns

1. Open R
1. Follow Code COVID19.DESEq.Function_trial_3.r 
  1. At the start of the code you will need to load packages required for analysis, for example:
```
library(vegan)
```
* If you don't have this package you will need to install it, for example:
```
install.pakcages("vegan")
```

## For the Phage Analysis
The Required Input Data
* Microbe Raw Abundance #1 (taxonomy_DNAbacteriaarchaeavirus_metagenome_S.txt):An abundance table with microbial species in rows and samples in columns.
* Microbe Raw Abundance #2 (counts_taxonomy_DNAbacteriaarchaeavirus_metagenome_extra_S.txt):An abundance table with microbial species in rows and samples in columns.
* Microbe Relative Abundance (relab_taxonomy_DNAbacteriaarchaeavirus_metagenome_S.txt):A relative abundance table with microbial species in rows and samples in columns.
* MetaData (201104.DNA.UPDATE.txt): A spreadsheet with samples in rows and metadata in columns

1. Open R
1. Follow Code COVID19.DESEq.Function_trial_3.r 
  1. At the start of the code you will need to load packages required for analysis, for example:
```
library(vegan)
```
* If you don't have this package you will need to install it, for example:
```
install.pakcages("vegan")
```
## For the Fungal Analysis
The Required Input Data
* Microbe Raw Abundance #1 (counts_taxonomy_fungi_metagenome_S.txt):An abundance table with microbial species in rows and samples in columns.
* Microbe Raw Abundance #2 (counts_taxonomy_fungi_metagenome_extra_S.txt):An abundance table with microbial species in rows and samples in columns.
* Microbe Relative Abundance (relab_taxonomy_fungi_metagenome_S.txt):A relative abundance table with microbial species in rows and samples in columns.
* MetaData (201104.WGS.UPDATE.txt): A spreadsheet with samples in rows and metadata in columns

1. Open R
1. Follow Code COVID19.DESEq.Function_trial_3.r 
  1. At the start of the code you will need to load packages required for analysis, for example:
```
library(vegan)
```
* If you don't have this package you will need to install it, for example:
```
install.pakcages("vegan")
```
## For the Resistome Analysis
The Required Input Data
* Microbe Raw Abundance (simplifiedcounts_megares_metagenome.txt):An abundance table with microbial species in rows and samples in columns.
* MetaData (201211.WGS.txt): A spreadsheet with samples in rows and metadata in columns

1. Open R
1. Follow Code COVID19.DESEq.Function_trial.r 
  1. At the start of the code you will need to load packages required for analysis, for example:
```
library(vegan)
```
* If you don't have this package you will need to install it, for example:
```
install.pakcages("vegan")
```
## For the Functional Analysis
The Required Input Data
* Microbe Raw Abundance (KO_Pathway_broken_Metagenome_sum_path03.txt):An abundance table with microbial species in rows and samples in columns.
* MetaData (201104.WGS.UPDATE): A spreadsheet with samples in rows and metadata in columns

1. Open R
1. Follow Code COVID19.DESEq.Function_trial.r 
  1. At the start of the code you will need to load packages required for analysis, for example:
```
library(vegan)
```
* If you don't have this package you will need to install it, for example:
```
install.pakcages("vegan")
```