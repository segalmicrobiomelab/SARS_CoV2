*Metatranscriptome Analysis*
## Running the Analysis
* [Bacterial](https://github.com/segalmicrobiomelab/SARS_CoV2/tree/main/Metatranscriptome/Bacterial)
* [DNA Virome](https://github.com/segalmicrobiomelab/SARS_CoV2/tree/main/Metatranscriptome/DNA_Virome)
* [RNA Virome](https://github.com/segalmicrobiomelab/SARS_CoV2/tree/main/Metatranscriptome/RNA_Virome)
* [Phages](https://github.com/segalmicrobiomelab/SARS_CoV2/tree/main/Metatranscriptome/Phages)
* [Fungal](https://github.com/segalmicrobiomelab/SARS_CoV2/tree/main/Metatranscriptome/Fungal)
* [Resistome](https://github.com/segalmicrobiomelab/SARS_CoV2/tree/main/Metatranscriptome/Resistome)
* [Functional](https://github.com/segalmicrobiomelab/SARS_CoV2/tree/main/Metatranscriptome/Functional)

## For the Bacterial Analysis
The Required Input Data
* Microbe Raw Abundance (taxonomy_bacteria_metatranscriptome_G.txt):An abundance table with microbial species in rows and samples in columns.
* Microbe Relative Abundance (relab_taxonomy_bacteria_metatranscriptome_S.txt):A relative abundance table with microbial species in rows and samples in columns.
* MetaData (201211.RNA.lns.txt): A spreadsheet with samples in rows and metadata in columns

1. Open R
1. Follow Code Bronchiectasis.Neutrophil.Elastase.r 
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
* Microbe Raw Abundance (taxonomy_DNAvertebratevirus_metatranscriptome_S.sarscov2.txt):An abundance table with microbial species in rows and samples in columns.
* Microbe Relative Abundance (relab_taxonomy_DNAvertebratevirus_metatranscriptome_S.txt):A relative abundance table with microbial species in rows and samples in columns.
* MetaData (201211.RNA.lns.txt): A spreadsheet with samples in rows and metadata in columns

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
## For the RNA Virome Analysis
The Required Input Data
* Microbe Raw Abundance (taxonomy_RNAvertebratevirus_metatranscriptome_S.sarscov2.txt):An abundance table with microbial species in rows and samples in columns.
* Microbe Relative Abundance (relab_taxonomy_RNAvertebratevirus_metatranscriptome_S.txt):A relative abundance table with microbial species in rows and samples in columns.
* MetaData (201211.RNA.lns.txt): A spreadsheet with samples in rows and metadata in columns

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
* Microbe Raw Abundance (taxonomy_DNAbacteriaarchaeavirus_metatranscriptome_S.sarscov2.txt):An abundance table with microbial species in rows and samples in columns.
* Microbe Relative Abundance (relab_taxonomy_DNAbacteriaarchaeavirus_metatranscriptome_S.txt):A relative abundance table with microbial species in rows and samples in columns.
* MetaData (201211.RNA.lns.txt): A spreadsheet with samples in rows and metadata in columns

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
* Microbe Raw Abundance (counts_taxonomy_fungi_metatranscriptome_S.txt):An abundance table with microbial species in rows and samples in columns.
* Microbe Relative Abundance (relab_taxonomy_fungi_metatranscriptome_S.txt):A relative abundance table with microbial species in rows and samples in columns.
* MetaData (201211.RNA.lns.txt): A spreadsheet with samples in rows and metadata in columns

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
## For the Resistome Analysis
The Required Input Data
* Microbe Raw Abundance (simplifiedcounts_megares_metatranscriptome.txt):An abundance table with microbial species in rows and samples in columns.
* MetaData (201211.RNA.lns.txt): A spreadsheet with samples in rows and metadata in columns

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
* Microbe Raw Abundance (KO_Pathway_broken_Metatranscriptome_sum_path03.txt):An abundance table with microbial species in rows and samples in columns.
* MetaData (201211.RNA.lns.txt): A spreadsheet with samples in rows and metadata in columns

1. Open R
1. Follow Code COVID19.DESEq.Function_trial.r 
  1. At the start of the code you will need to load packages required for analysis, for example:
```
library(vegan)
```
* If you don't have this package you will need to install it, for example:
```
install.pakcages("vegan")