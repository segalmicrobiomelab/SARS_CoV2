# Microbial signatures associated with poor clinical outcome in the lower airways of mechanically ventilated COVID19 patients 

Matthew Chung  
2021/01/11

<!-- MarkdownTOC autolink="true" levels="1,2,3,4" -->

- [Analysis setup](#analysis-setup)
	- [Set directories](#set-directories)
	- [Create directories](#create-directories)
	- [Download references](#download-references)
- [Pre-process reads](#pre-process-reads)
	- [Trim adaptors](#trim-adaptors)
		- [Run Trimmomatic to remove adaptor sequences](#run-trimmomatic-to-remove-adaptor-sequences)
		- [Combine single-end FASTQs for each technical replicate](#combine-single-end-fastqs-for-each-technical-replicate)
	- [Remove rRNA reads from RNA samples](#remove-rrna-reads-from-rna-samples)
		- [Split all FASTQs into small subsets](#split-all-fastqs-into-small-subsets)
		- [Remove rRNA reads from subset FASTQs](#remove-rrna-reads-from-subset-fastqs)
		- [Combine rRNA subset FASTQs from SortMeRNA for each sample](#combine-rrna-subset-fastqs-from-sortmerna-for-each-sample)
		- [Remove subset FASTQ files](#remove-subset-fastq-files)
		- [Compress FASTQ files from SortMeRNA](#compress-fastq-files-from-sortmerna)
	- [Remove human-mapping reads from all samples by Bowtie2 mapping](#remove-human-mapping-reads-from-all-samples-by-bowtie2-mapping)
		- [Create Bowtie2 index of human genome](#create-bowtie2-index-of-human-genome)
		- [Align reads to human genome](#align-reads-to-human-genome)
		- [Create FASTQ of reads that do not map to human genome](#create-fastq-of-reads-that-do-not-map-to-human-genome)
	- [Calculate read statistics from each technical replicate](#calculate-read-statistics-from-each-technical-replicate)
	- [Prepare pooled FASTQ files for analysis](#prepare-pooled-fastq-files-for-analysis)
		- [Pool technical replicates into final FASTQ files](#pool-technical-replicates-into-final-fastq-files)
		- [Compress final FASTQ files](#compress-final-fastq-files)
- [Taxonomically classify and analyze metagenomic and metatranscriptomic sequences](#taxonomically-classify-and-analyze-metagenomic-and-metatranscriptomic-sequences)
	- [Find taxonomic profiles for each sample using Kraken2](#find-taxonomic-profiles-for-each-sample-using-kraken2)
		- [Download or build Kraken databases](#download-or-build-kraken-databases)
			- [Human, Bacteria, Archaea](#human-bacteria-archaea)
			- [Virus](#virus)
			- [Fungi](#fungi)
		- [Run Kraken](#run-kraken)
			- [Human and Bacteria](#human-and-bacteria)
			- [Viruses](#viruses)
			- [Fungi](#fungi-1)
	- [Construct data frames of raw read counts for different taxa](#construct-data-frames-of-raw-read-counts-for-different-taxa)
		- [Bacteria](#bacteria)
			- [Set R inputs](#set-r-inputs)
			- [Load packages and view sessionInfo](#load-packages-and-view-sessioninfo)
			- [Set list of taxonomy ranks of interest](#set-list-of-taxonomy-ranks-of-interest)
			- [Read in Kraken2 outputs for metagenome and metatranscriptome samples](#read-in-kraken2-outputs-for-metagenome-and-metatranscriptome-samples)
			- [Exclude metagenome samples with low read counts](#exclude-metagenome-samples-with-low-read-counts)
			- [Compile list of detected taxonomic identifiers for metagenome and metatranscriptome samples](#compile-list-of-detected-taxonomic-identifiers-for-metagenome-and-metatranscriptome-samples)
			- [Construct data frame of the raw reads for each taxonomic identifier](#construct-data-frame-of-the-raw-reads-for-each-taxonomic-identifier)
			- [Output data frames containing the raw reads for each taxonomic identifier](#output-data-frames-containing-the-raw-reads-for-each-taxonomic-identifier)
		- [Virus](#virus-1)
			- [Set R inputs](#set-r-inputs-1)
			- [Load packages and view sessionInfo](#load-packages-and-view-sessioninfo-1)
			- [Set list of taxonomy ranks of interest](#set-list-of-taxonomy-ranks-of-interest-1)
			- [Read in Kraken2 outputs for metagenome and metatranscriptome samples](#read-in-kraken2-outputs-for-metagenome-and-metatranscriptome-samples-1)
			- [Exclude metagenome samples with low read counts](#exclude-metagenome-samples-with-low-read-counts-1)
			- [Compile list of detected taxonomic identifiers for metagenome and metatranscriptome samples](#compile-list-of-detected-taxonomic-identifiers-for-metagenome-and-metatranscriptome-samples-1)
			- [Construct data frame of the raw reads for each taxonomic identifier](#construct-data-frame-of-the-raw-reads-for-each-taxonomic-identifier-1)
			- [Output data frames containing the raw reads for each taxonomic identifier](#output-data-frames-containing-the-raw-reads-for-each-taxonomic-identifier-1)
		- [Fungi](#fungi-2)
			- [Set R inputs](#set-r-inputs-2)
			- [Load packages and view sessionInfo](#load-packages-and-view-sessioninfo-2)
			- [Set list of taxonomy ranks of interest](#set-list-of-taxonomy-ranks-of-interest-2)
			- [Read in Kraken2 outputs for metagenome and metatranscriptome samples](#read-in-kraken2-outputs-for-metagenome-and-metatranscriptome-samples-2)
			- [Exclude metagenome samples with low read counts](#exclude-metagenome-samples-with-low-read-counts-2)
			- [Compile list of detected taxonomic identifiers for metagenome and metatranscriptome samples](#compile-list-of-detected-taxonomic-identifiers-for-metagenome-and-metatranscriptome-samples-2)
			- [Construct data frame of the raw reads for each taxonomic identifier](#construct-data-frame-of-the-raw-reads-for-each-taxonomic-identifier-2)
			- [Output data frames containing the raw reads for each taxonomic identifier](#output-data-frames-containing-the-raw-reads-for-each-taxonomic-identifier-2)
	- [Create genera and species lists for DNA and RNA viruses with vertebrate, bacteria/archaea, and other hosts](#create-genera-and-species-lists-for-dna-and-rna-viruses-with-vertebrate-bacteriaarchaea-and-other-hosts)
		- [Set R inputs](#set-r-inputs-3)
		- [Load packages and view sessionInfo](#load-packages-and-view-sessioninfo-3)
		- [Read in list of mammalian virus taxonomy families](#read-in-list-of-mammalian-virus-taxonomy-families)
		- [Read in list of kraken reports](#read-in-list-of-kraken-reports)
		- [Parse all kraken outputs to map mammalian viral genera and species taxa to their families](#parse-all-kraken-outputs-to-map-mammalian-viral-genera-and-species-taxa-to-their-families)
		- [Manually remove duplicate genera entries incorrectly binned](#manually-remove-duplicate-genera-entries-incorrectly-binned)
		- [Output a list containing all mammalian virus genera and species split into DNA, RNA viruses and host categories](#output-a-list-containing-all-mammalian-virus-genera-and-species-split-into-dna-rna-viruses-and-host-categories)
	- [Combine taxonomy raw counts data frames into one data frame](#combine-taxonomy-raw-counts-data-frames-into-one-data-frame)
		- [Genera](#genera)
		- [Species](#species)
	- [Add species categorizations to counts data frames and create combined relative abundance dataframes](#add-species-categorizations-to-counts-data-frames-and-create-combined-relative-abundance-dataframes)
		- [Set R inputs](#set-r-inputs-4)
		- [Load packages and view sessionInfo](#load-packages-and-view-sessioninfo-4)
		- [Compile list of taxonomic counts files](#compile-list-of-taxonomic-counts-files)
		- [Compile taxa lists](#compile-taxa-lists)
		- [Manually curate categorizations of viruses not categorized by Kraken2](#manually-curate-categorizations-of-viruses-not-categorized-by-kraken2)
		- [Set function for calculating relative abundance](#set-function-for-calculating-relative-abundance)
		- [Set list of samples to exclude due to IRB issues](#set-list-of-samples-to-exclude-due-to-irb-issues)
		- [Write counts and relative abundance tables split by taxa categories](#write-counts-and-relative-abundance-tables-split-by-taxa-categories)
		- [Write combined counts and relative abundance table](#write-combined-counts-and-relative-abundance-table)
	- [Plot box plot of most abundant OTUs](#plot-box-plot-of-most-abundant-otus)
		- [Set R inputs](#set-r-inputs-5)
		- [Load packages and view sessionInfo](#load-packages-and-view-sessioninfo-5)
		- [Read in taxonomic relative abundance files](#read-in-taxonomic-relative-abundance-files)
		- [Read in meta-data maps](#read-in-meta-data-maps)
		- [Read in curated list of vertebrate viruses](#read-in-curated-list-of-vertebrate-viruses)
		- [Set function for creating box plots](#set-function-for-creating-box-plots)
		- [Ranks all taxa in each of the different datasets](#ranks-all-taxa-in-each-of-the-different-datasets)
		- [Writes the top 50 taxa for BAL and UA samples separately](#writes-the-top-50-taxa-for-bal-and-ua-samples-separately)
		- [Create ordering data frames](#create-ordering-data-frames)
		- [Plot OTU relative abundance box plots](#plot-otu-relative-abundance-box-plots)
- [Deplete Kraken human-assigned reads from final FASTQ files](#deplete-kraken-human-assigned-reads-from-final-fastq-files)
	- [Split Kraken outputs](#split-kraken-outputs)
	- [Split FASTQ files](#split-fastq-files)
	- [Create subset FASTQs without human-assigned reads](#create-subset-fastqs-without-human-assigned-reads)
	- [Combine human-depleted FASTQs for each sample](#combine-human-depleted-fastqs-for-each-sample)
	- [Compress combined human-depleted FASTQs](#compress-combined-human-depleted-fastqs)
- [Conduct functional profiling analysis of metagenome and metatranscriptome sequences](#conduct-functional-profiling-analysis-of-metagenome-and-metatranscriptome-sequences)
	- [Conduct BLAST searches of human-depleted reads against FMAP KO reference database](#conduct-blast-searches-of-human-depleted-reads-against-fmap-ko-reference-database)
	- [Derive count values from FMAP BLASTX searches](#derive-count-values-from-fmap-blastx-searches)
		- [Split FMAP BLASTX outputs](#split-fmap-blastx-outputs)
		- [Create abundance dataframes for KO terms](#create-abundance-dataframes-for-ko-terms)
	- [Construct KO counts data frames from FMAP output](#construct-ko-counts-data-frames-from-fmap-output)
		- [Set R inputs](#set-r-inputs-6)
		- [Load packages and view sessionInfo](#load-packages-and-view-sessioninfo-6)
		- [Read in FMAP abundance files](#read-in-fmap-abundance-files)
		- [Compile list of KO terms and identifiers quantified across all FMAP abundance files](#compile-list-of-ko-terms-and-identifiers-quantified-across-all-fmap-abundance-files)
		- [Compile counts for identified KOs for all samples](#compile-counts-for-identified-kos-for-all-samples)
		- [Identify pathway KOs for each KO term identified by FMAP](#identify-pathway-kos-for-each-ko-term-identified-by-fmap)
		- [Construct table of KO pathway counts](#construct-table-of-ko-pathway-counts)
		- [Append rownames of KO terms counts data frame with descriptors](#append-rownames-of-ko-terms-counts-data-frame-with-descriptors)
		- [Output counts data frames for KO terms and pathway KO terms](#output-counts-data-frames-for-ko-terms-and-pathway-ko-terms)
- [Conduct resistome analysis of metagenome and metatranscriptome sequences](#conduct-resistome-analysis-of-metagenome-and-metatranscriptome-sequences)
	- [Download MEGAres database of antibiotic resistance genes](#download-megares-database-of-antibiotic-resistance-genes)
	- [Quantify antibiotic resistance genes using Salmon](#quantify-antibiotic-resistance-genes-using-salmon)
		- [Construct index for reference MEGAres database](#construct-index-for-reference-megares-database)
		- [Run Salmon to quantify genes](#run-salmon-to-quantify-genes)
	- [Construct resistome count tables](#construct-resistome-count-tables)
		- [Set R inputs](#set-r-inputs-7)
		- [Load packages and view sessionInfo](#load-packages-and-view-sessioninfo-7)
		- [Set list of samples to exclude due to IRB issues](#set-list-of-samples-to-exclude-due-to-irb-issues-1)
		- [Read input count files of quantified MEGAres genes](#read-input-count-files-of-quantified-megares-genes)
		- [Filter out MEGAres genes so that only genes that actively confer resistance are kept](#filter-out-megares-genes-so-that-only-genes-that-actively-confer-resistance-are-kept)
		- [Simplify MEGAres counts dataframes](#simplify-megares-counts-dataframes)
			- [Collapse genes by orthologous gene families](#collapse-genes-by-orthologous-gene-families)
			- [Collapse genes by antibiotic classes](#collapse-genes-by-antibiotic-classes)
		- [Write output MEGAres quantified antibiotic resistance gene count tables](#write-output-megares-quantified-antibiotic-resistance-gene-count-tables)
- [Characterize CRISPR arrays in metagenome data](#characterize-crispr-arrays-in-metagenome-data)
	- [Extract bacterial reads for each metagenome sample](#extract-bacterial-reads-for-each-metagenome-sample)
		- [Split Kraken outputs](#split-kraken-outputs-1)
		- [Split FASTQ files](#split-fastq-files-1)
		- [Create subset metagenome FASTQs with only bacteria assigned reads](#create-subset-metagenome-fastqs-with-only-bacteria-assigned-reads)
		- [Combine bacteria-selected FASTQs for each metagenome sample](#combine-bacteria-selected-fastqs-for-each-metagenome-sample)
		- [Compress combined bacteria-selected FASTQs](#compress-combined-bacteria-selected-fastqs)
	- [Identify spacers and direct repeats from metagenome bacterial reads using Crass](#identify-spacers-and-direct-repeats-from-metagenome-bacterial-reads-using-crass)
		- [Run Crass](#run-crass)
		- [Combine Crass FASTA files for each sample](#combine-crass-fasta-files-for-each-sample)
	- [Conduct BLAST search on Crass spacers](#conduct-blast-search-on-crass-spacers)
		- [Bacterial DB](#bacterial-db)
		- [Viral DB](#viral-db)
	- [Conduct analysis on bacterial-viral interactions present in patient samples](#conduct-analysis-on-bacterial-viral-interactions-present-in-patient-samples)
		- [Set R inputs](#set-r-inputs-8)
		- [Load packages and view sessionInfo](#load-packages-and-view-sessioninfo-8)
		- [Parse BLASTN output files for Crass hits into a data frame](#parse-blastn-output-files-for-crass-hits-into-a-data-frame)
		- [Map taxonomic identifers to accession ids](#map-taxonomic-identifers-to-accession-ids)
		- [Read in meta-data maps](#read-in-meta-data-maps-1)
		- [Read in list of clinically relevant meta-data and remove continuous meta-data variables](#read-in-list-of-clinically-relevant-meta-data-and-remove-continuous-meta-data-variables)
		- [Remove bad viral hits](#remove-bad-viral-hits)
		- [Check for significantly over-represented bacterial-viral interactions in data set](#check-for-significantly-over-represented-bacterial-viral-interactions-in-data-set)
			- [Set function for conducting Fisher's exact test on Crass data](#set-function-for-conducting-fishers-exact-test-on-crass-data)
			- [Checking for significantly over-represented bacteria-virus interactions](#checking-for-significantly-over-represented-bacteria-virus-interactions)
			- [Checking for significantly over-represented bacteria present in bacteria-virus interactions](#checking-for-significantly-over-represented-bacteria-present-in-bacteria-virus-interactions)
			- [Checking for significantly over-represented viruses present in bacteria-virus interactions](#checking-for-significantly-over-represented-viruses-present-in-bacteria-virus-interactions)
		- [Creating plotting data frames for bacterial and viral interactions](#creating-plotting-data-frames-for-bacterial-and-viral-interactions)
		- [Plot heatmaps summarizing the relative abundance and presence of CRISPR array data](#plot-heatmaps-summarizing-the-relative-abundance-and-presence-of-crispr-array-data)
			- [Set function for plotting a binary heatmap of the presence of Crass interactions](#set-function-for-plotting-a-binary-heatmap-of-the-presence-of-crass-interactions)
			- [Set an aesthetic map to color significant bacterial and viral members and interactions between BAL and UA samples](#set-an-aesthetic-map-to-color-significant-bacterial-and-viral-members-and-interactions-between-bal-and-ua-samples)
			- [Plot binary heatmap of Crass interactions](#plot-binary-heatmap-of-crass-interactions)

<!-- /MarkdownTOC -->

# Analysis setup

## Set directories 

```{bash, eval = F}
WORKING_DIR=/hpcdata/lpd_sg/mchung/covid_bal_metagenome_metatranscriptome/
SCRIPTS_DIR=~/scripts/
TEMP_DIR=/hpcdata/scratch/mchung

FMAP_BIN_DIR=/hpcdata/lpd_sg/mchung/covid_bal_metagenome_metatranscriptome/db/FMAP

HUMANN_DB_DIR=/hpcdata/bio_data/biobakery/
BLAST_DB_DIR=/hpcdata/bio_data/blast_db_29SEP2020/
KRAKEN2_DB_DIR=/hpcdata/bio_data/kraken_db/
KRAKEN2_FUNGI_DB_DIR=/hpcdata/lpd_sg/mchung/covid_bal_metagenome_metatranscriptome/db/kraken_fungi
KRAKEN2_VIRUS_DB_DIR=/hpcdata/lpd_sg/mchung/covid_bal_metagenome_metatranscriptome/db/kraken/
METAPHLAN_DB_DIR=/hpcdata/bio_data/biobakery/
SORTMERNA_DB_DIR=/hpcdata/bio_data/sortmerna/data/rRNA_databases/
THREADS=4
```

## Create directories

```{bash, eval = F}
mkdir -p "$WORKING_DIR"/final_reads/metagenome
mkdir -p "$WORKING_DIR"/final_reads/metatranscriptome
mkdir -p "$WORKING_DIR"/reads
mkdir -p "$WORKING_DIR"/references

mkdir -p "$TEMP_DIR"/metagenome
mkdir -p "$TEMP_DIR"/metatranscriptome
```

## Download references
```{bash, eval = F}
wget https://urldefense.proofpoint.com/v2/url?u=ftp-3A__ftp.ncbi.nlm.nih.gov_genomes_all_GCF_009_858_895_GCF-5F009858895.2-5FASM985889v3_GCF-5F009858895.2-5FASM985889v3-5Fgenomic.fna.gz&d=DwIGAg&c=j5oPpO0eBH1iio48DtsedeElZfc04rx3ExJHeIIZuCs&r=ufTTecG4vkmes8mF_HId8gR_hbrb8-dKDFW3-voxs2A&m=9Ldy7Wd7XZ7NQzW7DY5909D3ElJvgX4saq7j5jlDI94&s=hIAbRFaq8l61FyWOCWVSQdUWuww-43niQUDcwAWvonQ&e=  -O "$WORKING_DIR"/references/covid19.fna.gz

wget https://urldefense.proofpoint.com/v2/url?u=https-3A__ftp.ncbi.nlm.nih.gov_genomes_all_GCF_000_013_425_GCF-5F000013425.1-5FASM1342v1_GCF-5F000013425.1-5FASM1342v1-5Fgenomic.fna.gz&d=DwIGAg&c=j5oPpO0eBH1iio48DtsedeElZfc04rx3ExJHeIIZuCs&r=ufTTecG4vkmes8mF_HId8gR_hbrb8-dKDFW3-voxs2A&m=9Ldy7Wd7XZ7NQzW7DY5909D3ElJvgX4saq7j5jlDI94&s=9dQH-OxiHE9T4kEReFusTbB-XRwIU4pBojq7POXynwk&e=  -O "$WORKING_DIR"/references/saureus.fna.gz
```

# Pre-process reads

## Trim adaptors

### Run Trimmomatic to remove adaptor sequences
```{bash, eval = F}
module load trimmomatic/0.36

find "$WORKING_DIR"/reads -name "*fastq.gz" | sed "s/_R[12]_001//g" | sed "s/[.].*//g" | sort -n | uniq -c | awk '$1 == 2 {print $NF}' | sort -n > "$WORKING_DIR"/trimmomatic.list

echo -e 'SAMPLE="$(sed -ne "$SGE_TASK_ID p" "$WORKING_DIR"/trimmomatic.list)" \n
java -jar "$EBROOTTRIMMOMATIC"/trimmomatic-0.36.jar PE \
	"$SAMPLE"_R1_001.fastq.gz "$SAMPLE"_R2_001.fastq.gz \
	"$SAMPLE"_R1_001.trimmed.fastq.gz "$SAMPLE"_R1_001.se.trimmed.fastq.gz \
	"$SAMPLE"_R2_001.trimmed.fastq.gz "$SAMPLE"_R2_001.se.trimmed.fastq.gz \
	ILLUMINACLIP:"$EBROOTTRIMMOMATIC"/adapters/TruSeq3-PE-2.fa:2:30:10:2:keepBothReads \
	LEADING:3 \
	TRAILING:3 \
	MINLEN:36' | qsub -v WORKING_DIR="$WORKING_DIR" -N trimmomatic -wd "$WORKING_DIR"/stderrout -t 1-"$(wc -l "$WORKING_DIR"/trimmomatic.list | awk '{print $1}')"
```

### Combine single-end FASTQs for each technical replicate
```{bash, eval = F}
find "$WORKING_DIR"/reads -name "*fastq.gz" | sed "s/_R[12]_001//g" | sed "s/[.].*//g" | sort -n | uniq -c | awk '$1 == 6 {print $NF}' | sort -n > "$WORKING_DIR"/reformat.list

echo -e 'SAMPLE="$(sed -ne "$SGE_TASK_ID p" "$WORKING_DIR"/reformat.list)" \n
	zcat "$SAMPLE"_R1_001.se.trimmed.fastq.gz "$SAMPLE"_R2_001.se.trimmed.fastq.gz | gzip > "$SAMPLE".se.trimmed.fastq.gz' | qsub -v WORKING_DIR="$WORKING_DIR" -N reformat -t 1-"$(wc -l "$WORKING_DIR"/reformat.list | awk '{print $1}')"
```


## Remove rRNA reads from RNA samples

### Split all FASTQs into small subsets
```{bash, eval = F}
find "$WORKING_DIR"/reads -name "*fastq.gz" | grep -v Plate | sed "s/_R[12]_001//g" | sed "s/[.].*//g" | sort -n | uniq -c | awk '$1 == 7 {print $NF}' | sort -n > "$WORKING_DIR"/sortmerna_split.list

echo -e 'SAMPLE="$(sed -ne "$SGE_TASK_ID p" "$WORKING_DIR"/sortmerna_split.list)" \n
	zcat "$SAMPLE"_R1_001.trimmed.fastq.gz | split -l 4000000 -d - "$SAMPLE"_R1_001.trimmed.fastq' | qsub -v WORKING_DIR="$WORKING_DIR" -N split -t 1-"$(wc -l "$WORKING_DIR"/sortmerna_split.list | awk '{print $1}')"
echo -e 'SAMPLE="$(sed -ne "$SGE_TASK_ID p" "$WORKING_DIR"/sortmerna_split.list)" \n
	zcat "$SAMPLE"_R2_001.trimmed.fastq.gz | split -l 4000000 -d - "$SAMPLE"_R2_001.trimmed.fastq' | qsub -v WORKING_DIR="$WORKING_DIR" -N split -t 1-"$(wc -l "$WORKING_DIR"/sortmerna_split.list | awk '{print $1}')"
echo -e 'SAMPLE="$(sed -ne "$SGE_TASK_ID p" "$WORKING_DIR"/sortmerna_split.list)" \n
	zcat "$SAMPLE".se.trimmed.fastq.gz | split -l 4000000 -d - "$SAMPLE".se.trimmed.fastq' | qsub -v WORKING_DIR="$WORKING_DIR" -N split -t 1-"$(wc -l "$WORKING_DIR"/sortmerna_split.list | awk '{print $1}')"
```

### Remove rRNA reads from subset FASTQs

##### Paired-end reads

```{bash, eval = F}
module load sortmerna/4.2.0

find "$WORKING_DIR"/reads -name "*fastq[0-9]*" | grep R1 | sort -n > "$WORKING_DIR"/sortmerna_pe.list

echo -e 'SAMPLE="$(sed -ne "$SGE_TASK_ID p" "$WORKING_DIR"/sortmerna_pe.list)" \n
sortmerna \
	--ref "$SORTMERNA_DB_DIR"/rfam-5s-database-id98.fasta \
	--ref "$SORTMERNA_DB_DIR"/rfam-5.8s-database-id98.fasta \
	--ref "$SORTMERNA_DB_DIR"/silva-arc-16s-id95.fasta \
	--ref "$SORTMERNA_DB_DIR"/silva-arc-23s-id98.fasta \
	--ref "$SORTMERNA_DB_DIR"/silva-euk-28s-id98.fasta \
	--ref "$SORTMERNA_DB_DIR"/silva-bac-23s-id98.fasta \
	--ref "$SORTMERNA_DB_DIR"/silva-euk-18s-id95.fasta \
	--ref "$SORTMERNA_DB_DIR"/silva-bac-16s-id90.fasta \
	--reads "$SAMPLE" \
	--reads "$(echo $SAMPLE | sed "s/R1/R2/g")" \
	--aligned "$SAMPLE".pe.trimmed.rRNA \
	--other "$SAMPLE".pe.trimmed.non-rRNA \
	--workdir "$TEMP_DIR"/$(basename $SAMPLE) \
	-m 100000 \
	--threads "$THREADS" \
	--fastx --paired_in --otu_map --out2' |  qsub -pe threaded "$THREADS" -l mem_free=100G -v SORTMERNA_DB_DIR="$SORTMERNA_DB_DIR",THREADS="$THREADS",TEMP_DIR="$TEMP_DIR",WORKING_DIR="$WORKING_DIR" -N sortmerna_pe -t 1-"$(wc -l "$WORKING_DIR"/sortmerna_pe.list | awk '{print $1}')"
```

##### Single-end reads

```{bash, eval = F}
module load sortmerna/4.2.0

find "$WORKING_DIR"/reads/ -name "*fastq[0-9]*" | grep se | sort -n > "$WORKING_DIR"/sortmerna_se.list

echo -e 'SAMPLE="$(sed -ne "$SGE_TASK_ID p" "$WORKING_DIR"/sortmerna_se.list)" \n
	sortmerna \
	--ref "$SORTMERNA_DB_DIR"/rfam-5s-database-id98.fasta \
	--ref "$SORTMERNA_DB_DIR"/rfam-5.8s-database-id98.fasta \
	--ref "$SORTMERNA_DB_DIR"/silva-arc-16s-id95.fasta \
	--ref "$SORTMERNA_DB_DIR"/silva-arc-23s-id98.fasta \
	--ref "$SORTMERNA_DB_DIR"/silva-euk-28s-id98.fasta \
	--ref "$SORTMERNA_DB_DIR"/silva-bac-23s-id98.fasta \
	--ref "$SORTMERNA_DB_DIR"/silva-euk-18s-id95.fasta \
	--ref "$SORTMERNA_DB_DIR"/silva-bac-16s-id90.fasta \
	--reads "$SAMPLE" \
	--aligned "$SAMPLE".se.trimmed.rRNA \
	--other "$SAMPLE".se.trimmed.non-rRNA \
	--workdir "$TEMP_DIR"/$(basename $SAMPLE) \
	-m 100000 \
	--threads "$THREADS" \
	--fastx --otu_map' |  qsub -pe threaded "$THREADS" -l mem_free=100G -v SORTMERNA_DB_DIR="$SORTMERNA_DB_DIR",THREADS="$THREADS",TEMP_DIR="$TEMP_DIR",WORKING_DIR="$WORKING_DIR" -N sortmerna_se -t 1-"$(wc -l "$WORKING_DIR"/sortmerna_se.list | awk '{print $1}')"
```

### Combine rRNA subset FASTQs from SortMeRNA for each sample

##### Paired-end reads

```{bash, eval = F}
for SAMPLE in $(find $PWD -name "*fastq[0-9]*" | grep -v "[.]se[.]" | grep rRNA | sed "s/.trimmed.fastq.*//g" | sort -n | uniq)
do
	for i in $(ls "$SAMPLE"* | grep -v "[.]se[.]" | grep [.]rRNA_ | grep -v otu | grep fwd | sort -n)
	do
		cat $i >> $(echo $SAMPLE | sed "s/_R1_001//g").pe.trimmed.rRNA_fwd.fastq
	done
	for i in $(ls "$SAMPLE"* | grep -v "[.]se[.]" | grep [.]rRNA_ | grep -v otu | grep rev | sort -n)
	do
		cat $i >> $(echo $SAMPLE | sed "s/_R1_001//g").pe.trimmed.rRNA_rev.fastq
	done
	for i in $(ls "$SAMPLE"* | grep -v "[.]se[.]" | grep [.]non-rRNA_ | grep -v otu | grep fwd | sort -n)
	do
		cat $i >> $(echo $SAMPLE | sed "s/_R1_001//g").pe.trimmed.non-rRNA_fwd.fastq
	done
	for i in $(ls "$SAMPLE"* | grep -v "[.]se[.]" | grep [.]non-rRNA_ | grep -v otu | grep rev | sort -n)
	do
		cat $i >> $(echo $SAMPLE | sed "s/_R1_001//g").pe.trimmed.non-rRNA_rev.fastq
	done
done
```

##### Single-end reads

```{bash, eval = F}
for SAMPLE in $(find $PWD -name "*fastq[0-9]*" | grep "[.]se[.]" | grep rRNA | sed "s/.trimmed.fastq.*//g" | sort -n | uniq)
do
	for i in $(ls "$SAMPLE"* | grep [.]rRNA | grep -v otu | grep -v log | sort -n)
	do
		cat $i >> "$SAMPLE".trimmed.rRNA.fastq
	done
	for i in $(ls "$SAMPLE"* | grep [.]non-rRNA | grep -v otu | grep -v log | sort -n)
	do
		cat $i >> "$SAMPLE".trimmed.non-rRNA.fastq
	done
done
```

### Remove subset FASTQ files

```{bash, eval = F}
rm *fastq[0-9]*
```

### Compress FASTQ files from SortMeRNA

```{bash, eval = F}
find "$WORKING_DIR"/reads -name "*fastq" > "$WORKING_DIR"/gzip.list

echo -e 'FASTQ="$(sed -ne "$SGE_TASK_ID p" "$WORKING_DIR"/gzip.list)" \n
gzip $FASTQ' | qsub -v WORKING_DIR="$WORKING_DIR" -N gzip -wd "$WORKING_DIR"/stderrout/ -t 1-"$(wc -l "$WORKING_DIR"/gzip.list | awk '{print $1}')"

```

## Remove human-mapping reads from all samples by Bowtie2 mapping

### Create Bowtie2 index of human genome

```{bash, eval = F}
module load bowtie2/2.3.4.1

bowtie2-build "$WORKING_DIR"/references/human.fna "$WORKING_DIR"/references/human.fna
```

### Align reads to human genome

##### Metagenome

```{bash, eval = F}
module load bowtie2/2.3.4.1
module load samtools/1.9-goolf-1.7.20

find "$WORKING_DIR"/reads -name "*fastq.gz" | grep Plate | sed "s/_R[12]_001//g" | sed "s/[.].*//g" | sort -n | uniq -c | awk '$1 == 7 {print $NF}' | sort -n > "$WORKING_DIR"/bowtie2.list

echo -e 'SAMPLE="$(sed -ne "$SGE_TASK_ID p" "$WORKING_DIR"/bowtie2.list)" \n
bowtie2 \
-x "$WORKING_DIR"/references/human.fna \
-1 "$SAMPLE"_R1_001.trimmed.fastq.gz \
-2 "$SAMPLE"_R2_001.trimmed.fastq.gz \
-U "$SAMPLE".se.trimmed.fastq.gz \
-p "$THREADS" | samtools sort -o "$SAMPLE".trimmed.human.bam' | qsub -pe threaded "$THREADS" -v THREADS="$THREADS",WORKING_DIR="$WORKING_DIR" -N bowtie2 -wd "$WORKING_DIR"/stderrout/ -t 1-"$(wc -l "$WORKING_DIR"/bowtie2.list | awk '{print $1}')"
```

##### Metatranscriptome

```{bash, eval = F}
module load bowtie2/2.3.4.1
module load samtools/1.9-goolf-1.7.20

find "$WORKING_DIR"/reads -name "*fastq.gz" | grep -v Plate | sed "s/_R[12]_001//g" | sed "s/[.].*//g" | sort -n | uniq -c | awk '$1 == 13 {print $NF}' | sort -n > "$WORKING_DIR"/bowtie2.list

echo -e 'SAMPLE="$(sed -ne "$SGE_TASK_ID p" "$WORKING_DIR"/bowtie2.list)" \n
bowtie2 \
-x "$WORKING_DIR"/references/human.fna \
-1 "$SAMPLE".pe.trimmed.non-rRNA_fwd.fastq.gz \
-2 "$SAMPLE".pe.trimmed.non-rRNA_rev.fastq.gz \
-U "$SAMPLE".se.trimmed.non-rRNA.fastq.gz \
-p "$THREADS" | samtools sort -o "$SAMPLE".trimmed.non-rRNA.human.bam' | qsub -pe threaded "$THREADS" -v THREADS="$THREADS",WORKING_DIR="$WORKING_DIR" -N bowtie2 -wd "$WORKING_DIR"/stderrout/ -t 1-"$(wc -l "$WORKING_DIR"/bowtie2.list | awk '{print $1}')"
```

### Create FASTQ of reads that do not map to human genome

##### Metagenome
```{bash, eval = F}
module load samtools/1.9-goolf-1.7.20

echo -e 'SAMPLE="$(sed -ne "$SGE_TASK_ID p" "$WORKING_DIR"/bowtie2.list)" \n
samtools fastq -f 77 "$SAMPLE".trimmed.human.bam | gzip > "$SAMPLE".trimmed.nohuman.1.fastq.gz' | qsub -v WORKING_DIR="$WORKING_DIR" -N samtools-fastq -t 1-"$(wc -l "$WORKING_DIR"/bowtie2.list| awk '{print $1}')"
echo -e 'SAMPLE="$(sed -ne "$SGE_TASK_ID p" "$WORKING_DIR"/bowtie2.list)" \n
samtools fastq -f 141 "$SAMPLE".trimmed.human.bam | gzip  > "$SAMPLE".trimmed.nohuman.2.fastq.gz' | qsub -v WORKING_DIR="$WORKING_DIR" -N samtools-fastq -t 1-"$(wc -l "$WORKING_DIR"/bowtie2.list| awk '{print $1}')"
echo -e 'SAMPLE="$(sed -ne "$SGE_TASK_ID p" "$WORKING_DIR"/bowtie2.list)" \n
samtools fastq -f 4 -F 1 "$SAMPLE".trimmed.human.bam | gzip  > "$SAMPLE".trimmed.nohuman.se.fastq.gz' | qsub -v WORKING_DIR="$WORKING_DIR" -N samtools-fastq -t 1-"$(wc -l "$WORKING_DIR"/bowtie2.list | awk '{print $1}')"
```

##### Metatranscriptome
```{bash, eval = F}
module load samtools/1.9-goolf-1.7.20

echo -e 'SAMPLE="$(sed -ne "$SGE_TASK_ID p" "$WORKING_DIR"/bowtie2.list)" \n
samtools fastq -f 77 "$SAMPLE".trimmed.non-rRNA.human.bam | gzip > "$SAMPLE".trimmed.non-rRNA.nohuman.1.fastq.gz' | qsub -v WORKING_DIR="$WORKING_DIR" -N samtools-fastq -t 1-"$(wc -l "$WORKING_DIR"/bowtie2.list| awk '{print $1}')"
echo -e 'SAMPLE="$(sed -ne "$SGE_TASK_ID p" "$WORKING_DIR"/bowtie2.list)" \n
samtools fastq -f 141 "$SAMPLE".trimmed.non-rRNA.human.bam | gzip  > "$SAMPLE".trimmed.non-rRNA.nohuman.2.fastq.gz' | qsub -v WORKING_DIR="$WORKING_DIR" -N samtools-fastq -t 1-"$(wc -l "$WORKING_DIR"/bowtie2.list| awk '{print $1}')"
echo -e 'SAMPLE="$(sed -ne "$SGE_TASK_ID p" "$WORKING_DIR"/bowtie2.list)" \n
samtools fastq -f 4 -F 1 "$SAMPLE".trimmed.non-rRNA.human.bam | gzip  > "$SAMPLE".trimmed.non-rRNA.nohuman.se.fastq.gz' | qsub -v WORKING_DIR="$WORKING_DIR" -N samtools-fastq -t 1-"$(wc -l "$WORKING_DIR"/bowtie2.list | awk '{print $1}')"
```

## Calculate read statistics from each technical replicate

##### Metagenome

Column 1: Replicate ID
Column 2: Reads in R1
Column 3: Reads in R2
Column 4: Reads in trimmed R1
Column 5: Reads in trimmed R2
Column 6: Reads in trimmed singletons
Column 7: Reads in trimmed, human-depleted R1
Column 8: Reads in trimmed, human-depleted R2
Column 9: Reads in trimmed, human-depleted singletons

```{bash, eval = F}
find "$WORKING_DIR"/reads/2020-11-23 -name "*fastq.gz" | grep Plate | sed "s/_R[12]_001//g" | sed "s/[.].*//g" | sort -n | uniq > "$WORKING_DIR"/metagenome_read_stats.list

echo -e 'SAMPLE="$(sed -ne "$SGE_TASK_ID p" "$WORKING_DIR"/metagenome_read_stats.list)" \n
	echo -e ""$SAMPLE"\t"$(expr $(zcat "$SAMPLE"_R1_001.fastq.gz | wc -l) / 4)"\t"$(expr $(zcat "$SAMPLE"_R2_001.fastq.gz | wc -l) / 4)"\t"$(expr $(zcat "$SAMPLE"_R1_001.trimmed.fastq.gz | wc -l) / 4)"\t"$(expr $(zcat "$SAMPLE"_R2_001.trimmed.fastq.gz | wc -l) / 4)"\t"$(expr $(zcat "$SAMPLE".se.trimmed.fastq.gz | wc -l) / 4)"\t"$(expr $(zcat "$SAMPLE".trimmed.nohuman.1.fastq.gz | wc -l) / 4)"\t"$(expr $(zcat "$SAMPLE".trimmed.nohuman.2.fastq.gz | wc -l) / 4)"\t"$(expr $(zcat "$SAMPLE".trimmed.nohuman.se.fastq.gz | wc -l) / 4)"" >> "$WORKING_DIR"/stats/metagenome_stats.tsv' | qsub -v WORKING_DIR="$WORKING_DIR" -N metagenome_stats -t 1-"$(wc -l "$WORKING_DIR"/metagenome_read_stats.list | awk '{print $1}')"
```

##### Metatranscriptome

Column 1: Replicate ID
Column 2: Reads in R1
Column 3: Reads in R2
Column 4: Reads in trimmed R1
Column 5: Reads in trimmed R2
Column 6: Reads in trimmed singletons
Column 7: Reads in trimmed, rRNA-depleted R1
Column 8: Reads in trimmed, rRNA-depleted R2
Column 9: Reads in trimmed, rRNA-depleted singletons
Column 10: Reads in trimmed, rRNA-depleted, human-depleted R1
Column 11: Reads in trimmed, rRNA-depleted, human-depleted R2
Column 12: Reads in trimmed, rRNA-depleted, human-depleted singletons

```{bash, eval = F}
find "$WORKING_DIR"/reads -name "*fastq.gz" | grep -v Plate | sed "s/_R[12]_001//g" | sed "s/[.].*//g" | sort -n | uniq > "$WORKING_DIR"/metatranscriptome_read_stats.list

echo -e 'SAMPLE="$(sed -ne "$SGE_TASK_ID p" "$WORKING_DIR"/metatranscriptome_read_stats.list)" \n
	echo -e ""$SAMPLE"\t"$(expr $(zcat "$SAMPLE"_R1_001.fastq.gz | wc -l) / 4)"\t"$(expr $(zcat "$SAMPLE"_R2_001.fastq.gz | wc -l) / 4)"\t"$(expr $(zcat "$SAMPLE"_R1_001.trimmed.fastq.gz | wc -l) / 4)"\t"$(expr $(zcat "$SAMPLE"_R2_001.trimmed.fastq.gz | wc -l) / 4)"\t"$(expr $(zcat "$SAMPLE".se.trimmed.fastq.gz | wc -l) / 4)"\t"$(expr $(zcat "$SAMPLE".pe.trimmed.non-rRNA_fwd.fastq.gz | wc -l) / 4)"\t"$(expr $(zcat "$SAMPLE".pe.trimmed.non-rRNA_rev.fastq.gz | wc -l) / 4)"\t"$(expr $(zcat "$SAMPLE".se.trimmed.non-rRNA.fastq.gz | wc -l) / 4)"\t"$(expr $(zcat "$SAMPLE".trimmed.non-rRNA.nohuman.1.fastq.gz | wc -l) / 4)"\t"$(expr $(zcat "$SAMPLE".trimmed.non-rRNA.nohuman.2.fastq.gz | wc -l) / 4)"\t"$(expr $(zcat "$SAMPLE".trimmed.non-rRNA.nohuman.se.fastq.gz | wc -l) / 4)"" >> "$WORKING_DIR"/stats/metatranscriptome_stats.tsv' | qsub -v WORKING_DIR="$WORKING_DIR" -N metatranscriptome_stats -t 1-"$(wc -l "$WORKING_DIR"/metatranscriptome_read_stats.list | awk '{print $1}')"
```

## Prepare pooled FASTQ files for analysis

### Pool technical replicates into final FASTQ files

##### Metagenome

```{bash, eval = F}
cat "$WORKING_DIR"/maps/COVID.WGS.Map.a1.txt | tail -n+2 | cut -f1 | while read SAMPLE
do
	rm "$WORKING_DIR"/final_reads/metagenome/"$(echo $SAMPLE | sed "s/[.]/_/g")".fastq
	rm "$WORKING_DIR"/final_reads/metagenome/"$(echo $SAMPLE | sed "s/[.]/_/g")".1.fastq
	rm "$WORKING_DIR"/final_reads/metagenome/"$(echo $SAMPLE | sed "s/[.]/_/g")".2.fastq
	rm "$WORKING_DIR"/final_reads/metagenome/"$(echo $SAMPLE | sed "s/[.]/_/g")".se.fastq

for i in $(find $WORKING_DIR | grep Plate | grep "$SAMPLE" | grep "nohuman" | grep "fastq.gz" | sort -n)
do
zcat "$i" >> "$WORKING_DIR"/final_reads/metagenome/"$(echo $SAMPLE | sed "s/[.]/_/g")".fastq
done

	for i in $(find $WORKING_DIR/reads | grep Plate | grep "$SAMPLE" | grep "nohuman" | grep "[.]1[.]fastq.gz" | sort -n)
	do
		zcat "$i" >> "$WORKING_DIR"/final_reads/metagenome/"$(echo $SAMPLE | sed "s/[.]/_/g")".1.fastq
	done
	for i in $(find $WORKING_DIR/reads | grep Plate | grep "$SAMPLE" | grep "nohuman" | grep "[.]2[.]fastq.gz" | sort -n)
	do
		zcat "$i" >> "$WORKING_DIR"/final_reads/metagenome/"$(echo $SAMPLE | sed "s/[.]/_/g")".2.fastq
	done
	for i in $(find $WORKING_DIR/reads | grep Plate | grep "$SAMPLE" | grep "nohuman" | grep "[.]se[.]fastq.gz" | sort -n)
	do
		zcat "$i" >> "$WORKING_DIR"/final_reads/metagenome/"$(echo $SAMPLE | sed "s/[.]/_/g")".se.fastq
	done
done
```

##### Metatranscriptome

```{bash, eval = F}
cat "$WORKING_DIR"/maps/COVID.RNA.Map.a1.txt | tail -n+2 | cut -f1 | while read SAMPLE
do
	rm "$WORKING_DIR"/final_reads/metatranscriptome/"$(echo $SAMPLE | sed "s/[.]/_/g")".fastq
	rm "$WORKING_DIR"/final_reads/metatranscriptome/"$(echo $SAMPLE | sed "s/[.]/_/g")".1.fastq
	rm "$WORKING_DIR"/final_reads/metatranscriptome/"$(echo $SAMPLE | sed "s/[.]/_/g")".2.fastq
	rm "$WORKING_DIR"/final_reads/metatranscriptome/"$(echo $SAMPLE | sed "s/[.]/_/g")".se.fastq

	for i in $(find $WORKING_DIR/reads | grep -v Plate | grep "$SAMPLE" | grep "nohuman" | grep "[.]1[.]fastq.gz" | sort -n)
	do
		zcat "$i" >> "$WORKING_DIR"/final_reads/metatranscriptome/"$(echo $SAMPLE | sed "s/[.]/_/g")".1.fastq
	done
	for i in $(find $WORKING_DIR/reads | grep -v Plate | grep "$SAMPLE" | grep "nohuman" | grep "[.]2[.]fastq.gz" | sort -n)
	do
		zcat "$i" >> "$WORKING_DIR"/final_reads/metatranscriptome/"$(echo $SAMPLE | sed "s/[.]/_/g")".2.fastq
	done
	for i in $(find $WORKING_DIR/reads | grep -v Plate | grep "$SAMPLE" | grep "nohuman" | grep "[.]se[.]fastq.gz" | sort -n)
	do
		zcat "$i" >> "$WORKING_DIR"/final_reads/metatranscriptome/"$(echo $SAMPLE | sed "s/[.]/_/g")".se.fastq
	done
	for i in $(find $WORKING_DIR | grep -v Plate | grep -v 2020-05-29 | grep "$SAMPLE" | grep "nohuman" | grep "fastq.gz" | sort -n)
	do
		zcat "$i" >> "$WORKING_DIR"/final_reads/metatranscriptome/"$(echo $SAMPLE | sed "s/[.]/_/g")".fastq
	done
done
```

### Compress final FASTQ files

```{bash, eval = F}
find "$WORKING_DIR"/final_reads -name "*fastq" > "$WORKING_DIR"/gzip.list

echo -e 'FASTQ="$(sed -ne "$SGE_TASK_ID p" "$WORKING_DIR"/gzip.list)" \n
gzip $FASTQ' | qsub -WORKING_DIR="$WORKING_DIR" -N gzip -wd "$WORKING_DIR"/stderrout/ -t 1-"$(wc -l "$WORKING_DIR"/gzip.list | awk '{print $1}')"
```

# Taxonomically classify and analyze metagenomic and metatranscriptomic sequences

## Find taxonomic profiles for each sample using Kraken2

### Download or build Kraken databases

#### Human, Bacteria, Archaea

##### Download and build database

```{bash, eval = F}
module load kraken/2.0.7-foss-2016b-Perl-5.24.0
module load blast+/2.9.0

kraken2-build --standard --db "$KRAKEN2_DB_DIR"
```

#### Virus

##### Download and build database

```{bash, eval = F}
module load kraken/2.0.7-foss-2016b-Perl-5.24.0
module load blast+/2.9.0

kraken2-build --download-library virus --db "$KRAKEN2_VIRUS_DB_DIR"
kraken2-build --build --db "$KRAKEN2_VIRUS_DB_DIR" --threads "$THREADS"
```

#### Fungi

##### Download database

```{bash, eval = F}
module load kraken/2.0.7-foss-2016b-Perl-5.24.0
module load blast+/2.9.0

kraken2-build --download-library fungi --db "$KRAKEN2_FUNGI_DB_DIR"
```

##### Download and add Candida auris sequences to database

```{bash, eval = F}
wget -P "$KRAKEN2_FUNGI_DB_DIR" https://urldefense.proofpoint.com/v2/url?u=ftp-3A__ftp.ncbi.nlm.nih.gov_genomes_all_GCA_003_013_715_GCA-5F003013715.2-5FASM301371v2_GCA-5F003013715.2-5FASM301371v2-5Fgenomic.fna.gz&d=DwIGAg&c=j5oPpO0eBH1iio48DtsedeElZfc04rx3ExJHeIIZuCs&r=ufTTecG4vkmes8mF_HId8gR_hbrb8-dKDFW3-voxs2A&m=9Ldy7Wd7XZ7NQzW7DY5909D3ElJvgX4saq7j5jlDI94&s=hg6onEQCfKjcBr623aYGQK_72Zu-YJlUeXQ4YDpDizk&e= 
wget -P "$KRAKEN2_FUNGI_DB_DIR" https://urldefense.proofpoint.com/v2/url?u=ftp-3A__ftp.ncbi.nlm.nih.gov_genomes_all_GCA_008_275_145_GCA-5F008275145.1-5FASM827514v1_GCA-5F008275145.1-5FASM827514v1-5Fgenomic.fna.gz&d=DwIGAg&c=j5oPpO0eBH1iio48DtsedeElZfc04rx3ExJHeIIZuCs&r=ufTTecG4vkmes8mF_HId8gR_hbrb8-dKDFW3-voxs2A&m=9Ldy7Wd7XZ7NQzW7DY5909D3ElJvgX4saq7j5jlDI94&s=INvRfCWc7ocA2QfV9hqShHonmLjUIPCF7D0ALJ5BJrY&e= 
gunzip "$KRAKEN2_FUNGI_DB_DIR"/GCA_003013715.2_ASM301371v2_genomic.fna.gz
gunzip "$KRAKEN2_FUNGI_DB_DIR"/GCA_008275145.1_ASM827514v1_genomic.fna.gz

sed -i "s/>/>kraken:taxid|498019|/g" "$KRAKEN2_FUNGI_DB_DIR"/GCA_003013715.2_ASM301371v2_genomic.fna
sed -i "s/>/>kraken:taxid|498019|/g" "$KRAKEN2_FUNGI_DB_DIR"/GCA_008275145.1_ASM827514v1_genomic.fna

kraken2-build --add-to-library "$KRAKEN2_FUNGI_DB_DIR"/GCA_003013715.2_ASM301371v2_genomic.fna --db "$KRAKEN2_FUNGI_DB_DIR"
kraken2-build --add-to-library "$KRAKEN2_FUNGI_DB_DIR"/GCA_008275145.1_ASM827514v1_genomic.fna --db "$KRAKEN2_FUNGI_DB_DIR"
```

##### Build database

```{bash, eval = F}
kraken2-build --build --db "$KRAKEN2_FUNGI_DB_DIR" --threads "$THREADS"
```

### Run Kraken

#### Human and Bacteria

##### Metagenome

```{bash, eval = F}
module load kraken/2.0.7-foss-2016b-Perl-5.24.0

find "$WORKING_DIR"/final_reads/metagenome/ -name "*fastq.gz" | sed "s/[.].*//g" | sort -n | uniq > "$WORKING_DIR"/kraken_metagenome.list

echo -e 'SAMPLE="$(sed -ne "$SGE_TASK_ID p" "$WORKING_DIR"/kraken_metagenome.list)" \n
kraken2 --gzip-compressed  \
--db "$KRAKEN2_DB_DIR" \
--threads "$THREADS" \
--classified-out "$WORKING_DIR"/kraken/metagenome/"$(basename "$SAMPLE")".classified.fastq \
--unclassified-out "$WORKING_DIR"/kraken/metagenome/"$(basename "$SAMPLE")".unclassified.fastq \
--report "$WORKING_DIR"/kraken/metagenome/"$(basename "$SAMPLE")".kraken_report.txt \
--output "$WORKING_DIR"/kraken/metagenome/"$(basename "$SAMPLE")".kraken_output.txt \
"$SAMPLE".fastq.gz' | qsub -pe threaded "$THREADS" -l mem_free=100G -v KRAKEN2_DB_DIR="$KRAKEN2_DB_DIR",THREADS="$THREADS",WORKING_DIR="$WORKING_DIR" -N kraken2 -t 1-"$(wc -l "$WORKING_DIR"/kraken_metagenome.list | awk '{print $1}')"
```

##### Metatranscriptome

```{bash, eval = F}
module load kraken/2.0.7-foss-2016b-Perl-5.24.0

find "$WORKING_DIR"/final_reads/metatranscriptome/ -name "*fastq.gz" | sed "s/[.].*//g" | sort -n | uniq > "$WORKING_DIR"/kraken_metatranscriptome.list

echo -e 'SAMPLE="$(sed -ne "$SGE_TASK_ID p" "$WORKING_DIR"/kraken_metatranscriptome.list)" \n
kraken2 --gzip-compressed  \
--db "$KRAKEN2_DB_DIR" \
--threads "$THREADS" \
--classified-out "$WORKING_DIR"/kraken/metatranscriptome/"$(basename "$SAMPLE")".classified.fastq \
--unclassified-out "$WORKING_DIR"/kraken/metatranscriptome/"$(basename "$SAMPLE")".unclassified.fastq \
--report "$WORKING_DIR"/kraken/metatranscriptome/"$(basename "$SAMPLE")".kraken_report.txt \
--output "$WORKING_DIR"/kraken/metatranscriptome/"$(basename "$SAMPLE")".kraken_output.txt \
"$SAMPLE".fastq.gz' | qsub -pe threaded "$THREADS" -l mem_free=100G -v KRAKEN2_DB_DIR="$KRAKEN2_DB_DIR",THREADS="$THREADS",WORKING_DIR="$WORKING_DIR" -N kraken2 -t 1-"$(wc -l "$WORKING_DIR"/kraken_metatranscriptome.list | awk '{print $1}')"
```

#### Viruses

##### Metagenome

```{bash, eval = F}
module load kraken/2.0.7-foss-2016b-Perl-5.24.0

find "$WORKING_DIR"/final_reads/metagenome/ -name "*fastq.gz" | sed "s/[.].*//g" | sort -n | uniq > "$WORKING_DIR"/kraken_metagenome.list

echo -e 'SAMPLE="$(sed -ne "$SGE_TASK_ID p" "$WORKING_DIR"/kraken_metagenome.list)" \n
kraken2 --gzip-compressed  \
--db "$KRAKEN2_VIRUS_DB_DIR" \
--threads "$THREADS" \
--classified-out "$WORKING_DIR"/kraken/metagenome/"$(basename "$SAMPLE")".virus.classified.fastq \
--unclassified-out "$WORKING_DIR"/kraken/metagenome/"$(basename "$SAMPLE")".virus.unclassified.fastq \
--report "$WORKING_DIR"/kraken/metagenome/"$(basename "$SAMPLE")".virus.kraken_report.txt \
--output "$WORKING_DIR"/kraken/metagenome/"$(basename "$SAMPLE")".virus.kraken_output.txt \
"$SAMPLE".fastq.gz' | qsub -pe threaded "$THREADS" -l mem_free=100G -v KRAKEN2_VIRUS_DB_DIR="$KRAKEN2_VIRUS_DB_DIR",THREADS="$THREADS",WORKING_DIR="$WORKING_DIR" -wd "$WORKING_DIR"/stderrout -N kraken2 -t 1-"$(wc -l "$WORKING_DIR"/kraken_metagenome.list | awk '{print $1}')"
```

##### Metatranscriptome

```{bash, eval = F}
module load kraken/2.0.7-foss-2016b-Perl-5.24.0

find "$WORKING_DIR"/final_reads/metatranscriptome/ -name "*fastq.gz" | sed "s/[.].*//g" | sort -n | uniq > "$WORKING_DIR"/kraken_metatranscriptome.list

echo -e 'SAMPLE="$(sed -ne "$SGE_TASK_ID p" "$WORKING_DIR"/kraken_metatranscriptome.list)" \n
kraken2 --gzip-compressed  \
--db "$KRAKEN2_VIRUS_DB_DIR" \
--threads "$THREADS" \
--classified-out "$WORKING_DIR"/kraken/metatranscriptome/"$(basename "$SAMPLE")".virus.classified.fastq \
--unclassified-out "$WORKING_DIR"/kraken/metatranscriptome/"$(basename "$SAMPLE")".virus.unclassified.fastq \
--report "$WORKING_DIR"/kraken/metatranscriptome/"$(basename "$SAMPLE")".virus.kraken_report.txt \
--output "$WORKING_DIR"/kraken/metatranscriptome/"$(basename "$SAMPLE")".virus.kraken_output.txt \
"$SAMPLE".fastq.gz' | qsub -pe threaded "$THREADS" -l mem_free=100G -v KRAKEN2_VIRUS_DB_DIR="$KRAKEN2_VIRUS_DB_DIR",THREADS="$THREADS",WORKING_DIR="$WORKING_DIR" -wd "$WORKING_DIR"/stderrout -N kraken2 -t 1-"$(wc -l "$WORKING_DIR"/kraken_metatranscriptome.list | awk '{print $1}')"
```

#### Fungi

##### Metagenome

```{bash, eval = F}
module load kraken/2.0.7-foss-2016b-Perl-5.24.0

find "$WORKING_DIR"/final_reads/metagenome/ -name "*fastq.gz" | sed "s/[.].*//g" | sort -n | uniq > "$WORKING_DIR"/kraken_metagenome.list

echo -e 'SAMPLE="$(sed -ne "$SGE_TASK_ID p" "$WORKING_DIR"/kraken_metagenome.list)" \n
kraken2 --gzip-compressed  \
--db "$KRAKEN2_FUNGI_DB_DIR" \
--threads "$THREADS" \
--classified-out "$WORKING_DIR"/kraken/metagenome/"$(basename "$SAMPLE")".fungi.classified.fastq \
--unclassified-out "$WORKING_DIR"/kraken/metagenome/"$(basename "$SAMPLE")".fungi.unclassified.fastq \
--report "$WORKING_DIR"/kraken/metagenome/"$(basename "$SAMPLE")".fungi.kraken_report.txt \
--output "$WORKING_DIR"/kraken/metagenome/"$(basename "$SAMPLE")".fungi.kraken_output.txt \
"$SAMPLE".fastq.gz' | qsub -pe threaded "$THREADS" -l mem_free=100G -v KRAKEN2_FUNGI_DB_DIR="$KRAKEN2_FUNGI_DB_DIR",THREADS="$THREADS",WORKING_DIR="$WORKING_DIR" -wd "$WORKING_DIR"/stderrout -N kraken2 -t 1-"$(wc -l "$WORKING_DIR"/kraken_metagenome.list | awk '{print $1}')"
```

##### Metatranscriptome

```{bash, eval = F}
module load kraken/2.0.7-foss-2016b-Perl-5.24.0

find "$WORKING_DIR"/final_reads/metatranscriptome/ -name "*fastq.gz" | sed "s/[.].*//g" | sort -n | uniq > "$WORKING_DIR"/kraken_metatranscriptome.list

echo -e 'SAMPLE="$(sed -ne "$SGE_TASK_ID p" "$WORKING_DIR"/kraken_metatranscriptome.list)" \n
kraken2 --gzip-compressed  \
--db "$KRAKEN2_FUNGI_DB_DIR" \
--threads "$THREADS" \
--classified-out "$WORKING_DIR"/kraken/metatranscriptome/"$(basename "$SAMPLE")".fungi.classified.fastq \
--unclassified-out "$WORKING_DIR"/kraken/metatranscriptome/"$(basename "$SAMPLE")".fungi.unclassified.fastq \
--report "$WORKING_DIR"/kraken/metatranscriptome/"$(basename "$SAMPLE")".fungi.kraken_report.txt \
--output "$WORKING_DIR"/kraken/metatranscriptome/"$(basename "$SAMPLE")".fungi.kraken_output.txt \
"$SAMPLE".fastq.gz' | qsub -pe threaded "$THREADS" -l mem_free=100G -v KRAKEN2_FUNGI_DB_DIR="$KRAKEN2_FUNGI_DB_DIR",THREADS="$THREADS",WORKING_DIR="$WORKING_DIR" -wd "$WORKING_DIR"/stderrout -N kraken2 -t 1-"$(wc -l "$WORKING_DIR"/kraken_metatranscriptome.list | awk '{print $1}')"
```

## Construct data frames of raw read counts for different taxa

### Bacteria

#### Set R inputs

```{R}
WORKING.DIR <- ""
```

#### Load packages and view sessionInfo

```{R}
sessionInfo()
```

```{R, eval = F}
R version 4.0.2 (2020-06-22)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 18363)

Matrix products: default

locale:
[1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252    LC_MONETARY=English_United States.1252
[4] LC_NUMERIC=C                           LC_TIME=English_United States.1252    

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

loaded via a namespace (and not attached):
[1] compiler_4.0.2 tools_4.0.2    knitr_1.29     xfun_0.16
```

#### Set list of taxonomy ranks of interest 

```{R}
taxonomy.list.names <- c("P","C","O","F","G","S")
```

#### Read in Kraken2 outputs for metagenome and metatranscriptome samples

```{R}
kraken_output.metagenome.list <- list.files(paste0(WORKING.DIR,"/kraken/metagenome/"),
											pattern="*kraken_report.txt")
kraken_output.metatranscriptome.list <- list.files(paste0(WORKING.DIR,"/kraken/metatranscriptome/"),
												   pattern="*kraken_report.txt")
```

#### Exclude metagenome samples with low read counts

The following DNA/metagenome samples do not have a sufficient number of reads for the analysis: 
Blank_234_1 (no sequencing data present)
Blank_234_2
UCS_0136_BKG

All RNA/metatranscriptome samples have a sufficient number of reads.

```{R}
kraken_output.metagenome.list <- grep("Blank_234_2",kraken_output.metagenome.list,value = T,invert = T)
kraken_output.metagenome.list <- grep("UCS_0136_BKG",kraken_output.metagenome.list,value = T,invert = T)
```

#### Compile list of detected taxonomic identifiers for metagenome and metatranscriptome samples

##### Metagenome

```{R}
taxonomy_bacteria.metagenome.list <- vector(mode = "list", length = length(taxonomy.list.names))

for(i in 1:length(kraken_output.metagenome.list)){
  kraken_output <- read.delim(paste0(WORKING.DIR,"/kraken/metagenome/",kraken_output.metagenome.list[i]),
                              header = F)
  kraken_output[,6] <- gsub("^ {2,}","",kraken_output[,6])
  
  d.coords <- which(kraken_output[,4] == "D")
  
  if(length(which(kraken_output[,6] == "Bacteria")) != 0){
    bacteria.coords1 <- which(kraken_output[,6] == "Bacteria")
    if(which(d.coords == bacteria.coords1) == length(d.coords)){
      bacteria.coords2 <- nrow(kraken_output)
    }else{
      bacteria.coords2 <- d.coords[which(d.coords == bacteria.coords1)+1]-1
    }
    kraken_output_bacteria <- kraken_output[bacteria.coords1:bacteria.coords2,]
  }

  for(j in 1:length(taxonomy.list.names)){
    taxonomy_bacteria.metagenome.list[[j]] <- unique(c(taxonomy_bacteria.metagenome.list[[j]],kraken_output_bacteria[kraken_output_bacteria[,4] == taxonomy.list.names[j],6]))
  }
}
```

##### Metatranscriptome

```{R}
taxonomy_bacteria.metatranscriptome.list <- vector(mode = "list", length = length(taxonomy.list.names))

for(i in 1:length(kraken_output.metatranscriptome.list)){
  kraken_output <- read.delim(paste0(WORKING.DIR,"/kraken/metatranscriptome/",kraken_output.metatranscriptome.list[i]),
                              header = F)
  kraken_output[,6] <- gsub("^ {2,}","",kraken_output[,6])
  
  d.coords <- which(kraken_output[,4] == "D")
  
  bacteria.coords1 <- which(kraken_output[,6] == "Bacteria")
  if(which(d.coords == bacteria.coords1) == length(d.coords)){
    bacteria.coords2 <- nrow(kraken_output)
  }else{
    bacteria.coords2 <- d.coords[which(d.coords == bacteria.coords1)+1]-1
  }
  
  kraken_output_bacteria <- kraken_output[bacteria.coords1:bacteria.coords2,]

  for(j in 1:length(taxonomy.list.names)){
    taxonomy_bacteria.metatranscriptome.list[[j]] <- unique(c(taxonomy_bacteria.metatranscriptome.list[[j]],kraken_output_bacteria[kraken_output_bacteria[,4] == taxonomy.list.names[j],6]))
  }
}
```

#### Construct data frame of the raw reads for each taxonomic identifier

##### Metagenome

```{R}
taxonomy_bacteria.metagenome.df <- vector(mode = "list", length = length(taxonomy.list.names))

for(i in 1:length(taxonomy.list.names)){
  taxonomy_bacteria.metagenome.df[[i]] <- as.data.frame(matrix(0,
                                                        nrow=length(taxonomy_bacteria.metagenome.list[[i]]),
                                                        ncol=length(kraken_output.metagenome.list)))
  rownames(taxonomy_bacteria.metagenome.df[[i]]) <- taxonomy_bacteria.metagenome.list[[i]]
  colnames(taxonomy_bacteria.metagenome.df[[i]]) <- gsub("[.].*","",kraken_output.metagenome.list)
  
  for(j in 1:length(kraken_output.metagenome.list)){
    kraken_output <- read.delim(paste0(WORKING.DIR,"/kraken/metagenome/",kraken_output.metagenome.list[j]),
                                header = F)
    kraken_output[,6] <- gsub("^ {2,}","",kraken_output[,6])
   
    kraken_output <- kraken_output[kraken_output[,4] == taxonomy.list.names[i],]
    
    taxonomy_bacteria.metagenome.df[[i]][,j] <- kraken_output[match(rownames(taxonomy_bacteria.metagenome.df[[i]]),kraken_output[,6]),2]
  }
  taxonomy_bacteria.metagenome.df[[i]][is.na(taxonomy_bacteria.metagenome.df[[i]])] <- 0
}
```

##### Metatranscriptome

```{R}
taxonomy_bacteria.metatranscriptome.df <- vector(mode = "list", length = length(taxonomy.list.names))

for(i in 1:length(taxonomy.list.names)){
  taxonomy_bacteria.metatranscriptome.df[[i]] <- as.data.frame(matrix(0,
                                                        nrow=length(taxonomy_bacteria.metatranscriptome.list[[i]]),
                                                        ncol=length(kraken_output.metatranscriptome.list)))
  rownames(taxonomy_bacteria.metatranscriptome.df[[i]]) <- taxonomy_bacteria.metatranscriptome.list[[i]]
  colnames(taxonomy_bacteria.metatranscriptome.df[[i]]) <- gsub("[.].*","",kraken_output.metatranscriptome.list)
  
  for(j in 1:length(kraken_output.metatranscriptome.list)){
    kraken_output <- read.delim(paste0(WORKING.DIR,"/kraken/metatranscriptome/",kraken_output.metatranscriptome.list[j]),
                                header = F)
    kraken_output[,6] <- gsub("^ {2,}","",kraken_output[,6])
    
    kraken_output <- kraken_output[kraken_output[,4] == taxonomy.list.names[i],]
    
    taxonomy_bacteria.metatranscriptome.df[[i]][,j] <- kraken_output[match(rownames(taxonomy_bacteria.metatranscriptome.df[[i]]),kraken_output[,6]),2]
  }
  taxonomy_bacteria.metatranscriptome.df[[i]][is.na(taxonomy_bacteria.metatranscriptome.df[[i]])] <- 0
}
```

#### Output data frames containing the raw reads for each taxonomic identifier

##### Metagenome

```{R}
for(i in 1:length(taxonomy_bacteria.metagenome.df)){
  write.table(taxonomy_bacteria.metagenome.df[[i]],
              paste0(WORKING.DIR,"/data_frames/taxonomy_bacteria_metagenome_",taxonomy.list.names[i],".tsv"),
              row.names = T,
              col.names = T,
              quote = F,
              sep = "\t")
}
```

##### Metatranscriptome

```{R}
for(i in 1:length(taxonomy_bacteria.metatranscriptome.df)){
  write.table(taxonomy_bacteria.metatranscriptome.df[[i]],
              paste0(WORKING.DIR,"/data_frames/taxonomy_bacteria_metatranscriptome_",taxonomy.list.names[i],".tsv"),
              row.names = T,
              col.names = T,
              quote = F,
              sep = "\t")
}
```

### Virus

#### Set R inputs

```{R}
WORKING.DIR <- ""
```

#### Load packages and view sessionInfo

```{R}
sessionInfo()
```

```{R, eval = F}
R version 4.0.2 (2020-06-22)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 18363)

Matrix products: default

locale:
[1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252    LC_MONETARY=English_United States.1252
[4] LC_NUMERIC=C                           LC_TIME=English_United States.1252    

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

loaded via a namespace (and not attached):
[1] compiler_4.0.2 tools_4.0.2    knitr_1.29     xfun_0.16
```

#### Set list of taxonomy ranks of interest 

```{R}
taxonomy.list.names <- c("P","C","O","F","G","S")
```

#### Read in Kraken2 outputs for metagenome and metatranscriptome samples

```{R}
kraken_output.metagenome.list <- list.files(paste0(WORKING.DIR,"/kraken/metagenome/"),
											pattern="*virus.kraken_report.txt")
kraken_output.metatranscriptome.list <- list.files(paste0(WORKING.DIR,"/kraken/metatranscriptome/"),
												   pattern="*virus.kraken_report.txt")
```

#### Exclude metagenome samples with low read counts

The following DNA/metagenome samples do not have a sufficient number of reads for the analysis: 
Blank_234_1 (no sequencing data present)
Blank_234_2
UCS_0136_BKG

All RNA/metatranscriptome samples have a sufficient number of reads.

```{R}
kraken_output.metagenome.list <- grep("Blank_234_2",kraken_output.metagenome.list,value = T,invert = T)
kraken_output.metagenome.list <- grep("UCS_0136_BKG",kraken_output.metagenome.list,value = T,invert = T)
```

#### Compile list of detected taxonomic identifiers for metagenome and metatranscriptome samples

##### Metagenome

```{R}
taxonomy_virus.metagenome.list <- vector(mode = "list", length = length(taxonomy.list.names))

for(i in 1:length(kraken_output.metagenome.list)){
  kraken_output <- read.delim(paste0(WORKING.DIR,"/kraken/metagenome/",kraken_output.metagenome.list[i]),
                              header = F)
  kraken_output[,6] <- gsub("^ {2,}","",kraken_output[,6])
  
  d.coords <- which(kraken_output[,4] == "D")
  
  if(length(which(kraken_output[,6] == "virus")) != 0){
    virus.coords1 <- which(kraken_output[,6] == "virus")
    if(which(d.coords == virus.coords1) == length(d.coords)){
      virus.coords2 <- nrow(kraken_output)
    }else{
      virus.coords2 <- d.coords[which(d.coords == virus.coords1)+1]-1
    }
    kraken_output_virus <- kraken_output[virus.coords1:virus.coords2,]
  }

  for(j in 1:length(taxonomy.list.names)){
    taxonomy_virus.metagenome.list[[j]] <- unique(c(taxonomy_virus.metagenome.list[[j]],kraken_output_virus[kraken_output_virus[,4] == taxonomy.list.names[j],6]))
  }
}
```

##### Metatranscriptome

```{R}
taxonomy_virus.metatranscriptome.list <- vector(mode = "list", length = length(taxonomy.list.names))

for(i in 1:length(kraken_output.metatranscriptome.list)){
  kraken_output <- read.delim(paste0(WORKING.DIR,"/kraken/metatranscriptome/",kraken_output.metatranscriptome.list[i]),
                              header = F)
  kraken_output[,6] <- gsub("^ {2,}","",kraken_output[,6])
  
  d.coords <- which(kraken_output[,4] == "D")
  
  virus.coords1 <- which(kraken_output[,6] == "virus")
  if(which(d.coords == virus.coords1) == length(d.coords)){
    virus.coords2 <- nrow(kraken_output)
  }else{
    virus.coords2 <- d.coords[which(d.coords == virus.coords1)+1]-1
  }
  
  kraken_output_virus <- kraken_output[virus.coords1:virus.coords2,]

  for(j in 1:length(taxonomy.list.names)){
    taxonomy_virus.metatranscriptome.list[[j]] <- unique(c(taxonomy_virus.metatranscriptome.list[[j]],kraken_output_virus[kraken_output_virus[,4] == taxonomy.list.names[j],6]))
  }
}
```

#### Construct data frame of the raw reads for each taxonomic identifier

##### Metagenome

```{R}
taxonomy_virus.metagenome.df <- vector(mode = "list", length = length(taxonomy.list.names))

for(i in 1:length(taxonomy.list.names)){
  taxonomy_virus.metagenome.df[[i]] <- as.data.frame(matrix(0,
                                                        nrow=length(taxonomy_virus.metagenome.list[[i]]),
                                                        ncol=length(kraken_output.metagenome.list)))
  rownames(taxonomy_virus.metagenome.df[[i]]) <- taxonomy_virus.metagenome.list[[i]]
  colnames(taxonomy_virus.metagenome.df[[i]]) <- gsub("[.].*","",kraken_output.metagenome.list)
  
  for(j in 1:length(kraken_output.metagenome.list)){
    kraken_output <- read.delim(paste0(WORKING.DIR,"/kraken/metagenome/",kraken_output.metagenome.list[j]),
                                header = F)
    kraken_output[,6] <- gsub("^ {2,}","",kraken_output[,6])
   
    kraken_output <- kraken_output[kraken_output[,4] == taxonomy.list.names[i],]
    
    taxonomy_virus.metagenome.df[[i]][,j] <- kraken_output[match(rownames(taxonomy_virus.metagenome.df[[i]]),kraken_output[,6]),2]
  }
  taxonomy_virus.metagenome.df[[i]][is.na(taxonomy_virus.metagenome.df[[i]])] <- 0
}
```

##### Metatranscriptome

```{R}
taxonomy_virus.metatranscriptome.df <- vector(mode = "list", length = length(taxonomy.list.names))

for(i in 1:length(taxonomy.list.names)){
  taxonomy_virus.metatranscriptome.df[[i]] <- as.data.frame(matrix(0,
                                                        nrow=length(taxonomy_virus.metatranscriptome.list[[i]]),
                                                        ncol=length(kraken_output.metatranscriptome.list)))
  rownames(taxonomy_virus.metatranscriptome.df[[i]]) <- taxonomy_virus.metatranscriptome.list[[i]]
  colnames(taxonomy_virus.metatranscriptome.df[[i]]) <- gsub("[.].*","",kraken_output.metatranscriptome.list)
  
  for(j in 1:length(kraken_output.metatranscriptome.list)){
    kraken_output <- read.delim(paste0(WORKING.DIR,"/kraken/metatranscriptome/",kraken_output.metatranscriptome.list[j]),
                                header = F)
    kraken_output[,6] <- gsub("^ {2,}","",kraken_output[,6])
    
    kraken_output <- kraken_output[kraken_output[,4] == taxonomy.list.names[i],]
    
    taxonomy_virus.metatranscriptome.df[[i]][,j] <- kraken_output[match(rownames(taxonomy_virus.metatranscriptome.df[[i]]),kraken_output[,6]),2]
  }
  taxonomy_virus.metatranscriptome.df[[i]][is.na(taxonomy_virus.metatranscriptome.df[[i]])] <- 0
}
```

#### Output data frames containing the raw reads for each taxonomic identifier

##### Metagenome

```{R}
for(i in 1:length(taxonomy_virus.metagenome.df)){
  write.table(taxonomy_virus.metagenome.df[[i]],
              paste0(WORKING.DIR,"/data_frames/taxonomy_virus_metagenome_",taxonomy.list.names[i],".tsv"),
              row.names = T,
              col.names = T,
              quote = F,
              sep = "\t")
}
```

##### Metatranscriptome

```{R}
for(i in 1:length(taxonomy_virus.metatranscriptome.df)){
  write.table(taxonomy_virus.metatranscriptome.df[[i]],
              paste0(WORKING.DIR,"/data_frames/taxonomy_virus_metatranscriptome_",taxonomy.list.names[i],".tsv"),
              row.names = T,
              col.names = T,
              quote = F,
              sep = "\t")
}
```

### Fungi

#### Set R inputs

```{R}
WORKING.DIR <- ""
```

#### Load packages and view sessionInfo

```{R}
sessionInfo()
```

```{R, eval = F}
R version 4.0.2 (2020-06-22)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 18363)

Matrix products: default

locale:
[1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252    LC_MONETARY=English_United States.1252
[4] LC_NUMERIC=C                           LC_TIME=English_United States.1252    

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

loaded via a namespace (and not attached):
[1] compiler_4.0.2 tools_4.0.2    knitr_1.29     xfun_0.16
```

#### Set list of taxonomy ranks of interest 

```{R}
taxonomy.list.names <- c("P","C","O","F","G","S")
```

#### Read in Kraken2 outputs for metagenome and metatranscriptome samples

```{R}
kraken_output.metagenome.list <- list.files(paste0(WORKING.DIR,"/kraken/metagenome/"),
											pattern="*fungi.kraken_report.txt")
kraken_output.metatranscriptome.list <- list.files(paste0(WORKING.DIR,"/kraken/metatranscriptome/"),
												   pattern="*fungi.kraken_report.txt")
```

#### Exclude metagenome samples with low read counts

The following DNA/metagenome samples do not have a sufficient number of reads for the analysis: 
Blank_234_1 (no sequencing data present)
Blank_234_2
UCS_0136_BKG

All RNA/metatranscriptome samples have a sufficient number of reads.

```{R}
kraken_output.metagenome.list <- grep("Blank_234_2",kraken_output.metagenome.list,value = T,invert = T)
kraken_output.metagenome.list <- grep("UCS_0136_BKG",kraken_output.metagenome.list,value = T,invert = T)
```

#### Compile list of detected taxonomic identifiers for metagenome and metatranscriptome samples

##### Metagenome

```{R}
taxonomy_fungi.metagenome.list <- vector(mode = "list", length = length(taxonomy.list.names))

for(i in 1:length(kraken_output.metagenome.list)){
  kraken_output <- read.delim(paste0(WORKING.DIR,"/kraken/metagenome/",kraken_output.metagenome.list[i]),
                              header = F)
  kraken_output[,6] <- gsub("^ {2,}","",kraken_output[,6])
  
  d.coords <- which(kraken_output[,4] == "D")
  
  if(length(which(kraken_output[,6] == "Eukaryota")) != 0){
    fungi.coords1 <- which(kraken_output[,6] == "Eukaryota")
    if(which(d.coords == fungi.coords1) == length(d.coords)){
      fungi.coords2 <- nrow(kraken_output)
    }else{
      fungi.coords2 <- d.coords[which(d.coords == fungi.coords1)+1]-1
    }
    kraken_output_fungi <- kraken_output[fungi.coords1:fungi.coords2,]
  }

  for(j in 1:length(taxonomy.list.names)){
    taxonomy_fungi.metagenome.list[[j]] <- unique(c(taxonomy_fungi.metagenome.list[[j]],kraken_output_fungi[kraken_output_fungi[,4] == taxonomy.list.names[j],6]))
  }
}
```

##### Metatranscriptome

```{R}
taxonomy_fungi.metatranscriptome.list <- vector(mode = "list", length = length(taxonomy.list.names))

for(i in 1:length(kraken_output.metatranscriptome.list)){
  kraken_output <- read.delim(paste0(WORKING.DIR,"/kraken/metatranscriptome/",kraken_output.metatranscriptome.list[i]),
                              header = F)
  kraken_output[,6] <- gsub("^ {2,}","",kraken_output[,6])
  
  d.coords <- which(kraken_output[,4] == "D")
  
  fungi.coords1 <- which(kraken_output[,6] == "Eukaryota")
  if(which(d.coords == fungi.coords1) == length(d.coords)){
    fungi.coords2 <- nrow(kraken_output)
  }else{
    fungi.coords2 <- d.coords[which(d.coords == fungi.coords1)+1]-1
  }
  
  kraken_output_fungi <- kraken_output[fungi.coords1:fungi.coords2,]

  for(j in 1:length(taxonomy.list.names)){
    taxonomy_fungi.metatranscriptome.list[[j]] <- unique(c(taxonomy_fungi.metatranscriptome.list[[j]],kraken_output_fungi[kraken_output_fungi[,4] == taxonomy.list.names[j],6]))
  }
}
```

#### Construct data frame of the raw reads for each taxonomic identifier

##### Metagenome

```{R}
taxonomy_fungi.metagenome.df <- vector(mode = "list", length = length(taxonomy.list.names))

for(i in 1:length(taxonomy.list.names)){
  taxonomy_fungi.metagenome.df[[i]] <- as.data.frame(matrix(0,
                                                        nrow=length(taxonomy_fungi.metagenome.list[[i]]),
                                                        ncol=length(kraken_output.metagenome.list)))
  rownames(taxonomy_fungi.metagenome.df[[i]]) <- taxonomy_fungi.metagenome.list[[i]]
  colnames(taxonomy_fungi.metagenome.df[[i]]) <- gsub("[.].*","",kraken_output.metagenome.list)
  
  for(j in 1:length(kraken_output.metagenome.list)){
    kraken_output <- read.delim(paste0(WORKING.DIR,"/kraken/metagenome/",kraken_output.metagenome.list[j]),
                                header = F)
    kraken_output[,6] <- gsub("^ {2,}","",kraken_output[,6])
   
    kraken_output <- kraken_output[kraken_output[,4] == taxonomy.list.names[i],]

    taxonomy_fungi.metagenome.df[[i]][,j] <- kraken_output[match(rownames(taxonomy_fungi.metagenome.df[[i]]),kraken_output[,6]),2]
  }
  taxonomy_fungi.metagenome.df[[i]][is.na(taxonomy_fungi.metagenome.df[[i]])] <- 0
}
```

##### Metatranscriptome

```{R}
taxonomy_fungi.metatranscriptome.df <- vector(mode = "list", length = length(taxonomy.list.names))

for(i in 1:length(taxonomy.list.names)){
  taxonomy_fungi.metatranscriptome.df[[i]] <- as.data.frame(matrix(0,
                                                        nrow=length(taxonomy_fungi.metatranscriptome.list[[i]]),
                                                        ncol=length(kraken_output.metatranscriptome.list)))
  rownames(taxonomy_fungi.metatranscriptome.df[[i]]) <- taxonomy_fungi.metatranscriptome.list[[i]]
  colnames(taxonomy_fungi.metatranscriptome.df[[i]]) <- gsub("[.].*","",kraken_output.metatranscriptome.list)
  
  for(j in 1:length(kraken_output.metatranscriptome.list)){
    kraken_output <- read.delim(paste0(WORKING.DIR,"/kraken/metatranscriptome/",kraken_output.metatranscriptome.list[j]),
                                header = F)
    kraken_output[,6] <- gsub("^ {2,}","",kraken_output[,6])
   
    kraken_output <- kraken_output[kraken_output[,4] == taxonomy.list.names[i],]

    taxonomy_fungi.metatranscriptome.df[[i]][,j] <- kraken_output[match(rownames(taxonomy_fungi.metatranscriptome.df[[i]]),kraken_output[,6]),2]
  }
  taxonomy_fungi.metatranscriptome.df[[i]][is.na(taxonomy_fungi.metatranscriptome.df[[i]])] <- 0
}
```

#### Output data frames containing the raw reads for each taxonomic identifier

##### Metagenome

```{R}
for(i in 1:length(taxonomy_fungi.metagenome.df)){
  write.table(taxonomy_fungi.metagenome.df[[i]],
              paste0(WORKING.DIR,"/data_frames/taxonomy_fungi_metagenome_",taxonomy.list.names[i],".tsv"),
              row.names = T,
              col.names = T,
              quote = F,
              sep = "\t")
}
```

##### Metatranscriptome

```{R}
for(i in 1:length(taxonomy_fungi.metatranscriptome.df)){
  write.table(taxonomy_fungi.metatranscriptome.df[[i]],
              paste0(WORKING.DIR,"/data_frames/taxonomy_fungi_metatranscriptome_",taxonomy.list.names[i],".tsv"),
              row.names = T,
              col.names = T,
              quote = F,
              sep = "\t")
}
```

## Create genera and species lists for DNA and RNA viruses with vertebrate, bacteria/archaea, and other hosts

### Set R inputs

```{R}
WORKING.DIR <- ""
```

### Load packages and view sessionInfo

```{R}
sessionInfo()
```

```{R, eval = F}
R version 4.0.2 (2020-06-22)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 18363)

Matrix products: default

locale:
[1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252    LC_MONETARY=English_United States.1252
[4] LC_NUMERIC=C                           LC_TIME=English_United States.1252    

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

loaded via a namespace (and not attached):
[1] compiler_4.0.2 tools_4.0.2    knitr_1.30     xfun_0.
```

### Read in list of mammalian virus taxonomy families

```{r}
virus_family.info <- read.delim(paste0(WORKING.DIR,"/maps/virus_family_info.tsv"),
                                header = F)
```

### Read in list of kraken reports

```{r}
kraken_report.file_list <- list.files(paste0(WORKING.DIR,"/kraken/"),
                                             recursive = T,
                                             pattern="kraken_report.txt")
kraken_report.file_list <- grep("antibiotic_resistance",kraken_report.file_list,value = T,invert = T)
kraken_report.file_list <- grep("metagenome/Blank_234_2",kraken_report.file_list,value = T,invert = T)
kraken_report.file_list <- grep("metagenome/UCS_0136_BKG",kraken_report.file_list,value = T,invert = T)
```

### Parse all kraken outputs to map mammalian viral genera and species taxa to their families

```{r}
dna_virus_genus.list <- list(vertebrates=list(),
                             bacteria_archaea=list(),
                             other=list())
rna_virus_genus.list <- list(vertebrates=list(),
                             bacteria_archaea=list(),
                             other=list())
dna_virus_species.list <- list(vertebrates=list(),
                             bacteria_archaea=list(),
                             other=list())
rna_virus_species.list <- list(vertebrates=list(),
                             bacteria_archaea=list(),
                             other=list())

for(i in 1:length(kraken_report.file_list)){
  kraken_report <- read.delim(paste0(WORKING.DIR,"/kraken/",kraken_report.file_list[i]))
  kraken_report[,6] <- gsub("^ {2,}","",kraken_report[,6])
  for(j in 1:nrow(virus_family.info)){
    if(nrow(kraken_report[kraken_report[,4] == "F" & kraken_report[,6] == virus_family.info[j,1],]) > 0){
      start_row <- which(kraken_report[,4] == "F" & kraken_report[,6] == virus_family.info[j,1])
      
      stop_row <- intersect(grep("^G",kraken_report[,4],invert = T),grep("^S",kraken_report[,4],invert = T))
      stop_row <- intersect(stop_row,grep("^F[0-9]",kraken_report[,4],invert = T))
      stop_row <- stop_row[stop_row > start_row][1]-1
      stop_row <- ifelse(is.na(stop_row),nrow(kraken_report),stop_row)
      kraken_report.subset <- kraken_report[start_row:stop_row,]
      
      if(virus_family.info[j,2] == "DNA"){
        if(virus_family.info[j,3] == "vertebrates"){
          dna_virus_genus.list$vertebrates = c(unlist(dna_virus_genus.list$vertebrates), kraken_report.subset[kraken_report.subset[,4] == "G",6])
          dna_virus_species.list$vertebrates = c(unlist(dna_virus_species.list$vertebrates), kraken_report.subset[kraken_report.subset[,4] == "S",6])
        }else if(virus_family.info[j,3] == "other"){
          dna_virus_genus.list$other = c(unlist(dna_virus_genus.list$other), kraken_report.subset[kraken_report.subset[,4] == "G",6])
          dna_virus_species.list$other = c(unlist(dna_virus_species.list$other), kraken_report.subset[kraken_report.subset[,4] == "S",6])
        }else{
          dna_virus_genus.list$bacteria_archaea = c(unlist(dna_virus_genus.list$bacteria_archaea), kraken_report.subset[kraken_report.subset[,4] == "G",6])
          dna_virus_species.list$bacteria_archaea = c(unlist(dna_virus_species.list$bacteria_archaea), kraken_report.subset[kraken_report.subset[,4] == "S",6])
        }
      }else{
        if(virus_family.info[j,3] == "vertebrates"){
          rna_virus_genus.list$vertebrates = c(unlist(rna_virus_genus.list$vertebrates), kraken_report.subset[kraken_report.subset[,4] == "G",6])
          rna_virus_species.list$vertebrates = c(unlist(rna_virus_species.list$vertebrates), kraken_report.subset[kraken_report.subset[,4] == "S",6])
        }else if(virus_family.info[j,3] == "other"){
          rna_virus_genus.list$other = c(unlist(rna_virus_genus.list$other), kraken_report.subset[kraken_report.subset[,4] == "G",6])
          rna_virus_species.list$other = c(unlist(rna_virus_species.list$other), kraken_report.subset[kraken_report.subset[,4] == "S",6])
        }else{
          rna_virus_genus.list$bacteria_archaea = c(unlist(rna_virus_genus.list$bacteria_archaea), kraken_report.subset[kraken_report.subset[,4] == "G",6])
          rna_virus_species.list$bacteria_archaea = c(unlist(rna_virus_species.list$bacteria_archaea), kraken_report.subset[kraken_report.subset[,4] == "S",6])
        }
      }
    }
  }
}
```

### Manually remove duplicate genera entries incorrectly binned

```{r}
dna_virus_genus.list$other <- dna_virus_genus.list$other[dna_virus_genus.list$other != "Salterprovirus"]
dna_virus_genus.list$other <- dna_virus_genus.list$other[dna_virus_genus.list$other != "Dinodnavirus"]

rna_virus_genus.list$other <- rna_virus_genus.list$other[rna_virus_genus.list$other != "Idaeovirus"]
```

### Output a list containing all mammalian virus genera and species split into DNA, RNA viruses and host categories

##### DNA

```{r}
write.table(as.data.frame(unique(dna_virus_genus.list$vertebrates)),
            paste0(WORKING.DIR,"/maps/dna_vertebrate_viral_genera.list"),
            col.names = F,
            row.names = F,
            quote = F,
            sep = "\t")
write.table(as.data.frame(unique(dna_virus_genus.list$bacteria_archaea)),
            paste0(WORKING.DIR,"/maps/dna_bacteria_archaea_viral_genera.list"),
            col.names = F,
            row.names = F,
            quote = F,
            sep = "\t")
write.table(as.data.frame(unique(dna_virus_genus.list$other)),
            paste0(WORKING.DIR,"/maps/dna_other_viral_genera.list"),
            col.names = F,
            row.names = F,
            quote = F,
            sep = "\t")

write.table(as.data.frame(unique(dna_virus_species.list$vertebrates)),
            paste0(WORKING.DIR,"/maps/dna_vertebrate_viral_species.list"),
            col.names = F,
            row.names = F,
            quote = F,
            sep = "\t")
write.table(as.data.frame(unique(dna_virus_species.list$bacteria_archaea)),
            paste0(WORKING.DIR,"/maps/dna_bacteria_archaea_viral_species.list"),
            col.names = F,
            row.names = F,
            quote = F,
            sep = "\t")
write.table(as.data.frame(unique(dna_virus_species.list$other)),
            paste0(WORKING.DIR,"/maps/dna_other_viral_species.list"),
            col.names = F,
            row.names = F,
            quote = F,
            sep = "\t")
```

##### RNA

```{r}
write.table(as.data.frame(unique(rna_virus_genus.list$vertebrates)),
            paste0(WORKING.DIR,"/maps/rna_vertebrate_viral_genera.list"),
            col.names = F,
            row.names = F,
            quote = F,
            sep = "\t")
write.table(as.data.frame(unique(rna_virus_genus.list$bacteria_archaea)),
            paste0(WORKING.DIR,"/maps/rna_bacteria_archaea_viral_genera.list"),
            col.names = F,
            row.names = F,
            quote = F,
            sep = "\t")
write.table(as.data.frame(unique(rna_virus_genus.list$other)),
            paste0(WORKING.DIR,"/maps/rna_other_viral_genera.list"),
            col.names = F,
            row.names = F,
            quote = F,
            sep = "\t")

write.table(as.data.frame(unique(rna_virus_species.list$vertebrates)),
            paste0(WORKING.DIR,"/maps/rna_vertebrate_viral_species.list"),
            col.names = F,
            row.names = F,
            quote = F,
            sep = "\t")
write.table(as.data.frame(unique(rna_virus_species.list$bacteria_archaea)),
            paste0(WORKING.DIR,"/maps/rna_bacteria_archaea_viral_species.list"),
            col.names = F,
            row.names = F,
            quote = F,
            sep = "\t")
write.table(as.data.frame(unique(rna_virus_species.list$other)),
            paste0(WORKING.DIR,"/maps/rna_other_viral_species.list"),
            col.names = F,
            row.names = F,
            quote = F,
            sep = "\t")
```

## Combine taxonomy raw counts data frames into one data frame

### Genera

##### Metagenome

```{bash, eval = F}
head -n1 "$WORKING_DIR"/data_frames/taxonomy_bacteria_metagenome_G.tsv > "$WORKING_DIR"/data_frames/metagenome_header.tsv

cat "$WORKING_DIR"/data_frames/taxonomy_bacteria_metagenome_G.tsv "$WORKING_DIR"/data_frames/taxonomy_DNAbacteriaarchaeavirus_metagenome_G.sarscov2.tsv "$WORKING_DIR"/data_frames/taxonomy_DNAothervirus_metagenome_G.sarscov2.tsv "$WORKING_DIR"/data_frames/taxonomy_DNAvertebratevirus_metagenome_G.sarscov2.tsv "$WORKING_DIR"/data_frames/taxonomy_RNAbacteriaarchaeavirus_metagenome_G.sarscov2.tsv "$WORKING_DIR"/data_frames/taxonomy_RNAothervirus_metagenome_G.sarscov2.tsv "$WORKING_DIR"/data_frames/taxonomy_RNAvertebratevirus_metagenome_G.sarscov2.tsv "$WORKING_DIR"/data_frames/taxonomy_fungi_metagenome_G.tsv | grep -v "^UCS" | sort -n | uniq > "$WORKING_DIR"/data_frames/metagenome_body.tsv

cat "$WORKING_DIR"/data_frames/metagenome_header.tsv "$WORKING_DIR"/data_frames/metagenome_body.tsv > counts_taxonomy_all_metagenome_G.tsv
rm "$WORKING_DIR"/data_frames/metagenome_header.tsv "$WORKING_DIR"/data_frames/metagenome_body.tsv
```

##### Metatranscriptome

```{bash, eval = F}
head -n1 "$WORKING_DIR"/data_frames/taxonomy_bacteria_metatranscriptome_G.tsv > "$WORKING_DIR"/data_frames/metatranscriptome_header.tsv

cat "$WORKING_DIR"/data_frames/taxonomy_bacteria_metatranscriptome_G.tsv "$WORKING_DIR"/data_frames/taxonomy_DNAbacteriaarchaeavirus_metatranscriptome_G.sarscov2.tsv "$WORKING_DIR"/data_frames/taxonomy_DNAothervirus_metatranscriptome_G.sarscov2.tsv "$WORKING_DIR"/data_frames/taxonomy_DNAvertebratevirus_metatranscriptome_G.sarscov2.tsv "$WORKING_DIR"/data_frames/taxonomy_RNAbacteriaarchaeavirus_metatranscriptome_G.sarscov2.tsv "$WORKING_DIR"/data_frames/taxonomy_RNAothervirus_metatranscriptome_G.sarscov2.tsv "$WORKING_DIR"/data_frames/taxonomy_RNAvertebratevirus_metatranscriptome_G.sarscov2.tsv "$WORKING_DIR"/data_frames/taxonomy_fungi_metatranscriptome_G.tsv | grep -v "^UCS" | sort -n | uniq > "$WORKING_DIR"/data_frames/metatranscriptome_body.tsv

cat "$WORKING_DIR"/data_frames/metatranscriptome_header.tsv "$WORKING_DIR"/data_frames/metatranscriptome_body.tsv > "$WORKING_DIR"/data_frames/counts_taxonomy_all_metatranscriptome_G.tsv
rm "$WORKING_DIR"/data_frames/metatranscriptome_header.tsv "$WORKING_DIR"/data_frames/metatranscriptome_body.tsv
```

### Species

##### Metagenome

```{bash, eval = F}
head -n1 "$WORKING_DIR"/data_frames/taxonomy_bacteria_metagenome_S.tsv > "$WORKING_DIR"/data_frames/metagenome_header.tsv

cat "$WORKING_DIR"/data_frames/taxonomy_bacteria_metagenome_S.tsv "$WORKING_DIR"/data_frames/taxonomy_DNAbacteriaarchaeavirus_metagenome_S.sarscov2.tsv "$WORKING_DIR"/data_frames/taxonomy_DNAothervirus_metagenome_S.sarscov2.tsv "$WORKING_DIR"/data_frames/taxonomy_DNAvertebratevirus_metagenome_S.sarscov2.tsv "$WORKING_DIR"/data_frames/taxonomy_RNAbacteriaarchaeavirus_metagenome_S.sarscov2.tsv "$WORKING_DIR"/data_frames/taxonomy_RNAothervirus_metagenome_S.sarscov2.tsv "$WORKING_DIR"/data_frames/taxonomy_RNAvertebratevirus_metagenome_S.sarscov2.tsv "$WORKING_DIR"/data_frames/taxonomy_fungi_metagenome_S.tsv | grep -v "^UCS" | sort -n | uniq > "$WORKING_DIR"/data_frames/metagenome_body.tsv

cat "$WORKING_DIR"/data_frames/metagenome_header.tsv "$WORKING_DIR"/data_frames/metagenome_body.tsv > "$WORKING_DIR"/data_frames/counts_taxonomy_all_metagenome_S.tsv
rm "$WORKING_DIR"/data_frames/metagenome_header.tsv "$WORKING_DIR"/data_frames/metagenome_body.tsv
```

##### Metatranscriptome

```{bash, eval = F}
head -n1 "$WORKING_DIR"/data_frames/taxonomy_bacteria_metatranscriptome_S.tsv > "$WORKING_DIR"/data_frames/metatranscriptome_header.tsv

cat "$WORKING_DIR"/data_frames/taxonomy_bacteria_metatranscriptome_S.tsv "$WORKING_DIR"/data_frames/taxonomy_DNAbacteriaarchaeavirus_metatranscriptome_S.sarscov2.tsv "$WORKING_DIR"/data_frames/taxonomy_DNAothervirus_metatranscriptome_S.sarscov2.tsv "$WORKING_DIR"/data_frames/taxonomy_DNAvertebratevirus_metatranscriptome_S.sarscov2.tsv "$WORKING_DIR"/data_frames/taxonomy_RNAbacteriaarchaeavirus_metatranscriptome_S.sarscov2.tsv "$WORKING_DIR"/data_frames/taxonomy_RNAothervirus_metatranscriptome_S.sarscov2.tsv "$WORKING_DIR"/data_frames/taxonomy_RNAvertebratevirus_metatranscriptome_S.sarscov2.tsv "$WORKING_DIR"/data_frames/taxonomy_fungi_metatranscriptome_S.tsv | grep -v "^UCS" | sort -n | uniq > "$WORKING_DIR"/data_frames/metatranscriptome_body.tsv

cat "$WORKING_DIR"/data_frames/metatranscriptome_header.tsv "$WORKING_DIR"/data_frames/metatranscriptome_body.tsv > counts_taxonomy_all_metatranscriptome_S.tsv
rm "$WORKING_DIR"/data_frames/metatranscriptome_header.tsv "$WORKING_DIR"/data_frames/metatranscriptome_body.tsv
```

## Add species categorizations to counts data frames and create combined relative abundance dataframes

### Set R inputs

```{r}
WORKING.DIR <- ""
```

### Load packages and view sessionInfo

```{R}
sessionInfo()
```

```{R, eval = F}
R version 4.0.2 (2020-06-22)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 18363)

Matrix products: default

locale:
[1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252    LC_MONETARY=English_United States.1252
[4] LC_NUMERIC=C                           LC_TIME=English_United States.1252    

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

loaded via a namespace (and not attached):
[1] compiler_4.0.2 tools_4.0.2    knitr_1.30     xfun_0.
```

### Compile list of taxonomic counts files

```{r}
counts.filelist <- grep("species",list.files(paste0(WORKING.DIR,"/data_frames/"),pattern="^counts_taxonomy_all_",full.names = T),value = T,invert=T)
#counts.filelist <- grep("_all_",list.files(paste0(WORKING.DIR,"/data_frames/"),pattern="^RAREFY",full.names = T),value = T)
```

### Compile taxa lists

##### Genera

```{r}
genus <- list()

genus$bacteria <- unique(c(rownames(read.delim(paste0(WORKING.DIR,"/data_frames/taxonomy_bacteria_metagenome_G.tsv"))),
                           rownames(read.delim(paste0(WORKING.DIR,"/data_frames/taxonomy_bacteria_metatranscriptome_G.tsv"))),
                           rownames(read.delim(paste0(WORKING.DIR,"/data_frames/taxonomy_bacteria_metagenome_extra_G.tsv")))))

genus$fungi <- unique(c(rownames(read.delim(paste0(WORKING.DIR,"/data_frames/taxonomy_fungi_metagenome_G.tsv"))),
                        rownames(read.delim(paste0(WORKING.DIR,"/data_frames/taxonomy_fungi_metatranscriptome_G.tsv"))),
                        rownames(read.delim(paste0(WORKING.DIR,"/data_frames/taxonomy_fungi_metagenome_extra_G.tsv")))))

genus$DNAvertebratevirus <- read.delim(paste0(WORKING.DIR,"/maps/dna_vertebrate_viral_genera.list"), header = F)[,1]
genus$DNAbacteriaarchaeavirus <- read.delim(paste0(WORKING.DIR,"/maps/dna_bacteria_archaea_viral_genera.list"), header = F)[,1]
genus$DNAothervirus <- read.delim(paste0(WORKING.DIR,"/maps/dna_other_viral_genera.list"), header = F)[,1]

genus$RNAvertebratevirus <- read.delim(paste0(WORKING.DIR,"/maps/rna_vertebrate_viral_genera.list"), header = F)[,1]
genus$RNAbacteriaarchaeavirus <- read.delim(paste0(WORKING.DIR,"/maps/rna_bacteria_archaea_viral_genera.list"), header = F)[,1]
genus$RNAothervirus <- read.delim(paste0(WORKING.DIR,"/maps/rna_other_viral_genera.list"), header = F)[,1]

genus_map.df <- as.data.frame(matrix(nrow=0,
                                     ncol=2))

for(i in 1:length(genus)){
  genus_map.df <- as.data.frame(rbind(genus_map.df,
                                      cbind(genus[[i]],names(genus)[i])))
}
```

##### Species

```{r}
species <- list()

species$bacteria <- unique(c(rownames(read.delim(paste0(WORKING.DIR,"/data_frames/taxonomy_bacteria_metagenome_S.tsv"))),
                             rownames(read.delim(paste0(WORKING.DIR,"/data_frames/taxonomy_bacteria_metatranscriptome_S.tsv"))),
                             rownames(read.delim(paste0(WORKING.DIR,"/data_frames/taxonomy_bacteria_metagenome_extra_S.tsv")))))

species$fungi <- unique(c(rownames(read.delim(paste0(WORKING.DIR,"/data_frames/taxonomy_fungi_metagenome_S.tsv"))),
                          rownames(read.delim(paste0(WORKING.DIR,"/data_frames/taxonomy_fungi_metatranscriptome_S.tsv"))),
                          rownames(read.delim(paste0(WORKING.DIR,"/data_frames/taxonomy_fungi_metagenome_extra_S.tsv")))))

species$DNAvertebratevirus <- read.delim(paste0(WORKING.DIR,"/maps/dna_vertebrate_viral_species.list"), header = F)[,1]
species$DNAbacteriaarchaeavirus <- read.delim(paste0(WORKING.DIR,"/maps/dna_bacteria_archaea_viral_species.list"), header = F)[,1]
species$DNAothervirus <- read.delim(paste0(WORKING.DIR,"/maps/dna_other_viral_species.list"), header = F)[,1]

species$RNAvertebratevirus <- read.delim(paste0(WORKING.DIR,"/maps/rna_vertebrate_viral_species.list"), header = F)[,1]
species$RNAbacteriaarchaeavirus <- read.delim(paste0(WORKING.DIR,"/maps/rna_bacteria_archaea_viral_species.list"), header = F)[,1]
species$RNAothervirus <- read.delim(paste0(WORKING.DIR,"/maps/rna_other_viral_species.list"), header = F)[,1]

species_map.df <- as.data.frame(matrix(nrow=0,
                                       ncol=2))

for(i in 1:length(species)){
  species_map.df <- as.data.frame(rbind(species_map.df,
                                        cbind(species[[i]],names(species)[i])))
}
```

### Manually curate categorizations of viruses not categorized by Kraken2

```{r}
manualcurate.rnaothervirus.list <- c("Hubei picorna-like virus 27","Hubei picorna-like virus 15","Wuhan insect virus 19",
                                     "Wenling zhaovirus-like virus 1","Hubei sobemo-like virus 11", "Wenling crustacean virus 3",
                                     "Wenling crustacean virus 14","Sanxia water strider virus 7","Beihai picorna-like virus 59",
                                     "Shuangao lacewing virus 2","Hubei endorna-like virus 1", "Wenzhou picorna-like virus 45",
                                     "Wuhan spider virus 4","Shuangao insect virus 12","Beihai mantis shrimp virus 5",
                                     "Hubei insect virus 2","Shayang ascaridia galli virus 2","Wuhan coneheads virus 1",
                                     "Nephila clavipes virus 2","Beihai mollusks virus 1","Hubei diptera virus 11",
                                     "Wuhan insect virus 13","Wuhan heteroptera virus 1","Beihai sipunculid worm virus 5",
                                     "Changjiang tombus-like virus 22","Wuhan nido-like virus 1","Xingshan nematode virus 2",
                                     "Changjiang picorna-like virus 10","Wenzhou hepe-like virus 1","Shayang spider virus 4",
                                     "Beihai tombus-like virus 18","Wenzhou shrimp virus 7","Beihai sesarmid crab virus 1",
                                     "Halastavi arva RNA virus","Hubei sobemo-like virus 42","Wenzhou gastropodes virus 2",
                                     "Hubei myriapoda virus 2","Wenzhou picorna-like virus 36","Drosophila subobscura Nora virus",
                                     "Wuhan pillworm virus 1","Wuhan house centipede virus 1","Wenzhou shrimp virus 5",
                                     "Sanxia tombus-like virus 6","Wenling crustacean virus 4","Shahe picorna-like virus 13",
                                     "Sanxia atyid shrimp virus 2","Beihai picorna-like virus 4","Xinzhou nematode virus 1",
                                     "Osedax japonicus RNA virus 1","Beihai hepe-like virus 9","Changjiang picorna-like virus 11",
                                     "Wuhan spider virus 6", "Hubei picorna-like virus 52","Hubei tetragnatha maxillosa virus 3",
                                     "Beihai barnacle viurs 1","Sanxia picorna-like virus 5","Beihai picorna-like virus 44",
                                     "Shahe isopoda virus 1","Changjiang picorna-like virus 6" ,"Hubei picorna-like virus 40",
                                     "Beihai tombus-like virus 9","Lonestar tick chuvirus 1","Beihai picorna-like virus 65",
                                     "Wenzhou picorna-like virus 52","Hubei tombus-like virus 13","Hubei picorna-like virus 26",
                                     "Hubei dimarhabdovirus virus 1","Adelphocoris suturalis virus","Ying Kou virus",
                                     "Hubei odonate virus 2","Wenzhou tombus-like virus 11","Hubei tombus-like virus 18",
                                     "Wallerfield virus","Hubei picorna-like virus 56","Hubei earwig virus 3",
                                     "Beihai picorna-like virus 87","Cordoba virus","Lodeiro virus",
                                     "Culex negev-like virus 3","Wabat virus")
                                     
manualcurate.dnaothervirus.list <-  c( "Bufonid herpesvirus 1","Cellulophaga phage phi48:2","Longjawed orbweaver circular virus 1",
                                       "Halovirus HCTV-2","Sewage-associated circular DNA virus-2","Halovirus HHTV-1","Halovirus HGTV-1",
                                       "Halovirus VNH-1","Halovirus HSTV-1","Thermococcus prieurii virus 1")

manualcurate.rnaothervirus.species_map.df <- as.data.frame(cbind(manualcurate.rnaothervirus.list,"RNAothervirus"))
manualcurate.dnaothervirus.species_map.df <- as.data.frame(cbind(manualcurate.dnaothervirus.list,"DNAothervirus"))

colnames(manualcurate.rnaothervirus.species_map.df) <- c("V1","V2")
colnames(manualcurate.dnaothervirus.species_map.df) <- c("V1","V2")

species_map.df <- as.data.frame(rbind(species_map.df,
                                      manualcurate.rnaothervirus.species_map.df,
                                      manualcurate.dnaothervirus.species_map.df))
```

### Set function for calculating relative abundance

```{r}
calculate_relab <- function(df){
  for(i in 1:ncol(df)){
    df[,i] <- df[,i]/sum(df[,i])
    
  }
  return(df)
}
```

### Set list of samples to exclude due to IRB issues

```{r}
sample_exclude.list <- c("UCS_0045","UCS_0052","UCS_0170","UCS_0171")
```

### Write counts and relative abundance tables split by taxa categories

```{r}
for(i in 1:length(counts.filelist)){
  counts.df <- read.delim(counts.filelist[i])
  counts.df <- counts.df[,!(gsub("_[BU].*","",colnames(counts.df)) %in% sample_exclude.list)]
  
  if(grepl("metagenome",counts.filelist[i],fixed=T)){
    category <- species_map.df[match(rownames(counts.df),species_map.df[,1]),2] 
    counts.df <- counts.df[!(grepl("RNA",category)),]
  }
  
  counts.relab.df <- calculate_relab(counts.df)
  if(grepl("_G.tsv",counts.filelist[i],fixed = T)){
    for(j in 1:length(genus)){
      counts.subset.df <- counts.df[rownames(counts.df) %in% genus[[j]],]
      write.table(counts.subset.df,
                  gsub("_all_",paste0("_",names(genus)[j],"_"),counts.filelist[i]),
                  row.names = T,
                  col.names = T,
                  quote = F,
                  sep = "\t")
      
      counts.relab.subset.df <- counts.relab.df[rownames(counts.relab.df) %in% genus[[j]],]
      write.table(counts.relab.subset.df,
                  gsub("counts_","relab_",gsub("_all_",paste0("_",names(genus)[j],"_"),counts.filelist[i])),
                  row.names = T,
                  col.names = T,
                  quote = F,
                  sep = "\t")
    }
  }
  else if(grepl("_S.tsv",counts.filelist[i],fixed = T)){
    for(j in 1:length(species)){
      counts.subset.df <- counts.df[rownames(counts.df) %in% species[[j]],]
      write.table(counts.subset.df,
                  gsub("_all_",paste0("_",names(species)[j],"_"),counts.filelist[i]),
                  row.names = T,
                  col.names = T,
                  quote = F,
                  sep = "\t")
      
      counts.relab.subset.df <- counts.relab.df[rownames(counts.relab.df) %in% species[[j]],]
      write.table(counts.relab.subset.df,
                  gsub("counts_","relab_",gsub("_all_",paste0("_",names(species)[j],"_"),counts.filelist[i])),
                  row.names = T,
                  col.names = T,
                  quote = F,
                  sep = "\t")
    }
  }
}
```

### Write combined counts and relative abundance table

```{r}
for(i in 1:length(counts.filelist)){
  counts.df <- read.delim(counts.filelist[i])
  counts.df <- counts.df[,!(gsub("_[BU].*","",colnames(counts.df)) %in% sample_exclude.list)]
  
  if(grepl("metagenome",counts.filelist[i],fixed=T)){
    category <- species_map.df[match(rownames(counts.df),species_map.df[,1]),2] 
    counts.df <- counts.df[!(grepl("RNA",category)),]
  }
  
  counts.relab.df <- calculate_relab(counts.df)
  
  if(grepl("_G.tsv",counts.filelist[i],fixed = T)){
    counts.output.df <- counts.df
    counts.relab.output.df <- counts.relab.df
    
    counts.output.df$category <- genus_map.df[match(rownames(counts.df),genus_map.df[,1]),2]
    counts.relab.output.df$category <- genus_map.df[match(rownames(counts.relab.df),genus_map.df[,1]),2]

    write.table(counts.output.df,
                gsub(".tsv$",".speciescategorized.tsv",counts.filelist[i]),
                row.names = T,
                col.names = T,
                quote = F,
                sep = "\t")
  
    write.table(counts.relab.output.df,
                gsub("counts_","relab_",gsub(".tsv$",".speciescategorized.tsv",counts.filelist[i])),
                row.names = T,
                col.names = T,
                quote = F,
                sep = "\t")
  }else if(grepl("_S.tsv",counts.filelist[i],fixed = T)){
    counts.output.df <- counts.df
    counts.relab.output.df <- counts.relab.df
    
    counts.output.df$category <- species_map.df[match(rownames(counts.df),species_map.df[,1]),2]
    counts.relab.output.df$category <- species_map.df[match(rownames(counts.relab.df),species_map.df[,1]),2]

    write.table(counts.output.df,
                gsub(".tsv$",".speciescategorized.tsv",counts.filelist[i]),
                row.names = T,
                col.names = T,
                quote = F,
                sep = "\t")
  
    write.table(counts.relab.output.df,
                gsub("counts_","relab_",gsub(".tsv$",".speciescategorized.tsv",counts.filelist[i])),
                row.names = T,
                col.names = T,
                quote = F,
                sep = "\t")
    
    write.table(counts.output.df[,c(grep("BAL",colnames(counts.output.df)),ncol(counts.output.df))],
                gsub(".tsv$",".speciescategorized.bal.tsv",counts.filelist[i]),
                row.names = T,
                col.names = T,
                quote = F,
                sep = "\t")
  
    write.table(counts.relab.output.df[,c(grep("BAL",colnames(counts.relab.output.df)),ncol(counts.relab.output.df))],
                gsub("counts_","relab_",gsub(".tsv$",".speciescategorized.bal.tsv",counts.filelist[i])),
                row.names = T,
                col.names = T,
                quote = F,
                sep = "\t")
  }
}
```

## Plot box plot of most abundant OTUs

### Set R inputs

```{r}
WORKING.DIR <- ""
```

### Load packages and view sessionInfo

```{r}
library(cowplot)
library(ggplot2)
library(matrixStats)
library(reshape2)

sessionInfo()
```

```{r,eval = F}
R version 4.0.2 (2020-06-22)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 18363)

Matrix products: default

locale:
[1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252    LC_MONETARY=English_United States.1252
[4] LC_NUMERIC=C                           LC_TIME=English_United States.1252    

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] reshape2_1.4.4     matrixStats_0.56.0 ggplot2_3.3.2      cowplot_1.1.0     

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.5       rstudioapi_0.11  knitr_1.30       magrittr_1.5     tidyselect_1.1.0 munsell_0.5.0    colorspace_1.4-1
 [8] R6_2.4.1         rlang_0.4.7      stringr_1.4.0    plyr_1.8.6       dplyr_1.0.2      tools_4.0.2      grid_4.0.2      
[15] gtable_0.3.0     xfun_0.17        withr_2.3.0      ellipsis_0.3.1   tibble_3.0.3     lifecycle_0.2.0  crayon_1.3.4    
[22] purrr_0.3.4      vctrs_0.3.4      glue_1.4.2       stringi_1.5.3    compiler_4.0.2   pillar_1.4.6     generics_0.0.2  
[29] scales_1.1.1     pkgconfig_2.0.3 
```

### Read in taxonomic relative abundance files

```{r}
taxonomy_relab.genus.metagenome.df <-  read.delim(paste0(WORKING.DIR,"/data_frames/relab_taxonomy_all_metagenome_G.speciescategorized.tsv"))
taxonomy_relab.genus.metatranscriptome.df <- read.delim(paste0(WORKING.DIR,"/data_frames/relab_taxonomy_all_metatranscriptome_G.speciescategorized.tsv"))

taxonomy_relab.species.metagenome.df <-  read.delim(paste0(WORKING.DIR,"/data_frames/relab_taxonomy_all_metagenome_S.speciescategorized.tsv"))
taxonomy_relab.species.metatranscriptome.df <- read.delim(paste0(WORKING.DIR,"/data_frames/relab_taxonomy_all_metatranscriptome_S.speciescategorized.tsv"))
```

### Read in meta-data maps

```{r}
metagenome.map <- read.delim(paste0(WORKING.DIR,"/maps/201015.WGS.txt"),
                             row.names = 1)
metatranscriptome.map <- read.delim(paste0(WORKING.DIR,"/maps/201015.RNA.txt"),
                                    row.names = 1)

rownames(metagenome.map) <- gsub("[.]","_",rownames(metagenome.map))
rownames(metatranscriptome.map) <- gsub("[.]","_",rownames(metatranscriptome.map))

metagenome.map <- metagenome.map[match(colnames(taxonomy_relab.species.metagenome.df),gsub("[.]","_",rownames(metagenome.map))),]
metatranscriptome.map <- metatranscriptome.map[match(colnames(taxonomy_relab.species.metatranscriptome.df),gsub("[.]","_",rownames(metatranscriptome.map))),]
```

### Read in curated list of vertebrate viruses

```{r}
dna_vertebrate_virus_species <- read.delim(paste0(WORKING.DIR,"/maps/dna_vertebrate_viral_species_Human_curated.list.txt"), header = F)[,1]
rna_vertebrate_virus_species <- read.delim(paste0(WORKING.DIR,"/maps/rna_vertebrate_viral_species_Human_curated.list.txt"), header = F)[,1]
```

### Set function for creating box plots

```{r}
plot_otu_box <- function(relab.df){
  n=20
  
  relab.df <- log2(relab.df*100)
  
  relab.df <- relab.df[order(rowMedians(as.matrix(relab.df)),decreasing=T),]
  relab.df <- relab.df[1:n,]
  
  plot.df <- melt(as.matrix(relab.df))
  plot.df[,1] <- factor(plot.df[,1],levels=rev(rownames(relab.df)[1:n]))

  boxplot <- ggplot(mapping=aes(x=!!plot.df[,1],
                                y=!!plot.df[,3]))+
  geom_boxplot()+
  labs(x="",y="relative abundance")+
  coord_flip()+
  theme_bw()

  return(boxplot)
}

plot_otu_box_setorder <- function(relab.df,order.df){
  n=50
  
  #relab.df <- log2(relab.df*100)
 # order.df <- log2(order.df*100)
  
  order.df <- order.df[order(rowMedians(as.matrix(order.df)),decreasing=T),]
  order.df <- order.df[1:n,]
  
  plot.df <- melt(as.matrix(relab.df))
  plot.df <- plot.df[plot.df[,1] %in% rownames(order.df),]
  plot.df[,1] <- factor(plot.df[,1],levels=rev(rownames(order.df)[1:n]))
  
  boxplot <- ggplot(mapping=aes(x=!!plot.df[,1],
                                y=!!plot.df[,3]))+
  geom_boxplot()+
  labs(x="",y="relative abundance")+
  coord_flip()+
  #coord_flip(ylim=c(0,1))+
  theme_bw()

  return(boxplot)
}

plot_otu_box_pullxlim <- function(order.df,categories,category.levels){
  
  ceiling_dec <- function(x, level=1) round(x + 5*10^(-level-1), level)
  
  categories.df <- as.data.frame(cbind(rownames(order.df),
                                       categories))
   
  order.df <- order.df[order(rowMedians(as.matrix(order.df)),decreasing=T),]
  
  xlim.vector <- c()
  for(i in 1:length(category.levels)){
    categories.df.subset <- categories.df[match(rownames(order.df),categories.df[,1]),]
    order.df.subset <- order.df[match(categories.df.subset[categories.df.subset[,2] == category.levels[i],1],rownames(order.df)),]
    
    xlim.vector[i] <- ceiling_dec(max(order.df.subset),1)
  }
  return(xlim.vector)
}

plot_otu_box_facetplot <- function(relab.df,order.df,categories,category.levels,color.levels,xaxis.lim){
  n=50
  
  #relab.df <- log2(relab.df*100)
  #order.df <- log2(order.df*100)
  
  categories <- categories[order(rowMedians(as.matrix(relab.df)),decreasing=T)]
  
  order.df <- order.df[order(rowMedians(as.matrix(order.df)),decreasing=T),]
  relab.df <- relab.df[order(rowMedians(as.matrix(relab.df)),decreasing=T),]
  #order.df <- order.df[1:n,]
  
  categories.df <- as.data.frame(cbind(rownames(relab.df),
                                       categories,
                                       seq(1,nrow(order.df))[match(rownames(relab.df),rownames(order.df))],
                                       seq(1,nrow(relab.df))))
  categories.df <- categories.df[order(as.numeric(as.character(categories.df[,3]))),]
  
  keeptaxa.vector <- c()
  for(i in 1:length(category.levels)){
    keeptaxa.vector <- c(keeptaxa.vector,categories.df[which(categories.df[,2] == category.levels[i])[1],1])
  }
  keeptaxa.vector <- unique(c(keeptaxa.vector,categories.df[1:n,1]))
  
  categories.df <- categories.df[categories.df[,1] %in% keeptaxa.vector,]
  
  plot.df <- melt(as.matrix(relab.df))
  plot.df <- plot.df[plot.df[,1] %in% keeptaxa.vector,]

  plot.df[,4] <- categories.df[match(plot.df[,1],categories.df[,1]),2]
  plot.df[,4] <- factor(plot.df[,4],levels=category.levels)
  
  levels <- rev(paste0(categories.df[,1]," (",categories.df[,4],")"))
  plot.df[,1] <- paste0(plot.df[,1]," (", categories.df[match(plot.df[,1],categories.df[,1]),4],")")
  plot.df[,1] <- factor(plot.df[,1],levels=rev(paste0(categories.df[,1]," (",categories.df[,4],")")))
  
  plot.list <- list()
  for(i in 1:length(category.levels)){
    plot.df.subset <- plot.df[plot.df[,4] == category.levels[i],]
    plot.list[[i]] <- ggplot(mapping=aes(x=!!plot.df.subset[,3]*100+1,
                                         y=!!plot.df.subset[,1]))+
      geom_boxplot(fill=color.levels[i])+
      #coord_cartesian(xlim = c(1,xaxis.lim[i]*100+1))+
      coord_cartesian(xlim = c(1,100))+
      scale_x_continuous(trans = "log10")+
      labs(x="",y="")+
      #coord_flip(ylim=c(0,xaxis.lim[i]))+
      #coord_flip(ylim=c(0,1))+
      guides(fill=F)+
      theme_bw()+
      theme(axis.title.y = element_blank())
  }
  return(plot.list)
}
```

### Ranks all taxa in each of the different datasets

```{r}
find_rank_by_median <- function(df,taxon.vector){
  df <- df[,colnames(df) != "category"]
  rank <- rank(-rowMedians(as.matrix(df)))
  return(rank[match(taxon.vector,rownames(df))])
}

find_median_relab <- function(df,taxon.vector){
  df <- df[,colnames(df) != "category"]
  relab <- rowMedians(as.matrix(df))
  return(relab[match(taxon.vector,rownames(df))])
}

 
rank.df <- as.data.frame(matrix(nrow=length(unique(c(rownames(taxonomy_relab.species.metagenome.df),
                                                     rownames(taxonomy_relab.species.metatranscriptome.df)))),
                                ncol=9))
rownames(rank.df) <- unique(c(rownames(taxonomy_relab.species.metagenome.df),
                              rownames(taxonomy_relab.species.metatranscriptome.df)))
colnames(rank.df) <- c("category",
                       "bal_metagenome_rank","ua_metagenome_rank",
                       "bal_metatranscriptome_rank","ua_metatranscriptome_rank",
                       "bal_metagenome_median_relab","ua_metagenome_median_relab",
                       "bal_metatranscriptome_median_relab","ua_metatranscriptome_median_relab")

categories <- unique(as.data.frame(cbind(c(rownames(taxonomy_relab.species.metagenome.df),
                                           rownames(taxonomy_relab.species.metatranscriptome.df)),
                                         c(taxonomy_relab.species.metagenome.df$category,
                                           taxonomy_relab.species.metatranscriptome.df$category))))

rank.df[,1] <- categories[match(rownames(rank.df),categories[,1]),2]

rank.df[,2] <- find_rank_by_median(taxonomy_relab.species.metagenome.df[,colnames(taxonomy_relab.species.metagenome.df) %in% rownames(metagenome.map)[metagenome.map$Sample.Type == "BAL"]],
                                   rownames(rank.df))
rank.df[,3] <- find_rank_by_median(taxonomy_relab.species.metagenome.df[,colnames(taxonomy_relab.species.metagenome.df) %in% rownames(metagenome.map)[metagenome.map$Sample.Type == "UA"]],
                                   rownames(rank.df))
rank.df[,4] <- find_rank_by_median(taxonomy_relab.species.metatranscriptome.df[,colnames(taxonomy_relab.species.metatranscriptome.df) %in% rownames(metatranscriptome.map)[metatranscriptome.map$Sample.Type == "BAL"]],
                                   rownames(rank.df))
rank.df[,5] <- find_rank_by_median(taxonomy_relab.species.metatranscriptome.df[,colnames(taxonomy_relab.species.metatranscriptome.df) %in% rownames(metatranscriptome.map)[metatranscriptome.map$Sample.Type == "UA"]],
                                   rownames(rank.df))

rank.df[,6] <- find_median_relab(taxonomy_relab.species.metagenome.df[,colnames(taxonomy_relab.species.metagenome.df) %in% rownames(metagenome.map)[metagenome.map$Sample.Type == "BAL"]],
                                   rownames(rank.df))
rank.df[,7] <- find_median_relab(taxonomy_relab.species.metagenome.df[,colnames(taxonomy_relab.species.metagenome.df) %in% rownames(metagenome.map)[metagenome.map$Sample.Type == "UA"]],
                                   rownames(rank.df))
rank.df[,8] <- find_median_relab(taxonomy_relab.species.metatranscriptome.df[,colnames(taxonomy_relab.species.metatranscriptome.df) %in% rownames(metatranscriptome.map)[metatranscriptome.map$Sample.Type == "BAL"]],
                                   rownames(rank.df))
rank.df[,9] <- find_median_relab(taxonomy_relab.species.metatranscriptome.df[,colnames(taxonomy_relab.species.metatranscriptome.df) %in% rownames(metatranscriptome.map)[metatranscriptome.map$Sample.Type == "UA"]],
                                   rownames(rank.df))

write.table(rank.df,
            paste0(WORKING.DIR,"/data_frames/relab_ranking.tsv"),
            row.names = T,
            col.names = T,
            quote = F,
            sep = "\t")
```

### Writes the top 50 taxa for BAL and UA samples separately

```{r}
order.df <- taxonomy_relab.species.metagenome.df[,1:ncol(taxonomy_relab.species.metagenome.df)-1]
order.df <- order.df[order(rowMedians(as.matrix(order.df)),decreasing=T),]
  
write.table(rownames(order.df)[1:50],
            paste0(WORKING.DIR,"/data_frames/top_50_taxa.bal.metagenome.list"),
            row.names = F,
            col.names = F,
            quote = F,
            sep = "\t")

order.df <- taxonomy_relab.species.metatranscriptome.df[,1:ncol(taxonomy_relab.species.metatranscriptome.df)-1]
order.df <- order.df[order(rowMedians(as.matrix(order.df)),decreasing=T),]
  
write.table(rownames(order.df)[1:50],
            paste0(WORKING.DIR,"/data_frames/top_50_taxa.bal.metatranscriptome.list"),
            row.names = F,
            col.names = F,
            quote = F,
            sep = "\t")
```

### Create ordering data frames

```{r}
order.bal.metagenome <- taxonomy_relab.species.metagenome.df[,colnames(taxonomy_relab.species.metagenome.df) %in% rownames(metagenome.map)[metagenome.map$Sample.Type == "BAL"]]
order.ua.metagenome <- taxonomy_relab.species.metagenome.df[,colnames(taxonomy_relab.species.metagenome.df) %in% rownames(metagenome.map)[metagenome.map$Sample.Type == "UA"]]
order.bal.metatranscriptome <- taxonomy_relab.species.metatranscriptome.df[,colnames(taxonomy_relab.species.metatranscriptome.df) %in% rownames(metatranscriptome.map)[metatranscriptome.map$Sample.Type == "BAL"]]
order.ua.metatranscriptome <- taxonomy_relab.species.metatranscriptome.df[,colnames(taxonomy_relab.species.metatranscriptome.df) %in% rownames(metatranscriptome.map)[metatranscriptome.map$Sample.Type == "UA"]]
```

### Plot OTU relative abundance box plots

category.levels <- c("RNAvertebratevirus","DNAvertebratevirus","DNAbacteriaarchaeavirus","bacteria","fungi")
color.levels <- ("#75d6fd","#ff7e78","#ff7e78","#ff2401","#d6fb78")

##### Top 50 Combined Taxa | Species | Metagenome | BAL Metagenome Ordered

```{r,fig.height=10,fig.width=12}
df <- taxonomy_relab.species.metagenome.df[match(rownames(order.bal.metagenome),rownames(taxonomy_relab.species.metagenome.df)),]

category.levels <- c("RNAvertebratevirus","DNAbacteriaarchaeavirus","bacteria","fungi")
color.levels <- c("#75d6fd","#ff7e78","#ff2401","#d6fb78")
category.subset <- df$category[match(rownames(order.bal.metagenome),rownames(df))]
xaxis.lim <- plot_otu_box_pullxlim(order.bal.metagenome,category.subset,category.levels)

relab.subset <- df[,colnames(df) %in% rownames(metagenome.map)[metagenome.map$Sample.Type == "BAL"]]
metagenome_bal.boxplot <- plot_otu_box_facetplot(relab.subset,order.bal.metagenome,category.subset,category.levels,color.levels,xaxis.lim)

relab.subset <- df[,colnames(df) %in% rownames(metagenome.map)[metagenome.map$Sample.Type == "UA"]]
metagenome_ua.boxplot <- plot_otu_box_facetplot(relab.subset,order.bal.metagenome,category.subset,category.levels,color.levels,xaxis.lim)


pdf(paste0(WORKING.DIR,"/plots/metagenome_taxa_boxplots.metagenomebalordered.pdf"),
    width=12, 
    height=10)
plot_grid(plotlist = list(metagenome_bal.boxplot[[3]],metagenome_ua.boxplot[[3]],
                          metagenome_bal.boxplot[[4]],metagenome_ua.boxplot[[4]]),
          ncol=2,
          align = 'v',
          rel_heights=c(3,6))
dev.off()

```
##### Top 50 Combined Taxa | Species | Metatranscriptome | BAL Metatranscriptome Ordered

```{r,fig.height=10,fig.width=12}
df <- taxonomy_relab.species.metatranscriptome.df[match(rownames(order.bal.metatranscriptome),rownames(taxonomy_relab.species.metatranscriptome.df)),]

category.levels <- c("RNAvertebratevirus","DNAbacteriaarchaeavirus","bacteria","fungi")
color.levels <- c("#75d6fd","#ff7e78","#ff2401","#d6fb78")
order.bal.metatranscriptome <- df[,colnames(df) %in% rownames(metatranscriptome.map)[metatranscriptome.map$Sample.Type == "BAL"]]
category.subset <- df$category[match(rownames(order.bal.metatranscriptome),rownames(df))]
xaxis.lim <- plot_otu_box_pullxlim(order.bal.metatranscriptome,category.subset,category.levels)

relab.subset <- df[,colnames(df) %in% rownames(metatranscriptome.map)[metatranscriptome.map$Sample.Type == "BAL"]]
metatranscriptome_bal.boxplot <- plot_otu_box_facetplot(relab.subset,order.bal.metatranscriptome,category.subset,category.levels,color.levels,xaxis.lim)

relab.subset <- df[,colnames(df) %in% rownames(metatranscriptome.map)[metatranscriptome.map$Sample.Type == "UA"]]
metatranscriptome_ua.boxplot <- plot_otu_box_facetplot(relab.subset,order.bal.metatranscriptome,category.subset,category.levels,color.levels,xaxis.lim)

pdf(paste0(WORKING.DIR,"/plots/metatranscriptome_taxa_boxplots.metatranscriptomebalordered.pdf"),
    width=12, 
    height=10)
plot_grid(plotlist = list(metatranscriptome_bal.boxplot[[1]],metatranscriptome_ua.boxplot[[1]],
                          metatranscriptome_bal.boxplot[[2]],metatranscriptome_ua.boxplot[[2]],
                          metatranscriptome_bal.boxplot[[3]],metatranscriptome_ua.boxplot[[3]],
                          metatranscriptome_bal.boxplot[[4]],metatranscriptome_ua.boxplot[[4]]),
          ncol=2,
          align = 'v',
          rel_heights=c(1.2,1.2,5,8))
dev.off()

```

# Deplete Kraken human-assigned reads from final FASTQ files

## Split Kraken outputs

##### Metagenome

```{bash, eval = F}
find "$WORKING_DIR"/kraken/metagenome -name "*kraken_output.txt" | grep -v sarscov2 | sort -n > "$WORKING_DIR"/kraken_output_split.list

echo -e 'KRAKEN_OUTPUT="$(sed -ne "$SGE_TASK_ID p" "$WORKING_DIR"/kraken_output_split.list)" \n
cat $KRAKEN_OUTPUT | split -l 20000000 -d - "$TEMP_DIR"/metagenome/"$(basename "$KRAKEN_OUTPUT")"' | qsub -v TEMP_DIR="$TEMP_DIR",WORKING_DIR="$WORKING_DIR" -wd "$WORKING_DIR"/stderrout -N split -t 1-"$(wc -l "$WORKING_DIR"/kraken_output_split.list | awk '{print $1}')"
```

##### Metatranscriptome

```{bash, eval = F}
find "$WORKING_DIR"/kraken/metatranscriptome -name "*kraken_output.txt" > "$WORKING_DIR"/kraken_output_split.list

echo -e 'KRAKEN_OUTPUT="$(sed -ne "$SGE_TASK_ID p" "$WORKING_DIR"/kraken_output_split.list)" \n
	cat $KRAKEN_OUTPUT | split -l 20000000 -d - "$TEMP_DIR"/metatranscriptome/"$(basename "$KRAKEN_OUTPUT")"' | qsub -v WORKING_DIR="$WORKING_DIR" -wd "$WORKING_DIR"/stderrout -N split -t 1-"$(wc -l "$WORKING_DIR"/kraken_output_split.list | awk '{print $1}')"
```

## Split FASTQ files

##### Metagenome

```{bash, eval = F}
find "$WORKING_DIR"/final_reads/metagenome/ -name "*fastq.gz" | sed "s/[.].*//g" | sort -n | uniq > "$WORKING_DIR"/metagenome_split.list

echo -e 'SAMPLE="$(sed -ne "$SGE_TASK_ID p" "$WORKING_DIR"/metagenome_split.list)" \n
	zcat "$SAMPLE".fastq.gz | split -l 120000000 -d - "$TEMP_DIR"/metagenome/"$(basename $SAMPLE)".fastq' | qsub -v TEMP_DIR="$TEMP_DIR",WORKING_DIR="$WORKING_DIR" -wd "$WORKING_DIR"/stderrout -N split -t 1-"$(wc -l "$WORKING_DIR"/metagenome_split.list | awk '{print $1}')"
```

##### Metatranscriptome

```{bash, eval = F}
find "$WORKING_DIR"/final_reads/metatranscriptome/ -name "*fastq.gz" | sed "s/[.].*//g" | sort -n | uniq > "$WORKING_DIR"/metatranscriptome_split.list

echo -e 'SAMPLE="$(sed -ne "$SGE_TASK_ID p" "$WORKING_DIR"/metatranscriptome_split.list)" \n
	zcat "$SAMPLE".fastq.gz | split -l 120000000 -d - "$TEMP_DIR"/metatranscriptome/"$(basename $SAMPLE)".fastq' | qsub -v WORKING_DIR="$WORKING_DIR" -wd "$WORKING_DIR"/stderrout -N split -t 1-"$(wc -l "$WORKING_DIR"/metatranscriptome_split.list | awk '{print $1}')"
```

## Create subset FASTQs without human-assigned reads

##### Metagenome

```{bash, eval = F}
module load KrakenTools/0.1-alpha
module load python/3.7.3-foss-2016b

for KRAKEN_OUTPUT in $(find "$TEMP_DIR"/metagenome/ -name "*kraken_output*" | grep -v sarscov2)
do
SAMPLE="$(basename "$KRAKEN_OUTPUT" | sed "s/[.]kraken_output.*//g")"

find "$TEMP_DIR"/metagenome/ | grep "$SAMPLE"[.]fastq | grep -v include | grep -v exclude > "$TEMP_DIR"/metagenome/"$SAMPLE"_fastq.list

echo -e 'FASTQ="$(sed -ne "$SGE_TASK_ID p" "$TEMP_DIR"/metagenome/"$SAMPLE"_fastq.list)" \n
~/scripts/KrakenTools/extract_kraken_reads.py --exclude --fastq-output -t 9606 \
--max 1000000000 \
-s "$FASTQ" \
-k "$KRAKEN_OUTPUT" \
-r "$WORKING_DIR"/kraken/metagenome/"$SAMPLE".kraken_report.txt \
-o "$(echo "$FASTQ" | sed "s/[.]fastq/.exclude.s_homo_sapiens.fastq/g")"_${KRAKEN_OUTPUT: -2}' | qsub -v SAMPLE="$SAMPLE",KRAKEN_OUTPUT="$KRAKEN_OUTPUT",TEMP_DIR="$TEMP_DIR",WORKING_DIR="$WORKING_DIR" -l mem_free=200G -wd "$WORKING_DIR"/stderrout -N extract_kraken_reads_nohuman -t 1-"$(wc -l "$TEMP_DIR"/metagenome/"$SAMPLE"_fastq.list | awk '{print $1}')"
done
```

##### Metatranscriptome

```{bash, eval = F}
module load KrakenTools/0.1-alpha
module load python/3.7.3-foss-2016b

for KRAKEN_OUTPUT in $(find "$TEMP_DIR"/metatranscriptome/ -name "*kraken_output*" | grep -v sarscov2)
do
SAMPLE="$(basename "$KRAKEN_OUTPUT" | sed "s/[.]kraken_output.*//g")"

find "$TEMP_DIR"/metatranscriptome/ | grep "$SAMPLE"[.]fastq | grep -v include | grep -v exclude > "$TEMP_DIR"/metatranscriptome/"$SAMPLE"_fastq.list

echo -e 'FASTQ="$(sed -ne "$SGE_TASK_ID p" "$TEMP_DIR"/metatranscriptome/"$SAMPLE"_fastq.list)" \n
~/scripts/KrakenTools/extract_kraken_reads.py --exclude --fastq-output -t 9606 \
--max 1000000000 \
-s "$FASTQ" \
-k "$KRAKEN_OUTPUT" \
-r "$WORKING_DIR"/kraken/metatranscriptome/"$SAMPLE".kraken_report.txt \
-o "$(echo "$FASTQ" | sed "s/[.]fastq/.exclude.s_homo_sapiens.fastq/g")"_${KRAKEN_OUTPUT: -2}' | qsub -v SAMPLE="$SAMPLE",KRAKEN_OUTPUT="$KRAKEN_OUTPUT",TEMP_DIR="$TEMP_DIR",WORKING_DIR="$WORKING_DIR" -l mem_free=200G -wd "$WORKING_DIR"/stderrout -N extract_kraken_reads_nohuman -t 1-"$(wc -l "$TEMP_DIR"/metatranscriptome/"$SAMPLE"_fastq.list | awk '{print $1}')"
done
```

## Combine human-depleted FASTQs for each sample

##### Metagenome

```{bash, eval = F}
for SAMPLE in $(find "$TEMP_DIR"/metagenome/ -name "*.exclude.s_homo_sapiens.fastq[0-9]*" | sed "s/[.].*//g" | sort -n | uniq)
do
	for i in $(ls "$SAMPLE"* | grep ".exclude.s_homo_sapiens.fastq[0-9]*" | sort -n)
	do
		cat $i >> "$WORKING_DIR"/kraken/metagenome/"$(basename "$SAMPLE")".exclude.s_homo_sapiens.fastq
	done
done
```

##### Metatranscriptome

```{bash, eval = F}
for SAMPLE in $(find "$TEMP_DIR"/metatranscriptome/ -name "*.exclude.s_homo_sapiens.fastq[0-9]*" | sed "s/[.].*//g" | sort -n | uniq)
do
	for i in $(ls "$SAMPLE"* | grep ".exclude.s_homo_sapiens.fastq[0-9]*" | sort -n)
	do
		cat $i >> "$WORKING_DIR"/kraken/metatranscriptome/"$(basename "$SAMPLE")".exclude.s_homo_sapiens.fastq
	done
done
```

## Compress combined human-depleted FASTQs
```{bash, eval = F}
find "$WORKING_DIR"/kraken/ -name "*fastq" > "$WORKING_DIR"/gzip.list

echo -e 'FASTQ="$(sed -ne "$SGE_TASK_ID p" "$WORKING_DIR"/gzip.list)" \n
gzip -f $FASTQ' | qsub -v WORKING_DIR="$WORKING_DIR" -N gzip -wd "$WORKING_DIR"/stderrout/ -t 1-"$(wc -l "$WORKING_DIR"/gzip.list | awk '{print $1}')"

```

# Conduct functional profiling analysis of metagenome and metatranscriptome sequences

## Conduct BLAST searches of human-depleted reads against FMAP KO reference database

##### Metagenome

```{bash, eval = F}
module load fmap/0.15-goolf-1.7.20

find "$WORKING_DIR"/kraken/metagenome -name "*.exclude.s_homo_sapiens.fastq.gz" | sed "s/[.].*//g" | sed "s/.*\\///g" | sort -n > "$WORKING_DIR"/fmap_metagenome.list

echo -e 'SAMPLE="$(sed -ne "$SGE_TASK_ID p" "$WORKING_DIR"/fmap_metagenome.list)" \n
perl "$FMAP_BIN_DIR"/FMAP_mapping.pl \
	-p "$THREADS" \
	-t "$TEMP_DIR"/fmap/metagenome/ \
	"$WORKING_DIR"/kraken/metagenome/"$SAMPLE".exclude.s_homo_sapiens.fastq.gz > "$WORKING_DIR"/fmap/metagenome/"$SAMPLE".blastx_hits.txt' | qsub -pe threaded "$THREADS" -l mem_free=20G -v FMAP_BIN_DIR="$FMAP_BIN_DIR",THREADS="$THREADS",TEMP_DIR="$TEMP_DIR",WORKING_DIR="$WORKING_DIR" -N fmap_metagenome -wd "$WORKING_DIR"/stderrout/ -t 1-"$(wc -l "$WORKING_DIR"/fmap_metagenome.list | awk '{print $1}')"
```

##### Metatranscriptome

```{bash, eval = F}
module load fmap/0.15-goolf-1.7.20

find "$WORKING_DIR"/kraken/metatranscriptome -name "*.exclude.s_homo_sapiens.fastq.gz" | sed "s/[.].*//g" | sed "s/.*\\///g" | sort -n > "$WORKING_DIR"/fmap_metatranscriptome.list

echo -e 'SAMPLE="$(sed -ne "$SGE_TASK_ID p" "$WORKING_DIR"/fmap_metatranscriptome.list)" \n
perl "$FMAP_BIN_DIR"/FMAP_mapping.pl \
	-p "$THREADS" \
	-t "$TEMP_DIR"/fmap/metatranscriptome/ \
	"$WORKING_DIR"/kraken/metatranscriptome/"$SAMPLE".exclude.s_homo_sapiens.fastq.gz > "$WORKING_DIR"/fmap/metatranscriptome/"$SAMPLE".blastx_hits.txt' | qsub -pe threaded "$THREADS" -l mem_free=20G -v FMAP_BIN_DIR="$FMAP_BIN_DIR",THREADS="$THREADS",TEMP_DIR="$TEMP_DIR",WORKING_DIR="$WORKING_DIR" -N fmap_metatranscriptome -wd "$WORKING_DIR"/stderrout/ -t 1-"$(wc -l "$WORKING_DIR"/fmap_metatranscriptome.list | awk '{print $1}')"
```

## Derive count values from FMAP BLASTX searches

### Split FMAP BLASTX outputs

##### Metagenome

```{bash, eval = F}
for FMAP_MAPPING_OUTPUT in $(find "$WORKING_DIR"/fmap/metagenome -name "*blastx_hits.txt")
do
split -l 10000000 -d "$FMAP_MAPPING_OUTPUT" "$FMAP_MAPPING_OUTPUT"
done
```

##### Metatranscriptome

```{bash, eval = F}
for FMAP_MAPPING_OUTPUT in $(find "$WORKING_DIR"/fmap/metatranscriptome -name "*blastx_hits.txt")
do
split -l 10000000 -d "$FMAP_MAPPING_OUTPUT" "$FMAP_MAPPING_OUTPUT"
done
```

### Create abundance dataframes for KO terms

##### Metagenome

```{bash, eval = F}
module load fmap/0.15-goolf-1.7.20

find "$WORKING_DIR"/fmap/metagenome -name "*blastx_hits.txt[0-9]*" > "$WORKING_DIR"/fmap_metagenome.list

echo -e 'FMAP_MAPPING_OUTPUT="$(sed -ne "$SGE_TASK_ID p" "$WORKING_DIR"/fmap_metagenome.list)" \n
perl "$FMAP_BIN_DIR"/FMAP_quantification.pl "$FMAP_MAPPING_OUTPUT" > "$(echo "$FMAP_MAPPING_OUTPUT" | sed "s/[.]blastx_hits.txt/.abundance.txt/g")"' | qsub -l mem_free=20G -v FMAP_BIN_DIR="$FMAP_BIN_DIR",WORKING_DIR="$WORKING_DIR" -N fmap_metagenome -wd "$WORKING_DIR"/stderrout/ -t 1-"$(wc -l "$WORKING_DIR"/fmap_metagenome.list | awk '{print $1}')"
```

##### Metatranscriptome

```{bash, eval = F}
module load fmap/0.15-goolf-1.7.20

find "$WORKING_DIR"/fmap/metatranscriptome -name "*blastx_hits.txt[0-9]*" > "$WORKING_DIR"/fmap_metatranscriptome.list

echo -e 'FMAP_MAPPING_OUTPUT="$(sed -ne "$SGE_TASK_ID p" "$WORKING_DIR"/fmap_metatranscriptome.list)" \n
perl "$FMAP_BIN_DIR"/FMAP_quantification.pl "$FMAP_MAPPING_OUTPUT" > "$(echo "$FMAP_MAPPING_OUTPUT" | sed "s/[.]blastx_hits.txt/.abundance.txt/g")"' | qsub -l mem_free=20G -v FMAP_BIN_DIR="$FMAP_BIN_DIR",WORKING_DIR="$WORKING_DIR" -N fmap_metatranscriptome -wd "$WORKING_DIR"/stderrout/ -t 1-"$(wc -l "$WORKING_DIR"/fmap_metatranscriptome.list | awk '{print $1}')"
```

## Construct KO counts data frames from FMAP output

### Set R inputs

```{R}
WORKING.DIR <- ""
```

### Load packages and view sessionInfo

```{r}
library(KEGGREST)

sessionInfo()
```
```{R, eval = F}
R version 4.0.2 (2020-06-22)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 18363)

Matrix products: default

locale:
[1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252    LC_MONETARY=English_United States.1252
[4] LC_NUMERIC=C                           LC_TIME=English_United States.1252    

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
[1] KEGGREST_1.28.0

loaded via a namespace (and not attached):
 [1] IRanges_2.22.2      png_0.1-7           Biostrings_2.56.0   digest_0.6.25       crayon_1.3.4        R6_2.4.1            stats4_4.0.2       
 [8] evaluate_0.14       httr_1.4.2          zlibbioc_1.34.0     rlang_0.4.7         XVector_0.28.0      rstudioapi_0.11     S4Vectors_0.26.1   
[15] rmarkdown_2.4       tools_4.0.2         tinytex_0.26        parallel_4.0.2      xfun_0.17           yaml_2.2.1          compiler_4.0.2     
[22] BiocGenerics_0.34.0 htmltools_0.5.0     knitr_1.30
```

### Read in FMAP abundance files

```{R}
abundance.files <- list.files(paste0(WORKING.DIR,"/fmap/"),
                              pattern="abundance.txt[0-9].*",
                              recursive = T)
metagenome.abundance.files <- grep("metagenome",abundance.files,value = T)
metatranscriptome.abundance.files <- grep("metatranscriptome",abundance.files,value = T)
```

### Compile list of KO terms and identifiers quantified across all FMAP abundance files

```{R}
ko_ids <- c()
ko_descriptors <- c()
for(i in 1:length(abundance.files)){
  abundance.file <- read.delim(paste0(WORKING.DIR,"/fmap/",abundance.files[i]),
                               quote="")
  ko_ids <- unique(c(ko_ids,abundance.file[,1]))
  ko_descriptors <- unique(c(ko_descriptors,paste0(abundance.file[,1],": ",abundance.file[,2])))
}
ko_ids <- sort(ko_ids) 
ko_descriptors <- sort(ko_descriptors)
```

### Compile counts for identified KOs for all samples

##### Metagenome

```{R}
metagenome.counts <- as.data.frame(matrix(0,
                                          nrow=length(ko_ids),
                                          ncol=length(unique(gsub("[.].*","",metagenome.abundance.files)))))
rownames(metagenome.counts) <- ko_ids
colnames(metagenome.counts) <- unique(basename(gsub("[.].*","",metagenome.abundance.files)))

for(i in 1:length(metagenome.abundance.files)){
  sample <- basename(gsub("[.].*","",metagenome.abundance.files[i]))
  abundance.file <- read.delim(paste0(WORKING.DIR,"/fmap/",metagenome.abundance.files[i]),
                               quote="")
  counts.subset <- metagenome.counts[,colnames(metagenome.counts) == sample] + abundance.file$count[match(rownames(metagenome.counts),abundance.file[,1]),drop = F]
  counts.subset[is.na(counts.subset)] <- 0
  metagenome.counts[,colnames(metagenome.counts) == sample] <- metagenome.counts[,colnames(metagenome.counts) == sample] + counts.subset
}
```

##### Metatranscriptome

```{R}
metatranscriptome.counts <- as.data.frame(matrix(0,
                                          nrow=length(ko_ids),
                                          ncol=length(unique(gsub("[.].*","",metatranscriptome.abundance.files)))))
rownames(metatranscriptome.counts) <- ko_ids
colnames(metatranscriptome.counts) <- unique(basename(gsub("[.].*","",metatranscriptome.abundance.files)))

for(i in 1:length(metatranscriptome.abundance.files)){
  sample <- basename(gsub("[.].*","",metatranscriptome.abundance.files[i]))
  abundance.file <- read.delim(paste0(WORKING.DIR,"/fmap/",metatranscriptome.abundance.files[i]),
                               quote="")
  counts.subset <- metatranscriptome.counts[,colnames(metatranscriptome.counts) == sample] + abundance.file$count[match(rownames(metatranscriptome.counts),abundance.file[,1]),drop = F]
  counts.subset[is.na(counts.subset)] <- 0
  metatranscriptome.counts[,colnames(metatranscriptome.counts) == sample] <- metatranscriptome.counts[,colnames(metatranscriptome.counts) == sample] + counts.subset
}
```

### Identify pathway KOs for each KO term identified by FMAP

```{r}
get_pathway_kos <- function(base_ko.vector){
  pathway_ko.list <- list()
  while(length(base_ko.vector) > 0){
    if(length(base_ko.vector)){
      query <- keggGet(base_ko.vector[1:10])
    }else{
      query <- keggGet(base_ko.vector)
    }
    for(i in 1:length(query)){
      for(j in 1:length(query[[i]]$PATHWAY)){
        pathway_ko.list[paste0(names(query[[i]]$PATHWAY[j]),"|",query[[i]]$PATHWAY[j])] <- paste(as.character(pathway_ko.list[paste0(names(query[[i]]$PATHWAY[j]),"|",query[[i]]$PATHWAY[j])]),query[[i]]$ENTRY,sep="|")
      }
    }
    if(length(base_ko.vector) > 10){
      base_ko.vector <- base_ko.vector[11:length(base_ko.vector)]
    }else{
      base_ko.vector <- c()
    }
    print(length(base_ko.vector))
  }
  return(pathway_ko.list)
}

pathway_ko.list <- get_pathway_kos(ko_ids)
pathway_ko.list <- lapply(pathway_ko.list,function(x){unique(unlist(strsplit(gsub("NULL[|]","",x),split="[|]")))})
pathway_ko.list <- pathway_ko.list[grep("ko",names(pathway_ko.list))]
```

### Construct table of KO pathway counts

##### Metagenome

```{r}
metagenome.pathwayko.counts <- as.data.frame(matrix(nrow=length(pathway_ko.list),
                                                    ncol=ncol(metagenome.counts)))
rownames(metagenome.pathwayko.counts) <- names(pathway_ko.list)
colnames(metagenome.pathwayko.counts) <- colnames(metagenome.counts)

for(i in 1:nrow(metagenome.pathwayko.counts)){
  counts.subset <- metagenome.counts[rownames(metagenome.counts) %in% pathway_ko.list[[i]],]
  metagenome.pathwayko.counts[i,] <- colSums(counts.subset)
}
```

##### Metatranscriptome

```{r}
metatranscriptome.pathwayko.counts <- as.data.frame(matrix(nrow=length(pathway_ko.list),
                                                    ncol=ncol(metatranscriptome.counts)))
rownames(metatranscriptome.pathwayko.counts) <- names(pathway_ko.list)
colnames(metatranscriptome.pathwayko.counts) <- colnames(metatranscriptome.counts)

for(i in 1:nrow(metatranscriptome.pathwayko.counts)){
  counts.subset <- metatranscriptome.counts[rownames(metatranscriptome.counts) %in% pathway_ko.list[[i]],]
  metatranscriptome.pathwayko.counts[i,] <- colSums(counts.subset)
}
```

### Append rownames of KO terms counts data frame with descriptors

```{r}
rownames(metagenome.counts) <- ko_descriptors
rownames(metatranscriptome.counts) <- ko_descriptors
```

### Output counts data frames for KO terms and pathway KO terms

##### Metagenome

```{R}
write.table(metagenome.counts,
            paste0(WORKING.DIR,"/data_frames/counts_ko_metagenome.tsv"),
            row.names = T,
            col.names = T,
            quote = F,
            sep = "\t")

write.table(metagenome.pathwayko.counts,
            paste0(WORKING.DIR,"/data_frames/counts_pathwayko_metagenome.tsv"),
            row.names = T,
            col.names = T,
            quote = F,
            sep = "\t")
```

##### Metatranscriptome

```{R}
write.table(metatranscriptome.counts,
            paste0(WORKING.DIR,"/data_frames/counts_ko_metatranscriptome.tsv"),
            row.names = T,
            col.names = T,
            quote = F,
            sep = "\t")

write.table(metatranscriptome.pathwayko.counts,
            paste0(WORKING.DIR,"/data_frames/counts_pathwayko_metatranscriptome.tsv"),
            row.names = T,
            col.names = T,
            quote = F,
            sep = "\t")
```

# Conduct resistome analysis of metagenome and metatranscriptome sequences

## Download MEGAres database of antibiotic resistance genes

```{bash, eval = F}
wget https://urldefense.proofpoint.com/v2/url?u=https-3A__megares.meglab.org_download_megares-5Fv2.00_megares-5Ffull-5Fdatabase-5Fv2.00.fasta&d=DwIGAg&c=j5oPpO0eBH1iio48DtsedeElZfc04rx3ExJHeIIZuCs&r=ufTTecG4vkmes8mF_HId8gR_hbrb8-dKDFW3-voxs2A&m=9Ldy7Wd7XZ7NQzW7DY5909D3ElJvgX4saq7j5jlDI94&s=wlo0CgSwNHfFRJ3Ofae5alp46TIYWO2JCxyVfwVSocQ&e=  -O "$WORKING_DIR"/references/megares_full_database_v2.00.fasta
```

## Quantify antibiotic resistance genes using Salmon

### Construct index for reference MEGAres database

```{bash, eval = F}
module load salmon/1.3.0

salmon index --keepDuplicates -t "$WORKING_DIR"/references/megares_full_database_v2.00.fasta -i "$WORKING_DIR"/references/megares_full_database_v2.00.fasta.salmon.index 
```

### Run Salmon to quantify genes

##### Metagenome

```{bash, eval = F}
module load salmon/1.3.0

find "$WORKING_DIR"/final_reads/metagenome/ -name "*fastq.gz" | sed "s/[.].*//g" | sort -n | uniq > "$WORKING_DIR"/salmon.list

echo -e 'SAMPLE="$(sed -ne "$SGE_TASK_ID p" "$WORKING_DIR"/salmon.list)" \n
salmon quant -i "$WORKING_DIR"/references/megares_full_database_v2.00.fasta.salmon.index  \
--libType A \
-1 "$SAMPLE".1.fastq.gz \
-2 "$SAMPLE".2.fastq.gz \
-o "$WORKING_DIR"/antibiotic_resistance/salmon/metagenome"$(basename "$SAMPLE")" \
--allowDovetail --meta' | qsub -v WORKING_DIR="$WORKING_DIR" -wd "$WORKING_DIR"/stderrout -N salmon -t 1-"$(wc -l "$WORKING_DIR"/salmon.list | awk '{print $1}')"
```

##### Metatranscriptome

```{bash, eval = F}
module load salmon/1.3.0

find "$WORKING_DIR"/final_reads/metatranscriptome/ -name "*fastq.gz" | sed "s/[.].*//g" | sort -n | uniq > "$WORKING_DIR"/salmon.list

echo -e 'SAMPLE="$(sed -ne "$SGE_TASK_ID p" "$WORKING_DIR"/salmon.list)" \n
salmon quant -i "$WORKING_DIR"/references/megares_full_database_v2.00.fasta.salmon.index  \
--libType A \
-1 "$SAMPLE".1.fastq.gz \
-2 "$SAMPLE".2.fastq.gz \
-o "$WORKING_DIR"/antibiotic_resistance/salmon/metatranscriptome"$(basename "$SAMPLE")" \
--allowDovetail --meta' | qsub -v WORKING_DIR="$WORKING_DIR" -wd "$WORKING_DIR"/stderrout -N salmon -t 1-"$(wc -l "$WORKING_DIR"/salmon.list | awk '{print $1}')"
```

## Construct resistome count tables

### Set R inputs

```{R}
WORKING.DIR <- ""
```

### Load packages and view sessionInfo

```{R}
library(DESeq2)
library(cowplot)
library(ggplot2)
library(gridExtra)
library(see)

sessionInfo()
```

```{R, eval = F}
R version 4.0.2 (2020-06-22)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 18363)

Matrix products: default

locale:
[1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252    LC_MONETARY=English_United States.1252
[4] LC_NUMERIC=C                           LC_TIME=English_United States.1252    

attached base packages:
[1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] see_0.6.1                   gridExtra_2.3               ggplot2_3.3.3               cowplot_1.1.1              
 [5] DESeq2_1.28.1               SummarizedExperiment_1.18.2 DelayedArray_0.14.1         matrixStats_0.57.0         
 [9] Biobase_2.48.0              GenomicRanges_1.40.0        GenomeInfoDb_1.24.2         IRanges_2.22.2             
[13] S4Vectors_0.26.1            BiocGenerics_0.34.0        

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.5             locfit_1.5-9.4         lattice_0.20-41        digest_0.6.25          plyr_1.8.6            
 [6] R6_2.5.0               ggridges_0.5.2         RSQLite_2.2.1          pillar_1.4.7           zlibbioc_1.34.0       
[11] rlang_0.4.7            rstudioapi_0.13        annotate_1.66.0        blob_1.2.1             Matrix_1.2-18         
[16] effectsize_0.4.1       splines_4.0.2          BiocParallel_1.22.0    geneplotter_1.66.0     RCurl_1.98-1.2        
[21] bit_4.0.4              munsell_0.5.0          compiler_4.0.2         xfun_0.17              parameters_0.10.1     
[26] pkgconfig_2.0.3        insight_0.11.1         tidyselect_1.1.0       tibble_3.0.4           GenomeInfoDbData_1.2.3
[31] XML_3.99-0.5           crayon_1.3.4           dplyr_1.0.2            withr_2.3.0            bitops_1.0-6          
[36] grid_4.0.2             xtable_1.8-4           gtable_0.3.0           lifecycle_0.2.0        DBI_1.1.0             
[41] magrittr_2.0.1         bayestestR_0.8.0       scales_1.1.1           XVector_0.28.0         genefilter_1.70.0     
[46] ellipsis_0.3.1         generics_0.1.0         vctrs_0.3.6            RColorBrewer_1.1-2     tools_4.0.2           
[51] bit64_4.0.5            glue_1.4.2             purrr_0.3.4            survival_3.1-12        AnnotationDbi_1.50.3  
[56] colorspace_2.0-0       memoise_1.1.0          knitr_1.30            
```

### Set list of samples to exclude due to IRB issues

```{r}
sample_exclude.list <- c("UCS_0045","UCS_0052","UCS_0170","UCS_0171")
```

### Read input count files of quantified MEGAres genes

##### Metagenome

```{R}
salmon.metagenome.files <- list.files(paste0(WORKING.DIR,"/antibiotic_resistance/salmon/metagenome/"),
                                      recursive = T,
                                      pattern = "quant.sf")

rawcounts.metagenome.df <- as.data.frame(matrix(nrow=nrow(read.delim(paste0(WORKING.DIR,"/antibiotic_resistance/salmon/metagenome/",salmon.metagenome.files[1]),
                                                                  header = T)),
                                  ncol=length(salmon.metagenome.files)))
rownames(rawcounts.metagenome.df) <- read.delim(paste0(WORKING.DIR,"/antibiotic_resistance/salmon/metagenome/",salmon.metagenome.files[1]))[,1]
colnames(rawcounts.metagenome.df) <- dirname(salmon.metagenome.files)
for(i in 1:ncol(rawcounts.metagenome.df)){
  salmon.output <- read.delim(paste0(WORKING.DIR,"/antibiotic_resistance/salmon/metagenome/",salmon.metagenome.files[i]),
                              header = T)
  rawcounts.metagenome.df[,i] <- salmon.output[,5]
}

rawcounts.metagenome.df <- rawcounts.metagenome.df[,!(gsub("_[BU].*","",colnames(rawcounts.metagenome.df)) %in% sample_exclude.list)]
rawcounts.metagenome.bal.df <- rawcounts.metagenome.df[,grep("BAL",colnames(rawcounts.metagenome.df))]
```

##### Metatranscriptome

```{R}
salmon.metatranscriptome.files <- list.files(paste0(WORKING.DIR,"/antibiotic_resistance/salmon/metatranscriptome/"),
                                      recursive = T,
                                      pattern = "quant.sf")

rawcounts.metatranscriptome.df <- as.data.frame(matrix(nrow=nrow(read.delim(paste0(WORKING.DIR,"/antibiotic_resistance/salmon/metatranscriptome/",salmon.metatranscriptome.files[1]),
                                                                  header = T)),
                                  ncol=length(salmon.metatranscriptome.files)))
rownames(rawcounts.metatranscriptome.df) <- read.delim(paste0(WORKING.DIR,"/antibiotic_resistance/salmon/metatranscriptome/",salmon.metatranscriptome.files[1]))[,1]
colnames(rawcounts.metatranscriptome.df) <- dirname(salmon.metatranscriptome.files)
for(i in 1:ncol(rawcounts.metatranscriptome.df)){
  salmon.output <- read.delim(paste0(WORKING.DIR,"/antibiotic_resistance/salmon/metatranscriptome/",salmon.metatranscriptome.files[i]),
                              header = T)
  rawcounts.metatranscriptome.df[,i] <- salmon.output[,5]
}

rawcounts.metatranscriptome.df <- rawcounts.metatranscriptome.df[,!(gsub("_[BU].*","",colnames(rawcounts.metatranscriptome.df)) %in% sample_exclude.list)]
rawcounts.metatranscriptome.bal.df <- rawcounts.metatranscriptome.df[,grep("BAL",colnames(rawcounts.metatranscriptome.df))]
```

### Filter out MEGAres genes so that only genes that actively confer resistance are kept

```{r}
rawcounts.metagenome.df <- rawcounts.metagenome.df[!(grepl("RequiresSNPConfirmation",rownames(rawcounts.metagenome.df))),]
rawcounts.metagenome.df <- rawcounts.metagenome.df[!(grepl("resistant",rownames(rawcounts.metagenome.df),ignore.case = T)),]

rawcounts.metatranscriptome.df <- rawcounts.metatranscriptome.df[match(rownames(rawcounts.metagenome.df),rownames(rawcounts.metatranscriptome.df)),]
```

### Simplify MEGAres counts dataframes

#### Collapse genes by orthologous gene families

```{r}
antibiotic_genes <- unique(paste(sapply(strsplit(rownames(rawcounts.metagenome.df), "[|]"), "[", 3),
                                 sapply(strsplit(rownames(rawcounts.metagenome.df), "[|]"), "[", 4),
                                 sapply(strsplit(rownames(rawcounts.metagenome.df), "[|]"), "[", 5),
                                 sep="|"))

counts.metagenome.df <- as.data.frame(matrix(nrow=length(antibiotic_genes),
                                             ncol=ncol(rawcounts.metagenome.df)))
rownames(counts.metagenome.df) <- antibiotic_genes
colnames(counts.metagenome.df) <- colnames(rawcounts.metagenome.df)

for(i in 1:nrow(counts.metagenome.df)){
   counts.metagenome.df[i,] <- colSums(rawcounts.metagenome.df[grep(paste0(gsub("[|]",".",rownames(counts.metagenome.df)[i]),"$"),
                                                                    rownames(rawcounts.metagenome.df)),])
}

counts.metatranscriptome.df <- as.data.frame(matrix(nrow=length(antibiotic_genes),
                                             ncol=ncol(rawcounts.metatranscriptome.df)))
rownames(counts.metatranscriptome.df) <- antibiotic_genes
colnames(counts.metatranscriptome.df) <- colnames(rawcounts.metatranscriptome.df)

for(i in 1:nrow(counts.metatranscriptome.df)){
   counts.metatranscriptome.df[i,] <- colSums(rawcounts.metatranscriptome.df[grep(paste0(gsub("[|]",".",rownames(counts.metatranscriptome.df)[i]),"$"),
                                                                    rownames(rawcounts.metatranscriptome.df)),])
}

counts.metagenome.bal.df <- counts.metagenome.df[,grep("BAL",colnames(counts.metagenome.df))]
counts.metatranscriptome.bal.df <- counts.metatranscriptome.df[,grep("BAL",colnames(counts.metatranscriptome.df))]
```

#### Collapse genes by antibiotic classes

```{r}
antibiotic_classes <- unique(paste(sapply(strsplit(rownames(rawcounts.metagenome.df), "[|]"), "[", 3)))
                                
simplifiedcounts.metagenome.df <- as.data.frame(matrix(nrow=length(antibiotic_classes),
                                                       ncol=ncol(rawcounts.metagenome.df)))
rownames(simplifiedcounts.metagenome.df) <- antibiotic_classes
colnames(simplifiedcounts.metagenome.df) <- colnames(rawcounts.metagenome.df)

for(i in 1:nrow(simplifiedcounts.metagenome.df)){
   simplifiedcounts.metagenome.df[i,] <- colSums(rawcounts.metagenome.df[grep(paste0("[|]",rownames(simplifiedcounts.metagenome.df)[i],"[|]"),
                                     rownames(rawcounts.metagenome.df)),])
}

simplifiedcounts.metatranscriptome.df <- as.data.frame(matrix(nrow=length(antibiotic_classes),
                                                       ncol=ncol(rawcounts.metatranscriptome.df)))
rownames(simplifiedcounts.metatranscriptome.df) <- antibiotic_classes
colnames(simplifiedcounts.metatranscriptome.df) <- colnames(rawcounts.metatranscriptome.df)

for(i in 1:nrow(simplifiedcounts.metatranscriptome.df)){
   simplifiedcounts.metatranscriptome.df[i,] <- colSums(rawcounts.metatranscriptome.df[grep(paste0("[|]",rownames(simplifiedcounts.metatranscriptome.df)[i],"[|]"),
                                     rownames(rawcounts.metatranscriptome.df)),])
}
simplifiedcounts.metagenome.bal.df <- simplifiedcounts.metagenome.df[,grep("BAL",colnames(simplifiedcounts.metagenome.df))]
simplifiedcounts.metatranscriptome.bal.df <- simplifiedcounts.metatranscriptome.df[,grep("BAL",colnames(simplifiedcounts.metatranscriptome.df))]
```

### Write output MEGAres quantified antibiotic resistance gene count tables

```{r}
write.table(rawcounts.metagenome.df,
            paste0(WORKING.DIR,"/data_frames/rawcounts_megares_metagenome.tsv"),
            row.names = T,
            col.names = T,
            quote = F,
            sep = "\t")

write.table(counts.metatranscriptome.df,
            paste0(WORKING.DIR,"/data_frames/rawcounts_megares_metatranscriptome.tsv"),
            row.names = T,
            col.names = T,
            quote = F,
            sep = "\t")

write.table(counts.metagenome.df,
            paste0(WORKING.DIR,"/data_frames/counts_megares_metagenome.tsv"),
            row.names = T,
            col.names = T,
            quote = F,
            sep = "\t")

write.table(counts.metatranscriptome.df,
            paste0(WORKING.DIR,"/data_frames/counts_megares_metatranscriptome.tsv"),
            row.names = T,
            col.names = T,
            quote = F,
            sep = "\t")

write.table(simplifiedcounts.metagenome.df,
            paste0(WORKING.DIR,"/data_frames/simplifiedcounts_megares_metagenome.tsv"),
            row.names = T,
            col.names = T,
            quote = F,
            sep = "\t")

write.table(simplifiedcounts.metatranscriptome.df,
            paste0(WORKING.DIR,"/data_frames/simplifiedcounts_megares_metatranscriptome.tsv"),
            row.names = T,
            col.names = T,
            quote = F,
            sep = "\t")
```

# Characterize CRISPR arrays in metagenome data

## Extract bacterial reads for each metagenome sample

### Split Kraken outputs

```{bash, eval = F}
find "$WORKING_DIR"/kraken/metagenome -name "*kraken_output.txt" > "$WORKING_DIR"/kraken_output_split.list

echo -e 'KRAKEN_OUTPUT="$(sed -ne "$SGE_TASK_ID p" "$WORKING_DIR"/kraken_output_split.list)" \n
	cat $KRAKEN_OUTPUT | split -l 20000000 -d - "$TEMP_DIR"/metagenome/"$(basename "$KRAKEN_OUTPUT")"' | qsub -v WORKING_DIR="$WORKING_DIR" -wd "$WORKING_DIR"/stderrout -N split -t 1-"$(wc -l "$WORKING_DIR"/kraken_output_split.list | awk '{print $1}')"
```

### Split FASTQ files

```{bash, eval = F}
find "$WORKING_DIR"/final_reads/metagenome/ -name "*fastq.gz" | sed "s/[.].*//g" | sort -n | uniq > "$WORKING_DIR"/metagenome_split.list

echo -e 'SAMPLE="$(sed -ne "$SGE_TASK_ID p" "$WORKING_DIR"/metagenome_split.list)" \n
	zcat "$SAMPLE".fastq.gz | split -l 120000000 -d - "$TEMP_DIR"/metagenome/"$(basename $SAMPLE)".fastq' | qsub -v WORKING_DIR="$WORKING_DIR" -wd "$WORKING_DIR"/stderrout -N split -t 1-"$(wc -l "$WORKING_DIR"/metagenome_split.list | awk '{print $1}')"
```

### Create subset metagenome FASTQs with only bacteria assigned reads

```{bash, eval = F}
module load KrakenTools/0.1-alpha
module load python/3.7.3-foss-2016b

for KRAKEN_OUTPUT in $(find "$TEMP_DIR"/metagenome/ -name "*kraken_output*")
do
	SAMPLE="$(basename "$KRAKEN_OUTPUT" | sed "s/[.]kraken_output.*//g")"

	find "$TEMP_DIR"/metagenome/ | grep "$SAMPLE" | grep -v include | grep -v exclude > "$TEMP_DIR"/metagenome/"$SAMPLE"_fastq.list

	echo -e 'FASTQ="$(sed -ne "$SGE_TASK_ID p" "$TEMP_DIR"/metagenome/"$SAMPLE"_fastq.list)" \n
		extract_kraken_reads.py --include-children --fastq-output -t 2 \
		--max 1000000000 \
		-s "$FASTQ" \
		-k "$KRAKEN_OUTPUT" \
		-r "$WORKING_DIR"/kraken/metagenome/"$SAMPLE".kraken_report.txt \
		-o "$(echo "$FASTQ" | sed "s/[.]fastq/.include.d_bacteria.fastq/g")"_${KRAKEN_OUTPUT: -2}' | qsub -v SAMPLE="$SAMPLE",KRAKEN_OUTPUT="$KRAKEN_OUTPUT",WORKING_DIR="$WORKING_DIR" -l mem_free=200G -wd "$WORKING_DIR"/stderrout -N extract_kraken_reads_bacteria -t 1-"$(wc -l "$TEMP_DIR"/metagenome/"$SAMPLE"_fastq.list | awk '{print $1}')"
done
```

### Combine bacteria-selected FASTQs for each metagenome sample

```{bash, eval = F}
for SAMPLE in $(find "$TEMP_DIR"/metagenome/ -name "*.include.d_bacteria.fastq[0-9]*" | sed "s/[.].*//g" | sort -n | uniq)
do
	for i in $(ls "$SAMPLE"* | grep ".include.d_bacteria.fastq[0-9]*" | sort -n)
	do
		cat $i >> "$WORKING_DIR"/kraken/metagenome/"$(basename "$SAMPLE")".include.d_bacteria.fastq
	done
done
```

### Compress combined bacteria-selected FASTQs
```{bash, eval = F}
find "$WORKING_DIR"/kraken/ -name "*fastq" > "$WORKING_DIR"/gzip.list

echo -e 'FASTQ="$(sed -ne "$SGE_TASK_ID p" "$WORKING_DIR"/gzip.list)" \n
	gzip $FASTQ' | qsub -v WORKING_DIR="$WORKING_DIR" -N gzip -wd "$WORKING_DIR"/stderrout/ -t 1-"$(wc -l "$WORKING_DIR"/gzip.list | awk '{print $1}')"

```

## Identify spacers and direct repeats from metagenome bacterial reads using Crass

### Run Crass

```{bash, eval = F}
module load crass/1.0.1

find "$WORKING_DIR"/kraken/metagenome -name "*include.d_bacteria.fastq.gz" > "$WORKING_DIR"/crass.list

echo -e 'FASTQ="$(sed -ne "$SGE_TASK_ID p" "$WORKING_DIR"/crass.list)" \n
	crass -o "$WORKING_DIR"/crass/"$(basename "$FASTQ" | sed "s/[.].*//g")" -l 4 "$FASTQ"' |  qsub -v WORKING_DIR="$WORKING_DIR"  -wd "$WORKING_DIR"/stderrout -N crass -t 1-"$(wc -l "$WORKING_DIR"/crass.list | awk '{print $1}')"
```

### Combine Crass FASTA files for each sample

The follow samples had no linkers or spacers assembled from Crass:

UCS_0002_BAL,UCS_0017_BAL,UCS_0018_BAL,UCS_0020_BAL,UCS_0021_BAL,
UCS_0022_BAL,UCS_0029_BAL,UCS_0031_BAL,UCS_0033_BAL,UCS_0042_BAL,
UCS_0051_BAL,UCS_0056_BAL,UCS_0059_BAL,UCS_0064_BAL,UCS_0073_BAL,
UCS_0082_BAL,UCS_0097_BAL,UCS_0100_BAL,UCS_0147_BAL,UCS_0153_BAL,
UCS_0159_BAL,UCS_0192_BAL

UCS_0040_BKG,UCS_0069_BKG,UCS_0070_BKG,UCS_0085_BKG

UCS_0083_UA,UCS_0192_UA

```{bash, eval = F}
for CRASS_DIR in $(find "$WORKING_DIR"/crass/ -maxdepth 1 -type d)
do
	cat "$CRASS_DIR"/*.fa > "$CRASS_DIR"/"$(basename "$CRASS_DIR")".crass.fna
done
```

## Conduct BLAST search on Crass spacers

### Bacterial DB

```{bash, eval = F}
module load blast+/2.9.0

find "$WORKING_DIR"/crass -name "*crass.fna" > "$WORKING_DIR"/crass_blastn.list

echo -e 'FASTA="$(sed -ne "$SGE_TASK_ID p" "$WORKING_DIR"/crass_blastn.list)" \n 
blastn \
-db nt \
-query "$FASTA" \
-out "$(echo "$FASTA" | sed "s/[.].*//g")".crass.bacteria.blastn.out \
-remote -entrez_query "Bacteria[Organism]" \
-max_target_seqs 1 -max_hsps 1 \
-evalue 1e-5 -outfmt 6' | qsub  -wd "$WORKING_DIR"/stderrout -v BLAST_DB_DIR="$BLAST_DB_DIR",WORKING_DIR="$WORKING_DIR" -N blast_crass_bacteria -t 1-"$(wc -l "$WORKING_DIR"/crass_blastn.list | awk '{print $1}')"
```

### Viral DB

```{bash, eval = F}
module load blast+/2.9.0

find "$WORKING_DIR"/crass -name "*crass.fna" > "$WORKING_DIR"/crass_blastn.list

echo -e 'FASTA="$(sed -ne "$SGE_TASK_ID p" "$WORKING_DIR"/crass_blastn.list)" \n 
blastn \
-db nt \
-query "$FASTA" \
-out "$(echo "$FASTA" | sed "s/[.].*//g")".crass.virus.blastn.out \
-remote -entrez_query "Viruses[Organism]" \
-max_target_seqs 1 -max_hsps 1 \
-evalue 1e-5 -outfmt 6' | qsub  -wd "$WORKING_DIR"/stderrout -v BLAST_DB_DIR="$BLAST_DB_DIR",WORKING_DIR="$WORKING_DIR" -N blast_crass_virus -t 1-"$(wc -l "$WORKING_DIR"/crass_blastn.list | awk '{print $1}')"
```

## Conduct analysis on bacterial-viral interactions present in patient samples

### Set R inputs

```{r}
WORKING.DIR <- ""
```

### Load packages and view sessionInfo

```{r}
library(bipartite)
library(ggplot2)
library(pvclust)
library(RColorBrewer)
library(reshape2)
library(taxonomizr)

sessionInfo()
```

```{r, eval = F}
R version 4.0.2 (2020-06-22)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 18363)

Matrix products: default

locale:
[1] LC_COLLATE=English_United States.1252  LC_CTYPE=English_United States.1252    LC_MONETARY=English_United States.1252
[4] LC_NUMERIC=C                           LC_TIME=English_United States.1252    

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] taxonomizr_0.5.3     reshape2_1.4.4       RColorBrewer_1.1-2   pvclust_2.2-0        ggplot2_3.3.3        bipartite_2.15      
 [7] sna_2.6              network_1.16.1       statnet.common_4.4.1 vegan_2.5-7          lattice_0.20-41      permute_0.9-5       

loaded via a namespace (and not attached):
 [1] Rcpp_1.0.5        plyr_1.8.6        pillar_1.4.7      compiler_4.0.2    tools_4.0.2       dotCall64_1.0-0   lifecycle_0.2.0  
 [8] tibble_3.0.4      rle_0.9.2         nlme_3.1-148      gtable_0.3.0      mgcv_1.8-31       pkgconfig_2.0.3   rlang_0.4.7      
[15] Matrix_1.2-18     igraph_1.2.6      rstudioapi_0.13   parallel_4.0.2    spam_2.6-0        xfun_0.17         coda_0.19-4      
[22] stringr_1.4.0     withr_2.3.0       dplyr_1.0.2       cluster_2.1.0     knitr_1.30        generics_0.1.0    vctrs_0.3.6      
[29] fields_11.6       maps_3.3.0        tidyselect_1.1.0  grid_4.0.2        data.table_1.13.6 glue_1.4.2        R6_2.5.0         
[36] purrr_0.3.4       magrittr_2.0.1    scales_1.1.1      ellipsis_0.3.1    MASS_7.3-51.6     splines_4.0.2     colorspace_2.0-0 
[43] stringi_1.5.3     munsell_0.5.0     crayon_1.3.4     
```

### Parse BLASTN output files for Crass hits into a data frame

```{r}
evalue.cutoff <- 1e-6

sample.list <- list.files(paste0(WORKING.DIR,"/crass/"),
                          recursive = T,
                          pattern = "*.crass.fna")
sample.list <- gsub("\\/.*","",sample.list)

irb_exclude_samples.list <- c("UCS_0045","UCS_0052","UCS_0170","UCS_0171")

exclude_samples1.list <- c("UCS_0002_BAL","UCS_0017_BAL","UCS_0018_BAL","UCS_0020_BAL","UCS_0021_BAL",
                          "UCS_0022_BAL","UCS_0029_BAL","UCS_0031_BAL","UCS_0033_BAL","UCS_0042_BAL",
                          "UCS_0051_BAL","UCS_0056_BAL","UCS_0059_BAL","UCS_0064_BAL","UCS_0073_BAL",
                          "UCS_0082_BAL","UCS_0097_BAL","UCS_0100_BAL","UCS_0147_BAL","UCS_0153_BAL",
                          "UCS_0159_BAL","UCS_0192_BAL","UCS_0040_BKG","UCS_0069_BKG","UCS_0070_BKG",
                          "UCS_0085_BKG","UCS_0083_UA","UCS_0192_UA")

exclude_samples2.list <- c("UCS_0005_BAL","UCS_0015_BAL","UCS_0040_BAL","UCS_0076_BAL","UCS_0103_BAL",
                           "UCS_0144_BAL","UCS_0191_BAL")

sample.list <- sample.list[!(gsub("_[BU].*","",sample.list) %in% irb_exclude_samples.list)]
sample.list <- sample.list[!(sample.list %in% exclude_samples1.list)]
sample.list <- sample.list[!(sample.list %in% exclude_samples2.list)]

crass.df <- as.data.frame(matrix(nrow=0,
                                 ncol=4))
crass.df.list <- list()
for(i in 1:length(sample.list)){
  bacteria_blastn.out <- read.delim(paste0(WORKING.DIR,"/crass/",sample.list[i],"/",sample.list[i],".crass.bacteria.blastn.out"),
                                    header = F)
  virus_blastn.out <- read.delim(paste0(WORKING.DIR,"/crass/",sample.list[i],"/",sample.list[i],".crass.virus.blastn.out"),
                                header = F)
  bacteria_blastn.out <- bacteria_blastn.out[bacteria_blastn.out[,11] < evalue.cutoff,]
  virus_blastn.out <- virus_blastn.out[virus_blastn.out[,11] < evalue.cutoff,]
  
  crass.reads <- unique(c(bacteria_blastn.out[,1],virus_blastn.out[,1]))
 
  crass.df.list[[i]] <- as.data.frame(matrix(nrow=length(crass.reads),
                                             ncol=3))
  crass.df.list[[i]][,1] <- crass.reads
  crass.df.list[[i]][,2] <- bacteria_blastn.out[match(crass.reads,bacteria_blastn.out[,1]),2]
  crass.df.list[[i]][,3] <- virus_blastn.out[match(crass.reads,virus_blastn.out[,1]),2]
 
  crass.df <- as.data.frame(rbind(crass.df,
                                  as.data.frame(cbind(sample.list[i],crass.df.list[[i]]))))
}

crass.df <- crass.df[rowSums(is.na(crass.df)) == 0,] 
```

### Map taxonomic identifers to accession ids

```{r}
taxid.vector <- unique(c(crass.df[,3],crass.df[,4]))
taxid.vector <- accessionToTaxa(taxid.vector,"nameNode.sqlite")
tax.df <- as.data.frame(getTaxonomy(taxid.vector,'nameNode.sqlite'))

tax.df$accession <- unique(c(crass.df[,3],crass.df[,4]))

crass.tax.df <- crass.df
crass.tax.df[,3] <- tax.df$species[match(crass.df[,3],tax.df$accession)]
crass.tax.df[,4] <- tax.df$species[match(crass.df[,4],tax.df$accession)]
crass.tax.df <- crass.tax.df[rowSums(is.na(crass.tax.df)) == 0,]
```

### Read in meta-data maps

```{r}
metagenome.map <- read.delim(paste0(WORKING.DIR,"/maps/201015.WGS.txt"),
                                    row.names = 1)
metagenome.map <- metagenome.map[match(unique(crass.df[,1]),gsub("[.]","_",rownames(metagenome.map))),]

metagenome.map$deseq2_group[metagenome.map$Composite.Month.LNS == "Good"] <- "<28 days"
metagenome.map$deseq2_group[metagenome.map$Composite.Month.LNS == "Bad" & metagenome.map$death_alive_kk != "death"] <- ">28 days"
metagenome.map$deseq2_group[metagenome.map$Composite.Month.LNS == "Bad" & metagenome.map$death_alive_kk == "death"] <- "death"
```

### Read in list of clinically relevant meta-data and remove continuous meta-data variables

```{r}
relevant_metadata.map <- read.delim(paste0(WORKING.DIR,"/maps/vars2.Univariate.Omic.201104.txt"))
relevant_metadata.map <- relevant_metadata.map[relevant_metadata.map$Categorical_Num != 2,]

relevant_metadata.map <- relevant_metadata.map[,]
metagenome.clinical.map <- metagenome.map[,colnames(metagenome.map) %in% relevant_metadata.map$Vars[relevant_metadata.map$Clinically_Relevant == 1],drop = F]
metagenome.topo.map <- metagenome.map[,colnames(metagenome.map) %in% relevant_metadata.map$Vars[relevant_metadata.map$Topographical_relevant == 1],drop = F]
metagenome.deseq2.map <- metagenome.map[,colnames(metagenome.map) == "deseq2_group",drop = F]
```

### Remove bad viral hits

```{r}
viral_bad_hits.vector <- c("Homo sapiens", "Staphylococcus aureus", "Staphylococcus epidermidis",
                           "Streptococcus pneumoniae","Streptococcus pyogenes","Saccharopolyspora erythraea",
                           "Pseudomonas aeruginosa","Lactobacillus johnsonii","Erwinia amylovora",
                           "Clostridioides difficile",
                           "Influenza A virus","Parvovirus NIH-CQV",
                           "Rabbit picobirnavirus","Pestivirus D","Pestivirus A",
                           "Norwalk virus")

crass.tax.df <- crass.tax.df[!(crass.tax.df[,4] %in% viral_bad_hits.vector),]
```

### Check for significantly over-represented bacterial-viral interactions in data set

#### Set function for conducting Fisher's exact test on Crass data

```{r}
crass_fishers_exact_test <- function(fishers.df,map){
  fishers.final.output.df <- as.data.frame(matrix(nrow = 0,
                                                  ncol = 7))
  for(i in 1:ncol(map)){
    for(j in 1:length(unique(map[,i]))){
      fishers.subset.df <- fishers.df[fishers.df[,1] %in% gsub("[.]","_",rownames(map)[map[,i] == unique(map[,i])[j]]),]
      terms <- unique(fishers.subset.df[,4])
      fishers.output.df <- as.data.frame(matrix(nrow = length(terms),
                                                ncol = 7))
      for(k in 1:length(terms)){
         a <- nrow(fishers.subset.df[fishers.subset.df[,4] == terms[k],])
         b <- nrow(fishers.df[fishers.df[,4] == terms[k],]) - a
         c <- nrow(fishers.subset.df) - a
         d <- nrow(fishers.df) - b - c
         
         fisherexact.matrix <- matrix(c(a,b,
                                        c,d),
                         					    nrow = 2,
                        					    ncol = 2)
         fishers.test <- fisher.test(fisherexact.matrix)
        
         fishers.output.df[k,] <- c(terms[k],
        	                          a,
        	                          a+b,
        	                          fishers.test$estimate,
                                    colnames(map)[i],
                                    unique(map[,i])[j],
                                    fishers.test$p.value)
      }
      fishers.final.output.df <- as.data.frame(rbind(fishers.final.output.df,
                                                     fishers.output.df))
    }
  }
  fishers.final.output.df[,8] <- p.adjust(fishers.final.output.df[,7], method="bonferroni",n=nrow(fishers.final.output.df))
  colnames(fishers.final.output.df) <- c("term",
                                         "cluster_occurences",
                                         "genome_occurences",
                                         "oddsratio",
                                         "category",
                                         "group",
                                         "p-value",
                                         "FDR")
  for(i in c(2,3,4,7,8)){
    fishers.final.output.df[,i] <- as.numeric(as.character(fishers.final.output.df[,i]))
  }
  return(fishers.final.output.df)
}
```

#### Checking for significantly over-represented bacteria-virus interactions

```{r}
fishers.df <- unique(crass.tax.df[,c(1,3,4)])
fishers.df[,4] <- paste(fishers.df[,2],fishers.df[,3],sep="|")
fishers.bal.df <- fishers.df[grep("BAL",fishers.df[,1]),]

fishers.topo.output.df <- crass_fishers_exact_test(fishers.df,metagenome.topo.map)
fishers.topo.output.sig.df <- fishers.topo.output.df[fishers.topo.output.df$oddsratio > 1 & fishers.topo.output.df$FDR < 0.05,]

fishers.clinical.output.df <- crass_fishers_exact_test(fishers.bal.df,metagenome.clinical.map)
fishers.clinical.output.sig.df <- fishers.clinical.output.df[fishers.clinical.output.df$oddsratio > 1 & fishers.clinical.output.df$FDR < 0.05,]

fishers.deseq2.output.df <- crass_fishers_exact_test(fishers.bal.df,metagenome.deseq2.map)
fishers.deseq2.output.sig.df <- fishers.deseq2.output.df[fishers.deseq2.output.df$oddsratio > 1 & fishers.deseq2.output.df$FDR < 0.05,]
```

#### Checking for significantly over-represented bacteria present in bacteria-virus interactions

```{r}
fishers.bacteria.df <- fishers.df
fishers.bacteria.df[,4] <- fishers.bacteria.df[,2]
fishers.bacteria.df <- fishers.bacteria.df[!(duplicated(paste(fishers.bacteria.df[,1],fishers.bacteria.df[,4]))),]

fishers.bal.df <- fishers.bacteria.df[grep("BAL",fishers.bacteria.df[,1]),]

fishers.topo.bacteria.output.df <- crass_fishers_exact_test(fishers.bacteria.df,metagenome.topo.map)
fishers.topo.bacteria.output.sig.df <- fishers.topo.bacteria.output.df[fishers.topo.bacteria.output.df$oddsratio > 1 & fishers.topo.bacteria.output.df$FDR < 0.05,]

fishers.clinical.bacteria.output.df <- crass_fishers_exact_test(fishers.bal.df,metagenome.clinical.map)
fishers.clinical.bacteria.output.sig.df <- fishers.clinical.bacteria.output.df[fishers.clinical.bacteria.output.df$oddsratio > 1 & fishers.clinical.bacteria.output.df$FDR < 0.05,]

fishers.deseq2.bacteria.output.df <- crass_fishers_exact_test(fishers.bal.df,metagenome.deseq2.map)
fishers.deseq2.bacteria.output.sig.df <- fishers.deseq2.bacteria.output.df[fishers.deseq2.bacteria.output.df$oddsratio > 1 & fishers.deseq2.bacteria.output.df$FDR < 0.05,]
```

#### Checking for significantly over-represented viruses present in bacteria-virus interactions

```{r}
fishers.virus.df <- fishers.df
fishers.virus.df[,4] <- fishers.virus.df[,3]
fishers.virus.df <- fishers.virus.df[!(duplicated(paste(fishers.virus.df[,1],fishers.virus.df[,4]))),]

fishers.bal.df <- fishers.virus.df[grep("BAL",fishers.virus.df[,1]),]

fishers.topo.virus.output.df <- crass_fishers_exact_test(fishers.virus.df,metagenome.topo.map)
fishers.topo.virus.output.sig.df <- fishers.topo.virus.output.df[fishers.topo.virus.output.df$oddsratio > 1 & fishers.topo.virus.output.df$FDR < 0.05,]

fishers.clinical.virus.output.df <- crass_fishers_exact_test(fishers.bal.df,metagenome.clinical.map)
fishers.clinical.virus.output.sig.df <- fishers.clinical.virus.output.df[fishers.clinical.virus.output.df$oddsratio > 1 & fishers.clinical.virus.output.df$FDR < 0.05,]

fishers.deseq2.virus.output.df <- crass_fishers_exact_test(fishers.bal.df,metagenome.deseq2.map)
fishers.deseq2.virus.output.sig.df <- fishers.deseq2.virus.output.df[fishers.deseq2.virus.output.df$oddsratio > 1 & fishers.deseq2.virus.output.df$FDR < 0.05,]
```

### Creating plotting data frames for bacterial and viral interactions

```{r}
crass.genera_v_species.df <- crass.df
crass.genera_v_species.df[,3] <- tax.df$genus[match(crass.genera_v_species.df[,3],tax.df$accession)]
crass.genera_v_species.df[,4] <- tax.df$species[match(crass.genera_v_species.df[,4],tax.df$accession)]
crass.genera_v_species.df <- crass.genera_v_species.df[rowSums(is.na(crass.genera_v_species.df)) == 0,] 
crass.genera_v_species.df <- unique(crass.genera_v_species.df[,-2])

crass.species_v_species.df <- crass.df
crass.species_v_species.df[,3] <- tax.df$species[match(crass.species_v_species.df[,3],tax.df$accession)]
crass.species_v_species.df[,4] <- tax.df$species[match(crass.species_v_species.df[,4],tax.df$accession)]
crass.species_v_species.df <- crass.species_v_species.df[rowSums(is.na(crass.species_v_species.df)) == 0,] 
crass.species_v_species.df <- unique(crass.species_v_species.df[,-2])

crass.genera_v_species.bal.df <- crass.genera_v_species.df[grep("BAL",crass.genera_v_species.df[,1]),]
crass.genera_v_species.ua.df <- crass.genera_v_species.df[grep("UA",crass.genera_v_species.df[,1]),]

crass.species_v_species.bal.df <- crass.species_v_species.df[grep("BAL",crass.species_v_species.df[,1]),]
crass.species_v_species.ua.df <- crass.species_v_species.df[grep("UA",crass.species_v_species.df[,1]),]
```

### Plot heatmaps summarizing the relative abundance and presence of CRISPR array data

#### Set function for plotting a binary heatmap of the presence of Crass interactions

```{r,fig.height=6,fig.width=11}
plot_crass_binaryheatmap <- function(crass.df,aes.map){
  set.seed(2)
  
  sample_cutoff <- 5

  counts.df <- as.data.frame(table(paste0(crass.df[,2],"|",crass.df[,3])))
  counts.df <- counts.df[counts.df[,2] >= sample_cutoff,] 
  
  plot.df <- as.data.frame(cbind(gsub("[|].*","",counts.df[,1]),
                                 gsub(".*[|]","",counts.df[,1]),
                                 1))
  plot.df <- plot.df[!(plot.df[,2] %in% viral_bad_hits.vector),]
  
  plot.df[,3] <- as.numeric(as.character(plot.df[,3]))
  
  colorder.df <- dcast(plot.df, V1 ~ V2)
  colorder.df[is.na(colorder.df)] <- 0
  rownames(colorder.df) <- colorder.df[,1]
  colorder.df <- colorder.df[,-1]
  
  bipartite.modules <- metaComputeModules(colorder.df)
  bipartite.modules.info <- listModuleInformation(bipartite.modules)
  
  taxa.vector <- c()
  modules.vector <- c()
  for(i in 1:length(bipartite.modules.info[[2]])){
    taxa.vector <- c(taxa.vector,bipartite.modules.info[[2]][[i]][[1]],bipartite.modules.info[[2]][[i]][[2]])
    modules.vector <- c(modules.vector,rep(i,length(c(bipartite.modules.info[[2]][[i]][[1]],bipartite.modules.info[[2]][[i]][[2]]))))
  }
  rowlevels <- c(rev(rownames(bipartite.modules@moduleWeb)), unique(plot.df[!(plot.df[,1] %in% rev(rownames(bipartite.modules@moduleWeb))),1]))
  collevels <- c(unique(plot.df[!(plot.df[,2] %in% colnames(bipartite.modules@moduleWeb)),2]),colnames(bipartite.modules@moduleWeb))
  
  plot.df[,1] <- factor(plot.df[,1],levels=rowlevels)
  plot.df[,2] <- factor(plot.df[,2],levels=collevels)
  plot.df <- plot.df[!(is.na(plot.df[,1])) & !(is.na(plot.df[,2])),]
  
  plot.df$xlab_col <- aes.map$xlab_col[match(plot.df[,1],aes.map[,2])]
  plot.df$ylab_col <- aes.map$ylab_col[match(plot.df[,2],aes.map[,3])]
  plot.df$cell_col <- aes.map$cell_col[match(paste0(plot.df[,1],"|",plot.df[,2]),paste0(aes.map[,2],"|",aes.map[,3]))]
  plot.df$cell_opa <- aes.map$cell_opa[match(paste0(plot.df[,1],"|",plot.df[,2]),paste0(aes.map[,2],"|",aes.map[,3]))]

  xlab_col <- plot.df$xlab_col[match(sort(unique(plot.df[,1])),plot.df[,1])]
  ylab_col <- plot.df$ylab_col[match(sort(unique(plot.df[,2])),plot.df[,2])]
  
  for(i in 1:nrow(plot.df)){
    if(length(modules.vector[taxa.vector == plot.df[i,1]]) ==  0){
      plot.df[i,4] <- max(modules.vector) + 1
    }else{
      if(modules.vector[taxa.vector == plot.df[i,1]] == modules.vector[taxa.vector == plot.df[i,2]]){
        plot.df[i,4] <- modules.vector[taxa.vector == plot.df[i,1]]
      }else{
        plot.df[i,4] <- 0
      }
    }
  }
  
  fillcol <- brewer.pal(length(unique(plot.df[,4])), "Set1")
  fillcol[1] <- "grey"
  
  heatmap <- ggplot(mapping=aes(x=!!plot.df[,1],
                                y=!!plot.df[,2],
                                fill=!!as.character(plot.df[,4])))+
    geom_tile(alpha=plot.df$cell_opa, color = plot.df$cell_col)+
    guides(fill=F)+
    labs(x="",y="")+
    theme_minimal()+
    scale_fill_manual(values = fillcol)+
    theme(axis.text.x = element_text(color=xlab_col,angle = 90,vjust = 0.5, hjust = 1,size=8),
          axis.text.y = element_text(color=ylab_col,size=8))
  
  return(heatmap)
}
```

#### Set an aesthetic map to color significant bacterial and viral members and interactions between BAL and UA samples

```{r}
aes.map <- crass.species_v_species.df
aes.map$xlab_col <- "black"
aes.map$xlab_col[which(aes.map[,2] %in% fishers.topo.bacteria.output.sig.df$term[fishers.topo.bacteria.output.sig.df$group == "BAL"])] <- "blue"
aes.map$xlab_col[which(aes.map[,2] %in% fishers.topo.bacteria.output.sig.df$term[fishers.topo.bacteria.output.sig.df$group == "UA"])] <- "orange"

aes.map$ylab_col <- "black"
aes.map$ylab_col[which(aes.map[,3] %in% fishers.topo.virus.output.sig.df$term[fishers.topo.virus.output.sig.df$group == "BAL"])] <- "blue"
aes.map$ylab_col[which(aes.map[,3] %in% fishers.topo.virus.output.sig.df$term[fishers.topo.virus.output.sig.df$group == "UA"])] <- "orange"

aes.map$cell_col <- ifelse(paste0(crass.species_v_species.df[,2],"|",crass.species_v_species.df[,3]) %in% fishers.topo.output.sig.df$term,
                                                     "black",NA)
aes.map$cell_opa <- ifelse(paste0(crass.species_v_species.df[,2],"|",crass.species_v_species.df[,3]) %in% fishers.topo.output.sig.df$term,
                                                     1,0.25)
```

#### Plot binary heatmap of Crass interactions

```{r,fig.height=6,fig.width=12}
bal.hm <- plot_crass_binaryheatmap(crass.species_v_species.bal.df,aes.map)
ua.hm <- plot_crass_binaryheatmap(crass.species_v_species.ua.df,aes.map)

pdf(paste0(WORKING.DIR,"/plots/crass_binaryhm.pdf"),
    width=12, 
    height=6)
plot_grid(bal.hm,
          ua.hm,
          ncol=2,
          align='h',
          labels=c("A","B"),
          rel_widths = c(1.1,2))
dev.off()
```
