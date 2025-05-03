# SquiDPipe: taxonomic classification and POD5 extraction pipeline for SquiDBase

<div align="center">
  <img src="squidpipe.png" alt="Description of Image" width="300">
</div>

## About
Nanopore sequencing generates POD5 files containing raw signal data. When dealing with multiplexed runs—where multiple species are sequenced on the same flowcell, you will end up with a set of POD5 and FASTQ files per barcode. In cases where a single barcode contains multiple species (e.g., human and Zika virus reads), SquiDPipe allows you to isolate and create a new POD5 file containing only the raw nanopore data for your species of interest.

This Nextflow pipeline automates this process by:

- Running Kraken2 for taxonomic classification on FASTQ files
- Identifying reads that match a specified taxonomic identifier
- Performing mapping and coverage analysis to validate the extracted data
- Extracting the corresponding raw signal data from POD5 files

Each step performed by the pipeline is explained at the end of this README under [Key pipeline processes](#key-pipeline-processes). 

SquiDpipe is specifically designed for microbial and viral sequencing applications, and was specifically developed to facilitate uploading your raw Nanopore data to [SquiDBase](https://squidbase.org/). 

## Table of contents
- [Installation](#installation)
- [Usage](#usage)
- [Output](#output)
- [CSV_file_specifications](#about-the-input-csv-file)
- [Key pipeline processes](#key-pipeline-processes)
- [Acknowledgements](#key-pipeline-processes)


## Installation

1. Install [Nextflow](https://www.nextflow.io/docs/latest/index.html) to run the pipeline.

2. Clone this repository:
```bash
git clone https://github.com/Cuypers-Wim/squidpipe.git
cd squidpipe
```

3. Configure the pipeline by editing the `nextflow.config` file to select from which environment you whish to run the pipeline.

### Evironment settings

- Run using Conda: set `conda.enabled = true`
- Run using Docker: set `docker.enabled = true`
- Run using Apptainer: set `apptainer.enabled = true`

Note that if you do not whish to use the above options, you can install the following dependencies:

- [Conda](https://docs.conda.io/en/latest/)
- [h5py](https://www.h5py.org/)
- [Kraken2](https://ccb.jhu.edu/software/kraken2/)
- [minimap2](https://github.com/lh3/minimap2)
- [ncbi-datasets-cli](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/command-line-tools/download-and-install/)
- [pandas](https://pandas.pydata.org/)
- [pigz](https://zlib.net/pigz/)
- [pod5](https://github.com/nanoporetech/pod5-file-format)
- [samtools](http://www.htslib.org/)
- [seqkit](https://bioinf.shenwei.me/seqkit/)
- [seqtk](https://github.com/lh3/seqtk)
- [vbz-h5py-plugin](https://github.com/nanoporetech/vbz-h5py-plugin)

## Usage 

Before running **SquiDPipe**, you need to set the required parameters. The easiest way to do this is by modifying the `nextflow.config` file. Below are the key parameters that should be adjusted based on your dataset and computing environment.  

### Required Parameters  

#### Input Files and Directories  
These parameters specify the locations of input data:  

- **`csv_file`** – Path to the CSV file containing metadata and sample information.  
  - Example: `'full.csv'`  
- **`fastq_dir`** – Directory containing FASTQ files used for taxonomic classification.  
  - Example: `'data/fastq'`  
- **`pod5_dir`** – Directory containing POD5 files with raw nanopore signals.  
  - Example: `'data/pod5'`  

#### Processing Options  
These options control the pipeline behavior:  

- **`csvMeta`** – If `true`, includes all metadata and generates an input file for Squidbase.  
  - Options: `true` or `false`  
  - Read more about the structure of the CSV file in the section - [CSV file](#about-the-input-csv-file)
- **`pod5_split`** – If `true`, generates subsetted POD5 files based on taxonomic classification.  
  - Options: `true` or `false`  

#### Database Configuration  
The pipeline requires a Kraken2 database for taxonomic classification:  
- **`kraken_db`** – Path to the Kraken2 database.  
  - Example: `'/home/the_darwim/databases/kraken2/refseq_virus'`  

In addition, you need to specify which reference genomes need to be included in addition to already specified taxonomic identifiers. For example, if you're mapping data originating from a human patient containing both virus and human reads, you want to include the human reference genome to improve the quality of SquiDPipe's mapping step:
- **`references`** – Path to the reference genome file (if needed for downstream analysis).  
  - Example: `'/some_location/GCF_000001405.40_GRCh38.p14_genomic.fna'`  

### Running the pipeline

```bash
nextflow run main.nf 
```

## Output

The pipeline outputs mapping statistics for the reference genome or the best matching reference genome. Additionally, if requested, the pipeline generates POD5 files containing the raw nanopore signal data (the so-called "squiggles") specific to the identified species. These POD5 files can be further uploaded to SquiDBase.

Running this pipeline generates a `results` folder containing multiple the output files, and  subdirectories containing important intermediate output.

### Output files

- **`coverage_metrics.csv`**: Contains coverage statistics for each sample contig (entire genome or chromosome). Columns include:  
  1. Contig name  
  2. Read coverage (percentage of positions covered by at least one read)  
  3. Depth of coverage (average number of reads covering a position)  
  4. Percentage of positions with ≥10 reads  
  5. Percentage of positions with ≥20 reads  
  6. Percentage of positions with ≥30 reads  

- **`final_metadata.csv`** *(optional)*: Generated if `csvMeta` is set to `true` and metadata is correctly provided. This file is intended for upload to SquiDBase and links POD5 files (column: `filename`) to metadata fields such as:  
  - `species_name`  
  - `species_taxid`  
  - `year_of_isolation`  
  - `country_of_isolation`  
  - `geographic_origin`  
  - `strain_lineage`  
  - `source_id`  
  - `host_taxid`  
  - `internal_lab_id`  
  - `diagnostic_method_id`  
  - `remarks`  

### Output directories

- `mapping`: Contains mapping statistics and BAM files related to the read alignment and filtering processes.
- `reads`: Contains the extracted reads in FASTQ format that correspond to the target species.
- `references`: Stores the reference genomes used for mapping.
- `POD5`: Includes POD5 files containing raw nanopore signal data ("squiggles") specific to the identified species.

## About the input CSV file

This pipeline is designed to process data generated from a nanopore sequencing run. Typically, such a run produces folders containing FASTQ and POD5 files, where each barcode folder may contain multiple read files. The pipeline is particularly useful when you need to identify and analyze a specific species within a sample, such as the Chikungunya virus with NCBI taxon identifier 37124, expected to be found in barcode one.

To run the pipeline, you need to provide a CSV file that specifies the expected species and their corresponding NCBI taxon identifiers for each barcode. The pipeline uses this information to extract reads corresponding to the target species from each barcode folder. It then performs mapping to confirm the presence of the species and eliminate any false positives.

###  Example input CSV file

#### Short version

The CSV file required for SquiDPipe must include, at a minimum, the barcode folder (e.g., barcode01), the corresponding taxon name (e.g., CHIKV-1), and its taxonomic identifier (e.g., 37124). An example CSV file is shown below:
|           |           |          |
|-----------|-----------|----------|
| barcode01 | CHIKV-1   | 37124    |
| barcode02 | Denv1-5   | 11053    |
| barcode03 | Denv2-7   | 11060    |
| barcode04 | Denv3-10  | 11069    |
| barcode05 | EEEV-16   | 11021    |


#### Long version including metadata

The full CSV format includes additional metadata. To use this option, the `csvMeta` option must be set to `true` in the NextFlow configuration file. An example CSV file is shown below:

| filename  | species_name    | species_taxid | year_of_isolation | country_of_isolation | geographic_origin | strain_lineage | source_id | host_taxid | internal_lab_id | diagnostic_method_id | remarks         |
|-----------|---------------|--------------|-------------------|---------------------|------------------|---------------|-----------|-----------|--------------|-------------------|----------------|
| barcode01 | DENV          | 11053        | 2014              | NA                  | NA               | ECSA          | NA        | NA        | NA           | NA                | Collection from Institute A |
| barcode02 | HIV           | 11676        | 2015              | BE                  | NI               | ECSA          | NA        | NA        | NA           | NA                | Collection from Institute A |
| barcode03 | ZIKV          | 64320        | 2015              | BE                  | ID               | NA            | NA        | NA        | NA           | NA                | Collection from Institute A |
| barcode04 | SARS-CoV-2_A  | 2697049      | 2018              | BE                  | PE               | NA            | NA        | NA        | NA           | NA                | Collection from Institute A |
| barcode05 | SARS-CoV-2_B  | 2697049      | 2018              | BE                  | PE               | NA            | NA        | NA        | NA           | NA                | Collection from Institute A |


## Key pipeline processes:

- `CONCATENATE_FASTQ`: Merges all FASTQ files in a barcode folder into a single FASTQ file for downstream analysis.
- `RUNKRAKEN2`: Executes Kraken2 on the concatenated FASTQ files for taxonomic classification, generating reports and extracting read IDs.
- `SUBSET_FASTQ`: Extracts reads of interest from the concatenated FASTQ files based on the Kraken2 classification results.
- `RETRIEVE_GENOMEREFS`: Downloads and formats reference genomes based on specified Taxon IDs.
- `MAP_READS`: Aligns the subsetted FASTQ reads to the reference genomes using minimap2 (to remove false-positives), followed by sorting the output BAM file.
- `BAM_PROCESSING`: Marks duplicates in the BAM file and creates an index for efficient access.
- `MAPPING_STATS_ALL`: Generates mapping statistics for the BAM files, including flagstat and idxstats reports.
- `EXTRACT_READIDS`: Filters the BAM file to retain reads that map to the target species and extracts the corresponding read IDs.
- `SAMTOOLS_DEPTH`: Computes the depth of coverage across the reference genomes using the filtered BAM files.
- `DEPTH_OF_COVERAGE`: Analyzes the depth of coverage data to calculate coverage statistics, outputting the results in CSV format.
- `SUBSET_POD5`: Filters POD5 files to retain only the raw signals/squiggles corresponding to the extracted read IDs.


## Citation  

Our publication is currently under review. In the meantime, please cite the SquiDBase [Pre-print](https://doi.org/10.1101/2025.04.28.650941).

SquiDBase: a community resource of raw nanopore data from microbes. 
Wim L. Cuypers, Halil Ceylan, Eline Turcksin, Laura Raes, Nicky de Vrij, Johan Michiels, Sandra Coppens, Tessa de Block, Daan Jansen, Kevin K. Arien, Philippe Selhorst, Koen Vercauteren, Julia M. Gauglitz, Wout Bittremieux, Kris Laukens. 
bioRxiv 2025.04.28.650941; doi: https://doi.org/10.1101/2025.04.28.650941

