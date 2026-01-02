# mLOY Caller

A high-performance CLI tool for detecting **Mosaic Loss of Y (mLOY)** from diverse genomic data formats. [cite_start]This tool implements the LRR (Log R Ratio) and Cell Fraction metrics defined in the [MoChA](https://github.com/freeseek/mocha) framework, with extended support for Illumina DNA methylation arrays and Whole Genome Sequencing.

## Features

* **Multi-Format Support**: Native processing for:
    * **Methylation Arrays**: Illumina IDAT (EPIC, 450k)
    * **Genotyping Arrays**: Illumina IDAT (GSA, Omni, etc.)
    * **Whole Genome Sequencing**: WGS BAM/CRAM and VCF/BCF
    * **Legacy**: Affymetrix CEL files
* **MoChA Compatible**: Outputs standard metrics (LRR, Cell Fraction, Classification).
* **Resource Efficient**: optimized for laptops using batch "chunking" and multicore parallelization (for macOS and Linux).
* **Build Support**: Defaults to `GRCh38` (hg38), with support for `GRCh37` (hg19).
* **Safe Environment**: Robust pre-flight checks and automated tool installation with interactive sudo-checks.


## Installation

1.  **Clone the repository:**
    ```bash
    git clone https://github.com/brianchengithub/mloy_caller.git
    cd mloy_caller
    ```

2.  **Check/Install R/Bioconductor Dependencies:**
    Run the setup script once to install required Bioconductor packages.
    ```bash
    Rscript install_deps.R
    ```

3.  **Setup the Tool:**
    Make script executable.
    ```bash
    chmod +x mloy_setup.sh mloy genome_audit.sh   
    ```

## Usage

The recommended way to run the tool is using the `mloy_setup.sh` script. 

Note: Windows users must install WSL2 (Ubuntu recommended). Open PowerShell as Administrator and run "wsl --install". After restarting, all commands below should be run inside the WSL terminal.

### Basic Syntax
```bash
./mloy_setup.sh <INPUT_DIRECTORY_OR_FILE> [EXTRA_OPTIONS]
```


### Examples by Data Type

### 1. Methylation Arrays (Illumina IDAT)
#### Recursively detects .idat files. No manifest required (uses Sesame annotation)
```bash
./mloy_setup.sh /path/to/idats/ --cores 6 --chunk 100
```

### 2. Genotyping Arrays (Illumina IDAT)
#### Critical: You must provide a CSV manifest mapping probes to chromosomes.
```bash
./mloy_setup.sh /path/to/geno_idats/ --manifest manifests/GSA_v3.csv
```

#### Manifest Format: Must contain columns Address, Chromosome, and Position.


### 3. WGS Data (BAM/CRAM)
```bash
./mloy_setup.sh /path/to/bams/ --build GRCh38
```

### 4. Variant Calls (VCF/BCF)
#### Extracts depth from DP or AD fields. Ideal for gVCFs.
```bash
./mloy_setup.sh /path/to/vcfs/  
``` 


## Options
| Option | Description | Default |
| :--- | :--- | :--- |
| `--build` | Genome build (`GRCh38` or `GRCh37`) | `GRCh38` |
| `--cores` | Number of CPU cores to use | `4` |
| `--chunk` | Number of files to process per RAM batch | `50` |
| `--manifest` | Path to CSV manifest (Required for Genotyping) | `NULL` |
| `--type` | Manually force file type (`meth`, `geno`, `bam`, `vcf`) | `Auto` |
| `--quiet` | Suppress progress bar | `False` |



## Output Format
### The output is a Tab-Separated Value (TSV) stream printed to stdout. You can redirect it to a file:
```bash
./mloy_setup.sh data/ > results.tsv
```

| Column | Description |
| :--- | :--- |
| `sample_id` | Filename or Sample ID |
| `computed_gender` | Inferred sex (**M**ale/**F**emale) based on Y-intensity |
| `rel_cov` | **Log R Ratio (LRR)**: Log2 ratio of observed Y vs. expected male Y |
| `type` | Classification: `Loss`, `Normal`, or `Undetermined` |
| `cf` | **Cell Fraction**: Estimate of cells with Y loss (0.0 - 1.0) |
| `beg_pos` / `end_pos` | Genomic coordinates of the analyzed nonPAR region |