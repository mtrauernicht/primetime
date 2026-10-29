<p align="center"><img src="misc/logo.png" alt="primetime" width="80%"></p>

# Table of contents

- [**Introduction**](#introduction)
- [**Installation**](#installation)
- [**Setting up your configuration file**](#setting-up-your-configuration-file)
  - [**Setting up the input files**](#setting-up-the-input-files)
  - [**Setting up the conditions for the comparative analysis**](#setting-up-the-conditions-for-the-comparative-analysis)
  - [**Setting up the output directory**](#setting-up-the-output-directory)
  - [**Optional: Changing the p-value threshold**](#optional-changing-the-p-value-threshold)
  - [**Optional: Changing the library information**](#optional-changing-the-library-information)
- [**Running primetime**](#running-primetime)
- [**Output Files**](#output-files)

# Introduction

**primetime** is a pipeline designed for the analysis of TF prime reporter data. It processes fastq files, counts barcodes, clusters them, annotates them, and performs a differential TF activity analysis. The pipeline compares different conditions of the samples based on the 'condition' field in the configuration file.

# Installation

To install and run **primetime**, follow these steps:

1. Clone the repository:
    ```sh
    git clone https://github.com/vansteensellab/primetime.git
    cd primetime
    ```

2. Make sure you have snakemake installed. If you don't have it, you can install it with conda:
    ```sh
    conda install -c bioconda snakemake
    ```

We recommend trying to run primetime with our test data to check if everything was installed correctly.
This should run without any errors.

```sh
snakemake --configfile config.yaml --use-conda --cores 10 --printshellcmds
```

# Setting up your configuration file

Before running the **primetime**, you need to change the parameters on the `config.yaml` according to your data. Below, we will use our test data (that is meant to access the changes on TF activity on U2OS cells uppon calcitriol treatment) to guide you through the different sections of the configuration file and how to set them up:

> Note: The configuration file is sensitive to the number of spaces before each field. Make sure the format of your file is the same as the examples we provide.

## Setting up the input files

The `INPUT_DATA` parameter specifies the input data for the pipeline. Each sample should have a unique identifier and include information about whether it is pDNA, the condition, and the paths to the fastq files.

Example:
```yaml
INPUT_DATA:
  U2OS_DMSO_12W:
    is_pdna: False
    condition: DMSO
    fastq: 
      - test_data/U2OS_DMSO_1.fastq.gz
      - test_data/U2OS_DMSO_2.fastq.gz
      - test_data/U2OS_DMSO_3.fastq.gz

  U2OS_Calcitriol_12W:
    is_pdna: False
    condition: Calcitriol
    fastq: 
      - test_data/U2OS_calcitriol_1.fastq.gz
      - test_data/U2OS_calcitriol_2.fastq.gz
      - test_data/U2OS_calcitriol_3.fastq.gz

  pDNA:
    is_pdna: True
    condition: None
    fastq: 
      - test_data/pDNA_1.fastq.gz
```

## Setting up the conditions for the comparative analysis

The `COMPARATIVE_ANALYSIS` parameter defines the reference and contrast conditions for the comparative analysis. For our test data, we are comparing the cells treated with calcitriol against the control cells, treated with DMSO.

Example:
```yaml
COMPARATIVE_ANALYSIS:
  REFERENCE_CONDITION: DMSO
  CONTRAST_CONDITION: Calcitriol
```

## Setting up the output directory

The `OUTPUT_DIRECTORY` parameter specifies the directory where the output files will be saved.

Example:
```yaml
OUTPUT_DIRECTORY: test_data_output
```

## Optional: Changing the p-value threshold

The `PVALUE_THRESHOLD` parameter sets the p-value threshold for the differential analysis.

Example:
```yaml
PVALUE_THRESHOLD: 0.05
```

## Optional: Turning normalization on or off

Set `NORMALIZE_COUNTS` to `True` to keep the promoter-based activity correction, or to `False` to skip it and run the comparative analysis on the raw cDNA-derived activities.

Example:
```yaml
NORMALIZE_COUNTS: True
```

`RANDOM_NORM_METHOD` selects how the shift is estimated from the RANDOM (no-TF) reporters:
- `median` (default): subtract the median of each promoter's RANDOM reporters.
- `shrunk`: reporters are precision-weighted (BCalm standard error plus reporter-to-reporter excess variance), a pooled shift is estimated over all promoters, and each promoter's deviation from it is shrunk by `tau^2 / (tau^2 + se_p^2)` (`tau^2` = between-promoter variance). This reduces the noise from a few noisy RANDOM reporters while keeping promoter-specific shifts that are clearly supported by the data.

```yaml
RANDOM_NORM_METHOD: shrunk
```

Set `DEPTH_NORMALIZATION` to `True` to fit a TF-specific loess between cDNA RPM and
sequencing depth, normalize cDNA counts toward the median sequencing depth, and
generate the faceted diagnostic plot
`primetime_QC/cDNA_depth_correlations.pdf`. The default is `False`.

```yaml
DEPTH_NORMALIZATION: False
```

Set `BARCODE_OFFSET_CORRECTION` to `True` to estimate stable barcode-specific
offsets within each TF/promoter reporter from centered log2(cDNA/pDNA) activity
and correct the cDNA counts while preserving each sample's overall activity.
The diagnostic plot is written to
`primetime_QC/barcode_offset_correction.pdf`. The default is `False`.

```yaml
BARCODE_OFFSET_CORRECTION: False
```

Set `SHRINKAGE` to `True` (screen pipeline only) to shrink noisy low-count TF
measurements. For every TF and cDNA replicate, the log2 RPM of the barcode-summed
count is compared with the TF's mean over the control replicates. A robust loess of
that difference against the control level gives the replicate's own systematic shift
(for example a lower floor for weakly active reporters). The shrinkage target is
control mean + this trend, so only TF-specific noise is shrunk and the replicate-wide
shift is left to the other corrections. The weight is `w = tv / (tv + nv)`, where
`nv = 1 / (ln(2)^2 * count)` is the count-noise variance and `tv` is the median
excess variance around the target across all TFs and replicates. Low-count
measurements are shrunk hardest, which also damps real fold-changes of low-count TFs.
Weights and factors are written to `tmp_primetime/activity/shrinkage_fits.tsv`.
The default is `False`.

```yaml
SHRINKAGE: False
```

For the screen pipeline (`Snakefile_screen`), depth, shrinkage and barcode-offset factors are
estimated independently from the same post-filter, uncorrected counts and then
combined. Banana correction is also estimated from that uncorrected baseline before
being applied to the corrected counts. This makes the QC activity outputs and
comparative BCalm results independent of the order in which these corrections are
applied in code. Because the estimates no longer depend on previously corrected
counts, results can differ from runs made with earlier sequential fitting.

The banana correction (`BANANA_CORRECTION`) uses as reference only the
`BANANA_REFERENCE_N` (default `3`) control/DMSO replicates with the lowest pDNA
bleed-through; all other samples, including the remaining DMSO replicates, are
corrected against the mean of these. The selected replicates are listed in
`tmp_primetime/activity/banana_reference_samples.tsv`. Set `BANANA_REFERENCE_N: 0`
to use all control replicates. `primetime_QC/bleedthrough_estimation_all_reporters.pdf`
shows the bleed-through fits using all reporters (the main
`bleedthrough_estimation.pdf` uses negative controls only).

Set `BANANA_PER_PLATE` to `True` to select the `BANANA_REFERENCE_N` anchors separately
within each plate (`PLATES_FILE` plate ID) and correct every sample against the anchors of its own
plate (in both the QC activity table and the BCalm analysis). Samples whose replicate name
contains no plate ID, or on a plate without control replicates, use the global anchors. The default is `False`.

```yaml
BANANA_PER_PLATE: False
```

## Optional: Changing the library information

There are some additional parameters on the config file related to the barcodes used in the analysis

> Be aware that changing this parameters will impact the counting of the barcodes in the input reads. Only change this parameters if you know what you are doing.

Example:
```yaml
BARCODE_LENGTH: 12
BARCODE_DOWNSTREAM_SEQUENCE: CATCGTCGCATCCAAGAGGCTAGCTAACTA 
MAX_MISMATCH_DOWNSTREAM_SEQ: 3
BARCODE_ANNOTATION_FILE: misc/bc_annotation_prime.csv
EXPECTED_PDNA_COUNTS: misc/expected_pDNA_counts.txt
```

# Running **primetime**

After setting up your configuration file, you can run **primetime** with snakemake. 

```sh
snakemake --configfile <your_config.yaml> --use-conda --cores 10 --printshellcmds
```


# Output Files

**primetime** generates several output files during its execution:

### 1. **Quality Check outputs**
Inside the `primetime_QC` folder, several QC plots will be placed: 
- `barcode_correlations.pdf`: correlation of Log2(cDNA/pDNA) for different barcodes of each replicate.
- `replicate_correlations.pdf`: correlation of Log2(cDNA/pDNA) -- after averaging the different barcodes -- for each sample.
- `bleedthrough_estimation.pdf`: estimation of the ammount of pDNA bleedthough (percentage of cDNA counts coming from pDNA) foeach replicate.
- `distribution_of_BC_counts.pdf`: distribution of the counts from all the barcodes of each replicate.
- `expected_vs_observed_pDNA_counts.pdf`: correlation of your pDNA counts (observed) with the ones of our lab (expected).
- `read_counts.pdf`: total amount of reads coming from each replicate.

### 2. **Main Results**
- `primetime_results/primetime_results.txt`: Main result file, containing the adjusted p-value and the fold-change values for eacTF, as well as the activity of the TFs for each condition.
- `primetime_results/primetime_volcano.pdf`: Volcano plot of the differential activity results.
- `primetime_results/primetime_lollipop.pdf`: Lollipop plot showing the activity of each TF for both conditions, highlighting thdifferentially active ones.

### Additional files
**primetime** also saves some additional files in the `tmp_primetime` folder, such as the barcode counts, and the results of the barcode clustering.

These files provide a comprehensive overview of the TF activity analysis and can be used for further downstream analysis.