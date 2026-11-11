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

Set `NORMALIZE_COUNTS` to `True` to apply RANDOM-reporter normalization, or to `False` to skip it.

Example:
```yaml
NORMALIZE_COUNTS: True
```

When enabled, normalization is always count-level: counts are scaled before BCalm separately for each promoter and replicate, using that promoter's RANDOM-reporter count sum as the reference.

BCalm's additional global library-size normalization is disabled so it does not undo the promoter-specific scaling. One QC PDF, `primetime_QC/count_level_normalization.pdf`, shows the log2 shift per condition and promoter relative to plain library-size normalization (log2 of the RANDOM share of pDNA over the RANDOM share of cDNA; positive = corrected up). The scaling target (1e7) cancels in this shift.

Both pipelines also produce `primetime_QC/raw_random_reporter_rpm.pdf`. For each comparison, the shaded reference condition is shown on the left, followed by its contrasts, with promoters in separate rows. Pages show at most 15 conditions, repeating the reference on each page. Each point represents a RANDOM negative-control barcode in one cDNA replicate, coloured by reporter. RPM is the original barcode count divided by the sum of all annotated barcode counts in that replicate, multiplied by one million. The plot is made before filtering or corrections, and therefore includes replicates subsequently excluded by `MIN_READS_PER_REPLICATE`.

For each barcode, `log2(RPM + 1)` is compared to the median across reference replicates of that same barcode, reporter and promoter. The black diamond is the median of these barcode-replicate differences per condition and promoter; zero means no change relative to the reference. The added 1 permits zero counts and is used only for this plot. These are descriptive cDNA RPM changes, not pDNA-normalized activities or statistical tests. Measurements without a matching reference barcode are omitted with a warning.

The regular pipeline writes barcode plots for each comparison under `primetime_results/top_tf_barcode_counts/before_normalization/` and `after_normalization/`. The before plot retains the existing log2(cDNA RPM / pDNA RPM) view from the shared count tables, after upstream barcode and replicate filtering but before comparison-specific banana and RANDOM count normalization. The after plot shows the actual per-barcode BCalm input, log2((cDNA + 1) / (pDNA + 1)), including enabled banana and RANDOM normalization, or BCalm's global count normalization when RANDOM normalization is disabled. Thus the two plots use different library-scaling and pseudocount conventions; their absolute y-values are not directly interchangeable. TF significance and logFC labels in both plots come from the same comparison model; shrinkage is applied to TF-level estimates, not barcode counts.

Set `DEPTH_NORMALIZATION` (screen pipeline only) to `True` to fit a TF-specific loess between cDNA RPM and
sequencing depth, normalize cDNA counts toward the median sequencing depth, and
generate the faceted diagnostic plot
`primetime_QC/cDNA_depth_correlations.pdf`. The default is `False`.

```yaml
DEPTH_NORMALIZATION: False
```

Set `BARCODE_OFFSET_CORRECTION` (screen pipeline only) to `True` to estimate stable barcode-specific
offsets within each TF/promoter reporter from centered log2(cDNA/pDNA) activity
and correct the cDNA counts while preserving each sample's overall activity.
The diagnostic plot is written to
`primetime_QC/barcode_offset_correction.pdf`. The default is `False`.

```yaml
BARCODE_OFFSET_CORRECTION: False
```

Set `SHRINKAGE` to `True` to add empirical-Bayes shrinkage of
the TF log2 fold-changes. It runs after BCalm and never touches the counts or the
statistics: p-values, FDR and the `Upregulated`/`Downregulated` calls stay on the
unshrunk logFC. For every contrast, a zero-centred mixture-of-normals prior
(adaptive-shrinkage style, learned from the logFC and standard errors of all non-RANDOM
TFs) gives a posterior-mean `logFC_shrunk`. Noisy estimates (large SE) are pulled towards 0,
precise ones are left almost unchanged. The `output_data` tables gain the columns `logFC_se`,
`logFC_shrunk`, `shrink_weight` (1 = no shrinkage) and `<contrast>_shrunk` (activity). In the
regular pipeline, standard plots remain based on unshrunk values; the screen pipeline also
produces shrunk plot variants. The BCalm standard errors are conservative (barcodes are shared
between conditions but modelled as independent groups), so the SEs are first rescaled by the spread of
the z-scores (logFC / SE) of the RANDOM reporters, which are true nulls (factor in `se_scale`, capped at 1).
`primetime_QC/shrinkage_effect.pdf` (and a per-contrast `.tsv`) summarises the effect: shrunk vs unshrunk
logFC, weight vs SE and |logFC|, weights and SE scale per contrast, and how many significant hits
fall below the logFC cutoff after shrinkage. The default is `False`.

```yaml
SHRINKAGE: False
```

For the screen pipeline (`Snakefile_screen`), depth and barcode-offset factors are
estimated independently from the same post-filter, uncorrected counts and then
combined. Banana correction is also estimated from that uncorrected baseline before
being applied to the corrected counts. This makes the QC activity outputs and
comparative BCalm results independent of the order in which these corrections are
applied in code. Because the estimates no longer depend on previously corrected
counts, results can differ from runs made with earlier sequential fitting.

The banana correction (`BANANA_CORRECTION`) is done separately for every comparison. Within the
reference condition of a comparison (first column of `COMPARISONS_FILE`), the
`BANANA_REFERENCE_N` (default `3`) replicates with the lowest pDNA bleed-through (per plate) are the
anchors; the contrast replicates and the remaining reference replicates are corrected against the mean
of these anchors, and the anchors themselves are left untouched. Set `BANANA_REFERENCE_N: 0` to use all
replicates of the reference condition. The bleed-through of every reference replicate and the selected
anchors are listed in `tmp_primetime/activity/banana_reference_samples.tsv`. The QC plots in
`primetime_QC/` compare each replicate with the anchors of its own comparison.
`primetime_QC/bleedthrough_estimation_all_reporters.pdf` shows the bleed-through fits using all reporters
(the main `bleedthrough_estimation.pdf` uses negative controls only).

Anchors are selected separately within each plate (`PLATES_FILE` plate ID) and every sample is corrected
against the anchors of its own plate. Samples whose replicate name contains no plate ID, or on a plate
without anchors, use the anchors of the whole reference condition.

When `BANANA_CORRECTION` is `True`, the banana is fitted on counts that already carry the depth,
barcode-offset corrections (not on the uncorrected baseline). Depth and banana both
remove the low-activity floor of high-depth samples, so this way the banana only removes what is
left and the same bias is not corrected twice. The result is independent of the order of the YAML flags.

`MIN_READS_PER_REPLICATE` (default `0`, i.e. off) excludes cDNA replicates with fewer reads than this value from all downstream analysis. `MAX_BLEEDTHROUGH` (default `0.5`, i.e. 50%; set `-1` to disable) excludes cDNA replicates whose estimated bleed-through slope is strictly greater than the cutoff. The slope is the linear fit of negative-control cDNA RPM against mean pDNA RPM, with an intercept, after the read-count filter but before barcode filtering or corrections. An unestimable slope stops the run. `tmp_primetime/activity/bleedthrough_filter_estimates.tsv` records all tested replicates and decisions; `excluded_high_bleedthrough_replicates.tsv` lists the exclusions. The later bleed-through diagnostic is fitted after barcode filtering and may differ from the filter estimate. Comparisons that lose every replicate of one side are written with empty outputs; if all cDNA replicates are excluded, the QC step stops with an error.

The bleed-through fit uses the pipeline's RPM convention, which adds one pseudocount per annotated barcode before calculating each library total. The cutoff is active by default, including for existing configs that do not specify it; set `MAX_BLEEDTHROUGH: -1` to retain the previous behavior.

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

`primetime_QC/banana_correction.pdf` shows the banana correction as one page per corrected replicate (M = replicate minus anchor activity vs anchor activity, before then after correction). Anchor replicates are not corrected and have no page.
