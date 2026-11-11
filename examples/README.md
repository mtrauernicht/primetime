# Primetime example data

Small, self-contained datasets to run both pipelines end to end. All paths in the configs are
relative, so run the commands from the respective example directory with the `primetime` conda
environment activated (`conda activate /DATA/shared/conda_envs/primetime`).

| Example | Pipeline | Size | Runtime (8 cores) |
|---|---|---|---|
| `regular/` | `Snakefile` | ~38 MB | ~13 min |
| `screen/`  | `Snakefile_screen` | ~14 MB | ~5 min |

## `regular/`: DOT1L WT vs KO T cells (two replicates)

Subset of the DOT1L comparison run: 6 conditions (WT/KO x unactivated, aCD3/aCD28, PMA/Ionomycin),
two replicates each, plus pDNA (subsampled to 200k reads). FASTQs are single-end reads of the reporter
barcodes.

```bash
cd examples/regular
snakemake --configfile config.yaml -s ../../Snakefile --cores 8
```

Files: `samples.txt` (condition, FASTQ), `comparisons.txt` (reference followed by contrasts, tab-separated),
`config.yaml`.

## `screen/`: K562 drug-screen subset

13 wells of one 384-well plate (4 DMSO, 3 PMA, 3 Forskolin, 3 Ionomycin; each well subsampled to 60k reads)
and pDNA (300k reads). The reads are already demultiplexed per well
(`demultiplexed_fastq/<plate>_<well>.fastq.gz`), so `RUN_DEMULTIPLEXING` is `False` and the raw plate
FASTQ listed in `plates_to_file.txt` is only a placeholder that is never read.

```bash
cd examples/screen
snakemake --configfile config_screen.yaml -s ../../Snakefile_screen --cores 8
```

Files: `plates_to_file.txt` (plate, FASTQ, well specification), `well_specification.csv`
(well, sample, replicate), `comparisons.txt`, `config_screen.yaml`.

## Expected output

Written to `OUTPUT/` in each example directory: `primetime_QC/` (QC plots), `primetime_results/`
(`output_data/*.txt`, volcano/lollipop plots, `heatmap_comparisons.pdf`, `umap_conditions.pdf`) and
`tmp_primetime/` (barcode counts and activity tables). Each comparison table has 131 TF rows.
Check with a dry run first (`-n`). Delete `OUTPUT/` and `.snakemake/` to start from scratch.
