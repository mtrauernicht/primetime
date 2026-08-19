import sys
import inspect
import os
import re

__version__ = "1.0.0"
##########################################################################################
# primetime: TF reporter pipeline ########################################################
# ########################################################################################
# Description:
#
# This pipeline performs comprehensive analysis of TF prime reporter data,
# supporting multiple comparisons across diverse experimental designs, including 
# 96-well plates, drug screens, and large-scale studies.
# It processes fastq files, counts barcodes, clusters them, annotates them, and 
# performs differential TF activity analysis for all specified condition contrasts.
#
# Functionality:
#
# The pipeline is divided into the following steps:
# 1) Get barcode counts
#   - Counts the number of times each barcode appears in fastq files, using a
#     constant downstream sequence and allowing a configurable number of mismatches.
# 2) Cluster barcodes
#   - Clusters barcodes based on a distance of 1 to account for sequencing errors.
# 3) Annotate the barcodes
#   - Annotates barcodes using a provided reference file containing barcode information.
# 4) Get activity and QC plots
#   - Calculates TF activity and generates QC plots after annotation.
# 5) Run comparative analysis
#   - Performs differential TF activity analysis using the BCalm package, 
#     automatically handling all specified condition comparisons, including 
#     high-throughput formats like 96-well plates and drug screens.
# 6) Save final results
#   - Outputs comprehensive results for each comparison, including corrected 
#     activity and p-values for differential activity.
#

##########################################################################################
# Preamble functions #####################################################################
##########################################################################################

# Setting up directories and paths
filename = inspect.getframeinfo(inspect.currentframe()).filename
base_dir = os.path.dirname(os.path.abspath(filename))
scripts_dir = os.path.join(base_dir, "scripts")
conda_envs_dir = os.path.join(base_dir, "conda_envs")
output_dir = config["OUTPUT_DIRECTORY"]
comparisons_file = config['COMPARISONS_FILE']
samples_file = config['SAMPLES_FILE']
normalize_counts = "TRUE" if config.get("NORMALIZE_COUNTS", True) else "FALSE"


# Goodbye message

def format_config(key, value, space=0):
    """
    Helper function to format the config file in a nice way for the message
    """
    lines = []
    # Check if value is a dictionary
    if isinstance(value, dict):
        lines.append(f"{' ' * space}\033[1m{key}:\033[0m")
        for k, v in value.items():
            if isinstance(v, dict):
                lines.extend(format_config(k, v, space + 2))
            else:
                if isinstance(v, list):
                    lines.append(f"{' ' * space} - \033[1m{k}:\033[0m")
                    for i in v:
                        # Escape curly braces for Snakemake message formatting
                        escaped_i = str(i).replace("{", "{{").replace("}", "}}")
                        lines.append(f"{'   ' * space}{escaped_i}")
                else:
                    # Escape curly braces for Snakemake message formatting
                    escaped_v = str(v).replace("{", "{{").replace("}", "}}")
                    lines.append(f"{' ' * space} - \033[1m{k}:\033[0m {escaped_v}")
    else:
        # Escape curly braces for Snakemake message formatting
        escaped_value = str(value).replace("{", "{{").replace("}", "}}")
        lines.append(f"\033[1m{key}:\033[0m {escaped_value}")
    return lines


def create_goodbye_message(sample_replicate_files, comparisons_dict):
    """
    Create a formatted goodbye message with pipeline info and config
    """
    lines = []
    lines.append("\n"+"\033[0m=" * 80)
    lines.append("primetime")
    lines.append(f"Version {__version__}")
    lines.append("-" * 80)
    
    for key, value in config.items():
        lines.extend(format_config(key, value))
    
    # Show all the samples that were found
    lines.append("\033[1mSAMPLES:\033[0m")
    for sample in sample_replicate_files:
        lines.append(f" - {sample}")

    # Show all the comparisons that were made
    lines.append("\033[1mCOMPARISONS:\033[0m")
    for ref, contrasts in comparisons_dict.items():
        for contrast in contrasts:
            lines.append(f" - \033[4m{ref}\033[0m vs {contrast}")

    lines.append("=" * 80)
    return "\n".join(lines)

# Parse the samples and replicates on the samples_file
sample_replicate_files = {}
with open(samples_file, "r") as f:
    for line in f:
        row = line.strip().split("\t")
        sample = row[0]
        fastq = row[1]
        if sample not in sample_replicate_files:
            sample_replicate_files[sample] = {'1': fastq}
        else:
            # Get the replicate number based on how many files of this sample are there
            replicate = f"{len(sample_replicate_files[sample]) + 1}"
            sample_replicate_files[sample][replicate] = fastq


# Parse all the comparisons in the comparisons file
comparisons_dict = {}
with open(comparisons_file, 'r') as f:
    for line in f:
        row = line.strip().split('\t')
        # Reference condition
        ref_name = row[0]
        # Add the list of contrasts
        comparisons_dict[ref_name] = row[1:]

# Get the paths to the annotated dataframes (used in rule 'get_activity_and_qc_plots')
path_to_annotated_dfs = [
    os.path.join(output_dir, f"tmp_primetime/bc_counts/{sample}_{replicate}.cluster.annotated.txt")
    for sample in sample_replicate_files
    for replicate in sample_replicate_files[sample]
]

##########################################################################################
# Rules ##################################################################################
##########################################################################################

rule all:
    input:
        cDNA=os.path.join(output_dir, "tmp_primetime/activity/cDNA_counts.txt"),
        pDNA=os.path.join(output_dir, "tmp_primetime/activity/pDNA_counts.txt"),
        qc=os.path.join(output_dir, "primetime_QC/distribution_of_BC_counts.pdf"),
        heatmap_comparisons=os.path.join(output_dir, "primetime_results/heatmap_comparisons.pdf"),
        umap_plot=os.path.join(output_dir, "primetime_results/umap_conditions.pdf"),
        all_comparisons=[
            os.path.join(output_dir, f"primetime_results/output_data/{contrast}_vs_{ref}.txt")
            for ref, contrasts in comparisons_dict.items()
            for contrast in contrasts
        ],
        all_volcanos=[
            os.path.join(output_dir, f"primetime_results/volcano_plots/{contrast}_vs_{ref}.pdf")
            for ref, contrasts in comparisons_dict.items()
            for contrast in contrasts
        ],
        all_lollipops=[
            os.path.join(output_dir, f"primetime_results/lollipop_plots/{contrast}_vs_{ref}.pdf")
            for ref, contrasts in comparisons_dict.items()
            for contrast in contrasts
        ],
        all_top_tf_plots=[
            os.path.join(output_dir, f"primetime_results/top_tf_barcode_counts/{contrast}_vs_{ref}.pdf")
            for ref, contrasts in comparisons_dict.items()
            for contrast in contrasts
        ]
    message:
        create_goodbye_message(sample_replicate_files, comparisons_dict)

##########################################################################################

rule get_barcode_counts:
    input:
        lambda wildcards: sample_replicate_files[wildcards.sample][wildcards.replicate]
    output:
        txt=os.path.join(output_dir, "tmp_primetime/bc_counts/{sample}_{replicate}.txt"),
        stats=os.path.join(
            output_dir, "tmp_primetime/bc_counts/{sample}_{replicate}.stats"
        ),
        invalid_bc_file=os.path.join(output_dir, "tmp_primetime/bc_counts/{sample}_{replicate}.invalid.txt")
    params:
        script=os.path.join(scripts_dir, "primetime_get_barcode_counts.py"),
        bc_length=config["BARCODE_LENGTH"],
        bc_downstream_seq=config["BARCODE_DOWNSTREAM_SEQUENCE"],
        max_mismatch=config["MAX_MISMATCH_DOWNSTREAM_SEQ"],
    threads: 20
    conda:
        os.path.join(conda_envs_dir, "test_bc_counts.yaml")
    shell:
        """
        python {params.script} \
--fastq {input} \
--bc_length {params.bc_length} \
--bc_downstream_seq {params.bc_downstream_seq} \
--max_mismatch {params.max_mismatch} \
--threads {threads} \
--invalid_bc_file {output.invalid_bc_file} \
> {output.txt} 2> {output.stats}
        """

##########################################################################################

# 2) Cluster barcodes
rule cluster_barcodes:
    input:
        os.path.join(output_dir, "tmp_primetime/bc_counts/{sample}_{replicate}.txt"),
    output:
        clustered=os.path.join(
            output_dir, "tmp_primetime/bc_counts/{sample}_{replicate}.cluster.txt"
        ),
        tmp=temp(os.path.join(
            output_dir, "tmp_primetime/bc_counts/{sample}_{replicate}.input"
        ))
    params:
        amount_spikes=100000,
        bc_annotation=config["BARCODE_ANNOTATION_FILE"],
    conda:
        os.path.join(conda_envs_dir, "starcode.yaml")
    threads: 10
    shell:
        """
        # First, repeat every BC on the annotation 100000 times
        tail -n +2 {params.bc_annotation} | awk -v FS=',' '{{for(i=1;i<={params.amount_spikes};i++) print $1}}' > {output.tmp}

        # Now, concatenate the input and this and call starcode with sphere clustering
        cat {input} {output.tmp} | starcode \
        --threads {threads} \
        --print-clusters \
        --quiet \
        --sphere \
        --dist 1 | \
        sort -k1,1 > {output.clustered}
        """
# rule cluster_barcodes:
#     input:
#         os.path.join(output_dir, "tmp_primetime/bc_counts/{sample}_{replicate}.txt"),
#     output:
#         os.path.join(
#             output_dir, "tmp_primetime/bc_counts/{sample}_{replicate}.cluster.txt"
#         ),
#     conda:
#         os.path.join(conda_envs_dir, "starcode.yaml")
#     threads: 10
#     shell:
#         """
#         starcode --threads {threads} --print-clusters -i {input} --dist 1 2> /dev/null | \
# sort -k1,1 > {output}
#         """

##########################################################################################

# 3) Annotate the barcodes
rule annotate_barcodes:
    input:
        clustered = os.path.join(
            output_dir, "tmp_primetime/bc_counts/{sample}_{replicate}.cluster.txt"
        ),
        stats = os.path.join(output_dir, "tmp_primetime/bc_counts/{sample}_{replicate}.stats"),
    output:
        os.path.join(
            output_dir,
            "tmp_primetime/bc_counts/{sample}_{replicate}.cluster.annotated.txt",
        ),
    params:
        script=os.path.join(scripts_dir, "primetime_annotate_barcodes.R"),
        bc_annotation=config["BARCODE_ANNOTATION_FILE"],
        amount_spikes=100000
    conda:
        os.path.join(conda_envs_dir, "r_plotting.yaml")
    shell:
        """
        Rscript {params.script} \
--bc_counts {input.clustered} \
--stats {input.stats} \
--amount_spikes {params.amount_spikes} \
--bc_annotation {params.bc_annotation} \
--sample {wildcards.sample} \
--replicate_number {wildcards.replicate} \
--output {output}
        """

##########################################################################################

rule get_activity_and_qc_plots:
    input:
        dfs=path_to_annotated_dfs
    output:
        bc_corr=os.path.join(output_dir, "primetime_QC/barcode_correlations.pdf"),
        replicate_corr=os.path.join(
            output_dir, "primetime_QC/replicate_correlations.pdf"
        ),
        expected_vs_obs=os.path.join(
            output_dir, "primetime_QC/expected_vs_observed_pDNA_counts.pdf"
        ),
        read_counts=os.path.join(output_dir, "primetime_QC/read_counts.pdf"),
        bc_counts=os.path.join(output_dir, "primetime_QC/distribution_of_BC_counts.pdf"),
        control_bc_corr=os.path.join(output_dir, "primetime_QC/control_barcode_correlations.pdf"),
        control_bc_neg_corr_stats=os.path.join(output_dir, "primetime_QC/control_barcode_residual_negative_correlation_stats.tsv"),
        bleedthrough=os.path.join(
            output_dir, "primetime_QC/bleedthrough_estimation.pdf"
        ),
        cDNA=os.path.join(output_dir, "tmp_primetime/activity/cDNA_counts.txt"),
        pDNA=os.path.join(output_dir, "tmp_primetime/activity/pDNA_counts.txt"),
        barcode_activity=os.path.join(output_dir, "tmp_primetime/activity/barcode_activity.txt"),
        read_count_summary=os.path.join(output_dir, "tmp_primetime/activity/read_count_per_sample.tsv"),
    params:
        script=os.path.join(scripts_dir, "primetime_get_activity_and_qc_plots.R"),
        df_basedir=os.path.join(output_dir, "tmp_primetime/bc_counts"),
        expected_pdna_counts=config["EXPECTED_PDNA_COUNTS"],
        plots_basedir=os.path.join(output_dir, "primetime_QC/"),
        activity_basedir=os.path.join(output_dir, "tmp_primetime/activity"),
        barcode_activity_output=os.path.join(output_dir, "tmp_primetime/activity/barcode_activity.txt"),
    conda:
        os.path.join(conda_envs_dir, "r_plotting.yaml")
    shell:
        """
        mkdir -p {params.plots_basedir}
        mkdir -p {params.activity_basedir}

        Rscript {params.script} \
        --list_of_annotated_files "{input}" \
        --design {samples_file} \
        --plots_basedir {params.plots_basedir} \
        --activity_basedir {params.activity_basedir} \
        --expected_pdna {params.expected_pdna_counts} \
        --cdna_output {output.cDNA} \
        --barcode_activity_output {params.barcode_activity_output}
        """

rule run_comparative_analysis:
    input:
        cDNA=os.path.join(output_dir, "tmp_primetime/activity/cDNA_counts.txt"),
        pDNA=os.path.join(output_dir, "tmp_primetime/activity/pDNA_counts.txt"),
    output:
        txt=os.path.join(output_dir, "primetime_results/output_data/{contrast}_vs_{ref}.txt"),
        plots=os.path.join(output_dir, "primetime_results/volcano_plots/{contrast}_vs_{ref}.pdf"),
        loli=os.path.join(output_dir, "primetime_results/lollipop_plots/{contrast}_vs_{ref}.pdf"),
        circular_loli=os.path.join(output_dir, "primetime_results/circular_lollipop_plots/{contrast}_vs_{ref}.pdf"),
    params:
        script=os.path.join(scripts_dir, "primetime_run_comparative_analysis_bcalm.R"),
        out_basedir=os.path.join(output_dir, "tmp_primetime/activity"),
        pval_threshold=config["PVALUE_THRESHOLD"],
        plot_output_dir=os.path.join(output_dir, "primetime_results"),
        contrast_condition=lambda wildcards: wildcards.contrast,
        num_replicates_contrast=lambda wildcards: len(sample_replicate_files[wildcards.contrast]),
        reference_condition=lambda wildcards: wildcards.ref,
        num_replicates_reference=lambda wildcards: len(sample_replicate_files[wildcards.ref]),
        split_by_promoter=config.get("SPLIT_COMPARATIVE_ANALYSIS_BY_PROMOTER", True),
        normalize=normalize_counts
    conda:
        os.path.join(conda_envs_dir, "comparative_analysis.yaml")
    threads: 1
    shell:
        """
        Rscript {params.script} \
        --cdna {input.cDNA} \
        --pdna {input.pDNA} \
        --output {params.out_basedir} \
        --pval_threshold {params.pval_threshold} \
        --contrast_condition {params.contrast_condition} \
        --num_replicates_contrast {params.num_replicates_contrast} \
        --reference_condition {params.reference_condition} \
        --num_replicates_reference {params.num_replicates_reference} \
        --plot_output {params.plot_output_dir} \
        --split_by_promoter {params.split_by_promoter} \
        --normalize {params.normalize}
        """

##########################################################################################

rule get_bc_counts_for_top_tfs:
    input:
        result=os.path.join(output_dir, "primetime_results/output_data/{contrast}_vs_{ref}.txt"),
        cdna=os.path.join(output_dir, "tmp_primetime/activity/cDNA_counts.txt"),
        pdna=os.path.join(output_dir, "tmp_primetime/activity/pDNA_counts.txt"),
    output:
        plot=os.path.join(output_dir, "primetime_results/top_tf_barcode_counts/{contrast}_vs_{ref}.pdf"),
        design=temp(os.path.join(output_dir, "tmp_primetime/activity/design/{contrast}_vs_{ref}.txt")),
    params:
        script=os.path.join(scripts_dir, "primetime_get_bc_counts_for_top_tfs.R"),
        reference_condition=lambda wildcards: wildcards.ref,
        contrast_condition=lambda wildcards: wildcards.contrast,
        num_replicates_reference=lambda wildcards: len(sample_replicate_files[wildcards.ref]),
        num_replicates_contrast=lambda wildcards: len(sample_replicate_files[wildcards.contrast]),
    conda:
        os.path.join(conda_envs_dir, "r_plotting.yaml")
    shell:
        """
        mkdir -p $(dirname {output.plot})
        mkdir -p $(dirname {output.design})

        echo -e "replicate\ttreatment\tsample\tpDNA" > {output.design}
        for i in $(seq 1 {params.num_replicates_reference}); do
            echo -e "{params.reference_condition}_$i\t{params.reference_condition}\t{params.reference_condition}\tFalse" >> {output.design}
        done
        for i in $(seq 1 {params.num_replicates_contrast}); do
            echo -e "{params.contrast_condition}_$i\t{params.contrast_condition}\t{params.contrast_condition}\tFalse" >> {output.design}
        done

        Rscript {params.script} \
        --result {input.result} \
        --cdna {input.cdna} \
        --pdna {input.pdna} \
        --output {output.plot} \
        --design {output.design} \
        --reference_condition {params.reference_condition}
        """

##########################################################################################

rule get_heatmap_of_conditions:
    input:
        results=expand(
            os.path.join(output_dir, f"primetime_results/output_data/{contrast}_vs_{ref}.txt")
            for ref, contrasts in comparisons_dict.items()
            for contrast in contrasts
        ),
        negative_correlation_stats=os.path.join(output_dir, "primetime_QC/control_barcode_residual_negative_correlation_stats.tsv"),
        read_count_summary=os.path.join(output_dir, "tmp_primetime/activity/read_count_per_sample.tsv")
    output:
        os.path.join(output_dir, "primetime_results/heatmap_comparisons.pdf"),
    params:
        script=os.path.join(scripts_dir, "primetime_get_heatmap_of_conditions.R"),
    conda:
        os.path.join(conda_envs_dir, "r_plotting.yaml")
    shell:
        """
        if [ $(echo {input.results} | wc -w) -le 1 ]; then
            echo "Not enough input files for heatmap, skipping."; touch {output};
        else
            results=$(echo {input.results} | tr ' ' ',')
            Rscript {params.script} \
            --results $results \
            --negative-correlation-stats {input.negative_correlation_stats} \
            --read-count-summary {input.read_count_summary} \
            --output {output} > /dev/null 2>&1
        fi
        """

##########################################################################################
# 8) Get UMAP plot of all conditions
rule get_umap_plot:
    input:
        results=expand(
            os.path.join(output_dir, f"primetime_results/output_data/{contrast}_vs_{ref}.txt")
            for ref, contrasts in comparisons_dict.items()
            for contrast in contrasts
        )
    output:
        os.path.join(output_dir, "primetime_results/umap_conditions.pdf"),
    params:
        r_script=os.path.join(scripts_dir, "primetime_get_pca_plot.R"),
    conda:
        os.path.join(conda_envs_dir, "r_plotting.yaml")
    shell:
        """
        FILELIST=$(mktemp)
        printf '%s\n' {input.results} > "$FILELIST"
        Rscript {params.r_script} --results-file "$FILELIST" --output {output}
        rm "$FILELIST"
        """