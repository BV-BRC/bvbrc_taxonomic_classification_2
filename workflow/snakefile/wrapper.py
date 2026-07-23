#!/opt/patric-common/runtime/bin/python3
import csv
import json 
import os
import re
import shlex
import shutil
import subprocess
import sys

def check_input_fastqs(input_path, sample_id, gzip_targets):
    # start new 
    sample_id = clean_sample_id(sample_id)

    if input_path.endswith("fastq.gz") or input_path.endswith("fq.gz"):
        return sample_id, ""
    
    elif input_path.endswith(".fastq") or input_path.endswith(".fq"):
        zipped_path = input_path + ".gz"
        # check if the zipped file already exists in the list of input files 
        if os.path.exists(zipped_path) == True:
            msg =f"zipped file exists using {input_path}.gz instead \n"
            sys.stderr.write(msg)
            return sample_id, ""

        else:
            msg = f"{input_path} identified as unzipped, zipping for analysis \n"
            sys.stderr.write(msg)
            gzip_targets.append(input_path)
            return sample_id, ".gz"
    else:
        msg = f"Warning: {input_path} does not end in 'fastq'/'fq' or 'fastq.gz'/'fq.gz'. \n"
        sys.stderr.write(msg)
        if "gz" in input_path:
            msg = f"{input_path} contains 'gz', treating as already zipped \n"
            sys.stderr.write(msg)
            return sample_id, ""
        else:
            msg = f"{input_path} identified as unzipped, zipping for analysis \n"
            sys.stderr.write(msg)
            gzip_targets.append(input_path)
            return sample_id, ".gz"


def clean_sample_id(s):
    return re.sub(r"[^a-zA-Z0-9_]", "", s)

def format_inputs(raw_params):
    input_dict = json.loads(raw_params)
    return input_dict


def load_hisat_indicies(input_dict, hisat_indicies_path, input_dir):
    # get the refrence genome from parameters
    if input_dict["host_genome"] == "no_host":
        msg = "No host genome selected. Not getting hisat indices \n"
        sys.stderr.write(msg)
    else:
        host_indicies_file = os.path.join(hisat_indicies_path, input_dict["host_genome"])
        staging_indices_dir = os.path.join(input_dir, "indices_dir")
        # shutil copy requires filename
        hisat_indicies= os.path.join(input_dir, "hisat_indicies")
        # cmd 1
        os.makedirs(staging_indices_dir, exist_ok=True )
        # cmd 2
        shutil.copy(host_indicies_file, hisat_indicies)
        # cmd 3
        tar_cmd = ["tar", "-xvf", "{}".format(hisat_indicies)]
        subprocess.run(tar_cmd)
    return


def get_paired_and_single_sample_ids(input_dict):
    paired_sample_ids = [clean_sample_id(item["sample_id"]) for item in input_dict.get("paired_end_libs", [])]
    single_sample_ids = [clean_sample_id(item["sample_id"]) for item in input_dict.get("single_end_libs", [])]
    return paired_sample_ids, single_sample_ids


# Build the per-sample/per-step (Complete/Fail) rows for JobStatus.html.
# Job-level steps (steps that produce one output for the whole job rather than
# per sample) are reported under the pseudo-sample "Multisample".
def build_job_status_rows(all_sample_ids, output_dir, input_dict, paired_sample_ids):
    analysis_type = input_dict["analysis_type"]
    host_genome = input_dict.get("host_genome", "no_host")
    multi_sample = len(all_sample_ids) > 1
    rows = []

    for sample_name in all_sample_ids:
        is_paired = sample_name in paired_sample_ids

        if is_paired:
            fastqc_ok = all(
                os.path.isfile(f"{output_dir}/{sample_name}/fastqc_results/raw_reads/raw_{sample_name}_R{read_num}_fastqc.html")
                for read_num in (1, 2)
            )
        else:
            # se_fastq_processing_16s uses "raw_read" (singular); every other pipeline uses "raw_reads"
            raw_reads_dir = "raw_read" if analysis_type == "16S" else "raw_reads"
            fastqc_ok = os.path.isfile(f"{output_dir}/{sample_name}/fastqc_results/{raw_reads_dir}/raw_{sample_name}_fastqc.html")
        rows.append((sample_name, "FastQC (raw reads)", fastqc_ok))

        if analysis_type in ("pathogen", "microbiome") and host_genome != "no_host":
            if is_paired:
                host_removal_ok = all(
                    os.path.isfile(f"{output_dir}/{sample_name}/hisat2_results/{sample_name}_host_removed_R{read_num}.fastq.gz")
                    for read_num in (1, 2)
                )
            else:
                host_removal_ok = os.path.isfile(f"{output_dir}/{sample_name}/hisat2_results/{sample_name}_host_removed.fastq.gz")
            rows.append((sample_name, "Host Removal (HISAT2)", host_removal_ok))

        if analysis_type == "16S":
            if is_paired:
                trim_ok = all(
                    os.path.isfile(f"{output_dir}/{sample_name}/trimmed_reads/{sample_name}_R{read_num}_trimmed.fastq.gz")
                    for read_num in (1, 2)
                )
            else:
                trim_ok = os.path.isfile(f"{output_dir}/{sample_name}/trimmed_read/{sample_name}_trimmed.fastq.gz")
            rows.append((sample_name, "Adapter Trimming (Trim Galore)", trim_ok))

        kraken_ok = os.path.isfile(f"{output_dir}/{sample_name}/kraken_output/{sample_name}_k2_report.txt")
        rows.append((sample_name, "Taxonomic Classification (Kraken2)", kraken_ok))

        if analysis_type in ("microbiome", "16S"):
            bracken_ok = os.path.isfile(f"{output_dir}/{sample_name}/bracken_output/{sample_name}_bracken_output.txt")
            rows.append((sample_name, "Abundance Estimation (Bracken)", bracken_ok))

            alpha_ok = os.path.isfile(f"{output_dir}/{sample_name}/{sample_name}_alpha_diversity.csv")
            rows.append((sample_name, "Alpha Diversity", alpha_ok))

        krona_ok = os.path.isfile(f"{output_dir}/{sample_name}_krona.html")
        rows.append((sample_name, "Krona Plot", krona_ok))

    if multi_sample:
        rows.append(("Multisample", "Multisample Comparison", os.path.isfile(f"{output_dir}/multisample_comparison.html")))
        rows.append(("Multisample", "Multisample Krona Plot", os.path.isfile(f"{output_dir}/multisample_krona.html")))
        if analysis_type in ("microbiome", "16S"):
            rows.append(("Multisample", "Beta Diversity", os.path.isfile(f"{output_dir}/beta_diversity.csv")))

    if analysis_type in ("microbiome", "16S"):
        rows.append(("Multisample", "Merged Alpha Diversity", os.path.isfile(f"{output_dir}/alpha_diversity.csv")))

    sankey_file = "multisample_sankey.html" if multi_sample else "sankey.html"
    rows.append(("Multisample", "Sankey Plot", os.path.isfile(f"{output_dir}/{sankey_file}")))

    multiqc_report = "Taxonomic-Classification-Service-BVBRC_multiqc_report.html"
    rows.append(("Multisample", "MultiQC Report", os.path.isfile(f"{output_dir}/{multiqc_report}")))

    return rows


def write_job_status_html(rows, output_dir):
    table_rows = "\n".join(
        f'<tr><td>{sample_name}</td><td>{step}</td>'
        f'<td><span class="badge {"complete" if ok else "fail"}">{"&#10003; Complete" if ok else "&#10007; Fail"}</span></td></tr>'
        for sample_name, step, ok in rows
    )
    num_failed = sum(1 for _, _, ok in rows if not ok)
    summary_class = "complete" if num_failed == 0 else "fail"
    summary_text = "All steps completed successfully" if num_failed == 0 else f"{num_failed} step(s) need review"

    html = f"""<html>
<head>
<title>Taxonomic Classification Job Status</title>
<style>
body {{
    font-family: -apple-system, "Segoe UI", Helvetica, Arial, sans-serif;
    background: #f4f6f8;
    color: #222;
    margin: 0;
    padding: 2.5em 1em;
}}
.container {{
    max-width: 900px;
    margin: 0 auto;
    background: #fff;
    border-radius: 10px;
    box-shadow: 0 2px 12px rgba(0, 0, 0, 0.08);
    overflow: hidden;
}}
.banner {{
    background: linear-gradient(135deg, #1f3a5f, #2f6690);
    color: #fff;
    padding: 1.5em 2em;
}}
.banner h1 {{
    margin: 0 0 0.2em 0;
    font-size: 1.5em;
}}
.banner .summary {{
    display: inline-block;
    margin-top: 0.5em;
    padding: 0.25em 0.8em;
    border-radius: 999px;
    font-size: 0.9em;
    font-weight: 600;
}}
.banner .summary.complete {{ background: #2e7d46; color: #fff; }}
.banner .summary.fail {{ background: #b3261e; color: #fff; }}
table {{
    border-collapse: collapse;
    width: 100%;
}}
th, td {{
    padding: 0.6em 1.2em;
    text-align: left;
}}
thead th {{
    background: #eef1f4;
    color: #333;
    font-size: 0.85em;
    text-transform: uppercase;
    letter-spacing: 0.03em;
    border-bottom: 2px solid #dde2e7;
}}
tbody tr:nth-child(odd) {{ background: #fafbfc; }}
tbody tr:hover {{ background: #eef4fa; }}
td {{ border-bottom: 1px solid #eee; }}
.badge {{
    display: inline-block;
    padding: 0.2em 0.7em;
    border-radius: 999px;
    font-size: 0.85em;
    font-weight: 600;
}}
.badge.complete {{ background: #e3f5e9; color: #1e7d3c; }}
.badge.fail {{ background: #fbe6e5; color: #b3261e; }}
.steps {{
    padding: 1.5em 2em 2em 2em;
    border-top: 1px solid #eee;
}}
.steps h2 {{
    font-size: 1.1em;
    margin-bottom: 0.8em;
}}
.steps dl {{
    margin: 0;
}}
.steps dt {{
    font-weight: 600;
    margin-top: 1em;
}}
.steps dd {{
    margin: 0.2em 0 0 0;
    color: #444;
    line-height: 1.4em;
}}
.steps dd.advice {{
    color: #444;
    background: #fdf3f2;
    border-left: 4px solid #b3261e;
    border-radius: 4px;
    padding: 0.6em 1em;
    margin-top: 0.5em;
}}
.steps dd.advice > *:first-child {{ margin-top: 0; }}
.steps dd.advice > *:last-child {{ margin-bottom: 0; }}
.steps dd h3 {{
    font-size: 1em;
    margin: 0.6em 0 0.3em 0;
    color: #222;
}}
.steps dd ul {{
    margin: 0.4em 0;
    padding-left: 1.2em;
}}
.steps dd li {{
    margin-bottom: 0.3em;
}}
a {{ color: #2f6690; }}
</style>
</head>
<body>
<div class="container">
<div class="banner">
<h1>Taxonomic Classification Job Status</h1>
<span class="summary {summary_class}">{summary_text}</span>
</div>
<table>
<thead><tr><th>Sample</th><th>Step</th><th>Status</th></tr></thead>
<tbody>
{table_rows}
</tbody>
</table>
<div class="steps">
<h2>Pipeline Overview</h2>
<dl>

<dt>FastQC (raw reads)</dt>
<dd>Runs quality control on the raw input reads.</dd>
<dd class="advice">If this step failed, please review: check that your reads were fully uploaded to the workspace (see this <a href="https://www.bv-brc.org/docs/quick_references/workspaces/data_upload.html" target="_blank">data upload guide</a>). A common cause of failure is uploading a FASTA file where a FASTQ file was expected &mdash; see this <a href="https://sequenceserver.com/blog/fasta-vs-fastq-formats/" target="_blank">guide on FASTA vs. FASTQ formats</a>.</dd>

<dt>Host Removal (HISAT2)</dt>
<dd>Removes host-derived reads by aligning to the selected host genome (only runs if a host genome was specified).</dd>
<dd class="advice">If this step failed, please review: confirm the host genome you selected actually matches the organism your reads were sequenced from.</dd>

<dt>Adapter Trimming (Trim Galore)</dt>
<dd>Trims sequencing adapters from the reads before classification. This step is specific to the 16S pipeline, which does not include a host-removal step.</dd>
<dd class="advice">If this step failed, please review: the FastQC results for the raw reads &mdash; trimming failures are usually caused by the same malformed or empty input issues that would show up there.</dd>

<dt>Taxonomic Classification (Kraken2)</dt>
<dd>Assigns taxonomic labels to each read against the reference database.</dd>
<dd class="advice">If this step failed, please review: the FastQC results for an explanation.</dd>

<dt>Abundance Estimation (Bracken)</dt>
<dd>Re-estimates species/genus-level abundance from the Kraken2 report.</dd>
<dd class="advice">If this step failed, please review: consider rerunning your samples with the Species Identification pipeline, which relies solely on the Kraken2 results for all of its outputs.</dd>

<dt>Alpha Diversity / Merged Alpha Diversity</dt>
<dd>Calculates within-sample diversity statistics from the Bracken results.</dd>
<dd class="advice">If this step failed, please review: the raw Kraken2 files available in your job results at &lt;sample_name&gt;/kraken_output/&lt;sample_name&gt;_k2_report.txt and &lt;sample_name&gt;/kraken_output/&lt;sample_name&gt;_k2_output.txt.</dd>

<dt>Beta Diversity</dt>
<dd>Calculates between-sample diversity comparisons (only runs when more than one sample is present).</dd>
<dd class="advice">If this step failed, please review: confirm that Bracken completed successfully for every sample in the job &mdash; this step depends on all of them.</dd>

<dt>Krona Plot / Multisample Krona Plot</dt>
<dd>Interactive charts visualizing the taxonomic composition of each sample (or all samples combined).</dd>
<dd class="advice">If this step failed, please review: the raw Kraken2 files available in your job results at &lt;sample_name&gt;/kraken_output/&lt;sample_name&gt;_k2_report.txt and &lt;sample_name&gt;/kraken_output/&lt;sample_name&gt;_k2_output.txt.</dd>

<dt>Sankey Plot</dt>
<dd>A flow diagram visualizing how reads are distributed across taxonomic ranks.</dd>
<dd class="advice">If this step failed, please review: the raw Kraken2 files available in your job results at &lt;sample_name&gt;/kraken_output/&lt;sample_name&gt;_k2_report.txt and &lt;sample_name&gt;/kraken_output/&lt;sample_name&gt;_k2_output.txt.</dd>

<dt>Multisample Comparison</dt>
<dd>A combined comparison table across all samples in the job (only runs when more than one sample is present).</dd>
<dd class="advice">If this step failed, please review: confirm that Kraken2 completed successfully for every sample in the job &mdash; this step depends on all of them.</dd>

<dt>MultiQC Report</dt>
<dd>Aggregates quality-control metrics from earlier steps into a single summary report.</dd>
<dd class="advice">
<h3>Using MultiQC for Troubleshooting</h3>
<p>If the classification results are unexpected or difficult to interpret, reviewing read quality is often the best first troubleshooting step. Problems such as poor sequencing quality, adapter contamination, low-complexity reads, or an unusually low number of reads can reduce the number of confidently classified sequences and may lead to poor abundance estimates.</p>
<p>The MultiQC report provides a convenient summary of these metrics across all samples. Pay particular attention to the following sections:</p>
<ul>
<li><strong>General Statistics</strong> &ndash; Verify the total number of reads and ensure the sample contains sufficient sequencing data.</li>
<li><strong>Per Base Sequence Quality</strong> &ndash; Look for consistently high quality scores (typically Phred scores above 30). A sharp decline in quality toward the ends of reads may indicate the need for additional trimming.</li>
<li><strong>Per Sequence Quality Scores</strong> &ndash; Check that most reads have high overall quality rather than a large fraction of low-quality reads.</li>
<li><strong>Per Base Sequence Content</strong> &ndash; Large deviations from the expected nucleotide composition can indicate library preparation issues or contamination.</li>
<li><strong>Adapter Content</strong> &ndash; High levels of adapter contamination suggest that trimming may be incomplete.</li>
<li><strong>Sequence Duplication Levels</strong> &ndash; Excessive duplication may indicate PCR bias or low library complexity.</li>
<li><strong>Overrepresented Sequences</strong> &ndash; Review any highly abundant sequences, which may represent contaminants, adapters, or other technical artifacts.</li>
</ul>
<p>While minor deviations are common, samples with low read counts, poor sequence quality, or substantial contamination may produce unreliable taxonomic classifications and abundance estimates. Reviewing the MultiQC report can help determine whether read quality, rather than the classification software, is responsible for unexpected results.</p>
</dd>

</dl>
</div>
</div>
</body>
</html>
"""
    with open(f"{output_dir}/JobStatus.html", "w") as job_status_file:
        job_status_file.write(html)


def post_processing_check(all_sample_ids, output_dir, input_dict):
    paired_sample_ids, _ = get_paired_and_single_sample_ids(input_dict)
    rows = build_job_status_rows(all_sample_ids, output_dir, input_dict, paired_sample_ids)
    write_job_status_html(rows, output_dir)

    complete = [f"{sample_name}: {step}" for sample_name, step, ok in rows if ok]
    incomplete = [f"{sample_name}: {step}" for sample_name, step, ok in rows if not ok]

    ## write the complete and incomplete steps to stderr for user to see ##
    if len(incomplete) != 0:
        msg = f"Reliable kraken results produced for the following samples: {complete}. \n \
                Review the following samples: **{incomplete}** \n"
        sys.stderr.write(msg)
    else:
        msg = f"Reliable kraken results produced for the following samples: {complete}. \n \
        Post Kraken processing complete"
        sys.stderr.write(msg)
        return True


def preprocessing_check(input_dir, output_dir, input_dict):
    ## Using the input dictionary instead of path from staging directory due to file name confusion
    paired_sample_ids, single_sample_ids = get_paired_and_single_sample_ids(input_dict)
    # Parse "sample_id" from "srr_libs"
    srr_sample_ids = [clean_sample_id(item["sample_id"]) for item in input_dict.get("srr_libs", [])]
    # Merge all sample IDs into one list
    all_sample_ids = paired_sample_ids + single_sample_ids + srr_sample_ids

    # file check for preprocessing/kraken run
    dict_samples = {}
    complete = []
    incomplete = []

    # check if Kraken report exists for each sample ID
    for sample_name in all_sample_ids:
        report_path = f"{output_dir}/{sample_name}/kraken_output/{sample_name}_k2_output.txt"
        if os.path.isfile(report_path) == True:
            complete.append(sample_name)
        else:
            dict_samples[sample_name] = False
            incomplete.append(sample_name)
    if len(incomplete) != 0:
        msg = f"Ending job. Not proceeding with the rest of the analysis due to errors in FASTQ proccessing \n \
                FASTQ Processing is complete for the following samples: {complete}. \n \
                FASTQ Processing is INCOMPLETE for the following samples: **{incomplete}** \n"
        sys.stderr.write(msg)
        ## return False
    else:
        msg = "preprocessing complete \n"
        sys.stderr.write(msg)
        return True, all_sample_ids

def run_16s_snakefile(input_dict, input_dir, output_dir,  config):
    # process any paired reads
    input_dir = config["input_data_dir"]
    SNAKEMAKE_PATH = config["snakemake"]
    SNAKEFILE_DIR = f"{config['workflow_dir']}/snakefile/"

    common_params = [
        SNAKEMAKE_PATH,
        "--cores", str(config["cores"]),
        "--use-singularity",
        "--verbose",
        "--printshellcmds",
        "--keep-going"
        ]
    if config["cores"] == 1:
        common_params.append("--debug")
        
    if os.path.exists(f"{input_dir}/pe_reads"):
        SNAKEFILE = os.path.join(SNAKEFILE_DIR, "pe_fastq_processing_16s")
        cmd = common_params + ["--snakefile",  SNAKEFILE]
        print(shlex.join(cmd), file=sys.stderr)
        subprocess.run(cmd)

    # Process any single end reads
    if os.path.exists(f"{input_dir}/se_reads"):
        SNAKEFILE = os.path.join(SNAKEFILE_DIR, "se_fastq_processing_16s")
        cmd = common_params + ["--snakefile",  SNAKEFILE]
        print(shlex.join(cmd), file=sys.stderr)
        subprocess.run(cmd)
    kraken_check_result = preprocessing_check(input_dir, output_dir, input_dict)
    kraken_check = kraken_check_result[0]
    all_sample_ids = kraken_check_result[1]
    
    if  kraken_check == True:
        # If all of the kraken report files exist then the final step will trigger.
        # Paired and single end reads will be processed together
        if input_dict["analysis_type"] == "16S":
            SNAKEFILE = os.path.join(SNAKEFILE_DIR, "16s_analysis")
            cmd = common_params + ["--snakefile",  SNAKEFILE]
            print(shlex.join(cmd), file=sys.stderr)
            subprocess.run(cmd)
    # File check if anything is wrong with the kraken outputs it should fail to 
    # Produce krona plot
    krona_check = post_processing_check(all_sample_ids, output_dir, input_dict)
    if krona_check == True:
        msg = "16s analysis is complete \n"
        sys.stderr.write(msg)
    return

# Run the snakefile command for WGS
def run_wgs_snakefile(input_dict, input_dir, output_dir,  config):

    input_dir = config["input_data_dir"]
    SNAKEMAKE_PATH = config["snakemake"]
    SNAKEFILE_DIR = f"{config['workflow_dir']}/snakefile/"

    common_params = [
        SNAKEMAKE_PATH,
        "--cores", str(config["cores"]),
        "--use-singularity",
        "--verbose",
        "--printshellcmds",
        "--keep-going"
        ]

    if config["cores"] == 1:
        common_params.append("--debug")
    
    # Process any paired reads
    if os.path.exists(f"{input_dir}/pe_reads"):
        SNAKEFILE = os.path.join(SNAKEFILE_DIR, "pe_fastq_processing")
        cmd = common_params + ["--snakefile",  SNAKEFILE]
        print(shlex.join(cmd), file=sys.stderr)
        subprocess.run(cmd)

    # Process any single reads
    if os.path.exists(f"{input_dir}/se_reads"):
        SNAKEFILE = os.path.join(SNAKEFILE_DIR, "se_fastq_processing")
        cmd = common_params + ["--snakefile",  SNAKEFILE]
        print(shlex.join(cmd), file=sys.stderr)
        subprocess.run(cmd)
    kraken_check_result = preprocessing_check(input_dir, output_dir, input_dict)
    ### check!
    kraken_check = kraken_check_result[0]
    all_sample_ids = kraken_check_result[1]

    if  kraken_check == True:
        # If all of the kraken report files exist then the final step will trigger.
        # Paired and single end reads will be processed together
        if input_dict["analysis_type"] == "pathogen":
            SNAKEFILE = os.path.join(SNAKEFILE_DIR, "pathogen_analysis")
            cmd = common_params + ["--snakefile",  SNAKEFILE]
            print(shlex.join(cmd), file=sys.stderr)
            subprocess.run(cmd)

        if input_dict["analysis_type"] == "microbiome":
            SNAKEFILE = os.path.join(SNAKEFILE_DIR, "microbiome_analysis")
            cmd = common_params + ["--snakefile",  SNAKEFILE]
            print(shlex.join(cmd), file=sys.stderr)
            subprocess.run(cmd)
    # Final file check
    krona_check = post_processing_check(all_sample_ids, output_dir, input_dict)
    if krona_check == True:
        msg = " WGS analysis is complete \n"
        sys.stderr.write(msg)
    return


def set_up_sample_dictionary(input_dir, input_dict, output_dir, cores):
    #### paired reads ####
    files_to_move = []
    gzip_targets = []

    if input_dict.get("paired_end_libs"):
        paired_sample_dict = {}
        pe_path= f"{input_dir}/pe_reads"
        os.makedirs(pe_path, exist_ok = True)

        # for item in ws_paired_reads:
        for item in input_dict["paired_end_libs"]:
            # Managing files
            sample_id_r1, file_extension_r1 = check_input_fastqs(item["read1"], item["sample_id"], gzip_targets)
            sample_id_r2, file_extension_r2 = check_input_fastqs(item["read2"], item["sample_id"], gzip_targets)
            # This will either use a zipped file or the file will already be zipped when the command runs
            files_to_move.append((item["read1"]+"{}".format(file_extension_r1), "{}/{}_R1.fastq.gz".format(pe_path, sample_id_r1)))
            files_to_move.append((item["read2"]+"{}".format(file_extension_r2), "{}/{}_R2.fastq.gz".format(pe_path, sample_id_r2)))
            
            # Logging input for sample key csv
            read1_filename  = item["read1"].split("/")[-1]
            paired_sample_dict[read1_filename] = "{}/{}_R1.fastq.gz".format(pe_path, sample_id_r1)
            read2_filename  = item["read2"].split("/")[-1]
            paired_sample_dict[read2_filename] = "{}/{}_R2.fastq.gz".format(pe_path, sample_id_r2)

    #### single reads ####
    if input_dict.get("single_end_libs"):
        single_end_sample_dict = {}
        se_path= f"{input_dir}/se_reads"
        os.makedirs(se_path, exist_ok = True)

        # for item in ws_single_end_reads:
        for item in input_dict["single_end_libs"]:
            # Managing files
            sample_id, file_extension = check_input_fastqs(item["read"], item["sample_id"], gzip_targets)
            # This will either use a zipped file or the file will already be zipped when the command runs
            files_to_move.append((item["read"]+"{}".format(file_extension), "{}/{}.fastq.gz".format(se_path, sample_id)))


            # Logging for sample key csv
            se_filename = item["read"].split("/")[-1]
            single_end_sample_dict[se_filename] = "{}/{}.fastq.gz".format(se_path, sample_id)

    # In parallel, gzip anything that needs to be gzipped
    # We deferred the move until after the zip
    #
    if gzip_targets:
        par_inp = "\n".join(gzip_targets + [""])
        subprocess.run(["parallel", "-j", str(cores), "gzip", "{}"], input=par_inp,text=True)
    

    for (src, dest) in files_to_move:
        # evaulate the paths for gz path
        # actual_src = src if src.endswith("gz") else src + ".gz"
        shutil.move(src, dest)

    # Export the sample dictionary to .CSV
    with open(f"{output_dir}/sample_key.csv", "w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(["User sample name", "Analysis sample name"])
        if len(input_dict["paired_end_libs"]) != 0:
            for key, value in paired_sample_dict.items():
                writer.writerow([key, value])
        if len(input_dict["single_end_libs"]) != 0:
            for key, value in single_end_sample_dict.items():
                writer.writerow([key, value])
    return


# run the script from service-script/app_taxonomic_classification perl script
# It takes a single argument which is the pathname of a config.json file
# This contains the app parameters in the params slot.
def main(argv):
    params_file = argv[0]
    try:
        fh = open(params_file, encoding="utf-8")
        config = json.load(fh)
        fh.close()
    except IOError:
        print(f"Could not read params file {params_file}")
        sys.exit(1)
    except json.JSONDecodeError as e:
        print(f"JSON parse error in pyfile {params_file}: {e}")
        sys.exit(1)

    input_dict = config["params"]
    input_dir = config["input_data_dir"]
    output_dir = config["output_data_dir"]
    hisat_indicies_path = config["hisat2_indicies_path"]
    set_up_sample_dictionary(input_dir, input_dict, output_dir, min(8, int(config["cores"])))
    try:
        load_hisat_indicies(input_dict, hisat_indicies_path, input_dir)
    except:
        pass
    if input_dict["analysis_type"] == "pathogen" or input_dict["analysis_type"] == "microbiome":
        run_wgs_snakefile(input_dict, input_dir, output_dir, config)
    elif input_dict["analysis_type"] == "16S":
        run_16s_snakefile(input_dict, input_dir, output_dir,  config)
    else:
        msg = f"unsuported analysis type {input_dict['analysis_type']} \n"
        sys.stderr.write(msg)
        sys.exit(1)

if __name__ == "__main__":
    main(sys.argv[1:])
