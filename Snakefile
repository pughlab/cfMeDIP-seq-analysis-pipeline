import yaml
import pandas as pd
import gzip
import statistics

configfile: "config.yml"

# ---------------------------- #
#  Important Shared Functions  #
# ---------------------------- #

fragment_length_windows = pd.read_csv('data/manual/fragment_length_windows.tsv', delimiter='\t')
fragment_length_windows_list = (fragment_length_windows.chr + '-' + fragment_length_windows.start.astype(str) + '-' + fragment_length_windows.end.astype(str)).tolist()

def get_active_cohorts():
    """Returns a list of active cohorts.
    Active cohorts are defined in config.yml data > cohorts > [cohortname] > active: True
    """
    return(c for c in config['data']['cohorts'] if config['data']['cohorts'][c]['active'])

def get_cohort_config(cohort):
    """Returns the configuration details for specific cohort.
    Importantly, retrieving a cohort using this function (rather than directly
    accessing it via the config dictionary) will also inject the default settings
    (config > data > defaults) wherever a specific setting is not specified.
    """
    config_data = dict(config['data']['defaults'])
    if 'settings' in config['data']['cohorts'][cohort]:
        for key in config['data']['cohorts'][cohort]['settings']:
            if key in config_data:
                config_data[key] = config['data']['cohorts'][cohort]['settings'][key]
            else:
                raise(Exception('Setting {} does not exist for cohort {}. Available settings: {}'.format(
                    key, cohort, ', '.join(config_data.keys())
                )))
    return(config_data)

def get_cohort_outdir(cohort):
    config_data = get_cohort_config(cohort)
    return(config_data['output_dir'])

def get_cohort_data(cohort):
    """Parses the samplesheet for a specific cohort.
    Also removes any excluded_cases from the samplesheet before it is returned.
    """
    cohort_config = get_cohort_config(cohort);
    samplesheet = pd.read_csv(config['data']['cohorts'][cohort]['samplesheet'], comment='#').drop_duplicates()
    if 'exclude' in cohort_config and cohort_config['exclude'] is not None:
        samplesheet = samplesheet[~samplesheet.sample_name.isin(cohort_config['exclude'])]
    return(samplesheet)

def get_fastq_path(cohort, sample, library, read_in_pair=1):
    """Retrieves the path to the fastq file.

    Keyword arguments:
        cohort -- name of the cohort whose samplesheet should be accessed.
        sample -- identifier of the sample as specified in the samplesheet.
        library -- integer representing the library index as specified in the samplesheet.
        read_in_pair -- 1 or 2 - representing read 1 or read 2 in paired end data.
    """
    library = int(library)
    cohort_data = get_cohort_data(cohort)
    sample_line = cohort_data[
        (cohort_data.sample_name == sample) &
        (cohort_data.library_index == library) &
        (cohort_data.read_in_pair == read_in_pair)
    ]
    return(sample_line.path.to_list()[0])

def get_all_samples(cohort=None):
    """Retrieves all samples to be processed.
    Does so by calling get_cohort_data, and therefore filters out excluded_cases.

    Keyword arguments:
        cohort -- Name of a cohort, OPTIONAL. If not specified, returns all samples
                  across all cohorts.
    """
    all_samples = pd.concat([
        get_cohort_data(cohort_name).assign(cohort_name = cohort_name)
        for cohort_name
        in config['data']['cohorts']
        if config['data']['cohorts'][cohort_name]['active']
        ])

    if cohort is None:
        return(all_samples)
    else:
        return(all_samples[all_samples.cohort_name == cohort])

def get_all_samples_list():
    """Returns all samples in list format.
    By calling get_all_samples()
    """
    return(get_all_samples().sample_name.unique().tolist())

def get_all_samples_with_cohorts():
    """Returns a list of tuples with cohort name and sample name."""
    samples = get_all_samples()[['cohort_name', 'sample_name']]
    return(zip(
        samples.drop_duplicates().cohort_name.tolist(),
        samples.drop_duplicates().sample_name.tolist()
    ))

def clean(command):
    """Cleans a snakemake command string by replacing whitespace with a single space.
    Useful for stripping down multiline shell commands so that they read more
    cleanly on snakemake's printed output.
    """
    command = command.replace('\n', ' ').replace('\t', ' ')
    while '  ' in command:
        command = command.replace('  ', ' ')
    return(command.strip())

def get_cohort_signature_pairs(bed=False):
    """Iterator returning every signature for every active cohort."""
    for active_cohort in get_active_cohorts():
        for signature in get_cohort_config(active_cohort)['signatures_bed' if bed else 'signatures']:
            yield (get_cohort_outdir(active_cohort), active_cohort, signature)

# ------------------------------ #
#  Beginning of Snakemake Rules  #
# ------------------------------ #

rule all:
    input:
        #['{output_dir}/pipeline/{cohort}/results/qc/{sample}_qc_complete'.format(
        #    output_dir = get_cohort_outdir(cohort),
        #    cohort = cohort,
        #    sample = sample
        #    ) for cohort,sample in get_all_samples_with_cohorts()
        #],
        #['{output_dir}/pipeline/{cohort}/results/aggregated/all_insertsize.tsv'.format(
        #    output_dir = get_cohort_outdir(cohort),
        #    cohort = cohort
        #    ) for cohort in get_active_cohorts()
        #],
        #[
        #    '{}/{}/results/aggregated/all_bedpe/{}.tsv'.format(output_dir, cohort, signature)
        #    for output_dir, cohort, signature in get_cohort_signature_pairs(bed=True)
        #],
        [
            '{output_dir}/{cohort}/results/aggregated/fragment_length_windows.feather'.format(
                output_dir = get_cohort_outdir(cohort),
                cohort = cohort
            ) for cohort in get_active_cohorts()
        ],
        #expand(
        #    [
        #    '{}/{}/results/extract_signature_matrix/{}-{}-{{valuetype}}.parquet'.format(output_dir, cohort, cohort, signature)
        #        for output_dir, cohort, signature in get_cohort_signature_pairs()
        #    ],
        #    valuetype = ['coverage', 'posterior']
        #    #valuetype = ['coverage', 'posterior', 'medestrand']
        #),

# ------------------------------------------ #
#  Install conda dependencies automatically  #
# ------------------------------------------ #
# This section is only necessary to support air-gapped environments,
# where the cluster nodes have no direct access to internet.
# In this setting, you must first load all necessary packages from
# a build node with internet access. This can be done simply
# by running 
# bash update_conda.sh
# from a node with internet access, which loads install_dependencies
# rule below to install all of the necessary conda dependencies
# in advance. Any time the git repo is pulled or if you modify
# any of the conda environments in conda_env, this needs to be
# re-run.

rule install_dependencies:
    input:
        'conda_env/biopython_installed',
        'conda_env/cfmedip_r_installed',
        'conda_env/fastqc_installed',
        'conda_env/qualimap_installed',
        'conda_env/samtools_installed',
        'conda_env/umitools_installed',
        'conda_env/trimgalore_installed'
    output:
        'conda_env/dependencies'
    shell:
        'touch {output}'

rule biopython:
    output:
        'conda_env/biopython_installed'
    conda: 'conda_env/biopython.yml'
    shell:
        'touch {output}'

rule cfmedip_r:
    output:
        'conda_env/cfmedip_r_installed'
    conda: 'conda_env/cfmedip_r.yml'
    shell:
        'touch {output}'

rule fastqc:
    output:
        'conda_env/fastqc_installed'
    conda: 'conda_env/fastqc.yml'
    shell:
        'touch {output}'

rule qualimap:
    output:
        'conda_env/qualimap_installed'
    conda: 'conda_env/qualimap.yml'
    shell:
        'touch {output}'

rule samtools:
    output:
        'conda_env/samtools_installed'
    conda: 'conda_env/samtools.yml'
    shell:
        'touch {output}'

rule umitools:
    output:
        'conda_env/umitools_installed'
    conda: 'conda_env/umitools.yml'
    shell:
        'touch {output}'

rule trimgalore:
    output:
        'conda_env/trimgalore_installed'
    conda: 'conda_env/trimgalore.yml'
    shell:
        'touch {output}'

# -------------------------- #
#  Pre-process input FASTQs  #
# -------------------------- #

rule gunzip_fastq:
    input:
        lambda wildcards: get_fastq_path(wildcards.cohort, wildcards.sample, int(wildcards.lib), int(wildcards.read)),
    output:
        temp('{output_dir}/pipeline/{cohort}/tmp/gunzip_fastq/{sample}_lib{lib}_R{read}.fastq')
    resources: cpus=1, mem_mb=8000, time_min='24:00:00'
    shell:
        'gunzip -dc {input} > {output}'

# QC of input FASTQs using FASTQC
rule fastqc_fastq:
    input:
        '{output_dir}/pipeline/{cohort}/tmp/gunzip_fastq/{sample}_lib{lib}_R{read}.fastq'
    output:
        html='{output_dir}/pipeline/{cohort}/results/qc_input/{sample}_lib{lib}_R{read}_fastqc.html',
        zipfile='{output_dir}/pipeline/{cohort}/results/qc_input/{sample}_lib{lib}_R{read}_fastqc.zip'
    resources: cpus=1, mem_mb=8000, time_min='24:00:00'
    params:
        outdir = lambda wildcards, output: '/'.join(output.html.split('/')[0:-1])
    conda: 'conda_env/fastqc.yml'
    shell:
        'fastqc --outdir {params.outdir} {input}'


# ---------------------------------- #
#  Extract Barcodes using UMI Tools  #
# ---------------------------------- #
#
# This has superseded the use of ConsensusCruncher.
# UMI barcodes are removed and placed into FASTQ headers.
rule umi_tools_extract:
    input:
        R1_qc = '{output_dir}/pipeline/{cohort}/results/qc_input/{sample}_lib{lib}_R1_fastqc.html',
        R2_qc = '{output_dir}/pipeline/{cohort}/results/qc_input/{sample}_lib{lib}_R2_fastqc.html',
        R1 = '{output_dir}/pipeline/{cohort}/tmp/gunzip_fastq/{sample}_lib{lib}_R1.fastq',
        R2 = '{output_dir}/pipeline/{cohort}/tmp/gunzip_fastq/{sample}_lib{lib}_R2.fastq'
    output:
        R1 = temp('{output_dir}/pipeline/{cohort}/tmp/umi_tools_extract/{sample}_lib{lib}_umitools_R1.fastq'),
        R2 = temp('{output_dir}/pipeline/{cohort}/tmp/umi_tools_extract/{sample}_lib{lib}_umitools_R2.fastq'),
        logfile = '{output_dir}/pipeline/{cohort}/results/umi_tools_extract/{sample}_lib{lib}_umitools_log.txt'
    params:
        barcodes = lambda wildcards: get_cohort_config(wildcards.cohort)['barcodes']
    resources: cpus=1, mem_mb=30000, time_min='5-00:00:00'
    conda: 'conda_env/umitools.yml'
    shell:
        'umi_tools extract '
        '--extract-method=regex '
        '--stdin={input.R1} '
        '--read2-in={input.R2} '
        '--bc-pattern="{params.barcodes}" '
        '--bc-pattern2="{params.barcodes}" '
        '--stdout={output.R1} '
        '--read2-out={output.R2} '
        '--log={output.logfile} '

# Extract Barcodes using ConsensusCruncher
# Pulls the path to extract_barcodes.py from config > paths > dependencies > extract_barcodes_path
#rule extract_barcodes:
#    input:
#        R1_qc = '{output_dir}/pipeline/{cohort}/results/qc_input/{sample}_lib{lib}_R1_fastqc.html',
#        R2_qc = '{output_dir}/pipeline/{cohort}/results/qc_input/{sample}_lib{lib}_R2_fastqc.html',
#        R1 = '{output_dir}/pipeline/{cohort}/tmp/gunzip_fastq/{sample}_lib{lib}_R1.fastq',
#        R2 = '{output_dir}/pipeline/{cohort}/tmp/gunzip_fastq/{sample}_lib{lib}_R2.fastq'
#    output:
#        R1 = temp('{output_dir}/pipeline/{cohort}/tmp/extract_barcodes/{sample}_lib{lib}_extract_barcode_R1.fastq'),
#        R2 = temp('{output_dir}/pipeline/{cohort}/tmp/extract_barcodes/{sample}_lib{lib}_extract_barcode_R2.fastq')
#    params:
#        outprefix = lambda wildcards, output: output.R1.split('_barcode_')[0],
#        barcodes = lambda wildcards: get_cohort_config(wildcards.cohort)['barcodes']
#    resources: cpus=1, mem_mb=30000, time_min='5-00:00:00'
#    conda: 'conda_env/biopython.yml'
#    shell:
#        clean(r'''
#        python {extract_barcodes}
#            --read1 {{input.R1}}
#            --read2 {{input.R2}}
#            --outfile {{params.outprefix}}
#            {{params.barcodes}}
#        '''.format(extract_barcodes = config['paths']['dependencies']['extract_barcodes_path']))

# Trims FASTQ using trimgalore to remove barcode sequences
# By default, trims 10 base pairs from the 5' end, which seems to be correct for OICR cfMeDIP-seq output.
# This can be configured in the config.yml under data > cohorts > settings > trimgalore.
rule trim_fastq:
    input:
        R1 = '{output_dir}/pipeline/{cohort}/tmp/umi_tools_extract/{sample}_lib{lib}_umitools_R1.fastq',
        R2 = '{output_dir}/pipeline/{cohort}/tmp/umi_tools_extract/{sample}_lib{lib}_umitools_R2.fastq',
    output:
        trimmed_1 = temp('{output_dir}/pipeline/{cohort}/tmp/trim_fastq/{sample}_lib{lib}_umitools_R1_val_1.fq'),
        trimmed_2 = temp('{output_dir}/pipeline/{cohort}/tmp/trim_fastq/{sample}_lib{lib}_umitools_R2_val_2.fq'),
        report_1 = '{output_dir}/pipeline/{cohort}/results/qc/{sample}_lib{lib}_umitools_R1.fastq_trimming_report.txt',
        report_2 = '{output_dir}/pipeline/{cohort}/results/qc/{sample}_lib{lib}_umitools_R2.fastq_trimming_report.txt'
    params:
        outdir = lambda wildcards, output: '/'.join(output.trimmed_1.split('/')[0:-1]),
        trimgalore_settings = lambda wildcards: get_cohort_config(wildcards.cohort)['trimgalore']
    resources: cpus=4, mem_mb=8000, time_min='24:00:00'
    conda: 'conda_env/trimgalore.yml'
    shell:
        'trim_galore --cores 4 --dont_gzip --paired '
        '{params.trimgalore_settings} '
        '--output_dir {params.outdir} '
        '{input.R1} {input.R2} '
        '&& mv {wildcards.output_dir}/pipeline/{wildcards.cohort}/tmp/trim_fastq/{wildcards.sample}_lib{wildcards.lib}_umitools_R*.fastq_trimming_report.txt {wildcards.output_dir}/pipeline/{wildcards.cohort}/results/qc '

# --------------------------- #
#  Align FASTQs to reference  #
# --------------------------- #

def get_read_group_from_fastq(fastq_file, sample_name):
    """Extracts the read group data from FASTQs automatically.
    This information is used in the shell script of rule bwa_mem.
    """
    with gzip.open(fastq_file, 'rt') as fastq:
        header = next(fastq)
        (instrument, run_number, flowcell, lane, tile, xpos, ypos) = header.split(' ')[0].split(':')
        lib_value = header.strip().split(' ')[1].split(':')[3]
        rg_line = r"@RG\tID:{flowcell}_{lane}\tSM:{sample}\tPL:Illumina\tPU:.\tLB:{lib_value}".format(
            flowcell = flowcell,
            lane = lane,
            sample = sample_name,
            lib_value = lib_value
        )
        return(rg_line)

# Run BWA mem on FASTQs after extracting barcodes.
rule bwa_mem:
    input:
        '{output_dir}/pipeline/{cohort}/tmp/trim_fastq/{sample}_lib{lib}_umitools_R1_val_1.fq',
        '{output_dir}/pipeline/{cohort}/tmp/trim_fastq/{sample}_lib{lib}_umitools_R2_val_2.fq',
    output:
        temp('{output_dir}/pipeline/{cohort}/tmp/bwa_mem/{sample}_lib{lib}.sam')
    resources: cpus=4, mem_mb=16000, time_min='72:00:00'
    params:
        read_group = lambda wildcards, input: get_read_group_from_fastq(
            fastq_file = get_fastq_path(wildcards.cohort, wildcards.sample, wildcards.lib),
            sample_name = wildcards.sample
        ),
        bwa_index = lambda wildcards: get_cohort_config(wildcards.cohort)['bwa_index']
    conda: 'conda_env/samtools.yml'
    shell:
        clean(r"""
        bwa mem -M -t4
        -R'{params.read_group}' 
        {params.bwa_index}
        {input} > {output}""")

# Converts SAM to BAM with sorting
rule sam_to_sorted_bam:
    input:
        '{output_dir}/pipeline/{cohort}/tmp/bwa_mem/{sample}_lib{lib}.sam'
    output:
        bam = temp('{output_dir}/pipeline/{cohort}/tmp/bwa_mem/{sample}_lib{lib}.sorted.bam'),
        index = temp('{output_dir}/pipeline/{cohort}/tmp/bwa_mem/{sample}_lib{lib}.sorted.bam.bai'),
    resources: cpus=32, mem_mb=30000, time_min='72:00:00'
    conda: 'conda_env/samtools.yml'
    shell:
        # Try it without fixmate for replication purposes
        #"samtools view -buS -f 2 -F 4 -@4 {input} | samtools fixmate -m - - | samtools sort -@4 -o {output.bam} && samtools index {output.bam}"
        'samtools view -buS -f 2 -F 4 -@32 {input} | '
        'samtools fixmate -m - - | '
        'samtools sort -@32 -o {output.bam} && samtools index {output.bam} '

def get_libraries_of_sample(sample):
    """Returns all library indices of a sample based on samplesheet."""
    filtered_table = get_all_samples()[get_all_samples().sample_name == sample]
    return(list(set(filtered_table.library_index.to_list())))

# If there are multiple libraries for a given sample, as specified in samplesheet,
# these libraries are automatically merged at this step into a single unified BAM.
rule merge_bam:
    input:
        lambda wildcards: expand(
                '{{output_dir}}/pipeline/{{cohort}}/tmp/bwa_mem/' + wildcards.sample + '_lib{lib}.sorted.bam',
                lib=get_libraries_of_sample(wildcards.sample)
        )
    output:
        temp('{output_dir}/pipeline/{cohort}/tmp/merge_bam/{sample}.aligned.sorted.bam'),
    resources: cpus=1, mem_mb=8000, time_min='24:00:00'
    conda: 'conda_env/samtools.yml'
    shell:
        'samtools merge {output} {input} '

rule index_merged_bam:
    input:
        '{output_dir}/pipeline/{cohort}/tmp/merge_bam/{sample}.aligned.sorted.bam',
    output:
        temp('{output_dir}/pipeline/{cohort}/tmp/merge_bam/{sample}.aligned.sorted.bam.bai')
    resources: cpus=1, mem_mb=8000, time_min='2:00:00'
    threads: 12
    conda: 'conda_env/samtools.yml'
    shell:
        'samtools index -@ {threads} {input} '

# We next deduplicate with UMItools.
rule bam_umitools_dedup:
    input:
        bam='{output_dir}/pipeline/{cohort}/tmp/merge_bam/{sample}.aligned.sorted.bam',
        index='{output_dir}/pipeline/{cohort}/tmp/merge_bam/{sample}.aligned.sorted.bam.bai'
    output:
        bam = temp('{output_dir}/pipeline/{cohort}/results/bam/{sample}_tmp.bam'),
    params:
        stat_prefix = '{output_dir}/pipeline/{cohort}/results/bam/{sample}.aligned.sorted.dedup',
    conda: 'conda_env/umitools.yml'
    resources: cpus=1, mem_mb=30000, time_min='5-00:00:00'
    shell:
        'umi_tools dedup --paired -I {input.bam} -S {output.bam} --umi-separator="_" --output-stats={params.stat_prefix}'

# We apply some filters to the reads, yielding final BAM file
rule finalize_bam:
    input:
        '{output_dir}/pipeline/{cohort}/results/bam/{sample}_tmp.bam',
    output:
        bam = '{output_dir}/pipeline/{cohort}/results/bam/{sample}.aligned.sorted.dedup.bam',
    conda: 'conda_env/samtools.yml'
    shell:
        'samtools view -b -f 2 -F 2828 --threads {threads} {input} > {output.bam} && '

rule index_final_bam:
    input:
        '{output_dir}/pipeline/{cohort}/results/bam/{sample}.aligned.sorted.dedup.bam',
    output:
        '{output_dir}/pipeline/{cohort}/results/bam/{sample}.aligned.sorted.dedup.bam.bai',
    resources: cpus=1, mem_mb=8000, time_min='2:00:00'
    threads: 12
    conda: 'conda_env/samtools.yml'
    shell:
        'samtools index -@ {threads} {input} '

# Extract fragment length data
# in 5 Mb windows
rule fragment_length_windows:
    input:
        bam = '{output_dir}/pipeline/{cohort}/results/bam/{sample}.aligned.sorted.dedup.bam',
        index = '{output_dir}/pipeline/{cohort}/results/bam/{sample}.aligned.sorted.dedup.bam.bai',
        windows = "data/manual/fragment_length_windows.tsv"
    output:
        '{output_dir}/pipeline/{cohort}/tmp/fragment_length_windows/{sample}.tsv',
    resources: cpus=1, mem_mb=8000, time_min='4:00:00'
    conda: 'conda_env/samtools.yml'
    shell:
        "python src/python/fragment_lengths.py "
        "-b {input.bam} "
        "-w {input.windows} "
        "-o {output} "

rule fragment_length_windows_paths:
    input:
        lambda wildcards: expand(
            '{{output_dir}}/pipeline/{{cohort}}/tmp/fragment_length_windows/{sample}.tsv',
            sample = get_cohort_data(wildcards.cohort).sample_name.unique().tolist(),
        )
    output:
        '{output_dir}/pipeline/{cohort}/tmp/fragment_length_windows_paths/paths.txt',
    resources: cpus=1, mem_mb=8000, time_min='1:00:00'
    run:
        with open(output[0], 'w') as out:
            out.write('\n'.join(input))

rule fragment_length_windows_aggregate:
    input:
        '{output_dir}/pipeline/{cohort}/tmp/fragment_length_windows_paths/paths.txt',
    output:
        '{output_dir}/pipeline/{cohort}/results/aggregated/fragment_length_windows.feather',
    resources: cpus=1, mem_mb=8000, time_min='1:00:00'
    conda: 'conda_env/cfmedip_r.yml'
    shell:
        'Rscript src/R/row_bind_tables.R -p {input} -o {output} --in-tsv --out-feather'

# Bam markdup and create index. 
# Now superseded by UMItools dedup above
# This step finalizes the definitive BAM file.
#rule bam_markdup:
#    input:
#        '{output_dir}/pipeline/{cohort}/tmp/merge_bam/{sample}.aligned.sorted.bam'
#    output:
#        bam = '{output_dir}/pipeline/{cohort}/results/bam_markdup/{sample}.aligned.sorted.markdup.bam',
#        index = '{output_dir}/pipeline/{cohort}/results/bam_markdup/{sample}.aligned.sorted.markdup.bam.bai'
#    resources: cpus=1, mem_mb=8000, time_min='24:00:00'
#    conda: 'conda_env/samtools.yml'
#    shell:
#        "samtools markdup -r {input} {output.bam} && samtools index {output.bam}"

# ----------------- #
#  QC of BAM files  #
# ----------------- #

# Run FASTQC on final BAM files.
rule fastqc_bam:
    input:
        bam='{output_dir}/pipeline/{cohort}/results/bam/{sample}.aligned.sorted.dedup.bam',
        index = '{output_dir}/pipeline/{cohort}/results/bam/{sample}.aligned.sorted.dedup.bam.bai',
    output:
        html = '{output_dir}/pipeline/{cohort}/results/qc/{sample}.aligned.sorted.dedup_fastqc.html',
        zipfile = '{output_dir}/pipeline/{cohort}/results/qc/{sample}.aligned.sorted.dedup_fastqc.zip',
    resources: cpus=1, mem_mb=8000, time_min='24:00:00'
    params:
        outdir = lambda wildcards, output: '/'.join(output.html.split('/')[0:-1])
    conda: 'conda_env/fastqc.yml'
    shell:
        'fastqc --outdir {params.outdir} {input.bam}'

# Run Qualimap QC metrics on final BAM files.
rule qualimap_bam:
    input:
        '{output_dir}/pipeline/{cohort}/results/bam/{sample}.aligned.sorted.dedup.bam',
        '{output_dir}/pipeline/{cohort}/results/bam/{sample}.aligned.sorted.dedup.bam.bai',
    output:
        html = '{output_dir}/pipeline/{cohort}/results/qc/qualimap_{sample}/qualimapReport.html',
    params:
        outdir = lambda wildcards, output: '/'.join(output.html.split('/')[0:-1])
    conda: 'conda_env/qualimap.yml'
    resources: cpus=1, mem_mb=10000, time_min='24:00:00'
    shell:
        'qualimap bamqc -bam {input[0]} -outdir {params.outdir} --java-mem-size=8G'

# Run Flagstat to get basic stats on final BAM files.
rule bam_flagstat:
    input:
        '{output_dir}/pipeline/{cohort}/results/bam/{sample}.aligned.sorted.dedup.bam',
        '{output_dir}/pipeline/{cohort}/results/bam/{sample}.aligned.sorted.dedup.bam.bai',
    output:
        '{output_dir}/pipeline/{cohort}/results/qc/{sample}.aligned.sorted.dedup.bam.flagstat',
    resources: cpus=1, mem_mb=8000, time_min='1:00:00'
    conda: 'conda_env/samtools.yml'
    shell:
        'samtools flagstat {input[0]} > {output}'

# Unified QC rule which runs all of the above QCs
rule bam_qc:
    input:
        fastqc = '{output_dir}/pipeline/{cohort}/results/qc/{sample}.aligned.sorted.dedup_fastqc.html',
        flagstat = '{output_dir}/pipeline/{cohort}/results/qc/{sample}.aligned.sorted.dedup.bam.flagstat',
        qualimap = '{output_dir}/pipeline/{cohort}/results/qc/qualimap_{sample}/qualimapReport.html',
    output:
        '{output_dir}/pipeline/{cohort}/results/qc/{sample}_qc_complete',
    resources: cpus=1, mem_mb=1000, time_min='00:01:00'
    conda: 'conda_env/samtools.yml'
    shell:
        'touch {output}'

# ------------------------------------------ #
#  Compute BAM bin stats including coverage  #
# ------------------------------------------ #

# Generate bin stats for each chromosome individually.
rule bam_bin_stats:
    input:
        '{output_dir}/pipeline/{cohort}/results/bam/{sample}.aligned.sorted.dedup.bam',
        '{output_dir}/pipeline/{cohort}/results/bam/{sample}.aligned.sorted.dedup.bam.bai',
    output:
        binstat = temp('{output_dir}/pipeline/{cohort}/tmp/bam_bin_stats/bin_stats_{sample}_{species}_{chrom}.tsv'),
        filtered = '{output_dir}/pipeline/{cohort}/results/bam_bin_stats/removed_bins_{sample}_{species}_{chrom}.tsv',
    params:
        bsgenome = lambda wildcards: get_cohort_config(wildcards.cohort)['bsgenome'][wildcards.species],
    resources: cpus=1, mem_mb=30000, time_min='24:00:00'
    conda: 'conda_env/cfmedip_r.yml'
    shell:
        clean('''
        Rscript src/R/bin_stats.R
            -b {input[0]}
            -g {params.bsgenome}
            -c {wildcards.chrom}
            -o {output.binstat}
            --filtered {output.filtered}
        ''')

# Preload chromosome values for below rule
chromosomes = {}
for c in config['data']['cohorts']:
    if config['data']['cohorts'][c]['active']:
        chr_data = get_cohort_config(c)['chromosomes']
        chromosome_tuples = [(species, chrom) for species in chr_data for chrom in chr_data[species].split(',')]
        chromosomes[c] = chromosome_tuples

# Merge the bin stats across all chromosomes.
rule merge_bin_stats:
    input:
        lambda wildcards: ['{{output_dir}}/pipeline/{{cohort}}/tmp/bam_bin_stats/bin_stats_{{sample}}_{species}_{chrom}.tsv'.format(species=a[0], chrom=a[1]) for a in chromosomes[wildcards.cohort]],
    output:
        '{output_dir}/pipeline/{cohort}/results/merge_bin_stats/bin_stats_{sample}.feather'
    params:
        paths = lambda wildcards, input: ','.join(input)
    resources: cpus=1, mem_mb=8000, time_min='24:00:00'
    conda: 'conda_env/cfmedip_r.yml'
    shell:
        'Rscript src/R/row_bind_tables.R -p "{params.paths}" -o {output} --in-tsv --out-feather --omit-paths'

# Converts bam to bedpe file - useful for fragmentomics work
rule bam_to_chr_signature_bedpe:
    input:
        signature = lambda wildcards: get_cohort_config(wildcards.cohort)['signatures_bed'][wildcards.signature],
        bam = '{output_dir}/pipeline/{cohort}/results/bam/{sample}.aligned.sorted.dedup.bam',
        index = '{output_dir}/pipeline/{cohort}/results/bam/{sample}.aligned.sorted.dedup.bam.bai',
    output: 
        temp('{output_dir}/pipeline/{cohort}/tmp/bam_to_chr_signature_bedpe/{sample}/{signature}/{chromosome}.bedpe'),
    resources: cpus=1, mem_mb=16000, time_min='24:00:00'
    conda: 'conda_env/samtools.yml'
    shell:
        'samtools view -b -F 0x900 {input.bam} {wildcards.chromosome} | samtools sort -T /tmp/ -n - | bedtools pairtobed -abam stdin -b {input.signature} -bedpe > {output}'

rule signature_bedpe:
    input:
        lambda wildcards: ['{{output_dir}}/pipeline/{{cohort}}/tmp/bam_to_chr_signature_bedpe/{{sample}}/{{signature}}/{chromosome}.bedpe'.format(chromosome=a[1]) for a in chromosomes[wildcards.cohort]],
    output:
        '{output_dir}/pipeline/{cohort}/results/signature_bedpe/{signature}/{sample}_{signature}.bedpe',
    resources: cpus=1, mem_mb=4000, time_min='24:00:00'
    shell:
        'cat {input} > {output}'

rule insert_size_distribution:
    input:
        bam='{output_dir}/pipeline/{cohort}/results/bam/{sample}.aligned.sorted.dedup.bam',
        index='{output_dir}/pipeline/{cohort}/results/bam/{sample}.aligned.sorted.dedup.bam.bai',
    output:
        txt = '{output_dir}/pipeline/{cohort}/results/insert_size_distribution/{sample}.picard.insertsize.txt',
        pdf = '{output_dir}/pipeline/{cohort}/results/insert_size_distribution/{sample}.picard.insertsize.pdf',
    conda: 'conda_env/samtools.yml'
    resources: cpus=1, mem_mb=8000, time_min='24:00:00'
    shell:
        'picard CollectInsertSizeMetrics I={input.bam} O={output.txt} H={output.pdf}'

rule all_bedpe:
    input:
        lambda wildcards: expand(
            '{{output_dir}}/pipeline/{{cohort}}/results/signature_bedpe/{{signature}}/{sample}_{{signature}}.bedpe',
            sample = get_cohort_data(wildcards.cohort).sample_name.unique().tolist()
        ),
        '{output_dir}/pipeline/{cohort}/results/aggregated/all_insertsize.tsv',
    output:
        '{output_dir}/pipeline/{cohort}/results/aggregated/all_bedpe/{signature}.tsv',
    run:
        with open(output[0], 'w') as outfile:
            for input_path in input:
                for line in open(input_path):
                    outfile.write('\t'.join(line.strip().split('\t') + [input_path.split('/')[-1]]) + '\n')

rule all_insertsize:
    input:
        lambda wildcards: expand(
            '{{output_dir}}/pipeline/{{cohort}}/results/insert_size_distribution/{sample}.picard.insertsize.txt',
            sample = get_cohort_data(wildcards.cohort).sample_name.unique().tolist()
        )
    output:
        '{output_dir}/pipeline/{cohort}/results/aggregated/all_insertsize.tsv',
    run:
        with open(output[0], 'w') as outfile:
            first_file = True
            for path in input:
                infile = open(path)
                for line in infile:
                    if line.startswith('## HISTOGRAM'):
                        break
                if first_file:
                    outfile.write('file\t' + next(infile))
                    first_file = False
                else:
                    next(infile)
                for line in infile:
                    outfile.write(path.split('/')[-1] + '\t' + line)

# ----------------------------------------------------------------------------- #
#  Fit cfMeDIP-seq coverage stats to infer absolute methylation using MedReMix  #
# ----------------------------------------------------------------------------- #

# This is the currently used negative binomial GLM approach to fitting
rule cfmedip_nbglm:
    input:
        '{output_dir}/pipeline/{cohort}/results/merge_bin_stats/bin_stats_{sample}.feather'
    output:
        fit='{output_dir}/pipeline/{cohort}/results/cfmedip_nbglm/{sample}_fit_nbglm.feather',
        model='{output_dir}/pipeline/{cohort}/results/cfmedip_nbglm/{sample}_fit_nbglm_model.Rds',
    resources: cpus=1, time_min='5-00:00:00', mem_mb=lambda wildcards, attempt: 16000 if attempt == 1 else 30000
    conda: 'conda_env/cfmedip_r.yml'
    shell:
        'Rscript src/R/cfmedip_nbglm.R -i {input} -o {output.fit} --modelout {output.model}'

# Below is the full bayesian approach for fitting, which remains under development
rule fit_bin_stats:
    input:
        '{output_dir}/pipeline/{cohort}/results/merge_bin_stats/bin_stats_{sample}.feather'
    output:
        '{output_dir}/pipeline/{cohort}/results/fit_bin_stats/{sample}_fit_{method}.tsv'
    resources: cpus=1, mem_mb=30000, time_min='5-00:00:00'
    conda: 'conda_env/cfmedip_r.yml'
    shell:
        'Rscript src/R/fit_cpg_bins.R -i {input} -o {output} --method {wildcards.method}'

# --------------------------------------------- #
#  MeDEStrand Method of methylation correction  #
# --------------------------------------------- #

rule medestrand:
    input:
        '{output_dir}/pipeline/{cohort}/results/bam/{sample}.aligned.sorted.dedup.bam',
    output:
        '{output_dir}/pipeline/{cohort}/results/medestrand/{sample}_medestrand_output.feather',
    resources: cpus=1, mem_mb=30000, time_min='5-00:00:00'
    conda: 'conda_env/cfmedip_r.yml'
    params: bsgenome = lambda wildcards: get_cohort_config(wildcards.cohort)['bsgenome']['human']
    shell:
        'Rscript src/R/run_medestrand.R -b {{input}} -o {{output}} -m {} --bsgenome {{params.bsgenome}}'.format(
            config['paths']['dependencies']['medestrand_path']
        )

# ----------------------------------------------------------- #
#  Assemble results overlapping specific regions of interest  #
# ----------------------------------------------------------- #
# Signatures define regions of interest from which we
# would like to capture quantitative methylation information.
#
# The results are similar to those from the matrices/binsets,
# except that signatures don't have to match bins in the data,
# and the scripts here will automatically find all bins
# overlapping the signature regions and output data from them.

# Each type of matrix is composed from specific files
# and specific column within each file category.
# For example, a matrix of posteriors pulls data from
# column 8 of the output files of cfmedip_nbglm.

data_sources = {
    'coverage': {
        'path': '{output_dir}/pipeline/{cohort}/results/cfmedip_nbglm/{sample}_fit_nbglm.feather',
        'filetype': 'feather',
        'colname': 'coverage',
        'clip': '_fit_nbglm.feather'
    },
    'posterior': {
        'path': '{output_dir}/pipeline/{cohort}/results/cfmedip_nbglm/{sample}_fit_nbglm.feather',
        'filetype': 'feather',
        'colname': 'methylated_posterior',
        'clip': '_fit_nbglm.feather'
    },
    'medestrand': {
        'path': '{output_dir}/pipeline/{cohort}/results/medestrand/{sample}_medestrand_output.feather',
        'filetype': 'feather',
        'colname': 'binMethyl',
        'clip': '_medestrand_output.feather'
    }
}

def get_matrix_inputs(wildcards):
    """Get paths of all files used to assemble a specific matrix.
    """
    samples = get_cohort_data(wildcards.cohort).sample_name.unique().tolist()
    paths = expand(data_sources[wildcards.valuetype]['path'], output_dir = wildcards.output_dir, cohort = wildcards.cohort, sample = samples)
    return(paths)

# Assembles paths of files for assembling a matrix and
# write them to a file, one path per line
rule signature_matrix_paths:
    input:
        get_matrix_inputs
    output:
        temp('{output_dir}/pipeline/{cohort}/tmp/signature_matrix_paths/matrix.{valuetype}.paths.txt')
    run:
        with open(output[0], 'w') as out:
            out.write('\n'.join(input))

# Using the files assembled in signature_matrix_paths,
# and signature regions specified in a file pointed to in
# config > data > cohorts > signatures or the defaults
# specified in config > data > defaults > signatures.
rule extract_signature_matrix:
    input:
        signature = lambda wildcards: get_cohort_config(wildcards.cohort)['signatures'][wildcards.signature],
        paths = '{output_dir}/pipeline/{cohort}/tmp/signature_matrix_paths/matrix.{valuetype}.paths.txt'
    output:
        '{output_dir}/pipeline/{cohort}/results/extract_signature_matrix/{cohort}-{signature}-{valuetype}.parquet'
    resources: cpus=1, mem_mb=30000, time_min='5-00:00:00'
    conda: 'conda_env/cfmedip_r.yml'
    params:
        colname = lambda wildcards: data_sources[wildcards.valuetype]['colname'],
        filetype = lambda wildcards: data_sources[wildcards.valuetype]['filetype'],
        clip = lambda wildcards: data_sources[wildcards.valuetype]['clip']
    shell:
        clean("""
        Rscript src/R/extract_signature_matrix.R 
            -s {input.signature} 
            -p {input.paths} 
            -c {params.colname} 
            -f {params.filetype}
            -o {output} 
            --clip '{params.clip}'
        """)

# -------------- #
#  Legacy rules  #
# -------------- #

#include: 'legacy.smk'
# This line includes some rules from legacy code, which are
# gradually being retired but still accessible if needed.


