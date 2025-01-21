import glob

# set project directory
work_dir = 'path/to/data'

# get sample list
samples = [s.split('/')[-1].split('_')[0] for s in glob.glob(work_dir+'/reads/raw/*_R1.fastq.gz')]
print("samples:", samples)

# HiSat2 human database
hisat2_index = 'path/to/directory/with/index'

# run preprocessing pipeline
rule run_preproc:
    input: expand("{work_dir}/reads/raw__filtered/{sample}_R1.fastq.gz", sample = samples, work_dir = work_dir), expand("{work_dir}/reads/raw__filtered/{sample}_R2.fastq.gz", sample = samples, work_dir = work_dir)

# gzip'ing obtained reads
rule pigz_reads:
    input:
        reads = ["{work_dir}/reads/raw__filtered/{sample}_R1.fastq", "{work_dir}/reads/raw__filtered/{sample}_R2.fastq"],
    output:
        gzipped = ["{work_dir}/reads/raw__filtered/{sample}_R1.fastq.gz", "{work_dir}/reads/raw__filtered/{sample}_R2.fastq.gz"],
    conda: "/path/to/envs/pigz.yaml"
    threads: 10
    shell:
        "pigz -p {threads} {input.reads}"

# get decontaminated reads 
rule get_unmapped_reads:
    input:
        "{work_dir}/reads/raw__filtered/{sample}_unmapped_sorted.bam",
    output:
        R1 = temp("{work_dir}/reads/raw__filtered/{sample}_R1.fastq"),
        R2 = temp("{work_dir}/reads/raw__filtered/{sample}_R2.fastq")
    conda: "/path/to/envs/bedtools.yaml"
    shell:
        "bedtools bamtofastq -i {input} -fq {output.R1} -fq2 {output.R2}"

# get unmmaped bam
rule get_unmapped_bam:
    input:
        bam = "{work_dir}/reads/raw__filtered/{sample}.bam"
    output: unmapped_bam = temp("{work_dir}/reads/raw__filtered/{sample}_unmapped_sorted.bam")
    conda: 
        "/path/to/envs/samtools.yaml"
    shell:
        "samtools view -b -f 12 -F 256 {input} | samtools sort > {output}"

# HiSat2: remove human
rule hisat2_remove_human:
    input: reads = ["{work_dir}/reads/raw__fastp/{sample}_R1.fastq.gz", "{work_dir}/reads/raw__fastp/{sample}_R2.fastq.gz"],
        idx = hisat2_index
    output:
        temp("{work_dir}/reads/raw__filtered/{sample}.bam"),
    log:
        "{work_dir}/reports/hisat2/{sample}.report",
    params:
        extra = "--very-sensitive",
    threads: 10
    wrapper:
        "v2.6.0/bio/hisat2/align"

# fastp: quality filtering with average Q > 30 and autodetection of adapter content
rule fastp_quality_filtering:
    input:
        sample = ["{work_dir}/reads/raw/{sample}_R1.fastq.gz", "{work_dir}/reads/raw/{sample}_R2.fastq.gz"]
    output:
        trimmed = temp(["{work_dir}/reads/raw__fastp/{sample}_R1.fastq.gz", "{work_dir}/reads/raw__fastp/{sample}_R2.fastq.gz"]),
        html = "{work_dir}/reports/fastp/{sample}_fastp.html",
        json = "{work_dir}/reports/fastp/{sample}_fastp.json"
    params:
        extra = "--detect_adapter_for_pe --overrepresentation_analysis --html {output.html} --correction --dup_calc_accuracy 6 --average_qual 30"
    threads: 10
    wrapper:
        "v2.6.0/bio/fastp"
