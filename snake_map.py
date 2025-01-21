import glob
import numpy as np
import re

work_dir='path/to/data'

samples = [s.split('/')[-1].split('_R1')[0] for s in glob.glob(work_dir + '/reads/raw__filtered/*_R1.fastq.gz')]
print(samples)

rule run_map:
        input: expand("{work_dir}/analysis/megahit/{sample}/{sample}_sorted.bam", sample = samples, work_dir=work_dir)

rule sort_unmapped_bam:
        input: "{work_dir}/analysis/megahit/{sample}/{sample}.bam"
        output: "{work_dir}/analysis/megahit/{sample}/{sample}_sorted.bam"
        conda: "/data12/bio/runs-jeniaole/snake-files/MAG_assembly_pipeline/envs/samtools.yaml"
        shell: "samtools sort {input} -o {output}"

rule sam_to_bam:
        input: "{work_dir}/analysis/megahit/{sample}/{sample}.sam"
        output: temp("{work_dir}/analysis/megahit/{sample}/{sample}.bam")
        conda: "/data12/bio/runs-jeniaole/snake-files/MAG_assembly_pipeline/envs/samtools.yaml"
        shell: "samtools view -Sb -F 0x4 {input} > {output}"

rule hisat2_map:
        input: R1 = "{work_dir}/reads/raw__filtered/{sample}_R1.fastq.gz", R2 = "{work_dir}/reads/raw__filtered/{sample}_R2.fastq.gz", index_file = "{work_dir}/analysis/megahit/{sample}/contigs.8.ht2"
        output: sam = temp("{work_dir}/analysis/megahit/{sample}/{sample}.sam"), statistics = "{work_dir}/analysis/megahit/{sample}/{sample}.hisat"
        threads: 10
        params: prefix = "{work_dir}/analysis/megahit/{sample}/contigs"
        conda: "/data12/bio/runs-jeniaole/snake-files/MAG_assembly_pipeline/envs/hisat2.yaml"
        shell: "hisat2 --very-sensitive -p {threads} -x {params.prefix} -1 {input.R1} -2 {input.R2} -S {output.sam} 2> {output.statistics}"

rule hisat2_index:
        input: "{work_dir}/analysis/megahit/{sample}/contigs.fa"
        output:
            "{work_dir}/analysis/megahit/{sample}/contigs.1.ht2",
            "{work_dir}/analysis/megahit/{sample}/contigs.2.ht2",
            "{work_dir}/analysis/megahit/{sample}/contigs.3.ht2",
            "{work_dir}/analysis/megahit/{sample}/contigs.4.ht2",
            "{work_dir}/analysis/megahit/{sample}/contigs.5.ht2",
            "{work_dir}/analysis/megahit/{sample}/contigs.6.ht2",
            "{work_dir}/analysis/megahit/{sample}/contigs.7.ht2",
            "{work_dir}/analysis/megahit/{sample}/contigs.8.ht2"
        params: prefix = "{work_dir}/analysis/megahit/{sample}/contigs"
        conda: "/data12/bio/runs-jeniaole/snake-files/MAG_assembly_pipeline/envs/hisat2.yaml"
        shell: "hisat2-build {input} {params.prefix}"

rule reformat_contigs:
        input: "{work_dir}/analysis/megahit/{sample}/final.contigs.fa"
        output: "{work_dir}/analysis/megahit/{sample}/contigs.fa"
        conda: "/data12/bio/runs-jeniaole/snake-files/MAG_assembly_pipeline/envs/bbmap.yaml"
        shell: "reformat.sh in={input} out={output} minlength=1000"
