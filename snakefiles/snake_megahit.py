import glob
import numpy as np

work_dir='path/to/data'

samples = [s.split('/')[-1].split('_R')[0] for s in glob.glob(work_dir + '/reads/raw__filtered/*_R1.fastq.gz')]
print(samples)

rule run_megahit:
	input: 
		expand("{work_dir}/analysis/megahit/{sample}/final.contigs.fa", sample = samples, work_dir = work_dir)

rule megahit:
	input: 
		R1 = "{work_dir}/reads/raw__filtered/{sample}_R1.fastq.gz", 
		R2 = "{work_dir}/reads/raw__filtered/{sample}_R2.fastq.gz"
	threads: 10
	output: "{work_dir}/analysis/megahit/{sample}/final.contigs.fa"
	conda: "/path/to/envs/megahit.yaml"
	shell: "megahit -f -t {threads} -m 1000e9 -1 {input.R1} -2 {input.R2} -o {work_dir}/analysis/megahit/{wildcards.sample}"
