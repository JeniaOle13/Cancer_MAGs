import glob

work_dir='path/to/data'

samples = [s.split('/')[-1].split('_R')[0] for s in glob.glob(work_dir +'/reads/raw__filtered/*_R1.fastq.gz')]
print("samples:", samples)

rule run_binning:
	input: expand("{work_dir}/analysis/metabat2/{sample}/{sample}.stat", sample = samples, work_dir = work_dir), expand("{work_dir}/analysis/semibin2/{sample}/contig_bins.tsv", sample = samples, work_dir = work_dir), expand("{work_dir}/analysis/maxbin2/{sample}/{sample}.log", sample = samples, work_dir = work_dir), expand("{work_dir}/analysis/dastool/{sample}.stat", sample = samples, work_dir = work_dir)

rule dastool:
    input: 
        metabat2 = "{work_dir}/analysis/metabat2/{sample}/{sample}.metabat2.contigs2bin.tsv",
        semibin2 = "{work_dir}/analysis/semibin2/{sample}/{sample}.semibin2.contigs2bin.tsv",
        maxbin2 = "{work_dir}/analysis/maxbin2/{sample}/{sample}.maxbin2.contigs2bin.tsv",
        contigs = "{work_dir}/analysis/megahit/{sample}/contigs.fa"
    output: 
        statistic = "{work_dir}/analysis/dastool/{sample}.stat"
    conda: "/data12/bio/runs-kanaevavera/envs/dastool.yaml"
    threads: 10
    shell: "DAS_Tool --write_bin_evals --write_bins --threads {threads} --score_threshold 0.5 -i {input.metabat2},{input.semibin2},{input.maxbin2} -l metabat,semibin,maxbin -c {input.contigs} -o {work_dir}/analysis/dastool/{wildcards.sample}/{wildcards.sample} &>> {output.statistic}"

rule get_maxbin2_contigs_stats:
    input: "{work_dir}/analysis/maxbin2/{sample}/{sample}.log"
    params: "{work_dir}/analysis/maxbin2/{sample}"
    output: "{work_dir}/analysis/maxbin2/{sample}/{sample}.maxbin2.contigs2bin.tsv"
    conda: "/data12/bio/runs-kanaevavera/envs/dastool.yaml"
    shell: "Fasta_to_Contig2Bin.sh -i {params} -e fasta > {output}"

rule get_metabat2_contigs_stats:
    input: "{work_dir}/analysis/metabat2/{sample}/{sample}.stat"
    params: "{work_dir}/analysis/metabat2/{sample}"
    output: "{work_dir}/analysis/metabat2/{sample}/{sample}.metabat2.contigs2bin.tsv"
    conda: "/data12/bio/runs-kanaevavera/envs/dastool.yaml"
    shell: "Fasta_to_Contig2Bin.sh -i {params} -e fa > {output}"

rule get_semibin2_contigs_stats:
    input: "{work_dir}/analysis/semibin2/{sample}/contig_bins.tsv"
    params: "{work_dir}/analysis/semibin2/{sample}/output_bins"
    output: "{work_dir}/analysis/semibin2/{sample}/{sample}.semibin2.contigs2bin.tsv"
    conda: "/data12/bio/runs-kanaevavera/envs/dastool.yaml"
    shell: "Fasta_to_Contig2Bin.sh -i {params} -e fa > {output}"

# running semibin2
rule semibin2:
	input: contigs = "{work_dir}/analysis/megahit/{sample}/contigs.fa", bam = "{work_dir}/analysis/megahit/{sample}/{sample}_sorted.bam"
	output: "{work_dir}/analysis/semibin2/{sample}/contig_bins.tsv"
	params: "{work_dir}/analysis/semibin2/{sample}"
	conda: "/data12/bio/runs-kanaevavera/envs/semibin2.yaml"
	threads: 20
	shell: "SemiBin2 single_easy_bin -i {input.contigs} -b {input.bam} -o {params} --environment global --compression none -p {threads}"

# running metabat2
rule metabat2:
	input: 
         depth = "{work_dir}/analysis/megahit/{sample}/{sample}.depth.txt", 
         contigs = "{work_dir}/analysis/megahit/{sample}/contigs.fa"
	output: statistic = "{work_dir}/analysis/metabat2/{sample}/{sample}.stat"
	conda: "/data12/bio/runs-kanaevavera/envs/metabat2.yaml"
	shell: "metabat2 -i {input.contigs} -a {input.depth} -o {work_dir}/analysis/metabat2/{wildcards.sample}/{wildcards.sample} -m 1500 &>> {output.statistic}"

rule get_depth_for_metabat2:
	input: "{work_dir}/analysis/megahit/{sample}/{sample}_sorted.bam"
	output: "{work_dir}/analysis/megahit/{sample}/{sample}.depth.txt"
	conda: "/data12/bio/runs-kanaevavera/envs/metabat2.yaml"
	shell: "jgi_summarize_bam_contig_depths --outputDepth {output} {input}"

# running maxbin2
rule maxbin2:
        input:
                abund = "{work_dir}/analysis/megahit/{sample}/{sample}.abundance.txt",
                contigs = "{work_dir}/analysis/megahit/{sample}/contigs.fa",
                R1 = "{work_dir}/reads/raw__filtered/{sample}_R1.fastq.gz",
                R2 = "{work_dir}/reads/raw__filtered/{sample}_R2.fastq.gz"
        threads: 20
        output: "{work_dir}/analysis/maxbin2/{sample}/{sample}.log",
        conda: "/data12/bio/runs-kanaevavera/envs/maxbin2.yaml"
        shell: "run_MaxBin.pl -reads {input.R1} -reads2 {input.R2} -thread {threads} -contig {input.contigs} -out {work_dir}/analysis/maxbin2/{wildcards.sample}/{wildcards.sample} -abund {input.abund} &>> {output}"

rule mkdir:
    input: "{work_dir}/analysis/megahit/{sample}/{sample}.abundance.txt"
    output: directory('{work_dir}/analysis/maxbin2/{wildcards.sample}/')
    shell: "mkdir {output}"

rule get_abundance_for_maxbin2:
        input: "{work_dir}/analysis/megahit/{sample}/{sample}.cov.txt"
        output: "{work_dir}/analysis/megahit/{sample}/{sample}.abundance.txt"
        shell: """awk '{{print $1"\\t"$5}}' {input} | grep -v '^#' > {output}"""

rule get_coverage_for_maxbin2:
    input: "{work_dir}/analysis/megahit/{sample}/{sample}_sorted.bam"
    output: "{work_dir}/analysis/megahit/{sample}/{sample}.cov.txt"
    conda: "/data12/bio/runs-kanaevavera/envs/bbmap.yaml"
    shell: "pileup.sh in={input} out={output} 32bit=t"
