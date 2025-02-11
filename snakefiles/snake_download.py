import pandas as pd

samples = pd.read_csv("/path/to/table/sampleid.txt")
samples = samples['accession'].tolist()
print(samples)

rule run_download:
    input: expand("/path/to/data/{sample}_2.fastq.gz", sample = samples)

rule get_fastq:
    output: "/path/to/data/{sample}_2.fastq.gz"
    conda: "/path/to/envs/kingfisher.yaml"
    threads: 10
    shell: "kingfisher get -r {wildcards.sample} -m aws-http -f fastq.gz --download-threads {threads} --output-directory /path/to/data"
