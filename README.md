# Catalog of stool metagenome-assembled genomes from patients with different cancer types

## Description
A non-redundant catalog of 3,816 genomes with average 91.9 ± 6.8% completeness and 2.13 ± 2.77% contamination assembled from metagenomes has been introduced. Samples of 976 metagenomes from 14 studies were obtained from patients receiving immunotherapy for the treatment of different types of cancers.

### Accession codes of the used stool metagenome sequences
[PRJNA397906](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA397906), [PRJEB22893](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJEB22893), [PRJNA399742](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA399742), [PRJNA678737](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA678737), [PRJNA672867](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA672867), [PRJNA770295](https://www.ncbi.nlm.nih.gov/bioproject/770295), [PRJEB43119](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJEB43119), [PRJNA762360](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA762360), [PRJNA1011235](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA1011235), [PRJNA928744](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA928744), [PRJNA615114](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA615114), [PRJNA866654](https://www.ncbi.nlm.nih.gov/bioproject/866654), [PRJNA494824](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA494824), [PRJEB49516](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJEB49516).

![](https://lh7-rt.googleusercontent.com/docsz/AD_4nXeJy0hY1mi8HOGaU633bZqRL4Laruz9jO__GE7g7XHasyAmuzqyUwLX-puD2bmTM6VaRFnqXBFjwauj3bEqOZLBm_F5oMPfhnmXdtSul6V-A1TEzylWpLMmc4m0r0xOHFgjh0EQ?key=1_Rc7IYQvnC-DPwR2bAX7g)
**Approximate maximum likelihood phylogenetic tree generated using CheckM with 43 AA marker sequences and 3,816 MAGs assembled from the stool metagenomes of cancer patients.**

## Dependencies for computational pipeline
| **Tool**  | **Publication**  |
|---|---|
|[conda](https://github.com/conda/conda) v23.7.4| - |
|[mamba](https://github.com/mamba-org/mamba) v1.5.1| - |
|[snakemake](https://github.com/snakemake/snakemake) v7.32.4|[https://doi.org/10.12688/f1000research.29032.2](https://doi.org/10.12688/f1000research.29032.2)|
|[kingfisher](https://github.com/wwood/kingfisher-download) v0.4.1| - |
|[samtools](https://github.com/samtools/samtools) v1.17|[https://doi.org/10.1093/gigascience/giab008](https://doi.org/10.1093/gigascience/giab008)|
|[bedtools](https://github.com/arq5x/bedtools2) v2.31.0|[https://doi.org/10.1093/bioinformatics/btq033](https://doi.org/10.1093/bioinformatics/btq033)|
|[pigz](https://github.com/madler/pigz) v2.6| - |
|[fastp](https://github.com/OpenGene/fastp) v0.23.4|[https://doi.org/10.1002/imt2.107](https://doi.org/10.1002/imt2.107)|
|[hisat2](https://github.com/DaehwanKimLab/hisat2) v2.2.1|[https://doi.org/10.1038/s41587-019-0201-4](https://doi.org/10.1038/s41587-019-0201-4)|
|[bbmap](https://sourceforge.net/projects/bbmap/) v39.06| - |
|[megahit](https://github.com/voutcn/megahit) v1.2.9|[https://doi.org/10.1093/bioinformatics/btv033](https://doi.org/10.1093/bioinformatics/btv033)|
|[metabat2](https://bitbucket.org/berkeleylab/metabat/src/master/) v2.12.1|[https://doi.org/10.7717/peerj.7359](https://doi.org/10.7717/peerj.7359)|
|[maxbin2](https://sourceforge.net/projects/maxbin2/) v2.2.7|[https://doi.org/10.1093/bioinformatics/btv638](https://doi.org/10.1093/bioinformatics/btv638)|
|[semibin2](https://github.com/BigDataBiology/SemiBin) v2.1.0|[https://doi.org/10.1093/bioinformatics/btad209](https://doi.org/10.1093/bioinformatics/btad209)|
|[dastool](https://github.com/cmks/DAS_Tool) v1.1.7|[https://doi.org/10.1038/s41564-018-0171-1](https://doi.org/10.1038/s41564-018-0171-1)|
|[drep](https://github.com/MrOlm/drep) v3.4.5|[https://doi.org/10.1038/ismej.2017.126](https://doi.org/10.1038/ismej.2017.126)|
|[gtdbtk](https://github.com/Ecogenomics/GTDBTk) v2.1.1|[https://doi.org/10.1093/bioinformatics/btz848](https://doi.org/10.1093/bioinformatics/btz848)|
|[fasttree](https://morgannprice.github.io/fasttree/) v2.1.11|[https://doi.org/10.1371/journal.pone.0009490](https://doi.org/10.1371/journal.pone.0009490)|
|[qiime2](https://view.qiime2.org/) v2022.8.3|[https://doi.org/10.1038/s41587-019-0209-9](https://doi.org/10.1038/s41587-019-0209-9)|
|[empress](https://github.com/biocore/empress) v1.2.0|[https://doi.org/10.1128/msystems.01216-20](https://doi.org/10.1128/msystems.01216-20)|

## Computational pipeline
The `conda` and/or `mamba` package management systems must be installed to run computational scripts. It is also required to create a snakemake working environment.

> conda create -n snakemake -c bioconda snakemake

> conda activate snakemake

Before running scripts, rewrite `work_dir` path the paths to your project in all snakemake scripts. Also you need to specify paths for conda environments (`/envs/*yaml files`).

The organization of folders inside the project folder should look like this:

> /path/to/data/reads/raw/{sample}_*.fastq.gz

All computations were performed on the Lopukhin FRCC PCM cluster with 2 x AMD EPYC 7763 and 2 Tb RAM.

### 1. Download samples from SRA
The IDs of the metagenomic samples used in the study are provided in the `sampleid.txt` file. Use `snake_download.py` to download files from the repository. To start the download with 60 threads, use following command after activating `snakemake` environment:

> snakemake -s snake_download_test.py -k -p --latency-wait 150 -j 10 --use-conda run_download

All downloaded files will be stored in `path/to/data/reads/raw` directory.

### 2. Pre-processing
Metagenomic reads will be filtered using fastp with the following parameters `--detect_adapter_for_pe` `--overrepresentation_analysis` `--correction --dup_calc_accuracy 6` `--average_qual 30`. Quality controlled reads are mapped to the human genome using `hisat2` tool with `--very-sensitive` parameter. After that, `samtools` were used for filtering unmapped reads, `bedtools` for converting unmapped `bam` files to the `fastq` format as well as `pigz` gzipping the resulting `fastq` files. The results of the module will be gzipped fastq files, while all intermediate files will be removed:

> /path/to/data/reads/raw__filtered/{sample}_R*.fastq.gz

The human reference genome [GRCh38](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.26/) required for this step. It is also required to create hisat2 index:

> hisat2-build genome.fa genome

Within the `snake_preprocessing.py` script you need to specify the directory where the index was created `hisat2_index=”path/to/directory/with/index”`.

The script is run from the command line with the activated snakemake environment using 10 threads:

> snakemake -s snake_preprocessing.py -k -p --latency-wait 150 -j 10 --use-conda run_preproc

### 3. Assembly
Assembly step was performed using `MEGAHIT` short reads assembler with default parameters. To perform a metagenomic assembly, run the script by command line with the snakemake environment activated:

> snakemake -s snake_megahit.py -k -p --latency-wait 150 -j 10 --use-conda run_megahit

Results were storage in separation folder:

> /path/to/data/analysis/megahit/{sample}/final.contigs.fa

### 4. Binning
Metagenomic binning is performed in two steps. The first step will involve mapping pre-processed reads obtained on the step 1 to the metagenomic assembly obtained on the step 2. Firstly, all contigs with length less than 1000 bp are removed from the `final.contig.fa` files using the `reformat.sh` script from the `bbmap` package. The resulting filtered file `contigs.fa` is used by the `hisat2-build` program to build a reference to perform short read mapping. After that, `hisat2` with `--very-sensitive` flag was used for mapping reads to reference as well as `samtools` used for converting resulting `sam` files to sorted `bam` files. All intermediate files were removed from hard disk space.

To perform short reads mapping to metagenomic contigs, execute on the command line with the `snakemake` environment activated:

> snakemake -s snake_map.py -k -p --latency-wait 150 -j 10 --use-conda run_map

The result of the computational step is the production of sorted mapping files in `bam` format:

> /path/to/data/analysis/megahit/{sample}_sorted.bam

The second step is to perform metagenomic binning using metagenomic contigs `contigs.fa` and `bam` files obtained in the previous steps of the computational steps. Firstly, `maxbin2` and `metabat2` tools required files with abundance statistics for each contig. For `maxbin2` `pileup.sh` script from `bbmap` package and `awk` were used for the obtained coverage file. Metabat2 required the coverage depth file which can be obtained using `jgi_summarize_bam_contig_depths` script. `SemiBin2` can be run using `contigs.fa` and bam files and requires no prior additional computation. After performing the binning procedure using the three programs, `dastool` is used to form an optimized, non-redundant set of bins. Firstly, contigs to bin mapping files for all binning programs must be generated using `Fasta_to_Contig2Bin.sh`. After that `dastool` runs with `--score_threshold 0.5`.

To run binning procedure described above execute following script with activated snakemake environment:

> snakemake -s snake_binning.py -k -p --latency-wait 150 -j 10 --use-conda run_binning

The result of the computations will be metagenomic bins generated by each of the programs:

> /path/to/data/analysis/maxbin2/{sample}/*fasta
> /path/to/data/analysis/metabat2/{sample}/*fa
> /path/to/data/analysis/semibin2/{sample}/output_bins/*fa

Also an optimized, non-redundant set of bins will be storage:

> /path/to/data/analysis/dastool/{sample}/{sample}_DASTool_bins/*fa

### 5. Dereplication:
In order to prepare the received bins for dereplication, it is necessary to run a script for renamed bins and forming links for `dRep`. The code is designed to create a folder named `/path/to/data/analysis/dastool_bins` and to establish symbolic links to bins that bear unified names. These names reflect the sample ID and the name of the binning program that generated the bin. Flag `--work_dir` accepts the full path to the `/path/to/data/analysis` folder of your project.

To run script execute in command line:

> python prepare_bins.py --work_dir /path/to/data/analysis

These required `dRep` tool which must be installed via `conda` or `mamba`:

> conda create -n drep -c bioconda drep

> conda activate drep

To run dereplication execute in command line with activated `dRep conda` environment:

> dRep dereplicate -g /path/to/data/analysis/dastool_bins/*fa --P_ani 0.9 --S_ani 0.98 -p 10 /path/to/data/analysis/drep

Results:

> /path/to/data/analysis/drep/dereplicated_genomes/*fa

### 6. Taxonomic annotation
These required `gtdb-tk` tool which must be installed via conda or mamba:

> conda create -n gtdbtk -c bioconda gtdbtk

> conda activate gtdbtk

Before running `gtdb-tk` taxonomic annotation, Genome Taxonomy Database (GTDB) can be downloaded from https://gtdb.ecogenomic.org/.

To run taxonomic annotation of resulting metagenome-assembled genome catalog execute in command line with gtdb-tk conda activated environment:

> gtdbtk classify_wf --genome_dir /path/to/data/analysis/drep/dereplicated_genomes/ -x fa --out_dir /path/to/data/analysis/gtdbtk --cpus 10 --pplacer_cpus 10 --mash_db /path/to/gtdbtk/database

Results will be storage in folder:
> /path/to/data/analysis/gtdbtk/classify/

### 7. Aggregation the MAGs metadata:
This item provides a Python script that summarizes the results of the previous steps.  `make_mags_meta.py` generates the metadata table for metagenome-assembled genomes (MAGs), including genome qualities and taxonomy annotations. Pass the full path to your project's `/path/to/data/analysis` folder to the `--work_dir` flag.

To execute run in command line:

> python make_mags_meta.py --work_dir /path/to/data/analysis

### 8. Phylogenetic tree construction and visualisation:
These require the `CheckM` tool which is included in the `dRep` `conda` environment. The `checkm_tree.sh` script is employed to generate protein translation to the dereplicated genomes, followed by the execution of multiple sequence alignment of reference genomes and MAGs sequences. These processes facilitate the inference of a *de novo* genome tree. It is imperative to provide the full path to the `/path/to/data/analysis` folder of the project for script execution.

To run placing MAGs in the reference genome tree described above execute in command line with activated `drep` `conda` environment:

> bash checkm_tree.sh /path/to/data/analysis

The subsequent step is essential for those who want to visualize only the MAGs without the reference genomes. The utilization of the `sed` command entails the elimination of the reference IMG bins from the `checkm_alignment` file, thereby giving rise to the `checkm_alignment_woRef` file:

> sed '/>IMG/,+1 d' /path/to/data/analysis/phylo_tree/checkm_alignment > /path/to/data/analysis/phylo_tree/checkm_alignment_woRef

The next step involves the de novo construction of an unrooted tree using the `fasttree` which must be installed via `conda` or `mamba`:

> conda create -n fasttree -c bioconda fasttree

> conda activate fasttree

To run phylogenetic tree construction execute in command line with activated fasttree conda environment:
> fasttree < /path/to/data/analysis/phylo_tree/checkm_alignment_woRef > /path/to/data/analysis/phylo_tree/fasttree_aa

The following script utilizes `Qiime2` and its `Empress` package to generate a tree visualization. This implementation facilitates the acquisition of a rooted tree file in `Qiime2` format and the visualization of the data in the `.qzv` format. The visualization can be accessed via the `Qiime2` viewer at https://view.qiime2.org/. `Qiime2` tool must be installed via `conda` or `mamba` (Instructions for installing latest versions of QIIME2 can be found at the following URL: https://docs.qiime2.org/):
> wget https://data.qiime2.org/distro/core/qiime2-2022.8-py38-linux-conda.yml

> conda env create -n qiime2-2022.8 --file qiime2-2022.2-py38-linux-conda.yml

> conda activate qiime2-2022.8

`Empress` should be installed using `pip` with activated `Qiime2` `conda` environment:
> pip install empress

The script is necessary for the metadata file from step 6. Additionally, it is necessary to provide the full path to the `/path/to/data/analysis` folder of the project for script execution. To execute run in command line:
> bash empress_vis.sh /path/to/data/analysis

## DATA
|Description|Size|Links|
|---|---|---|
|all MAGs sequences (n=13,227)|11 GB|[Lopukhin FRCC PCM FTP](http://download.ripcm.com/mag_catalog/)|
|high-quality MAGs sequences (n=9,156)|7.4 GB|[Lopukhin FRCC PCM FTP](http://download.ripcm.com/mag_catalog/)|
|non-redundant MAGs catalog (n = 3,816)|2.7 GB|[Lopukhin FRCC PCM FTP](http://download.ripcm.com/mag_catalog/)|
|Merged non-redundant MAGs catalog|2.7 GB|[Zenodo](https://zenodo.org/records/14701025/files/dereplicated.tar.gz?download=1); [Figshare](https://figshare.com/ndownloader/files/51718619)|
|Mapping file (MAG>contig)|43 MB|[Zenodo](https://zenodo.org/records/14701025/files/mapping.tsv?download=1); [Figshare](https://figshare.com/ndownloader/files/51718547)|
|Dereplication statistic|940 KB|[Zenodo](https://zenodo.org/records/14358262/files/dereplication_stats.tsv?download=1); [Figshare](https://figshare.com/ndownloader/files/51079292)|
|Taxonomic annotation and assembly quality statistic|645 KB|[Zenodo](https://zenodo.org/records/14358262/files/statistics.tsv?download=1); [Figshare](https://figshare.com/ndownloader/files/51062807)|
|DNA extraction and sequencing methods|13 KB|[Zenodo](https://zenodo.org/records/14358262/files/methods.docx?download=1); [Figshare](https://figshare.com/ndownloader/files/51081881)|

All MAGs sequences were deposited to NCBI under BioProject [PRJNA1196825](https://www.ncbi.nlm.nih.gov/bioproject/1196825).

## Studies links:
1) [Frankel, Arthur E., et al. "Metagenomic shotgun sequencing and unbiased metabolomic profiling identify specific human gut microbiota and metabolites associated with immune checkpoint therapy efficacy in melanoma patients." Neoplasia 19.10 (2017): 848-855.](https://www.sciencedirect.com/science/article/pii/S1476558617302385)
2) [Peng, Zhi, et al. "The gut microbiome is associated with clinical response to anti–PD-1/PD-L1 immunotherapy in gastrointestinal cancer." Cancer Immunology Research 8.10 (2020): 1251-1261.](https://aacrjournals.org/cancerimmunolres/article/8/10/1251/466881/The-Gut-Microbiome-Is-Associated-with-Clinical)
3) [Gopalakrishnan, Vancheswaran, et al. "Gut microbiome modulates response to anti–PD-1 immunotherapy in melanoma patients." Science 359.6371 (2018): 97-103.](https://www.science.org/doi/10.1126/science.aan4236)
4) [Matson, Vyara, et al. "The commensal microbiome is associated with anti–PD-1 efficacy in metastatic melanoma patients." Science 359.6371 (2018): 104-108.](https://www.science.org/doi/10.1126/science.aao3290)
5) [Baruch, Erez N., et al. "Fecal microbiota transplant promotes response in immunotherapy-refractory melanoma patients." Science 371.6529 (2021): 602-609.](https://www.science.org/doi/10.1126/science.abb5920)
6) [Davar, Diwakar, et al. "Fecal microbiota transplant overcomes resistance to anti–PD-1 therapy in melanoma patients." Science 371.6529 (2021): 595-602.](https://www.science.org/doi/10.1126/science.abf3363)
7) [Spencer, Christine N., et al. "Dietary fiber and probiotics influence the gut microbiome and melanoma immunotherapy response." Science 374.6575 (2021): 1632-1640.](https://www.science.org/doi/10.1126/science.aaz7015?url_ver=Z39.88-2003&rfr_id=ori:rid:crossref.org&rfr_dat=cr_pub%20%200pubmed)
8) [Lee, Karla A., et al. "Cross-cohort gut microbiome associations with immune checkpoint inhibitor response in advanced melanoma." Nature Medicine 28.3 (2022): 535-544.](https://www.nature.com/articles/s41591-022-01695-5)
9) [McCulloch, John A., et al. "Intestinal microbiota signatures of clinical response and immune-related adverse events in melanoma patients treated with anti-PD-1." Nature Medicine 28.3 (2022): 545-556.](https://www.nature.com/articles/s41591-022-01698-2)
10) [Tsakmaklis, Anastasia, et al. "TIGIT+ NK cells in combination with specific gut microbiota features predict response to checkpoint inhibitor therapy in melanoma patients." BMC Cancer 23.1 (2023): 1160.](https://bmccancer.biomedcentral.com/articles/10.1186/s12885-023-11551-5)
11) [Routy, Bertrand, et al. "Fecal microbiota transplantation plus anti-PD-1 immunotherapy in advanced melanoma: a phase I trial." Nature Medicine 29.8 (2023): 2121-2132.](https://www.nature.com/articles/s41591-023-02453-x)
12) [Liu, Ben, et al. "Exploring gut microbiome in predicting the efficacy of immunotherapy in non-small cell lung cancer." Cancers 14.21 (2022): 5401.](https://www.mdpi.com/2072-6694/14/21/5401)
13) [Heshiki, Yoshitaro, et al. "Predictable modulation of cancer treatment outcomes by the gut microbiota." Microbiome 8 (2020): 1-14.](https://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-020-00811-2)
14) [Gunjur, Ashray, et al. "A gut microbial signature for combination immune checkpoint blockade across cancer types." Nature Medicine 30.3 (2024): 797-809.](https://www.nature.com/articles/s41591-024-02823-z)

## Data links:
1) Lopukhin FRCC PCM FTP: [http://download.ripcm.com/mag_catalog](http://download.ripcm.com/mag_catalog/).
2) Zenodo: [https://doi.org/10.5281/zenodo.14358262](https://doi.org/10.5281/zenodo.14358262).
3) Figshare: [https://doi.org/10.6084/m9.figshare.27993506.v4](https://doi.org/10.6084/m9.figshare.27993506.v4).

## Contacts:
Also, you can send your feedback to jeniaole01@gmail.com, vera.a.kanaeva@gmail.com.

## Authors:
MAGs assembly, pipeline developer: [Vera Kanaeva](https://scholar.google.ru/citations?hl=ru&user=Ie7RMLAAAAAJ) (Lopukhin FRCC PCM, MIPT).

Idea, supervisor, pipeline developer: [Evgenii Olekhnovich](https://scholar.google.ru/citations?user=RA9ItlsAAAAJ&hl=ru) (Lopukhin FRCC PCM).

This approach was first used for strain profiling in our publication:
[Zakharevich, Natalia V., et al. "Systemic metabolic depletion of gut microbiome undermines responsiveness to melanoma immunotherapy." Life Science Alliance 7.5 (2024).](https://www.life-science-alliance.org/content/7/5/e202302480/tab-figures-data)

## Fundings:
Financial support for this study was provided by the Russian Science Foundation under the grant #22-75-10029 (https://rscf.ru/project/22-75-10029/).
