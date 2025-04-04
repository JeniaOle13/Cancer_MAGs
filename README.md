# Catalog of stool metagenome-assembled genomes from patients with different cancer types

## Motivation
Immune checkpoint inhibitors are a cancer treatment approach that aims to activate the body's own immunity against the tumor. This approach can fight cancer even in advanced stages. However, not every patient can benefit from immunotherapy. Numerous studies of patients undergoing cancer immunotherapy describe the influence of the gut microbiota on anti-tumor immunity and immunotherapy efficacy, as well as changes in the microbial profile of patients who respond to treatment and those who do not.

Our computational pipeline makes a notable contribution to the expanding field of research investigating the role of the gut microbiome in cancer immunotherapy. Despite extensive studies of the microbiome in the context of melanoma, the gut microbiota of patients with other types of cancer remains poorly understood. The catalog of high-quality metagenome-assembled genomes (MAGs) we present here helps to fill this gap.

The resulting catalog is valuable for studying the composition and functional characteristics of the microbiome, which may have significant potential implications for future microbiome research and the development of personalized approaches to cancer immunotherapy.

## Description
A non-redundant catalog of 3,816 genomes with average 91.9 ± 6.8% completeness and 2.13 ± 2.77% contamination assembled from metagenomes has been introduced. Samples of 976 metagenomes from 14 studies were obtained from patients receiving immunotherapy for the treatment of different types of cancers.

### Accession codes of the used stool metagenome sequences
[PRJNA397906](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA397906), [PRJEB22893](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJEB22893), [PRJNA399742](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA399742), [PRJNA678737](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA678737), [PRJNA672867](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA672867), [PRJNA770295](https://www.ncbi.nlm.nih.gov/bioproject/770295), [PRJEB43119](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJEB43119), [PRJNA762360](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA762360), [PRJNA1011235](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA1011235), [PRJNA928744](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA928744), [PRJNA615114](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA615114), [PRJNA866654](https://www.ncbi.nlm.nih.gov/bioproject/866654), [PRJNA494824](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJNA494824), [PRJEB49516](https://www.ncbi.nlm.nih.gov/bioproject/?term=PRJEB49516).

### Computational pipeline
Description of computational steps is available at https://github.com/JeniaOle13/Cancer_MAGs/blob/main/PIPELINE.md.

![](https://github.com/kanaevavera/Cancer_MAGs/blob/main/catalog_picture.png)
**Approximate maximum likelihood phylogenetic tree generated using CheckM with 43 AA marker sequences and 3,816 MAGs assembled from the stool metagenomes of cancer patients.**

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
