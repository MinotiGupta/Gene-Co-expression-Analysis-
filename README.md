# Gene Co-expression Network Analysis 

## Overview
This project, conducted at **Amrita School of Artificial Intelligence, Amrita Vishwa Vidyapeetham, Coimbatore**, focuses on **Gene Co-expression Network Analysis (GCNA)** to investigate gene expression patterns related to **COVID-19**. The analysis uses RNA-sequencing data from the **NCBI Gene Expression Omnibus (GEO) dataset GSE151418** to construct co-expression networks, identify key gene modules, and perform functional enrichment analysis. The project also addresses **ethical considerations** and **intellectual property rights** in genomic research.

## Abstract
Gene Co-expression Network Analysis (GCNA) is employed to uncover relationships among genes and their roles in biological processes and diseases, specifically targeting **COVID-19**. Using **RNA-sequencing data** from GSE151418, the project constructs co-expression networks and conducts **enrichment analysis** to identify **key gene modules** (green and turquoise) and **hub genes** like **H2BC9**, **ITGB3**, and **IGHV2026**. These findings highlight potential **biomarkers** and **therapeutic targets** for COVID-19, with **GO**, **KEGG**, **Reactome**, and **WikiPathways** revealing processes such as **nucleosome assembly**, **immune response**, and **blood clotting**. Ethical issues, including **informed consent** and **health disparities**, are also explored alongside **intellectual property rights** in genomic research.

## Introduction
GCNA is a systems biology tool that maps **coordinated gene expression patterns** to reveal biological processes and disease mechanisms. Using **Weighted Gene Co-expression Network Analysis (WGCNA)**, the project analyzes **RNA-seq data** to identify **gene modules** and **hub genes** critical to **COVID-19 pathogenesis**. **Enrichment analysis** with databases like **GO**, **KEGG**, **Reactome**, and **WikiPathways** annotates these modules, linking them to functions such as **immune pathways** and **DNA organization**. The integrative approach bridges **transcriptomic data** with **functional interpretation**, aiding in **biomarker discovery** and **therapeutic development**.

## Literature Review
The literature review covers studies and patents relevant to GCNA and COVID-19:
1. **Weighted Gene Co-Expression Network Analysis for COVID-19 and Stroke** (Cen et al., 2021): Identified **ITGA2B** and **ITGB3** as key genes linking COVID-19 to stroke via **platelet activation** and **integrin signaling**.
2. **Differential Co-Expression Network Analysis for COVID-19** (Hasankhani et al., 2021): Highlighted **hub genes** with high connectivity as potential **therapeutic targets** using **RNA-seq** and **PPI networks**.
3. **WGCNA and Machine Learning for SARS-CoV-2** (Karami et al., 2021): Identified **PGLYRP4** and **HEPHL1** in **immune response pathways** and proposed **17 FDA-approved drugs** for repurposing.
4. **Drug Repurposing for SARS-CoV-2** (MotieGhader et al., 2021): Proposed drugs like **FLUOROURACIL** and **CISPLATIN** based on **TF-miRNA-Target Gene networks**.
5. **Ethical Papers**: Discussed **informed consent**, **data privacy**, **health disparities**, and **big data ethics** in genomic research.
6. **Patents**: Included methods for **lung cancer diagnosis** (US20200191791A1) and **colon cancer recurrence** (CN109872772B) using GCNA.

## Ethics Related to Gene Expression Analysis
The project examines ethical challenges in genomic research:
- **Informed Consent**: Broad consent is needed for future unspecified studies, but it risks **autonomy**. Dynamic consent platforms are proposed.
- **Data Privacy**: Genomic data is **identifiable** and requires robust **privacy protections** to prevent discrimination.
- **Health Disparities**: Over 80% of genomic data comes from **European descent**, leading to biased datasets. **Diverse representation** and **community engagement** are critical.
- **Big Data Ethics**: Issues include **unconsented data use**, **algorithmic bias**, and **re-identification risks**, necessitating **transparency** and **explainable AI**.

## Methodology - Gene Co-expression Analysis
The methodology involves:
1. **Dataset Collection**: RNA-seq data from **GSE151418** (17 COVID-19 patients, 17 healthy controls) focusing on **PBMCs** and **immune responses**.
2. **Data Preprocessing**: Using **R** with **GEOquery**, **WGCNA**, and **DESeq2** for quality control, **outlier removal** (e.g., GSM4615000), and **normalization**.
3. **Network Construction**: **WGCNA** with a **soft-threshold power of 18** to build a **scale-free network**, identifying **turquoise** and **green modules**.
4. **Module-Trait Correlation**: **Heatmaps** and **p-values** link modules to **COVID-19 severity**.
5. **Hub Gene Identification**: Genes with high **module membership** (e.g., **H2BC9**, **IGHV2026**) are prioritized.
6. **Module Extraction**: **Turquoise** and **green** module genes exported for enrichment analysis.

## Methodology - Enrichment Analysis
Enrichment analysis is performed in **Python** using:
- **Gene Annotation**: **MyGene.info** API maps **Ensembl IDs** to **gene symbols**, **Entrez IDs**, and **GO terms**.
- **Enrichment Analysis**: **Enrichr** queries **GO Biological Process 2023**, **KEGG 2021 Human**, **Reactome 2022**, and **WikiPathways 2019 Human** to identify enriched **biological processes** (e.g., **nucleosome assembly**, **immune response**) and **pathways** (e.g., **neutrophil extracellular trap formation**).
- **Visualization**: **Dotplots** and **heatmaps** display significant terms (FDR < 0.05).

## Results
Key findings include:
- **Cluster Dendrogram**: Identified outliers (e.g., **GSM4614993**) for removal.
- **PCA**: Confirmed sample variability and **outlier exclusion**.
- **Module Significance**: **Turquoise** module (largest, immune-focused) and **green** module (DNA organization, blood clotting) are highly correlated with **COVID-19 severity**.
- **Enrichment Results**:
  - **Green Module**: Enriched in **nucleosome assembly** (**H2BC9**), **megakaryocyte development** (**GP9**), and **immune response**.
  - **Turquoise Module**: Enriched in **immune response** (**IGHV2026**), **neutrophil extracellular trap formation**, and **systemic lupus erythematosus**.
- **Visualizations**: **Dotplots** and **heatmaps** highlight significant pathways (FDR < 0.05).

## Discussion
The analysis reveals:
- **Green Module**: Genes like **H2BC9** and **ITGB3** are linked to **DNA organization** and **blood clotting**, critical for **COVID-19 coagulopathy**.
- **Turquoise Module**: **Immunoglobulin genes** (**IGHV2026**) drive **immune responses**, suggesting roles in **cytokine storms**.
- **Ethical Implications**: Findings underscore the need for **equitable data sharing** and **privacy safeguards** in genomic research.
- **Clinical Relevance**: Identified genes are potential **biomarkers** or **therapeutic targets** for **COVID-19**.

## Conclusion
The project successfully identified **turquoise** and **green** gene modules as key players in **COVID-19 pathogenesis**, with **enrichment analysis** linking them to **immune**, **clotting**, and **DNA organization** pathways. Ethical considerations and **intellectual property rights** were thoroughly explored, emphasizing **transparent governance** and **equitable access** to genomic innovations.

## Acknowledgment
Gratitude is extended to **Dr. Soman K.P. (Dean)**, **Mrs. Reshma Sanal**, and **Dr. Neelesh Ashok** for their guidance, and to **Amrita Vishwa Vidyapeetham** for providing resources.

## Requirements
- **R**: Packages: **WGCNA**, **DESeq2**, **GEOquery**, **tidyverse**, **ggplot2**, **corrplot**.
- **Python**: Packages: **gseapy**, **mygene**, **pandas**, **matplotlib**, **openpyxl**.
- **Dataset**: **GSE151418** from NCBI GEO.

## Usage
1. **Gene Co-expression Analysis**:
   - Load **GSE151418** data in **R**.
   - Run preprocessing and **WGCNA** scripts to construct networks and identify modules.
   - Export **turquoise** and **green** module genes to CSV.
2. **Enrichment Analysis**:
   - Use **Python** scripts to annotate genes with **MyGene.info**.
   - Perform **Enrichr** analysis and generate **dotplots** for visualization.
3. **Outputs**:
   - CSV files for annotated genes and enrichment results.
   - Visualizations (dendrograms, heatmaps, dotplots) in the output directory.
