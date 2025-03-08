# Gene Co-Expression Network Analysis (GCNA)

This project explores gene-gene interactions in disease conditions using Gene Co-Expression Network Analysis (GCNA). We focus on analyzing gene expression data to identify gene modules, hub genes, and perform enrichment analysis.

## Project Overview
- **Dataset**: NCBI GEO (GSE152418) - Systems biological assessment of immunity to severe and mild COVID-19 infection.
- **Objective**:
  - Identify gene clusters and hub genes in COVID-19-related datasets.
  - Perform enrichment analysis on identified modules.
  - Compare co-expression networks with similar diseases.
  - Explore ethical concerns, intellectual property rights (IPR), and patent-related issues.

## Methodology
1. **Data Collection & Preprocessing**:
   - Downloaded gene expression data from NCBI GEO.
   - Merged expression data with metadata.
   - Performed quality control using hierarchical clustering, PCA, and outlier detection.
   - Normalized data using DESeq2.

2. **Network Construction**:
   - Built a weighted gene co-expression network (WGCNA).
   - Identified gene modules and correlated them with disease traits.

3. **Analysis & Results**:
   - Generated cluster dendrograms for hierarchical clustering.
   - Visualized gene module correlations using heatmaps.
   - Selected an appropriate soft threshold for network construction.

4. **Ethics & Intellectual Property Considerations**:
   - Discussed legal and ethical concerns regarding genetic data usage.
   - Reviewed patent-related issues in gene expression analysis.

## Key Findings
- Gene modules and hub genes associated with COVID-19 were identified.
- Ethical concerns like genetic data privacy and accessibility were analyzed.
- Intellectual property rights on gene sequencing technologies were explored.

## References
- Weighted Gene Co-Expression Network Analysis (2021, 2022)
- Drug Repurposing for COVID-19 using GCNA (2021)
- Intellectual Property Rights in Genetic Research (2000)
- World Intellectual Property Organization 
---

### Contributors
- **Team 11**: Minoti K, Nivedhitha S, Jaison John Samuel J, Sanchita B  
- **Faculty In-Charge**: Kelath Murali Manoj, Mrs. Reshma Sanal, Neelesh Ashok
