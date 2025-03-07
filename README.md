# Analysis of a Dataset About the Effects of COVID19 on Lung Tissue
This Jupyter Notebook presents an analysis of an scRNA-seq dataset exploring the effects of COVID-19 on lung tissue. It covers preprocessing, integration, cell annotation, clustering, data analysis, and visualization.
### Analysis Overview:
- **Preprocessing**: Doublet removal, mitochondrial and ribosomal gene labeling, quality control (QC), and normalization.
- **Integration**: Performed using scVI.
- **Clustering**: Applied the Leiden algorithm.
- **Cell Annotation**: Manually curated using PanglaoDB and DAVID.
- **Data Analysis & Visualization**: Conducted with Scanpy.

# Requirements:
- matplotlib
- numpy
- pandas
- seaborn
- scanpy
- scipy
- scVI
