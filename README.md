# Analysis of a Dataset About the Effects of COVID19 on Lung Tissue
This Jupyter Notebook presents an analysis of an scRNA-seq dataset exploring the effects of COVID-19 on lung tissue. It covers preprocessing, integration, cell annotation, clustering, data analysis, and visualization.
The data used in this analysis comes from [this study](https://www.nature.com/articles/s41586-021-03569-1) and is available [here](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE171524).
### Analysis Overview:
- **Preprocessing**: Doublet removal, mitochondrial and ribosomal gene labeling, quality control (QC), and normalization.
- **Integration**: Performed using scVI.
- **Clustering**: Applied the Leiden algorithm.
- **Cell Annotation**: Manually curated using [PanglaoDB](panglaodb.se) and [DAVID](davidbioinformatics.nih.gov).
- **Data Analysis & Visualization**: Conducted with Scanpy.

# Requirements:
- matplotlib
- numpy
- pandas
- seaborn
- scanpy
- scipy
- scVI

# References
- [A molecular single-cell lung atlas of lethal COVID-19](https://www.nature.com/articles/s41586-021-03569-1).
- Dataset [GSE171524](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE171524).
- [Sambomics](https://www.youtube.com/watch?v=uvyG9yLuNSE).
- [List of Ribosome Genes](https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/KEGG_RIBOSOME.html).
- [scanpy.readthedocs.io/en/stable/tutorials/basics/clustering.html](scanpy.readthedocs.io/en/stable/tutorials/basics/clustering.html)
- [sc-best-practices.org/preprocessing_visualization/quality_control.html](sc-best-practices.org/preprocessing_visualization/quality_control.html)
- [docs.scvi-tools.org/en/stable/tutorials/notebooks/scrna/scVI_DE_worm.html](docs.scvi-tools.org/en/stable/tutorials/notebooks/scrna/scVI_DE_worm.html)
- [sc-best-practices.org/cellular_structure/annotation.html](sc-best-practices.org/cellular_structure/annotation.html)
- [panglaodb.se](panglaodb.se)
- [davidbioinformatics.nih.gov](davidbioinformatics.nih.gov)
- [scanpy.readthedocs.io/en/stable/tutorials/plotting/core.html](scanpy.readthedocs.io/en/stable/tutorials/plotting/core.html)
