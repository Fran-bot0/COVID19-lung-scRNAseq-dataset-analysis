import scanpy as sc
import scvi
import pandas as pd
import numpy as np
from scipy.sparse import csr_matrix


def create_adata(path):
    '''Creates and Anndata object and annotates the sample name.'''

    # Create AnnData object.
    new_adata = sc.read_csv(path).T

    
    new_adata.obs['Sample'] = path[28:-15]
    
    # Converting the expression matrix of the AnnData object into a sparse matrix format.
    new_adata.X = csr_matrix(new_adata.X)

    return new_adata


def filter_genes(unfilt_adata, min_cell_num=10, top_var_genes=2000):
    '''Filters genes from the input AnnData object based on minimum cell presence and variability.'''

    out_adata = unfilt_adata.copy()
    
    # Only keep genes which are found in atleast 10 of the cells.
    sc.pp.filter_genes(out_adata, min_cells = min_cell_num)

    # Only keep the 2000 most variable genes.
    sc.pp.highly_variable_genes(out_adata, n_top_genes = top_var_genes, subset = True, flavor = 'seurat_v3')
    
    return out_adata


def identify_doublets(filt_adata, score_threshold=1):
    '''Trains a Variational Autoencoder (VAE) model, and applies the SOLO tool to classify cells as doublets or singlets.
    It calculates the score difference between these predictions and filters out low-certainty doublet predictions based on the provided threshold.'''
    
    # Preparing the AnnData object for use with scVI.
    scvi.model.SCVI.setup_anndata(filt_adata)

    # Initializing a Variational Autoencoder (VAE) model using the prepared data.
    vae = scvi.model.SCVI(filt_adata)
    
    # Training the VAE on the RNA-seq data.
    vae.train()

    # SOLO is a specialized tool built on scVI for identifying doublets in single-cell RNA-seq data.
    solo = scvi.external.SOLO.from_scvi_model(vae)
    
    # Trains the SOLO model for doublet classification.
    solo.train()

    # Creating a column with the prediction.
    filter_df = solo.predict()

    # Creating a dataframe with the predicted scores for doublet or singlet.
    filter_df['prediction'] = solo.predict(soft = False)

    # Creating a column in which the value is the difference in the scores.
    filter_df['score_difference'] = filter_df['doublet'] - filter_df['singlet']

    # Eliminating the cells predicted as doublets with low scores (low certainty).
    doublets_df = filter_df[(filter_df['prediction'] == 'doublet') & (filter_df['score_difference'] > score_threshold)]

    return doublets_df


def remove_doublets(wd_adata, doublets):
    '''Removes cells identified as doublets from an AnnData object.'''
    
    # Adding a column to the new AnnData object which stores a False value if the cell is a singlet and a True value if its a doublet.
    wd_adata.obs['doublet'] = wd_adata.obs.index.isin(doublets.index)
    
    # Keep only the rows (cells) which are not doublets (labled as True).
    wd_adata = wd_adata[~wd_adata.obs.doublet]

    return wd_adata


def label_mito(m_adata):
    '''Labels mitochondrial genes in an AnnData object.'''
    
    for gene in m_adata.var.index:
        if 'mt-' in gene or 'Mt-' in gene or 'MT-' in gene:
            mito_id = gene[:3]
            break

    # Adding a new annotation to the variables (genes) which stores a Boolean value for each gene, indicating whether it is a mito gene.
    m_adata.var['mito'] = m_adata.var.index.str.startswith(mito_id)

    return m_adata


def label_ribo(r_adata, path):
    '''Labels ribosomal genes in an AnnData object.'''

    # Getting a dataframe with the ribosomal genes.
    ribo_genes = pd.read_table(path, skiprows=2, header = None)

    # Adding a column to the AnnData object which stores a True value if the gene is a ribosomal gene and False if not.
    r_adata.var['ribo'] = r_adata.var_names.isin(ribo_genes[0].values)

    return r_adata


def quality_control(qc_adata, min_gene_num=200, qt=0.98, max_mito_pct_count=10, max_ribo_pct_count=30):
    '''Performs quality control on an AnnData object by filtering low-quality cells based on minimum gene count per cell, hihg gene count outliers,
    mitochondrial RNA content, and ribosomal RNA content.'''
    
    # Calculate the Quality Control (QC) metrics
    sc.pp.calculate_qc_metrics(qc_adata, qc_vars=['mito', 'ribo'], percent_top=None, log1p=False, inplace=True)

    # Remove the cells with less than a 200 gene count.
    sc.pp.filter_cells(qc_adata, min_genes=min_gene_num)

    # Determining the value below which 99% of the n_genes_by_counts are.
    upper_lim = np.quantile(qc_adata.obs.n_genes_by_counts.values, qt)
    
    # Removing n_genes_by_counts outliers.
    qc_adata = qc_adata[qc_adata.obs.n_genes_by_counts < upper_lim]
    
    # Removing cells with more than 10% mtRNA.
    qc_adata = qc_adata[qc_adata.obs.pct_counts_mito < max_mito_pct_count]
    
    # Removing cells with more than 30% rbRNA.
    qc_adata = qc_adata[qc_adata.obs.pct_counts_ribo < max_ribo_pct_count]

    return qc_adata


def normalization(norm_adata):
    '''Normalizes the raw counts and converts them to log-transformed values.'''
    
    # Normalize the counts
    sc.pp.normalize_total(norm_adata, target_sum=1e4)
    
    # Convert the normalized counts to log counts
    sc.pp.log1p(norm_adata)

    # Assigning the current state of the AnnData object to its .raw attribute.
    norm_adata.raw = norm_adata

    return norm_adata

