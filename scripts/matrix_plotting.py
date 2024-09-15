import scanpy as sc
import numpy as np
import argparse

def parse_args():
    """
    Parse command-line arguments.
    """
    parser = argparse.ArgumentParser(description="Process single-cell RNA-seq data.")
    parser.add_argument('--input_dir', type=str, required=True, help='Path to the 10X Genomics directory of filtered feature-barcode matrix.')
    return parser.parse_args()

def read_data(input_dir):
    """
    Read the 10X Genomics data.
    """
    ann_data = sc.read_10x_mtx(input_dir, var_names='gene_symbols', cache=True)
    #print(ann_data.obs)  # Metadata about cells
    #print(ann_data.var)  # Metadata about genes
    return ann_data

def calculate_qc_metrics(ann_data):
    """
    Calculate quality control metrics for the data.
    """
    # Total counts for each cell
    ann_data.obs['total_counts'] = ann_data.X.sum(axis=1)
    # Number of genes expressed in each cell
    ann_data.obs['n_genes_by_counts'] = (ann_data.X > 0).sum(axis=1)
    
    # Compute mitochondrial gene percentage
    mito_genes = ann_data.var_names.str.startswith('MT-')
    ann_data.obs['mito_percent'] = np.sum(ann_data[:, mito_genes].X, axis=1) / np.sum(ann_data.X, axis=1)
    
    # Check that these metrics have been added
    #print(ann_data.obs.head())

def normalize_data(ann_data):
    """
    Normalize the data.
    """
    sc.pp.normalize_total(ann_data, target_sum=1e4)
    sc.pp.log1p(ann_data)

def perform_dimensionality_reduction_and_clustering(ann_data):
    """
    Perform PCA, compute neighbors, run UMAP and Leiden clustering.
    """
    sc.pp.pca(ann_data)
    sc.pp.neighbors(ann_data)
    sc.tl.umap(ann_data)
    sc.tl.leiden(ann_data, flavor="igraph", n_iterations=2)

def visualize_data(ann_data):
    """
    Generate and save PCA and UMAP plots.
    """
    # Save PCA plot
    sc.pl.pca(ann_data, color='total_counts', save='.png')
    # Save UMAP plots
    sc.pl.umap(ann_data, color=['total_counts', 'n_genes_by_counts', 'mito_percent'], save='.png')
    sc.pl.umap(ann_data, color='leiden', save='.pdf')

def main():

    args = parse_args()

    # Set custom directory where images will be saved
    sc.settings.figdir = f'{args.input_dir.split('/')[0]}/analysis/plots' # args.input_dir.split('/')[0] = sample_name (args.input_dir is defined in Snakefile)
    
    # 1. Read Data
    ann_data = read_data(args.input_dir)
    
    # 2. Calculate Quality Control Metrics
    calculate_qc_metrics(ann_data)
    
    # 3. Normalize Data
    normalize_data(ann_data)
    
    # 4. Dimensionality Reduction and Clustering
    perform_dimensionality_reduction_and_clustering(ann_data)
    
    # 5. Visualization
    visualize_data(ann_data)

if __name__ == "__main__":
    main()
