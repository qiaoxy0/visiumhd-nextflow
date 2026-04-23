import os
os.environ["OPENCV_IO_MAX_IMAGE_PIXELS"] = str(2**40)

import logging
import matplotlib.pyplot as plt
import scanpy as sc
import numpy as np
import pandas as pd
import tensorflow as tf
from tensorflow.python.client import device_lib
import bin2cell as b2c
import pySTIM as pst
from skimage.measure import regionprops
import scipy
from scipy.io import mmwrite
import seaborn as sns
import anndata as ad
from copy import deepcopy
from tifffile import tifffile
import skimage.segmentation


# Setup logging
logging.basicConfig(filename='pipeline.log', level=logging.INFO, format='%(asctime)s - %(message)s')


def main(path, source_image_path, spaceranger_image_path, mpp, out_dir):
    """
    Main function to process and analyze spatial transcriptomic data.
    
    Parameters:
    path (str): Path to the directory containing the binned output files.
    source_image_path (str): Path to the source image file.
    spaceranger_image_path (str): Path to the spaceranger image directory.
    mpp (float): Microns per pixel (mpp) for scaling.
    out_dir (str): Output directory to save the results.
    """
    logging.info("Starting main function")
    os.makedirs(out_dir, exist_ok=True)
    
    # Load and preprocess spatial transcriptomic data
    logging.info("Loading spatial transcriptomic data")
    adata = b2c.read_visium(path, 
                            source_image_path=source_image_path, 
                            spaceranger_image_path=spaceranger_image_path)
    adata.var_names_make_unique()
    sc.pp.filter_genes(adata, min_cells=3)
    sc.pp.filter_cells(adata, min_counts=1)
    logging.info(f"Shape of adata after filtering: {adata.shape}")

    # Generate and save scaled H&E image
    logging.info("Generating and saving scaled H&E image")
    b2c.scaled_he_image(adata, mpp=mpp, save_path=f"{out_dir}/HE.tiff", crop=False)
    
    # Run Stardist segmentation to generate label masks
    logging.info("Running Stardist segmentation")
    b2c.stardist(image_path=f"{out_dir}/HE.tiff", 
                 labels_npz_path=f"{out_dir}/HE.npz", 
                 stardist_model="2D_versatile_he", 
                 prob_thresh=0.1)
    
    border_color=[255, 0, 0]

    image_path = f"{out_dir}/HE.tiff"
    labels_path = f"{out_dir}/HE.npz"
    img = tifffile.imread(image_path)
    labels = scipy.sparse.load_npz(labels_path)

    labels_sparse = scipy.sparse.coo_matrix(labels)
    border_sparse = scipy.sparse.coo_matrix(
                skimage.segmentation.find_boundaries(
                    np.array(labels_sparse.todense())
                )
            )
    
    img[border_sparse.row, border_sparse.col, :] = border_color

    fig, ax = plt.subplots(figsize = (20,20), dpi=100)
    ax.imshow(img)
    ax.axis('off')
    plt.tight_layout()
    plt.savefig(f'{out_dir}/segmentation_whole.png', dpi=300, bbox_inches='tight')
    plt.close(fig)
    
    # Insert label data into AnnData object
    logging.info("Inserting label data into AnnData object")
    b2c.insert_labels(adata, 
                      labels_npz_path=f"{out_dir}/HE.npz", 
                      basis="spatial", 
                      spatial_key="spatial", 
                      mpp=mpp, 
                      labels_key="labels_he")
    
    # Expand the label annotations
    logging.info("Expanding label annotations")
    b2c.expand_labels(adata, 
                      labels_key='labels_he', 
                      expanded_labels_key="labels_he_expanded", 
                      max_bin_distance=2)
    
    logging.info(f"Non-zero labels in 'labels_he': {len(adata.obs['labels_he'][adata.obs['labels_he'] != 0])}")
    logging.info(f"Non-zero labels in 'labels_he_expanded': {len(adata.obs['labels_he_expanded'][adata.obs['labels_he_expanded'] != 0])}")
    
    labels_path = adata.uns["bin2cell"]["labels_npz_paths"].get('labels_he', None)
    if labels_path:
        logging.info("Processing nuclei centroids")
        labels_sparse = scipy.sparse.load_npz(labels_path)
        nuclei_centroids = extract_nuclei_centroids(labels_sparse)
        adata.uns["nuclei_centroid"] = nuclei_centroids
        
        labels_key = "labels_he_expanded"
        spatial_keys = ["spatial"]
        
        adata = adata[adata.obs[labels_key] != 0]
        logging.info(f"Filtered adata to non-zero labels, new shape: {adata.shape}")
        
        cell_to_bin = pd.get_dummies(adata.obs[labels_key], sparse=True).sparse.to_coo().tocsr().T
        cell_names = [str(i) for i in range(cell_to_bin.shape[0])]
        
        X = cell_to_bin.dot(adata.X).tocsr()
        cell_adata = ad.AnnData(X, var=adata.var)
        cell_adata.obs_names = cell_names
        cell_adata.obs['object_id'] = [int(i) for i in cell_names]
        
        bin_count = np.asarray(cell_to_bin.sum(axis=1)).flatten()
        row_means = scipy.sparse.diags(1 / np.maximum(bin_count, 1))
        cell_adata.obs['bin_count'] = bin_count
        
        cell_adata.obs['nuclei_row'] = [nuclei_centroids.get(int(cell_name), (np.nan, np.nan))[0] for cell_name in cell_names]
        cell_adata.obs['nuclei_col'] = [nuclei_centroids.get(int(cell_name), (np.nan, np.nan))[1] for cell_name in cell_names]
        
        for spatial_key in spatial_keys:
            cell_adata.obsm[spatial_key] = row_means.dot(cell_to_bin).dot(adata.obsm[spatial_key])
        
        diameter_scale_factor = np.sqrt(np.mean(bin_count))
        library = list(adata.uns['spatial'].keys())[0]
        
        cell_adata.uns['spatial'] = deepcopy(adata.uns['spatial'])
        cell_adata.uns['spatial'][library]['scalefactors']['spot_diameter_fullres'] *= diameter_scale_factor
        
        visualize_qc_metrics(cell_adata, out_dir)
        
        cell_adata.write_h5ad(f"{out_dir}/cell_adata.h5ad")

        cell_adata.obsm['spatial'] = cell_adata.obs[['nuclei_col', 'nuclei_row']].values
        
        mmwrite(f"{out_dir}/adata_count.mtx", cell_adata.X)
        pd.DataFrame(cell_adata.obs_names, columns=["cells"]).to_csv(f"{out_dir}/adata_cells.csv")
        pd.DataFrame(cell_adata.var_names, columns=["genes"]).to_csv(f"{out_dir}/adata_genes.csv")
        cell_adata.obs.to_csv(f"{out_dir}/adata_meta.csv")

def visualize_qc_metrics(cell_adata, out_dir):
    logging.info("Generating QC metrics")
    mean_bin_count = cell_adata.obs["bin_count"].mean()
    
    fig, ax = plt.subplots(dpi = 150, figsize = (4, 3))
    ax = sns.histplot(cell_adata.obs["bin_count"], bins=120, kde=False)
    ax.set_xlim(0, 100)
    ax.axvline(mean_bin_count, color='red', linestyle='--', linewidth=2, label=f'Mean = {mean_bin_count:.2f}')
    ax.text(mean_bin_count + 0.5, ax.get_ylim()[1] * 0.9, f'Mean = {mean_bin_count:.2f}', color='red', fontsize=10)
    ax.set_xlabel('Bin counts')
    ax.set_ylabel('Frequency')
    plt.savefig(f"{out_dir}/hist_bin.png", bbox_inches='tight', dpi=300)
    logging.info("Saved bin count histogram")
    
    sc.pp.calculate_qc_metrics(cell_adata, inplace=True)
    mean_cell_count = cell_adata.obs["total_counts"].mean()
    fig, ax = plt.subplots(dpi = 150, figsize = (4, 3))
    ax = sns.histplot(cell_adata.obs["total_counts"], bins=200, kde=False, ax = ax)
    ax.axvline(mean_cell_count, color='red', linestyle='--', linewidth=2, label=f'Mean = {mean_cell_count:.2f}')
    ax.text(mean_cell_count + 150, ax.get_ylim()[1] * 0.9, f'Mean = {mean_cell_count:.2f}', color='red', fontsize=10)
    ax.set_xlim(-20, 6000)
    plt.savefig(f"{out_dir}/total_count.png", bbox_inches = 'tight', dpi = 300)
    logging.info("Saved total count histogram")
    


def extract_nuclei_centroids(labels_sparse):
    labels_dense = labels_sparse.toarray()
    props = regionprops(labels_dense)
    centroids = {prop.label: prop.centroid for prop in props}
    return centroids


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Process spatial transcriptomic data.')
    parser.add_argument('--path', required=True, help='Path to binned output files')
    parser.add_argument('--source_image_path', required=True, help='Path to source image file')
    parser.add_argument('--spaceranger_image_path', required=True, help='Path to spaceranger image directory')
    parser.add_argument('--mpp', type=float, required=True, help='Microns per pixel (mpp) for scaling')
    parser.add_argument('--out_dir', required=True, help='Output directory to save the results')
    args = parser.parse_args()
    main(args.path, args.source_image_path, args.spaceranger_image_path, args.mpp, args.out_dir)
  