# Load modules
from dataclasses import dataclass
from unicodedata import name
import scanpy
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy.stats as pearsonr

import bulk2space
from bulk2space import Bulk2Space

model = Bulk2Space()

# Specify file paths
bulk_data = "Path/to/file.csv"
ref_sc_data = "Path/to/file.csv"
ref_sc_meta = "Path/to/file.csv"
st_meta = "Path/to/file.csv"
st_data = "Path/to/file.csv"

# Generate and train vae model
generate_sc_meta, generate_sc_data = model.train_vae_and_generate(
    input_bulk_path=bulk_data,
    input_sc_data_path=ref_sc_data,
    input_sc_meta_path=ref_sc_meta,
    input_st_data_path=st_data,
    input_st_meta_path=st_meta,
    ratio_num=1,
    top_marker_num=500,
    gpu=0,
    batch_size=512,
    learning_rate=1e-4,
    hidden_size=256,
    epoch_num=3500,
    vae_save_dir="save model",
    vae_save_name="trained_model",
    generate_save_dir="output",
    generate_save_name="results",
)

# Load the deconvoluted scRNA-Seq data
generate_sc_data = pd.read_csv(
    "output_no_doublets/SN_human_from_mouse_no_doublets_sc_data.csv")
generate_sc_meta = pd.read_csv(
    "output_no_doublets/SN_human_from_mouse_no_doublets_sc_celltype.csv")

# Calculate cell type number in deconvoluted bulk-seq data
ct_stat = pd.DataFrame(generate_sc_meta["Cell_type"].value_counts())
ct_name = list(ct_stat.index)
ct_num = list(ct_stat["Cell_type"])
color = [
    "#EF4444",
    "#14B8A6",
    "#F97316",
    "#7C3AED",
    "#0EA5E9",
    "#3B82F6",
    "#000000",
    "#10B981",
    "#EC4899",
    "#7DD3FC",
    "#F472B6",
    "#D946EF",
]

# Plot the results
plt.barh(ct_name, ct_num, color=color)
plt.yticks(ct_name, ct_name, fontsize=20)
plt.xlabel("Cell number", fontsize=18)
plt.xticks(fontsize=18)
plt.ylabel("Cell type")
plt.show()

# load reference sc data
input_data = bulk2space.utils.load_data(
    input_bulk_path=bulk_data,
    input_sc_data_path=ref_sc_data,
    input_sc_meta_path=ref_sc_meta,
    input_st_data_path=st_data,
    input_st_meta_path=st_meta,
)

# Calculate 200 marker genes of each cell type
sc = scanpy.AnnData(input_data["input_sc_data"].T)
sc.obs = input_data["input_sc_meta"][["Cell_type"]]
scanpy.tl.rank_genes_groups(sc, "Cell_type", method="wilcoxon")
marker_df = pd.DataFrame(sc.uns["rank_genes_groups"]["names"]).head(200)
marker = list(np.unique(np.ravel(np.array(marker_df))))

# the mean expression of 200 marker genes of input sc data
sc_marker = input_data["input_sc_data"].loc[marker, :].T
sc_marker["Cell_type"] = input_data["input_sc_meta"]["Cell_type"]
sc_marker_mean = sc_marker.groupby("Cell_type")[marker].mean()

# the mean expression of 200 marker genes of deconvoluted bulk-seq data
generate_sc_meta.index = list(generate_sc_meta["Cell"])
generate_sc_data_new = generate_sc_data.T
new_header = generate_sc_data_new.iloc[0]
generate_sc_data_new = generate_sc_data_new[1:]
generate_sc_data_new.columns = new_header
# ---------------------------------------------------------
# Added this line to fix code breaking at groupby
generate_sc_data_new = generate_sc_data_new.astype(float)
# ---------------------------------------------------------
generate_sc_data_new["Cell_type"] = generate_sc_meta["Cell_type"]
generate_sc_marker_mean = generate_sc_data_new.groupby(["Cell_type"])[
    marker].mean()

intersect_cell = list(
    set(sc_marker_mean.index).intersection(set(generate_sc_marker_mean.index))
)

generate_sc_marker_mean = generate_sc_marker_mean.loc[intersect_cell]
sc_marker_mean = sc_marker_mean.loc[intersect_cell]

# Calculate expression correlation
sc_marker_mean = sc_marker_mean.T
generate_sc_marker_mean = generate_sc_marker_mean.T
coeffmat = np.zeros(
    (sc_marker_mean.shape[1], generate_sc_marker_mean.shape[1]))
for i in range(sc_marker_mean.shape[1]):
    for j in range(generate_sc_marker_mean.shape[1]):
        corrtest = pearsonr.pearsonr(
            sc_marker_mean[sc_marker_mean.columns[i]],
            generate_sc_marker_mean[generate_sc_marker_mean.columns[j]],
        )
        coeffmat[i, j] = corrtest[0]

rf_ct = list(sc_marker_mean.columns)
generate_ct = list(generate_sc_marker_mean.columns)

# Plot the results
fig, ax = plt.subplots()
im = ax.imshow(coeffmat, cmap="RdBu_r")
ax.set_xticks(np.arange(len(rf_ct)))
ax.set_xticklabels(rf_ct, fontsize=14)
ax.set_yticks(np.arange(len(generate_ct)))
ax.set_yticklabels(generate_ct, fontsize=14)
plt.xlabel("scRNA-seq reference", fontweight='bold', fontsize=15)
plt.ylabel("Deconvoluted bulk data", fontweight='bold', fontsize=15)
plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor")
plt.colorbar(im)
ax.set_title("Expression correlation", fontsize=15)
fig.tight_layout(pad=0.1)
plt.show()
fig.savefig('Correlation matrix.png', dpi=600)
