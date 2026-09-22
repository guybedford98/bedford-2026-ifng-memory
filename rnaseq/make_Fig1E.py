"""Fig 1E: per-gene Rechallenged vs Stimulated mean TPM, linear axes 1-2500.

Points are coloured by where they fall relative to the log2FC +/-0.5 guides (slopes
1.414 and 0.7072); the green dashed line is Rechallenged = Stimulated. Labelled genes
come from the three lists in data/fig1e_labels.
"""

import sys
from pathlib import Path

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "shared"))
import config as cfg
import figstyle as fs

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8")

LABELS = Path(__file__).resolve().parent / "data" / "fig1e_labels"

# Set the font to Arial for all elements
plt.rcParams['font.family'] = 'Arial'

tpms_df = pd.read_csv(cfg.path("rnaseq", "scatter_mean_tpm"))

genes_to_label_df_1 = pd.read_csv(LABELS / "Genes to label Up 250+.csv")
genes_to_label_df_2 = pd.read_csv(LABELS / "Genes to label nonmem 250+.csv")
genes_to_label_df_3 = pd.read_csv(LABELS / "Genes to label DOWN 250+.csv")

fig = plt.figure(figsize=(2.75, 2.75))  # Adjust figure size to square
sns.scatterplot(x='IFNg AVG', y='Rechallenge AVG', data=tpms_df, color='grey', s=50)  # Plot all genes in grey
point_size = 50


# Function to determine the color of the point based on its position relative to the lines
def determine_color(x, y):
    if y > 1.414 * x:
        return 'red'  # Above the upper line
    elif y < 0.7072 * x:
        return 'blue'  # Below the lower line
    return 'grey'  # For all other points


# Overlay points with conditional colors and outlines
for index, row in tpms_df.iterrows():
    color = determine_color(row['IFNg AVG'], row['Rechallenge AVG'])
    # Only plot with outlines if not grey to avoid redundancy
    if color != 'grey':
        plt.scatter(row['IFNg AVG'], row['Rechallenge AVG'], color=color, edgecolor='white', linewidth=0.5, s=point_size)

# Plotting and labeling genes from all three groups
offset = 0.03
for df, color, ha in [(genes_to_label_df_1, 'red', 'right'),
                      (genes_to_label_df_2, '#cc00c9', 'right'),
                      (genes_to_label_df_3, 'blue', 'left')]:  # Third group in blue
    for gene in df['Gene']:
        gene_data = tpms_df[tpms_df['Gene'] == gene]
        if not gene_data.empty:
            plt.scatter(gene_data['IFNg AVG'], gene_data['Rechallenge AVG'], color=color, edgecolors='white', s=point_size, zorder=3)
            plt.text(float(gene_data['IFNg AVG'].iloc[0]) * (1 + offset),
                     float(gene_data['Rechallenge AVG'].iloc[0]) * (1 + offset),
                     gene, fontsize=9, ha=ha, color='black')

plt.xlim(1, 2500)
plt.ylim(1, 2500)

plt.plot([1, 2500], [1.414, 2500 * 1.414], color='purple', linestyle=':', linewidth=3, label='log2(FC)>0.5')
plt.plot([1, 2500], [1, 2500], color='green', linestyle='--', linewidth=3, label='Trained+IFNg=Naive+IFNg')
plt.plot([1, 2500], [0.7072, (0.7072) * 2500], color='purple', linestyle=':', linewidth=3, label='log2(FC)<-0.5')

plt.gca().set_aspect('equal', adjustable='box')

plt.xlabel('Naive+IFNg (TPM)', fontsize=9)
plt.ylabel('Trained+IFNg (TPM)', fontsize=9)

plt.xticks(fontsize=9)
plt.yticks(fontsize=9)

# Adjust the legend to be outside and to the right of the graph
plt.legend(loc='upper left', bbox_to_anchor=(1, 1), fontsize='9')

fs.save(fig, cfg.out_dir("Fig1E") / "Fig1E")
