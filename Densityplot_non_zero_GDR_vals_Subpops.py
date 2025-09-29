import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

# Load your data
file_path = '/home/abera/2024/All_GDRs/updated_Sub_gdr.tsv'
df = pd.read_csv(file_path, sep='\t', index_col=0).reset_index()

# Assign the first column to 'GeneName'
df.columns = ['GeneName'] + list(df.columns[1:])

# Simplify gene names
def simplify_gene_name(gene_str):
    if isinstance(gene_str, str) and gene_str.count('___') >= 2:
        parts = gene_str.split('___')
        return '___'.join(parts[:3])
    return gene_str

df['GeneName'] = df['GeneName'].apply(simplify_gene_name)

# Reshape data
long_df = df.melt(id_vars=['GeneName'], var_name='Sub', value_name='GDR')
long_df = long_df.dropna(subset=['GDR'])
long_df['GDR'] = pd.to_numeric(long_df['GDR'], errors='coerce')

#  Exclude GDR == 0 values
long_df = long_df[long_df['GDR'] > 0]

# Map sub-populations to super-populations
super_pop_map = {
    'GWD': 'AFR', 'LWK': 'AFR', 'YRI': 'AFR', 'ESN': 'AFR', 'MSL': 'AFR', 'ASW': 'AFR', 'ACB': 'AFR',
    'CEU': 'EUR', 'FIN': 'EUR', 'GBR': 'EUR', 'IBS': 'EUR', 'TSI': 'EUR',
    'CHB': 'EAS', 'JPT': 'EAS', 'CHS': 'EAS', 'CDX': 'EAS', 'KHV': 'EAS',
    'GIH': 'SAS', 'PJL': 'SAS', 'STU': 'SAS', 'ITU': 'SAS', 'BEB': 'SAS',
    'MXL': 'AMR', 'PUR': 'AMR', 'CLM': 'AMR', 'PEL': 'AMR'
}
long_df['SuperPop'] = long_df['Sub'].map(super_pop_map)

# Create PDF output
pdf_output_file = '/home/abera/2024/Results/Sub_pops/All_18174_GDR_Superpops_DensityLines.pdf'

with PdfPages(pdf_output_file) as pdf:
    sns.set(style="whitegrid")
    superpops = long_df['SuperPop'].dropna().unique()

    for superpop in superpops:
        subset = long_df[long_df['SuperPop'] == superpop]

        plt.figure(figsize=(10, 6))
        ax = sns.kdeplot(
            data=subset,
            x='GDR',
            bw_adjust=1.2,   # smoothness factor
            hue='Sub',
            common_norm=False,
            fill=False,  # no solid fill
            lw=2,  # line width
            legend=True  #  let seaborn handle legend
        )

        ax.set_xlim(long_df['GDR'].min(), long_df['GDR'].max())

        ax.set_title(f"GDR Distribution : {superpop}", fontsize=16)
        ax.set_xlabel("GDR Value", fontsize=12)
        ax.set_ylabel("Density", fontsize=12)

        plt.tight_layout()
        pdf.savefig()
        output_file = f"/home/abera/2024/Results/Sub_pops/18174_GDR_DensityLine_{superpop}.png"
        plt.savefig(output_file, dpi=600, bbox_inches='tight')
        plt.close()
