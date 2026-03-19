import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib

matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

df = pd.read_csv("Brugge_5_prime_merged.csv", index_col=0)

genes = ['CDKN2A-p16', 'CDKN2A-p14', 'CDKN2A-p14+p16', 'CDKN2A', 'dialup-p16', 'dialup-p14']

total_cells = len(df)
percentages = [np.sum(df[gene] > 0) / total_cells * 100 for gene in genes]
sems = [df.groupby('sample')[gene].apply(lambda x: (x > 0).mean()).sem() * 100 for gene in genes]

x_labels = [
    r'p16$^{INK4A}$' + '\n10x 3\'RNA-seq',
    r'p14$^{ARF}$' + '\n10x 3\'RNA-seq',
    'CDKN2A\nexon 2+3 GTF',
    'CDKN2A\nstandard GTF',
    r'p16$^{INK4A}$' + '\nenrichment',
    r'p14$^{ARF}$' + '\nenrichment',
]

colors = ['blue', 'green', 'gray', 'dimgray', 'blue', 'green']

fig, ax = plt.subplots(figsize=(8, 6))
bars = ax.bar(range(len(genes)), percentages, yerr=sems, capsize=5, color=colors, width=0.6)

# Add bracket connecting the two CDKN2A bars to highlight they are the same
bracket_y = max(percentages[2], percentages[3]) + max(sems[2], sems[3]) + 1.5
ax.plot([2, 2, 3, 3], [bracket_y - 0.3, bracket_y, bracket_y, bracket_y - 0.3], color='black', linewidth=1)
ax.text(2.5, bracket_y + 0.3, 'n.s.', ha='center', va='bottom', fontsize=10)

ax.set_xticks(range(len(genes)))
ax.set_xticklabels(x_labels, rotation=45, ha='right', fontsize=9)
ax.set_ylabel('Percentage of cells (%)', fontsize=12)
ax.set_ylim(0, 45)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
plt.savefig("bar_plot_CDKN2A_exp_expanded.png", dpi=300, bbox_inches='tight')
plt.savefig("bar_plot_CDKN2A_exp_expanded.pdf", dpi=300, bbox_inches='tight')
plt.show()
print("Saved bar_plot_CDKN2A_exp_expanded.png and .pdf")
