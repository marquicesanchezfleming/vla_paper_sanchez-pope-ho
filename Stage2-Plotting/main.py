import pandas as pd
import matplotlib.pyplot as plt
import re
from matplotlib.lines import Line2D

df = pd.read_csv('/Users/Djslime07/RecreatingStroh/25A-386 Project Data - Light Curve Data.csv')

classification_styles = {
    'archival_none': {'color': '#660099', 'marker': 'o'},
    'VLASS_none': {'color': "#660099", 'marker': 's'},
    'archival_poor': {'color': '#04e3ff', 'marker': 'o'},
    'VLASS_poor': {'color': '#04e3ff', 'marker': 's'},
    'new_VLASS': {'color': '#fb19dd', 'marker': '*'},
}

label_these = ['SN2024ehs', 'SN2019pqo', 'SN2021bmf', 'SN2021tcv', 'SN2022cb', 
               'PTF11qcj', 'SN2004c', 'SN2003bg', 'SN2004dk', 'SN2016coi', '1986J', 'SN2012au']

fig, ax = plt.subplots(figsize=(20, 12))
texts = []

for sn, group in df.groupby('SN'):
    group_sorted = group.sort_values(by='dt (years)')
    group_sorted = group_sorted[group_sorted['dt (years)'] >= 0] 
    if group_sorted.empty:
        continue

    line_cls = str(group_sorted.iloc[0]['classification'])
    line_style = classification_styles.get(line_cls, {'color': 'gray', 'marker': 'o'})

    if sn in label_these:
        ax.plot(group_sorted['dt (years)'], group_sorted['Luminosity'],
                linestyle='-', color='black', linewidth=4.5, alpha=0.6, zorder=0)

    ax.plot(group_sorted['dt (years)'], group_sorted['Luminosity'],
            linestyle='-', color=line_style['color'], linewidth=1.8, alpha=0.8, zorder=1)

    for _, row in group_sorted.iterrows():
        cls = str(row['classification'])
        status = row['status']
        style = classification_styles.get(cls, {'color': 'gray', 'marker': 'o'})
        marker = style['marker'] if status == 'real' else 'v'

        ax.scatter(row['dt (years)'], row['Luminosity'],
                   color=style['color'], marker=marker,
                   s=100, edgecolor='black', linewidth=0.8, zorder=2)

    if sn in label_these:
        final_row = group_sorted.iloc[-1]
        short_label = re.sub(r'^SN20(\d{2})', r'\1', sn)
        x_offset = final_row['dt (years)'] * 1.04
        y_coord = final_row['Luminosity']

        text = ax.text(x_offset, y_coord, short_label,
                       fontsize=14, color=line_style['color'],
                       ha='left', va='center')
        texts.append(text)

legend_elements = [
    Line2D([0], [0], marker='o', color='#660099', label='Archival Observations', markersize=10, linestyle='None'),      #none
    Line2D([0], [0], marker='s', color='#660099', label='VLASS Observations', markersize=10, linestyle='None'),         #none
    Line2D([0], [0], marker='o', color='#04e3ff', label='Archival Observations', markersize=10, linestyle='None'),      #poor
    Line2D([0], [0], marker='s', color='#04e3ff', label='VLASS Observations', markersize=10, linestyle='None'),         #poor
    Line2D([0], [0], marker='*', color="#fb19dd", label='VLASS 25A-386 Observations', markersize=14, linestyle='None'), #new
]

ax.legend(handles=legend_elements, loc='lower left', fontsize=14, frameon=True)

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel('Time Since Explosion (years)', fontsize=18)
ax.set_ylabel('Luminosity (erg/s/Hz)', fontsize=18)
plt.tick_params(axis='both', which='major', labelsize=20)
plt.tight_layout()
plt.show()