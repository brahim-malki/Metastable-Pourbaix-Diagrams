"""
Created on Sun Oct  6 12:15:43 2019

@author: brahim
"""

# from pymatgen import MPRester modified in this pymatgen version as:
from pymatgen.ext.matproj import MPRester    

from pymatgen.analysis.pourbaix_diagram import PourbaixDiagram, generate_entry_label, PourbaixPlotter
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import seaborn as sb
import numpy as np

    
from random import seed
from random import random

import pandas as pd


# Initialize the MP Rester

#if __name__ == "__main__":
#    MAPI_KEY = ""  # You must change this to your Materials API key! (or set MAPI_KEY env variable)  
#mpr = MPRester(MAPI_KEY)  # object for connecting to MP Rest interface

# use new api key
if __name__ == "__main__":
    MAPI_KEY = ""  # You must change this to your Materials API key! (or set MAPI_KEY env variable)  
mpr = MPRester(MAPI_KEY)  # object for connecting to MP Rest interface




# Get all pourbaix entries corresponding to the Fe-Cr-water chemical system.
entries = mpr.get_pourbaix_entries(["Fe", "Cr"])
concentration = 1e-6

# with fixed concentration
pbx = PourbaixDiagram(entries, comp_dict={"Fe": 0.83, "Cr": 0.17  }, conc_dict={"Fe": concentration, "Cr": concentration})
#plotter = PourbaixPlotter(pbx)


# return the list of stable and unstable entries
stable_entries = pbx.stable_entries
unstable_entries = pbx.unstable_entries

# define (pH, phi) limits 
limits = [[-0, 14],[-1.1, 1.1]]


stable_domain   = pbx.get_pourbaix_domains(stable_entries, limits)
list_domains = list(stable_domain[1].keys())

vertices_list=list(stable_domain[1].values())

# Add coloring.
#colors = sb.color_palette("Set2", len(stable_entries))
colors2D = sb.color_palette("Blues", len(list_domains))


# Plotting details...
font = "Times new roman"
fig = plt.figure(figsize=(14, 9))
ax1 = fig.gca()
ax1.set_xlim([-0, 14])
ax1.set_xticklabels([int(t) for t in ax1.get_xticks()], fontname=font, fontsize=38)
ax1.set_ylim(-1.1, 1.1)
yticks = np.arange(-1.1, 1.1, 0.5)
ax1.set_yticks(yticks)

#ax1.set_yticklabels(ax1.get_yticks(), fontname=font, fontsize=38)
ax1.set_yticklabels([round(t, 1) for t in ax1.get_yticks()], fontname=font, fontsize=38)

#ax1.set_xlabel("pH", fontname=font, fontsize=18)
#ax1.set_ylabel("Potential vs. SHE (V)", fontname=font, fontsize=18)

# Outline water's stability range.
ax1.plot([-2, 16], [0, -0.829], color="gray", linestyle="--", alpha=0.7, linewidth=2)
ax1.plot([-2, 16], [1.229, 0.401], color="gray", linestyle="--", alpha=0.7, linewidth=2)


i = 0
j = 0
for vertices in vertices_list:
    label =  generate_entry_label(list_domains[j])
    j += 1
    col = colors2D[i]
    i += 1
    center_x = sum([v[0] for v in vertices])/len(vertices)
    center_y = sum([v[1] for v in vertices])/len(vertices)
    patch = Polygon(vertices, closed=True, fill=True, facecolor=col, linewidth=1, edgecolor="black")
#    ax1.text(center_x, center_y, label, verticalalignment="center", horizontalalignment="center", fontname=font, fontsize=11)
    ax1.add_patch(patch)

# Display plot
plt.subplots_adjust(left=0.08, right=0.95, top=0.95, bottom=0.08)
plt.savefig('Figure-c=10-6-Fe17Cr-stable-PD-slabel.png', dpi=800)
plt.show()

