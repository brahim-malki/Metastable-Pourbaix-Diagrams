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
import copy

from pymatgen.core.ion import Ion
from pymatgen.entries.computed_entries import ComputedEntry
from pymatgen.core.composition import Composition
from pymatgen.analysis.pourbaix_diagram import PourbaixDiagram, PourbaixPlotter, PourbaixEntry, ComputedEntry, IonEntry



# Initialize the MP Rester

#if __name__ == "__main__":
#    MAPI_KEY = "9TnOADgQ0sFZ6AAdt"  # You must change this to your Materials API key! (or set MAPI_KEY env variable)  
#mpr = MPRester(MAPI_KEY)  # object for connecting to MP Rest interface

# use new api key
if __name__ == "__main__":
    MAPI_KEY = "SNcD3tHhDv9Y6aOKDq0kKv7acAkmKeD8"  # You must change this to your Materials API key! (or set MAPI_KEY env variable)  
mpr = MPRester(MAPI_KEY)  # object for connecting to MP Rest interface

# Get all pourbaix entries corresponding to the Fe-Cr-water chemical system.
entries_init = mpr.get_pourbaix_entries(["Fe", "Cr"])
concentration = 1e-6


#  fonction de changement d'energie qui me manquait tant !

def modify_entries(entries_init, list_index_phase_to_modify, energy_increment):
    modified_entries = []
    for i, entry in enumerate(entries_init):
        if i in list_index_phase_to_modify:
            initial_entry = entry.entry  # Récupérer l'entrée initiale à partir de l'objet PourbaixEntry
#            print(f"Original energy: {initial_entry.energy_per_atom}")
            modified_energy = initial_entry.energy_per_atom + energy_increment
#            print(f"Modified energy: {modified_energy}")
            modified_entry = ComputedEntry(initial_entry.composition, modified_energy * initial_entry.composition.num_atoms, parameters=initial_entry.data)
            if isinstance(initial_entry, IonEntry):
                modified_entry = IonEntry(modified_entry)
            modified_entry = PourbaixEntry(modified_entry)
        else:
            modified_entry = entry

        # Skip the entry if normalization_factor would result in a division by zero
        if modified_entry.num_atoms - modified_entry.composition.get("H", 0) - modified_entry.composition.get("O", 0) == 0:
            continue

        modified_entries.append(modified_entry)

    return modified_entries


# with fixed concentration
entries_init_copy = copy.deepcopy(entries_init)

pbx = PourbaixDiagram(entries_init_copy, comp_dict={"Fe": 0.83, "Cr": 0.17  }, conc_dict={"Fe": concentration, "Cr": concentration})

# return the list of stable and unstable entries
#unstable_entries = pbx.unstable_entries
stable_entries = pbx.stable_entries

#determination des index des phases à changer : list_index
list_index_MP_id = []
list_index       = []
list_name        = []

for entry in stable_entries:
    for i, member in enumerate(entry.phase_type):        
        if "Solid" in member:
            if(entry.entry_id[i] not in list_index_MP_id):
               list_index_MP_id.append(entry.entry_id[i])
            
for id in list_index_MP_id:
    for i, entry in enumerate(entries_init_copy):
        if entry.entry_id == id:
            print(entry.entry.name)
            list_index.append(i)
            list_name.append(entry.entry.name)

#list_index_phase_to_modify = list_index
list_index_phase_to_modify = [114]  # au cas ou il faut changer qune seule phase





limits = [[-0, 14],[-1.1, 1.1]]
delta_energy = -0.2
#recalcul du diagramme de pourbaix modifié
modified_entries = modify_entries(entries_init_copy, list_index_phase_to_modify, delta_energy)
pbx = PourbaixDiagram(modified_entries, comp_dict={"Fe": 0.83, "Cr": 0.17}, conc_dict={"Fe": concentration, "Cr": concentration})
stable_entries = pbx.stable_entries

stable_domain   = pbx.get_pourbaix_domains(stable_entries, limits)
list_domains = list(stable_domain[1].keys())
vertices_list=list(stable_domain[1].values())




# construction de la liste de référence àexcuter une seule fois au début pour delta_energy = 0.0
reference_domains = []
for domain in list_domains:
    name = generate_entry_label(domain)
    if name not in reference_domains:
        reference_domains.append(name)

# Add coloring.
#colors = sb.color_palette("Set2", len(stable_entries))
color2D = sb.color_palette("Blues", len(reference_domains))
color_phase = dict(zip(reference_domains, color2D))

#color_phase["Fe(s) + Cr2O3(s)"] = "#FFA500"
color_phase["Cr2O3(s) + Fe(s)"] = "#FFA500"


   
         
# Plotting details...
font = "Times new roman"
fig = plt.figure(figsize=(14, 9))
ax1 = fig.gca()
ax1.set_xlim([-0.1, 14.1])
ax1.set_xticklabels([int(t) for t in ax1.get_xticks()], fontname=font, fontsize=38)
ax1.set_ylim(-1.2, 1.2)
ax1.set_yticklabels(ax1.get_yticks(), fontname=font, fontsize=38)
#ax1.set_xlabel("pH", fontname=font, fontsize=20)
#ax1.set_ylabel("Potential vs. SHE (V)", fontname=font, fontsize=20)

# Outline water's stability range.
ax1.plot([-2, 16], [0, -0.829], color="red", linestyle="--", alpha=0.7, linewidth=2)
ax1.plot([-2, 16], [1.229, 0.401], color="red", linestyle="--", alpha=0.7, linewidth=2)

label_phase = "Fe(s) + Cr2O3(s)"
i = 0
j = 0

for vertices in vertices_list:
    name = generate_entry_label(list_domains[j])
    j += 1
    center_x = sum([v[0] for v in vertices])/len(vertices)
    center_y = sum([v[1] for v in vertices])/len(vertices)

    matched = False
    reversed_matched = False
    for ref_name in reference_domains:
        if name == ref_name:
            matched = True
            label = name
            break
        elif " + ".join(name.split(" + ")[::-1]) == ref_name:
            reversed_matched = True
            label = ref_name
            break

    if not matched and not reversed_matched:
        reference_domains.append(name)
        color_phase[name] = "#CCFFCC"
        label = name

    if label == label_phase:
        patch = Polygon(vertices, closed=True, fill=True, facecolor="#FFA500", linewidth=2, edgecolor="black")
    else:
        col = color=color_phase[label]
        i += 1
        patch = Polygon(vertices, closed=True, fill=True, facecolor=col, linewidth=2, edgecolor="black")
    
#    ax1.text(center_x, center_y, label, verticalalignment="center", horizontalalignment="center", fontname=font, fontsize=18)
    ax1.add_patch(patch)



# Display plot
plt.subplots_adjust(left=0.06, right=0.99, top=0.99, bottom=0.06)
plt.savefig('Figure-c=10-1-Fe17Cr-effet-germination-Cr2O3-DP-2D-modifié-0.2eV-slab.png', dpi=800)
plt.show()


