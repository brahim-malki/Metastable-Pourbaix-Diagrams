# -*- coding: utf-8 -*-
"""
Created on Fri Apr  7 17:55:06 2023

@author: Brahim
"""

from mp_api.client import MPRester
from pymatgen.analysis.pourbaix_diagram import PourbaixDiagram, PourbaixPlotter, generate_entry_label
import matplotlib.pyplot as plt
from matplotlib.patches import Polygon
import seaborn as sb
import numpy as np
from scipy.interpolate import interp1d
from random import seed
from random import random
import pandas as pd
from scipy.interpolate import interp1d

from pymatgen.core.ion import Ion
from pymatgen.entries.computed_entries import ComputedEntry
from pymatgen.core.composition import Composition

from pymatgen.entries.compatibility import MaterialsProjectCompatibility

import copy

from pymatgen.analysis.pourbaix_diagram import PourbaixDiagram, PourbaixPlotter, PourbaixEntry, ComputedEntry, IonEntry


# Initialize the MP Rester
if __name__ == "__main__":
    MAPI_KEY = ""  # You must change this to your Materials API key! (or set MAPI_KEY env variable)  
mpr = MPRester(MAPI_KEY)  # object for connecting to MP Rest interface

# Get all pourbaix entries corresponding to the Ti-O-H chemical system.
entries_init = mpr.get_pourbaix_entries(["Fe", "Cr"])

entries = entries_init
concentration = 1e-8
# tracé d'bord du diagramme de pourbaix 2D
pbx = PourbaixDiagram(entries, comp_dict={"Fe": 0.83, "Cr": 0.17  }, conc_dict={"Fe": concentration, "Cr": concentration})


# return the list of stable and unstable entries
stable_entries = pbx.stable_entries
unstable_entries = pbx.unstable_entries
# define (pH, phi) limits 
limits = [[4.0, 9.0],[-0.4, 0.5]]

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
ax1.set_xlim([4.0, 9.0])
ax1.set_xticklabels([int(t) for t in ax1.get_xticks()], fontname=font, fontsize=38)
ax1.set_ylim(-0.4, 0.5)
yticks = np.arange(-0.4, 0.5, 0.3)
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
#plt.savefig('Figure-c=10-6-Fe17Cr-stable-PD-slabel.png', dpi=800)
plt.show()





###############################################################################
################## calcul du diagramme métastable)   ##########################
###############################################################################


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
            print(i)
            list_index.append(i)
            list_name.append(entry.entry.name)

#list_index_phase_to_modify = list_index
list_index_phase_to_modify = [328]  # au cas ou il faut changer qune seule phase


# list of surface energies, shape factors and molar volumes of entries to modify 
surface_energies = dict()
gamma_Fe2O3 = 0.1

surface_energies = { \
#        compound to add
#        315: [0.75, 5.84, 26.0, entries[315].uncorrected_energy],   \
#        223: [0.57, 5.84, 26.0, entries[223].uncorrected_energy],   \
#        244: [0.79, 5.84, 26.0, entries[244].uncorrected_energy],   \
#        107: [0.99, 5.84, 26.0, entries[107].uncorrected_energy],   \
#        186: [0.40, 5.84, 26.0, entries[186].uncorrected_energy],   \
#        210: [0.34, 5.84, 26.0, entries[210].uncorrected_energy],   \ 
#
#        active compound: Fe2O3 compounds mp-id : 77119
         328: [0.75, 5.84, 26.0, gamma_Fe2O3]    \
       }
                          

    #definition de la fonction d'interpolation en points proches
def interpolate_adjacent_points(x1, y1, y2, num_points=50):
        new_x = np.linspace(x1, x1 + 1, num_points)
        interp_func_y1 = interp1d([x1, x1 + 1], [y1, y2], kind='linear')
        new_y1 = interp_func_y1(new_x)
        return new_x, new_y1

    
    
    
#  fonction de changement d'energie qui me manquait tant !

def modify_entries(entries_init, list_index_phase_to_modify, radius):
    factor = 6.24e-3 # unit conversion factor
    modified_entries = []
    for i, entry in enumerate(entries_init):
        if i in list_index_phase_to_modify:
            initial_entry = entry.entry  # Récupérer l'entrée initiale à partir de l'objet PourbaixEntry
#           print(f"Original energy: {initial_entry.energy_per_atom}")
            delta_energy = ((surface_energies.get(i)[0]*surface_energies.get(i)[1]*surface_energies.get(i)[2])*factor)/radius
#            print(delta_energy)
            modified_energy = initial_entry.energy_per_atom + delta_energy
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

# return the list of stable and unstable entries
stable_entries = pbx.stable_entries
#unstable_entries = pbx.unstable_entries


# define (pH, phi) limits et initialisation des couleurs et du dictionnaire phase_type_scatter, dt et interateur l = 0
limits = [[6.4, 6.8],[-0.5, 0.5]]
colors = sb.color_palette("tab20c", len(stable_entries))
color_phase =dict()

# initialisation de phase_type_scatter
phase_type_scatter = {}

# set radius range in nanoetre from 200 nm to 1nm  
d_min =  2.0
d_max =  -0.5
d_num =  800

#dt = np.linspace(d_min, d_max, d_num)
dt = np.logspace(d_min, d_max, d_num)

l = 0
   

# boucle modifiant progressivement les énergies d'une liste de phases pour voir l'effet sur le diagramme metastable dans le plan (dt[i], V)

for i in range(len(dt)):
    # modify entries
#    increment = (d_max - (d_min)) / (d_num - 1)

#    Créez une copie des entrées initiales pour les conserver inchangées.
    entries_init_copy = copy.deepcopy(entries_init)

    # Appelez la fonction modify_entries avec la copie des entrées initiales et dt[i].
    modified_entries = modify_entries(entries_init_copy, list_index_phase_to_modify, dt[i])          
                   
    pbx1 = PourbaixDiagram(modified_entries, comp_dict={"Fe": 0.83, "Cr": 0.17}, conc_dict={"Fe": concentration, "Cr": concentration})
    stable_entries = pbx1.stable_entries
    stable_domain = pbx1.get_pourbaix_domains(stable_entries, limits)
    list_domains = list(stable_domain[1].keys())
    vertices_list = list(stable_domain[1].values())

    for w in range(0, len(list_domains)):
        list_n = list(color_phase.keys())
        name = list_domains[w].name
        if not (name in list_n):  # needed to fix stable domain color
            a = name.split()
            b = a[::-1]
            c = " ".join(b)
            if not (c in list_n):
                color_phase[name] = colors[l]
                l += 1
            else:
                color_phase[name] = color_phase[list_n[list_n.index(c)]]

        list_p = list(phase_type_scatter.keys())
        if not (name in list_p):
            a_ = name.split()
            b_ = a_[::-1]
            c_ = " ".join(b_)
            if not (c_ in list_p):
                phase_type_scatter[list_domains[w].name] = {'x1': [], 'y1': [], 'y2': [], 'mp-id': []}  # fill phase_type_scatter dict
                dict_tmp = phase_type_scatter[list_domains[w].name]
            else:
                name = c_
                dict_tmp = phase_type_scatter[c_]
        else:
            dict_tmp = phase_type_scatter[name]

#        if (len(vertices_list[w]) == 4):  # save domain coordinates from small polygons
#        if (len(vertices_list[w]) == 4) or (len(vertices_list[w]) == 5) or (len(vertices_list[w]) == 3): # Modifiez cette condition
        if (len(vertices_list[w]) == 4) or (len(vertices_list[w]) == 5) or (len(vertices_list[w]) == 3) or (len(vertices_list[w]) == 6): # Modifiez cette condition

            tmp_x = sorted(vertices_list[w][:, 0])
            tmp_y = sorted(vertices_list[w][:, 1])

            dict_tmp['x1'].append(dt[i])

            if len(vertices_list[w]) == 3:
                dict_tmp['y1'].append((tmp_y[0] + tmp_y[1]) / 2)
                dict_tmp['y2'].append(tmp_y[2])
            else:
                dict_tmp['y1'].append((tmp_y[0] + tmp_y[1]) / 2)
                dict_tmp['y2'].append((tmp_y[2] + tmp_y[3]) / 2)
            dict_tmp['mp-id'].append(list_domains[w].entry_id)
            phase_type_scatter[name].update(dict_tmp)

        else:
            print(w)
            print('vertices_list[w]')
            print(vertices_list[w])
# Ajoutez ces lignes pour imprimer les informations sur les polygones ignorés
            print(f"Ignored polygon at dt[{i}] ({dt[i]}):")
            print(f"  Phase: {name}")
            print(f"  Vertices: {vertices_list[w]}")
            print("\n")
            phase_type_scatter[name].update(dict_tmp)



# tracé du diagramme métastable et affichage des labels des phases directement sur le graphique 

import matplotlib.cm as cm
import seaborn as sb
import numpy as np
from matplotlib.ticker import FormatStrFormatter

unique_phases = list(phase_type_scatter.keys())
colors = sb.color_palette("Paired", len(unique_phases))
#colors = sb.color_palette("Blues", len(unique_phases))
color_phase = dict(zip(unique_phases, colors))
#color_phase["Cr2O3(s) + Fe(s)"] = "yellow"
#color_phase['Fe(s) + Cr2O3(s)'] = (1.0, 1.0, 0.75)
#color_phase["Fe(s) + Cr2O3(s)"] = "#FFA500"

fig, ax = plt.subplots()

num_points = 400


#utilisation de scipy pour interpoler des valeurs nonlineaires

from scipy.interpolate import UnivariateSpline
from scipy.interpolate import CubicSpline, PchipInterpolator

for phase_name, coords_dict in phase_type_scatter.items():
    x_list  = coords_dict['x1']
    y1_list = coords_dict['y1']
    y2_list = coords_dict['y2']

    if len(x_list) > 0 and len(y1_list) > 0 and len(y2_list) > 0:
        # Tri des données en x et y
        sorted_indices = np.argsort(x_list)
        x_list_sorted = np.array(x_list)[sorted_indices]
        y1_list_sorted = np.array(y1_list)[sorted_indices]
        y2_list_sorted = np.array(y2_list)[sorted_indices]
        
        x_new = np.linspace(min(x_list_sorted), max(x_list_sorted), num_points)
        
        # Utilisation de l'interpolation spline
#        spline_y1 = UnivariateSpline(x_list_sorted, y1_list_sorted, s=0.1)
#        spline_y2 = UnivariateSpline(x_list_sorted, y2_list_sorted, s=0.1)

        # ou Utilisation de PchipInterpolator
        spline_y1 = PchipInterpolator(x_list_sorted, y1_list_sorted, extrapolate = 'yes')
        spline_y2 = PchipInterpolator(x_list_sorted, y2_list_sorted, extrapolate = 'yes')

        # Utilisation de CubicSpline
#       spline_y1 = CubicSpline(x_list_sorted, y1_list_sorted)
#        spline_y2 = CubicSpline(x_list_sorted, y2_list_sorted)

        y1_new = spline_y1(x_new)
        y2_new = spline_y2(x_new)

        ax.plot(x_new, y1_new, color='black', linewidth=0.5)
        ax.plot(x_new, y2_new, color='black', linewidth=0.5)
        
        ax.fill_between(x_new, y1_new, y2_new, alpha=0.7, color=color_phase[phase_name])
        
        # Find the center position for the label
        x_center = np.mean(x_new)
        y_center = np.mean([np.mean(y1_new), np.mean(y2_new)])
        
        # Add the label using the `annotate` function
        ax.annotate(phase_name, (x_center, y_center), fontsize=8, ha='center', va='center')



# Plotting details...
font = "serif"
plt.xlim(dt[0], dt[-1])
plt.gca().set_xlim([100, 0.5])
#plt.ylim(limits[1][0], limits[1][1])
plt.ylim(-0.1, 0.3)

plt.xlabel("Stabilisation factor", fontname=font, fontsize=12)
plt.ylabel("Potential vs. SHE (V)", fontname=font, fontsize=12)

plt.gca().set_xscale("log") # Correction ici
# Add more ticks on x axis
xticks = np.linspace(dt[0], dt[-1], num=11)
plt.xticks(xticks)

# Format y axis labels to show one decimal place
plt.gca().yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

#plt.savefig('Figure-4-c=10-1-all-solids-avec-label.png', dpi=800) 

plt.show()
