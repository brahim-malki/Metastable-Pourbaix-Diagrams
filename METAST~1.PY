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
    MAPI_KEY = "SNcD3tHhDv9Y6aOKDq0kKv7acAkmKeD8"  # You must change this to your Materials API key! (or set MAPI_KEY env variable)  
mpr = MPRester(MAPI_KEY)  # object for connecting to MP Rest interface

# Get all pourbaix entries corresponding to the Ti-O-H chemical system.
entries_init = mpr.get_pourbaix_entries(["Fe", "Cr"])

entries = entries_init
concentration = 1e-1

# with fixed concentration
pbx = PourbaixDiagram(entries, comp_dict={"Fe": 0.83, "Cr": 0.17  }, conc_dict={"Fe": concentration, "Cr": concentration})


# return the list of stable and unstable entries
stable_entries = pbx.stable_entries
unstable_entries = pbx.unstable_entries


#definition de la fonction d'interpolation en points proches
def interpolate_adjacent_points(x1, y1, y2, num_points=50):
    new_x = np.linspace(x1, x1 + 1, num_points)
    interp_func_y1 = interp1d([x1, x1 + 1], [y1, y2], kind='linear')
    new_y1 = interp_func_y1(new_x)
    return new_x, new_y1


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
    for i, entry in enumerate(entries):
        if entry.entry_id == id:
            print(entry.entry.name)
            list_index.append(i)
            list_name.append(entry.entry.name)

list_index_phase_to_modify = list_index
#list_index_phase_to_modify = [64]  # au cas ou il faut changer qune seule phase


# define (pH, phi) limits et initialisation des couleurs et du dictionnaire phase_type_scatter, dt et interateur l = 0
pH = 7.0
limits = [[pH-0.1, pH+0.1],[-1.1, 1.1]]

colors = sb.color_palette("tab20c", len(stable_entries))
color_phase =dict()

# initialisation de phase_type_scatter
phase_type_scatter = {}

# set radius range in nanoetre from 200 nm to 1nm  
d_min = -2.0
d_max =  2.0
d_num =  400

dt = np.linspace(d_min, d_max, d_num)
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
        if (len(vertices_list[w]) == 4) or (len(vertices_list[w]) == 5) or (len(vertices_list[w]) == 3): # Modifiez cette condition
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
# Ajoutez ces lignes pour imprimer les informations sur les polygones ignorés
            print(f"Ignored polygon at dt[{i}] ({dt[i]}):")
            print(f"  Phase: {name}")
            print(f"  Vertices: {vertices_list[w]}")
            print("\n")
            phase_type_scatter[name].update(dict_tmp)

 


# methode 1 : tracé du diagramme métastable et affichage des labels des phases directement sur le graphique 

import matplotlib.cm as cm
import seaborn as sb
import numpy as np
from matplotlib.ticker import FormatStrFormatter

unique_phases = list(phase_type_scatter.keys())
#colors = sb.color_palette("Paired", len(unique_phases))
colors = sb.color_palette("Blues", len(unique_phases))
color_phase = dict(zip(unique_phases, colors))
#color_phase["Cr2O3(s) + Fe(s)"] = "yellow"
#color_phase['Fe(s) + Cr2O3(s)'] = (1.0, 1.0, 0.75)
color_phase["Fe(s) + Cr2O3(s)"] = "#FFA500"

fig, ax = plt.subplots()

num_points = 200

for phase_name, coords_dict in phase_type_scatter.items():
    x_list  = coords_dict['x1']
    y1_list = coords_dict['y1']
    y2_list = coords_dict['y2']

    if len(x_list) > 0 and len(y1_list) > 0 and len(y2_list) > 0:
        x_new = np.linspace(min(x_list), max(x_list), num_points)
        y1_new = np.interp(x_new, x_list, y1_list)
        y2_new = np.interp(x_new, x_list, y2_list)
        
        ax.plot(x_new, y1_new, color='black', linewidth=1.0)
        ax.plot(x_new, y2_new, color='black', linewidth=1.0)
        
        ax.fill_between(x_new, y1_new, y2_new, alpha=0.8, color=color_phase[phase_name])
        
        # Find the center position for the label
        x_center = np.mean(x_new)
        y_center = np.mean([np.mean(y1_new), np.mean(y2_new)])
        
        # Add the label using the `annotate` function
        ax.annotate(phase_name, (x_center, y_center), fontsize=8, ha='center', va='center')



# Plotting details...
font = "serif"

# Plotting details...
font = "serif"
plt.xlim(dt[0], dt[-1])
plt.ylim(limits[1][0], limits[1][1])
plt.xlabel("Stabilisation factor", fontname=font, fontsize=12)
plt.ylabel("Potential vs. SHE (V)", fontname=font, fontsize=12)

# Add more ticks on x axis
xticks = np.linspace(dt[0], dt[-1], num=11)
plt.xticks(xticks)

# Format y axis labels to show one decimal place
plt.gca().yaxis.set_major_formatter(FormatStrFormatter('%.1f'))

plt.savefig('Figure-4-c=10-1-all-solids-avec-label.png', dpi=800) 

plt.show()








# methode 2 : tracé du diagramme métastable et affichage sans labels des phases

import matplotlib.cm as cm
import seaborn as sb
import numpy as np
from matplotlib.ticker import FormatStrFormatter

unique_phases = list(phase_type_scatter.keys())
#colors = sb.color_palette("Paired", len(unique_phases))
colors = sb.color_palette("Blues", len(unique_phases))
color_phase = dict(zip(unique_phases, colors))
color_phase = dict(zip(unique_phases, colors))
#color_phase["Cr2O3(s) + Fe(s)"] = "yellow"
#color_phase["Cr2O3(s) + Fe(s)"] = (1.0, 1.0, 0.75)
color_phase["Fe(s) + Cr2O3(s)"] = "#FFA500"



fig, ax = plt.subplots()

num_points = 100

for phase_name, coords_dict in phase_type_scatter.items():
    x_list  = coords_dict['x1']
    y1_list = coords_dict['y1']
    y2_list = coords_dict['y2']

    if len(x_list) > 0 and len(y1_list) > 0 and len(y2_list) > 0:
        x_new = np.linspace(min(x_list), max(x_list), num_points)
        y1_new = np.interp(x_new, x_list, y1_list)
        y2_new = np.interp(x_new, x_list, y2_list)
        
        ax.plot(x_new, y1_new, color='black', linewidth=1.0)
        ax.plot(x_new, y2_new, color='black', linewidth=1.0)
        
        ax.fill_between(x_new, y1_new, y2_new, alpha=0.8, color=color_phase[phase_name], label=phase_name)


# Plotting details...
font = "serif"
# Plotting details...
font = "serif"
plt.xlim(dt[0], dt[-1])
plt.ylim(limits[1][0], limits[1][1])
plt.xlabel("Stabilisation factor", fontname=font, fontsize=12)
plt.ylabel("Potential vs. SHE (V)", fontname=font, fontsize=12)

# Add more ticks on x axis
xticks = np.linspace(dt[0], dt[-1], num=11)
plt.xticks(xticks)

# Format y axis labels to show one decimal place
plt.gca().yaxis.set_major_formatter(FormatStrFormatter('%.1f'))
plt.subplots_adjust(left=0.06, right=0.95, top=0.95, bottom=0.06)

plt.savefig('Figure-4-c=10-1-all-solids-Slabel.png', dpi=800) 

plt.show()


# pour afficher la bare des labels Avoid duplicate labels in the legend
#handles, labels = plt.gca().get_legend_handles_labels()
#by_label = dict(zip(labels, handles))
#plt.legend(by_label.values(), by_label.keys())





# methode 2 tracé du diagramme métastable et interpolation avec cubispline et affichage des labels des phases directement sur le graphique 

import matplotlib.cm as cm
import seaborn as sb
import numpy as np
from scipy.interpolate import CubicSpline
from matplotlib.ticker import FormatStrFormatter



unique_phases = list(phase_type_scatter.keys())
#colors = sb.color_palette("Paired", len(unique_phases))
colors = sb.color_palette("Blues", len(unique_phases))
color_phase = dict(zip(unique_phases, colors))
color_phase = dict(zip(unique_phases, colors))
color_phase["Cr2O3(s) + Fe(s)"] = "yellow"

fig, ax = plt.subplots()

num_points = 200

for phase_name, coords_dict in phase_type_scatter.items():
    x_list  = coords_dict['x1']
    y1_list = coords_dict['y1']
    y2_list = coords_dict['y2']

    if len(x_list) > 0 and len(y1_list) > 0 and len(y2_list) > 0:
        x_new = np.linspace(min(x_list), max(x_list), num_points)
        
        # Utiliser la méthode CubicSpline pour l'interpolation
        spline_y1 = CubicSpline(x_list, y1_list)
        spline_y2 = CubicSpline(x_list, y2_list)
        y1_new = spline_y1(x_new)
        y2_new = spline_y2(x_new)
        
        ax.plot(x_new, y1_new, color='black', linewidth=1.0)
        ax.plot(x_new, y2_new, color='black', linewidth=1.0)
        
        ax.fill_between(x_new, y1_new, y2_new, alpha=0.8, color=color_phase[phase_name])
        
        # Find the center position for the label
        x_center = np.mean(x_new)
        y_center = np.mean([np.mean(y1_new), np.mean(y2_new)])
        
        # Add the label using the `annotate` function
        ax.annotate(phase_name, (x_center, y_center), fontsize=8, ha='center', va='center')



# Plotting details...
font = "serif"
plt.xlim(dt[0], dt[-1])
plt.ylim(limits[1][0], limits[1][1])
plt.xlabel("Stabilisation factor", fontname=font, fontsize=12)
plt.ylabel("Potential vs. SHE (V)", fontname=font, fontsize=12)
plt.show()




###############################################################################
#########################            Annexes              #####################
###############################################################################


#obtenir l'énergie d'un matériau à partir de Materials Project et ensuite créer une nouvelle IonEntry avec l'énergie modifiée

from mp_api.client import MPRester
from pymatgen.analysis.pourbaix_diagram import IonEntry
from pymatgen.core.ion import Ion

MAPI_KEY = "SNcD3tHhDv9Y6aOKDq0kKv7acAkmKeD8"  # Remplacez par votre clé API
material_id = 'mp-19399'  # Remplacez par l'ID du matériau correspondant à Cr2O3

# Obtenez l'énergie du matériau à partir de Materials Project
with MPRester(MAPI_KEY) as mpr:
    entry = mpr.get_entry_by_material_id(material_id, compatible_only=True)
    entry = entries[0]
    energy = entry.energy_per_atom
    
# Créez une IonEntry existante avec l'énergie obtenue
ion = Ion.from_formula('Cr2O3')
entry = IonEntry(ion, energy)

# Modifiez l'énergie
new_energy = energy + 0.1

# Créez une nouvelle IonEntry avec l'énergie modifiée
modified_entry = IonEntry(ion, new_energy)

print(modified_entry)



# modifier l'energie d'un composé Cr2O3 et voir l'effet sur le diagramme de pourbaix 

from mp_api.client import MPRester
from pymatgen.analysis.pourbaix_diagram import PourbaixDiagram, PourbaixPlotter, IonEntry, PourbaixEntry
from pymatgen.core.ion import Ion

MAPI_KEY = "SNcD3tHhDv9Y6aOKDq0kKv7acAkmKeD8"  # Remplacez par votre clé API

# Obtenez toutes les entrées pour les diagrammes de Pourbaix pour les éléments donnés (Fer et Chrome)
with MPRester(MAPI_KEY) as mpr:
    fe_entries = mpr.get_pourbaix_entries(["Fe"])
    cr_entries = mpr.get_pourbaix_entries(["Cr"])

# Identifiez l'entrée correspondant à Cr2O3 et modifiez son énergie

modified_entries = []
for entry in cr_entries:
    if entry.composition.reduced_formula == "Cr2O3":
        new_energy = entry.energy + 1.0  # Modifiez l'énergie ici
        modified_entry = PourbaixEntry(entry.entry, new_energy)  # Créez une nouvelle instance de PourbaixEntry avec la nouvelle énergie
        modified_entries.append(modified_entry)
    else:
        modified_entries.append(entry)
        
# Combine Fe et Cr entries
modified_entries.extend(fe_entries)

# Définissez les concentrations pour Fe et Cr
concentration = 1e-8  # Remplacez cette valeur par la concentration souhaitée

# Créez le diagramme de Pourbaix avec les entrées modifiées pour le système Fe-Cr
pbx = PourbaixDiagram(modified_entries, comp_dict={"Fe": 0.83, "Cr": 0.17}, conc_dict={"Fe": concentration, "Cr": concentration})
plotter = PourbaixPlotter(pbx)
plt = plotter.get_pourbaix_plot()
plt.show()


# un code similaire avec une liste de phase et sur une boucle en new_energy
# pour l'instant la méthode utilisée ne permet pas de changer 
   
        
        
import numpy as np
from pymatgen.analysis.pourbaix_diagram import PourbaixDiagram, PourbaixPlotter, PourbaixEntry
from pymatgen.ext.matproj import MPRester
from pymatgen.core import Element

# Configuration du MPRester et récupération des données pour les éléments Fe et Cr
MAPI_KEY = "SNcD3tHhDv9Y6aOKDq0kKv7acAkmKeD8"  # Remplacez par votre clé API
mpr = MPRester(MAPI_KEY)

# Concentration des éléments
concentration = 1e-8

# Paramètres pour la boucle sur new_energy
d_min = -1
d_max = 1
d_num = 10

# Liste des phases à modifier
list_index_phase_to_modify = [64]

# Récupération des entrées initiales de Pourbaix
entries = mpr.get_pourbaix_entries(["Fe", "Cr"])

def get_modified_pourbaix(entries, list_index_phase_to_modify, new_energy_offsets):
    # Création d'une copie des entrées pour les modifier
    modified_entries = entries.copy()

    # Modification de l'énergie des phases sélectionnées
    for index_phase, new_energy_offset in zip(list_index_phase_to_modify, new_energy_offsets):
        entry = entries[index_phase]
        new_energy = entry.energy + new_energy_offset
        modified_entry = PourbaixEntry(entry.entry, new_energy)  # Créez une nouvelle instance de PourbaixEntry avec la nouvelle énergie

        # Remplacement de l'entrée initiale par l'entrée modifiée
        modified_entries[index_phase] = modified_entry

    # Création du diagramme de Pourbaix modifié
    pbx = PourbaixDiagram(modified_entries, comp_dict={"Fe": 0.83, "Cr": 0.17},
                          conc_dict={"Fe": concentration, "Cr": concentration})
    return pbx

      
# Boucle sur les valeurs de new_energy
for new_energy_offset in np.linspace(d_min, d_max, d_num):
    # Création et tracé du diagramme de Pourbaix modifié
    new_energy_offsets = [new_energy_offset] * len(list_index_phase_to_modify)
    pbx = get_modified_pourbaix(entries, list_index_phase_to_modify, new_energy_offsets)
    plotter = PourbaixPlotter(pbx)
    plt = plotter.get_pourbaix_plot()
    plt.title(f"Diagramme de Pourbaix Fe-Cr modifié (dE = {new_energy_offset:.2f})")
    plt.show()









from mp_api.client import MPRester
from pymatgen.analysis.pourbaix_diagram import PourbaixDiagram, PourbaixPlotter, IonEntry, PourbaixEntry, ComputedEntry
from pymatgen.core.ion import Ion

MAPI_KEY = "SNcD3tHhDv9Y6aOKDq0kKv7acAkmKeD8"  # Remplacez par votre clé API


from mp_api.client import MPRester
from pymatgen.analysis.pourbaix_diagram import PourbaixDiagram, PourbaixPlotter, PourbaixEntry, ComputedEntry, IonEntry

# Get energy formation entries
MAPI_KEY = "SNcD3tHhDv9Y6aOKDq0kKv7acAkmKeD8"  # Remplacez par votre clé API

mpr = MPRester(MAPI_KEY)  # object for connecting to MP Rest interface
entries = mpr.get_entries_in_chemsys(["Fe", "Cr", "O", "H"])

# Modify energy of the desired entry
modified_entries = []
for entry in entries:
    if entry.composition.reduced_formula == "Cr2O3":
        print(f"Original energy: {entry.energy_per_atom}")
        modified_energy = entry.energy_per_atom + 1.0
        print(f"Modified energy: {modified_energy}")
        modified_entry = ComputedEntry(entry.composition, modified_energy * entry.composition.num_atoms, parameters=entry.parameters)
        if isinstance(entry, IonEntry):
            modified_entry = IonEntry(modified_entry)
        modified_entry = PourbaixEntry(modified_entry)
    else:
        modified_entry = PourbaixEntry(entry)
    
    # Skip the entry if normalization_factor would result in a division by zero
    if modified_entry.num_atoms - modified_entry.composition.get("H", 0) - modified_entry.composition.get("O", 0) == 0:
        continue
        
    modified_entries.append(modified_entry)





# Generate Pourbaix diagram
conc_dict = {"Fe": 1e-8, "Cr": 1e-8}
pbx = PourbaixDiagram(modified_entries, conc_dict)
plotter = PourbaixPlotter(pbx)
plotter.get_pourbaix_plot().show()





# autre proposition 



from mp_api.client import MPRester
from pymatgen.analysis.pourbaix_diagram import PourbaixDiagram, PourbaixPlotter, PourbaixEntry, ComputedEntry, IonEntry

MAPI_KEY = "SNcD3tHhDv9Y6aOKDq0kKv7acAkmKeD8"  # Remplacez par votre clé API
mpr = MPRester(MAPI_KEY)  # object for connecting to MP Rest interface
entries_init = mpr.get_pourbaix_entries(["Fe", "Cr"])



conc_dict = {"Fe": 1e-8, "Cr": 1e-8}
pbx1 = PourbaixDiagram(entries_init, conc_dict)
plotter = PourbaixPlotter(pbx1)
plotter.get_pourbaix_plot().show()


# Modify energy of the desired entry
modified_entries = []
for entry in entries_init:
    initial_entry = entry.entry  # Récupérer l'entrée initiale à partir de l'objet PourbaixEntry
    if initial_entry.composition.reduced_formula == "Cr2O3":
        print(f"Original energy: {initial_entry.energy_per_atom}")
        modified_energy = initial_entry.energy_per_atom + 0.0
        print(f"Modified energy: {modified_energy}")
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

# Generate Pourbaix diagram
conc_dict = {"Fe": 1e-1, "Cr": 1e-1}
pbx = PourbaixDiagram(modified_entries, conc_dict)
plotter = PourbaixPlotter(pbx)
plotter.get_pourbaix_plot().show()

print("Concentrations utilisées: ", conc_dict)

for entry in modified_entries:
    if entry.composition.reduced_formula == "Cr2O3":
        print(f"Énergie modifiée pour Cr2O3: {entry.energy_per_atom}")














