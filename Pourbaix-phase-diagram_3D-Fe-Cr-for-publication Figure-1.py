# -*- coding: utf-8 -*-T
"""
Created on Sun Oct  6 12:15:43 2019

@author: brahim
"""

###############################################################################
############## version 1 : using collection3d method   ########################
###############################################################################


# from pymatgen import MPRester modified in this pymatgen version as:
from pymatgen.ext.matproj import MPRester 
# or 
from mp_api.client import MPRester

from pymatgen.analysis.pourbaix_diagram import PourbaixDiagram, PourbaixPlotter, generate_entry_label

import numpy as np

import matplotlib.pyplot as plt

import seaborn as sb

import mpl_toolkits.mplot3d as a3
import matplotlib.colors as colors


# Initialize the MP Rester

#if __name__ == "__main__":
#    MAPI_KEY = "9TnOADgQ0sFZ6AAdt"  #  (set MAPI_KEY env variable)
#    system = ["Mn"]  # system we want to get open PD for    
#mpr = MPRester(MAPI_KEY)  # object for connecting to MP Rest interface
 
# use new api key
if __name__ == "__main__":
    MAPI_KEY = "SNcD3tHhDv9Y6aOKDq0kKv7acAkmKeD8"  # You must change this to your Materials API key! (or set MAPI_KEY env variable)  
mpr = MPRester(MAPI_KEY)  # object for connecting to MP Rest interface
   
# Get all pourbaix entries corresponding to the Ti-O-H chemical system.
entries = mpr.get_pourbaix_entries(["Fe", "Cr"])


# with fixed concentration
concentration = 1e-1

# with fixed concentration
pbx = PourbaixDiagram(entries, comp_dict={"Fe": 0.83, "Cr": 0.17}, conc_dict={"Fe": concentration, "Cr": concentration})
#plotter = PourbaixPlotter(pbx)


# return the list of stable and unstable entries
stable_entries = pbx.stable_entries
unstable_entries = pbx.unstable_entries

limits = [[-0.1, 14],[-1.1, 1.1]]
list_domain   = {}

# get st
stable_domain   = pbx.get_pourbaix_domains(stable_entries, limits)
list_domains = list(stable_domain[1].keys())
vertices_list=list(stable_domain[1].values())


collection_2D = []
collection_3D = []
vecteur_2D = np.zeros((1,3))
vecteur_3D = np.zeros((1,3))

for i in range(0, len(list_domains)):
    test = vertices_list[i]
    n = len(test)
    vecteur_2D = np.zeros((n,3))
    vecteur_3D = np.zeros((n,3))
    for j in range(0,len(test)):
        he = pbx.get_hull_energy(test[j][0],test[j][1])
        vecteur_2D[j,:] = [test[j][0], test[j][1], -15.0]
        vecteur_3D[j,:] = [test[j][0], test[j][1], he] 
    collection_2D.append(vecteur_2D)
    collection_3D.append(vecteur_3D)
  
def create_figure():
    font = "serif"  
    ax.view_init(20, 35)
    ax.set_xlim(-0.1, 14.1)
    ax.set_xticks([0, 5, 10, 14])
    ax.set_xticklabels(ax.get_xticks(), fontname=font, fontsize=14)
    ax.set_ylim(-1.1, 1.1)
    ax.set_yticks([-1, 0, 1])
    ax.set_yticklabels(ax.get_yticks(), fontname=font, fontsize=14)
    ax.set_zlim(-15, 2)
    ax.set_zticks([-15, -10, -5, 0, 2])
    ax.set_zticklabels(ax.get_zticks(), fontname=font, fontsize=14)
#    ax.set_xlabel("pH", fontname=font, fontsize=14)
#    ax.set_ylabel("Potential vs. SHE (V)", fontname=font, fontsize=14)
#    ax.set_zlabel("Pourbaix Potential (eV/Mn)", fontname=font, fontsize=14)
    return ax    

fig = plt.figure(figsize=(14, 9))
ax = a3.Axes3D(fig)
create_figure()

colors2D = sb.color_palette("Blues", len(list_domains))
colors3D = sb.color_palette("Blues", len(list_domains))

for i in range(0, len(collection_2D)): 
    vertices_xy = vertices_list[i]
    vertices_z  = collection_3D[i]
    tri_2D = a3.art3d.Poly3DCollection([collection_2D[i]])
    tri_3D = a3.art3d.Poly3DCollection([collection_3D[i]])
    color2D = colors2D[i]  #   or color = colors.rgb2hex(sp.rand(3))
    color3D = colors3D[i]  #   or color = colors.rgb2hex(sp.rand(3)
    tri_2D.set_color(color2D)
    tri_3D.set_color(color3D)
    tri_3D.set_edgecolor('k')
    ax.add_collection3d(tri_2D)
    ax.add_collection3d(tri_3D)
plt.show()


###############################################################################
############## version 1 : using collection3d method   ########################
##############               adding unstable entries   ########################
###############################################################################


# return the list of stable and unstable entries
stable_entries = pbx.stable_entries
unstable_entries = pbx.unstable_entries

limits = [[-0.1, 14],[-1.1, 1.1]]
list_domain   = {}

# get st
stable_domain   = pbx.get_pourbaix_domains(stable_entries, limits)
list_domains = list(stable_domain[1].keys())
vertices_list=list(stable_domain[1].values())

unstable_domain_   = pbx.get_pourbaix_domains(unstable_entries, limits)
list_domains_ = list(unstable_domain_[1].keys())
vertices_list_=list(unstable_domain_[1].values())


collection_3D  = []
collection_3D_ = []

for i in range(0, len(list_domains)):
    test = vertices_list[i]
    n = len(test)
    vecteur_3D = np.zeros((n,3))
    for j in range(0,len(test)):
        he = pbx.get_hull_energy(test[j][0],test[j][1])
        vecteur_3D[j,:] = [test[j][0], test[j][1], he] 
    collection_3D.append(vecteur_3D)

for i in range(0, len(list_domains_)):
    test = vertices_list_[i]
    n = len(test)
    vecteur_3D = np.zeros((n,3))
    for j in range(0,len(test)):
        he = pbx.get_hull_energy(test[j][0],test[j][1])
        vecteur_3D[j,:] = [test[j][0], test[j][1], he] 
    collection_3D_.append(vecteur_3D)
    
    
  
def create_figure():
    font = "serif"  
    ax.view_init(20, 35)
    ax.set_xlim(-0.1, 14.1)
    ax.set_xticks([0, 5, 10, 14])
    ax.set_xticklabels(ax.get_xticks(), fontname=font, fontsize=14)
    ax.set_ylim(-1.1, 1.1)
    ax.set_yticks([-1, 0, 1])
    ax.set_yticklabels(ax.get_yticks(), fontname=font, fontsize=14)
    ax.set_zlim(-15, 2)
    ax.set_zticks([-15, -10, -5, 0, 2])
    ax.set_zticklabels(ax.get_zticks(), fontname=font, fontsize=14)
#    ax.set_xlabel("pH", fontname=font, fontsize=14)
#    ax.set_ylabel("Potential vs. SHE (V)", fontname=font, fontsize=14)
#    ax.set_zlabel("Pourbaix Potential (eV/Mn)", fontname=font, fontsize=14)
    return ax    

ax = a3.Axes3D(plt.figure(figsize=(14, 9)))
create_figure()

colors = sb.color_palette("Blues", len(collection_3D))
#colors = sb.color_palette("Set2", len(collection_3D))
  


for i in range(0, 28): 
    tri_3D = a3.art3d.Poly3DCollection([collection_3D[i]])
    if (i >=11):
        color = (0.5, 0.5, 0.5)
    else:
        color = colors[i]  #   or color = colors.rgb2hex(sp.rand(3))
    tri_3D.set_color(color)
    tri_3D.set_edgecolor('k')
    ax.add_collection3d(tri_3D)
#    ax.set_animated(True)
plt.show()




