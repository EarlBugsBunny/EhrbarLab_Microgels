# -*- coding: utf-8 -*-
"""
Created on Thu Aug  4 14:11:22 2022

@author: lukas
"""

import os
import cv2
import numpy as np

overall_folder = {}

for folder in os.listdir():
    overall_folder[folder] = {}
    for subfolder in os.listdir(folder):
        overall_folder[folder][subfolder] = {}
        same_spot = {}
        for file in os.listdir(folder+'/'+subfolder):
            spot = file[:-10]
            if spot in same_spot.keys():
                same_spot[spot].append(file)
            else:
                same_spot[spot] = [file]
        for key in same_spot.keys():
            all_files_same_spot = same_spot[key]
            all_matrix = []
            for direction in all_files_same_spot:
                matrix = cv2.imread(folder+'/'+subfolder+'/'+direction)[:,:,0]
                all_matrix.append(matrix)
            max_projection = np.maximum.reduce(all_matrix)
            overall_folder[folder][subfolder][key] = max_projection

import matplotlib.pyplot as plt
import numpy as np
from scipy import ndimage as ndi
import cv2
from skimage.filters import threshold_otsu
import os
import matplotlib.pyplot as plt
import numpy as np
from scipy import ndimage as ndi
import cv2
from skimage.filters import threshold_otsu
import pandas as pd

from skimage import (
    color, feature, filters, measure, morphology, segmentation, util
)

def intensity_per_cell(image_dapi,image_yap,size_pixel = 200, min_distance=20):
    thresh = threshold_otsu(image_dapi)
    cells = image_dapi > thresh

    distance = ndi.distance_transform_edt(cells)
    local_max_coords = feature.peak_local_max(distance, min_distance=20)
    '''plat with the figure min_distance!!'''
    local_max_mask = np.zeros(distance.shape, dtype=bool)
    local_max_mask[tuple(local_max_coords.T)] = True
    markers = measure.label(local_max_mask)
    segmented_cells = segmentation.watershed(-distance, markers, mask=cells)
    
    dict_all= {}
    
    for label in np.unique(segmented_cells):
        if label == 0:
            continue
        area = np.count_nonzero(segmented_cells == label)
        if area > size_pixel:
            mask = (segmented_cells == label)
            average_dapi = np.mean(image_dapi, where=mask)
            average_yap = np.mean(image_yap, where=mask)
            dict_all[label] = [area, average_dapi, average_yap]
    
    return dict_all

column_names = ['condition', 'image_name','label','area','dapi_intensity', 'YAP_intensity']
df = pd.DataFrame(columns = column_names)

for folder in overall_folder.keys():
    for list_images in overall_folder[folder]['dapi'].keys():
        matrix_dapi = overall_folder[folder]['dapi'][list_images]
        matrix_yap = overall_folder[folder]['yap'][list_images]
        dict_cells = intensity_per_cell(matrix_dapi, matrix_yap)
        for key in dict_cells.keys():
            new_row = {'condition': folder, 'image_name': list_images,'label': key,
                       'area': dict_cells[key][0],
                       'dapi_intensity': dict_cells[key][1],
                       'YAP_intensity': dict_cells[key][2]}
            df = df.append(new_row, ignore_index = True)
            
        
        
#%%

import seaborn as sns
df["area"] = df.area.astype(float)
df["label"] = df.label.astype(float)
df["dapi_intensity"] = df.dapi_intensity.astype(float)
df["YAP_intensity"] = df.YAP_intensity.astype(float)

df["YAP_normalized"] = df['YAP_intensity'] / df['dapi_intensity']
#%%
df.to_csv('int_per_cell.csv')
#%%

# ax = sns.violinplot(x="condition", y="YAP_intensity", data=df,
#                     order=[ "soft", "medium", "stiff"])
# ax = sns.violinplot(x="condition", y="YAP_normalized", data=df, 
#                     order=[ "soft", "medium", "stiff"])
# ax = sns.violinplot(x="condition", y="YAP_intensity", data=df,
#                     order=[ "small", "medium", "large"])
ax = sns.violinplot(x="condition", y="YAP_normalized", data=df, 
                    order=[ "small", "medium", "large"])

                
#%%

import seaborn as sns
import matplotlib.pyplot as plt

plt.figure(figsize=(14,8))

cat = ['soft', 'medium', 'stiff']
cat = ['small', 'medium', 'large']
sns.color_palette('colorblind')

df['condition'] = pd.Categorical(df['condition'], cat)

ax = sns.histplot(data=df, x="YAP_normalized", hue="condition", alpha=0.5, kde=True, 
             line_kws={'color': 'crimson', 'lw':5})

ax.set_xlabel(xlabel='Intensity of YAP in the Nucleus (Normalized by DAPI)', size = 25)
ax.set_ylabel(ylabel='Counted Nuclei', size = 25)
ax.set_xlim(-0.0,2)
ax.tick_params(axis='both', which='major', labelsize=15)
plt.setp(ax.get_legend().get_texts(), fontsize='30') # for legend text
plt.setp(ax.get_legend().get_title(), fontsize='35') # for legend title

plt.show()
#%%
ax = sns.kdeplot(data=df, x="YAP_normalized", hue="condition",lw=3)

            
            