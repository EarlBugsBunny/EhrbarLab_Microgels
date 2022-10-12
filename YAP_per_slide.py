# -*- coding: utf-8 -*-
"""
Created on Thu Oct 28 15:47:53 2021
@author: Lukas
"""
import cv2
import numpy as np
import os

'With induvidual threshold'
#colors
light_red = '#CA6573'
dark_red = '#791422'
light_violet = '#625192'
dark_violet = '#281758'
light_green = '#9DC462'
dark_green = '#4F7514'
light_yellow = '#D6D147'
dark_yellow = '#99940E'
my_pal = [light_red,light_green, light_violet, light_yellow]

def get_ratio_yap(yap_picture, dapi_image, actin_image):
    threshold_image_dapi = cv2.threshold(dapi_image, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)[1]
    boolean_matrix_dapi = np.array(threshold_image_dapi, dtype=bool)
    count_dapi = np.count_nonzero(boolean_matrix_dapi)
    if count_dapi < 50:
        return float('NaN')
    average_nucleus = np.average(yap_picture, weights = boolean_matrix_dapi)
    # yap_for_mean = yap_picture[boolean_matrix_dapi]
    # average_nucleus = np.median(yap_for_mean)
    
    
    threshold_image_actin = cv2.threshold(actin_image, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)[1]
    boolean_matrix_actin = np.array(threshold_image_actin, dtype=bool)
    count_actin = np.count_nonzero(boolean_matrix_actin)
    if count_actin < 50:
        return float('NaN')
    
    not_nucleus = np.logical_not(boolean_matrix_dapi)
    boolean_matrix_cytosol = np.logical_and(boolean_matrix_actin,not_nucleus)
    count_cytosol = np.count_nonzero(boolean_matrix_cytosol)
    if count_cytosol < 50:
        return float('NaN')
    
    average_cytosol = np.average(yap_picture, weights = boolean_matrix_cytosol)
    # yap_for_mean_actin = yap_picture[boolean_matrix_cytosol]
    # average_cytosol = np.median(yap_for_mean_actin)
    ratio = average_nucleus / average_cytosol
    
    return ratio


overall_folder = os.listdir()
dic_all_values = {}

for folder in overall_folder:
    dic_all_values[folder] = []
    list_all_images = []
    folder_list = os.listdir(folder)
    for subfolder in folder_list:
        liste_images =  os.listdir(folder+'/'+subfolder)
        liste_images = [subfolder+'/'+ s  for s in liste_images]
        list_all_images.append(liste_images)
    for count, actin_image_name in enumerate(list_all_images[0]):
        actin_image = cv2.imread(folder+'/'+actin_image_name)[:,:,0]
        dapi_image = cv2.imread(folder+'/'+list_all_images[1][count])[:,:,0]
        yap_image = cv2.imread(folder+'/'+list_all_images[2][count])[:,:,0]
        ratio = get_ratio_yap(yap_image, dapi_image, actin_image)
        dic_all_values[folder].append(ratio)
        
for key in dic_all_values.keys():
    liste = dic_all_values[key]
    liste = np.array(liste)
    liste = liste[~np.isnan(liste)]
    liste = liste[liste < 1E308]
    print(key)
    print(np.mean(liste))
    print(np.median(liste))

#Plotting of a figure to see result
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
sns.set_style("darkgrid")

over_all_list = []
percent = [5,7.5,10,12.5]
i = 0

for key in dic_all_values.keys():
    liste = dic_all_values[key]
    liste = np.array(liste)
    liste = liste[~np.isnan(liste)]
    liste = liste[liste < 8]
    for num in liste:
        over_all_list.append([percent[i],num])
    i += 1
    
df = pd.DataFrame(over_all_list,columns=['Percentage of PEG-NB Gel','Ratio (Intensity of Nucleus/Cytosol)'])
  
 
p = sns.violinplot(x='Percentage of PEG-NB Gel', y='Ratio (Intensity of Nucleus/Cytosol)', data=df, palette=my_pal)
l = sns.stripplot(x='Percentage of PEG-NB Gel', y='Ratio (Intensity of Nucleus/Cytosol)', data=df, jitter=True, color='black', s= 7)
l.set_xlabel('Percentage of PEG-NB Gel', fontsize=30)
l.set_ylabel('Ratio (Intensity of Nucleus/Cytosol)', fontsize=30)
l.tick_params(axis='both', which='major', labelsize=30)
l.set_ylim(bottom=None, top=10)
l.set_title('Yap location based on intensity of antibody staining', fontsize=50)
plt.show()



""""Alternatives which were tried and gave worse results:"""

"""
#%%
import cv2
import numpy as np
import os
import skimage
from skimage import filters

'''With all images the same threshold (but induvidual channels different ones)'''


overall_folder = os.listdir()
dic_all_values = {}

actin_all = []
dapi_all = []
yap_all = []

overall_flat_actin = np.asarray([])
overall_flat_dapi = np.asarray([])
overall_flat_yap = np.asarray([])

for folder in overall_folder:
    folder_list = os.listdir(folder)
    for sub_folder in folder_list:
        if 'actin' in sub_folder:
            print(sub_folder)
            pictures = os.listdir(folder+'/'+sub_folder)
            for picture in pictures:
                image = cv2.imread(folder+'/'+sub_folder+'/'+picture)[:,:,0]
                overall_flat_actin = np.append(np.copy(overall_flat_actin), np.copy(image.flatten()))
        if 'dapi' in sub_folder:
            print(sub_folder)
            pictures = os.listdir(folder+'/'+sub_folder)
            for picture in pictures:
                image = cv2.imread(folder+'/'+sub_folder+'/'+picture)[:,:,0]
                overall_flat_dapi = np.append(np.copy(overall_flat_dapi), np.copy(image.flatten()))
        if 'yap' in sub_folder:
            print(sub_folder)
            pictures = os.listdir(folder+'/'+sub_folder)
            for picture in pictures:
                image = cv2.imread(folder+'/'+sub_folder+'/'+picture)[:,:,0]
                overall_flat_yap = np.append(np.copy(overall_flat_yap), np.copy(image.flatten()))

otsu_actin = skimage.filters.threshold_otsu(overall_flat_actin)
otsu_dapi = skimage.filters.threshold_otsu(overall_flat_dapi)
otsu_yap = skimage.filters.threshold_otsu(overall_flat_yap)   



def get_ratio_yap(yap_picture, dapi_image, actin_image):
    boolean_matrix_dapi = dapi_image > otsu_dapi
    count_dapi = np.count_nonzero(boolean_matrix_dapi)
    if count_dapi < 50:
        return float('NaN')
    
    average_nucleus = np.average(yap_picture, weights = boolean_matrix_dapi)
    
    boolean_matrix_actin = actin_image > otsu_actin
    count_actin = np.count_nonzero(boolean_matrix_actin)
    if count_actin < 50:
        return float('NaN')
    
    not_nucleus = np.logical_not(boolean_matrix_dapi)
    boolean_matrix_cytosol = np.logical_and(boolean_matrix_actin,not_nucleus)
    
    average_cytosol = np.average(yap_picture, weights = boolean_matrix_cytosol)
    
    ratio = average_nucleus / average_cytosol
    
    return ratio


overall_folder = os.listdir()
dic_all_values = {}

for folder in overall_folder:
    dic_all_values[folder] = []
    list_all_images = []
    folder_list = os.listdir(folder)
    for subfolder in folder_list:
        liste_images =  os.listdir(folder+'/'+subfolder)
        liste_images = [subfolder+'/'+ s  for s in liste_images]
        list_all_images.append(liste_images)
    for count, actin_image_name in enumerate(list_all_images[0]):
        actin_image = cv2.imread(folder+'/'+actin_image_name)[:,:,0]
        dapi_image = cv2.imread(folder+'/'+list_all_images[1][count])[:,:,0]
        yap_image = cv2.imread(folder+'/'+list_all_images[2][count])[:,:,0]
        ratio = get_ratio_yap(yap_image, dapi_image, actin_image)
        dic_all_values[folder].append(ratio)
        
for key in dic_all_values.keys():
    liste = dic_all_values[key]
    liste = np.array(liste)
    liste = liste[~np.isnan(liste)]
    liste = liste[liste < 1E308]
    print(key)
    print(np.mean(liste))
    print(np.median(liste))

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
sns.set_style("darkgrid")

over_all_list = []
percent = [5,7.5,10,12.5]
i = 0

for key in dic_all_values.keys():
    liste = dic_all_values[key]
    liste = np.array(liste)
    liste = liste[~np.isnan(liste)]
    liste = liste[liste < 8]
    for num in liste:
        over_all_list.append([percent[i],num])
    i += 1
    
df = pd.DataFrame(over_all_list,columns=['Percentage of PEG-NB Gel','Ratio (Intensity of Nucleus/Cytosol)'])
    
p = sns.violinplot(x='Percentage of PEG-NB Gel', y='Ratio (Intensity of Nucleus/Cytosol)', data=df)
l = sns.stripplot(x='Percentage of PEG-NB Gel', y='Ratio (Intensity of Nucleus/Cytosol)', data=df, jitter=True, color='black')
l.set_xlabel('Percentage of PEG-NB Gel', fontsize=30)
l.set_ylabel('Ratio (Intensity of Nucleus/Cytosol)', fontsize=30)
l.tick_params(axis='both', which='major', labelsize=30)
l.set_ylim(bottom=None, top=8)
l.set_title('Yap location based on intensity of antibody staining*', fontsize=50)
plt.show()


#%%

'with induvidiiuall threshold but median filter before'
import cv2
import numpy as np
import os
from scipy import ndimage

'With induvidual threshold'

def get_ratio_yap(yap_picture, dapi_image, actin_image):
    dapi_image = ndimage.median_filter(dapi_image, size=20)
    threshold_image_dapi = cv2.threshold(dapi_image, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)[1]
    boolean_matrix_dapi = np.array(threshold_image_dapi, dtype=bool)
    count_dapi = np.count_nonzero(boolean_matrix_dapi)
    if count_dapi < 50:
        return float('NaN')
    average_nucleus = np.average(yap_picture, weights = boolean_matrix_dapi)
    
    actin_image = ndimage.median_filter(actin_image, size=20)
    threshold_image_actin = cv2.threshold(actin_image, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)[1]
    boolean_matrix_actin = np.array(threshold_image_actin, dtype=bool)
    count_actin = np.count_nonzero(boolean_matrix_actin)
    if count_actin < 50:
        return float('NaN')
    
    not_nucleus = np.logical_not(boolean_matrix_dapi)
    boolean_matrix_cytosol = np.logical_and(boolean_matrix_actin,not_nucleus)
    count_cytosol = np.count_nonzero(boolean_matrix_cytosol)
    if count_cytosol < 50:
        return float('NaN')
    
    yap_picture = ndimage.median_filter(yap_picture, size=20)
    average_cytosol = np.average(yap_picture, weights = boolean_matrix_cytosol)
    
    ratio = average_nucleus / average_cytosol
    
    return ratio


overall_folder = os.listdir()
dic_all_values = {}

for folder in overall_folder:
    dic_all_values[folder] = []
    list_all_images = []
    folder_list = os.listdir(folder)
    for subfolder in folder_list:
        liste_images =  os.listdir(folder+'/'+subfolder)
        liste_images = [subfolder+'/'+ s  for s in liste_images]
        list_all_images.append(liste_images)
    for count, actin_image_name in enumerate(list_all_images[0]):
        actin_image = cv2.imread(folder+'/'+actin_image_name)[:,:,0]
        dapi_image = cv2.imread(folder+'/'+list_all_images[1][count])[:,:,0]
        yap_image = cv2.imread(folder+'/'+list_all_images[2][count])[:,:,0]
        ratio = get_ratio_yap(yap_image, dapi_image, actin_image)
        dic_all_values[folder].append(ratio)
        
for key in dic_all_values.keys():
    liste = dic_all_values[key]
    liste = np.array(liste)
    liste = liste[~np.isnan(liste)]
    liste = liste[liste < 1E308]
    print(key)
    print(np.mean(liste))
    print(np.median(liste))

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
sns.set_style("darkgrid")

over_all_list = []
percent = [5,7.5,10,12.5]
i = 0

for key in dic_all_values.keys():
    liste = dic_all_values[key]
    liste = np.array(liste)
    liste = liste[~np.isnan(liste)]
    liste = liste[liste < 8]
    for num in liste:
        over_all_list.append([percent[i],num])
    i += 1
    
df = pd.DataFrame(over_all_list,columns=['Percentage of PEG-NB Gel','Ratio (Intensity of Nucleus/Cytosol)'])
    
p = sns.violinplot(x='Percentage of PEG-NB Gel', y='Ratio (Intensity of Nucleus/Cytosol)', data=df)
l = sns.stripplot(x='Percentage of PEG-NB Gel', y='Ratio (Intensity of Nucleus/Cytosol)', data=df, jitter=True, color='black')
l.set_xlabel('Percentage of PEG-NB Gel', fontsize=30)
l.set_ylabel('Ratio (Intensity of Nucleus/Cytosol)', fontsize=30)
l.tick_params(axis='both', which='major', labelsize=30)
l.set_ylim(bottom=None, top=20)
l.set_title('Yap location based on intensity of antibody staining*', fontsize=50)
plt.show()

#%%
import cv2
'With induvidual threshold but local theshold'

def get_ratio_yap(yap_picture, dapi_image, actin_image):
    threshold_image_dapi = cv2.adaptiveThreshold(dapi_image, 255, cv2.ADAPTIVE_THRESH_MEAN_C,
                                          cv2.THRESH_BINARY, 199, 5)
    # threshold_image_dapi = cv2.threshold(dapi_image, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)[1]
    boolean_matrix_dapi = np.array(threshold_image_dapi, dtype=bool)
    count_dapi = np.count_nonzero(boolean_matrix_dapi)
    if count_dapi < 50:
        return float('NaN')
    average_nucleus = np.average(yap_picture, weights = boolean_matrix_dapi)
    
    threshold_image_actin = cv2.adaptiveThreshold(actin_image, 255, cv2.ADAPTIVE_THRESH_MEAN_C,
                                          cv2.THRESH_BINARY, 199, 5)
    # threshold_image_actin = cv2.threshold(actin_image, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU)[1]
    boolean_matrix_actin = np.array(threshold_image_actin, dtype=bool)
    count_actin = np.count_nonzero(boolean_matrix_actin)
    if count_actin < 50:
        return float('NaN')
    
    not_nucleus = np.logical_not(boolean_matrix_dapi)
    boolean_matrix_cytosol = np.logical_and(boolean_matrix_actin,not_nucleus)
    count_cytosol = np.count_nonzero(boolean_matrix_cytosol)
    if count_cytosol < 50:
        return float('NaN')
    average_cytosol = np.average(yap_picture, weights = boolean_matrix_cytosol)
    
    ratio = average_nucleus / average_cytosol
    
    return ratio


overall_folder = os.listdir()
dic_all_values = {}

for folder in overall_folder:
    dic_all_values[folder] = []
    list_all_images = []
    folder_list = os.listdir(folder)
    for subfolder in folder_list:
        liste_images =  os.listdir(folder+'/'+subfolder)
        liste_images = [subfolder+'/'+ s  for s in liste_images]
        list_all_images.append(liste_images)
    for count, actin_image_name in enumerate(list_all_images[0]):
        actin_image = cv2.imread(folder+'/'+actin_image_name)[:,:,0]
        dapi_image = cv2.imread(folder+'/'+list_all_images[1][count])[:,:,0]
        yap_image = cv2.imread(folder+'/'+list_all_images[2][count])[:,:,0]
        ratio = get_ratio_yap(yap_image, dapi_image, actin_image)
        dic_all_values[folder].append(ratio)
        
for key in dic_all_values.keys():
    liste = dic_all_values[key]
    liste = np.array(liste)
    liste = liste[~np.isnan(liste)]
    liste = liste[liste < 1E308]
    print(key)
    print(np.mean(liste))
    print(np.median(liste))

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
sns.set_style("darkgrid")

over_all_list = []
percent = [5,7.5,10,12.5]
i = 0

for key in dic_all_values.keys():
    liste = dic_all_values[key]
    liste = np.array(liste)
    liste = liste[~np.isnan(liste)]
    liste = liste[liste < 8]
    for num in liste:
        over_all_list.append([percent[i],num])
    i += 1
    
df = pd.DataFrame(over_all_list,columns=['Percentage of PEG-NB Gel','Ratio (Intensity of Nucleus/Cytosol)'])
    
p = sns.violinplot(x='Percentage of PEG-NB Gel', y='Ratio (Intensity of Nucleus/Cytosol)', data=df)
l = sns.stripplot(x='Percentage of PEG-NB Gel', y='Ratio (Intensity of Nucleus/Cytosol)', data=df, jitter=True, color='black')
l.set_xlabel('Percentage of PEG-NB Gel', fontsize=30)
l.set_ylabel('Ratio (Intensity of Nucleus/Cytosol)', fontsize=30)
l.tick_params(axis='both', which='major', labelsize=30)
l.set_ylim(bottom=None, top=15)
l.set_title('Yap location based on intensity of antibody staining*', fontsize=50)
plt.show()


"""














