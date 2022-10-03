# Analyzing microgel and cell fluorescent images
- - - 

At the moment there are two scripts written to analyze the YAP content of confocal images to see the distribution of YAP content between different cells. Both are based on image with three channels (actin stained channel, Nuclei stained channel (DAPI) and YAP stained channel. The input images are separated images from a confocal microsope Z-Stack .

## Script 1
In this approach the images were analyzed per slide. The main principle is the get the ratio between the YAP intensity in the nucleus over the YAP intensity in the cytosol per image. Nucleus is defined as the region were there is DAPI signal, cytosol is defined as the region where there is actin signal but no DAPI signal.

More details are following.

## Script 2

In the second approach the YAP intensity of only the nuclei was measured (again defined by de DAPI signal). This signal was then normalized by the DAPI signal.


---
All rights 


