import healpy as hp
import pandas as pd
import numpy as np
import os


nside = 256
ra = 249.7119627755216
dec = -44.60816045259986
radius = 2

# Get the list of pixels within the radius
target_pix = hp.query_disc(nside, hp.pixelfunc.ang2vec(ra,dec,lonlat=True), np.radians(radius), nest=True)

# 如果范围太小，可能只有一个pix，换用单pix函数
if len(target_pix) == 0:
    target_pix = [hp.ang2pix(nside,ra,dec,nest=True,lonlat=True)]
print(target_pix)

target_file = []
for pix in target_pix:
    for i in os.listdir('/media/hyz/dwarfcave/data/gaiaDR3/ms'):
        if i.endswith('.csv'):
            if int(i.split('-')[1].split('.csv')[0]) >= pix and int(i.split('-')[0].split('_')[1]) <= pix:
                target_file.append(i)
                break
target_file = np.array(target_file)
target_file = np.unique(target_file)
print(target_file)

gaiatb = pd.DataFrame()
for f in target_file:
    newtb = pd.read_csv('/media/hyz/dwarfcave/data/gaiaDR3/ms/%s'%f,comment='#')
    gaiatb = pd.concat((gaiatb,newtb),axis=0)
selection_g = gaiatb[(gaiatb['ra']-ra)**2+(gaiatb['dec']-dec)**2<=radius**2]
selection_g.to_csv('./selection20230525_gaia.csv')
print(selection_g[:1000:100])