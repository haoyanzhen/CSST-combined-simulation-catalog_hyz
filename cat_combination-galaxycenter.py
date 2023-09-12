import pandas as pd
import numpy as np
import os
from astropy.coordinates import SkyCoord
from astropy.table import Table
import astropy.units as u
import healpy as hp
import cat2hdf5

# 计算银心位置对应的仿真和Gaia星表文件
radius = 2
nside = 256
gall,galb = 0,0
print(gall,galb)
cencoor = SkyCoord(l=gall,b=galb,unit='deg',frame='galactic')
icrsc = cencoor.transform_to('icrs')
ra, dec = icrsc.ra.degree, icrsc.dec.degree
print(ra,dec)
# Get the list of pixels within the radius
ipix_list = hp.query_disc(nside, hp.pixelfunc.ang2vec(ra,dec,lonlat=True), np.radians(radius), nest=True)

# 如果范围太小，可能只有一个pix，换用单pix函数
if len(ipix_list) == 0:
    ipix_list = [hp.ang2pix(nside,gall,galb,nest=True,lonlat=True)]
ipix_gaia_l,ipix_gaia_r = min(ipix_list), max(ipix_list)

# 从gaia星表里得到包含这些区域的文件
gaia_root = '/home/haoyanzhen_shao/project/simulation_work/data_dir/catalog/CSST_trilegal/'
gaia_files = []
for gf in os.listdir(gaia_root):
    if gf.endswith('.csv'):
        left = int(gf.split('_')[1].split('-')[0])
        right = int(gf.split('-')[1].split('.c')[0])
        if left <= ipix_gaia_r and right >= ipix_gaia_l:
            gaia_files.append(gf)
len(gaia_files)

# 提取Gaia星表
def gaia_stractor(gaia_files,ra,dec,radius):
    selection = pd.DataFrame()
    columns = ['RA','Dec','app_sdss_g','teff','grav','feh','pmra','pmdec','RV','parallax','label']
    for file in gaia_files:
        path = os.path.join(gaia_root,file)
        cat = pd.read_csv(path,comment='#',usecols=['ra','dec','phot_g_mean_mag','mh_gspphot','logg_gspphot','teff_gspphot','pmra','pmdec','radial_velocity','parallax'])
        cat = cat[(cat['ra']-ra)**2+(cat['dec']-dec)**2<=radius**2]
        cat['RA'] = cat['ra']
        cat['Dec'] = cat['dec']
        cat['app_sdss_g'] = cat['phot_g_mean_mag']
        cat['teff'] = cat['teff_gspphot']
        cat['grav'] = cat['logg_gspphot']
        cat['feh'] = cat['mh_gspphot']
        cat['RV'] = cat['radial_velocity']
        cat['label'] = 'Gaia'
        selection = pd.concat([selection,cat[columns]])
    return selection

catgaia = gaia_stractor(gaia_files,ra,dec,radius)

# 提取仿真星表
nside = 128

root = '/home/haoyanzhen_shao/project/simulation_work/data_dir/catalog/CSST_trilegal/'
ipix_list = hp.query_disc(nside, hp.pixelfunc.ang2vec(gall,galb,lonlat=True), np.radians(radius), nest=True)

# 如果范围太小，可能只有一个pix，换用单pix函数
if len(ipix_list) == 0:
    ipix_list = [hp.ang2pix(nside,gall,galb,nest=True,lonlat=True)]
ipix_list

def simulation_stractor(ipix_list,ra,dec,radius):
    target = pd.DataFrame()
    columns = ['RA','Dec','app_sdss_g','teff','grav','feh','pmra','pmdec','RV','parallax','label']
    for id in ipix_list:
        f = f'triout_128_{id}.fits'
        t = Table.read(os.path.join(root,f))
        t = t.to_pandas()
        # print(t.columns)
        # print(t['gall'])
        # print(t['galb'])
        _coor = SkyCoord(l=t['gall'],b=t['galb'],unit='deg',frame='galactic')
        # print('probe')
        _coor = _coor.transform_to('icrs')
        _ra, _dec = _coor.ra.degree, _coor.dec.degree
        t['RA'] = _ra
        t['Dec'] = _dec
        print('ra,dec:\n',t['RA'][::10000],t['Dec'][::10000])
        t = t[(t['RA']-ra)**2+(t['Dec']-dec)**2<=radius**2]
        t['app_sdss_g'] = t['gSmag']
        t['teff'] = 10**t['logTe']
        t['grav'] = t['logg']
        t['feh'] = t['M_H']
        t['pmra'] = t['PMracos']/np.cos(t['Dec'])
        t['pmdec'] = t['PMdec']
        t['RV'] = t['Vrad']
        t['parallax'] = 2*u.au/(10**(0.2*t['mu0']+1) * (1*u.pc).to(u.au).value)
        print('parallax',t['parallax'][::10000])
        t['label'] = 'sim'
        target = pd.concat([target,t[columns]])
    return target

catsim = simulation_stractor(ipix_list,ra,dec,radius)


# 将仿真星表中较亮的星用Gaia完全代替
catsim = catsim.loc[catsim['app_sdss_g']>=catgaia['app_sdss_g'].max(),:]
concat = pd.concat([catsim,catgaia])
concat = Table.from_pandas(concat)

# 保存文件
save_dir = '/home/haoyanzhen_shao/project/simulation_work/data_dir/catalog/'
concat.write(root+'galactic_center_20230423.fits', format='fits',overwrite=True)
cat2hdf5.fits_to_hdf5(root+'galactic_center_20230423.fits',root+'galactic_center_20230423.hdf5')
print('test succ')