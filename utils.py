import healpy as hp
from astropy.coordinates import SkyCoord
import numpy as np
from astropy.table import Table
import os
import astropy.units as u
import pandas as pd
import matplotlib.pyplot as plt

#####################################################
 # 检验
def readfits(id,root):
    f = f'triout_128_{id}.fits'
    t = Table.read(os.path.join(root,f))
    t = t.to_pandas()
    t = t[['gall','galb']]
    t = t[0::1000]
    _coor = SkyCoord(l=t['gall']*u.degree,b=t['galb']*u.degree,unit='deg',frame='galactic')
    _coor = _coor.transform_to('icrs')
    _ra, _dec = _coor.ra.degree, _coor.dec.degree 
    

    t['ra'] = _ra
    t['dec'] = _dec
    return t
    
def simulation_stractor_single(ipix_list,ra,dec,distance,root):
    target = pd.DataFrame()
    columns = ['M_H', 'logTe', 'logg', 'mu0',
            'NUVmag', 'umag', 'gmag', 'rmag', 'imag', 'zmag', 'ymag', 
            'uSmag', 'gSmag', 'rSmag', 'iSmag', 'zSmag', 
            'PMracos', 'PMdec', 'Vrad', 'ra', 'dec', 'gall', 'galb']
    for id in ipix_list:
        f = f'triout_128_{id}.fits'
        t = Table.read(os.path.join(root,f))
        t = t.to_pandas()
        _coor = SkyCoord(l=t['gall']*u.degree,b=t['galb']*u.degree,unit='deg',frame='galactic')
        _coor = _coor.transform_to('icrs')
        _ra, _dec = _coor.ra.degree, _coor.dec.degree
        t['ra'] = _ra
        t['dec'] = _dec
        # t['parallax'] = 10**(0.2*t['mu0']+1) * u.pc
        # t['parallax'] = 2*u.au/t['parallax'].to(u.au)
        t = t[(t['ra']-ra)**2+(t['dec']-dec)**2<=distance**2]
        # target = pd.concat([target,t])
        target = pd.concat([target,t[columns]])
    return target
    
def sim_stractor(ra,dec,radius,save_path):
    nside = 128

    # root path
    root = '/media/hyz/travaller/CSSTsim/CSSTsim.20211214.pSDSSpPS1'

    # turn icrs to gal
    cencoor = SkyCoord(ra=ra,dec=dec,unit='deg',frame='icrs')
    galc = cencoor.transform_to('galactic')
    gall, galb = galc.l.degree, galc.b.degree

    ipix_list = hp.query_disc(nside, hp.pixelfunc.ang2vec(gall,galb,lonlat=True), np.radians(radius), nest=True)

    # 如果范围太小，可能只有一个pix，换用单pix函数
    if len(ipix_list) == 0:
        ipix_list = [hp.ang2pix(nside,gall,galb,nest=True,lonlat=True)]
    print(ipix_list)
    
    # 检验
    print(readfits(ipix_list[0],root))    

    selection = simulation_stractor_single(ipix_list,ra,dec,radius,root)

    selection.to_csv(save_path[:-4]+'_sim'+save_path[-4:])

#####################################################
def gaia_stractor_single(ipix_list,ra,dec,distance):
        target = pd.DataFrame()
        columns = ['phot_g_mean_mag', 'phot_rp_mean_mag', 'phot_bp_mean_mag',
                'l', 'b', 'ra', 'ra_error', 'dec', 'dec_error',
                'pmra', 'pmra_error', 'pmdec', 'pmdec_error', 
                'parallax', 'parallax_error', 'radial_velocity', 'radial_velocity_error',
                'teff_gspphot', 'logg_gspphot', 'mh_gspphot']
        for id in ipix_list:
            t = pd.read_csv('/media/hyz/dwarfcave/data/gaiaDR3/ms/%s'%id,comment='#',usecols=columns)
            t = t[(t['ra']-ra)**2+(t['dec']-dec)**2<=distance**2]
            # target = pd.concat([target,t])
            target = pd.concat([target,t[columns]])
        return target

def gaia_stractor(ra,dec,radius,save_path):
    nside = 256

    target_pix = hp.query_disc(nside, hp.pixelfunc.ang2vec(ra,dec,lonlat=True), np.radians(radius), nest=True)

    if len(target_pix) == 0:
        target_pix = [hp.ang2pix(nside,ra,dec,nest=True,lonlat=True)]
    plt.figure(figsize=(3,3))
    plt.hist(target_pix)
    print(len(target_pix))
    
    target_file = []
    for pix in target_pix:
        for i in os.listdir('/media/hyz/dwarfcave/data/gaiaDR3/ms'):
            if i.endswith('.csv') and i not in target_file:
                if int(i.split('-')[1].split('.csv')[0]) >= pix and int(i.split('-')[0].split('_')[1]) <= pix:
                    target_file.append(i)
                    break
    print(target_file)

    selection_g = gaia_stractor_single(target_file,ra,dec,radius)
    selection_g.to_csv(save_path[:-4]+'_gaia'+save_path[-4:])
    print(len(selection_g))
    

#####################################################
def cat_merge_padding(save_path, J2000=False):
    sim = pd.read_csv(save_path[:-4]+'_sim'+save_path[-4:])
    if J2000:
        gaia = pd.read_csv(save_path[:-4]+'_gaia_j2000'+save_path[-4:])
    else:
        gaia = pd.read_csv(save_path[:-4]+'_gaia'+save_path[-4:])
    sim = sim.sort_values('gmag')
    gaia = gaia.sort_values('phot_g_mean_mag')
    sim.index = np.arange(len(sim))
    gaia.index = np.arange(len(gaia))
    sim_columns_raw = sim.columns
    print(len(sim), len(gaia))
    
    gi_sdss = sim['gSmag'] - sim['iSmag']
    ri_sdss = sim['rSmag'] - sim['iSmag']
    sim['gmag_gaia'] = -0.074189 - 0.51409*gi_sdss - 0.080607*gi_sdss**2 + 0.0016001*gi_sdss**3 + sim['gSmag']  # Gaia-SDSS12 relationship
    sim['rmag_gaia'] = -0.029869 - 1.1303*ri_sdss + sim['rSmag']
    sim['bmag_gaia'] = 0.1463 + 1.7244*ri_sdss - 1.1912*ri_sdss**2 + 0.22004*ri_sdss**3 + sim['rSmag']
    sim['pmra_gaia'] = sim['PMracos']
    sim['parallax_gaia'] = u.AU.to(u.pc)/np.power(10, 0.2*sim['mu0']+1)/(1/3600*0.01*u.degree.to(u.rad))
    sim['radial_velocity_gaia'] = sim['Vrad']
    
    ssim_list = ['ra_error_gaia','dec_error_gaia','pmra_error_gaia','pmdec_error_gaia','parallax_error_gaia','radial_velocity_error_gaia']
    num_sim, num_gaia = len(sim), len(gaia)
    for isim in ssim_list:
        _hist,_edge = np.histogram(gaia[isim[:-5]],bins=100,range=[gaia[isim[:-5]].min(),gaia[isim[:-5]].max()])
        sample = np.empty(1)
        for i in range(len(_hist)):
            sample = np.concatenate([sample, np.random.rand(int(num_sim*_hist[i]/num_gaia)) * (_edge[i+1]-_edge[i]) + _edge[i]], axis=0)
        if len(sample)<num_sim:
            sample = np.concatenate([sample,np.ones(num_sim-len(sample))*np.nan])
        np.random.shuffle(sample)
        print(len(sim),len(sample))
        sim[isim] = sample[:num_sim]
    
    _hist,_edge = np.histogram(gaia['phot_g_mean_mag'].dropna(),bins=50)
    truncate_mag = _edge[(_hist[1:]-_hist[:-1]).argmin()+1]
    sim.drop(sim[sim['gmag_gaia']<truncate_mag].index,inplace=True)
    
    columns = ['gmag_gaia','rmag_gaia','bmag_gaia',
               'gall_gaia','galb_gaia','ra_gaia','ra_error_gaia','dec_gaia','dec_error_gaia',
               'pmra_gaia','pmra_error_gaia','pmdec_gaia','pmdec_error_gaia','parallax_gaia','parallax_error_gaia','radial_velocity_gaia','radial_velocity_error_gaia']
    columns_gaia = ['phot_g_mean_mag','phot_rp_mean_mag','phot_bp_mean_mag',
                    'l','b','ra','ra_error','dec','dec_error',
                    'pmra','pmra_error','pmdec','pmdec_error','parallax','parallax_error','radial_velocity','radial_velocity_error']
    columns_sim = ['gmag_gaia','rmag_gaia','bmag_gaia',
                   'gall','galb','ra','ra_error_gaia','dec','dec_error_gaia',
                   'pmra_gaia','pmra_error_gaia','PMdec','pmdec_error_gaia','parallax_gaia','parallax_error_gaia','radial_velocity_gaia','radial_velocity_error_gaia']
    
    target = pd.DataFrame()
    for (_,e) in enumerate(columns):
        # print(_,e,columns_sim[_],columns_gaia[_])
        target[e] = np.concatenate([sim[columns_sim[_]].to_numpy(),gaia[columns_gaia[_]].to_numpy()])
        # print(sim[columns_sim[_]][::10000],gaia[columns_gaia[_]][::10000],target[e][::10000])
    for e in sim_columns_raw:
        target[e] = np.concatenate([sim[e].to_numpy(),np.ones(num_gaia)*np.nan])


    target.loc[:len(sim),'if_gaia_true'] = int(0)
    target.loc[len(sim):,'if_gaia_true'] = int(1)
    print('len(sim_truncated),len(gaia),len(target)',len(sim),len(gaia),len(target))

    def insert_na(t):
        for col in t.columns:
            if not col.endswith('_gaia'):
                continue
            data = t[col]
            num_na = data.isna().sum()
            idx_na = data[data.isna()].index
            if num_na == 0:
                continue
            if 'mag' in col:
                insert_data = t[['gmag_gaia','rmag_gaia','bmag_gaia']].mean(axis=1)[idx_na]
                insert_data_ = pd.Series(np.random.rand(insert_data.isna().sum())+t['gmag_gaia'].max())
                insert_data_.index = insert_data[insert_data.isna()].index
                insert_data = insert_data.fillna(insert_data_)
                insert_data = np.random.rand(num_na) * 2 + insert_data
            elif 'error' in col:
                insert_data = np.random.rand(num_na) * data.std() + data.max()
            else:
                insert_data = np.random.randn(num_na) + data.mean()
            insert_data = pd.Series(insert_data)
            insert_data.index = idx_na
            t[col] = data.fillna(insert_data)
        return t

    insert_na(target)
    # insert_na(sim)[::20,'gmag_gaia':]

    # sim[['gmag_gaia','rmag_gaia','bmag_gaia']][-20:]

    target.to_csv(save_path[:-4]+'_concat'+save_path[-4:])
    
    return target
    
    
##############################################################
def cat_merge_forsim(save_path, J2000=False):
    sim = pd.read_csv(save_path[:-4]+'_sim'+save_path[-4:])
    if J2000:
        gaia = pd.read_csv(save_path[:-4]+'_gaia_j2000'+save_path[-4:])
    else:
        gaia = pd.read_csv(save_path[:-4]+'_gaia'+save_path[-4:])
    sim = sim.sort_values('gmag')
    gaia = gaia.sort_values('phot_g_mean_mag')
    sim.index = np.arange(len(sim))
    gaia.index = np.arange(len(gaia))
    print(len(sim), len(gaia))
    sim['RA'] = sim['ra']
    sim['Dec'] = sim['dec']
    sim['app_sdss_g'] = sim['gSmag']
    sim['teff'] = np.power(10,sim['logTe'])
    sim['feh'] = sim['M_H']
    sim['grav'] = sim['logg']
    sim['pmra'] = sim['PMracos']
    sim['pmdec'] = sim['PMdec']
    sim['parallax'] = u.AU.to(u.pc)/np.power(10, 0.2*sim['mu0']+1)/(1/3600*0.01*u.degree.to(u.rad))
    sim['RV'] = sim['Vrad']
    
    gaia_bprp = gaia['phot_bp_mean_mag'] - gaia['phot_rp_mean_mag']
    gaia['sdss_g'] = gaia['phot_g_mean_mag']-(0.13518-0.46245*gaia_bprp-0.25171*gaia_bprp**2+0.021349*gaia_bprp**3)
    
    _hist,_edge = np.histogram(gaia['sdss_g'].dropna(),bins=50)
    truncate_mag = _edge[(_hist[1:]-_hist[:-1]).argmin()+1]
    sim.drop(sim[sim['app_sdss_g']<truncate_mag].index,inplace=True)
    # num_gaia = len(gaia)
    # sim.loc[:num_gaia,'RA'] = gaia.loc[:num_gaia,'ra']
    # sim.loc[:num_gaia,'Dec'] = gaia.loc[:num_gaia,'dec']
    # sim.loc[:num_gaia,'g'] = gaia.loc[:num_gaia,'phot_g_mean_mag']
    # sim.loc[:num_gaia,'teff'] = gaia.loc[:num_gaia,'teff_gspphot']
    # sim.loc[:num_gaia,'grav'] = gaia.loc[:num_gaia,'logg_gspphot']
    # sim.loc[:num_gaia,'feh'] = gaia.loc[:num_gaia,'mh_gspphot']
    # sim.loc[:num_gaia,'pmra'] = gaia.loc[:num_gaia,'pmra']
    # sim.loc[:num_gaia,'pmdec'] = gaia.loc[:num_gaia,'pmdec']
    # sim.loc[:num_gaia,'parallax'] = gaia.loc[:num_gaia,'parallax']
    # sim.loc[:num_gaia,'rv'] = gaia.loc[:num_gaia,'radial_velocity']
    
    columns = ['RA','Dec','app_sdss_g',
               'teff','grav','feh',
               'pmra','pmdec','parallax','RV']
    columns_gaia = ['ra','dec','sdss_g',
                    'teff_gspphot','logg_gspphot','mh_gspphot',
                    'pmra','pmdec','parallax','radial_velocity']
    target = pd.DataFrame()
    for (_,e) in enumerate(columns):
        target[e] = np.concatenate([sim[e].to_numpy(),gaia[columns_gaia[_]].to_numpy()])
    
    # sim = sim[columns]
    
    # sim.to_csv(save_path[:-4]+'_merge_forsim'+save_path[-4:])
    target.to_csv(save_path[:-4]+'_merge_forsim'+save_path[-4:])
    return target




##############################################################
import numpy as np
from astropy.time import Time
from astropy.table import Table, hstack
from astropy.coordinates import SkyCoord,Distance 
import astropy.units as u

def epoch_trans(save_path,t0,t1):
    data_c = Table.read(save_path[:-4]+'_gaia'+save_path[-4:])
    #天体的基本信息

    Alpha_t1 = np.array(data_c['ra'])
    Delta_t1 = np.array(data_c['dec'])
    pmra1 = np.array(data_c['pmra'])
    pmdec1 = np.array(data_c['pmdec']) 
    para1 = np.array(data_c['parallax'])
    radial_velocity1 = np.zeros(len(data_c))

    ra_2000 = []
    dec_2000 = []

    c = SkyCoord(ra=Alpha_t1 * u.deg,
                dec=Delta_t1 * u.deg,
                distance=Distance(parallax=abs(para1) * u.mas),
                pm_ra_cosdec=pmra1 * u.mas/u.yr,   #np.cos(np.radians(table['dec'])
                radial_velocity=radial_velocity1*u.km/u.s,
                pm_dec=pmdec1 * u.mas/u.yr,
                obstime=Time(2016, format='jyear',
                            scale='tcb'),  frame="icrs" )


    epoch_now = Time(t0, format='jyear', scale='tcb')
    obstime=Time(t1, format='jyear', scale='tcb')
    c_epoch_now = c.apply_space_motion(epoch_now)

    ra_2000 = c_epoch_now.ra.degree
    dec_2000 = c_epoch_now.dec.degree
        

    RA_l = Table([ra_2000],names='q')
    DEC_l = Table([dec_2000],names='w')

    RA_l['q'].name='ra_j2000'
    DEC_l['w'].name='dec_j2000'

    All_data = hstack([data_c, RA_l, DEC_l])

    All_data.write(save_path[:-4]+'_gaia_j2000'+save_path[-4:],format="ascii.csv",overwrite=True) 
