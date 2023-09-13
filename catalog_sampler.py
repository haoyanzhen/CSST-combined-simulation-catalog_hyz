from astropy.table import Table
import pandas as pd
import numpy as np
from mpi4py import MPI
import os

comm = MPI.COMM_WORLD
csize = comm.size
crank = comm.rank

cat_path = '/home/haoyanzhen_shao/project/simulation_work/data_dir/reference_catalog/CT1.5-N6397-GDR3-G26.0-rev_J2000.fits'
cat_raw = Table.read(cat_path)
cat_raw = cat_raw.to_pandas()
cat_raw = cat_raw[['RA','Dec','grav','teff','feh','app_sdss_g','pmra','pmdec','parallax','RV']]
# name_selection = ['Dist','grav','teff','feh','app_Gaia_G']

name_all = ['grav','teff','feh','app_sdss_g','pmra','pmdec','parallax','RV']
cat_raw_desc = cat_raw[name_all].describe()

stamp_nstar = '5e6'
nstar = int(5e6)
nstar = int(nstar/csize)
ra = np.random.rand(nstar) * 1.4 + 270.3
dec = np.random.rand(nstar) * 1.4 - 41.7

_hist,_edges = np.histogram(cat_raw['RV'],bins=90,range=(-280,440))
np.save('/home/haoyanzhen_shao/project/simulation_work/data_dir/catalog/crowdedField_sampler0627/histgram_rv.npy',np.array([_hist,_edges],dtype=object))
_hist = _hist/len(cat_raw['RV'])
_nedges = _edges + (_edges[1]-_edges[0])/2
_nedges = _nedges[:-1]
from scipy.interpolate import CubicSpline
interp = CubicSpline(_nedges,_hist)
a = np.random.rand(nstar)*(438+279) - 279
p = interp(a)
p[p<0] = 0
p = p/np.sum(p)
rv = np.random.choice(a,size=nstar,p=p)
del _hist,_edges,interp,_nedges,a,p

def stat(df):
    df_selection_len = len(df.columns)
    df_len = len(df)
    df_desc = df.describe()
    bins = np.array([np.arange(-1.6,0.9,0.1)**3+4.8,
                    np.arange(2000,8700,670),
                    np.arange(-3.72,0.54,0.426),
                    np.arange(-2.7,1.4,0.1)**3+25],dtype=object)  # 5-27.8
    range_ = np.array([df_desc.min()[:],df_desc.iloc[3:].max()[:]]).T
    df_hist,_edges = np.histogramdd(df.to_numpy(),bins,range_)
    np.save('/home/haoyanzhen_shao/project/simulation_work/data_dir/catalog/crowdedField_sampler0627/histgram_phy.npy',np.array([df_hist,_edges],dtype=object))

    # hist调整
    _lens = []
    for i in range(4):
        _lens.append(_edges[i][1:] - _edges[i][:-1])
        _edges[i] = (_edges[i][:-1]+_edges[i][1:])/2
    df_hist = df_hist/_lens[0][:,None,None,None]/_lens[1][None,:,None,None]/_lens[2][None,None,:,None]/_lens[3][None,None,None,:]
    df_hist = df_hist/df_hist.sum()

    # 概率分布函数插值
    from scipy.interpolate import RegularGridInterpolator

    interp = RegularGridInterpolator(_edges,df_hist,method='cubic')
    return interp, _edges

# def stat2(df):
#     df_len = len(df)
#     df_desc = df.describe()
#     bins = np.array([np.arange(-8,8,0.5)**3,
#                     np.arange(-8,8,0.5)**3,
#                     np.float_power(10,np.arange(-8,2,0.25))],dtype=object)
#     range_ = np.array([df_desc.min()[:],df_desc.iloc[3:].max()[:]]).T
#     df_hist,_edges = np.histogramdd(df.to_numpy(),bins,range_)
#     np.save('/home/haoyanzhen_shao/project/simulation_work/data_dir/catalog/crowdedField_sampler0627/histgram_ast.npy',np.array([df_hist,_edges],dtype=object))

#     # hist调整
#     _lens = []
#     for i in range(3):
#         _lens.append(_edges[i][1:] - _edges[i][:-1])
#         # _edges[i] = (_edges[i][:-1]+_edges[i][1:])/2
#         _edges[i] = _edges[i][:-1]
#     df_hist = df_hist/_lens[0][:,None,None]/_lens[1][None,:,None]/_lens[2][None,None,:]
#     df_hist = df_hist/df_hist.sum()

#     # 概率分布函数插值
#     from scipy.interpolate import RegularGridInterpolator

#     interp = RegularGridInterpolator(_edges,df_hist,method='nearest')
#     return interp, _edges

name_groups = []
name_groups.append(name_all[:4])
name_groups.append(name_all[4:-1])
bins = np.array([40,20,20,100,20,20,20])
interp_phy, _edge_phy = stat(cat_raw[name_groups[0]])
# interp_ast, _edge_ast = stat2(cat_raw[name_groups[1]])

def starcount_phy_log(x):
    ndim = 4
    try:
        p = interp_phy(x)
        if p <= 0:
            p = 1e-20
        return np.log10(p)
    except Exception:
        p = 0
        for i in range(ndim):
            if x[i] > _edge_phy[i][-1]:
                p += -30 - 30*(x[i]-_edge_phy[i][-1])/(_edge_phy[i][-1]-_edge_phy[i][0])
            elif x[i] < _edge_phy[i][0]:
                p += -30 - 30*(_edge_phy[i][0]-x[i])/(_edge_phy[i][-1]-_edge_phy[i][0])
        if p > -10:
            p -= 20
        return p

# def starcount_ast_log(x):
#     ndim = 3
#     try:
#         p = interp_ast(x)
#         if p <= 0:
#             p = 1e-20
#         return np.log10(p)
#     except Exception:
#         p = 0
#         for i in range(ndim):
#             if x[i] > _edge_ast[i][-1]:
#                 p += -30 - 100*(x[i]-_edge_ast[i][-1])/(_edge_ast[i][-1]-_edge_ast[i][0])
#             elif x[i] < _edge_ast[i][0]:
#                 p += -30 - 100*(_edge_ast[i][0]-x[i])/(_edge_ast[i][-1]-_edge_ast[i][0])
#         if p > -20:
#             p -= 20
#         return p
    
    
from emcee import EnsembleSampler

nwalkers = 100
ndim = 4
psteps = 100
nsteps = int(nstar/nwalkers)
# x[grav,teff,feh,mag_g]
# [ 3.76236096e-01,  5.35259056e+00],
# [ 6.57588529e+02,  8.69925293e+03],
# [-3.72083282e+00,  5.36125898e-01],
# [ 1.68997382e+00,  2.74582329e+01]
x0 = np.empty((nwalkers,ndim))
x0[:,0] = np.random.rand(nwalkers) * 1.6 + 3.7
x0[:,1] = np.random.rand(nwalkers) * 8000 + 660
x0[:,2] = np.random.rand(nwalkers) * 4 - 3.7
x0[:,3] = np.random.rand(nwalkers) * 25 + 2

sampler = EnsembleSampler(nwalkers,ndim,starcount_phy_log)
x = sampler.run_mcmc(x0, psteps, progress=True)
sampler.reset()
sampler.run_mcmc(x, nsteps, progress=True)
samples_phy = sampler.get_chain().reshape((-1,ndim))

# nwalkers = 100
# ndim = 3
# psteps = 1000
# nsteps = int(nstar/nwalkers)
# # x[mag_g,pmra,pmdec,parallax]
# # [ 5.21828318e+00,  2.74582329e+01],
# # [-3.84212126e+02,  1.68569655e+02],
# # [-4.81517687e+02,  1.01266962e+02],
# # [ 8.70646100e-04,  7.40291565e+01]
# x0 = np.empty((nwalkers,ndim))
# x0[:,0] = np.random.rand(nwalkers) * 500 - 380
# x0[:,1] = np.random.rand(nwalkers) * 550 - 480
# x0[:,2] = np.random.rand(nwalkers) * 70 + 1

# # 不确定mag对pm和pa产生影响的原因，可能是基于星系动力学，由于采样较难，目前假设为均匀分布情况下无影响。
# # x_mag = samples_phy[:,-1][:,None]

# sampler = EnsembleSampler(nwalkers,ndim,starcount_ast_log)
# sampler.reset()
# # samplel = []
# # for i in range(nsamples):
# #     x0 = sampler.sample(x0)
# #     samplel.append(x0)
# x = sampler.run_mcmc(x0, psteps, progress=True)
# sampler.reset()
# sampler.run_mcmc(x, nsteps, progress=True)
# samples_ast = sampler.get_chain().reshape((-1,ndim))

grav = samples_phy[:,0].copy()
teff = samples_phy[:,1].copy()
feh = samples_phy[:,2].copy()
gaia_g = samples_phy[:,3].copy()
# pmra = samples_ast[:,0]
# pmdec = samples_ast[:,1]
# parallax = samples_ast[:,2]
del interp_phy,_edge_phy,sampler,samples_phy

from scipy.optimize import curve_fit

def norm(x,mu,var):
    return 1/((2*np.pi)**0.5 * var) * np.exp(-(x-mu)**2/var**2)

x = cat_raw['pmra']
bins = np.arange(-400,400,0.01)
_hist, _edge = np.histogram(x,bins)
# popt,pcov,infodict,mesg,flg = curve_fit(norm,_edge[:-1],_hist)
popt,pcov = curve_fit(norm,_edge[:-1],_hist)
pmra = np.random.normal(loc=popt[0],scale=popt[1],size=nstar)

x = cat_raw['pmdec']
_hist, _edge = np.histogram(x,bins)
popt,pcov = curve_fit(norm,_edge[:-1],_hist)
pmdec = np.random.normal(loc=popt[0],scale=popt[1],size=nstar)

def exp(x,k,e):
    return k * 2 ** (-e*x)
x = cat_raw['parallax']
bins = np.arange(0,40,0.01)
_hist, _edge = np.histogram(x,bins)
popt,pcov,info,mesg,flg= curve_fit(exp,_edge[:-1],_hist,full_output=True)
a = np.random.rand(nstar)*10
p = exp(a,popt[0],popt[1])
p /= p.sum()
parallax = np.random.choice(a,size=nstar,p=p)
del a,p

cat_len = nstar
for _data in [ra,dec,gaia_g,teff,grav,feh,pmra,pmdec,rv,parallax]:
    cat_len = (len(_data) if len(_data)<cat_len else cat_len)
for _data in [ra,dec,gaia_g,teff,grav,feh,pmra,pmdec,rv,parallax]:
    _data = _data[:cat_len]

cat_new = Table([ra,dec,gaia_g,teff,grav,feh,pmra,pmdec,rv,parallax],
                names=('RA','Dec','app_sdss_g','teff','grav','feh','pmra','pmdec','RV','parallax'))
cat_new.write('/home/haoyanzhen_shao/project/simulation_work/data_dir/catalog/crowdedField_sampler0627/from_N6397_%s_%s.fits'%(stamp_nstar,crank), format='fits',overwrite=True)

ssss=0
comm.gather(ssss,root=0)
if crank != 0:
    pass
else:
    from astropy.table import vstack
    for i in range(csize-1):
        i+=1
        cati = Table.read('/home/haoyanzhen_shao/project/simulation_work/data_dir/catalog/crowdedField_sampler0627/from_N6397_%s_%s.fits'%(stamp_nstar,i))
        cat_new = vstack((cat_new,cati))
        del cati
    cat_new = cat_new.to_pandas()
    cat_new = cat_new.drop(cat_new[(cat_new['grav']<1.3)|(cat_new['grav']>5.5)|
                                   (cat_new['teff']<2000)|(cat_new['teff']>10000)|
                                   (cat_new['feh']<-4)|(cat_new['feh']>0.5)].index)
    cat_new = cat_new.dropna(axis=0)
    cat_new = Table.from_pandas(cat_new)
    cat_new.write('/home/haoyanzhen_shao/project/simulation_work/data_dir/catalog/crowdedField_sampler0627/from_N6397_%s.fits'%stamp_nstar, format='fits',overwrite=True)

    from self_test.cat2hdf5 import fits_to_hdf5 as f2h
    f2h("/home/haoyanzhen_shao/project/simulation_work/data_dir/catalog/crowdedField_sampler0627/from_N6397_%s.fits"%stamp_nstar,
        "/home/haoyanzhen_shao/project/simulation_work/data_dir/catalog/crowdedField_sampler0627/from_N6397_%s.hdf5"%stamp_nstar)
