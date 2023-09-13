# csst_simulation_catalog-hyz

*This is the first time that I pull codes to github. If there is any question, please confirmly let me know. As a beginner of coding, I will be appreciated if any advice could help.*

*Simulation catalog generation created by hyz, containing different ways to generate catalogs for CSST simulation observation.*

*The codes were not carefully revised, remaining a complete clumsy development process.*

> There are three catalogs used in these codes: [Gaia](https://gea.esac.esa.int/archive/), [CSST-Trilegal](https://nadc.china-vo.org/data/data/csst-trilegal/f), simulation catalog provided by Shao Zhengyi, [CosmoDC2](https://data.lsstdesc.org/doc/cosmodc2).

## cat_combination.ipynb
A serious test and method of combination. The chapter 银心仿真数据生成 refers to [cat_combination-galaxycenter.py](./cat_combination-galaxycenter.py), and the chapter 优化合并星表 is the latest version. [utils.py](./utils.py) is needed in it. The latest version containing:
1. Extraction catalog from Trilegal and Gaia within a certain area. Parameters includes astrometry, radial velocity and three atmospheric parameters.
2. Transferming the epoch of Gaia from 2016 to 2000 for matching.
3. Importing [magnitude conversion system](https://gea.esac.esa.int/archive/documentation/GDR3/Data_processing/chap_cu5pho/cu5pho_sec_photSystem/cu5pho_ssec_photRelations.html). Replace the brightest stars in Trilegal with Gaia. 
4. Combining Gaia catalog and Trilegal catalog. Several changes have been made for requirement of CSST simulation:
    1. Empty errors of astrometry parameters and radial velocity are sampled based on attribution of Gaia.
    2. Nan value in Gaia is inserted by a formula:
       ```python
       insert_data = np.random.rand(num_na) * data.std() + data.max()
       ```
    3. Names of columns are changed to CSST-simulation format.

## cat_combination-galaxycenter.py
An early version of [cat_combination.ipynb](./cat_combination.ipynb).

## catalog_generator.ipynb
Contents:
- Early unsuccessful trials.
- A simple method assumes that parameters of stars follow simple distribution laws.
- A procedure from raw data to pointing file, stellar catalog and galaxy catalog in format of csst-simulation.

## catalog_sampler.ipynb
A method of generating catalog based on statistic and sampling. Template catalog is provided by Shao Zhengyi.

## catalog_sampler.py
A .py version of [catalog_sampler.ipynb](./catalog_sampler.py).

## utils.ipynb
Wraps some functions used in [cat_combination.ipynb](./cat_combination.ipynb).

## healpy_learning-Trilegal_statistics.ipynb
A basically learning code with healpy and a simple statistical analyze of Trilegal catalog.
