# csst_simulation_catalog-hyz

*This is the first time that I pull codes to github. If there is any question, please confirmly let me know. As a beginner of coding, I will be appreciated if any advice could help.*

*Simulation catalog generation created by hyz, containing different ways to generate catalogs for CSST simulation observation.*

*The codes were not carefully revised, remaining a complete clumsy development process.*

> There are three catalogs used in these codes: [Gaia](https://gea.esac.esa.int/archive/), [CSST-Trilegal](https://nadc.china-vo.org/data/data/csst-trilegal/f), simulation catalog provided by Shao Zhengyi.

## cat_extraction-galaxycenter.py
A combination method of Gaia and CSST-Trilegal catalog. Extractions from specific sky area are firstly done by healpix. Next combination simply replaces the brightest stars in Trilegal catalog by Gaia catalog, according to Gaia's limiting magnitude. Notice in this code, magnitude system has not been been converted, which will leading to a serious error in replacement. A better version is at the end of [cat_extraction.ipynb](./cat_extraction.ipynb).


## cat_combination.ipynb
A serious test and method of combination. The chapter 银心仿真数据生成 refers to [cat_combination-galaxycenter.py](./cat_combination-galaxycenter.py), and the chapter 优化合并星表 is the latest version. The latest version containing:
1. Extraction catalog from Trilegal and Gaia within a certain area. Parameters includes astrometry, radial velocity and three atmospheric parameters.
2. Transferming the epoch of Gaia from 2016 to 2000 for matching.
3. Importing [magnitude conversion system](https://gea.esac.esa.int/archive/documentation/GDR3/Data_processing/chap_cu5pho/cu5pho_sec_photSystem/cu5pho_ssec_photRelations.html).
4. Combining Gaia catalog and Trilegal catalog. Several changes have been made for requirement of CSST simulation:
    1. Empty errors of astrometry parameters and radial velocity are sampled based on attribution of Gaia.
    2. Nan value in Gaia is inserted by a formula:
       ```python
       insert_data = np.random.rand(num_na) * data.std() + data.max()
       ```
    3. Names of columns are changed to CSST-simulation format.
