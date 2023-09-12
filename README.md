# csst_simulation_catalog-hyz

*This is the first time that I pull codes to github. If there is any question, please confirmly let me know. As a beginner of coding, I will be appreciated if any advice could help.*

*Simulation catalog generation created by hyz, containing different ways to generate catalogs for CSST simulation observation.*

*The codes were not carefully revised, remaining a complete clumsy development process.*

> There are three catalogs used in these codes: [Gaia](https://gea.esac.esa.int/archive/), [CSST-Trilegal], simulation catalog provided by xxx.

## cat_extraction-galaxycenter
A combination method of Gaia and CSST-Trilegal catalog. Extractions from specific sky area are firstly done by healpix. Next combination simply replaces the brightest stars in Trilegal catalog by Gaia catalog, according to Gaia's limiting magnitude. Notice in this code, magnitude system has not been been converted, which will leading to a serious error in replacement. A better version is in [cat_extraction.ipynb](./cat_extraction.ipynb).
