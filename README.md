# README

Author: Andrew Eden

E-Mail: aeden2019@my.fit.edu

---

# Contents

* [Summary](#summary)
* [Requirements](#requirements)
    * [Gizmo Analysis](#requirements-gizmo-analysis)
    * [Py Ananke](#requirements-py-ananke)
    * [Nba & Pygadgetreader](#requirements-nba-pygadgetreader)
* [Ananke Parameters](#ananke-parameters)
    * [Input Parameters](#ananke-parameters-input) 
    * [Output Parameters](#ananke-parameters-output)
* [Configuraton Files](#configuration-files)
* [Future Work](#future-work)

---

## Summary <a name="summary"></a>

This code utilizes the [py-ananke](https://arxiv.org/abs/2312.02268) library to generate synthetic star catalogs from FIRE simulated data. In particular, it applies the constraints of the Vera Rubin telescope to generate a mock catalog, and then provides visualization via mollweide plots. 

To obtain the mock catalog, which will be stored in `survey.sim.h5`: 

```Bash
python src/mock/make_mock.py
```

To obtain a mollweide plot visualzation, which will be stored in `ananke_mollweide_plot.png`:

```Bash
python src/visualization/generate_plot_from_catalog.py
```

Key parameters can be modified in `src/mock/config.yaml`. 

You can also select a different yaml file when you the python script in the command line:

```Bash
python src/mock/make_mock.py --param foo.yaml
```

Additional pre-made configuration files used in notebooks can be found in `src/mock/`, such as `m12b_600_inner.yaml`. For more information about these files, refer to the Configuration Files section of the README. 


---

## Requirements <a name="requirements"></a>

General requirements (installable through pip):
* python
* numpy
* astropy
* numpy
* vaex
* args
* healpy

Specific requirements (see explanations below): 
* gizmo_analysis
* ananke
* nba
* pygadgetreader

For an easy installation of these packages, you can also use the provided shell script as follows:
```Bash 
bash dependencies.sh
```

### Gizmo Analysis <a name="requirements-gizmo-analysis"></a> 

This module is part of the processing for FIRE. Written by Andrew Wetzel, it can be found here: [https://bitbucket.org/awetzel/gizmo_analysis/src/master/](https://bitbucket.org/awetzel/gizmo_analysis/src/master/). Note that it also requires the [\utilities](https://bitbucket.org/awetzel/utilities/src/master/) package.  

The easiest installation is as follows:

```Bash 
git clone https://bitbucket.org/awetzel/gizmo_analysis.git
bash ./gizmo_analysis/install_helper.sh
```

### Py Ananke <a name="requirements-py-ananke"></a> 

This module is foundational to our code.  Written by Adrien Thob, it can be found here: [https://github.com/athob/py-ananke](https://github.com/athob/py-ananke). 

One of the possible installation methods is as follows:
```Bash
git clone https://github.com/athob/py-ananke
cd py-ananke
pip install .

```

### Nba & Pygadgetreader <a name="requirements-nba-pygadgetreader"></a> 

The modules are used for the visualization process. They can be cloned and installed from the github of [Nicolas Garavito Camargo](https://github.com/jngaravitoc) as follows:

```Bash
git clone https://github.com/jngaravitoc/nba
cd nba
pip install .
```

```Bash 
git clone https://github.com/jngaravitoc/pygadgetreader
cd pygadgetreader
pip install .
```

---

## Ananke Parameters <a name="ananke-parameters"></a> 

In this section you can find details regarding the input and ouput parameters for `py-ananke`. For additional information, refer to the `testing_ananke.ipynb` notebook located at [https://github.com/athob/py-ananke/blob/main/jupyter/testing_ananke.ipynb](https://github.com/athob/py-ananke/blob/main/jupyter/testing_ananke.ipynb). 

### Input Parameters <a name="ananke-parameters-input"></a> 

Below is a table of the parameters used with `py-ananke` to best simulate the characteristics of the Vera Rubin telescope. These values can be changes in `src/mock/make_mock.py` during the py-ananke process `ananke = an.Ananke()`.

| Parameter      | Value Used       | Description                                                   |
| -------------- | ---------------- | ------------------------------------------                    |
| particles      | p                | Dictionary of particle data from FIRE.                        |
| photo_sys      | 'padova/LSST'    | Photometric system Galaxia should use to generate the survey. |
| cmd_magnames   | ‘rmag,gmag-rmag’ | Names of the filters Galaxia should use for the color-magnitude diagram box selection. First argument is used for appMagLimits and the second for colorLimits.                             |
| app_mag_lim_lo | 17               | Lower limit in apparent mag. Related to cmd_magnames.         |
| app_mag_lim_hi | 27.5             | Upper limit in apparent mag. Related to cmd_magnames.         |
| abs_mag_lim_lo | -7.0             | Lower limit in absolute mag.                                  |
| abs_mag_lim_hi | 10.0             | Upper limit in absolute mag.                                  |

### Output Parameters <a name="ananke-parameters-output"></a> 

Below is a table of the parameters outputed by `py-ananke`. 

| Ouput Parameter | Description                   | Unit          |
| --------------- | ----------------------------- | ------------- |
| age             | Age                           | log (age/yr)  |
| alpha           | Alpha abundance               | [alpha/Fe]    |
| dec             | Declination                   | degree        |
| dmod            | Distance modulus              | NA            |
| feh             | Metallicity                   | [Fe/H]        |
| glat            | Galactic latitude             | degree        |
| glon            | Galactic longitude            | degree        |
| grav            | Surface gravity               | log(gravity)  |
| lsst_gmag       | Magnitude in g band           | mag           |
| lsst_imag       | Magnitude in i band           | mag           |
| lsst_rmag       | Magnitude in r band           | mag           |
| lsst_umag       | Magnitude in u band           | mag           |
| lsst_ymag       | Magnitude in y band           | mag           |
| lsst_zmag       | Magnitude in z band           | mag           |
| lum             | Luminosity                    | L_solar       |
| mact            | Actual solar mass             | M_solar       |
| mtip            | Mass at giant branch tip      | M_solar       |
| parentid        | Parent particle #             | NA            |
| partid          | 0 if at parent coords, else 1 | NA            |
| px              | Position x                    | kpc           |
| py              | Position y                    | kpc           |
| pz              | Position z                    | kpc           |
| ra              | Right Ascension               | degree        |
| rad             | Radial distance               | kpc           |
| smass           | Initial stellar mass          | M_solar       |
| teff            | Effective Temperature         | log(T/Kelvin) |
| vx              | Velocity x                    | km/s          |
| vy              | Velocity y                    | km/s          |
| vz              | Velocity z                    | km/s          |

---

## Configuraton Files <a name="configuration-files"></a>

The configuration files below are used in `population_plots.ipynb` and are all configured for LSST. They are labeled inner since they contain data from a radius of 5 kpc to a radius of 300 kpc, and require that the disk be removed in code. 

Simulation m12b:
* m12b_inner_356: snapshot taken before pericenter
* m12b_inner_385: snapshot at pericenter
* m12b_inner_396: snapshot taken after pericenter
* m12b_inner_600: current-day snapshot 

Simulation m12i:
* m12i_inner_600: current-day snapshot 

---

## Future Work <a name="future-work"></a>

Below is a list of the work that remains to be done:
* Fit the power law to the plots in `population_plots.ipynb`
* In `population_plots.ipynb` add region quadrant density plots for populations (BHB, MSTO, Kgiant)
* Update the README to describe each notebook in depth
* Add additional errors to our ananke output such as distance and proper motion errors (currently we only have the photometric errors)
* Rotate the outputed data so that the m12b satellite lines up with the LMC in the Milky Way
* Manually remove satellite objects using Rockstar data
* Create all-sky plots in each photometric band 
* Clean up the plots using matplotlib configuration and colormaps
* Plot metalicity distribution of the stars from FIRE
* Add extinction values using automated script 
