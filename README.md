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

Additionally, key parameters can be modified in `src/mock/config.yaml`.

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

In this section you can find details regarding the input and ouput parameters for `py-ananke`. 

### Input Parameters <a name="ananke-parameters-input"></a> 

Below is a table of the parameters used with `py-ananke` to best simulate the characteristics of the Vera Rubin telescope. These values can be changes in `src/mock/make_mock.py` during the py-ananke process `ananke = an.Ananke()`.

| Parameter   | Value Used  | Description    |
| ----------- | ----------- | -------------- |
| particles   | p           | Dictionary of particle data from FIRE. |
| photo_sys   | 'padova/LSST' | Photometric system Galaxia should use to generate the survey. |
| cmd_magnames| ‘rmag,gmag-rmag’ | Names of the filters Galaxia should use for the color-magnitude diagram box selection. First argument is used for appMagLimits and the second for colorLimits. |
| app_mag_lim_lo | 17 | Lower limit in apparent mag. Related to cmd_magnames. |
| app_mag_lim_hi | 27.5 | Upper limit in apparent mag. Related to cmd_magnames. |
| abs_mag_lim_lo | -7.0 | Lower limit in absolute mag. |
| abs_mag_lim_hi | 10.0 | Upper limit in absolute mag. |

### Output Parameters <a name="ananke-parameters-output"></a> 

Below is a table of the parameters outputed by `py-ananke`. This is a work in progress, as indicated by the parameters with `Unkown` units. 

| Ouput Parameter | Description                   | Unit          |
| --------------- | ----------------------------- | ------------- |
| age             | Age                           | log (age/yr)  |
| alpha           | Alpha abundance               | [alpha/Fe]    |
| dec             | Declination                   | degree        |
| dmod            | Distance modulus              | Unknown       |
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
| lum             | Luminosity                    | Unknown       |
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
| vx              | Velocity x                    | Unknown       |
| vy              | Velocity y                    | Unknown       |
| vz              | Velocity z                    | Unknown       |

---
