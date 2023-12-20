# README

Author: Andrew Eden

E-Mail: aeden2019@my.fit.edu

---

## Summary

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

## Requirements

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

### Gizmo Analysis

This module is part of the processing for FIRE. Written by Andrew Wetzel, it can be found here: [https://bitbucket.org/awetzel/gizmo_analysis/src/master/](https://bitbucket.org/awetzel/gizmo_analysis/src/master/). Note that it also requires the [\utilities](https://bitbucket.org/awetzel/utilities/src/master/) package.  

The easiest installation is as follows:

```Bash 
git clone https://bitbucket.org/awetzel/gizmo_analysis.git
bash ./gizmo_analysis/install_helper.sh
```

### Py Ananke 

This module is foundational to our code.  Written by Adrien Thob, it can be found here: [https://github.com/athob/py-ananke](https://github.com/athob/py-ananke). 

One of the possible installation methods is as follows:
```Bash
git clone https://github.com/athob/py-ananke
cd py-ananke
pip install .

```

### Nba & Pygadgetreader

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

## Ananke Parameters

Work in progress

---
