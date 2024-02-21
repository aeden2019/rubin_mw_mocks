# Make a packages directory
mkdir packages
cd packages/

# Install modules available through pip
pip install args
pip install pynbody

# Clone nba repo and install via pip
git clone https://github.com/jngaravitoc/nba
cd nba/
python -m pip install .
cd -

# Clone pygadgetreader repo and install via pip
git clone https://github.com/jngaravitoc/pygadgetreader
cd pygadgetreader/
python -m pip install .
cd -

# Clone py-ananke repo and install via pip
git clone https://github.com/athob/py-ananke
cd py-ananke
python -m pip install .
cd -

# Clone halo_analysis repo and install via pip
git clone https://bitbucket.org/awetzel/halo_analysis.git
cd halo_analysis
python -m pip install .
cd -
