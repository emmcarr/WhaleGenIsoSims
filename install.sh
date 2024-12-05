#!/usr/bin/env sh

# set up virtual environment
python -m venv env
source env/bin/activate

# install prerequisites 
pip install --upgrade pip
pip install -r requirements.txt

# get simuPOP installed
git clone https://github.com/BoPeng/simuPOP.git simuPOP
cd simuPOP
python setup.py install --verbose
cd ..
