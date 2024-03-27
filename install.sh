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
git fetch -v origin pull/115/head:pr-115
git checkout pr-115
python setup.py install --verbose
cd ..
