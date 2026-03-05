# Introduction 

Install Miniforge
#open miniforge prompt
conda create --name geoenv python=3.10  # Use a recent Python version
conda activate geoenv
conda install gdal
conda install -c conda-forge grass
# add additional required packages via requirements.txt
# install grass 7.8 and update variables in grass_functions.py