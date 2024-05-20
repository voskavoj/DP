This repostitory contains the source code created as a part of the diploma thesis **Signals of opportunity positioning using LEO satellites** by Vojtěch Voska at the CTU FEE.

## How to install
1. Install Python 3.10 (https://www.python.org/)
1. Install miniconda (https://docs.anaconda.com/free/miniconda/index.html)
1. Create and activate a virtual environment in miniconda
	1. conda create NAME
	1. conda activate NAME
1. From Conda, install (conda install \<package\>):
	1. GNURadio 3.10 - `gnuradio` (https://anaconda.org/conda-forge/gnuradio)
	1. SoapySDR - `soapysdr` (https://anaconda.org/conda-forge/soapysdr)
	1. SoapySDR PlutoSDR - `soapysdr-module-plutosdr` (https://anaconda.org/conda-forge/soapysdr-module-plutosdr), for PlutoSDR
	1. gr-iridium - `gnuradio-iridium` (https://anaconda.org/conda-forge/gnuradio-iridium)
1. From Git, download iridium-toolkit (https://github.com/muccc/iridium-toolkit) to /External
1. Run pip install -r requirements.txt to install the required packages


## How to run
1. Activate the virtual environment (see How to install)
1. power up radio equipment
2. capture frames into file: `iridium-extractor -D 4 CONFIG >> Data/DATA/output.bits` where CONFIG is a path to configuration file, e.g. "\src\config\plutosdr-soapy.conf", and DATA is the name of the data folder (the folder must exist!)
2. either
	1. in src/config/setup.py, set EXP_NAME to DATA
	1. run `python External/iridium-toolkit/iridium-parser.py -p Data/DATA/output.bits --harder >> Data/DATA/decoded.txt`
	1. run src/app/process_offline_data.py
	1. run src/app/find_position_offline.py
3. or run src/app/find_position_from_offline_data.py