This repostitory contains the source code created as a part of the diploma thesis **Signals of opportunity positioning using LEO satellites** by Vojtěch Voska at the CTU FEE.

**The version valid at thesis submission is available in the Releases (https://github.com/voskavoj/DP/releases/tag/submission) or under the tag *submission* (https://github.com/voskavoj/DP/tree/submission)**

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


## Code structure
- **External** houses external files, namely it should include the iridium-toolkit
- **Data** contains folders with captured data
- **Plots** contains some plotting utilities, but most are in src/app
- **src** contains the source code
	- **app** contains applications and program entry points. Only these scripts should be run.
	- **config** contains configuration files and setup, parameters, locations, ...
	- **navigation** contains navigation calculations and data processing tools, and is (mostly) constellation-independed
	- **radio** contains files specific to radio support and constellation support, such as satellite identification, timing, frame parsing, ...
	- **satellites** contains the Satellite class, TLE download and operations, SGP4 predictions and generally anything to do with the actual satellites, (mostly) constellation-independed
	- **utils** contains utilities such as plotting and data dumping, ...
- **Utils** contains utilities such as bat files and command line instructions

## Code flow
**Signal capture, decoding**: Done separately, by gr-iridium and iridium-toolkit respectivelly.

**Processing**: Done by the script app/process_offline_data.
1. load TLEs (satellites/download_tles) from the data folder (argument offline_dir)
1. Load demodulated frames from Data directory
1. Find start time by the use IBC frame - radio/iridium_start_time
1. Process frames into array of time, f, fb, sat_id - radio/iridium_offline_radio
1. Process frames into navigation data - navigation/data_processing
1. save data as pickle

**Calculation**: Done by app/find_position_offline
1. load TLEs (satellites/download_tles) from the data folder (argument offline_dir)
1. Load data frames from Data directory
1. run solve from navigation/curve_fitting_method


## Extending the code
To implement a new position estimation method, put it in src/navigation. It must implement the *solve* function and it should be called from find_position_offline. Existing radio and processing chain can stay mostly intact.

To implement a new method of signal capture, a new radio chain in src/radio needs to be implemented, mainly the **frame_parsing** part. If IRA frames will not be decoded, the **channel identification** needs to change as well.

To support a new constellation, a new radio chain in src/radio needs to be implemented, mainly the **channel identification** part. In *navigation/data_processing*, the proper identification function needs to be called, and in *config/setup* the constellation needs to be added. Otherwise, the code should be able to remain intact.