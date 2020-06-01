# optofluidics
Optofluidics is a Python library for the analysis of spectroscopic data stored in the HDF5 file format.

The HDF5 must contain measurement groups containing a set of timelapses. Each timelapse contains a collection of timestamped spectral data.

Each spectra is expected to have a background spectra and reference spectra defined stored as an attribute.

HDFS File
- measurement
  - timelapse_0
    - spectra_0
    - spectra_1
    - ...
  - timelapse_1
    - spectra_0
    - spectra_1
    - ...

## Installation
This package an be installed by cloning this repo. It will be for installation via PyPi in due course.

```bash
pip install optofluidics
```

## Usage
Optofluidics defines three object types: Datafiles, Datasets and Reactions.

 - A Datafile is the HDF5 file that may contain multiple Datasets.
 - A Dataset is one timelapse containing multiple spectra at different times.
 - A Reaction is a Dataset that has been processed to give concentration profiles.

```python
import optofluidics as of

datafile_1=of.Datafile('file_path') # loads the datafile
datafile_1.list_groups() # returns a list of measurements in the datafile
dataset_1=of.Dataset(datafile_1, group_key) # loads a specific dataset
dataset_1.pre_process() # returns a pd Dataframe with background-correction
dataset_1.calculate_abs() # returns a pd Dataframe with absorbances (calculated from reference spectra)
reaction_1=of.Reaction(dataset_1,wav_centre,epsilon,path_length) # initialises concentration profile
reaction_1.calculate_conc(wav_range) # calculates concentration using absorbance values for wav_centre +- wav_range
reaction_1.linear_drift(times_arr) # applies a linear drift correction by fitting to nil absorption points specified in times_arr
```

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.

## License
[MIT](https://choosealicense.com/licenses/mit/)
