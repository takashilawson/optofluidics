# optofluidics
Optofluidics is a Python library for the analysis of spectroscopic data stored in the HDF5 file format.

The HDF5 must contain measurement groups containing a set of timelapses. Each timelapse contains a collection of timestamped spectral data.

Each spectrum is expected to have a background spectrum and reference spectrum defined stored as an attribute.

Modelling functions assume a photoreduction process with carbon dots as the photosensitiser.

Carbon dots are modelled with an error function to describe their activation period.

HDFS File
- measurement ABC
  - timelapse_0
    - spectrum_0
    - spectrum_1
    - ...
  - timelapse_1
    - spectrum_0
    - spectrum_1
    - ...
- measurement XYZ

## Installation
This package an be installed by cloning this repo. It will be available for installation via PyPi in due course.

```bash
pip install optofluidics
```

## Usage
Optofluidics defines three object types: Datafiles, Datasets and Reactions.

 - A Datafile is the HDF5 file that may contain multiple Datasets.
 - A Dataset is one timelapse containing multiple spectra at different times.
 - A Reaction is a Dataset that has been processed to give concentration profiles.

 Model creation expects certain params.

 params must consist of
     {
     t_erf representing the delay time in Erf model
     kc representing the complexation rate constant
     kbr representing the product bleaching rate constant
     k representing the saturation rate constant in the Erf model
     kr representing the back-reaction rate constant
     o representing the standard deviation in the Erf model
     c0 representing the initial concentration of reactant
     K representing the equilibrium constant for complexation
     end representing the total time for the model
     }

```python
from pre_proc import Datafile, Dataset
from reaction import Reaction
from standard_functions import *
from lmfit import Parameters

datafile_1=Datafile('file_path') # loads the datafile
datafile_1.list_groups() # returns a list of measurements in the datafile

dataset_1=Dataset(datafile_1, group_key, timelapse_key) # loads a specific dataset
dataset_1.pre_process() # returns a Pandas DataFrame with background-correction
dataset_1.calculate_abs() # returns a Pandas DataFrame with absorbances (calculated from reference spectra)

reaction_1=Reaction(dataset_1,wav_centre,epsilon,path_length) # initialises concentration profile
reaction_1.calculate_conc(wav_range) # calculates concentration using absorbance values for wav_centre +- wav_range/2
reaction_1.linear_drift(times_arr) # applies a linear drift correction by fitting to nil absorption points specified in times_arr
reaction_1.find_turn_points(first_point, prominence) # finds turning points (except the first one which you must specify)
reaction_1.create_model(params) # returns a Pandas DataFrame with model
reaction_1.fit_model(params,method) # runs a lmfit fitting routine to optimise the model parameters
reaction_1.model2csv(file_name, file_path) # saves model to csv file
```

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.

## License
[MIT](https://choosealicense.com/licenses/mit/)
