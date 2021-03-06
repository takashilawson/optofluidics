# optofluidics
Optofluidics is a Python library for the analysis of spectroscopic data stored in the HDF5 file format.

The HDF5 must contain measurement groups containing a set of timelapses. Each timelapse contains a collection of timestamped spectral data.

Each spectrum is expected to have a background spectrum stored as an attribute.

Absorbance calculations also assume a reference spectrum stored as an attribute.

Modelling functions assume a methyl viologen photoreduction process with carbon dots as the photosensitiser and EDTA as the sacrificial electron donor. Carbon dots are modelled with an error function to describe their delay/activation period. The bleaching of the photoreduced methyl viologen radical, and complexation of methyl viologen with EDTA, are both accounted for in the model.

HDF5 File Structure
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
This package can be installed using pip.

```bash
pip install --index-url https://test.pypi.org/simple/ optofluidics
```

## Usage
Optofluidics defines four object types: Datafiles, Datasets, Reactions and SV.

 - A Datafile is the HDF5 file that may contain multiple Datasets.
 - A Dataset is one timelapse containing multiple spectra at different times.
 - A Reaction is a Dataset that has been processed to give concentration profiles.
 - A SV is a Dataset that has been processed to give fluorescent counts.

```python
import optofluidics as of
from lmfit import Parameters

datafile_1=of.Datafile('file_path') # loads the datafile
datafile_1.list_groups() # returns a list of measurements in the datafile

dataset_1=of.Dataset(datafile_1, group_key, timelapse_key) # loads a specific dataset
dataset_1.pre_process() # returns a Pandas DataFrame with background-correction
dataset_1.calculate_abs() # returns a Pandas DataFrame with absorbances (calculated from reference spectra)

reaction_1=of.Reaction(dataset_1,wav_centre,epsilon,path_length) # initialises concentration profile
reaction_1.calculate_conc(wav_range) # calculates concentration using absorbance values for wav_centre +- wav_range/2
reaction_1.linear_drift(times_arr) # applies a linear drift correction by fitting to nil absorption points specified in times_arr
reaction_1.find_turn_points(first_point, prominence) # finds turning points (except the first one which you must specify)
reaction_1.create_model(params) # returns a Pandas DataFrame with model
reaction_1.fit_model(params,method) # runs a lmfit fitting routine to optimise the model parameters
reaction_1.model2csv(file_name, file_path) # saves model to csv file
```

Model creation expects certain params (lmfit parameter object).
- t_erf representing the delay time in the error function model
- kc representing the complexation rate constant
- kbr representing the product bleaching rate constant
- k representing the saturation rate constant in the Erf model
- kr representing the back-reaction rate constant
- o representing the standard deviation in the Erf model
- c0 representing the initial concentration of reactant
- K representing the equilibrium constant for complexation
- end representing the total time for the model

## Contributing
Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Please make sure to update tests as appropriate.

## License
[MIT](https://choosealicense.com/licenses/mit/)
