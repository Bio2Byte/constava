<a id="readme-top"></a>

<!-- PROJECT LOGO -->
<br /> 
<div align="center">
  <a href="bio2byte.be/b2btools" target="_blank" ref="noreferrer noopener">
  <img src="https://pbs.twimg.com/profile_images/1247824923546079232/B9b_Yg7n_400x400.jpg" width="224px"/>
   </a>


# Constava
</div>

## Description
This software is used to calculate conformational states propensities & conformational state variability from a 
protein structure ensemble.
This is done by calculating the propensities for each conformational state for each residue in a protein ensemble. 
Then, conformational state variability is calculated from the change among conformational states, inferred from trained 
kernel density estimators (KDEs).

By default, this code retrains the conformational states propensities KDEs with the data set which we
provide, as described in the associated publication. This will generate KDEs that are compatible with your current
SciKit-learn version.

If the user wishes to train KDEs with a different set of dihedrals, a new set of dihedrals can be employed.
This set must be provided in a json file, with the name of the conformational states as keys and a list of lists [[phi, psi], [phi, psi], ...] as values.

<p align="right">(<a href="#readme-top">back to top</a>)</p>

## Getting Started

### Shell execution as a python package
#### Installation

1. Clone the repository:
   ```
   git clone https://bitbucket.org/bio2byte/constava/
   cd constava
   ```

2. Create a virtual environment (optional but recommended):
   ```
   python3 -m venv venv
   source venv/bin/activate
   ```

3. Install the project dependencies:
   ```
   pip install -r requirements.txt
   ```

Then, from the root directory run:

- Build: `make build`
- Install locally: `make install`

If the package requires to be uninstalled, run `make uninstall` in the terminal from the root directory. 

#### Shell execution usage

Once you installed the Python package locally (`make install`), the software's usage is as follows:

```
usage: constava [-h] -i INFILE -o OUTFILE [--input-format {auto,xvg,csv}] [--output-format {auto,csv,json}] [-k KDE] 
[-d TRAINING_DATA] [--kde-dump KDE_DUMP] [--use-publication-kdes USE_PUBLICATION_KDES]
[--window WINDOW [WINDOW ...]] [--bootstrap BOOTSTRAP [BOOTSTRAP ...]] [--bootstrap-samples BOOTSTRAP_SAMPLES] 
[--quick] [--degrees-to-radians]
```

with the flags described as: 

```shell
optional arguments:
  -h, --help            show this help message and exit

Input/Output Options:
  -i INPUT_FILE [INPUT_FILE ...], --input-file INPUT_FILE [INPUT_FILE ...]
                        Input file with dihedral angles (default: None)
  -o OUTPUT_FILE, --output-file OUTPUT_FILE
                        Output file (default: None)
  --input-format {auto,xvg,csv}
                        Format of input file (default: auto)
  --input-degrees       Add this flag if input is provided in degrees (instead of radians) (default: False)
  --output-format {auto,csv,json}
                        Format of output file (default: auto)

KDE options:
  -k <file.pkl>, --kdes <file.pkl>
                        Load KDEs from the given file (default: None)
  -d <data.json>, --kde-from-data <data.json>
                        For KDEs from the given data. The data is provided in a json file, with the name of the 
                        conformational states as keys and a list of lists [[phi, psi], [phi, psi], ...] as values (default: None)
  --kde-from-degrees    Add this flag if the data to fit the KDEs is provided in degrees (instead of radians) (default: False)
  --dump-kdes <file.pkl>
                        Dump the fitted KDEs as a pickled file, so that they can be reused later on using the --kdes flag. (default: None)
  --kde-bandwidth <float>

Miscellaneous Options:
  --window <int> [<int> ...]
                        Do inference using a moving reading-frame of <int> consecutive samples.
  --bootstrap <int> [<int> ...]
                        Do inference using <Int> samples obtained through bootstrapping. (By default a run with 3 and 25 is performed.)
  --bootstrap-samples <int>
                        If bootstrap, sample <Int> times from the input data (default: 500)
  --seed <int>          Set random seed for bootstrap sampling (default: None)
  --precision PRECISION <int> 
                        Sets de number of decimals in the output files (default: 4)
```



#### Shell execution example
An example command to run Constava would be: 

```shell
constava.py -i input_file.csv -o output_file.csv --bootstrap 5 10 15 --input-degrees
```

This example would run do the following things:
- Train the KDEs with the default data, since no additional training data was provided with the flag ```--training-data```. 
- Process ```input_file.csv``` and automatically detect the format, since no format was specified with the flag 
  ```--input-format```. It also converts the degrees of the phi-psi angles to radians because the flag 
  ```input-degrees``` is provided.
- Calculate conformational state variability for bootstrap (5, 500), (10, 500) and (15, 500), since the bootstrap size 
  was provided with the 
  flag ```--bootstrap``` and 500 is the default number of samples, which was not modified with the flag 
  ```--bootstrap-samples```.
- Output the results in ```output_file.csv```

<p align="right">(<a href="#readme-top">back to top</a>)</p>

#### Prerequisites
- Python 3.6 or higher
- pip

### Execution as a python library
#### Installation
We recommend installation via PyPI:

1. Create a virtual environment (optional but recommended):
   ```
   python3 -m venv venv
   source venv/bin/activate
   ```

2. Install the python module
    ```
   pip install constava
   ```
   
#### Python execution usage

This example code will generate an output for a protein:

```python
import constava
cons = constava.ConStaVa()
cons.train_kde()

# The argument infile accepts both shell syntax or a list of paths
# All the files in infile should be part of the same protein ensemble
cons.read_input_files(degrees=True, infile="../md_simulations/1akp/*.xvg")

cons.calculate_results(window=[1,2,4], bootstrap=[2,3,5,25])
cons.save_results(output_file='the_results.csv')
```

<p align="right">(<a href="#readme-top">back to top</a>)</p>


<!-- LICENSE -->
## License

Distributed under the GNU General Public License v3 (GPLv3) License.

<p align="right">(<a href="#readme-top">back to top</a>)</p>

<!-- AUTHORS -->
## Authors

- Jose Gavalda-Garcia<sup>&spades;</sup> 
[![ORCID](https://orcid.org/sites/default/files/images/orcid_16x16.png)](https://orcid.org/0000-0001-6431-3442) - 
[jose.gavalda.garcia@vub.be](mailto:jose.gavalda.garcia@vub.be)

- David Bickel<sup>&spades;</sup> 
[![ORCID](https://orcid.org/sites/default/files/images/orcid_16x16.png)](https://orcid.org/0000-0003-0332-8338) - 
[david.bickel@vub.be](mailto:david.bickel@vub.be)

- Joel Roca-Martinez 
[![ORCID](https://orcid.org/sites/default/files/images/orcid_16x16.png)](
https://orcid.org/0000-0002-4313-3845) - 
[joel.roca.martinez@vub.be](mailto:joel.roca.martinez@vub.be)

- Daniele Raimondi -
[![ORCID](https://orcid.org/sites/default/files/images/orcid_16x16.png)](https://orcid.org/0000-0003-1157-1899) - 
[daniele.raimondi@kuleuven.be](mailto:daniele.raimondi@kuleuven.be)

- Gabriele Orlando -
[![ORCID](https://orcid.org/sites/default/files/images/orcid_16x16.png)](https://orcid.org/0000-0002-5935-5258) - 
[gabriele.orlando@kuleuven.be](mailto:gabriele.orlando@kuleuven.be)

- Wim Vranken - 
[![ORCID](https://orcid.org/sites/default/files/images/orcid_16x16.png)](https://orcid.org/0000-0001-7470-4324) - 
[Personal page](https://researchportal.vub.be/en/persons/wim-vranken) - 
[wim.
  vranken@vub.be](mailto:wim.vranken@vub.be)

<sup>&spades;</sup> Authors contributed equally to this work.

<p align="right">(<a href="#readme-top">back to top</a>)</p>


<!-- ACKNOWLEDGMENTS -->
## Acknowledgments

We thank Adrian Diaz [![ORCID](https://orcid.org/sites/default/files/images/orcid_16x16.png)](https://orcid.org/0000-0003-0165-1318) for the invaluable help in the distribution of this software. 

<!-- CONTACT -->
## Contact

Wim Vranken - [wim.vranken@vub.be](mailto:wim.vranken@vub.be)

Bio2Byte website: [https://bio2byte.be/](https://bio2byte.be/)

<p align="right">(<a href="#readme-top">back to top</a>)</p>
