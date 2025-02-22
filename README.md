[![DOI](https://zenodo.org/badge/756268318.svg)](https://zenodo.org/badge/latestdoi/756268318)
[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/constava)](https://pypi.org/project/constava/)
[![PyPI - Version](https://img.shields.io/pypi/v/constava)](https://pypi.org/project/constava/)
[![Anaconda - Version](https://anaconda.org/bioconda/constava/badges/version.svg)](https://anaconda.org/bioconda/constava)
[![License](https://anaconda.org/bioconda/constava/badges/license.svg)](https://www.gnu.org/licenses/gpl-3.0.en.html)

<!-- PROJECT LOGO -->
<br /> 
<div align="center">
  <a href="bio2byte.be/b2btools" target="_blank" ref="noreferrer noopener">
  <img src="https://pbs.twimg.com/profile_images/1247824923546079232/B9b_Yg7n_400x400.jpg" width="80px"/>
  </a>

  # Constava
</div>

## Table of content

* [Table of content](#table-of-content)
* [Description](#description)
* [Installation](#installation)
  * [Prerequisites](#prerequisites)
  * [Installation through PyPI](#installation-through-pypi)
  * [Installation through conda](#installation-through-conda)
  * [Installation from the source](#installation-from-the-source)
* [Running constava from a containter (Docker)](#running-constava-from-a-container-docker)
  * [Using constava as a command line tool inside a Docker container](#using-constava-as-a-command-line-tool-inside-a-docker-container)
  * [Using constava as a library inside a Docker container](#using-constava-as-a-library-inside-a-docker-container)
* [Usage](#usage)
  * [Execution from the command line](#execution-from-the-command-line)
    * [Extracting backbone dihedrals from a trajectory](#extracting-backbone-dihedrals-from-a-trajectory)
    * [Analyzing a conformational ensemble](#analyzing-a-conformational-ensemble)
    * [Generating custom conformational state models](#generating-custom-conformational-state-models)
  * [Execution as a python library](#execution-as-a-python-library)
    * [Extracting backbone dihedrals as a DataFrame](#extracting-backbone-dihedrals-as-a-dataframe)
    * [Setting parameters and analyzing a conformational ensemble](#setting-parameters-and-analyzing-a-conformational-ensemble)
    * [Generating and loading conformational state models](#generating-and-loading-conformational-state-models)
    * [Constava-class parameters vs. command line arguments](#constava-class-parameters-vs-command-line-arguments)
* [License](#license)
* [Citation](#citation)
* [Authors](#authors)
* [Acknowledgements](#acknowledgments)
* [Contact](#contact)

[<Go to top>](#constava)

## Description

Constava analyzes conformational ensembles calculating conformational state 
propensities and conformational state variability. The conformational state 
propensities indicate the likelihood of a residue residing in a given 
conformational state, while the conformational state variability is a measure 
of the residues ability to transiton between conformational states.

Each conformational state is a statistical model of based on the backbone 
dihedrals (phi, psi). The default models were derived from an analysis of NMR
ensembles and chemical shifts. To analyze a conformational ensemble, the phi- 
and psi-angles for each conformational state in the ensemble need to be 
provided. 

As input data Constava needs the backbone dihedral angles extracted from the 
conformational ensemble. These dihedrals can be obtained using GROMACS' 
`gmx chi` module (set `--input-format=xvg`) or using the  `constava dihedrals` 
submodule, which supports a wide range of MD and structure formats.

[<Go to top>](#constava)

## Installation

### Prerequisites
- Python 3.8 or higher
- pip

[<Go to top>](#constava)

### Installation through PyPI

We recommend this installation for most users.

```sh
# Create a virtual environment (optional but recommended):
python3 -m venv constava
source constava/bin/activate

# Install the python module:
pip install constava

# Run tests to ensure the successful installation (optional but recommended):
constava test
```

If the package requires to be uninstalled, run `pip uninstall constava`. 

[<Go to top>](#constava)

### Installation through conda

To install constava through conda please follow the instructions below (both [Conda-Forge](https://conda-forge.org) and [Bioconda](https://bioconda.github.io) channels are required to install Constava dependencies).

```sh
# Create a conda environment (optional but recommended):
conda create -n constava python=3.12
conde activate constava

# Install constava
conda install -c bioconda -c conda-forge constava

# Run tests to ensure the successful installation (optional but recommended):
constava test
```

If the package requires to be uninstalled, run `conda remove constava`. 

[<Go to top>](#constava)

### Installation from the source

To download and install the latest version of the software from the source code
follow the instructions below.

```sh
# Clone the repository:
git clone https://bitbucket.org/bio2byte/constava/
cd constava

# Create a virtual environment (optional but recommended):
python3 -m venv constava
source constava/bin/activate

# Build and install the package from the packages root directory:
# ... build package from source
make build
# ... install it locally
make install
# ... test the installation
make test
```

If the package requires to be uninstalled, run `make uninstall` in the terminal 
from the package's root directory. 

[<Go to top>](#constava)

### Troubleshooting

#### Libtiff issues

If you run `constava` and see an error related to the library `libtiff` such as `libtiff.5.dylib' (no such file)`, you can try to fix it by installing `libtiff`. For instance, using conda: 

```sh
conda install libtiff
```

## Running constava from a container (Docker)

### Using constava as a command line tool inside a Docker container

To use constava's Docker image generated by the Biocontainers project based on 
the Bioconda package, follow the instructions below. You can find the container 
tags on https://quay.io/repository/biocontainers/constava?tab=tags. In this 
example, the latest tag is `1.1.0--pyhdfd78af_0`:

```sh
# Pull the constava image from quay.io
docker pull quay.io/biocontainers/constava:1.1.0--pyhdfd78af_0

# Run a container with the constava image
docker run \
  -it quay.io/biocontainers/constava:1.1.0--pyhdfd78af_0 \
  constava <COMMAND-LINE-OPTIONS>

# Optionally, you can mount a local directory to the container for accessing your data
docker run \
  -it -v /path/to/your/data:/data quay.io/biocontainers/constava:1.1.0--pyhdfd78af_0 \
  constava <COMMAND-LINE-OPTIONS>
```

To stop and remove the container, use the following commands:

```sh
# List all running containers
docker ps

# Stop a running container (replace <container_id> with the actual container ID)
docker stop <container_id>

# Remove the stopped container (replace <container_id> with the actual container ID)
docker rm <container_id>
```

If the image requires to be removed, run `docker rmi -f quay.io/biocontainers/constava:1.1.0--pyhdfd78af_0`.

### Using constava as a library inside a Docker container

To use constava as a library inside the Docker container, follow the instructions 
below. This allows you to interact with the constava library directly within a 
Python session inside the Docker container.

```sh
# Start an interactive Python session inside the constava container
docker run \
  --rm -it quay.io/biocontainers/constava:1.1.0--pyhdfd78af_0 \
  python
# This will start a Python shell where you can import constava
# >>> import constava
# >>>

# Alternatively, execute a python script inside the constava container
docker run \
  --rm -it quay.io/biocontainers/constava:1.1.0--pyhdfd78af_0 \
  python <python-script>
```

If the image requires to be removed, run `docker rmi -f quay.io/biocontainers/constava:1.1.0--pyhdfd78af_0`.

[<Go to top>](#constava)

## Usage

The software provides two modes of interaction. Shell user may use the software
from the command line, while users skilled in Python can import it as a module.
We provide a couple of usage examples in a [Colab notebook](https://colab.research.google.com/github/Bio2Byte/public_notebooks/blob/main/constava_examples.ipynb).

[<Go to top>](#constava)

### Execution from the command line

The software is subdivided in **three submodules**:

The `constava dihedrals` submodule provides a simple way to extract backbone 
dihedral angles from MD simulations or PDB ensembles. For more information
run: `constava dihedrals -h`. Alternatively, the backbone dihedrals may be
extracted with GROMACS' `gmx chi` module.

The `constava analyze` submodule analyzes the provided backbone dihedral angles
and infers the propensities for each residue to reside in a given 
conformational state. For more information run: `constava analyze -h`.

The `constava fit-model` can be used to train a custom probabilistic model of
confromational states.  The default models were derived from an analysis of NMR
ensembles and chemical shifts; they cover six conformational states:

* Core Helix - Exclusively alpha-helical, low backbone dynamics
* Surrounding Helix - Mostly alpha-helical, high backbone dynamics
* Core Sheet - Exclusively beta-sheet, low backbone dynamics
* Surrounding Sheet - Mostly extended conformation, high backbone dynamics
* Turn - Mostly turn, high backbone dynamics
* Other - Mostly coil, high backbone dynamics

[<Go to top>](#constava)

#### Extracting backbone dihedrals from a trajectory

To extract dihedral angles from a trajectory the `constava dihedrals` submodule 
is used.

```
usage: constava dihedrals [-h] [-s <file.pdb>] [-f <file.xtc> [<file.xtc> ...]] [-o OUTPUT] [--selection SELECTION] [--precision PRECISION] [--degrees] [-O]

The `constava dihedrals` submodule is used to extract the backbone dihedrals
needed for the analysis from confromational ensembles. By default the results
are written out in radians as this is the preferred format for
`constava analyze`.

Note: For the first and last residue in a protein only one backbone dihedral
can be extracted. Thus, those residues are omitted by default.

optional arguments:
  -h, --help            Show this help message and exit

Input & output options:
  -s <file.pdb>, --structure <file.pdb>
                        Structure file with atomic information: [pdb, gro, tpr]
  -f <file.xtc> [<file.xtc> ...], --trajectory <file.xtc> [<file.xtc> ...]
                        Trajectory file with coordinates: [pdb, gro, trr, xtc, crd, nc]
  -o OUTPUT, --output OUTPUT
                        CSV file to write dihedral information to. (default: dihedrals.csv)

Input & output options:
  --selection SELECTION
                        Selection for the dihedral calculation. (default: 'protein')
  --precision PRECISION
                        Defines the number of decimals written for the dihedrals. (default: 5)
  --degrees             If set results are written in degrees instead of radians.
  -O, --overwrite       If set any previously generated output will be overwritten.
```

An example:

```sh
# Obtain backbone dihedrals (overwriting any existing files)
constava dihedrals -O -s "2mkx.gro" -f "2mkx.xtc" -o "2mkx_dihedrals.csv"
```

[<Go to top>](#constava)

#### Analyzing a conformational ensemble

To analyze the backbone dihedral angles extracted from a confromational ensemble,
the `constava analyze` submodule is used.

```
usage: constava analyze [-h] [-i <file.csv> [<file.csv> ...]] [--input-format {auto,xvg,csv}] [-o <file.csv>] [--output-format {auto,csv,json,tsv}] [-m <file.pkl>] [--window <int> [<int> ...]]
                        [--window-series <int> [<int> ...]] [--bootstrap <int> [<int> ...]] [--bootstrap-series <int> [<int> ...]] [--bootstrap-samples <int>] [--degrees] [--precision <int>] [--seed <int>] [-v]

The `constava analyze` submodule analyzes the provided backbone dihedral angles
and infers the propensities for each residue to reside in a given 
conformational state. 

Each conformational state is a statistical model of based on the backbone 
dihedrals (phi, psi). The default models were derived from an analysis of NMR
ensembles and chemical shifts. To analyze a conformational ensemble, the phi- 
and psi-angles for each conformational state in the ensemble need to be 
provided. 

As input data the backbone dihedral angles extracted from the conformational 
ensemble need to be provided. Those can be generated using the 
`constava dihedrals` submodule (`--input-format csv`) or GROMACS'
`gmx chi` module (`--input-format xvg`).

optional arguments:
  -h, --help            Show this help message and exit

Input & output options:
  -i <file.csv> [<file.csv> ...], --input <file.csv> [<file.csv> ...]
                        Input file(s) that contain the dihedral angles.
  --input-format {auto,xvg,csv}
                        Format of the input file: {'auto', 'csv', 'xvg'}
  -o <file.csv>, --output <file.csv>
                        The file to write the results to.
  --output-format {auto,csv,json,tsv}
                        Format of output file: {'csv', 'json', 'tsv'}. (default: 'auto')

Conformational state model options:
  -m <file.pkl>, --load-model <file.pkl>
                        Load a conformational state model from the given pickled 
                        file. If not provided, the default model will be used.

Subsampling options:
  --window <int> [<int> ...]
                        Do inference using a moving reading-frame. Each reading 
                        frame consists of <int> consecutive samples. Multiple 
                        values can be provided.
  --window-series <int> [<int> ...]
                        Do inference using a moving reading-frame. Each reading 
                        frame consists of <int> consecutive samples. Return the 
                        results for every window rather than the average. This can
                        result in very large output files. Multiple values can be 
                        provided.
  --bootstrap <int> [<int> ...]
                        Do inference using <Int> samples obtained through 
                        bootstrapping. Multiple values can be provided.
  --bootstrap-series <int> [<int> ...]
                        Do inference using <Int> samples obtained through 
                        bootstrapping. Return the results for every subsample
                        rather than the average. This can result in very 
                        large output files. Multiple values can be provided.
  --bootstrap-samples <int>
                        When bootstrapping, sample <Int> times from the input data.
                        (default: 500)

Miscellaneous options:
  --degrees             Set this flag, if dihedrals in the input files are in 
                        degrees.
  --precision <int>     Sets the number of decimals in the output files.
  --seed <int>          Set random seed for bootstrap sampling
  -v, --verbose         Set verbosity level of screen output. Flag can be given 
                        multiple times (up to 2) to gradually increase output to 
                        debugging mode.
```

An example:

```sh
# Run constava with debug-level output
constava analyze \
    -i "2mkx_dihedrals.csv" \
    -o "2mkx_constava.json" --output-format json \
    --window 3 5 25 \
    -vv
```

[<Go to top>](#constava)

#### Generating custom conformational state models

To  train a custom probabilistic model of confromational states, the `constava fit-model` 
submodule is used. 

```
usage: constava fit-model [-h] [-i <file.json>] -o <file.pkl> [--model-type {kde,grid}] [--kde-bandwidth <float>] [--grid-points <int>] [--degrees] [-v]

The `constava fit-model` submodule is used to generate the probabilistic
conformational state models used in the analysis. By default, when running
`constava analyze` these models are generated on-the-fly. In selected cases 
generating a model beforehand and loading it can be useful, though.

We provide two model types. kde-Models are the default. They are fast to fit
but may be slow in the inference in large conformational ensembles (e.g., 
long-timescale MD simulations). The idea of grid-Models is, to replace
the continuous probability density function of the kde-Model by a fixed set
of grid-points. The PDF for any sample is then estimated by linear 
interpolation between the nearest grid points. This is slightly less
accurate than the kde-Model but speeds up inference significantly.

optional arguments:
  -h, --help            Show this help message and exit

Input and output options:
  -i <file.json>, --input <file.json>
                        The data to which the new conformational state models will
                        be fitted. It should be provided as a JSON file. The 
                        top-most key should indicate the names of the 
                        conformational states. On the level below, lists of phi-/
                        psi pairs for each stat should be provided. If not provided 
                        the default data from the publication will be used.
  -o <file.pkl>, --output <file.pkl>
                        Write the generated model to a pickled file, that can be
                        loaded gain using `constava analyze --load-model`

Conformational state model options:
  --model-type {kde,grid}
                        The probabilistic conformational state model used. The 
                        default is `kde`. The alternative `grid` runs significantly
                        faster while slightly sacrificing accuracy: {'kde', 'grid'}
                        (default: 'kde')
  --kde-bandwidth <float>
                        This flag controls the bandwidth of the Gaussian kernel 
                        density estimator. (default: 0.13)
  --grid-points <int>   This flag controls how many grid points are used to 
                        describe the probability density function. Only applies if
                        `--model-type` is set to `grid`. (default: 10000)

Miscellaneous options:
  --degrees             Set this flag, if dihedrals in `model-data` are in degrees 
                        instead of radians.
  -v, --verbose         Set verbosity level of screen output. Flag can be given 
                        multiple times (up to 2) to gradually increase output to 
                        debugging mode.
```

An example:

```sh
# Generates a faster 'grid-interpolation model' using the default dataset
constava fit-model -v \
    -o default_grid.pkl \
    --model-type grid \
    --kde-bandwidth 0.13 \
    --grid-points 6400
```

[<Go to top>](#constava)

### Execution as a python library

The module provides the `Constava` class a general interface to software's 
features. The only notable exception is the extraction of dihedrals,
which is done through a separate function.

[<Go to top>](#constava)

#### Extracting backbone dihedrals as a DataFrame

```python
import pandas as pd
from constava.utils.dihedrals import calculate_dihedrals

# Calculate dihedrals as a DataFrame
dihedrals = calculate_dihedrals(structure="./2mkx.pdb", trajectory="2mkx.xtc")

# Write dihedrals out as a csv
dihedrals.to_csv("2mkx_dihedrals.csv", index=False, float_format="%.4f")
```

[<Go to top>](#constava)

#### Setting parameters and analyzing a conformational ensemble

This example code will generate an output for a protein:

```python
# Initialize Constava Python interface with parameters
import glob
from constava import Constava

# Define input and output files
PDBID = "2mkx"
input_files = glob.glob(f"./{PDBID}/ramaPhiPsi*.xvg")
output_file = f"./{PDBID}_constava.csv"

# Initialize Constava Python interface with parameters
c = Constava(
    input_files = input_files,
    output_file = output_file,
    bootstrap = [3,5,10,25],
    input_degrees = True,
    verbose = 2)

# Alter parameters after initialization
c.set_param("window", [1,3,5])

# Run the calculation and write results
c.run()
```
This protein, with 48 residues and 100 frames per residue runs in about 1 minute.

The original MD ensembles from the manuscript can be found in 
[https://doi.org/10.5281/zenodo.8160755](https://doi.org/10.5281/zenodo.8160755).

[<Go to top>](#constava)

#### Generating and loading conformational state models

Conformational state models are usually fitted at runtime. This is usually the
safest option to retain compatibility. For `kde` models, refitting usually
takes less than a second and is almost neglectable. However, `grid` interpolation
models take longer to generate. Thus, it makes sense to store them when 
running multiple predictions on the same model.

**Note:** Conformational state model-pickles are intended for quickly rerunning 
simulations. They are **not for storing or sharing your conformational state models**. 
When you need to store or share a custom conformational state model, provide 
the training data and and model-fitting parameters.

```python
from constava import Constava

# Fit the grid-interpolation model
c = Constava(verbose = 1)
csmodel = c.fit_csmodel(model_type = "grid",
                        kde_bandwidth = .13,
                        grid_points = 10_201)

# Write the fitted model out as a pickle
csmodel.dump_pickle("grid_model.pkl")

# Use the new model to analyze a confromational ensemble
PDBID = "2mkx"
input_files = glob.glob(f"./{PDBID}_dihedrals.csv")
output_file = f"./{PDBID}_constava.csv"
c = Constava(
    input_files = input_files,
    output_file = output_file,
    model_load = "grid_model.pkl",
    input_degrees=True,
    window = [1, 5, 10, 25],
    verbose = 1)
c.run()
```

[<Go to top>](#constava)

#### Constava-class parameters vs. command line arguments

In the following table, all available parameters of the Python interface (`Constava` 
class) and their corresponding command line arguments are listed. The defaults for 
parameters in Python and command line are the same.

| Python parameter                      | Command line argument                                    | Description                                                                                                                                                                  |
|---------------------------------------|----------------------------------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `input_files : List[str] or str`      | `constava analyze --input <file> [<file> ...]`           | Input file(s) that contain the dihedral angles.                                                                                                                              |
| `input_format : str`                  | `constava analyze --input-format <enum>`                 | Format of the input file: `{'auto', 'csv', 'xvg'}`                                                                                                                           |
| `output_file : str`                   | `constava analyze --output <file>`                       | The file to write the output to.                                                                                                                                             |
| `output_format : str`                 | `constava analyze --output-format <enum>`                | Format of output file: `{'auto', 'csv', 'json', 'tsv'}`                                                                                                                      |
|                                       |                                                          |                                                                                                                                                                              |
| `model_type : str`                    | `constava fit-model --model-type <enum>`                 | The probabilistic conformational state model used. Default is `kde`. The alternative `grid` runs significantly faster while slightly sacrificing accuracy: `{'kde', 'grid'}` |
| `model_load : str`                    | `constava analyze --load-model <file>`                   | Load a conformational state model from the given pickled file.                                                                                                               |
| `model_data : str`                    | `constava fit-model --input <file>`                      | Fit conformational state models to data provided in the given file.                                                                                                          |
| `model_dump : str`                    | `constava fit-model --output <file>`                     | Write the generated model to a pickled file, that can be loaded again using `model_load`.                                                                                    |
|                                       |                                                          |                                                                                                                                                                              |
| `window : List[int] or int`           | `constava analyze --window <Int>  [<Int> ...]`           | Do inference using a moving reading-frame of <int> consecutive samples. Multiple values can be given as a list.                                                              |
| `window_series : List[int] or int`    | `constava analyze --window-series <Int>  [<Int> ...]`    | Do inference using a moving reading-frame of <int> consecutive samples. Return the results for every window rather than the average. Multiple values can be given as a list. |
| `bootstrap : List[int] or int`        | `constava analyze --bootstrap <Int> [<Int>  ...]`        | Do inference using <Int> samples obtained through bootstrapping. Multiple values can be given as a list.                                                                     |
| `bootstrap_series : List[int] or int` | `constava analyze --bootstrap-series <Int>  [<Int> ...]` | Do inference using <Int> samples obtained through bootstrapping. Return the results for every bootstrap rather than the average. Multiple values can be given as a list.     |
| `bootstrap_samples : int`             | `constava analyze --bootstrap-samples <Int> `            | When bootstrapping, sample <Int> times from the input data.                                                                                                                  |
|                                       |                                                          |                                                                                                                                                                              |
| `input_degrees : bool`                | `constava analyze --degrees`                             | Set `True` if input files are in degrees.                                                                                                                                    |
| `model_data_degrees : bool`           | `constava fit-model --degrees`                           | Set `True` if the data given under `model_data` to is given in degrees.                                                                                                      |
| `precision : int`                     | `constava analyze --precision <int> `                    | Sets the number of decimals in the output files. By default, 4 decimals.                                                                                                     |
| `kde_bandwidth : float`               | `constava fit-model --kde-bandwidth <float>`             | This controls the bandwidth of the Gaussian kernel density estimator.                                                                                                        |
| `grid_points : int`                   | `constava analyze --grid-points <int>`                   | When `model_type` equals 'grid', this controls how many grid points are used to describe the probability density function.                                                   |
| `seed : int`                          | `constava analyze --seed <int>`                          | Set the random seed especially for bootstrapping.                                                                                                                            |
| `verbose : int`                       | `constava <...> -v [-v] `                                | Set verbosity level of screen output.                                                                                                                                        |

[<Go to top>](#constava)

## License

Distributed under the GNU General Public License v3 (GPLv3) License.

[<Go to top>](#constava)

## Citation

> Gavalda-Garcia, J., Bickel, D., Roca-Martinez, J., Raimondi, D., Orlando, G., & Vranken, W. (2024). *Data-driven probabilistic definition of the low energy conformational states of protein residues*. NAR Genomics and Bioinformatics, 6(3), lqae082. [https://doi.org/10.1093/nargab/lqae082](https://doi.org/10.1093/nargab/lqae082)

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
[wim.vranken@vub.be](mailto:wim.vranken@vub.be)

<sup>&spades;</sup> Authors contributed equally to this work.

[<Go to top>](#constava)

## Acknowledgments

We thank Adrian Diaz [![ORCID](https://orcid.org/sites/default/files/images/orcid_16x16.png)](https://orcid.org/0000-0003-0165-1318) for the invaluable help in the distribution of this software. 

[<Go to top>](#constava)

## Contact

Wim Vranken - [wim.vranken@vub.be](mailto:wim.vranken@vub.be)

Bio2Byte website: [https://bio2byte.be/](https://bio2byte.be/)

[<Go to top>](#constava)
