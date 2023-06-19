# Constava

## Description
Provide a brief description of your project here.

## Getting Started

### Quick start

**To setup the source code:**

- Build: `make build`
- Install locally: `make install`
- Uninstall locally: `make uninstall`
- Publish on PyPI: `make publish`
- Clean up build and dist files: `make clean`

**To run the template example:**

Either as a Python module:

```
python -m template -m "Hello scientific community"
```

Or once you installed the Python package locally (`make install`):

```
constava_template -m "Hello scientific community"
```

**To add more modules:**
Feel free to add new directories and add the command line inside `setup.py` as it was made for the template example.

### Prerequisites
- Python 3.x
- pip

### Installation
To install the project, follow these steps:

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

### Usage

#### Build
To build the project, run the following command:

```
make build
```

This command will create the source distribution and wheel distribution files using `python3 setup.py`.

#### Install
To install the project, use the following command:

```
make install
```

This command will install the latest version of the project from the `dist` directory using `pip`.

#### Uninstall
To uninstall the project, execute the following command:

```
make uninstall
```

This command will uninstall the project package using `pip`.

#### Publish
To publish the project to a package repository, run the following command:

```
make publish
```

This command will upload the distribution files to the package repository using `twine`.

#### Clean
To clean the project by removing the build artifacts, use the following command:

```
make clean
```

This command will remove the `dist`, `build`, and `constava.egg-info` directories.

## Contributing

If you would like to contribute to this project, please follow these guidelines:
- Fork the repository
- Create a new branch
- Make your changes
- Open a pull request

## License
Include information about the license used for the project.

## Acknowledgments
Mention any acknowledgments or references used in the project.

## Contact
Provide your contact information if users have questions or want to reach out for support.
