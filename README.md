# CHAMP Helpers

A set of helper modules and scripts for running, generating inputs and
analyzing the outputs of [CHAMP](https://github.com/QMC-Cornell/CHAMP) program. 
Currently, they are used and tested specifically for 2D graphene-like systems 
(using [CHAMP Fork](https://github.com/g1obal/CHAMP)). However, basic outputs
of CHAMP's variational and diffusion Monte Carlo calculations (e.g. total energy
or correlation time) can be analyzed for all systems. 

The scripts and modules included in this repository can also be useful examples
for other projects.

## Requirements
numpy <br />
pandas <br />
matplotlib <br />
networkx <br />
scikit-learn <br />

This package is expected to work with various versions of the dependencies. 
However, if you notice any inconsistencies, see *requirements.txt* for tested
versions.

## Installation
Installation is not required since all scripts work as long as *champio*
package is located in the same directory with the corresponding script (or
included in PYTHONPATH environment variable), and the environment meets the 
library/package requirements.

Download or clone this repository.
```
$ git clone https://github.com/g1obal/champ-helpers.git
```

## Usage
The scripts located in the root directory are usable, and there are useful
scripts in *tools* directory.

For .py files:
```
$ python3 ScriptName.py
```

For .sh files:
```
$ bash ScriptName.sh
```

## Contact
For questions or bug reports, please create an issue on this repository or 
contact me via [e-mail](mailto:gooztarhan@gmail.com).
