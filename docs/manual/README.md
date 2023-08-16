# Manual: Hybdrid modelling with BMTK and COMSOL 

## Table of contents

<!---
TODO:
- github
- vsc
- nomachine
-->

- [Background](./background/)
  - [ESAT server](./background/esat.md)
  - [Terminal](./background/terminal.md)
  - [Python package management](./background/packages.md)
- [BMTK](./bmtk/)
  - [Installation guide](./bmtk/installation.md)
- [COMSOL](./comsol/)

## The Brain Modeling Toolkit

A software development package for building, simulating, and analyzing large-scale networks of different levels of resolution.

- [GitHub repo](https://github.com/AllenInstitute/bmtk)
- [Main paper](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1008386)
- [User guide](https://alleninstitute.github.io/bmtk/)


### [Installation guide](https://alleninstitute.github.io/bmtk/installation.html)

bmtk requires Python 2.7 or 3.5+, plus [additional python dependencies](#dependencies). There are three ways to install with base requirements from a command-line:

- Using your favourite python package manager
    ```bash
    $ pip install bmtk
    ```
     OR
    ```bash
    $ conda install -c kaeldai bmtk
    ```
  - Both pip and conda should automatcally download the necessary python [dependencies](#dependencies).
  - However, you will have to download any tutorials, examples, documentation... separately from the [GitHub repo](https://github.com/AllenInstitute/bmtk).
- Installing from the source
  ```bash
  $ git clone https://github.com/AllenInstitute/bmtk.git
  $ cd bmtk
  $ python setup.py install
  ```
  - This method will download create a copy of the github repo on your computer, giving you access to all tutorials, examples, documentation...
  - You will probably need to install python [dependencies](#dependencies) manually.
  - There are examples of building models and running simulations located in docs/examples/. Some of the simulation engines may require additional requirements to run.
- BMTK will work as is on your computer for simple examples, but if you want to take advantage of parallel computing on an HPC cluster, you also need the mpi4py python package.

#### Dependencies

- numpy
- h5py
- pandas
- matplotlib
- jsonschema
- pytest (optional for running unit tests)

BMTK will work as is on your computer for simple examples, but if you want to take advantage of parallel computing on an HPC cluster, you also need
- mpi4py

