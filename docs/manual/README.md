# Manual: Hybdrid modelling with BMTK and COMSOL 

## The Brain Modeling Toolkit

A software development package for building, simulating, and analyzing large-scale networks of different levels of resolution.

- [GitHub repo](https://github.com/AllenInstitute/bmtk)
- [Main paper](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1008386)
- [User guide](https://alleninstitute.github.io/bmtk/)


### Quickstart

bmtk requires Python 2.7 or 3.5+, plus [additional python dependicies](https://alleninstitute.github.io/bmtk/index.html#base-installation). To install with
base requirements from a command-line:

```bash
 $ git clone https://github.com/AllenInstitute/bmtk.git
 $ cd bmtk
 $ python setup.py install
```

There are examples of building models and running simulations located in docs/examples/. Some of the simulation engines may require additional requirements to run.

##### Tests

There are a collection of unit tests in `bmtk.tests` which can be run using pytest

```bash
  $ cd bmtk
  $ py.test
```

### Documentation


- [Building network models](https://alleninstitute.github.io/bmtk/builder.html)
- [Running biophysical simulations](https://alleninstitute.github.io/bmtk/bionet.html)
- [Running point-neuron simulations](https://alleninstitute.github.io/bmtk/pointnet.html)
- [Running population-level simulations](https://alleninstitute.github.io/bmtk/popnet.html)

Copyright 2017 Allen Institute
