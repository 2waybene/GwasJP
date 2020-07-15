GwasJP
===============================

version number: 0.0.1
author: Jianying Li

Overview
--------

This project is to help Dr. Motsinger-Reif for her accord clinical trial project pipeline automation.

Requirements
--------

GwasJP was developed using Python 3.7.1. In addition to requirements specified in setup.py, GwasJP requires installation of the following tools:

* [plink](http://zzz.bwh.harvard.edu/plink/)
* [smartpca](https://doi.org/10.1371/journal.pgen.0020190)
* [GCTA](http://cnsgenomics.com/software/gcta)

Most of the computational jobs need to be done on the high performance computer cluster. [SLURM](https://slurm.schedmd.com/documentation.html) cluster management is required and needs to be set up properly.

Prerequisit for the "accord" project
-------------------
The current development was designated for a historical clinical trial project ["the accord project"](https://github.com/2waybene/GwasJP/) -- need a real "accord" project link here

There are a few files need to be prepared outside this project

### Model step 1

- **/rundir/rhtn/ pheno_data_rhtn.txt**
- **/rundir/rhtn/rhtn_combined/forced_covars.txt**
- **/rundir/rhtn/rhtn_combined/starting_covars.txt**
- **/rundir/rhtn/rhtn_combined/phenotypes.txt**

Model step 2

In this package, a set of R codes are included and need to be downloaded and put in a separate bin folder. This "bin" folder needs to be set up properly

To successfully run this project, a set of genotype data is needed. They are developed under Dr. Motsinger-Reif's supervision and can be requested. These genotype data need to be placed in a "data" folder.

Example layout will be somthing like:

/homeDir/local/accord/bin/
<br>
/homeDir/local/accord/data/


Prerequisit for the EPR project
-------------------

> Here are some important dependency for the EPR project

### Third-party software

- **[plink](http://zzz.bwh.harvard.edu/plink/)**
- **[smartpca](https://doi.org/10.1371/journal.pgen.0020190)**
- **[GCTA](http://cnsgenomics.com/software/gcta)**
- **king**
- **etc.**


### Execution software (R codes etc.)
### Phenotype data preparation
### Genotype data
### Statistical modeling
### Main analysis work flow


Installation / Usage
--------------------

To install use pip:

    $ pip install GwasJP


Or clone the repo:

    $ git clone https://github.com/2waybene/GwasJP.git
    $ python setup.py install
    
Contributing
------------

TBD

Example
-------

TBD
