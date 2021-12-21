# FEM wrapper for Solid Mechanics simulations on Open-Source platforms
Yenan Shen, Anvitha Sudhakar, and Zhiyu Shi <br>
Dec. 20, 2021

## Introduction
In this project, we designed an wrapper for using open-source FEM pakages in solving continuum mechanics problem. There are six main folders in this repository: `mesh`, `material`, `residual`, `annealer`, `postprocess`, and `optimization`. `mesh` is resposible for generating mesh from the given yaml file. `material` assigns material to each discrete node and calculates the Lame constants. `residual` calculates the total energy and sets the gradient of the energy to be 0, solving the energy minimization problem. `annealer` does time discretization and gradually increases traction to avoid numerical error. `postprocess` outputs the result to .pvd (ParaView), .txt, or .png plot. `optimization` calculates the optimized mesh size for efficient computation.

We hope this package can allow newusers to use python open source FEM solves without learning their specific syntax. Forthe experienced users, we hope they don’t need to change the code every time when the boundary geometry, physical parameter, or elasticity model is changed.

## Installation
Firstly, using git clone command to download this pakage from remote repository to local computer. Once you have the repository, you can use Anaconda to create a new environment. The required packages are in `environment.yml` file. Use
```console
$ conda env create -f environment.yml
```
to create the conda environment with required packages. Then, use 
```console
$ conda acativate femtastic
```
to activate the environment. <br>
**Note that FEniCS pakage only works on Mac or Linux system. If you are working in windows, connect to a computer cluster to install the pakages**

## Run
To run the package, you first need to create a yaml file (\*.yml) that contains the parameters for this test. A sample yaml file can be found as `test.yml`. Users can change parameters in the yaml file for their own applications. In this yaml file, user need to specify mesh, material constant, equilibrium, annealer, and post-processing parameters. Once your own test.yml file is created, use
```console
$ python driver.py test.yml
```
to run the whole program. Depending on users' parameter, they will get .pvd, .txt, or .png files, along with the optimal mesh size calculation (depending on whether turning on the optimizer).

## Building Documentation
We used Sphinx for automated documentation. To build documents, if you want a pdf file , use
```console
$ make latexpdf
```
to get LaTex pdf version of the documentation. If you want an html file, use
```console
$ make html
```
to get html version of the documents.
