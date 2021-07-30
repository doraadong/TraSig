# Introduction
TraSig (**Tra**jectory-based **Sig**nalling genes inference) identifies interacting cell types pairs and significant ligand-receptors based on the expression of genes as well as the *pseudo-time ordering* of cells. For any two groups of cells that are expected to overlap in time, TraSig takes the pseudo-time ordering for each group and the expression of genes along the trajectory as input and then outputs an interaction score and p-value for each possible ligand-receptor pair. It also outputs a summary score for cell type pairs by combining individual ligand-receptors' scores. 

![flowchart](./method_diagram.png)

## Table of Contents
- [Get started](#Get&nbsp;started)
- [Command-line tools](#Command-line)
- [Tutorials](#Tutorials)
- [Learn more](#Learn&nbsp;more)
- [Credits](#Credits)

# Get started 
## Prerequisites 
* Python >= 3.6
* Python side-packages:   
-- numpy >= 1.19.5  
-- pandas >= 0.23.4  
-- Bottleneck >= 1.3.2  
-- statsmodels >= 0.12.1 (required for post-analysis only)  
-- scipy >= 1.5.4 (required for post-analysis only)  
-- matplotlib >= 3.3.4 (required for post-analysis only)  
-- seaborn >= 0.11.0 (required for post-analysis only)

## Installation 

### Install within a virtual environment 

It is recommended to use a virtural environment/pacakges manager such as [Anaconda](https://www.anaconda.com/). After successfully installing Anaconda/Miniconda, create an environment by following: 

```shell
conda create -n myenv python=3.6
```

You can then install and run the package in the virtual environment. Activate the virtural environment by: 

```shell
conda activate myenv
```

Make sure you have **pip** installed in your environment. You may check by 

```shell
conda list
```

If not installed, then: 

```shell
conda install pip
```

Then install TraSig, together with all its dependencies by: 

```shell
pip install git+https://github.com/doraadong/TraSig.git
```

If you want to upgrade TraSig to the newest version, then first uninstall it by:

```shell
pip uninstall trasig
```
And then just run the pip install command again. 

### Not using virtural environment

If you prefer not to use a virtual envrionment, then you may install TraSig and its dependencies by (may need to use **sudo**): 

```shell
pip3 install git+https://github.com/doraadong/TraSig.git
```

You may find where the package is installed by:
 
```shell
pip show trasig
```

# Command-line 

## Run TraSig

Run TraSig by (arguments are taken for example): 

```shell
main.py -i input -o output -d oligodendrocyte-differentiation-clusters_marques -g None -b ti_slingshot -n 1000 -s smallerWindow
```
The usage of this command is listed as follows:  

```shell
usage: main.py [-h] -i INPUT -o OUTPUT -d PROJECT -g PREPROCESS -b MODELNAME
               [-t LISTTYPE] [-l NLAP] [-m METRIC] [-z NAN2ZERO] [-n NUMPERMS]
               [-p MULTIPROCESS] [-c NCORES] [-s STARTINGTREATMENT]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        string, folder to find inputs
  -o OUTPUT, --output OUTPUT
                        string, folder to put outputs
  -d PROJECT, --project PROJECT
                        string, project name
  -g PREPROCESS, --preprocess PREPROCESS
                        string, preprocessing steps applied to the data /
                        project, default None
  -b MODELNAME, --modelName MODELNAME
                        string, name of the trajectory model
  -t LISTTYPE, --listType LISTTYPE
                        string, optional, interaction list type, default
                        ligand_receptor
  -l NLAP, --nLap NLAP  integer, optional, sliding window size, default 20
  -m METRIC, --metric METRIC
                        string, optional, scoring metric, default dot
  -z NAN2ZERO, --nan2zero NAN2ZERO
                        boolean, optional, if treat nan as zero, default True
  -n NUMPERMS, --numPerms NUMPERMS
                        integer, optional, number of permutations, default
                        10000
  -p MULTIPROCESS, --multiProcess MULTIPROCESS
                        boolean, optional, if use multi-processing, default
                        True
  -c NCORES, --ncores NCORES
                        integer, optional, number of cores to use for multi-
                        processing, default 4
  -s STARTINGTREATMENT, --startingTreatment STARTINGTREATMENT
                        string, optional, way to treat values at the beginning
                        of an edge with sliding window size smaller than nLap,
                        None/parent/discard/smallerWindow, default
                        smallerWindow, need to provide an extra input
                        'path_info.pickle' for 'parent' option
```

## Prepare inputs for TraSig (from dynverse outputs)

Prepare inputs by (arguments are taken for example): 

```shell
python prepare_inputs.py -i ../trajectory/input -o ../example/input -d oligodendrocyte-differentiation-clusters_marques -t ../trajectory/output/output.h5 -g None -b ti_slingshot -e None
```
The usage of this command is listed as follows:  

```shell
usage: prepare_inputs.py [-h] -i INPUT -o OUTPUT -d PROJECT -t TRAJECTORYFILE
                         -g PREPROCESS -b MODELNAME [-e OTHERIDENTIFIER]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        string, folder to find inputs for trajectory inference
  -o OUTPUT, --output OUTPUT
                        string, folder to save inputs for TraSig
  -d PROJECT, --project PROJECT
                        string, project name
  -t TRAJECTORYFILE, --trajectoryFile TRAJECTORYFILE
                        string, trajectory output file from dynverse, default
                        ../trajectory/output/output.h5
  -g PREPROCESS, --preprocess PREPROCESS
                        string, preprocessing steps applied to the data /
                        project, default None
  -b MODELNAME, --modelName MODELNAME
                        string, name of the trajectory model
  -e OTHERIDENTIFIER, --otherIdentifier OTHERIDENTIFIER
                        string, optional, other identifier for the output,
                        default None

```

## Analyze outputs from TraSig

Analyze outputs by (arguments are taken for example): 

```shell
python analyze_outputs.py -i ../example/input -o ../example/output -d oligodendrocyte-differentiation-clusters_marques -g None -b ti_slingshot -e None -n 100000 -s smallerWindow
```
The usage of this command is listed as follows:  

```shell
usage: analyze_outputs.py [-h] -i INPUT -o OUTPUT -d PROJECT -g PREPROCESS -b
                          MODELNAME [-t LISTTYPE] [-e OTHERIDENTIFIER]
                          [-l NLAP] [-m METRIC] [-z NAN2ZERO] [-n NUMPERMS]
                          [-s STARTINGTREATMENT]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        string, folder to find TraSig's inputs
  -o OUTPUT, --output OUTPUT
                        string, folder to find TraSig's outputs
  -d PROJECT, --project PROJECT
                        string, project name
  -g PREPROCESS, --preprocess PREPROCESS
                        string, preprocessing steps applied to the data /
                        project, default None
  -b MODELNAME, --modelName MODELNAME
                        string, name of the trajectory model
  -t LISTTYPE, --listType LISTTYPE
                        string, optional, interaction list type, default
                        ligand_receptor
  -e OTHERIDENTIFIER, --otherIdentifier OTHERIDENTIFIER
                        string, optional, other identifier for the output,
                        default None
  -l NLAP, --nLap NLAP  integer, optional, sliding window size, default 20
  -m METRIC, --metric METRIC
                        string, optional, scoring metric, default dot
  -z NAN2ZERO, --nan2zero NAN2ZERO
                        boolean, optional, if treat nan as zero, default True
  -n NUMPERMS, --numPerms NUMPERMS
                        integer, optional, number of permutations, default
                        10000
  -s STARTINGTREATMENT, --startingTreatment STARTINGTREATMENT
                        string, optional, way to treat values at the beginning
                        of an edge with sliding window size smaller than nLap,
                        None/parent/discard/smallerWindow, default
                        smallerWindow, need to provide an extra input
                        'path_info.pickle' for 'parent' option

```

# Tutorials

Github rendering disables some functionalities of Jupyter notebooks. We recommend using [nbviewer](https://nbviewer.jupyter.org/) to view the following tutorials. 

## Run TraSig on example data and analyze outputs 
The example inputs and outputs can be found under the folder [example](example). You may follow the [tutorial](tutorials/Run_TraSig_on_example_data.ipynb) to run TraSig on the example data and analyze the outputs. You may also obtain the analysis outputs by running the aforementioned script [analyze_outputs](tutorails/analyze_outputs.py) using command-line. See the tutorial for more details. 

## Prepare inputs 
To run TraSig, we need to have 4 input files. Here is a [tutorial](tutorials/Prepare_input_from_dynverse_ti_methods.ipynb), showing how to prepare these files from the inference results of any trajectory inference method included in [dynverse](https://dynverse.org/). You can find the example expression data (input) and trajectory inference result (output) under the folder [trajectory](trajectory). You may also prepare the inputs for TraSig by running the aforementioned script [prepare_inputs](tutorails/prepare_inputs.py) using command-line. See the tutorial for more details. 

# Learn more
Check our preprint at [biorxiv](https://www.biorxiv.org/content/10.1101/2021.07.28.454054v1). 

# Credits
The software is an implementation of the method TraSig, jointly developed by [Dora Li](https://github.com/doraadong), [Jun Ding](https://github.com/phoenixding) and Ziv Bar-Joseph from [System Biology Group @ Carnegie Mellon University](http://sb.cs.cmu.edu/). We also acknowledge Jeremy J. Velazquez, Joshua Hislop and Mo R. Ebrahimkhani from University of Pittsburgh for the fruitful discussions on method development. 

# Contacts
* dongshul at andrew.cmu.edu 

# License 
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details

