MPC-PolyPhen
==============================

Python module to get MPC values and PolyPhen label predictions for a set of patients and mutations.
*Last update: Sep 02, 2019*

Project Organization
------------

    ├── README.md          <- The top-level README for developers using this project.
    │
    ├── analysis           <- The results and files created by the code, cytoscape sessions, etc.
    │
    ├── data
    │   ├── processed      <- The final, canonical data sets for modeling.
    │   └── raw            <- The original, immutable data dump.
    │
    ├── environment.yml   <- The requirements file for reproducing the analysis environment, e.g.
    │                         generated with `conda env export --no-builds > env.yml`
    │
    ├── notebooks          <- Jupyter notebooks. Naming convention is a number (for ordering),
    │                         the creator's initials, and a short `-` delimited description, e.g.
    │                         `1.0-jqp-initial-data-exploration`.
    │
    ├── references         <- Data dictionaries, manuals, and all other explanatory materials.
    │
    ├── reports            <- Generated analysis as HTML, PDF, LaTeX, etc.
    │   └── figures        <- Generated graphics and figures to be used in reporting
    │
    ├── src                <- Source code for use in this project.
    │   └── __init__.py    <- Makes src a Python module


### Description

This repo contains a Python module to annotate MPC values and pph2 prediction labels to a set 
of patients and mutations.

### Interface

This project contains a module called [getMPC.py](https://github.com/mgmartinezl/MPC-pph2/blob/master/getMPC.py), 
which contains a set of functions integrated and called by the [main-getMPC.py](https://github.com/mgmartinezl/MPC-pph2/blob/master/main-getMPC.py) 
script. 

#### Positional parameters

* **inputFile:** it is mandatory to specify the absolute path to the input file containing the 
patient mutations.
    - Example: ```$ python main-getMPC.py /home/Docs/mutations-file.txt```
* **pathwaysDirectory:**  it is mandatory to specify the absolute path to the directory that 
contains the pathway files with gene annotations.
    - Example: ```$ python main-getMPC.py /home/Docs/mutations-file.txt /home/pathways/```
* **inputMPC:**  it is mandatory to specify the absolute path to the directory that 
contains the MPC official values file. It is also possible to set a path to a directory that
contains chunks of the original file when it is partitioned.
    - Example 1: ```$ python main-getMPC.py /home/Docs/mutations-file.txt /home/pathways/ home/path-to-mpc-official-values-file```
    - Example 2: ```$ python main-getMPC.py /home/Docs/mutations-file.txt /home/pathways/ home/path-to-mpc-official-values-directory```

**IMPORTANT:** if an MPC official values file is provided, it must be a comma separated file.

#### Optional parameters

* **-pathway:** optional argument to process only specific pathways from the
mutations file. If no setting is provided,
all available pathways will be extracted by default.
    - Example 1: it will run only for pathway R-HSA-69620 \
    ```$ python main-getMPC.py /home/Docs/mutations-file.txt /home/pathways/ -pathway R-HSA-69620 ```
    - Example 2: it will run for all possible pathways \
    ```$ python main-getMPC.py /home/Docs/mutations-file.txt /home/pathways/ ```    
    
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; For more than one pathway, you can specify one of the following: \
\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; * Subset of pathways separated by comma (without spaces):  \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; ```-pathway R-HSA-69620,0051705 ``` \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; * A txt tab delimited file with no headers and the desired pathways to filter written in the first column: \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; ```-pathway path-to-my-file-containing-pathways.txt``` \

* **-gene:** optional argument to process only specific genes from the
mutations file. If no setting is provided,
all available genes will be extracted by default.
    - Example: it will run only for gene CTR9 \
    ```$ python main-getMPC.py /home/Docs/mutations-file.txt /home/pathways/ -gene CTR9 ```
   
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; For more than one gene, you can specify one of the following: \
\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; * Subset of genes separated by comma (without spaces):  \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; ```-gene CTR9,NOCL2 ``` \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; * A txt tab delimited file with no headers and the desired genes to filter written in the first column: \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; ```-gene path-to-my-file-containing-genes.txt``` \

* **-patient:** optional argument to process only specific patient IDs. If no setting is provided,
all available patients will be extracted by default.
    - Example: it will run only for patient with ID 1 \
    ```$ python main-getMPC.py /home/Docs/mutations-file.txt /home/pathways/ -patient Patient_1 ```
   
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; For more than one patient ID, you can specify one of the following: \
\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; * Subset of patients separated by comma (without spaces):  \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; ```-patient Patient_X,Patient_Y ``` \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; * A txt tab delimited file with no headers and the desired patient IDs to filter written in the first column: \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; ```-patient path-to-my-file-containing-patientIDs.txt``` \

* **-mutation:** optional argument to process only specific mutations. If no setting is provided,
all available consequences will be processed by default.
    - Example: it will run only for mutations of type 'missense_variant' \
    ```$ python main-getMPC.py /home/Docs/mutations-file.txt /home/pathways/ -mutation missense_variant```
   
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; For more than one mutation, you can specify one of the following: \
\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; * Subset of mutations separated by comma (without spaces):  \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; ```-mutation missense_variant,Intron ``` \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; * A txt tab delimited file with no headers and the desired mutations to filter written in the first column: \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; ```-mutation path-to-my-file-containing-mutations.txt``` \

**Note:** help() is available for all the parameters via the command line. 

### Input files
A sample of the file required by the first positional argument (the input file containing mutations and patients)
can be found in the folder **InputFiles**. For direct access, click [here.](https://github.com/mgmartinezl/MPC-pph2/blob/master/InputFiles/original-mutations-file.txt) 

Similarly, in the **InputFiles** folder, a sample of a directory containing pathways can be found, which is a mandatory
parameter for _pathwaysDirectory_ entry. See it directly [here.](https://github.com/mgmartinezl/MPC-pph2/tree/master/InputFiles/Pathways)

Additional example text files to filter 
[genes](https://github.com/mgmartinezl/MPC-pph2/tree/master/InputFiles/my-file-containing-genes.txt), 
[patients](https://github.com/mgmartinezl/MPC-pph2/tree/master/InputFiles/my-file-containing-patients.txt) and 
[mutations](https://github.com/mgmartinezl/MPC-pph2/tree/master/InputFiles/my-file-containing-mutations.txt) can be found 
in the folder as well.

## Running tests

In the folders [Logs]() and 
[Output](), three different running
examples can be found. Each of them generates a log containing the parameters set to run the
program, as well as the desired output.

## How to run this script

The scripts [getMPC.py](https://github.com/mgmartinezl/MPC-pph2/blob/master/getMPC.py), 
 and [main-getMPC.py](https://github.com/mgmartinezl/MPC-pph2/blob/master/main-getMPC.py) 
 are written in Python 3, which uses up-to-date libraries for this version as well.
 
To run the main module in a linux environment, simply call the script and the
arguments it needs:

```python3 main-getMPC.py inputFile pathwaysDirectory inputMPC -pathway -gene -patient -mutation ```

For any additional information, contact me: 

*Gabriela Martinez* <br>
*airamgabriela17@gmail.com*

--------

<p><small>Project based on the <a target="_blank" href="https://drivendata.github.io/cookiecutter-data-science/">cookiecutter data science project template</a>. #cookiecutterdatascience</small></p>
