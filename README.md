MPC-PolyPhen
==============================

Python module to get MPC values and PolyPhen label predictions for a set of patients and mutations. \
*Last update: Sep 03, 2019*

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

This project contains a module called [getMPC.py](https://github.com/mgmartinezl/Stalicla-MPC-PolyPhen/blob/master/src/getMPC.py), 
which contains a set of functions integrated and called by the [main.py](https://github.com/mgmartinezl/MPC-pph2/blob/master/main.py) 
script. 

#### Positional parameters

* **inputFile:** it is mandatory to specify the absolute path to the input file containing the 
patient mutations.
    - Example: ```$ python3 main.py ~/data/raw/original-mutations-file.txt```
* **pathwaysDirectory:**  it is mandatory to specify the absolute path to the directory that 
contains the pathway files with gene annotations.
    - Example: ```$ python3 main.py ~/data/raw/original-mutations-file.txt /data/raw/pathways/```
* **inputMPC:**  it is mandatory to specify the absolute path to the directory that 
contains the MPC official values file. It is also possible to set a path to a directory that
contains chunks of the original file when it is partitioned.
    - Example 1: ```$ python3 main.py ~/data/raw/original-mutations-file.txt /data/raw/pathways/ path-to-mpc-official-values-file```
    - Example 2: ```$ python3 main.py ~/data/raw/original-mutations-file.txt /data/raw/pathways/ path-to-mpc-official-values-directory```

**IMPORTANT:** if an MPC official values file is provided, it must be a comma separated file.

#### Optional parameters

* **-pathway:** optional argument to process only specific pathways from the
mutations file. If no setting is provided,
all available pathways will be extracted by default.
    - Example 1: it will run only for pathway R-HSA-69620 \
    ```$ python3 main.py ~/data/raw/original-mutations-file.txt /data/raw/pathways/ -pathway R-HSA-69620 ```
    - Example 2: it will run for all possible pathways \
    ```$ python3 main.py ~/data/raw/original-mutations-file.txt /data/raw/pathways/ ```    
    
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
    ```$ python3 main.py ~/data/raw/original-mutations-file.txt /data/raw/pathways/ -gene CTR9 ```
   
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; For more than one gene, you can specify one of the following: \
\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; * Subset of genes separated by comma (without spaces):  \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; ```-gene CTR9,NOCL2 ``` \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; * A txt tab delimited file with no headers and the desired genes to filter written in the first column: \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; ```-gene path-to-my-file-containing-genes.txt``` \

* **-patient:** optional argument to process only specific patient IDs. If no setting is provided,
all available patients will be extracted by default.
    - Example: it will run only for patient with ID 1 \
    ```$ python3 main.py ~/data/raw/original-mutations-file.txt /data/raw/pathways/ -patient Patient_1 ```
   
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; For more than one patient ID, you can specify one of the following: \
\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; * Subset of patients separated by comma (without spaces):  \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; ```-patient Patient_X,Patient_Y ``` \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; * A txt tab delimited file with no headers and the desired patient IDs to filter written in the first column: \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; ```-patient path-to-my-file-containing-patientIDs.txt``` \

* **-mutation:** optional argument to process only specific mutations. If no setting is provided,
all available consequences will be processed by default.
    - Example: it will run only for mutations of type 'missense_variant' \
    ```$ python3 main.py ~/data/raw/original-mutations-file.txt /data/raw/pathways/ -mutation missense_variant```
   
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; For more than one mutation, you can specify one of the following: \
\
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; * Subset of mutations separated by comma (without spaces):  \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; ```-mutation missense_variant,Intron ``` \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; * A txt tab delimited file with no headers and the desired mutations to filter written in the first column: \
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; ```-mutation path-to-my-file-containing-mutations.txt``` \

**Note:** help() is available for all the parameters via the command line. 


### Input files
A sample of the file required by the first argument (the input file containing mutations and 
patients) can be found in the folder **data/raw** of the repo. The name of the file is
 **original-mutations-file.txt**. Also, a file called **dummy 

Similarly, in the [**data/raw/pathways**] folder, a sample of a directory containing pathways 
can be found, which is a mandatory parameter for _pathwaysDirectory_ entry. 

Additional example text files to filter can be found in the folder **data/raw/filters**.

## Running tests

In the folder [reports](https://github.com/mgmartinezl/Stalicla-MPC-PolyPhen/tree/master/reports), 
one example of the output generated by the script can be found, as well as a log containing the 
parameters set to run the program. To see the logs please visit **~/MPC-PolyPhen/analysis**.

## How to run this script

The scripts getMPC.py, and main.py are written in Python 3, 
which uses up-to-date libraries for this version as well.
 
To run the main module in a linux environment, simply call the script and the
arguments it needs:
 
To run the main module in a linux environment, simply call the script and the
arguments it needs:

```python3 main.py inputFile pathwaysDirectory inputMPC -pathway -gene -patient -mutation ```

For any additional information, contact me: 

*Gabriela Martinez* <br>
*airamgabriela17@gmail.com*

--------

<p><small>Project based on the <a target="_blank" href="https://drivendata.github.io/cookiecutter-data-science/">cookiecutter data science project template</a>. #cookiecutterdatascience</small></p>
