## [pymgpipe](https://korem-lab.github.io/pymgpipe/) | [![Python package](https://github.com/korem-lab/pymgpipe/actions/workflows/python-package.yml/badge.svg?branch=main)](https://github.com/korem-lab/pymgpipe/actions/workflows/python-package.yml) [![Tests](https://github.com/korem-lab/pymgpipe/actions/workflows/tests.yml/badge.svg?branch=main)](https://github.com/korem-lab/pymgpipe/actions/workflows/tests.yml) [![Docs](https://github.com/korem-lab/pymgpipe/actions/workflows/docs.yml/badge.svg)](https://github.com/korem-lab/pymgpipe/actions/workflows/docs.yml)  
<!-- Pytest Coverage Comment:Begin -->
<!-- Pytest Coverage Comment:End -->

### API Docss
https://korem-lab.github.io/pymgpipe/

### Installation
Run `pip install 'pymgpipe @ git+https://github.com/korem-lab/pymgpipe'`<br/><br/>
You should now just be able to `import pymgpipe`. If you get a **ModuleNotFoundError**, make sure you're using the same python environment you used with `pip`. I highly recommend using the `pyenv` package to manage python versions & environments.

### Additional Dependencies
Need at least one of the following solvers-

-  [cplex](<https://www-01.ibm.com/software/commerce/optimization/cplex-optimizer/>)
-  [gurobipy](<http://www.gurobi.com>)

If you're running python 3.8.*, you should be able to just use `pip install cplex` or `pip install gurobipy`. For later versions of python, this might not work. Also, this does not actually create a license, it just installs the python interface to interact with these solvers. Both gurobipy and cplex offer free academic licenses. 

### Inputs
To create multi-species community models with **pymgpipe**, you need two things to start-

-  Folder with individual taxa models (either in `.mat` or `.xml` format)
-  Relative abundance matrix (as a `.csv`) with samples as **columns** and taxa as **rows**. Taxa names should correspond to file names within taxa folder (excluding extension)

Examples of both can be found in the  `examples/` folder

### Outputs
The exact location and names of output files will vary depending on the parameters you pass into each function. However, the default output for pymgpipe's **build_models** function will look something like this-

```

* pymgpipe input*
.
├── taxaModels/
│   ├── taxa1.xml
│   ├── taxa2.xml
│   └── ...
├── abundances.csv

# pymgpipe output*
.
├── out/
    ├── reaction_content.csv
    ├── reaction_abundance.csv
    ├── sample_label_conversion.csv
    ├── metabolic_diversity.png
    ├── problems/
    │   ├── mc1.mps.gz
    │   ├── mc2.mps.gz
    │   └── ...
    └── models/
        ├── mc1.xml.gz
        ├── mc2.xml.gz
        └── ...
```

Here is a breakdown of each output files/directory and their descriptions-

| File | Type | Description |  
|---|---|---|
| reaction_content | CSV | Matrix showing binary presence/absence of all reactions within each sample | 
| reaction_abundance | CSV | Matrix showing scaled abundance (between 0 and 1) of all reactions within each sample  |  
| sample_label_conversion | CSV | Dictionary with conversion between original sample names and model names (default `sample_prefix` is 'mc') | 
| metabolic_diversity | PNG | Plot depicting # of unique reactions & taxa present within each sample | 
| problems | dir | Directory containing LP problems (default format is .mps, with `compressed` set to True) |  
| models | dir | Directory containing COBRA moddels (default format is .xml, with `compressed` set to True) | 

### Examples
Clone and run through `workflow.ipynb` in the **examples/** folder (see below)

### Creating a local copy (optional)
If you want to mess around with this code, clone repository into desired directory. Run `pip install .` from within directory, or `pip install path/to/directory` (if you've previously installed using github link above, you should run `pip uninstall pymgpipe` first to avoid conflicting versions). If you're having issues with dependencies, try running `pip install -r requirements.txt`

