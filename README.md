## pymgpipe

### Installation
Run `pip install git+https://github.com/korem-lab/pymgpipe`<br/><br/>
You should now just be able to `import pymgpipe`. If you get a **ModuleNotFoundError**, make sure you're using the same python environment you used with `pip`. I highly recommend using the `pyenv` package to manage python versions & environments.


### Additional Dependencies
Need at least one of the following solvers-

-  [cplex](<https://www-01.ibm.com/software/commerce/optimization/cplex-optimizer/>)
-  [gurobipy](<http://www.gurobi.com>)

If you're running python 3.8.*, you should be able to just use `pip install cplex` or `pip install gurobipy`. For later versions of python, this might not work. Also, this does not actually create a license, it just installs the python interface to interact with these solvers. Both gurobipy and cplex offer free academic licenses. 

### Starting materials
To create multi-species community models with **pymgpipe**, you need two things to start-

-  Folder with individual taxa models (either in `.mat` or `.xml` format)
-  Relative abundance matrix (as a `.csv`) with samples as **columns** and taxa as **rows**. Taxa names should correspond to file names within taxa folder.

Examples of both can be found in the  `examples/` folder

### Testing solver
To ensure you have your solver setup correctly, you can run the following code snippet

```
from pymgpipe.test import test_pymgpipe
test_pymgpipe(solver='gurobi')
```
If everything is properly configured, this will run with no fatal exceptions.

### Testing model creation
In order to test model creation using provided examples, you'll have to first create a local copy of this repository (see below). After you've done that, you can follow the example shown at `examples/workflow.ipynb`. Change the `solver` param to correspond to your solver of choice and run through the shown steps.

After running all three steps, your directory should look like this

```
examples
│   workflow.ipynb
│   normCoverage.csv
│   normCoverage_formatted.csv
│   sample_label_conversion.csv
└───problems
│   │   ...
└───models
│   │   ...
└───fva
│   │   ...
└───panModels
│   │   ...
```

### Creating a local copy (optional)
If you want to mess around with this code, clone repository into desired directory. Run `pip install .` from within directory, or `pip install path/to/directory`. If you're having issues with dependencies, try running `pip install -r requirements.txt`
