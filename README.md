## pymgpipe

### Installation
Run `pip install git+https://github.com/korem-lab/pymgpipe.git`<br/><br/>
You should now just be able to `import pymgpipe`. If you get a **ModuleNotFoundError**, make sure you're using the same python environment you used with `pip`. I highly recommend using the `pyenv` package to manage python versions & environments.


### Additional Dependencies
Need at least one of the following solvers-

-  [cplex](<https://www-01.ibm.com/software/commerce/optimization/cplex-optimizer/>)
-  [gurobipy](<http://www.gurobi.com>)

If you're running python 3.8.*, you should be able to just use `pip install cplex` or `pip install gurobipy`. For later versions of python, this might not work. Also, this does not actually create a license, it just installs the python interface to interact with these solvers. Both gurobipy and cplex offer free academic licenses. 

### Creating a local copy
Clone repository into desired directory. Run `pip install .` from within directory, or `pip install path/to/directory`<br/>
If you're having issues with dependencies, try running `pip install -r requirements.txt`

### Testing installation
To ensure you have **pymgpipe** running correctly, follow the example shown at `examples/workflow.ipynb`. Change the `solver` param to correspond to your solver of choice and run through the steps.

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
