## pymgpipe

### Installation
Run `pip install git+https://github.com/korem-lab/pymgpipe.git`<br/><br/>
You should now just be able to `import pymgpipe`. If you get a **ModuleNotFoundError**, make sure you're using the same python environment you used with `pip`. I highly recommend using the `pyenv` package to manage python versions & environments.


### Additional Dependencies
Need at least one of the following solvers-

-  [cplex](<https://www-01.ibm.com/software/commerce/optimization/cplex-optimizer/>)
-  [gurobipy](<http://www.gurobi.com>)

### Creating a local copy
Clone repository into desired directory. Run `pip install .` from within directory, or `pip install path/to/directory`<br/>
If you're having issues with dependencies, try running `pip install -r requirements.txt`
