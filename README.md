# Simultaenous Nitrogen Aggregation and Cooling (SNAC)
A package for modelling simultaneous nitrogen aggregation and cooling in diamond. See also Wincott et al. 2026[^1].

## Requirements
This package was developed under Python v. 3.12.9. All required packages are specified in requirements.txt

## Insallation
1) Install Python v. 3.12 or higher
2) Download or clone the SNAC package and place into a directory of your choice
3) Create a virtual environment in your working directory, then activate the newly created virtual environment:
    
    ```python -m venv \path\to\new\virtual\environment```

    or, in VSCode, open the Command Palette (Ctrl+Shift+P), search for the Python: Create Environment command, and select it. Then select a Python interpreter.

    ```.\venv\Scripts\activate```

4) install dependencies

    ```pip install -r requirements.txt```

5) install the snac package

    ```pip install -e .```

## Quick Start
Refer to Ipython notebooks in **documentation**:

automatedSNAC.ipynb

eventfulSNAC.ipynb

manualSNAC.ipynb


## Overview of modules and and their purpose
### snac

This is the main package. It contains the following modules.

#### diamond
This module defines the class Diamond. It is used for storing relevant information about a diamond (core, rim and kimberlite ages as well as nitrogen aggregation data). Diamond objects can be created from json files (Diamond.from_json() method) and stored as json files (Diamond.to_json() method). A Diamond object instance is passed to the SNACmodel (see below).

Example:

    from snac.diamond import Diamond
    diamond = Diamond(
        age_core=3520,
        age_rim=1860,
        age_kimberlite=0,
        c_NT=625,
        c_agg=0.863,
        r_NT=801,
        r_agg=0.197,
    )

#### SNACmodel
This module defines the class AggregationModel, which is used for the forward modelling of simultaneous nitrogen aggregation and cooling. After instantiating an AggregationModel object, its .run() method (`AggregationModel.run()`) can be used to optimise the initial temperature and cooling rate so the predicted aggregation state matches the measured data. Upon initialisation, the parameters cooling_rate0 and T_start0 are set, which will be used as first guesses once the run() method is used. The `AggregationModel.plot_T_history()` and `AggregationModel.plot_aggregation_history()` methods can be used to produce diagrams from the output. The temperature and aggregation history can be saved as csv via `AggregationModel.save_history(filename)` (see automatedSNAC.ipynb in **documentation**).

`AggregationModel.run()` accesses the function `aggregation.aggregate_and_cool()`, which is the computational foundation of the SNAC model. Beyond its basic use for predicting nitrogen aggregation during continuous cooling, a number of temperature scenarios can be modelled. Examples are provided in eventfulSNAC.ipynb within **documentation**.

Basic example:
    
    from SNACmodel import AggregtionModel
    model = AggregationModel(
        diamond=diamond,
        cooling_rate0=0.01,
        T_start0=1200,
        rate_bounds=(0.001, 0.12),
        T_bounds=(1000, 1450),
        dt=1
    )
    # fit model to measured N aggregation state of the diamond
    model.run()

    # plot results
    model.plot_T_history()
    model.plot_aggregation_history()

    # save results
    model.save_history('SNACoutput.csv')

Example with temperature scenario:

    from SNACmodel import AggregtionModel
    model = AggregationModel(
        diamond=diamond,
        cooling_rate0=0.01,
        T_start0=1200,
        rate_bounds=(0.001, 0.12),
        T_bounds=(1000, 1450),
        dt=1,
        T_scenario='hot_spike',
        scenario_params=(50, 1000, 25)
    )
    # fit model to measured N aggregation state of the diamond
    model.run()

    # plot results
    model.plot_T_history()
    model.plot_aggregation_history()

    # save results
    model.save_history('SNACoutput.csv')

**Note**: using `AggregationModel.get_history()` or `AggregationModel.plot_aggregation_history()` without prior fitting (i.e. without using `AggregationModel.run()` first) will use the initial guesses forward model nitrogen aggregation and cooling. This can be useful to explore the effect of different parameters on the outcome. See also manualSNAC.ipynb in **documentation**.

#### aggregation
This module contains a number of functions, including the main function used to model simultaneous nitrogen aggregation:
- **aggregate()**: this function is used to calculate the concentration of nitrogen in A-centres after a given duration spent at a given temperature
- **Temp_N()**: this function can be used to calculate a conventional model temperature based on the measured nitrogen aggregation state and diamond age
- **aggregate_and_cool()**: the computational heart of the SNAC model. This function models simultaneous nitrogen aggregation and cooling, including temperature "scenarios".

#### cooling
This module provides function describing cooling over time

### documentation
This directory contains Ipython notebooks that show how the snac package can be used. Each notebook contains a worked example.

- **automatedSNAC**: Example of using the SNAC model for automated modelling

- **eventfulSNAC**: Examples of using temperature scenarios for nitrogen aggregation modelling

- **manualSNAC**: Example of using the SNAC model to manually explore simultaneous cooling and nitrogen aggregation

### autoSNAC.py
This script allows calling the SNAC model from another script or, from the command line, e.g.:
```
python autoSNAC.py --save_dir "C:/User/SNAC/"  --model_file "C:/GitHub/SNAC/documentation/example_model.json" 
```

Note: You may have to activate your virtual Python environment before executing the script, e.g.
```
.\.venv\Scripts\activate
```


[^1]: Matt's paper!
