# Optimal Design of Experiments for Chemical Experiments

This repository contains MATLAB scripts for creating optimal experimental designs for chemical experiments. The main functionality is provided by two scripts:

1. `calculateOptimalReducedDesign.m`
2. `createExperimentDesign.m`

## Functions

### calculateOptimalReducedDesign.m

This function constructs model-variance-optimized experimental designs using various optimality criteria:

- A-optimality
- D-optimality
- E-optimality
- G-optimality

#### Usage:
```
[ChosenDesign, CandidateAssessment] = calculateOptimalReducedDesign(levels, varargin)
```

- levels: A 1xN vector containing the factors and levels. Each number represents a different factor, with the magnitude of the number representing the number of levels for that particular factor.
- Optional parameters can be passed as name-value pairs, e.g., 'OptiCond', 'startValue', etc.

### createExperimentDesign.m
This function extends calculateOptimalReducedDesign with a user interface specifically for chemical experiments.
Usage:
```
FinalDesign = createExperimentDesign()
```

The function guides the user through an interactive process to input experiment parameters and then generates an optimal experimental design.

### Instructions

1. Start MATLAB and navigate to the directory containing the scripts.
2. Run createExperimentDesign:
3. Follow the prompts in the dialog boxes:

- Enter the maximum number of levels for any factor
- Enter the total number of factors
- Name the factors (e.g., Temperature, Time, Amount, Species)
- Fill in the table with specific values for each level
- Select which factors should be retained

4. The script will then calculate an optimal experimental design and return it as FinalDesign.

### Notes

- Calculations may take some time depending on the complexity of the experimental design.
- Ensure you provide valid inputs to avoid errors.
- The generated experimental designs are optimized to minimize variance and maximize efficiency.
