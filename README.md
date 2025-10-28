# Contextual-adaptation-provides-a-unified-explanation-for-multidimensional-prosociality

# DESCRIPTION
This project is designed for the production of the paper titled Contextual adaptation provides a unified explanation for multidimensional prosociality.

## Implementation Details
- Data generation codes were implemented in **Julia 1.9.4**
- Figure plotting codes were implemented in **Python 3.11.5**

# FILE STRUCTURE

## Core Analysis Scripts
- `frequencyPrediction.jl` - Computes structure-coefficients for given networks
- `networkModification.jl` - Generates group partitions for given networks
- `demo.jl` - An example

## Visualization Resources
- `Figure2/` - Data and code for Figure 2 
- `Figure3/` - Data and code for Figure 3
- `ExtendedDataFig1/` - Data and code for Extended Data Figure 1
- `ExtendedDataFig2/` - Data and code for Extended Data Figure 2
- `ExtendedDataFig3/` - Data and code for Extended Data Figure 3
- `ExtendedDataFig4/` - Data and code for Extended Data Figure 4
- `ExtendedDataFig5/` - Data and code for Extended Data Figure 5
- `ExtendedDataFig6/` - Data and code for Extended Data Figure 6
- `ExtendedDataFig7/` - Data and code for Extended Data Figure 7
- `ExtendedDataFig8/` - Data and code for Extended Data Figure 8
- `ExtendedDataFig9/` - Data and code for Extended Data Figure 9
- `ExtendedDataFig10/` - Data and code for Extended Data Figure 10
- `ExtendedDataFig13/` - Data and code for Extended Data Figure 13

## Empirical Data
- `EmpiricalNetwork_edgelist/` - Edge lists of four empirical social networks analyzed in the study

# INSTALLATION

## Julia Installation
Download and install Julia by following these [instructions](https://julialang.org/downloads/), or by running:

> Linux and MacOS:
>
> ```console
> $ curl -fsSL https://install.julialang.org | sh
> ```
>
> Windows:
>
> ```console
> > winget install julia -s msstore
> ```

## Run scripts

After installation, Julia will be accessible through your command-line interface. Execute scripts using:

```console
$ julia demo.jl
```

# Support
For any question about this program, please contact <br>
Yu Xia, Email: sjtuxy2019@sjtu.edu.cn
