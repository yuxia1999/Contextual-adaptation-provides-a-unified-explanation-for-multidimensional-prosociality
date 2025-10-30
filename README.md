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
- `frequencySimulation.jl` - Monte Carlo simulation of strategy frequencies in given networks.
- `demo.jl` - An example

**Note**: Input network nodes must be labeled consecutively starting from `0` when using this code for calculations.

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

All data are provided as ZIP archives. **Please unzip the corresponding data archive into the same directory before running the plotting code.**

## Empirical Data
- `EmpiricalNetwork_edgelist/` - Edge lists of four empirical social networks analyzed in the study

## Synthetic Data
- `SF_50_degree4_1.txt` – Edge list of a scale-free network with 50 nodes and an average degree of 4

# SYSTEM REQUIREMENTS
## Operating Systems

- **Tested on:** Windows 11 
- **Hardware:** Standard desktop computer (e.g., Intel i7, 16GB RAM)  

## Software Dependencies
- **Required Python Packages:**
```
json
matplotlib
numpy
```

- **Required Julia Packages**
```
ArgParse
DataFrames
Distributions
IterativeSolvers
LinearAlgebra
Printf
Random
SparseArrays
StatsBase
```

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

### Expected Output

The structure-coefficients lambda1, lambda2, lambda3, K1 and K2 of the network will be printed.
  
### Runtime for Demo

~20 seconds on a standard desktop computer.

# PSEUDOCODE FOR SIMULATION
```
## Inputs
- A network `G` (nodes represent individuals, edges represent potential interactions)
- For each individual `i`, a partition of its neighbors into groups:  
  `groups[i] = [G_i¹, G_i², ..., G_i^{m_i}]`
- Strategy set `S = {1, 2, ..., K}`
- Payoff matrix `P ∈ ℝ^{K×K}`, where `P[s, s′]` is the payoff to an individual using strategy `s` against a neighbor using `s′`
- Selection strength `δ ≥ 0`
- Mutation (exploration) rate `μ ∈ [0, 1]`
- Total number of update steps `T`

## Initialization
For each individual `i`:  
 For each group `g` in `groups[i]`:  
  `strategy[i][g] ←` uniformly random strategy from `S`

## Simulation Loop
For `t = 1` to `T`:
1. **Select an update target**:  
 Choose an individual `i` uniformly at random.  
 If `groups[i]` is non-empty, choose a group index `g` uniformly at random from `groups[i]`.

2. **Gather candidate strategies for imitation**:  
 Initialize an empty list `candidates`.  
 For each neighbor `j` of `i`:
  - Let `g_j` be the group in `groups[j]` that contains `i`.  
  - Let `s_j→i = strategy[j][g_j]` (the strategy `j` uses toward `i`).  
  - Compute `j`’s total payoff `u_j`:  
    `u_j = Σ_{k ∈ N(j)} P[ strategy[j][g_{jk}], strategy[k][g_{kj}] ]`  
    where `g_{jk}` is the group in `groups[j]` containing `k`,  
    and `g_{kj}` is the group in `groups[k]` containing `j`.  
  - Compute fitness: `f_j = exp(δ ⋅ u_j)`  
  - Append `(s_j→i, f_j)` to `candidates`.

3. **Update the selected strategy**:  
 If `random() < μ`:  
  `strategy[i][g] ←` random strategy from `S`  // exploration  
 Else if `candidates` is non-empty:  
  Select a strategy `s` from `candidates` with probability proportional to its associated `f_j`  
  `strategy[i][g] ← s`  // imitation  
 // Otherwise, retain the current strategy.

## Output
Compute the frequency of each strategy across all `(individual, group)` pairs over the simulation (or at the end), and report the distribution.
```

# SUPPORT
For any question about this program, please contact <br>
Yu Xia, Email: sjtuxy2019@sjtu.edu.cn
