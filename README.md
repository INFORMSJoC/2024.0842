[![INFORMS Journal on Computing Logo](https://INFORMSJoC.github.io/logos/INFORMS_Journal_on_Computing_Header.jpg)](https://pubsonline.informs.org/journal/ijoc)

# Exact MILP Models for Redundancy Allocation with Mixed Components

This archive is distributed in association with the [INFORMS Journal on
Computing](https://pubsonline.informs.org/journal/ijoc) under the [MIT License](LICENSE).

The software and data in this repository are a snapshot of the software and data that were used in the research reported on in the paper [Exact MILP Models for Redundancy Allocation with Mixed Components](https://doi.org/10.1287/ijoc.2024.0842) by Jianguang Feng, Zhi-Long Chen, Ada Che, Chengbin Chu. 


## Cite

To cite the contents of this repository, please cite both the paper and this repo, using their respective DOIs.

https://doi.org/10.1287/ijoc.2024.0842

https://doi.org/10.1287/ijoc.2024.0842.cd

Below is the BibTex for citing this snapshot of the repository.

```
@misc{feng2025exact,
  author =        {Jianguang Feng and Zhi-Long Chen and Ada Che and Chengbin Chu},
  publisher =     {INFORMS Journal on Computing},
  title =         {Exact {MILP} models for redundancy allocation with mixed components},
  year =          {2025},
  doi =           {10.1287/ijoc.2024.0842.cd},
  url =           {https://github.com/INFORMSJoC/2024.0842},
  note =          {Available for download at https://github.com/INFORMSJoC/2024.0842},
}  
```

## Description

This software provides the implementation of different exact solution methods introduced for solving redundancy allocation problems with mixed components described in the paper, including several different MILP models and a benchmark branch-and-bound method.


## Repository Structure

This repository is organized as follows:
- `src/`— Contains the source code
  - Files prefixed with `mm_*.h/.cpp` implement the MILP models.
  - Files prefixed with `bb_*.h/.cpp` implement the branch-and-bound method.
- `data`/— Contains benchmark instances for testing the proposed approaches. See the `README` file in this directory for data format details.
  - `data/H_2_4_Gamma_0_1.0`/ Data used in  used in Sections 4.2.1-4.2.3, the data with 5, 6, and 7 subsystems are from [1].
  - `data/H_5_6_Gamma_0_1.0`/ Data with 5 or 6 component types,  used in Sections 4.2.5. 
  - `data/H_2_4_Gamma_0_0.5`/ Data with component reliability drawn from [0, 0.5], also used in Sections 4.2.5. 
- `results`— Contains the results of computational experiments. See the `README` file in this directory for more details.
  - `results/H_2_4_Gamma_0_1.0`/ Results corresponding to Sections 4.2.1-4.2.3
  - `results/H_5_6_Gamma_0_1.0`/  Results for instances with 5 or 6 component types in Sections 4.2.5. 
  - `results/H_2_4_Gamma_0_0.5`/ Results with random reliability values in [0, 0.5], from Sections 4.2.5. 
- `cmake`/ — Includes custom CMake modules. 
- `CMakeList.txt` — Defines the build configuration for generating executables.

## Building

This project requires the following dependencies:
- [Gurobi](www.gurobi.com) (version 11.0 and later)
- [CMake](https://cmake.org/) (version 3.24 and later)
- [cxxopts](https://github.com/jarro2783/cxxopts) (version 3.1.1 and later)

To build the code for the MILP models: 

```shell
mkdir build
cmake -B build
cmake --build build --target ra_mm
```
To build the code for the branch-and-bound method:

```shell
mkdir build
cmake -B build
cmake --build build --target ra_bb
```

## Replicating

### Running `ra_mm` (MILP Models)

```shell
ra_mm [OPTION...]
```
Available options:

```
-h, --help              display usage information
-s, --system arg        system identifier associated with the instance: 1, 2, ..., 12
-f, --file arg          path to the instance data file
-m, --method arg        model to use: 
                        	110: [Q-BASIC], 120: [Q-TREE], 130: [Q-OB-TREE]
                        	210: [BASIC], 220: [TREE], 230: [OB-TREE]
-x, --time-limit arg    time limit in seconds, 0 for unlimited (default: 0)
-t, --test arg          test mode: 0, 1, 2, 3, 4 (default: 0)
```

> **Note:** Only the default value `0` for the `--test` option is relevant to the experiments reported in the paper.

Example:

```shell
ra_mm -s 1 -f rrap_ns5_nh2_m2_seed1.txt -m 210 -x 3600
```
This command runs the [BASIC] model (`210`) on system `1` with data from `rrap_ns5_nh2_m2_seed1.txt` and a time limit of 3600 seconds.

The output will be written to the **same directory** as the input data file. The output filename is constructed by adding a **prefix** and **suffix** to the original data filename, in the following format: 

```shell
s_<system identifier>_<original_filename>_<model name>_sol.txt
```

For the example above, the output file would be:

``````
s1_rrap_ns5_nh2_m2_seed1_BASIC_sol.txt
``````

------

### Running `ra_bb` (Branch-and-Bound)

```shell
ra_bb [OPTION...]
```
Available options:

```
-h, --help            display usage information
-s, --system arg      system identifier associated with the instance [1,2,...,12]
-f, --file arg        path to the instance data file
-x, --time-limit arg  time limit in seconds, 0 for unlimited (default: 0)
-t, --test arg        test mode: 0, 1, 2, 3, 4 (default: 0)
```

> **Note:** Only the default value `0` for the `--test` option is relevant to the experiments reported in the paper.

Example:

```shell
ra_bb -s 1 -f rrap_ns5_nh2_m2_seed1.txt -x 3600
```
This solves the instance for system `1` using the data in `rrap_ns5_nh2_m2_seed1.txt`, with a 3600-second time limit.

The output will be written to the **same directory** as the input data file. The output filename is constructed by adding a **prefix** and **suffix** to the original data filename, in the following format: 

```shell
s_<system identifier>_<original_filename>_BB_sol.txt
```

For the example above, the output file would be:

``````
s1_rrap_ns5_nh2_m2_seed1_BB_sol.txt
``````

## References

[1] [Young Woong Park](https://pubsonline.informs.org/action/doSearch?text1=Park%2C+Young+Woong&field1=Contrib) (2020) MILP Models for Complex System Reliability Redundancy Allocation with Mixed Components. INFORMS Journal on Computing 32(3):600-619. 
