# GetBlunted
For bluntifying overlapped GFAs

## Method
Overlaps are resolved using POA, as in the following diagram:
![POA alignment](https://github.com/rlorigro/GetBlunted/blob/dev/images/overlap_poa_diagram.svg)

The adjacency of nodes (as described by the Link lines in the GFA) is preserved. In particular, cases in which adjacency components are not fully connected (case C), new paths are **not** introduced between unlinked nodes. This is achieved by computing a biclique cover for each adjacency component and creating branches in the terminus of any node which participates in multiple bicliques.
![Diploid examples](https://github.com/rlorigro/GetBlunted/blob/dev/images/example_bluntification_cases.svg)

## Usage

```./GetBlunted /path/to/input.gfa```

For a typical phased human assembly GFA (5.1GB) about 8GB of RAM are used, and run time is 1m 30s.

## Installation

GetBlunted is compatible with Ubuntu (tested on 18.04) and MacOS

External dependencies: `libomp`, `cmake`

1. Clone the repo
2. `mkdir build`
3. `cd build`
4. `cmake ..`
5. `make -j [n_threads]`

Currently there is no make install step, so the executable is found in`GetBlunted/build/`
