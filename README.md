# HyperFlowCutter
Prototype implementation of the HyperFlowCutter algorithm and its refinement variant for hypergraph bipartitioning, accompanying our publication to appear at ESA19. arXiv preprint available soon.

# Setup
Get a recent C++17 ready compiler.
We use PaToH to obtain initial partitions as well as terminal pairs. Get PaToH at https://www.cc.gatech.edu/~umit/software.html and move the library (libpatoh.a) and header file (patoh.h) to the directory extern/
We use boost::dynamic_bitset for bitvectors. Get boost at https://www.boost.org/ or via your package manager.
We use tlx for parsing command line parameters. It is included as a submodule.

To compile, navigate to the top level directory and

```
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make HyperFlowCutter
make ReBaHFC
```
build ReBaHFC for the refinement algorithm, and build HyperFlowCutter for the plain HyperFlowCutter algorithm.

We read hypergraphs in hMETIS format. You can get the benchmark instances we used in the paper at https://algo2.iti.kit.edu/schlag/sea2017/benchmark_set.tar

To run the refinement algorithm with standard parameters on an initial partition obtained with PaToH's quality preset

```
./ReBaHFC -g <path to hypergraph> --eps <imbalance> --ip-eps <imbalance of initial partition> --patoh-preset Q
```
We recommend setting the same value for `--eps` and `--ip-eps`.
Use `--help` for an overview of the available options.

