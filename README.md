# Simulation using Best-fit Algorithm
SimBA (**Sim**ulation busing **B**est-fit **A**lgorithms) 


# Get Started

Please download a precompiled version of SimBA from the [releases page](https://github.com/ComputationalGenomics/SimBA/releases/). We provide binaries for Linux and macOS.

## How to run?

Here we provide a short description on how to use it from the command line.

The following command simulates 40 tetraploid samples based on 30 founders. Input markers are read from file ```input.vcf```. Output population is written in file ```output.vcf```:

```sh
$ simba-hap --founders 30 --samples 40 --ploidy 4 --input-vcf input.vcf --output-vcf output.vcf
```

## How to build packages

Go to the build folder:

```sh
$ cd build
```

Then run cmake, specifying the path to IBM CPLEX. For example, on Linux:

```sh
$ cmake .. -DSTATIC_BUILDS=ON \ 
           -DILOG_ROOT_DIR=/opt/ibm/ILOG/CPLEX_Studio127/cplex \ 
           -DILOG_CPLEX_LIBRARY=/opt/ibm/ILOG/CPLEX_Studio127/cplex/lib/x86-64_linux/static_pic/libcplex.a \
           -DILOG_CONCERT_LIBRARY=/opt/ibm/ILOG/CPLEX_Studio127/concert/lib/x86-64_linux/static_pic/libconcert.a
```

Or on macOS: 

```sh
$ cmake .. -DSTATIC_BUILDS=ON \ 
           -DCMAKE_CXX_COMPILER=/opt/local/bin/g++-mp-4.9 \
           -DILOG_ROOT_DIR=/opt/ibm/ILOG/CPLEX_Studio126/cplex/ \
           -DILOG_CPLEX_LIBRARY=/opt/ibm/ILOG/CPLEX_Studio126/cplex/lib/x86-64_osx/static_pic/libcplex.a \
           -DILOG_CONCERT_LIBRARY=/opt/ibm/ILOG/CPLEX_Studio126/concert/lib/x86-64_osx/static_pic/libconcert.a
```

Finally build the package by invoking make as follows: 

```sh
$ make package
```

# Citation

Please cite the following article if you use SimBA in your research:

E. Siragusa, N. Haiminen, F. Utro, L. Parida. Linear time algorithms to construct populations fitting multiple constraint distributions at genomic scales. IEEE/ACM Transactions on Computational Biology and Bioinformatics. 2017

# Contributing

We welcome contributions but request that you follow the guidelines indicated [here](https://github.com/ComputationalGenomics/SimBA/blob/master/Contributing/Contributing.md).

# Apache License v. 2.0
Copyright 2016 IBM Corporation

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

# Contact

If you have suggestions, questions or comments regarding SimBA, please email us to: 

pa_ri_da (at) us.ibm.com  (remove the underscores)

# Open Source @ IBM

Find more open source projects on the [IBM Github Page](http://ibm.github.io/)
