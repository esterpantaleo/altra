###altra: a Bayesian method for simultaneous transcript reconstruction and abundance estimation with RNA-Seq data in multiple samples

This repository contains **altra**, a software implementing a Bayesian method
for simultaneous transcript reconstruction and abundance estimation with 
RNA-Seq data in multiple samples

**altra** is free software, you can redistribute it and/or modify it under
the terms of the GNU General Public License.

The GNU General Public License does not permit this software to be
redistributed in proprietary programs.

This library is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE
.
####Dependences

Installation instructions below assume that the *bedtools* and *samtools* are in your path and that the *gsl* and *boost dynamic_bitset* are installed on your system. Also the software uses *R* and the R package [*FLLat*](http://cran.r-project.org/web/packages/FLLat/index.html).

####Installation

To install the package, clone the repository in *~/src/altra*, then:

    cd ~/src/altra
    make all
    make check

and put folder *~/src/altra/scripts* in your `PATH` by adding the following lines to your *~/.bashrc* file:

    export PATH=$PATH":$HOME/src/altra/scripts/"

After adding this line to your *.bashrc*, remember to either login again, or do

    source ~/.bashrc


####Usage

To print the usage run **altra** with option `-h`:

    altra -h

