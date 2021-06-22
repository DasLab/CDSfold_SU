# CDSfold SU (Stanford Update)
This repo is an update the of the CDSfold mRNA design algorithm to modern C++ while maintaining performance.

In addition, the repo fixes up some compilation issues with recent versions of ViennaRNA folding and explores use of higher temperature folding and 'jittered' energy functions to increase diversity of mRNA solutions.

The original CDSfold codebase is now available [on Github](https://github.com/gterai/CDSfold), but is not, to our understanding, under active development. As a result, a group of developers based at Stanford has chosen to implement this update as a separate, detached repository rather than a conventional fork.

## Install:
Before compiling the source codes, you need to install 
the Vienna RNA package (version 2.4.9), which is available  
from the [University of Vienna](http://www.tbi.univie.ac.at/RNA/)

Next, add the following lines to your ```~/.bashrc``` file:

```bash
# environment variables for CDSfold
export VIENNA_INSTALL="BUILD_DIR/ViennaRNA-2.4.9/src/ViennaRNA"
export CDS_HOME="CDS_PATH/CDSfold_SU"
```

where ```BUILD_DIR``` is the same path where ViennaRNA-2.4.9 is located and
```CDS_PATH``` is the path to where the ```CDSfold_SU``` directory is. Both
of these environment variables should be absolute paths.

Source your ~/.bashrc script:

```
source ~/.bashrc
```

Finally, run the following commands to build CDSfold:
``` 
cd CDSfold_SU/src
make
```

The executable CDSfold will showup in CDSFold_SU/src.  

## Usage:

```CDSfold [options] seqfile```

seqfile:

    amino acid sequence file in fasta format.

options:

```   
-w<integer>
   set the maximum distance between base-paired nucleotides.

-e<string>
   set codons to be avoided. When you set multiple codons, they 
   are specified by a comma-separated string (e.g. e=ACU,GCU).

-r
   perform a heuristic pseudo-MFE maximization.

-f<integer>
   design a CDS with partially stable or unstable secondary structure.
   The -f option indicates the start position of an interval in which the
   secondary structure should be stable or unstable. The option must
   be used together with -t option. See, examples below.

-t<integer>
   design a CDS with partially stable or unstable secondary structure.
   The -t option indicates the ending position of an interval in which the
   secondary structure should be stable or unstabl. See, examples below.

-R
   perform a random traceback

-C<float>
   temperature (in Celcius) the Vienna energy parameters are calculated at;
   default is 37 C.

-j<float>
   Enable jitter mode and specify the range of variance for Vienna energy
   parameters. This argument can only be passed when CDSfold is built with 
   the JitteredViennaEnergyModel using `make clean` followed by 
   `make jitter`. The argument must be in the range from 0 to 1. 
   In jitter mode, the energy parameters are multiplied by a random 
   number in the range [1-n, 1+n] where n is the value passed by the -j argument.
 ```

## Examples:
The following command finds a CDS with the most stable secondary structure
among all possible CDSs translated into the amino acid sequence in test.faa.

```./src/CDSfold ./example/test.faa```

To limit the maximum distance between base-pairs to be 20-bp, please try:

```./src/CDSfold -w 20 ./example/test.faa```

To avoid some codons (e.g. GUA and GUC), please try:

```./src/CDSfold -e GUA,GUC ./example/test.faa```

To design a CDS with unstable secondary structure by a heuristic pseudo-MFE
maximization, please try:

```./src/CDSfold -r ./example/test.faa```

To design a CDS of which the secondary strucrure from the first codon to 10-th
codon is unstable, please try:

```./src/CDSfold -r -f 1 -t 10 ./example/test.faa```

To design a CDS of which the secondary strucrure from the first codon to 10-th
codon is stable, please try:

```./src/CDSfold -f 1 -t 10 ./example/test.faa```

## Continuous Integration
Some basic CI is set up via GitHub, but we're seeing issues due to randomness in running test executables, so the CI is not very useful at the moment. 

Please see [issue #32](https://github.com/eternagame/CDSfold_SU/issues/32), and please file a pull request if you see the fix.


