# CDSfold Scripts
This directory contains several scripts for running and testing CDSfold. The test
scripts in the ```tests``` directory also import scripts in this directory.

## multirun.py
```multirun.py``` is a script for running the ```CDSfold``` binary many times, 
potentially with different command line arguments. The arguments for this script
are:

```
usage: multirun.py [-h] [--repeat REPEAT] [--out OUT] [-l]
                   [cds_args] [faa_path]

Script for running CDSfold multiple times.

positional arguments:
  cds_args         quote enclosed CDSfold options
  faa_path         path to faa file relative to CDS_HOME

optional arguments:
  -h, --help       show this help message and exit
  --repeat REPEAT  number of times to run CDSfold
  --out OUT        name of output file - requires .json extension
  -l               use list specified in the script instead of cmd line
```

The arguments for CDSfold can either be supplied from the command line or in the
python script itself. If the ```-l``` flag is not used, then the first argument
is the flags for CDSfold, enclosed in quotes. The second argument is the path
relative to ```CDS_HOME``` where the ```.faa``` file is. There are optional arguments for
specifying how many times to run CDSfold with the given arguments, and where to
save the output. Below is an example of running the script using CDS args specified
in the command line:

```
python3 multirun.py "-w 20 -e UAC, CCA" example/mev.faa --repeat 10 --out test.json
```

If the ```-l``` argument is used, then the arguments ```cds_args``` and ```faa_path```
is ignored. Instead, the arguments for CDSfold are specified by the ```RUN_SET```
variable in ```multirun.py```. See the documentation at the top of the script
for a description of how the list specifying command line args is formatted.
Below is an example of running with the ```-l``` option. The CDSfold arguments
and path to the ```.faa``` file are specified by the python variable ```RUN_SET```.

```
python3 multirun.py -l
```

The ```--out``` option argument specifies a path relative to the current
directory that the results are saved in. Two ```.json``` files are saved. One
file with ```_full``` appended includes run time, amino acid sequence, fold pattern,
and other information grouped by the ```.faa``` file. The second json file
groups sequences by common amino acid sequences in schema format for eterna upload.
Note that the ```puzzleID``` field is assigned starting at zero and incremented
for subsequent sequences; this field will need to be manually modified before
uploading to eterna.

## util.py
This file contains several utility functions used by the test script and by
```multirun.py```. This file is not a standalone script that can be run from
the command line.
