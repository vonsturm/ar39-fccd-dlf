
# Ar39stat

Toy Monte Carlo sampler of simulated <sup>39</sup>Ar GERDA spectra. 
Different `FCCD` and `DLF` models and energy ranges can be specified
in a `.json` config file to &chi;<sup>2</sup>-test against the toy experiments.
All &chi;<sup>2</sup>-distributions are saved to file.

## Download the code

Download including submodules

```
git clone --recurse-submodules git@github.com:vonsturm/ar39-fccd-dlf.git
```

## Compilation 

A `Makefile` is present that saves the executables in `./bin`. Just invoke
`make` and you are all set.

## Usage

The sampler takes input arguments via command line and a master-config-file in `.json` format.
Consult the example master-config-file `./settings/ar39conf.json` and use `./bin/ar39stat --help`
to display the help menu

```
Chi2 test statistic sampling with Ar39 GERDA simulation

USAGE : ./ar39-stat --json config.json <options>

OPTIONS :
    -h --help          : print this help text
    --json <opt>       : master config file [conf.json]
    -c --channel <opt> : channel [0-40]
    --fccd <opt>       : fccd value in um [450-3000:50] must be available
    --dlf <opt>        : dlf value as fraction [0.00-1.00:0.05] must be available
    -s --stat <opt>    : statistics to be sampled in each toy experiment
    --emin <opt>       : minimum energy for chi2 test [45-100]
    --emax <opt>       : maximum energy for chi2 test [100-200]
    -t --toys <opt>    : number of toy experiments
    -b <opt>           : rebin
    -v                 : more output
```

Note:
- all options are optional, if not given the default value is used
- commmand line options can be given in arbitrary order 
- options given via a config.json file are overridden by all other command line options

Default values:
```
--json         : none
-c (--channel) : 0
--fccd         : 1000
--dlf          : 0.5
-s (--stat)    : 10000
--emin         : 45
--emax         : 150
-t (--toys)    : 100
-b             : 1
-v             : false
```


## Pre-requisits and dependencies

The code depends on the [ROOT](https://root.cern/) libraries to be
compiled using the following options (standard `c++17`)

```
cmake -Dhttp="ON"     \
      -Dmathmore="ON" \
      -Dminuit2="ON"  \
      -Droofit="ON"   \
      -DCMAKE_CXX_STANDARD=17 \
      -Dimt="ON"      \
      -Dvdt="ON"      \
      path/to/source
```

The filter algorithm is based on the [gram savitzky golay
library](https://github.com/arntanguy/gram_savitzky_golay) which should
be included as submodule at some point and compiled automatically

```
git clone --recursive
git@github.com:arntanguy/gram_savitzky_golay.git mkdir build && cd build
cmake -DCMAKE_BUILD_TYPE=Release .. make sudo make install
```

You need to install `Eigen3` (usually found in the repository of your
distribution of choice) for GSG to work and you may need to relink the
Eigen3 include files if they are installed in an unexpected place like a
subfolder in /usr/include or similar.

The [progressbar](https://github.com/gipert/progressbar) utility
by [gipert](https://github.com/gipert) is included as submodule.

Interaction with `.json` files is done via
[nlohmann/json](https://github.com/nlohmann/json). Consult the
documentation for details.


### Auxillary files

Download from mpik

Data
```
cd ar39-fccd-dlf 
mkdir -p data
cd data
rsync -av mpik:/lfs/l3/gerda/pertoldi/gerda/gerda-fitter/data/*v07.00.root
```

MC files
```
cd ar39-fccd-dlf 
rsync -avrz --copy-links --exclude edep --exclude prod-settings mpik:/lfs/l3/gerda/sturm/gerda/gerda-pdfs/releases/ph2p-ar39
```
