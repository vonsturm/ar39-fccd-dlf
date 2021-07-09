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
	--emin <opt>       : minimum energy for chi2 test [45-100]
	--emax <opt>       : maximum energy for chi2 test [100-200]
	-c --channel <opt> : channel
	-r <opt>           : rebin
	-o <opt>           : output directory
	--test <opt>       : test statistics even number plain, odd number delta ts
	                   : (0|1 = Chi2Test, 2|3 = Chi2Test/NDF, 4|5 = Kolmogorov, 6/7 Chi2 by-hand)
	-i --interpolate   : use gerda-ar39-pdf to interpolate between discrete pdfs in fccd and dlf
	-v                 : more output
```

Note:
- all options are optional, if not given the default value is used
- commmand line options can be given in arbitrary order
- options given via a config.json file are overridden by all other command line options

## Dependencies

The code depends on the `Root CERN` libraries to be compiled using the
following options (standard `c++17`)

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

and on the [progressbar](https://github.com/gipert/progressbar) utility by [gipert](https://github.com/gipert) included as submodule.

Interaction with `.json` files is done via [nlohmann/json](https://github.com/nlohmann/json). Consult the documentation for details.

## Auxillary files

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
