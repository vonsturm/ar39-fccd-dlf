
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

## Compilation 

A `Makefile` is present that saves the executables in `./bin`. Just invoke
`make` and you are all set.

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
