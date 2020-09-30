
# Ar39stat

Toy Monte Carlo sampler of simulated $^{39}$Ar GERDA spectra. 
Different `FCCD` and `DLF` models and energy ranges can be specified
in a `.json` config file to $\chi^2$-test against the toy experiments.
All $\chi^2$-distributions are saved to file.

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

and on the `progressbar` by `gipert` utility included as submodule to the code.

## Compilation 

A `Makefile` is present that saves the executables in `./bin`. Just invoke
`make` and you are all set.

## Auxillary files

### Data

Download from mpik

```
rsync -av mpik:/lfs/l3/gerda/pertoldi/gerda/gerda-fitter/data/*v07.00.root
```

### MC files

Download from mpik

```
rsync -avrz --copy-links --exclude edep --exclude prod-settings mpik:/lfs/l3/gerda/sturm/gerda/gerda-pdfs/releases/ph2p-ar39
```
