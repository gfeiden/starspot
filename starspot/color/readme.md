# Bolometric Corrections for Stellar Models 


## Compiling

Fortran 95 code included in this sub-package is written to be compatible with 
the NumPy/SciPy package `f2py`, allowing the Fortran code to be wrapped natively
into Python scripts. To compile and wrap the Fortran modules, execute the
following command:

```bash 
f2py -c -m bolcor utils.f95 interpolate.f95 bolcor.f95 
```

With installations of NumPy/SciPy on Mac with either Homebrew or MacPort,
`f2py` may need to be called with the Python version number, such as

```bash
f2py-2.7 -c -m bolcor utils.f95 interpolate.f95 bolcor.f95 
```

This creates a library called `bolcor.so` that can be imported into Python
routines.

## Usage

Within a Python script, the Fortran routines can be called after importing
the `bolcor` library. Modules and global variables are then accessible by 
calling `bolcor.module_name` or `bolcor.global_var_name`, but often one wants 
to access subroutines in the separate modules. This is achieved by calling 
`bolcor.module_name.subroutine(args)`. Here is a quick example to derive 
magnitudes for a single set of input parameters:

```python
import bolcor as bc

# initialize log file
bc.utils.log_init('example.log')

# transform stellar parameters into UBVRIJHK magnitudes
filters = ['U', 'B', 'V', 'R', 'I', 'J', 'H', 'K']
bc.bolcorrection.bc_init(0.0, 0.0, 'marcs', filters)

magnitudes = bc.bolcorrection.bc_eval(3000.0, 5.0, -2.65, 8)

bc.bolcorrection.bc_clean()
bc.utils.log_close()

```

For the moment, you *must* initialize the log file. This will hopefully 
be updated soon so the code recognizes whether this step has been performed.
It is also advisable to execute the clean and close routines to ensure 
that memory is properly released and that the log file stream is properly 
closed and the log data saved.

## References

MARCS bolometric correction tables computed for the entire grid of MARCS 
models (Gustafsson et al. [2008, A&A, 486, 951](http://adsabs.harvard.edu/abs/2008A%26A...486..951G)) 
by Casagrande & VandenBerg ([2014, MNRAS, 444, 392](http://adsabs.harvard.edu/abs/2014MNRAS.444..392C)). 
Full routines and tables are available online at

ftp://cdsarc.u-strasbg.fr/pub/cats/J/MNRAS/444/392/

PHOENIX bolometric corrections and color transformations computed for custom
grid of PHOENIX AMES-COND atmosphere model fluxes (Hauschildt et al. [1999,
ApJ, 512, 377](http://adsabs.harvard.edu/abs/1999ApJ...512..377H); [1999, ApJ,
525, 871](http://adsabs.harvard.edu/abs/1999ApJ...525..871H)), as described 
in Dotter et al. ([2007, PhDT](http://adsabs.harvard.edu/abs/2007PhDT........17D); 
[2008, ApJS, 178, 89](http://adsabs.harvard.edu/abs/2008ApJS..178...89D)). See

http://stellar.dartmouth.edu/models/

