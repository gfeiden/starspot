# Bolometric Corrections for Stellar Models 

Stellar evolution models yield predictions of stellar fundamental properties 
(radius, effective temperature, luminosity) as a function of time. While 
this is useful data, most observations do not yield this information directly.
Instead, observations measure the brightness of a star, or groups of stars,
in particular wavelength regimes, isolated using filters that permit only
certain wavelengths of light to pass through to the detector. It is therefore
advantageous to convert stellar fundamental properties to direct photometric
observables. This set of Fortran routines allows one to do just that.

## Compiling

To comiple directly for use in other Fortran programs, simply compile in
the following order,

```bash
gfortran -o program utils.f95 interpolation.f95 bolcor.f95 example.f95 
```

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

Bolometric corrections may be obtained for Johnson _UBV_, Cousins
_RI_, 2MASS _JHK_, and SDSS _ugriz_ bandpasses. HST magnitudes will be
available in the near future.

Here is a quick example of how to use these routines in a Fortran program
to compute bolometric correction for a single set of stellar parameters.

```fortran
program example
    use utils
    use bolcorrection
    implicit none
    
    real(dp) :: feh = 0.0_dp, afe = 0.0_dp
    real(dp) :: teff = 3000.0_dp, logg = 5.0_dp, logl = -2.65_dp
    real(dp), dimension(8) :: magnitudes
    
    character(len=15) :: brand
    character(len=1), dimension(8) :: filters = ['U', 'B', 'V', 'R', 'I', 'J', 'H', 'K']
    
    ! initialize log file
    call log_init('example.log')
    
    ! initialize bolometric correction table at fixed [Fe/H] and [a/Fe]
    call bc_init(feh, afe, brand, filters)
    
    ! transform stellar parameters to UBVRIJHK magnitudes
    call bc_eval(teff, logg, logl, size(filters, 1), magnitudes)
    
    ! release allocated memory and close log file
    call bc_clean()
    call log_close()

end program example
```

This example assumes that no changes are made to the Fortran routines. However,
if one is opting for Fortran over Python (see below), then it is possible
to simplify the call and allow magnitudes to be an allocatable array. Not 
crucial, it can just make life easier and the code neater. To do so, look 
into `bolcor` for segements of commented code referring to `f2py` not 
playing nicely.

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

# initialize bolometric correction table at fixed [Fe/H] and [a/Fe]
Fe_H = 0.0          # iron abundance relative to solar
A_Fe = 0.0          # alpha abundance relative to solar
brand = 'marcs'     # type of bolometric correction
filters = ['U', 'B', 'V', 'R', 'I', 'J', 'H', 'K']
bc.bolcorrection.bc_init(Fe_H, A_Fe, brand, filters)

# transform stellar parameters into UBVRIJHK magnitudes
Teff = 3000.0       # effective temperature, in K
logg = 5.0          # log10(gravity) in cgs units 
logL = -2.65        # log10(Luminosity / Luminosity_sun)
magnitudes = bc.bolcorrection.bc_eval(Teff, logg, logL, len(filters))

# release allocated memory and close log file
bc.bolcorrection.bc_clean()
bc.utils.log_close()

```

For the moment, you *must* initialize a log file. This will hopefully 
be updated soon so the code recognizes whether this step has been performed.
It is also advisable to execute the clean and close routines to ensure 
that memory is properly released and that the log file stream is properly 
closed and the log data saved.

Note that, for the moment, arguments for `bc_eval` should be scalars. In 
the future this will be updated to handle arrays for both spotted and un-spotted
stars.

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

Semi-empirical bolometric corrections were derived by VandenBerg & Clem 
([2003, AJ, 126, 778](http://adsabs.harvard.edu/abs/2003AJ....126..778V))
by empirically correcting theoretical bolometric corrections derived from 
MARCS and ATLAS9 model atmospheres. See their publication for details.

Empirical bolometric corrections for main-sequence stars at solar metallicity
are taken from Pecaut & Mamajek ([2013, ApJS, 208, 9](http://adsabs.harvard.edu/abs/2013ApJS..208....9P))
and maintained by Eric Mamajek (Version 2015.07.03). 
Data tables packaged with this software are modified compared to the original
tables. See [Eric Mamajek's webpage](http://www.pas.rochester.edu/~emamajek) 
for the original table "A Modern Stellar Color and Effective Temperature Sequence
for O9V - Y0V Dwarf Stars."

# Citations

If you make use of this code, please cite the relevant publications and 
note their contributions. 

* Feiden (in prep.) --- creation and distribution of Dartmouth model isochrones for young stars, including colors and magnitudes.

* Feiden & Christophe (in prep.) --- starspot modeling procedures.

* Casagrande & VandenBerg ([2014, MNRAS, 444, 392](http://adsabs.harvard.edu/abs/2014MNRAS.444..392C)) --- Bolometric correction tables using MARCS atmosphere fluxes.

* Dotter et al. ([2008, ApJS, 178, 89](http://adsabs.harvard.edu/abs/2008ApJS..178...89D)) --- Bolometric correction tables using PHOENIX AMES-COND atmosphere fluxes.

* Gustafsson et al. ([2008, A&A, 486, 951](http://adsabs.harvard.edu/abs/2008A%26A...486..951G)) --- Creation and distribution of the MARCS model atmosphere code and data products.

* Hauschildt et al. ([1999, ApJ, 512, 377](http://adsabs.harvard.edu/abs/1999ApJ...512..377H); [1999, ApJ,
525, 871](http://adsabs.harvard.edu/abs/1999ApJ...525..871H)) --- Creation and distribution of the PHOENIX AMES-COND model atmosphere code and data products.

* Pecaut & Mamajek ([2013, ApJS, 208, 9](http://adsabs.harvard.edu/abs/2013ApJS..208....9P)) --- Empirical, main-sequence, solar metallicity bolometric corrections.

Once the code is in a better and more user-friendly form, it will be provided
with an official release version number (via GitHub) and a doi number through
[Zenodo](http://www.zenodo.org). Stay tuned!
