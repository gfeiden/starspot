# Bolometric Corrections and Color Transformations

MARCS bolometric corrections and color transformation tables computed for the
entire grid by Casagrande & VandenBerg (2014, MNRAS, 444, 392). Full routines
and tables available at

    ftp://cdsarc.u-strasbg.fr/pub/cats/J/MNRAS/444/392/

PHOENIX bolometric corrections and color transformations computed for custom
grid of PHOENIX atmosphere model fluxes, as described in Dotter et al. (2007,
PhDT; 2008, ApJS, 178, 89). See

    http://stellar.dartmouth.edu/models/

MARCS
-----

This file provides a detailed explanation of the content of BCtables.tar.gz and 
BCcodes.tar.gz, along with some details concerning the execution of their 
various programs. A quick step-by-step set of instructions 
(UNIX/Linux-oriented) is contained in the file INSTRUCTIONS.txt.

The BCtables.tar.gz file contains the following folders, each of which
contains many sub-directories:
1) alpha_m04   2) alpha_p00   3) alpha_p04   4) alpha_std

The BCcodes.tar.gz file contains the following Fortran computer programs:
1) getinputbcs.for  ... code to retrieve BC data for several fixed E(B-V) values
2) getbctable.for   ... code to generate BC transformation tables for any E(B-V)
3) gettestbcs.for   ... code to retrieve original MARCS data (see below)
4) chkbctable.for   ... verification program (see below)
5) bcutil.for       ... subroutines to interpolate in BC transformation tables
6) colstar.for      ... program to illustrate the use of bcinterp.for
7) bcstars.for      ... code to generate BCs for stars of known Teff,logg,[Fe/H]

In addition, BCcodes.tar.gz also contains the following files:
8) selectbc.data    ... to select photometric systems and MARCS models to use
9) input.sample     ... an example of input file for bcstars.for
10)bcgo             ... an executable file to compile, link and run bcstars.for

Selectbc.data (used to select the filters of interest), getinputbcs.for, and 
getbctable.for are described in considerable detail in the Appendix A of 
Casagrande & VandenBerg (2014), so we include just a few additional details 
concerning them here. Files 3-6 are fully described in sections C to E in the 
following documentation: they provide the means to verify that the program 
"getbctable" is working correctly. A detailed example is given in section F.

A few points worth mentioning/emphasizing:
1) The various codes have been set up for execution on Linux (or UNIX)
   computers, which use a forward slash to separate sub-directories, as opposed
   to a backslash that is used for the same purpose on WINDOWS machines.  If
   this package is implemented on a WINDOWS computer, it is necessary to edit
   "getinputbcs.for" and "gettestbcs.for" to change the definition of "slash"
   in the relevant DATA statements: comments are contained in these programs
   to help identify which statement must be corrected.
2) The interpolation subroutines (e.g., "getbc_std.for") will halt their
   execution with a STOP code if the input [Fe/H] or Teff values are outside
   the ranges spanned by the interpolation tables.  If the input log g value
   is too high or too low, the output bolometric corrections are set to a
   value of 9.9999.  It is unlikely that this will occur, as interpolations or
   extrapolations are performed at log g values from -0.5 to +5.5 at the
   coolest temperatures.
3) Only the "bc_std.data" file includes BCs for log g = -0.5.  The other 
   bc*.data files contain BCs for log g values from 0.0 to 5.5, though
   extrapolations to log g = -0.5 are performed when the data are read.
4) The interpolation subroutines send a message to unit 6 (normally associated
   with the monitor) when the bc_*.data files are read to inform the user
   of the E(B-V) value, the input [Fe/H] value, and the filters for which
   BC transformations are provided.  That (and any other) "write(6,*" 
   statement may be commented out if the user prefers to have no such output.
5) The output tables contain bolometric corrections, with or without the
   effects of reddening taken into account, *NOT* colors.  To obtain colors,
   the usual procedure is to calculate M_bol = 4.75 - 2.5 log L/Lsun, then
   M_i = M_bol - BC_i and M_j = M_bol - BC_j (for the i, j filters) and finally 
   the color i - j = M_i - M_j = BC_j - BC_i (where i is a bluer filter than j).
   Note that we have assumed M_bol = 4.75 for the Sun throughout this paper;
   i.e., all of the BC tables have been generated assuming this value.  If
   anyone wishes to adopt a different value, all of the BC tables that are
   generated using our software should be adjusted accordingly.


SOME DETAILS CONCERNING THE EXECUTION OF THE VARIOUS PROGRAMS FOLLOW:

*** Steps A and B are normally all that are needed to obtain the desired
*** BC tables.

A. "getinputbcs.for"

When executing "getinputbcs", a brief summary is sent to the monitor.  For
instance, if requesting transformations for the ubvri90 B and V filters,
"getinputbcs" will generate the following summary:
  files inputbc_r00.data, inputbc_r12.data, etc., have been created for
           ubvri90 : jc_B
           ubvri90 : jc_V
("r00" indicates tables for E(B-V) = 0.0, "r12" is attached to the name of
tables for E(B-V) = 0.12, etc.)  The "inputbc_r00.data", etc. files are
rewritten each time that this program is executed - and they can be deleted
when they are no longer needed.  Normally, "getbctable" should be executed
immediately after running "getinputbcs" in order to produce tables in the
form used by the interpolation subroutines.  

B. "getbctable.for"

This code uses Akima splines to interpolate in the output files produced by
executing "getinputbcs" (i.e., "inputbc_r00.data", etc.) for the value of
E(B-V) that is specified by the user.  The code asks the user to enter the
value of this quantity, which can be any number between 0.0 and 0.72.  The
output file will have one of 4 possible names (see below) - "bc_std.data",
"bc_p04.data", "bc_p00.data", or "bc_m04.data" depending on the value of
the parameter "ialf" in "selectbc.data".  If a file with the same name is
already resident in the home directory, it will be rewritten; otherwise,
a new file with this name will be created.  If, e.g., ialf and E(B-V) were
set to values of 1 and 0.0, respectively, the successful execution of
"getbctable" would generate the one-line summary:
   bc_std.data has been created assuming E(B-V) = 0.000
NOTE that "selectbc.data" is read by several of the codes, so it should
not be edited between the execution of "getinputbcs" and "getbctable" (for
instance).

*** The following provides the means to verify that everything is working
*** correctly.  Once sufficient familiarity with this package has been gained,
*** steps C to E can be bypassed.

C. "gettestbcs.for"

This code retrieves from the relevant directories the BC values derived
directly from MARCS model atmospheres to provide data sets that can be used to
check the results of the interpolation code.  These data are organized in data
files for fixed values of log g (e.g., "sphm05.data" contains BC values for
log g = -0.5, whereas "sph15.data" and "ppl40.data" contain such data for
log g = 1.5 and 4.0, respectively).  At log g = 3.5 and less, the BC values
were derived from spherical model atmospheres (as indicated by "sph" in the
file name), while they are based on plane-parallel atmospheres at higher
gravities (hence "ppl" in the file name).  This code asks the user to choose
one of the following six E(B-V) values: 0.0, 0.08, 0.14, 0.28, 0.44, and 0.56
(an integer must be entered to select one of these choices).  (These particular
values were chosen arbitrarily (i.e., at random) simply to provide values of
E(B-V) that differ from those assumed in the tables that are read by
"getinputbcs.for".)  Once this selection is made, "getbctable" must be executed
for the same reddening value so that the interpolation tables and the
"sphm05.data", "sph00.data", ..., "ppl55.data" files all assume exactly the
same value of E(B-V) that was assumed by the previous code.  For the example
used above (i.e., ubvri90 B and V filters), the execution of "gettestbcs" will
generate the following message (which is sent to the monitor):
  files sphm05.data, sph00.data, sph05.data, etc., have been created
  for E(B-V) = 0.00 to test
           ubvri90 : jc_B
           ubvri90 : jc_V
The "sphm05.data", etc. files are rewritten each time that this program is
executed- and they can be deleted when they are no longer needed.

D. "getbctable.for"

Any execution of "gettestbcs" should be *immediately* followed by the running
of "getbctable" so that the output table is generated for exactly the same
E(B-V) value that was assumed by the former code.

E. "chkbctable.for"

This code interpolates in the output file generated by "getbctable" for the
specific values of log g, Teff, and [Fe/H] that are listed in "sphm05.data",
"sph00.data", ..., and compares the BC values in the latter with those obtained
by interpolation in the former.  Any discrepancies larger than 0.001 are
flagged.  If the correct files have been created and the codes execute
properly, no such discrepancies should be found.  (The largest differences
are found when comparing the original and interpolated transformations for
ultraviolet filters, but even for them, the maximum differences are always
< 0.001.)  Zero differences will be found if the selected value of E(B-V) is
0.0 since the tabulated/interpolated BC values will be identical with those
derived directly from the MARCS atmospheres.


F. an example

First compile, in turn, the 4 FORTRAN codes (specifically, "getinputbcs.for",
"gettestbcs.for", "getbctable.for", and "chkbctable.for") that constitute
this software package.  If the Intel FORTRAN compiler on a Linux or WINDOWS
computer is used, this is accomplished by issuing the commands
                 ifort getinputbcs.for 
                 ifort gettestbcs.for
                       etc.
(Whatever commands are appropriate for different compilers and operating
systems obviously need to be employed.)  On a Linux machine, the executable
module is normally named "a.out": this should be renamed after each compilation
to be, e.g., "getinputbcs.exe" so that all of the executables have unique names.
Alternatively, the desired executable module can be obtained in one step via
the command:     ifort getinputbcs.for -o getinputbcs.exe       (etc.)
Once executable modules have been generated for them, there will be no further
need to recompile these codes.


Suppose transformations to the "hst_vega" system and the filters "f390m",
"f606w", and "f814w" are needed on the assumption of the standard variation
of [alpha/Fe] with [Fe/H].

Step 1: edit "selectbc.data" so that the integers on the first 5 lines have
        the following values:

  1  = ialf (= [alpha/Fe] variation: select from choices listed below)
  3  = nfil (= number of filter bandpasses to be considered; maximum = 5)
  4 14  =  photometric system and filter  { 4 selects the hst_vega photometric
  4 20  =  photometric system and filter  { system and 14, 20, 23 select the 
  4 23  =  photometric system and filter  { f390m, f606w, and f814w filters.
  7 37  ... not used (line can be left as is)
  7 39  ... not used (line can be left as is)

(The sixth and seventh lines would be used if nfil = 4 or 5.  In this example,
these lines will not be read, but they should be left as is in case they are
needed in the future.)

Step 2: execute "getinputbcs" (i.e., issue the command "getinputbcs.exe").  
        (One may then choose to jump directly to step 6 to generate the BC
        tables for the selected filters, with the effects of reddening taken
        into account.  The next 3 steps, which enable the user to verify that
        everything is working correctly and to check the accuracy of the
        interpolated tables for a few representative E(B-V) values, will
        normally be bypassed when one has gained sufficient familiarity with
        this package.

Step 3: execute "gettestbcs".  When the value of E(B-V) is requested, enter 3
        - one of the six permitted choices (so that, in this case, MARCS BCs
        for E(B-V) = 0.14 will be retrieved from the relevant directory).
        The data in the sph*.data and ppl*.data files are compared with the
        interpolated BCs in step 5.  As they are used only for testing purposes,
        the sph*.data and ppl*.data files may be deleted after completing
        step 5.

Step 4: execute "getbctable".  When asked to do so, enter 0.14 so that the
        tables created by "getinputbcs" are interpolated to the same value of
        E(B-V) that was selected in Step 3.

Step 5: execute "chkbctable".  The result of this execution will be the 
        following table, which lists for each gravity, the number of MARCS
        BC values, the number of interpolations that were performed, and in the
        last 6 columns, the BC values for the selected filter and the maximum
        differences found between the MARCS and interpolated BC values.  The
        largest of these differences occurred at the Teff and [Fe/H] values
        listed in the 3rd and 4th columns.  The agreement is clearly excellent
        given that the differences between the "input" and interpolated BC
        values is at the +/- 0.0001 level or better.  (Only in the case of 
        the HST uv filters (e.g., F218W) do the differences rise to as much
        as 0.0007 mag - and only for some reddening choices.  Thus, the 
        interpolation codes *always* reproduce the input MARCS values to
        within +/- 0.001 mag ... and usually much better than this.)

        input/
 log g checked  Teff  [Fe/H]   f390m   delta   f606w   delta   f814w   delta
 -0.5   57/ 57  4000   0.50  -6.8054  0.0001 -1.0218  0.0001  0.3335  0.0000
  0.0  161/161  5000  -0.50  -3.3588  0.0001 -0.4103  0.0000  0.4819  0.0000
  0.5  237/237  3100  -2.00  -6.4001  0.0000 -2.3103  0.0001  0.0608  0.0000
  1.0  278/278  3100  -0.75  -5.9080  0.0000 -3.8406  0.0001 -0.3713  0.0000
  1.5  318/318  6500   1.00  -1.7774  0.0001 -0.1344  0.0001  0.3585  0.0000
  2.0  332/332  3100  -1.00  -5.1269  0.0000 -3.1481  0.0001 -0.0354  0.0000
  2.5  400/400  6500   0.75  -1.5884  0.0001 -0.1650  0.0000  0.3372  0.0000
  3.0  411/411  3100  -2.00  -5.4495  0.0000 -1.6838  0.0001  0.2652  0.0000
  3.5  408/408  3100  -1.00  -5.0410  0.0001 -2.4185  0.0001  0.1899  0.0000
  4.0  418/418  3100  -0.25  -5.2412  0.0001 -2.6951  0.0000  0.0824  0.0000
  4.5  412/412  6250   0.00  -1.5128  0.0001 -0.2700  0.0001  0.3082  0.0001
  5.0  415/415  3100  -0.25  -5.7941  0.0001 -2.5948  0.0000  0.1139  0.0000
  5.5  169/169  3100   0.75  -5.5821  0.0000 -2.9712  0.0001 -0.0598  0.0000

(Steps 3, 4, and 5 may be repeated for other choices of the reddenings to
ensure that there are no problems at any of the 6 E(B-V) values for which 
the original MARCS results have been stored.)

Step 6: execute "getbctable" for the particular E(B-V) value of interest: the
        maximum permitted value of E(B-V) is 0.72.  The selected value of
        E(B-V) is recorded in the first line of the output "bc_*.data"
        file.  It is this file that should be used to transform absolute
        bolometric magnitudes to (in this example) absolute f390m, f606w,
        and f814w magnitudes.

In practice, many users will likely adopt the BC transformations for the
"standard variation" of [alpha/Fe] with [Fe/H], which assumes [alpha/Fe] =
+0.4 for [Fe/H] from -4.0 to -1.0, inclusive; 0.0 for [Fe/H] from +0.0 to
+1.0, inclusive; and [alpha/Fe] = -0.4 [Fe/H] at intermediate metallicities.
Transformations for constant values of [alpha/Fe] may be obtained for smaller
ranges of [Fe/H], as noted in "selectbc.data".  There is substantial overlap
of the assumed ranges in [Fe/H], so that, e.g., transformations for [Fe/H] =
0.0 and [alpha/Fe] = 0.0 can be obtained by running either the ialf=1 or 3
options.  Note that only in the case of ialf=1 are transformations provided
for log g = -0.5; normally, the range in log g is from 0.0 to 5.0 (but with 
extrapolations to log g = -0.5).

In any case, the summary line which is sent to the monitor when "getbctable"
is executed will say whether "bc_std.data" (ialf=1), "bc_p04.data", (ialf=2),
"bc_p00.data" (ialf=3), or "bc_m04.data" (ialf=4) has been created.  Separate
interpolation codes are available for these four cases; namely,
"marcsbc_std.for", "marcsbc_p04.for", "marcsbc_p00.for", and "marcsbc_m04.for",
respectively.  These subroutines are contained in the "bcutil.for" library.
Note that the "bc_*.data" files which are read and interpolated in by these
subroutines are attached to unit 30 (see the OPEN and CLOSE statements in the
latter).  In order to obtain the BC values appropriate to, e.g., [alpha/Fe] =
0.24, one will first have to interpolate in the tables for [alpha/Fe] = -0.4,
0.0, and 0.4 to obtain the BCs for the input values of the temperature, gravity,
and [Fe/H], and then perform a final interpolation in these BCs to obtain the
desired result.

"colstar.for" (also part of the package) is a simple code that calls all 
four subroutines to provide an example of how the latter are used.  To execute
"colstar", it is necessary to cycle through the ialf = 1, 2, 3, and 4 cases
in "selectbc.data" to generate the bc_std.data, bc_p04.data, bc_p00.data, and
bc_m04.data files (i.e., the BC interpolation tables for the four different
variations of [alpha/Fe] vs [Fe/H]).  Note that "bcutil.for" must be compiled
with "colstar.for" to create an executable module - by, e.g., issuing the
command:     ifort colstar.for bcutil.for             (on a WINDOWS computer)
or           ifort colstar.for bcutil.for -o colstar.exe   (Linux, UNIX)
The results generated by "colstar.exe" should match those reported in the
beginning commented section of "colstar.for".
 
"bcstars.for" is a simple code that given an input sample of stars of known 
logg, [Fe/H] and Teff ("input.sample") determines BCs (up to 5 filters, as 
selected in selectbc.data). This program can be compiled, linked and executed 
running the shell script "bcgo".
