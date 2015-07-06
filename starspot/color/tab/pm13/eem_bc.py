# Quick script to reformat origianl Pecaut and Mamajek table
# for use in bolometric corrections routines.
#
import numpy as np

data = np.genfromtxt('EEM_dwarf_UBVIJHK_colors_Teff.txt')

# ascending Teffs
data = np.flipud(data)

# filters
filters = ['U', 'B', 'V', 'Rc', 'Ic', 'J', 'H', 'Ks']

# UBV(RI)c JHKs magnitudes
corrections = np.empty((len(data), len(filters)))
#--
corrections[:, 2] = data[:, 2]                       # V
corrections[:, 1] = corrections[:, 2] - data[:, 5]   # B
corrections[:, 0] = corrections[:, 1] - data[:, 7]   # U
corrections[:, 3] = corrections[:, 2] + data[:, 8]   # Rc
corrections[:, 4] = corrections[:, 2] + data[:, 9]   # Ic
corrections[:, 7] = corrections[:, 2] + data[:,10]   # Ks
corrections[:, 6] = corrections[:, 7] - data[:,12]   # H
corrections[:, 5] = corrections[:, 6] - data[:,11]   # J

table = np.column_stack((data[:, 0], corrections))

np.savetxt('Pecaut_Mamajek_2013_BCs.txt', table, fmt="%16.8f")
