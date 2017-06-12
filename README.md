# CLAMATO2017

Code for pre-processing data required for tomographic reconstruction of Ly-a forest data from the 2017 CLAMATO Survey

GEN_DACH_INPUT.PRO: From the observed spectra, fits continua, applies pixel masks, and extracts Ly-a forest pixels to create a binary file in a format appropriate for tomographic reconstruction using the Dachshund code (Stark et al 2015, MNRAS, 453, 311)

INPUT.CFG: Input ASCII file for Dachshund code

GEN_CLAMATO_VEC.PRO: Similar to GEN_DACH_INPUT.PRO, but splits the data files into separate files for positions, fluxes, sigmas etc

PLOT_CL2017_SLICE_LONG_SMOOTHED.PRO: Generate 'slice maps' of the tomographic reconstruction.

PLOT_CL2017_SLICE_LONG_SMOOTHED_SKYMAP.PRO: Generate 'slice maps' of the tomographic reconstruction, but with little [RA/Dec] maps showing the position of each slicex