## BPASS grid creation

This directory contains Python and Fortran scripts for creating new BPASS grids for arbitrary modelling choices (e.g. Initial Mass Function). Both the FOrtran and Python scripts should produce identical results.

## Set up

First download the spectra files from google drive (linked [here](https://bpass.auckland.ac.nz/9.html)) and `gunzip` them in the `bpass_files` folder.

## Python Scripts

You'll need to install [hoki](https://github.com/HeloiseS/hoki) for the file handling, and [spectres](https://github.com/ACCarnall/SpectRes) for resampling. Update `combine_bpass_spec.py` to point to the chosen files in the `glob` call, choose a suitable file name, then just run 

    python combine_bpass_spec.py
    python mass_file.py
 
This should create the binary grid file, as well as the mass table. Finally, you just need to update the file name (`spec_prep`, line 232), spectral resolution (`nspec`, line 234), and size of the metallicity grid (`nz`, line 224) in `sps_vars.f90`. Then re-install FSPS as normal (`make clean && make`).

## Fortran Scripts

`bpass_bin.f90` is a Fortran script for creating the same grids. You can compile it within the FSPS makefile if you add the relevant syntax, but need to include the BPASS compiler flag within `sps_vars.f90`.



