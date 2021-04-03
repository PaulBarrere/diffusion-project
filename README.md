M2R AAIS numerical method project: Scattering and Monte-Carlo method
===============================================================

![Mie diffusion](https://www.iap.uni-jena.de/iapmedia/de/Gruppe+Wyro/2018/Pic_Mie+Scattering-width-500-height-158.jpg)

*Credits: Institute of Aplied Physics (Universität Jena)*

Foreword
----------------

Light scattering by a single particule is well analitically described (Rayleigh and Mie theories). Nevertheless describing the scattering phenomenon through a layer composed of several particles, numerical methods are needed like the Monte-Carlo one using randomness. 

This project gathers FORTRAN 90 codes that simulate light scattering through a layer of spherical particles between 2 boundless parallel planes. First, isotropic and then Rayleigh scattering was implemented. Finally, light polarisation was considered in the frame of Mie formalism to simulate scattering of light by bigger particules.

Post-processing Python routines and the Makefile are also saved in this github project. They are easily editable according to what data is meant to be plotted. These are here just exemples of what can be done with the results of the simulations.

----------------------------------------
## Modules in FORTRAN 90
  * mod_param.f90
  * mod_random.f90
  * mod_fzero.f90
  * mod_proj_stereo.f90
  * mod_enveloppe.f90
  * mod_iso.f90
  * mod_diff.f90
  * mod_mie.f90
  * mie.f90

## Program in FORTRAN 90
### Isotropic scattering
  * iso_diff.f90
  * iso_diff_abs.f90

### Rayleigh scattering
  * rayl_diff.f90
  * rayl_RGB.f90
  * rayl_stereo_RGB.f90

### Mie scattering
  * mie_diff.f90
  * mie_diff_CloudySky.f90
  * mie_RGB.f90
  * mie_stereo_RGB.f90

### Other scattering models
  * hg_diff.f90
  * polyn_diff.f90

### Some tests
  * random_test.f90
  * random_test2.f90
  * env_eff.f90
  * test_bhmie.f90
  * test_bhmie_rb.f90


----------------------------------------
## Compilation with *make*
To compile and execute the most of the programs, this simple command is used:

    $ make program.exe && ./Bin/program.exe

The two exepctions are rayl_stereo_RGB.f90

    $ make rayl_RGB.exe && ./Bin/rayl_RGB.exe
    $ make rayl_stereo_RGB.exe && ./Bin/rayl_stereo_RGB.exe

and mie_stereo_RGB.f90

    $ make mie_RGB.exe && ./Bin/mie_RGB.exe
    $ make mie_stereo_RGB.exe && ./Bin/mie_stereo_RGB.exe


---------------------------------------
## Attached files 
### Post-processing routines in Python
  * deg_pol.py
  * diag_pol_diff.py
  * diag_pol_diff_bhmie.py
  * diag_pol_diff_RGB.py
  * diag_pol_mie.py
  * diag_pol_mie_RGB.py
  * env_eff.py
  * f_phase.py
  * iso_diff_abs.py
  * test_bhmie.py
  * test_bhmie_rb.py

### Makefile
  * Makefile


----------------------------------------
Plots and codes can also be found on GitLab <https://gitlab.obspm.fr/pbarrere/diffusion2>
