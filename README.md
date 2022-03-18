# README for LO_roms_user

### This repo is a place for user versions of ROMS code associated with the git repository LO_roms_source.

---

#### npzd_banas
This folder started as copies of the Fennel code in LO_roms_source/ROMS/Nonlinear/Biology. It also includes External/bio_Fennel.in. Then these files were edited to retain the Fennel code, including NH4 and Chl variables, but modifying the parameters in the .in and the equations in fennel.h to reproduce the Banas/Siedlecki/Davis model as closely as possible, while allowing a separate NH4 pool. All changes in the .h code are denoted with:
```
! PM Edit
[new code]
! End PM Edit
```
The first test of using this will be the uu0kb executable.
