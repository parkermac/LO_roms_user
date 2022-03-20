# README for LO_roms_user

### This repo is a place for user versions of ROMS code associated with the git repository LO_roms_source.

NOTE: Currently these directions are sort of specific to Parker, but I'll make them more generic soon.

---

(1) Put the ROMS source code on the machine you will be compiling on, e.g. klone:

git clone https://pmaccc@www.myroms.org/git/src LO_roms_source

Newer information about the ROMS repo and documentation seems to be going here: see https://www.myroms.org/doxygen/.

Below we will refer to the directory LO_roms_source as (*)

You can then do git pull anytime - it will ask for your ROMS password.

NOTE: when I do git pull I get this warning:
```
[pmacc@klone1 test0]$ git pull
warning: Pulling without specifying how to reconcile divergent branches is
discouraged. You can squelch this message by running one of the following
commands sometime before your next pull:

  git config pull.rebase false  # merge (the default strategy)
  git config pull.rebase true   # rebase
  git config pull.ff only       # fast-forward only
```
So I should figure out what these mean and then choose one.

---

#### test0

This is the upwelling test case that come with ROMS.  It is always the first thing you should try to run when moving to a new version of ROMS or a new machine.

I have created a few files to run it on klone:
- build_roms.sh is copied from (*)/ROMS/Bin/build_roms.sh and you can see the few edits I made by using diff.
- upwelling.h
- roms_upwelling.in
- klone_batch0.sh

In the directory test0 do:
```
pmsrun
./build_roms.sh
logout
```
This will spew several minutes of stuff to the screen and eventually result in the executable...

Then to run ROMS do:
```
./klone_batch.sh
```

NOTE: when I ran it today, 2022.03.19, it generated these messages but appeared to run fine.
```
[LOG_CAT_SBGP] libnuma.so: cannot open shared object file: No such file or directory
[LOG_CAT_SBGP] Failed to dlopen libnuma.so. Fallback to GROUP_BY_SOCKET manual.
```

---

#### npzd_banas

This folder started as copies of the Fennel code in LO_roms_source/ROMS/Nonlinear/Biology. It also includes External/bio_Fennel.in. Then these files were edited to retain the Fennel code, including NH4 and Chl variables, but modifying the parameters in the .in and the equations in fennel.h to reproduce the Banas/Siedlecki/Davis model as closely as possible, while allowing a separate NH4 pool. All changes in the .h code are denoted with:
```
! PM Edit
[new code]
! End PM Edit
```
The first test of using this will be the uu0kb executable.
