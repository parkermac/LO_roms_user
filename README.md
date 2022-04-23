# README for LO_roms_user

### This repo is a place for user versions of ROMS code associated with the git repository LO_roms_source.

These notes are written for klone, but should work with only minor changes for mox.

---

(1) Put the ROMS source code on the machine you will be compiling on, e.g. klone:

```
git clone https://[your ROMS username]@www.myroms.org/git/src LO_roms_source

```
Newer information about the ROMS repo and documentation seems to be going here: see https://www.myroms.org/doxygen/.

Below we will refer to the directory LO_roms_source as (*)

You can then do git pull anytime - it will ask for your ROMS password.

By issuing the command while in the directory on klone
```
git config pull.ff only
```
it means you will always use fast-forward for reconciling branches. If you don't do this you get a warning message from klone. This seems to work with my usual one-way git workflow: edit on my mac => commit and push to cloud => then pull to klone.

---

(2) This directory, LO_roms_user, is a git repo that you clone to klone at the same level as LO, LO_roms_source, etc. It is a place for copies of ROMS code that I need to edit. There are a few categories:
- The compiler instructions, copied from (*)/Compilers.
- The files needed for making a specific executable, like test0, uu0k, etc., see below.
- The custom bio code, which goes in npzd_banas. This has its own README.

For all the edited code you can use diff to compare it to the original version in (*).

---

#### test0

This is the upwelling test case that come with ROMS.  It is always the first thing you should try to run when moving to a new version of ROMS or a new machine.

I have created a few files to run it on klone:
- build_roms.sh from (*)/ROMS/Bin
- upwelling.h from (*)/ROMS/Include
- roms_upwelling.in from (*)/ROMS/External
- klone_batch0.sh created from scratch

In the directory test0 do:
```
pmsrun
./build_roms.sh < /dev/null > bld.log &
logout
```
This will take a few minutes and result in the executable romsM. It also makes a folder Build_roms full of intermediate files such as the .f90 that result from the preprocessing of the original .F files. You issue these three commands anytime you are compiling. The /dev/null input is to avoid having the process stopped by any keyboard input.  The large amount of standard output will end up in bld.log, and & just escapes you back to shell.

pmsrun is an alias on klone, which appears as a line in my .bashcr
```
alias pmsrun='srun -p compute -A macc --pty bash -l'

```
The purpose of pmsrun is to log you onto one of our compute nodes because in the hyak system you are supposed to compile on a compute node, leaving the head node for stuff like running our drivers and moving files around.

The build_roms.sh script orchestrated the actual compiling, like you may have done in the past with a makefile.

When you logout after compiling it leaves you back on the klone head node.

Then to run ROMS do:
```
./klone_batch.sh
```

If it ran correctly it will create a log file roms_log.txt and NetCDf output: roms_[his, dia, avg, rst].nc

NOTE: when I ran it today, 2022.03.19, it generated these messages, but appeared to run fine.
```
[LOG_CAT_SBGP] libnuma.so: cannot open shared object file: No such file or directory
[LOG_CAT_SBGP] Failed to dlopen libnuma.so. Fallback to GROUP_BY_SOCKET manual.
```

---

#### uu0k

This is meant to exactly reproduce a physics-only version of the current LiveOcean forecast (cas6_v0_u0kb).

NOTE: the naming of the folder is important.  We currently rely on a repeated first letter, "uu" in this case, so that LO/driver/driver_roms1.py knows this is using the LO_roms_user system and not the old LiveOcean_roms system.

Because this is an LO-style run, it has more external gizmos than test0:
- It is meant to be run by driver_roms1.py
- It will look on apogee for forcing files, e.g. in LO_output/cas6_v0
- It will look on klone for a history file from the previous day to start from (in this case we can copy one of the current cas6_v0_u0kb files)
- It needs a dot_in instance such as cas6_v0_uu0k

Then when all these are in place I run in LO/driver with a command like:
```
python3 driver_roms1.py -g cas6 -t v0 -x uu0k -r backfill -0 2021.11.10 -np 400 -N 40 --move_his False --short_roms True --get_forcing False < /dev/null > uu0k.log &
```

---

#### uu1k

This is much like uu0k except it drops the cppdefs flags associated with atm forcing.  This makes it useful for analytical runs that don't have atm forcing.

I run this with `driver_roms2.py` because that is where I fixed a bug in the start_type=new logic.

```
python3 driver_roms2.py -g ae0 -t v0 -x uu1k -r backfill -s new -0 2020.01.01 -1 2020.01.02 -np 40 -N 40 < /dev/null > ae.log &
```
