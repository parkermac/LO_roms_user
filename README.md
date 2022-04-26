# README for LO_roms_user

### This repo is a place for user versions of code used to compile ROMS code associated with the git repository LO_roms_source.

These notes are written for klone, but should work with only minor changes for mox.

---

Once you have gotten a klone account from our system administrator, you have two directories to be aware of.

**First directory:** In your home directory (~) you will need to add some lines to your .bashrc using vi or whatever your favorite command line text editor is.  Here is what mine looks like:
```
# .bashrc

# Source global definitions
if [ -f /etc/bashrc ]; then
	. /etc/bashrc
fi

# User specific environment
if ! [[ "$PATH" =~ "$HOME/.local/bin:$HOME/bin:" ]]
then
    PATH="$HOME/.local/bin:$HOME/bin:$PATH"
fi
export PATH

module load intel/oneAPI
LODIR=/gscratch/macc/local
#OMPI=${LODIR}/openmpi-ifort
NFDIR=${LODIR}/netcdf-ifort
NCDIR=${LODIR}//netcdf-icc
export LD_LIBRARY_PATH=${NFDIR}/lib:${NCDIR}/lib:${LD_LIBRARY_PATH}
export PATH=/gscratch/macc/local/netcdf-ifort/bin:$PATH
export PATH=/gscratch/macc/local/netcdf-icc/bin:$PATH
#export PATH=/gscratch/macc/local/openmpi-ifort/bin:$PATH


# Uncomment the following line if you don't like systemctl's auto-paging feature:
# export SYSTEMD_PAGER=

# User specific aliases and functions
alias cdpm='cd /gscratch/macc/parker'
alias cdLo='cd /gscratch/macc/parker/LO'
alias cdLu='cd /gscratch/macc/parker/LO_user'
alias cdLoo='cd /gscratch/macc/parker/LO_output'
alias cdLor='cd /gscratch/macc/parker/LO_roms'
alias cdLru='cd /gscratch/macc/parker/LO_roms_user'
alias cdLrs='cd /gscratch/macc/parker/LO_roms_source'
alias cdLra='cd /gscratch/macc/parker/LO_roms_source_alt'
alias cdLod='cd /gscratch/macc/parker/LO_data'
alias pmsrun='srun -p compute -A macc --pty bash -l'
```
In particular you will need to copy and paste in the section with all the module and export lines.  These make sure you are using the right NetCDF and MPI libraries.

The section of aliases are what I use to help move around quickly.  You might want similar aliases but be sure to substitute the name of your working directory for `parker`.

**Second directory:** The main place where you will install, compile, and run ROMS is your working directory: /gscratch/macc/[whatever].  We call this (+) below.

---

Working from (+), clone the LO, and LO_roms_source_alt repos:
```
git clone https://github.com/parkermac/LO.git
```
```
git clone https://github.com/parkermac/LO_roms_source_alt.git
```
Also clone your own LO_user repo.

---

Before you get the ROMS code repo you need to get a ROMS account.  See the first bullet link below.

Places for ROMS info:
- https://www.myroms.org/ Main page.  Click the "Register" tab to get an account.
- https://www.myroms.org/wiki/Documentation_Portal Main documentation portal.  Tons of info about everything.
- https://www.myroms.org/doxygen/ Alternate version of the portal.
- https://www.myroms.org/forum/index.php The ROMS Forum, where you can search for discussions of many issues, and post your own questions.  This is incredibly useful.  You will get speedy answers from the experts!

---

Then put the ROMS source code on klone, again working in (+).  Do this using svn (subversion, similar to git).  Just type this command (substituting in your ROMS username, no []).  This will create a folder LO_roms_source with all the ROMS code.
```
svn checkout --username [your ROMS username] https://www.myroms.org/svn/src/trunk LO_roms_source

```
You can bring the repo up to date anytime from inside LO_roms_source by typing `svn update`.  It will ask for your ROMS password.  You can also ask for a specific revision of the code, and many other things.  Type `svn help` to find out more.

NOTE: the LO_roms_source_alt repo that you cloned above has versions of bits of the ROMS source code that we have edited, such as the modified biogeochemical code.  We keep these in a separate repo to they do not conflict with the original source code.  We are about to point our compiler to look in this directory for what it needs.

---

Next, create (on your personal computer) a git repo called LO_roms_user, and publish it to your account on GitHub.

Copy some of my code from https://github.com/parkermac/LO_roms_user into your LO_roms_user.  Specifically you want to get the folder "upwelling".

This is the upwelling test case that comes with ROMS.  It is always the first thing you should try to run when moving to a new version of ROMS or a new machine.

I have created a few files to run it on klone:
- `build_roms.sh` modified from LO_roms_source/ROMS/Bin.  **You need to edit line 112 so that MY_ROOT_DIR is equal to your (+).**
- `upwelling.h` modified from LO_roms_source/ROMS/Include.  No need to edit.
- `roms_upwelling.in` modified from LO_roms_source/ROMS/External.  **You will need to edit line 77 so that the path to varinfo.yaml points to (+).**
- `klone_batch0.sh` created from scratch.  **You will need to edit line 23 so that RUN_DIR points to (+).**

After you have edited everything on your personal computer, push it to GitHub, and clone it to (+) on klone.

---

**Now you are ready to compile and run ROMS (in parallel) for the first time!**

Working on klone **in the directory LO_roms_user/upwelling**, do these steps, waiting for each to finish, to compile ROMS:
```
srun -p compute -A macc --pty bash -l
```
The purpose of this is to log you onto one of our compute nodes because in the hyak system you are supposed to compile on a compute node, leaving the head node for stuff like running our drivers and moving files around.  You should notice that your prompt changes, now showing which node number you are on. Any user in the LiveOcean group should be able to use this command as-is because "macc" refers to our group ownership of nodes, not a single user.  Note that in my .bashrc I made an alias `pmsrun` for this hard-to-remember command.

Then to do the actual compiling, do this:

```
./build_roms.sh -j 10
```

This will take about six minutes, spew a lot of text to your screen, and result in the executable `romsM`. It also makes a folder `Build_roms` full of intermediate things such as the .f90 files that result from the preprocessing of the original .F files.

The `-j 10` argument means that we use 10 cores to compile, which is faster.  Note that each node on klone had 40 cores.

If you don't want all the screen output, do something like.
```
./build_roms.sh -j 10 > bld.log &
```
This redirects (>) the screen output to a file bld.log, and escapes to shell (&).

On occasion I have things where keyboard input (like hitting Return because you are impatient) causes the job to stop.  If you find this happening you can add a `< /dev/null` thing like this:
```
./build_roms.sh -j 10 < /dev/null > bld.log &
```
After compiling is done, do:
```
logout
```
to get off of the compute node and back to the head node.

Then to run ROMS do:
```
./klone_batch.sh
```
This will run the ROMS upwelling test case on 4 cores.  It should take a couple of minutes.  You can use the < > & things to not have to wait for it to finish.

If it ran correctly it will create a log file roms_log.txt and NetCDf output: roms_[his, dia, avg, rst].nc

---

## LO Compiler Configurations

Below we list the current folders where we define LO-specific compiling choices.  The name of each folder refers to [ex_name] in the LO run naming system.  Before compiling, each contains only:
- `build_roms.sh` Which can be copied directly from your `upwelling` folder, without need to edit.
- `[ex_name].h` This has configuration specific compiler flags.  You can explore the full range of choices and their meanings in `LO_roms_source/ROMS/External/cppdefs.h`.

NOTE: to run any of these, or your own versions, you have to make the LO_data folder in (+) and use scp to get your grid folder from perigee or apogee.

---

#### uu0k

This is meant to exactly reproduce a physics-only version of the current LiveOcean forecast (cas6_v0_u0kb), but updated to the newest ROMS, and using forcing files that use the new varinfo.yaml to automate the naming of things in the NetCDF forcing files.

Because this is an LO-style run, it has more external gizmos than test0:
- It is meant to be run by driver_roms2.py
- It will look on apogee for forcing files in LO_output/cas6_v0
- It will look on klone for a history file from the previous day to start from (in this case we can copy one of the current cas6_v0_u0kb files)
- It needs a dot_in instance such as cas6_v0_uu0k

Then when all these are in place I run in LO/driver with a command like:
```
python3 driver_roms2.py -g cas6 -t v0 -x uu0k -r backfill -0 2021.11.10 -np 40 -N 40 --move_his False --short_roms True --get_forcing False < /dev/null > uu0k.log &
```
Check out the comments in `driver_roms2.py` to see what all the commend line arguments do.

---

#### uu1k

This is much like uu0k except it drops the cppdefs flags associated with atm forcing.  This makes it useful for analytical runs that don't have atm forcing.

Example command to run it:

```
python3 driver_roms2.py -g ae0 -t v0 -x uu1k -r backfill -s new -0 2020.01.01 -1 2020.01.02 -np 40 -N 40 < /dev/null > ae.log &
```
