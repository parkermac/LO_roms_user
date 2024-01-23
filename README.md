# README for LO_roms_user

### This repo is a place for user versions of code used to compile ROMS code associated with the git repository LO_roms_source_git.

These notes are written for klone, but should work with only minor changes for mox.

---

#### Overview: klone and mox are UW supercomputers in the hyak system.

Here are examples of aliases I have on my mac ~/.bash_profile (equivalent to ~/.bashrc on the linux machines) to quickly get to my machines
```
alias klo='ssh pmacc@klone1.hyak.uw.edu'
alias mox1='ssh pmacc@mox1.hyak.uw.edu'
alias mox2='ssh pmacc@mox2.hyak.uw.edu'
alias pgee='ssh parker@perigee.ocean.washington.edu'
alias agee='ssh parker@apogee.ocean.washington.edu'
```
Note: klone1 is the same as klone.  If you just ssh to mox you end up randomly at either mox1 or mox2, which are the same machine except that they keep separate crontabs. I always ssh to mox1 to avoid confusion.

---

#### Tools to control jobs running on klone

`/gscratch/macc` is our working directory on both mox and klone because on hyak we are the "macc" group. I have created my own directory inside that: "parker", where all my code for running ROMS is stored.

When you have a job running on klone you can check on it using:
```
squeue -A macc
```
or on mox:
```
squeue -p macc
```
If you want to stap a running job, find the job ID (the number to the left in the squeue listing) and issue the command:
```
scancel [job ID]
```
Since your job will typically have been launched by a python driver you will also want to stop tat driver. Use "top" to find the associated job ID, and then use the "kill" command.

---

#### Getting resource info

`hyakstorage` will give info about storage on klone.  Use `hyakstorage --help` to get more info on command options. This is not yet working on mox.

To check on our disk allocation on mox you can also look in the file `/gscratch/macc/usage_report.txt` although this will be phased out soon.

`hyakalloc` will give info on the nodes we own.

**mox**: we own 148 cores (3 nodes with 28 cores each = 84, plus another 2 nodes with 32 cores each = 64). Unfortunately you cannot mix-and-match with these, you have to use nodes with the same number of cores on a given job.

**klone**: we own 600 cores (15 nodes with 40 cores each). We are allocated 1 TB of storage for each node, so 15 TB total.

---

#### Once you have gotten a klone account from our system administrator, you have two directories to be aware of.

**First directory:** In your home directory (~) you will need to add some lines to your .bashrc using vi or whatever your favorite command line text editor is.

Here is my .bashrc on klone:

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

#module load intel/oneAPI
LODIR=/gscratch/macc/local
#OMPI=${LODIR}/openmpi-ifort
NFDIR=${LODIR}/netcdf-ifort
NCDIR=${LODIR}//netcdf-icc
PIODIR=${LODIR}/pio
PNDIR=${PNDIR}/pnetcdf
export LD_LIBRARY_PATH=${PIODIR}/lib:${PNDIR}:${NFDIR}/lib:${NCDIR}/lib:${LD_LIBRARY_PATH}
#export LD_LIBRARY_PATH=${NFDIR}/lib:${NCDIR}/lib:${LD_LIBRARY_PATH}
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
alias cdLrs='cd /gscratch/macc/parker/LO_roms_source_git'
alias cdLod='cd /gscratch/macc/parker/LO_data'
alias pmsrun='srun -p compute -A macc --pty bash -l'
alias buildit='./build_roms.sh -j 10 < /dev/null > bld.log &'
alias mli='module load intel/oneAPI'
```

The section of aliases are what I use to help move around quickly.  You might want similar aliases but be sure to substitute the name of your working directory for "parker".

In particular you will need to copy and paste in the section with all the module and export lines.  These make sure you are using the right NetCDF and MPI libraries.

Note: I need to clean this up by getting rid of obsolete export calls, and setting the base working directory as a variable.

**Second directory:** The main place where you will install, compile, and run ROMS is your working directory:

**/gscratch/macc/[your directory name]**  We call this **(+)** below.

Note: Even though my username on klone is "pmacc" my main directory is "parker". This implies that there is less restriction in naming things on klone compared to apogee and perigee. I don't recall who set up my initial directory. Either David Darr or I did it.

---

#### Set up ssh-keygen to apogee

The LO ROMS driver system tries to minimize the files we store on hyak, because the ROMS output files could quickly exceed our quotas.  To do this the drivers (e.d. LO/driver/driver_roms3.py) uses scp to copy forcing files and ROMS output files to apogee or perigee where we have lots of storage.  Then the driver automatically deletes unneeded files on hyak after each day it runs.  To allow the driver to do this automatically you have to grant it access to your account on perigee or apogee, using the ssh-keygen steps described here.

Log onto klone1 and do:
```
ssh-keygen
```
and hit return for most everything.  However, you may encounter a prompt like this:
```
Enter file in which to save the key (/mmfs1/home/pmacc/.ssh/id_rsa):
/mmfs1/home/pmacc/.ssh/id_rsa already exists.
Overwrite (y/n)?
```
Looking [HERE](https://www.hostdime.com/kb/hd/linux-server/the-guide-to-generating-and-uploading-ssh-keys), I found out that id_rsa is the default name that it looks for automatically. You can name the key anything and then just refer to it when using ssh and etc. like:
```
ssh parker@apogee.ocean.washington.edu -i /path/to/ssh/key
```
In the interests of tidying up I chose to **overwrite** in the above. When I did this it asked for a passphrase and I hit return (no passphrase).

Then I did:
```
ssh-copy-id parker@apogee.ocean.washington.edu
```
(it asks for my apogee password)

And now I can ssh and scp from klone to apogee without a password, and on apogee it added a key with pmacc@klone1.hyak.local at the end to my ~/.ssh/authorized_keys.

Similarly, on klone there is now an entry in ~/.ssh/known_hosts for apogee.ocean.washington.edu.

So, in summary: for going from klone1 to apogee it added to:
- ~/.ssh/known_hosts on klone (boiler and mox1 are also there), and
- ~/.ssh/authorized_keys on apogee

Now I can run `ssh-copy-id` again for other computers, without having to do the `ssh-keygen` step.

Don't worry if things get messed up. Just delete the related entries in the .ssh files and start again. This is a good place to remind yourself that you need to be able to edit text files from the command line on remote machines, e.g. using vi.

---

#### Working from (+), clone the LO repo:
```
git clone https://github.com/parkermac/LO.git
```
Also clone your own LO_user repo.

---

#### Before you start using ROMS you should get a ROMS account.  See the first bullet link below.

Places for ROMS info:
- https://www.myroms.org/ Main page.  Click the "Register" tab to get an account.
- https://www.myroms.org/wiki/Documentation_Portal Main documentation portal.  Tons of info about everything.
- https://www.myroms.org/doxygen/ Alternate version of the portal.
- https://www.myroms.org/forum/index.php The ROMS Forum, where you can search for discussions of many issues, and post your own questions.  This is incredibly useful.  You will get speedy answers from the experts!

---

#### Get the ROMS source code

Then put the ROMS source code on klone, again working in (+).  Do this using git.  Just type this command.  This will create a folder LO_roms_source_git with all the ROMS code.
```
git clone https://github.com/myroms/roms.git LO_roms_source_git
```
You can bring the repo up to date anytime from inside LO_roms_source_git by typing `git pull`.

---

#### Next, create (on your personal computer) a git repo called LO_roms_user, and publish it to your account on GitHub.

Copy some of my code from https://github.com/parkermac/LO_roms_user into your LO_roms_user.  Specifically you want to get the folder "upwelling".

This is the upwelling test case that comes with ROMS.  It is always the first thing you should try to run when moving to a new version of ROMS or a new machine.

I have created a few files to run it on klone:
- `build_roms.sh` modified from LO_roms_source_git/ROMS/Bin.  **You need to edit line 152 so that MY_ROOT_DIR is equal to your (+).**
- `upwelling.h` copied from LO_roms_source_git/ROMS/Include.  No need to edit.
- `roms_upwelling.in` modified from LO_roms_source_git/ROMS/External.  **You will need to edit line 78 so that the path to varinfo.yaml points to (+).**
- `klone_batch0.sh` created from scratch.  **You will need to edit line 24 so that RUN_DIR points to (+).**

After you have edited everything on your personal computer, push it to GitHub, and clone it to (+) on klone.

---

#### Now you are ready to compile and run ROMS (in parallel) for the first time!

Working on klone **in the directory LO_roms_user/upwelling**, do these steps, waiting for each to finish, to compile ROMS:
```
srun -p compute -A macc --pty bash -l
```
or on mox:
```
srun -p macc -A macc --pty bash -l
```
The purpose of this is to log you onto one of our compute nodes because in the hyak system you are supposed to compile on a compute node, leaving the head node for stuff like running our drivers and moving files around.  You should notice that your prompt changes, now showing which node number you are on. Any user in the LiveOcean group should be able to use this command as-is because "macc" refers to our group ownership of nodes, not a single user.  Note that in my .bashrc I made an alias `pmsrun` for this hard-to-remember command.

Then before you can do the compiling on klone (ignore this step on mox) you have to do:
```
module load intel/oneAPI
```
I have this aliased to `mli` in my .bashrc.

Then to actually compile you do:
```
./build_roms.sh -j 10 < /dev/null > bld.log &
```

This will take about six minutes, spew a lot of text to bld.log, and result in the executable `romsM`. It also makes a folder `Build_roms` full of intermediate things such as the .f90 files that result from the preprocessing of the original .F files. I have this aliased as `buildit` in my .bashrc.

The `-j 10` argument means that we use 10 cores to compile, which is faster.  Note that each node on klone had 40 cores.

On occasion I have a problem where keyboard input (like hitting Return because you are impatient) causes the job to stop.  That is why I added the `< /dev/null` thing to this command.

#### >>> After compiling is done, DO NOT FORGET TO: <<<
```
logout
```
to get off of the compute node and back to the head node. If I forget to do logout and instead try to run ROMS from the compute node it will appear to be working but not make any progress.

Then to run ROMS do (from the klone head node, meaning after you logged out of the compute node):
```
sbatch -p compute -A macc klone_batch0.sh
```
or if you are working on mox the command is:
```
sbatch -p macc -A macc mox_batch0.sh
```
This will run the ROMS upwelling test case on 4 cores.  It should take a couple of minutes.  You can use the < > & things to not have to wait for it to finish.

If it ran correctly it will create a log file roms_log.txt and NetCDf output: roms_[his, dia, avg, rst].nc

---

#### Running things by cron

These are mainly used by the daily forecast but can also be helpful for checking on long hindcasts and sending you an email. See LO/driver/crontabs for my current versions.  These are discussed more in LO/README.md.

---

## LO Compiler Configurations

Below we list the current folders where we define LO-specific compiling choices.  The name of each folder refers to [ex_name] in the LO run naming system.  Before compiling, each contains:
- `build_roms.sh` Which can be copied directly from your `upwelling` folder, without need to edit.
- `[ex_name].h` This has configuration specific compiler flags.  You can explore the full range of choices and their meanings in `LO_roms_source/ROMS/External/cppdefs.h`.
- `fennel.h` if this is a run with biology.

NOTE: to run any of these, or your own versions, you have to make the LO_data folder in (+) and use scp to get your grid folder from perigee or apogee.

**NOTE: the ex_name can have numbers, but no underscores, and all letters MUST be lowercase.**

---

## CURRENT

#### x4b

Like x2b but modified by Aurora Leeson (her meV00) to increase the light attenuation by a factor of three for the Salish Sea. In her tests she also used the full TRAPS forcing (WWTP bug fixed) and MPDATA for bio tracer advection in the dot_in. This should be the default code for the long hindcast and for whenever we update the forecast. 2023.11.05

---

#### xa0

Meant for an analytical run. Basically identical to x4b but with the atmospheric forcing set to zero, and biology turned off. This replaces uu1k which did the same thing previously.

---

## OBSOLETE BUT RECENT

Mostly I call these _obsolete_ becasue they use the somewhat older ROMS we had from svn, and they rely on varinfo.yaml in LO_roms_source_alt. But some of them are being used for the current daily forecast, specifically x2b and xn0b.

#### uu0mb

This is a major step in the ROMS update process.
- It uses the near-latest version of ROMS.
- It is meant to be run using `driver_roms3.py`. Please look carefully at the top of that code to see all the command line arguments.
- It uses the PERFECT_RESTART cpp flag. The leads to a smoother run and fewer blow-ups. It also means that it no longer writes an ocean_his_0001.nc file. This would be identical to the 0025 file from the previous day. This change is accounted for in `Lfun.get_fn_list()`.
- It incorporates rain (EMINUSP).
- It assumes that the forcing was created using driver_forcing3.py. This uses the new organizational structure where forcing is put in a [gridname] folder, not [gridname_tag] or [gtag].
- See `LO/dot_in/cas6_v00_uu0mb` for an example dot_in that runs this.

For the bio code:
- It uses my edited version of the fennel bio code, which I keep in `LO_roms_source_alt/npzd_banas`.
- We correct att and opt in the bio code, the match BSD as written.
- Better atm CO2.

---

#### uu0m

This is just like uu0mb except without biology.

---

#### x0mb

Like uu0mb but with the rOxN* ratios set back to the original Fennel values (instead of the larger Stock values). Also some changes to the benthic remin: (i) fixed a bug in the if statement to test if the aerobic flux would pull DO negative, and (ii) a simpler handling of denitrification from benthic flux, but ensuring it does not pull NO3 negative.

I introduced a new name here because I had been recycling uu0mb to many times!

---

#### x1b

Like x0mb but I edited the bio code to include the "optimum uptake" form of nutrient limitation for NH4. It was already in NO3. Created 2023.04.08.

It is poor design to have the bio code in a separate folder. For example, if I now recompiled x0mb I would get code that reflected x1mb. **So I am going to put fennel.h in this folder and then set
MY_ANALYTICAL_DIR=${MY_PROJECT_DIR} in `build_roms.sh`.**

I am also dropping the "m" for mox. There unless I was running parallel forecasts on both mox and klone (as I was once) there is no reason for this.

---

#### x2b

This starts from the fennel.h code in x1b and modifies it so the the benthic flux conforms more closely to Siedlecki et al. (2015) except with the necessary change that remineralization goes into NH4 instead of NO3. Denitrification still comes out of NO3. I also turn off the light limitation in Nitrification.

---

#### x3b

An experiment using the fennel.h code from x2b but modifying light attenuation to be what it is in the current forecast. Note that the current forecast has bugs in this part of the code that make it different from Davis et al. (2014) as written.

---

#### uu1k

This is much like uu0mb except it drops the cppdefs flags associated with atm forcing and biology. This makes it useful for analytical runs that don't have atm forcing. Note carefully the ANA flags used in the cpp file. Like uu0mb, it makes use of forcing files that use the new varinfo.yaml to automate the naming of things in the NetCDF forcing files (the "A0" sequence).

---

#### xn0

Designed to run a nested model. Omits tidal forcing. No biology. Otherwise based on x2b.

---

#### xn0b

Designed to run a nested model. Omits tidal forcing. Has biology from x2b. Otherwise based on xn0.

---

## OBSOLETE AND VERY OLD

### These notes are only relevant to the old ROMS installation used in the LiveOcean (not LO) system

#### From David Darr: klone requires only modest changes to Linux-ifort_mox.mk, renamed Linux-ifort_klone.mk (in LiveOcean_roms/LO_ROMS/Compilers).  The only difference actually is that the two instances of:
```
NC_CONFIG ?= nc-config
```
are now:
```
NC_CONFIG ?= nf-config
```

#### Steps to compile:

On klone:

```
cd LiveOcean_roms/LO_ROMS
srun -p compute -A macc --pty bash -l
make clean
make -f /gscratch/macc/parker/LiveOcean_roms/makefiles/[ex_name]/makefile
```
Then `logout` to get back to the usual shell.  You have to do this because the `srun` command logged you onto one of the compute nodes.

On mox the steps are only slightly different. The `compute` in the srun command is `macc` on **mox**:
```
cd LiveOcean_roms/LO_ROMS
srun -p macc -A macc --pty bash -l
make clean
make -f /gscratch/macc/parker/LiveOcean_roms/makefiles/[ex_name]/makefile
```
