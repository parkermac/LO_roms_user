# README for LO_roms_user

### This repo is a place for user versions of code used to compile ROMS code associated with the git repository LO_roms_source_git.

These notes are written for klone.

---

#### Overview: klone is a UW supercomputer in the hyak system. It is what we use for all our ROMS simulations.

Here are examples of aliases I have on my mac ~/.bash_profile (equivalent to ~/.bashrc on the linux machines) to quickly get to my machines
```
alias klo='ssh pmacc@klone1.hyak.uw.edu'
alias pgee='ssh parker@perigee.ocean.washington.edu'
alias agee='ssh parker@apogee.ocean.washington.edu'
```
Note: klone1 is the same as klone.

---

#### Tools to control jobs running on klone

There is excellent documentation of the UW hyak system, for example starting here:

https://hyak.uw.edu/docs/compute/start-here

I encourage yo to explore the tabs on the left of that page to answer any questions you have.

`/gscratch/macc` is our working directory on klone because on hyak we are the "macc" group. I have created my own directory inside that: "parker", where all my code for running ROMS is stored. TO DO: How do I create directories for new members of our group?

Here are a few useful commands.

When you have a job running on klone you can check on it using:
```
squeue -A macc
```
-A refers to the account (macc in our case). We also have resources in "coenv".

You could also use other command line arguments to get specific info:

-p [compute, cgu-g2, ckpt-g2] for different "partitions" which is a fancy word meaning which type of computer is your job running on.

-u [a UW NetID] to see info on the job(s) for a specific user

If you want to stop a running job, find the job ID (the number to the left in the squeue listing) and issue the command:
```
scancel [job ID]
```
Since your job will typically have been launched by a python driver you will also want to stop that driver. Use "top" to find the associated job ID, and then use the "kill" command.

---

#### Getting resource info

`hyakstorage` will give info about storage on klone.  Use `hyakstorage --help` to get more info on command options.

`hyakalloc` will give info on the nodes we have priority access to.

More specifics about nodes we can use:

**-A macc -p compute**: These are the original klone nodes. We own 600 cores (15 nodes with 40 cores each). We are allocated 1 TB of storage for each node, so 15 TB total.

**-A macc -p cpu-g2**: These are the new klone nodes. We own 480 cores (15 "slices" with 32 cores each). Each node consists of 6 slices, so we own 2.5 nodes. The advantage of running on the these slices is that it is easier for the scheduler to allocate resources because they are all on one node. They are also faster. **Currently 6 slices are reserved for the daily forecast system.**

**-A coenv -p cpu-g2**: We own 320 cores (10 slices with 32 cores each). This is in a separate account because of the history of how they were purchased.

**-p ckpt-g2**: These are cpu-g2 nodes that are available to anyone in the UW system, no -A account needed. They have proven to be useful even for long, large runs.

--- 

### Notes on using these resources

#### There are some specific ways to make better use of these nodes, and knowing this will make your jobs start faster, and run more reliably.

**For the old "compute" nodes** you would typically need to ask for more than one node. So for example to use 200 cores you would run with a command like:
```
python3 driver_roms4.py -np 200 -N 40 --cpu_choice compute --group_choice macc
```
And it will assign your job to run on 4 compute nodes (confirm using squeue).

**For the new cpu-g2 nodes** the strategy is different, because each node has 192 cores. **So in almost all cases you really want to run on one node.** If you used the node strategy that we used for compute nodes then your job would be forced to be spread out accoss separate nodes, which makes it much harder for the sbatch system to allocate and run. You could do this using driver_roms4.py with a command like:
```
python3 driver_roms4.py -np 160 -N 160 --cpu_choice cpu-g2 --group_choice coenv
```
At least I think this should work - I haven't tested it yet. This should run your job on most of 1 node (again conform using squeue).

I have also created a special driver that is tailored for the daily forecast, which I use with a command like:
```
python3 driver_roms5.py --group_choice macc --cpu_choice cpu-g2 -tpn 192
```
This uses a different batch script that is hard-coded to just use one node, and to use special sbatch commands: --exclusive and --mem=0. These force the scheduler to only allow you to use a node, and to use all the memory on the node.  I also created driver_roms5a.py that is for more general use, that does not have these special spatch commands.

### In general everyone except for Parker should use `driver_roms4.py` using the guidelines above, and stick to the coenv/cpu-g2 or macc/compute or ckpt-g2 unless I specifically give you the okay to work on macc/cpu-g2.

---

#### Once you have gotten a klone account from our system administrator, you have two directories to be aware of.

**First directory:** In your home directory (~) you will need to add some lines to your .bashrc using vi or whatever your favorite command line text editor is.

Here is my .bashrc on klone as of 6/28/2025:

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
NCDIR=${LODIR}/netcdf-icc
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
LOd=/gscratch/macc/parker/LO/driver
alias cdpm='cd /gscratch/macc/parker'
alias cdLo='cd /gscratch/macc/parker/LO'
alias cdLu='cd /gscratch/macc/parker/LO_user'
alias cdLoo='cd /gscratch/macc/parker/LO_output'
alias cdLor='cd /gscratch/macc/parker/LO_roms'
alias cdLru='cd /gscratch/macc/parker/LO_roms_user'
alias cdLrs='cd /gscratch/macc/parker/LO_roms_source_git'
alias cdLod='cd /gscratch/macc/parker/LO_data'
alias pmsrun='srun -p compute -A macc --pty bash -l'
alias pmsrun2='srun -p cpu-g2 -A macc --pty bash -l'
alias buildit='./build_roms.sh -j 10 < /dev/null > bld.log &'
alias buildit_dev='./build_roms.sh -j 10 -b develop < /dev/null > bld.log &'
alias mli='module load intel/oneAPI'

# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/gscratch/macc/parker/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/gscratch/macc/parker/miniconda3/etc/profile.d/conda.sh" ]; then
        . "/gscratch/macc/parker/miniconda3/etc/profile.d/conda.sh"
    else
        export PATH="/gscratch/macc/parker/miniconda3/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<
```

The section of aliases are what I use to help move around quickly.  You might want similar aliases but be sure to substitute the name of your working directory for "parker".

In particular you will need to copy and paste in the section with all the module and export lines.  These make sure you are using the right NetCDF and MPI libraries.

The conda part was added automatically when I set up a python environment on klone. At this point, however you DO NOT need to create a new python environment on klone. The one that is already there is enough to do all our model runs.

TO DO: I need to clean this up by getting rid of obsolete export calls, and setting the base working directory as a variable.

**Second directory:** The main place where you will install, compile, and run ROMS is your working directory:

**/gscratch/macc/[your directory name]**  We call this **(+)** below.

Note: Even though my username on klone is "pmacc" my main directory is "parker". This implies that there is less restriction in naming things on klone compared to apogee and perigee. I don't recall who set up my initial directory. Either David Darr or I did it.

---

#### Set up ssh-keygen to apogee

The LO ROMS driver system tries to minimize the files we store on hyak, because the ROMS output files could quickly exceed our quotas.  To do this the drivers (e.d. LO/driver/driver_roms4.py) uses scp to copy forcing files and ROMS output files from apogee or perigee where we have lots of storage.  Then the driver automatically deletes unneeded files on hyak after each day it runs.  To allow the driver to do this automatically you have to grant it access to your account on perigee or apogee, using the ssh-keygen steps described here.

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
- ~/.ssh/known_hosts on klone, and
- ~/.ssh/authorized_keys on apogee

Now I can run `ssh-copy-id` again for other computers, without having to do the `ssh-keygen` step.

Don't worry if things get messed up. Just delete the related entries in the .ssh files and start again. This is a good place to remind yourself that you need to be able to edit text files from the command line on remote machines, e.g. using vi.

---

#### Working from (+), clone the LO repo:
```
git clone https://github.com/parkermac/LO.git
```
Also clone your own LO_user repo. Note that you do not have to install the "loenv" python environment. All the code we run on klone is designed to work with the default python installation that is already there.

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
- `build_roms.sh` modified from LO_roms_source_git/ROMS/Bin.  **You need to edit line 173 so that MY_ROOT_DIR is equal to your (+).**
- You may also need to issue the command `chmod u+x build_roms.sh` so that you have permisssion to execute that script.
- `upwelling.h` copied from LO_roms_source_git/ROMS/Include.  No need to edit.
- `roms_upwelling.in` modified from LO_roms_source_git/ROMS/External.  **You will need to edit line 76 so that the path to varinfo.yaml points to (+).**
- `klone_batch0.sh` created from scratch.  **You will need to edit line 24 so that RUN_DIR points to (+).**

After you have edited everything on your personal computer, push it to GitHub, and clone it to (+) on klone.

---

#### Now you are ready to compile and run ROMS (in parallel) for the first time!

Working on klone **in the directory LO_roms_user/upwelling**, do these steps, waiting for each to finish, to compile ROMS:
```
srun -p compute -A macc --pty bash -l
```
The purpose of this is to log you onto one of our compute nodes because in the hyak system you are supposed to compile on a compute node, leaving the head node for stuff like running our drivers and moving files around.  You should notice that your prompt changes, now showing which node number you are on. Any user in the LiveOcean group should be able to use this command as-is because "macc" refers to our group ownership of nodes, not a single user.  Note that in my .bashrc I made an alias `pmsrun` for this hard-to-remember command. I also have `pmsrun2` to use "-p cpu-g2", the next-generation nodes. Don't use these without checking with Parker; some are reserved for the daily forecast system!

Then before you can do the compiling on klone you have to do:
```
module load intel/oneAPI
```
I have this aliased to `mli` in my .bashrc.

Then to actually compile you do:
```
./build_roms.sh -j 10 < /dev/null > bld.log &
```

This will take about ten minutes, spew a lot of text to bld.log, and result in the executable `romsM`. It also makes a folder `Build_romsM` full of intermediate things such as the .f90 files that result from the preprocessing of the original .F files. I have this aliased as `buildit` in my .bashrc.

The `-j 10` argument means that we use 10 cores to compile, which is faster.  Note that each node on klone had 40 cores (or 32 if we were using the cpu-g2 partition).

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
This will run the ROMS upwelling test case on 4 cores.  It should take a couple of minutes.  You can add the < > & things to the sbatch command line not have to wait for it to finish.

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

Naming conventions: There are not formal naming conventions, but I typically start with a letter like "x", or "xn" if it is for a nested (no tides) case. Then a number like "4" to give some indication of where it is in our development. I typically append "b" if the run includes biology.

---

## CURRENT

#### x11b

This is the current primary forecast executable, as of 6/22/2025. It is like x10ab but with updated ROMS version (4.3, as of 2025.06.12) and compiled using the -b develop branch (probably not needed if you did git pull recently in LO_roms_source_git), and defining OMEGA_IMPLICIT. This has the Harcourt turbulence improvements as the default advection scheme. I also created a new build_roms.sh, to keep it up to date with the one in the ROMS repo.

#### x11ab

Like x11ab but defining averages. Intended for daily saves.

#### xn11b

Like x11b but without tides. Intended for nested runs.

---

#### x10ab

Like x4b but with 50% burial of organic particulate N and C in the Salish Sea. It also saves averages. It is designed to run only saving two history files per day and an average file. No PERFECT_RESTART. This was a big step in the development leading from x4b to x11b. It uses new flags in the .h to turn on and off the bgc edits that apply only to the Salish Sea.

#### x4b

The old default code used for the long hindcast and daily forecast cas7_t0_x4b. The fennel.h code has lines to increase the light attenuation by a factor of three for the Salish Sea. It allows for vertical point sources (like wastewater treatment plants) which requires a more recent ROMS repo (~January 2024). It uses MPDATA for bio tracer advection in the dot_in.

#### xn4b

Like x4b but without tides, for nested runs.

---

#### x4a

Like x4b but no biology, no perfect restart, but defining AVERAGES. This is for the new experiment with GLORYS forcing (May 2025) but could be used for any physics-only experiment where we want to run fast.

---

#### x4

Like x4b but without biology, for testing physics changes like tidal forcing.

---

#### x4tf

Like x4 but with the tidal tractive force turned on.

---

#### xa0

Meant for an analytical run. Basically identical to x4b but with the atmospheric forcing set to zero, and biology turned off. This replaces uu1k which did the same thing previously.

---

## OBSOLETE

Mostly I call these _obsolete_ becasue they use the somewhat older ROMS we had from svn, and they rely on varinfo.yaml in LO_roms_source_alt. But some of them are being used for the current daily forecast, specifically xn0b.

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
