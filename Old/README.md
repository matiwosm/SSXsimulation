This README explains the useful files that I've worked on, and how to use Bridges and Dedalus for simulation and VAPOR for simulation.

Fist, you need to install Dedalus on Bridges. Transfer the cluster_install.sh file to your Bridges directory and run it. It will take care of everything. To actually run a simulation, I highly recommend you transfer

  dedalusBatch.job, mergeBatch.job, SSX_model_A, SSX_model_A_load, spheromak.py, merge1.py, and merge2.py

to your $SCRATCH directory on Bridges. Running them from my $HOME directory caused problems.

SSX_model_A.py will output field data and load_data. The field data (fields) contains what you may want to visualize, while the load_data is everything that Dedalus needs to load a run. There is a lot of redundancy between these two, so only using one might be an important optimization.

Use the dedalusBatch script to run an initial conditions run. Use the loadBatch script to run a loaded run.
The paths are specific to my Bridges account, so you will need to change the "source" and "cd" calls in both. cd to your $SCRATCH directory, and source from your $HOME directory (where you installed dedalus).

You may need to change the batch scripts to set the amount of time you run for, the number of cores you use, or the allocation you use. Changing the number of cores is a bit tricky, because the total number of cores needs to be a multiple of the mesh dimensions in SSX_model_A*.py. A full explanation of the batch scripts is pretty easy to understand and on the bridges website.

After a run, the output file folders will appear in your $SCRATCH directory.
A slightly more in-depth explanation of running can be found on the document "Using Bridges".

To visualize with VAPOR: download VAPOR (not VAPOR3) from ucar.
Before running you must source VAPOR and Dedalus.
To source VAPOR: "source /Applications/VAPOR/VAPOR.app/Contents/MacOS/vapor-setup.sh".
To source Dedalus "source /Users/dedalus/dedalus/bin/activate".
Your download locations may vary. You will need to do this every time you open a new terminal window.

Then run toVapor.py, taking in the merged h5 file you want to visualize the name of the file you want to output (exclude the file extension for this argument), and optionally dimension scales you want to visualize. (See below). This will output a .vdf file.
In VAPOR, go to Data->Load a Dataset into Current Session (or just command+D) and find your .vdf file. It will then load in the data.

To summarize: Things you probably need to to:
1. Install Dedalus on your Bridges account
2. Copy the necessary files to your $SCRATCH directory on Bridges
3. Modify the batch scripts to source from your account
4. Make any changes you need to the SSX_model_A*.py files (like changing resolution)
5. Install VAPOR (not VAPOR3!) on your local machine.

Files:

===dedalusBatch.job===
A batch script used on Bridges. Calls SSX_model_A.py to simulate spheromak evolution, then calls merge1.py and merge2.py to merge data together into an h5 file  named merged.h5.

===mergeBatch.job===
A batch script used on Bridges. Calls merge1.py and merge2.py to merge data into a merged file. Use in case dedalusBatch.job fails to merge (out of time). --cleanup is true by default.

===SSX_model_A.py===
The Dedalus script used to set up and simulate spheromak evolution. Requires spheromak.py to be in same directory.
Parameters--
Eta: Magnetic diffusivity
Mu: Viscosity*rho0
Kappa: Temperature diffusivity*rho0
Gamma: Ideal gas adiabatic index

===SSX_model_A_load.py===
The Dedalus script used to load a previous run and resume simulation.
Takes the name of an h5 file and will load the last time-step (you can change this in the code). Also takes an initial timestep.
It's just like SSX_model_A.py but it loads initial conditions instead of creating them.

===Spheromak.py===
Dependency of SSX_model_A.py. Sets spheromak initial current density. On top of the regular spheromak functions there's a tanh mask, because the functions are periodic, and we only want one spheromak.

===merge1.py===
Merges data from seperate processessors into coherent, time ordered file.
Usage:
    merge1.py <base_path> [--cleanup]
Options:
    --cleanup   Delete distributed files after merging

base_path is a directory containing the files to be merged

===merge2.py===
Merges data from seperate time-ordered h5 files (output of merge1.py), and creates one large, merged file.
Usage:
    merge2.py <joint_path> <set_path> [--cleanup]

Options:
    --cleanup   Delete distributed files after merging

joint_path is file name of finished, merged file
set_path is a directory containing files to be merged

***Merging***
The two merging scripts work, but there is now dedalus wrappers that run the merging:
dedalus merge_procs <base_path> [--cleanup]
dedalus merge_sets <joint_path> <set_path>... [--cleanup]

===toVapor.py===
Script to create a VAPOR-readable .vdf file from an h5 file produced by Dedalus. Uses netCDF as an intermediate format.
Pulls all variables output by Dedalus and captures time evolution if present.
Usage:
    toVapor.py <fileIn> <fileOutName> [<dimensionRatio>]

fileIn is the relevant h5 file
fileOutName is the name (no extension) of the resultant file (eg. "ssx")
dimensionRatio is an optional argument specifying the desired visualized dimensions. Syntax: "minZ:minY:minX:maxZ:maxY:maxX"

It's possible that the VAPOR call and Dedalus call will not work. In this case, source VAPOR and Dedalus directly by executing the command manually. Something like: source /Applications/VAPOR/VAPOR.app/Contents/MacOS/vapor-setup.sh and source /Users/dedalus/dedalus/bin/activate

FOR HELP regarding getting data into VAPOR more generally: see Eric Hester's comment here https://groups.google.com/forum/#!topic/dedalus-users/2tS6PS-zKLM

===particlePuncture.py===
Makes a puncture plot of a particle trajectory given information about where to take the slice, and with what precision.

===particleStaticFieldSpheromak.py===
Produces netCDF file with data about particle orbit and fields. Use ncdfvdfcreate and ncdf2vdf to make VAPOR-readable vdf files.
For more info, see the above link.

===cluster_install.sh===
Script to install Dedalus on Bridges. Trying to do it manually is truly a pain.
