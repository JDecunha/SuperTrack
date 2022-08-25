# SuperTrack (v1.0)
An application to superimpose tracks of ionizing radiation on to volumes for the determination of microdosimetric spectra.
## Installation and Requirements
The required software and libraries to compile SuperTrack include:
1. CERN ROOT 6.24 or greater. Must be built from source rather than installed from a binary.
2. CMake 3.3 or greater.
3. An NVIDIA CUDA capable card with compute capability > 6.0 and CUDA Runtime version > 10.1.
4. Any required dependancies of CERN ROOT.

SuperTrack is built by:
```
#from SuperTrack directory, in BASH, on Linux
mkdir build
cd build
cmake ..
make
```
## Running the Software
The software is typically run from a subdirectory of the project main directory. Any command which takes the following form will be appropriate for running the software:
```
../build/SuperTrack ../macros/[macro name and extension here]
```
In the instance that a path to a macro file is not provided the software will search for a macro file at "../macros/default.ini". If a default .ini macro file is not at that location the program will throw an error.
### Input Parameters
**Command Line**

There are no optional command line parameters. The only command line parameter is the mandatory macro file parameter described in the above section.

**Macro File**

**[Run] Type** 

Specifies the type of track library to be importing in to SuperTrack. Currently, "Folder" is the only option.

**[Run] Path**

Path to the input folder of tracks.
 
**[Run] FirstFile**

If using Run type = Folder, this specifies the first track in the folder (in alphanumeric order) to start analyzing. If this parameter and LastFile are both set to -1 all files in the folder will be analyzed.

**[Run] LastFile**

If using Run type = Folder, this specifies the last track in the folder (in alphanumeric order) to start analyzing. If this parameter and FirstFile are both set to -1 all files in the folder will be analyzed.

**[Run] Threads**

Number of CPU threads to run the program on.

**[Run] Oversamples**

Number of times to re-sample each track.

**[Simulation] Method**

Describe the geometry and track sampling technique for superimposing. Currently only "VoxelConstrainedSphere" has been implemented.

Note that the src/Simulation_method/VoxelConstrainedSphereMethod folder contains its own README.

**[VoxelConstrainedSphere] ScoringRegionHalfLength**

Set the half-length of the scoring "voxel" in nanometers.

**[VoxelConstrainedSphere] ScoringSphereDiameter**

Set the scoring sphere diameter in nanometers.

**[VoxelConstrainedSphere] SuggestedCudaBlocks**

This is the number of CUDA blocks the geometry releted kernels will run on. Tweak this parameter for your hardware. 

**[VoxelConstrainedSphere] SuggestedCudaThreads**

This is the number of CUDA threads the geometry releted kernels will run on. Tweak this parameter for your hardware.

**[Histogram] Type**

Set either linear ("lin") or logarithmic ("log") bins for output histogram.

**[Histogram] NBins**

Number of histogram bins.

**[Histogram] BinLower**

Lower end of the smallest histogram bin.

**[Histogram] BinUpper**

Upper end of the greatest histogram bin.

**[Histogram] SuggestedCudaAccumulateBlocks**

This is the number of CUDA blocks the histogram releted kernels will run on. Tweak this parameter for your hardware. Compared to the blocks and threads used for the SimulationMethod, the histogram block and thread configuration affects performance much less.

**[Histogram] SuggestedCudaAccumulateThreads**

This is the number of CUDA threads the histogram releted kernels will run on. Tweak this parameter for your hardware. Compared to the blocks and threads used for the SimulationMethod, the histogram block and thread configuration affects performance much less.

**[Output] Path**

Path to where the output microdosimetric spectra should be saved.

**[Output] NamePrefix**

This string will preface the output file names.

### Example macro file

A macro file for SuperTrack to superimpose 50 MeV proton tracks on a 5 um cubic side length voxel is provided:
```
; Example macro file for SuperTrack

[Run]
Type = Folder
Path = [Path to your tracks]/MicroTrackGenerator/output/proton/50.0MeV/ #this is just an example, you can store your track library anywhere
FirstFile = -1 ; #-1 is the default case, will analyze whole folder
LastFile = -1 ;
Threads = 10
Oversamples = 10

[Simulation]             
Method = VoxelConstrainedSphere

[VoxelConstrainedSphere]
;Geometry parameters
ScoringRegionHalfLength = 2.5e3    ; side length in nanometers
ScoringSphereDiameter = 10        ; sphere diameter in nanometers

;Meta parameters
SuggestedCudaBlocks = 256
SuggestedCudaThreads = 32


[Histogram]
Type = lin 
NBins = 3000
BinLower = 0.01
BinUpper = 300.0

;Meta Parameters
SuggestedCudaAccumulateBlocks = 64
SuggestedCudaAccumulateThreads = 32

[Output]
Path = ../output/proton/
NamePrefix = proton_50.0MeV_10nm_diameter_x10_oversamples_
```

## Using the runTools

Using the methods described in the previous section, a user can run the application by making their own macro file and invoking the application from the command line or a script. We highly recommend however, that users planning to run MicroTrackGenerator use the bundled runTools to automatically create their own macro files and run scripts. The runTools are a Python2 based command line interface application. A Python2 installation with numpy is all that should be required to make use of the runTools. The runTools are designed to rapidly generate 1.) macro files and 2.) run files for computing cluster schedulers. The runTools are invoked on the command line from the runTools folder by: 
```
python2 SuperTrack_runtools.py
```
The commands available from the main shell of the runTools are _build_, _configure_, _help_, and _quit_. The use of the _configure_ and _build_ commands will be explained in the next sections.
### Configuring the runTools
When first using the runTools the _configure_ command should be invoked. When _configure_ is called, the terminal will prompt you to provide a path to a file containing a template of the run file for your cluster. You will have to make a template appropriate for the computing cluster you intend to run the software on. Bundled with the software, in the runTools/templates folder is an example template which is for use on the MD Anderson Cancer Center Seadragon computing cluster. The template is as follows:
```
#!/bin/bash

#BSUB -W {walltime_request}
#BSUB -o /rsrch3/home/radphys_rsch/jdecunha/SuperTrack/run_logfiles
#BSUB -cwd /rsrch3/home/radphys_rsch/jdecunha/SuperTrack/runFiles
#BSUB -q gpu
#BSUB -m gdragon[003:004]
#BSUB -gpu num=1:gmem=16
#BSUB -M 48
#BSUB -R rusage[mem=48]
#BSUB -n 10
#BSUB -u jdecunha@mdanderson.org
#BSUB -J {job_name}

source /rsrch3/home/radphys_rsch/jdecunha/configure_GPU.sh

{run_command}

```
In addition to all of the necessary commands for your scheduler, the template should include: `{walltime_request}` `{job_name}` and `{run_command}` in the appropriate locations. The configure command will also ask you for the file extension you desire for any generated runfiles. Make sure to include the "." preceeding your file extension.
### Building with the runTools
Once your template has been created and _configure_ has been invoked, you are ready to call _build_. The build procedure will prompt you to give the required inputs one at a time, in order to generate the macro and run files. _build_ allows you to generate a single macro and run file pair, or a series of macro files for a single geometry and many different particle energies. This is useful if you desire to generate a library of microdosimetric spectra at various energies. If you call _build_ then _series_ again after already building a series, you will be given the option to change the target size, for generating microdosimetric spectra in various different target sizes.

## Software License

Those wishing to use, modify, or reproduce this software must contact Joseph DeCunha at jdecunha@mdanderson.org to discuss an appropriate collaboration agreement. Copyright is claimed by Joseph M. DeCunha, 2022. All rights not expressly granted under this license are reserved.

This software is provided by the copyright holder "As is" and any express or implied warranties, including, but not limited to, implied warranties of merchantability, of satisfactory quality, and fitness for a particular purpose or use are disclaimed. The copyright holder makes no representation that the software and modifications thereof, will not infringe any patent, copyright, trade secret or other proprietary right.

The copyright holder shall have no liability for direct, indirect, special, incidental, consequential, exemplary, or punitive damages of any character including, without limitation, procurement of substitute goods or services, loss of use, data or profits, or business interruption, however caused and on any theory of contract, warranty, tort (including negligence), product liability or otherwise, arising in any way out of the use of this software, even if advised of the possibility of such damages.
