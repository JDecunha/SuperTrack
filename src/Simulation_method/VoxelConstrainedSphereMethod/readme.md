# VoxelConstrainedSphereMethod
## Geometry 
The geometry of the VoxelConstrainedSphereMethod, consists of a series of spherical targets confined to a cubic voxel. The end user specifies the target sphere size as well as the desired voxel side halflength in the input .ini file. The targets span the voxel without any spacing between adjacent targets and without spacing between them and the edge of the voxel.

The voxel side length itself is malleable, and if the sphere size you've requested does not evenly tile the voxel, the voxel side length will reduce itself until an integer number of targets span it. There is no guarantee whether there will be a target at (0,0,0). If an even number of targets spans the length of the voxel, then there will be no target at the origin. If an odd number of targets spans the voxel, then there will be a target situated at the origin, with an equal number of targets on each side of it. 

## How Tracks are Handled and Superimposed

The tracks for this method must be developed from the sister project MicroTrackGenerator. Given the side length of the cube you requested an X-Y shift which spans the side length of the cube will be generated. The purpose is to allow the track to originate from any location on the negative Z surface of the cube (the tracks move towards positive Z).

## Output

Currently the only available output information is lineal energy in keV/micron, output through the standard SuperTrack Histogram functionality. See the documentation for Histogram in this project for further information.

## Potential Future Functionality if Desired by Community:

**Other output data (energy imparted and specific energy):** A functor could be passed to the existing CUDA kernels to allow for different values to be passed out. However since the targets are all the same size, you can already convert to specific energy or energy imparted in post-processing.

**List of Volume-Edep pair as output:** For trouble shooting or visualization purposes it might be interesting to output the volumes and associated energy deposition which occured in each sphere. 
Currently only a Histogram is possible.

**Place start of track explicitly on voxel edge:** See above section on word of caution to see what this is referring to. 

**Efficiency improvements:**
  
  - *Malloc GPU memory only once, based on the longest track to be analyzed*. Currently every time a new track is analyzed, new malloc calls are made on the GPU for the track, the histogram, the lists that get compacted etc. If we looked at all the tracks and found the biggest one, we could just malloc once. This would gain some performance benefits.
  - *Implement track pre-fetching, and only run the CUDA context from one thread*. My knowledge of CUDA has evolved from when I started developing this project. I am now of the belief that I could have all of the CPU threads but one focused on "feeding" tracks to the GPU, and a single thread manages all the GPU calls. Currently the CUDA context gets handed off rapidly between the different threads at some performance cost. More than that, my fetching of tracks is not parallelized at all. I suspect this will yield measurable performance improvements.
  - *Migrate to std::async*. This topic would be a modification to SuperTrackManager as a whole. This would be preferred because debugging the code with GDB would be simplified. I also believe the overhead of std::async is much reduced compared to the current use of fork as well. 
  - *Implement clear distinction between main thread, and sub-process parts of the method*. Currently many things are done on each sub-process that the VoxelConstrainedSphere method is produced on (i.e. the macro file is read in nThreads number of times). This class should be re-structured such that it can be created on the main thread, and some of its values are copied over to the new process. This would be more efficient.
  - *Figure out a method to increase the persistency of the GPU histogram*. Currently the GPU histogram is accumulated after every oversample. Surely we could do that accumulation less often to improve performance somehow right?
  - *Allow single precision floating point operations.* This would just be implemented by templating a bunch of the kernels, but also the track would have to be read in as float as well.

## Notes on the CUDA Kernel Algorithms

The following comments concern themselves with the inner workings of the algorithm. If you only intend to be an end-user of the code you can safely ignore this section. If you intend to extend or modify the code then this will be of interest to you.

**FilterInScoringBox:** A kernel which applies a random shift to the x,y coordinates of the track. Checks which events are in the box. Performs stream compation on the events in the box. The stream-compaction is accomplished with shared atomics so the code, while efficient, is hard to understand if you're new to CUDA. 
Outside of the shared atomics there is nothing I think would be very puzzling about this kernel.

**FilterTrackInSphere:** A kernel which determines which randomly shifted track step points are in a sphere. Performs stream compaction using shared atomics on those events. The tricky part of this kernel is where the distance from nearest sphere is calculated. Essnetially, if you take the floor of position/targetDiameter you get the index of the sphere nearest the point. You then add 0.5 to that index to bring you to the center of the sphere. Multiply by the diameter to get the location of that sphere center. Subtract the center of the sphere from the step position to get the distance in that axis.

**ScoreTrackInSphere:** A kernel which, with the compacted list from FilterTrackInSphere, generates a VolumeEdepPair (a series of pairs of volume identifiers and the amount of edep there). This volume edep pair still has to be reduced to combine all of the energy depositions which occured in the same sphere. The only tricky part of this code is calculating the indices. Basically you can get the index from the position/targetDiameter, however we subtract the position from the edge of the box to make it so that all positions are taken relative to the box edge. Taking the position relative to the box edge makes all of this positions positive, so that all of the indices are positive too. This ensures that the sphere at the extreme negative edge of the voxel has index (0,0,0).
