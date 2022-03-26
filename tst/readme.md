# SuperTrack Tests

A series of tests using the gtest framework have been implemented for SuperTrack. Specifically they test for the correct output of the Histogram class and SuperTrack as whole.

## How to run the tests

**The tests  must be run in independent processes.** There are two ways to achieve this.

1.) Through invoking the dedicated runTestsSeparate.sh script.

2.) Through filtering for the specific test:
> ../build/test/SuperTrack_tst --gtest_filter=*A3

<details>
  <summary>Developer's Note: reason why tests must be run separately</summary>
  
  >Currently there are two identified issues which require tests be run separately. 
  >
  >The SuperTrack threading paradigm involves the launching of multiple processes using fork(). This can run into issues if the Histogram tests are run first, because the Histogram tests create a CUDA context in the main process. Later when fork() is called and SuperTrack tries to access the CUDA context (generated in the main process) from the sub-processes this leads to a segfault. 
  >
  >Secondly. The SuperTrackManager is implemented as a singleton and has state which can extend between tests. The primary issue that you'll run into with this, is the SuperTrackManager will try to analyze tracks which have been fed to it from other tests, in addition to the tracks being analyzed in the current test. For both of these reasons tests must be run in independent processes.
</details>

Tests are run from the SuperTrack/tst folder, attempting to invoke the tests from another directory will cause a test failure as SuperTrack relies on relative directories to reach the test macro files.

## Generating the Test Tracks
A version of CERN ROOT newer than v. 6.24 is recommended. From a bash terminal invoke the command `root tstTrackGenerator.cc` from the SuperTrack/tst/testTracks folder.

## Appendix: List of Tests 
### Histogram Tests
**1.) HandlesVectorOfOneInDifferentVolumes**
Verifies the ability of Histogram to process a VolumeEdepPair consisting of 1 eV each in a different volume.

**2.) HandlesVectorOfFiveInDifferentVolumes**
Verifies the ability of Histogram to process a VolumeEdepPair consisting of 5 eV each in a different volume.

**3.) HandlesVectorOfOnesInSameVolume**
Verifies the ability of Histogram to process a VolumeEdepPair consisting of 1 eV all in the same volume.

**4.) HandlesVectorOfThousandInSameVolume**
Verifies the ability of Histogram to process a VolumeEdepPair consisting of 1000 eV all in the same volume.

**Same as above but logarithmic histogram:**

**5.) HandlesVectorOfOneInDifferentVolumesLogarithmic**

**6.) HandlesVectorOfFiveInDifferentVolumesLogarithmic**

**7.) HandlesVectorOfOnesInSameVolumeLogarithmic**

**8.) HandlesVectorOfThousandInSameVolumeLogarithmic**

### SuperTrack Tests
**Test A1: A single energy deposition event**
Verifies the ability of SuperTrack to quantify a single energy deposition event leading to a 1 keV/um lineal energy deposit in a single 1 um diameter sphere.

**Test A2: Energy deposition events in a line**
Verifies the abiliy of SuperTrack to quantify a series of energy deposition events, in each of the 3000, 1 um diameter spheres which span a 3 mm voxel.

**Test A3: Energy deposition events in a plane**
Similar to above, except an entire plane of spheres is spanned. Each containing a singular energy deposition event.

**Test A4: Testing of oversamples and multiple tracks**
Performs a test which resembles Test A2 (energy depositions in a line), except 5 identical tracks are considered which are sampled 1000 times each.

**Test A5: Multiple energy depositions in each volume**
Builds upon test A4, except rather than a singular energy deposition point there are 1000 depositions in every sphere. Together they accumulate to a 1000 keV/um lineal energy deposit. The test analyzes 5 identical tracks which are sampled 1000 times each.

**Test A6: Testing the bounding geometry**
A series of energy deposition points are placed on the X and Y edges of the 3 mm bounding voxel. Testing along the z-edge is not conducted because the program does not filter tracks which are out of range in the z-axis, as it assumes that tracks have been cleaved on the z-axis prior to use in SuperTrack. This test is passed if none of the energy depositions on the edge of the box appear in the output.

**Test A7: Verifying the randomess of track generation**
This test involves sampling of a track which traverses a straight line along the z-axis. The track is sampled 10000 times each with a random shift in the x-y plane. The number of energy deposition events which land inside the spheres compared to outside is used to make a prediction of pi. If the prediction of pi is accurate for 10000 samples, within ~5% I consider this evidence that the random track shifts are working appropriately.
