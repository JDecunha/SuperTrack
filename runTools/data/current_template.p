(lp0
S'\xef\xbb\xbf#!/bin/bash\n\n#BSUB -W {walltime_request}\n#BSUB -o /rsrch3/home/imag_phy/jdecunha/SuperTrack/run_logfiles\n#BSUB -cwd /rsrch3/home/imag_phy/jdecunha/SuperTrack/runFiles\n#BSUB -q short\n#BSUB -M 165\n#BSUB -R rusage [mem=165]\n#BSUB -n 28\n#BSUB -u jdecunha@mdanderson.org\n#BSUB -J {job_name}\n\nsource /rsrch3/home/radphys_rsch/jdecunha/configure.sh\n\n{run_command}\n'
p1
aS'.lsf'
p2
a.