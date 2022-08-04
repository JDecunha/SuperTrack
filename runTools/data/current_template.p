(lp0
S'\xef\xbb\xbf#!/bin/bash\n\n#BSUB -W {walltime_request}\n#BSUB -o /rsrch3/home/radphys_rsch/jdecunha/SuperTrack/run_logfiles\n#BSUB -cwd /rsrch3/home/radphys_rsch/jdecunha/SuperTrack/runFiles\n#BSUB -q gpu\n#BSUB -m gdragon004\n#BSUB -gpu num=1:gmem=16\n#BSUB -M 48\n#BSUB -R rusage[mem=48]\n#BSUB -n 10\n#BSUB -u jdecunha@mdanderson.org\n#BSUB -J {job_name}\n\nsource /rsrch3/home/radphys_rsch/jdecunha/configure_GPU.sh\n\n{run_command}\n'
p1
aS'.lsf'
p2
a.