#!/bin/bash

/home/ubuntu/MATLAB/SEED/md5ToEC.awk /mnt/extra/SEED/MGMData/function_map /mnt/extra/SEED/MGMData/md5_protein_map > /home/ubuntu/MATLAB/SEED/md5ToEC.map
/home/ubuntu/MATLAB/SEED/mgmMd5ToEC.awk /home/ubuntu/MATLAB/SEED/md5ToEC.map /mnt/extra/SEED/MGMData/mgm4448044.3.650.protein.sims.head > /home/ubuntu/MATLAB/SEED/testEC.txt
grep -P "(\d|-)+\.(\d|-)+\.(\d|-)+\.(\d|-)+" /home/ubuntu/MATLAB/SEED/testEC.txt > /home/ubuntu/MATLAB/SEED/testEC2.txt
