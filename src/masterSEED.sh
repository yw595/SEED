#!/bin/bash

/home/fs01/yw595/MATLAB/SEED/src/md5ToEC.awk /home/fs01/yw595/MATLAB/SEED/input/MGMData/function_map /home/fs01/yw595/MATLAB/SEED/input/MGMData/md5_protein_map > /home/fs01/yw595/MATLAB/SEED/md5ToEC.map
/home/fs01/yw595/MATLAB/SEED/src/mgmMd5ToEC.awk /home/fs01/yw595/MATLAB/SEED/md5ToEC.map /home/fs01/yw595/MATLAB/SEED/input/MGMData/mgm4448044.3.650.protein.sims > /home/fs01/yw595/MATLAB/SEED/testEC.txt
grep -P "(\d|-)+\.(\d|-)+\.(\d|-)+\.(\d|-)+" /home/fs01/yw595/MATLAB/SEED/testEC.txt > /home/fs01/yw595/MATLAB/SEED/testEC2.txt
