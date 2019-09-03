# SUMMARY:      extract_pfotran tecplot output files
# USAGE:        simple scripts to extract single h5 from pflotran output
# ORG:          Pacific Northwest National Laboratory
# AUTHOR:       Xuehang Song
# E-MAIL:       xuehang.song@pnnl.gov
# ORIG-DATE:    August-2019
# DESCRIPTION:
# DESCRIPTION-END

import glob
import re
import numpy as np

simu_dir = "/global/cscratch1/sd/xhsong/John_case_optim_5/"

# get file names
tec_files = glob.glob(simu_dir+"*.tec")[0:2]

# read first observation file
with open(tec_files[0], "r") as f:
    header = [x.replace('"', "").lower().lstrip()
              for x in f.readline().split(',')[:-1]]
data = np.genfromtxt(tec_files[0], skip_header=1)

# read rest of the files
for itec in tec_files[1:]:
    data = np.hstack((data, np.genfromtxt(itec, skip_header=1)))
    with open(itec, "r") as f:
        header.append([x.replace('"', "").lower().lstrip()
                       for x in f.readline().split(',')[:-1]])

header_0 = header

xxx = [[x] if x == header[0]
       else [x.split("(")[0].split(" ")[-2]] for x in header]
# x.split("(")[-2].split(")")[0],
# x.split("(")[-1].split(")")[0]]
