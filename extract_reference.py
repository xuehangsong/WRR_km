# SUMMARY:      extract_pfotran tecplot output files
# USAGE:        simple scripts to extract single h5 from pflotran output
# ORG:          Pacific Northwest National Laboratory
# AUTHOR:       Xuehang Song
# E-MAIL:       xuehang.song@pnnl.gov
# ORIG-DATE:    August-2019
# DESCRIPTION:
# DESCRIPTION-END

import glob
import numpy as np
import json
import argparse


def read_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('-dir',
                        '--simu_dir',
                        type=str,
                        default="")
    args = vars(parser.parse_args())
    return(args)


# read parser
args = read_args()
for ikey, ivari in args.items():
    exec(ikey + '=ivari')


# get file names
tec_files = glob.glob(simu_dir+"*.tec")

# read first observation file
with open(tec_files[0], "r") as f:
    header = [x.replace('"', "").lower().lstrip()
              for x in f.readline().split("\n")[0].split(',')]
data = np.genfromtxt(tec_files[0], skip_header=1)

# # read rest of the files
for itec in tec_files[1:]:
    data = np.hstack((data, np.genfromtxt(itec, skip_header=1)))
    with open(itec, "r") as f:
        header += [x.replace('"', "").lower().lstrip()
                   for x in f.readline().split("\n")[0].split(',')]

# break header
header = [[x] if x == header[0]
          else [x.split(x.split("(")[0].split(" ")[-2])[0].strip(),
                x.split("(")[0].split(" ")[-2].strip(),
                x.split("(")[-2].split(")")[0].strip(),
                x.split("(")[-1].split(")")[0].strip()] for x in header]

# remove dumplicated time column
time_index = [i for i, x in enumerate(header) if x == header[0]]
time = data[:, 0]
time_unit = header[0][0].split("[")[1].split("]")[0]
header = [x for i, x in enumerate(header) if i not in time_index]
data = np.delete(data, time_index, axis=1)

# form the dict
obs_dict = {"header": header,
            "time_unit": time_unit,
            "time": time.tolist(),
            "data": data.tolist()}

# dump jason file
with open(simu_dir+"obs.json", "w") as f:
    json.dump(obs_dict, f, indent=4)
