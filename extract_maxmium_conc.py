# SUMMARY:      extract_maximum_conc.py
# USAGE:        simple scripts to extract single h5 from pflotran output
# ORG:          Pacific Northwest National Laboratory
# AUTHOR:       Xuehang Song
# E-MAIL:       xuehang.song@pnnl.gov
# ORIG-DATE:    sep-2019
# DESCRIPTION:
# DESCRIPTION-END

import glob
import numpy as np
import h5py as h5
import copy
import joblib

multi_dir = "/media/sf_e/john/optim_5/single_h5/"
results_dir = "/media/sf_e/john/optim_5/"

# list files and read times
snap_files = glob.glob(multi_dir+"Time*.h5")
times = [float(s.split(" ")[-2]) for s in snap_files]
snap_files = [snap_files[x] for x in np.argsort(times)]
times = [times[x] for x in np.argsort(times)]

# collect tracer data
river_pattern = dict()
for i_index, ifile in enumerate(snap_files[:]):
    # read data
    f = h5.File(ifile, "r")
#    print(list(f[list(f.keys())[0]]))
    print(ifile)
    tracer = f[list(f.keys())[0]]["Total_Tracer_river_middle [M]"][:] + \
        f[list(f.keys())[0]]["Total_Tracer_river_north [M]"][:] + \
        f[list(f.keys())[0]]["Total_Tracer_river_south [M]"][:] + \
        f[list(f.keys())[0]]["Total_Tracer_river_upstream [M]"][:]
    pressure = f[list(f.keys())[0]]["Liquid_Pressure [Pa]"][:]
    f.close()

    # get maxium river in each column
    max_river = copy.deepcopy(pressure)
    max_river[pressure < 101325] = 0
    max_river[pressure >= 101325] = 1
    max_river = max_river*tracer
    max_river = np.max(max_river, 2)*1000
    river_pattern[i_index] = max_river

river_pattern["times"] = times
joblib.dump(river_pattern, results_dir+"river_pattern.joblib")
