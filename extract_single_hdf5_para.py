# SUMMARY:      extract_single_hdf5_para.py
# USAGE:        simple scripts to extract single h5 from pflotran output
# ORG:          Pacific Northwest National Laboratory
# AUTHOR:       Xuehang Song
# E-MAIL:       xuehang.song@pnnl.gov
# ORIG-DATE:    August-2019
# DESCRIPTION:
# DESCRIPTION-END


import h5py as h5
import os
import numpy as np
import glob
from shutil import copyfile
import multiprocessing as mp


def create_single_h5(ifile):
    print(ifile)
    print(os.path.isfile(ifile))
    datafile = h5.File(ifile, "r")
    if ifile == single_simu[0]:
        groups = list(datafile.keys())  # [0:2]
    else:
        groups = [T for T in list(datafile.keys()) if "Time:" in T]  # [0:2]
    for igroup in groups:
        if "Time:" in igroup:
            jgroup = ("Time:  " +
                      "{:.5E}".format(float(igroup.split(" ")[2]) +
                                      time_increment)+" h.h5")
        else:
            jgroup = igroup+".h5"
        print(jgroup)
#        jgroup=jgroup+ifile.split("/")[-2]
        output_file = h5.File(multi_dir+jgroup, "w")
        datafile.copy(igroup, output_file)
        output_file.close()
    datafile.close()


multi_dir = "/global/cscratch1/sd/xhsong/John_case_optim_5/single_h5/"
single_simu = [
    "/global/cscratch1/sd/xhsong/John_case_optim_5/pflotran_bigplume-000.h5",
    "/global/cscratch1/sd/xhsong/John_case_optim_5/pflotran_bigplume-001.h5",
    "/global/cscratch1/sd/xhsong/John_case_optim_5/pflotran_bigplume-002.h5",
    "/global/cscratch1/sd/xhsong/John_case_optim_5/pflotran_bigplume-003.h5",
    "/global/cscratch1/sd/xhsong/John_case_optim_5/pflotran_bigplume-004.h5",
    "/global/cscratch1/sd/xhsong/John_case_optim_5/pflotran_bigplume-005.h5"]

time_increment = 0
ncore = 1
pool = mp.Pool(processes=ncore)
pool.map(create_single_h5, single_simu)
