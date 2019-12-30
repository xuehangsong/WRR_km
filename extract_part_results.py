import h5py as h5
import os
import numpy as np


# define input and output locations
simu_dir = "/global/cscratch1/sd/chenxy/Xuehang/npt/simu_base/"
simu_dir = "/global/cscratch1/sd/renh686/Xuehang/npt/simu_homo/"
simu_dir = '/global/cscratch1/sd/hd09/Xuehang/npt/simu_base_smooth/'
simu_dir = "/global/cscratch1/sd/xhsong/npt/base/flow/"

# selected time
mar_first = np.sum([366, 365*3, 366, 365])*24+np.sum([31, 28, 31, 30]*24)
output_times = range(mar_first, mar_first+31*24, 1)

# files names
input_name = simu_dir + "pflotran_bigplume.h5"
output_name = simu_dir + "pflotran_2014.h5"

# remove files
if os.path.isfile(output_name):
    os.remove(output_name)

# open file and write
input_file = h5.File(input_name, "r")
output_file = h5.File(output_name, "w")

output_group = ["Coordinates", "Provenance"]
for itime in output_times:
    groupname = "Time:  "+"{0:.5E}".format(itime)+" h"
    output_group = np.append(output_group, groupname)
for igroup in output_group:
    print(igroup)
    input_file.copy(igroup, output_file)

input_file.close()
output_file.close()

#output_file = h5.File(case_name+"/"+output_name+".h5","r")
# if os.path.isfile(case_name+"/"+output_name+".h5"):
#     os.remove(case_name+"/"+output_name+".h5")
#groupnames = list(input_file.keys())


print("Hello World!")
