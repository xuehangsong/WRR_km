import h5py as h5
import os
import numpy as np


input_name = "/mnt/e/john/vec/pflotran_2014.h5"
output_name = "/mnt/e/john/vec/pflotran_2014_center.h5"

input_file = h5.File(input_name, "r")
output_file = h5.File(output_name, "w")
groups = list(input_file.keys())
times = [x for x in groups if "Time" in x]


for igroup in ["Coordinates", "Provenance"]:
    print(igroup)
    input_file.copy(igroup, output_file)

for itime in times:
    print(itime)
    output_file.create_group(itime)
        
    output_file[itime].create_dataset(
        "Material_ID",
        data=input_file[itime]["Material_ID"])

    
    xflux_face = np.array(input_file[itime]["Liquid X-Flux Velocities"])
    xflux_update = np.concatenate(
        (xflux_face,np.expand_dims(xflux_face[-1,:,:],0)),
        axis=0)
    output_file[itime].create_dataset(
        "Liquid X-Velocity [m_per_h]",
        data=xflux_update)

    yflux_face = np.array(input_file[itime]["Liquid Y-Flux Velocities"])
    yflux_update = np.concatenate(
        (yflux_face,np.expand_dims(yflux_face[:,-1,:],1)),
        axis=1)
    output_file[itime].create_dataset(
        "Liquid Y-Velocity [m_per_h]",
        data=yflux_update)

    
    zflux_face = np.array(input_file[itime]["Liquid Z-Flux Velocities"])
    zflux_update = np.concatenate(
        (zflux_face,np.expand_dims(zflux_face[:,:,-1],2)),
        axis=2)    
    output_file[itime].create_dataset(
        "Liquid Z-Velocity [m_per_h]",
        data=zflux_update)
    
input_file.close()
output_file.close()

