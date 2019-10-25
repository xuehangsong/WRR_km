# SUMMARY:      map_particle_tracking.py
# USAGE:        simple scripts to extract single h5 from pflotran output
# ORG:          Pacific Northwest National Laboratory
# AUTHOR:       Xuehang Song
# E-MAIL:       xuehang.song@pnnl.gov
# ORIG-DATE:    Oct-2019
# DESCRIPTION:  map 1.6km model domain results to 3km model domain
# DESCRIPTION-END

import numpy as np
import h5py as h5


def proj_to_model(origin, coord):
    """
    convert projected coordinates to model coordinates
    """
    output = np.zeros(coord.shape)
    output[:, 0] = (coord[:, 0]-origin[0])*np.cos(origin[2]) + \
        (coord[:, 1]-origin[1])*np.sin(origin[2])
    output[:, 1] = (coord[:, 1]-origin[1])*np.cos(origin[2]) - \
        (coord[:, 0]-origin[0])*np.sin(origin[2])
    return(output)


def model_to_proj(origin, coord):
    """
    convert model coordinates to projected coordinates 
    """
    output = np.zeros(coord.shape)
    output[:, 0] = origin[0]+coord[:, 0] * \
        np.cos(origin[2])-coord[:, 1]*np.sin(origin[2])
    output[:, 1] = origin[1]+coord[:, 0] * \
        np.sin(origin[2])+coord[:, 1]*np.cos(origin[2])
    return(output)


origin_file = "/media/sf_e/john/particles/mapped/data.h5"
model_origin_1 = [594186, 115943, 15/180*np.pi]
model_origin_2 = [593000, 114500, 0]

data = h5.File(origin_file, "r+")
dgrp_names = [x for x in list(data.keys()) if "." in x]
for iname in dgrp_names:
    print(iname)
    coord_data = np.array(data[iname]["Coordinates"])
    # coord_data[:, 0:2] = proj_to_model(
    #     model_origin_2, model_to_proj(model_origin_1, coord_data[:, 0:2]))
    coord_data[:, 0:2] = model_to_proj(model_origin_1, coord_data[:, 0:2])
    data[iname]["Coordinates"][...] = coord_data
    # change particle index
    if iname == dgrp_names[0]:
        y_index = np.argsort(np.argsort(coord_data[:, 1])[::-1])
        # del data["Particle index"]
        # data.create_dataset("Particle index", data=y_index)
        data["Particle index"][...] = y_index
        data["Particle start time"][...] = y_index.astype("float")
data.close()
# coord = np.ones((5, 2))
# coord = proj_to_model(model_origin_1, coord)
# coord = model_to_proj(model_origin_1, coord)
