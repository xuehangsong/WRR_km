import json
import argparse
import numpy as np

simu_dir = "/media/sf_e/john/optim_5/"

# load jason file
with open(simu_dir+"obs.json", "r") as f:
    obs_dict = json.load(f)
header = obs_dict["header"]
data = np.array(obs_dict["data"])
time = np.array(obs_dict["time"])
time_unit = obs_dict["time_unit"]

# tracer name
tracer = ['total tracer_northbc [m]',
          'total tracer_river_middle [m]',
          'total tracer_river_north [m]',
          'total tracer_river_south [m]',
          'total tracer_river_upstream [m]']
max_tracer = 0.001

varis = np.unique([x[0] for x in header]).tolist()
loc_name = np.unique([x[1] for x in header]).tolist()

# wells
wells = ["399-2-2",
         "399-2-32",
         "399-1-21a"]

# get flux averaged tracer
simu_tracer = dict()
for iwell in wells:
    simu_tracer[iwell] = dict()

    # use flux in xy plane as weight
    qlx_index = [i for i, x in enumerate(
        header) if "qlx" in x[0] and iwell in x[1]]
    qly_index = [i for i, x in enumerate(
        header) if "qly" in x[0] and iwell in x[1]]
    weight = (data[:, qlx_index]**2+data[:, qly_index]**2)**0.5

    # extract weight
    for itracer in tracer:
        loc_index = [i for i, x in enumerate(
            header) if x[0] == itracer and iwell in x[1]]
        tracer_data = data[:, loc_index]
        simu_tracer[iwell][itracer] = np.append(0, np.average(
            tracer_data[1:, ],
            weights=weight[1:, ], axis=1))/max_tracer

# Get total river tracer
for iwell in wells:
    simu_tracer[iwell]["river"] = simu_tracer[iwell][
        "total tracer_river_middle [m]"] + \
        simu_tracer[iwell]["total tracer_river_north [m]"] + \
        simu_tracer[iwell]["total tracer_river_south [m]"] + \
        simu_tracer[iwell]["total tracer_river_upstream [m]"]
tracer = tracer+["river"]

# get well level
well_level = dict()
for iwell in wells:

    # locate pressure
    pressure_index = [i for i, x in enumerate(
        header) if "pressure" in x[0] and iwell in x[1]]
    liquid_data = data[:, pressure_index]
    liquid_data[liquid_data < 101325] = np.nan
    screen_ele = [float(header[x][3].split(" ")[-1])
                  for x in pressure_index]

    # caculate well well
    well_level[iwell] = np.nanmean((liquid_data-101325) /
                                   9.8068/998.2+screen_ele, axis=1)
