# SUMMARY:      spatial cluster.py
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
import matplotlib.pyplot as plt
import random
from datetime import datetime, timedelta
from scipy.spatial import distance
from scipy.cluster import hierarchy
from scipy.cluster.hierarchy import dendrogram
from matplotlib.pyplot import cm
from fastdtw import fastdtw
from shapely.geometry import MultiPoint
from descartes import PolygonPatch
from matplotlib.patches import Patch
from matplotlib.lines import Line2D
from pylab import *
import matplotlib.gridspec as gridspec
from scipy.spatial.distance import squareform


def DTWDistance(s1, s2):
    """
    The basic to caculate norm2 DTW.
    reference:"https://nbviewer.jupyter.org/github/alexminnaar/time-series-classification-and-clustering/tree/master/"
    For large dataset, we can use the dtw package in python
    Xuehang Song
    07/19/2019
    """
    DTW = {}

    for i in range(len(s1)):
        DTW[(i, -1)] = float('inf')
    for i in range(len(s2)):
        DTW[(-1, i)] = float('inf')
    DTW[(-1, -1)] = 0

    for i in range(len(s1)):
        for j in range(len(s2)):
            dist = (s1[i]-s2[j])**2
            DTW[(i, j)] = dist + min(DTW[(i-1, j)],
                                     DTW[(i, j-1)], DTW[(i-1, j-1)])
    return (DTW[len(s1)-1, len(s2)-1])**0.5


def dtw():
    # dtw (slow version)
    # ri = random.choices(np.arange(len(river_tracer[0])), k=50)
    # dist_spatial = [DTWDistance(river_tracer[:, ri[i]], river_tracer[:, ri[j]])
    #              for i in range(len(ri)) for j in range(len(ri))]
    ri = random.choices(np.arange(len(river_tracer[0])), k=400)
    dist_spatial_tri = [[fastdtw(river_tracer[:, ri[i]], river_tracer[:, ri[j]])[0]
                         for i in range(j+1)] for j in range(len(ri))]
    dist_spatial = np.zeros((len(ri), len(ri)))
    for i, j in enumerate(dist_spatial_tri):
        dist_spatial[i][0:len(j)] = j
    dist_spatial = dist_spatial+np.transpose(dist_spatial)
    dist_spatial = distance.squareform(
        dist_spatial, force="tovector", checks=True)


def batch_delta_to_time(origin, x, time_format, delta_format):
    nx = len(x)
    y = []
    for ix in range(nx):
        if delta_format == "hours":
            temp_y = origin + timedelta(hours=x[ix])
        elif delta_format == "days":
            temp_y = origin + timedelta(days=x[ix])
        elif delta_format == "minutes":
            temp_y = origin + timedelta(minutes=x[ix])
        elif delta_format == "weeks":
            temp_y = origin + timedelta(weeks=x[ix])
        elif delta_format == "seconds":
            temp_y = origin + timedelta(seconds=x[ix])
        elif delta_format == "microseconds":
            temp_y = origin + timedelta(microseconds=x[ix])
        elif delta_format == "milliseconds":
            temp_y = origin + timedelta(milliseconds=x[ix])
        else:
            print("Sorry, this naive program only solve single time unit")
        y.append(temp_y.strftime(time_format))
    y = np.asarray(y)
    return(y)


def batch_time_to_delta(origin, x, time_format):
    nx = len(x)
    y = []
    for ix in range(nx):
        temp_y = abs(datetime.strptime(
            x[ix], time_format) - origin).total_seconds()
        y.append(temp_y)
    y = np.asarray(y)


multi_dir = "/media/sf_e/john/optim_5/single_h5/"
results_dir = "/media/sf_e/john/optim_5/"
img_dir = "/media/sf_e/john/optim_5/figures/"
material_file = "/media/sf_e/john/optim_5/300A_material.h5"
river_stage = "/media/sf_e/john/optim_5/DatumH_River_filtered_6h_321.txt"
sampling_file = "/media/sf_e/john/optim_5/Sample_Data_2015_U.csv"
wells_dir = "/media/sf_e/john/wells/"

colors = cm.Set1(np.linspace(0, 1, 9))[
    np.array([6, 2, 1, 3, 4, 0, 7, 8, 5])]
#colors = np.array(["brown", "green", "blue", "purple", "orange", "red"])

# cluster information from Chris
well = ["399-1-10A",
        "399-1-21A",
        "399-2-1",
        "399-2-2",
        "399-2-23",
        "399-2-3",
        "399-2-32",
        "399-2-33",
        "399-2-5",
        "399-2-7",
        "399-3-19",
        "399-3-9",
        "399-4-9"]
well_name = ["1-10A",
             "1-21A",
             "2-01",
             "2-02",
             "2-23",
             "2-03",
             "2-32",
             "2-33",
             "2-05",
             "2-07",
             "3-19",
             "3-09",
             "4-09"]
well_cluster = [1, 3, 4, 4, 2, 4, 2, 2, 2, 2, 3, 4, 1]
well_coord = np.array([[594346.57, 116733.99],
                       [594160.78, 116183.88],
                       [594467.25, 116121.22],
                       [594385.8, 116282.6],
                       [594272.27, 116057.59],
                       [594377.49, 116220.4],
                       [594284.62, 116195.13],
                       [594263.56, 116145.76],
                       [594287.74, 116068.8],
                       [594235.2, 116084.49],
                       [594071.94, 116030.22],
                       [594504.51, 115917.93],
                       [594537.9, 115741.5]])
well_marker = ["o",  "^", "X", "*", "s", "v", "<", ">"]
line_type = ["-", ":"]

# define key date
date_origin = datetime.strptime("2010-01-01 00:00:00", "%Y-%m-%d %H:%M:%S")
date_start = datetime.strptime("2013-03-01 00:00:00", "%Y-%m-%d %H:%M:%S")
date_end = datetime.strptime("2015-01-01 00:00:00", "%Y-%m-%d %H:%M:%S")
date_mark = [datetime.strptime(x, "%Y-%m-%d %H:%M:%S") for x in
             ["2013-03-01 00:00:00",
              "2013-09-01 00:00:00",
              "2014-03-01 00:00:00",
              "2014-09-01 00:00:00",
              ]]

# load spatial data
river_pattern = joblib.load(results_dir+"river_pattern.joblib")
river_tracer = np.array([river_pattern[x].flatten(order="F")
                         for x in np.arange(len(river_pattern)-1)])
times = river_pattern["times"][:len(river_tracer)]
ntime = len(times)
ncell = len(river_tracer[0])

# read model  dimension
f = h5.File(multi_dir+"Coordinates.h5", "r")
x = np.array(f['Coordinates']["X [m]"])+593000
y = np.array(f['Coordinates']["Y [m]"])+114500
z = np.array(f['Coordinates']["Z [m]"])
ox = x[0]
oy = y[0]
oz = z[0]
ex = x[-1]
ey = y[-1]
ez = z[-1]
dx = np.diff(x)
dy = np.diff(y)
dz = np.diff(z)
nx = len(dx)
ny = len(dy)
nz = len(dz)
x = ox+np.cumsum(dx)-0.5*dx
y = oy+np.cumsum(dy)-0.5*dy
z = oz+np.cumsum(dz)-0.5*dz
x_array = np.array([x]*ny).flatten(order="C")
y_array = np.array([y]*nx).flatten(order="F")
f.close()

# read mass1 data
mass1 = np.genfromtxt(river_stage)
mass1_time = mass1[:, 0]/3600
mass1_date = batch_delta_to_time(date_origin,
                                 mass1_time,
                                 "%Y-%m-%d %H:%M:%S",
                                 "hours")
mass1_date = np.array([datetime.strptime(x, "%Y-%m-%d %H:%M:%S")
                       for x in mass1_date])
mass1_level = mass1[:, -1]

# read groundwater sample
sample_data = np.genfromtxt(sampling_file,
                            skip_header=1,
                            delimiter=",",
                            dtype="str")
sample_well = np.array(sample_data[:, 0])
sample_time = np.array([datetime.strptime(x, "%d-%b-%Y %H:%M:%S")
                        for x in sample_data[:, 1]])
sample_spc = np.array(sample_data[:, 3], dtype=np.float)
sample_u = np.array(sample_data[:, 3], dtype=np.float)
sample_data = dict()
for iwell in np.unique(sample_well):
    sample_data[iwell] = dict()
    sample_data[iwell]["time"] = sample_time[sample_well == iwell]
    sample_data[iwell]["spc"] = sample_spc[sample_well == iwell]
    sample_data[iwell]["u"] = sample_u[sample_well == iwell]
max_spc = np.max(sample_spc)
min_spc = np.min(sample_spc)
max_u = np.max(sample_u)
min_u = np.min(sample_u)
# representive wells
cluster_well = [[],
                ["399-2-32"],
                ["399-1-21A"],
                ["399-2-2"],
                ["399-1-10A"]]
cluster_well_name = [[],
                     ["Well 2-32"],
                     ["Well 1-21A"],
                     ["Well 2-02"],
                     ["Well 1-10A"]]

cluster_well = [[],
                ["399-1-21A"],
                ["399-2-32", "399-2-2"],
                [],
                [],
                ["399-1-10A"]]
cluster_well_name = [[],
                     ["Well 1-21A"],
                     ["Well 2-32", "Well 2-02"],
                     [],
                     [],
                     ["Well 1-10A"]]

# read material information
f = h5.File(material_file, "r")
material = f["Materials"]["Material Ids"][:].reshape((nx, ny, nz), order="F")
f.close()
river_bed = np.array([[z[np.where(material[ix, iy, :] > 0)[0][-1]]
                       for iy in range(ny)] for ix in range(nx)])
ringold_top = np.array([[z[np.where((material[ix, iy, :] != 1) *
                                    (material[ix, iy, :] != 0))[0][-1]]
                         for iy in range(ny)] for ix in range(nx)])
river_bed = river_bed  # .flatten(order="F")

# limit data to date range
time_x = batch_delta_to_time(date_origin,
                             times,
                             "%Y-%m-%d %H:%M:%S",
                             "hours")
time_date = np.array([datetime.strptime(x, "%Y-%m-%d %H:%M:%S")
                      for x in time_x])
time_index = (time_date >= date_start)*(time_date <= date_end)
times = np.array(times)[time_index]
time_x = time_x[time_index]
time_date = time_date[time_index]
river_tracer = river_tracer[time_index]
mean_tracer = np.mean(river_tracer, 0)
river_tracer_ori = copy.deepcopy(river_tracer)

# interp well data to simualtion data
interp_well = dict()
for iwell in well:
    interp_well[iwell] = np.interp(
        [(x-date_origin).total_seconds() for x in time_date],
        [(x-date_origin).total_seconds() for x in sample_data[iwell]["time"]],
        1-(sample_data[iwell]["spc"]-min_spc)/(max_spc-min_spc),
        left=np.nan,
        right=np.nan,)
cluster_data = np.transpose(
    np.array([interp_well[x] for x in well]))
dist_obs = np.ma.corrcoef(
    np.ma.masked_invalid(np.transpose(cluster_data)))
dist_obs = squareform(1-dist_obs, force="tovector", checks=True)
obs_cluster = hierarchy.ward(dist_obs)

# exclude low concentration locations from clustering
high_conc = np.where(np.max(river_tracer, 0) > 0.1)[0]
# exclude the second chanel part
river_index = copy.deepcopy(river_bed)
for iy in range(ny):
    left_bank = np.where(river_bed[:, iy] < 104)[0][0]
    river_index[left_bank:, iy] = 0
high_conc = high_conc[river_index.flatten(order="F")[high_conc] > 104]
river_tracer = river_tracer[:, high_conc]

# find poly_gon boundary
poly_bc = []
for iy in y:
    ix = x_array[high_conc][y_array[high_conc] == iy]
    poly_bc += [[min(ix), iy], [max(ix), iy]]
poly_bc = np.array(poly_bc)
# high_shape = MultiPoint(poly_bc).convex_hull
# poly_patch = PolygonPatch(high_shape)
# ax.add_patch(poly_patch)
# poly_bc = np.array(poly_bc)


# random.seed(1)
# ri = random.choices(np.arange(len(river_tracer[0])), k=100)
# imgfile = img_dir+"river_series.png"
# fig = plt.figure()
# ax = plt.subplot(111)
# for i in ri:
#     ax.plot(times, river_tracer[:, i])
# fig.set_size_inches(16, 4)
# plt.tight_layout()
# fig.savefig(imgfile, bbox_inches=0, density=600)
# plt.close(fig)

# # only take sample from spatial simualtions
random.seed(100)
spatial_index = random.choices(np.arange(len(river_tracer[0])), k=40000)
# compute distance for spatial clusters
dist_spatial = distance.pdist(X=np.transpose(
    river_tracer[:, spatial_index]), metric='correlation')
#    river_tracer[:, spatial_index]), metric='euclidean')
spatial_cluster = hierarchy.ward(dist_spatial)

# compute distance for spatial clusters
# random.seed(100)
# spatial_index_1 = random.choices(np.arange(len(river_tracer[0])), k=10000)
# cluster_data = river_tracer[:, spatial_index_1]
# dist_combined_1 = distance.pdist(X=np.transpose(
#     cluster_data), metric='correlation')
# cluster_data = np.hstack((river_tracer[:, spatial_index_1], np.transpose(
#     np.array([interp_well[x] for x in well]))))
# # dist_combined_1 = np.ma.corrcoef(
# #     np.ma.masked_invalid(np.transpose(cluster_data)))
# combined_cluster_1 = hierarchy.ward(dist_combined_1)
# # compute distance for spatial clusters

# compute distance for temporal clusters
dis_time = distance.pdist(X=river_tracer, metric='correlation')
time_cluster = hierarchy.ward(dis_time)


def plot_spatial_distance():
    """
    plot cluster dendrogram map
    """
    ncluster = 6
    custom_colors = [colors[0].tolist()]*2+colors[0:6].tolist()
#    custom_colors = custom_colors[::-1]
    hierarchy.set_link_color_palette([matplotlib.colors.rgb2hex(x)
                                      for x in custom_colors])
    imgfile = img_dir+"spatial_distance.png"
    fig = plt.figure()
    ax = plt.subplot(111)
    dendrogram(spatial_cluster,
               #               truncate_mode='lastp',
               truncate_mode='level',
               count_sort=False,
               p=5,
               #               show_contracted=True,
               color_threshold=12,
               #               link_color_func=lambda k: colors[k],
               above_threshold_color='black',
               orientation="top")
#    ax.set_ylim(0, 60)
    legend_elements = list()
    for icluster in range(ncluster):
        legend_elements.append(
            Line2D([0], [0],
                   color=colors[icluster],
                   label='Spatial cluster #'+str(icluster+1)))
    legend1 = plt.legend(handles=legend_elements,
                         fontsize=8,

                         loc='upper left')
    ax.add_artist(legend1)
    # fig.text(0, 0.5, "Height", rotation="90")
    # fig.text(0.45, 0.01, "Number of cells")
    # plt.subplots_adjust(top=0.98, bottom=0.18, left=0.08,
    #                     right=0.98, wspace=0.0, hspace=0.05)

    ax.set_xlabel("Number of cells")
    ax.set_ylabel("Height")
    fig.set_size_inches(8, 2.5)
    plt.subplots_adjust(top=0.98, bottom=0.25, left=0.065,
                        right=0.98, wspace=0.0, hspace=0.05)

#    plt.tight_layout()
    fig.savefig(imgfile, bbox_inches=0, dpi=300)
    plt.close(fig)


def plot_obs_distance():
    """
    plot cluster dendrogram map
    """
    colors = np.array(["#0000ff", "#ff0000", "#00ff00"])
    markers = np.array(["s", "^", "o"])
    hierarchy.set_link_color_palette([matplotlib.colors.rgb2hex(x)
                                      for x in colors])

    imgfile = img_dir+"obs_distance.png"
    fig = plt.figure()
    ax = plt.subplot(111)
    dendrogram(obs_cluster,
               truncate_mode='level',
               count_sort=False,
               p=5,
               above_threshold_color='black',
               color_threshold=0.8,
               labels=well_name,
               leaf_rotation=0,
               leaf_font_size=8,
               orientation="top")
    legend_elements = list()
    for icluster in range(ncluster):
        legend_elements.append(
            Line2D([0], [0],
                   color=colors[icluster],
                   label='Well cluster #'+str(icluster+1)))
    legend1 = plt.legend(handles=legend_elements,
                         fontsize=8,
                         loc='upper right')
    ax.add_artist(legend1)
    # fig.text(0.015, 0.5, "Height", rotation="90")
    # fig.text(0.45, 0.01, "Wells")
    ax.set_xlabel("Wells")
    ax.set_ylabel("Height")
    fig.set_size_inches(8, 2.5)
    plt.subplots_adjust(top=0.98, bottom=0.25, left=0.065,
                        right=0.98, wspace=0.0, hspace=0.05)
    fig.savefig(imgfile, bbox_inches=0, dpi=300)
    plt.close(fig)

# ncluster = 5
# cluster_ids = hierarchy.fcluster(
#     spatial_cluster, t=ncluster, criterion="maxclust")
# imgfile = img_dir+"subset_point.png"
# fig = plt.figure()
# ax = plt.subplot(111)
# ax.scatter(x_array[high_conc[ri]],
#            y_array[high_conc[ri]],
#            color=colors[cluster_ids-1])
# ax.set_xlim(ox, ex)
# ax.set_ylim(oy, ey)
# ax.set_xlabel("Easting (m)")
# ax.set_ylabel("Northing (m)")
# ax.set_aspect(1)
# legend_elements = list()
# for icluster in range(ncluster):
#     legend_elements.append(Patch(facecolor=colors[icluster],
#                                  edgecolor=colors[icluster],
#                                  label='Cluster '+str(icluster+1)))
# ax.legend(handles=legend_elements, loc='upper left')
# fig.set_size_inches(6, 9)
# plt.tight_layout()
# fig.savefig(imgfile, bbox_inches=0, density=600)
# plt.close(fig)


def plot_cluster_time_series():
    """
    plot time series of spatial clusters
    """

    ncluster = 5
    cluster_ids = hierarchy.fcluster(
        spatial_cluster, t=ncluster, criterion="maxclust")
    imgfile = img_dir+"cluster_time_series.png"
    fig, axs = plt.subplots(5, 1)
    for i, ax in enumerate(fig.axes[0:ncluster]):
        for i_line in np.where(cluster_ids == i+1)[0][0:500]:
            ax.plot(time_date,
                    river_tracer[:, spatial_index[i_line]],
                    color=colors[i],
                    linewidth=0.1,
                    )
        for iwell in cluster_well[i]:
            ax.plot(time_date,
                    interp_well[iwell],
                    color="black")
            # ax.plot(sample_data[iwell]["time"],
            #         1-(sample_data[iwell]["spc"]-min_spc)/(max_spc-min_spc),
            #         color="black")
        ax.set_ylim(0, 1)
        ax.set_xticks([])
#        ax.set_yticks([])
        ax.set_title("Spatial cluster #"+str(i+1), fontsize=12)
        if len(cluster_well[i]) > 0:
            legend_elements = list()
            for iwell_index, iwell in enumerate(cluster_well[i]):
                legend_elements.append(
                    Line2D([0], [0], color='black',
                           label=cluster_well_name[i][iwell_index]))
            ax.legend(handles=legend_elements,
                      loc="upper right",
                      fancybox=True, framealpha=0.4)
        ax.set_ylabel("River water fraction(-)")
#    ax.set_ylabel("Normalized concentration (-)")
    ax.set_ylabel("River water fraction(-)")
    ax.set_xticks(date_mark)  # rotation=45)
    ax.set_xticklabels([x.strftime("%Y-%m") for x in date_mark],
                       rotation=0)
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1])

    for i, ax in enumerate(fig.axes[ncluster:]):
        ax.set_visible(False)
    plt.subplots_adjust(top=0.97, bottom=0.04, left=0.15,
                        right=0.98, wspace=0.12, hspace=0.2)
    fig.set_size_inches(5, 9)
    fig.savefig(imgfile, bbox_inches=0, dpi=300, transparent=True)
    plt.close(fig)


def plot_spatial_clusters(ncluster):
    """
    plot ringole top
    """
#    ncluster = 7
    cluster_ids = hierarchy.fcluster(
        spatial_cluster, t=ncluster, criterion="maxclust")

    imgfile = img_dir+"spatial_cluster_"+str(ncluster)+".png"
    fig = plt.figure()
    ax = fig.add_subplot(111)
    cf = ax.contourf(x, y,
                     np.transpose(ringold_top),
                     cmap=plt.cm.Greys_r)
    ax.scatter(poly_bc[:, 0],
               poly_bc[:, 1],
               color="black",
               s=0.5)
    ax.scatter(x_array[high_conc[spatial_index]],
               y_array[high_conc[spatial_index]],
               linewidth=0,
               color=colors[cluster_ids-1], s=0.8)
    # for i, icluster in enumerate(well_cluster):
    #     ax.scatter(well_coord[i, 0],
    #                well_coord[i, 1],
    #                color="black",
    #                marker=well_marker[icluster-1],
    #                s=50)
    ax.set_xlim(ox, ex)
    ax.set_ylim(oy, ey)
    ax.set_xlabel("Easting (m)")
    ax.set_ylabel("Northing (m)")
    ax.set_aspect(1)
    cb = plt.colorbar(cf,
                      panchor=(0.1, 0.1),
                      fraction=0.03, pad=0.04)
    cb.set_label('Ringold elevation (m)', labelpad=13, rotation=270)

    legend_elements = list()
    for icluster in range(ncluster):
        legend_elements.append(Patch(facecolor=colors[icluster],
                                     edgecolor=colors[icluster],
                                     label='Spatial cluster #'+str(icluster+1)))
    legend1 = plt.legend(handles=legend_elements,
                         loc='upper left')

    legend_elements = list()
    for i, icluster in enumerate(well_marker[0:(ncluster+1)]):
        legend_elements.append(
            Line2D([0], [0], marker=icluster, color='w',
                   label='Well cluster #'+str(i+1),
                   markerfacecolor='black', markersize=10))
    # legend2 = plt.legend(handles=legend_elements,
    #                      loc='lower left')

    ax.add_artist(legend1)
#    ax.add_artist(legend2)

    fig.set_size_inches(6.7, 9)
    plt.subplots_adjust(top=0.98, bottom=0.05, left=0.1,
                        right=0.9, wspace=0.0, hspace=0.05)
    fig.savefig(imgfile, bbox_inches=0,
                dpi=300, transparent=True)
    plt.close(fig)


def plot_spatial_clusters_with_obs():
    """
    plot ringole top
    """
    ncluster = 3
    cluster_ids = hierarchy.fcluster(
        spatial_cluster, t=ncluster, criterion="maxclust")

    # get subclusters
    n_sub_cluster = 8
    sub_cluster_ids = hierarchy.fcluster(
        spatial_cluster, t=n_sub_cluster, criterion="maxclust")
    sub_cluster_ids = sub_cluster_ids-2
    sub_cluster_ids[cluster_ids == 1] = 1
    # sub_cluster_ids = -sub_cluster_ids+n_sub_cluster+1
    sub_clusters = np.unique(sub_cluster_ids)

    # generate sub boundaries
    # sub_bc = []
    # id_array = np.full(nx*ny, np.nan)
    # id_array[high_conc[spatial_index]] = sub_cluster_ids
    # id_array = id_array.reshape((nx, ny), order="F")
    # for ix in range(nx):
    #     iys = np.where(~np.isnan(id_array[ix, :]))[0]
    #     if len(iys) > 1:
    #         bc_iys = iys[np.where(np.diff(subset_bc[ix, iys]) != 0)]
    #         sub_bc += [[x[ix], y[iy]] for iy in bc_iys]
    # for iy in range(ny):
    #     ixs = np.where(~np.isnan(id_array[:, iy]))[0]
    #     if len(ixs) > 1:
    #         bc_ixs = ixs[np.where(np.diff(subset_bc[ixs, iy]) != 0)]
    #         sub_bc += [[x[ix], y[iy]] for ix in bc_ixs]

    # sub_bc = np.array(sub_bc)

    # all_xy = np.transpose(np.vstack((x_array[high_conc[spatial_index]],
    #                                  y_array[high_conc[spatial_index]])))
    # bc_xy = np.transpose(np.vstack((bc_point_x, bc_point_y)))
    # bc_xy = np.array([x for x in set(tuple(x) for x in bc_xy)
    #                   & set(tuple(x) for x in all_xy)])

    nobs_cluster = 3
    obs_ids = hierarchy.fcluster(
        obs_cluster, t=nobs_cluster, criterion="maxclust")

    imgfile = img_dir+"spatial_cluster_with_obs.png"
    fig = plt.figure()
    ax = fig.add_subplot(111)
    cf = ax.contourf(x, y,
                     np.transpose(ringold_top),
                     cmap=plt.cm.Greys_r)
    ax.scatter(poly_bc[:, 0],
               poly_bc[:, 1],
               color="black",
               s=0.5)
    # ax.scatter(sub_bc[:, 0],
    #            sub_bc[:, 1],
    #            color="black",
    #            zorder=1000,
    #            s=0.5)
    # for sub_index, sub_wells in enumerate(cluster_well):
    #     for name_index, iwell in enumerate(sub_wells):
    #         well_index = np.where(np.array(well) == iwell)[0][0]
    #         ax.text(well_coord[well_index][0]-50,
    #                 well_coord[well_index][1]+20,
    #                 cluster_well_name[
    #                     sub_index][name_index].split("Well ")[-1],
    #                 fontsize=10)

    ax.scatter(x_array[high_conc[spatial_index]],
               y_array[high_conc[spatial_index]],
               linewidth=0,
               color=colors[sub_cluster_ids-1], s=0.8)
    for i, icluster in enumerate(obs_ids):
        ax.scatter(well_coord[i, 0],
                   well_coord[i, 1],
                   color="black",
                   marker=well_marker[icluster-1],
                   s=40)
    ax.set_xlim(ox, ex)
    ax.set_ylim(oy, ey)
    ax.set_xlabel("Easting (m)")
    ax.set_ylabel("Northing (m)")
    ax.set_aspect(1)
    cb = plt.colorbar(cf,
                      panchor=(0.1, 0.1),
                      fraction=0.03, pad=0.04)
    cb.set_label('Ringold elevation (m)', labelpad=13, rotation=270)

    legend_elements = list()
    for icluster in range(len(sub_clusters)):
        legend_elements.append(Patch(facecolor=colors[icluster],
                                     edgecolor=colors[icluster],
                                     label='Spatial cluster #'+str(icluster+1)))
    legend1 = plt.legend(handles=legend_elements,
                         loc='upper left')

    legend_elements = list()
    for i, icluster in enumerate(well_marker[0:nobs_cluster]):
        legend_elements.append(
            Line2D([0], [0], marker=icluster, color='w',
                   label='Well cluster #'+str(i+1),
                   markerfacecolor='black', markersize=10))
    legend2 = plt.legend(handles=legend_elements,
                         loc='lower left')

    ax.add_artist(legend1)
    ax.add_artist(legend2)

    fig.set_size_inches(6.7, 9)
    plt.subplots_adjust(top=0.98, bottom=0.05, left=0.1,
                        right=0.9, wspace=0.0, hspace=0.05)
    fig.savefig(imgfile, bbox_inches=0,
                dpi=300, transparent=True)
    plt.close(fig)


def plot_sub_cluster_time_series():
    """
    plot time series of spatial clusters
    """
    ncluster = 3
    cluster_ids = hierarchy.fcluster(
        spatial_cluster, t=ncluster, criterion="maxclust")

    # get subclusters
    n_sub_cluster = 8
    sub_cluster_ids = hierarchy.fcluster(
        spatial_cluster, t=n_sub_cluster, criterion="maxclust")
    sub_cluster_ids = sub_cluster_ids-2
    sub_cluster_ids[cluster_ids == 1] = 1
    # sub_cluster_ids = -sub_cluster_ids+n_sub_cluster+1
    sub_clusters = np.unique(sub_cluster_ids)

    ncluster = len(sub_clusters)

    imgfile = img_dir+"cluster_time_series.png"
    fig, axs = plt.subplots(ncluster, 1)
    for i, ax in enumerate(fig.axes[0:ncluster]):
        for i_line in np.where(sub_cluster_ids == i+1)[0][0:500]:
            ax.plot(time_date,
                    river_tracer[:, spatial_index[i_line]],
                    color=colors[i],
                    linewidth=0.1,
                    )
        for well_index, iwell in enumerate(cluster_well[i]):
            # ax.plot(time_date,
            #         interp_well[iwell],
            #         color="black")
            ax.plot(sample_data[iwell]["time"],
                    (sample_data[iwell]["spc"]-min_spc)/(max_spc-min_spc),
                    linestyle=line_type[well_index],
                    color="black")
        ax.set_ylim(0, 1)
        ax.set_xticks([])
#        ax.set_yticks([])
        ax.set_title("Spatial cluster #"+str(i+1), fontsize=11)
        if len(cluster_well[i]) > 0:
            legend_elements = list()
            for iwell_index, iwell in enumerate(cluster_well[i]):
                legend_elements.append(
                    Line2D([0], [0], color='black',
                           linestyle=line_type[iwell_index],
                           label=cluster_well_name[i][iwell_index]))
            ax.legend(handles=legend_elements,
                      loc="lower right", fontsize=8,
                      fancybox=True, framealpha=0.4)
#        ax.set_ylabel("River water fraction(-)")
#    ax.set_ylabel("Normalized concentration (-)")
#    ax.set_ylabel("River water fraction(-)")
    ax.set_xticks(date_mark)  # rotation=45)
    ax.set_xticklabels([x.strftime("%Y-%m") for x in date_mark],
                       rotation=0)
    ax.set_yticks([0, 0.2, 0.4, 0.6, 0.8, 1])
    fig.text(0.02, 0.41, "Normalized concentration(-)", rotation=90)
    for i, ax in enumerate(fig.axes[ncluster:]):
        ax.set_visible(False)
    plt.subplots_adjust(top=0.97, bottom=0.04, left=0.15,
                        right=0.98, wspace=0.12, hspace=0.22)
    fig.set_size_inches(4.5, 9)
    fig.savefig(imgfile, bbox_inches=0, dpi=300, transparent=True)
    plt.close(fig)


# poly_bc = np.array(poly_bc)
# imgfile = img_dir+"ringold_top_2.png"
# fig = plt.figure()
# ax = plt.subplot(111)
# cf = ax.contourf(x, y,
#                  np.transpose(ringold_top),
#                  cmap=plt.cm.jet)
# ax.scatter(poly_bc[:, 0],
#            poly_bc[:, 1],
#            color="black",
#            s=0.5)
# ax.set_xlim(ox, ex)
# ax.set_ylim(oy, ey)
# ax.set_xlabel("Easting (m)")
# ax.set_ylabel("Northing (m)")
# ax.set_aspect(1)
# fig.set_size_inches(6, 9)
# plt.tight_layout()
# fig.savefig(imgfile, bbox_inches=0, dpi=300)
# plt.close(fig)


# poly_bc = np.array(poly_bc)
# imgfile = img_dir+"mean_tracer.png"
# fig = plt.figure()
# ax = plt.subplot(111)
# ax.contourf(x, y,
#             np.transpose(mean_tracer.reshape((nx, ny), order="F")),
#             cmap=plt.cm.jet)
# ax.scatter(poly_bc[:, 0],
#            poly_bc[:, 1],
#            color="black",
#            s=0.5)
# ax.set_xlim(ox, ex)
# ax.set_ylim(oy, ey)
# ax.set_aspect(1)
# fig.set_size_inches(6, 9)
# plt.tight_layout()
# fig.savefig(imgfile, bbox_inches=0, density=600)
# plt.close(fig)


def plot_temporal_distance():
    """
    plot cluster dendrogram map
    """

    colors = np.array(["#0000ff", "#ff0000", "#00ff00"])
    markers = np.array(["s", "^", "o"])
    ncluster = 3

    hierarchy.set_link_color_palette([matplotlib.colors.rgb2hex(x)
                                      for x in colors])

    imgfile = img_dir+"distance_time.png"
    fig = plt.figure()
    ax = plt.subplot(111)
    dendrogram(time_cluster,
               truncate_mode='lastp',
               count_sort=False,
               p=30,
               above_threshold_color='black',
               color_threshold=0.66,
               show_leaf_counts=True,
               #               labels=well_name,
               leaf_rotation=90,
               leaf_font_size=8,
               orientation="top")
    legend_elements = list()
    for icluster in range(ncluster):
        legend_elements.append(
            Line2D([0], [0],
                   color=colors[icluster],
                   label='Temporal cluster #'+str(icluster+1)))
    legend1 = plt.legend(handles=legend_elements,
                         fontsize=8,
                         loc='upper right')
    ax.add_artist(legend1)
    # fig.text(0.0, 0.47, "Height", rotation="90")
    # fig.text(0.4, 0.01, "Number of snapshots")
    ax.set_xlabel("Number of snapshots")
    ax.set_ylabel("Height")
    # plt.subplots_adjust(top=1, bottom=0.1, left=0.08,
    #                     right=0.98, wspace=0.0, hspace=0.05)
#    plt.tight_layout()
    fig.set_size_inches(8, 2.5)
    plt.subplots_adjust(top=0.98, bottom=0.25, left=0.065,
                        right=0.98, wspace=0.0, hspace=0.05)
    fig.savefig(imgfile, bbox_inches=0, dpi=300)
    plt.close(fig)


def plot_spatial_map_of_temporal_cluster():
    """
    plot spatial map of time clusters
    """

    ncluster = 3
    time_ids = hierarchy.fcluster(
        time_cluster, t=ncluster, criterion="maxclust")
    mass1_index = [np.argmin(abs(mass1_time-x)) for x in times]
    river_tracer_plot = np.full((len(times), nx*ny), np.nan)
    river_tracer_plot[:, high_conc] = river_tracer
#    colors = cm.Set1(np.linspace(0, 1, 9))
    imgfile = img_dir+"temporal_cluster_map.png"
    fig, axs = plt.subplots(1, ncluster+1,
                            gridspec_kw={"width_ratios": [1, 1, 1, 0.05]})
    for i, ax in enumerate(fig.axes[0:ncluster]):
        time_subset = (time_ids == i+1)
        mean_subset = (np.mean(river_tracer_plot[time_subset, ], axis=0)).reshape(
            (nx, ny), order="F")
        cf = ax.contourf(x, y,
                         np.transpose(mean_subset),
                         levels=np.arange(0, 1.01, 0.01),
                         cmap=plt.cm.jet)
        ax.scatter(poly_bc[:, 0],
                   poly_bc[:, 1],
                   color="black",
                   s=0.5)
        ax.set_xlabel("Easting (m)")
        ax.set_ylabel("Northing (m)")
        ax.set_title("Temporal cluster #"+str(i+1), fontsize=12)
        ax.set_xlim(ox, ex)
        ax.set_ylim(oy, ey)
        ax.set_aspect(1)
    ax = fig.axes[-1]
    ax.axis("off")
    cb = fig.colorbar(cf, ax=ax, fraction=1, pad=0,
                      ticks=np.arange(0, 1.01, 0.1))
    cb.set_label('Mean normalized concentration (-)',
                 labelpad=13, rotation=270)
    fig.set_size_inches(11.5, 5)
    plt.tight_layout()
    fig.savefig(imgfile, bbox_inches=0, dpi=300, transparent=True)
    plt.close(fig)


def plot_temporal_clusters():
    colors = np.array(["#0000ff", "#ff0000", "#00ff00"])
    markers = np.array(["s", "^", "o"])
    # colors = cm.Set1(np.linspace(0, 1, 9))[
    #     np.array([6, 2, 1, 3, 4, 0, 7, 8, 5])]
    # markers = ["o",  "^", "X", "*", "s", "v", "<", ">"]

    ncluster = 3
    time_ids = hierarchy.fcluster(
        time_cluster, t=ncluster, criterion="maxclust")

#    colors = cm.Set1(np.linspace(0, 1, 9))

    imgfile = img_dir+"time_cluster.png"
    fig = plt.figure()
    ax = plt.subplot(111)
    for cluster_id in np.unique(time_ids):
        mass1_index = [np.argmin(abs(mass1_time-x))
                       for x in times[time_ids == cluster_id]]
        ax.scatter(mass1_date[mass1_index],
                   mass1_level[mass1_index],
                   color=colors[cluster_id-1],
                   s=20,
                   marker=markers[cluster_id-1],
                   zorder=1000)
    ax.plot(mass1_date,  # [mass1_index],
            mass1_level,  # [mass1_index],
            linewidth=0.5,
            color="black",
            zorder=100)
    mass1_index = [np.argmin(abs(mass1_time-x))
                   for x in times]
    ax.set_xlim(mass1_date[mass1_index[0]],
                mass1_date[mass1_index[-1]])
    ax.set_ylim(104, 108)
    ax.set_ylabel("River stage (m)")
    legend_elements = list()
    for icluster in range(ncluster):
        legend_elements.append(
            Line2D([0], [0], marker=markers[icluster], color="w",
                   label='Temporal cluster #'+str(icluster+1),
                   markerfacecolor=colors[icluster], markersize=6))

    ax.legend(handles=legend_elements, loc='upper right')
    fig.set_size_inches(11, 2.5)
    plt.tight_layout()
    fig.savefig(imgfile, bbox_inches=0, dpi=300, transparent=True)
    plt.close(fig)


def plot_combined_clusters():
    """
    plot ringole top
    """
    ncluster = 5
    cluster_ids = hierarchy.fcluster(
        combined_cluster, t=ncluster, criterion="maxclust")

    imgfile = img_dir+"combined_cluster.png"
    fig = plt.figure()
    ax = fig.add_subplot(111)
    cf = ax.contourf(x, y,
                     np.transpose(ringold_top),
                     cmap=plt.cm.Greys_r)
    ax.scatter(poly_bc[:, 0],
               poly_bc[:, 1],
               color="black",
               s=0.5)
    ax.scatter(x_array[high_conc[spatial_index]],
               y_array[high_conc[spatial_index]],
               color=colors[cluster_ids[0:-len(well)]-1], s=0.1)
    for iwell in range(len(well)):
        ax.scatter(well_coord[iwell, 0],
                   well_coord[iwell, 1],
                   color="black",
                   marker=well_marker[cluster_ids[-len(well):][iwell]-1],
                   s=50)
    ax.set_xlim(ox, ex)
    ax.set_ylim(oy, ey)
    ax.set_xlabel("Easting (m)")
    ax.set_ylabel("Northing (m)")
    ax.set_aspect(1)
    cb = plt.colorbar(cf,
                      panchor=(0.1, 0.1),
                      fraction=0.03, pad=0.04)
    cb.set_label('Ringold elevation (m)', labelpad=13, rotation=270)

    legend_elements = list()
    for icluster in range(ncluster):
        legend_elements.append(Patch(facecolor=colors[icluster],
                                     edgecolor=colors[icluster],
                                     label='Spatial cluster #'+str(icluster+1)))
    legend1 = plt.legend(handles=legend_elements,
                         loc='upper left')

    legend_elements = list()
    for i, icluster in enumerate(well_marker):
        legend_elements.append(
            Line2D([0], [0], marker=icluster, color='w',
                   label='Well cluster #'+str(i+1),
                   markerfacecolor='black', markersize=10))
    legend2 = plt.legend(handles=legend_elements,
                         loc='lower left')

    ax.add_artist(legend1)
    ax.add_artist(legend2)

    fig.set_size_inches(6.7, 9)
    plt.subplots_adjust(top=0.98, bottom=0.05, left=0.1,
                        right=0.9, wspace=0.0, hspace=0.05)
    fig.savefig(imgfile, bbox_inches=0,
                dpi=300, transparent=True)
    plt.close(fig)


# plot_spatial_distance()
# plot_cluster_time_series()
# plot_spatial_clusters()
# plot_temporal_distance()
# plot_spatial_cluster()
# plot_temporal_clusters()

# plot_spatial_clusters()


# combined_cluster = combined_cluster_2
# spatial_index = spatial_index_1
# plot_combined_clusters()

# for ncluster in np.arange(1, 10):
#     plot_spatial_clusters(ncluster)
