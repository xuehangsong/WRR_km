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
import matplotlib.gridspec as gridspec


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
    # dis_river = [DTWDistance(river_tracer[:, ri[i]], river_tracer[:, ri[j]])
    #              for i in range(len(ri)) for j in range(len(ri))]
    ri = random.choices(np.arange(len(river_tracer[0])), k=400)
    dis_river_tri = [[fastdtw(river_tracer[:, ri[i]], river_tracer[:, ri[j]])[0]
                      for i in range(j+1)] for j in range(len(ri))]
    dis_river = np.zeros((len(ri), len(ri)))
    for i, j in enumerate(dis_river_tri):
        dis_river[i][0:len(j)] = j
    dis_river = dis_river+np.transpose(dis_river)
    dis_river = distance.squareform(dis_river, force="tovector", checks=True)


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


multi_dir = "/media/sf_e/john/optim_5/single_h5/"
results_dir = "/media/sf_e/john/optim_5/"
img_dir = "/media/sf_e/john/optim_5/figures/"
material_file = "/media/sf_e/john/optim_5/300A_material.h5"
river_stage = "/media/sf_e/john/optim_5/DatumH_River_filtered_6h_321.txt"

we1ll_name = ["399-1-10A",
              "399-1-21A",
              "399-2-01",
              "399-2-02",
              "399-2-23",
              "399-2-03",
              "399-2-32",
              "399-2-33",
              "399-2-05",
              "399-2-07",
              "399-3-19",
              "399-3-09",
              "399-4-09"]
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
well_marker = ["o", "*", "^", "X"]


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

# load data
river_pattern = joblib.load(results_dir+"river_pattern.joblib")
river_tracer = np.array([river_pattern[x].flatten(order="F")
                         for x in np.arange(len(river_pattern)-1)])
times = river_pattern["times"][:len(river_tracer)]
ntime = len(times)
ncell = len(river_tracer[0])

# define dimension
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

# read mass1
mass1 = np.genfromtxt(river_stage)
mass1_time = mass1[:, 0]/3600
mass1_date = batch_delta_to_time(date_origin,
                                 mass1_time,
                                 "%Y-%m-%d %H:%M:%S",
                                 "hours")
mass1_date = np.array([datetime.strptime(x, "%Y-%m-%d %H:%M:%S")
                       for x in mass1_date])
mass1_level = mass1[:, -1]

# read material
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

# exclude low concentration locations from clustering
# low_conc = np.where(np.max(river_tracer, 0) < 0.1)[0]
high_conc = np.where(np.max(river_tracer, 0) > 0.1)[0]


# exclude part below river
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

random.seed(1)
ri = random.choices(np.arange(len(river_tracer[0])), k=100)
imgfile = img_dir+"river_series.png"
fig = plt.figure()
ax = plt.subplot(111)
for i in ri:
    ax.plot(times, river_tracer[:, i])
fig.set_size_inches(16, 4)
plt.tight_layout()
fig.savefig(imgfile, bbox_inches=0, density=600)
plt.close(fig)

random.seed(1)
ri = random.choices(np.arange(len(river_tracer[0])), k=10000)

# direct compute distance
dis_river = distance.pdist(X=np.transpose(
    #    river_tracer[:, ri]), metric='euclidean')
    river_tracer[:, ri]), metric='correlation')


# def plot_distance_cluster():
#     """
#     plot cluster dendrogram map
#     """
subset_cluster = hierarchy.ward(dis_river)
imgfile = img_dir+"subset_distance_cluster.png"
fig = plt.figure()
ax = plt.subplot(111)
dendrogram(subset_cluster,
           orientation="right")
fig.set_size_inches(6, 15)
plt.tight_layout()
fig.savefig(imgfile, bbox_inches=0, density=600)
plt.close(fig)

ncluster = 5
cluster_ids = hierarchy.fcluster(
    subset_cluster, t=ncluster, criterion="maxclust")
colors = cm.Set1(np.linspace(0, 1, 9))


imgfile = img_dir+"subset_point.png"
fig = plt.figure()
ax = plt.subplot(111)
ax.scatter(x_array[high_conc[ri]],
           y_array[high_conc[ri]],
           color=colors[cluster_ids-1])
ax.set_xlim(ox, ex)
ax.set_ylim(oy, ey)
ax.set_xlabel("Easting (m)")
ax.set_ylabel("Northing (m)")
ax.set_aspect(1)
legend_elements = list()
for icluster in range(ncluster):
    legend_elements.append(Patch(facecolor=colors[icluster],
                                 edgecolor=colors[icluster],
                                 label='Cluster '+str(icluster+1)))
ax.legend(handles=legend_elements, loc='upper left')
fig.set_size_inches(6, 9)
plt.tight_layout()
fig.savefig(imgfile, bbox_inches=0, density=600)
plt.close(fig)


def test2():
    ncluster = 5
    imgfile = img_dir+"cluster_series.png"
    fig, axs = plt.subplots(5, 1)
    for i, ax in enumerate(fig.axes[0:ncluster]):
        for i_line in np.where(cluster_ids == i+1)[0][0:500]:
            ax.plot(time_date,
                    river_tracer[:, ri[i_line]],
                    color=colors[i],
                    linewidth=0.1,
                    )
        ax.set_ylim(0, 1)
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title("Spatial cluster #"+str(i+1), fontsize=12)
    ax.set_ylabel("Normalized concentration (-)")
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

# high_shape = MultiPoint(poly_bc).convex_hull
# poly_patch = PolygonPatch(high_shape)
# ax.add_patch(poly_patch)
# poly_bc = np.array(poly_bc)


def test():
    imgfile = img_dir+"ringold_top.png"
    # gs = gridspec.GridSpec(1, 2,
    #                        width_ratios=[10, 1])
    fig = plt.figure()
#    ax = fig.add_subplot(gs[0])
    ax = fig.add_subplot(111)
    cf = ax.contourf(x, y,
                     np.transpose(ringold_top),
                     cmap=plt.cm.Greys_r)
    ax.scatter(poly_bc[:, 0],
               poly_bc[:, 1],
               color="black",
               s=0.5)
    ax.scatter(x_array[high_conc[ri]],
               y_array[high_conc[ri]],
               color=colors[cluster_ids-1], s=0.3)
    for i, icluster in enumerate(well_cluster):
        ax.scatter(well_coord[i, 0],
                   well_coord[i, 1],
                   color="black",
                   marker=well_marker[icluster-1],
                   s=100)
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


poly_bc = np.array(poly_bc)
imgfile = img_dir+"ringold_top_2.png"
fig = plt.figure()
ax = plt.subplot(111)
cf = ax.contourf(x, y,
                 np.transpose(ringold_top),
                 cmap=plt.cm.jet)
ax.scatter(poly_bc[:, 0],
           poly_bc[:, 1],
           color="black",
           s=0.5)
ax.set_xlim(ox, ex)
ax.set_ylim(oy, ey)
ax.set_xlabel("Easting (m)")
ax.set_ylabel("Northing (m)")
ax.set_aspect(1)
fig.set_size_inches(6, 9)
plt.tight_layout()
fig.savefig(imgfile, bbox_inches=0, dpi=300)
plt.close(fig)


poly_bc = np.array(poly_bc)
imgfile = img_dir+"mean_tracer.png"
fig = plt.figure()
ax = plt.subplot(111)
ax.contourf(x, y,
            np.transpose(mean_tracer.reshape((nx, ny), order="F")),
            cmap=plt.cm.jet)
ax.scatter(poly_bc[:, 0],
           poly_bc[:, 1],
           color="black",
           s=0.5)
ax.set_xlim(ox, ex)
ax.set_ylim(oy, ey)
ax.set_aspect(1)
fig.set_size_inches(6, 9)
plt.tight_layout()
fig.savefig(imgfile, bbox_inches=0, density=600)
plt.close(fig)


dis_time = distance.pdist(X=river_tracer, metric='correlation')
time_cluster = hierarchy.ward(dis_time)
imgfile = img_dir+"distance_time.png"
fig = plt.figure()
ax = plt.subplot(111)
dendrogram(time_cluster,
           orientation="right")
fig.set_size_inches(6, 15)
plt.tight_layout()
fig.savefig(imgfile, bbox_inches=0, density=600)
plt.close(fig)


def test3():
    ncluster = 3
    time_ids = hierarchy.fcluster(
        time_cluster, t=ncluster, criterion="maxclust")
    mass1_index = [np.argmin(abs(mass1_time-x)) for x in times]
    river_tracer_plot = np.full((len(times), nx*ny), np.nan)
    river_tracer_plot[:, high_conc] = river_tracer
    colors = cm.Set1(np.linspace(0, 1, 9))
    imgfile = img_dir+"time_map.png"
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


def test4():
    ncluster = 3
    time_ids = hierarchy.fcluster(
        time_cluster, t=ncluster, criterion="maxclust")
    mass1_index = [np.argmin(abs(mass1_time-x)) for x in times]
    colors = cm.Set1(np.linspace(0, 1, 9))

    imgfile = img_dir+"time_point.png"
    fig = plt.figure()
    ax = plt.subplot(111)
    ax.scatter(mass1_date[mass1_index],
               mass1_level[mass1_index],
               color=colors[time_ids-1],
               s=50,
               zorder=1000)
    ax.plot(mass1_date,  # [mass1_index],
            mass1_level,  # [mass1_index],
            linewidth=0.5,
            color="black",
            zorder=100)
    ax.set_xlim(mass1_date[mass1_index[0]],
                mass1_date[mass1_index[-1]])
    ax.set_ylim(104, 108)
#    ax.set_xlabel("Easting (m)")
    ax.set_ylabel("River stage (m)")
    legend_elements = list()
    for icluster in range(ncluster):
        legend_elements.append(Patch(facecolor=colors[icluster],
                                     edgecolor=colors[icluster],
                                     label='Temporal cluster #'+str(icluster+1)))
    ax.legend(handles=legend_elements, loc='upper right')
    fig.set_size_inches(11, 2.5)
    plt.tight_layout()
    fig.savefig(imgfile, bbox_inches=0, dpi=300, transparent=True)
    plt.close(fig)
