import pickle as pk
import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
import matplotlib.gridspec as gridspec


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
    return(y)


img_dir = "/mnt/e/john/vec/stage/"

date_origin = datetime.strptime("2010-01-01 00:00:00", "%Y-%m-%d %H:%M:%S")
year_label = ["2009",
              "2010",
              "2011",
              "2012",
              "2013",
              "2014",
              "2015"]
year_label_loc = batch_time_to_delta(date_origin, year_label,
                                     "%Y") / 3600

river_file = "/mnt/e/john/optim_5/" + \
    "DatumH_River_filtered_6h_322.txt"
river_level = np.loadtxt(river_file)
river_level = np.asarray(river_level)
river_level[:, 0] = river_level[:, 0] / 3600 / 24+(365+366)


date_origin = datetime.strptime("2008-01-01 00:00:00", "%Y-%m-%d %H:%M:%S")
year_label = [
    "2014-04-25",
    "2014-05-01",
    "2014-05-10",
    "2014-05-20",
    "2014-05-31",
    "2014-06-05",
]

year_label_loc = batch_time_to_delta(date_origin, year_label,
                                     "%Y-%m-%d") / 3600

fill_start = np.sum([366, 365*3, 366, 365])*24+np.sum([31, 28, 31, 30]*24)
all_end = range(fill_start, fill_start+31*24, 1)

end_label = batch_delta_to_time(
    date_origin, all_end, "%Y-%m-%d %H:%M:%S", "hours")
end_label = ["Time: "+x for x in end_label]
batch_delta_to_time(date_origin, [100], "%Y-%m-%d", "hours")
# level_index = np.where(((river_level[:, 0]*24) > year_label_loc[0])
#                        * ((river_level[:, 0]*24) < year_label_loc[-1]))[0]
for ifig,(fill_end, ititle) in enumerate(zip(all_end, end_label)):
    fig_name = img_dir + str(ifig)+".png"
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.plot((river_level[:, 0] * 24),
             (river_level[:, 3]),
             zorder=100,
             color="black", linewidth=1)
    ax1.plot([fill_end]*2, [104, 109],color="red",linewidth=2)
    ax1.set_title(ititle)
    ax1.set_xlabel('')
    ax1.set_ylabel('River stage (m)')
    ax1.set_xlim(year_label_loc[0], year_label_loc[-1])
    ax1.set_ylim(105.5, 107.5)
    ax1.set_xticks(year_label_loc[1:-1])
    ax1.set_xticklabels(year_label[1:-1])
    fig.set_size_inches(7, 2.5)
    fig.tight_layout()
    fig.savefig(fig_name, dpi=160, transparent=False)
    plt.close(fig)
