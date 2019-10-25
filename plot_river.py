import pickle as pk
import numpy as np
import matplotlib.pyplot as plt
from sklearn.externals import joblib
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


img_dir = "/media/sf_e/john/optim_5/figures/"

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

river_file = "/media/sf_e/john/optim_5/" + \
    "DatumH_River_filtered_6h_322.txt"
river_level = np.loadtxt(river_file)
river_level = np.asarray(river_level)
river_level[:, 0] = river_level[:, 0] / 3600 / 24+(365+366)

# fig_name = img_dir + "river_level.png"
# fig = plt.figure()
# ax1 = fig.add_subplot(111)
# ax1.plot((river_level[:, 0] * 24),
#          (river_level[:, 3]),
#          color="black", linewidth=0.5)
# #ax1.set_title("River stage")
# ax1.set_xlabel('Time (year)')
# ax1.set_ylabel('River stage (m)')
# #ax1.set_xlim(8784, river_level[-1, 0]*24)
# ax1.set_xticks(year_label_loc)
# ax1.set_xticklabels(year_label)
# fig.set_size_inches(10, 3)
# fig.tight_layout()
# fig.savefig(fig_name, dpi=600, transparent=False)
# plt.close(fig)


date_origin = datetime.strptime("2008-01-01 00:00:00", "%Y-%m-%d %H:%M:%S")
year_label = [
    "2012-10",
    #    "2013-01",
    "2013-04",
    #    "2013-07",
    "2013-10",
    #    "2014-01",
    "2014-04",
    "2014-10",
    "2015-04"
]

year_label = [
    "2010-01",
    #    "2013-01",
    "2011-01",
    #    "2013-07",
    "2012-01",
    #    "2014-01",
    "2013-01",
    "2014-01",
    "2015-01"
]


year_label_loc = batch_time_to_delta(date_origin, year_label,
                                     "%Y-%m") / 3600

fill_start = 43245.73774269
all_end = (fill_start+np.arange(100)*100).tolist()
end_label = batch_delta_to_time(date_origin, all_end, "%Y-%m-%d", "hours")
end_label = ["Time: "+x for x in end_label]
batch_delta_to_time(date_origin, [100], "%Y-%m-%d", "hours")
level_index = np.where(((river_level[:, 0]*24) > year_label_loc[0])
                       * ((river_level[:, 0]*24) < year_label_loc[-1]))[0]
for fill_end, ititle in zip(all_end[-2:], end_label[-2:]):
    fig_name = img_dir + str(fill_end)+".png"
    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    ax1.plot((river_level[level_index, 0] * 24),
             (river_level[level_index, 3]),
             zorder=100,
             color="black", linewidth=0.5)
    ax1.fill_between([fill_start, fill_end],
                     [104]*2,
                     [109]*2,
                     zorder=1,
                     color="lightgrey", linewidth=0.5)
#    ax1.set_title(ititle)
    ax1.set_xlabel('Date')
    ax1.set_ylabel('River stage (m)')
#    ax1.set_xlim(year_label_loc[0], year_label_loc[-1])
    ax1.set_ylim(104, 109)
    ax1.set_xticks(year_label_loc[0:-1])
    ax1.set_xticklabels(year_label[0:-1])
    fig.set_size_inches(7, 2.5)
    fig.tight_layout()
    fig.savefig(fig_name, dpi=300, transparent=False)
    plt.close(fig)
