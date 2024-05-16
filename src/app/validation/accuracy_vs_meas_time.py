import numpy as np
import matplotlib.pyplot as plt

from src.utils.data import get_fig_filename
from src.utils.run_for_all_data import run_for_all_data, load_results
from src.config.locations import LOCATIONS
from src.navigation.calculations import latlon_distance
from src.navigation.curve_fit_method import solve
from src.config.parameters import CurveFitMethodParameters

# first: 5, alt
# second: 5, no alt
# third: 5, alt, no eststate
RUN_NEW_DATA = False
TIME_STEP = 5  # min
DATA_IDX = -1
EXCLUDED_DATA = ["val07",
                 "val05",
                 # "val10",
                 "val04",
                 # "val08"
                 ]

home_lon, home_lat, home_alt = LOCATIONS["HOME"][0], LOCATIONS["HOME"][1], LOCATIONS["HOME"][2]
LON_HOME, LAT_HOME, ALT_HOME = LOCATIONS["HOME"][0], LOCATIONS["HOME"][1], LOCATIONS["HOME"][2]
plt.style.use("Plots/ctuthesis.mplstyle")


# section: ------------------------------------------------------------------------- SOLVING
def slice_data_into_time_chunks(data, chunk_time):
    if chunk_time is None:
        return [data]

    time_chunks = list()
    start_time, end_time = data[0, 0], data[-1, 0]
    i = 0

    while start_time + chunk_time * i < end_time:
        mask = (data[:, 0] >= start_time + chunk_time * i) & (data[:, 0] < start_time + chunk_time * (i + 1))
        time_chunks.append(data[mask])
        i += 1

    return time_chunks


def cb_solve(nav_data, satellites, default_parameters: CurveFitMethodParameters, *_):
    results = list()
    default_parameters.min_curve_length = 0
    # default_parameters.iteration.alt.init_step = 0
    est_state = None

    cummulative_data = None
    for data_chunk in slice_data_into_time_chunks(nav_data, TIME_STEP * 60):
        if cummulative_data is None:
            cummulative_data = data_chunk
        else:
            cummulative_data = np.vstack((cummulative_data, data_chunk))
        print(cummulative_data.shape)
        res = solve(cummulative_data, satellites, default_parameters, init_state=est_state)
        est_state = res
        results.append((*res, cummulative_data.shape[0]))
    return results


# section: ------------------------------------------------------------------------- ANALYSIS
results = load_results("accuracy_vs_meas_time", index=DATA_IDX)
if RUN_NEW_DATA or results is None:
    run_for_all_data(cb_solve, "accuracy_vs_meas_time")
    results = load_results("accuracy_vs_meas_time")

# section: data
for excl in EXCLUDED_DATA:
    if excl in results:
        print(f"Excluding {excl}")
        results.pop(excl)

parsed_results = dict()
for res in results.values():
    print(len(res), res)

parsed_results = list()
parsed_counts = list()
for exp_name, res_exp in results.items():
    ress = list()
    count = list()
    for i, res in enumerate(res_exp):
        lat, lon, alt, off, dft, cnt = res
        hor_dist = latlon_distance(LAT_HOME, lat, LON_HOME, lon)
        abs_dist = latlon_distance(LAT_HOME, lat, LON_HOME, lon, ALT_HOME, alt)
        #                    0    1    2    3    4    5         6
        ress.append(hor_dist)
        count.append(cnt)
    parsed_results.append(ress)
    parsed_counts.append(count)


# clip all results to minimum common length
def clip_to_min_len(results, counts, percentile):
    mean_dist_arr = np.empty((len(results)))
    for i, exp_res in enumerate(results):
        mean_dist_arr[i] = np.array([exp_res]).mean()
    mask = mean_dist_arr <= np.percentile(mean_dist_arr, percentile)
    mask_idx = np.where(mask == np.True_)[0]
    good_parsed_results = [results[i] for i in mask_idx]
    good_counts = [counts[i] for i in mask_idx]
    min_common_length = min([len(v) for v in good_parsed_results])
    parsed_results_comm, counts_comm = list(), list()
    for exp_res, cnts in zip(good_parsed_results, good_counts):
        parsed_results_comm.append([exp_res[i] for i in range(min_common_length)])
        counts_comm.append([cnts[i] for i in range(min_common_length)])
    return np.array(parsed_results_comm), np.array(counts_comm)


res_arr_all, cnt_arr_all = clip_to_min_len(parsed_results, parsed_counts, 100)
res_arr_cep, cnt_arr_cep = clip_to_min_len(parsed_results, parsed_counts, 50)
res_arr_95 , cnt_arr_95  = clip_to_min_len(parsed_results, parsed_counts, 95)


# section: plotting
# all results, just for debug
plt.figure()
for i, exp_res in enumerate(parsed_results):
    plt.plot(exp_res, label=list(results.keys())[i])
# plt.yscale('log')
plt.legend()

# all results, just for debug
plt.figure()
for i, exp_res in enumerate(parsed_counts):
    plt.plot(exp_res, label=list(results.keys())[i])
plt.legend()

for res_arr, cnt_arr in zip([res_arr_all], [cnt_arr_all]):
    fig, (ax1, axB) = plt.subplots(2, figsize=(5, 6))
    res_arr /= 1000
    max_arr, mean_arr, min_arr = res_arr.max(axis=0), res_arr.mean(axis=0), res_arr.min(axis=0)
    max_cnt_arr, mean_cnt_arr, min_cnt_arr = cnt_arr.max(axis=0), cnt_arr.mean(axis=0), cnt_arr.min(axis=0)
    times = np.array((list(range(1, len(mean_arr) + 1)))) * TIME_STEP
    eb = ax1.errorbar(times, mean_arr, yerr=[min_arr, max_arr], fmt=".-", label="Accuracy")
    eb[-1][0].set_linestyle('--')
    eb[-1][0].set_linewidth(0.8)
    ax1.set_ylabel("2D error (km)")

    ax2 = ax1.twinx()
    eb2 = ax2.plot(times, mean_cnt_arr, ".-", color="orange", label="Mean frame count")
    ax2.set_ylabel("Frame count")
    ax2.grid(False)
    ax1.set_xlabel("Time (min)")

    eb = axB.errorbar(times, mean_arr, yerr=[min_arr, max_arr], fmt=".-", label="__nolabel__")
    eb[-1][0].set_linestyle('--')
    eb[-1][0].set_linewidth(0.8)
    axB.set_xlim([14, ax1.get_xlim()[1]])
    axB.set_ylim([0, 8])
    # axB.set(xlim=[14, ax1.get_xlim()[1]], ylim=[0, 8], aspect=1)
    axB.set_ylabel("2D error (km)")
    axB.set_xlabel("Time (min) (close-up)")

    fig.legend()
    fig.tight_layout()

    plt.savefig(get_fig_filename("validation\\" + f"exp_accuracy_vs_meas_time", idx=False), dpi=600)

plt.show()
