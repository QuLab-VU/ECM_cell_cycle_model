import numpy as np
import matplotlib.pyplot as plt

def plot_cell_count(cell_count, run_time):
    t_sim = np.arange(0, run_time, 1 / 3)
    cell_count_t = np.zeros(len(t_sim))
    time_tmp_m = 0
    for time_tmp_p, cell_count_time_tmp in cell_count:
        cell_count_t[(time_tmp_m <= t_sim) & (t_sim < time_tmp_p)] = cell_count_time_tmp
        time_tmp_m = time_tmp_p
    DIP_sim = np.log2(cell_count_t / cell_count_t[0])
    fig, ax = plt.subplots()
    ax.scatter(t_sim, DIP_sim, color="black", s=1)
    ax.set_xlabel('Time (h)')
    ax.set_ylabel('DIP')
    plt.show()

    return fig, ax
