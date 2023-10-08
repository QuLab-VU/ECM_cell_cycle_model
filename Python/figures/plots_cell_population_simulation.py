import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
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


def break_down_cell_life(df): # Doesn't work yet
    cell_life_breakdown = []
    old_id_to_new_id = {}  # Maps old Cell ID to new Cell ID
    new_id_counter = 0

    ID_list = df['ID'].tolist()
    new_ID_list = []

    for _, row in df.iterrows():
        cell_id = row['ID']
        lineage = row['Lineage']
        events = row['event']
        divisions = row['division']
        birth_time = row['Birth time']
        division_history = row['divisions history']
        lineage_history = row['division lineage history']
        end_of_recording = row['end of recording']

        life_breakdown = []
        start_time = birth_time  # Initialize start_time here

        if len(division_history) == 0:
            # Handle the case when division_history is empty
            for i in range(len(events)):
                if i == len(events) - 1:
                    end_time = end_of_recording
                else:
                    end_time = end_of_recording if i == len(events) - 1 else birth_time

                life_breakdown.append({
                    'Old ID': cell_id,
                    'New ID': new_cell_id,
                    'Lineage': new_cell_id,
                    'Event': events[i],
                    'Start Time': start_time,
                    'End Time': end_time
                })
                start_time = end_time  # Update start_time

        for i in range(len(division_history)):
            if cell_id in new_ID_list:
                new_cell_id = ID_list[-1] + 1
                ID_list.append(new_cell_id)
                new_ID_list.append(new_cell_id)
            else:
                new_cell_id = cell_id
                new_ID_list.append(cell_id)

            if isinstance(division_history[i], list) and division_history[i]:  # Check if it's a non-empty list
                end_time = division_history[i][0]
                life_breakdown.append({
                    'Old ID': cell_id,
                    'New ID': new_cell_id,
                    'Lineage': new_cell_id if i == 0 else lineage,
                    'Event': 'Mitosis',
                    'Start Time': start_time,
                    'End Time': end_time
                })

                # Update Lineage for related cells
                for cell_lineage, cell_id_list in zip(lineage_history[i], division_history[i]):
                    if isinstance(cell_id_list, list):
                        for old_id in cell_id_list:
                            old_id_to_new_id[old_id] = new_cell_id
                start_time = end_time  # Update start_time

            if i < len(events):
                if i == len(events) - 1:
                    end_time = end_of_recording
                elif isinstance(division_history[i], list) and division_history[i]:
                    end_time = division_history[i][0]
                else:
                    end_time = end_of_recording

                life_breakdown.append({
                    'Old ID': cell_id,
                    'New ID': new_cell_id,
                    'Lineage': new_cell_id,
                    'Event': events[i],
                    'Start Time': start_time,
                    'End Time': end_time
                })
                start_time = end_time  # Update start_time

        cell_life_breakdown.extend(life_breakdown)

    return pd.DataFrame(cell_life_breakdown)


if __name__=="__main__":
    # Sample data
    data = {
        "ID": [0, 1, 2, 3, 4, 5, 6, 7],
        "Lineage": ['Initial cell', 0, 0, 0, 1, 2, 1, 4],
        "event": ['End of experiment', 'End of experiment', 'End of experiment', 'End of experiment',
                  'End of experiment', 'End of experiment', 'End of experiment', 'End of experiment'],
        "division": [3, 3, 2, 1, 2, 1, 1, 1],
        "Birth time": ['Initial cell', 0, 21.67224080268, 38.138755980861, 20.86956521739, 42.96499748019,
                       43.966507177029996, 40.133779264210006],
        "divisions history": [
            [0, 21.67224080268, 38.138755980861], [20.86956521739, 43.966507177029996], [42.96499748019], [],
            [40.133779264210006], [], [], []],
        "division lineage history": [[1, 2, 3], [4, 6], [5], [], [7], [], [], []],
        "end of recording": [50.0956937799, 50.0956937799, 50.09606221652, 50.088573258490996, 50.066889632110005,
                             50.084619966098, 50.082587892655994, 50.086018816616004]
    }

    df = pd.DataFrame(data)

    # Break down cell life
    cell_life_breakdown = break_down_cell_life(df)
    print(cell_life_breakdown)