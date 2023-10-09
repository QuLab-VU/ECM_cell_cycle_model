import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import networkx as nx
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

    # While df is not empty
    while len(df) > 0:
        # Get the first row
        row = df.iloc[0]
        # Remove the first row from df
        df = df.iloc[1:]

        cell_id = row['ID']
        lineage = row['Lineage']
        event = row['event']
        birth_time = row['Birth time']
        division_history = row['divisions history']
        lineage_history = row['division lineage history']
        end_of_recording = row['end of recording']

        life_breakdown = []

        if len(division_history) == 0:
            # Handle the case when division_history is empty
            if lineage == 'Initial cell':
                starting_cell = True
                start_time = 0
            else:
                starting_cell = False
                start_time = birth_time

            life_breakdown.append({
                'Old ID': cell_id,
                'New ID': cell_id,
                'Lineage': lineage,
                'Event': event,
                'Start Time': start_time,
                'Birth Time': birth_time,
                'End Time': end_of_recording,
                'Starting Cell': starting_cell
            })

        else: # Handle the case when division_history is not empty
            # From birth to first mitosis
            new_cell_id = cell_id
            new_ID_list.append(cell_id)
            end_time = division_history[0]
            event_new = 'Mitosis'
            birth_time_new = birth_time
            if birth_time_new == 'Initial cell':
                start_time = 0
                starting_cell = True
            else:
                start_time = birth_time_new
                starting_cell = False
            new_lineage = lineage

            life_breakdown.append({
                'Old ID': cell_id,
                'New ID': new_cell_id,
                'Lineage': new_lineage,
                'Event': event_new,
                'Start Time': start_time,
                'Birth Time': birth_time_new,
                'End Time': end_time,
                'Starting Cell': starting_cell
            })
            start_time = end_time
            new_lineage = new_cell_id

            # First occurence of lineage_history[i] in df["ID"]
            df.loc[df["ID"] == lineage_history[0], "Lineage"] = new_cell_id


            for i in range(len(division_history)):
                # From birth to first division
                new_cell_id = ID_list[-1] + 1
                ID_list.append(new_cell_id)
                new_ID_list.append(new_cell_id)

                if i == len(division_history) - 1: # From last mitosis to end of recording
                    end_time = end_of_recording
                    event_new = event
                    birth_time_new = division_history[i-1]

                else: # From mitosis to mitosis
                    end_time = division_history[i+1]
                    event_new = 'Mitosis'
                    birth_time_new = division_history[i-1]

                life_breakdown.append({
                    'Old ID': cell_id,
                    'New ID': new_cell_id,
                    'Lineage': new_lineage,
                    'Event': event_new,
                    'Start Time': start_time,
                    'Birth Time': start_time,
                    'End Time': end_time,
                    'Starting Cell': False
                })

                start_time = end_time
                new_lineage = new_cell_id

                # First occurence of lineage_history[i] in df["ID"]
                if i < len(division_history) - 1:
                    df.loc[df["ID"] == lineage_history[i+1], "Lineage"] = new_cell_id



        cell_life_breakdown.extend(life_breakdown)

    population_history = pd.DataFrame(cell_life_breakdown)
    population_history["Lifetime"] = population_history["End Time"] - population_history["Start Time"]


    return population_history


# Source: https://stackoverflow.com/questions/29586520/can-one-get-hierarchical-graphs-from-networkx-with-python-3


def plot_digraph_population_history(population_history):
    # Create a directed graph
    G = nx.DiGraph()

    # Define colors for different fates
    fate_colors = {
        'Mitosis': 'lightblue',
        'Apoptosis': 'lightcoral',
        'End of experiment': 'lightgreen'
    }

    # Define shapes for starting and non-starting cells
    shape_starting = 's'  # Square for starting cells
    shape_non_starting = 'o'  # Circle for non-starting cells

    # Add nodes to the graph
    for _, row in population_history.iterrows():
        node_id = row['New ID']
        node_label = f"Cell {node_id}"
        node_color = fate_colors.get(row['Event'], 'lightgray')  # Default to lightgray for unknown events
        node_size = row['Lifetime'] * 10  # Adjust the multiplier as needed
        node_shape = shape_starting if row['Starting Cell'] else shape_non_starting
        G.add_node(node_id, label=node_label, color=node_color, size=node_size, shape=node_shape)

    # Add edges to represent lineage relationships
    for _, row in population_history.iterrows():
        lineage = row['Lineage']
        if lineage != 'Initial cell':
            G.add_edge(lineage, row['New ID'])

    # Create a layout for the graph
    # pos = nx.planar_layout(G)
    # pos = nx.fruchterman_reingold_layout(G)
    pos = nx.kamada_kawai_layout(G)
    # pos = nx.shell_layout(G)
    # pos = nx.spectral_layout(G)

    # Extract node colors, sizes, and shapes for plotting
    node_colors = [G.nodes[node]['color'] for node in G.nodes]
    node_sizes = [G.nodes[node]['size'] for node in G.nodes]
    node_shapes = [G.nodes[node]['shape'] for node in G.nodes]

    # Plot the lineage graph with custom shapes and sizes
    plt.figure(figsize=(12, 8), dpi=400)
    labels = nx.get_node_attributes(G, 'label')
    for i, node in enumerate(G.nodes):
        nx.draw_networkx_labels(G, pos, labels={node: labels[node]}, font_size=8, font_color='black', font_weight='bold', verticalalignment='bottom')
        nx.draw_networkx_nodes(G, pos, nodelist=[node], node_color=node_colors[i], node_size=node_sizes[i], node_shape=node_shapes[i])
    nx.draw_networkx_edges(G, pos)
    plt.title('Cell Lineage Graph')
    plt.show()



if __name__=="__main__":
    # Sample data
    data = {
        "ID": [0, 1, 2, 3, 4, 5],
        "Lineage": ['Initial cell', 'Initial cell', 0, 1, 1, 3],
        "event": ['End of experiment', 'End of experiment', 'End of experiment', 'End of experiment', 'Apoptosis', 'End of experiment'],
        "division": [1, 2, 1, 2, 1, 2],
        "Birth time": ['Initial cell', 'Initial cell', 12.84280936455, 0.1003344481605, 36.05504587156, 18.4615384615405],
        "divisions history": [[12.84280936455], [0.1003344481605, 36.05504587156], [], [18.4615384615405], [], []],
        "division lineage history": [[2], [3, 4], [], [5], [], []],
        "end of recording": [40.09174311927, 40.09174311927, 40.09494744386, 40.0832057012475, 39.072821100917, 40.0944368131905]
    }

    df = pd.DataFrame(data)

    # Break down cell life
    population_history = break_down_cell_life(df)
    print(population_history)

    plot_digraph_population_history(population_history)