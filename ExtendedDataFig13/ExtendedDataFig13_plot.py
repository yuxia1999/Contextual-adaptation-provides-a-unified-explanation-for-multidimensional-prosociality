import numpy as np
import matplotlib.pylab as plt 
import matplotlib.patches as patches


def plot_gamma_bc(overlineM, mutation_rate, index):
    plt.rcParams['font.family'] = 'Arial'
    plt.rc('font', size=12.5)

    if index == 0 and mutation_rate == 0.01:
        bc_initial = 52.2579
    elif index == 1 and mutation_rate == 0.01:
        bc_initial = 7.9707
    elif index == 2 and mutation_rate == 0.01:
        bc_initial = 465.2082
    elif index == 3 and mutation_rate == 0.01:
        bc_initial = 6.2735

    fig, ax = plt.subplots(figsize=(3.25, 3.25))

    data_name = ["FourthGradeStudent",  "KarateClub", "SouthernWomen","ModeMath"]
    title_name = ["Students",  "Karate Club", "Women", "Superintendents" ]
    num_nodes_set = [24, 34, 18, 30]
    filename = f"EDF13_data/{data_name[index]}_K1K2.txt"

    data_matrix = np.loadtxt(filename)

    u_set = [mutation_rate]

    colors = ['#c1272c', '#458c6e', '#0070bb']  
    plt.axhline(y=bc_initial, linewidth = 2.0, linestyle='-', color=colors[1], label = f'$\overline{{M}}=1.0$')

    for i, u in enumerate(u_set):
        rows = data_matrix[np.isclose(data_matrix[:, 2], u) & np.isclose(data_matrix[:, 0], overlineM[0])]

        gamma_values = rows[:, 1]
        y_values = rows[:, 5]

        plt.plot(gamma_values, y_values, linestyle='-', linewidth = 2, color=colors[0], label = f'$\overline{{M}}={overlineM[0]}$')
    
    plt.axvline(x=1, linewidth = 1.5, color='grey', linestyle='--')

    plt.xlabel(r'Adaptation scaling exponent, $\gamma$')
    plt.ylabel(r'Critical threshold, $(b/c)^{\ast}$')
    plt.title(f"{title_name[index]}", fontsize=14)
    plt.legend()

    plt.xlim(-10,10)

    # plt.savefig(f'ExtendedFig_EN_optimal_{data_name[index]}_mutationRate{mutation_rate}.svg')

    plt.show()

for index in range(4):
    plot_gamma_bc([2.0],0.01,index)


