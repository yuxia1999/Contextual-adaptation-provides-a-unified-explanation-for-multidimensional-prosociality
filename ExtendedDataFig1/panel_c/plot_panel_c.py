import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from matplotlib.ticker import ScalarFormatter


def plot_mutationRate_bc(M, degree, index):
    fig, ax = plt.subplots(figsize=(4.3,4.3))

    filename = f'DonationGame_SF_DB_bcCritical_50_degree{degree}_{index}_M50.txt'
    data_matrix = np.loadtxt(filename)
    gamma_set = [1.0]
    colors = ['#458c6e']  
    for i, gamma in enumerate(gamma_set):
        rows = data_matrix[np.isclose(data_matrix[:, 0], gamma)]
        u_values = rows[:, 1]
        bc_values = rows[:, 4]

        positive_mask = bc_values >= 0
        negative_mask = bc_values < 0

        plt.scatter(u_values[positive_mask], bc_values[positive_mask], marker='o', s =20, color=colors[i], label=r'$\overline{M} = 1$')
        
        plt.scatter(u_values[negative_mask], bc_values[negative_mask], marker='o', s =20, color=colors[i])

        if np.any(positive_mask) and np.any(negative_mask):
            last_positive_x = u_values[positive_mask][-1]
            first_negative_x = u_values[negative_mask][0]

            transition_x = (last_positive_x + first_negative_x) / 2

            plt.axvline(x=transition_x, linestyle='--', color=colors[i])

    filename = f'DonationGame_SF_DB_bcCritical_50_degree{degree}_{index}_M{M}.txt'
    data_matrix = np.loadtxt(filename)
    gamma_set = [1.0, 0.0]
    colors = ['#c1272c', '#0070bb']  
    for i, gamma in enumerate(gamma_set):
        rows = data_matrix[np.isclose(data_matrix[:, 0], gamma)]
        u_values = rows[:, 1]
        bc_values = rows[:, 4]

        positive_mask = bc_values >= 0
        negative_mask = bc_values < 0

        plt.scatter(u_values[positive_mask], bc_values[positive_mask], marker='o', s =20, color=colors[i], label=fr'$\gamma = {gamma}$')
        
        plt.scatter(u_values[negative_mask], bc_values[negative_mask], marker='o', s =20, color=colors[i])

        if np.any(positive_mask) and np.any(negative_mask):
            last_positive_x = u_values[positive_mask][-1]
            first_negative_x = u_values[negative_mask][0]

            transition_x = (last_positive_x + first_negative_x) / 2

            plt.axvline(x=transition_x, linestyle='--', color=colors[i])
    
    plt.axhline(y=0.0, linestyle='--', color='grey')
    plt.ticklabel_format(axis = 'y',style='sci',  scilimits=(0,0))

    plt.xlabel(r'Mutation rate, $\mu$')
    plt.ylabel(r'Critical benefit-to-cost threshold, $(b/c)^{\ast}$')
    plt.title(rf'$\bar{{d}} = {degree}, \overline{{M}}=2$')
    plt.legend(loc='upper left')
    plt.xlim(1e-3,1)
    plt.ylim(-5.3e3,5.3e3)
    plt.ylim(-2.13e3,2.13e3)
    plt.yticks([-2e3,-1e3,0,1e3,2e3])
    
    plt.xscale('log')
    ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax.ticklabel_format(style='scientific', axis='y', scilimits=(0, 0))

    # plt.savefig(f'ExtendedDataFig1_panel_c_u_bc_degree{degree}_{index}_M{M}.svg')

    plt.show()



plot_mutationRate_bc(100, 30, 1)
