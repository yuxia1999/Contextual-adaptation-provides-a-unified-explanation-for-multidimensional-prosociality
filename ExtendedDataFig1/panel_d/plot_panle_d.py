import numpy as np
import matplotlib.pylab as plt 
import matplotlib.patches as patches

def plot_gamma_bc(M, degree, index):
    filename = f'DonationGame_SF_DB_bcCritical_50_degree{degree}_{index}_M{M}_withInfiniteGamma.txt'
    data_matrix = np.loadtxt(filename)
    u_set = [0.001, 0.01, 0.1]
    colors = ['#c1272c', '#458c6e', '#0070bb'] 
    plt.figure(figsize=(4.2,4.2))

    for i, u in enumerate(u_set):
        rows = data_matrix[np.isclose(data_matrix[:, 1], u)]
        gamma_values = rows[:, 0]
        y_values = rows[:, 4]
        plt.plot(gamma_values, y_values, linestyle='-', marker='o', markersize = 4,linewidth = 2.5, color=colors[i], label=f'u = {u}')

    plt.xlabel(r'Adaptation scaling exponent, $\gamma$')
    plt.ylabel(r'Critical benefit-to-cost threshold, $(b/c)^{\ast}$')
    plt.title(r'$\overline{M}=2$')
    plt.legend(loc = 'lower left')

    plt.xlim(-12.8,12.8)
    plt.xticks([-12.8,-10, -5, 0, 5, 10, 12.8],[r'$-\infty$',r'$-$10', r'$-$5', r'0', r'5', r'10', r'$+\infty$'])
    plt.ylim([2.7, 8.3])
    # plt.savefig(f'panel_d_M{M}_degree{degree}_index{index}_line.svg')

    plt.show()

plot_gamma_bc(100, 4, 1)
