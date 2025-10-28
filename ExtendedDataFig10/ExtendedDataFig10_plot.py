import matplotlib.pyplot as plt
import numpy as np

def extendedDataFig_variance_plot(num_nodes, overlineM, degree):
    plt.rcParams['font.family'] = 'Arial'
    plt.rc('font', size=12.5)
    fig,ax = plt.subplots(figsize=(4, 4))

    colors = ['#c1272c', '#0070bb', '#458c6e', '#e69542']

    for index in range(1,101):

        file_name = f"EDF10_data/SF_{num_nodes}_overlineM{overlineM}_degree{degree}_{index}_variance.txt"
        data = np.loadtxt(file_name, delimiter=' ')
        gamma = data[:, 0]
        original_variance = data[:, 1]
        weighted_variance = data[:, 2]

        ax.scatter(gamma, weighted_variance, color=colors[int(overlineM-1)], alpha=0.05)

    plt.axvline(x=1, linewidth = 1.5, color='grey', linestyle='--')
    plt.xlabel(r'$\gamma$')
    plt.ylabel('Variance of $d_i/L_i$')

    plt.title(rf'$\overline{{M}}$ = {overlineM:.2f}')

    plt.xlim(-4,4)

    plt.show()

def extendedDataFig_variance_plot_combined(num_nodes, degree):
    plt.rcParams['font.family'] = 'Arial'
    plt.rc('font', size=12.5)
    fig, ax = plt.subplots(figsize=(4, 4))

    colors = ['#c1272c', '#0070bb', '#458c6e', '#e69542']
    overlineM_list = [1, 2, 3, 4]

    for overlineM in overlineM_list:
        all_gamma = None
        all_weighted_variances = []

        for index in range(1, 101):
            file_name = f"EDF10_data/SF_{num_nodes}_overlineM{overlineM}_degree{degree}_{index}_variance.txt"
            data = np.loadtxt(file_name, delimiter=' ')
            gamma = data[:, 0]
            weighted_variance = data[:, 2]

            if all_gamma is None:
                all_gamma = gamma  

            all_weighted_variances.append(weighted_variance)

        all_weighted_variances = np.array(all_weighted_variances)  

        mean_weighted_variance = np.mean(all_weighted_variances, axis=0)
        std_weighted_variance = np.std(all_weighted_variances, axis=0)

        ax.plot(all_gamma, mean_weighted_variance, color=colors[overlineM-1], label=rf'$\overline{{M}}$ = {overlineM}')

        ax.fill_between(all_gamma,
                        mean_weighted_variance - std_weighted_variance,
                        mean_weighted_variance + std_weighted_variance,
                        color=colors[overlineM-1],
                        alpha=0.2)

    plt.axvline(x=1, linewidth=1.5, color='grey', linestyle='--')
    plt.xlabel(r'$\gamma$')
    plt.ylabel('Variance of $d_i/L_i$')
    plt.title(rf'$\overline{{d}}$ = {degree}')
    plt.xlim(-4, 4)

    plt.legend(frameon=False)

    # plt.savefig(f"variance_SF{num_nodes}_deg{degree}.svg")
    plt.show()
for deg in [4,6,8,10]:
    extendedDataFig_variance_plot_combined(50, deg)