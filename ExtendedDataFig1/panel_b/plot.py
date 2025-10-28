import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter


def simulation_process(selection_intensity):
    results_M50 = np.zeros(5)
    results_M100_positive_1 = np.zeros(5)
    results_M100_negative_1 = np.zeros(5)
    results_M100_0 = np.zeros(5)

    ite = 50
    time = 100000000

    for index in range(1,ite+1):
        temp_M50 = np.loadtxt(f'data_EDF1/data_figure2_panel_b_M50_selectionIntenstiy{selection_intensity}_{time}times_{index}.txt')
        temp_M100_positive_1 = np.loadtxt(f'data_EDF1/data_figure2_panel_b_M100_gamma1.0_selectionIntenstiy{selection_intensity}_{time}times_{index}.txt')
        temp_M100_negative_1 = np.loadtxt(f'data_EDF1/data_figure2_panel_b_M100_gamma-1.0_selectionIntenstiy{selection_intensity}_{time}times_{index}.txt')
        temp_M100_0 = np.loadtxt(f'data_EDF1/data_figure2_panel_b_M100_gamma0.0_selectionIntenstiy{selection_intensity}_{time}times_{index}.txt') 
        #print(index)
        results_M100_0 += temp_M100_0[:,1]
        results_M50 += temp_M50[:,1]
        results_M100_positive_1 += temp_M100_positive_1[:,1]
        results_M100_negative_1 += temp_M100_negative_1[:,1]

    results_M100_0 /= (time*ite*100)
    results_M50 /= (time*ite*50)
    results_M100_positive_1 /= (time*ite*100)
    results_M100_negative_1 /= (time*ite*100)

    b_M50 = temp_M50[:,0]
    b_M100_positive_1 = temp_M100_positive_1[:,0]
    b_M100_negative_1 = temp_M100_negative_1[:,0]
    b_M100_0 = temp_M100_0[:,0]

    return b_M50, b_M100_0, b_M100_negative_1, b_M100_positive_1, results_M50, results_M100_positive_1, results_M100_negative_1, results_M100_0



def panel_b_cooperation_defference(selection_intensity):
    plt.rcParams['font.family'] = 'Arial'
    plt.rc('font', size=11)

    colors = [ '#458c6e', '#c1272c', '#e69542', '#0070bb']
    #            green       red      orange      blue
    
    b_M50, b_M100_0, b_M100_negative_1, b_M100_positive_1, results_M50, results_M100_positive_1, results_M100_negative_1, results_M100_0 = simulation_process(selection_intensity)
    b_M50_theoretical = np.array([value * 0.001 + 5.48 for value in range(841)])
    b_M100_negative_1_theoretical = np.array([value * 0.001 + 6.16 for value in range(841)])
    b_M100_positive_1_theoretical = np.array([value * 0.001 + 2.67 for value in range(861)])
    b_M100_0_theoretical = np.array([value * 0.001 + 4.0 for value in range(801)])

    k1_M50 = 124.78607381671657
    k2_M50 = 21.024593151714843

    k1_M100_0 = 109.44588795165677
    k2_M100_0 = 24.605459908226607

    k1_M100_negative_1 = 153.90399849262238
    k2_M100_negative_1 = 23.45769829843473

    k1_M100_positive_1 = 89.73796920376826
    k2_M100_positive_1 = 28.5816942413923

    y_M50_theoretical = selection_intensity /2 * (-k1_M50 + k2_M50 * b_M50_theoretical)
    y_M100_negative_1_theoretical = selection_intensity /2 * (-k1_M100_negative_1 + k2_M100_negative_1 * b_M100_negative_1_theoretical)
    y_M100_positive_1_theoretical = selection_intensity /2 * (-k1_M100_positive_1 + k2_M100_positive_1 * b_M100_positive_1_theoretical)
    y_M100_0_theoretical = selection_intensity /2 * (-k1_M100_0 + k2_M100_0 * b_M100_0_theoretical)

    fig, ax = plt.subplots(figsize=(4, 4))

    critical_b_values = [5.935243, 3.139701,  6.560916, 4.448033,]

    for i,b in enumerate(critical_b_values):
        ax.axvline(x=b, color = colors[i], linestyle='--', linewidth = 1.5)
        #ax.text(b, ax.get_ylim()[1], rf'$(b/c)^*=${b:.2f}', ha='center', va='top')
        #ax.text(b, ax.get_ylim()[1], f'{b:.2f}', ha='center', va='top')

    ax.axhline(y=0.0, color='gray', linestyle='--', linewidth = 1.5)


    ax.plot(b_M50_theoretical, y_M50_theoretical, color=colors[0], linewidth = 1.5, label=r'$\overline{M}=1$')
    ax.plot(b_M100_positive_1_theoretical, y_M100_positive_1_theoretical, color=colors[1], linewidth = 1.5, label=r'$\gamma = 1.0$')
    ax.plot(b_M100_negative_1_theoretical, y_M100_negative_1_theoretical, color=colors[2], linewidth = 1.5, label=r'$\gamma = -1.0$')
    ax.plot(b_M100_0_theoretical, y_M100_0_theoretical, color=colors[3], linewidth = 1.5, label=r'$\gamma = 0.0$')

    ax.scatter(b_M50, 2*results_M50-1, facecolor=colors[0], marker='s', s=50, edgecolor=colors[0], linewidths=1.5)
    ax.scatter(b_M100_positive_1, 2*results_M100_positive_1-1, facecolor=colors[1], marker='s', s=50, edgecolor=colors[1], linewidths=1.5)
    ax.scatter(b_M100_negative_1, 2*results_M100_negative_1-1, facecolor=colors[2], marker='s', s=50, edgecolor=colors[2], linewidths=1.5)
    ax.scatter(b_M100_0, 2*results_M100_0-1, facecolor=colors[3], marker='s', s=50, edgecolor=colors[3], linewidths=1.5)


    ax.set_xlim(2.5,7.2)
    
    ax.set_ylim(-3.5e-2,3.5e-2)
    ax.set_xlabel(r'Benefit, $b$', fontsize=11.5)
    ax.set_ylabel(r'Cooperation minus defection rate, $\langle x_C - x_D \rangle$', fontsize=11.5)

    ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax.ticklabel_format(style='scientific', axis='y', scilimits=(0, 0))

    plt.legend(loc='best')
    plt.title(r'$\overline{M}=2$')
    #plt.savefig(f'ExtendedDataFig1_panelb_delta{selection_intensity}.svg')
    plt.show()




panel_b_cooperation_defference(0.005)