import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import matplotlib.patches as patches

def dictatorGame_pew_frequency_continual(lambda1, lambda2, lambda3, selection_intensity):
    p = 1/2 + selection_intensity * (- lambda2 * 1/6 - lambda3 * 1/12)
    return p

def plot_pq_vs_u_withSimulation(file_pattern, deg_values, M_numbers, gamma, selection_intensity=0.001):
    plt.rcParams['font.family'] = 'Arial'
    plt.rc('font', size=12.5)

    fig, ax = plt.subplots(figsize=(3.25, 3.25))
    colors = ['#c1272c', '#0070bb']  

    for i, M_number in enumerate(M_numbers):
        for deg in deg_values:
            file_name = file_pattern.format(deg=deg, M_number=M_number)
            data = np.loadtxt(file_name)
            gamma_values, u_values, lambda1_values, lambda2_values, lambda3_values = data.T

            filtered_indices = np.where(gamma_values == gamma)
            u_values_filtered = u_values[filtered_indices]
            lambda1_filtered = lambda1_values[filtered_indices]
            lambda2_filtered = lambda2_values[filtered_indices]
            lambda3_filtered = lambda3_values[filtered_indices]

            p_values = []

            for l1, l2, l3 in zip(lambda1_filtered, lambda2_filtered, lambda3_filtered):
                p = dictatorGame_pew_frequency_continual(l1, l2, l3, selection_intensity)
                p_values.append(p)

            plt.plot(u_values_filtered, p_values, linewidth=2.0, label=rf'$\overline{{M}}={M_number/30:.2f}$', color=colors[i], linestyle='-')
            

    ite = 30
    data30 = np.zeros(5)
    data50 = np.zeros(5)

    for index in range(1,ite+1):
        filename30 = f'EDF7_data/DTG_M30_gamma1.0_selectionIntenstiy0.005_100000000times_{index}.txt'
        filename50 = f'EDF7_data/DTG_M50_gamma1.0_selectionIntenstiy0.005_100000000times_{index}.txt'

        temp30 = np.loadtxt(filename30)
        temp50 = np.loadtxt(filename50)

        data30 += temp30
        data50 += temp50

    data30 /= ite * 100000000 * 30
    data50 /= ite * 100000000 * 50

    u_set = [10**(-2.5), 10**(-2), 10**(-1.5), 10**(-1), 10**(-0.5)]
    plt.scatter(u_set, data30, marker='s', color=colors[0])
    plt.scatter(u_set, data50, marker='s', color=colors[1])

    plt.xlim(1e-3,1)
    plt.ylim(0.462, 0.503)
    plt.xscale('log')
    plt.xlabel(r'Mutation rate, $\mu$', fontsize = 14)
    plt.ylabel('Mean value', fontsize = 14)
    plt.title(fr'$\bar{{d}} = {deg}, \gamma = {gamma}$', fontsize = 14)
    plt.legend(loc='upper right')

    ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax.ticklabel_format(style='scientific', axis='y', scilimits=(0, 0))

    # plt.savefig(rf"ExtendedDataFig7_DTG_pqFrequency_withSimulation_SF_degree{deg_values[0]}_new.svg")
    plt.show()

def plot_pq_vs_gamma_variousMutationRate(file_pattern, deg_values, M_numbers, mutation_rates,selection_intensity=0.005):
    plt.rcParams['font.family'] = 'Arial'
    plt.rc('font', size=12.5)

    fig, ax = plt.subplots(figsize=(3.25, 3.25))

    colors = ['#c1272c', '#0070bb'] 

    M_number = M_numbers[0]  

    for i, mutation_rate in enumerate(mutation_rates):
        for deg in deg_values:
            file_name = file_pattern.format(deg=deg, M_number=M_number)
            data = np.loadtxt(file_name)
            gamma_values, u_values, lambda1_values, lambda2_values, lambda3_values = data.T
            filtered_indices = np.where(u_values == mutation_rate)
            gamma_values_filtered = gamma_values[filtered_indices]
            lambda1_filtered = lambda1_values[filtered_indices]
            lambda2_filtered = lambda2_values[filtered_indices]
            lambda3_filtered = lambda3_values[filtered_indices]

            phi_values = []
            for l1, l2, l3 in zip(lambda1_filtered, lambda2_filtered, lambda3_filtered):
                phi = dictatorGame_pew_frequency_continual(l1, l2, l3, selection_intensity)
                phi_values.append(phi)

            plt.plot(gamma_values_filtered, phi_values, linewidth=2.0, label=rf'$\mu$={mutation_rate}, $\langle p \rangle$', color=colors[i], linestyle='-')

            if M_number == 30:
                mean_phi_value = np.mean(phi_values)
                plt.axhline(mean_phi_value, color=colors[i], linestyle='-', linewidth=2.0)

    plt.xlim(-4,4)
    ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax.ticklabel_format(style='scientific', axis='y', scilimits=(0, 0))

    plt.axvline(x=1, linewidth = 1.5, color='grey', linestyle='--')

    plt.xlabel(r'$\gamma$', fontsize = 14)
    plt.ylabel('Mean value', fontsize = 14)
    plt.title(fr'$\overline{{M}}={M_number}$', fontsize = 14)
    plt.legend(fontsize = 12)

    # plt.savefig(rf"ExtendedDataFig7_DTG_pqFrequency_vs_gamma_verious_u_SF_degree{deg_values}_M{M_number}_delta{selection_intensity}.svg")
    plt.show()

def plot_pq_vs_M_SF_RR_variousMutationRate_line(deg, M_value, gamma, selection_intensity=0.005):
    plt.rcParams['font.family'] = 'Arial'
    plt.rc('font', size=12.5)
    fig, ax = plt.subplots(figsize=(3.25, 3.25))
    colors = ['#c1272c', '#0070bb','orange','black']

    accumulated_p_RR = np.zeros(50)
    accumulated_p_SF = np.zeros(50)

    p_SF_all = []
    p_RR_all = []

    ite_number = 10
    for index in range(1, ite_number+1):
        file_name = f'EDF7_data/StructureCoefficient_RR_DB_30_degree{deg}_{index}_M{M_value}_u_variation.txt'
        data = np.loadtxt(file_name)
        gamma_values, u_values, lambda1_values, lambda2_values, lambda3_values = data.T

        filtered_indices = np.where((gamma_values == 1.0))
        u_values_filtered = u_values[filtered_indices]
        lambda1_filtered = lambda1_values[filtered_indices]
        lambda2_filtered = lambda2_values[filtered_indices]
        lambda3_filtered = lambda3_values[filtered_indices]

        p_values = []
        for l1, l2, l3 in zip(lambda1_filtered, lambda2_filtered, lambda3_filtered):
            p = dictatorGame_pew_frequency_continual(l1, l2, l3, selection_intensity)
            p_values.append(p)

        accumulated_p_RR += p_values

        p_RR_all.append(p_values)
        file_name = f'EDF7_data/StructureCoefficient_SF_DB_30_degree{deg}_{index}_M{M_value}_u_variation.txt'
        data = np.loadtxt(file_name)
        gamma_values, u_values, lambda1_values, lambda2_values, lambda3_values = data.T

        filtered_indices = np.where((gamma_values == gamma))
        u_values_filtered_SF = u_values[filtered_indices]
        lambda1_filtered = lambda1_values[filtered_indices]
        lambda2_filtered = lambda2_values[filtered_indices]
        lambda3_filtered = lambda3_values[filtered_indices]

        p_values_SF = []
        for l1, l2, l3 in zip(lambda1_filtered, lambda2_filtered, lambda3_filtered):
            p = dictatorGame_pew_frequency_continual(l1, l2, l3, selection_intensity)
            p_values_SF.append(p)

        accumulated_p_SF += p_values_SF
        p_SF_all.append(p_values_SF)

    accumulated_p_RR /= ite_number
    accumulated_p_SF /= ite_number

    p_RR_all = np.array(p_RR_all)
    p_SF_all = np.array(p_SF_all)

    error_margin_RR_p = np.std(p_RR_all,axis=0)  
    error_margin_SF_p = np.std(p_SF_all,axis=0)  

    plt.fill_between(u_values_filtered, accumulated_p_RR - error_margin_RR_p, accumulated_p_RR + error_margin_RR_p,
                     color=colors[1], alpha=0.2)

    plt.fill_between(u_values_filtered, accumulated_p_SF - error_margin_SF_p, accumulated_p_SF + error_margin_SF_p,
                     color=colors[0], alpha=0.2)

    plt.plot(u_values_filtered, accumulated_p_RR, label=rf'RR, $\langle p \rangle$', color=colors[1], linewidth=2.0)
    plt.plot(u_values_filtered, accumulated_p_SF, label=rf'BA, $\langle p \rangle$', color=colors[0], linewidth=2.0)

    plt.xlim(1e-3, 1)
    plt.xscale('log')

    plt.xlabel(r'Mutation rate, $\mu$', fontsize = 14)
    plt.ylabel('Mean value', fontsize = 14)
    plt.title(fr'$\overline{{M}}={M_value}, \gamma={gamma}$')
    plt.legend(loc='best', fontsize = 12)

    ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax.ticklabel_format(style='scientific', axis='y', scilimits=(0, 0))

    # plt.savefig(rf"ExtendedDataFig7_DTG_prFrequency_vs_M_SF_RR_M{M_value}_gamma{gamma}_delta{selection_intensity}_line.svg")
    plt.show()


def plot_pq_vs_gamma_variousNetwork(file_pattern, deg_values, M_numbers,  mutation_rates,selection_intensity=0.005):
    plt.rcParams['font.family'] = 'Arial'
    plt.rc('font', size=12.5)
    fig, ax = plt.subplots(figsize=(3.25, 3.25))
    colors = ['#c1272c', '#0070bb']  
    M_number = M_numbers[0]  

    for i, mutation_rate in enumerate(mutation_rates):
        for deg in deg_values:
            file_name = file_pattern.format(deg=deg, M_number=M_number)
            data = np.loadtxt(file_name)
            gamma_values, u_values, lambda1_values, lambda2_values, lambda3_values = data.T

            filtered_indices = np.where(u_values == mutation_rate)
            gamma_values_filtered = gamma_values[filtered_indices]
            lambda1_filtered = lambda1_values[filtered_indices]
            lambda2_filtered = lambda2_values[filtered_indices]
            lambda3_filtered = lambda3_values[filtered_indices]

            phi_values = []
            for l1, l2, l3 in zip(lambda1_filtered, lambda2_filtered, lambda3_filtered):
                phi = dictatorGame_pew_frequency_continual(l1, l2, l3, selection_intensity)
                phi_values.append(phi)

            plt.plot(gamma_values_filtered, phi_values, linewidth=2.0, label=rf'SF, $\langle p \rangle$', color=colors[i], linestyle='-')

    file_name = "EDF7_data/SC_RR_DB_30_degree4_1_M_variation.txt"
    data = np.loadtxt(file_name)
    gamma_values, u_values, M_values, lambda1_values, lambda2_values, lambda3_values = data.T

    filtered_indices = np.where((u_values == mutation_rate) & (M_values == M_number))
    lambda1_filtered = lambda1_values[filtered_indices]
    lambda2_filtered = lambda2_values[filtered_indices]
    lambda3_filtered = lambda3_values[filtered_indices]

    phi_values = []

    for l1, l2, l3 in zip(lambda1_filtered, lambda2_filtered, lambda3_filtered):
        phi = dictatorGame_pew_frequency_continual(l1, l2, l3, selection_intensity)
        phi_values.append(phi)

    phi_values = np.repeat(phi_values, len(gamma_values_filtered))
    plt.plot(gamma_values_filtered, phi_values, linewidth=2.0, label=rf'RR, $\langle p \rangle$', color=colors[1], linestyle='-')

    plt.xlim(-4,4)

    ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax.ticklabel_format(style='scientific', axis='y', scilimits=(0, 0))

    plt.axvline(x=1, linewidth = 1.5, color='grey', linestyle='--')

    plt.xlabel(r'$\gamma$', fontsize = 14)
    plt.ylabel('Mean value', fontsize = 14)
    plt.title(fr'$\overline{{M}}={M_number}$', fontsize = 14)
    plt.legend(fontsize = 12)

    # plt.savefig(rf"ExtendedDataFig7_DTG_gamma_heterogeneity_degree{deg_values}_M{M_number}_delta{selection_intensity}.svg")
    plt.show()


plot_pq_vs_u_withSimulation("EDF7_data/StructureCoefficient_SF_DB_30_degree{deg}_1_M{M_number}_u_variation.txt", [4], [30,  50,], 1.0,0.005)
plot_pq_vs_gamma_variousMutationRate("EDF7_data/StructureCoefficient_SF_DB_30_degree{deg}_1_M{M_number}_gamma_variation.txt", [4], [50], [0.01,0.1], 0.005)
plot_pq_vs_M_SF_RR_variousMutationRate_line(4, 50, 1.0, 0.005)
plot_pq_vs_gamma_variousNetwork("EDF7_data/StructureCoefficient_SF_DB_30_degree{deg}_1_M{M_number}_gamma_variation.txt", [4], [50], [0.01], 0.005)