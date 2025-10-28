import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import matplotlib.patches as patches

def trustGame_pr_frequency_continual(lambda1, lambda2, lambda3, selection_intensity, b):
    p = 1/2 + selection_intensity * ( lambda1 * 1/12 * (b-1)     - lambda2 * 1/12 + lambda3 * (b-2)/24)
    r = 1/2 + selection_intensity * ( lambda2 * (-b/12)  - lambda3*(b/24) )
    return p,r

def plot_pq_vs_u_withSimulation(file_pattern, deg_values, M_numbers, b, gamma, selection_intensity=0.005):
    plt.rcParams['font.family'] = 'Arial'
    plt.rc('font', size=12.5)

    fig, ax = plt.subplots(figsize=(3.25, 3.25))
    colors = ['#c1272c', '#4758A2', '#0070bb', '#c1272c']  

    for i, M_number in enumerate(M_numbers):
        if (M_number == 40) or (M_number == 60):
            continue
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
            r_values = []
            for l1, l2, l3 in zip(lambda1_filtered, lambda2_filtered, lambda3_filtered):
                p, r = trustGame_pr_frequency_continual(l1, l2, l3, selection_intensity, b)
                p_values.append(p)
                r_values.append(r)

            plt.plot(u_values_filtered, p_values, linewidth=2.0, label=rf'M={M_number}, $\langle p \rangle$', color=colors[i], linestyle='-')
            plt.plot(u_values_filtered, r_values, linewidth=2.0, label=rf'M={M_number}, $\langle r \rangle$', color=colors[i], linestyle='--')

    ite = 50
    data30 = np.zeros((5,2))
    data50 = np.zeros((5,2))

    for index in range(1,ite+1):
        filename30 = f'EDF3_data/data_figure4_panel_a_M30_gamma1.0_selectionIntenstiy0.005_100000000times_{index}.txt'
        filename50 = f'EDF3_data/data_figure4_panel_a_M50_gamma1.0_selectionIntenstiy0.005_100000000times_{index}.txt'
        temp30 = np.loadtxt(filename30)
        temp50 = np.loadtxt(filename50)
        data30 += temp30
        data50 += temp50

    for index in range(1,6):
        filename30_ori = f'EDF3_data/data_figure4_panel_a_M30_gamma1.0_selectionIntenstiy0.005_1000000000times_{index}.txt'
        filename50_ori = f'EDF3_data/data_figure4_panel_a_M50_gamma1.0_selectionIntenstiy0.005_1000000000times_{index}.txt'
        temp30_ori = np.loadtxt(filename30_ori)
        temp50_ori = np.loadtxt(filename50_ori)
        data30 += temp30_ori
        data50 += temp50_ori

    data30 /= ite * 100000000 * 30 * 2
    data50 /= ite * 100000000 * 50 * 2 

    u_set = [10**(-2.5), 10**(-2), 10**(-1.5), 10**(-1), 10**(-0.5)]
    plt.scatter(u_set, data30[:,0], marker='s', color=colors[0])
    plt.scatter(u_set, data30[:,1], marker='o', color=colors[0])
    plt.scatter(u_set, data50[:,0], marker='s', color=colors[2])
    plt.scatter(u_set, data50[:,1], marker='o', color=colors[2])

    ite_addition=10 + 10
    data30_addition = np.zeros((10,2))
    data50_addition = np.zeros((10,2))

    for index in range(1,ite_addition+1):
        filename30_addition = f'EDF3_data/data_figure4_panel_a_M30_gamma1.0_selectionIntenstiy0.005_500000000times_{index}.txt'
        filename50_addition = f'EDF3_data/data_figure4_panel_a_M50_gamma1.0_selectionIntenstiy0.005_500000000times_{index}.txt'

        temp30_addition = np.loadtxt(filename30_addition)
        temp50_addition = np.loadtxt(filename50_addition)

        data30_addition += temp30_addition
        data50_addition += temp50_addition

    data30_addition /= ite_addition * 500000000 * 30
    data50_addition /= ite_addition * 500000000 * 50

    u_set_addition = [10**(-2.5-0.5/3), 10**(-2.5+0.5/3), 10**(-2.5 + 1.0/3), 10**(-2+0.5/3), 10**(-2+1.0/3), 10**(-1.5+0.5/3), 10**(-1.5+1.0/3), 10**(-1+0.5/3), 10**(-1+1.0/3), 10**(-0.5+0.5/3)]
    plt.scatter(u_set_addition, data30_addition[:,0], marker='s', color=colors[0])
    plt.scatter(u_set_addition, data30_addition[:,1], marker='o', color=colors[0])
    plt.scatter(u_set_addition, data50_addition[:,0], marker='s', color=colors[2])
    plt.scatter(u_set_addition, data50_addition[:,1], marker='o', color=colors[2])

    ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax.ticklabel_format(style='scientific', axis='y', scilimits=(0, 0))

    plt.ylim(0.4735,0.5015)

    plt.xlim(1e-3, 1)

    plt.xscale('log')
    plt.xlabel(r'Mutation rate, $\mu$', fontsize = 14)
    plt.ylabel('Mean value', fontsize = 14)
    plt.title(fr'$\gamma={gamma}$', fontsize = 14)
    plt.legend(loc='lower right')
    # plt.savefig(rf"ExtendedDataFig3_TG_prFrequency_vs_mutationRate_SF_degree{deg_values}_b{b}_gamma{gamma}_delta{selection_intensity}_new.svg")
    plt.show()

def plot_pq_vs_gamma_variousMutationRate(file_pattern, deg_values, M_numbers, b, mutation_rates,selection_intensity=0.005):
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
            omega_values = []
            for l1, l2, l3 in zip(lambda1_filtered, lambda2_filtered, lambda3_filtered):
                phi,omega = trustGame_pr_frequency_continual(l1, l2, l3, selection_intensity, b)
                phi_values.append(phi)
                omega_values.append(omega)

            plt.plot(gamma_values_filtered, phi_values, linewidth=2.0, label=rf'$\mu$={mutation_rate}, $\langle \varphi \rangle$', color=colors[i], linestyle='-')
            plt.plot(gamma_values_filtered, omega_values, linewidth=2.0, label=rf'$\mu$={mutation_rate}, $\langle \omega \rangle$', color=colors[i], linestyle='--')

            if M_number == 30:
                mean_phi_value = np.mean(phi_values)
                mean_omega_value = np.mean(omega_values)
                plt.axhline(mean_phi_value, color=colors[i], linestyle='-', linewidth=2.0)
                plt.axhline(mean_omega_value, color=colors[i], linestyle='--', linewidth=2.0)

    plt.xlim(-4,4)

    ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax.ticklabel_format(style='scientific', axis='y', scilimits=(0, 0))

    plt.axvline(x=1, linewidth = 1.5, color='grey', linestyle='--')

    plt.xlabel(r'$\gamma$', fontsize = 14)
    plt.ylabel('Mean value', fontsize = 14)
    plt.title(fr'$\overline{{M}}={M_number}$', fontsize = 14)
    plt.legend(fontsize = 12)

    # plt.savefig(rf"ExtendedDataFig3_TG_pqFrequency_vs_gamma_verious_u_SF_degree{deg_values}_b{b}_M{M_number}_delta{selection_intensity}.svg")
    plt.show()

def plot_pq_vs_M_SF_RR_variousMutationRate_line(deg, b, M_value, gamma, selection_intensity=0.005):
    plt.rcParams['font.family'] = 'Arial'
    plt.rc('font', size=12.5)

    fig, ax = plt.subplots(figsize=(3.25, 3.25))

    colors = ['#c1272c', '#0070bb','orange','black']

    accumulated_p_RR = np.zeros(50)
    accumulated_p_SF = np.zeros(50)
    accumulated_r_RR = np.zeros(50)
    accumulated_r_SF = np.zeros(50)

    p_SF_all = []
    r_SF_all = []
    p_RR_all = []
    r_RR_all = []

    ite_number = 10
    for index in range(1, ite_number+1):
        file_name = f'EDF3_data/StructureCoefficient_RR_DB_30_degree{deg}_{index}_M{M_value}_u_variation.txt'
        data = np.loadtxt(file_name)
        gamma_values, u_values, lambda1_values, lambda2_values, lambda3_values = data.T

        filtered_indices = np.where((gamma_values == 1.0))
        u_values_filtered = u_values[filtered_indices]
        lambda1_filtered = lambda1_values[filtered_indices]
        lambda2_filtered = lambda2_values[filtered_indices]
        lambda3_filtered = lambda3_values[filtered_indices]

        p_values, r_values = [], []
        for l1, l2, l3 in zip(lambda1_filtered, lambda2_filtered, lambda3_filtered):
            p, r = trustGame_pr_frequency_continual(l1, l2, l3, selection_intensity, b)
            p_values.append(p)
            r_values.append(r)

        accumulated_p_RR += p_values
        accumulated_r_RR += r_values

        p_RR_all.append(p_values)
        r_RR_all.append(r_values)

        file_name = f'EDF3_data/StructureCoefficient_SF_DB_30_degree{deg}_{index}_M{M_value}_u_variation.txt'
        data = np.loadtxt(file_name)
        gamma_values, u_values, lambda1_values, lambda2_values, lambda3_values = data.T

        filtered_indices = np.where((gamma_values == gamma))
        u_values_filtered_SF = u_values[filtered_indices]
        lambda1_filtered = lambda1_values[filtered_indices]
        lambda2_filtered = lambda2_values[filtered_indices]
        lambda3_filtered = lambda3_values[filtered_indices]

        p_values_SF, r_values_SF = [], []
        for l1, l2, l3 in zip(lambda1_filtered, lambda2_filtered, lambda3_filtered):
            p, r = trustGame_pr_frequency_continual(l1, l2, l3, selection_intensity, b)
            p_values_SF.append(p)
            r_values_SF.append(r)

        accumulated_p_SF += p_values_SF
        accumulated_r_SF += r_values_SF

        p_SF_all.append(p_values_SF)
        r_SF_all.append(r_values_SF)

    accumulated_p_RR /= ite_number
    accumulated_p_SF /= ite_number
    accumulated_r_RR /= ite_number
    accumulated_r_SF /= ite_number

    p_SF_all = np.array(p_SF_all)
    r_SF_all = np.array(r_SF_all)
    p_RR_all = np.array(p_RR_all)
    r_RR_all = np.array(r_RR_all)

    error_margin_RR_p = np.std(p_RR_all, axis=0)  
    error_margin_SF_p = np.std(p_SF_all, axis=0) 
    error_margin_RR_r = np.std(r_RR_all, axis=0) 
    error_margin_SF_r = np.std(r_SF_all, axis=0)  

    plt.fill_between(u_values_filtered, accumulated_p_RR - error_margin_RR_p, accumulated_p_RR + error_margin_RR_p,
                     color=colors[1], alpha=0.2)
    plt.fill_between(u_values_filtered, accumulated_r_RR - error_margin_RR_r, accumulated_r_RR + error_margin_RR_r,
                     color=colors[1], alpha=0.2)

    plt.fill_between(u_values_filtered, accumulated_p_SF - error_margin_SF_p, accumulated_p_SF + error_margin_SF_p,
                     color=colors[0], alpha=0.2)
    plt.fill_between(u_values_filtered, accumulated_r_SF - error_margin_SF_r, accumulated_r_SF + error_margin_SF_r,
                     color=colors[0], alpha=0.2)

    plt.plot(u_values_filtered, accumulated_p_RR, label=rf'RR, $\langle p \rangle$', color=colors[1], linewidth=2.0)
    plt.plot(u_values_filtered, accumulated_r_RR, label=rf'RR, $\langle r \rangle$', color=colors[1], linewidth=2.0, linestyle='--')
    plt.plot(u_values_filtered, accumulated_p_SF, label=rf'BA, $\langle p \rangle$', color=colors[0], linewidth=2.0)
    plt.plot(u_values_filtered, accumulated_r_SF, label=rf'BA, $\langle r \rangle$', color=colors[0], linewidth=2.0, linestyle='--')

    plt.xlim(1e-3, 1)
    plt.xscale('log')

    plt.yticks([0.47, 0.48, 0.49, 0.50])
    plt.xlabel(r'Mutation rate, $\mu$', fontsize = 14)
    plt.ylabel('Mean value', fontsize = 14)
    plt.title(fr'$\overline{{M}}={M_value}, \gamma={gamma}$')
    plt.legend(loc='best', fontsize = 12)

    ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax.ticklabel_format(style='scientific', axis='y', scilimits=(0, 0))

    # plt.savefig(rf"ExtendedDataFig3_TG_prFrequency_vs_M_SF_RR_b{b}_M{M_value}_gamma{gamma}_delta{selection_intensity}_line.svg")
    plt.show()

def plot_pq_vs_gamma_variousNetwork(file_pattern, deg_values, M_numbers, b, mutation_rates,selection_intensity=0.005):
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
            omega_values = []
            for l1, l2, l3 in zip(lambda1_filtered, lambda2_filtered, lambda3_filtered):
                phi,omega = trustGame_pr_frequency_continual(l1, l2, l3, selection_intensity, b)
                phi_values.append(phi)
                omega_values.append(omega)

            plt.plot(gamma_values_filtered, phi_values, linewidth=2.0, label=rf'SF, $\langle p \rangle$', color=colors[i], linestyle='-')
            plt.plot(gamma_values_filtered, omega_values, linewidth=2.0, label=rf'SF, $\langle r \rangle$', color=colors[i], linestyle='--')

    file_name = "EDF3_data/SC_RR_DB_30_degree4_1_M_variation.txt"
    data = np.loadtxt(file_name)
    gamma_values, u_values, M_values, lambda1_values, lambda2_values, lambda3_values = data.T

    filtered_indices = np.where((u_values == mutation_rate) & (M_values == M_number))
    lambda1_filtered = lambda1_values[filtered_indices]
    lambda2_filtered = lambda2_values[filtered_indices]
    lambda3_filtered = lambda3_values[filtered_indices]

    phi_values = []
    omega_values = []
    for l1, l2, l3 in zip(lambda1_filtered, lambda2_filtered, lambda3_filtered):
        phi,omega = trustGame_pr_frequency_continual(l1, l2, l3, selection_intensity, b)
        phi_values.append(phi)
        omega_values.append(omega)

    phi_values = np.repeat(phi_values, len(gamma_values_filtered))
    omega_values = np.repeat(omega_values, len(gamma_values_filtered))

    plt.plot(gamma_values_filtered, phi_values, linewidth=2.0, label=rf'RR, $\langle p \rangle$', color=colors[1], linestyle='-')
    plt.plot(gamma_values_filtered, omega_values, linewidth=2.0, label=rf'RR, $\langle r \rangle$', color=colors[1], linestyle='--')

    plt.xlim(-4,4)

    ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax.ticklabel_format(style='scientific', axis='y', scilimits=(0, 0))

    plt.axvline(x=1, linewidth = 1.5, color='grey', linestyle='--')

    plt.xlabel(r'$\gamma$', fontsize = 14)
    plt.ylabel('Mean value', fontsize = 14)
    plt.title(fr'$\overline{{M}}={M_number}$', fontsize = 14)
    plt.legend(fontsize = 12)
    # plt.savefig(rf"ExtendedDataFig3_TG_gamma_heterogeneity_degree{deg_values}_b{b}_M{M_number}_delta{selection_intensity}.svg")
    plt.show()

plot_pq_vs_u_withSimulation("EDF3_data/StructureCoefficient_SF_DB_30_degree{deg}_1_M{M_number}_u_variation.txt", [4], [30, 40, 50, 60], 1.5 ,1.0,0.005) # panel a

plot_pq_vs_gamma_variousMutationRate("EDF3_data/StructureCoefficient_SF_DB_30_degree{deg}_1_M{M_number}_gamma_variation.txt", [4], [50], 1.5 ,[0.01,0.1], 0.005) # panel b

plot_pq_vs_M_SF_RR_variousMutationRate_line(4, 1.5, 50, 1.0, 0.005) # panel c

plot_pq_vs_gamma_variousNetwork("EDF3_data/StructureCoefficient_SF_DB_30_degree{deg}_1_M{M_number}_gamma_variation.txt", [4], [50], 1.5 ,[0.01], 0.005) # panel d


