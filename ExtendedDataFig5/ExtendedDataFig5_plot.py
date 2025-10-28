import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import matplotlib.patches as patches

def dictatorGame_pew_frequency_continual(lambda1, lambda2, lambda3, selection_intensity):
    p = 1/2 + selection_intensity * (- lambda2 * 1/6 - lambda3 * 1/12)
    return p

def plot_pq_vs_u_withSimulation(file_pattern, deg_values, M_numbers, gamma, sigma, selection_intensity=0.005):
    plt.rcParams['font.family'] = 'Arial'
    plt.rc('font', size=12.5)

    fig, ax = plt.subplots(figsize=(3.25, 3.25))
    colors = ['#c1272c',  '#0070bb']  

    for i, M_number in enumerate(M_numbers):
        for deg in deg_values:
            file_name = file_pattern.format(deg=deg, M_number=M_number)
            data = np.loadtxt(file_name)
            gamma_values, u_values, K1_values, K2_values, bc_values = data.T

            filtered_indices = np.where((gamma_values == 1.0))
            u_values_filtered = u_values[filtered_indices]
            K1_filtered = K1_values[filtered_indices]
            K2_filtered = K2_values[filtered_indices]
            
            R1 = (K1_filtered + K2_filtered)/4
            R2 = (K1_filtered - K2_filtered)/4
            p_values = 1/2 + selection_intensity /2 *1*( (sigma-1)*R1 -R2 )

            plt.plot(u_values_filtered, p_values, linewidth=2.0, label=rf'$\overline{{M}}={M_number/50}$', color=colors[i], linestyle='-')
        
    ite = 30
    data50 = np.zeros((5,1))
    data100 = np.zeros((5,1))

    for index in range(1,ite+1):
        filename50 = f'EDF5_data/SHG_M50_selectionIntenstiy0.005_100000000times_{index}.txt'
        filename100 = f'EDF5_data/SHG_M100_gamma1.0_selectionIntenstiy0.005_100000000times_{index}.txt'

        temp50 = np.loadtxt(filename50)
        temp100 = np.loadtxt(filename100)
        
        data50 += temp50[:,1][:, np.newaxis]
        data100 += temp100[:,1][:, np.newaxis]

    data50 /= ite * 100000000 * 50
    data100 /= ite * 100000000 * 100

    u_set = [10**(-2.5), 10**(-2), 10**(-1.5), 10**(-1), 10**(-0.5)]
    plt.scatter(u_set, data50[:,0], marker='s', color=colors[0])
    plt.scatter(u_set, data100[:,0], marker='s', color=colors[1])

    plt.xlim(1e-3,1)
    plt.xscale('log')
    plt.xlabel(r'Mutation rate, $\mu$', fontsize = 14)
    plt.ylabel('Mean value', fontsize = 14)
    plt.title(fr'$\bar{{d}} = {deg}, \gamma = {gamma}$', fontsize = 14)
    plt.legend(loc='upper right')

    ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax.ticklabel_format(style='scientific', axis='y', scilimits=(0, 0))

    # plt.savefig(rf"ExtendedDataFig5_SHG_pqFrequency_Simulation_sigma{sigma}.svg")
    plt.show()

def plot_pq_vs_gamma_variousMutationRate(file_name, deg_values, M_numbers,  mutation_rates,sigma,selection_intensity=0.005):
    plt.rcParams['font.family'] = 'Arial'
    plt.rc('font', size=12.5)

    fig, ax = plt.subplots(figsize=(3.25, 3.25))

    colors = ['#c1272c', '#0070bb']  
    M_number = M_numbers[0]  

    for i, mutation_rate in enumerate(mutation_rates):
        data = np.loadtxt(file_name)
        gamma_values, u_values, K1_values, K2_values, bc_values = data.T

        filtered_indices = np.where(u_values == mutation_rate)
        gamma_values_filtered = gamma_values[filtered_indices]
        K1_filtered = K1_values[filtered_indices]
        K2_filtered = K2_values[filtered_indices]

        R1 = (K1_filtered + K2_filtered)/4
        R2 = (K1_filtered - K2_filtered)/4

        phi_values = 1/2 + selection_intensity /2 *1*( (sigma-1)*R1 -R2 )

        plt.plot(gamma_values_filtered, phi_values, linewidth=2.0, label=rf'$\mu$={mutation_rate}', color=colors[i], linestyle='-')

        if M_number == 30:
            mean_phi_value = np.mean(phi_values)
            plt.axhline(mean_phi_value, color=colors[i], linestyle='-', linewidth=2.0)

    plt.xlim(-4,4)

    ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax.ticklabel_format(style='scientific', axis='y', scilimits=(0, 0))

    plt.axvline(x=1, linewidth = 1.5, color='grey', linestyle='--')

    plt.xlabel(r'$\gamma$', fontsize = 14)
    plt.ylabel(r'$\langle h \rangle$', fontsize = 14)
    plt.title(fr'$\overline{{M}}={M_number}$', fontsize = 14)
    plt.legend(fontsize = 12)

    # # plt.savefig(rf"ExtendedDataFig5_SHG_pqFrequency_vs_gamma_verious_u_SF_degree{deg_values}_M{M_number}_delta{selection_intensity}.svg")
    plt.show()

def plot_pq_vs_M_SF_RR_variousMutationRate_line(deg, M_value, gamma, sigma, selection_intensity=0.005):
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
        file_name = f'EDF5_data/DonationGame_RR_DB_bcCritical_50_degree{deg}_{index}_M{M_value}.txt'
        data = np.loadtxt(file_name)
        gamma_values, u_values, K1_values, K2_values, bc_values = data.T

        filtered_indices = np.where((gamma_values == 1.0))
        u_values_filtered = u_values[filtered_indices]
        K1_filtered = K1_values[filtered_indices]
        K2_filtered = K2_values[filtered_indices]
        

        R1 = (K1_filtered + K2_filtered)/4
        R2 = (K1_filtered - K2_filtered)/4
        p_values = 1/2 + selection_intensity /2 *1*( (sigma-1)*R1 -R2 )

        p_RR_all.append(p_values)

        accumulated_p_RR += p_values

        file_name = f'EDF5_data/DonationGame_SF_DB_bcCritical_50_degree{deg}_{index}_M{M_value}.txt'
        data = np.loadtxt(file_name)
        gamma_values, u_values, K1_values, K2_values, bc_values = data.T

        filtered_indices = np.where((gamma_values == gamma))
        u_values_filtered = u_values[filtered_indices]
        K1_filtered = K1_values[filtered_indices]
        K2_filtered = K2_values[filtered_indices]
        

        R1 = (K1_filtered + K2_filtered)/4
        R2 = (K1_filtered - K2_filtered)/4
        p_values_SF = 1/2 + selection_intensity /2 *1*( (sigma-1)*R1 -R2 )

        accumulated_p_SF += p_values_SF

        p_SF_all.append(p_values_SF)

    p_values = np.array(p_RR_all)
    p_values_SF = np.array(p_SF_all)

    accumulated_p_RR /= ite_number
    accumulated_p_SF /= ite_number

    error_margin_RR_p = np.std(p_RR_all, axis=0)  
    error_margin_SF_p = np.std(p_SF_all, axis=0)  

    plt.fill_between(u_values_filtered, accumulated_p_RR - error_margin_RR_p, accumulated_p_RR + error_margin_RR_p,
                     color=colors[1], alpha=0.2)

    plt.fill_between(u_values_filtered, accumulated_p_SF - error_margin_SF_p, accumulated_p_SF + error_margin_SF_p,
                     color=colors[0], alpha=0.2)

    plt.plot(u_values_filtered, accumulated_p_RR, label=rf'RR', color=colors[1], linewidth=2.0)
    plt.plot(u_values_filtered, accumulated_p_SF, label=rf'SF', color=colors[0], linewidth=2.0)

    plt.xlim(1e-3, 1)
    plt.xscale('log')

    plt.xlabel(r'Mutation rate, $\mu$', fontsize = 14)
    plt.ylabel(r'$\langle h \rangle$', fontsize = 14)
    plt.title(fr'$\overline{{M}}={M_value}, \gamma={gamma}$')
    plt.legend(loc='best', fontsize = 12)

    ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax.ticklabel_format(style='scientific', axis='y', scilimits=(0, 0))

    # plt.savefig(rf"ExtendedDataFig5_SHG_prFrequency_vs_M_SF_RR_M{M_value}_gamma{gamma}_delta{selection_intensity}_line.svg")
    plt.show()

def plot_pq_vs_gamma_variousNetwork(file_name, deg_values, M_numbers,  mutation_rates,sigma,selection_intensity=0.005):
    plt.rcParams['font.family'] = 'Arial'
    plt.rc('font', size=12.5)
    fig, ax = plt.subplots(figsize=(3.25, 3.25))
    colors = ['#c1272c', '#0070bb']  
    M_number = M_numbers[0]  

    for i, mutation_rate in enumerate(mutation_rates):
        data = np.loadtxt(file_name)
        gamma_values, u_values, K1_values, K2_values, bc_values = data.T

        filtered_indices = np.where(u_values == mutation_rate)
        gamma_values_filtered = gamma_values[filtered_indices]
        K1_filtered = K1_values[filtered_indices]
        K2_filtered = K2_values[filtered_indices]

        R1 = (K1_filtered + K2_filtered)/4
        R2 = (K1_filtered - K2_filtered)/4

        phi_values = 1/2 + selection_intensity /2 *1*( (sigma-1)*R1 -R2 )

        plt.plot(gamma_values_filtered, phi_values, linewidth=2.0, label=rf'SF', color=colors[i], linestyle='-')

    k1_rr = 94.3289017244 
    k2_rr = 28.3222623592

    r1_RR = (k1_rr+k2_rr)/4
    r2_RR = (k1_rr-k2_rr)/4

    values_RR = 1/2+selection_intensity /2 *1*( (sigma-1)*r1_RR -r2_RR )

    plt.axhline(y=values_RR, color=colors[1], linestyle='-', linewidth=2.0, label=rf'RR')

    plt.xlim(-4,4)

    ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax.ticklabel_format(style='scientific', axis='y', scilimits=(0, 0))

    plt.axvline(x=1, linewidth = 1.5, color='grey', linestyle='--')

    plt.xlabel(r'$\gamma$', fontsize = 14)
    plt.ylabel(r'$\langle h \rangle$', fontsize = 14)
    plt.title(fr'$\overline{{M}}={M_number}$', fontsize = 14)
    plt.legend(fontsize = 12)

    # plt.savefig(rf"ExtendedDataFig5_SHG_gamma_heterogeneity_degree{deg_values}_M{M_number}_delta{selection_intensity}.svg")
    plt.show()

sigma = 1.5

plot_pq_vs_u_withSimulation("EDF5_data/SF_DB_bcCritical_50_degree{deg}_1_M{M_number}.txt", [4], [50, 100], 1.0,sigma, 0.005)

plot_pq_vs_gamma_variousMutationRate("EDF5_data/SF_DB_bcCritical_50_degree4_1_M100_gammaVariation.txt", [4],[100], [0.01,0.1], sigma,  0.005)

plot_pq_vs_M_SF_RR_variousMutationRate_line(4, 100, 1.0, sigma,  0.005)

plot_pq_vs_gamma_variousNetwork("EDF5_data/SF_DB_bcCritical_50_degree4_1_M100_gammaVariation.txt", [4],[100], [0.01], sigma,  0.005)
