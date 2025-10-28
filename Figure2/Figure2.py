import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import matplotlib.patches as patches

def trustGame_pr_frequency_continual(lambda1, lambda2, lambda3, selection_intensity, b):
    p = 1/2 + selection_intensity * ( lambda1 * 1/12 * (b-1)     - lambda2 * 1/12 + lambda3 * (b-2)/24)
    r = 1/2 + selection_intensity * ( lambda2 * (-b/12)  - lambda3*(b/24) )
    return p,r
def plot_M_trustGame(deg, b, mutation_rate, gamma, selection_intensity=0.005):
    plt.rcParams['font.family'] = 'Arial'
    plt.rc('font', size=12.5)

    fig, ax = plt.subplots(figsize=(3.25, 3.25))

    colors = ['#c1272c', '#0070bb','orange','black']

    accumulated_p_RR = np.zeros(151)
    accumulated_p_SF = np.zeros(151)
    accumulated_r_RR = np.zeros(151)
    accumulated_r_SF = np.zeros(151)

    p_SF_all = []
    r_SF_all = []
    p_RR_all = []
    r_RR_all = []

    ite_number = 100
    for index in range(1, ite_number+1):
        # RR Data
        file_name = f'data_M_vs_lambda1_2_3/SC_RR_DB_50_degree{deg}_{index}_M_variation.txt'
        data = np.loadtxt(file_name)
        gamma_values, u_values, M_values, lambda1_values, lambda2_values, lambda3_values = data.T

        filtered_indices = np.where((u_values == mutation_rate) & (gamma_values == 1.0))
        M_values_filtered = M_values[filtered_indices]
        lambda1_filtered = lambda1_values[filtered_indices]
        lambda2_filtered = lambda2_values[filtered_indices]
        lambda3_filtered = lambda3_values[filtered_indices]

        p_values, r_values = [], []
        for l1, l2, l3 in zip(lambda1_filtered, lambda2_filtered, lambda3_filtered):
            p, r = trustGame_pr_frequency_continual(l1, l2, l3, selection_intensity, b)
            p_values.append(p)
            r_values.append(r)

        # Accumulate values
        accumulated_p_RR += p_values
        accumulated_r_RR += r_values

        p_RR_all.append(p_values)
        r_RR_all.append(r_values)

        # SF Data
        file_name = f'data_M_vs_lambda1_2_3/SC_SF_DB_50_degree{deg}_{index}_gamma1.0_M_variation.txt'
        data = np.loadtxt(file_name)
        gamma_values, u_values, M_values, lambda1_values, lambda2_values, lambda3_values = data.T

        filtered_indices = np.where((u_values == mutation_rate) & (gamma_values == gamma))
        M_values_filtered_SF = M_values[filtered_indices]
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

    # Average values
    accumulated_p_RR /= ite_number
    accumulated_p_SF /= ite_number
    accumulated_r_RR /= ite_number
    accumulated_r_SF /= ite_number

    p_RR_all = np.array(p_RR_all)
    r_RR_all = np.array(r_RR_all)
    p_SF_all = np.array(p_SF_all)
    r_SF_all = np.array(r_SF_all)

    # Error margins (for demonstration, adjust accordingly)
    error_margin_RR_p = np.std(p_RR_all, axis=0)  # Example for error calculation
    error_margin_SF_p = np.std(p_SF_all, axis=0)  # Example for error calculation
    error_margin_RR_r = np.std(r_RR_all, axis=0)  # Example for error calculation
    error_margin_SF_r = np.std(r_SF_all, axis=0)  # Example for error calculation

    M_values_filtered /= 50

    # Plotting means with error regions
    plt.fill_between(M_values_filtered, accumulated_p_RR - error_margin_RR_p, accumulated_p_RR + error_margin_RR_p,
                     color=colors[1], alpha=0.2)
    plt.fill_between(M_values_filtered, accumulated_r_RR - error_margin_RR_r, accumulated_r_RR + error_margin_RR_r,
                     color=colors[1], alpha=0.2)

    plt.fill_between(M_values_filtered, accumulated_p_SF - error_margin_SF_p, accumulated_p_SF + error_margin_SF_p,
                     color=colors[0], alpha=0.2)
    plt.fill_between(M_values_filtered, accumulated_r_SF - error_margin_SF_r, accumulated_r_SF + error_margin_SF_r,
                     color=colors[0], alpha=0.2)

    # Solid lines for mean values
    plt.plot(M_values_filtered, accumulated_p_RR, label=rf'RR, $\langle p \rangle$', color=colors[1], linewidth=2.0)
    plt.plot(M_values_filtered, accumulated_r_RR, label=rf'RR, $\langle r \rangle$', color=colors[1], linewidth=2.0, linestyle='--')
    plt.plot(M_values_filtered, accumulated_p_SF, label=rf'SF, $\langle p \rangle$', color=colors[0], linewidth=2.0)
    plt.plot(M_values_filtered, accumulated_r_SF, label=rf'SF, $\langle r \rangle$', color=colors[0], linewidth=2.0, linestyle='--')

    # Customizing the plot
    plt.xlim(1, 4)
    # plt.xticks([1, 1.5, 2])
    plt.ylim(0.4665, 0.5035)
    #plt.yticks([0.48, 0.49, 0.50])
    plt.xlabel(r'$\overline{M}$', fontsize = 14)
    plt.ylabel('Mean value', fontsize = 14)
    plt.title(fr'$u={mutation_rate}, \gamma={gamma}$')
    plt.legend(loc='best', fontsize = 12)

    ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax.ticklabel_format(style='scientific', axis='y', scilimits=(0, 0))

    # plt.savefig(rf"Fig_sixGamex_M_trustGame_prFrequency_vs_M_SF_RR_b{b}_u{mutation_rate}_gamma{gamma}_delta{selection_intensity}_newColor.svg")
    plt.show()


def ultimatumGame_pq_frequency_continual(lambda1, lambda2, lambda3, selection_intensity):
    p = 1/2 + selection_intensity * ( lambda1*1/12     - lambda2*1/12)
    q = 1/2 + selection_intensity * ( lambda1*(-1/12)  - lambda3*(1/24) )
    return p,q

def plot_M_ultimatumGame(deg,  mutation_rate, gamma, selection_intensity=0.005):
    plt.rcParams['font.family'] = 'Arial'
    plt.rc('font', size=12.5)

    fig, ax = plt.subplots(figsize=(3.25, 3.25))
    # fig, ax = plt.subplots(figsize=(1.75, 1.75))

    colors = ['#c1272c', '#0070bb','orange','black']
    #colors = ['#c1272c', '#0070bb','orange','black']

    accumulated_p_RR = np.zeros(151)
    accumulated_p_SF = np.zeros(151)
    accumulated_r_RR = np.zeros(151)
    accumulated_r_SF = np.zeros(151)

    p_RR_all = []
    r_RR_all = []
    p_SF_all = []
    r_SF_all = []

    ite_number = 100
    for index in range(1, ite_number+1):
        # RR Data
        file_name = f'data_M_vs_lambda1_2_3/SC_RR_DB_50_degree{deg}_{index}_M_variation.txt'
        data = np.loadtxt(file_name)
        gamma_values, u_values, M_values, lambda1_values, lambda2_values, lambda3_values = data.T

        filtered_indices = np.where((u_values == mutation_rate) & (gamma_values == 1.0))
        M_values_filtered = M_values[filtered_indices]
        lambda1_filtered = lambda1_values[filtered_indices]
        lambda2_filtered = lambda2_values[filtered_indices]
        lambda3_filtered = lambda3_values[filtered_indices]

        p_values, r_values = [], []
        for l1, l2, l3 in zip(lambda1_filtered, lambda2_filtered, lambda3_filtered):
            p, r = ultimatumGame_pq_frequency_continual(l1, l2, l3, selection_intensity)
            p_values.append(p)
            r_values.append(r)

        # Accumulate values
        accumulated_p_RR += p_values
        accumulated_r_RR += r_values

        p_RR_all.append(p_values)
        r_RR_all.append(r_values)

        # SF Data
        file_name = f'data_M_vs_lambda1_2_3/SC_SF_DB_50_degree{deg}_{index}_gamma1.0_M_variation.txt'
        data = np.loadtxt(file_name)
        gamma_values, u_values, M_values, lambda1_values, lambda2_values, lambda3_values = data.T

        filtered_indices = np.where((u_values == mutation_rate) & (gamma_values == gamma))
        M_values_filtered_SF = M_values[filtered_indices]
        lambda1_filtered = lambda1_values[filtered_indices]
        lambda2_filtered = lambda2_values[filtered_indices]
        lambda3_filtered = lambda3_values[filtered_indices]

        p_values_SF, r_values_SF = [], []
        for l1, l2, l3 in zip(lambda1_filtered, lambda2_filtered, lambda3_filtered):
            p, r = ultimatumGame_pq_frequency_continual(l1, l2, l3, selection_intensity)
            p_values_SF.append(p)
            r_values_SF.append(r)

        accumulated_p_SF += p_values_SF
        accumulated_r_SF += r_values_SF

        p_SF_all.append(p_values_SF)
        r_SF_all.append(r_values_SF)

    # Average values
    accumulated_p_RR /= ite_number
    accumulated_p_SF /= ite_number
    accumulated_r_RR /= ite_number
    accumulated_r_SF /= ite_number

    p_RR_all = np.array(p_RR_all)
    r_RR_all = np.array(r_RR_all)
    p_SF_all = np.array(p_SF_all)
    r_SF_all = np.array(r_SF_all)


    # Error margins (for demonstration, adjust accordingly)
    error_margin_RR_p = np.std(p_RR_all, axis=0)  # Example for error calculation
    error_margin_SF_p = np.std(p_SF_all, axis=0)  # Example for error calculation
    error_margin_RR_q = np.std(r_RR_all, axis=0)  # Example for error calculation
    error_margin_SF_q = np.std(r_SF_all, axis=0)  # Example for error calculation

    M_values_filtered /= 50

    # Plotting means with error regions
    plt.fill_between(M_values_filtered, accumulated_p_RR - error_margin_RR_p, accumulated_p_RR + error_margin_RR_p,
                     color=colors[1], alpha=0.2)
    plt.fill_between(M_values_filtered, accumulated_r_RR - error_margin_RR_q, accumulated_r_RR + error_margin_RR_q,
                     color=colors[1], alpha=0.2)

    plt.fill_between(M_values_filtered, accumulated_p_SF - error_margin_SF_p, accumulated_p_SF + error_margin_SF_p,
                     color=colors[0], alpha=0.2)
    plt.fill_between(M_values_filtered, accumulated_r_SF - error_margin_SF_q, accumulated_r_SF + error_margin_SF_q,
                     color=colors[0], alpha=0.2)

    # Solid lines for mean values
    plt.plot(M_values_filtered, accumulated_p_RR, label=rf'RR, $\langle p \rangle$', color=colors[1], linewidth=2.0)
    plt.plot(M_values_filtered, accumulated_r_RR, label=rf'RR, $\langle q \rangle$', color=colors[1], linewidth=2.0, linestyle='--')
    plt.plot(M_values_filtered, accumulated_p_SF, label=rf'SF, $\langle p \rangle$', color=colors[0], linewidth=2.0)
    plt.plot(M_values_filtered, accumulated_r_SF, label=rf'SF, $\langle q \rangle$', color=colors[0], linewidth=2.0, linestyle='--')


    rect1 = patches.Rectangle((2,0.511), 1.95, 0.004, linewidth=2, edgecolor='black', facecolor='none', linestyle='-')
    # 添加矩形框到图上
    plt.gca().add_patch(rect1)

    # Customizing the plot
    plt.xlim(1, 4)
    # plt.ylim(0.490-0.001, 0.505+0.001)
    # plt.yticks([0.49, 0.495, 0.50, 0.505])

    # plt.xlim(3.8, 4)
    # plt.ylim(0.503, 0.505)

    plt.xlabel(r'$\overline{M}$', fontsize = 14)
    plt.ylabel('Mean value', fontsize = 14)
    plt.title(fr'$u={mutation_rate}, \gamma={gamma}$')
    # plt.legend(loc='best', fontsize = 12)

    ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax.ticklabel_format(style='scientific', axis='y', scilimits=(0, 0))

    # plt.savefig(rf"Fig_sixGamex_M_ultimatumGame_pqFrequency_SF_RR_u{mutation_rate}_gamma{gamma}_delta{selection_intensity}_line_newColor.svg")
    plt.show()

def plot_M_ultimatumGame_subfig(deg,  mutation_rate, gamma, selection_intensity=0.005):
    plt.rcParams['font.family'] = 'Arial'
    plt.rc('font', size=12.5)

    # fig, ax = plt.subplots(figsize=(3.25, 3.25))
    fig, ax = plt.subplots(figsize=(1.75, 1.75))

    colors = ['#c1272c', '#0070bb','orange','black']
    #colors = ['#c1272c', '#0070bb','orange','black']

    accumulated_p_RR = np.zeros(151)
    accumulated_p_SF = np.zeros(151)
    accumulated_r_RR = np.zeros(151)
    accumulated_r_SF = np.zeros(151)

    p_RR_all = []
    r_RR_all = []
    p_SF_all = []
    r_SF_all = []

    ite_number = 100
    for index in range(1, ite_number+1):
        # RR Data
        file_name = f'data_M_vs_lambda1_2_3/SC_RR_DB_50_degree{deg}_{index}_M_variation.txt'
        data = np.loadtxt(file_name)
        gamma_values, u_values, M_values, lambda1_values, lambda2_values, lambda3_values = data.T

        filtered_indices = np.where((u_values == mutation_rate) & (gamma_values == 1.0))
        M_values_filtered = M_values[filtered_indices]
        lambda1_filtered = lambda1_values[filtered_indices]
        lambda2_filtered = lambda2_values[filtered_indices]
        lambda3_filtered = lambda3_values[filtered_indices]

        p_values, r_values = [], []
        for l1, l2, l3 in zip(lambda1_filtered, lambda2_filtered, lambda3_filtered):
            p, r = ultimatumGame_pq_frequency_continual(l1, l2, l3, selection_intensity)
            p_values.append(p)
            r_values.append(r)

        # Accumulate values
        accumulated_p_RR += p_values
        accumulated_r_RR += r_values

        p_RR_all.append(p_values)
        r_RR_all.append(r_values)

        # SF Data
        file_name = f'data_M_vs_lambda1_2_3/SC_SF_DB_50_degree{deg}_{index}_gamma1.0_M_variation.txt'
        data = np.loadtxt(file_name)
        gamma_values, u_values, M_values, lambda1_values, lambda2_values, lambda3_values = data.T

        filtered_indices = np.where((u_values == mutation_rate) & (gamma_values == gamma))
        M_values_filtered_SF = M_values[filtered_indices]
        lambda1_filtered = lambda1_values[filtered_indices]
        lambda2_filtered = lambda2_values[filtered_indices]
        lambda3_filtered = lambda3_values[filtered_indices]

        p_values_SF, r_values_SF = [], []
        for l1, l2, l3 in zip(lambda1_filtered, lambda2_filtered, lambda3_filtered):
            p, r = ultimatumGame_pq_frequency_continual(l1, l2, l3, selection_intensity)
            p_values_SF.append(p)
            r_values_SF.append(r)

        accumulated_p_SF += p_values_SF
        accumulated_r_SF += r_values_SF

        p_SF_all.append(p_values_SF)
        r_SF_all.append(r_values_SF)

    # Average values
    accumulated_p_RR /= ite_number
    accumulated_p_SF /= ite_number
    accumulated_r_RR /= ite_number
    accumulated_r_SF /= ite_number

    p_RR_all = np.array(p_RR_all)
    r_RR_all = np.array(r_RR_all)
    p_SF_all = np.array(p_SF_all)
    r_SF_all = np.array(r_SF_all)


    # Error margins (for demonstration, adjust accordingly)
    error_margin_RR_p = np.std(p_RR_all, axis=0)  # Example for error calculation
    error_margin_SF_p = np.std(p_SF_all, axis=0)  # Example for error calculation
    error_margin_RR_q = np.std(r_RR_all, axis=0)  # Example for error calculation
    error_margin_SF_q = np.std(r_SF_all, axis=0)  # Example for error calculation

    M_values_filtered /= 50

    # Plotting means with error regions
    plt.fill_between(M_values_filtered, accumulated_p_RR - error_margin_RR_p, accumulated_p_RR + error_margin_RR_p,
                     color=colors[1], alpha=0.2)
    plt.fill_between(M_values_filtered, accumulated_r_RR - error_margin_RR_q, accumulated_r_RR + error_margin_RR_q,
                     color=colors[1], alpha=0.2)

    plt.fill_between(M_values_filtered, accumulated_p_SF - error_margin_SF_p, accumulated_p_SF + error_margin_SF_p,
                     color=colors[0], alpha=0.2)
    plt.fill_between(M_values_filtered, accumulated_r_SF - error_margin_SF_q, accumulated_r_SF + error_margin_SF_q,
                     color=colors[0], alpha=0.2)

    # Solid lines for mean values
    plt.plot(M_values_filtered, accumulated_p_RR, label=rf'RR, $\langle p \rangle$', color=colors[1], linewidth=2.0)
    plt.plot(M_values_filtered, accumulated_r_RR, label=rf'RR, $\langle q \rangle$', color=colors[1], linewidth=2.0, linestyle='--')
    plt.plot(M_values_filtered, accumulated_p_SF, label=rf'SF, $\langle p \rangle$', color=colors[0], linewidth=2.0)
    plt.plot(M_values_filtered, accumulated_r_SF, label=rf'SF, $\langle q \rangle$', color=colors[0], linewidth=2.0, linestyle='--')


    # rect1 = patches.Rectangle((69,0.5027), 10, 0.0025, linewidth=2, edgecolor='black', facecolor='none', linestyle='-')
    # # 添加矩形框到图上
    # plt.gca().add_patch(rect1)

    # Customizing the plot
    plt.xlim(2, 4)
    # plt.ylim(0.490-0.001, 0.505+0.001)
    # plt.yticks([0.49, 0.495, 0.50, 0.505])

    # plt.xlim(3.8, 4)
    plt.ylim(0.5116, 0.5143)

    # plt.xlabel(r'$\overline{M}$', fontsize = 14)
    # plt.ylabel('Mean value', fontsize = 14)
    # plt.title(fr'$u={mutation_rate}, \gamma={gamma}$')
    # plt.legend(loc='best', fontsize = 12)

    ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax.ticklabel_format(style='scientific', axis='y', scilimits=(0, 0))

    # plt.savefig(rf"subFig_sixGames_M_ultimatumGame_pqFrequency_SF_RR_u{mutation_rate}_gamma{gamma}_delta{selection_intensity}_line_newColor.svg")
    plt.show()

def plot_M_donationGame(gamma, mut, degree):
    plt.rcParams['font.family'] = 'Arial'
    plt.rc('font', size=12.5)

    fig, ax = plt.subplots(figsize=(3.25, 3.25))

    accumulated_RR = np.zeros((100, 151))
    accumulated_SF = np.zeros((100, 151))
    M_values = None

    ite_number = 100
    for index in range(1, ite_number+1):
        filenameSF = f'data_M_vs_lambda1_2_3/SC_SF_DB_50_degree{degree}_{index}_gamma1.0_M_variation.txt'
        data_matrixSF = np.loadtxt(filenameSF)
        rows = data_matrixSF[np.isclose(data_matrixSF[:, 0], gamma) & np.isclose(data_matrixSF[:, 1], mut)]
        M_values = rows[:, 2]
        lambda1_values = rows[:, 3]
        lambda2_values = rows[:, 4]
        lambda3_values = rows[:, 5]
        bc_values = (lambda1_values + lambda2_values + lambda3_values) / (lambda1_values- lambda2_values)
        accumulated_SF[index-1, :] = bc_values

        filenameRR = f'data_M_vs_lambda1_2_3/SC_RR_DB_50_degree{degree}_{index}_M_variation.txt'
        data_matrixRR = np.loadtxt(filenameRR)
        rows = data_matrixRR[np.isclose(data_matrixRR[:, 0], 1.0) & np.isclose(data_matrixRR[:, 1], mut)]
        lambda1_values = rows[:, 3]
        lambda2_values = rows[:, 4]
        lambda3_values = rows[:, 5]
        bc_values = (lambda1_values + lambda2_values + lambda3_values) / (lambda1_values- lambda2_values)
        accumulated_RR[index-1, :] = bc_values

    
    # Calculate the mean and standard deviation for error bars
    mean_SF = np.mean(accumulated_SF, axis=0)
    std_SF = np.std(accumulated_SF, axis=0)
    mean_RR = np.mean(accumulated_RR, axis=0)
    std_RR = np.std(accumulated_RR, axis=0)

    M_values /= 50

    # Plot the mean lines with shaded error bars
    plt.plot(M_values, mean_SF, color='#c1272c', label='SF', linewidth=2.0)
    plt.fill_between(M_values, mean_SF - std_SF, mean_SF + std_SF, color='#c1272c', alpha=0.2)

    plt.plot(M_values, mean_RR, color='#0070bb', label='RR', linewidth=2.0)
    plt.fill_between(M_values, mean_RR - std_RR, mean_RR + std_RR, color='#0070bb', alpha=0.2)

    

    # Add horizontal lines for the initial value
    #plt.axhline(y=mean_SF[0], linestyle='--', color='#c1272c')
    #plt.axhline(y=mean_RR[0], linestyle='--', color='#0070bb')
    #plt.axhline(y=mean_SF[0], xmin=0.5/21, xmax=0.5/21+20/21, linestyle='--', color='#c1272c')
    #plt.axhline(y=mean_RR[0], xmin=0.5/21, xmax=0.5/21+20/21, linestyle='--', color='#0070bb')

    # Add labels and title
    plt.xlabel(r'$\overline{M}$', fontsize = 14)
    plt.ylabel(r'Critical ratio, $(b/c)^{\ast}$', fontsize = 14)
    plt.title(rf'$\bar{{d}} = {degree} , \gamma={gamma}, \mu = {mut}$')
    plt.legend()
    plt.xlim(1,4)
    plt.ylim(1.9, 6.1)
    plt.yticks([2,3,4,5,6])

    # plt.savefig(f'Fig_sixGames_M_bc_SF_RR_deg{degree}_gamma{gamma}_mutation{mut}_1000graph_line_newColor.svg')
    plt.show()



def plot_M_stagHuntGame(gamma, mut, degree, sigma, b, delta=0.001):
    accumulated_RR = np.zeros((100, 151))
    accumulated_SF = np.zeros((100, 151))
    M_values = None

    # Create a new figure
    fig, ax = plt.subplots(figsize=(3.25, 3.25))

    for index in range(1,101):
        filenameSF = f'data_M_vs_lambda1_2_3/SC_SF_DB_50_degree{degree}_{index}_gamma1.0_M_variation.txt'
        data_matrixSF = np.loadtxt(filenameSF)
        rows = data_matrixSF[np.isclose(data_matrixSF[:, 0], gamma) & np.isclose(data_matrixSF[:, 1], mut)]
        M_values = rows[:, 2]
        lambda1_values = rows[:, 3]
        lambda2_values = rows[:, 4]
        lambda3_values = rows[:, 5]
        K1_values = lambda1_values + lambda2_values + lambda3_values
        K2_values = lambda1_values - lambda2_values
        #plt.scatter(M_values, bc_values, marker='o', s=10, color='#c1272c', alpha=0.08)
        R1 = (K1_values + K2_values)/4
        R2 = (K1_values - K2_values)/4
        accumulated_SF[index-1, :] = 1/2 + delta /2 *b*( (sigma-1)*R1 -R2 )

        filenameRR = f'data_M_vs_lambda1_2_3/SC_RR_DB_50_degree{degree}_{index}_M_variation.txt'
        data_matrixRR = np.loadtxt(filenameRR)
        rows = data_matrixRR[np.isclose(data_matrixRR[:, 0], 1.0) & np.isclose(data_matrixRR[:, 1], mut)]
        lambda1_values = rows[:, 3]
        lambda2_values = rows[:, 4]
        lambda3_values = rows[:, 5]
        K1_values = lambda1_values + lambda2_values + lambda3_values
        K2_values = lambda1_values - lambda2_values
        R1 = (K1_values + K2_values)/4
        R2 = (K1_values - K2_values)/4
        #plt.scatter(M_values, bc_values, marker='o', s=10, color='#0070bb', alpha=0.08)
        accumulated_RR[index-1, :] = 1/2 + delta /2 *b*( (sigma-1)*R1 -R2 )

    # Calculate the mean and standard deviation for error bars
    mean_SF = np.mean(accumulated_SF, axis=0)
    std_SF = np.std(accumulated_SF, axis=0)
    mean_RR = np.mean(accumulated_RR, axis=0)
    std_RR = np.std(accumulated_RR, axis=0)

    M_values /= 50

    # Plot the mean lines with shaded error bars
    plt.plot(M_values, mean_SF, color='#c1272c', label='SF', linewidth=2.0)
    plt.fill_between(M_values, mean_SF - std_SF, mean_SF + std_SF, color='#c1272c', alpha=0.2)

    plt.plot(M_values, mean_RR, color='#0070bb', label='RR', linewidth=2.0)
    plt.fill_between(M_values, mean_RR - std_RR, mean_RR + std_RR, color='#0070bb', alpha=0.2)

    # Add horizontal lines for the initial value
    #plt.axhline(y=mean_SF[0], linestyle='--', color='#c1272c')
    #plt.axhline(y=mean_RR[0], linestyle='--', color='#0070bb')
    #plt.axhline(y=mean_SF[0], xmin=0.5/21, xmax=0.5/21+20/21, linestyle='--', color='#c1272c')
    #plt.axhline(y=mean_RR[0], xmin=0.5/21, xmax=0.5/21+20/21, linestyle='--', color='#0070bb')

    ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax.ticklabel_format(style='scientific', axis='y', scilimits=(0, 0))

    # Add labels and title
    plt.xlabel(r'$\overline{M}$', fontsize = 14)
    plt.ylabel('Mean value', fontsize = 14)
    plt.title(rf'$\gamma={gamma}, \mu = {mut}$')
    plt.legend()
    plt.xlim(1,4)
    plt.ylim(0.496, 0.502)
    # Show the figure
    # plt.savefig(f'Fig_sixGames_M_stagHuntGame_SF_RR_deg{degree}_gamma{gamma}_mutation{mut}_1000graph_newColor.svg')
    plt.show()


def GEG_pew_frequency_continual(lambda1, lambda2, lambda3, selection_intensity, b):
    phi = 1/2 + selection_intensity * (- lambda2 * 1/6 - lambda3 * 1/12)
    omega = 1/2 + selection_intensity * ((b-1)/12 *lambda1 - lambda2 * (b+1)/12 - lambda3 * 1/12 )
    return phi,omega

def plot_M_giftExchangeGame(deg, g, mutation_rate, gamma, selection_intensity=0.005):
    plt.rcParams['font.family'] = 'Arial'
    plt.rc('font', size=12.5)

    fig, ax = plt.subplots(figsize=(3.25, 3.25))

    colors = ['#c1272c', '#0070bb','orange','black']

    accumulated_phi_RR = np.zeros(151)
    accumulated_phi_SF = np.zeros(151)
    accumulated_omega_RR = np.zeros(151)
    accumulated_omega_SF = np.zeros(151)

    phi_SF_all = []
    omega_SF_all = []
    phi_RR_all = []
    omega_RR_all = []


    ite_number = 100
    for index in range(1, ite_number+1):
        # RR Data
        file_name = f'data_M_vs_lambda1_2_3/SC_RR_DB_50_degree{deg}_{index}_M_variation.txt'
        data = np.loadtxt(file_name)
        gamma_values, u_values, M_values, lambda1_values, lambda2_values, lambda3_values = data.T

        filtered_indices = np.where((u_values == mutation_rate) & (gamma_values == 1.0))
        M_values_filtered = M_values[filtered_indices]
        lambda1_filtered = lambda1_values[filtered_indices]
        lambda2_filtered = lambda2_values[filtered_indices]
        lambda3_filtered = lambda3_values[filtered_indices]

        phi_values = []
        omega_values = []
        for l1, l2, l3 in zip(lambda1_filtered, lambda2_filtered, lambda3_filtered):
            phi,omega = GEG_pew_frequency_continual(l1, l2, l3, selection_intensity, g)
            phi_values.append(phi)
            omega_values.append(omega)

        phi_RR_all.append(phi_values)
        omega_RR_all.append(omega_values)

        # Accumulate values
        accumulated_phi_RR += phi_values
        accumulated_omega_RR += omega_values

        # SF Data
        file_name = f'data_M_vs_lambda1_2_3/SC_SF_DB_50_degree{deg}_{index}_gamma1.0_M_variation.txt'
        data = np.loadtxt(file_name)
        gamma_values, u_values, M_values, lambda1_values, lambda2_values, lambda3_values = data.T

        filtered_indices = np.where((u_values == mutation_rate) & (gamma_values == gamma))
        M_values_filtered_SF = M_values[filtered_indices]
        lambda1_filtered = lambda1_values[filtered_indices]
        lambda2_filtered = lambda2_values[filtered_indices]
        lambda3_filtered = lambda3_values[filtered_indices]

        phi_values_SF, omega_values_SF = [], []
        for l1, l2, l3 in zip(lambda1_filtered, lambda2_filtered, lambda3_filtered):
            phi,omega = GEG_pew_frequency_continual(l1, l2, l3, selection_intensity, g)
            phi_values_SF.append(phi)
            omega_values_SF.append(omega)

        accumulated_phi_SF += phi_values_SF
        accumulated_omega_SF += omega_values_SF

        phi_SF_all.append(phi_values_SF)
        omega_SF_all.append(omega_values_SF)

    # Average values
    accumulated_phi_RR /= ite_number
    accumulated_phi_SF /= ite_number
    accumulated_omega_RR /= ite_number
    accumulated_omega_SF /= ite_number

    phi_SF_all = np.array(phi_SF_all)
    omega_SF_all = np.array(omega_SF_all)
    phi_RR_all = np.array(phi_RR_all)
    omega_RR_all = np.array(omega_RR_all)

    # Error margins 
    error_margin_RR_phi = np.std(phi_RR_all, axis=0)  
    error_margin_SF_phi = np.std(phi_SF_all, axis=0)  
    error_margin_RR_omega = np.std(omega_RR_all, axis=0)  
    error_margin_SF_omega = np.std(omega_SF_all, axis=0)  

    #print(f"error_SF_phi: {error_margin_SF_phi}")

    M_values_filtered /= 50

    # Plotting means with error regions
    plt.fill_between(M_values_filtered, accumulated_phi_RR - error_margin_RR_phi, accumulated_phi_RR + error_margin_RR_phi,
                     color=colors[1], alpha=0.2)
    plt.fill_between(M_values_filtered, accumulated_omega_RR - error_margin_RR_omega, accumulated_omega_RR + error_margin_RR_omega,
                     color=colors[1], alpha=0.2)

    plt.fill_between(M_values_filtered, accumulated_phi_SF - error_margin_SF_phi, accumulated_phi_SF + error_margin_SF_phi,
                     color=colors[0], alpha=0.2)
    plt.fill_between(M_values_filtered, accumulated_omega_SF - error_margin_SF_omega, accumulated_omega_SF + error_margin_SF_omega,
                     color=colors[0], alpha=0.2)

    # Solid lines for mean values
    plt.plot(M_values_filtered, accumulated_phi_RR, label=rf'RR, $\langle \varphi \rangle$', color=colors[1], linewidth=2.0)
    plt.plot(M_values_filtered, accumulated_omega_RR, label=rf'RR, $\langle \omega \rangle$', color=colors[1], linewidth=2.0, linestyle='--')
    plt.plot(M_values_filtered, accumulated_phi_SF, label=rf'SF, $\langle \varphi \rangle$', color=colors[0], linewidth=2.0)
    plt.plot(M_values_filtered, accumulated_omega_SF, label=rf'SF, $\langle \omega \rangle$', color=colors[0], linewidth=2.0, linestyle='--')

    # Customizing the plot
    plt.xlim(1, 4)
    # plt.xticks([1, 1.5, 2])
    plt.ylim(0.455, 0.497)

    plt.xlabel(r'$\overline{M}$', fontsize = 14)
    plt.ylabel('Mean value', fontsize = 14)
    plt.title(fr'$u={mutation_rate}$', fontsize = 14)
    plt.legend(loc='best', fontsize = 12)

    ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax.ticklabel_format(style='scientific', axis='y', scilimits=(0, 0))

    # plt.savefig(rf"Fig_sixGames_M_GiftExchangeGame_prFrequency_vs_SF_RR_g{g}_u{mutation_rate}_gamma{gamma}_delta{selection_intensity}_line_newColor.svg")
    plt.show()


def dictatorGame_pew_frequency_continual(lambda1, lambda2, lambda3, selection_intensity):
    p = 1/2 + selection_intensity * (- lambda2 * 1/6 - lambda3 * 1/12)
    return p
def plot_M_dictatorGame(deg, mutation_rate, gamma, selection_intensity=0.005):
    plt.rcParams['font.family'] = 'Arial'
    plt.rc('font', size=12.5)

    fig, ax = plt.subplots(figsize=(3.25, 3.25))

    colors = ['#c1272c', '#0070bb','orange','black']

    accumulated_phi_RR = np.zeros(151)
    accumulated_phi_SF = np.zeros(151)

    phi_SF_all = []
    phi_RR_all = []

    ite_number = 100
    for index in range(1, ite_number+1):
        # RR Data
        file_name = f'data_M_vs_lambda1_2_3/SC_RR_DB_50_degree{deg}_{index}_M_variation.txt'
        data = np.loadtxt(file_name)
        gamma_values, u_values, M_values, lambda1_values, lambda2_values, lambda3_values = data.T

        filtered_indices = np.where((u_values == mutation_rate) & (gamma_values == 1.0))
        M_values_filtered = M_values[filtered_indices]
        lambda1_filtered = lambda1_values[filtered_indices]
        lambda2_filtered = lambda2_values[filtered_indices]
        lambda3_filtered = lambda3_values[filtered_indices]

        phi_values = []
        for l1, l2, l3 in zip(lambda1_filtered, lambda2_filtered, lambda3_filtered):
            phi = dictatorGame_pew_frequency_continual(l1, l2, l3, selection_intensity)
            phi_values.append(phi)

        # Accumulate values
        accumulated_phi_RR += phi_values
        phi_RR_all.append(phi_values)

        # SF Data
        file_name = f'data_M_vs_lambda1_2_3/SC_SF_DB_50_degree{deg}_{index}_gamma1.0_M_variation.txt'
        data = np.loadtxt(file_name)
        gamma_values, u_values, M_values, lambda1_values, lambda2_values, lambda3_values = data.T

        filtered_indices = np.where((u_values == mutation_rate) & (gamma_values == gamma))
        M_values_filtered_SF = M_values[filtered_indices]
        lambda1_filtered = lambda1_values[filtered_indices]
        lambda2_filtered = lambda2_values[filtered_indices]
        lambda3_filtered = lambda3_values[filtered_indices]

        phi_values_SF= []
        for l1, l2, l3 in zip(lambda1_filtered, lambda2_filtered, lambda3_filtered):
            phi = dictatorGame_pew_frequency_continual(l1, l2, l3, selection_intensity)
            phi_values_SF.append(phi)

        accumulated_phi_SF += phi_values_SF
        phi_SF_all.append(phi_values_SF)

    # Average values
    accumulated_phi_RR /= ite_number
    accumulated_phi_SF /= ite_number


    phi_SF_all = np.array(phi_SF_all)
    phi_RR_all = np.array(phi_RR_all)

    # Error margins
    error_margin_RR_phi = np.std(phi_RR_all, axis=0)
    error_margin_SF_phi = np.std(phi_SF_all, axis=0)

    # print(f"error_SF_phi: {error_margin_SF_phi}")

    M_values_filtered /= 50

    # Plotting means with error regions
    plt.fill_between(M_values_filtered, accumulated_phi_RR - error_margin_RR_phi, accumulated_phi_RR + error_margin_RR_phi,
                     color=colors[1], alpha=0.2)

    plt.fill_between(M_values_filtered, accumulated_phi_SF - error_margin_SF_phi, accumulated_phi_SF + error_margin_SF_phi,
                     color=colors[0], alpha=0.2)

    # Solid lines for mean values
    plt.plot(M_values_filtered, accumulated_phi_RR, label=rf'RR', color=colors[1], linewidth=2.0)
    plt.plot(M_values_filtered, accumulated_phi_SF, label=rf'SF', color=colors[0], linewidth=2.0)

    # Customizing the plot
    plt.xlim(1, 4)
    # plt.xticks([1, 1.5, 2])
    plt.ylim(0.455, 0.485)

    # plt.yticks([0.47, 0.475, 0.48])

    plt.xlabel(r'$\overline{M}$', fontsize = 14)
    plt.ylabel('Mean value', fontsize = 14)
    plt.title(fr'$u={mutation_rate}$', fontsize = 14)
    plt.legend(loc='best', fontsize = 12)

    ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax.ticklabel_format(style='scientific', axis='y', scilimits=(0, 0))

    # plt.savefig(rf"Fig_sixGames_M_dictatorGame_prFrequency_SF_RR_u{mutation_rate}_gamma{gamma}_delta{selection_intensity}_line_newColor.svg")
    plt.show()


mutation_rate = 0.01
mutation_rate_UG = 0.01
gamma = 1.0
selection_intensity = 0.005
degree = 4

plot_M_donationGame(gamma, mutation_rate, degree)

plot_M_trustGame(degree, 1.5, mutation_rate, gamma, selection_intensity)

plot_M_ultimatumGame(degree, mutation_rate_UG, gamma, selection_intensity)
plot_M_ultimatumGame_subfig(degree, mutation_rate_UG, gamma, selection_intensity)

plot_M_stagHuntGame(gamma, mutation_rate, degree, 1.5, 1)

plot_M_giftExchangeGame(degree, 2, mutation_rate, gamma, selection_intensity)

plot_M_dictatorGame(degree, mutation_rate, gamma, selection_intensity)