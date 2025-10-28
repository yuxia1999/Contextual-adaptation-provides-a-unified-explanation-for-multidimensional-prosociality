import numpy as np
import matplotlib.pyplot as plt
from matplotlib.transforms import Bbox
import matplotlib.patches as patches
from matplotlib.ticker import ScalarFormatter

def trustGame_pr_frequency_continual(lambda1, lambda2, lambda3, selection_intensity, b):
    p = 1/2 + selection_intensity * ( lambda1 * 1/12 * (b-1)     - lambda2 * 1/12 + lambda3 * (b-2)/24)
    r = 1/2 + selection_intensity * ( lambda2 * (-b/12)  - lambda3*(b/24) )
    return p,r
def ultimatumGame_pq_frequency_continual(lambda1, lambda2, lambda3, selection_intensity):
    p = 1/2 + selection_intensity * ( lambda1*1/12     - lambda2*1/12)
    q = 1/2 + selection_intensity * ( lambda1*(-1/12)  - lambda3*(1/24) )
    return p,q
def GEG_pew_frequency_continual(lambda1, lambda2, lambda3, selection_intensity, b):
    phi = 1/2 + selection_intensity * (- lambda2 * 1/6 - lambda3 * 1/12)
    omega = 1/2 + selection_intensity * ((b-1)/12 *lambda1 - lambda2 * (b+1)/12 - lambda3 * 1/12 )
    return phi,omega
def dictatorGame_pew_frequency_continual(lambda1, lambda2, lambda3, selection_intensity):
    p = 1/2 + selection_intensity * (- lambda2 * 1/6 - lambda3 * 1/12)
    return p

def DG_u_bc_bigdata_line(index, M, degree):

    plt.rcParams['font.family'] = 'Arial'
    plt.rc('font', size=12.5)

    plt.figure(figsize=(3.25, 3.25))

    gamma_1_positive_data = {'positive': [], 'm_values': []}
    gamma_1_negative_data = {'negative': [], 'm_values': []}
    gamma_minus_1_positive_data = {'positive': [], 'm_values': []}
    gamma_minus_1_negative_data = {'negative': [], 'm_values': []}
    gamma_0_positive_data = {'positive': [], 'm_values': []}
    gamma_0_negative_data = {'negative': [], 'm_values': []}

    for idx in range(1, 2):
        filenameSF1 = f'EDF9_data/DonationGame_SF_BD_bcCritical_50_gamma1.0_degree{degree}_{idx}_u_variation.txt'
        try:
            data_matrixSF = np.loadtxt(filenameSF1)
            rows = data_matrixSF[np.isclose(data_matrixSF[:, 2], M)]
            for m_value in np.unique(rows[:, 1]):
                m_rows = rows[rows[:, 1] == m_value]
                positive_values = m_rows[m_rows[:, 5] >= 0, 5]
                negative_values = m_rows[m_rows[:, 5] < 0, 5]
                if len(positive_values) > 0:
                    gamma_1_positive_data['positive'].append(np.mean(positive_values))
                    gamma_1_positive_data['m_values'].append(m_value)
                if len(negative_values) > 0:
                    gamma_1_negative_data['negative'].append(np.mean(negative_values))
                    gamma_1_negative_data['m_values'].append(m_value)
        except OSError:
            continue

        filenameSF2 = f'EDF9_data/DonationGame_SF_BD_bcCritical_50_M50_degree{degree}_{idx}_u_variation.txt'
        try:
            data_matrixSF = np.loadtxt(filenameSF2)
            rows = data_matrixSF[np.isclose(data_matrixSF[:, 2], 50)]
            for m_value in np.unique(rows[:, 1]):
                m_rows = rows[rows[:, 1] == m_value]
                positive_values = m_rows[m_rows[:, 5] >= 0, 5]
                negative_values = m_rows[m_rows[:, 5] < 0, 5]
                if len(positive_values) > 0:
                    gamma_minus_1_positive_data['positive'].append(np.mean(positive_values))
                    gamma_minus_1_positive_data['m_values'].append(m_value)
                if len(negative_values) > 0:
                    gamma_minus_1_negative_data['negative'].append(np.mean(negative_values))
                    gamma_minus_1_negative_data['m_values'].append(m_value)
        except OSError:
            continue

        filenameSF3 = f'EDF9_data/DonationGame_SF_BD_bcCritical_50_gamma0.0_degree{degree}_{idx}_u_variation.txt'
        try:
            data_matrixSF = np.loadtxt(filenameSF3)
            rows = data_matrixSF[np.isclose(data_matrixSF[:, 2], M)]
            for m_value in np.unique(rows[:, 1]):
                m_rows = rows[rows[:, 1] == m_value]
                positive_values = m_rows[m_rows[:, 5] >= 0, 5]
                negative_values = m_rows[m_rows[:, 5] < 0, 5]
                if len(positive_values) > 0:
                    gamma_0_positive_data['positive'].append(np.mean(positive_values))
                    gamma_0_positive_data['m_values'].append(m_value)
                if len(negative_values) > 0:
                    gamma_0_negative_data['negative'].append(np.mean(negative_values))
                    gamma_0_negative_data['m_values'].append(m_value)
        except OSError:
            continue

    def safe_plot(data, label, color):
        positive_values = data.get('positive', [])
        negative_values = data.get('negative', [])

        if len(data['m_values']) > 0:
            if len(positive_values) > 0:
                plt.plot(data['m_values'], positive_values, color=color,  linestyle='-', alpha=1, linewidth=2.0)
            if len(negative_values) > 0:
                plt.plot(data['m_values'], negative_values, color=color,  linestyle='-', alpha=1, label=label, linewidth=2.0)

    def plot_vertical_lines(positive_data, negative_data, color):
        if len(positive_data) > 0 and len(negative_data) > 0:
            max_positive_x = max(positive_data['m_values'])
            min_negative_x = min(negative_data['m_values'])
            vertical_x = (max_positive_x + min_negative_x) / 2
            plt.axvline(x=vertical_x, color=color, linestyle='--', linewidth = 1.5)

    def safe_plot_vertical_lines(data, label, color, negative_data=None):
        if len(data['m_values']) > 0:
            if negative_data and len(negative_data['m_values']) > 0:  
                plot_vertical_lines(data, negative_data, color) 

    # M = N
    safe_plot(gamma_minus_1_positive_data, 'Ungrouped', '#458c6e')
    safe_plot(gamma_minus_1_negative_data, r'$M = 1.0$', '#458c6e')
    safe_plot_vertical_lines(gamma_minus_1_positive_data, 'Ungrouped', '#458c6e',gamma_minus_1_negative_data)
    # gamma = 1.0
    safe_plot(gamma_1_positive_data, r'$\gamma = 1.0$', '#0070bb')
    safe_plot(gamma_1_negative_data, r'$\gamma = 1.0$', '#0070bb')
    safe_plot_vertical_lines(gamma_1_positive_data, 'Ungrouped', '#0070bb', gamma_1_negative_data)

    # gamma = 0.0
    safe_plot(gamma_0_positive_data, r'$\gamma = 0.0$', '#c1272c')
    safe_plot(gamma_0_negative_data, r'$\gamma = 0.0$', '#c1272c')
    safe_plot_vertical_lines(gamma_0_positive_data, 'Ungrouped','#c1272c', gamma_0_negative_data)

    rect1 = patches.Rectangle((0.0011,1e1), 0.01-0.0011, 1e2-1e1, linewidth=1.5, edgecolor='black', facecolor='none', linestyle='-')
    plt.gca().add_patch(rect1)

    plt.xlabel(r'Mutation rate, $\mu$', fontsize=14)
    plt.ylabel(r'$(b/c)^{\ast}$', fontsize=14)
    plt.title(rf'$\bar{{d}} = {degree} , M = {M}$', fontsize=14)
    plt.legend(loc='lower right', markerscale=0.01, frameon=True, fontsize = 11)
    plt.xlim(1e-3, 1)
    plt.xscale('log')
    plt.ylim(-3e3, 3e3)
    plt.axhline(y=0, color='grey', linestyle='--', linewidth = 1.5)
    plt.yscale('symlog', linthresh=10)

    # plt.savefig(f'ExtendedDataFig9_largeNetwork_BD_DG_u_vs_bcCritical_degree{degree}_1_line.svg')
    plt.show()    

def TG_pq_vs_u_for_degrees(file_pattern, deg_values, M_numbers, b, selection_intensity=0.001):
    plt.rcParams['font.family'] = 'Arial'
    plt.rc('font', size=12.5)

    plt.figure(figsize=(3.25, 3.25))

    colors = ['#c1272c',  '#0070bb']  

    for i, M_number in enumerate(M_numbers):
        for deg in deg_values:
            file_name = file_pattern.format(deg=deg, M_number=M_number)
            data = np.loadtxt(file_name)
            u_values, lambda1_values, lambda2_values, lambda3_values = data.T

            p_values = []
            r_values = []
            for l1, l2, l3 in zip(lambda1_values, lambda2_values, lambda3_values):
                p, r = trustGame_pr_frequency_continual(l1, l2, l3, selection_intensity,b)
                p_values.append(p)
                r_values.append(r)

            plt.plot(u_values, p_values, linewidth=1.8 , label=rf'$\overline{{M}}=${M_number/30:.2f}, $\langle p \rangle$', color=colors[i], linestyle='-')
            plt.plot(u_values, r_values, linewidth=1.8 , label=rf'$\overline{{M}}=${M_number/30:.2f}, $\langle r \rangle$', color=colors[i], linestyle='--')

    plt.xscale('log')
    plt.xlim(1e-3,1)
    if deg_values[0] >=16:
        plt.xlabel(r'Mutation rate, $\mu$', fontsize = 13.5)
    if deg_values[0] in [4,10,16]:
        plt.ylabel('Mean value', fontsize = 13.5)

    plt.yticks([0.488, 0.492, 0.496, 0.5])

    plt.gca().yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))

    plt.title(fr'$\bar{{d}} = {deg}$', fontsize = 13.5)
    plt.legend(loc='lower right', fontsize = 11)

    # plt.savefig(rf"ExtendedDataFig_largeNetwork_BD_TG_prFrequency_SF_degree{deg_values}_b{b}_gamma1.0.svg")
    plt.show()
 
def UG_pq_vs_u_for_degrees(file_pattern, deg_values, M_numbers, selection_intensity=0.001):
    plt.rcParams['font.family'] = 'Arial'
    plt.rc('font', size=12.5)

    plt.figure(figsize=(3.25, 3.25))
    colors = ['#c1272c',  '#0070bb']  

    for i, M_number in enumerate(M_numbers):
        for deg in deg_values:
            file_name = file_pattern.format(deg=deg, M_number=M_number)
            data = np.loadtxt(file_name)
            u_values, lambda1_values, lambda2_values, lambda3_values = data.T

            p_values = []
            q_values = []
            for l1, l2, l3 in zip(lambda1_values, lambda2_values, lambda3_values):
                p, q = ultimatumGame_pq_frequency_continual(l1, l2, l3, selection_intensity)
                p_values.append(p)
                q_values.append(q)

            plt.plot(u_values, p_values, linewidth=1.8 , label=rf'$\overline{{M}}=${M_number/30:.2f}, $\langle p \rangle$', color=colors[i], linestyle='-')
            plt.plot(u_values, q_values, linewidth=1.8 , label=rf'$\overline{{M}}=${M_number/30:.2f}, $\langle q \rangle$', color=colors[i], linestyle='--')

    plt.xlim(1e-3,1)
    plt.xscale('log')
    if deg_values[0] >=16:
        plt.xlabel(r'Mutation rate, $\mu$', fontsize = 13.5)
    if deg_values[0] in [4,10,16]:
        plt.ylabel('Mean value', fontsize = 13.5)

    plt.gca().yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))


    plt.title(fr'$\bar{{d}} = {deg}$', fontsize = 13.5)
    plt.legend(loc='lower right', fontsize = 11)

    # plt.savefig(rf"ExtendedDataFig9_largeNetwork_BD_ultimatumGame_pqFrequency_SF_degree{deg_values}_gamma1.0.svg")
    plt.show()

def GEG_pq_vs_u_for_degrees(file_pattern, deg_values, M_numbers, b, selection_intensity=0.001):
    plt.rcParams['font.family'] = 'Arial'
    plt.rc('font', size=12.5)

    plt.figure(figsize=(3.25, 3.25))

    colors = ['#c1272c',  '#0070bb'] 

    for i, M_number in enumerate(M_numbers):
        for deg in deg_values:
            file_name = file_pattern.format(deg=deg, M_number=M_number)
            data = np.loadtxt(file_name)
            u_values, lambda1_values, lambda2_values, lambda3_values = data.T

            p_values = []
            r_values = []
            for l1, l2, l3 in zip(lambda1_values, lambda2_values, lambda3_values):
                p, r = GEG_pew_frequency_continual(l1, l2, l3, selection_intensity, b)
                p_values.append(p)
                r_values.append(r)

            plt.plot(u_values, p_values, linewidth=1.8 , label=rf'$\overline{{M}}=${M_number/30:.2f}, $\langle \varphi \rangle$', color=colors[i], linestyle='-')
            plt.plot(u_values, r_values, linewidth=1.8 , label=rf'$\overline{{M}}=${M_number/30:.2f}, $\langle \omega \rangle$', color=colors[i], linestyle='--')

    plt.xscale('log')
    plt.xlim(1e-3,1)
    if deg_values[0] >=16:
        plt.xlabel(r'Mutation rate, $\mu$', fontsize = 13.5)
    if deg_values[0] in [4,10,16]:
        plt.ylabel('Mean value', fontsize = 13.5)

    plt.ylim(0.479, 0.501)

    plt.gca().yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))

    plt.title(fr'$\bar{{d}} = {deg}$', fontsize = 13.5)
    plt.legend(loc='lower right', fontsize = 11)

    # plt.savefig(rf"ExtendedDataFig9_largeNetwork_BD_GEG_prFrequency_SF_degree{deg_values}_b{b}_gamma1.0.svg")
    plt.show()


def DTG_pq_vs_u_for_degrees(file_pattern, deg_values, M_numbers, selection_intensity=0.001):
    plt.rcParams['font.family'] = 'Arial'
    plt.rc('font', size=12.5)

    plt.figure(figsize=(3.25, 3.25))
    colors = ['#c1272c',  '#0070bb']  

    for i, M_number in enumerate(M_numbers):
        for deg in deg_values:
            file_name = file_pattern.format(deg=deg, M_number=M_number)
            data = np.loadtxt(file_name)
            u_values, lambda1_values, lambda2_values, lambda3_values = data.T

            p_values = []
            for l1, l2, l3 in zip(lambda1_values, lambda2_values, lambda3_values):
                p = dictatorGame_pew_frequency_continual(l1, l2, l3, selection_intensity)
                p_values.append(p)

            plt.plot(u_values, p_values, linewidth=1.8 , label=rf'$\overline{{M}}=${M_number/30:.2f}, $\langle p \rangle$', color=colors[i], linestyle='-')

    plt.xlim(1e-3,1)
    plt.xscale('log')
    if deg_values[0] >=16:
        plt.xlabel(r'Mutation rate, $\mu$', fontsize = 13.5)
    if deg_values[0] in [4,10,16]:
        plt.ylabel('Mean value', fontsize = 13.5)
    plt.yticks([0.485, 0.490, 0.495, 0.5])
    plt.gca().yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))

    plt.title(fr'$\bar{{d}} = {deg}$', fontsize = 13.5)
    plt.legend(loc='lower right', fontsize = 11)
    # plt.savefig(rf"ExtendedDataFig9_largeNetwork_BD_DTG_pqFrequency_SF_degree{deg_values}_gamma1.0.svg")
    plt.show()

def SHG_pq_vs_u_for_degrees(degree, selection_intensity=0.001):
    plt.rcParams['font.family'] = 'Arial'
    plt.rc('font', size=12.5)

    plt.figure(figsize=(3.25, 3.25))
    colors = ['#458c6e', '#c1272c',  '#0070bb']  

    file_name = f'EDF9_data/DonationGame_SF_BD_bcCritical_50_M50_degree4_1_u_variation.txt'
    data = np.loadtxt(file_name)
    gamma_values, u_values, M_values, K1_values, K2_values, bc_values = data.T

    R1 = (K1_values + K2_values)/4
    R2 = (K1_values - K2_values)/4

    sigma = 1.5
    b = 1.0

    p_values = 1/2 + selection_intensity /2 *b*( (sigma-1)*R1 -R2 )

    M_number = M_values[0]
    plt.plot(u_values, p_values, linewidth=1.8 , label=rf'$\overline{{M}}=${M_number/50:.2f}', color=colors[0], linestyle='-')

    for gamma in [0.0, 1.0]:
        file_name = f'EDF9_data/DonationGame_SF_BD_bcCritical_50_gamma{gamma}_degree{degree}_1_u_variation.txt'
        data = np.loadtxt(file_name)
        gamma_values, u_values, M_values, K1_values, K2_values, bc_values = data.T

        R1 = (K1_values + K2_values)/4
        R2 = (K1_values - K2_values)/4

        sigma = 1.5
        b = 1.0

        p_values = 1/2 + selection_intensity /2 *b*( (sigma-1)*R1 -R2 )

        M_number = M_values[0]
        plt.plot(u_values, p_values, linewidth=1.8 , label=rf'$\overline{{M}}=${M_number/50:.2f}, $\gamma = {gamma}$', color=colors[int(gamma+1)], linestyle='-')

    plt.xlim(1e-3,1)
    plt.xscale('log')

    plt.gca().yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))

    plt.title(fr'$\bar{{d}} = {degree}$', fontsize = 13.5)
    plt.legend(loc='lower right', fontsize = 11)

    # plt.savefig(rf"ExtendedDataFig9_largeNetwork_BD_SHG_pqFrequency_SF_degree{degree}_gamma1.0.svg")
    plt.show()

DG_u_bc_bigdata_line(1, 100, 4)
TG_pq_vs_u_for_degrees("EDF9_data/MS_BD_SF_30_degree{deg}_1_M{M_number}_gamma1.0_structureCoefficient.txt", [4], [30,  50], b=1.5)
UG_pq_vs_u_for_degrees("EDF9_data/MS_BD_SF_30_degree{deg}_1_M{M_number}_gamma1.0_structureCoefficient.txt", [4], [30,  50])
SHG_pq_vs_u_for_degrees(4)
GEG_pq_vs_u_for_degrees("EDF9_data/MS_BD_SF_30_degree{deg}_1_M{M_number}_gamma1.0_structureCoefficient.txt", [4], [30,  50], b=5)
DTG_pq_vs_u_for_degrees("EDF9_data/MS_BD_SF_30_degree{deg}_1_M{M_number}_gamma1.0_structureCoefficient.txt", [4], [30,  50])
