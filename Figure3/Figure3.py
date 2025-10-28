import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import matplotlib.patches as patches

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
def donationGame_bcCritical(lambda1, lambda2, lambda3):
    return (lambda1 + lambda2 + lambda3) / (lambda1-lambda2)
def stagHuntGame_frequency(lambda1, lambda2, lambda3, delta, b, sigma):
    k1 = lambda1 + lambda2 + lambda3
    k2 = lambda1 - lambda2
    R1 = (k1 + k2) / 4
    R2 = (k1 - k2) / 4
    return 1/2 + delta / 2 * b * ((sigma - 1) * R1 - R2)

def get_lambdaValues(ite_number, mutation_rate):
    lambda_matrix = np.zeros((1959, 5))
    for index in range(1, ite_number + 1):
        filename = f"data_gammaVariance/StructureCoefficient_SF_DB_50_degree4_1_M100_gamma_variation_repeatation{index}.txt"
        data = np.loadtxt(filename)
        lambda_matrix += data
    lambda_matrix /= ite_number
    gamma_values, u_values, lambda1_values, lambda2_values, lambda3_values = lambda_matrix.T
    tol = 1e-12                       
    mask = np.isclose(u_values, mutation_rate, atol=tol)
    filtered_indices = np.where(mask)[0]

    gamma_values_filtered = gamma_values[filtered_indices]
    lambda1_filtered = lambda1_values[filtered_indices]
    lambda2_filtered = lambda2_values[filtered_indices]
    lambda3_filtered = lambda3_values[filtered_indices]

    return gamma_values_filtered, lambda1_filtered, lambda2_filtered, lambda3_filtered


def plot_donationGame_gamma(u, degree, index, ite_number=80):
    plt.rcParams['font.family'] = 'Arial'
    plt.rc('font', size=12.5)

    fig, ax = plt.subplots(figsize=(3.25, 3.25))

    if u == 0.01:
        lambda1_rigid = 65.6056489063
        lambda2_rigid = 44.5810555190 
        lambda3_rigid = 14.5993692353
    elif u == 0.03:
        lambda1_rigid = 35.8337674633
        lambda2_rigid = 22.7767219359
        lambda3_rigid = 21.9177730571
    elif u == 0.1:
        lambda1_rigid = 10.6739398911
        lambda2_rigid = 5.8233208336
        lambda3_rigid = 17.6548683620

    hline = donationGame_bcCritical(lambda1_rigid, lambda2_rigid, lambda3_rigid)

    plt.axhline(y=hline, linestyle='-', color='#c1272c', linewidth = 2.0, label=r'$\overline{M}=1$')

    gamma_values_filtered, lambda1_filtered, lambda2_filtered, lambda3_filtered = get_lambdaValues(ite_number, u)
    y_values = donationGame_bcCritical(lambda1_filtered, lambda2_filtered, lambda3_filtered)
    plt.plot(gamma_values_filtered, y_values, linestyle='-', linewidth = 2.0, color='#0070bb', label=r'$\overline{M}=2$')

    x_val = 1.0
    plt.axvline(x=x_val, linewidth = 1.5, color='grey', linestyle='--')
    y_pos = 4.96  
    plt.text(x_val - 1, y_pos, rf'$\gamma$ = {x_val}', color='black',
            bbox=dict(facecolor='white', edgecolor='none', pad=1.0))

    plt.xlabel(r'$\gamma$', fontsize = 13)
    plt.ylabel(r'Critical ratio, $(b/c)^{\ast}$', fontsize = 13)
    plt.title('Donation Game', fontsize = 13)
    # plt.legend(loc = 'lower left')

    plt.xlim(-10,10)
    plt.ylim([2.7, 7.3])
    # plt.savefig(f'Fig_sixGames_gamma_donationGame_u{u}_degree{degree}_index{index}_line.svg')

    plt.show()

def plot_trustGame_gamma(u, degree,b, ite_number=80, selection_intensity=0.005):
    plt.rcParams['font.family'] = 'Arial'
    plt.rc('font', size=12.5)

    fig, ax = plt.subplots(figsize=(3.25, 3.25))
    mutation_rate = u

    colors = ['#c1272c', '#0070bb']  

    if u == 0.01:
        lambda1_rigid = 65.6056489063
        lambda2_rigid = 44.5810555190 
        lambda3_rigid = 14.5993692353
    elif u == 0.03:
        lambda1_rigid = 35.8337674633
        lambda2_rigid = 22.7767219359
        lambda3_rigid = 21.9177730571
    elif u == 0.1:
        lambda1_rigid = 10.6739398911
        lambda2_rigid = 5.8233208336
        lambda3_rigid = 17.6548683620

    hline_p, hline_r = trustGame_pr_frequency_continual(lambda1_rigid, lambda2_rigid, lambda3_rigid, selection_intensity, b)    
    plt.axhline(y=hline_p, linestyle='-', color=colors[0], linewidth = 2.0, label=r'$\overline{M}=1$, $\langle p \rangle$')
    plt.axhline(y=hline_r, linestyle='--', color=colors[0], linewidth = 2.0, label=r'$\overline{M}=1$, $\langle r \rangle$')

    gamma_values_filtered, lambda1_filtered, lambda2_filtered, lambda3_filtered = get_lambdaValues(ite_number, mutation_rate)
    p_values, r_values = trustGame_pr_frequency_continual(lambda1_filtered, lambda2_filtered, lambda3_filtered, selection_intensity, b)
    plt.plot(gamma_values_filtered, p_values, linestyle='-', linewidth = 2.0, color=colors[1], label=r'$\overline{M}=2$, $\langle p \rangle$')
    plt.plot(gamma_values_filtered, r_values, linestyle='--', linewidth = 2.0, color=colors[1], label=r'$\overline{M}=2$, $\langle r \rangle$')


    plt.xlim(-10,10)
    # plt.ylim(0.469,0.499+0.002)

    x_val = 1.0
    plt.axvline(x=x_val, linewidth = 1.5, color='grey', linestyle='--')
    # y_pos = 4.96  
    # plt.text(x_val - 1, y_pos, rf'$\gamma$ = {x_val}', color='black',
    #         bbox=dict(facecolor='white', edgecolor='none', pad=1.0))

    ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax.ticklabel_format(style='scientific', axis='y', scilimits=(0, 0))


    plt.xlabel(r'$\gamma$', fontsize = 13)
    plt.ylabel('Mean value', fontsize = 13)
    plt.title('Trust Game', fontsize = 13)
    # plt.legend()
    ##plt.grid(True)
    # plt.savefig(rf"Fig_sixGames_gamma_TG_degree{degree}_b{b}_u{mutation_rate}_delta{selection_intensity}.svg")
    plt.show()


def plot_ultimatumGame_gamma(u, degree, ite_number=80, selection_intensity=0.005):
    plt.rcParams['font.family'] = 'Arial'
    plt.rc('font', size=12.5)

    mutation_rate = u
    fig, ax = plt.subplots(figsize=(3.25, 3.25))
    colors = ['#c1272c', '#0070bb']  

    if u == 0.01:
        lambda1_rigid = 65.6056489063
        lambda2_rigid = 44.5810555190 
        lambda3_rigid = 14.5993692353
    elif u == 0.03:
        lambda1_rigid = 35.8337674633
        lambda2_rigid = 22.7767219359
        lambda3_rigid = 21.9177730571
    elif u == 0.1:
        lambda1_rigid = 10.6739398911
        lambda2_rigid = 5.8233208336
        lambda3_rigid = 17.6548683620

    hline_p, hline_q = ultimatumGame_pq_frequency_continual(lambda1_rigid, lambda2_rigid, lambda3_rigid, selection_intensity)
    plt.axhline(y=hline_p, linestyle='-', color=colors[0], linewidth = 2.0, label=r'$\overline{M}=1$, $\langle p \rangle$')
    plt.axhline(y=hline_q, linestyle='--', color=colors[0], linewidth = 2.0, label=r'$\overline{M}=1$, $\langle q \rangle$')
    gamma_values_filtered, lambda1_filtered, lambda2_filtered, lambda3_filtered = get_lambdaValues(ite_number, mutation_rate)
    p_values, q_values = ultimatumGame_pq_frequency_continual(lambda1_filtered, lambda2_filtered, lambda3_filtered, selection_intensity)
    plt.plot(gamma_values_filtered, p_values, linestyle='-', linewidth = 2.0, color=colors[1], label=r'$\overline{M}=2$, $\langle p \rangle$')
    plt.plot(gamma_values_filtered, q_values, linestyle='--', linewidth = 2.0, color=colors[1], label=r'$\overline{M}=2$, $\langle q \rangle$')


    plt.xlim(-10,10)
    plt.ylim(0.46,0.52)

    x_val = 1.0
    plt.axvline(x=x_val, linewidth = 1.5, color='grey', linestyle='--')


    plt.xlabel(r'$\gamma$', fontsize = 13)
    plt.ylabel('Mean value', fontsize = 13)
    plt.title('Ultimatum Game', fontsize = 13)

    ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax.ticklabel_format(style='scientific', axis='y', scilimits=(0, 0))

    #plt.legend()

    #plt.grid(True)
    # plt.savefig(rf"Fig_sixGames_gamma_UG_degree{degree}_u{mutation_rate}_delta{selection_intensity}.svg")
    plt.show()

def plot_stagHuntGame_gamma(u, degree, index, sigma, b, ite_number=80,delta=0.001):
    plt.rcParams['font.family'] = 'Arial'
    plt.rc('font', size=12.5)

    fig, ax = plt.subplots(figsize=(3.25, 3.25))

    if u == 0.01:
        lambda1_rigid = 65.6056489063
        lambda2_rigid = 44.5810555190 
        lambda3_rigid = 14.5993692353
    elif u == 0.03:
        lambda1_rigid = 35.8337674633
        lambda2_rigid = 22.7767219359
        lambda3_rigid = 21.9177730571
    elif u == 0.1:
        lambda1_rigid = 10.6739398911
        lambda2_rigid = 5.8233208336
        lambda3_rigid = 17.6548683620

    hline = stagHuntGame_frequency(lambda1_rigid, lambda2_rigid, lambda3_rigid, delta, b, sigma)

    plt.axhline(y=hline, linestyle='-', color='#c1272c', linewidth = 2.0, label=r'$\overline{M}=1$')

    gamma_values_filtered, lambda1_filtered, lambda2_filtered, lambda3_filtered = get_lambdaValues(ite_number, u)
    y_values = stagHuntGame_frequency(lambda1_filtered, lambda2_filtered, lambda3_filtered, delta, b, sigma)
    plt.plot(gamma_values_filtered, y_values, linestyle='-', linewidth = 2.0, color='#0070bb', label=r'$\overline{M}=2$')

    x_val = 1.0
    plt.axvline(x=x_val, linewidth = 1.5, color='grey', linestyle='--')
    # y_pos = 0.496  
    # plt.text(x_val - 1, y_pos, rf'$\gamma$ = {x_val}', color='black',
    #         bbox=dict(facecolor='white', edgecolor='none', pad=1.0))

    ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax.ticklabel_format(style='scientific', axis='y', scilimits=(0, 0))

    plt.xlabel(r'$\gamma$', fontsize = 13)
    plt.ylabel('Mean value', fontsize = 13)
    plt.title('Stag-Hunt Game', fontsize = 13)
    # plt.legend(loc = 'lower left')

    plt.xlim(-10,10)
    # plt.ylim([0.4935, 0.5005])
    # plt.yticks([0.494, 0.496, 0.498, 0.500])

    # plt.savefig(f'Fig_sixGames_gamma_StaghuntGame_u{u}_degree{degree}_index{index}_line.svg')
    plt.show()

def plot_giftExchangeGame_gamma(mutation_rate, degree, g, ite_number=80, selection_intensity=0.005):
    plt.rcParams['font.family'] = 'Arial'
    plt.rc('font', size=12.5)

    fig, ax = plt.subplots(figsize=(3.25, 3.25))
    colors = ['#c1272c', '#0070bb']  
    u = mutation_rate

    if u == 0.01:
        lambda1_rigid = 65.6056489063
        lambda2_rigid = 44.5810555190 
        lambda3_rigid = 14.5993692353
    elif u == 0.03:
        lambda1_rigid = 35.8337674633
        lambda2_rigid = 22.7767219359
        lambda3_rigid = 21.9177730571
    elif u == 0.1:
        lambda1_rigid = 10.6739398911
        lambda2_rigid = 5.8233208336
        lambda3_rigid = 17.6548683620

    phi_rigid, omega_rigid = GEG_pew_frequency_continual(lambda1_rigid, lambda2_rigid, lambda3_rigid, selection_intensity, g)
    plt.axhline(y=phi_rigid, linestyle='-', color=colors[0], linewidth = 2.0, label=r'$\overline{M}=1$, $\langle \phi \rangle$')
    plt.axhline(y=omega_rigid, linestyle='--', color=colors[0], linewidth = 2.0, label=r'$\overline{M}=1$, $\langle \omega \rangle$')
    gamma_values_filtered, lambda1_filtered, lambda2_filtered, lambda3_filtered = get_lambdaValues(ite_number, u)
    phi_values, omega_values = GEG_pew_frequency_continual(lambda1_filtered, lambda2_filtered, lambda3_filtered, selection_intensity, g)
    plt.plot(gamma_values_filtered, phi_values, linestyle='-', linewidth = 2.0, color=colors[1], label=r'$\overline{M}=2$, $\langle \phi \rangle$')
    plt.plot(gamma_values_filtered, omega_values, linestyle='--', linewidth = 2.0, color=colors[1], label=r'$\overline{M}=2$, $\langle \omega \rangle$')

    plt.xlim(-10,10)
    plt.ylim(0.44,0.49)

    x_val = 1.0
    plt.axvline(x=x_val, linewidth = 1.5, color='grey', linestyle='--')

    ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax.ticklabel_format(style='scientific', axis='y', scilimits=(0, 0))

    plt.xlabel(r'$\gamma$', fontsize = 13)
    plt.ylabel('Mean value', fontsize = 13)
    plt.title('Gift-Exchange Game', fontsize = 13)
    # plt.legend(fontsize = 12)

    # plt.savefig(rf"Fig_sixGames_gamma_GEG_degree{degree}_g{g}_u{mutation_rate}_delta{selection_intensity}.svg")
    plt.show()

def plot_dictatorGame_gamma(mutation_rate, degree, ite_number=80, selection_intensity=0.005):
    plt.rcParams['font.family'] = 'Arial'
    plt.rc('font', size=12.5)

    fig, ax = plt.subplots(figsize=(3.25, 3.25))
    colors = ['#c1272c', '#0070bb']  

    u = mutation_rate

    if u == 0.01:
        lambda1_rigid = 65.6056489063
        lambda2_rigid = 44.5810555190 
        lambda3_rigid = 14.5993692353
    elif u == 0.03:
        lambda1_rigid = 35.8337674633
        lambda2_rigid = 22.7767219359
        lambda3_rigid = 21.9177730571
    elif u == 0.1:
        lambda1_rigid = 10.6739398911
        lambda2_rigid = 5.8233208336
        lambda3_rigid = 17.6548683620

    p_rigid = dictatorGame_pew_frequency_continual(lambda1_rigid, lambda2_rigid, lambda3_rigid, selection_intensity)
    plt.axhline(y=p_rigid, linestyle='-', color=colors[0], linewidth = 2.0, label=r'$\overline{M}=1$, $\langle p \rangle$')
    gamma_values_filtered, lambda1_filtered, lambda2_filtered, lambda3_filtered = get_lambdaValues(ite_number, u)  
    p_values = dictatorGame_pew_frequency_continual(lambda1_filtered, lambda2_filtered, lambda3_filtered, selection_intensity)
    plt.plot(gamma_values_filtered, p_values, linestyle='-', linewidth = 2.0, color=colors[1], label=r'$\overline{M}=2$, $\langle p \rangle$')

    plt.xlim(-10,10)
    plt.ylim(0.44, 0.48)

    x_val = 1.0
    plt.axvline(x=x_val, linewidth = 1.5, color='grey', linestyle='--')

    ax.yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax.ticklabel_format(style='scientific', axis='y', scilimits=(0, 0))


    plt.xlabel(r'$\gamma$', fontsize = 13)
    plt.ylabel('Mean value', fontsize = 13)
    plt.title('Dictator game', fontsize = 13)
    # plt.legend(fontsize = 12)

    # plt.savefig(rf"Fig_sixGames_gamma_dictatorGame_degree{degree}_u{mutation_rate}_delta{selection_intensity}.svg")
    plt.show()

selection_intensity = 0.005
ite_number = 80
mutation_rate = 0.01
degree = 4
net_index = 1

plot_donationGame_gamma(mutation_rate, degree, net_index, ite_number)

plot_trustGame_gamma(mutation_rate, degree, 1.5, ite_number, selection_intensity)

plot_stagHuntGame_gamma(mutation_rate, degree, net_index, 1, 1.5, ite_number, selection_intensity)

plot_ultimatumGame_gamma(mutation_rate, degree, ite_number, selection_intensity)

g=2
plot_giftExchangeGame_gamma(mutation_rate, degree, g, ite_number, selection_intensity)
plot_dictatorGame_gamma(mutation_rate, degree, ite_number, selection_intensity)
