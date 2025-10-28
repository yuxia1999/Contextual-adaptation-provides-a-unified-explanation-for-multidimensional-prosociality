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

def DG_bc_barchart_resort(p_ungrouped_file, p_grouped_file):
    plt.rcParams['font.family'] = 'Arial'
    plt.rc('font', size=12.5)

    with open(p_ungrouped_file, 'r') as f:
        p_ungrouped = np.array([float(line.strip()) for line in f.readlines()])

    with open(p_grouped_file, 'r') as f:
        p_grouped = np.array([float(line.strip()) for line in f.readlines()])

    ungrouped_color = '#F1BFB577'  
    grouped_color = '#788AB6'  
    ungrouped_edgecolor = '#ED3B0E77'  
    grouped_edgecolor = '#5B679B'  

    positive_indices = np.where(p_grouped >= 0)[0]
    negative_indices = np.where(p_grouped < 0)[0]

    sorted_pos_indices = positive_indices[np.argsort(p_grouped[positive_indices])]
    sorted_neg_indices = negative_indices[np.argsort(p_grouped[negative_indices])]

    sorted_indices = np.concatenate((sorted_pos_indices, sorted_neg_indices))

    p_ungrouped_sorted = p_ungrouped[sorted_indices]
    p_grouped_sorted = p_grouped[sorted_indices]

    bar_width = 0.99
    plt.rcParams['font.family'] = 'Arial'
    plt.rc('font', size=12.5)

    fig, ax = plt.subplots(figsize=(3.25, 3.25))

    plt.bar(np.arange(len(p_grouped_sorted)), p_grouped_sorted, width=bar_width, label=r'$\overline{{M}}=\bar{d}$', 
            color=grouped_color, alpha=0.4, edgecolor=grouped_edgecolor)
    plt.bar(np.arange(len(p_ungrouped_sorted)), p_ungrouped_sorted, width=bar_width, label=r'$\overline{{M}}=1.0$', 
            color=ungrouped_color, edgecolor=ungrouped_edgecolor)

    plt.xlabel('Graph', fontsize = 14)
    plt.ylabel(r'$(b/c)^{\ast}$', fontsize = 14)
    plt.title(r'$\mu = 0.1$', fontsize = 14)

    plt.xticks([0, 29, 59, 89, 111], ['1', '30', '60', '90','112'])
    plt.axvline(x=60.5, color='gray', linestyle='--', linewidth=0.8)
    plt.yscale('symlog', linthresh=10)

    plt.legend(fontsize = 12)
    plt.xlim(-0.5, 111.5)
    plt.ylim(-7e2,7e2)

    # plt.savefig("ExtendedDataFig_smallNetwork_DG_BD_bcCritical_u0.1.svg")
    plt.show()

def UG_pq_barchart(p_ungrouped_file, p_grouped_file, q_ungrouped_file, q_grouped_file):
    with open(p_ungrouped_file, 'r') as f:
        p_ungrouped = np.array([float(line.strip()) for line in f.readlines()])

    with open(p_grouped_file, 'r') as f:
        p_grouped = np.array([float(line.strip()) for line in f.readlines()])

    with open(q_ungrouped_file, 'r') as f:
        q_ungrouped = np.array([float(line.strip()) for line in f.readlines()])

    with open(q_grouped_file, 'r') as f:
        q_grouped = np.array([float(line.strip()) for line in f.readlines()])

    ungrouped_color = '#F1BFB577'
    grouped_color = '#788AB6'
    ungrouped_edgecolor = '#ED3B0E77'
    grouped_edgecolor = '#5B679B'

    sorted_indices = np.argsort(p_ungrouped)
    p_ungrouped_sorted = p_ungrouped[sorted_indices]
    p_grouped_sorted = p_grouped[sorted_indices]
    q_ungrouped_sorted = q_ungrouped[sorted_indices]
    q_grouped_sorted = q_grouped[sorted_indices]

    bar_width = 0.99

    plt.figure(figsize=(5.2/1.6, 2.5/1.6))
    plt.bar(np.arange(len(p_grouped_sorted)), p_grouped_sorted, width=bar_width, label=r'$\overline{{M}}=\bar{d}$', 
            color=grouped_color, alpha=0.4, edgecolor=grouped_edgecolor)
    plt.bar(np.arange(len(p_ungrouped_sorted)), p_ungrouped_sorted, width=bar_width, label=r'$\overline{{M}}=1.0$', 
            color=ungrouped_color, edgecolor=ungrouped_edgecolor)

    plt.xlabel('Graph', fontsize = 14)
    plt.ylabel(r'$\langle p \rangle$', fontsize = 14)
    plt.title(r'$\delta = 0.001, \mu = 0.1$', fontsize = 14)
    plt.xlim(-0.5, 111.5)

    plt.xticks([0, 29, 59, 89, 111], ['1', '30', '60', '90','112'])

    plt.ylim(0.4983 , 0.501)
    plt.gca().yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))

    # plt.savefig('ExtendedDataFig_smallNetwork_BD_UG_p_u0.1_delta0.001.svg')
    plt.show()

    plt.figure(figsize=(5.2/1.6, 2.5/1.6))
    plt.bar(np.arange(len(q_grouped_sorted)), q_grouped_sorted, width=bar_width, label=r'$\overline{{M}}=\bar{d}$', 
            color=grouped_color, alpha=0.4, edgecolor=grouped_edgecolor)
    plt.bar(np.arange(len(q_ungrouped_sorted)), q_ungrouped_sorted, width=bar_width, label=r'$\overline{{M}}=1.0$', 
            color=ungrouped_color, edgecolor=ungrouped_edgecolor)

    plt.xlabel('Graph', fontsize = 14)
    plt.ylabel(r'$\langle q \rangle$', fontsize = 14)
    plt.title(r'$\delta = 0.001, \mu = 0.1$', fontsize = 14)
    plt.xlim(-0.5, 111.5)
    plt.ylim(0.497, 0.500)
    plt.yticks([0.497, 0.498, 0.4990, 0.5])

    plt.xticks([0, 29, 59, 89, 111], ['1', '30', '60', '90','112'])
    plt.gca().yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))

    # plt.savefig('ExtendedDataFig_smallNetwork_BD_UG_q_u0.1_delta0.001.svg')
    plt.show()

def TG_pq_barchart(p_ungrouped_file, p_grouped_file, q_ungrouped_file, q_grouped_file):
    with open(p_ungrouped_file, 'r') as f:
        p_ungrouped = np.array([float(line.strip()) for line in f.readlines()])

    with open(p_grouped_file, 'r') as f:
        p_grouped = np.array([float(line.strip()) for line in f.readlines()])

    with open(q_ungrouped_file, 'r') as f:
        q_ungrouped = np.array([float(line.strip()) for line in f.readlines()])

    with open(q_grouped_file, 'r') as f:
        q_grouped = np.array([float(line.strip()) for line in f.readlines()])

    ungrouped_color = '#F1BFB577'
    grouped_color = '#788AB6'
    ungrouped_edgecolor = '#ED3B0E77'
    grouped_edgecolor = '#5B679B'

    sorted_indices = np.argsort(p_ungrouped)
    p_ungrouped_sorted = p_ungrouped[sorted_indices]
    p_grouped_sorted = p_grouped[sorted_indices]
    q_ungrouped_sorted = q_ungrouped[sorted_indices]
    q_grouped_sorted = q_grouped[sorted_indices]

    bar_width = 0.99

    plt.figure(figsize=(5.2/1.6, 2.5/1.6))
    plt.bar(np.arange(len(p_grouped_sorted)), p_grouped_sorted, width=bar_width, label=r'$\overline{{M}}=\bar{d}$', 
            color=grouped_color, alpha=0.4, edgecolor=grouped_edgecolor)
    plt.bar(np.arange(len(p_ungrouped_sorted)), p_ungrouped_sorted, width=bar_width, label=r'$\overline{{M}}=1.0$', 
            color=ungrouped_color, edgecolor=ungrouped_edgecolor)

    plt.xlabel('Graph', fontsize = 14)
    plt.ylabel(r'$\langle p \rangle$', fontsize = 14)
    plt.title(r'$\delta = 0.001, \mu = 0.1$', fontsize = 14)
    plt.xlim(-0.5, 111.5)

    plt.xticks([0, 29, 59, 89, 111], ['1', '30', '60', '90','112'])
    plt.ylim(0.4972, 0.5)
    plt.yticks([0.498, 0.499, 0.50])
    plt.gca().yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))

    # plt.savefig('ExtendedDataFig_smallNetwork_BD_TG_s_p_u0.1_delta0.001.svg')
    plt.show()

    plt.figure(figsize=(5.2/1.6, 2.5/1.6))
    plt.bar(np.arange(len(q_grouped_sorted)), q_grouped_sorted, width=bar_width, label=r'$\overline{{M}}=\bar{d}$', 
            color=grouped_color, alpha=0.4, edgecolor=grouped_edgecolor)
    plt.bar(np.arange(len(q_ungrouped_sorted)), q_ungrouped_sorted, width=bar_width, label=r'$\overline{{M}}=1.0$', 
            color=ungrouped_color, edgecolor=ungrouped_edgecolor)

    plt.xlabel('Graph', fontsize = 14)
    plt.ylabel(r'$\langle r \rangle$', fontsize = 14)
    plt.title(r'$\delta = 0.001, \mu = 0.1$', fontsize = 14)
    plt.xlim(-0.5, 111.5)
    plt.ylim(0.493, 0.499)
    plt.yticks([0.493, 0.495, 0.497, 0.499])

    plt.xticks([0, 29, 59, 89, 111], ['1', '30', '60', '90','112'])
    plt.gca().yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))

    # plt.savefig('ExtendedDataFig_smallNetwork_BD_TG_q_u0.1_delta0.001.svg')
    plt.show()

def SHG_bc_barchart_resort(p_ungrouped_file, p_grouped_file):
    plt.rcParams['font.family'] = 'Arial'
    plt.rc('font', size=12.5)

    with open(p_ungrouped_file, 'r') as f:
        p_ungrouped = np.array([float(line.strip()) for line in f.readlines()])

    with open(p_grouped_file, 'r') as f:
        p_grouped = np.array([float(line.strip()) for line in f.readlines()])

    ungrouped_color = '#F1BFB577'  
    grouped_color = '#788AB6'  
    ungrouped_edgecolor = '#ED3B0E77' 
    grouped_edgecolor = '#5B679B'  

    positive_indices = np.where(p_ungrouped >= 0)[0]
    negative_indices = np.where(p_ungrouped < 0)[0]

    sorted_pos_indices = positive_indices[np.argsort(p_ungrouped[positive_indices])]
    sorted_neg_indices = negative_indices[np.argsort(p_ungrouped[negative_indices])]

    sorted_indices = np.concatenate((sorted_pos_indices, sorted_neg_indices))

    p_ungrouped_sorted = p_ungrouped[sorted_indices]
    p_grouped_sorted = p_grouped[sorted_indices]

    bar_width = 0.99
    plt.rcParams['font.family'] = 'Arial'
    plt.rc('font', size=12.5)

    fig, ax = plt.subplots(figsize=(3.25, 3.25))

    plt.bar(np.arange(len(p_grouped_sorted)), p_grouped_sorted, width=bar_width, label=r'$\overline{{M}}=\bar{d}$', 
            color=grouped_color, alpha=0.4, edgecolor=grouped_edgecolor)
    plt.bar(np.arange(len(p_ungrouped_sorted)), p_ungrouped_sorted, width=bar_width, label=r'$\overline{{M}}=1.0$', 
            color=ungrouped_color, edgecolor=ungrouped_edgecolor)

    plt.xlabel('Graph', fontsize = 14)
    plt.ylabel(r'$\langle h  \rangle$', fontsize = 14)
    plt.title(r'$\mu = 0.1$', fontsize = 14)

    plt.xticks([0, 29, 59, 89, 111], ['1', '30', '60', '90','112'])

    plt.ylim(0.498, 0.5)
    plt.yticks([0.498, 0.499, 0.50])
    plt.gca().yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))

    plt.legend(fontsize = 12)
    plt.xlim(-0.5, 111.5)
    # plt.savefig("ExtendedDataFig_smallNetwork_SHG_BD_bcCritical_u0.1.svg")
    plt.show()

def GEG_pq_barchart(p_ungrouped_file, p_grouped_file, q_ungrouped_file, q_grouped_file):
    with open(p_ungrouped_file, 'r') as f:
        p_ungrouped = np.array([float(line.strip()) for line in f.readlines()])

    with open(p_grouped_file, 'r') as f:
        p_grouped = np.array([float(line.strip()) for line in f.readlines()])

    with open(q_ungrouped_file, 'r') as f:
        q_ungrouped = np.array([float(line.strip()) for line in f.readlines()])

    with open(q_grouped_file, 'r') as f:
        q_grouped = np.array([float(line.strip()) for line in f.readlines()])

    ungrouped_color = '#F1BFB577'
    grouped_color = '#788AB6'
    ungrouped_edgecolor = '#ED3B0E77'
    grouped_edgecolor = '#5B679B'

    sorted_indices = np.argsort(p_ungrouped)
    p_ungrouped_sorted = p_ungrouped[sorted_indices]
    p_grouped_sorted = p_grouped[sorted_indices]
    q_ungrouped_sorted = q_ungrouped[sorted_indices]
    q_grouped_sorted = q_grouped[sorted_indices]

    bar_width = 0.99

    plt.figure(figsize=(5.2/1.6, 2.5/1.6))
    plt.bar(np.arange(len(p_grouped_sorted)), p_grouped_sorted, width=bar_width, label=r'$\overline{{M}}=\bar{d}$', 
            color=grouped_color, alpha=0.4, edgecolor=grouped_edgecolor)
    plt.bar(np.arange(len(p_ungrouped_sorted)), p_ungrouped_sorted, width=bar_width, label=r'$\overline{{M}}=1.0$', 
            color=ungrouped_color, edgecolor=ungrouped_edgecolor)

    plt.xlabel('Graph', fontsize = 14)
    plt.ylabel(r'$\langle \varphi \rangle$', fontsize = 14)
    plt.title(r'$\delta = 0.001, \mu = 0.1$', fontsize = 14)
    plt.xlim(-0.5, 111.5)

    plt.xticks([0, 29, 59, 89, 111], ['1', '30', '60', '90','112'])

    plt.ylim(0.498, 0.5)
    plt.yticks([0.498, 0.499, 0.50])
    plt.gca().yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))

    plt.show()

    plt.figure(figsize=(5.2/1.6, 2.5/1.6))
    plt.bar(np.arange(len(q_grouped_sorted)), q_grouped_sorted, width=bar_width, label=r'$\overline{{M}}=\bar{d}$', 
            color=grouped_color, alpha=0.4, edgecolor=grouped_edgecolor)
    plt.bar(np.arange(len(q_ungrouped_sorted)), q_ungrouped_sorted, width=bar_width, label=r'$\overline{{M}}=1.0$', 
            color=ungrouped_color, edgecolor=ungrouped_edgecolor)

    plt.xlabel('Graph', fontsize = 14)
    plt.ylabel(r'$\langle \omega \rangle$', fontsize = 14)
    plt.title(r'$\delta = 0.001, \mu = 0.1$', fontsize = 14)
    plt.xlim(-0.5, 111.5)

    plt.ylim(0.498, 0.5)

    plt.xticks([0, 29, 59, 89, 111], ['1', '30', '60', '90','112'])
    plt.gca().yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))

    # plt.savefig('ExtendedDataFig_smallNetwork_BD_GEG_q_u0.1_delta0.001.svg')
    plt.show()


def DTG_bc_barchart_resort(p_ungrouped_file, p_grouped_file):
    plt.rcParams['font.family'] = 'Arial'
    plt.rc('font', size=12.5)

    with open(p_ungrouped_file, 'r') as f:
        p_ungrouped = np.array([float(line.strip()) for line in f.readlines()])

    with open(p_grouped_file, 'r') as f:
        p_grouped = np.array([float(line.strip()) for line in f.readlines()])

    ungrouped_color = '#F1BFB577'  
    grouped_color = '#788AB6' 
    ungrouped_edgecolor = '#ED3B0E77' 
    grouped_edgecolor = '#5B679B'  

    positive_indices = np.where(p_ungrouped >= 0)[0]
    negative_indices = np.where(p_ungrouped < 0)[0]

    sorted_pos_indices = positive_indices[np.argsort(p_ungrouped[positive_indices])]
    sorted_neg_indices = negative_indices[np.argsort(p_ungrouped[negative_indices])]

    sorted_indices = np.concatenate((sorted_pos_indices, sorted_neg_indices))

    p_ungrouped_sorted = p_ungrouped[sorted_indices]
    p_grouped_sorted = p_grouped[sorted_indices]

    bar_width = 0.99
    plt.rcParams['font.family'] = 'Arial'
    plt.rc('font', size=12.5)

    fig, ax = plt.subplots(figsize=(3.25, 3.25))

    plt.bar(np.arange(len(p_grouped_sorted)), p_grouped_sorted, width=bar_width, label=r'$\overline{{M}}=\bar{d}$', 
            color=grouped_color, alpha=0.4, edgecolor=grouped_edgecolor)
    plt.bar(np.arange(len(p_ungrouped_sorted)), p_ungrouped_sorted, width=bar_width, label=r'$\overline{{M}}=1.0$', 
            color=ungrouped_color, edgecolor=ungrouped_edgecolor)

    plt.xlabel('Graph', fontsize = 14)
    plt.ylabel(r'$\langle p  \rangle$', fontsize = 14)
    plt.title(r'$\mu = 0.1$', fontsize = 14)

    plt.xticks([0, 29, 59, 89, 111], ['1', '30', '60', '90','112'])
    
    plt.axvline(x=60.5, color='gray', linestyle='--', linewidth=0.8)

    plt.ylim(0.498, 0.5)
    plt.yticks([0.498, 0.499, 0.50])
    plt.gca().yaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))

    plt.legend(fontsize = 12)
    plt.xlim(-0.5, 111.5)

    # plt.savefig("ExtendedDataFig_smallNetwork_DTG_BD_bcCritical_u0.1.svg")
    plt.show()

DG_bc_barchart_resort(
    'EDF8_data/donationGame_networkOfSize6_bc_ungrouped_u0.1.txt',
    'EDF8_data/donationGame_networkOfSize6_bc_grouped_u0.1.txt'
)

TG_pq_barchart(
    'EDF8_data/trustGame_networkOfSize6_p_ungrouped_u0.1.txt',
    'EDF8_data/trustGame_networkOfSize6_p_grouped_u0.1.txt',
    'EDF8_data/trustGame_networkOfSize6_q_ungrouped_u0.1.txt',
    'EDF8_data/trustGame_networkOfSize6_q_grouped_u0.1.txt')

UG_pq_barchart(
    'EDF8_data/ultimatumGame_networkOfSize6_p_ungrouped_u0.1.txt',
    'EDF8_data/ultimatumGame_networkOfSize6_p_grouped_u0.1.txt',
    'EDF8_data/ultimatumGame_networkOfSize6_q_ungrouped_u0.1.txt',
    'EDF8_data/ultimatumGame_networkOfSize6_q_grouped_u0.1.txt')

SHG_bc_barchart_resort(
    'EDF8_data/SHG_networkOfSize6_bc_ungrouped_u0.1.txt',
    'EDF8_data/SHG_networkOfSize6_bc_grouped_u0.1.txt'
)

GEG_pq_barchart(
    'EDF8_data/GEG_networkOfSize6_p_ungrouped_u0.1.txt',
    'EDF8_data/GEG_networkOfSize6_p_grouped_u0.1.txt',
    'EDF8_data/GEG_networkOfSize6_q_ungrouped_u0.1.txt',
    'EDF8_data/GEG_networkOfSize6_q_grouped_u0.1.txt')

DTG_bc_barchart_resort(
    'EDF8_data/DTG_networkOfSize6_bc_ungrouped_u0.1.txt',
    'EDF8_data/DTG_networkOfSize6_bc_grouped_u0.1.txt'
)