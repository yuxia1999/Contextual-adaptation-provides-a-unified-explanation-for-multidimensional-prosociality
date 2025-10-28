import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches

def plot_M_bc_new_withErrorBar(gamma, mut, degree):
    import numpy as np
    import matplotlib.pyplot as plt

    accumulated_RR = np.zeros((100, 151))
    accumulated_SF = np.zeros((100, 151))
    M_values = None

    plt.figure(figsize=(4.3, 4.3))

    for index in range(100):
        filenameSF = f'EDF2_data/DonationGame_SF_DB_bcCritical_50_degree{degree}_{index + 1}_M_variation.txt'
        data_matrixSF = np.loadtxt(filenameSF)
        rows = data_matrixSF[np.isclose(data_matrixSF[:, 0], gamma) & np.isclose(data_matrixSF[:, 1], mut)]
        M_values = rows[:, 2]
        bc_values = rows[:, 5]

        accumulated_SF[index, :] = bc_values

        filenameRR = f'EDF2_data/DonationGame_RR_DB_bcCritical_50_degree{degree}_{index + 1}_M_variation.txt'
        data_matrixRR = np.loadtxt(filenameRR)
        rows = data_matrixRR[np.isclose(data_matrixRR[:, 0], 0.0) & np.isclose(data_matrixRR[:, 1], mut)]
        bc_values = rows[:, 5]

        accumulated_RR[index, :] = bc_values

    M_values /= 50 

    mean_SF = np.mean(accumulated_SF, axis=0)
    std_SF = np.std(accumulated_SF, axis=0)
    mean_RR = np.mean(accumulated_RR, axis=0)
    std_RR = np.std(accumulated_RR, axis=0)

    plt.plot(M_values, mean_SF, color='#c1272c', label='BA')
    plt.fill_between(M_values, mean_SF - std_SF, mean_SF + std_SF, color='#c1272c', alpha=0.2)

    plt.plot(M_values, mean_RR, color='#0070bb', label='RR')
    plt.fill_between(M_values, mean_RR - std_RR, mean_RR + std_RR, color='#0070bb', alpha=0.2)

    plt.axhline(y=mean_SF[0], xmin=0.5/21, xmax=0.5/21+20/21, linestyle='--', color='#c1272c')
    plt.axhline(y=mean_RR[0], xmin=0.5/21, xmax=0.5/21+20/21, linestyle='--', color='#0070bb')

    plt.xlabel(r'$M$')
    plt.ylabel(r'$(b/c)^{\ast}$')
    plt.title(rf'$\bar{{d}} = {degree} , \gamma={gamma}, \mu = {mut}$')
    plt.legend()
    plt.xlim((50-7.5/2)/50, (200+7.5/2)/50)
    plt.ylim(1.9, 6.1)
    plt.yticks([2,3,4,5,6])

    # plt.savefig(f'figure3_M_vs_bc_SF_RR_deg{degree}_gamma{gamma}_mutation{mut}_1000graph_line.svg')
    plt.show()

def plot_u_bc_new_line(gamma, degree,M):
    accumulated_RR = np.zeros(50)
    accumulated_SF = np.zeros(50)

    plt.figure(figsize=(4.3,4.3))

    ite_number = 100
    accumulated_SF = np.zeros((ite_number,50))
    accumulated_RR = np.zeros((ite_number,50))

    for index in range(1,ite_number+1):
        filenameSF = f'EDF2_data/DonationGame_SF_DB_bcCritical_50_degree{degree}_{index}_M{M}.txt'
        data_matrixSF = np.loadtxt(filenameSF)

        rows = data_matrixSF[np.isclose(data_matrixSF[:, 0], gamma) ]
        u_values = rows[:, 1]
        bc_values = rows[:, 4]

        accumulated_SF[index-1, :] = bc_values

        filenameRR = f'EDF2_data/DonationGame_RR_DB_bcCritical_50_degree{degree}_{index}_M{M}.txt'
        data_matrixRR = np.loadtxt(filenameRR)

        rows = data_matrixRR[np.isclose(data_matrixRR[:, 0], 1.0) ]
        u_values = rows[:, 1]
        bc_values = rows[:, 4]

        accumulated_RR[index-1, :] = bc_values

    std_SF = np.std(accumulated_SF, axis=0)
    std_RR = np.std(accumulated_RR, axis=0)
    mean_SF = np.mean(accumulated_SF, axis=0)
    mean_RR = np.mean(accumulated_RR, axis=0)

    plt.fill_between(u_values[:48], mean_SF[:48] - std_SF[:48], mean_SF[:48] + std_SF[:48], color='#c1272c', alpha=0.2)
    plt.fill_between(u_values[:48], mean_RR[:48] - std_RR[:48], mean_RR[:48] + std_RR[:48], color='#0070bb', alpha=0.2)
    plt.plot(u_values[:48], mean_SF[:48], linewidth = 2, color='#c1272c', label=fr'BA')
    plt.plot(u_values[:48], mean_RR[:48], linewidth = 2, color= '#0070bb', label=fr'RR')

    plt.xlabel(r'Mutation rate, $\mu$')
    plt.ylabel(r'$(b/c)^{\ast}$')
    plt.title(rf'$\bar{{d}} = {degree} , \gamma={gamma}, \overline{{M}} = {M/50}$')
    plt.legend()

    x_max = 10**(3*0.5/20)
    x_min = 10**(-3*0.5/20-3)
    plt.xlim(x_min, x_max)

    plt.ylim(1.5,22.5)
    plt.yticks([2, 6, 10,14, 18, 22])
    plt.xscale('log')

    # plt.savefig(f'figure3_u_vs_bc_SF_RR_deg{degree}_gamma{gamma}_M{M}_line.svg')

    plt.show()


plot_M_bc_new_withErrorBar(0.0,0.01, 4) # panel a
plot_u_bc_new_line(0.0, 4, 100) # panel b
plot_M_bc_new_withErrorBar(1.0,0.01, 4) # panel c
plot_u_bc_new_line(1.0, 4, 100) # panel d
