import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from read_data import hist
from choose_sequences import sequences
from utils import fit_peak

combine = 1
relative = False
colors = ['r','b','g','y']


def sequence_plot(parameter,ax,amu_min,amu_max,title,legend=True,relative=False):
    new = False
    if parameter == 'voltage':
            delta_ms = np.array(['50V','70V','90V'])
            path = "restgasspektrum_FARx7_20prozent_5res_{}.csv"
    if parameter == 'detector':
            delta_ms = np.array(['SEM', 'FAR'])
            path = "restgasspektrum_{}.csv"
    if parameter == 'baseline':
            path = 'xenonbaseline_highres.csv'
            delta_ms = ['']
            new = True
    for i, c in zip(delta_ms,colors):
        amu, p, err = sequences(path.format(i), False, combine, amu_min, amu_max,
                            relative, False,new)
        ax.plot(amu, p,color=c,label='{}'.format(i))
        ax.errorbar(amu, p, err, capsize=3, capthick=0.4, ecolor="black", elinewidth=0.4, fmt='none')
        if parameter == 'baseline':
            popt = fit_peak(amu, p, ax=ax)
        ax.set_title(title)
        if legend:
            ax.legend()
        if relative:
            ax.set_ylabel(r"$p_{part}$ / $p_{tot}$ [%]")
        else:
            ax.set_ylabel(r"$p$ [Torr]")
        ax.set_xlabel("amu")
        
def double_plot(parameter,title,amu_min1,amu_max1,title1,amu_min2,amu_max2,title2,relative=False):
    fig, ax = plt.subplots(1, 2)
    ax1 = ax[0]
    sequence_plot(parameter,ax1,amu_min1,amu_max1,title1,False,relative)
    ax2 = ax[1]  
    sequence_plot(parameter,ax2,amu_min2,amu_max2,title2,True,relative)
    fig.suptitle(r'variation of {}'.format(parameter))
    plt.tight_layout()
    plt.savefig("Report/DataResultsPlots/{}_variation_{}_and_{}.pdf".format(parameter,title1,title2))
    plt.plot()
    
if __name__ == '__main__':
    title = r'Acceleration voltage'
    param = 'voltage'

    double_plot(param,title,0,2,r'$H_2$',16,19,r'$H_2O$')
    
    double_plot(param,title,26,28,r'$N_2$',30,32,r'$O_2$')

    fig, ax = plt.subplots(1, 1)
    sequence_plot(param,ax,42,44,r'$CO_2$',True)
    fig.suptitle(r'variation of {}'.format(param))
    plt.tight_layout()
    plt.savefig("Report/DataResultsPlots/{}_variation_CO2.pdf".format(param))
    plt.show()

