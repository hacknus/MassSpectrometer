import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from read_data import hist
from choose_sequences import sequences
from utils import fit_peak
import matplotlib.ticker as ticker

combine = 1
relative = False
colors = ['r','b','g','y']


def sequence_plot(parameter,ax,amu_min,amu_max,title=False,legend=True,relative=False):
    new = False
    width = 0.15
    if parameter == 'voltage':
            delta_ms = np.array(['50V','70V','90V'])
            path = "restgasspektrum_FARx7_20prozent_5res_{}.csv"
            label = delta_ms
            pos = [-width/3,0,width/3]
    if parameter == 'detector':
            delta_ms = np.array(['SEM', 'FAR'])
            path = "restgasspektrum_{}.csv"
            label = delta_ms
            pos = [-width/4, width/4]
    if parameter == 'baseline':
            path = 'xenonbaseline_highres.csv'
            delta_ms = ['']
            label = delta_ms
            new = True
            pos = 0
    if parameter == 'delta m':
            path = "restgasspektrum_FARx7_{}prozent.csv"
            delta_ms = np.array([-10, 5, 20, 35])
            label = np.array(['-10%', '5%', '20%', '35%'])
            pos = [-3*width/8,-1/8*width,1/8*width,3/8*width]
    if parameter == 'resolution':
            path = "restgasspektrum_FARx7_20prozent_{}res.csv"
            delta_ms = np.array([-5, 5, 10])
            label = np.array(['-5%', '5%', '10%'])
            pos = pos = [-width/3,0,width/3]
    
    for i, c, l, j in zip(delta_ms,colors, label, pos):
        amu, p, err = sequences(path.format(i), False, combine, amu_min, amu_max,
                            relative, False,new)
        #width = 0.9/len(delta_ms)
        ax.bar(amu+j, p,color=c,label='{}'.format(l), width = width/len(delta_ms))
        ax.errorbar(amu+j, p, err, capsize=3, capthick=0.4, color='black', ecolor='black', elinewidth=0.4, fmt='none')
        ax.set_xticks(np.arange(min(amu),max(amu)+0.4,0.4))
        ax.xaxis.set_minor_locator(ticker.MultipleLocator(0.2))
        if parameter == 'baseline':
            popt = fit_peak(amu, p, ax=ax)
       #     ax.set_title(title)
        if legend:
            ax.legend()
        if relative:
            ax.set_ylabel(r"$p_{part}$ / $p_{tot}$ [%]")
        else:
            ax.set_ylabel(r"$p$ [Torr]")
        ax.set_xlabel("m/q [amu/e]")
        ax.tick_params(axis='both', which='major', labelsize=8)
def double_plot(parameter,title,amu_min1,amu_max1,title1,amu_min2,amu_max2,title2,relative=False):
    fig, ax = plt.subplots(1, 2)
    ax1 = ax[0]
   
    sequence_plot(parameter,ax1,amu_min1,amu_max1,False,False,relative)
    ax2 = ax[1]  

    sequence_plot(parameter,ax2,amu_min2,amu_max2,False,True,relative)
    fig.suptitle(r'variation of {}'.format(parameter))
    
    plt.tight_layout()
    plt.savefig("Report/DataResultsPlots/{}_variation_{}_and_{}.pdf".format(parameter,title1,title2))
    plt.show()
    
if __name__ == '__main__':
    
    params = ['resolution','delta m','voltage','detector']
    for param in params:
        double_plot(param,'title',0,2,'h2',16,18.4,'h2o')
        double_plot(param,'title',26,28,'n2',30,32,'o2')

        fig, ax = plt.subplots(1, 1)
        sequence_plot(param,ax,42,44,r'$CO_2$',True)
        fig.suptitle(r'variation of {}'.format(param))
        plt.tight_layout()
        plt.savefig("Report/DataResultsPlots/{}_variation_CO2.pdf".format(param))
        plt.show()
        
        double_plot(param,'title',0,2,'h2',42,44,'co2')
    
