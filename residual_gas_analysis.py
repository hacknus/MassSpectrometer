import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from read_data import hist
from choose_sequences import sequences
from voltage_variation import double_plot, sequence_plot
from utils import fit_peak
from nist import get_nist_peaks, nist_aprox





#select significant peaks
combine=1
amu_min=0
amu_max= 140
relative=False
amu, p, err = sequences('xenonbaseline_highres.csv', False, combine, amu_min, amu_max,
                            relative, False,new=True)
atom=[1]
for i in np.arange(len(amu)-1):
    if err[i]<p[i] and err[i-1]<p[i-1]:
        peak=int(np.round(amu[i]))
        if np.array(atom)[-1]!=peak:
            atom.append(peak)
atom = np.array(atom)
amu_min = atom-1.3
amu_min[0]=0
amu_max = atom-0.3
integr_p = []
integr_p_err = []

#Gaus fit over selected peaks
for i,j,k in zip(amu_min,amu_max,atom):
    fig, ax = plt.subplots(1, 1)
    amu, p, err = sequences('xenonbaseline_highres.csv', False, combine, i, j, relative,plot=False,new= True)
    ax.plot(amu,p)
    ax.errorbar(amu, p, err, capsize=3, capthick=0.4, ecolor="black", elinewidth=0.4, fmt='none')
    ax.set_title('peak at {} amu'.format(k))
    ax.errorbar(amu, p, err, capsize=3, capthick=0.4, ecolor="black", elinewidth=0.4, fmt='none')
    ax.set_ylabel(r"$p$ [Torr]")
    ax.set_xlabel("amu")
    plt.tight_layout()
    plt.show()
    popt, pcov = fit_peak(amu, p, ax=ax)
    integr_p.append(popt[0]*popt[1]*np.sqrt(2*np.pi))
    integr_p_err.append(np.sqrt(2*np.pi*(popt[0]**2*pcov[1][1]+popt[1]**2*pcov[0][0])))
integr_p=np.array(integr_p)
integr_p_err= np.array(integr_p_err)

nist_aprox(atom,integr_p,integr_p_err)