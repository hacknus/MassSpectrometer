import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from read_data import hist
from choose_sequences import sequences
from voltage_variation import double_plot, sequence_plot
from utils import fit_peak
from nist import get_nist_peaks, nist_aprox




def gas_analysis(file,baseline_file=False):
    #select significant peaks
    combine=1
    amu_min=0
    amu_max= 140
    relative=False
    amu, p, err = sequences(file, baseline_file, combine, amu_min, amu_max,relative, False,new=True)
    plt.plot(amu,p)
    plt.show
    atom=[1]
    for i in np.arange(0,len(amu)-1,10):
        if amu[i] < 50:
            if  err[i]<p[i] and err[i+1]<p[i+1] and err[i-1]<p[i-1]:
                peak=int(np.round(amu[i]))
                if np.array(atom)[-1] < peak:
                    atom.append(peak)
                   
                  
              
        if amu[i] > 50:                                                     #amu verschiebt sich leicht, daher korrektur
            if  err[i]<p[i] and err[i+1]<p[i+1] and err[i+2]<p[i+2]:
                peak=int(np.round(amu[i]))
                if np.array(atom)[-1] < peak:
                    atom.append(peak)
                    if peak == 64:
                        atom.append(64.5)
                    if peak == 65:
                        atom.append(65.5)
                 
    atom = np.array(atom)
    amu_min = atom-1.3
    amu_min[0]=0
    amu_max = atom-0.3
    if max(amu) >60:
        for  i in np.arange(len(atom+2)):
            if atom[i]==64: 
                amu_min[i]+=+0.1
                amu_max[i]+=-0.7
            if atom[i]==64.5: 
                amu_min[i]+=0.25
                amu_max[i]-=0.25
            
            if atom[i]==65:
                amu_min[i]+=+0.25
                amu_max[i]+=-0.55
                
            if atom[i]==65.5: 
                amu_min[i]+=0.25
                amu_max[i]-=0.25
                
            if atom[i]==66:
                amu_min[i]+=0.1
                amu_max[i]+=-0.2
            if atom[i]==128:
                amu_min[i]+=0
                amu_max[i]+=-0.3
            if atom[i]==98:
                amu_min[i]+=0.2
                amu_max[i]+=-0.4
    integr_p = []
    integr_p_err = []

    #Gaus fit over selected peaks
    for i,j,k in zip(amu_min,amu_max,atom):
        #if k == 64 and file == 'xenon_highres.csv':
        #    integr_p.append(0)
        #    integr_p_err.append(0)
        #    continue
        fig, ax = plt.subplots(1, 1)
        amu, p, err = sequences(file, baseline_file, combine, i, j, relative,plot=False,new= True)
        ax.plot(amu,p)
        ax.errorbar(amu, p, err, capsize=3, capthick=0.4, ecolor="black", elinewidth=0.4, fmt='none')
        ax.set_title('peak at {} amu'.format(k))
        ax.errorbar(amu, p, err, capsize=3, capthick=0.4, ecolor="black", elinewidth=0.4, fmt='none')
        ax.set_ylabel(r"$p$ [Torr]")
        ax.set_xlabel("amu")
        plt.tight_layout()
        popt, pcov = fit_peak(amu, p, ax=ax)
        plt.show()
        if popt[2] > 0.4:
            popt[0] = 0
            pcov[0][0] = 0
        integr_p.append(popt[0]*abs(popt[2])*np.sqrt(2*np.pi))
        integr_p_err.append(np.sqrt(2*np.pi*(popt[0]**2*pcov[2][2]+popt[2]**2*pcov[0][0])))
    integr_p=np.array(integr_p)
    integr_p_err= np.array(integr_p_err)
    
    return atom, integr_p, integr_p_err
    