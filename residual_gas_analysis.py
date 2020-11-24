import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from read_data import hist
from choose_sequences import sequences
from voltage_variation import double_plot, sequence_plot
from utils import fit_peak
from nist import get_nist_peaks


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
  





#Nist data
integr_p=np.array(integr_p)
integr_p_err= np.array(integr_p_err)

n_co2 = get_nist_peaks("co2", p_number=6)
amu_co2 = np.array(n_co2.m)
p_co2 = np.array(n_co2.y)
p_co2_scaled = integr_p[atom==44]/p_co2[amu_co2==44]*p_co2

n_o2 = get_nist_peaks("oxygen", p_number=6)
amu_o2 = np.array(n_o2.m)
p_o2 = np.array(n_o2.y)
p_o2_scaled = integr_p[atom==32]/p_o2[amu_o2==32]*p_o2
for i in list(set(amu_o2) & set(amu_co2)):
    p_o2_scaled[amu_o2==i] = p_o2_scaled[amu_o2==i]+p_co2_scaled[amu_co2==i]

n_n2 = get_nist_peaks("nitrogen", p_number=6)
amu_n2 = np.array(n_n2.m)
p_n2 = np.array(n_n2.y)
p_n2_scaled =  (integr_p[atom==28]-p_co2_scaled[amu_co2==28])/p_n2[amu_n2==28]*p_n2


n_h2o = get_nist_peaks("water", p_number=6)
amu_h2o = np.array(n_h2o.m)
p_h2o = np.array(n_h2o.y)
p_h2o_scaled = integr_p[atom==18]/p_h2o[amu_h2o==18]*p_h2o
#for i in list(set(amu_n2) & set(np.cascatenate(amu_co2,amu_o2)):
#    p_o2_scaled[amu_o2==i] = p_o2_scaled[amu_o2==i]+p_co2_scaled[amu_co2==i]

amu_h2=np.array([1,2])
p_h2 = np.array([2,100])       #mouse pad
p_h2_scaled = integr_p[atom==2]/100*p_h2

#plot
for i in amu_h2:
    if i in amu_h2o:
        p_h2_scaled[amu_h2==i]=p_h2_scaled[amu_h2==i]+p_h2o_scaled[amu_h2o==i]
    if i in amu_n2:
        p_h2_scaled[amu_h2==i]=p_h2_scaled[amu_h2==i]+p_n2_scaled[amu_n2==i]
    if i in amu_o2:
        p_h2_scaled[amu_h2==i]=p_h2_scaled[amu_h2==i]+p_o2_scaled[amu_o2==i]
    if i in amu_co2:
        p_h2_scaled[amu_h2==i]=p_h2_scaled[amu_h2==i]+p_co2_scaled[amu_co2==i]
        
for i in amu_h2o:
    if i in amu_n2:
        p_h2o_scaled[amu_h2o==i]=p_h2o_scaled[amu_h2o==i]+p_n2_scaled[amu_n2==i]
    if i in amu_o2:
        p_h2o_scaled[amu_h2o==i]=p_h2o_scaled[amu_h2o==i]+p_o2_scaled[amu_o2==i]
    if i in amu_co2:
        p_h2o_scaled[amu_h2o==i]=p_h2o_scaled[amu_h2o==i]+p_co2_scaled[amu_co2==i]
        
for i in amu_n2:
    if i in amu_o2:
        p_n2_scaled[amu_n2==i]=p_n2_scaled[amu_n2==i]+p_o2_scaled[amu_o2==i]
    if i in amu_co2:
        p_n2_scaled[amu_n2==i]=p_n2_scaled[amu_n2==i]+p_co2_scaled[amu_co2==i]
for i in amu_o2:
    if i in amu_co2:
        p_o2_scaled[amu_o2==i]=p_o2_scaled[amu_o2==i]+p_co2_scaled[amu_co2==i]
    

width=0.9
fig, ax = plt.subplots(1, 1)
ax.set_title('residual gas')
ax.bar(atom+width/4, integr_p, width=width/2, color="black",label='measured')
ax.bar(amu_h2-width/4,p_h2_scaled, width=width/2, color="m",label=r'nist $H_2$')
ax.bar(amu_h2o-width/4,p_h2o_scaled, width=width/2, color="red",label=r'nist $H_2O$')
ax.bar(amu_n2-width/4,p_n2_scaled, width=width/2, color="y",label=r'nist $N_2$')
ax.bar(amu_o2-width/4,p_o2_scaled, width=width/2, color="blue",label=r'nist $O_2$')
ax.bar(amu_co2-width/4,p_co2_scaled, width=width/2, color="green",label=r'nist $CO_2$')


ax.legend()
ax.errorbar(atom+width/4, integr_p, integr_p_err, capsize=3, capthick=0.4, ecolor="black", elinewidth=0.4, fmt='none')
if relative:
            ax.set_ylabel(r"$p_{part}$ / $p_{tot}$ [%]")
else:
    ax.set_ylabel(r"$p$ [Torr]")
    ax.set_xlabel("amu")
    ax.set_xticks(np.arange(0,max(atom)+1,2))
plt.savefig("Report/DataResultsPlots/residual_gas.pdf")
plt.show()


