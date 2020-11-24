import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from read_data import hist
from choose_sequences import sequences
from voltage_variation import double_plot, sequence_plot

colors = ['r']
title = 'base line'
param = 'baseline'

double_plot(param,title,0,2.2,r'$H_2$',15,20,r'$H_2O$')

double_plot(param,title,26.8,27.5,r'$N_2$',30,33,r'$O_2$')

fig, ax = plt.subplots(1, 1)
sequence_plot(param,ax,42.5,43.8,r'$CO_2$',True)
fig.suptitle(r'variation of {}'.format(param))
plt.tight_layout()
plt.savefig("Report/DataResultsPlots/{}_variation_CO2.pdf".format(param))
plt.show()


combine=1
amu_min=0
amu_max= 140
relative=False
fig, ax = plt.subplots(1, 1)
amu, p, err = sequences('xenonbaseline_highres.csv', False, combine, amu_min, amu_max,
                            relative, False,new=True)
for i in np.arange(len(amu)-1):
    if err[i]>p[i]:
        p[i] = 0
        err[i] = 0
    else:    print(amu[i])
        
ax.plot(amu, p)
ax.errorbar(amu, p, err, capsize=3, capthick=0.4, ecolor="black", elinewidth=0.4, fmt='none')
ax.set_title(title)
plt.tight_layout()
#plt.savefig("Report/DataResultsPlots/residual_gas_an.pdf".format(param))
plt.show()