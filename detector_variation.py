import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os
from read_data import hist
from choose_sequences import sequences

combine = 1
relative = False
delta_ms = np.array(['_SEMx10', '_FARx7_20prozent_10res'])
# delta_ms= np.array(['_SEMx10_50prozent_10res','_FARx7_50prozent_10res'])
detector = np.array(['SEM', 'FAR'])

fig, ax = plt.subplots(1, 2)
ax1 = ax[0]
ax2 = ax[1]
amu_min = 0
amu_max = 2
for i in [0, 1]:
    amu, p, err = sequences("restgasspektrum{}.csv".format(delta_ms[i]), False, combine, amu_min, amu_max, relative,
                            False)
    ax1.plot(amu, p)
    ax1.errorbar(amu, p, err, capsize=3, capthick=0.4, ecolor="black", elinewidth=0.4, fmt='none')
ax1.set_title(r'$H_2$')
if relative:
    ax1.set_ylabel(r"$p_{part}$ / $p_{tot}$ [%]")
else:
    ax1.set_ylabel(r"$p$ [Torr]")
ax1.set_xlabel("amu")

amu_min = 16
amu_max = 19
for i in [0, 1]:
    amu, p, err = sequences("restgasspektrum{}.csv".format(delta_ms[i]), False, combine, amu_min, amu_max, relative,
                            False)
    ax2.plot(amu, p, label=detector[i])
    ax2.errorbar(amu, p, err, capsize=3, capthick=0.4, ecolor="black", elinewidth=0.4, fmt='none')
ax2.set_title(r'$H_2O$')
ax2.legend()
if relative:
    ax2.set_ylabel(r"$p_{part}$ / $p_{tot}$ [%]")
else:
    ax1.set_ylabel(r"$p$ [Torr]")
ax2.set_xlabel("amu")

fig.suptitle(r'variation of detector')
plt.tight_layout()
plt.savefig("Report/DataResultsPlots/delta_m_variation_H2_and_H20.pdf")
plt.plot()

fig, ax = plt.subplots(1, 2)
ax3 = ax[0]
ax4 = ax[1]
amu_min = 27
amu_max = 28
for i in [0, 1]:
    amu, p, err = sequences("restgasspektrum{}.csv".format(delta_ms[i]), False, combine, amu_min, amu_max, relative,
                            False)
    ax3.plot(amu, p, label=detector[i])
    ax3.errorbar(amu, p, err, capsize=3, capthick=0.4, ecolor="black", elinewidth=0.4, fmt='none')
ax3.set_title(r'$N_2$')

if relative:
    ax3.set_ylabel(r"$p_{part}$ / $p_{tot}$ [%]")
else:
    ax3.set_ylabel(r"$p$ [Torr]")
ax3.set_xlabel("amu")

amu_min = 31
amu_max = 32
for i in [0, 1]:
    amu, p, err = sequences("restgasspektrum{}.csv".format(delta_ms[i]), False, combine, amu_min, amu_max, relative,
                            False)
    ax4.plot(amu, p, label=detector[i])
    ax4.errorbar(amu, p, err, capsize=3, capthick=0.4, ecolor="black", elinewidth=0.4, fmt='none')
ax4.set_title(r'$O_2$')
ax4.legend()
# if relative:
#    ax2.set_ylabel(r"$p_{part}$ / $p_{tot}$ [%]")
# else: ax2.set_ylabel(r"$p$ [Torr]")
ax4.set_xlabel("amu")

fig.suptitle(r'variation of detector')
plt.tight_layout()
plt.savefig("Report/DataResultsPlots/detector_variation_N2_and_O2.pdf")
plt.plot()

fig, ax = plt.subplots(1, 1)
ax5 = ax
amu_min = 43
amu_max = 44
for i in [0, 1]:
    amu, p, err = sequences("restgasspektrum{}.csv".format(delta_ms[i]), False, combine, amu_min, amu_max, relative,
                            False)
    ax5.plot(amu, p, label=detector[i])
    ax5.errorbar(amu, p, err, capsize=3, capthick=0.4, ecolor="black", elinewidth=0.4, fmt='none')
ax5.set_title(r'$CO_2$')
ax5.legend()
if relative:
    ax5.set_ylabel(r"$p_{part}$ / $p_{tot}$ [%]")
else:
    ax5.set_ylabel(r"$p$ [Torr]")
ax5.set_xlabel("amu")

fig.suptitle(r'variation of detector')
plt.tight_layout()
plt.savefig("Report/DataResultsPlots/detector_variation_CO2.pdf")
plt.show()
