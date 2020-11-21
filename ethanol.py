from choose_sequences import sequences
import numpy as np
import matplotlib.pyplot as plt
from utils import fit_peak
from nist import get_nist_peaks


relative = False
combine = 1
start = 0
end = 50

n = get_nist_peaks("ethanol", p_number=3)

amu, p, err = sequences('ethanol.csv', False, combine, start, end, relative, False, new=True)
fig, ax = plt.subplots(1, 1)

ax.bar(n.m, n.y * 1.1 * np.max(p) / 100000, width=0.9, alpha=0.5, color="orange", label="NIST Ethanol")

ax.bar(amu, p, width=0.1 * combine, color="red", label="ethanol")
ax.errorbar(amu, p, err, capsize=3, capthick=0.4, ecolor="black", elinewidth=0.4, fmt='none')
# popt, pcov = fit_peak(amu, p, m1=start, m2=end, ax=ax)
ax.legend()
if relative:
    ax.set_ylabel(r"$p_{part}$ / $p_{tot}$ [%]")
else:
    ax.set_ylabel(r"$p$ [Torr]")
ax.set_xlabel("amu")
plt.show()
