from choose_sequences import sequences
import numpy as np
import matplotlib.pyplot as plt
from utils import fit_peak
from nist import get_nist_peaks


relative = False
combine = 1
start = 0
end = 50

n_butane = get_nist_peaks("butane", p_number=3)
n_propane = get_nist_peaks("propane", p_number=3)

amu, p, err = sequences('deo.csv', False, combine, start, end, relative, False, new=True)
fig, ax = plt.subplots(1, 1)

ax.bar(n_propane.m, n_propane.y * 1.1 * np.max(p) / 10000, width=0.9, alpha=0.5, color="blue", label="NIST Propane")
ax.bar(n_butane.m, n_butane.y * 1.1 * np.max(p) / 10000, width=0.9, alpha=0.5, color="orange", label="NIST Butane")

ax.bar(amu, p, width=0.1 * combine, color="red", label="deo")
ax.errorbar(amu, p, err, capsize=3, capthick=0.4, ecolor="black", elinewidth=0.4, fmt='none')
# popt, pcov = fit_peak(amu, p, m1=start, m2=end, ax=ax)
ax.legend()
if relative:
    ax.set_ylabel(r"$p_{part}$ / $p_{tot}$ [%]")
else:
    ax.set_ylabel(r"$p$ [Torr]")
ax.set_xlabel("amu")
plt.show()
