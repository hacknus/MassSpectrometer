from choose_sequences import sequences
import numpy as np
import matplotlib.pyplot as plt
from utils import fit_peak
from nist import get_nist_peaks

relative = False
combine = 1
start = 0
end = 50

n = get_nist_peaks("ethanol", p_number=10)
nn = get_nist_peaks("nitrogen", p_number=4)
nco2 = get_nist_peaks("co2", p_number=4)
no2 = get_nist_peaks("oxygen", p_number=4)
nag = get_nist_peaks("argon", p_number=4)

amu, p, err = sequences('ethanol.csv', False, combine, start, end, relative, False, new=True)
fig, ax = plt.subplots(1, 1)

scale = np.max(p[(amu > 30) & (amu < 32)]) / 10000
scalen = np.max(p[(amu > 27) & (amu < 29)]) / 10000
scaleco2 = np.max(p[(amu > 43) & (amu < 45)]) / 10000
scaleo2 = np.max(p[(amu > 31.5) & (amu < 33)]) / 10000
scaleag = np.max(p[(amu > 39) & (amu < 41)]) / 10000

ax.bar(n.m, n.y * scale, width=0.9, alpha=0.5, color="red", label="NIST Ethanol")
ax.bar(nco2.m, nco2.y * scaleco2, width=0.9, alpha=0.5, color="orange", label=r"NIST CO$_2$")
ax.bar(nn.m, nn.y * scalen, width=0.9, alpha=0.5, color="green", label="NIST N$_2$")
ax.bar(no2.m, no2.y * scaleo2, width=0.9, alpha=0.5, color="yellow", label="NIST O$_2$")
ax.bar(nag.m, nag.y * scaleag, width=0.9, alpha=0.5, color="brown", label="NIST Ag")

ax.bar(amu, p, width=0.1 * combine, color="blue", label="ethanol + air")
ax.errorbar(amu, p, err, capsize=3, capthick=0.4, ecolor="black", elinewidth=0.4, fmt='none')
# popt, pcov = fit_peak(amu, p, m1=start, m2=end, ax=ax)
ax.legend()
if relative:
    ax.set_ylabel(r"$p_{part}$ / $p_{tot}$ [%]")
else:
    ax.set_ylabel(r"$p$ [Torr]")
ax.set_xlabel("amu")
ax.set_ylim(0, 1e-6)
ax.set_xlim(25, 47)
ax.set_xticks(range(25, 47))
plt.savefig("Report/DataResultsPlots/ethanol2.pdf")
plt.show()
