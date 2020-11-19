import numpy as np
import matplotlib.pyplot as plt
from utils import fit_peak
from utils import fit_xenon
import pandas as pd
from nist import get_nist_peaks


argon_path = "Data2Avg/argon.csv"
xenon_path = "Data2Avg/xenon.csv"

argon = pd.read_csv(argon_path)
xenon = pd.read_csv(xenon_path)
print(xenon.head())
amu = np.array(xenon.amu)
p = np.array(xenon.p)
err = np.array(xenon.err)

n = get_nist_peaks("xenon", 3)

fig, ax = plt.subplots(1, 1)
ax.bar(n.m, n.y*np.sum(p)/3500000, width=0.9, alpha=0.5, color="orange", label="NIST")
ax.bar(amu, p, width=0.1, color="blue", label="measurements")
ax.errorbar(amu, p, err, capsize=3, capthick=0.4, ecolor="black", elinewidth=0.4, fmt='none')
popt = fit_xenon(amu, p, m1=128, m2=136.5, ax=ax)
ax.legend()
ax.set_xlim(128,137)
ax.set_ylim(0,3e-8)
plt.show()


