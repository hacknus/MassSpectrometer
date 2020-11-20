from choose_sequences import sequences
import numpy as np
import matplotlib.pyplot as plt
from utils import fit_peak

relative = False
combine = 1
start = 0
end = 50

o2 = np.arange(8, dtype=np.float)
co2 = np.arange(8, dtype=np.float)
breath = np.arange(1, 7)
colors = ["blue", "orange", "red", "green", "black", "yellow", "brown"]
amu0, p0, err0 = sequences('air0.csv', False, combine, start, end, relative, False, new=True)
fig, ax = plt.subplots(1, 1)
o2 = []
co2 = []
for i, c in zip(breath, colors):
    amu, p, err = sequences(f"air{i}.csv", False, combine, start, end, relative, False, new=True)
    popt = fit_peak(amu, p, err, m1=30, m2=34, ax=ax)
    o2.append(popt[0])
    popt = fit_peak(amu, p, err, m1=42, m2=46, ax=ax)
    co2.append(popt[0])
    # print(popt[0])
    ax.plot(amu, p, color=c, label=f"breath {i}")
    # exit()
    # p = p * p0[int(5 * 28 / combine - 1)] / p[int(5 * 28 / combine - 1)]
    # o2[i] = p[32 - 1]
    # co2[i] = p[44 - 1]
    s = np.sum(p)
    co2[-1] /= s
    o2[-1] /= s
# plt.plot(amu,p,label='breathing {}'.format(i))
# plt.errorbar(amu, p, err, capsize=3, capthick=0.4 ,ecolor="black", elinewidth=0.4 ,fmt ='none')
#
plt.xlim(30, 50)
plt.ylim(1e-9, 1e-6)
plt.legend()
plt.show()
plt.plot(np.arange(6), o2, label=r'$O_2$')
plt.plot(np.arange(6), co2, label=r'$CO_2$')
plt.errorbar(np.arange(6), co2, co2err, capsize=3, capthick=0.4 ,ecolor="black", elinewidth=0.4 ,fmt ='o',label=r'$CO_2$')
plt.legend()
if relative:
    plt.ylabel(r"$p_{part}$ / $p_{tot}$ [%]")
else:
    plt.ylabel(r"$p$ [Torr]")
# plt.xlabel('amu')
plt.xlabel('breath')
plt.savefig("Report/DataResultsPlots/air.pdf")
plt.show()
