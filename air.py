from choose_sequences import sequences
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


def gauss(m, a, mu, sigma):
    return a * np.exp(- ((m - mu) ** 2) / (2 * sigma ** 2))


relative = False
combine = 1
start = 0
end = 50

o2 = np.arange(8, dtype=np.float)
co2 = np.arange(8, dtype=np.float)
breath = np.arange(1, 8)
colors = ["blue", "orange", "red", "green", "black", "yellow", "brown"]
amu0, p0, err0 = sequences('air_background.csv', False, combine, start, end, relative, False)
o2 = []
co2 = []
for i, c in zip(breath, colors):
    amu, p, err = sequences(f"air{i}.csv", False, combine, start, end, relative, False)
    mask_o2 = (amu > 31) & (amu < 33)
    mask_co2 = (amu > 43) & (amu < 45)
    popt, pcov = curve_fit(gauss, amu[mask_o2], p[mask_o2], p0=[max(p[mask_o2]), 32, 1])
    o2.append(popt[0])
    o2popt = popt
    popt, pcov = curve_fit(gauss, amu[mask_co2], p[mask_co2], p0=[max(p[mask_co2]), 44, 1])
    co2.append(popt[0])
    # print(popt[0])
    plt.plot(amu[mask_co2], p[mask_co2], color=c, label=f"breath {i}")
    plt.plot(np.linspace(43, 45, 100), gauss(np.linspace(43, 45, 100), *popt), ls="--", color=c)
    plt.plot(amu[mask_o2], p[mask_o2], color=c)
    plt.plot(np.linspace(31, 33, 100), gauss(np.linspace(31, 33, 100), *o2popt), ls="--", color=c)
    # exit()
    # p = p * p0[int(5 * 28 / combine - 1)] / p[int(5 * 28 / combine - 1)]
    # o2[i] = p[32 - 1]
    # co2[i] = p[44 - 1]
# plt.plot(amu,p,label='breathing {}'.format(i))
# plt.errorbar(amu, p, err, capsize=3, capthick=0.4 ,ecolor="black", elinewidth=0.4 ,fmt ='none')
#
plt.legend()
plt.show()
plt.plot(breath, o2, label=r'$O_2$')
plt.plot(breath, co2, label=r'$CO_2$')
plt.legend()
if relative:
    plt.ylabel(r"$p_{part}$ / $p_{tot}$ [%]")
else:
    plt.ylabel(r"$p$ [Torr]")
# plt.xlabel('amu')
plt.xlabel('breath')
plt.savefig("Report/DataResultsPlots/air.pdf")
plt.show()
