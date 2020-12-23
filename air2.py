from choose_sequences import sequences
import numpy as np
import matplotlib.pyplot as plt
from utils import fit_peak, gauss
from scipy.integrate import quad
from scipy.optimize import curve_fit
from nist import get_nist_peaks, nist_aprox
from gas_analysis import gas_analysis


def oxygen(n, n0, p):
    return n0 * p ** n


co2 = np.zeros(6, float)
co2err = np.zeros(6, float)
o2 = np.zeros(6, float)
o2err = np.zeros(6, float)
frac = np.zeros(6, float)
frac_err = np.zeros(6, float)
p = np.zeros(6, float)
perr = np.zeros(6, float)
co2rel = np.zeros(6, float)
co2relerr = np.zeros(6, float)
o2rel = np.zeros(6, float)
o2relerr = np.zeros(6, float)

for i in np.arange(6):
    if i != 0:
        air = nist_aprox(*gas_analysis('air{}.csv'.format(i), 'airbaseline.csv'), 'air{}'.format(i))
    else:
        air = nist_aprox(*gas_analysis('air00.csv', 'airbaseline.csv'), 'air00')
    co2[i] = air[1][11]
    co2err[i] = air[2][11]
    o2[i] = air[1][10]
    o2err[i] = air[2][10]
    p[i] = sum(air[1]) - air[1][5]  # -air[1][9]-air[1][8]
    perr[i] = np.sqrt(sum(air[2] ** 2))
    # perr[i] = np.sqrt(co2err[i]**2 + o2err[i]**2)

frac = co2 / o2
frac_err = np.sqrt((1 / o2) ** 2 * co2err ** 2 + co2 ** 2 / o2 ** 4 * o2err ** 2)
co2rel = co2 / p
co2relerr = np.sqrt(((1 / p) ** 2 - co2 / p ** 3) * co2err ** 2 + (co2 / p ** 2) ** 2 * perr ** 2)
o2rel = o2 / p
o2relerr = np.sqrt(((1 / p) ** 2 - o2 / p ** 3) * o2err ** 2 + (o2 / p ** 2) ** 2 * perr ** 2)

popt, pcov = curve_fit(oxygen, np.arange(0, 6), o2rel, sigma=o2relerr)
n = np.linspace(0, 10, 100)
print(f"a = {popt[1]:.4f} +/- {np.sqrt(pcov[1][1]):.4f}")
print(f"n0 = {100 * popt[0]:.4f} +/- {100 * np.sqrt(pcov[0][0]):.4f}")
plt.plot(n, 100 * oxygen(n, *popt), ls="--", color="red")


def carbondioxide(n, n0):
    return popt[0] - n0 * popt[1] ** n


co2relerr[0] = co2relerr[1]

poptco2, pcovco2 = curve_fit(carbondioxide, np.arange(0, 6), co2rel, sigma=co2relerr)
n = np.linspace(0, 10, 100)
print(f"n0 = {100 * poptco2[0]:.4f} +/- {100 * np.sqrt(pcovco2[0][0]):.4f}")
plt.plot(n, 100 * carbondioxide(n, *poptco2), ls="--", color="blue")

plt.errorbar(np.arange(0, 6), 100 * np.array(co2rel), np.array(co2relerr) * 100, color="blue", capsize=3, capthick=0.4,
             ecolor="black",
             elinewidth=0.4,
             fmt='.',
             label=r'CO$_2$')
plt.errorbar(np.arange(0, 6), 100 * np.array(o2rel), np.array(o2relerr) * 100, color="red", capsize=3, capthick=0.4,
             ecolor="black",
             elinewidth=0.4,
             fmt='.',
             label=r'O$_2$')

plt.xticks(np.arange(0, 6))
plt.xlim(-0.5, 5.5)
plt.legend()

plt.ylabel(r"amount [%]")
plt.xlabel('breath')
plt.savefig("Report/DataResultsPlots/air.pdf")
plt.show()

# fraction
plt.errorbar(np.arange(0, 6), frac, frac_err, color="red", capsize=3, capthick=0.4,
             ecolor="black",
             elinewidth=0.4,
             fmt='.', )
plt.xticks(np.arange(0, 6))
plt.xlim(-0.5, 5.5)
plt.ylabel(r"ratio  $p_{CO_2}/p_{O_2}$")
plt.xlabel('breath')
plt.savefig("Report/DataResultsPlots/air_relative.pdf")
plt.show()
frac = np.delete(frac, [2, 3])
frac_err = np.delete(frac_err, [2, 3])

plt.errorbar(np.arange(0, 4), frac, frac_err, color="red", capsize=3, capthick=0.4, ecolor="black", elinewidth=0.4,
             fmt='.', label=r'O$_2$')
plt.xticks([0, 1, 2, 3], [0, 1, 4, 5])
plt.xlim(-0.5, 3.5)
plt.ylabel(r"ratio  $p_{CO_2}/p_{O_2}$")
plt.xlabel('breath')
plt.savefig("Report/DataResultsPlots/air_relative_cuted.pdf")
plt.show()
print(p)
