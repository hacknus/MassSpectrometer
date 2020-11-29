from choose_sequences import sequences
import numpy as np
import matplotlib.pyplot as plt
from utils import fit_peak, gauss
from scipy.integrate import quad
from scipy.optimize import curve_fit

relative = False
combine = 1
start = 0
end = 50

breath = np.arange(1, 6)
colors = ["blue", "orange", "red", "green", "black", "yellow", "brown"]
amu0, p0, err0 = sequences('air0.csv', False, combine, start, end, relative, False, new=True)
fig, ax = plt.subplots(1, 1)
o2 = []
co2 = []
o2err = []
co2err = []
for i, c in zip(breath, colors):
    amu, p, err = sequences(f"air{i + 1}.csv", False, combine, start, end, relative, False, new=True)

    p -= 1e-10

    popt, pcov = fit_peak(amu, p, err, m1=30, m2=34, ax=ax)
    a, a_err = quad(gauss, popt[1] - 3 * popt[2], popt[1] + 3 * popt[2], args=tuple(popt))
    o2.append(a)
    o2err.append(2 * np.pi * np.sqrt(popt[0] ** 2 * pcov[2][2] + popt[2] ** 2 * pcov[0][0]))
    popt, pcov = fit_peak(amu, p, err, m1=42, m2=46, ax=ax)
    a, a_err = quad(gauss, popt[1] - 3 * popt[2], popt[1] + 3 * popt[2], args=tuple(popt))
    co2.append(a)
    co2err.append(2 * np.pi * np.sqrt(popt[0] ** 2 * pcov[2][2] + popt[2] ** 2 * pcov[0][0]))

    ax.plot(amu, p, color=c, label=f"breath {i}")
    s = np.sum(p)
    co2[-1] /= s
    o2[-1] /= s

# plt.xlim(30, 50)
plt.semilogy()
plt.ylim(1e-9, 1e-5)
plt.legend()
plt.show()


def oxygen(n, n0, p):
    return n0 * p ** n


def carbondioxide(n, n0, p):
    return n0 * (1 + p) ** n


popt, pcov = curve_fit(oxygen, np.arange(1, 6), o2)
n = np.linspace(0, 6, 100)
print(popt)
plt.plot(n, 100 * oxygen(n, *popt), ls="--", color="red")

popt, pcov = curve_fit(carbondioxide, np.arange(1, 6), co2)
n = np.linspace(0, 6, 100)
print(popt)
plt.plot(n, 100 * carbondioxide(n, *popt), ls="--", color="blue")

plt.errorbar(np.arange(1, 6), 100 * np.array(co2), co2err, color="blue", capsize=3, capthick=0.4, ecolor="black",
             elinewidth=0.4,
             fmt='.',
             label=r'$CO_2$')
plt.errorbar(np.arange(1, 6), 100 * np.array(o2), o2err, color="red", capsize=3, capthick=0.4, ecolor="black",
             elinewidth=0.4,
             fmt='.',
             label=r'$O_2$')

plt.xticks(np.arange(0, 6))
plt.xlim(0.5, 5.5)
plt.legend()
if relative:
    plt.ylabel(r"$p_{part}$ / $p_{tot}$ [%]")
else:
    plt.ylabel(r"$amount$ [%]")
plt.xlabel('breath')
plt.savefig("Report/DataResultsPlots/air.pdf")
plt.show()
