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


def carbondioxide(n, n0, p):
    return n0 * (1 - p ** n)


co2 = np.zeros(6,float)
co2err = np.zeros(6,float)
o2 = np.zeros(6,float)
o2err = np.zeros(6,float)
frac = np.zeros(6,float)
frac_err = np.zeros(6,float)

for i in np.arange(6):
    air=nist_aprox(*gas_analysis('air{}.csv'.format(i),'airbaseline.csv'),'air{}'.format(i))
    co2[i] = air[1][11]
    co2err[i] = air[2][11]
    o2[i] = air[1][10]
    o2err[i] = air[2][10]
frac = co2/o2
frac_err = np.sqrt((1/o2)**2*co2err**2 + co2**2/o2**4 * o2err**2)


#popt, pcov = curve_fit(carbondioxide, np.arange(0, 6), frac, sigma=frac_err, p0 = [1,1])
#n = np.linspace(0, 6, 100)
#print(f"a = {popt[1]:.4f} +/- {np.sqrt(pcov[1][1]):.4f}")
#print(f"n0 = {100 * popt[0]:.4f} +/- {100 * np.sqrt(pcov[0][0]):.4f}")
#plt.plot(n, 100 * carbondioxide(n, *popt), ls="--", color="blue")

plt.errorbar(np.arange(0, 6), 100 * frac, frac_err * 100, color="red", capsize=3, capthick=0.4,
             ecolor="black",
             elinewidth=0.4,
             fmt='.',
             label=r'O$_2$')



plt.xticks(np.arange(0, 6))
plt.xlim(-0.5, 5.5)
plt.ylabel(r"partial pressure CO$_2$/O$_2$ [%]")
plt.xlabel('breath')
plt.savefig("Report/DataResultsPlots/air_relative.pdf")
plt.show()
frac =np.delete(frac,[2,3])
frac_err = np.delete(frac_err,[2,3])
plt.errorbar(np.arange(0,4), 100 * frac, frac_err * 100, color="red", capsize=3, capthick=0.4, ecolor="black", elinewidth=0.4, fmt='.', label=r'O$_2$')

#popt, pcov = curve_fit(carbondioxide, np.arange(0, 4), frac, sigma=frac_err, p0 = [1,1])
#n = np.linspace(0, 4, 100)
#print(f"a = {popt[1]:.4f} +/- {np.sqrt(pcov[1][1]):.4f}")
#print(f"n0 = {100 * popt[0]:.4f} +/- {100 * np.sqrt(pcov[0][0]):.4f}")
#plt.plot(n, 100 * carbondioxide(n, *popt), ls="--", color="blue")


plt.xticks([0,1,2,3],[0,1,4,5 ])
plt.xlim(-0.5, 3.5)
plt.ylabel(r"partial pressure CO$_2$/O$_2$ [%]")
plt.xlabel('breath')
plt.savefig("Report/DataResultsPlots/air_relative_cuted.pdf")
plt.show()