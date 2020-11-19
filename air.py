from choose_sequences import sequences
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


def gauss(m, a, mu, sigma):
    return a * np.exp(- (m - mu) ** 2 / (2 * sigma ** 2))


relative = False
combine = 1
start = 0
end = 50
gause=False

o2 = np.arange(8, dtype=np.float)
co2 = np.arange(8, dtype=np.float)
breath = np.arange(8)
#amu_base, p_base, err_base = sequences('air1.csv', False, combine, start, end, relative, False)
amu0,p0,err0 = sequences('air_background.csv', False, combine, start, end, relative, False)
print(sequences('xenon.csv', False, combine, start, end, relative, False)[1])
o2[0] = p0[int(32*5/combine - 3)]/sum(p0)
co2[0] = p0[int(44*5/combine - 3)]/sum(p0)
#n2[0]=p0[]

for i in np.arange(1, 8):
    amu, p, err = sequences("air{}.csv".format(i), False, combine, start, end, relative, False)
    p = p/sum(p)
    if gause==True:
        mask_o2 = (amu > 31) & (amu < 33)
        mask_co2 = (amu > 43) & (amu < 45)
        popt, pcov = curve_fit(gauss, amu[mask_o2], p[mask_o2], p0=[max(p[mask_o2]), 32, 1])
        o2[i] = popt[0]
        o2popt = popt
        popt, pcov = curve_fit(gauss, amu[mask_co2], p[mask_co2], p0=[max(p[mask_co2]), 44, 1])
        co2[i] = popt[0]
        print(popt[0])
        plt.plot(amu[mask_co2], p[mask_co2], color="orange")
        plt.plot(np.linspace(43, 45, 100), gauss(np.linspace(43, 45, 100), *popt), ls="--", color="red")
        plt.plot(amu[mask_o2], p[mask_o2], color="blue")
        plt.plot(np.linspace(31, 33, 100), gauss(np.linspace(31, 33, 100), *o2popt), ls="--", color="red")
        plt.show()
    else:
        o2[i] = p[int(32*5/combine )-3]
        co2[i] = p[int(44*5/combine)-3]
# plt.plot(amu,p,label='breathing {}'.format(i))
# plt.errorbar(amu, p, err, capsize=3, capthick=0.4 ,ecolor="black", elinewidth=0.4 ,fmt ='none')
#
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
