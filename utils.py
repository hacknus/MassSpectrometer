from scipy.optimize import curve_fit
import numpy as np


def gauss(m, a, mu, sigma):
    return a * np.exp(- ((m - mu) ** 2) / (2 * sigma ** 2))


def fit_peak(m, p, m1=None, m2=None, ax=None):
    if m1 and m2:
        mask = (m > m1) & (m < m2)
        p = p[mask]
        m = m[mask]
    peak = (m[-1] - m[0]) / 2 + m[0]
    popt, pcov = curve_fit(gauss, m, p, p0=[max(p), peak, 1])
    print(f"found peak of pressure p={popt[0]} at m={popt[1]} with width={popt[2]}")
    if ax:
        m_lin = np.linspace(m[0], m[-1], 100)
        ax.plot(m_lin, gauss(m_lin, *popt), ls="--", color="red")
    return popt
