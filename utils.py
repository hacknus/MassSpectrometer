from scipy.optimize import curve_fit
import numpy as np


def gauss(m, a, mu, sigma):
    return a * np.exp(- ((m - mu) ** 2) / (2 * sigma ** 2))


def xenon(m, a, a2, a3, mu, sigma, sigma2, sigma3):
    return gauss(m, a, mu, sigma) + gauss(m, a2, mu-1, sigma2) + gauss(m, a3, mu-3, sigma3)

def xenon_half(m, a, a2, a3, mu, sigma, sigma2, sigma3):
    return gauss(m, a, mu, sigma) + gauss(m, a2, mu+1, sigma2)+gauss(m, a3, mu+1.5, sigma3)




def fit_peak(m, p, err=None, m1=None, m2=None, ax=None):
    if m1 and m2:
        mask = (m > m1) & (m < m2)
        p = p[mask]
        m = m[mask]
    peak = (m[-1] - m[0]) / 2 + m[0]
    if type(err) != None:
        popt, pcov = curve_fit(gauss, m, p, p0=[max(p), peak, 1], maxfev=2000)
    else:
        popt, pcov = curve_fit(gauss, m, p, p0=[max(p), peak, 1], sigma=err, maxfev=2000)
    print(f"found peak of pressure p={popt[0]} at m={popt[1]:.5f} +/- {pcov[1][1]:.5f} with width={popt[2]}")
    if ax:
        m_lin = np.linspace(m[0], m[-1], 1000)
        ax.plot(m_lin, gauss(m_lin, *popt), ls="--", color="red")
        ax.fill_between(m_lin,gauss(m_lin, *popt),0,color="red",alpha=0.5)
    return popt, pcov


def fit_xenon(m, p, err=None, m1=None, m2=None, ax=None):
    if m1 and m2:
        mask = (m > m1) & (m < m2)
        p = p[mask]
        m = m[mask]
    peak = (m[-1] - m[0]) / 2 + m[0]
    try:
        if type(err) != None:
            popt, pcov = curve_fit(xenon, m, p, p0=[max(p), max(p), max(p), peak, 0.1, 0.1, 0.1], maxfev=10000)
        else:
            popt, pcov = curve_fit(xenon, m, p, p0=[max(p), max(p), max(p), peak, 0.1, 0.1, 0.1], sigma=err, maxfev=10000)
        print(f"found first peak of pressure p={popt[0]} at m={popt[3]:.5f} +/- {pcov[3][3]:.5f} with width={popt[3]:.2f}")
        print(f"found second peak of pressure p={popt[1]} at m={popt[3]-1:.5f} +/- {pcov[3][3]:.5f} with width={popt[4]:.2f}")
        print(f"found thirs peak of pressure p={popt[2]} at m={popt[3]-3:.5f} +/- {pcov[3][3]:.5f} with width={popt[5]:.2f}")
    except:
        popt = [max(p), max(p), max(p), peak, 0.1, 0.1, 0.1]
    if ax:
        m_lin = np.linspace(m[0], m[-1], 100)
        ax.plot(m_lin, xenon(m_lin, *popt), ls="--", color="red")
    return popt, pcov

def fit_xenon_half(m, p, err=None, m1=None, m2=None, ax=None):
    if m1 and m2:
        mask = (m > m1) & (m < m2)
        p = p[mask]
        m = m[mask]
    peak = (m[-1] - m[0]) / 2 + m[0]
    try:
        if type(err) != None:
            popt, pcov = curve_fit(xenon_half, m, p, p0=[max(p), max(p), max(p), peak, 0.1, 0.1, 0.1], maxfev=10000)
        else:
            popt, pcov = curve_fit(xenon_half, m, p, p0=[max(p), max(p), max(p), peak, 0.1, 0.1, 0.1], sigma=err, maxfev=10000)
        print(f"found first peak of pressure p={popt[0]} at m={popt[3]:.5f} +/- {pcov[3][3]:.5f} with width={popt[3]:.2f}")
        print(f"found second peak of pressure p={popt[1]} at m={popt[3]-1:.5f} +/- {pcov[3][3]:.5f} with width={popt[4]:.2f}")
        print(f"found thirs peak of pressure p={popt[2]} at m={popt[3]-3:.5f} +/- {pcov[3][3]:.5f} with width={popt[5]:.2f}")
    except:
        popt = [max(p),max(p),max(p), peak, 0.1, 0.1, 0.1]
    if ax:
        m_lin = np.linspace(m[0], m[-1], 100)
        ax.plot(m_lin, xenon_half(m_lin, *popt), ls="--", color="red")
    return popt, pcov