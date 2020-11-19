import numpy as np
import matplotlib.pyplot as plt
from utils import fit_peak
from utils import fit_xenon
import pandas as pd
from nist import get_nist_peaks


def make_plot_xenon(plot=True):
    xenon_path = "Data2Avg/xenon.csv"
    df = pd.read_csv(xenon_path)
    n = get_nist_peaks("xenon", 3)
    amu = np.array(df.amu)
    p = np.array(df.p)
    err = np.array(df.err)
    fig, ax = plt.subplots(1, 1)
    ax.bar(n.m, n.y * np.sum(p) / 3500000, width=0.9, alpha=0.5, color="orange", label="NIST")
    ax.bar(amu, p, width=0.1, color="blue", label="measurements")
    ax.errorbar(amu, p, err, capsize=3, capthick=0.4, ecolor="black", elinewidth=0.4, fmt='none')
    popt = fit_xenon(amu, p, m1=128, m2=136.5, ax=ax)
    ax.legend()
    ax.set_xlim(128, 137)
    ax.set_ylim(0, 3e-8)
    if plot:
        plt.show()
    return popt[3], np.array(n.m)[0]


def make_plot_argon(plot=True):
    argon_path = "Data2Avg/argon.csv"
    df = pd.read_csv(argon_path)
    n = get_nist_peaks("argon", 1)
    amu = np.array(df.amu)
    p = np.array(df.p)
    err = np.array(df.err)
    fig, ax = plt.subplots(1, 1)
    ax.bar(n.m, n.y * 1.1*np.max(p)/10000, width=0.9, alpha=0.5, color="orange", label="NIST")
    ax.bar(amu, p, width=0.1, color="blue", label="measurements")
    ax.errorbar(amu, p, err, capsize=3, capthick=0.4, ecolor="black", elinewidth=0.4, fmt='none')
    popt = fit_peak(amu, p, m1=38, m2=42, ax=ax)
    ax.legend()
    ax.set_xlim(38, 42)
    ax.set_ylim(0, 3e-7)
    if plot:
        plt.show()
    return popt[1], np.array(n.m)[0]


def calibrate_dataset(df):
    m_x, m_x_true = make_plot_xenon(True)
    m_a, m_a_true = make_plot_argon(True)
    df.amu = m_a_true + (m_x_true - m_a_true) / (m_x - m_a) * (df.amu - m_a)
    return df
