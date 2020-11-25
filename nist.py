from jcamp import JCAMP_reader
import pandas as pd
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from read_data import hist
from choose_sequences import sequences
from voltage_variation import double_plot, sequence_plot
from utils import fit_peak


def get_nist_peaks(name, p_number=1):
    path = "NIST/" + name + ".jdx"
    d = JCAMP_reader(path)
    m = d["x"]
    y = d["y"]
    nist = {"m": m,
            "y": y}
    n = pd.DataFrame(data=nist)
    n.sort_values("y", inplace=True, ascending=False)
    n.reset_index()
    n = n[:p_number]
    return n

<<<<<<< HEAD
def scale_nist(atom,integr_p,gas,number=12):
=======

def scale_nist(atom, integr_p, gas, number=12):
>>>>>>> 202ac3ab2ceb119e8e6fc80624c1a0e57aaed4bc
    if gas == 'h2':
        amu = np.array([1, 2])
        p = np.array([2, 100])
    else:
        nist = get_nist_peaks(gas, p_number=number)
        amu = np.array(nist.m)
        p = np.array(nist.y)
    amu_max = amu[max(p) == p]
    if max(amu_max) in atom:
        return amu, integr_p[atom == amu_max] / p[amu == amu_max] * p
    else:
        return amu, np.zeros(len(amu), int)


def nist_aprox(atom, integr_p, integr_p_err, plot_name, ethanol=False):
    # scale nist data
    amu_co2, p_co2_scaled = scale_nist(atom, integr_p, 'co2')
    amu_o2, p_o2_scaled = scale_nist(atom, integr_p, 'oxygen')
    p_o2_0 = p_o2_scaled
    amu_argon, p_argon_scaled = scale_nist(atom, integr_p, 'argon')
    amu_xenon, p_xenon_scaled = scale_nist(atom, integr_p, 'xenon')
    amu_krypton, p_krypton_scaled = scale_nist(atom, integr_p, 'krypton')
    amu_h2o, p_h2o_scaled = scale_nist(atom, integr_p, 'water')
    p_h2o_0 = p_h2o_scaled
    amu_butane, p_butane_scaled = scale_nist(atom, integr_p, 'butane')
    p_butane_0 = p_butane_scaled
    amu_h2, p_h2_scaled = scale_nist(atom, integr_p, 'h2')
    amu_ethanol, p_ethanol_scaled = scale_nist(atom, integr_p, 'ethanol')
    p_ethanol_0 = p_ethanol_scaled

    if ethanol:
        p_butane_scaled = np.zeros(len(p_butane_scaled))
        p_butane_0 = p_butane_scaled

    n_propane = get_nist_peaks("propane", p_number=6)
    amu_propane = np.array(n_propane.m)
    p_propane = np.array(n_propane.y)
    if 29 in atom:
        p_propane_scaled = (integr_p[atom == 29] - p_butane_scaled[amu_butane == 29]) / p_propane[
            amu_propane == 29] * p_propane
    else:
        p_propane_scaled = np.zeros(len(amu_propane), int)
    p_propane_0 = p_propane_scaled
    if ethanol:
        p_propane_scaled = np.zeros(len(p_propane_scaled))
        p_propane_0 = p_propane_scaled

    n_n2 = get_nist_peaks("nitrogen", p_number=6)
    amu_n2 = np.array(n_n2.m)
    p_n2 = np.array(n_n2.y)
    if 28 in atom:
        p_n2_scaled = (integr_p[atom == 28] - p_co2_scaled[amu_co2 == 28] - p_butane_scaled[amu_butane == 28] -
                       p_propane_scaled[amu_propane == 28]) / p_n2[amu_n2 == 28] * p_n2
    else:
        p_n2_scaled = np.zeros(len(amu_n2), int)
    p_n2_0 = p_n2_scaled

    # make sum of the bars
    if ethanol != True:
        for i in amu_propane:
            if i in amu_butane:
                p_propane_scaled[amu_propane == i] = p_propane_scaled[amu_propane == i] + p_butane_scaled[
                    amu_butane == i]
            if i in amu_h2o:
                p_propane_scaled[amu_propane == i] = p_propane_scaled[amu_propane == i] + p_h2o_scaled[amu_h2o == i]
            if i in amu_n2:
                p_propane_scaled[amu_propane == i] = p_propane_scaled[amu_propane == i] + p_n2_scaled[amu_n2 == i]
            if i in amu_o2:
                p_propane_scaled[amu_propane == i] = p_propane_scaled[amu_propane == i] + p_o2_scaled[amu_o2 == i]
            if i in amu_co2:
                p_propane_scaled[amu_propane == i] = p_propane_scaled[amu_propane == i] + p_co2_scaled[amu_co2 == i]

        for i in amu_butane:
            if i in amu_h2o:
                p_butane_scaled[amu_butane == i] = p_butane_scaled[amu_butane == i] + p_h2o_scaled[amu_h2o == i]
            if i in amu_n2:
                p_butane_scaled[amu_butane == i] = p_butane_scaled[amu_butane == i] + p_n2_scaled[amu_n2 == i]
            if i in amu_o2:
                p_butane_scaled[amu_butane == i] = p_butane_scaled[amu_butane == i] + p_o2_scaled[amu_o2 == i]
            if i in amu_co2:
                p_butane_scaled[amu_butane == i] = p_butane_scaled[amu_butane == i] + p_co2_scaled[amu_co2 == i]
    else:
        for i in amu_ethanol:
            if i in amu_h2o:
                p_ethanol_scaled[amu_ethanol == i] = p_ethanol_scaled[amu_ethanol == i] + p_h2o_scaled[amu_h2o == i]
            if i in amu_n2:
                p_ethanol_scaled[amu_ethanol == i] = p_ethanol_scaled[amu_ethanol == i] + p_n2_scaled[amu_n2 == i]
            if i in amu_o2:
                p_ethanol_scaled[amu_ethanol == i] = p_ethanol_scaled[amu_ethanol == i] + p_o2_scaled[amu_o2 == i]
            if i in amu_co2:
                p_ethanol_scaled[amu_ethanol == i] = p_ethanol_scaled[amu_ethanol == i] + p_co2_scaled[amu_co2 == i]

    for i in amu_h2o:
        if i in amu_n2:
            p_h2o_scaled[amu_h2o == i] = p_h2o_scaled[amu_h2o == i] + p_n2_scaled[amu_n2 == i]
        if i in amu_o2:
            p_h2o_scaled[amu_h2o == i] = p_h2o_scaled[amu_h2o == i] + p_o2_scaled[amu_o2 == i]
        if i in amu_co2:
            p_h2o_scaled[amu_h2o == i] = p_h2o_scaled[amu_h2o == i] + p_co2_scaled[amu_co2 == i]

    for i in amu_n2:
        if i in amu_o2:
            p_n2_scaled[amu_n2 == i] = p_n2_scaled[amu_n2 == i] + p_o2_scaled[amu_o2 == i]
        if i in amu_co2:
            p_n2_scaled[amu_n2 == i] = p_n2_scaled[amu_n2 == i] + p_co2_scaled[amu_co2 == i]
    for i in amu_o2:
        if i in amu_co2:
            p_o2_scaled[amu_o2 == i] = p_o2_scaled[amu_o2 == i] + p_co2_scaled[amu_co2 == i]

    # plot bars
    width = 0.9
    fig, ax = plt.subplots(1, 1)
    ax.set_title(plot_name)

    ax.bar(atom + width / 4, integr_p, width=width / 2, color="black", label='measured')
    if sum(p_h2_scaled) != 0:
        ax.bar(amu_h2 - width / 4, p_h2_scaled, width=width / 2, color="m", label=r'nist $H_2$')
    if sum(p_krypton_scaled) != 0:
        ax.bar(amu_krypton - width / 4, p_krypton_scaled, width=width / 2, color="cyan", label=r'nist Kr')
    if sum(p_argon_scaled) != 0:
        ax.bar(amu_argon - width / 4, p_argon_scaled, width=width / 2, color="steelblue", label=r'nist Ar')
    if sum(p_xenon_scaled) != 0:
        ax.bar(amu_xenon - width / 4, p_xenon_scaled, width=width / 2, color="deepskyblue", label=r'nist Xe')
    if sum(p_ethanol_scaled) != 0:
        ax.bar(amu_ethanol - width / 4, p_ethanol_scaled, width=width / 2, color="deepskyblue", label=r'nist ethanol')
    if sum(p_propane_0) != 0:
        ax.bar(amu_propane - width / 4, p_propane_scaled, width=width / 2, color="pink", label=r'nist propane')
    if sum(p_butane_0) != 0:
        ax.bar(amu_butane - width / 4, p_butane_scaled, width=width / 2, color="orange", label=r'nist butane')
    if sum(p_h2o_0) != 0:
        ax.bar(amu_h2o - width / 4, p_h2o_scaled, width=width / 2, color="red", label=r'nist $H_2O$')
    if sum(p_n2_0) != 0:
        ax.bar(amu_n2 - width / 4, p_n2_scaled, width=width / 2, color="y", label=r'nist $N_2$')
    if sum(p_o2_0) != 0:
        ax.bar(amu_o2 - width / 4, p_o2_scaled, width=width / 2, color="blue", label=r'nist $O_2$')
    if sum(p_co2_scaled) != 0:
        ax.bar(amu_co2 - width / 4, p_co2_scaled, width=width / 2, color="green", label=r'nist $CO_2$')
    ax.legend()
    ax.errorbar(atom + width / 4, integr_p, integr_p_err, capsize=3, capthick=0.4, ecolor="black", elinewidth=0.4,
                fmt='none')

    ax.set_ylabel(r"$p$ [Torr]")
    ax.set_xlabel("amu")
    if max(atom) < 50:
        ax.set_xticks(np.arange(0, max(atom) + 1, 2))
    if max(atom) > 50:
        ax.set_xticks(np.arange(0, max(atom) + 1, 10))
    plt.savefig("Report/DataResultsPlots/{}.pdf".format(plot_name))
    plt.show()


if __name__ == "__main__":
    print(get_nist_peaks("oxygen", 2).head())
