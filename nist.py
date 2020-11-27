from jcamp import JCAMP_reader
import pandas as pd
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
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

def scale_nist(atom, integr_p, integr_p_err, gas, number=12):
    nist = get_nist_peaks(gas, p_number=number)
    amu = np.array(nist.m)
    p = np.array(nist.y)
    amu_max = amu[max(p) == p]
    if max(amu_max) in atom:
        scaled = integr_p[atom == amu_max] / p[amu == amu_max] * p
        scaled_err = p/ p[amu == amu_max]*integr_p_err  
        return amu, scaled, scaled_err
    else:
        return amu, np.zeros(len(amu), int), np.zeros(len(amu), int)


def nist_aprox(atom, integr_p, integr_p_err, plot_name, ethanol=False):
    amu_o2, p_o2_scaled, p_o2_scaled_err  = scale_nist(atom, integr_p,integr_p_err, 'oxygen')
    p_o2_0 = p_o2_scaled
    amu_argon, p_argon_scaled ,p_argon_scaled_err = scale_nist(atom, integr_p, integr_p_err, 'argon')
    amu_xenon, p_xenon_scaled, p_xenon_scaled_err = scale_nist(atom, integr_p,  integr_p_err,'xenon')
    amu_krypton, p_krypton_scaled, p_krypton_scaled_err = scale_nist(atom, integr_p,  integr_p_err,'krypton')
    amu_h2o, p_h2o_scaled , p_h2o_scaled_err= scale_nist(atom, integr_p,  integr_p_err,'water')
    p_h2o_0 = p_h2o_scaled
    amu_h2, p_h2_scaled, ph2_scaled_err = scale_nist(atom, integr_p, integr_p_err, 'hydrogen')
    amu_ethanol, p_ethanol_scaled, p_ethanol_scaled_err= scale_nist(atom, integr_p,  integr_p_err,'ethanol')
    p_ethanol_0 = p_ethanol_scaled

# cupled par (propane and butane)
    n_propane = get_nist_peaks("propane", p_number=6)
    amu_propane = np.array(n_propane.m)
    p_propane = np.array(n_propane.y)

    
    n_butane = get_nist_peaks("butane", p_number=6)
    amu_butane = np.array(n_butane.m)
    p_butane = np.array(n_butane.y)
    

    if 29 and 43 in atom:
        for i in np.arange(len(atom)):
            if atom[i] == 29:
                p29=integr_p[i]
                p29_err = integr_p_err[i]
            if atom[i] == 43:
                p43=integr_p[i]
                p43_err = integr_p_err[i]
        for i in np.arange(len(amu_butane)):
            if amu_butane[i] == 43:
                but43 = p_butane[i]
            if amu_butane[i] == 29:
                but29 = p_butane[i]
        for i in np.arange(len(amu_propane)):
            if amu_propane[i] == 43:
                pro43 = p_propane[i]
            if amu_propane[i] == 29:
                pro29 = p_propane[i]

        k_propane = (p29 - p43 / but43 * but29) / (pro29 - pro43/but43*but29) 
        k_propane_err = np.sqrt((p29_err / (pro29 - pro43/but43*but29))**2 + ((p43_err/ but43 * but29) / (pro29 - pro43/but43*but29))**2) 
        k_butane = (p43 - k_propane * pro43)/ but43     
        k_butane_err = np.sqrt((p43_err / but43)**2 + ((k_propane_err * pro43)/but43)**2)
        p_propane_scaled = k_propane * p_propane
        p_propane_scaled_err = k_propane_err * p_propane
        p_butane_scaled = k_butane * p_butane
        p_butane_scaled_err = k_butane_err * p_butane
        

    

    else:
        p_propane_scaled = np.zeros(len(amu_propane), int)
        p_butane_scaled = np.zeros(len(amu_butane), int)
        p_propane_scaled_err = np.zeros(len(amu_propane), int)
        p_butane_scaled_err = np.zeros(len(amu_butane), int)
    
    if ethanol:
        p_propane_scaled = np.zeros(len(amu_propane), int)
        p_butane_scaled = np.zeros(len(amu_butane), int)
        p_propane_scaled_err = np.zeros(len(amu_propane), int)
        p_butane_scaled_err = np.zeros(len(amu_butane), int)
    p_butane_0 = p_butane_scaled
    p_propane_0 = p_propane_scaled

#dependent maxima for co2
    n_co2 = get_nist_peaks("co2", p_number=6)
    amu_co2 = np.array(n_co2.m)
    p_co2 = np.array(n_co2.y)
    if 44 in atom:
        p_co2_scaled = (integr_p[atom == 44] - p_propane_scaled[amu_propane == 44]) / p_co2[amu_co2 == 44] * p_co2
        p_co2_scaled_err = np.sqrt((integr_p_err[atom == 44]/p_co2[amu_co2 == 44] * p_co2)**2 + (p_propane_scaled_err[amu_propane == 44] / p_co2[amu_co2 == 44] * p_co2)**2)
    else:
        p_co2_scaled = np.zeros(len(amu_co2), int)
    p_co2_0 = p_co2_scaled

#dependent maxima for nitrogen 
    n_n2 = get_nist_peaks("nitrogen", p_number=6)
    amu_n2 = np.array(n_n2.m)
    p_n2 = np.array(n_n2.y)
    if 28 in atom:
        p_n2_scaled = (integr_p[atom == 28] - p_co2_scaled[amu_co2 == 28] - p_butane_scaled[amu_butane == 28] -
                       p_propane_scaled[amu_propane == 28] - p_ethanol_scaled[amu_ethanol == 28]) / p_n2[amu_n2 == 28] * p_n2
        p_n2_scaled_err = np.sqrt(integr_p_err[atom == 28]**2 + p_co2_scale_errd[amu_co2 == 28]**2 + p_butane_scaled_err[amu_butane == 28]**2 + p_propane_scaled_err[amu_propane == 28]**2 + p_ethanol_scaled[amu_ethanol == 28]**2) / p_n2[amu_n2 == 28] * p_n2
    else:
        p_n2_scaled = np.zeros(len(amu_n2), int)
        p_n2_scaled_err = np.zeros(len(amu_n2), int)
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

#delet half masses:
    if 64.5 in atom:
        integr_p[atom==64.5] = 0
        integr_p_err[atom==64.5] = 0
    if 64.5 in atom:
        integr_p[atom==65.5] = 0
        integr_p_err[atom==65.5] = 0

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
    if max(atom) < 60:
        ax.set_xticks(np.arange(0, max(atom) + 1, 5))
    if max(atom) > 60:
        ax.set_xticks(np.arange(0, max(atom) + 1, 10))
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(1))
    
    plt.savefig("Report/DataResultsPlots/{}.pdf".format(plot_name))
    plt.show()


if __name__ == "__main__":
    print(get_nist_peaks("oxygen", 2).head())
