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
        scaled_err = integr_p_err[atom == amu_max]/ p[amu == amu_max]*p  
        return amu, scaled, scaled_err
    else:
        return amu, np.zeros(len(amu), float), np.zeros(len(amu), float)


n_isobutane = get_nist_peaks("isobutane", p_number=12)
amu_isobutane = np.array(n_isobutane.m)
p_isobutane = np.array(n_isobutane.y)
#print('isobutane')
#print(n_isobutane)

n_isopropanol = get_nist_peaks("isopropanol", p_number=12)
amu_isopropanol = np.array(n_isopropanol.m)
p_isopropanol = np.array(n_isopropanol.y)
#print('isobropanol')
#print(n_isopropanol)

n_butane = get_nist_peaks("butane", p_number=12)
#print('butane')
#print(n_butane)
n_propane = get_nist_peaks("propane", p_number=12)
#print('propane')
#print(n_propane)
n_n2 = get_nist_peaks("nitrogen", p_number=12)
#print(n_n2)
n_ethanol = get_nist_peaks("ethanol", p_number=12)
#print('ethanol')
#print(n_ethanol)
#print(get_nist_peaks('oxygen',12))
print('acetanaldeyde')
print(get_nist_peaks('acetaldehyde',12))



def nist_aprox(atom, integr_p, integr_p_err, plot_name, ethanol=False):
    amu_o2, p_o2_scaled, p_o2_scaled_err  = scale_nist(atom, integr_p,integr_p_err, 'oxygen')
    p_o2_0 = p_o2_scaled
    amu_argon, p_argon_scaled ,p_argon_scaled_err = scale_nist(atom, integr_p, integr_p_err, 'argon')
    amu_xenon, p_xenon_scaled, p_xenon_scaled_err = scale_nist(atom, integr_p,  integr_p_err,'xenon')
    amu_krypton, p_krypton_scaled, p_krypton_scaled_err = scale_nist(atom, integr_p,  integr_p_err,'krypton')
    amu_h2o, p_h2o_scaled , p_h2o_scaled_err= scale_nist(atom, integr_p,  integr_p_err,'water')
    p_h2o_0 = p_h2o_scaled
    amu_h2, p_h2_scaled, p_h2_scaled_err = scale_nist(atom, integr_p, integr_p_err, 'hydrogen')
    amu_ethanol, p_ethanol_scaled, p_ethanol_scaled_err= scale_nist(atom, integr_p,  integr_p_err,'ethanol')
    p_ethanol_0 = sum(p_ethanol_scaled)

# cupled par (propane and butane)
    n_propane = get_nist_peaks("propane", p_number=7)
    amu_propane = np.array(n_propane.m)
    p_propane = np.array(n_propane.y)

    
    n_butane = get_nist_peaks("butane", p_number=6)
    amu_butane = np.array(n_butane.m)
    p_butane = np.array(n_butane.y)
    
    n_isobutane = get_nist_peaks("isobutane", p_number=12)
    amu_isobutane = np.array(n_isobutane.m)
    p_isobutane = np.array(n_isobutane.y)
    
    n_isobutanol = get_nist_peaks("isopropanol", p_number=12)
    amu_isopropanol = np.array(n_isopropanol.m)
    p_isopropanol = np.array(n_isopropanol.y)
    
    

    if 29 and 43 and 41 in atom and ethanol != True:
        for i in np.arange(len(amu_ethanol)):
            if amu_ethanol[i] == 29:
                eth29 = p_ethanol_scaled[i]
                eth29_err = p_ethanol_scaled_err[i]
            if amu_ethanol[i] == 43:
                eth43 = p_ethanol_scaled[i]
                eth43_err = p_ethanol_scaled_err[i]
        for i in np.arange(len(atom)):
            if atom[i] == 29:
                p29=integr_p[i] - eth29
                p29_err = np.sqrt(integr_p_err[i]**2 + eth29_err**2)
            if atom[i] == 43:
                p43=integr_p[i] - eth43
                p43_err = np.sqrt(integr_p_err[i]**2 + eth43**2)
            if atom[i] == 41:
                p41=integr_p[i]
                p41_err = integr_p_err[i]
            
        for i in np.arange(len(amu_butane)):
            if amu_butane[i] == 43:
                but43 = p_butane[i]
            if amu_butane[i] == 29:
                but29 = p_butane[i]
            if amu_butane[i] == 41:
                but41 = p_butane[i]
            
                
        for i in np.arange(len(amu_propane)):
            if amu_propane[i] == 43:
                pro43 = p_propane[i]
            if amu_propane[i] == 29:
                pro29 = p_propane[i]
            if amu_propane[i] == 41:
                pro41 = p_propane[i]
        
        for i in np.arange(len(amu_isobutane)):
            if amu_isobutane[i] == 43:
                iso43 = p_isobutane[i]
            if amu_isobutane[i] == 29:
                iso29 = p_isobutane[i]
            if amu_isobutane[i] == 41:
                iso41 = p_isobutane[i]
            
        
        div = (-but29 * iso41 * pro43 + but29 * iso43 * pro41 + but41 * iso29 * pro43 - but41 * iso43 * pro29 - but43 * iso29 * pro41 + but43 * iso41 * pro29)         
        
        k_propane = -(but29 * iso41 * p43 - but29 * iso43 * p41 - but41 * iso29 * p43 + but41 * iso43 * p29 + but43 * iso29 * p41 - but43 * iso41 *  p29)/div
        k_propane_err = np.sqrt((but29 * iso41 * p43_err / div)**2 +  (but29 * iso43 * p41_err / div)**2 + ( but41 * iso29 * p43_err/div)**2 + ( but41 * iso43 * p29_err/div)**2 +  (but43 * iso29 * p41_err/div)**2 + (but43 * iso41 *  p29_err/div)**2)
        k_butane = -(iso29 * p43 * pro41 - iso29 * p41 * pro43 - iso41 * p43 * pro29 + iso41 *  p29 * pro43 + iso43 * p41 * pro29 - iso43 *  p29 * pro41)/div
        k_butane_err = np.sqrt((iso29 * p43_err * pro41/div)**2 + (iso29 * p41_err * pro43/div)**2 + (iso41 * p43_err * pro29/div)**2 + (iso41 *  p29_err * pro43/div)**2 + (iso43 * p41_err * pro29/div)**2 + (iso43 *  p29_err * pro41/div)**2)
        
        
        #k_iso = -(-but29 * p43 * pro41 + but29 * p41 * p3 + b2 * p43 * p1 - b2 * p29 * p3 - b3 *p41* p1 + b3 * p29 * pro41)
        k_isobutane = -(-but29 * p43 * pro41 + but29 * p41 * pro43 + but41 * p43 * pro29 - but41 * p29 * pro43 - but43 * p41 * pro29 + but43 * p29 * pro41)/div
        k_isobutane_err = np.sqrt((but29 * p43_err * pro41/div)**2 + (but29 * p41_err * pro43/div)**2 + (but41 * p43_err * pro29/div)**2 + (but41 *  p29_err * pro43/div)**2 + (but43 * p41_err * pro29/div)**2 + (but43 *  p29_err * pro41/div)**2)
        
        
            

        #k_propane = (p29 - p43 / but43 * but29) / (pro29 - pro43/but43*but29) 
        #k_propane_err = np.sqrt((p29_err / (pro29 - pro43/but43*but29))**2 + ((p43_err/ but43 * but29) / (pro29 - pro43/but43*but29))**2) 
        #k_butane = (p43 - k_propane * pro43)/ but43     
        #k_butane_err = np.sqrt((p43_err / but43)**2 + ((k_propane_err * pro43)/but43)**2)
        p_butane_scaled = k_butane * p_butane
        p_butane_scaled_err = k_butane_err * p_butane
        p_propane_scaled = k_propane * p_propane
        p_propane_scaled_err = k_propane_err * p_propane
        p_isobutane_scaled = k_isobutane * p_isobutane
        p_isobutane_scaled_err = k_isobutane_err * p_isobutane
        
        

    

    else:
        p_propane_scaled = np.zeros(len(amu_propane), float)
        p_butane_scaled = np.zeros(len(amu_butane), float)
        p_propane_scaled_err = np.zeros(len(amu_propane), float)
        p_butane_scaled_err = np.zeros(len(amu_butane), float)
        p_isobutane_scaled = np.zeros(len(amu_isobutane), float)
        p_isobutane_scaled_err = np.zeros(len(amu_isobutane), float)
    
#    if ethanol:
#        p_propane_scaled = np.zeros(len(amu_propane), int)
#        p_butane_scaled = np.zeros(len(amu_butane), int)
#        p_propane_scaled_err = np.zeros(len(amu_propane), int)
#        p_butane_scaled_err = np.zeros(len(amu_butane), int)
#        p_isobutane_scaled = np.zeros(len(amu_isobutane), int)
#        p_isobutane_scaled_err = np.zeros(len(amu_isobutane), int)
    
#    if ethanol:
#        p_propane_scaled = np.zeros(len(amu_propane), int)
#        p_butane_scaled = np.zeros(len(amu_butane), int)
#        p_propane_scaled_err = np.zeros(len(amu_propane), int)
#        p_butane_scaled_err = np.zeros(len(amu_butane), int)
    p_butane_0 = sum(p_butane_scaled)
    p_isobutane_0 = sum(p_isobutane_scaled)
    p_propane_0 = sum(p_propane_scaled)
    

#dependent maxima for co2
    n_co2 = get_nist_peaks("co2", p_number=6)
    amu_co2 = np.array(n_co2.m)
    p_co2 = np.array(n_co2.y)
    if 44 in atom:
        p_co2_scaled = (integr_p[atom == 44] - p_propane_scaled[amu_propane == 44] - p_isobutane_scaled[amu_isobutane == 44]) / p_co2[amu_co2 == 44] * p_co2
        p_co2_scaled_err = np.sqrt((integr_p_err[atom == 44]/p_co2[amu_co2 == 44] * p_co2)**2 + (p_propane_scaled_err[amu_propane == 44] / p_co2[amu_co2 == 44] * p_co2)**2 + (p_isobutane_scaled_err[amu_isobutane == 44] / p_co2[amu_co2 == 44] * p_co2)**2)
        if p_co2_scaled[0] <= 0:
            p_co2_scaled = p_co2_scaled * 0
            p_co2_scaled_err = p_co2_scaled * 0
            
    else:
        p_co2_scaled = np.zeros(len(amu_co2), float)
        p_co2_scaled_err = np.zeros(len(amu_co2), float)      
    p_co2_0 = p_co2_scaled

#dependent maxima for nitrogen 
    n_n2 = get_nist_peaks("nitrogen", p_number=2)
    amu_n2 = np.array(n_n2.m)
    p_n2 = np.array(n_n2.y)
    if 28 in atom:
        p_n2_scaled = (integr_p[atom == 28] - p_co2_scaled[amu_co2 == 28] - p_butane_scaled[amu_butane == 28] - p_ethanol_scaled[amu_ethanol == 28] -
                       p_propane_scaled[amu_propane == 28]  - p_isobutane_scaled[amu_isobutane == 28]) / p_n2[amu_n2 == 28] * p_n2
        p_n2_scaled_err = np.sqrt(integr_p_err[atom == 28]**2 + p_co2_scaled_err[amu_co2 == 28]**2 + p_butane_scaled_err[amu_butane == 28]**2 + p_propane_scaled_err[amu_propane == 28]**2 + p_ethanol_scaled[amu_ethanol == 28]**2 + p_isobutane_scaled[amu_isobutane == 28]**2) / p_n2[amu_n2 == 28] * p_n2
    else:
        p_n2_scaled = np.zeros(len(amu_n2), float)
        p_n2_scaled_err = np.zeros(len(amu_n2), float)
    p_n2_0 = p_n2_scaled
    
 # p  matrix:
    amu_gases = [amu_ethanol, amu_isobutane, amu_butane, amu_propane, amu_xenon, amu_argon, amu_krypton, amu_h2, amu_n2, amu_h2o, amu_o2, amu_co2]
    p_gases = [p_ethanol_scaled, p_isobutane_scaled, p_butane_scaled, p_propane_scaled, p_xenon_scaled, p_argon_scaled, p_krypton_scaled, p_h2_scaled, p_n2_scaled, p_h2o_scaled, p_o2_scaled, p_co2_scaled]
    
    p_gases_err = [p_ethanol_scaled_err, p_isobutane_scaled_err, p_butane_scaled_err, p_propane_scaled_err, p_xenon_scaled_err, p_argon_scaled_err, p_krypton_scaled_err, p_h2_scaled_err, p_n2_scaled_err, p_h2o_scaled_err, p_o2_scaled_err, p_co2_scaled_err]
    M = np.zeros((len(atom),len(amu_gases)),float)
    M_err = np.zeros((len(atom),len(amu_gases)),float)
    M_rel = np.zeros((len(atom),len(amu_gases)),float)
    M_rel_err = np.zeros((len(atom),len(amu_gases)),float)
    for amu_gas, p_gas, p_gas_err, j in zip(amu_gases, p_gases,p_gases_err, np.arange(len(amu_gases))):
        for amu, i in zip(atom, np.arange(len(atom))):
            if amu in amu_gas:
                M[i][j] = p_gas[amu_gas == amu]
                M_err[i][j] =  p_gas_err[amu_gas == amu]
            else:  
                M[i][j] = 0
                M_err[i][j] = 0

   # p relative matrix:
    for i in np.arange(len(M)):        
        if sum(M[i]) != 0:
            M_rel[i] = M[i]/sum(M[i])
            if 1 not in M_rel[i]:
                M_rel_err[i] = np.sqrt(( 1/sum(M[i])**2 - 2 * M[i] / sum(M[i])**3) * M_err[i]**2 + (M[i]**2 / sum(M[i])**4) * sum(M_err[i]**2))
    
    p_element = np.zeros(len(amu_gases)) 
    p_element_err = np.zeros(len(amu_gases)) 
    for i in np.arange(len(p_element)):
        p_element[i] = sum(integr_p * M_rel[:,i])
        p_element_err[i] = np.sqrt(sum((integr_p_err * M_rel[:,i])**2 + (integr_p * M_rel_err[:,i])**2))
    element = ['ethanol', 'isobutane', 'butane', 'propane', 'Xe', 'Ar', 'Kr', r'H$_2$', r'N$_2$', r'H$_2$O', r'O$_2$', r'CO$_2$']
    
     
    #plot :  make sum of the bars
    if ethanol != True:
        for i in amu_ethanol:
            if i in amu_isobutane:
                p_ethanol_scaled[amu_ethanol == i] += p_isobutane_scaled[amu_isobutane == i]
            if i in amu_propane:
                p_ethanol_scaled[amu_ethanol == i] += p_propane_scaled[amu_propane == i]
            if i in amu_butane:
                p_ethanol_scaled[amu_ethanol == i] += p_butane_scaled[amu_butane == i]
            if i in amu_h2o:
                p_ethanol_scaled[amu_ethanol == i] += p_h2o_scaled[amu_h2o == i]
            if i in amu_n2:
                p_ethanol_scaled[amu_ethanol == i] += p_n2_scaled[amu_n2 == i]
            if i in amu_o2:
                p_ethanol_scaled[amu_ethanol == i] += p_o2_scaled[amu_o2 == i]
            if i in amu_co2:
                p_ethanol_scaled[amu_ethanol == i] += p_co2_scaled[amu_co2 == i]
                
        for i in amu_isobutane:
            if i in amu_propane:
                p_isobutane_scaled[amu_isobutane == i] += p_propane_scaled[amu_propane == i]
            if i in amu_butane:
                p_isobutane_scaled[amu_isobutane == i] += p_butane_scaled[amu_butane == i]
            if i in amu_h2o:
                p_isobutane_scaled[amu_isobutane == i] += p_h2o_scaled[amu_h2o == i]
            if i in amu_n2:
                p_isobutane_scaled[amu_isobutane == i] += p_n2_scaled[amu_n2 == i]
            if i in amu_o2:
                p_isobutane_scaled[amu_isobutane == i] += p_o2_scaled[amu_o2 == i]
            if i in amu_co2:
                p_isobutane_scaled[amu_isobutane == i] += p_co2_scaled[amu_co2 == i]
        
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
        ax.bar(amu_h2 - width / 4, p_h2_scaled, width=width / 2, color="m", label=r'NIST H$_2$')
    if sum(p_krypton_scaled) != 0:
        ax.bar(amu_krypton - width / 4, p_krypton_scaled, width=width / 2, color="cyan", label=r'NIST Kr')
    if sum(p_argon_scaled) != 0:
        ax.bar(amu_argon - width / 4, p_argon_scaled, width=width / 2, color="steelblue", label=r'NIST Ar')
    if sum(p_xenon_scaled) != 0:
        ax.bar(amu_xenon - width / 4, p_xenon_scaled, width=width / 2, color="deepskyblue", label=r'NIST Xe')
    if p_ethanol_0 != 0:
        ax.bar(amu_ethanol - width / 4, p_ethanol_scaled, width=width / 2, color="deepskyblue", label=r'NIST ethanol')
    if p_isobutane_0 != 0:
        ax.bar(amu_isobutane - width / 4, p_isobutane_scaled, width=width / 2, color="purple", label=r'NIST isobutane')
    if p_propane_0 != 0:
        ax.bar(amu_propane - width / 4, p_propane_scaled, width=width / 2, color="pink", label=r'NIST propane')
    if p_butane_0 != 0:
        ax.bar(amu_butane - width / 4, p_butane_scaled, width=width / 2, color="orange", label=r'NIST butane')
    if sum(p_h2o_0) != 0:
        ax.bar(amu_h2o - width / 4, p_h2o_scaled, width=width / 2, color="red", label=r'NIST H$_2$O')
    if sum(p_n2_0) != 0:
        ax.bar(amu_n2 - width / 4, p_n2_scaled, width=width / 2, color="y", label=r'NIST N$_2$')
    if sum(p_o2_0) != 0:
        ax.bar(amu_o2 - width / 4, p_o2_scaled, width=width / 2, color="blue", label=r'NIST O$_2$')
    if sum(p_co2_scaled) != 0:
        ax.bar(amu_co2 - width / 4, p_co2_scaled, width=width / 2, color="green", label=r'NIST CO$_2$')
    ax.legend()
    ax.errorbar(atom + width / 4, integr_p, integr_p_err, capsize=3, capthick=0.4, ecolor="black", elinewidth=0.4,
                fmt='none')

    ax.set_ylabel(r"$p$ [Torr]")
    ax.set_xlabel("m/q [amu/e]")
    if max(atom) < 60:
        ax.set_xticks(np.arange(0, max(atom) + 1, 5))
    if max(atom) > 60:
        ax.set_xticks(np.arange(0, max(atom) + 1, 10))
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(1))
    
    plt.savefig("Report/DataResultsPlots/{}.pdf".format(plot_name))
    plt.show()
    
    return element, p_element, p_element_err 

if __name__ == "__main__":
    print(get_nist_peaks("oxygen", 2).head())
