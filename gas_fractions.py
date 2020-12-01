import numpy as np
from nist import get_nist_peaks, nist_aprox
from gas_analysis import gas_analysis
import csv

def Fraction_Table(comp, p, err, exp_comp,measured):
    exp_comp=np.array(exp_comp)
    with open('komp_{}.csv'.format(measured),'w', newline='') as csvfile:
        fieldnames = ['substance', 'fraction', 'err']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for i in exp_comp:
            #
            a =np.in1d(comp,i)#.reshape(comp.shape)
            b =np.in1d(comp,exp_comp)#.reshape(comp.shape)
            frac = np.round(p[a]/sum(p[b])*100,2)[0]
            frac_err = np.round(np.sqrt(( 1/sum(p[b])**2 - 2 * p[a] / sum(p[b])**3) * err[a]**2 + (p[a]**2 / sum(p[b])**4) * sum(err[b]**2)) *100,2)[0]
            writer.writerow({'substance': i, 'fraction': frac, 'err': frac_err})
    return exp_comp, frac, frac_err
  

deo = nist_aprox(*gas_analysis('deo.csv','airbaseline.csv'),'Deo')
exp_comp = ['ethanol','propane','butane','isobutane']
Fraction_Table(*deo,exp_comp,'deo')

residual = nist_aprox(*gas_analysis('xenonbaseline_highres.csv',False),'resiudal gas')
exp_comp = [r'H$_2$O',r'CO$_2$',r'N$_2$',r'O$_2$']
Fraction_Table(*residual,exp_comp,'residual')

nobel_gas = nist_aprox(*gas_analysis('mix2.csv','mix_baseline2.csv'),'mix')
exp_comp = ['Ar','Xe','Kr']
Fraction_Table(*nobel_gas,exp_comp,'nobel_gas')

air = nist_aprox(*gas_analysis('air00.csv','airbaseline.csv'),'air')
exp_comp = [r'N$_2$',r'O$_2$']
Fraction_Table(*air,exp_comp,'air_o2_vs_n2')
