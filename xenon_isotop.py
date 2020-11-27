from gas_analysis import gas_analysis
import matplotlib.pyplot as plt
import numpy as np
import csv

amu, p, err = gas_analysis('xenon_highres.csv','xenonbaseline_highres.csv')

xe_iso = [128,129,130,131,132,134,136]
xe_half = [64,64.5,65,65.5,66,67,68]

p_xe_iso = []
p_xe_iso_err = []
for i in xe_iso:
    p_xe_iso.append(*p[amu==i])
    p_xe_iso_err.append(*err[amu==i])

p_xe_half = []
p_xe_half_err = []
for i in xe_half:
    p_xe_half.append(*p[amu==i])
    p_xe_half_err.append(*err[amu==i])
    
p_xe =np.array(np.array(p_xe_iso) + np.array(p_xe_half))
p_xe_err = np.array(np.sqrt(np.array(p_xe_iso_err)**2 + np.array(p_xe_half_err)**2))

print(p_xe)
print(p_xe_err)

frac_xe = p_xe/sum(p_xe)

s = sum(p_xe)
s_err = np.sqrt(sum(p_xe_err**2))




frac_xe_err = np.sqrt((1/s)**2 * p_xe_err**2 + (p_xe/s**2)**2 * s_err**2)


print(frac_xe)
print(frac_xe_err)

with open('xenon_isotop.csv','w', newline='') as csvfile:
    fieldnames = ['Isotop', 'fraction []', 'err']
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()
    for i in np.arange(len(p_xe)):
        writer.writerow({'Isotop': '$Xe^{}$'.format(xe_iso[i]), 'fraction []': frac_xe[i], 'err': frac_xe_err[i]})