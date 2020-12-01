from gas_analysis import gas_analysis
import matplotlib.pyplot as plt
import numpy as np
import csv

#xenon with xenon data
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


s = sum(p_xe)
s_err = np.sqrt(sum(p_xe_err**2))
frac_xe = np.round(p_xe/sum(p_xe)*100,2)
frac_xe_err = np.round(np.sqrt(((1/s)**2  - p_xe/s**3 )* p_xe_err**2 + (p_xe/s**2)**2 * s_err**2)*100,2)


with open('isotop_xenon.csv','w', newline='') as csvfile:
    fieldnames = ['Isotop', 'fraction', 'err']
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()
    for i in np.arange(len(p_xe)):
        writer.writerow({'Isotop': r'Xe$^{}$'.format(xe_iso[i]), 'fraction': frac_xe[i], 'err': frac_xe_err[i]})


#krpyton mix data
amu, p, err = gas_analysis('mix2.csv','mix_baseline2.csv')
kr = [84,86]     
p_kr = []
p_kr_err = []
for i in kr:
    p_kr.append(*p[amu==i])
    p_kr_err.append(*err[amu==i])
    
p_kr = np.array(p_kr)
p_kr_err = np.array(p_kr_err)
s = sum(p_kr)
s_err = np.sqrt(sum(p_kr_err**2))

frac_kr = np.round(p_kr/sum(p_kr)*100,2)
frac_kr_err = np.round(np.sqrt(((1/s)**2  - p_kr/s**3 )* p_kr_err**2  + (p_kr/s**2)**2 * s_err**2)*100,2)


with open('isotop_krypton.csv','w', newline='') as csvfile:
    fieldnames = ['Isotop', 'fraction', 'err']
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()
    for i in np.arange(len(p_kr)):
        writer.writerow({'Isotop': r'Kr$^{}$'.format(kr[i]), 'fraction': frac_kr[i], 'err': frac_kr_err[i]})



#argon mix data
amu, p, err = gas_analysis('mix2.csv','mix_baseline2.csv')

ar_half_mix = np.array([20])
ar_iso_mix = np.array([40])

p_ar_iso_mix = p[amu==40]
p_ar_iso_err_mix = err[amu==40]  

p_ar_half_mix = p[amu==20]
p_ar_half_err_mix = err[amu==20]  

p_ar_mix = p_ar_iso_mix + p_ar_half_mix
p_ar_err_mix = np.sqrt(p_ar_iso_err_mix**2+p_ar_half_err_mix**2)

############################################################################ double Ionization plot 
#xenon with mix data
amu, p, err = gas_analysis('mix2.csv','mix_baseline2.csv')

xe_half_mix = np.array([66])
xe_iso_mix = np.array([132])

p_xe_iso_mix = p[amu==132]
p_xe_iso_err_mix = err[amu==132]  

p_xe_half_mix = p[amu==66]
p_xe_half_err_mix = err[amu==66]  

p_xe_mix = p_xe_iso_mix + p_xe_half_mix
p_xe_err_mix = np.sqrt(p_xe_iso_err_mix**2+p_xe_half_err_mix**2)


#argon with argon data
amu, p, err = gas_analysis('argon2.csv','argonbaseline.csv')

ar_half = np.array([20])
ar_iso = np.array([40])

p_ar_iso = p[amu==40]
p_ar_iso_err = err[amu==40]  

p_ar_half = p[amu==20]
p_ar_half_err = err[amu==20]  

p_ar = p_ar_iso + p_ar_half
p_ar_err = np.sqrt(p_ar_iso_err**2+p_ar_half_err**2)

#Doppeljonisierung in abhängigkeit von massenzahl:
dop_xe = p_xe_half/p_xe
dop_xe_err = np.sqrt((1/p_xe*p_xe_half_err)**2 + (p_xe_half/p_xe**2*p_xe_err)**2)
dop_xe_mix = p_xe_half_mix/p_xe_mix
dop_xe_err_mix = np.sqrt((1/p_xe_mix*p_xe_half_err_mix)**2 + (p_xe_half_mix/p_xe_mix**2*p_xe_err_mix)**2)

dop_ar_mix = p_ar_half_mix/p_ar_mix
dop_ar_err_mix = np.sqrt((1/p_ar_mix*p_ar_half_err_mix)**2 + (p_ar_half_mix/p_ar_mix**2*p_ar_err_mix)**2)
dop_ar = p_ar_half/p_ar
dop_ar_err = np.sqrt((1/p_ar*p_ar_half_err)**2 + (p_ar_half/p_ar**2*p_ar_err)**2)




plt.plot(dop_ar,ar_iso,label = 'gas mix')

dop_mix = np.concatenate((dop_ar_mix,dop_xe_mix))
dop_err_mix = np.concatenate((dop_ar_err_mix,dop_xe_err_mix))
dop_amu_mix = np.concatenate((ar_iso_mix,xe_iso_mix))

fig, ax = plt.subplots(1, 1)
ax.set_title('double ionization')
ax.plot(ar_iso, dop_ar*100,'o', label='Ar data')
ax.plot(xe_iso, dop_xe*100,'o', label='Xe data')
ax.plot(dop_amu_mix, dop_mix*100, 'o',label='mix data')
ax.legend()
plt.show()