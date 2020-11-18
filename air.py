from choose_sequences import sequences
import numpy as np
import matplotlib.pyplot as plt
relative=False
combine=5 
start=0
end = 50

o2 = np.arange(8,dtype= np.float)
co2 = np.arange(8,dtype= np.float)
breath= np.arange(8)
amu0,p0,err0= sequences('air_background.csv',False,combine,start,end,relative,False)
o2[0]=p0[32-1]
co2[0]=p0[44-1]
for i in np.arange(1,8):
    amu, p, err =sequences("air{}.csv".format(i),False,combine,start,end,relative,False)
    p = p*p0[int(5*28/combine-1)]/p[int(5*28/combine-1)]
    o2[i]=p[32-1]
    co2[i]=p[44-1]
#plt.plot(amu,p,label='breathing {}'.format(i))
#plt.errorbar(amu, p, err, capsize=3, capthick=0.4 ,ecolor="black", elinewidth=0.4 ,fmt ='none')
#
plt.plot(breath,o2,label=r'$O_2$')
plt.plot(breath,co2,label=r'$CO_2$')
plt.legend()
if relative:
    plt.ylabel(r"$p_{part}$ / $p_{tot}$ [%]")
else: plt.ylabel(r"$p$ [Torr]")
#plt.xlabel('amu')
plt.xlabel('breath')
plt.savefig("Report/DataResultsPlots/air.pdf")
plt.show()