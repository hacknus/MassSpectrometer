from choose_sequences import sequences
import numpy as np
import matplotlib.pyplot as plt
relative=True
combine=5 
start=0
end = 50

amu0,p0,err0= sequences('air2.csv',False,combine,start,end,relative,False)
for i in [3,4,6,7]:
    amu, p, err =sequences("air{}.csv".format(i),False,combine,start,end,relative,False)
    p = p*p0[int(5*28/combine-1)]/p[int(5*28/combine-1)]-p0
    plt.plot(amu,p,label='breathing {}'.format(i))
    plt.errorbar(amu, p, err, capsize=3, capthick=0.4 ,ecolor="black", elinewidth=0.4 ,fmt ='none')
plt.legend()
if relative:
    plt.ylabel(r"$p_{part}$ / $p_{tot}$ [%]")
else: plt.ylabel(r"$p$ [Torr]")
plt.xlabel('amu')
plt.savefig("Report/DataResultsPlots/air.pdf")
plt.plot()