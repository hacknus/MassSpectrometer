import matplotlib.pyplot as plt
import numpy as np

amu = np.array([1,2,3,4])
p_A = np.array([0,1.5,2,4])*10**-7
p_B = np.array([1,2,3,0])*10**-7
p = np.array([2,3,5,4])*10**-7

width = 0.6
fig, ax = plt.subplots(1, 1)
ax.bar(amu+width/4,p,width=width/2,color='black',label= 'measured')
ax.bar(amu-width/4,p_A+p_B,width=width/2,label = 'gas j')
ax.bar(amu-width/4,p_A,width=width/2, label = 'gas i')
ax.legend()
ax.set_xlabel('amu')
ax.set_xticklabels([r'A$_{max~b}-2$',r'A$_{max~b}-1$',r'A$_{max~b}$',r'A$_{max~a}$'])
ax.set_ylabel('partial pressure [Torr]' )
ax.set_xticks(amu)
ax.set_xlabel('mass number')
plt.savefig('Report\DataResultsPlots\ilustration_nist_approx1.pdf')
plt.show()

width = 0.6
fig, ax = plt.subplots(1, 1)
ax.bar(amu+width/4,p,width=width/2,color='black',label= 'measured')
ax.bar(amu-width/4,p,width=width/2,label = 'gas j')
ax.bar(amu-width/4,p_A/(p_A+p_B)*p,width=width/2, label = 'gas i')
ax.legend()
ax.set_xlabel('mass number')

ax.set_xticks(amu)
ax.set_xticklabels([r'A$_{max~b}-2$',r'A$_{max~b}-1$',r'A$_{max~b}$',r'A$_{max~a}$'])
ax.set_ylabel('partial pressure [Torr]' )
plt.savefig('Report\DataResultsPlots\ilustration_nist_approx2.pdf')
plt.show()