"""
La routine env_eff.py trace l'eficacité des enveloppes en fonction du facteur de forme.
"""
import numpy as np
import matplotlib.pyplot as plt

color=['red','green','blue']
label=['Rayleigh','Henyey-Greenstein','Poly3']

def efficacite(data_list):
    fig, ax= plt.subplots(figsize=(6.5,9))
    for i in range(0,len(data_list)):
        tab=data_list[i]
        ax.plot(tab[:,0],tab[:,1],color=color[i],label=label[i])

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel(r'$x$')
    ax.set_ylabel(r'Efficacité')
    ax.grid()
    
    plt.legend()
    plt.show()

data1=open('Data/env_eff_rayl.dat','r')
data2=open('Data/env_eff_hg.dat','r')
data3=open('Data/env_eff_polyn.dat','r')
tab1=np.array([s.strip().split() for s in data1],dtype='f')
tab2=np.array([s.strip().split() for s in data2],dtype='f')
tab3=np.array([s.strip().split() for s in data3],dtype='f')
data_list=(tab1,tab2,tab3)

efficacite(data_list)
