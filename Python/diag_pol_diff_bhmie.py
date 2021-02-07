"""
La routine diag_pol_diff_bhmie.py compare la fonction de phase 
sortie de la subroutine bhmie et un type d'enveloppe dans un 
diagramme polaire.
"""
import numpy as np
import matplotlib.pyplot as plt

env='rayl' #Type d'enveloppe

data1=open('Data/bhmie_diff.dat','r')
if (env=='rayl'):
    data2=open('Data/rayl_diff.dat','r')
elif (env=='hg'):
    data2=open('Data/hg_diff.dat','r')
elif (env=='polyn'):
    data2=open('Data/polyn_diff.dat','r')
tab1=np.array([s.strip().split() for s in data1],dtype='f')
tab2=np.array([s.strip().split() for s in data2],dtype='f')
data_list=(tab1,tab2)

def diag_pol(data_list):
    fig, ax= plt.subplots(figsize=(6.5,9))
    for i in range(0,len(data_list)):
        tab=data_list[i]

        diag_pol=np.zeros((tab.shape[0],2))
        diag_pol[:,0]=-tab[:,1]*np.cos(tab[:,2])
        diag_pol[:,1]=tab[:,1]*np.sin(tab[:,2])

        ax.plot(diag_pol[:,0],diag_pol[:,1],color=color[i])
        ax.plot(diag_pol[:,0],-diag_pol[:,1],linestyle='--',color=color[i])

    ax.set_xlabel(r'$P_{esc}cos(\theta)$')
    ax.set_ylabel(r'$P_{esc}sin(\theta)$')
    ax.set_ylim([-1,1])
    ax.set_xlim([-1,1])
    ax.grid()
    
    plt.legend()
    plt.show()
diag_pol(data_list)
