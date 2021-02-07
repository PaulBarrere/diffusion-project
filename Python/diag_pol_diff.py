"""
La routine diag_pol_diff.py affiche les diagrammes polaires des
fonctions de phase.
"""
import numpy as np
import matplotlib.pyplot as plt

diff='rayl' #Type de diffusion

if (diff=='iso'):
    data1=open('Data/iso_diff.dat','r')
elif (diff=='rayl'):
    data1=open('Data/rayl_diff.dat','r')
elif (diff=='hg'):
    data1=open('Data/hg_diff.dat','r')
elif (diff=='polyn'):
    data1=open('Data/polyn_diff.dat','r')
elif (diff=='bhmie'):
    data1=open('Data/bhmie_diff.dat','r')
tab1=np.array([s.strip().split() for s in data1],dtype='f')
data_list=(tab1,)
color=['black']

def diag_pol(data_list):
    fig, ax= plt.subplots(figsize=(6.5,9))
    for i in range(0,len(data_list)):
        tab=data_list[i]

        diag_pol=np.zeros((tab.shape[0],2))
        diag_pol[:,0]=-tab[:,1]*np.cos(tab[:,2])
        diag_pol[:,1]=tab[:,1]*np.sin(tab[:,2])

        ax.plot(diag_pol[:,0],diag_pol[:,1],color=color[i])
        ax.plot(diag_pol[:,0],-diag_pol[:,1],linestyle='--',color=color[i])

    ax.text(-0.95,0.90,r'$\tau_{max}=0$')
    ax.set_xlabel(r'$P_{esc}cos(\theta)$')
    ax.set_ylabel(r'$P_{esc}sin(\theta)$')
    ax.set_ylim([-1,1])
    ax.set_xlim([-1,1])
    ax.grid()
    
    plt.show()
diag_pol(data_list)
