"""
La routine diag_pol_mie.py affiche le diagramme dipolaire de la diffusion
de Mie.
"""
import numpy as np
import matplotlib.pyplot as plt

color=['black']

data1=open('Data/theta_mie_diff.dat')
tab1=np.array([s.strip().split() for s in data1],dtype='f')
print(tab1.shape)
data_list=(tab1,)

def diag_pol(data_list):
    fig, ax= plt.subplots(figsize=(6.5,9))
    for i in range(0,len(data_list)):
        tab=data_list[i]

        diag_pol=np.zeros((tab.shape[0],2))
        diag_pol[:,0]=tab[:,1]*np.cos(tab[:,5])
        diag_pol[:,1]=tab[:,1]*np.sin(tab[:,5])

        ax.plot(diag_pol[:,0],diag_pol[:,1],color=color[i])
        ax.plot(diag_pol[:,0],-diag_pol[:,1],linestyle='--',color=color[i])

    ax.text(-0.95,0.90,r'$\tau_{max}=1$ '+'et'+r'$x=0.01$')
    ax.set_xlabel(r'$P_{esc}cos(\theta)$')
    ax.set_ylabel(r'$P_{esc}sin(\theta)$')
    ax.set_ylim([-1,1])
    ax.set_xlim([-1,1])
    ax.grid()
    
    plt.show()
diag_pol(data_list)
