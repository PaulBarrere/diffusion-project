"""
La routine daig_pol_mie_RGB.py affiche les diagrammes polaires 
des fonction de phase de la diffusion de Mie dans les couleurs RGB.
"""
import numpy as np
import matplotlib.pyplot as plt

color=['red','green','blue']

nb=open("Data/nb_photons_RGB_rayl.dat",'r')
data1=open('Data/theta_mie_diff_R.dat','r')
data2=open('Data/theta_mie_diff_G.dat','r')
data3=open('Data/theta_mie_diff_B.dat','r')
tab_nb=np.array([s.strip().split() for s in nb],dtype='f')
tab1=np.array([s.strip().split() for s in data1],dtype='f')
tab2=np.array([s.strip().split() for s in data2],dtype='f')
tab3=np.array([s.strip().split() for s in data3],dtype='f')
data_list=(tab1,tab2,tab3)

def diag_pol(data_list):
    fig, ax= plt.subplots(figsize=(6.5,9))
    for i in range(0,len(data_list)):
        tab=data_list[i]

        diag_pol=np.zeros((tab.shape[0],2))
        diag_pol[:,0]=-tab[:,1]*(tab_nb[0,i]/tab_nb[0,-1])*np.cos(tab[:,5])
        diag_pol[:,1]=tab[:,1]*(tab_nb[0,i]/tab_nb[0,-1])*np.sin(tab[:,5])

        ax.plot(diag_pol[:,0],diag_pol[:,1],color=color[i])
        ax.plot(diag_pol[:,0],-diag_pol[:,1],linestyle='--',color=color[i])

    ax.text(-0.95,0.90,r'$D=0.4$')
    ax.set_xlabel(r'$P_{esc}cos(\theta)$')
    ax.set_ylabel(r'$P_{esc}sin(\theta)$')
    ax.set_ylim([-1,1])
    ax.set_xlim([-1,1])
    ax.grid()
    
    plt.legend()
    plt.show()
diag_pol(data_list)
