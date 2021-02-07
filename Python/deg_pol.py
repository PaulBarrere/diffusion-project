"""
La routine deg_pol.py trace le degré de polarisation en fonction
de l'angle de diffusion.
"""
import numpy as np
import matplotlib.pyplot as plt

degpol='mie_low_tmax'

if (degpol=='mie_high_tmax'):
    data1=open('Data/theta_mie_diff_1.dat')
    data2=open('Data/theta_mie_diff_2.dat')
    data3=open('Data/theta_mie_diff_5.dat')
    data4=open('Data/theta_mie_diff_10.dat')
    tau=np.array([1,2,5,10])
    #N=np.array([10**7,10**6,10**7,5*10**6])
elif (degpol=='mie_low_tmax'):
    data1=open('Data/theta_mie_diff_005.dat')
    data2=open('Data/theta_mie_diff_001.dat')
    data3=open('Data/theta_mie_diff_0005.dat')
    data4=open('Data/theta_mie_diff_0001.dat')
    tau=np.array([0.05,0.01,0.005,0.001])
elif (degpol=='mie_CS_tmax'):
    data1=open('Data/theta_mie_diff_CloudySky_05.dat')
    data2=open('Data/theta_mie_diff_CloudySky_1.dat')
    data3=open('Data/theta_mie_diff_CloudySky_3.dat')
    data4=open('Data/theta_mie_diff_CloudySky_5.dat')
    tau=np.array([0.5,1,3,5])

tab1=np.array([s.strip().split() for s in data1],dtype='f')
tab2=np.array([s.strip().split() for s in data2],dtype='f')
tab3=np.array([s.strip().split() for s in data3],dtype='f')
tab4=np.array([s.strip().split() for s in data4],dtype='f')
data_list=(tab1,tab2,tab3,tab4)

if (degpol=='mie'):
    data1=open('Data/theta_mie_diff_CloudySky.dat')
elif (degpol=='mie_CS'):
    data1=open('Data/theta_mie_diff.dat')

tab1=np.array([s.strip().split() for s in data1],dtype='f')
data_list=(tab1,)

def deg_pol(data_list):
    fig, ax= plt.subplots(figsize=(6.5,9))
    for i in range(0,len(data_list)):
        tab=data_list[i]
        P=np.sqrt(tab[:,2]*tab[:,2]+tab[:,3]*tab[:,3]+tab[:,4]*tab[:,4])/tab[:,1]
        ax.plot(tab[:,5],P)#,label=r'$\tau=$'+"{:.0e}".format(tau[i]))#' ,$N_{part}=$'+"{:.0e}".format(N[i]))

    P_th=(1-np.cos(tab[:-1,5])**2)/(1+np.cos(tab[:-1,5])**2)
    ax.plot(tab[:-1,5],P_th,label='P théorique')
    ax.set_xlabel(r'$\theta$')
    ax.set_ylabel(r'P')
    ax.grid()

    plt.legend()
    plt.show()
deg_pol(data_list)
