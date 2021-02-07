"""
La routine test_bhmie.py trace la dépendance de l'intensité sortie
de la subroutine bhmie en fonction du rayon 'a' de la particule ou
de la longueur d'onde 'lambda'.
"""
import numpy as np
import matplotlib.pyplot as plt

dpd='a'

if (dpd=='lambda'):
    data1=open('Data/Ivlambda_rayl_test_bhmie.dat','r')
    tab1=np.array([s.strip().split() for s in data1],dtype='f')

    print(tab1.shape)
    lamb=np.zeros(len(tab1[:,0]))
    for i in range(len(tab1[:,0])):
        lamb[i]=1.25*2**i

    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel(r'$\lambda\:(\mu m)$')
    plt.ylabel(r'$I_s\:(I_ir^{-2})$')
    plt.title(r'I vs $\lambda$ ($a=2\:\mu m$)')
    plt.axvline(x=125,color='black',linestyle='--')
    plt.plot(lamb,tab1[:,0],label=r'$\theta=0$')
    plt.plot(lamb,tab1[:,50],label=r'$\theta=\pi/4$')
    plt.plot(lamb,tab1[:,100],label=r'$\theta=\pi/2$')
    plt.plot(lamb,tab1[:,150],label=r'$\theta=3\pi/4$')
    plt.plot(lamb,tab1[:,-1],label=r'$\theta=\pi$')
    plt.legend()
    plt.grid()
    plt.show()   

elif (dpd=='a'):
    data2=open('Data/Iva_rayl_test_bhmie.dat','r')
    tab2=np.array([s.strip().split() for s in data2],dtype='f')

    print(tab2.shape)
    lamb=np.zeros(len(tab2[:,0]))
    for i in range(len(tab2[:,0])):
        lamb[i]=7.95*10**(-6)*2**i

    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel(r'$a$ ($\mu m$)')
    plt.ylabel(r'$I_s\:(I_ir^{-2})$')
    plt.title(r'I vs $a$ ($\lambda=0.5\:\mu m$)')
    plt.axvline(x=7.95*10**(-3),color='black',linestyle='--')
    plt.plot(lamb,tab2[:,0],label=r'$\theta=0$')
    plt.plot(lamb,tab2[:,50],label=r'$\theta=\pi/4$')
    plt.plot(lamb,tab2[:,100],label=r'$\theta=\pi/2$')
    plt.plot(lamb,tab2[:,150],label=r'$\theta=3\pi/4$')
    plt.plot(lamb,tab2[:,-1],label=r'$\theta=\pi$')
    plt.legend()
    plt.grid()
    plt.show()

