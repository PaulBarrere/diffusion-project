"""
La routine f_phase.py permet de vérifier le fit fait par la subroutine
dgels de la libraire LAPACK.
"""
import numpy as np
import matplotlib.pyplot as plt

data1=open('Data/rayl_diff.dat')
tab1=np.array([s.strip().split() for s in data1],dtype='f')

# On utilise le nombre de photons dans le cas où on utilise les 
# fonctions de phase pour la projection stéréographique RGB.
#N=open('Data/nb_photons_RGB_mie.dat')
#tab_N=np.array([s.strip().split() for s in N],dtype='f')[0,:]

coeff=open('Data/coefficients')
tab_coeff=np.array([s.strip().split() for s in coeff],dtype='f')[0,:]
tab2=np.zeros(len(tab1[:,0]))

for i in range(0,len(tab_coeff)):
    tab2+=tab_coeff[i]*np.cos(tab1[:,2])**i

#plt.figure(figsize=(6.5,9))
plt.xlabel(r'$\theta$')
plt.ylabel('P')
plt.plot(tab1[1:,2],tab2[1:],label=r'Fit',color='green')
plt.scatter(tab1[1:,2],tab1[1:,1],label=r'Données MC',marker='+')
plt.legend()
plt.grid()
plt.show()    
