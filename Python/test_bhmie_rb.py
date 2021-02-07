"""
La routine test_bhmie_rb.py trace l'intensité sortie de la subroutine 
bhmie en fonction de l'angle de rétrodiffusion pour mettre en avant
le phénomène d'arc-en-ciel.
"""
import numpy as np
import matplotlib.pyplot as plt

data=open('Data/Ivtheta_test_bhmie_rb.dat','r')
tab=np.array([s.strip().split() for s in data],dtype='f')

theta=np.linspace(0,180,len(tab[0,:]))

plt.xlabel(r'$\gamma\:(°)$')
plt.ylabel(r'$I_s\:(I_ir^{-2})$')
plt.plot(theta,tab[0,::-1],color='red')
plt.plot(theta,tab[1,::-1],color='green')
plt.plot(theta,tab[2,::-1],color='blue')
plt.xlim([35,60])
plt.ylim([0,8*10**6])
plt.grid()
plt.show()

