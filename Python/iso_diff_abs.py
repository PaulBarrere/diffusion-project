"""
La routine iso_diff_abs.py trace la réflectance de la diffusion isotrope
en fonction de l'albédo.
"""
import numpy as np
import matplotlib.pyplot as plt

data1=open('iso_diff_abs.dat','r')
tab=np.array([s.strip().split() for s in data1],dtype='f')
print(tab.shape)

fig, ax= plt.subplots(figsize=(5,5))
ax.plot(tab[:,0],tab[:,1])
ax.set_xlabel(r'$a$')
ax.set_ylabel(r'$R$')
ax.grid()
plt.show()

