from matplotlib import pyplot as plt
import numpy as np
import matplotlib
matplotlib.use('TkAgg')

A_Ec = []
A_Ep = []
A_Diff = []
fichier = open("Valeur.txt",'r')
Lignes = fichier.readlines()
for i in range(2392):
    ligne = Lignes[i].split(' ')
    A_Ec.append(float(ligne[0]))
    A_Ep.append(float(ligne[1]))
    A_Diff.append(float(ligne[0])-float(ligne[1]))

plt.close('all')
T = [k for k in range(2392)]
plt.plot(T,A_Ec,label='dh*f')
plt.plot(T,A_Ep,label='h*df')
plt.plot(T,A_Diff,label='Soustraction')
plt.legend()
plt.show()