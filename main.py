## Import de librairies
import numpy as np
from matplotlib import pyplot as plt
import matplotlib
matplotlib.use('TkAgg')
from mpl_toolkits.mplot3d import Axes3D

################################################## Traitement du texte pour le pendule simple en 3D #####################################
fichier = open('res_1pend.txt','r')

L_T = [] 
L_Teta = []
L_Phi = []
Lignes = fichier.readlines()
for i in range(len(Lignes)):
    ligne = Lignes[i].split(' ')
    L_T.append(float(ligne[0]))
    L_Teta.append(float(ligne[1]))
    L_Phi.append(float(ligne[3]))

L_X = []
L_Y = []
L_Z = []

Rho = 2 #Longueur du pendule

for i in range(len(L_Teta)): #From spherical to cartesian
    Teta = L_Teta[i]
    Phi = L_Phi[i]
    x = Rho*np.cos(Phi)*np.sin(Teta)
    y = Rho*np.sin(Phi)*np.sin(Teta)
    z = Rho*np.cos(Teta)
    L_X.append(-x)
    L_Y.append(-y)
    L_Z.append(-z+Rho)


plt.figure('Affichage de la trajectoire du pendule')
axes = plt.axes(projection="3d")
axes.set_xlabel("X")
axes.set_ylabel("Y")
axes.set_zlabel("Z")
axes.set_xlim(-Rho,Rho)
axes.set_ylim(-Rho,Rho)
axes.set_zlim(0,Rho)
axes.plot(L_X,L_Y,L_Z,label = 'Position du pendule')
plt.show()

# For an animation of the simple pendulum, uncomment it
# for i in range(0,len(L_X),2):
#     plt.figure('Affichage des coordonnées sphériques dans le temps')
#     axes = plt.axes(projection="3d")
#     axes.set_xlabel("X")
#     axes.set_ylabel("Y")
#     axes.set_zlabel("Z")
#     axes.set_xlim(-Rho,Rho)
#     axes.set_ylim(-Rho,Rho)
#     axes.set_zlim(0,Rho)
#     axes.plot(L_X[:i],L_Y[:i],L_Z[:i],label = 'Position du pendule')
#     plt.pause(0.0001)
#     plt.clf()
# plt.show()

################################################# Traitement de texte pour le double pendule en 3D ########################################
plt.close('all')
fichier = open('res_2pend.txt','r')

l1,l2 = 1,1 # Longueurs respectives du premier et second pendule en m

L_T = []

L_Q1 = []
L_Q2 = []
L_Q3 = []
L_Q4 = []

Lignes = fichier.readlines()
for i in range(len(Lignes)):
    ligne = Lignes[i].split(' ')
    L_T.append(float(ligne[0]))

    L_Q1.append(float(ligne[1]))
    L_Q2.append(float(ligne[2]))
    L_Q3.append(float(ligne[3]))
    L_Q4.append(float(ligne[4]))

fichier.close()
    
#Passe les coordonnées stéréographiques en cartésiennes
A_Q1 = np.array(L_Q1)
A_Q2 = np.array(L_Q2)
A_Q3 = np.array(L_Q3)
A_Q4 = np.array(L_Q4)

""" plt.figure("Coordonnées stéréographiques des pendules")
plt.plot(A_Q1,A_Q2,label="Pendule 1")
plt.plot(A_Q3,A_Q4,label="Pendule 2")
plt.legend()
plt.show() """

#On récupère les coordonnées cartésiennes

L_X1 = []
L_Y1 = []
L_Z1 = []
L_X2 = []
L_Y2 = []
L_Z2 = []

fichier = open("cart_2pend.txt","r")
Lignes = fichier.readlines()
for i in range(len(Lignes)):    #Récupère les coordonnées cartésiennes
    ligne = Lignes[i].split(' ')

    L_X1.append(float(ligne[0]))
    L_Y1.append(float(ligne[1]))
    L_Z1.append(float(ligne[2]))
    L_X2.append(float(ligne[3]))
    L_Y2.append(float(ligne[4]))
    L_Z2.append(float(ligne[5]))
fichier.close()    

plt.figure('Affichage de la trajectoire des deux pendules')
axes = plt.axes(projection="3d")
axes.set_xlabel("X")
axes.set_ylabel("Y")
axes.set_zlabel("Z")
axes.set_xlim(-l1-l2,l1+l2)
axes.set_ylim(-l1-l2,l1+l2)
axes.set_zlim(0,2*l1+2*l2)
axes.plot([0],[0],[l1+l2],marker = 'o',color = 'black')
axes.plot(L_X1[:2000],L_Y1[:2000],L_Z1[:2000],label = '1st trajectory',color='yellow',marker = 'o')
axes.plot(L_X2[:2000],L_Y2[:2000],L_Z2[:2000],label = '2nd trajectory',color='cyan',marker = 'o')

axes.plot(L_X1[0],L_Y1[0],L_Z1[0],label = '1st end',color='orange',marker = 'o')
axes.plot(L_X2[0],L_Y2[0],L_Z2[0],label = '2nd end',color='blue',marker = 'o')
plt.legend()
plt.show() 

# Trajectory with a tail, uncomment for the animation
Int_ind = 3
A_Affichage = np.arange(0,len(L_X1),Int_ind)
Ind_trainée = 80
# for i in A_Affichage:
#     plt.figure("Trajectoire des deux pendules animé")
#     axes = plt.axes(projection="3d")
#     axes.set_xlabel("X")
#     axes.set_ylabel("Y")
#     axes.set_zlabel("Z")
#     axes.set_xlim(-l1-l2,l1+l2)
#     axes.set_ylim(-l1-l2,l1+l2)
#     axes.set_zlim(0,2*l1+2*l2)
#     axes.plot([L_X1[i]],[L_Y1[i]],[L_Z1[i]],marker = 'o',color = 'blue')
#     axes.plot([L_X2[i]],[L_Y2[i]],[L_Z2[i]],marker = 'o',color = 'red')
#     if i >= Ind_trainée:
#         axes.plot([0],[0],[l1+l2],marker = 'o',color = 'black',label = "Origine")

#         axes.plot(L_X1[i-Ind_trainée:i],L_Y1[i-Ind_trainée:i],L_Z1[i-Ind_trainée:i],color='blue',label="Trajectoire 1")
        
#         axes.plot(L_X2[i-Ind_trainée:i],L_Y2[i-Ind_trainée:i],L_Z2[i-Ind_trainée:i],color='red',label = "Trajectoire 2")
#     else:
#         axes.plot([0],[0],[l1+l2],marker = 'o',color = 'black',label = "Origine")

#         axes.plot(L_X1[:i],L_Y1[:i],L_Z1[:i],color='blue',label="Trajectoire 1")
        
#         axes.plot(L_X2[:i],L_Y2[:i],L_Z2[:i],color='red',label = "Trajectoire 2")
#     plt.legend()
#     plt.pause(0.0001)
#     plt.clf()
# plt.show()

plt.close('all')

## Etude de deux systèmes très proches afin de voir si le caractère chaotique est vérifié

""" Le premier système est parti avec q1 = 0.7,q2 = 0,q3 = 0.9,q4 = 0.9,p1 = 0,p2 = 0,p3 = 0,p4 = 0 """
""" Le second part avec q1 = 0.7,q2 = 0,q3 = 0.9,q4 = 0.89,p1 = 0,p2 = 0,p3 = 0,p4 = 0 """
# For an animation with 2 double pendulum, care you need to run 2 simulatio before
# L_X1_2 = []
# L_Y1_2 = []
# L_Z1_2 = []
# L_X2_2 = []
# L_Y2_2 = []
# L_Z2_2 = []

# fichier = open("Coord_cart_2.txt","r")
# Lignes = fichier.readlines()
# for i in range(len(Lignes)):
#     ligne = Lignes[i].split(' ')

#     L_X1_2.append(float(ligne[0]))
#     L_Y1_2.append(float(ligne[1]))
#     L_Z1_2.append(float(ligne[2]))
#     L_X2_2.append(float(ligne[3]))
#     L_Y2_2.append(float(ligne[4]))
#     L_Z2_2.append(float(ligne[5]))
# fichier.close()


# for i in A_Affichage:
#     plt.figure('Comparatison entre deux systèmes')
#     axes = plt.axes(projection="3d")
#     axes.set_facecolor("black")
#     axes.view_init(elev = 10,azim = i*np.pi/(10*Int_ind))
#     axes.set_xlabel("X")
#     axes.set_ylabel("Y")
#     axes.set_zlabel("Z")
#     axes.set_xlim(-l1-l2,l1+l2)
#     axes.set_ylim(-l1-l2,l1+l2)
#     axes.set_zlim(0,2*l1+2*l2)
#     axes.quiver([0,0,0],[0,0,0],[l1+l2,l1+l2,l1+l2],[0.5,0,0],[0,0.5,0],[0,0,0.5],color="Black")
#     axes.plot([0],[0],[l1+l2],marker = 'o',color = 'white')

#     axes.plot([L_X1[i]],[L_Y1[i]],[L_Z1[i]],marker = 'o',color = 'cyan',label = 'Système 1')
#     axes.plot([0,L_X1[i]],[0,L_Y1[i]],[l1+l2,L_Z1[i]],color='cyan')
    
#     axes.plot([L_X2[i]],[L_Y2[i]],[L_Z2[i]],marker = 'o',color = 'blue')
#     axes.plot([L_X1[i],L_X2[i]],[L_Y1[i],L_Y2[i]],[L_Z1[i],L_Z2[i]],color='blue')
    
#     axes.plot([L_X1_2[i]],[L_Y1_2[i]],[L_Z1_2[i]],marker = 'o',color = 'orange',label = 'Système 2')
#     axes.plot([0,L_X1_2[i]],[0,L_Y1_2[i]],[l1+l2,L_Z1_2[i]],color='orange')
    
#     axes.plot([L_X2_2[i]],[L_Y2_2[i]],[L_Z2_2[i]],marker = 'o',color = 'red')
#     axes.plot([L_X1_2[i],L_X2_2[i]],[L_Y1_2[i],L_Y2_2[i]],[L_Z1_2[i],L_Z2_2[i]],color='red')
#     plt.legend()
#     plt.savefig("F:\Mag 1\Image_info\image_"+str(i)+".png")
#     plt.pause(0.0001)
#     plt.clf()
# plt.show(block=False)

## Enregistrement et création de la vidéo selon des images triées dans un ordre numérique
""" import cv2
import glob
import re

img_array = []
numbers = re.compile(r'(\d+)')
def numericalSort(value):
    parts = numbers.split(value)
    parts[1::2] = map(int, parts[1::2])
    return parts

for filename in sorted(glob.glob("F:\\Mag 1\\Image_info\\image_*.png") , key=numericalSort):
    img = cv2.imread(filename)
    height, width, layers = img.shape
    size = (width,height)
    img_array.append(img)

out = cv2.VideoWriter('Comparaison_pendule.mp4',cv2.VideoWriter_fourcc(*'DIVX'), 30, size)

for i in range(len(img_array)):
    out.write(img_array[i])
out.release()
 """


