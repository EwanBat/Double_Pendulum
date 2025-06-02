// Import de librairies
#include<iostream>
#include<fstream>
#include<math.h>
#define M_PI           3.14159265358979323846
#include<stdlib.h>
using namespace std;

const double g = 9.81; // Constantante gravitationnelle terrestre en m.s-2
double m = 10; // Masse du point du pendule en kg
double Rho = 2; // Longueur du pendule en m 
int Nb_var = 4; // Nombre de variable pour la méthode de Runge-Kutta

// Création de la fonction rk4

void Ini_Mat(void (*f)(double*,double*,double),double*** Mat_Pond,double* Coord,double t,double dt,int Nb){
  int i,j,k;           // Fonction qui construit un tableau contenant les valeurs p1,p2,p3,p4 de la méthode rk4
  for (j=0;j<4;j++){ // Selon i, on fixe soit pj soit le terme dans le système qui permet de calculer pj
        for (i=0;i<2;i++){
            if (i==0 && j==0){
                for (k=0;k<Nb;k++){
		  Mat_Pond[0][j][k] = Coord[k]; // Définit le paramètre pour calculer p1
                }
            }
            else if (i==0 && j>0 && j<3){
                for (k=0;k<Nb;k++){
                    Mat_Pond[0][j][k] = Coord[k] + dt*Mat_Pond[1][j-1][k]/2;  // Définit le paramètre pour calculer p2 et p3 dépendant de p1 et p2
                }
            }
            else if (i==0 && j==3){
                for (k=0;k<Nb;k++){ // Définit le paramètre pour calculer p4 dépendant de p3
                    Mat_Pond[0][j][k] = Coord[k] + dt*Mat_Pond[1][j-1][k];
                }
            }
            else if (i==1 && j==0){
	            f(Mat_Pond[0][j],Mat_Pond[1][j],t); // Définit p1
            }
            else if (i==1 && j>0 && j<3){
                f(Mat_Pond[0][j],Mat_Pond[1][j],t+dt/2); // Définit p2 et p3
            }
            else if (i==1 && j==3){
	            f(Mat_Pond[0][j],Mat_Pond[1][j],t+dt); // Définit p4
            }
        }// k porte sur les différentes variables du système
    }
}

void rk_4(void (*f)(double*,double*,double),double* Coord,double tf,double dt,int Nb,const char* Nom_fichier)
{
    double t=0;
    double ***Mat_Pond;
    Mat_Pond = (double***)malloc(2*sizeof(double**)); // Définit la matrice 2*4*Nb qui portera les variables du
    int i,j,k; // système ainsi que les paramètres utiles de la méthode de Runge-Kutta
    for (i=0;i<2;i++){
        Mat_Pond[i] = (double**)malloc(4*sizeof(double*));
        for (j=0;j<4;j++){
            Mat_Pond[i][j] = (double*)malloc(Nb*sizeof(double));
        }
    }
    Ini_Mat(f,Mat_Pond,Coord,t,dt,Nb); // On initialise les pi pour t=0
    
    fstream fich;
    fich.open(Nom_fichier,ios::out);
    fich << t; // On affiche les conditions initiales dans un fichier texte pour l'afficher plus tard sous Python
    for (k=0;k<Nb_var;k++){
      fich <<' '<<Coord[k];
    }
    fich <<endl;
    while (t<=tf){ // Avancée temporelle
        t += dt;
	fich <<t;
        for (k=0;k<Nb_var;k++){ // Pour chaques coordonnées on avance par développement
            Coord[k] = Coord[k] + dt*(Mat_Pond[1][0][k] + 2*Mat_Pond[1][1][k] + 2*Mat_Pond[1][2][k] + Mat_Pond[1][3][k])/6;
	    fich << ' ' << Coord[k];
        }
        fich <<endl;
        Ini_Mat(f,Mat_Pond,Coord,t,dt,Nb); // On fait évoluer les pi à t+dt pour la prochaine récurrence du système
    }
    fich.close();
    for (i=0;i<2;i++){ // On libère la mémoire
        for (j=0;j<4;j++){
	  free(Mat_Pond[i][j]);
	}
	free(Mat_Pond[i]);
    }
    free(Mat_Pond);
}

// Mise en équation

void pendule_3D(double* Coord,double* dCoord,double t){ //Equation du pendule simple en 3D
  double d2Teta,d2Phi;
  double Teta = Coord[0],dTeta = Coord[1],/*Phi = Coord[2],*/dPhi = Coord[3];
  d2Teta = (Rho*dPhi*dPhi*sin(Teta)*cos(Teta)-g*sin(Teta))/Rho;
  d2Phi = -(2*dTeta*dPhi*cos(Teta))/(sin(Teta));
  dCoord[0] = Coord[1]; dCoord[1] = d2Teta; dCoord[2] = Coord[3]; dCoord[3] = d2Phi;
}

int main(){
    double* Coord;
    double tf=100,dt=0.01;
    const char* Nom_fichier = "res_1pend.txt";
    Coord = (double*)malloc(Nb_var*sizeof(double));
    Coord[0] = M_PI/3;Coord[1] = 0;Coord[2] = 0;Coord[3] = M_PI/10;
    rk_4(pendule_3D,Coord,tf,dt,Nb_var,Nom_fichier);
    return 0;
}
