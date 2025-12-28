#include "../include/simple_pendulum.hpp"

// Création de la fonction rk4

void SimplePendulum::Ini_Mat(double*** Mat_Pond, double t){
  int i,j,k;           // Fonction qui construit un tableau contenant les valeurs p1,p2,p3,p4 de la méthode rk4
  for (j=0;j<4;j++){ // Selon i, on fixe soit pj soit le terme dans le système qui permet de calculer pj
        for (i=0;i<2;i++){
            if (i==0 && j==0){
                for (k=0;k<Nb_var;k++){
		  Mat_Pond[0][j][k] = Coord[k]; // Définit le paramètre pour calculer p1
                }
            }
            else if (i==0 && j>0 && j<3){
                for (k=0;k<Nb_var;k++){
                    Mat_Pond[0][j][k] = Coord[k] + dt*Mat_Pond[1][j-1][k]/2;  // Définit le paramètre pour calculer p2 et p3 dépendant de p1 et p2
                }
            }
            else if (i==0 && j==3){
                for (k=0;k<Nb_var;k++){ // Définit le paramètre pour calculer p4 dépendant de p3
                    Mat_Pond[0][j][k] = Coord[k] + dt*Mat_Pond[1][j-1][k];
                }
            }
            else if (i==1 && j==0){
	            pendule_3D(Mat_Pond[0][j],Mat_Pond[1][j],t); // Définit p1
            }
            else if (i==1 && j>0 && j<3){
                pendule_3D(Mat_Pond[0][j],Mat_Pond[1][j],t+dt/2); // Définit p2 et p3
            }
            else if (i==1 && j==3){
	            pendule_3D(Mat_Pond[0][j],Mat_Pond[1][j],t+dt); // Définit p4
            }
        }// k porte sur les différentes variables du système
    }
}

void SimplePendulum::rk_4(const char* Nom_fichier)
{
    double t=0;
    double ***Mat_Pond;
    Mat_Pond = (double***)malloc(2*sizeof(double**)); // Définit la matrice 2*4*Nb qui portera les variables du
    int i,j,k; // système ainsi que les paramètres utiles de la méthode de Runge-Kutta
    for (i=0;i<2;i++){
        Mat_Pond[i] = (double**)malloc(4*sizeof(double*));
        for (j=0;j<4;j++){
            Mat_Pond[i][j] = (double*)malloc(Nb_var*sizeof(double));
        }
    }
    Ini_Mat(Mat_Pond, t); // On initialise les pi pour t=0
    
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
        Ini_Mat(Mat_Pond, t); // On fait évoluer les pi à t+dt pour la prochaine récurrence du système
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

void SimplePendulum::pendule_3D(double* Coord,double* dCoord,double t){ //Equation du pendule simple en 3D
  double d2Teta,d2Phi;
  double Teta = Coord[0],dTeta = Coord[1],/*Phi = Coord[2],*/dPhi = Coord[3];
  d2Teta = (l*dPhi*dPhi*sin(Teta)*cos(Teta)-g*sin(Teta))/l;
  d2Phi = -(2*dTeta*dPhi*cos(Teta))/(sin(Teta));
  dCoord[0] = Coord[1]; dCoord[1] = d2Teta; dCoord[2] = Coord[3]; dCoord[3] = d2Phi;
}
