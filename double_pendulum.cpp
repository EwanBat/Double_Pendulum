// Import de librairies et définition des constantes importantes
#include<iostream>
#include<fstream>
#include<math.h>
#define M_PI           3.14159265358979323846
#include<stdlib.h>
using namespace std;

const double g = 9.81; // Constante gravitationnelle terrestre normalisée en m.s-2
double m1 = 1,m2 = 1; // Masses respectives du premier et second pendule en kg
double l1 = 1,l2 = 1; // Longueurs respectives du premier et second pendule en m
int Nb_var = 8; // Nombre totale de variable pour la méthode de Runge-Kutta
double tf=100,dt=0.01; // Temps final et pas de temps en s

int Delta(int i,int j){ // Delta de Dirac, on fera attention à donner i+1,j car les listes commencent à 0
  if (i == j){
    return 1;
  }
  else if (i != j){
    return 0;
  }
  return 0;
}

// Création de la fonction rk4

void Ini_Mat(void (*f)(double*,double*,double),double*** Mat_Pond,double* Coord,double t,double dt,int Nb){
  int i,j,k;           // Fonction qui construit un tableau contenant les valeurs p1,p2,p3,p4 de la méthode rk4
  for (j=0;j<4;j++){ // Selon i, on fixe soit pj soit le terme dans le système qui permet de calculer pj
        for (i=0;i<2;i++){
            if (i==0 && j==0){
	      for (k=0;k<Nb;k++){ // Pour chaques coordonnées
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

void rk_4(void (*f)(double*,double*,double),double* Coord,double tf,double dt,int Nb,const char* Nom_fichier )
{ //Fonction qui a pour variables un système d'équation, les coordonnées initiales, le temps final, un pas de temps, le nombre d'inconnues et le fichier pour les valeurs
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
    fich << t; // On enregistre les conditions initiales dans un fichier texte pour l'afficher plus tard sous Python
    for (k=0;k<Nb_var;k++){
      fich <<' '<<Coord[k];
    }
    fich <<endl;

    while (t<=tf){ // Avancée temporelle
        t += dt;
	fich <<t;
        for (k=0;k<Nb_var;k++){ // Pour chaques coordonnées on avance par développement limité
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
/* On note q3 et q4 les coordonnées stéréographiques du pendule 1 et q1,q2 les coordonnées stéréographiques du pendule 2*/
/* On note p3 et p4 les impulsions stéréographiques du pendule 1 et q1,q2 les impulsions stéréographiques du pendule 2*/

////////////////  On utilise des fonctions tiers afin de "simplifier" la dérivation de l'hamiltonien ///////
double f_a(double* Coord){
  double q1 = Coord[0],q2 = Coord[1];
  return 1 + q1*q1 + q2*q2;
}

double f_b(double* Coord){
  double q3 = Coord[2],q4 = Coord[3];
  return 1 + q3*q3 + q4*q4;
}

double f_k(double* Coord){
  double q1 = Coord[0],q2 = Coord[1],q3 = Coord[2],q4 = Coord[3];
  return (q1-q3)*(q1-q3) + (q2-q4)*(q2-q4);
}

double f_A(double* Coord){
  double q2 = Coord[1],q4 = Coord[3];
  double p1 = Coord[4];
  double a = f_a(Coord), k = f_k(Coord);
  return p1*l2*a*(a*(q4-q2) + q2*k);
}

double f_B(double* Coord){
  double q1 = Coord[0],q3 = Coord[2];
  double p2 = Coord[5];
  double a = f_a(Coord), k = f_k(Coord);
  return p2*l2*a*(a*(q3-q1) + q1*k);
}

double f_C(double* Coord){
  double q2 = Coord[1],q4 = Coord[3];
  double p3 = Coord[6];
  double b = f_b(Coord), k = f_k(Coord);
  return p3*l1*b*(b*(q2-q4) + q4*k);
}

double f_D(double* Coord){
  double q1 = Coord[0],q3 = Coord[2];
  double p4 = Coord[7];
  double b = f_b(Coord), k = f_k(Coord);
  return p4*l1*b*(b*(q1-q3) + q3*k);
}
  
double f_y1(double* Coord){
  double a = f_a(Coord),b = f_b(Coord);
  double p1 = Coord[4];
  return p1*l2*a*a*b/2;
}

double f_y2(double* Coord){
  double a = f_a(Coord),b = f_b(Coord);
  double p2 = Coord[5];
  return p2*l2*a*a*b/2;
}

double f_y3(double* Coord){
  double a = f_a(Coord),b = f_b(Coord);
  double p3 = Coord[6];
  return p3*l1*b*b*a/2;
}

double f_y4(double* Coord){
  double a = f_a(Coord),b = f_b(Coord);
  double p4 = Coord[7];
  return p4*l1*b*b*a/2;
}
	      
double f_F(double* Coord){
  double q1 = Coord[0],q2 = Coord[1];
  double p1 = Coord[4],p2 = Coord[5];
  double a = f_a(Coord),b = f_b(Coord);
  return a*b*l2*(q1*p1 + q2*p2);
}

double f_G(double* Coord){
  double q3 = Coord[2],q4 = Coord[3];
  double p3 = Coord[6],p4 = Coord[7];
  double a = f_a(Coord),b = f_b(Coord);
  return a*b*l1*(q4*p4 + q3*p3);
}

////////////////////////////// On effectue leurs dérivées selon les qi et pi ////////////////////////

double f_da(double* Coord,int i){
  if (i==0 || i==1){ // Dérive selon q1 ou q2
    return 2*Coord[i];
  }
  return 0;
}

double f_db(double* Coord,int i){
  if (i==2 || i==3){ // Dérive selon q3 ou q4
    return 2*Coord[i];
  }
  return 0;
}

double f_dk(double* Coord,int i){
  double q1 = Coord[0],q2 = Coord[1],q3 = Coord[2],q4 = Coord[3];
  if (i == 0) return 2*(q1-q3); // Dérive selon q1
  if (i == 1) return 2*(q2-q4); // Dérive selon q2
  if (i == 2) return 2*(q3-q1); // Dérive selon q3
  if (i == 3) return 2*(q4-q2); // Dérive selon q4
  return 0;
}

double f_dA(double* Coord,int i){
  double q2 = Coord[1],q4 = Coord[3];
  double p1 = Coord[4];
  double a = f_a(Coord), k = f_k(Coord);
  double da = f_da(Coord,i),dk = f_dk(Coord,i);
  if (0<=i && i<=3){ // Dérive selon qi
    return p1*l2*(da*(a*(q4-q2)+q2*k) + a*(da*(q4-q2) + a*(Delta(i+1,4) - Delta(i+1,2)) + Delta(i+1,2)*k + q2*dk));
  }
  if (4<=i && i<=7){ // Dérive selon pi
    return Delta(i+1-4,1)*l2*a*(a*(q4-q2) + q2*k);
  }
  return 0;
}

double f_dB(double* Coord,int i){
  double q1 = Coord[0],q3 = Coord[2];
  double p2 = Coord[5];
  double a = f_a(Coord), k = f_k(Coord);
  double da = f_da(Coord,i),dk = f_dk(Coord,i);
  if (0<=i && i<=3){ // Dérive selon qi
    return p2*l2*(da*(a*(q3-q1)+q1*k) + a*(da*(q3-q1) + a*(Delta(i+1,3) - Delta(i+1,1)) + Delta(i+1,1)*k + q1*dk));
  }
  if (4<=i && i<=7){ // Dérive selon pi
    return Delta(i+1-4,2)*l2*a*(a*(q3-q1) + q1*k);
  }
  return 0;
}

double f_dC(double* Coord,int i){
  double q2 = Coord[1],q4 = Coord[3];
  double p3 = Coord[6];
  double b = f_b(Coord), k = f_k(Coord);
  double db = f_db(Coord,i),dk = f_dk(Coord,i);
  if (0<=i && i<=3){ // Dérive selon qi
    return p3*l1*(db*(b*(q2-q4)+q4*k) + b*(db*(q2-q4) + b*(Delta(i+1,2) - Delta(i+1,4)) + Delta(i+1,4)*k + q4*dk));
  }
  if (4<=i && i<=7){ // Dérive selon pi
    return Delta(i+1-4,3)*l1*b*(b*(q2-q4) + q4*k);
  }
  return 0;
}

double f_dD(double* Coord,int i){
  double q1 = Coord[0],q3 = Coord[2];
  double p4 = Coord[7];
  double b = f_b(Coord), k = f_k(Coord);
  double db = f_db(Coord,i),dk = f_dk(Coord,i);
  if (0<=i && i<=3){ // Dérive selon qi
    return p4*l1*(db*(b*(q1-q3)+q3*k) + b*(db*(q1-q3) + b*(Delta(i+1,1) - Delta(i+1,3)) + Delta(i+1,3)*k + q3*dk));
  }
  if (4<=i && i<=7){ // Dérive selon pi
    return Delta(i+1-4,4)*l1*b*(b*(q1-q3) + q3*k);
  }
  return 0;
}

double f_dy1(double* Coord,int i){
  double a = f_a(Coord),b = f_b(Coord);
  double da = f_da(Coord,i),db = f_db(Coord,i);
  double p1 = Coord[4];
  if (0<=i && i<=3){ // Dérive selon qi
    return p1*l2*(2*a*da*b + db*a*a)/2;
  }
  if (i == 4){ // Dérive selon p1
    return l2*a*a*b/2;
  }
  return 0;
}

double f_dy2(double* Coord,int i){
  double a = f_a(Coord),b = f_b(Coord);
  double da = f_da(Coord,i),db = f_db(Coord,i);
  double p2 = Coord[5];
  if (0<=i && i<=3){ // Dérive selon qi
    return p2*l2*(2*a*da*b + db*a*a)/2;
  }
  if (i == 5){ // Dérive selon p2
    return l2*a*a*b/2;
  }
  return 0;
}

double f_dy3(double* Coord,int i){
  double a = f_a(Coord),b = f_b(Coord);
  double da = f_da(Coord,i),db = f_db(Coord,i);
  double p3 = Coord[6];
  if (0<=i && i<=3){ // Dérive selon qi
    return p3*l1*(2*a*db*b + da*b*b)/2;
  }
  if (i == 6){ // Dérive selon p3
    return l1*b*b*a/2;
  }
  return 0;
}

double f_dy4(double* Coord,int i){
  double a = f_a(Coord),b = f_b(Coord);
  double da = f_da(Coord,i),db = f_db(Coord,i);
  double p4 = Coord[7];
  if (0<=i && i<=3){ // Dérive selon qi 
    return p4*l1*(2*a*db*b + da*b*b)/2;
  }
  if (i == 7){ // Dérive selon p4
    return l1*b*b*a/2;
  }
  return 0;
}

double f_dF(double* Coord,int i){
  double q1 = Coord[0],q2 = Coord[1];
  double p1 = Coord[4],p2 = Coord[5];
  double a = f_a(Coord),b = f_b(Coord);
  double da = f_da(Coord,i),db = f_db(Coord,i);
  if (0<=i && i<=3){ // Dérive selon qi
    return l2*((da*b+a*db)*(q1*p1+q2*p2) + a*b*(Delta(i+1,1)*p1 + Delta(i+1,2)*p2));
  }
  if (4<=i && i<=7){ // Dérive selon pi
    return a*b*l2*(q1*Delta(i+1-4,1) + q2*Delta(i+1-4,2));
  }
  return 0;
}

double f_dG(double* Coord,int i){
  double q3 = Coord[2],q4 = Coord[3];
  double p3 = Coord[6],p4 = Coord[7];
  double a = f_a(Coord),b = f_b(Coord);
  double da = f_da(Coord,i),db = f_db(Coord,i);
  if (0<=i && i<=3){ // Dérive selon qi
    return l1*((da*b+a*db)*(q4*p4+q3*p3) + a*b*(Delta(i+1,3)*p3 + Delta(i+1,4)*p4));
  }
  if (4<=i && i<=7){ // Dérive selon pi
    return a*b*l1*(q4*Delta(i+1-4,4) + q3*Delta(i+1-4,3));
  }
  return 0;
}

///////////////// On écrit encore d'autres fonctions pour simplifier en ayant H = h/f + l ////////

double f_h(double* Coord){
  double q1 = Coord[0],q2 = Coord[1],q3 = Coord[2],q4 = Coord[3];
  double A = f_A(Coord),B = f_B(Coord),C = f_C(Coord),D = f_D(Coord);
  double y1 = f_y1(Coord),y2 = f_y2(Coord),y3 = f_y3(Coord),y4 = f_y4(Coord);
  double F = f_F(Coord),G = f_G(Coord);
  return (A - B + C - D)*(A - B + C - D)*m1/m2 + (A-B)*(A-B) + m2/m1*(y1*y1 + y2*y2) + ((y1 - y3)*(y1 - y3) + (y2 - y4)*(y2 - y4)) + 2*F*(q1*y3 + q2*y4) + 2*G*(q3*y1 + q4*y2) - 2*F*G*(1 + q1*q3 + q2*q4);
}

double f_f(double* Coord){
  double a = f_a(Coord),b = f_b(Coord), k = f_k(Coord);
  return 2*l1*l1*l2*l2*(m2*a*a*b*b + 4*m1*k*(a*b - k));
}

double f_dh(double* Coord,int i){ // Différentielle de h
  double q1 = Coord[0],q2 = Coord[1],q3 = Coord[2],q4 = Coord[3];
  double A = f_A(Coord),B = f_B(Coord),C = f_C(Coord),D = f_D(Coord);
  double dA = f_dA(Coord,i),dB = f_dB(Coord,i),dC = f_dC(Coord,i),dD = f_dD(Coord,i);
  double y1 = f_y1(Coord),y2 = f_y2(Coord),y3 = f_y3(Coord),y4 = f_y4(Coord);
  double dy1 = f_dy1(Coord,i),dy2 = f_dy2(Coord,i),dy3 = f_dy3(Coord,i),dy4 = f_dy4(Coord,i);
  double F = f_F(Coord),G = f_G(Coord),dF = f_dF(Coord,i),dG = f_dG(Coord,i);
  if (0<=i && i<=3){ // Dérive selon qi
    return 2*(dA-dB+dC-dD)*(A-B+C-D)*m1/m2 + 2*dA*(A-B) - 2*dB*(A-B)  + 2*m2/m1*(dy1*y1+dy2*y2) + 2*(dy1-dy3)*(y1-y3) + 2*(dy2-dy4)*(y2-y4) + 2*dF*(q1*y3+q2*y4) + 2*F*(Delta(i+1,1)*y3 + q1*dy3 + q2*dy4 + y4*Delta(i+1,2)) + 2*G*(Delta(i+1,3)*y1 + q3*dy1 + q4*dy2 + y2*Delta(i+1,4)) + 2*dG*(q3*y1 + q4*y2) - 2*F*G*(Delta(i+1+2,3)*q3 + Delta(i+1+2,4)*q4 + Delta(i+1-2,1)*q1 + Delta(i+1-2,2)*q2) - 2*(dF*G + dG*F)*(1+q1*q3+q2*q4);
  }
  if (4<=i && i<=7){ // Dérive selon pi
    return 2*(dA-dB+dC-dD)*(A-B+C-D)*m1/m2 + 2*dA*(A-B) - 2*dB*(A-B) + 2*m2/m1*(dy1*y1 + dy2*y2) + 2*(dy1-dy3)*(y1-y3) + 2*(dy2-dy4)*(y2-y4) + 2*dF*(q1*y3+q2*y4) + 2*F*(q1*dy3 + q2*dy4) + 2*G*(q3*dy1 + q4*dy2) + 2*dG*(q3*y1+q4*y2) - 2*(dF*G + F*dG)*(1+q1*q3+q2*q4);
  }
  return 0;
}

double f_df(double* Coord,int i){ // Différentielle de f
  double a = f_a(Coord),b = f_b(Coord), k = f_k(Coord);
  double da = f_da(Coord,i),db = f_db(Coord,i),dk = f_dk(Coord,i);
  return 2*l1*l1*l2*l2*(m2*(2*a*da*b*b + 2*b*db*a*a) + 4*m1*dk*(a*b-k) + 4*m1*k*(da*b+a*db-dk));
}

double f_dl(double* Coord,int i){ // Différentielle de l'énergie potentielle
  double a = f_a(Coord),b = f_b(Coord);
  double da = f_da(Coord,i),db = f_db(Coord,i);
  return g*((m1+m2)*(2*db/b/b)*l2 + m1*(2*da/a/a)*l1);
}

///////////////// On écrit le système et on résout ///////////////////////////////

void double_pendule_3D(double* Coord,double* dCoord,double t){
  int i;
  double h = f_h(Coord),f = f_f(Coord);
  double dh,df,dl;
  for (i=0;i<=7;i+=1){ // On applique les équations d'Hamilton pour q et p
    if (0<=i && i<=3){ // On dérive selon les pi pour calculer les qi point
      dh = f_dh(Coord,i+4);
      df = f_df(Coord,i+4);
      dl = f_dl(Coord,i+4);
      dCoord[i] = (dh*f - h*df)/(f*f) + dl;
    }
    if (4<= i && i<=7){ // On dérive selon les qi pour calculer les pi point
      dh = f_dh(Coord,i-4);
      df = f_df(Coord,i-4);
      dl = f_dl(Coord,i-4);
      dCoord[i] = -(dh*f - h*df)/(f*f) - dl;
    }
  }
}
  
int main(){
    double* Coord;
    double q1 = 0.7,q2 = 0,q3 = 0.9,q4 = 0.89;
    double p1 = 0,p2 = 0,p3 = 0,p4 = 0;
    const char* Nom_fichier = "res_2pend.txt"; //Fichier de stockage des résultats
    
    Coord = (double*)malloc(Nb_var*sizeof(double));
    Coord[0] = q1,Coord[1] = q2,Coord[2] = q3, Coord[3] = q4; // Coordonnées initiales
    Coord[4] = p1,Coord[5] = p2,Coord[6] = p3,Coord[7] = p4;
    rk_4(double_pendule_3D,Coord,tf,dt,Nb_var,Nom_fichier);

    fstream f_stereo,f_cart;
    f_stereo.open(Nom_fichier,ios :: in); // On récupère les coordonnées stéréographiques
    f_cart.open("cart_2pend.txt",ios :: out); // On va sortir les coordonnées cartésiennes
    double t,q1_t,q2_t,q3_t,q4_t,p1_t,p2_t,p3_t,p4_t;
    double x1,y1,z1,x2,y2,z2;
    while (f_stereo >> t >> q1_t >> q2_t >> q3_t >> q4_t >> p1_t >> p2_t >> p3_t >> p4_t) { // On passe des coordonnées stéréographiques en cartésiennes
      x2 = l2*2*q1_t/(1+q1_t*q1_t+q2_t*q2_t);
      y2 = l2*2*q2_t/(1+q1_t*q1_t+q2_t*q2_t);
      z2 = l2*(-1 + q1_t*q1_t + q2_t*q2_t)/(1 + q1_t*q1_t + q2_t*q2_t); // Formules données en annexe

      x1 = l1*2*q3_t/(1+q3_t*q3_t+q4_t*q4_t);
      y1 = l1*2*q4_t/(1+q3_t*q3_t+q4_t*q4_t);
      z1 = l1*(-1 + q3_t*q3_t + q4_t*q4_t)/(1 + q3_t*q3_t + q4_t*q4_t)+ l1 + l2 ;
      f_cart << x1 << " " << y1 << " " << z1 << " " << x1+x2 << " " << y1+y2 << " " << z1+z2 << " " <<endl;
      }
    f_stereo.close();
    f_cart.close();
    return 0;
}
