#pragma once

// Import de librairies
#include<iostream>
#include<fstream>
#include<math.h>
#define M_PI           3.14159265358979323846
#include<stdlib.h>
using namespace std;

class DoublePendulum {
public:
    DoublePendulum(double length1, double length2, double mass1, double mass2, double time_step, double final_time)
        {
            this->m_l1 = length1;
            this->m_l2 = length2;
            this->m_m1 = mass1;
            this->m_m2 = mass2;
            this->m_dt = time_step;
            this->m_tf = final_time;
            this->m_Coord = new double[m_Nb_var];
        };
    
    ~DoublePendulum() {
        delete[] m_Coord;
    };
    
    void initial_conditions(double q1, double q2, double q3, double q4, 
                           double p1, double p2, double p3, double p4)
        {
            this->m_q1 = q1;
            this->m_q2 = q2;
            this->m_q3 = q3;
            this->m_q4 = q4;
            this->m_p1 = p1;
            this->m_p2 = p2;
            this->m_p3 = p3;
            this->m_p4 = p4;
            // Initialiser Coord avec les conditions initiales
            m_Coord[0] = q1;
            m_Coord[1] = q2;
            m_Coord[2] = q3;
            m_Coord[3] = q4;
            m_Coord[4] = p1;
            m_Coord[5] = p2;
            m_Coord[6] = p3;
            m_Coord[7] = p4;
        };
    
    void rk_4(const char* Nom_fichier);

private:
    double m_l1;      // Longueur du premier pendule en m
    double m_l2;      // Longueur du second pendule en m
    double m_m1;      // Masse du premier pendule en kg
    double m_m2;      // Masse du second pendule en kg
    double m_dt;      // Pas de temps en s
    double m_tf;      // Temps final en s
    
    // Conditions initiales
    double m_q1, m_q2, m_q3, m_q4;  // Coordonnées stéréographiques
    double m_p1, m_p2, m_p3, m_p4;  // Impulsions stéréographiques

    double* m_Coord; // Current coordinates

    const double m_g = 9.81; // Constante gravitationnelle terrestre en m.s-2
    const int m_Nb_var = 8; // Nombre de variable pour la méthode de Runge-Kutta

    // Méthodes privées
    int Delta(int i, int j);
    void Ini_Mat(double*** Mat_Pond, double t);
    void double_pendule_3D(double* Coord, double* dCoord, double t);
    
    // Fonctions auxiliaires
    double f_a(double* Coord);
    double f_b(double* Coord);
    double f_k(double* Coord);
    double f_A(double* Coord);
    double f_B(double* Coord);
    double f_C(double* Coord);
    double f_D(double* Coord);
    double f_y1(double* Coord);
    double f_y2(double* Coord);
    double f_y3(double* Coord);
    double f_y4(double* Coord);
    double f_F(double* Coord);
    double f_G(double* Coord);
    
    // Dérivées
    double f_da(double* Coord, int i);
    double f_db(double* Coord, int i);
    double f_dk(double* Coord, int i);
    double f_dA(double* Coord, int i);
    double f_dB(double* Coord, int i);
    double f_dC(double* Coord, int i);
    double f_dD(double* Coord, int i);
    double f_dy1(double* Coord, int i);
    double f_dy2(double* Coord, int i);
    double f_dy3(double* Coord, int i);
    double f_dy4(double* Coord, int i);
    double f_dF(double* Coord, int i);
    double f_dG(double* Coord, int i);
    
    // Hamiltonien
    double f_h(double* Coord);
    double f_f(double* Coord);
    double f_dh(double* Coord, int i);
    double f_df(double* Coord, int i);
    double f_dl(double* Coord, int i);
};
