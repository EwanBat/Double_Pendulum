#pragma once

// Import de librairies
#include<iostream>
#include<fstream>
#include<math.h>
#define M_PI           3.14159265358979323846
#include<stdlib.h>
using namespace std;

class SimplePendulum {
public:
    SimplePendulum(double length, double mass, double time_step, double final_time)
        {
            this->m_l = length;
            this->m_mass = mass;
            this->m_dt = time_step;
            this->m_tf = final_time;
            this->m_Coord = new double[m_Nb_var];
        };
    
    ~SimplePendulum() {
        delete[] m_Coord;
    };
    
    void initial_conditions(double theta0, double dtheta0, double phi0, double dphi0)
        {
            this->m_theta0 = theta0;
            this->m_dtheta0 = dtheta0;
            this->m_phi0 = phi0;
            this->m_dphi0 = dphi0;
            // Initialiser Coord avec les conditions initiales
            m_Coord[0] = theta0;
            m_Coord[1] = dtheta0;
            m_Coord[2] = phi0;
            m_Coord[3] = dphi0;
        };
    
    void rk_4(const char* Nom_fichier);

private:
    double m_l;      // Length of the pendulum
    double m_mass;      // Mass of the pendulum bob
    double m_dt;     // Time step for the simulation
    double m_tf;     // Final time for the simulation
    double m_theta0;  // Angle from the vertical (in radians)
    double m_dtheta0;  // Angular velocity (in radians per second)
    double m_phi0;    // Azimuthal angle (in radians)
    double m_dphi0;   // Azimuthal angular velocity (in radians per second)

    double* m_Coord; // Current coordinates

    const double m_g = 9.81; // Constantante gravitationnelle terrestre en m.s-2
    const int m_Nb_var = 4; // Nombre de variable pour la méthode de Runge-Kutta

    void Ini_Mat(double*** Mat_Pond, double t);
    
    void pendule_3D(double* Coord,double* dCoord,double t); //Equation du pendule simple en 3D
};