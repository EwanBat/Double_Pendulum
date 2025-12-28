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
            this->l = length;
            this->m = mass;
            this->dt = time_step;
            this->tf = final_time;
        };
    
    void initial_conditions(double theta0, double dtheta0, double phi0, double dphi0)
        {
            this->theta0 = theta0;
            this->dtheta0 = dtheta0;
            this->phi0 = phi0;
            this->dphi0 = dphi0;
        };
    
    void rk_4(const char* Nom_fichier);

private:
    double l;      // Length of the pendulum
    double m;      // Mass of the pendulum bob
    double dt;     // Time step for the simulation
    double tf;     // Final time for the simulation
    double theta0;  // Angle from the vertical (in radians)
    double dtheta0;  // Angular velocity (in radians per second)
    double phi0;    // Azimuthal angle (in radians)
    double dphi0;   // Azimuthal angular velocity (in radians per second)

    double* Coord; // Current coordinates

    const double g = 9.81; // Constantante gravitationnelle terrestre en m.s-2
    const int Nb_var = 4; // Nombre de variable pour la méthode de Runge-Kutta

    void Ini_Mat(double*** Mat_Pond, double t);
    
    void pendule_3D(double* Coord,double* dCoord,double t); //Equation du pendule simple en 3D
};