#include "../include/simple_pendulum.hpp"

int main(){
    // Définition des paramètres du pendule simple
    double length = 1.0;      // Longueur du pendule en mètres
    double mass = 1.0;        // Masse du pendule en kg
    double time_step = 0.01;  // Pas de temps en secondes
    double final_time = 100.0; // Temps final en secondes

    // Création de l'objet SimplePendulum
    SimplePendulum pendulum(length, mass, time_step, final_time);

    // Définition des conditions initiales
    double theta0 = M_PI / 4;   // Angle initial (45 degrés)
    double dtheta0 = 0.0;       // Vitesse angulaire initiale
    double phi0 = 0.0;          // Angle azimutal initial
    double dphi0 = 0.0;         // Vitesse angulaire azimutale initiale

    pendulum.initial_conditions(theta0, dtheta0, phi0, dphi0);

    // Exécution de la simulation et enregistrement des résultats
    const char* output_file = "simple_pendulum_output.txt";
    pendulum.rk_4(output_file);

    return 0;
}