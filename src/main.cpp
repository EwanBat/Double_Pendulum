#include "../include/simple_pendulum.hpp"
#include "../include/double_pendulum.hpp"

int main(){
    //////////////////////////////// Simple Pendulum Simulation ////////////////////////////////
    // Définition des paramètres du pendule simple
    double length = 1.0;      // Longueur du pendule en mètres
    double mass = 1.0;        // Masse du pendule en kg
    double time_step = 0.001;  // Pas de temps en secondes
    double final_time = 10; // Temps final en secondes

    // Création de l'objet SimplePendulum
    SimplePendulum pendulum(length, mass, time_step, final_time);

    // Définition des conditions initiales
    double theta0 = M_PI/4;   // Angle initial (45 degrés) (Ne doit pas être égal à 0 pour observer le mouvement)
    double dtheta0 = 0.0;       // Vitesse angulaire initiale
    double phi0 = 0.0;          // Angle azimutal initial
    double dphi0 = -0.1;         // Vitesse angulaire azimutale initiale

    pendulum.initial_conditions(theta0, dtheta0, phi0, dphi0);

    // Exécution de la simulation et enregistrement des résultats
    const char* output_file = "../data/simple_pendulum_output.txt";
    pendulum.rk_4(output_file);

    // Run plotting script
    std::string cmd = "python3 ../src/simple_pendulum_plot.py";
    int result = std::system(cmd.c_str());
    if (result != 0) {
        std::cerr << "Error executing plotting script." << std::endl;
    }

    //////////////////////////////// Double Pendulum Simulation ////////////////////////////////
    // Définition des paramètres du double pendule
    double length1 = 1.0;      // Longueur du premier pendule en mètres
    double length2 = 1.0;      // Longueur du second pendule en mètres
    double mass1 = 1.0;        // Masse du premier pendule en kg
    double mass2 = 1.0;        // Masse du second pendule en kg
    time_step = 0.001;  // Pas de temps en secondes
    final_time = 100; // Temps final en secondes

    // Création de l'objet DoublePendulum
    DoublePendulum double_pendulum(length1, length2, mass1, mass2, time_step, final_time);

    // Définition des conditions initiales
    double q1 = 0.0;   // Coordonnée stéréographique du premier pendule
    double q2 = -0.1;   // Coordonnée stéréographique du premier pendule
    double q3 = 0.1;   // Coordonnée stéréographique du second pendule
    double q4 = 0.0;   // Coordonnée stéréographique du second pendule
    double p1 = 0.1;   // Impulsion stéréographique du premier pendule
    double p2 = 0.0;   // Impulsion stéréographique du premier pendule
    double p3 = 0.0;   // Impulsion stéréographique du second pendule
    double p4 = 0.0;   // Impulsion stéréographique du second pendule

    double_pendulum.initial_conditions(q1, q2, q3, q4, p1, p2, p3, p4);

    // Exécution de la simulation et enregistrement des résultats
    const char* double_output_file = "../data/double_pendulum_output_bis.txt";
    double_pendulum.rk_4(double_output_file);
    // Run plotting script
    std::string double_cmd = "python3 ../src/double_pendulum_plot.py";
    result = std::system(double_cmd.c_str());
    if (result != 0) {
        std::cerr << "Error executing double pendulum plotting script." << std::endl;
    }

    return 0;
}