/*
single_particle_code.cpp:

This version simulates the trajectory of a single electron inside a magnetic bottle.
*/

#include <iostream>
#include <numbers>
#include <vector>
#include <string>
#include <fstream>
#include <cmath>

#include "domain.h"
#include "electron.h"


/* function declarations */

double index_from_position(double x, double dx, double x0);
double interpolate(std::vector<double>& X, std::vector<double>& Y, double x);

// Runge Kutta 4th order solver for 2nd order ODE
bool RungeKutta4(double dt, std::vector<double>& s, std::vector<double>& B, std::vector<double>& gradB, Electron& e);
double f1(double x1, double x2);
double f2(double x1, double x2, Electron& e, std::vector<double>& s, std::vector<double>& gradB);

int main()
{
    // some annoying math constants the GUI won't figure out -_-
    const double pi = 3.14159265358979323846;
    const double to_rad = pi / 180.0;

    /* ----- simulation_environment_setup ----- */

    // magnetic field line. Value must be > 0.
    const double L = 4.0;

    // number of spatial cells on each hemisphere. Value must be > 0.
    const int ns = 500;

    // spatial domain object
    Domain sim(L, ns);


    // reference vector for the spatial domain constructed in our Domain object
    std::vector<double> s = sim.get_space();
    std::vector<double> magField = sim.get_magnetic_field();
    std::vector<double> magGrad = sim.get_magnetic_gradient();

    /* --------------------- put some electrons at the origin and follow its orbits ------------------------- */
    std::vector<Electron> particles;

    // keep track of the number of created particles
    int e_count = 0;

    // all electrons are initially located at the equator of the magnetic bottle.
    int eqtr{ int( (0 - s[0])/(s[1] - s[0]) ) };
    double speed, energy, v_par, v_per, mu_e;
    double me = sim.get_me();
    double qe = sim.get_qe();

    // loop for energies
    for (int i = 0; i <= 3; i++)

        // loop for pitch angles
        for (int p = 15; p <= 75; p += 15) {

            e_count += 1;

            energy = std::abs(qe) * std::pow(10.0, double(i));
            speed = std::sqrt(2.0 * energy / me);
            v_par = speed * std::cos(to_rad * p);
            v_per = speed * std::sin(to_rad * p);
            mu_e = energy * std::pow(std::sin(to_rad * p), 2.0) / magField[eqtr];

            Electron temp(e_count, s[eqtr], v_par, v_per, mu_e);
            particles.push_back(temp);

        }

    std::cout << "done creating " << e_count << " electrons\n";

    // advance the simulation up to 20.0 s in intervals 'dt'.
    double dt{ .01 };    // time discretization unit
    double tc{ .0 };     // current time into simulation

    // the simulation runs until 20.0 s, or until all particles abandon the spatial grid
    while (tc <= 20.0 && !particles.empty()) {

        std::cout << "tc = " << tc << " s\n";

        for (auto itr = particles.begin(); itr != particles.end(); ) {

            //std::cout << "iterating particle " << (*itr).get_id() << '\n';

            if (!RungeKutta4(dt, s, magField, magGrad, *itr)) {

                // output data to file before calling destructor on electron
                std::ofstream write_output("electron" + std::to_string((*itr).get_id()) + ".csv");

                if (write_output.is_open()) write_output << *itr << std::endl;

                write_output.close();

                itr = particles.erase(itr);
            }
            else ++itr;
        }

        tc += dt;
    }

    // output the trajectories of the remaining particles
    for (auto& electron : particles ) {

        // output data to file before calling destructor on electron
        std::ofstream write_output("electron" + std::to_string(electron.get_id()) + ".csv");

        if (write_output.is_open()) write_output << electron << std::endl;

        write_output.close();
    }

}


/*  function index from position: double index(x)
    determines the cell number on the positions array, where a given coordinate 'x' is located
    No tests implemented to ensure valid outputs yet. Value must be higher than zero. */
double index_from_position(double x, double dx, double x0) { return ((x - x0) / dx); }

/*  function interpolate:
    this function takes an array Y, and its corresponding spatial domain X, and interpolates
    the value 'y' that corresponds to its coordinate 'x'. */
double interpolate(std::vector<double> &X, std::vector<double>& Y, double x) {

    // find the cell number, nx, where 'x' is located in the array X
    double nx{ (x - X[0]) / (X[1] - X[0]) };

    // locate the indexes closer to that position
    int n1{ int(std::floor(nx)) };
    int n2{ n1 + 1 };

    // interpolate
    return Y[n2] * (nx - n1) + Y[n1] * (n2 - nx);
}


// Ejecutar programa: Ctrl + F5 o menú Depurar > Iniciar sin depurar
// Depurar programa: F5 o menú Depurar > Iniciar depuración

// Sugerencias para primeros pasos: 1. Use la ventana del Explorador de soluciones para agregar y administrar archivos
//   2. Use la ventana de Team Explorer para conectar con el control de código fuente
//   3. Use la ventana de salida para ver la salida de compilación y otros mensajes
//   4. Use la ventana Lista de errores para ver los errores
//   5. Vaya a Proyecto > Agregar nuevo elemento para crear nuevos archivos de código, o a Proyecto > Agregar elemento existente para agregar archivos de código existentes al proyecto
//   6. En el futuro, para volver a abrir este proyecto, vaya a Archivo > Abrir > Proyecto y seleccione el archivo .sln


/* Runge Kutta 4th order solver for 2nd order ODE

http://lampx.tugraz.at/~hadley/num/ch8/rk4ode2.php

Any second order differential equation can be written as two coupled first order equations,
dx/dt = f1(x, v, t)
dv/dt = f2(x1, x2, t).

The relevant equations for an electron moving inside a magnetic bottle take the form:
ds/dt = v					--->	f1(x, v, t) = v
dv/dt = - mue * gradB / me	--->	f2(x, v, t) = - mue * gradB(x) / me

with intial conditions x = xo, v = vo.
*/
double f1(double x, double v) { return v; };

double f2(double x, double v, Electron& e, std::vector<double>& s, std::vector<double>& gradB) {
    return -e.get_magneticMoment() / e.Me * interpolate(s, gradB, x);
};


// this function advances the motion in time of the electron 'e'. It returns TRUE if after iteration the particle
// remains within the simulation domain, and returns FALSE otherwise.
bool RungeKutta4(double dt, std::vector<double>& s, std::vector<double>& B, std::vector<double>& gradB, Electron& e) {

    double x1 = e.get_position();
    double x2 = e.get_parallelSpeed();
    double x3 = e.get_perpendicularSpeed();

    // estimate the magnetic intensity at the beginning of the iteration
    double B_initial = interpolate(s, B, x1);

    double ks1 = dt * f1(x1, x2);
    double kv1 = dt * f2(x1, x2, e, s, gradB);
    if (ks1 < s[0] || ks1 > s[s.size() - 1]) return false;

    double ks2 = dt * f1(x1 + 0.5 * ks1, x2 + 0.5 * kv1);
    double kv2 = dt * f2(x1 + 0.5 * ks1, x2 + 0.5 * kv1, e, s, gradB);
    if (ks2 < s[0] || ks2 > s[s.size() - 1]) return false;

    double ks3 = dt * f1(x1 + 0.5 * ks2, x2 + 0.5 * kv2);
    double kv3 = dt * f2(x1 + 0.5 * ks2, x2 + 0.5 * kv2, e, s, gradB);
    if (ks3 < s[0] || ks3 > s[s.size() - 1]) return false;

    double ks4 = dt * f1(x1 + ks3, x2 + kv3);
    double kv4 = dt * f2(x1 + ks3, x2 + kv3, e, s, gradB);
    if (ks4 < s[0] || ks4 > s[s.size() - 1]) return false;

    x1 += (ks1 + 2.0 * (ks2 + ks3) + ks4) / 6.0;
    x2 += (kv1 + 2.0 * (kv2 + kv3) + kv4) / 6.0;

    // update the electrons trajectory
    double B_final = interpolate(s, B, x1);
    x3 *= std::sqrt(std::abs(B_final / B_initial));

    e.set_position(x1);
    e.set_parallelSpeed(x2);
    e.set_perpendicularSpeed(x3);

    return true;
}
