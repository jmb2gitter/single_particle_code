#include "domain.h"
#include <cmath>
#include <iostream>

Domain::Domain(double Lshell, int number_cells) : L{ Lshell } {

    // latitude at which the magnetic field line intersects the earths surface (rad)
    // Basic Space Plasma Physics - Baumjohann - eq 3.9
    const double li = std::acos(1.0 / std::sqrt(L));

    // si: the distance[m] between the equatorial plane and the surface of the earth along the magnetic field line L
    const double si = L * Re * ( 3.0 * std::sin(li) * std::sqrt( 1.0 + 3.0 * std::pow(std::sin(li), 2.0) ) +
        std::sqrt(3.0) * std::asinh( std::sqrt(3.0) * std::sin(li) ) ) / 6.0;

    // spatial discretization interval. In this case the spatial domain covers the range [-si:si)
    delta_s = 2.0 * si / double(number_cells);

    // the spatial domain coordinate 'space [m]' vector initialization
    for (int i = 0; i < number_cells; i++) space.push_back(double(i) * delta_s - si);


    /* magnetic field array - magnetic bottle B = Beq(1 + (s / so) ^ 2) [T]. (eq. 2.74, pág 41) */

    // magnetic field intensity at que equator:     Beq = Be/L^3
    double Beq = Be / std::pow(L, 3);
    //std::cout << "Beq = " << Beq << std::endl;

    // magnetic field intensity at the ionosphere:  Bi = Beq * sqrt(1 + 3*sin(li)^2)/cos(li)^6
    // Basic Space Plasma Physics - Baumjohann - eq 3.8
    double Bi = Beq * std::sqrt(1.0 + 3.0 * std::pow(std::sin(li), 2)) / std::pow(std::cos(li), 6);
    //std::cout << "Bi = " << Bi << std::endl;

    // define a spatial parameter for numerical convenience: so = si / sqrt(Bi/Beq - 1)
    double so = si / std::sqrt(Bi / Beq - 1.0);
    //std::cout << "so = " << so << std::endl;

    // magnetic intensity (and magnetic gradient) vector declaration & initialization:
    for (int i = 0; i < space.size(); i++) {
        magnetic_field.push_back(Beq * (1.0 + std::pow(space[i] / so, 2)));
        magnetic_gradient.push_back(2.0 * Beq * space[i] / std::pow(so, 2));
    }

    // outputting spatial profiles for validation purposes
    std::ofstream write_output("spatial_profiles.csv");
    if (write_output.is_open()) {
        write_output << "position (Re),magnetic field (T),magnetic gradient (T/m)" << '\n';
        for (int i = 0; i < space.size(); i++) write_output << space[i] << ',' << magnetic_field[i] << ',' << magnetic_gradient[i] << '\n';
    }
    write_output.close();

}

const std::vector<double> Domain::get_space() { return space; }
const std::vector<double> Domain::get_magnetic_field() { return magnetic_field; }
const std::vector<double> Domain::get_magnetic_gradient() { return magnetic_gradient; }
