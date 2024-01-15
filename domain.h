#ifndef DOMAIN_H
#define DOMAIN_H

#include <vector>
#include <fstream>

class Domain
{
	std::vector<double> space;					// spatial grid vector
	std::vector<double> magnetic_field;			// magnetic field vector
	std::vector<double> magnetic_gradient;		// gradient of the magnetic field

	/* ----- physical constants and parameters ----- */
	const double pi = 3.14159265358979323846;
	const double to_rad = pi / 180.0;

	const double qe = -1.602176462e-19;     // electron charge (C)
	const double me = 9.10938188e-31;       // electron mass (kg)
	const double mp = 1.67262192e-27;       // proton mass (kg)
	const double mu0 = 4.0e-7 * pi;         // vacuum magnetic permeability
	const double Be = 3.11e-5;              // equatorial magnetic field (T)
	const double Re = 6.371e6;              // earth radius (m)
	const double vc = 299792458;            // light speed (m/s)

	const double L;
	double delta_s;

public:

	// Class constructor.
	Domain(double Lshell, int number_cells);

	// inline access functions to some physical quantites
	double get_qe() const { return qe; };
	double get_me() const { return me; };

	// access functions to spatial profiles
	const std::vector<double> get_space();
	const std::vector<double> get_magnetic_field();
	const std::vector<double> get_magnetic_gradient();

};

#endif