#include "electron.h"

#include <fstream>

// Class constructor. All necessary calculations and tests are left to the main loop.
Electron::Electron(int ID, double s0, double vpar0, double vper0, double mue)
	: id{ ID }, position{ s0 }, parallelSpeed{ vpar0 }, perpendicularSpeed{ vper0 }, magneticMoment{ mue }
{};

int Electron::get_id() { return id; };

double Electron::get_position() {
	return position[position.size()-1];
};
double Electron::get_parallelSpeed() {
	return parallelSpeed[position.size()-1];
};
double Electron::get_perpendicularSpeed() {
	return perpendicularSpeed[position.size()-1];
};
double Electron::get_magneticMoment() {
	return magneticMoment;
};

void Electron::set_position(double s) {
	position.push_back(s);
};
void Electron::set_parallelSpeed(double vpar) {
	parallelSpeed.push_back(vpar);
};
void Electron::set_perpendicularSpeed(double vper) {
	perpendicularSpeed.push_back(vper);
};


std::ostream& operator<<(std::ostream& os, const Electron& e)
{
	os << "s(m), v_par(m/s), v_perp(m/s)\n";

	for (int i = 0; i < e.position.size(); i++)
		os << e.position[i] << ", " << e.parallelSpeed[i] << ", " << e.perpendicularSpeed[i] << '\n';

	return os;
}
