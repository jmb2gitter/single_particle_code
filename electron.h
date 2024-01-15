#ifndef ELECTRON_H
#define ELECTRON_H

#include <vector>
#include <fstream>

class Electron
{
	int id;		// each electron is assigned an ID, used to name the trajectory file where data is to be stored.

	std::vector<double> position;
	std::vector<double> parallelSpeed;
	std::vector<double> perpendicularSpeed;
	double magneticMoment;

public:
	static constexpr double Qe{ -1.602176462e-19 };			// electron charge (C)
	static constexpr double Me{ 9.10938188e-31 };			// electron mass (kg)

	// Class constuctor
	Electron(int ID, double s0, double vpar0, double vper0, double mue);

	int get_id();
	double get_position();
	double get_parallelSpeed();
	double get_perpendicularSpeed();
	double get_magneticMoment();

	void set_position(double s);
	void set_parallelSpeed(double vpar);
	void set_perpendicularSpeed(double mue);

	friend std::ostream& operator<<(std::ostream& os, const Electron& e);

};

#endif