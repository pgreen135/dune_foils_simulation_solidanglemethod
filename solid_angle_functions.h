#ifndef SOLID_ANGLE_H
#define SOLID_ANGLE_H

// solid angle functions
// calculates number of pmt_hits for VUV light using solid angle method

#include <vector>
#include <iostream>
#include <cmath>

#include "TVector3.h"
#include "TRandom3.h"

namespace solid_angle{

	std::vector<double> SolidAngleAnalyzer(int Nphotons_created, const double quantum_efficiency, const double catcov, const double vuvfrac, const double visfrac, int pmt_number, TVector3 ScintPoint, TVector3 OpDetPoint);

	// solid angle of rectangle calculation
	// data structure
	struct acc{
		// ax,ay,az = centre of rectangle; w = width; h = height
  		double ax, ay, az, w, h; 
	};
 	// functions
	double omega(double a, double b, double d);
	double solid(acc& out, TVector3 v);

}

#endif