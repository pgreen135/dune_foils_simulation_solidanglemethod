#ifndef SOLID_ANGLE_H
#define SOLID_ANGLE_H

// class containing functions for calculating number of hits on each optical channel via the solid angle method
// calculates number of VUV hits, using gaisser-hillas corrections for Rayleigh scattering
// calculated number of visible hits [preliminary]

#include <vector>
#include <string>

#include "TF1.h"
#include "TVector3.h"

class solid_angle {

private:

	// useful constants
	const double pi = 3.141592653589793;

	// detector type flag
	const std::string flagDetector;


	// *************************************************************************************************
	//                    	OPTICAL DETECTOR SHAPE/SIZE
	// *************************************************************************************************	
	const int optical_detector_type;
	// Arapucas: detector_type = 0;	[matching size used in Diego's .gdml for "small_bars"]
	double height = 10.16;	// cm
	double width = 10.16;	// cm
	// PMTs: decector_type = 1; 
	double raidus = 8*2.54/2.;	//	8" PMT diameter  // cm	// not in use

	// structure definition for solid angle of rectangle function
	struct acc{
		// ax,ay,az = centre of rectangle; w = width; h = height
  		double ax, ay, az, w, h; 
	};


	// *************************************************************************************************
	//                    NUMBER OF VUV HITS PARAMETRIZATION
	// *************************************************************************************************
	// VUV hits Gaisser-Hillas Rayleigh scattering correction
	TF1* GH[9];
	//  SBN Gaisser-Hillas
	const int flagRS;
	const double GH_RS60cm_SBN[4][9] = { {1.31326, 1.28293, 1.24304, 1.13739, 1.03286, 0.904908, 0.779762, 0.654461, 0.525548},
					     {87.5397, 95.0615, 89.6917, 94.9943, 93.9111, 111.708, 114.998, 132.535, 144.225},
					     {57.3686, 59.9412, 53.8011, 56.1887, 63.5782, 61.0104, 66.3173, 63.379, 64.4428},
					     {-200, -200, -200, -200, -200, -200, -200, -200, -200} };
	const double GH_RS120cm_SBN[4][9] = { {1.16866, 1.13776, 1.08677, 0.993735, 0.885896, 0.760942, 0.628574, 0.498151, 0.37109}, 
					      {86.5078, 96.3383, 90.7074, 97.8305, 99.7487, 119.343, 129.554, 148.349, 174.775}, 
					      {88.0653, 89.4535, 82.9928, 84.9811, 93.5736, 89.2581, 93.4002, 91.709, 86.177}, 
					      {-200, -200, -200, -200, -200, -200, -200, -200, -200} };
	const double GH_RS180cm_SBN[4][9] = { {1.11436, 1.08245, 1.02718, 0.939834, 0.83028, 0.704374, 0.566687, 0.430942, 0.299605}, 
					      {85.1942, 95.5351, 90.6834, 98.2426, 102.794, 123.284, 135.601, 153.886, 189.781}, 
					      {118.062, 119.348, 111.604, 112.464, 121.064, 114.846, 120.385, 119.249, 106.702}, 
					      {-200, -200, -200, -200, -200, -200, -200, -200, -200} };
	//  DUNE-SP Gaisser-Hillas
	const double GH_RS60cm_SP[4][9] = { {1.37378, 1.3634, 1.31054, 1.23488, 1.14697, 1.01977, 0.886863, 0.751005, 0.592496}, 
					     {113.764, 128.753, 122.512, 141.309, 140.16, 153.797, 170.915, 184.999, 199.248}, 
					     {81.3747, 78.791, 87.2706, 81.9593, 92.3303, 102.592, 110.304, 112.577, 107.575}, 
					     {-200, -200, -200, -200, -200, -200, -200, -200, -200} }; 
	const double GH_RS120cm_SP[4][9] = { {1.22881, 1.20776, 1.15355, 1.08087, 0.988751, 0.868487, 0.736578, 0.604445, 0.465248}, 
					     {120.126, 137.211, 129.695, 150.215, 151.926, 168.741, 199.556, 223.586, 260.437}, 
					     {120.445, 115.844, 127.995, 114.96, 130.093, 141.39, 147.55, 154.139, 136.948}, 
					     {-200, -200, -200, -200, -200, -200, -200, -200, -200} }; 
	const double GH_RS180cm_SP[4][9] = { {1.16447, 1.14188, 1.08141, 1.00912, 0.911832, 0.793711, 0.656118, 0.517022, 0.38575}, 
					     {120.862, 138.321, 129.506, 152.468, 154.87, 171.04, 210.579, 240.266, 297.42}, 
					     {156.572, 146.229, 173.181, 147.513, 165.223, 175.133, 182.79, 211.805, 173.369}, 
					     {-200, -200, -200, -200, -200, -200, -200, -200, -200} };
	//  DUNE-DP Gaisser-Hillas
	const double GH_RS60cm_DP[4][9] = { {1.2378, 1.24291, 1.20084, 1.13647, 1.04805, 0.928209, 0.81468, 0.687154, 0.538787}, 
					     {95.9886, 105.046, 114.902, 121.08, 126.533, 142.666, 143.314, 156.796, 159.649}, 
					     {170.762, 161.485, 146.444, 136.313, 128.357, 112.543, 106.582, 88.2847, 81.1439}, 
					     {-200, -200, -200, -200, -200, -200, -200, -200, -200} };
	const double GH_RS120cm_DP[4][9] = { {1.15393, 1.13664, 1.09137, 1.02059, 0.927886, 0.800404, 0.675405, 0.544569, 0.410451}, 
					     {120.77, 131.444, 136.192, 143.061, 153.195, 168.244, 176.448, 189.741, 209.645}, 
					     {250.962, 234.726, 217.918, 202.817, 185.121, 167.992, 156.201, 133.598, 111.581}, 
					     {-200, -200, -200, -200, -200, -200, -200, -200, -200} }; 
	const double GH_RS180cm_DP[4][9] = { {1.11356, 1.09376, 1.04515, 0.969378, 0.873762, 0.741592, 0.613477, 0.476078, 0.344799}, 
					     {129.151, 146.954, 147.726, 155.173, 166.108, 180.721, 193.948, 206.119, 240.13}, 
					     {321.246, 291.538, 278.321, 258.852, 234.141, 219.926, 197.511, 175.507, 139.981}, 
					     {-200, -200, -200, -200, -200, -200, -200, -200, -200} };

	// bin size of the offset angle used in the parametrization in degrees
	const double delta_angulo = 10.;
	// LAr absorption length in cm
	const double L_abs = 2000.;


	// *************************************************************************************************
	//                    NUMBER OF VIS HITS PARAMETRIZATION
	// *************************************************************************************************
	// properties of LAr for reflected light path calculation
	// refractive indices in LAr
	const double n_LAr_VUV = 2.632;     // effective index due to group vel.
	const double n_LAr_vis = 1.23;

	// properties of detector
	// plane depth
	// Dune
	const double plane_depth = 363.38405;	// cm
	// Cathode plane covered by foils
	// Dune [ UPDATE ACTUAL VALUES ]
	// size
	double height_foils = 1200;	// cm
	double width_foils = 1400;	// cm
	// centre coordinates
	double x_foils = 363.38405; double y_foils = 0; double z_foils = 700;	// cm

	
public:	
	// constructor 
	solid_angle(const int detector_type, const int flagRS, const std::string flagDetector);

	// destructor
	~solid_angle(){};

	// hits calculating functions
	int VUVHits(int Nphotons_created, TVector3 ScintPoint, TVector3 OpDetPoint);
	int VisHits(int Nphotons_created, TVector3 ScintPoint, TVector3 OpDetPoint);

	// gaisser-hillas function
	static Double_t GaisserHillas(double *x, double *par);

	// solid angle of rectangular aperture calculation functions
	double omega(double a, double b, double d);
	double solid(acc& out, TVector3 v);

	// solid angle of spherical aperture calculation functions		TO DO

};

#endif