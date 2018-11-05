#include "solid_angle_functions.h"

// solid angle functions
// calculates number of pmt_hits for VUV light using solid angle method

std::vector<double> solid_angle::SolidAngleAnalyzer(int Nphotons_created, const double quantum_efficiency, const double catcov, const double vuvfrac, const double visfrac, int pmt_number, TVector3 ScintPoint, TVector3 OpDetPoint) {
	
	std::vector<double> pmt_hits;

	const double pi = 3.14159265358979323846;

	// -------- VUV photons --------

	// needs Rayleigh scattering correction

	// set up arapuca info struct for solid angle function
	solid_angle::acc detPoint; 
	detPoint.ax = OpDetPoint[0]; detPoint.ay = OpDetPoint[1]; detPoint.az = OpDetPoint[2];	// centre coordinates
	// dimensions set to match those in diegos gdml file
  detPoint.w = 10.16; detPoint.h = 10.16;	// width and height in cm of arapuca active window     ARAPUCA WINDOW DIMENSIONS: width = 7.8cm, height = 9.8cm

	// get scintillation point coordinates relative to arapuca window centre
	TVector3 ScintPoint_rel = ScintPoint - OpDetPoint;	

	// calculate solid angle
	double solid_angle = solid_angle::solid(detPoint, ScintPoint_rel);

	// calculate number of photons incident on pmt
	double hits_vuv = (solid_angle / (4*pi)) * Nphotons_created;
  int int_hits_vuv = hits_vuv;

	// determine number of incident photons detected 
  int total_hits_vuv = 0;
  for(int i = 0; i < int_hits_vuv; i++){
    if(gRandom->Uniform(1.) <= vuvfrac*quantum_efficiency*0.5){total_hits_vuv++;}		// extra factor 1/2 due to only half vuv photons incident on TPB coated detector are detected, other half remitted in opposite direction as visible photons
  }
  
	// push result back to vector 
	pmt_hits.push_back(total_hits_vuv);

	// -------- VIS photons --------

	// to do

	int total_hits_vis = 0;
	pmt_hits.push_back(total_hits_vis);

	return pmt_hits;
}



// solid angle of rectanglular aperture - arapuca - calculation functions
// written by Franciole Marinho
double solid_angle::omega(double a, double b, double d){

  double aa = a/(2.0*d);
  double bb = b/(2.0*d);
  double aux = (1+aa*aa+bb*bb)/((1.+aa*aa)*(1.+bb*bb));
  return 4*std::acos(std::sqrt(aux));

}

double solid_angle::solid(solid_angle::acc& out, TVector3 v){

  //This function gives the solid angle for a segment adapted to tallbo geometry

  //v is the position of the track segment with respect to 
  //the center position of the arapuca window 
  //That is the: segment position wrt tallbo frame - center position of the 
  //arapuca window wrt to tallbo frame

  //out is only a struct that I use to pass info regarding the arapuca 
  //dimensions and position

  // arapuca plane fixed in x direction	

  if( v.Y()==0.0 && v.Z()==0.0){
    return omega(out.w,out.h,v.X());
  }
  
  if( (std::abs(v.Y()) > out.w/2.0) && (std::abs(v.Z()) > out.h/2.0)){
    double A, B, a, b, d;
    A = std::abs(v.Y())-out.w/2.0;
    B = std::abs(v.Z())-out.h/2.0;
    a = out.w;
    b = out.h;
    d = abs(v.X());
    double to_return = (omega(2*(A+a),2*(B+b),d)-omega(2*A,2*(B+b),d)-omega(2*(A+a),2*B,d)+omega(2*A,2*B,d))/4.0;
    return to_return;
  }
  
  if( (std::abs(v.Y()) <= out.w/2.0) && (std::abs(v.Z()) <= out.h/2.0)){
    double A, B, a, b, d;
    A = -std::abs(v.Y())+out.w/2.0;
    B = -std::abs(v.Z())+out.h/2.0;
    a = out.w;
    b = out.h;
    d = abs(v.X());
    double to_return = (omega(2*(a-A),2*(b-B),d)+omega(2*A,2*(b-B),d)+omega(2*(a-A),2*B,d)+omega(2*A,2*B,d))/4.0;
    return to_return;
  }

  if( (std::abs(v.Y()) > out.w/2.0) && (std::abs(v.Z()) <= out.h/2.0)){
    double A, B, a, b, d;
    A = std::abs(v.Y())-out.w/2.0;
    B = -std::abs(v.Z())+out.h/2.0;
    a = out.w;
    b = out.h;
    d = abs(v.X());
    double to_return = (omega(2*(A+a),2*(b-B),d)-omega(2*A,2*(b-B),d)+omega(2*(A+a),2*B,d)-omega(2*A,2*B,d))/4.0;
    return to_return;
  }

  if( (std::abs(v.Y()) <= out.w/2.0) && (std::abs(v.Z()) > out.h/2.0)){
    double A, B, a, b, d;
    A = -std::abs(v.Y())+out.w/2.0;
    B = std::abs(v.Z())-out.h/2.0;
    a = out.w;
    b = out.h;
    d = abs(v.X());
    double to_return = (omega(2*(a-A),2*(B+b),d)-omega(2*(a-A),2*B,d)+omega(2*A,2*(B+b),d)-omega(2*A,2*B,d))/4.0;
    return to_return;
  }
  // error message if none of these cases, i.e. something has gone wrong!
  std::cout << "Warning: invalid solid angle call." << std::endl;
  return 0.0;
}
