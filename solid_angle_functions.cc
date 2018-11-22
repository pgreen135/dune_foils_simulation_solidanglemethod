#include "solid_angle_functions.h"

// solid angle class implementation

#include <iostream>
#include <cmath>

#include "TRandom.h"

using namespace std;

// constructor
solid_angle::solid_angle(const int optical_detector, const int scattering, const string flagDet): optical_detector_type{optical_detector}, flagRS{scattering}, flagDetector{flagDet} {

  // initialise gaisser hillas functions for VUV Rayleigh scattering correction
  std::cout<<"Loading the parameters for the geometry corrections depending on the detector and raileigh Scattering length ..."<<std::endl;
  std::cout<<"Simulation made with a <RS> = "<<60.*(flagRS)<<" cm"<<std::endl;
  double pars_ini[4] = {0,0,0,0};
  // For SBN-like (size) detectors: SBND, MicroBooNE and ICARUS
  if(flagDetector == "SBN") {
    std::cout<<"Light simulation for SBN detector (SBND, MicroBooNE and ICARUS)"<<std::endl;
    for(int bin = 0; bin < 9; bin++) {
      GH[bin] =  new TF1("GH",GaisserHillas,0.,2000,4);      
      if(flagRS == 1) {
        for(int j=0; j < 4; j++) {
          pars_ini[j] = GH_RS60cm_SBN[j][bin];
        }
      }
      else if(flagRS == 2) { 
        for(int j=0; j < 4; j++) {
          pars_ini[j] = GH_RS120cm_SBN[j][bin];
        }
      }
      else if(flagRS == 3) {
        for(int j=0; j < 4; j++) {
          pars_ini[j] = GH_RS180cm_SBN[j][bin];
        }
      }  
      else {
        cout << "Error in Gaissr-Hillas settings!"<<endl; 
        break;
      }  
      GH[bin]->SetParameters(pars_ini);
    }
  }
  // For DUNE Single Phase like detector
  if(flagDetector == "SP") {
    std::cout<<"Light simulation for DUNE Single Phase detector"<<std::endl;
    for(int bin = 0; bin < 9; bin++) {
      GH[bin] =  new TF1("GH",GaisserHillas,0.,2000,4);      
      if(flagRS == 1) {
        for(int j=0; j < 4; j++) {
          pars_ini[j] = GH_RS60cm_SP[j][bin];
        }
      }
      else if(flagRS == 2) {
        for(int j=0; j < 4; j++) {
          pars_ini[j] = GH_RS120cm_SP[j][bin];
        }
      }
      else if(flagRS == 3) {
        for(int j=0; j < 4; j++) {
          pars_ini[j] = GH_RS180cm_SP[j][bin];  
        }
      }
      else {
        cout << "Error in Gaissr-Hillas settings!"<<endl; 
        break;
      }  
      GH[bin]->SetParameters(pars_ini);
    }
  }
  // For DUNE Dual Phase like detector
  if(flagDetector == "DP") {
    std::cout<<"Light simulation for DUNE Dual Phase detector"<<std::endl;
    for(int bin = 0; bin < 9; bin++) {
      GH[bin] =  new TF1("GH",GaisserHillas,0.,2000,4);      
      if(flagRS == 1) {
        for(int j=0; j < 4; j++) {
          pars_ini[j] = GH_RS60cm_DP[j][bin];
        }
      }
      else if(flagRS == 2) {
        for(int j=0; j < 4; j++) {
          pars_ini[j] = GH_RS120cm_DP[j][bin];
        }
      }
      else if(flagRS == 3) {
        for(int j=0; j < 4; j++) {
          pars_ini[j] = GH_RS180cm_DP[j][bin];  
        }
      }
      else {
        cout << "Error in Gaissr-Hillas settings!"<<endl; 
        break;
      }  
      GH[bin]->SetParameters(pars_ini);
    }
  }

}


// VUV hits calculation
int solid_angle::VUVHits(int Nphotons_created, TVector3 ScintPoint, TVector3 OpDetPoint) {
  gRandom->SetSeed(0);

  // set Arapuca geometry struct for solid angle function
  acc detPoint; 
  detPoint.ax = OpDetPoint[0]; detPoint.ay = OpDetPoint[1]; detPoint.az = OpDetPoint[2];  // centre coordinates of optical detector
  detPoint.w = width; detPoint.h = height; // width and height in cm of arapuca active window

  // distance and angle between ScintPoint and OpDetPoint
  double distance = sqrt(pow(ScintPoint[0] - OpDetPoint[0],2) + pow(ScintPoint[1] - OpDetPoint[1],2) + pow(ScintPoint[2] - OpDetPoint[2],2));
  double cosine = sqrt(pow(ScintPoint[0] - OpDetPoint[0],2)) / distance;
  double theta = acos(cosine)*180./pi;

  // get scintillation point coordinates relative to arapuca window centre
  TVector3 ScintPoint_rel = ScintPoint - OpDetPoint;  

  // calculate solid angle
  double solid_angle = solid_angle::solid(detPoint, ScintPoint_rel);

  // calculate number of photons hits by geometric acceptance: accounting for solid angle and LAr absorbtion length
  double hits_geo = exp(-1.*distance/L_abs) * (solid_angle / (4*pi)) * Nphotons_created;

  // apply Gaisser-Hillas correction for Rayleigh scattering distance and angular dependence
  // offset angle bin
  int j = std::round(theta/delta_angulo);
  double hits_rec = gRandom->Poisson( GH[j]->Eval(distance)*hits_geo/cosine );

  // round to integer value, cannot have non-integer number of hits
  int hits_vuv = std::round(hits_rec);

  return hits_vuv;
}

// Visible hits calculation
int solid_angle::VisHits(int Nphotons_created, TVector3 ScintPoint, TVector3 OpDetPoint) {
  gRandom->SetSeed(0);

  // 1). calculate total number of hits of VUV photons on reflective foils via solid angle + Gaisser-Hillas corrections:

  // set cathode plane struct for solid angle function
  acc cathode_plane; 
  cathode_plane.ax = x_foils; cathode_plane.ay = y_foils; cathode_plane.az = z_foils;   // centre coordinates of cathode plane
  cathode_plane.w = width_foils; cathode_plane.h = height_foils;                        // width and height in cm

  // get scintpoint coords relative to centre of cathode plane
  TVector3 cathodeCentrePoint(x_foils,y_foils,z_foils);
  TVector3 ScintPoint_relative = ScintPoint - cathodeCentrePoint; 

  // calculate distance and angle between ScintPoint and cathode centre point
  double distance_cathode = sqrt(pow(ScintPoint[0] - cathodeCentrePoint[0],2) + pow(ScintPoint[1] - cathodeCentrePoint[1],2) + pow(ScintPoint[2] - cathodeCentrePoint[2],2));
  double cosine_cathode = sqrt(pow(ScintPoint[0] - cathodeCentrePoint[0],2)) / distance_cathode;
  double theta_cathode = acos(cosine_cathode)*180./pi;

  // calculate solid angle
  double solid_angle_cathode = solid_angle::solid(cathode_plane, ScintPoint_relative);

  // calculate hits on cathode plane via geometric acceptance
  double cathode_hits_geo = exp(-1.*distance_cathode/L_abs) * (solid_angle_cathode / (4*pi)) * Nphotons_created;

  // apply Gaisser-Hillas correction for Rayleigh scattering distance and angular dependence
  // offset angle bin
  int j = std::round(theta_cathode/delta_angulo);
  double cathode_hits_rec = GH[j]->Eval(distance_cathode)*cathode_hits_geo/cosine_cathode;
  //double cathode_hits_rec = gRandom->Poisson( GH[j]->Eval(distance)*cathode_hits_geo/cosine_cathode );      // SWAP TO THIS ?


  // 2). calculate number of these hits which reach the Arapuca via solid angle 

  // set Arapuca geometry struct for solid angle function
  acc detPoint; 
  detPoint.ax = OpDetPoint[0]; detPoint.ay = OpDetPoint[1]; detPoint.az = OpDetPoint[2];  // centre coordinates of optical detector
  detPoint.w = width; detPoint.h = height; // width and height in cm of arapuca active window

  // calculate hotspot location  
  TVector3 v_to_wall(plane_depth - ScintPoint[0],0,0);        
  TVector3 hotspot = ScintPoint + v_to_wall;

  // get hotspot coordinates relative to detpoint
  TVector3 emission_relative = hotspot - OpDetPoint;

  // calculate solid angle
  double solid_angle_detector = solid_angle::solid(detPoint, emission_relative);

  // calculate number of hits via geometeric acceptance
  double hits_geo = (solid_angle_detector / (4*pi)) * cathode_hits_rec;

  // apply some form of correction

  // round final result
  int hits_vis = std::round(hits_geo);

  return hits_vis;
}


// gaisser-hillas function definition
Double_t solid_angle::GaisserHillas(double *x,double *par) {
  //This is the Gaisser-Hillas function
  Double_t X_mu_0=par[3];
  Double_t Normalization=par[0];
  Double_t Diff=par[1]-X_mu_0;
  Double_t Term=pow((*x-X_mu_0)/Diff,Diff/par[2]);
  Double_t Exponential=TMath::Exp((par[1]-*x)/par[2]);
  
  return ( Normalization*Term*Exponential);
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
