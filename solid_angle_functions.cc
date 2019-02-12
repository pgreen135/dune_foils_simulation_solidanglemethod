#include "solid_angle_functions.h"

// solid angle class implementation

#include <iostream>
#include <cmath>

#include "TRandom.h"
#include "TSystem.h"
#include "TMath.h"
#include "TFormula.h"
#include "Math/SpecFuncMathMore.h"

using namespace std;

// constructor
solid_angle::solid_angle(const int optical_detector, const int scattering, const string flagDet): optical_detector_type{optical_detector}, flagRS{scattering}, flagDetector{flagDet} {

  // load mathmore library
  gSystem->Load("libMathMore.so");
  if(gSystem->Load("libMathMore.so") < 0) {
      throw(std::runtime_error("Unable to load MathMore library"));
    }
  _mathmore_loaded_ = true;

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
      // uses the values for Arapucas, as these replace pmt values in solid_angle_functions.h; would want to add switch between PMT or Arapuca configs
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


  // initialise pol5 functions for visible hits correction -- arapucas, RS60cm only [preliminary]
  double pars_ini_vis[6] = {0,0,0,0,0,0};
  for (int bin = 0; bin < 9; bin++) {
    VIS_pol[bin] = new TF1 ("pol", "pol5", 0, 2000);
    for (int j = 0; j < 6; j++){
      pars_ini_vis[j] = VIS_RS60cm_SP[j][bin];
    }
    VIS_pol[bin]->SetParameters(pars_ini_vis);
  }


}


// VUV hits calculation
int solid_angle::VUVHits(int Nphotons_created, TVector3 ScintPoint, TVector3 OpDetPoint) {
  gRandom->SetSeed(0);

  // distance and angle between ScintPoint and OpDetPoint
  double distance = sqrt(pow(ScintPoint[0] - OpDetPoint[0],2) + pow(ScintPoint[1] - OpDetPoint[1],2) + pow(ScintPoint[2] - OpDetPoint[2],2));
  double cosine = sqrt(pow(ScintPoint[0] - OpDetPoint[0],2)) / distance;
  double theta = acos(cosine)*180./pi;

  // calculate solid angle:
  double solid_angle;
  // Arapucas
  if (optical_detector_type == 0) {
    // set Arapuca geometry struct for solid angle function
    acc detPoint; 
    detPoint.ax = OpDetPoint[0]; detPoint.ay = OpDetPoint[1]; detPoint.az = OpDetPoint[2];  // centre coordinates of optical detector
    detPoint.w = width; detPoint.h = height; // width and height in cm of arapuca active window

    // get scintillation point coordinates relative to arapuca window centre
    TVector3 ScintPoint_rel = ScintPoint - OpDetPoint;  

    // calculate solid angle
    solid_angle = solid_angle::solid(detPoint, ScintPoint_rel);
  }
  // PMTs
  else if (optical_detector_type == 1) {
    // offset in z-y plane
    double d = sqrt(pow(ScintPoint[1] - OpDetPoint[1],2) + pow(ScintPoint[2] - OpDetPoint[2],2));
    // drift distance (in x)
    double h =  sqrt(pow(ScintPoint[0] - OpDetPoint[0],2));
    // Solid angle of a disk
    solid_angle = Disk_SolidAngle(d, h, radius);
  }
  else {
    std::cout << "Error: Invalid optical detector type." << endl;
    exit(1);
  }  

  // calculate number of photons hits by geometric acceptance: accounting for solid angle and LAr absorbtion length
  double hits_geo = exp(-1.*distance/L_abs) * (solid_angle / (4*pi)) * Nphotons_created;

  // apply Gaisser-Hillas correction for Rayleigh scattering distance and angular dependence
  // offset angle bin
  int j = (theta/delta_angulo);
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

  // calculate solid angle
  double solid_angle_cathode = solid_angle::solid(cathode_plane, ScintPoint_relative);

  // calculate distance and angle between ScintPoint and hotspot
  // vast majority of hits in hotspot region directly infront of scintpoint,therefore consider attenuation for this distance and on axis GH instead of for the centre coordinate
  double distance_cathode = plane_depth - ScintPoint[0];
  double cosine_cathode = 1;
  double theta_cathode = 0;

  // calculate hits on cathode plane via geometric acceptance
  double cathode_hits_geo = exp(-1.*distance_cathode/L_abs) * (solid_angle_cathode / (4.*pi)) * Nphotons_created;

  // apply Gaisser-Hillas correction for Rayleigh scattering distance and angular dependence
  // offset angle bin
  int j = (theta_cathode/delta_angulo);
  double cathode_hits_rec = GH[j]->Eval(distance_cathode)*cathode_hits_geo/cosine_cathode;
  //double cathode_hits_rec = gRandom->Poisson( GH[j]->Eval(distance)*cathode_hits_geo/cosine_cathode );      // Swap to this? 


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
  // ##################################################################
  // NOTE: there was an error here in previous version (if you had already copied this code), was using 4pi instead of 2pi i.e. forgetting about the vm2000 foils (doh!)
  // required corrections are now factor 2 less due to this!

  double hits_geo = (solid_angle_detector / (2*pi)) * cathode_hits_rec;

  // ##################################################################

  // distance to hotspot
  double distance_vuv = sqrt(pow(ScintPoint[0] - hotspot[0],2) + pow(ScintPoint[1] - hotspot[1],2) + pow(ScintPoint[2] - hotspot[2],2));
  // distance from hotspot to arapuca
  double distance_vis = sqrt(pow(hotspot[0] - OpDetPoint[0],2) + pow(hotspot[1] - OpDetPoint[1],2) + pow(hotspot[2] - OpDetPoint[2],2));
  // angle between hotspot and arapuca
  double cosine_vis = sqrt(pow(hotspot[0] - OpDetPoint[0],2)) / distance_vis;
  double theta_vis = acos(cosine_vis)*180./pi;

  // apply correction
  int k = (theta_vis/delta_angle);
  double hits_rec = gRandom->Poisson(VIS_pol[k]->Eval(distance_vuv)*hits_geo/cosine_vis);

  // round final result
  int hits_vis = std::round(hits_rec);
  
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


// solid angle of rectanglular aperture calculation functions

double solid_angle::omega(double a, double b, double d){

  double aa = a/(2.0*d);
  double bb = b/(2.0*d);
  double aux = (1+aa*aa+bb*bb)/((1.+aa*aa)*(1.+bb*bb));
  return 4*std::acos(std::sqrt(aux));

}

double solid_angle::solid(solid_angle::acc& out, TVector3 v){

  //v is the position of the track segment with respect to 
  //the center position of the arapuca window 
 
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


// solid angle of circular aperture
double solid_angle::Disk_SolidAngle(double* x, double *p) {
  const double d = x[0];
  const double h = x[1];
  const double b = p[0];
  if(b <= 0. || d < 0. || h <= 0.) return 0.; 
  const double aa = TMath::Sqrt(h*h/(h*h+(b+d)*(b+d)));
  if(d == 0) {
    return 2.*TMath::Pi()*(1.-aa);
  }
  const double bb = TMath::Sqrt(4*b*d/(h*h+(b+d)*(b+d)));
  const double cc = 4*b*d/((b+d)*(b+d));

  if(!_mathmore_loaded_) {
    if(gSystem->Load("libMathMore.so") < 0) {
      throw(std::runtime_error("Unable to load MathMore library"));
    }
    _mathmore_loaded_ = true;
  }
  if(TMath::Abs(ROOT::Math::comp_ellint_1(bb) - bb) < 1e-10 && TMath::Abs(ROOT::Math::comp_ellint_3(cc,bb) - cc) <1e-10) {
    throw(std::runtime_error("please do gSystem->Load(\"libMathMore.so\") before running Disk_SolidAngle for the first time!"));
  }
  if(d < b) {
    return 2.*TMath::Pi() - 2.*aa*(ROOT::Math::comp_ellint_1(bb) + TMath::Sqrt(1.-cc)*ROOT::Math::comp_ellint_3(cc,bb));
  }
  if(d == b) {
    return TMath::Pi() - 2.*aa*ROOT::Math::comp_ellint_1(bb);
  }
  if(d > b) {
    return 2.*aa*(TMath::Sqrt(1.-cc)*ROOT::Math::comp_ellint_3(cc,bb) - ROOT::Math::comp_ellint_1(bb));
  }

  return 0.;
}

double solid_angle::Disk_SolidAngle(double d, double h, double b) {
  double x[2] = { d, h };
  double p[1] = { b };
  if(!_mathmore_loaded_) {
    if(gSystem->Load("libMathMore.so") < 0) {
      throw(std::runtime_error("Unable to load MathMore library"));
    }
    _mathmore_loaded_ = true;
  }
  return Disk_SolidAngle(x,p);
}