#include "CRYGenerator.h"
#include "CRYSetup.h"
#include "CRYParticle.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <string>

// ROOT
#include "TFile.h"
#include "TTree.h"
#include "TRandom.h"

double interp1d(const std::vector<double>&x,const std::vector<double>& y,double xq);

int main() {

    std::ifstream input("setup.file");

    std::stringstream buffer;
    buffer << input.rdbuf();

    CRYSetup setup(buffer.str(),
                   "/Users/changhyonha/Software/cry_v1.7/data");
    CRYGenerator gen(&setup);

    // tunnel profile
    std::vector<double> station;
    std::vector<double> mwe;
    std::vector<double> rock_m;
    std::vector<double> water_m;
    
    std::ifstream prof("tunnel_profile.txt");

    if(!prof.is_open()) {
      std::cerr << "ERROR: Cannot open tunnel_profile.txt" << std::endl;
      return 1;
    }
    
    std::string line;
    double s_val, mwe_val, rock_val, water_val;
    
    while (std::getline(prof, line)) {
      
      if(line.empty()) continue;
      if(line[0] == '#') continue;
      
      std::stringstream ss(line);
      
      if(ss >> s_val >> mwe_val >> rock_val >> water_val) {
        station.push_back(s_val);
        mwe.push_back(mwe_val);
        rock_m.push_back(rock_val);
        water_m.push_back(water_val);
      }
    }
    
    std::cout << "Loaded profile points: "
	      << station.size() << std::endl;

    if(station.size() == 0) {
      std::cerr << "ERROR: No profile points loaded." << std::endl;
      return 1;
    }

    std::vector<CRYParticle*> particles;

    TFile *fout = new TFile("muons.root", "RECREATE");
    TTree *tree = new TTree("mu", "CRY muons");

    // variables
    double E;
    double theta;
    double phi;

    double u,v,w;

    double Emin;
    int survive;


    double x1, y1;
    double x2, y2;
    int coincidence;
    int detected;

    // tunnel variables
    double s_m;
    double X0;
    double Xslant;

    double rock;
    double water;

 
    // branches
    tree->Branch("E", &E, "E/D");
    tree->Branch("theta", &theta, "theta/D");
    tree->Branch("phi", &phi, "phi/D");

    tree->Branch("u", &u, "u/D");
    tree->Branch("v", &v, "v/D");
    tree->Branch("w", &w, "w/D");

    tree->Branch("Emin", &Emin, "Emin/D");
    tree->Branch("survive", &survive, "survive/I");

    tree->Branch("x1", &x1, "x1/D");
    tree->Branch("y1", &y1, "y1/D");
    tree->Branch("x2", &x2, "x2/D");
    tree->Branch("y2", &y2, "y2/D");
    tree->Branch("coincidence", &coincidence, "coincidence/I");
    tree->Branch("detected", &detected, "detected/I");

    tree->Branch("s_m", &s_m, "s_m/D");
    tree->Branch("X0", &X0, "X0/D");
    tree->Branch("Xslant", &Xslant, "Xslant/D");

    tree->Branch("rock_m", &rock, "rock_m/D");
    tree->Branch("water_m", &water, "water_m/D");

    // detector geometry
    double Lx  = 1.0;   // m, along tunnel
    double Ly  = 0.5;   // m, across tunnel
    double gap = 0.5;   // m, vertical separation

    //    int Nevents = 100000;
    int Nevents = 5000;

    std::ofstream foutRate("rate_profile.txt");
    foutRate << "# s_m X0_mwe X0_shifted_mean_mwe ratio rate_Hz Ncoin_surface Ncoin_underground\n";

    for(size_t ip=0; ip<mwe.size(); ip++) {
      
      s_m = station[ip];
      X0  = mwe[ip];
      rock = rock_m[ip];
      water = water_m[ip];

      long Ncoin_surface = 0;
      long Ncoin_underground = 0;
      
      double X0_shifted_sum = 0.0;
      long Nshift = 0;
 
      for(int iev=0; iev<Nevents; iev++) {
	
        particles.clear();
	
        gen.genEvent(&particles);
	
        for(auto p : particles){
	  
	  E = p->ke()/1000.0; // GeV
	  
	  u = p->u();
	  v = p->v();
	  w = p->w();
	  
	  // downward only
	  if(w >= 0) continue;
	  
	  double costh = -w;
	  theta = acos(costh);
	  phi = atan2(v,u);
	  
	  theta *= 180.0/M_PI;
	  phi *= 180.0/M_PI;
	  
	  // sample hit position on upper panel
	  x1 = (gRandom->Rndm() - 0.5) * Lx;
	  y1 = (gRandom->Rndm() - 0.5) * Ly;
	  
	  
	  x2 = x1 + gap * u / costh;
	  y2 = y1 + gap * v / costh;
	  
	  // check lower-panel hit
	  coincidence = (fabs(x2) < Lx/2.0) && (fabs(y2) < Ly/2.0);
	  if (coincidence) Ncoin_surface++;
	  
	  //	  Xslant = X0 / costh;
	  double h_eff = rock_m[ip] + water_m[ip];   // meters, physical height
	  double ds = h_eff * u / costh;
	  double s_shift = s_m + ds;
	  double X0_shifted = interp1d(station, mwe, s_shift);
	  Xslant = X0_shifted / costh;
	  Emin = 0.2 * Xslant;
	  
	  survive = (E > Emin);

	  if(coincidence){
	    X0_shifted_sum += X0_shifted;
	    Nshift++;
	  }
	  // underground survived AND two-panel coincidence
	  detected = survive && coincidence;
	  if (detected) Ncoin_underground++;
	  
	  tree->Fill();
        }
      }
      double X0_shifted_mean = 0.0;
      if(Nshift > 0) {
	X0_shifted_mean = X0_shifted_sum / Nshift;
      }

      double ratio = 0.0;
      if(Ncoin_surface > 0) {
        ratio = double(Ncoin_underground) / double(Ncoin_surface);
      }
      
      double rateHz = 36.5 * ratio;
      
      foutRate
        << s_m << " "
        << X0 << " "
        << X0_shifted_mean << " "
        << ratio << " "
        << rateHz << " "
        << Ncoin_surface << " "
        << Ncoin_underground << "\n";

      std::cout
        << "s = " << s_m
        << " m, X0 = " << X0
        << " m.w.e., rate = " << rateHz
        << " Hz" << std::endl;
    }
    
    
    tree->Write();
    fout->Close();
    std::cout << "Saved muons.root" << std::endl;

    foutRate.close();
    
    return 0;
}

double interp1d(const std::vector<double>& x,const std::vector<double>& y,double xq){
  if(xq <= x.front()) return y.front();
  if(xq >= x.back())  return y.back();
  
  for(size_t i=0; i<x.size()-1; i++) {
    if(xq >= x[i] && xq < x[i+1]) {
      double t = (xq - x[i]) / (x[i+1] - x[i]);
      return y[i]*(1.0-t) + y[i+1]*t;
    }
  }
  return y.back();
}
