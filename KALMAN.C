#include "TCanvas.h"
#include "TGraph.h"
#include <iostream>

#include <vector>
#include <TMath.h>
#include <math.h>

#include <fstream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <string>

#include <TMatrixD.h>

#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TClonesArray.h"

#include <algorithm>
#include <limits>


using namespace std;

// Element-wise addition of two vectors
std::vector<double> vectorAdd(const std::vector<double>& a, const std::vector<double>& b) {
    std::vector<double> result(a.size());
    for (size_t i = 0; i < a.size(); i++) {
        result[i] = a[i] + b[i];
    }
    return result;
}

// Scalar multiplication of a vector
std::vector<double> vectorScalarMultiply(const std::vector<double>& vec, double scalar) {
    std::vector<double> result(vec.size());
    for (size_t i = 0; i < vec.size(); i++) {
        result[i] = vec[i] * scalar;
    }
    return result;
}


// Function to calculate the derivatives

std::vector<double> derivatives(const std::vector<double>& x, double z, double Bx, double By, double Bz, double q_over_p, double kappa) {
    std::vector<double> dxdt(5);

    dxdt[0] = x[2]; // dx/dz = tx = px/pz
    dxdt[1] = x[3]; // dy/dz = ty = py/pz
    dxdt[2] = kappa * q_over_p * sqrt(x[2]*x[2] + x[3]*x[3] + 1) * x[3] * Bz - (1 + x[2]*x[2]) * By + x[2] * x[3] * Bx;
    dxdt[3] = kappa * q_over_p * sqrt(x[2]*x[2] + x[3]*x[3] + 1) * (-x[2]) * Bz - (1 + x[3]*x[3]) * Bx - x[2] * x[3] * By;
    dxdt[4] = 0.0; // d(q/p)/dz = 0
   // std::cout<<dxdt[3]<<endl;
    return dxdt;
}



void KALMAN() {

     TCanvas *c1 = new TCanvas("c1", "My Graph", 800, 600);
     c1->Divide(2, 2);
     c1->Draw();


     TH2F *rx_vs_ry = new TH2F("rx_vs_ry", "rungex_vs_rungey", 720, 0, -3, 1000, 0, 2.0);

    // Initial conditions and parameters
     Double_t px,py,pz;
     Double_t q = 1.6022*TMath::Power(10,-19);   //charge of the particle in C
     px = 1;
     py = 1;
     pz = 2.0;

     // Calculate the momentum magnitude
     Double_t momentum = sqrt(px * px + py * py + pz * pz);                     //momentum in MeV/c.

     // Print the momentum magnitude
     //std::cout << "Momentum Magnitude: " << momentum << " MeV/c" << std::endl;


     std::vector<double> x0 = {0.0, 0.0, px/pz, py/pz, q/momentum};  // initial values for x, y, tx, ty, q/p
     Double_t z0 = 0.0; // initial z position
     Double_t Bx = 0.0; // magnetic field component Bx
     Double_t By = 0.0; // magnetic field component By
     Double_t Bz = 3.0;  // magnetic field component Bz in Tesla
     Double_t q_over_p = q/momentum; // q/p value
     Double_t kappa = 0.299792458 * 0.001; // kappa value in MeV/c.m.T
     Double_t z_target = 200.0; // target z position in m
     Double_t h = 1; // step size     in m. 


  //   std::cout<< h  <<  endl;

    // Perform RK4 integration
     while (z0 < z_target){
           // Calculate the RK4 intermediate steps
           std::vector<double> k1 = derivatives(x0, z0, Bx, By, Bz, q_over_p, kappa);
           std::vector<double> k2 = derivatives(vectorAdd(x0, vectorScalarMultiply(k1, 0.5 * h)), z0 + 0.5 * h, Bx, By, Bz, q_over_p, kappa);
           std::vector<double> k3 = derivatives(vectorAdd(x0, vectorScalarMultiply(k2, 0.5 * h)), z0 + 0.5 * h, Bx, By, Bz, q_over_p, kappa);
           std::vector<double> k4 = derivatives(vectorAdd(x0, vectorScalarMultiply(k3, h)), z0 + h, Bx, By, Bz, q_over_p, kappa);

           //std::cout<< k1[0] << " " << k4[1] << " " << k4[2] << " " <<k4[3]  << " " << k4[4] << endl;

           // Update the state vector for the next step
           for (int i = 0; i < 5; i++) {
               x0[i] = x0[i] + (h / 6.0) * (k1[i] + 2.0 * k2[i] + 2.0 * k3[i] + k4[i]);
              // std::cout << x0[i] << endl << endl;
              // rx_vs_ry->Fill(x[0],x[1]);
           }
          // std::cout << x0[4] <<" " <<  x0[3]  << endl;
           rx_vs_ry->Fill(x0[0],x0[1]);

           z0 +=h;
     }

     // Plot the histogram
     rx_vs_ry->SetMarkerStyle(22);
     rx_vs_ry->SetMarkerSize(0.5);
     rx_vs_ry->Draw();

     return 0;
}

