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
std::vector<Double_t> vectorAdd(const std::vector<Double_t>& a, const std::vector<Double_t>& b) {
    std::vector<Double_t> result(a.size());
    for (size_t t = 0; t < a.size(); t++) {
        result[t] = a[t] + b[t];
    }
    return result;
}

// Scalar multiplication of a vector
std::vector<Double_t> vectorScalarMultiply(const std::vector<Double_t>& vec, Double_t scalar) {
    std::vector<Double_t> result(vec.size());
    for (size_t t = 0; t < vec.size(); t++) {
        result[t] = vec[t] * scalar;
    }
    return result;
}

//Function for energy Loss Calculation
void energyloss(std::string file, TGraph* eLossCurve){
     std::string eLossFileName_;
     eLossFileName_ = file;
     std::ifstream elossfile;
     std::string line;

     Float_t ener = 0;
     Float_t stp = 0;

     try {

	 elossfile.open(eLossFileName_,ios::in);
	 if ((elossfile.peek() == std::ifstream::traits_type::eof())) {
            std::cout << " Error: Energy loss file not found! Exiting..."<< "\n";
      	    std::exit(EXIT_FAILURE);
    	 }

	//std::cout << " Processing energy loss data file " << eLossFileName_ << "\n";
	 for (auto i = 0; i < 2; ++i) {
             std::getline(elossfile,line);
    	 }
	 while (std::getline(elossfile, line)) {
	       std::istringstream data(line);
	       data >> ener >> stp ;
      	       eLossCurve->SetPoint(eLossCurve->GetN(),ener,stp);
    	 }

     } catch (...) {
     }

}

Double_t GetEnergy(Double_t vx, Double_t vy, Double_t vz){
         Double_t vv, bet, gamma1, Energy;
         Double_t c = 3.0*TMath::Power(10, 8);

         vv = TMath::Sqrt(TMath::Power((vx),2)+TMath::Power((vy),2)+TMath::Power((vz),2));
         bet = vv/c;
//cout<<bet<<endl;
         gamma1 = TMath::Sqrt(1/(1-TMath::Power(bet,2)));
         Energy = (gamma1-1)*931.494028;
         return Energy;
}

Double_t StoppingPower(Double_t Energy){
         Double_t stpSim;
         TGraph* eLossCurve = new TGraph();
         energyloss("stpHydrogen.txt", eLossCurve);
         Double_t  gasMediumDensity_ = 8.988e-5;   //g/cm3 at 1 bar

         stpSim = eLossCurve->Eval(Energy);
         stpSim*=1.6021773349e-13;
         stpSim*= gasMediumDensity_*100;
         stpSim/= 1.6726*TMath::Power(10,-27);
         return stpSim;
}

// Function to calculate the derivatives
std::vector<double> derivatives (Double_t t,const std::vector<double>& x ,double Ex, double Ey, double Ez, double Bx, double By, double Bz, double q_over_m) {
    Double_t Energy = GetEnergy(x[3], x[4], x[5]);
    Double_t st = StoppingPower(Energy);
  //  cout<<x[3]<<endl<<endl;
    std::vector<double> dxdt(6);
    Double_t rr,az,po;

    rr = TMath::Sqrt(TMath::Power(x[3],2)+TMath::Power(x[4],2)+TMath::Power(x[5],2));
    az = TMath::ATan2(x[4] , x[3]);
    po = TMath::ACos(x[5]/rr);

    dxdt[0] = x[3];
    dxdt[1] = x[4];
    dxdt[2] = x[5];
    dxdt[3] = q_over_m *(Ex+x[4]*Bz-x[5]*By) - st*TMath::Sin(po)*TMath::Cos(az); //dxdt with energyloss compensation.;
    dxdt[4] = q_over_m *(Ey+x[5]*Bx-x[3]*Bz) - st*TMath::Sin(po)*TMath::Sin(az);
    dxdt[5] = q_over_m *(Ez+x[3]*By-x[4]*Bx) - st*TMath::Cos(po);
   // std::cout<<dxdt[3]<<endl;
    return dxdt;
}


void KALMAN() {

     TCanvas *c1 = new TCanvas("c1", "My Graph", 800, 600);
     c1->Divide(2, 2);
     c1->Draw();

     TH2F *rx_vs_ry = new TH2F("rx_vs_ry", "rungex_vs_rungey", 720, 0, -3, 1000, 0, 2.0);
     TH2F *propagatorx_vs_propagatory = new TH2F("propagatorx_vs_propagatory", "propagatorx_vs_propagatory", 720, 0, -3, 1000, 0, 2.0);
     TH3F *F_projection = new TH3F("F_projection", "propagator_projection", 720, 5.0, -3, 100, 0, 10.0,100, -5, 10.0);

     Double_t Ex,Ey,Ez;              //Electric field in V/m.
     Double_t Bx,By,Bz;              // magnetic field in Tesla
     Double_t t0=0.0;                  //Initial time.
     Int_t n = 1000;

     Double_t x_pos[n], y_pos[n],z_pos[n];
     Double_t vx[n],vy[n],vz[n];

     Double_t q = 1.6022*TMath::Power(10,-19);   //charge of the particle in C
     Double_t m = 1.6726*TMath::Power(10,-27); 	// mass of the particle in kg
     Double_t q_over_m = q/m;                   // q/m value
     //cout << q_over_m << endl;
     Double_t B= 3.0;                 // Applied magnetic field (T).
     Double_t E=TMath::Cos((q*B)/m) * 500;   // Applied Electric field(V/m).
     Bx = 0;
     By = 0;                	// magnetic field in x and y  direction in Tesla.
     Bz = B;          // magnetic field in the z direction in Tesla.
     Ex = 0;
     Ey = 0;			// Electric field in the x and  direction in V/m.
     Ez = -E;                        // Electric field in the z direction.


// Define the initial conditions.
     std::vector<double> x0 = {0, 0, 0, 10e6, 10e6,10e6};  // initial values for x, y, z,vx,vy,vz
     std::vector<std::vector<double>> state_vectors;

//Define the time step
     Double_t h = TMath::Power(10,-10);     //in seconds.
     Double_t Energy = GetEnergy(x0[3],x0[4],x0[5]);
     Int_t i = 0;




/*----------------------------------------*/
cout<<endl;
cout<<endl;
cout<<"+----------------+"<<endl;
cout<<" Running kalman.C"<<endl;
cout<<"+----------------+"<<endl;

cout<<"Initial Energy value: "<< Energy << " MeV" <<endl;
cout<<endl;
cout<<endl;
/*----------------------------------------*/

    while (Energy > 0.1) {
          std::vector<double> k1 = derivatives(t0,x0, Ex, Ey, Ez, Bx, By, Bz, q_over_m);
          std::vector<double> k2 = derivatives(t0+0.5*h,vectorAdd(x0, vectorScalarMultiply(k1, 0.5 * h)), Ex, Ey, Ez, Bx, By, Bz, q_over_m);
          std::vector<double> k3 = derivatives(t0+0.5*h,vectorAdd(x0, vectorScalarMultiply(k2, 0.5 * h)), Ex, Ey, Ez, Bx, By, Bz, q_over_m);
          std::vector<double> k4 = derivatives(t0+h,vectorAdd(x0, vectorScalarMultiply(k3, h)), Ex, Ey, Ez, Bx, By, Bz, q_over_m);

          // Update the state vector for the next step
          x0[0] = x0[0] + (h / 6) * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]);
          x0[1] = x0[1] + (h / 6) * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1]);
          x0[2] = x0[2] + (h / 6) * (k1[2] + 2 * k2[2] + 2 * k3[2] + k4[2]);
          x0[3] = x0[3] + (h / 6) * (k1[3] + 2 * k2[3] + 2 * k3[3] + k4[3]);
          x0[4] = x0[4] + (h / 6) * (k1[4] + 2 * k2[4] + 2 * k3[4] + k4[4]);
          x0[5] = x0[5] + (h / 6) * (k1[5] + 2 * k2[5] + 2 * k3[5] + k4[5]);

          // Update the time
          t0 = t0 + h;
          //std::cout<<x0[0]<<" "<<x0[1]<<" "<<x0[2]<<endl;
          state_vectors.push_back(x0);
          rx_vs_ry->Fill(x0[0], x0[1]);
          Energy = GetEnergy(x0[3], x0[4], x0[5]);

          i += 1;
          //std::cout << Energy << endl << endl;
    }

    Int_t numStates = state_vectors.size(); // Number of state vectors
    Int_t numDimensions = state_vectors[0].size(); // Number of dimensions in each state vector

    std::cout<< numStates << endl;
    TMatrixD state_matrix(numDimensions, numStates); // Create the state matrix

    // Populate the state matrix using the values from the state_vectors
    for (int i = 0; i < numStates; i++) {
        for (int j = 0; j < numDimensions; j++) {
            state_matrix(j, i) = state_vectors[i][j];
        }
    }

    // Initialize the Jacobian matrix
    TMatrixD jacobian(6, 6);

    // Calculate the derivatives at the initial state
    std::vector<double> dxdt = derivatives(t0, x0, Ex, Ey, Ez, Bx, By, Bz, q_over_m);

    // Calculate the Jacobian matrix
    for (int i = 0; i < 6; i++) {
        // Perturb the i-th component of x0 by a small amount
        double epsilon = 1e-6;
        std::vector<double> x0_perturbed = x0;
        x0_perturbed[i] += epsilon;

        // Calculate the perturbed derivatives
        std::vector<double> dxdt_perturbed = derivatives(t0, x0_perturbed, Ex, Ey, Ez, Bx, By, Bz, q_over_m);

        // Calculate the i-th column of the Jacobian matrix
        for (int j = 0; j < 6; j++) {
            jacobian(j, i) = (dxdt_perturbed[j] - dxdt[j]) / epsilon;
        }
    }
    //jacobian.Print();

    TMatrixD identity(6, 6);
    identity.UnitMatrix(); // Create an identity matrix

    // Calculate propagator matrix using intermediate matrices
    // Initialize F as identity matrix
    TMatrixD F(6,6);
    F.UnitMatrix();
    // Compute F1, F2, F3, F4
    Double_t dt_1 = 1e-10;
    TMatrixD IplusFi (6,6);
    IplusFi.Zero();
    TMatrixD Fi(6,6);

    for (int j = 0; j < 4; j++) {
        IplusFi = F* dt_1;
        Fi  = (1e-10)*jacobian*IplusFi;
        F     += Fi *((j==0 || j==3) ?  1.0/6.0 : 1.0/3.0);
        //F.Print();
    }

    // Define matrix to hold time derivatives of state vectors
    TMatrixD  state_dot_matrix(numDimensions,numStates-1);

    // Calculate time derivatives of state vectors
    for (int k = 0; k < numStates; k++) {
    // Extract state vector at time t.
        TMatrixD state_vector(6,1);
        for (int j = 0; j < 6; j++) {
            state_vector(j,0) = state_matrix(j,k);
        }

    // Multiply propagator matrix with state vector to get time derivative of state vector
    TMatrixD state_dot_vector = F * state_vector;
    propagatorx_vs_propagatory->Fill(state_dot_vector(0,0), state_dot_vector(1,0));
    F_projection->Fill(state_dot_vector(0,0),state_dot_vector(1,0),state_dot_vector(2,0));
    }


    // Plot the histogram
    c1->cd(1);
        rx_vs_ry->Draw();
        rx_vs_ry->SetMarkerStyle(21);
        rx_vs_ry->SetMarkerSize(0.3);
    c1->cd(2);
        propagatorx_vs_propagatory->Draw();
        propagatorx_vs_propagatory->SetMarkerStyle(21);
        propagatorx_vs_propagatory->SetMarkerSize(0.3);

    c1->cd(4);
        F_projection->Draw();
        F_projection->SetMarkerStyle(21);
        F_projection->SetMarkerSize(0.3);


    return 0;
}
