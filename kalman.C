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



void Get_Energy(Double_t M,Double_t IZ,Double_t BRHO,Double_t &E);


std::vector<AtHitCluster> fHitClusterArray; //< Clusterized hits container

std::vector<AtHitCluster> *GetHitClusterArray() {


         return &fHitClusterArray;


 }


std::vector<AtHit> fHitArray;
std::vector<AtHit> *GetHitArray() {


         return &fHitArray;


 }




void  Get_Energy(Double_t M,Double_t IZ,Double_t BRHO,Double_t &E){

  //Energy per nucleon

  Float_t  AM=931.5;
  Float_t X=BRHO/0.1439*IZ/M;
  X=pow(X,2);
  X=2.*AM*X;
  X=X+pow(AM,2);
  E=TMath::Sqrt(X)-AM;


  }



void energyloss(std::string file, TGraph* eLossCurve){

	std::string eLossFileName_;
	eLossFileName_ = file;

        Float_t ener = 0;
        Float_t stp = 0;

	std::ifstream elossfile;


	try {

		elossfile.open(eLossFileName_,ios::in);

		if ((elossfile.peek() == std::ifstream::traits_type::eof())) {

      			std::cout << " Error: Energy loss file not found! Exiting..."<< "\n";
      			std::exit(EXIT_FAILURE);
    		}


	//std::cout << " Processing energy loss data file " << eLossFileName_ << "\n";
	std::string line;

	for (auto i = 0; i < 2; ++i) {
        	std::getline(elossfile,line);
        //	std::cout << line << "\n";
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

// Function that calculates the derivative of the position and velocity of the particle

Double_t  dxdt(double t,double vx, double vy, Double_t vz){

        Double_t Energy = GetEnergy(vx, vy, vz);
        Double_t st = StoppingPower(Energy);
	Double_t Ex,Ey,Ez,f1;            	//Electric field in V/m. 
	Double_t Bx,By,Bz;	      	// magnetic field in Tesla
	Double_t rr,az,po;
	Double_t q = 1.6022*TMath::Power(10,-19); 	//charge of the particle(proton) in (C)
	Double_t m = 1.6726*TMath::Power(10,-27); 	// mass of the particle in kg
	Double_t B=3.0;                 // Applied magnetic field (T).
        Double_t E=TMath::Cos((q*B)/m) * 500;   // Applied Electric field(V/m).
        Bx = 0;
	By = 0;                	// magnetic field in x and y  direction in Tesla.
	Bz = B;          // magnetic field in the z direction in Tesla.
        Ex = 0;
	Ey = 0;			// Electric field in the x and  direction in V/m.
	Ez = -E;                        // Electric field in the z direction. 


        rr = TMath::Sqrt(TMath::Power(vx,2)+TMath::Power(vy,2)+TMath::Power(vz,2));
        az = TMath::ATan2(vy,vx);
        po = TMath::ACos(vz/rr);

	f1 =  q/m * (Ex + vy*Bz-vz*By) - st*TMath::Sin(po)*TMath::Cos(az); //-s*TMath::Sin(po)*TMath::Cos(az) ;                                    //dxdt with energyloss compensation.

        double bro = Bz * rr / TMath::Sin(po)/ 1000.0;
//cout<<Energy<<endl;
//cout<<bro<<endl;
	//std::cout<<Ex <<" "<< vy<<" "<<vz << " " <<By<< " "<< f1<<std::endl; 
        return f1;

        }

double  dydt(double t,double vx, double vy, Double_t vz){
Double_t Energy = GetEnergy(vx, vy, vz);
        Double_t st = StoppingPower(Energy);
	Double_t Ex,Ey,Ez,f2;              //Electric field in V/m. 
        Double_t Bx,By,Bz;              // magnetic field in Tesla
        Double_t q = 1.6022*TMath::Power(10,-19);   //charge of the particle in C
        Double_t B=3.0;                 // Applied magnetic field.
        Double_t rr,az,po;
        Double_t m = 1.6726*TMath::Power(10,-27);  // mass of the particle in kg
	Double_t E = TMath::Cos((q*B)/m) * 500 ;
        Bx = 0;                      // magnetic field in x and y  direction in Tesla.
	By = 0;
        Bz = B;          // magnetic field in the z direction in Tesla.
        Ex = 0;                      // Electric field in the x and  direction in V/m.
	Ey = 0;
        Ez = -E;                        // Electric field in the z direction.

	rr = TMath::Sqrt(TMath::Power(vx,2)+TMath::Power(vy,2)+TMath::Power(vz,2));
        az = TMath::ATan(vy/vx);
        po = TMath::ACos(vz/rr);


       f2 =  q/m * (Ey + vz*Bx - vx*Bz) - st*TMath::Sin(po)*TMath::Sin(az); 

        double bro = Bz * rr / TMath::Sin(po)/ 1000.0;
//cout<<bro<<endl;
//cout<<Energy<<endl;
	//std::cout<<f2<<std::endl;
        return f2;
        }

double  dzdt(double t,double vx, double vy,Double_t vz){
Double_t Energy = GetEnergy(vx, vy, vz);
        Double_t st = StoppingPower(Energy);
        Double_t Ex,Ey,Ez,f3;              //Electric field in V/m. 
        Double_t Bx,By,Bz;              // magnetic field in Tesla
        Double_t q = 1.6022*pow(10,-19);   //charge of the particle in eV
        Double_t B=3.0;                 // Applied magnetic field.
        Double_t rr,az,po;
        Double_t m = 1.6726*TMath::Power(10,-27);  // mass of the particle in kg
	Double_t E = TMath::Cos((q*B)/m) * 500 ; 
        Bx = 0;                      // magnetic field in x and y  direction in Tesla.
	By = 0;
        Bz = B;          // magnetic field in the z direction in Tesla.
        Ex = 0;                      // Electric field in the x and  direction in V/m. 
	Ey = 0;
        Ez = -E;                        // Electric field in the z direction.

	rr = TMath::Sqrt(TMath::Power(vx,2)+TMath::Power(vy,2)+TMath::Power(vz,2));
        az = TMath::ATan(vy/vx);
	//std::cout<<az<<endl;
        po = TMath::ACos(vz/rr);


        f3 =  q/m * (Ez + vx*By - vy*Bx) - st*TMath::Cos(po);
        double bro = Bz * rr / TMath::Sin(po)/ 1000.0;
//cout<<bro<<endl;
//cout<<Energy<<endl<<endl;
	//std::cout <<r<< " " << TMath::Power(vy,2) + TMath::Power(vx,2) + TMath::Power(vz,2)<< " " << TMath::Power(vz,2) <<  std::endl; 
        return f3;
        }

//  Functions for positions. 

Double_t fx(Double_t t, Double_t vx){


       Double_t f1x =  vx; 
        return f1x;


        }

Double_t fy(Double_t t,Double_t vy){


        Double_t f2y = vy;
        return f2y;


        }

Double_t fz(Double_t t,Double_t vz){

        Double_t f3z =  vz; 
        return f3z;

        }

// Define a functions to calculate the Jacobian matrix of the equations.
void Jacobi_matrice(TMatrixD &Jacobi_matrix){

 // Define the size of the matrix
        Int_t rows = 6; // the number of rows. In my case I have a 6*6 matrices for the state vectors
        Int_t cols = 6; // the number of cols.
        Double_t h = 7.49 *TMath::Power(10,-10) ; //in seconds.

        // Define the jacobi matrix to hold state vectors
        Double_t Ex,Ey,Ez,f1;                   //Electric field in V/m. 
        Double_t Bx,By,Bz;              // magnetic field in Tesla
        Double_t rr,az,po;
        Double_t q = 1.6022*TMath::Power(10,-19);       //charge of the particle(proton) in (C)
        Double_t m = 1.6726*TMath::Power(10,-27);       // mass of the particle in kg
        Double_t B=3.0;                 // Applied magnetic field (T).
        Double_t E=TMath::Cos((q*B)/m) * 500;   // Applied Electric field(V/m).
        Bx = 0;
        By = 0;                 // magnetic field in x and y  direction in Tesla.
        Bz = B;          // magnetic field in the z direction in Tesla.
        Ex = 0;
        Ey = 0;                 // Electric field in the x and  direction in V/m.
        Ez = -E;


        Jacobi_matrix.ResizeTo(rows,cols);
        Jacobi_matrix.Zero();

        Jacobi_matrix(0,3) = 1*h;

        Jacobi_matrix[1][4] = 1*h;

        Jacobi_matrix(2,5) = 1*h;


        Jacobi_matrix(3,4) = q/m*(Ex+Bz)*h;
        Jacobi_matrix(3,5) = q/m *(Ex-By)*h;


        Jacobi_matrix(4,3) = q/m * (Ey-Bz)*h;
        Jacobi_matrix(4,5) = q/m*(Ey+Bx)*h;


        Jacobi_matrix(5,3) = q/m * (Ez+By)*h;
        Jacobi_matrix(5,4) = q/m * (Ez-Bx)*h;




}

// Define a functions to calculate the process noise matrix Q .
//this is a 6*6 matrice for my case. 
void generateGaussianNoise(TMatrixD &Q, Double_t mean, Double_t stdDev)
{
    // Define the size of the matrix
    Int_t rows = 6; // the number of rows. In my case I have a 6*6 matrices for the state vectors
    Int_t cols = 6; // the number of cols.

    // Define the matrix to hold process noise. 
    Q.ResizeTo(rows, cols);
    Q.Zero();


}



void generateMeasurementNoise(TMatrixD& R, Double_t xvariance, Double_t yvariance)
{
    // Define the size of the matrix
    Int_t rows = 3;
    Int_t cols = 3;

    // Define the matrix to hold measurement noise
    R.ResizeTo(rows, cols);
    R.Zero();

    R(0,0) = 1e-6;

    R(1,1) = 1e-6;

    R(2,2) = 1e-6;

}



// Define a functions to calculate the Initial Covariance P .
//this is a 6*6 matrice for my case. 
void Ini_P(TMatrixD &P){


 // Define the size of the matrix
        Int_t rows = 6; // the number of rows. In my case I have a 6*6 matrices for the state vectors
        Int_t cols = 6; // the number of cols.

        // Define the noise matrix to hold process noise. 
        P.ResizeTo(rows,cols);
        P.Zero();


        // Define the  matrixes to hold covariance matrix.

        P(0,0) = 1e-4;

        P(1,1) = 1e-4;

        P(2,2) = 1e-4;

        P(3,3) = 1e-7;

        P(4,4) = 1e-7;

        P(5,5) = 1e-7;



}

//Define the Identity matrix
void I_matrix(TMatrixD &I){


 // Define the size of the matrix
        Int_t rows = 6; // the number of rows. In my case I have a 6*6 matrices for the state vectors
        Int_t cols = 6; // the number of cols.

        // Define the noise matrix to hold process noise. 
        I.ResizeTo(rows,cols);
        I.Zero();


        // Define the  matrixes to hold Identity matrix.

        I(0,0) = 1;

        I(1,1) = 1;

        I(2,2) = 1;

        I(3,3) = 1;

        I(4,4) = 1;

        I(5,5) = 1;



}

//Define the initial predicted state for the kalman state vectors. 

void statepred(TMatrixD &Initial_state){


// Define the size of the matrix
        Int_t rows = 6; // the number of rows. 
        Int_t cols = 1; // the number of cols.
        Double_t Am = 1.007; // atomic mass of proton in u. 

 // Define the initial conditions
         Initial_state(0,0) = 0.0;      // initial x position   unit in m.
         Initial_state(1,0) = 0.0;      // initial y position  unit in m
         Initial_state(2,0) = 0.0;      // initial z position  unit in m
         Initial_state(3,0) = 10e6 * Am;     // initial x momentum unit in u m/s
         Initial_state(4,0) = 10e6 * Am;     // initial y momentum  unit in u m/s
         Initial_state(5,0) = 10e6 * Am ;     // initial z momentum unit in u m/s.

}

Double_t GetPhi(Double_t x1, Double_t y1, Double_t x2, Double_t y2) {
         Double_t phi = TMath::ATan2(y2 - y1, -(x2 - x1));
         Double_t phiDeg;
         if (phi < 0){
            phiDeg = 360+ phi*TMath::RadToDeg();
            }  else {
                    phiDeg = phi*TMath::RadToDeg();
                    }
         return phiDeg;
}



// Function to calculate the virtual plane parameters
// Plane equation: ax + by + cz + d = 0
// Inputs:
//    hitPos: Hit position
//    px, py, pz: Track parameters(direction of the track)
// Output:
//     d: Virtual plane parameters
Double_t  GetVirtualPlane( Double_t px, Double_t py, Double_t pz, Double_t x, Double_t y, Double_t z) {


          Double_t a = px;
          Double_t b = py;
          Double_t c = pz;
          Double_t d = a * x + b * y + c * z;
          return d;

}


void GetPOCA(Double_t x1, Double_t y1, Double_t z1, Double_t px, Double_t py, Double_t pz, Double_t& xi, Double_t& yi, Double_t& zi) {
     Double_t d = GetVirtualPlane(px,py,pz,x1,y1,z1);
     Double_t dx = px;
     Double_t dy = py;
     Double_t dz = pz;
     Double_t t = d - (px*x1 + py*y1 + pz*z1)  / (px*dx + py*dy + pz*dz);                       // calculate the intersection
     xi = x1 + t * dx;
     yi = y1 + t * dy;
     zi = z1 + t * dz;

     std::cout << "t value: " << t << endl;
     std::cout << "Intersection point: (" << xi << ", " << yi << ", " << zi << ")" << std::endl;
     std::cout <<endl;
}


void GetMeasurementMatrix(TMatrixD& H_k, Double_t px, Double_t py, Double_t pz, Double_t x1, Double_t y1, Double_t z1) {
    // Define the size of the matrix
    Int_t rows = 3;
    Int_t cols = 6;

    // Define the matrix to hold the measurement matrix
    H_k.ResizeTo(rows, cols);
    H_k.Zero();

    Double_t d = GetVirtualPlane(px,py,pz,x1,y1,z1);
    Double_t dx = px;
    Double_t dy = py;
    Double_t dz = pz;
    Double_t t = d - (px*x1 + py*y1 + pz*z1)  / (px*dx + py*dy + pz*dz);                       // calculate the intersection


    // Set the values of the measurement matrix
    H_k(0, 0) = 1;
    H_k(1, 1) = 1;
    H_k(2, 2) = 1;

}

void ProjectMeasurement(TMatrixD& h_k, TMatrixD state_k, Double_t px, Double_t py, Double_t pz, Double_t x1, Double_t y1, Double_t z1) {
    // Define the measurement matrix
    TMatrixD H_k(3, 6);

    // Get the measurement matrix
    GetMeasurementMatrix(H_k, px, py, pz, x1, y1, z1);

    // Project the measurement
    h_k.ResizeTo(3, 1);
    h_k = H_k * state_k;
}





void kalman(){

         TH2F *rx_vs_ry = new TH2F("rx_vs_ry", "rungex_vs_rungey", 720, 0, -3, 1000, 0, 2.0);
         TH2F *kx_vs_ky = new TH2F("kx_vs_ky", "kalmanx_vs_kalmany", 720, 0, -3, 1000, 0, 2.0);

         TH2F *vx_vs_vy = new TH2F("vx_vs_rvy", "rungevx_vs_rungevy", 720, 0, -3, 1000, 0, 2.0);
         TH2F *Intx_vs_Inty = new TH2F("Intx_vs_Inty", "Intersectionx_vs_Intersectiony", 720, 0, -3, 1000, 0, 2.0);

         TH2F *propagatorx_vs_propagatory = new TH2F("propagatorx_vs_propagatory", "propagatorx_vs_propagatory", 720, 0, -3, 1000, 0, 2.0);
  //       TH2F *kx_vs_ky = new TH2F("kx_vs_ky", "kalmanx_vs_kalmany", 720, 0, -3, 100, 0, 2.0);

        TH3F *R_projection = new TH3F("R_projection", "runge_projection", 720, 5.0, -3, 100, 0, 10.0,100, -5, 10.0);
        TH3F *F_projection = new TH3F("F_projection", "propagator_projection", 720, 5.0, -3, 100, 0, 10.0,100, -5, 10.0);
        TH3F *Intersection = new TH3F("Intersection", "Intersection", 720, 5.0, -3, 100, 0, 10.0,100, -5, 10.0);



         TCanvas *c1 = new TCanvas();
         c1->Divide(2, 2);
         c1->Draw();


 	//TCanvas *c1 = new TCanvas("c1","Particle in a magnetic field",500,600);
        //c1->Divide(2,2);

        const Int_t n = 1000;

	Double_t t[n],vx[n],vy[n],vz[n];
	Double_t x[n],y[n],z[n],z0,zf;

	Double_t k1x[n],k2x[n],k3x[n],k4x[n];
        Double_t k1y[n],k2y[n],k3y[n],k4y[n];
	Double_t k1z[n],k2z[n],k3z[n],k4z[n];

	Double_t k1vx[n],k2vx[n],k3vx[n],k4vx[n];
        Double_t k1vy[n],k2vy[n],k3vy[n],k4vy[n];
        Double_t k1vz[n],k2vz[n],k3vz[n],k4vz[n];

       // Define the initial conditions
         x[0] = 0.0;      // initial x position   unit in m.
         y[0] = 0.0;      // initial y position  unit in m
         z[0] = 0.0;      // initial z position  unit in m
         vx[0] = 10e6;     // initial x velocity unit in m/s
         vy[0] = 10e6;     // initial y velocity unit in m/s
         vz[0] = 10e6;     // initial z velocity unit in m/s.
        

        // Define the time step 

        Double_t h = TMath::Power(10,-10) ; //in seconds.
	//double tf=8.18*TMath::Power(10,-7);
	//double t0=0.0;


	//Graph to evaluate Energy loss
	//Double_t  gasMediumDensity_ = 8.375e-5;   //g/cm3
	Double_t p[n],s[n],vv[n];
	TGraph* eLossCurve = new TGraph();

       energyloss("stpHydrogen.txt", eLossCurve);

	//eLossCurve->Draw();  // Plotting energy vs stopping power

	//std::cout<<x[0] << " " << y[0] << " " << z[0] << std::endl;
	//Start Rk4. 


	//for(Int_t i=0;i <n-1;  i++){
        Double_t Energy = GetEnergy(vx[0],vy[0],vz[0]);
        Int_t i = 0;

/*----------------------------------------*/
cout<<endl;
cout<<endl;
cout<<"+----------------+"<<endl;
cout<<" Running kalman.C"<<endl;
cout<<"+----------------+"<<endl;

//cout<<"Initial Energy value: "<< Energy << " MeV" <<endl;
//cout<<endl;
//cout<<endl;
/*----------------------------------------*/



        while (Energy>0.1){
	   k1x[i] = fx(t[i],vx[i]);
           k1y[i] = fy(t[i],vy[i]);
           k1z[i] = fz(t[i],vz[i]);

	   k1vx[i] = dxdt(t[i],vx[i],vy[i],vz[i]);
           k1vy[i] = dydt(t[i],vx[i],vy[i],vz[i]);
           k1vz[i] = dzdt(t[i],vx[i],vy[i],vz[i]);

           k2x[i] = fx(t[i]+h*0.5,vx[i]+0.5*h*k1x[i]);
           k2y[i] = fy(t[i]+h*0.5,vy[i]+0.5*h*k1y[i]);
           k2z[i] = fz(t[i]+h*0.5,vz[i]+0.5*h*k1z[i]);

           k2vx[i] = dxdt(t[i]+h*0.5,vx[i]+k1vx[i]*0.5*h,vy[i]+k1vy[i]*h*0.5,vz[i]+k1vz[i]*h*0.5);
           k2vy[i] = dydt(t[i]+h*0.5,vx[i]+k1vx[i]*h*0.5, vy[i]+k1vy[i]*h*0.5,vz[i]+k1vz[i]*h*0.5);
           k2vz[i] = dzdt(t[i]+h*0.5,vx[i]+k1vx[i]*h*0.5, vy[i]+k1vy[i]*h*0.5, vz[i]+k1vz[i]*h*0.5);


	   k3x[i] = fx(t[i]+h*0.5,vx[i]+0.5*h*k2x[i]);
           k3y[i] = fy(t[i]+h*0.5,vy[i]+0.5*h*k2y[i]);
           k3z[i] = fz(t[i]+h*0.5,vz[i]+0.5*h*k2z[i]);

           k3vx[i] = dxdt(t[i]+h*0.5,vx[i]+k2vx[i]*0.5*h,vy[i]+k2vy[i]*h*0.5,vz[i] + k2vz[i]*h*0.5);
           k3vy[i] = dydt(t[i]+h*0.5,vx[i]+k2vx[i]*h*0.5, vy[i]+k2vy[i]*h*0.5,vz[i]+k2vz[i]*h*0.5);
           k3vz[i] = dzdt(t[i]+h*0.5,vx[i]+k2vx[i]*h*0.5, vy[i]+k2vy[i]*h*0.5, vz[i]+k2vz[i]*h*0.5);

           k4x[i] = fx(t[i]+h,vx[i]+k3x[i]*h);
           k4y[i] = fy(t[i]+h,vy[i]+k3y[i]*h);
           k4z[i] = fz(t[i]+h,vz[i]+k3z[i]*h);

           k4vx[i] = dxdt(t[i]+h,vx[i]+k3vx[i]*h,vy[i]+k3vy[i]*h,vz[i]+k3vz[i]*h);
           k4vy[i] = dydt(t[i]+h,vx[i]+k3vx[i]*h,vy[i]+k3vy[i]*h,vz[i]+k3vz[i]*h);
           k4vz[i] = dzdt(t[i]+h,vx[i]+k3vx[i]*h,vy[i]+k3vy[i]*h,vz[i]+k3vz[i]*h);


          vx[i+1] = vx[i] + h/6 *( k1vx[i] + 2*k2vx[i] + 2*k3vx[i] + k4vx[i]);
          vy[i+1]  = vy[i] + h/6 *( k1vy[i] + 2*k2vy[i] + 2*k3vy[i] + k4vy[i]);
	  vz[i+1]  = vz[i] + h/6 *( k1vz[i] + 2*k2vz[i] + 2*k3vz[i] + k4vz[i]);

	  x[i+1] = x[i] + h/6 *( k1x[i] + 2*k2x[i] + 2*k3x[i] + k4x[i]);
          y[i+1]  = y[i] + h/6 *( k1y[i] + 2*k2y[i] + 2*k3y[i] + k4y[i]);
          z[i+1]  = z[i] + h/6 *( k1z[i] + 2*k2z[i] + 2*k3z[i] + k4z[i]);

          t[i+1]= t[i] + h;

	 // std::cout<<vx[i]<< " "<<vy[i]<<" " << vz[i] << " " << h << endl;

          vx[i] = vx[i+1];
	  vz[i] = vz[i+1];
          vy[i] = vy[i+1];

          x[i] = x[i+1];
          z[i] = z[i+1];
          y[i] = y[i+1];

          t[i] = t[i+1];

        // std::cout<<x[i+1]<<" "<<y[i+1]<<" "<<vz[i+1]<<endl;

        rx_vs_ry->Fill(x[i],y[i]);
        R_projection->Fill(x[i],y[i],z[i]);

        Energy = GetEnergy(vx[i],vy[i],vz[i]);
        i+=1;
        }

        // Define the size of the matrix
        Int_t rows = 6; // the number of state variables
        Int_t cols = i; // the number of steptime.

	// Define matrix to hold state vectors
        TMatrixD  state_matrix(rows,cols);

        // Fill in matrix with state vectors
        for (int j = 0; j < i; j++) {
            state_matrix(0,j) = x[j];
            state_matrix(1,j) = y[j];
            state_matrix(2,j) = z[j];
            state_matrix(3,j) = vx[j];
            state_matrix(4,j) = vy[j];
            state_matrix(5,j) = vz[j];
        }


	TMatrixD Jacobi_matrix(6,6);
	Jacobi_matrice(Jacobi_matrix);
        //Jacobi_matrix.Print();


       //Identity matrix.
       TMatrixD I(6,6);
       I_matrix(I);
       //I.Print();

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
             Fi  = (1e-10)*Jacobi_matrix*IplusFi;
             F     += Fi *((j==0 || j==3) ?  1.0/6.0 : 1.0/3.0);
          //F.Print();
        }

        // Define matrix to hold time derivatives of state vectors
        TMatrixD  state_dot_matrix(rows,cols-1);

     // Calculate time derivatives of state vectors
        for (int k = 0; k < cols; k++) {
    // Extract state vector at time t. 
            TMatrixD state_vector(6,1);
            for (int j = 0; j < 6; j++) {
                state_vector(j,0) = state_matrix(j,k);
            }

    // Multiply propagator matrix with state vector to get time derivative of state vector
            TMatrixD state_dot_vector = F * state_vector;
            propagatorx_vs_propagatory->Fill(state_dot_vector(0,0), state_dot_vector(1,0));
            F_projection->Fill(state_dot_vector(0,0),state_dot_vector(1,0),state_dot_vector(2,0));
          // state_vector.Print();
        }
/*

c1->cd(1);
        rx_vs_ry->Draw();
        rx_vs_ry->SetMarkerStyle(21);
        rx_vs_ry->SetMarkerSize(0.3);
        c1->cd(2);
        propagatorx_vs_propagatory->Draw();
        propagatorx_vs_propagatory->SetMarkerStyle(21);
        propagatorx_vs_propagatory->SetMarkerSize(0.3);
        c1->cd(3);
        R_projection->Draw();
        R_projection->SetMarkerStyle(21);
        R_projection->SetMarkerSize(0.3);
        c1->cd(4);
        F_projection->Draw();
        F_projection->SetMarkerStyle(21);
        F_projection->SetMarkerSize(0.3);
*/


// Needs rethink!!!!

        TCanvas *c2 = new TCanvas();
        c2->Divide(2, 2);
        c2->Draw();
        TH2F *angle_vs_energy = new TH2F("angle_vs_energy", "angle_vs_energy", 720, 0, 179, 1000, 0, 100.0);
        TH2F *Z_vs_Y = new TH2F("Z_vs_Y", "Z_vs_Y", 720, 0, -3, 1000, 0, 2.0);
        TH2F *energy_vs_Zorb = new TH2F("energy_vs_Zorb", "energy_vs_Zorb", 720, 0, -5, 100, 0, 5.0);
        TH1F *phi_pattern = new TH1F("phi_pattern", "phi_pattern", 1440, -50, 400);
        TH2F *angle_vs_energy_pattern =new TH2F("angle_vs_energy_pattern", "angle_vs_energy_pattern", 720, 0, 179, 1000, 0, 100.0);


 //to read the data file.

/*----------------------------------------*/
cout<<endl;
cout<<endl;
cout<<"+----------------+"<<endl;
cout<<" Testing the kalman Filter "<<endl;
cout<<"+----------------+"<<endl;

//cout<<"Initial Energy value: "<< Energy << " MeV" <<endl;
//cout<<endl;
//cout<<endl;
/*----------------------------------------*/



        FairRunAna *run = new FairRunAna();
        std::vector<TString> files = {"data/output_digi"};
        Int_t nEvents = 1000;
        Int_t first_Event = 0;
        Int_t last_Event = nEvents;
       // TMatrixD plane1 = ();
        //read_file(files, nEvents);


        for (auto iFile = 0; iFile < files.size(); ++iFile) {

            TString mcFileNameHead = files[iFile];
            TString mcFileNameTail = ".root";
            TString mcFileName = mcFileNameHead + mcFileNameTail;
            std::cout << " Analysis of simulation file  " << mcFileName << endl;

            TFile *file = new TFile(mcFileName.Data(), "READ");
            TTree *tree = (TTree *)file->Get("cbmsim");

            tree = (TTree *)file->Get("cbmsim");
            TTreeReader Reader1("cbmsim", file);
            TTreeReaderValue<TClonesArray> eventArray(Reader1, "AtPatternEvent");

            ROOT::Math::XYZPoint iniPos;
            ROOT::Math::XYZPoint secPos;
            //Loop over  the events in the simulation. 
            for (Int_t i = first_Event; i < last_Event; i++) {
               // std::cout << " Event Number : " << i << " has " << "\n";
                Reader1.Next();

                AtPatternEvent *patternEvent = (AtPatternEvent *)eventArray->At(0);

                if (patternEvent) {
                    std::vector<AtTrack> &patternTrackCand = patternEvent->GetTrackCand();
                 //   std::cout <<  patternTrackCand.size() << " Number of pattern tracks "<< "\n";

                    for (auto track : patternTrackCand) {

                     //   std::cout << " with === Track " << track.GetTrackID() << " also having: "<< track.GetHitClusterArray()->size() << " clusters "  << "\n";
                       // std::cout <<endl;

                        if ( track.GetHitClusterArray()->size() < 5) {
                           //std::cout << " Track is noise or has less than 5 clusters! "  << "\n";
                           //std::cout<<endl;
                           continue;
                        }

                        Double_t theta = track.GetGeoTheta();
                        Double_t rad   = track.GetGeoRadius();
                        Double_t B_f = 3.0;                         //in Tesla.
                        Double_t brho = B_f * rad / TMath::Sin(theta) / 1000.0;     //in Tm.
                        Double_t ener = 0;
                        Double_t Am = 1.007; // atomic mass of proton in u. 
                        Double_t Z = 1.0;
                        Get_Energy(Am, Z, brho, ener);
                        angle_vs_energy_pattern->Fill(theta * TMath::RadToDeg(), ener);
                        Double_t p = brho * Z * 2.99792458*100.0;      //Magnitude of the momentum In MeV/c.
                       // std::cout << "This particle has:";
                       // std:: cout << "Ener:"<<ener<<" Momentum:" << p << " and"   <<  " Magnetic Rigidity: " << brho <<  endl;
                       // std::cout<<endl;
                        auto hitClusterArray = track.GetHitClusterArray();

                        AtHitCluster iniCluster;
                        AtHitCluster SecCluster;

                        TMatrixD x_pred(6,1);
                        TMatrixD P_pred(6,6);
                        TMatrixD K(6,3);
                        TMatrixD Y_1(3,1); //Observation matrix
                        TMatrixD C_1(3,3);
                        TMatrixD Z_1(3,1);


                        TMatrixD Q(6,6);
                        generateGaussianNoise(Q, 0.0, 1.0); // Generates a 6x6 matrix of Gaussian noise with mean 0 and standard deviation 1
                         Q.Print();

                        TMatrixD P(6,6);
                        Ini_P(P);

                        TMatrixD X_1(6,1);
        		statepred(X_1);

                        // Error in Measurement.
                        TMatrixD R_1(3,3);
                        generateMeasurementNoise(R_1,0.0,1.0);

                        Double_t Matrix_C[9] =  {1,0,0,0,1,0,0,0,1};
                        C_1.Use(C_1.GetNrows(), C_1.GetNcols(), Matrix_C);


                        Double_t xaxis[hitClusterArray->size()],yaxis[hitClusterArray->size()],zaxis[hitClusterArray->size()];

                        std::vector<int> eventNumbers = {395};
                        // Iterate over event numbers and access the corresponding events
                        if (std::find(eventNumbers.begin(), eventNumbers.end(), i) != eventNumbers.end()) {

                           std::cout<< "Processing event " << i  << "with " << track.GetHitClusterArray()->size() << " clusters" << endl;

                           for(auto iclus = 0; iclus < hitClusterArray->size()-1; ++iclus){
                            //  std::cout << "Loop iteration: " << iclus << std::endl;
                              //std::cout << "Vector size: " << hitClusterArray->size() << std::endl;
                              auto Cluster1 = hitClusterArray->at(iclus+1);
                              auto inipos = Cluster1.GetPosition();
                              Double_t x1 = inipos.X();
                              Double_t y1 = inipos.Y();
                              Double_t z1 = inipos.Z();
                            //  std::cout << x1<<","<<y1<< ","<< z1<<endl;
                              Double_t distance = TMath::Sqrt(x1 * x1 + y1 * y1);
                              Z_vs_Y->Fill(x1,y1);

                              // Select another cluster to compare with Cluster1

                             auto Cluster2 = hitClusterArray->at(iclus);
                             auto inipos2 = Cluster2.GetPosition();
                             Double_t x2 = inipos2.X();
                             Double_t y2 = inipos2.Y();
                             Double_t z2 = inipos2.Z();

                             // Calculate phi angle between the two points
                             Double_t phiDeg = GetPhi(x1,y1,x2,y2);

                             phi_pattern->Fill(phiDeg);

                             auto px = p * cos(phiDeg) * sin(theta);
                             auto py = p * sin(phiDeg) * sin(theta);
                             auto pz = p * cos(theta);
                             //std::cout << px << " ," << py << " ," << pz << std::endl;
                             // Plane equation: ax + by + cz = d
                             Double_t plane = GetVirtualPlane(px,py,pz,x1,y1,z1);
                             Double_t xi,yi,zi = 0.0;
                             GetPOCA(x1,y1,z1,px,py,pz,xi,yi,zi);
                             //std::cout << "xi:" << xi << "yi:" << yi <<  "zi:" << zi <<  endl; 
                             //Intx_vs_Inty->Fill(xi,yi);



                            //Perform Kalman here.
                            //Initial state predictions for protons.

                             TMatrixD H_1(3,6);
                             // Get the measurement matrix
                             GetMeasurementMatrix(H_1, px, py, pz, x1, y1, z1);



                             x_pred = (F * X_1) ;
                             P_pred =  (F *TMatrixD(P, TMatrixD::kMultTranspose,F))  + Q;
                             //updates
                             K =  TMatrixD(P_pred, TMatrixD::kMultTranspose,H_1) *  (H_1 * TMatrixD(P_pred, TMatrixD::kMultTranspose,H_1) + R_1).Invert();

                             xaxis[iclus] = xi;
                             yaxis[iclus] = yi;
                             zaxis[iclus] = zi;
                             //std::cout << "x:" <<xaxis[iclus] <<"y:" << yaxis[iclus] << "z:"  << zaxis[iclus]<< endl;
                             Double_t Matrix_Y[3] = {xaxis[iclus],yaxis[iclus],zaxis[iclus]};
                             Y_1.Use(Y_1.GetNrows(), Y_1.GetNcols(), Matrix_Y);

                             Z_1 = C_1* Y_1;
                             X_1 = x_pred + (K *(Z_1-(H_1*x_pred)));
                             P =(I-K*H_1)*P_pred;

                         //    std::cout << "the predicted state:" <<std::endl;
                            // x_pred.Print(); 
                           //  std::cout << "the predicted covariance:" <<std::endl; 
                           //  P_pred.Print();
                             kx_vs_ky->Fill(X_1(0,0), X_1(1,0)); 


                            //Y_1.Print();
                           }

                        }

                    }

                }
            }


        }

c2->cd(1);
        angle_vs_energy_pattern->SetMarkerStyle(20);
        angle_vs_energy_pattern->SetLineWidth(1);
        angle_vs_energy_pattern->Draw();
c2->cd(2);
        Z_vs_Y->SetMarkerStyle(20);
        Z_vs_Y->SetLineWidth(1);
        Z_vs_Y->Draw();

c2->cd(3);
        phi_pattern->SetMarkerStyle(20);
        phi_pattern->SetLineWidth(1);
        phi_pattern->Draw();

c2->cd(4);
        kx_vs_ky->SetMarkerStyle(21);
        kx_vs_ky->SetLineWidth(1);
        kx_vs_ky->Draw();


/*

// We basically just transfer the cordinates into array..
        Int_t n2 = tcor.size();
        Double_t xaxis[n],yaxis[n],zaxis[n];
  


        //kalman storage
        TMatrixD x_pred(6,1);
        TMatrixD P_pred(6,6);
        TMatrixD K(6,3);
  //      TMatrixD H(6,6);


       // Initial Values 
        // Define the initial conditions
         x[0] = 0.0;      // initial x position
         y[0] = 0.0;      // initial y position
         z[0] = 0.0;      // initial z position
         vx[0] = 1.0;     // initial x velocity
         vy[0] = 1.0;     // initial y velocity
         vz[0] = 1.0;     // initial z velocity


        TMatrixD X_1(6,1);
        Double_t Matrix_X[6] = {x[0],y[0],z[0],vx[0],vy[0],vz[0]};
        X_1.Use(6, 1, Matrix_X);
   
        TMatrixD Q(6,6);
        Process_noise(Q);

        TMatrixD P(6,6);
        Ini_P(P);


        //Measurement matrice. 
        //Values of the measurements. Should be modified with real data. this is just my example. 
        //initial values.
      //  Double_t x1[11] = {0,0.03,0.004,-0.02,-04.93,-05.04,-08.96,-099.35,-073.36,-045.89,-022.58};
        //Double_t y1[11] ={0,099.6,02.39,05.04,01.09,-02.72,-08.61,-02.64,-01.88,-1.27,02.98};
        //Double_t z1[11] ={0,49.6,32.39,25.04,30.09,-24.72,-28.61,-24.64,-28.88,24.27,22.98};

       //Beacuse only x, yand z are observed.
        TMatrixD Z(3,1);
        TMatrixD Y_1(3,1); //Observation matrix
        TMatrixD H_1(3,6);
        TMatrixD R_1(3,3);
        TMatrixD C_1(3,3);

       //fill the observation matrice.
         Double_t Matrix_Y[3] ;
         Y_1.Use(Y_1.GetNrows(), Y_1.GetNcols(), Matrix_Y);


         Double_t Matrix_H[18] =  {1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0};
         H_1.Use(H_1.GetNrows(), H_1.GetNcols(), Matrix_H);

        // Error in Measurement.
         Double_t Matrix_R[9] = {1e-4,0,0,0,1e-4,0,0,0,1e-4};
         R_1.Use(R_1.GetNrows(), R_1.GetNcols(), Matrix_R);

         Double_t Matrix_C[9] =  {1,0,0,0,1,0,0,0,1};
         C_1.Use(C_1.GetNrows(), C_1.GetNcols(), Matrix_C);


        //start kalman
        for (Int_t i=0;i<n2;++i){
            xaxis[i] = xcor.at(i);
            zaxis[i] = zcor.at(i);
            yaxis[i]  = ycor.at(i);

            Double_t Matrix_Y[3] = {xaxis[i],yaxis[i],zaxis[i]};
            Y_1.Use(Y_1.GetNrows(), Y_1.GetNcols(), Matrix_Y);

            x_pred = (F * X_1) ; 
            P_pred =  (F *TMatrixD(P, TMatrixD::kMultTranspose,F))  + Q;
        

      //updates



            K =  TMatrixD(P_pred, TMatrixD::kMultTranspose,H_1) *  (H_1 * TMatrixD(P_pred, TMatrixD::kMultTranspose,H_1) + R_1).Invert();

            Z = C_1* Y_1;
            X_1 = x_pred + (K *(Z-(H_1*x_pred)));
            P =(I-K*H_1)*P_pred;

            std::cout << "the predicted state:" <<std::endl;
            x_pred.Print(); 
            std::cout << "the predicted covariance:" <<std::endl; 
            P_pred.Print();
             kx_vs_ky->Fill(x_pred(0,0), x_pred(1,0)); 
            //Z.Print();
           // std::cout<< "measurement" << Z << "EStimate" << x_pred <<  std::endl;

        }


        c2->cd(1);
        kx_vs_ky->Draw();
        kx_vs_ky->SetMarkerStyle(21);
        kx_vs_ky->SetLineWidth(1);
        kx_vs_ky->Draw();
*/
  //      file.close();
    //    writefile.close();
        return(0);


}


