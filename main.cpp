///
/// Copyright 2022 Mohanad Youssef (Al-khwarizmi)
///
/// Use of this source code is governed by an GPL-3.0 - style
/// license that can be found in the LICENSE file or at
/// https://opensource.org/licenses/GPL-3.0
///
/// @author Mohanad Youssef <mohanad.magdy.hammad@gmail.com>
/// @file main.cpp
///

#include <iostream>
#include <vector>
#include <TMath.h>
#include <math.h>
#include <string>
#include <TMatrixD.h>
#include <algorithm>
#include <limits>
#include <fstream>
#include "TGraph.h"


#include "kalman_filter/types.h"
#include "kalman_filter/unscented_kalman_filter.h"

using namespace std;

static constexpr size_t DIM_X{ 6 };
static constexpr size_t DIM_V{ 6 };
static constexpr size_t DIM_Z{ 3 };
static constexpr size_t DIM_N{ 3 };

void runExample1();

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
cout<<"this is the values of bet:"<<endl;
cout<<bet<<endl;
cout<<"this is the value for vv:" << endl;
cout<<vv<<endl;
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



kf::Vector<DIM_X> funcF(const kf::Vector<DIM_X> & x, const kf::Vector<DIM_V> & v)
{
    Double_t Ex,Ey,Ez;              //Electric field in V/m.
    Double_t Bx,By,Bz;              // magnetic field in Tesla
    Double_t Energy = GetEnergy(x[3], x[4], x[5]);
    Double_t st = StoppingPower(Energy);

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

    Double_t rr,az,po;
    rr = TMath::Sqrt(TMath::Power(x[3],2)+TMath::Power(x[4],2)+TMath::Power(x[5],2));
    az = TMath::ATan2(x[4] , x[3]);
    po = TMath::ACos(x[5]/rr);

    kf::Vector<DIM_X> y;
    y[0] = x[3];
    y[1] = x[4];
    y[2] = x[5];
    y[3] = q_over_m *(Ex+x[4]*Bz-x[5]*By) - st*TMath::Sin(po)*TMath::Cos(az); //dxdt with energyloss compensation.;
    y[4] = q_over_m *(Ey+x[5]*Bx-x[3]*Bz) - st*TMath::Sin(po)*TMath::Sin(az);
    y[5] = q_over_m *(Ez+x[3]*By-x[4]*Bx) - st*TMath::Cos(po);

    return y;
}




kf::Vector<DIM_Z> funcH(const kf::Vector<DIM_X> & x, const kf::Vector<DIM_N> & n)
{
    kf::Vector<DIM_Z> y;

    y[0] = x[0];
    y[1] = x[1];
    y[2] = x[2];
    return y;
}


int main(int argc, char ** argv)
{
    // example 1
    runExample1();

    return 0;
}


void runExample1()
{


    kf::Vector<DIM_X> x;
    x << 0.0, 0.0, 0.0, 10e6, 10e6, 10e6; //initial values for x,y,z,vx,vy and vz.
    Double_t Energy = GetEnergy(x[3],x[4],x[5]);


/*----------------------------------------*/
    cout<<endl;
    cout<<endl;
    cout<<"+----------------+"<<endl;
    cout<<" Running kalman Filter"<<endl;
    cout<<"+----------------+"<<endl;
    cout<<endl;
/*----------------------------------------*/


    kf::Matrix<DIM_X, DIM_X> P; //in meters
    P << 0.01, 0.0, 0.0, 0.0, 0.0, 0.0,
         0.0, 0.01, 0.0, 0.0, 0.0, 0.0,
         0.0, 0.0, 0.01, 0.0, 0.0, 0.0,
         0.0, 0.0, 0.0, 0.01, 0.0, 0.0,
         0.0, 0.0, 0.0, 0.0, 0.01, 0.0,
         0.0, 0.0, 0.0, 0.0, 0.0, 0.01;

    kf::Matrix<DIM_V, DIM_V> Q;
    Q << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
         0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
         0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
         0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
         0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
         0.0, 0.0, 0.0, 0.0, 0.0, 0.0;


    kf::Matrix<DIM_N, DIM_N> R;
    R << 0.01, 0.0, 0.0,  0.0, 0.01,0.0, 0.0, 0.0, 0.01;

    kf::Vector<DIM_Z> z;
    z << 2.5, 0.05,0.05;

    kf::UnscentedKalmanFilter<DIM_X, DIM_Z, DIM_V, DIM_N> ukf;

    ukf.vecX() = x;
    ukf.matP() = P;

    ukf.setCovarianceQ(Q);
    ukf.setCovarianceR(R);

    ukf.predictUKF(funcF);

    std::cout << "x = \n" << ukf.vecX() << std::endl;
    std::cout << "P = \n" << ukf.matP() << std::endl;



    ukf.correctUKF(funcH, z);

    std::cout << "x = \n" << ukf.vecX() << std::endl;
    std::cout << "P = \n" << ukf.matP() << std::endl;


    std::cout << " End of Kalman: ===========================" << std::endl;

}

