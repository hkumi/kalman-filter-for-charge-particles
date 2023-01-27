#include "TCanvas.h"
#include "TGraph.h"
#include <iostream>

#include <vector>
#include <TMath.h>
#include <math.h>
#include <TMatrixD.h>
#include "TVectorD.h"
#include <TMatrixT.h>

#include <cmath>

#include <fstream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <string>



using namespace std;




//Kalman Filter for a radar

void kalman1 () {
	
	// using TMATRICES FROM CERN.
	//DEfine all matrices

	TMatrixD A_1(6,6);
//	TMatrixD B_1(2,1);
	TMatrixD X_1(6,1);
        TMatrixD Y_1(2,1); //Observation matrix
	TMatrixD H_1(2,6);
        TMatrixD R_1(2,2);
	TMatrixD C_1(2,2);
	TMatrixD ax_1(6,1);
        TMatrixD P_1(6,6);
	TMatrixD Q_1(6,6);
	TMatrixD I_1(6,6);

	//kalman storage
	TMatrixD x_pred(6,1);
	TMatrixD P_pred(6,6);
	TMatrixD K(6,2);
	TMatrixD Z(2,1);
	

	Int_t n = 5;
	Double_t dt;
	Double_t ax;
	Double_t v_x[n],a_x[n],v_y[n],a_y[n],p[n],Px[n],Py[n];

	//initial values.
	Double_t x[35] = {-393.66,-375.93,-351.04,-328.96,-299.35,-273.36,-245.89,-222.58,-198.03,-174.17,-146.32,-123.72,-103.47,-78.23,-52.63,-23.34,25.96,49.72,76.94,95.38,119.83,144.01,161.84,180.56,201.42,222.62,239.4,252.51,266.26,271.75,277.4,294.12,301.23,291.8,299.89};
	Double_t y[35] ={300.4,301.78,295.1,305.19,301.06,302.05,300,303.57,296.33,297.65,297.41,299.61,299.6,302.39,295.04,300.09,294.72,298.61,294.64,284.88,272.82,264.93,251.46,241.27,222.98,203.73,184.1,166.12,138.71,119.71,100.41,79.76,50.62,32.99,2.14};
	
	Px[0]=0.0;
	Py[0] = 0.0;

	v_x[0] = 0.0; //  Units in meters.
	a_x[0] = 0.0;
	v_y[0] = 0.0;
	a_y[0] = 0.0;

	//Initial state vector guess.
	p[0] = 500.0; 
	//ax= 2.0; 	//in m/sec².
	dt = 1.0;

	
	//filling the matrices
	Double_t Matrix_A[36] = {1,dt,0.5*TMath::Power(dt,2),0,0,0,0,1,dt,0,0,0,0,0,1,0,0,0,0,0,0,1,dt,0.5*TMath::Power(dt,2),0,0,0,0,1,dt,0,0,0,0,0,1};
	A_1.Use(A_1.GetNrows(), A_1.GetNcols(), Matrix_A);


	//Double_t Matrix_B[2] = {0.5*TMath::Power(dt,2),dt};
	//B_1.Use(B_1.GetNrows(), B_1.GetNcols(), Matrix_B);

	Double_t Matrix_X[6] = {Px[0],v_x[0],a_x[0],Py[0],v_y[0],a_y[0]};
	X_1.Use(X_1.GetNrows(), X_1.GetNcols(), Matrix_X);

	Double_t Matrix_Y[6] ;
	Y_1.Use(Y_1.GetNrows(), Y_1.GetNcols(), Matrix_Y);

	Double_t Matrix_H[12] =  {1,0,0,0,0,0,0,0,0,1,0,0};
	H_1.Use(H_1.GetNrows(), H_1.GetNcols(), Matrix_H);

	Double_t Matrix_R[4] = {9,0,0,9};
	R_1.Use(R_1.GetNrows(), R_1.GetNcols(), Matrix_R);

	Double_t Matrix_C[4] =  {1,0,0,1};
	C_1.Use(C_1.GetNrows(), C_1.GetNcols(), Matrix_C);

	Double_t I[36]  =  {1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1};
	I_1.Use(I_1.GetNrows(), I_1.GetNcols(), I);

	Double_t Matrix_P[36] =  {p[0],0,0,0,0,0,0,p[0],0,0,0,0,0,0,p[0],0,0,0,0,0,0,p[0],0,0,0,0,0,0,p[0],0,0,0,0,0,0,p[0]}; // generated errors.
	P_1.Use(P_1.GetNrows(), P_1.GetNcols(), Matrix_P);

	Double_t Matrix_Q[36] =  {pow(0.2,2)*TMath::Power(dt,4)/4,pow(0.2,2)*TMath::Power(dt,3)/2,pow(0.2,2)*TMath::Power(dt,2)/2,0,0,0,pow(0.2,2)*TMath::Power(dt,3)/2,pow(0.2,2)*TMath::Power(dt,2),pow(0.2,2)*dt,0,0,0,pow(0.2,2)*pow(dt,2)/2,pow(0.2,2)*dt,pow(0.2,2)*1,0,0,0,0,0,0,pow(0.2,2)*TMath::Power(dt,4)/4,pow(0.2,2)*TMath::Power(dt,3)/2,pow(0.2,2)*TMath::Power(dt,2)/2,0,0,0,pow(0.2,2)*pow(dt,3)/2,pow(0.2,2)*pow(dt,2),pow(0.2,2)*dt,0,0,0,pow(0.2,2)*pow(dt,2)/2,pow(0.2,2)*dt,pow(0.2,2)*1}; // generated errors.
        Q_1.Use(Q_1.GetNrows(), Q_1.GetNcols(), Matrix_Q);


	  





	for(Int_t i=0; i<2; i++){

	Double_t Matrix_Y[2] = {x[i],y[i],};
        Y_1.Use(Y_1.GetNrows(), Y_1.GetNcols(), Matrix_Y);

	//start kalman

	x_pred = (A_1 * X_1) ;//+( (B_1*ax) + Q) ;
	P_pred =  (A_1 *TMatrixD(P_1, TMatrixD::kMultTranspose,A_1))  + Q_1;


	//updates

	K =  TMatrixD(P_pred, TMatrixD::kMultTranspose,H_1) *  (H_1 * TMatrixD(P_pred, TMatrixD::kMultTranspose,H_1) + R_1).Invert();
	Z = C_1*Y_1;
	X_1 = x_pred + (K *(Z-(H_1*x_pred)));
	P_1 =(I_1-K*H_1)*P_pred;

	x_pred.Print();
        //Z.Print();
	//std::cout<< z[n] << " " << x_pred << std::endl;
	}


	//std::cout<<Q<<endl;
	P_1.Print();
	//X_1.Print();
	return 0;


}
