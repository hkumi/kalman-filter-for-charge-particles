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

void macro6 () {
	
	// using TMATRICES FROM CERN.
	//DEfine all matrices

	TMatrixD A_1(2,2);
	TMatrixD B_1(2,1);
	TMatrixD X_1(2,1);
        TMatrixD Y_1(2,1); //Observation matrix
	TMatrixD H_1(2,2);
        TMatrixD R_1(2,2);
	TMatrixD C_1(2,2);
	TMatrixD ax_1(2,1);
        TMatrixD P_1(2,2);
	TMatrixD Q_1(2,2);
	TMatrixD I_1(2,2);

	//kalman storage
	TMatrixD x_pred(2,1);
	TMatrixD P_pred(2,2);
	TMatrixD K(2,2);
	TMatrixD Z(2,1);
	

	Int_t n = 5;
	Double_t Q, E_Est[n],E_CE[n],E_Mea,KG[n],q,mn[n],C_E[n];
	Double_t dt,w,ay;
	Double_t dx,dv;
	Double_t ax;
	//initial values.
	dx = 20.0; // meters.
	dv = 5; //m/s.
	Q = 0.0; // Error in the process.
	Double_t x[5] = {4000.0,4260.0,4550.0,4860.0,5110.0};
	Double_t v[5] ={280.0,282.0,285.0,286.0,290.0};
	//x[0] = 4000.0; //  Units in meters.
	//v[0] = 280.0; //    Units in meter/seconds.
	ax= 2.0; 	//in m/sec².
	dt = 1.0;
	C_E[0] = 10; // unit in m.
	Double_t z[10] = {49.95,49.967,50.10,50.106,49.992,49.819,49.933,50.007,50.023,49.99};  // units in m. 
	
	//filling the matrices
	Double_t Matrix_A[4] = {1,dt,0,1};
	A_1.Use(A_1.GetNrows(), A_1.GetNcols(), Matrix_A);

	Double_t Matrix_B[2] = {0.5*TMath::Power(dt,2),dt};
	B_1.Use(B_1.GetNrows(), B_1.GetNcols(), Matrix_B);

	Double_t Matrix_X[2] = {x[0],v[0]};
	X_1.Use(X_1.GetNrows(), X_1.GetNcols(), Matrix_X);

	Double_t Matrix_Y[2] ;
	Y_1.Use(Y_1.GetNrows(), Y_1.GetNcols(), Matrix_Y);

	Double_t Matrix_H[4] = {1,0,0,1} ;
	H_1.Use(H_1.GetNrows(), H_1.GetNcols(), Matrix_H);

	Double_t Matrix_R[4] ={625,0,0,36};
	R_1.Use(R_1.GetNrows(), R_1.GetNcols(), Matrix_R);

	Double_t Matrix_C[4] = {1,0,0,1};
	C_1.Use(C_1.GetNrows(), C_1.GetNcols(), Matrix_C);

	Double_t I[4]  = {1,0,0,1};
	I_1.Use(I_1.GetNrows(), I_1.GetNcols(), I);

	Double_t Matrix_P[4] = {TMath::Power(dx,2),0.0,0.0,TMath::Power(dv,2)} ; // generated errors.
	P_1.Use(P_1.GetNrows(), P_1.GetNcols(), Matrix_P);



	  





	for(Int_t i=0; i<n; i++){

	Double_t Matrix_Y[2] = {x[i],v[i]};
        Y_1.Use(Y_1.GetNrows(), Y_1.GetNcols(), Matrix_Y);

	//start kalman

	x_pred = (A_1 * X_1) +( (B_1*ax) + Q) ;
	P_pred =  (A_1 *TMatrixD(P_1, TMatrixD::kMultTranspose,A_1))  + Q;


	//updates

	K =  TMatrixD(P_pred, TMatrixD::kMultTranspose,H_1) *  ((H_1 * TMatrixD(P_pred, TMatrixD::kMultTranspose,H_1)) + R_1).Invert();
	Z = C_1*Y_1;
	X_1 = x_pred + K *(Z-H_1*x_pred);
	P_1 =(I_1-K*H_1)*P_pred;

	X_1.Print();
        //Z.Print();
	//std::cout<< z[n] << " " << x_pred << std::endl;
	}
	//std::cout<<Q<<endl;
	K.Print();
	Z.Print();
	return 0;


}
