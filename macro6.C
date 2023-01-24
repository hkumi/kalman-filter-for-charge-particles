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

	TCanvas *c1 = new TCanvas("c1", "Radial", 500,600);
// KG= KALMAN GAIN , E_Est = ERROR IN ESTIMATE ,E_Mea = ERROR IN MEASUREMENT,EST_t = CURRENT ESTIMATE,EST_t-1 = PREVIOUS ESTIMATE,MEA MEASUREMENT;
	Int_t n =5; 
	Double_t alpha,beta;
	alpha = 0.8;
	beta = 0.5;

	
	// using TMATRICES FROM CERN.
	//DEfine all matrices

	TMatrixD A_1(2,2);
	TMatrixD B_1(2,1);
	TMatrixD X_1(2,1);
        TMatrixD Y_1(2,1);
	TMatrixD H_1(2,2);
        TMatrixD R_1(2,2);
	TMatrixD C_1(2,2);
	TMatrixD ax_1(2,1);
        TMatrixD P_1(2,2);

	//kalman storage
	TMatrixD x_pred(2,1);
	TMatrixD P_pred(2,2);
	TMatrixD K(2,2);
	//TMatrixD H(2,1);
	//TMatrixD P(2,2);



//	A_1.SetMatrixArray(data.GetArray());


	Double_t Q, E_Est[n],E_CE[n],E_Mea,KG[n],q,mn[n],C_E[n];
	Double_t x[n],dt,w,y[n],ay,v[n];
	Double_t dx,dv;

	//initial values.
	dx = 20.0; // meters.
	dv = 5; //m/s.
	Q = 0.0;
	x[0] = 4000.0; //  Units in meters.
	v[0] = 280.0; //    Units in meter/seconds.
	//ax= 2.0; 	//in m/sec².
	dt = 1.0;
	C_E[0] = 10; // unit in m.
	Double_t z[10] = {49.95,49.967,50.10,50.106,49.992,49.819,49.933,50.007,50.023,49.99};  // units in m. 
	
	//filling the matrices
	Double_t Matrix_A[4] = {1,0,dt,1};
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

	Double_t ax[2]  = {2.0,0.0};
	ax_1.Use(ax_1.GetNrows(), ax_1.GetNcols(), ax);

	Double_t Matrix_P[4] = {TMath::Power(dx,2),0.0,0.0,TMath::Power(dv,2)} ; // generated errors.
	P_1.Use(P_1.GetNrows(), P_1.GetNcols(), Matrix_P);


	for(Int_t i=0; i<2; i++){
	//start kalman
	x_pred = (A_1 * X_1) + (B_1*ax_1) ;
	P_pred = A_1 * P_1 * A_1 + Q;

	//updates
	K =  P_pred * H_1/(H_1 * P_pred * H_1 + R_1);
	X_1 = x_pred + K *(z[n]*x_pred);
	P=(1-K*H_1)*P_pred;
	std::cout<< z[n] << " " << x_pred << std::endl;
	}
	return 0;

      //for printing the TMatrix.



	const Int_t *rIndex = A_1.GetRowIndexArray();
	const Int_t *cIndex = A_1.GetColIndexArray();
	const Double_t *pData = A_1.GetMatrixArray();
	for (Int_t irow = 0; irow < A_1.GetNrows(); irow++) {
		const Int_t sIndex = rIndex[irow];
		const Int_t eIndex = rIndex[irow+1];

		for (Int_t index = sIndex; index < eIndex; index++) {
			const Int_t icol = cIndex[index];
			const Double_t data = pData[index];
			printf("data(%d,%d) = %.4en",irow+A_1.GetRowLwb(),
			icol+A_1.GetColLwb(),data);
		}
	}

    
//	return 0;
/*
	TGraph *gr1 = new TGraph(n,mn,z);
   	TGraph *gr2 = new TGraphErrors(n,mn,C_E);

   // create a multigraph and draw it
   TMultiGraph  *mg  = new TMultiGraph();
   mg->Add(gr1);
   mg->Add(gr2);
   mg->Draw("ALP");
*/
//	TGraph *gr1 = new TGraph(n,mn,C_E);
//	gr1->GetXaxis()->SetTitle("mn");
///	gr1->GetYaxis()->SetTitle("temperature");
//	gr1->Draw("ALP");

}
