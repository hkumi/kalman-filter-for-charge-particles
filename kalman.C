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


using namespace std;


// RK4 for alpha particle moving in electric and magnetic field.
//declaration of the functions.


//Function for the stopping power .
void energyloss(std::string file, TGraph* eLossCurve){

	std::string eLossFileName_;
	eLossFileName_ = file;

	Float_t ener = 0;
	TString enerUnit = "";
	Float_t dEdx_elec = 0;
	Float_t dEdx_nucl = 0;
	Float_t range = 0;
	TString rangeUnit = "";
	Float_t lonStra = 0;
	TString lonStraUnit = "";
	Float_t latStra = 0;
	TString latStraUnit = "";

	std::ifstream elossfile;


	try {

		elossfile.open(eLossFileName_,ios::in);

		if ((elossfile.peek() == std::ifstream::traits_type::eof())) {

      			std::cout << " Error: Energy loss file not found! Exiting..."<< "\n";
      			std::exit(EXIT_FAILURE);
    		}


	//std::cout << " Processing energy loss data file " << eLossFileName_ << "\n";
	std::string line;

	for (auto i = 0; i < 3; ++i) {
        	std::getline(elossfile,line);
        //	std::cout << line << "\n";
    	}


	while (std::getline(elossfile, line)) {
		std::istringstream data(line);
		data >> ener >> enerUnit >> dEdx_elec >> dEdx_nucl >> range >> rangeUnit >> lonStra >> lonStraUnit >> latStra >> latStraUnit ;
      		if(enerUnit.Contains("keV"))
        	ener/=1000.0;

 //     std::cout<<ener<<" "<<enerUnit.Data()<<" "<<dEdx_elec<<" "<<dEdx_nucl<<" "<<range<<" "<<rangeUnit.Data()<<" "<<lonStra<<" "<<lonStraUnit.Data()<<" "<<latStra<<" "<<latStraUnit<<"\n";
      	//	if (ener > maxKinEnergy_)
        //	maxKinEnergy_ = ener;
      		eLossCurve->SetPoint(eLossCurve->GetN(),ener,dEdx_elec+dEdx_nucl);
      	//	elosscurve->SetPoint(elosscurve->GetN(),ener,dEdx_elec+dEdx_nucl);
      //std::cout<<eLossCurve_->GetN()<<" "<<elosscurve->GetN()<<"\n";//NB: eLossCurve_ only valid if only one AtFitter::AtGenfit object is created.
         //	if(elossfile->eof()) break;
    	}


  } catch (...) {
  }

}


// Function that calculates the derivative of the position and velocity of the particle

Double_t  dxdt(double t,double vx, double vy, Double_t vz){

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


	f1 =  q/m * (Ex + vy*Bz-vz*By); //-s*TMath::Sin(po)*TMath::Cos(az) ;                                    //dxdt with energyloss compensation.

	//std::cout<<Ex <<" "<< vy<<" "<<vz << " " <<By<< " "<< f1<<std::endl; 
        return f1;

        }

double  dydt(double t,double vx, double vy, Double_t vz){
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
 

        f2 =  q/m * (Ey + vz*Bx - vx*Bz);// - s*TMath::Sin(po)*TMath::Sin(az); 

	//std::cout<<f2<<std::endl;
        return f2;
        }

double  dzdt(double t,double vx, double vy,Double_t vz){
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
 

        f3 =  q/m * (Ez + vx*By - vy*Bx);// - s*TMath::Cos(po);
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
void Jacobi_matrice(Double_t Jacobi_matrix[][6]){

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


        Jacobi_matrix[0][0] = 0;
        Jacobi_matrix[0][1] = 0;
        Jacobi_matrix[0][2] = 0;
        Jacobi_matrix[0][3] = 1*h;
        Jacobi_matrix[0][4] = 0;
        Jacobi_matrix[0][5] = 0;

        Jacobi_matrix[1][0] = 0;
        Jacobi_matrix[1][1] = 0;
        Jacobi_matrix[1][2] = 0;
        Jacobi_matrix[1][3] = 0;
        Jacobi_matrix[1][4] = 1*h;
        Jacobi_matrix[1][5] = 0;

	Jacobi_matrix[2][0] = 0;
        Jacobi_matrix[2][1] = 0;
        Jacobi_matrix[2][2] = 0;
        Jacobi_matrix[2][3] = 0;
        Jacobi_matrix[2][4] = 0;
        Jacobi_matrix[2][5] = 1*h;

        Jacobi_matrix[3][0] = 0;
        Jacobi_matrix[3][1] = 0;
        Jacobi_matrix[3][2] = 0;
        Jacobi_matrix[3][3] = 0;
        Jacobi_matrix[3][4] = q/m*(Ex+Bz)*h;
        Jacobi_matrix[3][5] = q/m *(Ex-By)*h;

        Jacobi_matrix[4][0] = 0;
        Jacobi_matrix[4][1] = 0;
        Jacobi_matrix[4][2] = 0;
        Jacobi_matrix[4][3] = q/m * (Ey-Bz)*h;
        Jacobi_matrix[4][4] = 0;
        Jacobi_matrix[4][5] = q/m*(Ey+Bx)*h;
 

	Jacobi_matrix[5][0] = 0;
        Jacobi_matrix[5][1] = 0;
        Jacobi_matrix[5][2] = 0;
        Jacobi_matrix[5][3] = q/m * (Ez+By)*h;
        Jacobi_matrix[5][4] = q/m * (Ez-Bx)*h;
        Jacobi_matrix[5][5] = 0;



}



void kalman(){

 	TCanvas *c1 = new TCanvas("c1","Particle in a magnetic field",500,600);
        c1->Divide(2,2);

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
         x[0] = 0.0;      // initial x position
         y[0] = 0.0;      // initial y position
         z[0] = 0.0;      // initial z position
         vx[0] = 1.0;     // initial x velocity
         vy[0] = 1.0;     // initial y velocity
         vz[0] = 1.0;     // initial z velocity
        

        // Define the time step and total time
        double dt = 0.01;     // time step size
        double T = 10.0;      // total time
        Double_t h = 7.49 *TMath::Power(10,-10) ; //in seconds.
	//double tf=8.18*TMath::Power(10,-7);
	//double t0=0.0;


	//Graph to evaluate Energy loss
	Double_t  gasMediumDensity_ = 0.0153236;   //g/cm3
	Double_t p[n],s[n],vv[n],bet[n],gamma1[n],Energy[n];
	TGraph* eLossCurve = new TGraph();

       energyloss("van.txt", eLossCurve);

	//eLossCurve->Draw();

	//std::cout<<x[0] << " " << y[0] << " " << z[0] << std::endl;
	//Start Rk4. 


	for(Int_t i=0;i <n-1;  i++){
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

           k3vx[i] = dxdt(t[i]+h*0.5,vx[i]+k2vx[i]*0.5*h,vy[i]+k2vy[i]*h*0.5,k2vz[i]*vz[i]*h*0.5);
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

        // std::cout<<x[i+1]<<" "<<y[i+1]<<" "<<z[i+1]<<endl;


        }
        // Define the size of the matrix
        Int_t rows = 6; // the number of state variables
        Int_t cols = n; // the number of steptime.

	// Define matrix to hold state vectors
        Double_t  state_matrix[rows][cols];

        // Fill in matrix with state vectors
        for (int i = 0; i < n-1; i++) {
            state_matrix[0][i] = x[i];
            state_matrix[1][i] = y[i];
            state_matrix[2][i] = z[i];
            state_matrix[3][i] = vx[i];
            state_matrix[4][i] = vy[i];
            state_matrix[5][i] = vz[i];
        }

	Double_t Jacobi_matrix[6][6];
	Jacobi_matrice(Jacobi_matrix);


     // Define matrix to hold time derivatives of state vectors
        Double_t  state_dot_matrix[rows][cols-1];

     // Calculate time derivatives of state vectors
        for (int i = 0; i < cols-1; i++) {
    // Extract state vector at time t
            Double_t state_vector[6];
            for (int j = 0; j < 6; j++) {
                state_vector[j] = state_matrix[j][i];
            }
    // Multiply Jacobian matrix with state vector to get time derivative of state vector
            for (int j = 0; j < 6; j++) {
                state_dot_matrix[j][i] = 0;
                for (int k = 0; k < 6; k++) {
                    state_dot_matrix[j][i] += Jacobi_matrix[j][k] * state_vector[k];
                }
            }
        }



        double I[6][6] = {{0}};
        for (int i = 0; i < n; i++) {
            I[i][i] = 1.0;
         }


        cout << "propagator state matrix:" << endl;
	for (int i = 0; i < 6; i++) {
  	    for (int j = 0; j < cols-1; j++) {
    		cout << state_dot_matrix[i][j] << " ";
  	    }
            cout << endl;
        }



/*
	// For starting  kalman filter.

	//sytem parameters.

	Double_t Ex,Ey,Ez,f3;              //Electric field in V/m. 
        Double_t Bx,By,Bz,m1,m2;              // magnetic field in Tesla
        Double_t q = 3.2*pow(10,-19);   //charge of the particle in eV
        Double_t B=2.0;                 // Applied magnetic field. 
        m1 = 6.64*TMath::Power(10,-27);  // mass of the particle in kg
        Double_t E = TMath::Cos((q*B)/m) * 500 ; 
        Bx = 0;                      // magnetic field in x and y  direction in Tesla.
        By = 0;
        Bz = B;          // magnetic field in the z direction in Tesla.
        Ex = 0;                      // Electric field in the x and  direction in V/m. 
        Ey = 0;
        Ez = -E; 
	m2=q/m1;

	//DEfine all matrices
	TMatrixD  height_1(6,1);
	TMatrixD Jacobian(6,6);
        TMatrixD F_1(3,3);      //state transition matrix.
        TMatrixD G_1(3,1);	
	TMatrixD U_1(3,1);  
        TMatrixD X_1(3,1);	//State Estimate matrix.
        TMatrixD Y_1(3,1);	//Observation matrix or measurements.
        TMatrixD H_1(3,3);	//Measurement matrix
        TMatrixD R_1(3,3);	//Measurement Noise covariance.
        TMatrixD C_1(3,3);	
        TMatrixD ax_1(3,1);	
        TMatrixD P_1(3,3);	//Estimate Error Covariance
        TMatrixD Q_1(3,3);	//Process Noise Covariance.
        TMatrixD I_1(6,6);	//Identity matrix.

        //kalman storage
        TMatrixD x_pred(3,1);
        TMatrixD P_pred(3,3);
        TMatrixD K(3,3);	//Kalman Gain.
        TMatrixD Z(3,1);
 
	//filling the matrices
	Double_t Matrix_height[6] = {h_z,h_z,h_z,h,h,h};
        height_1.Use(height_1.GetNrows(), height_1.GetNcols(), Matrix_height);

	Double_t Matrix_JA[36] = {0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,(q/m)*Bz,(-q/m)*By,0,0,0,(-q/m)*Bz,0,(q/m)*Bx,0,0,0,(q/m)*(Ez+By),(q/m)*(Ez-Bx),0};
	Jacobian.Use(Jacobian.GetNrows(),Jacobian.GetNcols(),Matrix_JA);

        Double_t Matrix_F[9] = {0,1,1,1,0,1,1,1,0};
        F_1.Use(F_1.GetNrows(), F_1.GetNcols(), Matrix_F);

        Double_t Matrix_G[3] = {1,1,1};
        G_1.Use(G_1.GetNrows(), G_1.GetNcols(), Matrix_G);

	Double_t Matrix_U[3] = {m1*Ex,m1*Ey,m1*Ez};
        U_1.Use(U_1.GetNrows(), U_1.GetNcols(), Matrix_U);

	// Initial predicted state
	Double_t x_1[n],y_1[n],z_1[n],vx_1[n],vy_1[n],vz_1[n];
	vx_1[0] = 0.0001;
	vy_1[0] = 0.0001;
	vz_1[0] = 0.0001; 

        Double_t Matrix_X[3] = {vx_1[0],vy_1[0],vz_1[0]};
        X_1.Use(X_1.GetNrows(), X_1.GetNcols(), Matrix_X);


        Double_t Matrix_Y[3] ;
        Y_1.Use(Y_1.GetNrows(), Y_1.GetNcols(), Matrix_Y);

        Double_t Matrix_H[9] = {1,0,0,0,1,0,0,0,1} ;
        H_1.Use(H_1.GetNrows(), H_1.GetNcols(), Matrix_H);


        Double_t Matrix_C[9] = {1,0,0,0,1,0,0,0,1};
        C_1.Use(C_1.GetNrows(), C_1.GetNcols(), Matrix_C);

        Double_t I[36]  = {1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,0,0,0,0,0,1};
        I_1.Use(I_1.GetNrows(), I_1.GetNcols(), I);
	
	Double_t Matrix_Q[9] = {0,0,0,0,0,0,0,0,0};
        Q_1.Use(Q_1.GetNrows(), Q_1.GetNcols(), Matrix_Q);


//	calculating the propagator matrix F. 
	TMatrixD F(6,1); 
	//F[0] = 0.0;  // Initial values of F. 
	
	for (Int_t i = 0; i<4;i++){	
		F = h*Jacobian; (I_1 + F[i]*(t[i]-t[0])/h);
		F.Print();
	}






	//Calculating the variance of the process and measured moise. 
	Double_t sum=0.0  ,mean, standardDeviation=0.0,st,variance,covariance;
	Double_t sum1=0.0 , mean1, standardDeviation1=0.0,st1,variance1,covariance1;
	Double_t sum2=0.0, mean2, standardDeviation2=0.0,st2,variance2,covariance2;

	Double_t P_n[n],I_vx[n],K_G[n],measurement[10],estimate[10];
	Double_t P_ny[n],I_vy[n],K_G1[n];
	Double_t P_nz[n],I_vz[n],K_G2[n];

	for(Int_t i = 0; i < n-1; ++i) {
		sum += vx[i];
		sum1 += vy[i];
		sum2 += vz[i];
	}

	mean = sum / n-1;
	mean1 = sum1/n-1;
	mean2 = sum2/n-1;

	for(Int_t i = 0; i < n; ++i) {
		standardDeviation += pow(vx[i] - mean, 2);
		standardDeviation1 += pow(vy[i] - mean1,2);
		standardDeviation2 += pow(vz[i] - mean2,2);
	}

  	st = sqrt(standardDeviation / n-1);
	st1 = sqrt(standardDeviation1 / n-1);
	st2 = sqrt(standardDeviation2 / n-1);

	variance = st*st;
	variance1 = st1*st1;
	variance2 = st2*st2;


	P_n[0] =variance*variance;  
        P_ny[0] = variance1*variance1;
        P_nz[0] = variance2*variance2;

        Double_t r_x[n],r_y[n],r_z[n];
        r_x[0] = variance;
        r_y[0] = variance1;
        r_z[0] = variance2;



	//std::cout<< I_vx[0] <<" " <<mean2 <<endl;
	//Calculating the covariance between vx,vy and vz. 

	 for(Int_t i = 0; i < n-1; ++i) {
                sum += vx[i];
                sum1 += vy[i];
                sum2 += vz[i];
        }

        mean = sum / n-1;
        mean1 = sum1/n-1;
        mean2 = sum2/n-1;

        for(Int_t i = 0; i < n; ++i) {
                standardDeviation += (vx[i] - mean)* (vy[i] - mean1);
                standardDeviation1 += (vx[i] - mean)* (vz[i] - mean2);
                standardDeviation2 += (vy[i] - mean)* (vz[i] - mean1);;
        }

        covariance = (standardDeviation / n-1);
        covariance1 = (standardDeviation1 / n-1);
        covariance2 = (standardDeviation2 / n-1);



        P_n[0] =covariance;  
        P_ny[0] = covariance1;
        P_nz[0] = covariance2;

	//std::cout<<r_x[0]<<" "<<P_n[0]<<std::endl;



	Double_t Matrix_P[9] = {r_x[0],P_n[0],P_ny[0],P_n[0],r_y[0],P_nz[0],P_ny[0],P_nz[0],r_z[0]} ; // generated errors.
        P_1.Use(P_1.GetNrows(), P_1.GetNcols(), Matrix_P);

	Double_t Matrix_R[9] ={r_x[0],0,0,0,r_y[0],0,0,0,r_z[0]};          // Error in measurements.
        R_1.Use(R_1.GetNrows(), R_1.GetNcols(), Matrix_R);
	
	P_1.Print();
	for(Int_t i=0; i<n; i++){

        Double_t Matrix_Y[3] = {vx[i],vy[i],vz[i]};
        Y_1.Use(Y_1.GetNrows(), Y_1.GetNcols(), Matrix_Y);

        //start kalman

        x_pred = (F_1 * X_1) ;//+ (B_1*U_1) ;
        P_pred =  (F_1 *TMatrixD(P_1, TMatrixD::kMultTranspose,F_1))  + Q_1;


        //updates



        K =  TMatrixD(P_pred, TMatrixD::kMultTranspose,H_1) *  (H_1 * TMatrixD(P_pred, TMatrixD::kMultTranspose,H_1) + R_1).Invert();

        Z = C_1*Y_1;
        X_1 = x_pred + (K *(Z-(H_1*x_pred)));
        P_1 =(I_1-K*H_1)*P_pred;

        //x_pred.Print();
        //Z.Print();
        //std::cout<< z[n] << " " << x_pred << std::endl;
        }


        x_pred.Print();
        P_pred.Print();
	

*/
//	I_1.Print();


/*

  c1->cd(1);
  auto r1=new TGraph2D(n,vx,vy,vz);
  r1->SetTitle("Particle motion; X axis;Y axis; Z axis");
//  gStyle->SetPalette(1);
  //gr1->SetMarkerStyle(20);
  r1->Draw(" LINE");


   c1->cd(2);

   TGraph *gr1 = new TGraph (n, vz,vx);
   gr1->GetXaxis()->SetTitle("z position");
   gr1->GetYaxis()->SetTitle("x Position");
   gr1->GetXaxis()->CenterTitle();
   gr1->GetYaxis()->CenterTitle();
   gr1->Draw("AL");


 c1->cd(3);

   TGraph *g = new TGraph (n, vz,vy);
   g->GetXaxis()->SetTitle("z position");
   g->GetYaxis()->SetTitle("y Position");
   g->GetXaxis()->CenterTitle();
   g->GetYaxis()->CenterTitle();
   g->Draw("AL");

 c1->cd(4);

   TGraph *gn = new TGraph (n, vx,vy);
   gn->GetXaxis()->SetTitle("x position");
   gn->GetYaxis()->SetTitle("y Position");
   gn->GetXaxis()->CenterTitle();
   gn->GetYaxis()->CenterTitle();
   gn->Draw("AL");


*/

//elossfile.close();
}


