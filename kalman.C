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



Double_t  dxdt(double t,double vx, double vy, Double_t vz,Double_t s){

	Double_t Ex,Ey,Ez,f1;            	//Electric field in V/m. 
	Double_t Bx,By,Bz;	      	// magnetic field in Tesla
	Double_t rr,az,po;
	Double_t q = 3.2*TMath::Power(10,-19); 	//charge of the particle in eV
	Double_t m = 6.64*TMath::Power(10,-27); 	// mass of the particle in kg
	Double_t B=2.0;                 // Applied magnetic field.
        Double_t E=TMath::Cos((q*B)/m) * 500;   // Applied Electric field.
        Bx = 0;
	By = 0;                	// magnetic field in x and y  direction in Tesla.
	Bz = B;          // magnetic field in the z direction in Tesla.
        Ex = 0;
	Ey = 0;			// Electric field in the x and  direction in V/m.
	Ez = -E;                        // Electric field in the z direction. 
	
	rr = TMath::Sqrt(TMath::Power(vx,2)+TMath::Power(vy,2)+TMath::Power(vz,2));
	az = TMath::ATan(vy/vx);
	po = TMath::ACos(vz/rr);


	f1 =  q/m * (Ex + vy*Bz-vz*By); //-s*TMath::Sin(po)*TMath::Cos(az) ;                                    //dxdt with energyloss compensation.

	//std::cout<<Ex <<" "<< vy<<" "<<vz << " " <<By<< " "<< f1<<std::endl; 
        return f1;

        }

double  dydt(double t,double vx, double vy, Double_t vz,Double_t s){
	Double_t Ex,Ey,Ez,f2;              //Electric field in V/m. 
        Double_t Bx,By,Bz;              // magnetic field in Tesla
        Double_t q = 3.2*TMath::Power(10,-19);   //charge of the particle in eV
        Double_t B=2.0;                 // Applied magnetic field.
        Double_t rr,az,po;   
        Double_t m = 6.64*TMath::Power(10,-27);  // mass of the particle in kg
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

double  dzdt(double t,double vx, double vy,Double_t vz, Double_t s){
	Double_t Ex,Ey,Ez,f3;              //Electric field in V/m. 
        Double_t Bx,By,Bz;              // magnetic field in Tesla
        Double_t q = 3.2*pow(10,-19);   //charge of the particle in eV
        Double_t B=2.0;                 // Applied magnetic field.
        Double_t rr,az,po;  
        Double_t m = 6.64*TMath::Power(10,-27);  // mass of the particle in kg
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
 

        f3 =  q/m * (Ez + vx*By - vy*Bx);// - s*TMath::Cos(po);
	//std::cout <<r<< " " << TMath::Power(vy,2) + TMath::Power(vx,2) + TMath::Power(vz,2)<< " " << TMath::Power(vz,2) <<  std::endl; 
        return f3;
        }

//  Functions for positions. 

Double_t fx( Double_t vx,Double_t vy, Double_t vz){
        return vx/vz;
        }

Double_t fy(Double_t vx,Double_t vy, Double_t vz){
        return vy/vz;
        }

Double_t  dt_xdz(Double_t vx, Double_t vy, Double_t vz){
	Double_t k,q,p,tx,ty,Bz,By,Bx,dt1;
	k = pow(2.9979248,17);
	Double_t m = 6.64*TMath::Power(10,-27);  // mass of the particle in kg
	q = 3.2*pow(10,-19);   //charge of the particle in eV
	p= 1.0;
	Double_t B=2.0;                 // Applied magnetic field.
	Bx = 0;                      // magnetic field in x and y  direction in Tesla.
        By = 0;
        Bz = B;          // magnetic field in the z direction in Tesla.
	tx = vx/vz;
	ty = vy/vz;

	dt1 = k*(q/m*vx)*TMath::Sqrt(pow(tx,2)+pow(ty,2)+1)*(ty*Bz-(1+pow(tx,2))*By+tx*ty*By);

        return dt1;
        }


Double_t  dt_ydz(Double_t vx, Double_t vy, Double_t vz){
        Double_t k,q,p,tx,ty,Bz,By,Bx,dt;
        k = pow(2.9979248,17);
        q = 3.2*pow(10,-19);   //charge of the particle in eV
	Double_t m = 6.64*TMath::Power(10,-27);  // mass of the particle in kg
        p= 1.0;
        Double_t B=2.0;                 // Applied magnetic field.
        Bx = 0;                      // magnetic field in x and y  direction in Tesla.
        By = 0;
        Bz = B;          // magnetic field in the z direction in Tesla.
        tx = vx/vz;
        ty = vy/vz;

        dt = k*(q/m*vy)*TMath::Sqrt(pow(tx,2)+pow(ty,2)+1)*(-tx*Bz-(1+pow(ty,2))*Bx-tx*ty*By);

        return dt;
        }




void kalman(){

 	TCanvas *c1 = new TCanvas("c1","Fast Fourier Transform",500,600);
        c1->Divide(2,2);

        const Int_t n = 1000;
	Double_t m = 6.64*TMath::Power(10,-27);
	Int_t v = 1;
	Double_t t[n],vx[n],vy[n],vz[n],c,beta,gamma;
	Double_t x[n],y[n],z[n],z0,zf;
	c = 2.99792*TMath::Power(10,8);
	gamma =1/TMath::Sqrt(1-TMath::Power(v/c,2));
	Double_t F1[n],F2[n],F3[n],F4[n];
        Double_t G1[n],G2[n],G3[n],G4[n];
	Double_t L1[n],L2[n],L3[n],L4[n];
	Double_t theta,phi;         // Using spherical coordinates.
        theta=2.83135;
	phi=-1.41593;

//Initial Conditions

        t[0]=0.0003;
        //vx[0]= gamma*c*TMath::Sin(theta/180*M_PI)*TMath::Cos(phi/180*M_PI);
	//vx[0] = 1.0;
	//vy[0] = 1.0;
	//vz[0] = 1.0;
	double tf=8.18*TMath::Power(10,-7);
	double t0=0.0;



        vy[0]=gamma*c*TMath::Sin(theta/180*M_PI)*TMath::Sin(phi/180*M_PI);
	vz[0]=gamma*c*TMath::Cos(theta/180*M_PI);
	Double_t  h=(tf-t0)/n;
                             // step_size in seconds
	x[0] = 0.000001;
	y[0] = 0.000001;
	z[0] = 0.000001;

	z0 = 0.000001;
	zf = 20.0;
	Double_t h_z = (zf-z0)/n;
	//Graph to evaluate Energy loss
	Double_t  gasMediumDensity_ = 0.0153236;   //g/cm3
	Double_t p[n],s[n],vv[n],bet[n],gamma1[n],Energy[n];
	TGraph* eLossCurve = new TGraph();

       energyloss("van.txt", eLossCurve);

	//eLossCurve->Draw();

	std::cout<<x[0] << " " << y[0] << " " << z[0] << std::endl;
	//Start Rk4. 
/*

	for(Int_t i=0;i <n-1;  i++){
  //       std::cout<<vx[i] << " "<<vy[i]<< " " <<vz[i]<<std::endl;
          F1[i] = dxdt(t[i],vx[i],vy[i],vz[i],s[i]);
          G1[i] = dydt(t[i],vx[i],vy[i],vz[i],s[i]);
	  L1[i] = dzdt(t[i],vx[i],vy[i],vz[i],s[i]);

          F2[i] = dxdt(t[i]+h*0.5,vx[i]+0.5*h*F1[i],vy[i]+G1[i]*h*0.5,vz[i]+0.5*h*L1[i],s[i]);
          G2[i] = dydt(t[i]+h*0.5,vx[i]+F1[i]*h*0.5,vy[i]+G1[i]*h*0.5,vz[i]+0.5*h*L1[i],s[i]);
	  L2[i] = dzdt(t[i]+h*0.5,vx[i]+F1[i]*h*0.5,vy[i]+G1[i]*h*0.5,vz[i]+0.5*h*L1[i],s[i]);

          F3[i] = dxdt(t[i]+h*0.5,vx[i]+F2[i]*h*0.5,vy[i]+G2[i]*h*0.5, vz[i]+0.5*h*L2[i],s[i]);
          G3[i] = dydt(t[i]+h*0.5,vx[i]+F2[i]*h*0.5,vy[i]+G2[i]*h*0.5, vz[i]+0.5*h*L2[i],s[i]);
	  L3[i] = dzdt(t[i]+h*0.5,vx[i]+F2[i]*h*0.5,vy[i]+G2[i]*h*0.5, vz[i]+0.5*h*L2[i],s[i]);

          F4[i] = dxdt(t[i]+h,vx[i]+F3[i]*h,vy[i]+h*G3[i],vz[i]+h*L3[i],s[i]);
          G4[i] = dydt(t[i]+h,vx[i]+h*F3[i],vy[i]+h*G3[i],vz[i]+h*L3[i],s[i]);
	  L4[i] = dzdt(t[i]+h,vx[i]+F3[i]*h,vy[i]+G3[i]*h,vz[i]+h*L3[i],s[i]);

          vx[i+1] = vx[i] + h/6 *( F1[i] + 2*F2[i] + 2*F3[i] + F4[i]);
          vy[i+1]  = vy[i] + h/6 *( G1[i] + 2*G2[i] + 2*G3[i] + G4[i]);
	  vz[i+1]  = vz[i] + h/6 *( L1[i] + 2*L2[i] + 2*L3[i] + L4[i]);
          t[i+1]= t[i] + h;



          vx[i] = vx[i+1];
	  vz[i] = vz[i+1];
          vy[i] = vy[i+1];
          t[i] = t[i+1];

	  //std::cout<<vx[i]<<endl;


	  vv[i] = TMath::Sqrt(TMath::Power(vx[i],2)+TMath::Power(vy[i],2)+TMath::Power(vz[i],2));
          bet[i] = vv[i]/c;
//          gamma1[i] = TMath::Sqrt(4/(TMath::Power(bet[i],2)));
  //        Energy[i] = (gamma1[i] - 1) *931.494028;

	  Energy[i] = 4 * 931.495*0.5*TMath::Power(bet[i],2);    // For non-relativity
          //std::cout<< vv[i] << " " << bet[i] << std:: endl;
          if (Energy[i] < 0.01)
         {

		p[i] = Energy[i];

	  return 1;


	}




          //p[i] = Energy[i];
          s[i] = gasMediumDensity_*eLossCurve->Eval(p[i]);
          s[i] = s[i]* 1.6021773349*TMath::Power(10,-13); // convert MeV to J (kg.m^2.s^-2)
          s[i] = s[i]* gasMediumDensity_*100; // multiply with density to get kg.m.s^-2 (Force)
          s[i] = s[i]/ m; // divide by mass to get homogeneous with acceleration (m.s^-2)
          //std::cout<<" Energy : "<<gamma1[i]<<" - dE/dx :     "<<eLossCurve->Eval(p[i])<<" - dE/dx (density) :"<<s[i]<<"\n";
         // std::cout<<" Energy (Curves) : "<<s[i]<<" - dE/dx :     "<<" - dE/dx (density) :"<<Energy[i]<<"\n";

          //if(elossfile->eof()) break;


        }
*/

//	TGraph* gr = new TGraph(100,p,s);
  //      gr->Draw();
	//RK4 for tracking the position of the particle.
	Double_t K1[n],K2[n],K3[n],K4[n];
        Double_t M1[n],M2[n],M3[n],M4[n];
        Double_t H1[n],H2[n],H3[n],H4[n];
	Double_t N1[n],N2[n],N3[n],N4[n];

	for(Int_t i=0;i <n-1;  i++){
  //       std::cout<<vx[i] << " "<<vy[i]<< " " <<vz[i]<<std::endl;
          K1[i] = fx(x[i],y[i],z[i]);
          M1[i] = fy(x[i],y[i],z[i]);
          H1[i] = dt_xdz(vx[i],vy[i],vz[i]);
	  N1[i] = dt_ydz(vx[i],vy[i],vz[i]);

          K2[i] = fx(x[i]+0.5*h_z*K1[i],y[i]+M1[i]*h_z*0.5,z[i]+0.5*h_z);
          M2[i] = fy(x[i]+K1[i]*h_z*0.5,y[i]+M1[i]*h_z*0.5,z[i]+0.5*h_z);
          H2[i] = dt_xdz(vx[i]+0.5*h_z*H1[i],vy[i]+N1[i]*h_z*0.5,vz[i]+0.5*h_z);
	  N2[i] = dt_ydz(vx[i]+0.5*h_z*H1[i],vy[i]+N1[i]*h_z*0.5,vz[i]+0.5*h_z);

          K3[i] = fx(x[i]+K2[i]*h_z*0.5, y[i]+M2[i]*h_z*0.5,z[i]+0.5*h_z);
          M3[i] = fy(x[i]+K2[i]*h_z*0.5,y[i]+M2[i]*h_z*0.5,z[i]+0.5*h_z);
          H3[i] = dt_xdz(vx[i]+0.5*h_z*H2[i],vy[i]+N2[i]*h_z*0.5, vz[i]+0.5*h_z);
	  N3[i] = dt_ydz(vx[i]+0.5*h_z*H2[i],vy[i]+N2[i]*h_z*0.5,vz[i]+0.5*h_z);

          K4[i] = fx(x[i]+K3[i]*h_z,y[i]+M3[i]*h_z,z[i]+h_z);
          M4[i] = fy(x[i]+h_z*K3[i],y[i]+M3[i]*h_z,z[i]+h_z);
          H4[i] = dt_xdz(vx[i]+h_z*H3[i],vy[i]+N3[i]*h_z,vz[i]+h_z);
	  N4[i] = dt_ydz(vx[i]+h_z*H3[i],vy[i]+N3[i]*h_z,vz[i]+h_z);	 


          x[i+1] = x[i] + h_z/6 *( K1[i] + 2*K2[i] + 2*K3[i] + K4[i]);
          y[i+1]  = y[i] + h_z/6 *( M1[i] + 2*M2[i] + 2*M3[i] + M4[i]);
          z[i+1]  = z[i] + h_z;
          //t[i+1]= t[i] + h;

	  vx[i+1] = vx[i] + h_z/6 *( H1[i] + 2*H2[i] + 2*H3[i] + H4[i]);
          vy[i+1]  = vy[i] + h_z/6 *( N1[i] + 2*N2[i] + 2*N3[i] + N4[i]);
          vz[i+1]  = vz[i] + h_z;


          x[i] = x[i+1];
          z[i] = z[i+1];
          y[i] = y[i+1];
          //t[i] = t[i+1];

	  vx[i] = vx[i+1];
          vz[i] = vz[i+1];
          vy[i] = vy[i+1];


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
        TMatrixD I_1(3,3);	//Identity matrix.

        //kalman storage
        TMatrixD x_pred(3,1);
        TMatrixD P_pred(3,3);
        TMatrixD K(3,3);	//Kalman Gain.
        TMatrixD Z(3,1);
 
	//filling the matrices
        Double_t Matrix_F[9] = {0,1,1,1,0,1,1,1,0};
        F_1.Use(F_1.GetNrows(), F_1.GetNcols(), Matrix_F);

        Double_t Matrix_G[3] = {1,1,1};
        G_1.Use(G_1.GetNrows(), G_1.GetNcols(), Matrix_G);

	Double_t Matrix_U[3] = {m1*Ex,m1*Ey,m1*Ez};
        U_1.Use(U_1.GetNrows(), U_1.GetNcols(), Matrix_U);

	// Initial predicted state
	Double_t vx_1[n],vy_1[n],vz_1[n];
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

        Double_t I[9]  = {1,0,0,0,1,0,0,0,1};
        I_1.Use(I_1.GetNrows(), I_1.GetNcols(), I);
	
	Double_t Matrix_Q[9] = {0,0,0,0,0,0,0,0,0};
        Q_1.Use(Q_1.GetNrows(), Q_1.GetNcols(), Matrix_Q);




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
	



	return 0;
	
*/


//  c1->cd(1);
  //auto r1=new TGraph2D(n,x,y,z);
  //r1->SetTitle("Particle motion; X axis;Y axis; Z axis");
//  gStyle->SetPalette(1);
//  gr1->SetMarkerStyle(20);
  //r1->Draw(" LINE");


   c1->cd(2);

   TGraph *gr1 = new TGraph (n, x,y);
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




//elossfile.close();
}


