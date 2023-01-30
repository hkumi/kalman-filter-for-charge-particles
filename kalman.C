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
 

        f2 =  q/m * (Ey + vz*Bx - vx*Bz) ;//- s*TMath::Sin(po)*TMath::Sin(az); 

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
 

        f3 =  q/m * (Ez + vx*By - vy*Bx) ;//- s*TMath::Cos(po);
	//std::cout <<r<< " " << TMath::Power(vy,2) + TMath::Power(vx,2) + TMath::Power(vz,2)<< " " << TMath::Power(vz,2) <<  std::endl; 
        return f3;
        }

//  Functions for positions. 

Double_t fx(Double_t t,Double_t vx){
        return vx;
        }

Double_t fy(Double_t t,Double_t vy){
        return vy;
        }

Double_t  fz(Double_t t,Double_t vz){
        return vz;
        }



void kalman(){

 	TCanvas *c1 = new TCanvas("c1","Fast Fourier Transform",500,600);
        c1->Divide(2,2);

        const Int_t n = 1000;
	Double_t m = 6.64*TMath::Power(10,-27);
	Int_t v = 1;
	Double_t t[n],vx[n],vy[n],vz[n],c,beta,gamma;
	Double_t x[n],y[n],z[n];
	c = 2.99792*TMath::Power(10,8);
	gamma =1/TMath::Sqrt(1-TMath::Power(v/c,2));
	Double_t F1,F2,F3,F4=0.0;
        Double_t G1,G2,G3,G4=0.0;
	Double_t L1,L2,L3,L4=0.0;
	Double_t theta,phi;         // Using spherical coordinates.
        theta=2.83135;
	phi=-1.41593;

//Initial Conditions

        t[0]=0.0003;
        vx[0]= gamma*c*TMath::Sin(theta/180*M_PI)*TMath::Cos(phi/180*M_PI);
	//vx[0] = 1.0;
	//vy[0] = 1.0;
	//vz[0] = 1.0;
	double tf=8.18*TMath::Power(10,-7);
	double t0=0.0;
        vy[0]=gamma*c*TMath::Sin(theta/180*M_PI)*TMath::Sin(phi/180*M_PI);
	vz[0]=gamma*c*TMath::Cos(theta/180*M_PI);
	Double_t  h=(tf-t0)/n;                             // step_size in seconds
	//x[0] = vx[0] ;
	//y[0] = vy[0];
	//z[0] = vz[0];
	//Graph to evaluate Energy loss
	Double_t  gasMediumDensity_ = 0.0153236;   //g/cm3
	Double_t p[n],s[n],vv[n],bet[n],gamma1[n],Energy[n];
	TGraph* eLossCurve = new TGraph();

       energyloss("van.txt", eLossCurve);

	//eLossCurve->Draw();

	for(Int_t i=0;i <n-1;  i++){
  //       std::cout<<vx[i] << " "<<vy[i]<< " " <<vz[i]<<std::endl;
          F1 = dxdt(t[i],vx[i],vy[i],vz[i],s[i]);
          G1 = dydt(t[i],vx[i],vy[i],vz[i],s[i]);
	  L1 = dzdt(t[i],vx[i],vy[i],vz[i],s[i]);

          F2=  dxdt(t[i]+h*0.5,vx[i]+0.5*h*F1,vy[i]+G1*h*0.5,vz[i]+0.5*h*L1,s[i]);
          G2=  dydt(t[i]+h*0.5,vx[i]+F1*h*0.5,vy[i]+G1*h*0.5,vz[i]+0.5*h*L1,s[i]);
	  L2=  dzdt(t[i]+h*0.5,vx[i]+F1*h*0.5,vy[i]+G1*h*0.5,vz[i]+0.5*h*L1,s[i]);

          F3=  dxdt(t[i]+h*0.5,vx[i]+F2*h*0.5,vy[i]+G2*h*0.5, vz[i]+0.5*h*L2,s[i]);
          G3=  dydt(t[i]+h*0.5,vx[i]+F2*h*0.5,vy[i]+G2*h*0.5, vz[i]+0.5*h*L2,s[i]);
	  L3=  dzdt(t[i]+h*0.5,vx[i]+F2*h*0.5,vy[i]+G2*h*0.5, vz[i]+0.5*h*L2,s[i]);

          F4=  dxdt(t[i]+h,vx[i]+F3*h,vy[i]+h*G3,vz[i]+h*L3,s[i]);
          G4=  dydt(t[i]+h,vx[i]+h*F3,vy[i]+h*G3,vz[i]+h*L3,s[i]);
	  L4=  dzdt(t[i]+h,vx[i]+F3*h,vy[i]+G3*h,vz[i]+h*L3,s[i]);

          vx[i+1] = vx[i] + h/6 *( F1 + 2*F2 + 2*F3 + F4);
          vy[i+1]  = vy[i] + h/6 *( G1 + 2*G2 + 2*G3 + G4);
	  vz[i+1]  = vz[i] + h/6 *( L1 + 2*L2 + 2*L3 + L4);
          t[i+1]= t[i] + h;



          vx[i] = vx[i+1];
	  vz[i] = vz[i+1];
          vy[i] = vy[i+1];
          t[i] = t[i+1];

	  //std::cout<<vx[i]<<endl;

/*
	  vv[i] = TMath::Sqrt(TMath::Power(vx[i],2)+TMath::Power(vy[i],2)+TMath::Power(vz[i],2));
          bet[i] = vv[i]/c;
//          gamma1[i] = TMath::Sqrt(4/(TMath::Power(bet[i],2)));
  //        Energy[i] = (gamma1[i] - 1) *931.494028;
	  Energy[i] = 4 * 931.495*0.5*TMath::Power(bet[i],2);    // For non-relativity
          std::cout<< vv[i] << " " << bet[i] << std:: endl;
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
*/
        }
//	TGraph* gr = new TGraph(100,p,s);
  //      gr->Draw();

	// For starting  kalman filter.

	//DEfine all matrices

        TMatrixD A_1(3,3);      //state transition matrix.
//        TMatrixD B_1(2,1);	
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
        Double_t Matrix_A[9] = {1,0,0,0,1,0,0,0,1};
        A_1.Use(A_1.GetNrows(), A_1.GetNcols(), Matrix_A);

  //      Double_t Matrix_B[2] = {0.5*TMath::Power(dt,2),dt};
    //    B_1.Use(B_1.GetNrows(), B_1.GetNcols(), Matrix_B);

	// Initial predicted state
	Double_t vx_1[n],vy_1[n],vz_1[n];
	vx_1[0] = vx[0];
	vy_1[0] = vy[0];
	vz_1[0] = vz[0]; 

        Double_t Matrix_X[3] = {vx_1[0],vy_1[0],vz_1[0]};
        X_1.Use(X_1.GetNrows(), X_1.GetNcols(), Matrix_X);


        Double_t Matrix_Y[3] ;
        Y_1.Use(Y_1.GetNrows(), Y_1.GetNcols(), Matrix_Y);

        Double_t Matrix_H[9] = {1,0,0,0,1,0,0,0,1} ;
        H_1.Use(H_1.GetNrows(), H_1.GetNcols(), Matrix_H);

        Double_t Matrix_R[9] ={0.001,0,0,0,0.001,0,0,0,0.001};
        R_1.Use(R_1.GetNrows(), R_1.GetNcols(), Matrix_R);

        Double_t Matrix_C[9] = {1,0,0,0,1,0,0,0,1};
        C_1.Use(C_1.GetNrows(), C_1.GetNcols(), Matrix_C);

        Double_t I[9]  = {1,0,0,0,1,0,0,0,1};
        I_1.Use(I_1.GetNrows(), I_1.GetNcols(), I);
	
	Double_t Matrix_Q[9] = {0,0,0,0,0,0,0,0,0};
        Q_1.Use(Q_1.GetNrows(), Q_1.GetNcols(), Matrix_Q);





	Double_t sum = 0.0, mean, standardDeviation = 0.0,st,variance;
	Double_t sum1 = 0.0, mean1, standardDeviation1 = 0.0,st1,variance1;
	Double_t sum2 = 0.0, mean2, standardDeviation2 = 0.0,st2,variance2;

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

	for(Int_t i = 0; i < 10; ++i) {
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



	//std::cout<< I_vx[0] <<" " <<mean2 <<endl;
	P_n[0] =variance*variance;  
	P_ny[0] = variance1*variance1;
	P_nz[0] = variance2*variance2;

	Double_t Matrix_P[9] = {P_n[0],0,0,0,P_ny[0],0,0,0,P_nz[0]} ; // generated errors.
        P_1.Use(P_1.GetNrows(), P_1.GetNcols(), Matrix_P);


	for(Int_t i=0; i<1; i++){

        Double_t Matrix_Y[3] = {vx[i],vy[i],vz[i]};
        Y_1.Use(Y_1.GetNrows(), Y_1.GetNcols(), Matrix_Y);

        //start kalman

        x_pred = (A_1 * X_1) ;//+( (B_1*ax) + Q) ;
        P_pred =  (A_1 *TMatrixD(P_1, TMatrixD::kMultTranspose,A_1))  + Q_1;


        //updates

        K =  TMatrixD(P_pred, TMatrixD::kMultTranspose,H_1) *  (H_1 * TMatrixD(P_pred, TMatrixD::kMultTranspose,H_1) + R_1).Invert();

        Z = C_1*Y_1;
        X_1 = x_pred + (K *(Z-(H_1*x_pred)));
        P_1 =(I_1-K*H_1)*P_pred;

        //x_pred.Print();
        //Z.Print();
        //std::cout<< z[n] << " " << x_pred << std::endl;
        }


        P_pred.Print();
        K.Print();
	



	return 0;
	


/*
  c1->cd(1);
  auto r1=new TGraph2D(n,I_vx,I_vy,I_vz);
  r1->SetTitle("Particle motion; X axis;Y axis; Z axis");
//  gStyle->SetPalette(1);
//  gr1->SetMarkerStyle(20);
  r1->Draw(" LINE");
*/

   c1->cd(2);

   TGraph *gr1 = new TGraph (10, t,estimate);
   gr1->GetXaxis()->SetTitle("vz position");
   gr1->GetYaxis()->SetTitle("vy Position");
   gr1->GetXaxis()->CenterTitle();
   gr1->GetYaxis()->CenterTitle();
   gr1->Draw("AL");
/*
 c1->cd(3);

   TGraph *g = new TGraph (n, vz,vy);
   g->GetXaxis()->SetTitle("vz position");
   g->GetYaxis()->SetTitle("vy Position");
   g->GetXaxis()->CenterTitle();
   g->GetYaxis()->CenterTitle();
   g->Draw("AL");

 c1->cd(4);

   TGraph *gn = new TGraph (n, vx,vy);
   gn->GetXaxis()->SetTitle("vx position");
   gn->GetYaxis()->SetTitle("vz Position");
   gn->GetXaxis()->CenterTitle();
   gn->GetYaxis()->CenterTitle();
   gn->Draw("AL");

*/

//elossfile.close();
}


