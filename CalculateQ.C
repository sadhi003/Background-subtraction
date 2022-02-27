// run this macro as root -l -b -q 'QValue_final.C("rootFile_40Cut_9.root","/Users/shankar/Expt_work/QValue", 1, 2)'

#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif
#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooDataHist.h"
#include "RooGaussian.h"
#include "RooConstVar.h"
#include "RooAddPdf.h"
#include "RooChebychev.h"
#include "RooFitResult.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "RooPlot.h"
#include "TFile.h"
#include "TStyle.h"
#include "TH2.h"
#include "TMatrixDSym.h"
#include "TMath.h"
#include "TRandom.h"
#include "TH1D.h"
#include "TTree.h"
#include <cmath>
#include <iostream>
#include <string>
#include <cstdlib>
#include <iomanip>
#include <fstream>
#include <map>
#include <sstream>
#include <stdexcept>
#include <vector>
#include "TMatrixT.h"
#include "TVectorD.h"
#include "TGraph.h"
#include "TFile.h"
#include "TH1.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include "TSpline.h"
#include "TMath.h"
using namespace RooFit ;


void CalculateQ(Int_t i_e,Int_t f_e)  //initial event 
{
	// read data file to get information for Q factor
    TString workdir ="./";

	TString outdir = workdir + "RootQ/";
	string makeoutdir("mkdir " +  outdir);
	system(makeoutdir.c_str());

	TString plotDir = workdir + "Plots/";
	string makeoutdir1("mkdir " + plotDir);
	system(makeoutdir1.c_str());


	TString indir1 ="./";
	TString infilename1= indir1 + "test.root";
	TChain *chain = new TChain("output");
	chain->Add(infilename1);

	Double_t pi = TMath::Pi();
	Int_t N_c = 500;

	Double_t i_f = 1.1;
	Double_t f_f = 1.5;
	Double_t f_Epx = 0.005;

	Int_t nEvent = chain->GetEntries();

	std::cout << "Number of events: " << nEvent << std::endl;

	// how many event need to run 
	TString s_nEvent = Form("%i",nEvent);
	Int_t digits = s_nEvent.Length();
	
	TString s_i_e = Form("%i",i_e);
	Int_t l_s_i = s_i_e.Length();
	while (l_s_i < digits)
	{
		s_i_e = "0"+s_i_e;
		l_s_i = s_i_e.Length();	  
	}
	
	TString s_f_e = Form("%i",f_e);
	Int_t l_s_f = s_f_e.Length();
	while (l_s_f < digits)
	{
		s_f_e = "0"+s_f_e;
		l_s_f = s_f_e.Length();	  
	}
	
	TString fileout = outdir + "CalculateQ" + "_" + s_i_e + "_" + s_f_e + ".root";


	TLorentzVector Fb;
	TLorentzVector Fpr;
	TLorentzVector Fk;
	TLorentzVector Fp;
	TLorentzVector Fpi;
	TLorentzVector Flam;
	TVector3 Tb, Tpr,Tp,Tpi,Tk;
	float Xpr, Ypr,Zpr;
	float beta1_x,beta1_y,beta1_z;
	float beta2_x,beta2_y,beta2_z;

	Float_t Mb = 0.0;
	Float_t Mpr = 0.938272;
	Float_t Mk = 0.493667;
	Float_t MDel = 1.220;
	Float_t Mpi = 0.13957;
	Float_t Mp = 0.938272;
	
	

/* #############################  descriptions of input data #########################*/

	chain->SetBranchStatus("*",0);

	// initializing tuples
	double ebeam,cosTheta_Ks_hel,phi_Ks_hel,cosTheta_Ks_cm,phi_Ks_cm,cosTheta_delta_cm,phi_delta_cm,cosTheta_Ks_gj,phi_Ks_gj;
	Double_t cosTheta_proton_hel, phi_proton_hel, cosTheta_pip_gj, phi_pip_gj;
	double men_t;
	double mppip1_kinfit, mKsKm_kinfit;
	double time_weight;
	

	int pol;

	float delta_rf;

	chain -> SetBranchStatus("ebeam",1);
	chain -> SetBranchStatus("mppip1_kinfit", 1);
	chain -> SetBranchStatus("mKsKm_kinfit", 1);
	chain -> SetBranchStatus("cosTheta_delta_cm", 1);
	chain -> SetBranchStatus("cosTheta_Ks_hel",1);
	chain -> SetBranchStatus("phi_Ks_hel",1);
	chain -> SetBranchStatus("men_t",1);
	chain -> SetBranchStatus("time_weight",1);
	chain -> SetBranchStatus("pol",1);
	chain -> SetBranchStatus("cosTheta_pip_gj",1);
	chain -> SetBranchStatus("phi_pip_gj",1);

	
	
    
	// getting branch info

	chain -> SetBranchAddress("ebeam", &ebeam);
	chain -> SetBranchAddress("cosTheta_delta_cm", &cosTheta_delta_cm);	
	chain -> SetBranchAddress("mppip1_kinfit", &mppip1_kinfit);
	chain -> SetBranchAddress("mKsKm_kinfit", &mKsKm_kinfit);
	chain -> SetBranchAddress("cosTheta_Ks_hel", &cosTheta_Ks_hel);
	chain -> SetBranchAddress("phi_Ks_hel", &phi_Ks_hel);
	chain -> SetBranchAddress("men_t", &men_t);
	chain -> SetBranchAddress("time_weight", &time_weight);
	chain -> SetBranchAddress("pol", &pol);
	chain -> SetBranchAddress("cosTheta_pip_gj", &cosTheta_pip_gj);
	chain -> SetBranchAddress("phi_pip_gj", &phi_pip_gj);





    /*######################  end input file description ######################*/

	//############################## output file description  #####################
	TFile fout(fileout,"recreate");
	TTree *tree = new TTree("output_q","Description of data set");


	

	

	float ebeam_q, mppip1_kinfit_q, mKsKm_kinfit_q, cosTheta_delta_cm_q, cosTheta_Ks_hel_q, phi_Ks_hel_q;
	float men_t_q;
	float time_weight_q;
	Int_t pol_q;
	float cosTheta_pip_gj_q,phi_pip_gj_q;
	

	tree->Branch("ebeam_q",&ebeam_q, "ebeam_q/F");
	tree->Branch("men_t_q",&men_t_q, "men_t_q/F");
	tree->Branch("mppip1_kinfit_q",&mppip1_kinfit_q, "mppip1_kinfit_q/F");
	tree->Branch("mKsKm_kinfit_q",&mKsKm_kinfit_q, "mKsKm_kinfit_q/F");
	tree->Branch("cosTheta_delta_cm_q",&cosTheta_delta_cm_q, "cosTheta_delta_cm_q/F");
	tree->Branch("cosTheta_Ks_hel_q",&cosTheta_Ks_hel_q, "cosTheta_Ks_hel_q/F");
	tree->Branch("phi_Ks_hel_q",&phi_Ks_hel_q, "phi_Ks_hel_q/F");
	tree->Branch("time_weight_q",&time_weight_q, "time_weight_q/F");
	tree->Branch("pol_q",&pol_q, "pol_q/I");
	tree->Branch("cosTheta_pip_gj_q",&cosTheta_pip_gj_q, "cosTheta_pip_gj_q/F");
	tree->Branch("phi_pip_gj_q",&phi_pip_gj_q, "phi_pip_gj_q/F");
	
	

	

	float d_mass[300],dpar_Mat[5];
	Double_t evnt_mass,sigVal,meanVal, Pol_a0,Pol_a1,q_fac,S_Fun,B_Fun,Total_Fun,dQ_dmean, dQ_dsig,dQ_da0,dQ_da1,dQ_df,Q_Err,f;

	tree->Branch("evnt_mass",&evnt_mass,"evnt_mass/D");
    tree->Branch("d_mass",&d_mass,"d_mass[300]/F");  // The closest 100 events mass to the seed event 
    tree->Branch("sigVal",&sigVal,"sigVal/D");
    tree->Branch("Pol_a0",&Pol_a0,"Pol_a0/D");
    tree->Branch("Pol_a1",&Pol_a1,"Pol_a1/D");         
    tree->Branch("S_Fun",&S_Fun,"S_Fun/D");            // Value of signal function at event mass bin
    tree->Branch("B_Fun",&B_Fun,"B_Fun/D");            // Value of background function at event mass bin
    tree->Branch("Total_Fun",&Total_Fun,"Total_Fun/D");
    tree->Branch("f",&f,"f/D");                        // Fraction of signal events
    tree->Branch("q_fac",&q_fac,"q_fac/D");            // Q-factor 
    tree->Branch("dpar_Mat",&dpar_Mat,"dpar_Mat[4]/F");// Derivatives of variables
    tree->Branch("Q_Err",&Q_Err,"Q_Err/D");



    float COSth_p_pi;

	TMatrixD Xi(0,nEvent-1,0,4);  // 1 -> rowlb, 2 -> rowub, 3 -> collb, 4 -> colub

	for (Int_t i = 0; i < nEvent ; i++)
	{
		chain->GetEntry(i); 
		if(!(i%(nEvent/10))) std::cout << "Tree done " << i << " out of " << nEvent-1 << " ==> " << TMath::Nint(double(i)*100.0/double(nEvent)) << "%" << std::endl;


			//cout << cosTheta_delta_cm << endl;


		   Xi(i,0) = cosTheta_delta_cm;
		   Xi(i,1) = cosTheta_pip_gj;
		   Xi(i,2) = phi_pip_gj;
		   //Xi(i,3) = COSth_p_pi;
		   Xi(i,3) = men_t;



	}

	

	TArrayD d2(nEvent);
	Int_t s_e = f_e-i_e+1;

	for (Int_t i = i_e; i <= f_e; i++)
	{
		chain->GetEntry(i);
		evnt_mass = mppip1_kinfit;

		std::cout << "Calculating " << i+1-i_e << " out of " << s_e << " ==> " << double(i+1-i_e)*100.0/double(s_e) << "%" << std::endl;
		
		TTree *taux = new TTree("Aux","Auxiliar tree");
		Double_t mppip1_new;
		
		taux->Branch("mppip1_new",&mppip1_new,"mppip1_new/D");

		for (Int_t j = 0; j < nEvent; j++)
		{
		  d2.AddAt(pow((Xi(i,0)-Xi(j,0))/2,2) + pow((Xi(i,1)-Xi(j,1))/2,2) + pow((Xi(i,2)-Xi(j,2))/(2*pi),2) + pow((Xi(i,3)-Xi(j,3))/0.5,2), j);  // divide by range
		

		}

		const Int_t n = d2.fN;

		cout << "d2 dimension " << n << endl;


		Int_t *index = new int[n];
		TMath::Sort(n,d2.fArray,index,0);

		for (Int_t k = 0; k < N_c; k++)
		{
			chain->GetEntry(index[k]);  // take how many nearest neighbor need to use

			//cout <<"which index is this " << index[k]<< endl;
			Double_t mppip1 = chain->GetLeaf("mppip1_kinfit")->GetValue(0);

			if(!((k+1)%(N_c/10))) std::cout << "Selecting events for event " << i << ", " << k+1  << " out of " << N_c << " ==> " << double(k+1)*100.0/double(N_c) << "%" << std::endl;
			mppip1_new = mppip1;
			//d_mass[i] = mm2KP_new;

			taux->Fill();

		}

		Double_t f_bins = 25;// (f_f-i_f)/f_Epx;


		TString *s_i = new TString(Form("%d", i));
		TString *s_Nc = new TString(Form("%d", N_c));

		TString *f_name = new TString("E"+*s_i+"-C"+*s_Nc);

		TH1F *hist = new TH1F("hist",*f_name,f_bins,i_f,f_f);
		TCanvas *f_canvas = new TCanvas(*f_name);
		taux->Draw("mppip1_new>>hist");

		TH1* hh = hist;
      	RooRealVar x("x","M(Delta)",i_f,f_f) ;
      //--------------------------------------------------------------------------------
      // Building composite models with fractions(or plain non-extended composite mode)
      // -------------------------------------------------------------------------------
      	
      RooDataHist data("data","data",x,Import(*hh));

      RooPlot* xframe = x.frame(Title("M(Delta)"));
	  data.plotOn(xframe);
	  double Mean1 = 1.220;
	  double Width1 = 0.006;
	  double Sigma = 0.120;
      // Create a Gaussian PDFs g1(x,mean,sigma1) parameters
	  RooRealVar mean("mean1","mean of gaussians",Mean1, 1.18, 1.26);
      RooRealVar sigma1("sigma1","sigma of gaussians",Sigma, 0.001, 0.2);
      RooRealVar width("width","width of bw",Width1);
      //  RooGaussian sig1("sig1","Signal component 1",x,mean1,sigma1);
      RooVoigtian signal("signal","voigt function", x, mean, width,sigma1);
      // Build c polynomial p.d.f.  
      RooRealVar a0("a0","a0",0.5,-10.0,10.0);
      //      RooRealVar a1("a1","a1",-0.5,-100.0,100.0);
      //RooRealVar a2("a2","a2",0.1, -100, 100);
      RooChebychev bkg("bkg","background maximum Likelihood fitp.d.f.",x,RooArgList(a0)) ;

	  // Sum the composite signal and background 
	  RooRealVar bkgfrac("bkgfrac","fraction of background",0.57,0.,1.);
	  RooAddPdf  model("model","g1+a",RooArgList(signal,bkg),bkgfrac);    

	  RooFitResult* r = model.fitTo(data,Save());
      
      meanVal=mean.getVal();
      sigVal=sigma1.getVal();    
      Pol_a0=a0.getVal();    
      //      Pol_a1=a1.getVal(); 
      
      RooAbsReal* d_a0=bkg.derivative(a0,1); 
      //      RooAbsReal* d_a1=bkg.derivative(a1,1); 
      RooAbsReal* d_sig=signal.derivative(sigma1,1);
      RooAbsReal* d_mean=signal.derivative(mean,1);    
      
      
      q_fac = 0.0;   
    
      if(evnt_mass > i_f && evnt_mass < f_f){

				x.setVal(evnt_mass);// At the seed mass

				f=bkgfrac.getVal();

				S_Fun=signal.getVal(x); 
				B_Fun=bkg.getVal(x); 

				float b = (1.- f) * B_Fun;
				float s =  f * S_Fun;

				

				q_fac= s/(s + b);

				if(s < 1e-4) q_fac = s;
				if(TMath::IsNaN(q_fac)) q_fac = -0.5;
      
      }

      Q_Err = 0.0;
           
      if(evnt_mass > i_f && evnt_mass < f_f){

		      const TMatrixDSym& cov = r->covarianceMatrix();

		      x.setVal(evnt_mass);// At the seed mass
      
		      Double_t D_a0=d_a0->getVal(x); 
		      // Double_t D_a1=d_a1->getVal(x); 
		      Double_t D_sig=d_sig->getVal(x);
		      Double_t D_mean=d_mean->getVal(x);

		  
		      dpar_Mat[0]=dQ_da0=(((1-f)*pow(q_fac,2))/(f*S_Fun))*D_a0 + pow(q_fac,2)/(pow(f,2)*S_Fun);
		      //	      dpar_Mat[1]=dQ_da1=(((1-f)*pow(q_fac,2))/(f*S_Fun))*D_a1 + pow(q_fac,2)/(pow(f,2)*S_Fun);      
		      dpar_Mat[1]=dQ_df=pow(q_fac,2)/(pow(f,2)*S_Fun);
		      dpar_Mat[3]=dQ_dsig=((1-f)*B_Fun*pow(q_fac,2))/(f*pow(S_Fun,2))*D_sig + pow(q_fac,2)/(pow(f,2)*S_Fun);
		      dpar_Mat[2]=dQ_dmean=((1-f)*B_Fun*pow(q_fac,2))/(f*pow(S_Fun,2))*D_mean + pow(q_fac,2)/(pow(f,2)*S_Fun);

		 
		      TMatrixD aM(1,4),bM(4,1),tempM(4,1),Q_errM(1,1);
		      bM[0][0]=aM[0][0]=dQ_da0;
		      //bM[1][0]=aM[0][1]=dQ_da1;
		      bM[1][0]=aM[0][1]=dQ_df;
		      bM[3][0]=aM[0][3]=dQ_dsig;
		      bM[2][0]=aM[0][2]=dQ_dmean;
		      
		      for(Int_t m=0;m<4;m++){  // To consider only the diagonal elements 
					for(Int_t n=0;n<4;n++)
			  			{
			    			if(m!=n)
			      			cov(m,n)==0.0;
			  			}
			  		}   

		      tempM.Mult(cov,bM);
		      Q_errM.Mult(aM,tempM);
		      Q_Err=sqrt(Q_errM[0][0]);

	}
      
      

      cout << "Q and Q_Error \t " << q_fac << "\t" <<  Q_Err << endl;

      
      //************************************************************************* // This part will help you to plot
      // P L O T T I N G                                                        *
      //----------------------------------------------------------------------  */
      
            
	  model.plotOn(xframe);
	  model.plotOn(xframe,Components(bkg),LineStyle(kDashed));
	  model.plotOn(xframe,Components(signal),LineStyle(kDashed));

	  //TString *p_status = new TString(status);
	  //hist->SetTitle(*f_name+"-"+*p_status);
	  

	  	  if (i%10 == 0){
	  	  			xframe->Draw();
	  				f_canvas->Print(plotDir+*f_name+".png");
		    }
		
	
      //*************************************************************************

      //__________________________________________________________________________
      //---------------------------Error Calculation--------------------------- 
      //__________________________________________________________________________
      
	
		
	//TH1D *h1 = new TH1D("h1", "  ", 100, 1.05, 1.19);
	//TH1D *h2 = new TH1D("h2", "  ", 100, 1.05, 1.19);

	//h1 -> Fill(mmkp_fit_q, q_fac);
	//h2 -> Fill(mmkp_fit_q,(1 - q_fac));	
      { 
			chain -> GetEntry(i);


			
			mppip1_kinfit_q = evnt_mass;
		    ebeam_q = ebeam;
		    cosTheta_delta_cm_q = Xi(i,0);
		    cosTheta_Ks_hel_q = cosTheta_Ks_hel;
		    phi_Ks_hel_q = phi_Ks_hel;
			time_weight_q = time_weight;
			mKsKm_kinfit_q = mKsKm_kinfit;
			men_t_q = Xi(i,3);
			pol_q = pol;
			cosTheta_pip_gj_q = Xi(i,1);
			phi_pip_gj_q = Xi(i,2);

			
			tree->Fill();

      }
		d2.Reset();
/***************** delete unnecessary junk ********************/

		delete [] index;
		delete taux;
		delete hist;
		
		gDirectory->Delete("Aux;*");
	}

	std::cout << "WRITING TREE..." << std::endl;
	tree->Write();
	fout.Close();





}






