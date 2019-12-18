
#include <iostream>
#include <cmath>
#include <cstdlib>
#include <unistd.h>
#include <ctype.h>

#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include "TCut.h"

#include "constants.h"

using namespace std;
int main(int argc, char ** argv){

	if (argc<3){
		cerr << "Wrong number of arguments used.\n\tPlease instead use: ./code [inputCalcFile] [OutputHistsFile]\n";
		return -1;
	}

	// Setup input tree
	TFile * inFile = new TFile(argv[1]);
	TTree * inTree = (TTree*)inFile->Get("skim");
	double p_e	=0;	
	double theta_e	=0;
	double phi_e	=0;	
	double p_p	=0;	
	double theta_p	=0;
	double phi_p	=0;	
	double p_n	=0;	
	double theta_n	=0;
	double phi_n	=0;	
	double q 	=0;	
	double theta_q 	=0;
	double phi_q 	=0;	
	double nu 	=0;	
	double Q2 	=0;	
	double xB 	=0;	
	double W2 	=0;	
	double phi_nq 	=0;
	double theta_nq  =0;	
	double CosTheta_nq=0; 	
	double E_n 	=0;	
	double Wp 	=0;	
	double Xp 	=0;	
	double As 	=0;	
	double weighting =0;
	int e_tagged 	=0;	
	int n_tagged 	=0;	
	inTree->Branch("p_e"		,&p_e		);
        inTree->Branch("theta_e"	,&theta_e	);
        inTree->Branch("phi_e"		,&phi_e		);
        inTree->Branch("p_p"		,&p_p		);
        inTree->Branch("theta_p"	,&theta_p	);
        inTree->Branch("phi_p"		,&phi_p		);
        inTree->Branch("p_n"		,&p_n		);
        inTree->Branch("theta_n"	,&theta_n	);
        inTree->Branch("phi_n"		,&phi_n		);
        inTree->Branch("q" 		,&q 		);
        inTree->Branch("theta_q" 	,&theta_q 	);
        inTree->Branch("phi_q" 		,&phi_q 	);
        inTree->Branch("nu"		,&nu 		);
        inTree->Branch("Q2"		,&Q2 		);
        inTree->Branch("xB"		,&xB 		);
        inTree->Branch("W2"		,&W2 		);
        inTree->Branch("phi_nq" 	,&phi_nq 	);
        inTree->Branch("theta_nq" 	,&theta_nq 	);
        inTree->Branch("CosTheta_nq" 	,&CosTheta_nq 	);
        inTree->Branch("E_n" 		,&E_n 		);
        inTree->Branch("Wp" 		,&Wp 		);
        inTree->Branch("Xp" 		,&Xp 		);
	inTree->Branch("As" 		,&As 		);	
        inTree->Branch("weighting"	,&weighting	);
	inTree->Branch("e_tagged"	,&e_tagged	);
	inTree->Branch("n_tagged"	,&n_tagged	);

	// Define histograms for the output
	TFile * outFile = new TFile(argv[2],"RECREATE");
	TH1D * h1_sim_p_e	= new TH1D("h1_sim_p_e",	"h1_sim_p_e",		20,2.5,4.5);
	TH1D * h1_sim_th_e	= new TH1D("h1_sim_th_e",	"h1_sim_th_e",		25,0,25);
	TH2D * h2_sim_th_ph_e	= new TH2D("h2_sim_th_ph_e",	"h2_sim_th_ph_e",	108,-180,180,60,0,30);
	TH1D * h1_sim_xB_e	= new TH1D("h1_sim_xB_e",	"h1_sim_xB_e",		10,0.5,1.5);
	TH1D * h1_sim_Q2_e	= new TH1D("h1_sim_Q2_e",	"h1_sim_Q2_e",		15,0,3);
	TH1D * h1_sim_W_e	= new TH1D("h1_sim_W_e",	"h1_sim_W_e",		20,0,2);
	TH1D * h1_sim_p_p	= new TH1D("h1_sim_p_p",	"h1_sim_p_p",		20,0,4);
	TH1D * h1_sim_th_p	= new TH1D("h1_sim_th_p",	"h1_sim_th_p",		20,20,80);
	TH2D * h2_sim_th_ph_p	= new TH2D("h2_sim_th_ph_p",	"h2_sim_th_ph_p",	36,-180,180,20,20,80);
	TH2D * h2_sim_th_th_pe	= new TH2D("h2_sim_th_th_pe",	"h2_sim_th_th_pe",	60,0,30,20,20,80);
	TH2D * h2_sim_p_p_pe	= new TH2D("h2_sim_p_p_pe",	"h2_sim_p_p_pe",	20,2.5,4.5,20,0,4);
	TH2D * h2_sim_ph_ph_pe	= new TH2D("h2_sim_ph_ph_pe",	"h2_sim_ph_ph_pe",	108,-180,180,36,-180,180);
	TH1D * h1_sim_mass_m	= new TH1D("h1_sim_mass_m",	"h1_sim_mass_m",	10,0,0.1);
	TH1D * h1_sim_p_m	= new TH1D("h1_sim_p_m",	"h1_sim_p_m",		6,0.2,0.8);
	TH1D * h1_sim_px_m	= new TH1D("h1_sim_px_m",	"h1_sim_px_m",		10,-0.5,0.5);
	TH1D * h1_sim_py_m	= new TH1D("h1_sim_py_m",	"h1_sim_py_m",		10,-0.5,0.5);
	TH2D * h2_sim_p_th_m	= new TH2D("h2_sim_p_th_m",	"h2_sim_p_th_m",	36,0,180,6,0.2,0.8);

	// Cuts used for the plots
	TCut inclusive = "weighting*(e_tagged == 1)";
	TCut tagged = "weighting*(e_tagged == 1 && n_tagged == 1)";

	inTree->Draw("p_e >> h1_sim_p_e"							,inclusive);
	inTree->Draw("theta_e*180./TMath::Pi() >> h1_sim_th_e"					,inclusive);
	inTree->Draw("theta_e*180./TMath::Pi() : phi_e*180./TMath::Pi() - 180>> h2_sim_th_ph_e"	,inclusive);
	inTree->Draw("xB >> h1_sim_xB_e"							,inclusive);
	inTree->Draw("Q2 >> h1_sim_Q2_e"							,inclusive);
	inTree->Draw("sqrt(W2) >> h1_sim_W_e"							,inclusive);
	inTree->Draw("p_p >> h1_sim_p_p"							,inclusive);
	inTree->Draw("theta_p*180./TMath::Pi() >> h1_sim_th_p"					,inclusive);
	inTree->Draw("theta_p*180./TMath::Pi() : phi_p*180./TMath::Pi() - 180>> h2_sim_th_ph_p"	,inclusive);
	inTree->Draw("theta_p*180./TMath::Pi() : theta_e*180./TMath::Pi()>> h2_sim_th_th_pe"	,inclusive);
	inTree->Draw("p_p: p_e >> h2_sim_p_p_pe"						,inclusive);
	inTree->Draw("phi_p*180./TMath::Pi() -180 : phi_e*180./TMath::Pi() - 180 >> h2_sim_ph_ph_pe"	,inclusive);
	inTree->Draw("m_miss >> h1_sim_mass_m"							,inclusive);
	inTree->Draw("p_miss >> h1_sim_p_m"							,inclusive);
	inTree->Draw("p_miss*sin(theta_miss)*cos(phi_miss) >> h1_sim_px_m"			,inclusive);
	inTree->Draw("p_miss*sin(theta_miss)*sin(phi_miss) >> h1_sim_py_m"			,inclusive);
	inTree->Draw("p_miss : theta_miss*180./TMath::Pi() >> h2_sim_p_th_m"			,inclusive);

	outFile->cd();
	h1_sim_p_e	-> Write();
	h1_sim_th_e	-> Write();
	h2_sim_th_ph_e	-> Write();
	h1_sim_xB_e	-> Write();
	h1_sim_Q2_e	-> Write();
	h1_sim_W_e	-> Write();
	h1_sim_p_p	-> Write();
	h1_sim_th_p	-> Write();
	h2_sim_th_ph_p	-> Write();
	h2_sim_th_th_pe	-> Write();
	h2_sim_p_p_pe	-> Write();
	h2_sim_ph_ph_pe	-> Write();
	h1_sim_mass_m	-> Write();
	h1_sim_p_m	-> Write();
	h1_sim_px_m	-> Write();
	h1_sim_py_m	-> Write();
	h2_sim_p_th_m	-> Write();
	outFile->Close();

	return 1;
}
