#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include "TH2.h"
#include "TH1.h"
#include "TVectorT.h"

#include "reader.h"
#include "bank.h"

#include "BParticle.h"
#include "BCalorimeter.h"
#include "BScintillator.h"
#include "BBand.h"
#include "BEvent.h"

#include "RCDB/Connection.h"

#include "constants.h"

using namespace std;

int getRunNumber( string filename );
double getBeamEnergy( int runNum );
void getEventInfo( BEvent eventInfo, double &integrated_charge, double &livetime, double &starttime );
void getElectronInfo( BParticle particles, int& pid, TVector3& momentum, TVector3& vertex, 
			double& time, int& charge, double& beta, double& chi2pid, int& status );
bool checkElectron( int pid, TVector3 momentum, TVector3 vertex, double time, int charge, double beta, double chi2pid, int status,
			double lV, double lW , double E_tot);
void getProtonInfo( BParticle particles, int& pid, TVector3& momentum, TVector3& vertex,
			double& time, int& charge, double& beta, double& chi2pid, int& status , int& multiplicity );
bool checkProton( int pid, TVector3 momentum, TVector3 del_vertex, double time, int charge, double beta, double chi2pid, int status, int mult );


int main(int argc, char** argv) {
	// check number of arguments
	if( argc < 3 ){
		cerr << "Incorrect number of arugments. Instead use:\n\t./code [outputFile] [inputFile]\n\n";
		return -1;
	}

	// Create output tree
	TFile * outFile = new TFile(argv[1],"RECREATE");
	TTree * outTree = new TTree("skim","CLAS and BAND Physics");
	double gated_charge	= 0;
	double eChi2pid		= 0;
	double lU		= 0;
	double lV		= 0;
	double lW		= 0;
	double E_tot		= 0;
	double p_e		= 0;
	double theta_e		= 0;
	double phi_e		= 0;
	double q		= 0;
	double theta_q		= 0;
	double phi_q		= 0;
	double nu		= 0;
	double Q2		= 0;
	double xB		= 0;
	double W2		= 0;
	double p_p		= 0;
	double theta_p		= 0;
	double phi_p		= 0;
	double p_miss		= 0;
	double m_miss		= 0;
	double theta_miss	= 0;
        double phi_miss		= 0;
	outTree->Branch("gated_charge"	,&gated_charge		);
	outTree->Branch("eChi2pid"	,&eChi2pid		);
	outTree->Branch("E_tot"		,&E_tot			);
	outTree->Branch("lU"		,&lU			);
	outTree->Branch("lV"		,&lV			);
	outTree->Branch("lW"		,&lW			);
	outTree->Branch("p_e"		,&p_e			);
	outTree->Branch("theta_e"	,&theta_e		);
	outTree->Branch("phi_e"		,&phi_e			);
	outTree->Branch("q"		,&q			);
	outTree->Branch("theta_q"	,&theta_q		);
	outTree->Branch("phi_q"		,&phi_q			);
	outTree->Branch("nu"		,&nu			);
	outTree->Branch("Q2"		,&Q2			);
	outTree->Branch("xB"		,&xB			);
	outTree->Branch("W2"		,&W2			);
	outTree->Branch("p_p"		,&p_p			);
	outTree->Branch("theta_p"	,&theta_p		);
	outTree->Branch("phi_p"		,&phi_p			);
	outTree->Branch("p_miss"	,&p_miss		);
	outTree->Branch("m_miss"	,&m_miss		);
	outTree->Branch("theta_miss"	,&theta_miss		);
	outTree->Branch("phi_miss"	,&phi_miss		);
	TH2D * h2_EoP_pe	= new TH2D("h2_EoP_pe",		"h2_EoP_pe",		400,0,10,400,0,1);
	TH1D * h1_vX_e		= new TH1D("h1_vX_e",		"h1_vX_e",		100,-10,10);
	TH1D * h1_vY_e		= new TH1D("h1_vY_e",		"h1_vY_e",		100,-10,10);
	TH2D * h2_vZ_ph_e	= new TH2D("h2_vZ_ph_e",	"h2_vZ_ph_e",		400,-180,180,200,-20,20);
	TH1D * h1_tof_e		= new TH1D("h1_tof_e",		"h1_tof_e",		300,0,30);
	TH1D * h1_chi2pid_e	= new TH1D("h1_chi2pid_e",	"h1_chi2pid_e",		100,-5,5);
	TH1D * h1_pcal_e	= new TH1D("h1_pcal_e",		"h1_pcal_e",		1000,0,10);
	TH2D * h2_EoP_lU	= new TH2D("h2_EoP_lU",		"h2_EoP_lU",		400,0,400,400,0,1);
	TH2D * h2_EoP_lV	= new TH2D("h2_EoP_lV",		"h2_EoP_lV",		400,0,400,400,0,1);
	TH2D * h2_EoP_lW	= new TH2D("h2_EoP_lW",		"h2_EoP_lW",		400,0,400,400,0,1);
	TH1D * h1_mult_p	= new TH1D("h1_mult_p",		"h1_mult_p",		10,0,10);
	TH1D * h1_vdelX_pe	= new TH1D("h1_vdelX_pe",	"h1_vdelX_pe",		200,-10,10);
        TH1D * h1_vdelY_pe	= new TH1D("h1_vdelY_pe",	"h1_vdelY_pe",		200,-10,10);
        TH2D * h2_vdelZ_pe	= new TH2D("h2_vdelZ_pe",	"h2_vdelZ_pe",		400,-180,180,300,-15,15);
	TH2D * h2_BvP_p		= new TH2D("h2_BvP_p",		"h2_BvP_p",		100,0,5,60,0,1.2);
	TH1D * h1_chi2pid_p	= new TH1D("h1_chi2pid_p",	"h1_chi2pid_p",		100,-5,5);
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

	// Connect to the RCDB
	rcdb::Connection connection("mysql://rcdb@clasdb.jlab.org/rcdb");

	gated_charge = 0;
	// Load input file
	for( int i = 2 ; i < argc ; i++ ){
		// Using run number of current file, grab the beam energy from RCDB
		int runNum = getRunNumber(argv[i]);
		auto cnd = connection.GetCondition(runNum, "beam_energy");
		double Ebeam = cnd->ToDouble() / 1000.; // [GeV]

		// Setup hipo reading for this file
		TString inputFile = argv[i];
		hipo::reader reader;
		reader.open(inputFile);
		hipo::dictionary  factory;      
		hipo::schema	  schema;
		reader.readDictionary(factory); 
		BEvent		event_info		(factory.getSchema("REC::Event"		));
		BParticle	particles		(factory.getSchema("REC::Particle"	));
		BCalorimeter	calorimeter		(factory.getSchema("REC::Calorimeter"	));
		BScintillator	scintillator		(factory.getSchema("REC::Scintillator"	));
		hipo::bank	scaler			(factory.getSchema("RUN::scaler"	));
		hipo::event 	readevent;
		
		// Loop over all events in file
		int event_counter = 0;
		while(reader.next()==true){
			eChi2pid 	= 0;
			lU		= 0;
			lV		= 0;
			lW		= 0;
			E_tot		= 0;
			p_e		= 0;
			theta_e		= 0;
			phi_e		= 0;
			q		= 0;
			theta_q		= 0;
			phi_q		= 0;
			nu		= 0;
			Q2		= 0;
			xB		= 0;
			W2		= 0;
			p_p		= 0;
			theta_p		= 0;
			phi_p		= 0;
			p_miss		= 0;
			m_miss		= 0;
			theta_miss	= 0;
			phi_miss	= 0;


			if(event_counter%10000==0) cout << "event: " << event_counter << endl;
			if(event_counter >= 1000000) break;
			event_counter++;

			// Load data structure for this event:
			reader.read(readevent);
			readevent.getStructure(event_info);
			readevent.getStructure(particles);
			readevent.getStructure(calorimeter);
			readevent.getStructure(scintillator);
			readevent.getStructure(scaler);
	
			// Currently, REC::Event has uncalibrated livetime / charge, so these will have to work
			double livetime 	= 	scaler.getFloat(2,0);
			gated_charge 		= 	scaler.getFloat(0,0) * 0.001; // [microC] -- this seems to be ~10-20% accurate

			// Get integrated charge, livetime and start-time from REC::Event
			double starttime = 0;
			if( event_info.getRows() == 0 ) continue;
			getEventInfo( event_info, gated_charge, livetime, starttime );

			// Get electron from particle bank REC::Particle
			TVector3 eVertex, eMomentum;
			int ePid = 0, eCharge = 0, eStatus = 0;
			double eTime = 0, eBeta = 0;
			getElectronInfo( particles, ePid, eMomentum, eVertex, eTime ,eCharge, eBeta, eChi2pid, eStatus );
			//	get electron information from scint and calo banks:
			double t_e 	= scintillator.getTime(0) - starttime;
			E_tot 	= calorimeter.getTotE(0);
			double E_pcal 	= calorimeter.getPcalE(0);
			lU	= calorimeter.getLU(0);
			lV	= calorimeter.getLV(0);
			lW	= calorimeter.getLW(0);
			//	Do electron PID cuts
			bool ePass = checkElectron( ePid, eMomentum, eVertex, eTime ,eCharge, eBeta, eChi2pid, eStatus , lV , lW , E_tot );
			if( !ePass ) continue;

			// From electron information and beam information, create kinematic variables
			TVector3 beamMomentum(0,0,Ebeam);
			TVector3 qMomentum; qMomentum = beamMomentum - eMomentum;
			p_e		= eMomentum.Mag();
			theta_e		= eMomentum.Theta();
			phi_e		= eMomentum.Phi();
			q		= qMomentum.Mag();
			theta_q		= qMomentum.Theta();
			phi_q		= qMomentum.Phi();
			nu		= Ebeam - p_e;
			Q2		= q*q - nu*nu;
			xB		= Q2 / (2.*mP*nu);
			W2		= mP*mP - Q2 + 2.*nu*mP;
			// Histograms to see quality of PID cuts made
			h2_EoP_pe	-> Fill( p_e , E_tot/p_e	);
			h1_vX_e		-> Fill( eVertex.X() 		);
			h1_vY_e		-> Fill( eVertex.Y() 		);
			h2_vZ_ph_e	-> Fill( phi_e * 180./M_PI, eVertex.Z() 	);
			h1_tof_e	-> Fill( t_e 			);
			h1_chi2pid_e	-> Fill( eChi2pid 		);
			h1_pcal_e	-> Fill( E_pcal 		);
			h2_EoP_lU	-> Fill( lU ,  E_tot / p_e 	);
			h2_EoP_lV	-> Fill( lV ,  E_tot / p_e 	);
			h2_EoP_lW	-> Fill( lW ,  E_tot / p_e 	);
			
			
			// Grab the proton information:
			TVector3 pVertex, pMomentum;
			int pPid = 0, pCharge = 0, pStatus = 0, pMult = 0;
			double pTime = 0, pBeta = 0, pChi2pid = 0;
			getProtonInfo( particles, pPid, pMomentum, pVertex, pTime ,pCharge, pBeta, pChi2pid, pStatus, pMult );
			//	do proton PID cuts		
			bool pPass = checkProton( pPid, pMomentum, pVertex-eVertex, pTime ,pCharge, pBeta, pChi2pid, pStatus , pMult );
			if( !pPass ) continue;
			// Create more kinematic variables
			p_p		= pMomentum.Mag();
			theta_p		= pMomentum.Theta();
			phi_p		= pMomentum.Phi();
			double E_p	= sqrt(p_p*p_p + mP*mP);
			// 	assuming d e -> e' p 
			TVector3 missMomentum; missMomentum = pMomentum - qMomentum;
			p_miss		= missMomentum.Mag();
			m_miss		= mP - mD + sqrt( pow(nu+mD-E_p,2) + pow(p_miss,2) );
			theta_miss	= missMomentum.Theta();
			phi_miss	= missMomentum.Phi();
			// Histograms to see quality of PID cuts made
			h1_mult_p	-> Fill( pMult );
			h1_vdelX_pe	-> Fill( pVertex.X() - eVertex.X() );
			h1_vdelY_pe	-> Fill( pVertex.Y() - eVertex.Y() );
			h2_vdelZ_pe	-> Fill( phi_e*180./M_PI , pVertex.Z() - eVertex.Z() );
			h2_BvP_p	-> Fill( pMomentum.Mag() , pBeta );
			h1_chi2pid_p	-> Fill( pChi2pid );


			// Histograms to compare with simulation for counting (e,e'p):
			h1_sim_p_e	-> Fill( p_e );
			h1_sim_th_e	-> Fill( theta_e*180./M_PI );
			h2_sim_th_ph_e	-> Fill( phi_e*180./M_PI , theta_e*180./M_PI );
			h1_sim_xB_e	-> Fill( xB );
			h1_sim_Q2_e	-> Fill( Q2 );
			h1_sim_W_e	-> Fill( sqrt(W2) );
			
			h1_sim_p_p	-> Fill( p_p );
			h1_sim_th_p	-> Fill( theta_p*180./M_PI );
			h2_sim_th_ph_p	-> Fill( phi_p*180./M_PI , theta_p*180./M_PI );
			h2_sim_th_th_pe	-> Fill( theta_e*180./M_PI , theta_p*180./M_PI );
			h2_sim_p_p_pe	-> Fill( p_e , p_p );
			h2_sim_ph_ph_pe	-> Fill( phi_e*180./M_PI , phi_p*180./M_PI );
			h1_sim_mass_m	-> Fill( m_miss );
			h1_sim_p_m	-> Fill( p_miss );
			h1_sim_px_m	-> Fill( p_miss*sin(theta_miss)*cos(phi_miss) );
			h1_sim_py_m	-> Fill( p_miss*sin(theta_miss)*sin(phi_miss) );
			h2_sim_p_th_m	-> Fill( theta_miss*180./M_PI , p_miss );

			// Fill tree to do any more plots on the fly
			outTree->Fill();

		} // end loop over events
	}// end loop over files
	cout << "Total charge collected in file: " << gated_charge << " [microC]\n";
	cout << "Total number of events (e,e'p): " << outTree->GetEntries() << "\n";
	

	outFile->cd();
	TVectorT<double> fileInfo(2);
	fileInfo(0) = gated_charge;
	fileInfo(1) = outTree->GetEntries();
	fileInfo.Write("fileInfo");
	h2_EoP_pe	->Write();
	h1_vX_e		->Write();
	h1_vY_e		->Write();
	h2_vZ_ph_e	->Write();
	h1_tof_e	->Write();
	h1_chi2pid_e	->Write();
	h1_pcal_e	->Write();
	h2_EoP_lU	->Write();
	h2_EoP_lV	->Write();
	h2_EoP_lW	->Write();
	h1_mult_p	->Write();
	h1_vdelX_pe	->Write();
	h1_vdelY_pe	->Write();
	h2_vdelZ_pe	->Write();
	h2_BvP_p	->Write();
	h1_chi2pid_p	->Write();
	h1_sim_p_e	->Write();
	h1_sim_th_e	->Write(); 
	h2_sim_th_ph_e	->Write(); 
	h1_sim_xB_e	->Write(); 
	h1_sim_Q2_e	->Write(); 
	h1_sim_W_e	->Write(); 
	h1_sim_p_p	->Write(); 
	h1_sim_th_p	->Write(); 
	h2_sim_th_ph_p	->Write(); 
	h2_sim_th_th_pe	->Write(); 
	h2_sim_p_p_pe	->Write(); 
	h2_sim_ph_ph_pe	->Write(); 
	h1_sim_mass_m	->Write(); 
	h1_sim_p_m	->Write(); 
	h1_sim_px_m	->Write(); 
	h1_sim_py_m	->Write(); 
	h2_sim_p_th_m	->Write();
	outTree->Write();
	outFile->Close();

	return 0;
}


int getRunNumber( string filename ){
	//string parsed = filename.substr( filename.find("inc") );
	string parsed = filename.substr( filename.find("_clas") );
	string moreparse = parsed.substr(6,8);
        return stoi(moreparse);
}

double getBeamEnergy( int runNum ){
        double thisEn = 0;

        if( runNum <= 6399 ) thisEn = 10.6;
        else{ thisEn = 10.2; }
        if( runNum == 6523 || runNum == 6524 || runNum == 6525 ) thisEn = 10.;
	

        return thisEn;
}

void getEventInfo( BEvent eventInfo, double &integrated_charge, double &livetime, double &starttime ){
	if( eventInfo.getRows() != 1 ){ 
		cerr << "getEventInfo::NotImplementedFunction\n"; 
		exit(-1); 
	}
	//integrated_charge       += (double)eventInfo.getBCG(0); 	// not calibrated currently
	//livetime 		= (double)eventInfo.getLT(0);		// not calibrated currently
	starttime		= (double)eventInfo.getSTTime(0);
	return;
}
void getElectronInfo( BParticle particles, int& pid, TVector3& momentum, TVector3& vertex,
			double& time, int& charge, double& beta, double& chi2pid, int& status ){
	pid 		= particles.getPid(0);
	momentum 	= particles.getV3P(0);
	vertex		= particles.getV3v(0);
	time		= particles.getVt(0);
	charge		= particles.getCharge(0);
	beta		= particles.getBeta(0);
	chi2pid		= particles.getChi2pid(0);
	status		= particles.getStatus(0);
	return;
}
void getProtonInfo( BParticle particles, int& pid, TVector3& momentum, TVector3& vertex,
			double& time, int& charge, double& beta, double& chi2pid, int& status , int& multiplicity ){
	// Takes the first proton in the bank
	multiplicity = 0;
	for( int row = 0 ; row < particles.getRows() ; row++ ){
		pid 		= particles.getPid(row);
		charge		= particles.getCharge(row);
		if( pid == 2212 && charge == 1 ){
			if( multiplicity == 0 ){ // if this is the first proton, save the information
				momentum 	= particles.getV3P(row);
				vertex		= particles.getV3v(row);
				time		= particles.getVt(row);
				beta		= particles.getBeta(row);
				chi2pid		= particles.getChi2pid(row);
				status		= particles.getStatus(row);
			}
			multiplicity ++;
		}
	}
}
bool checkElectron( int pid, TVector3 momentum, TVector3 vertex, double time, int charge, double beta, double chi2pid, int status,
			double lV, double lW , double E_tot){
			// Anymore fiducial cuts that are needed for the electron:
			// 	E/p cut
			// 	vertex cut
			// 	minimum momentum  cut / maximum  momentum  cut
			// 	electron ToF cut
			// 	minimum W cut
			//	chi2 cut
			//	pcal energy cut?
	if( pid != 11 || charge != -1 ) return false;
	if( vertex.X() < -2 || vertex.X() > 2) return false;
	if( vertex.Y() < -2 || vertex.Y() > 2) return false;
	if( vertex.Z() < -7 || vertex.Z() > 2) return false;
	if( time < 15 ) return false;
	if( lV < 15 || lW < 15 ) return false;
	if( chi2pid < -3 || chi2pid > 2 ) return false;
	if( momentum.Mag() < 2 || momentum.Mag() > 10.6 ) return false;
	if( E_tot / momentum.Mag() < 0.15 || E_tot / momentum.Mag() > 0.3 ) return false;

	return true;
}
bool checkProton( int pid, TVector3 momentum, TVector3 del_vertex, double time, int charge, double beta, double chi2pid, int status, int mult ){
	if( momentum.Mag() == 0 || momentum.Mag() > 10.6 ) return false;
	if( beta <= 0 || beta > 1 ) return false;
	if( mult != 1 ) return false;
	if( del_vertex.Z() < -5 || del_vertex.Z() > 5 ) return false;
	if( chi2pid < -3 || chi2pid > 3 ) return false;
	
	
	return true;
}
