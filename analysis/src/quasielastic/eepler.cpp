#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"

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
bool checkElectron( int pid, TVector3 momentum, TVector3 vertex, double time, int charge, double beta, double chi2pid, int status );

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
	outTree->Branch("gated_charge"	,&gated_charge		);
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
	TH2D * h2_EoP_pe	= new TH2D("h2_EoP_pe",		"h2_EoP_pe",		1000,0,10,10);
	TH1D * h1_vX_e		= new TH1D("h1_vX_e",		"h1_vX_e",		);
	TH1D * h1_vY_e		= new TH1D("h1_vY_e",		"h1_vY_e",		);
	TH1D * h1_vZ_e		= new TH1D("h1_vZ_e",		"h1_vZ_e",		);
	TH1D * h1_tof_e		= new TH1D("h1_tof_e",		"h1_tof_e",		);
	TH1D * h1_chi2pid_e	= new TH1D("h1_chi2pid_e",	"h1_chi2pid_e",		);
	TH1D * h1_pcal_e	= new TH1D("h1_pcal_e",		"h1_pcal_e",		);
	TH2D * h2_EoP_lU	= new TH2D("h2_EoP_lU",		"h2_EoP_lU",		);
	TH2D * h2_EoP_lV	= new TH2D("h2_EoP_lV",		"h2_EoP_lV",		);
	TH2D * h2_EoP_lW	= new TH2D("h2_EoP_lW",		"h2_EoP_lW",		);




	// Connect to the RCDB
	rcdb::Connection connection("mysql://rcdb@clasdb.jlab.org/rcdb");

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
		gated_charge = 0;
		while(reader.next()==true){
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
			if(event_counter%10000==0) cout << "event: " << event_counter << endl;
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
			double eTime = 0, eBeta = 0, eChi2pid = 0;
			getElectronInfo( particles, ePid, eMomentum, eVertex, eTime ,eCharge, eBeta, eChi2pid, eStatus );
			//	do basic electron cuts
			bool ePass = checkElectron( ePid, eMomentum, eVertex, eTime ,eCharge, eBeta, eChi2pid, eStatus );
			if( !ePass ) continue;
			//	get electron information from scint and calo banks:
			double t_e 	= scintillator.getTime(0) - starttime;
			double E_tot 	= calorimeter.getTotE(0);
			double E_pcal 	= calorimeter.getPcalE(0);
			double lU	= calorimeter.getLU(0);
			double lV	= calorimeter.getLV(0);
			double lW	= calorimeter.getLW(0);

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
			W2		= mP*mP - Q2 - 2.*nu*mP;

			// Anymore fiducial cuts that are needed for the electron:
			// 	E/p cut
			// 	vertex cut
			// 	minimum momentum  cut / maximum  momentum  cut
			// 	electron ToF cut
			// 	minimum W cut
			//	chi2 cut
			//	pcal energy cut?
			//	uvw cut?
			h2_EoP_pe	-> Fill( E_tot / p_e , p_e 	);
			h1_vX_e		-> Fill( eVertex.X() 		);
			h1_vY_e		-> Fill( eVertex.Y() 		);
			h1_vZ_e		-> Fill( eVertex.Z() 		);
			h1_tof_e	-> Fill( t_e 			);
			h1_chi2pid_e	-> Fill( eChi2pid 		);
			h1_pcal_e	-> Fill( E_pcal 		);
			h2_EoP_lU	-> Fill( E_tot / p_e , lU 	);
			h2_EoP_lV	-> Fill( E_tot / p_e , lV 	);
			h2_EoP_lW	-> Fill( E_tot / p_e , lW 	);

			
			// Grab the proton information:
			TVector3 pVertex, pMomentum;
			int pPid = 0, pCharge = 0, pStatus = 0;
			double pTime = 0, pBeta = 0, pChi2pid = 0;
			getProtonInfo( particles, pPid, pMomentum, pVertex, pTime ,pCharge, pBeta, pChi2pid, pStatus );

			// Proton PID cuts:
			//	beta vs p cut
			//	vertex cut on diff of vertex of electron & proton
			//	missing mass cut
			//	missing momentum
			//	# protons cut
			//	pe vs pp
			//	theta p vs p_p
			//	phi_p vs phi_e
			//	deltaT vs p cut

	
			// Fill tree if cuts pass to give us counts:
			outTree->Fill();

		} // end loop over events
		cout << "Total charge collected in file: " << gated_charge << " [microC]\n";
	}// end loop over files
	

	outFile->cd();
	h2_EoP_pe	->Write();
	h1_vX_e		->Write();
	h1_vY_e		->Write();
	h1_vZ_e		->Write();
	h1_tof_e	->Write();
	h1_chi2pid_e	->Write();
	h1_pcal_e	->Write();
	h2_EoP_lU	->Write();
	h2_EoP_lV	->Write();
	h2_EoP_lW	->Write();
	outTree->Write();
	outFile->Close();

	return 0;
}


int getRunNumber( string filename ){
	string parsed = filename.substr( filename.find("inc") );
	string moreparse = parsed.substr(4,6);
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
bool getProtonInfo( BParticle particles, int& pid, TVector3& momentum, TVector3& vertex,
			double& time, int& charge, double& beta, double& chi2pid, int& status ){
	// Takes the first proton in the bank
	for( int row = 0 ; row < particles->getRows() ; row++ ){
		pid 		= particles.getPid(row);
		momentum 	= particles.getV3P(row);
		vertex		= particles.getV3v(row);
		time		= particles.getVt(row);
		charge		= particles.getCharge(row);
		beta		= particles.getBeta(row);
		chi2pid		= particles.getChi2pid(row);
		status		= particles.getStatus(row);
			
		if( pid == 2212 && charge == 1 ) return true;
	}
	return false;
}
bool checkElectron( int pid, TVector3 momentum, TVector3 vertex, double time, int charge, double beta, double chi2pid, int status ){
	if( pid != 11 || charge != -1 ) return false;
	return true;
}
