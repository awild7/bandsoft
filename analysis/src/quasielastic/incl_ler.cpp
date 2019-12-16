#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>

#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include "TH1.h"
#include "TH2.h"

#include "reader.h"
#include "bank.h"

#include "BParticle.h"
#include "BCalorimeter.h"
#include "BScintillator.h"
#include "BBand.h"
#include "BEvent.h"

#include "RCDB/Connection.h"

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
	TH2D * phi_theta = new TH2D("phi_v_theta","phiE_v_thetaE;phi;theta;Counts",100,0,2*M_PI,100,0.1,0.7);
	TH2D * eop_p = new TH2D("eop_v_p","eE/pE_v_pE;eE/pE;pE;Counts",110,0,11,150,0,1.5);

	vector<TH1*> hist_list_1D;
	TH1D * hist_xB =  new TH1D("hist_xB_Incl" ,"hist;xB;Counts",50,0,1.5);
	hist_list_1D.push_back(hist_xB);
	TH1D * hist_Q2 =  new TH1D("hist_Q2_Incl" ,"hist;Q2;Counts",100,0,5);
	hist_list_1D.push_back(hist_Q2);
	TH1D * hist_pE =  new TH1D("hist_pE_Incl" ,"hist;pE;Counts",110,0,11);
	hist_list_1D.push_back(hist_pE);
	TH1D * hist_W2 =  new TH1D("hist_W2_Incl" ,"hist;W2;Counts",100,0,8);
	hist_list_1D.push_back(hist_W2);
	TH1D * hist_chi =  new TH1D("hist_chi_Incl" ,"hist;chi;Counts",50,0,1.5);
	hist_list_1D.push_back(hist_chi);

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
		BCalorimeter    calo                    (factory.getSchema("REC::Calorimeter"   ));
		hipo::event 	readevent;
		
		// Loop over all events in file
		int event_counter = 0;
		double int_charge = 0;
		while(reader.next()==true){
			if(event_counter%10000==0) cout << "event: " << event_counter << endl;
			event_counter++;

			// Load data structure for this event:
			reader.read(readevent);
			readevent.getStructure(event_info);
			readevent.getStructure(particles);
			readevent.getStructure(calo);
			
			// Get integrated charge, livetime and start-time from REC::Event
			double livetime = 0 , starttime = 0;
			if( event_info.getRows() == 0 ) continue;
			getEventInfo( event_info, int_charge, livetime, starttime );

			// Get electron from particle bank REC::Particle
			TVector3 eVertex, eMomentum;
			int ePid = 0, eCharge = 0, eStatus = 0;
			double eTime = 0, eBeta = 0, eChi2pid = 0;
			getElectronInfo( particles, ePid, eMomentum, eVertex, eTime ,eCharge, eBeta, eChi2pid, eStatus );
			//	do electron cuts
			bool ePass = checkElectron( ePid, eMomentum, eVertex, eTime ,eCharge, eBeta, eChi2pid, eStatus );
			if( !ePass ) continue;
			double Ee = calo.getTotE (0);
			// Anymore fiducial cuts that are needed

			TVector3 beamVec(0,0,fixed_Ebeam);
			TVector3 qVec = eMomentum - beamVec;
			// From electron information and beam information, create kinematic variables
			double p_e 		= eMomentum.Mag();
			double theta_e 	        = eMomentum.Theta();
			double phi_e		= eMomentum.Phi();
			double q		= qVec.Mag();
			double theta_q		= qVec.Theta();
			double phi_q		= qVec.Phi();
			double nu 		= fixed_Ebeam - sqrt( p_e*p_e + mE*mE );
			double Q2 		= q*q - nu*nu;
			double xB 		= Q2 / (2.*mP*nu);
			double W2     		= mP*mP - Q2 + 2*nu*mP;	
			double EoP		= Ee / p_e;
			//			double sector_e	        = calo.getSector(0);	// technically this is a bit wrong, but it should all have same sector...
			double vrt_x_e		= eVertex.X();
			double vrt_y_e		= eVertex.Y();
			double vrt_z_e		= eVertex.Z();
			//			double t_e       	= scintillator.getTime(0);
			
			phi_theta->Fill(phi_e,theta_e);
			eop_p->Fill(p_e,EoP);
			hist_xB->Fill(xB);
			hist_Q2->Fill(Q2);
			hist_pE->Fill(p_e);
			hist_W2->Fill(W2);
			hist_chi->Fill(eChi2pid);


		} // end loop over events
		cout << int_charge << "\n";
	}// end loop over files
	

	outFile->cd();
	outTree->Write();
	for(int k=0; k<hist_list_1D.size(); k++){
	  hist_list_1D[k]->Write();
	}
	phi_theta->Write();
	eop_p->Write();			 
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
	integrated_charge       += (double)eventInfo.getBCG(0);
	livetime 		= (double)eventInfo.getLT(0);
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
bool checkElectron( int pid, TVector3 momentum, TVector3 vertex, double time, int charge, double beta, double chi2pid, int status ){
	if( pid != 11 || charge != -1 ) return false;
	return true;
}
