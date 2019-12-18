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

	vector<TH1*> hist_list_2D;	

	TH2D * eop_U = new TH2D("eop_v_U","eE/pE_v_U;eE/pE;U;Counts",40,0,300,100,0,0.5);
	hist_list_2D.push_back(eop_U);
	TH2D * eop_V = new TH2D("eop_v_V","eE/pE_v_V;eE/pE;V;Counts",40,0,300,100,0,0.5);
	hist_list_2D.push_back(eop_V);
	TH2D * eop_W = new TH2D("eop_v_W","eE/pE_v_W;eE/pE;W;Counts",40,0,300,100,0,0.5);
	hist_list_2D.push_back(eop_W);
	TH2D * chi_U = new TH2D("chi_v_U","chi_v_U;chi;U;Counts",40,0,300,140,-7,7);
	hist_list_2D.push_back(chi_U);
	TH2D * chi_V = new TH2D("chi_v_V","chi_v_V;chi;V;Counts",40,0,300,140,-7,7);
	hist_list_2D.push_back(chi_V);
	TH2D * chi_W = new TH2D("chi_v_W","chi_v_W;chi;W;Counts",40,0,300,140,-7,7);
	hist_list_2D.push_back(chi_W);

	vector<TH1*> hist_list_1D;
	int lengA = 6;
	double ULower[] = {0,10,20,30,40,50};
	double VWLower[] = {0,5,10,15,20,25};
	TH1D * hist_chi[lengA][lengA];
	TH1D * hist_eop[lengA][lengA];
	for(int k=0; k<lengA; k++){
	  for(int l=0; l<lengA; l++){
	  char temp1[100];
	  sprintf(temp1,"hist_chi_%f_%f",ULower[k],VWLower[l]);
	  hist_chi[k][l] = new TH1D(temp1,"hist;chi;Counts",140,-7,7);
	  hist_list_1D.push_back(hist_chi[k][l]);

	  char temp2[100];
	  sprintf(temp2,"hist_eop_%f_%f",ULower[k],VWLower[l]);
	  hist_eop[k][l] = new TH1D(temp2,"hist;eop;Counts",100,0,0.5);
	  hist_list_1D.push_back(hist_eop[k][l]);
	  }
	}


	
	TH1D * hist_xB =  new TH1D("hist_xB_Incl" ,"hist;xB;Counts",50,0,1.5);
	hist_list_1D.push_back(hist_xB);
	TH1D * hist_Q2 =  new TH1D("hist_Q2_Incl" ,"hist;Q2;Counts",100,0,8);
	hist_list_1D.push_back(hist_Q2);
	TH1D * hist_pE =  new TH1D("hist_pE_Incl" ,"hist;pE;Counts",110,0,11);
	hist_list_1D.push_back(hist_pE);
	TH1D * hist_W2 =  new TH1D("hist_W2_Incl" ,"hist;W2;Counts",100,0,30);
	hist_list_1D.push_back(hist_W2);
	
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
		BScintillator	scintillator		(factory.getSchema("REC::Scintillator"	));
		hipo::bank	scaler			(factory.getSchema("RUN::scaler"	));
		hipo::event 	readevent;
		
		// Loop over all events in file
		int event_counter = 0;
		double int_charge = 0;
		while( (reader.next()==true) && (event_counter<10000) ){
			if(event_counter%10000==0) cout << "event: " << event_counter << endl;
			event_counter++;

			// Load data structure for this event:
			reader.read(readevent);
			readevent.getStructure(event_info);
			readevent.getStructure(particles);
			readevent.getStructure(calo);
			readevent.getStructure(scintillator);
			readevent.getStructure(scaler);
			
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
			double t_e 	= scintillator.getTime(0) - starttime;
			double E_tot 	= calo.getTotE(0);
			double E_pcal 	= calo.getPcalE(0);
			double lU	= calo.getLU(0);
			double lV	= calo.getLV(0);
			double lW	= calo.getLW(0);
			// Anymore fiducial cuts that are needed

			TVector3 beamVec(0,0,Ebeam);
			TVector3 qVec = beamVec - eMomentum;
			// From electron information and beam information, create kinematic variables
			double p_e 		= eMomentum.Mag();
			double theta_e 	        = eMomentum.Theta();
			double phi_e		= eMomentum.Phi();
			double q		= qVec.Mag();
			double theta_q		= qVec.Theta();
			double phi_q		= qVec.Phi();
			double nu 		= Ebeam - sqrt( p_e*p_e + mE*mE );
			double Q2 		= q*q - nu*nu;
			double xB 		= Q2 / (2.*mP*nu);
			double W2     		= mP*mP - Q2 + 2*nu*mP;	
			double EoP		= E_tot / p_e;
			//			double sector_e	        = calo.getSector(0);	// technically this is a bit wrong, but it should all have same sector...
			double vrt_x_e		= eVertex.X();
			double vrt_y_e		= eVertex.Y();
			double vrt_z_e		= eVertex.Z();
			//			double t_e       	= scintillator.getTime(0);
			
			eop_U        	->Fill ( lU    		,EoP      	);
			eop_V        	->Fill ( lV    		,EoP      	);
			eop_W        	->Fill ( lW    		,EoP      	);
			chi_U        	->Fill ( lU    		,eChi2pid      	);
			chi_V        	->Fill ( lV    		,eChi2pid      	);
			chi_W        	->Fill ( lW    		,eChi2pid      	);

			for(int k=0; k<lengA; k++){
			  for(int l=0; l<lengA; l++){
			    if(lU > ULower[k]){
			      if( (lV > VWLower[l]) && (lW > VWLower[l]) ){
				hist_eop[k][l]->Fill( EoP );
				hist_chi[k][l]->Fill( eChi2pid );
			      } 
			    } 
			  }
			}

			
			hist_xB      	->Fill ( xB           	);
			hist_Q2      	->Fill ( Q2           	);
			hist_pE      	->Fill ( p_e          	);
			hist_W2      	->Fill ( W2           	);

		} // end loop over events
		cout << int_charge << "\n";
	}// end loop over files
	

	outFile->cd();
	outTree->Write();
	for(int k=0; k<hist_list_2D.size(); k++){
	  hist_list_2D[k]->Write();
	}	
	for(int k=0; k<hist_list_1D.size(); k++){
	  hist_list_1D[k]->Write();
	}
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

	/*TH2D * phi_theta = new TH2D("phi_v_theta","phiE_v_thetaE;phi;theta;Counts",100,-M_PI,M_PI,100,0.1,0.7);
	hist_list_2D.push_back(phi_theta);
	TH2D * eop_p = new TH2D("eop_v_p","eE/pE_v_pE;eE/pE;pE;Counts",110,0,11,100,0,0.5);
	hist_list_2D.push_back(eop_p);
	TH2D * eop_W2 = new TH2D("eop_v_W2","eE/pE_v_W2;eE/pE;W2;Counts",100,0,30,100,0,0.5);
	hist_list_2D.push_back(eop_W2);
	TH2D * eop_xB = new TH2D("eop_v_xB","eE/pE_v_xB;eE/pE;xB;Counts",50,0,1.5,100,0,0.5);
	hist_list_2D.push_back(eop_xB);
	TH2D * eop_vtx_X = new TH2D("eop_v_vtx_X","eE/pE_v_vtx_X;eE/pE;vtx_X;Counts",40,-10,10,100,0,0.5);
	hist_list_2D.push_back(eop_vtx_X);
	TH2D * eop_vtx_Y = new TH2D("eop_v_vtx_Y","eE/pE_v_vtx_Y;eE/pE;vtx_Y;Counts",40,-10,10,100,0,0.5);
	hist_list_2D.push_back(eop_vtx_Y);
	TH2D * eop_vtx_Z = new TH2D("eop_v_vtx_Z","eE/pE_v_vtx_Z;eE/pE;vtx_Z;Counts",35,-20,15,100,0,0.5);
	hist_list_2D.push_back(eop_vtx_Z);
	TH2D * chi_p = new TH2D("chi_v_p","chi_v_pE;chi;pE;Counts",110,0,11,140,-7,7);
	hist_list_2D.push_back(chi_p);
	TH2D * chi_W2 = new TH2D("chi_v_W2","chi_v_W2;chi;W2;Counts",100,0,30,140,-7,7);
	hist_list_2D.push_back(chi_W2);
	TH2D * chi_xB = new TH2D("chi_v_xB","chi_v_xB;chi;xB;Counts",50,0,1.5,140,-7,7);
	hist_list_2D.push_back(chi_xB);
	TH2D * chi_vtx_X = new TH2D("chi_v_vtx_X","chi_v_vtx_X;chi;vtx_X;Counts",40,-10,10,140,-7,7);
	hist_list_2D.push_back(chi_vtx_X);
	TH2D * chi_vtx_Y = new TH2D("chi_v_vtx_Y","chi_v_vtx_Y;chi;vtx_Y;Counts",40,-10,10,140,-7,7);
	hist_list_2D.push_back(chi_vtx_Y);
	TH2D * chi_vtx_Z = new TH2D("chi_v_vtx_Z","chi_v_vtx_Z;chi;vtx_Z;Counts",35,-20,15,140,-7,7);
	hist_list_2D.push_back(chi_vtx_Z);
	*/

	/*
	TH1D * hist_U =  new TH1D("hist_U_Incl" ,"hist;U;Counts",40,0,300);
	hist_list_1D.push_back(hist_U);
	TH1D * hist_V =  new TH1D("hist_V_Incl" ,"hist;V;Counts",40,0,300);
	hist_list_1D.push_back(hist_V);
	TH1D * hist_W =  new TH1D("hist_W_Incl" ,"hist;W;Counts",40,0,300);
	hist_list_1D.push_back(hist_W);

	TH1D * hist_vtx_X =  new TH1D("hist_vtx_X_Incl" ,"hist;vtx_X;Counts",40,-10,10);
	hist_list_1D.push_back(hist_vtx_X);
	TH1D * hist_vtx_Y =  new TH1D("hist_vtx_Y_Incl" ,"hist;vtx_Y;Counts",40,-10,10);
	hist_list_1D.push_back(hist_vtx_Y);
	TH1D * hist_vtx_Z =  new TH1D("hist_vtx_Z_Incl" ,"hist;vtx_Z;Counts",35,-20,15);
	hist_list_1D.push_back(hist_vtx_Z);
	*/

/*
			phi_theta    	->Fill ( phi_e		,theta_e	);
			eop_p        	->Fill ( p_e	       	,EoP      	);
			eop_W2        	->Fill ( W2    		,EoP      	);
			eop_xB        	->Fill ( xB    		,EoP      	);
			eop_vtx_X      	->Fill ( vrt_x_e       	,EoP      	);
			eop_vtx_Y     	->Fill ( vrt_y_e       	,EoP      	);
			eop_vtx_Z      	->Fill ( vrt_z_e       	,EoP      	);
			chi_p        	->Fill ( p_e	       	,eChi2pid      	);
			chi_W2        	->Fill ( W2    		,eChi2pid      	);
			chi_xB        	->Fill ( xB    		,eChi2pid      	);
			chi_vtx_X      	->Fill ( vrt_x_e       	,eChi2pid      	);
			chi_vtx_Y     	->Fill ( vrt_y_e       	,eChi2pid      	);
			chi_vtx_Z      	->Fill ( vrt_z_e       	,eChi2pid      	);
*/
			/*
			hist_U       	->Fill ( lU           	);
			hist_V       	->Fill ( lV           	);
			hist_W       	->Fill ( lW           	);
			hist_vtx_X     	->Fill ( vrt_x_e       	);
			hist_vtx_Y     	->Fill ( vrt_y_e       	);
			hist_vtx_Z     	->Fill ( vrt_z_e       	);
			*/
	/*
	for(int k=0; k<lengA; k++){
	  for(int l=0; l<lengA; l++){
	    if(lU > ULower[k]){
	      if( (lV > VWLower[l]) && (lW > VWLower[l]) ){
		TF1 * myChiFit = new TF1("myChiFit","gaus",-7,7);   //Define Fit
		myChiFit->SetParameter(0,hist_chi[k][l]->GetMaximum());
		myADCfit->SetParameter(1,hist_chi[k][l]->GetMean());
		myADCfit->SetParameter(2,hist_chi[k][l]-GetStd());
		TFitResultPtr pointChi = hist_chi[k][l]->Fit(myChiFit,"qeSrn","",-7,7);
		cout<<"U Cut: "<<ULower[k]<<"   V & W Cut: "<<VWLower[k]
		    <<"Fit Mean: "<<pointChi->Parameter(1)<<"   Fit Mean Error: "<<pointChi->ParError(1)
		    <<"Fit Sigma: "<<pointChi->Parameter(2)<<"   Fit Sigma Error: "<<pointChi->ParError(2);

	      } 
	    } 
	  }
	}
	*/
