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
#include "TCanvas.h"
#include "TPaveText.h"
#include "TVectorT.h"

#include "constants.h"

using namespace std;
void TH2D_Style( TH2D * hist, TString xLabel, TString yLabel, int line_color );
void TH1D_Style( TH1D * hist, TString xLabel, TString yLabel, int line_color );
TCanvas * print1Hists( TH1D * dat , TH1D * sim , int idx , TCanvas * saveon );
TCanvas * print2Hists( TH2D * dat , TH2D * sim , int idx , TCanvas * saveon );
int main(int argc, char ** argv){

	if (argc<3){
		cerr << "Wrong number of arguments used.\n\tPlease instead use: ./code [SimulationHists] [DataHists]\n";
		return -1;
	}

	TFile *inFileSim 	= new TFile(argv[1]);
	TFile *inFileData 	= new TFile(argv[2]);
	

	cout<<"Loading all the simulation histograms.\n";

	TH1D * sim_h1_sim_p_e		= (TH1D*)inFileSim->Get("h1_sim_p_e");
	TH1D * sim_h1_sim_th_e		= (TH1D*)inFileSim->Get("h1_sim_th_e");	
	TH2D * sim_h2_sim_th_ph_e	= (TH2D*)inFileSim->Get("h2_sim_th_ph_e");
	TH1D * sim_h1_sim_xB_e		= (TH1D*)inFileSim->Get("h1_sim_xB_e");
	TH1D * sim_h1_sim_Q2_e		= (TH1D*)inFileSim->Get("h1_sim_Q2_e");	
	TH1D * sim_h1_sim_W_e		= (TH1D*)inFileSim->Get("h1_sim_W_e");	
	TH1D_Style( sim_h1_sim_p_e		, "Electron Momentum [GeV]"	, "Counts [a.u.]"		, 2 );
	TH1D_Style( sim_h1_sim_th_e		, "Electron Theta [deg.]"	, "Counts [a.u.]"       	, 2 );
	TH2D_Style( sim_h2_sim_th_ph_e		, "Electron Phi [deg.]"		, "Electron Theta [deg.]"	, 2 );
	TH1D_Style( sim_h1_sim_xB_e		, "Bjorken-x"    		, "Counts [a.u.]"       	, 2 );
	TH1D_Style( sim_h1_sim_Q2_e		, "Q2 [GeV2]" 		    	, "Counts [a.u.]"       	, 2 );
	TH1D_Style( sim_h1_sim_W_e		, "W [GeV]"     		, "Counts [a.u.]"       	, 2 );

	cout<<"Loading all the data histograms.\n";

	TH1D * dat_h1_sim_p_e		= (TH1D*)inFileData->Get("h1_sim_p_e");
	TH1D * dat_h1_sim_th_e		= (TH1D*)inFileData->Get("h1_sim_th_e");	
	TH2D * dat_h2_sim_th_ph_e	= (TH2D*)inFileData->Get("h2_sim_th_ph_e");
	TH1D * dat_h1_sim_xB_e		= (TH1D*)inFileData->Get("h1_sim_xB_e");
	TH1D * dat_h1_sim_Q2_e		= (TH1D*)inFileData->Get("h1_sim_Q2_e");	
	TH1D * dat_h1_sim_W_e		= (TH1D*)inFileData->Get("h1_sim_W_e");	
	TH1D_Style( dat_h1_sim_p_e		, "Electron Momentum [GeV]"	, "Counts [a.u.]"		, 8 );
	TH1D_Style( dat_h1_sim_th_e		, "Electron Theta [deg.]"	, "Counts [a.u.]"       	, 8 );
	TH2D_Style( dat_h2_sim_th_ph_e		, "Electron Phi [deg.]"		, "Electron Theta [deg.]"	, 8 );
	TH1D_Style( dat_h1_sim_xB_e		, "Bjorken-x"    		, "Counts [a.u.]"       	, 8 );
	TH1D_Style( dat_h1_sim_Q2_e		, "Q2 [GeV2]" 		    	, "Counts [a.u.]"       	, 8 );
	TH1D_Style( dat_h1_sim_W_e		, "W [GeV]"     		, "Counts [a.u.]"       	, 8 );
	
	cout<<"Making TCanvas.\n";
	
	TCanvas * c = new TCanvas("c","c",700,900);
	TPaveText *pt = new TPaveText(.05,.1,.95,.8);
	pt->AddText(Form("Gated Charge Collected: %f [microC]", 0));
	pt->AddText(Form("(e,e) events: %i", 0));
	pt->Draw();
	c->Modified() ; c->Update();
	c->Print("ee_counts.pdf(");
	
	cout<<"Printing to TCanvas.\n";

	TCanvas ** c1 = new TCanvas*[50];
	c1[0] = print1Hists( dat_h1_sim_p_e			, sim_h1_sim_p_e 			, 0 , c );
	c1[1] = print1Hists( dat_h1_sim_th_e			, sim_h1_sim_th_e	 		, 1 , c );
	c1[2] = print2Hists( dat_h2_sim_th_ph_e			, sim_h2_sim_th_ph_e	 		, 2 , c );
	c1[3] = print1Hists( dat_h1_sim_xB_e			, sim_h1_sim_xB_e	 		, 3 , c );
	c1[4] = print1Hists( dat_h1_sim_Q2_e			, sim_h1_sim_Q2_e	 		, 4 , c );
	c1[5] = print1Hists( dat_h1_sim_W_e			, sim_h1_sim_W_e	 		, 5 , c );
	for( int i = 0 ; i < 6 ; i++ )
		c1[i]->Print("ee_counts.pdf");

	cout<<"Writing out PDF.\n";

	c->Print("ee_counts.pdf)");


	return 1;
}

void TH1D_Style( TH1D * hist, TString xLabel, TString yLabel, int line_color ){
	hist->GetXaxis()->SetTitle(xLabel);
	hist->GetYaxis()->SetTitle(yLabel);
	hist->SetTitle("");
	hist->SetLineColor(line_color);
	hist->SetLineWidth(3);
	hist->SetStats(0);
	hist->SetMinimum(0);
}
void TH2D_Style( TH2D * hist, TString xLabel, TString yLabel, int line_color ){
	hist->GetXaxis()->SetTitle(xLabel);
	hist->GetYaxis()->SetTitle(yLabel);
	hist->SetTitle("");
	hist->SetLineColor(line_color);
	hist->SetLineWidth(3);
	hist->SetStats(0);
}

TCanvas * print1Hists( TH1D * dat , TH1D * sim , int idx , TCanvas * saveon ){
	TCanvas * printon = new TCanvas(Form("c1_%i",idx),Form("c1_%i",idx),700,900);
	printon->Divide(1,2);

	double scale_dat = (double)sim->Integral() / (double)dat->Integral();
	dat->Scale(scale_dat);
	dat->SetMinimum(0);
	dat->SetTitle("Data");
	sim->SetTitle("Simulation");

	printon->cd(1);
	dat->Draw("H,COLZ");
	printon->cd(2);
	sim->Draw("H,COLZ");

	printon->Modified(); printon->Update();

	return printon;
}
TCanvas * print2Hists( TH2D * dat , TH2D * sim , int idx , TCanvas * saveon ){
	TCanvas * printon = new TCanvas(Form("c1_%i",idx),Form("c1_%i",idx),700,900);
	printon->Divide(1,2);

	dat->SetTitle("Data");
	sim->SetTitle("Simulation");

	printon->cd(1);
	dat->Draw("H,COLZ");
	printon->cd(2);
	sim->Draw("H,COLZ");

	printon->Modified(); printon->Update();

	return printon;
}
