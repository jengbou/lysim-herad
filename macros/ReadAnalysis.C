#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <cmath>
#include <vector>
#include <map>
#include <TROOT.h>
#include <TMath.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TGraphErrors.h>
#include <TLegend.h>
#include <TAxis.h>
#include <TImageDump.h>
#include <TFile.h>
#include <TTree.h>
#include "Math/Interpolator.h"

using TMath::Exp;
using TMath::Power;
using TMath::Log;

void readRatios(string inputFile);

// Process one line of the GEANT output. 
void processLine(string line, double& absorption, double& efficiency)
{
	stringstream ss(line);
	string field;
	int ctr = 0; //column index of tab-delimited row
	while(getline(ss, field, '\t'))
	{
		if(ctr==0) {stringstream temp(field); temp >> absorption;}
		else if (ctr==1) {stringstream temp(field); temp >> efficiency;}
		else {}
		
		ctr++;
	}
}

// Read one GEANT output file. 
// Produce a TGraphErrors of efficiency vs absorption coefficient
// with vertical error bars calculated assuming 1e5 events.
TGraphErrors* readAnalysisFile(string inputFile)
{
	vector<double> absorptionV;
	vector<double> efficiencyV;
	vector<double> errorV;
	TGraphErrors* gr;
	
	ifstream inputStream(inputFile.c_str());
	string line;
	if(inputStream.is_open())
	{
		while(getline(inputStream, line))
		{
			if(line[0]=='#') continue; //Skip comment lines
			else 
			{
				//Process the line to get beam position and efficiency
				double absorption;
				double efficiency;
				double errorEfficiency;
				processLine(line, absorption, efficiency);
				absorptionV.push_back(absorption);
				efficiencyV.push_back(efficiency);
				// Binomial error from 1e5 events.
				errorEfficiency = sqrt(efficiency * (1 - efficiency) / 1e5);
				errorV.push_back(errorEfficiency);
			}
		}
		Int_t n = absorptionV.size();
		gr = new TGraphErrors(n, &absorptionV[0], &efficiencyV[0], NULL, &errorV[0]);
		//cout << gr << endl;
		//gr->Draw("ALP");
	}
	else {
		cout<<"Error: Input file not found." <<endl;
		gr = NULL;
	}
	//cout << gr << endl;
	return gr;
}

string constructFileName(const char* runName, int layerNo, int ieta)
{
	TString fileName("../RUNNAME/Analysis-RUNNAME_layerLAYERNUMBER_ietaIETA.txt");
	fileName.ReplaceAll("RUNNAME", runName);
	fileName.ReplaceAll("LAYERNUMBER", TString::Itoa(layerNo, 10));
	fileName.ReplaceAll("IETA", TString::Itoa(ieta, 10));
	string fileNameString(fileName.Data());
	return  fileNameString;
}

void processDoseLine(string line, int& ieta, double& dose, double& fluence)
{
	stringstream ss(line);
	string field;
	int ctr = 0; //column index of tab-delimited row
	while(getline(ss, field, '\t'))
	{
		if(ctr==0) {stringstream temp(field); temp >> ieta;}
		else if (ctr==4) {stringstream temp(field); temp >> dose;}
		else if (ctr==5) {stringstream temp(field); temp >> fluence;}
		else {}
		
		ctr++;
	}
}

void readDoseFile(string inputFile, std::map<int, double>& doseMap, std::map<int, double>& fluenceMap)
{
	ifstream inputStream(inputFile.c_str());
	string line;
	if(inputStream.is_open())
	{
		while(getline(inputStream, line))
		{
			if(line[0]=='#') continue; //Skip comment lines
			else 
			{
				//Process the line to get ieta and dose
				int ieta;
				double dose, fluence;
				processDoseLine(line, ieta, dose, fluence);
				doseMap.insert(std::pair<int, double>(ieta, dose));
				fluenceMap.insert(std::pair<int, double>(ieta, fluence / 1e12));
			}
		}
	}
	else {cout <<"Error: Input file not found." <<endl;}
}

void dosePerTile()
{
	////////////////////////////////////////////////////////
	// Calculate dose per 100fb^-1 for each tile
	////////////////////////////////////////////////////////
	std::map<int, double> doseMapLayer1, fluenceMapLayer1;
	readDoseFile("doseLayer1.txt", doseMapLayer1, fluenceMapLayer1);
	cout << "ieta \t dose" << endl;
	for (int i=17; i<=29; i++) {
		cout << i << "\t" << doseMapLayer1[i] << endl;
	}
}

// Function to model efficiency vs attenuation length.
Double_t myfunction(const Double_t *x, Double_t *par)
{
	Double_t xx =x[0];
	Double_t f = par[1] * Exp(-par[0]/(xx));
	return f;
}

// Function to model efficiency vs 1 / attenuation length.
Double_t myfunction2(const Double_t *x, Double_t *par)
{
	Double_t xx =x[0];
	Double_t f = par[1] * Exp(-par[0]*(xx));
	return f;
}

// Function to model efficiency vs 1 / attenuation length.
Double_t EfficiencyVsMu(const Double_t *x, Double_t *par)
{
	Double_t xx =x[0];
	Double_t f = Exp(par[3]*pow(xx, 3) + par[2]*pow(xx, 2) + par[1]*(xx) + par[0]);
	return f;
}

// Inverse of EfficiencyVsMu
/*
Double_t MuVsEfficiency(const Double_t *x, Double_t *par)
{
	Double_t xx =x[0];
	Double_t f = (-par[1] - sqrt(pow(par[1], 2) - 4*par[2]*(par[0] - Log(xx)))) / (2*par[2]);
	return f;
}
*/

// Test EfficiencyVsMu and inverse
/*
void test()
{
	double par[3] = {-1.5, -46, 254};
	double x[1], y[1], z[1];
	x[0] = 0.0101;
	y[0] = EfficiencyVsMu(x, par);
	z[0] = MuVsEfficiency(y, par);
	cout << y[0] << endl;
	cout << z[0] << endl;
}
*/

// Function to model efficiency vs 1 / attenuation length.
Double_t LogEfficiencyVsMu(const Double_t *x, Double_t *par)
{
	Double_t xx = x[0];
	Double_t f = par[2] * pow(xx, 2) + par[1] * xx + par[0];
	return f;
}

// Inverse of myfunction 
Double_t invmyfunction(const Double_t *f, Double_t *par)
{
	Double_t ff = f[0];
	Double_t x = par[0] / Log( par[1]/ff );
	return x;
}

// Inverse of myfunction2 
Double_t invmyfunction2(const Double_t *f, Double_t *par)
{
	Double_t ff = f[0];
	Double_t x = - Log(ff / par[1]) / par[0];
	return x;
}

// Function to map dose in kGy to mu(Tile).
// mu(Fiber) is implicitly set in simulation. 
Double_t DoseToMu(const Double_t *x, Double_t *par)
{
	Double_t xx = x[0];
	Double_t f = par[1] * xx + par[0];
	return f;
}

// Function to map (dose, dose rate) in kGy to mu(Tile).
// mu(Fiber) is implicitly set in simulation. 
Double_t CalculateMu(const Double_t *x, Double_t *par)
{
	Double_t dose = x[0];
	Double_t doseRate = x[1];
	Double_t doseRateBase = par[2];
	Double_t A = par[3];
	Double_t f = ( 1 - A * log10(doseRate/doseRateBase) ) * par[1] * dose + par[0];
	return f;
}

void myfunc()
{
//	TF1 *f1 = new TF1("myfunc", myfunction, 0.1, 100.0, 2);
//	Double_t par[] = {1.0, 1.0};
//	f1->SetParameters(par);
//	f1->SetParNames("average track length", "constant");
//	TGraphErrors* gr = longrun6test();
//	gr->GetXaxis()->SetTitle("Absorption length (cm)");
//	gr->GetYaxis()->SetTitle("Efficiency");
//	gr->GetXaxis()->SetLimits(0.1,300);
//	gr->GetYaxis()->SetRangeUser(1e-7,1e-0);
//	gr->Draw("AP");
//	gr->Fit(f1, "", "", 0.5, 100);
//	f1->GetParameters(par);
//	cout << par[0] << "\t" << par[1] << endl;
}

void processLaserDataLine(string line, double& lumi, double* ratios)
{
	// ratio[0] through ratio[12] corresponds to ieta 17 to 29.
	// Make sure "ratios" has 13 elements.
	stringstream ss(line);
	string field;
	int ctr = 0; //column index of tab-delimited row
	while(getline(ss, field, '\t'))
	{
		if(ctr==0) {
			stringstream temp(field);
			temp >> lumi;
			//cout << lumi;
		}
		else if (true) {//ctr >= 1 && ctr <= 13) {
			stringstream temp(field); 
			double ratio;
			temp >> ratio;
			//cout << ratio << "\t";
			ratios[ctr - 1] = ratio;
		}
		else {
			cout << endl;
			break;
		}
		
		ctr++;
	}
}

void readLaserDataFile(string inputFile, vector<double>& lumis, std::map< int, vector<double> > & ratiosMap)
{
	ifstream inputStream(inputFile.c_str());
	string line;
	if(inputStream.is_open())
	{
		while(getline(inputStream, line))
		{
			if(line[0]=='#') continue; //Skip comment lines
			else 
			{
				//Process the line to get ieta and dose
				double lumi;
				double ratios[13];
				processLaserDataLine(line, lumi, ratios);
				lumis.push_back(lumi);
				for (int i=0; i < 13; i++) {
					ratiosMap[i+17].push_back(ratios[i]);
					//cout << "ratios[i] " <<ratios[i] << endl;
					//cout << "ratiosMap[i+17] " << (ratiosMap[i+17])[1] << endl;
				}
			}
		}
	}
	else {cout <<"Error: Input file not found." <<endl;}
}

void readLaserDataFile2(string inputFile, vector<double>& lumis, std::map< int, vector<double> > & ratiosMap)
{
	ifstream inputStream(inputFile.c_str());
	string line;
	if(inputStream.is_open())
	{
		while(getline(inputStream, line))
		{
			if(line[0]=='#') continue; //Skip comment lines
			else 
			{
				//Process the line to get ieta and dose
				double lumi;
				double ratios[13];
				processLaserDataLine(line, lumi, ratios);
				lumis.push_back(lumi);
				for (int i=0; i < 13; i++) {
					ratiosMap[i+17].push_back(ratios[i]);
					//cout << "ratios[i] " <<ratios[i] << endl;
					//cout << "ratiosMap[i+17] " << (ratiosMap[i+17])[1] << endl;
				}
			}
		}
	}
	else {cout <<"Error: Input file not found." <<endl;}
}

void ratiosMap(int ieta = 29) //bool doPlot = false)
{
	vector<double> lumis;
	std::map<int, vector<double> > ratiosMap;

	readLaserDataFile("../HELaserData.txt", lumis, ratiosMap);

	TCanvas* c1 = new TCanvas("c1","canvas1",800,600);
	TGraph* graph;
	graph = new TGraph(lumis.size(), &lumis[0], &((ratiosMap[ieta])[0]));
	c1->cd();
	graph->Draw("AP");
}

// Function to model attenuation length vs dose 
Double_t LambdaVsDose(const Double_t *x, Double_t *par)
{
	Double_t xx =x[0];
	Double_t f = 1 / (par[0] * xx + par[1]);
	return f;
}

// Function to model 1/attenuation length vs dose 
Double_t LambdaVsDose2(const Double_t *x, Double_t *par)
{
	Double_t xx =x[0];
	Double_t f = par[0] * (xx - 0) + par[1];
	return f;
}

// Deprecated code
/*
void temp(int minLumiIndex, int maxLumiIndex, double refAttLength, double lumioffset, bool doPlot = false, bool doSavePlots = false, bool doDoseOrFluence = true) 
{
	// Vector of colors.
	vector<Color_t> colors(18, kBlack);
	{
		colors[0] = kBlue;
		colors[1] = kOrange;
		colors[2] = kGreen;
		colors[3] = kRed;
		colors[4] = kCyan;
		colors[5] = kMagenta;
		// Define colors[6] through colors[11]
		for (int i = 0; i < 6; i++ ) {
			colors[i+6] = colors[i] + 1;
		}
		// Define colors[12] through colors[17]
		for (int i = 0; i < 6; i++ ) {
			colors[i+12] = colors[i] + 2;
		}
	}
	
	TCanvas* c1;
	c1 = new TCanvas("c1", "Efficiency vs Mu", 800, 600);
	//gPad->SetLogx();
	//gPad->SetLogy();
	
	// graphsAll contains Efficiency vs Mu data for all tiles from the Geant4 simulations.
	// namesAll contains names for all graphs in graphsAll.
	// graphsTypeX contains Efficiency vs Lambda data for tiles of type X from the Geant4 simulations.
	vector<TGraphErrors*> graphsLayer1, graphsLayer7, graphsAll, graphsAllLog;
	vector<string> namesAll;

	// Read from Analysis files and add to graphsTypeX, and graphsAll.
	// (Only reading type 1 tiles.)
	{
		graphsLayer1.push_back(readAnalysisFile("../run5b/Analysis-run5b_layer1_ieta22.txt"));
		graphsLayer1.push_back(readAnalysisFile("../run5b/Analysis-run5b_layer1_ieta23.txt"));
		graphsLayer1.push_back(readAnalysisFile("../run5b/Analysis-run5b_layer1_ieta24.txt"));
		graphsLayer1.push_back(readAnalysisFile("../run5b/Analysis-run5b_layer1_ieta25.txt"));
		graphsLayer1.push_back(readAnalysisFile("../run5b/Analysis-run5b_layer1_ieta26.txt"));
		graphsLayer1.push_back(readAnalysisFile("../run5b/Analysis-run5b_layer1_ieta27.txt"));
		graphsLayer1.push_back(readAnalysisFile("../run5b/Analysis-run5b_layer1_ieta28.txt"));
		graphsLayer1.push_back(readAnalysisFile("../run5b/Analysis-run5b_layer1_ieta29.txt"));
		
	//	graphsType2.push_back(readAnalysisFile("../run5/Analysis-run5_layer1_ieta17.txt"));
	//	graphsType2.push_back(readAnalysisFile("../run5/Analysis-run5_layer1_ieta18.txt"));
	//	graphsType2.push_back(readAnalysisFile("../run5/Analysis-run5_layer1_ieta19.txt"));
	//	graphsType2.push_back(readAnalysisFile("../run5/Analysis-run5_layer1_ieta20.txt"));
	//	
	//	graphsType3.push_back(readAnalysisFile("../run5/Analysis-run5_layer1_ieta17.txt"));
	//	graphsType3.push_back(readAnalysisFile("../run5/Analysis-run5_layer1_ieta18.txt"));
	//	graphsType3.push_back(readAnalysisFile("../run5/Analysis-run5_layer1_ieta19.txt"));
	//	graphsType3.push_back(readAnalysisFile("../run5/Analysis-run5_layer1_ieta20.txt"));
		
		// Add all graphs from graphsLayer1 to graphsAll. Add all legend names to namesAll.
		for (unsigned int i=0; i < graphsLayer1.size(); i++) {
			graphsAll.push_back(graphsLayer1[i]);
			stringstream name;
			name << "ieta " << (i+22) << "Layer1";
			namesAll.push_back(name.str());
		}
	}

	TLegend* leg = new TLegend(0.7, 0.6, 0.9, 0.85);
	// Run through all graphs in graphsAll and modify error bars if necessary.
	// Push log version of graphsAll into graphsAllLog.
	// Set graph properties of all graphs in graphsAll and draw to c1 pad 1.
	// Add all legend names from namesAll to leg.
	for (unsigned int i=0; i < graphsAll.size(); i++) {
		//cout << graphsAll.size() << ": " << graphsAll[i] << endl;
		
		int n = graphsAll[i]->GetN();
		double* xV = graphsAll[i]->GetX();
		// Offset x-values for clarity.
		// (Disabled)
		for (int j=0; j < n; j++) {
			//xV[j] += i*0.0002;
		}
		double* yV = graphsAll[i]->GetY();
		double* exV = graphsAll[i]->GetEX();
		double* eyV = graphsAll[i]->GetEY();
	
		graphsAll[i] = new TGraphErrors(n, xV, yV, exV, eyV);
		graphsAll[i]->SetTitle("");//("Light collection efficiency vs Absorption Length");
		graphsAll[i]->GetXaxis()->SetTitle("Absorption coefficient [cm^{-1}]");
		graphsAll[i]->GetYaxis()->SetTitle("Efficiency");
		//graphsAll[i]->GetXaxis()->SetLimits(0.0, 0.05);
		graphsAll[i]->GetYaxis()->SetRangeUser(1e-2, 2e-1);
		graphsAll[i]->SetMarkerColor(colors[i]);
		graphsAll[i]->SetMarkerSize(0.5);
		graphsAll[i]->SetMarkerStyle(21);

		graphsAllLog.push_back( (TGraphErrors*) graphsAll[i]->Clone() );

		c1->cd(1);
		c1->SetLogy();
		if(i==0) {
			graphsAll[i]->Draw("AP");
		}
		else {
			graphsAll[i]->Draw("P same");
		}
	
		leg->AddEntry(graphsAll[i], namesAll[i].c_str(), "p");
	}

	vector<TF1*> functionV (graphsAll.size(), new TF1("EfficiencyVsMu", EfficiencyVsMu, 0.0, 0.0, 3));
	// Fit all graphs in graphsAllLog to myfunction.
	for (unsigned int i=0; i < graphsAll.size(); i++) {
		Double_t par[] = {1.0, 1.0, 1.0, 1.0};
		functionV[i]->SetParameters(par);
		functionV[i]->SetLineColor(colors[i]);
		functionV[i]->SetParNames("p0", "p1", "p2");
		graphsAll[i]->Fit(functionV[i], "", "", 0.01, 0.045);
	}

	// Make log version of graphsAllLog
	for (unsigned int i=0; i < graphsAllLog.size(); i++) {
		//cout << graphsAllLogLog.size() << ": " << graphsAllLog[i] << endl;
		
		int n = graphsAllLog[i]->GetN();
		double* xV = graphsAllLog[i]->GetX();
		double* yV = graphsAllLog[i]->GetY();
		for (int j=0; j < n; j++) {
			yV[j] = log10( yV[j] );
		}
		double* exV = graphsAllLog[i]->GetEX();
		double* eyV = graphsAllLog[i]->GetEY();

		graphsAllLog[i] = new TGraphErrors(n, xV, yV, exV, eyV);
		graphsAllLog[i]->SetTitle("");//("Light collection efficiency vs Absorption Length");
		graphsAllLog[i]->GetXaxis()->SetTitle("Absorption coefficient [cm^{-1}]");
		graphsAllLog[i]->GetYaxis()->SetTitle("Log(Efficiency)");
		//graphsAllLog[i]->GetXaxis()->SetLimits(-0.005, 0.045);
		//graphsAllLog[i]->GetYaxis()->SetRangeUser(1e-6, 1e0);
		graphsAllLog[i]->SetMarkerColor(colors[i]);
		graphsAllLog[i]->SetMarkerSize(0.5);
		graphsAllLog[i]->SetMarkerStyle(21);

		//c1->cd(1);
		//if(i==0) {
		//	graphsAllLog[i]->Draw("AP");
		//}
		//else {
		//	graphsAllLog[i]->Draw("P same");
		//}
	}

	// Draw leg to c1 pad 1.
	{
		c1->cd(1);
		leg->SetFillColor(0);
		leg->SetBorderSize(0);
		leg->SetTextSize(0.03);
		leg->SetTextFont(42);
		leg->Draw();
	}

	vector<TF1*> functionLogV (graphsAllLog.size(), new TF1("LogEfficiencyVsMu", LogEfficiencyVsMu, 0.0, 0.0, 3));
	// Fit all graphs in graphsAllLog to myfunction.
	for (unsigned int i=0; i < graphsAllLog.size(); i++) {
		Double_t par[] = {1.0, 1.0, 1.0, 1.0};
		functionLogV[i]->SetParameters(par);
		functionLogV[i]->SetLineColor(colors[i]);
		functionLogV[i]->SetParNames("p0", "p1", "p2");
		graphsAllLog[i]->Fit(functionLogV[i], "", "", 0.01, 0.045);
	}

//	vector<TF1*> function2V (graphsAll.size(), new TF1("myfunction2", myfunction2, 0.01, 100.0, 2));
//	// Create graphs of Efficiency vs Mu from graphsAll and add to graphsAll2.
//	// Run through all graphs in graphsAll2 and modify error bars if necessary.
//	// Set graph properties of all graphs in graphsAll2 and draw to c1 pad 2.
//	for (unsigned int i=0; i < graphsAll.size(); i++) {
//		graphsAll2.push_back( (TGraphErrors*) (graphsAll[i])->Clone() );
//		
//		int n = graphsAll2[i]->GetN();
//
//		double* xV = graphsAll2[i]->GetX();
//		vector<double> xVector;
//		for (int j=0; j < n; j++) {
//			xVector.push_back(1 / xV[j]);
//		}
//		double* yV = graphsAll2[i]->GetY();
//		double* exV = graphsAll2[i]->GetEX();
//		double* eyV = graphsAll2[i]->GetEY();
//
//		graphsAll2[i] = new TGraphErrors(n, &xVector[0], yV, exV, eyV);
//		graphsAll2[i]->SetTitle("");//("Light collection efficiency vs Absorption Length");
//		graphsAll2[i]->GetXaxis()->SetTitle("1 / Absorption length (cm)");
//		graphsAll2[i]->GetYaxis()->SetTitle("Efficiency");
//		graphsAll2[i]->GetXaxis()->SetLimits(0.005, 0.25);
//		graphsAll2[i]->GetYaxis()->SetRangeUser(-0.02, 0.14);
//		graphsAll2[i]->SetMarkerColor(colors[i]);
//		graphsAll2[i]->SetMarkerSize(0.5);
//		graphsAll2[i]->SetMarkerStyle(21);
//
//		c1->cd(2);
//		if(i==0) {
//			graphsAll2[i]->Draw("AP");
//		}
//		else {
//			graphsAll2[i]->Draw("P same");
//		}
//	}
//
//	// Save c1 as image.
//	{
//		c1->SaveAs("Efficiency vs Lambda and Mu.png");
//	}
//
//	// Fit all graphs in graphsAll2 to myfunction2.
//	for (unsigned int i=0; i < graphsAll2.size(); i++) {
//		Double_t par[] = {1.0, 0.5};
//		function2V[i]->SetParameters(par);
//		function2V[i]->SetLineColor(colors[i]);
//		function2V[i]->SetParNames("average track length", "constant", "offset");
//		graphsAll2[i]->Fit(function2V[i], "", "", 0.005, 10);
//	}
//
//	vector<double> trackLengths, ietas;
//	// Create vectors of average track lengths and ieta from graphsAll2.
//	for (unsigned int i=0; i < function2V.size(); i++) {
//		trackLengths.push_back( (function2V[i]->GetParameters())[0] );
//		ietas.push_back( i+21 );
//	}
//
//	// Create graph of average track lengths vs ieta.
//	// Draw to new canvas.
//	if (doPlot) {
//		TCanvas* trackLengthCanvas = new TCanvas("track lengths", "track lengths", 800, 600);
//
//		TGraph* graphTrackLength = new TGraph( trackLengths.size(), &ietas[0], &trackLengths[0] );
//		graphTrackLength->SetTitle("");
//		graphTrackLength->GetXaxis()->SetTitle("ieta");
//		graphTrackLength->GetYaxis()->SetTitle("Average track length (cm)");
//		graphTrackLength->GetXaxis()->SetLabelSize(0.035);
//		graphTrackLength->GetYaxis()->SetLabelSize(0.035);
//	//	graphTrackLength->GetXaxis()->SetLimits(0.3,200);
//	//	graphTrackLength->GetYaxis()->SetRangeUser(1e-6,3e-1);
//		graphTrackLength->SetMarkerColor(colors[0]);
//		graphTrackLength->SetMarkerSize(0.5);
//		graphTrackLength->SetMarkerStyle(21);
//
//		trackLengthCanvas->cd();
//		graphTrackLength->Draw("AP");
//	}
//
	// Read laser data file.
	vector<double> lumis;
	std::map<int, vector<double> > ratiosMap;
	readLaserDataFile("../HELaserData.txt", lumis, ratiosMap);

//	// Read dose and fluence map data files.
//	std::map<int, double> doseMapLayer1;
//	std::map<int, double> fluenceMapLayer1;
//	readDoseFile("doseLayer1.txt", doseMapLayer1, fluenceMapLayer1);
//
//	TLegend* legEfficiencyDose = new TLegend(0.7, 0.2, 0.9, 0.6);
//	TLegend* legEfficiencyFluence = new TLegend(0.7, 0.2, 0.9, 0.6);
//	//vector<TGraph*> graphsEfficiencyDoseLinear, graphsEfficiencyDoseLog;
//	vector<double>  kParameterDoseVector, kParameterFluenceVector, ietaVector;
//	
//	// Run through all ietas. 
//	for (int i=29; i>=21; i--)
//	{
//		stringstream name;
//		name << "ieta " << (i);
//		TCanvas *canvasDose, *canvasFluence;
//		// Create canvases with 2 pads each. One for Dose and one for Fluence.
//		if (doPlot) {
//			canvasDose = new TCanvas(name.str().c_str(), name.str().c_str(), 1200, 600);
//			canvasDose->SetLogx(0);
//			canvasDose->SetLogy(0);
//			canvasDose->Divide(2,0);
//			canvasFluence = new TCanvas(name.str().c_str(), name.str().c_str(), 1200, 600);
//			canvasFluence->SetLogx(0);
//			canvasFluence->SetLogy(0);
//			canvasFluence->Divide(2,0);
//		}
//
//		vector<double> xVector, xVector2, yVector, yVector2;
//		// Calculate measured Lambda, Mu, Dose, Fluence and add to vectors	
//		for (int j = minLumiIndex; j < maxLumiIndex + 1; j++)
//		{
//			double* fitParameters = functionV[i - 21]->GetParameters();
//			double* fitParameters2 = function2V[i - 21]->GetParameters();
//			double refEfficiency = myfunction( &refAttLength, fitParameters );
//			double measuredEfficiency = ratiosMap[i][j] * refEfficiency;
//			double measuredLambda = invmyfunction( &measuredEfficiency, fitParameters );
//			double measuredMu = invmyfunction2( &measuredEfficiency, fitParameters2 );
//
//			double refDose = doseMapLayer1[i];
//			double measuredDose = refDose * (lumis[j] - lumioffset) / 100;
//
//			double refFluence = fluenceMapLayer1[i];
//			double measuredFluence = refFluence * (lumis[j] - lumioffset) / 100;
//
//			xVector.push_back(measuredDose);
//			xVector2.push_back(measuredFluence);
//			yVector.push_back(measuredLambda);
//			yVector2.push_back(measuredMu);
//			
//			//cout << "i, j = " << i << ", " << j << endl;
//		}
//
//		TGraph *graphDose, *graphFluence, *graphDose2, *graphFluence2;
//		// Create graphs of Lambda vs Dose and Lambda vs Fluence.
//		{
//	   		graphDose = new TGraph( (maxLumiIndex - minLumiIndex + 1), &xVector[0], &yVector[0] );
//			graphDose->GetXaxis()->SetTitle("Dose (kGy)");
//			graphDose->SetTitle("");
//			graphDose->GetYaxis()->SetTitle("Extracted absorption length (cm)");
//			graphDose->GetXaxis()->SetLabelSize(0.035);
//			graphDose->GetYaxis()->SetLabelSize(0.035);
//	//		graphDose->GetXaxis()->SetLimits(0.3,200);
//	//		graphDose->GetYaxis()->SetRangeUser(1e-6,3e-1);
//			graphDose->SetMarkerColor(colors[0]);
//			graphDose->SetMarkerSize(0.5);
//			graphDose->SetMarkerStyle(21);
//
//	   		graphFluence = new TGraph( (maxLumiIndex - minLumiIndex + 1), &xVector2[0], &yVector[0] );
//			graphFluence->GetXaxis()->SetTitle("Neutral hadron fluence (10^12 cm^-1)");
//			graphFluence->SetTitle("");
//			graphFluence->GetYaxis()->SetTitle("Extracted absorption length (cm)");
//			graphFluence->GetXaxis()->SetLabelSize(0.035);
//			graphFluence->GetYaxis()->SetLabelSize(0.035);
//	//		graphFluence->GetXaxis()->SetLimits(0.3,200);
//	//		graphFluence->GetYaxis()->SetRangeUser(1e-6,3e-1);
//			graphFluence->SetMarkerColor(colors[0]);
//			graphFluence->SetMarkerSize(0.5);
//			graphFluence->SetMarkerStyle(21);
//		}
//
//		// Fit graphDose and graphFluence to LambdaVsDose.
//		{
//			TF1* functionDose = new TF1("LambdaVsDose", LambdaVsDose, 0.0, 1e+5, 2);
//			double parameters[] = {0.1, 0.1, 1.0};
//			functionDose->SetParameters(parameters);
//			functionDose->SetLineColor(colors[1]);
//			functionDose->SetParNames("p0", "p1", "p2");
//			graphDose->Fit(functionDose, "", "", 0.0, 1e+5);
//
//			TF1* functionFluence = (TF1*) functionDose->Clone();
//			graphFluence->Fit(functionFluence, "", "", 0.0, 1e+5);
//		}
//
//		legEfficiencyDose->AddEntry(graphDose, name.str().c_str(), "p");
//		legEfficiencyFluence->AddEntry(graphFluence, name.str().c_str(), "p");
//
//		// Draw graphDose and graphFluence to canvasDose and canvasFluence (pad 1).
//		if (canvasDose && canvasFluence && doPlot) {
//			canvasDose->cd(1);
//			graphDose->Draw("AP");
//			canvasFluence->cd(1);
//			graphFluence->Draw("AP");
//		}
//
//		// Create graphs of Mu vs Dose and Mu vs Fluence.
//		{
//	   		graphDose2 = new TGraph( (maxLumiIndex - minLumiIndex + 1), &xVector[0], &yVector2[0] );
//			graphDose2->GetXaxis()->SetTitle("Dose (kGy)");
//			graphDose2->SetTitle("");
//			graphDose2->GetYaxis()->SetTitle("1 / Extracted asorption length (cm^-1)");
//			graphDose2->GetXaxis()->SetLabelSize(0.035);
//			graphDose2->GetYaxis()->SetLabelSize(0.035);
//		//	graphDose2->GetXaxis()->SetLimits(0.3,200);
//		//	graphDose2->GetYaxis()->SetRangeUser(1e-6,3e-1);
//			graphDose2->SetMarkerColor(colors[0]);
//			graphDose2->SetMarkerSize(0.5);
//			graphDose2->SetMarkerStyle(21);
//
//	   		graphFluence2 = new TGraph( (maxLumiIndex - minLumiIndex + 1), &xVector2[0], &yVector2[0] );
//			graphFluence2->GetXaxis()->SetTitle("Neutral hadron fluence (10^12 cm^-1)");
//			graphFluence2->SetTitle("");
//			graphFluence2->GetYaxis()->SetTitle("1 / Extracted asorption length (cm^-1)");
//			graphFluence2->GetXaxis()->SetLabelSize(0.035);
//			graphFluence2->GetYaxis()->SetLabelSize(0.035);
//		//	graphFluence2->GetXaxis()->SetLimits(0.3,200);
//		//	graphFluence2->GetYaxis()->SetRangeUser(1e-6,3e-1);
//			graphFluence2->SetMarkerColor(colors[0]);
//			graphFluence2->SetMarkerSize(0.5);
//			graphFluence2->SetMarkerStyle(21);
//		}
//
//		// Fit graphDose2 and graphFluence2 to LambdaVsDose2.
//		// Add k-parameter of Mu vs Dose and Mu vs Fluence to vectors.
//		{
//			TF1* functionDose2 = new TF1("LambdaVsDose", LambdaVsDose2, 0.0, 1e+5, 2);
//			double parameters[] = {0.1, 0.1, 1.0};
//			functionDose2->SetParameters(parameters);
//			functionDose2->SetLineColor(colors[1]);
//			functionDose2->SetParNames("p0", "p1", "p2");
//			graphDose2->Fit(functionDose2, "", "", 0.0, 1e+5);
//
//			TF1* functionFluence2 = (TF1*) functionDose2->Clone();
//			graphFluence2->Fit(functionFluence2, "", "", 0.0, 1e+5);
//
//			kParameterDoseVector.push_back( (functionDose2->GetParameters())[0] );
//			kParameterFluenceVector.push_back( (functionFluence2->GetParameters())[0] );
//			ietaVector.push_back(i);
//		}
//
//		// Draw graphDose2 and graphFluence2 to canvasDose and canvasFluence (pad 2).
//		if (canvasDose && canvasFluence && doPlot) {
//			canvasDose->cd(2);
//			graphDose2->Draw("AP");
//			canvasFluence->cd(2);
//			graphFluence2->Draw("AP");
//		}
//
//	//legEfficiencyDose->SetFillColor(0);
//	//legEfficiencyDose->SetBorderSize(0);
//	//legEfficiencyDose->SetTextSize(0.03);
//	//legEfficiencyDose->SetTextFont(42);
//	//legEfficiencyDose->Draw();
//
//		// Save Lambda (Mu) vs Dose and Lambda (Mu) vs Fluence plots.
//		if (doSavePlots && doPlot) {
//			stringstream filenameDose, filenameFluence;
//			filenameDose << "plots/LambdaVSDose_ieta" << (i) << ".png";
//			filenameFluence << "plots/LambdaVSFluence_ieta" << (i) << ".png";
//			TImageDump *imgdumpDose = new TImageDump(filenameDose.str().c_str());
//			canvasDose->Paint();
//			imgdumpDose->Close();
//			TImageDump *imgdumpFluence = new TImageDump(filenameFluence.str().c_str());
//			canvasFluence->Paint();
//			imgdumpFluence->Close();
//		}
//	}
//
//	// Plot k-parameter vs ieta. k is slope of Mu vs Dose, or Mu vs Fluence.
//	{
//		string titleDose, ytitleDose;
//		titleDose = "k-parameter [Dose]";
//		ytitleDose = "k-parameter [cm^{-1} / kGy]";
//		string titleFluence, ytitleFluence;
//		titleFluence = "k-parameter [Fluence]";
//		ytitleFluence = "k-parameter [cm^{-1} / cm^{-2}]";
//
//		TCanvas* canvasDose = new TCanvas(titleDose.c_str(), titleDose.c_str(), 800, 600);
//
//		TGraph* graphDose = new TGraph( kParameterDoseVector.size(), &ietaVector[0], &kParameterDoseVector[0] );
//		graphDose->SetTitle("");
//		graphDose->GetXaxis()->SetTitle("ieta");
//		graphDose->GetYaxis()->SetTitle(ytitleDose.c_str());
//		graphDose->GetXaxis()->SetLabelSize(0.035);
//		graphDose->GetYaxis()->SetLabelSize(0.035);
//	//	graphDose->GetXaxis()->SetLimits(0.3,200);
//	//	graphDose->GetYaxis()->SetRangeUser(1e-6,3e-1);
//		graphDose->SetMarkerColor(colors[0]);
//		graphDose->SetMarkerSize(1.0);
//		graphDose->SetMarkerStyle(21);
//
//		canvasDose->cd();
//		graphDose->Draw("AP");
//		canvasDose->SaveAs("k-parameter (Dose).png");
//
//		TCanvas* canvasFluence = new TCanvas(titleFluence.c_str(), titleFluence.c_str(), 800, 600);
//
//		TGraph* graphFluence = new TGraph( kParameterFluenceVector.size(), &ietaVector[0], &kParameterFluenceVector[0] );
//		graphFluence->SetTitle("");
//		graphFluence->GetXaxis()->SetTitle("ieta");
//		graphFluence->GetYaxis()->SetTitle(ytitleFluence.c_str());
//		graphFluence->GetXaxis()->SetLabelSize(0.035);
//		graphFluence->GetYaxis()->SetLabelSize(0.035);
//	//	graphFluence->GetXaxis()->SetLimits(0.3,200);
//	//	graphFluence->GetYaxis()->SetRangeUser(1e-6,3e-1);
//		graphFluence->SetMarkerColor(colors[0]);
//		graphFluence->SetMarkerSize(1.0);
//		graphFluence->SetMarkerStyle(21);
//
//		canvasFluence->cd();
//		graphFluence->Draw("AP");
//		canvasFluence->SaveAs("k-parameter (Fluence).png");
//	}
}
*/

void MakeTrees(const char* runName, double doseFactor, const char* treeLabel) 
{
	// Vector of colors.
	vector<Color_t> colors(18, kBlack);
	{
		colors[0] = kBlue;
		colors[1] = kOrange;
		colors[2] = kGreen;
		colors[3] = kRed;
		colors[4] = kCyan;
		colors[5] = kMagenta;
		// Define colors[6] through colors[11]
		for (int i = 0; i < 6; i++ ) {
			colors[i+6] = colors[i] + 1;
		}
		// Define colors[12] through colors[17]
		for (int i = 0; i < 6; i++ ) {
			colors[i+12] = colors[i] + 2;
		}
	}
	
	TCanvas* c1;
	c1 = new TCanvas("c1", "Efficiency vs Mu", 800, 600);
	//gPad->SetLogx();
	//gPad->SetLogy();
	
	// graphsAll contains Efficiency vs Mu data for all tiles from the Geant4 simulations.
	// namesAll contains names for all graphs in graphsAll.
	// graphsTypeX contains Efficiency vs Lambda data for tiles of type X from the Geant4 simulations.
	vector<TGraphErrors*> graphsLayer1, graphsLayer7, graphsAll, graphsAllLog;
	vector<string> namesAll;

	// Read from Analysis files and add to graphsTypeX, and graphsAll.
	// (Only reading type 1 tiles.)
	{
		graphsLayer1.push_back(readAnalysisFile(constructFileName(runName, 1, 22)));
		graphsLayer1.push_back(readAnalysisFile(constructFileName(runName, 1, 23)));
		graphsLayer1.push_back(readAnalysisFile(constructFileName(runName, 1, 24)));
		graphsLayer1.push_back(readAnalysisFile(constructFileName(runName, 1, 25)));
		graphsLayer1.push_back(readAnalysisFile(constructFileName(runName, 1, 26)));
		graphsLayer1.push_back(readAnalysisFile(constructFileName(runName, 1, 27)));
		graphsLayer1.push_back(readAnalysisFile(constructFileName(runName, 1, 28)));
		graphsLayer1.push_back(readAnalysisFile(constructFileName(runName, 1, 29)));
		
	//	graphsType2.push_back(readAnalysisFile("../run5/Analysis-run5_layer1_ieta17.txt"));
	//	graphsType2.push_back(readAnalysisFile("../run5/Analysis-run5_layer1_ieta18.txt"));
	//	graphsType2.push_back(readAnalysisFile("../run5/Analysis-run5_layer1_ieta19.txt"));
	//	graphsType2.push_back(readAnalysisFile("../run5/Analysis-run5_layer1_ieta20.txt"));
	//	
	//	graphsType3.push_back(readAnalysisFile("../run5/Analysis-run5_layer1_ieta17.txt"));
	//	graphsType3.push_back(readAnalysisFile("../run5/Analysis-run5_layer1_ieta18.txt"));
	//	graphsType3.push_back(readAnalysisFile("../run5/Analysis-run5_layer1_ieta19.txt"));
	//	graphsType3.push_back(readAnalysisFile("../run5/Analysis-run5_layer1_ieta20.txt"));
		
		// Add all graphs from graphsLayer1 to graphsAll. Add all legend names to namesAll.
		for (unsigned int i=0; i < graphsLayer1.size(); i++) {
			graphsAll.push_back(graphsLayer1[i]);
			stringstream name;
			name << "ieta " << (i+22) << "Layer1";
			namesAll.push_back(name.str());
		}
	}

	TLegend* leg = new TLegend(0.7, 0.6, 0.9, 0.85);
	// Run through all graphs in graphsAll and modify error bars if necessary.
	// Push log version of graphsAll into graphsAllLog.
	// Set graph properties of all graphs in graphsAll and draw to c1 pad 1.
	// Add all legend names from namesAll to leg.
	for (unsigned int i=0; i < graphsAll.size(); i++) {
		//cout << graphsAll.size() << ": " << graphsAll[i] << endl;
		
		int n = graphsAll[i]->GetN();
		double* xV = graphsAll[i]->GetX();
		// Offset x-values for clarity.
		// (Disabled)
		for (int j=0; j < n; j++) {
			//xV[j] += i*0.0002;
		}
		double* yV = graphsAll[i]->GetY();
		double* exV = graphsAll[i]->GetEX();
		double* eyV = graphsAll[i]->GetEY();
	
		graphsAll[i] = new TGraphErrors(n, xV, yV, exV, eyV);
		graphsAll[i]->SetTitle("");//("Light collection efficiency vs Absorption Length");
		graphsAll[i]->GetXaxis()->SetTitle("Absorption coefficient [cm^{-1}]");
		graphsAll[i]->GetYaxis()->SetTitle("Efficiency");
		//graphsAll[i]->GetXaxis()->SetLimits(0.0, 0.05);
		graphsAll[i]->GetYaxis()->SetRangeUser(1e-5, 2e-1);
		graphsAll[i]->SetMarkerColor(colors[i]);
		graphsAll[i]->SetMarkerSize(0.5);
		graphsAll[i]->SetMarkerStyle(21);

		graphsAllLog.push_back( (TGraphErrors*) graphsAll[i]->Clone() );

		c1->cd(1);
		c1->SetLogy();
		if(i==0) {
			graphsAll[i]->Draw("AP");
		}
		else {
			graphsAll[i]->Draw("P same");
		}
	
		leg->AddEntry(graphsAll[i], namesAll[i].c_str(), "p");
	}

	vector<TF1*> functionV (graphsAll.size());
	for (unsigned int i=0; i < graphsAll.size(); i++) {
		functionV[i] = new TF1("EfficiencyVsMu", EfficiencyVsMu, 0.0, 10.0, 4);
	}
	// Fit all graphs in graphsAll to EfficiencyVsMu.
	for (unsigned int i=0; i < graphsAll.size(); i++) {
		Double_t par[] = {1.0, 1.0, 1.0, -1.0};
		cout << endl << "ieta: " << i+22 << endl; //*-*
		functionV[i]->SetParameters(par);
		//functionV[i]->FixParameter(2, 0.0); //*-* Only fit to exp(linear function).
		functionV[i]->SetLineColor(colors[i]);
		functionV[i]->SetParNames("p0", "p1", "p2", "p3");
		graphsAll[i]->Fit(functionV[i], "", "", 0.01, 0.12);
	}

	// Make log version of graphsAllLog
	for (unsigned int i=0; i < graphsAllLog.size(); i++) {
		//cout << graphsAllLogLog.size() << ": " << graphsAllLog[i] << endl;
		
		int n = graphsAllLog[i]->GetN();
		double* xV = graphsAllLog[i]->GetX();
		double* yV = graphsAllLog[i]->GetY();
		for (int j=0; j < n; j++) {
			yV[j] = log10( yV[j] );
		}
		double* exV = graphsAllLog[i]->GetEX();
		double* eyV = graphsAllLog[i]->GetEY();

		graphsAllLog[i] = new TGraphErrors(n, xV, yV, exV, eyV);
		graphsAllLog[i]->SetTitle("");//("Light collection efficiency vs Absorption Length");
		graphsAllLog[i]->GetXaxis()->SetTitle("Absorption coefficient [cm^{-1}]");
		graphsAllLog[i]->GetYaxis()->SetTitle("Log(Efficiency)");
		//graphsAllLog[i]->GetXaxis()->SetLimits(-0.005, 0.045);
		//graphsAllLog[i]->GetYaxis()->SetRangeUser(1e-6, 1e0);
		graphsAllLog[i]->SetMarkerColor(colors[i]);
		graphsAllLog[i]->SetMarkerSize(0.5);
		graphsAllLog[i]->SetMarkerStyle(21);

		//c1->cd(1);
		//if(i==0) {
		//	graphsAllLog[i]->Draw("AP");
		//}
		//else {
		//	graphsAllLog[i]->Draw("P same");
		//}
	}

	// Draw leg to c1 pad 1.
	{
		c1->cd(1);
		leg->SetFillColor(0);
		leg->SetBorderSize(0);
		leg->SetTextSize(0.03);
		leg->SetTextFont(42);
		leg->Draw();
	}

	vector<TF1*> functionLogV (graphsAllLog.size(), new TF1("LogEfficiencyVsMu", LogEfficiencyVsMu, 0.0, 0.0, 3));
	// Fit all graphs in graphsAllLog to myfunction.
	for (unsigned int i=0; i < graphsAllLog.size(); i++) {
		Double_t par[] = {1.0, 1.0, 1.0};
		functionLogV[i]->SetParameters(par);
		functionLogV[i]->SetLineColor(colors[i]);
		functionLogV[i]->SetParNames("p0", "p1", "p2");
		//*-*graphsAllLog[i]->Fit(functionLogV[i], "", "", 0.01, 0.045);
	}

	// Read laser data file.
	// ratiosMap[ieta] contains ratios for ieta at luminosity values contained in lumis.
	vector<double> lumis;
	std::map<int, vector<double> > ratiosMapData;
	readLaserDataFile("../HELaserData.txt", lumis, ratiosMapData);

	// Read dose and fluence map data files.
	std::map<int, double> doseMapLayer1;
	std::map<int, double> fluenceMapLayer1;
	readDoseFile("../doseLayer1Mars.txt", doseMapLayer1, fluenceMapLayer1); //*-*

	// Variables to read and write from Trees.
	Int_t ieta, lumiIndex;
	Double_t lumi, ratio;
	ieta = 0;
	lumiIndex = 0;
	lumi = 0.;
	ratio = 0.;

	// Prepare to read dataRatios Tree.
	TFile* dataFile = new TFile("dataRatios.root");
	if (!dataFile->IsOpen()) {cout << "dataRatios.root not found"; return;}
	TTree* dataRatios = (TTree*)dataFile->Get("dataRatios");
	dataRatios->Show(10);
	dataRatios->SetBranchAddress("ieta", &ieta);
	dataRatios->SetBranchAddress("lumiIndex", &lumiIndex);
	dataRatios->SetBranchAddress("lumi", &lumi);
	dataRatios->SetBranchAddress("ratio", &ratio);
	
	// Prepare mcRatios Tree to be written. mcRatios shares branch addresses with dataRatios.
	TString mcFileName(runName);
	mcFileName.Prepend("mcRatios-").Append("-").Append(treeLabel).Append(".root");
	TFile* file = new TFile(mcFileName.Data(), "recreate");
	if (!file->IsOpen()) {cout << "mcRatios.root not found"; return;}
	TTree* mcRatios = new TTree("mcRatios", "ratios from Geant4");
	mcRatios->Branch("ieta", &ieta, "ieta/I");
	mcRatios->Branch("lumiIndex", &lumiIndex, "lumiIndex/I");
	mcRatios->Branch("lumi", &lumi, "lumi/D");
	mcRatios->Branch("ratio", &ratio, "ratio/D");

	for (unsigned int i=0; i < graphsAll.size(); i++) {
		// Define parameters for DoseToMu and calculate undamaged efficiency.
		//double par[2] = {0.02, 0.002};
		double doseRateBase = 2.67e-3;
		double A = 0.5;
		double par[2] = {0.023, 0.001};
		double par2[4] = {0.023, 0.001, doseRateBase, A};
		double muTileUndamaged = 0.023;

		for (unsigned int j=0; j < lumis.size(); j++) {
			dataRatios->GetEntryWithIndex(i+22, j); // Get entry with arbitrary ieta and lumiIndex == j.
			//stringstream cutFormula; cutFormula << "ieta == " << i << " && lumiIndex == " << j;
			//cout << cutFormula.str() << endl;
			//dataRatios->Draw("ieta:lumi:ratio:lumiIndex", cutFormula.str().c_str(), "goff");
			//ieta = (dataRatios->GetVal(1))[0];

			//cout << "ieta: " << ieta << endl; 
			//cout << "lumiIndex: " << lumiIndex << endl; 
			//cout << "lumi: " << lumi << endl; 
			//cout << "ratio: " << ratio << endl; 

			double dose = doseFactor * (doseMapLayer1[ieta] / 100) * lumi; // Dose[kGy] per fb-1.
			double doseRate = doseMapLayer1[ieta] * 3.6e-5; // Dose rate [kGy/hr]
			double xx[2] = {dose, doseRate};
			double muTile = DoseToMu(&dose, par);
			//double muTile = CalculateMu(xx, par2);
			//double efficiency = functionV[i]->Eval(muTile);
			double efficiency = graphsAll[i]->Eval(muTile); //*-*
			if(j==0) {
				muTileUndamaged = muTile;
			}
			//double efficiencyUndamaged = functionV[i]->Eval(muTileUndamaged);
			double efficiencyUndamaged = graphsAll[i]->Eval(muTileUndamaged); //*-*
			ratio = efficiency / efficiencyUndamaged;
			double ratioFromData = (ratiosMapData[ieta])[j];

//			if(j == lumis.size()-1) {
//				cout << "ieta: " << i+22 << endl;
//				cout << "dose: " << dose << endl; 
//				cout << "muTile: " << muTile << endl; 
//				cout << "efficiency: " << efficiency << endl; 
//				cout << "efficiencyUndamaged: " << efficiencyUndamaged << endl; 
//				cout << "ratio: " << ratio << endl; 
//				cout << "ratioFromData: " << ratioFromData << endl; 
//				cout << endl;
//			}
			mcRatios->Fill();
		}
	}

	int nEntries = dataRatios->GetEntries();
	for(int i=0; i < nEntries; i++) {
//		dataRatios->GetEntry(i);
//
//		cout << "ieta: " << ieta << endl; 
//		cout << "lumiIndex: " << lumiIndex << endl; 
//		cout << "lumi: " << lumi << endl; 
//		cout << "ratio: " << ratio << endl; 
	}
	
	mcRatios->Write();

}

//Test
void test2(const char* runName) 
{
	// Vector of colors.
	vector<Color_t> colors(18, kBlack);
	{
		colors[0] = kBlue;
		colors[1] = kOrange;
		colors[2] = kGreen;
		colors[3] = kRed;
		colors[4] = kCyan;
		colors[5] = kMagenta;
		// Define colors[6] through colors[11]
		for (int i = 0; i < 6; i++ ) {
			colors[i+6] = colors[i] + 1;
		}
		// Define colors[12] through colors[17]
		for (int i = 0; i < 6; i++ ) {
			colors[i+12] = colors[i] + 2;
		}
	}
	
	TCanvas* c1;
	c1 = new TCanvas("c1", "Efficiency vs Mu", 800, 600);
	//gPad->SetLogx();
	//gPad->SetLogy();
	
	// graphsAll contains Efficiency vs Mu data for all tiles from the Geant4 simulations.
	// namesAll contains names for all graphs in graphsAll.
	// graphsTypeX contains Efficiency vs Lambda data for tiles of type X from the Geant4 simulations.
	vector<TGraphErrors*> graphsLayer1, graphsLayer7, graphsAll, graphsAllLog;
	vector<string> namesAll;

	// Read from Analysis files and add to graphsTypeX, and graphsAll.
	// (Only reading type 1 tiles.)
	{
		graphsLayer1.push_back(readAnalysisFile(constructFileName(runName, 1, 22)));
		graphsLayer1.push_back(readAnalysisFile(constructFileName(runName, 1, 23)));
		graphsLayer1.push_back(readAnalysisFile(constructFileName(runName, 1, 24)));
		graphsLayer1.push_back(readAnalysisFile(constructFileName(runName, 1, 25)));
		graphsLayer1.push_back(readAnalysisFile(constructFileName(runName, 1, 26)));
		graphsLayer1.push_back(readAnalysisFile(constructFileName(runName, 1, 27)));
		graphsLayer1.push_back(readAnalysisFile(constructFileName(runName, 1, 28)));
		graphsLayer1.push_back(readAnalysisFile(constructFileName(runName, 1, 29)));
		
	//	graphsType2.push_back(readAnalysisFile("../run5/Analysis-run5_layer1_ieta17.txt"));
	//	graphsType2.push_back(readAnalysisFile("../run5/Analysis-run5_layer1_ieta18.txt"));
	//	graphsType2.push_back(readAnalysisFile("../run5/Analysis-run5_layer1_ieta19.txt"));
	//	graphsType2.push_back(readAnalysisFile("../run5/Analysis-run5_layer1_ieta20.txt"));
	//	
	//	graphsType3.push_back(readAnalysisFile("../run5/Analysis-run5_layer1_ieta17.txt"));
	//	graphsType3.push_back(readAnalysisFile("../run5/Analysis-run5_layer1_ieta18.txt"));
	//	graphsType3.push_back(readAnalysisFile("../run5/Analysis-run5_layer1_ieta19.txt"));
	//	graphsType3.push_back(readAnalysisFile("../run5/Analysis-run5_layer1_ieta20.txt"));
		
		// Add all graphs from graphsLayer1 to graphsAll. Add all legend names to namesAll.
		for (unsigned int i=0; i < graphsLayer1.size(); i++) {
			graphsAll.push_back(graphsLayer1[i]);
			stringstream name;
			name << "ieta " << (i+22) << "Layer1";
			namesAll.push_back(name.str());
		}
	}

	TLegend* leg = new TLegend(0.7, 0.6, 0.9, 0.85);
	// Run through all graphs in graphsAll and modify error bars if necessary.
	// Push log version of graphsAll into graphsAllLog.
	// Set graph properties of all graphs in graphsAll and draw to c1 pad 1.
	// Add all legend names from namesAll to leg.
	for (unsigned int i=0; i < graphsAll.size(); i++) {
		//cout << graphsAll.size() << ": " << graphsAll[i] << endl;
		
		int n = graphsAll[i]->GetN();
		double* xV = graphsAll[i]->GetX();
		// Offset x-values for clarity.
		// (Disabled)
		for (int j=0; j < n; j++) {
			//xV[j] += i*0.0002;
		}
		double* yV = graphsAll[i]->GetY();
		double* exV = graphsAll[i]->GetEX();
		double* eyV = graphsAll[i]->GetEY();
	
		graphsAll[i] = new TGraphErrors(n, xV, yV, exV, eyV);
		graphsAll[i]->SetTitle("");//("Light collection efficiency vs Absorption Length");
		graphsAll[i]->GetXaxis()->SetTitle("Absorption coefficient [cm^{-1}]");
		graphsAll[i]->GetYaxis()->SetTitle("Efficiency");
		//graphsAll[i]->GetXaxis()->SetLimits(0.0, 0.05);
		graphsAll[i]->GetYaxis()->SetRangeUser(1e-5, 2e-1);
		graphsAll[i]->SetMarkerColor(colors[i]);
		graphsAll[i]->SetMarkerSize(0.5);
		graphsAll[i]->SetMarkerStyle(21);

		graphsAllLog.push_back( (TGraphErrors*) graphsAll[i]->Clone() );

		c1->cd(1);
		c1->SetLogy();
		if(i==0) {
			graphsAll[i]->Draw("AP");
		}
		else {
			graphsAll[i]->Draw("P same");
		}
	
		leg->AddEntry(graphsAll[i], namesAll[i].c_str(), "p");
	}

	vector<TF1*> functionV (graphsAll.size());
	for (unsigned int i=0; i < graphsAll.size(); i++) {
		functionV[i] = new TF1("EfficiencyVsMu", EfficiencyVsMu, 0.0, 10.0, 4);
	}
	for (unsigned int i=0; i==0; i++) {
		double x[1000], y[1000];
		for (int j = 0; j<1000; j++) {
			x[j] = 0.0 + (0.15 - 0.0) * j/1000;
			y[j] = graphsAll[i]->Eval(x[j], 0, "S");
		}
		TGraph* fit = new TGraph(1000, x, y);
		fit->Draw("C same");
	}

	// Make log version of graphsAllLog
	for (unsigned int i=0; i < graphsAllLog.size(); i++) {
		//cout << graphsAllLogLog.size() << ": " << graphsAllLog[i] << endl;
		
		int n = graphsAllLog[i]->GetN();
		double* xV = graphsAllLog[i]->GetX();
		double* yV = graphsAllLog[i]->GetY();
		for (int j=0; j < n; j++) {
			yV[j] = log10( yV[j] );
		}
		double* exV = graphsAllLog[i]->GetEX();
		double* eyV = graphsAllLog[i]->GetEY();

		graphsAllLog[i] = new TGraphErrors(n, xV, yV, exV, eyV);
		graphsAllLog[i]->SetTitle("");//("Light collection efficiency vs Absorption Length");
		graphsAllLog[i]->GetXaxis()->SetTitle("Absorption coefficient [cm^{-1}]");
		graphsAllLog[i]->GetYaxis()->SetTitle("Log(Efficiency)");
		//graphsAllLog[i]->GetXaxis()->SetLimits(-0.005, 0.045);
		//graphsAllLog[i]->GetYaxis()->SetRangeUser(1e-6, 1e0);
		graphsAllLog[i]->SetMarkerColor(colors[i]);
		graphsAllLog[i]->SetMarkerSize(0.5);
		graphsAllLog[i]->SetMarkerStyle(21);

		//c1->cd(1);
		//if(i==0) {
		//	graphsAllLog[i]->Draw("AP");
		//}
		//else {
		//	graphsAllLog[i]->Draw("P same");
		//}
	}

	// Draw leg to c1 pad 1.
	{
		c1->cd(1);
		leg->SetFillColor(0);
		leg->SetBorderSize(0);
		leg->SetTextSize(0.03);
		leg->SetTextFont(42);
		leg->Draw();
	}

	vector<TF1*> functionLogV (graphsAllLog.size(), new TF1("LogEfficiencyVsMu", LogEfficiencyVsMu, 0.0, 0.0, 3));
	// Fit all graphs in graphsAllLog to myfunction.
	for (unsigned int i=0; i < graphsAllLog.size(); i++) {
		Double_t par[] = {1.0, 1.0, 1.0};
		functionLogV[i]->SetParameters(par);
		functionLogV[i]->SetLineColor(colors[i]);
		functionLogV[i]->SetParNames("p0", "p1", "p2");
		//*-*graphsAllLog[i]->Fit(functionLogV[i], "", "", 0.01, 0.045);
	}

}

void readRatios(string inputFile)
{
	TFile* file = new TFile("dataRatios.root", "recreate");
	TTree* dataRatios = new TTree("dataRatios", "ratios from laser data");
	Int_t ieta, lumiIndex;
	lumiIndex = 0;
	Double_t lumi, ratio;
	dataRatios->Branch("ieta", &ieta, "ieta/I");
	dataRatios->Branch("lumiIndex", &lumiIndex, "lumiIndex/I");
	dataRatios->Branch("lumi", &lumi, "lumi/D");
	dataRatios->Branch("ratio", &ratio, "ratio/D");

	ifstream inputStream(inputFile.c_str());
	string line;
	if(inputStream.is_open())
	{
		while(getline(inputStream, line))
		{
			if(line[0]=='#') continue; //Skip comment lines
			else {
				// Process the line to get ieta and dose.
				//double lumi;
				double ratios[13];
				processLaserDataLine(line, lumi, ratios);
				for (int i=0; i < 13; i++) {
					ieta = i + 17;
					ratio = ratios[i];

					cout << "ieta: " << ieta << endl; 
					cout << "lumiIndex: " << lumiIndex << endl; 
					cout << "lumi: " << lumi << endl; 
					cout << "ratio: " << ratio << endl << endl; 

					dataRatios->Fill();
					//cout << "ratios[i] " <<ratios[i] << endl;
				}
				lumiIndex = lumiIndex +1; 
			}
		}
	}
	else {cout <<"Error: Input file not found." <<endl;}
	dataRatios->BuildIndex("ieta", "lumiIndex");
	dataRatios->Write();
}

void ratiosPlot() //bool doPlot = false)
{
	readRatios("../HELaserData.txt");

//	TCanvas* c1 = new TCanvas("c1","canvas1",800,600);
//	TGraph* graph;
//	graph = new TGraph(lumis.size(), &lumis[0], &((ratiosMap[ieta])[0]));
//	c1->cd();
//	graph->Draw("AP");
}

void test(const char* input)
{
	TString* fileName;
	int layer = 1;
	fileName = new TString("../RUNNAME/Analysis-RUNNAME_layerLAYERNUMBER_ietaIETA.txt");
	fileName->ReplaceAll("RUNNAME", input);
	fileName->Append(TString::Itoa(layer,10));
	cout << fileName->Data() <<endl;
	readAnalysisFile("../run5bc/Analysis-run5bc_layer1_ieta22.txt");
}

// Copy of MakeTrees to do enable fitting of dose factor.
void FitEffVsMu(const char* runName, double doseFactor)//, const char* treeLabel) 
{
	// Vector of colors.
	vector<Color_t> colors(18, kBlack);
	{
		colors[0] = kBlue;
		colors[1] = kOrange;
		colors[2] = kGreen;
		colors[3] = kRed;
		colors[4] = kCyan;
		colors[5] = kMagenta;
		// Define colors[6] through colors[11]
		for (int i = 0; i < 6; i++ ) {
			colors[i+6] = colors[i] + 1;
		}
		// Define colors[12] through colors[17]
		for (int i = 0; i < 6; i++ ) {
			colors[i+12] = colors[i] + 2;
		}
	}
	
	TCanvas* c1;
	c1 = new TCanvas("c1", "Efficiency vs Mu", 800, 600);
	//gPad->SetLogx();
	//gPad->SetLogy();
	
	// graphsAll contains Efficiency vs Mu data for all tiles from the Geant4 simulations.
	// namesAll contains names for all graphs in graphsAll.
	// graphsTypeX contains Efficiency vs Lambda data for tiles of type X from the Geant4 simulations.
	vector<TGraphErrors*> graphsLayer1, graphsLayer7, graphsAll, graphsAllLog;
	vector<string> namesAll;

	// Read from Analysis files and add to graphsTypeX, and graphsAll.
	// (Only reading type 1 tiles.)
	{
		graphsLayer1.push_back(readAnalysisFile(constructFileName(runName, 1, 22)));
		graphsLayer1.push_back(readAnalysisFile(constructFileName(runName, 1, 23)));
		graphsLayer1.push_back(readAnalysisFile(constructFileName(runName, 1, 24)));
		graphsLayer1.push_back(readAnalysisFile(constructFileName(runName, 1, 25)));
		graphsLayer1.push_back(readAnalysisFile(constructFileName(runName, 1, 26)));
		graphsLayer1.push_back(readAnalysisFile(constructFileName(runName, 1, 27)));
		graphsLayer1.push_back(readAnalysisFile(constructFileName(runName, 1, 28)));
		graphsLayer1.push_back(readAnalysisFile(constructFileName(runName, 1, 29)));

		graphsLayer1.push_back(readAnalysisFile(constructFileName(runName, 7, 22)));
		graphsLayer1.push_back(readAnalysisFile(constructFileName(runName, 7, 23)));
		graphsLayer1.push_back(readAnalysisFile(constructFileName(runName, 7, 24)));
		graphsLayer1.push_back(readAnalysisFile(constructFileName(runName, 7, 25)));
		graphsLayer1.push_back(readAnalysisFile(constructFileName(runName, 7, 26)));
		graphsLayer1.push_back(readAnalysisFile(constructFileName(runName, 7, 27)));
		graphsLayer1.push_back(readAnalysisFile(constructFileName(runName, 7, 28)));
		graphsLayer1.push_back(readAnalysisFile(constructFileName(runName, 7, 29)));
		
	//	graphsType2.push_back(readAnalysisFile("../run5/Analysis-run5_layer1_ieta17.txt"));
	//	graphsType2.push_back(readAnalysisFile("../run5/Analysis-run5_layer1_ieta18.txt"));
	//	graphsType2.push_back(readAnalysisFile("../run5/Analysis-run5_layer1_ieta19.txt"));
	//	graphsType2.push_back(readAnalysisFile("../run5/Analysis-run5_layer1_ieta20.txt"));
	//	
	//	graphsType3.push_back(readAnalysisFile("../run5/Analysis-run5_layer1_ieta17.txt"));
	//	graphsType3.push_back(readAnalysisFile("../run5/Analysis-run5_layer1_ieta18.txt"));
	//	graphsType3.push_back(readAnalysisFile("../run5/Analysis-run5_layer1_ieta19.txt"));
	//	graphsType3.push_back(readAnalysisFile("../run5/Analysis-run5_layer1_ieta20.txt"));
		
		// Add all graphs from graphsLayer1 to graphsAll. Add all legend names to namesAll.
		for (unsigned int i=0; i < graphsLayer1.size(); i++) {
			graphsAll.push_back(graphsLayer1[i]);
			stringstream name;
			name << "ieta " << (i+22) << "Layer1";
			namesAll.push_back(name.str());
		}
	}

	TLegend* leg = new TLegend(0.7, 0.6, 0.9, 0.85);
	// Run through all graphs in graphsAll and modify error bars if necessary.
	// Push log version of graphsAll into graphsAllLog.
	// Set graph properties of all graphs in graphsAll and draw to c1 pad 1.
	// Add all legend names from namesAll to leg.
	for (unsigned int i=0; i < graphsAll.size(); i++) {
		//cout << graphsAll.size() << ": " << graphsAll[i] << endl;
		
		int n = graphsAll[i]->GetN();
		double* xV = graphsAll[i]->GetX();
		// Offset x-values for clarity.
		// (Disabled)
		for (int j=0; j < n; j++) {
			//xV[j] += i*0.0002;
		}
		double* yV = graphsAll[i]->GetY();
		double* exV = graphsAll[i]->GetEX();
		double* eyV = graphsAll[i]->GetEY();
	
		graphsAll[i] = new TGraphErrors(n, xV, yV, exV, eyV);
		graphsAll[i]->SetTitle("");//("Light collection efficiency vs Absorption Length");
		graphsAll[i]->GetXaxis()->SetTitle("Absorption coefficient [cm^{-1}]");
		graphsAll[i]->GetYaxis()->SetTitle("Efficiency");
		//graphsAll[i]->GetXaxis()->SetLimits(0.0, 0.05);
		graphsAll[i]->GetYaxis()->SetRangeUser(1e-5, 2e-1);
		graphsAll[i]->SetMarkerColor(colors[i]);
		graphsAll[i]->SetMarkerSize(0.5);
		graphsAll[i]->SetMarkerStyle(21);

		graphsAllLog.push_back( (TGraphErrors*) graphsAll[i]->Clone() );

		c1->cd(1);
		c1->SetLogy();
		if(i==0) {
			graphsAll[i]->Draw("AP");
		}
		else {
			graphsAll[i]->Draw("P same");
		}
	
		leg->AddEntry(graphsAll[i], namesAll[i].c_str(), "p");
	}

	ofstream outfile;
	outfile.open("EffVsMuFitParametersRun13.txt");

	vector<TF1*> functionV (graphsAll.size());
	for (unsigned int i=0; i < graphsAll.size(); i++) {
		functionV[i] = new TF1("EfficiencyVsMu", EfficiencyVsMu, 0.0, 10.0, 4);
	}
	// Fit all graphs in graphsAll to EfficiencyVsMu.
	outfile << "layer\tieta\tp1\tp0" << endl;
	for (unsigned int i=0; i < graphsAll.size(); i++) {
		Double_t par[] = {1.0, 1.0, 1.0, -1.0};
		functionV[i]->SetParameters(par);
		functionV[i]->FixParameter(3, 0.0); //*-* Only fit to exp(linear function).
		functionV[i]->FixParameter(2, 0.0); //*-* Only fit to exp(linear function).
		functionV[i]->SetLineColor(colors[i]);
		functionV[i]->SetParNames("p0", "p1", "p2", "p3");
		graphsAll[i]->Fit(functionV[i], "", "", 0.04, 0.12);

		int ieta = (i % 8) + 22;
		int layer = i<8 ? 1 : 7;
		double p1 = functionV[i]->GetParameter("p1");
		double p0 = functionV[i]->GetParameter("p0");
		outfile << layer << "\t" << ieta << "\t";//*-*
		outfile << p1 << "\t" << p0 << endl;
	}

	// Make log version of graphsAllLog
	for (unsigned int i=0; i < graphsAllLog.size(); i++) {
		//cout << graphsAllLogLog.size() << ": " << graphsAllLog[i] << endl;
		
		int n = graphsAllLog[i]->GetN();
		double* xV = graphsAllLog[i]->GetX();
		double* yV = graphsAllLog[i]->GetY();
		for (int j=0; j < n; j++) {
			yV[j] = log10( yV[j] );
		}
		double* exV = graphsAllLog[i]->GetEX();
		double* eyV = graphsAllLog[i]->GetEY();

		graphsAllLog[i] = new TGraphErrors(n, xV, yV, exV, eyV);
		graphsAllLog[i]->SetTitle("");//("Light collection efficiency vs Absorption Length");
		graphsAllLog[i]->GetXaxis()->SetTitle("Absorption coefficient [cm^{-1}]");
		graphsAllLog[i]->GetYaxis()->SetTitle("Log(Efficiency)");
		//graphsAllLog[i]->GetXaxis()->SetLimits(-0.005, 0.045);
		//graphsAllLog[i]->GetYaxis()->SetRangeUser(1e-6, 1e0);
		graphsAllLog[i]->SetMarkerColor(colors[i]);
		graphsAllLog[i]->SetMarkerSize(0.5);
		graphsAllLog[i]->SetMarkerStyle(21);

		//c1->cd(1);
		//if(i==0) {
		//	graphsAllLog[i]->Draw("AP");
		//}
		//else {
		//	graphsAllLog[i]->Draw("P same");
		//}
	}

	// Draw leg to c1 pad 1.
	{
		c1->cd(1);
		leg->SetFillColor(0);
		leg->SetBorderSize(0);
		leg->SetTextSize(0.03);
		leg->SetTextFont(42);
		leg->Draw();
	}

	/*
	vector<TF1*> functionLogV (graphsAllLog.size(), new TF1("LogEfficiencyVsMu", LogEfficiencyVsMu, 0.0, 0.0, 3));
	// Fit all graphs in graphsAllLog to myfunction.
	for (unsigned int i=0; i < graphsAllLog.size(); i++) {
		Double_t par[] = {1.0, 1.0, 1.0};
		functionLogV[i]->SetParameters(par);
		functionLogV[i]->SetLineColor(colors[i]);
		functionLogV[i]->SetParNames("p0", "p1", "p2");
		//*-*graphsAllLog[i]->Fit(functionLogV[i], "", "", 0.01, 0.045);
	}

	// Read laser data file.
	// ratiosMap[ieta] contains ratios for ieta at luminosity values contained in lumis.
	vector<double> lumis;
	std::map<int, vector<double> > ratiosMapData;
	readLaserDataFile("../HELaserData.txt", lumis, ratiosMapData);

	// Read dose and fluence map data files.
	std::map<int, double> doseMapLayer1;
	std::map<int, double> fluenceMapLayer1;
	readDoseFile("../doseLayer1Mars.txt", doseMapLayer1, fluenceMapLayer1); //*-*

	// Variables to read and write from Trees.
	Int_t ieta, lumiIndex;
	Double_t lumi, ratio;
	ieta = 0;
	lumiIndex = 0;
	lumi = 0.;
	ratio = 0.;

	// Prepare to read dataRatios Tree.
	TFile* dataFile = new TFile("dataRatios.root");
	if (!dataFile->IsOpen()) {cout << "dataRatios.root not found"; return;}
	TTree* dataRatios = (TTree*)dataFile->Get("dataRatios");
	dataRatios->Show(10);
	dataRatios->SetBranchAddress("ieta", &ieta);
	dataRatios->SetBranchAddress("lumiIndex", &lumiIndex);
	dataRatios->SetBranchAddress("lumi", &lumi);
	dataRatios->SetBranchAddress("ratio", &ratio);
	
	// Prepare mcRatios Tree to be written. mcRatios shares branch addresses with dataRatios.
	TString mcFileName(runName);
	mcFileName.Prepend("mcRatios-").Append("-").Append(treeLabel).Append(".root");
	TFile* file = new TFile(mcFileName.Data(), "recreate");
	if (!file->IsOpen()) {cout << "mcRatios.root not found"; return;}
	TTree* mcRatios = new TTree("mcRatios", "ratios from Geant4");
	mcRatios->Branch("ieta", &ieta, "ieta/I");
	mcRatios->Branch("lumiIndex", &lumiIndex, "lumiIndex/I");
	mcRatios->Branch("lumi", &lumi, "lumi/D");
	mcRatios->Branch("ratio", &ratio, "ratio/D");

	for (unsigned int i=0; i < graphsAll.size(); i++) {
		// Define parameters for DoseToMu and calculate undamaged efficiency.
		//double par[2] = {0.02, 0.002};
		double doseRateBase = 2.67e-3;
		double A = 0.5;
		double par[2] = {0.023, 0.001};
		double par2[4] = {0.023, 0.001, doseRateBase, A};
		double muTileUndamaged = 0.023;

		for (unsigned int j=0; j < lumis.size(); j++) {
			dataRatios->GetEntryWithIndex(i+22, j); // Get entry with arbitrary ieta and lumiIndex == j.
			//stringstream cutFormula; cutFormula << "ieta == " << i << " && lumiIndex == " << j;
			//cout << cutFormula.str() << endl;
			//dataRatios->Draw("ieta:lumi:ratio:lumiIndex", cutFormula.str().c_str(), "goff");
			//ieta = (dataRatios->GetVal(1))[0];

			//cout << "ieta: " << ieta << endl; 
			//cout << "lumiIndex: " << lumiIndex << endl; 
			//cout << "lumi: " << lumi << endl; 
			//cout << "ratio: " << ratio << endl; 

			double dose = doseFactor * (doseMapLayer1[ieta] / 100) * lumi; // Dose[kGy] per fb-1.
			double doseRate = doseMapLayer1[ieta] * 3.6e-5; // Dose rate [kGy/hr]
			double xx[2] = {dose, doseRate};
			double muTile = DoseToMu(&dose, par);
			//double muTile = CalculateMu(xx, par2);
			//double efficiency = functionV[i]->Eval(muTile);
			double efficiency = graphsAll[i]->Eval(muTile); //*-*
			if(j==0) {
				muTileUndamaged = muTile;
			}
			//double efficiencyUndamaged = functionV[i]->Eval(muTileUndamaged);
			double efficiencyUndamaged = graphsAll[i]->Eval(muTileUndamaged); //*-*
			ratio = efficiency / efficiencyUndamaged;
			double ratioFromData = (ratiosMapData[ieta])[j];

//			if(j == lumis.size()-1) {
//				cout << "ieta: " << i+22 << endl;
//				cout << "dose: " << dose << endl; 
//				cout << "muTile: " << muTile << endl; 
//				cout << "efficiency: " << efficiency << endl; 
//				cout << "efficiencyUndamaged: " << efficiencyUndamaged << endl; 
//				cout << "ratio: " << ratio << endl; 
//				cout << "ratioFromData: " << ratioFromData << endl; 
//				cout << endl;
//			}
			mcRatios->Fill();
		}
	}

	int nEntries = dataRatios->GetEntries();
	for(int i=0; i < nEntries; i++) {
//		dataRatios->GetEntry(i);
//
//		cout << "ieta: " << ieta << endl; 
//		cout << "lumiIndex: " << lumiIndex << endl; 
//		cout << "lumi: " << lumi << endl; 
//		cout << "ratio: " << ratio << endl; 
	}
	
	mcRatios->Write();
*/
}

