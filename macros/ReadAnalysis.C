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

using TMath::Exp;
using TMath::Power;
using TMath::Log;

void processLine(string line, double& absLength, double& efficiency)
{
	stringstream ss(line);
	string field;
	int ctr = 0; //column index of tab-delimited row
	while(getline(ss, field, '\t'))
	{
		if(ctr==0) {stringstream temp(field); temp >> absLength;}
		else if (ctr==1) {stringstream temp(field); temp >> efficiency;}
		else {}
		
		ctr++;
	}
}

TGraphErrors* readAnalysisFile(string inputFile)
{
	vector<double> absLengthV;
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
				double absLength;
				double efficiency;
				double errorEfficiency;
				processLine(line, absLength, efficiency);
				absLengthV.push_back(absLength);
				efficiencyV.push_back(efficiency);
				// Binomial error from 1e6 events.
				errorEfficiency = sqrt(efficiency * (1 - efficiency) / 1e6);
				errorV.push_back(errorEfficiency);
			}
		}
		Int_t n = absLengthV.size();
		gr = new TGraphErrors(n, &absLengthV[0], &efficiencyV[0], NULL, &errorV[0]);
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

void longrun6()
{
	TCanvas* c1;
	// Vector of colors.
	vector<Color_t> colors(18, kBlack);
	colors[0] = kBlue;
	colors[1] = kOrange;
	colors[2] = kGreen;
	colors[3] = kRed;
	colors[4] = kCyan;
	colors[5] = kMagenta;
	for (int i = 0; i < 6; i++ ) {
		colors[i+6] = colors[i] + 1;
	}
	for (int i = 0; i < 6; i++ ) {
		colors[i+12] = colors[i] + 2;
	}
	
	
	c1 = new TCanvas("c1","canvas1",800,600);
	c1->SetLogx();
	c1->SetLogy();
	vector<TGraphErrors*> graphsType1, graphsType2, graphsType3, graphsAll;
	vector<string> namesAll;
	TLegend* leg = new TLegend(0.7, 0.2, 0.9, 0.6);

	graphsType1.push_back(readAnalysisFile("longrun6_ieta21_type1.txt"));
	graphsType1.push_back(readAnalysisFile("longrun6_ieta22_type1.txt"));
	graphsType1.push_back(readAnalysisFile("longrun6_ieta23_type1.txt"));
	graphsType1.push_back(readAnalysisFile("longrun6_ieta24_type1.txt"));
	graphsType1.push_back(readAnalysisFile("longrun6_ieta25_type1.txt"));
	graphsType1.push_back(readAnalysisFile("longrun6_ieta26_type1.txt"));
	graphsType1.push_back(readAnalysisFile("longrun6_ieta27_type1.txt"));
	graphsType1.push_back(readAnalysisFile("longrun6_ieta28_type1.txt"));
	graphsType1.push_back(readAnalysisFile("longrun6_ieta29_type1.txt"));
	
	graphsType2.push_back(readAnalysisFile("longrun6_ieta17_type2.txt"));
	graphsType2.push_back(readAnalysisFile("longrun6_ieta18_type2.txt"));
	graphsType2.push_back(readAnalysisFile("longrun6_ieta19_type2.txt"));
	graphsType2.push_back(readAnalysisFile("longrun6_ieta20_type2.txt"));
	
	graphsType3.push_back(readAnalysisFile("longrun6_ieta17_type3.txt"));
	graphsType3.push_back(readAnalysisFile("longrun6_ieta18_type3.txt"));
	graphsType3.push_back(readAnalysisFile("longrun6_ieta19_type3.txt"));
	graphsType3.push_back(readAnalysisFile("longrun6_ieta20_type3.txt"));
	
	// Add all graphs to graphsAll. Add all legend names to namesAll;
	for (unsigned int i=0; i < graphsType3.size(); i++) {
		graphsAll.push_back(graphsType3[i]);
		stringstream name;
		name << "ieta " << (i+17) << " type3";
		namesAll.push_back(name.str());
	}
	for (unsigned int i = 0; i < graphsType2.size(); i++) {
		graphsAll.push_back(graphsType2[i]);
		stringstream name;
		name << "ieta " << (i+17) << " type2";
		namesAll.push_back(name.str());
	}
	for (unsigned int i = 0; i < graphsType1.size(); i++) {
		graphsAll.push_back(graphsType1[i]);
		stringstream name;
		name << "ieta " << (i+21) << " type1";
		namesAll.push_back(name.str());
	}


	for (unsigned int i=0; i < graphsAll.size(); i++) {
		//cout << graphsAll.size() << ": " << graphsAll[i] << endl;
		
		// Offset x-values for clarity.
		int n = graphsAll[i]->GetN();
		
		double* xV = graphsAll[i]->GetX();
		for (int j=0; j < n; j++) {
			xV[j] *= pow(1.02, i);
		}
		double* yV = graphsAll[i]->GetY();
		double* exV = graphsAll[i]->GetEX();
		double* eyV = graphsAll[i]->GetEY();
		
		graphsAll[i] = new TGraphErrors(n, xV, yV, exV, eyV);
		graphsAll[i]->SetTitle("");//("Light collection efficiency vs Absorption Length");

		graphsAll[i]->GetXaxis()->SetTitle("Absorption length (cm)");
		graphsAll[i]->GetYaxis()->SetTitle("Efficiency");
		graphsAll[i]->GetXaxis()->SetLimits(0.1,200);
		graphsAll[i]->GetYaxis()->SetRangeUser(1e-7,3e-1);
		graphsAll[i]->SetMarkerColor(colors[i]);
		graphsAll[i]->SetMarkerSize(0.5);
		graphsAll[i]->SetMarkerStyle(21);
		if(i==0) {
			graphsAll[i]->Draw("AP");
		}
		else {
			graphsAll[i]->Draw("P same");
		}
		
		leg->AddEntry(graphsAll[i], namesAll[i].c_str(), "p");
	}

	leg->SetFillColor(0);
	leg->SetBorderSize(0);
	leg->SetTextSize(0.03);
	leg->SetTextFont(42);
	leg->Draw();
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

TGraphErrors* longrun6test()
{
	TCanvas* c1;
	// Vector of colors.
	vector<Color_t> colors(18, kBlack);
	colors[0] = kBlue;
	colors[1] = kOrange;
	colors[2] = kGreen;
	colors[3] = kRed;
	colors[4] = kCyan;
	colors[5] = kMagenta;
	for (int i = 0; i < 6; i++ ) {
		colors[i+6] = colors[i] + 1;
	}
	for (int i = 0; i < 6; i++ ) {
		colors[i+12] = colors[i] + 2;
	}
	
	
	c1 = new TCanvas("c1","canvas1",800,600);
	c1->SetLogx();
	c1->SetLogy();
	vector<TGraphErrors*> graphsType1, graphsType2, graphsType3, graphsAll;
	vector<string> namesAll;
	TLegend* leg = new TLegend(0.7, 0.2, 0.9, 0.6);

	graphsType1.push_back(readAnalysisFile("longrun6_ieta21_type1.txt"));
	graphsType1.push_back(readAnalysisFile("longrun6_ieta22_type1.txt"));
	graphsType1.push_back(readAnalysisFile("longrun6_ieta23_type1.txt"));
	graphsType1.push_back(readAnalysisFile("longrun6_ieta24_type1.txt"));
	graphsType1.push_back(readAnalysisFile("longrun6_ieta25_type1.txt"));
	graphsType1.push_back(readAnalysisFile("longrun6_ieta26_type1.txt"));
	graphsType1.push_back(readAnalysisFile("longrun6_ieta27_type1.txt"));
	graphsType1.push_back(readAnalysisFile("longrun6_ieta28_type1.txt"));
	graphsType1.push_back(readAnalysisFile("longrun6_ieta29_type1.txt"));
	
//	graphsType2.push_back(readAnalysisFile("longrun6_ieta17_type2.txt"));
//	graphsType2.push_back(readAnalysisFile("longrun6_ieta18_type2.txt"));
//	graphsType2.push_back(readAnalysisFile("longrun6_ieta19_type2.txt"));
//	graphsType2.push_back(readAnalysisFile("longrun6_ieta20_type2.txt"));
//	
//	graphsType3.push_back(readAnalysisFile("longrun6_ieta17_type3.txt"));
//	graphsType3.push_back(readAnalysisFile("longrun6_ieta18_type3.txt"));
//	graphsType3.push_back(readAnalysisFile("longrun6_ieta19_type3.txt"));
//	graphsType3.push_back(readAnalysisFile("longrun6_ieta20_type3.txt"));
	
	// Add all graphs to graphsAll. Add all legend names to namesAll;
	for (unsigned int i=0; i < graphsType3.size(); i++) {
		graphsAll.push_back(graphsType3[i]);
		stringstream name;
		name << "ieta " << (i+17) << " type3";
		namesAll.push_back(name.str());
	}
	for (unsigned int i = 0; i < graphsType2.size(); i++) {
		graphsAll.push_back(graphsType2[i]);
		stringstream name;
		name << "ieta " << (i+17) << " type2";
		namesAll.push_back(name.str());
	}
	for (unsigned int i = 0; i < graphsType1.size(); i++) {
		graphsAll.push_back(graphsType1[i]);
		stringstream name;
		name << "ieta " << (i+21) << " type1";
		namesAll.push_back(name.str());
	}

	vector<TF1*> functionV (graphsAll.size(), new TF1("myfunc", myfunction, 0.1, 100.0, 2));

	for (unsigned int i=0; i < graphsAll.size(); i++) {
		//cout << graphsAll.size() << ": " << graphsAll[i] << endl;
		
		// Offset x-values for clarity.
		int n = graphsAll[i]->GetN();
		
		double* xV = graphsAll[i]->GetX();
		for (int j=0; j < n; j++) {
			xV[j] *= pow(1.02, i);
		}
		double* yV = graphsAll[i]->GetY();
		double* exV = graphsAll[i]->GetEX();
		double* eyV = graphsAll[i]->GetEY();
	
		graphsAll[i] = new TGraphErrors(n, xV, yV, exV, eyV);
		graphsAll[i]->SetTitle("");//("Light collection efficiency vs Absorption Length");

		graphsAll[i]->GetXaxis()->SetTitle("Absorption length (cm)");
		graphsAll[i]->GetYaxis()->SetTitle("Efficiency");
		graphsAll[i]->GetXaxis()->SetLimits(0.3,200);
		graphsAll[i]->GetYaxis()->SetRangeUser(1e-6,3e-1);
		graphsAll[i]->SetMarkerColor(colors[i]);
		graphsAll[i]->SetMarkerSize(0.5);
		graphsAll[i]->SetMarkerStyle(21);
		if(i==0) {
			graphsAll[i]->Draw("AP");
		}
		else {
			graphsAll[i]->Draw("P same");
		}
	
		leg->AddEntry(graphsAll[i], namesAll[i].c_str(), "p");
		
		Double_t par[] = {1.0, 0.5};
		functionV[i]->SetParameters(par);
		functionV[i]->SetLineColor(colors[i]);
		functionV[i]->SetParNames("average track length", "constant", "offset");
		cout << "Fit number: " << i << endl;
		graphsAll[i]->Fit(functionV[i], "", "", 0.1, 200);
	}
	leg->SetFillColor(0);
	leg->SetBorderSize(0);
	leg->SetTextSize(0.03);
	leg->SetTextFont(42);
	leg->Draw();
	return graphsAll[0];
}

void myfunc()
{
	TF1 *f1 = new TF1("myfunc", myfunction, 0.1, 100.0, 2);
	Double_t par[] = {1.0, 1.0};
	f1->SetParameters(par);
	f1->SetParNames("average track length", "constant");
	TGraphErrors* gr = longrun6test();
	gr->GetXaxis()->SetTitle("Absorption length (cm)");
	gr->GetYaxis()->SetTitle("Efficiency");
	gr->GetXaxis()->SetLimits(0.1,300);
	gr->GetYaxis()->SetRangeUser(1e-7,1e-0);
	gr->Draw("AP");
	gr->Fit(f1, "", "", 0.5, 100);
	f1->GetParameters(par);
	cout << par[0] << "\t" << par[1] << endl;
}

void processLaserDataLine(string line, double& lumi, double* ratios)
{
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

void ratiosMap(int ieta = 29) //bool doPlot = false)
{
	vector<double> lumis;
	std::map<int, vector<double> > ratiosMap;

	readLaserDataFile("HELaserData.txt", lumis, ratiosMap);


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


void MakePlots(int minLumiIndex, int maxLumiIndex, double refAttLength, double lumioffset, bool doPlot = false, bool doSavePlots = false, bool doDoseOrFluence = true) 
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
	c1 = new TCanvas("c1", "Efficiency vs Lambda and Mu", 1200, 600);
	c1->Divide(2,0);
	c1->cd(1);
	gPad->SetLogx();
	gPad->SetLogy();
	c1->cd(2);
	// graphsAll contains Efficiency vs Lambda data for all tiles from the Geant4 simulations.
	// namesAll contains names for all graphs in graphsAll.
	// graphsTypeX contains Efficiency vs Lambda data for tiles of type X from the Geant4 simulations.
	vector<TGraphErrors*> graphsType1, graphsType2, graphsType3, graphsAll, graphsAll2;
	vector<string> namesAll;

	// Read from Analysis files and add to graphsTypeX, and graphsAll.
	// (Only reading type 1 tiles.)
	{
		graphsType1.push_back(readAnalysisFile("longrun6_ieta21_type1.txt"));
		graphsType1.push_back(readAnalysisFile("longrun6_ieta22_type1.txt"));
		graphsType1.push_back(readAnalysisFile("longrun6_ieta23_type1.txt"));
		graphsType1.push_back(readAnalysisFile("longrun6_ieta24_type1.txt"));
		graphsType1.push_back(readAnalysisFile("longrun6_ieta25_type1.txt"));
		graphsType1.push_back(readAnalysisFile("longrun6_ieta26_type1.txt"));
		graphsType1.push_back(readAnalysisFile("longrun6_ieta27_type1.txt"));
		graphsType1.push_back(readAnalysisFile("longrun6_ieta28_type1.txt"));
		graphsType1.push_back(readAnalysisFile("longrun6_ieta29_type1.txt"));
		
	//	graphsType2.push_back(readAnalysisFile("longrun6_ieta17_type2.txt"));
	//	graphsType2.push_back(readAnalysisFile("longrun6_ieta18_type2.txt"));
	//	graphsType2.push_back(readAnalysisFile("longrun6_ieta19_type2.txt"));
	//	graphsType2.push_back(readAnalysisFile("longrun6_ieta20_type2.txt"));
	//	
	//	graphsType3.push_back(readAnalysisFile("longrun6_ieta17_type3.txt"));
	//	graphsType3.push_back(readAnalysisFile("longrun6_ieta18_type3.txt"));
	//	graphsType3.push_back(readAnalysisFile("longrun6_ieta19_type3.txt"));
	//	graphsType3.push_back(readAnalysisFile("longrun6_ieta20_type3.txt"));
		
		// Add all graphs from graphsType3 to graphsAll. Add all legend names to namesAll.
		for (unsigned int i=0; i < graphsType3.size(); i++) {
			graphsAll.push_back(graphsType3[i]);
			stringstream name;
			name << "ieta " << (i+17) << " type3";
			namesAll.push_back(name.str());
		}
		// Add all graphs from graphsType2 to graphsAll. Add all legend names to namesAll;
		for (unsigned int i = 0; i < graphsType2.size(); i++) {
			graphsAll.push_back(graphsType2[i]);
			stringstream name;
			name << "ieta " << (i+17) << " type2";
			namesAll.push_back(name.str());
		}
		// Add all graphs from graphsType1 to graphsAll. Add all legend names to namesAll;
		for (unsigned int i = 0; i < graphsType1.size(); i++) {
			graphsAll.push_back(graphsType1[i]);
			stringstream name;
			name << "ieta " << (i+21) << " type1";
			namesAll.push_back(name.str());
		}
	}

	TLegend* leg = new TLegend(0.7, 0.2, 0.9, 0.6);
	// Run through all graphs in graphsAll and modify error bars if necessary.
	// Set graph properties of all graphs in graphsAll and draw to c1 pad 1.
	// Add all legend names from namesAll to leg.
	for (unsigned int i=0; i < graphsAll.size(); i++) {
		//cout << graphsAll.size() << ": " << graphsAll[i] << endl;
		
		int n = graphsAll[i]->GetN();
		double* xV = graphsAll[i]->GetX();
		// Offset x-values for clarity.
		// (Disabled)
		for (int j=0; j < n; j++) {
			//xV[j] *= pow(1.02, i);
		}

		double* yV = graphsAll[i]->GetY();
		double* exV = graphsAll[i]->GetEX();
		double* eyV = graphsAll[i]->GetEY();
	
		graphsAll[i] = new TGraphErrors(n, xV, yV, exV, eyV);
		graphsAll[i]->SetTitle("");//("Light collection efficiency vs Absorption Length");
		graphsAll[i]->GetXaxis()->SetTitle("Absorption length (cm)");
		graphsAll[i]->GetYaxis()->SetTitle("Efficiency");
		graphsAll[i]->GetXaxis()->SetLimits(4, 200);
		graphsAll[i]->GetYaxis()->SetRangeUser(3e-3, 2e-1);
		graphsAll[i]->SetMarkerColor(colors[i]);
		graphsAll[i]->SetMarkerSize(0.5);
		graphsAll[i]->SetMarkerStyle(21);

		c1->cd(1);
		if(i==0) {
			graphsAll[i]->Draw("AP");
		}
		else {
			graphsAll[i]->Draw("P same");
		}
	
		leg->AddEntry(graphsAll[i], namesAll[i].c_str(), "p");
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

	vector<TF1*> functionV (graphsAll.size(), new TF1("myfunction", myfunction, 0.01, 100.0, 2));
	// Fit all graphs in graphsAll to myfunction.
	for (unsigned int i=0; i < graphsAll.size(); i++) {
		Double_t par[] = {1.0, 0.5};
		functionV[i]->SetParameters(par);
		functionV[i]->SetLineColor(colors[i]);
		functionV[i]->SetParNames("average track length", "constant", "offset");
		graphsAll[i]->Fit(functionV[i], "", "", 0.1, 200);
	}

	vector<TF1*> function2V (graphsAll.size(), new TF1("myfunction2", myfunction2, 0.01, 100.0, 2));
	// Create graphs of Efficiency vs Mu from graphsAll and add to graphsAll2.
	// Run through all graphs in graphsAll2 and modify error bars if necessary.
	// Set graph properties of all graphs in graphsAll2 and draw to c1 pad 2.
	for (unsigned int i=0; i < graphsAll.size(); i++) {
		graphsAll2.push_back( (TGraphErrors*) (graphsAll[i])->Clone() );
		
		int n = graphsAll2[i]->GetN();

		double* xV = graphsAll2[i]->GetX();
		vector<double> xVector;
		for (int j=0; j < n; j++) {
			xVector.push_back(1 / xV[j]);
		}
		double* yV = graphsAll2[i]->GetY();
		double* exV = graphsAll2[i]->GetEX();
		double* eyV = graphsAll2[i]->GetEY();

		graphsAll2[i] = new TGraphErrors(n, &xVector[0], yV, exV, eyV);
		graphsAll2[i]->SetTitle("");//("Light collection efficiency vs Absorption Length");
		graphsAll2[i]->GetXaxis()->SetTitle("1 / Absorption length (cm)");
		graphsAll2[i]->GetYaxis()->SetTitle("Efficiency");
		graphsAll2[i]->GetXaxis()->SetLimits(0.005, 0.25);
		graphsAll2[i]->GetYaxis()->SetRangeUser(-0.02, 0.14);
		graphsAll2[i]->SetMarkerColor(colors[i]);
		graphsAll2[i]->SetMarkerSize(0.5);
		graphsAll2[i]->SetMarkerStyle(21);

		c1->cd(2);
		if(i==0) {
			graphsAll2[i]->Draw("AP");
		}
		else {
			graphsAll2[i]->Draw("P same");
		}
	}

	// Save c1 as image.
	{
		c1->SaveAs("Efficiency vs Lambda and Mu.png");
	}

	// Fit all graphs in graphsAll2 to myfunction2.
	for (unsigned int i=0; i < graphsAll2.size(); i++) {
		Double_t par[] = {1.0, 0.5};
		function2V[i]->SetParameters(par);
		function2V[i]->SetLineColor(colors[i]);
		function2V[i]->SetParNames("average track length", "constant", "offset");
		graphsAll2[i]->Fit(function2V[i], "", "", 0.005, 10);
	}

	vector<double> trackLengths, ietas;
	// Create vectors of average track lengths and ieta from graphsAll2.
	for (unsigned int i=0; i < function2V.size(); i++) {
		trackLengths.push_back( (function2V[i]->GetParameters())[0] );
		ietas.push_back( i+21 );
	}

	// Create graph of average track lengths vs ieta.
	// Draw to new canvas.
	if (doPlot) {
		TCanvas* trackLengthCanvas = new TCanvas("track lengths", "track lengths", 800, 600);

		TGraph* graphTrackLength = new TGraph( trackLengths.size(), &ietas[0], &trackLengths[0] );
		graphTrackLength->SetTitle("");
		graphTrackLength->GetXaxis()->SetTitle("ieta");
		graphTrackLength->GetYaxis()->SetTitle("Average track length (cm)");
		graphTrackLength->GetXaxis()->SetLabelSize(0.035);
		graphTrackLength->GetYaxis()->SetLabelSize(0.035);
	//	graphTrackLength->GetXaxis()->SetLimits(0.3,200);
	//	graphTrackLength->GetYaxis()->SetRangeUser(1e-6,3e-1);
		graphTrackLength->SetMarkerColor(colors[0]);
		graphTrackLength->SetMarkerSize(0.5);
		graphTrackLength->SetMarkerStyle(21);

		trackLengthCanvas->cd();
		graphTrackLength->Draw("AP");
	}

	// Read laser data file.
	vector<double> lumis;
	std::map<int, vector<double> > ratiosMap;
	readLaserDataFile("HELaserData.txt", lumis, ratiosMap);

	// Read dose and fluence map data files.
	std::map<int, double> doseMapLayer1;
	std::map<int, double> fluenceMapLayer1;
	readDoseFile("doseLayer1.txt", doseMapLayer1, fluenceMapLayer1);

	TLegend* legEfficiencyDose = new TLegend(0.7, 0.2, 0.9, 0.6);
	TLegend* legEfficiencyFluence = new TLegend(0.7, 0.2, 0.9, 0.6);
	//vector<TGraph*> graphsEfficiencyDoseLinear, graphsEfficiencyDoseLog;
	vector<double>  kParameterDoseVector, kParameterFluenceVector, ietaVector;
	
	// Run through all ietas. 
	for (int i=29; i>=21; i--)
	{
		stringstream name;
		name << "ieta " << (i);
		TCanvas *canvasDose, *canvasFluence;
		// Create canvases with 2 pads each. One for Dose and one for Fluence.
		if (doPlot) {
			canvasDose = new TCanvas(name.str().c_str(), name.str().c_str(), 1200, 600);
			canvasDose->SetLogx(0);
			canvasDose->SetLogy(0);
			canvasDose->Divide(2,0);
			canvasFluence = new TCanvas(name.str().c_str(), name.str().c_str(), 1200, 600);
			canvasFluence->SetLogx(0);
			canvasFluence->SetLogy(0);
			canvasFluence->Divide(2,0);
		}

		vector<double> xVector, xVector2, yVector, yVector2;
		// Calculate measured Lambda, Mu, Dose, Fluence and add to vectors	
		for (int j = minLumiIndex; j < maxLumiIndex + 1; j++)
		{
			double* fitParameters = functionV[i - 21]->GetParameters();
			double* fitParameters2 = function2V[i - 21]->GetParameters();
			double refEfficiency = myfunction( &refAttLength, fitParameters );
			double measuredEfficiency = ratiosMap[i][j] * refEfficiency;
			double measuredLambda = invmyfunction( &measuredEfficiency, fitParameters );
			double measuredMu = invmyfunction2( &measuredEfficiency, fitParameters2 );

			double refDose = doseMapLayer1[i];
			double measuredDose = refDose * (lumis[j] - lumioffset) / 100;

			double refFluence = fluenceMapLayer1[i];
			double measuredFluence = refFluence * (lumis[j] - lumioffset) / 100;

			xVector.push_back(measuredDose);
			xVector2.push_back(measuredFluence);
			yVector.push_back(measuredLambda);
			yVector2.push_back(measuredMu);
			
			//cout << "i, j = " << i << ", " << j << endl;
		}

		TGraph *graphDose, *graphFluence, *graphDose2, *graphFluence2;
		// Create graphs of Lambda vs Dose and Lambda vs Fluence.
		{
	   		graphDose = new TGraph( (maxLumiIndex - minLumiIndex + 1), &xVector[0], &yVector[0] );
			graphDose->GetXaxis()->SetTitle("Dose (kGy)");
			graphDose->SetTitle("");
			graphDose->GetYaxis()->SetTitle("Extracted absorption length (cm)");
			graphDose->GetXaxis()->SetLabelSize(0.035);
			graphDose->GetYaxis()->SetLabelSize(0.035);
	//		graphDose->GetXaxis()->SetLimits(0.3,200);
	//		graphDose->GetYaxis()->SetRangeUser(1e-6,3e-1);
			graphDose->SetMarkerColor(colors[0]);
			graphDose->SetMarkerSize(0.5);
			graphDose->SetMarkerStyle(21);

	   		graphFluence = new TGraph( (maxLumiIndex - minLumiIndex + 1), &xVector2[0], &yVector[0] );
			graphFluence->GetXaxis()->SetTitle("Neutral hadron fluence (10^12 cm^-1)");
			graphFluence->SetTitle("");
			graphFluence->GetYaxis()->SetTitle("Extracted absorption length (cm)");
			graphFluence->GetXaxis()->SetLabelSize(0.035);
			graphFluence->GetYaxis()->SetLabelSize(0.035);
	//		graphFluence->GetXaxis()->SetLimits(0.3,200);
	//		graphFluence->GetYaxis()->SetRangeUser(1e-6,3e-1);
			graphFluence->SetMarkerColor(colors[0]);
			graphFluence->SetMarkerSize(0.5);
			graphFluence->SetMarkerStyle(21);
		}

		// Fit graphDose and graphFluence to LambdaVsDose.
		{
			TF1* functionDose = new TF1("LambdaVsDose", LambdaVsDose, 0.0, 1e+5, 2);
			double parameters[] = {0.1, 0.1, 1.0};
			functionDose->SetParameters(parameters);
			functionDose->SetLineColor(colors[1]);
			functionDose->SetParNames("p0", "p1", "p2");
			graphDose->Fit(functionDose, "", "", 0.0, 1e+5);

			TF1* functionFluence = (TF1*) functionDose->Clone();
			graphFluence->Fit(functionFluence, "", "", 0.0, 1e+5);
		}

		legEfficiencyDose->AddEntry(graphDose, name.str().c_str(), "p");
		legEfficiencyFluence->AddEntry(graphFluence, name.str().c_str(), "p");

		// Draw graphDose and graphFluence to canvasDose and canvasFluence (pad 1).
		if (canvasDose && canvasFluence && doPlot) {
			canvasDose->cd(1);
			graphDose->Draw("AP");
			canvasFluence->cd(1);
			graphFluence->Draw("AP");
		}

		// Create graphs of Mu vs Dose and Mu vs Fluence.
		{
	   		graphDose2 = new TGraph( (maxLumiIndex - minLumiIndex + 1), &xVector[0], &yVector2[0] );
			graphDose2->GetXaxis()->SetTitle("Dose (kGy)");
			graphDose2->SetTitle("");
			graphDose2->GetYaxis()->SetTitle("1 / Extracted asorption length (cm^-1)");
			graphDose2->GetXaxis()->SetLabelSize(0.035);
			graphDose2->GetYaxis()->SetLabelSize(0.035);
		//	graphDose2->GetXaxis()->SetLimits(0.3,200);
		//	graphDose2->GetYaxis()->SetRangeUser(1e-6,3e-1);
			graphDose2->SetMarkerColor(colors[0]);
			graphDose2->SetMarkerSize(0.5);
			graphDose2->SetMarkerStyle(21);

	   		graphFluence2 = new TGraph( (maxLumiIndex - minLumiIndex + 1), &xVector2[0], &yVector2[0] );
			graphFluence2->GetXaxis()->SetTitle("Neutral hadron fluence (10^12 cm^-1)");
			graphFluence2->SetTitle("");
			graphFluence2->GetYaxis()->SetTitle("1 / Extracted asorption length (cm^-1)");
			graphFluence2->GetXaxis()->SetLabelSize(0.035);
			graphFluence2->GetYaxis()->SetLabelSize(0.035);
		//	graphFluence2->GetXaxis()->SetLimits(0.3,200);
		//	graphFluence2->GetYaxis()->SetRangeUser(1e-6,3e-1);
			graphFluence2->SetMarkerColor(colors[0]);
			graphFluence2->SetMarkerSize(0.5);
			graphFluence2->SetMarkerStyle(21);
		}

		// Fit graphDose2 and graphFluence2 to LambdaVsDose2.
		// Add k-parameter of Mu vs Dose and Mu vs Fluence to vectors.
		{
			TF1* functionDose2 = new TF1("LambdaVsDose", LambdaVsDose2, 0.0, 1e+5, 2);
			double parameters[] = {0.1, 0.1, 1.0};
			functionDose2->SetParameters(parameters);
			functionDose2->SetLineColor(colors[1]);
			functionDose2->SetParNames("p0", "p1", "p2");
			graphDose2->Fit(functionDose2, "", "", 0.0, 1e+5);

			TF1* functionFluence2 = (TF1*) functionDose2->Clone();
			graphFluence2->Fit(functionFluence2, "", "", 0.0, 1e+5);

			kParameterDoseVector.push_back( (functionDose2->GetParameters())[0] );
			kParameterFluenceVector.push_back( (functionFluence2->GetParameters())[0] );
			ietaVector.push_back(i);
		}

		// Draw graphDose2 and graphFluence2 to canvasDose and canvasFluence (pad 2).
		if (canvasDose && canvasFluence && doPlot) {
			canvasDose->cd(2);
			graphDose2->Draw("AP");
			canvasFluence->cd(2);
			graphFluence2->Draw("AP");
		}

	//legEfficiencyDose->SetFillColor(0);
	//legEfficiencyDose->SetBorderSize(0);
	//legEfficiencyDose->SetTextSize(0.03);
	//legEfficiencyDose->SetTextFont(42);
	//legEfficiencyDose->Draw();

		// Save Lambda (Mu) vs Dose and Lambda (Mu) vs Fluence plots.
		if (doSavePlots && doPlot) {
			stringstream filenameDose, filenameFluence;
			filenameDose << "plots/LambdaVSDose_ieta" << (i) << ".png";
			filenameFluence << "plots/LambdaVSFluence_ieta" << (i) << ".png";
			TImageDump *imgdumpDose = new TImageDump(filenameDose.str().c_str());
			canvasDose->Paint();
			imgdumpDose->Close();
			TImageDump *imgdumpFluence = new TImageDump(filenameFluence.str().c_str());
			canvasFluence->Paint();
			imgdumpFluence->Close();
		}
	}

	// Plot k-parameter vs ieta. k is slope of Mu vs Dose, or Mu vs Fluence.
	{
		string titleDose, ytitleDose;
		titleDose = "k-parameter [Dose]";
		ytitleDose = "k-parameter [cm^{-1} / kGy]";
		string titleFluence, ytitleFluence;
		titleFluence = "k-parameter [Fluence]";
		ytitleFluence = "k-parameter [cm^{-1} / cm^{-2}]";

		TCanvas* canvasDose = new TCanvas(titleDose.c_str(), titleDose.c_str(), 800, 600);

		TGraph* graphDose = new TGraph( kParameterDoseVector.size(), &ietaVector[0], &kParameterDoseVector[0] );
		graphDose->SetTitle("");
		graphDose->GetXaxis()->SetTitle("ieta");
		graphDose->GetYaxis()->SetTitle(ytitleDose.c_str());
		graphDose->GetXaxis()->SetLabelSize(0.035);
		graphDose->GetYaxis()->SetLabelSize(0.035);
	//	graphDose->GetXaxis()->SetLimits(0.3,200);
	//	graphDose->GetYaxis()->SetRangeUser(1e-6,3e-1);
		graphDose->SetMarkerColor(colors[0]);
		graphDose->SetMarkerSize(1.0);
		graphDose->SetMarkerStyle(21);

		canvasDose->cd();
		graphDose->Draw("AP");
		canvasDose->SaveAs("k-parameter (Dose).png");

		TCanvas* canvasFluence = new TCanvas(titleFluence.c_str(), titleFluence.c_str(), 800, 600);

		TGraph* graphFluence = new TGraph( kParameterFluenceVector.size(), &ietaVector[0], &kParameterFluenceVector[0] );
		graphFluence->SetTitle("");
		graphFluence->GetXaxis()->SetTitle("ieta");
		graphFluence->GetYaxis()->SetTitle(ytitleFluence.c_str());
		graphFluence->GetXaxis()->SetLabelSize(0.035);
		graphFluence->GetYaxis()->SetLabelSize(0.035);
	//	graphFluence->GetXaxis()->SetLimits(0.3,200);
	//	graphFluence->GetYaxis()->SetRangeUser(1e-6,3e-1);
		graphFluence->SetMarkerColor(colors[0]);
		graphFluence->SetMarkerSize(1.0);
		graphFluence->SetMarkerStyle(21);

		canvasFluence->cd();
		graphFluence->Draw("AP");
		canvasFluence->SaveAs("k-parameter (Fluence).png");
	}
}


