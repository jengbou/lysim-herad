#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <fstream>
#include <TROOT.h>

using namespace std;

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

TGraph* readAnalysisFile(string inputFile)
{
	vector<double> absLengthV;
	vector<double> efficiencyV;
	TGraph* gr;
	TCanvas* c1;
	
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
				processLine(line, absLength, efficiency);
				absLengthV.push_back(absLength);
				efficiencyV.push_back(efficiency);
			}
		}
		Int_t n = absLengthV.size();
		gr = new TGraph(n, &absLengthV[0], &efficiencyV[0]);
		gr->SetTitle("");//("Light collection efficiency vs Absorption Length");
		gr->GetXaxis()->SetTitle("Absorption length (cm)");
		gr->GetYaxis()->SetTitle("Efficiency");
		gr->SetMarkerSize(1.0);
		gr->SetMarkerStyle(21);
		cout << gr << endl;
		gr->Draw("ALP");
	}
	else {cout<<"Error: Input file not found." <<endl;}
	//cout << gr << endl;
	return gr;
}

void longrun3()
{
	// Vector of colors.
	vector<Color_t> colors(36, kBlack);
	colors[0] = kBlue;
	colors[1] = kOrange;
	colors[2] = kGreen;
	colors[3] = kRed;
	colors[4] = kCyan;
	colors[5] = kMagenta;
	for (int i = 0; i < 6; i++ ) {
		colors[i+6] = colors[i] + 2;
	}
	
	
	c1 = new TCanvas("c1","canvas1",700,500);
	c1->SetLogx();
	c1->SetLogy();
	vector<TGraph*> graphs;
	TLegend* leg = new TLegend(0.7, 0.2, 0.9, 0.5);

	graphs.push_back(readAnalysisFile("Analysis-longrun3_ieta22.txt"));
	graphs.push_back(readAnalysisFile("Analysis-longrun3_ieta23.txt"));
	graphs.push_back(readAnalysisFile("Analysis-longrun3_ieta24.txt"));
	graphs.push_back(readAnalysisFile("Analysis-longrun3_ieta25.txt"));
	graphs.push_back(readAnalysisFile("Analysis-longrun3_ieta26.txt"));
	graphs.push_back(readAnalysisFile("Analysis-longrun3_ieta27.txt"));
	graphs.push_back(readAnalysisFile("Analysis-longrun3_ieta28.txt"));
	graphs.push_back(readAnalysisFile("Analysis-longrun3_ieta29.txt"));
	for (int i = 0; i < 8; i++) {
		cout<< graphs[i] << endl;
 		graphs[i]->GetYaxis()->SetRangeUser(1e-5,5e-1);
 		graphs[i]->SetMarkerColor(colors[i]);
		if(i==0) {
			graphs[i]->Draw("ALP");
		}
		else {
			graphs[i]->Draw("LP same");
		}
 		
 		stringstream legendName;
 		legendName << "ieta " << (i+22);
 		cout<< "ieta " << (i+22) << endl;
 		leg->AddEntry(graphs[i], legendName.str().c_str(), "p");
	}
	leg->Draw();
// 	TGraph* grIeta22 = readAnalysisFile("Analysis-longrun3_ieta22.txt");
// 	TGraph* grIeta23 = readAnalysisFile("Analysis-longrun3_ieta23.txt");
// 	TGraph* grIeta24 = readAnalysisFile("Analysis-longrun3_ieta24.txt");
// 	TGraph* grIeta25 = readAnalysisFile("Analysis-longrun3_ieta25.txt");
// 	TGraph* grIeta26 = readAnalysisFile("Analysis-longrun3_ieta26.txt");
// 	TGraph* grIeta27 = readAnalysisFile("Analysis-longrun3_ieta27.txt");
// 	TGraph* grIeta28 = readAnalysisFile("Analysis-longrun3_ieta28.txt");
// 	TGraph* grIeta29 = readAnalysisFile("Analysis-longrun3_ieta29.txt");
// 	grIeta22->GetYaxis()->SetRangeUser(1e-5,5e-1);
// 	grIeta22->SetMarkerColor(colors[0]);
// 	grIeta23->SetMarkerColor(colors[1]);
// 	grIeta24->SetMarkerColor(colors[2]);
// 	grIeta25->SetMarkerColor(colors[3]);
// 	grIeta26->SetMarkerColor(colors[4]);
// 	grIeta27->SetMarkerColor(colors[5]);
// 	grIeta28->SetMarkerColor(colors[6]);
// 	grIeta29->SetMarkerColor(colors[7]);
// 	grIeta22->Draw("ALP");
// 	grIeta23->Draw("LP same");
// 	grIeta24->Draw("LP same");
// 	grIeta25->Draw("LP same");
// 	grIeta26->Draw("LP same");
// 	grIeta27->Draw("LP same");
// 	grIeta28->Draw("LP same");
// 	grIeta29->Draw("LP same");

// 	leg->AddEntry(grIeta22, "ieta 22", "p");	
// 	leg->AddEntry(grIeta23, "ieta 23", "p");
// 	leg->AddEntry(grIeta24, "ieta 24", "p");
// 	leg->AddEntry(grIeta25, "ieta 25", "p");
// 	leg->AddEntry(grIeta26, "ieta 26", "p");
// 	leg->AddEntry(grIeta27, "ieta 27", "p");
// 	leg->AddEntry(grIeta28, "ieta 28", "p");
// 	leg->AddEntry(grIeta29, "ieta 29", "p");
// 	leg->Draw();
}
