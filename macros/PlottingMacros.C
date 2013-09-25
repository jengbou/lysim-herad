void GetTrees(TTree* &mc, TTree* &data);

void Plot()
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

	TTree *mcRatios, *dataRatios;
	GetTrees(mcRatios, dataRatios);

	mcRatios->Draw("ieta:lumi:ratio:lumiIndex", "ieta==29", "goff");
	int nMC = mcRatios->GetSelectedRows();
	double* temp; 

	temp = mcRatios->GetVal(0);
	double* ietasMC;
	ietasMC = new double[nMC];
	for(int i=0; i<nMC; i++) {ietasMC[i] = temp[i];}

	temp = mcRatios->GetVal(1);
	double* lumisMC;
	lumisMC = new double[nMC];
	for(int i=0; i<nMC; i++) {lumisMC[i] = temp[i];}

	temp = mcRatios->GetVal(2);
	double* ratiosMC;
	ratiosMC = new double[nMC];
	for(int i=0; i<nMC; i++) {ratiosMC[i] = temp[i];}

	dataRatios->Draw("ieta:lumi:ratio:lumiIndex", "ieta==29", "goff");
	int nData = dataRatios->GetSelectedRows();
	double* ietasData = dataRatios->GetVal(1);
	double* lumisData = dataRatios->GetVal(2);
	double* ratiosData = dataRatios->GetVal(3);

	temp = dataRatios->GetVal(0);
	double* ietasData;
	ietasData = new double[nData];
	for(int i=0; i<nData; i++) {ietasData[i] = temp[i];}

	temp = dataRatios->GetVal(1);
	double* lumisData;
	lumisData = new double[nData];
	for(int i=0; i<nData; i++) {lumisData[i] = temp[i];}

	temp = dataRatios->GetVal(2);
	double* ratiosData;
	ratiosData = new double[nData];
	for(int i=0; i<nData; i++) {ratiosData[i] = temp[i];}

	TCanvas* canvas = new TCanvas("canvas","canvas title",800,600);
	//canvas->Divide(2,1);

	TGraph* graphData;
	graphData = new TGraph(nData, lumisData, ratiosData);
	graphData->SetTitle("");
	graphData->GetXaxis()->SetTitle("lumi [fb^{-1}]");
	graphData->GetYaxis()->SetTitle("damaged/undamaged ratio");
	graphData->GetYaxis()->SetRangeUser(0.5, 1.2);
	graphData->SetMarkerColor(colors[0]);
	canvas->cd(1);
	graphData->Draw("AP");

	TGraph* graphMC;
	graphMC = new TGraph(nMC, lumisMC, ratiosMC);
	graphMC->SetMarkerColor(colors[5]);
	canvas->cd(1);
	graphMC->Draw("P same");

	TLegend* leg = new TLegend(0.2, 0.2, 0.4, 0.3);
	leg->AddEntry(graphData, "Data", "p");
	leg->AddEntry(graphMC, "Simulation", "p");

	leg->SetFillColor(0);
	leg->SetBorderSize(0);
//	leg->SetTextSize(0.03);
	leg->SetTextFont(42);
	leg->Draw();

//	delete lumisData;
//	cout<<ietas<<endl;
}

void GetTrees(TTree* &mc, TTree* &data)
{
	TFile* dataFile = new TFile("dataRatios.root");
	data = (TTree*)dataFile->Get("dataRatios");

	TFile* mcFile = new TFile("mcRatios2.root");
	mc =  (TTree*)mcFile->Get("mcRatios");
}
