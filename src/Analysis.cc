#include "Analysis.hh"
#include "AnalysisMessenger.hh"
#include "G4UnitsTable.hh"
#include "G4SDManager.hh"
#include "G4PhysicalConstants.hh"
#include "LYSimDetectorConstruction.hh"
#include "G4ios.hh"
#include <iostream>
#include <fstream>
using namespace std;

//ROOT Stuff
//#include "TProfile.h"
#include "TFile.h"
#include "TH1D.h"

Analysis* Analysis::singleton = 0;

//Constructor
Analysis::Analysis()
{
	fMessenger = new AnalysisMessenger(this);
}

Analysis::~Analysis()
{
	if(fMessenger) delete fMessenger;
}


void Analysis::PrepareNewEvent(const G4Event* /*anEvent*/)
{
}

void Analysis::EndOfEvent(const G4Event* anEvent)
{
	G4String hitCollName = "PMTHitsCollection";
	G4SDManager* SDman = G4SDManager::GetSDMpointer();	
	static G4int hitCollID = -1;
	if ( hitCollID < 0 )
	hitCollID = SDman->GetCollectionID( hitCollName );
	
	G4HCofThisEvent* hitsCollections = 0;
	hitsCollections = anEvent->GetHCofThisEvent();
	
	LYSimPMTHitsCollection* hits = 0;
	
	if ( hitsCollections )
	{
		hits = static_cast<LYSimPMTHitsCollection*> ( hitsCollections->GetHC(hitCollID) );
	}
	else
	{
		G4cerr << "hitsCollection not found" << G4endl;
		return;
	}
	LYSimPMTHit* hit = NULL;
	G4double EventEnergy = 0;
	G4int EventPhotonCount = 0;
	PhotonCount++;
	if(hits->entries()==1) //if there is an entry in hits collection
	{
		hit = (*hits)[0];
		EventEnergy = hit->GetEnergy();
		EventPhotonCount = hit->GetPhotonCount();
		HitCount++;
		G4double trackLength = hit->GetTrackLength();
		TrackLengthHisto->Fill(trackLength/cm);
	}
	//G4cout	<< "Energy deposited in PMT this event is " << G4BestUnit(EventEnergy,"Energy") << G4endl
	//		<< "# of photon hits in PMT this event is " << EventPhotonCount << G4endl;
	//EnergyHist->Fill(EventEnergy/eV);
	//PhotonHitsHist->Fill(EventPhotonCount);
	
	//// Add TrackLength to histogram.
	//int trajectorySize = anEvent->GetTrajectoryContainer()->size();	
	//cout << "#############" << endl << "Number of trajectories is " << trajectorySize << endl;

}


void Analysis::PrepareNewRun(const G4Run* /*aRun*/ )
{
	G4cout << "Analysis::PrepareNewRun() begin" << G4endl;
	//Reset variables relative to the run
	PhotonCount = 0;
	HitCount = 0; //
	RootFile = new TFile("test.root", "RECREATE");
	TrackLengthHisto = new TH1D("test", "test histogram of tracklengths", 50, 0., 5000.);
}


void Analysis::EndOfRun(const G4Run* aRun)
{
	G4cout << "Analysis::EndOfRun() begin" << G4endl;
	
	//Write to file
	//~ TFile* rootfile = TFile::Open(fROOTFileName.c_str(),"recreate");
	//~ EnergyHist->Write();
	//~ PhotonHitsHist->Write();
	//~ rootfile->Close();
	//~ 
	std::ofstream outputfile;
	outputfile.open(fOutputFileName.c_str(), ofstream::out | ofstream::app);
	G4cout	<< "Efficiency in this run is " << (G4double)HitCount/(G4double)PhotonCount << G4endl;
	if (outputfile.is_open())
	{
		outputfile << "#Mu [cm^-1] \t Efficiency" << G4endl;
		outputfile << inducedMuTile << "\t"<< (G4double)HitCount/(G4double)PhotonCount << G4endl;
	}
	else
	{
		G4cout << "Output file not open" << G4endl;
		G4cout << "#Mu [cm^-1] \t Efficiency" << G4endl;
		G4cout << inducedMuTile << "\t"<< (G4double)HitCount/(G4double)PhotonCount << G4endl;
	}
	outputfile.close();

	if(TrackLengthHisto) {TrackLengthHisto->Write();}
	RootFile->Close();
	G4cout << "Analysis::EndOfRun() end" << G4endl;
}
