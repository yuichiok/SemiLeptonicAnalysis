#include "ttbarTreeWriter.h"
#include <iostream>
#include <math.h>
#include <sstream>

#include <LCIOSTLTypes.h>
#include <EVENT/LCCollection.h>
#include <EVENT/ReconstructedParticle.h>
#include <UTIL/PIDHandler.h>
#include <EVENT/MCParticle.h>
#include <EVENT/Vertex.h>
#include <HelixClass.h>
#include <UTIL/LCRelationNavigator.h>
// ----- include for verbosity dependend logging ---------
#include "marlin/VerbosityLevels.h"

// GEAR include files
#include <marlin/Global.h>
#include <gear/GEAR.h>
#include <gear/TPCParameters.h>
#include <gear/CalorimeterParameters.h>
#include <gear/BField.h>

using namespace lcio;
using namespace marlin;

ttbarTreeWriter attbarTreeWriter;

ttbarTreeWriter::ttbarTreeWriter() : Processor("ttbarTreeWriter")
{
  
  // processor description

  _description = "Analyse ttbar to semileptonic mode at Reconstruction and MC level";
  
  
  // register steering parameters: name, description, class-variable, default value

  registerProcessorParameter( "ROOTFileName",
                              "Output ROOT File Name",
                              _hfilename,
                              std::string("ttbarTreeWriter.root") );

  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			   "PandoraCollection" ,
			   "Output collection given by Pandora" ,
			   _pfosCollection ,
			   std::string("PandoraPFOs") );

  registerInputCollection( LCIO::MCPARTICLE,
			   "MCCollection" ,
			   "MC collection (MCParticle for Mokka, MCParticlesSkimmed for DST)" ,
			   _mcCollection ,
			   std::string("MCParticle") );
  
  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
                           "MyPFOsCollection" ,
                           "Name of my collection"  ,
                           _myPfosCollection ,
                           std::string("PhD_PFOs") );

  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
                           "4Jets" ,
                           "Name of the 4 Jets collection"  ,
                           _myJetsCollection ,
                           std::string("PhD_Durham4Jets") );

  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
                           "Lepton" ,
                           "Name of the lepton collection"  ,
                           _myLeptonCollection ,
                           std::string("PhD_SelectedLepton") );

}


void ttbarTreeWriter::init()
{ 
  
  _nRun = 0;
  _nEvt = 0;

  ROOTfile = new TFile(_hfilename.c_str(), "RECREATE", _hfilename.c_str());
  
  _tree = new TTree("tree", "tree");
  
  _tree->Branch("runNumber", &_nRun, "runNumber/I");
  _tree->Branch("eventNumber", &_nEvt, "eventNumber/I");
  _tree->Branch("crossSection", &_xsec, "crossSection/F");

  _Process = new TString();
  _tree->Branch("Process","TString",&_Process,16000,0);
  _tree->Branch("isTTbarEvent", &_isttbarevent, "isTTbarEvent/O");
  _tree->Branch("isTTbarSL", &_isttsl, "isTTbarSL/O");
  _tree->Branch("isZWWSL", &_iszwwsl, "isZWWSL/O");
  _tree->Branch("isBackground", &_isbkg, "isBackground/O");


  //  Reconstructed variables 

  for(int i = 0 ; i < 4 ; i++)
    {
      _jets[i] = new TLorentzVector();
      std::stringstream name;
      name << "jet" << i;
      _tree->Branch(name.str().c_str(),"TLorentzVector",&_jets[i],16000,0);
      _tree->Branch((name.str()+"_btag").c_str(), &_jetsBtag[i], (name.str()+"_btag/F").c_str());
      _tree->Branch((name.str()+"_ctag").c_str(), &_jetsCtag[i], (name.str()+"_ctag/F").c_str());
      _tree->Branch((name.str()+"_VertexCharge").c_str(), &_jetsVertexCharge[i], (name.str()+"_VertexCharge/I").c_str());
      _tree->Branch((name.str()+"_VertexChargeTer").c_str(), &_jetsVertexChargeTer[i], (name.str()+"_VertexChargeTer/I").c_str());
      _tree->Branch((name.str()+"_TrueFlavour").c_str(), &_jetsTrueFlavour[i], (name.str()+"_TrueFlavour/I").c_str());
      _tree->Branch((name.str()+"_QuarkOrigin").c_str(), &_jetsQuarkOrigin[i], (name.str()+"_QuarkOrigin/F").c_str());
      _tree->Branch((name.str()+"_TrueCharge").c_str(), &_jetsTrueCharge[i], (name.str()+"_TrueCharge/I").c_str());
      _tree->Branch((name.str()+"_TrueBmesonCharge").c_str(), &_jetsTrueBmesonCharge[i], (name.str()+"_TrueBmesonCharge/I").c_str());
      _tree->Branch((name.str()+"_DecayLength").c_str(), &_jetsDecayLength[i], (name.str()+"_DecayLength/F").c_str());
      _tree->Branch((name.str()+"_DecayLengthSig").c_str(), &_jetsDecayLengthSig[i], (name.str()+"_DecayLengthSig/F").c_str());
      _tree->Branch((name.str()+"_NumofTracks").c_str(), &_jetsNumofTracks[i], (name.str()+"_NumofTracks/I").c_str());
      _tree->Branch((name.str()+"_NumVertices").c_str(), &_jetsNumVertices[i], (name.str()+"_NumVertices/I").c_str());
      _tree->Branch((name.str()+"_VtxMom").c_str(), &_jetsVtxMom[i], (name.str()+"_VtxMom/F").c_str());
      _tree->Branch((name.str()+"_VtxPtmass").c_str(), &_jetsVtxPtmass[i], (name.str()+"_VtxPtmass/F").c_str());
      _tree->Branch((name.str()+"_NotIsoMuFound").c_str(), &_NotIsoMuonFound[i], (name.str()+"_NotIsoMuFound/O").c_str());
      _tree->Branch((name.str()+"_NotIsoMuEnergy").c_str(), &_NotIsoMuons_energy[i], (name.str()+"_NotIsoMuEnergy/F").c_str());
      _tree->Branch((name.str()+"_NotIsoMuz4").c_str(), &_NotIsoMuons_z4[i], (name.str()+"_NotIsoMuz4/F").c_str());
      _tree->Branch((name.str()+"_NotIsoMupT").c_str(), &_NotIsoMuons_pT[i], (name.str()+"_NotIsoMupT/F").c_str());
      _tree->Branch((name.str()+"_NotIsoMuxT").c_str(), &_NotIsoMuons_xT[i], (name.str()+"_NotIsoMuxT/F").c_str());
      _tree->Branch((name.str()+"_NotIsoMucharge").c_str(), &_NotIsoMuons_charge[i], (name.str()+"_NotIsoMucharge/I").c_str());
    }

  _tree->Branch("YPlus", &_yplus, "YPlus/F");
  _tree->Branch("YMinus", &_yminus, "YMinus/F");
  _tree->Branch("isGoodBjets", &_isGoodBs, "isGoodBjets/O");
  _tree->Branch("isGoodCombination", &_isGoodComb, "isGoodCombination/O");

  _myLepton = new TLorentzVector();
  _tree->Branch("lepton","TLorentzVector",&_myLepton,16000,0);
  _tree->Branch("leptonPDG", &_lepPDG, "leptonPDG/I");
  _tree->Branch("leptonCharge", &_lepCharge, "leptonCharge/F");
  _tree->Branch("NbCandidates", &_nbOfLep, "NbCandidates/I");
  _tree->Branch("isGoodLepton", &_isGoodLep, "isGoodLepton/O");
  _tree->Branch("isLepFromTau", &_isLepFromTau, "isLepFromTau/O");
  _totalWOlepton = new TLorentzVector();
  _tree->Branch("allWOlepton","TLorentzVector",&_totalWOlepton,16000,0);
  _Whad = new TLorentzVector();
  _tree->Branch("Whad","TLorentzVector",&_Whad,16000,0);
  _tree->Branch("SelectedBjet", &_selectedBjet, "SelectedBjet/I");
  _top_had = new TLorentzVector();
  _tree->Branch("top_had","TLorentzVector",&_top_had,16000,0);
  _tree->Branch("chi2DBD", &_chi2DBD, "chi2DBD/F");
  _tree->Branch("chi2new", &_chi2new, "chi2new/F");
  _tree->Branch("EventCharge", &_eventcharge, "EventCharge/I");
  _tree->Branch("isLeptonVtxChargeOK", &_isLeptonVtxChargeOK, "isLeptonVtxChargeOK/O");
  
  _tree->Branch("numberOfPFOs", &_numberOfPFOs, "numberOfPFOs/I");
  _tree->Branch("pfosFoxWolfram", _pfosFoxWolfram, "pfosFoxWolfram[7]/F");
  _tree->Branch("pfosSphericity", &_pfosSphericity, "pfosSphericity/F");
  _tree->Branch("pfosAplanarity", &_pfosAplanarity, "pfosAplanarity/F");
  _tree->Branch("pfosOblateness", &_pfosOblateness, "pfosOblateness/F");
  _tree->Branch("pfosThrust", &_pfosThrust, "pfosThrust/F");


  // Monte Carlo variables 

  _MCtotal = new TLorentzVector();
  _tree->Branch("MCtotal","TLorentzVector",&_MCtotal,16000,0);
  
  _MCphoton[0] = new TLorentzVector();
  _tree->Branch("MCphoton1","TLorentzVector",&_MCphoton[0],16000,0);
  _MCphoton[1] = new TLorentzVector();
  _tree->Branch("MCphoton2","TLorentzVector",&_MCphoton[1],16000,0);
  
  
  for(int i = 0 ; i < 6 ; i++)
    {
      _MCf[i] = new TLorentzVector();
      std::stringstream name;
      name << "MCfermion" << i;
      _tree->Branch(name.str().c_str(),"TLorentzVector",&_MCf[i],16000,0);
      name << "_PDG";
      _tree->Branch(name.str().c_str(), &_MCfpdg[i], (name.str()+"/I").c_str());
    }
  
  _MCb = new TLorentzVector();
  _tree->Branch("MCb","TLorentzVector",&_MCb,16000,0);

  _MCbbar = new TLorentzVector();
  _tree->Branch("MCbbar","TLorentzVector",&_MCbbar,16000,0);
  
  _MClepton = new TLorentzVector();
  _tree->Branch("MClepton","TLorentzVector",&_MClepton,16000,0);
  _tree->Branch("MClepton_PDG", &_MClepton_PDG, "MClepton_PDG/I");
  
  _MCneutrino = new TLorentzVector();
  _tree->Branch("MCneutrino","TLorentzVector",&_MCneutrino,16000,0);
  _tree->Branch("MCneutrino_PDG", &_MCneutrino_PDG, "MCneutrino_PDG/I");

  _MCq[0] = new TLorentzVector();
  _tree->Branch("MCq1","TLorentzVector",&_MCq[0],16000,0);
  _tree->Branch("MCq1_PDG", &_MCqpdg[0], "MCq1_PDG/I");  

  _MCq[1] = new TLorentzVector();
  _tree->Branch("MCq2","TLorentzVector",&_MCq[1],16000,0);
  _tree->Branch("MCq2_PDG", &_MCqpdg[1], "MCq2_PDG/I");  
  
  _MCW_jet = new TLorentzVector();
  _tree->Branch("MCW_jet","TLorentzVector",&_MCW_jet,16000,0);

  _MCW_lepton = new TLorentzVector();
  _tree->Branch("MCW_lepton","TLorentzVector",&_MCW_lepton,16000,0);

  _MCtop_jet = new TLorentzVector();
  _tree->Branch("MCtop_jet","TLorentzVector",&_MCtop_jet,16000,0);
  _tree->Branch("MCtop_jet_sign", &_MCtop_jet_Q, "MCtop_jet_sign/I");

  _MCtop_lepton = new TLorentzVector();
  _tree->Branch("MCtop_lepton","TLorentzVector",&_MCtop_lepton,16000,0);
  _tree->Branch("MCtop_lepton_sign", &_MCtop_lepton_Q, "MCtop_lepton_sign/I");
 
  _MCbbbar = new TLorentzVector();
  _tree->Branch("MCbbbar","TLorentzVector",&_MCbbbar,16000,0);
  
  
}


void ttbarTreeWriter::processRunHeader(LCRunHeader* run)
{ 
  _nRun++;
} 


void ttbarTreeWriter::processEvent(LCEvent * evt) 
{ 
  
  _nEvt ++;
  if(_nEvt%50 == 0) std::cout << "****************   Event " << _nEvt << "   **********************" << std::endl;
  ROOTfile->cd();
  
  try {

    // Get Process name and cross section
    *_Process = evt->getParameters().getStringVal("Process");
    _xsec = evt->getParameters().getFloatVal("CrossSection_fb");
    
    _isttsl = false;
    _iszwwsl = false;
    _isbkg = false;
    bool havetop = false;
    bool havetopbar = false;
    _isttbarevent = false;
 
    // Get information on the MCParticles
    LCCollection* mccol = evt->getCollection(_mcCollection);
	

    // ********** MONTE CARLO STUFF **********

    // Initialise
    _MCtotal->SetPxPyPzE(0.,0.,0.,0.);
    _MCtop_jet->SetPxPyPzE(0.,0.,0.,0.);
    _MCtop_jet_Q = 0;
    _MCtop_lepton->SetPxPyPzE(0.,0.,0.,0.);
    _MCtop_lepton_Q = 0;
    _MCphoton[0]->SetPxPyPzE(0.,0.,0.,0.);
    _MCphoton[1]->SetPxPyPzE(0.,0.,0.,0.);
    _MCq[0]->SetPxPyPzE(0.,0.,0.,0.) ; _MCqpdg[0] = 0;
    _MCq[1]->SetPxPyPzE(0.,0.,0.,0.) ; _MCqpdg[1] = 0;
    for(int i = 0 ; i < 6 ; i++)
      { 
	_MCf[i]->SetPxPyPzE(0.,0.,0.,0.);
	_MCfpdg[i] = 0; 
      }
    _MCW_jet->SetPxPyPzE(0.,0.,0.,0.);
    _MCW_lepton->SetPxPyPzE(0.,0.,0.,0.);
    _MCb->SetPxPyPzE(0.,0.,0.,0.);
    _MCbbar->SetPxPyPzE(0.,0.,0.,0.);
    _MClepton->SetPxPyPzE(0.,0.,0.,0.);
    _MClepton_PDG = 0;
    _MCneutrino->SetPxPyPzE(0.,0.,0.,0.);
    _MCneutrino_PDG = 0;
	 

    if(mccol->getNumberOfElements() > 11)
      {
	for(int i = 0; i <  mccol->getNumberOfElements(); i++)
	  {
	    MCParticle *mcpart = dynamic_cast<MCParticle*>(mccol->getElementAt(i));
	    TLorentzVector mcVec(TVector3(mcpart->getMomentum()),mcpart->getEnergy());

	    if(mcpart->getGeneratorStatus() != 2) continue;
	
	    if(mcpart->getPDG() == 6) havetop = true;
	    if(mcpart->getPDG() == -6) havetopbar = true;	
	  }
	if(havetop == true && havetopbar == true) _isttbarevent = true;

       	// 2 photons + the 6 leptons
	for(int i = 4; i < 12; i++)
	  {
	    MCParticle *mcpart = dynamic_cast<MCParticle*>(mccol->getElementAt(i));
	    TLorentzVector mcVec(TVector3(mcpart->getMomentum()),mcpart->getEnergy());

	    if(mcpart->getGeneratorStatus() != 2) continue;

	    *_MCtotal += mcVec;
	
	    if(i < 6) *_MCphoton[i] = mcVec;	
	    else 
	      {
		*_MCf[i-6] = mcVec;
		_MCfpdg[i-6] = mcpart->getPDG();
	      }
	  }
	    
	// get b quarks, lepton and neutrino information
	for(int i = 0 ; i < 6 ; i++)
	  {
	    if(_MCfpdg[i] == 5) *_MCb = *_MCf[i]; // b quark has charge -1/3

	    else if(_MCfpdg[i] == -5) *_MCbbar = *_MCf[i]; // bbar quark has charge +1/3
	
	    else if(abs(_MCfpdg[i]) == 11 || abs(_MCfpdg[i]) == 13 || abs(_MCfpdg[i]) == 15 ) // lepton
	      {
		*_MClepton = *_MCf[i];
		_MClepton_PDG = _MCfpdg[i];
	      }
	
	    else if(abs(_MCfpdg[i]) == 12 || abs(_MCfpdg[i]) == 14 || abs(_MCfpdg[i]) == 16 ) // neutrino
	      {
		*_MCneutrino = *_MCf[i];
		_MCneutrino_PDG = _MCfpdg[i];
	      }

	    else if(abs(_MCfpdg[i]) <= 4 && abs(_MCfpdg[i])%2 == 0) // first non-b quark
	      {
		*_MCq[0] = *_MCf[i];
		_MCqpdg[0] = _MCfpdg[i];
	      }

	    else if(abs(_MCfpdg[i]) <= 4 && abs(_MCfpdg[i])%2 == 1) // second non-b quark
	      {
		*_MCq[1] = *_MCf[i];
		_MCqpdg[1] = _MCfpdg[i];
	      }
	  }

	// Reconstruct the b-bbar system
	*_MCbbbar = *_MCb + *_MCbbar;
    
	// reconstruct the two W's
	*_MCW_jet = *_MCq[0] + *_MCq[1];
	*_MCW_lepton = *_MClepton + *_MCneutrino;

	// reconstruct tops
	if(_MCqpdg[0] == 0 &&_MCqpdg[1] == 0) // i.e. two lepton 
	  {
	    _isbkg = true;
	  }
	else if(_MClepton_PDG > 0) // i.e. lepton charge = -1
	  {
	    *_MCtop_lepton = *_MCbbar + *_MCW_lepton;
	    _MCtop_lepton_Q = -1;
	    *_MCtop_jet = *_MCb + *_MCW_jet;
	    _MCtop_jet_Q = +1;
	  }

	else if(_MClepton_PDG < 0) // i.e. lepton charge = +1
	  {
	    *_MCtop_lepton = *_MCb + *_MCW_lepton;
	    _MCtop_lepton_Q = +1;
	    *_MCtop_jet = *_MCbbar + *_MCW_jet;
	    _MCtop_jet_Q = -1;
	  }
	else if(_MClepton_PDG == 0) // i.e. no lepton
	  {
	    _isbkg = true;
	  }
	   
	if(!_isbkg)
	  { 
	    if((TMath::Abs(_MCtop_lepton->M()-174)<5*1.51 || TMath::Abs(_MCtop_jet->M()-174)<5*1.51 || TMath::Sqrt(pow(_MCtop_jet->M()-174,2)+pow(_MCtop_lepton->M()-174,2))<(15*1.51)) && (_MCtop_jet->M()>100 && _MCtop_lepton->M()>100))
	      {
		_isttsl = true; 
	      }
	    else _iszwwsl = true;
	  }
      }
      
    // ********** END OF MONTE CARLO STUFF **********
 

    // Write PFOs information
    _numberOfPFOs = 0;
    LCCollection* col = evt->getCollection(_pfosCollection);
    _numberOfPFOs = col->getNumberOfElements();
    for(int i = 0 ; i < 7 ; i++)
      {
	_pfosFoxWolfram[i] = 0.;
	std::stringstream name;
	name << "FoxWolfram_moment(" << i << ")";
	_pfosFoxWolfram[i] = col->getParameters().getFloatVal(name.str());
      }
    _pfosSphericity = col->getParameters().getFloatVal("sphericity");
    _pfosAplanarity = col->getParameters().getFloatVal("aplanarity");
    _pfosOblateness = col->getParameters().getFloatVal("Oblateness");
    _pfosThrust = col->getParameters().getFloatVal("principleThrustValue");

    // Write myPFOs information
    col = evt->getCollection(_myPfosCollection);
    // Calculate total 4-momentum (WO lepton)
    _totalWOlepton->SetPxPyPzE(0.,0.,0.,0.);
    for(int i = 0 ; i < col->getNumberOfElements() ; i++)
      {
	ReconstructedParticle *a_part = dynamic_cast<ReconstructedParticle*>(col->getElementAt(i));
	TLorentzVector Partvec(TVector3(a_part->getMomentum()),a_part->getEnergy());
	*_totalWOlepton += Partvec;
      }

    // Initialise and write : LEPTON
    _myLepton->SetPxPyPzE(0.,0.,0.,0.);
    LCCollection* leptonCol = evt->getCollection(_myLeptonCollection);
    _nbOfLep = 0;
    _nbOfLep = leptonCol->getParameters().getIntVal("NumberOfLeptonsCandidates");
    _lepPDG = 0 ;
    _lepPDG = leptonCol->getParameters().getIntVal("SelectedLeptonPDG");
    _lepCharge = 0.;
    _isGoodLep = false;
    _isLepFromTau = false;
    if(leptonCol->getParameters().getIntVal("isGoodLeptonFound") == 1) _isGoodLep = true;
    if(leptonCol->getParameters().getIntVal("isLeptonFromTau") == 1) _isLepFromTau = true;
    for(int i = 0 ; i < (int) leptonCol->getNumberOfElements() ; i++)
      {
	ReconstructedParticle *mylepton = dynamic_cast<ReconstructedParticle*>(leptonCol->getElementAt(i));
	TLorentzVector lepVec(TVector3(mylepton->getMomentum()),mylepton->getEnergy());
	_lepCharge += mylepton->getCharge();
	*_myLepton += lepVec;
      }

    // Initialise and write : JETS
    LCCollection* jetCol = evt->getCollection(_myJetsCollection);
    for(int i = 0 ; i < 4 ; i++)
      {
	_jets[i]->SetPxPyPzE(0.,0.,0.,0.);
	_jetsBtag[i] = 0.;
	_jetsCtag[i] = 0.;
	_jetsVertexCharge[i] = 10; // 10 means no vertex found
	_jetsVertexChargeTer[i] = 10; // 10 means no vertex found
	_jetsTrueFlavour[i] = 0;
	_jetsQuarkOrigin[i] = 0;  // 0 for uds, 1 for c and 2 for b
	_jetsTrueCharge[i] = 0;
	_jetsTrueBmesonCharge[i] = 3; // 3 means not defined charge
	_jetsDecayLength[i] = -1; // -1 means no vertex to calculate the decay length
	_jetsDecayLengthSig[i] = -1; // -1 means no vertex to calculate the decay length significance
	_jetsNumofTracks[i] = -1; // -1 means no vertex to calculate the number of tracks
	_jetsNumVertices[i] = 0;
	_jetsVtxMom[i] = -1; // -1 means no vertex
	_jetsVtxPtmass[i] = -1; // -1 means no vertex
	_NotIsoMuonFound[i] = false;
	_NotIsoMuons_energy[i] = 0;
	_NotIsoMuons_z4[i] = 0;
	_NotIsoMuons_pT[i] = 0;
	_NotIsoMuons_xT[i] = 0;
	_NotIsoMuons_charge[i] = 0;
      }
    _yplus = 0.;
    _yminus = 0.;
    _isGoodBs = false;
    _isGoodComb = false;
    _chi2DBD = 0;
    _chi2new = 0;
    _eventcharge = 0;
    _isLeptonVtxChargeOK = false;

    // Links PFOs - MCtruth
    LCCollection *mctruth = evt->getCollection("RecoMCTruthLink");
    LCRelationNavigator findPFOs(mctruth);
     
    if(jetCol->getNumberOfElements() == 4)
      {
	_yplus = jetCol->getParameters().getFloatVal("YPlus");
	_yminus = jetCol->getParameters().getFloatVal("YMinus");
	//get the PIDHandler fot this jet collection
	PIDHandler pidh(jetCol);
	//get the algorithm id for the Flavour Tag info
	int alid = pidh.getAlgorithmID("lcfiplus");
	int ibtag = pidh.getParameterIndex(alid,"BTag");
	int ictag = pidh.getParameterIndex(alid,"CTag");
	int iDLtag = pidh.getParameterIndex(alid,"vtxlen1"); //decay length of the first vertex in the jet (zero if no vertex is found) 
	int iDLStag = pidh.getParameterIndex(alid,"vtxsig1"); //decay length significance of the first vertex in the jet (zero if no vertex is found) 
	int iMulttag = pidh.getParameterIndex(alid,"vtxmult1"); //number of tracks included in the first vertex (zero if no vertex is found) 
	int iNVtxtag = pidh.getParameterIndex(alid,"nvtxall"); // number of vertices in a jet
	int iVtxMtag = pidh.getParameterIndex(alid,"vtxmom1"); // magnitude of the vector sum of the momenta of all tracks combined into the first vertex
	int iVtxPttag = pidh.getParameterIndex(alid,"vtxmasspc"); // mass of the vtx with minimum pt correction allowed by the error matrices of the primary and secondary vtx
	// intermediate variables
	TLorentzVector jets[4];
	float bTag[4], cTag[4];
	float jetBangle[4], jetBbarangle[4], jetDecayLength[4], jetDecayLengthSig[4], jetVtxMom[4], jetVtxPtmass[4], jetQuarkOrigin[4];
	int jetType[4], jetCharge[4], VertexCharge[4], VertexChargeTer[4], jetBmesonCharge[4], jetNumofTracks[4], jetNumVertices[4], notIsoMuons_charge[4];
	float notIsoMuons_energy[4], notIsoMuons_z4[4], notIsoMuons_pT[4], notIsoMuons_xT[4];
	bool notIsoMuonFound[4];
	int Nmeson_found = 0;
	bool twoMesons_found = false;

	for(int i = 0 ; i < 4 ; i++)
	  {
	    ReconstructedParticle *a_jet = dynamic_cast<ReconstructedParticle*>(jetCol->getElementAt(i));
	    TLorentzVector Jet4vec(TVector3(a_jet->getMomentum()),a_jet->getEnergy());
	    jets[i] = Jet4vec;
	    notIsoMuons_energy[i] = 0;
	    notIsoMuons_z4[i] = 0;
	    notIsoMuons_pT[i] = 0;
	    notIsoMuons_xT[i] = 0;
	    notIsoMuonFound[i] = false;
	    notIsoMuons_charge[i] = 0;
	    float jet_QuarkOrigin = 0; // 0 for uds, 1 for c and 2 for b
	    ReconstructedParticleVec jetParticles = a_jet->getParticles();
	    for(int k = 0; k < (int)jetParticles.size(); k++)
	      {
		TLorentzVector thisPFO(TVector3(jetParticles.at(k)->getMomentum()),jetParticles.at(k)->getEnergy());
		
		int object_QuarkOrigin = -1; // 0 for uds, 1 for c and 2 for b
		if(findPFOs.getRelatedToObjects(jetParticles.at(k)).size() > 0)
		  {
		    double maxWeight = 0;
		    int good_relatedObject = 0;	
		    for(unsigned int j = 0; j < findPFOs.getRelatedToObjects(jetParticles.at(k)).size(); j++)
		      {
			if(findPFOs.getRelatedToWeights(jetParticles.at(k)).at(j) > maxWeight)
			  {
			    maxWeight = findPFOs.getRelatedToWeights(jetParticles.at(k)).at(j);
			    good_relatedObject = j;
			  }
		      }
		    MCParticle *mcObj = dynamic_cast<MCParticle*>( findPFOs.getRelatedToObjects(jetParticles.at(k)).at(good_relatedObject) );
		    
		    while(mcObj->getParents().size() != 0 && object_QuarkOrigin == -1)
		      {
			mcObj = mcObj->getParents()[0];
			if(mcObj->getParents().size() != 0)
			  {
			    for(unsigned int parent = 0; parent < mcObj->getParents().size(); parent++)
			      { 
				MCParticle *mcparents = mcObj->getParents()[parent];
			
				if(abs(mcparents->getPDG()) <= 3) object_QuarkOrigin = 0; 
				if(abs(mcparents->getPDG()) == 4) object_QuarkOrigin = 1; 
				if(abs(mcparents->getPDG()) == 5) object_QuarkOrigin = 2; 
			      }
			  }
		      }
		    //std::cout << "Jet " << i << " object_QuarkOrigin = " << object_QuarkOrigin << std::endl;

		    if( object_QuarkOrigin != -1) jet_QuarkOrigin += object_QuarkOrigin;
		  }

		// ***********  Apply our LeptonID on all objects with a minimum energy = _leptonEmin to find not isolated muons  *******************
		// get clusters
		float Eecal = 0; // ecal energy deposition
		float Etotal = 0;  // total (ecal+hcal) energy deposition
		ClusterVec aclustervec = jetParticles.at(k)->getClusters();
		for(int cluster_i = 0; cluster_i < (int)aclustervec.size(); cluster_i++)
		  {
		    Cluster* acluster = aclustervec.at(cluster_i);
		    Eecal += acluster->getSubdetectorEnergies()[0];
		    Eecal += acluster->getSubdetectorEnergies()[3];
		    Etotal += acluster->getEnergy();
		  }

		// get tracks momentum - if no track associated : TrackMom = 0
		float TrackMom = 0; // track momentum
		TrackVec atrackvec = jetParticles.at(k)->getTracks();
		if(atrackvec.size() >= 1)
		  {
		    Track *a_track = atrackvec[0];
		    float tkd0 = a_track->getD0();
		    float tkphi = a_track->getPhi();
		    float tkomega = a_track->getOmega();
		    float tkz0 = a_track->getZ0();
		    float tktanlambda = a_track->getTanLambda();
		    float bField = Global::GEAR->getBField().at( gear::Vector3D( 0., 0., 0.) ).z();
		    // fit with an helix
		    HelixClass* helix = new HelixClass();
		    helix->Initialize_Canonical( tkphi , tkd0 , tkz0 , tkomega , tktanlambda , bField );
		    TVector3 ParTrack( helix->getMomentum() );
		    // Use real momentum instead of energy in the case of particles with one track
		    TrackMom = ParTrack.Mag();
		  }

		// Change PDG if the particle satisfies the criteria :
		// 0 < Eecal/Etotal < 0.6 && 0 < Etotal/TrackMom < 0.6 for muons
		// 0.9 < Eecal/Etotal < 1 && 0.2 < Etotal/TrackMom < 1.5 for electrons
		int object_PDG = 0;
		object_PDG = jetParticles.at(k)->getType();
		if(TrackMom > 5)
		  {
		    int sign = 1;
		    if(jetParticles.at(k)->getCharge() > 0) sign = -1;
		    if(Etotal/TrackMom > 0 && Etotal/TrackMom < 0.6 && Eecal/Etotal > 0 && Eecal/Etotal < 0.6 ) object_PDG = 13*sign;		    
		    else if(Etotal/TrackMom > 0.2 && Etotal/TrackMom < 1.5 && Eecal/Etotal > 0.9 && Eecal/Etotal < 1) object_PDG = 11*sign;
		  }
		
		float object_energy = thisPFO.E();
		// Calculate z4, pT4 and xT
		float z4 = thisPFO.E()/Jet4vec.E();
		float theta = TMath::Abs(thisPFO.Angle(Jet4vec.Vect()));
		float pT = thisPFO.P()*sin(theta);
		float xT = pT/Jet4vec.M();
	
		// ********** ISOLATION **********
	
		if(z4<0.6 && xT<0.25 &&  object_energy > 5 && abs(object_PDG) == 13 && object_energy > notIsoMuons_energy[i])
		  {
		    std::cout << " *********************** Not isolated muons found in jet ***************************** " << i << std::endl; 
		    notIsoMuonFound[i] = true;
		    notIsoMuons_energy[i] = object_energy;
		    notIsoMuons_z4[i] = z4;
		    notIsoMuons_pT[i] = pT;
		    notIsoMuons_xT[i] = xT;
		    notIsoMuons_charge[i] = jetParticles.at(k)->getCharge();
		  }
			
	      }
	    
	    jet_QuarkOrigin = jet_QuarkOrigin/(double)jetParticles.size();
	    std::cout << "Jet " << i << " QuarkOrigin = " << jet_QuarkOrigin << std::endl;
	    jetQuarkOrigin[i] = jet_QuarkOrigin;

	    //get the Particle id object containing the Flavour Tag info
	    const ParticleID& pid = pidh.getParticleID(a_jet,alid);
	    //get the actual flavour tags
	    bTag[i] = pid.getParameters()[ibtag];
	    cTag[i] = pid.getParameters()[ictag];
	    jetDecayLength[i] = pid.getParameters()[iDLtag];
	    jetDecayLengthSig[i] = pid.getParameters()[iDLStag];
	    jetNumofTracks[i] = (int)pid.getParameters()[iMulttag];
	    jetNumVertices[i] = (int)pid.getParameters()[iNVtxtag];
	    jetVtxMom[i] = pid.getParameters()[iVtxMtag];
	    jetVtxPtmass[i] = pid.getParameters()[iVtxPttag];
	    //get the angles to know the true flavour
	    jetBangle[i] = fabs(Jet4vec.Angle(_MCb->Vect()));
	    jetBbarangle[i] = fabs(Jet4vec.Angle(_MCbbar->Vect()));


	    // Vertex charge
	    //Link vertex - jets collection
	    LCCollection *vtxjetrelCol = evt->getCollection("FinalJets_rel");
	    LCRelationNavigator findVertex(vtxjetrelCol);
	  
	    VertexCharge[i] = 10;
	    jetBmesonCharge[i] = 3;
	    if(findVertex.getRelatedToObjects(a_jet).size() > 0)
	      {
		double maxWeight = 0;
		int good_relatedObject = 0;	
		for(unsigned int j = 0; j < findVertex.getRelatedToObjects(a_jet).size(); j++)
		  {
		    if(findVertex.getRelatedToWeights(a_jet).at(j) > maxWeight)
		      {
			maxWeight = findVertex.getRelatedToWeights(a_jet).at(j);
			good_relatedObject = j;
		      }
		  }
		Vertex *jetVertex = dynamic_cast<Vertex*>( findVertex.getRelatedToObjects(a_jet).at(good_relatedObject) );
		ReconstructedParticle *Vertexpart = dynamic_cast<ReconstructedParticle*> (jetVertex->getAssociatedParticle());
		VertexCharge[i] = (int)Vertexpart->getCharge();
		if(findVertex.getRelatedToObjects(a_jet).size() == 2)
		  {
		    Vertex *jetVertex2 = dynamic_cast<Vertex*>( findVertex.getRelatedToObjects(a_jet).at(1) );
		    ReconstructedParticle *Vertexpart2 = dynamic_cast<ReconstructedParticle*> (jetVertex2->getAssociatedParticle());
		    VertexChargeTer[i] = (int)Vertexpart2->getCharge();
		  }
		std::cout << " Jet " << i << " : Btag = " << bTag[i] << " Charge = " << VertexCharge[i] << std::endl;
		bool meson_found = false;
		for(unsigned int k = 0; k < Vertexpart->getParticles().size(); k++)
		  {   
		    if(findPFOs.getRelatedToObjects(Vertexpart->getParticles()[k]).size() > 0 && meson_found == false)
		      {
			maxWeight = 0;
			good_relatedObject = 0;	
			for(unsigned int j = 0; j < findPFOs.getRelatedToObjects(Vertexpart->getParticles()[k]).size(); j++)
			  {
			    if(findPFOs.getRelatedToWeights(Vertexpart->getParticles()[k]).at(j) > maxWeight)
			      {
				maxWeight = findPFOs.getRelatedToWeights(Vertexpart->getParticles()[k]).at(j);
				good_relatedObject = j;
			      }
			  }
			MCParticle *mcObj = dynamic_cast<MCParticle*>(findPFOs.getRelatedToObjects(Vertexpart->getParticles()[k]).at(good_relatedObject));
			MCParticle *mcParent = mcObj;
			if(twoMesons_found) std::cout << "Already two  B mesons found" << std::endl; 
			while(!meson_found  && !twoMesons_found)
			  {
			    if(mcObj->getParents().size() != 0)
			      {
				for(unsigned int parent = 0; parent < mcObj->getParents().size(); parent++)
				  {
				    mcParent = dynamic_cast<MCParticle*>(mcObj->getParents()[parent]);
				    int PDG = mcParent->getPDG();
				    if(abs(PDG)==511||abs(PDG)==10511||abs(PDG)==513||abs(PDG)==10513||abs(PDG)==20513||abs(PDG)==515||abs(PDG)==531||abs(PDG)==10531||abs(PDG)==533||abs(PDG)==10533||abs(PDG)==20533 || abs(PDG)==535 || abs(PDG)==5122 || abs(PDG)==5232)
				      {
					std::cout << "B0 meson jet" << std::endl;
					meson_found = true;
					jetBmesonCharge[i] = 0;
				      }
				    if(abs(PDG)==521||abs(PDG)==10521||abs(PDG)==523||abs(PDG)==10523||abs(PDG)==20523||abs(PDG)==525||abs(PDG)==541||abs(PDG)==10541||abs(PDG)==543||abs(PDG)==10543||abs(PDG)==545 || abs(PDG)==20543 || abs(PDG)==5132 || abs(PDG)==5332)
				      {
					if(PDG > 0) 
					  {
					    std::cout << "B+ meson jet" << std::endl;
					    meson_found = true;
					    jetBmesonCharge[i] = 1;
					  }
					if(PDG < 0)
					  {
					    std::cout << "B- meson jet" << std::endl;
					    meson_found = true;
					    jetBmesonCharge[i] = -1;
					  }
				      }
				  }
				if(!meson_found) mcObj = dynamic_cast<MCParticle*>(mcObj->getParents()[0]);
				if(meson_found) Nmeson_found++;
				if(Nmeson_found == 2) twoMesons_found = true;
			      }
			    if(mcObj->getParents().size() == 0)
			      {
				std::cout << "No B meson found" << std::endl;
				break;
			      }
			  }				  
		      }
		  }
	      }
	  }   // end loop on the 4 jets

	// select the smallest jetBangle
	int smallestjetBanglePosition = 0;
	int smallestjetBbaranglePosition = 0;
	for(int i = 1; i < 4; i++)
	  {
	    if(jetBangle[i] < jetBangle[smallestjetBanglePosition]) smallestjetBanglePosition = i;
	  }
	// select the smallest jetBbarangle among 3 other jets
	for(int i = 1; i < 4; i++)
	  {
	    if(i == smallestjetBanglePosition) continue;
	    if(jetBbarangle[i] < jetBbarangle[smallestjetBbaranglePosition]) smallestjetBbaranglePosition = i;
	  }

	if(jetBbarangle[smallestjetBanglePosition] < jetBbarangle[smallestjetBbaranglePosition]) // do the opposite order
	  {
	    smallestjetBanglePosition = 0;
	    smallestjetBbaranglePosition = 0;
	    // select the smallest jetBbarangle
	    for(int i = 1; i < 4; i++)
	      {
		if(jetBbarangle[i] < jetBbarangle[smallestjetBbaranglePosition]) smallestjetBbaranglePosition = i;
	      }
	    // select the smallest jetBangle among 3 other jets
	    for(int i = 1; i < 4; i++)
	      {
		if(i == smallestjetBbaranglePosition) continue;
		if(jetBangle[i] < jetBangle[smallestjetBanglePosition]) smallestjetBanglePosition = i;
	      }
	  }
 

	// Set the "almost" true flavour and charge to the B jets
	for(int i = 0 ; i < 4 ; i++) 
	  {
	    if(i == smallestjetBanglePosition)
	      {
		jetType[i] = 5;
		jetCharge[i] = -1;
	      }
	    else if(i == smallestjetBbaranglePosition)
	      {
		jetType[i] = 5;
		jetCharge[i] = 1;
	      }
	    else 
	      {
		jetType[i] = 0;
		jetCharge[i] = 0;
	      }
	  }


	// sort jets by btag
	int sortPositions[4];
	float tempBtag[4];
	for(int i = 0 ; i < 4 ; i++) 
	  {
	    tempBtag[i] = bTag[i];
	    sortPositions[i] = i;
	  }
	for(int i = 0 ; i < 4 ; i++)
	  {
	    int highestBtagPosition = i;
	    for(int j = i ; j < 4 ; j++)
	      {
		if(tempBtag[j] > tempBtag[highestBtagPosition]) highestBtagPosition = j;
	      }
	    float temp = tempBtag[i];
	    tempBtag[i] = tempBtag[highestBtagPosition];
	    tempBtag[highestBtagPosition] = temp;
	    int tempInt = sortPositions[i];
	    sortPositions[i] = sortPositions[highestBtagPosition];
	    sortPositions[highestBtagPosition] = tempInt;
	  }
	
	// write sorted jets in the tree
	for(int i = 0 ; i < 4 ; i++)
	  {
	    *_jets[i] = jets[sortPositions[i]];
	    _jetsBtag[i] = bTag[sortPositions[i]];
	    _jetsCtag[i] = cTag[sortPositions[i]];
	    _jetsVertexCharge[i] = VertexCharge[sortPositions[i]];
	    _jetsVertexChargeTer[i] = VertexChargeTer[sortPositions[i]];
	    _jetsTrueFlavour[i] = jetType[sortPositions[i]];
	    _jetsQuarkOrigin[i] = jetQuarkOrigin[sortPositions[i]];
	    _jetsTrueCharge[i] = jetCharge[sortPositions[i]];
	    _jetsTrueBmesonCharge[i] = jetBmesonCharge[sortPositions[i]];
	    _jetsDecayLength[i] = jetDecayLength[sortPositions[i]];
	    _jetsDecayLengthSig[i] = jetDecayLengthSig[sortPositions[i]];
	    _jetsNumofTracks[i] = jetNumofTracks[sortPositions[i]];
	    _jetsNumVertices[i] = jetNumVertices[sortPositions[i]];
	    _jetsVtxMom[i] = jetVtxMom[sortPositions[i]];
	    _jetsVtxPtmass[i] = jetVtxPtmass[sortPositions[i]];
	    _NotIsoMuonFound[i] = notIsoMuonFound[sortPositions[i]];
	    _NotIsoMuons_energy[i] = notIsoMuons_energy[sortPositions[i]];
	    _NotIsoMuons_z4[i] = notIsoMuons_z4[sortPositions[i]];
	    _NotIsoMuons_pT[i] = notIsoMuons_pT[sortPositions[i]];
	    _NotIsoMuons_xT[i] = notIsoMuons_xT[sortPositions[i]];
	    _NotIsoMuons_charge[i] = notIsoMuons_charge[sortPositions[i]];
	  }
  
	if(_jetsTrueFlavour[0] == 5 && _jetsTrueFlavour[1] == 5) _isGoodBs = true;

	// Reconstruct the top
	*_Whad = *_jets[2] + *_jets[3];
	_selectedBjet = -1;
	_top_had->SetPxPyPzE(0.,0.,0.,0.);
	TLorentzVector top1 = *_Whad + *_jets[0];
	TLorentzVector top2 = *_Whad + *_jets[1];
	Double_t bpstar1 = top1.Gamma()*_jets[0]->P()*(1 - top1.Beta()*TMath::Cos(top1.Angle(_jets[0]->Vect())));
	Double_t bpstar2 = top2.Gamma()*_jets[1]->P()*(1 - top2.Beta()*TMath::Cos(top2.Angle(_jets[1]->Vect())));
	Double_t WBangle1 = TMath::Cos(_Whad->Angle(_jets[0]->Vect()));
	Double_t WBangle2 = TMath::Cos(_Whad->Angle(_jets[1]->Vect()));
	Double_t EstarW1 = (top1.M2() + _Whad->M2()) / (2 * top1.M()); 
	Double_t EstarW2 = (top2.M2() + _Whad->M2()) / (2 * top2.M()); 

	if(pow(top1.M()-174,2)/pow(6.3,2) + pow(top1.E()-250,2)/pow(8.0,2) + pow(bpstar1-68,2)/pow(5,2) + pow(WBangle1-0.23,2)/pow(0.14,2) < pow(top2.M()-174,2)/pow(6.3,2) + pow(top2.E()-250,2)/pow(8.0,2) + pow(bpstar2-68,2)/pow(5,2) + pow(WBangle2-0.23,2)/pow(0.14,2))
	  {
	    *_top_had = top1;
	    _selectedBjet = 0;
	    if(_isGoodBs && _jetsTrueCharge[0] == -_MCtop_jet_Q) _isGoodComb = true;
	    _chi2DBD = pow(bpstar1-68,2)/pow(5,2) + pow(top1.Gamma()-1.435,2)/pow(0.05,2) + pow(WBangle1-0.23,2)/pow(0.14,2);
	    _chi2new = pow(bpstar1-68,2)/pow(8.5,2) + pow(top1.M()-174,2)/pow(15,2) + pow(EstarW1-105,2)/pow(10,2);
	  }
	else
	  {
	    *_top_had = top2;
	    _selectedBjet = 1;
	    if(_isGoodBs && _jetsTrueCharge[1] == -_MCtop_jet_Q) _isGoodComb = true;
	    _chi2DBD = pow(bpstar2-68,2)/pow(5,2) + pow(top2.Gamma()-1.435,2)/pow(0.05,2) + pow(WBangle2-0.23,2)/pow(0.14,2);
	    _chi2new = pow(bpstar2-68,2)/pow(8.5,2) + pow(top2.M()-174,2)/pow(15,2) + pow(EstarW2-105,2)/pow(10,2);
	  }
	_eventcharge = _jetsVertexCharge[0] - _jetsVertexCharge[1];
	if((_selectedBjet == 0 && _eventcharge/_lepCharge > 0) || (_selectedBjet == 1 && _eventcharge/_lepCharge < 0)) _isLeptonVtxChargeOK = true;
      }
    _tree->Fill();
  }
  
  catch(DataNotAvailableException &e) { }
  
}



void ttbarTreeWriter::check(LCEvent * evt)
{ 
  
}


void ttbarTreeWriter::end()
{ 
  ROOTfile->Write();
  ROOTfile->Close();
  delete ROOTfile;
}
