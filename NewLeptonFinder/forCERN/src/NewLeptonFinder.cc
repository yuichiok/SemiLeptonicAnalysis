#include "NewLeptonFinder.h"
#include <iostream>
#include <math.h>

#include <IMPL/LCCollectionVec.h>
#include <EVENT/ReconstructedParticle.h>
#include <IMPL/ReconstructedParticleImpl.h>
#include <EVENT/ParticleID.h>
#include <IMPL/ParticleIDImpl.h>
#include <UTIL/PIDHandler.h>
#include <LCIOSTLTypes.h>
#include <EVENT/MCParticle.h>
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

NewLeptonFinder aNewLeptonFinder;

NewLeptonFinder::NewLeptonFinder() : Processor("NewLeptonFinder")
{
  
  // modify processor description

  _description = "Processor made for leptonic W decays : it takes the best PFO lepton (electron/muon) according to isolation criteria";
  
  
  // register steering parameters: name, description, class-variable, default value

  registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
                           "ReconstructedParticle" ,
                           "Name of the PFOs collection"  ,
                           _InPFOscolName ,
                           std::string("PandoraPFOs") );

  registerProcessorParameter( "4jetsCollection" ,
			      "Name of the 4 jets collection" ,
			      _col4jetsName ,
			      std::string("RefinedJets") );
  
  registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			    "PFOsCollection",
			    "Pandora PFOs Collection Name (without the selected lepton)",
			    _OutPFOscolName,
			    std::string("NewPandoraPFOs"));
  
  registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE,
			    "LeptonCollection",
			    "Isolated lepton Collection Name",
			    _Lepton,
			    std::string("SelectedLepton"));

  registerProcessorParameter( "ROOTFileName",
                              "Output ROOT File Name",
                              _hfilename,
                              std::string("LeptonFinder.root") );

  registerProcessorParameter( "LeptonEnergyCut",
                              "Cut on the lepton minimal energy",
                              _leptonEmin,
                              float(5) );

  registerProcessorParameter( "UseOwnPID",
                              "Replace PFA ID by our own lepton ID",
                              _useLeptonID,
                              bool(true) );
 
}


void NewLeptonFinder::init()
{ 
  
  _nRun = 0;
  _nEvt = 0;

  ROOTfile = new TFile( _hfilename.c_str(), "RECREATE", _hfilename.c_str() );
  
  _tree = new TTree( "tree", "tree" );
  
  _tree->Branch("eventNumber", &_nEvt, "eventNumber/I");

  _MClepton = new TLorentzVector();
  _tree->Branch("MClepton","TLorentzVector",&_MClepton,16000,0);
  _tree->Branch("MCz4", &_MCz4 , "MCz4/F");
  _tree->Branch("MCpT", &_MCpT , "MCpT/F");
  _tree->Branch("MCxT", &_MCxT , "MCxT/F");
  _tree->Branch("isMCLeptonFound", &_isMCLeptonFound , "isMCLeptonFound/O");
  _tree->Branch("isMCLeptonFromTau", &_isMCLeptonFromTau , "isMCLeptonFromTau/O");
  _tree->Branch("MCleptonPosition", &_posMCpfo , "MCleptonPosition/I");
  _tree->Branch("MCleptonPDG", &_MCpdg , "MCleptonPDG/I");

  _tree->Branch("numberOfPFOs", &_nObj, "numberOfPFOs/I");
  _tree->Branch("PFO_energy", _object_energy, "PFO_energy[numberOfPFOs]/F");
  _tree->Branch("PFO_P", _object_momentum, "PFO_P[numberOfPFOs]/F");
  _tree->Branch("PFO_px", _object_px, "PFO_px[numberOfPFOs]/F");
  _tree->Branch("PFO_py", _object_py, "PFO_py[numberOfPFOs]/F");
  _tree->Branch("PFO_pz", _object_pz, "PFO_pz[numberOfPFOs]/F");
  _tree->Branch("PFO_phi", _object_phi, "PFO_phi[numberOfPFOs]/F");
  _tree->Branch("PFO_theta", _object_theta, "PFO_theta[numberOfPFOs]/F");
  _tree->Branch("PFO_PDG", _object_PDG, "PFO_PDG[numberOfPFOs]/I");
  _tree->Branch("PFO_mcPDG", _object_mcPDG, "PFO_mcPDG[numberOfPFOs]/I");
  _tree->Branch("PFO_jetNb", _object_jetNb, "PFO_jetNb[numberOfPFOs]/I");
  _tree->Branch("Eecal", _Eecal, "Eecal[numberOfPFOs]/F");
  _tree->Branch("Ehcal", _Ehcal, "Ehcal[numberOfPFOs]/F");
  _tree->Branch("Eyoke", _Eyoke, "Eyoke[numberOfPFOs]/F");
  _tree->Branch("Ebcal", _Ebcal, "Ebcal[numberOfPFOs]/F");
  _tree->Branch("Etotal", _Etotal, "Etotal[numberOfPFOs]/F");
  _tree->Branch("Momentum", _Momentum, "Momentum[numberOfPFOs]/F");
  _tree->Branch("TrackMom", _TrackMom, "TrackMom[numberOfPFOs]/F");
  _tree->Branch("z4", _z4 , "z4[numberOfPFOs]/F");
  _tree->Branch("pT", _pT , "pT[numberOfPFOs]/F");
  _tree->Branch("xT", _xT , "xT[numberOfPFOs]/F");
  _tree->Branch("PFO_quarkOrigin", _object_QuarkOrigin, "PFO_quarkOrigin[numberOfPFOs]/I");
  _tree->Branch("PFOfromtau", _objectFromtau , "PFOfromtau[numberOfPFOs]/O");
  _tree->Branch("PFOfromB", _objectFromBmeson , "PFOfromB[numberOfPFOs]/O");
  _jet0 = new TLorentzVector();
  _jet1 = new TLorentzVector();
  _jet2 = new TLorentzVector();
  _jet3 = new TLorentzVector();
  _tree->Branch("jet0","TLorentzVector",&_jet0,16000,0);
  _tree->Branch("jet1","TLorentzVector",&_jet1,16000,0);
  _tree->Branch("jet2","TLorentzVector",&_jet2,16000,0);
  _tree->Branch("jet3","TLorentzVector",&_jet3,16000,0);
  
  _tree->Branch("selectedLeptonPosition", &_leptonNb, "selectedLeptonPosition/I");
  _Recolepton = new TLorentzVector();
  _tree->Branch("Recolepton","TLorentzVector",&_Recolepton,16000,0);
  _tree->Branch("Recoz4", &_Recoz4 , "Recoz4/F");
  _tree->Branch("RecopT", &_RecopT , "RecopT/F");
  _tree->Branch("RecoxT", &_RecoxT , "RecoxT/F");

  _tree->Branch("NbOfLeptons", &_nbOfLep , "NbOfLeptons/I");
  _tree->Branch("candidate_PFO", _candidate_PFO , "candidate_PFO[NbOfLeptons]/I");
  _tree->Branch("NbOfPhotons", &_nbPhotons , "NbOfPhotons/I");
  _tree->Branch("isGoodLeptonFound", &_isGoodLeptonFound , "isGoodLeptonFound/I");
  _tree->Branch("lepFromtau", &_lepFromtau , "lepFromtau/I");
  _Hadrons = new TLorentzVector();
  _tree->Branch("Hadrons","TLorentzVector",&_Hadrons,16000,0);
  _tree->Branch("NbOfNotIsoLeptons", &_nbNotIsoLep , "NbOfNotIsoLeptons/I");
  _tree->Branch("NotIsoLepcharge", _NotIsoLepcharge , "NotIsoLepcharge[NbOfNotIsoLeptons]/I");
  _tree->Branch("NotIsoLepPDG", _NotIsoLepPDG , "NotIsoLepPDG[NbOfNotIsoLeptons]/I");
}


void NewLeptonFinder::processRunHeader(LCRunHeader* run)
{ 
  _nRun++;
} 


void NewLeptonFinder::processEvent(LCEvent * evt)
{  
  _nEvt ++;
  if(_nEvt%50 == 0) std::cout << "*************  Event " << _nEvt << "  *******************" << std::endl;
  ROOTfile->cd();

  try {
    
    // Getting the information on the Reconstructed Objects
    LCCollection* recobj = evt->getCollection( _InPFOscolName );
    _nObj = recobj->getNumberOfElements();
    if(_nObj > 1000) std::cout << "Event : " << _nEvt << " -> " << _nObj << " objects !!!" << std::endl;
    // Get 4-jets reconstructed collection
    LCCollection *fourjetscol = evt->getCollection(_col4jetsName); 
    // Links to the corresponding PFOs
    LCCollection *mctruth = evt->getCollection("RecoMCTruthLink");
    LCRelationNavigator findPFOs(mctruth);

   
    // ********** MONTE CARLO STUFF **********
    
    // Initialisation
    _MClepton->SetPxPyPzE(0.,0.,0.,0.);
    _MCz4 = 0; _MCpT = 0; _MCxT = 0;
    _posMCpfo = -2; // -1 is for _leptonNb init
    _MCpdg = 0;
    _isMCLeptonFound = false;
    _isMCLeptonFromTau = false;
    
    // Get information on the MCParticle
    LCCollection* mccol = evt->getCollection("MCParticlesSkimmed");
    bool leptonSeen = false;
    bool leptonTwoSeen = false;
    MCParticle *relationlepton;
  
    if(mccol->getNumberOfElements() > 11)
      {
	int mcobj_num = 6;
	while(!leptonSeen && mcobj_num < mccol->getNumberOfElements())
	  {
	    MCParticle *mcpart = dynamic_cast<MCParticle*>( mccol->getElementAt(mcobj_num) );
	    //std::cout << "PDG: " <<(mcpart->getPDG())  << " charge: " << mcpart->getCharge() << "\n";
	    if((abs(mcpart->getPDG()) == 11 || abs(mcpart->getPDG()) == 13 || abs(mcpart->getPDG()) == 15) && mcpart->getGeneratorStatus() == 2)
	      {
		if(mcpart->getParents().size() == 2 && mcpart->getDaughters().size() == 1)
		  {
		    if(abs(mcpart->getParents()[0]->getPDG())  == 11 && abs(mcpart->getParents()[1]->getPDG())  == 11)
		      {
			TLorentzVector mcVec(TVector3(mcpart->getMomentum()),mcpart->getEnergy());
			*_MClepton = mcVec;
			leptonSeen = true;
			_MCpdg = mcpart->getPDG();
			if(abs(mcpart->getPDG()) == 15) _isMCLeptonFromTau = true;   
		      }
		  }
	      }
	    mcobj_num++;
	    if(mcobj_num == mccol->getNumberOfElements()) std::cout << "Warning no MC lepton found" << std::endl;
	  }
	
	if(leptonSeen)
	  {
	    _isMCLeptonFound = true;

	    for(int i = 12 ; i < mccol->getNumberOfElements() ; i++)
		  {
		    if(!leptonTwoSeen)
		      {
			MCParticle *mcpart = dynamic_cast<MCParticle*>( mccol->getElementAt(i) );
			if(mcpart->getPDG() == _MCpdg && (mcpart->hasLeftDetector() || mcpart->isDecayedInCalorimeter() || mcpart->isDecayedInTracker()))
			  {
			    leptonTwoSeen = true;
			    relationlepton = dynamic_cast<MCParticle*>( mccol->getElementAt(i) );
			  }
		      }
		  }
	    if(leptonTwoSeen)
	      {
		if(findPFOs.getRelatedFromObjects(relationlepton).size() != 0)
		  {
		    double maxWeight = 0;
		    int good_relatedObject = 0;	
		    for(unsigned int i = 0; i < findPFOs.getRelatedFromObjects(relationlepton).size(); i++)
		      {
			if(findPFOs.getRelatedFromWeights(relationlepton).at(i) > maxWeight)
			  {
			    maxWeight = findPFOs.getRelatedFromWeights(relationlepton).at(i);
			    good_relatedObject = i;
			  }
		      }
		    int pos = 0;
		    ReconstructedParticle *recolepton = dynamic_cast<ReconstructedParticle*>( findPFOs.getRelatedFromObjects(relationlepton).at(good_relatedObject) );
		    for(int j = 0 ; j < fourjetscol->getNumberOfElements() ; j++)
		      {
			ReconstructedParticle *recojet = dynamic_cast<ReconstructedParticle*>( fourjetscol->getElementAt(j) );
			ReconstructedParticleVec jetParticles = recojet->getParticles();
			for(int k = 0 ; k < (int)jetParticles.size() ; k++)
			  {
			    if(jetParticles.at(k)->id() == recolepton->id())
			      {
				_posMCpfo = pos;
				TLorentzVector jetVec(TVector3(recojet->getMomentum()),recojet->getEnergy());
				_MCz4 = _MClepton->E()/jetVec.E();
				float theta = TMath::Abs(_MClepton->Angle(jetVec.Vect()));
				_MCpT = _MClepton->P()*sin(theta);
				_MCxT = _MCpT/jetVec.M();
			      }
			    pos++;
			  }
		      }
		  }
	      }
	  }
      }
    
  
    // ********** END OF MONTE CARLO STUFF **********
        
    /*
***************
Find the lepton
***************
*/    
    
    // Initialise the reconstructed lepton
    _Recolepton->SetPxPyPzE(0.,0.,0.,0.);
    _Recoz4 = 0.;
    _RecopT = 0.;
    _RecoxT = 0.;
    _leptonNb = -1;
    _nbPhotons = 0;
    _isGoodLeptonFound = 0;
    _lepFromtau = 0;
    _Hadrons->SetPxPyPzE(0.,0.,0.,0.);
    _nbNotIsoLep = 0;

    // Looks for the best suitable reconstructed lepton (electron/muon)
    // Then copies information into ROOT Tree
    std::vector<int> leptonsPositions;
      
    // Loop on Jets and their inner Reconstructed objects
    _jet0->SetPxPyPzE(0.,0.,0.,0.);
    _jet1->SetPxPyPzE(0.,0.,0.,0.);
    _jet2->SetPxPyPzE(0.,0.,0.,0.);
    _jet3->SetPxPyPzE(0.,0.,0.,0.);
  
    int pfopos = 0;
    int thispos = 0;
    bool usedphoton[nObjects]; // to tag the already used photons to cure bremsstrahlung and FSR 
    for(int i = 0; i < nObjects; i++)
      {
	usedphoton[i] = false;
      }

    for(int i = 0; i < fourjetscol->getNumberOfElements(); i++)
      {
	ReconstructedParticle *recojet = dynamic_cast<ReconstructedParticle*>( fourjetscol->getElementAt(i) );
	TLorentzVector Jet4vec(TVector3(recojet->getMomentum()),recojet->getEnergy());
	if(i == 0) *_jet0 = Jet4vec; 
	else if(i == 1) *_jet1 = Jet4vec;
	else if(i == 2) *_jet2 = Jet4vec;
	else if(i == 3) *_jet3 = Jet4vec;

	ReconstructedParticleVec jetParticles = recojet->getParticles();
	for(int k = 0; k < (int)jetParticles.size(); k++)
	  {
	    _lepton_pos[pfopos] = -1;
	    for(int phNb = 0; phNb < 5 ; phNb++)
	      { 
		_photon_pos[pfopos][phNb] = -1;
	      }

	    TLorentzVector thisPFO(TVector3(jetParticles.at(k)->getMomentum()),jetParticles.at(k)->getEnergy());
	  
	    // ***************************  Get mcPDG *****************************
	    _object_mcPDG[pfopos] = 0;
	    _object_QuarkOrigin[pfopos] = -1; // -1 if not comming from a quark, 0 for uds, 1 for c and 2 for b
	    _objectFromtau[pfopos] = false;
	    _objectFromBmeson[pfopos] = false;
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
		_object_mcPDG[pfopos] = mcObj->getPDG();
	
		while(mcObj->getParents().size() != 0)
		  {
		    mcObj = mcObj->getParents()[0];
		    if(mcObj->getParents().size() != 0)
		      {
			for(unsigned int parent = 0; parent < mcObj->getParents().size(); parent++)
			  { 
			    MCParticle *mcparents = mcObj->getParents()[parent];
			    if(_objectFromtau[pfopos] == false && abs(mcparents->getPDG()) == 15) _objectFromtau[pfopos] = true;
			    if(_objectFromBmeson[pfopos] == false && (abs(mcparents->getPDG()) == 511 || abs(mcparents->getPDG()) == 521 || abs(mcparents->getPDG()) == 531 || abs(mcparents->getPDG()) == 541 || abs(mcparents->getPDG()) == 5122 || abs(mcparents->getPDG()) == 5332 || abs(mcparents->getPDG()) == 5132 || abs(mcparents->getPDG()) == 5232))
				  {
				    _objectFromBmeson[pfopos]=true;
				  } 
			    if(_object_QuarkOrigin[pfopos] == -1) 
			      {
				if(abs(mcparents->getPDG()) <= 3) _object_QuarkOrigin[pfopos] = 0; 
				if(abs(mcparents->getPDG()) == 4) _object_QuarkOrigin[pfopos] = 1; 
				if(abs(mcparents->getPDG()) == 5) _object_QuarkOrigin[pfopos] = 2; 
			      }
			  }
		      }
		  }
	      }
	    
	    
	    // ******** Apply our LeptonID on all objects with a minimum energy = _leptonEmin
	    // get clusters
	    _Eecal[pfopos] = 0; // ecal energy deposition
	    _Ehcal[pfopos] = 0; // hcal energy deposition
	    _Eyoke[pfopos] = 0; // yoke energy deposition
	    _Ebcal[pfopos] = 0; // bcal energy deposition
	    _Etotal[pfopos] = 0;  // total (ecal+hcal) energy deposition
	    ClusterVec aclustervec = jetParticles.at(k)->getClusters();
	    for(int cluster_i = 0; cluster_i < (int)aclustervec.size(); cluster_i++)
	      {
		Cluster* acluster = aclustervec.at(cluster_i);
		_Eecal[pfopos] += acluster->getSubdetectorEnergies()[0];
		_Eecal[pfopos] += acluster->getSubdetectorEnergies()[3];
		_Ehcal[pfopos] += acluster->getSubdetectorEnergies()[1];
		_Ehcal[pfopos] += acluster->getSubdetectorEnergies()[4];
		_Eyoke[pfopos] += acluster->getSubdetectorEnergies()[2];
		_Ebcal[pfopos] += acluster->getSubdetectorEnergies()[5];
		_Etotal[pfopos] += acluster->getEnergy();
	      }

	    // get particle momentum
	    _Momentum[pfopos] = TVector3(jetParticles.at(k)->getMomentum()).Mag(); 

	    // get tracks momentum - if no track associated : TrackMom = 0
	    _TrackMom[pfopos] = 0; // track momentum
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
		_TrackMom[pfopos] = ParTrack.Mag();
	      }



	    _object_PDG[pfopos] = jetParticles.at(k)->getType();// Use normal PFOID

	    if(_useLeptonID)
	      {
		// Change PDG if the particle satisfies the criteria :
		// 0 < Eecal/Etotal < 0.6 && 0 < Etotal/TrackMom < 0.6 for muons
		// 0.9 < Eecal/Etotal < 1 && 0.2 < Etotal/TrackMom < 1.5 for electrons
		if(_TrackMom[pfopos] > _leptonEmin)
		  {
		    int sign = 1;
		    if(jetParticles.at(k)->getCharge() > 0) sign = -1;
		    if(_Etotal[pfopos]/_TrackMom[pfopos] > 0 && _Etotal[pfopos]/_TrackMom[pfopos] < 0.6 &&
		       _Eecal[pfopos]/_Etotal[pfopos] > 0 && _Eecal[pfopos]/_Etotal[pfopos] < 0.6 ) _object_PDG[pfopos] = 13*sign;

		    else if(_Etotal[pfopos]/_TrackMom[pfopos] > 0.2 && _Etotal[pfopos]/_TrackMom[pfopos] < 1.5 &&
			    _Eecal[pfopos]/_Etotal[pfopos] > 0.9 && _Eecal[pfopos]/_Etotal[pfopos] < 1) _object_PDG[pfopos] = 11*sign;
		  }
	      }
	    

	    // ********** CURE BREMSSTRAHLUNG and/or FSR **********
	    // Look for bremsstralung and FSR photons to correct the lepton energy
	    int photonNb = 0;
	    if((abs(_object_PDG[pfopos]) == 11 || abs(_object_PDG[pfopos]) == 13) && thisPFO.E() > _leptonEmin)
	      {
		_lepton_pos[thispos] = pfopos;
		TVector3 leptonCluster;
		if(jetParticles.at(k)->getClusters().size() > 0) leptonCluster = TVector3(jetParticles.at(k)->getClusters().at(0)->getPosition());
		int pfopos2 = 0;
		photonNb = 0;
		for(int j = 0; j < fourjetscol->getNumberOfElements(); j++)
		  {
		    ReconstructedParticle *recojet2 = dynamic_cast<ReconstructedParticle*>( fourjetscol->getElementAt(j) );
		    ReconstructedParticleVec jetParticles2 = recojet2->getParticles();
		    for(int l = 0; l < (int)jetParticles2.size(); l++)
		      {
			if(usedphoton[pfopos2] == false)
			  {
			    TLorentzVector thisPFO2(TVector3(jetParticles2.at(l)->getMomentum()),jetParticles2.at(l)->getEnergy());
			    if(jetParticles2.at(l)->getType() == 22 && thisPFO2.E() > _leptonEmin && photonNb < 5)
			      {
				float dist = 999.;
				if(jetParticles2.at(l)->getClusters().size() > 0 && jetParticles.at(k)->getClusters().size() > 0)
				  {
				    TVector3 photonCluster(jetParticles2.at(l)->getClusters().at(0)->getPosition());
				    TVector3 vect = photonCluster.Cross(leptonCluster);
				    dist = vect.Mag()/photonCluster.Mag();
				  }
				if((fabs(cos(thisPFO2.Angle(thisPFO.Vect()))) > 0.999 && j == i) // case for FSRs
				   || ( dist < 20. &&  j == i && abs(_object_PDG[pfopos]) == 11)) // simple case for bremsstrahlung
				  {
				    thisPFO += thisPFO2;
				    _photon_pos[thispos][photonNb] = pfopos2;
				    usedphoton[pfopos2] = true;
				    photonNb++;
				  }
			      }
			  }
			pfopos2++;
		      }
		  }
		thispos++;
	      }
	   
	    _object_energy[pfopos] = thisPFO.E();
	    _object_px[pfopos] = thisPFO.Px();
	    _object_py[pfopos] = thisPFO.Py();
	    _object_pz[pfopos] = thisPFO.Pz();
	    _object_momentum[pfopos] = thisPFO.P();
	    _object_phi[pfopos] = thisPFO.Phi();
	    _object_theta[pfopos] = thisPFO.Theta();
	    _object_jetNb[pfopos] = i;
	    
	    // Calculate z4, pT4 and for each PFO
	    _z4[pfopos] = thisPFO.E()/Jet4vec.E();
	    float theta = TMath::Abs(thisPFO.Angle(Jet4vec.Vect()));
	    _pT[pfopos] = thisPFO.P()*sin(theta);
	    _xT[pfopos] = _pT[pfopos]/Jet4vec.M();

	    pfopos++;
	  }
      } // Now all information needed on PFOs is written
    
    
    // ********** ISOLATION **********
    // Identify leptons with cuts on isolation and minimum energy
    for(int mypos = 0; mypos < _nObj; mypos++)
      {
	if(!(_z4[mypos]<0.6 && _xT[mypos]<0.25) &&  _object_energy[mypos] > _leptonEmin)
	  {
	    if(abs(_object_PDG[mypos]) == 11 || abs(_object_PDG[mypos]) == 13) leptonsPositions.push_back(mypos);
	  }
	if(_z4[mypos]<0.6 && _xT[mypos]<0.25 &&  _object_energy[mypos] > _leptonEmin)
	  {
	    /*if(abs(_object_PDG[mypos]) == 11)
	      {
		_nbNotIsoLep++;
		_NotIsoLepcharge[_nbNotIsoLep -1] = -_object_PDG[mypos]/11;
		_NotIsoLepPDG[_nbNotIsoLep -1] = _object_PDG[mypos];
		}*/
	    if(abs(_object_PDG[mypos]) == 13)
	      {
		_nbNotIsoLep++;
		_NotIsoLepcharge[_nbNotIsoLep -1] = -_object_PDG[mypos]/13;
		_NotIsoLepPDG[_nbNotIsoLep -1] = _object_PDG[mypos];
	      }
	  }
      }	
    
    std::cout << "Number of not isolated muons = " << _nbNotIsoLep << std::endl; 

    if(_nObj != pfopos) std::cout << "Number of Objects = " << _nObj << " - Number of PFOs in jets = " << pfopos << std::endl; 

    // ********** SELECT THE BEST LEPTON **********
    // Number of leptons found
    _nbOfLep = leptonsPositions.size();
    for(int ilep = 0; ilep < _nbOfLep; ilep++)
      {
	_candidate_PFO[ilep] = leptonsPositions.at(ilep);
      }
    
    // Apply selection criteria
    int part_nb = -1;
    
    // ********** First : Only one lepton or no lepton found
    if(leptonsPositions.size() == 1) part_nb = leptonsPositions.at(0);
    
    // ********** Second : if more leptons,
    // Highest z lepton
    else if(leptonsPositions.size() > 1)
      {
	float zz = 0.;
	for(unsigned int i = 0; i < leptonsPositions.size(); i++)
	  {
	    int j = leptonsPositions.at(i);
	    float zobj =  _z4[j];
	    // Find best lepton (highest energy)
	    if(zobj > zz)
	      {
		zz = zobj;
		part_nb = j;
	      }
	  } 
      }
    
    // final position of the selected PFO
    _leptonNb = part_nb;
    if(_leptonNb == _posMCpfo) _isGoodLeptonFound = true;
    if(_isMCLeptonFromTau) _lepFromtau = 1;
    int photonIndex = -1;
    for(int i = 0; i < _nObj; i++)
      {
	if(_lepton_pos[i] == _leptonNb) photonIndex = i;
      }
    
    /*
***********************************
Copy of the Reconstructed Particles
***********************************
    */
    // Now that one lepton was selected (with index part_nb), separate the PFOs in two collections
    IMPL::LCCollectionVec * myLepton = new IMPL::LCCollectionVec(EVENT::LCIO::RECONSTRUCTEDPARTICLE);
    myLepton->setSubset(true);
    IMPL::LCCollectionVec * otherPFOs = new IMPL::LCCollectionVec(EVENT::LCIO::RECONSTRUCTEDPARTICLE);
    otherPFOs->setSubset(true);
    
    myLepton->parameters().setValue("NumberOfLeptonsCandidates",_nbOfLep);
    myLepton->parameters().setValue("SelectedLeptonPDG",_object_PDG[_leptonNb]);
    myLepton->parameters().setValue("isGoodLeptonFound",_isGoodLeptonFound);
    myLepton->parameters().setValue("isLeptonFromTau",_lepFromtau);
    
    // copy the RecoParticles
    int pos = 0;
    for(int i = 0 ; i < fourjetscol->getNumberOfElements(); i++)
      {
        ReconstructedParticle *recojet = dynamic_cast<ReconstructedParticle*>( fourjetscol->getElementAt(i) );
        ReconstructedParticleVec jetParticles = recojet->getParticles();
        for(int j = 0; j < (int)jetParticles.size(); j++)
          {
	    ReconstructedParticle* aPFO = dynamic_cast<ReconstructedParticle*>( jetParticles.at(j) );
	  
	    if(pos == _leptonNb)
	      {
		// ROOT stuff
		TLorentzVector lepVec(TVector3(aPFO->getMomentum()),aPFO->getEnergy());
		*_Recolepton = lepVec;
		_Recoz4 = _z4[pos];
		_RecopT = _pT[pos];
		_RecoxT = _xT[pos];
		//
		myLepton->addElement( aPFO );
	      }
	    else if( pos == _photon_pos[photonIndex][0]
		  || pos == _photon_pos[photonIndex][1]
	          || pos == _photon_pos[photonIndex][2]
	          || pos == _photon_pos[photonIndex][3]
	          || pos == _photon_pos[photonIndex][4]
	           )
	      {
		_nbPhotons++;
		myLepton->addElement( aPFO );
	      }
	    else
	      {
		TLorentzVector pfoVec(TVector3(aPFO->getMomentum()),aPFO->getEnergy());
		*_Hadrons += pfoVec;
		otherPFOs->addElement( aPFO );
	      }
	    pos ++;
	  }
      }

    evt->addCollection( otherPFOs , _OutPFOscolName );
    evt->addCollection( myLepton , _Lepton );
    
    // Fill ROOT tree
    _tree->Fill();
    // Return status false if no lepton is found
    if(_nbOfLep == 0) setReturnValue( false );
    else setReturnValue( true );
  }
  
  
  catch( DataNotAvailableException &e) {}
  
}


void NewLeptonFinder::check( LCEvent * evt )
{   
}


void NewLeptonFinder::end()
{ 
  ROOTfile->Write();
  ROOTfile->Close();
  delete ROOTfile;
}
