#ifndef LeptonFinder_h
#define LeptonFinder_h 1

#include "marlin/Processor.h"
#include "lcio.h"
#include <string>
#include "TFile.h"
#include "TTree.h"
#include "TLorentzVector.h"
#include "TMath.h"

#include <vector>

using namespace lcio;
using namespace marlin;


/** 
 *  <h4>Input - Prerequisites</h4>
 *  Needs the collection of PFOs : PandoraPFO
 *
 *  <h4>Output</h4> 
 *  A lepton and a collection of PFOs separated + Number of candidates
 * 
 * @param CollectionName: Name of the PFOs collection
 * @param CollectionName: Name of the 4 jets collection
 * 
 * @author Philippe Doublet, LAL
 * @version $Id: LeptonFinder.h,v 2.0.1 2011/01/19 $ 
 * New in version 2.1.1 :
 * + returnValue(false) if no lepton is found + initial 4 jets
 *
 * New version 3.0 2012/09/27
 * Update for the DBD
 * Now use collection from LCFIPlus
 * 
 */

class LeptonFinder : public Processor {
  
 public:
  
  virtual Processor*  newProcessor() { return new LeptonFinder; }
  
  
  LeptonFinder();

  virtual void init();

  virtual void processRunHeader( LCRunHeader* run );

  virtual void processEvent( LCEvent * evt ); 
    
  virtual void check( LCEvent * evt ); 

  virtual void end();


 protected:
  
  std::string _InPFOscolName;
  std::string _leptonType;
  std::string _OutPFOscolName;
  std::string _Lepton;
  std::string _col4jetsName;
  float _leptonEmin;

  int _nRun;
  int _nEvt;

  TFile* ROOTfile;
  std::string _hfilename;
  TTree* _tree;

  // To store PFOs information
  static const int nObjects = 1000;
  int _nObj;
  float _object_phi[nObjects];
  float _object_theta[nObjects];
  float _object_energy[nObjects];
  float _object_momentum[nObjects];
  float _object_px[nObjects];
  float _object_py[nObjects];
  float _object_pz[nObjects];
  int _object_PDG[nObjects];
  int _object_mcPDG[nObjects];
  int _object_jetNb[nObjects];
  bool _isZombieObject[nObjects]; // zombies are leading-isolated particles (non-leptons)
  bool _objectFromtau[nObjects];
  // PFOID
  float _Eecal[nObjects];
  float _Ehcal[nObjects];
  float _Eyoke[nObjects];
  float _Ebcal[nObjects];
  float _Etotal[nObjects];
  float _TrackMom[nObjects];
  float _object_deltaP[nObjects];
  // Separation criteria + jets
  float _z4[nObjects];
  float _pT[nObjects];
  float _xT[nObjects];
  TLorentzVector *_jet0;
  TLorentzVector *_jet1;
  TLorentzVector *_jet2;
  TLorentzVector *_jet3;
  // For Bremsstrahlung
  int _lepton_pos[nObjects];
  int _photon_pos[nObjects][5];
  // MC information of the lepton
  TLorentzVector *_MClepton;
  float _MCz4;
  float _MCpT;
  float _MCxT;
  bool _isMCLeptonFound;
  bool _isMCLeptonFromTau;
  int _posMCpfo;
  int _MCpdg;
  int _isGoodLeptonFound;
  int _lepFromtau;
  int _isZombieEvent; // tells if zombies appeared

  // Reconstructed lepton
  TLorentzVector *_Recolepton;
  float _Recoz4;
  float _RecopT;
  float _RecoxT;
  int _nbOfLep;
  int _candidate_PFO[nObjects];
  int _leptonNb; // id of the selected lepton among PFOs
  int _nbPhotons; // number of photons added to the lepton
  TLorentzVector *_Hadrons;


};

#endif
