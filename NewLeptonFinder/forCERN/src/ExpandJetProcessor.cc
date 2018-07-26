#include "ExpandJetProcessor.h"

#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>

// ----- include for verbosity dependend logging ---------
#include <marlin/VerbosityLevels.h>

using namespace lcio ;
using namespace marlin ;

ExpandJetProcessor aExpandJetProcessor ;

ExpandJetProcessor::ExpandJetProcessor()
	: Processor("ExpandJetProcessor") {

		// Processor description
		_description = "" ;

		// register steering parameters: name, description, class-variable, default value
		registerInputCollection( LCIO::RECONSTRUCTEDPARTICLE,
				"InputCollection" ,
				"Input collection of ReconstructedParticles (jets)",
				_inputCollection,
				std::string("InputJets"));

		registerOutputCollection( LCIO::RECONSTRUCTEDPARTICLE,
				"OutputCollection",
				"Output collection of particles that didn't pass the selection",
				_outputCollection,
				std::string("JetPFOs") );
	}


void ExpandJetProcessor::init() { 
	streamlog_out(DEBUG) << "   init called  " << std::endl ;
	printParameters() ;
}

void ExpandJetProcessor::processRunHeader( LCRunHeader* run) { 
} 

void ExpandJetProcessor::processEvent( LCEvent * evt ) { 

	LCCollection* jetCol = evt->getCollection( _inputCollection ) ;

	// Output PFOs
	LCCollectionVec* outCol = new LCCollectionVec( LCIO::RECONSTRUCTEDPARTICLE ) ;
	outCol->setSubset(true) ;

	// jet loop
	unsigned int njet = jetCol->getNumberOfElements();
	for (unsigned int i = 0; i<njet; i++ ) {
		ReconstructedParticle* jet = dynamic_cast<ReconstructedParticle*>( jetCol->getElementAt(i) );
		const ReconstructedParticleVec& pfovec = jet->getParticles();
		for (unsigned int j=0; j<pfovec.size(); ++j) {
			outCol->addElement( pfovec[j] );
		}
	}

	streamlog_out(DEBUG) << "   processing event: " << evt->getEventNumber() 
		<< "   in run:  " << evt->getRunNumber() 
		<< std::endl ;

	// Add PFOs to new collection
	evt->addCollection( outCol, _outputCollection.c_str() );
}

void ExpandJetProcessor::check( LCEvent * evt ) { 
}

void ExpandJetProcessor::end() { 
}
